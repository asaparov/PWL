#ifndef NATURAL_DEDUCTION_MH_H_
#define NATURAL_DEDUCTION_MH_H_

#include "theory.h"
#include "natural_deduction.h"

#include <core/random.h>
#include <math/log.h>

constexpr double LOG_2 = 0.693147180559945309417232121458176568075500134360255254120;

using namespace core;

template<typename Formula>
struct proof_transformations {
	array_map<nd_step<Formula, true>*, proof_substitution<Formula, true>> transformed_proofs;

	static inline void free(proof_transformations<Formula>& t) {
		for (auto entry : t.transformed_proofs)
			core::free(entry.value);
		core::free(t.transformed_proofs);
	}
};

template<typename Formula>
bool propose_transformation(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		nd_step<Formula, true>* old_step, nd_step<Formula, true>* new_step)
{
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	array<Proof*> stack(16);
	hash_set<Proof*> visited(32);
	stack[0] = old_step; stack.length++;
	while (stack.length > 0) {
		Proof* step = stack.pop();
		if (!visited.add(step)) return false;

		if (T.observations.contains(step)) {
			/* `step` is an observed logical form in the theory `T` */
			if (!proposed_proofs.transformed_proofs.check_size())
				return false;

			unsigned int index = proposed_proofs.transformed_proofs.index_of(step);
			if (index == proposed_proofs.transformed_proofs.size) {
				proposed_proofs.transformed_proofs.keys[index] = step;
				if (!init(proposed_proofs.transformed_proofs.values[index]))
					return false;
				proposed_proofs.transformed_proofs.size++;
			}

			proof_substitution<Formula, true>& value = proposed_proofs.transformed_proofs.values[index];
			value.map.put(old_step, new_step);
		}

		if (!stack.ensure_capacity(stack.length + step->children.length))
			return false;
		for (Proof* child : step->children) {
			if (visited.contains(child)) continue;
			stack[stack.length++] = child;
		}
	}
}

template<bool FirstSample, typename Formula>
void select_axiom(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		const concept<natural_deduction<Formula, true>>& c,
		array<pair<uint_fast8_t, unsigned int>>& axiom_indices,
		array<unsigned int>& selected_types, array<unsigned int>& selected_negated_types,
		array<relation>& selected_relations, array<relation>& selected_negated_relations,
		array<unsigned int>& intersection)
{
	unsigned int type_value; relation relation_value;
	/* TODO: This should be more likely if the resulting intersecting subset is
	   closer to the set of all concepts with the given `type`. `A` being
	   "closer" to being a subset of `B` can be measured in different ways:
		(1) if `A_1` has fewer elements not in `B` than `A_2` (i.e.
		    |A_1 \ B| < |A_2 \ B|), then `A_1`is closer than `A_2`,
		(2) if `A_1` has more elements in `B` than `A_2` (i.e.
			|A_1 intersect B| > |A_2 intersect B|), then `A_1` is closer than
			`A_2`. */
	unsigned int selected_index = sample_uniform(axiom_indices.length);
	switch (axiom_indices[selected_index].key) {
	case 0:
		type_value = c.types.keys[axiom_indices[selected_index].value];
		if (!selected_types.add(type_value)) return false;
		if (FirstSample) intersection.append(T.types.get(type_value).key);
		else set_intersect(intersection, T.types.get(type_value).key);
	case 1:
		type_value = c.negated_types.keys[axiom_indices[selected_index].value];
		if (!selected_negated_types.add(type_value)) return false;
		if (FirstSample) intersection.append(T.types.get(type_value).value);
		else set_intersect(intersection, T.types.get(type_value).value);
	case 2:
		relation_value = c.relations.keys[axiom_indices[selected_index].value];
		if (!selected_relations.add(relation_value)) return false;
		if (FirstSample) intersection.append(T.relations.get(relation_value).key);
		else set_intersect(intersection, T.relations.get(relation_value).key);
	case 3:
		relation_value = c.negated_relations.keys[axiom_indices[selected_index].value];
		if (!selected_negated_relations.add(relation_value)) return false;
		if (FirstSample) intersection.append(T.relations.get(relation_value).value);
		else set_intersect(intersection, T.relations.get(relation_value).value);
	}

	axiom_indices.remove(selected_index);
}

template<typename Formula, bool Negated>
inline bool get_axiom(
		const array_map<unsigned int, nd_step<Formula, true>*>& types,
		unsigned int predicate, array<nd_step<Formula, true>*>& conjunct_steps)
{
#if !defined(NDEBUG)
	if (!types.contains(predicate))
		fprintf(stderr, "get_axiom WARNING: The given `predicate` does not exist in the map `types`.\n");
#endif
	return conjunct_steps.add(types.get(predicate));
}

template<typename Formula, bool Negated>
inline bool get_axiom(
		const array_map<relation, nd_step<Formula, true>*>& relations,
		relation rel, array<nd_step<Formula, true>*>& conjunct_steps)
{
#if !defined(NDEBUG)
	if (!relations.contains(rel))
		fprintf(stderr, "get_axiom WARNING: The given `rel` does not exist in the map `relations`.\n");
#endif
	return conjunct_steps.add(relations.get(rel));
}

template<bool Negated, unsigned int Arity>
struct atom {
	unsigned int predicate;
	unsigned int args[Arity];
};

template<typename Formula, unsigned int Index, bool Negated, unsigned int Arity, typename... Formulas>
inline Formula* _helper(const atom<Negated, Arity> a, Formulas&&... args) {
	if (Index == 0)
		return Formula::new_atom(a.predicate, std::forward<Formulas>(args)...);
	else return _helper<Formula, Index - 1>(a,
			(a.args[Index - 1] == 0 ? Formula::new_variable(1) : Formula::new_constant(a.args[Index - 1])), args);
}

template<typename Formula, bool Negated, unsigned int Arity>
inline Formula* (const atom<Negated, Arity> a) {
	return _helper<Formula, Arity>(a);
}

template<bool Negated, typename Formula>
inline const array<unsigned int>& get_negated_set(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		const atom<Negated, 1> a)
{
	return (Negated ? T.types.get(a.predicate).key : T.types.get(a.predicate).value);
}

template<bool Negated, typename Formula>
const array<unsigned int>& get_negated_set(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		const atom<Negated, 2> a)
{
	relation r = { a.predicate, a.args[0], a.args[1] };
	return (Negated ? T.relations.get(r).key : T.relations.get(r).value);
}

template<bool Negated, typename Formula>
const nd_step<Formula, true>* get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		const atom<Negated, 1> a)
{
	array_map<unsigned int, nd_step<Formula, true>*>& proofs = Negated ? instance.negated_types : instance.types;
	unsigned int index = proofs.index_of(a.predicate);
#if !defined(NDEBUG)
	if (index == proofs.size) {
		fprintf(stderr, "propose_atom_generalization WARNING: Theory is invalid.\n");
		return NULL;
	}
#endif
	return proofs.values[index];
}

template<bool Negated, typename Formula>
const nd_step<Formula, true>* get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		const atom<Negated, 2> a)
{
	array_map<relation, nd_step<Formula, true>*>& proofs = Negated ? instance.negated_relations : instance.relations;
	relation r = { a.predicate, a.args[0], a.args[1] };
	unsigned int index = proofs.index_of(r);
#if !defined(NDEBUG)
	if (index == proofs.size) {
		fprintf(stderr, "propose_atom_generalization WARNING: Theory is invalid.\n");
		return NULL;
	}
#endif
	return proofs.values[index];
}

template<bool FirstSample, bool Negated, typename Formula>
inline void get_satisfying_concepts_helper(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		const Formula* literal, array<unsigned int>& intersection)
{
	typedef typename Formula::TermType TermType;

	if (literal->atom.arg2.type == TermType::NONE) {
		/* this is a unary atom */
		const pair<array<unsigned int>, array<unsigned int>>& list = T.types.get(literal->atom.predicate);
		const array<unsigned int>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist);
		else set_intersect(intersection, sublist);
	} else {
		/* this is a binary atom */
		relation r = { literal->atom.predicate,
				(literal->atom.arg1.type == TermType::VARIABLE ? 0 : literal->atom.arg1.constant),
				(literal->atom.arg2.type == TermType::VARIABLE ? 0 : literal->atom.arg2.constant) };
		const pair<array<unsigned int>, array<unsigned int>>& list = T.relations.get(r);
		const array<unsigned int>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist);
		else set_intersect(intersection, sublist);
	}
}

template<bool FirstSample, typename Formula>
bool get_satisfying_concepts(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		const Formula* formula, array<unsigned int>& intersection)
{
	typedef typename Formula::Type FormulaType;

	if (formula->type == FormulaType::ATOM) {
		get_satisfying_concepts_helper<FirstSample, false>(T, formula, intersection);
	} else if (formula->type == FormulaType::NOT && formula->unary.operand->type == FormulaType::ATOM) {
		get_satisfying_concepts_helper<FirstSample, true>(T, formula->unary->operand, intersection);
	} else if (formula->type == FormulaType::AND) {
		for (unsigned int i = 0; i < antecedent->array.length; i++) {
			Formula* conjunct = antecedent->array.operand[i];
			if ((i == 0 && !get_satisfying_concepts<FirstSample>(T, conjunct, antecedent_concepts))
			 || (i > 0 && !get_satisfying_concepts<false>(T, conjunct, antecedent_concepts)))
				return false;
		}
	} else {
		fprintf(stderr, "get_satisfying_concepts ERROR: Expected a literal.\n");
		return false;
	}
	return true;
}

template<typename Formula>
bool is_subset(
		const Formula** first, unsigned int first_length,
		const Formula** second, unsigned int second_length)
{
	unsigned int i = 0, j = 0;
	while (i < first_length && j < second_length)
	{
		if (first[i] == second[j] || *first[i] == *second[j]) {
			i++; j++;
		} else if (*first[i] < *second[j]) {
			return false;
		} else {
			j++;
		}
	}
	return (i == first_length);
}

template<typename Formula>
bool is_subset(const Formula* first, const Formula* second)
{
	typedef typename Formula::Type FormulaType;

	if (first->type == FormulaType::AND) {
		if (second->type == Formula::AND) {
			return is_subset(first->array.operands, first->array.length, second->array.operands, second->array.length);
		} else {
			return false;
		}
	} else {
		if (second->type == Formula::AND) {
			return is_subset(&first, 1, second->array.operands, second->array.length);
		} else {
			return (first == second || *first == *second);
		}
	}
}

template<bool Negated, unsigned int Arity, typename Formula>
bool propose_universal_intro(
		theory<Formula, natural_deduction<Formula, true>>& T,
		const atom<Negated, Arity> a,
		unsigned int concept,
		double& log_proposal_probability_ratio)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula, true>> ConceptType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	proof_transformations<Formula> proposed_proofs;

	/* NOTE: We can't check if there is an existing universally-quantified
			 formula that implies the new formula (that we propose in this
			 function), or vice versa, since the inverse proposal has zero
			 probability. Rather, we constrain the theory such that for any
			 literal observation, if it can be proved by a universally-
			 quantified axiom, we initialize its proof to use that axiom,
			 rather than adding a trivial proof of the literal. */

	/* compute the set of constants for which the selected axiom is true */
	const array<unsigned int>& negated_set = get_negated_set(T, a);

	/* find other ground axioms that are connected to this constant */
	const ConceptType& c = T.ground_concepts.get(concept);
	array<pair<uint_fast8_t, unsigned int>> axiom_indices(
			c.types.size + c.negated_types.size + c.relations.size + c.negated_relations.size);
	for (unsigned int i = 0; i < c.types.size; i++) {
		if (!Negated && Arity == 1 && c.types.keys[i] == a.predicate) continue; /* avoid selecting the consequent as a conjunct in the antecedent */
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(0, i));
	} for (unsigned int i = 0; i < c.negated_types.size; i++) {
		if (Negated && Arity == 1 && c.negated_types.keys[i] == a.predicate) continue; /* avoid selecting the consequent as a conjunct in the antecedent */
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(1, i));
	} for (unsigned int i = 0; i < c.relations.size; i++) {
		/* avoid selecting the consequent as a conjunct in the antecedent */
		if (!Negated && Arity == 2 && c.relations.keys[i].type == a.predicate
		 && c.relations.keys[i].arg1 == a.args[0]
		 && c.relations.keys[i].arg2 == a.args[1]) continue;
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(2, i));
	} for (unsigned int i = 0; i < c.negated_relations.size; i++) {
		/* avoid selecting the consequent as a conjunct in the antecedent */
		if (Negated && Arity == 2 && c.negated_relations.keys[i].type == a.predicate
		 && c.negated_relations.keys[i].arg1 == a.args[0]
		 && c.negated_relations.keys[i].arg2 == a.args[1]) continue;
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(3, i));
	}

	if (axiom_indices.length == 0) {
		/* no candidate axioms are available */
		return true;
	}

	array<unsigned int> selected_types(8);
	array<unsigned int> selected_negated_types(8);
	array<relation> selected_relations(8);
	array<relation> selected_negated_relations(8);

	unsigned int count = sample_uniform(axiom_indices.length - 1) + 1;
	array<unsigned int> intersection(32);
	log_proposal_probability_ratio -= -log_cache<double>::instance().get(axiom_indices.length - 1)
			- lgamma(axiom_indicies.length + 1) + lgamma(count + 1) + lgamma(axiom_indicies.length - count + 1);
	for (unsigned int i = 0; i < count; i++) {
		/* select a candidate axiom at random, and compute the set of constants
		   for which it's true and intersect it with `intersection` */
		if (i == 0) {
			select_axiom<true>(T, c, axiom_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		} else {
			select_axiom<false>(T, c, axiom_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		}
	}

	if (has_intersection(intersection.data, intersection.length, negated_set.data, negated_set.length))
		return true;

	/* compute the probability of the inverse proposal */
	unsigned int new_axiom_count = T.ground_axiom_count + T.universal_quantifications.length + 1;
	if (Arity == 2 && a.args[0] == a.args[1]) new_axiom_count -= 3 * intersection.length;
	else if (Arity == 2) new_axiom_count -= 2 * intersection.length;
	else new_axiom_count -= intersection.length;
	if (!log_cache<double>::instance().ensure_size(new_axiom_count)) return false;
	log_proposal_probability_ratio += -log_cache<double>::instance().get(new_axiom_count);

	/* we've finished constructing the new universally-quantified formula,
	   and now we need to transform all the relevant proofs */
	array<Formula*> conjuncts(32);
	for (unsigned int type : selected_types) {
		Formula* conjunct = Formula::new_atom(type, Formula::new_variable(1));
		if (conjunct == NULL || !conjuncts.add(conjunct)) {
			free_formulas(conjuncts); return false;
		}
	}
	for (unsigned int type : selected_negated_types) {
		Formula* conjunct = Formula::new_not(Formula::new_atom(type, Formula::new_variable(1))));
		if (conjunct == NULL || !conjuncts.add(conjunct)) {
			free_formulas(conjuncts); return false;
		}
	}
	for (const relation& r : selected_relations) {
		Formula* conjunct = Formula::new_atom(
				(r.type == 0) ? Formula::new_variable(1) : Formula::new_constant(r.type),
				(r.arg1 == 0) ? Formula::new_variable(1) : Formula::new_constant(r.arg1),
				(r.arg2 == 0) ? Formula::new_variable(1) : Formula::new_constant(r.arg2));
		if (conjunct == NULL || !conjuncts.add(conjunct)) {
			free_formulas(conjuncts); return false;
		}
	}
	for (const relation& r : selected_negated_relations) {
		Formula* conjunct = Formula::new_not(Formula::new_atom(
				(r.type == 0) ? Formula::new_variable(1) : Formula::new_constant(r.type),
				(r.arg1 == 0) ? Formula::new_variable(1) : Formula::new_constant(r.arg1),
				(r.arg2 == 0) ? Formula::new_variable(1) : fol_forFormulamula::new_constant(r.arg2)));
		if (conjunct == NULL || !conjuncts.add(conjunct)) {
			free_formulas(conjuncts); return false;
		}
	}

	Formula* axiom = Formula::new_for_all(1, Formula::new_if_then(
			Formula::new_and(conjuncts), <Formula>(a)));
	free_formulas(conjuncts);
	if (axiom == NULL) return false;

	Formula* canonicalized = canonicalize(axiom);
	free(*axiom); if (axiom->reference_count == 0) free(axiom);
	if (canonicalized == NULL)
		return false;

	nd_step<Formula, true>* axiom_step = ProofCalculus::new_axiom(canonicalized);
	free(*canonicalized);
	if (canonicalized->reference_count == 0)
		free(canonicalized);
	if (axiom_step == NULL)
		return false;

	/* check if the new axiom already exists in the theory */
	for (Proof* old_axiom : T.universal_quantifications) {
		if (*old_axiom->formula == *canonicalized) {
			free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
			return true; /* we fail here so that the inverse proposal is easier to implement */
		}
	}

	array<nd_step<Formula, Canonical>*> conjunct_steps(antecedent->array.length);
	for (unsigned int concept : intersection) {
		const ConceptType& instance = T.ground_concepts.get(concept);

		for (unsigned int i = 0; i < antecedent->array.length; i++) {
			Formula* atom = antecedent->array.operands[i];
			if (atom->type == FormulaType::NOT) {
				atom = atom->unary.operand;
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, true>(instance.negated_types, atom->predicate, conjunct_steps)) {
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
						return false;
					}
				} else {
					const relation rel = { atom->predicate,
							atom->arg1.type == TermType::CONSTANT ? atom->arg1.constant : 0,
							atom->arg2.type == TermType::CONSTANT ? atom->arg2.constant : 0 };
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, conjunct_steps)) {
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
						return false;
					}
				}
			} else {
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, false>(instance.types, atom->predicate, conjunct_steps)) {
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
						return false;
					}
				} else {
					const relation rel = { atom->predicate,
							atom->arg1.type == TermType::CONSTANT ? atom->arg1.constant : 0,
							atom->arg2.type == TermType::CONSTANT ? atom->arg2.constant : 0 };
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, conjunct_steps)) {
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
						return false;
					}
				}
			}
		}
	}

	for (unsigned int k = 0; k < intersection.length; k++) {
		const unsigned int concept = intersection[k];
		ConceptType& instance = T.ground_concepts.get(concept);

		Proof* old_step = get_proof(instance, a);
#if !defined(NDEBUG)
		else if (old_step->type != nd_step_type::AXIOM)
			fprintf(stderr, "propose_atom_generalization WARNING: Expected an axiom.\n");
#endif

		Proof* new_step = ProofCalculus::new_implication_elim(
				ProofCalculus::new_universal_elim(axiom_step, Formula::new_constant(concept)),
				ProofCalculus::new_conjunction_intro(conjunct_steps));
		if (new_step == NULL) {
			for (unsigned int i = k; k < conjunct_steps.length; k++) {
				free(*conjunct_steps[k]);
				if (conjunct_steps[k]->reference_count == 0)
					free(conjunct_steps[k]);
			}
			free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
			return false;
		}
		new_step->reference_count++;

		if (!propose_transformation(T, proposed_proofs, old_step, new_step)) {
			free(*new_step);
			if (new_step->reference_count == 0)
				free(new_step);
			for (unsigned int i = k; k < conjunct_steps.length; k++) {
				free(*conjunct_steps[k]);
				if (conjunct_steps[k]->reference_count == 0)
					free(conjunct_steps[k]);
			}
			free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
			return false;
		}
	}

	return true;
}

template<typename Formula>
void free_proofs(nd_step<Formula, true>** proofs, unsigned int count) {
	for (unsigned int i = 0; i < count; i++) {
		if (proofs[i] == NULL) continue;
		free(*proofs[i]);
		if (proofs[i]->reference_count == 0)
			free(proofs[i]);
	}
}

template<typename Formula>
inline bool make_grounded_conjunct(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		Formula* lifted_conjunct, unsigned int variable,
		typename Formula::Term constant,
		nd_step<Formula, true>** new_axioms,
		unsigned int i,
		hash_map<unsigned int, nd_step<Formula, true>*>& new_grounded_axioms)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, true> ProofCalculus;

	if (new_axioms[i] != NULL) return true;

	/* first check if the grounded conjunct already exists as an axiom in the theory */
	const Formula* operand = lifted_conjunct;
	if (operand->type == FormulaType::ATOM) {
		if (operand->atom.arg2.type == TermType::NONE) {
			/* atom is unary */
			Proof* proof = c.types.get(operand->atom.predicate, contains);
			if (contains)
				/* we have to fail here as the transformation is no longer reversible */
				return true;
		} else {
			/* atom is binary */
			relation r = { operand->atom.predicate,
					(operand->atom.arg1.type == TermType::VARIABLE ? 0 : operand->atom.arg1.constant),
					(operand->atom.arg2.type == TermType::VARIABLE ? 0 : operand->atom.arg2.constant) };
			Proof* proof = c.relations.get(r, contains);
			if (contains)
				/* we have to fail here as the transformation is no longer reversible */
				return true;
		}
	} else if (operand->type == FormulaType::NOT && operand->unary.operand->type == FormulaType::ATOM) {
		operand = operand->unary.operand;
		if (operand->atom.arg2.type == TermType::NONE) {
			/* atom is unary */
			Proof* proof = c.negated_types.get(operand->atom.predicate, contains);
			if (contains)
				/* we have to fail here as the transformation is no longer reversible */
				return true;
		} else {
			/* atom is binary */
			relation r = { operand->atom.predicate,
					(operand->atom.arg1.type == TermType::VARIABLE ? 0 : operand->atom.arg1.constant),
					(operand->atom.arg2.type == TermType::VARIABLE ? 0 : operand->atom.arg2.constant) };
			Proof* proof = c.negated_relations.get(r, contains);
			if (contains)
				/* we have to fail here as the transformation is no longer reversible */
				return true;
		}
	} else {
		fprintf(stderr, "make_grounded_conjunct ERROR: Expected a conjunction of literals.\n");
		return false;
	}

	new_axioms[i] = ProofCalculus::new_axiom(
			substitute(*lifted_conjunct, Formula::new_variable(variable), constant));
	if (new_axioms[i] == NULL)
		return false;
	new_axioms[i]->reference_count++;

	if (!new_grounded_axioms.put(constant.constant, new_axioms[i])) {
		free(*new_axioms[i]); free(new_axioms[i]);
		new_axioms[i] = NULL; return false;
	}
	return true;
}

template<typename Formula>
bool make_grounded_conjunction(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		Formula** conjuncts, unsigned int conjunct_count,
		unsigned int variable, typename Formula::Term constant,
		nd_step<Formula, true>** new_axioms,
		hash_map<unsigned int, nd_step<Formula, true>*>& new_grounded_axioms,
		nd_step<Formula, true>*& new_step)
{
	typedef natural_deduction<Formula, true> ProofCalculus;

	bool contains;
	const concept<ProofCalculus>& c = T.ground_concepts.get(constant.constant, contains);
#if !defined(NDEBUG)
	if (!contains)
		fprintf(stderr, "make_grounded_conjunction WARNING: Theory does not contain a concept for ID %u.\n", constant.constant);
#endif
	for (unsigned int i = 0; i < consequent_count; i++) {
		if (!make_grounded_conjunct(T, conjuncts[i], variable, constant, new_axioms, i, new_grounded_axioms)) {
			free_proofs(new_axioms, i);
			return false;
		} if (new_axioms[i] == NULL) {
			new_step = NULL;
			return true;
		}
	}

	if (consequent_count == 1)
		return new_axioms[0];

	new_step = ProofCalculus::new_conjunction_intro(array_view(new_axioms, consequent_count));
	if (new_step == NULL) {
		free_proofs(new_axioms, consequent_count);
		return false;
	}
	new_step->reference_count++;
	return true;
}

template<typename Formula,
	typename TheoryPrior, typename AxiomPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior>
bool propose_universal_elim(
		theory<Formula, natural_deduction<Formula, true>>& T,
		unsigned int axiom_index,
		double& log_proposal_probability_ratio,
		double log_proof_stop_probability,
		double log_proof_continue_probability,
		TheoryPrior& theory_prior, AxiomPrior& axiom_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula, true>> ConceptType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	proof_transformations<Formula> proposed_proofs;

	/* first check that all uses of this axiom eliminate the universal quantification; otherwise, we can't remove it */
	Proof* step = T.universal_quantifications[axiom_index];
	if (T.observations.contains(axiom->formula)) return true;

	/* this proposal is not invertible if the consequent is not atomic */
	bool is_consequent_binary; bool is_consequent_symmetric;
	unsigned int consequent_length;
	Formula* consequent = axiom->formula->quantifier.operand->binary.right;
	if (consequent->type == FormulaType::ATOM) {
		consequent_length = 1;
		Formula* atom = consequent;
		is_consequent_binary = (atom->atom.arg2.type != TermType::NONE);
		is_consequent_symmetric = false;
	} else if (consequent->type == FormulaType::NOT && consequent->unary.operand->type == FormulaType::ATOM) {
		consequent_length = 1;
		Formula* atom = consequent->unary.operand;
		is_consequent_binary = (atom->atom.arg2.type != TermType::NONE);
		is_consequent_symmetric = (atom->atom.arg1.type == TermType::VARIABLE && atom->atom.arg2.type == TermType::VARIABLE);
	} else if (consequent->type == FormulaType::AND) {
		consequent_length = consequent->array.length;
		return true;
	} else {
		return true;
	}

	unsigned int antecedent_length;
	Formula* antecedent = axiom->formula->quantifier.operand->binary.left;
	if (antecedent->type == FormulaType::ATOM {
		antecedent_length = 1;
	} else if (antecedent->type == FormulaType::NOT && antecedent->unary.operand->type == FormulaType::ATOM) {
		antecedent_length = 1;
	} else if (antecedent->type == FormulaType::AND) {
		antecedent_length = antecedent->array.length;
	} else {
		return true;
	}

	hash_map<unsigned int, Proof*> new_grounded_axioms(16);
	for (Proof* child : axiom->children) {
		if (child->type != nd_step_type::UNIVERSAL_ELIMINATION
		 || child->operands[1]->term.type == TermType::CONSTANT)
			return true;

		for (Proof* grandchild : child->children) {
			if (grandchild->type != nd_step_type::IMPLICATION_ELIMINATION)
				return true;

			Proof* new_step = NULL;
			Proof** new_axioms = (Proof**) calloc(consequent_length, sizeof(Proof*));
			if (new_axioms == NULL) return false;

			/* check if `grandchild` is observed */
			if (T.observations.contains(grandchild)) {
				if (!make_grounded_conjunction(T,
						(consequent->type == FormulaType::AND ? consequent->array.operands : &consequent),
						(consequent->type == FormulaType::AND ? consequent_length : 1),
						formula->quantifier.variable, child->operands[1]->term, new_axioms, new_grounded_axioms, new_step)
				 || new_step == NULL)
				{
					bool fail = (new_step == NULL);
					free_proofs(new_axioms);
					free(new_axioms); return !fail;
				}

				/* propose `new_step` to substitute `grandchild` */
				if (!propose_transformation(T, proposed_proofs, grandchild, new_step)) {
					free(*new_step); if (new_step->reference_count == 0) free(new_step);
					free_proofs(new_axioms); free(new_axioms); return false;
				}
			}
			for (Proof* descendant : grandchild->children) {
				if (descendant->type != nd_step_type::CONJUNCTION_ELIMINATION) {
					/* we need to prove all conjuncts in `grandchild` */
					if (new_step == NULL) {
						if (!make_grounded_conjunction(T,
								(consequent->type == FormulaType::AND ? consequent->array.operands : &consequent),
								(consequent->type == FormulaType::AND ? consequent_length : 1),
								formula->quantifier.variable, child->operands[1]->term, new_axioms, new_grounded_axioms)
						 || new_step == NULL)
						{
							bool fail = (new_step == NULL);
							free_proofs(new_axioms);
							free(new_axioms); return !fail;
						}

						/* propose `new_step` to substitute `grandchild` if we haven't already */
						if (!propose_transformation(T, proposed_proofs, grandchild, new_step)) {
							free(*new_step); if (new_step->reference_count == 0) free(new_step);
							free_proofs(new_axioms); free(new_axioms); return false;
						}
					}
				} else {
					const array<unsigned int>& indices = descendant->operands[1]->parameters;
					for (unsigned int index : indices) {
						if (!make_grounded_conjunct(T, consequent->array.operands[index],
								formula->quantifier.variable, child->operands[1]->term, new_axioms, index, new_grounded_axioms)
						 || new_axioms[index] == NULL)
						{
							bool fail = (new_axioms[index] == NULL);
							if (new_step != NULL) {
								free(*new_step); if (new_step->reference_count == 0) free(new_step);
							}
							free_proofs(new_axioms); free(new_axioms); return !fail;
						}
					}

					Proof* new_indexed_step = ProofCalculus::new_conjunction_intro(
							indexed_array_view(new_axioms, indices.data, indices.length));
					if (new_indexed_step == NULL) {
						if (new_step != NULL) {
							free(*new_step); if (new_step->reference_count == 0) free(new_step);
						}
						free_proofs(new_axioms); free(new_axioms); return false;
					}

					/* propose `new_indexed_step` to substitute `descendant` */
					if (!propose_transformation(T, proposed_proofs, descendant, new_indexed_step)) {
						if (new_step != NULL) {
							free(*new_step); if (new_step->reference_count == 0) free(new_step);
						}
						free_proofs(new_axioms); free(new_axioms); return false;
					}
				}
			}
			free_proofs(new_axioms);
			free(new_axioms);
		}
	}

	/* Compute the probability of inverting this proposal. NOTE: the change to
	   `T` that is being proposed here is removing an element from
	   `T.universal_quantifications` and adding an axiom for every conjunct in
	   the consequent and for every constant that is used to eliminate the
	   universal quantification in the aforementioned axiom. */
	array<unsigned int> intersection(64);
	if (!get_satisfying_concepts<true>(T, antecedent, intersection))
		return false;
	/* if the antecedent conjuncts aren't grounded, then we can't propose the inverse transformation */
	if (intersection.length == 0) return true;
	unsigned int new_axiom_count = T.ground_axiom_count + T.universal_quantifications.length - 1;
	if (is_consequent_symmetric) new_axiom_count += new_grounded_axioms.size * 3;
	else if (is_consequent_binary) new_axiom_count += new_grounded_axioms.size * 2;
	else new_axiom_count += new_grounded_axioms.size;
	if (!log_cache<double>::instance().ensure_size(new_axiom_count)) return false;
	double* log_probabilities = (double*) calloc(intersection.length, sizeof(double));
	for (unsigned int i = 0; i < intersection.length; i++) {
		const ConceptType& c = T.ground_concepts.get(intersection[i]);
		unsigned int count = c.types.length + c.negated_types.length + c.relations.length + c.negated_relations.length;
		if (new_grounded_axioms.table.contains(intersection[i])) count--;

		log_probabilities[i] += -log_cache<double>::instance().get(count - 1)
				- lgamma(count + 1) + lgamma(antecedent_length + 1) + lgamma(count - antecedent_length + 1)
				- log_cache<double>::instance().get(new_axiom_count);
	}
	log_proposal_probability_ratio += logsumexp(log_probabilities, intersection.length);
	free(log_probabilities);

	/* go ahead and compute the probability ratios and perform the accept/reject step */
	if (consequent->type == FormulaType::ATOM) {
		Formula* atom = consequent;
		if (atom->atom.arg2.type == TermType::NONE) {
			atom<false, 1> consequent_atom;
			consequent_atom.predicate = atom->atom.predicate;
			consequent_atom.args[0] = 0;
			return do_mh_universal_elim(T, proposed_proofs, axiom_index, new_grounded_axioms, consequent_atom,
				log_proposal_probability_ratio, log_proof_stop_probability, log_proof_continue_probability,
				theory_prior, axiom_prior, universal_introduction_prior, universal_elimination_prior);
		} else {
			atom<false, 2> consequent_atom;
			consequent_atom.predicate = atom->atom.predicate;
			consequent_atom.args[0] = (atom->atom.arg1.type == TermType::VARIABLE ? 0 : atom->atom.arg1.constant);
			consequent_atom.args[1] = (atom->atom.arg2.type == TermType::VARIABLE ? 0 : atom->atom.arg2.constant);
			return do_mh_universal_elim(T, proposed_proofs, axiom_index, new_grounded_axioms, consequent_atom,
				log_proposal_probability_ratio, log_proof_stop_probability, log_proof_continue_probability,
				theory_prior, axiom_prior, universal_introduction_prior, universal_elimination_prior);
		}
	} else if (consequent->type == FormulaType::NOT && consequent->unary.operand->type == FormulaType::ATOM) {
		Formula* atom = consequent->unary.operand;
		if (atom->atom.arg2.type == TermType::NONE) {
			atom<true, 1> consequent_atom;
			consequent_atom.predicate = atom->atom.predicate;
			consequent_atom.args[0] = 0;
			return do_mh_universal_elim(T, proposed_proofs, axiom_index, new_grounded_axioms, consequent_atom,
				log_proposal_probability_ratio, log_proof_stop_probability, log_proof_continue_probability,
				theory_prior, axiom_prior, universal_introduction_prior, universal_elimination_prior);
		} else {
			atom<true, 2> consequent_atom;
			consequent_atom.predicate = atom->atom.predicate;
			consequent_atom.args[0] = (atom->atom.arg1.type == TermType::VARIABLE ? 0 : atom->atom.arg1.constant);
			consequent_atom.args[1] = (atom->atom.arg2.type == TermType::VARIABLE ? 0 : atom->atom.arg2.constant);
			return do_mh_universal_elim(T, proposed_proofs, axiom_index, new_grounded_axioms, consequent_atom,
				log_proposal_probability_ratio, log_proof_stop_probability, log_proof_continue_probability,
				theory_prior, axiom_prior, universal_introduction_prior, universal_elimination_prior);
		}
	}
}

template<typename Formula>
bool transform_proofs(const proof_transformations<Formula>& proposed_proofs)
{
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	hash_map<Proof*, Proof*> transformations(32);
	for (auto entry : proposed_proofs.transformed_proofs) {
		if (!transformations.check_size(transformations.size + entry.value.map.size))
			return false;
		for (auto transformation : entry.value.map)
			transformations.put(transformation.key, transformation.value);
	}

	for (auto entry : transformations) {
		if (!entry.value.children.ensure_capacity(entry.key.children.length))
			return false;

		/* change the operand of any child to the new proof step */
		for (unsigned int i = entry.key.children.length; i > 0; i--) {
			Proof* child = entry.key.children[i - 1];

			bool error = false;
			switch (child->type) {
			case nd_step_type::CONJUNCTION_INTRODUCTION:
				for (unsigned int j = 0; j < child->operand_array.length; j++) {
					if (child->operand_array[j] == entry.key) {
						child->operand_array[j] = entry.value;
						entry.key->reference_count--;
						if (entry.key->reference_count == 0) free(entry.key);
						entry.value->reference_count++;
					}
				}
				break;
			case nd_step_type::CONJUNCTION_ELIMINATION:
			case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
			case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
			case nd_step_type::DISJUNCTION_INTRODUCTION:
			case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
			case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
			case nd_step_type::DISJUNCTION_ELIMINATION:
			case nd_step_type::IMPLICATION_INTRODUCTION:
			case nd_step_type::IMPLICATION_ELIMINATION:
			case nd_step_type::BICONDITIONAL_INTRODUCTION:
			case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT:
			case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT:
			case nd_step_type::PROOF_BY_CONTRADICTION:
			case nd_step_type::NEGATION_ELIMINATION:
			case nd_step_type::UNIVERSAL_INTRODUCTION:
			case nd_step_type::UNIVERSAL_ELIMINATION:
			case nd_step_type::EXISTENTIAL_INTRODUCTION:
			case nd_step_type::EXISTENTIAL_ELIMINATION:
				for (unsigned int j = 0; j < ND_OPERAND_COUNT; j++) {
					if (child->operands[j] == entry.key) {
						child->operands[j] = entry.value;
						entry.key->reference_count--;
						if (entry.key->reference_count == 0) free(entry.key);
						entry.value->reference_count++;
					}
				}
				break;
			case nd_step_type::TERM_PARAMETER:
			case nd_step_type::ARRAY_PARAMETER:
			case nd_step_type::AXIOM:
			case nd_step_type::FORMULA_PARAMETER:
				error = true; break;
			}
			if (error) {
				fprintf(stderr, "transform_proofs ERROR: Invalid nd_step_type.\n");
				return false;
			}

			entry.value.children[entry.value.children.length++] = child;
			entry.key.children.length--;
		}
	}
	return true;
}

template<typename Formula,
	bool Negated, unsigned int Arity,
	typename TheoryPrior, typename AxiomPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior>
bool do_mh_universal_intro(
		theory<Formula, natural_deduction<Formula, true>>& T,
		const proof_transformations<Formula>& proposed_proofs,
		nd_step<Formula, true>* new_universal_quantification,
		double log_proposal_probability_ratio,
		double log_proof_stop_probability,
		double log_proof_continue_probability,
		TheoryPrior& theory_prior, AxiomPrior& axiom_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior)
{
	/* compute the proof portion of the prior for both current and proposed theories */
	for (const auto& entry : proposed_proofs.transformed_proofs) {
		/* contribution from the old proof */
		log_proposal_probability_ratio -= log_probability(
				*entry.key, log_proof_stop_probability,
				log_proof_continue_probability, axiom_prior,
				universal_introduction_prior, universal_elimination_prior);

		/* contribution from the new proof */
		log_proposal_probability_ratio += log_probability(
				*entry.key, log_proof_stop_probability,
				log_proof_continue_probability, axiom_prior,
				universal_introduction_prior, universal_elimination_prior,
				entry.value);
	}

	/* compute the prior of the current and proposed theories and add them to `log_proposal_probability_ratio` */
	log_proposal_probability_ratio += log_probability_ratio(T, theory_prior /* add arguments here */);

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we've accepted the proposal */
		if (!transform_proofs(proposed_proofs)) return false;

		if (!T.add_universal_quantification(new_universal_quantification)) return false;
		for (auto entry : new_grounded_axioms)
			if (!add_ground_axiom(T, consequent_atom, entry.key, entry.value)) return false;
	}
	return true;
}

template<bool Negated, typename Formula>
inline bool add_ground_axiom(
		theory<Formula, natural_deduction<Formula, true>>& T,
		const atom<Negated, 1> consequent_atom,
		unsigned int constant, nd_step<Formula, true>* axiom)
{
	return T.add_unary_axiom<Negated>(consequent_atom.predicate, constant, axiom);
}

template<bool Negated, typename Formula>
inline bool add_ground_axiom(
		theory<Formula, natural_deduction<Formula, true>>& T,
		const atom<Negated, 2> consequent_atom,
		unsigned int constant, nd_step<Formula, true>* axiom)
{
	relation rel = { consequent_atom.predicate,
			(consequent_atom.args[0] == 0 ? constant : consequent_atom.args[0]),
			(consequent_atom.args[1] == 0 ? constant : consequent_atom.args[1]) };
	return T.add_binary_axiom<Negated>(rel, axiom);
}

template<typename Formula,
	bool Negated, unsigned int Arity,
	typename TheoryPrior, typename AxiomPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior>
bool do_mh_universal_elim(
		theory<Formula, natural_deduction<Formula, true>>& T,
		const proof_transformations<Formula>& proposed_proofs,
		unsigned int axiom_index,
		const hash_map<unsigned int, nd_step<Formula, true>*>& new_grounded_axioms,
		const atom<Negated, Arity> consequent_atom,
		double log_proposal_probability_ratio,
		double log_proof_stop_probability,
		double log_proof_continue_probability,
		TheoryPrior& theory_prior, AxiomPrior& axiom_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior)
{
	/* compute the proof portion of the prior for both current and proposed theories */
	for (const auto& entry : proposed_proofs.transformed_proofs) {
		/* contribution from the old proof */
		log_proposal_probability_ratio -= log_probability(
				*entry.key, log_proof_stop_probability,
				log_proof_continue_probability, axiom_prior,
				universal_introduction_prior, universal_elimination_prior);

		/* contribution from the new proof */
		log_proposal_probability_ratio += log_probability(
				*entry.key, log_proof_stop_probability,
				log_proof_continue_probability, axiom_prior,
				universal_introduction_prior, universal_elimination_prior,
				entry.value);
	}

	/* compute the prior of the current and proposed theories and add them to `log_proposal_probability_ratio` */
	log_proposal_probability_ratio += log_probability_ratio(T, theory_prior, axiom_index, new_grounded_axioms);

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we've accepted the proposal */
		if (!transform_proofs(proposed_proofs)) return false;

		T.remove_universal_quantification(axiom_index);
		for (auto entry : new_grounded_axioms)
			if (!add_ground_axiom(T, consequent_atom, entry.key, entry.value)) return false;
	}
	return true;
}

template<typename Formula, typename AxiomPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior>
bool do_mh_step(
		const theory<Formula, natural_deduction<Formula, true>>& T,
		double log_proof_stop_probability,
		double log_proof_continue_probability
		AxiomPrior& axiom_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, true> ProofCalculus;

	double log_proposal_probability_ratio = 0.0;

	/* TODO: select an axiom from `T` uniformly at random */
	unsigned int axiom_count = T.ground_axiom_count + T.universal_quantifications.length;
	unsigned int random = sample_uniform(axiom_count);
	if (!log_cache<double>::instance().ensure_size(axiom_count)) return false;
	log_proposal_probability_ratio -= -log_cache<double>::instance().get(axiom_count);
	if (random < T.ground_axiom_count) {
		/* we've selected a grounded axiom */
		for (auto entry : T.ground_concepts) {
			concept<ProofCalculus>& c = entry.value;
			if (random < c.types.size) {
				atom<false, 1> a = {c.types[random], {0}};
				return propose_universal_intro(T, a, log_proposal_probability_ratio);
			}
			random -= c.types.size;
			if (random < c.negated_types.size) {
				atom<true, 1> a = {c.negated_types[random], {0}};
				return propose_universal_intro(T, a, log_proposal_probability_ratio);
			}
			random -= r.negated_types.size;
			if (random < c.relations.size) {
				relation rel = c.relations[random];
				atom<false, 2> a = {rel.predicate, {rel.arg1, rel.arg2};
				return propose_universal_intro(T, a, log_proposal_probability_ratio);
			}
			random -= c.relations.size;
			if (random < c.negated_relations.size) {
				relation rel = c.relations[random];
				atom<true, 2> a = {rel.predicate, {rel.arg1, rel.arg2}};
				return propose_universal_intro(T, a, log_proposal_probability_ratio);
			}
			random -= c.negated_relations.size;
		}
	}
	random -= T.ground_axiom_count;

	if (random < T.universal_quantifications.length) {
		/* we've selected a universally-quantified formula */
		return propose_universal_elim(T, random,
				log_proposal_probability_ratio, log_proof_stop_probability,
				log_proof_continue_probability, theory_prior, axiom_prior,
				universal_introduction_prior, universal_elimination_prior);
	}

	fprintf(stderr, "propose ERROR: Unable to select axiom.\n");
	return false;
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
