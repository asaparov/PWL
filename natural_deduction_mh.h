#ifndef NATURAL_DEDUCTION_MH_H_
#define NATURAL_DEDUCTION_MH_H_

#include "theory.h"
#include "natural_deduction.h"

#include <core/random.h>

constexpr double LOG_2 = 0.693147180559945309417232121458176568075500134360255254120;

using namespace core;

template<typename Formula>
struct proof_transformation {
	array_map<nd_step<Formula, true>*, nd_step<Formula, true>*> transformations;

	static inline void free(proof_transformation<Formula>& t) {
		for (auto entry : t.transformations) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}
		core::free(t.transformations);
	}
};

template<typename Formula>
inline bool init(proof_transformation<Formula>& t) {
	return array_map_init(t.transformations, 4);
}

template<typename Formula>
struct proof_transformations {
	array_map<nd_step<Formula, true>*, proof_transformation<Formula>> transformed_proofs;

	static inline void free(proof_transformations<Formula>& t) {
		for (auto entry : t.transformed_proofs)
			core::free(entry.value);
		core::free(t.transformed_proofs);
	}
};

bool propose_transformation(
		const theory<Formula, nd_step<Formula, true>>& T,
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
			if (!proposed_proofs.transformed_proofs.check_size())
				return false;

			unsigned int index = proposed_proofs.transformed_proofs.index_of(step);
			if (index == proposed_proofs.transformed_proofs.size) {
				proposed_proofs.transformed_proofs.keys[index] = step;
				if (!init(proposed_proofs.transformed_proofs.values[index]))
					return false;
			}
			proposed_proofs.transformed_proofs.size++;

			proof_transformation<Formula>& value = proposed_proofs.transformed_proofs.values[index];
			value.transformations.put(old_step, new_step);
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
void select_axiom(const theory<Formula, nd_step<Formula, true>>& T,
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
		if (FirstSample) intersection.append(T.negated_relations.get(relation_value).value);
		else set_intersect(intersection, T.negated_relations.get(relation_value).value);
	}

	axiom_indices.remove(selected_index);
}

template<typename Formula, bool Negated>
bool get_axiom(
		array_map<unsigned int, nd_step<Formula, true>*>& types,
		unsigned int predicate, unsigned int arg,
		array<nd_step<Formula, true>*>& conjunct_steps)
{
	if (!types.ensure_capacity(types.size + 1))
		return false;
	unsigned int index = linear_search(types.keys, predicate, 0, types.size);
	if (index == types.size || types.keys[index] != predicate) {
		Proof* new_axiom = Negated ?
				Proof::new_axiom(Formula::new_not(Formula::new_atom(predicate, Formula::new_constant(arg)))) :
				Proof::new_axiom(Formula::new_atom(predicate, Formula::new_constant(arg)));
		if (new_axiom == NULL) return false;
		new_axiom->reference_count++;

		if (!conjunct_steps.add(new_axiom)) {
			free(*new_axiom);
			if (new_axiom->reference_count == 0)
				free(new_axiom);
		}

		shift_right(types.keys, types.size, index);
		shift_right(types.values, types.size, index);
		types.keys[index] = predicate;
		types.values[index] = new_axiom;
		types.size++;
	} else {
		if (!conjunct_steps.add(types.values[index]))
			return false;
	}
	return true;
}

template<typename Formula, bool Negated>
bool get_axiom(
		array_map<relation, nd_step<Formula, true>*>& relations,
		relation rel, unsigned int arg,
		array<nd_step<Formula, true>*>& conjunct_steps)
{
	typedef natural_deduction<Formula, true> Proof;

	if (!relations.ensure_capacity(relations.size + 1))
		return false;
	unsigned int index = linear_search(relations.keys, rel, 0, relations.size);
	if (index == relations.size || relations.keys[index] != rel) {
		Proof* atom = Formula::new_atom((rel.type == 0 ? arg : rel.type),
				Formula::new_constant(rel.arg1 == 0 ? arg : rel.arg1),
				Formula::new_constant(rel.arg2 == 0 ? arg : rel.arg2));
		if (atom == NULL) return false;
		atom->reference_count++;

		Proof* new_axiom = Negated ? Proof::new_axiom(Formula::new_not(atom)) : Proof::new_axiom(atom);
		if (new_axiom == NULL) {
			free(*atom); if (atom->reference_count == 0) free(atom);
			return false;
		}

		if (!conjunct_steps.add(new_axiom)) {
			free(*new_axiom);
			if (new_axiom->reference_count == 0)
				free(new_axiom);
		}

		shift_right(relations.keys, relations.size, index);
		shift_right(relations.values, relations.size, index);
		relations.keys[index] = rel;
		relations.values[index] = new_axiom;
		relations.size++;
	} else {
		if (!conjunct_steps.add(relations.values[index]))
			return false;
	}
	return true;
}

template<bool Negated, unsigned int Arity, unsigned int LiftedArgIndex>
struct atom {
	static_assert(LiftedArgIndex < Arity, "LiftedArgIndex is out of bounds");

	unsigned int predicate;
	unsigned int args[Arity];
};

template<typename Formula, unsigned int Index, bool Negated, unsigned int Arity, unsigned int LiftedArgIndex, typename... Formulas>
inline Formula* new_lifted_atom_helper(const atom<Negated, Arity, LiftedArgIndex> a, Formulas&&... args) {
	if (Index == 0)
		return Formula::new_atom(a.predicate, std::forward<Formulas>(args)...);
	else return new_lifted_atom_helper<Formula, Index - 1>(a,
			(Index - 1 == LiftedArgIndex ? Formula::new_variable(1) : Formula::new_constant(a.args[Index - 1])), args);
}

template<typename Formula, bool Negated, unsigned int Arity, unsigned int LiftedArgIndex>
inline Formula* new_lifted_atom(const atom<Negated, Arity, LiftedArgIndex> a) {
	return new_lifted_atom_helper<Formula, Arity>(a);
}

template<bool Negated, typename Formula>
inline const array<unsigned int>& get_negated_set(
		const theory<Formula, nd_step<Formula, true>>& T,
		const atom<Negated, 1, 0> a)
{
	return (Negated ? T.types.get(a.predicate).key : T.types.get(a.predicate).value);
}

template<bool Negated, unsigned int LiftedArgIndex, typename Formula>
const array<unsigned int>& get_negated_set(
		const theory<Formula, nd_step<Formula, true>>& T,
		const atom<Negated, 2, LiftedArgIndex> a)
{
	relation r = { a.predicate, (LiftedArgIndex == 0 ? 0 : a.args[0]), (LiftedArgIndex == 1 ? 0 : a.args[1]) };
	return (Negated ? T.relations.get(r).key : T.relations.get(r).value);
}

template<bool Negated, typename Formula>
const nd_step<Formula, true>* get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		const atom<Negated, 1, 0> a)
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

template<bool Negated, unsigned int LiftedArgIndex, typename Formula>
const nd_step<Formula, true>* get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		const atom<Negated, 2, LiftedArgIndex> a)
{
	array_map<relation, nd_step<Formula, true>*>& proofs = Negated ? instance.negated_relations : instance.relations;
	relation r = { a.predicate, (LiftedArgIndex == 0 ? 0 : a.args[0]), (LiftedArgIndex == 1 ? 0 : a.args[1]) };
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
inline void get_satisfying_concepts(
		const theory<Formula, nd_step<Formula, true>>& T,
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
		const theory<Formula, nd_step<Formula, true>>& T,
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

template<bool Negated, unsigned int Arity, unsigned int LiftedArgIndex, typename Formula>
bool propose_universal_intro(
		const theory<Formula, nd_step<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		const atom<Negated, Arity, LiftedArgIndex> a, double stop_probability)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula, true>> ConceptType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

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
	const ConceptType& c = T.ground_concepts.get(a.args[LiftedArgIndex]);
	array<pair<uint_fast8_t, unsigned int>> axiom_indices(
			c.types.length + c.negated_types.length + c.relations.length + c.negated_relations.length);
	for (unsigned int i = 0; i < c.types.length; i++)
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(0, i));
	for (unsigned int i = 0; i < c.negated_types.length; i++)
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(1, i));
	for (unsigned int i = 0; i < c.relations.length; i++)
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(2, i));
	for (unsigned int i = 0; i < c.negated_relations.length; i++)
		axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(3, i));

	if (axiom_indices.length == 0) {
		/* no candidate axioms are available */
		return true;
	}

	array<unsigned int> selected_types(8);
	array<unsigned int> selected_negated_types(8);
	array<relation> selected_relations(8);
	array<relation> selected_negated_relations(8);

	bool first_sample = true;
	unsigned int count = sample_geometric(stop_probability);
	array<unsigned int> intersection(32);
	while (true) {
		/* select a candidate axiom at random, and compute the set of constants
		for which it's true and intersect it with `intersection` */
		if (first_sample) {
			select_axiom<true>(T, c, axiom_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		} else {
			select_axiom<false>(T, c, axiom_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		}
		if (axiom_indices.length == 0) {
			/* it is impossible to create a universally-quantified rule using
				only ground concepts and relations */
			return true;
		}

		if (!has_intersection(intersection.data, intersection.length, negated_set.data, negated_set.length))
			break;
		first_sample = false;
	}

	while (true) {
		if (selected_types.length + selected_negated_types.length
		+ selected_relations.length + selected_negated_relations.length >= count)
			break;

		/* select a candidate axiom at random */
		select_axiom<false>(T, c, axiom_indices, selected_types, selected_negated_types,
				selected_relations, selected_negated_relations, intersection);
		if (axiom_indices.length == 0)
			break;
	}

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
			Formula::new_and(conjuncts), new_lifted_atom<Formula>(a)));
	free_formulas(conjuncts);
	if (axiom == NULL) return false;

	Formula* canonicalized = canonicalize(axiom);
	free(*axiom); if (axiom->reference_count == 0) free(axiom);
	if (canonicalized == NULL)
		return false;

	nd_step<Formula, true>* axiom_step = ProofCalculus::new_axiom(canonicalized);
	if (axiom_step == NULL) {
		free(*canonicalized);
		if (canonicalized->reference_count == 0)
			free(canonicalized);
		return false;
	}

	array<nd_step<Formula, Canonical>*> conjunct_steps(antecedent->array.length);
	for (unsigned int concept : intersection) {
		const ConceptType& instance = T.ground_concepts.get(concept);

		for (unsigned int i = 0; i < antecedent->array.length; i++) {
			Formula* atom = antecedent->array.operands[i];
			if (atom->type == FormulaType::NOT) {
				atom = atom->unary.operand;
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, true>(instance.negated_types, atom->predicate, concept, conjunct_steps)) {
						free(*canonicalized);
						if (canonicalized->reference_count == 0)
							free(canonicalized);
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
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, concept, conjunct_steps)) {
						free(*canonicalized);
						if (canonicalized->reference_count == 0)
							free(canonicalized);
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); if (axiom_step->reference_count == 0) free(axiom_step);
						return false;
					}
				}
			} else {
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, false>(instance.types, atom->predicate, concept, conjunct_steps)) {
						free(*canonicalized);
						if (canonicalized->reference_count == 0)
							free(canonicalized);
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
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, concept, conjunct_steps)) {
						free(*canonicalized);
						if (canonicalized->reference_count == 0)
							free(canonicalized);
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
	free(*canonicalized);
	if (canonicalized->reference_count == 0)
		free(canonicalized);

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
		const theory<Formula, nd_step<Formula, true>>& T,
		Formula* consequent, unsigned int variable,
		typename Formula::Term constant,
		nd_step<Formula, true>** new_axioms,
		unsigned int i)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, true> ProofCalculus;

	if (new_axioms[i] != NULL) return true;

	/* first check if the grounded conjunct already exists as an axiom in the theory */
	const Formula* operand = consequent->array.operands[i];
	if (operand->type == FormulaType::ATOM) {
		if (operand->atom.arg2.type == TermType::NONE) {
			/* atom is unary */
			Proof* proof = c.types.get(operand->atom.predicate, contains);
			if (contains) {
				new_axioms[i] = proof;
				proof->reference_count++;
				return true;
			}
		} else {
			/* atom is binary */
			relation r = { operand->atom.predicate,
					(operand->atom.arg1.type == TermType::VARIABLE ? 0 : operand->atom.arg1.constant),
					(operand->atom.arg2.type == TermType::VARIABLE ? 0 : operand->atom.arg2.constant) };
			Proof* proof = c.relations.get(r, contains);
			if (contains) {
				new_axioms[i] = proof;
				proof->reference_count++;
				return true;
			}
		}
	} else if (operand->type == FormulaType::NOT && operand->unary.operand->type == FormulaType::ATOM) {
		operand = operand->unary.operand;
		if (operand->atom.arg2.type == TermType::NONE) {
			/* atom is unary */
			Proof* proof = c.negated_types.get(operand->atom.predicate, contains);
			if (contains) {
				new_axioms[i] = proof;
				proof->reference_count++;
				return true;
			}
		} else {
			/* atom is binary */
			relation r = { operand->atom.predicate,
					(operand->atom.arg1.type == TermType::VARIABLE ? 0 : operand->atom.arg1.constant),
					(operand->atom.arg2.type == TermType::VARIABLE ? 0 : operand->atom.arg2.constant) };
			Proof* proof = c.negated_relations.get(r, contains);
			if (contains) {
				new_axioms[i] = proof;
				proof->reference_count++;
				return true;
			}
		}
	} else {
		fprintf(stderr, "make_grounded_conjunct ERROR: Expected a conjunction of literals.\n");
		return false;
	}

	new_axioms[i] = ProofCalculus::new_axiom(
			substitute(*consequent->array.operands[i],
			Formula::new_variable(variable), constant));
	if (new_axioms[i] == NULL)
		return false;
	new_axioms[i]->reference_count++;
	return true;
}

template<typename Formula>
nd_step<Formula, true>* make_grounded_conjunction(
		const theory<Formula, nd_step<Formula, true>>& T,
		Formula* consequent, unsigned int variable,
		typename Formula::Term constant,
		nd_step<Formula, true>** new_axioms)
{
	typedef natural_deduction<Formula, true> ProofCalculus;

	bool contains;
	const concept<ProofCalculus>& c = T.ground_concepts.get(constant.constant, contains);
#if !defined(NDEBUG)
	if (!contains)
		fprintf(stderr, "make_grounded_conjunction WARNING: Theory does not contain a concept for ID %u.\n", constant.constant);
#endif
	for (unsigned int i = 0; i < consequent->array.length; i++) {
		if (!make_grounded_conjunct(T, consequent, variable, constant, new_axioms, i)) {
			free_proofs(new_axioms, i);
			return NULL;
		}
	}
	new_step = ProofCalculus::new_conjunction_intro(array_view(new_axioms, consequent->array.length));
	if (new_step == NULL) {
		free_proofs(new_axioms, consequent->array.length);
		return NULL;
	}
	new_step->reference_count++;
	return new_step;
}

template<typename Formula>
bool propose_universal_elim(
		const theory<Formula, nd_step<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		const nd_step<Formula, true>* axiom)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	/* first check that all uses of this axiom eliminate the universal quantification; otherwise, we can't remove it */
	if (T.observations.contains(axiom->formula)) return false;

	Formula* consequent = formula->quantifier.operand->binary.right;
	for (Proof* child : axiom->children) {
		if (child->type != nd_step_type::UNIVERSAL_ELIMINATION
		 || child->operands[1]->term.type == TermType::CONSTANT)
			return true;

		for (Proof* grandchild : child->children) {
			if (grandchild->type != nd_step_type::IMPLICATION_ELIMINATION)
				return true;

			Proof* new_step = NULL;
			/* TODO: there is a potential problem; if a new grounded atom is
					 proveable by another univerally-quantified axiom, the
					 currently implemented algorithm will make that grounded
					 atom as an axiom */
			Proof** new_axioms = (Proof**) calloc(consequent->array.length, sizeof(Proof*));
			if (new_axioms == NULL) return false;

			/* check if `grandchild` is observed */
			if (T.observations.contains(grandchild)) {
				new_step = make_grounded_conjunction(consequent,
						formula->quantifier.variable, child->operands[1]->term, new_axioms);
				if (new_step == NULL) {
					free_proofs(new_axioms);
					free(new_axioms); return false;
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
						new_step = make_grounded_conjunction(consequent,
								formula->quantifier.variable, child->operands[1]->term, new_axioms);
						if (new_step == NULL) {
							free_proofs(new_axioms);
							free(new_axioms); return false;
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
						if (!make_grounded_conjunct(T, consequent, formula->quantifier.variable, child->operands[1]->term, new_axioms, index)) {
							if (new_step != NULL) {
								free(*new_step); if (new_step->reference_count == 0) free(new_step);
							}
							free_proofs(new_axioms); free(new_axioms); return false;
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
	return true;
}

template<typename Formula>
bool propose(const theory<Formula, nd_step<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		bool new_antecedent_stop_probability,
		double& log_proposal_probability_ratio)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	/* TODO: select an axiom from `T` uniformly at random */
	unsigned int random = sample_uniform(T.type_count + T.relation_count + T.universal_quantifications.length);
	log_proposal_probability_ratio -= log(T.type_count + T.relation_count + T.universal_quantifications.length);
	if (random < T.type_count) {
		/* we've selected a formula of the form `t(c)` */
		for (const auto& entry : T.types) {
			if (random < entry.value.key.length) {
				atom<false, 1, 0> a;
				a.predicate = entry.key;
				a.args[0] = entry.value.key[random];
				return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
			}
			random -= entry.value.key.length;
			if (random < entry.value.value.length) {
				atom<true, 1, 0> a;
				a.predicate = entry.key;
				a.args[0] = entry.value.key[random];
				return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
			}
			random -= entry.value.value.length;
		}
	}
	random -= T.type_count;

	if (random < T.relation_count) {
		/* we've selected a formula of the form `r(a,b)` */
		for (const auto& entry : T.relations) {
			if (random < entry.value.key.length) {
				log_proposal_probability_ratio -= LOG_2;
				if (sample_uniform(2) == 0) {
					atom<false, 2, 0> a;
					a.predicate = axiom->atom.predicate;
					a.args[0] = axiom->atom.arg1;
					a.args[1] = axiom->atom.arg2;
					return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
				} else {
					atom<false, 2, 1> a;
					a.predicate = axiom->atom.predicate;
					a.args[0] = axiom->atom.arg1;
					a.args[1] = axiom->atom.arg2;
					return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
				}
			}
			random -= entry.value.key.length;
			if (random < entry.value.value.length) {
				log_proposal_probability_ratio -= LOG_2;
				if (sample_uniform(2) == 0) {
					atom<true, 2, 0> a;
					a.predicate = axiom->atom.predicate;
					a.args[0] = axiom->atom.arg1;
					a.args[1] = axiom->atom.arg2;
					return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
				} else {
					atom<false, 2, 2> a;
					a.predicate = axiom->atom.predicate;
					a.args[0] = axiom->atom.arg1;
					a.args[1] = axiom->atom.arg2;
					return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
				}
			}
			random -= entry.value.value.length;
		}
	}
	random -= T.relation_count;

	if (random < T.universal_quantifications.length) {
		/* we've selected a universally-quantified formula */
		return propose_universal_elim(T, proposed_proofs, T.universal_quantifications[random]->formula);
	}

	fprintf(stderr, "propose ERROR: Unable to select axiom.\n");
	return false;
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
