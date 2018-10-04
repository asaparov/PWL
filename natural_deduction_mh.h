#ifndef NATURAL_DEDUCTION_MH_H_
#define NATURAL_DEDUCTION_MH_H_

#include "theory.h"
#include "natural_deduction.h"

#include <core/random.h>

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

template<typename Formula, unsigned int Index, bool Negated, unsigned int Arity, unsigned int LiftedArgIndex>
inline Formula* new_lifted_atom_helper(atom<Negated, Arity, LiftedArgIndex> a, Formula*... args) {
	if (Index == 0)
		return Formula::new_atom(a.predicate, args...);
	else return new_lifted_atom_helper<Formula, Index - 1>(a,
			(Index - 1 == LiftedArgIndex ? Formula::new_variable(1) : Formula::new_constant(a.args[Index - 1])), args);
}

template<typename Formula, bool Negated, unsigned int Arity, unsigned int LiftedArgIndex>
inline Formula* new_lifted_atom(atom<Negated, Arity, LiftedArgIndex> a) {
	return new_lifted_atom_helper<Formula, Arity>(a);
}

template<bool Negated, typename Formula>
inline const array<unsigned int>& get_negated_set(
		const theory<Formula, nd_step<Formula, true>>& T,
		atom<Negated, 1, 0> a)
{
	return (Negated ? T.types.get(a.predicate).key : T.types.get(a.predicate).value);
}

template<bool Negated, unsigned int LiftedArgIndex, typename Formula>
const array<unsigned int>& get_negated_set(
		const theory<Formula, nd_step<Formula, true>>& T,
		atom<Negated, 2, LiftedArgIndex> a)
{
	relation r = { a.predicate, (LiftedArgIndex == 0 ? 0 : a.args[0]), (LiftedArgIndex == 1 ? 0 : a.args[1]) };
	return (Negated ? T.relations.get(r).key : T.relations.get(r).value);
}

template<bool Negated, typename Formula>
const proof<ProofCalculus>& get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		atom<Negated, 1, 0> a)
{
	typedef natural_deduction<Formula, true> ProofCalculus;

	array_map<unsigned int, proof<ProofCalculus>>& proofs = Negated ? instance.negated_types : instance.types;
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
const proof<ProofCalculus>& get_proof(
		const concept<natural_deduction<Formula, true>>& c,
		atom<Negated, 2, LiftedArgIndex> a)
{
	typedef natural_deduction<Formula, true> ProofCalculus;

	array_map<relation, proof<ProofCalculus>>& proofs = Negated ? instance.negated_relations : instance.relations;
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

template<bool Negated, unsigned int Arity, unsigned int LiftedArgIndex, typename Formula>
bool propose_universal_intro(
		const theory<Formula, nd_step<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		atom<Negated, Arity, LiftedArgIndex> a, double stop_probability)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula, true>> ConceptType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

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
	array<unsigned int> intersection(32);
	unsigned int count = sample_geometric(stop_probability);
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
	if (selected_types.length > 1)
		sort(selected_types);
	if (selected_negated_types.length > 1)
		sort(selected_negated_types);
	if (selected_relations.length > 1)
		sort(selected_relations);
	if (selected_negated_relations.length > 1)
		sort(selected_negated_relations);
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

	Formula* canonicalized = canonicalize(*axiom->quantifier.operand->binary.left);
	if (canonicalized == NULL) {
		free(*axiom); free(axiom);
		return false;
	}

	nd_step<Formula, Canonical>* axiom_step = ProofCalculus::new_axiom(axiom);
	free(*axiom) if (axiom->reference_count == 0) free(axiom);
	if (axiom_step == NULL) {
		free(*canonicalized); free(canonicalized);
		return false;
	}

	array<nd_step<Formula, Canonical>*> conjunct_steps(canonicalized->array.length);
	for (unsigned int concept : intersection) {
		const ConceptType& instance = T.ground_concepts.get(concept);

		for (unsigned int i = 0; i < canonicalized->array.length; i++) {
			Formula* atom = canonicalized->array.operands[i];
			if (atom->type == FormulaType::NOT) {
				atom = atom->unary.operand;
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, true>(instance.negated_types, atom->predicate, concept, conjunct_steps)) {
						free(*canonicalized); free(canonicalized);
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); free(axiom_step);
						return false;
					}
				} else {
					const relation rel = { atom->predicate,
							atom->arg1.type == TermType::CONSTANT ? atom->arg1.constant : 0,
							atom->arg2.type == TermType::CONSTANT ? atom->arg2.constant : 0 };
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, concept, conjunct_steps)) {
						free(*canonicalized); free(canonicalized);
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); free(axiom_step);
						return false;
					}
				}
			} else {
				if (atom->atom.arg2.type == TermType::NONE) {
					if (!get_axiom<Formula, false>(instance.types, atom->predicate, concept, conjunct_steps)) {
						free(*canonicalized); free(canonicalized);
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); free(axiom_step);
						return false;
					}
				} else {
					const relation rel = { atom->predicate,
							atom->arg1.type == TermType::CONSTANT ? atom->arg1.constant : 0,
							atom->arg2.type == TermType::CONSTANT ? atom->arg2.constant : 0 };
					if (!get_axiom<Formula, true>(instance.negated_relations, rel, concept, conjunct_steps)) {
						free(*canonicalized); free(canonicalized);
						for (auto step : conjunct_steps) {
							free(*step); if (step->reference_count == 0) free(step);
						}
						free(*axiom_step); free(axiom_step);
						return false;
					}
				}
			}
		}
	}
	free(*canonicalized); free(canonicalized);

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
			free(*axiom_step); free(axiom_step);
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
			free(*axiom_step); free(axiom_step);
			return false;
		}
	}

	return true;
}

void free_proofs(nd_step<Formula, true>** proofs, unsigned int count) {
	for (unsigned int i = 0; i < count; i++) {
		if (proofs[i] == NULL) continue;
		free(*proofs[i]);
		if (proofs[i]->reference_count == 0)
			free(proofs[i]);
	}
}

template<typename Formula>
nd_step<Formula, true>* make_grounded_conjunction(
		Formula* consequent, unsigned int variable,
		typename Formula::Term constant,
		nd_step<Formula, true>** new_axioms)
{
	typedef natural_deduction<Formula, true> ProofCalculus;

	for (unsigned int i = 0; i < consequent->array.length; i++) {
		if (new_axioms[i] != NULL) continue;

		new_axioms[i] = ProofCalculus::new_axiom(
				substitute(*consequent->array.operands[i],
				Formula::new_variable(variable), constant));
		if (new_axioms[i] == NULL) {
			free_proofs(new_axioms);
			return NULL;
		}
		new_axioms[i]->reference_count++;
	}
	new_step = ProofCalculus::new_conjunction_intro(array_view(new_axioms, consequent->array.length));
	if (new_step == NULL) {
		free_proofs(new_axioms);
		return NULL;
	}
	new_step->reference_count++;
	return new_step;
}

template<typename Formula>
bool propose_universal_elim(
		const theory<Formula, nd_step<Formula, true>>& T,
		proof_transformations<Formula>& proposed_proofs,
		const Formula* formula)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, true> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	/* find the proof step for the universally-quantified axiom */
	unsigned int i = 0;
	for (; i < T.universal_quantifications.length; i++)
		if (T.universal_quantifications[i].axiom->formula == formula) break;
#if !defined(NDEBUG)
	if (i == T.universal_quantifications.length || T.universal_quantifications[i].axiom->type != nd_step_type::AXIOM) {
		fprintf(stderr, "propose_universal_elim ERROR: Unable to find axiom in theory.\n");
		return false;
	}
#endif

	/* first check that all uses of this axiom eliminate the universal quantification; otherwise, we can't remove it */
	const proof<ProofCalculus>& pf = T.universal_quantifications[i];
	for (const Proof* proof : pf.proofs)
		if (proof->type == nd_step_type::AXIOM && proof->formula == formula) return false;

	Formula* consequent = formula->quantifier.operand->binary.right;
	for (Proof* child : pf.axiom->children) {
		if (child->type != nd_step_type::UNIVERSAL_ELIMINATION
		 || child->operands[1]->term.type == TermType::CONSTANT)
			return true;

		for (Proof* grandchild : child->children) {
			if (grandchild->type != nd_step_type::IMPLICATION_ELIMINATION)
				return true;

			Proof* new_step = NULL;
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
						if (new_axioms[index] != NULL) continue;

						new_axioms[index] = ProofCalculus::new_axiom(
								substitute(*consequent->array.operands[index],
								Formula::new_variable(variable), child->operands[1]->term));
						if (new_axioms[index] == NULL) {
							if (new_step != NULL) {
								free(*new_step); if (new_step->reference_count == 0) free(new_step);
							}
							free_proofs(new_axioms); free(new_axioms); return false;
						}
						new_axioms[index]->reference_count++;
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
		array<pair<nd_step<Formula, true>*, nd_step<Formula, true>*>>& proposed_proofs,
		bool new_antecedent_stop_probability)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	/* TODO: select an axiom from `T` uniformly at random */
	Formula* axiom;

	if (axiom->type == FormulaType::FOR_ALL) {
		return propose_universal_elim(T, axiom->quantifier.operand, proposed_proofs);
	} else if (axiom->type == FormulaType::ATOM) {
		if (axiom->atom.arg2.type == TermType::NONE) {
			atom<false, 1, 0> a;
			a.predicate = axiom->atom.predicate;
			a.args[0] = axiom->atom.arg1;
			return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
		} else if (sample_uniform(2) == 0) {
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
	} else if (axiom->type == FormulaType::NOT) {
		axiom = axiom->unary.operand;
		if (axiom->type == FormulaType::ATOM) {
			if (axiom->atom.arg2.type == TermType::NONE) {
				atom<true, 1, 0> a;
				a.predicate = axiom->atom.predicate;
				a.args[0] = axiom->atom.arg1;
				return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
			} else if (sample_uniform(2) == 0) {
				atom<true, 2, 0> a;
				a.predicate = axiom->atom.predicate;
				a.args[0] = axiom->atom.arg1;
				a.args[1] = axiom->atom.arg2;
				return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
			} else {
				atom<true, 2, 1> a;
				a.predicate = axiom->atom.predicate;
				a.args[0] = axiom->atom.arg1;
				a.args[1] = axiom->atom.arg2;
				return propose_universal_intro(T, proposed_proofs, a, new_antecedent_stop_probability);
			}
		}
	}
	fprintf(stderr, "propose ERROR: Selected an axiom with unsupported form.\n");
	return false;
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
