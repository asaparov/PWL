#ifndef NATURAL_DEDUCTION_MH_H_
#define NATURAL_DEDUCTION_MH_H_

#include "theory.h"
#include "natural_deduction.h"

#include <core/random.h>

using namespace core;

template<bool FirstSample, typename Formula>
void select_axiom(const theory<Formula, nd_step<Formula, true>>& T,
		const concept<natural_deduction<Formula, true>>& c,
		array<pair<uint_fast8_t, unsigned int>>& concept_indices,
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
	unsigned int selected_index = sample_uniform(concept_indices.length);
	switch (concept_indices[selected_index].key) {
	case 0:
		type_value = c.types.keys[concept_indices[selected_index].value];
		if (!selected_types.add(type_value)) return false;
		if (FirstSample) intersection.append(T.types.get(type_value).key);
		else set_intersect(intersection, T.types.get(type_value).key);
	case 1:
		type_value = c.negated_types.keys[concept_indices[selected_index].value];
		if (!selected_negated_types.add(type_value)) return false;
		if (FirstSample) intersection.append(T.types.get(type_value).value);
		else set_intersect(intersection, T.types.get(type_value).value);
	case 2:
		relation_value = c.relations.keys[concept_indices[selected_index].value];
		if (!selected_relations.add(relation_value)) return false;
		if (FirstSample) intersection.append(T.relations.get(relation_value).key);
		else set_intersect(intersection, T.relations.get(relation_value).key);
	case 3:
		relation_value = c.negated_relations.keys[concept_indices[selected_index].value];
		if (!selected_negated_relations.add(relation_value)) return false;
		if (FirstSample) intersection.append(T.negated_relations.get(relation_value).value);
		else set_intersect(intersection, T.negated_relations.get(relation_value).value);
	}

	concept_indices.remove(selected_index);
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
	if (!relations.ensure_capacity(relations.size + 1))
		return false;
	unsigned int index = linear_search(relations.keys, rel, 0, relations.size);
	if (index == relations.size || relations.keys[index] != rel) {
		Proof* atom = Formula::new_atom((rel.type == 0 ? arg : rel.type),
				Formula::new_constant(rel.arg1 == 0 ? arg : rel.arg1),
				Formula::new_constant(rel.arg2 == 0 ? arg : rel.arg2));
		if (atom == NULL) return false;

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

template<typename Formula, bool Negated>
bool propose_atom_generalization(
		const theory<Formula, nd_step<Formula, true>>& T,
		unsigned int constant, unsigned int type,
		array<pair<nd_step<Formula, true>*, nd_step<Formula, true>*>>& proposed_proofs,
		double stop_probability)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula, true>> ConceptType;
	typedef natural_deduction<Formula, true> Proof;

	/* compute the set of constants for which the selected axiom is true */
	const array<unsigned int>& proposed_set = (Negated ? T.types.get(type).value : T.types.get(type).key);

	/* find other ground axioms that are connected to this constant */
	const ConceptType& c = T.ground_concepts.get(constant);
	array<pair<uint_fast8_t, unsigned int>> concept_indices(
			c.types.length + c.negated_types.length + c.relations.length + c.negated_relations.length);
	for (unsigned int i = 0; i < c.types.length; i++)
		concept_indices.add(make_pair<uint_fast8_t, unsigned int>(0, i));
	for (unsigned int i = 0; i < c.negated_types.length; i++)
		concept_indices.add(make_pair<uint_fast8_t, unsigned int>(1, i));
	for (unsigned int i = 0; i < c.relations.length; i++)
		concept_indices.add(make_pair<uint_fast8_t, unsigned int>(2, i));
	for (unsigned int i = 0; i < c.negated_relations.length; i++)
		concept_indices.add(make_pair<uint_fast8_t, unsigned int>(3, i));

	if (concept_indices.length == 0) {
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
			select_axiom<true>(T, c, concept_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		} else {
			select_axiom<false>(T, c, concept_indices, selected_types, selected_negated_types,
					selected_relations, selected_negated_relations, intersection);
		}
		if (concept_indices.length == 0) {
			/* it is impossible to create a universally-quantified rule using
				only ground concepts and relations */
			return true;
		}

		if (is_subset(intersection.data, intersection.length, proposed_set.data, proposed_set.length))
			break;
		first_sample = false;
	}

	while (true) {
		if (selected_types.length + selected_negated_types.length
		  + selected_relations.length + selected_negated_relations.length >= count)
			break;

		/* select a candidate axiom at random */
		select_axiom<false>(T, c, concept_indices, selected_types, selected_negated_types,
				selected_relations, selected_negated_relations, intersection);
		if (concept_indices.length == 0)
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

	if (!T.universal_quantications.ensure_capacity(T.universal_quantications.length + 1))
		return false;

	Formula* axiom = Formula::new_for_all(1, Formula::new_if_then(
			Formula::new_and(conjuncts), Formula::new_atom(type, Formula::new_variable(1))));
	free_formulas(conjuncts);
	if (axiom == NULL) return false;

	Formula* canonicalized = canonicalize(*axiom->quantifier.operand->binary.left);
	if (canonicalized == NULL) {
		free(*axiom); free(axiom);
		return false;
	}

	nd_step<Formula, Canonical>* axiom_step = Proof::new_axiom(axiom);
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
						return false;
					}
				}
			}
		}
	}
	free(*canonicalized); free(canonicalized);

	T.universal_quantications.add(axiom);
	for (unsigned int k = 0; k < intersection.length; k++) {
		const unsigned int concept = intersection[k];
		ConceptType& instance = T.ground_concepts.get(concept);

		array_map<unsigned int, nd_step<Formula, Canonical>*>& proofs = Negated ? instance.negated_types : instance.types;
		unsigned int index = proofs.index_of(type);
#if !defined(NDEBUG)
		if (index == proofs.size)
			fprintf(stderr, "propose_atom_generalization WARNING: Theory is invalid.\n");
		else if (proofs.values[index]->type != nd_step_type::AXIOM)
			fprintf(stderr, "propose_atom_generalization WARNING: Expected an axiom.\n");
#endif

		Proof* implication = Proof::new_universal_elim(axiom_step, Formula::new_constant(concept));
		if (implication == NULL) {
			for (unsigned int i = k; k < conjunct_steps.length; k++) {
				free(*conjunct_steps[k]);
				if (conjunct_steps[k]->reference_count == 0)
					free(conjunct_steps[k]);
			}
			return false;
		}

		Proof* antecedent = Proof::new_conjunction_intro(conjunct_steps);
		if (antecedent == NULL) {
			free(*implication); free(implication);
			for (unsigned int i = k; k < conjunct_steps.length; k++) {
				free(*conjunct_steps[k]);
				if (conjunct_steps[k]->reference_count == 0)
					free(conjunct_steps[k]);
			}
			return false;
		}

		nd_step<Formula, Canonical>* step = proofs.values[index];
		step->type = nd_step_type::IMPLICATION_ELIMINATION;
		Formula* old_formula = step->formula;
		step->operands[0] = implication;
		step->operands[1] = antecedent;
		for (unsigned int i = 2; i < ND_OPERAND_COUNT; i++)
			step->operands[i] = NULL;
		step->operands[0]->children.add(step); /* NOTE: `children` is initialized with positive capacity */
		step->operands[1]->children.add(step);

		/* remove the now redundant axiom from the theory */
		shift_left(proofs.keys + index, proofs.size - index);
		shift_left(proofs.values + index, proofs.size - index);
		proofs.size--;
	}
}

template<typename Formula>
bool propose(const theory<Formula, nd_step<Formula, true>>& T,
		array<pair<nd_step<Formula, true>*, nd_step<Formula, true>*>>& proposed_proofs)
{
	typedef typename Formula::Type FormulaType;

	/* TODO: select an axiom from `T` uniformly at random */
	Formula* axiom;

	if (axiom->type == FormulaType::FOR_ALL) {
		return propose_universal_elimination(T, axiom, proposed_proofs);
	} else if (axiom->type == FormulaType::ATOM) {
		return propose_atom_generalization(T, axiom, proposed_proofs);
	} else if (axiom->type == FormulaType::NOT) {
		if (axiom->unary.operand->type == FormulaType::ATOM)
			return propose_atom_generalization(T, axiom, proposed_proofs);
		else if (axiom->unary.operand.type == FormulaType::EXISTS)	
			return propose_exists_generalization(T, axiom, proposed_proofs);
	}
	fprintf(stderr, "propose ERROR: Selected an axiom with unsupported form.\n");
	return false;
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
