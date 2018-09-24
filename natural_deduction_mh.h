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
bool propose_atom_generalization(
		const theory<Formula, nd_step<Formula, true>>& T,
		unsigned int constant, unsigned int type,
		array<pair<nd_step<Formula, true>*, nd_step<Formula, true>*>>& proposed_proofs,
		double stop_probability)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef concept<natural_deduction<Formula>> ConceptType;
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
	if (axiom_step == NULL) {
		free(*canonicalized); free(canonicalized);
		free(*axiom) if (axiom->reference_count == 0) free(axiom);
		return false;
	}

	array<nd_step<Formula, Canonical>*> conjunct_steps(canonicalized->array.length);
	for (unsigned int concept : intersection) {
		const ConceptType& instance = T.ground_concepts.get(concept);

		array_map<unsigned int, nd_step<Formula, Canonical>*>& proofs = Negated ? instance.negated_types : instance.types;
		unsigned int index = proofs.index_of(type);
#if !defined(NDEBUG)
		if (index == proofs.size)
			fprintf(stderr, "propose_atom_generalization WARNING: Theory is invalid.\n");
		else if (proofs.values[index]->type != nd_step_type::AXIOM)
			fprintf(stderr, "propose_atom_generalization WARNING: Expected an axiom.\n");
#endif

		bool contains;
		for (unsigned int i = 0; i < canonicalized->array.length; i++) {
			Formula* atom = canonicalized->array.operands[i];
			if (atom->type == FormulaType::NOT) {
				atom = atom->unary.operand;
				if (atom->atom.arg2.type == TermType::NONE) {
					Proof* axiom = instance.negated_types.get(atom->atom.predicate, contains);
					/* TODO: implement this */
				} else {
					/* TODO: implement this */
				}
			} else {
				if (atom->atom.arg2.type == TermType::NONE) {
					/* TODO: implement this */
				} else {
					/* TODO: implement this */
				}
			}
		}

		nd_step<Formula, Canonical>* step = proofs.values[index];
		step->type = nd_step_type::IMPLICATION_ELIMINATION;
		free(*step->formula);
		if (step->formula->reference_count == 0)
			free(step->formula);
		step->operands[0] = Proof::new_universal_elim(axiom_step, Formula::new_constant(concept));
		step->operands[1] = Proof::new_conjunction_intro(conjunct_steps);
		for (unsigned int i = 2; i < ND_OPERAND_COUNT; i++)
			step->operands[i] = NULL;
		if (!step->operands[0]->children.add(step) || !step->operands[1]->children.add(step)) {
			/* TODO: implement this */
		}
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
