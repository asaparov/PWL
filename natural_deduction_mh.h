#ifndef NATURAL_DEDUCTION_MH_H_
#define NATURAL_DEDUCTION_MH_H_

#include "natural_deduction.h"
#include "theory.h"

#include <core/random.h>
#include <math/log.h>

constexpr double LOG_2 = 0.693147180559945309417232121458176568075500134360255254120;

using namespace core;

template<typename Proof>
constexpr inline bool accept(const hash_set<Proof*>& proofs, const array<typename Proof::FormulaType&>& extra_axioms, double proof_prior_diff) {
	return true;
}

template<typename T, bool AutomaticallyFree, bool OtherAutomaticallyFree>
inline bool add(hash_multiset<T, AutomaticallyFree>& multiset, const array_multiset<T, OtherAutomaticallyFree>& items)
{
	return multiset.add(items);
}

template<typename T, bool AutomaticallyFree, bool OtherAutomaticallyFree>
void subtract(hash_multiset<T, AutomaticallyFree>& multiset, const array_multiset<T, OtherAutomaticallyFree>& items)
{
	for (unsigned int i = 0; i < items.counts.size; i++) {
		bool contains; unsigned int bucket;
		unsigned int& count = multiset.counts.get(items.counts.keys[i], contains, bucket);
#if !defined(NDEBUG)
		if (count < items.counts.values[i]) {
			fprintf(stderr, "subtract WARNING: Attempted to remove more items from a bin than it contains.\n");
			count = 0;
		} else count -= items.counts.values[i];
#else
		count -= items.counts.values[i];
#endif
		if (count == 0)
			multiset.counts.remove_at(bucket);
	}
	multiset.sum -= items.sum;
}

template<typename Formula>
struct proof_transformations {
	array_map<nd_step<Formula>*, proof_substitution<Formula>> transformed_proofs;

	static inline void free(proof_transformations<Formula>& t) {
		t.free();
		core::free(t.transformed_proofs);
	}

private:
	proof_transformations() { } /* delete default constructor */

	inline void free() {
		for (auto entry : transformed_proofs)
			core::free(entry.value);
	}
};

template<typename Formula>
inline bool init(proof_transformations<Formula>& transformations) {
	return array_map_init(transformations.transformed_proofs, 16);
}

template<typename Formula, typename Canonicalizer>
bool propose_transformation(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		proof_transformations<Formula>& proposed_proofs,
		array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
		nd_step<Formula>* old_step, nd_step<Formula>* new_step)
{
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	array<Proof*> stack(16);
	hash_set<Proof*> visited(32);
	stack[0] = old_step; stack.length++;
	while (stack.length > 0) {
		Proof* step = stack.pop();
		if (!visited.add(step)) return false;

		if (T.observations.contains(step)) {
			/* `step` is an observed logical form in the theory `T` */
			if (!proposed_proofs.transformed_proofs.ensure_capacity(proposed_proofs.transformed_proofs.size + 1))
				return false;

			unsigned int index = proposed_proofs.transformed_proofs.index_of(step);
			if (index == proposed_proofs.transformed_proofs.size) {
				proposed_proofs.transformed_proofs.keys[index] = step;
				if (!init(proposed_proofs.transformed_proofs.values[index]))
					return false;
				proposed_proofs.transformed_proofs.size++;
			}

			proof_substitution<Formula>& value = proposed_proofs.transformed_proofs.values[index];
			value.map.put(old_step, new_step);
			new_step->reference_count++;
		}

		if (!stack.ensure_capacity(stack.length + step->children.length))
			return false;
		for (Proof* child : step->children) {
			if (visited.contains(child)) continue;
			stack[stack.length++] = child;
		}
	}

	/* check if `old_step` was an observation; if so, add `new_step` as a potential observation */
	if (T.observations.contains(old_step))
		return observation_changes.add({old_step, new_step});
	return true;
}

template<bool FirstSample, typename Formula, typename Canonicalizer>
bool select_axiom(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		const concept<natural_deduction<Formula>>& c,
		array<pair<uint_fast8_t, unsigned int>>& axiom_indices,
		array<typename Formula::Term*>& selected_types,
		array<typename Formula::Term*>& selected_negated_types,
		array<relation>& selected_relations, array<relation>& selected_negated_relations,
		array<unsigned int>& intersection)
{
	typedef typename Formula::Term Term;

	relation relation_value;
	/* TODO: This should be more likely if the resulting intersecting subset is
	   closer to the set of all concepts with the given `type`. `A` being
	   "closer" to being a subset of `B` can be measured in different ways:
		(1) if `A_1` has fewer elements not in `B` than `A_2` (i.e.
		    |A_1 \ B| < |A_2 \ B|), then `A_1`is closer than `A_2`,
		(2) if `A_1` has more elements in `B` than `A_2` (i.e.
			|A_1 intersect B| > |A_2 intersect B|), then `A_1` is closer than
			`A_2`. */
	unsigned int selected_index = sample_uniform(axiom_indices.length);
	if (axiom_indices[selected_index].key == 0) {
		Term& type_value = c.types.keys[axiom_indices[selected_index].value];
		if (!selected_types.add(&type_value)) return false;
		const array<unsigned int>& ids = T.atoms.get(type_value).key;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 1) {
		Term& type_value = c.negated_types.keys[axiom_indices[selected_index].value];
		if (!selected_negated_types.add(&type_value)) return false;
		const array<unsigned int>& ids = T.atoms.get(type_value).value;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 2) {
		relation_value = c.relations.keys[axiom_indices[selected_index].value];
		if (!selected_relations.add(relation_value)) return false;
		const array<unsigned int>& ids = T.relations.get(relation_value).key;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 3) {
		relation_value = c.negated_relations.keys[axiom_indices[selected_index].value];
		if (!selected_negated_relations.add(relation_value)) return false;
		const array<unsigned int>& ids = T.relations.get(relation_value).value;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	}

	axiom_indices.remove(selected_index);
	return true;
}

template<typename Formula, bool Negated>
inline bool get_axiom(
		const array_map<typename Formula::Term, nd_step<Formula>*>& types,
		const typename Formula::Term& atom,
		array<nd_step<Formula>*>& conjunct_steps)
{
#if !defined(NDEBUG)
	if (!types.contains(atom))
		fprintf(stderr, "get_axiom WARNING: The given `atom` does not exist in the map `types`.\n");
#endif
	return conjunct_steps.add(types.get(atom));
}

template<typename Formula, bool Negated>
inline bool get_axiom(
		const array_map<relation, nd_step<Formula>*>& relations,
		relation rel, array<nd_step<Formula>*>& conjunct_steps)
{
#if !defined(NDEBUG)
	if (!relations.contains(rel))
		fprintf(stderr, "get_axiom WARNING: The given `rel` does not exist in the map `relations`.\n");
#endif
	return conjunct_steps.add(relations.get(rel));
}

template<bool Negated, typename Formula, typename Canonicalizer>
inline const array<unsigned int>& get_concept_set(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		const typename Formula::Term& a)
{
	typedef typename Formula::TermType TermType;
	if (a.type == TermType::UNARY_APPLICATION) {
		return (Negated ? T.atoms.get(a).value : T.atoms.get(a).key);
	} else {
		relation r = { a.ternary.first->constant,
				a.ternary.second->type == TermType::VARIABLE ? 0 : a.ternary.second->constant,
				a.ternary.third->type == TermType::VARIABLE ? 0 : a.ternary.third->constant };
		return (Negated ? T.relations.get(r).value : T.relations.get(r).key);
	}
}

template<bool Negated, typename Formula>
nd_step<Formula>* get_proof(
		concept<natural_deduction<Formula>>& c,
		const typename Formula::Term& a)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	if (a.type == TermType::UNARY_APPLICATION) {
		array_map<Term, nd_step<Formula>*>& proofs = Negated ? c.negated_types : c.types;
		unsigned int index = proofs.index_of(a);
#if !defined(NDEBUG)
		if (index == proofs.size) {
			fprintf(stderr, "get_proof WARNING: Theory is invalid.\n");
			return NULL;
		}
#endif
		return proofs.values[index];
	} else {
		array_map<relation, nd_step<Formula>*>& proofs = Negated ? c.negated_relations : c.relations;
		relation r = { a.ternary.first->constant,
				a.ternary.second->type == TermType::VARIABLE ? 0 : a.ternary.second->constant,
				a.ternary.third->type == TermType::VARIABLE ? 0 : a.ternary.third->constant };
		unsigned int index = proofs.index_of(r);
#if !defined(NDEBUG)
		if (index == proofs.size) {
			fprintf(stderr, "get_proof WARNING: Theory is invalid.\n");
			return NULL;
		}
#endif
		return proofs.values[index];
	}
}

template<bool FirstSample, bool Negated,
	typename Formula, typename Canonicalizer>
inline void get_satisfying_concepts_helper(
		const theory<natural_deduction<Formula>, Canonicalizer>& T, const typename Formula::Term* predicate,
		const typename Formula::Term* arg1, const typename Formula::Term* arg2, array<unsigned int>& intersection)
{
	typedef typename Formula::TermType TermType;

	if (arg2 == NULL) {
		/* this is a unary atom */
		const pair<array<unsigned int>, array<unsigned int>>& list = T.atoms.get(*predicate);
		const array<unsigned int>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist.data, sublist.length);
		else set_intersect(intersection, sublist);
	} else {
		/* this is a binary atom */
		relation r = { predicate,
				(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
				(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
		const pair<array<unsigned int>, array<unsigned int>>& list = T.relations.get(r);
		const array<unsigned int>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist.data, sublist.length);
		else set_intersect(intersection, sublist);
	}
}

template<bool FirstSample, typename Formula, typename Canonicalizer>
bool get_satisfying_concepts(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		const Formula* formula, array<unsigned int>& intersection)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;

	Term const* predicate; Term const* arg1; Term const* arg2;
	if (is_atomic(*formula, predicate, arg1, arg2)) {
		get_satisfying_concepts_helper<FirstSample, false>(T, predicate, arg1, arg2, intersection);
	} else if (formula->type == FormulaType::NOT && is_atomic(*formula->unary.operand, predicate, arg1, arg2)) {
		get_satisfying_concepts_helper<FirstSample, true>(T, predicate, arg1, arg2, intersection);
	} else if (formula->type == FormulaType::AND) {
		for (unsigned int i = 0; i < formula->array.length; i++) {
			Formula* conjunct = formula->array.operands[i];
			if ((i == 0 && !get_satisfying_concepts<FirstSample>(T, conjunct, intersection))
			 || (i > 0 && !get_satisfying_concepts<false>(T, conjunct, intersection)))
				return false;
		}
	} else {
		fprintf(stderr, "get_satisfying_concepts ERROR: Expected a literal.\n");
		return false;
	}
	return true;
}

template<bool Negated, typename Formula>
struct universal_intro_proposal {
	typedef typename Formula::Term Term;

	unsigned int concept_id;
	nd_step<Formula>* old_axiom;
	nd_step<Formula>* new_axiom;
	Term& consequent_atom;
	Formula* new_antecedent_set_size_axiom;
	Formula* new_consequent_set_size_axiom;

	~universal_intro_proposal() {
		free(*old_axiom); if (old_axiom->reference_count == 0) free(old_axiom);
		free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
	}
};

template<bool Negated, typename Formula>
inline universal_intro_proposal<Negated, Formula> make_universal_intro_proposal(
		unsigned int concept_id,
		nd_step<Formula>* old_axiom,
		nd_step<Formula>* new_axiom,
		typename Formula::Term& consequent_atom,
		Formula* new_antecedent_set_size_axiom,
		Formula* new_consequent_set_size_axiom)
{
	old_axiom->reference_count++;
	new_axiom->reference_count++;
	return {concept_id, old_axiom, new_axiom, consequent_atom, new_antecedent_set_size_axiom, new_consequent_set_size_axiom};
}

template<typename Formula>
inline void free_formulas(array<Formula*>& formulas) {
	for (Formula* formula : formulas) {
		free(*formula);
		if (formula->reference_count == 0)
			free(formula);
	}
}

template<typename Formula>
bool get_conjunct_step(const Formula* atom,
		const concept<natural_deduction<Formula>>& instance,
		array<nd_step<Formula>*>& conjunct_steps)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	Term const* predicate = nullptr; Term const* arg1 = nullptr; Term const* arg2 = nullptr;
	if (atom->type == FormulaType::NOT) {
		atom = atom->unary.operand;
		is_atomic(*atom, predicate, arg1, arg2);
		if (arg2 == NULL) {
			if (!get_axiom<Formula, true>(instance.negated_types, *atom, conjunct_steps)) {
				for (auto step : conjunct_steps) {
					free(*step); if (step->reference_count == 0) free(step);
				}
				return false;
			}
		} else {
			if (predicate->type != TermType::CONSTANT)
				return false;
			const relation rel = { predicate->constant,
					arg1->type == TermType::CONSTANT ? arg1->constant : 0,
					arg2->type == TermType::CONSTANT ? arg2->constant : 0 };
			if (!get_axiom<Formula, true>(instance.negated_relations, rel, conjunct_steps)) {
				for (auto step : conjunct_steps) {
					free(*step); if (step->reference_count == 0) free(step);
				}
				return false;
			}
		}
	} else {
		is_atomic(*atom, predicate, arg1, arg2);
		if (arg2 == NULL) {
			if (!get_axiom<Formula, false>(instance.types, *atom, conjunct_steps)) {
				for (auto step : conjunct_steps) {
					free(*step); if (step->reference_count == 0) free(step);
				}
				return false;
			}
		} else {
			if (predicate->type != TermType::CONSTANT)
				return false;
			const relation rel = { predicate->constant,
					arg1->type == TermType::CONSTANT ? arg1->constant : 0,
					arg2->type == TermType::CONSTANT ? arg2->constant : 0 };
			if (!get_axiom<Formula, true>(instance.negated_relations, rel, conjunct_steps)) {
				for (auto step : conjunct_steps) {
					free(*step); if (step->reference_count == 0) free(step);
				}
				return false;
			}
		}
	}
	return true;
}

template<bool Negated, typename ProofCalculus, typename Formula>
bool intensional_element_of(
		const concept<ProofCalculus>& c, const Formula* set_formula)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef typename ProofCalculus::Proof Proof;

	const array_map<Term, Proof*>& types = (Negated ? c.negated_types : c.types);
	const array_map<relation, Proof*>& relations = (Negated ? c.negated_relations : c.relations);
	switch (set_formula->type) {
	case FormulaType::UNARY_APPLICATION:
		if (set_formula->binary.right->type == TermType::VARIABLE
		 && set_formula->binary.right->variable == 1
		 && set_formula->binary.left->type == TermType::CONSTANT) {
			return types.contains(*set_formula);
		} else {
			return false;
		}
	case FormulaType::BINARY_APPLICATION:
		if ((set_formula->ternary.first->type != TermType::CONSTANT && (set_formula->ternary.first->type != TermType::VARIABLE || set_formula->ternary.first->variable != 1))
		 || (set_formula->ternary.second->type != TermType::CONSTANT && (set_formula->ternary.second->type != TermType::VARIABLE || set_formula->ternary.second->variable != 1))
		 || (set_formula->ternary.third->type != TermType::CONSTANT && (set_formula->ternary.third->type != TermType::VARIABLE || set_formula->ternary.third->variable != 1)))
			return false;
		return relations.contains({
				set_formula->ternary.first->type == TermType::VARIABLE ? 0 : set_formula->ternary.first->constant,
				set_formula->ternary.second->type == TermType::VARIABLE ? 0 : set_formula->ternary.second->constant,
				set_formula->ternary.third->type == TermType::VARIABLE ? 0 : set_formula->ternary.third->constant });
	case FormulaType::NOT:
		return intensional_element_of<!Negated>(c, set_formula->unary.operand);
	case FormulaType::AND:
		for (unsigned int i = 0; i < set_formula->array.length; i++)
			if (!intensional_element_of<Negated>(c, set_formula->array.operands[i])) return false;
		return true;
	case FormulaType::OR:
		for (unsigned int i = 0; i < set_formula->array.length; i++)
			if (intensional_element_of<Negated>(c, set_formula->array.operands[i])) return true;
		return false;
	case FormulaType::IF_THEN:
		return intensional_element_of<Negated>(c, set_formula->binary.left)
			&& intensional_element_of<!Negated>(c, set_formula->binary.right);
	case FormulaType::TRUE:
		return !Negated;
	case FormulaType::FALSE:
		return Negated;
	case FormulaType::EXISTS:
		/* since the proposal that uses this function currently only supports conjunctions of atoms, we return false here */
		return false;
	case FormulaType::EQUALS:
	case FormulaType::IFF:
	case FormulaType::FOR_ALL:
	case FormulaType::LAMBDA:
		fprintf(stderr, "intensional_element_of ERROR: Not implemented.\n");
		return false;
	case FormulaType::VARIABLE:
	case FormulaType::CONSTANT:
	case FormulaType::PARAMETER:
	case FormulaType::NUMBER:
	case FormulaType::STRING:
	case FormulaType::UINT_LIST:
	case FormulaType::ANY:
	case FormulaType::ANY_RIGHT:
	case FormulaType::ANY_RIGHT_ONLY:
	case FormulaType::ANY_ARRAY:
	case FormulaType::ANY_CONSTANT:
	case FormulaType::ANY_CONSTANT_EXCEPT:
	case FormulaType::ANY_QUANTIFIER:
	case FormulaType::VARIABLE_PREIMAGE:
		break;
	}
	fprintf(stderr, "intensional_element_of ERROR: Unrecognized formula type.\n");
	return false;
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool is_satisfying_antecedent(
		unsigned int concept_id, unsigned int antecedent_id,
		const theory<ProofCalculus, Canonicalizer>& T)
{
	/* check if `concept_id` cannot belong to the set given by `antecedent_id` (if there is a contradiction) */
	/* first check if we can prove that the entry belongs to any descendant of `antecedent_id` */
	bool provably_consistent = false;
	array<unsigned int> stack(8);
	hash_set<unsigned int> visited(16);
	stack[stack.length++] = antecedent_id;
	visited.add(antecedent_id);
	while (stack.length > 0) {
		unsigned int descendant = stack.pop();
		if (intensional_element_of<false>(T.ground_concepts[concept_id - T.new_constant_offset], T.sets.sets[descendant].set_formula())) {
			provably_consistent = true;
			break;
		}

		for (unsigned int child : T.sets.intensional_graph.vertices[descendant].children) {
			if (visited.contains(child)) continue;
			if (!visited.add(child) || !stack.add(child)) return false;
		} /*for (const auto& entry : T.sets.extensional_graph.vertices[descendant].children) {
			if (visited.contains(entry.key)) continue;
			if (!visited.add(entry.key) || !stack.add(entry.key)) return false;
		}*/
	}

	return provably_consistent;
}

template<bool Negated, typename Formula, typename Canonicalizer>
bool get_satisfying_extensional_edges(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		const typename Formula::Term& a,
		unsigned int concept_id,
		array<pair<unsigned int, unsigned int>>& satisfying_extensional_edges)
{
	typedef typename Formula::Type FormulaType;

	for (unsigned int i = 1; i < T.sets.set_count + 1; i++) {
		if (T.sets.sets[i].size_axioms.data == nullptr) continue;
		Formula* set_formula = T.sets.sets[i].set_formula();
		if (Negated) {
			if (set_formula->type != FormulaType::NOT) continue;
			set_formula = set_formula->unary.operand;
		}
		if (a == *set_formula) {
			for (const auto& entry : T.sets.extensional_graph.vertices[i].children) {
				if (!is_satisfying_antecedent(concept_id, entry.key, T)) continue;

				/* as far as we can tell, we can't prove that `concept_id` cannot belong to the antecedent set */
				for (unsigned int j = 0; j < entry.value.length; j++) {
					Formula* formula = entry.value[j]->formula->quantifier.operand;
					while (formula->type == FormulaType::FOR_ALL)
						formula = formula->quantifier.operand;
					if (*T.sets.sets[entry.key].set_formula() != *formula->binary.left || a != *formula->binary.right) continue;
					if (!satisfying_extensional_edges.add({i, entry.key})) return false;
				}
			}
		}
	}
	return true;
}

constexpr double SET_SIZE_PROPOSAL_P = 0.1;
double SET_SIZE_PROPOSAL_LOG_P = log(SET_SIZE_PROPOSAL_P);
double SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P = log(1.0 - SET_SIZE_PROPOSAL_P);

struct unfixed_set_counter {
	int change;
	double log_probability;
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		unfixed_set_counter& counter)
{
	if (max_set_size == UINT_MAX) {
		out = min_set_size + sample_geometric(0.1);
		counter.log_probability += (out - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		out = min_set_size + sample_uniform(max_set_size - min_set_size + 1);
		if (!log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2))
			return false;
		counter.log_probability += -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
	counter.change++;
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const unfixed_set_counter& visitor)
{ }

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const unfixed_set_counter& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const unfixed_set_counter& visitor)
{ }

template<
	bool Negated, typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_universal_intro(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		typename Formula::Term& a,
		unsigned int concept_id,
		double& log_proposal_probability_ratio,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		const ProposalDistribution& proposal_distribution)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	/* look for existing universally-quantified formulas that imply this formula */
	array<pair<unsigned int, unsigned int>> satisfying_extensional_edges(16);
	if (!get_satisfying_extensional_edges<Negated>(T, a, concept_id, satisfying_extensional_edges)) return false;

	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) return false;

	Formula* antecedent; Formula* consequent;
	Formula* new_antecedent_set_size_axiom = nullptr;
	Formula* new_consequent_set_size_axiom = nullptr;
	unsigned int antecedent_set, consequent_set;
	nd_step<Formula>* axiom_step = nullptr;
	unsigned int random = sample_uniform(satisfying_extensional_edges.length + 1);
	log_proposal_probability_ratio -= -log_cache<double>::instance().get(satisfying_extensional_edges.length + 1);
	unfixed_set_counter unfixed_set_count_change = {0};
	if (random < satisfying_extensional_edges.length) {
		const pair<unsigned int, unsigned int>& selected = satisfying_extensional_edges[random];
		consequent_set = selected.key;
		antecedent_set = selected.value;
		consequent = T.sets.sets[consequent_set].set_formula();
		antecedent = T.sets.sets[antecedent_set].set_formula();
		array<nd_step<Formula>*> proofs = T.sets.extensional_graph.vertices[selected.key].children.get(selected.value);
		for (nd_step<Formula>* proof : proofs) {
			Formula* formula = proof->formula->quantifier.operand;
			if (formula->type == FormulaType::FOR_ALL)
				formula = formula->quantifier.operand;
			if (*antecedent == *formula->binary.left && *consequent == *formula->binary.right) {
				axiom_step = proof;
				break;
			}
		}
		axiom_step->reference_count++;

	} else {

		/* compute the set of constants for which the selected axiom is true */
		const array<unsigned int>& negated_set = get_concept_set<!Negated>(T, a);

		/* find other ground axioms that are connected to this constant */
		const concept<ProofCalculus>& c = T.ground_concepts[concept_id - T.new_constant_offset];
		array<pair<uint_fast8_t, unsigned int>> axiom_indices(
				c.types.size + c.negated_types.size + c.relations.size + c.negated_relations.size);
		for (unsigned int i = 0; i < c.types.size; i++) {
			if (!Negated && a.type == TermType::UNARY_APPLICATION && *c.types.keys[i].binary.left == *a.binary.left) continue; /* avoid selecting the consequent as a conjunct in the antecedent */
			axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(0, i));
		} for (unsigned int i = 0; i < c.negated_types.size; i++) {
			if (Negated && a.type == TermType::UNARY_APPLICATION && *c.negated_types.keys[i].binary.left == *a.binary.left) continue; /* avoid selecting the consequent as a conjunct in the antecedent */
			axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(1, i));
		} for (unsigned int i = 0; i < c.relations.size; i++) {
			/* avoid selecting the consequent as a conjunct in the antecedent */
			if (!Negated && a.type == TermType::BINARY_APPLICATION
			 && a.ternary.first->type == TermType::CONSTANT && c.relations.keys[i].predicate == a.ternary.first->constant
			 && ((c.relations.keys[i].arg1 == 0 && a.ternary.second->type == TermType::VARIABLE && a.ternary.second->variable == 0)
			  || (c.relations.keys[i].arg1 == a.ternary.second->constant && a.ternary.second->type == TermType::CONSTANT))
			 && ((c.relations.keys[i].arg2 == 0 && a.ternary.third->type == TermType::VARIABLE && a.ternary.third->variable == 0)
			  || (c.relations.keys[i].arg2 == a.ternary.third->constant && a.ternary.third->type == TermType::CONSTANT))) continue;
			axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(2, i));
		} for (unsigned int i = 0; i < c.negated_relations.size; i++) {
			/* avoid selecting the consequent as a conjunct in the antecedent */
			if (Negated && a.type == TermType::BINARY_APPLICATION
			 && a.ternary.first->type == TermType::CONSTANT && c.negated_relations.keys[i].predicate == a.ternary.first->constant
			 && ((c.negated_relations.keys[i].arg1 == 0 && a.ternary.second->type == TermType::VARIABLE && a.ternary.second->variable == 0)
			  || (c.negated_relations.keys[i].arg1 == a.ternary.second->constant && a.ternary.second->type == TermType::CONSTANT))
			 && ((c.negated_relations.keys[i].arg2 == 0 && a.ternary.third->type == TermType::VARIABLE && a.ternary.third->variable == 0)
			  || (c.negated_relations.keys[i].arg2 == a.ternary.third->constant && a.ternary.third->type == TermType::CONSTANT))) continue;
			axiom_indices.add(make_pair<uint_fast8_t, unsigned int>(3, i));
		}

		if (axiom_indices.length == 0) {
			/* no candidate axioms are available */
			free(proposed_proofs); return true;
		}

		array<Term*> selected_types(8);
		array<Term*> selected_negated_types(8);
		array<relation> selected_relations(8);
		array<relation> selected_negated_relations(8);

		array<unsigned int> intersection(32);
		unsigned int count = sample_uniform(axiom_indices.length) + 1;
		log_proposal_probability_ratio -= -log_cache<double>::instance().get(axiom_indices.length)
				- lgamma(axiom_indices.length + 1) + lgamma(count + 1) + lgamma(axiom_indices.length - count + 1);
		for (unsigned int i = 0; i < count; i++) {
			/* select a candidate axiom at random, and compute the set of constants
			for which it's true and intersect it with `intersection` */
			if (i == 0) {
				if (!select_axiom<true>(T, c, axiom_indices, selected_types, selected_negated_types,
						selected_relations, selected_negated_relations, intersection)) { free(proposed_proofs); return false; }
			} else {
				if (!select_axiom<false>(T, c, axiom_indices, selected_types, selected_negated_types,
						selected_relations, selected_negated_relations, intersection)) { free(proposed_proofs); return false; }
			}
		}

		if (has_intersection(intersection.data, intersection.length, negated_set.data, negated_set.length)) {
			free(proposed_proofs);
			return true;
		}

		/* we only consider concepts that have the axiom `a` */
		set_intersect(intersection, get_concept_set<Negated>(T, a));

		/* we've finished constructing the new universally-quantified formula,
		   and now we need to transform all the relevant proofs */
		array<Formula*> conjuncts(32);
		for (Term* type : selected_types) {
			Formula* conjunct = Formula::new_apply(type->binary.left, &Term::template variables<1>::value);
			if (conjunct == NULL || !conjuncts.add(conjunct)) {
				free_formulas(conjuncts); free(proposed_proofs); return false;
			}
			Term::template variables<1>::value.reference_count++;
			type->binary.left->reference_count++;
		}
		for (Term* type : selected_negated_types) {
			Formula* conjunct = Formula::new_not(Formula::new_apply(type->binary.left, &Term::template variables<1>::value));
			if (conjunct == NULL || !conjuncts.add(conjunct)) {
				free_formulas(conjuncts); free(proposed_proofs); return false;
			}
			Term::template variables<1>::value.reference_count++;
			type->binary.left->reference_count++;
		}
		for (const relation& r : selected_relations) {
			Formula* conjunct = Formula::new_atom(r.predicate,
					(r.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(r.arg1),
					(r.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(r.arg2));
			if (conjunct == NULL || !conjuncts.add(conjunct)) {
				free_formulas(conjuncts); free(proposed_proofs); return false;
			}
			if (r.arg1 == 0) Term::template variables<1>::value.reference_count++;
			if (r.arg2 == 0) Term::template variables<1>::value.reference_count++;
		}
		for (const relation& r : selected_negated_relations) {
			Formula* conjunct = Formula::new_not(Formula::new_atom(r.predicate,
					(r.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(r.arg1),
					(r.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(r.arg2)));
			if (conjunct == NULL || !conjuncts.add(conjunct)) {
				free_formulas(conjuncts); free(proposed_proofs); return false;
			}
			if (r.arg1 == 0) Term::template variables<1>::value.reference_count++;
			if (r.arg2 == 0) Term::template variables<1>::value.reference_count++;
		}

		Formula* axiom;
		if (conjuncts.length == 1) {
			axiom = Formula::new_for_all(1, Formula::new_if_then(conjuncts[0], (Negated ? Formula::new_not(&a) : &a)));
		} else {
			axiom = Formula::new_for_all(1, Formula::new_if_then(
					Formula::new_and(make_array_view(conjuncts.data, conjuncts.length)), (Negated ? Formula::new_not(&a) : &a)));
		}
		if (axiom == NULL) {
			free_formulas(conjuncts);
			free(proposed_proofs); return false;
		} else {
			for (Formula* formula : conjuncts)
				formula->reference_count++;
			free_formulas(conjuncts);
			a.reference_count++;
		}

		Formula* canonicalized = Canonicalizer::canonicalize(*axiom);
		free(*axiom); if (axiom->reference_count == 0) free(axiom);
		if (canonicalized == NULL) {
			free(proposed_proofs); return false;
		}

		bool is_antecedent_new, is_consequent_new;
		antecedent = canonicalized->quantifier.operand->binary.left;
		consequent = canonicalized->quantifier.operand->binary.right;
		axiom_step = T.template get_subset_axiom<false>(antecedent, consequent, 1, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, unfixed_set_count_change);
		free(*canonicalized);
		if (canonicalized->reference_count == 0)
			free(canonicalized);
		if (axiom_step == NULL) {
			free(proposed_proofs); return true; /* the proposed axiom is inconsistent */
		}
		axiom_step->reference_count++;

		antecedent = axiom_step->formula->quantifier.operand->binary.left;
		consequent = axiom_step->formula->quantifier.operand->binary.right;

		if (is_antecedent_new)
			new_antecedent_set_size_axiom = T.sets.sets[antecedent_set].size_axioms[0]->formula;
		if (is_consequent_new)
			new_consequent_set_size_axiom = T.sets.sets[consequent_set].size_axioms[0]->formula;
	}

	/* compute the probability of the inverse proposal */
	int delta_atom_count = 0;
	if (a.type == TermType::BINARY_APPLICATION && *a.ternary.second == *a.ternary.third) delta_atom_count -= 3;
	else if (a.type == TermType::BINARY_APPLICATION) delta_atom_count -= 2;
	else delta_atom_count -= 1;
	log_proposal_probability_ratio += log_probability(proposal_distribution, delta_atom_count, (satisfying_extensional_edges.length == 0 ? 1 : 0), unfixed_set_count_change.change);
	log_proposal_probability_ratio += unfixed_set_count_change.log_probability;

	unsigned int antecedent_length;
	if (antecedent->type != FormulaType::AND) {
		antecedent_length = 1;
	} else {
		antecedent_length = antecedent->array.length;
	}

	array<nd_step<Formula>*> conjunct_steps(antecedent_length);
	concept<ProofCalculus>& instance = T.ground_concepts[concept_id - T.new_constant_offset];

	if (antecedent->type != FormulaType::AND) {
		if (!get_conjunct_step(antecedent, instance, conjunct_steps)) {
			free(proposed_proofs); free(*axiom_step);
			if (axiom_step->reference_count == 1)
				T.free_subset_axiom(antecedent, consequent, 1, antecedent_set, consequent_set);
			return false;
		}
	} else {
		for (unsigned int i = 0; i < antecedent->array.length; i++) {
			Formula* atom = antecedent->array.operands[i];
			if (!get_conjunct_step(atom, instance, conjunct_steps)) {
				free(proposed_proofs); free(*axiom_step);
				if (axiom_step->reference_count == 1)
					T.free_subset_axiom(antecedent, consequent, 1, antecedent_set, consequent_set);
				return false;
			}
		}
	}

	array<pair<Proof*, Proof*>> observation_changes(1);
	Proof* old_step = get_proof<Negated>(instance, a);
#if !defined(NDEBUG)
	if (old_step->type != nd_step_type::AXIOM)
		fprintf(stderr, "propose_universal_intro WARNING: Expected an axiom.\n");
#endif

	Proof* antecedent_proof = (antecedent_length == 1 ? conjunct_steps[0] : ProofCalculus::new_conjunction_intro(make_array_view(conjunct_steps.data, antecedent_length)));
	Term* constant = Formula::new_constant(concept_id);
	Proof* new_step = ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom_step, constant), antecedent_proof);
	free(*constant); if (constant->reference_count == 0) free(constant);
	if (new_step == NULL) {
		for (unsigned int i = 0; i < conjunct_steps.length; i++) {
			free(*conjunct_steps[i]);
			if (conjunct_steps[i]->reference_count == 0)
				free(conjunct_steps[i]);
		}
		free(proposed_proofs); free(*axiom_step);
		if (axiom_step->reference_count == 1)
			T.free_subset_axiom(antecedent, consequent, 1, antecedent_set, consequent_set);
		return false;
	}

	if (!propose_transformation(T, proposed_proofs, observation_changes, old_step, new_step)) {
		new_step->reference_count++;
		free(*new_step); if (new_step->reference_count == 0) free(new_step);
		for (unsigned int i = 0; i < conjunct_steps.length; i++) {
			free(*conjunct_steps[i]);
			if (conjunct_steps[i]->reference_count == 0)
				free(conjunct_steps[i]);
		}
		free(proposed_proofs); free(*axiom_step);
		if (axiom_step->reference_count == 1)
			T.free_subset_axiom(antecedent, consequent, 1, antecedent_set, consequent_set);
		return false;
	}

	old_step->reference_count++;
	bool success = do_mh_universal_intro(T, proposed_proofs, observation_changes,
			make_universal_intro_proposal<Negated>(concept_id, old_step, axiom_step, a, new_antecedent_set_size_axiom, new_consequent_set_size_axiom),
			log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector);
	free(proposed_proofs);
	free(*old_step); if (old_step->reference_count == 0) free(old_step);
	free(*axiom_step);
	if (axiom_step->reference_count == 1)
		T.free_subset_axiom(antecedent, consequent, 1, antecedent_set, consequent_set);
	return success;
}

template<typename Formula>
void free_proofs(nd_step<Formula>** proofs, unsigned int count) {
	for (unsigned int i = 0; i < count; i++) {
		if (proofs[i] == NULL) continue;
		free(*proofs[i]);
		if (proofs[i]->reference_count == 0)
			free(proofs[i]);
	}
}

template<typename Formula, typename Canonicalizer>
inline nd_step<Formula>* make_grounded_conjunct(
		const theory<natural_deduction<Formula>, Canonicalizer>& T,
		const concept<natural_deduction<Formula>>& c,
		Formula* lifted_conjunct, unsigned int variable,
		typename Formula::Term* constant)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula> ProofCalculus;

	/* first check if the grounded conjunct already exists as an axiom in the theory */
	bool contains;
	const Formula* operand = lifted_conjunct;
	Term const* predicate; Term const* arg1; Term const* arg2;
	if (is_atomic(*operand, predicate, arg1, arg2)) {
		if (arg2 == NULL) {
			/* atom is unary */
			c.types.get(*lifted_conjunct, contains);
			if (contains) {
				/* we have to fail here as the transformation is no longer reversible */
				fprintf(stderr, "make_grounded_conjunct WARNING: The theory is in a state where this proposed transformation is not invertible.\n");
				return NULL;
			}
		} else {
			/* atom is binary */
			if (predicate->type != TermType::CONSTANT)
				return nullptr;
			relation r = { predicate->constant,
					(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
					(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
			c.relations.get(r, contains);
			if (contains) {
				/* we have to fail here as the transformation is no longer reversible */
				fprintf(stderr, "make_grounded_conjunct WARNING: The theory is in a state where this proposed transformation is not invertible.\n");
				return NULL;
			}
		}
	} else if (operand->type == FormulaType::NOT && is_atomic(*operand->unary.operand, predicate, arg1, arg2)) {
		if (arg2 == NULL) {
			/* atom is unary */
			c.negated_types.get(*lifted_conjunct->unary.operand, contains);
			if (contains) {
				/* we have to fail here as the transformation is no longer reversible */
				fprintf(stderr, "make_grounded_conjunct WARNING: The theory is in a state where this proposed transformation is not invertible.\n");
				return NULL;
			}
		} else {
			/* atom is binary */
			if (predicate->type != TermType::CONSTANT)
				return nullptr;
			relation r = { predicate->constant,
					(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
					(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
			c.negated_relations.get(r, contains);
			if (contains) {
				/* we have to fail here as the transformation is no longer reversible */
				fprintf(stderr, "make_grounded_conjunct WARNING: The theory is in a state where this proposed transformation is not invertible.\n");
				return NULL;
			}
		}
	} else {
		fprintf(stderr, "make_grounded_conjunct ERROR: Expected a literal.\n");
		return NULL;
	}

	Term* new_variable = Formula::new_variable(variable);
	nd_step<Formula>* new_axiom = ProofCalculus::new_axiom(substitute<TermType::VARIABLE>(lifted_conjunct, new_variable, constant));
	free(*new_variable); if (new_variable->reference_count == 0) free(new_variable);
	if (new_axiom == NULL) return NULL;
	free(*new_axiom->formula);
	new_axiom->reference_count++;
	return new_axiom;
}

template<typename Formula>
struct extensional_edge {
	unsigned int consequent_set;
	unsigned int antecedent_set;
	nd_step<Formula>* axiom;
	nd_step<Formula>* child;
	nd_step<Formula>* grandchild;
};

template<bool Negated, typename Formula>
struct universal_elim_proposal {
	typedef typename Formula::Term Term;

	extensional_edge<Formula> edge;
	unsigned int constant;
	nd_step<Formula>* new_axiom;
	Term& consequent_atom;
	Formula* old_antecedent_set_size_axiom;
	Formula* old_consequent_set_size_axiom;
};

template<bool Negated, typename Formula>
inline universal_elim_proposal<Negated, Formula> make_universal_elim_proposal(
		extensional_edge<Formula> edge,
		unsigned int constant,
		nd_step<Formula>* new_axiom,
		typename Formula::Term& consequent_atom,
		Formula* old_antecedent_set_size_axiom,
		Formula* old_consequent_set_size_axiom)
{
	return {edge, constant, new_axiom, consequent_atom, old_antecedent_set_size_axiom, old_consequent_set_size_axiom};
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_universal_elim(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const extensional_edge<Formula>& selected_edge,
		double& log_proposal_probability_ratio,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		const ProposalDistribution& proposal_distribution)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) return false;

	/* this proposal is not invertible if the consequent is not atomic */
	bool is_consequent_binary; bool is_consequent_symmetric;
	Formula* consequent = selected_edge.axiom->formula->quantifier.operand->binary.right;
	Term const* arg1; Term const* arg2;
	if (is_atomic(*consequent, arg1, arg2)) {
		is_consequent_binary = (arg2 != NULL);
		is_consequent_symmetric = false;
	} else if (consequent->type == FormulaType::NOT && is_atomic(*consequent->unary.operand, arg1, arg2)) {
		is_consequent_binary = (arg2 != NULL);
		is_consequent_symmetric = (arg1->type == TermType::VARIABLE && arg2->type == TermType::VARIABLE);
	} else {
		free(proposed_proofs);
		return true;
	}

	/* this proposal is not invertible if the antecedent cannot be selected in the inverse step */
#if !defined(NDEBUG)
	if (selected_edge.grandchild->type != nd_step_type::IMPLICATION_ELIMINATION)
		fprintf(stderr, "propose_universal_elim WARNING: Expected an implication elimination step.\n");
#endif
	if (selected_edge.grandchild->operands[1]->type == nd_step_type::AXIOM && selected_edge.grandchild->operands[1]->reference_count == 2) {
		free(proposed_proofs);
		return true;
	}

	unsigned int antecedent_length;
	Formula* antecedent = selected_edge.axiom->formula->quantifier.operand->binary.left;
	if (is_atomic(*antecedent)) {
		antecedent_length = 1;
	} else if (antecedent->type == FormulaType::NOT && is_atomic(*antecedent->unary.operand)) {
		antecedent_length = 1;
	} else if (antecedent->type == FormulaType::AND) {
		antecedent_length = antecedent->array.length;
	} else {
		free(proposed_proofs); return true;
	}

	array<pair<Proof*, Proof*>> observation_changes(8);
	if (selected_edge.child->type != nd_step_type::UNIVERSAL_ELIMINATION
	 || selected_edge.child->operands[1]->term->type != TermType::CONSTANT) {
		free(proposed_proofs); return true;
	}

	unsigned int constant = selected_edge.child->operands[1]->term->constant;
	const concept<ProofCalculus>& c = T.ground_concepts[constant - T.new_constant_offset];
#if !defined(NDEBUG)
	if (c.types.keys == NULL)
		fprintf(stderr, "propose_universal_elim WARNING: Theory does not contain a concept for ID %u.\n", constant);
#endif

	/* check if `grandchild` is observed */
	Proof* new_axiom = make_grounded_conjunct(T, c, consequent,
				selected_edge.axiom->formula->quantifier.variable,
				selected_edge.child->operands[1]->term);
	if (new_axiom == NULL) {
		free(proposed_proofs); return false;
	}

	/* propose `new_step` to substitute `grandchild` */
	if (!propose_transformation(T, proposed_proofs, observation_changes, selected_edge.grandchild, new_axiom)) {
		free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
		free(proposed_proofs); return false;
	}

	/* compute the probability of inverting this proposal */
	array<pair<unsigned int, unsigned int>> satisfying_extensional_edges(16);
	bool success = false;
	if (consequent->type == FormulaType::NOT) {
		success = get_satisfying_extensional_edges<true>(T, *consequent->unary.operand, constant, satisfying_extensional_edges);
	} else {
		success = get_satisfying_extensional_edges<false>(T, *consequent, constant, satisfying_extensional_edges);
	}
	if (!success) {
		free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
		free(proposed_proofs); return false;
	}
	if (satisfying_extensional_edges.length == 0) {
		/* this could happen if the antecedent is satisfied by another
		   implication elimination (another extentional edge, for example), in
		   which case, this transformation is not invertible */
		free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
		free(proposed_proofs); return true;
	}

	unfixed_set_counter unfixed_set_count_change = {0};
	Formula* old_antecedent_set_size_axiom = nullptr;
	Formula* old_consequent_set_size_axiom = nullptr;
	if (selected_edge.axiom->reference_count == 3) {
		/* if this edge is removed from the set graph, two unfixed sets may also be freed */
		bool size_axiom_is_used_in_proof = false;
		for (Proof* size_axiom : T.sets.sets[selected_edge.consequent_set].size_axioms) {
			if (size_axiom->children.length != 0) {
				size_axiom_is_used_in_proof = true;
				break;
			}
		}
		if (T.sets.is_unfixed(selected_edge.consequent_set, T.observations)
		 && T.sets.extensional_graph.vertices[selected_edge.consequent_set].parents.size
		  + T.sets.extensional_graph.vertices[selected_edge.consequent_set].children.size == 1
		 && !size_axiom_is_used_in_proof)
		{
			double log_prob = 0.0;
			set_size_proposal_log_probability(selected_edge.consequent_set, T.sets, log_prob);
			unfixed_set_count_change.log_probability -= log_prob;
			unfixed_set_count_change.change--;
			old_consequent_set_size_axiom = T.sets.sets[selected_edge.consequent_set].size_axioms[0]->formula;
		}

		size_axiom_is_used_in_proof = false;
		for (Proof* size_axiom : T.sets.sets[selected_edge.antecedent_set].size_axioms) {
			if (size_axiom->children.length != 0) {
				size_axiom_is_used_in_proof = true;
				break;
			}
		}
		if (T.sets.is_unfixed(selected_edge.antecedent_set, T.observations)
		 && T.sets.extensional_graph.vertices[selected_edge.antecedent_set].parents.size
		  + T.sets.extensional_graph.vertices[selected_edge.antecedent_set].children.size == 1
		 && !size_axiom_is_used_in_proof)
		{
			double log_prob = 0.0;
			set_size_proposal_log_probability(selected_edge.antecedent_set, T.sets, log_prob);
			unfixed_set_count_change.log_probability -= log_prob;
			unfixed_set_count_change.change--;
			old_antecedent_set_size_axiom = T.sets.sets[selected_edge.antecedent_set].size_axioms[0]->formula;
		}
	}

	unsigned int new_axiom_indices_length = c.types.size + c.negated_types.size + c.relations.size + c.negated_relations.size;
	unsigned int new_satisfying_extensional_edges_length = satisfying_extensional_edges.length;
	if (is_consequent_symmetric) {
		new_axiom_indices_length += 3;
	} else if (is_consequent_binary) {
		new_axiom_indices_length += 2;
	} else {
		new_axiom_indices_length += 1;
	}
	if (selected_edge.axiom->children.length == 1 && selected_edge.child->children.length == 1) {
		/* this is the last reference to this axiom */
		new_satisfying_extensional_edges_length--;
	}
	if (!log_cache<double>::instance().ensure_size(max(new_axiom_indices_length, new_satisfying_extensional_edges_length) + 1)) {
		free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
		free(proposed_proofs); return false;
	}

	double log_rebuild_probability = 0.0;
	new_axiom_indices_length--;
	/* this is the probability of selecting the length of the antecedent */
	log_rebuild_probability += -log_cache<double>::instance().get(new_axiom_indices_length);
	/* this is the probability of selecting the conjuncts in the antecedent */
	log_rebuild_probability += -lgamma(new_axiom_indices_length + 1)
			+ lgamma(antecedent_length + 1) + lgamma(new_axiom_indices_length - antecedent_length + 1);

	log_proposal_probability_ratio += logsumexp(
			-log_cache<double>::instance().get(new_satisfying_extensional_edges_length + 1), /* log probability of relifting the grounded consequent to this extensional edge */
			-log_cache<double>::instance().get(new_satisfying_extensional_edges_length + 1) + log_rebuild_probability) /* log probability of rebuilding this extensional edge */
			+ unfixed_set_count_change.log_probability; /* log probability of selecting the set size when creating sets */

	/* go ahead and compute the probability ratios and perform the accept/reject step */
	selected_edge.axiom->reference_count++;
	if (consequent->type == FormulaType::NOT) {
		success = do_mh_universal_elim(T, proposed_proofs, observation_changes,
				make_universal_elim_proposal<true>(selected_edge, constant, new_axiom, *consequent->unary.operand, old_antecedent_set_size_axiom, old_consequent_set_size_axiom),
				log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector);
	} else {
		success = do_mh_universal_elim(T, proposed_proofs, observation_changes,
				make_universal_elim_proposal<false>(selected_edge, constant, new_axiom, *consequent, old_antecedent_set_size_axiom, old_consequent_set_size_axiom),
				log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector);
	}
	free(*new_axiom); if (new_axiom->reference_count == 0) free(new_axiom);
	free(proposed_proofs); free(*selected_edge.axiom);
	if (selected_edge.axiom->reference_count == 1)
		T.sets.free_subset_axiom(selected_edge.axiom);
	return success;
}

template<typename Formula, typename Canonicalizer,
	typename SizePrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_change_set_size(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		unsigned int selected_set,
		double& log_proposal_probability_ratio,
		SizePrior& set_size_prior,
		TheorySampleCollector& sample_collector,
		const ProposalDistribution& proposal_distribution)
{
	/* compute the upper and lower bounds of the size of this set */
	unsigned int lower_bound = 0, upper_bound = 0;
	if (!T.sets.get_size_lower_bound(selected_set, lower_bound)
	 || !T.sets.get_size_upper_bound(selected_set, upper_bound)) return false;

	unsigned int old_size = T.sets.sets[selected_set].set_size;
	unsigned int new_size = sample(set_size_prior, lower_bound, upper_bound);
	log_proposal_probability_ratio += log_probability(proposal_distribution, 0, 0, 0);
	double proof_prior_diff = log_probability(new_size, set_size_prior) - log_probability(old_size, set_size_prior);
	log_proposal_probability_ratio += proof_prior_diff;

#if !defined(NDEBUG)
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "propose_change_set_size WARNING: The computed log probability ratio is NaN.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		T.sets.template set_size_axiom<false>(selected_set, new_size);

		bool is_consistent = true;
		if (T.sets.sets[selected_set].element_count() == new_size) {
			/* if the set becomes full, check for consistency by calling `check_set_membership_after_addition` */
			hol_term* set_formula = T.sets.sets[selected_set].set_formula();
			array<hol_term*> conjuncts(8);
			if (set_formula->type == hol_term_type::AND) {
				if (!conjuncts.append(set_formula->array.operands, set_formula->array.length))
					return false;
			} else {
				conjuncts[conjuncts.length++] = set_formula;
			}

			for (hol_term* conjunct : conjuncts) {
				if (conjunct->type == hol_term_type::UNARY_APPLICATION
				 && conjunct->binary.right->type == hol_term_type::VARIABLE
				 && conjunct->binary.right->variable == 1
				 && !T.template check_set_membership_after_addition<false>(conjunct))
				{
					is_consistent = false;
					break;
				}
			}
		} else if (T.sets.sets[selected_set].element_count() == old_size) {
			/* if the set becomes not full, check for consistency by calling `check_set_membership_after_subtraction` */
			hol_term* set_formula = T.sets.sets[selected_set].set_formula();
			array<hol_term*> conjuncts(8);
			if (set_formula->type == hol_term_type::AND) {
				if (!conjuncts.append(set_formula->array.operands, set_formula->array.length))
					return false;
			} else {
				conjuncts[conjuncts.length++] = set_formula;
			}

			for (hol_term* conjunct : conjuncts) {
				if (conjunct->type == hol_term_type::UNARY_APPLICATION
				 && conjunct->binary.right->type == hol_term_type::VARIABLE
				 && conjunct->binary.right->variable == 1
				 && !T.template check_set_membership_after_subtraction(conjunct))
				{
					is_consistent = false;
					break;
				}
			}
		}

		if (!is_consistent) {
			T.sets.template set_size_axiom<false>(selected_set, old_size);
		} else {
			/* TODO: we could keep track of the extra axioms within `theory` */
			array<hol_term*> extra_axioms(16);
			T.get_extra_axioms(extra_axioms);
			bool is_axiom = (T.sets.sets[selected_set].size_axioms[0]->children.length != 0
						  || extra_axioms.contains(T.sets.sets[selected_set].size_axioms[0]->formula));
			if (!sample_collector.accept(T.observations, extra_axioms, is_axiom ? proof_prior_diff : 0.0))
				return false;
		}
	}
	return true;
}

template<unsigned int NewConceptProbabilityAlpha>
struct proof_sampler {
	/* this tracking of inconsistent constants may too liberally exclude
	   constants if the expression contains a set size term, since we don't
	   keep track of inconsistent set sizes here */
	/*hash_map<Formula, hash_set<unsigned int>> inconsistencies;*/
	array<unsigned int> removed_set_sizes;
	/*bool all_descendants_inconsistent;*/
	double log_probability;
	double set_size_log_probability;
	bool undo;

	proof_sampler() : /*inconsistencies(64),*/ removed_set_sizes(4), log_probability(0.0), set_size_log_probability(0.0), undo(false) { }
	~proof_sampler() {
		/*for (auto entry : inconsistencies) {
			free(entry.key);
			free(entry.value);
		}*/
	}
};

template<typename Proof, unsigned int NewConceptProbabilityAlpha>
inline void visit_node(const Proof& proof, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { }

template<bool Negated, typename Term, unsigned int NewConceptProbabilityAlpha> constexpr bool visit_unary_atom(const Term* term, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<bool Negated, unsigned int NewConceptProbabilityAlpha> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<typename Proof, unsigned int NewConceptProbabilityAlpha> constexpr bool visit_subset_axiom(const Proof& proof, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<unsigned int NewConceptProbabilityAlpha> constexpr bool visit_existential_intro(const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<unsigned int NewConceptProbabilityAlpha> constexpr bool visit_negated_universal_intro(const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<unsigned int NewConceptProbabilityAlpha> constexpr bool visit_negated_conjunction(const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }
template<unsigned int NewConceptProbabilityAlpha> constexpr bool visit_disjunction_intro(const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }

template<unsigned int NewConceptProbabilityAlpha>
inline void on_subtract_changes(proof_sampler<NewConceptProbabilityAlpha>& visitor) {
	visitor.undo = true;
}

template<typename Formula, unsigned int NewConceptProbabilityAlpha>
constexpr bool on_undo_filter_operands(const Formula* formula, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }

template<typename Theory, typename Formula, unsigned int NewConceptProbabilityAlpha>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable, const proof_sampler<NewConceptProbabilityAlpha>& visitor) { return true; }

template<typename Formula, unsigned int NewConceptProbabilityAlpha>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, proof_sampler<NewConceptProbabilityAlpha>& sampler)
{
	if (!filter_operands(formula, indices))
		return false;

	if (!log_cache<double>::instance().ensure_size(indices.length + 1)) return false;
	sampler.log_probability += -log_cache<double>::instance().get(indices.length);

	indices.length = 1;
	return true;
}

template<typename ProofCalculus, typename Canonicalizer, unsigned int NewConceptProbabilityAlpha>
inline bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		proof_sampler<NewConceptProbabilityAlpha>& sampler)
{
	if (!filter_constants_helper(T, formula, variable, constants))
		return false;

	if (NewConceptProbabilityAlpha == 1) {
		if (!log_cache<double>::instance().ensure_size(constants.length + 1)) return false;
		sampler.log_probability += -log_cache<double>::instance().get(constants.length);
	} else {
		unsigned int index = 0;
		for (; index < constants.length; index++)
			if (constants[index].type == instance_type::ANY) break;
		if (index < constants.length) {
			double pi = (NewConceptProbabilityAlpha - 1.0) / (constants.length - 1 + NewConceptProbabilityAlpha);
			if (sample_uniform<double>() < pi)
				swap(constants[0], constants[index]);
			if (constants[0].type == instance_type::ANY) {
				sampler.log_probability += log((double) NewConceptProbabilityAlpha / (constants.length - 1 + NewConceptProbabilityAlpha));
			} else {
				sampler.log_probability += -log(constants.length - 1 + NewConceptProbabilityAlpha);
			}
		} else {
			if (!log_cache<double>::instance().ensure_size(constants.length + 1)) return false;
			sampler.log_probability += -log_cache<double>::instance().get(constants.length);
		}
	}

	/*if (!sampler.inconsistencies.check_size()) return false;
	bool contains; unsigned int bucket;
	hash_set<unsigned int>& inconsistent_constants = sampler.inconsistencies.get(*formula, contains, bucket);
	if (!contains) {
		if (!hash_set_init(inconsistent_constants, 8)) {
			return false;
		} else if (!init(sampler.inconsistencies.table.keys[bucket], *formula)) {
			free(inconsistent_constants);
			return false;
		}
		sampler.inconsistencies.table.size++;
	}
	unsigned int original_constant_count = constants.length;
	constants.length = 0;
	for (unsigned int i = 0; i < original_constant_count; i++) {
		if (!inconsistent_constants.contains(constants[i])) {
			constants[0] = constants[i];
			constants.length = 1;
			break;
		}
	}
	sampler.all_descendants_inconsistent = true;*/
	constants.length = 1;
	return true;
}

template<typename Formula, unsigned int NewConceptProbabilityOffset>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, proof_sampler<NewConceptProbabilityOffset>& sampler) { return true; }

template<typename Formula, unsigned int NewConceptProbabilityOffset>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, proof_sampler<NewConceptProbabilityOffset>& sampler) { return true; }

template<typename Formula, unsigned int NewConceptProbabilityOffset>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count, proof_sampler<NewConceptProbabilityOffset>& sampler) {
	/*bool contains; unsigned int bucket;
	hash_set<unsigned int>& inconsistent_constants = sampler.inconsistencies.get(*formula, contains, bucket);
	sampler.all_descendants_inconsistent = (inconsistent_constants.size == original_constant_count);
	if (sampler.all_descendants_inconsistent) {
		free(sampler.inconsistencies.table.keys[bucket]);
		core::free(inconsistent_constants);
		sampler.inconsistencies.remove_at(bucket);
	}*/
}

struct inverse_proof_sampler {
	array<unsigned int> removed_set_sizes;
	double log_probability;
	double set_size_log_probability;

	inverse_proof_sampler() : removed_set_sizes(4), log_probability(0.0), set_size_log_probability(0.0) { }
};

template<typename Proof>
inline void visit_node(const Proof& proof, const inverse_proof_sampler& visitor) { }

template<bool Negated, typename Term> constexpr bool visit_unary_atom(const Term* term, const inverse_proof_sampler& visitor) { return true; }
template<bool Negated> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const inverse_proof_sampler& visitor) { return true; }
template<typename Proof> constexpr bool visit_subset_axiom(const Proof& proof, const inverse_proof_sampler& visitor) { return true; }
constexpr bool visit_existential_intro(const inverse_proof_sampler& visitor) { return true; }
constexpr bool visit_negated_universal_intro(const inverse_proof_sampler& visitor) { return true; }
constexpr bool visit_negated_conjunction(const inverse_proof_sampler& visitor) { return true; }
constexpr bool visit_disjunction_intro(const inverse_proof_sampler& visitor) { return true; }

inline void on_subtract_changes(inverse_proof_sampler& visitor) { }

template<typename Formula>
constexpr bool filter_operands(const Formula* formula, array<unsigned int>& indices, inverse_proof_sampler& sampler) { return true; }

template<typename ProofCalculus, typename Canonicalizer>
constexpr bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		inverse_proof_sampler& sampler)
{
	return true;
}

template<typename Formula>
inline bool on_undo_filter_operands(Formula* formula, inverse_proof_sampler& sampler) {
	typedef typename Formula::Type FormulaType;

	array<unsigned int> indices(8);
	if (formula->type == FormulaType::NOT && formula->unary.operand->type == FormulaType::AND) {
		if (!indices.ensure_capacity(formula->unary.operand->array.length))
			return false;
		for (unsigned int i = 0; i < formula->unary.operand->array.length; i++)
			indices[i] = i + 1;
		indices.length = formula->unary.operand->array.length;
	} else if (formula->type == FormulaType::IF_THEN) {
		indices[0] = 1; indices[1] = 2;
		indices.length = 2;
	} else if (formula->type == FormulaType::OR) {
		if (!indices.ensure_capacity(formula->array.length))
			return false;
		for (unsigned int i = 0; i < formula->array.length; i++)
			indices[i] = i + 1;
		indices.length = formula->array.length;
	}

	if (!filter_operands(formula, indices))
		return false;

	if (!log_cache<double>::instance().ensure_size(indices.length + 1)) return false;
	sampler.log_probability += -log_cache<double>::instance().get(indices.length);
	return true;
}

template<typename Theory, typename Formula>
inline bool on_undo_filter_constants(Theory& T, Formula* quantified, unsigned int variable, inverse_proof_sampler& sampler)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	array<instance> constants(T.ground_concept_capacity + 1);
	array<hol_number> numbers(64); array<string*> strings(64);
	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys != NULL) {
			constants[constants.length].type = instance_type::CONSTANT;
			constants[constants.length++].constant = T.new_constant_offset + i;

			for (const auto& entry : T.ground_concepts[i].function_values) {
				Term* constant = entry.value->formula->binary.right;
				if (constant->type == TermType::NUMBER) {
					if (!numbers.contains(constant->number) && !numbers.add(constant->number)) return false;
				} else if (constant->type == TermType::STRING) {
					bool contains = false;
					for (const string* str : strings)
						if (str == &constant->str || *str == constant->str) { contains = true; break; }
					if (!contains && !strings.add(&constant->str)) return false;
				}
			}
		}
	} for (unsigned int i = 1; i < T.sets.set_count + 1; i++) {
		if (T.sets.sets[i].size_axioms.data == nullptr) continue;
		hol_number number;
		number.integer = T.sets.sets[i].set_size;
		number.decimal = 0;
		if (!numbers.contains(number) && !numbers.add(number))
			return false;
	}
	constants[constants.length++].type = instance_type::ANY;
	if (!constants.ensure_capacity(constants.length + numbers.length + strings.length))
		return false;
	for (hol_number number : numbers) {
		constants[constants.length].type = instance_type::NUMBER;
		constants[constants.length++].number = number;
	} for (string* str : strings) {
		constants[constants.length].type = instance_type::STRING;
		constants[constants.length++].str = str;
	}

	if (!filter_constants_helper(T, quantified, variable, constants))
		return false;

	if (!log_cache<double>::instance().ensure_size(constants.length + 1)) return false;
	sampler.log_probability += -log_cache<double>::instance().get(constants.length);
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, inverse_proof_sampler& sampler) { return true; }

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, inverse_proof_sampler& sampler) { return true; }

template<typename Formula>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count, inverse_proof_sampler& sampler) { }

struct undo_remove_sets {
	array<unsigned int>& removed_set_sizes;

	undo_remove_sets(array<unsigned int>& removed_set_sizes) : removed_set_sizes(removed_set_sizes) { }
};

inline void on_subtract_changes(const undo_remove_sets& visitor) { }

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const undo_remove_sets& visitor) { return true; }

template<typename Theory, typename Formula>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable, const undo_remove_sets& visitor) { return true; }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		const undo_remove_sets& visitor)
{
	out = visitor.removed_set_sizes.last();
#if !defined(NDEBUG)
	if (visitor.removed_set_sizes.length == 0)
		fprintf(stderr, "compute_new_set_size WARNING: The recorded array of set sizes is empty.\n");
	if (out < min_set_size || out > max_set_size)
		fprintf(stderr, "compute_new_set_size WARNING: The recorded set size in `undo_remove_sets.removed_set_sizes` is outside the bounds of `[min_set_size, max_set_size]`.\n");
#endif
	visitor.removed_set_sizes.length--;
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const undo_remove_sets& visitor)
{ }

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const undo_remove_sets& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const undo_remove_sets& visitor)
{ }

template<typename Formula, typename Canonicalizer>
bool undo_proof_changes(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		typename theory<natural_deduction<Formula>, Canonicalizer>::changes& old_proof_changes,
		typename theory<natural_deduction<Formula>, Canonicalizer>::changes& new_proof_changes,
		nd_step<Formula>* old_proof, nd_step<Formula>* new_proof,
		const undo_remove_sets& old_sets, const undo_remove_sets& new_sets)
{
	typedef nd_step<Formula> Proof;

	array<pair<Proof*, Proof*>> observation_changes(8);
	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) {
		return false;
	} else if (!propose_transformation(T, proposed_proofs, observation_changes, new_proof, old_proof)) {
		free(proposed_proofs); return false;
	} else if (!transform_proofs(proposed_proofs)) {
		free(proposed_proofs); return false;
	}

	for (auto& entry : observation_changes) {
		T.observations.remove(entry.key);
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		T.observations.add(entry.value);
		entry.value->reference_count++;
	}
	free(proposed_proofs);

	set_changes<Formula> set_diff;
	if (!T.subtract_changes(new_proof_changes, set_diff, new_sets)) return false;

	/* add the changes back into T from `old_proof` */
	return T.add_changes(old_proof_changes, set_diff, old_sets);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void set_size_proposal_log_probability(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		double& log_probability_value)
{
	unsigned int min_set_size, max_set_size;
	if (!sets.set_size_bounds(set_id, min_set_size, max_set_size)) {
		fprintf(stderr, "set_size_proposal_log_probability ERROR: `set_size_bounds` failed.\n");
		return;
	}

#if !defined(NDEBUG)
	if (sets.sets[set_id].set_size < min_set_size || sets.sets[set_id].set_size > max_set_size)
		fprintf(stderr, "set_size_proposal_log_probability WARNING: The set with ID %u has size outside the bounds computed by `set_reasoning.set_size_bounds`.\n", set_id);
#endif

	if (max_set_size == UINT_MAX) {
		log_probability_value += (sets.sets[set_id].set_size - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2);
		log_probability_value += -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
}

struct inverse_set_size_log_probability {
	double value;
	array<unsigned int> removed_set_sizes;

	inverse_set_size_log_probability() : value(0), removed_set_sizes(4) { }
};

inline void on_subtract_changes(const inverse_set_size_log_probability& visitor) { }

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const inverse_set_size_log_probability& visitor) { return true; }

template<typename Theory, typename Formula>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable, const inverse_set_size_log_probability& visitor) { return true; }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		inverse_set_size_log_probability& visitor)
{
	if (max_set_size == UINT_MAX) {
		out = min_set_size + sample_geometric(0.1);
		visitor.value -= (out - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		out = min_set_size + sample_uniform(max_set_size - min_set_size + 1);
		if (!log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2))
			return false;
		visitor.value -= -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		inverse_set_size_log_probability& visitor)
{
	set_size_proposal_log_probability(set_id, sets, visitor.value);
	visitor.removed_set_sizes.add(sets.sets[set_id].set_size);
}

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const inverse_set_size_log_probability& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const inverse_set_size_log_probability& visitor)
{ }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, unsigned int NewConceptProbabilityOffset>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		proof_sampler<NewConceptProbabilityOffset>& sampler)
{
	if (sampler.undo) {
		out = sampler.removed_set_sizes.last();
#if !defined(NDEBUG)
		if (sampler.removed_set_sizes.length == 0)
			fprintf(stderr, "compute_new_set_size WARNING: The recorded array of set sizes is empty.\n");
		if (out < min_set_size || out > max_set_size)
			fprintf(stderr, "compute_new_set_size WARNING: The recorded set size in `proof_sampler.removed_set_sizes` is outside the bounds of `[min_set_size, max_set_size]`.\n");
#endif
		sampler.removed_set_sizes.length--;
		return true;
	}

	if (max_set_size == UINT_MAX) {
		out = min_set_size + sample_geometric(0.1);
		sampler.set_size_log_probability += (out - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		out = min_set_size + sample_uniform(max_set_size - min_set_size + 1);
		if (!log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2))
			return false;
		sampler.set_size_log_probability += -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, unsigned int NewConceptProbabilityAlpha>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		proof_sampler<NewConceptProbabilityAlpha>& sampler)
{
	if (sampler.undo) return;
	sampler.removed_set_sizes.add(sets.sets[set_id].set_size);
}

template<typename Proof, unsigned int NewConceptProbabilityAlpha>
constexpr bool on_new_size_axiom(Proof* new_size_axiom,
		const proof_sampler<NewConceptProbabilityAlpha>& visitor)
{
	return true;
}

template<typename Proof, unsigned int NewConceptProbabilityAlpha>
inline void on_old_size_axiom(Proof* old_size_axiom,
		const proof_sampler<NewConceptProbabilityAlpha>& visitor)
{ }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		inverse_proof_sampler& visitor)
{
	if (max_set_size == UINT_MAX) {
		out = min_set_size + sample_geometric(0.1);
		visitor.set_size_log_probability -= (out - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		out = min_set_size + sample_uniform(max_set_size - min_set_size + 1);
		if (!log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2))
			return false;
		visitor.set_size_log_probability -= -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		inverse_proof_sampler& visitor)
{
	set_size_proposal_log_probability(set_id, sets, visitor.set_size_log_probability);
	visitor.removed_set_sizes.add(sets.sets[set_id].set_size);
}

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const inverse_proof_sampler& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const inverse_proof_sampler& visitor)
{ }

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_disjunction_intro(
	theory<natural_deduction<Formula>, Canonicalizer>& T,
	pair<Formula*, nd_step<Formula>*> selected_step,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	typedef nd_step<Formula> Proof;
	typedef theory<natural_deduction<Formula>, Canonicalizer> Theory;

	set_changes<Formula> set_diff;
	inverse_proof_sampler inverse_sampler;
	typename Theory::changes old_proof_changes;
	if (!T.get_theory_changes(*selected_step.value, old_proof_changes))
		return false;
	/* check to make sure this proof wouldn't remove any subset edges, since we cannot make the inverse proposal */
	/*for (const typename Theory::change& c : old_proof_changes.list)
		if (c.type == Theory::change_type::SUBSET_AXIOM) return true;*/
	T.subtract_changes(old_proof_changes, set_diff, inverse_sampler);
	/* some removed set size axioms may actually become extra axioms */
	for (const typename Theory::change& change : old_proof_changes.list) {
		if (change.type == Theory::change_type::SET_SIZE_AXIOM) {
			bool was_removed = false;
			for (Formula* old_axiom : set_diff.old_set_axioms) {
				if (old_axiom == change.axiom->formula) {
					was_removed = true;
					break;
				}
			}
			if (!was_removed)
				set_diff.new_set(change.axiom->formula);
		}
	}

	nd_step<Formula>* new_proof;
	proof_sampler<ProposalDistribution::NewConceptProbabilityAlpha> sampler;
	set_changes<Formula> new_set_diff;
unsigned int debug = 0;
bool debug_flag = false;
if (debug_flag) {
print(*selected_step.key, stderr); print('\n', stderr);
}
	while (true) {
		/* sample a new proof, only selecting one path at every branch,
		   and avoiding paths that we've previously proved to be inconsistent,
		   also compute the log probability of the new path */
		unsigned int new_constant = 0;
		sampler.log_probability = 0.0;
		sampler.set_size_log_probability = 0.0;
		sampler.removed_set_sizes.clear();
		sampler.undo = false;
if (debug_flag) {
fprintf(stderr, "INNER DEBUG: %u\n", debug);
T.print_axioms(stderr);
}
debug++;
		new_set_diff.clear();
		new_proof = T.template make_proof<false, true, false>(selected_step.key, new_set_diff, new_constant, sampler);
		if (new_proof != NULL)
			break;
	}
if (debug_flag) {
fprintf(stderr, "INNER DEBUG: %u (loop broken)\n", debug);
T.print_axioms(stderr);
}

	log_proposal_probability_ratio -= sampler.log_probability;
	log_proposal_probability_ratio -= sampler.set_size_log_probability;
	log_proposal_probability_ratio += inverse_sampler.log_probability;
	log_proposal_probability_ratio += inverse_sampler.set_size_log_probability;

	typename Theory::changes new_proof_changes;
	if (!T.get_theory_changes(*new_proof, new_proof_changes))
		return false;

	/* check if the proposed proof is the same as the original proof */
	if (*selected_step.value == *new_proof) {
		undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, new_proof, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		return true;
	}

	/* some new set size axioms may actually have already been extra axioms */
	for (const typename Theory::change& change : new_proof_changes.list) {
		if (change.type == Theory::change_type::SET_SIZE_AXIOM) {
			bool is_new = false;
			for (Formula* new_axiom : new_set_diff.new_set_axioms) {
				if (new_axiom == change.axiom->formula) {
					is_new = true;
					break;
				}
			}
			if (!is_new)
				new_set_diff.old_set(change.axiom->formula);
		}
	}
	for (Formula* old_axiom : new_set_diff.old_set_axioms)
		set_diff.old_set(old_axiom);
	for (Formula* new_axiom : new_set_diff.new_set_axioms)
		set_diff.new_set(new_axiom);

	/* propose `new_proof` to substitute `selected_step.value` */
	array<pair<Proof*, Proof*>> observation_changes(8);
	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) {
		undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, new_proof, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		return false;
	} else if (!propose_transformation(T, proposed_proofs, observation_changes, selected_step.value, new_proof)) {
		undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, new_proof, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(proposed_proofs); return false;
	}

	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			set_diff.old_set_axioms, set_diff.new_set_axioms,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

	if (!transform_proofs(proposed_proofs)) {
		undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, new_proof, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(*selected_step.value); if (selected_step.value->reference_count == 0) free(selected_step.value);
		free(proposed_proofs); return false;
	}

	for (auto& entry : observation_changes) {
		T.observations.remove(entry.key);
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		if (!T.observations.add(entry.value)) {
			undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, new_proof, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
			free(*selected_step.value); if (selected_step.value->reference_count == 0) free(selected_step.value);
			free(proposed_proofs); return false;
		}
		entry.value->reference_count++;
	}

	/* count all unfixed sets */
	array<unsigned int> unfixed_sets(8);
	if (!T.sets.get_unfixed_sets(unfixed_sets, T.observations)) return false;

	/* count all eliminable extensional edges */
	array<extensional_edge<Formula>> eliminable_extensional_edges(8);
	get_eliminable_extensional_edges(T, eliminable_extensional_edges);

	array<pair<relation, relation>> mergeable_events(8);
	get_mergeable_events(T, mergeable_events);

	array<relation> splittable_events(8);
	get_splittable_events(T, splittable_events);

	log_proposal_probability_ratio += log_probability(proposal_distribution, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_step);

	bool success = do_mh_disjunction_intro(T, selected_step, new_proof, observation_changes,
			old_proof_changes, new_proof_changes, proof_axioms, old_axioms, new_axioms, log_proposal_probability_ratio,
			log_proposal_probability_ratio + sampler.set_size_log_probability - inverse_sampler.set_size_log_probability,
			undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes),
			proof_prior_diff, sample_collector);
	free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
	free(proposed_proofs);
	return success;
}

struct proof_initializer {
	array<instance> expected_constants;
	array<unsigned int> expected_operand_indices;
	unsigned int constant_position;
	unsigned int operand_position;
	array<unsigned int> removed_set_sizes;
	double set_size_log_probability;
	bool undo;

	proof_initializer() : expected_constants(8), expected_operand_indices(4), constant_position(0), operand_position(0), removed_set_sizes(4), set_size_log_probability(0.0), undo(false) { }

	static inline void free(proof_initializer& initializer) {
		core::free(initializer.expected_constants);
		core::free(initializer.expected_operand_indices);
		core::free(initializer.removed_set_sizes);
	}
};

inline bool init(proof_initializer& initializer) {
	if (!array_init(initializer.expected_constants, 8)) {
		return false;
	} else if (!array_init(initializer.expected_operand_indices, 8)) {
		free(initializer.expected_constants);
		return false;
	} else if (!array_init(initializer.removed_set_sizes, 4)) {
		free(initializer.expected_constants);
		free(initializer.expected_operand_indices);
		return false;
	}
	initializer.constant_position = 0;
	initializer.operand_position = 0;
	initializer.set_size_log_probability = 0.0;
	initializer.undo = false;
	return true;
}

template<typename Proof>
inline void visit_node(const Proof& proof, const proof_initializer& visitor) { }

template<bool Negated, typename Term> constexpr bool visit_unary_atom(const Term* term, const proof_initializer& visitor) { return true; }
template<bool Negated> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const proof_initializer& visitor) { return true; }
template<typename Proof> constexpr bool visit_subset_axiom(const Proof& proof, const proof_initializer& visitor) { return true; }
constexpr bool visit_existential_intro(const proof_initializer& visitor) { return true; }
constexpr bool visit_negated_universal_intro(const proof_initializer& visitor) { return true; }
constexpr bool visit_negated_conjunction(const proof_initializer& visitor) { return true; }
constexpr bool visit_disjunction_intro(const proof_initializer& visitor) { return true; }

inline void on_subtract_changes(proof_initializer& visitor) {
	visitor.undo = true;
}

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const proof_initializer& visitor) { return true; }

template<typename Theory, typename Formula>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable, const proof_initializer& visitor) { return true; }

template<typename Formula>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, proof_initializer& initializer)
{
	if (!filter_operands(formula, indices))
		return false;

	if (initializer.operand_position == initializer.expected_operand_indices.length) {
		fprintf(stderr, "filter_operands ERROR: `proof_initializer` has no further `expected_operand_indices`.\n");
		exit(EXIT_FAILURE);
	} else if (!indices.contains(initializer.expected_operand_indices[initializer.operand_position])) {
		indices.length = 0;
		return false;
	}
	indices[0] = initializer.expected_operand_indices[initializer.operand_position++];
	indices.length = 1;
	return true;
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		proof_initializer& initializer)
{
	if (!filter_constants_helper(T, formula, variable, constants))
		return false;

	if (initializer.constant_position == initializer.expected_constants.length) {
		fprintf(stderr, "filter_constants ERROR: `proof_initializer` has no further `expected_constants`.\n");
		exit(EXIT_FAILURE);
	}
	instance expected_constant = initializer.expected_constants[initializer.constant_position];
	if (expected_constant.type == instance_type::CONSTANT && T.ground_concepts[expected_constant.constant - T.new_constant_offset].types.keys == nullptr)
		expected_constant = instance_any();
	if (!constants.contains(expected_constant)) {
		constants.length = 0;
		return false;
	}
	constants[0] = expected_constant;
	constants.length = 1;
	initializer.constant_position++;
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, proof_initializer& initializer) { return true; }

template<typename Formula>
inline bool inconsistent_constant(const Formula* formula, const instance& constant, proof_initializer& initializer) {
	initializer.constant_position--;
	return true;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count, proof_initializer& initializer) { }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(
		unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		proof_initializer& initializer)
{
	if (initializer.undo) {
		out = initializer.removed_set_sizes.last();
#if !defined(NDEBUG)
		if (initializer.removed_set_sizes.length == 0)
			fprintf(stderr, "compute_new_set_size WARNING: The recorded array of set sizes is empty.\n");
		if (out < min_set_size || out > max_set_size)
			fprintf(stderr, "compute_new_set_size WARNING: The recorded set size in `proof_initializer.removed_set_sizes` is outside the bounds of `[min_set_size, max_set_size]`.\n");
#endif
		initializer.removed_set_sizes.length--;
		return true;
	}

	if (max_set_size == UINT_MAX) {
		out = min_set_size + sample_geometric(0.1);
		initializer.set_size_log_probability += (out - min_set_size) * SET_SIZE_PROPOSAL_LOG_ONE_MINUS_P + SET_SIZE_PROPOSAL_LOG_P;
	} else {
		out = min_set_size + sample_uniform(max_set_size - min_set_size + 1);
		if (!log_cache<double>::instance().ensure_size(max_set_size - min_set_size + 2))
			return false;
		initializer.set_size_log_probability += -log_cache<double>::instance().get(max_set_size - min_set_size + 1);
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		proof_initializer& initializer)
{
	if (initializer.undo) return;
	initializer.removed_set_sizes.add(sets.sets[set_id].set_size);
}

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const proof_initializer& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const proof_initializer& visitor)
{ }

template<typename Formula>
inline bool get_proof_initializer(nd_step<Formula>* proof, proof_initializer& initializer)
{
	typedef typename Formula::TermType TermType;

	switch (proof->type) {
	case nd_step_type::IMPLICATION_INTRODUCTION:
		if (proof->operands[0]->type == nd_step_type::FALSITY_ELIMINATION
		 && proof->operands[0]->operands[0]->type == nd_step_type::NEGATION_ELIMINATION
		 && proof->operands[0]->operands[0]->operands[1] == proof->operands[1])
		{
			if (!initializer.expected_operand_indices.add(1)) return false;
		} else {
			if (!initializer.expected_operand_indices.add(0)) return false;
		}
		break;
	case nd_step_type::PROOF_BY_CONTRADICTION:
		if (proof->operands[0]->type == nd_step_type::NEGATION_ELIMINATION
		 && proof->operands[0]->operands[0]->type == nd_step_type::CONJUNCTION_ELIMINATION
		 && proof->operands[0]->operands[0]->operands[0] == proof->operands[1])
		{
			unsigned int index = proof->operands[0]->operands[0]->operands[1]->parameters[0] + 1;
			if (!initializer.expected_operand_indices.add(index)) return false;
		}
		break;

	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		if (proof->operands[2]->term->type == TermType::CONSTANT) {
			if (!initializer.expected_constants.add(instance_constant(proof->operands[2]->term->constant)))
				return false;
		} else if (proof->operands[2]->term->type == TermType::NUMBER) {
			if (!initializer.expected_constants.add(instance_number(proof->operands[2]->term->number)))
				return false;
		} else if (proof->operands[2]->term->type == TermType::STRING) {
			if (!initializer.expected_constants.add(instance_string(&proof->operands[2]->term->str)))
				return false;
		}
		break;

	case nd_step_type::DISJUNCTION_INTRODUCTION:
		if (!initializer.expected_operand_indices.add(proof->operands[2]->parameter + 1))
			return false;
		break;

	case nd_step_type::AXIOM:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_ELIMINATION:
	case nd_step_type::IMPLICATION_ELIMINATION:
	case nd_step_type::BICONDITIONAL_INTRODUCTION:
	case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT:
	case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT:
	case nd_step_type::NEGATION_ELIMINATION:
	case nd_step_type::FALSITY_ELIMINATION:
	case nd_step_type::UNIVERSAL_INTRODUCTION:
	case nd_step_type::UNIVERSAL_ELIMINATION:
	case nd_step_type::EXISTENTIAL_ELIMINATION:
	case nd_step_type::EQUALITY_ELIMINATION:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
	case nd_step_type::COUNT:
		break;
	}

	unsigned int operand_count;
	nd_step<Formula>* const* operands;
	proof->get_subproofs(operands, operand_count);
	for (unsigned int i = 0; i < operand_count; i++) {
		if (operands[i] == nullptr) continue;
		if (!get_proof_initializer(operands[i], initializer))
			return false;
	}
	return true;
}

template<typename Formula, typename Canonicalizer,
	typename DstEventType,
	typename MapConstantsFunction,
	typename ComputeLogProbabilityFunction,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
inline bool do_split_merge(
	theory<natural_deduction<Formula>, Canonicalizer>& T,
	const relation& src_event,
	const DstEventType& dst_event,
	MapConstantsFunction map_constants,
	ComputeLogProbabilityFunction compute_log_probability,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef theory<natural_deduction<Formula>, Canonicalizer> Theory;
	typedef typename Theory::changes TheoryChanges;
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	array<pair<Formula*, Proof*>> old_proofs(8);
	for (unsigned int i = 0; i < T.existential_intro_nodes.length; i++) {
		pair<Formula*, Proof*>& entry = T.existential_intro_nodes[i];
		Term* term = entry.value->operands[2]->term;
		if (term->type != TermType::CONSTANT) continue;
		if (term->constant != src_event.predicate && term->constant != src_event.arg1 && term->constant != src_event.arg2)
			continue;

		bool is_root = true;
		array<Proof*> stack(max(8, entry.value->children.length));
		for (Proof* child : entry.value->children)
			stack[stack.length++] = child;
		while (stack.length != 0) {
			Proof* current = stack.pop();
			if (current->type == nd_step_type::EXISTENTIAL_INTRODUCTION) {
				Term* term = current->operands[2]->term;
				if (term->type == TermType::CONSTANT && (term->constant == src_event.predicate || term->constant == src_event.arg1 || term->constant == src_event.arg2)) {
					is_root = false;
					break;
				}
			}

			for (Proof* child : current->children)
				if (!stack.add(child)) return false;
		}

		if (is_root && !old_proofs.add(entry))
			return false;
	}

	TheoryChanges* old_proof_changes = (TheoryChanges*) malloc(sizeof(TheoryChanges) * old_proofs.length);
	if (old_proof_changes == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		return false;
	}

	TheoryChanges* new_proof_changes = (TheoryChanges*) malloc(sizeof(TheoryChanges) * old_proofs.length);
	if (new_proof_changes == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		free(old_proof_changes); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!Theory::init(old_proof_changes[i])) {
			for (unsigned int j = 0; j < i; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); return false;
		}
	} for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!Theory::init(new_proof_changes[i])) {
			for (unsigned int j = 0; j < i; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); return false;
		}
	}

	proof_initializer* initializers = (proof_initializer*) malloc(sizeof(proof_initializer) * old_proofs.length);
	if (initializers == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		free(old_proof_changes); free(new_proof_changes); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!init(initializers[i])) {
			for (unsigned int j = 0; j < i; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); free(initializers); return false;
		}
	}

	Proof** new_proofs = (Proof**) malloc(sizeof(Proof*) * old_proofs.length);
	if (new_proofs == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		free(old_proof_changes); free(new_proof_changes); free(initializers); return false;
	}

	set_changes<Formula> set_diff;
	inverse_set_size_log_probability set_size_log_probability;
	array_map<const Proof*, unsigned int> reference_counts(32);
	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!get_proof_initializer(old_proofs[i].value, initializers[i])
		 || !reference_counts.ensure_capacity(reference_counts.size + 1))
		{
			set_changes<Formula> dummy;
			for (unsigned int j = i; j > 0; j--)
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes);
			free(initializers); free(new_proofs); return false;
		}

		for (instance& inst : initializers[i].expected_constants)
			map_constants(inst);

		array<const Proof*> discharged_axioms(16);
		unsigned int index = reference_counts.index_of(old_proofs[i].value);
		if (index == reference_counts.size) {
			reference_counts.keys[index] = old_proofs[i].value;
			reference_counts.values[index] = old_proofs[i].value->reference_count;
			reference_counts.size++;
		}
		reference_counts.values[index]--;
		if (!T.get_theory_changes(*old_proofs[i].value, discharged_axioms, reference_counts, old_proof_changes[i]))
			return false;
		/* check to make sure this proof wouldn't remove any subset edges, since we cannot make the inverse proposal */
		/*for (const typename Theory::change& c : old_proof_changes[i].list) {
			if (c.type == Theory::change_type::SUBSET_AXIOM) {
				set_changes<Formula> dummy;
				for (unsigned int j = i; j > 0; j--)
					T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
				for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
				for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
				for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
				free(old_proof_changes); free(new_proof_changes);
				free(initializers); free(new_proofs); return false;
			}
		}*/
		T.subtract_changes(old_proof_changes[i], set_diff, set_size_log_probability);
		/* some removed set size axioms may actually become extra axioms */
		for (const typename Theory::change& change : old_proof_changes[i].list) {
			if (change.type == Theory::change_type::SET_SIZE_AXIOM) {
				bool was_removed = false;
				for (Formula* old_axiom : set_diff.old_set_axioms) {
					if (old_axiom == change.axiom->formula) {
						was_removed = true;
						break;
					}
				}
				if (!was_removed)
					set_diff.new_set(change.axiom->formula);
			}
		}
	}
	log_proposal_probability_ratio += set_size_log_probability.value;

	array<pair<Proof*, Proof*>> observation_changes(8);
	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) {
		set_changes<Formula> dummy;
		for (unsigned int j = old_proofs.length; j > 0; j--)
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		free(old_proof_changes); free(new_proof_changes);
		free(initializers); free(new_proofs); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		unsigned int new_constant = 0;
		set_changes<Formula> new_set_diff;
		new_proofs[i] = T.template make_proof<false, true, false>(old_proofs[i].key, new_set_diff, new_constant, initializers[i]);
		if (new_proofs[i] == nullptr) {
			set_changes<Formula> dummy;
			for (unsigned int j = i; j > 0; j--) {
				T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
				free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
			} for (unsigned int j = old_proofs.length; j > 0; j--)
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
			free(initializers); free(new_proofs); return true;
		}

		log_proposal_probability_ratio -= initializers[i].set_size_log_probability;

		if (!T.get_theory_changes(*new_proofs[i], new_proof_changes[i])) {
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
			free(initializers); free(new_proofs); return true;
		}
		/* some new set size axioms may actually have already been extra axioms */
		for (const typename Theory::change& change : new_proof_changes[i].list) {
			if (change.type == Theory::change_type::SET_SIZE_AXIOM) {
				bool is_new = false;
				for (Formula* new_axiom : new_set_diff.new_set_axioms) {
					if (new_axiom == change.axiom->formula) {
						is_new = true;
						break;
					}
				}
				if (!is_new)
					new_set_diff.old_set(change.axiom->formula);
			}
		}

		/* propose `new_proofs[i]` to substitute `old_proofs[i].value` */
		if (!propose_transformation(T, proposed_proofs, observation_changes, old_proofs[i].value, new_proofs[i])) {
			set_changes<Formula> dummy;
			for (unsigned int j = i + 1; j > 0; j--) {
				T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
				free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
			} for (unsigned int j = old_proofs.length; j > 0; j--)
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
			free(initializers); free(new_proofs); return true;
		}

		for (Formula* old_axiom : new_set_diff.old_set_axioms)
			set_diff.old_set(old_axiom);
		for (Formula* new_axiom : new_set_diff.new_set_axioms)
			set_diff.new_set(new_axiom);
	}

	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			set_diff.old_set_axioms, set_diff.new_set_axioms,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

	if (!transform_proofs(proposed_proofs)) {
		set_changes<Formula> dummy;
		for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
			free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
		} for (unsigned int j = old_proofs.length; j > 0; j--)
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
		free(initializers); free(new_proofs); return false;
	}

	auto undo_proof_changes = [&]() {
		array<pair<Proof*, Proof*>> observation_changes(8);
		proof_transformations<Formula>& inverse_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
		if (!init(inverse_proofs))
			return false;
		for (unsigned int i = 0; i < old_proofs.length; i++) {
			if (!propose_transformation(T, inverse_proofs, observation_changes, new_proofs[i], old_proofs[i].value)) {
				free(inverse_proofs);
				return false;
			}
		}
		if (!transform_proofs(inverse_proofs)) {
			free(inverse_proofs);
			return false;
		}

		for (auto& entry : observation_changes) {
			T.observations.remove(entry.key);
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
			T.observations.add(entry.value);
			entry.value->reference_count++;
		}
		free(inverse_proofs);

		set_changes<Formula> dummy;
		for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
			free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
		} for (unsigned int j = old_proofs.length; j > 0; j--)
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
		free(initializers); free(new_proofs);
		return true;
	};

	for (auto& entry : observation_changes) {
		T.observations.remove(entry.key);
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		if (!T.observations.add(entry.value)) {
			undo_proof_changes();
			return false;
		}
		entry.value->reference_count++;
	}

	/* count all unfixed sets */
	array<unsigned int> unfixed_sets(8);
	if (!T.sets.get_unfixed_sets(unfixed_sets, T.observations)) return false;

	/* count all eliminable extensional edges */
	array<extensional_edge<Formula>> eliminable_extensional_edges(8);
	get_eliminable_extensional_edges(T, eliminable_extensional_edges);

	array<pair<relation, relation>> mergeable_events(8);
	get_mergeable_events(T, mergeable_events);

	array<relation> splittable_events(8);
	get_splittable_events(T, splittable_events);

	log_proposal_probability_ratio += log_probability(proposal_distribution, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, dst_event);
	compute_log_probability(log_proposal_probability_ratio);

#if !defined(NDEBUG)
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "do_split_merge WARNING: The computed log probability ratio is NaN.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we accepted the new proofs */
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	} else {
		return undo_proof_changes();
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		free(*new_proofs[i]); if (new_proofs[i]->reference_count == 0) free(new_proofs[i]);
	}
	for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
	for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
	for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
	free(old_proof_changes); free(new_proof_changes); free(proposed_proofs);
	free(initializers); free(new_proofs);
	return true;
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_merge_events(
	theory<natural_deduction<Formula>, Canonicalizer>& T,
	const pair<relation, relation>& events,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	auto map_constants = [events](instance& inst) {
		if (inst.type != instance_type::CONSTANT)
			return;
		if (inst.constant == events.value.predicate)
			inst.constant = events.key.predicate;
		else if (inst.constant == events.value.arg1)
			inst.constant = events.key.arg1;
		else if (inst.constant == events.value.arg2)
			inst.constant = events.key.arg2;
	};
	auto compute_log_probability = [events,&T](double& log_proposal_probability_ratio) {
		/* count the number of Bernoulli samples in the inverse proposal */
		unsigned int num_bernoulli_samples = 0;
		if (events.key.arg1 != 0)
			num_bernoulli_samples++;
		if (events.key.arg2 != 0)
			num_bernoulli_samples++;
		for (const pair<Formula*, nd_step<Formula>*>& entry : T.existential_intro_nodes) {
			Term* term = entry.value->operands[2]->term;
			if (term->type != TermType::CONSTANT) continue;
			if (term->constant == events.key.predicate)
				num_bernoulli_samples++;
			else if (term->constant == events.key.arg1)
				num_bernoulli_samples++;
			else if (term->constant == events.key.arg2)
				num_bernoulli_samples++;
		}
		log_proposal_probability_ratio += -log_cache<double>::instance().get(2) * num_bernoulli_samples;
	};
	return do_split_merge(T, events.value, events.key, map_constants, compute_log_probability, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool propose_split_event(
	theory<natural_deduction<Formula>, Canonicalizer>& T,
	const relation& event,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	unsigned int num_bernoulli_samples = 0;
	unsigned int new_predicate = T.get_free_concept_id();
	unsigned int new_arg1 = 0;
	unsigned int new_arg2 = 0;
	if (event.arg1 != 0) {
		if (event.arg2 == 0 || sample_uniform(1) == 1)
			new_arg1 = T.get_free_concept_id(new_predicate - T.new_constant_offset + 1);
		if (event.arg2 != 0)
			num_bernoulli_samples++;
	} if (event.arg2 != 0) {
		if (event.arg1 == 0 || sample_uniform(1) == 1)
			new_arg2 = T.get_free_concept_id(max(new_predicate, new_arg1) - T.new_constant_offset + 1);
		if (event.arg1 != 0)
			num_bernoulli_samples++;
	}

	unsigned int num_predicates = 0;
	unsigned int num_arg1 = 0;
	unsigned int num_arg2 = 0;
	unsigned int num_new_predicates = 0;
	unsigned int num_new_arg1 = 0;
	unsigned int num_new_arg2 = 0;
	auto map_constants = [event,new_predicate,new_arg1,new_arg2,&num_bernoulli_samples,&num_predicates,&num_arg1,&num_arg2,&num_new_predicates,&num_new_arg1,&num_new_arg2](instance& inst) {
		if (inst.type != instance_type::CONSTANT)
			return;
		if (inst.constant == event.predicate) {
			if (sample_uniform(1) == 1) {
				inst.constant = new_predicate;
				num_new_predicates++;
			}
			num_bernoulli_samples++;
			num_predicates++;
		} else if (new_arg1 != 0 && inst.constant == event.arg1) {
			if (sample_uniform(1) == 1) {
				inst.constant = new_arg1;
				num_new_arg1++;
			}
			num_bernoulli_samples++;
			num_arg1++;
		} else if (new_arg2 != 0 && inst.constant == event.arg2) {
			if (sample_uniform(1) == 1) {
				inst.constant = new_arg2;
				num_new_arg2++;
			}
			num_bernoulli_samples++;
			num_arg2++;
		}
	};
	auto compute_log_probability = [&num_bernoulli_samples,num_predicates,num_arg1,num_arg2,num_new_predicates,num_new_arg1,num_new_arg2](double& log_proposal_probability_ratio) {
		if (num_predicates == num_new_predicates || (num_arg1 != 0 && num_arg1 == num_new_arg1) || (num_arg2 != 0 && num_arg2 == num_new_arg2))
			log_proposal_probability_ratio = -std::numeric_limits<double>::infinity();
		else log_proposal_probability_ratio -= -log_cache<double>::instance().get(2) * num_bernoulli_samples;
	};
	pair<relation, relation> dst_event = {event, {new_predicate, new_arg1, new_arg2}};
	return do_split_merge(T, event, dst_event, map_constants, compute_log_probability, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

template<typename Formula>
bool transform_proofs(const proof_transformations<Formula>& proposed_proofs)
{
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	hash_map<Proof*, Proof*> transformations(32);
	for (auto entry : proposed_proofs.transformed_proofs) {
		if (!transformations.check_size(transformations.table.size + entry.value.map.size))
			return false;
		for (auto transformation : entry.value.map)
			transformations.put(transformation.key, transformation.value);
	}

	for (auto entry : transformations) {
		if (!entry.value->children.ensure_capacity(entry.key->children.length))
			return false;

		/* this is to avoid freeing `entry.key` while we are inside the next loop */
		entry.key->reference_count++;

		/* change the operand of any child to the new proof step */
		for (unsigned int i = entry.key->children.length; i > 0; i--) {
			Proof* child = entry.key->children[i - 1];

			bool error = false;
			switch (child->type) {
			case nd_step_type::CONJUNCTION_INTRODUCTION:
				for (unsigned int j = 0; j < child->operand_array.length; j++) {
					if (child->operand_array[j] == entry.key) {
						child->operand_array[j] = entry.value;
						entry.key->reference_count--;
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
			case nd_step_type::FALSITY_ELIMINATION:
			case nd_step_type::UNIVERSAL_INTRODUCTION:
			case nd_step_type::UNIVERSAL_ELIMINATION:
			case nd_step_type::EXISTENTIAL_INTRODUCTION:
			case nd_step_type::EXISTENTIAL_ELIMINATION:
			case nd_step_type::EQUALITY_ELIMINATION:
				for (unsigned int j = 0; j < ND_OPERAND_COUNT; j++) {
					if (child->operands[j] == entry.key) {
						child->operands[j] = entry.value;
						entry.key->reference_count--;
						entry.value->reference_count++;
					}
				}
				break;
			case nd_step_type::AXIOM:
			case nd_step_type::COMPARISON_INTRODUCTION:
			case nd_step_type::PARAMETER:
			case nd_step_type::TERM_PARAMETER:
			case nd_step_type::ARRAY_PARAMETER:
			case nd_step_type::FORMULA_PARAMETER:
			case nd_step_type::BETA_EQUIVALENCE:
			case nd_step_type::COUNT:
				error = true; break;
			}
			if (error) {
				fprintf(stderr, "transform_proofs ERROR: Invalid nd_step_type.\n");
				return false;
			}

			entry.value->children[entry.value->children.length++] = child;
			entry.key->children.length--;
		}

		free(*entry.key);
		if (entry.key->reference_count == 0)
			free(entry.key);
	}
	return true;
}

template<bool Negated, typename Formula, typename Canonicalizer>
inline bool add_ground_axiom(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const typename Formula::Term& consequent_atom,
		unsigned int constant, nd_step<Formula>* axiom)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	if (consequent_atom.type == TermType::UNARY_APPLICATION) {
		Term* atom = Term::new_apply(consequent_atom.binary.left, Term::new_constant(constant));
		if (atom == nullptr) return false;
		consequent_atom.binary.left->reference_count++;
		bool result = T.template add_unary_atom<Negated, true>(*atom, axiom);
		free(*atom); free(atom);
		return result;
	} else {
		relation rel = { consequent_atom.ternary.first->constant,
				(consequent_atom.ternary.second->type == TermType::VARIABLE ? constant : consequent_atom.ternary.second->constant),
				(consequent_atom.ternary.third->type == TermType::VARIABLE ? constant : consequent_atom.ternary.third->constant) };
		return T.template add_binary_atom<Negated, true>(rel, axiom);
	}
}

template<bool Negated, typename Formula, typename Canonicalizer>
inline bool remove_ground_axiom(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const typename Formula::Term& consequent_atom, unsigned int constant)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	if (consequent_atom.type == TermType::UNARY_APPLICATION) {
		Term* atom = Term::new_apply(consequent_atom.binary.left, Term::new_constant(constant));
		if (atom == nullptr) return false;
		consequent_atom.binary.left->reference_count++;
		bool result = T.template remove_unary_atom<Negated>(*atom);
		free(*atom); free(atom);
		return result;
	} else {
		relation rel = { consequent_atom.ternary.first->constant,
				(consequent_atom.ternary.second->type == TermType::VARIABLE ? constant : consequent_atom.ternary.second->constant),
				(consequent_atom.ternary.third->type == TermType::VARIABLE ? constant : consequent_atom.ternary.third->constant) };
		return T.template remove_binary_atom<Negated>(rel);
	}
}

template<
	typename Formula, typename Canonicalizer, bool Negated,
	typename ProofPrior, typename TheorySampleCollector>
bool do_mh_universal_intro(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const proof_transformations<Formula>& proposed_proofs,
		const array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
		const universal_intro_proposal<Negated, Formula>& proposal,
		double log_proposal_probability_ratio, ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	array<Formula*> old_extra_observations(1);
	array<Formula*> new_extra_observations(2);
	if (proposal.new_antecedent_set_size_axiom != nullptr)
		new_extra_observations[new_extra_observations.length++] = proposal.new_antecedent_set_size_axiom;
	if (proposal.new_consequent_set_size_axiom != nullptr)
		new_extra_observations[new_extra_observations.length++] = proposal.new_consequent_set_size_axiom;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			old_extra_observations, new_extra_observations,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

#if !defined(NDEBUG)
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "do_mh_universal_intro WARNING: The computed log probability ratio is NaN.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we've accepted the proposal */
		if (!remove_ground_axiom<Negated>(T, proposal.consequent_atom, proposal.concept_id))
			return false;
		if (!transform_proofs(proposed_proofs)) {
			add_ground_axiom<Negated>(T, proposal.consequent_atom, proposal.concept_id, proposal.old_axiom);
			return false;
		}
		for (auto& entry : observation_changes) {
			T.observations.remove_at(T.observations.index_of(entry.key));
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
			T.observations.add(entry.value);
			entry.value->reference_count++;
		}
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	}
	return true;
}

template<
	typename Formula, typename Canonicalizer, bool Negated,
	typename ProofPrior, typename TheorySampleCollector>
bool do_mh_universal_elim(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const proof_transformations<Formula>& proposed_proofs,
		const array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
		const universal_elim_proposal<Negated, Formula>& proposal,
		double log_proposal_probability_ratio, ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	array<Formula*> old_extra_observations(2);
	array<Formula*> new_extra_observations(1);
	if (proposal.old_antecedent_set_size_axiom != nullptr)
		old_extra_observations[old_extra_observations.length++] = proposal.old_antecedent_set_size_axiom;
	if (proposal.old_consequent_set_size_axiom != nullptr)
		old_extra_observations[old_extra_observations.length++] = proposal.old_consequent_set_size_axiom;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			old_extra_observations, new_extra_observations,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

#if !defined(NDEBUG)
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "do_mh_universal_elim WARNING: The computed log probability ratio is NaN.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we've accepted the proposal */
		if (!add_ground_axiom<Negated>(T, proposal.consequent_atom, proposal.constant, proposal.new_axiom))
			return false;
		if (!transform_proofs(proposed_proofs)) return false;

		for (auto& entry : observation_changes) {
			T.observations.remove(entry.key);
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
			T.observations.add(entry.value);
			entry.value->reference_count++;
		}
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	}
	return true;
}

template<typename Formula, typename Canonicalizer,
	typename PriorState, typename PriorStateChanges,
	typename TheorySampleCollector>
inline bool do_mh_disjunction_intro(
	theory<natural_deduction<Formula>, Canonicalizer>& T,
	pair<Formula*, nd_step<Formula>*>& selected_step,
	nd_step<Formula>* proposed_proof,
	const array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
	typename theory<natural_deduction<Formula>, Canonicalizer>::changes& old_proof_changes,
	typename theory<natural_deduction<Formula>, Canonicalizer>::changes& new_proof_changes,
	PriorState& proof_axioms,
	const PriorStateChanges& old_axioms,
	const PriorStateChanges& new_axioms,
	double log_proposal_probability_ratio,
	double log_proposal_probability_ratio_without_set_sizes,
	const undo_remove_sets& old_sets,
	const undo_remove_sets& new_sets,
	double proof_prior_diff,
	TheorySampleCollector& sample_collector)
{
#if !defined(NDEBUG)
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "do_mh_disjunction_intro WARNING: The computed log probability ratio is NaN.\n");
	if (*selected_step.value == *proposed_proof && fabs(log_proposal_probability_ratio_without_set_sizes) > 1.0e-12)
		fprintf(stderr, "do_mh_disjunction_intro WARNING: This identity proposal does not have probability ratio 1.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we accepted the new proof */
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	} else {
		return undo_proof_changes(T, old_proof_changes, new_proof_changes, selected_step.value, proposed_proof, old_sets, new_sets);
	}
	return true;
}

template<typename Formula, typename Canonicalizer>
bool is_eliminable_extensional_edge(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const nd_step<Formula>* axiom)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef nd_step<Formula> Proof;

	const Formula* formula = axiom->formula;
	if (formula->type != FormulaType::FOR_ALL
	 || formula->quantifier.operand->type != FormulaType::IF_THEN)
		return false;

	if (T.observations.contains(axiom)) return false;

	const Formula* right = formula->quantifier.operand->binary.right;

	const Formula* atom = right;
	if (atom->type == FormulaType::NOT)
		atom = atom->unary.operand;
	if (!is_atomic(*atom)) return false;

#if !defined(NDEBUG)
	bool contains;
	unsigned int set_id = T.sets.set_ids.get(*right, contains);
	if (!contains)
		fprintf(stderr, "set_reasoning.is_eliminable_extensional_edge WARNING: The given set formula is not in `set_ids`.\n");
#else
	unsigned int set_id = T.sets.set_ids.get(*right);
#endif

	for (auto entry : T.sets.extensional_graph.vertices[set_id].children) {
		if (axiom == entry.value) {
			for (Proof* child : entry.value->children) {
				if (child->type != nd_step_type::UNIVERSAL_ELIMINATION
				 || child->operands[1]->type != nd_step_type::TERM_PARAMETER
				 || child->operands[1]->term->type != TermType::CONSTANT) continue;
				for (Proof* grandchild : child->children) {
					if (grandchild->type != nd_step_type::IMPLICATION_ELIMINATION) continue;
					return true;
				}
			}
		}
	}
	return false;
}

template<typename Formula, typename Canonicalizer>
bool get_eliminable_extensional_edges(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	for (unsigned int i = 1; i < T.sets.set_count + 1; i++) {
		if (T.sets.sets[i].size_axioms.data == nullptr) continue;

		const Formula* set_formula = T.sets.sets[i].set_formula();
		if (set_formula->type == FormulaType::NOT)
			set_formula = set_formula->unary.operand;
		if (!is_atomic(*set_formula)) continue;
		for (auto entry : T.sets.extensional_graph.vertices[i].children) {
			for (Proof* axiom : entry.value) {
				/* check that this universally-quantified axiom can be removed */
				if (T.observations.contains(axiom)) continue;
				for (Proof* child : axiom->children) {
					if (child->type != nd_step_type::UNIVERSAL_ELIMINATION
					|| child->operands[1]->type != nd_step_type::TERM_PARAMETER
					|| child->operands[1]->term->type != TermType::CONSTANT) continue;
					for (Proof* grandchild : child->children) {
						if (grandchild->type != nd_step_type::IMPLICATION_ELIMINATION) continue;
						if (!eliminable_extensional_edges.add({i, entry.key, axiom, child, grandchild}))
							return false;
					}
				}
			}
		}
	}

	return true;
}

template<typename Formula, typename Canonicalizer>
bool are_mergeable(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		const typename Formula::Term* first,
		const typename Formula::Term* second,
		const relation& first_fragment,
		const relation& second_fragment)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	if (first == nullptr || second == nullptr) {
		return true;
	} else if (first->type == TermType::CONSTANT) {
		if (second->type != TermType::CONSTANT)
			return false;
		if (second->constant == first->constant)
			return true;
		Term* first_arg1 = T.get_arg(first->constant, (unsigned int) built_in_predicates::ARG1);
		Term* first_arg2 = T.get_arg(first->constant, (unsigned int) built_in_predicates::ARG2);
		Term* second_arg1 = T.get_arg(second->constant, (unsigned int) built_in_predicates::ARG1);
		Term* second_arg2 = T.get_arg(second->constant, (unsigned int) built_in_predicates::ARG2);

		if (first_arg1 != nullptr && second_arg1 != nullptr) {
			if (first_arg1->type == TermType::CONSTANT) {
				if (second_arg1->type != TermType::CONSTANT) return false;
				unsigned int first_constant = first_arg1->constant;
				if (first_arg1->constant == first_fragment.predicate)
					first_constant = second_fragment.predicate;
				else if (first_arg1->constant == first_fragment.arg1)
					first_constant = second_fragment.arg1;
				else if (first_arg1->constant == first_fragment.arg2)
					first_constant = second_fragment.arg2;
				if (first_constant != second->constant)
					return false;
			} else {
				if (first_arg1 != second_arg1 && *first_arg1 != *second_arg1)
					return false;
			}
		} if (first_arg2 != nullptr && second_arg2 != nullptr) {
			if (first_arg2->type == TermType::CONSTANT) {
				if (second_arg2->type != TermType::CONSTANT) return false;
				unsigned int first_constant = first_arg2->constant;
				if (first_arg2->constant == first_fragment.predicate)
					first_constant = second_fragment.predicate;
				else if (first_arg2->constant == first_fragment.arg2)
					first_constant = second_fragment.arg2;
				else if (first_arg2->constant == first_fragment.arg2)
					first_constant = second_fragment.arg2;
				if (first_constant != second->constant)
					return false;
			} else {
				if (first_arg2 != second_arg2 && *first_arg2 != *second_arg2)
					return false;
			}
		}
		return true;
	} else {
		return (first == second || *first == *second);
	}
}

template<typename Formula, typename Canonicalizer>
bool get_mergeable_events(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<pair<relation, relation>>& mergeable_events)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;

		Term* first_arg1 = T.get_arg(i + T.new_constant_offset, (unsigned int) built_in_predicates::ARG1);
		Term* first_arg2 = T.get_arg(i + T.new_constant_offset, (unsigned int) built_in_predicates::ARG2);

		for (unsigned int j = i + 1; j < T.ground_concept_capacity; j++) {
			if (T.ground_concepts[j].types.keys == nullptr) continue;

			/* make sure concepts `i` and `j` share at least one type */
			bool found_common_type = false;
			for (const auto& first_type : T.ground_concepts[i].types) {
				for (const auto& second_type : T.ground_concepts[j].types) {
					if (first_type.key == second_type.key) {
						found_common_type = true;
						break;
					}
				}
				if (found_common_type) break;
			}
			if (!found_common_type) continue;

			Term* second_arg1 = T.get_arg(j + T.new_constant_offset, (unsigned int) built_in_predicates::ARG1);
			Term* second_arg2 = T.get_arg(j + T.new_constant_offset, (unsigned int) built_in_predicates::ARG2);

			relation first, second;
			first.predicate = T.new_constant_offset + i;
			first.arg1 = ((first_arg1 != nullptr && first_arg1->type == TermType::CONSTANT) ? first_arg1->constant : 0);
			first.arg2 = ((first_arg2 != nullptr && first_arg2->type == TermType::CONSTANT) ? first_arg2->constant : 0);
			second.predicate = T.new_constant_offset + j;
			second.arg1 = ((second_arg1 != nullptr && second_arg1->type == TermType::CONSTANT) ? second_arg1->constant : 0);
			second.arg2 = ((second_arg2 != nullptr && second_arg2->type == TermType::CONSTANT) ? second_arg2->constant : 0);

			/* we only consider fragments with at least two vertices */
			if ((first.arg1 == 0 && first.arg2 == 0)
			 || first.arg1 == first.predicate || first.arg2 == first.predicate
			 || (first.arg1 == first.arg2 && first.arg1 != 0)
			 || !are_mergeable(T, first_arg1, second_arg1, first, second)
			 || !are_mergeable(T, first_arg2, second_arg2, first, second))
				continue;

			if (first.arg1 == 0 && second.arg1 != 0)
				first.arg1 = second.arg1;
			if (first.arg2 == 0 && second.arg2 != 0)
				first.arg2 = second.arg2;

			if (!mergeable_events.add(make_pair(first, second)))
				return false;
		}
	}
	return true;
}

template<typename Formula, typename Canonicalizer>
bool get_splittable_events(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<relation>& splittable_events)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;

		Term* arg1 = T.get_arg(i + T.new_constant_offset, (unsigned int) built_in_predicates::ARG1);
		Term* arg2 = T.get_arg(i + T.new_constant_offset, (unsigned int) built_in_predicates::ARG2);

		relation fragment;
		fragment.predicate = T.new_constant_offset + i;
		fragment.arg1 = ((arg1 != nullptr && arg1->type == TermType::CONSTANT) ? arg1->constant : 0);
		fragment.arg2 = ((arg2 != nullptr && arg2->type == TermType::CONSTANT) ? arg2->constant : 0);

		/* we only consider fragments with at least two vertices */
		if ((fragment.arg1 == 0 && fragment.arg2 == 0)
		 || fragment.arg1 == fragment.predicate
		 || fragment.arg2 == fragment.predicate
		 || (fragment.arg1 == fragment.arg2 && fragment.arg1 != 0))
			continue;

		if (!splittable_events.add(fragment))
			return false;

		array<nd_step<Formula>*> proof_roots(8);
		unsigned int predicate_instantiation_count = 0;
		unsigned int arg1_instantiation_count = 0;
		unsigned int arg2_instantiation_count = 0;
		for (const pair<Formula*, nd_step<Formula>*>& entry : T.existential_intro_nodes) {
			Term* term = entry.value->operands[2]->term;
			if (term->type != TermType::CONSTANT) continue;
			if (term->constant == fragment.predicate)
				predicate_instantiation_count++;
			if (term->constant == fragment.arg1)
				arg1_instantiation_count++;
			if (term->constant == fragment.arg2)
				arg2_instantiation_count++;
		}

		if (predicate_instantiation_count < 2) continue;
		if (fragment.arg1 != 0 && arg1_instantiation_count < 2) continue;
		if (fragment.arg2 != 0 && arg2_instantiation_count < 2) continue;
	}
	return true;
}

struct uniform_proposal {
	unsigned int mergeable_event_count;
	unsigned int splittable_event_count;
	double merge_weight;
	double split_weight;
	double log_merge_weight;
	double log_split_weight;

	double old_normalization;

	static constexpr unsigned int NewConceptProbabilityAlpha = 1;

	uniform_proposal(double merge_weight = 1.0, double split_weight = 1.0) :
			merge_weight(merge_weight), split_weight(split_weight),
			log_merge_weight(log(merge_weight)), log_split_weight(log(split_weight))
	{ }
};

template<typename Formula, typename Canonicalizer>
inline unsigned int sample(
		uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		double& log_proposal_probability_ratio)
{
	unsigned int axiom_count = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length == 0)
		fprintf(stderr, "sample WARNING: `axiom_count + mergeable_events.length + splittable_events.length` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length;
	double random = sample_uniform<double>() * normalization;
	proposal_distribution.old_normalization = normalization;
	if (random < axiom_count) {
		log_proposal_probability_ratio -= -log(normalization);
		return (unsigned int) random;
	}
	random -= axiom_count;

	if (random < proposal_distribution.merge_weight * mergeable_events.length) {
		log_proposal_probability_ratio -= proposal_distribution.log_merge_weight - log(normalization);
		return axiom_count + (unsigned int) (random / proposal_distribution.merge_weight);
	}
	random -= proposal_distribution.merge_weight * mergeable_events.length;

	if (random < proposal_distribution.split_weight * splittable_events.length) {
		log_proposal_probability_ratio -= proposal_distribution.log_split_weight - log(normalization);
		return axiom_count + mergeable_events.length + (unsigned int) (random / proposal_distribution.split_weight);
	}
	random -= proposal_distribution.split_weight * splittable_events.length;

	fprintf(stderr, "sample ERROR: Unable to sample Metropolis-Hastings proposal.\n");
	return UINT_MAX;
}

inline double log_probability(
		const uniform_proposal& proposal_distribution,
		int delta_atom_count, int delta_extensional_edges,
		int delta_unfixed_set_count)
{
	double new_normalization = proposal_distribution.old_normalization + (delta_atom_count + delta_extensional_edges + delta_unfixed_set_count);
	return -log(new_normalization);
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<Formula*, nd_step<Formula>*>& selected_disjunction_intro)
{
	unsigned int axiom_count = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length;
	return -log(normalization);
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<relation, relation>& selected_merge_event)
{
	unsigned int axiom_count = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length;
	return proposal_distribution.log_merge_weight - log(normalization);
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const relation& selected_split_event)
{
	unsigned int axiom_count = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length;
	return proposal_distribution.log_split_weight - log(normalization);
}

struct exploration_proposal {
	static constexpr unsigned int NewConceptProbabilityAlpha = 10;
};

template<typename Formula, typename Canonicalizer>
inline unsigned int sample(
		exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		double& log_proposal_probability_ratio)
{
	unsigned int axiom_count = unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
#if !defined(NDEBUG)
	if (axiom_count == 0)
		fprintf(stderr, "sample WARNING: `axiom_count` is 0.\n");
#endif
	return T.ground_axiom_count + eliminable_extensional_edges.length + sample_uniform(axiom_count);
}

inline double log_probability(
		const exploration_proposal& proposal_distribution,
		int delta_atom_count, int delta_extensional_edges,
		int delta_unfixed_set_count)
{
	return 1.0e100;
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<Formula*, nd_step<Formula>*>& selected_disjunction_intro)
{
	return 1.0e100;
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<relation, relation>& selected_merge_event)
{
	return 1.0e100;
}

template<typename Formula, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const relation& selected_split_event)
{
	return 1.0e100;
}

template<typename Formula, typename FallbackProposal>
struct query_proposal {
	double sample_query_probability;
	double log_sample_query_probability;
	double log_one_minus_sample_query_probability;
	const nd_step<Formula>* query_proof;
	bool selected_query;
	FallbackProposal& fallback_proposal;

	static constexpr unsigned int NewConceptProbabilityAlpha = 1;

	query_proposal(double sample_query_probability, const nd_step<Formula>* query_proof, FallbackProposal& fallback_proposal) :
		sample_query_probability(sample_query_probability),
		log_sample_query_probability(log(sample_query_probability)),
		log_one_minus_sample_query_probability(log(1.0 - sample_query_probability)),
		query_proof(query_proof), fallback_proposal(fallback_proposal)
	{ }
};

template<typename Formula, typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline unsigned int sample(
		query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		double& log_proposal_probability_ratio)
{
	unsigned int offset = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length;

	unsigned int query_index;
	for (query_index = 0; query_index < T.existential_intro_nodes.length; query_index++)
		if (T.existential_intro_nodes[query_index].value == proposal_distribution.query_proof) break;

	if (sample_uniform<double>() < proposal_distribution.sample_query_probability) {
		log_proposal_probability_ratio -= proposal_distribution.log_sample_query_probability;
		proposal_distribution.selected_query = true;
		return offset + query_index;
	} else {
		proposal_distribution.selected_query = false;
		return sample(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, log_proposal_probability_ratio);
	}
}

template<typename Formula, typename FallbackProposal>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		int delta_atom_count, int delta_extensional_edges, int delta_unfixed_set_count)
{
	if (proposal_distribution.selected_query)
		return proposal_distribution.log_sample_query_probability;
	else return logsumexp(proposal_distribution.log_one_minus_sample_query_probability,
			log_probability(proposal_distribution.fallback_proposal, delta_atom_count, delta_extensional_edges, delta_unfixed_set_count));
}

template<typename Formula, typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<Formula*, nd_step<Formula>*>& selected_disjunction_intro)
{
	if (selected_disjunction_intro.value == proposal_distribution.query_proof) {
		return logsumexp(proposal_distribution.log_sample_query_probability,
				log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_disjunction_intro));
	} else {
		return proposal_distribution.log_one_minus_sample_query_probability
			 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_disjunction_intro);
	}
}

template<typename Formula, typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<relation, relation>& selected_merge_event)
{
	return proposal_distribution.log_one_minus_sample_query_probability
		 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_merge_event);
}

template<typename Formula, typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const relation& selected_split_event)
{
	return proposal_distribution.log_one_minus_sample_query_probability
		 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_split_event);
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
bool do_mh_step(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		ProposalDistribution& proposal_distribution)
{
	typedef natural_deduction<Formula> ProofCalculus;
	typedef typename Formula::Term Term;

	array<unsigned int> unfixed_sets(8);
	if (!T.sets.get_unfixed_sets(unfixed_sets, T.observations)) return false;

	/* count all eliminable extensional edges */
	array<extensional_edge<Formula>> eliminable_extensional_edges(8);
	get_eliminable_extensional_edges(T, eliminable_extensional_edges);

	array<pair<relation, relation>> mergeable_events(8);
	get_mergeable_events(T, mergeable_events);

	array<relation> splittable_events(8);
	get_splittable_events(T, splittable_events);

	double log_proposal_probability_ratio = 0.0;

	/* select an axiom from `T` uniformly at random */
	unsigned int random = sample(proposal_distribution, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, log_proposal_probability_ratio);
	if (random < T.ground_axiom_count) {
		/* we've selected a grounded axiom */
		for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
			if (T.ground_concepts[i].types.keys == NULL) continue;
			unsigned int concept_id = T.new_constant_offset + i;
			concept<ProofCalculus>& c = T.ground_concepts[i];
			if (random < c.types.size)
				return propose_universal_intro<false>(T, c.types.keys[random], concept_id, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
			random -= c.types.size;
			if (random < c.negated_types.size)
				return propose_universal_intro<true>(T, c.negated_types.keys[random], concept_id, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
			random -= c.negated_types.size;
			if (random < c.relations.size) {
				relation rel = c.relations.keys[random];
				Term* atom = Term::new_apply(Term::new_constant(rel.predicate),
						(rel.arg1 == 0 ? &Term::template variables<1>::value : Term::new_constant(rel.arg1)),
						(rel.arg2 == 0 ? &Term::template variables<1>::value : Term::new_constant(rel.arg2)));
				if (atom == nullptr) return false;
				if (rel.arg1 == 0) Term::template variables<1>::value.reference_count++;
				if (rel.arg2 == 0) Term::template variables<1>::value.reference_count++;
				return propose_universal_intro<false>(T, *atom, concept_id, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
			}
			random -= c.relations.size;
			if (random < c.negated_relations.size) {
				relation rel = c.relations.keys[random];
				Term* atom = Term::new_apply(Term::new_constant(rel.predicate),
						(rel.arg1 == 0 ? &Term::template variables<1>::value : Term::new_constant(rel.arg1)),
						(rel.arg2 == 0 ? &Term::template variables<1>::value : Term::new_constant(rel.arg2)));
				if (atom == nullptr) return false;
				if (rel.arg1 == 0) Term::template variables<1>::value.reference_count++;
				if (rel.arg2 == 0) Term::template variables<1>::value.reference_count++;
				return propose_universal_intro<true>(T, *atom, concept_id, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
			}
			random -= c.negated_relations.size;
		}
	}
	random -= T.ground_axiom_count;

	if (random < eliminable_extensional_edges.length) {
		/* we've selected a universally-quantified formula */
		return propose_universal_elim(T, eliminable_extensional_edges[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= eliminable_extensional_edges.length;

	if (random < unfixed_sets.length) {
		/* we've selected to resample the size of a set */
		return propose_change_set_size(T, unfixed_sets[random], log_proposal_probability_ratio, proof_prior.axiom_prior.base_distribution.set_size_distribution, sample_collector, proposal_distribution);
	}
	random -= unfixed_sets.length;

	if (random < T.disjunction_intro_nodes.length) {
		return propose_disjunction_intro(T, T.disjunction_intro_nodes[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= T.disjunction_intro_nodes.length;
	if (random < T.negated_conjunction_nodes.length) {
		return propose_disjunction_intro(T, T.negated_conjunction_nodes[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= T.negated_conjunction_nodes.length;
	if (random < T.implication_intro_nodes.length) {
		return propose_disjunction_intro(T, T.implication_intro_nodes[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= T.implication_intro_nodes.length;
	if (random < T.existential_intro_nodes.length) {
		return propose_disjunction_intro(T, T.existential_intro_nodes[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= T.existential_intro_nodes.length;

	if (random < mergeable_events.length) {
		return propose_merge_events(T, mergeable_events[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= mergeable_events.length;
	if (random < splittable_events.length) {
		return propose_split_event(T, splittable_events[random], log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
	}
	random -= splittable_events.length;

	fprintf(stderr, "propose ERROR: Unable to select axiom.\n");
	return false;
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector>
inline bool do_mh_step(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	static thread_local uniform_proposal default_proposal(2.0, 0.001);
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, default_proposal);
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector>
inline bool do_exploratory_mh_step(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	static thread_local exploration_proposal proposal_distribution;
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

template<typename Formula, typename Canonicalizer,
	typename ProofPrior, typename TheorySampleCollector>
inline bool do_mh_step(
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		nd_step<Formula>* query_proof)
{
	static thread_local uniform_proposal default_proposal(2.0, 0.001);
	query_proposal<Formula, uniform_proposal> proposal_distribution(0.01, query_proof, default_proposal);
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
