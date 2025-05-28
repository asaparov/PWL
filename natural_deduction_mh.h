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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool propose_transformation(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		proof_transformations<Formula>& proposed_proofs,
		array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
		nd_step<Formula>* old_step, nd_step<Formula>* new_step)
{
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
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

template<bool FirstSample, typename Formula, bool Intuitionistic, typename Canonicalizer>
bool select_axiom(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		const concept<natural_deduction<Formula, Intuitionistic>>& c,
		array<pair<uint_fast8_t, unsigned int>>& axiom_indices,
		array<typename Formula::Term*>& selected_types,
		array<typename Formula::Term*>& selected_negated_types,
		array<relation>& selected_relations, array<relation>& selected_negated_relations,
		array<instance>& intersection)
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
		const array<instance>& ids = T.atoms.get(type_value).key;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 1) {
		Term& type_value = c.negated_types.keys[axiom_indices[selected_index].value];
		if (!selected_negated_types.add(&type_value)) return false;
		const array<instance>& ids = T.atoms.get(type_value).value;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 2) {
		relation_value = c.relations.keys[axiom_indices[selected_index].value];
		if (!selected_relations.add(relation_value)) return false;
		const array<instance>& ids = T.relations.get(relation_value).key;
		if (FirstSample) intersection.append(ids.data, ids.length);
		else set_intersect(intersection, ids);
	} else if (axiom_indices[selected_index].key == 3) {
		relation_value = c.negated_relations.keys[axiom_indices[selected_index].value];
		if (!selected_negated_relations.add(relation_value)) return false;
		const array<instance>& ids = T.relations.get(relation_value).value;
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

template<bool Negated, typename Formula, bool Intuitionistic, typename Canonicalizer>
inline const array<instance>& get_concept_set(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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

template<bool Negated, typename Formula, bool Intuitionistic>
nd_step<Formula>* get_proof(
		concept<natural_deduction<Formula, Intuitionistic>>& c,
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
	typename Formula, bool Intuitionistic, typename Canonicalizer>
inline void get_satisfying_concepts_helper(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		const typename Formula::Term* predicate, const typename Formula::Term* arg1,
		const typename Formula::Term* arg2, array<instance>& intersection)
{
	typedef typename Formula::TermType TermType;

	if (arg2 == NULL) {
		/* this is a unary atom */
		const pair<array<instance>, array<instance>>& list = T.atoms.get(*predicate);
		const array<instance>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist.data, sublist.length);
		else set_intersect(intersection, sublist);
	} else {
		/* this is a binary atom */
		relation r = { predicate,
				(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
				(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
		const pair<array<instance>, array<instance>>& list = T.relations.get(r);
		const array<instance>& sublist = Negated ? list.value : list.key;
		if (FirstSample) intersection.append(sublist.data, sublist.length);
		else set_intersect(intersection, sublist);
	}
}

template<bool FirstSample, typename Formula, bool Intuitionistic, typename Canonicalizer>
bool get_satisfying_concepts(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		const Formula* formula, array<instance>& intersection)
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

template<typename Formula, bool Intuitionistic>
bool get_conjunct_step(const Formula* atom,
		const concept<natural_deduction<Formula, Intuitionistic>>& instance,
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

template<bool Negated, typename Formula, bool Intuitionistic, typename Canonicalizer>
bool get_satisfying_extensional_edges(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
		unsigned int min_set_size, unsigned int max_set_size,
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
	bool Negated, typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_universal_intro(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
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
		array<nd_step<Formula>*>& proofs = T.sets.extensional_graph.vertices[selected.key].children.get(selected.value);
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
		const array<instance>& negated_set = get_concept_set<!Negated>(T, a);

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

		array<instance> intersection(32);
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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline nd_step<Formula>* make_grounded_conjunct(
		const theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		const concept<natural_deduction<Formula, Intuitionistic>>& c,
		Formula* lifted_conjunct, unsigned int variable,
		typename Formula::Term* constant)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;

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

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_universal_elim(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
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
	if (selected_edge.axiom->reference_count == 2) {
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
			set_size_proposal_log_probability(selected_edge.consequent_set, T.sets, log_prob, 0, UINT_MAX);
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
			set_size_proposal_log_probability(selected_edge.antecedent_set, T.sets, log_prob, 0, UINT_MAX);
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
	if (consequent->type == FormulaType::NOT) {
		return do_mh_universal_elim(T, proposed_proofs, observation_changes,
				make_universal_elim_proposal<true>(selected_edge, constant, new_axiom, *consequent->unary.operand, old_antecedent_set_size_axiom, old_consequent_set_size_axiom),
				log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector);
	} else {
		return do_mh_universal_elim(T, proposed_proofs, observation_changes,
				make_universal_elim_proposal<false>(selected_edge, constant, new_axiom, *consequent, old_antecedent_set_size_axiom, old_consequent_set_size_axiom),
				log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector);
	}
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename SizePrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_change_set_size(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
		T.sets.sets[selected_set].change_size(new_size);

		bool is_consistent = true;
		if (T.sets.sets[selected_set].provable_elements.length == new_size) {
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
				if (!T.template check_set_membership_after_addition<false>(conjunct)) {
					is_consistent = false;
					break;
				}
			}
		} else if (T.sets.sets[selected_set].provable_elements.length == old_size) {
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
				if (!T.check_set_membership_after_subtraction(conjunct, 0)) {
					is_consistent = false;
					break;
				}
			}
		}

		if (!is_consistent) {
			T.sets.sets[selected_set].change_size(old_size);
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

template<bool IsExploratory>
struct proof_sampler {
	array<unsigned int> removed_set_sizes;
	double log_probability;
	double set_size_log_probability;
	proof_disjunction_nodes old_proof;
	proof_disjunction_nodes new_proof;
	bool is_old_proof;
	bool undo;

	proof_sampler() : removed_set_sizes(4), log_probability(0.0), set_size_log_probability(0.0), is_old_proof(true), undo(false) { }

	inline void clear() {
		old_proof.constant_position = 0;
		old_proof.operand_position = 0;
		new_proof.clear();
		is_old_proof = true;
		log_probability = 0.0;
		set_size_log_probability = 0.0;
		removed_set_sizes.clear();
		undo = false;
	}
};

template<typename Proof, bool IsExploratory>
inline void visit_node(const Proof& proof, const proof_sampler<IsExploratory>& visitor) { }

template<bool Negated, typename Term, bool IsExploratory> constexpr bool visit_unary_atom(const Term* term, const proof_sampler<IsExploratory>& visitor) { return true; }
template<bool Negated, bool IsExploratory> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const proof_sampler<IsExploratory>& visitor) { return true; }
template<typename Proof, bool IsExploratory> constexpr bool visit_subset_axiom(const Proof& proof, const proof_sampler<IsExploratory>& visitor) { return true; }
template<bool IsExploratory> constexpr bool visit_existential_intro(const proof_sampler<IsExploratory>& visitor) { return true; }
template<bool IsExploratory> constexpr bool visit_negated_universal_intro(const proof_sampler<IsExploratory>& visitor) { return true; }
template<bool IsExploratory> constexpr bool visit_negated_conjunction(const proof_sampler<IsExploratory>& visitor) { return true; }
template<bool IsExploratory> constexpr bool visit_disjunction_intro(const proof_sampler<IsExploratory>& visitor) { return true; }

template<bool IsExploratory>
inline void on_subtract_changes(proof_sampler<IsExploratory>& visitor) {
	visitor.undo = true;
}

template<typename Formula, bool IsExploratory>
constexpr bool on_undo_filter_operands(const Formula* formula, const proof_sampler<IsExploratory>& visitor) { return true; }

template<typename Theory, typename Formula, typename Proof, bool IsExploratory>
constexpr bool on_undo_filter_constants(
		const Theory& T, const Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, const proof_sampler<IsExploratory>& visitor)
{ return true; }

template<typename Formula, bool IsExploratory>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, proof_sampler<IsExploratory>& sampler)
{
	if (!filter_operands(formula, indices))
		return false;

	if (sampler.is_old_proof) {
		if (sampler.old_proof.operand_position == sampler.old_proof.expected_operand_indices.length)
			/* this is currently possible when the new proof has more disjunction
			   introductions than the old proof, for example if the proof of the
			   expression `a(b)` chose a different set definition of `a` than in
			   the old proof */
			return false;

		unsigned int index = indices.index_of(sampler.old_proof.expected_operand_indices[sampler.old_proof.operand_position]);
		if (index != 0) {
			sampler.is_old_proof = false;
		} else {
			sampler.old_proof.operand_position++;
		}
	}

	if (!log_cache<double>::instance().ensure_size(indices.length + 1)) return false;
	sampler.log_probability += -log_cache<double>::instance().get(indices.length);

	if (!sampler.new_proof.expected_operand_indices.add(indices[0]))
		return false;

	indices.length = 1;
	return true;
}

inline double* compute_constant_probabilities(const array<instance>& constants, double& sum)
{
	sum = 0.0;
	double* probabilities = (double*) malloc(sizeof(double) * constants.length);
	if (probabilities == nullptr) {
		fprintf(stderr, "compute_constant_probabilities ERROR: Out of memory.\n");
		return nullptr;
	}
	for (unsigned int i = 0; i < constants.length; i++) {
		probabilities[i] = exp(constants[i].matching_types * max(2.0, constants.length / 40.0) - constants[i].mismatching_types * max(2.0, constants.length / 40.0));
		sum += probabilities[i];
	}

	double max_probability = 0.0;
	size_t any_index = constants.length;
	for (size_t i = 0; i < constants.length; i++) {
		max_probability = max(max_probability, probabilities[i]);
		if (constants[i].type == instance_type::ANY)
			any_index = i;
	}
	if (any_index != constants.length && probabilities[any_index] * 3 < max_probability) {
		sum -= probabilities[any_index];
		probabilities[any_index] = max_probability / 3;
		sum += probabilities[any_index];
	}
	return probabilities;
}

template<typename ProofCalculus, typename Canonicalizer>
inline double* compute_constant_probabilities(
		const array<instance>& constants, double& sum,
		const theory<ProofCalculus, Canonicalizer>& T,
		proof_disjunction_nodes& cached_proof,
		bool& is_cached_proof)
{
	double* probabilities = compute_constant_probabilities(constants, sum);

	if (is_cached_proof) {
		unsigned int index;
		const instance& key = cached_proof.expected_constants[cached_proof.constant_position];
		for (index = 0; index < constants.length; index++) {
			if (key.type == instance_type::CONSTANT && key.constant >= T.new_constant_offset
			 && T.ground_concepts[key.constant - T.new_constant_offset].types.keys == nullptr
			 && constants[index].type == instance_type::ANY)
			{
				break;
			} else if (key == constants[index]) {
				break;
			}
		}
		if (index != constants.length) {
			double min_probability = max(1.0, sum * (cached_proof.constant_position + cached_proof.operand_position));
			if (probabilities[index] < min_probability) {
				sum -= probabilities[index];
				probabilities[index] = min_probability;
				sum += probabilities[index];
			}
			cached_proof.constant_position++;
		}
	}

	return probabilities;
}

template<typename ProofCalculus, typename Canonicalizer, bool IsExploratory>
inline bool get_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions,
		proof_sampler<IsExploratory>& sampler)
{
	if (!get_possible_constants(T, constants))
		return false;

	array<pair<unsigned int, bool>> new_name_events(8);
	array<const typename ProofCalculus::Language*> arg1_of(4);
	array<const typename ProofCalculus::Language*> arg2_of(4);
	if (!filter_constants_helper<false>(T, formula, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of))
		return false;

	if (IsExploratory) {
		unsigned int index = index_of_any(constants);
		if (index < constants.length)
			constants[index].matching_types += 1;
	}

	double sum;
	double* probabilities = compute_constant_probabilities(constants, sum, T, sampler.old_proof, sampler.is_old_proof);
	if (probabilities == nullptr)
		return false;

	unsigned int random = sample_categorical(probabilities, sum, constants.length);
	sampler.log_probability += log(probabilities[random] / sum);
	free(probabilities);
	swap(constants[0], constants[random]);

	if (!sampler.new_proof.expected_constants.add(constants[0]))
		return false;

	constants.length = 1;
	return true;
}

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer, bool IsExploratory>
inline constexpr bool is_impossible(
		const typename ProofCalculus::Language* formula,
		const theory<ProofCalculus, Canonicalizer>& T,
		const proof_sampler<IsExploratory>& sampler)
{
	return false;
}

template<typename Formula, bool IsExploratory>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, proof_sampler<IsExploratory>& sampler) { return true; }

template<typename Formula, bool IsExploratory>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, proof_sampler<IsExploratory>& sampler) { return true; }

template<typename Formula, bool IsExploratory>
inline void finished_constants(const Formula* formula, proof_sampler<IsExploratory>& sampler) { }

template<typename Formula, bool IsExploratory>
inline void finished_operand_indices(const Formula* formula, proof_sampler<IsExploratory>& sampler) { }

struct inverse_proof_sampler {
	array<unsigned int> removed_set_sizes;
	double set_size_log_probability;
	array<array<unsigned int>> operand_indices;
	array<pair<array<instance>, unsigned int>> constants;

	inverse_proof_sampler() : removed_set_sizes(4), set_size_log_probability(0.0), operand_indices(4), constants(4) { }

	~inverse_proof_sampler() {
		for (array<unsigned int>& indices : operand_indices)
			free(indices);
		for (pair<array<instance>, unsigned int>& entry : constants)
			free(entry.key);
	}
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
inline bool get_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions,
		inverse_proof_sampler& sampler)
{
	if (!get_possible_constants(T, constants))
		return false;
	shuffle(constants);
	return true;
}

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer>
inline constexpr bool is_impossible(
		const typename ProofCalculus::Language* formula,
		const theory<ProofCalculus, Canonicalizer>& T,
		const inverse_proof_sampler& sampler)
{
	return false;
}

template<typename Formula>
inline bool on_undo_filter_operands(Formula* formula, inverse_proof_sampler& sampler) {
	typedef typename Formula::Type FormulaType;

	if (!sampler.operand_indices.ensure_capacity(sampler.operand_indices.length + 1)
	 || !array_init(sampler.operand_indices[sampler.operand_indices.length], 8))
		return false;
	array<unsigned int>& indices = sampler.operand_indices[sampler.operand_indices.length];
	sampler.operand_indices.length++;

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

	return filter_operands(formula, indices);
}

template<typename Theory, typename Formula, typename Proof>
inline bool on_undo_filter_constants(
		Theory& T, Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, inverse_proof_sampler& sampler)
{
	typedef typename Formula::TermType TermType;

	if (!sampler.constants.ensure_capacity(sampler.constants.length + 1)
	 || !array_init(sampler.constants[sampler.constants.length].key, T.ground_concept_capacity + T.constant_types.size + T.constant_negated_types.size + 1))
		return false;
	array<instance>& constants = sampler.constants[sampler.constants.length].key;
	unsigned int& index = sampler.constants[sampler.constants.length].value;
	sampler.constants.length++;
	if (!get_possible_constants(T, constants))
		return false;

	array<pair<unsigned int, bool>> new_name_events(8);
	array<const Formula*> arg1_of(4);
	array<const Formula*> arg2_of(4);
	if (!filter_constants_helper<false>(T, quantified, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of))
		return false;

	bool is_new_concept = (term->type == TermType::CONSTANT && term->constant >= T.new_constant_offset && T.ground_concepts[term->constant - T.new_constant_offset].types.keys == nullptr);
	for (index = 0; index < constants.length; index++) {
		if (is_new_concept && constants[index].type == instance_type::ANY)
			break;
		if ((constants[index].type == instance_type::CONSTANT && term->type == TermType::CONSTANT && constants[index].constant == term->constant)
		 || (constants[index].type == instance_type::NUMBER && term->type == TermType::NUMBER && constants[index].number == term->number)
		 || (constants[index].type == instance_type::STRING && term->type == TermType::STRING && *constants[index].str == term->str))
		{
			break;
		}
	}

	if (index == constants.length) {
		/* NOTE: this could happen if an existential quantifier is instantiated
		   as an integer, which was previously the size of a set, but a
		   subsequent MH transition has changed the size of that set, and so
		   the original integer is not a member of `constants` */
		for (index = 0; index < constants.length; index++)
			if (constants[index].type == instance_type::ANY) break;
	}
	return true;
}

template<typename Theory>
double proof_sample_log_probability(const inverse_proof_sampler& sampler, const Theory& T, proof_disjunction_nodes& new_proof)
{
	double log_probability = 0.0;
	bool is_new_proof = true;
	for (unsigned int i = sampler.operand_indices.length; i > 0; i--) {
		const array<unsigned int>& indices = sampler.operand_indices[i - 1];

		if (is_new_proof) {
			unsigned int index = indices.index_of(new_proof.expected_operand_indices[new_proof.operand_position]);
			if (index != 0) {
				is_new_proof = false;
			} else {
				new_proof.operand_position++;
			}
		}

		log_cache<double>::instance().ensure_size(indices.length + 1);
		log_probability += -log_cache<double>::instance().get(indices.length);
	}
	for (unsigned int i = sampler.constants.length; i > 0; i--) {
		const array<instance>& constants = sampler.constants[i - 1].key;

		double sum;
		double* probabilities = compute_constant_probabilities(constants, sum, T, new_proof, is_new_proof);
		if (probabilities == nullptr)
			return false;

		log_probability += log(probabilities[sampler.constants[i - 1].value] / sum);
		free(probabilities);
	}
	return log_probability;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, inverse_proof_sampler& sampler) { return true; }

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, inverse_proof_sampler& sampler) { return true; }

template<typename Formula>
inline void finished_constants(const Formula* formula, inverse_proof_sampler& sampler) { }

template<typename Formula>
inline void finished_operand_indices(const Formula* formula, inverse_proof_sampler& sampler) { }

struct undo_remove_sets {
	array<unsigned int>& removed_set_sizes;

	undo_remove_sets(array<unsigned int>& removed_set_sizes) : removed_set_sizes(removed_set_sizes) { }
};

inline void on_subtract_changes(const undo_remove_sets& visitor) { }

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const undo_remove_sets& visitor) { return true; }

template<typename Theory, typename Formula, typename Proof>
constexpr bool on_undo_filter_constants(const Theory& T,
		const Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, const undo_remove_sets& visitor)
{ return true; }

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
		unsigned int min_set_size, unsigned int max_set_size,
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

template<typename ProofCalculus>
bool transform_proofs(const proof_transformations<typename ProofCalculus::Language>& proposed_proofs)
{
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
			case nd_step_type::INEQUALITY_INTRODUCTION:
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

template<bool FreeProposedProofs, typename Formula, bool Intuitionistic, typename Canonicalizer>
bool undo_proof_changes(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::changes& old_proof_changes,
		typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::changes& new_proof_changes,
		nd_step<Formula>* old_proof, nd_step<Formula>* new_proof,
		proof_transformations<Formula>& proposed_proofs,
		const undo_remove_sets& old_sets, const undo_remove_sets& new_sets)
{
	typedef nd_step<Formula> Proof;

	array<pair<Proof*, Proof*>> observation_changes(8);
	proof_transformations<Formula>& inverse_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(inverse_proofs)) {
		if (FreeProposedProofs) free(proposed_proofs);
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
	} else if (!propose_transformation(T, inverse_proofs, observation_changes, new_proof, old_proof)) {
		free(inverse_proofs);
		if (FreeProposedProofs) free(proposed_proofs);
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
	} else if (!transform_proofs<natural_deduction<Formula, Intuitionistic>>(inverse_proofs)) {
		free(inverse_proofs);
		if (FreeProposedProofs) free(proposed_proofs);
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
	}

	for (auto& entry : observation_changes) {
		unsigned int index = T.observations.index_of(entry.key);
		T.observations[index] = entry.value;
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		entry.value->reference_count++;
	}
	free(inverse_proofs);

	if (FreeProposedProofs)
		free(proposed_proofs);

	set_changes<Formula> set_diff;
	if (!T.subtract_changes(new_proof_changes, set_diff, new_sets)) {
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
	}
	free(new_proof_changes);
	free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);

	/* add the changes back into T from `old_proof` */
	if (!T.add_changes(old_proof_changes, set_diff, old_sets)) {
		free(old_proof_changes);
		fprintf(stderr, "undo_proof_changes ERROR: Failed to add changes back to theory.\n");
		exit(0);
		return false;
	}
	free(old_proof_changes);
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void set_size_proposal_log_probability(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		double& log_probability_value, unsigned int min_set_size, unsigned int max_set_size)
{
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

template<typename Theory, typename Formula, typename Proof>
constexpr bool on_undo_filter_constants(
		const Theory& T, const Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, const inverse_set_size_log_probability& visitor)
{ return true; }

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
		unsigned int min_set_size, unsigned int max_set_size,
		inverse_set_size_log_probability& visitor)
{
	set_size_proposal_log_probability(set_id, sets, visitor.value, min_set_size, max_set_size);
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

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, bool IsExploratory>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		proof_sampler<IsExploratory>& sampler)
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

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, bool IsExploratory>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int min_set_size, unsigned int max_set_size,
		proof_sampler<IsExploratory>& sampler)
{
	if (sampler.undo) return;
	sampler.removed_set_sizes.add(sets.sets[set_id].set_size);
}

template<typename Proof, bool IsExploratory>
constexpr bool on_new_size_axiom(Proof* new_size_axiom,
		const proof_sampler<IsExploratory>& visitor)
{
	return true;
}

template<typename Proof, bool IsExploratory>
inline void on_old_size_axiom(Proof* old_size_axiom,
		const proof_sampler<IsExploratory>& visitor)
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
		unsigned int min_set_size, unsigned int max_set_size,
		inverse_proof_sampler& visitor)
{
	set_size_proposal_log_probability(set_id, sets, visitor.set_size_log_probability, min_set_size, max_set_size);
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

template<bool IsExploratory, typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename PriorState,
	typename PriorStateChanges, typename TheorySampleCollector>
inline bool do_mh_disjunction_intro(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
	typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::proof_node& selected_step,
	nd_step<Formula>* proposed_proof,
	proof_transformations<Formula>& proposed_proofs,
	const array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
	typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::changes& old_proof_changes,
	typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::changes& new_proof_changes,
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
	if (!IsExploratory && *selected_step.proof == *proposed_proof && fabs(log_proposal_probability_ratio_without_set_sizes) > 1.0e-12)
		fprintf(stderr, "do_mh_disjunction_intro WARNING: This identity proposal does not have probability ratio 1.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we accepted the new proof */
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		free(proposed_proofs);
		free(new_proof_changes);
		free(*proposed_proof); if (proposed_proof->reference_count == 0) free(proposed_proof);
		free(old_proof_changes);
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	} else {
		return undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, selected_step.proof, proposed_proof, proposed_proofs, old_sets, new_sets);
	}
	return true;
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_disjunction_intro(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
	typename theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>::proof_node& selected_step,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	typedef nd_step<Formula> Proof;
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
	typedef theory<ProofCalculus, Canonicalizer> Theory;
	typedef typename Theory::proof_node ProofNode;

	/* check to make sure this proof wouldn't remove any subset edges, since we cannot make the inverse proposal */
	if (has_proof_step<nd_step_type::IMPLICATION_ELIMINATION>(selected_step.proof))
		return true;

	ProofNode selected_proof_step(selected_step);

	set_changes<Formula> set_diff;
	inverse_proof_sampler inverse_sampler;
	typename Theory::changes& old_proof_changes = *((typename Theory::changes*) alloca(sizeof(typename Theory::changes)));
	if (!Theory::init(old_proof_changes)) {
		return false;
	} else if (!T.get_theory_changes(*selected_proof_step.proof, old_proof_changes)) {
		free(old_proof_changes);
		return false;
	}
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
	proof_sampler<ProposalDistribution::IsExploratory> sampler;
	set_changes<Formula> new_set_diff;
	get_proof_disjunction_nodes(selected_proof_step.proof, sampler.old_proof);
unsigned int debug = 0;
if (debug_flag) {
print(*selected_proof_step.formula, stderr); print('\n', stderr);
}
	while (true) {
		/* sample a new proof, only selecting one path at every branch,
		   and avoiding paths that we've previously proved to be inconsistent,
		   also compute the log probability of the new path */
		unsigned int new_constant = 0;
		sampler.clear();
if (debug_flag) {
fprintf(stderr, "INNER DEBUG: %u\n", debug);
T.template print_axioms<true>(stderr, *debug_terminal_printer);
}
debug++;
		new_set_diff.clear();
		array_map<unsigned int, unsigned int> requested_set_sizes(4);
		unsigned int old_set_definition_count = selected_proof_step.set_definitions.size;
		new_proof = T.template make_proof<false, true, false>(selected_proof_step.formula, selected_proof_step.set_definitions, requested_set_sizes, new_set_diff, new_constant, sampler);
		selected_proof_step.set_definitions.size = old_set_definition_count;
		if (new_proof != NULL)
			break;
	}
if (debug_flag) {
fprintf(stderr, "INNER DEBUG: %u (loop broken)\n", debug);
T.template print_axioms<true>(stderr, *debug_terminal_printer);
}

	log_proposal_probability_ratio -= sampler.log_probability;
	log_proposal_probability_ratio -= sampler.set_size_log_probability;
	log_proposal_probability_ratio += proof_sample_log_probability(inverse_sampler, T, sampler.new_proof);
	log_proposal_probability_ratio += inverse_sampler.set_size_log_probability;

	typename Theory::changes& new_proof_changes = *((typename Theory::changes*) alloca(sizeof(typename Theory::changes)));
	if (!Theory::init(new_proof_changes)) {
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes); return false;
	} else if (!T.get_theory_changes(*new_proof, new_proof_changes)) {
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
	}

	/* check if the proposed proof is the same as the original proof */
	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (*selected_proof_step.proof == *new_proof) {
		undo_proof_changes<false>(T, old_proof_changes, new_proof_changes, selected_proof_step.proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
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

	/* propose `new_proof` to substitute `selected_proof_step.proof` */
	array<pair<Proof*, Proof*>> observation_changes(8);
	if (!init(proposed_proofs)) {
		proposed_proofs.transformed_proofs.keys = nullptr;
		undo_proof_changes<false>(T, old_proof_changes, new_proof_changes, selected_proof_step.proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		return false;
	} else if (!propose_transformation(T, proposed_proofs, observation_changes, selected_proof_step.proof, new_proof)) {
		undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, selected_proof_step.proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(proposed_proofs); return false;
	}

	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			set_diff.old_set_axioms, set_diff.new_set_axioms,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

	if (!transform_proofs<ProofCalculus>(proposed_proofs)) {
		undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, selected_proof_step.proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*selected_proof_step.proof); if (selected_proof_step.proof->reference_count == 0) free(selected_proof_step.proof);
		return false;
	}

	for (auto& entry : observation_changes) {
		unsigned int index = T.observations.index_of(entry.key);
		T.observations[index] = entry.value;
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
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

	log_proposal_probability_ratio += log_probability(proposal_distribution, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_proof_step.proof);

	return do_mh_disjunction_intro<ProposalDistribution::IsExploratory>(
			T, selected_proof_step, new_proof, proposed_proofs, observation_changes, old_proof_changes,
			new_proof_changes, proof_axioms, old_axioms, new_axioms, log_proposal_probability_ratio,
			log_proposal_probability_ratio + sampler.set_size_log_probability - inverse_sampler.set_size_log_probability,
			undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes),
			proof_prior_diff, sample_collector);
}

struct proof_initializer {
	proof_disjunction_nodes nodes;
	unsigned int operand_position;
	array<unsigned int> removed_set_sizes;
	double set_size_log_probability;
	bool undo;

	proof_initializer() : removed_set_sizes(4), set_size_log_probability(0.0), undo(false) { }

	static inline void free(proof_initializer& initializer) {
		core::free(initializer.nodes);
		core::free(initializer.removed_set_sizes);
	}
};

inline bool init(proof_initializer& initializer) {
	if (!init(initializer.nodes)) {
		return false;
	} else if (!array_init(initializer.removed_set_sizes, 4)) {
		free(initializer.nodes);
		return false;
	}
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

template<typename Theory, typename Formula, typename Proof>
constexpr bool on_undo_filter_constants(
		const Theory& T, const Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, const proof_initializer& visitor)
{ return true; }

template<typename Formula>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, proof_initializer& initializer)
{
	if (!filter_operands(formula, indices))
		return false;

	if (initializer.nodes.operand_position == initializer.nodes.expected_operand_indices.length) {
		fprintf(stderr, "filter_operands ERROR: `proof_initializer` has no further `expected_operand_indices`.\n");
		exit(EXIT_FAILURE);
	} else if (!indices.contains(initializer.nodes.expected_operand_indices[initializer.nodes.operand_position])) {
		indices.length = 0;
		return false;
	}
	indices[0] = initializer.nodes.expected_operand_indices[initializer.nodes.operand_position++];
	indices.length = 1;
	return true;
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool get_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions,
		proof_initializer& initializer)
{
	if (initializer.nodes.constant_position == initializer.nodes.expected_constants.length)
		/* this is currently possible when the new proof has more existential
		   introductions than the old proof, for example if the proof of the
		   expression `a(b)` chose a different set definition of `a` than in
		   the old proof */
		return false;

	constants[0] = initializer.nodes.expected_constants[initializer.nodes.constant_position];
	constants.length = 1;

	array<pair<unsigned int, bool>> new_name_events(8);
	array<const typename ProofCalculus::Language*> arg1_of(4);
	array<const typename ProofCalculus::Language*> arg2_of(4);
	if (!filter_constants_helper<false>(T, formula, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of))
		return false;

	initializer.nodes.constant_position++;
	return true;
}

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer>
inline constexpr bool is_impossible(
		const typename ProofCalculus::Language* formula,
		const theory<ProofCalculus, Canonicalizer>& T,
		const proof_initializer& initializer)
{
	return false;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, proof_initializer& initializer) {
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, proof_initializer& initializer) {
	return true;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, proof_initializer& initializer) {
	initializer.nodes.constant_position--;
}

template<typename Formula>
inline void finished_operand_indices(const Formula* formula, proof_initializer& initializer) {
	initializer.nodes.operand_position--;
}

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
		unsigned int min_set_size, unsigned int max_set_size,
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

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename DstEventType,
	typename MapConstantsFunction,
	typename ComputeLogProbabilityFunction,
	typename ProofPrior, typename TheorySampleCollector,
	typename ProposalDistribution>
inline bool do_split_merge(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	typedef theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer> Theory;
	typedef typename Theory::changes TheoryChanges;
	typedef typename Theory::proof_node ProofNode;
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
	typedef typename ProofCalculus::Proof Proof;

	array<ProofNode> old_proofs(8);
	for (unsigned int i = 0; i < T.existential_intro_nodes.length; i++) {
		ProofNode& entry = T.existential_intro_nodes[i];
		Term* term = entry.proof->operands[2]->term;
		if (term->type != TermType::CONSTANT) continue;
		if (term->constant != src_event.predicate && term->constant != src_event.arg1 && term->constant != src_event.arg2)
			continue;

		bool is_root = true;
		array<Proof*> stack(max(8, entry.proof->children.length));
		for (Proof* child : entry.proof->children)
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

			for (Proof* child : current->children) {
				if (!stack.add(child)) {
					for (ProofNode& node : old_proofs) node.~ProofNode();
					return false;
				}
			}
		}

		if (!is_root) continue;
		if (!old_proofs.ensure_capacity(old_proofs.length + 1)) {
			for (ProofNode& node : old_proofs) node.~ProofNode();
			return false;
		}
		new (&old_proofs[old_proofs.length]) ProofNode(entry);
		old_proofs.length++;
	}

	/* check to make sure these proofs wouldn't remove any subset edges, since we cannot make the inverse proposal */
	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (has_proof_step<nd_step_type::IMPLICATION_ELIMINATION>(old_proofs[i].proof)) {
			for (ProofNode& node : old_proofs) node.~ProofNode();
			return true;
		}
	}

	TheoryChanges* old_proof_changes = (TheoryChanges*) malloc(sizeof(TheoryChanges) * old_proofs.length);
	if (old_proof_changes == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (ProofNode& node : old_proofs) node.~ProofNode();
		return false;
	}

	TheoryChanges* new_proof_changes = (TheoryChanges*) malloc(sizeof(TheoryChanges) * old_proofs.length);
	if (new_proof_changes == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!Theory::init(old_proof_changes[i])) {
			for (unsigned int j = 0; j < i; j++) free(old_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes); return false;
		}
	} for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!Theory::init(new_proof_changes[i])) {
			for (unsigned int j = 0; j < i; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes); return false;
		}
	}

	proof_initializer* initializers = (proof_initializer*) malloc(sizeof(proof_initializer) * old_proofs.length);
	if (initializers == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); free(new_proof_changes); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!init(initializers[i])) {
			for (unsigned int j = 0; j < i; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes); free(initializers); return false;
		}
	}

	Proof** new_proofs = (Proof**) malloc(sizeof(Proof*) * old_proofs.length);
	if (new_proofs == nullptr) {
		fprintf(stderr, "do_split_merge ERROR: Out of memory.\n");
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); free(new_proof_changes); free(initializers); return false;
	}

	set_changes<Formula> set_diff;
	inverse_set_size_log_probability set_size_log_probability;
	array_map<const Proof*, unsigned int> reference_counts(32);
	for (unsigned int i = 0; i < old_proofs.length; i++) {
		if (!get_proof_disjunction_nodes(old_proofs[i].proof, initializers[i].nodes)
		 || !reference_counts.ensure_capacity(reference_counts.size + 1))
		{
			set_changes<Formula> dummy;
			for (unsigned int j = i; j > 0; j--) {
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
				free(old_proof_changes[j - 1]);
			}
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (unsigned int j = i; j < old_proofs.length; j++) free(old_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes);
			free(initializers); free(new_proofs); return false;
		}

		for (instance& inst : initializers[i].nodes.expected_constants)
			map_constants(inst);

		array<const Proof*> discharged_axioms(16);
		unsigned int index = reference_counts.index_of(old_proofs[i].proof);
		if (index == reference_counts.size) {
			reference_counts.keys[index] = old_proofs[i].proof;
			reference_counts.values[index] = old_proofs[i].proof->reference_count;
			reference_counts.size++;
		}
		reference_counts.values[index]--;
		array_map<unsigned int, unsigned int> requested_set_sizes(4);
		if (!T.get_theory_changes(*old_proofs[i].proof, discharged_axioms, reference_counts, requested_set_sizes, old_proof_changes[i])) {
			for (ProofNode& node : old_proofs) node.~ProofNode();
			return false;
		}
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
		for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			free(old_proof_changes[j - 1]);
		}
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); free(new_proof_changes);
		free(initializers); free(new_proofs); return false;
	}

	for (unsigned int i = 0; i < old_proofs.length; i++) {
		unsigned int new_constant = 0;
		set_changes<Formula> new_set_diff;
		array_map<unsigned int, unsigned int> requested_set_sizes(4);
		unsigned int old_set_definition_count = old_proofs[i].set_definitions.size;
		new_proofs[i] = T.template make_proof<false, true, false>(old_proofs[i].formula, old_proofs[i].set_definitions, requested_set_sizes, new_set_diff, new_constant, initializers[i]);
		old_proofs[i].set_definitions.size = old_set_definition_count;
		if (new_proofs[i] == nullptr) {
			free(proposed_proofs);
			set_changes<Formula> dummy;
			for (unsigned int j = i; j > 0; j--) {
				T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
				free(new_proof_changes[j - 1]);
				free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
			} for (unsigned int j = old_proofs.length; j > 0; j--) {
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
				free(old_proof_changes[j - 1]);
			}
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = i; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes);
			free(initializers); free(new_proofs); return true;
		}

		log_proposal_probability_ratio -= initializers[i].set_size_log_probability;

		if (!T.get_theory_changes(*new_proofs[i], new_proof_changes[i])) {
			free(proposed_proofs);
			set_changes<Formula> dummy;
			free(*new_proofs[i]); if (new_proofs[i]->reference_count == 0) free(new_proofs[i]);
			for (unsigned int j = i; j > 0; j--) {
				T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
				free(new_proof_changes[j - 1]);
				free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
			} for (unsigned int j = old_proofs.length; j > 0; j--) {
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
				free(old_proof_changes[j - 1]);
			}
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = i; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes);
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

		/* propose `new_proofs[i]` to substitute `old_proofs[i].proof` */
		if (!propose_transformation(T, proposed_proofs, observation_changes, old_proofs[i].proof, new_proofs[i])) {
			free(proposed_proofs);
			set_changes<Formula> dummy;
			for (unsigned int j = i + 1; j > 0; j--) {
				T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
				free(new_proof_changes[j - 1]);
				free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
			} for (unsigned int j = old_proofs.length; j > 0; j--) {
				T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
				free(old_proof_changes[j - 1]);
			}
			for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
			for (unsigned int j = i + 1; j < old_proofs.length; j++) free(new_proof_changes[j]);
			for (ProofNode& node : old_proofs) node.~ProofNode();
			free(old_proof_changes); free(new_proof_changes);
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

	if (!transform_proofs<ProofCalculus>(proposed_proofs)) {
		free(proposed_proofs);
		set_changes<Formula> dummy;
		for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
			free(new_proof_changes[j - 1]);
			free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
		} for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			free(old_proof_changes[j]);
		}
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); free(new_proof_changes);
		free(initializers); free(new_proofs); return false;
	}

	auto undo_proof_changes = [&]() {
		array<pair<Proof*, Proof*>> observation_changes(8);
		proof_transformations<Formula>& inverse_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
		if (!init(inverse_proofs))
			return false;
		for (unsigned int i = 0; i < old_proofs.length; i++) {
			if (!propose_transformation(T, inverse_proofs, observation_changes, new_proofs[i], old_proofs[i].proof)) {
				free(inverse_proofs);
				return false;
			}
		}
		if (!transform_proofs<ProofCalculus>(inverse_proofs)) {
			free(inverse_proofs);
			return false;
		}

		for (auto& entry : observation_changes) {
			unsigned int index = T.observations.index_of(entry.key);
			T.observations[index] = entry.value;
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
			entry.value->reference_count++;
		}
		free(inverse_proofs);

		free(proposed_proofs);
		set_changes<Formula> dummy;
		for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.subtract_changes(new_proof_changes[j - 1], dummy, undo_remove_sets(initializers[j - 1].removed_set_sizes));
			free(new_proof_changes[j - 1]);
			free(*new_proofs[j - 1]); if (new_proofs[j - 1]->reference_count == 0) free(new_proofs[j - 1]);
		} for (unsigned int j = old_proofs.length; j > 0; j--) {
			T.add_changes(old_proof_changes[j - 1], dummy, undo_remove_sets(set_size_log_probability.removed_set_sizes));
			free(old_proof_changes[j - 1]);
		}
		for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
		for (ProofNode& node : old_proofs) node.~ProofNode();
		free(old_proof_changes); free(new_proof_changes);
		free(initializers); free(new_proofs);
		return true;
	};

	for (auto& entry : observation_changes) {
		unsigned int index = T.observations.index_of(entry.key);
		T.observations[index] = entry.value;
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		entry.value->reference_count++;
	}

	/* count all unfixed sets */
	array<unsigned int> unfixed_sets(8);
	if (!T.sets.get_unfixed_sets(unfixed_sets, T.observations)) {
		for (ProofNode& node : old_proofs) node.~ProofNode();
		return false;
	}

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
		if (!proof_axioms.add(new_axioms)) {
			for (ProofNode& node : old_proofs) node.~ProofNode();
			return false;
		}
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes)) {
			for (ProofNode& node : old_proofs) node.~ProofNode();
			return false;
		}
	} else {
		return undo_proof_changes();
	}

	free(proposed_proofs);
	for (unsigned int i = 0; i < old_proofs.length; i++) {
		free(*new_proofs[i]); if (new_proofs[i]->reference_count == 0) free(new_proofs[i]);
	}
	for (unsigned int j = 0; j < old_proofs.length; j++) free(initializers[j]);
	for (unsigned int j = 0; j < old_proofs.length; j++) free(new_proof_changes[j]);
	for (unsigned int j = 0; j < old_proofs.length; j++) free(old_proof_changes[j]);
	for (ProofNode& node : old_proofs) node.~ProofNode();
	free(old_proof_changes); free(new_proof_changes);
	free(initializers); free(new_proofs);
	return true;
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_merge_events(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
		for (const auto& node : T.existential_intro_nodes) {
			Term* term = node.proof->operands[2]->term;
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

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_split_event(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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

template<bool Negated, typename Formula,
	bool Intuitionistic, typename Canonicalizer>
inline bool add_ground_axiom(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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

template<bool Negated, typename Formula, bool Intuitionistic, typename Canonicalizer>
inline bool remove_ground_axiom(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	typename Formula, bool Intuitionistic,
	typename Canonicalizer, bool Negated,
	typename ProofPrior, typename TheorySampleCollector>
bool do_mh_universal_intro(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
		if (!transform_proofs<natural_deduction<Formula, Intuitionistic>>(proposed_proofs)) {
			add_ground_axiom<Negated>(T, proposal.consequent_atom, proposal.concept_id, proposal.old_axiom);
			return false;
		}
		for (auto& entry : observation_changes) {
			unsigned int index = T.observations.index_of(entry.key);
			T.observations[index] = entry.value;
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
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

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, bool Negated,
	typename ProofPrior, typename TheorySampleCollector>
bool do_mh_universal_elim(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		proof_transformations<Formula>& proposed_proofs,
		const array<pair<nd_step<Formula>*, nd_step<Formula>*>>& observation_changes,
		const universal_elim_proposal<Negated, Formula>& proposal,
		double log_proposal_probability_ratio, ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	typedef typename Formula::Type FormulaType;

	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	array<Formula*> old_extra_observations(4);
	array<Formula*> new_extra_observations(1);
	if (proposal.old_antecedent_set_size_axiom != nullptr)
		old_extra_observations[old_extra_observations.length++] = proposal.old_antecedent_set_size_axiom;
	if (proposal.old_consequent_set_size_axiom != nullptr)
		old_extra_observations[old_extra_observations.length++] = proposal.old_consequent_set_size_axiom;

	if (proposal.edge.axiom->reference_count == 2) {
		/* removing this universal quantifier could also remove axioms from the proof of the antecedent */
		if (proposal.edge.grandchild->operands[1]->type == nd_step_type::CONJUNCTION_INTRODUCTION) {
			for (nd_step<Formula>* conjunct : proposal.edge.grandchild->operands[1]->operand_array) {
				if (conjunct->type == nd_step_type::AXIOM && conjunct->reference_count == 2
				 && conjunct->formula->type == FormulaType::UNARY_APPLICATION
				 && conjunct->formula->binary.right->type == FormulaType::CONSTANT
				 && conjunct->formula->binary.right->constant >= T.new_constant_offset)
				{
					old_extra_observations.add(conjunct->formula);
				}
			}
		} else {
			nd_step<Formula>* conjunct = proposal.edge.grandchild->operands[1];
			if (conjunct->type == nd_step_type::AXIOM && conjunct->reference_count == 2
			 && conjunct->formula->type == FormulaType::UNARY_APPLICATION
			 && conjunct->formula->binary.right->type == FormulaType::CONSTANT
			 && conjunct->formula->binary.right->constant >= T.new_constant_offset)
			{
				old_extra_observations.add(conjunct->formula);
			}
		}
	}

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
		if (!add_ground_axiom<Negated>(T, proposal.consequent_atom, proposal.constant, proposal.new_axiom)
		 || !transform_proofs<natural_deduction<Formula, Intuitionistic>>(proposed_proofs))
		{
			free(*proposal.new_axiom); if (proposal.new_axiom->reference_count == 0) free(proposal.new_axiom);
			free(proposed_proofs);
			if (proposal.edge.axiom->reference_count == 1)
				T.sets.free_subset_axiom(proposal.edge.axiom);
			return false;
		}

		for (Formula* old_formula : old_extra_observations) {
			if (old_formula->type == FormulaType::UNARY_APPLICATION) {
				old_formula->reference_count++;
				T.template remove_unary_atom<false>(*old_formula);
				free(*old_formula); if (old_formula->reference_count == 0) free(old_formula);
			} else if (old_formula->type == FormulaType::NOT && old_formula->unary.operand->type == FormulaType::UNARY_APPLICATION) {
				old_formula->reference_count++;
				T.template remove_unary_atom<true>(*old_formula->unary.operand);
				free(*old_formula); if (old_formula->reference_count == 0) free(old_formula);
			}
		}

		for (auto& entry : observation_changes) {
			unsigned int index = T.observations.index_of(entry.key);
			T.observations[index] = entry.value;
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
			entry.value->reference_count++;
		}
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms)) {
			free(*proposal.new_axiom); if (proposal.new_axiom->reference_count == 0) free(proposal.new_axiom);
			free(proposed_proofs);
			if (proposal.edge.axiom->reference_count == 1)
				T.sets.free_subset_axiom(proposal.edge.axiom);
			return false;
		}
		free(*proposal.new_axiom); if (proposal.new_axiom->reference_count == 0) free(proposal.new_axiom);
		free(proposed_proofs);
		if (proposal.edge.axiom->reference_count == 1)
			T.sets.free_subset_axiom(proposal.edge.axiom);

		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		return sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes);
	} else {
		free(*proposal.new_axiom); if (proposal.new_axiom->reference_count == 0) free(proposal.new_axiom);
		free(proposed_proofs);
		if (proposal.edge.axiom->reference_count == 1)
			T.sets.free_subset_axiom(proposal.edge.axiom);
		return true;
	}
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool propose_rebind_anaphora(
	theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
	unsigned int sentence_id, unsigned int anaphora_id,
	double& log_proposal_probability_ratio,
	ProofPrior& proof_prior,
	typename ProofPrior::PriorState& proof_axioms,
	TheorySampleCollector& sample_collector,
	const ProposalDistribution& proposal_distribution)
{
	typedef nd_step<Formula> Proof;
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
	typedef theory<ProofCalculus, Canonicalizer> Theory;

	unsigned int anaphora_start = T.ctx.referent_iterators[sentence_id].anaphora.values[anaphora_id];
	if (anaphora_start <= 1)
		/* there is no more than one possible referent for this anaphora */
		return true;

	/* check to make sure this proof wouldn't remove any subset edges, since we cannot make the inverse proposal */
	Proof* proof = T.observations[sentence_id];
	if (has_proof_step<nd_step_type::IMPLICATION_ELIMINATION>(proof))
		return true;

	/* first remove the observation with ID given by `sentence_id` */
	set_changes<Formula> set_diff;
	inverse_proof_sampler inverse_sampler;
	typename Theory::changes& old_proof_changes = *((typename Theory::changes*) alloca(sizeof(typename Theory::changes)));
	if (!Theory::init(old_proof_changes)) {
		return false;
	} else if (!T.get_theory_changes(*proof, old_proof_changes)) {
		free(old_proof_changes);
		return false;
	}
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

	/* try re-binding the selected anaphora */
	Proof* new_proof;
	proof_sampler<ProposalDistribution::IsExploratory> sampler;
	set_changes<Formula> new_set_diff;
	unsigned int old_referent_id = T.ctx.bindings[sentence_id].indices[anaphora_id];
	while (true) {
		unsigned int new_referent_id = sample_uniform(anaphora_start);
		T.ctx.bindings[sentence_id].indices[anaphora_id] = anaphora_start - new_referent_id - 1;
		if (T.ctx.bindings[sentence_id].indices[anaphora_id] == old_referent_id) {
			/* for performance, we assume that if the proposed anaphora binding
			   is the same as the old one, we just propose the identity transformation */
			set_changes<Formula> set_diff;
			if (!T.add_changes(old_proof_changes, set_diff, undo_remove_sets(inverse_sampler.removed_set_sizes))) {
				free(old_proof_changes);
				fprintf(stderr, "propose_rebind_anaphora ERROR: Failed to add changes back to theory.\n");
				exit(0);
				return false;
			}
			free(old_proof_changes);
			return true;
		}
		if (!is_anaphora_binding_legal(T.ctx.referent_iterators[sentence_id], T.ctx.bindings[sentence_id].indices))
			continue;

		hol_term* resolved_logical_form = bind_anaphora(T.ctx.referent_iterators[sentence_id], T.ctx.bindings[sentence_id].indices);
		if (resolved_logical_form == nullptr)
			continue;

		unsigned int new_constant = 0;
		sampler.clear();
		new_set_diff.clear();
		new_proof = T.add_formula_helper(resolved_logical_form, new_set_diff, new_constant, sampler);
		core::free(*resolved_logical_form);
		if (resolved_logical_form->reference_count == 0)
			core::free(resolved_logical_form);
		if (new_proof != nullptr)
			break;
	}

	/* compute the log probability of the new theory */
	log_proposal_probability_ratio -= sampler.log_probability;
	log_proposal_probability_ratio -= sampler.set_size_log_probability;
	log_proposal_probability_ratio += proof_sample_log_probability(inverse_sampler, T, sampler.new_proof);
	log_proposal_probability_ratio += inverse_sampler.set_size_log_probability;

	typename Theory::changes& new_proof_changes = *((typename Theory::changes*) alloca(sizeof(typename Theory::changes)));
	if (!Theory::init(new_proof_changes)) {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes); return false;
	} else if (!T.get_theory_changes(*new_proof, new_proof_changes)) {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		return false;
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

	/* propose `new_proof` to substitute `proof` */
	array<pair<Proof*, Proof*>> observation_changes(8);
	proof_transformations<Formula>& proposed_proofs = *((proof_transformations<Formula>*) alloca(sizeof(proof_transformations<Formula>)));
	if (!init(proposed_proofs)) {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		proposed_proofs.transformed_proofs.keys = nullptr;
		undo_proof_changes<false>(T, old_proof_changes, new_proof_changes, proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		return false;
	} else if (!propose_transformation(T, proposed_proofs, observation_changes, proof, new_proof)) {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(proposed_proofs); return false;
	}

	/* compute the proof portion of the prior for both current and proposed theories */
	typename ProofPrior::PriorStateChanges old_axioms, new_axioms;
	double proof_prior_diff = log_probability_ratio(proposed_proofs.transformed_proofs,
			set_diff.old_set_axioms, set_diff.new_set_axioms,
			proof_prior, proof_axioms, old_axioms, new_axioms, sample_collector);
	log_proposal_probability_ratio += proof_prior_diff;

	if (!transform_proofs<ProofCalculus>(proposed_proofs)) {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
		free(*proof); if (proof->reference_count == 0) free(proof);
		return false;
	}

	for (auto& entry : observation_changes) {
		unsigned int index = T.observations.index_of(entry.key);
		T.observations[index] = entry.value;
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
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

	log_proposal_probability_ratio += log_probability(proposal_distribution, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, proof);

#if !defined(NDEBUG)
	double log_proposal_probability_ratio_without_set_sizes = log_proposal_probability_ratio + sampler.set_size_log_probability - inverse_sampler.set_size_log_probability;
	if (isnan(log_proposal_probability_ratio))
		fprintf(stderr, "propose_rebind_anaphora WARNING: The computed log probability ratio is NaN.\n");
	if (!ProposalDistribution::IsExploratory && set_diff.old_set_axioms.length == 0 && set_diff.new_set_axioms.length == 0 && *proof == *new_proof && fabs(log_proposal_probability_ratio_without_set_sizes) > 1.0e-12)
		fprintf(stderr, "propose_rebind_anaphora WARNING: This identity proposal does not have probability ratio 1.\n");
#endif

	if (sample_uniform<double>() < exp(log_proposal_probability_ratio)) {
		/* we accepted the new proof */
		proof_axioms.subtract(old_axioms);
		if (!proof_axioms.add(new_axioms))
			return false;
		free(proposed_proofs);
		free(new_proof_changes);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		free(old_proof_changes);
		/* TODO: we could keep track of the extra axioms within `theory` */
		array<hol_term*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		if (!sample_collector.accept_with_observation_changes(T.observations, extra_axioms, proof_prior_diff, observation_changes))
			return false;
	} else {
		T.ctx.bindings[sentence_id].indices[anaphora_id] = old_referent_id;
		return undo_proof_changes<true>(T, old_proof_changes, new_proof_changes, proof, new_proof, proposed_proofs, undo_remove_sets(inverse_sampler.removed_set_sizes), undo_remove_sets(sampler.removed_set_sizes));
	}
	return true;
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool is_eliminable_extensional_edge(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool get_eliminable_extensional_edges(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool are_mergeable(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		const typename Formula::Term* first,
		const typename Formula::Term* second,
		const relation& first_fragment,
		const relation& second_fragment,
		const array<typename Formula::Term*>& first_name_terms,
		const array<typename Formula::Term*>& second_name_terms,
		const typename Formula::Term* first_arg1,
		const typename Formula::Term* first_arg2)
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
		Term* second_arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(second->constant);
		Term* second_arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(second->constant);

		if ((first_arg1 != nullptr && *first_arg1 == *second)
		 || (first_arg2 != nullptr && *first_arg2 == *second)
		 || (second_arg1 != nullptr && *second_arg1 == *first)
		 || (second_arg2 != nullptr && *second_arg2 == *first))
		{
			/* prevent cases where an event is an argument of itself */
			return false;
		}

		/* do not consider merging concepts if it will create a concept with more names */
		unsigned int a = 0, b = 0;
		while (a < first_name_terms.length && b < second_name_terms.length) {
			if (*first_name_terms[a] == *second_name_terms[b]) {
				a++; b++;
			} else {
				break;
			}
		}
		if (a < first_name_terms.length || b < second_name_terms.length)
			return false;

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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool get_mergeable_events(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<pair<relation, relation>>& mergeable_events)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;

		Term* first_arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(i + T.new_constant_offset);
		Term* first_arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(i + T.new_constant_offset);

		array<Term*> first_name_terms(4);
		array<Term*> first_arg1_name_terms(4);
		array<Term*> first_arg2_name_terms(4);
		if (!T.get_concept_names(i + T.new_constant_offset, first_name_terms)
		 || (first_arg1 != nullptr && first_arg1->type == TermType::CONSTANT && !T.get_concept_names(first_arg1->constant, first_arg1_name_terms))
		 || (first_arg2 != nullptr && first_arg2->type == TermType::CONSTANT && !T.get_concept_names(first_arg2->constant, first_arg2_name_terms)))
			return false;
		if (first_name_terms.length > 1)
			sort(first_name_terms, pointer_sorter());
		if (first_arg1_name_terms.length > 1)
			sort(first_arg1_name_terms, pointer_sorter());
		if (first_arg2_name_terms.length > 1)
			sort(first_arg2_name_terms, pointer_sorter());

		array<unsigned int> other_concepts(8);
		for (const auto& first_type : T.ground_concepts[i].types) {
			const array<instance>& constants = T.atoms.get(first_type.key).key;
			for (const instance& constant : constants) {
				if (constant.type != instance_type::CONSTANT || constant.constant < T.new_constant_offset || constant.constant == i + T.new_constant_offset)
					continue;

				/* only consider pairs of concepts that share at least one type */
				if (!other_concepts.add(constant.constant))
					return false;
			}
		}

		if (other_concepts.length > 1) {
			sort(other_concepts);
			unique(other_concepts);
		}

		Term* arg1_of_first_arg1;
		Term* arg2_of_first_arg1;
		Term* arg1_of_first_arg2;
		Term* arg2_of_first_arg2;
		if (first_arg1 != nullptr && first_arg1->type == TermType::CONSTANT) {
			arg1_of_first_arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(first_arg1->constant);
			arg2_of_first_arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(first_arg1->constant);
		} else {
			arg1_of_first_arg1 = nullptr;
			arg2_of_first_arg1 = nullptr;
		} if (first_arg2 != nullptr && first_arg2->type == TermType::CONSTANT) {
			arg1_of_first_arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(first_arg2->constant);
			arg2_of_first_arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(first_arg2->constant);
		} else {
			arg1_of_first_arg2 = nullptr;
			arg2_of_first_arg2 = nullptr;
		}

		for (unsigned int j : other_concepts) {
			/* do not consider other concepts `j` if it results in a new concept with more names */
			array<Term*> second_name_terms(4);
			if (!T.get_concept_names(j, second_name_terms))
				return false;
			if (second_name_terms.length > 1)
				sort(second_name_terms, pointer_sorter());
			unsigned int a = 0, b = 0;
			while (a < first_name_terms.length && b < second_name_terms.length) {
				if (*first_name_terms[a] == *second_name_terms[b]) {
					a++; b++;
				} else {
					break;
				}
			}
			if (a < first_name_terms.length || b < second_name_terms.length)
				continue;

			Term* second_arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(j);
			Term* second_arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(j);

			array<Term*> second_arg1_name_terms(4);
			array<Term*> second_arg2_name_terms(4);
			if ((second_arg1 != nullptr && second_arg1->type == TermType::CONSTANT && !T.get_concept_names(second_arg1->constant, second_arg1_name_terms))
			 || (second_arg2 != nullptr && second_arg2->type == TermType::CONSTANT && !T.get_concept_names(second_arg2->constant, second_arg2_name_terms)))
				return false;
			if (second_arg1_name_terms.length > 1)
				sort(second_arg1_name_terms, pointer_sorter());
			if (second_arg2_name_terms.length > 1)
				sort(second_arg2_name_terms, pointer_sorter());

			relation first, second;
			first.predicate = T.new_constant_offset + i;
			first.arg1 = ((first_arg1 != nullptr && first_arg1->type == TermType::CONSTANT) ? first_arg1->constant : 0);
			first.arg2 = ((first_arg2 != nullptr && first_arg2->type == TermType::CONSTANT) ? first_arg2->constant : 0);
			second.predicate = j;
			second.arg1 = ((second_arg1 != nullptr && second_arg1->type == TermType::CONSTANT) ? second_arg1->constant : 0);
			second.arg2 = ((second_arg2 != nullptr && second_arg2->type == TermType::CONSTANT) ? second_arg2->constant : 0);

			/* we only consider fragments with at least two vertices */
			if ((first.arg1 == 0 && first.arg2 == 0)
			 || first.arg1 == first.predicate || first.arg2 == first.predicate
			 || (first.arg1 == first.arg2 && first.arg1 != 0)
			 || !are_mergeable(T, first_arg1, second_arg1, first, second, first_arg1_name_terms, second_arg1_name_terms, arg1_of_first_arg1, arg2_of_first_arg1)
			 || !are_mergeable(T, first_arg2, second_arg2, first, second, first_arg2_name_terms, second_arg2_name_terms, arg1_of_first_arg2, arg2_of_first_arg2))
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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
bool get_splittable_events(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<relation>& splittable_events)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;

		Term* arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(i + T.new_constant_offset);
		Term* arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(i + T.new_constant_offset);

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
		for (const auto& node : T.existential_intro_nodes) {
			Term* term = node.proof->operands[2]->term;
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

	static constexpr bool IsExploratory = false;

	uniform_proposal(double merge_weight = 1.0, double split_weight = 1.0) :
			merge_weight(merge_weight), split_weight(split_weight),
			log_merge_weight(log(merge_weight)), log_split_weight(log(split_weight))
	{ }
};

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline unsigned int sample(
		uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	unsigned int anaphora_count = 0;
	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++)
		anaphora_count += T.ctx.referent_iterators[i].anaphora.size;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length + anaphora_count == 0)
		fprintf(stderr, "sample WARNING: `axiom_count + mergeable_events.length + splittable_events.length + anaphora_count` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length + anaphora_count;
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

	if (random < anaphora_count) {
		log_proposal_probability_ratio -= -log(normalization);
		return axiom_count + mergeable_events.length + splittable_events.length + (unsigned int) random;
	}
	random -= anaphora_count;

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

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const nd_step<Formula>* selected_disjunction_intro)
{
	unsigned int axiom_count = T.ground_axiom_count
			+ eliminable_extensional_edges.length + unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;
	unsigned int anaphora_count = 0;
	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++)
		anaphora_count += T.ctx.referent_iterators[i].anaphora.size;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length + anaphora_count == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length + anaphora_count` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length + anaphora_count;
	return -log(normalization);
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	unsigned int anaphora_count = 0;
	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++)
		anaphora_count += T.ctx.referent_iterators[i].anaphora.size;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length + anaphora_count == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length + anaphora_count` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length + anaphora_count;
	return proposal_distribution.log_merge_weight - log(normalization);
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const uniform_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	unsigned int anaphora_count = 0;
	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++)
		anaphora_count += T.ctx.referent_iterators[i].anaphora.size;
#if !defined(NDEBUG)
	if (axiom_count + mergeable_events.length + splittable_events.length + anaphora_count == 0)
		fprintf(stderr, "log_probability WARNING: `axiom_count + mergeable_events.length + splittable_events.length + anaphora_count` is 0.\n");
#endif
	double normalization = axiom_count + proposal_distribution.merge_weight * mergeable_events.length + proposal_distribution.split_weight * splittable_events.length + anaphora_count;
	return proposal_distribution.log_split_weight - log(normalization);
}

struct exploration_proposal {
	static constexpr bool IsExploratory = true;
};

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline unsigned int sample(
		exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		double& log_proposal_probability_ratio)
{
	unsigned int axiom_count = unfixed_sets.length
			+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
			+ T.implication_intro_nodes.length + T.existential_intro_nodes.length;

	unsigned int anaphora_count = 0;
	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++)
		anaphora_count += T.ctx.referent_iterators[i].anaphora.size;

#if !defined(NDEBUG)
	if (axiom_count + anaphora_count == 0)
		fprintf(stderr, "sample WARNING: `axiom_count + anaphora_count` is 0.\n");
#endif

	unsigned int random = sample_uniform(axiom_count + anaphora_count);
	if (random < axiom_count) {
		return T.ground_axiom_count + eliminable_extensional_edges.length + sample_uniform(axiom_count);
	} else {
		random -= axiom_count;
		unsigned int anaphora_offset = T.ground_axiom_count
				+ eliminable_extensional_edges.length + unfixed_sets.length
				+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
				+ T.implication_intro_nodes.length + T.existential_intro_nodes.length
				+ mergeable_events.length + splittable_events.length;
		return anaphora_offset + random;
	}
}

inline double log_probability(
		const exploration_proposal& proposal_distribution,
		int delta_atom_count, int delta_extensional_edges,
		int delta_unfixed_set_count)
{
	return 1.0e100;
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const nd_step<Formula>* selected_disjunction_intro)
{
	return 1.0e100;
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<extensional_edge<Formula>>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<relation, relation>& selected_merge_event)
{
	return 1.0e100;
}

template<typename Formula, bool Intuitionistic, typename Canonicalizer>
inline double log_probability(
		const exploration_proposal& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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
	unsigned int query_step_count;
	FallbackProposal& fallback_proposal;

	static constexpr bool IsExploratory = false;

	query_proposal(double sample_query_probability, const nd_step<Formula>* query_proof, FallbackProposal& fallback_proposal) :
		sample_query_probability(sample_query_probability),
		log_sample_query_probability(log(sample_query_probability)),
		log_one_minus_sample_query_probability(log(1.0 - sample_query_probability)),
		query_proof(query_proof), fallback_proposal(fallback_proposal)
	{ }
};

template<typename Formula, bool Intuitionistic,
	typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline unsigned int sample(
		query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
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

	unsigned int observation_index;
	for (observation_index = 0; observation_index < T.observations.length; observation_index++)
		if (T.observations[observation_index] == proposal_distribution.query_proof) break;
	unsigned int num_anaphora = T.ctx.referent_iterators[observation_index].anaphora.size;

	array<const nd_step<Formula>*> proof_steps(4);
	get_proof_steps<nd_step_type::EXISTENTIAL_INTRODUCTION>(proposal_distribution.query_proof, proof_steps);
	proposal_distribution.query_step_count = proof_steps.length + num_anaphora;

	if (proof_steps.length != 0 && sample_uniform<double>() < proposal_distribution.sample_query_probability) {
		log_cache<double>::instance().ensure_size(proof_steps.length + num_anaphora + 1);
		log_proposal_probability_ratio -= proposal_distribution.log_sample_query_probability - log_cache<double>::instance().get(proof_steps.length + num_anaphora);
		proposal_distribution.selected_query = true;
		unsigned int random = sample_uniform(proof_steps.length + num_anaphora);
		if (random < proof_steps.length) {
			const nd_step<Formula>* selected_step = proof_steps[random];

			unsigned int query_index;
			for (query_index = 0; query_index < T.existential_intro_nodes.length; query_index++)
				if (T.existential_intro_nodes[query_index].proof == selected_step) break;
			return offset + query_index;
		} else {
			unsigned int anaphora_offset = T.ground_axiom_count
					+ eliminable_extensional_edges.length + unfixed_sets.length
					+ T.disjunction_intro_nodes.length + T.negated_conjunction_nodes.length
					+ T.implication_intro_nodes.length + T.existential_intro_nodes.length
					+ mergeable_events.length + splittable_events.length;
			for (unsigned int i = 0; i < observation_index; i++)
				anaphora_offset += T.ctx.referent_iterators[i].anaphora.size;
			random -= proof_steps.length;
			return anaphora_offset + random;
		}
	} else {
		proposal_distribution.selected_query = false;
		return proposal_distribution.log_one_minus_sample_query_probability + sample(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, log_proposal_probability_ratio);
	}
}

template<typename Formula, typename FallbackProposal>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		int delta_atom_count, int delta_extensional_edges, int delta_unfixed_set_count)
{
	if (proposal_distribution.selected_query)
		return proposal_distribution.log_sample_query_probability - log_cache<double>::instance().get(proposal_distribution.query_step_count);
	else return logsumexp(proposal_distribution.log_one_minus_sample_query_probability,
			log_probability(proposal_distribution.fallback_proposal, delta_atom_count, delta_extensional_edges, delta_unfixed_set_count));
}

template<typename Formula, bool Intuitionistic,
	typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const nd_step<Formula>* selected_disjunction_intro)
{
	if (proposal_distribution.selected_query) {
		return logsumexp(proposal_distribution.log_sample_query_probability - log_cache<double>::instance().get(proposal_distribution.query_step_count),
				log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_disjunction_intro));
	} else {
		return proposal_distribution.log_one_minus_sample_query_probability
			 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_disjunction_intro);
	}
}

template<typename Formula, bool Intuitionistic,
	typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const pair<relation, relation>& selected_merge_event)
{
	return proposal_distribution.log_one_minus_sample_query_probability
		 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_merge_event);
}

template<typename Formula, bool Intuitionistic,
	typename FallbackProposal, typename Canonicalizer, typename ExtensionalEdge>
inline double log_probability(
		const query_proposal<Formula, FallbackProposal>& proposal_distribution,
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		array<ExtensionalEdge>& eliminable_extensional_edges,
		array<unsigned int>& unfixed_sets,
		array<pair<relation, relation>>& mergeable_events,
		array<relation>& splittable_events,
		const relation& selected_split_event)
{
	return proposal_distribution.log_one_minus_sample_query_probability
		 + log_probability(proposal_distribution.fallback_proposal, T, eliminable_extensional_edges, unfixed_sets, mergeable_events, splittable_events, selected_split_event);
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector, typename ProposalDistribution>
bool do_mh_step(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		ProposalDistribution& proposal_distribution)
{
	typedef natural_deduction<Formula, Intuitionistic> ProofCalculus;
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

	for (unsigned int i = 0; i < T.ctx.referent_iterators.length; i++) {
		if (random < T.ctx.referent_iterators[i].anaphora.size) {
			/* propose resampling the j-th anaphora in the i-th sentence, where `j = random` */
			return propose_rebind_anaphora(T, i, random, log_proposal_probability_ratio, proof_prior, proof_axioms, sample_collector, proposal_distribution);
		}
		random -= T.ctx.referent_iterators[i].anaphora.size;
	}

	fprintf(stderr, "propose ERROR: Unable to select axiom.\n");
	return false;
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector>
inline bool do_mh_step(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	static thread_local uniform_proposal default_proposal(2.0, 0.001);
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, default_proposal);
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector>
inline bool do_exploratory_mh_step(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector)
{
	static thread_local exploration_proposal proposal_distribution;
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector>
inline bool do_mh_step(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		nd_step<Formula>* query_proof,
		double query_proposal_probability = 0.01)
{
	static thread_local uniform_proposal default_proposal(2.0, 0.001);
	query_proposal<Formula, uniform_proposal> proposal_distribution(query_proposal_probability, query_proof, default_proposal);
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

template<typename Formula, bool Intuitionistic,
	typename Canonicalizer, typename ProofPrior,
	typename TheorySampleCollector>
inline bool do_exploratory_mh_step(
		theory<natural_deduction<Formula, Intuitionistic>, Canonicalizer>& T,
		ProofPrior& proof_prior,
		typename ProofPrior::PriorState& proof_axioms,
		TheorySampleCollector& sample_collector,
		nd_step<Formula>* query_proof,
		double query_proposal_probability = 0.01)
{
	static thread_local exploration_proposal default_proposal;
	query_proposal<Formula, exploration_proposal> proposal_distribution(query_proposal_probability, query_proof, default_proposal);
	return do_mh_step(T, proof_prior, proof_axioms, sample_collector, proposal_distribution);
}

#endif /* NATURAL_DEDUCTION_MH_H_ */
