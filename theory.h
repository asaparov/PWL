#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <math/multiset.h>

#include "array_view.h"
#include "set_reasoning.h"

using namespace core;


enum class built_in_predicates : unsigned int {
	ZERO = 0,
	UNKNOWN,
	ARG1,
	ARG2,
	ARG3,
	ARG1_OF,
	ARG2_OF,
	ARG3_OF,
	SIZE,
	INVERSE,
	OWN,
	EXIST,
	PRESENT,
	PRESENT_PROGRESSIVE,
	PRESENT_PERFECT,
	PRESENT_PERFECT_PROGRESSIVE,
	PAST,
	PAST_PROGRESSIVE,
	PAST_PERFECT,
	PAST_PERFECT_PROGRESSIVE,
	FUTURE,
	FUTURE_PROGRESSIVE,
	FUTURE_PERFECT,
	FUTURE_PERFECT_PROGRESSIVE,
	EMPTY,
	EMPTY_REF,
	WIDE_SCOPE,

	EQUALS, /* this is only used to refer to set definitions of the form `A=^[x]:f(x)` */
	SUBSET,
	SAME,

	COUNT
};

/* WARNING: The below should preserve the order of the entries in the enum. */

unsigned int PAST_OR_PRESENT[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PAST
};

unsigned int PRESENT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE
};

unsigned int FUTURE_PREDICATES[] = {
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int PERFECT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int NON_PERFECT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE
};

unsigned int PROGRESSIVE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int NON_PROGRESSIVE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT
};

unsigned int TENSE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("<ZERO>", (unsigned int) built_in_predicates::ZERO)
		&& names.put("unknown", (unsigned int) built_in_predicates::UNKNOWN)
		&& names.put("arg1", (unsigned int) built_in_predicates::ARG1)
		&& names.put("arg2", (unsigned int) built_in_predicates::ARG2)
		&& names.put("arg3", (unsigned int) built_in_predicates::ARG3)
		&& names.put("arg1_of", (unsigned int) built_in_predicates::ARG1_OF)
		&& names.put("arg2_of", (unsigned int) built_in_predicates::ARG2_OF)
		&& names.put("arg3_of", (unsigned int) built_in_predicates::ARG3_OF)
		&& names.put("inverse", (unsigned int) built_in_predicates::INVERSE)
		&& names.put("own", (unsigned int) built_in_predicates::OWN)
		&& names.put("exist", (unsigned int) built_in_predicates::EXIST)
		&& names.put("size", (unsigned int) built_in_predicates::SIZE)
		&& names.put("present", (unsigned int) built_in_predicates::PRESENT)
		&& names.put("present_progressive", (unsigned int) built_in_predicates::PRESENT_PROGRESSIVE)
		&& names.put("present_perfect", (unsigned int) built_in_predicates::PRESENT_PERFECT)
		&& names.put("present_perfect_progressive", (unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE)
		&& names.put("past", (unsigned int) built_in_predicates::PAST)
		&& names.put("past_progressive", (unsigned int) built_in_predicates::PAST_PROGRESSIVE)
		&& names.put("past_perfect", (unsigned int) built_in_predicates::PAST_PERFECT)
		&& names.put("past_perfect_progressive", (unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE)
		&& names.put("future", (unsigned int) built_in_predicates::FUTURE)
		&& names.put("future_progressive", (unsigned int) built_in_predicates::FUTURE_PROGRESSIVE)
		&& names.put("future_perfect", (unsigned int) built_in_predicates::FUTURE_PERFECT)
		&& names.put("future_perfect_progressive", (unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE)
		&& names.put("empty", (unsigned int) built_in_predicates::EMPTY)
		&& names.put("empty_ref", (unsigned int) built_in_predicates::EMPTY_REF)
		&& names.put("W", (unsigned int) built_in_predicates::WIDE_SCOPE)
		&& names.put("=", (unsigned int) built_in_predicates::EQUALS)
		&& names.put("subset", (unsigned int) built_in_predicates::SUBSET)
		&& names.put("same", (unsigned int) built_in_predicates::SAME);
}

template<typename Proof>
inline void visit_node(const Proof& proof) { }

template<bool Negated>
constexpr bool visit_unary_atom(unsigned int predicate, unsigned int arg) { return true; }

template<bool Negated>
constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2) { return true; }

template<typename Proof>
constexpr bool visit_subset_axiom(const Proof& proof) { return true; }

constexpr bool visit_existential_intro() { return true; }
constexpr bool visit_negated_conjunction() { return true; }
constexpr bool visit_negated_universal_intro() { return true; }
constexpr bool visit_disjunction_intro() { return true; }

struct relation {
	unsigned int predicate;
	unsigned int arg1; /* `0` here indicates the source vertex */
	unsigned int arg2; /* `0` here indicates the source vertex */

	relation() { }

	relation(unsigned int predicate, unsigned int arg1, unsigned int arg2) :
		predicate(predicate), arg1(arg1), arg2(arg2)
	{ }

	static inline bool is_empty(const relation& key) {
		return key.predicate == 0;
	}

	static inline unsigned int hash(const relation& key) {
		return default_hash(key.predicate) ^ default_hash(key.arg1) ^ default_hash(key.arg2);
	}

	static inline void move(const relation& src, relation& dst) {
		dst.predicate = src.predicate;
		dst.arg1 = src.arg1;
		dst.arg2 = src.arg2;
	}
};

inline bool operator == (const relation& first, const relation& second) {
	return first.predicate == second.predicate
		&& first.arg1 == second.arg1
		&& first.arg2 == second.arg2;
}

inline bool operator != (const relation& first, const relation& second) {
	return first.predicate != second.predicate
		|| first.arg1 != second.arg1
		|| first.arg2 != second.arg2;
}

template<typename ProofCalculus>
struct concept
{
	typedef typename ProofCalculus::Proof Proof;

	array_map<unsigned int, Proof*> types;
	array_map<unsigned int, Proof*> negated_types;
	array_map<relation, Proof*> relations;
	array_map<relation, Proof*> negated_relations;

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (auto entry : types) {
			if (!print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_types) {
			if (!print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : relations) {
			if (!print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_relations) {
			if (!print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		}
		return true;
	}

	static inline void move(const concept<ProofCalculus>& src, concept<ProofCalculus>& dst) {
		core::move(src.types, dst.types);
		core::move(src.negated_types, dst.negated_types);
		core::move(src.relations, dst.relations);
		core::move(src.negated_relations, dst.negated_relations);
	}

	static inline void free(concept<ProofCalculus>& c) {
		for (auto entry : c.types) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.negated_types) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.relations) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.negated_relations) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}
		core::free(c.types);
		core::free(c.negated_types);
		core::free(c.relations);
		core::free(c.negated_relations);
	}
};

template<typename ProofCalculus>
inline bool init(concept<ProofCalculus>& c) {
	if (!array_map_init(c.types, 8)) {
		return false;
	} else if (!array_map_init(c.negated_types, 8)) {
		free(c.types); return false;
	} else if (!array_map_init(c.relations, 8)) {
		free(c.negated_types);
		free(c.types); return false;
	} else if (!array_map_init(c.negated_relations, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); return false;
	}
	return true;
}

template<typename Formula>
inline Formula* preprocess_formula(Formula* src) {
	Formula* temp = same_to_equals(src);
	if (temp == nullptr) return nullptr;

	Formula* next = remove_exists(temp);
	free(*temp); if (temp->reference_count == 0) free(temp);
	return next;
}

template<typename Formula,
	typename ProofCalculus,
	typename Canonicalizer>
struct theory
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;

	unsigned int new_constant_offset;

	/* A map from `x` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `x(y_i)` and for any `z_i` there is an axiom in the theory
	   `~x(z_i)`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `x(u)` or `~x(u)`
	   are in the theory. */
	hash_map<unsigned int, pair<array<unsigned int>, array<unsigned int>>> types;

	/* A map from `R` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]R` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]R`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]R` or `~[u/0]R`
	   are in the theory. */
	hash_map<relation, pair<array<unsigned int>, array<unsigned int>>> relations;

	concept<ProofCalculus>* ground_concepts;
	unsigned int ground_concept_capacity;
	unsigned int ground_axiom_count;

	hash_set<Proof*> observations;
	hash_multiset<Formula*, false> proof_axioms;
	set_reasoning<built_in_predicates, Formula, ProofCalculus> sets;

	array<pair<Formula*, Proof*>> disjunction_intro_nodes;
	array<pair<Formula*, Proof*>> negated_conjunction_nodes;
	array<pair<Formula*, Proof*>> implication_intro_nodes;
	array<pair<Formula*, Proof*>> existential_intro_nodes;

	Proof* empty_set_axiom;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), types(64), relations(64),
			ground_concept_capacity(64), ground_axiom_count(0), observations(64), proof_axioms(64),
			disjunction_intro_nodes(16), negated_conjunction_nodes(16),
			implication_intro_nodes(16), existential_intro_nodes(16)
	{
		ground_concepts = (concept<ProofCalculus>*) malloc(sizeof(concept<ProofCalculus>) * ground_concept_capacity);
		if (ground_concepts == NULL) {
			fprintf(stderr, "theory ERROR: Insufficient memory for `ground_concepts`.\n");
			exit(EXIT_FAILURE);
		}
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			ground_concepts[i].types.keys = NULL; /* this is used to indicate that this concept is uninitialized */

		Formula* empty_set_formula = Formula::new_for_all(1, Formula::new_equals(
			Formula::new_equals(Formula::new_atom((unsigned int) built_in_predicates::SIZE, Formula::new_variable(1)), Formula::new_int(0)),
			Formula::new_not(Formula::new_exists(2, Formula::new_apply(Formula::new_variable(1), Formula::new_variable(2))))
		));
		if (empty_set_formula == NULL) {
			core::free(ground_concepts);
			exit(EXIT_FAILURE);
		}
		empty_set_axiom = ProofCalculus::new_axiom(empty_set_formula);
		core::free(*empty_set_formula);
		if (empty_set_formula->reference_count == 0)
			core::free(empty_set_formula);
		if (empty_set_axiom == NULL) exit(EXIT_FAILURE);
		empty_set_axiom->reference_count++;
	}

	~theory() {
		for (auto entry : types) {
			core::free(entry.value.key);
			core::free(entry.value.value);
		} for (auto entry : relations) {
			core::free(entry.value.key);
			core::free(entry.value.value);
		} for (Proof* proof : observations) {
			core::free(*proof);
			if (proof->reference_count == 0)
				core::free(proof);
		} for (auto entry : disjunction_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : negated_conjunction_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : implication_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : existential_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		}

		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			if (ground_concepts[i].types.keys != NULL) core::free(ground_concepts[i]);
		core::free(ground_concepts);

		core::free(*empty_set_axiom);
		if (empty_set_axiom->reference_count == 0)
			core::free(empty_set_axiom);
	}

	unsigned int get_free_concept_id() {
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			if (ground_concepts[i].types.keys == NULL) return i + new_constant_offset;

		unsigned int new_capacity = ground_concept_capacity;
		expand_capacity(new_capacity, ground_concept_capacity + 1);

		if (!resize(ground_concepts, new_capacity))
			return EXIT_FAILURE;
		for (unsigned int i = ground_concept_capacity; i < new_capacity; i++)
			ground_concepts[i].types.keys = NULL; /* this is used to indicate that this concept is uninitialized */
		ground_concept_capacity = new_capacity;
		return ground_concept_capacity + new_constant_offset;
	}

	inline bool try_init_concept(unsigned int id) {
		if (id < new_constant_offset) return true;
		concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.keys != NULL)
			return true;
		else return init(c);
	}

	void free_concept_id(unsigned int id) {
#if !defined(NDEBUG)
		if (id < new_constant_offset)
			fprintf(stderr, "theory.free_concept_id WARNING: The given `id` is less than `new_constant_offset`.\n");
#endif
		core::free(ground_concepts[id - new_constant_offset]);
		ground_concepts[id - new_constant_offset].types.keys = NULL;
	}

	void try_free_concept_id(unsigned int id) {
		if (id < new_constant_offset) return;
		const concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.size != 0 || c.negated_types.size != 0 || c.relations.size != 0 || c.negated_relations.size != 0)
			return;

		bool contains;
		unsigned int count = sets.symbols_in_formulas.counts.get(id, contains);
		if (contains && count > 0) return;

		if (sets.element_map.table.contains(id)) return;

		free_concept_id(id);
	}

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != NULL
			 && !ground_concepts[i].print_axioms(out, std::forward<Printer>(printer)...)) return false;
		}
		return sets.print_axioms(out, std::forward<Printer>(printer)...);
	}

	bool add_formula(Formula* formula, unsigned int& new_constant)
	{
		new_constant = 0;

		Formula* new_formula = preprocess_formula(formula);
		if (new_formula == NULL) return false;

		Formula* canonicalized = Canonicalizer::canonicalize(*new_formula);
		free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
		if (canonicalized == NULL) return false;

/* TODO: for debugging; delete this */
print("canonicalized: ", stderr); print(*canonicalized, stderr); print('\n', stderr);
		Proof* new_proof = make_proof<false, true, true>(canonicalized, new_constant);
if (new_proof != NULL) {
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put((unsigned int) built_in_predicates::UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof<built_in_predicates, typename ProofCalculus::ProofCanonicalizer>(*new_proof, expected_conclusion))
fprintf(stderr, "add_formula WARNING: `check_proof` failed.\n");
free(*expected_conclusion); if (expected_conclusion->reference_count == 0) free(expected_conclusion);
}
		core::free(*canonicalized);
		if (canonicalized->reference_count == 0)
			core::free(canonicalized);
		if (new_proof == NULL) {
			return false;
		} else if (!observations.add(new_proof)) {
			free_proof(new_proof);
			return false;
		}

		/* add the axioms in the new proof to `proof_axioms` */
		if (!get_axioms(new_proof, proof_axioms)) {
			observations.remove(new_proof);
			free_proof(new_proof);
			return false;
		}
		return true;
	}

	bool check_proof_axioms() const
	{
		bool success = true;
		hash_multiset<Formula*, false> computed_axioms(64);
		for (const Proof* proof : observations) {
			array_multiset<Formula*, false> axioms(16);
			if (!get_axioms(proof, axioms)) {
				fprintf(stderr, "theory.check_proof_axioms ERROR: get_axioms failed.\n");
				return false;
			}

			for (const auto& entry : axioms.counts) {
				bool contains;
				unsigned int count = proof_axioms.counts.get(entry.key, contains);
				if (!contains || count == 0) {
					print("theory.check_proof_axioms WARNING: Found axiom of a proof in"
							" `observations` that is not in `proof_axioms`.\n", stderr);
					print("  Axiom: ", stderr); print(*entry.key, stderr); print('\n', stderr);
					success = false;
				}
			}

			if (!computed_axioms.add(axioms)) return false;
		}

		for (const auto& entry : proof_axioms.counts) {
			bool contains;
			unsigned int count = computed_axioms.counts.get(entry.key, contains);
			if (!contains) {
				print("theory.check_proof_axioms WARNING: `proof_axioms` contains an "
						"axiom that does not belong to a proof in `observations`.\n", stderr);
				print("  Axiom: ", stderr); print(*entry.key, stderr); print('\n', stderr);
				success = false;
			} else if (entry.value != count) {
				print("theory.check_proof_axioms WARNING: The axiom '", stderr);
				print(*entry.key, stderr); print("' has expected frequency ", stderr);
				print(count, stderr); print(" but has computed frequency ", stderr);
				print(entry.value, stderr); print('\n', stderr);
				success = false;
			}
		}

		return success;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_unary_atom(unsigned int predicate, unsigned int arg, Proof* axiom, Args&&... visitor)
	{
		Formula* lifted_literal;
		Formula* lifted_atom = Formula::new_atom(predicate, Formula::new_variable(1));
		if (lifted_atom == NULL) return false;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom); return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}
		if (!move_element_to_subset<ResolveInconsistencies>(arg, lifted_literal, std::forward<Args>(visitor)...)) {
			free(*lifted_literal); free(lifted_literal);
			return false;
		}
		free(*lifted_literal); free(lifted_literal);

#if !defined(NDEBUG)
		if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_unary_atom WARNING: `ground_concepts` does not contain the concept %u.\n", arg);
#endif
		bool contains; unsigned int bucket;
		if (!types.check_size()) return false;
		pair<array<unsigned int>, array<unsigned int>>& instance_pair = types.get(predicate, contains, bucket);
		if (!contains) {
			if (!array_init(instance_pair.key, 8)) {
				return false;
			} else if (!array_init(instance_pair.value, 8)) {
				core::free(instance_pair.key); return false;
			}
			types.table.keys[bucket] = predicate;
			types.table.size++;
		}

		array<unsigned int>& instances = (Negated ? instance_pair.value : instance_pair.key);
		array_map<unsigned int, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);
		if (!instances.ensure_capacity(instances.length + 1)
		 || !ground_types.ensure_capacity(ground_types.size + 1)) return false;

		add_sorted<false>(instances, arg);
		ground_types.keys[ground_types.size] = predicate;
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		ground_axiom_count++;
		return true;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_binary_atom(relation rel, Proof* axiom, Args&&... visitor)
	{
		/* check if this observation is inconsistent with the theory (can we prove its negation?) */
		if (rel.arg1 == rel.arg2) {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_and(
					Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_constant(rel.arg2)),
					Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_variable(1)),
					Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), Formula::new_variable(1)));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_subset<ResolveInconsistencies>(rel.arg1, lifted_literal, std::forward<Args>(visitor)...)) {
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			free(*lifted_literal); free(lifted_literal);

		} else {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_constant(rel.arg2));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_subset<ResolveInconsistencies>(rel.arg1, lifted_literal, std::forward<Args>(visitor)...)) {
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			free(*lifted_literal); free(lifted_literal);

			lifted_atom = Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), Formula::new_variable(1));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_subset<ResolveInconsistencies>(rel.arg2, lifted_literal, std::forward<Args>(visitor)...)) {
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			free(*lifted_literal); free(lifted_literal);
		}

#if !defined(NDEBUG)
		if (rel.arg1 < new_constant_offset || ground_concepts[rel.arg1 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (rel.arg2 < new_constant_offset || ground_concepts[rel.arg2 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg2);
#endif

		bool contains; unsigned int bucket;
		if (!relations.check_size(relations.table.size + 4)) return false;
		pair<array<unsigned int>, array<unsigned int>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(predicate_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(predicate_instance_pair.value, 8)) {
				core::free(predicate_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {0, rel.arg1, rel.arg2};
			relations.table.size++;
		}

		pair<array<unsigned int>, array<unsigned int>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(arg1_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(arg1_instance_pair.value, 8)) {
				core::free(arg1_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, 0, rel.arg2};
			relations.table.size++;
		}

		pair<array<unsigned int>, array<unsigned int>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0}, contains, bucket);
		if (!contains) {
			if (!array_init(arg2_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(arg2_instance_pair.value, 8)) {
				core::free(arg2_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, rel.arg1, 0};
			relations.table.size++;
		}

		array<unsigned int>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<unsigned int>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<unsigned int>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);
		if (!predicate_instances.ensure_capacity(predicate_instances.length + 1)
		 || !arg1_instances.ensure_capacity(arg1_instances.length + 1)
		 || !arg2_instances.ensure_capacity(arg2_instances.length + 1)
		 || !ground_arg1.ensure_capacity(ground_arg1.size + 3)
		 || !ground_arg2.ensure_capacity(ground_arg2.size + 3)) return false;

		if (rel.arg1 == rel.arg2) {
			pair<array<unsigned int>, array<unsigned int>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0}, contains, bucket);
			if (!contains) {
				if (!array_init(both_arg_instance_pair.key, 8)) {
					return false;
				} else if (!array_init(both_arg_instance_pair.value, 8)) {
					core::free(both_arg_instance_pair.key);
					return false;
				}
				relations.table.keys[bucket] = {rel.predicate, 0, 0};
				relations.table.size++;
			}

			array<unsigned int>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			if (!both_arg_instances.ensure_capacity(both_arg_instances.length + 1))
				return false;
			both_arg_instances[both_arg_instances.length++] = rel.arg1;
			insertion_sort(both_arg_instances);
		}

		add_sorted<false>(predicate_instances, rel.predicate);
		add_sorted<false>(arg1_instances, rel.arg1);
		add_sorted<false>(arg2_instances, rel.arg2);
		ground_arg1.keys[ground_arg1.size] = {rel.predicate, 0, rel.arg2};
		ground_arg1.values[ground_arg1.size++] = axiom;
		ground_arg2.keys[ground_arg2.size] = {rel.predicate, rel.arg1, 0};
		ground_arg2.values[ground_arg2.size++] = axiom;
		axiom->reference_count += 2;
		ground_axiom_count += 2;
		if (rel.arg1 == rel.arg2) {
			/* in this case, `ground_arg1` and `ground_arg2` are the same */
			ground_arg1.keys[ground_arg1.size] = {rel.predicate, 0, 0};
			ground_arg1.values[ground_arg1.size++] = axiom;
			axiom->reference_count++;
			ground_axiom_count++;
		}
		return true;
	}

	template<bool Negated, typename... Args>
	inline bool remove_unary_atom(unsigned int predicate, unsigned int arg, Args&&... visitor)
	{
		Formula* lifted_literal;
		Formula* lifted_atom = Formula::new_atom(predicate, Formula::new_variable(1));
		if (lifted_atom == NULL) return false;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom); return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}
		if (!move_element_to_superset(arg, lifted_literal, std::forward<Args>(visitor)...)) {
			fprintf(stderr, "theory.remove_unary_atom ERROR: `move_element_to_superset` failed.\n");
			if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
			return false;
		}
		free(*lifted_literal); free(lifted_literal);

#if !defined(NDEBUG)
		if (!types.table.contains(predicate))
			fprintf(stderr, "theory.remove_unary_atom WARNING: `types` does not contain the key %u.\n", predicate);
		if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", arg);
#endif

		array<unsigned int>& instances = (Negated ? types.get(predicate).value : types.get(predicate).key);
		array_map<unsigned int, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);

		unsigned int index = instances.index_of(arg);
#if !defined(NDEBUG)
		if (index == instances.length)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `instances` does not contain %u.\n", arg);
#endif
		shift_left(instances.data + index, instances.length - index - 1);
		instances.length--;

		index = ground_types.index_of(predicate);
#if !defined(NDEBUG)
		if (index == ground_types.size)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_types` does not contain %u.\n", predicate);
#endif
		Proof* axiom = ground_types.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		ground_types.remove_at(index);
		ground_axiom_count--;
		return true;
	}

	template<bool Negated, typename... Args>
	inline bool remove_binary_atom(relation rel, Args&&... visitor)
	{
		if (rel.arg1 == rel.arg2) {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_and(
					Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), Formula::new_variable(1)),
					Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_variable(1)),
					Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_constant(rel.arg2)));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_superset(rel.arg1, lifted_literal, std::forward<Args>(visitor)...)) {
				fprintf(stderr, "theory.remove_binary_atom ERROR: `move_element_to_superset` failed.\n");
				if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
				return false;
			}
			free(*lifted_literal); free(lifted_literal);

		} else {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_atom(rel.predicate, Formula::new_variable(1), Formula::new_constant(rel.arg2));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_superset(rel.arg1, lifted_literal, std::forward<Args>(visitor)...)) {
				fprintf(stderr, "theory.remove_binary_atom ERROR: `move_element_to_superset` failed.\n");
				if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
				return false;
			}
			free(*lifted_literal); free(lifted_literal);

			lifted_atom = Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), Formula::new_variable(1));
			if (lifted_atom == NULL) return false;
			if (Negated) {
				lifted_literal = Formula::new_not(lifted_atom);
				if (lifted_literal == NULL) {
					free(*lifted_atom); free(lifted_atom); return false;
				}
			} else {
				lifted_literal = lifted_atom;
			}
			if (!move_element_to_superset(rel.arg2, lifted_literal, std::forward<Args>(visitor)...)) {
				fprintf(stderr, "theory.remove_binary_atom ERROR: `move_element_to_superset` failed.\n");
				if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
				return false;
			}
			free(*lifted_literal); free(lifted_literal);
		}

#if !defined(NDEBUG)
		if (!relations.table.contains({0, rel.arg1, rel.arg2})
		 || !relations.table.contains({rel.predicate, 0, rel.arg2})
		 || !relations.table.contains({rel.predicate, rel.arg1, 0}))
			fprintf(stderr, "theory.add_binary_atom WARNING: `relations` does not contain the necessary relations.\n");
		if (rel.arg1 < new_constant_offset || ground_concepts[rel.arg1 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (rel.arg2 < new_constant_offset || ground_concepts[rel.arg2 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg2);
#endif

		pair<array<unsigned int>, array<unsigned int>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0});

		array<unsigned int>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<unsigned int>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<unsigned int>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);

		unsigned int index;
		if (rel.arg1 == rel.arg2) {
			pair<array<unsigned int>, array<unsigned int>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0});
			array<unsigned int>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			index = both_arg_instances.index_of(rel.arg1);
#if !defined(NDEBUG)
			if (index == both_arg_instances.length)
				fprintf(stderr, "theory.add_binary_atom WARNING: `both_arg_instances` does not contain %u.\n", rel.arg1);
#endif
			shift_left(both_arg_instances.data + index, both_arg_instances.length - index - 1);
			both_arg_instances.length--;
		}

		index = predicate_instances.index_of(rel.predicate);
#if !defined(NDEBUG)
		if (index == predicate_instances.length)
			fprintf(stderr, "theory.add_binary_atom WARNING: `predicate_instances` does not contain %u.\n", rel.predicate);
#endif
		shift_left(predicate_instances.data + index, predicate_instances.length - index - 1);
		predicate_instances.length--;

		index = arg1_instances.index_of(rel.arg1);
#if !defined(NDEBUG)
		if (index == arg1_instances.length)
			fprintf(stderr, "theory.add_binary_atom WARNING: `arg1_instances` does not contain %u.\n", rel.arg1);
#endif
		shift_left(arg1_instances.data + index, arg1_instances.length - index - 1);
		arg1_instances.length--;

		index = arg2_instances.index_of(rel.arg2);
#if !defined(NDEBUG)
		if (index == arg2_instances.length)
			fprintf(stderr, "theory.add_binary_atom WARNING: `arg2_instances` does not contain %u.\n", rel.arg2);
#endif
		shift_left(arg2_instances.data + index, arg2_instances.length - index - 1);
		arg2_instances.length--;

		index = ground_arg1.index_of(relation(rel.predicate, 0, rel.arg2));
#if !defined(NDEBUG)
		if (index == ground_arg1.size)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_arg1` does not contain the requested predicate.\n");
#endif
		Proof* axiom = ground_arg1.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		ground_arg1.remove_at(index);

		index = ground_arg2.index_of(relation(rel.predicate, rel.arg1, 0));
#if !defined(NDEBUG)
		if (index == ground_arg2.size)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_arg2` does not contain the requested predicate.\n");
#endif
		axiom = ground_arg2.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		ground_arg2.remove_at(index);
		ground_axiom_count -= 2;

		if (rel.arg1 == rel.arg2) {
			/* in this case, `ground_arg1` and `ground_arg2` are the same */
			index = ground_arg1.index_of(relation(rel.predicate, 0, 0));
#if !defined(NDEBUG)
			if (index == ground_arg1.size)
				fprintf(stderr, "theory.add_binary_atom WARNING: `ground_arg1` does not contain the requested predicate.\n");
#endif
			axiom = ground_arg1.values[index];
			core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
			ground_arg1.remove_at(index);
			ground_axiom_count--;
		}

		return true;
	}

	unsigned int index_of(const array<pair<Formula*, Proof*>>& elements, Proof* proof) const {
		unsigned int index = elements.length;
		for (unsigned int i = 0; i < elements.length; i++)
			if (elements[i].value == proof) { index = i; break; }
#if !defined(NDEBUG)
		if (index == elements.length)
			fprintf(stderr, "theory.index_of WARNING: `elements` does not contain `proof`.\n");
#endif
		return index;
	}

	struct changes {
		array<pair<pair<unsigned int, unsigned int>, Proof*>> unary_atoms;
		array<pair<pair<unsigned int, unsigned int>, Proof*>> negated_unary_atoms;
		array<pair<relation, Proof*>> binary_atoms;
		array<pair<relation, Proof*>> negated_binary_atoms;

		array<Proof*> subset_axioms;
		array<Proof*> set_size_axioms;

		array<pair<Formula*, Proof*>> implication_intro_nodes;
		array<pair<Formula*, Proof*>> negated_conjunction_nodes;
		array<pair<Formula*, Proof*>> existential_intro_nodes;
		array<pair<Formula*, Proof*>> disjunction_intro_nodes;

		changes() : unary_atoms(4), negated_unary_atoms(4),
				binary_atoms(4), negated_binary_atoms(4),
				subset_axioms(4), set_size_axioms(4),
				implication_intro_nodes(4), negated_conjunction_nodes(4),
				existential_intro_nodes(4), disjunction_intro_nodes(4)
		{ }
	};

	template<typename... Args>
	bool add_changes(const theory::changes& changes, Args&&... visitor)
	{
		for (const auto& entry : changes.unary_atoms) {
			if (!try_init_concept(entry.key.key) || !try_init_concept(entry.key.value)
			 || !add_unary_atom<false, false>(entry.key.key, entry.key.value, entry.value, std::forward<Args>(visitor)...)) return false;
		} for (const auto& entry : changes.negated_unary_atoms) {
			if (!try_init_concept(entry.key.key) || !try_init_concept(entry.key.value)
			 || !add_unary_atom<true, false>(entry.key.key, entry.key.value, entry.value, std::forward<Args>(visitor)...)) return false;
		} for (const auto& entry : changes.binary_atoms) {
			if (!try_init_concept(entry.key.predicate)
			 || !try_init_concept(entry.key.arg1) || !try_init_concept(entry.key.arg2)
			 || !add_binary_atom<false, false>(entry.key, entry.value, std::forward<Args>(visitor)...)) return false;
		} for (const auto& entry : changes.negated_binary_atoms) {
			if (!try_init_concept(entry.key.predicate)
			 || !try_init_concept(entry.key.arg1) || !try_init_concept(entry.key.arg2)
			 || !add_binary_atom<true, false>(entry.key, entry.value, std::forward<Args>(visitor)...)) return false;
		} for (Proof* subset_axiom : changes.subset_axioms) {
			unsigned int variable = subset_axiom->formula->quantifier.variable;
			if (subset_axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
				unsigned int predicate; Term const* arg1; Term const* arg2;
				Formula* left = subset_axiom->formula->quantifier.operand->binary.left;
				if (is_atomic(*left, predicate, arg1, arg2)
				 && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL)
				{
					if (!try_init_concept(predicate)) return false;
				}
			}
			if (!sets.add_subset_axiom(subset_axiom)) return false;
		} for (Proof* axiom : changes.set_size_axioms) {
			unsigned int set_id;
			Formula* set_formula = axiom->formula->binary.left->binary.right->quantifier.operand;
			if (!sets.get_set_id(set_formula, set_id)
			 || !sets.sets[set_id].set_size_axiom(axiom)) return false;
		} for (const auto& entry : changes.implication_intro_nodes) {
			if (!implication_intro_nodes.add(entry)) return false;
			entry.key->reference_count++;
		} for (const auto& entry : changes.negated_conjunction_nodes) {
			if (!negated_conjunction_nodes.add(entry)) return false;
			entry.key->reference_count++;
		} for (const auto& entry : changes.existential_intro_nodes) {
			if (!existential_intro_nodes.add(entry)) return false;
			entry.key->reference_count++;
		} for (const auto& entry : changes.disjunction_intro_nodes) {
			if (!disjunction_intro_nodes.add(entry)) return false;
			entry.key->reference_count++;
		}
		return true;
	}

	template<typename... Args>
	bool subtract_changes(const theory::changes& changes, Args&&... visitor)
	{
		for (const auto& entry : changes.unary_atoms) {
			remove_unary_atom<false>(entry.key.key, entry.key.value, std::forward<Args>(visitor)...);
			try_free_concept_id(entry.key.key);
			try_free_concept_id(entry.key.value);
		} for (const auto& entry : changes.negated_unary_atoms) {
			remove_unary_atom<true>(entry.key.key, entry.key.value, std::forward<Args>(visitor)...);
			try_free_concept_id(entry.key.key);
			try_free_concept_id(entry.key.value);
		} for (const auto& entry : changes.binary_atoms) {
			remove_binary_atom<false>(entry.key, std::forward<Args>(visitor)...);
			try_free_concept_id(entry.key.predicate);
			try_free_concept_id(entry.key.arg1);
			try_free_concept_id(entry.key.arg2);
		} for (const auto& entry : changes.negated_binary_atoms) {
			remove_binary_atom<true>(entry.key, std::forward<Args>(visitor)...);
			try_free_concept_id(entry.key.predicate);
			try_free_concept_id(entry.key.arg1);
			try_free_concept_id(entry.key.arg2);
		} for (Proof* subset_axiom : changes.subset_axioms) {
			unsigned int variable = subset_axiom->formula->quantifier.variable;
			if (subset_axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
				unsigned int predicate; Term const* arg1; Term const* arg2;
				Formula* left = subset_axiom->formula->quantifier.operand->binary.left;
				if (is_atomic(*left, predicate, arg1, arg2)
				 && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL)
				{
					try_free_concept_id(predicate);
				}
			}
			sets.free_subset_axiom(subset_axiom);
		} for (Proof* axiom : changes.set_size_axioms) {
			unsigned int set_id;
			Formula* set_formula = axiom->formula->binary.left->binary.right->quantifier.operand;
			if (!sets.get_set_id(set_formula, set_id)) return false;
			axiom->reference_count--;
			bool is_freeable = sets.is_freeable(set_id);
			axiom->reference_count++;
			if (is_freeable)
				sets.free_set(set_formula, set_id);
		} for (const auto& entry : changes.implication_intro_nodes) {
			implication_intro_nodes.remove(index_of(implication_intro_nodes, entry.value));
		} for (const auto& entry : changes.negated_conjunction_nodes) {
			negated_conjunction_nodes.remove(index_of(negated_conjunction_nodes, entry.value));
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		} for (const auto& entry : changes.existential_intro_nodes) {
			existential_intro_nodes.remove(index_of(existential_intro_nodes, entry.value));
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		} for (const auto& entry : changes.disjunction_intro_nodes) {
			disjunction_intro_nodes.remove(index_of(disjunction_intro_nodes, entry.value));
			free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
		}
		return true;
	}

	template<typename... Visitor>
	bool get_theory_changes(
			Proof& proof, array<const Proof*>& discharged_axioms,
			array_map<const Proof*, unsigned int>& reference_counts,
			theory::changes& changes, Visitor&&... visitor) const
	{
		visit_node(proof, std::forward<Visitor>(visitor)...);
#if !defined(NDEBUG)
		bool contains;
		unsigned int reference_count = reference_counts.get(&proof, contains);
		if (!contains) fprintf(stderr, "theory.remove_proof WARNING: The given proof is not in the map `reference_counts`.\n");
#else
		unsigned int reference_count = reference_counts.get(&proof);
#endif

		unsigned int old_discharged_axiom_count = discharged_axioms.length;

		unsigned int predicate = 0; Term* arg1 = nullptr; Term* arg2 = nullptr;
		switch (proof.type) {
		case ProofType::AXIOM:
			if (discharged_axioms.contains(&proof)) break;
			if (is_atomic(*proof.formula, predicate, arg1, arg2)) {
				if (reference_count != 1) return true;
				if (arg1->type == TermType::CONSTANT && arg2 == NULL) {
					if (!visit_unary_atom<false>(predicate, arg1->constant, std::forward<Visitor>(visitor)...)
					 || !changes.unary_atoms.add({{predicate, arg1->constant}, &proof})) return false;
				} else if (arg1->type == TermType::CONSTANT && arg2->type == TermType::CONSTANT) {
					if (!visit_binary_atom<false>(predicate, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
					 || !changes.binary_atoms.add({{predicate, arg1->constant, arg2->constant}, &proof})) return false;
				} else {
					fprintf(stderr, "free WARNING: Found unexpected literal axiom.\n");
				}
			} else if (proof.formula->type == FormulaType::NOT) {
				if (is_atomic(*proof.formula->unary.operand, predicate, arg1, arg2)) {
					if (reference_count != 1) return true;
					if (arg1->type == TermType::CONSTANT && arg2 == NULL) {
						if (!visit_unary_atom<true>(predicate, arg1->constant, std::forward<Visitor>(visitor)...)
						 || !changes.negated_unary_atoms.add({{predicate, arg1->constant}, &proof})) return false;
					} else if (arg1->type == TermType::CONSTANT && arg2->type == TermType::CONSTANT) {
						if (!visit_binary_atom<true>(predicate, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
						 || !changes.negated_binary_atoms.add({{predicate, arg1->constant, arg2->constant}, &proof})) return false;
					} else {
						fprintf(stderr, "free WARNING: Found unexpected literal axiom.\n");
					}
				} else {
					fprintf(stderr, "free WARNING: Found unexpected axiom.\n");
				}
			} else if (proof.formula->type == FormulaType::FOR_ALL
					&& proof.formula->quantifier.operand->type == FormulaType::IF_THEN)
			{
				if (reference_count == 1)
					return visit_subset_axiom(proof, std::forward<Visitor>(visitor)...)
						&& changes.subset_axioms.add(&proof);
				return true;
			} else if (proof.formula->type == FormulaType::EQUALS) {
				if (is_atomic(*proof.formula->binary.left, predicate, arg1, arg2)
				 && predicate == (unsigned int) built_in_predicates::SIZE
				 && arg1->type == TermType::LAMBDA && arg2 == NULL
				 && proof.formula->binary.right->type == TermType::INTEGER)
				{
					if (reference_count == 1)
						return changes.set_size_axioms.add(&proof);
					return true;
				} else {
					fprintf(stderr, "free WARNING: Found unexpected equals axiom.\n");
				}
			} else {
				fprintf(stderr, "free WARNING: Found unexpected axiom.\n");
			}
			break;

		case ProofType::IMPLICATION_INTRODUCTION:
			if (reference_count != 0) return true;
			if (!changes.implication_intro_nodes.add({implication_intro_nodes[index_of(implication_intro_nodes, &proof)].key, &proof})
			 || !discharged_axioms.add(proof.operands[1]))
				return false;
			break;

		case ProofType::PROOF_BY_CONTRADICTION:
			if (reference_count != 0) return true;
			discharged_axioms.add(proof.operands[1]);
			if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
			  && proof.operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
			  && proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				if (!visit_negated_conjunction(std::forward<Visitor>(visitor)...)
				 || !changes.negated_conjunction_nodes.add(
					{negated_conjunction_nodes[index_of(negated_conjunction_nodes, &proof)].key, &proof})) return false;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::UNIVERSAL_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				if (!visit_negated_universal_intro(std::forward<Visitor>(visitor)...)
				 || !changes.existential_intro_nodes.add(
					{existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key, &proof})) return false;
			} else if (proof.operands[0]->type == ProofType::DISJUNCTION_ELIMINATION
					&& proof.operands[0]->operands[0] == proof.operands[1])
			{
				break;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::IMPLICATION_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				break;
			} else {
				fprintf(stderr, "free WARNING: Found unexpected proof by contradiction step.\n");
			}
			break;

		case ProofType::EXISTENTIAL_INTRODUCTION:
			if (reference_count != 0) return true;
			if (!visit_existential_intro(std::forward<Visitor>(visitor)...)
			 || !changes.existential_intro_nodes.add(
				{existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key, &proof})) return false;
			break;

		case ProofType::DISJUNCTION_INTRODUCTION:
		case ProofType::DISJUNCTION_INTRODUCTION_LEFT:
		case ProofType::DISJUNCTION_INTRODUCTION_RIGHT:
			if (reference_count != 0) return true;
			if (!visit_disjunction_intro(std::forward<Visitor>(visitor)...)
			 || !changes.disjunction_intro_nodes.add(
				{disjunction_intro_nodes[index_of(disjunction_intro_nodes, &proof)].key, &proof})) return false;
			break;

		case ProofType::CONJUNCTION_INTRODUCTION:
		case ProofType::BETA_EQUIVALENCE:
		case ProofType::CONJUNCTION_ELIMINATION:
		case ProofType::CONJUNCTION_ELIMINATION_LEFT:
		case ProofType::CONJUNCTION_ELIMINATION_RIGHT:
		case ProofType::DISJUNCTION_ELIMINATION:
		case ProofType::IMPLICATION_ELIMINATION:
		case ProofType::BICONDITIONAL_INTRODUCTION:
		case ProofType::BICONDITIONAL_ELIMINATION_LEFT:
		case ProofType::BICONDITIONAL_ELIMINATION_RIGHT:
		case ProofType::NEGATION_ELIMINATION:
		case ProofType::FALSITY_ELIMINATION:
		case ProofType::UNIVERSAL_INTRODUCTION:
		case ProofType::UNIVERSAL_ELIMINATION:
		case ProofType::EXISTENTIAL_ELIMINATION:
		case ProofType::EQUALITY_ELIMINATION:
		case ProofType::PARAMETER:
		case ProofType::TERM_PARAMETER:
		case ProofType::ARRAY_PARAMETER:
		case ProofType::FORMULA_PARAMETER:
		case ProofType::COUNT:
			if (reference_count != 0) return true;
			break;
		}

		unsigned int operand_count;
		Proof* const* operands;
		proof.get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			if (!reference_counts.ensure_capacity(reference_counts.size + 1))
				return false;
			unsigned int index = reference_counts.index_of(operands[i]);
			if (index == reference_counts.size) {
				reference_counts.keys[index] = operands[i];
				reference_counts.values[index] = operands[i]->reference_count;
				reference_counts.size++;
			}
			reference_counts.values[index]--;
			if (!get_theory_changes(*operands[i], discharged_axioms, reference_counts, changes, std::forward<Visitor>(visitor)...))
				return false;
		}
		discharged_axioms.length = old_discharged_axiom_count;
		return true;
	}

	template<typename... Visitor>
	inline bool get_theory_changes(Proof& proof, theory::changes& changes, Visitor&&... visitor) const
	{
		array<const Proof*> discharged_axioms(16);
		array_map<const Proof*, unsigned int> reference_counts(32);
		reference_counts.keys[0] = &proof;
		reference_counts.values[0] = proof.reference_count - 1;
		reference_counts.size++;
		return get_theory_changes(proof, discharged_axioms, reference_counts, changes, std::forward<Visitor>(visitor)...);
	}

	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_proof(Formula* canonicalized, unsigned int& new_constant, Args&&... args)
	{
		unsigned int predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<DefinitionsAllowed, Contradiction, ResolveInconsistencies>(predicate, arg1, arg2, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::NOT) {
			return make_proof<!Contradiction, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->unary.operand, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
					return NULL;

				Term* constant;
				Proof* exists_not_proof = make_exists_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (exists_not_proof == NULL) return NULL;
				Proof* assumption = ProofCalculus::new_axiom(canonicalized);
				if (assumption == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof); return NULL;
				}
				assumption->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_universal_elim(assumption, constant), exists_not_proof), assumption);
				free(*assumption); if (assumption->reference_count == 0) free(assumption);
				if (proof == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof); return NULL;
				}
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				proof->reference_count++;
				/* record the formula with this existential */
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				return proof;
			}

			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == FormulaType::IF_THEN) {
				Formula* left = canonicalized->quantifier.operand->binary.left;
				Formula* right = canonicalized->quantifier.operand->binary.right;
				Term const* arg1; Term const* arg2;
				bool atomic = is_atomic(*left, predicate, arg1, arg2);
				if (atomic && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL
				 && predicate == (unsigned int) built_in_predicates::UNKNOWN)
				{
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_proof ERROR: Definitions are not allowed in this context.\n");
						return NULL;
					}

					/* this is a definition of a type */
					/* check the right-hand side is a valid definition */
					if (!valid_definition(right, variable)) {
						fprintf(stderr, "theory.make_proof ERROR: This is not a valid type definition.\n");
						return NULL;
					}

					new_constant = get_free_concept_id();

					if (!try_init_concept(new_constant)) return NULL;

					/* TODO: definitions of new concepts should be biconditionals */
					Formula* new_left = Formula::new_atom(new_constant, Formula::new_variable(variable));
					if (new_left == NULL) { free_concept_id(new_constant); return NULL; }

					Proof* new_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(new_left, right);
					new_axiom->reference_count++;
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					return new_axiom;
				} else {
					/* make sure there are no unknown predicates */
					if (contains_constant(*left, (unsigned int) built_in_predicates::UNKNOWN)
					 || contains_constant(*right, (unsigned int) built_in_predicates::UNKNOWN)) {
						fprintf(stderr, "theory.make_proof ERROR: Universally-quantified statement has unknown constants.\n");
						return NULL;
					}

					/* this is a formula of form `![x]:(t(x) => f(x))` */
					Proof* new_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(left, right);
					new_axiom->reference_count++;
					return new_axiom;
				}
			}

		} else if (canonicalized->type == FormulaType::AND) {
			if (Contradiction) {
				if (!negated_conjunction_nodes.ensure_capacity(negated_conjunction_nodes.length + 1))
					return NULL;

				array<unsigned int> indices(canonicalized->array.length);
				for (unsigned int i = 0; i < canonicalized->array.length; i++)
					indices[i] = i + 1;
				indices.length = canonicalized->array.length;
				shuffle(indices);

				unsigned int old_index_count = indices.length;
				if (!filter_constants(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

				Formula* negated_conjunction = Formula::new_not(canonicalized);
				if (negated_conjunction == NULL) return NULL;
				canonicalized->reference_count++;

				for (unsigned int index : indices) {
					unsigned int index_minus_one = index - 1;
					Proof* operand = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index_minus_one], new_constant, std::forward<Args>(args)...);
					if (operand != NULL) {
						/* we found a disproof of a conjunct */
						Proof* axiom = ProofCalculus::new_axiom(canonicalized);
						if (axiom == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
							free_proof(operand); return NULL;
						}
						axiom->reference_count++;
						Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
								ProofCalculus::new_conjunction_elim(axiom, make_array_view(&index_minus_one, 1)), operand), axiom);
						free(*axiom); if (axiom->reference_count == 0) free(axiom);
						if (proof == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
							free_proof(operand); return NULL;
						}
						free(*operand); if (operand->reference_count == 0) free(operand);
						proof->reference_count++;
						/* record the formula with this negated conjunction node */
						negated_conjunction_nodes[negated_conjunction_nodes.length++] = { negated_conjunction, proof };
						return proof;
					}

					if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
				}

				/* we couldn't find a disproof of any of the conjuncts */
				finished_constants(canonicalized, old_index_count, std::forward<Args>(args)...);
				free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
				return NULL;
			}

			Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
			if (operands == NULL) {
				fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
				return NULL;
			}
			for (unsigned int i = 0; i < canonicalized->array.length; i++) {
				if (new_constant == 0)
					operands[i] = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);
				else operands[i] = make_proof<false, false, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);

				if (operands[i] == NULL) {
					for (unsigned int j = 0; j < i; j++)
						/* undo the changes made by the recursive calls to `make_proof` */
						free_proof(operands[j]);
					free(operands); return NULL;
				}
			}
			Proof* conjunction = ProofCalculus::new_conjunction_intro(make_array_view(operands, canonicalized->array.length));
			if (conjunction == NULL) {
				for (unsigned int j = 0; j < canonicalized->array.length; j++)
					/* undo the changes made by the recursive calls to `make_proof` */
					free_proof(operands[j]);
				free(operands); return NULL;
			}
			for (unsigned int j = 0; j < canonicalized->array.length; j++) {
				free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
			}
			free(operands);
			conjunction->reference_count++;
			return conjunction;

		} else if (canonicalized->type == FormulaType::OR) {
			if (Contradiction) {
				Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
				if (operands == NULL) {
					fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
					return NULL;
				}
				for (unsigned int i = 0; i < canonicalized->array.length; i++) {
					if (new_constant == 0)
						operands[i] = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);
					else operands[i] = make_proof<true, false, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);

					if (operands[i] == NULL) {
						for (unsigned int j = 0; j < i; j++)
							/* undo the changes made by recursive calls to `make_proof` */
							free_proof(operands[j]);
						free(operands); return NULL;
					} else if (operands[i]->type == ProofType::PROOF_BY_CONTRADICTION) {
						Proof* absurdity = operands[i]->operands[0];
						absurdity->reference_count++;
						free(*operands[i]); if (operands[i]->reference_count == 0) free(operands[i]);
						operands[i] = absurdity;
					} else {
						Proof* absurdity = ProofCalculus::new_negation_elim(
								ProofCalculus::new_axiom(canonicalized->array.operands[i]), operands[i]);
						if (absurdity == NULL) {
							for (unsigned int j = 0; j < i; j++)
								/* undo the changes made by recursive calls to `make_proof` */
								free_proof(operands[j]);
							free(operands); return NULL;
						}
						absurdity->reference_count++;
						free(*operands[i]); if (operands[i]->reference_count == 0) free(operands[i]);
						operands[i] = absurdity;
					}
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					for (unsigned int j = 0; j < canonicalized->array.length; j++)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j]);
					free(operands); return NULL;
				}
				axiom->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(
						ProofCalculus::new_disjunction_elim(axiom, make_array_view(operands, canonicalized->array.length)), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					for (unsigned int j = 0; j < canonicalized->array.length; j++)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j]);
					free(operands); return NULL;
				}
				proof->reference_count++;
				for (unsigned int j = 0; j < canonicalized->array.length; j++) {
					free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
				}
				free(operands);
				return proof;
			}

			if (!disjunction_intro_nodes.ensure_capacity(disjunction_intro_nodes.length + 1))
				return NULL;

			array<unsigned int> indices(canonicalized->array.length);
			for (unsigned int i = 0; i < canonicalized->array.length; i++)
				indices[i] = i + 1;
			indices.length = canonicalized->array.length;
			shuffle(indices);

			unsigned int old_index_count = indices.length;
			if (!filter_constants(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

			for (unsigned int index : indices) {
				Proof* operand = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index - 1], new_constant, std::forward<Args>(args)...);
				if (operand != NULL) {
					/* we found a proof of a disjunct */
					array<Formula*> other_disjuncts(max(1u, canonicalized->array.length - 1));
					for (unsigned int i = 0; i < canonicalized->array.length; i++) {
						if (i == index - 1) continue;
						other_disjuncts[other_disjuncts.length++] = canonicalized->array.operands[i];
					}
					Formula* other_disjunction;
					if (other_disjuncts.length == 1) {
						other_disjunction = other_disjuncts[0];
						other_disjuncts[0]->reference_count++;
					} else {
						other_disjunction = Formula::new_or(make_array_view(other_disjuncts.data, other_disjuncts.length));
						if (other_disjunction == NULL) {
							/* undo changes made by the recursive call to `make_proof` */
							free_proof(operand); return NULL;
						}
						for (Formula* other_disjunct : other_disjuncts)
							other_disjunct->reference_count++;
					}
					Proof* proof = ProofCalculus::new_disjunction_intro(operand, other_disjunction, index - 1);
					free(*other_disjunction); if (other_disjunction->reference_count == 0) free(other_disjunction);
					if (proof == NULL) {
						/* undo changes made by the recursive call to `make_proof` */
						free_proof(operand); return NULL;
					}
					free(*operand); if (operand->reference_count == 0) free(operand);
					proof->reference_count++;
					/* record the formula with this disjunction intro node */
					disjunction_intro_nodes[disjunction_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					return proof;
				}

				if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
			}

			/* we couldn't find a proof of any of the disjuncts */
			finished_constants(canonicalized, old_index_count, std::forward<Args>(args)...);
			return NULL;

		} else if (canonicalized->type == FormulaType::IF_THEN) {
			if (Contradiction) {
				Proof* left = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, new_constant, std::forward<Args>(args)...);
				if (left == NULL) return NULL;

				Proof* right;
				if (new_constant == 0 && DefinitionsAllowed)
					right = make_proof<true, true, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);
				else right = make_proof<true, false, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);

				if (right == NULL) {
					free_proof(left); return NULL;
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					free_proof(left); free_proof(right);
					return NULL;
				}
				axiom->reference_count++;

				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_implication_elim(axiom, left), right), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					free_proof(left); free_proof(right);
					return NULL;
				}
				free(*left); if (left->reference_count == 0) free(left);
				free(*right); if (right->reference_count == 0) free(right);
				proof->reference_count++;
				return proof;
			}

			if (!implication_intro_nodes.ensure_capacity(implication_intro_nodes.length + 1))
				return NULL;

			array<unsigned int> indices(2);
			indices[0] = 1; indices[1] = 2;
			indices.length = 2;
			if (sample_uniform(2) == 2)
				swap(indices[0], indices[1]);
			if (!filter_constants(canonicalized, indices, std::forward<Args>(args)...)) return NULL;
			for (unsigned int i = 0; i < 2; i++) {
				if (indices[i] == 1) {
					Proof* left = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, new_constant, std::forward<Args>(args)...);
					if (left == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* axiom = ProofCalculus::new_axiom(canonicalized->binary.left);
					if (axiom == NULL) {
						free_proof(left); return NULL;
					}
					axiom->reference_count++;

					Proof* proof = ProofCalculus::new_implication_intro(
							ProofCalculus::new_falsity_elim(
								ProofCalculus::new_negation_elim(left, axiom),
								canonicalized->binary.right),
							axiom);
					free(*axiom); if (axiom->reference_count == 0) free(axiom);
					if (proof == NULL) {
						free_proof(left); return NULL;
					}
					free(*left); if (left->reference_count == 0) free(left);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				} else {
					Proof* right = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);
					if (right == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* proof = ProofCalculus::new_implication_intro(right, ProofCalculus::new_axiom(canonicalized->binary.left));
					if (proof == NULL) {
						free_proof(right); return NULL;
					}
					free(*right); if (right->reference_count == 0) free(right);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				}
			}

			finished_constants(canonicalized, 2, std::forward<Args>(args)...);
			return NULL;

		} else if (canonicalized->type == FormulaType::EXISTS) {
			if (Contradiction) {
				/* get empty set size axiom, forcing inconsistencies to be resolved */
				Formula* set_formula = canonicalized->quantifier.operand;
				Formula* lambda_formula = Formula::new_lambda(1, set_formula);
				if (lambda_formula == NULL) return NULL;
				set_formula->reference_count++;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(set_formula, 0);
				if (set_size_axiom == NULL) {
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					return NULL;
				}
				Formula* shifted_lambda_formula = shift_bound_variables(lambda_formula, 1);
				if (shifted_lambda_formula == NULL) {
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					sets.try_free_set(set_formula); return NULL;
				}

				unsigned int zero = 0;
				Formula* beta_left = Formula::new_not(Formula::new_exists(1, Formula::new_apply(shifted_lambda_formula, Formula::new_variable(1))));
				Formula* beta_right = Formula::new_not(Formula::new_exists(1, set_formula));
				if (beta_left == NULL || beta_right == NULL) {
					if (beta_left != NULL) { free(*beta_left); if (beta_left->reference_count == 0) free(beta_left); }
					if (beta_right != NULL) { free(*beta_right); if (beta_right->reference_count == 0) free(beta_right); }
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					sets.try_free_set(set_formula); return NULL;
				}
				set_formula->reference_count++;
				Proof* proof = ProofCalculus::new_equality_elim(
						ProofCalculus::new_beta(beta_left, beta_right),
						ProofCalculus::new_equality_elim(ProofCalculus::new_universal_elim(empty_set_axiom, lambda_formula), set_size_axiom, make_array_view(&zero, 1)),
						make_array_view(&zero, 1));
				free(*beta_left); if (beta_left->reference_count == 0) free(beta_left);
				free(*beta_right); if (beta_right->reference_count == 0) free(beta_right);
				free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
				if (proof == NULL) { sets.try_free_set(set_formula); return NULL; }
				proof->reference_count++;
				return proof;
			} else {
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
					return NULL;

				Term* variable = Formula::new_variable(canonicalized->quantifier.variable);
				if (variable == NULL) return NULL;

				Term* constant;
				Proof* operand = make_exists_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (operand == NULL) {
					free(*variable); if (variable->reference_count == 0) free(variable);
					return NULL;
				}
				free(*constant); if (constant->reference_count == 0) free(constant);

				array<unsigned int> indices(8);
				if (!compute_indices(*canonicalized->quantifier.operand, *variable, indices)) {
					free(*variable); if (variable->reference_count == 0) free(variable);
					free_proof(operand); return NULL;
				}
				free(*variable); if (variable->reference_count == 0) free(variable);
				Proof* proof = ProofCalculus::new_existential_intro(operand, indices.data, indices.length);
				if (proof == NULL) {
					free_proof(operand); return NULL;
				}
				free(*operand); if (operand->reference_count == 0) free(operand);
				proof->reference_count++;
				/* record the formula with this existential */
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				return proof;
			}

		} else if (canonicalized->type == FormulaType::EQUALS) {
			Formula* left = canonicalized->binary.left;
			Formula* right = canonicalized->binary.right;
			Term* arg1; Term* arg2;
			bool atomic = is_atomic(*left, predicate, arg1, arg2);
			if (atomic && predicate == (unsigned int) built_in_predicates::SIZE
			 && arg1->type == TermType::LAMBDA && arg2 == NULL
			 && right->type == TermType::INTEGER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(arg1->quantifier.operand, right->integer);
				set_size_axiom->reference_count++;
				return set_size_axiom;
			} else if (left == right || *left == *right) {
				Proof* proof = ProofCalculus::new_beta(left, left);
				if (proof == NULL) return NULL;
				proof->reference_count++;
				return proof;
			} else if (left->type == hol_term_type::CONSTANT && right->type == hol_term_type::CONSTANT) {
				if (left->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return NULL;
					} else {
						new_constant = right->constant;
						Proof* proof = ProofCalculus::new_beta(right, right);
						if (proof == NULL) return NULL;
						proof->reference_count++;
						return proof;
					}
				} else if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					new_constant = left->constant;
					Proof* proof = ProofCalculus::new_beta(left, left);
					if (proof == NULL) return NULL;
					proof->reference_count++;
					return proof;
				} else {
					if (Contradiction) {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return NULL;
					} else {
						/* this is impossible */
						return NULL;
					}
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

private:
	/* NOTE: this function finds a constant that proves `quantified`, and not ?[x]:`quantified` */
	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_exists_proof(Formula* quantified, unsigned int variable, Term*& constant, unsigned int& new_constant, Args&&... args)
	{
		Formula* var = Formula::new_variable(variable);
		if (var == NULL) return NULL;

		array<unsigned int> constants(ground_concept_capacity + 1);
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			if (ground_concepts[i].types.keys != NULL) constants[constants.length++] = new_constant_offset + i;
		constants[constants.length++] = UINT_MAX;
		shuffle(constants);

		unsigned int original_constant_count = constants.length;
		if (!filter_constants(quantified, constants, std::forward<Args>(args)...)) {
			free(*var); if (var->reference_count == 0) free(var);
			return NULL;
		}

		for (unsigned int id : constants) {
			if (id != UINT_MAX) {
				constant = Formula::new_constant(id);
				if (constant == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}

				Proof* proof = make_proof<Contradiction, DefinitionsAllowed, ResolveInconsistencies>(substituted, new_constant, std::forward<Args>(args)...);
				free(*substituted); if (substituted->reference_count == 0) free(substituted);
				if (proof != NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return proof;
				}
				free(*constant); if (constant->reference_count == 0) free(constant);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
			} else {
				unsigned int constant_id = get_free_concept_id();

				if (!try_init_concept(constant_id)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}

				constant = Formula::new_constant(constant_id);
				if (constant == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free_concept_id(constant_id); return NULL;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_concept_id(constant_id); return NULL;
				}

				Proof* proof = make_proof<Contradiction, false, ResolveInconsistencies>(substituted, new_constant, std::forward<Args>(args)...);
				free(*substituted); if (substituted->reference_count == 0) free(substituted);
				if (proof != NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return proof;
				}
				free(*constant); if (constant->reference_count == 0) free(constant);
				if (ground_concepts[constant_id - new_constant_offset].types.keys != NULL)
					/* `make_proof` could have already freed the new concept */
					free_concept_id(constant_id);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
			}
		}

		finished_constants(quantified, original_constant_count, std::forward<Args>(args)...);
		free(*var); if (var->reference_count == 0) free(var);
		return NULL;
	}

	template<bool DefinitionsAllowed, bool Negated, bool ResolveInconsistencies, typename... Args>
	Proof* make_atom_proof(
			unsigned int predicate,
			Term* arg1, Term* arg2,
			unsigned int& new_constant,
			Args&&... args)
	{
		if (arg1->type == TermType::CONSTANT && arg2 == NULL)
		{
			/* this is a unary formula */
			if (predicate != (unsigned int) built_in_predicates::UNKNOWN) {
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					/* this is a definition of an object */
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_atom_proof ERROR: Definitions are not allowed in this context.\n");
						return NULL;
					}
					new_constant = get_free_concept_id();

					Formula* formula = Negated ?
							Formula::new_not(Formula::new_atom(predicate, Formula::new_constant(new_constant))) :
							Formula::new_atom(predicate, Formula::new_constant(new_constant));
					if (formula == NULL) return NULL;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;

					if (!try_init_concept(new_constant)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
					} if (!add_unary_atom<Negated, ResolveInconsistencies>(predicate, new_constant, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom);
						free_concept_id(new_constant); return NULL;
					}
					return new_axiom;
				} else {
					/* this is a formula of form `t(c)` */

					/* first check if there is already an extensional edge that proves this formula (for this particular instantiation) */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;
						const Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type == FormulaType::NOT) {
							if (Negated) set_formula = set_formula->unary.operand;
							else continue;
						} else if (Negated) continue;
						if (set_formula->type != FormulaType::UNARY_APPLICATION)
							continue;

						unsigned int expected_constant_eliminator = 0;
						if (set_formula->binary.left->type == TermType::VARIABLE) {
							expected_constant_eliminator = predicate;
						} else if (set_formula->binary.left->type != TermType::CONSTANT || set_formula->binary.left->constant != predicate) {
							continue;
						}

						if (set_formula->binary.right->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == 0)
								expected_constant_eliminator = arg1->constant;
							else if (expected_constant_eliminator != arg1->constant)
								continue;
						} else if (*set_formula->binary.right != *arg1) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							if (entry.value->type == ProofType::UNIVERSAL_ELIMINATION
							 && entry.value->operands[1]->type == ProofType::TERM_PARAMETER
							 && entry.value->operands[1]->term->type == TermType::CONSTANT
							 && entry.value->operands[1]->term->constant == expected_constant_eliminator)
							{
								for (Proof* grandchild : entry.value->children) {
									if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
										/* we found a proof of the atomic formula */
										grandchild->reference_count++;
										return grandchild;
									}
								}
							}
						}
					}

					/* there is no extensional edge that proves this atomic formula,
					  so check if the axiom already exists, and if not, create a new axiom */
					bool contains;
					concept<ProofCalculus>& c = ground_concepts[arg1->constant - new_constant_offset];
					Proof* axiom = (Negated ?
							c.negated_types.get(predicate, contains) :
							c.types.get(predicate, contains));
					if (contains) {
						axiom->reference_count++;
						return axiom;
					}
					Formula* formula = Negated ?
							Formula::new_not(Formula::new_atom(predicate, arg1)) :
							Formula::new_atom(predicate, arg1);
					if (formula == NULL) return NULL;
					arg1->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated, ResolveInconsistencies>(predicate, arg1->constant, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else if (arg1->type == TermType::CONSTANT
				&& arg2->type == TermType::CONSTANT)
		{
			/* this is a binary formula */
			if (predicate != (unsigned int) built_in_predicates::UNKNOWN)
			{
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN
				 || arg2->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
					return NULL;
				} else {
					/* this is a formula of form `r(c_1,c_2)` */

					/* first check if there is already an extensional edge that proves this formula (for this particular instantiation) */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;
						const Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type == FormulaType::NOT) {
							if (Negated) set_formula = set_formula->unary.operand;
							else continue;
						} else if (Negated) continue;
						if (set_formula->type != FormulaType::BINARY_APPLICATION)
							continue;

						unsigned int expected_constant_eliminator = 0;
						if (set_formula->ternary.first->type == TermType::VARIABLE) {
							expected_constant_eliminator = predicate;
						} else if (set_formula->ternary.first->type != TermType::CONSTANT || set_formula->ternary.first->constant != predicate) {
							continue;
						}

						if (set_formula->ternary.second->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == 0)
								expected_constant_eliminator = arg1->constant;
							else if (expected_constant_eliminator != arg1->constant)
								continue;
						} else if (*set_formula->ternary.second != *arg1) {
							continue;
						}

						if (set_formula->ternary.third->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == 0)
								expected_constant_eliminator = arg2->constant;
							else if (expected_constant_eliminator != arg2->constant)
								continue;
						} else if (*set_formula->ternary.third != *arg2) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							if (entry.value->type == ProofType::UNIVERSAL_ELIMINATION
							 && entry.value->operands[1]->type == ProofType::TERM_PARAMETER
							 && entry.value->operands[1]->term->type == TermType::CONSTANT
							 && entry.value->operands[1]->term->constant == expected_constant_eliminator)
							{
								for (Proof* grandchild : entry.value->children) {
									if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
										/* we found a proof of the atomic formula */
										grandchild->reference_count++;
										return grandchild;
									}
								}
							}
						}
					}

					/* there is no extensional edge that proves this atomic formula,
					   so check if the axiom already exists, and if not, create a new axiom */
					bool contains;
					concept<ProofCalculus>& c = ground_concepts[arg1->constant - new_constant_offset];
					Proof* axiom = (Negated ?
							c.negated_relations.get(relation(predicate, 0, arg2->constant), contains) :
							c.relations.get(relation(predicate, 0, arg2->constant), contains));
					if (contains) {
						axiom->reference_count++;
						return axiom;
					}
					Formula* formula = Negated ?
							Formula::new_not(Formula::new_atom(predicate, arg1, arg2)) :
							Formula::new_atom(predicate, arg1, arg2);
					if (formula == NULL) return NULL;
					arg1->reference_count++; arg2->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated, ResolveInconsistencies>({predicate, arg1->constant, arg2->constant}, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else {
			fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
			return NULL;
		}
	}

	void free_proof(Proof* proof) {
		theory::changes changes;
		if (!get_theory_changes(*proof, changes)) return;
		subtract_changes(changes);
		free(*proof); if (proof->reference_count == 0) free(proof);
	}

	inline bool valid_definition(const Formula* right, unsigned int quantified_variable)
	{
		unsigned int predicate; Term const* arg1; Term const* arg2;
		if (is_atomic(*right, predicate, arg1, arg2)) {
			return predicate != (unsigned int) built_in_predicates::UNKNOWN
				&& arg1->type == TermType::VARIABLE
				&& arg1->variable == quantified_variable
				&& arg2 == NULL;
		} else if (right->type == FormulaType::AND) {
			for (unsigned int i = 0; i < right->array.length; i++)
				if (valid_definition(right->array.operands[i], quantified_variable)) return true;
			return false;
		} else {
			return false;
		}
	}

	template<bool Negated>
	bool make_conjunct_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		unsigned int predicate; Term* arg1; Term* arg2;
		if (is_atomic(*constraint, predicate, arg1, arg2)) {
			bool contains;
			if (arg2 == NULL) {
				/* `constraint` is a unary literal */
				proof = (Negated ? c.negated_types.get(predicate, contains) : c.types.get(predicate, contains));
				if (!contains) proof = NULL;
				return true;
			} else {
				/* `constraint` is a binary literal */
				relation r = { predicate,
						(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
						(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
				proof = (Negated ? c.negated_relations.get(r, contains) : c.relations.get(r, contains));
				if (!contains) proof = NULL;
				return true;
			}
		} else if (constraint->type == FormulaType::NOT) {
			return make_conjunct_proof<true>(constraint->unary.operand, c, proof);
		} else {
			return false;
		}
	}

	bool make_universal_elim_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		if (constraint->type == FormulaType::AND) {
			Proof** conjunct_proofs = (Proof**) malloc(sizeof(Proof*) * constraint->array.length);
			if (conjunct_proofs == NULL) {
				fprintf(stderr, "theory.make_universal_elim_proof ERROR: Out of memory.\n");
				return false;
			}
			for (unsigned int i = 0; i < constraint->array.length; i++) {
				Formula* conjunct = constraint->array.operands[i];
				if (!make_conjunct_proof<false>(conjunct, c, conjunct_proofs[i])) return false;
				if (conjunct_proofs[i] == NULL) {
					proof = NULL; free(conjunct_proofs);
					return true;
				}
			}
			proof = ProofCalculus::new_conjunction_intro(make_array_view(conjunct_proofs, constraint->array.length));
			free(conjunct_proofs);
			return proof != NULL;
		} else {
			return make_conjunct_proof<false>(constraint, c, proof);
		}
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool move_element_to_subset(unsigned int element, Formula* lifted_literal, Args&&... visitor)
	{
		Formula* new_set_formula; bool contains;
		unsigned int old_set_id = sets.element_map.get(element, contains);
		Formula* old_set_formula = (contains ? sets.sets[old_set_id].set_formula() : NULL);
		if (old_set_formula == NULL) {
			/* there are no other atomic formulas for this concept */
			new_set_formula = lifted_literal;
			new_set_formula->reference_count++;
		} else {
			new_set_formula = Formula::new_and(old_set_formula, lifted_literal);
			old_set_formula->reference_count++;
			lifted_literal->reference_count++;
		}
		if (new_set_formula == NULL)
			return false;
		Formula* canonicalized = Canonicalizer::canonicalize(*new_set_formula);
		free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		if (canonicalized == NULL)
			return false;
		if (!sets.template move_element_to_set<ResolveInconsistencies>(element, canonicalized, std::forward<Args>(visitor)...)) {
			free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
			return false;
		}

		/* make sure we can't currently prove the negation of this new axiom
		   via set containment (can we prove that the instance does not belong to any ancestor?) */
		unsigned int new_set_id = sets.set_ids.get(*canonicalized);
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == NULL) continue;
			Formula* set_formula = sets.sets[i].set_formula();
			Formula* negated_set_formula = Formula::new_not(set_formula);
			if (negated_set_formula == NULL) {
				if (old_set_formula == NULL) sets.remove_element(element);
				else sets.move_element_to_superset(element, old_set_formula);
				return false;
			}
			set_formula->reference_count++;
			bool contradiction = is_subset(lifted_literal, negated_set_formula) && sets.sets[i].descendants.contains(new_set_id);
			free(*negated_set_formula); free(negated_set_formula);
			if (contradiction) {
				if (old_set_formula == NULL) sets.remove_element(element);
				else sets.move_element_to_superset(element, old_set_formula);
				return false;
			}
		}

		return true;
	}

	bool subtract_conjuncts(
		Formula** first, unsigned int first_count,
		Formula** second, unsigned int second_count,
		array<Formula*>& difference)
	{
		unsigned int i = 0, j = 0;
		while (i < first_count && j < second_count)
		{
			if (*first[i] == *second[j]) {
				i++; j++;
			} else if (*first[i] < *second[j]) {
				if (!difference.add(first[i])) return false;
				i++;
			} else {
				j++;
			}
		}

		while (i < first_count) {
			if (!difference.add(first[i])) return false;
			i++;
		}
		return true;
	}

	template<typename... Args>
	bool move_element_to_superset(unsigned int element, Formula* lifted_literal, Args&&... visitor)
	{
		Formula* old_set_formula = sets.sets[sets.element_map.get(element)].set_formula();

		Formula* canonicalized = Canonicalizer::canonicalize(*lifted_literal);
		array<Formula*> difference(8);
		if (canonicalized == NULL) {
			return false;
		} else if (old_set_formula->type != FormulaType::AND) {
			if (canonicalized->type == FormulaType::AND) {
				if (!subtract_conjuncts(&old_set_formula, 1, canonicalized->array.operands, canonicalized->array.length, difference)) {
					free(*canonicalized); free(canonicalized); return false;
				}
			} else {
				if (!subtract_conjuncts(&old_set_formula, 1, &canonicalized, 1, difference)) {
					free(*canonicalized); free(canonicalized); return false;
				}
			}
		} else if (canonicalized->type == FormulaType::AND) {
			if (!subtract_conjuncts(old_set_formula->array.operands, old_set_formula->array.length, canonicalized->array.operands, canonicalized->array.length, difference)) {
				free(*canonicalized); free(canonicalized); return false;
			}
		} else {
			if (!subtract_conjuncts(old_set_formula->array.operands, old_set_formula->array.length, &canonicalized, 1, difference)) {
				free(*canonicalized); free(canonicalized); return false;
			}
		}
		free(*canonicalized); free(canonicalized);

		Formula* new_set_formula;
		if (difference.length == 0)
			new_set_formula = NULL;
		else if (difference.length == 1)
			new_set_formula = difference[0];
		else new_set_formula = Formula::new_and(make_array_view(difference.data, difference.length));
		if (difference.length != 0 && new_set_formula == NULL) {
			return false;
		}
		for (Formula* conjunct : difference)
			conjunct->reference_count++;

		bool success;
		if (new_set_formula == NULL) {
			sets.remove_element(element, std::forward<Args>(visitor)...);
			success = true;
		} else {
			success = sets.move_element_to_superset(element, new_set_formula, std::forward<Args>(visitor)...);
			free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		}
		return success;
	}
};

template<typename Formula>
constexpr bool filter_constants(const Formula* formula, const array<unsigned int>& constants) {
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int constant) {
	return true;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count) { }

template<bool Negated, typename Formula>
bool intensional_element_of(
		unsigned int predicate,
		const typename Formula::Term* arg1,
		const typename Formula::Term* arg2,
		const Formula* set_formula,
		typename Formula::Term*& unifying_substitution)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	switch (set_formula->type) {
	case FormulaType::UNARY_APPLICATION:
		if (Negated || arg1 == NULL || arg2 != NULL) return false;
		if (set_formula->binary.left->type == TermType::VARIABLE && set_formula->binary.left->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = Term::new_constant(predicate);
				if (unifying_substitution == NULL) return false;
			} else if (unifying_substitution->type != TermType::CONSTANT || unifying_substitution->constant != predicate) {
				return false;
			}
		} else if (set_formula->binary.left->type != TermType::CONSTANT || set_formula->binary.left->constant != predicate) {
			return false;
		}

		if (set_formula->binary.right->type == TermType::VARIABLE && set_formula->binary.right->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg1;
				arg1->reference_count++;
			} else if (*unifying_substitution != *arg1) {
				return false;
			}
		} else if (*set_formula->binary.right != *arg1) {
			return false;
		}
		return true;
	case FormulaType::BINARY_APPLICATION:
		if (Negated || arg1 == NULL || arg2 == NULL) return false;
		if (set_formula->ternary.first->type == TermType::VARIABLE && set_formula->ternary.first->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = Term::new_constant(predicate);
				if (unifying_substitution == NULL) return false;
			} else if (unifying_substitution->type != TermType::CONSTANT || unifying_substitution->constant != predicate) {
				return false;
			}
		} else if (set_formula->ternary.first->type != TermType::CONSTANT || set_formula->ternary.first->constant != predicate) {
			return false;
		}

		if (set_formula->ternary.second->type == TermType::VARIABLE && set_formula->ternary.second->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg1;
				arg1->reference_count++;
			} else if (*unifying_substitution != *arg1) {
				return false;
			}
		} else if (*set_formula->ternary.second != *arg1) {
			return false;
		}

		if (set_formula->ternary.third->type == TermType::VARIABLE && set_formula->ternary.third->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg2;
				arg2->reference_count++;
			} else if (*unifying_substitution != *arg2) {
				return false;
			}
		} else if (*set_formula->ternary.third != *arg2) {
			return false;
		}
		return true;
	case FormulaType::NOT:
		return intensional_element_of<!Negated>(predicate, arg1, arg2, set_formula, unifying_substitution);
	case FormulaType::AND:
		for (unsigned int i = 0; i < set_formula->array.length; i++)
			if (!intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->array.operands[i], unifying_substitution)) return false;
		return true;
	case FormulaType::OR:
		for (unsigned int i = 0; i < set_formula->array.length; i++) {
			if (intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->array.operands[i], unifying_substitution))
				return true;
			if (unifying_substitution != NULL) {
				free(*unifying_substitution);
				if (unifying_substitution->reference_count == 0)
					free(unifying_substitution);
				unifying_substitution = NULL;
			}
		}
		return false;
	case FormulaType::IF_THEN:
		return intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->binary.left, unifying_substitution)
			&& intensional_element_of<!Negated>(predicate, arg1, arg2, set_formula->binary.right, unifying_substitution);
	case FormulaType::TRUE:
		return !Negated;
	case FormulaType::FALSE:
		return Negated;
	case FormulaType::EQUALS:
	case FormulaType::IFF:
	case FormulaType::FOR_ALL:
	case FormulaType::EXISTS:
	case FormulaType::LAMBDA:
		fprintf(stderr, "intensional_element_of ERROR: Not implemented.\n");
		return false;
	case FormulaType::VARIABLE:
	case FormulaType::CONSTANT:
	case FormulaType::PARAMETER:
	case FormulaType::INTEGER:
		break;
	}
	fprintf(stderr, "intensional_element_of ERROR: Unrecognized formula type.\n");
	return false;
}

template<typename Formula, typename ProofCalculus, typename Canonicalizer>
Formula* make_lifted_conjunction(unsigned int concept_id,
		const theory<Formula, ProofCalculus, Canonicalizer>& T)
{
	const concept<ProofCalculus>& c = T.ground_concepts.get(concept_id);
	array<Formula*> conjuncts(16);
	for (const auto& entry : c.types) {
		Formula* new_conjunct = Formula::new_atom(entry.key, Formula::new_variable(1));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
	} for (const auto& entry : c.negated_types) {
		Formula* new_conjunct = Formula::new_not(Formula::new_atom(entry.key, Formula::new_variable(1)));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
	} for (const auto& entry : c.relations) {
		Formula* new_conjunct = Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? Formula::new_variable(1) : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? Formula::new_variable(1) : Formula::new_constant(entry.key.arg2));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
	} for (const auto& entry : c.negated_relations) {
		Formula* new_conjunct = Formula::new_not(Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? Formula::new_variable(1) : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? Formula::new_variable(1) : Formula::new_constant(entry.key.arg2)));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
	}
	Formula* conjunction = Formula::new_and(conjuncts);
	return conjunction;
}

#endif /* THEORY_H_ */
