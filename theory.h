#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <math/multiset.h>
#include <stdexcept>

#if !defined(NDEBUG)
#include <functional>
#endif

#include "array_view.h"
#include "set_reasoning.h"
#include "built_in_predicates.h"

using namespace core;


template<typename Proof>
inline void visit_node(const Proof& proof) { }

template<bool Negated, typename Term>
constexpr bool visit_unary_atom(const Term* term) { return true; }

template<bool Negated>
constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2) { return true; }

template<typename Proof>
constexpr bool visit_subset_axiom(const Proof& proof) { return true; }

constexpr bool visit_existential_intro() { return true; }
constexpr bool visit_negated_conjunction() { return true; }
constexpr bool visit_negated_universal_intro() { return true; }
constexpr bool visit_disjunction_intro() { return true; }
inline void on_subtract_changes() { }

enum class instance_type {
	ANY,
	CONSTANT,
	INTEGER,
	STRING
};

struct instance {
	instance_type type;
	union {
		unsigned int constant;
		int64_t integer;
		string* str;
	};

	static inline void swap(instance& first, instance& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(instance); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void move(const instance& src, instance& dst) {
		dst.type = src.type;
		switch (src.type) {
		case instance_type::ANY: return;
		case instance_type::CONSTANT: dst.constant = src.constant; return;
		case instance_type::INTEGER: dst.integer = src.integer; return;
		case instance_type::STRING: dst.str = src.str; return;
		}
		fprintf(stderr, "instance.move ERROR: Unrecognized instance_type.\n");
		exit(EXIT_FAILURE);
	}
};

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
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Term Term;

	array_map<Term, Proof*> types;
	array_map<Term, Proof*> negated_types;
	array_map<relation, Proof*> relations;
	array_map<relation, Proof*> negated_relations;
	array<Proof*> definitions;
	array<Proof*> existential_intro_nodes;

	array_map<unsigned int, Proof*> function_values;

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, const char* prefix, Printer&&... printer) const {
		for (auto entry : types) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_types) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : relations) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_relations) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (Proof* definition : definitions) {
			if (!print(prefix, out) || !print(*definition->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : function_values) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		}
		return true;
	}

	bool check_axioms() const {
		bool success = true;
		for (auto entry : types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : negated_types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : relations) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : negated_relations) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		}
		return success;
	}

	static inline void move(const concept<ProofCalculus>& src, concept<ProofCalculus>& dst) {
		core::move(src.types, dst.types);
		core::move(src.negated_types, dst.negated_types);
		core::move(src.relations, dst.relations);
		core::move(src.negated_relations, dst.negated_relations);
		core::move(src.definitions, dst.definitions);
		core::move(src.existential_intro_nodes, dst.existential_intro_nodes);
		core::move(src.function_values, dst.function_values);
	}

	static inline void free(concept<ProofCalculus>& c) {
		/* we set the initial reference_count to 2 since `theory.free_proof`
		   will free definitions when their reference_count is 1 */
		core::free(*c.definitions[0]);
		for (auto entry : c.types) {
			core::free(entry.key);
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.negated_types) {
			core::free(entry.key);
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
		} for (Proof* definition : c.definitions) {
			core::free(*definition);
			if (definition->reference_count == 0)
				core::free(definition);
		} for (auto entry : c.function_values) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}
		core::free(c.types);
		core::free(c.negated_types);
		core::free(c.relations);
		core::free(c.negated_relations);
		core::free(c.definitions);
		core::free(c.existential_intro_nodes);
		core::free(c.function_values);
	}
};

template<typename ProofCalculus>
inline bool init(concept<ProofCalculus>& c, unsigned int constant) {
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;

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
	} else if (!array_init(c.definitions, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		return false;
	} else if (!array_init(c.existential_intro_nodes, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); return false;
	} else if (!array_map_init(c.function_values, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.existential_intro_nodes);
		return false;
	}

	Formula* constant_expr = Formula::new_constant(constant);
	if (constant_expr == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(c.existential_intro_nodes); return false;
	}
	Formula* axiom = Formula::new_equals(constant_expr,
			Formula::new_lambda(1, Formula::new_apply(constant_expr, &Term::template variables<1>::value)));
	if (constant_expr == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(*constant_expr); free(constant_expr);
		free(c.existential_intro_nodes); return false;
	}
	Term::template variables<1>::value.reference_count++;
	constant_expr->reference_count += 2 - 1;
	c.definitions[0] = ProofCalculus::new_axiom(axiom);
	if (c.definitions[0] == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(*axiom); free(axiom);
		free(c.existential_intro_nodes); return false;
	}
	free(*axiom);
	/* we set the initial reference_count to 2 since `theory.free_proof`
	   will free definitions when their reference_count is 1 */
	c.definitions[0]->reference_count += 2;
	c.definitions.length++;
	return true;
}

template<typename Formula>
inline Formula* preprocess_formula(Formula* src) {
	Formula* first = remove_inverse(src);
	if (first == nullptr) return nullptr;

	Formula* second = same_to_equals(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = remove_exists(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = arg_of_to_arg(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = remove_tense(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = remove_object(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = normalize_set_operations(second);
	free(*second); if (second->reference_count == 0) free(second);
	return first;
}

/* this is useful in `theory.add_definition` where if any sets are created in
   `set_reasoning.get_subset_axiom`, we need the sets to have a specific size */
struct required_set_size {
	unsigned int set_size;
	required_set_size(unsigned int set_size) : set_size(set_size) { }
};

template<typename... Args>
inline bool compute_new_set_size(
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		required_set_size& required,
		Args&&... visitor)
{
	if (required.set_size == UINT_MAX) {
		if (!compute_new_set_size(out, min_set_size, max_set_size, std::forward<Args>(visitor)...))
			return false;
		required.set_size = out;
		return true;
	} else {
		return required.set_size >= min_set_size && required.set_size <= max_set_size
			&& compute_new_set_size(out, required.set_size, required.set_size, std::forward<Args>(visitor)...);
	}
}

template<typename Proof, typename... Args>
inline bool compute_new_set_size(unsigned int& out,
		unsigned int min_set_size, unsigned int max_set_size,
		const array<Proof*>& freeable_axioms,
		Args&&... visitor)
{
	return compute_new_set_size(out, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus>& sets,
		const array<typename ProofCalculus::Proof*>& freeable_axioms,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
}

template<typename ProofCalculus, typename Canonicalizer>
struct theory
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	template<unsigned int Value> using Constants = typename Term::template constants<Value>;
	template<unsigned int Value> using Variables = typename Term::template variables<Value>;

	unsigned int new_constant_offset;

	/* A map from `t` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]t` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]t`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]t` or `~[u/0]t`
	   are in the theory. */
	hash_map<Term, pair<array<unsigned int>, array<unsigned int>>> atoms;

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
	set_reasoning<built_in_predicates, ProofCalculus> sets;

	array<pair<Formula*, Proof*>> disjunction_intro_nodes;
	array<pair<Formula*, Proof*>> negated_conjunction_nodes;
	array<pair<Formula*, Proof*>> implication_intro_nodes;
	array<pair<Formula*, Proof*>> existential_intro_nodes;

	Proof* empty_set_axiom;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), atoms(64), relations(64),
			ground_concept_capacity(64), ground_axiom_count(0), observations(64),
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
			Formula::new_equals(Formula::new_atom((unsigned int) built_in_predicates::SIZE, &Variables<1>::value), Formula::new_int(0)),
			Formula::new_not(Formula::new_exists(2, Formula::new_apply(&Variables<1>::value, &Variables<2>::value)))
		));
		if (empty_set_formula == NULL) {
			core::free(ground_concepts);
			exit(EXIT_FAILURE);
		}
		Variables<1>::value.reference_count += 2;
		Variables<2>::value.reference_count++;
		empty_set_axiom = ProofCalculus::new_axiom(empty_set_formula);
		core::free(*empty_set_formula);
		if (empty_set_formula->reference_count == 0)
			core::free(empty_set_formula);
		if (empty_set_axiom == NULL) exit(EXIT_FAILURE);
		empty_set_axiom->reference_count++;
	}

	~theory() {
		for (auto entry : atoms) {
			core::free(entry.key);
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

		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;

			auto& c = ground_concepts[i];
			for (unsigned int i = 0; i < c.definitions.length; i++) {
				Proof* definition = c.definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					for (unsigned int j = i + 1; j < c.definitions.length; j++) {
						Proof* other_definition = c.definitions[j];
						if (other_definition->formula->binary.right->type != FormulaType::LAMBDA) continue;

						Proof* axiom = sets.template get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, 1);
						free(*axiom);
						if (axiom->reference_count == 1)
							sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand);

						axiom = sets.template get_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1);
						free(*axiom);
						if (axiom->reference_count == 1)
							sets.free_subset_axiom(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand);
					}
				}
			}

			core::free(c);
		}
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
		unsigned int old_capacity = ground_concept_capacity;
		ground_concept_capacity = new_capacity;
		return new_constant_offset + old_capacity;
	}

	inline bool try_init_concept(unsigned int id) {
		if (id < new_constant_offset) return true;
		concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.keys != NULL)
			return true;
		else return init(c, id);
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
		if (c.types.size != 0 || c.negated_types.size != 0 || c.relations.size != 0
		 || c.negated_relations.size != 0 || c.definitions.length > 1 || c.definitions[0]->reference_count > 2
		 || c.existential_intro_nodes.length != 0 || c.function_values.size != 0)
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
			if (ground_concepts[i].types.keys == NULL) continue;
			if (!print("Concept ", out) || !print(i + new_constant_offset, out, std::forward<Printer>(printer)...) || !print(":\n", out)
			 || !ground_concepts[i].print_axioms(out, "  ", std::forward<Printer>(printer)...)) return false;
		}
		return sets.print_axioms(out, std::forward<Printer>(printer)...);
	}

	Proof* add_formula(Formula* formula, unsigned int& new_constant)
	{
		new_constant = 0;

		Formula* new_formula = preprocess_formula(formula);
		if (new_formula == NULL) return nullptr;

		Formula* canonicalized = Canonicalizer::canonicalize(*new_formula);
		free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
		if (canonicalized == NULL) return nullptr;

/* TODO: for debugging; delete this */
print("canonicalized: ", stderr); print(*canonicalized, stderr); print('\n', stderr);
		Proof* new_proof = make_proof<false, true, true>(canonicalized, new_constant);
print_axioms(stderr);
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
			return nullptr;
		} else if (!observations.add(new_proof)) {
			free_proof(new_proof);
			return nullptr;
		}
		return new_proof;
	}

	void remove_formula(Proof* proof) {
		observations.remove(proof);
		free_proof(proof);
	}

	template<ProofType Type>
	bool check_disjunction_introductions(const array<pair<Formula*, Proof*>>& proof_steps) const
	{
		bool success = true;
		array<const Proof*> computed_steps(64);
		for (const Proof* proof : observations) {
			array<const Proof*> steps(16);
			if (!get_proof_steps<Type>(proof, steps)) {
				fprintf(stderr, "theory.check_proof_disjunctions ERROR: `get_proof_steps` failed.\n");
				return false;
			}

			for (const Proof* step : steps) {
				bool contains = false;
				for (const auto& entry : proof_steps) {
					if (*step == *entry.value) {
						contains = true;
						break;
					}
				}
				if (!contains) {
					print("theory.check_proof_disjunctions WARNING: Found step "
							"of a proof in `observations` that is not in `proof_steps`.\n", stderr);
					success = false;
				}
			}

			for (const Proof* step : steps) {
				if (!computed_steps.contains(step) && !computed_steps.add(step))
					return false;
			}
		}

		for (const auto& entry : proof_steps) {
			if (!computed_steps.contains(entry.value)) {
				print("theory.check_proof_disjunctions WARNING: `proof_steps` "
						"contains a step that does not belong to a proof in `observations`.\n", stderr);
				print("  Formula: ", stderr); print(*entry.key, stderr); print('\n', stderr);
				success = false;
			}
		}

		return success;
	}

	inline bool check_disjunction_introductions() const {
		return check_disjunction_introductions<ProofType::DISJUNCTION_INTRODUCTION>(disjunction_intro_nodes)
			&& check_disjunction_introductions<ProofType::EXISTENTIAL_INTRODUCTION>(existential_intro_nodes)
			&& check_disjunction_introductions<ProofType::IMPLICATION_INTRODUCTION>(implication_intro_nodes);
	}

	inline bool check_concept_axioms() const {
		bool success = true;
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;
			success &= ground_concepts[i].check_axioms();
		}
		return success;
	}

	inline bool is_provably_not_a_set(unsigned int constant) const {
		if (constant < new_constant_offset || ground_concepts[constant - new_constant_offset].types.keys == NULL)
			return false;
		return ground_concepts[constant - new_constant_offset].types.size != 0
			|| ground_concepts[constant - new_constant_offset].function_values.size != 0;
	}

	inline bool is_provably_a_set(unsigned int constant) const {
		if (constant < new_constant_offset || ground_concepts[constant - new_constant_offset].types.keys == NULL)
			return false;
		const Formula* first_set_definition = ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.right->quantifier.operand;
		if (sets.set_ids.table.contains(*first_set_definition))
			return true;

		Term* atom = Term::new_apply(Term::new_constant(constant), &Variables<1>::value);
		if (atom == nullptr)
			return false;
		Variables<1>::value.reference_count++;

		bool contains;
		pair<array<unsigned int>, array<unsigned int>>& type_instances = atoms.get(*atom, contains);
		free(*atom); free(atom);
		if (contains && type_instances.key.length != 0)
			return true;
		return false;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_unary_atom(const Term& atom, Proof* axiom, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (atom.binary.right->type != TermType::CONSTANT)
			fprintf(stderr, "add_unary_atom WARNING: The operand of this application is not constant.\n");
#endif
		unsigned int arg = atom.binary.right->constant;

		Formula* lifted_literal;
		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom);
				return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}
		if (!move_element_to_subset<ResolveInconsistencies>(arg, lifted_literal, std::forward<Args>(visitor)...)) {
			free(*lifted_literal); free(lifted_literal);
			return false;
		}

#if !defined(NDEBUG)
		if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_unary_atom WARNING: `ground_concepts` does not contain the concept %u.\n", arg);
#endif
		bool contains; unsigned int bucket;
		if (!atoms.check_size()) {
			free(*lifted_literal); free(lifted_literal);
			return false;
		}

		pair<array<unsigned int>, array<unsigned int>>& instance_pair = atoms.get(*lifted_literal, contains, bucket);
		if (!contains) {
			if (!array_init(instance_pair.key, 8)) {
				free(*lifted_literal); free(lifted_literal);
				return false;
			} else if (!array_init(instance_pair.value, 8)) {
				core::free(instance_pair.key);
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			atoms.table.keys[bucket] = *lifted_literal;
			atoms.table.size++;
		}

		array<unsigned int>& instances = (Negated ? instance_pair.value : instance_pair.key);
		array_map<Term, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);
		if (!instances.ensure_capacity(instances.length + 1)
		 || !ground_types.ensure_capacity(ground_types.size + 1))
		{
			free(*lifted_literal); free(lifted_literal);
			return false;
		}

		add_sorted<false>(instances, arg);
		ground_types.keys[ground_types.size] = *lifted_literal;
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		ground_axiom_count++;
		free(*lifted_literal); free(lifted_literal);
		return true;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_binary_atom(relation rel, Proof* axiom, Args&&... visitor)
	{
		/* check if this observation is inconsistent with the theory (can we prove its negation?) */
		if (rel.arg1 == rel.arg2) {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_and(
					Formula::new_atom(rel.predicate, &Variables<1>::value, Formula::new_constant(rel.arg2)),
					Formula::new_atom(rel.predicate, &Variables<1>::value, &Variables<1>::value),
					Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), &Variables<1>::value));
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count += 4;
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
			Formula* lifted_atom = Formula::new_atom(rel.predicate, &Variables<1>::value, Formula::new_constant(rel.arg2));
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count++;
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

			lifted_atom = Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), &Variables<1>::value);
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count++;
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
	inline bool remove_unary_atom(const Term& atom, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (atom.binary.right->type != TermType::CONSTANT)
			fprintf(stderr, "remove_unary_atom WARNING: The operand of this application is not constant.\n");
#endif
		unsigned int arg = atom.binary.right->constant;

		Formula* lifted_literal;
		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom);
				return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}
		if (!move_element_to_superset(arg, lifted_literal, std::forward<Args>(visitor)...)) {
			fprintf(stderr, "theory.remove_unary_atom ERROR: `move_element_to_superset` failed.\n");
			if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
			return false;
		}

#if !defined(NDEBUG)
		if (!atoms.table.contains(*lifted_literal)) {
			print("theory.remove_unary_atom WARNING: `atoms` does not contain the key ", stderr);
			print(*lifted_literal, stderr); print(".\n", stderr);
		} if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", arg);
#endif

		array<unsigned int>& instances = (Negated ? atoms.get(*lifted_literal).value : atoms.get(*lifted_literal).key);
		array_map<Term, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);

		unsigned int index = instances.index_of(arg);
#if !defined(NDEBUG)
		if (index == instances.length)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `instances` does not contain %u.\n", arg);
#endif
		shift_left(instances.data + index, instances.length - index - 1);
		instances.length--;

		index = ground_types.index_of(*lifted_literal);
#if !defined(NDEBUG)
		if (index == ground_types.size) {
			print("theory.remove_unary_atom WARNING: `ground_types` does not contain ", stderr);
			print(*lifted_literal, stderr); print(".\n", stderr);
		}
#endif
		Proof* axiom = ground_types.values[index];
		free(*axiom); if (axiom->reference_count == 0) free(axiom);
		ground_types.remove_at(index);
		ground_axiom_count--;
		free(*lifted_literal); free(lifted_literal);
		return true;
	}

	template<bool Negated, typename... Args>
	inline bool remove_binary_atom(relation rel, Args&&... visitor)
	{
		if (rel.arg1 == rel.arg2) {
			Formula* lifted_literal;
			Formula* lifted_atom = Formula::new_and(
					Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), &Variables<1>::value),
					Formula::new_atom(rel.predicate, &Variables<1>::value, &Variables<1>::value),
					Formula::new_atom(rel.predicate, &Variables<1>::value, Formula::new_constant(rel.arg2)));
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count += 4;
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
			Formula* lifted_atom = Formula::new_atom(rel.predicate, &Variables<1>::value, Formula::new_constant(rel.arg2));
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count++;
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

			lifted_atom = Formula::new_atom(rel.predicate, Formula::new_constant(rel.arg1), &Variables<1>::value);
			if (lifted_atom == NULL) return false;
			Variables<1>::value.reference_count++;
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

	enum class change_type {
		UNARY_ATOM,
		NEGATED_UNARY_ATOM,
		BINARY_ATOM,
		NEGATED_BINARY_ATOM,
		SUBSET_AXIOM,
		SET_SIZE_AXIOM,
		DEFINITION,
		FUNCTION_VALUE,
		IMPLICATION_INTRO_NODE,
		NEGATED_CONJUNCTION_NODE,
		EXISTENTIAL_INTRO_NODE,
		DISJUNCTION_INTRO_NODE
	};

	struct change {
		change_type type;
		union {
			pair<Term, Proof*> unary_atom;
			pair<relation, Proof*> binary_atom;
			Proof* axiom;
			pair<Formula*, Proof*> intro_node;
		};

		static inline void free(change& c) {
			switch (c.type) {
			case change_type::UNARY_ATOM:
			case change_type::NEGATED_UNARY_ATOM:
				core::free(c.unary_atom.key);
				core::free(*c.unary_atom.value); if (c.unary_atom.value->reference_count == 0) core::free(c.unary_atom.value);
				return;
			case change_type::BINARY_ATOM:
			case change_type::NEGATED_BINARY_ATOM:
				core::free(*c.binary_atom.value); if (c.binary_atom.value->reference_count == 0) core::free(c.binary_atom.value);
				return;
			case change_type::SUBSET_AXIOM:
			case change_type::SET_SIZE_AXIOM:
			case change_type::DEFINITION:
			case change_type::FUNCTION_VALUE:
				core::free(*c.axiom); if (c.axiom->reference_count == 0) core::free(c.axiom);
				return;
			case change_type::IMPLICATION_INTRO_NODE:
			case change_type::NEGATED_CONJUNCTION_NODE:
			case change_type::EXISTENTIAL_INTRO_NODE:
			case change_type::DISJUNCTION_INTRO_NODE:
				core::free(*c.intro_node.key); if (c.intro_node.key->reference_count == 0) core::free(c.intro_node.key);
				core::free(*c.intro_node.value); if (c.intro_node.value->reference_count == 0) core::free(c.intro_node.value);
				return;
			}
			fprintf(stderr, "theory.change.free ERROR: Unrecognized `change_type`.\n");
			exit(EXIT_FAILURE);
		}
	};

	struct changes {
		array<change> list;

		changes() : list(8) { }

		~changes() {
			for (change& c : list)
				core::free(c);
		}

		inline bool add(change_type type, const pair<Term, Proof*>& unary_atom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length].unary_atom.key = unary_atom.key;
			list[list.length++].unary_atom.value = unary_atom.value;
			unary_atom.value->reference_count++;
			return true;
		}

		inline bool add(change_type type, const pair<relation, Proof*>& binary_atom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].binary_atom = binary_atom;
			binary_atom.value->reference_count++;
			return true;
		}

		inline bool add(change_type type, Proof* axiom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].axiom = axiom;
			axiom->reference_count++;
			return true;
		}

		inline bool add(change_type type, const pair<Formula*, Proof*>& intro_node) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].intro_node = intro_node;
			intro_node.key->reference_count++;
			intro_node.value->reference_count++;
			return true;
		}
	};

	template<typename... Args>
	bool add_changes(const theory::changes& changes, Args&&... visitor)
	{
		for (const change& c : changes.list) {
			array_multiset<unsigned int> constants(8);
			switch (c.type) {
			case change_type::UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<false, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::NEGATED_UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<true, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<false, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<true, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::SUBSET_AXIOM:
				{
					unsigned int variable = c.axiom->formula->quantifier.variable;
					if (c.axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
						unsigned int predicate; Term const* arg1; Term const* arg2;
						Formula* left = c.axiom->formula->quantifier.operand->binary.left;
						if (is_atomic(*left, predicate, arg1, arg2)
						 && arg1->type == TermType::VARIABLE
						 && arg1->variable == variable && arg2 == NULL)
						{
							if (!try_init_concept(predicate)) return false;
						}
					}
					if (!sets.add_subset_axiom(c.axiom, std::forward<Args>(visitor)...)) return false;
				}
				continue;
			case change_type::SET_SIZE_AXIOM:
				{
					unsigned int set_id;
					Formula* set_formula = c.axiom->formula->binary.left->binary.right->quantifier.operand;
					if (!sets.get_set_id(set_formula, set_id, std::forward<Args>(visitor)...)
					 || !sets.sets[set_id].set_size_axiom(c.axiom)) return false;
				}
				continue;
			case change_type::DEFINITION:
				if (!add_definition<false>(c.axiom, std::forward<Args>(visitor)...))
					return false;
				continue;
			case change_type::FUNCTION_VALUE:
				if (!add_function_value<false>(c.axiom))
					return false;
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				if (!implication_intro_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				if (!negated_conjunction_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					if (!existential_intro_nodes.add(c.intro_node)) return false;

					Term* term;
					if (c.intro_node.value->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.value->operands[2]->term;
					else term = c.intro_node.value->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						if (!try_init_concept(term->constant)
						 || !ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.add(c.intro_node.value)) {
							unsigned int index = existential_intro_nodes.index_of(c.intro_node);
							existential_intro_nodes.remove(index);
							return false;
						}
					}
					c.intro_node.key->reference_count++;
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				if (!disjunction_intro_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			}
			fprintf(stderr, "theory.add_changes ERROR: Unrecognized `change_type`.\n");
			return false;
		}
		return true;
	}

	template<typename... Args>
	bool subtract_changes(const theory::changes& changes, Args&&... visitor)
	{
		on_subtract_changes(std::forward<Args>(visitor)...);
		array<Proof*> freeable_set_size_axioms(8);
		for (unsigned int i = changes.list.length; i > 0; i--) {
			const change& c = changes.list[i - 1];
			array_multiset<unsigned int> constants(8);
			switch (c.type) {
			case change_type::UNARY_ATOM:
				remove_unary_atom<false>(c.unary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					try_free_concept_id(constants.counts.keys[i]);
				continue;
			case change_type::NEGATED_UNARY_ATOM:
				remove_unary_atom<true>(c.unary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					try_free_concept_id(constants.counts.keys[i]);
				continue;
			case change_type::BINARY_ATOM:
				remove_binary_atom<false>(c.binary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				try_free_concept_id(c.binary_atom.key.predicate);
				try_free_concept_id(c.binary_atom.key.arg1);
				try_free_concept_id(c.binary_atom.key.arg2);
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				remove_binary_atom<true>(c.binary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				try_free_concept_id(c.binary_atom.key.predicate);
				try_free_concept_id(c.binary_atom.key.arg1);
				try_free_concept_id(c.binary_atom.key.arg2);
				continue;
			case change_type::SUBSET_AXIOM:
				{
					unsigned int variable = c.axiom->formula->quantifier.variable;
					if (c.axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
						unsigned int predicate; Term const* arg1; Term const* arg2;
						Formula* left = c.axiom->formula->quantifier.operand->binary.left;
						if (is_atomic(*left, predicate, arg1, arg2)
						 && arg1->type == TermType::VARIABLE
						 && arg1->variable == variable && arg2 == NULL)
						{
							try_free_concept_id(predicate);
						}
					}
					unsigned int child_count = c.axiom->children.length;
					if (changes.list.length == 1) child_count++;
					if (c.axiom->reference_count == child_count + 2) {
						sets.free_subset_axiom(c.axiom, freeable_set_size_axioms, std::forward<Args>(visitor)...);
					} else {
						free(*c.axiom);
					}
					continue;
				}
			case change_type::SET_SIZE_AXIOM:
				{
					unsigned int set_id;
					Formula* set_formula = c.axiom->formula->binary.left->binary.right->quantifier.operand;
					if (!sets.get_set_id(set_formula, set_id)
					 || !freeable_set_size_axioms.add(c.axiom)) return false;
					unsigned int old_ref_count = c.axiom->reference_count;
					c.axiom->reference_count = 1;
					bool is_freeable = sets.is_freeable(set_id);
					c.axiom->reference_count = old_ref_count;
					if (is_freeable) {
						on_free_set(set_id, sets, std::forward<Args>(visitor)...);
						sets.free_set(set_formula, set_id);
					}
					continue;
				}
			case change_type::DEFINITION:
				remove_definition(c.axiom, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				continue;
			case change_type::FUNCTION_VALUE:
				remove_function_value(c.axiom);
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				implication_intro_nodes.remove(index_of(implication_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				negated_conjunction_nodes.remove(index_of(negated_conjunction_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					Term* term;
					if (c.intro_node.value->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.value->operands[2]->term;
					else term = c.intro_node.value->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						unsigned int index = ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.index_of(c.intro_node.value);
						ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.remove(index);
						try_free_concept_id(term->constant);
					}

					existential_intro_nodes.remove(index_of(existential_intro_nodes, c.intro_node.value));
					c.intro_node.key->reference_count--;
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				disjunction_intro_nodes.remove(index_of(disjunction_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			}
			fprintf(stderr, "theory.subtract_changes ERROR: Unrecognized `change_type`.\n");
			return false;
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

		Formula* formula = nullptr;
		Term* predicate = nullptr;
		Term* arg1 = nullptr;
		Term* arg2 = nullptr;
		switch (proof.type) {
		case ProofType::AXIOM:
			if (discharged_axioms.contains(&proof)) break;
			if (is_atomic(*proof.formula, predicate, arg1, arg2)) {
				if (reference_count != 1) return true;
				if (arg2 == NULL) {
					if (!visit_unary_atom<false>(proof.formula, std::forward<Visitor>(visitor)...)
					 || !changes.add(change_type::UNARY_ATOM, {*proof.formula, &proof})) return false;
				} else {
					if (predicate->type != TermType::CONSTANT) return false;
					if (!visit_binary_atom<false>(predicate->constant, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
					 || !changes.add(change_type::BINARY_ATOM, {{predicate->constant, arg1->constant, arg2->constant}, &proof})) return false;
				}
			} else if (proof.formula->type == FormulaType::NOT) {
				if (is_atomic(*proof.formula->unary.operand, predicate, arg1, arg2)) {
					if (reference_count != 1) return true;
					if (arg2 == NULL) {
						if (!visit_unary_atom<true>(proof.formula->unary.operand, std::forward<Visitor>(visitor)...)
						 || !changes.add(change_type::NEGATED_UNARY_ATOM, {*proof.formula->unary.operand, &proof})) return false;
					} else {
						if (predicate->type != TermType::CONSTANT) return false;
						if (!visit_binary_atom<true>(predicate->constant, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
						 || !changes.add(change_type::NEGATED_BINARY_ATOM, {{predicate->constant, arg1->constant, arg2->constant}, &proof})) return false;
					}
				} else {
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected axiom.\n");
				}
			} else if (proof.formula->type == FormulaType::FOR_ALL
					&& proof.formula->quantifier.operand->type == FormulaType::IF_THEN)
			{
				if (!visit_subset_axiom(proof, std::forward<Visitor>(visitor)...))
					return false;
				if (reference_count != 1) return true;
				return changes.add(change_type::SUBSET_AXIOM, &proof);
			} else if (proof.formula->type == FormulaType::EQUALS) {
				bool atomic = is_atomic(*proof.formula->binary.left, predicate, arg1, arg2);
				if (atomic && predicate->type == TermType::CONSTANT
				 && predicate->constant == (unsigned int) built_in_predicates::SIZE
				 && arg1->type == TermType::LAMBDA && arg2 == NULL
				 && proof.formula->binary.right->type == TermType::INTEGER)
				{
					if (reference_count != 1) return true;
					return changes.add(change_type::SET_SIZE_AXIOM, &proof);
				} else if (proof.formula->binary.left->type == TermType::CONSTANT) {
					if (reference_count != 1) return true;
					/* make sure we keep track of decrementing the reference counts of bidirectional subset edges */
					if (proof.formula->binary.right->type == TermType::LAMBDA) {
						unsigned int concept_id = proof.formula->binary.left->constant;
						for (Proof* other_definition : ground_concepts[concept_id - new_constant_offset].definitions) {
							if (proof.formula->binary.right == other_definition->formula->binary.right
							 || other_definition->formula->binary.right->type != FormulaType::LAMBDA)
								continue;
							if (!reference_counts.ensure_capacity(reference_counts.size + 2))
								return false;

							Proof* axiom = sets.template get_existing_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, proof.formula->binary.right->quantifier.operand);
							unsigned int index = reference_counts.index_of(axiom);
							if (index == reference_counts.size) {
								reference_counts.keys[index] = axiom;
								reference_counts.values[index] = axiom->reference_count;
								reference_counts.size++;
							}
							reference_counts.values[index]--;

							axiom = sets.template get_existing_subset_axiom<false>(proof.formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand);
							index = reference_counts.index_of(axiom);
							if (index == reference_counts.size) {
								reference_counts.keys[index] = axiom;
								reference_counts.values[index] = axiom->reference_count;
								reference_counts.size++;
							}
							reference_counts.values[index]--;
						}
					}
					return changes.add(change_type::DEFINITION, &proof);
				} else if (proof.formula->binary.left->type == TermType::UNARY_APPLICATION) {
					if (reference_count != 1) return true;
					return changes.add(change_type::FUNCTION_VALUE, &proof);
				} else {
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected equals axiom.\n");
				}
			} else {
				fprintf(stderr, "get_theory_changes WARNING: Found unexpected axiom.\n");
			}
			break;

		case ProofType::IMPLICATION_INTRODUCTION:
			if (reference_count != 0) return true;
			formula = implication_intro_nodes[index_of(implication_intro_nodes, &proof)].key;
			if (!changes.add(change_type::IMPLICATION_INTRO_NODE, {formula, &proof})
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
				formula = negated_conjunction_nodes[index_of(negated_conjunction_nodes, &proof)].key;
				if (!visit_negated_conjunction(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::NEGATED_CONJUNCTION_NODE, {formula, &proof})) return false;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::UNIVERSAL_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				formula = existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key;
				if (!visit_negated_universal_intro(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, {formula, &proof})) return false;
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
				fprintf(stderr, "get_theory_changes WARNING: Found unexpected proof by contradiction step.\n");
			}
			break;

		case ProofType::EXISTENTIAL_INTRODUCTION:
			if (reference_count != 0) return true;
			formula = existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key;
			if (!visit_existential_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, {formula, &proof})) return false;
			break;

		case ProofType::DISJUNCTION_INTRODUCTION:
		case ProofType::DISJUNCTION_INTRODUCTION_LEFT:
		case ProofType::DISJUNCTION_INTRODUCTION_RIGHT:
			if (reference_count != 0) return true;
			formula = disjunction_intro_nodes[index_of(disjunction_intro_nodes, &proof)].key;
			if (!visit_disjunction_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::DISJUNCTION_INTRO_NODE, {formula, &proof})) return false;
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
		Term* predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<DefinitionsAllowed, Contradiction, ResolveInconsistencies>(canonicalized, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::NOT) {
			return make_proof<!Contradiction, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->unary.operand, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
					return NULL;

				Term* constant;
				Proof* exists_not_proof = make_exists_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (exists_not_proof == NULL) return NULL;
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				 && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				{
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				Proof* assumption = ProofCalculus::new_axiom(canonicalized);
				if (assumption == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				assumption->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_universal_elim(assumption, constant), exists_not_proof), assumption);
				free(*assumption); if (assumption->reference_count == 0) free(assumption);
				if (proof == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				free(*constant); if (constant->reference_count == 0) free(constant);
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				return proof;
			}

			Formula* operand = canonicalized->quantifier.operand;
			array_map<unsigned int, unsigned int> variable_map(8);
			variable_map.keys[0] = canonicalized->quantifier.variable;
			variable_map.values[0] = 1;
			variable_map.size = 1;
			while (operand->type == FormulaType::FOR_ALL) {
				if (!variable_map.ensure_capacity(variable_map.size + 1))
					return nullptr;
				variable_map.keys[variable_map.size] = operand->quantifier.variable;
				variable_map.values[variable_map.size] = variable_map.size + 1;
				variable_map.size++;
				operand = operand->quantifier.operand;
			}

			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == FormulaType::IF_THEN) {
				Proof* new_axiom;
				Formula* left = canonicalized->quantifier.operand->binary.left;
				Formula* right = canonicalized->quantifier.operand->binary.right;

				Formula* new_left = map_variables(left, variable_map);
				if (new_left == nullptr) return nullptr;

				Formula* new_right = map_variables(right, variable_map);
				if (new_right == nullptr) {
					free(*new_left); free(new_left);
					return nullptr;
				}

				Term const* predicate; Term const* arg1; Term const* arg2;
				bool atomic = is_atomic(*left, predicate, arg1, arg2);
				if (variable_map.size == 1 && atomic && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL
				 && predicate->type == TermType::CONSTANT
				 && predicate->constant == (unsigned int) built_in_predicates::UNKNOWN)
				{
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_proof ERROR: Definitions are not allowed in this context.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					/* this is a definition of a type */
					/* check the right-hand side is a valid definition */
					if (!valid_definition(right, variable)) {
						fprintf(stderr, "theory.make_proof ERROR: This is not a valid type definition.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					new_constant = get_free_concept_id();

					/* TODO: definitions of new concepts should be biconditionals */
					if (!try_init_concept(new_constant)) return NULL;

					new_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(new_left, new_right, variable_map.size, std::forward<Args>(args)...);
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				} else {
					/* make sure there are no unknown predicates */
					if (contains_constant(*left, (unsigned int) built_in_predicates::UNKNOWN)
					 || contains_constant(*right, (unsigned int) built_in_predicates::UNKNOWN)) {
						fprintf(stderr, "theory.make_proof ERROR: Universally-quantified statement has unknown constants.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					/* this is a formula of form `![x]:(t(x) => f(x))` */
					new_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(new_left, new_right, variable_map.size, std::forward<Args>(args)...);
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				}
				return new_axiom;
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
				if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

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
							free_proof(operand, std::forward<Args>(args)...); return NULL;
						}
						axiom->reference_count++;
						Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
								ProofCalculus::new_conjunction_elim(axiom, make_array_view(&index_minus_one, 1)), operand), axiom);
						free(*axiom); if (axiom->reference_count == 0) free(axiom);
						if (proof == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
							free_proof(operand, std::forward<Args>(args)...); return NULL;
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
					for (unsigned int j = i; j > 0; j--)
						/* undo the changes made by the recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
					free(operands); return NULL;
				}
			}
			Proof* conjunction = ProofCalculus::new_conjunction_intro(make_array_view(operands, canonicalized->array.length));
			if (conjunction == NULL) {
				for (unsigned int j = canonicalized->array.length; j > 0; j--)
					/* undo the changes made by the recursive calls to `make_proof` */
					free_proof(operands[j - 1], std::forward<Args>(args)...);
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
						for (unsigned int j = i; j > 0; j--)
							/* undo the changes made by recursive calls to `make_proof` */
							free_proof(operands[j - 1], std::forward<Args>(args)...);
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
							for (unsigned int j = i; j > 0; j--)
								/* undo the changes made by recursive calls to `make_proof` */
								free_proof(operands[j - 1], std::forward<Args>(args)...);
							free(operands); return NULL;
						}
						absurdity->reference_count++;
						free(*operands[i]); if (operands[i]->reference_count == 0) free(operands[i]);
						operands[i] = absurdity;
					}
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
					free(operands); return NULL;
				}
				axiom->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(
						ProofCalculus::new_disjunction_elim(axiom, make_array_view(operands, canonicalized->array.length)), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
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
			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

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
							free_proof(operand, std::forward<Args>(args)...); return NULL;
						}
						for (Formula* other_disjunct : other_disjuncts)
							other_disjunct->reference_count++;
					}
					Proof* proof = ProofCalculus::new_disjunction_intro(operand, other_disjunction, index - 1);
					free(*other_disjunction); if (other_disjunction->reference_count == 0) free(other_disjunction);
					if (proof == NULL) {
						/* undo changes made by the recursive call to `make_proof` */
						free_proof(operand, std::forward<Args>(args)...); return NULL;
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
					free_proof(left, std::forward<Args>(args)...); return NULL;
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					free_proof(right, std::forward<Args>(args)...);
					free_proof(left, std::forward<Args>(args)...);
					return NULL;
				}
				axiom->reference_count++;

				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_implication_elim(axiom, left), right), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					free_proof(right, std::forward<Args>(args)...);
					free_proof(left, std::forward<Args>(args)...);
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
			if (sample_uniform(2) == 1)
				swap(indices[0], indices[1]);
			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;
			for (unsigned int i = 0; i < 2; i++) {
				if (indices[i] == 1) {
					Proof* left = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, new_constant, std::forward<Args>(args)...);
					if (left == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* axiom = ProofCalculus::new_axiom(canonicalized->binary.left);
					if (axiom == NULL) {
						free_proof(left, std::forward<Args>(args)...);
						return NULL;
					}
					axiom->reference_count++;

					Proof* proof = ProofCalculus::new_implication_intro(
							ProofCalculus::new_falsity_elim(
								ProofCalculus::new_negation_elim(left, axiom),
								canonicalized->binary.right),
							axiom);
					free(*axiom); if (axiom->reference_count == 0) free(axiom);
					if (proof == NULL) {
						free_proof(left, std::forward<Args>(args)...);
						return NULL;
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
						free_proof(right, std::forward<Args>(args)...);
						return NULL;
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
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(set_formula, 0, std::forward<Args>(args)...);
				if (set_size_axiom == NULL) {
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					return NULL;
				}

				Formula* beta_left = Formula::new_not(Formula::new_exists(1, Formula::new_apply(lambda_formula, &Variables<1>::value)));
				Formula* beta_right = Formula::new_not(Formula::new_exists(1, set_formula));
				if (beta_left == NULL || beta_right == NULL) {
					if (beta_left != NULL) { free(*beta_left); if (beta_left->reference_count == 0) free(beta_left); }
					if (beta_right != NULL) { free(*beta_right); if (beta_right->reference_count == 0) free(beta_right); }
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					sets.try_free_set(set_formula, std::forward<Args>(args)...); return NULL;
				}
				Variables<1>::value.reference_count++;
				set_formula->reference_count++;
				Proof* proof = ProofCalculus::new_equality_elim(
						ProofCalculus::new_beta(beta_left, beta_right),
						ProofCalculus::new_equality_elim(ProofCalculus::new_universal_elim(empty_set_axiom, lambda_formula), set_size_axiom, make_repeated_array_view(0u, 1)),
						make_repeated_array_view(0u, 1));
				free(*beta_left); if (beta_left->reference_count == 0) free(beta_left);
				free(*beta_right); if (beta_right->reference_count == 0) free(beta_right);
				free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
				if (proof == NULL) { sets.try_free_set(set_formula, std::forward<Args>(args)...); return NULL; }
				lambda_formula->reference_count++;
				proof->reference_count++;
				return proof;
			} else {
				Term* variable = Formula::new_variable(canonicalized->quantifier.variable);
				if (variable == NULL) return NULL;

				Term* constant;
				Proof* operand = make_exists_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (operand == NULL) {
					free(*variable); if (variable->reference_count == 0) free(variable);
					return NULL;
				}

				array<unsigned int> indices(8);
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1)
				 || (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				  && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				 || !compute_indices(*canonicalized->quantifier.operand, *variable, indices))
				{
					free(*variable); if (variable->reference_count == 0) free(variable);
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_proof(operand, std::forward<Args>(args)...);
					return NULL;
				}
				free(*variable); if (variable->reference_count == 0) free(variable);
				Proof* proof = ProofCalculus::new_existential_intro(operand, indices.data, indices.length, constant);
				if (proof == NULL) {
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_proof(operand, std::forward<Args>(args)...);
					return NULL;
				}
				free(*operand); if (operand->reference_count == 0) free(operand);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				free(*constant); if (constant->reference_count == 0) free(constant);
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				return proof;
			}

		} else if (canonicalized->type == FormulaType::EQUALS) {
			Formula* left = canonicalized->binary.left;
			Formula* right = canonicalized->binary.right;
			Term* arg1; Term* arg2;
			bool atomic = is_atomic(*left, predicate, arg1, arg2);
			if (atomic && predicate->type == TermType::CONSTANT
			 && predicate->constant == (unsigned int) built_in_predicates::SIZE
			 && arg1->type == TermType::LAMBDA && arg2 == NULL
			 && right->type == TermType::INTEGER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(arg1->quantifier.operand, right->integer, std::forward<Args>(args)...);
				if (set_size_axiom == nullptr) return nullptr;
				set_size_axiom->reference_count++;
				return set_size_axiom;

			} else if (atomic && predicate->type == TermType::CONSTANT
					&& predicate->constant == (unsigned int) built_in_predicates::SIZE
			 		&& arg1->type == TermType::CONSTANT && arg2 == NULL
			 		&& right->type == TermType::INTEGER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				/* this is a statement on the size of a set */
				Proof* definition = ground_concepts[arg1->constant - new_constant_offset].definitions[0];
#if !defined(NDEBUG)
				if (ground_concepts[arg1->constant - new_constant_offset].definitions.length == 0
				 || definition->formula->binary.right->type != FormulaType::LAMBDA
				 || definition->formula->binary.right->quantifier.operand->type != FormulaType::UNARY_APPLICATION
				 || definition->formula->binary.right->quantifier.operand->binary.left->type != TermType::CONSTANT
				 || definition->formula->binary.right->quantifier.operand->binary.left->constant != arg1->constant
				 || definition->formula->binary.right->quantifier.operand->binary.right->type != TermType::VARIABLE
				 || definition->formula->binary.right->quantifier.operand->binary.right->variable != definition->formula->binary.right->quantifier.variable)
					fprintf(stderr, "theory.make_proof ERROR: Expected a set definition of the form c=λx.c(x).\n");
#endif
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(definition->formula->binary.right->quantifier.operand, right->integer, std::forward<Args>(args)...);
				Proof* proof = ProofCalculus::new_equality_elim(definition, set_size_axiom, make_repeated_array_view(3u, 1));
				if (proof == NULL)
					return NULL;
				proof->reference_count++;
				return proof;

			} else if (left == right || *left == *right) {
				Proof* proof = ProofCalculus::new_beta(left, left);
				if (proof == NULL) return NULL;
				proof->reference_count++;
				return proof;

			} else if (left->type == TermType::CONSTANT && right->type == TermType::CONSTANT) {
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
			} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT) {
				if (right->type == TermType::CONSTANT && right->constant != (unsigned int) built_in_predicates::UNKNOWN)
					swap(left, right);

				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				} else {
					if (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT
					 && right->binary.right->type == TermType::CONSTANT
					 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
					{
						/* we disallow objects to be arguments of themselves */
						if (right->binary.right->constant == left->constant)
							return NULL;

						/* check if the other object has this function value */
						if (right->binary.right->constant >= new_constant_offset
						 && ground_concepts[right->binary.right->constant - new_constant_offset].function_values.contains((unsigned int) right->binary.left->constant))
							return NULL;
					}

					/* check that this constant could be a set */
					if (right->type == FormulaType::LAMBDA && is_provably_not_a_set(left->constant))
						return NULL;

					unsigned int min_variable = UINT_MAX;
					min_bound_variable(*right, min_variable);
					Formula* new_right;
					if (min_variable != UINT_MAX) {
						new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
					} else {
						new_right = right;
						right->reference_count++;
					}

					/* check if anything else has this definition */
					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						if (ground_concepts[i].types.keys == NULL || i == left->constant - new_constant_offset) continue;
						for (Proof* proof : ground_concepts[i].definitions) {
							if (proof->formula->binary.right == new_right || *proof->formula->binary.right == *new_right) {
								/* we found a different concept with this definition */
								free(*new_right); if (new_right->reference_count == 0) free(new_right);
								return NULL;
							}
						}
					}

					Formula* new_formula = Formula::new_equals(left, new_right);
					if (new_formula == NULL) {
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						return NULL;
					}
					left->reference_count++;

					Proof* new_proof = ProofCalculus::new_axiom(new_formula);
					free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
					if (new_proof == NULL)
						return NULL;
					new_proof->reference_count++;

					Proof* definition = add_definition<ResolveInconsistencies>(new_proof, std::forward<Args>(args)...);
					if (definition != new_proof) {
						free(*new_proof); free(new_proof);
					}
					return definition;
				}
			} else if ((left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT && left->binary.right->type == TermType::CONSTANT)
					|| (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT && right->binary.right->type == TermType::CONSTANT))
			{
				bool swap_order = (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT && right->binary.right->type == TermType::CONSTANT);
				if (swap_order) swap(left, right);

				if (is_provably_a_set(left->binary.right->constant))
					return NULL;

				/* check if anything else has this as a definition */
				for (unsigned int i = 0; i < ground_concept_capacity; i++) {
					if (ground_concepts[i].types.keys == NULL) continue;
					for (Proof* proof : ground_concepts[i].definitions) {
						if (proof->formula->binary.right == left || *proof->formula->binary.right == *left)
							return NULL;
					}
				}

				/* this is a function value definition */
				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				Formula* new_formula = Formula::new_equals(left, new_right);
				if (new_formula == NULL) {
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					return NULL;
				}
				left->reference_count++;

				Proof* new_proof = ProofCalculus::new_axiom(new_formula);
				free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
				if (new_proof == NULL)
					return NULL;
				new_proof->reference_count++;

				Proof* function_value = add_function_value<ResolveInconsistencies>(new_proof);
				if (function_value != new_proof) {
					free(*new_proof); free(new_proof);
				} if (function_value == nullptr) {
					return nullptr;
				}

				if (swap_order) {
					Proof* swapped_proof = ProofCalculus::new_equality_elim(
							function_value, ProofCalculus::new_beta(right, right), make_repeated_array_view(2u, 1));
					if (swapped_proof == NULL) {
						free_proof(function_value);
						return NULL;
					}
					free(*function_value);
					swapped_proof->reference_count++;
					return swapped_proof;
				} else {
					return function_value;
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
	template<bool ResolveInconsistencies, typename... Args>
	Proof* add_definition(Proof* definition, Args&&... args)
	{
		Formula* constant = definition->formula->binary.left;
		Formula* new_definition = definition->formula->binary.right;

		if (!try_init_concept(constant->constant))
			return NULL;

		/* check if the axiom already exists */
		for (Proof* definition : ground_concepts[constant->constant - new_constant_offset].definitions) {
			if ((definition->formula->binary.left == constant || *definition->formula->binary.left == *constant)
			 && (definition->formula->binary.right == new_definition || *definition->formula->binary.right == *new_definition))
			{
				definition->reference_count++;
				return definition;
			}
		}

		if (new_definition->type == FormulaType::LAMBDA) {
			/* check if this constant defines any other sets, and indicate to the set reasoning module that they are the same set */
			bool contains;
			unsigned int set_size = UINT_MAX;
			unsigned int set_id = sets.set_ids.get(*new_definition->quantifier.operand, contains);
			if (contains) {
				set_size = sets.sets[set_id].set_size;
			} else if (ground_concepts[constant->constant - new_constant_offset].definitions.length != 0) {
				for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
					Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
					if (definition->formula->binary.right->type != FormulaType::LAMBDA)
						continue;
					set_id = sets.set_ids.get(*definition->formula->binary.right->quantifier.operand, contains);
					if (contains) {
						set_size = sets.sets[set_id].set_size;
						break;
					}
				}
			}

			for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					required_set_size set_size_computer(set_size);
					Proof* first_subset_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, set_size_computer, std::forward<Args>(args)...);
					if (first_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Proof* axiom = sets.template get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

							axiom = sets.template get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
						}
						return NULL;
					}
					first_subset_axiom->reference_count++;

					Proof* second_subset_axiom = sets.template get_subset_axiom<ResolveInconsistencies>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
					if (second_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						free(*first_subset_axiom);
						if (first_subset_axiom->reference_count == 1)
							sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Proof* axiom = sets.template get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

							axiom = sets.template get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
						}
						return NULL;
					}
					second_subset_axiom->reference_count++;
				}
			}
		}

		if (!ground_concepts[constant->constant - new_constant_offset].definitions.add(definition)) {
			/* undo the changes we've made so far */
			on_subtract_changes(std::forward<Args>(args)...);
			for (unsigned int j = 0; j < ground_concepts[constant->constant - new_constant_offset].definitions.length; j++) {
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
				Proof* axiom = sets.template get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == 1)
					sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

				axiom = sets.template get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == 1)
					sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
			}
			return NULL;
		}
		definition->reference_count++;
		return definition;
	}

	template<typename... Args>
	void remove_definition(Proof* definition, Args&&... args) {
		unsigned int concept_id = definition->formula->binary.left->constant;

		/* remove subset edges from other set definitions for `concept_id` */
		unsigned int index = ground_concepts[concept_id - new_constant_offset].definitions.length;
		for (unsigned int i = ground_concepts[concept_id - new_constant_offset].definitions.length; i > 0; i--) {
			Proof* other_definition = ground_concepts[concept_id - new_constant_offset].definitions[i - 1];
			if (definition->formula->binary.right == other_definition->formula->binary.right) {
				index = i - 1;
				continue;
			}
			if (other_definition->formula->binary.right->type != FormulaType::LAMBDA)
				continue;
			if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
				Proof* axiom = sets.template get_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					sets.free_subset_axiom(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

				axiom = sets.template get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);
			}
		}

#if !defined(NDEBUG)
		if (index == ground_concepts[concept_id - new_constant_offset].definitions.length)
			fprintf(stderr, "remove_definition WARNING: Unable to find definition to remove.\n");
#endif

		ground_concepts[concept_id - new_constant_offset].definitions.remove(index);
		try_free_concept_id(concept_id);
		free(*definition); if (definition->reference_count == 0) free(definition);
	}

	template<bool ResolveInconsistencies>
	Proof* add_function_value(Proof* function_value_axiom)
	{
		Formula* constant = function_value_axiom->formula->binary.left->binary.right;
		if (!try_init_concept(constant->constant))
			return NULL;

		/* first check if this function value is already defined */
		array_map<unsigned int, Proof*>& function_values = ground_concepts[constant->constant - new_constant_offset].function_values;
		if (!function_values.ensure_capacity(function_values.size + 1))
			return NULL;
		unsigned int index = function_values.index_of(function_value_axiom->formula->binary.left->binary.left->constant);
		if (index < function_values.size) {
			if ((function_value_axiom->formula->binary.right == function_values.values[index]->formula->binary.right)
			 || (*function_value_axiom->formula->binary.right == *function_values.values[index]->formula->binary.right))
			{
				/* this definition is already an axiom, so return it */
				function_values.values[index]->reference_count++;
				return function_values.values[index];
			} else {
				return NULL;
			}
		}

		function_values.keys[index] = function_value_axiom->formula->binary.left->binary.left->constant;
		function_values.values[index] = function_value_axiom;
		function_values.size++;
		function_value_axiom->reference_count++;
		return function_value_axiom;
	}

	void remove_function_value(Proof* function_value_axiom) {
		unsigned int concept_id = function_value_axiom->formula->binary.left->binary.right->constant;

		array_map<unsigned int, Proof*>& function_values = ground_concepts[concept_id - new_constant_offset].function_values;
		unsigned int index = function_values.index_of(function_value_axiom->formula->binary.left->binary.left->constant);
#if !defined(NDEBUG)
		if (index == function_values.size)
			fprintf(stderr, "remove_function_value WARNING: Unable to find axiom to remove.\n");
#endif

		function_values.remove_at(index);
		try_free_concept_id(concept_id);
		free(*function_value_axiom); if (function_value_axiom->reference_count == 0) free(function_value_axiom);
	}

	/* NOTE: this function finds a constant that proves `quantified`, and not ?[x]:`quantified` */
	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_exists_proof(Formula* quantified, unsigned int variable, Term*& constant, unsigned int& new_constant, Args&&... args)
	{
		Formula* var = Formula::new_variable(variable);
		if (var == NULL) return NULL;

		array<instance> constants(ground_concept_capacity + 1);
		hash_set<int64_t> integers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != NULL) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length++].constant = new_constant_offset + i;

				for (const auto& entry : ground_concepts[i].function_values) {
					Term* constant = entry.value->formula->binary.right;
					if (constant->type == TermType::INTEGER) {
						if (!integers.add(constant->integer)) return NULL;
					} else if (constant->type == TermType::STRING) {
						bool contains = false;
						for (const string* str : strings)
							if (str == &constant->str || *str == constant->str) { contains = true; break; }
						if (!contains && !strings.add(&constant->str)) return NULL;
					}
				}
			}
		}
		constants[constants.length++].type = instance_type::ANY;
		if (!constants.ensure_capacity(constants.length + integers.size + strings.length))
			return NULL;
		for (int64_t integer : integers) {
			constants[constants.length].type = instance_type::INTEGER;
			constants[constants.length++].integer = integer;
		} for (string* str : strings) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length++].str = str;
		}
		shuffle(constants);

		unsigned int original_constant_count = constants.length;
		if (!filter_constants(*this, quantified, variable, constants, std::forward<Args>(args)...)) {
			free(*var); if (var->reference_count == 0) free(var);
			return NULL;
		}

		for (const instance& id : constants) {
			if (id.type == instance_type::ANY || (id.type == instance_type::CONSTANT && ground_concepts[id.constant - new_constant_offset].types.keys == nullptr)) {
				unsigned int constant_id;
				if (id.type == instance_type::ANY)
					constant_id = get_free_concept_id();
				else constant_id = id.constant;

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
			} else {
				if (id.type == instance_type::CONSTANT) {
					constant = Formula::new_constant(id.constant);
				} else if (id.type == instance_type::INTEGER) {
					constant = Formula::new_int(id.integer);
				} else if (id.type == instance_type::STRING) {
					constant = Formula::new_string(*id.str);
				}
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
			}
		}

		finished_constants(quantified, original_constant_count, std::forward<Args>(args)...);
		free(*var); if (var->reference_count == 0) free(var);
		return NULL;
	}

	template<bool DefinitionsAllowed, bool Negated, bool ResolveInconsistencies, typename... Args>
	Proof* make_atom_proof(
			Term* atom,
			unsigned int& new_constant,
			Args&&... args)
	{
		if (atom->type == TermType::UNARY_APPLICATION && atom->binary.right->type == TermType::CONSTANT)
		{
			/* this is a unary formula */
			Term* arg1 = atom->binary.right;
			if (atom->binary.left->type != TermType::CONSTANT || atom->binary.left->constant != (unsigned int) built_in_predicates::UNKNOWN) {
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					/* make sure that `predicate` could be a set, and that `arg` could be not a set */
					if (!Negated && atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant))
						return NULL;

					/* this is a definition of an object */
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_atom_proof ERROR: Definitions are not allowed in this context.\n");
						return NULL;
					}
					new_constant = get_free_concept_id();

					Formula* formula = Negated ?
							Formula::new_not(Term::new_apply(atom->binary.left, Term::new_constant(new_constant))) :
							Term::new_apply(atom->binary.left, Term::new_constant(new_constant));
					if (formula == NULL) return NULL;
					atom->binary.left->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;

					if (!try_init_concept(new_constant)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
					} if (!add_unary_atom<Negated, ResolveInconsistencies>(Negated ? *formula->unary.operand : *formula, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom);
						free_concept_id(new_constant); return NULL;
					}
					return new_axiom;
				} else {
					/* this is a formula of form `t(c)` */

					/* we do not allow sets to contain themselves */
					if (!Negated && *atom->binary.left == *arg1)
						return NULL;

					/* make sure that `predicate` could be a set, and that `arg` could be not a set */
					if (!Negated && ((atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant)) || is_provably_a_set(arg1->constant)))
						return NULL;

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

						Term* expected_constant_eliminator = nullptr;
						if (set_formula->binary.left->type == TermType::VARIABLE) {
							expected_constant_eliminator = atom->binary.left;
						} else if (set_formula->binary.left->type != TermType::CONSTANT || *set_formula->binary.left != *atom->binary.left) {
							continue;
						}

						if (set_formula->binary.right->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg1;
							else if (*expected_constant_eliminator != *arg1)
								continue;
						} else if (*set_formula->binary.right != *arg1) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							for (Proof* child : entry.value->children) {
								if (child->type == ProofType::UNIVERSAL_ELIMINATION
								 && child->operands[1]->type == ProofType::TERM_PARAMETER
								 && *child->operands[1]->term == *expected_constant_eliminator)
								{
									for (Proof* grandchild : child->children) {
										if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
											/* we found a proof of the atomic formula */
											grandchild->reference_count++;
											return grandchild;
										}
									}
								}
							}
						}
					}

					/* there is no extensional edge that proves this atomic formula,
					  so check if the axiom already exists, and if not, create a new axiom */
					Term* lifted_atom = Term::new_apply(atom->binary.left, &Variables<1>::value);
					if (lifted_atom == nullptr) return nullptr;
					atom->binary.left->reference_count++;
					Variables<1>::value.reference_count++;

					bool contains;
					concept<ProofCalculus>& c = ground_concepts[arg1->constant - new_constant_offset];
					Proof* axiom = (Negated ?
							c.negated_types.get(*lifted_atom, contains) :
							c.types.get(*lifted_atom, contains));
					if (contains) {
						axiom->reference_count++;
						free(*lifted_atom); free(lifted_atom);
						return axiom;
					}
					free(*lifted_atom); free(lifted_atom);
					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) return NULL;
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated, ResolveInconsistencies>(*atom, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_atom_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else if (atom->type == TermType::BINARY_APPLICATION
				&& atom->ternary.second->type == TermType::CONSTANT
				&& atom->ternary.third->type == TermType::CONSTANT)
		{
			/* this is a binary formula */
			Term* arg1 = atom->ternary.second;
			Term* arg2 = atom->ternary.third;
			if (atom->ternary.first->type != TermType::CONSTANT || atom->ternary.first->constant != (unsigned int) built_in_predicates::UNKNOWN)
			{
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN
				 || arg2->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
					return NULL;
				} else {
					/* this is a formula of form `r(c_1,c_2)` */
					if (atom->ternary.first->type != TermType::CONSTANT) {
						fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
						return nullptr;
					}

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

						Term* expected_constant_eliminator = nullptr;
						if (set_formula->ternary.first->type == TermType::VARIABLE) {
							expected_constant_eliminator = atom->ternary.first;
						} else if (set_formula->ternary.first->type != TermType::CONSTANT || *set_formula->ternary.first != *atom->ternary.first) {
							continue;
						}

						if (set_formula->ternary.second->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg1;
							else if (*expected_constant_eliminator != *arg1)
								continue;
						} else if (*set_formula->ternary.second != *arg1) {
							continue;
						}

						if (set_formula->ternary.third->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg2;
							else if (*expected_constant_eliminator != *arg2)
								continue;
						} else if (*set_formula->ternary.third != *arg2) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							for (Proof* child : entry.value->children) {
								if (child->type == ProofType::UNIVERSAL_ELIMINATION
								 && child->operands[1]->type == ProofType::TERM_PARAMETER
								 && *child->operands[1]->term == *expected_constant_eliminator)
								{
									for (Proof* grandchild : child->children) {
										if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
											/* we found a proof of the atomic formula */
											grandchild->reference_count++;
											return grandchild;
										}
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
							c.negated_relations.get(relation(atom->ternary.first->constant, 0, arg2->constant), contains) :
							c.relations.get(relation(atom->ternary.first->constant, 0, arg2->constant), contains));
					if (contains) {
						axiom->reference_count++;
						return axiom;
					}
					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) return NULL;
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated, ResolveInconsistencies>({atom->ternary.first->constant, arg1->constant, arg2->constant}, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_atom_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else {
			fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
			return NULL;
		}
	}

	template<typename... Args>
	void free_proof(Proof* proof, Args&&... args) {
		theory::changes changes;
		if (!get_theory_changes(*proof, changes, std::forward<Args>(args)...)) return;
		subtract_changes(changes, std::forward<Args>(args)...);
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
				Term* atom = Term::new_apply(Term::new_constant(predicate), &Variables<1>::value);
				if (atom == nullptr) return false;
				Variables<1>::value.reference_count++;
				proof = (Negated ? c.negated_types.get(*atom, contains) : c.types.get(*atom, contains));
				free(*atom); free(atom);
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
constexpr bool filter_operands(const Formula* formula, const array<unsigned int>& constants) { return true; }

inline unsigned int index_of_any(const array<instance>& constants) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::ANY)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_constant(const array<instance>& constants, unsigned int constant) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::CONSTANT && constants[i].constant == constant)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_integer(const array<instance>& constants, int64_t integer) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::INTEGER && constants[i].integer == integer)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_string(const array<instance>& constants, string* str) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::STRING && (constants[i].str == str || *constants[i].str == *str))
			return i;
	}
	return constants.length;
}

inline bool contains_any(const array<instance>& constants) {
	return index_of_any(constants) < constants.length;
}

inline bool contains_constant(const array<instance>& constants, unsigned int constant) {
	return index_of_constant(constants, constant) < constants.length;
}

inline bool contains_integer(const array<instance>& constants, int64_t integer) {
	return index_of_integer(constants, integer) < constants.length;
}

inline bool contains_string(const array<instance>& constants, string* str) {
	return index_of_string(constants, str) < constants.length;
}

template<typename ProofCalculus, typename Canonicalizer>
bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Term Term;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	if (formula->type == FormulaType::EQUALS) {
		Term* left = formula->binary.left;
		Term* right = formula->binary.right;
		if (left == right || *left == *right) {
			/* no-op */
		} else if ((left->type == TermType::VARIABLE && left->variable == variable)
				|| (right->type == TermType::VARIABLE && right->variable == variable))
		{
			if (right->type == TermType::VARIABLE && right->variable == variable)
				swap(left, right);
			if (right->type == TermType::CONSTANT) {
				unsigned int constant_id = right->constant;
				if (!contains_constant(constants, constant_id) && (!contains_any(constants) || T.ground_concepts[constant_id - T.new_constant_offset].types.keys != nullptr))
					return false;
				constants[0].type = instance_type::CONSTANT;
				constants[0].constant = right->constant;
				constants.length = 1;
			} else if (right->type == TermType::INTEGER) {
				int64_t integer = right->integer;
				if (!contains_integer(constants, integer) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::INTEGER;
				constants[0].integer = right->integer;
				constants.length = 1;
			} else if (right->type == TermType::STRING) {
				string* str = &right->str;
				if (!contains_string(constants, str) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::STRING;
				constants[0].str = &right->str;
				constants.length = 1;
			}
			/* TODO: handle the case where `right` is an integer or a string */
			else if (right->type == TermType::VARIABLE) {
				/* no-op */
			} else {
				if (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT
				 && right->binary.right->type == TermType::CONSTANT
				 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
				{
					/* we disallow objects to be arguments of themselves */
					unsigned int index = index_of_constant(constants, right->binary.right->constant);
					if (index < constants.length)
						constants.remove(index);

					/* check if the other object has this function value */
					if (right->binary.right->constant >= T.new_constant_offset
					 && T.ground_concepts[right->binary.right->constant - T.new_constant_offset].function_values.contains((unsigned int) right->binary.left->constant))
						return false;
				}

				for (unsigned int i = 0; right->type == FormulaType::LAMBDA && i < constants.length; i++) {
					/* check that this constant could be a set */
					if (constants[i].type == instance_type::ANY) {
						continue;
					} else if (constants[i].type != instance_type::CONSTANT) {
						constants.remove(i--);
					} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
						continue;
					} else if (T.is_provably_not_a_set(constants[i].constant)) {
						constants.remove(i--);
					}
				}
				if (constants.length == 0) return false;

				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				/* check if anything else has this definition */
				for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
					if (T.ground_concepts[i].types.keys == nullptr) continue;
					for (Proof* proof : T.ground_concepts[i].definitions) {
						if (proof->formula->binary.right == new_right || *proof->formula->binary.right == *new_right) {
							/* we found a different concept with this definition */
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							unsigned int constant_id = T.new_constant_offset + i;
							if (!contains_constant(constants, constant_id) && (!contains_any(constants) || T.ground_concepts[i].types.keys != nullptr))
								return false;
							constants[0].type = instance_type::CONSTANT;
							constants[0].constant = constant_id;
							constants.length = 1;
							return true;
						}
					}
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
			}
		} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT) {
			if (right->type == TermType::CONSTANT && right->constant != (unsigned int) built_in_predicates::UNKNOWN)
				swap(left, right);

			/* check that this constant could be a set */
			if (right->type == FormulaType::LAMBDA && T.is_provably_not_a_set(left->constant))
				return false;

			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY || constants[i].type != instance_type::CONSTANT) {
					continue;
				} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
					continue;
				}

				/* substitute `variable` in the right-hand side with `constants[i]` */
				Term* src_var = Term::new_variable(variable);
				if (src_var == nullptr) return false;
				Term* dst_term = Term::new_constant(constants[i].constant);
				if (dst_term == nullptr) {
					free(*src_var); free(src_var);
					return false;
				}
				Formula* new_right = substitute(right, src_var, dst_term);
				free(*src_var); free(src_var);
				free(*dst_term); if (dst_term->reference_count == 0) free(dst_term);
				if (new_right == nullptr)
					return false;

				if (*left == *new_right) {
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					return true;
				} else if (new_right->type != TermType::CONSTANT) {
					if (new_right->type == TermType::UNARY_APPLICATION && new_right->binary.left->type == TermType::CONSTANT
					 && new_right->binary.right->type == TermType::CONSTANT
					 && (new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					  || new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					  || new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
					{
						/* we disallow objects to be arguments of themselves */
						if (new_right->binary.right->constant == left->constant) {
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							constants.remove(i--); continue;
						}

						/* check if the other object has this function value */
						if (new_right->binary.right->constant >= T.new_constant_offset
						 && T.ground_concepts[new_right->binary.right->constant - T.new_constant_offset].function_values.contains((unsigned int) new_right->binary.left->constant))
						{
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							constants.remove(i--); continue;
						}
					}

					unsigned int min_variable = UINT_MAX;
					min_bound_variable(*new_right, min_variable);
					Formula* new_new_right;
					if (min_variable != UINT_MAX) {
						new_new_right = shift_bound_variables(new_right, -((int) (min_variable - 1)));
					} else {
						new_new_right = new_right;
						new_right->reference_count++;
					}

					/* check if anything else has this definition */
					bool found_conflicting_definition = false;
					for (unsigned int j = 0; j < T.ground_concept_capacity && !found_conflicting_definition; j++) {
						if (T.ground_concepts[j].types.keys == NULL || j == left->constant - T.new_constant_offset) continue;
						for (Proof* proof : T.ground_concepts[j].definitions) {
							if (proof->formula->binary.right == new_new_right || *proof->formula->binary.right == *new_new_right) {
								/* we found a different concept with this definition */
								found_conflicting_definition = true;
								break;
							}
						}
					}
					free(*new_new_right); if (new_new_right->reference_count == 0) free(new_new_right);
					if (found_conflicting_definition)
						constants.remove(i--);
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
			}
		} else if (left->type == TermType::UNARY_APPLICATION
				&& ((left->binary.left->type == TermType::CONSTANT && left->binary.right->type == TermType::VARIABLE && left->binary.right->variable == variable)
				 || (left->binary.right->type == TermType::CONSTANT && left->binary.left->type == TermType::VARIABLE && left->binary.left->variable == variable)))
		{
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT) {
					constants.remove(i--); continue;
				} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
					continue;
				}

				/* substitute `variable` in the left-hand side with `constants[i]` */
				Term* src_var = Term::new_variable(variable);
				if (src_var == nullptr) return false;
				Term* dst_term = Term::new_constant(constants[i].constant);
				if (dst_term == nullptr) {
					free(*src_var); free(src_var);
					return false;
				}
				Formula* new_left = substitute(left, src_var, dst_term);
				free(*src_var); free(src_var);
				free(*dst_term); if (dst_term->reference_count == 0) free(dst_term);
				if (new_left == nullptr)
					return false;

				/* check if anything else has this as a definition */
				bool found_conflicting_definition = false;
				for (unsigned int i = 0; i < T.ground_concept_capacity && !found_conflicting_definition; i++) {
					if (T.ground_concepts[i].types.keys == NULL) continue;
					for (Proof* proof : T.ground_concepts[i].definitions) {
						if (*proof->formula->binary.right == *new_left) {
							found_conflicting_definition = true;
							break;
						}
					}
				}
				if (found_conflicting_definition) {
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					constants.remove(i--); continue;
				}

				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				/* first check if this function value is already defined */
				const array_map<unsigned int, Proof*>& function_values = T.ground_concepts[new_left->binary.right->constant - T.new_constant_offset].function_values;
				unsigned int index = function_values.index_of(new_left->binary.left->constant);
				if (index < function_values.size) {
					if ((new_right != function_values.values[index]->formula->binary.right)
					 && (*new_right != *function_values.values[index]->formula->binary.right))
					{
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						free(*new_left); if (new_left->reference_count == 0) free(new_left);
						constants.remove(i--); continue;
					}
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
			}
		}
	} else if (formula->type == FormulaType::UNARY_APPLICATION) {
		Term* left = formula->binary.left;
		Term* right = formula->binary.right;
		if (left->type == TermType::VARIABLE && left->variable == variable) {
			if (right->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, right->constant);
				if (index < constants.length)
					constants.remove(index);
			}
			/* make sure x could be a set in `x(y)` */
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT || T.is_provably_not_a_set(constants[i].constant)) {
					constants.remove(i);
					i--;
				}
			}
		} if (right->type == TermType::VARIABLE && right->variable == variable) {
			if (left->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, left->constant);
				if (index < constants.length)
					constants.remove(index);
			}
			/* make sure y could not be a set in `x(y)` */
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT || T.is_provably_a_set(constants[i].constant)) {
					constants.remove(i);
					i--;
				}
			}
		}
	} else if (formula->type == FormulaType::AND) {
		for (unsigned int i = 0; i < formula->array.length; i++)
			if (!filter_constants(T, formula->array.operands[i], variable, constants)) return false;
	} else if (formula->type == FormulaType::EXISTS || formula->type == FormulaType::FOR_ALL) {
		return filter_constants(T, formula->quantifier.operand, variable, constants);
	}
	return (constants.length > 0);
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant) {
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index) {
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

template<typename ProofCalculus, typename Canonicalizer>
typename ProofCalculus::Language* make_lifted_conjunction(unsigned int concept_id, const theory<ProofCalculus, Canonicalizer>& T)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;

	const concept<ProofCalculus>& c = T.ground_concepts.get(concept_id);
	array<Formula*> conjuncts(16);
	for (const auto& entry : c.types) {
		Formula* new_conjunct = Formula::new_apply(&entry.key, &Term::template variables<1>::value);
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		entry.key.reference_count++;
		Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.negated_types) {
		Formula* new_conjunct = Formula::new_not(Formula::new_apply(&entry.key, &Term::template variables<1>::value));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		entry.key.reference_count++;
		Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.relations) {
		Formula* new_conjunct = Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg2));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		if (entry.key.arg1 == 0) Term::template variables<1>::value.reference_count++;
		if (entry.key.arg2 == 0) Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.negated_relations) {
		Formula* new_conjunct = Formula::new_not(Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg2)));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		if (entry.key.arg1 == 0) Term::template variables<1>::value.reference_count++;
		if (entry.key.arg2 == 0) Term::template variables<1>::value.reference_count++;
	}
	Formula* conjunction = Formula::new_and(conjuncts);
	return conjunction;
}

template<typename Proof>
struct theory_sample {
	Proof** proofs;
	unsigned int count;
	double log_probability;
#if !defined(NDEBUG)
	/* a unique identifier for debugging */
	unsigned int id;
#endif

	static inline unsigned int hash(const theory_sample<Proof>& key) {
		unsigned int hash_value = 0;
		for (unsigned int i = 0; i < key.count; i++)
			hash_value ^= Proof::hash(*key.proofs[i]);
		return hash_value;
	}

	static inline bool is_empty(const theory_sample<Proof>& key) {
		return key.proofs == nullptr;
	}

	static inline void move(const theory_sample<Proof>& src, theory_sample<Proof>& dst) {
		dst.proofs = src.proofs;
		dst.count = src.count;
		dst.log_probability = src.log_probability;
#if !defined(NDEBUG)
		dst.id = src.id;
#endif
	}

	static inline void free(theory_sample<Proof>& sample) { sample.free(); }

private:
	inline void free() {
		for (unsigned int i = 0; i < count; i++) { core::free(*proofs[i]); if (proofs[i]->reference_count == 0) core::free(proofs[i]); }
		core::free(proofs);
	}
};

template<typename Proof>
bool init(theory_sample<Proof>& sample, const hash_set<Proof*>& proofs, double log_probability) {
	sample.log_probability = log_probability;
	sample.proofs = (Proof**) malloc(max((size_t) 1, sizeof(Proof*) * proofs.size));
	if (sample.proofs == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `theory_sample.proofs`.\n");
		return false;
	}

	sample.count = 0;
	for (Proof* proof : proofs) {
		sample.proofs[sample.count] = (Proof*) malloc(sizeof(Proof));
		if (sample.proofs[sample.count] == nullptr
		 || !Proof::clone(*proof, *sample.proofs[sample.count]))
		{
			if (sample.proofs[sample.count] != nullptr) free(sample.proofs[sample.count]);
			for (unsigned int j = 0; j < sample.count; j++)
				free(sample.proofs[j]);
			free(sample.proofs);
			return false;
		}
		sample.count++;
	}

	/* sort the proofs in canonical order */
	if (sample.count > 1)
		sort(sample.proofs, sample.count, pointer_sorter());
	return true;
}

struct null_collector {
	template<typename Proof>
	constexpr inline bool has_prior(const Proof* proof) const {
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept(const hash_set<Proof*>& sample, double proof_prior_diff) const {
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept_with_observation_changes(const hash_set<Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes) const
	{
		return true;
	}
};

template<typename Proof>
inline bool operator == (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr || first.count != second.count)
		return false;
	for (unsigned int i = 0; i < first.count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return false;
	}
	return true;
}

template<typename Proof>
inline bool operator != (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr || first.count != second.count)
		return true;
	for (unsigned int i = 0; i < first.count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return true;
	}
	return false;
}

template<typename Proof>
inline double log_probability(const theory_sample<Proof>& sample) {
	return sample.log_probability;
}

template<typename ProofCalculus, typename Canonicalizer, typename OnProofSampleFunction = no_op>
struct model_evidence_collector
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;

	hash_set<theory_sample<Proof>> samples;
	unsigned int observation_count;
	double current_log_probability;
	Proof* test_proof;
	OnProofSampleFunction on_new_proof_sample;

#if !defined(NDEBUG)
	std::function<double(void)> compute_current_log_probability;
#endif

	template<typename ProofPrior>
	model_evidence_collector(const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior, Proof* test_proof) :
		model_evidence_collector(T, proof_prior, test_proof, *this)
	{ }

	template<typename ProofPrior, typename TheorySampleCollector>
	model_evidence_collector(const theory<ProofCalculus, Canonicalizer>& T,
			ProofPrior& proof_prior, Proof* test_proof, TheorySampleCollector& sample_collector,
			OnProofSampleFunction on_new_proof_sample = no_op()) :
		samples(1024), observation_count(T.observations.size), test_proof(test_proof), on_new_proof_sample(on_new_proof_sample)
	{
#if !defined(NDEBUG)
		if (test_proof->type != ProofType::EXISTENTIAL_INTRODUCTION) {
			fprintf(stderr, "model_evidence_collector ERROR: `test_proof` is not an existential introduction.\n");
			throw std::runtime_error("`test_proof` is not an existential introduction.");
		}
#endif

		/* initialize `current_log_probability` */
		current_log_probability = log_probability(T.observations, proof_prior, sample_collector);
#if !defined(NDEBUG)
		compute_current_log_probability = [&]() { return log_probability(T.observations, proof_prior, sample_collector); };
#endif

		/* add the first sample */
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!init(new_sample, T.observations, current_log_probability))
			throw std::runtime_error("Failed to initialize first theory_sample.");
		on_new_proof_sample(test_proof, current_log_probability);

#if !defined(NDEBUG)
		new_sample.id = samples.size;
#endif
		unsigned int bucket = samples.index_of(new_sample);
		move(new_sample, samples.keys[bucket]);
		samples.size++;
	}

	~model_evidence_collector() { free(); }

	constexpr inline bool has_prior(const Proof* proof) const {
		return (proof != test_proof);
	}

	bool accept(const hash_set<typename ProofCalculus::Proof*>& sample, double proof_prior_diff)
	{
		typedef typename ProofCalculus::Proof Proof;

		current_log_probability += proof_prior_diff;
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!samples.check_size()
		 || !init(new_sample, sample, current_log_probability))
			return false;

		bool contains;
		unsigned int bucket = samples.index_of(new_sample, contains);
		if (contains) {
			/* we've already seen this sample before */
#if !defined(NDEBUG)
			if (fabs(new_sample.log_probability - samples.keys[bucket].log_probability) > 1.0e-9)
				fprintf(stderr, "model_evidence_collector WARNING: Found an"
						" old sample with a different computed log probability."
						" Old log probability: %lf, new log probability: %lf.\n",
						samples.keys[bucket].log_probability, new_sample.log_probability);
#endif
			core::free(new_sample);
			current_log_probability = samples.keys[bucket].log_probability;
			return true;
		}

		/* we've never seen this sample before */
#if !defined(NDEBUG)
		new_sample.id = samples.size;
		double expected_log_probability = compute_current_log_probability();
		if (fabs(expected_log_probability - new_sample.log_probability) > 1.0e-9) {
			fprintf(stderr, "model_evidence_collector WARNING: The computed"
					" log probability of the sample (%lf) differs from the expected log probability (%lf).\n",
					new_sample.log_probability, expected_log_probability);
		}
#endif
		on_new_proof_sample(test_proof, current_log_probability);
		move(new_sample, samples.keys[bucket]);
		samples.size++;
		return true;
	}

	inline bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes)
	{
		if (test_proof != nullptr) {
			for (unsigned int i = 0; i < observation_changes.length; i++) {
				if (observation_changes[i].key == test_proof) {
					test_proof = observation_changes[i].value;
					break;
				}
			}
		}

		return accept(sample, proof_prior_diff);
	}

	inline double total_log_probability() const {
		return logsumexp(samples);
	}

private:
	void free() {
		for (auto& sample : samples)
			core::free(sample);
	}
};

template<typename ProofCalculus, typename Canonicalizer, typename OnProofSampleFunction = no_op>
struct provability_collector
{
	typedef typename ProofCalculus::Proof Proof;

	model_evidence_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction> internal_collector;

	template<typename ProofPrior>
	provability_collector(const theory<ProofCalculus, Canonicalizer>& T,
			ProofPrior& proof_prior, Proof* test_proof,
			OnProofSampleFunction on_new_proof_sample = no_op()) :
		internal_collector(T, proof_prior, test_proof, *this, on_new_proof_sample)
	{ }

	inline bool has_prior(const Proof* proof) const {
		return (proof != internal_collector.test_proof);
	}

	inline bool accept(const hash_set<typename ProofCalculus::Proof*>& sample, double proof_prior_diff)
	{
		return internal_collector.accept(sample, proof_prior_diff);
	}

	bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes)
	{
		return internal_collector.accept_with_observation_changes(sample, proof_prior_diff, observation_changes);
	}

	inline double total_log_probability() const {
		return internal_collector.total_log_probability();
	}
};

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction = no_op>
provability_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction> make_provability_collector(
		const theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofCalculus::Proof* test_proof,
		OnProofSampleFunction on_new_proof_sample = no_op())
{
	return provability_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction>(T, proof_prior, test_proof, on_new_proof_sample);
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, unsigned int num_samples)
{
	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, nullptr);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);
	return collector.total_log_probability();
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability_of_observation(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples)
{
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	Proof* new_proof = T.add_formula(logical_form, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!add(new_proof, proof_prior, proof_axioms)) {
		T.remove_formula(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);

	subtract(collector.test_proof, proof_prior, proof_axioms);
	T.remove_formula(collector.test_proof);
	return collector.total_log_probability();
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability_of_truth(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples)
{
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	Proof* new_proof = T.add_formula(logical_form, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!add(new_proof, proof_prior, proof_axioms)) {
		T.remove_formula(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	provability_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);

	subtract(collector.internal_collector.test_proof, proof_prior, proof_axioms);
	T.remove_formula(collector.internal_collector.test_proof);
	return collector.total_log_probability();
}

template<typename ProofCanonicalizer, typename OnProofSampleFunction>
struct lambda_proof_sample_delegate
{
	OnProofSampleFunction on_new_proof_sample;

	lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) : on_new_proof_sample(on_new_proof_sample) { }

	template<typename Proof>
	inline void operator() (const Proof* test_proof, double log_probability)
	{
		typedef typename Proof::FormulaType Formula;
		typedef typename Formula::Term Term;

		/* get the current value of the term used to introduce the existential quantifier */
		Formula* instantiated_formula = compute_proof_conclusion<built_in_predicates, ProofCanonicalizer>(*test_proof->operands[0]);
		if (test_proof->operands[1]->parameters.length == 0) {
			fprintf(stderr, "lambda_proof_sample_delegate.operator () ERROR: Formula does not depend on the lambda variable.\n");
			core::free(*instantiated_formula); if (instantiated_formula->reference_count == 0) core::free(instantiated_formula);
			return;
		}
		Term* current_term = get_term_at_index(*instantiated_formula, test_proof->operands[1]->parameters[0]);
		on_new_proof_sample(current_term, log_probability);
		free(*instantiated_formula); if (instantiated_formula->reference_count == 0) free(instantiated_formula);
	}
};

template<typename ProofCanonicalizer, typename OnProofSampleFunction>
lambda_proof_sample_delegate<ProofCanonicalizer, OnProofSampleFunction> make_lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) {
	return lambda_proof_sample_delegate<ProofCanonicalizer, OnProofSampleFunction>(on_new_proof_sample);
}

/* TODO: for debugging; delete this */
#include <core/random.h>

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction>
bool log_joint_probability_of_lambda(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		OnProofSampleFunction on_new_proof_sample)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

#if !defined(NDEBUG)
	typedef typename Formula::Type FormulaType;
	if (logical_form->type != FormulaType::LAMBDA)
		fprintf(stderr, "log_joint_probability_of_lambda WARNING: `logical_form` is not a lambda expression.\n");
#endif

	Formula* existential = Formula::new_exists(logical_form->quantifier.variable, logical_form->quantifier.operand);
	if (existential == nullptr)
		return false;
	existential->quantifier.operand->reference_count++;

	unsigned int new_constant;
extern const string_map_scribe* debug_terminal_printer;
T.print_axioms(stderr, *debug_terminal_printer);
	Proof* new_proof = T.add_formula(existential, new_constant);
	free(*existential); if (existential->reference_count == 0) free(existential);
	if (new_proof == nullptr) {
		return false;
	} else if (!proof_axioms.template add<false>(new_proof, proof_prior)) {
		T.remove_formula(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	auto new_proof_sample_delegate = make_lambda_proof_sample_delegate<typename ProofCalculus::ProofCanonicalizer>(on_new_proof_sample);
	auto collector = make_provability_collector(T, proof_prior, new_proof, new_proof_sample_delegate);
	for (unsigned int t = 0; t < num_samples; t++)
{
fprintf(stderr, "DEBUG: t = %u\n", t);
proof_axioms.check_proof_axioms(T);
proof_axioms.check_universal_eliminations(T, collector);
T.check_concept_axioms();
T.check_disjunction_introductions();
T.sets.are_elements_unique();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
T.sets.are_set_sizes_valid();
T.sets.check_set_ids();
if (!T.observations.contains(collector.internal_collector.test_proof))
	fprintf(stderr, "log_joint_probability_of_lambda WARNING: `provability_collector.internal_collector.test_proof` is not an observation in the theory.\n");
/*if (t == 34)
fprintf(stderr, "DEBUG: BREAKPOINT\n");*/
T.print_axioms(stderr, *debug_terminal_printer);
		do_mh_step(T, proof_prior, proof_axioms, collector);
		if (t % 10 == 0)
		{
			/* make sure we are sampling the query logical form with sufficient frequency */
			unsigned int index;
			for (index = 0; index < T.existential_intro_nodes.length; index++)
				if (T.existential_intro_nodes[index].value == collector.internal_collector.test_proof) break;

		}

if (t % 1000 == 0)
	fprintf(stdout, "(seed = %u)\n", get_seed());
}

	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, proof_prior);
	T.remove_formula(collector.internal_collector.test_proof);
	return true;
}

#endif /* THEORY_H_ */
