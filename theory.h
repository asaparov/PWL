#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <math/multiset.h>

#include "array_view.h"
#include "set_reasoning.h"

using namespace core;


enum class built_in_predicates : unsigned int {
	UNKNOWN = 1,
	SET = 2,
	ALL = 3,
	SIZE = 4,
	COUNT
};

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("unknown",	(unsigned int) built_in_predicates::UNKNOWN)
		&& names.put("set",		(unsigned int) built_in_predicates::SET)
		&& names.put("all",		(unsigned int) built_in_predicates::ALL)
		&& names.put("size",	(unsigned int) built_in_predicates::SIZE);
}

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

	hash_map<unsigned int, concept<ProofCalculus>> ground_concepts;
	unsigned int ground_axiom_count;

	hash_set<Proof*> observations;
	hash_multiset<Formula*, false> proof_axioms;
	set_reasoning<built_in_predicates, Formula, ProofCalculus> sets;

	Proof* empty_set_axiom;

	Canonicalizer canonicalizer;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), types(64), relations(64),
			ground_concepts(64), ground_axiom_count(0), observations(64), proof_axioms(64)
	{
		Formula* empty_set_formula = Formula::new_for_all(1, Formula::new_equals(
			Formula::new_equals(Formula::new_atom((unsigned int) built_in_predicates::SIZE, Formula::new_lambda(2, Formula::new_variable(1))), Formula::new_int(0)),
			Formula::new_not(Formula::new_exists(2, Formula::new_variable(1)))
		));
		if (empty_set_formula == NULL) exit(EXIT_FAILURE);
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
		} for (auto entry : ground_concepts) {
			core::free(entry.value);
		}
		core::free(*empty_set_axiom);
		if (empty_set_axiom->reference_count == 0)
			core::free(empty_set_axiom);
	}

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (auto entry : ground_concepts) {
			if (!entry.value.print_axioms(out, std::forward<Printer>(printer)...)) return false;
		}
		return sets.print_axioms(out, std::forward<Printer>(printer)...);
	}

	bool add_formula(Formula* formula, unsigned int& new_constant)
	{
		new_constant = 0;
		Formula* canonicalized = canonicalize(*formula, canonicalizer);
		if (canonicalized == NULL) return false;

		Proof* new_proof = make_proof<true>(canonicalized, new_constant);
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put((unsigned int) built_in_predicates::UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof<built_in_predicates>(*new_proof, expected_conclusion, canonicalizer))
fprintf(stderr, "add_formula WARNING: `check_proof` failed.\n");
free(*expected_conclusion); if (expected_conclusion->reference_count == 0) free(expected_conclusion);
		core::free(*canonicalized);
		if (canonicalized->reference_count == 0)
			core::free(canonicalized);
		if (new_proof == NULL) {
			return false;
		} else if (!observations.add(new_proof)) {
			core::free(*new_proof);
			if (new_proof->reference_count == 0)
				core::free(new_proof);
			return false;
		}

		/* add the axioms in the new proof to `proof_axioms` */
		if (!get_axioms(new_proof, proof_axioms)) {
			observations.remove(new_proof);
			core::free(*new_proof);
			if (new_proof->reference_count == 0)
				core::free(new_proof);
			return false;
		}
		return true;
	}

	template<bool Negated>
	inline bool add_unary_atom(unsigned int predicate, unsigned int arg, Proof* axiom)
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
		if (!move_element_to_subset(arg, lifted_literal)) {
			fprintf(stderr, "theory.add_unary_atom ERROR: `move_element_to_subset` failed.\n");
			free(*lifted_literal); free(lifted_literal);
			return false;
		}
		free(*lifted_literal); free(lifted_literal);

#if !defined(NDEBUG)
		if (!ground_concepts.table.contains(arg))
			fprintf(stderr, "theory.add_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", arg);
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
		array_map<unsigned int, Proof*>& ground_types = (Negated ? ground_concepts.get(arg).negated_types : ground_concepts.get(arg).types);
		if (!instances.ensure_capacity(instances.length + 1)
		 || !ground_types.ensure_capacity(ground_types.size + 1)) return false;

		add_sorted<false>(instances, arg);
		ground_types.keys[ground_types.size] = predicate;
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		ground_axiom_count++;
		return true;
	}

	template<bool Negated>
	inline bool add_binary_atom(relation rel, Proof* axiom)
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
			if (!move_element_to_subset(rel.arg1, lifted_literal)) {
				fprintf(stderr, "theory.add_binary_atom ERROR: `move_element_to_subset` failed.\n");
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
			if (!move_element_to_subset(rel.arg1, lifted_literal)) {
				fprintf(stderr, "theory.add_binary_atom ERROR: `move_element_to_subset` failed.\n");
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
			if (!move_element_to_subset(rel.arg2, lifted_literal)) {
				fprintf(stderr, "theory.add_binary_atom ERROR: `move_element_to_subset` failed.\n");
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			free(*lifted_literal); free(lifted_literal);
		}

#if !defined(NDEBUG)
		if (!ground_concepts.table.contains(rel.arg1))
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (!ground_concepts.table.contains(rel.arg2))
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
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts.get(rel.arg1).negated_relations : ground_concepts.get(rel.arg1).relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts.get(rel.arg2).negated_relations : ground_concepts.get(rel.arg2).relations);
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

	template<bool Negated>
	inline bool remove_unary_atom(unsigned int predicate, unsigned int arg)
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
		if (!move_element_to_superset(arg, lifted_literal)) {
			fprintf(stderr, "theory.remove_unary_atom ERROR: `move_element_to_superset` failed.\n");
			if (lifted_literal != NULL) { free(*lifted_literal); free(lifted_literal); }
			return false;
		}
		free(*lifted_literal); free(lifted_literal);

#if !defined(NDEBUG)
		if (!types.table.contains(predicate))
			fprintf(stderr, "theory.remove_unary_atom WARNING: `types` does not contain the key %u.\n", predicate);
		if (!ground_concepts.table.contains(arg))
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", arg);
#endif

		array<unsigned int>& instances = (Negated ? types.get(predicate).value : types.get(predicate).key);
		array_map<unsigned int, Proof*>& ground_types = (Negated ? ground_concepts.get(arg).negated_types : ground_concepts.get(arg).types);

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

	template<bool Negated>
	inline bool remove_binary_atom(relation rel)
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
			if (!move_element_to_superset(rel.arg1, lifted_literal)) {
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
			if (!move_element_to_superset(rel.arg1, lifted_literal)) {
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
			if (!move_element_to_superset(rel.arg2, lifted_literal)) {
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
		if (!ground_concepts.table.contains(rel.arg1))
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (!ground_concepts.table.contains(rel.arg2))
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg2);
#endif

		pair<array<unsigned int>, array<unsigned int>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0});

		array<unsigned int>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<unsigned int>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<unsigned int>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts.get(rel.arg1).negated_relations : ground_concepts.get(rel.arg1).relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts.get(rel.arg2).negated_relations : ground_concepts.get(rel.arg2).relations);

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

private:
	template<bool DefinitionsAllowed>
	Proof* make_proof(Formula* canonicalized, unsigned int& new_constant)
	{
		unsigned int predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<DefinitionsAllowed, false>(canonicalized, predicate, arg1, arg2, new_constant);
		} else if (canonicalized->type == FormulaType::NOT) {
			if (is_atomic(*canonicalized->unary.operand, predicate, arg1, arg2)) {
				return make_atom_proof<DefinitionsAllowed, true>(canonicalized, predicate, arg1, arg2, new_constant);
			} else if (canonicalized->type == FormulaType::EXISTS) {
				/* get empty set size axiom, forcing inconsistencies to be resolved */
				Formula* set_formula = canonicalized->quantifier.operand;
				Proof* set_size_axiom = sets.template get_size_axiom<true>(set_formula, 0);
				if (set_size_axiom == NULL) return NULL;

				/* TODO: what are the semantics of the third argument to new_equality_elim? */
				unsigned int zero = 0;
				Proof* proof = ProofCalculus::new_equality_elim(set_size_axiom,
						ProofCalculus::new_universal_elim(empty_set_axiom, set_formula), make_array_view(&zero, 1));
				if (proof == NULL) return NULL;
				proof->reference_count++;
				return proof;
			}
		} else if (canonicalized->type == FormulaType::FOR_ALL) {
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

					if (!ground_concepts.check_size()) return NULL;
					new_constant = new_constant_offset + ground_concepts.table.size;

					bool contains; unsigned int bucket;
					concept<ProofCalculus>& c = ground_concepts.get(new_constant, contains, bucket);
					if (!contains) {
						if (!init(c)) return NULL;
						ground_concepts.table.keys[bucket] = new_constant;
						ground_concepts.table.size++;
					}

					/* TODO: definitions of new concepts should be biconditionals */
					Formula* new_left = Formula::new_atom(new_constant, Formula::new_variable(variable));
					if (new_left == NULL) return NULL;

					Proof* new_axiom = sets.template get_subset_axiom<true>(new_left, right);
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
					Proof* new_axiom = sets.template get_subset_axiom<true>(left, right);
					new_axiom->reference_count++;
					return new_axiom;
				}
			}
		} else if (canonicalized->type == FormulaType::AND) {
			Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
			if (operands == NULL) {
				fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
				return NULL;
			}
			for (unsigned int i = 0; i < canonicalized->array.length; i++) {
				if (new_constant == 0)
					operands[i] = make_proof<true>(canonicalized->array.operands[i], new_constant);
				else operands[i] = make_proof<false>(canonicalized->array.operands[i], new_constant);

				if (operands[i] == NULL) {
					for (unsigned int j = 0; j < i; j++) {
						free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
					}
				}
			}
			Proof* conjunction = ProofCalculus::new_conjunction_intro(make_array_view(operands, canonicalized->array.length));
			for (unsigned int j = 0; j < canonicalized->array.length; j++) {
				free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
			}
			return conjunction;
		} else if (canonicalized->type == FormulaType::EXISTS) {
			/* TODO: implement this */
			fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
			return NULL;
		} else if (canonicalized->type == FormulaType::EQUALS) {
			Formula* left = canonicalized->binary.left;
			Formula* right = canonicalized->binary.right;
			Term* arg1; Term* arg2;
			bool atomic = is_atomic(*left, predicate, arg1, arg2);
			if (atomic && predicate == (unsigned int) built_in_predicates::SIZE
			 && arg1->type == TermType::LAMBDA && arg2 == NULL
			 && right->type == TermType::INTEGER)
			{
				Proof* set_size_axiom = sets.template get_size_axiom<true>(arg1->quantifier.operand, 0);
				set_size_axiom->reference_count++;
				return set_size_axiom;
			}
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

	template<bool DefinitionsAllowed, bool Negated>
	Proof* make_atom_proof(Formula* canonicalized, unsigned int predicate,
			Term* arg1, Term* arg2, unsigned int& new_constant)
	{
		bool contains; unsigned int bucket;
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
					if (!ground_concepts.check_size()) return NULL;
					new_constant = new_constant_offset + ground_concepts.table.size;
					arg1->constant = new_constant;

					Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;

					concept<ProofCalculus>& c = ground_concepts.get(new_constant, contains, bucket);
					if (!contains) {
						if (!init(c)) {
							free(*new_axiom); free(new_axiom);
							return NULL;
						}
						ground_concepts.table.keys[bucket] = new_constant;
						ground_concepts.table.size++;
					} if (!add_unary_atom<Negated>(predicate, new_constant, new_axiom)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
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

					/* there is no extensional edge that proves this atomic formula, so create a new axiom */
					Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated>(predicate, arg1->constant, new_axiom)) {
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

					/* there is no extensional edge that proves this atomic formula, so create a new axiom */
					Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated>({predicate, arg1->constant, arg2->constant}, new_axiom)) {
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

	bool move_element_to_subset(unsigned int element, Formula* lifted_literal)
	{
		Formula* new_set_formula;
		Formula* old_set_formula = make_lifted_conjunction(element, *this);
		if (old_set_formula == NULL) {
			return false;
		} else if (old_set_formula->array.length == 0) {
			/* there are no other atomic formulas for this concept */
			free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			old_set_formula = NULL;
			new_set_formula = lifted_literal;
			new_set_formula->reference_count++;
		} else {
			if (old_set_formula->array.length == 1) {
				Formula* singleton = old_set_formula->array.operands[0];
				singleton->reference_count++;
				free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
				old_set_formula = singleton;
			}
			new_set_formula = Formula::new_and(old_set_formula, lifted_literal);
			old_set_formula->reference_count++;
			lifted_literal->reference_count++;
		}
		if (new_set_formula == NULL) {
			if (old_set_formula != NULL) {
				free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			}
			return false;
		}
		Formula* canonicalized = canonicalize(*new_set_formula, canonicalizer);
		free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		if (canonicalized == NULL) {
			if (old_set_formula != NULL) {
				free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			}
			return false;
		}
		if (!sets.template move_element_to_set<true>(element, old_set_formula, canonicalized)) {
			free(*canonicalized); free(canonicalized);
			if (old_set_formula != NULL) {
				free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			}
			return false;
		}

		/* make sure we can't currently prove the negation of this new axiom
		   via set containment (can we prove that the instance does not belong to any ancestor?) */
		unsigned int new_set_id = sets.set_ids.get(*canonicalized);
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == NULL) continue;
			Formula* set_formula = sets.sets[i].set_formula();
			Formula* negated_set_formula = Formula::new_not(set_formula);
			if (negated_set_formula == NULL) {
				if (old_set_formula == NULL) sets.remove_element_from_set(element, canonicalized);
				else sets.move_element_to_superset(element, canonicalized, old_set_formula);
				if (old_set_formula != NULL) { free(*old_set_formula); free(old_set_formula); }
				free(*canonicalized); free(canonicalized);
				return false;
			}
			set_formula->reference_count++;
			bool contradiction = is_subset(lifted_literal, negated_set_formula) && sets.sets[i].descendants.contains(new_set_id);
			free(*negated_set_formula); free(negated_set_formula);
			if (contradiction) {
				if (old_set_formula == NULL) sets.remove_element_from_set(element, canonicalized);
				else sets.move_element_to_superset(element, canonicalized, old_set_formula);
				if (old_set_formula != NULL) { free(*old_set_formula); free(old_set_formula); }
				free(*canonicalized); free(canonicalized);
				return false;
			}
		}

		if (old_set_formula != NULL) { free(*old_set_formula); free(old_set_formula); }
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		return true;
	}

	bool move_element_to_superset(unsigned int element, Formula* lifted_literal)
	{
		Formula* old_set_formula = make_lifted_conjunction(element, *this);
		if (old_set_formula == NULL) {
			return false;
		} else if (old_set_formula->array.length == 0) {
			/* there are no other atomic formulas for this concept */
			free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			return false;
		}

		Formula* canonicalized_old_set_formula = canonicalize(*old_set_formula, canonicalizer);
		if (canonicalized_old_set_formula == NULL) {
			free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);
			return false;
		}
		free(*old_set_formula); if (old_set_formula->reference_count == 0) free(old_set_formula);

		Formula* canonicalized = canonicalize(*lifted_literal, canonicalizer);
		array<Formula*> difference(8);
		if (canonicalized == NULL) {
			free(*canonicalized_old_set_formula);
			if (canonicalized_old_set_formula->reference_count == 0)
				free(canonicalized_old_set_formula);
			return false;
		} else if (canonicalized->type == FormulaType::AND) {
			unsigned int i = 0, j = 0;
			while (i < canonicalized_old_set_formula->array.length && j < canonicalized->array.length)
			{
				if (*canonicalized_old_set_formula->array.operands[i] == *canonicalized->array.operands[j]) {
					i++; j++;
				} else if (*canonicalized_old_set_formula->array.operands[i] < *canonicalized->array.operands[j]) {
					if (!difference.add(canonicalized_old_set_formula->array.operands[i])) {
						free(*canonicalized_old_set_formula);
						if (canonicalized_old_set_formula->reference_count == 0)
							free(canonicalized_old_set_formula);
						free(*canonicalized); free(canonicalized); return false;
					}
					i++;
				} else {
					j++;
				}
			}

			while (i < canonicalized_old_set_formula->array.length) {
				if (!difference.add(canonicalized_old_set_formula->array.operands[i])) {
					free(*canonicalized_old_set_formula);
					if (canonicalized_old_set_formula->reference_count == 0)
						free(canonicalized_old_set_formula);
					free(*canonicalized); free(canonicalized); return false;
				}
				i++;
			}
		} else {
			for (unsigned int i = 0; i < canonicalized_old_set_formula->array.length; i++) {
				if (*canonicalized_old_set_formula->array.operands[i] != *canonicalized
				 && !difference.add(canonicalized_old_set_formula->array.operands[i]))
				{
					free(*canonicalized_old_set_formula);
					if (canonicalized_old_set_formula->reference_count == 0)
						free(canonicalized_old_set_formula);
					free(*canonicalized); free(canonicalized); return false;
				}
			}
		}
		free(*canonicalized); free(canonicalized);

		Formula* new_set_formula;
		if (difference.length == 0)
			new_set_formula = NULL;
		else if (difference.length == 1)
			new_set_formula = difference[0];
		else new_set_formula = Formula::new_and(difference);
		if (difference.length != 0 && new_set_formula == NULL) {
			free(*canonicalized_old_set_formula);
			if (canonicalized_old_set_formula->reference_count == 0)
				free(canonicalized_old_set_formula);
			return false;
		}
		for (Formula* conjunct : difference)
			conjunct->reference_count++;

		bool success;
		if (new_set_formula == NULL) {
			sets.remove_element_from_set(element, canonicalized_old_set_formula);
			success = true;
		} else {
			success = sets.move_element_to_superset(element, canonicalized_old_set_formula, new_set_formula);
			free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		}
		free(*canonicalized_old_set_formula);
		if (canonicalized_old_set_formula->reference_count == 0)
			free(canonicalized_old_set_formula);
		return success;
	}

	inline bool add_new_element_to_set(unsigned int element, Formula* lifted_literal) {
		return sets.move_element_to_set<true>(element, NULL, lifted_literal);
	}
};

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
