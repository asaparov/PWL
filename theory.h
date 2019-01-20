#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <math/multiset.h>

#include "array_view.h"

using namespace core;


enum built_in_predicates : unsigned int {
	PREDICATE_UNKNOWN = 1,
	PREDICATE_COUNT
};

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("unknown", PREDICATE_UNKNOWN);
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

template<bool Unique, typename T>
inline void add_sorted(array<T>& list, const T& element) {
	unsigned int index = linear_search(list.data, element, 0, list.length);
	if (Unique && index < list.length && list[index] == element) return;
	shift_right(list.data, list.length, index);
	list[index] = element;
	list.length++;
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

	array<Proof*> universal_quantifications;

	hash_set<Proof*> observations;
	hash_multiset<Formula*, false> proof_axioms;

	Canonicalizer canonicalizer;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), types(64), relations(64),
			ground_concepts(64), ground_axiom_count(0), universal_quantifications(32),
			observations(64), proof_axioms(64) { }

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
		} for (Proof* proof : universal_quantifications) {
			core::free(*proof);
			if (proof->reference_count == 0)
				core::free(proof);
		} for (auto entry : ground_concepts) {
			core::free(entry.value);
		}
	}

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (Proof* axiom : universal_quantifications) {
			if (!print(*axiom->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : ground_concepts) {
			if (!entry.value.print_axioms(out, std::forward<Printer>(printer)...)) return false;
		}
		return true;
	}

	inline unsigned int axiom_count() const {
		return universal_quantifications.length + ground_axiom_count;
	}

	bool add_formula(Formula* formula, unsigned int& new_constant)
	{
		Formula* canonicalized = canonicalize(*formula, canonicalizer);
		if (canonicalized == NULL) return false;

		Proof* new_proof = make_proof(canonicalized, new_constant);
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put(PREDICATE_UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof(*new_proof, expected_conclusion, canonicalizer))
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

	inline bool add_universal_quantification(Proof* axiom) {
		if (!universal_quantifications.add(axiom))
			return false;
		axiom->reference_count++;
		return true;
	}

	inline void remove_universal_quantification(unsigned int axiom_index) {
		Proof* axiom = universal_quantifications[axiom_index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		universal_quantifications.remove(axiom_index);
	}

	template<bool Negated>
	inline bool add_unary_atom(unsigned int predicate, unsigned int arg, Proof* axiom)
	{
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
	inline void remove_unary_atom(unsigned int predicate, unsigned int arg)
	{
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
	}

	template<bool Negated>
	inline void remove_binary_atom(relation rel)
	{
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
	}

private:
	Proof* make_proof(Formula* canonicalized, unsigned int& new_constant)
	{
		unsigned int predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<false>(canonicalized, predicate, arg1, arg2, new_constant);
		} else if (canonicalized->type == FormulaType::NOT) {
			if (is_atomic(*canonicalized->unary.operand, predicate, arg1, arg2)) {
				return make_atom_proof<true>(canonicalized, predicate, arg1, arg2, new_constant);
			}
		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == FormulaType::IF_THEN) {
				const Formula* left = canonicalized->quantifier.operand->binary.left;
				Term const* arg1; Term const* arg2;
				bool atomic = is_atomic(*left, predicate, arg1, arg2);
				if (atomic && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 ==  NULL)
				{
					if (predicate == PREDICATE_UNKNOWN) {
						/* this is a definition of a type */
						Formula* right = canonicalized->quantifier.operand->binary.right;

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
						Formula* definition = Formula::new_for_all(variable, Formula::new_if_then(
								Formula::new_atom(new_constant, Formula::new_variable(variable)), right));
						if (definition == NULL) return NULL;
						right->reference_count++;

						Proof* new_axiom = ProofCalculus::new_axiom(definition);
						free(*definition); if (definition->reference_count == 0) free(definition);
						if (new_axiom == NULL) return NULL;
						new_axiom->reference_count++;
						if (!add_universal_quantification(new_axiom)) {
							free(*new_axiom); free(new_axiom);
							return NULL;
						}
						return new_axiom;

						/*Formula* definition = Formula::new_for_all(variable, Formula::new_iff(
								Formula::new_atom(new_constant, Formula::new_variable(variable)), right));
						if (definition == NULL) return NULL;
						right->reference_count++;

						Proof* proof = ProofCalculus::new_universal_intro(
							ProofCalculus::new_implication_intro(
								ProofCalculus::new_biconditional_elim_left(
									ProofCalculus::new_universal_elim(ProofCalculus::new_axiom(definition), Formula::new_parameter(1)),
									ProofCalculus::new_axiom(left)),
								ProofCalculus::new_axiom(left)), 1);
						if (proof == NULL) {
							free(*definition); free(definition);
						}
						return proof;*/
					} else {
						/* this is a formula of form `![x]:(t(x) => f(x))` */
						/* TODO: check that this is not implied by an existing univerally-quantified axiom */
						Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
						if (new_axiom == NULL) return NULL;
						new_axiom->reference_count++;
						if (!add_universal_quantification(new_axiom)) {
							free(*new_axiom); free(new_axiom);
							return NULL;
						}
						return new_axiom;
					}
				} else if (is_atomic(*left)) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				} else if (left->type == FormulaType::AND) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}
			}
		} else if (canonicalized->type == FormulaType::AND) {
			return ProofCalculus::new_conjunction_intro(
					make_proof(canonicalized->binary.left, new_constant),
					make_proof(canonicalized->binary.right, new_constant));
		} else if (canonicalized->type == FormulaType::EXISTS) {
			/* TODO: implement this */
			fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
			return NULL;
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

	template<bool Negated>
	Proof* make_atom_proof(Formula* canonicalized, unsigned int predicate,
			Term* arg1, Term* arg2, unsigned int& new_constant)
	{
		bool contains; unsigned int bucket;
		if (arg1->type == TermType::CONSTANT && arg2 == NULL)
		{
			/* this is a unary formula */
			if (predicate != PREDICATE_UNKNOWN) {
				if (arg1->constant == PREDICATE_UNKNOWN) {
					/* this is a definition of an object */
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

					/* check that this is not implied by an existing univerally-quantified axiom */
					for (unsigned int k = 0; k < universal_quantifications.length; k++) {
						Proof* axiom = universal_quantifications[k];
						Formula* consequent = axiom->formula->quantifier.operand->binary.right;
						unsigned int i = 0;
						if (consequent->type == FormulaType::AND) {
							for (; i < consequent->array.length; i++) {
								Term const* dst = NULL;
								Term* variable = Formula::new_variable(axiom->formula->quantifier.variable);
								if (unify(*consequent->array.operands[i], *canonicalized, variable, dst)) {
									free(*variable); if (variable->reference_count == 0) free(variable);
									break;
								}
								free(*variable); if (variable->reference_count == 0) free(variable);
							}
							if (i == consequent->array.length) continue;
						} else if (*consequent != *canonicalized) {
							i = consequent->array.length;
							continue;
						}

						/* check if the antecedent is satisfied by `c` */
						Proof* proof;
						if (!make_universal_elim_proof(consequent, ground_concepts.get(arg1->constant), proof))
							return NULL;
						if (proof != NULL) {
							proof->reference_count++;
							Proof* new_proof;
							if (i == consequent->array.length) {
								new_proof = ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, arg1), proof);
							} else {
								new_proof = ProofCalculus::new_conjunction_elim(
										ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, arg1), proof),
										make_array_view(&i, 1));
							}
							free(*proof); if (proof->reference_count == 0) free(proof);
							if (new_proof == NULL) return NULL;
							return new_proof;
						}
					}

					/* no existing universally-quantified axiom implies this observation */
					Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated>(predicate, arg1->constant, new_axiom)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
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
			if (predicate != PREDICATE_UNKNOWN)
			{
				if (arg1->constant == PREDICATE_UNKNOWN
				 || arg2->constant == PREDICATE_UNKNOWN) {
					fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
					return NULL;
				} else {
					/* this is a formula of form `r(c_1,c_2)` */

					/* check that this is not implied by an existing univerally-quantified axiom */
					for (unsigned int k = 0; k < universal_quantifications.length; k++) {
						Proof* axiom = universal_quantifications[k];
						Formula* consequent = axiom->formula->quantifier.operand->binary.right;
						Term const* unifying_term = NULL;
						unsigned int i = 0;
						if (consequent->type == FormulaType::AND) {
							for (; i < consequent->array.length; i++) {
								Term* variable = Formula::new_variable(axiom->formula->quantifier.variable);
								if (unify(*consequent->array.operands[i], *canonicalized, variable, unifying_term)) {
									free(*variable); if (variable->reference_count == 0) free(variable);
									break;
								}
								free(*variable); if (variable->reference_count == 0) free(variable);
							}
							if (i == consequent->array.length) continue;
						} else if (*consequent != *canonicalized) {
							i = consequent->array.length;
							continue;
						}

						/* check if the antecedent is satisfied by `c_1` or `c_2` (whichever unifies with the consequent) */
						Proof* proof;
						if (!make_universal_elim_proof(consequent, ground_concepts.get(unifying_term->constant), proof))
							return NULL;
						if (proof != NULL) {
							proof->reference_count++;
							Proof* new_proof;
							if (i == consequent->array.length) {
								new_proof = ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, arg1), proof);
							} else {
								new_proof = ProofCalculus::new_conjunction_elim(
										ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, arg1), proof),
										make_array_view(&i, 1));
							}
							free(*proof); if (proof->reference_count == 0) free(proof);
							if (new_proof == NULL) return NULL;
							return new_proof;
						}
					}

					/* no existing universally-quantified axiom implies this observation */
					Proof* new_axiom = ProofCalculus::new_axiom(canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated>({predicate, arg1->constant, arg2->constant}, new_axiom)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
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
			return predicate != PREDICATE_UNKNOWN
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
};

#endif /* THEORY_H_ */
