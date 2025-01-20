#ifndef NATURAL_DEDUCTION_H_
#define NATURAL_DEDUCTION_H_

#include <core/array.h>
#include <math/log.h>
#include <math/multiset.h>

#include "array_view.h"

using namespace core;

typedef uint_fast16_t nd_step_type_specifier;
enum class nd_step_type : nd_step_type_specifier
{
	AXIOM = 0,
	PARAMETER,
	ARRAY_PARAMETER,
	TERM_PARAMETER,
	FORMULA_PARAMETER,

	BETA_EQUIVALENCE,

	CONJUNCTION_INTRODUCTION,
	CONJUNCTION_ELIMINATION, /* for the canonicalizing proof calculus */
	CONJUNCTION_ELIMINATION_LEFT,
	CONJUNCTION_ELIMINATION_RIGHT,
	DISJUNCTION_INTRODUCTION, /* for the canonicalizing proof calculus */
	DISJUNCTION_INTRODUCTION_LEFT,
	DISJUNCTION_INTRODUCTION_RIGHT,
	DISJUNCTION_ELIMINATION,
	IMPLICATION_INTRODUCTION,
	IMPLICATION_ELIMINATION,
	BICONDITIONAL_INTRODUCTION,
	BICONDITIONAL_ELIMINATION_LEFT,
	BICONDITIONAL_ELIMINATION_RIGHT,
	EQUALITY_ELIMINATION,
	PROOF_BY_CONTRADICTION,
	NEGATION_ELIMINATION,
	FALSITY_ELIMINATION,
	COMPARISON_INTRODUCTION,
	INEQUALITY_INTRODUCTION,

	UNIVERSAL_INTRODUCTION,
	UNIVERSAL_ELIMINATION,
	EXISTENTIAL_INTRODUCTION,
	EXISTENTIAL_ELIMINATION,

	COUNT
};

template<bool Intuitionistic, typename Stream>
inline bool print(nd_step_type type, Stream& out) {
	switch (type) {
	case nd_step_type::AXIOM: return print("Ax", out);
	case nd_step_type::PARAMETER: return print("PARAMETER", out);
	case nd_step_type::ARRAY_PARAMETER: return print("ARRAY_PARAMETER", out);
	case nd_step_type::TERM_PARAMETER: return print("TERM_PARAMETER", out);
	case nd_step_type::FORMULA_PARAMETER: return print("FORMULA_PARAMETER", out);
	case nd_step_type::BETA_EQUIVALENCE: return print("β", out);
	case nd_step_type::CONJUNCTION_INTRODUCTION: return print("∧I", out);
	case nd_step_type::CONJUNCTION_ELIMINATION: return print("∧E", out);
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT: return print("∧Eᴸ", out);
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT: return print("∧Eᴿ", out);
	case nd_step_type::DISJUNCTION_INTRODUCTION: return print("∨I", out);
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT: return print("∨Iᴸ", out);
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT: return print("∨Iᴿ", out);
	case nd_step_type::DISJUNCTION_ELIMINATION: return print("∨E", out);
	case nd_step_type::IMPLICATION_INTRODUCTION: return print("→I", out);
	case nd_step_type::IMPLICATION_ELIMINATION: return print("→E", out);
	case nd_step_type::BICONDITIONAL_INTRODUCTION: return print("↔I", out);
	case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT: return print("↔Eᴸ", out);
	case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT: return print("↔Eᴿ", out);
	case nd_step_type::EQUALITY_ELIMINATION: return print("=E", out);
	case nd_step_type::PROOF_BY_CONTRADICTION: return print(Intuitionistic ? "¬I" : "⊥ᶜ", out);
	case nd_step_type::NEGATION_ELIMINATION: return print("¬E", out);
	case nd_step_type::FALSITY_ELIMINATION: return print("⊥E", out);
	case nd_step_type::COMPARISON_INTRODUCTION: return print("≥I", out);
	case nd_step_type::INEQUALITY_INTRODUCTION: return print("≠I", out);
	case nd_step_type::UNIVERSAL_INTRODUCTION: return print("∀I", out);
	case nd_step_type::UNIVERSAL_ELIMINATION: return print("∀E", out);
	case nd_step_type::EXISTENTIAL_INTRODUCTION: return print("∃I", out);
	case nd_step_type::EXISTENTIAL_ELIMINATION: return print("∃E", out);
	case nd_step_type::COUNT: break;
	}
	fprintf(stderr, "print ERROR: Unrecognized `nd_step_type`.\n");
	return false;
}

constexpr static unsigned int ND_OPERAND_COUNT = 3;
static double LOG_ND_RULE_COUNT = log((double) nd_step_type::COUNT);

template<typename Formula>
struct nd_step
{
	typedef Formula FormulaType;
	typedef typename Formula::Term Term;

	nd_step_type type;
	unsigned int reference_count;
	union {
		Term* term;
		Formula* formula;
		unsigned int parameter;
		array<unsigned int> parameters;

		nd_step<Formula>* operands[ND_OPERAND_COUNT];
		array<nd_step<Formula>*> operand_array;
	};

	/* list of proof steps that use the subproof rooted at this rule instance */
	array<nd_step<Formula>*> children;

	static inline void free(nd_step<Formula>& step) {
		step.reference_count--;
		if (step.reference_count == 0) {
			step.free();
			core::free(step.children);
		}
	}

	inline bool is_parameter() const {
		switch (type) {
		case nd_step_type::PARAMETER:
		case nd_step_type::TERM_PARAMETER:
		case nd_step_type::ARRAY_PARAMETER:
		case nd_step_type::FORMULA_PARAMETER:
			return true;
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::BETA_EQUIVALENCE:
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
			return false;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.is_parameter ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	inline void get_subproofs(nd_step<Formula>* const*& subproofs, unsigned int& length) {
		switch (type) {
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::PARAMETER:
		case nd_step_type::TERM_PARAMETER:
		case nd_step_type::ARRAY_PARAMETER:
		case nd_step_type::FORMULA_PARAMETER:
			subproofs = NULL; length = 0;
			return;
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			subproofs = operand_array.data;
			length = operand_array.length;
			return;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			subproofs = operands;
			length = ND_OPERAND_COUNT;
			return;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.get_subproofs ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	inline void get_subproofs(const nd_step<Formula>* const*& subproofs, unsigned int& length) const {
		switch (type) {
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::PARAMETER:
		case nd_step_type::TERM_PARAMETER:
		case nd_step_type::ARRAY_PARAMETER:
		case nd_step_type::FORMULA_PARAMETER:
			subproofs = NULL; length = 0;
			return;
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			subproofs = operand_array.data;
			length = operand_array.length;
			return;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			subproofs = operands;
			length = ND_OPERAND_COUNT;
			return;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.get_subproofs ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	template<typename ProofStep>
	inline void remove_child(const ProofStep& child) {
		unsigned int index = children.index_of(child);
#if !defined(NDEBUG)
		if (index == children.length) {
			fprintf(stderr, "nd_step.remove_child WARNING: Index out of bounds.\n");
			return;
		}
#endif
		children.remove(index);
	}

	inline void free() {
		switch (type) {
		case nd_step_type::PARAMETER:
			return;
		case nd_step_type::ARRAY_PARAMETER:
			core::free(parameters); return;
		case nd_step_type::TERM_PARAMETER:
			core::free(*term);
			if (term->reference_count == 0)
				core::free(term);
			return;
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::FORMULA_PARAMETER:
			core::free(*formula);
#if !defined(NDEBUG)
			if (children.length > 0 && formula->reference_count == 0)
				fprintf(stderr, "nd_step.free WARNING: Detected double free.\n");
#endif
			if (formula->reference_count == 0)
				core::free(formula);
			return;
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			for (unsigned int i = 0; i < operand_array.length; i++) {
				operand_array[i]->remove_child(this);
				core::free(*operand_array[i]);
				if (operand_array[i]->reference_count == 0)
					core::free(operand_array[i]);
			}
			core::free(operand_array);
			return;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
				if (operands[i] == NULL) break;
				operands[i]->remove_child(this);
				core::free(*operands[i]);
				if (operands[i]->reference_count == 0)
					core::free(operands[i]);
			}
			return;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.free ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline unsigned int hash(const nd_step<Formula>& key)
	{
		unsigned int hash_value = 0;
		/* TODO: precompute these and store them in a table for faster access */
		unsigned int type_hash = default_hash<nd_step_type, 0x25a5dd87>(key.type);
		switch (key.type) {
		case nd_step_type::PARAMETER:
			return type_hash ^ default_hash(key.parameter);
		case nd_step_type::TERM_PARAMETER:
			return type_hash ^ Term::hash(*key.term);
		case nd_step_type::ARRAY_PARAMETER:
			return type_hash ^ default_hash(key.parameters.data, key.parameters.length);
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::FORMULA_PARAMETER:
			return type_hash ^ Formula::hash(*key.formula);
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			hash_value = type_hash ^ default_hash(key.operand_array.length);
			for (unsigned int i = 0; i < key.operand_array.length; i++)
				hash_value ^= hash(*key.operand_array[i]);
			return hash_value;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			hash_value = type_hash;
			for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
				if (key.operands[i] != NULL)
					hash_value ^= hash(*key.operands[i]);
			}
			return hash_value;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.hash ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline bool clone(const nd_step<Formula>* src, nd_step<Formula>*& dst,
			array_map<const nd_step<Formula>*, nd_step<Formula>*>& pointer_map,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		unsigned int index = pointer_map.index_of(src);
		if (index < pointer_map.size) {
			dst = pointer_map.values[index];
			dst->reference_count++;
		} else {
			dst = (nd_step<Formula>*) malloc(sizeof(nd_step<Formula>));
			if (dst == nullptr) {
				fprintf(stderr, "nd_step.clone ERROR: Out of memory.\n");
				return false;
			} else if (!clone(*src, *dst, pointer_map, formula_map)) {
				core::free(dst);
				return false;
			} else if (!pointer_map.ensure_capacity(pointer_map.size + 1)) {
				core::free(*dst); core::free(dst);
				return false;
			}
			pointer_map.keys[pointer_map.size] = src;
			pointer_map.values[pointer_map.size] = dst;
			pointer_map.size++;
		}
		return true;
	}

	static inline bool clone(const nd_step<Formula>& src, nd_step<Formula>& dst,
			array_map<const nd_step<Formula>*, nd_step<Formula>*>& pointer_map,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		dst.type = src.type;
		dst.reference_count = 1;
		if (!array_init(dst.children, max((size_t) 1, src.children.length)))
			return false;
		switch (src.type) {
		case nd_step_type::PARAMETER:
			dst.parameter = src.parameter;
			return true;
		case nd_step_type::ARRAY_PARAMETER:
			if (!array_init(dst.parameters, src.parameters.length)) {
				core::free(dst.children);
				return false;
			}
			for (unsigned int i = 0; i < src.parameters.length; i++)
				dst.parameters[i] = src.parameters[i];
			dst.parameters.length = src.parameters.length;
			return true;
		case nd_step_type::TERM_PARAMETER:
			return ::clone(src.term, dst.term, formula_map);
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::FORMULA_PARAMETER:
			return ::clone(src.formula, dst.formula, formula_map);
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			if (!array_init(dst.operand_array, src.operand_array.length)) {
				core::free(dst.children);
				return false;
			}
			for (unsigned int i = 0; i < src.operand_array.length; i++) {
				if (!clone(src.operand_array[i], dst.operand_array[i], pointer_map, formula_map)) {
					for (unsigned int j = 0; j < i; j++) {
						core::free(*dst.operand_array[j]);
						if (dst.operand_array[j]->reference_count == 0)
							core::free(dst.operand_array[j]);
					}
					core::free(dst.children); return false;
				}
				dst.operand_array[i]->children.add(&dst);
			}
			dst.operand_array.length = src.operand_array.length;
			return true;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
				if (src.operands[i] == nullptr) {
					dst.operands[i] = nullptr;
				} else {
					if (!clone(src.operands[i], dst.operands[i], pointer_map, formula_map)) {
						for (unsigned int j = 0; j < i; j++) {
							if (dst.operands[i] != nullptr) {
								core::free(*dst.operands[j]);
								if (dst.operands[j]->reference_count == 0)
									core::free(dst.operands[j]);
							}
						}
						core::free(dst.children); return false;
					}
					dst.operands[i]->children.add(&dst);
				}
			}
			return true;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.clone ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void move(const nd_step<Formula>& src, nd_step<Formula>& dst) {
		dst.type = src.type;
		dst.reference_count = src.reference_count;
		core::move(src.children, dst.children);
		switch (src.type) {
		case nd_step_type::PARAMETER:
			dst.parameter = src.parameter; return;
		case nd_step_type::ARRAY_PARAMETER:
			core::move(src.parameters, dst.parameters); return;
		case nd_step_type::TERM_PARAMETER:
			dst.term = src.term; return;
		case nd_step_type::AXIOM:
		case nd_step_type::COMPARISON_INTRODUCTION:
		case nd_step_type::INEQUALITY_INTRODUCTION:
		case nd_step_type::FORMULA_PARAMETER:
			dst.formula = src.formula; return;
		case nd_step_type::CONJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_ELIMINATION:
			core::move(src.operand_array, dst.operand_array); return;
		case nd_step_type::BETA_EQUIVALENCE:
		case nd_step_type::CONJUNCTION_ELIMINATION:
		case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		case nd_step_type::DISJUNCTION_INTRODUCTION:
		case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
		case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
			for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++)
				dst.operands[i] = src.operands[i];
			return;
		case nd_step_type::COUNT:
			break;
		}
		fprintf(stderr, "nd_step.move ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void swap(nd_step<Formula>& first, nd_step<Formula>& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(nd_step<Formula>); i++)
			core::swap(first_data[i], second_data[i]);
	}
};

template<typename Formula>
bool get_proof_map(const nd_step<Formula>* proof,
		hash_map<const nd_step<Formula>*, unsigned int>& pointer_map,
		hash_map<const Formula*, unsigned int>& formula_map)
{
	if (!pointer_map.check_size())
		return false;
	bool contains; unsigned int index;
	pointer_map.get(proof, contains, index);
	if (contains) return true;
	pointer_map.table.keys[index] = proof;
	pointer_map.values[index] = pointer_map.table.size;
	pointer_map.table.size++;

	switch (proof->type) {
	case nd_step_type::PARAMETER:
	case nd_step_type::ARRAY_PARAMETER:
		return true;
	case nd_step_type::TERM_PARAMETER:
		return get_formula_map(proof->term, formula_map);
	case nd_step_type::AXIOM:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::FORMULA_PARAMETER:
		return get_formula_map(proof->formula, formula_map);
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_ELIMINATION:
		for (const nd_step<Formula>* operand : proof->operand_array)
			if (!get_proof_map(operand, pointer_map, formula_map)) return false;
		return true;
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
		for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
			if (proof->operands[i] == nullptr) break;
			if (!get_proof_map(proof->operands[i], pointer_map, formula_map))
				return false;
		}
		return true;
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "get_proof_map ERROR: Unrecognized nd_step_type.\n");
	return false;
}

template<typename Formula>
inline bool init_proof_array(nd_step<Formula>**& proofs, size_t length)
{
	proofs = (nd_step<Formula>**) malloc(max(1, sizeof(nd_step<Formula>*) * length));
	if (proofs == nullptr) {
		fprintf(stderr, "init_proof_array ERROR: Insufficient memory for `proofs` array.\n");
		return false;
	}

	for (size_t i = 0; i < length; i++) {
		if (!new_nd_step(proofs[i], nd_step_type::PARAMETER)) {
			fprintf(stderr, "init_proof_array ERROR: Insufficient memory for `proofs` array.\n");
			for (size_t j = 0; j < i; j++) free(proofs[j]);
			free(proofs); return false;
		}
		proofs[i]->reference_count = 1;
	}
	return true;
}

/* NOTE: this function assumes that every element in `proofs` has the fields
   `reference_count` and `children` already initialized, and that every element
   in `formulas` has the field `reference_count` initialized */
template<typename Formula, typename Stream>
bool read(nd_step<Formula>& proof, Stream& in,
		nd_step<Formula>** proofs, Formula** formulas)
{
	nd_step_type_specifier type;
	if (!read(type, in))
		return false;
	proof.type = (nd_step_type) type;

	unsigned int index;
	switch (proof.type) {
	case nd_step_type::PARAMETER:
		return read(proof.parameter, in);
	case nd_step_type::ARRAY_PARAMETER:
		if (!read(proof.parameters, in)) {
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		return true;
	case nd_step_type::TERM_PARAMETER:
		if (!read(index, in)) {
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		proof.term = formulas[index];
		formulas[index]->reference_count++;
		return true;
	case nd_step_type::AXIOM:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::FORMULA_PARAMETER:
		if (!read(index, in)) {
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		proof.formula = formulas[index];
		formulas[index]->reference_count++;
		return true;
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (!read(proof.operand_array.length, in)) {
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		proof.operand_array.data = (nd_step<Formula>**) malloc(max(1, sizeof(nd_step<Formula>*) * proof.operand_array.length));
		if (proof.operand_array.data == nullptr) {
			fprintf(stderr, "read ERROR: Insufficient memory for `nd_step.operand_array.data`.\n");
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		for (size_t i = 0; i < proof.operand_array.length; i++) {
			unsigned int index;
			if (!read(index, in)
			 || !proofs[index]->children.add(&proof))
			{
				proof.operand_array.length = i;
				return false;
			}
			proof.operand_array[i] = proofs[index];
			proofs[index]->reference_count++;
		}
		return true;
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
		if (!read(index, in)) {
			/* this is to make sure that when this proof is freed, no field other than `children` is freed */
			proof.type = nd_step_type::PARAMETER; return false;
		}
		for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++)
			proof.operands[i] = nullptr;
		for (unsigned int i = 0; i < index; i++) {
			unsigned int ind;
			if (!read(ind, in)
			 || !proofs[ind]->children.add(&proof))
			{
				return false;
			}
			proof.operands[i] = proofs[ind];
			proofs[ind]->reference_count++;
		}
		return true;
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "read ERROR: Unrecognized nd_step_type.\n");
	return false;
}

template<typename Formula, typename Stream>
bool write(const nd_step<Formula>& proof, Stream& out,
		const hash_map<const nd_step<Formula>*, unsigned int>& pointer_map,
		const hash_map<const Formula*, unsigned int>& formula_map)
{
	if (!write((nd_step_type_specifier) proof.type, out))
		return false;

	unsigned int operand_count;
	switch (proof.type) {
	case nd_step_type::PARAMETER:
		return write(proof.parameter, out);
	case nd_step_type::ARRAY_PARAMETER:
		return write(proof.parameters, out);
	case nd_step_type::TERM_PARAMETER:
		return write(formula_map.get(proof.term), out);
	case nd_step_type::AXIOM:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::FORMULA_PARAMETER:
		return write(formula_map.get(proof.formula), out);
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (!write(proof.operand_array.length, out))
			return false;
		for (size_t i = 0; i < proof.operand_array.length; i++)
			if (!write(pointer_map.get(proof.operand_array[i]), out)) return false;
		return true;
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
		operand_count = 0;
		while (operand_count < ND_OPERAND_COUNT && proof.operands[operand_count] != nullptr) operand_count++;
		if (!write(operand_count, out)) return false;
		for (unsigned int i = 0; i < operand_count; i++)
			if (!write(pointer_map.get(proof.operands[i]), out)) return false;
		return true;
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "write ERROR: Unrecognized nd_step_type.\n");
	return false;
}

template<typename Formula>
inline int_fast8_t compare(const nd_step<Formula>& first, const nd_step<Formula>& second)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;

	int_fast8_t result;
	switch (first.type) {
	case nd_step_type::PARAMETER:
		if (first.parameter < second.parameter) return -1;
		else if (first.parameter > second.parameter) return 1;
		else return 0;
	case nd_step_type::TERM_PARAMETER:
		return compare(first.term, second.term);
	case nd_step_type::ARRAY_PARAMETER:
		if (first.parameters.length < second.parameters.length) return -1;
		else if (first.parameters.length > second.parameters.length) return 1;
		for (unsigned int i = 0; i < first.parameters.length; i++) {
			if (first.parameters[i] < second.parameters[i]) return -1;
			else if (first.parameters[i] > second.parameters[i]) return 1;
		}
		return 0;
	case nd_step_type::AXIOM:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::FORMULA_PARAMETER:
		return compare(*first.formula, *second.formula);
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (first.operand_array.length < second.operand_array.length) return -1;
		else if (first.operand_array.length > second.operand_array.length) return 1;
		for (unsigned int i = 0; i < first.operand_array.length; i++) {
			int_fast8_t result = compare(*first.operand_array[i], *second.operand_array[i]);
			if (result != 0) return result;
		}
		return 0;
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
		for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
			if (first.operands[i] == NULL) {
				if (second.operands[i] == NULL)
					continue;
				else return -1;
			} else {
				if (second.operands[i] == NULL)
					return 1;
				result = compare(*first.operands[i], *second.operands[i]);
				if (result != 0) return result;
			}
		}
		return 0;
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "compare ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula>
inline bool operator < (const nd_step<Formula>& first, const nd_step<Formula>& second) {
	return compare(first, second) < 0;
}

template<typename Formula>
inline bool operator > (const nd_step<Formula>& first, const nd_step<Formula>& second) {
	return compare(first, second) > 0;
}

template<typename Formula>
inline bool operator == (const nd_step<Formula>& first, const nd_step<Formula>& second)
{
	if (first.type != second.type) return false;

	switch (first.type) {
	case nd_step_type::PARAMETER:
		return first.parameter == second.parameter;
	case nd_step_type::TERM_PARAMETER:
		return *first.term == *second.term;
	case nd_step_type::ARRAY_PARAMETER:
		if (first.parameters.length != second.parameters.length) return false;
		for (unsigned int i = 0; i < first.parameters.length; i++)
			if (first.parameters[i] != second.parameters[i]) return false;
		return true;
	case nd_step_type::AXIOM:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::FORMULA_PARAMETER:
		return *first.formula == *second.formula;
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (first.operand_array.length != second.operand_array.length) return false;
		for (unsigned int i = 0; i < first.operand_array.length; i++)
			if (*first.operand_array[i] != *second.operand_array[i]) return false;
		return true;
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
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
		for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
			if (first.operands[i] == NULL) {
				if (second.operands[i] == NULL)
					continue;
				else return false;
			} else {
				if (second.operands[i] == NULL)
					return false;
				else if (*first.operands[i] != *second.operands[i])
					return false;
			}
		}
		return true;
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula>
inline bool operator != (const nd_step<Formula>& first, const nd_step<Formula>& second) {
	return !(first == second);
}

template<typename Formula>
struct proof_state {
	array<const Formula*> assumptions;
	Formula* formula;

	proof_state() : assumptions(8) { }
	~proof_state() { free(); }

	inline unsigned int index_of_assumption(const Formula* assumption) const {
		for (unsigned int i = 0; i < assumptions.length; i++)
			if (*assumptions[i] == *assumption) return i;
		return assumptions.length;
	}

	inline bool assumptions_have_parameter(unsigned int parameter) const {
		for (const Formula* assumption : assumptions)
			if (contains_parameter(*assumption, parameter)) return true;
		return false;
	}

	static inline void move(const proof_state<Formula>& src, proof_state<Formula>& dst) {
		core::move(src.assumptions, dst.assumptions);
		core::move(src.formula, dst.formula);
	}

	static inline void free(proof_state<Formula>& state) {
		state.free();
		core::free(state.assumptions);
	}

private:
	inline void free() {
		if (formula != NULL) {
			core::free(*formula);
			if (formula->reference_count == 0)
				core::free(formula);
		}
	}
};

template<typename Formula>
inline bool init(proof_state<Formula>& state) {
	state.formula = NULL;
	return array_init(state.assumptions, 8);
}

template<typename K, typename Formula>
inline void free_proof_states(hash_map<K, proof_state<Formula>>& map) {
	for (auto entry : map)
		free(entry.value);
}

template<typename K, typename Formula>
inline void free_proof_states(hash_map<K, pair<proof_state<Formula>, unsigned int>>& map) {
	for (auto entry : map)
		free(entry.value.key);
}

template<typename FormulaPtr>
struct proof_state_formulas {
	typedef typename std::remove_pointer<FormulaPtr>::type Formula;

	proof_state<Formula>** states;
	unsigned int length;

	proof_state_formulas(proof_state<Formula>** states, unsigned int length) : states(states), length(length) { }

	inline FormulaPtr operator[] (size_t index) const {
		return states[index]->formula;
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename Formula>
struct proof_state_assumptions {
	proof_state<Formula>** states;
	unsigned int length;

	proof_state_assumptions(proof_state<Formula>** states, unsigned int length) : states(states), length(length) { }

	inline const array<const Formula*>& operator[] (size_t index) const {
		return states[index]->assumptions;
	}
};

template<typename Formula>
proof_state_assumptions<Formula> make_proof_state_assumptions(proof_state<Formula>** states, unsigned int length) {
	return proof_state_assumptions<Formula>(states, length);
}

template<typename T, typename ExcludedType>
struct exclude {
	T* elements;
	unsigned int length;
	ExcludedType excluded;

	exclude(T* elements, unsigned int length, ExcludedType& excluded) :
			elements(elements), length(length), excluded(excluded) { }

	inline void operator = (const exclude<T, ExcludedType>& src) {
		elements = src.elements;
		length = src.length;
		excluded = src.excluded;
	}
};

template<typename T, typename ExcludedType>
inline exclude<T, ExcludedType> make_exclude(T* elements, unsigned int length, ExcludedType& excluded) {
	return exclude<T, ExcludedType>(elements, length, excluded);
}

template<typename T>
bool pass_hypotheses(array<T>& dst, const array<T>& first) {
	return dst.append(first.data, first.length);
}

template<typename T>
bool pass_hypotheses(array<T>& dst, const array<T>& first, const array<T>& second) {
	return set_union(dst, first, second);
}

template<typename T, typename ExcludedType>
bool pass_hypotheses(array<T>& dst, const exclude<T, ExcludedType>& first)
{
	if (!dst.ensure_capacity(dst.length + first.length))
		return false;

	for (unsigned int i = 0; i < first.length; i++) {
		if (*first.elements[i] == *first.excluded) { continue; }
		dst[dst.length] = first.elements[i];
		dst.length++;
	}
	return true;
}

template<typename T, typename ExcludedType, bool RemoveDuplicates = true>
bool pass_hypotheses(array<T>& dst, const array<T>& first, const exclude<T, ExcludedType>& second)
{
	if (!dst.ensure_capacity(first.length + second.length))
		return false;

	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		if (*second.excluded == *second.elements[j]) { j++; continue; }

		if (first[i] == second.elements[j]) {
			set_union_helper<RemoveDuplicates>(dst.data, dst.length, first[i]);
			i++; j++;
		} else if (first[i] < second.elements[j]) {
			set_union_helper<RemoveDuplicates>(dst.data, dst.length, first[i]);
			i++;
		} else {
			set_union_helper<RemoveDuplicates>(dst.data, dst.length, second.elements[j]);
			j++;
		}
	}

	while (i < first.length) {
		set_union_helper<RemoveDuplicates>(dst.data, dst.length, first[i]);
		i++;
	} while (j < second.length) {
		if (*second.excluded == *second.elements[j]) { j++; continue; }
		set_union_helper<RemoveDuplicates>(dst.data, dst.length, second.elements[j]);
		j++;
	}

	return true;
}

template<typename T, typename... Lists>
bool pass_hypotheses(array<T>& out, const array<T>& first, Lists&&... lists) {
	if (!pass_hypotheses(out, std::forward<Lists>(lists)...)) return false;
	array<T> temp = array<T>(out.length + first.length);
	swap(temp, out);
	return set_union(out, temp, first);
}

template<typename T, typename ExcludedType, typename... Lists>
bool pass_hypotheses(array<T>& out, const exclude<T, ExcludedType>& first, Lists&&... lists) {
	if (!pass_hypotheses(out, std::forward<Lists>(lists)...)) return false;
	array<T> temp = array<T>(out.length + first.length);
	swap(temp, out);
	return pass_hypotheses(out, temp, first);
}

template<typename T, typename Arrays>
bool pass_hypotheses(array<T>& out, const Arrays& lists) {
	if (lists.length == 0)
		return true;
	if (!pass_hypotheses(out, lists[0])) return false;

	array<T> temp = array<T>(16);
	for (unsigned int i = 1; i < lists.length; i++) {
		if (!temp.ensure_capacity(out.length + lists[i].length))
			return false;
		swap(temp, out);
		if (!pass_hypotheses(out, temp, lists[i]))
			return false;
		temp.clear();
	}
	return true;
}

template<typename Formula>
bool subtract_unifying_hypotheses(
	array<const Formula*>& dst_hypotheses,
	const array<const Formula*>& hypotheses,
	const Formula& src, const Formula& C)
{
	typedef typename Formula::Term Term;

	for (const Formula* hypothesis : hypotheses) {
		unsigned int parameter;
		Term* variable = Formula::new_variable(src.quantifier.variable);
		if (unifies_parameter(*src.quantifier.operand, *hypothesis, variable, parameter)) {
			free(*variable); if (variable->reference_count == 0) free(variable);
			if (contains_parameter(C, parameter)) {
				if (!dst_hypotheses.add(hypothesis)) return false;
			}
		} else {
			free(*variable); if (variable->reference_count == 0) free(variable);
			if (!dst_hypotheses.add(hypothesis)) return false;
		}
	}
	return true;
}

template<typename Canonicalizer, typename Formula>
inline Formula* try_canonicalize(Formula* formula) {
	Formula* out = Canonicalizer::canonicalize(*formula);
#if !defined(NDEBUG)
	if (out == nullptr) {
		print("try_canonicalize ERROR: Failed to canonicalize ", stderr);
		print(*formula, stderr, *debug_terminal_printer); print('\n', stderr);
	}
#endif
	free(*formula); if (formula->reference_count == 0) free(formula);
	return out;
}

template<typename Formula>
struct proof_substitution {
	array_map<nd_step<Formula>*, nd_step<Formula>*> map;

	static inline void free(proof_substitution<Formula>& s) {
		for (auto entry : s.map) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}
		core::free(s.map);
	}
};

template<typename Formula>
inline bool init(proof_substitution<Formula>& substitution) {
	return array_map_init(substitution.map, 8);
}

template<typename Formula>
inline nd_step<Formula>* map(
		nd_step<Formula>* const step)
{
	return step;
}

template<typename Formula>
inline const nd_step<Formula>* map_const(
		const nd_step<Formula>* const step)
{
	return step;
}

template<typename Formula, typename ProofStep,
	typename std::enable_if<!std::is_const<typename std::remove_pointer<ProofStep>::type>::value>::type* = nullptr>
inline nd_step<Formula>* map(ProofStep step,
		const proof_substitution<Formula>& substitution)
{
	bool contains;
	auto value = substitution.map.get(step, contains);
	nd_step<Formula>* new_step = value;
	if (contains) return new_step;
	else return step;
}

template<typename Formula, typename ProofStep>
inline const nd_step<Formula>* map_const(
		const ProofStep& step,
		const proof_substitution<Formula>& substitution)
{
	bool contains;
	auto value = substitution.map.get(step, contains);
	const nd_step<Formula>* new_step = value;
	if (contains) return new_step;
	else return step;
}

template<typename BuiltInPredicates, typename Canonicalizer,
	bool Intuitionistic, typename Formula, typename... ProofMap>
bool check_proof(proof_state<Formula>& out,
		const nd_step<Formula>& proof,
		proof_state<Formula>** operand_states,
		unsigned int operand_count,
		ProofMap&&... proof_map)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	Formula* formula;
	const nd_step<Formula>* second_operand;
	const nd_step<Formula>* third_operand;
	Term* parameter; Term* variable;
	array<unsigned int> constants = array<unsigned int>(8);
	exclude<const Formula*, Formula*>* assumptions;
	switch (proof.type) {
	case nd_step_type::AXIOM:
		if (!is_canonical<Canonicalizer>(*proof.formula)) {
			print("check_proof ERROR: Axiom '", stderr);
			print(*proof.formula, stderr, *debug_terminal_printer);
			print("' is not in canonical form.\n", stderr);
			return false;
		}
		out.formula = proof.formula;
		out.formula->reference_count++;
		return out.assumptions.add(out.formula);
	case nd_step_type::COMPARISON_INTRODUCTION:
		if (proof.formula->type == FormulaType::NOT) {
			Formula* operand = proof.formula->unary.operand;
			if (operand->type != FormulaType::BINARY_APPLICATION
			 || operand->ternary.first->type != TermType::CONSTANT
			 || operand->ternary.first->constant != (unsigned int) BuiltInPredicates::GREATER_THAN_OR_EQUAL
			 || operand->ternary.second->type != TermType::NUMBER
			 || operand->ternary.third->type != TermType::NUMBER
			 || operand->ternary.third->number <= operand->ternary.second->number)
				return false;
		} else {
			Formula* operand = proof.formula;
			if (operand->type != FormulaType::BINARY_APPLICATION
			 || operand->ternary.first->type != TermType::CONSTANT
			 || operand->ternary.first->constant != (unsigned int) BuiltInPredicates::GREATER_THAN_OR_EQUAL
			 || operand->ternary.second->type != TermType::NUMBER
			 || operand->ternary.third->type != TermType::NUMBER
			 || operand->ternary.second->number < operand->ternary.third->number)
				return false;
		}
		out.formula = proof.formula;
		out.formula->reference_count++;
		return true;
	case nd_step_type::INEQUALITY_INTRODUCTION:
		if (proof.formula->type != FormulaType::NOT || proof.formula->unary.operand->type != FormulaType::EQUALS)
			return false;
		formula = proof.formula->unary.operand;
		if (formula->binary.right->type != TermType::CONSTANT && formula->binary.right->type != TermType::NUMBER && formula->binary.right->type != TermType::STRING)
			return false;
		if (formula->binary.left->type == TermType::CONSTANT) {
			if (formula->binary.right->type == TermType::CONSTANT && formula->binary.left->constant == formula->binary.right->constant)
				return false;
		} else if (formula->binary.left->type == TermType::NUMBER) {
			if (formula->binary.right->type == TermType::NUMBER && formula->binary.left->number == formula->binary.right->number)
				return false;
		} else if (formula->binary.left->type == TermType::STRING) {
			if (formula->binary.right->type == TermType::STRING && formula->binary.left->str == formula->binary.right->str)
				return false;
		} else {
			return false;
		}
		out.formula = proof.formula;
		out.formula->reference_count++;
		return true;
	case nd_step_type::BETA_EQUIVALENCE:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || proof.operands[0]->type != nd_step_type::FORMULA_PARAMETER || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		formula = beta_reduce(proof.operands[0]->formula);
		if (formula != nullptr) {
			Formula* relabeled = relabel_variables(formula);
			free(*formula); if (formula->reference_count == 0) free(formula);
			if (relabeled == nullptr)
				return false;
			formula = relabeled;
		} else {
			return false;
		}
		out.formula = beta_reduce(second_operand->formula);
		if (out.formula != nullptr) {
			Formula* relabeled = relabel_variables(out.formula);
			free(*out.formula); if (out.formula->reference_count == 0) free(out.formula);
			if (relabeled == nullptr) {
				free(*formula); if (formula->reference_count == 0) free(formula);
				return false;
			}
			out.formula = relabeled;
		} else {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		if (*formula != *out.formula) {
			fprintf(stderr, "check_proof ERROR: The two formula are not beta-equivalent.\n");
			free(*formula); if (formula->reference_count == 0) free(formula);
			free(*out.formula); if (out.formula->reference_count == 0) free(out.formula);
			out.formula = NULL;
			return false;
		}
		free(*formula); if (formula->reference_count == 0) free(formula);
		free(*out.formula); if (out.formula->reference_count == 0) free(out.formula);

		out.formula = Formula::new_equals(proof.operands[0]->formula, second_operand->formula);
		if (out.formula == NULL) return false;
		proof.operands[0]->formula->reference_count++;
		second_operand->formula->reference_count++;
		return true;
	case nd_step_type::CONJUNCTION_INTRODUCTION:
		if (operand_count < 2) return false;
		out.formula = Formula::new_and(proof_state_formulas<Formula*>(operand_states, operand_count));
		if (out.formula == NULL) return false;
		for (unsigned int i = 0; i < operand_count; i++)
			operand_states[i]->formula->reference_count++;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, make_proof_state_assumptions(operand_states, operand_count));
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if ((proof.type == nd_step_type::CONJUNCTION_ELIMINATION && (operand_count != 2 || second_operand->type != nd_step_type::ARRAY_PARAMETER))
		 || (proof.type != nd_step_type::CONJUNCTION_ELIMINATION && operand_count != 1))
			return false;
		formula = operand_states[0]->formula;
		if (formula == NULL) {
			return false;
		} else if (formula->type != FormulaType::AND || formula->array.length < 2) {
			fprintf(stderr, "check_proof ERROR: Expected a conjunction.\n");
			return false;
		}
		if (proof.type == nd_step_type::CONJUNCTION_ELIMINATION_LEFT) {
			out.formula = formula->array.operands[0];
			out.formula->reference_count++;
		} else if (proof.type == nd_step_type::CONJUNCTION_ELIMINATION_RIGHT) {
			out.formula = formula->array.operands[1];
			out.formula->reference_count++;
		} else if (proof.type == nd_step_type::CONJUNCTION_ELIMINATION) {
			const indexed_array_view<Formula*> conjuncts(formula->array.operands, second_operand->parameters.data, second_operand->parameters.length);
			if (*std::max_element(conjuncts.indices, conjuncts.indices + conjuncts.length) >= formula->array.length
			 || !std::is_sorted(conjuncts.indices, conjuncts.indices + conjuncts.length)
			 || std::adjacent_find(conjuncts.indices, conjuncts.indices + conjuncts.length) != conjuncts.indices + conjuncts.length)
				return false;
			if (conjuncts.length == 1) {
				out.formula = conjuncts[0];
			} else {
				out.formula = Formula::new_and(conjuncts);
				if (out.formula == NULL) return false;
			}
			for (unsigned int i = 0; i < conjuncts.length; i++)
				conjuncts[i]->reference_count++;
		}
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_INTRODUCTION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		third_operand = map_const(proof.operands[2], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 3 || second_operand->type != nd_step_type::FORMULA_PARAMETER || third_operand->type != nd_step_type::PARAMETER)
			return false;
		if (second_operand->formula->type == FormulaType::OR) {
			out.formula = Formula::new_or(make_included_array_view(
					second_operand->formula->array.operands,
					second_operand->formula->array.length,
					operand_states[0]->formula, third_operand->parameter));
		} else if (third_operand->parameter == 0) {
			out.formula = Formula::new_or(operand_states[0]->formula, second_operand->formula);
		} else {
			out.formula = Formula::new_or(second_operand->formula, operand_states[0]->formula);
		}
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		second_operand->formula->reference_count++;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 3 || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_LEFT)
			out.formula = Formula::new_or(operand_states[0]->formula, second_operand->formula);
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT)
			out.formula = Formula::new_or(second_operand->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		second_operand->formula->reference_count++;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (operand_count < 3 || operand_states[0]->formula->type != FormulaType::OR)
			return false;
		for (unsigned int i = 2; i < operand_count; i++)
			if (*operand_states[i]->formula != *operand_states[1]->formula) return false;
		out.formula = operand_states[1]->formula;
		out.formula->reference_count++;
		assumptions = (exclude<const Formula*, Formula*>*)
				malloc(sizeof(exclude<const Formula*, Formula*>) * (operand_count - 1));
		for (unsigned int i = 1; i < operand_count; i++) {
			assumptions[i - 1] = make_exclude(operand_states[i]->assumptions.data,
					operand_states[i]->assumptions.length, operand_states[0]->formula->array.operands[i - 1]);
		}
		if (!pass_hypotheses(out.assumptions, operand_states[0]->assumptions, make_array_view(assumptions, operand_count - 1))) {
			free(assumptions);
			return false;
		}
		free(assumptions);
		return true;
	case nd_step_type::IMPLICATION_INTRODUCTION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || second_operand->type != nd_step_type::AXIOM)
			return false;
		out.formula = Formula::new_if_then(operand_states[1]->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		operand_states[1]->formula->reference_count++;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions,
				make_exclude(operand_states[0]->assumptions.data, operand_states[0]->assumptions.length, operand_states[1]->formula));
	case nd_step_type::IMPLICATION_ELIMINATION:
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::IF_THEN
		 || *operand_states[0]->formula->binary.left != *operand_states[1]->formula)
			return false;
		out.formula = operand_states[0]->formula->binary.right;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::BICONDITIONAL_INTRODUCTION:
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::IF_THEN
		 || operand_states[1]->formula->type != FormulaType::IF_THEN
		 || *operand_states[0]->formula->binary.left == *operand_states[1]->formula->binary.right
		 || *operand_states[0]->formula->binary.right == *operand_states[1]->formula->binary.left)
			return false;
		out.formula = Formula::new_iff(operand_states[0]->formula->binary.left, operand_states[0]->formula->binary.right);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->binary.left->reference_count++;
		operand_states[0]->formula->binary.right->reference_count++;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT:
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::IFF
		 || *operand_states[0]->formula->binary.left != *operand_states[1]->formula)
			return false;
		out.formula = operand_states[0]->formula->binary.right;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT:
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::IFF
		 || *operand_states[0]->formula->binary.right != *operand_states[1]->formula)
			return false;
		out.formula = operand_states[0]->formula->binary.left;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::PROOF_BY_CONTRADICTION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FALSE
		 || second_operand->type != nd_step_type::AXIOM)
			return false;
		if (!Intuitionistic && operand_states[1]->formula->type == FormulaType::NOT) {
			out.formula = second_operand->formula->unary.operand;
			out.formula->reference_count++;
		} else {
			out.formula = Formula::new_not(second_operand->formula);
			second_operand->formula->reference_count++;
		}
		return pass_hypotheses(out.assumptions,
				make_exclude(operand_states[0]->assumptions.data, operand_states[0]->assumptions.length, second_operand->formula));
	case nd_step_type::NEGATION_ELIMINATION:
		if (operand_count != 2) {
			return false;
		} else if (operand_states[1]->formula->type == FormulaType::NOT
				&& *operand_states[0]->formula == *operand_states[1]->formula->unary.operand)
		{
			out.formula = Formula::new_false();
		} else if (operand_states[0]->formula->type == FormulaType::NOT
				&& *operand_states[1]->formula == *operand_states[0]->formula->unary.operand)
		{
			out.formula = Formula::new_false();
		} else { return false; }
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::FALSITY_ELIMINATION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FALSE
		 || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		out.formula = second_operand->formula;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::UNIVERSAL_INTRODUCTION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || operand_states[0]->formula == NULL || second_operand->type != nd_step_type::PARAMETER) {
			return false;
		} else if (operand_states[0]->assumptions_have_parameter(second_operand->parameter)) {
			/* the parameter is not allowed to occur free in the assumptions */
			return false;
		}

		parameter = Formula::new_parameter(second_operand->parameter);
		formula = substitute<TermType::PARAMETER, 1>(operand_states[0]->formula, parameter, &Term::template variables<1>::value);
		free(*parameter); if (parameter->reference_count == 0) free(parameter);
		if (formula == NULL) return false;
		out.formula = Formula::new_for_all(1, formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::UNIVERSAL_ELIMINATION:
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FOR_ALL
		 || second_operand->type != nd_step_type::TERM_PARAMETER)
			return false;
		variable = Formula::new_variable(operand_states[0]->formula->quantifier.variable);
		{
			unsigned int max_variable = 0;
			Formula* to_substitute = second_operand->term;
			if (max_bound_variable(*operand_states[0]->formula, max_variable)) {
				to_substitute = shift_bound_variables(to_substitute, max_variable - 1);
			} else {
				to_substitute->reference_count++;
			}
			Formula* substituted = substitute<TermType::VARIABLE, -1>(operand_states[0]->formula->quantifier.operand, variable, to_substitute);
			free(*to_substitute); if (to_substitute->reference_count == 0) free(to_substitute);
			out.formula = relabel_variables(substituted);
			free(*substituted); if (substituted->reference_count == 0) free(substituted);
		}
		free(*variable); if (variable->reference_count == 0) free(variable);
		if (out.formula == NULL) return false;
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		if (operand_count != 3) return false;
		second_operand = map_const(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (second_operand->type == nd_step_type::ARRAY_PARAMETER && proof.operands[2]->type == nd_step_type::TERM_PARAMETER) {
			Term* replaced_term = get_term_at_index(*operand_states[0]->formula, second_operand->parameters.data[0]);
			if (replaced_term == nullptr || *replaced_term != *proof.operands[2]->term)
				return false;

			Formula* temp = shift_bound_variables(operand_states[0]->formula, 1);
			if (temp == NULL) return false;

			formula = substitute(temp, second_operand->parameters.data, second_operand->parameters.length, &Term::template variables<1>::value);
			free(*temp); if (temp->reference_count == 0) free(temp);
		} else {
			return false;
		}

		if (formula == NULL) {
			fprintf(stderr, "check_proof ERROR: Invalid term indices for substitution.\n");
			return false;
		}
		out.formula = Formula::new_exists(1, formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::EXISTENTIAL_ELIMINATION:
		if (operand_count != 2 || operand_states[0]->formula->type != FormulaType::EXISTS)
			return false;
		if (!subtract_unifying_hypotheses(out.assumptions, operand_states[1]->assumptions, *operand_states[0]->formula, *operand_states[1]->formula))
			return false;
		out.formula = operand_states[1]->formula;
		out.formula->reference_count++;
		return true;
	case nd_step_type::EQUALITY_ELIMINATION:
		if (operand_count != 3 || operand_states[0]->formula->type != FormulaType::EQUALS) return false;
		if (proof.operands[2]->type == nd_step_type::ARRAY_PARAMETER) {
			out.formula = substitute(operand_states[1]->formula,
					proof.operands[2]->parameters.data, proof.operands[2]->parameters.length,
					operand_states[0]->formula->binary.left, operand_states[0]->formula->binary.right);
			if (out.formula == NULL) {
				out.formula = substitute(operand_states[1]->formula,
						proof.operands[2]->parameters.data, proof.operands[2]->parameters.length,
						operand_states[0]->formula->binary.right, operand_states[0]->formula->binary.left);
				if (out.formula == NULL) return false;
			}
		} else {
			return false;
		}

		out.formula = try_canonicalize<Canonicalizer>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::PARAMETER:
	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "check_proof ERROR: Unexpected nd_step_type.\n");
	return false;
}

template<bool RevisitVisitedNodes = false, typename Formula, typename Visitor>
bool visit(const nd_step<Formula>* proof, Visitor& visitor)
{
	array<const nd_step<Formula>*> stack(64);
	hash_set<const nd_step<Formula>*> visited(128);
	if (!stack.add(proof)) return false;
	while (stack.length > 0)
	{
		const nd_step<Formula>* node = stack.pop();
		if (!RevisitVisitedNodes && !visited.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);

		if (!visit_node(node, visitor)) return false;

		if (operand_count == 0) continue;
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;

			if (!RevisitVisitedNodes && visited.contains(operands[i]))
				continue;
			if (!stack.add(operands[i]))
				return false;
		}
	}
	return true;
}

template<nd_step_type Type, typename Collection>
struct proof_step_collector {
	Collection& steps;
};

template<nd_step_type Type, typename Formula, typename Collection>
inline bool visit_node(const nd_step<Formula>* proof, proof_step_collector<Type, Collection>& collector) {
	if (proof->type == Type)
		return collector.steps.add(proof);
	return true;
}

template<nd_step_type Type, typename Formula, typename Collection>
inline bool get_proof_steps(const nd_step<Formula>* proof, Collection& steps) {
	proof_step_collector<Type, Collection> collector = {steps};
	return visit<true>(proof, collector);
}

template<nd_step_type Type>
struct static_proof_step_finder { };

template<nd_step_type Type, typename Formula>
inline bool visit_node(const nd_step<Formula>* proof, static_proof_step_finder<Type>& visitor) {
	if (proof->type == Type)
		return false;
	return true;
}

template<nd_step_type Type, typename Formula>
inline bool has_proof_step(const nd_step<Formula>* proof) {
	static_proof_step_finder<Type> visitor;
	return !visit<false>(proof, visitor);
}

template<bool RevisitVisitedNodes = false, typename Formula, typename Visitor>
bool visit_with_pruning(const nd_step<Formula>* proof, Visitor& visitor)
{
	array<const nd_step<Formula>*> stack(64);
	hash_set<const nd_step<Formula>*> visited(128);
	if (!stack.add(proof)) return false;
	while (stack.length > 0)
	{
		const nd_step<Formula>* node = stack.pop();
		if (!RevisitVisitedNodes && !visited.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);

		bool recurse_visit = true;
		if (!visit_node(node, visitor, recurse_visit)) return false;

		if (operand_count == 0 || !recurse_visit) continue;
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;

			if (!RevisitVisitedNodes && visited.contains(operands[i]))
				continue;
			if (!stack.add(operands[i]))
				return false;
		}
	}
	return true;
}

template<nd_step_type Type, typename Formula>
struct most_recent_proof_step_finder {
	array<const nd_step<Formula>*>& steps;

	most_recent_proof_step_finder(array<const nd_step<Formula>*>& steps) : steps(steps) { }
};

template<nd_step_type Type, typename Formula>
inline bool visit_node(const nd_step<Formula>* proof, most_recent_proof_step_finder<Type, Formula>& visitor, bool& recurse_visit) {
	if (proof->type == Type) {
		if (!visitor.steps.add(proof))
			return false;
		recurse_visit = false;
	}
	return true;
}

template<nd_step_type Type, typename Formula>
inline bool get_most_recent_proof_steps(const nd_step<Formula>* proof, array<const nd_step<Formula>*>& steps)
{
	most_recent_proof_step_finder<Type, Formula> visitor(steps);
	return !visit_with_pruning<true>(proof, visitor);
}

template<typename Formula, typename... ProofMap>
bool compute_in_degrees(const nd_step<Formula>* proof,
		hash_map<const nd_step<Formula>*, unsigned int>& in_degrees,
		ProofMap&&... proof_map)
{
	array<const nd_step<Formula>*> stack(64);
	hash_set<const nd_step<Formula>*> visited(128);
	if (!stack.add(map_const(proof, std::forward<ProofMap>(proof_map)...)))
		return false;
	while (stack.length > 0)
	{
		const nd_step<Formula>* node = stack.pop();
		if (!visited.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		if (!in_degrees.check_size(in_degrees.table.size + operand_count + 1))
			return false;

		bool contains; unsigned int bucket;
		unsigned int& degree = in_degrees.get(node, contains, bucket);
		if (!contains) {
			in_degrees.table.keys[bucket] = node;
			in_degrees.table.size++;
			degree = 0;
		}

		if (operand_count == 0) continue;
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;

			const nd_step<Formula>* operand = map_const(
					operands[i], std::forward<ProofMap>(proof_map)...);
			unsigned int& degree = in_degrees.get(operand, contains, bucket);
			if (!contains) {
				in_degrees.table.keys[bucket] = operand;
				in_degrees.table.size++;
				degree = 1;
			} else {
				degree++;
			}

			if (visited.contains(operand))
				continue;
			if (!stack.add(operand))
				return false;
		}
	}

	return true;
}

template<typename BuiltInPredicates, typename Canonicalizer, bool Intuitionistic, typename Formula, typename... ProofMap>
Formula* compute_proof_conclusion(const nd_step<Formula>& proof, ProofMap&&... proof_map)
{
	/* first list the proof steps in reverse topological order */
	hash_map<const nd_step<Formula>*, unsigned int> in_degrees(128);
	if (!compute_in_degrees(&proof, in_degrees, std::forward<ProofMap>(proof_map)...)) return NULL;

	array<const nd_step<Formula>*> stack(32);
	for (const auto& entry : in_degrees) {
		if (entry.value == 0 && !stack.add(entry.key))
			return NULL;
	}

	array<const nd_step<Formula>*> topological_order(64);
	while (stack.length > 0) {
		const nd_step<Formula>* node = stack.pop();
		if (!node->is_parameter() && !topological_order.add(node)) return NULL;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			const nd_step<Formula>* operand = map_const(
					operands[i], std::forward<ProofMap>(proof_map)...);
			unsigned int& degree = in_degrees.get(operand);
			degree--;

			if (degree == 0 && !stack.add(operand))
				return NULL;
		}
	}

	/* process the proof steps in reverse topological order */
	hash_map<const nd_step<Formula>*, proof_state<Formula>> proof_states(128);
	array<proof_state<Formula>*> operand_states(16);
	for (unsigned int k = topological_order.length; k > 0; k--) {
		const nd_step<Formula>* node = topological_order[k - 1];
		if (!proof_states.check_size()) return NULL;

		bool contains; unsigned int bucket;
		proof_state<Formula>& state = proof_states.get(node, contains, bucket);
		if (contains) {
			fprintf(stderr, "compute_proof_conclusion ERROR: The proof state at this node should be uninitialized.\n");
			free_proof_states(proof_states); return NULL;
		} else if (!init(state)) {
			fprintf(stderr, "compute_proof_conclusion ERROR: Unable to initialize new proof_state.\n");
			free_proof_states(proof_states); return NULL;
		}
		proof_states.table.keys[bucket] = node;
		proof_states.table.size++;

		/* get the proof states of the operands */
		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		if (!operand_states.ensure_capacity(operand_count)) {
			free_proof_states(proof_states); return NULL;
		}
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) break;
			const nd_step<Formula>* operand = map_const(
					operands[i], std::forward<ProofMap>(proof_map)...);
			operand_states[i] = &proof_states.get(operand, contains);
			if (!operand->is_parameter() && !contains) {
				fprintf(stderr, "check_proof ERROR: The proof is not topologically ordered.\n");
				free_proof_states(proof_states); return NULL;
			}
			operand_states.length++;
		}

		/* check this proof step */
		if (!check_proof<BuiltInPredicates, Canonicalizer, Intuitionistic>(state, *node, operand_states.data, operand_states.length, std::forward<ProofMap>(proof_map)...)) {
			free_proof_states(proof_states);
			return NULL;
		}
		operand_states.clear();
	}

	/* get the proof state of the last deduction step */
	bool contains;
	const proof_state<Formula>& root_state = proof_states.get(
			map_const(&proof, std::forward<ProofMap>(proof_map)...), contains);
	if (!contains) {
		fprintf(stderr, "check_proof ERROR: Unable to find proof state of root.\n");
		free_proof_states(proof_states); return NULL;
	}
	Formula* formula = root_state.formula;
	formula->reference_count++;
	free_proof_states(proof_states);
	return formula;
}

template<typename BuiltInPredicates, typename Canonicalizer,
	bool Intuitionistic, typename Formula, typename... ProofMap>
bool check_proof(const nd_step<Formula>& proof,
		const Formula* expected_conclusion,
		ProofMap&&... proof_map)
{
	Formula* actual_conclusion = compute_proof_conclusion<BuiltInPredicates, Canonicalizer, Intuitionistic>(proof, std::forward<ProofMap>(proof_map)...);
	if (actual_conclusion == NULL) return false;
	bool success = (*actual_conclusion == *expected_conclusion);
/* TODO: for debugging; delete this */
//print("actual_conclusion:   ", stderr); print(*actual_conclusion, stderr, *debug_terminal_printer); print('\n', stderr);
//print("expected_conclusion: ", stderr); print(*expected_conclusion, stderr, *debug_terminal_printer); print('\n', stderr);
	if (!success)
		fprintf(stderr, "check_proof ERROR: Actual concluding formula does not match the expected formula.\n");
	free(*actual_conclusion);
	if (actual_conclusion->reference_count == 0)
		free(actual_conclusion);
	return success;
}

template<typename Formula>
bool new_nd_step(nd_step<Formula>*& step, nd_step_type type)
{
	step = (nd_step<Formula>*) malloc(sizeof(nd_step<Formula>));
	if (step == NULL) {
		fprintf(stderr, "new_nd_step ERROR: Out of memory.\n");
		return false;
	} else if (!array_init(step->children, 4)) {
		free(step); return false;
	}
	step->type = type;
	return true;
}

template<typename BuiltInPredicates, typename Canonicalizer, bool Intuitionistic, typename Formula, typename Stream, typename... Printer>
bool print(const nd_step<Formula>& proof, Stream& out, Printer&&... printer)
{
	/* first list the proof steps in reverse topological order */
	hash_map<const nd_step<Formula>*, unsigned int> in_degrees(128);
	if (!compute_in_degrees(&proof, in_degrees)) return false;

	array<const nd_step<Formula>*> stack(32);
	for (const auto& entry : in_degrees) {
		if (entry.value == 0 && !stack.add(entry.key))
			return false;
	}

	array<const nd_step<Formula>*> topological_order(64);
	while (stack.length > 0) {
		const nd_step<Formula>* node = stack.pop();
		if (!node->is_parameter() && !topological_order.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			const nd_step<Formula>* operand = operands[i];
			unsigned int& degree = in_degrees.get(operand);
			degree--;

			if (degree == 0 && !stack.add(operand))
				return false;
		}
	}

	/* process the proof steps in reverse topological order */
	hash_map<const nd_step<Formula>*, pair<proof_state<Formula>, unsigned int>> proof_states(128);
	array<proof_state<Formula>*> operand_states(16);
	array<unsigned int> operand_indices(16);
	unsigned int proof_step_index = 0;
	for (unsigned int k = topological_order.length; k > 0; k--) {
		const nd_step<Formula>* node = topological_order[k - 1];
		if (!proof_states.check_size()) return false;

		bool contains; unsigned int bucket;
		pair<proof_state<Formula>, unsigned int>& state = proof_states.get(node, contains, bucket);
		if (contains) {
			fprintf(stderr, "print ERROR: The proof state at this node should be uninitialized.\n");
			free_proof_states(proof_states); return false;
		} else if (!init(state.key)) {
			fprintf(stderr, "print ERROR: Unable to initialize new proof_state.\n");
			free_proof_states(proof_states); return false;
		}
		if (!node->is_parameter())
			state.value = ++proof_step_index;
		else state.value = 0;
		proof_states.table.keys[bucket] = node;
		proof_states.table.size++;

		/* get the proof states of the operands */
		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		if (!operand_states.ensure_capacity(operand_count)
		 || !operand_indices.ensure_capacity(operand_count))
		{
			free_proof_states(proof_states);
			return false;
		}
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) break;
			const nd_step<Formula>* operand = operands[i];
			pair<proof_state<Formula>, unsigned int>& operand_state = proof_states.get(operand, contains);
			operand_states[i] = &operand_state.key;
			if (!operand->is_parameter() && !contains) {
				fprintf(stderr, "print ERROR: The proof is not topologically ordered.\n");
				free_proof_states(proof_states); return false;
			}
			operand_states.length++;
			if (!operand->is_parameter()) {
				operand_indices[i] = operand_state.value;
				operand_indices.length++;
			}
		}

		/* compute the conclusion at this proof step */
		if (!check_proof<BuiltInPredicates, Canonicalizer, Intuitionistic>(state.key, *node, operand_states.data, operand_states.length)) {
			print("print ERROR: `check_proof` failed to compute conclusion of proof step ", stderr);
			print(state.value, stderr); print(".\n", stderr);
			free_proof_states(proof_states);
			return false;
		}

		if (!node->is_parameter()) {
			/* print this proof step */
			if (!print("  [", out) || !print(state.value, out) || !print("] ", out)
			 || !print(*state.key.formula, out, std::forward<Printer>(printer)...)
			 || !print(" by ", out) || !print<Intuitionistic>(node->type, out))
				return false;
			if (operand_indices.length != 0) {
				if (!print(", Premises: [", out) || !print(operand_indices[0], out)) return false;
				for (unsigned int i = 1; i < operand_indices.length; i++) {
					if (!print("], [", out) || !print(operand_indices[i], out)) return false;
				}
				if (!print("]", out)) return false;
			}
			if (!print('\n', out))
				return false;
		}

		operand_states.clear();
		operand_indices.clear();
	}
	free_proof_states(proof_states);
	return true;
}

template<typename Formula, bool IsIntuitionistic, typename Canonicalizer = identity_canonicalizer>
struct natural_deduction
{
	typedef Formula Language;
	typedef nd_step<Formula> Proof;
	typedef nd_step_type ProofType;
	typedef typename Formula::Term Term;
	typedef Canonicalizer ProofCanonicalizer;

	static constexpr bool Intuitionistic = IsIntuitionistic;

	static inline Proof* new_axiom(Formula* axiom) {
		return new_parameterized_step<nd_step_type::AXIOM>(axiom);
	}

	static inline Proof* new_beta(Formula* first, Formula* second) {
		return new_binary_step<nd_step_type::BETA_EQUIVALENCE>(new_formula_parameter(first), new_formula_parameter(second));
	}

	static inline Proof* new_comparison_introduction(Formula* comparison) {
		return new_parameterized_step<nd_step_type::COMPARISON_INTRODUCTION>(comparison);
	}

	static inline Proof* new_inequality_introduction(Formula* inequality) {
		return new_parameterized_step<nd_step_type::INEQUALITY_INTRODUCTION>(inequality);
	}

	template<typename... Proofs>
	static inline Proof* new_conjunction_intro(Proof* proof, Proofs&&... other_proofs) {
		return new_array_step<nd_step_type::CONJUNCTION_INTRODUCTION, 2>(proof, std::forward<Proofs>(other_proofs)...);
	}

	template<template<typename> class Array>
	static inline Proof* new_conjunction_intro(const Array<Proof*>& operands) {
		return new_array_step<nd_step_type::CONJUNCTION_INTRODUCTION, 2>(operands);
	}

	template<template<typename> class Array>
	static inline Proof* new_conjunction_elim(Proof* proof, const Array<unsigned int>& indices) {
		return new_binary_step<nd_step_type::CONJUNCTION_ELIMINATION>(proof, new_array_parameter(indices));
	}

	static inline Proof* new_conjunction_elim_left(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_LEFT>(proof);
	}

	static inline Proof* new_conjunction_elim_right(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_RIGHT>(proof);
	}

	static inline Proof* new_disjunction_intro(Proof* proof, Formula* parameter, unsigned int index_to_insert) {
		return new_ternary_step<nd_step_type::DISJUNCTION_INTRODUCTION>(proof, new_formula_parameter(parameter), new_uint_parameter(index_to_insert));
	}

	static inline Proof* new_disjunction_intro_left(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_LEFT>(proof, new_formula_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_right(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT>(proof, new_formula_parameter(parameter));
	}

	template<template<typename> class Array>
	static inline Proof* new_disjunction_elim(Proof* disjunction, Array<Proof*> operands) {
		return new_array_step<nd_step_type::DISJUNCTION_ELIMINATION, 3>(make_prepended_array_view(disjunction, operands));
	}

	static inline Proof* new_implication_intro(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_step_type::AXIOM) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(proof, assumption);
			return NULL;
		}
		return new_binary_step<nd_step_type::IMPLICATION_INTRODUCTION>(proof, assumption);
	}

	static inline Proof* new_implication_elim(Proof* implication, Proof* antecedent) {
		return new_binary_step<nd_step_type::IMPLICATION_ELIMINATION>(implication, antecedent);
	}

	static inline Proof* new_biconditional_intro(Proof* forward, Proof* backward) {
		return new_binary_step<nd_step_type::BICONDITIONAL_INTRODUCTION>(forward, backward);
	}

	static inline Proof* new_biconditional_elim_left(Proof* biconditional, Proof* left) {
		return new_binary_step<nd_step_type::BICONDITIONAL_ELIMINATION_LEFT>(biconditional, left);
	}

	static inline Proof* new_biconditional_elim_right(Proof* biconditional, Proof* right) {
		return new_binary_step<nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT>(biconditional, right);
	}

	template<template<typename> class Array>
	static inline Proof* new_equality_elim(Proof* equality, Proof* formula, const Array<unsigned int>& term_indices) {
		return new_ternary_step<nd_step_type::EQUALITY_ELIMINATION>(equality, formula, new_array_parameter(term_indices));
	}

	static inline Proof* new_proof_by_contradiction(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_step_type::AXIOM) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(proof, assumption);
			return NULL;
		}
		return new_binary_step<nd_step_type::PROOF_BY_CONTRADICTION>(proof, assumption);
	}

	static inline Proof* new_negation_elim(Proof* proof, Proof* negated) {
		return new_binary_step<nd_step_type::NEGATION_ELIMINATION>(proof, negated);
	}

	static inline Proof* new_falsity_elim(Proof* proof, Formula* formula) {
		return new_binary_step<nd_step_type::FALSITY_ELIMINATION>(proof, new_formula_parameter(formula));
	}

	static inline Proof* new_universal_intro(Proof* proof, unsigned int parameter) {
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_uint_parameter(parameter));
	}

	static inline Proof* new_universal_elim(Proof* proof, Term* term) {
		return new_binary_step<nd_step_type::UNIVERSAL_ELIMINATION>(proof, new_term_parameter(term));
	}

	static inline Proof* new_existential_intro(Proof* proof, unsigned int* term_indices, unsigned int term_index_count, Term* term) {
		return new_ternary_step<nd_step_type::EXISTENTIAL_INTRODUCTION>(proof, new_array_parameter(make_array_view(term_indices, term_index_count)), new_term_parameter(term));
	}

	static inline Proof* new_existential_elim(Proof* existential, Proof* proof) {
		return new_binary_step<nd_step_type::EXISTENTIAL_ELIMINATION>(existential, proof);
	}

private:
	static inline Proof* new_formula_parameter(Formula* parameter) {
		return new_parameterized_step<nd_step_type::FORMULA_PARAMETER>(parameter);
	}

	static inline Proof* new_uint_parameter(unsigned int parameter) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::PARAMETER)) return NULL;
		step->reference_count = 0;
		step->parameter = parameter;
		return step;
	}

	static inline Proof* new_term_parameter(Term* term) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::TERM_PARAMETER)) return NULL;
		step->reference_count = 0;
		step->term = term;
		term->reference_count++;
		return step;
	}

	template<template<typename> class Array>
	static inline Proof* new_array_parameter(const Array<unsigned int>& parameters)
	{
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::ARRAY_PARAMETER)) {
			return NULL;
		} else if (!array_init(step->parameters, max((decltype(parameters.length)) 1, parameters.length))) {
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < parameters.length; i++)
			step->parameters[i] = parameters[i];
		step->parameters.length = parameters.length;
		step->reference_count = 0;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_parameterized_step(Formula* parameter) {
		if (parameter == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		step->formula = parameter;
		parameter->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_unary_step(Proof* proof) {
		if (proof == NULL) {
			free_proof_operands(proof);
			return NULL;
		}

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) {
			free_proof_operands(proof);
			return NULL;
		}
		step->reference_count = 0;
		step->operands[0] = proof;
		step->operands[1] = NULL;
		step->operands[2] = NULL;
		proof->reference_count++;
		if (!proof->children.add(step)) {
			free(*step); free(step);
			free_proof_operands(proof);
			return NULL;
		}
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_binary_step(Proof* left, Proof* right) {
		if (left == NULL || right == NULL) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(left, right);
			return NULL;
		}

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(left, right);
			return NULL;
		}
		step->reference_count = 0;
		step->operands[0] = left;
		step->operands[1] = right;
		step->operands[2] = NULL;
		left->reference_count++;
		right->reference_count++;
		if (!left->children.add(step) || !right->children.add(step)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(left, right);
			free(*step); free(step); return NULL;
		}
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_ternary_step(Proof* first, Proof* second, Proof* third) {
		if (first == NULL || second == NULL || third == NULL) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(first, second, third);
			return NULL;
		}

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(first, second, third);
			return NULL;
		}
		step->reference_count = 0;
		step->operands[0] = first;
		step->operands[1] = second;
		step->operands[2] = third;
		first->reference_count++;
		second->reference_count++;
		third->reference_count++;
		if (!first->children.add(step) || !second->children.add(step) || !third->children.add(step)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(first, second, third);
			free(*step); free(step); return NULL;
		}
		return step;
	}

	template<unsigned int Index>
	static constexpr bool new_array_step_helper(nd_step<Formula>* step) {
		return true;
	}

	template<unsigned int Index, typename... Args>
	static inline bool new_array_step_helper(nd_step<Formula>* step, Proof* arg, Args&&... args) {
		step->operand_array[Index] = arg;
		arg->reference_count++;
		if (!arg->children.add(step)) {
			free(*arg); return false;
		} else if (!new_array_step_helper<Index + 1>(step, std::forward<Args>(args)...)) {
			arg->remove_child(step);
			free(*arg); return false;
		}
		return true;
	}

	template<nd_step_type Type, unsigned int MinOperandCount, typename... Args,
			typename std::enable_if<sizeof...(Args) >= MinOperandCount>::type* = nullptr>
	static inline Proof* new_array_step(Args&&... args) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(std::forward<Args>(args)...);
			return NULL;
		}
		step->reference_count = 0;
		if (!array_init(step->operand_array, sizeof...(Args))) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(std::forward<Args>(args)...);
			free(step); return NULL;
		}
		if (!new_array_step_helper<0>(step, std::forward<Args>(args)...)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(std::forward<Args>(args)...);
			free(step->operand_array); free(step); return NULL;
		}
		step->operand_array.length = sizeof...(Args);
		return step;
	}

	template<nd_step_type Type, unsigned int MinOperandCount, typename Array>
	static inline Proof* new_array_step(const Array& operands) {
		if (operands.size() < MinOperandCount) {
			fprintf(stderr, "natural_deduction.new_array_step ERROR: "
					"This proof step requires at least two operands.\n");
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(operands);
			return NULL;
		}

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(operands);
			return NULL;
		}
		step->reference_count = 0;
		if (!array_init(step->operand_array, operands.size())) {
			natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(operands);
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < operands.size(); i++) {
			step->operand_array[i] = operands[i];
			operands[i]->reference_count++;
			if (!operands[i]->children.add(step)) {
				for (unsigned int j = 0; j < i; j++) {
					operands[j]->remove_child(step);
					free(*operands[j]);
				}
				free(step->operand_array); free(step);
				natural_deduction<Formula, Intuitionistic, Canonicalizer>::free_proof_operands(operands);
				return NULL;
			}
		}
		step->operand_array.length = operands.size();
		return step;
	}

	static inline void free_proof_operands() { }

	template<typename... Args>
	static inline void free_proof_operands(Proof* operand, Args&&... operands) {
		if (operand != NULL && operand->reference_count == 0) { operand->reference_count = 1; free(*operand); free(operand); }
		free_proof_operands(std::forward<Args>(operands)...);
	}

	template<typename Array>
	static inline void free_proof_operands(const Array& operands) {
		for (unsigned int i = 0; i < operands.size(); i++)
			if (operands[i] != NULL && operands[i]->reference_count == 0) { operands[i]->reference_count = 1; free(*operands[i]); free(operands[i]); }
	}
};

template<typename Formula>
struct nd_canonicalizer {
	inline bool operator() (const nd_step<Formula>* first, const nd_step<Formula>* second) {
		return *first > *second;
	}
};

template<typename Formula, typename... ProofMap>
bool canonicalize(const nd_step<Formula>& proof,
		array<const nd_step<Formula>*>& canonical_order,
		ProofMap&&... proof_map)
{
	hash_map<const nd_step<Formula>*, unsigned int> in_degrees(128);
	if (!compute_in_degrees(&proof, in_degrees, std::forward<ProofMap>(proof_map)...)) return false;

	array<const nd_step<Formula>*> heap(32);
	for (const auto& entry : in_degrees) {
		if (entry.value == 0 && !heap.add(entry.key))
			return false;
	}

	std::make_heap(heap.begin(), heap.end(), nd_canonicalizer<Formula>());
	while (heap.length > 0) {
		std::pop_heap(heap.begin(), heap.end(), nd_canonicalizer<Formula>());
		const nd_step<Formula>* node = heap.pop();
		if (!canonical_order.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			const nd_step<Formula>* operand = map_const(
					operands[i], std::forward<ProofMap>(proof_map)...);
			unsigned int& degree = in_degrees.get(operand);
			degree--;

			if (degree == 0) {
				if (!heap.add(operand))
					return false;
				std::push_heap(heap.begin(), heap.end(), nd_canonicalizer<Formula>());
			}
		}
	}

	reverse(canonical_order);
	return true;
}

template<typename Formula,
	typename ConjunctionIntroductions,
	typename ConjunctionEliminations,
	typename UniversalIntroductions,
	typename UniversalEliminations,
	typename TermIndices,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula>* const* proof,
		unsigned int step_index,
		unsigned int& formula_counter,
		array<unsigned int>& available_parameters,
		ConjunctionIntroductions& conjunction_introductions,
		ConjunctionEliminations& conjunction_eliminations,
		UniversalIntroductions& universal_introductions,
		UniversalEliminations& universal_eliminations,
		TermIndices& term_indices, ProofMap&&... proof_map)
{
	typedef typename Formula::TermType TermType;

	const nd_step<Formula>& current_step = *proof[step_index];
	unsigned int index;
	const nd_step<Formula>* operand;
	const nd_step<Formula>* const* operands;
	switch (current_step.type) {
	case nd_step_type::PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
		/* these aren't actual proof steps in the calculus */
		return 0.0;
	case nd_step_type::AXIOM:
		formula_counter++;
		get_parameters(*current_step.formula, available_parameters);
		if (available_parameters.length > 1) {
			sort(available_parameters); unique(available_parameters);
		}
		/* the contribution from the axiom_prior is computed at the end */
		return 0.0;
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
		formula_counter++;
		return -LOG_ND_RULE_COUNT;
	case nd_step_type::BETA_EQUIVALENCE:
		/* TODO: this is a hack; correct it */
		formula_counter++;
		return -LOG_ND_RULE_COUNT;
	case nd_step_type::FALSITY_ELIMINATION:
		/* TODO: this is a hack; correct it */
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::CONJUNCTION_ELIMINATION:
		operand = map_const(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		if (!conjunction_eliminations.add(make_array_view(operand->parameters.data, operand->parameters.length)))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::CONJUNCTION_INTRODUCTION:
		operands = current_step.operand_array.data;
		if (!conjunction_introductions.add(make_pair(
				make_array_view(operands, current_step.operand_array.length), make_array_view(proof, step_index))))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT;
	case nd_step_type::IMPLICATION_INTRODUCTION: /* TODO: is this correct? */
	case nd_step_type::IMPLICATION_ELIMINATION:
	case nd_step_type::BICONDITIONAL_INTRODUCTION:
	case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT:
	case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT:
	case nd_step_type::PROOF_BY_CONTRADICTION: /* TODO: is this correct? */
	case nd_step_type::NEGATION_ELIMINATION:
	case nd_step_type::EXISTENTIAL_ELIMINATION:
		return -LOG_ND_RULE_COUNT - 2*log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::DISJUNCTION_ELIMINATION:
		return -LOG_ND_RULE_COUNT - (current_step.operand_array.length - 1) * log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
		/* TODO: we need to compute the prior on the new formula */
		/* TODO: this is a hack; correct it */
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::UNIVERSAL_INTRODUCTION:
		operand = map_const(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		/* we need to compute the prior on the parameter */
		if (!universal_introductions.add(make_pair(operand->parameter, available_parameters)))
			exit(EXIT_FAILURE);
		
		index = available_parameters.index_of(operand->parameter);
		if (index < available_parameters.length)
			shift_left(available_parameters.data + index, available_parameters.length - index - 1);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::UNIVERSAL_ELIMINATION:
		operand = map_const(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand->term->type == TermType::PARAMETER) {
			if (available_parameters.ensure_capacity(available_parameters.length + 1)) exit(EXIT_FAILURE);
			add_sorted<true>(available_parameters, operand->term->parameter);
		}
		/* we need to compute the prior on the term */
		if (!universal_eliminations.add(*operand->term)) exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		/* we need to compute the prior on the parameter (it's a list of term indices) */
		operand = map_const(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		if (!term_indices.add(make_array_view(operand->parameters.data, operand->parameters.length)))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::EQUALITY_ELIMINATION:
		/* we need to compute the prior on the parameter (it's a list of term indices) */
		operand = map_const(current_step.operands[2], std::forward<ProofMap>(proof_map)...);
		if (!term_indices.add(make_array_view(operand->parameters.data, operand->parameters.length)))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - 2*log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "log_probability ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula>
struct changes_collector {
	typedef typename Formula::Term Term;

	array_multiset<Formula*, false>& axioms;
	array_multiset<Term, true>& universal_eliminations;
};

template<typename Formula>
inline bool visit_node(const nd_step<Formula>* proof, changes_collector<Formula>& collector) {
	if (proof->type == nd_step_type::AXIOM)
		return collector.axioms.add(proof->formula);
	else if (proof->type == nd_step_type::UNIVERSAL_ELIMINATION)
		return collector.universal_eliminations.add(*proof->operands[1]->term);
	return true;
}

template<typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename ConjunctionEliminationPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior>
struct canonicalized_proof_prior
{
	typedef typename std::remove_pointer<typename AxiomPrior::ObservationType>::type Formula;
	typedef typename Formula::Term Term;

	struct prior_state_changes {
		array_multiset<Formula*, false> proof_axioms;
		array_multiset<Term, true> universal_eliminations;
		typename AxiomPrior::PriorStateChanges axiom_prior_state;

		prior_state_changes() : proof_axioms(16), universal_eliminations(8) { }
	};

	struct prior_state {
		hash_multiset<Formula*, false> proof_axioms;
		hash_multiset<Term, true> universal_eliminations;
		typename AxiomPrior::PriorState axiom_prior_state;

		prior_state() : proof_axioms(64), universal_eliminations(16) { }

		prior_state(const prior_state& src, const hash_map<const Formula*, Formula*>& formula_map) :
			proof_axioms(src.proof_axioms.counts.table.capacity),
			universal_eliminations(src.universal_eliminations.counts.table.capacity),
			axiom_prior_state(src.axiom_prior_state, formula_map)
		{
			if (!init_helper(src, formula_map))
				exit(EXIT_FAILURE);
		}

		inline bool add(const prior_state_changes& changes)
		{
			if (!proof_axioms.add(changes.proof_axioms)) {
				return false;
			} else if (!universal_eliminations.add(changes.universal_eliminations)) {
				proof_axioms.template subtract<true>(changes.proof_axioms);
				return false;
			} else if (!axiom_prior_state.add(changes.axiom_prior_state)) {
				proof_axioms.template subtract<true>(changes.proof_axioms);
				universal_eliminations.template subtract<true>(changes.universal_eliminations);
				return false;
			}
			return true;
		}

		template<bool HasPrior = true>
		inline bool add(const nd_step<Formula>* proof, const array<Formula*>& extra_axioms,
				const canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior)
		{
			array_multiset<Formula*, false> proof_axiom_changes(16);
			array_multiset<Term, true> universal_elimination_changes(8);
			changes_collector<Formula> collector = {proof_axiom_changes, universal_elimination_changes};
			if (!visit<true>(proof, collector)
			 || !proof_axioms.counts.check_size(proof_axioms.counts.table.size + proof_axiom_changes.counts.size)
			 || (HasPrior && !universal_eliminations.add(universal_elimination_changes)))
				return false;
			array<Formula*> new_axioms(max((size_t) 1, proof_axiom_changes.counts.size + extra_axioms.length));
			for (const auto& pair : proof_axiom_changes.counts) {
				bool contains; unsigned int bucket;
				unsigned int& count = proof_axioms.counts.get(pair.key, contains, bucket);
				if (!contains) {
					proof_axioms.counts.table.keys[bucket] = pair.key;
					proof_axioms.counts.values[bucket] = pair.value;
					proof_axioms.counts.table.size++;
					new_axioms.add(pair.key);
				} else {
					count += pair.value;
				}
			}
			proof_axioms.sum += proof_axiom_changes.sum;

			for (Formula* extra_axiom : extra_axioms) {
				bool contains = false;
				for (Formula* existing_axiom : new_axioms) {
					if (existing_axiom == extra_axiom || *existing_axiom == *extra_axiom) {
						contains = true;
						break;
					}
				}
				if (!contains)
					new_axioms.add(extra_axiom);
			}

			for (unsigned int i = 0; i < new_axioms.length; i++) {
				if (!axiom_prior_state.add(new_axioms[i], prior.axiom_prior)) {
					for (unsigned int j = 0; j < i; j++)
						axiom_prior_state.subtract(new_axioms[j], prior.axiom_prior);
					if (HasPrior) universal_eliminations.template subtract<true>(universal_elimination_changes);
					proof_axioms.template subtract<true>(proof_axiom_changes);
					return false;
				}
			}
			return true;
		}

		inline void subtract(const prior_state_changes& changes) {
			proof_axioms.template subtract<true>(changes.proof_axioms);
			universal_eliminations.template subtract<true>(changes.universal_eliminations);
			axiom_prior_state.subtract(changes.axiom_prior_state);
		}

		template<bool HasPrior = true>
		inline void subtract(const nd_step<Formula>* proof, const array<Formula*>& extra_axioms,
				const canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior)
		{
			array_multiset<Formula*, false> proof_axiom_changes(16);
			array_multiset<Term, true> universal_elimination_changes(8);
			changes_collector<Formula> collector = {proof_axiom_changes, universal_elimination_changes};
			if (!visit<true>(proof, collector)) {
				fprintf(stderr, "canonicalized_proof_prior.prior_state.subtract ERROR: Unable to collect proof changes.\n");
				return;
			}
			array<Formula*> old_axioms(max((size_t) 1, proof_axiom_changes.counts.size + extra_axioms.length));
			for (const auto& pair : proof_axiom_changes.counts) {
				bool contains; unsigned int bucket;
				unsigned int& count = proof_axioms.counts.get(pair.key, contains, bucket);
#if !defined(NDEBUG)
				if (!contains || count < pair.value)
					fprintf(stderr, "canonicalized_proof_prior.prior_state.subtract WARNING:"
							" Attempted to remove more items from a bin than it contains.\n");
#endif
				count -= pair.value;
				if (count == 0) {
					old_axioms.add(pair.key);
					proof_axioms.counts.remove_at(bucket);
				}
			}
			proof_axioms.sum -= proof_axiom_changes.sum;

			for (Formula* extra_axiom : extra_axioms) {
				bool contains = false;
				for (Formula* existing_axiom : old_axioms) {
					if (existing_axiom == extra_axiom || *existing_axiom == *extra_axiom) {
						contains = true;
						break;
					}
				}
				if (!contains)
					old_axioms.add(extra_axiom);
			}

			for (unsigned int i = 0; i < old_axioms.length; i++)
				axiom_prior_state.subtract(old_axioms[i], prior.axiom_prior);
			if (HasPrior)
				universal_eliminations.template subtract<true>(universal_elimination_changes);
		}

		template<typename Theory>
		bool check_proof_axioms(const Theory& T) const
		{
			typedef typename Formula::TermType TermType;

			bool success = true;
			array<const nd_step<Formula>*> axiom_collection(64);
			array_multiset<Formula*, false> computed_axioms(64);
			for (const nd_step<Formula>* proof : T.observations) {
				array_multiset<const nd_step<Formula>*, false> axioms(16);
				if (!get_proof_steps<nd_step_type::AXIOM>(proof, axioms)) {
					fprintf(stderr, "canonicalized_proof_prior.prior_state.check_proof_axioms ERROR: `get_proof_steps` failed.\n");
					return false;
				}

				for (const auto& entry : axioms.counts) {
					bool contains;
					unsigned int count = proof_axioms.counts.get(entry.key->formula, contains);
					if (!contains || count == 0) {
						print("canonicalized_proof_prior.prior_state.check_proof_axioms WARNING: "
								"Found axiom of a proof in `observations` that is not in `proof_axioms`.\n", stderr);
						print("  Axiom: ", stderr); print(*entry.key->formula, stderr); print('\n', stderr);
						success = false;
					}

					if (!computed_axioms.add(entry.key->formula, entry.value)
					 || !axiom_collection.add(entry.key)) return false;
				}
			}

			for (const auto& entry : proof_axioms.counts) {
				bool contains;
				unsigned int count = computed_axioms.counts.get(entry.key, contains);
				if (!contains) {
					print("canonicalized_proof_prior.prior_state.check_proof_axioms WARNING: "
							"`proof_axioms` contains an axiom that does not belong to a proof"
							" in `observations`.\n", stderr);
					print("  Axiom: ", stderr); print(*entry.key, stderr); print('\n', stderr);
					success = false;
				} else if (entry.value != count) {
					print("canonicalized_proof_prior.prior_state.check_proof_axioms WARNING: The axiom '", stderr);
					print(*entry.key, stderr); print("' has expected frequency ", stderr);
					print(count, stderr); print(" but has computed frequency ", stderr);
					print(entry.value, stderr); print('\n', stderr);
					success = false;
				}
			}

			/* check that definition axioms are correctly stored in the appropriate `definitions` array */
			for (const nd_step<Formula>* axiom : axiom_collection) {
				if (axiom->formula->type == TermType::EQUALS && axiom->formula->binary.left->type == TermType::CONSTANT) {
					unsigned int constant = axiom->formula->binary.left->constant;
					if (constant < T.new_constant_offset) continue;
					if (!T.ground_concepts[constant - T.new_constant_offset].definitions.contains(axiom)) {
						print("canonicalized_proof_prior.prior_state.check_proof_axioms WARNING: The axiom '", stderr);
						print(*axiom->formula, stderr);
						print("' does not exist in the appropriate `definitions` array.\n", stderr);
						success = false;
					}
				}
			}
			return success;
		}

		template<typename Theory, typename TheorySampleCollector>
		bool check_universal_eliminations(const Theory& T, TheorySampleCollector& theory_sample_collector) const
		{
			array_multiset<Formula*, false> proof_axiom_changes(16);
			array_multiset<Term, true> universal_elimination_changes(8);
			changes_collector<Formula> collector = {proof_axiom_changes, universal_elimination_changes};
			for (const nd_step<Formula>* proof : T.observations) {
				if (!theory_sample_collector.has_prior(proof)) continue;
				if (!visit<true>(proof, collector))
					return false;
			}

			bool success = true;
			for (const auto& entry : universal_elimination_changes.counts) {
				bool contains;
				unsigned int count = universal_eliminations.counts.get(entry.key, contains);
				if (!contains) {
					print("canonicalized_proof_prior.prior_state.check_universal_eliminations WARNING:"
							" The observation's proofs contain the term '", stderr);
					print(entry.key, stderr); print("' but `universal_eliminations` does not.\n", stderr);
					success = false;
				} else if (count != entry.value) {
					fprintf(stderr, "canonicalized_proof_prior.prior_state.check_universal_eliminations WARNING:"
							" The observation's proofs contain %u instances of the term '", entry.value);
					print(entry.key, stderr); fprintf(stderr, "' but `universal_eliminations` contains %u.\n", count);
					success = false;
				}
			} for (const auto& entry : universal_eliminations.counts) {
				if (!universal_eliminations.counts.table.contains(entry.key)) {
					print("canonicalized_proof_prior.prior_state.check_universal_eliminations WARNING:"
							" `universal_eliminations` contains the term '", stderr);
					print(entry.key, stderr); print("' but the observation's proofs do not.\n", stderr);
					success = false;
				}
			}
			return success;
		}

		static inline bool get_formula_map(const prior_state& state,
				hash_map<const Formula*, unsigned int>& formula_map)
		{
			for (const auto& entry : state.proof_axioms.counts) {
				if (!::get_formula_map(entry.key, formula_map))
					return false;
			} for (const auto& entry : state.universal_eliminations.counts) {
				if (!::get_formula_map(entry.key, formula_map))
					return false;
			}
			return AxiomPrior::PriorState::get_formula_map(state.axiom_prior_state, formula_map);
		}

		template<typename Stream>
		static bool read(prior_state& state,
				Stream& in, Formula** formulas)
		{
			typedef typename Formula::Term Term;

			decltype(state.proof_axioms.counts.table.size) proof_axiom_count;
			if (!core::read(proof_axiom_count, in)
			 || !init(state.proof_axioms, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (proof_axiom_count == 0 ? 1 : proof_axiom_count)) + 1)))
			{
				return false;
			}
			unsigned int index;
			for (unsigned int i = 0; i < proof_axiom_count; i++) {
				if (!core::read(index, in)) {
					core::free(state.proof_axioms);
					return false;
				}
				Formula* formula = formulas[index];
				unsigned int bucket = state.proof_axioms.counts.table.index_to_insert(formula);
				if (!core::read(state.proof_axioms.counts.values[bucket], in)) {
					core::free(state.proof_axioms);
					return false;
				}
				state.proof_axioms.counts.table.keys[bucket] = formula;
				state.proof_axioms.counts.table.size++;
				state.proof_axioms.sum += state.proof_axioms.counts.values[bucket];
			}

			decltype(state.universal_eliminations.counts.table.size) universal_elimination_count;
			if (!core::read(universal_elimination_count, in)
			 || !init(state.universal_eliminations, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (universal_elimination_count == 0 ? 1 : universal_elimination_count)) + 1)))
			{
				core::free(state.proof_axioms);
				return false;
			}
			Term& term = *((Term*) alloca(sizeof(Term)));
			for (unsigned int i = 0; i < universal_elimination_count; i++) {
				if (!::read(term, in, formulas)) {
					core::free(state.proof_axioms);
					core::free(state.universal_eliminations);
					return false;
				}
				term.reference_count = 1;
				unsigned int bucket = state.universal_eliminations.counts.table.index_to_insert(term);
				if (!core::read(state.universal_eliminations.counts.values[bucket], in)) {
					core::free(state.proof_axioms);
					core::free(state.universal_eliminations);
					return false;
				}
				core::move(term, state.universal_eliminations.counts.table.keys[bucket]);
				state.universal_eliminations.counts.table.size++;
				state.universal_eliminations.sum += state.universal_eliminations.counts.values[bucket];
			}

			if (!AxiomPrior::PriorState::read(state.axiom_prior_state, in, formulas)) {
				core::free(state.proof_axioms);
				core::free(state.universal_eliminations);
				return false;
			}
			return true;
		}

		template<typename Stream>
		static bool write(const prior_state& state, Stream& out,
				const hash_map<const Formula*, unsigned int>& formula_map)
		{
			if (!core::write(state.proof_axioms.counts.table.size, out))
				return false;
			for (const auto& entry : state.proof_axioms.counts) {
				if (!core::write(formula_map.get(entry.key), out)
				 || !core::write(entry.value, out))
					return false;
			}

			if (!core::write(state.universal_eliminations.counts.table.size, out))
				return false;
			for (const auto& entry : state.universal_eliminations.counts) {
				if (!::write(entry.key, out, formula_map)
				 || !core::write(entry.value, out))
					return false;
			}

			return AxiomPrior::PriorState::write(state.axiom_prior_state, out, formula_map);
		}

		static bool clone(const prior_state& src, prior_state& dst, const hash_map<const Formula*, Formula*>& formula_map)
		{
			if (!init(dst.proof_axioms, src.proof_axioms.counts.table.capacity)) {
				return false;
			} else if (!init(dst.universal_eliminations, src.universal_eliminations.counts.table.capacity)) {
				core::free(dst.proof_axioms);
				return false;
			} else if (!AxiomPrior::PriorState::clone(src.axiom_prior_state, dst.axiom_prior_state, formula_map)) {
				core::free(dst.proof_axioms);
				core::free(dst.universal_eliminations);
				return false;
			} else if (!dst.init_helper(src, formula_map)) {
				core::free(dst.proof_axioms);
				core::free(dst.universal_eliminations);
				core::free(dst.axiom_prior_state);
				return false;
			}
			return true;
		}

		static inline void free(prior_state& state) {
			core::free(state.proof_axioms);
			core::free(state.universal_eliminations);
			core::free(state.axiom_prior_state);
		}

		inline bool init_helper(const prior_state& src, const hash_map<const Formula*, Formula*>& formula_map)
		{
			for (const auto& entry : src.proof_axioms.counts) {
#if !defined(NDEBUG)
				bool contains;
				Formula* formula = formula_map.get(entry.key, contains);
				if (!contains)
					fprintf(stderr, "canonicalized_proof_prior.prior_state.init_helper WARNING: Formula doesn't exist in `formula_map`.\n");
				proof_axioms.counts.put(formula, entry.value);
#else
				proof_axioms.counts.put(formula_map.get(entry.key), entry.value);
#endif
			}
			proof_axioms.sum = src.proof_axioms.sum;

			for (const auto& entry : src.universal_eliminations.counts) {
				unsigned int index = universal_eliminations.counts.table.index_to_insert(entry.key);
				if (!::clone(entry.key, universal_eliminations.counts.table.keys[index], formula_map))
					return false;
				universal_eliminations.counts.values[index] = entry.value;
				universal_eliminations.counts.table.size++;
			}
			universal_eliminations.sum = src.universal_eliminations.sum;
			return true;
		}
	};

	typedef prior_state PriorState;
	typedef prior_state_changes PriorStateChanges;

	AxiomPrior axiom_prior;
	ConjunctionIntroductionPrior conjunction_introduction_prior;
	ConjunctionEliminationPrior conjunction_elimination_prior;
	UniversalIntroductionPrior universal_introduction_prior;
	UniversalEliminationPrior universal_elimination_prior;
	TermIndicesPrior term_indices_prior;
	ProofLengthPrior proof_length_prior;

	double log_negated_antecedent_probability;
	double log_consequent_probability;

	canonicalized_proof_prior(
			const AxiomPrior& axiom_prior,
			const ConjunctionIntroductionPrior& conjunction_introduction_prior,
			const ConjunctionEliminationPrior& conjunction_elimination_prior,
			const UniversalIntroductionPrior& universal_introduction_prior,
			const UniversalEliminationPrior& universal_elimination_prior,
			const TermIndicesPrior& term_indices_prior,
			const ProofLengthPrior& proof_length_prior,
			double consequent_probability) :
		axiom_prior(axiom_prior),
		conjunction_introduction_prior(conjunction_introduction_prior),
		conjunction_elimination_prior(conjunction_elimination_prior),
		universal_introduction_prior(universal_introduction_prior),
		universal_elimination_prior(universal_elimination_prior),
		term_indices_prior(term_indices_prior),
		proof_length_prior(proof_length_prior),
		log_negated_antecedent_probability(log(1.0 - consequent_probability)),
		log_consequent_probability(log(consequent_probability))
	{ }

	canonicalized_proof_prior(const canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& src) :
		axiom_prior(src.axiom_prior),
		conjunction_introduction_prior(src.conjunction_introduction_prior),
		conjunction_elimination_prior(src.conjunction_elimination_prior),
		universal_introduction_prior(src.universal_introduction_prior),
		universal_elimination_prior(src.universal_elimination_prior),
		term_indices_prior(src.term_indices_prior),
		proof_length_prior(src.proof_length_prior),
		log_negated_antecedent_probability(src.log_negated_antecedent_probability),
		log_consequent_probability(src.log_consequent_probability)
	{ }

	template<typename Formula>
	inline bool add_axiom(Formula* axiom) {
		return axiom_prior.add(axiom);
	}

	template<typename Formula>
	inline void remove_axiom(Formula* axiom) {
		axiom_prior.remove(axiom);
	}
};

template<typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename ConjunctionEliminationPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior>
inline canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>
make_canonicalized_proof_prior(
		const AxiomPrior& axiom_prior,
		const ConjunctionIntroductionPrior& conjunction_introduction_prior,
		const ConjunctionEliminationPrior& conjunction_elimination_prior,
		const UniversalIntroductionPrior& universal_introduction_prior,
		const UniversalEliminationPrior& universal_elimination_prior,
		const TermIndicesPrior& term_indices_prior,
		const ProofLengthPrior& proof_length_prior,
		double consequent_probability)
{
	return canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>(
			axiom_prior, conjunction_introduction_prior, conjunction_elimination_prior, universal_introduction_prior, universal_elimination_prior, term_indices_prior, proof_length_prior, consequent_probability);
}

template<typename Formula, typename... ProofMap>
inline bool count_axioms(nd_step<Formula>& proof,
		array_multiset<Formula*, false>& axioms,
		ProofMap&&... proof_map)
{
	/* count the axioms */
	array<nd_step<Formula>*> stack(8);
	stack[stack.length++] = map(&proof, std::forward<ProofMap>(proof_map)...);
	while (stack.length > 0) {
		nd_step<Formula>* node = stack.pop();
		if (node->type == nd_step_type::AXIOM
		 && !axioms.add_unsorted(node->formula))
			return false;

		unsigned int operand_count;
		nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			nd_step<Formula>* operand = map(
					operands[i], std::forward<ProofMap>(proof_map)...);

			if (!stack.add(operand))
				return false;
		}
	}
	return true;
}

template<typename Formula,
	typename ConjunctionIntroductions,
	typename ConjunctionEliminations,
	typename UniversalIntroductions,
	typename UniversalEliminations,
	typename TermIndices,
	typename ProofLengthPrior,
	typename... ProofMap>
inline double log_probability_helper(
		nd_step<Formula>& proof,
		array_multiset<Formula*, false>& axioms,
		ConjunctionIntroductions& conjunction_introductions,
		ConjunctionEliminations& conjunction_eliminations,
		UniversalIntroductions& universal_introductions,
		UniversalEliminations& universal_eliminations,
		TermIndices& term_indices,
		ProofLengthPrior& proof_length_prior,
		double log_negated_antecedent_probability,
		double log_consequent_probability,
		ProofMap&&... proof_map)
{
	if (!count_axioms(proof, axioms, std::forward<ProofMap>(proof_map)...))
		exit(EXIT_FAILURE);

	array<const nd_step<Formula>*> canonical_order(64);
	if (!canonicalize(proof, canonical_order, std::forward<ProofMap>(proof_map)...)) {
		fprintf(stderr, "log_probability ERROR: Unable to canonicalize proof.\n");
		exit(EXIT_FAILURE);
	}

	double value = log_probability(canonical_order.length, proof_length_prior);
	unsigned int formula_counter = 0;
	array<unsigned int> available_parameters(16);
	log_cache<double>::instance().ensure_size(canonical_order.length + 1);
	for (unsigned int i = 0; i < canonical_order.length; i++) {
		if (canonical_order[i]->type == nd_step_type::IMPLICATION_INTRODUCTION) {
			const nd_step<Formula>* operand = map_const(
					canonical_order[i]->operands[0], std::forward<ProofMap>(proof_map)...);
			if (operand->type == nd_step_type::FALSITY_ELIMINATION) {
				value += log_negated_antecedent_probability;
			} else {
				value += log_consequent_probability;
			}
		}
		value += log_probability(canonical_order.data, i, formula_counter, available_parameters,
				conjunction_introductions, conjunction_eliminations, universal_introductions,
				universal_eliminations, term_indices, std::forward<ProofMap>(proof_map)...);
	}
	return value;
}

template<
	typename Formula, typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename ConjunctionEliminationPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula>& proof,
		canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior,
		ProofMap&&... proof_map)
{
	typedef typename ConjunctionIntroductionPrior::ObservationCollection ConjunctionIntroductions;
	typedef typename ConjunctionEliminationPrior::ObservationCollection ConjunctionEliminations;
	typedef typename UniversalIntroductionPrior::ObservationCollection UniversalIntroductions;
	typedef typename UniversalEliminationPrior::ObservationCollection UniversalEliminations;
	typedef typename TermIndicesPrior::ObservationCollection TermIndices;

	array_multiset<Formula*, false> axioms(16);
	ConjunctionIntroductions conjunction_introductions;
	ConjunctionEliminations conjunction_eliminations;
	UniversalIntroductions universal_introductions;
	UniversalEliminations universal_eliminations;
	TermIndices term_indices;
	double value = log_probability_helper(proof, axioms, conjunction_introductions, conjunction_eliminations,
			universal_introductions, universal_eliminations, prior.log_continue_probability,
			prior.log_stop_probability, term_indices, prior.proof_length_prior, prior.log_negated_antecedent_probability,
			prior.log_consequent_probability, std::forward<ProofMap>(proof_map)...);
	return value + log_probability(axioms, prior.axiom_prior)
		 + log_probability(conjunction_introductions, prior.conjunction_introduction_prior)
		 + log_probability(conjunction_eliminations, prior.conjunction_elimination_prior)
		 + log_probability(universal_introductions, prior.universal_introduction_prior)
		 + log_probability(universal_eliminations, prior.universal_elimination_prior)
		 + log_probability(term_indices, prior.term_indices_prior);
}

template<
	typename Formula, typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename ConjunctionEliminationPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior,
	typename TheorySampleCollector>
double log_probability(
		const array<nd_step<Formula>*>& proofs,
		const array<Formula*>& extra_observations,
		canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior,
		TheorySampleCollector& theory_sample_collector)
{
	typedef typename ConjunctionIntroductionPrior::ObservationCollection ConjunctionIntroductions;
	typedef typename ConjunctionEliminationPrior::ObservationCollection ConjunctionEliminations;
	typedef typename UniversalIntroductionPrior::ObservationCollection UniversalIntroductions;
	typedef typename UniversalEliminationPrior::ObservationCollection UniversalEliminations;
	typedef typename TermIndicesPrior::ObservationCollection TermIndices;

#if defined(DEBUG_LOG_PROBABILITY)
	unsigned int proof_index = 0;
#endif

	double value = 0.0;
	array_multiset<Formula*, false> axioms(16);
	UniversalEliminations universal_eliminations;
	for (nd_step<Formula>* entry : proofs) {
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "Proof with index %u:\n", proof_index);
		proof_index++;
#endif
		double proof_log_probability = 0.0;
		if (theory_sample_collector.has_prior(entry)) {
			ConjunctionIntroductions conjunction_introductions;
			ConjunctionEliminations conjunction_eliminations;
			UniversalIntroductions universal_introductions;
			TermIndices term_indices;
			double value = log_probability_helper(
					*entry, axioms, conjunction_introductions, conjunction_eliminations,
					universal_introductions, universal_eliminations, term_indices, prior.proof_length_prior,
					prior.log_negated_antecedent_probability, prior.log_consequent_probability);
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  log_probability_helper returned %lf.\n", value);
#endif
			proof_log_probability += value;

			value = log_probability(conjunction_introductions, prior.conjunction_introduction_prior);
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  log probability of `conjunction_introductions`: %lf.\n", value);
#endif
			proof_log_probability += value;
			value = log_probability(conjunction_eliminations, prior.conjunction_elimination_prior);
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  log probability of `conjunction_eliminations`: %lf.\n", value);
#endif
			proof_log_probability += value;
			value = log_probability(universal_introductions, prior.universal_introduction_prior);
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  log probability of `universal_introductions`: %lf.\n", value);
#endif
			proof_log_probability += value;
			value = log_probability(term_indices, prior.term_indices_prior);
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  log probability of `term_indices`: %lf.\n", value);
#endif
			proof_log_probability += value;
		} else {
#if defined(DEBUG_LOG_PROBABILITY)
			fprintf(stderr, "  theory_sample_collector.has_prior returned true.\n");
#endif
			if (!count_axioms(*entry, axioms))
				exit(EXIT_FAILURE);
		}
		value += proof_log_probability;
	}

	double current_value = log_probability(universal_eliminations, prior.universal_elimination_prior);
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  log probability of `universal_eliminations`: %lf.\n", current_value);
#endif
	value += current_value;

#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "Total log likelihood of `canonicalized_proof_prior`: %lf.\n", value);
#endif

	if (axioms.counts.size > 1)
		sort(axioms.counts.keys, axioms.counts.values, axioms.counts.size);
	return value + log_probability(axioms, extra_observations, prior.axiom_prior);
}

template<
	typename Formula, typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename ConjunctionEliminationPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior,
	typename PriorState,
	typename PriorStateChanges,
	typename TheorySampleCollector>
double log_probability_ratio(
		const array_map<nd_step<Formula>*, proof_substitution<Formula>>& proofs,
		const array<Formula*>& old_extra_axioms, const array<Formula*>& new_extra_axioms,
		canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, ConjunctionEliminationPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior,
		const PriorState& prior_state, PriorStateChanges& old_axioms, PriorStateChanges& new_axioms,
		TheorySampleCollector& theory_sample_collector)
{
	typedef typename ConjunctionIntroductionPrior::ObservationCollection ConjunctionIntroductions;
	typedef typename ConjunctionEliminationPrior::ObservationCollection ConjunctionEliminations;
	typedef typename UniversalIntroductionPrior::ObservationCollection UniversalIntroductions;
	typedef typename TermIndicesPrior::ObservationCollection TermIndices;

	double value = 0.0;
	for (const auto& entry : proofs) {
		if (theory_sample_collector.has_prior(entry.key)) {
			ConjunctionIntroductions old_conjunction_introductions, new_conjunction_introductions;
			ConjunctionEliminations old_conjunction_eliminations, new_conjunction_eliminations;
			UniversalIntroductions old_universal_introductions, new_universal_introductions;
			TermIndices old_term_indices, new_term_indices;
			value -= log_probability_helper(
					*entry.key, old_axioms.proof_axioms, old_conjunction_introductions, old_conjunction_eliminations,
					old_universal_introductions, old_axioms.universal_eliminations, old_term_indices, prior.proof_length_prior,
					prior.log_negated_antecedent_probability, prior.log_consequent_probability);
			value += log_probability_helper(
					*entry.key, new_axioms.proof_axioms, new_conjunction_introductions, new_conjunction_eliminations,
					new_universal_introductions, new_axioms.universal_eliminations, new_term_indices, prior.proof_length_prior,
					prior.log_negated_antecedent_probability, prior.log_consequent_probability, entry.value);

			value += log_probability_ratio(old_conjunction_introductions, new_conjunction_introductions, prior.conjunction_introduction_prior);
			value += log_probability_ratio(old_conjunction_eliminations, new_conjunction_eliminations, prior.conjunction_elimination_prior);
			value += log_probability_ratio(old_universal_introductions, new_universal_introductions, prior.universal_introduction_prior);
			value += log_probability_ratio(old_term_indices, new_term_indices, prior.term_indices_prior);
		} else {
			if (!count_axioms(*entry.key, old_axioms.proof_axioms)
			 || !count_axioms(*entry.key, new_axioms.proof_axioms, entry.value))
				exit(EXIT_FAILURE);
		}
	}

	value += log_probability_ratio(prior_state.universal_eliminations, old_axioms.universal_eliminations, new_axioms.universal_eliminations, prior.universal_elimination_prior);

	if (old_axioms.proof_axioms.counts.size > 1)
		sort(old_axioms.proof_axioms.counts.keys, old_axioms.proof_axioms.counts.values, old_axioms.proof_axioms.counts.size);
	if (new_axioms.proof_axioms.counts.size > 1)
		sort(new_axioms.proof_axioms.counts.keys, new_axioms.proof_axioms.counts.values, new_axioms.proof_axioms.counts.size);
	return value + log_probability_ratio(prior_state.proof_axioms,
			old_axioms.proof_axioms, new_axioms.proof_axioms,
			old_extra_axioms, new_extra_axioms,
			prior.axiom_prior, prior_state.axiom_prior_state,
			old_axioms.axiom_prior_state, new_axioms.axiom_prior_state);
}

#endif /* NATURAL_DEDUCTION_H_ */
