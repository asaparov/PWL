#ifndef NATURAL_DEDUCTION_H_
#define NATURAL_DEDUCTION_H_

#include <core/array.h>
#include <math/log.h>
#include <math/multiset.h>

#include "array_view.h"

using namespace core;

enum class nd_step_type : uint_fast16_t
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

	UNIVERSAL_INTRODUCTION,
	UNIVERSAL_ELIMINATION,
	EXISTENTIAL_INTRODUCTION,
	EXISTENTIAL_ELIMINATION,

	COUNT
};

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
};

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

template<typename FormulaPtr>
struct proof_state_formulas {
	typedef typename std::remove_pointer<FormulaPtr>::type Formula;

	proof_state<Formula>** states;
	unsigned int length;

	proof_state_formulas(proof_state<Formula>** states, unsigned int length) : states(states), length(length) { }

	inline FormulaPtr operator[] (size_t index) const {
		return states[index]->formula;
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

	bool discharged = false;
	for (unsigned int i = 0; i < first.length; i++) {
		if (*first.elements[i] == *first.excluded) { discharged = true; continue; }
		dst[dst.length] = first.elements[i];
		dst.length++;
	}
	return discharged;
}

template<typename T, typename ExcludedType, bool RemoveDuplicates = true>
bool pass_hypotheses(array<T>& dst, const array<T>& first, const exclude<T, ExcludedType>& second)
{
	if (!dst.ensure_capacity(first.length + second.length))
		return false;

	bool discharged = false;
	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		if (*second.excluded == *second.elements[j]) { j++; discharged = true; continue; }

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
		if (*second.excluded == *second.elements[j]) { j++; discharged = true; continue; }
		set_union_helper<RemoveDuplicates>(dst.data, dst.length, second.elements[j]);
		j++;
	}

	return discharged;
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

template<typename Formula, typename Canonicalizer>
inline Formula* try_canonicalize(Formula* formula, Canonicalizer& canonicalizer) {
	Formula* out = canonicalize(*formula, canonicalizer);
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
inline const nd_step<Formula>* map(
		const nd_step<Formula>* const step)
{
	return step;
}

template<typename Formula, typename ProofStep>
inline const nd_step<Formula>* map(
		const ProofStep& step,
		const proof_substitution<Formula>& substitution)
{
	bool contains;
	auto value = substitution.map.get(step, contains);
	const nd_step<Formula>* new_step = value;
	if (contains) return new_step;
	else return step;
}

template<typename BuiltInPredicates, typename Formula, typename Canonicalizer, typename... ProofMap>
bool check_proof(proof_state<Formula>& out,
		const nd_step<Formula>& proof,
		proof_state<Formula>** operand_states,
		unsigned int operand_count,
		Canonicalizer& canonicalizer,
		ProofMap&&... proof_map)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	Formula* formula;
	const nd_step<Formula>* second_operand;
	Term* parameter; Term* variable;
	array<unsigned int> constants = array<unsigned int>(8);
	exclude<const Formula*, Formula*>* assumptions;
	unsigned int max_variable;
	switch (proof.type) {
	case nd_step_type::AXIOM:
		if (!is_canonical(*proof.formula, canonicalizer)) {
			fprintf(stderr, "check_proof ERROR: Axiom is not in canonical form.\n");
			return false;
		}
		out.formula = proof.formula;
		out.formula->reference_count++;
		return out.assumptions.add(out.formula);
	case nd_step_type::BETA_EQUIVALENCE:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || proof.operands[0]->type != nd_step_type::FORMULA_PARAMETER || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		formula = beta_reduce(proof.operands[0]->formula);
		if (formula == NULL) return false;
		out.formula = beta_reduce(second_operand->formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		} else if (*formula != *out.formula) {
			fprintf(stderr, "check_proof ERROR: The two formula are not beta-equivalent.\n");
			free(*formula); if (formula->reference_count == 0) free(formula);
			free(*out.formula); if (out.formula->reference_count == 0) free(out.formula);
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
		out.formula = try_canonicalize(out.formula, canonicalizer);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, make_proof_state_assumptions(operand_states, operand_count));
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if ((proof.type == nd_step_type::CONJUNCTION_ELIMINATION && (operand_count != 2 || second_operand->type != nd_step_type::ARRAY_PARAMETER))
		 || (proof.type != nd_step_type::CONJUNCTION_ELIMINATION && operand_count != 1))
			return false;
		formula = operand_states[0]->formula;
		if (formula == NULL) {
			return false;
		} else if (formula->type != FormulaType::AND || formula->array.length < 2) {
			fprintf(stderr, "check_proof ERROR: Expected a conjunction.\n");
			free(*formula);
			if (formula->reference_count == 0)
				free(formula);
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
		free(*formula);
		if (formula->reference_count == 0)
			free(formula);
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_LEFT || proof.type == nd_step_type::DISJUNCTION_INTRODUCTION)
			out.formula = Formula::new_or(operand_states[0]->formula, second_operand->formula);
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT)
			out.formula = Formula::new_or(second_operand->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		second_operand->formula->reference_count++;
		out.formula = try_canonicalize(out.formula, canonicalizer);
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
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || second_operand->type != nd_step_type::AXIOM)
			return false;
		out.formula = Formula::new_if_then(operand_states[1]->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		operand_states[1]->formula->reference_count++;
		out.formula = try_canonicalize(out.formula, canonicalizer);
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
		out.formula = try_canonicalize(out.formula, canonicalizer);
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
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FALSE
		 || second_operand->type != nd_step_type::AXIOM)
			return false;
		if (operand_states[1]->formula->type == FormulaType::NOT) {
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
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FALSE
		 || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		out.formula = second_operand->formula;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::UNIVERSAL_INTRODUCTION:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || operand_states[0]->formula == NULL || second_operand->type != nd_step_type::PARAMETER) {
			return false;
		} else if (operand_states[0]->assumptions_have_parameter(second_operand->parameter)) {
			/* the parameter is not allowed to occur free in the assumptions */
			return false;
		}

		parameter = Formula::new_parameter(second_operand->parameter);
		variable = Formula::new_variable(1);
		formula = substitute<TermType::PARAMETER, 1>(operand_states[0]->formula, parameter, variable);
		free(*parameter); if (parameter->reference_count == 0) free(parameter);
		free(*variable); if (variable->reference_count == 0) free(variable);
		if (formula == NULL) return false;
		out.formula = Formula::new_for_all(1, formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize(out.formula, canonicalizer);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::UNIVERSAL_ELIMINATION:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FOR_ALL
		 || second_operand->type != nd_step_type::TERM_PARAMETER)
			return false;
		variable = Formula::new_variable(operand_states[0]->formula->quantifier.variable);
		max_variable = max_bound_variable(*operand_states[0]->formula->quantifier.operand);
		if (max_variable == 0) {
			formula = second_operand->term;
			formula->reference_count++;
		} else {
			formula = shift_bound_variables(second_operand->term, max_variable);
			if (formula == NULL) {
				free(*variable); if (variable->reference_count == 0) free(variable);
				return false;
			}
		}
		out.formula = substitute<TermType::VARIABLE, -1>(operand_states[0]->formula->quantifier.operand, variable, formula);
		free(*formula); if (formula->reference_count == 0) free(formula);
		free(*variable); if (variable->reference_count == 0) free(variable);
		if (out.formula == NULL) return false;
		out.formula = try_canonicalize(out.formula, canonicalizer);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		if (operand_count != 2) return false;
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (second_operand->type == nd_step_type::ARRAY_PARAMETER) {
			variable = Formula::new_variable(1);
			formula = substitute<1>(operand_states[0]->formula, second_operand->parameters.data, second_operand->parameters.length, variable);
			free(*variable); if (variable->reference_count == 0) free(variable);
		} else {
			return false;
		}

		if (formula == NULL) return false;
		out.formula = Formula::new_exists(1, formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize(out.formula, canonicalizer);
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
			max_variable = max_bound_variable(*operand_states[1]->formula);
			if (max_variable == 0) {
				formula = operand_states[0]->formula->binary.left;
				formula->reference_count++;
			} else {
				formula = shift_bound_variables(operand_states[0]->formula->binary.left, max_variable);
				if (formula == NULL) return false;
			}
			out.formula = substitute(operand_states[1]->formula,
					proof.operands[2]->parameters.data, proof.operands[2]->parameters.length,
					formula, operand_states[0]->formula->binary.right);
			free(*formula); if (formula->reference_count == 0) free(formula);
			if (out.formula == NULL) {
				if (max_variable == 0) {
					formula = operand_states[0]->formula->binary.right;
					formula->reference_count++;
				} else {
					formula = shift_bound_variables(operand_states[0]->formula->binary.right, max_variable);
					if (formula == NULL) return false;
				}
				out.formula = substitute(operand_states[1]->formula,
						proof.operands[2]->parameters.data, proof.operands[2]->parameters.length,
						formula, operand_states[0]->formula->binary.left);
				free(*formula); if (formula->reference_count == 0) free(formula);
				if (out.formula == NULL) return false;
			}
		} else {
			return false;
		}

		out.formula = try_canonicalize(out.formula, canonicalizer);
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

template<typename Formula, typename Visitor>
bool visit(const nd_step<Formula>* proof, Visitor& visitor)
{
	array<const nd_step<Formula>*> stack(64);
	hash_set<const nd_step<Formula>*> visited(128);
	if (!stack.add(proof)) return false;
	while (stack.length > 0)
	{
		const nd_step<Formula>* node = stack.pop();
		if (!visited.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula>* const* operands;
		node->get_subproofs(operands, operand_count);

		if (!visit_node(node, visitor)) return false;

		if (operand_count == 0) continue;
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;

			if (visited.contains(operands[i]))
				continue;
			if (!stack.add(operands[i]))
				return false;
		}
	}
	return true;
}

template<typename Collection>
struct axiom_collector {
	Collection& axioms;
};

template<typename Formula, typename Collection>
inline bool visit_node(const nd_step<Formula>* proof, axiom_collector<Collection>& collector) {
	if (proof->type == nd_step_type::AXIOM)
		return collector.axioms.add(proof->formula);
	return true;
}

template<typename Formula, typename Collection>
inline bool get_axioms(const nd_step<Formula>* proof, Collection& axioms) {
	axiom_collector<Collection> collector = {axioms};
	return visit(proof, collector);
}

template<typename Formula, typename... ProofMap>
bool compute_in_degrees(const nd_step<Formula>* proof,
		hash_map<const nd_step<Formula>*, unsigned int>& in_degrees,
		ProofMap&&... proof_map)
{
	array<const nd_step<Formula>*> stack(64);
	hash_set<const nd_step<Formula>*> visited(128);
	if (!stack.add(map(proof, std::forward<ProofMap>(proof_map)...)))
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

			const nd_step<Formula>* operand = map(
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

template<typename BuiltInPredicates, typename Formula, typename Canonicalizer, typename... ProofMap>
Formula* compute_proof_conclusion(const nd_step<Formula>& proof,
		Canonicalizer& canonicalizer, ProofMap&&... proof_map)
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
			const nd_step<Formula>* operand = map(
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
			const nd_step<Formula>* operand = map(
					operands[i], std::forward<ProofMap>(proof_map)...);
			operand_states[i] = &proof_states.get(operand, contains);
			if (!operand->is_parameter() && !contains) {
				fprintf(stderr, "check_proof ERROR: The proof is not topologically ordered.\n");
				free_proof_states(proof_states); return NULL;
			}
			operand_states.length++;
		}

		/* check this proof step */
		if (!check_proof<BuiltInPredicates>(state, *node, operand_states.data, operand_states.length, canonicalizer, std::forward<ProofMap>(proof_map)...)) {
			free_proof_states(proof_states);
			return NULL;
		}
		operand_states.clear();
	}

	/* get the proof state of the last deduction step */
	bool contains;
	const proof_state<Formula>& root_state = proof_states.get(
			map(&proof, std::forward<ProofMap>(proof_map)...), contains);
	if (!contains) {
		fprintf(stderr, "check_proof ERROR: Unable to find proof state of root.\n");
		free_proof_states(proof_states); return NULL;
	}
	Formula* formula = root_state.formula;
	formula->reference_count++;
	free_proof_states(proof_states);
	return formula;
}

template<typename BuiltInPredicates, typename Formula, typename Canonicalizer, typename... ProofMap>
bool check_proof(const nd_step<Formula>& proof,
		const Formula* expected_conclusion,
		Canonicalizer& canonicalizer, ProofMap&&... proof_map)
{
	Formula* actual_conclusion = compute_proof_conclusion<BuiltInPredicates>(proof, canonicalizer, std::forward<ProofMap>(proof_map)...);
	if (actual_conclusion == NULL) return false;
	bool success = (*actual_conclusion == *expected_conclusion);
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

template<typename Formula>
struct natural_deduction
{
	typedef nd_step<Formula> Proof;
	typedef nd_step_type ProofType;
	typedef typename Formula::Term Term;

	static inline Proof* new_axiom(Formula* axiom) {
		return new_parameterized_step<nd_step_type::AXIOM>(axiom);
	}

	static inline Proof* new_beta(Formula* first, Formula* second) {
		return new_binary_step<nd_step_type::BETA_EQUIVALENCE>(new_formula_parameter(first), new_formula_parameter(second));
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

	static inline Proof* new_disjunction_intro(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION>(proof, new_formula_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_left(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_LEFT>(proof, new_formula_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_right(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT>(proof, new_formula_parameter(parameter));
	}

	template<template<typename> class Array>
	static inline Proof* new_disjunction_elim(Proof* disjunction, Array<Proof*> operands) {
		return new_array_step<nd_step_type::DISJUNCTION_ELIMINATION, 3>(make_composed_array_view(disjunction, operands));
	}

	static inline Proof* new_implication_intro(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_step_type::AXIOM) return NULL;
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
		if (assumption == NULL || assumption->type != nd_step_type::AXIOM) return NULL;
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

	static inline Proof* new_existential_intro(Proof* proof, unsigned int* term_indices, unsigned int term_index_count) {
		return new_binary_step<nd_step_type::EXISTENTIAL_INTRODUCTION>(proof, new_array_parameter(make_array_view(term_indices, term_index_count)));
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
		if (proof == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		step->operands[0] = proof;
		step->operands[1] = NULL;
		step->operands[2] = NULL;
		proof->reference_count++;
		if (!proof->children.add(step)) {
			free(*step); free(step);
			return NULL;
		}
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_binary_step(Proof* left, Proof* right) {
		if (left == NULL || right == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		step->operands[0] = left;
		step->operands[1] = right;
		step->operands[2] = NULL;
		left->reference_count++;
		right->reference_count++;
		if (!left->children.add(step) || !right->children.add(step)) {
			free(*step); free(step);
			return NULL;
		}
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_ternary_step(Proof* first, Proof* second, Proof* third) {
		if (first == NULL || second == NULL || third == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		step->operands[0] = first;
		step->operands[1] = second;
		step->operands[2] = third;
		first->reference_count++;
		second->reference_count++;
		third->reference_count++;
		if (!first->children.add(step) || !second->children.add(step) || !third->children.add(step)) {
			free(*step); free(step);
			return NULL;
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
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		if (!array_init(step->operand_array, sizeof...(Args))) {
			free(step); return NULL;
		}
		if (!new_array_step_helper<0>(step, std::forward<Args>(args)...)) {
			free(step->operand_array); free(step);
			return NULL;
		}
		step->operand_array.length = sizeof...(Args);
		return step;
	}

	template<nd_step_type Type, unsigned int MinOperandCount, typename Array>
	static inline Proof* new_array_step(const Array& operands) {
		if (size(operands) < MinOperandCount) {
			fprintf(stderr, "natural_deduction.new_array_step ERROR: "
					"This proof step requires at least two operands.\n");
			return NULL;
		}

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		if (!array_init(step->operand_array, size(operands))) {
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < size(operands); i++) {
			step->operand_array[i] = operands[i];
			operands[i]->reference_count++;
			if (!operands[i]->children.add(step)) {
				for (unsigned int j = 0; j < i; j++) {
					operands[j]->remove_child(step);
					free(*operands[j]);
				}
				free(step->operand_array); free(step);
				return NULL;
			}
		}
		step->operand_array.length = size(operands);
		return step;
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
			const nd_step<Formula>* operand = map(
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
	typename UniversalIntroductions,
	typename UniversalEliminations,
	typename TermIndices,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula>* const* proof,
		unsigned int step_index,
		unsigned int& formula_counter,
		array<unsigned int>& available_parameters,
		array_multiset<Formula*, false>& axioms,
		ConjunctionIntroductions& conjunction_introductions,
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
		if (!axioms.add_unsorted(current_step.formula)) exit(EXIT_FAILURE);
		return 0.0;
	case nd_step_type::BETA_EQUIVALENCE:
		/* TODO: this is a hack; correct it */
		formula_counter++;
		return -LOG_ND_RULE_COUNT;
	case nd_step_type::FALSITY_ELIMINATION:
		/* TODO: this is a hack; correct it */
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		/* TODO: implement this */
		/* TODO: this is a hack; correct it */
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
		operand = map(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		/* we need to compute the prior on the parameter */
		if (!universal_introductions.add(make_pair(operand->parameter, available_parameters)))
			exit(EXIT_FAILURE);
		
		index = available_parameters.index_of(operand->parameter);
		if (index < available_parameters.length)
			shift_left(available_parameters.data + index, available_parameters.length - index - 1);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::UNIVERSAL_ELIMINATION:
		operand = map(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand->term->type == TermType::PARAMETER) {
			if (available_parameters.ensure_capacity(available_parameters.length + 1)) exit(EXIT_FAILURE);
			add_sorted<true>(available_parameters, operand->term->parameter);
		}
		/* we need to compute the prior on the term */
		if (!universal_eliminations.add(*operand->term)) exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		/* we need to compute the prior on the parameter (it's a list of term indices) */
		operand = map(current_step.operands[1], std::forward<ProofMap>(proof_map)...);
		if (!term_indices.add(make_array_view(operand->parameters.data, operand->parameters.length)))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::EQUALITY_ELIMINATION:
		/* we need to compute the prior on the parameter (it's a list of term indices) */
		operand = map(current_step.operands[2], std::forward<ProofMap>(proof_map)...);
		if (!term_indices.add(make_array_view(operand->parameters.data, operand->parameters.length)))
			exit(EXIT_FAILURE);
		return -LOG_ND_RULE_COUNT - 2*log_cache<double>::instance().get(formula_counter++);
	case nd_step_type::COUNT:
		break;
	}
	fprintf(stderr, "log_probability ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior>
struct canonicalized_proof_prior
{
	AxiomPrior axiom_prior;
	ConjunctionIntroductionPrior conjunction_introduction_prior;
	UniversalIntroductionPrior universal_introduction_prior;
	UniversalEliminationPrior universal_elimination_prior;
	TermIndicesPrior term_indices_prior;
	ProofLengthPrior proof_length_prior;

	canonicalized_proof_prior(
			const AxiomPrior& axiom_prior,
			const ConjunctionIntroductionPrior& conjunction_introduction_prior,
			const UniversalIntroductionPrior& universal_introduction_prior,
			const UniversalEliminationPrior& universal_elimination_prior,
			const TermIndicesPrior& term_indices_prior,
			const ProofLengthPrior& proof_length_prior) :
		axiom_prior(axiom_prior),
		conjunction_introduction_prior(conjunction_introduction_prior),
		universal_introduction_prior(universal_introduction_prior),
		universal_elimination_prior(universal_elimination_prior),
		term_indices_prior(term_indices_prior),
		proof_length_prior(proof_length_prior)
	{ }

	canonicalized_proof_prior(const canonicalized_proof_prior& src) :
		axiom_prior(src.axiom_prior),
		conjunction_introduction_prior(src.conjunction_introduction_prior),
		universal_introduction_prior(src.universal_introduction_prior),
		universal_elimination_prior(src.universal_elimination_prior),
		term_indices_prior(src.term_indices_prior),
		proof_length_prior(src.proof_length_prior)
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
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior>
inline canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>
make_canonicalized_proof_prior(
		const AxiomPrior& axiom_prior,
		const ConjunctionIntroductionPrior& conjunction_introduction_prior,
		const UniversalIntroductionPrior& universal_introduction_prior,
		const UniversalEliminationPrior& universal_elimination_prior,
		const TermIndicesPrior& term_indices_prior,
		const ProofLengthPrior& proof_length_prior)
{
	return canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>(
			axiom_prior, conjunction_introduction_prior, universal_introduction_prior, universal_elimination_prior, term_indices_prior, proof_length_prior);
}

template<typename Formula,
	typename ConjunctionIntroductions,
	typename UniversalIntroductions,
	typename UniversalEliminations,
	typename TermIndices,
	typename ProofLengthPrior,
	typename... ProofMap>
inline double log_probability_helper(
		const nd_step<Formula>& proof,
		array_multiset<Formula*, false>& axioms,
		ConjunctionIntroductions& conjunction_introductions,
		UniversalIntroductions& universal_introductions,
		UniversalEliminations& universal_eliminations,
		TermIndices& term_indices,
		ProofLengthPrior& proof_length_prior,
		ProofMap&&... proof_map)
{
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
		value += log_probability(canonical_order.data, i, formula_counter,
				available_parameters, axioms, conjunction_introductions, universal_introductions,
				universal_eliminations, term_indices, std::forward<ProofMap>(proof_map)...);
	}
	return value;
}

template<
	typename Formula, typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula>& proof,
		canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior,
		ProofMap&&... proof_map)
{
	typedef typename ConjunctionIntroductionPrior::ObservationCollection ConjunctionIntroductions;
	typedef typename UniversalIntroductionPrior::ObservationCollection UniversalIntroductions;
	typedef typename UniversalEliminationPrior::ObservationCollection UniversalEliminations;
	typedef typename TermIndicesPrior::ObservationCollection TermIndices;

	array_multiset<Formula*, false> axioms(16);
	ConjunctionIntroductions conjunction_introductions;
	UniversalIntroductions universal_introductions;
	UniversalEliminations universal_eliminations;
	TermIndices term_indices;
	double value = log_probability_helper(proof, axioms,
			conjunction_introductions, universal_introductions, universal_eliminations, prior.log_continue_probability,
			prior.log_stop_probability, term_indices, prior.proof_length_prior, std::forward<ProofMap>(proof_map)...);
	return value + log_probability(axioms, prior.axiom_prior)
		 + log_probability(conjunction_introductions, prior.conjunction_introduction_prior)
		 + log_probability(universal_introductions, prior.universal_introduction_prior)
		 + log_probability(universal_eliminations, prior.universal_elimination_prior)
		 + log_probability(term_indices, prior.term_indices_prior);
}

template<
	typename Formula, typename AxiomPrior,
	typename ConjunctionIntroductionPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename TermIndicesPrior,
	typename ProofLengthPrior,
	typename MultisetType,
	typename... AxiomPriorParameters>
double log_probability_ratio(
		const array_map<nd_step<Formula>*, proof_substitution<Formula>>& proofs,
		canonicalized_proof_prior<AxiomPrior, ConjunctionIntroductionPrior, UniversalIntroductionPrior, UniversalEliminationPrior, TermIndicesPrior, ProofLengthPrior>& prior,
		const MultisetType& proof_axioms,
		array_multiset<Formula*, false>& old_axioms,
		array_multiset<Formula*, false>& new_axioms,
		AxiomPriorParameters&&... axiom_prior_parameters)
{
	typedef typename ConjunctionIntroductionPrior::ObservationCollection ConjunctionIntroductions;
	typedef typename UniversalIntroductionPrior::ObservationCollection UniversalIntroductions;
	typedef typename UniversalEliminationPrior::ObservationCollection UniversalEliminations;
	typedef typename TermIndicesPrior::ObservationCollection TermIndices;

	double value = 0.0;
	for (const auto& entry : proofs) {
		ConjunctionIntroductions old_conjunction_introductions, new_conjunction_introductions;
		UniversalIntroductions old_universal_introductions, new_universal_introductions;
		UniversalEliminations old_universal_eliminations, new_universal_eliminations;
		TermIndices old_term_indices, new_term_indices;
		value -= log_probability_helper(*entry.key, old_axioms, old_conjunction_introductions,
				old_universal_introductions, old_universal_eliminations, old_term_indices, prior.proof_length_prior);
		value += log_probability_helper(*entry.key, new_axioms, new_conjunction_introductions,
				new_universal_introductions, new_universal_eliminations, new_term_indices, prior.proof_length_prior, entry.value);

		value += log_probability_ratio(old_conjunction_introductions, new_conjunction_introductions, prior.conjunction_introduction_prior);
		value += log_probability_ratio(old_universal_introductions, new_universal_introductions, prior.universal_introduction_prior);
		value += log_probability_ratio(old_universal_eliminations, new_universal_eliminations, prior.universal_elimination_prior);
		value += log_probability_ratio(old_term_indices, new_term_indices, prior.term_indices_prior);
	}

	sort(old_axioms.counts.keys, old_axioms.counts.values, old_axioms.counts.size);
	sort(new_axioms.counts.keys, new_axioms.counts.values, new_axioms.counts.size);
	return value + log_probability_ratio(proof_axioms, old_axioms, new_axioms,
			prior.axiom_prior, std::forward<AxiomPriorParameters>(axiom_prior_parameters)...);
}

#endif /* NATURAL_DEDUCTION_H_ */
