#ifndef NATURAL_DEDUCTION_H_
#define NATURAL_DEDUCTION_H_

#include <core/array.h>
#include <math/log.h>

using namespace core;

enum class nd_step_type : uint_fast16_t
{
	AXIOM = 0,
	PARAMETER,
	ARRAY_PARAMETER,
	TERM_PARAMETER,
	FORMULA_PARAMETER,

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
	PROOF_BY_CONTRADICTION,
	NEGATION_ELIMINATION,

	UNIVERSAL_INTRODUCTION,
	UNIVERSAL_ELIMINATION,
	EXISTENTIAL_INTRODUCTION,
	EXISTENTIAL_ELIMINATION,

	COUNT
};

constexpr static unsigned int ND_OPERAND_COUNT = 3;
constexpr static double LOG_ND_RULE_COUNT = log((double) nd_step_type::COUNT);

template<typename Formula, bool Canonical>
struct nd_step
{
	typedef typename Formula::Term Term;

	nd_step_type type;
	unsigned int reference_count;
	union {
		Term term;
		Formula* formula;
		unsigned int parameter;
		array<unsigned int> parameters;

		nd_step<Formula, Canonical>* operands[ND_OPERAND_COUNT];
		array<nd_step<Formula, Canonical>*> operand_array;
	};

	/* list of proof steps that use the subproof rooted at this rule instance */
	array<nd_step<Formula, Canonical>*> children;

	static inline void free(nd_step<Formula, Canonical>& step) {
		core::free(step.children);
		step.free();
	}

	inline void get_subproofs(nd_step<Formula, Canonical>**& subproofs, unsigned int& length) {
		switch (type) {
		case TERM_PARAMETER:
		case ARRAY_PARAMETER:
		case AXIOM:
		case FORMULA_PARAMETER:
			subproofs = NULL; length = 0;
			return;
		case CONJUNCTION_INTRODUCTION:
			subproofs = operand_array.data;
			length = operand_array.length;
			return;
		case CONJUNCTION_ELIMINATION:
		case CONJUNCTION_ELIMINATION_LEFT:
		case CONJUNCTION_ELIMINATION_RIGHT:
		case DISJUNCTION_INTRODUCTION:
		case DISJUNCTION_INTRODUCTION_LEFT:
		case DISJUNCTION_INTRODUCTION_RIGHT:
		case DISJUNCTION_ELIMINATION:
		case IMPLICATION_INTRODUCTION:
		case IMPLICATION_ELIMINATION:
		case BICONDITIONAL_INTRODUCTION:
		case BICONDITIONAL_ELIMINATION_LEFT:
		case BICONDITIONAL_ELIMINATION_RIGHT:
		case PROOF_BY_CONTRADICTION:
		case NEGATION_ELIMINATION:
		case UNIVERSAL_INTRODUCTION:
		case UNIVERSAL_ELIMINATION:
		case EXISTENTIAL_INTRODUCTION:
		case EXISTENTIAL_ELIMINATION:
			subproofs = operands;
			length = ND_OPERAND_COUNT;
			return true;
		}
		fprintf(stderr, "nd_step.get_subproofs ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

private:
	inline void free() {
		switch (type) {
		case TERM_PARAMETER:
			free(term); return;
		case ARRAY_PARAMETER:
			free(parameters); return;
		case AXIOM:
		case FORMULA_PARAMETER:
			core::free(*formula);
			if (formula->reference_count == 0)
				core::free(formula);
			return;
		case CONJUNCTION_INTRODUCTION:
			for (unsigned int i = 0; i < operand_array.length; i++) {
				operand_array[i]->remove_child(this);
				core::free(*operand_array[i]);
				if (operand_array[i]->reference_count == 0)
					core::free(operand_array[i]);
			}
			core::free(operand_array);
			return;
		case CONJUNCTION_ELIMINATION:
		case CONJUNCTION_ELIMINATION_LEFT:
		case CONJUNCTION_ELIMINATION_RIGHT:
		case DISJUNCTION_INTRODUCTION:
		case DISJUNCTION_INTRODUCTION_LEFT:
		case DISJUNCTION_INTRODUCTION_RIGHT:
		case DISJUNCTION_ELIMINATION:
		case IMPLICATION_INTRODUCTION:
		case IMPLICATION_ELIMINATION:
		case BICONDITIONAL_INTRODUCTION:
		case BICONDITIONAL_ELIMINATION_LEFT:
		case BICONDITIONAL_ELIMINATION_RIGHT:
		case PROOF_BY_CONTRADICTION:
		case NEGATION_ELIMINATION:
		case UNIVERSAL_INTRODUCTION:
		case UNIVERSAL_ELIMINATION:
		case EXISTENTIAL_INTRODUCTION:
		case EXISTENTIAL_ELIMINATION:
			for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
				if (operands[i] == NULL) break;
				operands[i]->remove_child(this);
				core::free(*operands[i]);
				if (operands[i]->reference_count == 0)
					core::free(operands[i]);
			}
			return;
		}
		fprintf(stderr, "nd_step.free ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}

	inline void remove_child(const nd_step<Formula, Canonical>* child) {
		unsigned int index = children.index_of(child);
#if !defined(NDEBUG)
		if (index == children.length) {
			fprintf(stderr, "nd_step.remove_child WARNING: Index out of bounds.\n");
			return;
		}
#endif
		children.remove(index);
	}
};

template<typename Formula, bool Canonical>
inline int_fast8_t compare(const nd_step<Formula, Canonical>& first, const nd_step<Formula, Canonical>& second)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;

	int_fast8_t result;
	switch (first.type) {
	case TERM_PARAMETER:
		return compare(first.term, second.term);
	case ARRAY_PARAMETER:
		if (first.parameters.length < second.parameters.length) return -1;
		else if (first.parameters.length > second.parameters.length) return 1;
		for (unsigned int i = 0; i < first.parameters.length; i++) {
			if (first.parameters[i] < second.parameters[i]) return -1;
			else if (first.parameters[i] > second.parameters[i]) return 1;
		}
		return 0;
	case AXIOM:
	case FORMULA_PARAMETER:
		return compare(*first.formula, *second.formula);
	case CONJUNCTION_INTRODUCTION:
		if (first.operand_array.length < second.operand_array.length) return -1;
		else if (first.operand_array.length > second.operand_array.length) return 1;
		for (unsigned int i = 0; i < first.operand_array.length; i++) {
			int_fast8_t result = compare(*first.operand_array[i], *second.operand_array[i]);
			if (result != 0) return result;
		}
		return 0;
	case CONJUNCTION_ELIMINATION:
	case CONJUNCTION_ELIMINATION_LEFT:
	case CONJUNCTION_ELIMINATION_RIGHT:
	case DISJUNCTION_INTRODUCTION:
	case DISJUNCTION_INTRODUCTION_LEFT:
	case DISJUNCTION_INTRODUCTION_RIGHT:
	case DISJUNCTION_ELIMINATION:
	case IMPLICATION_INTRODUCTION:
	case IMPLICATION_ELIMINATION:
	case BICONDITIONAL_INTRODUCTION:
	case BICONDITIONAL_ELIMINATION_LEFT:
	case BICONDITIONAL_ELIMINATION_RIGHT:
	case PROOF_BY_CONTRADICTION:
	case NEGATION_ELIMINATION:
	case UNIVERSAL_INTRODUCTION:
	case UNIVERSAL_ELIMINATION:
	case EXISTENTIAL_INTRODUCTION:
	case EXISTENTIAL_ELIMINATION:
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
	}
	fprintf(stderr, "compare ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula, bool Canonical>
inline bool operator < (const nd_step<Formula, Canonical>& first, const nd_step<Formula, Canonical>& second) {
	return compare(first, second) < 0;
}

template<typename Formula>
struct proof_state {
	array<Formula*> assumptions;
	Formula* formula;

	proof_state() : assumptions(8) { }
	~proof_state() { free(); }

	inline unsigned int index_of_assumption(Formula* assumption) const {
		for (unsigned int i = 0; i < assumptions.length; i++)
			if (*assumptions[i] == *assumption) return i;
		return assumptions.length;
	}

	inline bool assumptions_have_parameter(unsigned int parameter) const {
		for (const Formula* assumption : assumptions)
			if (contains_parameter(*assumption, parameter)) return true;
		return false;
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

template<template<typename> class Array, typename T>
struct exclude {
	Array<T>& elements;
	T& excluded;

	exclude(Array<T>& elements, T& excluded) : elements(elements), excluded(excluded) { }
};

template<template<typename> class Array, typename T>
inline exclude<Array<T>, T> make_exclude(Array<T>& elements, T& excluded) {
	return exclude<Array<T>, T>(elements, excluded);
} 

template<typename T>
bool pass_hypotheses(array<T>& dst, const array<T>& first) {
	return dst.append(first.data, first.length);
}

template<template<typename> class Array, typename T>
bool pass_hypotheses(array<T>& dst, const exclude<Array, T>& first)
{
	if (!dst.ensure_capacity(dst.length + first.elements.length))
		return false;

	bool discharged = false;
	for (unsigned int i = 0; i < first.elements.length; i++) {
		if (*first.elements[i] == *first.excluded) { discharged = true; continue; }
		dst[dst.length] = first.elements[i];
		dst.length++;
	}
	return discharged;
}

template<template<typename> class Array, typename T, bool RemoveDuplicates = true>
bool pass_hypotheses(array<T>& dst, const array<T>& first, const exclude<Array, T>& second)
{
	if (!dst.ensure_capacity(first.length + second.elements.length))
		return false;

	bool discharged = false;
	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.elements.length)
	{
		if (*second.excluded == *second.elements[j]) { j++; discharged = true; continue; }

		if (first[i] == second.elements[j]) {
			set_union_helper<RemoveDuplicates>(dst, dst.length, first[i]);
			i++; j++;
		} else if (first[i] < second.elements[j]) {
			set_union_helper<RemoveDuplicates>(dst, dst.length, first[i]);
			i++;
		} else {
			set_union_helper<RemoveDuplicates>(dst, dst.length, second.elements[j]);
			j++;
		}
	}

	while (i < first.length) {
		set_union_helper<RemoveDuplicates>(dst, dst.length, first[i]);
		i++;
	} while (j < second.elements.length) {
		if (*discharge_second == *second.elements[j]) { j++; discharged = true; continue; }
		set_union_helper<RemoveDuplicates>(dst, dst.length, second.elements[j]);
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

template<template<typename> typename Array, typename T, typename... Lists>
bool pass_hypotheses(array<T>& out, const exclude<Array, T>& first, Lists&&... lists) {
	if (!pass_hypotheses(out, std::forward<Lists>(lists)...)) return false;
	array<T> temp = array<T>(out.length + first.elements.length);
	swap(temp, out);
	return pass_hypotheses(out, temp, first);
}

template<template<typename> typename Array, typename T>
bool pass_hypotheses(array<T>& out, const Array<array<T>>& lists) {
	if (lists.length == 0)
		return true;
	if (!pass_hypotheses(out, lists[0])) return false;

	array<T> temp = array<T>(16);
	for (unsigned int i = 1; i < lists.length; i++) {
		if (!temp.ensure_capacity(out.length + lists[i].length))
			return false;
		swap(temp, out);
		if (!set_union(out, temp, lists[i]))
			return false;
		temp.clear();
	}
	return true;
}

template<typename Formula>
bool subtract_unifying_hypotheses(
	array<Formula*>& dst_hypotheses,
	const array<Formula*>& hypotheses,
	const Formula& src, const Formula& C)
{
	for (Formula* hypothesis : hypotheses) {
		unsigned int parameter;
		if (unifies_parameter(*src.quantifier.operand, *hypothesis, Formula::new_variable(src.quantifier.variable), parameter)) {
			if (contains_parameter(C, parameter)) {
				if (!dst_hypotheses.add(hypothesis)) return false;
			}
		} else {
			if (!dst_hypotheses.add(hypothesis)) return false;
		}
	}
	return true;
}

template<bool Canonical, typename Formula>
Formula* try_canonicalize(Formula* formula) {
	if (Canonical) {
		Formula* out = canonicalize(*formula);
		free(*formula); if (formula->reference_count == 0) free(formula);
		return out;
	} else {
		return formula;
	}
}

template<typename Formula>
struct proof_state_formulas {
	proof_state<Formula>** states;
	unsigned int length;

	proof_state_formulas(proof_state<Formula>** states, unsigned int length) : states(states), length(length) { }

	inline Formula* operator[] (size_t index) {
		return states[index]->formula;
	}
};

template<typename Formula>
struct proof_state_assumptions {
	proof_state<Formula>** states;
	unsigned int length;

	proof_state_assumptions(proof_state<Formula>** states, unsigned int length) : states(states), length(length) { }

	inline array<Formula*>& operator[] (size_t index) {
		return states[index]->assumptions;
	}
};

template<typename T>
struct array_view {
	T* array;
	unsigned int length;

	array_view(T* array, unsigned int length) : array(array), length(length) { }

	inline T& operator[] (size_t index) {
		return array[index];
	}
};

template<typename T>
struct indexed_array_view {
	T* array;
	unsigned int* indices;
	unsigned int length;

	indexed_array_view(T* array, unsigned int* indices, unsigned int length) : array(array), indices(indices), length(length) { }

	inline T& operator[] (size_t index) {
		return array[indices[index]];
	}
};

template<typename Formula, bool Canonical>
struct proof_substitution {
	const nd_step<Formula, Canonical>* old_step;
	const nd_step<Formula, Canonical>* new_step;
};

template<typename Formula, bool Canonical>
inline const nd_step<Formula, Canonical>* map(const nd_step<Formula, Canonical>* step) {
	return step;
}

template<typename Formula, bool Canonical>
inline const nd_step<Formula, Canonical>* map(const nd_step<Formula, Canonical>* step,
		const proof_substitution<Formula, Canonical>& substitution)
{
	if (step == substitution.old_step) return new_step;
	else return step;
}

template<typename Formula, bool Canonical, typename... ProofMap>
bool check_proof(proof_state<Formula>& out,
		const nd_step<Formula, Canonical>& proof,
		const proof_state<Formula>** operand_states,
		unsigned int operand_count,
		ProofMap&&... proof_map)
{
	typedef typename Formula::Type FormulaType;

	Formula* formula; Formula* second_operand; unsigned int parameter;
	switch (proof.type) {
	case nd_step_type::AXIOM:
		if (Canonical && !is_canonical(*proof.formula)) {
			fprintf(stderr, "check_proof ERROR: Axiom is not in canonical form.\n");
			return false;
		}
		out.formula = proof.formula;
		out.formula->reference_count++;
		return out.assumptions.add(out.formula);
	case nd_step_type::CONJUNCTION_INTRODUCTION:
		if (operand_count < 2) return false;
		out.formula = Formula::new_and(proof_state_formulas(operand_states, operand_count));
		if (out.formula == NULL) return false;
		for (unsigned int i = 0; i < operand_count; i++)
			operand_states[i]->formula->reference_count++;
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, proof_state_assumptions(operand_states, operand_count));
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
		if (!Canonical && proof.type == nd_step_type::CONJUNCTION_ELIMINATION) return false;
		if (Canonical && proof.type == nd_step_type::CONJUNCTION_ELIMINATION_LEFT && nd_step_type::CONJUNCTION_ELIMINATION_RIGHT) return false;
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
		if (proof.type == CONJUNCTION_ELIMINATION_LEFT) {
			out.formula = formula->array.operands[0];
			out.formula->reference_count++;
		} else if (proof.type == CONJUNCTION_ELIMINATION_RIGHT) {
			out.formula = formula->array.operands[1];
			out.formula->reference_count++;
		} else if (proof.type == CONJUNCTION_ELIMINATION) {
			const indexed_array_view<Formula*> conjuncts(formula->array.operands, second_operand->parameters.data, second_operand->parameters.length);
			if (std::max_element(conjuncts.indices, conjuncts.indices + conjuncts.length) >= formula->array.length
			 || !std::is_sorted(conjuncts.indices, conjuncts.indices + conjuncts.length)
			 || std::adjacent_find(conjuncts.indices, conjuncts.indices + conjuncts.length) != conjuncts.indices + conjuncts.length)
				return false;
			out.formula = Formula::new_and(conjuncts);
			if (out.formula == NULL) return false;
		}
		free(*formula);
		if (formula->reference_count == 0)
			free(formula);
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_INTRODUCTION:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
		if (!Canonical && proof.type == nd_step_type::DISJUNCTION_INTRODUCTION) return false;
		if (Canonical && proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_LEFT && nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT) return false;
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || second_operand->type != nd_step_type::FORMULA_PARAMETER)
			return false;
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_LEFT || proof.type == nd_step_type::DISJUNCTION_INTRODUCTION)
			out.formula = Formula::new_or(operand_states[0]->formula, second_operand->formula);
		if (proof.type == nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT)
			out.formula = Formula::new_or(second_operand->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		operand_states[1]->formula->reference_count++;
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::DISJUNCTION_ELIMINATION:
		if (operand_count != 3 || operand_states[0]->formula->type != FormulaType::OR
		 || *operand_states[1]->formula != *operand_states[2]->formula)
			return false;
		out.formula = operand_states[0]->formula;
		out.formula.reference_count++;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions,
				make_exclude(operand_states[1]->assumptions, operand_states[0]->formula.binary.left),
				make_exclude(operand_states[2]->assumptions, operand_states[0]->formula.binary.right));
	case nd_step_type::IMPLICATION_INTRODUCTION:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || second_operand->type != nd_step_type::AXIOM)
			return false;
		out.formula = Formula::new_if_then(operand_states[1]->formula, operand_states[0]->formula);
		if (out.formula == NULL) return false;
		operand_states[0]->formula->reference_count++;
		operand_states[1]->formula->reference_count++;
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, make_exclude(operand_states[0]->assumptions, operand_states[1]->formula));
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
		out.formula = try_canonicalize<Canonical>(out.formula);
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
		 || second_operand->type != nd_proof_type::axiom
		 || operand_states[1]->formula->type != FormulaType::NOT)
			return false;
		out.formula = second_operand->formula->unary;
		out.formula->reference_count++;
		return pass_hypotheses(out.assumptions, make_exclude(operand_states[0]->assumptions, second_operand->formula));
	case nd_step_type::NEGATION_ELIMINATION:
		if (operand_count != 2
		 || operand_states[1]->formula->type != FormulaType::NOT
		 || *operand_states[0]->formula != *operand_states[1]->formula->unary)
			return false;
		out.formula = Formula::new_false();
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions, operand_states[1]->assumptions);
	case nd_step_type::UNIVERSAL_INTRODUCTION:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2 || operand_states[0]->formula == NULL || second_operand->type != nd_proof_type::PARAMETER) {
			return false;
		} else if (operand_states[0]->assumptions_have_parameter(second_operand->parameter)) {
			/* the parameter is not allowed to occur free in the assumptions */
			return false;
		}

		formula = substitute<1>(*operand_states[0]->formula, Formula::new_parameter(second_operand->parameter), Formula::new_variable(1));
		if (formula == NULL) return false;
		out.formula = Formula::new_for_all(1, formula);
		if (out.formula == NULL) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::UNIVERSAL_ELIMINATION:
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (operand_count != 2
		 || operand_states[0]->formula->type != FormulaType::FOR_ALL
		 || second_operand->type != nd_proof_type::TERM_PARAMETER)
			return false;
		out.formula = substitute<-1>(*operand_states[0]->formula->quantifier.operand,
				Formula::new_variable(operand_states[0]->formula->quantifier.variable),
				second_operand->term);
		if (out.formula == NULL) return false;
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		if (operand_count != 2) return false;
		second_operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		if (second_operand->type == nd_proof_type::ARRAY_PARAMETER) {
			formula = substitute<1>(*operand_states[0]->formula, second_operand->parameters.data, second_operand->parameters.length, Formula::new_variable(1));
		} else {
			return false;
		}

		if (formula == NULL) return false;
		out.formula = Formula::new_exists(1, formula);
		if (out.formula == NULL) {
			free(*formula) if (formula->reference_count == 0) free(formula);
			return false;
		}
		out.formula = try_canonicalize<Canonical>(out.formula);
		if (out.formula == NULL) return false;
		return pass_hypotheses(out.assumptions, operand_states[0]->assumptions);
	case nd_step_type::EXISTENTIAL_ELIMINATION:
		if (operand_count != 2 || operand_states[0]->formula->type != FormulaType::EXISTS)
			return false;
		if (!subtract_unifying_hypotheses(out.assumptions, operand_states[1]->assumptions, *operand_states[1]->formula))
			return false;
		out.formula = operand_states[1]->formula;
		out.formula->reference_count++;
		return true;
	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
		break;
	}
	fprintf(stderr, "check_proof ERROR: Unexpected nd_step_type.\n");
	return false;
}

template<typename Formula, bool Canonical, typename... ProofMap>
bool compute_in_degrees(const nd_step<Formula, Canonical>* proof,
		hash_map<const nd_step<Formula, Canonical>*, unsigned int>& in_degrees,
		ProofMap&&... proof_map)
{
	array<const nd_step<Formula, Canonical>*> stack(64);
	hash_set<const nd_step<Formula, Canonical>*> visited(128);
	if (!stack.add(map(proof, std::forward<ProofMap>(proof_map)...)))
		return false;
	while (stack.length > 0)
	{
		const nd_step<Formula, Canonical>* node = stack.pop();
		if (!visited.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula, Canonical>* const* operands;
		node->get_subproofs(operands, operand_count);
		if (operand_count == 0) continue;
		if (!in_degrees.check_size(in_degrees.table.size + operand_count + 1))
			return false;

		bool contains; unsigned int bucket;
		unsigned int& degree = in_degrees.get(node, contains, bucket);
		if (!contains) {
			in_degrees.table.keys[bucket] = node;
			in_degrees.table.size++;
			degree = 0;
		}

		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;

			const nd_step<Formula, Canonical>* operand = map(
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
				return NULL;
		}
	}

	return true;
}

template<typename Formula, bool Canonical, typename... ProofMap>
Formula* check_proof(const nd_step<Formula, Canonical>& proof, ProofMap&&... proof_map)
{
	/* first list the proof steps in reverse topological order */
	hash_map<const nd_step<Formula, Canonical>*, unsigned int> in_degrees(128);
	if (!compute_in_degrees(&proof, in_degrees, std::forward<ProofMap>(proof_map)...)) return NULL;

	array<const nd_step<Formula, Canonical>*> stack(32);
	for (const auto& entry : in_degrees) {
		if (entry.value == 0 && !stack.add(entry.key))
			return NULL;
	}

	array<const nd_step<Formula, Canonical>*> topological_order(64);
	while (stack.length > 0) {
		const nd_step<Formula, Canonical>* node = stack.pop();
		if (!topological_order.add(node)) return NULL;

		unsigned int operand_count;
		const nd_step<Formula, Canonical>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			const nd_step<Formula, Canonical>* operand = map(
					operands[i], std::forward<ProofMap>(proof_map)...);
			unsigned int& degree = in_degrees.get(operand);
			degree--;

			if (degree == 0 && !stack.add(operand))
				return false;
		}
	}

	/* process the proof steps in reverse topological order */
	hash_map<const nd_step<Formula, Canonical>*, proof_state<Formula>> proof_states(128);
	array<proof_state<Formula>*> operand_states(16);
	for (unsigned int k = topological_order.length; k > 0; k--) {
		const nd_step<Formula, Canonical>* node = topological_order[k - 1];
		if (!proof_states.check_size()) return NULL;

		bool contains; unsigned int bucket;
		proof_state<Formula>& state = proof_states.get(node, contains, bucket);
		if (contains) {
			fprintf(stderr, "check_proof ERROR: The proof state at this node should be uninitialized.\n");
			free_proof_states(proof_states); return NULL;
		} else if (!init(state)) {
			fprintf(stderr, "check_proof ERROR: Unable to initialize new proof_state.\n");
			free_proof_states(proof_states); return NULL;
		}
		state.table.keys[bucket] = node;
		state.table.size++;

		/* get the proof states of the operands */
		unsigned int operand_count;
		const nd_step<Formula, Canonical>* const* operands;
		node->get_subproofs(operands, operand_count);
		if (!operand_states.ensure_capacity(operand_count)) {
			free_proof_states(proof_states); return NULL;
		}
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) break;
			const nd_step<Formula, Canonical>* operand = map(
					operands[i], std::forward<ProofMap>(proof_map)...);
			operand_states[i] = &proof_states.get(operand, contains);
			if (!contains) {
				fprintf(stderr, "check_proof ERROR: The proof is not topologically ordered.\n");
				free_proof_states(proof_states); return NULL;
			}
			operand_states.length++;
		}

		/* check this proof step */
		if (!check_proof(state, *node, operand_states.data, operand_states.length, std::forward<ProofMap>(proof_map)...)) {
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

template<typename Formula, bool Canonical, typename... ProofMap>
bool check_proof(const nd_step<Formula, Canonical>& proof,
		const Formula* expected_conclusion, ProofMap&&... proof_map)
{
	Formula* actual_conclusion = check_proof(proof, std::forward<ProofMap>(proof_map)...);
	bool success = (*actual_conclusion != *expected_conclusion)
	if (!success)
		fprintf(stderr, "check_proof ERROR: Actual concluding formula does not match the expected formula.\n");
	free(*actual_conclusion)
	if (actual_conclusion->reference_count == 0)
		free(actual_conclusion);
	return success;
}

template<typename Formula, bool Canonical>
bool new_nd_step(nd_step<Formula, Canonical>*& step, nd_step_type type)
{
	step = (nd_step<Formula, Canonical>*) malloc(sizeof(nd_step<Formula, Canonical>));
	if (step == NULL) {
		fprintf(stderr, "new_nd_step ERROR: Out of memory.\n");
		return false;
	} else if (!array_init(step->children, 4)) {
		free(step); return false;
	}
	step->type = type;
	return true;
}

template<typename Formula, bool Canonical>
struct natural_deduction
{
	typedef nd_step<Formula, Canonical> Proof;
	typedef typename Formula::Term Term;

	static inline Proof* new_axiom(Formula* axiom) {
		new_parameterized_step<nd_step_type::AXIOM>(axiom);
	}

	template<typename... Args>
	static inline Proof* new_conjunction_intro(Args&&... args) {
		return new_array_step<nd_step_type::CONJUNCTION_INTRODUCTION>(std::forward<Args>(args)...);
	}

	template<template<typename> class Array>
	static inline Proof* new_conjunction_intro(const Array<Proof>& operands) {
		if (!Canonical) return NULL;
		return new_array_step<nd_step_type::CONJUNCTION_INTRODUCTION>(operands);
	}

	template<template<typename> class Array>
	static inline Proof* new_conjunction_elim(Proof* proof, const Array<unsigned int>& indices) {
		if (!Canonical) return NULL;
		return new_binary_step<nd_step_type::CONJUNCTION_ELIMINATION>(proof, new_parameter(indices));
	}

	static inline Proof* new_conjunction_elim_left(Proof* proof) {
		if (Canonical) return NULL;
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_LEFT>(proof);
	}

	static inline Proof* new_conjunction_elim_right(Proof* proof) {
		if (Canonical) return NULL;
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_RIGHT>(proof);
	}

	static inline Proof* new_disjunction_intro(Proof* proof, Formula* parameter) {
		if (!Canonical) return NULL;
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION>(proof, new_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_left(Proof* proof, Formula* parameter) {
		if (Canonical) return NULL;
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_LEFT>(proof, new_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_right(Proof* proof, Formula* parameter) {
		if (Canonical) return NULL;
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT>(proof, new_parameter(parameter));
	}

	static inline Proof* new_disjunction_elim(Proof* disjunction, Proof* left, Proof* right) {
		return new_ternary_step<nd_step_type::DISJUNCTION_ELIMINATION>(disjunction, left, right);
	}

	static inline Proof* new_implication_intro(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_proof_type::AXIOM) return NULL;
		return new_binary_step<nd_step_type::IMPLICATION_INTRODUCTION>(proof, assumption);
	}

	static inline Proof* new_implication_elim(Proof* implication, Proof* antecedent) {
		return new_unary_step<nd_step_type::IMPLICATION_ELIMINATION>(implication, antecedent);
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

	static inline Proof* new_proof_by_contradiction(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_proof_type::AXIOM) return NULL;
		return new_binary_step<nd_step_type::PROOF_BY_CONTRADICTION>(proof, assumption);
	}

	static inline Proof* new_negation_elim(Proof* proof, Proof* negated) {
		return new_binary_step<nd_step_type::NEGATION_ELIMINATION>(proof, negated);
	}

	static inline Proof* new_universal_intro(Proof* proof, unsigned int parameter) {
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_parameter(parameter));
	}

	static inline Proof* new_universal_elim(Proof* proof, const Term& term) {
		return new_binary_step<nd_step_type::UNIVERSAL_ELIMINATION>(proof, new_parameter(term));
	}

	static inline Proof* new_existential_intro(Proof* proof, const unsigned int* term_indices, unsigned int term_index_count) {
		return new_binary_step<nd_step_type::EXISTENTIAL_INTRODUCTION>(proof, new_parameter(term_indices, term_index_count));
	}

	static inline Proof* new_existential_elim(Proof* existential, Proof* proof) {
		return new_binary_step<nd_step_type::EXISTENTIAL_ELIMINATION>(existential, proof);
	}

private:
	static inline Proof* new_parameter(Formula* parameter) {
		new_formula_parameterized_step<nd_step_type::FORMULA_PARAMETER>(parameter);
	}

	static inline Proof* new_parameter(unsigned int parameter) {
		nd_step<Formula, Canonical>* step;
		if (!new_nd_step(step, nd_step_type::PARAMETER)) return NULL;
		step->reference_count = 0;
		step->parameter = parameter;
		axiom->reference_count++;
		return step;
	}

	static inline Proof* new_parameter(const Term& term) {
		nd_step<Formula, Canonical>* step;
		if (!new_nd_step(step, nd_step_type::TERM_PARAMETER)) return NULL;
		step->reference_count = 0;
		step->term = term;
		axiom->reference_count++;
		return step;
	}

	template<template<typename> class Array>
	static inline Proof* new_parameter(const Array<unsigned int>& parameters)
	{
		nd_step<Formula, Canonical>* step;
		if (!new_nd_step(step, nd_step_type::ARRAY_PARAMETER)) {
			return NULL;
		} else if (!array_init(step->parameters, max((size_t) 1, parameters.length))) {
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < parameters.length; i++)
			step->parameters[i] = parameters[i];
		step->reference_count = 0;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_parameterized_step(Formula* parameter) {
		if (parameter == NULL) return NULL;

		nd_step<Formula, Canonical>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		step->formula = parameter;
		parameter->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_unary_step(Proof* proof) {
		if (proof == NULL) return NULL;

		nd_step<Formula, Canonical>* step;
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

		nd_step<Formula, Canonical>* step;
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

		nd_step<Formula, Canonical>* step;
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
	static inline bool new_array_step_helper(nd_step<Formula, Canonical>* step, Proof* arg) {
		step.operand_array[Index] = arg;
		arg->reference_count++;
		if (!arg->children.add(step)) {
			free(*arg); return false;
		}
		return true;
	}

	template<unsigned int Index, typename... Args>
	static inline bool new_array_step_helper(nd_step<Formula, Canonical>* step, Proof* arg, Args&&... args) {
		step.operand_array[Index] = arg;
		arg->reference_count++;
		if (!arg->children.add(step)) {
			free(*arg); return false;
		} else if (!new_array_step<Index + 1>(step, std::forward<Args>(args)...)) {
			arg->remove_child(step);
			free(*arg); return false;
		}
		return true;
	}

	template<nd_step_type Type, typename... Args>
	static inline Proof* new_array_step(Args&&... args) {
		static_assert(sizeof...(Args) >= 2, "This proof step requires at least two operands.");

		nd_step<Formula, Canonical>* step;
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

	template<nd_step_type Type, template<typename> class Array>
	static inline Proof* new_array_step(const Array<Proof>& operands) {
		if (operands.length < 2) return NULL;

		nd_step<Formula, Canonical>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 0;
		if (!array_init(step->operand_array, operands.length)) {
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < operands.length; i++) {
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
		step->operand_array.length = operand_count;
		return step;
	}
};

template<typename Formula, bool Canonical>
struct nd_canonicalizer {
	inline bool operator() (const nd_step<Formula, Canonical>* first, const nd_step<Formula, Canonical>* second) {
		return *first > *second;
	}
};

template<typename Formula, bool Canonical, typename... ProofMap>
bool canonicalize(const nd_step<Formula, Canonical>& proof,
		array<const nd_step<Formula, Canonical>*>& canonical_order,
		ProofMap&&... proof_map)
{
	hash_map<const nd_step<Formula, Canonical>*, unsigned int> in_degrees(128);
	if (!compute_in_degrees(&proof, in_degrees, std::forward<ProofMap>(proof_map)...)) return false;

	array<const nd_step<Formula, Canonical>*> heap(32);
	for (const auto& entry : in_degrees) {
		if (entry.value == 0 && !heap.add(entry.key))
			return false;
	}

	std::make_heap(heap.begin(), heap.end(), nd_canonicalizer<Formula>());
	while (heap.length > 0) {
		std::pop_heap(heap.begin(), heap.end());
		const nd_step<Formula, Canonical>* node = heap.pop();
		if (!canonical_order.add(node)) return false;

		unsigned int operand_count;
		const nd_step<Formula, Canonical>* const* operands;
		node->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			const nd_step<Formula, Canonical>* operand = map(
					operands[i], std::forward<ProofMap>(proof_map)...);
			unsigned int& degree = in_degrees.get(operand);
			degree--;

			if (degree == 0) {
				if (!heap.add(operand))
					return false;
				std::push_heap(heap.begin(), heap.end());
			}
		}
	}

	return true;
}

template<typename Formula,
	bool Canonical, typename FormulaPrior,
	typename ConjunctionIntroductionPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula, Canonical>& proof,
		unsigned int& formula_counter,
		array<typename Formula::Term>& introduced_terms,
		array<unsigned int>& available_parameters,
		FormulaPrior& formula_prior,
		ConjunctionIntroductionPrior& conjunction_introduction_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior,
		ProofMap&&... proof_map)
{
	typedef typename Formula::TermType TermType;

	double value; unsigned int index;
	switch (proof.type) {
	case TERM_PARAMETER:
	case ARRAY_PARAMETER:
	case FORMULA_PARAMETER:
		/* these aren't actual proof steps in the calculus */
		return 0.0;
	case AXIOM:
		/* TODO: we need to compute the prior */
		formula_counter++;
		get_parameters(*proof->formula, available_parameters);
		if (available_parameters.length > 1) {
			sort(available_parameters); unique(available_parameters);
		}
		return log_probability(*proof->formula, formula_prior);
	case CONJUNCTION_ELIMINATION:
		/* TODO: implement this */
		fprintf(stderr, "log_probability ERROR: Not implemented.\n");
		exit(EXIT_FAILURE);
	case CONJUNCTION_INTRODUCTION:
		return -LOG_ND_RULE_COUNT + log_probability(proof.operand_array.length, conjunction_introduction_prior)
			  + lgamma(formula_counter - proof.operand_array.length - 1) - lgamma(formula_counter++);
	case IMPLICATION_INTRODUCTION: /* TODO: is this correct? */
	case IMPLICATION_ELIMINATION:
	case BICONDITIONAL_INTRODUCTION:
	case BICONDITIONAL_ELIMINATION_LEFT:
	case BICONDITIONAL_ELIMINATION_RIGHT:
	case PROOF_BY_CONTRADICTION: /* TODO: is this correct? */
	case NEGATION_ELIMINATION:
	case EXISTENTIAL_ELIMINATION:
		return -LOG_ND_RULE_COUNT - 2*log_cache<V>::instance().get(formula_counter++);
	case DISJUNCTION_ELIMINATION:
		return -LOG_ND_RULE_COUNT - 3*log_cache<V>::instance().get(formula_counter++);
	case DISJUNCTION_INTRODUCTION:
	case DISJUNCTION_INTRODUCTION_LEFT:
	case DISJUNCTION_INTRODUCTION_RIGHT:
		/* TODO: we need to compute the prior on the new formula */
		formula_counter++;
		fprintf(stderr, "log_probability ERROR: Not implemented.\n"); exit(EXIT_FAILURE);
	case UNIVERSAL_INTRODUCTION:
		/* TODO: we need to compute the prior on the parameter */
		formula_counter++;
		operand = map(proof.operands[1], std::forward<ProofMap>(proof_map)...);
		value = log_probability(operand->parameter, universal_introduction_prior, available_parameters);
		index = available_parameters.index_of(operand->parameter);
		if (index < available_parameters.length)
			shift_left(available_parameters.data + index, available_parameters.length - index - 1);
		return value;
	case UNIVERSAL_ELIMINATION:
		/* TODO: we need to compute the prior on the term */
		formula_counter++;
		operand = map(proof.operands[2], std::forward<ProofMap>(proof_map)...);
		if (operand->term.type == TermType::PARAMETER) {
			available_parameters.add(operand->term.parameter);
			insertion_sort(available_parameters); unique(available_parameters);
		}
		return log_probability(operand->term, universal_elimination_prior);
	case EXISTENTIAL_INTRODUCTION:
		/* TODO: we need to compute the prior on the parameter (it can be a term or a list of term indices) */
		formula_counter++;
		fprintf(stderr, "log_probability ERROR: Not implemented.\n"); exit(EXIT_FAILURE);
	}
	fprintf(stderr, "log_probability ERROR: Unrecognized nd_step_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula,
	bool Canonical, typename FormulaPrior,
	typename UniversalIntroductionPrior,
	typename UniversalEliminationPrior,
	typename... ProofMap>
double log_probability(
		const nd_step<Formula, Canonical>& proof,
		double log_stop_probability,
		double log_continue_probability,
		FormulaPrior& formula_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior,
		ProofMap&&... proof_map)
{
	array<const nd_step<Formula, Canonical>*> canonical_order(64);
	if (!canonicalize(proof, canonical_order, std::forward<ProofMap>(proof_map)...)) {
		fprintf(stderr, "log_probability ERROR: Unable to canonicalize proof.\n");
		exit(EXIT_FAILURE);
	}

	universal_introduction_prior.clear();
	universal_elimination_prior.clear();

	double value = (canonical_order.length - 1) * log_continue_probability + log_stop_probability;
	unsigned int formula_counter = 0;
	array<typename Formula::Term>& introduced_terms(16);
	log_cache<V>::instance().ensure_size(canonical_order.length);
	for (const nd_step<Formula, Canonical>* step : canonical_order)
		value += log_probability(*step, counter, introduced_terms,
				log_max_parameter_count, formula_prior, universal_introduction_prior,
				universal_elimination_prior, std::forward<ProofMap>(proof_map)...);
	return value;
}

#endif /* NATURAL_DEDUCTION_H_ */
