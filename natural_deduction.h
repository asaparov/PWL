#ifndef NATURAL_DEDUCTION_H_
#define NATURAL_DEDUCTION_H_

#include <core/array.h>

using namespace core;

enum class nd_step_type
{
	AXIOM,
	ARRAY_PARAMETER,
	TERM_PARAMETER,
	FORMULA_PARAMETER,

	CONJUNCTION_INTRODUCTION,
	CONJUNCTION_ELIMINATION_LEFT,
	CONJUNCTION_ELIMINATION_RIGHT,
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
	EXISTENTIAL_ELIMINATION
};

constexpr static unsigned int ND_OPERAND_COUNT = 3;

template<typename Formula>
struct nd_step
{
	typedef typename Formula::Term Term;

	nd_step_type type;
	unsigned int reference_count;
	union {
		Term term;
		Formula* formula;
		array<unsigned int> parameters;

		nd_step<Formula>* operands[ND_OPERAND_COUNT];
	};

	static inline void free(nd_step<Formula>& step) { step.free(); }

	inline bool has_subproofs() {
		switch (type) {
		case TERM_PARAMETER:
		case ARRAY_PARAMETER:
		case AXIOM:
		case FORMULA_PARAMETER:
			return false;
		case CONJUNCTION_INTRODUCTION:
		case CONJUNCTION_ELIMINATION_LEFT:
		case CONJUNCTION_ELIMINATION_RIGHT:
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
			return false;
		}
		fprintf(stderr, "nd_step.has_subproofs ERROR: Unrecognized nd_step_type.\n");
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
		case CONJUNCTION_ELIMINATION_LEFT:
		case CONJUNCTION_ELIMINATION_RIGHT:
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
				core::free(*operands[i]);
				if (operands[i]->reference_count == 0)
					core::free(operands[i]);
			}
			return;
		}
		fprintf(stderr, "nd_step.free ERROR: Unrecognized nd_step_type.\n");
		exit(EXIT_FAILURE);
	}
};

template<typename Formula>
Formula* compute_conclusion(const nd_step<Formula>& proof)
{
	typedef typename Formula::Type FormulaType;

	Formula* formula; Formula* operand;
	switch (proof.type) {
	case nd_step_type::AXIOM:
		proof.formula->reference_count++;
		return proof.formula;
	case nd_step_type::CONJUNCTION_INTRODUCTION:
		if (proof.operands[0] == NULL || proof.operands[1] == NULL) return NULL;
		return Formula::new_and(compute_conclusion(*proof.operands[0], *proof.operands[1]));
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
		if (proof.operands[0] == NULL) return NULL;
		formula = compute_conclusion(*proof.operands[0]);
		if (formula == NULL) {
			return NULL;
		} else if (formula->type != FormulaType) {
			free(*formula);
			if (formula->reference_count == 0)
				free(formula);
			return NULL;
		}
		operand = formula->binary.left;
		operand->reference_count++;
		free(*formula);
		if (formula->reference_count == 0)
			free(formula);
		return operand;
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
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

	case nd_step_type::UNIVERSAL_INTRODUCTION:
	case nd_step_type::UNIVERSAL_ELIMINATION:
	case nd_step_type::EXISTENTIAL_INTRODUCTION:
	case nd_step_type::EXISTENTIAL_ELIMINATION:

	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
	}
}

template<typename Formula>
struct proof_checker {
	const nd_step<Formula>* step;
	Formula* operands[ND_OPERAND_COUNT];
	array<Formula*> assumptions;
};

template<typename Formula>
Formula* compute_conclusion(const nd_step<Formula>& proof)
{
	/* first list the proof steps in reverse topological order */
	array<const nd_step<Formula>*> first_stack(64);
	

	array<nd_step<Formula>*> first_stack(32);
	array<proof_checker> second_stack(32);
	if (!first_stack.add(&proof)) {
		fprintf(stderr, "compute_conclusion ERROR: Unable to add root node to first stack.\n");
		exit(EXIT_FAILURE);
	}

	while (first_stack.length > 0) {
		nd_step<Formula>* node = first_stack.last();
		if (!second_stack.add(node)) {
			fprintf(stderr, "compute_conclusion ERROR: Unable to add node to second stack.\n");
			exit(EXIT_FAILURE);
		}

		if (!node->has_subproofs()) continue;
		for (unsigned int i = 0; i < ND_OPERAND_COUNT; i++) {
			if (node->operands[i] == NULL)
				continue;
			if (!first_stack.add(node->operands[i])) { 
				fprintf(stderr, "compute_conclusion ERROR: Unable to add node to first stack.\n");
				exit(EXIT_FAILURE);
			}
		}
	}
}

template<typename Formula>
bool check_proof(
		const nd_step<Formula>& proof,
		const Formula* expected_conclusion)
{
	array<const nd_step<formula>*> stack(64);
	const nd_step<Formula>* node = &proof;
	while (stack.length > 0 || node != NULL) {
		if (node != NULL) {
			if (!stack.add()) {
				fprintf(stderr, "check_proof ERROR: Unable to add to stack.\n");
				exit(EXIT_FAILURE);
			}
		}


	}
}

template<typename Formula>
bool new_nd_step(nd_step<Formula>*& step, nd_step_type type) {
	step = (nd_step<Formula>*) malloc(sizeof(nd_step<Formula>));
	if (step == NULL) {
		fprintf(stderr, "new_nd_step ERROR: Out of memory.\n");
		return false;
	}
	step->type = type;
	return true;
}

template<typename Formula>
struct natural_deduction
{
	typedef nd_step<Formula> Proof;
	typedef typename Formula::Term Term;

	static inline Proof* new_axiom(Formula* axiom) {
		new_parameterized_step<nd_step_type::AXIOM>(axiom);
	}

	static inline Proof* new_conjunction_intro(Proof* left, Proof* right) {
		return new_binary_step<nd_step_type::CONJUNCTION_INTRODUCTION>(left, right);
	}

	static inline Proof* new_conjunction_elim_left(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_LEFT>(proof);
	}

	static inline Proof* new_conjunction_elim_right(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_RIGHT>(proof);
	}

	static inline Proof* new_disjunction_intro_left(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_LEFT>(proof, new_parameter(parameter));
	}

	static inline Proof* new_disjunction_intro_right(Proof* proof, Formula* parameter) {
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

	static inline Proof* new_universal_intro(Proof* proof, const Term& term) {
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_parameter(term));
	}

	static inline Proof* new_universal_intro(Proof* proof, const array<unsigned int>& variable_indices) {
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_parameter(variable_indices));
	}

	static inline Proof* new_universal_elim(Proof* proof, const Term& term) {
		return new_binary_step<nd_step_type::UNIVERSAL_ELIMINATION>(proof, new_parameter(term));
	}

	static inline Proof* new_existential_intro(Proof* proof, const Term& term) {
		return new_binary_step<nd_step_type::EXISTENTIAL_INTRODUCTION>(proof, new_parameter(term));
	}

	static inline Proof* new_existential_intro(Proof* proof, const array<unsigned int>& term_indices) {
		return new_binary_step<nd_step_type::EXISTENTIAL_INTRODUCTION>(proof, new_parameter(term_indices));
	}

	static inline Proof* new_existential_elim(Proof* existential, Proof* proof) {
		return new_binary_step<nd_step_type::EXISTENTIAL_ELIMINATION>(existential, proof);
	}

private:
	static inline Proof* new_parameter(Formula* parameter) {
		new_formula_parameterized_step<nd_step_type::FORMULA_PARAMETER>(parameter);
	}

	static inline Proof* new_parameter(const Term& term) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::TERM_PARAMETER)) return NULL;
		step->reference_count = 1;
		step->term = term;
		axiom->reference_count++;
		return step;
	}

	static inline Proof* new_parameter(const array<unsigned int>& parameters) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::ARRAY_PARAMETER)) {
			return NULL;
		} else if (!array_init(step->parameters, max((size_t) 1, parameters.length))) {
			free(step); return NULL;
		}
		for (unsigned int i = 0; i < parameters.length; i++)
			step->parameters[i] = parameters[i];
		step->reference_count = 1;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_parameterized_step(Formula* parameter) {
		if (parameter == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 1;
		step->formula = axiom;
		axiom->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_unary_step(Proof* proof) {
		if (proof == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 1;
		step->operands[0] = proof;
		step->operands[1] = NULL;
		step->operands[2] = NULL;
		proof->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_binary_step(Proof* left, Proof* right) {
		if (left == NULL || right == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 1;
		step->operands[0] = left;
		step->operands[1] = right;
		step->operands[2] = NULL;
		left->reference_count++;
		right->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	static inline Proof* new_ternary_step(Proof* first, Proof* second, Proof* third) {
		if (first == NULL || second == NULL || third == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 1;
		step->operands[0] = first;
		step->operands[1] = second;
		step->operands[2] = third;
		first->reference_count++;
		second->reference_count++;
		third->reference_count++;
		return step;
	}
};

#endif /* NATURAL_DEDUCTION_H_ */
