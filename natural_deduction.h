#ifndef NATURAL_DEDUCTION_H_
#define NATURAL_DEDUCTION_H_

enum class nd_step_type
{
	AXIOM,
	PARAMETER,
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
		unsigned int parameter;
		nd_step<Formula>* operands[ND_OPERAND_COUNT];
	};

	static inline void free(nd_step<Formula>& step) { step.free(); }

private:
	inline void free() {
		switch (type) {
		case PARAMETER:
			return;
		case TERM_PARAMETER:
			free(term); return;
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

	inline Proof* new_axiom(Formula* axiom) {
		new_parameterized_step<nd_step_type::AXIOM>(axiom);
	}

	inline Proof* new_conjunction_intro(Proof* left, Proof* right) {
		return new_binary_step<nd_step_type::CONJUNCTION_INTRODUCTION>(left, right);
	}

	inline Proof* new_conjunction_elim_left(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_LEFT>(proof);
	}

	inline Proof* new_conjunction_elim_right(Proof* proof) {
		return new_unary_step<nd_step_type::CONJUNCTION_ELIMINATION_RIGHT>(proof);
	}

	inline Proof* new_disjunction_intro_left(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_LEFT>(proof, new_parameter(parameter));
	}

	inline Proof* new_disjunction_intro_right(Proof* proof, Formula* parameter) {
		return new_binary_step<nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT>(proof, new_parameter(parameter));
	}

	inline Proof* new_disjunction_elim(Proof* disjunction, Proof* left, Proof* right) {
		return new_ternary_step<nd_step_type::DISJUNCTION_ELIMINATION>(disjunction, left, right);
	}

	inline Proof* new_implication_intro(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_proof_type::AXIOM) return NULL;
		return new_binary_step<nd_step_type::IMPLICATION_INTRODUCTION>(proof, assumption);
	}

	inline Proof* new_implication_elim(Proof* implication, Proof* antecedent) {
		return new_unary_step<nd_step_type::IMPLICATION_ELIMINATION>(implication, antecedent);
	}

	inline Proof* new_biconditional_intro(Proof* forward, Proof* backward) {
		return new_binary_step<nd_step_type::BICONDITIONAL_INTRODUCTION>(forward, backward);
	}

	inline Proof* new_biconditional_elim_left(Proof* biconditional, Proof* left) {
		return new_binary_step<nd_step_type::BICONDITIONAL_ELIMINATION_LEFT>(biconditional, left);
	}

	inline Proof* new_biconditional_elim_right(Proof* biconditional, Proof* right) {
		return new_binary_step<nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT>(biconditional, right);
	}

	inline Proof* new_proof_by_contradiction(Proof* proof, Proof* assumption) {
		if (assumption == NULL || assumption->type != nd_proof_type::AXIOM) return NULL;
		return new_binary_step<nd_step_type::PROOF_BY_CONTRADICTION>(proof, assumption);
	}

	inline Proof* new_negation_elim(Proof* proof, Proof* negated) {
		return new_binary_step<nd_step_type::NEGATION_ELIMINATION>(proof, negated);
	}

	inline Proof* new_universal_intro(Proof* proof, unsigned int parameter) {
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_parameter(negated));
	}

	inline Proof* new_universal_elim(Proof* proof, const Term& term) {
		return new_binary_step<nd_step_type::UNIVERSAL_ELIMINATION>(proof, new_parameter(term));
	}

	inline Proof* new_existential_intro(Proof* proof, unsigned int parameter) {
		/* TODO: change the second argument to be a list of term indices which indicate which terms will be replaced by a new variable */
		return new_binary_step<nd_step_type::UNIVERSAL_INTRODUCTION>(proof, new_parameter(negated));
	}

private:
	inline Proof* new_parameter(Formula* parameter) {
		new_formula_parameterized_step<nd_step_type::FORMULA_PARAMETER>(parameter);
	}

	inline Proof* new_parameter(unsigned int parameter) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::PARAMETER)) return NULL;
		step->reference_count = 1;
		step->parameter = parameter;
		axiom->reference_count++;
		return step;
	}

	inline Proof* new_parameter(const Term& term) {
		nd_step<Formula>* step;
		if (!new_nd_step(step, nd_step_type::TERM_PARAMETER)) return NULL;
		step->reference_count = 1;
		step->term = term;
		axiom->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	inline Proof* new_parameterized_step(Formula* parameter) {
		if (parameter == NULL) return NULL;

		nd_step<Formula>* step;
		if (!new_nd_step(step, Type)) return NULL;
		step->reference_count = 1;
		step->formula = axiom;
		axiom->reference_count++;
		return step;
	}

	template<nd_step_type Type>
	inline Proof* new_unary_step(Proof* proof) {
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
	inline Proof* new_binary_step(Proof* left, Proof* right) {
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
	inline Proof* new_ternary_step(Proof* first, Proof* second, Proof* third) {
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
