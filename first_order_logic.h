/**
 * first_order_logic.h
 *
 *  Created on: Jun 26, 2018
 *      Author: asaparov
 */

#ifndef FIRST_ORDER_LOGIC_H_
#define FIRST_ORDER_LOGIC_H_

enum class fol_term_type {
	VARIABLE,
	CONSTANT
};

struct fol_term {
	fol_term_type type;
	union {
		unsigned int variable;
		unsigned int constant;
	};
};

enum class fol_expression_type {
	ATOM,

	AND,
	OR,
	IF_THEN,
	IFF,
	NOT,

	FOR_ALL,
	EXISTS
};

struct fol_atom {
	unsigned int predicate;
	fol_term arg1;
	fol_term arg2;
};

struct fol_and_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_or_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_if_then_expression {
	fol_expression* antecedent;
	fol_expression* consequent;
};

struct fol_iff_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_not_expression {
	fol_expression* operand;
};

struct fol_for_all_expression {
	unsigned int variable;
	fol_expression* operand;
};

struct fol_exists_expression {
	unsigned int variable;
	fol_expression* operand;
};

struct fol_expression {
	fol_expression_type type;
	unsigned int reference_count;
	union {
		fol_atom atom;
		fol_and_expression conj;
		fol_or_expression disj;
		fol_if_then_expression if_then;
		fol_iff_expression iff;
		fol_for_all_expression for_all;
		fol_exists_expression exists;
	};
};

#endif /* FIRST_ORDER_LOGIC_H_ */
