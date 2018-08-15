#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>

#include "first_order_logic.h"

using namespace core;


enum built_in_predicates : unsigned int {
	PREDICATE_TYPE = 1,
	PREDICATE_UNKNOWN,
	PREDICATE_ARG1,
	PREDICATE_ARG2,
	PREDICATE_COUNT
};

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("type", PREDICATE_TYPE)
		&& names.put("unknown", PREDICATE_UNKNOWN)
		&& names.put("arg1", PREDICATE_ARG1)
		&& names.put("arg2", PREDICATE_ARG2);
}

struct theory {
	array<fol_formula*> definitions;
	array<fol_formula*> formulas;

	bool add_formula(fol_formula* formula)
	{
		fol_formula* canonicalized = canonicalize(*formula);
		if (!add_canonicalized_formula(canonicalized)) {
			free(*canonicalized);
			if (canonicalized->reference_count == 0)
				free(canonicalized);
			return false;
		}
		return true;
	}

private:
	bool add_canonicalized_formula(fol_formula* canonicalized)
	{
		if (canonicalized->type == fol_formula_type::ATOM) {
			if (canonicalized->atom.predicate == PREDICATE_TYPE
			 && canonicalized->atom.arg1.type == fol_term_type::CONSTANT)
			{
				if (canonicalized->atom.arg1.constant == PREDICATE_UNKNOWN) {
					/* this is a definition of an object */
					canonicalized->atom.arg1.constant = PREDICATE_COUNT + definitions.length;
					if (!definitions.add(canonicalized)) return false;
				} else {
					/* this is a formula of form `type(c,t)` */
					return add_formula_helper(canonicalized);
				}
			}
		} else if (canonicalized->type == fol_formula_type::FOR_ALL) {
			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == fol_formula_type::IF_THEN) {
				const fol_formula* left = canonicalized->quantifier.operand->binary.left;
				if (left->type == fol_formula_type::ATOM
				 && left->atom.predicate == PREDICATE_TYPE
				 && left->atom.arg1.type == fol_term_type::VARIABLE
				 && left->atom.arg1.variable == variable
				 && left->atom.arg2.type == fol_term_type::CONSTANT)
				{
					if (left->atom.arg2.constant == PREDICATE_UNKNOWN) {
						/* this is a definition of a type */
						fol_formula* right = canonicalized->quantifier.operand.binary.right;

						/* check the right-hand side is a valid definition */
						if (!valid_definition(right, variable)) {
							fprintf(stderr, "add_canonicalized_formula ERROR: This is not a valid type definition.\n");
							return false;
						}

						right->reference_count++;
						unsigned int new_type = PREDICATE_COUNT + definitions.length;
						fol_formula* definition = make_fol_for_all(variable, make_fol_iff(make_fol_atom(PREDICATE_TYPE, make_fol_variable(variable), make_fol_constant(new_type)), right));
						if (definition == NULL || !definitions.add(definition))
							return false;
					} else {
						/* this is a formula of form `![x]:(type(x,t) => f(x))` */
						return add_formula_helper(canonicalized);
					}
				}
			}
		} else if (canonicalized->type == fol_formula_type::EXISTS) {
			return add_formula_helper(canonicalized);
		} else if (canonicalized->type == fol_formula_type::AND) {
			return add_canonicalized_formula(canonicalized->binary.left)
				&& add_canonicalized_formula(canonicalized->binary.right);
		} else {
			fprintf(stderr, "theory.add_canonicalized_formula ERROR: Unsupported fol_formula_type.\n");
			return false;
		}
		return true;
	}

	inline bool add_formula_helper(fol_formula* formula) {
		for (const fol_formula* existing : formulas) {
			if (*existing == *formula) return true;
		}
		if (!formulas.add(formula)) return false;
		formula->reference_count++;
		return true;
	}

	inline bool valid_definition(const fol_formula* right,
			unsigned int quantified_variable)
	{
		if (right->type == fol_formula_type::ATOM) {
			return right->atom.predicate == PREDICATE_TYPE
				&& right->atom.arg1.type == fol_term_type::VARIABLE
				&& right->atom.arg1.variable == quantified_variable
				&& right->atom.arg2.type == fol_term_type::CONSTANT);
		} else if (right->type == fol_formula_type::AND) {
			return valid_definition(right->binary.left, quantified_variable)
				|| valid_definition(right->binary.right, quantified_variable);
		} else {
			return false;
		}
	}
};

#endif /* THEORY_H_ */
