#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>

#include "first_order_logic.h"

using namespace core;


enum built_in_predicates : unsigned int {
	PREDICATE_TYPE = 1,
	PREDICATE_UNKNOWN,
	PREDICATE_ARG1,
	PREDICATE_ARG2
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

	bool add_formula(fol_formula* formula) {
		/* TODO: fill this in */
	}
};

#endif /* THEORY_H_ */
