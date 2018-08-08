#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>

#include "first_order_logic.h"

using namespace core;

struct theory {
	array<fol_formula*> definitions;
	array<fol_formula*> formulas;

	bool add_formula(fol_formula* formula) {
		/* TODO: fill this in */
	}
};

#endif /* THEORY_H_ */
