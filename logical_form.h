/**
 * logical_form.h
 *
 *  Created on: Jun 22, 2018
 *      Author: asaparov
 */

#ifndef LOGICAL_FORM_H_
#define LOGICAL_FORM_H_

#include <core/map.h>
#include <core/utility.h>

using namespace core;

/* TODO: make this better (so we don't have to spend an hour renumbering
		 everything when we add/remove predicates) */
constexpr unsigned int PREDICATE_TYPE = 1;
constexpr unsigned int PREDICATE_AND = 2;
constexpr unsigned int PREDICATE_OR = 3;

bool populate_name_map(hash_map<string, unsigned int>& names) {
	bool success = true;
	success &= names.put("type", PREDICATE_TYPE);
	success &= names.put("and", PREDICATE_AND);
	success &= names.put("or", PREDICATE_OR);
	return success;
}

struct quantification {
	unsigned int group_count;
	unsigned int* group_sizes;
	unsigned int* variables;
};

enum class lf_expression_type {
	CONSTANT,
	VERTEX
};

struct lf_literal {
	unsigned int label;
};

struct lf_vertex {
	unsigned int variable;
	lf_expression* type;
	lf_expression* arg1;
	lf_expression* arg2;

	unsigned int raised;
	quantification indices;
};

struct lf_expression {
	lf_expression_type type;
	unsigned int reference_count;
	union {
		lf_literal constant;
		lf_vertex vertex;
	};
};

#endif /* LOGICAL_FORM_H_ */
