/**
 * first_order_logic.h
 *
 *  Created on: Jun 26, 2018
 *      Author: asaparov
 */

#ifndef FIRST_ORDER_LOGIC_H_
#define FIRST_ORDER_LOGIC_H_

#include <core/lex.h>
#include <cstdint>

using namespace core;


/* forward declarations */

struct fol_formula;

enum class fol_term_type {
	NONE,
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

template<typename Stream>
inline bool print_variable(unsigned int variable, Stream& out) {
	return print('$', out) + print(variable, out);
}

template<typename Stream, typename... Printer>
bool print(const fol_term& term, Stream& out, Printer&&... printer) {
	switch (term.type) {
	case fol_term_type::VARIABLE:
		return print_variable(term.variable, out);
	case fol_term_type::CONSTANT:
		return print(term.constant, out, std::forward<Printer>(printer)...);
	case fol_term_type::NONE:
		break;
	}

	fprintf(stderr, "print ERROR: Unexpected fol_term_type.\n");
	return false;
}

enum class fol_formula_type {
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

	static inline void move(const fol_atom& src, fol_atom& dst) {
		dst.predicate = src.predicate;
		dst.arg1 = src.arg1;
		dst.arg2 = src.arg2;
	}

	static inline void free(fol_atom& formula) { }
};

struct fol_unary_formula {
	fol_formula* operand;

	static inline void move(const fol_unary_formula& src, fol_unary_formula& dst);
	static inline void free(fol_unary_formula& formula);
};

struct fol_binary_formula {
	fol_formula* left;
	fol_formula* right;

	static inline void move(const fol_binary_formula& src, fol_binary_formula& dst);
	static inline void free(fol_binary_formula& formula);
};

struct fol_quantifier {
	unsigned int variable;
	fol_formula* operand;

	static inline void move(const fol_quantifier& src, fol_quantifier& dst);
	static inline void free(fol_quantifier& formula);
};

struct fol_formula {
	fol_formula_type type;
	unsigned int reference_count;
	union {
		fol_atom atom;
		fol_unary_formula unary;
		fol_binary_formula binary;
		fol_quantifier quantifier;
	};

	static inline void move(const fol_formula& src, fol_formula& dst);
	static inline void free(fol_formula& formula);
};

inline void fol_formula::move(const fol_formula& src, fol_formula& dst) {
	dst.type = src.type;
	dst.reference_count = src.reference_count;
	switch (src.type) {
	case fol_formula_type::ATOM:
		core::move(src.atom, dst.atom); return;
	case fol_formula_type::NOT:
		core::move(src.unary, dst.unary); return;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
		core::move(src.binary, dst.binary); return;
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		core::move(src.quantifier, dst.quantifier); return;
	}
	fprintf(stderr, "fol_formula.move ERROR: Unrecognized fol_formula_type.\n");
}

inline void fol_unary_formula::move(const fol_unary_formula& src, fol_unary_formula& dst) {
	dst.operand = src.operand;
}

inline void fol_binary_formula::move(const fol_binary_formula& src, fol_binary_formula& dst) {
	dst.left = src.left;
	dst.right = src.right;
}

inline void fol_quantifier::move(const fol_quantifier& src, fol_quantifier& dst) {
	dst.variable = src.variable;
	dst.operand = src.operand;
}

inline void fol_formula::free(fol_formula& formula) {
	formula.reference_count--;
	if (formula.reference_count == 0) {
		switch (formula.type) {
		case fol_formula_type::ATOM:
			core::free(formula.atom); return;
		case fol_formula_type::NOT:
			core::free(formula.unary); return;
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IF_THEN:
		case fol_formula_type::IFF:
			core::free(formula.binary); return;
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			core::free(formula.quantifier); return;
		}
		fprintf(stderr, "fol_formula.free ERROR: Unrecognized fol_formula_type.\n");
	}
}

inline void fol_unary_formula::free(fol_unary_formula& formula) {
	core::free(*formula.operand);
	if (formula.operand->reference_count == 0)
		core::free(formula.operand);
}

inline void fol_binary_formula::free(fol_binary_formula& formula) {
	core::free(*formula.left);
	if (formula.left->reference_count == 0)
		core::free(formula.left);

	core::free(*formula.right);
	if (formula.right->reference_count == 0)
		core::free(formula.right);
}

inline void fol_quantifier::free(fol_quantifier& formula) {
	core::free(*formula.operand);
	if (formula.operand->reference_count == 0)
		core::free(formula.operand);
}

struct canonlicalizer { };

int_fast8_t compare(
		const fol_formula&,
		const fol_formula&,
		const canonlicalizer&);

inline int_fast8_t compare(
		const fol_term& first,
		const fol_term& second,
		const canonlicalizer& sorter)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;
	switch (first.type) {
	case fol_term_type::CONSTANT:
		if (first.constant < second.constant) return -1;
		else if (first.constant > second.constant) return 1;
		else return 0;
	case fol_term_type::VARIABLE:
		/* TODO: how do we compare variables? */
		return -1;
	case fol_term_type::NONE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline int_fast8_t compare(
		const fol_atom& first,
		const fol_atom& second,
		const canonlicalizer& sorter)
{
	if (first.predicate < second.predicate) return -1;
	else if (first.predicate > second.predicate) return 1;

	int_fast8_t result = compare(first.arg1, second.arg1, sorter);
	if (result != 0) return result;

	return compare(first.arg2, second.arg2, sorter);
}

inline int_fast8_t compare(
		const fol_unary_formula& first,
		const fol_unary_formula& second,
		const canonlicalizer& sorter)
{
	return compare(*first.operand, *second.operand, sorter);
}

inline int_fast8_t compare(
		const fol_binary_formula& first,
		const fol_binary_formula& second,
		const canonlicalizer& sorter)
{
	int_fast8_t result = compare(*first.left, *second.left, sorter);
	if (result != 0) return result;
	return compare(*first.right, *second.right, sorter);
}

inline int_fast8_t compare(
		const fol_quantifier& first,
		const fol_quantifier& second,
		const canonlicalizer& sorter)
{
	/* TODO: how do we compare variables? */
	return compare(*first.operand, *second.operand, sorter);
}

int_fast8_t compare(
		const fol_formula& first,
		const fol_formula& second,
		const canonlicalizer& sorter)
{
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return compare(first.atom, second.atom, sorter);
	case fol_formula_type::NOT:
		return compare(first.unary, second.unary, sorter);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
		return compare(first.binary, second.binary, sorter);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return compare(first.quantifier, second.quantifier, sorter);
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

bool less_than(
		const fol_formula& first,
		const fol_formula& second,
		const canonlicalizer& sorter)
{
	return compare(first, second, sorter) < 0;
}

struct fol_commutative_scope {
	array<const fol_formula*> children;

	fol_commutative_scope() : children(8) { }
};

struct fol_noncommutative_scope {
	array<const fol_formula*> left;
	array<const fol_formula*> right;
	bool is_left;

	fol_noncommutative_scope() : left(8), right(8), is_left(true) { }
};

struct fol_quantifier_scope {
	const fol_formula* operand;
	unsigned int variable;
};

struct fol_scope {
	fol_scope* parent;
	fol_formula_type type;
	union {
		const fol_formula* unary;
		fol_commutative_scope commutative;
		fol_noncommutative_scope noncommutative;
		fol_quantifier_scope quantifier;
	};

	fol_scope(fol_scope* parent, fol_formula_type type) : parent(parent), type(type) {
		switch (type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			new (&commutative) fol_commutative_scope();
		case fol_formula_type::IF_THEN:
			new (&noncommutative) fol_noncommutative_scope();
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			new (&quantifier) fol_quantifier_scope();
		case fol_formula_type::NOT:
		case fol_formula_type::ATOM:
			unary = NULL;
		}
		fprintf(stderr, "fol_scope ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}

	~fol_scope() {
		switch (type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			commutative.~fol_commutative_scope();
		case fol_formula_type::IF_THEN:
			noncommutative.~fol_noncommutative_scope();
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			quantifier.~fol_quantifier_scope();
		case fol_formula_type::NOT:
		case fol_formula_type::ATOM:
			/* TODO: can unary be null here? */
			free(*unary);
			if (unary->reference_count == 0)
				free(unary);
		}
		fprintf(stderr, "~fol_scope ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}
};

template<fol_formula_type ScopeType>
inline fol_formula* canonicalize_scope(array<fol_formula*>& scope)
{
	sort(scope, canonicalizer());

	/* remove duplicate elements */
	unsigned int dst_index = 0;
	for (unsigned int i = 1; i < scope.length; i++) {
		if (scope[dst_index] != scope[i]) {
			scope[++dst_index] = scope[i];
		} else {
			free(*scope[i]);
			free(scope[i]);
		}
	}
	scope.length = dst_index + 1;

	fol_formula* inner = scope.last();
	inner->reference_count++;
	for (unsigned int i = scope.length - 1; i > 0; i--) {
		fol_formula* new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "canonicalize_scope: Out of memory.\n");
			free(*inner); if (inner->reference_count == 0) free(inner);
			return NULL;
		}
		new_formula->type = ScopeType;
		new_formula->binary.left = scope[i - 1];
		new_formula->binary.right = inner;
		new_formula->binary.left->reference_count++;
		/* the reference count of 'inner' is incremented and decremented
		   simultaneously here, so we don't change it */
		inner = new_formula;
	}

	return inner;
}

inline bool move_antecedent_to_highest_scope(
		fol_scope& scope, const fol_formula* formula,
		const array<unsigned int>& variables)
{
	fol_scope* prev = &scope;
	for (fol_scope* curr = scope.parent; curr != NULL; prev = curr, curr = curr->parent) {
		if (curr->type == fol_formula_type::IF_THEN && curr->noncommutative.is_left) continue;
		prev = curr;
	}
}

inline bool move_and_to_highest_scope(
		fol_scope& scope, const fol_formula* formula,
		const array<unsigned int>& variables)
{
	fol_scope* prev = &scope;
	fol_scope* curr = scope.parent;
	for (; curr != NULL; prev = curr, curr = curr->parent) {
		bool highest_scope = false;
		switch (curr->type) {
		case fol_formula_type::AND: continue;
		case fol_formula_type::IF_THEN:
			if (curr->noncommutative.is_left)
				return move_antecedent_to_highest_scope(*curr, formula, variables);
			else return curr->noncommutative.right.add(formula);
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:

		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			return prev->commutative.children.add(formula);
		case fol_formula_type::NOT:
		case fol_formula_type::ATOM:
		}
	}
}

bool move_to_highest_scope(
		fol_scope& scope, const fol_formula* formula,
		const array<unsigned int>& variables)
{
	switch (scope.type) {
	case fol_formula_type::AND:
		return move_and_to_highest_scope(scope, formula, variables);
	case fol_formula_type::IF_THEN:
		if (scope.noncommutative.is_left)
			return move_antecendent_to_highest_scope(scope, formula, variables);
		else return move_consequent_to_highest_scope(scope, formula, variables);
	}
}

template<fol_formula_type ScopeType>
inline bool new_commutative_scope(
		const fol_formula& src, fol_scope& parent_scope,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope new_scope(&parent_scope, ScopeType);
	if (!canonicalize(*src.binary.left, new_scope, variable_map)
	 || !canonicalize(*src.binary.right, new_scope, variable_map))
		return false;

	/* construct the canonicalized node */
	if (new_scope.commutative.children.length > 1) {
		fol_formula* canonicalized = canonicalize_scope<ScopeType>(new_scope.commutative.children);
		if (canonicalized == NULL) return false;

		/* TODO: move 'canonicalized' to the appropriate ancestor scope */

	} else {
		/* TODO: what do we do if there is less than 2 operands? is this even possible? */
	}
	return true;
}

template<fol_formula_type ScopeType>
inline bool new_noncommutative_scope(
		const fol_formula& src, fol_scope& parent_scope,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope new_scope(&parent_scope, ScopeType);
	if (!canonicalize(*src.binary.left, new_scope, variable_map)) return false;
	new_scope.noncommutative.is_left = false;
	if (!canonicalize(*src.binary.right, new_scope, variable_map)) return false;

	/* construct the canonicalized node */
	fol_formula* canonicalized_left;
	if (new_scope.noncommutative.left.length > 1) {
		canonicalized_left = canonicalize_scope<fol_formula_type::AND>(new_scope.noncommutative.left);
		if (canonicalized_left == NULL) return false;
	} else if (new_scope.noncommutative.left.length == 1) {
		canonicalized_left = new_scope.noncommutative.left[0];
		canonicalized_left->reference_count++;
	} else {
		/* TODO: can new_scope.noncommutative.left be empty? */
	}

	fol_formula* canonicalized_right;
	if (new_scope.noncommutative.right.length > 1) {
		canonicalized_right = canonicalize_scope<fol_formula_type::OR>(new_scope.noncommutative.right);
		if (canonicalized_right == NULL) {
			free(*canonicalized_left);
			if (canonicalized_left->reference_count == 0)
				free(canonicalized_left);
			return false;
		}
	} else if (new_scope.noncommutative.right.length == 1) {
		canonicalized_right = new_scope.noncommutative.right[0];
		canonicalized_right->reference_count++;
	} else {
		/* TODO: can new_scope.noncommutative.right be empty? */
	}

	fol_formula* canonicalized = (fol_formula*) malloc(sizeof(fol_formula));
	if (canonicalized == NULL) {
		fprintf(stderr, "new_noncommutative_scope ERROR: Out of memory.\n");
		free(*canonicalized_left); free(*canonicalized_right);
		if (canonicalized_left->reference_count == 0) free(canonicalized_left);
		if (canonicalized_right->reference_count == 0) free(canonicalized_right);
		return false;
	}
	canonicalized->type = ScopeType;
	canonicalized->binary.left = canonicalized_left;
	canonicalized->binary.right = canonicalized_right;

	/* TODO: move 'canonicalized' to the appropriate ancestor scope */

	return true;
}

inline bool canonicalize_negation(
		const fol_formula& src, fol_scope& parent_scope,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope new_scope(&parent_scope, fol_formula_type::NOT);
	/* TODO: do we need to initialize new_scope.unary? */
	if (!canonicalize(*src.unary.operand, new_scope, variable_map))
		return false;

	/* negation blocks all movement, so new_scope.unary->type cannot be NOT */
	/* TODO: move 'canonicalized' to the appropriate ancestor scope */
	const fol_formula* canonicalized = new_scope.unary;

	return true;
}

template<fol_formula_type QuantifierType>
inline bool new_quantifier_scope(
		const fol_formula& src, fol_scope& parent_scope,
		array_map<unsigned int, unsigned int>& variable_map,
		array<unsigned int>& variables)
{
	if (!variable_map.ensure_capacity(variable_map.size + 1))
		return false;
	unsigned int index = variable_map.index_of(src.quantifier.variable);
	if (index < variable_map.size) {
		fprintf(stderr, "new_quantifier_scope ERROR: Multiple declaration of variable %u.\n", src.quantifier.variable);
		return false;
	}
	variable_map.keys[index] = src.quantifier.variable;
	variable_map.values[index] = variable_map.size + 1;
	variable_map.size++;

	fol_scope new_scope(parent_scope, QuantifierType);
	new_scope.quantifier.variable = variable_map.values[index];
	array<unsigned int> child_variables = array<unsigned int>(8);
	if (!canonicalize(*src.quantifier.operand, new_scope, variable_map, child_variables))
		return false;

	variable_map.size = index;
	unsigned int var_index = child_variables.index_of(new_scope.quantifier.variable);
	if (var_index < child_variables.length) {
		shift_left(child_variables.data + var_index, child_variables.length - var_index - 1);
		child_variables.length--;
	}

	/* construct the canonicalized node */
	fol_formula* canonicalized = (fol_formula*) malloc(sizeof(fol_formula));
	if (canonicalized == NULL) {
		fprintf(stderr, "new_quantifier_scope ERROR: Out of memory.\n");
		return false;
	}
	canonicalized->type = QuantifierType;
	canonicalized->quantifier.variable = new_scope.quantifier.variable;
	canonicalized->quantifier.operand = new_scope.quantifier.operand;
	canonicalized->quantifier.operand->reference_count++;

	/* TODO: move 'canonicalized' to the appropriate ancestor scope */

	return true;
}

inline bool canonicalize_term(const fol_term& src, fol_term& dst,
		array_map<unsigned int, unsigned int>& variable_map,
		array<unsigned int>& variables)
{
	dst.type = src.type;

	unsigned int index;
	switch (src.type) {
	case fol_term_type::CONSTANT:
		dst.constant = src.constant; return true;
	case fol_term_type::VARIABLE:
		index = variable_map.index_of(src.variable);
		if (index < variable_map.size) {
			dst.variable = variable_map.values[index];
			return variables.add(dst.variable);
		} else {
			fprintf(stderr, "canonicalize_term ERROR: Undefined variable.\n");
			return false;
		}
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "canonicalize_term ERROR: Unrecognized fol_term_type.\n");
	return false;
}

inline bool canonicalize_atom(
		const fol_atom& src, fol_scope& parent_scope,
		array_map<unsigned int, unsigned int>& variable_map,
		array<unsigned int>& variables)
{
	fol_formula* canonicalized = (fol_formula*) malloc(sizeof(fol_formula));
	if (canonicalized == NULL) {
		fprintf(stderr, "canonicalize_atom ERROR: Out of memory.\n");
		return false;
	}
	canonicalized->type = fol_formula_type::ATOM;

	if (!canonicalize_term(src.arg1, canonicalized->atom.arg1, variable_map, variables)
	 || !canonicalize_term(src.arg2, canonicalized->atom.arg2, variable_map, variables)) {
		free(canonicalized); return false;
	}

	/* TODO: move 'canonicalized' to the appropriate ancestor scope */

	return true;
}

bool canonicalize(const fol_formula& src, fol_scope& scope,
		array_map<unsigned int, unsigned int>& variable_map,
		array<unsigned int>& variables)
{
	switch (src.type) {
	case fol_formula_type::AND:
		if (scope.type == fol_formula_type::AND) {
			return canonicalize(*src.binary.left, scope, variable_map, variables)
				&& canonicalize(*src.binary.right, scope, variable_map, variables);
		} else if (scope.type == fol_formula_type::IF_THEN && scope.noncommutative.is_left) {
			return canonicalize(*src.binary.left, scope, variable_map, variables)
				&& canonicalize(*src.binary.right, scope, variable_map, variables);
		} else {
			return new_commutative_scope<fol_formula_type::AND>(src, scope, variable_map);
		}

	case fol_formula_type::OR:
		if (scope.type == fol_formula_type::OR) {
			return canonicalize(*src.binary.left, scope, variable_map, variables)
				&& canonicalize(*src.binary.right, scope, variable_map, variables);
		} else {
			return new_commutative_scope<fol_formula_type::OR>(src, scope, variable_map);
		}

	case fol_formula_type::IF_THEN:
		if (scope.type == fol_formula_type::IF_THEN && !scope.noncommutative.is_left) {
			scope.noncommutative.is_left = true;
			if (!canonicalize(*src.binary.left, scope, variable_map, variables)) return false;
			scope.noncommutative.is_left = false;
			return canonicalize(*src.binary.right, scope, variable_map, variables);
		} else {
			return new_noncommutative_scope<fol_formula_type::IF_THEN>(src, scope, variable_map);
		}

	case fol_formula_type::IFF:
		if (scope.type == fol_formula_type::IFF) {
			return canonicalize(*src.binary.left, scope, variable_map, variables)
				&& canonicalize(*src.binary.right, scope, variable_map, variables);
		} else {
			return new_commutative_scope<fol_formula_type::IFF>(src, scope, variable_map);
		}

	case fol_formula_type::NOT:
		if (src.unary.operand->type == fol_formula_type::NOT) {
			/* cancel this negation with the next one */
			return canonicalize(*src.unary.operand->unary.operand, scope, variable_map, variables);
		} else {
			return canonicalize_negation(src, scope, variable_map);
		}

	case fol_formula_type::FOR_ALL:
		return new_quantifier_scope<fol_formula_type::FOR_ALL>(src, scope, variable_map);

	case fol_formula_type::EXISTS:
		return new_quantifier_scope<fol_formula_type::EXISTS>(src, scope, variable_map);

	case fol_formula_type::ATOM:
		return canonicalize_atom(src.atom, scope, variable_map, variables);
	}
	fprintf(stderr, "canonicalize ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

template<typename Stream, typename... Printer>
bool print(const fol_formula& formula, Stream& out, Printer&&... printer)
{
	switch (formula.type) {
	case fol_formula_type::ATOM:
		if (!print(formula.atom.predicate, out, std::forward<Printer>(printer)...) || !print('(', out))
			return false;

		if (formula.atom.arg1.type != fol_term_type::NONE) {
			if (!print(formula.atom.arg1, out, std::forward<Printer>(printer)...))
				return false;
			if (formula.atom.arg2.type != fol_term_type::NONE) {
				if (!print(',', out) || !print(formula.atom.arg2, out, std::forward<Printer>(printer)...))
					return false;
			}
		}

		return print(')', out);

	case fol_formula_type::NOT:
		return print('~', out) && print(*formula.unary.operand, out, std::forward<Printer>(printer)...);

	case fol_formula_type::AND:
		return print('(', out) && print(*formula.binary.left, out, std::forward<Printer>(printer)...)
			&& print(" & ", out) && print(*formula.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case fol_formula_type::OR:
		return print('(', out) && print(*formula.binary.left, out, std::forward<Printer>(printer)...)
			&& print(" | ", out) && print(*formula.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case fol_formula_type::IF_THEN:
		return print('(', out) && print(*formula.binary.left, out, std::forward<Printer>(printer)...)
			&& print(" => ", out) && print(*formula.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case fol_formula_type::IFF:
		return print('(', out) && print(*formula.binary.left, out, std::forward<Printer>(printer)...)
			&& print(" <=> ", out) && print(*formula.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case fol_formula_type::FOR_ALL:
		return print("![", out) && print_variable(formula.quantifier.variable, out) && print("]:", out)
			&& print(*formula.quantifier.operand, out, std::forward<Printer>(printer)...);

	case fol_formula_type::EXISTS:
		return print("?[", out) && print_variable(formula.quantifier.variable, out) && print("]:", out)
			&& print(*formula.quantifier.operand, out, std::forward<Printer>(printer)...);
	}

	fprintf(stderr, "print ERROR: Unrecognized fol_formula_type.\n");
	return false;
}


/**
 * Code for tokenizing/lexing first-order logic formulas in TPTP-like format.
 */

enum class tptp_token_type {
	LBRACKET,
	RBRACKET,
	LPAREN,
	RPAREN,
	COMMA,
	COLON,

	AND,
	OR,
	NOT,
	IF_THEN,
	IFF,
	FOR_ALL,
	EXISTS,

	IDENTIFIER
};

typedef lexical_token<tptp_token_type> tptp_token;

template<typename Stream>
inline bool print(tptp_token_type type, Stream& stream) {
	switch (type) {
	case tptp_token_type::LBRACKET:
		return print('[', stream);
	case tptp_token_type::RBRACKET:
		return print(']', stream);
	case tptp_token_type::LPAREN:
		return print('(', stream);
	case tptp_token_type::RPAREN:
		return print(')', stream);
	case tptp_token_type::COMMA:
		return print(',', stream);
	case tptp_token_type::COLON:
		return print(':', stream);
	case tptp_token_type::AND:
		return print('&', stream);
	case tptp_token_type::OR:
		return print('|', stream);
	case tptp_token_type::NOT:
		return print('~', stream);
	case tptp_token_type::IF_THEN:
		return print("=>", stream);
	case tptp_token_type::IFF:
		return print("<=>", stream);
	case tptp_token_type::FOR_ALL:
		return print('!', stream);
	case tptp_token_type::EXISTS:
		return print('?', stream);
	case tptp_token_type::IDENTIFIER:
		return print("IDENTIFIER", stream);
	}
	fprintf(stderr, "print ERROR: Unknown tptp_token_type.\n");
	return false;
}

enum class tptp_lexer_state {
	DEFAULT,
	IDENTIFIER,
};

bool tptp_emit_symbol(array<tptp_token>& tokens, const position& start, char symbol) {
	switch (symbol) {
	case ',':
		return emit_token(tokens, start, start + 1, tptp_token_type::COMMA);
	case ':':
		return emit_token(tokens, start, start + 1, tptp_token_type::COLON);
	case '(':
		return emit_token(tokens, start, start + 1, tptp_token_type::LPAREN);
	case ')':
		return emit_token(tokens, start, start + 1, tptp_token_type::RPAREN);
	case '[':
		return emit_token(tokens, start, start + 1, tptp_token_type::LBRACKET);
	case ']':
		return emit_token(tokens, start, start + 1, tptp_token_type::RBRACKET);
	case '&':
		return emit_token(tokens, start, start + 1, tptp_token_type::AND);
	case '|':
		return emit_token(tokens, start, start + 1, tptp_token_type::OR);
	case '~':
		return emit_token(tokens, start, start + 1, tptp_token_type::NOT);
	case '!':
		return emit_token(tokens, start, start + 1, tptp_token_type::FOR_ALL);
	case '?':
		return emit_token(tokens, start, start + 1, tptp_token_type::EXISTS);
	default:
		fprintf(stderr, "tptp_emit_symbol ERROR: Unexpected symbol.\n");
		return false;
	}
}

template<typename Stream>
inline bool tptp_lex_symbol(array<tptp_token>& tokens, Stream& input, wint_t next, position& current)
{
	if (next == ',' || next == ':' || next == '(' || next == ')'
	 || next == '[' || next == ']' || next == '&' || next == '|'
	 || next == '~' || next == '!' || next == '?')
	{
		return tptp_emit_symbol(tokens, current, next);
	} else if (next == '=') {
		next = fgetwc(input);
		if (next != '>') {
			read_error("Expected '>' after '='", current);
			return false;
		} if (!emit_token(tokens, current, current + 2, tptp_token_type::IF_THEN))
			return false;
		current.column++;
	} else if (next == '<') {
		next = fgetwc(input);
		if (next != '=') {
			read_error("Expected '=' after '<'", current);
			return false;
		}
		next = fgetwc(input);
		if (next != '>') {
			read_error("Expected '>' after '='", current);
			return false;
		} if (!emit_token(tokens, current, current + 3, tptp_token_type::IF_THEN))
			return false;
		current.column += 2;
	} else {
		fprintf(stderr, "tptp_lex_symbol ERROR: Unrecognized symbol.\n");
		return false;
	}
	return true;
}

inline bool append_to_token(
	array<char>& token, wint_t next, std::mbstate_t& shift)
{
	if (!token.ensure_capacity(token.length + MB_CUR_MAX))
		return false;
	size_t written = wcrtomb(token.data + token.length, next, &shift);
	if (written == static_cast<std::size_t>(-1))
		return false;
	token.length += written;
	return true;
}

template<typename Stream>
bool tptp_lex(array<tptp_token>& tokens, Stream& input) {
	position start = position(1, 1);
	position current = position(1, 1);
	tptp_lexer_state state = tptp_lexer_state::DEFAULT;
	array<char> token = array<char>(1024);

	std::mbstate_t shift = {0};
	wint_t next = fgetwc(input);
	bool new_line = false;
	while (next != WEOF) {
		switch (state) {
		case tptp_lexer_state::IDENTIFIER:
			if (next == ',' || next == ':' || next == '(' || next == ')'
			 || next == '[' || next == ']' || next == '&' || next == '|'
			 || next == '~' || next == '!' || next == '?' || next == '='
			 || next == '<')
			{
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER)
				 || !tptp_lex_symbol(tokens, input, next, current))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear(); shift = {0};
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear(); shift = {0};
				new_line = (next == '\n');
			} else {
				if (!append_to_token(token, next, shift)) return false;
			}
			break;

		case tptp_lexer_state::DEFAULT:
			if (next == ',' || next == ':' || next == '(' || next == ')'
			 || next == '[' || next == ']' || next == '&' || next == '|'
			 || next == '~' || next == '!' || next == '?' || next == '='
			 || next == '<')
			{
				if (!tptp_lex_symbol(tokens, input, next, current))
					return false;
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
			} else {
				if (!append_to_token(token, next, shift)) return false;
				state = tptp_lexer_state::IDENTIFIER;
				start = current;
			}
			break;
		}

		if (new_line) {
			current.line++;
			current.column = 1;
			new_line = false;
		} else current.column++;
		next = fgetwc(input);
	}

	if (state == tptp_lexer_state::IDENTIFIER)
		return emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER);
	return true;
}


/**
 * Recursive-descent parser for first-order logic formulas in TPTP-like format.
 */

bool tptp_interpret_unary_formula(
	const array<tptp_token>&,
	unsigned int&, fol_formula&,
	hash_map<string, unsigned int>&,
	array_map<string, unsigned int>&);
bool tptp_interpret(
	const array<tptp_token>&,
	unsigned int&, fol_formula&,
	hash_map<string, unsigned int>&,
	array_map<string, unsigned int>&);

bool tptp_interpret_argument_list(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	array<fol_term>& terms)
{
	if (!expect_token(tokens, index, tptp_token_type::LPAREN,
			"opening parenthesis for list of arguments in atomic formula"))
		return false;
	index++;

	while (true) {
		if (!terms.ensure_capacity(terms.length + 1)
		 || !expect_token(tokens, index, tptp_token_type::IDENTIFIER,
				"identifier in list of arguments in atomic formula"))
			return false;

		bool contains;
		fol_term& next_term = terms[terms.length];
		unsigned int variable = variables.get(tokens[index].text, contains);
		if (contains) {
			/* this argument is a variable */
			next_term.variable = variable;
			next_term.type = fol_term_type::VARIABLE;
		} else {
			if (!get_token(tokens[index].text, next_term.constant, names))
				return false;
			next_term.type = fol_term_type::CONSTANT;
		}
		terms.length++;
		index++;

		if (index >= tokens.length) {
			read_error("Unexpected end of input", tokens.last().end);
			return false;
		} else if (tokens[index].type == tptp_token_type::RPAREN) {
			index++;
			return true;
		} else if (tokens[index].type != tptp_token_type::COMMA) {
			read_error("Unexpected symbol. Expected a comma", tokens[index].start);
			return false;
		}
		index++;
	}
}

bool tptp_interpret_variable_list(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	if (!expect_token(tokens, index, tptp_token_type::LBRACKET, "left bracket for list of quantified variables"))
		return false;
	index++;

	while (true) {
		if (!expect_token(tokens, index, tptp_token_type::IDENTIFIER, "variable in list of quantified variables"))
			return false;
		if (names.table.contains(tokens[index].text)) {
			fprintf(stderr, "WARNING at %d:%d: Variable '", tokens[index].start.line, tokens[index].start.column);
			print(tokens[index].text, stderr); print("' shadows previously declared identifier.\n", stderr);
		} if (variables.contains(tokens[index].text)) {
			read_error("Variable redeclared", tokens[index].start);
			return false;
		} if (!variables.ensure_capacity(variables.size + 1)) {
			return false;
		}
		variables.keys[variables.size] = tokens[index].text;
		variables.values[variables.size] = variables.size + 1;
		variables.size++;
		index++;

		if (index >= tokens.length) {
			read_error("Unexpected end of input", tokens.last().end);
			return false;
		} else if (tokens[index].type == tptp_token_type::RBRACKET) {
			index++;
			return true;
		} else if (tokens[index].type != tptp_token_type::COMMA) {
			read_error("Unexpected symbol. Expected a comma", tokens[index].start);
			return false;
		}
		index++;
	}
}

template<fol_formula_type QuantifierType>
bool tptp_interpret_quantifier(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	unsigned int old_variable_count = variables.size;
	if (!tptp_interpret_variable_list(tokens, index, names, variables)
	 || !expect_token(tokens, index, tptp_token_type::COLON, "colon for quantified formula"))
		return false;
	index++;

	fol_formula* operand = (fol_formula*) malloc(sizeof(fol_formula));
	operand->reference_count = 1;
	if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
		free(operand); return false;
	}

	fol_formula* inner = operand;
	for (unsigned int i = variables.size - 1; i > old_variable_count; i--) {
		fol_formula* quantified = (fol_formula*) malloc(sizeof(fol_formula));
		if (quantified == NULL) {
			fprintf(stderr, "tptp_interpret_unary_formula ERROR: Out of memory.\n");
			free(*inner); free(inner);
		}
		quantified->quantifier.variable = variables.values[i];
		quantified->quantifier.operand = inner;
		quantified->type = QuantifierType;
		quantified->reference_count = 1;
		inner = quantified;
	}

	formula.quantifier.variable = variables.values[old_variable_count];
	formula.quantifier.operand = inner;
	formula.type = QuantifierType;
	formula.reference_count = 1;

	for (unsigned int i = old_variable_count; i < variables.size; i++)
		free(variables.keys[i]);
	variables.size = old_variable_count;
	return true;
}

bool tptp_interpret_unary_formula(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	if (index >= tokens.length) {
		fprintf(stderr, "ERROR: Unexpected end of input.\n");
		return false;

	} else if (tokens[index].type == tptp_token_type::NOT) {
		/* this is a negation of the form ~U */
		index++;
		fol_formula* operand = (fol_formula*) malloc(sizeof(fol_formula));
		operand->reference_count = 1;
		if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
			free(operand); return false;
		}

		formula.unary.operand = operand;
		formula.type = fol_formula_type::NOT;
		formula.reference_count = 1;

	} else if (tokens[index].type == tptp_token_type::LPAREN) {
		/* these are just grouping parenthesis of the form (F) */
		index++;
		if (!tptp_interpret(tokens, index, formula, names, variables)) {
			return false;
		} if (!expect_token(tokens, index, tptp_token_type::RPAREN, "closing parenthesis")) {
			free(formula); return false;
		}
		index++;

	} else if (tokens[index].type == tptp_token_type::FOR_ALL) {
		/* this is a universal quantifier of the form ![v_1,...,v_n]:U */
		index++;
		if (!tptp_interpret_quantifier<fol_formula_type::FOR_ALL>(tokens, index, formula, names, variables))
			return false;

	} else if (tokens[index].type == tptp_token_type::EXISTS) {
		/* this is an existential quantifier of the form ?[v_1,...,v_n]:U */
		index++;
		if (!tptp_interpret_quantifier<fol_formula_type::EXISTS>(tokens, index, formula, names, variables))
			return false;

	} else if (tokens[index].type == tptp_token_type::IDENTIFIER) {
		/* this is an atomic formula of the form P(T_1,...,T_n) */
		unsigned int predicate;
		if (variables.contains(tokens[index].text)) {
			fprintf(stderr, "WARNING at %d:%d: Predicate '", tokens[index].start.line, tokens[index].start.column);
			print(tokens[index].text, stderr); print("' is also a variable.\n", stderr);
		} if (!get_token(tokens[index].text, predicate, names))
			return false;
		index++;

		array<fol_term> terms = array<fol_term>(2);
		if (!tptp_interpret_argument_list(tokens, index, names, variables, terms))
			return false;
		if (terms.length == 0) {
			formula.atom.arg1.type = fol_term_type::NONE;
			formula.atom.arg2.type = fol_term_type::NONE;
		} else if (terms.length == 1) {
			formula.atom.arg1 = terms[0];
			formula.atom.arg2.type = fol_term_type::NONE;
		} else if (terms.length == 2) {
			formula.atom.arg1 = terms[0];
			formula.atom.arg2 = terms[1];
		} else {
			read_error("Atomic formulas with arity greater than 2 are not supported", tokens[index - 1].end);
			return false;
		}
		formula.atom.predicate = predicate;
		formula.type = fol_formula_type::ATOM;
		formula.reference_count = 1;

	} else {
		read_error("Unexpected symbol. Expected a unary formula", tokens[index].start);
		return false;
	}
	return true;
}

template<fol_formula_type OperatorType>
inline bool tptp_interpret_binary_formula(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	fol_formula* left)
{
	fol_formula* right = (fol_formula*) malloc(sizeof(fol_formula));
	if (right == NULL) {
		fprintf(stderr, "tptp_interpret_binary_formula ERROR: Out of memory.\n");
		free(*left); free(left); return false;
	} else if (!tptp_interpret_unary_formula(tokens, index, *right, names, variables)) {
		free(*left); free(left); free(right); return false;
	}

	formula.binary.left = left;
	formula.binary.right = right;
	formula.type = OperatorType;
	formula.reference_count = 1;
	return true;
}

template<fol_formula_type OperatorType> struct tptp_operator_type { };
template<> struct tptp_operator_type<fol_formula_type::AND> { static constexpr tptp_token_type type = tptp_token_type::AND; };
template<> struct tptp_operator_type<fol_formula_type::OR> { static constexpr tptp_token_type type = tptp_token_type::OR; };

template<fol_formula_type OperatorType>
bool tptp_interpret_binary_sequence(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	fol_formula* left)
{
	while (true) {
		fol_formula* next = (fol_formula*) malloc(sizeof(fol_formula));
		if (next == NULL) {
			fprintf(stderr, "tptp_interpret_binary_sequence ERROR: Out of memory.\n");
			free(*left); free(left); return false;
		} else if (!tptp_interpret_unary_formula(tokens, index, *next, names, variables)) {
			free(*left); free(left); free(next); return false;
		}

		if (index < tokens.length && tokens[index].type == tptp_operator_type<OperatorType>::type) {
			index++;
			fol_formula* parent = (fol_formula*) malloc(sizeof(fol_formula));
			if (parent == NULL) {
				fprintf(stderr, "tptp_interpret_binary_sequence ERROR: Out of memory.\n");
				free(*left); free(left); free(*next); free(next); return false;
			}

			parent->binary.left = left;
			parent->binary.right = next;
			parent->type = OperatorType;
			parent->reference_count = 1;
			left = parent;
		} else {
			formula.binary.left = left;
			formula.binary.right = next;
			formula.type = OperatorType;
			formula.reference_count = 1;
			return true;
		}
	}
}

bool tptp_interpret(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	fol_formula* left = (fol_formula*) malloc(sizeof(fol_formula));
	if (left == NULL) {
		fprintf(stderr, "tptp_interpret ERROR: Out of memory.\n");
		return false;
	} if (!tptp_interpret_unary_formula(tokens, index, *left, names, variables)) {
		return false;
	}

	if (index >= tokens.length) {
		move(*left, formula); free(left);
		return true;
	} else if (tokens[index].type == tptp_token_type::AND) {
		index++;
		if (!tptp_interpret_binary_sequence<fol_formula_type::AND>(tokens, index, formula, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::OR) {
		index++;
		if (!tptp_interpret_binary_sequence<fol_formula_type::OR>(tokens, index, formula, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::IF_THEN) {
		index++;
		if (!tptp_interpret_binary_formula<fol_formula_type::IF_THEN>(tokens, index, formula, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::IFF) {
		index++;
		if (!tptp_interpret_binary_formula<fol_formula_type::IFF>(tokens, index, formula, names, variables, left))
			return false;
	} else {
		move(*left, formula); free(left);
	}
	return true;
}

#endif /* FIRST_ORDER_LOGIC_H_ */
