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
bool operator == (const fol_formula&, const fol_formula&);
bool operator == (const fol_formula&, const fol_formula&);

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

inline bool operator == (const fol_term& first, const fol_term& second) {
	if (first.type != second.type) return false;
	switch (first.type) {
	case fol_term_type::VARIABLE:
		return first.variable == second.variable;
	case fol_term_type::CONSTANT:
		return first.constant == second.constant;
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

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
	EXISTS,

	TRUE,
	FALSE
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

inline bool operator == (const fol_atom& first, const fol_atom& second) {
	return first.predicate == second.predicate
		&& first.arg1 == second.arg1
		&& first.arg2 == second.arg2;
}

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

	fol_formula(fol_formula_type type) : type(type), reference_count(1) { }
	~fol_formula() { free_helper(); }

	static inline void move(const fol_formula& src, fol_formula& dst);
	static inline void free(fol_formula& formula) { formula.free_helper(); }

	void free_helper();
};

thread_local fol_formula FOL_TRUE(fol_formula_type::TRUE);
thread_local fol_formula FOL_FALSE(fol_formula_type::FALSE);

inline bool operator == (const fol_unary_formula& first, const fol_unary_formula& second) {
	return first.operand == second.operand
		|| *first.operand == *second.operand;
}

inline bool operator == (const fol_binary_formula& first, const fol_binary_formula& second) {
	return (first.left == second.left || *first.left == *second.left)
		&& (first.right == second.right || *first.right == *second.right);
}

inline bool operator == (const fol_quantifier& first, const fol_quantifier& second) {
	return first.variable == second.variable
		&& (first.operand == second.operand || *first.operand == *second.operand);
}

bool operator == (const fol_formula& first, const fol_formula& second)
{
	if (first.type != second.type) return false;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return first.atom == second.atom;
	case fol_formula_type::NOT:
		return first.unary == second.unary;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
		return first.binary == second.binary;
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return first.quantifier == second.quantifier;
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return true;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

inline bool operator != (const fol_formula& first, const fol_formula& second) {
	return !(first == second);
}

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
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return;
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

inline void fol_formula::free_helper() {
	reference_count--;
	if (reference_count == 0) {
		switch (type) {
		case fol_formula_type::ATOM:
			core::free(atom); return;
		case fol_formula_type::NOT:
			core::free(unary); return;
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IF_THEN:
		case fol_formula_type::IFF:
			core::free(binary); return;
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			core::free(quantifier); return;
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return;
		}
		fprintf(stderr, "fol_formula.free_helper ERROR: Unrecognized fol_formula_type.\n");
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

	case fol_formula_type::TRUE:
		return print('T', out);

	case fol_formula_type::FALSE:
		return print('F', out);

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
 * Below is code for canonicalizing first-order formulas.
 */


struct canonicalizer { };

int_fast8_t compare(
		const fol_formula&,
		const fol_formula&,
		const canonicalizer&);

inline int_fast8_t compare(
		const fol_term& first,
		const fol_term& second,
		const canonicalizer& sorter)
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
		const canonicalizer& sorter)
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
		const canonicalizer& sorter)
{
	return compare(*first.operand, *second.operand, sorter);
}

inline int_fast8_t compare(
		const fol_binary_formula& first,
		const fol_binary_formula& second,
		const canonicalizer& sorter)
{
	int_fast8_t result = compare(*first.left, *second.left, sorter);
	if (result != 0) return result;
	return compare(*first.right, *second.right, sorter);
}

inline int_fast8_t compare(
		const fol_quantifier& first,
		const fol_quantifier& second,
		const canonicalizer& sorter)
{
	/* TODO: how do we compare variables? */
	return compare(*first.operand, *second.operand, sorter);
}

int_fast8_t compare(
		const fol_formula& first,
		const fol_formula& second,
		const canonicalizer& sorter)
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
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

bool less_than(
		const fol_formula* first,
		const fol_formula* second,
		const canonicalizer& sorter)
{
	return compare(*first, *second, sorter) < 0;
}

inline bool relabel_variables(fol_term& term,
		array_map<unsigned int, unsigned int>& variable_map)
{
	unsigned int index;
	switch (term.type) {
	case fol_term_type::CONSTANT:
	case fol_term_type::NONE:
		return true;
	case fol_term_type::VARIABLE:
		index = variable_map.index_of(term.variable);
		if (index < variable_map.size) {
			term.variable = variable_map.values[index];
			return true;
		} else {
			fprintf(stderr, "relabel_variables ERROR: Undefined variable.\n");
			return false;
		}
	}
	fprintf(stderr, "relabel_variables ERROR: Unrecognized fol_term_type.\n");
	return false;
}


/* forward declarations */
bool relabel_variables(fol_formula&, array_map<unsigned int, unsigned int>&);


inline bool new_variable(unsigned int src, unsigned int& dst,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!variable_map.ensure_capacity(variable_map.size + 1))
		return false;
	unsigned int index = variable_map.index_of(src);
	if (index < variable_map.size) {
		fprintf(stderr, "new_variable ERROR: Multiple declaration of variable %u.\n", src);
		return false;
	}
	variable_map.keys[index] = src;
	dst = variable_map.size + 1;
	variable_map.values[index] = dst;
	variable_map.size++;
	return true;
}

inline bool relabel_variables(fol_atom& atom,
		array_map<unsigned int, unsigned int>& variable_map)
{
	return relabel_variables(atom.arg1, variable_map)
		&& relabel_variables(atom.arg2, variable_map);
}

inline bool relabel_variables(fol_quantifier& quantifier,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!new_variable(quantifier.variable, quantifier.variable, variable_map)
	 || !relabel_variables(*quantifier.operand, variable_map))
		return false;

	variable_map.size--;
	return true;
}

bool relabel_variables(fol_formula& formula,
		array_map<unsigned int, unsigned int>& variable_map)
{
	switch (formula.type) {
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
		return relabel_variables(*formula.binary.left, variable_map)
			&& relabel_variables(*formula.binary.right, variable_map);
	case fol_formula_type::NOT:
		return relabel_variables(*formula.unary.operand, variable_map);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return relabel_variables(formula.quantifier, variable_map);
	case fol_formula_type::ATOM:
		return relabel_variables(formula.atom, variable_map);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return true;
	}
	fprintf(stderr, "relabel_variables ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

inline bool relabel_variables(fol_formula& formula) {
	array_map<unsigned int, unsigned int> variable_map(16);
	return relabel_variables(formula, variable_map);
}


/* forward declarations */
struct fol_scope;


template<unsigned int Depth>
using fol_scope_child = typename std::conditional<Depth == 0, fol_formula_type*, fol_scope<Depth - 1>>::type;

template<unsigned int Depth>
using fol_scope_ptr = typename std::conditional<Depth == 0, fol_formula_type*, fol_scope<Depth - 1>*>::type;

template<unsigned int Depth>
struct fol_commutative_scope {
	array<fol_scope_child<Depth>> children;

	fol_commutative_scope() : children(8) { }

	static inline void move(const fol_commutative_scope& src, fol_commutative_scope& dst) {
		core::move(src.children, dst.children);
	}
};

template<unsigned int Depth>
struct fol_noncommutative_scope {
	array<fol_scope_child<Depth>> left;
	array<fol_scope_child<Depth>> right;
	bool is_left;

	fol_noncommutative_scope() : left(8), right(8), is_left(true) { }

	static inline void move(const fol_noncommutative_scope& src, fol_noncommutative_scope& dst) {
		core::move(src.left, dst.left);
		core::move(src.right, dst.right);
		dst.is_left = src.is_left;
	}
};

template<unsigned int Depth>
struct fol_quantifier_scope {
	fol_scope_ptr<Depth> operand;
	unsigned int variable;

	static inline void move(const fol_quantifier_scope& src, fol_quantifier_scope& dst) {
		dst.operand = src.operand;
		dst.variable = src.variable;
	}
};

template<unsigned int Depth>
struct fol_scope {
	fol_formula_type type;
	array<unsigned int> variables;
	union {
		fol_atom atom;
		fol_scope_ptr unary;
		fol_commutative_scope commutative;
		fol_noncommutative_scope noncommutative;
		fol_quantifier_scope quantifier;
	};

	fol_scope() : fol_scope(fol_formula_type::TRUE) { }

	fol_scope(fol_formula_type type) : variables(8) {
		if (!init_helper(type))
			exit(EXIT_FAILURE);
	}

	~fol_scope() { free_helper(); }

	static inline void move(const fol_scope<Depth>& src, fol_scope<Depth>& dst)
	{
		switch (src.type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			fol_commutative_scope::move(src.commutative, dst.commutative);
		case fol_formula_type::IF_THEN:
			fol_commutative_scope::move(src.noncommutative, dst.noncommutative);
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			fol_commutative_scope::move(src.quantifier, dst.quantifier);
		case fol_formula_type::ATOM:
			fol_atom::move(src.atom, dst.atom);
		case fol_formula_type::NOT:
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return;
		}
		fprintf(stderr, "fol_scope.move ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void free(fol_scope<Depth>& scope) {
		scope.free_helper();
		core::free(scope.variables);
	}

private:
	bool init_helper(fol_formula_type scope_type) {
		type = scope_type;
		switch (type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			new (&commutative) fol_commutative_scope(); return true;
		case fol_formula_type::IF_THEN:
			new (&noncommutative) fol_noncommutative_scope(); return true;
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			new (&quantifier) fol_quantifier_scope(); return true;
		case fol_formula_type::ATOM:
			new (&atom) fol_atom(); return true;
		case fol_formula_type::NOT:
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return true;
		}
		fprintf(stderr, "fol_scope.init_helper ERROR: Unrecognized fol_formula_type.\n");
		return false;
	}

	bool free_helper() {
		switch (type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			commutative.~fol_commutative_scope(); return true;
		case fol_formula_type::IF_THEN:
			noncommutative.~fol_noncommutative_scope(); return true;
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			quantifier.~fol_quantifier_scope(); return true;
		case fol_formula_type::ATOM:
			atom.~fol_atom(); return true;
		case fol_formula_type::NOT:
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return true;
		}
		fprintf(stderr, "fol_scope.free_helper ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}

	friend bool init(fol_scope<Depth>&, fol_formula_type);
};

inline fol_formula* scope_to_formula(fol_formula* scope) {
	return scope;
}

template<unsigned int Depth>
inline fol_formula* scope_to_formula(fol_scope<Depth>* scope) {
	return scope_to_formula(*scope);
}

template<fol_formula_type ScopeType, unsigned int Depth>
inline fol_formula* scope_to_formula(array<fol_scope_child<Depth>>& scope)
{
	fol_formula* inner = scope.first();
	for (unsigned int i = 1; i < scope.length; i++) {
		fol_formula* right = scope_to_formula(scope[i]);
		if (right == NULL) {
			free(*inner); if (inner->reference_count == 0); free(inner);
			return NULL;
		}

		fol_formula* new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "canonicalize_scope: Out of memory.\n");
			free(*inner); if (inner->reference_count == 0) free(inner);
			free(*right); if (right->reference_count == 0) free(right);
			return NULL;
		}
		new_formula->type = ScopeType;
		new_formula->binary.left = inner;
		new_formula->binary.right = right;
		new_formula->reference_count = 1;
		/* the reference count of `inner` is incremented and decremented
		   simultaneously here, so we don't change it */
		inner = new_formula;
	}

	return inner;
}

template<fol_formula_type QuantifierType, unsigned int Depth>
inline fol_formula* scope_to_formula(fol_quantifier_scope<Depth>& scope)
{
	static_assert(ScopeType == fol_formula_type::FOR_ALL
			   || ScopeType == fol_formula_type::EXISTS,
			"ScopeType is not a quantifier.");

	fol_formula* operand = scope_to_formula(scope.operand);
	if (operand == NULL) return NULL;

	new_formula = (fol_formula*) malloc(sizeof(fol_formula));
	if (new_formula == NULL) {
		fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
		return NULL;
	}
	new_formula->type = fol_formula_type::FOR_ALL;
	new_formula->reference_count = 1;
	new_formula->quantifier.variable = quantifier.variable;
	new_formula->quantifier.operand = operand;
	return new_formula;
}

template<unsigned int Depth>
inline fol_formula* scope_to_formula(fol_scope<Depth>& scope)
{
	fol_formula* new_formula;
	fol_formula* left; fol_formula* right;
	switch (scope.type) {
	case fol_formula_type::AND:
		return scope_to_formula<fol_formula_type::AND>(scope.commutative.children);
	case fol_formula_type::OR:
		return scope_to_formula<fol_formula_type::OR>(scope.commutative.children);
	case fol_formula_type::IFF:
		return scope_to_formula<fol_formula_type::IFF>(scope.commutative.children);
	case fol_formula_type::IF_THEN:
		left = scope_to_formula<fol_formula_type::AND>(scope.noncommutative.left);
		if (left == NULL) return NULL;
		right = scope_to_formula<fol_formula_type::AND>(scope.noncommutative.right);
		if (right == NULL) {
			free(*left); if (left->reference_count == 0) free(left);
			return NULL;
		}

		new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			return NULL;
		}
		new_formula->type = fol_formula_type::IF_THEN;
		new_formula->reference_count = 1;
		new_formula->binary.left = left;
		new_formula->binary.right = right;
		return new_formula;
	case fol_formula_type::NOT:
		left = scope_to_formula(scope.unary);
		if (left == NULL) return NULL;

		new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			return NULL;
		}
		new_formula->type = fol_formula_type::NOT;
		new_formula->reference_count = 1;
		new_formula->unary.operand = left;
		return new_formula;
	case fol_formula_type::FOR_ALL:
		return scope_to_formula<fol_formula_type::FOR_ALL>(scope.quantifier);
	case fol_formula_type::EXISTS:
		return scope_to_formula<fol_formula_type::EXISTS>(scope.quantifier);
	case fol_formula_type::ATOM:
		new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			return NULL;
		}
		new_formula->type = fol_formula_type::ATOM;
		new_formula->reference_count = 1;
		new_formula->atom = scope.atom;
		return new_formula;
	case fol_formula_type::TRUE:
		FOL_TRUE.reference_count++;
		return &FOL_TRUE;
	case fol_formula_type::FALSE:
		FOL_FALSE.reference_count++;
		return &FOL_FALSE;
	}
	fprintf(stderr, "scope_to_formula ERROR: Unrecognized fol_formula_type.\n");
	return NULL;
}

inline void move_variables(
		const array<unsigned int>& formula_variables,
		array<unsigned int>& scope_variables)
{
	array<unsigned int> variable_union = array<unsigned int>(
			max(scope_variables.capacity, formula_variables.length + scope_variables.length));
	set_union(variable_union.data, variable_union.length,
			formula_variables.data, formula_variables.length,
			scope_variables.data, scope_variables.length);
	swap(variable_union, scope_variables);
}

template<unsigned int Depth>
bool add_to_commutative_scope(fol_scope_child<Depth>& formula, fol_commutative_scope<Depth>& scope)
{
	/* add `formula` into the correct (sorted) position in `scope.children` */
	unsigned int i = 0;
	for (; i < scope.children.length; i++) {
		auto result = compare(formula, scope.children[i], canonicalizer());
		if (result < 0) {
			break;
		} else if (result == 0) {
			free(*formula); if (formula->reference_count == 0) free(formula);
			return true;
		}
	}

	/* `*formula` is unique, so add it at index `i` */
	if (scope.children.ensure_capacity(scope.children.length + 1))
		return false;
	shift_right(scope.children.data, scope.children.length, i);
	move(formula, scope.children[i]);
	scope.children.length++;
	return true;
}

template<unsigned int Depth>
void merge_scopes(array<fol_scope_child<Depth>>& dst,
	const array<fol_scope_child<Depth>>& first,
	const array<fol_scope_child<Depth>>& second)
{
	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[i], canonicalizer());
		if (result == 0) {
			dst[dst.length] = first[i];
			free(*second[j]); if (second[j]->reference_count == 0) free(second[j]); /* TODO: this needs to be a deep free */
			dst.length++; i++; j++;
		} else if (result < 0) {
			dst[dst.length] = first[i];
			dst.length++; i++;
		} else {
			dst[dst.length] = second[i];
			dst.length++; j++;
		}
	}

	while (i < first.length) {
		dst[dst.length] = first[i];
		dst.length++; i++;
	} while (j < second.length) {
		dst[dst.length] = second[j];
		dst.length++; j++;
	}
}

bool canonicalize_commutative(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope right;
	if (!canonicalize(*src.left, out, variable_map)
	 || !canonicalize(*src.right, right, variable_map))
		return false;

	/* consider the case where `*src.left` and `*src.right` are identical */
	if (out == right) return true;

	if (out.type == fol_formula_type::AND) {
		if (right.type == fol_formula_type::AND) {
			/* merge the two scopes */
			array<fol_formula*> both = array<fol_formula*>(out.commutative.children.length + right.commutative.children.length);
			merge_scopes(both, out.commutative.children, right.commutative.children);
			swap(both, out.commutative.children);
			move_variables(right.variables, out.variables);
		} else if (right.type == fol_formula_type::NOT && *right.unary == out) {
			/* we have the case `A & ~A` */
			free(*out); free(*right);
			out.type = fol_formula_type::FALSE;
		} else {
			fol_formula* right_formula = scope_to_formula(right);
			if (right_formula == NULL || !add_to_commutative_scope(right_formula, out.commutative))
				return false;
			move_variables(right.variables, out.variables);
		}

	} else if (out.type == fol_formula_type::NOT && *out.unary == right) {
		/* we have the case `A & ~A` */
		free(*out); free(*right);
		out.type = fol_formula_type::FALSE;
	} else {
		fol_formula* left_formula = scope_to_formula(out);
		if (left_formula == NULL) return false;

		if (right.type == fol_formula_type::AND) {
			if (!add_to_commutative_scope(left_formula, right.commutative)) return false;
			move_variables(out.variables, right.variables);
			swap(right, out);
		} else {
			fol_formula* right_formula = scope_to_formula(right);
			if (right_formula == NULL) {
				free(*left_formula); if (left_formula->reference_count == 0) free(left_formula);
				return false;
			}

			if (!init(out, fol_formula_type::AND)) {
				free(*left_formula); if (left_formula->reference_count == 0) free(left_formula);
				free(*right_formula); if (right_formula->reference_count == 0) free(right_formula);
				return false;
			} else if (less_than(*left_formula, *right_formula, canonicalizer())) {
				out.commutative.children[0] = left_formula;
				out.commutative.children[1] = right_formula;
			} else {
				out.commutative.children[0] = right_formula;
				out.commutative.children[1] = left_formula;
			}
			out.commutative.children.length = 2;
			move_variables(right.variables, out.variables);
		}
	}
	return true;
}

bool canonicalize(const fol_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	switch (src.type) {
	case fol_formula_type::AND:
		canonicalize(*src.binary.left, 
	}
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
		} if (!emit_token(tokens, current, current + 3, tptp_token_type::IFF))
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
		if (tokens[index].text == "T") {
			/* this the constant true */
			formula.type = fol_formula_type::TRUE;
			formula.reference_count = 1;
			index++;
		} else if (tokens[index].text == "F") {
			/* this the constant false */
			formula.type = fol_formula_type::FALSE;
			formula.reference_count = 1;
			index++;
		} else {
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
		}

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
		free(left); return false;
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
