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
	FALSE /* our canonicalization code assumes that FALSE is the last element of this enum */
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


struct fol_commutative_scope {
	array<fol_scope> children;
	array<fol_scope> negated;

	fol_commutative_scope() : children(4), negated(4) { }

	static inline void move(const fol_commutative_scope& src, fol_commutative_scope& dst) {
		core::move(src.children, dst.children);
		core::move(src.negated, dst.negated);
	}
};

struct fol_noncommutative_scope {
	array<fol_scope> left, left_negated;
	array<fol_scope> right, right_negated;

	fol_noncommutative_scope() : left(4), left_negated(4), right(4), right_negated(4) { }

	static inline void move(const fol_noncommutative_scope& src, fol_noncommutative_scope& dst) {
		core::move(src.left, dst.left);
		core::move(src.left_negated, dst.left_negated);
		core::move(src.right, dst.right);
		core::move(src.right_negated, dst.right_negated);
	}
};

struct fol_quantifier_scope {
	fol_scope* operand;
	unsigned int variable;

	static inline void move(const fol_quantifier_scope& src, fol_quantifier_scope& dst) {
		dst.operand = src.operand;
		dst.variable = src.variable;
	}
};

struct fol_scope {
	fol_formula_type type;
	array<unsigned int> variables;
	union {
		fol_atom atom;
		fol_scope* unary;
		fol_commutative_scope commutative;
		fol_noncommutative_scope noncommutative;
		fol_quantifier_scope quantifier;
	};

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

	static inline void free(fol_scope& scope) {
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

	friend bool init(fol_scope&, fol_formula_type);
};

template<bool Negated>
inline fol_formula* scope_to_formula(const fol_scope& scope);

template<>
inline fol_formula* scope_to_formula<false>(const fol_scope& scope) {
	return scope_to_formula(scope);
}

template<>
inline fol_formula* scope_to_formula<true>(const fol_scope& negated)
{
	fol_formula* not = (fol_formula*) malloc(sizeof(fol_formula));
	if (not == NULL) {
		fprintf(stderr, "negated_scope_to_formula ERROR: Out of memory.\n");
		return false;
	}
	not->type = fol_formula_type::NOT;
	not->reference_count = 1;

	not->unary.operand = scope_to_formula(negated);
	if (not->unary.operand == NULL) {
		free(not); return NULL;
	}
	return not;
}

template<fol_formula_type ScopeType, bool Negated>
inline fol_formula* scope_to_formula(const fol_scope* scope,
		unsigned int scope_length, fol_formula* first)
{
	fol_formula* left = first;

	for (unsigned int i = 0; i < scope_length; i++) {
		fol_formula* right = scope_to_formula<Negated>(scope[i]);
		if (right == NULL) {
			free(*left); if (left->reference_count == 0) free(left);
			return NULL;
		}

		fol_formula* new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			free(*left); if (left->reference_count == 0) free(left);
			free(*right); if (right->reference_count == 0) free(right);
			return NULL;
		}
		new_formula->type = ScopeType;
		new_formula->binary.left = left;
		new_formula->binary.right = right;
		new_formula->reference_count = 1;
		/* the reference counts of `left` and `right` are incremented and
		   decremented simultaneously here, so we don't change them */
		left = new_formula;
	}

	return left;
}

template<fol_formula_type ScopeType>
inline fol_formula* scope_to_formula(
		const array<fol_scope>& scope,
		const array<fol_scope>& negated)
{
	if (negated.length == 0) {
		fol_formula* first = scope_to_formula<false>(scope[0]);
		if (first == NULL) return NULL;
		return scope_to_formula<ScopeType, false>(scope.data + 1, scope.length - 1, first);
	} else {
		fol_formula* first = scope_to_formula<true>(negated[0]);
		if (first == NULL) return NULL;
		first = scope_to_formula<ScopeType, true>(negated.data + 1, negated.length - 1, first);
		if (first == NULL) return NULL;
		return scope_to_formula<ScopeType, false>(scope.data, scope.length, first);
	}
}

template<fol_formula_type QuantifierType>
inline fol_formula* scope_to_formula(const fol_quantifier_scope& scope)
{
	static_assert(ScopeType == fol_formula_type::FOR_ALL
			   || ScopeType == fol_formula_type::EXISTS,
			"ScopeType is not a quantifier.");

	fol_formula* operand = scope_to_formula(*scope.operand);
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

inline fol_formula* scope_to_formula(const fol_scope& scope)
{
	fol_formula* new_formula;
	fol_formula* left; fol_formula* right;
	switch (scope.type) {
	case fol_formula_type::AND:
		return scope_to_formula<fol_formula_type::AND>(scope.commutative.children, scope.commutative.negated);
	case fol_formula_type::OR:
		return scope_to_formula<fol_formula_type::OR>(scope.commutative.children, scope.commutative.negated);
	case fol_formula_type::IFF:
		left = scope_to_formula<fol_formula_type::IFF>(scope.commutative.children, scope.commutative.negated);
		if (left == NULL) return NULL;

		if (scope.commutative.children.length > 0 && scope.commutative.children.last().type == fol_formula_type::FALSE) {
			new_formula = (fol_formula*) malloc(sizeof(fol_formula));
			if (new_formula == NULL) {
				fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			}
			new_formula->type = fol_formula_type::NOT;
			new_formula->reference_count = 1;
			new_formula->unary.opearnd = left;
			return new_formula;
		} else {
			return left;
		}
	case fol_formula_type::IF_THEN:
		left = scope_to_formula<fol_formula_type::AND>(scope.noncommutative.left, scope.noncommutative.left_negated);
		if (left == NULL) return NULL;
		right = scope_to_formula<fol_formula_type::AND>(scope.noncommutative.right, scope.noncommutative.right_negated);
		if (right == NULL) {
			free(*left); if (left->reference_count == 0) free(left);
			return NULL;
		}

		new_formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (new_formula == NULL) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			free(*left); if (left->reference_count == 0) free(left);
			free(*right); if (right->reference_count == 0) free(right);
			return NULL;
		}
		new_formula->type = fol_formula_type::IF_THEN;
		new_formula->reference_count = 1;
		new_formula->binary.left = left;
		new_formula->binary.right = right;
		return new_formula;
	case fol_formula_type::NOT:
		left = scope_to_formula(*scope.unary);
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

inline bool scope_contains(const<fol_scope>& scope, const fol_scope& subscope, unsigned int& index)
{
	for (index = 0; index < scope.children.length; index++) {
		auto result = compare(subscope, scope.children[index], canonicalizer());
		if (result < 0) {
			break;
		} else if (result == 0) {
			return true;
		}
	}
	return false;
}

inline bool scope_contains(const<fol_scope>& scope, const fol_scope& subscope) {
	unsigned int index;
	return scope_contains(scope, subscope, index);
}

template<fol_formula_type Operator, bool Antecedent = true>
inline void add_true_to_scope(fol_scope& scope)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IF_THEN
			   || Operator == fol_formula_type::IFF,
			"Operator is not a binary operator.");

	if (Operator == fol_formula_type::OR) {
		/* TODO: deep free `scope` */
		scope.type = fol_formula_type::TRUE;
	} else if (Operator == fol_formula_type::IF_THEN && !Antecedent) {
		/* TODO: deep free `scope` */
		scope.type = fol_formula_type::TRUE;
	}
}

/* NOTE: this function assumes `scope.commutative.children` has capacity at
		 least `scope.commutative.children.length + 1` */
template<fol_formula_type Operator, bool Antecedent = true>
inline void add_false_to_scope(fol_scope& scope)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IF_THEN
			   || Operator == fol_formula_type::IFF,
			"Operator is not a binary operator.");

	if (Operator == fol_formula_type::AND) {
		/* TODO: deep free `scope` */
		scope.type = fol_formula_type::FALSE;
	} else if (Operator == fol_formula_type::IF_THEN && Antecedent) {
		/* TODO: deep free `scope` */
		scope.type = fol_formula_type::TRUE;
	} else if (Operator == fol_formula_type::IFF) {
		if (scope.commutative.children.last().type == fol_formula_type::FALSE) {
			scope.commutative.children.length--;
		} else {
			scope.commutative.children.length++;
			scope.commutative.children.last().type = fol_formula_type::FALSE;
		}
	}
}

template<fol_formula_type Operator>
bool add_to_scope(
		const fol_scope& subscope,
		array<fol_scope>& children,
		array<fol_scope>& negated,
		bool& found_negated)
{
	unsigned int i;
	if (scope_contains(subscope, negated, i)) {
		found_negated = true;
		if (Operator == fol_formula_type::IFF) {
			free(negated[i]);
			shift_left(negated.data + i, negated.length - i - 1);
			negated.length--;
		} else {
			return true;
		}
	}
	found_negated = false;

	/* add `subscope` into the correct (sorted) position in `children` */
	if (scope_contains(subscope, children, i)) {
		/* we found an operand in `children` that is identical to `subscope` */
		if (Operator == fol_formula_type::IFF) {
			free(children[i]);
			shift_left(children.data + i, children.length - i - 1);
			children.length--;
		}
		return true;
	}

	/* `subscope` is unique, so insert it at index `i` */
	if (children.ensure_capacity(children.length + 1))
		return false;
	shift_right(children.data, children.length, i);
	move(subscope, children[i]);
	children.length++;
	return true;
}

template<fol_formula_type Operator>
bool add_to_scope(
		const fol_scope& subscope,
		array<fol_scope>& children,
		array<fol_scope>& negated,
		bool& found_negated)
{
	/* check if `subscope` is the negation of any operand in `scope.commutative.children` */
	if (subscope.type == fol_formula_type::NOT) {
		return add_to_scope<Operator>(*subscope.unary,
				negated, children, found_negated);
	} else {
		return add_to_scope<Operator>(subscope,
				children, negated, found_negated);
	}
}

template<fol_formula_type Operator>
unsigned int intersection_size(
	array<fol_scope>& first,
	array<fol_scope>& second)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int intersection_count = 0;
	unsigned int i = 0, j = 0, first_index = 0, second_index = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j], canonicalizer());
		if (result == 0) {
			if (Operator == fol_formula_type::AND) {
				success = add_false_to_scope<Operator>(scope); return 1;
			} else if (Operator == fol_formula_type::OR) {
				success = add_true_to_scope<Operator>(scope); return 1;
			} else if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	while (i < first.length) {
		move(first[i], first[first_index]);
		i++; first_index++;
	} while (j < second.length) {
		move(second[j], second[second_index]);
		j++; second_index++;
	}
	first.length = first_index;
	second.length = second_index;
	return intersection_count;
}

template<fol_formula_type Operator>
unsigned int intersection_size(
	array<fol_scope>& first,
	array<fol_scope>& second,
	unsigned int skip_second_index)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int intersection_count = 0;
	unsigned int i = 0, j = 0, first_index = 0, second_index = 0;
	while (i < first.length && j < second.length)
	{
		if (j == skip_second_index) j++;

		auto result = compare(first[i], second[j], canonicalizer());
		if (result == 0) {
			if (Operator == fol_formula_type::AND) {
				success = add_false_to_scope<Operator>(scope); return 1;
			} else if (Operator == fol_formula_type::OR) {
				success = add_true_to_scope<Operator>(scope); return 1;
			} else if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	while (i < first.length) {
		move(first[i], first[first_index]);
		i++; first_index++;
	} while (j < second.length) {
		if (j == skip_second_index) j++;
		move(second[j], second[second_index]);
		j++; second_index++;
	}
	first.length = first_index;
	second.length = second_index;
	return intersection_count;
}

template<fol_formula_type Operator>
void merge_scopes(array<fol_scope>& dst,
	const array<fol_scope>& first,
	const array<fol_scope>& second)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j], canonicalizer());
		if (result == 0) {
			if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[i]);
				i++; j++;
			} else {
				move(first[i], dst[dst.length]);
				free(second[i]);
				dst.length++; i++; j++;
			}
		} else if (result < 0) {
			move(first[i], dst[dst.length]);
			dst.length++; i++;
		} else {
			move(second[j], dst[dst.length]);
			dst.length++; j++;
		}
	}

	while (i < first.length) {
		move(first[i], dst[dst.length]);
		dst.length++; i++;
	} while (j < second.length) {
		move(second[j], dst[dst.length]);
		dst.length++; j++;
	}
}

template<fol_formula_type Operator>
void merge_scopes(array<fol_scope>& dst,
	const array<fol_scope>& first,
	const array<fol_scope>& second,
	unsigned int skip_second_index,
	unsigned int& new_second_index)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		if (j == skip_second_index) { j++; new_second_index = dst.length; }

		auto result = compare(first[i], second[j], canonicalizer());
		if (result == 0) {
			if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[i]);
				i++; j++;
			} else {
				move(first[i], dst[dst.length]);
				free(second[i]);
				dst.length++; i++; j++;
			}
		} else if (result < 0) {
			move(first[i], dst[dst.length]);
			dst.length++; i++;
		} else {
			move(second[j], dst[dst.length]);
			dst.length++; j++;
		}
	}

	while (i < first.length) {
		move(first[i], dst[dst.length]);
		dst.length++; i++;
	} while (j < second.length) {
		if (j == skip_second_index) { j++; new_second_index = dst.length; }
		move(second[j], dst[dst.length]);
		dst.length++; j++;
	}
}

template<fol_formula_type Operator>
inline void merge_scopes(
		const array<fol_scope>& src, array<fol_scope>& dst,
		const array<fol_scope>& src_negated, array<fol_scope>& dst_negated,
		bool& found_negation)
{
	unsigned int intersection_count = intersection_size(src, dst_negated);
	if (Operator == fol_formula_type::AND && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::OR && intersection_count > 0) {
		found_negation = true; return;
	}

	intersection_count += intersection_size(negated, dst);
	if (Operator == fol_formula_type::AND && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::OR && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	/* merge the two scopes */
	array<fol_scope> both = array<fol_scope>(src.length + dst.length);
	merge_scopes(both, src, dst);
	swap(both, dst); both.clear();

	merge_scopes(both, negated, dst_negated);
	swap(both, dst_negated);
}

template<fol_formula_type Operator>
inline void merge_scopes(
		const array<fol_scope>& src, array<fol_scope>& dst,
		const array<fol_scope>& src_negated, array<fol_scope>& dst_negated,
		bool& found_negation, unsigned int skip_dst_index, unsigned int& new_dst_index)
{
	unsigned int intersection_count = intersection_size(src, dst_negated);
	if (Operator == fol_formula_type::AND && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::OR && intersection_count > 0) {
		found_negation = true; return;
	}

	intersection_count += intersection_size(negated, dst, skip_dst_index);
	if (Operator == fol_formula_type::AND && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::OR && intersection_count > 0) {
		found_negation = true; return;
	} else if (Operator == fol_formula_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	/* merge the two scopes */
	array<fol_scope> both = array<fol_scope>(src.length + dst.length);
	merge_scopes(both, src, dst, skip_dst_index, new_dst_index);
	swap(both, dst); both.clear();

	merge_scopes(both, negated, dst_negated);
	swap(both, dst_negated);
}

bool negate_scope(const fol_scope& src, fol_scope& dst)
{
	if (src.type == fol_formula_type::NOT) {
		move(*src.unary, dst);
	} else if (src.type == fol_formula_type::IFF) {
		move(src, dst);
		if (!dst.commutative.children.ensure_capacity(dst.commutative.children.length + 1))
			return false;
		add_false_to_scope<fol_formula_type::IFF>(dst);
	} else {
		dst.type = fol_formula_type::NOT;
		dst.unary = (fol_scope*) malloc(sizeof(fol_scope));
		if (dst.unary == NULL) {
			fprintf(stderr, "negate_scope ERROR: Out of memory.\n");
			return false;
		}
		move(src, *dst.unary);
	}
	return true;
}

template<fol_formula_type Operator>
bool canonicalize_commutative_scope(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	if (!canonicalize(*src.left, out, variable_map)) return false;

	fol_scope& right = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (out.type == fol_formula_type::FALSE) {
		if (Operator == fol_formula_type::AND) {
			return true;
		} else if (Operator == fol_formula_type::OR) {
			return canonicalize(*src.right, out, variable_map);
		} else if (Operator == fol_formula_type::IFF) {
			if (!canonicalize(*src.right, right, variable_map)) {
				return false;
			} else if (!negate_scope(right, out)) {
				free(right);
				return false;
			}
		}
	} else if (out.type == fol_formula_type::TRUE) {
		if (Operator == fol_formula_type::OR) {
			return true;
		} else {
			return canonicalize(*src.right, right, variable_map);
		}
	}

	if (!canonicalize(*src.right, right, variable_map))
		return false;

	if (out == right) {
		/* consider the case where `*src.left` and `*src.right` are identical */
		if (Operator == fol_formula_type::IFF) {
			free(out); free(right);
			out.type = fol_formula_type::TRUE;
		} else {
			free(right);
		}
	} else if (right.type == fol_formula_type::FALSE) {
		if (Operator == fol_formula_type::AND) {
			free(out);
			out.type = fol_formula_type::FALSE;
		} else if (Operator == fol_formula_type::IFF) {
			if (!negate_scope(out, right)) {
				free(out); return false;
			}
			move(right, out);
		}
	} else if (right.type == fol_formula_type::TRUE) {
		if (Operator == fol_formula_type::OR) {
			free(out);
			out.type = fol_formula_type::TRUE;
			return true;
		}
	} else if ((out.type == fol_formula_type::NOT && *out.unary == right)
			|| (right.type == fol_formula_type::NOT && *right.unary == out)) {
		/* we have the case `A op ~A` */
		free(right); free(out);
		if (Operator == fol_formula_type::AND || Operator == fol_formula_type::IFF)
			out.type = fol_formula_type::FALSE;
		else if (Operator == fol_formula_type::OR)
			out.type = fol_formula_type::TRUE;
		return true;
	} else if (out.type == Operator) {
		if (right.type == Operator) {
			bool found_negation;
			merge_scopes<Operator>(right.commutative.children, out.commutative.children,
					right.commutative.negated, out.commutative.negated, found_negation);
			move_variables(right.variables, out.variables);
			free(right);
			if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(out); out.type = fol_formula_type::FALSE;
				} else if (Operator == fol_formula_type::OR) {
					free(out); out.type = fol_formula_type::TRUE;
				} else if (Operator == fol_formula_type::IFF) {
					if (!out.commutative.children.ensure_capacity(out.commutative.children.length + 1)) {
						free(out); return false;
					}
					add_false_to_scope<Operator>(out);
				}
			}
		} else {
			bool found_negated;
			if (!add_to_scope(right, out.commutative.children, out.commutative.negated, found_negated)) {
				free(out); free(right); return false;
			} else if (found_negated) {
				free(right);
				if (Operator == fol_formula_type::AND) {
					free(out); return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(out); return init(out, fol_formula_type::TRUE);
				} else {
					add_false_to_scope<Operator>(out);
				}
			} else {
				move_variables(right.variables, out.variables);
			}
		}
	} else {
		if (right.type == Operator) {
			bool found_negated;
			if (!add_to_scope(out, right.commutative.children, right.commutative.negated, found_negated)) {
				free(out); free(right); return false;
			} else if (found_negated) {
				free(out);
				if (Operator == fol_formula_type::AND) {
					free(right); return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(right); return init(out, fol_formula_type::TRUE);
				} else {
					add_false_to_scope<Operator>(right);
				}
			} else {
				move_variables(out.variables, right.variables);
			}
			move(right, out);
		} else {
			fol_scope both = fol_scope(Operator);

			if (less_than(out, right, canonicalizer())) {
				move(out, both.commutative.children[0]);
				move(right, both.commutative.children[1]);
			} else {
				move(right, out.commutative.children[0]);
				move(out, out.commutative.children[1]);
			}
			out.commutative.children.length = 2;
			move_variables(right.variables, out.variables);
			swap(out, both);
		}
	}
	return true;
}

bool canonicalize_conditional_scope(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope& left = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (!canonicalize(*src.left, left, variable_map)) return false;

	if (left.type == fol_formula_type::FALSE) {
		return init(out, fol_formula_type::TRUE);
	} else if (left.type == fol_formula_type::TRUE) {
		return canonicalize(*src.right, out, variable_map);
	}

	if (!canonicalize(*src.right, out, variable_map))
		return false;

	if (out == left) {
		/* consider the case where `*src.left` and `*src.right` are identical */
		free(out); free(left);
		return init(out, fol_formula_type::TRUE);
	} else if (out.type == fol_formula_type::FALSE) {
		if (!negate_scope(out, left)) {
			free(out); return false;
		}
		move(left, out);
	} else if (out.type == fol_formula_type::TRUE) {
		/* this is a no-op */
	} else if ((out.type == fol_formula_type::NOT && *out.unary == left)
			|| (left.type == fol_formula_type::NOT && *left.unary == out)) {
		/* we have the case `A => ~A` or `~A => A`, which is also a no-op */
	} else {
		/* first construct the conditional */
		if (out.type == fol_formula_type::OR) {
			fol_scope temp = fol_scope(fol_formula_type::IF_THEN);
			swap(out.commutative.children, temp.noncommutative.right);
			swap(out.commutative.negated, temp.noncommutative.right_negated);
			swap(out.variables, temp.variables);
			swap(out, temp);

			/* check if any operands in the OR can be raised into this IF_THEN consequent scope */
			array<fol_scope> to_merge = array<fol_scope>(8);
			for (unsigned int i = 0; i < out.noncommutative.right.length; i++) {
				fol_scope& child = out.noncommutative.right[i];
				if (child.type != fol_formula_type::IF_THEN) continue;

				bool found_negation;
				fol_scope& temp = *((fol_scope) alloca(sizeof(fol_scope)));
				move(child, temp);
				merge_scopes<fol_formula_type::AND>(temp.noncommutative.left, out.noncommutative.left,
						temp.noncommutative.left_negated, out.noncommutative.left_negated, found_negation);
				temp.noncommutative.left.clear();
				temp.noncommutative.left_negated.clear();
				if (found_negation) {
					free(out); free(left); free(temp);
					return init(out, fol_formula_type::TRUE);
				}

				unsigned int new_index;
				merge_scopes<fol_formula_type::OR>(temp.noncommutative.right, out.noncommutative.right,
						temp.noncommutative.right_negated, out.noncommutative.right_negated, found_negation, i, new_index);
				temp.noncommutative.right.clear();
				temp.noncommutative.right_negated.clear();
				free(temp);
				if (found_negation) {
					free(out); free(left);
					return init(out, fol_formula_type::TRUE);
				}
				i = new_index;
			}
		} else if (out.type == fol_formula_type::NOT) {
			fol_scope temp = fol_scope(fol_formula_type::IF_THEN);
			if (!temp.variables.append(out.variables.data, out.variables.length)) {
				free(left); free(out);
				return false;
			}
			move(out, temp.noncommutative.right_negated[0]);
			temp.noncommutative.right_negated.length++;
			swap(out, temp);
		} else if (out.type != fol_formula_type::IF_THEN) {
			fol_scope temp = fol_scope(fol_formula_type::IF_THEN);
			if (!temp.variables.append(out.variables.data, out.variables.length)) {
				free(left); free(out);
				return false;
			}
			move(out, temp.noncommutative.right[0]);
			temp.noncommutative.right.length++;
			swap(out, temp);
		}

		/* now try merging `left` with the conditional */
		if (left.type == fol_formula_type::AND) {
			bool found_negation;
			merge_scopes<fol_formula_type::AND>(left.commutative.children, out.noncommutative.left,
					left.commutative.negated, out.noncommutative.left_negated, found_negation);
			left.commutative.children.clear();
			left.commutative.negated.clear();
			if (found_negation) {
				free(out); free(left);
				return init(out, fol_formula_type::TRUE);
			}
			move_variables(left.variables, out.variables);
			free(left);
		} else {
			bool found_negated;
			if (!add_to_scope(left, out.noncommutative.left, out.noncommutative.left_negated, found_negated)) {
				free(out); free(left); return false;
			} else if (found_negated) {
				free(out); free(left);
				return init(out, fol_formula_type::TRUE);
			}
			move_variables(left.variables, out.variables);
		}
	}
	return true;
}

inline bool canonicalize_negation_scope(
		const fol_unary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope& temp = *((fol_scope*) malloc(sizeof(fol_scope)));
	if (!canonicalize(*src.operand, temp, variable_map)) return false;
	return negate_scope(temp, out);
}

inline bool canonicalize_term(
		const fol_term& src, fol_term& dst,
		array<unsigned int>& variables)
{
	dst.type = src.type;

	switch (src.type) {
	case fol_term_type::CONSTANT:
		dst.constant = src.constant; return true;
	case fol_term_type::VARIABLE:
		dst.variable = src.variable;
		return variables.add(dst.variable);
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "canonicalize_term ERROR: Unrecognized fol_term_type.\n");
	return false;
}

inline bool canonicalize_atom(
		const fol_atom& src, fol_scope& out)
{
	if (!init(out, fol_formula_type::ATOM)) return false;

	out.atom.predicate = src.predicate;
	if (!canonicalize_term(src.arg1, out.atom.arg1, out.variables)
	 || !canonicalize_term(src.arg2, out.atom.arg2, out.variables)) {
		return false;
	}
	if (out.variables.length > 1) {
		insertion_sort(out.variables);
		unique(out.variables);
	}

	return true;
}

bool canonicalize_scope(const fol_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	switch (src.type) {
	case fol_formula_type::AND:
		return canonicalize_commutative_scope<fol_formula_type::AND>(src.binary, out, variable_map);
	case fol_formula_type::OR:
		return canonicalize_commutative_scope<fol_formula_type::OR>(src.binary, out, variable_map);
	case fol_formula_type::IFF:
		return canonicalize_commutative_scope<fol_formula_type::IFF>(src.binary, out, variable_map);
	case fol_formula_type::IF_THEN:
		return canonicalize_conditional_scope(src.binary, out, variable_map);
	case fol_formula_type::NOT:
		return canonicalize_negation_scope(src.unary, out, variable_map);
	case fol_formula_type::ATOM:
		return canonicalize_atom(src.atom, out);
	case fol_formula_type::FOR_ALL:
		/* TODO: implement this */
	case fol_formula_type::EXISTS:
		/* TODO: implement this */
	case fol_formula_type::TRUE:
		return init(out, fol_formula_type::TRUE);
	case fol_formula_type::FALSE:
		return init(out, fol_formula_type::FALSE);
	}
	fprintf(stderr, "canonicalize_scope ERROR: Unrecognized fol_formula_type.\n");
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
