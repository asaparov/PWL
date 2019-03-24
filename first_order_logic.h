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
bool operator != (const fol_formula&, const fol_formula&);

enum class fol_term_type {
	NONE = 0,
	VARIABLE,
	CONSTANT,
	PARAMETER
};

struct fol_term {
	fol_term_type type;
	unsigned int reference_count;
	union {
		unsigned int variable;
		unsigned int constant;
		unsigned int parameter;
	};

	static inline bool is_empty(const fol_term& key) {
		return key.type == fol_term_type::NONE;
	}

	static inline void set_empty(fol_term& key) {
		key.type = fol_term_type::NONE;
	}

	static inline bool copy(const fol_term& src, fol_term& dst) {
		dst = src;
		return true;
	}

	static inline void move(const fol_term& src, fol_term& dst) {
		dst = src;
	}

	static inline void swap(fol_term& first, fol_term& second) {
		fol_term temp = first;
		first = second;
		second = temp;
	}

	static inline unsigned int hash(const fol_term& key) {
		switch (key.type) {
		case fol_term_type::VARIABLE:
			return 0 + 3 * default_hash(key.variable);
		case fol_term_type::CONSTANT:
			return 1 + 3 * default_hash(key.constant);
		case fol_term_type::PARAMETER:
			return 2 + 3 * default_hash(key.parameter);
		case fol_term_type::NONE:
			break;
		}
		fprintf(stderr, "fol_term.hash ERROR: Unrecognized fol_term_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void free(fol_term& term) { }
};

inline bool operator == (const fol_term& first, const fol_term& second) {
	if (first.type != second.type) return false;
	switch (first.type) {
	case fol_term_type::VARIABLE:
		return first.variable == second.variable;
	case fol_term_type::CONSTANT:
		return first.constant == second.constant;
	case fol_term_type::PARAMETER:
		return first.parameter == second.parameter;
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline bool operator != (const fol_term& first, const fol_term& second) {
	return !(first == second);
}

inline bool operator < (const fol_term& first, const fol_term& second) {
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case fol_term_type::VARIABLE:
		return first.variable < second.variable;
	case fol_term_type::CONSTANT:
		return first.constant < second.constant;
	case fol_term_type::PARAMETER:
		return first.parameter < second.parameter;
	case fol_term_type::NONE:
		return false;
	}
	fprintf(stderr, "operator < ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

template<typename Stream>
inline bool print_subscript(unsigned int number, Stream& out) {
	static const char* subscripts[] = { "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉" };
	while (number > 0) {
		if (!print(subscripts[number % 10], out)) return false;
		number /= 10;
	}
	return true;
}

enum class fol_formula_syntax {
	TPTP,
	CLASSIC
};

template<fol_formula_syntax Syntax, typename Stream>
inline bool print_variable(unsigned int variable, Stream& out) {
	switch (Syntax) {
	case fol_formula_syntax::TPTP:
		return print('$', out) + print(variable, out);
	case fol_formula_syntax::CLASSIC:
		return print('x', out) && print_subscript(variable, out);
	}
	fprintf(stderr, "print_variable ERROR: Unrecognized fol_formula_syntax.\n");
	return false;
}

template<fol_formula_syntax Syntax, typename Stream>
inline bool print_parameter(unsigned int parameter, Stream& out) {
	switch (Syntax) {
	case fol_formula_syntax::TPTP:
		return print('#', out) + print(parameter, out);
	case fol_formula_syntax::CLASSIC:
		return print('a', out) && print_subscript(parameter, out);
	}
	fprintf(stderr, "print_parameter ERROR: Unrecognized fol_formula_syntax.\n");
	return false;
}

template<fol_formula_syntax Syntax = fol_formula_syntax::CLASSIC, typename Stream, typename... Printer>
bool print(const fol_term& term, Stream& out, Printer&&... printer) {
	switch (term.type) {
	case fol_term_type::VARIABLE:
		return print_variable<Syntax>(term.variable, out);
	case fol_term_type::CONSTANT:
		return print(term.constant, out, std::forward<Printer>(printer)...);
	case fol_term_type::PARAMETER:
		return print_parameter<Syntax>(term.parameter, out);
	case fol_term_type::NONE:
		break;
	}

	fprintf(stderr, "print ERROR: Unexpected fol_term_type.\n");
	return false;
}

enum class fol_formula_type {
	ATOM,
	EQUALS,

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

struct fol_equals {
	fol_term arg1;
	fol_term arg2;

	static inline void move(const fol_equals& src, fol_equals& dst) {
		dst.arg1 = src.arg1;
		dst.arg2 = src.arg2;
	}

	static inline void free(fol_equals& equals) { }
};

inline bool operator == (const fol_equals& first, const fol_equals& second) {
	return first.arg1 == second.arg1
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

struct fol_array_formula {
	fol_formula** operands;
	unsigned int length;

	static inline void move(const fol_array_formula& src, fol_array_formula& dst);
	static inline void free(fol_array_formula& formula);
};

struct fol_quantifier {
	unsigned int variable;
	fol_formula* operand;

	static inline void move(const fol_quantifier& src, fol_quantifier& dst);
	static inline void free(fol_quantifier& formula);
};

struct fol_formula
{
	typedef fol_formula_type Type;
	typedef fol_term Term;
	typedef fol_term_type TermType;

	static constexpr fol_term EMPTY_TERM = {fol_term_type::NONE, 1, {0}};

	fol_formula_type type;
	unsigned int reference_count;
	union {
		fol_atom atom;
		fol_equals equals;
		fol_unary_formula unary;
		fol_binary_formula binary;
		fol_array_formula array;
		fol_quantifier quantifier;
	};

	fol_formula(fol_formula_type type) : type(type), reference_count(1) { }
	~fol_formula() { free_helper(); }

	static fol_term* new_variable(unsigned int variable);
	static fol_term* new_constant(unsigned int constant);
	static fol_term* new_parameter(unsigned int parameter);
	static fol_formula* new_atom(unsigned int predicate, fol_term* arg1, fol_term* arg2);
	static inline fol_formula* new_atom(unsigned int predicate, fol_term* arg1);
	static inline fol_formula* new_atom(unsigned int predicate);
	static inline fol_formula* new_equals(fol_term arg1, fol_term arg2);
	static fol_formula* new_true();
	static fol_formula* new_false();
	template<typename... Args> static inline fol_formula* new_and(Args&&... args);
	template<typename... Args> static inline fol_formula* new_or(Args&&... args);
	template<typename... Args> static inline fol_formula* new_iff(Args&&... args);
	static fol_formula* new_if_then(fol_formula* first, fol_formula* second);
	static fol_formula* new_not(fol_formula* operand);
	static inline fol_formula* new_for_all(unsigned int variable, fol_formula* operand);
	static inline fol_formula* new_exists(unsigned int variable, fol_formula* operand);

	static inline void move(const fol_formula& src, fol_formula& dst);
	static inline void free(fol_formula& formula) { formula.free_helper(); }

private:
	void free_helper();
};

constexpr fol_term fol_formula::EMPTY_TERM;

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

inline bool operator == (const fol_array_formula& first, const fol_array_formula& second) {
	if (first.length != second.length) return false;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.operands[i] != second.operands[i] && *first.operands[i] != *second.operands[i]) return false;
	return true;
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
	case fol_formula_type::EQUALS:
		return first.equals == second.equals;
	case fol_formula_type::NOT:
		return first.unary == second.unary;
	case fol_formula_type::IF_THEN:
		return first.binary == second.binary;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return first.array == second.array;
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
	case fol_formula_type::EQUALS:
		core::move(src.equals, dst.equals); return;
	case fol_formula_type::NOT:
		core::move(src.unary, dst.unary); return;
	case fol_formula_type::IF_THEN:
		core::move(src.binary, dst.binary); return;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		core::move(src.array, dst.array); return;
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

inline void fol_array_formula::move(const fol_array_formula& src, fol_array_formula& dst) {
	dst.length = src.length;
	dst.operands = src.operands;
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
		case fol_formula_type::EQUALS:
			core::free(equals); return;
		case fol_formula_type::NOT:
			core::free(unary); return;
		case fol_formula_type::IF_THEN:
			core::free(binary); return;
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			core::free(array); return;
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

inline void fol_array_formula::free(fol_array_formula& formula) {
	for (unsigned int i = 0; i < formula.length; i++) {
		core::free(*formula.operands[i]);
		if (formula.operands[i]->reference_count == 0)
			core::free(formula.operands[i]);
	}
	core::free(formula.operands);
}

inline void fol_quantifier::free(fol_quantifier& formula) {
	core::free(*formula.operand);
	if (formula.operand->reference_count == 0)
		core::free(formula.operand);
}

inline bool is_atomic(
		const fol_formula& term, unsigned int& predicate,
		fol_term const*& arg1, fol_term const*& arg2)
{
	if (term.type != fol_formula_type::ATOM) return false;
	predicate = term.atom.predicate;
	arg1 = (term.atom.arg1.type == fol_term_type::NONE ? NULL : &term.atom.arg1);
	arg2 = (term.atom.arg2.type == fol_term_type::NONE ? NULL : &term.atom.arg2);
	return true;
}

inline bool is_atomic(const fol_formula& term,
		fol_term const*& arg1, fol_term const*& arg2)
{
	unsigned int predicate;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(const fol_formula& term) {
	unsigned int predicate;
	fol_term const* arg1; fol_term const* arg2;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(
		fol_formula& term, unsigned int& predicate,
		fol_term*& arg1, fol_term*& arg2)
{
	if (term.type != fol_formula_type::ATOM) return false;
	predicate = term.atom.predicate;
	arg1 = (term.atom.arg1.type == fol_term_type::NONE ? NULL : &term.atom.arg1);
	arg2 = (term.atom.arg2.type == fol_term_type::NONE ? NULL : &term.atom.arg2);
	return true;
}

inline bool is_atomic(fol_formula& term,
		fol_term*& arg1, fol_term*& arg2)
{
	unsigned int predicate;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(fol_formula& term) {
	unsigned int predicate;
	fol_term* arg1; fol_term* arg2;
	return is_atomic(term, predicate, arg1, arg2);
}

template<fol_formula_syntax Syntax> struct and_symbol;
template<fol_formula_syntax Syntax> struct or_symbol;
template<fol_formula_syntax Syntax> struct if_then_symbol;
template<fol_formula_syntax Syntax> struct iff_symbol;
template<fol_formula_syntax Syntax> struct not_symbol;
template<fol_formula_syntax Syntax> struct equals_symbol;
template<fol_formula_syntax Syntax> struct true_symbol;
template<fol_formula_syntax Syntax> struct false_symbol;

template<> struct and_symbol<fol_formula_syntax::TPTP> { static const char symbol[]; };
template<> struct or_symbol<fol_formula_syntax::TPTP> { static const char symbol[]; };
template<> struct iff_symbol<fol_formula_syntax::TPTP> { static const char symbol[]; };
template<> struct if_then_symbol<fol_formula_syntax::TPTP> { static const char symbol[]; };
template<> struct not_symbol<fol_formula_syntax::TPTP> { static const char symbol; };
template<> struct equals_symbol<fol_formula_syntax::TPTP> { static const char symbol; };
template<> struct true_symbol<fol_formula_syntax::TPTP> { static const char symbol; };
template<> struct false_symbol<fol_formula_syntax::TPTP> { static const char symbol; };

template<> struct and_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct or_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct iff_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct if_then_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct not_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct equals_symbol<fol_formula_syntax::CLASSIC> { static const char symbol; };
template<> struct true_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };
template<> struct false_symbol<fol_formula_syntax::CLASSIC> { static const char symbol[]; };

const char and_symbol<fol_formula_syntax::TPTP>::symbol[] = " & ";
const char or_symbol<fol_formula_syntax::TPTP>::symbol[] = " | ";
const char iff_symbol<fol_formula_syntax::TPTP>::symbol[] = " <=> ";
const char if_then_symbol<fol_formula_syntax::TPTP>::symbol[] = " => ";
const char not_symbol<fol_formula_syntax::TPTP>::symbol = '~';
const char equals_symbol<fol_formula_syntax::TPTP>::symbol = '=';
const char true_symbol<fol_formula_syntax::TPTP>::symbol = 'T';
const char false_symbol<fol_formula_syntax::TPTP>::symbol = 'F';

const char and_symbol<fol_formula_syntax::CLASSIC>::symbol[] = " ∧ ";
const char or_symbol<fol_formula_syntax::CLASSIC>::symbol[] = " ∨ ";
const char iff_symbol<fol_formula_syntax::CLASSIC>::symbol[] = " ↔ ";
const char if_then_symbol<fol_formula_syntax::CLASSIC>::symbol[] = " → ";
const char not_symbol<fol_formula_syntax::CLASSIC>::symbol[] = "¬";
const char equals_symbol<fol_formula_syntax::CLASSIC>::symbol = '=';
const char true_symbol<fol_formula_syntax::CLASSIC>::symbol[] = "⊤";
const char false_symbol<fol_formula_syntax::CLASSIC>::symbol[] = "⊥";

const char left_parens[] = "(";
const char right_parens[] = ")";

template<fol_formula_syntax Syntax, typename Stream>
inline bool print_for_all(unsigned int quantified_variable, Stream& out) {
	switch (Syntax) {
	case fol_formula_syntax::TPTP:
		return print("![", out) && print_variable<Syntax>(quantified_variable, out) && print("]:", out);
	case fol_formula_syntax::CLASSIC:
		return print("∀", out) && print_variable<Syntax>(quantified_variable, out);
	}
	fprintf(stderr, "print_for_all ERROR: Unrecognized fol_formula_syntax.\n");
	return false;
}

template<fol_formula_syntax Syntax, typename Stream>
inline bool print_exists(unsigned int quantified_variable, Stream& out) {
	switch (Syntax) {
	case fol_formula_syntax::TPTP:
		return print("?[", out) && print_variable<Syntax>(quantified_variable, out) && print("]:", out);
	case fol_formula_syntax::CLASSIC:
		return print("∃", out) && print_variable<Syntax>(quantified_variable, out);
	}
	fprintf(stderr, "print_exists ERROR: Unrecognized fol_formula_syntax.\n");
	return false;
}

template<fol_formula_syntax Syntax = fol_formula_syntax::CLASSIC, typename Stream, typename... Printer>
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

	case fol_formula_type::EQUALS:
		return print(formula.equals.arg1, out, std::forward<Printer>(printer)...)
			&& print(equals_symbol<Syntax>::symbol, out)
			&& print(formula.equals.arg2, out, std::forward<Printer>(printer)...);

	case fol_formula_type::TRUE:
		return print(true_symbol<Syntax>::symbol, out);

	case fol_formula_type::FALSE:
		return print(false_symbol<Syntax>::symbol, out);

	case fol_formula_type::NOT:
		return print(not_symbol<Syntax>::symbol, out) && print(*formula.unary.operand, out, std::forward<Printer>(printer)...);

	case fol_formula_type::AND:
		return print<fol_formula*, left_parens, right_parens, and_symbol<Syntax>::symbol>(formula.array.operands, formula.array.length, out, pointer_scribe(), std::forward<Printer>(printer)...);

	case fol_formula_type::OR:
		return print<fol_formula*, left_parens, right_parens, or_symbol<Syntax>::symbol>(formula.array.operands, formula.array.length, out, pointer_scribe(), std::forward<Printer>(printer)...);

	case fol_formula_type::IF_THEN:
		return print('(', out) && print(*formula.binary.left, out, std::forward<Printer>(printer)...)
			&& print(if_then_symbol<Syntax>::symbol, out) && print(*formula.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case fol_formula_type::IFF:
		return print<fol_formula*, left_parens, right_parens, iff_symbol<Syntax>::symbol>(formula.array.operands, formula.array.length, out, pointer_scribe(), std::forward<Printer>(printer)...);

	case fol_formula_type::FOR_ALL:
		return print_for_all<Syntax>(formula.quantifier.variable, out)
			&& print(*formula.quantifier.operand, out, std::forward<Printer>(printer)...);

	case fol_formula_type::EXISTS:
		return print_exists<Syntax>(formula.quantifier.variable, out)
			&& print(*formula.quantifier.operand, out, std::forward<Printer>(printer)...);
	}

	fprintf(stderr, "print ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

inline bool new_fol_formula(fol_formula*& new_formula) {
	new_formula = (fol_formula*) malloc(sizeof(fol_formula));
	if (new_formula == NULL) {
		fprintf(stderr, "new_fol_formula ERROR: Out of memory.\n");
		return false;
	}
	return true;
}

constexpr bool visit_constant(unsigned int constant) { return true; }
constexpr bool visit_equals(const fol_formula& formula) { return true; }
constexpr bool visit_variable(unsigned int variable) { return true; }
constexpr bool visit_parameter(unsigned int parameter) { return true; }
constexpr bool visit_true(const fol_formula& formula) { return true; }
constexpr bool visit_false(const fol_formula& formula) { return true; }

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula& formula) { return true; }

template<typename Term, typename... Visitor,
	typename std::enable_if<std::is_same<typename std::remove_cv<typename std::remove_reference<Term>::type>::type, fol_term>::value>::type* = nullptr>
inline bool visit(Term&& term, Visitor&&... visitor) {
	switch (term.type) {
	case fol_term_type::CONSTANT:
		return visit_constant(term.constant, std::forward<Visitor>(visitor)...);
	case fol_term_type::VARIABLE:
		return visit_variable(term.variable, std::forward<Visitor>(visitor)...);
	case fol_term_type::PARAMETER:
		return visit_parameter(term.parameter, std::forward<Visitor>(visitor)...);
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "visit ERROR: Unrecognized fol_term_type.\n");
	return false;
}

template<typename Atom, typename... Visitor,
	typename std::enable_if<std::is_same<typename std::remove_cv<typename std::remove_reference<Atom>::type>::type, fol_atom>::value>::type* = nullptr>
inline bool visit(Atom&& atom, Visitor&&... visitor) {
	return visit_constant(atom.predicate, std::forward<Visitor>(visitor)...)
		&& visit(atom.arg1, std::forward<Visitor>(visitor)...)
		&& visit(atom.arg2, std::forward<Visitor>(visitor)...);
}

template<typename Formula, typename... Visitor,
	typename std::enable_if<std::is_same<typename std::remove_cv<typename std::remove_reference<Formula>::type>::type, fol_formula>::value>::type* = nullptr>
bool visit(Formula&& formula, Visitor&&... visitor)
{
	switch (formula.type) {
	case fol_formula_type::ATOM:
		return visit(formula.atom, std::forward<Visitor>(visitor)...);
	case fol_formula_type::EQUALS:
		return visit_equals(formula, std::forward<Visitor>(visitor)...)
			&& visit(formula.equals.arg1, std::forward<Visitor>(visitor)...)
			&& visit(formula.equals.arg2, std::forward<Visitor>(visitor)...);
	case fol_formula_type::NOT:
		return visit_operator<fol_formula_type::NOT>(formula, std::forward<Visitor>(visitor)...)
			&& visit(*formula.unary.operand, std::forward<Visitor>(visitor)...);
	case fol_formula_type::AND:
		if (!visit_operator<fol_formula_type::AND>(formula, std::forward<Visitor>(visitor)...)) return false;
		for (unsigned int i = 0; i < formula.array.length; i++)
			if (!visit(*formula.array.operands[i], std::forward<Visitor>(visitor)...)) return false;
		return true;
	case fol_formula_type::OR:
		if (!visit_operator<fol_formula_type::OR>(formula, std::forward<Visitor>(visitor)...)) return false;
		for (unsigned int i = 0; i < formula.array.length; i++)
			if (!visit(*formula.array.operands[i], std::forward<Visitor>(visitor)...)) return false;
		return true;
	case fol_formula_type::IF_THEN:
		return visit_operator<fol_formula_type::IF_THEN>(formula, std::forward<Visitor>(visitor)...)
			&& visit(*formula.binary.left, std::forward<Visitor>(visitor)...)
			&& visit(*formula.binary.right, std::forward<Visitor>(visitor)...);
	case fol_formula_type::IFF:
		if (!visit_operator<fol_formula_type::IFF>(formula, std::forward<Visitor>(visitor)...)) return false;
		for (unsigned int i = 0; i < formula.array.length; i++)
			if (!visit(*formula.array.operands[i], std::forward<Visitor>(visitor)...)) return false;
		return true;
	case fol_formula_type::FOR_ALL:
		return visit_operator<fol_formula_type::FOR_ALL>(formula, std::forward<Visitor>(visitor)...)
			&& visit_variable(formula.quantifier.variable, std::forward<Visitor>(visitor)...)
			&& visit(*formula.quantifier.operand, std::forward<Visitor>(visitor)...);
	case fol_formula_type::EXISTS:
		return visit_operator<fol_formula_type::EXISTS>(formula, std::forward<Visitor>(visitor)...)
			&& visit_variable(formula.quantifier.variable, std::forward<Visitor>(visitor)...)
			&& visit(*formula.quantifier.operand, std::forward<Visitor>(visitor)...);
	case fol_formula_type::TRUE:
		return visit_true(formula, std::forward<Visitor>(visitor)...);
	case fol_formula_type::FALSE:
		return visit_false(formula, std::forward<Visitor>(visitor)...);
	}
	fprintf(stderr, "visit ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

struct parameter_comparator {
	unsigned int parameter;
};

constexpr bool visit_constant(unsigned int constant, const parameter_comparator& visitor) { return true; }
constexpr bool visit_variable(unsigned int variable, const parameter_comparator& visitor) { return true; }
constexpr bool visit_equals(const fol_formula& formula, const parameter_comparator& visitor) { return true; }
constexpr bool visit_true(const fol_formula& formula, const parameter_comparator& visitor) { return true; }
constexpr bool visit_false(const fol_formula& formula, const parameter_comparator& visitor) { return true; }

inline bool visit_parameter(unsigned int parameter, const parameter_comparator& visitor) {
	return visitor.parameter == parameter;
}

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula& formula, const parameter_comparator& visitor) { return true; }

inline bool contains_parameter(const fol_formula& src, unsigned int parameter) {
	parameter_comparator visitor = {parameter};
	return !visit(src, visitor);
}

struct parameter_collector {
	array<unsigned int>& parameters;
};

constexpr bool visit_constant(unsigned int constant, const parameter_collector& visitor) { return true; }
constexpr bool visit_variable(unsigned int variable, const parameter_collector& visitor) { return true; }
constexpr bool visit_equals(const fol_formula& formula, const parameter_collector& visitor) { return true; }
constexpr bool visit_true(const fol_formula& formula, const parameter_collector& visitor) { return true; }
constexpr bool visit_false(const fol_formula& formula, const parameter_collector& visitor) { return true; }

inline bool visit_parameter(unsigned int parameter, const parameter_collector& visitor) {
	return visitor.parameters.add(parameter);
}

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula& formula, const parameter_collector& visitor) { return true; }

inline bool get_parameters(const fol_formula& src, array<unsigned int>& parameters) {
	parameter_collector visitor = {parameters};
	return !visit(src, visitor);
}

inline bool clone_constant(unsigned int src_constant, unsigned int& dst_constant) {
	dst_constant = src_constant;
	return true;
}

inline bool clone_predicate(unsigned int src_predicate, unsigned int& dst_predicate) {
	dst_predicate = src_predicate;
	return true;
}

inline bool clone_variable(unsigned int src_variable, unsigned int& dst_variable) {
	dst_variable = src_variable;
	return true;
}

inline bool clone_parameter(unsigned int src_parameter, unsigned int& dst_parameter) {
	dst_parameter = src_parameter;
	return true;
}

template<typename... Cloner>
inline bool clone(const fol_term& src, fol_term& dst, Cloner&&... cloner) {
	dst.type = src.type;
	switch (src.type) {
	case fol_term_type::CONSTANT:
		return clone_constant(src.constant, dst.constant, std::forward<Cloner>(cloner)...);
	case fol_term_type::VARIABLE:
		return clone_variable(src.variable, dst.variable, std::forward<Cloner>(cloner)...);
	case fol_term_type::PARAMETER:
		return clone_parameter(src.parameter, dst.parameter, std::forward<Cloner>(cloner)...);
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "clone ERROR: Unrecognized fol_term_type.\n");
	return false;
}

template<typename... Cloner>
inline bool clone(const fol_atom& src, fol_atom& dst, Cloner&&... cloner) {
	return clone_predicate(src.predicate, dst.predicate, std::forward<Cloner>(cloner)...)
		&& clone(src.arg1, dst.arg1, std::forward<Cloner>(cloner)...)
		&& clone(src.arg2, dst.arg2, std::forward<Cloner>(cloner)...);
}

template<typename... Cloner>
bool clone(const fol_formula& src, fol_formula& dst, Cloner&&... cloner)
{
	dst.type = src.type;
	dst.reference_count = 1;
	switch (src.type) {
	case fol_formula_type::ATOM:
		return clone(src.atom, dst.atom, std::forward<Cloner>(cloner)...);
	case fol_formula_type::EQUALS:
		return clone(src.equals.arg1, dst.equals.arg1, std::forward<Cloner>(cloner)...)
			&& clone(src.equals.arg2, dst.equals.arg2, std::forward<Cloner>(cloner)...);
	case fol_formula_type::NOT:
		if (!new_fol_formula(dst.unary.operand)) return false;
		if (!clone(*src.unary.operand, *dst.unary.operand, std::forward<Cloner>(cloner)...)) {
			free(dst.unary.operand);
			return false;
		}
		return true;
	case fol_formula_type::IF_THEN:
		if (!new_fol_formula(dst.binary.left)) {
			return false;
		} else if (!new_fol_formula(dst.binary.right)) {
			free(dst.binary.left); return false;
		} else if (!clone(*src.binary.left, *dst.binary.left, std::forward<Cloner>(cloner)...)) {
			free(dst.binary.left); free(dst.binary.right);
			return false;
		} else if (!clone(*src.binary.right, *dst.binary.right, std::forward<Cloner>(cloner)...)) {
			free(*dst.binary.left); free(dst.binary.left);
			free(dst.binary.right); return false;
		}
		return true;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		dst.array.operands = (fol_formula**) malloc(sizeof(fol_formula*) * src.array.length);
		if (dst.array.operands == NULL) return false;
		for (unsigned int i = 0; i < src.array.length; i++) {
			if (!new_fol_formula(dst.array.operands[i])
			 || !clone(*src.array.operands[i], *dst.array.operands[i], std::forward<Cloner>(cloner)...)) {
				for (unsigned int j = 0; j < i; j++) {
					free(*dst.array.operands[j]); free(dst.array.operands[j]);
				}
				if (dst.array.operands[i] != NULL) free(dst.array.operands[i]);
				free(dst.array.operands); return false;
			}
		}
		dst.array.length = src.array.length;
		return true;
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		if (!clone_variable(src.quantifier.variable, dst.quantifier.variable, std::forward<Cloner>(cloner)...)
		 || !new_fol_formula(dst.quantifier.operand))
			return false;
		if (!clone(*src.quantifier.operand, *dst.quantifier.operand, std::forward<Cloner>(cloner)...)) {
			free(dst.quantifier.operand); return false;
		}
		return true;
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return true;
	}
	fprintf(stderr, "clone ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

template<typename... Cloner>
inline bool clone(const fol_formula* src, fol_formula* dst, Cloner&&... cloner)
{
	if (!new_fol_formula(dst)) return false;
	return clone(*src, *dst, std::forward<Cloner>(cloner)...);
}

template<typename... Function>
inline bool apply_to_terms(const fol_atom& src, fol_atom& dst, Function&&... function)
{
	dst.predicate = src.predicate;
	return apply(src.arg1, dst.arg1, std::forward<Function>(function)...)
		&& apply(src.arg2, dst.arg2, std::forward<Function>(function)...);
}

template<typename... Function>
inline bool apply_to_terms(const fol_equals& src, fol_equals& dst, Function&&... function)
{
	return apply(src.arg1, dst.arg1, std::forward<Function>(function)...)
		&& apply(src.arg2, dst.arg2, std::forward<Function>(function)...);
}

template<typename... Function>
fol_formula* apply_to_terms(fol_formula& src, Function&&... function)
{
	fol_formula* new_formula;
	fol_formula** new_formulas;
	fol_formula* left; fol_formula* right;
	fol_atom atom; fol_equals equals;
	bool changed;
	switch (src.type) {
	case fol_formula_type::ATOM:
		if (!apply_to_terms(src.atom, atom, std::forward<Function>(function)...)) {
			return NULL;
		} else if (atom == src.atom) {
			free(atom);
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				free(atom); return NULL;
			}
			new_formula->atom = atom;
			new_formula->type = fol_formula_type::ATOM;
			new_formula->reference_count = 1;
			return new_formula;
		}
	case fol_formula_type::EQUALS:
		if (!apply_to_terms(src.equals, equals, std::forward<Function>(function)...)) {
			return NULL;
		} else if (equals == src.equals) {
			free(equals);
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				free(equals); return NULL;
			}
			new_formula->equals = equals;
			new_formula->type = fol_formula_type::EQUALS;
			new_formula->reference_count = 1;
			return new_formula;
		}

	case fol_formula_type::NOT:
		left = apply_to_terms(*src.unary.operand, std::forward<Function>(function)...);
		if (left == NULL) {
			return NULL;
		} else if (left == src.unary.operand) {
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			}
			new_formula->unary.operand = left;
			new_formula->type = fol_formula_type::NOT;
			new_formula->reference_count = 1;
			return new_formula;
		}
	case fol_formula_type::IF_THEN:
		left = apply_to_terms(*src.binary.left, std::forward<Function>(function)...);
		if (left == NULL) return NULL;
		right = apply_to_terms(*src.binary.right, std::forward<Function>(function)...);
		if (right == NULL) {
			if (left != src.binary.left) {
				free(*left); if (left->reference_count == 0) free(left);
			}
			return NULL;
		} else if (left == src.binary.left && right == src.binary.right) {
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				if (left != src.binary.left) {
					free(*left); if (left->reference_count == 0) free(left);
				} if (right != src.binary.right) {
					free(*right); if (right->reference_count == 0) free(right);
				}
				return NULL;
			}
			new_formula->binary.left = left;
			new_formula->binary.right = right;
			if (left == src.binary.left) left->reference_count++;
			if (right == src.binary.right) right->reference_count++;
			new_formula->type = src.type;
			new_formula->reference_count = 1;
			return new_formula;
		}
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		new_formulas = (fol_formula**) malloc(sizeof(fol_formula*) * src.array.length);
		if (new_formulas == NULL) return NULL;
		changed = false;
		for (unsigned int i = 0; i < src.array.length; i++) {
			new_formulas[i] = apply_to_terms(*src.array.operands[i], std::forward<Function>(function)...);
			if (new_formulas[i] == NULL) {
				for (unsigned int j = 0; j < i; j++) {
					if (new_formulas[j] != src.array.operands[j]) {
						free(*new_formulas[j]); if (new_formulas[j]->reference_count == 0) free(new_formulas[j]);
					}
				}
				free(new_formulas); return NULL;
			} else if (new_formulas[i] != src.array.operands[i])
				changed = true;
		}

		if (!changed) {
			free(new_formulas);
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				for (unsigned int j = 0; j < src.array.length; j++) {
					if (new_formulas[j] != src.array.operands[j]) {
						free(*new_formulas[j]); if (new_formulas[j]->reference_count == 0) free(new_formulas[j]);
					}
				}
				free(new_formulas); return NULL;
			}
			new_formula->array.operands = new_formulas;
			new_formula->array.length = src.array.length;
			for (unsigned int i = 0; i < src.array.length; i++)
				if (new_formula->array.operands[i] == src.array.operands[i]) src.array.operands[i]->reference_count++;
			new_formula->type = src.type;
			new_formula->reference_count = 1;
			return new_formula;
		}
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		left = apply_to_terms(*src.quantifier.operand, std::forward<Function>(function)...);
		if (left == NULL) {
			return NULL;
		} else if (left == src.quantifier.operand) {
			return &src;
		} else {
			if (!new_fol_formula(new_formula)) {
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			}
			new_formula->quantifier.variable = src.quantifier.variable;
			new_formula->quantifier.operand = left;
			new_formula->type = src.type;
			new_formula->reference_count = 1;
			return new_formula;
		}
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return &src;
	}
	fprintf(stderr, "apply_to_terms ERROR: Unrecognized fol_formula_type.\n");
	return NULL;
}

template<int VariableShift>
struct term_substituter {
	fol_term src;
	fol_term dst;
};

template<int VariableShift>
inline bool apply(const fol_term& src, fol_term& dst, const term_substituter<VariableShift>& substituter) {
	if (src == substituter.src) {
		dst = substituter.dst;
	} else if (src.type == fol_term_type::VARIABLE) {
		dst.type = src.type;
		dst.variable = src.variable + VariableShift;
	} else {
		dst = src;
	}
	return true;
}

template<fol_term_type SrcTermType, int VariableShift = 0>
inline fol_formula* substitute(fol_formula* src,
		const fol_term* src_term, const fol_term* dst_term)
{
	const term_substituter<VariableShift> substituter = {*src_term, *dst_term};
	fol_formula* formula = apply_to_terms(*src, substituter);
	if (formula == src)
		formula->reference_count++;
	return formula;
}

template<int VariableShift>
struct index_substituter {
	fol_term src;
	fol_term dst;
	const unsigned int* term_indices;
	unsigned int term_index_count;
	unsigned int current_term_index;
};

template<int VariableShift>
inline bool apply(const fol_term& src, fol_term& dst, index_substituter<VariableShift>& substituter)
{
	if (substituter.term_index_count > 0 && *substituter.term_indices == substituter.current_term_index) {
		if (substituter.src.type == fol_term_type::NONE) {
			substituter.src = src;
		} else if (substituter.src != src) {
			/* this term is not identical to other substituted terms, which should not happen */
			return false;
		}
		dst = substituter.dst;
		substituter.term_indices++;
		substituter.term_index_count--;
	}
	substituter.current_term_index++;
	return true;
}

template<int VariableShift>
inline fol_formula* substitute(
		fol_formula* src, const unsigned int* term_indices,
		unsigned int term_index_count, const fol_term* dst_term)
{
	index_substituter<VariableShift> substituter = {fol_formula::EMPTY_TERM, *dst_term, term_indices, term_index_count, 0};
	fol_formula* formula = apply_to_terms(*src, substituter);
	if (formula == src)
		formula->reference_count++;
	return formula;
}

inline bool unify(
		const fol_term& first, const fol_term& second,
		const fol_term* src_term, fol_term const*& dst_term)
{
	if (first == *src_term) {
		if (dst_term == NULL) {
			dst_term = &second;
		} else if (second != *dst_term) {
			return false;
		}
	}
	return true;
}

inline bool unify(
		const fol_atom& first, const fol_atom& second,
		const fol_term* src_term, fol_term const*& dst_term)
{
	if (first.predicate != second.predicate) return false;
	return unify(first.arg1, second.arg1, src_term, dst_term)
		&& unify(first.arg2, second.arg2, src_term, dst_term);
}

inline bool unify(
		const fol_equals& first, const fol_equals& second,
		const fol_term* src_term, fol_term const*& dst_term)
{
	return unify(first.arg1, second.arg1, src_term, dst_term)
		&& unify(first.arg2, second.arg2, src_term, dst_term);
}

bool unify(
		const fol_formula& first, const fol_formula& second,
		const fol_term* src_term, fol_term const*& dst_term)
{
	if (first.type != second.type) return false;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return unify(first.atom, second.atom, src_term, dst_term);
	case fol_formula_type::EQUALS:
		return unify(first.equals, second.equals, src_term, dst_term);
	case fol_formula_type::NOT:
		return unify(*first.unary.operand, *second.unary.operand, src_term, dst_term);
	case fol_formula_type::IF_THEN:
		return unify(*first.binary.left, *second.binary.left, src_term, dst_term)
			&& unify(*first.binary.right, *second.binary.right, src_term, dst_term);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		if (first.array.length != second.array.length) return false;
		for (unsigned int i = 0; i < first.array.length; i++)
			if (!unify(*first.array.operands[i], *second.array.operands[i], src_term, dst_term)) return false;
		return true;
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		if (first.quantifier.variable != second.quantifier.variable) return false;
		return unify(*first.quantifier.operand, *second.quantifier.operand, src_term, dst_term);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return true;
	}
	fprintf(stderr, "unify ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

bool unifies_parameter(
		const fol_formula& first, const fol_formula& second,
		const fol_term* src_term, unsigned int& parameter)
{
	const fol_term* dst = NULL;
	if (!unify(first, second, src_term, dst)
	 || dst == NULL || dst->type != fol_term_type::PARAMETER)
		return false;
	parameter = dst->parameter;
	return true;
}


/**
 * Functions for easily constructing first-order logic expressions in code.
 */


fol_term* fol_formula::new_variable(unsigned int variable) {
	fol_term* term = (fol_term*) malloc(sizeof(fol_term));
	if (term == NULL) {
		fprintf(stderr, "fol_formula.new_variable ERROR: Out of memory.\n");
		return NULL;
	}
	term->type = fol_term_type::VARIABLE;
	term->variable = variable;
	return term;
}

fol_term* fol_formula::new_constant(unsigned int constant) {
	fol_term* term = (fol_term*) malloc(sizeof(fol_term));
	if (term == NULL) {
		fprintf(stderr, "fol_formula.new_constant ERROR: Out of memory.\n");
		return NULL;
	}
	term->type = fol_term_type::CONSTANT;
	term->constant = constant;
	return term;
}

fol_term* fol_formula::new_parameter(unsigned int parameter) {
	fol_term* term = (fol_term*) malloc(sizeof(fol_term));
	if (term == NULL) {
		fprintf(stderr, "fol_formula.new_parameter ERROR: Out of memory.\n");
		return NULL;
	}
	term->type = fol_term_type::PARAMETER;
	term->parameter = parameter;
	return term;
}

fol_formula* new_atom_helper(unsigned int predicate,
		const fol_term& arg1, const fol_term& arg2)
{
	fol_formula* atom;
	if (!new_fol_formula(atom)) return NULL;
	atom->reference_count = 1;
	atom->type = fol_formula_type::ATOM;
	atom->atom.predicate = predicate;
	atom->atom.arg1 = arg1;
	atom->atom.arg2 = arg2;
	return atom;
}

inline fol_formula* fol_formula::new_atom(
		unsigned int predicate,
		fol_term* arg1, fol_term* arg2)
{
	return new_atom_helper(predicate, *arg1, *arg2);
}

inline fol_formula* fol_formula::new_atom(unsigned int predicate, fol_term* arg1) {
	return new_atom_helper(predicate, *arg1, fol_formula::EMPTY_TERM);
}

inline fol_formula* fol_formula::new_atom(unsigned int predicate) {
	return new_atom_helper(predicate, fol_formula::EMPTY_TERM, fol_formula::EMPTY_TERM);
}

inline fol_formula* fol_formula::new_equals(
		fol_term arg1, fol_term arg2)
{
	fol_formula* equals;
	if (!new_fol_formula(equals)) return NULL;
	equals->reference_count = 1;
	equals->type = fol_formula_type::EQUALS;
	equals->equals.arg1 = arg1;
	equals->equals.arg2 = arg2;
	return equals;
}

inline fol_formula* fol_formula::new_true() {
	FOL_TRUE.reference_count++;
	return &FOL_TRUE;
}

inline fol_formula* fol_formula::new_false() {
	FOL_FALSE.reference_count++;
	return &FOL_FALSE;
}

template<fol_formula_type Operator, unsigned int Index>
inline void new_fol_array_helper(fol_formula** operands, fol_formula* arg)
{
	operands[Index] = arg;
}

template<fol_formula_type Operator, unsigned int Index, typename... Args>
inline void new_fol_array_helper(fol_formula** operands, fol_formula* arg, Args&&... args)
{
	operands[Index] = arg;
	new_fol_array_helper<Operator, Index + 1>(operands, std::forward<Args>(args)...);
}

template<fol_formula_type Operator, typename... Args>
inline fol_formula* new_fol_array(fol_formula* arg, Args&&... args)
{
	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = Operator;
	formula->array.length = 1 + sizeof...(Args);
	formula->array.operands = (fol_formula**) malloc(sizeof(fol_formula*) * (1 + sizeof...(Args)));
	if (formula->array.operands == NULL) {
		free(formula); return NULL;
	}
	new_fol_array_helper<Operator, 0>(formula->array.operands, arg, std::forward<Args>(args)...);
	return formula;
}

template<fol_formula_type Operator, template<typename> class Array>
inline fol_formula* new_fol_array(const Array<fol_formula*>& operands)
{
	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = Operator;
	formula->array.length = operands.length;
	formula->array.operands = (fol_formula**) malloc(sizeof(fol_formula*) * operands.length);
	if (formula->array.operands == NULL) {
		free(formula); return NULL;
	}
	for (unsigned int i = 0; i < operands.length; i++)
		formula->array.operands[i] = operands[i];
	return formula;
}

template<typename... Args>
inline fol_formula* fol_formula::new_and(Args&&... args) {
	return new_fol_array<fol_formula_type::AND>(std::forward<Args>(args)...);
}

template<typename... Args>
inline fol_formula* fol_formula::new_or(Args&&... args) {
	return new_fol_array<fol_formula_type::OR>(std::forward<Args>(args)...);
}

template<typename... Args>
inline fol_formula* fol_formula::new_iff(Args&&... args) {
	return new_fol_array<fol_formula_type::IFF>(std::forward<Args>(args)...);
}

fol_formula* fol_formula::new_if_then(fol_formula* first, fol_formula* second)
{
	if (first == NULL || second == NULL)
		return NULL;

	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = fol_formula_type::IF_THEN;
	formula->binary.left = first;
	formula->binary.right = second;
	return formula;
}

fol_formula* fol_formula::new_not(fol_formula* operand)
{
	if (operand == NULL) return NULL;

	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = fol_formula_type::NOT;
	formula->unary.operand = operand;
	return formula;
}

template<fol_formula_type QuantifierType>
fol_formula* new_fol_quantifier(unsigned int variable, fol_formula* operand)
{
	if (operand == NULL) return NULL;

	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = QuantifierType;
	formula->quantifier.variable = variable;
	formula->quantifier.operand = operand;
	return formula;
}

inline fol_formula* fol_formula::new_for_all(unsigned int variable, fol_formula* operand) {
	return new_fol_quantifier<fol_formula_type::FOR_ALL>(variable, operand);
}

inline fol_formula* fol_formula::new_exists(unsigned int variable, fol_formula* operand) {
	return new_fol_quantifier<fol_formula_type::EXISTS>(variable, operand);
}



/**
 * Below is code for canonicalizing first-order formulas.
 */


int_fast8_t compare(
		const fol_formula&,
		const fol_formula&);

inline int_fast8_t compare(
		const fol_term& first,
		const fol_term& second)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;
	switch (first.type) {
	case fol_term_type::CONSTANT:
		if (first.constant < second.constant) return -1;
		else if (first.constant > second.constant) return 1;
		else return 0;
	case fol_term_type::VARIABLE:
		if (first.variable < second.variable) return -1;
		else if (first.variable > second.variable) return 1;
		else return 0;
	case fol_term_type::PARAMETER:
		if (first.parameter < second.parameter) return -1;
		else if (first.parameter > second.parameter) return 1;
		else return 0;
	case fol_term_type::NONE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline int_fast8_t compare(
		const fol_term* first,
		const fol_term* second)
{
	return compare(*first, *second);
}

inline int_fast8_t compare_terms(
		const fol_term& first_arg1, const fol_term& first_arg2,
		const fol_term& second_arg1, const fol_term& second_arg2)
{
	int_fast8_t result = compare(first_arg1, second_arg1);
	if (result != 0) return result;

	return compare(first_arg2, second_arg2);
}

inline int_fast8_t compare(
		const fol_atom& first,
		const fol_atom& second)
{
	if (first.predicate < second.predicate) return -1;
	else if (first.predicate > second.predicate) return 1;

	return compare_terms(first.arg1, first.arg2, second.arg1, second.arg2);
}

inline int_fast8_t compare(
		const fol_equals& first,
		const fol_equals& second)
{
	return compare_terms(first.arg1, first.arg2, second.arg1, second.arg2);
}

inline int_fast8_t compare(
		const fol_unary_formula& first,
		const fol_unary_formula& second)
{
	return compare(*first.operand, *second.operand);
}

inline int_fast8_t compare(
		const fol_binary_formula& first,
		const fol_binary_formula& second)
{
	int_fast8_t result = compare(*first.left, *second.left);
	if (result != 0) return result;
	return compare(*first.right, *second.right);
}

inline int_fast8_t compare(
		const fol_array_formula& first,
		const fol_array_formula& second)
{
	if (first.length < second.length) return -1;
	else if (first.length > second.length) return 1;
	for (unsigned int i = 0; i < first.length; i++) {
		int_fast8_t result = compare(*first.operands[i], *second.operands[i]);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const fol_quantifier& first,
		const fol_quantifier& second)
{
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
	return compare(*first.operand, *second.operand);
}

int_fast8_t compare(
		const fol_formula& first,
		const fol_formula& second)
{
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return compare(first.atom, second.atom);
	case fol_formula_type::EQUALS:
		return compare(first.equals, second.equals);
	case fol_formula_type::NOT:
		return compare(first.unary, second.unary);
	case fol_formula_type::IF_THEN:
		return compare(first.binary, second.binary);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return compare(first.array, second.array);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return compare(first.quantifier, second.quantifier);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

inline bool operator < (
		const fol_formula& first,
		const fol_formula& second)
{
	return compare(first, second) < 0;
}


/* forward declarations */
bool relabel_variables(fol_formula&, array_map<unsigned int, unsigned int>&);


inline bool relabel_variables(fol_term& term,
		array_map<unsigned int, unsigned int>& variable_map)
{
	unsigned int index;
	switch (term.type) {
	case fol_term_type::CONSTANT:
	case fol_term_type::PARAMETER:
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

inline bool relabel_variables(fol_equals& equals,
		array_map<unsigned int, unsigned int>& variable_map)
{
	return relabel_variables(equals.arg1, variable_map)
		&& relabel_variables(equals.arg2, variable_map);
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
	case fol_formula_type::IF_THEN:
		return relabel_variables(*formula.binary.left, variable_map)
			&& relabel_variables(*formula.binary.right, variable_map);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		for (unsigned int i = 0; i < formula.array.length; i++)
			if (!relabel_variables(*formula.array.operands[i], variable_map)) return false;
		return true;
	case fol_formula_type::NOT:
		return relabel_variables(*formula.unary.operand, variable_map);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return relabel_variables(formula.quantifier, variable_map);
	case fol_formula_type::ATOM:
		return relabel_variables(formula.atom, variable_map);
	case fol_formula_type::EQUALS:
		return relabel_variables(formula.equals, variable_map);
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
bool operator == (const fol_scope&, const fol_scope&);
int_fast8_t compare(const fol_scope&, const fol_scope&);
void shift_variables(fol_scope&, unsigned int);

template<bool AllConstantsDistinct>
bool canonicalize_scope(const fol_formula&, fol_scope&, array_map<unsigned int, unsigned int>&);


struct fol_commutative_scope {
	array<fol_scope> children;
	array<fol_scope> negated;

	fol_commutative_scope() : children(4), negated(4) { }

	~fol_commutative_scope() {
		for (unsigned int i = 0; i < children.length; i++) free(children[i]);
		for (unsigned int i = 0; i < negated.length; i++) free(negated[i]);
	}

	static inline void move(const fol_commutative_scope& src, fol_commutative_scope& dst) {
		core::move(src.children, dst.children);
		core::move(src.negated, dst.negated);
	}
};

struct fol_noncommutative_scope {
	array<fol_scope> left, left_negated;
	array<fol_scope> right, right_negated;

	fol_noncommutative_scope() : left(4), left_negated(4), right(4), right_negated(4) { }

	~fol_noncommutative_scope() {
		for (unsigned int i = 0; i < left.length; i++) free(left[i]);
		for (unsigned int i = 0; i < left_negated.length; i++) free(left_negated[i]);
		for (unsigned int i = 0; i < right.length; i++) free(right[i]);
		for (unsigned int i = 0; i < right_negated.length; i++) free(right_negated[i]);
	}

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

	~fol_quantifier_scope() {
		free(*operand); free(operand);
	}

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
		fol_equals equals;
		fol_scope* unary;
		fol_commutative_scope commutative;
		fol_noncommutative_scope noncommutative;
		fol_quantifier_scope quantifier;
	};

	fol_scope(fol_formula_type type) : variables(8) {
		if (!init_helper(type))
			exit(EXIT_FAILURE);
	}

	fol_scope(fol_formula_type type, const array<unsigned int>& src_variables) : variables(src_variables.capacity) {
		if (!init_helper(type, src_variables))
			exit(EXIT_FAILURE);
	}

	~fol_scope() { free_helper(); }

	static inline void move(const fol_scope& src, fol_scope& dst)
	{
		dst.type = src.type;
		core::move(src.variables, dst.variables);
		switch (src.type) {
		case fol_formula_type::AND:
		case fol_formula_type::OR:
		case fol_formula_type::IFF:
			fol_commutative_scope::move(src.commutative, dst.commutative); return;
		case fol_formula_type::IF_THEN:
			fol_noncommutative_scope::move(src.noncommutative, dst.noncommutative); return;
		case fol_formula_type::FOR_ALL:
		case fol_formula_type::EXISTS:
			fol_quantifier_scope::move(src.quantifier, dst.quantifier); return;
		case fol_formula_type::ATOM:
			fol_atom::move(src.atom, dst.atom); return;
		case fol_formula_type::EQUALS:
			fol_equals::move(src.equals, dst.equals); return;
		case fol_formula_type::NOT:
			dst.unary = src.unary; return;
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return;
		}
		fprintf(stderr, "fol_scope.move ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void swap(fol_scope& first, fol_scope& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(fol_scope); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void free(fol_scope& scope) {
		scope.free_helper();
		core::free(scope.variables);
	}

private:
	inline bool init_helper(fol_formula_type scope_type) {
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
		case fol_formula_type::EQUALS:
			new (&equals) fol_equals(); return true;
		case fol_formula_type::NOT:
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return true;
		}
		fprintf(stderr, "fol_scope.init_helper ERROR: Unrecognized fol_formula_type.\n");
		return false;
	}

	inline bool init_helper(fol_formula_type scope_type, const array<unsigned int>& src_variables) {
		for (unsigned int i = 0; i < src_variables.length; i++)
			variables[i] = src_variables[i];
		variables.length = src_variables.length;
		return init_helper(scope_type);
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
		case fol_formula_type::EQUALS:
			equals.~fol_equals(); return true;
		case fol_formula_type::NOT:
			core::free(*unary); core::free(unary); return true;
		case fol_formula_type::TRUE:
		case fol_formula_type::FALSE:
			return true;
		}
		fprintf(stderr, "fol_scope.free_helper ERROR: Unrecognized fol_formula_type.\n");
		exit(EXIT_FAILURE);
	}

	friend bool init(fol_scope&, fol_formula_type, unsigned int);
	friend bool init(fol_scope&, fol_formula_type, const array<unsigned int>&);
};

inline bool init(fol_scope& scope, fol_formula_type type, unsigned int variable_capacity = 8) {
	if (!array_init(scope.variables, variable_capacity)) {
		return false;
	} else if (!scope.init_helper(type)) {
		free(scope.variables);
		return false;
	}
	return true;
}

inline bool init(fol_scope& scope, fol_formula_type type,
		const array<unsigned int>& src_variables)
{
	if (!array_init(scope.variables, src_variables.capacity)) {
		return false;
	} else if (!scope.init_helper(type, src_variables)) {
		free(scope.variables);
		return false;
	}
	return true;
}

inline bool operator != (const fol_scope& first, const fol_scope& second) {
	return !(first == second);
}

inline bool operator == (const fol_commutative_scope& first, const fol_commutative_scope& second)
{
	if (first.children.length != second.children.length
	 || first.negated.length != second.negated.length)
		return false;
	for (unsigned int i = 0; i < first.children.length; i++)
		if (first.children[i] != second.children[i]) return false;
	for (unsigned int i = 0; i < first.negated.length; i++)
		if (first.negated[i] != second.negated[i]) return false;
	return true;
}

inline bool operator == (const fol_noncommutative_scope& first, const fol_noncommutative_scope& second)
{
	if (first.left.length != second.left.length || first.left_negated.length != second.left_negated.length
	 || first.right.length != second.right.length || first.right_negated.length != second.right_negated.length)
		return false;
	for (unsigned int i = 0; i < first.left.length; i++)
		if (first.left[i] != second.left[i]) return false;
	for (unsigned int i = 0; i < first.left_negated.length; i++)
		if (first.left_negated[i] != second.left_negated[i]) return false;
	for (unsigned int i = 0; i < first.right.length; i++)
		if (first.right[i] != second.right[i]) return false;
	for (unsigned int i = 0; i < first.right_negated.length; i++)
		if (first.right_negated[i] != second.right_negated[i]) return false;
	return true;
}

inline bool operator == (const fol_quantifier_scope& first, const fol_quantifier_scope& second) {
	return first.variable == second.variable
		&& *first.operand == *second.operand;
}

inline bool operator == (const fol_scope& first, const fol_scope& second)
{
	if (first.type != second.type)
		return false;
	switch (first.type) {
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return first.commutative == second.commutative;
	case fol_formula_type::IF_THEN:
		return first.noncommutative == second.noncommutative;
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return first.quantifier == second.quantifier;
	case fol_formula_type::ATOM:
		return first.atom == second.atom;
	case fol_formula_type::EQUALS:
		return first.equals == second.equals;
	case fol_formula_type::NOT:
		return *first.unary == *second.unary;
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return true;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized fol_formula_type when comparing fol_scopes.\n");
	exit(EXIT_FAILURE);
}

inline int_fast8_t compare(
		const fol_commutative_scope& first,
		const fol_commutative_scope& second)
{
	if (first.children.length < second.children.length) return -1;
	else if (first.children.length > second.children.length) return 1;
	else if (first.negated.length < second.negated.length) return -1;
	else if (first.negated.length > second.negated.length) return 1;

	for (unsigned int i = 0; i < first.children.length; i++) {
		int_fast8_t result = compare(first.children[i], second.children[i]);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.negated.length; i++) {
		int_fast8_t result = compare(first.negated[i], second.negated[i]);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const fol_noncommutative_scope& first,
		const fol_noncommutative_scope& second)
{
	if (first.left.length < second.left.length) return -1;
	else if (first.left.length > second.left.length) return 1;
	else if (first.left_negated.length < second.left_negated.length) return -1;
	else if (first.left_negated.length > second.left_negated.length) return 1;
	else if (first.right.length < second.right.length) return -1;
	else if (first.right.length > second.right.length) return 1;
	else if (first.right_negated.length < second.right_negated.length) return -1;
	else if (first.right_negated.length > second.right_negated.length) return 1;

	for (unsigned int i = 0; i < first.left.length; i++) {
		int_fast8_t result = compare(first.left[i], second.left[i]);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.left_negated.length; i++) {
		int_fast8_t result = compare(first.left_negated[i], second.left_negated[i]);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.right.length; i++) {
		int_fast8_t result = compare(first.right[i], second.right[i]);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.right_negated.length; i++) {
		int_fast8_t result = compare(first.right_negated[i], second.right_negated[i]);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const fol_quantifier_scope& first,
		const fol_quantifier_scope& second)
{
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
	return compare(*first.operand, *second.operand);
}

int_fast8_t compare(
		const fol_scope& first,
		const fol_scope& second)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return compare(first.atom, second.atom);
	case fol_formula_type::EQUALS:
		return compare(first.equals, second.equals);
	case fol_formula_type::NOT:
		return compare(*first.unary, *second.unary);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return compare(first.commutative, second.commutative);
	case fol_formula_type::IF_THEN:
		return compare(first.noncommutative, second.noncommutative);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return compare(first.quantifier, second.quantifier);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_formula_type when comparing fol_scopes.\n");
	exit(EXIT_FAILURE);
}

struct fol_scope_canonicalizer { };

inline bool less_than(
		const fol_scope& first,
		const fol_scope& second,
		const fol_scope_canonicalizer& sorter)
{
	return compare(first, second) < 0;
}

inline void shift_variables(fol_term& term, unsigned int removed_variable) {
	switch (term.type) {
	case fol_term_type::VARIABLE:
		if (term.variable > removed_variable)
			term.variable--;
		return;
	case fol_term_type::CONSTANT:
	case fol_term_type::PARAMETER:
	case fol_term_type::NONE:
		return;
	}
	fprintf(stderr, "shift_variables ERROR: Unrecognized fol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline void shift_variables(fol_atom& atom, unsigned int removed_variable) {
	shift_variables(atom.arg1, removed_variable);
	shift_variables(atom.arg2, removed_variable);
}

inline void shift_variables(fol_equals& equals, unsigned int removed_variable) {
	shift_variables(equals.arg1, removed_variable);
	shift_variables(equals.arg2, removed_variable);
}

inline void shift_variables(fol_commutative_scope& scope, unsigned int removed_variable) {
	for (fol_scope& child : scope.children)
		shift_variables(child, removed_variable);
	for (fol_scope& child : scope.negated)
		shift_variables(child, removed_variable);
}

inline void shift_variables(fol_noncommutative_scope& scope, unsigned int removed_variable) {
	for (fol_scope& child : scope.left)
		shift_variables(child, removed_variable);
	for (fol_scope& child : scope.left_negated)
		shift_variables(child, removed_variable);
	for (fol_scope& child : scope.right)
		shift_variables(child, removed_variable);
	for (fol_scope& child : scope.right_negated)
		shift_variables(child, removed_variable);
}

void shift_variables(fol_scope& scope, unsigned int removed_variable) {
	for (unsigned int i = 0; i < scope.variables.length; i++) {
		if (scope.variables[i] > removed_variable)
			scope.variables[i]--;
	}
	switch (scope.type) {
	case fol_formula_type::ATOM:
		return shift_variables(scope.atom, removed_variable);
	case fol_formula_type::EQUALS:
		return shift_variables(scope.equals, removed_variable);
	case fol_formula_type::NOT:
		return shift_variables(*scope.unary, removed_variable);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return shift_variables(scope.commutative, removed_variable);
	case fol_formula_type::IF_THEN:
		return shift_variables(scope.noncommutative, removed_variable);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		if (scope.quantifier.variable > removed_variable)
			scope.quantifier.variable--;
		return shift_variables(*scope.quantifier.operand, removed_variable);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return;
	}
	fprintf(stderr, "shift_variables ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

fol_formula* scope_to_formula(const fol_scope& scope);

template<bool Negated>
inline fol_formula* scope_to_formula(const fol_scope& scope);

template<>
inline fol_formula* scope_to_formula<false>(const fol_scope& scope) {
	return scope_to_formula(scope);
}

template<>
inline fol_formula* scope_to_formula<true>(const fol_scope& negated)
{
	fol_formula* negation;
	if (!new_fol_formula(negation)) return NULL;
	negation->type = fol_formula_type::NOT;
	negation->reference_count = 1;

	negation->unary.operand = scope_to_formula(negated);
	if (negation->unary.operand == NULL) {
		free(negation); return NULL;
	}
	return negation;
}

template<fol_formula_type ScopeType, bool Negated>
inline fol_formula* scope_to_formula(const fol_scope* scope,
		unsigned int scope_length, fol_formula* first)
{
	static_assert(ScopeType == fol_formula_type::AND
			   || ScopeType == fol_formula_type::OR
			   || ScopeType == fol_formula_type::IFF,
			"ScopeType must be either: AND, OR, IFF.");

	if (scope_length == 0) return first;

	fol_formula* new_formula;
	if (!new_fol_formula(new_formula)) {
		free(*first); if (first->reference_count == 0) free(first);
		return NULL;
	}
	new_formula->type = ScopeType;
	new_formula->reference_count = 1;
	new_formula->array.operands = (fol_formula**) malloc(sizeof(fol_formula*) * (scope_length + 1));
	new_formula->array.length = scope_length + 1;
	if (new_formula->array.operands == NULL) {
		free(*first); if (first->reference_count == 0) free(first);
		free(new_formula); return NULL;
	}
	new_formula->array.operands[0] = first;
	for (unsigned int i = 0; i < scope_length; i++) {
		new_formula->array.operands[i + 1] = scope_to_formula<Negated>(scope[i]);
		if (new_formula->array.operands[i + 1] == NULL) {
			for (unsigned int j = 0; j < i + 1; j++) {
				free(*new_formula->array.operands[j]);
				if (new_formula->array.operands[j]->reference_count == 0)
					free(new_formula->array.operands[j]);
			}
			free(new_formula->array.operands); free(new_formula);
			return NULL;
		}
	}
	return new_formula;
}

template<fol_formula_type ScopeType>
inline fol_formula* scope_to_formula(
		const fol_scope* scope, unsigned int scope_length,
		const fol_scope* negated, unsigned int negated_length)
{
	static_assert(ScopeType == fol_formula_type::AND
			   || ScopeType == fol_formula_type::OR
			   || ScopeType == fol_formula_type::IFF,
			"ScopeType must be either: AND, OR, IFF.");

	if (scope_length == 1 && negated_length == 0)
		return scope_to_formula<false>(scope[0]);
	else if (scope_length == 0 && negated_length == 1)
		return scope_to_formula<true>(negated[0]);

	fol_formula* new_formula;
	if (!new_fol_formula(new_formula)) return NULL;
	new_formula->type = ScopeType;
	new_formula->reference_count = 1;
	new_formula->array.operands = (fol_formula**) malloc(sizeof(fol_formula*) * (scope_length + negated_length));
	new_formula->array.length = scope_length + negated_length;
	if (new_formula->array.operands == NULL) {
		free(new_formula); return NULL;
	}
	for (unsigned int i = 0; i < scope_length; i++) {
		new_formula->array.operands[i] = scope_to_formula<false>(scope[i]);
		if (new_formula->array.operands[i] == NULL) {
			for (unsigned int j = 0; j < i; j++) {
				free(*new_formula->array.operands[j]);
				if (new_formula->array.operands[j]->reference_count == 0)
					free(new_formula->array.operands[j]);
			}
			free(new_formula->array.operands); free(new_formula);
			return NULL;
		}
	} for (unsigned int i = 0; i < negated_length; i++) {
		new_formula->array.operands[scope_length + i] = scope_to_formula<true>(negated[i]);
		if (new_formula->array.operands[scope_length + i] == NULL) {
			for (unsigned int j = 0; j < scope_length + i; j++) {
				free(*new_formula->array.operands[j]);
				if (new_formula->array.operands[j]->reference_count == 0)
					free(new_formula->array.operands[j]);
			}
			free(new_formula->array.operands); free(new_formula);
			return NULL;
		}
	}
	return new_formula;
}

template<fol_formula_type QuantifierType>
inline fol_formula* scope_to_formula(const fol_quantifier_scope& scope)
{
	static_assert(QuantifierType == fol_formula_type::FOR_ALL
			   || QuantifierType == fol_formula_type::EXISTS,
			"QuantifierType is not a quantifier.");

	fol_formula* operand = scope_to_formula(*scope.operand);
	if (operand == NULL) return NULL;

	fol_formula* new_formula;
	if (!new_fol_formula(new_formula)) return NULL;
	new_formula->type = QuantifierType;
	new_formula->reference_count = 1;
	new_formula->quantifier.variable = scope.variable;
	new_formula->quantifier.operand = operand;
	return new_formula;
}

inline fol_formula* scope_to_formula(const fol_scope& scope)
{
	fol_formula* new_formula;
	fol_formula* left; fol_formula* right;
	bool negated = false;

	switch (scope.type) {
	case fol_formula_type::AND:
		return scope_to_formula<fol_formula_type::AND>(
				scope.commutative.children.data, scope.commutative.children.length,
				scope.commutative.negated.data, scope.commutative.negated.length);
	case fol_formula_type::OR:
		return scope_to_formula<fol_formula_type::OR>(
				scope.commutative.children.data, scope.commutative.children.length,
				scope.commutative.negated.data, scope.commutative.negated.length);
	case fol_formula_type::IFF:
		if (scope.commutative.children.length > 0 && scope.commutative.children.last().type == fol_formula_type::FALSE)
			negated = true;
		left = scope_to_formula<fol_formula_type::IFF>(
				scope.commutative.children.data, scope.commutative.children.length - (negated ? 1 : 0),
				scope.commutative.negated.data, scope.commutative.negated.length);
		if (left == NULL) return NULL;

		if (negated) {
			if (!new_fol_formula(new_formula)) {
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			}
			new_formula->type = fol_formula_type::NOT;
			new_formula->reference_count = 1;
			new_formula->unary.operand = left;
			return new_formula;
		} else {
			return left;
		}
	case fol_formula_type::IF_THEN:
		left = scope_to_formula<fol_formula_type::AND>(
				scope.noncommutative.left.data, scope.noncommutative.left.length,
				scope.noncommutative.left_negated.data, scope.noncommutative.left_negated.length);
		if (left == NULL) return NULL;
		right = scope_to_formula<fol_formula_type::OR>(
				scope.noncommutative.right.data, scope.noncommutative.right.length,
				scope.noncommutative.right_negated.data, scope.noncommutative.right_negated.length);
		if (right == NULL) {
			free(*left); if (left->reference_count == 0) free(left);
			return NULL;
		}

		if (!new_fol_formula(new_formula)) {
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

		if (!new_fol_formula(new_formula)) return NULL;
		new_formula->type = fol_formula_type::NOT;
		new_formula->reference_count = 1;
		new_formula->unary.operand = left;
		return new_formula;
	case fol_formula_type::FOR_ALL:
		return scope_to_formula<fol_formula_type::FOR_ALL>(scope.quantifier);
	case fol_formula_type::EXISTS:
		return scope_to_formula<fol_formula_type::EXISTS>(scope.quantifier);
	case fol_formula_type::ATOM:
		if (!new_fol_formula(new_formula)) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			return NULL;
		}
		new_formula->type = fol_formula_type::ATOM;
		new_formula->reference_count = 1;
		new_formula->atom = scope.atom;
		return new_formula;
	case fol_formula_type::EQUALS:
		if (!new_fol_formula(new_formula)) {
			fprintf(stderr, "scope_to_formula ERROR: Out of memory.\n");
			return NULL;
		}
		new_formula->type = fol_formula_type::EQUALS;
		new_formula->reference_count = 1;
		new_formula->equals = scope.equals;
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
	array<unsigned int> variable_union(max(scope_variables.capacity, formula_variables.length + scope_variables.length));
	set_union(variable_union.data, variable_union.length,
			formula_variables.data, formula_variables.length,
			scope_variables.data, scope_variables.length);
	swap(variable_union, scope_variables);
}

inline void recompute_variables(
		const array<fol_scope>& children,
		const array<fol_scope>& negated,
		array<unsigned int>& variables)
{
	variables.clear();
	for (const fol_scope& child : children)
		move_variables(child.variables, variables);
	for (const fol_scope& child : negated)
		move_variables(child.variables, variables);
}

inline void recompute_variables(
		const array<fol_scope>& left, const array<fol_scope>& left_negated,
		const array<fol_scope>& right, const array<fol_scope>& right_negated,
		array<unsigned int>& variables)
{
	variables.clear();
	for (const fol_scope& child : left)
		move_variables(child.variables, variables);
	for (const fol_scope& child : left_negated)
		move_variables(child.variables, variables);
	for (const fol_scope& child : right)
		move_variables(child.variables, variables);
	for (const fol_scope& child : right_negated)
		move_variables(child.variables, variables);
}

inline bool scope_contains(const fol_scope& subscope, const array<fol_scope>& scope, unsigned int& index)
{
	for (index = 0; index < scope.length; index++) {
		auto result = compare(subscope, scope[index]);
		if (result < 0) {
			break;
		} else if (result == 0) {
			return true;
		}
	}
	return false;
}

inline bool scope_contains(const fol_scope& subscope, const array<fol_scope>& scope) {
	unsigned int index;
	return scope_contains(subscope, scope, index);
}

template<fol_formula_type Operator>
bool add_to_scope_helper(
		fol_scope& subscope,
		array<fol_scope>& children,
		array<fol_scope>& negated,
		bool& found_negation)
{
	unsigned int i;
	if (scope_contains(subscope, negated, i)) {
		found_negation = true;
		if (Operator == fol_formula_type::IFF) {
			free(negated[i]);
			shift_left(negated.data + i, negated.length - i - 1);
			negated.length--;
		}
		free(subscope);
		return true;
	}
	found_negation = false;

	/* add `subscope` into the correct (sorted) position in `children` */
	if (scope_contains(subscope, children, i)) {
		/* we found an operand in `children` that is identical to `subscope` */
		if (Operator == fol_formula_type::IFF) {
			free(children[i]);
			shift_left(children.data + i, children.length - i - 1);
			children.length--;
		}
		free(subscope);
		return true;
	}

	/* `subscope` is unique, so insert it at index `i` */
	if (!children.ensure_capacity(children.length + 1))
		return false;
	shift_right(children.data, children.length, i);
	move(subscope, children[i]);
	children.length++;
	return true;
}

template<fol_formula_type Operator>
bool add_to_scope(
		fol_scope& subscope,
		array<fol_scope>& children,
		array<fol_scope>& negated,
		bool& found_negation)
{
	/* check if `subscope` is the negation of any operand in `scope.commutative.children` */
	if (subscope.type == fol_formula_type::NOT) {
		if (!add_to_scope_helper<Operator>(*subscope.unary, negated, children, found_negation)) return false;
		free(subscope.variables);
		free(subscope.unary);
		return true;
	} else if (subscope.type == fol_formula_type::IFF && subscope.commutative.children.length > 0
			&& subscope.commutative.children.last().type == fol_formula_type::FALSE)
	{
		free(subscope.commutative.children.last());
		subscope.commutative.children.length--;
		return add_to_scope_helper<Operator>(subscope, negated, children, found_negation);
	} else {
		return add_to_scope_helper<Operator>(subscope, children, negated, found_negation);
	}
}

template<fol_formula_type Operator>
bool add_to_scope(
		fol_scope& subscope,
		array<fol_scope>& children,
		array<fol_scope>& negated,
		array<unsigned int>& variables,
		bool& found_negation)
{
	/* check if `subscope` is the negation of any operand in `scope.commutative.children` */
	if (subscope.type == fol_formula_type::NOT) {
		if (Operator != fol_formula_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(*subscope.unary, negated, children, found_negation)) return false;
		if (Operator == fol_formula_type::IFF)
			recompute_variables(children, negated, variables);
		free(subscope.variables);
		free(subscope.unary);
	} else if (subscope.type == fol_formula_type::IFF && subscope.commutative.children.length > 0
			&& subscope.commutative.children.last().type == fol_formula_type::FALSE)
	{
		free(subscope.commutative.children.last());
		subscope.commutative.children.length--;
		if (Operator != fol_formula_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(subscope, negated, children, found_negation)) return false;
		if (Operator == fol_formula_type::IFF)
			recompute_variables(children, negated, variables);
	} else {
		if (Operator != fol_formula_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(subscope, children, negated, found_negation)) return false;
		if (Operator == fol_formula_type::IFF)
			recompute_variables(children, negated, variables);
	}
	return true;
}

template<fol_formula_type Operator>
unsigned int intersection_size(
	array<fol_scope>& first,
	array<fol_scope>& second)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IF_THEN
			   || Operator == fol_formula_type::IFF,
			"Operator is not a binary operator.");

	unsigned int intersection_count = 0;
	unsigned int i = 0, j = 0, first_index = 0, second_index = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR || Operator == fol_formula_type::IF_THEN) {
				return 1;
			} else if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			if (Operator == fol_formula_type::IFF) move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			if (Operator == fol_formula_type::IFF) move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	if (Operator != fol_formula_type::IFF) return intersection_count;

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
			   || Operator == fol_formula_type::IF_THEN
			   || Operator == fol_formula_type::IFF,
			"Operator is not a binary operator.");

	unsigned int intersection_count = 0;
	unsigned int i = 0, j = 0, first_index = 0, second_index = 0;
	while (i < first.length && j < second.length)
	{
		if (j == skip_second_index) {
			move(second[j], second[second_index]);
			j++; second_index++; continue;
		}

		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR || Operator == fol_formula_type::IF_THEN) {
				return 1;
			} else if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			if (Operator == fol_formula_type::IFF) move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	while (i < first.length) {
		if (Operator == fol_formula_type::IFF) move(first[i], first[first_index]);
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
void merge_scopes(array<fol_scope>& dst,
	array<fol_scope>& first,
	array<fol_scope>& second)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++;
			} else {
				move(first[i], dst[dst.length]);
				free(second[j]);
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
	array<fol_scope>& first,
	array<fol_scope>& second,
	unsigned int skip_second_index,
	unsigned int& new_second_index)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	unsigned int i = 0, j = 0;
	new_second_index = skip_second_index;
	while (i < first.length && j < second.length)
	{
		if (j == skip_second_index) {
			j++;
			new_second_index = dst.length - 1;
			if (j == second.length) break;
		}

		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == fol_formula_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++;
			} else {
				move(first[i], dst[dst.length]);
				free(second[j]);
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
		if (j == skip_second_index) {
			j++;
			new_second_index = dst.length - 1;
			if (j == second.length) break;
		}
		move(second[j], dst[dst.length]);
		dst.length++; j++;
	}
}

template<fol_formula_type Operator>
inline void merge_scopes(
		array<fol_scope>& src, array<fol_scope>& dst,
		array<fol_scope>& src_negated, array<fol_scope>& dst_negated,
		bool& found_negation)
{
	unsigned int intersection_count = intersection_size<Operator>(src, dst_negated);
	if (intersection_count > 0 && (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR)) {
		for (fol_scope& child : src) free(child);
		for (fol_scope& child : src_negated) free(child);
		found_negation = true; return;
	}

	intersection_count += intersection_size<Operator>(src_negated, dst);
	if (intersection_count > 0 && (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR)) {
		for (fol_scope& child : src) free(child);
		for (fol_scope& child : src_negated) free(child);
		found_negation = true; return;
	} else if (Operator == fol_formula_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	if (Operator == fol_formula_type::IFF
	 && src.length == 0 && dst.length == 0
	 && src_negated.length == 0 && dst_negated.length == 0)
		return; /* this happens if the elements of `src` are all negations of the elements of `dst` (and vice versa) */

	/* merge the two scopes */
	array<fol_scope> both = array<fol_scope>(max((size_t) 1, src.length + dst.length));
	merge_scopes<Operator>(both, src, dst);
	swap(both, dst); both.clear();

	array<fol_scope> both_negated = array<fol_scope>(max((size_t) 1, src_negated.length + dst_negated.length));
	merge_scopes<Operator>(both_negated, src_negated, dst_negated);
	swap(both_negated, dst_negated);
}

template<fol_formula_type Operator>
inline void merge_scopes(
		array<fol_scope>& src, array<fol_scope>& dst,
		array<fol_scope>& src_negated, array<fol_scope>& dst_negated,
		bool& found_negation, unsigned int skip_dst_index, unsigned int& new_dst_index)
{
	unsigned int intersection_count = intersection_size<Operator>(src, dst_negated);
	if (intersection_count > 0 && (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR)) {
		for (fol_scope& child : src) free(child);
		for (fol_scope& child : src_negated) free(child);
		found_negation = true; return;
	}

	intersection_count += intersection_size<Operator>(src_negated, dst, skip_dst_index);
	if (intersection_count > 0 && (Operator == fol_formula_type::AND || Operator == fol_formula_type::OR)) {
		for (fol_scope& child : src) free(child);
		for (fol_scope& child : src_negated) free(child);
		found_negation = true; return;
	} else if (Operator == fol_formula_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	if (Operator == fol_formula_type::IFF
	 && src.length == 0 && dst.length == 0
	 && src_negated.length == 0 && dst_negated.length == 0)
		return; /* this happens if the elements of `src` are all negations of the elements of `dst` (and vice versa) */

	/* merge the two scopes */
	array<fol_scope> both = array<fol_scope>(max((size_t) 1, src.length + dst.length));
	merge_scopes<Operator>(both, src, dst, skip_dst_index, new_dst_index);
	swap(both, dst); both.clear();

	array<fol_scope> both_negated = array<fol_scope>(max((size_t) 1, src_negated.length + dst_negated.length));
	merge_scopes<Operator>(both_negated, src_negated, dst_negated);
	swap(both_negated, dst_negated);
}

inline bool negate_iff(fol_scope& scope)
{
	if (!scope.commutative.children.ensure_capacity(scope.commutative.children.length + 1)) {
		return false;
	} else if (scope.commutative.children.length > 0 && scope.commutative.children.last().type == fol_formula_type::FALSE) {
		free(scope.commutative.children.last());
		scope.commutative.children.length--;
	} else {
		if (!init(scope.commutative.children[scope.commutative.children.length], fol_formula_type::FALSE))
			return false;
		scope.commutative.children.length++;
	}
	return true;
}

bool negate_scope(fol_scope& scope)
{
	if (scope.type == fol_formula_type::TRUE) {
		scope.type = fol_formula_type::FALSE;
	} else if (scope.type == fol_formula_type::FALSE) {
		scope.type = fol_formula_type::TRUE;
	} else if (scope.type == fol_formula_type::NOT) {
		fol_scope* operand = scope.unary;
		free(scope.variables);
		move(*operand, scope);
		free(operand);
	} else if (scope.type == fol_formula_type::IFF) {
		return negate_iff(scope);
	} else {
		fol_scope* operand = (fol_scope*) malloc(sizeof(fol_scope));
		if (operand == NULL) {
			fprintf(stderr, "negate_scope ERROR: Out of memory.\n");
			return false;
		}

		move(scope, *operand);
		if (!init(scope, fol_formula_type::NOT, operand->variables)) {
			move(*operand, scope); free(operand);
			return false;
		}
		scope.unary = operand;
	}
	return true;
}

inline bool are_negations(fol_scope& left, fol_scope& right) {
	if ((left.type == fol_formula_type::NOT && *left.unary == right)
	 || (right.type == fol_formula_type::NOT && *right.unary == left))
		return true;
	if (left.type == fol_formula_type::IFF && right.type == fol_formula_type::IFF) {
		if (left.commutative.children.length > 0 && left.commutative.children.last().type == fol_formula_type::FALSE) {
			/* `left` is a negated IFF expression */
			if (right.commutative.children.length > 0 && right.commutative.children.last().type == fol_formula_type::FALSE) {
				/* `right` is a negated IFF expression */
				return false;
			} else {
				left.commutative.children.length--;
				bool negated = (left == right);
				left.commutative.children.length++;
				return negated;
			}
		} else {
			if (right.commutative.children.length > 0 && right.commutative.children.last().type == fol_formula_type::FALSE) {
				/* `right` is a negated IFF expression */
				right.commutative.children.length--;
				bool negated = (left == right);
				right.commutative.children.length++;
				return negated;
			} else {
				return false;
			}
		}
	}
	return false;
}

template<fol_formula_type Operator, bool AllConstantsDistinct>
bool canonicalize_commutative_scope(
		const fol_array_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	if (!init(out, Operator)) return false;

	fol_scope& next = *((fol_scope*) alloca(sizeof(fol_scope)));
	for (unsigned int i = 0; i < src.length; i++) {
		if (!canonicalize_scope<AllConstantsDistinct>(*src.operands[i], next, variable_map)) {
			free(out); return false;
		}

		if (next.type == fol_formula_type::FALSE) {
			free(next);
			if (Operator == fol_formula_type::AND) {
				free(out);
				return init(out, fol_formula_type::FALSE);
			} else if (Operator == fol_formula_type::IFF) {
				if (!negate_iff(out)) return false;
			}
		} else if (next.type == fol_formula_type::TRUE) {
			free(next);
			if (Operator == fol_formula_type::OR) {
				free(out);
				return init(out, fol_formula_type::TRUE);
			}
		} else if (next.type == Operator) {
			bool found_negation;
			merge_scopes<Operator>(next.commutative.children, out.commutative.children,
					next.commutative.negated, out.commutative.negated, found_negation);
			free(next.commutative.children);
			free(next.commutative.negated);
			if (Operator == fol_formula_type::IFF)
				recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
			else move_variables(next.variables, out.variables);
			free(next.variables);
			if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(out);
					return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(out);
					return init(out, fol_formula_type::TRUE);
				} else if (Operator == fol_formula_type::IFF) {
					recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
					if (!negate_iff(out)) return false;
				}
			}
		} else {
			bool found_negation;
			if (!add_to_scope<Operator>(next, out.commutative.children, out.commutative.negated, out.variables, found_negation)) {
				free(out); free(next); return false;
			} else if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(out); return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(out); return init(out, fol_formula_type::TRUE);
				} else if (Operator == fol_formula_type::IFF) {
					if (!negate_iff(out)) return false;
				}
			}
		}
	}

	if (out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
		free(out);
		if (Operator == fol_formula_type::AND || Operator == fol_formula_type::IFF)
			return init(out, fol_formula_type::TRUE);
		else return init(out, fol_formula_type::FALSE);
	} else if (out.commutative.children.length == 1 && out.commutative.negated.length == 0) {
		move(out.commutative.children[0], next);
		out.commutative.children.clear();
		free(out); move(next, out);
	} else if (out.commutative.children.length == 0 && out.commutative.negated.length == 1) {
		move(out.commutative.negated[0], next);
		out.commutative.negated.clear();
		free(out); move(next, out);
		if (!negate_scope(out)) {
			free(out); return false;
		}
	}
	return true;
}

template<bool AllConstantsDistinct>
bool canonicalize_conditional_scope(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope& left = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (!canonicalize_scope<AllConstantsDistinct>(*src.left, left, variable_map)) return false;

	if (left.type == fol_formula_type::FALSE) {
		free(left);
		return init(out, fol_formula_type::TRUE);
	} else if (left.type == fol_formula_type::TRUE) {
		free(left);
		return canonicalize_scope<AllConstantsDistinct>(*src.right, out, variable_map);
	}

	if (!canonicalize_scope<AllConstantsDistinct>(*src.right, out, variable_map))
		return false;

	if (out == left) {
		/* consider the case where `*src.left` and `*src.right` are identical */
		free(out); free(left);
		return init(out, fol_formula_type::TRUE);
	} else if (out.type == fol_formula_type::FALSE) {
		free(out); move(left, out);
		if (!negate_scope(out)) {
			free(out); return false;
		}
	} else if (out.type == fol_formula_type::TRUE) {
		free(left); /* this is a no-op */
	} else if (are_negations(out, left)) {
		free(left); /* we have the case `A => ~A` or `~A => A`, which is also a no-op */
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
				fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
				move(child, temp);
				merge_scopes<fol_formula_type::AND>(temp.noncommutative.left, out.noncommutative.left,
						temp.noncommutative.left_negated, out.noncommutative.left_negated, found_negation);
				temp.noncommutative.left.clear();
				temp.noncommutative.left_negated.clear();
				if (found_negation) {
					child.noncommutative.left.clear();
					child.noncommutative.left_negated.clear();
					free(out); free(left);
					return init(out, fol_formula_type::TRUE);
				}

				unsigned int new_index;
				merge_scopes<fol_formula_type::OR>(temp.noncommutative.right, out.noncommutative.right,
						temp.noncommutative.right_negated, out.noncommutative.right_negated, found_negation, i, new_index);
				temp.noncommutative.right.clear();
				temp.noncommutative.right_negated.clear();
				if (found_negation) {
					child.noncommutative.left.clear();
					child.noncommutative.left_negated.clear();
					child.noncommutative.right.clear();
					child.noncommutative.right_negated.clear();
					free(out); free(left);
					return init(out, fol_formula_type::TRUE);
				} else {
					free(temp);
				}
				i = new_index;
			}
		} else if (out.type == fol_formula_type::NOT) {
			fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (!init(temp, fol_formula_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			move(*out.unary, temp.noncommutative.right_negated[0]);
			temp.noncommutative.right_negated.length++;
			free(out.variables); free(out.unary);
			move(temp, out);
		} else if (out.type == fol_formula_type::IFF && out.commutative.children.length > 0
				&& out.commutative.children.last().type == fol_formula_type::FALSE)
		{
			fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (!init(temp, fol_formula_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			free(out.commutative.children.last());
			out.commutative.children.length--;
			move(out, temp.noncommutative.right_negated[0]);
			temp.noncommutative.right_negated.length++;
			move(temp, out);
		} else if (out.type != fol_formula_type::IF_THEN) {
			fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (!init(temp, fol_formula_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			move(out, temp.noncommutative.right[0]);
			temp.noncommutative.right.length++;
			move(temp, out);
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
			bool found_negation;
			if (!add_to_scope<fol_formula_type::AND>(left, out.noncommutative.left, out.noncommutative.left_negated, out.variables, found_negation)) {
				free(out); free(left); return false;
			} else if (found_negation) {
				free(out);
				return init(out, fol_formula_type::TRUE);
			}
		}

		/* check if antecedent and consequent have any common operands */
		if (intersection_size<fol_formula_type::IF_THEN>(out.noncommutative.left, out.noncommutative.right) > 0
		 || intersection_size<fol_formula_type::IF_THEN>(out.noncommutative.left_negated, out.noncommutative.right_negated) > 0)
		{
			free(out);
			return init(out, fol_formula_type::TRUE);
		}
	}
	return true;
}

inline bool promote_from_quantifier_scope(
		array<fol_scope>& quantifier_operand,
		array<fol_scope>& dst,
		unsigned int quantifier_variable)
{
	unsigned int next = 0;
	for (unsigned int i = 0; i < quantifier_operand.length; i++) {
		fol_scope& child = quantifier_operand[i];
		if (!child.variables.contains(quantifier_variable)) {
			if (!dst.ensure_capacity(dst.length + 1)) {
				/* finish shifting the remainder of `quantifier_operand` so all
				   of its elements are valid (and we avoid double freeing, for example) */
				while (i < quantifier_operand.length) {
					move(quantifier_operand[i], quantifier_operand[next]);
					next++; i++;
				}
				quantifier_operand.length = next;
				return false;
			}
			shift_variables(child, quantifier_variable);
			move(child, dst[dst.length]);
			dst.length++;

		} else {
			move(quantifier_operand[i], quantifier_operand[next++]);
		}
	}
	quantifier_operand.length = next;
	return true;
}

template<fol_formula_type QuantifierType>
inline bool make_quantifier_scope(fol_scope& out, fol_scope* operand, unsigned int quantifier_variable)
{
	if (!init(out, QuantifierType, operand->variables)) return false;
	out.quantifier.variable = quantifier_variable;
	out.quantifier.operand = operand;

	unsigned int index = out.variables.index_of(quantifier_variable);
	if (index < out.variables.length) {
		shift_left(out.variables.data + index, out.variables.length - index - 1);
		out.variables.length--;
	}
	return true;
}


/* forward declarations */
template<fol_formula_type QuantifierType>
bool process_conditional_quantifier_scope(fol_scope&, fol_scope*, unsigned int);


template<fol_formula_type QuantifierType>
bool process_commutative_quantifier_scope(
		fol_scope& out, fol_scope* operand,
		unsigned int quantifier_variable)
{
	if (!init(out, operand->type)) {
		free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->commutative.children, out.commutative.children, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->commutative.negated, out.commutative.negated, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	}

	if (operand->commutative.children.length == 0 && operand->commutative.negated.length == 0) {
		/* we've moved all children out of the quantifier */
		move_variables(operand->variables, out.variables);
		free(*operand); free(operand);
	} else {
		fol_scope* quantifier_operand = (fol_scope*) malloc(sizeof(fol_scope));
		if (quantifier_operand == NULL) {
			fprintf(stderr, "process_commutative_quantifier_scope ERROR: Out of memory.\n");
			free(out); free(*operand); free(operand); return false;
		}

		if (operand->commutative.children.length == 1 && operand->commutative.negated.length == 0) {
			fol_scope& inner_operand = operand->commutative.children[0];
			move(inner_operand, *quantifier_operand);
			operand->commutative.children.clear();
			free(*operand); free(operand);
		} else if (operand->commutative.children.length == 0 && operand->commutative.negated.length == 1) {
			move(operand->commutative.negated[0], *quantifier_operand);
			operand->commutative.negated.clear();
			if (!negate_scope(*quantifier_operand)) {
				free(*quantifier_operand); free(quantifier_operand);
				free(out); free(*operand); free(operand);
				return false;
			}
			free(*operand); free(operand);
		} else {
			recompute_variables(operand->commutative.children, operand->commutative.negated, operand->variables);
			move(*operand, *quantifier_operand); free(operand);
		}

		/* removing an operand can allow movement from its children to the
		   parent, so check if the new operand allows movement */
		fol_scope& quantifier = *((fol_scope*) alloca(sizeof(fol_scope)));
		if (quantifier_operand->type != out.type && (quantifier_operand->type == fol_formula_type::AND || quantifier_operand->type == fol_formula_type::OR)) {
			if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
				return false;
		} else if (quantifier_operand->type == fol_formula_type::IF_THEN) {
			if (!process_conditional_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
				return false;
		} else if (!make_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable)) {
			free(out); return false;
		}

		if (out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
			free(out);
			move(quantifier, out);
		} else {
			bool found_negation;
			if (out.type == fol_formula_type::AND) {
				if (!add_to_scope<fol_formula_type::AND>(quantifier, out.commutative.children, out.commutative.negated, found_negation)) {
					free(out); free(quantifier); return false;
				} else if (found_negation) {
					free(out);
					return init(out, fol_formula_type::FALSE);
				}
			} else if (out.type == fol_formula_type::OR) {
				if (!add_to_scope<fol_formula_type::OR>(quantifier, out.commutative.children, out.commutative.negated, found_negation)) {
					free(out); free(quantifier);
				} else if (found_negation) {
					free(out);
					return init(out, fol_formula_type::TRUE);
				}
			}
			recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
		}
	}
	return true;
}

template<fol_formula_type QuantifierType>
bool process_conditional_quantifier_scope(
		fol_scope& out, fol_scope* operand,
		unsigned int quantifier_variable)
{
	if (!init(out, fol_formula_type::IF_THEN)) {
		free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->noncommutative.left, out.noncommutative.left, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->noncommutative.left_negated, out.noncommutative.left_negated, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->noncommutative.right, out.noncommutative.right, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	} else if (!promote_from_quantifier_scope(operand->noncommutative.right_negated, out.noncommutative.right_negated, quantifier_variable)) {
		free(out); free(*operand); free(operand);
		return false;
	}

	if (operand->noncommutative.left.length == 0 && operand->noncommutative.left_negated.length == 0) {
		/* we've moved all children out of the antecedent */
		if (operand->noncommutative.right.length == 0 && operand->noncommutative.right_negated.length == 0) {
			/* we've moved all children out of the consequent */
			free(*operand); free(operand);
		} else {
			fol_scope* quantifier_operand = (fol_scope*) malloc(sizeof(fol_scope));
			if (quantifier_operand == NULL) {
				fprintf(stderr, "process_conditional_quantifier_scope ERROR: Out of memory.\n");
				free(out); free(*operand); free(operand); return false;
			}

			if (operand->noncommutative.right.length == 1 && operand->noncommutative.right_negated.length == 0) {
				move(operand->noncommutative.right[0], *quantifier_operand);
				operand->noncommutative.right.clear();
				free(*operand); free(operand);
			} else if (operand->noncommutative.right.length == 0 && operand->noncommutative.right_negated.length == 1) {
				move(operand->noncommutative.right_negated[0], *quantifier_operand);
				operand->noncommutative.right_negated.clear();
				if (!negate_scope(*quantifier_operand)) {
					free(*quantifier_operand); free(quantifier_operand);
					free(out); free(*operand); free(operand);
					return false;
				}
				free(*operand); free(operand);
			} else {
				if (!init(*quantifier_operand, fol_formula_type::OR)) {
					fprintf(stderr, "process_conditional_quantifier_scope ERROR: Out of memory.\n");
					free(quantifier_operand); free(out); free(*operand); free(operand); return false;
				}
				swap(operand->noncommutative.right, quantifier_operand->commutative.children);
				swap(operand->noncommutative.right_negated, quantifier_operand->commutative.negated);
				swap(operand->variables, quantifier_operand->variables);
				recompute_variables(quantifier_operand->commutative.children, quantifier_operand->commutative.negated, quantifier_operand->variables);
				free(*operand); free(operand);
			}


			/* removing an operand can allow movement from its children to the
			   parent, so check if the new operand allows movement */
			fol_scope& quantifier = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (quantifier_operand->type == fol_formula_type::AND || quantifier_operand->type == fol_formula_type::OR) {
				if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (!make_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<fol_formula_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, fol_formula_type::TRUE);
			}
			recompute_variables(out.noncommutative.left, out.noncommutative.left_negated,
					out.noncommutative.right, out.noncommutative.right_negated, out.variables);
		}
	} else {
		/* the antecedent is non-empty */
		if (operand->noncommutative.right.length == 0 && operand->noncommutative.right_negated.length == 0) {
			/* we've moved all children out of the consequent */
			fol_scope* quantifier_operand = (fol_scope*) malloc(sizeof(fol_scope));
			if (quantifier_operand == NULL) {
				fprintf(stderr, "process_conditional_quantifier_scope ERROR: Out of memory.\n");
				free(out); free(*operand); free(operand); return false;
			}

			if (operand->noncommutative.left.length == 1 && operand->noncommutative.left_negated.length == 0) {
				move(operand->noncommutative.left[0], *quantifier_operand);
				operand->noncommutative.left.clear();
				if (!negate_scope(*quantifier_operand)) {
					free(*quantifier_operand); free(quantifier_operand);
					free(out); free(*operand); free(operand);
					return false;
				}
				free(*operand); free(operand);
			} else if (operand->noncommutative.left.length == 0 && operand->noncommutative.left_negated.length == 1) {
				move(operand->noncommutative.left_negated[0], *quantifier_operand);
				operand->noncommutative.left_negated.clear();
				free(*operand); free(operand);
			} else {
				fol_scope* conjunction = (fol_scope*) malloc(sizeof(fol_scope));
				if (conjunction == NULL) {
					fprintf(stderr, "process_conditional_quantifier_scope ERROR: Out of memory.\n");
					free(quantifier_operand); free(out); free(*operand); free(operand);
				}
				conjunction->type = fol_formula_type::AND;
				move(operand->noncommutative.left, conjunction->commutative.children);
				move(operand->noncommutative.left_negated, conjunction->commutative.negated);
				move(operand->variables, conjunction->variables);
				recompute_variables(conjunction->commutative.children, conjunction->commutative.negated, conjunction->variables);
				free(operand->noncommutative.right); free(operand->noncommutative.right_negated);
				free(operand);

				if (!init(*quantifier_operand, fol_formula_type::NOT)) {
					free(quantifier_operand); free(out);
					return false;
				}
				quantifier_operand->unary = conjunction;
				move_variables(conjunction->variables, quantifier_operand->variables);
			}


			/* removing an operand can allow movement from its children to the
			   parent, so check if the new operand allows movement */
			fol_scope& quantifier = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (quantifier_operand->type == fol_formula_type::AND || quantifier_operand->type == fol_formula_type::OR) {
				if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (quantifier_operand->type == fol_formula_type::IF_THEN) {
				if (!process_conditional_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (!make_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<fol_formula_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, fol_formula_type::TRUE);
			}
			recompute_variables(out.noncommutative.left, out.noncommutative.left_negated,
					out.noncommutative.right, out.noncommutative.right_negated, out.variables);
		} else {
			fol_scope& quantifier = *((fol_scope*) alloca(sizeof(fol_scope)));
			recompute_variables(operand->noncommutative.left, operand->noncommutative.left_negated,
					operand->noncommutative.right, operand->noncommutative.right_negated, operand->variables);
			if (!make_quantifier_scope<QuantifierType>(quantifier, operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<fol_formula_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, fol_formula_type::TRUE);
			}
		}

		if (out.noncommutative.left.length == 0 && out.noncommutative.left_negated.length == 0) {
			/* the antecendent of the new (parent) conditional is empty, so change the node into a disjunction */
			if (out.noncommutative.right.length == 1 && out.noncommutative.right_negated.length == 0) {
				fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
				move(out.noncommutative.right[0], temp);
				out.noncommutative.right.clear();
				free(out); move(temp, out);
			} else if (out.noncommutative.right.length == 0 && out.noncommutative.right_negated.length == 1) {
				fol_scope& temp = *((fol_scope*) alloca(sizeof(fol_scope)));
				move(out.noncommutative.right_negated[0], temp);
				out.noncommutative.right.clear();
				free(out); move(temp, out);
				if (!negate_scope(out)) {
					free(out); return false;
				}
			} else {
				array<fol_scope>& temp = *((array<fol_scope>*) alloca(sizeof(array<fol_scope>)));
				array<fol_scope>& temp_negated = *((array<fol_scope>*) alloca(sizeof(array<fol_scope>)));
				free(out.noncommutative.left);
				free(out.noncommutative.left_negated);
				move(out.noncommutative.right, temp);
				move(out.noncommutative.right_negated, temp_negated);
				out.type = fol_formula_type::OR;
				move(temp, out.commutative.children);
				move(temp_negated, out.commutative.negated);
				recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
			}
		} else {
			recompute_variables(out.noncommutative.left, out.noncommutative.left_negated,
					out.noncommutative.right, out.noncommutative.right_negated, out.variables);
		}
	}
	return true;
}

template<fol_formula_type QuantifierType, bool AllConstantsDistinct>
bool canonicalize_quantifier_scope(
		const fol_quantifier& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	static_assert(QuantifierType == fol_formula_type::FOR_ALL
			   || QuantifierType == fol_formula_type::EXISTS,
			"QuantifierType is not a quantifier.");

	unsigned int quantifier_variable;
	fol_scope* operand = (fol_scope*) malloc(sizeof(fol_scope));
	if (operand == NULL) {
		fprintf(stderr, "canonicalize_quantifier_scope ERROR: Out of memory.\n");
		return false;
	} else if (!new_variable(src.variable, quantifier_variable, variable_map)
			|| !canonicalize_scope<AllConstantsDistinct>(*src.operand, *operand, variable_map))
	{
		free(operand); return false;
	}

	variable_map.size--;

	/* check if the operand has any instances of the quantified variable */
	if (!operand->variables.contains(quantifier_variable)) {
		move(*operand, out); free(operand);
		return true;
	}

	if (operand->type == fol_formula_type::AND || operand->type == fol_formula_type::OR) {
		return process_commutative_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	} else if (operand->type == fol_formula_type::IF_THEN) {
		return process_conditional_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	} else {
		return make_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	}
	return true;
}

template<bool AllConstantsDistinct>
inline bool canonicalize_negation_scope(
		const fol_unary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!canonicalize_scope<AllConstantsDistinct>(*src.operand, out, variable_map))
		return false;
	if (!negate_scope(out)) {
		free(out); return false;
	}
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
	case fol_term_type::PARAMETER:
		dst.parameter = src.parameter; return true;
	case fol_term_type::VARIABLE:
		index = variable_map.index_of(src.variable);
		if (index < variable_map.size) {
			dst.variable = variable_map.values[index];
		} else {
			fprintf(stderr, "relabel_variables ERROR: Undefined variable.\n");
			return false;
		}

		index = linear_search(variables.data, dst.variable, 0, variables.length);
		if (index == variables.length) {
			if (!variables.ensure_capacity(variables.length + 1)) return false;
			variables[variables.length] = dst.variable;
			variables.length++;
		} else if (variables[index] != dst.variable) {
			if (!variables.ensure_capacity(variables.length + 1)) return false;
			shift_right(variables.data, variables.length, index);
			variables[index] = dst.variable;
			variables.length++;
		}
		return true;
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "canonicalize_term ERROR: Unrecognized fol_term_type.\n");
	return false;
}

inline bool canonicalize_atom(const fol_atom& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!init(out, fol_formula_type::ATOM)) return false;

	out.atom.predicate = src.predicate;
	if (!canonicalize_term(src.arg1, out.atom.arg1, variable_map, out.variables)
	 || !canonicalize_term(src.arg2, out.atom.arg2, variable_map, out.variables)) {
		return false;
	}

	return true;
}

template<bool AllConstantsDistinct>
inline bool canonicalize_equals(const fol_equals& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!init(out, fol_formula_type::EQUALS)) return false;

	if (!canonicalize_term(src.arg1, out.equals.arg1, variable_map, out.variables)
	 || !canonicalize_term(src.arg2, out.equals.arg2, variable_map, out.variables)) {
		return false;
	}

	if (out.equals.arg2 < out.equals.arg1)
		swap(out.equals.arg1, out.equals.arg2);

	if (out.equals.arg1 == out.equals.arg2) {
		free(out);
		return init(out, fol_formula_type::TRUE);
	} else if (AllConstantsDistinct) {
		switch (out.equals.arg1.type) {
		case fol_term_type::CONSTANT:
			if (out.equals.arg2.type == fol_term_type::CONSTANT) {
				free(out);
				/* since we consider the case `out.equals.arg1 == out.equals.arg2` earlier, the two constants must differ */
				return init(out, fol_formula_type::FALSE);
			} else {
				return true;
			}
		case fol_term_type::PARAMETER:
		case fol_term_type::VARIABLE:
		case fol_term_type::NONE:
			return true;
		}
		fprintf(stderr, "canonicalize_equals ERROR: Unrecognized fol_term_type.\n");
		return false;
	}

	return true;
}

template<bool AllConstantsDistinct>
bool canonicalize_scope(const fol_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	switch (src.type) {
	case fol_formula_type::AND:
		return canonicalize_commutative_scope<fol_formula_type::AND, AllConstantsDistinct>(src.array, out, variable_map);
	case fol_formula_type::OR:
		return canonicalize_commutative_scope<fol_formula_type::OR, AllConstantsDistinct>(src.array, out, variable_map);
	case fol_formula_type::IFF:
		return canonicalize_commutative_scope<fol_formula_type::IFF, AllConstantsDistinct>(src.array, out, variable_map);
	case fol_formula_type::IF_THEN:
		return canonicalize_conditional_scope<AllConstantsDistinct>(src.binary, out, variable_map);
	case fol_formula_type::FOR_ALL:
		return canonicalize_quantifier_scope<fol_formula_type::FOR_ALL, AllConstantsDistinct>(src.quantifier, out, variable_map);
	case fol_formula_type::EXISTS:
		return canonicalize_quantifier_scope<fol_formula_type::EXISTS, AllConstantsDistinct>(src.quantifier, out, variable_map);
	case fol_formula_type::NOT:
		return canonicalize_negation_scope<AllConstantsDistinct>(src.unary, out, variable_map);
	case fol_formula_type::ATOM:
		return canonicalize_atom(src.atom, out, variable_map);
	case fol_formula_type::EQUALS:
		return canonicalize_equals<AllConstantsDistinct>(src.equals, out, variable_map);
	case fol_formula_type::TRUE:
		return init(out, fol_formula_type::TRUE);
	case fol_formula_type::FALSE:
		return init(out, fol_formula_type::FALSE);
	}
	fprintf(stderr, "canonicalize_scope ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

struct identity_canonicalizer { };

inline fol_formula* canonicalize(fol_formula& src,
		const identity_canonicalizer& canonicalizer)
{
	src.reference_count++;
	return &src;
}

template<bool AllConstantsDistinct>
struct standard_canonicalizer { };

template<bool AllConstantsDistinct>
inline fol_formula* canonicalize(const fol_formula& src,
		const standard_canonicalizer<AllConstantsDistinct>& canonicalizer)
{
	array_map<unsigned int, unsigned int> variable_map(16);
	fol_scope& scope = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (!canonicalize_scope<AllConstantsDistinct>(src, scope, variable_map))
		return NULL;
	fol_formula* canonicalized = scope_to_formula(scope);
	free(scope);
	return canonicalized;
}

template<typename Canonicalizer>
bool is_canonical(const fol_formula& src, Canonicalizer& canonicalizer) {
	fol_formula* canonicalized = canonicalize(src, canonicalizer);
	if (canonicalized == NULL) {
		fprintf(stderr, "is_canonical ERROR: Unable to canonicalize formula.\n");
		exit(EXIT_FAILURE);
	}
	bool canonical = (src == *canonicalized);
	free(*canonicalized);
	if (canonicalized->reference_count == 0)
		free(canonicalized);
	return canonical;
}

template<>
constexpr bool is_canonical<identity_canonicalizer>(
		const fol_formula& src, identity_canonicalizer& canonicalizer)
{
	return true;
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
	EQUALS,

	IDENTIFIER,
	SEMICOLON
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
	case tptp_token_type::EQUALS:
		return print('=', stream);
	case tptp_token_type::SEMICOLON:
		return print(';', stream);
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
	case ';':
		return emit_token(tokens, start, start + 1, tptp_token_type::SEMICOLON);
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
	 || next == '~' || next == '!' || next == '?' || next == ';')
	{
		return tptp_emit_symbol(tokens, current, next);
	} else if (next == '=') {
		fpos_t pos;
		fgetpos(input, &pos);
		next = fgetwc(input);
		if (next != '>') {
			fsetpos(input, &pos);
			if (!emit_token(tokens, current, current + 1, tptp_token_type::EQUALS)) return false;
		} else {
			if (!emit_token(tokens, current, current + 2, tptp_token_type::IF_THEN)) return false;
		}
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

template<typename Stream>
bool tptp_lex(array<tptp_token>& tokens, Stream& input, position start = position(1, 1)) {
	position current = start;
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
			 || next == '<' || next == ';')
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
			 || next == '<' || next == ';')
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

inline bool tptp_interpret_argument(
	const string& identifier, fol_term& next_term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	bool contains;
	unsigned int variable = variables.get(identifier, contains);
	if (contains) {
		/* this argument is a variable */
		next_term.variable = variable;
		next_term.type = fol_term_type::VARIABLE;
	} else {
		if (!get_token(identifier, next_term.constant, names))
			return false;
		next_term.type = fol_term_type::CONSTANT;
	}
	return true;
}

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

		if (!tptp_interpret_argument(tokens[index].text, terms[terms.length], names, variables))
			return false;
		index++; terms.length++;

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

	fol_formula* operand;
	if (!new_fol_formula(operand)) return false;
	operand->reference_count = 1;
	if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
		free(operand); return false;
	}

	fol_formula* inner = operand;
	for (unsigned int i = variables.size - 1; i > old_variable_count; i--) {
		fol_formula* quantified;
		if (!new_fol_formula(quantified)) {
			free(*inner); free(inner);
			return false;
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
		fol_formula* operand;
		if (!new_fol_formula(operand)) return false;
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
			/* this is an atomic formula */
			const tptp_token& identifier = tokens[index];
			index++;
			if (index < tokens.length && tokens[index].type == tptp_token_type::EQUALS) {
				/* this is an atomic formula of the form x = y */
				if (!tptp_interpret_argument(identifier.text, formula.equals.arg1, names, variables))
					return false;
				index++;
				if (!expect_token(tokens, index, tptp_token_type::IDENTIFIER, "constant or variable on right-side of equality")
				 || !tptp_interpret_argument(tokens[index].text, formula.equals.arg2, names, variables))
					return false;
				index++;
				formula.type = fol_formula_type::EQUALS;
				formula.reference_count = 1;
			} else {
				/* this is an atomic formula of the form P(T_1,...,T_n) */
				unsigned int predicate;
				if (variables.contains(identifier.text)) {
					fprintf(stderr, "WARNING at %d:%d: Predicate '", identifier.start.line, identifier.start.column);
					print(identifier.text, stderr); print("' is also a variable.\n", stderr);
				} if (!get_token(identifier.text, predicate, names))
					return false;

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
	fol_formula* right;
	if (!new_fol_formula(right)) {
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
template<> struct tptp_operator_type<fol_formula_type::IFF> { static constexpr tptp_token_type type = tptp_token_type::IFF; };

template<fol_formula_type OperatorType>
bool tptp_interpret_binary_sequence(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	fol_formula* left)
{
	array<fol_formula*>& operands = *((array<fol_formula*>*) alloca(sizeof(array<fol_formula*>)));
	if (!array_init(operands, 8)) {
		free(*left); free(left); return false;
	}
	operands[0] = left;
	operands.length = 1;

	while (true) {
		fol_formula* next;
		if (!new_fol_formula(next) || !operands.ensure_capacity(operands.length + 1)
		 || !tptp_interpret_unary_formula(tokens, index, *next, names, variables))
		{
			if (next != NULL) free(next);
			for (unsigned int i = 0; i < operands.length; i++) {
				free(*operands[i]); free(operands[i]);
			}
			free(operands); return false;
		}
		operands[operands.length] = next;
		operands.length++;

		if (index < tokens.length && tokens[index].type == tptp_operator_type<OperatorType>::type)
			index++;
		else break;
	}

	formula.type = OperatorType;
	formula.reference_count = 1;
	formula.array.operands = operands.data;
	formula.array.length = operands.length;
	return true;
}

bool tptp_interpret(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	fol_formula* left;
	if (!new_fol_formula(left)) {
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
		if (!tptp_interpret_binary_sequence<fol_formula_type::IFF>(tokens, index, formula, names, variables, left))
			return false;
	} else {
		move(*left, formula); free(left);
	}
	return true;
}

template<typename Stream>
inline bool parse(Stream& in, fol_formula& formula,
		hash_map<string, unsigned int>& names,
		position start = position(1, 1))
{
	array<tptp_token> tokens = array<tptp_token>(128);
	if (!tptp_lex(tokens, in, start)) {
		read_error("Unable to parse first-order formula (lexical analysis failed)", start);
		free_tokens(tokens); return false;
	}

	unsigned int index = 0;
	array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
	if (!tptp_interpret(tokens, index, formula, names, variables)) {
		read_error("Unable to parse first-order formula", start);
		for (auto entry : variables) free(entry.key);
		free_tokens(tokens); return false;
	}
	free_tokens(tokens);
	return true;
}

#endif /* FIRST_ORDER_LOGIC_H_ */
