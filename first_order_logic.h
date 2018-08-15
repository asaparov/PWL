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

fol_term EMPTY_FOL_TERM = {fol_term_type::NONE, 0};

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

inline bool new_fol_formula(fol_formula*& new_formula) {
	new_formula = (fol_formula*) malloc(sizeof(fol_formula));
	if (new_formula == NULL) {
		fprintf(stderr, "new_fol_formula ERROR: Out of memory.\n");
		return false;
	}
	return true;
}

constexpr bool visit_constant(unsigned int constant) { return true; }
constexpr bool visit_predicate(unsigned int predicate) { return true; }
constexpr bool visit_variable(unsigned int variable) { return true; }

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula_type& formula) { return true; }

template<typename... Visiter>
inline bool visit(fol_term& term, Visiter&&... visiter) {
	switch (term.type) {
	case fol_term_type::CONSTANT:
		return visit_constant(term.constant, std::forward<Visiter>(visiter)...);
	case fol_term_type::VARIABLE:
		return visit_variable(term.variable, std::forward<Visiter>(visiter)...);
	case fol_term_type::NONE:
		return true;
	}
	fprintf(stderr, "visit ERROR: Unrecognized fol_term_type.\n");
	return false;
}

template<typename... Visiter>
inline bool visit(fol_atom& atom, Visiter&&... visiter) {
	return visit_predicate(atom.predicate, std::forward<Visiter>(visiter)...)
		&& visit(atom.arg1, std::forward<Visiter>(visiter)...)
		&& visit(atom.arg2, std::forward<Visiter>(visiter)...);
}

template<typename... Visiter>
bool visit(fol_formula& formula, Visiter&&... visiter)
{
	switch (formula.type) {
	case fol_formula_type::ATOM:
		return visit(formula.atom, std::forward<Visiter>(visiter)...);
	case fol_formula_type::NOT:
		return visit_operator<fol_formula_type::NOT>(formula, std::forward<Visiter>(visiter)...)
			&& visit(*formula.unary.operand, std::forward<Visiter>(visiter)...);
	case fol_formula_type::AND:
		return visit_operator<fol_formula_type::AND>(formula, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.left, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.right, std::forward<Visiter>(visiter)...);
	case fol_formula_type::OR:
		return visit_operator<fol_formula_type::OR>(formula, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.left, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.right, std::forward<Visiter>(visiter)...);
	case fol_formula_type::IF_THEN:
		return visit_operator<fol_formula_type::IF_THEN>(formula, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.left, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.right, std::forward<Visiter>(visiter)...);
	case fol_formula_type::IFF:
		return visit_operator<fol_formula_type::IFF>(formula, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.left, std::forward<Visiter>(visiter)...)
			&& visit(*formula.binary.right, std::forward<Visiter>(visiter)...);
	case fol_formula_type::FOR_ALL:
		return visit_operator<fol_formula_type::FOR_ALL>(formula, std::forward<Visiter>(visiter)...)
			&& visit_variable(formula.quantifier.variable, std::forward<Visiter>(visiter)...)
			&& visit(*formula.quantifier.operand, std::forward<Visiter>(visiter)...);
	case fol_formula_type::EXISTS:
		return visit_operator<fol_formula_type::EXISTS>(formula, std::forward<Visiter>(visiter)...)
			&& visit_variable(formula.quantifier.variable, std::forward<Visiter>(visiter)...)
			&& visit(*formula.quantifier.operand, std::forward<Visiter>(visiter)...);
	case fol_formula_type::TRUE:
		return visit_true(formula, std::forward<Visiter>(visiter)...);
	case fol_formula_type::FALSE:
		return visit_false(formula, std::forward<Visiter>(visiter)...);
	}
	fprintf(stderr, "visit ERROR: Unrecognized fol_formula_type.\n");
	return false;
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

template<typename... Cloner>
inline bool clone(const fol_term& src, fol_term& dst, Cloner&&... cloner) {
	dst.type = src.type;
	switch (src.type) {
	case fol_term_type::CONSTANT:
		return clone_constant(src.constant, dst.constant, std::forward<Cloner>(cloner)...);
	case fol_term_type::VARIABLE:
		return clone_variable(src.variable, dst.variable, std::forward<Cloner>(cloner)...);
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
	case fol_formula_type::NOT:
		if (!new_fol_formula(dst.unary.operand)) return false;
		if (!clone(*src.unary.operand, *dst.unary.operand, std::forward<Cloner>(cloner)...)) {
			free(dst.unary.operand);
			return false;
		}
		return true;
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
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


/**
 * Functions for easily constructing first-order logic expressions in code.
 */


fol_term make_fol_variable(unsigned int variable) {
	fol_term term;
	term.type = fol_term_type::VARIABLE;
	term.variable = variable;
	return term;
}

fol_term make_fol_constant(unsigned int constant) {
	fol_term term;
	term.type = fol_term_type::CONSTANT;
	term.constant = constant;
	return term;
}

fol_formula* make_fol_atom(
		unsigned int predicate,
		fol_term arg1, fol_term arg2)
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

inline fol_formula* make_fol_atom(unsigned int predicate, fol_term arg1) {
	return make_fol_atom(predicate, arg1, EMPTY_FOL_TERM);
}

inline fol_formula* make_fol_atom(unsigned int predicate) {
	return make_fol_atom(predicate, EMPTY_FOL_TERM, EMPTY_FOL_TERM);
}

#include <tuple>

template<unsigned int... I>
struct static_sequence { };

template<unsigned int N, unsigned int... I>
struct reverse_static_sequence : reverse_static_sequence<N - 1, I..., N - 1> { };

template<unsigned int... I>
struct reverse_static_sequence<0, I...> : static_sequence<I...> { };

template<fol_formula_type Operator, typename... Args>
fol_formula* make_fol_binary_helper(fol_formula* arg, Args&&... args)
{
	if (arg == NULL) return NULL;

	fol_formula* formula;
	if (!new_fol_formula(formula)) return NULL;
	formula->reference_count = 1;
	formula->type = Operator;
	formula->binary.right = arg;
	formula->binary.left = make_fol_and_helper(std::forward<Args>(args)...);
	if (formula->binary.left == NULL) {
		free(formula); free(*arg);
		if (arg->reference_count == 0) free(arg);
		return NULL;
	}
	return formula;
}

template<fol_formula_type Operator, typename Tuple, unsigned int... I>
inline fol_formula* make_fol_binary_helper(Tuple&& args, const static_sequence<I...>& seq)
{
	return make_fol_binary_helper<Operator>(std::get<I>(args)...);
}

template<fol_formula_type Operator, typename... Args>
inline fol_formula* make_fol_binary(Args&&... args)
{
	constexpr auto seq = reverse_static_sequence<sizeof...(Args)>();
	return make_fol_binary_helper<Operator>(std::forward_as_tuple(args...), seq);
}

template<typename... Args>
inline fol_formula* make_fol_and(Args&&... args) {
	return make_fol_binary<fol_formula_type::AND>(std::forward<Args>(args)...);
}

template<typename... Args>
inline fol_formula* make_fol_or(Args&&... args) {
	return make_fol_binary<fol_formula_type::OR>(std::forward<Args>(args)...);
}

template<typename... Args>
inline fol_formula* make_fol_iff(Args&&... args) {
	return make_fol_binary<fol_formula_type::IFF>(std::forward<Args>(args)...);
}

fol_formula* make_fol_if_then(fol_formula* first, fol_formula* second)
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

fol_formula* make_fol_not(fol_formula* operand)
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
fol_formula* make_fol_quantifier(unsigned int variable, fol_formula* operand)
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

inline fol_formula* make_fol_for_all(unsigned int variable, fol_formula* operand) {
	return make_fol_quantifier<fol_formula_type::FOR_ALL>(variable, operand);
}

inline fol_formula* make_fol_exists(unsigned int variable, fol_formula* operand) {
	return make_fol_quantifier<fol_formula_type::EXISTS>(variable, operand);
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
		if (first.variable < second.variable) return -1;
		else if (first.variable > second.variable) return 1;
		else return 0;
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
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
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


/* forward declarations */
bool relabel_variables(fol_formula&, array_map<unsigned int, unsigned int>&);


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
bool operator == (const fol_scope&, const fol_scope&);
int_fast8_t compare(const fol_scope&, const fol_scope&, const canonicalizer&);
void shift_variables(fol_scope&, unsigned int);
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
		const fol_commutative_scope& second,
		const canonicalizer& sorter)
{
	if (first.children.length < second.children.length) return -1;
	else if (first.children.length > second.children.length) return 1;
	else if (first.negated.length < second.negated.length) return -1;
	else if (first.negated.length > second.negated.length) return 1;

	for (unsigned int i = 0; i < first.children.length; i++) {
		int_fast8_t result = compare(first.children[i], second.children[i], sorter);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.negated.length; i++) {
		int_fast8_t result = compare(first.negated[i], second.negated[i], sorter);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const fol_noncommutative_scope& first,
		const fol_noncommutative_scope& second,
		const canonicalizer& sorter)
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
		int_fast8_t result = compare(first.left[i], second.left[i], sorter);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.left_negated.length; i++) {
		int_fast8_t result = compare(first.left_negated[i], second.left_negated[i], sorter);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.right.length; i++) {
		int_fast8_t result = compare(first.right[i], second.right[i], sorter);
		if (result != 0) return result;
	} for (unsigned int i = 0; i < first.right_negated.length; i++) {
		int_fast8_t result = compare(first.right_negated[i], second.right_negated[i], sorter);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const fol_quantifier_scope& first,
		const fol_quantifier_scope& second,
		const canonicalizer& sorter)
{
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
	return compare(*first.operand, *second.operand, sorter);
}

int_fast8_t compare(
		const fol_scope& first,
		const fol_scope& second,
		const canonicalizer& sorter)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;
	switch (first.type) {
	case fol_formula_type::ATOM:
		return compare(first.atom, second.atom, sorter);
	case fol_formula_type::NOT:
		return compare(*first.unary, *second.unary, sorter);
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IFF:
		return compare(first.commutative, second.commutative, sorter);
	case fol_formula_type::IF_THEN:
		return compare(first.noncommutative, second.noncommutative, sorter);
	case fol_formula_type::FOR_ALL:
	case fol_formula_type::EXISTS:
		return compare(first.quantifier, second.quantifier, sorter);
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized fol_formula_type when comparing fol_scopes.\n");
	exit(EXIT_FAILURE);
}

inline bool less_than(
		const fol_scope& first,
		const fol_scope& second,
		const canonicalizer& sorter)
{
	return compare(first, second, sorter) < 0;
}

inline bool less_than(
		const fol_scope* first,
		const fol_scope* second,
		const canonicalizer& sorter)
{
	return less_than(*first, *second, sorter);
}

inline void shift_variables(fol_term& term, unsigned int removed_variable) {
	switch (term.type) {
	case fol_term_type::VARIABLE:
		if (term.variable > removed_variable)
			term.variable--;
		return;
	case fol_term_type::CONSTANT:
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
	fol_formula* left = first;

	for (unsigned int i = 0; i < scope_length; i++) {
		fol_formula* right = scope_to_formula<Negated>(scope[i]);
		if (right == NULL) {
			free(*left); if (left->reference_count == 0) free(left);
			return NULL;
		}

		fol_formula* new_formula;
		if (!new_fol_formula(new_formula)) {
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
		const fol_scope* scope, unsigned int scope_length,
		const fol_scope* negated, unsigned int negated_length)
{
	if (scope_length == 0) {
		fol_formula* first = scope_to_formula<true>(negated[0]);
		if (first == NULL) return NULL;
		return scope_to_formula<ScopeType, true>(negated + 1, negated_length - 1, first);
	} else {
		fol_formula* first = scope_to_formula<false>(scope[0]);
		if (first == NULL) return NULL;
		first = scope_to_formula<ScopeType, false>(scope + 1, scope_length - 1, first);
		if (first == NULL) return NULL;
		return scope_to_formula<ScopeType, true>(negated, negated_length, first);
	}
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
		auto result = compare(subscope, scope[index], canonicalizer());
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
		auto result = compare(first[i], second[j], canonicalizer());
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

		auto result = compare(first[i], second[j], canonicalizer());
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
		auto result = compare(first[i], second[j], canonicalizer());
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

		auto result = compare(first[i], second[j], canonicalizer());
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
	} else if (scope.commutative.children.last().type == fol_formula_type::FALSE) {
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

inline bool make_negated_iff(fol_scope& out, fol_scope& operand)
{
	if (operand.type == fol_formula_type::TRUE) {
		free(operand);
		return init(out, fol_formula_type::FALSE);
	} else if (operand.type == fol_formula_type::FALSE) {
		free(operand);
		return init(out, fol_formula_type::TRUE);
	} else {
		if (!init(out, fol_formula_type::IFF)) {
			free(operand); return false;
		}
		if (operand.type == fol_formula_type::NOT) {
			move(*operand.unary, out.commutative.negated[0]);
			out.commutative.negated.length++;
			move_variables(operand.variables, out.variables);
			free(operand.unary);
			free(operand.variables);
		} else {
			move(operand, out.commutative.children[0]);
			out.commutative.children.length++;
			move_variables(operand.variables, out.variables);
		}
		if (!init(out.commutative.children[out.commutative.children.length], fol_formula_type::FALSE)) {
			free(out); return false;
		}
		out.commutative.children.length++;
		return true;
	}
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

template<fol_formula_type Operator>
bool canonicalize_commutative_scope(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	static_assert(Operator == fol_formula_type::AND
			   || Operator == fol_formula_type::OR
			   || Operator == fol_formula_type::IFF,
			"Operator is not a commutative operator.");

	if (!canonicalize_scope(*src.left, out, variable_map)) return false;

	fol_scope& right = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (out.type == fol_formula_type::FALSE) {
		if (Operator == fol_formula_type::AND) {
			return true;
		} else if (Operator == fol_formula_type::OR) {
			free(out);
			return canonicalize_scope(*src.right, out, variable_map);
		} else if (Operator == fol_formula_type::IFF) {
			free(out);
			if (!canonicalize_scope(*src.right, right, variable_map))
				return false;
			return make_negated_iff(out, right);
		}
	} else if (out.type == fol_formula_type::TRUE) {
		if (Operator == fol_formula_type::OR) {
			return true;
		} else {
			free(out);
			return canonicalize_scope(*src.right, out, variable_map);
		}
	}

	if (!canonicalize_scope(*src.right, right, variable_map))
		return false;

	if (out == right) {
		/* consider the case where `*src.left` and `*src.right` are identical */
		free(right);
		if (Operator == fol_formula_type::IFF) {
			free(out);
			return init(out, fol_formula_type::TRUE);
		}
	} else if (right.type == fol_formula_type::FALSE) {
		free(right);
		if (Operator == fol_formula_type::AND) {
			free(out);
			return init(out, fol_formula_type::FALSE);
		} else if (Operator == fol_formula_type::IFF) {
			if (!make_negated_iff(right, out)) return false;
			move(right, out);
			return true;
		}
	} else if (right.type == fol_formula_type::TRUE) {
		free(right);
		if (Operator == fol_formula_type::OR) {
			free(out);
			return init(out, fol_formula_type::TRUE);
		}
	} else if (are_negations(out, right)) {
		/* we have the case `A op ~A` */
		free(right); free(out);
		if (Operator == fol_formula_type::AND || Operator == fol_formula_type::IFF)
			return init(out, fol_formula_type::FALSE);
		else if (Operator == fol_formula_type::OR)
			return init(out, fol_formula_type::TRUE);
	} else if (out.type == Operator) {
		if (right.type == Operator) {
			bool found_negation;
			merge_scopes<Operator>(right.commutative.children, out.commutative.children,
					right.commutative.negated, out.commutative.negated, found_negation);
			free(right.commutative.children);
			free(right.commutative.negated);
			if (Operator == fol_formula_type::IFF)
				recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
			else move_variables(right.variables, out.variables);
			free(right.variables);
			if (Operator == fol_formula_type::IFF && out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
				free(out);
				return init(out, found_negation ? fol_formula_type::FALSE : fol_formula_type::TRUE);
			} else if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(out);
					return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(out);
					return init(out, fol_formula_type::TRUE);
				} else if (Operator == fol_formula_type::IFF) {
					recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
					return negate_iff(out);
				}
			}
		} else {
			bool found_negation;
			if (!add_to_scope<Operator>(right, out.commutative.children, out.commutative.negated, out.variables, found_negation)) {
				free(out); free(right); return false;
			} else if (Operator == fol_formula_type::IFF && out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
				free(out);
				return init(out, found_negation ? fol_formula_type::FALSE : fol_formula_type::TRUE);
			} else if (Operator == fol_formula_type::IFF && out.commutative.children.length == 1 && out.commutative.negated.length == 0) {
				move(out.commutative.children[0], right);
				out.commutative.children.clear();
				free(out); move(right, out);
				if (found_negation && !negate_scope(out)) {
					free(out); return false;
				}
			} else if (Operator == fol_formula_type::IFF && out.commutative.children.length == 0 && out.commutative.negated.length == 1) {
				move(out.commutative.negated[0], right);
				out.commutative.negated.clear();
				free(out);
				if (!found_negation && !negate_scope(right)) {
					free(right); return false;
				}
				move(right, out);
			} else if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(out); return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(out); return init(out, fol_formula_type::TRUE);
				} else if (Operator == fol_formula_type::IFF) {
					return negate_iff(out);
				}
			}
		}
	} else {
		if (right.type == Operator) {
			bool found_negation;
			if (!add_to_scope<Operator>(out, right.commutative.children, right.commutative.negated, right.variables, found_negation)) {
				free(out); free(right); return false;
			} else if (Operator == fol_formula_type::IFF && right.commutative.children.length == 0 && right.commutative.negated.length == 0) {
				free(right);
				return init(out, found_negation ? fol_formula_type::FALSE : fol_formula_type::TRUE);
			} else if (Operator == fol_formula_type::IFF && right.commutative.children.length == 1 && right.commutative.negated.length == 0) {
				move(right.commutative.children[0], out);
				right.commutative.children.clear();
				free(right);
				if (found_negation && !negate_scope(out)) {
					free(out); return false;
				}
			} else if (Operator == fol_formula_type::IFF && right.commutative.children.length == 0 && right.commutative.negated.length == 1) {
				move(right.commutative.negated[0], out);
				right.commutative.negated.clear();
				free(right);
				if (!found_negation && !negate_scope(out)) {
					free(out); return false;
				}
			} else if (found_negation) {
				if (Operator == fol_formula_type::AND) {
					free(right); return init(out, fol_formula_type::FALSE);
				} else if (Operator == fol_formula_type::OR) {
					free(right); return init(out, fol_formula_type::TRUE);
				} else if (Operator == fol_formula_type::IFF) {
					move(right, out);
					return negate_iff(out);
				}
			} else {
				move(right, out);
			}
		} else {
			fol_scope& both = *((fol_scope*) alloca(sizeof(fol_scope)));
			if (!init(both, Operator, out.variables)) {
				free(out); free(right);
				return false;
			}
			move_variables(right.variables, both.variables);

			if (out.type == fol_formula_type::NOT) {
				move(*out.unary, both.commutative.negated[both.commutative.negated.length]);
				both.commutative.negated.length++;
				free(out.variables);
				free(out.unary);
			} else if (out.type == fol_formula_type::IFF && out.commutative.children.length > 0
					&& out.commutative.children.last().type == fol_formula_type::FALSE)
			{
				free(out.commutative.children.last());
				out.commutative.children.length--;
				move(out, both.commutative.negated[both.commutative.negated.length]);
				both.commutative.negated.length++;
			} else {
				move(out, both.commutative.children[both.commutative.children.length]);
				both.commutative.children.length++;
			}
			if (right.type == fol_formula_type::NOT) {
				move(*right.unary, both.commutative.negated[both.commutative.negated.length]);
				both.commutative.negated.length++;
				free(right.variables);
				free(right.unary);
			} else if (right.type == fol_formula_type::IFF && right.commutative.children.length > 0
					&& right.commutative.children.last().type == fol_formula_type::FALSE)
			{
				free(right.commutative.children.last());
				right.commutative.children.length--;
				move(right, both.commutative.negated[both.commutative.negated.length]);
				both.commutative.negated.length++;
			} else {
				move(right, both.commutative.children[both.commutative.children.length]);
				both.commutative.children.length++;
			}
			if (both.commutative.children.length > 1)
				insertion_sort(both.commutative.children, canonicalizer());
			if (both.commutative.negated.length > 1)
				insertion_sort(both.commutative.negated, canonicalizer());
			move(both, out);
		}
	}
	return true;
}

bool canonicalize_conditional_scope(
		const fol_binary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	fol_scope& left = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (!canonicalize_scope(*src.left, left, variable_map)) return false;

	if (left.type == fol_formula_type::FALSE) {
		free(left);
		return init(out, fol_formula_type::TRUE);
	} else if (left.type == fol_formula_type::TRUE) {
		free(left);
		return canonicalize_scope(*src.right, out, variable_map);
	}

	if (!canonicalize_scope(*src.right, out, variable_map))
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

template<fol_formula_type QuantifierType>
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
			|| !canonicalize_scope(*src.operand, *operand, variable_map))
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

inline bool canonicalize_negation_scope(
		const fol_unary_formula& src, fol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map)
{
	if (!canonicalize_scope(*src.operand, out, variable_map))
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
	case fol_formula_type::FOR_ALL:
		return canonicalize_quantifier_scope<fol_formula_type::FOR_ALL>(src.quantifier, out, variable_map);
	case fol_formula_type::EXISTS:
		return canonicalize_quantifier_scope<fol_formula_type::EXISTS>(src.quantifier, out, variable_map);
	case fol_formula_type::NOT:
		return canonicalize_negation_scope(src.unary, out, variable_map);
	case fol_formula_type::ATOM:
		return canonicalize_atom(src.atom, out, variable_map);
	case fol_formula_type::TRUE:
		return init(out, fol_formula_type::TRUE);
	case fol_formula_type::FALSE:
		return init(out, fol_formula_type::FALSE);
	}
	fprintf(stderr, "canonicalize_scope ERROR: Unrecognized fol_formula_type.\n");
	return false;
}

fol_formula* canonicalize(const fol_formula& src)
{
	array_map<unsigned int, unsigned int> variable_map(16);
	fol_scope& scope = *((fol_scope*) alloca(sizeof(fol_scope)));
	if (!canonicalize_scope(src, scope, variable_map))
		return NULL;
	fol_formula* canonicalized = scope_to_formula(scope);
	free(scope);
	return canonicalized;
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

template<fol_formula_type OperatorType>
bool tptp_interpret_binary_sequence(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	fol_formula* left)
{
	while (true) {
		fol_formula* next;
		if (!new_fol_formula(next)) {
			free(*left); free(left); return false;
		} else if (!tptp_interpret_unary_formula(tokens, index, *next, names, variables)) {
			free(*left); free(left); free(next); return false;
		}

		if (index < tokens.length && tokens[index].type == tptp_operator_type<OperatorType>::type) {
			index++;
			fol_formula* parent;
			if (!new_fol_formula(parent)) {
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
		if (!tptp_interpret_binary_formula<fol_formula_type::IFF>(tokens, index, formula, names, variables, left))
			return false;
	} else {
		move(*left, formula); free(left);
	}
	return true;
}

#endif /* FIRST_ORDER_LOGIC_H_ */
