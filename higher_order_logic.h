/**
 * higher_order_logic.h
 *
 *  Created on: Dec 31, 2018
 *      Author: asaparov
 */

#ifndef HIGHER_ORDER_LOGIC_H_
#define HIGHER_ORDER_LOGIC_H_

#include <core/lex.h>
#include <math/multiset.h>
#include <cstdint>

#include "array_view.h"

using namespace core;

constexpr unsigned int ANY_VARIABLE = UINT_MAX;


/* forward declarations */
struct hol_term;
bool operator == (const hol_term&, const hol_term&);
bool operator != (const hol_term&, const hol_term&);

enum class hol_term_type {
	VARIABLE = 1,
	CONSTANT,
	PARAMETER,

	UNARY_APPLICATION,
	BINARY_APPLICATION,

	AND,
	OR,
	IF_THEN,
	EQUALS,
	IFF, /* this is only used in canonicalization */
	NOT,

	FOR_ALL,
	EXISTS,
	LAMBDA,

	INTEGER,
	STRING,
	UINT_LIST,

	ANY, /* represents a set of all logical forms that contain a subtree within a specified set, and do not contain trees in other specified sets */
	ANY_RIGHT, /* represents a set of all logical forms that contain a subtree within a specified set that is on the right-leaning branch from the root */
	ANY_ARRAY, /* represents the set of all logical forms with an "array" type (e.g. `AND`, `OR`, `IFF`)  */
	ANY_CONSTANT, /* represents a set of logical forms of type `CONSTANT` */
	ANY_QUANTIFIER, /* represents a node that is either `FOR_ALL`, `EXISTS`, or `LAMBDA` with any variable */

	TRUE,
	FALSE /* our canonicalization code assumes that FALSE is the last element of this enum */
};

template<typename Stream>
inline bool print_subscript(unsigned int number, Stream& out) {
	static const char* subscripts[] = { "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉" };
	if (number == 0)
		return print(subscripts[0], out);
	array<uint8_t> digits(32);
	while (number > 0) {
		if (!digits.add(number % 10)) return false;
		number /= 10;
	}
	for (unsigned int i = digits.length; i > 0; i--)
		if (!print(subscripts[digits[i - 1]], out)) return false;
	return true;
}

enum class hol_term_syntax {
	TPTP,
	CLASSIC
};

template<hol_term_syntax Syntax, typename Stream>
inline bool print_variable(unsigned int variable, Stream& out) {
	if (variable == ANY_VARIABLE)
		return print('*', out);
	switch (Syntax) {
	case hol_term_syntax::TPTP:
		return print('$', out) + print(variable, out);
	case hol_term_syntax::CLASSIC:
		return print('x', out) && print_subscript(variable, out);
	}
	fprintf(stderr, "print_variable ERROR: Unrecognized hol_term_syntax.\n");
	return false;
}

template<hol_term_syntax Syntax, typename Stream>
inline bool print_parameter(unsigned int parameter, Stream& out) {
	switch (Syntax) {
	case hol_term_syntax::TPTP:
		return print('#', out) + print(parameter, out);
	case hol_term_syntax::CLASSIC:
		return print('a', out) && print_subscript(parameter, out);
	}
	fprintf(stderr, "print_parameter ERROR: Unrecognized hol_term_syntax.\n");
	return false;
}

struct hol_unary_term {
	hol_term* operand;

	static inline unsigned int hash(const hol_unary_term& key);
	static inline void move(const hol_unary_term& src, hol_unary_term& dst);
	static inline void free(hol_unary_term& term);
};

struct hol_binary_term {
	hol_term* left;
	hol_term* right;

	static inline unsigned int hash(const hol_binary_term& key);
	static inline void move(const hol_binary_term& src, hol_binary_term& dst);
	static inline void free(hol_binary_term& term);
};

struct hol_ternary_term {
	hol_term* first;
	hol_term* second;
	hol_term* third;

	static inline unsigned int hash(const hol_ternary_term& key);
	static inline void move(const hol_ternary_term& src, hol_ternary_term& dst);
	static inline void free(hol_ternary_term& term);
};

struct hol_array_term {
	hol_term** operands;
	unsigned int length;

	static inline unsigned int hash(const hol_array_term& key);
	static inline void move(const hol_array_term& src, hol_array_term& dst);
	static inline void free(hol_array_term& term);

private:
	inline bool init_helper(const hol_array_term& src);

	friend struct hol_term;
};

struct hol_quantifier {
	unsigned int variable;
	hol_term* operand;

	static inline unsigned int hash(const hol_quantifier& key);
	static inline void move(const hol_quantifier& src, hol_quantifier& dst);
	static inline void free(hol_quantifier& term);
};

struct hol_any {
	hol_term* included;
	hol_term** excluded_trees;
	unsigned int excluded_tree_count;

	static inline unsigned int hash(const hol_any& key);
	static inline void move(const hol_any& src, hol_any& dst);
	static inline void free(hol_any& term);
};

struct hol_any_array {
	hol_term_type oper;
	hol_term* all;
	hol_array_term any; /* this could overlap with `left` or `right` */
	hol_array_term left;
	hol_array_term right;

	static inline unsigned int hash(const hol_any_array& key);
	static inline void move(const hol_any_array& src, hol_any_array& dst);
	static inline void free(hol_any_array& term);
};

struct hol_any_constant {
	unsigned int* constants;
	unsigned int length;

	static inline unsigned int hash(const hol_any_constant& key);
	static inline void move(const hol_any_constant& src, hol_any_constant& dst);
	static inline void free(hol_any_constant& term);
};

enum class hol_quantifier_type {
	FOR_ALL = (int) hol_term_type::FOR_ALL,
	EXISTS = (int) hol_term_type::EXISTS,
	LAMBDA = (int) hol_term_type::LAMBDA,

	FOR_ALL_OR_EXISTS,
	FOR_ALL_OR_LAMBDA,
	EXISTS_OR_LAMBDA,
	ANY
};

struct hol_any_quantifier {
	hol_quantifier_type quantifier;
	hol_term* operand;

	static inline unsigned int hash(const hol_any_quantifier& key);
	static inline void move(const hol_any_quantifier& src, hol_any_quantifier& dst);
	static inline void free(hol_any_quantifier& term);
};

struct hol_term
{
	typedef hol_term_type Type;
	typedef hol_term Term;
	typedef hol_term_type TermType;

	hol_term_type type;
	unsigned int reference_count;
	union {
		unsigned int variable;
		unsigned int constant;
		unsigned int parameter;
		int integer;
		string str;
		sequence uint_list;
		hol_unary_term unary;
		hol_binary_term binary;
		hol_ternary_term ternary;
		hol_array_term array;
		hol_quantifier quantifier;
		hol_any any;
		hol_any_array any_array;
		hol_any_constant any_constant;
		hol_any_quantifier any_quantifier;
	};

	hol_term() : type(hol_term_type::ANY), reference_count(1) {
		any.included = nullptr;
		any.excluded_trees = nullptr;
		any.excluded_tree_count = 0;
	}

	hol_term(hol_term_type type) : type(type), reference_count(1) {
#if !defined(NDEBUG)
		if (type != hol_term_type::TRUE && type != hol_term_type::FALSE && type != hol_term_type::ANY)
			fprintf(stderr, "hol_term WARNING: This constructor should only be used if `type` is `TRUE`, `FALSE`, or `ANY`.\n");
#endif
	}

	hol_term(const hol_term& src) : type(src.type), reference_count(1) {
		init_helper(src);
	}

	~hol_term() { free_helper(); }

	inline void operator = (const hol_term& src) {
		type = src.type;
		reference_count = 1;
		init_helper(src);
	}

	static hol_term* new_variable(unsigned int variable);
	static hol_term* new_constant(unsigned int constant);
	static hol_term* new_parameter(unsigned int parameter);
	static hol_term* new_int(int integer);
	static hol_term* new_string(const string& str);
	static hol_term* new_uint_list(const sequence& list);
	static hol_term* new_atom(unsigned int predicate, hol_term* arg1, hol_term* arg2);
	static inline hol_term* new_atom(unsigned int predicate, hol_term* arg1);
	static inline hol_term* new_atom(unsigned int predicate);
	static hol_term* new_true();
	static hol_term* new_false();
	static hol_term* new_apply(hol_term* function, hol_term* arg);
	static hol_term* new_apply(hol_term* function, hol_term* arg1, hol_term* arg2);
	template<typename... Args> static inline hol_term* new_and(Args&&... args);
	template<typename... Args> static inline hol_term* new_or(Args&&... args);
	template<typename... Args> static inline hol_term* new_iff(Args&&... args);
	static inline hol_term* new_equals(hol_term* first, hol_term* second);
	template<typename Array, typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type** = nullptr> static inline hol_term* new_and(Array& args);
	template<typename Array, typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type** = nullptr> static inline hol_term* new_or(Array& args);
	static hol_term* new_if_then(hol_term* first, hol_term* second);
	static hol_term* new_not(hol_term* operand);
	static inline hol_term* new_for_all(unsigned int variable, hol_term* operand);
	static inline hol_term* new_exists(unsigned int variable, hol_term* operand);
	static inline hol_term* new_lambda(unsigned int variable, hol_term* operand);
	static hol_term* new_any(hol_term* included);
	static hol_term* new_any(hol_term* included, hol_term** excluded_trees, unsigned int excluded_tree_count);
	static hol_term* new_any_right(hol_term* included);
	static hol_term* new_any_right(hol_term* included, hol_term** excluded_trees, unsigned int excluded_tree_count);
	template<typename AnyArray, typename LeftArray, typename RightArray,
		typename std::enable_if<has_index_operator<AnyArray, hol_term*>::value>::type** = nullptr,
		typename std::enable_if<has_index_operator<LeftArray, hol_term*>::value>::type** = nullptr,
		typename std::enable_if<has_index_operator<RightArray, hol_term*>::value>::type** = nullptr>
	static hol_term* new_any_array(hol_term_type oper, hol_term* all, const AnyArray& any, const LeftArray& left, const RightArray& right);
	template<typename... Args> static hol_term* new_any_constant(unsigned int constant, Args&&... args);
	template<typename Array, typename std::enable_if<has_index_operator<Array, unsigned int>::value>::type** = nullptr> static hol_term* new_any_constant(const Array& args);
	static hol_term* new_any_quantifier(hol_quantifier_type quantifier_type, hol_term* operand);

	static inline unsigned int hash(const hol_term& key);
	static inline bool is_empty(const hol_term& key);
	static inline void set_empty(hol_term& key);
	static inline void move(const hol_term& src, hol_term& dst);
	static inline void swap(hol_term& first, hol_term& second);
	static inline void free(hol_term& term) { term.free_helper(); }

private:
	void free_helper();

	inline bool init_helper(const hol_term& src) {
		switch (src.type) {
		case hol_term_type::VARIABLE:
			variable = src.variable; return true;
		case hol_term_type::CONSTANT:
			constant = src.constant; return true;
		case hol_term_type::PARAMETER:
			parameter = src.parameter; return true;
		case hol_term_type::INTEGER:
			integer = src.integer; return true;
		case hol_term_type::STRING:
			str = src.str; return true;
		case hol_term_type::UINT_LIST:
			uint_list = src.uint_list; return true;
		case hol_term_type::NOT:
			unary.operand = src.unary.operand;
			unary.operand->reference_count++;
			return true;
		case hol_term_type::IF_THEN:
		case hol_term_type::EQUALS:
		case hol_term_type::UNARY_APPLICATION:
			binary.left = src.binary.left;
			binary.right = src.binary.right;
			binary.left->reference_count++;
			binary.right->reference_count++;
			return true;
		case hol_term_type::BINARY_APPLICATION:
			ternary.first = src.ternary.first;
			ternary.second = src.ternary.second;
			ternary.third = src.ternary.third;
			ternary.first->reference_count++;
			ternary.second->reference_count++;
			ternary.third->reference_count++;
			return true;
		case hol_term_type::AND:
		case hol_term_type::OR:
		case hol_term_type::IFF:
			return array.init_helper(src.array);
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			quantifier.variable = src.quantifier.variable;
			quantifier.operand = src.quantifier.operand;
			quantifier.operand->reference_count++;
			return true;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
			any.included = src.any.included;
			if (any.included != nullptr)
				any.included->reference_count++;
			any.excluded_tree_count = src.any.excluded_tree_count;
			if (src.any.excluded_tree_count == 0) {
				any.excluded_trees = nullptr;
			} else {
				any.excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * src.any.excluded_tree_count);
				if (any.excluded_trees == nullptr) {
					fprintf(stderr, "hol_term.init_helper ERROR: Insufficient memory for `any.excluded_trees`.\n");
					if (any.included != nullptr)
						core::free(any.included);
					return false;
				}
			}
			for (unsigned int i = 0; i < src.any.excluded_tree_count; i++) {
				any.excluded_trees[i] = src.any.excluded_trees[i];
				any.excluded_trees[i]->reference_count++;
			}
			return true;
		case hol_term_type::ANY_ARRAY:
			any_array.oper = src.any_array.oper;
			if (!any_array.any.init_helper(src.any_array.any)) {
				return false;
			} else if (!any_array.left.init_helper(src.any_array.left)) {
				core::free(any_array.any);
				return false;
			} else if (!any_array.right.init_helper(src.any_array.right)) {
				core::free(any_array.any);
				core::free(any_array.left);
				return false;
			}
			any_array.all = src.any_array.all;
			any_array.all->reference_count++;
			return true;
		case hol_term_type::ANY_CONSTANT:
			any_constant.length = src.any_constant.length;
			any_constant.constants = (unsigned int*) malloc(sizeof(unsigned int) * src.any_constant.length);
			if (any_constant.constants == nullptr) {
				fprintf(stderr, "hol_term.init_helper ERROR: Insufficient memory for `any_constant.excluded_trees`.\n");
				return false;
			}
			for (unsigned int i = 0; i < any_constant.length; i++)
				any_constant.constants[i] = src.any_constant.constants[i];
			return true;
		case hol_term_type::ANY_QUANTIFIER:
			any_quantifier.quantifier = src.any_quantifier.quantifier;
			any_quantifier.operand = src.any_quantifier.operand;
			any_quantifier.operand->reference_count++;
			return true;
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
			return true;
		}
		fprintf(stderr, "hol_term.init_helper ERROR: Unrecognized hol_term_type.\n");
		return false;
	}

	friend bool init(hol_term&, const hol_term&);
};

thread_local hol_term HOL_TRUE(hol_term_type::TRUE);
thread_local hol_term HOL_FALSE(hol_term_type::FALSE);
thread_local hol_term HOL_ANY(hol_term_type::ANY);

inline bool init(hol_term& term, const hol_term& src) {
	term.type = src.type;
	term.reference_count = 1;
	return term.init_helper(src);
}

inline void initialize_any(hol_term& term) {
	term.type = hol_term_type::ANY;
	term.any.included = nullptr;
	term.any.excluded_trees = nullptr;
	term.any.excluded_tree_count = 0;
	term.reference_count = 1;
}

inline bool hol_array_term::init_helper(const hol_array_term& src)
{
	length = src.length;
	if (src.length == 0) {
		operands = nullptr;
		return true;
	}

	operands = (hol_term**) malloc(sizeof(hol_term*) * src.length);
	if (operands == NULL) {
		fprintf(stderr, "hol_array_term.init_helper ERROR: Insufficient memory for `operands`.\n");
		return false;
	}
	for (unsigned int i = 0; i < src.length; i++) {
		operands[i] = src.operands[i];
		operands[i]->reference_count++;
	}
	return true;
}

inline bool operator == (const hol_unary_term& first, const hol_unary_term& second) {
	return first.operand == second.operand
		|| *first.operand == *second.operand;
}

inline bool operator == (const hol_binary_term& first, const hol_binary_term& second) {
	return (first.left == second.left || *first.left == *second.left)
		&& (first.right == second.right || *first.right == *second.right);
}

inline bool operator == (const hol_ternary_term& first, const hol_ternary_term& second) {
	return (first.first == second.first || *first.first == *second.first)
		&& (first.second == second.second || *first.second == *second.second)
		&& (first.third == second.third || *first.third == *second.third);
}

inline bool operator == (const hol_array_term& first, const hol_array_term& second) {
	if (first.length != second.length) return false;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.operands[i] != second.operands[i] && *first.operands[i] != *second.operands[i]) return false;
	return true;
}

inline bool operator == (const hol_quantifier& first, const hol_quantifier& second) {
	return first.variable == second.variable
		&& (first.operand == second.operand || *first.operand == *second.operand);
}

inline bool operator == (const hol_any& first, const hol_any& second) {
	if (first.included == nullptr) {
		if (second.included != nullptr) return false;
	} else if (second.included == nullptr) {
		return false;
	} else if (first.included != second.included && *first.included != *second.included) {
		return false;
	}

	if (first.excluded_tree_count != second.excluded_tree_count) return false;
	for (unsigned int i = 0; i < first.excluded_tree_count; i++) {
		if (first.excluded_trees[i] != second.excluded_trees[i] && *first.excluded_trees[i] != *second.excluded_trees[i])
			return false;
	}
	return true;
}

inline bool operator == (const hol_any_array& first, const hol_any_array& second) {
	return (first.oper == second.oper)
		&& (first.all == second.all || *first.all == *second.all)
		&& (first.any == second.any)
		&& (first.left == second.left)
		&& (first.right == second.right);
}

inline bool operator == (const hol_any_constant& first, const hol_any_constant& second) {
	if (first.length != second.length) return false;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.constants[i] != second.constants[i]) return false;
	return true;
}

inline bool operator == (const hol_any_quantifier& first, const hol_any_quantifier& second) {
	return (first.quantifier == second.quantifier)
		&& (first.operand == second.operand || *first.operand == *second.operand);
}

bool operator == (const hol_term& first, const hol_term& second)
{
	if (hol_term::is_empty(first)) return false;
	if (first.type != second.type) return false;
	switch (first.type) {
	case hol_term_type::VARIABLE:
		return first.variable == second.variable;
	case hol_term_type::CONSTANT:
		return first.constant == second.constant;
	case hol_term_type::PARAMETER:
		return first.parameter == second.parameter;
	case hol_term_type::INTEGER:
		return first.integer == second.integer;
	case hol_term_type::STRING:
		return first.str == second.str;
	case hol_term_type::UINT_LIST:
		return first.uint_list == second.uint_list;
	case hol_term_type::NOT:
		return first.unary == second.unary;
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return first.binary == second.binary;
	case hol_term_type::BINARY_APPLICATION:
		return first.ternary == second.ternary;
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return first.array == second.array;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return first.quantifier == second.quantifier;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		return first.any == second.any;
	case hol_term_type::ANY_ARRAY:
		return first.any_array == second.any_array;
	case hol_term_type::ANY_CONSTANT:
		return first.any_constant == second.any_constant;
	case hol_term_type::ANY_QUANTIFIER:
		return first.any_quantifier == second.any_quantifier;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline bool operator != (const hol_term& first, const hol_term& second) {
	return !(first == second);
}

inline bool hol_term::is_empty(const hol_term& key) {
	return key.reference_count == 0;
}

inline void hol_term::set_empty(hol_term& key) {
	key.reference_count = 0;
}

inline unsigned int hol_unary_term::hash(const hol_unary_term& key) {
	return hol_term::hash(*key.operand);
}

inline unsigned int hol_binary_term::hash(const hol_binary_term& key) {
	return hol_term::hash(*key.left) + hol_term::hash(*key.right) * 131071;
}

inline unsigned int hol_ternary_term::hash(const hol_ternary_term& key) {
	return hol_term::hash(*key.first) + hol_term::hash(*key.second) * 127 + hol_term::hash(*key.third) * 524287;
}

inline unsigned int hol_array_term::hash(const hol_array_term& key) {
	unsigned int hash_value = default_hash(key.length);
	for (unsigned int i = 0; i < key.length; i++)
		hash_value ^= hol_term::hash(*key.operands[i]);
	return hash_value;
}

inline unsigned int hol_quantifier::hash(const hol_quantifier& key) {
	return default_hash(key.variable) ^ hol_term::hash(*key.operand);
}

inline unsigned int hol_any::hash(const hol_any& key) {
	unsigned int hash_value = default_hash(key.excluded_tree_count);
	if (key.included != nullptr)
		hash_value ^= hol_term::hash(*key.included);
	for (unsigned int i = 0; i < key.excluded_tree_count; i++)
		hash_value ^= hol_term::hash(*key.excluded_trees[i]);
	return hash_value;
}

inline unsigned int hol_any_array::hash(const hol_any_array& key) {
	return default_hash(key.oper)
		 ^ hol_term::hash(*key.all)
		 ^ hol_array_term::hash(key.any)
		 ^ hol_array_term::hash(key.left)
		 ^ hol_array_term::hash(key.right);
}

inline unsigned int hol_any_constant::hash(const hol_any_constant& key) {
	return default_hash(key.constants, key.length);
}

inline unsigned int hol_any_quantifier::hash(const hol_any_quantifier& key) {
	return default_hash(key.quantifier) ^ hol_term::hash(*key.operand);
}

inline unsigned int hol_term::hash(const hol_term& key) {
	/* TODO: precompute these and store them in a table for faster access */
	unsigned int type_hash = default_hash<hol_term_type, 571290832>(key.type);
	switch (key.type) {
	case hol_term_type::VARIABLE:
		return type_hash ^ default_hash(key.variable);
	case hol_term_type::CONSTANT:
		return type_hash ^ default_hash(key.constant);
	case hol_term_type::PARAMETER:
		return type_hash ^ default_hash(key.parameter);
	case hol_term_type::INTEGER:
		return type_hash ^ default_hash(key.integer);
	case hol_term_type::STRING:
		return type_hash ^ string::hash(key.str);
	case hol_term_type::UINT_LIST:
		return type_hash ^ sequence::hash(key.uint_list);
	case hol_term_type::NOT:
		return type_hash ^ hol_unary_term::hash(key.unary);
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return type_hash ^ hol_binary_term::hash(key.binary);
	case hol_term_type::BINARY_APPLICATION:
		return type_hash ^ hol_ternary_term::hash(key.ternary);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return type_hash ^ hol_array_term::hash(key.array);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return type_hash ^ hol_quantifier::hash(key.quantifier);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		return type_hash ^ hol_any::hash(key.any);
	case hol_term_type::ANY_ARRAY:
		return type_hash ^ hol_any_array::hash(key.any_array);
	case hol_term_type::ANY_CONSTANT:
		return type_hash ^ hol_any_constant::hash(key.any_constant);
	case hol_term_type::ANY_QUANTIFIER:
		return type_hash ^ hol_any_quantifier::hash(key.any_quantifier);
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return type_hash;
	}
	fprintf(stderr, "hol_term.hash ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline void hol_term::move(const hol_term& src, hol_term& dst) {
	dst.type = src.type;
	dst.reference_count = src.reference_count;
	switch (src.type) {
	case hol_term_type::VARIABLE:
		dst.variable = src.variable; return;
	case hol_term_type::CONSTANT:
		dst.constant = src.constant; return;
	case hol_term_type::PARAMETER:
		dst.parameter = src.parameter; return;
	case hol_term_type::INTEGER:
		dst.integer = src.integer; return;
	case hol_term_type::STRING:
		string::move(src.str, dst.str); return;
	case hol_term_type::UINT_LIST:
		sequence::move(src.uint_list, dst.uint_list); return;
	case hol_term_type::NOT:
		core::move(src.unary, dst.unary); return;
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		core::move(src.binary, dst.binary); return;
	case hol_term_type::BINARY_APPLICATION:
		core::move(src.ternary, dst.ternary); return;
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		core::move(src.array, dst.array); return;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		core::move(src.quantifier, dst.quantifier); return;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		core::move(src.any, dst.any); return;
	case hol_term_type::ANY_ARRAY:
		core::move(src.any_array, dst.any_array); return;
	case hol_term_type::ANY_CONSTANT:
		core::move(src.any_constant, dst.any_constant); return;
	case hol_term_type::ANY_QUANTIFIER:
		core::move(src.any_quantifier, dst.any_quantifier); return;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return;
	}
	fprintf(stderr, "hol_term.move ERROR: Unrecognized hol_term_type.\n");
}

inline void hol_term::swap(hol_term& first, hol_term& second) {
	char* first_data = (char*) &first;
	char* second_data = (char*) &second;
	for (unsigned int i = 0; i < sizeof(hol_term); i++)
		core::swap(first_data[i], second_data[i]);
}

inline void hol_unary_term::move(const hol_unary_term& src, hol_unary_term& dst) {
	dst.operand = src.operand;
}

inline void hol_binary_term::move(const hol_binary_term& src, hol_binary_term& dst) {
	dst.left = src.left;
	dst.right = src.right;
}

inline void hol_ternary_term::move(const hol_ternary_term& src, hol_ternary_term& dst) {
	dst.first = src.first;
	dst.second = src.second;
	dst.third = src.third;
}

inline void hol_array_term::move(const hol_array_term& src, hol_array_term& dst) {
	dst.length = src.length;
	dst.operands = src.operands;
}

inline void hol_quantifier::move(const hol_quantifier& src, hol_quantifier& dst) {
	dst.variable = src.variable;
	dst.operand = src.operand;
}

inline void hol_any::move(const hol_any& src, hol_any& dst) {
	dst.included = src.included;
	dst.excluded_trees = src.excluded_trees;
	dst.excluded_tree_count = src.excluded_tree_count;
}

inline void hol_any_array::move(const hol_any_array& src, hol_any_array& dst) {
	dst.oper = src.oper;
	dst.all = src.all;
	core::move(src.any, dst.any);
	core::move(src.left, dst.left);
	core::move(src.right, dst.right);
}

inline void hol_any_constant::move(const hol_any_constant& src, hol_any_constant& dst) {
	dst.length = src.length;
	dst.constants = src.constants;
}

inline void hol_any_quantifier::move(const hol_any_quantifier& src, hol_any_quantifier& dst) {
	dst.quantifier = src.quantifier;
	dst.operand = src.operand;
}

inline void hol_term::free_helper() {
	reference_count--;
	if (reference_count == 0) {
		switch (type) {
		case hol_term_type::NOT:
			core::free(unary); return;
		case hol_term_type::IF_THEN:
		case hol_term_type::EQUALS:
		case hol_term_type::UNARY_APPLICATION:
			core::free(binary); return;
		case hol_term_type::BINARY_APPLICATION:
			core::free(ternary); return;
		case hol_term_type::AND:
		case hol_term_type::OR:
		case hol_term_type::IFF:
			core::free(array); return;
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			core::free(quantifier); return;
		case hol_term_type::STRING:
			core::free(str); return;
		case hol_term_type::UINT_LIST:
			core::free(uint_list); return;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
			core::free(any); return;
		case hol_term_type::ANY_ARRAY:
			core::free(any_array); return;
		case hol_term_type::ANY_CONSTANT:
			core::free(any_constant); return;
		case hol_term_type::ANY_QUANTIFIER:
			core::free(any_quantifier); return;
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
		case hol_term_type::VARIABLE:
		case hol_term_type::CONSTANT:
		case hol_term_type::PARAMETER:
		case hol_term_type::INTEGER:
			return;
		}
		fprintf(stderr, "hol_term.free_helper ERROR: Unrecognized hol_term_type.\n");
	}
}

inline void hol_unary_term::free(hol_unary_term& term) {
	core::free(*term.operand);
	if (term.operand->reference_count == 0)
		core::free(term.operand);
}

inline void hol_binary_term::free(hol_binary_term& term) {
	core::free(*term.left);
	if (term.left->reference_count == 0)
		core::free(term.left);

	core::free(*term.right);
	if (term.right->reference_count == 0)
		core::free(term.right);
}

inline void hol_ternary_term::free(hol_ternary_term& term) {
	core::free(*term.first);
	if (term.first->reference_count == 0)
		core::free(term.first);

	core::free(*term.second);
	if (term.second->reference_count == 0)
		core::free(term.second);

	core::free(*term.third);
	if (term.third->reference_count == 0)
		core::free(term.third);
}

inline void hol_array_term::free(hol_array_term& term) {
	if (term.operands == nullptr) return;
	for (unsigned int i = 0; i < term.length; i++) {
		core::free(*term.operands[i]);
		if (term.operands[i]->reference_count == 0)
			core::free(term.operands[i]);
	}
	core::free(term.operands);
}

inline void hol_quantifier::free(hol_quantifier& term) {
	core::free(*term.operand);
	if (term.operand->reference_count == 0)
		core::free(term.operand);
}

inline void hol_any::free(hol_any& term) {
	if (term.included != nullptr) {
		core::free(*term.included);
		if (term.included->reference_count == 0)
			core::free(term.included);
	} if (term.excluded_trees != nullptr) {
		for (unsigned int i = 0; i < term.excluded_tree_count; i++) {
			core::free(*term.excluded_trees[i]);
			if (term.excluded_trees[i]->reference_count == 0)
				core::free(term.excluded_trees[i]);
		}
		core::free(term.excluded_trees);
	}
}

inline void hol_any_array::free(hol_any_array& term) {
	core::free(*term.all);
	if (term.all->reference_count == 0)
		core::free(term.all);
	core::free(term.any);
	core::free(term.left);
	core::free(term.right);
}

inline void hol_any_constant::free(hol_any_constant& term) {
	core::free(term.constants);
}

inline void hol_any_quantifier::free(hol_any_quantifier& term) {
	core::free(*term.operand);
	if (term.operand->reference_count == 0)
		core::free(term.operand);
}

inline bool is_atomic(
		const hol_term& term, unsigned int& predicate,
		hol_term const*& arg1, hol_term const*& arg2)
{
	if (term.type == hol_term_type::UNARY_APPLICATION) {
		if (term.binary.left->type != hol_term_type::CONSTANT) return false;
		predicate = term.binary.left->constant;
		arg1 = term.binary.right;
		arg2 = NULL;
		return true;
	} else if (term.type == hol_term_type::BINARY_APPLICATION) {
		if (term.ternary.first->type != hol_term_type::CONSTANT) return false;
		predicate = term.ternary.first->constant;
		arg1 = term.ternary.second;
		arg2 = term.ternary.third;
		return true;
	}
	return false;
}

inline bool is_atomic(const hol_term& term,
		hol_term const*& arg1, hol_term const*& arg2)
{
	unsigned int predicate;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(const hol_term& term) {
	unsigned int predicate;
	hol_term const* arg1; hol_term const* arg2;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(
		hol_term& term, unsigned int& predicate,
		hol_term*& arg1, hol_term*& arg2)
{
	if (term.type == hol_term_type::UNARY_APPLICATION) {
		if (term.binary.left->type != hol_term_type::CONSTANT) return false;
		predicate = term.binary.left->constant;
		arg1 = term.binary.right;
		arg2 = NULL;
		return true;
	} else if (term.type == hol_term_type::BINARY_APPLICATION) {
		if (term.ternary.first->type != hol_term_type::CONSTANT) return false;
		predicate = term.ternary.first->constant;
		arg1 = term.ternary.second;
		arg2 = term.ternary.third;
		return true;
	}
	return false;
}

inline bool is_atomic(hol_term& term,
		hol_term*& arg1, hol_term*& arg2)
{
	unsigned int predicate;
	return is_atomic(term, predicate, arg1, arg2);
}

inline bool is_atomic(hol_term& term) {
	unsigned int predicate;
	hol_term* arg1; hol_term* arg2;
	return is_atomic(term, predicate, arg1, arg2);
}


/* forward declarations for hol_term printing */

template<hol_term_syntax Syntax = hol_term_syntax::CLASSIC, typename Stream, typename... Printer>
bool print(const hol_term&, Stream&, Printer&&...);


template<hol_term_syntax Syntax> struct and_symbol;
template<hol_term_syntax Syntax> struct or_symbol;
template<hol_term_syntax Syntax> struct if_then_symbol;
template<hol_term_syntax Syntax> struct not_symbol;
template<hol_term_syntax Syntax> struct equals_symbol;
template<hol_term_syntax Syntax> struct for_all_symbol;
template<hol_term_syntax Syntax> struct exists_symbol;
template<hol_term_syntax Syntax> struct lambda_symbol;
template<hol_term_syntax Syntax> struct true_symbol;
template<hol_term_syntax Syntax> struct false_symbol;

template<> struct and_symbol<hol_term_syntax::TPTP> { static const char symbol[]; };
template<> struct or_symbol<hol_term_syntax::TPTP> { static const char symbol[]; };
template<> struct if_then_symbol<hol_term_syntax::TPTP> { static const char symbol[]; };
template<> struct not_symbol<hol_term_syntax::TPTP> { static const char symbol; };
template<> struct equals_symbol<hol_term_syntax::TPTP> { static const char symbol[]; };
template<> struct for_all_symbol<hol_term_syntax::TPTP> { static const char symbol; };
template<> struct exists_symbol<hol_term_syntax::TPTP> { static const char symbol; };
template<> struct lambda_symbol<hol_term_syntax::TPTP> { static const char symbol; };
template<> struct true_symbol<hol_term_syntax::TPTP> { static const char symbol; };
template<> struct false_symbol<hol_term_syntax::TPTP> { static const char symbol; };

template<> struct and_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct or_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct if_then_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct not_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct equals_symbol<hol_term_syntax::CLASSIC> { static const char symbol; };
template<> struct for_all_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct exists_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct lambda_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct true_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };
template<> struct false_symbol<hol_term_syntax::CLASSIC> { static const char symbol[]; };

const char and_symbol<hol_term_syntax::TPTP>::symbol[] = " & ";
const char or_symbol<hol_term_syntax::TPTP>::symbol[] = " | ";
const char if_then_symbol<hol_term_syntax::TPTP>::symbol[] = " => ";
const char not_symbol<hol_term_syntax::TPTP>::symbol = '~';
const char equals_symbol<hol_term_syntax::TPTP>::symbol[] = " = ";
const char for_all_symbol<hol_term_syntax::TPTP>::symbol = '!';
const char exists_symbol<hol_term_syntax::TPTP>::symbol = '?';
const char lambda_symbol<hol_term_syntax::TPTP>::symbol = '^';
const char true_symbol<hol_term_syntax::TPTP>::symbol = 'T';
const char false_symbol<hol_term_syntax::TPTP>::symbol = 'F';

const char and_symbol<hol_term_syntax::CLASSIC>::symbol[] = " ∧ ";
const char or_symbol<hol_term_syntax::CLASSIC>::symbol[] = " ∨ ";
const char if_then_symbol<hol_term_syntax::CLASSIC>::symbol[] = " → ";
const char not_symbol<hol_term_syntax::CLASSIC>::symbol[] = "¬";
const char equals_symbol<hol_term_syntax::CLASSIC>::symbol = '=';
const char for_all_symbol<hol_term_syntax::CLASSIC>::symbol[] = "∀";
const char exists_symbol<hol_term_syntax::CLASSIC>::symbol[] = "∃";
const char lambda_symbol<hol_term_syntax::CLASSIC>::symbol[] = "λ";
const char true_symbol<hol_term_syntax::CLASSIC>::symbol[] = "⊤";
const char false_symbol<hol_term_syntax::CLASSIC>::symbol[] = "⊥";

const char empty_string[] = "";
const char left_curly_brace[] = "{";
const char right_curly_brace[] = "}";
const char comma[] = ",";

template<hol_term_syntax Syntax, typename Stream, typename... Printer>
inline bool print_iff(const hol_array_term& term, Stream& out, Printer&&... printer) {
	if (term.length < 2) {
		fprintf(stderr, "print_iff ERROR: IFF term has fewer than two operands.\n");
		return false;
	}

	for (unsigned int i = 0; i < term.length - 1; i++) {
		if (!print('(', out) || !print<Syntax>(*term.operands[i], out, std::forward<Printer>(printer)...)
		 || !print(equals_symbol<Syntax>::symbol, out))
			return false;
	}
	if (!print<Syntax>(*term.operands[term.length - 1], out, std::forward<Printer>(printer)...)) return false;
	for (unsigned int i = 0; i < term.length - 1; i++)
		if (!print(')', out)) return false;
	return true;
}

template<hol_term_syntax Syntax, typename Stream>
inline bool print_quantifier(hol_quantifier_type quantifier, Stream& out) {
	switch (quantifier) {
	case hol_quantifier_type::FOR_ALL:
		return print(for_all_symbol<Syntax>::symbol, out);
	case hol_quantifier_type::EXISTS:
		return print(exists_symbol<Syntax>::symbol, out);
	case hol_quantifier_type::LAMBDA:
		return print(lambda_symbol<Syntax>::symbol, out);
	case hol_quantifier_type::FOR_ALL_OR_EXISTS:
		return print('{', out) && print(for_all_symbol<Syntax>::symbol, out) && print(',', out) && print(exists_symbol<Syntax>::symbol, out) && print('}', out);
	case hol_quantifier_type::FOR_ALL_OR_LAMBDA:
		return print('{', out) && print(for_all_symbol<Syntax>::symbol, out) && print(',', out) && print(lambda_symbol<Syntax>::symbol, out) && print('}', out);
	case hol_quantifier_type::EXISTS_OR_LAMBDA:
		return print('{', out) && print(exists_symbol<Syntax>::symbol, out) && print(',', out) && print(lambda_symbol<Syntax>::symbol, out) && print('}', out);
	case hol_quantifier_type::ANY:
		return print('*', out);
	}
	fprintf(stderr, "print_quantifier ERROR: Unrecognized hol_quantifier_type.\n");
	return false;
}

template<hol_term_syntax Syntax, typename Stream>
inline bool print_quantifier(hol_quantifier_type quantifier, unsigned int quantified_variable, Stream& out) {
	if (!print_quantifier<Syntax>(quantifier, out)) return false;
	switch (Syntax) {
	case hol_term_syntax::TPTP:
		return print('[', out) && print_variable<Syntax>(quantified_variable, out) && print("]:", out);
	case hol_term_syntax::CLASSIC:
		return print_variable<Syntax>(quantified_variable, out);
	}
	fprintf(stderr, "print_quantifier ERROR: Unrecognized hol_term_syntax.\n");
	return false;
}

template<hol_term_syntax Syntax, typename Stream, typename... Printer>
bool print(const hol_term& term, Stream& out, Printer&&... printer)
{
	bool first;
	switch (term.type) {
	case hol_term_type::VARIABLE:
		return print_variable<Syntax>(term.variable, out);

	case hol_term_type::CONSTANT:
		return print(term.constant, out, std::forward<Printer>(printer)...);

	case hol_term_type::PARAMETER:
		return print_parameter<Syntax>(term.parameter, out);

	case hol_term_type::INTEGER:
		return print(term.integer, out);

	case hol_term_type::STRING:
		return print('"', out) && print(term.str, out) && print('"', out);

	case hol_term_type::UINT_LIST:
		return print(term.uint_list.tokens, term.uint_list.length, out, std::forward<Printer>(printer)...);

	case hol_term_type::TRUE:
		return print(true_symbol<Syntax>::symbol, out);

	case hol_term_type::FALSE:
		return print(false_symbol<Syntax>::symbol, out);

	case hol_term_type::NOT:
		return print(not_symbol<Syntax>::symbol, out) && print(*term.unary.operand, out, std::forward<Printer>(printer)...);

	case hol_term_type::AND:
		return print<hol_term*, empty_string, empty_string, and_symbol<Syntax>::symbol>(term.array.operands, term.array.length, out, pointer_scribe(), std::forward<Printer>(printer)...);

	case hol_term_type::OR:
		return print<hol_term*, empty_string, empty_string, or_symbol<Syntax>::symbol>(term.array.operands, term.array.length, out, pointer_scribe(), std::forward<Printer>(printer)...);

	case hol_term_type::IFF:
		return print_iff<Syntax>(term.array, out, std::forward<Printer>(printer)...);

	case hol_term_type::IF_THEN:
		return print('(', out) && print(*term.binary.left, out, std::forward<Printer>(printer)...)
			&& print(if_then_symbol<Syntax>::symbol, out) && print(*term.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case hol_term_type::EQUALS:
		return print(*term.binary.left, out, std::forward<Printer>(printer)...)
			&& print(equals_symbol<Syntax>::symbol, out)
			&& print(*term.binary.right, out, std::forward<Printer>(printer)...);

	case hol_term_type::UNARY_APPLICATION:
		return print(*term.binary.left, out, std::forward<Printer>(printer)...) && print('(', out)
			&& print(*term.binary.right, out, std::forward<Printer>(printer)...) && print(')', out);

	case hol_term_type::BINARY_APPLICATION:
		return print(*term.ternary.first, out, std::forward<Printer>(printer)...) && print('(', out)
			&& print(*term.ternary.second, out, std::forward<Printer>(printer)...) && print(',', out)
			&& print(*term.ternary.third, out, std::forward<Printer>(printer)...) && print(')', out);

	case hol_term_type::FOR_ALL:
		if (!print_quantifier<Syntax>(hol_quantifier_type::FOR_ALL, term.quantifier.variable, out)) return false;
		if (term.quantifier.operand->type == hol_term_type::AND
		 || term.quantifier.operand->type == hol_term_type::EQUALS
		 || term.quantifier.operand->type == hol_term_type::IFF
		 || term.quantifier.operand->type == hol_term_type::OR)
		{
			return print('(', out) && print(*term.quantifier.operand, out, std::forward<Printer>(printer)...) && print(')', out);
		} else {
			return print(*term.quantifier.operand, out, std::forward<Printer>(printer)...);
		}

	case hol_term_type::EXISTS:
		if (!print_quantifier<Syntax>(hol_quantifier_type::EXISTS, term.quantifier.variable, out)) return false;
		if (term.quantifier.operand->type == hol_term_type::AND
		 || term.quantifier.operand->type == hol_term_type::EQUALS
		 || term.quantifier.operand->type == hol_term_type::IFF
		 || term.quantifier.operand->type == hol_term_type::OR)
		{
			return print('(', out) && print(*term.quantifier.operand, out, std::forward<Printer>(printer)...) && print(')', out);
		} else {
			return print(*term.quantifier.operand, out, std::forward<Printer>(printer)...);
		}

	case hol_term_type::LAMBDA:
		if (!print_quantifier<Syntax>(hol_quantifier_type::LAMBDA, term.quantifier.variable, out)) return false;
		if (term.quantifier.operand->type == hol_term_type::AND
		 || term.quantifier.operand->type == hol_term_type::EQUALS
		 || term.quantifier.operand->type == hol_term_type::IFF
		 || term.quantifier.operand->type == hol_term_type::OR)
		{
			return print('(', out) && print(*term.quantifier.operand, out, std::forward<Printer>(printer)...) && print(')', out);
		} else {
			return print(*term.quantifier.operand, out, std::forward<Printer>(printer)...);
		}

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		if (term.any.included == nullptr) {
			if (!print('*', out)) return false;
		} else {
			const char* any_string = (term.type == hol_term_type::ANY ? "C(" : "R(");
			if (!print(any_string, out) || !print(*term.any.included, out, std::forward<Printer>(printer)...) || !print(')', out))
				return false;
		} if (term.any.excluded_tree_count != 0) {
			if (!print("∖", out)) return false;
		} if (term.any.excluded_tree_count > 1) {
			if (!print('(', out)) return false;
		}
		first = true;
		for (unsigned int i = 0; i < term.any.excluded_tree_count; i++) {
			if (!first && !print(" ⋃ ", out)) return false;
			if (first) first = false;
			if (!print(*term.any.excluded_trees[i], out, std::forward<Printer>(printer)...))
				return false;
		}
		if (term.any.excluded_tree_count > 1) {
			if (!print(')', out)) return false;
		}
		return true;

	case hol_term_type::ANY_ARRAY:
		if (term.any_array.oper == hol_term_type::AND) {
			if (!print("AND array", out)) return false;
		} else if (term.any_array.oper == hol_term_type::OR) {
			if (!print("OR array", out)) return false;
		} else if (term.any_array.oper == hol_term_type::IFF) {
			if (!print("IFF array", out)) return false;
		} else if (term.any_array.oper == hol_term_type::ANY_ARRAY) {
			if (!print("* array", out)) return false;
		} else {
			fprintf(stderr, "print ERROR: Unexpected operator type in ANY_ARRAY.\n");
			return false;
		}
		if (!print(", all: ", out) || !print(*term.any_array.all, out, std::forward<Printer>(printer)...)) return false;
		if (term.any_array.any.length > 0) {
			if (!print(", any: ", out) || !print<hol_term*>(term.any_array.any.operands, term.any_array.any.length, out, pointer_scribe(), std::forward<Printer>(printer)...)) return false;
		} if (term.any_array.left.length > 0) {
			if (!print(", left: ", out) || !print<hol_term*>(term.any_array.left.operands, term.any_array.left.length, out, pointer_scribe(), std::forward<Printer>(printer)...)) return false;
		} if (term.any_array.right.length > 0) {
			if (!print(", right: ", out) || !print<hol_term*>(term.any_array.right.operands, term.any_array.right.length, out, pointer_scribe(), std::forward<Printer>(printer)...)) return false;
		}
		return true;

	case hol_term_type::ANY_CONSTANT:
		return print("(* in ", out)
			&& print<unsigned int, left_curly_brace, right_curly_brace, comma>(term.any_constant.constants, term.any_constant.length, out, std::forward<Printer>(printer)...)
			&& print(')', out);

	case hol_term_type::ANY_QUANTIFIER:
		if (!print_quantifier<Syntax>(term.any_quantifier.quantifier, ANY_VARIABLE, out)) return false;
		if (term.any_quantifier.operand->type == hol_term_type::AND
		 || term.any_quantifier.operand->type == hol_term_type::EQUALS
		 || term.any_quantifier.operand->type == hol_term_type::IF_THEN
		 || term.any_quantifier.operand->type == hol_term_type::IFF
		 || term.any_quantifier.operand->type == hol_term_type::OR)
		{
			return print('(', out) && print(*term.any_quantifier.operand, out, std::forward<Printer>(printer)...) && print(')', out);
		} else {
			return print(*term.any_quantifier.operand, out, std::forward<Printer>(printer)...);
		}
	}

	fprintf(stderr, "print ERROR: Unrecognized hol_term_type.\n");
	return false;
}

inline bool new_hol_term(hol_term*& new_term) {
	new_term = (hol_term*) malloc(sizeof(hol_term));
	if (new_term == NULL) {
		fprintf(stderr, "new_hol_term ERROR: Out of memory.\n");
		return false;
	}
	return true;
}

template<hol_term_type Type>
constexpr bool visit(const hol_term& term) { return true; }

template<hol_term_type Type, typename... Visitor>
constexpr bool end_visit(const hol_term& term, Visitor&&... visitor) { return true; }

template<typename... Visitor>
inline bool visit(hol_array_term& array, Visitor&&... visitor) {
	for (unsigned int i = 0; i < array.length; i++)
		if (!visit(*array.operands[i], std::forward<Visitor>(visitor)...)) return false;
	return true;
}

template<typename... Visitor>
inline bool visit(const hol_array_term& array, Visitor&&... visitor) {
	for (unsigned int i = 0; i < array.length; i++)
		if (!visit(*array.operands[i], std::forward<Visitor>(visitor)...)) return false;
	return true;
}

template<typename Term, typename... Visitor,
	typename std::enable_if<std::is_same<typename std::remove_cv<typename std::remove_reference<Term>::type>::type, hol_term>::value>::type* = nullptr>
bool visit(Term&& term, Visitor&&... visitor)
{
	switch (term.type) {
	case hol_term_type::CONSTANT:
		return visit<hol_term_type::CONSTANT>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::VARIABLE:
		return visit<hol_term_type::VARIABLE>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::PARAMETER:
		return visit<hol_term_type::PARAMETER>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::INTEGER:
		return visit<hol_term_type::INTEGER>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::STRING:
		return visit<hol_term_type::STRING>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::UINT_LIST:
		return visit<hol_term_type::UINT_LIST>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::NOT:
		return visit<hol_term_type::NOT>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.unary.operand, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::NOT>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::AND:
		return visit<hol_term_type::AND>(term, std::forward<Visitor>(visitor)...)
			&& visit(term.array, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::AND>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::OR:
		return visit<hol_term_type::OR>(term, std::forward<Visitor>(visitor)...)
			&& visit(term.array, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::OR>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::IFF:
		return visit<hol_term_type::IFF>(term, std::forward<Visitor>(visitor)...)
			&& visit(term.array, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::IFF>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::IF_THEN:
		return visit<hol_term_type::IF_THEN>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.left, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.right, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::IF_THEN>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::EQUALS:
		return visit<hol_term_type::EQUALS>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.left, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.right, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::EQUALS>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::UNARY_APPLICATION:
		return visit<hol_term_type::UNARY_APPLICATION>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.left, std::forward<Visitor>(visitor)...)
			&& visit(*term.binary.right, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::UNARY_APPLICATION>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::BINARY_APPLICATION:
		return visit<hol_term_type::BINARY_APPLICATION>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.ternary.first, std::forward<Visitor>(visitor)...)
			&& visit(*term.ternary.second, std::forward<Visitor>(visitor)...)
			&& visit(*term.ternary.third, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::BINARY_APPLICATION>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::FOR_ALL:
		return visit<hol_term_type::FOR_ALL>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.quantifier.operand, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::FOR_ALL>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::EXISTS:
		return visit<hol_term_type::EXISTS>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.quantifier.operand, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::EXISTS>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::LAMBDA:
		return visit<hol_term_type::LAMBDA>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.quantifier.operand, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::LAMBDA>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::TRUE:
		return visit<hol_term_type::TRUE>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::FALSE:
		return visit<hol_term_type::FALSE>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::ANY:
		/* NOTE: by default, `visit` does not traverse excluded subtrees in `hol_any` */
		return visit<hol_term_type::ANY>(term, std::forward<Visitor>(visitor)...)
			&& (term.any.included == nullptr || visit(*term.any.included, std::forward<Visitor>(visitor)...))
			&& end_visit<hol_term_type::ANY>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::ANY_RIGHT:
		return visit<hol_term_type::ANY_RIGHT>(term, std::forward<Visitor>(visitor)...)
			&& (term.any.included == nullptr || visit(*term.any.included, std::forward<Visitor>(visitor)...))
			&& end_visit<hol_term_type::ANY_RIGHT>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::ANY_ARRAY:
		return visit<hol_term_type::ANY_ARRAY>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.any_array.all, std::forward<Visitor>(visitor)...)
			&& visit(term.any_array.any, std::forward<Visitor>(visitor)...)
			&& visit(term.any_array.left, std::forward<Visitor>(visitor)...)
			&& visit(term.any_array.right, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::ANY_ARRAY>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::ANY_CONSTANT:
		return visit<hol_term_type::ANY_CONSTANT>(term, std::forward<Visitor>(visitor)...);
	case hol_term_type::ANY_QUANTIFIER:
		return visit<hol_term_type::ANY_QUANTIFIER>(term, std::forward<Visitor>(visitor)...)
			&& visit(*term.any_quantifier.operand, std::forward<Visitor>(visitor)...)
			&& end_visit<hol_term_type::ANY_QUANTIFIER>(term, std::forward<Visitor>(visitor)...);
	}
	fprintf(stderr, "visit ERROR: Unrecognized hol_term_type.\n");
	return false;
}

struct parameter_comparator {
	unsigned int parameter;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, const parameter_comparator& visitor) {
	if (Type == hol_term_type::PARAMETER)
		return visitor.parameter != term.parameter;
	else return true;
}

inline bool contains_parameter(const hol_term& src, unsigned int parameter) {
	parameter_comparator visitor = {parameter};
	return !visit(src, visitor);
}

struct constant_comparator {
	unsigned int constant;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, const constant_comparator& visitor) {
	if (Type == hol_term_type::CONSTANT)
		return visitor.constant != term.constant;
	else return true;
}

inline bool contains_constant(const hol_term& src, unsigned int constant) {
	constant_comparator visitor = {constant};
	return !visit(src, visitor);
}

struct variable_collector {
	array<unsigned int>& variables;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, const variable_collector& visitor) {
	if (Type == hol_term_type::VARIABLE) {
		if (!visitor.variables.contains(term.variable))
			return visitor.variables.add(term.variable);
	}
	return true;
}

inline bool get_variables(const hol_term& src, array<unsigned int>& variables) {
	variable_collector visitor = {variables};
	return !visit(src, visitor);
}

struct free_variable_collector {
	array<unsigned int>& variables;
	array<unsigned int> bound_variables;

	free_variable_collector(array<unsigned int>& variables) : variables(variables), bound_variables(8) { }
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, free_variable_collector& visitor) {
	if (Type == hol_term_type::VARIABLE) {
		if (!visitor.bound_variables.contains(term.variable) && !visitor.variables.contains(term.variable))
			return visitor.variables.add(term.variable);
	} else if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		return visitor.bound_variables.add(term.quantifier.variable);
	}
	return true;
}

template<hol_term_type Type>
inline bool end_visit(const hol_term& term, free_variable_collector& visitor) {
	if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA)
		visitor.bound_variables.remove(visitor.bound_variables.index_of(term.quantifier.variable));
	return true;
}

inline bool get_free_variables(const hol_term& src, array<unsigned int>& variables) {
	free_variable_collector visitor(variables);
	return !visit(src, visitor);
}

struct max_bound_variable_collector {
	bool has_variable;
	unsigned int variable;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, max_bound_variable_collector& visitor) {
	if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		visitor.variable = max(visitor.variable, term.quantifier.variable);
		visitor.has_variable = true;
	}
	return true;
}

inline bool max_bound_variable(const hol_term& src, unsigned int& max_variable) {
	max_bound_variable_collector visitor = {false, 0};
	visit(src, visitor);
	if (!visitor.has_variable)
		return false;
	max_variable = max(max_variable, visitor.variable);
	return true;
}

struct parameter_collector {
	array<unsigned int>& parameters;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, const parameter_collector& visitor) {
	if (Type == hol_term_type::PARAMETER) {
		if (!visitor.parameters.contains(term.parameter))
			return visitor.parameters.add(term.parameter);
	}
	return true;
}

inline bool get_parameters(const hol_term& src, array<unsigned int>& parameters) {
	parameter_collector visitor = {parameters};
	return !visit(src, visitor);
}

struct index_computer {
	const hol_term& term;
	array<unsigned int>& indices;
	unsigned int index;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, index_computer& visitor) {
	if (term == visitor.term) {
		if (!visitor.indices.add(visitor.index)) return false;
	}
	visitor.index++;
	return true;
}

inline bool compute_indices(const hol_term& src,
		const hol_term& term, array<unsigned int>& indices)
{
	index_computer visitor = {term, indices, 0};
	return visit(src, visitor);
}

struct constant_collector {
	array_multiset<unsigned int>& constants;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, constant_collector& visitor) {
	if (Type == hol_term_type::CONSTANT)
		return visitor.constants.add(term.constant);
	return true;
}

inline bool get_constants(const hol_term& src,
		array_multiset<unsigned int>& constants)
{
	constant_collector visitor = {constants};
	return visit(src, visitor);
}

struct unambiguity_visitor { };

template<hol_term_type Type>
inline bool visit(const hol_term& term, unambiguity_visitor& visitor) {
	if (Type == hol_term_type::ANY || Type == hol_term_type::ANY_RIGHT || Type == hol_term_type::ANY_ARRAY
	 || Type == hol_term_type::ANY_CONSTANT || Type == hol_term_type::ANY_QUANTIFIER)
		return false;
	return true;
}

inline bool is_ambiguous(const hol_term& src) {
	unambiguity_visitor visitor;
	return !visit(src, visitor);
}

inline bool clone_constant(unsigned int src_constant, unsigned int& dst_constant) {
	dst_constant = src_constant;
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

inline bool clone_integer(int src_integer, int& dst_integer) {
	dst_integer = src_integer;
	return true;
}

inline bool clone_string(const string& src_string, string& dst_string) {
	dst_string = src_string;
	return true;
}

inline bool clone_uint_list(const sequence& src_list, sequence& dst_list) {
	dst_list = src_list;
	return true;
}

template<typename... Cloner>
inline bool clone(const hol_array_term& src, hol_array_term& dst, Cloner&&... cloner)
{
	dst.length = src.length;
	if (src.length == 0) {
		dst.operands = nullptr;
		return true;
	}
	dst.operands = (hol_term**) malloc(sizeof(hol_term*) * src.length);
	if (dst.operands == NULL) return false;
	for (unsigned int i = 0; i < src.length; i++) {
		if (!new_hol_term(dst.operands[i])
		 || !clone(*src.operands[i], *dst.operands[i], std::forward<Cloner>(cloner)...))
		{
			for (unsigned int j = 0; j < i; j++) {
				free(*dst.operands[j]); free(dst.operands[j]);
			}
			if (dst.operands[i] != NULL) free(dst.operands[i]);
			free(dst.operands); return false;
		}
	}
	return true;
}

template<typename... Cloner>
bool clone(const hol_term& src, hol_term& dst, Cloner&&... cloner)
{
	dst.type = src.type;
	dst.reference_count = 1;
	switch (src.type) {
	case hol_term_type::CONSTANT:
		return clone_constant(src.constant, dst.constant, std::forward<Cloner>(cloner)...);
	case hol_term_type::VARIABLE:
		return clone_variable(src.variable, dst.variable, std::forward<Cloner>(cloner)...);
	case hol_term_type::PARAMETER:
		return clone_parameter(src.parameter, dst.parameter, std::forward<Cloner>(cloner)...);
	case hol_term_type::INTEGER:
		return clone_integer(src.integer, dst.integer, std::forward<Cloner>(cloner)...);
	case hol_term_type::STRING:
		return clone_string(src.str, dst.str, std::forward<Cloner>(cloner)...);
	case hol_term_type::UINT_LIST:
		return clone_uint_list(src.uint_list, dst.uint_list, std::forward<Cloner>(cloner)...);
	case hol_term_type::NOT:
		if (!new_hol_term(dst.unary.operand)) return false;
		if (!clone(*src.unary.operand, *dst.unary.operand, std::forward<Cloner>(cloner)...)) {
			free(dst.unary.operand);
			return false;
		}
		return true;
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		if (!new_hol_term(dst.binary.left)) {
			return false;
		} else if (!new_hol_term(dst.binary.right)) {
			free(dst.binary.left); return false;
		} else if (!clone(*src.binary.left, *dst.binary.left, std::forward<Cloner>(cloner)...)) {
			free(dst.binary.left); free(dst.binary.right);
			return false;
		} else if (!clone(*src.binary.right, *dst.binary.right, std::forward<Cloner>(cloner)...)) {
			free(*dst.binary.left); free(dst.binary.left);
			free(dst.binary.right); return false;
		}
		return true;
	case hol_term_type::BINARY_APPLICATION:
		if (!new_hol_term(dst.ternary.first)) {
			return false;
		} else if (!new_hol_term(dst.ternary.second)) {
			free(dst.ternary.first); return false;
		} else if (!new_hol_term(dst.ternary.third)) {
			free(dst.ternary.first); free(dst.ternary.second); return false;
		} else if (!clone(*src.ternary.first, *dst.ternary.first, std::forward<Cloner>(cloner)...)) {
			free(dst.ternary.first); free(dst.ternary.second);
			free(dst.ternary.third); return false;
		} else if (!clone(*src.ternary.second, *dst.ternary.second, std::forward<Cloner>(cloner)...)) {
			free(*dst.ternary.first); free(dst.ternary.first);
			free(dst.ternary.second); free(dst.ternary.third); return false;
		} else if (!clone(*src.ternary.third, *dst.ternary.third, std::forward<Cloner>(cloner)...)) {
			free(*dst.ternary.first); free(dst.ternary.first);
			free(*dst.ternary.second); free(dst.ternary.second);
			free(dst.ternary.third); return false;
		}
		return true;
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return clone(src.array, dst.array, std::forward<Cloner>(cloner)...);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!clone_variable(src.quantifier.variable, dst.quantifier.variable, std::forward<Cloner>(cloner)...)
		 || !new_hol_term(dst.quantifier.operand))
			return false;
		if (!clone(*src.quantifier.operand, *dst.quantifier.operand, std::forward<Cloner>(cloner)...)) {
			free(dst.quantifier.operand); return false;
		}
		return true;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		dst.any.excluded_tree_count = src.any.excluded_tree_count;
		if (src.any.included != nullptr) {
			if (!new_hol_term(dst.any.included)
			 || !clone(*src.any.included, *dst.any.included, std::forward<Cloner>(cloner)...))
			{
				if (dst.any.included != nullptr) free(dst.any.included);
				return false;
			}
		} else {
			dst.any.included = nullptr;
		}
		if (src.any.excluded_tree_count != 0) {
			dst.any.excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * src.any.excluded_tree_count);
			if (dst.any.excluded_trees == nullptr) {
				if (dst.any.included != nullptr) {
					core::free(*dst.any.included);
					core::free(dst.any.included);
				}
				return false;
			}
			for (unsigned int i = 0; i < dst.any.excluded_tree_count; i++) {
				if (!new_hol_term(dst.any.excluded_trees[i])
				 || !clone(*src.any.excluded_trees[i], *dst.any.excluded_trees[i], std::forward<Cloner>(cloner)...))
				{
					if (dst.any.excluded_trees[i] != nullptr) core::free(dst.any.excluded_trees[i]);
					for (unsigned int j = 0; j < i; j++) {
						core::free(*dst.any.excluded_trees[j]); core::free(dst.any.excluded_trees[j]);
					}
					free(dst.any.excluded_trees);
					if (dst.any.included != nullptr) {
						core::free(*dst.any.included);
						core::free(dst.any.included);
					}
					return false;
				}
			}
		} else {
			dst.any.excluded_trees = nullptr;
		}
		return true;
	case hol_term_type::ANY_ARRAY:
		dst.any_array.oper = src.any_array.oper;
		if (!new_hol_term(dst.any_array.all)
		 || !clone(*src.any_array.all, *dst.any_array.all, std::forward<Cloner>(cloner)...))
		{
			if (dst.any_array.all != nullptr) core::free(dst.any_array.all);
			return false;
		} if (!clone(src.any_array.left, dst.any_array.left, std::forward<Cloner>(cloner)...)) {
			core::free(*dst.any_array.all); if (dst.any_array.all->reference_count == 0) core::free(dst.any_array.all);
			return false;
		} if (!clone(src.any_array.right, dst.any_array.right, std::forward<Cloner>(cloner)...)) {
			core::free(*dst.any_array.all); if (dst.any_array.all->reference_count == 0) core::free(dst.any_array.all);
			core::free(dst.any_array.left);
			return false;
		} if (!clone(src.any_array.any, dst.any_array.any, std::forward<Cloner>(cloner)...)) {
			core::free(*dst.any_array.all); if (dst.any_array.all->reference_count == 0) core::free(dst.any_array.all);
			core::free(dst.any_array.left); core::free(dst.any_array.right);
			return false;
		}
		return true;
	case hol_term_type::ANY_CONSTANT:
		dst.any_constant.length = src.any_constant.length;
		dst.any_constant.constants = (unsigned int*) malloc(sizeof(unsigned int) * src.any_constant.length);
		if (dst.any_constant.constants == nullptr) return false;
		for (unsigned int i = 0; i < src.any_constant.length; i++) {
			if (!clone_constant(src.any_constant.constants[i], dst.any_constant.constants[i], std::forward<Cloner>(cloner)...)) {
				free(dst.any_constant.constants);
				return false;
			}
		}
		return true;
	case hol_term_type::ANY_QUANTIFIER:
		dst.any_quantifier.quantifier = src.any_quantifier.quantifier;
		if (!new_hol_term(dst.any_quantifier.operand)
		 || !clone(*src.any_quantifier.operand, *dst.any_quantifier.operand, std::forward<Cloner>(cloner)...))
		{
			if (dst.any_quantifier.operand != nullptr) core::free(dst.any_quantifier.operand);
			return false;
		}
		return true;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	}
	fprintf(stderr, "clone ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<typename... Cloner>
inline bool clone(const hol_term* src, hol_term* dst, Cloner&&... cloner)
{
	if (!new_hol_term(dst)) return false;
	return clone(*src, *dst, std::forward<Cloner>(cloner)...);
}

struct constant_relabeler {
	const array_map<unsigned int, unsigned int>& map;
};

inline bool clone_constant(unsigned int src_constant, unsigned int& dst_constant, constant_relabeler& relabeler) {
	bool contains;
	const unsigned int& dst = relabeler.map.get(src_constant, contains);
	if (contains) {
		dst_constant = dst;
		return true;
	} else {
		dst_constant = src_constant;
		return true;
	}
}

inline bool clone_predicate(unsigned int src_predicate, unsigned int& dst_predicate, constant_relabeler& relabeler) {
	bool contains;
	const unsigned int& dst = relabeler.map.get(src_predicate, contains);
	if (contains) {
		dst_predicate = dst;
		return true;
	} else {
		dst_predicate = src_predicate;
		return true;
	}
}

inline bool clone_variable(unsigned int src_variable, unsigned int& dst_variable, constant_relabeler& relabeler) {
	return clone_variable(src_variable, dst_variable);
}

inline bool clone_parameter(unsigned int src_parameter, unsigned int& dst_parameter, constant_relabeler& relabeler) {
	return clone_parameter(src_parameter, dst_parameter);
}

inline bool clone_integer(int src_integer, int& dst_integer, constant_relabeler& relabeler) {
	return clone_integer(src_integer, dst_integer);
}

inline bool clone_string(const string& src_string, string& dst_string, constant_relabeler& relabeler) {
	return clone_string(src_string, dst_string);
}

inline bool clone_uint_list(const sequence& src_list, sequence& dst_list, constant_relabeler& relabeler) {
	return clone_uint_list(src_list, dst_list);
}

template<typename Formula>
inline Formula* relabel_constants(const Formula* src,
		const array_map<unsigned int, unsigned int>& constant_map)
{
	Formula* dst = (Formula*) malloc(sizeof(Formula));
	if (dst == NULL) {
		fprintf(stderr, "relabel_constants ERROR: Out of memory.\n");
		return NULL;
	}

	constant_relabeler relabeler = {constant_map};
	if (!clone(*src, *dst, relabeler)) {
		free(dst); return NULL;
	}
	return dst;
}

template<typename... Function>
hol_term** default_apply_array(const hol_array_term& src, bool& changed, Function&&... function)
{
	hol_term** new_terms = (hol_term**) malloc(sizeof(hol_term*) * src.length);
	if (new_terms == nullptr) return nullptr;
	for (unsigned int i = 0; i < src.length; i++) {
		new_terms[i] = apply(src.operands[i], std::forward<Function>(function)...);
		if (new_terms[i] == nullptr) {
			for (unsigned int j = 0; j < i; j++) {
				if (new_terms[j] != src.operands[j]) {
					free(*new_terms[j]); if (new_terms[j]->reference_count == 0) free(new_terms[j]);
				}
			}
			free(new_terms); return nullptr;
		} else if (new_terms[i] != src.operands[i])
			changed = true;
	}
	return new_terms;
}

template<typename... Function>
hol_term* default_apply_any_array(hol_term* src, Function&&... function)
{
	bool changed = false;
	hol_term* all = apply(src->any_array.all, std::forward<Function>(function)...);
	if (all != src->any_array.all)
		changed = true;

	hol_term** any = default_apply_array(src->any_array.any, changed, std::forward<Function>(function)...);
	if (any == nullptr) {
		if (all != src->any_array.all) { free(*all); if (all->reference_count == 0) free(all); }
		return nullptr;
	}

	hol_term** left = default_apply_array(src->any_array.left, changed, std::forward<Function>(function)...);
	if (left == nullptr) {
		if (all != src->any_array.all) { free(*all); if (all->reference_count == 0) free(all); }
		for (unsigned int j = 0; j < src->any_array.any.length; j++) {
			if (any[j] != src->any_array.any.operands[j]) { free(*any[j]); if (any[j]->reference_count == 0) free(any[j]); }
		}
		free(any);
		return nullptr;
	}

	hol_term** right = default_apply_array(src->any_array.right, changed, std::forward<Function>(function)...);
	if (right == nullptr) {
		if (all != src->any_array.all) { free(*all); if (all->reference_count == 0) free(all); }
		for (unsigned int j = 0; j < src->any_array.any.length; j++) {
			if (any[j] != src->any_array.any.operands[j]) { free(*any[j]); if (any[j]->reference_count == 0) free(any[j]); }
		} for (unsigned int j = 0; j < src->any_array.left.length; j++) {
			if (left[j] != src->any_array.left.operands[j]) { free(*left[j]); if (left[j]->reference_count == 0) free(left[j]); }
		}
		free(any); free(left);
		return nullptr;
	}

	if (!changed) {
		free(any); free(left); free(right);
		return src;
	}

	hol_term* new_term;
	if (!new_hol_term(new_term)) {
		if (all != src->any_array.all) { free(*all); if (all->reference_count == 0) free(all); }
		for (unsigned int j = 0; j < src->any_array.any.length; j++) {
			if (any[j] != src->any_array.any.operands[j]) { free(*any[j]); if (any[j]->reference_count == 0) free(any[j]); }
		} for (unsigned int j = 0; j < src->any_array.left.length; j++) {
			if (left[j] != src->any_array.left.operands[j]) { free(*left[j]); if (left[j]->reference_count == 0) free(left[j]); }
		} for (unsigned int j = 0; j < src->any_array.right.length; j++) {
			if (right[j] != src->any_array.right.operands[j]) { free(*right[j]); if (right[j]->reference_count == 0) free(right[j]); }
		}
		free(any); free(left); free(right);
		return nullptr;
	}
	new_term->type = hol_term_type::ANY_ARRAY;
	new_term->reference_count = 1;
	new_term->any_array.oper = src->any_array.oper;
	new_term->any_array.all = all;
	new_term->any_array.left.operands = left;
	new_term->any_array.left.length = src->any_array.left.length;
	new_term->any_array.right.operands = right;
	new_term->any_array.right.length = src->any_array.right.length;
	new_term->any_array.any.operands = any;
	new_term->any_array.any.length = src->any_array.any.length;
	if (all == src->any_array.all) all->reference_count++;
	for (unsigned int j = 0; j < src->any_array.any.length; j++) {
		if (any[j] == src->any_array.any.operands[j]) any[j]->reference_count++;
	} for (unsigned int j = 0; j < src->any_array.left.length; j++) {
		if (left[j] == src->any_array.left.operands[j]) left[j]->reference_count++;
	} for (unsigned int j = 0; j < src->any_array.right.length; j++) {
		if (right[j] == src->any_array.right.operands[j]) right[j]->reference_count++;
	}
	return new_term;
}

/* NOTE: this function assumes `src.type == Type` */
template<hol_term_type Type, typename... Function>
hol_term* default_apply(hol_term* src, Function&&... function)
{
	hol_term* new_term;
	hol_term** new_terms;
	hol_term** other_terms;
	hol_term* first; hol_term* second; hol_term* third;
	bool changed;
	switch (Type) {
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
		return src;

	case hol_term_type::NOT:
		first = apply(src->unary.operand, std::forward<Function>(function)...);
		if (first == NULL) {
			return NULL;
		} else if (first == src->unary.operand) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				free(*first); if (first->reference_count == 0) free(first);
				return NULL;
			}
			new_term->unary.operand = first;
			new_term->type = hol_term_type::NOT;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		first = apply(src->binary.left, std::forward<Function>(function)...);
		if (first == NULL) return NULL;
		second = apply(src->binary.right, std::forward<Function>(function)...);
		if (second == NULL) {
			if (first != src->binary.left) {
				free(*first); if (first->reference_count == 0) free(first);
			}
			return NULL;
		} else if (first == src->binary.left && second == src->binary.right) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->binary.left) {
					free(*first); if (first->reference_count == 0) free(first);
				} if (second != src->binary.right) {
					free(*second); if (second->reference_count == 0) free(second);
				}
				return NULL;
			}
			new_term->binary.left = first;
			new_term->binary.right = second;
			if (first == src->binary.left) first->reference_count++;
			if (second == src->binary.right) second->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::BINARY_APPLICATION:
		first = apply(src->ternary.first, std::forward<Function>(function)...);
		if (first == NULL) return NULL;
		second = apply(src->ternary.second, std::forward<Function>(function)...);
		if (second == NULL) {
			if (first != src->ternary.first) {
				free(*first); if (first->reference_count == 0) free(first);
			}
			return NULL;
		}
		third = apply(src->ternary.third, std::forward<Function>(function)...);
		if (third == NULL) {
			if (first != src->ternary.first) {
				free(*first); if (first->reference_count == 0) free(first);
			} if (second != src->ternary.second) {
				free(*second); if (second->reference_count == 0) free(second);
			}
			return NULL;
		} else if (first == src->ternary.first && second == src->ternary.second && third == src->ternary.third) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->ternary.first) {
					free(*first); if (first->reference_count == 0) free(first);
				} if (second != src->ternary.second) {
					free(*second); if (second->reference_count == 0) free(second);
				} if (third != src->ternary.third) {
					free(*third); if (third->reference_count == 0) free(third);
				}
				return NULL;
			}
			new_term->ternary.first = first;
			new_term->ternary.second = second;
			new_term->ternary.third = third;
			if (first == src->ternary.first) first->reference_count++;
			if (second == src->ternary.second) second->reference_count++;
			if (third == src->ternary.third) third->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		changed = false;
		new_terms = default_apply_array(src->array, changed, std::forward<Function>(function)...);

		if (!changed) {
			free(new_terms);
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				for (unsigned int j = 0; j < src->array.length; j++) {
					if (new_terms[j] != src->array.operands[j]) {
						free(*new_terms[j]); if (new_terms[j]->reference_count == 0) free(new_terms[j]);
					}
				}
				free(new_terms); return NULL;
			}
			new_term->array.operands = new_terms;
			new_term->array.length = src->array.length;
			for (unsigned int i = 0; i < src->array.length; i++)
				if (new_term->array.operands[i] == src->array.operands[i]) src->array.operands[i]->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		first = apply(src->quantifier.operand, std::forward<Function>(function)...);
		if (first == NULL) {
			return NULL;
		} else if (first == src->quantifier.operand) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				free(*first); if (first->reference_count == 0) free(first);
				return NULL;
			}
			new_term->quantifier.variable = src->quantifier.variable;
			new_term->quantifier.operand = first;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		changed = false;
		if (src->any.included != 0) {
			first = apply(src->any.included, std::forward<Function>(function)...);
			if (first == nullptr) {
				return nullptr;
			} else if (first != src->any.included)
				changed = true;
		} else {
			first = nullptr;
		}
		if (src->any.excluded_tree_count != 0) {
			other_terms = (hol_term**) malloc(sizeof(hol_term*) * src->any.excluded_tree_count);
			if (other_terms == nullptr) {
				if (first != src->any.included) {
					free(*first); if (first->reference_count == 0) free(first);
				}
				return nullptr;
			}
			for (unsigned int i = 0; i < src->any.excluded_tree_count; i++) {
				other_terms[i] = apply(src->any.excluded_trees[i], std::forward<Function>(function)...);
				if (other_terms[i] == nullptr) {
					if (first != src->any.included) {
						free(*first); if (first->reference_count == 0) free(first);
					}
					for (unsigned int j = 0; j < i; j++) {
						if (other_terms[j] != src->any.excluded_trees[j]) {
							free(*other_terms[j]); if (other_terms[j]->reference_count == 0) free(other_terms[j]);
						}
					}
					free(other_terms);
					return nullptr;
				} else if (other_terms[i] != src->any.excluded_trees[i])
					changed = true;
			}
		} else {
			other_terms = nullptr;
		}

		if (!changed) {
			if (other_terms != nullptr) free(other_terms);
			return src;
		}

		if (!new_hol_term(new_term)) {
			if (first != src->any.included) {
				free(*first); if (first->reference_count == 0) free(first);
			} if (other_terms != nullptr) {
				for (unsigned int j = 0; j < src->any.excluded_tree_count; j++) {
					if (other_terms[j] != src->any.excluded_trees[j]) {
						free(*other_terms[j]); if (other_terms[j]->reference_count == 0) free(other_terms[j]);
					}
				}
				free(other_terms);
			}
			return nullptr;
		}
		new_term->any.included = first;
		new_term->any.excluded_trees = other_terms;
		new_term->any.excluded_tree_count = src->any.excluded_tree_count;
		if (first != nullptr && first == src->any.included)
			first->reference_count++;
		if (other_terms != nullptr)
			for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
				if (other_terms[i] == src->any.excluded_trees[i]) other_terms[i]->reference_count++;
		new_term->type = Type;
		new_term->reference_count = 1;
		return new_term;

	case hol_term_type::ANY_ARRAY:
		return default_apply_any_array(src, std::forward<Function>(function)...);

	case hol_term_type::ANY_QUANTIFIER:
		first = apply(src->any_quantifier.operand, std::forward<Function>(function)...);
		if (first == nullptr) {
			return nullptr;
		} else if (first == src->any_quantifier.operand) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				free(*first); if (first->reference_count == 0) free(first);
				return NULL;
			}
			new_term->any_quantifier.quantifier = src->any_quantifier.quantifier;
			new_term->any_quantifier.operand = first;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return src;
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type.\n");
	return NULL;
}

template<hol_term_type Type, typename... Function>
inline hol_term* apply(hol_term* src, Function&&... function)
{
	return default_apply<Type>(src, std::forward<Function>(function)...);
}

template<typename... Function>
inline hol_term* apply(hol_term* src, Function&&... function)
{
	switch (src->type) {
	case hol_term_type::CONSTANT:
		return apply<hol_term_type::CONSTANT>(src, std::forward<Function>(function)...);
	case hol_term_type::VARIABLE:
		return apply<hol_term_type::VARIABLE>(src, std::forward<Function>(function)...);
	case hol_term_type::PARAMETER:
		return apply<hol_term_type::PARAMETER>(src, std::forward<Function>(function)...);
	case hol_term_type::INTEGER:
		return apply<hol_term_type::INTEGER>(src, std::forward<Function>(function)...);
	case hol_term_type::STRING:
		return apply<hol_term_type::STRING>(src, std::forward<Function>(function)...);
	case hol_term_type::UINT_LIST:
		return apply<hol_term_type::UINT_LIST>(src, std::forward<Function>(function)...);
	case hol_term_type::NOT:
		return apply<hol_term_type::NOT>(src, std::forward<Function>(function)...);
	case hol_term_type::IF_THEN:
		return apply<hol_term_type::IF_THEN>(src, std::forward<Function>(function)...);
	case hol_term_type::EQUALS:
		return apply<hol_term_type::EQUALS>(src, std::forward<Function>(function)...);
	case hol_term_type::UNARY_APPLICATION:
		return apply<hol_term_type::UNARY_APPLICATION>(src, std::forward<Function>(function)...);
	case hol_term_type::BINARY_APPLICATION:
		return apply<hol_term_type::BINARY_APPLICATION>(src, std::forward<Function>(function)...);
	case hol_term_type::AND:
		return apply<hol_term_type::AND>(src, std::forward<Function>(function)...);
	case hol_term_type::OR:
		return apply<hol_term_type::OR>(src, std::forward<Function>(function)...);
	case hol_term_type::IFF:
		return apply<hol_term_type::IFF>(src, std::forward<Function>(function)...);
	case hol_term_type::FOR_ALL:
		return apply<hol_term_type::FOR_ALL>(src, std::forward<Function>(function)...);
	case hol_term_type::EXISTS:
		return apply<hol_term_type::EXISTS>(src, std::forward<Function>(function)...);
	case hol_term_type::LAMBDA:
		return apply<hol_term_type::LAMBDA>(src, std::forward<Function>(function)...);
	case hol_term_type::TRUE:
		return apply<hol_term_type::TRUE>(src, std::forward<Function>(function)...);
	case hol_term_type::FALSE:
		return apply<hol_term_type::FALSE>(src, std::forward<Function>(function)...);
	case hol_term_type::ANY:
		return apply<hol_term_type::ANY>(src, std::forward<Function>(function)...);
	case hol_term_type::ANY_RIGHT:
		return apply<hol_term_type::ANY_RIGHT>(src, std::forward<Function>(function)...);
	case hol_term_type::ANY_ARRAY:
		return apply<hol_term_type::ANY_ARRAY>(src, std::forward<Function>(function)...);
	case hol_term_type::ANY_CONSTANT:
		return apply<hol_term_type::ANY_CONSTANT>(src, std::forward<Function>(function)...);
	case hol_term_type::ANY_QUANTIFIER:
		return apply<hol_term_type::ANY_QUANTIFIER>(src, std::forward<Function>(function)...);
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type.\n");
	return NULL;
}

struct bound_variable_shifter {
	int shift;
	array<unsigned int> bound_variables;

	bound_variable_shifter(int shift) : shift(shift), bound_variables(8) { }
};

template<hol_term_type Type, typename std::enable_if<
	Type == hol_term_type::VARIABLE || Type == hol_term_type::FOR_ALL
 || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, bound_variable_shifter& shifter) {
	if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		if (!shifter.bound_variables.add(src->quantifier.variable)) return NULL;
		hol_term* operand = apply(src->quantifier.operand, shifter);
		if (operand == NULL) return NULL;
		else if (operand == src) operand->reference_count++;

		unsigned int new_variable = src->quantifier.variable + shifter.shift;
		hol_term* dst;
		if (Type == hol_term_type::LAMBDA) dst = hol_term::new_lambda(new_variable, operand);
		else if (Type == hol_term_type::FOR_ALL) dst = hol_term::new_for_all(new_variable, operand);
		else if (Type == hol_term_type::EXISTS) dst = hol_term::new_exists(new_variable, operand);
		if (dst == NULL) {
			free(*operand); if (operand->reference_count == 0) free(operand);
			return NULL;
		}
		return dst;
	}

	if (shifter.bound_variables.contains(src->variable)) {
		hol_term* dst;
		if (!new_hol_term(dst)) return NULL;
		dst->type = Type;
		dst->variable = src->variable + shifter.shift;
		dst->reference_count = 1;
		return dst;
	} else {
		return src;
	}
}

inline hol_term* shift_bound_variables(hol_term* src, int shift)
{
	bound_variable_shifter shifter(shift);
	hol_term* dst = apply(src, shifter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct variable_mapper {
	const array_map<unsigned int, unsigned int>& variable_map;
};

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, const variable_mapper& mapper) {
	if (Type == hol_term_type::VARIABLE) {
		unsigned int new_variable; bool contains;
		new_variable = mapper.variable_map.get(src->variable, contains);
		if (!contains) {
			return src;
		} else {
			return hol_term::new_variable(new_variable);
		}
	} else if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		unsigned int new_variable; bool contains;
		new_variable = mapper.variable_map.get(src->quantifier.variable, contains);
		if (!contains) {
			return default_apply<Type>(src, mapper);
		} else {
			hol_term* new_operand = apply(src->quantifier.operand, mapper);
			if (new_operand == nullptr)
				return nullptr;

			hol_term* new_quantifier;
			if (!new_hol_term(new_quantifier)) {
				free(*new_operand); if (new_operand->reference_count == 0) free(new_operand);
				return nullptr;
			}
			new_quantifier->type = Type;
			new_quantifier->reference_count = 1;
			new_quantifier->quantifier.variable = new_variable;
			new_quantifier->quantifier.operand = new_operand;
			return new_quantifier;
		}
	} else {
		return default_apply<Type>(src, mapper);
	}
}

inline hol_term* map_variables(hol_term* src,
		const array_map<unsigned int, unsigned int>& variable_map)
{
	const variable_mapper mapper = {variable_map};
	hol_term* dst = apply(src, mapper);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct variable_relabeler {
	array_map<unsigned int, unsigned int> variable_map;

	variable_relabeler() : variable_map(8) { }
};

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

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, variable_relabeler& relabeler) {
	if (Type == hol_term_type::VARIABLE) {
		unsigned int index = relabeler.variable_map.index_of(src->variable);
		unsigned int variable;
		if (index < relabeler.variable_map.size) {
			variable = relabeler.variable_map.values[index];
		} else if (!new_variable(src->variable, variable, relabeler.variable_map)) {
			return nullptr;
		}

		if (variable == src->variable) {
			return src;
		} else {
			return hol_term::new_variable(variable);
		}
	} else if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		unsigned int index = relabeler.variable_map.index_of(src->variable);
		unsigned int variable;
		if (index < relabeler.variable_map.size) {
			variable = relabeler.variable_map.values[index];
		} else if (!new_variable(src->quantifier.variable, variable, relabeler.variable_map)) {
			return nullptr;
		}

		hol_term* new_operand = apply(src->quantifier.operand, relabeler);
		relabeler.variable_map.size--;

		if (new_operand == nullptr)
			return nullptr;
		if (variable == src->quantifier.variable && new_operand == src->quantifier.operand)
			return src;

		hol_term* new_quantifier;
		if (Type == hol_term_type::FOR_ALL) {
			new_quantifier = hol_term::new_for_all(variable, new_operand);
		} else if (Type == hol_term_type::EXISTS) {
			new_quantifier = hol_term::new_exists(variable, new_operand);
		} else if (Type == hol_term_type::LAMBDA) {
			new_quantifier = hol_term::new_lambda(variable, new_operand);
		}
		if (new_operand == src->quantifier.operand)
			new_operand->reference_count++;
		return new_quantifier;
	} else {
		return default_apply<Type>(src, relabeler);
	}
}

inline hol_term* relabel_variables(hol_term* src)
{
	variable_relabeler relabeler;
	hol_term* dst = apply(src, relabeler);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

template<hol_term_type SrcTermType, int VariableShift>
struct static_term_substituter {
	const hol_term* src;
	hol_term* dst;
};

template<hol_term_type Type, hol_term_type SrcTermType, int VariableShift>
inline hol_term* apply(hol_term* src, const static_term_substituter<SrcTermType, VariableShift>& substituter) {
	if (Type == SrcTermType && *src == *substituter.src) {
		substituter.dst->reference_count++;
		return substituter.dst;
	} else if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		hol_term* operand = apply(src->quantifier.operand, substituter);
		if (operand == NULL) return NULL;
		else if (operand == src) operand->reference_count++;

		unsigned int new_variable = src->quantifier.variable + VariableShift;
		hol_term* dst;
		if (Type == hol_term_type::LAMBDA) dst = hol_term::new_lambda(new_variable, operand);
		else if (Type == hol_term_type::FOR_ALL) dst = hol_term::new_for_all(new_variable, operand);
		else if (Type == hol_term_type::EXISTS) dst = hol_term::new_exists(new_variable, operand);
		if (dst == NULL) {
			free(*operand); if (operand->reference_count == 0) free(operand);
			return NULL;
		}
		return dst;
		
	} else if (Type == hol_term_type::VARIABLE) {
		hol_term* dst;
		if (!new_hol_term(dst)) return NULL;
		dst->type = Type;
		dst->variable = src->variable + VariableShift;
		dst->reference_count = 1;
		return dst;
	} else {
		return default_apply<Type>(src, substituter);
	}
}

/* NOTE: this function assumes `src_term.type == SrcTermType` */
template<hol_term_type SrcTermType, int VariableShift = 0>
inline hol_term* substitute(hol_term* src,
		const hol_term* src_term, hol_term* dst_term)
{
	const static_term_substituter<SrcTermType, VariableShift> substituter = {src_term, dst_term};
	hol_term* dst = apply(src, substituter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct term_substituter {
	const hol_term* src;
	hol_term* dst;
};

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, const term_substituter& substituter) {
	if (*src == *substituter.src) {
		substituter.dst->reference_count++;
		return substituter.dst;
	} else {
		return default_apply<Type>(src, substituter);
	}
}

inline hol_term* substitute(hol_term* src,
		const hol_term* src_term, hol_term* dst_term)
{
	const term_substituter substituter = {src_term, dst_term};
	hol_term* dst = apply(src, substituter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

template<int VariableShift>
struct index_substituter {
	const hol_term* src;
	hol_term* dst;
	const unsigned int* term_indices;
	unsigned int term_index_count;
	unsigned int current_term_index;
};

template<hol_term_type Type, int VariableShift>
inline hol_term* apply(hol_term* src, index_substituter<VariableShift>& substituter)
{
	if (substituter.term_index_count > 0 && *substituter.term_indices == substituter.current_term_index) {
		if (substituter.src == NULL) {
			substituter.src = src;
		} else if (*substituter.src != *src) {
			/* this term is not identical to other substituted terms, which should not happen */
			return NULL;
		}
		hol_term* dst = substituter.dst;
		dst->reference_count++;
		substituter.term_indices++;
		substituter.term_index_count--;
		substituter.current_term_index++;
		return dst;
	} else {
		substituter.current_term_index++;
		return default_apply<Type>(src, substituter);
	}
}

/* NOTE: this function assumes `src.type == SrcType` */
template<int VariableShift = 0>
inline hol_term* substitute(
		hol_term* src, const unsigned int* term_indices,
		unsigned int term_index_count, hol_term* dst_term,
		const hol_term* expected_src_term = NULL)
{
	index_substituter<VariableShift> substituter = {expected_src_term, dst_term, term_indices, term_index_count, 0};
	hol_term* term = apply(src, substituter);
	if (term == src)
		term->reference_count++;
	return term;
}

struct beta_reducer {
	unsigned int max_variable;
	array_map<unsigned int, hol_term*> substitutions;

	beta_reducer(unsigned int max_variable, unsigned int initial_substitution_capacity) :
		max_variable(max_variable), substitutions(initial_substitution_capacity) { }
};

inline hol_term* beta_reduce(hol_term* left_src, hol_term* right_src, beta_reducer& reducer)
{
	hol_term* right = apply(right_src, reducer);
	if (right == NULL) return NULL;
	else if (right == right_src)
		right->reference_count++;

	hol_term* result;
	if (left_src->type == hol_term_type::LAMBDA) {
		reducer.max_variable = max(reducer.max_variable, left_src->quantifier.variable);
		if (!reducer.substitutions.put(left_src->quantifier.variable, right)) {
			free(*right); if (right->reference_count == 0) free(right);
			return NULL;
		}

		result = apply(left_src->quantifier.operand, reducer);
		free(*right); if (right->reference_count == 0) free(right);
		if (result == NULL) return NULL;
		else if (result == left_src->quantifier.operand)
			result->reference_count++;
#if !defined(NDEBUG)
		if (!reducer.substitutions.remove(left_src->quantifier.variable))
			fprintf(stderr, "beta_reduce WARNING: The variable in the lambda expression does not exist in `reducer.substitutions`.\n");
#else
		reducer.substitutions.remove(left_src->quantifier.variable);
#endif
	} else {
		hol_term* left = apply(left_src, reducer);
		if (left == NULL) {
			free(*right); if (right->reference_count == 0) free(right);
			return NULL;
		} else if (left == left_src)
			left->reference_count++;

		if (left->type == hol_term_type::LAMBDA) {
			reducer.max_variable = max(reducer.max_variable, left->quantifier.variable);
			beta_reducer new_reducer = beta_reducer(reducer.max_variable, 1);
			new_reducer.substitutions.keys[0] = left->quantifier.variable;
			new_reducer.substitutions.values[0] = right;
			new_reducer.substitutions.size++;

			result = apply(left->quantifier.operand, new_reducer);
			free(*right); if (right->reference_count == 0) free(right);
			if (result == NULL) {
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			} else if (result == left->quantifier.operand)
				result->reference_count++;
			free(*left); if (left->reference_count == 0) free(left);
		} else {
			result = hol_term::new_apply(left, right);
			if (result == NULL) {
				free(*right); if (right->reference_count == 0) free(right);
				free(*left); if (left->reference_count == 0) free(left);
				return NULL;
			}
		}
	}
	return result;
}

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, beta_reducer& reducer)
{
	if (Type == hol_term_type::UNARY_APPLICATION)
	{
		unsigned int old_max_variable = reducer.max_variable;
		hol_term* result = beta_reduce(src->binary.left, src->binary.right, reducer);
		reducer.max_variable = old_max_variable;
		return result;

	} else if (Type == hol_term_type::BINARY_APPLICATION) {
		unsigned int old_max_variable = reducer.max_variable;
		hol_term* reduced = beta_reduce(src->ternary.first, src->ternary.second, reducer);
		if (reduced == NULL) return NULL;
		hol_term* result = beta_reduce(reduced, src->ternary.third, reducer);
		free(*reduced); if (reduced->reference_count == 0) free(reduced);
		reducer.max_variable = old_max_variable;
		return result;
		
	} else if (Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::LAMBDA) {
		hol_term* result = default_apply<Type>(src, reducer);
		reducer.max_variable = max(reducer.max_variable, src->quantifier.variable);
		return result;

	} else if (Type == hol_term_type::VARIABLE) {
		bool contains;
		hol_term* dst = reducer.substitutions.get(src->variable, contains);
		if (contains) {
			hol_term* result;
			if (reducer.max_variable == 0) {
				result = dst;
				result->reference_count++;
			} else {
				result = shift_bound_variables(dst, reducer.max_variable);
			}
			return result;
		} else {
			return default_apply<Type>(src, reducer);
		}

	} else {
		return default_apply<Type>(src, reducer);
	}
}

inline hol_term* beta_reduce(hol_term* src)
{
	beta_reducer reducer(0, 8);
	hol_term* dst = apply(src, reducer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

bool unify(
		const hol_term& first, const hol_term& second,
		const hol_term* src_term, hol_term const*& dst_term)
{
	if (first.type != second.type) {
		return false;
	} else if (first == *src_term) {
		if (dst_term == NULL) {
			dst_term = &second;
		} else if (second != *dst_term) {
			return false;
		} else {
			return true;
		}
	}

	switch (first.type) {
	case hol_term_type::CONSTANT:
		return first.constant == second.constant;
	case hol_term_type::VARIABLE:
		return first.variable == second.variable;
	case hol_term_type::PARAMETER:
		return first.parameter == second.parameter;
	case hol_term_type::INTEGER:
		return first.integer == second.integer;
	case hol_term_type::STRING:
		return first.str == second.str;
	case hol_term_type::UINT_LIST:
		return first.uint_list == second.uint_list;
	case hol_term_type::NOT:
		return unify(*first.unary.operand, *second.unary.operand, src_term, dst_term);
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return unify(*first.binary.left, *second.binary.left, src_term, dst_term)
			&& unify(*first.binary.right, *second.binary.right, src_term, dst_term);
	case hol_term_type::BINARY_APPLICATION:
		return unify(*first.ternary.first, *second.ternary.first, src_term, dst_term)
			&& unify(*first.ternary.second, *second.ternary.second, src_term, dst_term)
			&& unify(*first.ternary.third, *second.ternary.third, src_term, dst_term);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		if (first.array.length != second.array.length) return false;
		for (unsigned int i = 0; i < first.array.length; i++)
			if (!unify(*first.array.operands[i], *second.array.operands[i], src_term, dst_term)) return false;
		return true;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (first.quantifier.variable != second.quantifier.variable) return false;
		return unify(*first.quantifier.operand, *second.quantifier.operand, src_term, dst_term);
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "unify ERROR: hol_term_types `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, `ANY_QUANTIFIER` unsupported.\n");
		return false;
	}
	fprintf(stderr, "unify ERROR: Unrecognized hol_term_type.\n");
	return false;
}

bool unifies_parameter(
		const hol_term& first, const hol_term& second,
		const hol_term* src_term, unsigned int& parameter)
{
	const hol_term* dst = NULL;
	if (!unify(first, second, src_term, dst)
	 || dst == NULL || dst->type != hol_term_type::PARAMETER)
		return false;
	parameter = dst->parameter;
	return true;
}


/**
 * Functions for easily constructing first-order logic expressions in code.
 */


hol_term* hol_term::new_variable(unsigned int variable) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::VARIABLE;
	term->variable = variable;
	return term;
}

hol_term* hol_term::new_constant(unsigned int constant) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::CONSTANT;
	term->constant = constant;
	return term;
}

hol_term* hol_term::new_parameter(unsigned int parameter) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::PARAMETER;
	term->parameter = parameter;
	return term;
}

hol_term* hol_term::new_int(int integer) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::INTEGER;
	term->integer = integer;
	return term;
}

hol_term* hol_term::new_string(const string& str) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::STRING;
	term->str = str;
	return term;
}

hol_term* hol_term::new_uint_list(const sequence& list) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::UINT_LIST;
	term->uint_list = list;
	return term;
}

hol_term* hol_term::new_atom(unsigned int predicate, hol_term* arg1, hol_term* arg2) {
	return hol_term::new_apply(hol_term::new_constant(predicate), arg1, arg2);
}

inline hol_term* hol_term::new_atom(unsigned int predicate, hol_term* arg1) {
	return hol_term::new_apply(hol_term::new_constant(predicate), arg1);

}

inline hol_term* hol_term::new_atom(unsigned int predicate) {
	return hol_term::new_constant(predicate);
}

inline hol_term* hol_term::new_true() {
	HOL_TRUE.reference_count++;
	return &HOL_TRUE;
}

inline hol_term* hol_term::new_false() {
	HOL_FALSE.reference_count++;
	return &HOL_FALSE;
}

template<hol_term_type Operator, unsigned int Index>
inline void new_hol_array_helper(hol_term** operands, hol_term* arg)
{
	operands[Index] = arg;
}

template<hol_term_type Operator, unsigned int Index, typename... Args>
inline void new_hol_array_helper(hol_term** operands, hol_term* arg, Args&&... args)
{
	operands[Index] = arg;
	new_hol_array_helper<Operator, Index + 1>(operands, std::forward<Args>(args)...);
}

template<hol_term_type Operator, typename... Args>
inline hol_term* new_hol_array(hol_term* arg, Args&&... args)
{
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = Operator;
	term->array.length = 1 + sizeof...(Args);
	term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * (1 + sizeof...(Args)));
	if (term->array.operands == NULL) {
		free(term); return NULL;
	}
	new_hol_array_helper<Operator, 0>(term->array.operands, arg, std::forward<Args>(args)...);
#if !defined(NDEBUG)
	if (term->array.length < 2)
		fprintf(stderr, "new_hol_array WARNING: Array length is not at least 2.\n");
#endif
	return term;
}

template<typename Array,
	typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type** = nullptr>
inline bool new_hol_array_term(hol_array_term& array, const Array& operands)
{
	array.length = operands.size();
	if (array.length == 0) {
		array.operands = nullptr;
		return true;
	}
	array.operands = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * operands.size()));
	if (array.operands == nullptr)
		return false;
	for (unsigned int i = 0; i < operands.size(); i++) {
		if (operands[i] == nullptr) return false;
		array.operands[i] = operands[i];
	}
	return true;
}

template<hol_term_type Operator, typename Array,
	typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type** = nullptr>
inline hol_term* new_hol_array(const Array& operands)
{
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = Operator;
	if (!new_hol_array_term(term->array, operands)) {
		free(term);
		return nullptr;
	}
#if !defined(NDEBUG)
	if (term->array.length < 2)
		fprintf(stderr, "new_hol_array WARNING: Array length is not at least 2.\n");
#endif
	return term;
}

template<typename... Args>
inline hol_term* hol_term::new_and(Args&&... args) {
	return new_hol_array<hol_term_type::AND>(std::forward<Args>(args)...);
}

template<typename... Args>
inline hol_term* hol_term::new_or(Args&&... args) {
	return new_hol_array<hol_term_type::OR>(std::forward<Args>(args)...);
}

template<typename Array,
	typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type**>
inline hol_term* hol_term::new_and(Array& args) {
	return new_hol_array<hol_term_type::AND>(args);
}

template<typename Array,
	typename std::enable_if<has_index_operator<Array, hol_term*>::value>::type**>
inline hol_term* hol_term::new_or(Array& args) {
	return new_hol_array<hol_term_type::OR>(args);
}

template<hol_term_type Type>
inline hol_term* new_hol_binary_term(hol_term* first, hol_term* second) {
	if (first == NULL || second == NULL)
		return NULL;

	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = Type;
	term->binary.left = first;
	term->binary.right = second;
	return term;
}

template<hol_term_type Type>
inline hol_term* new_hol_ternary_term(hol_term* first, hol_term* second, hol_term* third) {
	if (first == NULL || second == NULL)
		return NULL;

	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = Type;
	term->ternary.first = first;
	term->ternary.second = second;
	term->ternary.third = third;
	return term;
}

hol_term* hol_term::new_if_then(hol_term* first, hol_term* second) {
	return new_hol_binary_term<hol_term_type::IF_THEN>(first, second);
}

hol_term* hol_term::new_equals(hol_term* first, hol_term* second) {
	return new_hol_binary_term<hol_term_type::EQUALS>(first, second);
}

inline hol_term* new_equals_helper(hol_term* first, hol_term* second) {
	return new_hol_binary_term<hol_term_type::EQUALS>(first, second);
}

template<typename... Args>
inline hol_term* new_equals_helper(hol_term* first, Args&&... args) {
	return new_hol_binary_term<hol_term_type::EQUALS>(first, new_equals_helper(std::forward<Args>(args)...));
}

template<typename... Args>
hol_term* hol_term::new_iff(Args&&... args) {
	return new_equals_helper(std::forward<Args>(args)...);
}

hol_term* hol_term::new_apply(hol_term* function, hol_term* arg) {
	return new_hol_binary_term<hol_term_type::UNARY_APPLICATION>(function, arg);
}

hol_term* hol_term::new_apply(hol_term* function, hol_term* arg1, hol_term* arg2) {
	return new_hol_ternary_term<hol_term_type::BINARY_APPLICATION>(function, arg1, arg2);
}

hol_term* hol_term::new_not(hol_term* operand)
{
	if (operand == NULL) return NULL;

	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::NOT;
	term->unary.operand = operand;
	return term;
}

template<hol_term_type QuantifierType>
hol_term* new_hol_quantifier(unsigned int variable, hol_term* operand)
{
	if (operand == NULL) return NULL;

	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = QuantifierType;
	term->quantifier.variable = variable;
	term->quantifier.operand = operand;
	return term;
}

inline hol_term* hol_term::new_for_all(unsigned int variable, hol_term* operand) {
	return new_hol_quantifier<hol_term_type::FOR_ALL>(variable, operand);
}

inline hol_term* hol_term::new_exists(unsigned int variable, hol_term* operand) {
	return new_hol_quantifier<hol_term_type::EXISTS>(variable, operand);
}

inline hol_term* hol_term::new_lambda(unsigned int variable, hol_term* operand) {
	return new_hol_quantifier<hol_term_type::LAMBDA>(variable, operand);
}

template<hol_term_type Type>
inline hol_term* new_hol_any(hol_term* included,
		hol_term** excluded_trees, unsigned int excluded_tree_count)
{
	hol_term* term;
	if (!new_hol_term(term)) return nullptr;
	term->reference_count = 1;
	term->type = Type;
	term->any.excluded_tree_count = excluded_tree_count;
	if (excluded_tree_count == 0) {
		term->any.excluded_trees = nullptr;
	} else {
		term->any.excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * excluded_tree_count);
		if (term->any.excluded_trees == nullptr) {
			fprintf(stderr, "hol_term.new_any ERROR: Insufficient memory for `excluded_trees`.\n");
			core::free(term); return nullptr;
		}
	}
	term->any.included = included;
	for (unsigned int i = 0; i < term->any.excluded_tree_count; i++)
		term->any.excluded_trees[i] = excluded_trees[i];
	return term;
}

hol_term* hol_term::new_any(hol_term* included) {
	return new_hol_any<hol_term_type::ANY>(included, nullptr, 0);
}

hol_term* hol_term::new_any(hol_term* included,
		hol_term** excluded_trees, unsigned int excluded_tree_count)
{
	return new_hol_any<hol_term_type::ANY>(included, excluded_trees, excluded_tree_count);
}

inline hol_term* hol_term::new_any_right(hol_term* included) {
	return new_hol_any<hol_term_type::ANY_RIGHT>(included, nullptr, 0);
}

hol_term* hol_term::new_any_right(hol_term* included,
		hol_term** excluded_trees, unsigned int excluded_tree_count)
{
	return new_hol_any<hol_term_type::ANY_RIGHT>(included, excluded_trees, excluded_tree_count);
}

template<typename AnyArray, typename LeftArray, typename RightArray,
	typename std::enable_if<has_index_operator<AnyArray, hol_term*>::value>::type**,
	typename std::enable_if<has_index_operator<LeftArray, hol_term*>::value>::type**,
	typename std::enable_if<has_index_operator<RightArray, hol_term*>::value>::type**>
hol_term* hol_term::new_any_array(hol_term_type oper, hol_term* all,
		const AnyArray& any, const LeftArray& left, const RightArray& right)
{
	if (all == nullptr)
		return nullptr;

	hol_term* term;
	if (!new_hol_term(term)) return nullptr;
	term->reference_count = 1;
	term->type = hol_term_type::ANY_ARRAY;
	term->any_array.oper = oper;
	term->any_array.all = all;
	if (!new_hol_array_term(term->any_array.any, any)
	 || !new_hol_array_term(term->any_array.left, left)
	 || !new_hol_array_term(term->any_array.right, right))
	{
		core::free(term);
		return nullptr;
	}
	return term;
}

template<unsigned int Index>
inline void new_any_constant_helper(unsigned int* constants, unsigned int constant)
{
	constants[Index] = constant;
}

template<unsigned int Index, typename... Args>
inline void new_any_constant_helper(unsigned int* constants, unsigned int constant, Args&&... args)
{
	constants[Index] = constant;
	new_any_constant_helper<Index + 1>(constants, std::forward<Args>(args)...);
}

template<typename... Args>
hol_term* hol_term::new_any_constant(unsigned int constant, Args&&... args) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::ANY_CONSTANT;
	term->any_constant.length = 1 + sizeof...(Args);
	term->any_constant.constants = (unsigned int*) malloc(sizeof(unsigned int) * (1 + sizeof...(Args)));
	if (term->any_constant.constants == nullptr) {
		core::free(term); return nullptr;
	}
	new_any_constant_helper<0>(term->any_constant.constants, constant, std::forward<Args>(args)...);
#if !defined(NDEBUG)
	if (!std::is_sorted(term->any_constant.constants, term->any_constant.constants + term->any_constant.length))
		fprintf(stderr, "hol_term.new_any_constant WARNING: The given constant array is not sorted.\n");
#endif
	return term;
}

template<typename Array,
	typename std::enable_if<has_index_operator<Array, unsigned int>::value>::type**>
hol_term* hol_term::new_any_constant(const Array& constants) {
	hol_term* term;
	if (!new_hol_term(term)) return NULL;
	term->reference_count = 1;
	term->type = hol_term_type::ANY_CONSTANT;
	term->any_constant.length = constants.size();
	term->any_constant.constants = (unsigned int*) malloc(sizeof(unsigned int) * constants.size());
	if (term->any_constant.constants == nullptr) {
		core::free(term); return nullptr;
	}
	for (unsigned int i = 0; i < term->any_constant.length; i++)
		term->any_constant.constants[i] = constants[i];
#if !defined(NDEBUG)
	if (!std::is_sorted(term->any_constant.constants, term->any_constant.constants + term->any_constant.length))
		fprintf(stderr, "hol_term.new_any_constant WARNING: The given constant array is not sorted.\n");
#endif
	return term;
}

hol_term* hol_term::new_any_quantifier(hol_quantifier_type quantifier_type, hol_term* operand)
{
	if (operand == nullptr) return nullptr;

	hol_term* term;
	if (!new_hol_term(term)) return nullptr;
	term->reference_count = 1;
	term->type = hol_term_type::ANY_QUANTIFIER;
	term->any_quantifier.quantifier = quantifier_type;
	term->any_quantifier.operand = operand;
	return term;
}



/**
 * Below is code for computing the type of higher-order terms.
 */


/* forward declarations */
struct hol_type;
bool unify(const hol_type&, const hol_type&, hol_type&, array<hol_type>&);
template<bool Root> bool flatten_type_variable(hol_type&,
		array<pair<unsigned int, bool>>&, array<hol_type>&);


enum class hol_type_kind {
	CONSTANT,
	FUNCTION,
	VARIABLE,
	ANY,
	NONE
};

enum class hol_constant_type {
	BOOLEAN,
	INDIVIDUAL
};

struct hol_function_type {
	hol_type* left;
	hol_type* right;
};

struct hol_type {
	hol_type_kind kind;
	union {
		hol_constant_type constant;
		hol_function_type function;
		unsigned int variable;
	};

	hol_type() { }

	explicit hol_type(hol_type_kind kind) : kind(kind) { }

	explicit hol_type(hol_constant_type constant) : kind(hol_type_kind::CONSTANT), constant(constant) { }

	explicit hol_type(unsigned int variable) : kind(hol_type_kind::VARIABLE), variable(variable) { }

	explicit hol_type(const hol_type& type) {
		if (!init_helper(type)) exit(EXIT_FAILURE);
	}

	explicit hol_type(const hol_type& left, const hol_type& right) : kind(hol_type_kind::FUNCTION) {
		if (!init_helper(left, right)) exit(EXIT_FAILURE);
	}

	~hol_type() { free_helper(); }

	inline void operator = (const hol_type& src) {
		init_helper(src);
	}

	static inline void move(const hol_type& src, hol_type& dst) {
		dst.kind = src.kind;
		switch (src.kind) {
		case hol_type_kind::CONSTANT:
			dst.constant = src.constant; return;
		case hol_type_kind::VARIABLE:
			dst.variable = src.variable; return;
		case hol_type_kind::FUNCTION:
			dst.function.left = src.function.left;
			dst.function.right = src.function.right;
			return;
		case hol_type_kind::ANY:
		case hol_type_kind::NONE:
			return;
		}
		fprintf(stderr, "hol_type.move ERROR: Unrecognized hol_type_kind.\n");
	}

	static inline void swap(hol_type& first, hol_type& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(hol_type); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void free(hol_type& type) {
		type.free_helper();
	}

private:
	inline bool init_helper(const hol_type& src) {
		kind = src.kind;
		switch (src.kind) {
		case hol_type_kind::CONSTANT:
			constant = src.constant; return true;
		case hol_type_kind::VARIABLE:
			variable = src.variable; return true;
		case hol_type_kind::FUNCTION:
			return init_helper(*src.function.left, *src.function.right);
		case hol_type_kind::ANY:
		case hol_type_kind::NONE:
			return true;
		}
		fprintf(stderr, "hol_type.init_helper ERROR: Unrecognized hol_type_kind.\n");
		return false;
	}

	inline bool init_helper(const hol_type& left, const hol_type& right) {
		function.left = (hol_type*) malloc(sizeof(hol_type));
		if (function.left == NULL) {
			fprintf(stderr, "hol_type.init_helper ERROR: Insufficient memory for `function.left`.\n");
			return false;
		}
		function.right = (hol_type*) malloc(sizeof(hol_type));
		if (function.right == NULL) {
			fprintf(stderr, "hol_type.init_helper ERROR: Insufficient memory for `function.right`.\n");
			core::free(function.left); return false;
		} else if (!init(*function.left, left)) {
			core::free(function.left); core::free(function.right);
			return false;
		} else if (!init(*function.right, right)) {
			core::free(*function.left); core::free(function.left);
			core::free(function.right); return false;
		}
		return true;
	}

	inline void free_helper() {
		switch (kind) {
		case hol_type_kind::CONSTANT:
		case hol_type_kind::VARIABLE:
		case hol_type_kind::ANY:
		case hol_type_kind::NONE:
			return;
		case hol_type_kind::FUNCTION:
			core::free(*function.left); core::free(function.left);
			core::free(*function.right); core::free(function.right);
			return;
		}
		fprintf(stderr, "hol_type.free_helper ERROR: Unrecognized hol_type_kind.\n");
	}

	friend bool init(hol_type&, const hol_type&);
	friend bool init(hol_type&, const hol_type&, const hol_type&);
};

inline bool init(hol_type& type, hol_type_kind kind) {
	type.kind = kind;
	return true;
}

inline bool init(hol_type& type, const hol_type& src) {
	return type.init_helper(src);
}

inline bool init(hol_type& type, const hol_type& left, const hol_type& right) {
	type.kind = hol_type_kind::FUNCTION;
	return type.init_helper(left, right);
}

inline bool init(hol_type& type, hol_constant_type constant) {
	type.kind = hol_type_kind::CONSTANT;
	type.constant = constant;
	return true;
}

inline bool init(hol_type& type, unsigned int variable) {
	type.kind = hol_type_kind::VARIABLE;
	type.variable = variable;
	return true;
}

const hol_type HOL_BOOLEAN_TYPE(hol_constant_type::BOOLEAN);
const hol_type HOL_INTEGER_TYPE(hol_constant_type::INDIVIDUAL);
const hol_type HOL_STRING_TYPE(hol_constant_type::INDIVIDUAL);
const hol_type HOL_UINT_LIST_TYPE(hol_constant_type::INDIVIDUAL);

inline bool operator == (const hol_type& first, const hol_type& second) {
	if (first.kind != second.kind) return false;
	switch (first.kind) {
	case hol_type_kind::ANY:
	case hol_type_kind::NONE:
		return true;
	case hol_type_kind::CONSTANT:
		return first.constant == second.constant;
	case hol_type_kind::VARIABLE:
		return first.variable == second.variable;
	case hol_type_kind::FUNCTION:
		return *first.function.left == *second.function.left
			&& *first.function.right == *second.function.right;
	}
	fprintf(stderr, "operator == ERROR: Unexpected hol_type_kind.\n");
	return false;
}

inline bool operator != (const hol_type& first, const hol_type& second) {
	return !(first == second);
}

template<typename Stream>
inline bool print(const hol_constant_type& constant, Stream& out) {
	switch (constant) {
	case hol_constant_type::BOOLEAN:
		return print("𝝄", out);
	case hol_constant_type::INDIVIDUAL:
		return print("𝜾", out);
	}
	fprintf(stderr, "print ERROR: Unrecognized hol_constant_type.\n");
	return false;
}

template<typename Stream>
inline bool print_variable(unsigned int variable, Stream& out) {
	return print_variable<hol_term_syntax::CLASSIC>(variable, out);
}

template<typename Stream>
inline bool print_variable(unsigned int variable, Stream& out, array<unsigned int>& type_variables) {
	return print_variable(variable, out)
		&& (type_variables.contains(variable) || type_variables.add(variable));
}

template<typename Stream, typename... Printer>
bool print(const hol_type& type, Stream& out, Printer&&... variable_printer)
{
	switch (type.kind) {
	case hol_type_kind::CONSTANT:
		return print(type.constant, out);
	case hol_type_kind::FUNCTION:
		return print('(', out) && print(*type.function.left, out, std::forward<Printer>(variable_printer)...)
			&& print(" → ", out) && print(*type.function.right, out, std::forward<Printer>(variable_printer)...)
			&& print(')', out);
	case hol_type_kind::VARIABLE:
		return print_variable(type.variable, out, std::forward<Printer>(variable_printer)...);
	case hol_type_kind::ANY:
		return print('*', out);
	case hol_type_kind::NONE:
		return print("NONE", out);
	}
	fprintf(stderr, "print ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

template<typename Stream>
bool print_type(const hol_type& type, Stream& out,
		const array<hol_type>& type_variables)
{
	array<unsigned int> variables(8);
	if (!print(type, out, variables)) return false;

	if (variables.length > 0 && !print(" where ", out))
		return false;

	for (unsigned int i = 0; i < variables.length; i++) {
		if (i > 0 && !print(", ", out))
			return false;
		if (!print_variable<hol_term_syntax::CLASSIC>(variables[i], out)
		 || !print(" = ", out) || !print(type_variables[variables[i]], out, variables))
		{
			return false;
		}
	}
	return true;
}

inline bool unify_constant(
		hol_constant_type first, const hol_type& second,
		hol_type& out, array<hol_type>& type_variables)
{
	if (second.kind == hol_type_kind::ANY) {
		return init(out, first);
	} else if (second.kind == hol_type_kind::CONSTANT && second.constant == first) {
		return init(out, first);
	} else if (second.kind == hol_type_kind::VARIABLE) {
		if (!unify_constant(first, type_variables[second.variable], out, type_variables))
			return false;
		if (out.kind != hol_type_kind::NONE) {
			free(type_variables[second.variable]);
			return init(type_variables[second.variable], out);
		}
	}

	return init(out, hol_type_kind::NONE);
}

inline bool unify_function(
		const hol_function_type& first, const hol_type& second,
		hol_type& out, array<hol_type>& type_variables)
{
	if (second.kind == hol_type_kind::ANY) {
		return init(out, *first.left, *first.right);
	} else if (second.kind == hol_type_kind::VARIABLE) {
		if (!unify_function(first, type_variables[second.variable], out, type_variables))
			return false;
		if (out.kind == hol_type_kind::NONE) return false;
		free(type_variables[second.variable]);
		return init(type_variables[second.variable], out);
	} else if (second.kind != hol_type_kind::FUNCTION) {
		return init(out, hol_type_kind::NONE);
	}

	out.function.left = (hol_type*) malloc(sizeof(hol_type));
	if (out.function.left == NULL) {
		fprintf(stderr, "unify ERROR: Insufficient memory for `out.function.left`.\n");
		return false;
	}
	out.function.right = (hol_type*) malloc(sizeof(hol_type));
	if (out.function.right == NULL) {
		fprintf(stderr, "unify ERROR: Insufficient memory for `out.function.right`.\n");
		free(out.function.left); return false;
	}

	if (!unify(*first.left, *second.function.left, *out.function.left, type_variables)) {
		free(out.function.left); free(out.function.right);
		return false;
	} else if (out.function.left->kind == hol_type_kind::NONE) {
		free(*out.function.left); free(out.function.left); free(out.function.right);
		out.kind = hol_type_kind::NONE; return true;
	}

	if (!unify(*first.right, *second.function.right, *out.function.right, type_variables)) {
		free(*out.function.left); free(out.function.left);
		free(out.function.right); return false;
	} else if (out.function.right->kind == hol_type_kind::NONE) {
		free(*out.function.left); free(out.function.left);
		free(*out.function.right); free(out.function.right);
		out.kind = hol_type_kind::NONE; return true;
	}

	out.kind = hol_type_kind::FUNCTION;
	return true;
}

inline bool unify_variable(
		unsigned int first, const hol_type& second,
		hol_type& out, array<hol_type>& type_variables)
{
	unsigned int var;
	switch (second.kind) {
	case hol_type_kind::ANY:
		return init(out, first);
	case hol_type_kind::NONE:
		return init(out, hol_type_kind::NONE);
	case hol_type_kind::CONSTANT:
		if (!unify_constant(second.constant, type_variables[first], out, type_variables))
			return false;
		if (out.kind == hol_type_kind::NONE)
			return init(out, hol_type_kind::NONE);
		free(type_variables[first]);
		return init(type_variables[first], out);
	case hol_type_kind::FUNCTION:
		if (!unify_function(second.function, type_variables[first], out, type_variables))
			return false;
		if (out.kind == hol_type_kind::NONE)
			return init(out, hol_type_kind::NONE);
		free(type_variables[first]);
		return init(type_variables[first], out);
	case hol_type_kind::VARIABLE:
		var = second.variable;
		if (first == var) return init(out, var);
		while (type_variables[var].kind == hol_type_kind::VARIABLE) {
			var = type_variables[var].variable;
			if (first == var) return init(out, var);
		}

		if (!unify_variable(first, type_variables[var], out, type_variables))
			return false;
		if (out.kind == hol_type_kind::NONE)
			return init(out, hol_type_kind::NONE);
		free(type_variables[var]);
		move(out, type_variables[var]);
		return init(out, var);
	}
	fprintf(stderr, "unify_variable ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

bool unify(const hol_type& first, const hol_type& second,
		hol_type& out, array<hol_type>& type_variables)
{
	switch (first.kind) {
	case hol_type_kind::ANY:
		return init(out, second);
	case hol_type_kind::NONE:
		return init(out, hol_type_kind::NONE);
	case hol_type_kind::CONSTANT:
		return unify_constant(first.constant, second, out, type_variables);
	case hol_type_kind::FUNCTION:
		return unify_function(first.function, second, out, type_variables);
	case hol_type_kind::VARIABLE:
		return unify_variable(first.variable, second, out, type_variables);
	}
	fprintf(stderr, "unify ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

inline bool expect_type(
		const hol_type& actual_type,
		hol_type& expected_type,
		array<hol_type>& type_variables)
{
	hol_type& temp = *((hol_type*) alloca(sizeof(hol_type)));
	if (!unify(actual_type, expected_type, temp, type_variables))
		return false;
	swap(expected_type, temp); free(temp);
	if (expected_type.kind == hol_type_kind::NONE) {
		print("ERROR: Term is not well-typed.\n", stderr);
		print("  Computed type: ", stderr); print_type(actual_type, stderr, type_variables); print('\n', stderr);
		print("  Expected type: ", stderr); print_type(expected_type, stderr, type_variables); print('\n', stderr);
		return false;
	}
	return true;
}

template<hol_term_type Type, typename ComputedTypes>
inline bool compute_type(
		unsigned int symbol, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& symbol_types,
		array<hol_type>& type_variables)
{
	static_assert(Type == hol_term_type::CONSTANT
			   || Type == hol_term_type::VARIABLE
			   || Type == hol_term_type::PARAMETER,
			"Type must be either: CONSTANT, VARIABLE, PARAMETER.");

	if (!types.template push<Type>(term)
	 || !symbol_types.ensure_capacity(symbol_types.size + 1))
		return false;

	unsigned int index = symbol_types.index_of(symbol);
	if (index < symbol_types.size) {
		hol_type& new_type = *((hol_type*) alloca(sizeof(hol_type)));
		if (!unify(expected_type, symbol_types.values[index], new_type, type_variables))
			return false;
		if (new_type.kind == hol_type_kind::NONE) {
			print("ERROR: Term is not well-typed. Symbol ", stderr); print(symbol, stderr);
			print(" has conflicting types: ", stderr);
			print("  Type computed from earlier instances of symbol: ", stderr);
			print_type(symbol_types.values[index], stderr, type_variables); print('\n', stderr);
			print("  Expected type: ", stderr); print_type(expected_type, stderr, type_variables); print('\n', stderr);
			free(new_type); return false;
		} else {
			swap(new_type, symbol_types.values[index]);
			free(new_type); free(expected_type);
			if (!init(expected_type, symbol_types.values[index])) {
				expected_type.kind = hol_type_kind::ANY; /* make sure `expected_type` is valid since its destructor will be called */
				return false;
			}
			return types.template add<Type>(term, symbol_types.values[index]);
		}
	} else {
		symbol_types.keys[index] = symbol;
		if (!init(symbol_types.values[index], expected_type))
			return false;
		symbol_types.size++;
		return types.template add<Type>(term, expected_type);
	}
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_array_type(const hol_array_term& array_term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	for (unsigned int i = 0; i < array_term.length; i++) {
		if (!compute_type<PolymorphicEquality>(*array_term.operands[i], types, expected_type,
				constant_types, variable_types, parameter_types, type_variables))
		{
			return false;
		}
	}
	return true;
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_equals_type(
		const hol_binary_term& equals, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	unsigned int type_variable = type_variables.length;
	if (!types.template push<hol_term_type::EQUALS>(term)
		|| !expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
		|| !type_variables.ensure_capacity(type_variables.length + 1)
		|| !init(type_variables[type_variable], hol_type_kind::ANY))
		return false;
	type_variables.length++;

	hol_type first_type(type_variable);
	if (!compute_type<PolymorphicEquality>(*equals.left, types,
			first_type, constant_types, variable_types, parameter_types, type_variables))
		return false;
	if (PolymorphicEquality) {
		type_variable = type_variables.length;
		if (!type_variables.ensure_capacity(type_variables.length + 1)
		 || !init(type_variables[type_variables.length], hol_type_kind::ANY))
			return false;
		type_variables.length++;
	}

	hol_type second_type(type_variable);
	return compute_type<PolymorphicEquality>(*equals.right, types,
			second_type, constant_types, variable_types, parameter_types, type_variables)
		&& types.template add<hol_term_type::EQUALS>(term, HOL_BOOLEAN_TYPE, first_type, second_type);
}

template<hol_term_type Type, bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_type(
		const hol_quantifier& quantifier, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	static_assert(Type == hol_term_type::FOR_ALL
			   || Type == hol_term_type::EXISTS,
			"Type must be either: FOR_ALL, EXISTS.");

	if (!types.template push<Type>(term)
	 || !variable_types.ensure_capacity(variable_types.size + 1)
	 || !type_variables.ensure_capacity(type_variables.length + 1)
	 || !expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables))
	{
		return false;
	}

	unsigned int new_type_variable = type_variables.length;
	if (!init(type_variables[new_type_variable], hol_type_kind::ANY))
		return false;
	type_variables.length++;

#if !defined(NDEBUG)
	if (variable_types.contains(quantifier.variable))
		fprintf(stderr, "compute_type WARNING: `variable_types` already contains key %u.\n", quantifier.variable);
	unsigned int old_variable_types_size = variable_types.size;
#endif
	variable_types.keys[variable_types.size] = quantifier.variable;
	if (!init(variable_types.values[variable_types.size], new_type_variable))
		return false;
	variable_types.size++;

	if (!compute_type<PolymorphicEquality>(*quantifier.operand, types, expected_type,
	 		constant_types, variable_types, parameter_types, type_variables))
	{
		return false;
	}

#if !defined(NDEBUG)
	if (old_variable_types_size + 1 != variable_types.size
	 || variable_types.keys[variable_types.size - 1] != quantifier.variable)
		fprintf(stderr, "compute_type WARNING: Quantified term is not well-formed.\n");
#endif
	variable_types.size--;
	if (!types.template add<Type>(term, HOL_BOOLEAN_TYPE, variable_types.values[variable_types.size])) {
		free(variable_types.values[variable_types.size]);
		return false;
	}
	free(variable_types.values[variable_types.size]);
	return true;
}

inline bool get_function_child_types(unsigned int type_variable,
		array<hol_type>& type_variables, hol_type& left, hol_type& right)
{
	switch (type_variables[type_variable].kind) {
	case hol_type_kind::ANY:
		if (!type_variables.ensure_capacity(type_variables.length + 2))
			return false;
		if (!init(type_variables[type_variables.length], hol_type_kind::ANY)) return false;
		type_variables.length++;
		if (!init(type_variables[type_variables.length], hol_type_kind::ANY)) return false;
		type_variables.length++;
		if (!init(left, type_variables.length - 2)) return false;
		if (!init(right, type_variables.length - 1)) { free(left); return false; }
		free(type_variables[type_variable]);
		if (!init(type_variables[type_variable], left, right)) {
			free(left); free(right);
			return false;
		}
		return true;
	case hol_type_kind::FUNCTION:
		if (!init(left, *type_variables[type_variable].function.left)) return false;
		if (!init(right, *type_variables[type_variable].function.right)) { free(left); return false; }
		return true;
	case hol_type_kind::VARIABLE:
		if (!get_function_child_types(type_variables[type_variable].variable, type_variables, left, right)) return false;
		free(type_variables[type_variable]);
		if (!init(type_variables[type_variable], left, right)) {
			free(left); free(right);
			return false;
		}
		return true;
	case hol_type_kind::NONE:
	case hol_type_kind::CONSTANT:
		left.kind = hol_type_kind::NONE;
		return true;
	}
	fprintf(stderr, "get_function_child_types ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

inline bool get_function_child_types(const hol_type& type,
		array<hol_type>& type_variables, hol_type& left, hol_type& right)
{
	switch (type.kind) {
	case hol_type_kind::ANY:
		fprintf(stderr, "get_function_child_types ERROR: `type` is not supposed to be ANY.\n");
		return false;
	case hol_type_kind::FUNCTION:
		if (!init(left, *type.function.left)) return false;
		if (!init(right, *type.function.right)) { free(left); return false; }
		return true;
	case hol_type_kind::VARIABLE:
		return get_function_child_types(type.variable, type_variables, left, right);
	case hol_type_kind::NONE:
	case hol_type_kind::CONSTANT:
		left.kind = hol_type_kind::NONE;
		return true;
	}
	fprintf(stderr, "get_function_child_types ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_lambda_type(
		const hol_quantifier& quantifier, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	if (!types.template push<hol_term_type::LAMBDA>(term)
	 || !variable_types.ensure_capacity(variable_types.size + 1)) return false;

	hol_type& left_type = *((hol_type*) alloca(sizeof(hol_type)));
	hol_type& right_type = *((hol_type*) alloca(sizeof(hol_type)));
	if (!get_function_child_types(expected_type, type_variables, left_type, right_type))
		return false;
	if (left_type.kind == hol_type_kind::NONE) {
		print("ERROR: Term is not well-typed. Lambda expression has a non-function expected type: ", stderr);
		print_type(expected_type, stderr, type_variables); print(".\n", stderr); free(left_type);
		return false;
	}

#if !defined(NDEBUG)
	if (variable_types.contains(quantifier.variable))
		fprintf(stderr, "compute_lambda_type WARNING: `variable_types` already contains key %u.\n", quantifier.variable);
	unsigned int old_variable_types_size = variable_types.size;
#endif
	variable_types.keys[variable_types.size] = quantifier.variable;
	move(left_type, variable_types.values[variable_types.size]);
	variable_types.size++;

	if (!compute_type<PolymorphicEquality>(*quantifier.operand, types, right_type, constant_types, variable_types, parameter_types, type_variables))
		return false;

#if !defined(NDEBUG)
	if (old_variable_types_size + 1 != variable_types.size
	 || variable_types.keys[variable_types.size - 1] != quantifier.variable)
		fprintf(stderr, "compute_lambda_type WARNING: Lambda term is not well-formed.\n");
#endif
	variable_types.size--;
	free(expected_type);
	if (!init(expected_type, variable_types.values[variable_types.size], right_type)) {
		expected_type.kind = hol_type_kind::ANY; /* make sure `expected_type` is valid since its destructor will be called */
		free(variable_types.values[variable_types.size]); free(right_type);
		return false;
	}
	free(variable_types.values[variable_types.size]); free(right_type);
	return types.template add<hol_term_type::LAMBDA>(term, expected_type);
}

template<bool PolymorphicEquality, hol_term_type AnyType, typename Any, typename ComputedTypes>
inline bool compute_any_type(
		const Any& any, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	if (!types.template push<AnyType>(term))
		return false;

	if (any.included != nullptr) {
		hol_type included_type(hol_type_kind::ANY);
		if (!compute_type<PolymorphicEquality>(*any.included, types, included_type, constant_types, variable_types, parameter_types, type_variables))
			return false;
	}

	return types.template add<AnyType>(term, expected_type);
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_apply_type(
		const hol_binary_term& apply, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	if (!types.template push<hol_term_type::UNARY_APPLICATION>(term)
	 || !type_variables.ensure_capacity(type_variables.length + 1)
	 || !init(type_variables[type_variables.length], hol_type_kind::ANY))
		return false;
	type_variables.length++;
	hol_type type(hol_type(type_variables.length - 1), expected_type);
	if (!compute_type<PolymorphicEquality>(*apply.left, types, type,
			constant_types, variable_types, parameter_types, type_variables))
	{
		return false;
	}

	hol_type& arg_type = *type.function.left;
	swap(expected_type, *type.function.right);

	return compute_type<PolymorphicEquality>(*apply.right, types, arg_type, constant_types, variable_types, parameter_types, type_variables)
		&& types.template add<hol_term_type::UNARY_APPLICATION>(term, expected_type);
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_apply_type(
		const hol_ternary_term& apply, const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	if (!types.template push<hol_term_type::BINARY_APPLICATION>(term)
	 || !type_variables.ensure_capacity(type_variables.length + 2)
	 || !init(type_variables[type_variables.length], hol_type_kind::ANY))
		return false;
	type_variables.length++;
	if (!init(type_variables[type_variables.length], hol_type_kind::ANY))
		return false;
	type_variables.length++;
	hol_type type(hol_type(type_variables.length - 1), hol_type(hol_type(type_variables.length - 2), expected_type));
	if (!compute_type<PolymorphicEquality>(*apply.first, types, type,
			constant_types, variable_types, parameter_types, type_variables))
	{
		return false;
	}

	hol_type& arg1_type = *type.function.left;
	hol_type& arg2_type = *type.function.right->function.left;
	swap(expected_type, *type.function.right->function.right);

	return compute_type<PolymorphicEquality>(*apply.second, types, arg1_type, constant_types, variable_types, parameter_types, type_variables)
		&& compute_type<PolymorphicEquality>(*apply.third, types, arg2_type, constant_types, variable_types, parameter_types, type_variables)
		&& types.template add<hol_term_type::BINARY_APPLICATION>(term, expected_type);
}

template<bool PolymorphicEquality, typename ComputedTypes>
bool compute_type(const hol_term& term,
		ComputedTypes& types, hol_type& expected_type,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types,
		array<hol_type>& type_variables)
{
	switch (term.type) {
	case hol_term_type::CONSTANT:
		return compute_type<hol_term_type::CONSTANT>(term.constant, term, types, expected_type, constant_types, type_variables);
	case hol_term_type::VARIABLE:
		return compute_type<hol_term_type::VARIABLE>(term.variable, term, types, expected_type, variable_types, type_variables);
	case hol_term_type::PARAMETER:
		return compute_type<hol_term_type::PARAMETER>(term.parameter, term, types, expected_type, parameter_types, type_variables);
	case hol_term_type::INTEGER:
		return types.template push<hol_term_type::INTEGER>(term)
			&& expect_type(HOL_INTEGER_TYPE, expected_type, type_variables)
			&& types.template add<hol_term_type::INTEGER>(term, HOL_INTEGER_TYPE);
	case hol_term_type::STRING:
		return types.template push<hol_term_type::STRING>(term)
			&& expect_type(HOL_STRING_TYPE, expected_type, type_variables)
			&& types.template add<hol_term_type::STRING>(term, HOL_STRING_TYPE);
	case hol_term_type::UINT_LIST:
		return types.template push<hol_term_type::UINT_LIST>(term)
			&& expect_type(HOL_UINT_LIST_TYPE, expected_type, type_variables)
			&& types.template add<hol_term_type::UINT_LIST>(term, HOL_UINT_LIST_TYPE);
	case hol_term_type::UNARY_APPLICATION:
		return compute_apply_type<PolymorphicEquality>(term.binary, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::BINARY_APPLICATION:
		return compute_apply_type<PolymorphicEquality>(term.ternary, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::NOT:
		return types.template push<hol_term_type::NOT>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_type<PolymorphicEquality>(*term.unary.operand, types,
					expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::NOT>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::IF_THEN:
		return types.template push<hol_term_type::IF_THEN>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_type<PolymorphicEquality>(*term.binary.left, types,
					expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& compute_type<PolymorphicEquality>(*term.binary.right, types,
					expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::IF_THEN>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::EQUALS:
		return compute_equals_type<PolymorphicEquality>(term.binary, term, types,
					expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::AND:
		return types.template push<hol_term_type::AND>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.array, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::AND>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::OR:
		return types.template push<hol_term_type::OR>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.array, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::OR>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::IFF:
		return types.template push<hol_term_type::IFF>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.array, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::IFF>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::FOR_ALL:
		return compute_type<hol_term_type::FOR_ALL, PolymorphicEquality>(term.quantifier, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::EXISTS:
		return compute_type<hol_term_type::EXISTS, PolymorphicEquality>(term.quantifier, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::LAMBDA:
		return compute_lambda_type<PolymorphicEquality>(term.quantifier, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::TRUE:
		return types.template push<hol_term_type::TRUE>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& types.template add<hol_term_type::TRUE>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::FALSE:
		return types.template push<hol_term_type::FALSE>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& types.template add<hol_term_type::FALSE>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::ANY:
		return compute_any_type<PolymorphicEquality, hol_term_type::ANY>(term.any, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::ANY_RIGHT:
		return compute_any_type<PolymorphicEquality, hol_term_type::ANY_RIGHT>(term.any, term, types,
				expected_type, constant_types, variable_types, parameter_types, type_variables);
	case hol_term_type::ANY_ARRAY:
		return types.template push<hol_term_type::ANY_ARRAY>(term)
			&& expect_type(HOL_BOOLEAN_TYPE, expected_type, type_variables)
			&& compute_type<PolymorphicEquality>(*term.any_array.all,
					types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.any_array.any, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.any_array.left, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& compute_array_type<PolymorphicEquality>(term.any_array.right, types, expected_type, constant_types, variable_types, parameter_types, type_variables)
			&& types.template add<hol_term_type::ANY_ARRAY>(term, HOL_BOOLEAN_TYPE);
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "compute_type ERROR: hol_term_types `ANY_CONSTANT`, `ANY_QUANTIFIER` unsupported.\n");
		return false;
	}
	fprintf(stderr, "compute_type ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<bool Root>
inline bool flatten_type_variable(
		hol_type& type_variable, unsigned int variable,
		array<pair<unsigned int, bool>>& visited_variables,
		array<hol_type>& type_variables)
{
	bool is_trivial_alias = Root;
	for (unsigned int i = visited_variables.length; i > 0; i--) {
		if (visited_variables[i - 1].key == variable) {
			if (is_trivial_alias) {
				/* we found a cycle of trivial variable references, so all these variables become `ANY` */
				for (unsigned int j = i - 1; j < visited_variables.length; j++) {
					free(type_variables[visited_variables[j].key]);
					if (!init(type_variables[visited_variables[j].key], hol_type_kind::ANY)) return false;
				}
				free(type_variables[variable]);
				if (!init(type_variables[variable], hol_type_kind::ANY)) return false;
				free(type_variable);
				return init(type_variable, hol_type_kind::ANY);
			} else {
				/* we found a non-trival cycle of variable references */
				print("flatten_type_variable ERROR: Found infinite type ", stderr);
				print_type(type_variable, stderr, type_variables); print('\n', stderr);
				return false;
			}
		}
		is_trivial_alias &= visited_variables[i - 1].value;
	}

#if !defined(NDEBUG)
	unsigned int old_visited_variable_count = visited_variables.length;
#endif
	if (!visited_variables.add(make_pair(variable, Root))
	 || !flatten_type_variable<true>(type_variables[variable], visited_variables, type_variables))
		return false;
	visited_variables.length--;
#if !defined(NDEBUG)
	if (old_visited_variable_count != visited_variables.length)
		fprintf(stderr, "flatten_type_variable ERROR: `visited_variables` is invalid.\n");
#endif

	/* replace the current variable with its value */
	free(type_variable);
	return init(type_variable, type_variables[variable]);
}

template<bool Root>
bool flatten_type_variable(hol_type& type_variable,
		array<pair<unsigned int, bool>>& visited_variables,
		array<hol_type>& type_variables)
{
	switch (type_variable.kind) {
	case hol_type_kind::ANY:
	case hol_type_kind::NONE:
	case hol_type_kind::CONSTANT:
		return true;
	case hol_type_kind::FUNCTION:
		return flatten_type_variable<false>(*type_variable.function.left, visited_variables, type_variables)
			&& flatten_type_variable<false>(*type_variable.function.right, visited_variables, type_variables);
	case hol_type_kind::VARIABLE:
		return flatten_type_variable<Root>(type_variable, type_variable.variable, visited_variables, type_variables);
	}
	fprintf(stderr, "flatten_type_variable ERROR: Unrecognized hol_type_kind.\n");
	return false;
}

template<typename T>
inline void free_elements(array<T>& list) {
	for (T& element : list)
		free(element);
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_type(
		const hol_term& term, ComputedTypes& types,
		array_map<unsigned int, hol_type>& constant_types,
		array_map<unsigned int, hol_type>& variable_types,
		array_map<unsigned int, hol_type>& parameter_types)
{
	array<hol_type> type_variables(8);
	if (!init(type_variables[0], hol_type_kind::ANY)) return false;
	type_variables.length++;
	hol_type type(0);
	if (!compute_type<PolymorphicEquality>(term, types, type, constant_types, variable_types, parameter_types, type_variables))
		return false;

	array<pair<unsigned int, bool>> visited_variables(8);
	if (!types.apply([&](hol_type& type){ return flatten_type_variable<true>(type, visited_variables, type_variables); })) {
		free_elements(type_variables); return false;
	}
	for (auto entry : constant_types) {
		if (!flatten_type_variable<true>(entry.value, visited_variables, type_variables)) {
			free_elements(type_variables); return false;
		}
	} for (auto entry : variable_types) {
		if (!flatten_type_variable<true>(entry.value, visited_variables, type_variables)) {
			free_elements(type_variables); return false;
		}
	} for (auto entry : parameter_types) {
		if (!flatten_type_variable<true>(entry.value, visited_variables, type_variables)) {
			free_elements(type_variables); return false;
		}
	}
	free_elements(type_variables);
	return true;
}

template<bool PolymorphicEquality, typename ComputedTypes>
inline bool compute_type(
		const hol_term& term, ComputedTypes& types)
{
	array_map<unsigned int, hol_type> constant_types(8);
	array_map<unsigned int, hol_type> variable_types(8);
	array_map<unsigned int, hol_type> parameter_types(8);
	bool success = compute_type<PolymorphicEquality>(term, types, constant_types, variable_types, parameter_types);
	for (unsigned int j = 0; j < constant_types.size; j++) free(constant_types.values[j]);
	for (unsigned int j = 0; j < variable_types.size; j++) free(variable_types.values[j]);
	for (unsigned int j = 0; j < parameter_types.size; j++) free(parameter_types.values[j]);
	return success;
}

struct type_map {
	array_map<const hol_term*, hol_type> types;

	type_map(unsigned int initial_capacity) : types(initial_capacity) { }

	~type_map() {
		for (auto entry : types) free(entry.value);
	}

	template<hol_term_type Type>
	constexpr bool push(const hol_term& term) const { return true; }

	template<hol_term_type Type, typename... Args>
	inline bool add(const hol_term& term, const hol_type& type, Args&&... extra_types) {
		return types.put(&term, type);
	}

	template<typename Function>
	inline bool apply(Function function) {
		for (auto entry : types)
			if (!function(entry.value)) return false;
		return true;
	}

	inline void clear() {
		for (auto entry : types) free(entry.value);
		types.clear();
	}
};

struct equals_arg_types {
	array_map<const hol_term*, pair<hol_type, hol_type>> types;

	equals_arg_types(unsigned int initial_capacity) : types(initial_capacity) { }

	~equals_arg_types() { free_helper(); }

	template<hol_term_type Type>
	constexpr bool push(const hol_term& term) const { return true; }

	template<hol_term_type Type, typename... Args,
		typename std::enable_if<Type != hol_term_type::EQUALS>::type* = nullptr>
	constexpr bool add(const hol_term& term, const hol_type& type, Args&&... extra_types) const {
		return true;
	}

	template<hol_term_type Type,
		typename std::enable_if<Type == hol_term_type::EQUALS>::type* = nullptr>
	inline bool add(const hol_term& term, const hol_type& type,
			const hol_type& first_arg_type, const hol_type& second_arg_type)
	{
		if (!types.ensure_capacity(types.size + 1)) return false;

		unsigned int index;
		pair<hol_type, hol_type>& arg_types = types.get(&term, index);
		if (index == types.size) {
			if (!init(arg_types.key, first_arg_type)) {
				return false;
			} else if (!init(arg_types.value, second_arg_type)) {
				free(arg_types.key);
				return false;
			}
			types.keys[index] = &term;
			types.size++;
			return true;
		} else {
			fprintf(stderr, "equals_arg_types.add ERROR: We've already seen this term.\n");
			return false;
		}
	}

	template<typename Function>
	inline bool apply(Function function) {
		for (auto entry : types)
			if (!function(entry.value.key)
			 || !function(entry.value.value)) return false;
		return true;
	}

	inline void clear() {
		free_helper();
		types.clear();
	}

private:
	void free_helper() {
		for (auto entry : types) {
			free(entry.value.key);
			free(entry.value.value);
		}
	}
};


/**
 * Below is code for canonicalizing higher-order formulas.
 */


int_fast8_t compare(
		const hol_term&,
		const hol_term&);


inline int_fast8_t compare(
		const hol_unary_term& first,
		const hol_unary_term& second)
{
	return compare(*first.operand, *second.operand);
}

inline int_fast8_t compare(
		const hol_binary_term& first,
		const hol_binary_term& second)
{
	int_fast8_t result = compare(*first.left, *second.left);
	if (result != 0) return result;
	return compare(*first.right, *second.right);
}

inline int_fast8_t compare(
		const hol_ternary_term& first,
		const hol_ternary_term& second)
{
	int_fast8_t result = compare(*first.first, *second.first);
	if (result != 0) return result;
	result = compare(*first.second, *second.second);
	if (result != 0) return result;
	return compare(*first.third, *second.third);
}

inline int_fast8_t compare(
		const hol_array_term& first,
		const hol_array_term& second)
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
		const hol_quantifier& first,
		const hol_quantifier& second)
{
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
	return compare(*first.operand, *second.operand);
}

inline int_fast8_t compare(
		const hol_any& first,
		const hol_any& second)
{
	if (first.included == nullptr) {
		if (second.included != nullptr) return -1;
	} else if (second.included == nullptr) {
		return 1;
	} else {
		int_fast8_t result = compare(*first.included, *second.included);
		if (result != 0) return result;
	}

	if (first.excluded_tree_count < second.excluded_tree_count) return -1;
	else if (first.excluded_tree_count > second.excluded_tree_count) return 1;

	for (unsigned int i = 0; i < first.excluded_tree_count; i++) {
		int_fast8_t result = compare(*first.excluded_trees[i], *second.excluded_trees[i]);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const hol_any_array& first,
		const hol_any_array& second)
{
	if (first.oper < second.oper) return -1;
	else if (first.oper > second.oper) return 1;

	int_fast8_t result = compare(*first.all, *second.all);
	if (result != 0) return result;

	result = compare(first.any, second.any);
	if (result != 0) return result;

	result = compare(first.left, second.left);
	if (result != 0) return result;

	result = compare(first.right, second.right);
	if (result != 0) return result;

	return 0;
}

inline int_fast8_t compare(
		const hol_any_constant& first,
		const hol_any_constant& second)
{
	if (first.length < second.length) return -1;
	else if (first.length > second.length) return 1;

	for (unsigned int i = 0; i < first.length; i++) {
		if (first.constants[i] < second.constants[i]) return -1;
		else if (first.constants[i] > second.constants[i]) return 1;
	}
	return 0;
}

inline int_fast8_t compare(
		const hol_any_quantifier& first,
		const hol_any_quantifier& second)
{
	if (first.quantifier < second.quantifier) return -1;
	else if (first.quantifier > second.quantifier) return 1;
	else return compare(*first.operand, *second.operand);
}

inline int_fast8_t compare(
		const string& first,
		const string& second)
{
	if (first.length < second.length) return -1;
	else if (first.length > second.length) return 1;
	for (unsigned int i = 0; i < first.length; i++) {
		if (first.data[i] < second.data[i]) return -1;
		else if (first.data[i] > second.data[i]) return 1;
	}
	return 0;
}

inline int_fast8_t compare(
		const sequence& first,
		const sequence& second)
{
	if (first.length < second.length) return -1;
	else if (first.length > second.length) return 1;
	for (unsigned int i = 0; i < first.length; i++) {
		if (first.tokens[i] < second.tokens[i]) return -1;
		else if (first.tokens[i] > second.tokens[i]) return 1;
	}
	return 0;
}

int_fast8_t compare(
		const hol_term& first,
		const hol_term& second)
{
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case hol_term_type::VARIABLE:
		if (first.variable < second.variable) return -1;
		else if (first.variable > second.variable) return 1;
		else return 0;
	case hol_term_type::CONSTANT:
		if (first.constant < second.constant) return -1;
		else if (first.constant > second.constant) return 1;
		else return 0;
	case hol_term_type::PARAMETER:
		if (first.parameter < second.parameter) return -1;
		else if (first.parameter > second.parameter) return 1;
		else return 0;
	case hol_term_type::INTEGER:
		if (first.integer < second.integer) return -1;
		else if (first.integer > second.integer) return 1;
		else return 0;
	case hol_term_type::STRING:
		return compare(first.str, second.str);
	case hol_term_type::UINT_LIST:
		return compare(first.uint_list, second.uint_list);
	case hol_term_type::NOT:
		return compare(first.unary, second.unary);
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return compare(first.binary, second.binary);
	case hol_term_type::BINARY_APPLICATION:
		return compare(first.ternary, second.ternary);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return compare(first.array, second.array);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return compare(first.quantifier, second.quantifier);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		return compare(first.any, second.any);
	case hol_term_type::ANY_ARRAY:
		return compare(first.any_array, second.any_array);
	case hol_term_type::ANY_CONSTANT:
		return compare(first.any_constant, second.any_constant);
	case hol_term_type::ANY_QUANTIFIER:
		return compare(first.any_quantifier, second.any_quantifier);
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return 0;
	}
	fprintf(stderr, "compare ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline int_fast8_t compare(
		const hol_term* first,
		const hol_term* second)
{
	return compare(*first, *second);
}

inline bool operator < (
		const hol_term& first,
		const hol_term& second)
{
	return compare(first, second) < 0;
}


/* forward declarations */
struct hol_scope;
bool operator == (const hol_scope&, const hol_scope&);
int_fast8_t compare(const hol_scope&, const hol_scope&);
void shift_variables(hol_scope&, unsigned int);

template<bool AllConstantsDistinct>
bool canonicalize_scope(const hol_term&, hol_scope&, array_map<unsigned int, unsigned int>&, const equals_arg_types&);


template<unsigned int Arity>
struct hol_nary_scope {
	hol_scope* operands[Arity];

	~hol_nary_scope() {
		for (unsigned int i = 0; i < Arity; i++) {
			free(*operands[i]); free(operands[i]);
		}
	}

	static inline void move(const hol_nary_scope<Arity>& src, hol_nary_scope<Arity>& dst) {
		for (unsigned int i = 0; i < Arity; i++)
			dst.operands[i] = src.operands[i];
	}
};

struct hol_commutative_scope {
	array<hol_scope> children;
	array<hol_scope> negated;

	hol_commutative_scope() : children(4), negated(4) { }

	~hol_commutative_scope() {
		for (unsigned int i = 0; i < children.length; i++) free(children[i]);
		for (unsigned int i = 0; i < negated.length; i++) free(negated[i]);
	}

	static inline void move(const hol_commutative_scope& src, hol_commutative_scope& dst) {
		core::move(src.children, dst.children);
		core::move(src.negated, dst.negated);
	}
};

struct hol_noncommutative_scope {
	array<hol_scope> left, left_negated;
	array<hol_scope> right, right_negated;

	hol_noncommutative_scope() : left(4), left_negated(4), right(4), right_negated(4) { }

	~hol_noncommutative_scope() {
		for (unsigned int i = 0; i < left.length; i++) free(left[i]);
		for (unsigned int i = 0; i < left_negated.length; i++) free(left_negated[i]);
		for (unsigned int i = 0; i < right.length; i++) free(right[i]);
		for (unsigned int i = 0; i < right_negated.length; i++) free(right_negated[i]);
	}

	static inline void move(const hol_noncommutative_scope& src, hol_noncommutative_scope& dst) {
		core::move(src.left, dst.left);
		core::move(src.left_negated, dst.left_negated);
		core::move(src.right, dst.right);
		core::move(src.right_negated, dst.right_negated);
	}
};

struct hol_quantifier_scope {
	hol_scope* operand;
	unsigned int variable;

	~hol_quantifier_scope() {
		free(*operand); free(operand);
	}

	static inline void move(const hol_quantifier_scope& src, hol_quantifier_scope& dst) {
		dst.operand = src.operand;
		dst.variable = src.variable;
	}
};

struct hol_scope {
	hol_term_type type;
	array<unsigned int> variables;
	union {
		unsigned int variable;
		unsigned int constant;
		unsigned int parameter;
		int integer;
		string str;
		sequence uint_list;
		hol_scope* unary;
		hol_nary_scope<2> binary;
		hol_nary_scope<3> ternary;
		hol_commutative_scope commutative;
		hol_noncommutative_scope noncommutative;
		hol_quantifier_scope quantifier;
	};

	hol_scope(hol_term_type type) : variables(8) {
		if (!init_helper(type))
			exit(EXIT_FAILURE);
	}

	hol_scope(hol_term_type type, const array<unsigned int>& src_variables) : variables(src_variables.capacity) {
		if (!init_helper(type, src_variables))
			exit(EXIT_FAILURE);
	}

	~hol_scope() { free_helper(); }

	static inline void move(const hol_scope& src, hol_scope& dst)
	{
		dst.type = src.type;
		core::move(src.variables, dst.variables);
		switch (src.type) {
		case hol_term_type::CONSTANT:
			dst.constant = src.constant; return;
		case hol_term_type::VARIABLE:
			dst.variable = src.variable; return;
		case hol_term_type::PARAMETER:
			dst.parameter = src.parameter; return;
		case hol_term_type::INTEGER:
			dst.integer = src.integer; return;
		case hol_term_type::STRING:
			string::move(src.str, dst.str); return;
		case hol_term_type::UINT_LIST:
			sequence::move(src.uint_list, dst.uint_list); return;
		case hol_term_type::AND:
		case hol_term_type::OR:
		case hol_term_type::IFF:
			hol_commutative_scope::move(src.commutative, dst.commutative); return;
		case hol_term_type::IF_THEN:
			hol_noncommutative_scope::move(src.noncommutative, dst.noncommutative); return;
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			hol_quantifier_scope::move(src.quantifier, dst.quantifier); return;
		case hol_term_type::NOT:
			dst.unary = src.unary; return;
		case hol_term_type::UNARY_APPLICATION:
		case hol_term_type::EQUALS:
			hol_nary_scope<2>::move(src.binary, dst.binary); return;
		case hol_term_type::BINARY_APPLICATION:
			hol_nary_scope<3>::move(src.ternary, dst.ternary); return;
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
			return;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
		case hol_term_type::ANY_ARRAY:
		case hol_term_type::ANY_CONSTANT:
		case hol_term_type::ANY_QUANTIFIER:
			fprintf(stderr, "hol_scope.move ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
			exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
		}
		fprintf(stderr, "hol_scope.move ERROR: Unrecognized hol_term_type.\n");
		exit(EXIT_FAILURE);
	}

	static inline void swap(hol_scope& first, hol_scope& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(hol_scope); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void free(hol_scope& scope) {
		scope.free_helper();
		core::free(scope.variables);
	}

private:
	inline bool init_helper(hol_term_type scope_type) {
		type = scope_type;
		switch (type) {
		case hol_term_type::AND:
		case hol_term_type::OR:
		case hol_term_type::IFF:
			new (&commutative) hol_commutative_scope(); return true;
		case hol_term_type::IF_THEN:
			new (&noncommutative) hol_noncommutative_scope(); return true;
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			new (&quantifier) hol_quantifier_scope(); return true;
		case hol_term_type::UNARY_APPLICATION:
		case hol_term_type::EQUALS:
			new (&binary) hol_nary_scope<2>(); return true;
		case hol_term_type::BINARY_APPLICATION:
			new (&ternary) hol_nary_scope<3>(); return true;
		case hol_term_type::NOT:
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
		case hol_term_type::CONSTANT:
		case hol_term_type::VARIABLE:
		case hol_term_type::PARAMETER:
		case hol_term_type::INTEGER:
		case hol_term_type::STRING:
		case hol_term_type::UINT_LIST:
			return true;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
		case hol_term_type::ANY_ARRAY:
		case hol_term_type::ANY_CONSTANT:
		case hol_term_type::ANY_QUANTIFIER:
			fprintf(stderr, "hol_scope.init_helper ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
			exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
		}
		fprintf(stderr, "hol_scope.init_helper ERROR: Unrecognized hol_term_type.\n");
		return false;
	}

	inline bool init_helper(hol_term_type scope_type, const array<unsigned int>& src_variables) {
		for (unsigned int i = 0; i < src_variables.length; i++)
			variables[i] = src_variables[i];
		variables.length = src_variables.length;
		return init_helper(scope_type);
	}

	void free_helper() {
		switch (type) {
		case hol_term_type::AND:
		case hol_term_type::OR:
		case hol_term_type::IFF:
			commutative.~hol_commutative_scope(); return;
		case hol_term_type::IF_THEN:
			noncommutative.~hol_noncommutative_scope(); return;
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			quantifier.~hol_quantifier_scope(); return;
		case hol_term_type::NOT:
			core::free(*unary); core::free(unary); return;
		case hol_term_type::UNARY_APPLICATION:
		case hol_term_type::EQUALS:
			binary.~hol_nary_scope(); return;
		case hol_term_type::BINARY_APPLICATION:
			ternary.~hol_nary_scope(); return;
		case hol_term_type::STRING:
			core::free(str); return;
		case hol_term_type::UINT_LIST:
			core::free(uint_list); return;
		case hol_term_type::CONSTANT:
		case hol_term_type::VARIABLE:
		case hol_term_type::PARAMETER:
		case hol_term_type::INTEGER:
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
			return;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
		case hol_term_type::ANY_ARRAY:
		case hol_term_type::ANY_CONSTANT:
		case hol_term_type::ANY_QUANTIFIER:
			fprintf(stderr, "hol_scope.free_helper ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
			exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
		}
		fprintf(stderr, "hol_scope.free_helper ERROR: Unrecognized hol_term_type.\n");
		exit(EXIT_FAILURE);
	}

	friend bool init(hol_scope&, hol_term_type, unsigned int);
	friend bool init(hol_scope&, hol_term_type, const array<unsigned int>&);
};

inline bool init(hol_scope& scope, hol_term_type type, unsigned int variable_capacity = 8) {
	if (!array_init(scope.variables, variable_capacity)) {
		return false;
	} else if (!scope.init_helper(type)) {
		free(scope.variables);
		return false;
	}
	return true;
}

inline bool init(hol_scope& scope, hol_term_type type,
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

inline bool operator != (const hol_scope& first, const hol_scope& second) {
	return !(first == second);
}

template<unsigned int Arity>
inline bool operator == (const hol_nary_scope<Arity>& first, const hol_nary_scope<Arity>& second) {
	for (unsigned int i = 0; i < Arity; i++)
		if (*first.operands[i] != *second.operands[i]) return false;
	return true;
}

inline bool operator == (const hol_commutative_scope& first, const hol_commutative_scope& second)
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

inline bool operator == (const hol_noncommutative_scope& first, const hol_noncommutative_scope& second)
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

inline bool operator == (const hol_quantifier_scope& first, const hol_quantifier_scope& second) {
	return first.variable == second.variable
		&& *first.operand == *second.operand;
}

inline bool operator == (const hol_scope& first, const hol_scope& second)
{
	if (first.type != second.type)
		return false;
	switch (first.type) {
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return first.commutative == second.commutative;
	case hol_term_type::IF_THEN:
		return first.noncommutative == second.noncommutative;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return first.quantifier == second.quantifier;
	case hol_term_type::NOT:
		return *first.unary == *second.unary;
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::EQUALS:
		return first.binary == second.binary;
	case hol_term_type::BINARY_APPLICATION:
		return first.ternary == second.ternary;
	case hol_term_type::CONSTANT:
		return first.constant == second.constant;
	case hol_term_type::VARIABLE:
		return first.variable == second.variable;
	case hol_term_type::PARAMETER:
		return first.parameter == second.parameter;
	case hol_term_type::INTEGER:
		return first.integer == second.integer;
	case hol_term_type::STRING:
		return first.str == second.str;
	case hol_term_type::UINT_LIST:
		return first.uint_list == second.uint_list;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "operator == ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
		exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
	}
	fprintf(stderr, "operator == ERROR: Unrecognized hol_term_type when comparing hol_scopes.\n");
	exit(EXIT_FAILURE);
}

template<unsigned int Arity>
inline int_fast8_t compare(
		const hol_nary_scope<Arity>& first,
		const hol_nary_scope<Arity>& second)
{
	int_fast8_t result;
	for (unsigned int i = 0; i < Arity; i++) {
		result = compare(*first.operands[i], *second.operands[i]);
		if (result != 0) return result;
	}
	return 0;
}

inline int_fast8_t compare(
		const hol_commutative_scope& first,
		const hol_commutative_scope& second)
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
		const hol_noncommutative_scope& first,
		const hol_noncommutative_scope& second)
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
		const hol_quantifier_scope& first,
		const hol_quantifier_scope& second)
{
	if (first.variable < second.variable) return -1;
	else if (first.variable > second.variable) return 1;
	return compare(*first.operand, *second.operand);
}

int_fast8_t compare(
		const hol_scope& first,
		const hol_scope& second)
{
	if (first.type < second.type) return -1;
	else if (first.type > second.type) return 1;
	switch (first.type) {
	case hol_term_type::NOT:
		return compare(*first.unary, *second.unary);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return compare(first.commutative, second.commutative);
	case hol_term_type::IF_THEN:
		return compare(first.noncommutative, second.noncommutative);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return compare(first.quantifier, second.quantifier);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::EQUALS:
		return compare(first.binary, second.binary);
	case hol_term_type::BINARY_APPLICATION:
		return compare(first.ternary, second.ternary);
	case hol_term_type::CONSTANT:
		if (first.constant < second.constant) return -1;
		else if (first.constant > second.constant) return 1;
		else return 0;
	case hol_term_type::VARIABLE:
		if (first.variable < second.variable) return -1;
		else if (first.variable > second.variable) return 1;
		else return 0;
	case hol_term_type::PARAMETER:
		if (first.parameter < second.parameter) return -1;
		else if (first.parameter > second.parameter) return 1;
		else return 0;
	case hol_term_type::INTEGER:
		if (first.integer < second.integer) return -1;
		else if (first.integer > second.integer) return 1;
		else return 0;
	case hol_term_type::STRING:
		return compare(first.str, second.str);
	case hol_term_type::UINT_LIST:
		return compare(first.uint_list, second.uint_list);
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return 0;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "compare ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
		exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
	}
	fprintf(stderr, "compare ERROR: Unrecognized hol_term_type when comparing hol_scopes.\n");
	exit(EXIT_FAILURE);
}

struct hol_scope_canonicalizer { };

inline bool less_than(
		const hol_scope& first,
		const hol_scope& second,
		const hol_scope_canonicalizer& sorter)
{
	return compare(first, second) < 0;
}

template<unsigned int Arity>
inline void shift_variables(hol_nary_scope<Arity>& scope, unsigned int removed_variable) {
	for (unsigned int i = 0; i < Arity; i++)
		shift_variables(*scope.operands[i], removed_variable);
}

inline void shift_variables(hol_commutative_scope& scope, unsigned int removed_variable) {
	for (hol_scope& child : scope.children)
		shift_variables(child, removed_variable);
	for (hol_scope& child : scope.negated)
		shift_variables(child, removed_variable);
}

inline void shift_variables(hol_noncommutative_scope& scope, unsigned int removed_variable) {
	for (hol_scope& child : scope.left)
		shift_variables(child, removed_variable);
	for (hol_scope& child : scope.left_negated)
		shift_variables(child, removed_variable);
	for (hol_scope& child : scope.right)
		shift_variables(child, removed_variable);
	for (hol_scope& child : scope.right_negated)
		shift_variables(child, removed_variable);
}

void shift_variables(hol_scope& scope, unsigned int removed_variable) {
	for (unsigned int i = 0; i < scope.variables.length; i++) {
		if (scope.variables[i] > removed_variable)
			scope.variables[i]--;
	}
	switch (scope.type) {
	case hol_term_type::VARIABLE:
		if (scope.variable > removed_variable)
			scope.variable--;
		return;
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
		return;
	case hol_term_type::NOT:
		shift_variables(*scope.unary, removed_variable); return;
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		shift_variables(scope.commutative, removed_variable); return;
	case hol_term_type::IF_THEN:
		shift_variables(scope.noncommutative, removed_variable); return;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (scope.quantifier.variable > removed_variable)
			scope.quantifier.variable--;
		return shift_variables(*scope.quantifier.operand, removed_variable);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::EQUALS:
		shift_variables(scope.binary, removed_variable); return;
	case hol_term_type::BINARY_APPLICATION:
		shift_variables(scope.ternary, removed_variable); return;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "shift_variables ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
		exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
	}
	fprintf(stderr, "shift_variables ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

hol_term* scope_to_term(const hol_scope& scope);

template<bool Negated>
inline hol_term* scope_to_term(const hol_scope& scope);

template<>
inline hol_term* scope_to_term<false>(const hol_scope& scope) {
	return scope_to_term(scope);
}

template<>
inline hol_term* scope_to_term<true>(const hol_scope& negated)
{
	hol_term* negation;
	if (!new_hol_term(negation)) return NULL;
	negation->type = hol_term_type::NOT;
	negation->reference_count = 1;

	negation->unary.operand = scope_to_term(negated);
	if (negation->unary.operand == NULL) {
		free(negation); return NULL;
	}
	return negation;
}

template<hol_term_type ScopeType, bool Negated>
inline hol_term* scope_to_term(const hol_scope* scope,
		unsigned int scope_length, hol_term* first)
{
	static_assert(ScopeType == hol_term_type::AND
			   || ScopeType == hol_term_type::OR,
			"ScopeType must be either: AND or OR.");

	if (scope_length == 0) return first;

	hol_term* new_term;
	if (!new_hol_term(new_term)) {
		free(*first); if (first->reference_count == 0) free(first);
		return NULL;
	}
	new_term->type = ScopeType;
	new_term->reference_count = 1;
	new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * (scope_length + 1));
	new_term->array.length = scope_length + 1;
	if (new_term->array.operands == NULL) {
		free(*first); if (first->reference_count == 0) free(first);
		free(new_term); return NULL;
	}
	new_term->array.operands[0] = first;
	for (unsigned int i = 0; i < scope_length; i++) {
		new_term->array.operands[i + 1] = scope_to_term<Negated>(scope[i]);
		if (new_term->array.operands[i + 1] == NULL) {
			for (unsigned int j = 0; j < i + 1; j++) {
				free(*new_term->array.operands[j]);
				if (new_term->array.operands[j]->reference_count == 0)
					free(new_term->array.operands[j]);
			}
			free(new_term->array.operands); free(new_term);
			return NULL;
		}
	}
	return new_term;
}

template<hol_term_type ScopeType>
inline hol_term* scope_to_term(
		const hol_scope* scope, unsigned int scope_length,
		const hol_scope* negated, unsigned int negated_length)
{
	static_assert(ScopeType == hol_term_type::AND
			   || ScopeType == hol_term_type::OR,
			"ScopeType must be either: AND or OR.");

	if (scope_length == 1 && negated_length == 0)
		return scope_to_term<false>(scope[0]);
	else if (scope_length == 0 && negated_length == 1)
		return scope_to_term<true>(negated[0]);

	hol_term* new_term;
	if (!new_hol_term(new_term)) return NULL;
	new_term->type = ScopeType;
	new_term->reference_count = 1;
	new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * (scope_length + negated_length));
	new_term->array.length = scope_length + negated_length;
	if (new_term->array.operands == NULL) {
		free(new_term); return NULL;
	}
	for (unsigned int i = 0; i < scope_length; i++) {
		new_term->array.operands[i] = scope_to_term<false>(scope[i]);
		if (new_term->array.operands[i] == NULL) {
			for (unsigned int j = 0; j < i; j++) {
				free(*new_term->array.operands[j]);
				if (new_term->array.operands[j]->reference_count == 0)
					free(new_term->array.operands[j]);
			}
			free(new_term->array.operands); free(new_term);
			return NULL;
		}
	} for (unsigned int i = 0; i < negated_length; i++) {
		new_term->array.operands[scope_length + i] = scope_to_term<true>(negated[i]);
		if (new_term->array.operands[scope_length + i] == NULL) {
			for (unsigned int j = 0; j < scope_length + i; j++) {
				free(*new_term->array.operands[j]);
				if (new_term->array.operands[j]->reference_count == 0)
					free(new_term->array.operands[j]);
			}
			free(new_term->array.operands); free(new_term);
			return NULL;
		}
	}
	return new_term;
}

template<bool Negated>
inline hol_term* iff_scope_to_term(const hol_scope* scope,
		unsigned int scope_length, hol_term* first)
{
	if (scope_length == 0) return first;

	for (unsigned int i = scope_length; i > 0; i--) {
		hol_term* next = scope_to_term<Negated>(scope[i - 1]);
		if (next == NULL) {
			free(*first); if (first->reference_count == 0) free(first);
			return NULL;
		}

		hol_term* new_right;
		if (!new_hol_term(new_right)) {
			free(*first); if (first->reference_count == 0) free(first);
			free(*next); if (next->reference_count == 0) free(next);
			return NULL;
		}
		new_right->type = hol_term_type::EQUALS;
		new_right->reference_count = 1;
		new_right->binary.left = next;
		new_right->binary.right = first;
		first = new_right;
	}
	return first;
}

inline hol_term* iff_scope_to_term(
		const hol_scope* scope, unsigned int scope_length,
		const hol_scope* negated, unsigned int negated_length)
{
	if (scope_length == 1 && negated_length == 0)
		return scope_to_term<false>(scope[0]);
	else if (scope_length == 0 && negated_length == 1)
		return scope_to_term<true>(negated[0]);

	hol_term* right;
	if (negated_length > 0) {
		right = scope_to_term<true>(negated[negated_length - 1]);
		negated_length--;
	} else {
		right = scope_to_term<false>(scope[scope_length - 1]);
		scope_length--;
	}
	if (right == NULL) return NULL;

	right = iff_scope_to_term<true>(negated, negated_length, right);
	if (right == NULL) return NULL;
	return iff_scope_to_term<false>(scope, scope_length, right);
}

template<hol_term_type QuantifierType>
inline hol_term* scope_to_term(const hol_quantifier_scope& scope)
{
	static_assert(QuantifierType == hol_term_type::FOR_ALL
			   || QuantifierType == hol_term_type::EXISTS
			   || QuantifierType == hol_term_type::LAMBDA,
			"QuantifierType is not a quantifier.");

	hol_term* operand = scope_to_term(*scope.operand);
	if (operand == NULL) return NULL;

	hol_term* new_term;
	if (!new_hol_term(new_term)) return NULL;
	new_term->type = QuantifierType;
	new_term->reference_count = 1;
	new_term->quantifier.variable = scope.variable;
	new_term->quantifier.operand = operand;
	return new_term;
}

template<hol_term_type Type>
hol_term* scope_to_term(const hol_nary_scope<2>& scope) {
	static_assert(Type == hol_term_type::EQUALS
			   || Type == hol_term_type::UNARY_APPLICATION,
			"Type must be either: EQUALS or UNARY_APPLICATION.\n");

	hol_term* left = scope_to_term(*scope.operands[0]);
	if (left == NULL) return NULL;
	hol_term* right = scope_to_term(*scope.operands[1]);
	if (right == NULL) {
		free(*left); if (left->reference_count == 0) free(left);
		return NULL;
	}

	hol_term* new_term;
	if (!new_hol_term(new_term)) return NULL;
	new_term->type = Type;
	new_term->reference_count = 1;
	new_term->binary.left = left;
	new_term->binary.right = right;
	return new_term;
}

template<hol_term_type Type>
hol_term* scope_to_term(const hol_nary_scope<3>& scope) {
	static_assert(Type == hol_term_type::BINARY_APPLICATION,
			"Type must be BINARY_APPLICATION.\n");

	hol_term* first = scope_to_term(*scope.operands[0]);
	if (first == NULL) return NULL;
	hol_term* second = scope_to_term(*scope.operands[1]);
	if (second == NULL) {
		free(*first); if (first->reference_count == 0) free(first);
		return NULL;
	}
	hol_term* third = scope_to_term(*scope.operands[2]);
	if (third == NULL) {
		free(*first); if (first->reference_count == 0) free(first);
		free(*second); if (second->reference_count == 0) free(second);
		return NULL;
	}

	hol_term* new_term;
	if (!new_hol_term(new_term)) return NULL;
	new_term->type = Type;
	new_term->reference_count = 1;
	new_term->ternary.first = first;
	new_term->ternary.second = second;
	new_term->ternary.third = third;
	return new_term;
}

inline hol_term* scope_to_term(const hol_scope& scope)
{
	hol_term* new_term;
	hol_term* first; hol_term* second;

	switch (scope.type) {
	case hol_term_type::AND:
		return scope_to_term<hol_term_type::AND>(
				scope.commutative.children.data, scope.commutative.children.length,
				scope.commutative.negated.data, scope.commutative.negated.length);
	case hol_term_type::OR:
		return scope_to_term<hol_term_type::OR>(
				scope.commutative.children.data, scope.commutative.children.length,
				scope.commutative.negated.data, scope.commutative.negated.length);
	case hol_term_type::IFF:
		if (scope.commutative.children.length > 0 && scope.commutative.children.last().type == hol_term_type::FALSE) {
			first = iff_scope_to_term(
					scope.commutative.children.data, scope.commutative.children.length - 1,
					scope.commutative.negated.data, scope.commutative.negated.length);
			if (!new_hol_term(new_term)) {
				free(*first); if (first->reference_count == 0) free(first);
				return NULL;
			}
			new_term->type = hol_term_type::NOT;
			new_term->reference_count = 1;
			new_term->unary.operand = first;
			return new_term;
		} else {
			return iff_scope_to_term(
					scope.commutative.children.data, scope.commutative.children.length,
					scope.commutative.negated.data, scope.commutative.negated.length);
		}
	case hol_term_type::IF_THEN:
		first = scope_to_term<hol_term_type::AND>(
				scope.noncommutative.left.data, scope.noncommutative.left.length,
				scope.noncommutative.left_negated.data, scope.noncommutative.left_negated.length);
		if (first == NULL) return NULL;
		second = scope_to_term<hol_term_type::OR>(
				scope.noncommutative.right.data, scope.noncommutative.right.length,
				scope.noncommutative.right_negated.data, scope.noncommutative.right_negated.length);
		if (second == NULL) {
			free(*first); if (first->reference_count == 0) free(first);
			return NULL;
		}

		if (!new_hol_term(new_term)) {
			free(*first); if (first->reference_count == 0) free(first);
			free(*second); if (second->reference_count == 0) free(second);
			return NULL;
		}
		new_term->type = hol_term_type::IF_THEN;
		new_term->reference_count = 1;
		new_term->binary.left = first;
		new_term->binary.right = second;
		return new_term;
	case hol_term_type::NOT:
		first = scope_to_term(*scope.unary);
		if (first == NULL) return NULL;

		if (!new_hol_term(new_term)) return NULL;
		new_term->type = hol_term_type::NOT;
		new_term->reference_count = 1;
		new_term->unary.operand = first;
		return new_term;
	case hol_term_type::FOR_ALL:
		return scope_to_term<hol_term_type::FOR_ALL>(scope.quantifier);
	case hol_term_type::EXISTS:
		return scope_to_term<hol_term_type::EXISTS>(scope.quantifier);
	case hol_term_type::LAMBDA:
		return scope_to_term<hol_term_type::LAMBDA>(scope.quantifier);
	case hol_term_type::EQUALS:
		return scope_to_term<hol_term_type::EQUALS>(scope.binary);
	case hol_term_type::UNARY_APPLICATION:
		return scope_to_term<hol_term_type::UNARY_APPLICATION>(scope.binary);
	case hol_term_type::BINARY_APPLICATION:
		return scope_to_term<hol_term_type::BINARY_APPLICATION>(scope.ternary);
	case hol_term_type::CONSTANT:
		return hol_term::new_constant(scope.constant);
	case hol_term_type::VARIABLE:
		return hol_term::new_variable(scope.variable);
	case hol_term_type::PARAMETER:
		return hol_term::new_parameter(scope.parameter);
	case hol_term_type::INTEGER:
		return hol_term::new_int(scope.integer);
	case hol_term_type::STRING:
		return hol_term::new_string(scope.str);
	case hol_term_type::UINT_LIST:
		return hol_term::new_uint_list(scope.uint_list);
	case hol_term_type::TRUE:
		HOL_TRUE.reference_count++;
		return &HOL_TRUE;
	case hol_term_type::FALSE:
		HOL_FALSE.reference_count++;
		return &HOL_FALSE;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "scope_to_term ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT,` `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
		exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT,` `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
	}
	fprintf(stderr, "scope_to_term ERROR: Unrecognized hol_term_type.\n");
	return nullptr;
}

inline void move_variables(
		const array<unsigned int>& term_variables,
		array<unsigned int>& scope_variables)
{
	array<unsigned int> variable_union(max(scope_variables.capacity, term_variables.length + scope_variables.length));
	set_union(variable_union.data, variable_union.length,
			term_variables.data, term_variables.length,
			scope_variables.data, scope_variables.length);
	swap(variable_union, scope_variables);
}

inline void recompute_variables(
		const array<hol_scope>& children,
		const array<hol_scope>& negated,
		array<unsigned int>& variables)
{
	variables.clear();
	for (const hol_scope& child : children)
		move_variables(child.variables, variables);
	for (const hol_scope& child : negated)
		move_variables(child.variables, variables);
}

inline void recompute_variables(
		const array<hol_scope>& left, const array<hol_scope>& left_negated,
		const array<hol_scope>& right, const array<hol_scope>& right_negated,
		array<unsigned int>& variables)
{
	variables.clear();
	for (const hol_scope& child : left)
		move_variables(child.variables, variables);
	for (const hol_scope& child : left_negated)
		move_variables(child.variables, variables);
	for (const hol_scope& child : right)
		move_variables(child.variables, variables);
	for (const hol_scope& child : right_negated)
		move_variables(child.variables, variables);
}

inline bool scope_contains(const hol_scope& subscope, const array<hol_scope>& scope, unsigned int& index)
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

inline bool scope_contains(const hol_scope& subscope, const array<hol_scope>& scope) {
	unsigned int index;
	return scope_contains(subscope, scope, index);
}

template<hol_term_type Operator>
bool add_to_scope_helper(
		hol_scope& subscope,
		array<hol_scope>& children,
		array<hol_scope>& negated,
		bool& found_negation)
{
	unsigned int i;
	if (scope_contains(subscope, negated, i)) {
		found_negation = true;
		if (Operator == hol_term_type::IFF) {
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
		if (Operator == hol_term_type::IFF) {
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

template<hol_term_type Operator>
bool add_to_scope(
		hol_scope& subscope,
		array<hol_scope>& children,
		array<hol_scope>& negated,
		bool& found_negation)
{
	/* check if `subscope` is the negation of any operand in `scope.commutative.children` */
	if (subscope.type == hol_term_type::NOT) {
		if (!add_to_scope_helper<Operator>(*subscope.unary, negated, children, found_negation)) return false;
		free(subscope.variables);
		free(subscope.unary);
		return true;
	} else if (subscope.type == hol_term_type::IFF && subscope.commutative.children.length > 0
			&& subscope.commutative.children.last().type == hol_term_type::FALSE)
	{
		free(subscope.commutative.children.last());
		subscope.commutative.children.length--;
		return add_to_scope_helper<Operator>(subscope, negated, children, found_negation);
	} else {
		return add_to_scope_helper<Operator>(subscope, children, negated, found_negation);
	}
}

template<hol_term_type Operator>
bool add_to_scope(
		hol_scope& subscope,
		array<hol_scope>& children,
		array<hol_scope>& negated,
		array<unsigned int>& variables,
		bool& found_negation)
{
	/* check if `subscope` is the negation of any operand in `scope.commutative.children` */
	if (subscope.type == hol_term_type::NOT) {
		if (Operator != hol_term_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(*subscope.unary, negated, children, found_negation)) return false;
		if (Operator == hol_term_type::IFF)
			recompute_variables(children, negated, variables);
		free(subscope.variables);
		free(subscope.unary);
	} else if (subscope.type == hol_term_type::IFF && subscope.commutative.children.length > 0
			&& subscope.commutative.children.last().type == hol_term_type::FALSE)
	{
		free(subscope.commutative.children.last());
		subscope.commutative.children.length--;
		if (Operator != hol_term_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(subscope, negated, children, found_negation)) return false;
		if (Operator == hol_term_type::IFF)
			recompute_variables(children, negated, variables);
	} else {
		if (Operator != hol_term_type::IFF)
			move_variables(subscope.variables, variables);
		if (!add_to_scope_helper<Operator>(subscope, children, negated, found_negation)) return false;
		if (Operator == hol_term_type::IFF)
			recompute_variables(children, negated, variables);
	}
	return true;
}

template<hol_term_type Operator>
unsigned int intersection_size(
	array<hol_scope>& first,
	array<hol_scope>& second)
{
	static_assert(Operator == hol_term_type::AND
			   || Operator == hol_term_type::OR
			   || Operator == hol_term_type::IF_THEN
			   || Operator == hol_term_type::IFF,
			"Operator must be either: AND, OR, IF_THEN, IFF.");

	unsigned int intersection_count = 0;
	unsigned int i = 0, j = 0, first_index = 0, second_index = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == hol_term_type::AND || Operator == hol_term_type::OR || Operator == hol_term_type::IF_THEN) {
				return 1;
			} else if (Operator == hol_term_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			if (Operator == hol_term_type::IFF) move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			if (Operator == hol_term_type::IFF) move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	if (Operator != hol_term_type::IFF) return intersection_count;

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

template<hol_term_type Operator>
unsigned int intersection_size(
	array<hol_scope>& first,
	array<hol_scope>& second,
	unsigned int skip_second_index)
{
	static_assert(Operator == hol_term_type::AND
			   || Operator == hol_term_type::OR
			   || Operator == hol_term_type::IF_THEN
			   || Operator == hol_term_type::IFF,
			"Operator must be either: AND, OR, IF_THEN, IFF.");

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
			if (Operator == hol_term_type::AND || Operator == hol_term_type::OR || Operator == hol_term_type::IF_THEN) {
				return 1;
			} else if (Operator == hol_term_type::IFF) {
				free(first[i]); free(second[j]);
				i++; j++; intersection_count++;
			}
		} else if (result < 0) {
			if (Operator == hol_term_type::IFF) move(first[i], first[first_index]);
			i++; first_index++;
		} else {
			move(second[j], second[second_index]);
			j++; second_index++;
		}
	}

	while (i < first.length) {
		if (Operator == hol_term_type::IFF) move(first[i], first[first_index]);
		i++; first_index++;
	} while (j < second.length) {
		move(second[j], second[second_index]);
		j++; second_index++;
	}
	first.length = first_index;
	second.length = second_index;
	return intersection_count;
}

template<hol_term_type Operator>
void merge_scopes(array<hol_scope>& dst,
	array<hol_scope>& first,
	array<hol_scope>& second)
{
	static_assert(Operator == hol_term_type::AND
			   || Operator == hol_term_type::OR
			   || Operator == hol_term_type::IFF,
			"Operator must be either: AND, OR, IFF.");

	unsigned int i = 0, j = 0;
	while (i < first.length && j < second.length)
	{
		auto result = compare(first[i], second[j]);
		if (result == 0) {
			if (Operator == hol_term_type::IFF) {
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

template<hol_term_type Operator>
void merge_scopes(array<hol_scope>& dst,
	array<hol_scope>& first,
	array<hol_scope>& second,
	unsigned int skip_second_index,
	unsigned int& new_second_index)
{
	static_assert(Operator == hol_term_type::AND
			   || Operator == hol_term_type::OR
			   || Operator == hol_term_type::IFF,
			"Operator must be either: AND, OR, IFF.");

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
			if (Operator == hol_term_type::IFF) {
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

template<hol_term_type Operator>
inline void merge_scopes(
		array<hol_scope>& src, array<hol_scope>& dst,
		array<hol_scope>& src_negated, array<hol_scope>& dst_negated,
		bool& found_negation)
{
	unsigned int intersection_count = intersection_size<Operator>(src, dst_negated);
	if (intersection_count > 0 && (Operator == hol_term_type::AND || Operator == hol_term_type::OR)) {
		for (hol_scope& child : src) free(child);
		for (hol_scope& child : src_negated) free(child);
		found_negation = true; return;
	}

	intersection_count += intersection_size<Operator>(src_negated, dst);
	if (intersection_count > 0 && (Operator == hol_term_type::AND || Operator == hol_term_type::OR)) {
		for (hol_scope& child : src) free(child);
		for (hol_scope& child : src_negated) free(child);
		found_negation = true; return;
	} else if (Operator == hol_term_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	if (Operator == hol_term_type::IFF
	 && src.length == 0 && dst.length == 0
	 && src_negated.length == 0 && dst_negated.length == 0)
		return; /* this happens if the elements of `src` are all negations of the elements of `dst` (and vice versa) */

	/* merge the two scopes */
	array<hol_scope> both = array<hol_scope>(max((size_t) 1, src.length + dst.length));
	merge_scopes<Operator>(both, src, dst);
	swap(both, dst); both.clear();

	array<hol_scope> both_negated = array<hol_scope>(max((size_t) 1, src_negated.length + dst_negated.length));
	merge_scopes<Operator>(both_negated, src_negated, dst_negated);
	swap(both_negated, dst_negated);
}

template<hol_term_type Operator>
inline void merge_scopes(
		array<hol_scope>& src, array<hol_scope>& dst,
		array<hol_scope>& src_negated, array<hol_scope>& dst_negated,
		bool& found_negation, unsigned int skip_dst_index, unsigned int& new_dst_index)
{
	unsigned int intersection_count = intersection_size<Operator>(src, dst_negated);
	if (intersection_count > 0 && (Operator == hol_term_type::AND || Operator == hol_term_type::OR)) {
		for (hol_scope& child : src) free(child);
		for (hol_scope& child : src_negated) free(child);
		found_negation = true; return;
	}

	intersection_count += intersection_size<Operator>(src_negated, dst, skip_dst_index);
	if (intersection_count > 0 && (Operator == hol_term_type::AND || Operator == hol_term_type::OR)) {
		for (hol_scope& child : src) free(child);
		for (hol_scope& child : src_negated) free(child);
		found_negation = true; return;
	} else if (Operator == hol_term_type::IFF && intersection_count % 2 == 1) {
		found_negation = true;
	} else {
		found_negation = false;
	}

	if (Operator == hol_term_type::IFF
	 && src.length == 0 && dst.length == 0
	 && src_negated.length == 0 && dst_negated.length == 0)
		return; /* this happens if the elements of `src` are all negations of the elements of `dst` (and vice versa) */

	/* merge the two scopes */
	array<hol_scope> both = array<hol_scope>(max((size_t) 1, src.length + dst.length));
	merge_scopes<Operator>(both, src, dst, skip_dst_index, new_dst_index);
	swap(both, dst); both.clear();

	array<hol_scope> both_negated = array<hol_scope>(max((size_t) 1, src_negated.length + dst_negated.length));
	merge_scopes<Operator>(both_negated, src_negated, dst_negated);
	swap(both_negated, dst_negated);
}

inline bool negate_iff(hol_scope& scope)
{
	if (!scope.commutative.children.ensure_capacity(scope.commutative.children.length + 1)) {
		return false;
	} else if (scope.commutative.children.length > 0 && scope.commutative.children.last().type == hol_term_type::FALSE) {
		free(scope.commutative.children.last());
		scope.commutative.children.length--;
	} else {
		if (!init(scope.commutative.children[scope.commutative.children.length], hol_term_type::FALSE))
			return false;
		scope.commutative.children.length++;
	}
	return true;
}

bool negate_scope(hol_scope& scope)
{
	if (scope.type == hol_term_type::TRUE) {
		scope.type = hol_term_type::FALSE;
	} else if (scope.type == hol_term_type::FALSE) {
		scope.type = hol_term_type::TRUE;
	} else if (scope.type == hol_term_type::NOT) {
		hol_scope* operand = scope.unary;
		free(scope.variables);
		move(*operand, scope);
		free(operand);
	} else if (scope.type == hol_term_type::IFF) {
		return negate_iff(scope);
	} else {
		hol_scope* operand = (hol_scope*) malloc(sizeof(hol_scope));
		if (operand == NULL) {
			fprintf(stderr, "negate_scope ERROR: Out of memory.\n");
			return false;
		}

		move(scope, *operand);
		if (!init(scope, hol_term_type::NOT, operand->variables)) {
			move(*operand, scope); free(operand);
			return false;
		}
		scope.unary = operand;
	}
	return true;
}

inline bool are_negations(hol_scope& left, hol_scope& right)
{
	if ((left.type == hol_term_type::NOT && *left.unary == right)
	 || (right.type == hol_term_type::NOT && *right.unary == left))
		return true;
	if (left.type == hol_term_type::IFF && right.type == hol_term_type::IFF) {
		if (left.commutative.children.length > 0 && left.commutative.children.last().type == hol_term_type::FALSE) {
			/* `left` is a negated IFF expression */
			if (right.commutative.children.length > 0 && right.commutative.children.last().type == hol_term_type::FALSE) {
				/* `right` is a negated IFF expression */
				return false;
			} else {
				left.commutative.children.length--;
				bool negated = (left == right);
				left.commutative.children.length++;
				return negated;
			}
		} else {
			if (right.commutative.children.length > 0 && right.commutative.children.last().type == hol_term_type::FALSE) {
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

template<hol_term_type Operator, bool AllConstantsDistinct>
bool canonicalize_commutative_scope(
		const hol_array_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	static_assert(Operator == hol_term_type::AND
			   || Operator == hol_term_type::OR
			   || Operator == hol_term_type::IFF,
			"Operator must be either: AND, OR, IFF.");

	if (!init(out, Operator)) return false;

	hol_scope& next = *((hol_scope*) alloca(sizeof(hol_scope)));
	for (unsigned int i = 0; i < src.length; i++) {
		if (!canonicalize_scope<AllConstantsDistinct>(*src.operands[i], next, variable_map, types)) {
			free(out); return false;
		}

		if (next.type == hol_term_type::FALSE) {
			free(next);
			if (Operator == hol_term_type::AND) {
				free(out);
				return init(out, hol_term_type::FALSE);
			} else if (Operator == hol_term_type::IFF) {
				if (!negate_iff(out)) return false;
			}
		} else if (next.type == hol_term_type::TRUE) {
			free(next);
			if (Operator == hol_term_type::OR) {
				free(out);
				return init(out, hol_term_type::TRUE);
			}
		} else if (next.type == Operator) {
			bool found_negation;
			merge_scopes<Operator>(next.commutative.children, out.commutative.children,
					next.commutative.negated, out.commutative.negated, found_negation);
			free(next.commutative.children);
			free(next.commutative.negated);
			if (Operator == hol_term_type::IFF)
				recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
			else move_variables(next.variables, out.variables);
			free(next.variables);
			if (found_negation) {
				if (Operator == hol_term_type::AND) {
					free(out);
					return init(out, hol_term_type::FALSE);
				} else if (Operator == hol_term_type::OR) {
					free(out);
					return init(out, hol_term_type::TRUE);
				} else if (Operator == hol_term_type::IFF) {
					recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
					if (!negate_iff(out)) return false;
				}
			}
		} else {
			bool found_negation;
			if (!add_to_scope<Operator>(next, out.commutative.children, out.commutative.negated, out.variables, found_negation)) {
				free(out); free(next); return false;
			} else if (found_negation) {
				if (Operator == hol_term_type::AND) {
					free(out); return init(out, hol_term_type::FALSE);
				} else if (Operator == hol_term_type::OR) {
					free(out); return init(out, hol_term_type::TRUE);
				} else if (Operator == hol_term_type::IFF) {
					if (!negate_iff(out)) return false;
				}
			}
		}
	}

	if (out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
		free(out);
		if (Operator == hol_term_type::AND || Operator == hol_term_type::IFF)
			return init(out, hol_term_type::TRUE);
		else return init(out, hol_term_type::FALSE);
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
		const hol_binary_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	hol_scope& left = *((hol_scope*) alloca(sizeof(hol_scope)));
	if (!canonicalize_scope<AllConstantsDistinct>(*src.left, left, variable_map, types)) return false;

	if (left.type == hol_term_type::FALSE) {
		free(left);
		return init(out, hol_term_type::TRUE);
	} else if (left.type == hol_term_type::TRUE) {
		free(left);
		return canonicalize_scope<AllConstantsDistinct>(*src.right, out, variable_map, types);
	}

	if (!canonicalize_scope<AllConstantsDistinct>(*src.right, out, variable_map, types))
		return false;

	if (out == left) {
		/* consider the case where `*src.left` and `*src.right` are identical */
		free(out); free(left);
		return init(out, hol_term_type::TRUE);
	} else if (out.type == hol_term_type::FALSE) {
		free(out); move(left, out);
		if (!negate_scope(out)) {
			free(out); return false;
		}
	} else if (out.type == hol_term_type::TRUE) {
		free(left); /* this is a no-op */
	} else if (are_negations(out, left)) {
		free(left); /* we have the case `A => ~A` or `~A => A`, which is also a no-op */
	} else {
		/* first construct the conditional */
		if (out.type == hol_term_type::OR) {
			hol_scope temp = hol_scope(hol_term_type::IF_THEN);
			swap(out.commutative.children, temp.noncommutative.right);
			swap(out.commutative.negated, temp.noncommutative.right_negated);
			swap(out.variables, temp.variables);
			swap(out, temp);

			/* check if any operands in the OR can be raised into this IF_THEN consequent scope */
			array<hol_scope> to_merge = array<hol_scope>(8);
			for (unsigned int i = 0; i < out.noncommutative.right.length; i++) {
				hol_scope& child = out.noncommutative.right[i];
				if (child.type != hol_term_type::IF_THEN) continue;

				bool found_negation;
				hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
				move(child, temp);
				merge_scopes<hol_term_type::AND>(temp.noncommutative.left, out.noncommutative.left,
						temp.noncommutative.left_negated, out.noncommutative.left_negated, found_negation);
				temp.noncommutative.left.clear();
				temp.noncommutative.left_negated.clear();
				if (found_negation) {
					child.noncommutative.left.clear();
					child.noncommutative.left_negated.clear();
					free(out); free(left);
					return init(out, hol_term_type::TRUE);
				}

				unsigned int new_index;
				merge_scopes<hol_term_type::OR>(temp.noncommutative.right, out.noncommutative.right,
						temp.noncommutative.right_negated, out.noncommutative.right_negated, found_negation, i, new_index);
				temp.noncommutative.right.clear();
				temp.noncommutative.right_negated.clear();
				if (found_negation) {
					child.noncommutative.left.clear();
					child.noncommutative.left_negated.clear();
					child.noncommutative.right.clear();
					child.noncommutative.right_negated.clear();
					free(out); free(left);
					return init(out, hol_term_type::TRUE);
				} else {
					free(temp);
				}
				i = new_index;
			}
		} else if (out.type == hol_term_type::NOT) {
			hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
			if (!init(temp, hol_term_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			move(*out.unary, temp.noncommutative.right_negated[0]);
			temp.noncommutative.right_negated.length++;
			free(out.variables); free(out.unary);
			move(temp, out);
		} else if (out.type == hol_term_type::IFF && out.commutative.children.length > 0
				&& out.commutative.children.last().type == hol_term_type::FALSE)
		{
			hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
			if (!init(temp, hol_term_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			free(out.commutative.children.last());
			out.commutative.children.length--;
			move(out, temp.noncommutative.right_negated[0]);
			temp.noncommutative.right_negated.length++;
			move(temp, out);
		} else if (out.type != hol_term_type::IF_THEN) {
			hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
			if (!init(temp, hol_term_type::IF_THEN, out.variables)) {
				free(out); free(left);
				return false;
			}
			move(out, temp.noncommutative.right[0]);
			temp.noncommutative.right.length++;
			move(temp, out);
		}

		/* now try merging `left` with the conditional */
		if (left.type == hol_term_type::AND) {
			bool found_negation;
			merge_scopes<hol_term_type::AND>(left.commutative.children, out.noncommutative.left,
					left.commutative.negated, out.noncommutative.left_negated, found_negation);
			left.commutative.children.clear();
			left.commutative.negated.clear();
			if (found_negation) {
				free(out); free(left);
				return init(out, hol_term_type::TRUE);
			}
			move_variables(left.variables, out.variables);
			free(left);
		} else {
			bool found_negation;
			if (!add_to_scope<hol_term_type::AND>(left, out.noncommutative.left, out.noncommutative.left_negated, out.variables, found_negation)) {
				free(out); free(left); return false;
			} else if (found_negation) {
				free(out);
				return init(out, hol_term_type::TRUE);
			}
		}

		/* check if antecedent and consequent have any common operands */
		if (intersection_size<hol_term_type::IF_THEN>(out.noncommutative.left, out.noncommutative.right) > 0
		 || intersection_size<hol_term_type::IF_THEN>(out.noncommutative.left_negated, out.noncommutative.right_negated) > 0)
		{
			free(out);
			return init(out, hol_term_type::TRUE);
		}
	}
	return true;
}

inline bool promote_from_quantifier_scope(
		array<hol_scope>& quantifier_operand,
		array<hol_scope>& dst,
		unsigned int quantifier_variable)
{
	unsigned int next = 0;
	for (unsigned int i = 0; i < quantifier_operand.length; i++) {
		hol_scope& child = quantifier_operand[i];
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

template<hol_term_type QuantifierType>
inline bool make_quantifier_scope(hol_scope& out, hol_scope* operand, unsigned int quantifier_variable)
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
template<hol_term_type QuantifierType>
bool process_conditional_quantifier_scope(hol_scope&, hol_scope*, unsigned int);


template<hol_term_type QuantifierType>
bool process_commutative_quantifier_scope(
		hol_scope& out, hol_scope* operand,
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
		hol_scope* quantifier_operand = (hol_scope*) malloc(sizeof(hol_scope));
		if (quantifier_operand == NULL) {
			fprintf(stderr, "process_commutative_quantifier_scope ERROR: Out of memory.\n");
			free(out); free(*operand); free(operand); return false;
		}

		if (operand->commutative.children.length == 1 && operand->commutative.negated.length == 0) {
			hol_scope& inner_operand = operand->commutative.children[0];
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
		hol_scope& quantifier = *((hol_scope*) alloca(sizeof(hol_scope)));
		if (quantifier_operand->type != out.type && (quantifier_operand->type == hol_term_type::AND || quantifier_operand->type == hol_term_type::OR)) {
			if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
				return false;
		} else if (quantifier_operand->type == hol_term_type::IF_THEN) {
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
			if (out.type == hol_term_type::AND) {
				if (!add_to_scope<hol_term_type::AND>(quantifier, out.commutative.children, out.commutative.negated, found_negation)) {
					free(out); free(quantifier); return false;
				} else if (found_negation) {
					free(out);
					return init(out, hol_term_type::FALSE);
				}
			} else if (out.type == hol_term_type::OR) {
				if (!add_to_scope<hol_term_type::OR>(quantifier, out.commutative.children, out.commutative.negated, found_negation)) {
					free(out); free(quantifier);
				} else if (found_negation) {
					free(out);
					return init(out, hol_term_type::TRUE);
				}
			}
			recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
		}
	}
	return true;
}

template<hol_term_type QuantifierType>
bool process_conditional_quantifier_scope(
		hol_scope& out, hol_scope* operand,
		unsigned int quantifier_variable)
{
	if (!init(out, hol_term_type::IF_THEN)) {
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
			hol_scope* quantifier_operand = (hol_scope*) malloc(sizeof(hol_scope));
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
				if (!init(*quantifier_operand, hol_term_type::OR)) {
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
			hol_scope& quantifier = *((hol_scope*) alloca(sizeof(hol_scope)));
			if (quantifier_operand->type == hol_term_type::AND || quantifier_operand->type == hol_term_type::OR) {
				if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (!make_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<hol_term_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, hol_term_type::TRUE);
			}
			recompute_variables(out.noncommutative.left, out.noncommutative.left_negated,
					out.noncommutative.right, out.noncommutative.right_negated, out.variables);
		}
	} else {
		/* the antecedent is non-empty */
		if (operand->noncommutative.right.length == 0 && operand->noncommutative.right_negated.length == 0) {
			/* we've moved all children out of the consequent */
			hol_scope* quantifier_operand = (hol_scope*) malloc(sizeof(hol_scope));
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
				hol_scope* conjunction = (hol_scope*) malloc(sizeof(hol_scope));
				if (conjunction == NULL) {
					fprintf(stderr, "process_conditional_quantifier_scope ERROR: Out of memory.\n");
					free(quantifier_operand); free(out); free(*operand); free(operand);
				}
				conjunction->type = hol_term_type::AND;
				move(operand->noncommutative.left, conjunction->commutative.children);
				move(operand->noncommutative.left_negated, conjunction->commutative.negated);
				move(operand->variables, conjunction->variables);
				recompute_variables(conjunction->commutative.children, conjunction->commutative.negated, conjunction->variables);
				free(operand->noncommutative.right); free(operand->noncommutative.right_negated);
				free(operand);

				if (!init(*quantifier_operand, hol_term_type::NOT)) {
					free(quantifier_operand); free(out);
					return false;
				}
				quantifier_operand->unary = conjunction;
				move_variables(conjunction->variables, quantifier_operand->variables);
			}


			/* removing an operand can allow movement from its children to the
			   parent, so check if the new operand allows movement */
			hol_scope& quantifier = *((hol_scope*) alloca(sizeof(hol_scope)));
			if (quantifier_operand->type == hol_term_type::AND || quantifier_operand->type == hol_term_type::OR) {
				if (!process_commutative_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (quantifier_operand->type == hol_term_type::IF_THEN) {
				if (!process_conditional_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable))
					return false;
			} else if (!make_quantifier_scope<QuantifierType>(quantifier, quantifier_operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<hol_term_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, hol_term_type::TRUE);
			}
			recompute_variables(out.noncommutative.left, out.noncommutative.left_negated,
					out.noncommutative.right, out.noncommutative.right_negated, out.variables);
		} else {
			hol_scope& quantifier = *((hol_scope*) alloca(sizeof(hol_scope)));
			recompute_variables(operand->noncommutative.left, operand->noncommutative.left_negated,
					operand->noncommutative.right, operand->noncommutative.right_negated, operand->variables);
			if (!make_quantifier_scope<QuantifierType>(quantifier, operand, quantifier_variable)) {
				free(out); return false;
			}

			bool found_negation;
			if (!add_to_scope<hol_term_type::OR>(quantifier, out.noncommutative.right, out.noncommutative.right_negated, found_negation)) {
				free(out); free(quantifier); return false;
			} else if (found_negation) {
				free(out);
				return init(out, hol_term_type::TRUE);
			}
		}

		if (out.noncommutative.left.length == 0 && out.noncommutative.left_negated.length == 0) {
			/* the antecendent of the new (parent) conditional is empty, so change the node into a disjunction */
			if (out.noncommutative.right.length == 1 && out.noncommutative.right_negated.length == 0) {
				hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
				move(out.noncommutative.right[0], temp);
				out.noncommutative.right.clear();
				free(out); move(temp, out);
			} else if (out.noncommutative.right.length == 0 && out.noncommutative.right_negated.length == 1) {
				hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
				move(out.noncommutative.right_negated[0], temp);
				out.noncommutative.right.clear();
				free(out); move(temp, out);
				if (!negate_scope(out)) {
					free(out); return false;
				}
			} else {
				array<hol_scope>& temp = *((array<hol_scope>*) alloca(sizeof(array<hol_scope>)));
				array<hol_scope>& temp_negated = *((array<hol_scope>*) alloca(sizeof(array<hol_scope>)));
				free(out.noncommutative.left);
				free(out.noncommutative.left_negated);
				move(out.noncommutative.right, temp);
				move(out.noncommutative.right_negated, temp_negated);
				out.type = hol_term_type::OR;
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

template<hol_term_type QuantifierType, bool AllConstantsDistinct>
bool canonicalize_quantifier_scope(
		const hol_quantifier& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	static_assert(QuantifierType == hol_term_type::FOR_ALL
			   || QuantifierType == hol_term_type::EXISTS
			   || QuantifierType == hol_term_type::LAMBDA,
			"QuantifierType is not a quantifier.");

	unsigned int quantifier_variable;
	hol_scope* operand = (hol_scope*) malloc(sizeof(hol_scope));
	if (operand == NULL) {
		fprintf(stderr, "canonicalize_quantifier_scope ERROR: Out of memory.\n");
		return false;
	} else if (!new_variable(src.variable, quantifier_variable, variable_map)
			|| !canonicalize_scope<AllConstantsDistinct>(*src.operand, *operand, variable_map, types))
	{
		free(operand); return false;
	}

	variable_map.size--;

	/* check if the operand has any instances of the quantified variable */
	if (!operand->variables.contains(quantifier_variable)) {
		move(*operand, out); free(operand);
		return true;
	}

	if (operand->type == hol_term_type::AND || operand->type == hol_term_type::OR) {
		return process_commutative_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	} else if (operand->type == hol_term_type::IF_THEN) {
		return process_conditional_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	} else {
		return make_quantifier_scope<QuantifierType>(out, operand, quantifier_variable);
	}
	return true;
}

template<bool AllConstantsDistinct>
inline bool canonicalize_negation_scope(
		const hol_unary_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	if (!canonicalize_scope<AllConstantsDistinct>(*src.operand, out, variable_map, types))
		return false;
	if (!negate_scope(out)) {
		free(out); return false;
	}
	return true;
}

template<hol_term_type Type, bool AllConstantsDistinct>
bool canonicalize_nary_scope(
		const hol_binary_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	hol_scope* left = (hol_scope*) malloc(sizeof(hol_scope));
	if (left == NULL) {
		fprintf(stderr, "canonicalize_nary_scope ERROR: Out of memory.\n");
		free(out); return false;
	}
	hol_scope* right = (hol_scope*) malloc(sizeof(hol_scope));
	if (right == NULL) {
		fprintf(stderr, "canonicalize_nary_scope ERROR: Out of memory.\n");
		free(left); return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.left, *left, variable_map, types)) {
		free(left); free(right); return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.right, *right, variable_map, types)) {
		free(*left); free(left);
		free(right); return false;
	}

	if (!init(out, Type)) {
		free(*left); free(left);
		free(*right); free(right);
		return false;
	}
	out.binary.operands[0] = left;
	out.binary.operands[1] = right;
	move_variables(left->variables, out.variables);
	move_variables(right->variables, out.variables);
	return true;
}

template<hol_term_type Type, bool AllConstantsDistinct>
bool canonicalize_nary_scope(
		const hol_ternary_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	hol_scope* first = (hol_scope*) malloc(sizeof(hol_scope));
	if (first == NULL) {
		fprintf(stderr, "canonicalize_nary_scope ERROR: Out of memory.\n");
		return false;
	}
	hol_scope* second = (hol_scope*) malloc(sizeof(hol_scope));
	if (second == NULL) {
		fprintf(stderr, "canonicalize_nary_scope ERROR: Out of memory.\n");
		free(first); return false;
	}
	hol_scope* third = (hol_scope*) malloc(sizeof(hol_scope));
	if (third == NULL) {
		fprintf(stderr, "canonicalize_nary_scope ERROR: Out of memory.\n");
		free(first); free(second); return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.first, *first, variable_map, types)) {
		free(first); free(second);
		free(third); return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.second, *second, variable_map, types)) {
		free(*first); free(first);
		free(second); free(third); return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.third, *third, variable_map, types)) {
		free(*first); free(first);
		free(*second); free(second);
		free(third); return false;
	}

	if (!init(out, Type)) {
		free(*first); free(first);
		free(*second); free(second);
		free(*third); free(third);
		return false;
	}
	out.ternary.operands[0] = first;
	out.ternary.operands[1] = second;
	out.ternary.operands[2] = third;
	move_variables(first->variables, out.variables);
	move_variables(second->variables, out.variables);
	move_variables(third->variables, out.variables);
	return true;
}

template<bool AllConstantsDistinct>
bool canonicalize_equals_scope(
		const hol_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	hol_scope* left = (hol_scope*) malloc(sizeof(hol_scope));
	if (left == NULL) {
		fprintf(stderr, "canonicalize_equals_scope ERROR: Out of memory.\n");
		return false;
	} else if (!canonicalize_scope<AllConstantsDistinct>(*src.binary.left, *left, variable_map, types)) {
		free(left); return false;
	}

	const pair<hol_type, hol_type>& arg_types = types.types.get(&src);
	bool is_left_boolean = (arg_types.key.kind == hol_type_kind::CONSTANT && arg_types.key.constant == hol_constant_type::BOOLEAN);
	bool is_right_boolean = (arg_types.value.kind == hol_type_kind::CONSTANT && arg_types.value.constant == hol_constant_type::BOOLEAN);
	if (is_right_boolean && left->type == hol_term_type::FALSE) {
		free(*left); free(left);
		if (!canonicalize_scope<AllConstantsDistinct>(*src.binary.right, out, variable_map, types))
			return false;
		if (!negate_scope(out)) { free(out); return false; }
		return true;
	} else if (is_right_boolean && left->type == hol_term_type::TRUE) {
		free(*left); free(left);
		return canonicalize_scope<AllConstantsDistinct>(*src.binary.right, out, variable_map, types);
	} else if (is_right_boolean && left->type == hol_term_type::IFF) {

		hol_scope* right = (hol_scope*) malloc(sizeof(hol_scope));
		if (right == NULL) {
			fprintf(stderr, "canonicalize_equals_scope ERROR: Out of memory.\n");
			free(*left); free(left); return false;
		} else if (!canonicalize_scope<AllConstantsDistinct>(*src.binary.right, *right, variable_map, types)) {
			free(*left); free(left);
			free(right); return false;
		}

		if (right->type == hol_term_type::FALSE) {
			move(*left, out); free(left);
			free(*right); free(right);
			if (!negate_iff(out)) { free(out); return false; }
			return true;
		} else if (right->type == hol_term_type::TRUE) {
			move(*left, out); free(left);
			free(*right); free(right);
			return true;
		} else if (right->type == hol_term_type::IFF) {
			bool found_negation;
			move(*left, out); free(left);
			merge_scopes<hol_term_type::IFF>(right->commutative.children, out.commutative.children,
					right->commutative.negated, out.commutative.negated, found_negation);
			free(right->commutative.children);
			free(right->commutative.negated);
			recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
			free(right->variables); free(right);
			if (found_negation) {
				recompute_variables(out.commutative.children, out.commutative.negated, out.variables);
				if (!negate_iff(out)) { free(out); return false; }
			}
		} else {
			bool found_negation;
			swap(*left, out); free(left);
			if (!add_to_scope<hol_term_type::IFF>(*right, out.commutative.children, out.commutative.negated, out.variables, found_negation)) {
				free(*right); free(right);
				free(out); return false;
			} else if (found_negation) {
				if (!negate_iff(out)) { free(out); return false; }
			}
			free(right);
		}
	} else {

		hol_scope* right = (hol_scope*) malloc(sizeof(hol_scope));
		if (right == NULL) {
			fprintf(stderr, "canonicalize_equals_scope ERROR: Out of memory.\n");
			free(*left); free(left); return false;
		} else if (!canonicalize_scope<AllConstantsDistinct>(*src.binary.right, *right, variable_map, types)) {
			free(*left); free(left);
			free(right); return false;
		}

		if (is_left_boolean && right->type == hol_term_type::FALSE) {
			move(*left, out); free(left);
			free(*right); free(right);
			if (!negate_scope(out)) { free(out); return false; }
			return true;
		} else if (is_left_boolean && right->type == hol_term_type::TRUE) {
			move(*left, out); free(left);
			free(*right); free(right);
			return true;
		} else if (is_left_boolean && right->type == hol_term_type::IFF) {
			bool found_negation;
			move(*right, out); free(right);
			if (!add_to_scope<hol_term_type::IFF>(*left, out.commutative.children, out.commutative.negated, out.variables, found_negation)) {
				free(*left); free(left);
				free(out); return false;
			} else if (found_negation) {
				if (!negate_iff(out)) { free(out); return false; }
			}
			free(left);
		} else if (*right == *left) {
			free(*left); free(left);
			free(*right); free(right);
			return init(out, hol_term_type::TRUE);
		} else if (AllConstantsDistinct && left->type == hol_term_type::CONSTANT
				&& right->type == hol_term_type::CONSTANT && left->constant != right->constant)
		{
			free(*left); free(left);
			free(*right); free(right);
			return init(out, hol_term_type::FALSE);
		} else {
			/* if the types of `left` and `right` are BOOLEAN, then construct
			   an IFF node containing them; otherwise, create an EQUALS node */
			if (is_left_boolean && is_right_boolean) {
				/* child types are BOOLEAN, so construct a IFF node */
				bool first_negation, second_negation;
				if (!init(out, hol_term_type::IFF)) {
					free(*left); free(left);
					free(*right); free(right); return false;
				} else if (!add_to_scope<hol_term_type::IFF>(*left, out.commutative.children, out.commutative.negated, out.variables, first_negation)) {
					free(*left); free(left);
					free(*right); free(right);
					free(out); return false;
				}
				free(left);
				if (!add_to_scope<hol_term_type::IFF>(*right, out.commutative.children, out.commutative.negated, out.variables, second_negation)) {
					free(*right); free(right);
					free(out); return false;
				}
				free(right);
				if (first_negation ^ second_negation) {
					if (!negate_iff(out)) { free(out); return false; }
				}
			} else {
				/* child types are not known to be BOOLEAN, so construct an EQUALS node */
				if (!init(out, hol_term_type::EQUALS)) {
					free(*left); free(left);
					free(*right); free(right); return false;
				} else if (compare(*left, *right) > 0) {
					swap(left, right);
				}
				out.binary.operands[0] = left;
				out.binary.operands[1] = right;
				move_variables(left->variables, out.variables);
				move_variables(right->variables, out.variables);
				return true;
			}
		}
	}

	if (out.commutative.children.length == 0 && out.commutative.negated.length == 0) {
		free(out);
		return init(out, hol_term_type::TRUE);
	} else if (out.commutative.children.length == 1 && out.commutative.negated.length == 0) {
		hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
		move(out.commutative.children[0], temp);
		out.commutative.children.clear();
		free(out); move(temp, out);
	} else if (out.commutative.children.length == 0 && out.commutative.negated.length == 1) {
		hol_scope& temp = *((hol_scope*) alloca(sizeof(hol_scope)));
		move(out.commutative.negated[0], temp);
		out.commutative.negated.clear();
		free(out); move(temp, out);
		if (!negate_scope(out)) {
			free(out); return false;
		}
	}
	return true;
}

template<bool AllConstantsDistinct>
bool canonicalize_scope(const hol_term& src, hol_scope& out,
		array_map<unsigned int, unsigned int>& variable_map,
		const equals_arg_types& types)
{
	unsigned int index;
	switch (src.type) {
	case hol_term_type::AND:
		return canonicalize_commutative_scope<hol_term_type::AND, AllConstantsDistinct>(src.array, out, variable_map, types);
	case hol_term_type::OR:
		return canonicalize_commutative_scope<hol_term_type::OR, AllConstantsDistinct>(src.array, out, variable_map, types);
	case hol_term_type::IFF:
		return canonicalize_commutative_scope<hol_term_type::IFF, AllConstantsDistinct>(src.array, out, variable_map, types);
	case hol_term_type::IF_THEN:
		return canonicalize_conditional_scope<AllConstantsDistinct>(src.binary, out, variable_map, types);
	case hol_term_type::FOR_ALL:
		return canonicalize_quantifier_scope<hol_term_type::FOR_ALL, AllConstantsDistinct>(src.quantifier, out, variable_map, types);
	case hol_term_type::EXISTS:
		return canonicalize_quantifier_scope<hol_term_type::EXISTS, AllConstantsDistinct>(src.quantifier, out, variable_map, types);
	case hol_term_type::LAMBDA:
		return canonicalize_quantifier_scope<hol_term_type::LAMBDA, AllConstantsDistinct>(src.quantifier, out, variable_map, types);
	case hol_term_type::NOT:
		return canonicalize_negation_scope<AllConstantsDistinct>(src.unary, out, variable_map, types);
	case hol_term_type::EQUALS:
		return canonicalize_equals_scope<AllConstantsDistinct>(src, out, variable_map, types);
	case hol_term_type::UNARY_APPLICATION:
		return canonicalize_nary_scope<hol_term_type::UNARY_APPLICATION, AllConstantsDistinct>(src.binary, out, variable_map, types);
	case hol_term_type::BINARY_APPLICATION:
		return canonicalize_nary_scope<hol_term_type::BINARY_APPLICATION, AllConstantsDistinct>(src.ternary, out, variable_map, types);
	case hol_term_type::CONSTANT:
		if (!init(out, hol_term_type::CONSTANT)) return false;
		out.constant = src.constant; return true;
	case hol_term_type::PARAMETER:
		if (!init(out, hol_term_type::PARAMETER)) return false;
		out.parameter = src.parameter; return true;
	case hol_term_type::INTEGER:
		if (!init(out, hol_term_type::INTEGER)) return false;
		out.integer = src.integer; return true;
	case hol_term_type::STRING:
		if (!init(out, hol_term_type::STRING)) return false;
		out.str = src.str; return true;
	case hol_term_type::UINT_LIST:
		if (!init(out, hol_term_type::UINT_LIST)) return false;
		out.uint_list = src.uint_list; return true;
	case hol_term_type::VARIABLE:
		if (!init(out, hol_term_type::VARIABLE)) return false;
		index = variable_map.index_of(src.variable);
		if (index < variable_map.size) {
			out.variable = variable_map.values[index];
		} else if (!new_variable(src.variable, out.variable, variable_map)) {
			return false;
		}

		index = linear_search(out.variables.data, out.variable, 0, out.variables.length);
		if (index == out.variables.length) {
			if (!out.variables.ensure_capacity(out.variables.length + 1)) return false;
			out.variables[out.variables.length] = out.variable;
			out.variables.length++;
		} else if (out.variables[index] != out.variable) {
			if (!out.variables.ensure_capacity(out.variables.length + 1)) return false;
			shift_right(out.variables.data, out.variables.length, index);
			out.variables[index] = out.variable;
			out.variables.length++;
		}
		return true;
	case hol_term_type::TRUE:
		return init(out, hol_term_type::TRUE);
	case hol_term_type::FALSE:
		return init(out, hol_term_type::FALSE);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		fprintf(stderr, "canonicalize_scope ERROR: Canonicalization of formulas with expressions of type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` are not supported.\n");
		exit(EXIT_FAILURE); /* we don't support canonicalization of expressions with type `ANY`, `ANY_RIGHT`, `ANY_ARRAY`, `ANY_CONSTANT`, or `ANY_QUANTIFIER` */
	}
	fprintf(stderr, "canonicalize_scope ERROR: Unrecognized hol_term_type.\n");
	return false;
}

struct identity_canonicalizer { };

inline hol_term* canonicalize(hol_term& src,
		const identity_canonicalizer& canonicalizer)
{
	src.reference_count++;
	return &src;
}

template<bool AllConstantsDistinct, bool PolymorphicEquality>
struct standard_canonicalizer { };

template<bool AllConstantsDistinct, bool PolymorphicEquality>
inline hol_term* canonicalize(const hol_term& src,
		const standard_canonicalizer<AllConstantsDistinct, PolymorphicEquality>& canonicalizer)
{
	equals_arg_types types(16);
	if (!compute_type<PolymorphicEquality>(src, types))
		return NULL;

	array_map<unsigned int, unsigned int> variable_map(16);
	hol_scope& scope = *((hol_scope*) alloca(sizeof(hol_scope)));
	if (!canonicalize_scope<AllConstantsDistinct>(src, scope, variable_map, types))
		return NULL;
	hol_term* canonicalized = scope_to_term(scope);
	free(scope);
	return canonicalized;
}

template<typename Canonicalizer>
bool is_canonical(const hol_term& src, Canonicalizer& canonicalizer) {
	hol_term* canonicalized = canonicalize(src, canonicalizer);
	if (canonicalized == NULL) {
		fprintf(stderr, "is_canonical ERROR: Unable to canonicalize term.\n");
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
		const hol_term& src, identity_canonicalizer& canonicalizer)
{
	return true;
}


/**
 * Code for determining set relations with sets of the form {x : A} where A is
 * a higher-order logic formula.
 */

/* forward declarations */
bool is_subset(const hol_term* first, const hol_term* second);


bool is_conjunction_subset(
		const hol_term* const* first, unsigned int first_length,
		const hol_term* const* second, unsigned int second_length)
{
	unsigned int i = 0, j = 0;
	while (i < first_length && j < second_length)
	{
		if (first[i] == second[j] || *first[i] == *second[j]) {
			i++; j++;
		} else if (*first[i] < *second[j]) {
			i++;
		} else {
			bool found_subset = false;
			for (unsigned int k = 0; k < first_length; k++) {
				if (is_subset(first[k], second[j])) {
					j++; found_subset = true;
					break;
				}
			}
			if (!found_subset) return false;
		}
	}

	while (j < second_length) {
		bool found_subset = false;
		for (unsigned int k = 0; k < first_length; k++) {
			if (is_subset(first[k], second[j])) {
				j++; found_subset = true;
				break;
			}
		}
		if (!found_subset) return false;
	}

	return true;
}

bool is_disjunction_subset(
		const hol_term* const* first, unsigned int first_length,
		const hol_term* const* second, unsigned int second_length)
{
	unsigned int i = 0, j = 0;
	while (i < first_length && j < second_length)
	{
		if (first[i] == second[j] || *first[i] == *second[j]) {
			i++; j++;
		} else if (*first[i] < *second[j]) {
			bool found_subset = false;
			for (unsigned int k = 0; k < second_length; k++) {
				if (is_subset(first[i], second[k])) {
					i++; found_subset = true;
					break;
				}
			}
			if (!found_subset) return false;
		} else {
			j++;
		}
	}

	while (i < first_length) {
		bool found_subset = false;
		for (unsigned int k = 0; k < second_length; k++) {
			if (is_subset(first[i], second[k])) {
				i++; found_subset = true;
				break;
			}
		}
		if (!found_subset) return false;
	}

	return true;
}

bool is_subset(const hol_term* first, const hol_term* second)
{
	if (first->type == hol_term_type::TRUE) {
		return (second->type == hol_term_type::TRUE);
	} else if (second->type == hol_term_type::TRUE) {
		return true;
	} else if (first->type == hol_term_type::FALSE) {
		return true;
	} else if (second->type == hol_term_type::FALSE) {
		return (first->type == hol_term_type::FALSE);
	} else if (first->type == hol_term_type::AND) {
		if (second->type == hol_term_type::AND) {
			return is_conjunction_subset(first->array.operands, first->array.length, second->array.operands, second->array.length);
		} else {
			return is_conjunction_subset(first->array.operands, first->array.length, &second, 1);
		}
	} else if (second->type == hol_term_type::AND) {
		return false;
	} else if (first->type == hol_term_type::OR) {
		if (second->type == hol_term_type::OR) {
			return is_disjunction_subset(first->array.operands, first->array.length, second->array.operands, second->array.length);
		} else {
			return false;
		}
	} else if (second->type == hol_term_type::OR) {
		return is_disjunction_subset(&first, 1, second->array.operands, second->array.length);
	}

	switch (first->type) {
	case hol_term_type::CONSTANT:
		return (second->type == hol_term_type::CONSTANT && first->constant == second->constant);
	case hol_term_type::VARIABLE:
		return (second->type == hol_term_type::VARIABLE && first->variable == second->variable);
	case hol_term_type::PARAMETER:
		return (second->type == hol_term_type::PARAMETER && first->parameter == second->parameter);
	case hol_term_type::NOT:
		if (second->type == hol_term_type::NOT) {
			return is_subset(second->unary.operand, first->unary.operand);
		} else {
			return false;
		}
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
		return *first == *second;
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::IFF:
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		/* TODO: finish implementing this */
		fprintf(stderr, "is_subset ERROR: Not implemented.\n");
		exit(EXIT_FAILURE);

	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
		fprintf(stderr, "is_subset ERROR: `first` does not have type proposition.\n");
		exit(EXIT_FAILURE);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		/* this should be unreachable */
		break;
	}
	fprintf(stderr, "is_subset ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

inline hol_term* intersect(hol_term* first, hol_term* second)
{
	standard_canonicalizer<false, false> canonicalizer;
	hol_term* conjunction = hol_term::new_and(first, second);
	first->reference_count++; second->reference_count++;
	hol_term* canonicalized = canonicalize(*conjunction, canonicalizer);
	free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
	return canonicalized;
}


/**
 * Code to facilitate semantic parsing and generation.
 */


/* forward declarations */

bool is_reduceable(hol_term*, hol_term*, hol_term*&, hol_term*&, unsigned int&);

template<typename BuiltInPredicates, typename VariableComparator>
bool is_subset(hol_term*, hol_term*, VariableComparator&);

template<typename BuiltInPredicates>
bool subtract(array<hol_term*>&, hol_term*, hol_term*);

template<typename BuiltInPredicates>
bool subtract_any(array<hol_term*>&, hol_term*, hol_term*);

template<typename BuiltInPredicates, bool ComputeIntersection = true>
bool intersect_with_any(array<hol_term*>&, hol_term*, hol_term*);

template<typename BuiltInPredicates, bool ComputeIntersection = true>
bool intersect(array<hol_term*>&, hol_term*, hol_term*);


inline bool any_number(const hol_term& src) {
	return src.type == hol_term_type::ANY;
}

/* NOTE: this function assumes src is not ANY */
inline bool get_number(const hol_term& src, int& value) {
	if (src.type != hol_term_type::INTEGER)
		return false;
	value = src.integer;
	return true;
}

inline bool set_number(hol_term& exp,
		const hol_term& set, int value)
{
	if (set.type != hol_term_type::ANY && set.type != hol_term_type::INTEGER)
		return false;
	exp.type = hol_term_type::INTEGER;
	exp.integer = value;
	exp.reference_count = 1;
	return true;
}

inline bool any_uint_list(const hol_term& src) {
	return src.type == hol_term_type::ANY;
}

/* NOTE: this function assumes src is not ANY */
inline bool get_uint_list(const hol_term& src, sequence& value) {
	if (src.type != hol_term_type::UINT_LIST)
		return false;
	return init(value, src.uint_list);
}

inline bool set_uint_list(hol_term& exp,
		const hol_term& set, const sequence& value)
{
	if (set.type != hol_term_type::ANY && set.type != hol_term_type::UINT_LIST)
		return false;
	exp.type = hol_term_type::UINT_LIST;
	exp.uint_list = value;
	exp.reference_count = 1;
	return true;
}

inline bool is_reduceable(
		hol_array_term& first_src,
		hol_array_term& second_src,
		hol_term*& first_node,
		hol_term*& second_node,
		unsigned int& prefix_index)
{
	if (first_src.length != second_src.length) return false;
	for (unsigned int i = 0; i < first_src.length; i++)
		if (!is_reduceable(first_src.operands[i], second_src.operands[i], first_node, second_node, prefix_index)) return false;
	return true;
}

bool intersect(hol_quantifier_type& dst, hol_quantifier_type first, hol_quantifier_type second) {
	if (first == hol_quantifier_type::ANY || first == second) {
		dst = second;
	} else if (second == hol_quantifier_type::ANY) {
		dst = first;
	} else if (first == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
		if (second == hol_quantifier_type::LAMBDA) {
			return false;
		} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else {
			dst = second;
		}
	} else if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
		if (first == hol_quantifier_type::LAMBDA) {
			return false;
		} else if (first == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else if (first == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
		if (second == hol_quantifier_type::EXISTS) {
			return false;
		} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::LAMBDA;
		} else {
			dst = second;
		}
	} else if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
		if (first == hol_quantifier_type::EXISTS) {
			return false;
		} else if (first == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::LAMBDA;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::EXISTS_OR_LAMBDA) {
		if (second == hol_quantifier_type::FOR_ALL) {
			return false;
		} else {
			dst = second;
		}
	} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
		if (first == hol_quantifier_type::FOR_ALL) {
			return false;
		} else {
			dst = first;
		}
	} else {
		return false;
	}
	return true;
}

inline bool has_intersection(hol_quantifier_type first, hol_quantifier_type second) {
	hol_quantifier_type dummy;
	return intersect(dummy, first, second);
}

bool is_subset(hol_quantifier_type first, hol_quantifier_type second) {
	if (first == second || second == hol_quantifier_type::ANY) {
		return true;
	} else if (first == hol_quantifier_type::ANY) {
		return false;
	} else if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
		return (first == hol_quantifier_type::FOR_ALL || first == hol_quantifier_type::EXISTS);
	} else if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
		return (first == hol_quantifier_type::FOR_ALL || first == hol_quantifier_type::LAMBDA);
	} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
		return (first == hol_quantifier_type::EXISTS || first == hol_quantifier_type::LAMBDA);
	} else {
		return false;
	}
}

bool subtract(hol_quantifier_type& dst, hol_quantifier_type first, hol_quantifier_type second) {
	if (first == second || second == hol_quantifier_type::ANY) {
		return false;
	} else if (first == hol_quantifier_type::ANY) {
		if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
			dst = hol_quantifier_type::LAMBDA;
		} else if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else if (second == hol_quantifier_type::FOR_ALL) {
			dst = hol_quantifier_type::EXISTS_OR_LAMBDA;
		} else if (second == hol_quantifier_type::EXISTS) {
			dst = hol_quantifier_type::FOR_ALL_OR_LAMBDA;
		} else if (second == hol_quantifier_type::LAMBDA) {
			dst = hol_quantifier_type::EXISTS_OR_LAMBDA;
		}
#if !defined(NDEBUG)
		else fprintf(stderr, "subtract WARNING: Unrecognized hol_quantifier_type.\n");
#endif
	} else if (first == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
		if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else if (second == hol_quantifier_type::FOR_ALL) {
			dst = hol_quantifier_type::EXISTS;
		} else if (second == hol_quantifier_type::EXISTS) {
			dst = hol_quantifier_type::FOR_ALL;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
		if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
			dst = hol_quantifier_type::LAMBDA;
		} else if (second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else if (second == hol_quantifier_type::FOR_ALL) {
			dst = hol_quantifier_type::LAMBDA;
		} else if (second == hol_quantifier_type::LAMBDA) {
			dst = hol_quantifier_type::FOR_ALL;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::EXISTS_OR_LAMBDA) {
		if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS) {
			dst = hol_quantifier_type::LAMBDA;
		} else if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else if (second == hol_quantifier_type::EXISTS) {
			dst = hol_quantifier_type::LAMBDA;
		} else if (second == hol_quantifier_type::LAMBDA) {
			dst = hol_quantifier_type::EXISTS;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::FOR_ALL) {
		if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS || second == hol_quantifier_type::FOR_ALL_OR_LAMBDA) {
			return false;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::EXISTS) {
		if (second == hol_quantifier_type::FOR_ALL_OR_EXISTS || second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			return false;
		} else {
			dst = first;
		}
	} else if (first == hol_quantifier_type::LAMBDA) {
		if (second == hol_quantifier_type::FOR_ALL_OR_LAMBDA || second == hol_quantifier_type::EXISTS_OR_LAMBDA) {
			return false;
		} else {
			dst = first;
		}
	}
	return true;
}

bool is_reduceable(
		hol_term* first_src,
		hol_term* second_src,
		hol_term*& first_node,
		hol_term*& second_node,
		unsigned int& prefix_index)
{
	if (first_src->type != second_src->type) {
		if (first_node != nullptr) {
			return false;
		} else if (first_src->type == hol_term_type::ANY || second_src->type == hol_term_type::ANY) {
			first_node = first_src;
			second_node = second_src;
			return true;
		} else {
			return false;
		}
	} else if (first_src->type == hol_term_type::ANY || first_src->type == hol_term_type::ANY_RIGHT) {
		if ((first_src->any.included == nullptr || second_src->any.included == nullptr)
		 && first_src->any.included != second_src->any.included)
		{
			if (first_node != nullptr) return false;
			first_node = first_src;
			second_node = second_src;
			return true;
		} if (first_src->any.excluded_tree_count != second_src->any.excluded_tree_count)
		{
			if (first_node != nullptr) return false;
			first_node = first_src;
			second_node = second_src;
			return true;
		}
		for (unsigned int i = 0; i < first_src->any.excluded_tree_count; i++) {
			if (first_src->any.excluded_trees[i] != second_src->any.excluded_trees[i]
			 && *first_src->any.excluded_trees[i] != *second_src->any.excluded_trees[i])
			{
				if (first_node != nullptr) return false;
				first_node = first_src;
				second_node = second_src;
				return true;
			}
		}
		hol_term* first_inner_node = nullptr;
		hol_term* second_inner_node = nullptr;
		unsigned int inner_prefix_index = 0;
		if (first_src->any.included != nullptr && !is_reduceable(first_src->any.included, second_src->any.included, first_inner_node, second_inner_node, inner_prefix_index)) {
			if (first_node != nullptr) return false;
			first_node = first_src;
			second_node = second_src;
		} else if (first_inner_node != nullptr) {
			first_node = first_inner_node;
			second_node = second_inner_node;
			prefix_index += inner_prefix_index;
		}
		return true;
	}

	if (first_node == nullptr) prefix_index++;
	switch (first_src->type) {
	case hol_term_type::CONSTANT:
		return first_src->constant == second_src->constant;
	case hol_term_type::VARIABLE:
		return first_src->variable == second_src->variable;
	case hol_term_type::PARAMETER:
		return first_src->parameter == second_src->parameter;
	case hol_term_type::INTEGER:
		return first_src->integer == second_src->integer;
	case hol_term_type::STRING:
		return first_src->str == second_src->str;
	case hol_term_type::UINT_LIST:
		return first_src->uint_list == second_src->uint_list;
	case hol_term_type::NOT:
		return is_reduceable(first_src->unary.operand, second_src->unary.operand, first_node, second_node, prefix_index);
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return is_reduceable(first_src->binary.left, second_src->binary.left, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->binary.right, second_src->binary.right, first_node, second_node, prefix_index);
	case hol_term_type::BINARY_APPLICATION:
		return is_reduceable(first_src->ternary.first, second_src->ternary.first, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->ternary.second, second_src->ternary.second, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->ternary.third, second_src->ternary.third, first_node, second_node, prefix_index);
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		return is_reduceable(first_src->array, second_src->array, first_node, second_node, prefix_index);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:		
	case hol_term_type::LAMBDA:
		return (first_src->quantifier.variable != second_src->quantifier.variable)
			&& is_reduceable(first_src->quantifier.operand, second_src->quantifier.operand, first_node, second_node, prefix_index);
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	case hol_term_type::ANY_ARRAY:
		return (first_src->any_array.oper == hol_term_type::ANY_ARRAY || second_src->any_array.oper == hol_term_type::ANY_ARRAY || first_src->any_array.oper == second_src->any_array.oper)
			&& is_reduceable(first_src->any_array.all, second_src->any_array.all, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->any_array.left, second_src->any_array.left, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->any_array.any, second_src->any_array.any, first_node, second_node, prefix_index)
			&& is_reduceable(first_src->any_array.right, second_src->any_array.right, first_node, second_node, prefix_index);
	case hol_term_type::ANY_CONSTANT:
		return first_src->any_constant == second_src->any_constant;
	case hol_term_type::ANY_QUANTIFIER:
		return has_intersection(first_src->any_quantifier.quantifier, second_src->any_quantifier.quantifier)
			&& is_reduceable(first_src->any_quantifier.operand, second_src->any_quantifier.operand, first_node, second_node, prefix_index);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		break; /* we already handle this case above */
	}
	fprintf(stderr, "is_reduceable ERROR: Unrecognized hol_term_type.\n");
	return false;
}

struct reducer {
	unsigned int prefix_index;
	const unsigned int target_prefix_index;
	const unsigned int* kept_tree_indices;
	const unsigned int kept_tree_count;
};

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, reducer& r) {
#if defined(NDEBUG)
	if (Type == hol_term_type::ANY && r.prefix_index == r.target_prefix_index) {
#else
	if (r.prefix_index == r.target_prefix_index) {
		if (Type != hol_term_type::ANY)
			fprintf(stderr, "apply WARNING: Expected an `ANY` type during reduction.\n");
#endif

		if (src->any.included == nullptr && r.kept_tree_count == 0) {
			HOL_ANY.reference_count++;
			return &HOL_ANY;
		}

		hol_term* new_node;
		if (!new_hol_term(new_node)) return nullptr;
		new_node->type = hol_term_type::ANY;
		new_node->reference_count = 1;
		new_node->any.excluded_tree_count = 0;
		if (r.kept_tree_count == 0) {
			new_node->any.excluded_trees = nullptr;
		} else {
			new_node->any.excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * r.kept_tree_count);
			if (new_node->any.excluded_trees == nullptr) {
				free(new_node); return nullptr;
			}
			for (unsigned int i = 0; i < r.kept_tree_count; i++) {
				new_node->any.excluded_trees[new_node->any.excluded_tree_count] = src->any.excluded_trees[r.kept_tree_indices[i]];
				new_node->any.excluded_trees[new_node->any.excluded_tree_count++]->reference_count++;
			}
		}
		new_node->any.included = src->any.included;
		new_node->any.included->reference_count++;
	}

	r.prefix_index++;
	return default_apply<Type>(src, r);
}

inline hol_term* reduce(hol_term* src, const unsigned int target_prefix_index,
		const unsigned int* kept_tree_indices, const unsigned int kept_tree_count)
{
	const reducer r = {0, target_prefix_index, kept_tree_indices, kept_tree_count};
	hol_term* dst = apply(src, r);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

/**
 * Returns whether `C(first_subtree)` is a subset of `second`.
 */
template<typename BuiltInPredicates>
inline bool any_is_subset(hol_term* first_subtree, hol_term* second)
{
	if (second->type == hol_term_type::ANY || second->type == hol_term_type::ANY_RIGHT) {
		for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
			array<hol_term*> dummy(1);
			hol_term* any = hol_term::new_any(first_subtree);
			if (any == nullptr) return false;
			first_subtree->reference_count++;
			bool non_empty_intersection = intersect_with_any<BuiltInPredicates, false>(dummy, second->any.excluded_trees[i], any);
			free(*any); if (any->reference_count == 0) free(any);
			if (non_empty_intersection) return false;
		}
		return is_subset<BuiltInPredicates>(first_subtree, second);
	} else if (second->type == hol_term_type::ANY_RIGHT) {
		if (second->any.included != nullptr)
			return false;
		for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
			array<hol_term*> dummy(1);
			hol_term* any = hol_term::new_any_right(first_subtree);
			if (any == nullptr) return false;
			first_subtree->reference_count++;
			bool non_empty_intersection = intersect_with_any<BuiltInPredicates, false>(dummy, second->any.excluded_trees[i], any);
			free(*any); if (any->reference_count == 0) free(any);
			if (non_empty_intersection) return false;
		}
		return true;
	} else {
		return false;
	}
}

/**
 * Returns whether `C(first_subtree)` is a subset of the union of `second[i]`.
 */
template<typename BuiltInPredicates>
inline bool any_is_subset(hol_term* first_subtree, hol_term** second, unsigned int second_length)
{
	/* NOTE: we assume that `second` is in reduced form */
	for (unsigned int i = 0; i < second_length; i++)
		if (any_is_subset<BuiltInPredicates>(first_subtree, second[i])) return true;
	return false;
}

template<typename BuiltInPredicates, bool ExpandedSetInFirstSet>
inline bool reduce_union(
		hol_term** first, unsigned int& first_length,
		hol_term** second, unsigned int& second_length,
		unsigned int& expanded_set_index)
{
	hol_term* expanded_set = (ExpandedSetInFirstSet ? first[expanded_set_index] : second[expanded_set_index]);
	for (unsigned int i = 0; i < first_length; i++) {
		if (ExpandedSetInFirstSet && expanded_set_index == i) continue;

		/* consider all possible reductions between the trees `first[i]` and `expanded_set` */
		/* first check to see that `first[i]` and `expanded_set` differ in all but one node */
		unsigned int prefix_index = 0;
		hol_term* node = nullptr; hol_term* expanded_node = nullptr;
		if (!is_reduceable(first[i], expanded_set, node, expanded_node, prefix_index)) {
			continue;
		} else if (node == nullptr) {
			/* `first[i]` and `expanded_set` are identical */
			free(*first[i]);
			if (first[i]->reference_count == 0)
				free(first[i]);
			first[i] = first[--first_length];
			return true;
		}

		/* they differ in only one node (and one of them is ANY), so try reducing `first[i]` and `expanded_set` */
		if (node->type == hol_term_type::ANY) {
			unsigned int kept_tree_count = 0;
			unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * node->any.excluded_tree_count));
			if (kept_tree_indices == nullptr)
				return false;
			for (unsigned int k = 0; k < node->any.excluded_tree_count; k++) {
				if (!any_is_subset<BuiltInPredicates>(node->any.excluded_trees[k], expanded_node))
					kept_tree_indices[kept_tree_count++] = k;
			}

			if (kept_tree_count != node->any.excluded_tree_count) {
				/* we've removed some excluded sets from `node` so create a new `first[i]` */
				hol_term* new_set = reduce(first[i], prefix_index, kept_tree_indices, kept_tree_count);
				free(kept_tree_indices);
				if (new_set == nullptr) return false;
				free(*first[i]);
				if (first[i]->reference_count == 0)
					free(first[i]);
				first[i] = new_set;

				/* since `first_node` has not become a larger set, it may now be reducible with sets that we've already considered */
				if (!reduce_union<BuiltInPredicates, true>(first, first_length, second, second_length, i)) return false;
			} else {
				free(kept_tree_indices);
			}
		}
	}
	for (unsigned int i = 0; i < second_length; i++) {
		if (!ExpandedSetInFirstSet && expanded_set_index == i) continue;

		/* consider all possible reductions between the trees `second[i]` and `expanded_set` */
		/* first check to see that `second[i]` and `expanded_set` differ in all but one node */
		unsigned int prefix_index = 0;
		hol_term* node = nullptr; hol_term* expanded_node = nullptr;
		if (!is_reduceable(second[i], expanded_set, node, expanded_node, prefix_index)) {
			continue;
		} else if (node == nullptr) {
			/* `second[i]` and `expanded_set` are identical */
			free(*second[i]);
			if (second[i]->reference_count == 0)
				free(second[i]);
			second[i] = second[--second_length];
			return true;
		}

		/* they differ in only one node (and one of them is ANY), so try reducing `second[i]` and `expanded_set` */
		if (node->type == hol_term_type::ANY) {
			unsigned int kept_tree_count = 0;
			unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * node->any.excluded_tree_count));
			if (kept_tree_indices == nullptr)
				return false;
			for (unsigned int k = 0; k < node->any.excluded_tree_count; k++) {
				if (!any_is_subset<BuiltInPredicates>(node->any.excluded_trees[k], expanded_node))
					kept_tree_indices[kept_tree_count++] = k;
			}

			if (kept_tree_count != node->any.excluded_tree_count) {
				/* we've removed some excluded sets from `node` so create a new `second[i]` */
				hol_term* new_set = reduce(second[i], prefix_index, kept_tree_indices, kept_tree_count);
				free(kept_tree_indices);
				if (new_set == nullptr) return false;
				free(*second[i]);
				if (second[i]->reference_count == 0)
					free(second[i]);
				second[i] = new_set;

				/* since `second_node` has not become a larger set, it may now be reducible with sets that we've already considered */
				if (!reduce_union<BuiltInPredicates, false>(first, first_length, second, second_length, i)) return false;
			} else {
				free(kept_tree_indices);
			}
		}
	}
	return true;
}

template<typename BuiltInPredicates>
bool reduce_union(
		hol_term** first, unsigned int& first_length,
		hol_term** second, unsigned int& second_length)
{
	for (unsigned int i = 0; i < second_length; i++) {
		for (unsigned int j = 0; j < first_length; j++) {
			/* consider all possible reductions between the trees `second[i]` and `first[j]` */
			/* first check to see that `second[i]` and `first[j]` differ in all but one node */
			unsigned int prefix_index = 0;
			hol_term* first_node = nullptr; hol_term* second_node = nullptr;
			if (!is_reduceable(second[i], first[j], second_node, first_node, prefix_index)) {
				continue;
			} else if (first_node == nullptr) {
				/* `second[i]` and `first[j]` are identical */
				free(*second[i]);
				if (second[i]->reference_count == 0)
					free(second[i]);
				second[i--] = second[--second_length];
				break;
			}

			/* they differ in only one node (and one of them is ANY), so try reducing `src[i]` and `dst[j]` */
			if (first_node->type == hol_term_type::ANY) {
				unsigned int kept_tree_count = 0;
				unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * first_node->any.excluded_tree_count));
				if (kept_tree_indices == nullptr)
					return false;
				for (unsigned int k = 0; k < first_node->any.excluded_tree_count; k++) {
					if (!any_is_subset<BuiltInPredicates>(first_node->any.excluded_trees[k], second_node))
						kept_tree_indices[kept_tree_count++] = k;
				}

				if (kept_tree_count != first_node->any.excluded_tree_count) {
					/* we've removed some excluded sets from `first_node` so create a new `first[j]` */
					hol_term* new_set = reduce(first[j], prefix_index, kept_tree_indices, kept_tree_count);
					free(kept_tree_indices);
					if (new_set == nullptr) return false;
					free(*first[j]);
					if (first[j]->reference_count == 0)
						free(first[j]);
					first[j] = new_set;

					/* since `first_node` has not become a larger set, it may now be reducible with sets that we've already considered */
					unsigned int new_second_length = i + 1;
					if (!reduce_union<BuiltInPredicates, true>(first, first_length, second, new_second_length, j)) return false;
					/* some elements in `second` need to be shifted to the left by `i + 1 - new_second_length` */
					for (unsigned int j = i + 1; j < second_length; j++)
						second[j - (i + 1 - new_second_length)] = second[j];
					second_length -= (i + 1 - new_second_length);
				} else {
					free(kept_tree_indices);
				}
			} if (second_node->type == hol_term_type::ANY) {
				unsigned int kept_tree_count = 0;
				unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * second_node->any.excluded_tree_count));
				if (kept_tree_indices == nullptr)
					return false;
				for (unsigned int k = 0; k < second_node->any.excluded_tree_count; k++) {
					if (!any_is_subset<BuiltInPredicates>(second_node->any.excluded_trees[k], first_node))
						kept_tree_indices[kept_tree_count++] = k;
				}

				if (kept_tree_count != second_node->any.excluded_tree_count) {
					/* we've removed some excluded sets from `second_node` so create a new `first[j]` */
					hol_term* new_set = reduce(second[i], prefix_index, kept_tree_indices, kept_tree_count);
					free(kept_tree_indices);
					if (new_set == nullptr) return false;
					free(*second[i]);
					if (second[i]->reference_count == 0)
						free(second[i]);
					second[i] = new_set;

					/* since `second_node` has not become a larger set, it may now be reducible with sets that we've already considered */
					unsigned int new_second_length = i + 1;
					if (!reduce_union<BuiltInPredicates, false>(first, first_length, second, new_second_length, i)) return false;
					/* some elements in `second` need to be shifted to the left by `i + 1 - new_second_length` */
					for (unsigned int j = i + 1; j < second_length; j++)
						second[j - (i + 1 - new_second_length)] = second[j];
					second_length -= (i + 1 - new_second_length);
				} else {
					free(kept_tree_indices);
				}
			}
		}
	}
	return true;
}

template<typename BuiltInPredicates>
bool reduce_union(
		hol_term** sets, unsigned int& length)
{
	for (unsigned int i = 1; i < length; i++) {
		for (unsigned int j = 0; j < i; j++) {
			/* consider all possible reductions between the trees `sets[i]` and `sets[j]` */
			/* first check to see that `sets[i]` and `sets[j]` differ in all but one node */
			unsigned int prefix_index = 0;
			hol_term* first_node = nullptr; hol_term* second_node = nullptr;
			if (!is_reduceable(sets[i], sets[j], first_node, second_node, prefix_index)) {
				continue;
			} else if (first_node == nullptr) {
				/* `sets[i]` and `sets[j]` are identical */
				free(*sets[i]);
				if (sets[i]->reference_count == 0)
					free(sets[i]);
				sets[i--] = sets[--length];
				break;
			}

			/* they differ in only one node (and one of them is ANY), so try reducing `sets[i]` and `sets[j]` */
			if (first_node->type == hol_term_type::ANY) {
				unsigned int kept_tree_count = 0;
				unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * first_node->any.excluded_tree_count));
				if (kept_tree_indices == nullptr)
					return false;
				for (unsigned int k = 0; k < first_node->any.excluded_tree_count; k++) {
					if (!any_is_subset<BuiltInPredicates>(first_node->any.excluded_trees[k], second_node))
						kept_tree_indices[kept_tree_count++] = k;
				}

				if (kept_tree_count != first_node->any.excluded_tree_count)
				{
					/* we've removed some excluded sets from `first_node` so create a new `sets[i]` */
					hol_term* new_set = reduce(sets[i], prefix_index, kept_tree_indices, kept_tree_count);
					free(kept_tree_indices);
					if (new_set == nullptr) return false;
					free(*sets[i]);
					if (sets[i]->reference_count == 0)
						free(sets[i]);
					sets[i] = new_set;

					/* since `first_node` has not become a larger set, it may now be reducible with sets that we've already considered */
					unsigned int dummy = 0;
					if (!reduce_union<BuiltInPredicates, true>(sets, length, nullptr, dummy, i)) return false;
				} else {
					free(kept_tree_indices);
				}
			} if (second_node->type == hol_term_type::ANY) {
				unsigned int kept_tree_count = 0;
				unsigned int* kept_tree_indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * second_node->any.excluded_tree_count));
				if (kept_tree_indices == nullptr)
					return false;
				for (unsigned int k = 0; k < second_node->any.excluded_tree_count; k++) {
					if (!any_is_subset<BuiltInPredicates>(second_node->any.excluded_trees[k], first_node))
						kept_tree_indices[kept_tree_count++] = k;
				}

				if (kept_tree_count != second_node->any.excluded_tree_count)
				{
					/* we've removed some excluded sets from `second_node` so create a new `sets[j]` */
					hol_term* new_set = reduce(sets[j], prefix_index, kept_tree_indices, kept_tree_count);
					free(kept_tree_indices);
					if (new_set == nullptr) return false;
					free(*sets[j]);
					if (sets[j]->reference_count == 0)
						free(sets[j]);
					sets[j] = new_set;

					/* since `second_node` has not become a larger set, it may now be reducible with sets that we've already considered */
					unsigned int dummy = 0;
					if (!reduce_union<BuiltInPredicates, true>(sets, length, nullptr, dummy, j)) return false;
				} else {
					free(kept_tree_indices);
				}
			}
		}
	}
	return true;
}

template<typename BuiltInPredicates>
bool is_subset(hol_term* first, hol_term* second)
{
	if (first == second) {
		return true;
	} else if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) {
		if (second->type == hol_term_type::ANY || second->type == hol_term_type::ANY_RIGHT) {
			if (first->type == hol_term_type::ANY && second->type == hol_term_type::ANY_RIGHT && second->any.included != nullptr)
				return false;
			unsigned int first_tree_union_length = first->any.excluded_tree_count;
			unsigned int second_tree_union_length = 1;
			hol_term** first_tree_union = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * first_tree_union_length));
			hol_term** second_tree_union = &second;
			if (first_tree_union == nullptr) {
				fprintf(stderr, "is_subset ERROR: Out of memory.\n");
				return false;
			}
			for (unsigned int i = 0; i < first_tree_union_length; i++) {
				first_tree_union[i] = first->any.excluded_trees[i];
				first_tree_union[i]->reference_count++;
			}
			second_tree_union[0]->reference_count++;
			if (!reduce_union<BuiltInPredicates>(first_tree_union, first_tree_union_length, second_tree_union, second_tree_union_length)) {
				fprintf(stderr, "is_subset ERROR: `reduce_union` failed.\n");
				for (unsigned int i = 0; i < first_tree_union_length; i++) {
					free(*first_tree_union[i]);
					if (first_tree_union[i]->reference_count == 0)
						free(first_tree_union[i]);
				}
				free(*second_tree_union[0]);
				if (second_tree_union[0]->reference_count == 0)
					free(second_tree_union[0]);
				free(first_tree_union);
				return false;
			}

			if (second_tree_union_length == 0) {
				for (unsigned int i = 0; i < first_tree_union_length; i++) {
					free(*first_tree_union[i]);
					if (first_tree_union[i]->reference_count == 0)
						free(first_tree_union[i]);
				}
				free(first_tree_union);
				return false;
			} else if (first->any.included == nullptr) {
				if (second_tree_union[0]->any.included != nullptr) {
					for (unsigned int i = 0; i < first_tree_union_length; i++) {
						free(*first_tree_union[i]);
						if (first_tree_union[i]->reference_count == 0)
							free(first_tree_union[i]);
					}
					free(*second_tree_union[0]);
					if (second_tree_union[0]->reference_count == 0)
						free(second_tree_union[0]);
					free(first_tree_union);
					return false;
				}
			} else if (!is_subset<BuiltInPredicates>(first->any.included, second_tree_union[0])) {
				for (unsigned int i = 0; i < first_tree_union_length; i++) {
					free(*first_tree_union[i]);
					if (first_tree_union[i]->reference_count == 0)
						free(first_tree_union[i]);
				}
				free(*second_tree_union[0]);
				if (second_tree_union[0]->reference_count == 0)
					free(second_tree_union[0]);
				free(first_tree_union);
				return false;
			}

			bool result = true;
			for (unsigned int i = 0; i < second_tree_union[0]->any.excluded_tree_count; i++) {
				if (!is_subset<BuiltInPredicates>(second_tree_union[0]->any.excluded_trees[i], first_tree_union, first_tree_union_length)) {
					result = false;
					break;
				}
			}
			for (unsigned int i = 0; i < first_tree_union_length; i++) {
				free(*first_tree_union[i]);
				if (first_tree_union[i]->reference_count == 0)
					free(first_tree_union[i]);
			}
			free(*second_tree_union[0]);
			if (second_tree_union[0]->reference_count == 0)
				free(second_tree_union[0]);
			free(first_tree_union);
			return result;
		} else {
			return false;
		}
	} else if (second->type == hol_term_type::ANY) {
		for (unsigned int i = 0; i < second->any.excluded_tree_count; i++)
			if (has_intersection<BuiltInPredicates>(first, second->any.excluded_trees[i])) return false;
		if (second->any.included == nullptr || is_subset<BuiltInPredicates>(first, second->any.included)) return true;

		hol_term included_any;
		included_any.any.included = second->any.included;
		included_any.any.included->reference_count++;
		switch (first->type) {
		case hol_term_type::NOT:
			return is_subset<BuiltInPredicates>(first->unary.operand, &included_any);
		case hol_term_type::UNARY_APPLICATION:
		case hol_term_type::IF_THEN:
		case hol_term_type::EQUALS:
			return is_subset<BuiltInPredicates>(first->binary.left, &included_any)
				|| is_subset<BuiltInPredicates>(first->binary.right, &included_any);
		case hol_term_type::BINARY_APPLICATION:
			return is_subset<BuiltInPredicates>(first->ternary.first, &included_any)
				|| is_subset<BuiltInPredicates>(first->ternary.second, &included_any)
				|| is_subset<BuiltInPredicates>(first->ternary.third, &included_any);
		case hol_term_type::IFF:
		case hol_term_type::AND:
		case hol_term_type::OR:
			for (unsigned int i = 0; i < first->array.length; i++)
				if (is_subset<BuiltInPredicates>(first->array.operands[i], &included_any)) return true;
			return false;
		case hol_term_type::ANY_ARRAY:
			for (unsigned int i = 0; i < first->any_array.left.length; i++)
				if (is_subset<BuiltInPredicates>(first->any_array.left.operands[i], &included_any)) return true;
			for (unsigned int i = 0; i < first->any_array.right.length; i++)
				if (is_subset<BuiltInPredicates>(first->any_array.right.operands[i], &included_any)) return true;
			for (unsigned int i = 0; i < first->any_array.any.length; i++)
				if (is_subset<BuiltInPredicates>(first->any_array.any.operands[i], &included_any)) return true;
			return false;
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			return is_subset<BuiltInPredicates>(first->quantifier.operand, &included_any);
		case hol_term_type::ANY_QUANTIFIER:
			return is_subset<BuiltInPredicates>(first->any_quantifier.operand, &included_any);
		case hol_term_type::INTEGER:
		case hol_term_type::STRING:
		case hol_term_type::UINT_LIST:
		case hol_term_type::CONSTANT:
		case hol_term_type::VARIABLE:
		case hol_term_type::PARAMETER:
		case hol_term_type::ANY_CONSTANT:
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
			return false;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
			/* we already handle this before the switch statement */
			break;
		}
		fprintf(stderr, "is_subset ERROR: Unrecognized hol_term_type.\n");
		return false;
	} else if (second->type == hol_term_type::ANY_RIGHT) {
		for (unsigned int i = 0; i < second->any.excluded_tree_count; i++)
			if (has_intersection<BuiltInPredicates>(first, second->any.excluded_trees[i])) return false;
		if (second->any.included == nullptr || is_subset<BuiltInPredicates>(first, second->any.included)) return true;

		hol_term included_any;
		included_any.type = hol_term_type::ANY_RIGHT;
		included_any.any.included = second->any.included;
		included_any.any.included->reference_count++;
		switch (first->type) {
		case hol_term_type::NOT:
			return is_subset<BuiltInPredicates>(first->unary.operand, &included_any);
		case hol_term_type::UNARY_APPLICATION:
		case hol_term_type::IF_THEN:
		case hol_term_type::EQUALS:
			return is_subset<BuiltInPredicates>(first->binary.right, &included_any);
		case hol_term_type::BINARY_APPLICATION:
			return is_subset<BuiltInPredicates>(first->ternary.third, &included_any);
		case hol_term_type::IFF:
		case hol_term_type::AND:
		case hol_term_type::OR:
			return is_subset<BuiltInPredicates>(first->array.operands[first->array.length - 1], &included_any);
		case hol_term_type::ANY_ARRAY:
			if (first->any_array.right.length == 0)
				return is_subset<BuiltInPredicates>(first->any_array.all, &included_any);
			else return is_subset<BuiltInPredicates>(first->any_array.right.operands[first->any_array.right.length - 1], &included_any);
		case hol_term_type::FOR_ALL:
		case hol_term_type::EXISTS:
		case hol_term_type::LAMBDA:
			return is_subset<BuiltInPredicates>(first->quantifier.operand, &included_any);
		case hol_term_type::ANY_QUANTIFIER:
			return is_subset<BuiltInPredicates>(first->any_quantifier.operand, &included_any);
		case hol_term_type::INTEGER:
		case hol_term_type::STRING:
		case hol_term_type::UINT_LIST:
		case hol_term_type::CONSTANT:
		case hol_term_type::VARIABLE:
		case hol_term_type::PARAMETER:
		case hol_term_type::ANY_CONSTANT:
		case hol_term_type::TRUE:
		case hol_term_type::FALSE:
			return false;
		case hol_term_type::ANY:
		case hol_term_type::ANY_RIGHT:
			/* we already handle this before the switch statement */
			break;
		}
		fprintf(stderr, "is_subset ERROR: Unrecognized hol_term_type.\n");
		return false;
	} else if (first->type == hol_term_type::ANY_ARRAY) {
		if (second->type == hol_term_type::ANY_ARRAY) {
			if ((second->any_array.oper != hol_term_type::ANY_ARRAY && first->any_array.oper != second->any_array.oper)
			 || second->any_array.left.length > first->any_array.left.length
			 || second->any_array.right.length > first->any_array.right.length
			 || second->any_array.any.length > first->any_array.any.length
			 || !is_subset<BuiltInPredicates>(first->any_array.all, second->any_array.all))
				return false;

			for (unsigned int i = 0; i < second->any_array.left.length; i++)
				if (!is_subset<BuiltInPredicates>(first->any_array.left.operands[i], second->any_array.left.operands[i])) return false;
			for (unsigned int i = 0; i < second->any_array.right.length; i++)
				if (!is_subset<BuiltInPredicates>(first->any_array.right.operands[first->any_array.right.length - i - 1], second->any_array.right.operands[second->any_array.right.length - i - 1])) return false;

			bool is_any_subset = false;
			for (unsigned int i = 0; i < first->any_array.any.length - second->any_array.any.length; i++) {
				bool is_any_subset_i = true;
				for (unsigned int j = 0; j < second->any_array.any.length; j++) {
					if (!is_subset<BuiltInPredicates>(first->any_array.any.operands[i + j], second->any_array.any.operands[j])) {
						is_any_subset_i = false;
						break;
					}
				}
				if (is_any_subset_i) {
					is_any_subset = true;
					break;
				}
			}
			return is_any_subset;
		} else {
			/* we already considered the case where `second` has type `ANY` */
			return false;
		}
	} else if (second->type == hol_term_type::ANY_ARRAY) {
		if (first->type == hol_term_type::AND || first->type == hol_term_type::OR || first->type == hol_term_type::IFF) {
			if ((second->any_array.oper != hol_term_type::ANY_ARRAY && second->any_array.oper != first->type)
			 || first->array.length < second->any_array.left.length
			 || first->array.length < second->any_array.right.length
			 || first->array.length < second->any_array.any.length)
				return false;

			for (unsigned int i = 0; i < first->array.length; i++)
				if (!is_subset<BuiltInPredicates>(first->array.operands[i], second->any_array.all)) return false;
			for (unsigned int i = 0; i < second->any_array.left.length; i++)
				if (!is_subset<BuiltInPredicates>(first->array.operands[i], second->any_array.left.operands[i])) return false;
			for (unsigned int i = 0; i < second->any_array.right.length; i++)
				if (!is_subset<BuiltInPredicates>(first->array.operands[first->array.length - i - 1], second->any_array.right.operands[second->any_array.right.length - i - 1])) return false;

			bool found_any = false;
			for (unsigned int i = 0; i < first->array.length - second->any_array.any.length; i++) {
				bool found_any_i = true;
				for (unsigned int j = 0; j < second->any_array.any.length; j++) {
					if (!is_subset<BuiltInPredicates>(first->array.operands[i + j], second->any_array.any.operands[j])) {
						found_any_i = false;
						break;
					}
				}
				if (found_any_i) {
					found_any = true;
					break;
				}
			}
			return found_any;
		}

		/* check for the case that the length of `second` is 1 and `first` is the only item */
		if (second->any_array.any.length > 1 || second->any_array.left.length > 1 || second->any_array.right.length > 1)
			return false;
		if (second->any_array.any.length == 1 && !is_subset<BuiltInPredicates>(first, second->any_array.any.operands[0]))
			return false;
		if (second->any_array.left.length == 1 && !is_subset<BuiltInPredicates>(first, second->any_array.left.operands[0]))
			return false;
		if (second->any_array.right.length == 1 && !is_subset<BuiltInPredicates>(first, second->any_array.right.operands[0]))
			return false;
		return is_subset<BuiltInPredicates>(first, second->any_array.all);
	} else if (first->type == hol_term_type::ANY_CONSTANT) {
		if (second->type != hol_term_type::ANY_CONSTANT)
			return false;
		return is_subset(first->any_constant.constants, first->any_constant.length, second->any_constant.constants, second->any_constant.length);
	} else if (second->type == hol_term_type::ANY_CONSTANT) {
		if (first->type != hol_term_type::CONSTANT)
			return false;
		return index_of(first->constant, second->any_constant.constants, second->any_constant.length) < second->any_constant.length;
	} else if (first->type == hol_term_type::ANY_QUANTIFIER) {
		if (second->type == hol_term_type::ANY_QUANTIFIER) {
			return is_subset(first->any_quantifier.quantifier, second->any_quantifier.quantifier)
				&& is_subset<BuiltInPredicates>(first->any_quantifier.operand, second->any_quantifier.operand);
		} else {
			return false;
		}
	} else if (second->type == hol_term_type::ANY_QUANTIFIER) {
		if (first->type == hol_term_type::FOR_ALL || first->type == hol_term_type::EXISTS || first->type == hol_term_type::LAMBDA) {
			return is_subset((hol_quantifier_type) first->type, second->any_quantifier.quantifier)
				&& is_subset<BuiltInPredicates>(first->quantifier.operand, second->any_quantifier.operand);
		} else {
			return false;
		}
	} else if (first->type != second->type) {
		return false;
	}

	switch (first->type) {
	case hol_term_type::NOT:
		return is_subset<BuiltInPredicates>(first->unary.operand, second->unary.operand);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		return is_subset<BuiltInPredicates>(first->binary.left, second->binary.left)
			&& is_subset<BuiltInPredicates>(first->binary.right, second->binary.right);
	case hol_term_type::BINARY_APPLICATION:
		return is_subset<BuiltInPredicates>(first->ternary.first, second->ternary.first)
			&& is_subset<BuiltInPredicates>(first->ternary.second, second->ternary.second)
			&& is_subset<BuiltInPredicates>(first->ternary.third, second->ternary.third);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		if (first->array.length != second->array.length) return false;
		for (unsigned int i = 0; i < first->array.length; i++)
			if (!is_subset<BuiltInPredicates>(first->array.operands[i], second->array.operands[i])) return false;
		return true;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return (first->quantifier.variable == second->quantifier.variable)
			&& is_subset<BuiltInPredicates>(first->quantifier.operand, second->quantifier.operand);
	case hol_term_type::INTEGER:
		return first->integer == second->integer;
	case hol_term_type::STRING:
		return first->str == second->str;
	case hol_term_type::UINT_LIST:
		return first->uint_list == second->uint_list;
	case hol_term_type::CONSTANT:
		return first->constant == second->constant;
	case hol_term_type::VARIABLE:
		return first->variable == second->variable;
	case hol_term_type::PARAMETER:
		return first->parameter == second->parameter;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		/* we already handle this before the switch statement */
		break;
	}
	fprintf(stderr, "is_subset ERROR: Unrecognized hol_term_type.\n");
	return false;
}

/**
 * Returns whether `first` is a subset of the union of `second_trees[i]`.
 */
template<typename BuiltInPredicates>
bool is_subset(hol_term* first, hol_term** second_trees, unsigned int second_tree_count)
{
	if (first->type == hol_term_type::ANY) {
		unsigned int first_tree_union_length = first->any.excluded_tree_count;
		unsigned int second_tree_union_length = second_tree_count;
		hol_term** first_tree_union = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * first_tree_union_length));
		hol_term** second_tree_union = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * second_tree_union_length));
		if (first_tree_union == nullptr || second_tree_union == nullptr) {
			if (first_tree_union != nullptr) free(first_tree_union);
			fprintf(stderr, "is_subset ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < first_tree_union_length; i++) {
			first_tree_union[i] = first->any.excluded_trees[i];
			first_tree_union[i]->reference_count++;
		} for (unsigned int i = 0; i < second_tree_union_length; i++) {
			second_tree_union[i] = second_trees[i];
			second_tree_union[i]->reference_count++;
		}
		if (!reduce_union<BuiltInPredicates>(first_tree_union, first_tree_union_length, second_tree_union, second_tree_union_length)) {
			fprintf(stderr, "is_subset ERROR: `reduce_union` failed.\n");
			for (unsigned int i = 0; i < first_tree_union_length; i++) {
				free(*first_tree_union[i]); if (first_tree_union[i]->reference_count == 0) free(first_tree_union[i]);
			} for (unsigned int i = 0; i < second_tree_union_length; i++) {
				free(*second_tree_union[i]); if (second_tree_union[i]->reference_count == 0) free(second_tree_union[i]);
			}
			free(first_tree_union); free(second_tree_union);
			return false;
		}

		for (unsigned int i = 0; i < second_tree_union_length; i++) {
			if (is_subset<BuiltInPredicates>(first, second_tree_union[i])) {
				for (unsigned int i = 0; i < first_tree_union_length; i++) {
					free(*first_tree_union[i]); if (first_tree_union[i]->reference_count == 0) free(first_tree_union[i]);
				} for (unsigned int i = 0; i < second_tree_union_length; i++) {
					free(*second_tree_union[i]); if (second_tree_union[i]->reference_count == 0) free(second_tree_union[i]);
				}
				free(first_tree_union); free(second_tree_union);
				return true;	
			}
		}

		for (unsigned int i = 0; i < first_tree_union_length; i++) {
			free(*first_tree_union[i]); if (first_tree_union[i]->reference_count == 0) free(first_tree_union[i]);
		} for (unsigned int i = 0; i < second_tree_union_length; i++) {
			free(*second_tree_union[i]); if (second_tree_union[i]->reference_count == 0) free(second_tree_union[i]);
		}
		free(first_tree_union); free(second_tree_union);
		return false;

	} else {
		/* NOTE: we assume `second_trees` is already in reduced form */
		for (unsigned int i = 0; i < second_tree_count; i++)
			if (is_subset<BuiltInPredicates>(first, second_trees[i])) return true;
		return false;
	}
}

template<typename T, typename Function>
bool apply_to_cartesian_product(array<T>* sets,
		unsigned int set_count, Function function)
{
	if (set_count == 0)
		return function(nullptr);
	unsigned int* index_array = (unsigned int*) calloc(set_count, sizeof(unsigned int));
	while (true) {
		if (!function(index_array)) {
			free(index_array);
			return false;
		}

		/* increment `index_array` */
		bool remaining = false;
		for (unsigned int j = set_count; j > 0; j--) {
			index_array[j - 1]++;
			if (index_array[j - 1] == sets[j - 1].length) {
				index_array[j - 1] = 0;
			} else {
				remaining = true;
				break;
			}
		}
		if (!remaining) break;
	}
	free(index_array);
	return true;
}

template<typename BuiltInPredicates, typename Function>
inline bool subtract_any_with_array(hol_array_term& array_term, hol_term* superset, hol_term* second, Function function)
{
	if (array_term.length == 0)
		return function(nullptr, nullptr);
	array<hol_term*>* difference_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * array_term.length);
	array<hol_term*>* intersection_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * array_term.length);
	if (difference_array == nullptr || intersection_array == nullptr) {
		fprintf(stderr, "subtract_any ERROR: Out of memory.\n");
		if (difference_array != nullptr) free(difference_array);
		if (intersection_array != nullptr) free(intersection_array);
		return false;
	}
	for (unsigned int i = 0; i < array_term.length; i++) {
		if (!array_init(difference_array[i], 8)) {
			for (unsigned int j = 0; j < i; j++)
				free(intersection_array[i]);
			free(intersection_array);
			for (unsigned int j = 0; j < i; j++) {
				for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(difference_array[j]);
			}
			free(difference_array);
			return false;
		} else if (!array_init(intersection_array[i], 8)) {
			for (unsigned int j = 0; j < i; j++)
				free(intersection_array[i]);
			free(intersection_array);
			for (unsigned int j = 0; j <= i; j++) {
				for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(difference_array[j]);
			}
			free(difference_array);
			return false;
		}
		subtract_any<BuiltInPredicates>(difference_array[i], array_term.operands[i], second);
		for (unsigned int j = 0; superset != nullptr && j < difference_array[i].length; j++) {
			if (!is_subset<BuiltInPredicates>(difference_array[i][j], superset)) {
				free(*difference_array[i][j]);
				if (difference_array[i][j]->reference_count == 0)
					free(difference_array[i][j]);
				difference_array[i].remove(j--);
			}
		}
		if (difference_array[i].length == 0) {
			for (unsigned int j = 0; j <= i; j++) {
				for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(difference_array[j]); free(intersection_array[j]);
			}
			free(difference_array); free(intersection_array);
			return true;
		}
	}
	apply_to_cartesian_product(difference_array, array_term.length, [function,&array_term,&difference_array,&intersection_array](const unsigned int* index_array)
	{
		for (unsigned int i = 0; i < array_term.length; i++) {
			intersect<BuiltInPredicates, true>(intersection_array[i], array_term.operands[i], difference_array[i][index_array[i]]);
			if (intersection_array[i].length == 0) {
				for (unsigned int i = 0; i < array_term.length; i++) {
					for (hol_term* term : intersection_array[i]) { free(*term); if (term->reference_count == 0) free(term); }
					intersection_array[i].clear();
				}
				return true;
			}
		}
		bool result = apply_to_cartesian_product(intersection_array, array_term.length, [function,&intersection_array](const unsigned int* intersection_index_array) { return function(intersection_array, intersection_index_array); });
		for (unsigned int i = 0; i < array_term.length; i++) {
			for (hol_term* term : intersection_array[i]) { free(*term); if (term->reference_count == 0) free(term); }
			intersection_array[i].clear();
		}
		return result;
	});
	for (unsigned int i = 0; i < array_term.length; i++) {
		for (hol_term* term : difference_array[i]) { free(*term); if (term->reference_count == 0) free(term); }
		free(difference_array[i]); free(intersection_array[i]);
	}
	free(difference_array); free(intersection_array);
	return true;
}

template<typename BuiltInPredicates>
inline bool subtract_any_with_any(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
	unsigned int excluded_tree_count = 0;
	bool same_as_first = true;
	bool same_as_second = true;
	hol_term** excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * (first->any.excluded_tree_count + 1));
	for (unsigned int i = 0; i < first->any.excluded_tree_count; i++) {
		bool irreducible = true;
		if (is_subset<BuiltInPredicates>(first->any.excluded_trees[i], second))
			irreducible = false;
		if (irreducible) {
			excluded_trees[excluded_tree_count++] = first->any.excluded_trees[i];
			same_as_second = false;
		} else {
			same_as_first = true;
		}
	}
	unsigned int old_excluded_tree_count = excluded_tree_count;

	bool irreducible = true;
	for (unsigned int j = 0; irreducible && j < old_excluded_tree_count; j++) {
		if (is_subset<BuiltInPredicates>(second, excluded_trees[j])) {
			irreducible = false;
			break;
		}
	}
	if (irreducible) {
		excluded_trees[excluded_tree_count++] = second;
		same_as_first = false;
	} else {
		same_as_second = false;
	}

	/* reduce the `excluded_trees` union */
	for (unsigned int i = 0; i < excluded_tree_count; i++)
		excluded_trees[i]->reference_count++;
	if (!same_as_first && !same_as_second)
		reduce_union<BuiltInPredicates>(excluded_trees, excluded_tree_count);

	if (same_as_first) {
		for (unsigned int i = 0; i < excluded_tree_count; i++) {
			free(*excluded_trees[i]);
			if (excluded_trees[i]->reference_count == 0)
				free(excluded_trees[i]);
		}
		free(excluded_trees);
		if (!dst.ensure_capacity(dst.length + 1)) return false;
		dst[dst.length] = first;
		dst[dst.length++]->reference_count++;
		return true;
	}

	/* check if the set difference is empty */
	for (unsigned int i = 0; i < excluded_tree_count; i++) {
		if (index_of(excluded_trees[i], first->any.excluded_trees, first->any.excluded_tree_count) < first->any.excluded_tree_count)
			continue;
#if !defined(NDEBUG)
		if (*excluded_trees[i] == HOL_ANY)
			fprintf(stderr, "subtract_any WARNING: `excluded_trees` contains the set of all logical forms.\n");
#endif
		if (first->any.included != nullptr && is_subset<BuiltInPredicates>(first->any.included, excluded_trees[i])) {
			/* the difference is empty */
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]);
				if (excluded_trees[i]->reference_count == 0)
					free(excluded_trees[i]);
			}
			free(excluded_trees);
			return false;
		}
	}

	if (!dst.ensure_capacity(dst.length + 1)
	 || !new_hol_term(dst[dst.length]))
	{
		for (unsigned int i = 0; i < excluded_tree_count; i++) {
			free(*excluded_trees[i]);
			if (excluded_trees[i]->reference_count == 0)
				free(excluded_trees[i]);
		}
		free(excluded_trees);
		return false;
	}
	dst[dst.length]->type = first->type;
	dst[dst.length]->reference_count = 1;
	dst[dst.length]->any.included = first->any.included;
	if (first->any.included != nullptr)
		dst[dst.length]->any.included->reference_count++;
	dst[dst.length]->any.excluded_trees = excluded_trees;
	dst[dst.length]->any.excluded_tree_count = excluded_tree_count;
	dst.length++;
	return true;
}

template<typename BuiltInPredicates>
bool subtract_any(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (second->type != hol_term_type::ANY || second->any.excluded_tree_count != 0 || second->any.included == nullptr) {
		fprintf(stderr, "subtract_any ERROR: Expected `second` to have type `ANY` and have no excluded sets.\n");
		return false;
	}
#endif
	size_t old_dst_length = dst.length;
	if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) {
		return subtract_any_with_any<BuiltInPredicates>(dst, first, second);

	} else if (first->type == hol_term_type::ANY_ARRAY) {
		array<hol_term*> differences(8);
		if (second->any.included->type == hol_term_type::ANY_ARRAY || second->any.included->type == hol_term_type::AND || second->any.included->type == hol_term_type::OR || second->any.included->type == hol_term_type::IFF) {
			subtract<BuiltInPredicates>(differences, first, second->any.included);
		} else {
			differences[differences.length++] = first;
			first->reference_count++;
		}

		array<hol_term*> first_differences(8);
		subtract_any<BuiltInPredicates>(first_differences, first->any_array.all, second);
		for (hol_term* root : differences) {
			for (hol_term* all : first_differences) {
				bool result = subtract_any_with_array<BuiltInPredicates>(root->any_array.left, all, second, [root,first,all,second,&dst](array<hol_term*>* left_array, const unsigned int* left_index_array) {
					return subtract_any_with_array<BuiltInPredicates>(root->any_array.any, all, second, [root,first,all,second,left_array,left_index_array,&dst](array<hol_term*>* any_array, const unsigned int* any_index_array) {
						return subtract_any_with_array<BuiltInPredicates>(root->any_array.right, all, second, [root,first,all,left_array,left_index_array,any_array,any_index_array,&dst](array<hol_term*>* right_array, const unsigned int* right_index_array) {
							if (!dst.ensure_capacity(dst.length + 1))
								return false;

							bool same_as_root = (all == root->any_array.all);
							for (unsigned int i = 0; same_as_root && i < root->any_array.left.length; i++)
								if (root->any_array.left.operands[i] != left_array[i][left_index_array[i]]) same_as_root = false;
							for (unsigned int i = 0; same_as_root && i < root->any_array.right.length; i++)
								if (root->any_array.right.operands[i] != right_array[i][right_index_array[i]]) same_as_root = false;
							for (unsigned int i = 0; same_as_root && i < root->any_array.any.length; i++)
								if (root->any_array.any.operands[i] != any_array[i][any_index_array[i]]) same_as_root = false;

							if (same_as_root) {
								dst[dst.length++] = root;
								root->reference_count++;
								return true;
							}

							if (!new_hol_term(dst[dst.length])) return false;
							dst[dst.length]->type = first->type;
							dst[dst.length]->reference_count = 1;
							dst[dst.length]->any_array.oper = root->any_array.oper;
							dst[dst.length]->any_array.all = all;
							all->reference_count++;
							dst[dst.length]->any_array.left.length = first->any_array.left.length;
							if (first->any_array.left.length == 0) {
								dst[dst.length]->any_array.left.operands = nullptr;
							} else {
								dst[dst.length]->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * first->any_array.left.length);
								if (dst[dst.length]->any_array.left.operands == nullptr) {
									fprintf(stderr, "subtract_any ERROR: Out of memory.\n");
									free(dst[dst.length]); return false;
								}
								for (unsigned int i = 0; i < dst[dst.length]->any_array.left.length; i++) {
									dst[dst.length]->any_array.left.operands[i] = left_array[i][left_index_array[i]];
									dst[dst.length]->any_array.left.operands[i]->reference_count++;
								}
							}
							dst[dst.length]->any_array.any.length = first->any_array.any.length;
							if (first->any_array.any.length == 0) {
								dst[dst.length]->any_array.any.operands = nullptr;
							} else {
								dst[dst.length]->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * first->any_array.any.length);
								if (dst[dst.length]->any_array.any.operands == nullptr) {
									fprintf(stderr, "subtract_any ERROR: Out of memory.\n");
									for (unsigned int j = 0; j < dst[dst.length]->any_array.left.length; j++) {
										free(*dst[dst.length]->any_array.left.operands[j]);
										if (dst[dst.length]->any_array.left.operands[j]->reference_count == 0)
											free(dst[dst.length]->any_array.left.operands[j]);
									}
									free(dst[dst.length]->any_array.left.operands);
									free(dst[dst.length]); return false;
								}
								for (unsigned int i = 0; i < dst[dst.length]->any_array.any.length; i++) {
									dst[dst.length]->any_array.any.operands[i] = any_array[i][any_index_array[i]];
									dst[dst.length]->any_array.any.operands[i]->reference_count++;
								}
							}
							dst[dst.length]->any_array.right.length = first->any_array.right.length;
							if (first->any_array.right.length == 0) {
								dst[dst.length]->any_array.right.operands = nullptr;
							} else {
								dst[dst.length]->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * first->any_array.right.length);
								if (dst[dst.length]->any_array.right.operands == nullptr) {
									fprintf(stderr, "subtract_any ERROR: Out of memory.\n");
									for (unsigned int j = 0; j < dst[dst.length]->any_array.left.length; j++) {
										free(*dst[dst.length]->any_array.left.operands[j]);
										if (dst[dst.length]->any_array.left.operands[j]->reference_count == 0)
											free(dst[dst.length]->any_array.left.operands[j]);
									}
									free(dst[dst.length]->any_array.left.operands);
									for (unsigned int j = 0; j < dst[dst.length]->any_array.any.length; j++) {
										free(*dst[dst.length]->any_array.any.operands[j]);
										if (dst[dst.length]->any_array.any.operands[j]->reference_count == 0)
											free(dst[dst.length]->any_array.any.operands[j]);
									}
									free(dst[dst.length]->any_array.any.operands);
									free(dst[dst.length]); return false;
								}
								for (unsigned int i = 0; i < dst[dst.length]->any_array.right.length; i++) {
									dst[dst.length]->any_array.right.operands[i] = right_array[i][right_index_array[i]];
									dst[dst.length]->any_array.right.operands[i]->reference_count++;
								}
							}
							dst.length++;
							return true;
						});
					});
				});
				if (!result) return false;
			}
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	}

	array<hol_term*> differences(4);
	if (!subtract<BuiltInPredicates>(differences, first, second->any.included) || differences.length == 0)
		return false;
#if !defined(NDEBUG)
	for (const hol_term* diff : differences)
		if (diff->type != first->type)
			fprintf(stderr, "subtract_any: The root of the difference changed.\n");
#endif

	array<hol_term*> first_differences(8);
	array<hol_term*> second_differences(8);
	array<hol_term*> third_differences(8);
	array<hol_term*> any_differences(8);
	switch (first->type) {
	case hol_term_type::NOT:
		subtract_any<BuiltInPredicates>(first_differences, first->unary.operand, second);
		for (hol_term* root : differences) {
			for (hol_term* first_child : first_differences) {
				array<hol_term*> intersection(8);
				intersect<BuiltInPredicates, true>(intersection, root->unary.operand, first_child);
				if (!dst.ensure_capacity(dst.length + intersection.length)) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				} else if (intersection.length == 1 && intersection[0] == root->unary.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
				for (hol_term* term : intersection) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length]->type = first->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->unary.operand = term;
					term->reference_count++;
					dst.length++;
				}
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		subtract_any<BuiltInPredicates>(first_differences, first->binary.left, second);
		subtract_any<BuiltInPredicates>(second_differences, first->binary.right, second);
		for (hol_term* root : differences) {
			for (hol_term* first_child : first_differences) {
				for (hol_term* second_child : second_differences) {
					array<hol_term*> first_intersection(8);
					array<hol_term*> second_intersection(8);
					intersect<BuiltInPredicates, true>(first_intersection, root->binary.left, first_child);
					intersect<BuiltInPredicates, true>(second_intersection, root->binary.right, second_child);
					unsigned int intersection_count = first_intersection.length * second_intersection.length;
					if (!dst.ensure_capacity(dst.length + intersection_count)) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					} else if (intersection_count == 1 && first_intersection[0] == root->binary.left && second_intersection[0] == root->binary.right) {
						dst[dst.length++] = root;
						root->reference_count++;
						for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						continue;
					}
					for (hol_term* first_term : first_intersection) {
						for (hol_term* second_term : second_intersection) {
							if (!new_hol_term(dst[dst.length])) {
								for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
								return false;
							}
							dst[dst.length]->type = first->type;
							dst[dst.length]->reference_count = 1;
							dst[dst.length]->binary.left = first_term;
							dst[dst.length]->binary.right = second_term;
							first_term->reference_count++;
							second_term->reference_count++;
							dst.length++;
						}
					}
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				}
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::BINARY_APPLICATION:
		subtract_any<BuiltInPredicates>(first_differences, first->ternary.first, second);
		subtract_any<BuiltInPredicates>(second_differences, first->ternary.second, second);
		subtract_any<BuiltInPredicates>(third_differences, first->ternary.third, second);
		for (hol_term* root : differences) {
			for (hol_term* first_child : first_differences) {
				for (hol_term* second_child : second_differences) {
					for (hol_term* third_child : third_differences) {
						array<hol_term*> first_intersection(8);
						array<hol_term*> second_intersection(8);
						array<hol_term*> third_intersection(8);
						intersect<BuiltInPredicates, true>(first_intersection, root->ternary.first, first_child);
						intersect<BuiltInPredicates, true>(second_intersection, root->ternary.second, second_child);
						intersect<BuiltInPredicates, true>(third_intersection, root->ternary.third, third_child);
						unsigned int intersection_count = first_intersection.length * second_intersection.length * third_intersection.length;
						if (!dst.ensure_capacity(dst.length + intersection_count)) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							return false;
						} else if (intersection_count == 1 && first_intersection[0] == root->ternary.first && second_intersection[0] == root->ternary.second && third_intersection[0] == root->ternary.third) {
							dst[dst.length++] = root;
							root->reference_count++;
							for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							continue;
						}
						for (hol_term* first_term : first_intersection) {
							for (hol_term* second_term : second_intersection) {
								for (hol_term* third_term : third_intersection) {
									if (!new_hol_term(dst[dst.length])) {
										for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
										return false;
									}
									dst[dst.length]->type = first->type;
									dst[dst.length]->reference_count = 1;
									dst[dst.length]->ternary.first = first_term;
									dst[dst.length]->ternary.second = second_term;
									dst[dst.length]->ternary.third = third_term;
									first_term->reference_count++;
									second_term->reference_count++;
									third_term->reference_count++;
									dst.length++;
								}
							}
						}
						for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					}
				}
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		for (hol_term* root : differences) {
			subtract_any_with_array<BuiltInPredicates>(root->array, nullptr, second, [root,first,&dst](array<hol_term*>* intersection_array, const unsigned int* intersection_index_array) {
				if (!dst.ensure_capacity(dst.length + 1))
					return false;
				bool same_as_root = true;
				for (unsigned int i = 0; same_as_root && i < root->array.length; i++)
					if (intersection_array[i][intersection_index_array[i]] != root->array.operands[i]) same_as_root = false;

				if (same_as_root) {
					dst[dst.length++] = root;
					root->reference_count++;
					return true;
				}

				if (!new_hol_term(dst[dst.length])) return false;
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->array.length = first->array.length;
				dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
				if (dst[dst.length]->array.operands == nullptr) {
					fprintf(stderr, "subtract_any ERROR: Out of memory.\n");
					free(dst[dst.length]); return false;
				}
				for (unsigned int i = 0; i < first->array.length; i++) {
					dst[dst.length]->array.operands[i] = intersection_array[i][intersection_index_array[i]];
					dst[dst.length]->array.operands[i]->reference_count++;
				}
				dst.length++;
				return true;
			});
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		subtract_any<BuiltInPredicates>(first_differences, first->quantifier.operand, second);
		for (hol_term* root : differences) {
			for (hol_term* first_child : first_differences) {
				array<hol_term*> intersection(8);
				intersect<BuiltInPredicates, true>(intersection, root->quantifier.operand, first_child);
				if (!dst.ensure_capacity(dst.length + intersection.length)) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				} else if (intersection.length == 1 && intersection[0] == root->quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
				for (hol_term* term : intersection) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length]->type = first->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->quantifier.variable = first->quantifier.variable;
					dst[dst.length]->quantifier.operand = term;
					term->reference_count++;
					dst.length++;
				}
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY_CONSTANT:
		if (!dst.append(differences.data, differences.length)) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		return (dst.length > old_dst_length);
	case hol_term_type::ANY_QUANTIFIER:
		subtract_any<BuiltInPredicates>(first_differences, first->any_quantifier.operand, second);
		for (hol_term* root : differences) {
			for (hol_term* first_child : first_differences) {
				array<hol_term*> intersection(8);
				intersect<BuiltInPredicates, true>(intersection, root->any_quantifier.operand, first_child);
				if (!dst.ensure_capacity(dst.length + intersection.length)) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				} else if (intersection.length == 1 && intersection[0] == root->any_quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
				for (hol_term* term : intersection) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length]->type = root->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->any_quantifier.quantifier = root->any_quantifier.quantifier;
					dst[dst.length]->any_quantifier.operand = term;
					term->reference_count++;
					dst.length++;
				}
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
		break; /* we already handle this case before the switch statement */
	}
	fprintf(stderr, "subtract_any ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<typename BuiltInPredicates>
bool subtract_any_right(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (second->type != hol_term_type::ANY_RIGHT || second->any.included == nullptr) {
		fprintf(stderr, "subtract_any_right ERROR: Expected `second` to have type `ANY_ARRAY`.\n");
		return false;
	}
#endif
	if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) {
		return subtract_any_with_any<BuiltInPredicates>(dst, first, second);

	} else if (first->type == hol_term_type::ANY_ARRAY) {
		array<hol_term*> differences(8);
		if (second->any.included->type == hol_term_type::ANY_ARRAY || second->any.included->type == hol_term_type::AND || second->any.included->type == hol_term_type::OR || second->any.included->type == hol_term_type::IFF) {
			subtract<BuiltInPredicates>(differences, first, second->any.included);
		} else {
			differences[differences.length++] = first;
			first->reference_count++;
		}

		array<hol_term*> first_differences(8);
		hol_term* right;
		if (first->any_array.right.length == 0)
			right = first->any_array.all;
		else right = first->any_array.right.operands[first->any_array.right.length - 1];
		subtract_any_right<BuiltInPredicates>(first_differences, right, second);
		if (!dst.ensure_capacity(dst.length + differences.length * first_differences.length)) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		size_t old_dst_length = dst.length;
		for (hol_term* root : differences) {
			for (hol_term* new_right : first_differences) {
				if (root == first && new_right == right) {
					dst[dst.length] = first;
					dst[dst.length++]->reference_count++;
					continue;
				}
				dst[dst.length] = hol_term::new_any_array(root->any_array.oper, root->any_array.all,
						make_array_view(root->any_array.any.operands, root->any_array.any.length),
						make_array_view(root->any_array.left.operands, root->any_array.left.length),
						make_appended_array_view(make_array_view(root->any_array.right.operands, max(1u, root->any_array.right.length) - 1), new_right));
				if (dst[dst.length] == nullptr) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->any_array.all->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.any.length; i++)
					dst[dst.length]->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.left.length; i++)
					dst[dst.length]->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.right.length; i++)
					dst[dst.length]->any_array.right.operands[i]->reference_count++;
				dst.length++;
			}
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	}

	array<hol_term*> differences(4);
	if (!subtract<BuiltInPredicates>(differences, first, second->any.included) || differences.length == 0)
		return false;
#if !defined(NDEBUG)
	for (const hol_term* diff : differences)
		if (diff->type != first->type)
			fprintf(stderr, "subtract_any_right: The root of the difference changed.\n");
#endif

	size_t old_dst_length = dst.length;
	switch (first->type) {
	case hol_term_type::NOT:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->unary.operand, second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->unary.operand) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->unary.operand = first_child;
				dst[dst.length]->unary.operand->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->binary.right, second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->binary.right) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->binary.left = root->binary.left;
				dst[dst.length]->binary.right = first_child;
				dst[dst.length]->binary.left->reference_count++;
				dst[dst.length]->binary.right->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::BINARY_APPLICATION:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->ternary.third, second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->ternary.third) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->ternary.first = root->ternary.first;
				dst[dst.length]->ternary.second = root->ternary.second;
				dst[dst.length]->ternary.third = first_child;
				dst[dst.length]->ternary.first->reference_count++;
				dst[dst.length]->ternary.second->reference_count++;
				dst[dst.length]->ternary.third->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->array.operands[root->array.length - 1], second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->array.operands[root->array.length - 1]) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->array.length = root->array.length;
				dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * root->array.length);
				if (dst[dst.length]->array.operands == nullptr) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(dst[dst.length]); return false;
				}
				for (unsigned int i = 0; i + 1 < root->array.length; i++) {
					dst[dst.length]->array.operands[i] = root->array.operands[i];
					dst[dst.length]->array.operands[i]->reference_count++;
				}
				dst[dst.length]->array.operands[root->array.length - 1] = first_child;
				dst[dst.length]->array.operands[root->array.length - 1]->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->quantifier.operand, second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->quantifier.operand) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->quantifier.variable = root->quantifier.variable;
				dst[dst.length]->quantifier.operand = first_child;
				dst[dst.length]->quantifier.operand->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY_CONSTANT:
		return dst.append(differences.data, differences.length) && dst.length > 0;
	case hol_term_type::ANY_QUANTIFIER:
		for (hol_term* root : differences) {
			array<hol_term*> first_differences(8);
			subtract_any_right<BuiltInPredicates>(first_differences, root->any_quantifier.operand, second);
			if (!dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (first_differences.length == 1 && first_differences[0] == root->any_quantifier.operand) {
				dst[dst.length++] = root;
				root->reference_count++;
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				continue;
			}
			for (hol_term* first_child : first_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->any_quantifier.quantifier = root->any_quantifier.quantifier;
				dst[dst.length]->any_quantifier.operand = first_child;
				dst[dst.length]->any_quantifier.operand->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
		break; /* we already handle this case before the switch statement */
	}
	fprintf(stderr, "subtract_any_right ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<typename BuiltInPredicates>
bool subtract(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
	if (first == second) {
		return false;

	} else if (second->type == hol_term_type::ANY) {
		size_t old_dst_length = dst.length;
		if (second->any.included != nullptr) {
			if (second->any.excluded_tree_count == 0) {
				subtract_any<BuiltInPredicates>(dst, first, second);
			} else {
				hol_term* term = hol_term::new_any(second->any.included, nullptr, 0);
				if (term == nullptr) return false;
				second->any.included->reference_count++;
				subtract_any<BuiltInPredicates>(dst, first, term);
				free(*term); if (term->reference_count == 0) free(term);
			}
		}

		for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
			array<hol_term*> intersections(8);
			if (!intersect<BuiltInPredicates, true>(intersections, first, second->any.excluded_trees[i]))
				continue;

			/* subtract `intersections` with `dst` so far, to ensure the sets in `dst` are disjoint */
			for (unsigned int j = 0; j < dst.length; j++) {
				array<hol_term*> new_intersections(8);
				for (hol_term* intersection : intersections)
					subtract<BuiltInPredicates>(new_intersections, intersection, dst[j]);
				for (hol_term* term : intersections) {
					free(*term);
					if (term->reference_count == 0)
						free(term);
				}
				swap(intersections, new_intersections);
			}

			/* add the disjoint sets to `dst` */
			if (!dst.append(intersections.data, intersections.length)) {
				for (hol_term* term : intersections) {
					free(*term);
					if (term->reference_count == 0)
						free(term);
				}
				return false;
			}
		}
		return (dst.length > old_dst_length);

	} else if (first->type == hol_term_type::ANY) {
		unsigned int excluded_tree_count = 0;
		bool same_as_first = true;
		bool same_as_second = true;
		hol_term** excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * (first->any.excluded_tree_count + 1));
		for (unsigned int i = 0; i < first->any.excluded_tree_count; i++) {
			bool irreducible = true;
			if (is_subset<BuiltInPredicates>(first->any.excluded_trees[i], second)) {
				irreducible = false;
				break;
			}
			if (irreducible) {
				excluded_trees[excluded_tree_count++] = first->any.excluded_trees[i];
				same_as_second = false;
			} else {
				same_as_first = false;
			}
		}
		unsigned int old_excluded_tree_count = excluded_tree_count;

		bool irreducible = true;
		for (unsigned int j = 0; j < old_excluded_tree_count; j++) {
			if (is_subset<BuiltInPredicates>(second, excluded_trees[j])) {
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			excluded_trees[excluded_tree_count++] = second;
			same_as_first = false;
		} else {
			same_as_second = false;
		}

		/* reduce the `excluded_trees` union */
		for (unsigned int i = 0; i < excluded_tree_count; i++)
			excluded_trees[i]->reference_count++;
		if (!same_as_first && !same_as_second)
			reduce_union<BuiltInPredicates>(excluded_trees, excluded_tree_count);

		if (same_as_first) {
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]);
				if (excluded_trees[i]->reference_count == 0)
					free(excluded_trees[i]);
			}
			free(excluded_trees);
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length++] = first;
			first->reference_count++;
			return true;
		}

		if (!dst.ensure_capacity(dst.length + 1)
		 || !new_hol_term(dst[dst.length]))
		{
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]);
				if (excluded_trees[i]->reference_count == 0)
					free(excluded_trees[i]);
			}
			free(excluded_trees);
			return false;
		}
		dst[dst.length]->type = hol_term_type::ANY;
		dst[dst.length]->reference_count = 1;
		dst[dst.length]->any.included = first->any.included;
		if (first->any.included != nullptr)
			dst[dst.length]->any.included->reference_count++;
		dst[dst.length]->any.excluded_trees = excluded_trees;
		dst[dst.length]->any.excluded_tree_count = excluded_tree_count;
		dst.length++;
		return true;

	} else if (second->type == hol_term_type::ANY_RIGHT) {
		if (second->any.included != nullptr)
			return subtract_any_right<BuiltInPredicates>(dst, first, second);
		else return false;

	} else if (first->type == hol_term_type::ANY_RIGHT) {
		if (!has_intersection<BuiltInPredicates>(second, first)) {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		} else {
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		}

	} else if (first->type == hol_term_type::ANY_ARRAY) {
		if (second->type == hol_term_type::ANY_ARRAY) {
			if (first->any_array.oper != hol_term_type::ANY_ARRAY && first->any_array.oper != second->any_array.oper) {
				if (!dst.ensure_capacity(dst.length + 1)) return false;
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			}
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		} else if (second->type == hol_term_type::AND || second->type == hol_term_type::OR || second->type == hol_term_type::IFF) {
			if (first->any_array.oper != hol_term_type::ANY_ARRAY && first->any_array.oper != second->type) {
				if (!dst.ensure_capacity(dst.length + 1)) return false;
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			}
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		} else {
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		}

	} else if (second->type == hol_term_type::ANY_ARRAY) {
		if (first->type == hol_term_type::AND || first->type == hol_term_type::OR || first->type == hol_term_type::IFF) {
			if (second->any_array.oper != hol_term_type::ANY_ARRAY && second->any_array.oper != first->type) {
				if (!dst.ensure_capacity(dst.length + 1)) return false;
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			} else if (is_subset<BuiltInPredicates>(first, second)) {
				return false;
			}
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		} else {
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		}

	} else if (first->type == hol_term_type::ANY_CONSTANT) {
		if (second->type == hol_term_type::ANY_CONSTANT || second->type == hol_term_type::CONSTANT) {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			unsigned int* new_constants = (unsigned int*) malloc(sizeof(unsigned int) * first->any_constant.length);
			if (new_constants == nullptr) {
				fprintf(stderr, "subtract ERROR: Insufficient memory for `new_constants`.\n");
				return false;
			}
			unsigned int new_length = 0;
			if (second->type == hol_term_type::ANY_CONSTANT)
				set_subtract(new_constants, new_length, first->any_constant.constants, first->any_constant.length, second->any_constant.constants, second->any_constant.length);
			else set_subtract(new_constants, new_length, first->any_constant.constants, first->any_constant.length, &second->constant, 1);
			if (new_length == first->any_constant.length) {
				free(new_constants);
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			} else if (new_length == 0) {
				free(new_constants);
				return false;
			} else if (new_length == 1) {
				dst[dst.length] = hol_term::new_constant(new_constants[0]);
				free(new_constants);
				if (dst[dst.length] == nullptr)
					return false;
				dst.length++;
				return true;
			} else if (!new_hol_term(dst[dst.length])) {
				free(new_constants);
				return false;
			}

			dst[dst.length]->type = hol_term_type::ANY_CONSTANT;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->any_constant.length = new_length;
			dst[dst.length]->any_constant.constants = new_constants;
			dst.length++;
			return true;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}

	} else if (second->type == hol_term_type::ANY_CONSTANT) {
		if (first->type == hol_term_type::CONSTANT && index_of(first->constant, second->any_constant.constants, second->any_constant.length) < second->any_constant.length) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}

	} else if (first->type == hol_term_type::ANY_QUANTIFIER) {
		hol_term* second_operand;
		hol_quantifier_type second_quantifier;
		if (second->type == hol_term_type::ANY_QUANTIFIER) {
			second_operand = second->any_quantifier.operand;
			second_quantifier = second->any_quantifier.quantifier;
		} else if (second->type == hol_term_type::FOR_ALL || second->type == hol_term_type::EXISTS || second->type == hol_term_type::LAMBDA) {
			fprintf(stderr, "subtract ERROR: Unclosed subtraction.\n");
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}

		hol_quantifier_type quantifier;
		if (!subtract(quantifier, first->any_quantifier.quantifier, second_quantifier))
			return false;
		array<hol_term*> differences(8);
		if (!subtract<BuiltInPredicates>(differences, first->any_quantifier.operand, second_operand))
			return false;
		if (!dst.ensure_capacity(dst.length + differences.length)) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (differences.length == 1 && quantifier == first->any_quantifier.quantifier && differences[0] == first->any_quantifier.operand) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}

		for (hol_term* difference : differences) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = hol_term_type::ANY_QUANTIFIER;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->any_quantifier.quantifier = quantifier;
			dst[dst.length]->any_quantifier.operand = difference;
			difference->reference_count++;
			dst.length++;
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		return true;

	} else if (second->type == hol_term_type::ANY_QUANTIFIER) {
		if ((first->type == hol_term_type::FOR_ALL || first->type == hol_term_type::EXISTS || first->type == hol_term_type::LAMBDA)
		 && is_subset((hol_quantifier_type) first->type, second->any_quantifier.quantifier))
		{
			array<hol_term*> differences(8);
			subtract<BuiltInPredicates>(differences, first->quantifier.operand, second->any_quantifier.operand);
			if (!dst.ensure_capacity(dst.length + differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			} else if (differences.length == 1 && first->quantifier.operand == differences[0]) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			}

			for (hol_term* difference : differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->quantifier.variable = first->quantifier.variable;
				dst[dst.length]->quantifier.operand = difference;
				difference->reference_count++;
				dst.length++;
			}
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}

	} else if (first->type != second->type) {
		if (!dst.ensure_capacity(dst.length + 1)) return false;
		dst[dst.length] = first;
		dst[dst.length++]->reference_count++;
		return true;
	}

	array<hol_term*> first_differences(8);
	array<hol_term*> second_differences(8);
	array<hol_term*> third_differences(8);
	array<hol_term*>* difference_array;
	unsigned int difference_count;
	size_t old_dst_length = dst.length;
	switch (first->type) {
	case hol_term_type::NOT:
		subtract<BuiltInPredicates>(first_differences, first->unary.operand, second->unary.operand);
		if (!dst.ensure_capacity(dst.length + first_differences.length)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_differences.length == 1 && first_differences[0] == first->unary.operand) {
			dst[dst.length++] = first;
			first->reference_count++;
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		}
		for (hol_term* first_child : first_differences) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->unary.operand = first_child;
			first_child->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		subtract<BuiltInPredicates>(first_differences, first->binary.left, second->binary.left);
		if (!dst.ensure_capacity(dst.length + first_differences.length)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_differences.length == 1 && first_differences[0] == first->binary.left) {
			dst[dst.length++] = first;
			first->reference_count++;
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		}
		for (hol_term* first_child : first_differences) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->binary.left = first_child;
			dst[dst.length]->binary.right = first->binary.right;
			dst[dst.length]->binary.left->reference_count++;
			dst[dst.length]->binary.right->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		first_differences.clear();

		intersect<BuiltInPredicates, true>(first_differences, first->binary.left, second->binary.left);
		subtract<BuiltInPredicates>(second_differences, first->binary.right, second->binary.right);
		difference_count = first_differences.length * second_differences.length;
		if (!dst.ensure_capacity(dst.length + difference_count)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		for (hol_term* first_child : first_differences) {
			for (hol_term* second_child : second_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->binary.left = first_child;
				dst[dst.length]->binary.right = second_child;
				dst[dst.length]->binary.left->reference_count++;
				dst[dst.length]->binary.right->reference_count++;
				dst.length++;
			}
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::BINARY_APPLICATION:
		subtract<BuiltInPredicates>(first_differences, first->ternary.first, second->ternary.first);
		if (!dst.ensure_capacity(dst.length + first_differences.length)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_differences.length == 1 && first_differences[0] == first->ternary.first) {
			dst[dst.length++] = first;
			first->reference_count++;
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		}
		for (hol_term* first_child : first_differences) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->ternary.first = first_child;
			dst[dst.length]->ternary.second = first->ternary.second;
			dst[dst.length]->ternary.third = first->ternary.third;
			dst[dst.length]->ternary.first->reference_count++;
			dst[dst.length]->ternary.second->reference_count++;
			dst[dst.length]->ternary.third->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		first_differences.clear();

		intersect<BuiltInPredicates, true>(first_differences, first->ternary.first, second->ternary.first);
		if (first_differences.length == 0)
			return (dst.length > old_dst_length);
		subtract<BuiltInPredicates>(second_differences, first->ternary.second, second->ternary.second);
		difference_count = first_differences.length * second_differences.length;
		if (!dst.ensure_capacity(dst.length + difference_count)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		for (hol_term* first_child : first_differences) {
			for (hol_term* second_child : second_differences) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->ternary.first = first_child;
				dst[dst.length]->ternary.second = second_child;
				dst[dst.length]->ternary.third = first->ternary.third;
				dst[dst.length]->ternary.first->reference_count++;
				dst[dst.length]->ternary.second->reference_count++;
				dst[dst.length]->ternary.third->reference_count++;
				dst.length++;
			}
		}
		for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		second_differences.clear();

		intersect<BuiltInPredicates, true>(second_differences, first->ternary.second, second->ternary.second);
		if (second_differences.length == 0) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return (dst.length > old_dst_length);
		}
		subtract<BuiltInPredicates>(third_differences, first->ternary.third, second->ternary.third);
		difference_count = first_differences.length * second_differences.length;
		if (!dst.ensure_capacity(dst.length + difference_count)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		for (hol_term* first_child : first_differences) {
			for (hol_term* second_child : second_differences) {
				for (hol_term* third_child : third_differences) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length]->type = first->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->ternary.first = first_child;
					dst[dst.length]->ternary.second = second_child;
					dst[dst.length]->ternary.third = third_child;
					dst[dst.length]->ternary.first->reference_count++;
					dst[dst.length]->ternary.second->reference_count++;
					dst[dst.length]->ternary.third->reference_count++;
					dst.length++;
				}
			}
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		if (first->array.length != second->array.length) {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
		difference_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * first->array.length);
		if (difference_array == nullptr) {
			fprintf(stderr, "subtract ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < first->array.length; i++) {
			if (!array_init(difference_array[i], 8)) {
				for (unsigned int j = 0; j < i; j++) {
					for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(difference_array[j]);
				}
				free(difference_array);
				return false;
			}
		}
		difference_count = 1;
		for (unsigned int i = 0; i < first->array.length; i++) {
			if (i > 0) {
				intersect<BuiltInPredicates, true>(difference_array[i - 1], first->array.operands[i - 1], second->array.operands[i - 1]);
				if (difference_array[i - 1].length == 0) {
					for (unsigned int j = 0; j < first->array.length; j++) {
						for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(difference_array[j]);
					}
					free(difference_array);
					return (dst.length > old_dst_length);
				}
				difference_count *= difference_array[i - 1].length;
			}
			subtract<BuiltInPredicates>(difference_array[i], first->array.operands[i], second->array.operands[i]);
			unsigned int count = difference_count * difference_array[i].length;
			if (!dst.ensure_capacity(dst.length + count)) {
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(difference_array[j]);
				}
				free(difference_array);
				return false;
			}
			bool result = (count == 0) || apply_to_cartesian_product(difference_array, i + 1, [&dst,i,first,difference_array](const unsigned int* indices) {
				if (!new_hol_term(dst[dst.length]))
					return false;
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->array.length = first->array.length;
				dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
				if (dst[dst.length]->array.operands == nullptr) {
					fprintf(stderr, "subtract ERROR: Out of memory.\n");
					free(dst[dst.length]); return false;
				}
				for (unsigned int j = 0; j <= i; j++) {
					dst[dst.length]->array.operands[j] = difference_array[j][indices[j]];
					dst[dst.length]->array.operands[j]->reference_count++;
				} for (unsigned int j = i + 1; j < first->array.length; j++) {
					dst[dst.length]->array.operands[j] = first->array.operands[j];
					dst[dst.length]->array.operands[j]->reference_count++;
				}
				dst.length++;
				return true;
			});
			if (!result) {
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(difference_array[j]);
				}
				free(difference_array);
				return false;
			}
			for (hol_term* term : difference_array[i]) { free(*term); if (term->reference_count == 0) free(term); }
			difference_array[i].clear();
		}
		for (unsigned int j = 0; j < first->array.length; j++) {
			for (hol_term* term : difference_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
			free(difference_array[j]);
		}
		free(difference_array);
		return (dst.length > old_dst_length);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (first->quantifier.variable != second->quantifier.variable) {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
		subtract<BuiltInPredicates>(first_differences, first->quantifier.operand, second->quantifier.operand);
		if (!dst.ensure_capacity(dst.length + first_differences.length)) {
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_differences.length == 1 && first_differences[0] == first->quantifier.operand) {
			dst[dst.length++] = first;
			first->reference_count++;
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		}
		for (hol_term* first_child : first_differences) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->quantifier.variable = first->quantifier.variable;
			dst[dst.length]->quantifier.operand = first_child;
			first_child->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	case hol_term_type::INTEGER:
		if (first->integer == second->integer) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::STRING:
		if (first->str == second->str) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::UINT_LIST:
		if (first->uint_list == second->uint_list) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::CONSTANT:
		if (first->constant == second->constant) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::VARIABLE:
		if (first->variable == second->variable) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::PARAMETER:
		if (first->parameter == second->parameter) {
			return false;
		} else {
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		}
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return false;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		break; /* we already handle this case before the switch statement */
	}
	fprintf(stderr, "subtract ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<typename BuiltInPredicates>
inline bool any_is_excluded(
		hol_term* term, hol_term* first, hol_term* second,
		hol_term** excluded_trees, unsigned int excluded_tree_count,
		bool same_trees_as_first, bool same_trees_as_second)
{
	if (term != nullptr && (term->type == hol_term_type::ANY || term->type == hol_term_type::ANY_RIGHT)) {
		unsigned int first_union_length = term->any.excluded_tree_count;
		unsigned int second_union_length = excluded_tree_count;
		hol_term** first_union = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * first_union_length));
		hol_term** second_union = (hol_term**) malloc(max((size_t) 1, sizeof(hol_term*) * second_union_length));
		if (first_union == nullptr || second_union == nullptr) {
			if (first_union != nullptr) free(first_union);
			fprintf(stderr, "any_is_excluded ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < first_union_length; i++) {
			first_union[i] = term->any.excluded_trees[i];
			first_union[i]->reference_count++;
		} for (unsigned int i = 0; i < second_union_length; i++) {
			second_union[i] = excluded_trees[i];
			second_union[i]->reference_count++;
		}
		if (!reduce_union<BuiltInPredicates>(first_union, first_union_length, second_union, second_union_length)) {
			fprintf(stderr, "any_is_excluded ERROR: `reduce_union` failed.\n");
			for (unsigned int i = 0; i < first_union_length; i++) {
				free(*first_union[i]); if (first_union[i]->reference_count == 0) free(first_union[i]);
			} for (unsigned int i = 0; i < second_union_length; i++) {
				free(*second_union[i]); if (second_union[i]->reference_count == 0) free(second_union[i]);
			}
			free(first_union); free(second_union);
			return false;
		}

		for (unsigned int i = 0; i < second_union_length; i++) {
			if (term == first->any.included && index_of(second_union[i], first->any.excluded_trees, first->any.excluded_tree_count) < first->any.excluded_tree_count)
				continue;
			if (term == second->any.included && index_of(second_union[i], second->any.excluded_trees, second->any.excluded_tree_count) < second->any.excluded_tree_count)
				continue;
			if (is_subset<BuiltInPredicates>(term, second_union[i])) {
				/* `term` is a subset of the union of `excluded_trees[i]` */
				for (unsigned int j = 0; j < first_union_length; j++) {
					free(*first_union[j]); if (first_union[j]->reference_count == 0) free(first_union[j]);
				} for (unsigned int j = 0; j < second_union_length; j++) {
					free(*second_union[j]); if (second_union[j]->reference_count == 0) free(second_union[j]);
				}
				free(first_union); free(second_union);
				return false;
			}
		}

		/* `term` is *not* a subset of the union of `excluded_trees[i]` */
		for (unsigned int i = 0; i < first_union_length; i++) {
			free(*first_union[i]);
			if (first_union[i]->reference_count == 0)
				free(first_union[i]);
		} for (unsigned int i = 0; i < second_union_length; i++) {
			free(*second_union[i]);
			if (second_union[i]->reference_count == 0)
				free(second_union[i]);
		}
		free(first_union); free(second_union);

	} else {
		if (term != nullptr
		&& (term != first->any.included || !same_trees_as_first)
		&& (term != second->any.included || !same_trees_as_second))
		{
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				if (term == first->any.included && index_of(excluded_trees[i], first->any.excluded_trees, first->any.excluded_tree_count) < first->any.excluded_tree_count)
					continue;
				if (term == second->any.included && index_of(excluded_trees[i], second->any.excluded_trees, second->any.excluded_tree_count) < second->any.excluded_tree_count)
					continue;
				if (is_subset<BuiltInPredicates>(term, excluded_trees[i])) {
					/* `term` is a subset of the union of `excluded_trees[i]` */
					return false;
				}
			}
		}
	}
	return true;
}

template<typename BuiltInPredicates, bool ComputeIntersection>
inline bool intersect_any_with_any(array<hol_term*>& dst,
		hol_term* first, hol_term* second, hol_term* (&included)[2],
		hol_term** additional_excluded, size_t additional_excluded_count)
{
	unsigned int excluded_tree_count = 0;
	bool same_trees_as_first = true;
	bool same_trees_as_second = true;
	hol_term** excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * (first->any.excluded_tree_count + second->any.excluded_tree_count + additional_excluded_count));
	for (unsigned int i = 0; i < first->any.excluded_tree_count; i++) {
		bool irreducible = true;
		for (unsigned int j = 0; j < second->any.excluded_tree_count; j++) {
			if (is_subset<BuiltInPredicates>(first->any.excluded_trees[i], second->any.excluded_trees[j])) {
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			excluded_trees[excluded_tree_count++] = first->any.excluded_trees[i];
			same_trees_as_second = false;
		} else {
			same_trees_as_first = false;
		}
	}
	unsigned int old_excluded_tree_count = excluded_tree_count;
	for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
		bool irreducible = true;
		for (unsigned int j = 0; j < old_excluded_tree_count; j++) {
			if (is_subset<BuiltInPredicates>(second->any.excluded_trees[i], excluded_trees[j])) {
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			excluded_trees[excluded_tree_count++] = second->any.excluded_trees[i];
			same_trees_as_first = false;
		} else {
			same_trees_as_second = false;
		}
	}
	old_excluded_tree_count = excluded_tree_count;
	for (unsigned int i = 0; i < additional_excluded_count; i++) {
		bool irreducible = true;
		for (unsigned int j = 0; j < old_excluded_tree_count; j++) {
			if (is_subset<BuiltInPredicates>(additional_excluded[i], excluded_trees[j])) {
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			excluded_trees[excluded_tree_count++] = additional_excluded[i];
			same_trees_as_first = false;
			same_trees_as_second = false;
		}
	}

	/* reduce the `excluded_trees` union */
	for (unsigned int i = 0; i < excluded_tree_count; i++)
		excluded_trees[i]->reference_count++;
	if (!same_trees_as_first && !same_trees_as_second) {
		reduce_union<BuiltInPredicates>(excluded_trees, excluded_tree_count);
		if (excluded_tree_count == 0) {
			free(excluded_trees);
			excluded_trees = nullptr;
		}
	}

	/* for each k, check whether `included[k]` is a subset of the union of `excluded_trees[i]` */
	for (unsigned int k = 0; k < array_length(included); k++) {
		if (!any_is_excluded<BuiltInPredicates>(included[k], first, second, excluded_trees, excluded_tree_count, same_trees_as_first, same_trees_as_second)) {
			for (unsigned int j = 0; j < excluded_tree_count; j++) {
				free(*excluded_trees[j]); if (excluded_trees[j]->reference_count == 0) free(excluded_trees[j]);
			}
			free(excluded_trees); return false;
		}
	}

	if (!ComputeIntersection) {
		/* we know now that the intersection is non-empty, so return */
		if (excluded_trees != nullptr) {
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]);
				if (excluded_trees[i]->reference_count == 0)
					free(excluded_trees[i]);
			}
			free(excluded_trees);
		}
		return true;
	} else if (included[1] != nullptr) {
		fprintf(stderr, "intersect_any_with_any ERROR: Unclosed intersection.\n");
		if (excluded_trees != nullptr) {
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]);
				if (excluded_trees[i]->reference_count == 0)
					free(excluded_trees[i]);
			}
			free(excluded_trees);
		}
		return false;
	}

	if (included[0] == first->any.included && same_trees_as_first) {
		if (excluded_trees != nullptr) {
			for (unsigned int j = 0; j < excluded_tree_count; j++) {
				free(*excluded_trees[j]);
				if (excluded_trees[j]->reference_count == 0)
					free(excluded_trees[j]);
			}
			free(excluded_trees);
		}
		if (!dst.ensure_capacity(dst.length + 1)) return false;
		dst[dst.length] = first;
		dst[dst.length++]->reference_count++;
		return true;
	} if (included[0] == second->any.included && same_trees_as_second) {
		if (excluded_trees != nullptr) {
			for (unsigned int j = 0; j < excluded_tree_count; j++) {
				free(*excluded_trees[j]);
				if (excluded_trees[j]->reference_count == 0)
					free(excluded_trees[j]);
			}
			free(excluded_trees);
		}
		if (!dst.ensure_capacity(dst.length + 1)) return false;
		dst[dst.length] = second;
		dst[dst.length++]->reference_count++;
		return true;
	}

	if (!dst.ensure_capacity(dst.length + 1)
	 || !new_hol_term(dst[dst.length]))
	{
		if (excluded_trees != nullptr) {
			for (unsigned int j = 0; j < excluded_tree_count; j++) {
				free(*excluded_trees[j]);
				if (excluded_trees[j]->reference_count == 0)
					free(excluded_trees[j]);
			}
			free(excluded_trees);
		}
		return false;
	}
	dst[dst.length]->type = (included[0] == first->any.included ? first->type : second->type);
	dst[dst.length]->reference_count = 1;
	dst[dst.length]->any.included = included[0];
	if (included[0] != nullptr) included[0]->reference_count++;
	dst[dst.length]->any.excluded_trees = excluded_trees;
	dst[dst.length]->any.excluded_tree_count = excluded_tree_count;
	dst.length++;
	return true;
}

template<typename BuiltInPredicates, bool ComputeIntersection>
inline bool intersect_any_with_any(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (first->type != hol_term_type::ANY && first->type != hol_term_type::ANY_RIGHT) {
		fprintf(stderr, "intersect_any_with_any ERROR: Expected `first` to be of type `ANY` or `ANY_RIGHT`.\n");
		return false;
	} if (second->type != hol_term_type::ANY && second->type != hol_term_type::ANY_RIGHT) {
		fprintf(stderr, "intersect_any_with_any ERROR: Expected `second` to be of type `ANY` or `ANY_RIGHT`.\n");
		return false;
	}
	if (first->any.included != nullptr && *first->any.included == HOL_ANY)
		fprintf(stderr, "intersect_with_any WARNING: `first->any.included` is the set of all logical forms.\n");
	if (second->any.included != nullptr && *second->any.included == HOL_ANY)
		fprintf(stderr, "intersect_with_any WARNING: `second->any.included` is the set of all logical forms.\n");
#endif

	/* first intersect the `included` subset */
	hol_term* included[2];
	array<hol_term*> intersection(2);
	size_t first_intersection_count;
	if (first->any.included == nullptr) {
		if (second->any.included == nullptr) {
			included[0] = nullptr;
			included[1] = nullptr;
		} else {
			included[0] = second->any.included;
			included[1] = nullptr;
		}
	} else if (second->any.included == nullptr) {
		included[0] = first->any.included;
		included[1] = nullptr;
	} else if (first->type == hol_term_type::ANY_RIGHT && second->type == hol_term_type::ANY_RIGHT) {
		if (second->any.excluded_tree_count == 0) {
			intersect<BuiltInPredicates, true>(intersection, first->any.included, second);
		} else {
			hol_term* term = hol_term::new_any_right(second->any.included);
			if (term == nullptr) return false;
			second->any.included->reference_count++;
			intersect<BuiltInPredicates, true>(intersection, first->any.included, term);
			free(*term); if (term->reference_count == 0) free(term);
		}
		first_intersection_count = intersection.length;

		if (first->any.excluded_tree_count == 0) {
			intersect<BuiltInPredicates, true>(intersection, second->any.included, first);
		} else {
			hol_term* term = hol_term::new_any_right(first->any.included);
			if (term == nullptr) return false;
			first->any.included->reference_count++;
			intersect<BuiltInPredicates, true>(intersection, second->any.included, term);
			free(*term); if (term->reference_count == 0) free(term);
		}
		if (intersection.length == 0)
			return false;
	} else if (!(first->type == hol_term_type::ANY && second->type == hol_term_type::ANY_RIGHT) && is_subset<BuiltInPredicates>(first->any.included, second)) {
		included[0] = first->any.included;
		included[1] = nullptr;
	} else if (!(second->type == hol_term_type::ANY && first->type == hol_term_type::ANY_RIGHT) && is_subset<BuiltInPredicates>(second->any.included, first)) {
		included[0] = second->any.included;
		included[1] = nullptr;
	} else {
		included[0] = first->any.included;
		included[1] = second->any.included;
	}

	if (intersection.length != 0) {
		array<hol_term*> additional_excluded(max((size_t) 1, first_intersection_count));
		for (size_t i = 0; i < first_intersection_count; i++) {
			if (intersection[i] == first->any.included) {
				additional_excluded[additional_excluded.length++] = first;
				first->reference_count++;
			} else {
				hol_term* term = hol_term::new_any_right(intersection[i]);
				if (term == nullptr) {
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : additional_excluded) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				intersection[i]->reference_count++;
				additional_excluded[additional_excluded.length++] = term;
			}
		}

		included[1] = nullptr;
		size_t old_dst_length = dst.length;
		for (size_t i = 0; i < intersection.length; i++) {
			included[0] = intersection[i];
			bool intersection_not_empty;
			if (i < first_intersection_count)
				intersection_not_empty = intersect_any_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second, included, nullptr, 0);
			else intersection_not_empty = intersect_any_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second, included, additional_excluded.data, additional_excluded.length);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : additional_excluded) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			}
		}
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : additional_excluded) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	} else {
		return intersect_any_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second, included, nullptr, 0);
	}
}

template<typename BuiltInPredicates, bool ComputeIntersection>
inline bool intersect_with_any(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (second->type != hol_term_type::ANY)
		fprintf(stderr, "intersect_with_any WARNING: Expected the type of `second` to be `ANY`.\n");
#endif

	if (first == second || (second->any.included == nullptr && second->any.excluded_tree_count == 0)) {
		if (ComputeIntersection) {
			if (!dst.add(first)) return false;
			first->reference_count++;
		}
		return true;
	} else if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) {
		return intersect_any_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	}

	/* first subtract `second->any.excluded[i]` from `first` */
	array<hol_term*> differences(8);
	differences[0] = first;
	differences[0]->reference_count++;
	differences.length++;
	for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
		array<hol_term*> new_differences(8);
		for (hol_term* prev : differences)
			subtract<BuiltInPredicates>(new_differences, prev, second->any.excluded_trees[i]);
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		swap(new_differences, differences);
		if (differences.length == 0) return false;
	}

	hol_term* second_any;
	if (second->any.included == nullptr) {
		if (!ComputeIntersection) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		} else if (!dst.append(differences.data, differences.length)) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		return true;
	} else if (second->any.excluded_tree_count == 0) {
		second_any = second;
		second->reference_count++;
	} else {
		second_any = hol_term::new_any(second->any.included);
		if (second_any == nullptr) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		if (second->any.included != nullptr)
			second->any.included->reference_count++;
	}

	size_t old_dst_length = dst.length;
	switch (first->type) {
	case hol_term_type::NOT:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->unary.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->unary.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->unary.operand = first_child;
				first_child->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> first_differences(8);
			subtract_any<BuiltInPredicates>(first_differences, root->unary.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->unary.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = first->type;
					term->reference_count = 1;
					term->unary.operand = first_difference;
					first_difference->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->binary.left, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->binary.left) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->binary.left = first_child;
				dst[dst.length]->binary.right = root->binary.right;
				dst[dst.length]->binary.left->reference_count++;
				dst[dst.length]->binary.right->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> second_intersection(8);
			array<hol_term*> first_differences(8);
			intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(second_intersection, root->binary.right, second_any);
			subtract_any<BuiltInPredicates>(first_differences, root->binary.left, second_any);
			if (!ComputeIntersection && intersection_not_empty && first_differences.length > 0) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			} if (ComputeIntersection && !dst.ensure_capacity(dst.length + (first_differences.length * second_intersection.length))) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				for (hol_term* second_child : second_intersection) {
					if (first_difference == root->binary.left && second_child == root->binary.right) {
						dst[dst.length] = root;
						root->reference_count++;
					} else {
						if (!new_hol_term(dst[dst.length])) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							return false;
						}
						dst[dst.length]->type = first->type;
						dst[dst.length]->reference_count = 1;
						dst[dst.length]->binary.left = first_difference;
						dst[dst.length]->binary.right = second_child;
						first_difference->reference_count++;
						second_child->reference_count++;
					}
					dst.length++;
				}
			}
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> second_differences(8);
			subtract_any<BuiltInPredicates>(second_differences, root->binary.right, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + (first_differences.length * second_differences.length))) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				for (hol_term* second_difference : second_differences) {
					hol_term* term;
					if (first_difference == root->binary.left && second_difference == root->binary.right) {
						term = root;
						root->reference_count++;
					} else {
						if (!new_hol_term(term)) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							return false;
						}
						term->type = first->type;
						term->reference_count = 1;
						term->binary.left = first_difference;
						term->binary.right = second_difference;
						first_difference->reference_count++;
						second_difference->reference_count++;
					}
					intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
					free(*term); if (term->reference_count == 0) free(term);
					if (!ComputeIntersection && intersection_not_empty) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return true;
					}
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::BINARY_APPLICATION:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->ternary.first, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->ternary.first) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->ternary.first = first_child;
				dst[dst.length]->ternary.second = root->ternary.second;
				dst[dst.length]->ternary.third = root->ternary.third;
				dst[dst.length]->ternary.first->reference_count++;
				dst[dst.length]->ternary.second->reference_count++;
				dst[dst.length]->ternary.third->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> second_intersection(8);
			array<hol_term*> first_differences(8);
			intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(second_intersection, root->ternary.second, second_any);
			subtract_any<BuiltInPredicates>(first_differences, root->ternary.first, second_any);
			if (!ComputeIntersection && intersection_not_empty && first_differences.length > 0) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			} if (ComputeIntersection && !dst.ensure_capacity(dst.length + (first_differences.length * second_intersection.length))) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				for (hol_term* second_child : second_intersection) {
					if (first_difference == root->ternary.first && second_child == root->ternary.second) {
						dst[dst.length] = root;
						root->reference_count++;
					} else {
						if (!new_hol_term(dst[dst.length])) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							return false;
						}
						dst[dst.length]->type = first->type;
						dst[dst.length]->reference_count = 1;
						dst[dst.length]->ternary.first = first_difference;
						dst[dst.length]->ternary.second = second_child;
						dst[dst.length]->ternary.third = root->ternary.third;
						dst[dst.length]->ternary.first->reference_count++;
						dst[dst.length]->ternary.second->reference_count++;
						dst[dst.length]->ternary.third->reference_count++;
					}
					dst.length++;
				}
			}
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> third_intersection(8);
			array<hol_term*> second_differences(8);
			intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(third_intersection, root->ternary.third, second_any);
			subtract_any<BuiltInPredicates>(second_differences, root->ternary.second, second_any);
			if (!ComputeIntersection && intersection_not_empty && second_differences.length > 0) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			} if (ComputeIntersection && !dst.ensure_capacity(dst.length + (first_differences.length * second_differences.length * third_intersection.length))) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				for (hol_term* second_difference : second_differences) {
					for (hol_term* third_child : third_intersection) {
						if (first_difference == root->ternary.first && second_difference == root->ternary.second && third_child == root->ternary.third) {
							dst[dst.length] = root;
							root->reference_count++;
						} else {
							if (!new_hol_term(dst[dst.length])) {
								for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
								free(*second_any); if (second_any->reference_count == 0) free(second_any);
								return false;
							}
							dst[dst.length]->type = first->type;
							dst[dst.length]->reference_count = 1;
							dst[dst.length]->ternary.first = first_difference;
							dst[dst.length]->ternary.second = second_difference;
							dst[dst.length]->ternary.third = third_child;
							dst[dst.length]->ternary.first->reference_count++;
							dst[dst.length]->ternary.second->reference_count++;
							dst[dst.length]->ternary.third->reference_count++;
						}
						dst.length++;
					}
				}
			}
			for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> third_differences(8);
			subtract_any<BuiltInPredicates>(third_differences, root->ternary.third, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + (first_differences.length * second_differences.length * third_differences.length))) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				for (hol_term* second_difference : third_differences) {
					for (hol_term* third_difference : third_differences) {
						hol_term* term;
						if (first_difference == root->ternary.first && second_difference == root->ternary.second && third_difference == root->ternary.third) {
							term = root;
							root->reference_count++;
						} else {
							if (!new_hol_term(term)) {
								for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
								free(*second_any); if (second_any->reference_count == 0) free(second_any);
								return false;
							}
							term->type = first->type;
							term->reference_count = 1;
							term->ternary.first = first_difference;
							term->ternary.second = second_difference;
							term->ternary.third = third_difference;
							first_difference->reference_count++;
							second_difference->reference_count++;
							third_difference->reference_count++;
						}
						intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
						free(*term); if (term->reference_count == 0) free(term);
						if (!ComputeIntersection && intersection_not_empty) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							return true;
						}
					}
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_differences) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : third_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		for (hol_term* root : differences) {
			array<hol_term*>* child_differences = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * root->array.length);
			if (child_differences == nullptr) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (unsigned int i = 0; i < root->array.length; i++) {
				if (!array_init(child_differences[i], 8)) {
					for (unsigned int j = 0; j < i; j++) free(child_differences[j]);
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					free(child_differences);
					return false;
				}
			}

			bool skip = false;
			unsigned int child_difference_count = 1;
			for (unsigned int i = 0; i < root->array.length && child_difference_count > 0; i++) {
				array<hol_term*> child_intersection(8);
				bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(child_intersection, root->array.operands[i], second_any);
				if (i > 0)
					subtract_any<BuiltInPredicates>(child_differences[i - 1], root->array.operands[i - 1], second_any);
				if (i > 0 && !ComputeIntersection && intersection_not_empty && child_differences[i - 1].length > 0) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
				if (i > 0) child_difference_count *= child_differences[i - 1].length;
				unsigned int intersection_count = child_intersection.length * child_difference_count;
				if (intersection_count == 0) {
					for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				} else if (ComputeIntersection && !dst.ensure_capacity(dst.length + intersection_count)) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				} else if (i == 0 && intersection_count == 1) {
					if (child_intersection[0] == root->array.operands[0]) {
						dst[dst.length++] = root;
						root->reference_count++;
						for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						continue;
					}
				}
				unsigned int* index_array = (unsigned int*) calloc(max(1u, i), sizeof(unsigned int));
				if (index_array == nullptr) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				while (true) {
					for (hol_term* child : child_intersection) {
						if (!new_hol_term(dst[dst.length])) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (unsigned int j = 0; j < root->array.length; j++) {
								for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
								free(child_differences[j]);
							}
							free(child_differences); free(index_array);
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							return false;
						}
						dst[dst.length]->type = first->type;
						dst[dst.length]->reference_count = 1;
						dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
						if (dst[dst.length]->array.operands == nullptr) {
							for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
							for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
							for (unsigned int j = 0; j < root->array.length; j++) {
								for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
								free(child_differences[j]);
							}
							free(child_differences); free(index_array);
							free(*second_any); if (second_any->reference_count == 0) free(second_any);
							free(dst[dst.length]); return false;
						}
						for (unsigned int j = 0; j < i; j++) {
							dst[dst.length]->array.operands[j] = child_differences[j][index_array[j]];
							dst[dst.length]->array.operands[j]->reference_count++;
						}
						dst[dst.length]->array.operands[i] = child;
						child->reference_count++;
						for (unsigned int j = i + 1; j < first->array.length; j++) {
							dst[dst.length]->array.operands[j] = root->array.operands[j];
							dst[dst.length]->array.operands[j]->reference_count++;
						}
						dst.length++;
					}

					/* increment `index_array` */
					bool remaining = false;
					for (unsigned int j = i; j > 0; j--) {
						index_array[j - 1]++;
						if (index_array[j - 1] == child_differences[j - 1].length) {
							index_array[j - 1] = 0;
						} else {
							remaining = true;
							break;
						}
					}
					if (!remaining) break;
				}
				for (hol_term* term : child_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(index_array);
			}
			if (skip || child_difference_count == 0) {
				for (unsigned int j = 0; j < root->array.length; j++) {
					for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(child_differences[j]);
				}
				free(child_differences);
				continue;
			}

			subtract_any<BuiltInPredicates>(child_differences[first->array.length - 1], root->array.operands[first->array.length - 1], second_any);
			unsigned int intersection_count = child_difference_count * child_differences[first->array.length - 1].length;
			if (intersection_count == 0) {
				for (unsigned int j = 0; j < root->array.length; j++) {
					for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(child_differences[j]);
				}
				free(child_differences);
				continue;
			} else if (ComputeIntersection && !dst.ensure_capacity(dst.length + intersection_count)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (unsigned int j = 0; j < root->array.length; j++) {
					for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(child_differences[j]);
				}
				free(child_differences);
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (intersection_count == 1) {
				bool same = false;
				for (unsigned int i = 0; i < first->array.length; i++)
					if (child_differences[i][0] != root->array.operands[i]) same = false;
				if (same) {
					bool intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, root);
					if (!ComputeIntersection && intersection_not_empty) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (unsigned int j = 0; j < root->array.length; j++) {
							for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
							free(child_differences[j]);
						}
						free(child_differences);
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return true;
					}
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences);
					continue;
				}
			}
			unsigned int* index_array = (unsigned int*) calloc(first->array.length, sizeof(unsigned int));
			if (index_array == nullptr) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (unsigned int j = 0; j < root->array.length; j++) {
					for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(child_differences[j]);
				}
				free(child_differences);
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			while (true) {
				hol_term* term;
				if (!new_hol_term(term)) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences); free(index_array);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				term->type = first->type;
				term->reference_count = 1;
				term->array.length = first->array.length;
				term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
				if (term->array.operands == nullptr) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences); free(index_array); free(term);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				for (unsigned int j = 0; j < first->array.length; j++) {
					term->array.operands[j] = child_differences[j][index_array[j]];
					term->array.operands[j]->reference_count++;
				}
				bool intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < root->array.length; j++) {
						for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(child_differences[j]);
					}
					free(child_differences); free(index_array);
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}

				/* increment `index_array` */
				bool remaining = false;
				for (unsigned int j = first->array.length; j > 0; j--) {
					index_array[j - 1]++;
					if (index_array[j - 1] == child_differences[j - 1].length) {
						index_array[j - 1] = 0;
					} else {
						remaining = true;
						break;
					}
				}
				if (!remaining) break;
			}
			for (unsigned int j = 0; j < root->array.length; j++) {
				for (hol_term* term : child_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(child_differences[j]);
			}
			free(child_differences); free(index_array);
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->quantifier.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->quantifier.variable = root->quantifier.variable;
				dst[dst.length]->quantifier.operand = first_child;
				first_child->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> first_differences(8);
			subtract_any<BuiltInPredicates>(first_differences, root->quantifier.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->quantifier.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = first->type;
					term->reference_count = 1;
					term->quantifier.variable = root->quantifier.variable;
					term->quantifier.operand = first_difference;
					first_difference->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		for (hol_term* root : differences) {
			bool intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, root, second->any.included);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY_ARRAY:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->any_array.all, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (ComputeIntersection && root->any_array.any.length != 0 && first_intersection.length != 0) {
				fprintf(stderr, "intersect_with_any ERROR: Unclosed intersection.\n");
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->any_array.all) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				dst[dst.length] = hol_term::new_any_array(root->any_array.oper, root->any_array.all, make_array_view(&first_child, 1),
						make_array_view(root->any_array.left.operands, root->any_array.left.length),
						make_array_view(root->any_array.right.operands, root->any_array.right.length));
				if (dst[dst.length] == 0) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->any_array.all->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.any.length; i++)
					dst[dst.length]->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.left.length; i++)
					dst[dst.length]->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.right.length; i++)
					dst[dst.length]->any_array.right.operands[i]->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> first_differences(8);
			subtract_any<BuiltInPredicates>(first_differences, root->any_array.all, second_any);
			for (hol_term* all : first_differences) {
				bool intersection_not_empty = false;
				bool result = subtract_any_with_array<BuiltInPredicates>(root->any_array.left, all, second_any, [&intersection_not_empty,root,first,all,second_any,&dst](array<hol_term*>* left_array, const unsigned int* left_index_array) {
					return subtract_any_with_array<BuiltInPredicates>(root->any_array.any, all, second_any, [&intersection_not_empty,root,first,all,second_any,left_array,left_index_array,&dst](array<hol_term*>* any_array, const unsigned int* any_index_array) {
						return subtract_any_with_array<BuiltInPredicates>(root->any_array.right, all, second_any, [&intersection_not_empty,root,first,all,second_any,left_array,left_index_array,any_array,any_index_array,&dst](array<hol_term*>* right_array, const unsigned int* right_index_array) {
							hol_term* term;
							bool same_as_root = (all == root->any_array.all);
							for (unsigned int i = 0; same_as_root && i < root->any_array.left.length; i++)
								if (root->any_array.left.operands[i] != left_array[i][left_index_array[i]]) same_as_root = false;
							for (unsigned int i = 0; same_as_root && i < root->any_array.right.length; i++)
								if (root->any_array.right.operands[i] != right_array[i][right_index_array[i]]) same_as_root = false;
							for (unsigned int i = 0; same_as_root && i < root->any_array.any.length; i++)
								if (root->any_array.any.operands[i] != any_array[i][any_index_array[i]]) same_as_root = false;
							if (same_as_root) {
								term = root;
								root->reference_count++;
							} else {
								term = hol_term::new_any_array(root->any_array.oper, all,
										lookup_table_array_view<array<hol_term*>, hol_term*>(any_array, any_index_array, root->any_array.any.length),
										lookup_table_array_view<array<hol_term*>, hol_term*>(left_array, left_index_array, root->any_array.left.length),
										lookup_table_array_view<array<hol_term*>, hol_term*>(right_array, right_index_array, root->any_array.right.length));
								if (term == nullptr) return false;
								term->any_array.all->reference_count++;
								for (unsigned int i = 0; i < term->any_array.any.length; i++)
									term->any_array.any.operands[i]->reference_count++;
								for (unsigned int i = 0; i < term->any_array.left.length; i++)
									term->any_array.left.operands[i]->reference_count++;
								for (unsigned int i = 0; i < term->any_array.right.length; i++)
									term->any_array.right.operands[i]->reference_count++;
							}

							intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second_any->any.included, term);
							free(*term); if (term->reference_count == 0) free(term);
							if (!ComputeIntersection && intersection_not_empty)
								return false;
							return true;
						});
					});
				});
				if (!result) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					if (!ComputeIntersection && intersection_not_empty)
						return true;
					else return false;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY_QUANTIFIER:
		for (hol_term* root : differences) {
			array<hol_term*> first_intersection(8);
			bool intersection_not_empty = intersect_with_any<BuiltInPredicates, ComputeIntersection>(first_intersection, root->any_quantifier.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->any_quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					continue;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->any_quantifier.quantifier = root->any_quantifier.quantifier;
				dst[dst.length]->any_quantifier.operand = first_child;
				first_child->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> first_differences(8);
			subtract_any<BuiltInPredicates>(first_differences, root->any_quantifier.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->any_quantifier.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = first->type;
					term->reference_count = 1;
					term->any_quantifier.quantifier = root->any_quantifier.quantifier;
					term->any_quantifier.operand = first_difference;
					first_difference->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		/* we assume `first` is not of type `ANY` */
		break;
	}
	fprintf(stderr, "intersect_with_any ERROR: Unrecognized hol_term_type.\n");
	for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
	free(*second_any); if (second_any->reference_count == 0) free(second_any);
	return false;
}

template<typename BuiltInPredicates, bool ComputeIntersection>
inline bool intersect_with_any_right(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (second->type != hol_term_type::ANY_RIGHT)
		fprintf(stderr, "intersect_with_any_right WARNING: Expected the type of `second` to be `ANY_RIGHT`.\n");
#endif

	if (first == second) {
		if (ComputeIntersection) {
			if (!dst.add(first)) return false;
			first->reference_count++;
		}
		return true;
	} else if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) {
		return intersect_any_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	}

	/* first subtract `second->any.excluded[i]` from `first` */
	array<hol_term*> differences(8);
	differences[0] = first;
	differences[0]->reference_count++;
	differences.length++;
	for (unsigned int i = 0; i < second->any.excluded_tree_count; i++) {
		array<hol_term*> new_differences(8);
		for (hol_term* prev : differences)
			subtract<BuiltInPredicates>(new_differences, prev, second->any.excluded_trees[i]);
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		swap(new_differences, differences);
		if (differences.length == 0) return false;
	}

	hol_term* second_any;
	if (second->any.included == nullptr) {
		if (!ComputeIntersection) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		} else if (!dst.append(differences.data, differences.length)) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		return true;
	} else if (second->any.excluded_tree_count == 0) {
		second_any = second;
		second->reference_count++;
	} else {
		second_any = hol_term::new_any_right(second->any.included);
		if (second_any == nullptr) {
			for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		if (second->any.included != nullptr)
			second->any.included->reference_count++;
	}

	bool intersection_not_empty;
	array<hol_term*> first_intersection(8);
	array<hol_term*> first_differences(8);
	size_t old_dst_length = dst.length;
	switch (first->type) {
	case hol_term_type::NOT:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->unary.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->unary.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->unary.operand = first_child;
				first_child->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			subtract_any_right<BuiltInPredicates>(first_differences, root->unary.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->unary.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = root->type;
					term->reference_count = 1;
					term->unary.operand = first_difference;
					first_difference->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->binary.right, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->binary.right) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->binary.left = root->binary.left;
				dst[dst.length]->binary.right = first_child;
				dst[dst.length]->binary.left->reference_count++;
				dst[dst.length]->binary.right->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			array<hol_term*> intersections(4);
			intersect<BuiltInPredicates, true>(intersections, second->any.included, root);
			for (hol_term* intersection : intersections) {
				subtract_any_right<BuiltInPredicates>(first_differences, intersection->binary.right, second_any);
				if (!ComputeIntersection && first_differences.length > 0) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				} else if (!dst.ensure_capacity(dst.length + first_differences.length)) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				} else if (first_differences.length == 1 && first_differences[0] == intersection->binary.right) {
					dst[dst.length] = intersection;
					dst[dst.length++]->reference_count++;
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					first_differences.clear();
					continue;
				}
				for (hol_term* first_difference : first_differences) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					dst[dst.length]->type = first->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->binary.left = intersection->binary.left;
					dst[dst.length]->binary.right = first_difference;
					dst[dst.length]->binary.left->reference_count++;
					dst[dst.length]->binary.right->reference_count++;
					dst.length++;
				}
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				first_differences.clear();
			}
			for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::BINARY_APPLICATION:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->ternary.third, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->ternary.third) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->ternary.first = root->ternary.first;
				dst[dst.length]->ternary.second = root->ternary.second;
				dst[dst.length]->ternary.third = first_child;
				dst[dst.length]->ternary.first->reference_count++;
				dst[dst.length]->ternary.second->reference_count++;
				dst[dst.length]->ternary.third->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			subtract_any_right<BuiltInPredicates>(first_differences, root->ternary.third, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->ternary.third) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = root->type;
					term->reference_count = 1;
					term->ternary.first = root->ternary.first;
					term->ternary.second = root->ternary.second;
					term->ternary.third = first_difference;
					term->ternary.first->reference_count++;
					term->ternary.second->reference_count++;
					term->ternary.third->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->array.operands[root->array.length - 1], second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->array.operands[root->array.length - 1]) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->array.length = root->array.length;
				dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * root->array.length);
				if (dst[dst.length]->array.operands == nullptr) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					free(dst[dst.length]); return false;
				}
				for (unsigned int i = 0; i + 1 < root->array.length; i++) {
					dst[dst.length]->array.operands[i] = root->array.operands[i];
					dst[dst.length]->array.operands[i]->reference_count++;
				}
				dst[dst.length]->array.operands[root->array.length - 1] = first_child;
				dst[dst.length]->array.operands[root->array.length - 1]->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			subtract_any_right<BuiltInPredicates>(first_differences, root->array.operands[root->array.length - 1], second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->array.operands[root->array.length - 1]) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = first->type;
					term->reference_count = 1;
					term->array.length = root->array.length;
					term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * root->array.length);
					if (term->array.operands == nullptr) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						free(term); return false;
					}
					for (unsigned int i = 0; i + 1 < root->array.length; i++) {
						term->array.operands[i] = root->array.operands[i];
						term->array.operands[i]->reference_count++;
					}
					term->array.operands[root->array.length - 1] = first_difference;
					term->array.operands[root->array.length - 1]->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->quantifier.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->quantifier.variable = root->quantifier.variable;
				dst[dst.length]->quantifier.operand = first_child;
				dst[dst.length]->quantifier.operand->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			subtract_any_right<BuiltInPredicates>(first_differences, root->quantifier.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->quantifier.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = first->type;
					term->reference_count = 1;
					term->quantifier.variable = root->quantifier.variable;
					term->quantifier.operand = first_difference;
					term->quantifier.operand->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, root, second->any.included);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY_ARRAY:
		for (hol_term* root : differences) {
			if (root->any_array.right.length != 0) {
				intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->any_array.right.operands[root->any_array.right.length - 1], second_any);
			} else {
				intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->any_array.all, second_any);
			}
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if ((root->any_array.right.length != 0 && first_intersection[0] == root->any_array.right.operands[root->any_array.right.length - 1])
				 || (root->any_array.right.length == 0 && first_intersection[0] == root->any_array.all))
				{
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				dst[dst.length] = hol_term::new_any_array(root->any_array.oper, root->any_array.all,
						make_array_view(root->any_array.any.operands, root->any_array.any.length),
						make_array_view(root->any_array.left.operands, root->any_array.left.length),
						make_appended_array_view(make_array_view(root->any_array.right.operands, max(1u, root->any_array.right.length) - 1), first_child));
				if (dst[dst.length] == nullptr) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->any_array.all->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.any.length; i++)
					dst[dst.length]->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.left.length; i++)
					dst[dst.length]->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst[dst.length]->any_array.right.length; i++)
					dst[dst.length]->any_array.right.operands[i]->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			if (root->any_array.right.length != 0) {
				subtract_any_right<BuiltInPredicates>(first_differences, root->any_array.right.operands[root->any_array.right.length - 1], second_any);
			} else {
				subtract_any_right<BuiltInPredicates>(first_differences, root->any_array.all, second_any);
			}
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if ((root->any_array.right.length != 0 && first_difference == root->any_array.right.operands[root->any_array.right.length - 1])
				 || (root->any_array.right.length == 0 && first_difference == root->any_array.all))
				{
					term = root;
					root->reference_count++;
				} else {
					term = hol_term::new_any_array(root->any_array.oper, root->any_array.all,
							make_array_view(root->any_array.any.operands, root->any_array.any.length),
							make_array_view(root->any_array.left.operands, root->any_array.left.length),
							make_appended_array_view(make_array_view(root->any_array.right.operands, max(1u, root->any_array.right.length) - 1), first_difference));
					if (term == nullptr) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->any_array.all->reference_count++;
					for (unsigned int i = 0; i < term->any_array.any.length; i++)
						term->any_array.any.operands[i]->reference_count++;
					for (unsigned int i = 0; i < term->any_array.left.length; i++)
						term->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < term->any_array.right.length; i++)
						term->any_array.right.operands[i]->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY_QUANTIFIER:
		for (hol_term* root : differences) {
			intersection_not_empty = intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(first_intersection, root->any_quantifier.operand, second_any);
			if (!ComputeIntersection && intersection_not_empty) {
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return true;
			}

			if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			} else if (first_intersection.length == 1) {
				if (first_intersection[0] == root->any_quantifier.operand) {
					dst[dst.length++] = root;
					root->reference_count++;
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* first_child : first_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->any_quantifier.quantifier = root->any_quantifier.quantifier;
				dst[dst.length]->any_quantifier.operand = first_child;
				dst[dst.length]->any_quantifier.operand->reference_count++;
				dst.length++;
			}
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }

			subtract_any_right<BuiltInPredicates>(first_differences, root->any_quantifier.operand, second_any);
			if (ComputeIntersection && !dst.ensure_capacity(dst.length + first_differences.length)) {
				for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*second_any); if (second_any->reference_count == 0) free(second_any);
				return false;
			}
			for (hol_term* first_difference : first_differences) {
				hol_term* term;
				if (first_difference == root->any_quantifier.operand) {
					term = root;
					root->reference_count++;
				} else {
					if (!new_hol_term(term)) {
						for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						free(*second_any); if (second_any->reference_count == 0) free(second_any);
						return false;
					}
					term->type = root->type;
					term->reference_count = 1;
					term->any_quantifier.quantifier = root->any_quantifier.quantifier;
					term->any_quantifier.operand = first_difference;
					term->any_quantifier.operand->reference_count++;
				}
				intersection_not_empty = intersect<BuiltInPredicates, ComputeIntersection>(dst, second->any.included, term);
				free(*term); if (term->reference_count == 0) free(term);
				if (!ComputeIntersection && intersection_not_empty) {
					for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
					free(*second_any); if (second_any->reference_count == 0) free(second_any);
					return true;
				}
			}
			for (hol_term* term : first_differences) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
		free(*second_any); if (second_any->reference_count == 0) free(second_any);
		return (dst.length > old_dst_length);
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		/* we assume `first` is not of type `ANY` or `ANY_ARRAY` */
		break;
	}
	fprintf(stderr, "intersect_with_any_right ERROR: Unrecognized hol_term_type.\n");
	for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
	free(*second_any); if (second_any->reference_count == 0) free(second_any);
	return false;
}

template<typename BuiltInPredicates, bool ComputeIntersection>
bool intersect_with_any_array(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
#if !defined(NDEBUG)
	if (second->type != hol_term_type::ANY_ARRAY)
		fprintf(stderr, "intersect_with_any_array WARNING: Expected the type of `second` to be `ANY_ARRAY`.\n");
	if (first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT)
		fprintf(stderr, "intersect_with_any_array WARNING: Expected the type of `first` to not be `ANY` or `ANY_RIGHT`.\n");
#endif

	if (first->type == hol_term_type::ANY_ARRAY) {
		bool same_as_first = true;
		bool same_as_second = true;
		hol_term_type oper;
		if (first->any_array.oper == hol_term_type::ANY_ARRAY) {
			if (second->any_array.oper == hol_term_type::ANY_ARRAY) {
				oper = hol_term_type::ANY_ARRAY;
			} else {
				oper = second->any_array.oper;
				same_as_first = false;
			}
		} else {
			if (second->any_array.oper == hol_term_type::ANY_ARRAY) {
				oper = first->any_array.oper;
				same_as_second = false;
			} else if (first->any_array.oper == second->any_array.oper) {
				oper = first->any_array.oper;
			} else {
				return false;
			}
		}

		hol_term** any;
		unsigned int any_length;
		if (first->any_array.any.length == 0) {
			if (second->any_array.any.length == 0) {
				any = nullptr;
				any_length = 0;
			} else {
				any = (hol_term**) malloc(sizeof(hol_term*) * second->any_array.any.length);
				if (any == nullptr) {
					fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
					return false;
				}
				any_length = second->any_array.any.length;
				for (unsigned int i = 0; i < any_length; i++) {
					any[i] = second->any_array.any.operands[i];
					any[i]->reference_count++;
				}
				same_as_first = false;
			}
		} else {
			if (second->any_array.any.length == 0) {
				any = (hol_term**) malloc(sizeof(hol_term*) * first->any_array.any.length);
				if (any == nullptr) {
					fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
					return false;
				}
				any_length = first->any_array.any.length;
				for (unsigned int i = 0; i < any_length; i++) {
					any[i] = first->any_array.any.operands[i];
					any[i]->reference_count++;
				}
				same_as_second = false;
			} else {
				fprintf(stderr, "intersect_with_any_array ERROR: Unclosed intersection.\n");
				return false;
			}
		}

		unsigned int new_left_length = max(first->any_array.left.length, second->any_array.left.length);
		unsigned int new_right_length = max(first->any_array.right.length, second->any_array.right.length);
		if (new_left_length != first->any_array.left.length || new_right_length != first->any_array.right.length)
			same_as_first = false;
		if (new_left_length != second->any_array.left.length || new_right_length != second->any_array.right.length)
			same_as_second = false;

		array<hol_term*> all_intersections(8);
		array<hol_term*>* left_intersections = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * new_left_length));
		array<hol_term*>* right_intersections = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * new_right_length));
		if (left_intersections == nullptr || right_intersections == nullptr) {
			fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
			if (left_intersections != nullptr) free(left_intersections);
			if (any != nullptr) {
				for (unsigned int j = 0; j < any_length; j++) {
					free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
				}
				free(any);
			}
			return false;
		} for (unsigned int i = 0; i < new_left_length; i++) {
			if (!array_init(left_intersections[i], 4)) {
				for (unsigned int j = 0; j < i; j++) free(left_intersections[j]);
				free(left_intersections); free(right_intersections);
				if (any != nullptr) {
					for (unsigned int j = 0; j < any_length; j++) {
						free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
					}
					free(any);
				}
				return false;
			}
		} for (unsigned int i = 0; i < new_right_length; i++) {
			if (!array_init(right_intersections[i], 4)) {
				for (unsigned int j = 0; j < new_left_length; j++) free(left_intersections[j]);
				for (unsigned int j = 0; j < i; j++) free(right_intersections[j]);
				free(left_intersections); free(right_intersections);
				if (any != nullptr) {
					for (unsigned int j = 0; j < any_length; j++) {
						free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
					}
					free(any);
				}
				return false;
			}
		}

		if (!intersect<BuiltInPredicates, ComputeIntersection>(all_intersections, first->any_array.all, second->any_array.all))
			return false;
		for (unsigned int i = 0; i < min(first->any_array.left.length, second->any_array.left.length); i++) {
			if (!intersect<BuiltInPredicates, ComputeIntersection>(left_intersections[i], first->any_array.left.operands[i], second->any_array.left.operands[i])) {
				for (hol_term* term : all_intersections) { free(*term); if (term->reference_count == 0) free(term); }
				for (unsigned int j = 0; j < new_left_length; j++) {
					for (hol_term* term : left_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(left_intersections[j]);
				} for (unsigned int j = 0; j < new_right_length; j++) {
					for (hol_term* term : right_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(right_intersections[j]);
				}
				free(left_intersections); free(right_intersections);
				if (any != nullptr) {
					for (unsigned int j = 0; j < any_length; j++) {
						free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
					}
					free(any);
				}
				return false;
			}
		} for (unsigned int i = first->any_array.left.length; i < new_left_length; i++) {
			left_intersections[i].add(second->any_array.right.operands[i]);
		} for (unsigned int i = second->any_array.left.length; i < new_left_length; i++) {
			left_intersections[i].add(first->any_array.right.operands[i]);
		}
		for (unsigned int i = 0; i < min(first->any_array.right.length, second->any_array.right.length); i++) {
			if (!intersect<BuiltInPredicates, ComputeIntersection>(right_intersections[i], first->any_array.right.operands[i], second->any_array.right.operands[i])) {
				for (hol_term* term : all_intersections) { free(*term); if (term->reference_count == 0) free(term); }
				for (unsigned int j = 0; j < new_left_length; j++) {
					for (hol_term* term : left_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(left_intersections[j]);
				} for (unsigned int j = 0; j < new_right_length; j++) {
					for (hol_term* term : right_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(right_intersections[j]);
				}
				free(left_intersections); free(right_intersections);
				if (any != nullptr) {
					for (unsigned int j = 0; j < any_length; j++) {
						free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
					}
					free(any);
				}
				return false;
			}
		} for (unsigned int i = first->any_array.right.length; i < new_right_length; i++) {
			right_intersections[i].add(second->any_array.right.operands[i]);
		} for (unsigned int i = second->any_array.right.length; i < new_right_length; i++) {
			right_intersections[i].add(first->any_array.right.operands[i]);
		}
		if (!ComputeIntersection) {
			for (unsigned int j = 0; j < new_left_length; j++) {
				free(left_intersections[j]);
			} for (unsigned int j = 0; j < new_right_length; j++) {
				free(right_intersections[j]);
			}
			free(left_intersections); free(right_intersections);
			if (any != nullptr) {
				for (unsigned int j = 0; j < any_length; j++) {
					free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
				}
				free(any);
			}
			return true;
		}

		bool success = true;
		for (hol_term* all : all_intersections) {
			array<hol_term*>* any_intersections = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * any_length));
			if (any_intersections == nullptr) {
				success = false;
				break;
			}
			for (unsigned int i = 0; i < any_length; i++) {
				if (!array_init(any_intersections[i], 4)) {
					success = false;
					for (unsigned int j = 0; j < i; j++) free(any_intersections[j]);
					free(any_intersections);
					break;
				}
			}
			if (!success) break;
			bool skip = false;
			for (unsigned int i = 0; i < any_length; i++) {
				intersect<BuiltInPredicates, ComputeIntersection>(any_intersections[i], any[i], all);
				if (any_intersections[i].length == 0) {
					for (unsigned int j = 0; j < any_length; j++) {
						for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(any_intersections[j]);
					}
					free(any_intersections);
					skip = true;
					break;
				}
			}
			if (skip) continue;

			array<hol_term*>* new_left_intersections = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * new_left_length);
			if (new_left_intersections == nullptr) {
				success = false;
				for (unsigned int j = 0; j < any_length; j++) {
					for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(any_intersections[j]);
				}
				free(any_intersections);
				break;
			}
			for (unsigned int i = 0; i < new_left_length; i++) {
				if (!array_init(new_left_intersections[i], left_intersections[i].length)) {
					success = false;
					for (unsigned int j = 0; j < any_length; j++) {
						for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(any_intersections[j]);
					}
					free(any_intersections);
					for (unsigned int j = 0; j < i; j++) free(new_left_intersections[j]);
					free(new_left_intersections);
					break;
				}
				for (hol_term* term : left_intersections[i])
					if (is_subset<BuiltInPredicates>(term, all))
						new_left_intersections[i].add(term);
				if (new_left_intersections[i].length == 0) {
					for (unsigned int j = 0; j < any_length; j++) {
						for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(any_intersections[j]);
					}
					free(any_intersections);
					for (unsigned int j = 0; j < i; j++) free(new_left_intersections[j]);
					free(new_left_intersections);
					skip = true;
					break;
				}
			}
			if (!success) break;
			if (skip) continue;

			array<hol_term*>* new_right_intersections = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * new_right_length);
			if (new_right_intersections == nullptr) {
				success = false;
				for (unsigned int j = 0; j < any_length; j++) {
					for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(any_intersections[j]);
				}
				free(any_intersections);
				for (unsigned int j = 0; j < new_left_length; j++) free(new_left_intersections[j]);
				free(new_left_intersections);
				break;
			}
			for (unsigned int i = 0; i < new_right_length; i++) {
				if (!array_init(new_right_intersections[i], right_intersections[i].length)) {
					success = false;
					for (unsigned int j = 0; j < any_length; j++) {
						for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(any_intersections[j]);
					}
					free(any_intersections);
					for (unsigned int j = 0; j < new_left_length; j++) free(new_left_intersections[j]);
					for (unsigned int j = 0; j < i; j++) free(new_right_intersections[j]);
					free(new_left_intersections); free(new_right_intersections);
					break;
				}
				for (hol_term* term : right_intersections[i])
					if (is_subset<BuiltInPredicates>(term, all))
						new_right_intersections[i].add(term);
				if (new_right_intersections[i].length == 0) {
					for (unsigned int j = 0; j < any_length; j++) {
						for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(any_intersections[j]);
					}
					free(any_intersections);
					for (unsigned int j = 0; j < new_left_length; j++) free(new_left_intersections[j]);
					for (unsigned int j = 0; j < i; j++) free(new_right_intersections[j]);
					free(new_left_intersections); free(new_right_intersections);
					skip = true;
					break;
				}
			}
			if (!success) break;
			if (skip) continue;

			bool result = apply_to_cartesian_product(new_left_intersections, new_left_length, [first,second,same_as_first,same_as_second,all,oper,&dst,&any_intersections,&new_left_intersections,&new_right_intersections,any_length,new_left_length,new_right_length](const unsigned int* left_indices) {
				return apply_to_cartesian_product(new_right_intersections, new_right_length, [first,second,same_as_first,same_as_second,all,oper,&dst,&any_intersections,&new_left_intersections,&new_right_intersections,any_length,new_left_length,new_right_length,left_indices](const unsigned int* right_indices) {
					return apply_to_cartesian_product(any_intersections, any_length, [first,second,same_as_first,same_as_second,all,oper,&dst,&any_intersections,&new_left_intersections,&new_right_intersections,any_length,new_left_length,new_right_length,left_indices,right_indices](const unsigned int* any_indices) {
						/* check if the new logical form is the same as `first` or `second` */
						if (same_as_first) {
							bool same = (all == first->any_array.all);
							for (unsigned int i = 0; same && i < new_left_length; i++)
								if (first->any_array.left.operands[i] != new_left_intersections[i][left_indices[i]]) same = false;
							for (unsigned int i = 0; same && i < new_right_length; i++)
								if (first->any_array.right.operands[i] != new_right_intersections[i][right_indices[i]]) same = false;
							for (unsigned int i = 0; same && i < any_length; i++)
								if (first->any_array.any.operands[i] != any_intersections[i][any_indices[i]]) same = false;
							if (same) {
								dst[dst.length] = first;
								first->reference_count++;
								return true;
							}
						} if (same_as_second) {
							bool same = (all == second->any_array.all);
							for (unsigned int i = 0; same && i < new_left_length; i++)
								if (second->any_array.left.operands[i] != new_left_intersections[i][left_indices[i]]) same = false;
							for (unsigned int i = 0; same && i < new_right_length; i++)
								if (second->any_array.right.operands[i] != new_right_intersections[i][right_indices[i]]) same = false;
							for (unsigned int i = 0; same && i < any_length; i++)
								if (second->any_array.any.operands[i] != any_intersections[i][any_indices[i]]) same = false;
							if (same) {
								dst[dst.length] = second;
								second->reference_count++;
								return true;
							}
						}

						if (!new_hol_term(dst[dst.length]))
							return false;
						dst[dst.length]->type = hol_term_type::ANY_ARRAY;
						dst[dst.length]->reference_count = 1;
						dst[dst.length]->any_array.oper = oper;
						dst[dst.length]->any_array.all = all;
						all->reference_count++;
						dst[dst.length]->any_array.any.length = 0;
						dst[dst.length]->any_array.any.operands = nullptr;
						dst[dst.length]->any_array.left.length = 0;
						dst[dst.length]->any_array.left.operands = nullptr;
						dst[dst.length]->any_array.right.length = 0;
						dst[dst.length]->any_array.right.operands = nullptr;

						dst[dst.length]->any_array.any.length = any_length;
						if (any_length == 0) {
							dst[dst.length]->any_array.any.operands = nullptr;
						} else {
							dst[dst.length]->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * any_length);
							if (dst[dst.length]->any_array.any.operands == nullptr) {
								fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
								free(dst[dst.length]); return false;
							}
							for (unsigned int i = 0; i < any_length; i++) {
								dst[dst.length]->any_array.any.operands[i] = any_intersections[i][any_indices[i]];
								dst[dst.length]->any_array.any.operands[i]->reference_count++;
							}
						}
						dst[dst.length]->any_array.left.length = new_left_length;
						if (new_left_length == 0) {
							dst[dst.length]->any_array.left.operands = nullptr;
						} else {
							dst[dst.length]->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * new_left_length);
							if (dst[dst.length]->any_array.left.operands == nullptr) {
								fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
								free(*dst[dst.length]); free(dst[dst.length]);
								return false;
							}
							for (unsigned int i = 0; i < new_left_length; i++) {
								dst[dst.length]->any_array.left.operands[i] = new_left_intersections[i][left_indices[i]];
								dst[dst.length]->any_array.left.operands[i]->reference_count++;
							}
						}
						dst[dst.length]->any_array.right.length = new_right_length;
						if (new_right_length == 0) {
							dst[dst.length]->any_array.right.operands = nullptr;
						} else {
							dst[dst.length]->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * new_right_length);
							if (dst[dst.length]->any_array.right.operands == nullptr) {
								fprintf(stderr, "intersect_with_any_array ERROR: Out of memory.\n");
								free(*dst[dst.length]); free(dst[dst.length]);
								return false;
							}
							for (unsigned int i = 0; i < new_right_length; i++) {
								dst[dst.length]->any_array.right.operands[i] = new_right_intersections[i][right_indices[i]];
								dst[dst.length]->any_array.right.operands[i]->reference_count++;
							}
						}
						dst.length++;
						return true;
					});
				});
			});

			for (unsigned int j = 0; j < any_length; j++) {
				for (hol_term* term : any_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(any_intersections[j]);
			}
			for (unsigned int j = 0; j < new_left_length; j++) free(new_left_intersections[j]);
			for (unsigned int j = 0; j < new_right_length; j++) free(new_right_intersections[j]);
			free(any_intersections); free(new_left_intersections); free(new_right_intersections);
			if (!result) {
				success = false;
				break;
			}
		}

		for (hol_term* term : all_intersections) { free(*term); if (term->reference_count == 0) free(term); }
		for (unsigned int j = 0; j < new_left_length; j++) {
			for (hol_term* term : left_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
			free(left_intersections[j]);
		} for (unsigned int j = 0; j < new_right_length; j++) {
			for (hol_term* term : right_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
			free(right_intersections[j]);
		}
		free(left_intersections); free(right_intersections);
		if (any != nullptr) {
			for (unsigned int j = 0; j < any_length; j++) {
				free(*any[j]); if (any[j]->reference_count == 0) free(any[j]);
			}
			free(any);
		}
		return success && (dst.length > 0);

	} else if (first->type == hol_term_type::AND || first->type == hol_term_type::OR || first->type == hol_term_type::IFF) {
		if (second->any_array.oper != hol_term_type::ANY_ARRAY && second->any_array.oper != first->type)
			return false;
		if (first->array.length < second->any_array.left.length || first->array.length < second->any_array.right.length || first->array.length < second->any_array.any.length)
			return false;
		array<hol_term*>* terms = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * first->array.length);
		if (terms == nullptr) {
			fprintf(stderr, "intersect ERROR: Insufficient memory for `terms`.\n");
			if (terms != nullptr) free(terms);
			return false;
		}
		for (unsigned int i = 0; i < first->array.length; i++) {
			if (!array_init(terms[i], 8)) {
				for (unsigned int j = 0; j < i; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms);
				return false;
			}

			array<hol_term*> second_terms(2);
			if (i < second->any_array.left.length) {
				if (first->array.length - i - 1 < second->any_array.right.length) {
					/* both `second->any_array.left` and `second->any_array.right` apply to this term */
					if (!intersect<BuiltInPredicates, true>(second_terms, second->any_array.left.operands[i], second->any_array.right.operands[second->any_array.right.length - first->array.length + i])) {
						for (unsigned int j = 0; j <= i; j++) {
							for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
							free(terms[j]);
						}
						free(terms);
						return false;
					}
				} else {
					/* only the term from `second->any_array.left` applies to this term */
					second_terms[second_terms.length] = second->any_array.left.operands[i];
					second_terms[second_terms.length++]->reference_count++;
				}
			} else if (first->array.length - i - 1 < second->any_array.right.length) {
				/* only the term from `second->any_array.right` applies to this term */
				second_terms[second_terms.length] = second->any_array.right.operands[second->any_array.right.length - first->array.length + i];
				second_terms[second_terms.length++]->reference_count++;
			} else {
				/* no term from `second->any_array.left` or `second->any_array.right` applies to this term */
				second_terms[second_terms.length] = second->any_array.all;
				second_terms[second_terms.length++]->reference_count++;
			}

			for (hol_term* second_term : second_terms) {
				if (!intersect<BuiltInPredicates, true>(terms[i], first->array.operands[i], second_term)) {
					for (unsigned int j = 0; j <= i; j++) {
						for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(terms[j]);
					}
					free(terms);
					for (hol_term* term : second_terms) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
			}
			for (hol_term* term : second_terms) { free(*term); if (term->reference_count == 0) free(term); }
		}

		bool non_empty_intersection = false;
		bool result = apply_to_cartesian_product(terms, first->array.length, [first,second,&terms,&dst,&non_empty_intersection](const unsigned int* index_array) {
			/* handle the case where `second->any_array.any.length` is 0 */
			if (second->any_array.any.length == 0) {
				if (!ComputeIntersection) {
					non_empty_intersection = true;
					return true;
				} else if (!dst.ensure_capacity(dst.length + 1)) {
					return false;
				}

				/* check if the new term is the same as `first` */
				bool same_as_first = true;
				for (unsigned int j = 0; j < first->array.length; j++)
					if (terms[j][index_array[j]] != first->array.operands[j]) same_as_first = false;
				if (same_as_first) {
					dst[dst.length++] = first;
					first->reference_count++;
					return true;
				}

				if (!new_hol_term(dst[dst.length])) return false;
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->array.length = first->array.length;
				dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
				if (dst[dst.length]->array.operands == nullptr) {
					free(dst[dst.length]);
					return false;
				}
				for (unsigned int j = 0; j < first->array.length; j++) {
					dst[dst.length]->array.operands[j] = terms[j][index_array[j]];
					dst[dst.length]->array.operands[j]->reference_count++;
				}
				dst.length++;
				return true;
			}

			/* so now we can handle the case where `second->any_array.any.length` > 0 */
			for (unsigned int i = 0; i < first->array.length - second->any_array.any.length + 1; i++) {
				array<hol_term*>* child_intersections = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * second->any_array.any.length);
				if (child_intersections == nullptr) return false;

				for (unsigned int k = 0; k < second->any_array.any.length; k++) {
					if (!array_init(child_intersections[k], 4)) {
						for (unsigned int j = 0; j < k; j++) free(child_intersections[j]);
						free(child_intersections);
						return false;
					}
				}

				non_empty_intersection = true;
				for (unsigned int k = 0; k < second->any_array.any.length; k++)
					non_empty_intersection &= intersect<BuiltInPredicates, ComputeIntersection>(child_intersections[k], terms[i + k][index_array[i + k]], second->any_array.any.operands[k]);

				if (!non_empty_intersection) {
					for (unsigned int j = 0; j < second->any_array.any.length; j++) free(child_intersections[j]);
					free(child_intersections);
					continue;
				} else if (non_empty_intersection && !ComputeIntersection) {
					for (unsigned int j = 0; j < second->any_array.any.length; j++) free(child_intersections[j]);
					free(child_intersections);
					return true;
				}

				bool result = apply_to_cartesian_product(child_intersections, second->any_array.any.length, [first,second,i,index_array,&terms,&child_intersections,&dst](const unsigned int* inner_index_array) {
					hol_term* new_term;
					if (!new_hol_term(new_term)) return false;
					new_term->type = first->type;
					new_term->reference_count = 1;
					new_term->array.length = first->array.length;
					new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
					if (new_term->array.operands == nullptr) {
						free(new_term);
						return false;
					}
					for (unsigned int j = 0; j < i; j++)
						new_term->array.operands[j] = terms[j][index_array[j]];
					for (unsigned int k = 0; k < second->any_array.any.length; k++)
						new_term->array.operands[i + k] = child_intersections[k][inner_index_array[k]];
					for (unsigned int j = i + second->any_array.any.length; j < first->array.length; j++)
						new_term->array.operands[j] = terms[j][index_array[j]];
					for (unsigned int j = 0; j < first->array.length; j++)
						new_term->array.operands[j]->reference_count++;

					/* check if `new_term` is the same as `first` */
					bool same_as_first = true;
					for (unsigned int j = 0; j < first->array.length; j++)
						if (new_term->array.operands[j] != first->array.operands[j]) same_as_first = false;
					if (same_as_first) {
						free(*new_term); free(new_term);
						new_term = first;
						first->reference_count++;
					}

					/* to ensure that `new_term` is disjoint with all other sets in `dst`, we just subtract them */
					array<hol_term*> differences(8);
					differences[differences.length++] = new_term;
					unsigned int old_length = dst.length;
					for (unsigned int j = 0; j < old_length && differences.length > 0; j++) {
						array<hol_term*> new_differences(8);
						for (hol_term* difference : differences) {
							subtract<BuiltInPredicates>(new_differences, difference, dst[j]);
							free(*difference); if (difference->reference_count == 0) free(difference);
						}
						swap(differences, new_differences);
					}

					if (!dst.append(differences.data, differences.length)) {
						for (hol_term* term : differences) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					return true;
				});

				for (unsigned int j = 0; j < second->any_array.any.length; j++) {
					for (hol_term* term : child_intersections[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(child_intersections[j]);
				}
				free(child_intersections);
				if (!result) return false;
			}
			return true;
		});
		for (unsigned int j = 0; j < first->array.length; j++) {
			for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
			free(terms[j]);
		}
		free(terms);
		if (!ComputeIntersection && non_empty_intersection)
			return true;
		return result && (dst.length > 0);

	} else {

		if (second->any_array.left.length > 1 || second->any_array.right.length > 1 || second->any_array.any.length > 1)
			return false;

		if (second->any_array.left.length == 0 && second->any_array.right.length == 0 && second->any_array.any.length == 0)
			return intersect<BuiltInPredicates, ComputeIntersection>(dst, first, second->any_array.all);

		/* `second` could be an array of length 1 */
		array<hol_term*> left_intersection(8);
		if (second->any_array.left.length == 0) {
			left_intersection[left_intersection.length++] = first;
			first->reference_count++;
		} else if (!intersect<BuiltInPredicates, true>(left_intersection, first, second->any_array.left.operands[0])) {
			return false;
		}
		size_t old_dst_length = dst.length;
		for (hol_term* left : left_intersection) {
			array<hol_term*> right_intersection(8);
			if (second->any_array.right.length == 0) {
				right_intersection[right_intersection.length++] = left;
				left->reference_count++;
			} else if (!intersect<BuiltInPredicates, true>(right_intersection, left, second->any_array.right.operands[0])) {
				continue;
			}
			for (hol_term* right : right_intersection) {
				if (second->any_array.any.length == 0) {
					if (!ComputeIntersection) {
						for (hol_term* term : right_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : left_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return true;
					} else if (!dst.ensure_capacity(dst.length + 1)) {
						for (hol_term* term : right_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : left_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length++] = right;
					right->reference_count++;
				} else if (!intersect<BuiltInPredicates, ComputeIntersection>(dst, right, second->any_array.any.operands[0])) {
					continue;
				}
			}
			for (hol_term* term : right_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		}
		for (hol_term* term : left_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return (dst.length > old_dst_length);
	}
}

template<typename BuiltInPredicates, bool ComputeIntersection>
inline bool intersect_any_quantifier(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
	hol_quantifier_type second_quantifier;
	hol_term* second_operand;
	unsigned int second_variable;
	if (second->type == hol_term_type::ANY_QUANTIFIER) {
		second_quantifier = second->any_quantifier.quantifier;
		second_operand = second->any_quantifier.operand;
		second_variable = ANY_VARIABLE;
	} else if (second->type == hol_term_type::FOR_ALL || second->type == hol_term_type::EXISTS || second->type == hol_term_type::LAMBDA) {
		second_quantifier = (hol_quantifier_type) second->type;
		second_operand = second->quantifier.operand;
		second_variable = second->quantifier.variable;
	} else {
		return false;
	}

	hol_quantifier_type quantifier;
	if (!intersect(quantifier, first->any_quantifier.quantifier, second_quantifier))
		return false;
	array<hol_term*> intersection(8);
	if (!intersect<BuiltInPredicates, ComputeIntersection>(intersection, first->any_quantifier.operand, second_operand))
		return false;
	if (!ComputeIntersection) {
		return true;
	} else if (!dst.ensure_capacity(dst.length + intersection.length)) {
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	} else if (intersection.length == 1) {
		if (first->any_quantifier.quantifier == quantifier && second_variable == ANY_VARIABLE && first->any_quantifier.operand == intersection[0]) {
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		} if (second_quantifier == quantifier && second_operand == intersection[0]) {
			dst[dst.length] = second;
			dst[dst.length++]->reference_count++;
			return true;
		}
	}

	for (hol_term* term : intersection) {
		if (!new_hol_term(dst[dst.length])) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		if (second_variable == ANY_VARIABLE) {
			dst[dst.length]->type = hol_term_type::ANY_QUANTIFIER;
			dst[dst.length]->any_quantifier.quantifier = quantifier;
			dst[dst.length]->any_quantifier.operand = term;
		} else {
			dst[dst.length]->type = (hol_term_type) quantifier;
			dst[dst.length]->quantifier.variable = second_variable;
			dst[dst.length]->quantifier.operand = term;
		}
		term->reference_count++;
		dst[dst.length]->reference_count = 1;
		dst.length++;
	}
	return true;
}

template<typename BuiltInPredicates, bool ComputeIntersection>
bool intersect(array<hol_term*>& dst, hol_term* first, hol_term* second)
{
	if (first == second) {
		if (ComputeIntersection) {
			if (!dst.add(first)) return false;
			first->reference_count++;
		}
		return true;
	} else if (first->type == hol_term_type::ANY) {
		return intersect_with_any<BuiltInPredicates, ComputeIntersection>(dst, second, first);
	} else if (second->type == hol_term_type::ANY) {
		return intersect_with_any<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	} else if (first->type == hol_term_type::ANY_RIGHT) {
		return intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(dst, second, first);
	} else if (second->type == hol_term_type::ANY_RIGHT) {
		return intersect_with_any_right<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	} else if (first->type == hol_term_type::ANY_ARRAY) {
		return intersect_with_any_array<BuiltInPredicates, ComputeIntersection>(dst, second, first);
	} else if (second->type == hol_term_type::ANY_ARRAY) {
		return intersect_with_any_array<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	} else if (first->type == hol_term_type::ANY_CONSTANT) {
		if (second->type == hol_term_type::ANY_CONSTANT || second->type == hol_term_type::CONSTANT) {
			unsigned int* second_constants; unsigned int second_constant_count;
			if (second->type == hol_term_type::CONSTANT) {
				second_constants = &second->constant;
				second_constant_count = 1;
			} else {
				second_constants = second->any_constant.constants;
				second_constant_count = second->any_constant.length;
			}

			if (!ComputeIntersection)
				return has_intersection(first->any_constant.constants, first->any_constant.length, second_constants, second_constant_count);
			if (!dst.ensure_capacity(dst.length + 1))
				return false;

			unsigned int* intersection = (unsigned int*) malloc(sizeof(unsigned int) * min(first->any_constant.length, second_constant_count));
			if (intersection == nullptr) {
				fprintf(stderr, "intersect ERROR: Insufficient memory for `intersection`.\n");
				return false;
			}
			unsigned int intersection_size = 0;
			set_intersect(intersection, intersection_size, first->any_constant.constants, first->any_constant.length, second_constants, second_constant_count);
			if (intersection_size == first->any_constant.length) {
				free(intersection);
				dst[dst.length] = first;
				dst[dst.length++]->reference_count++;
				return true;
			} else if (intersection_size == second_constant_count) {
				free(intersection);
				dst[dst.length] = second;
				dst[dst.length++]->reference_count++;
				return true;
			} else if (intersection_size == 1) {
				dst[dst.length] = hol_term::new_constant(intersection[0]);
				free(intersection);
				if (dst[dst.length] == nullptr)
					return false;
				dst.length++;
				return true;
			} else if (!new_hol_term(dst[dst.length])) {
				free(intersection);
				return false;
			}
			dst[dst.length]->type = hol_term_type::ANY_CONSTANT;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->any_constant.constants = intersection;
			dst[dst.length]->any_constant.length = intersection_size;
			dst.length++;
			return true;
		} else {
			return false;
		}
	} else if (second->type == hol_term_type::ANY_CONSTANT) {
		if (first->type == hol_term_type::CONSTANT && index_of(first->constant, second->any_constant.constants, second->any_constant.length) < second->any_constant.length) {
			if (!ComputeIntersection) return true;
			if (!dst.ensure_capacity(dst.length + 1)) return false;
			dst[dst.length] = first;
			dst[dst.length++]->reference_count++;
			return true;
		} else {
			return false;
		}
	} else if (first->type == hol_term_type::ANY_QUANTIFIER) {
		return intersect_any_quantifier<BuiltInPredicates, ComputeIntersection>(dst, first, second);
	} else if (second->type == hol_term_type::ANY_QUANTIFIER) {
		return intersect_any_quantifier<BuiltInPredicates, ComputeIntersection>(dst, second, first);
	}

	array<hol_term*> first_intersection(8);
	array<hol_term*> second_intersection(8);
	array<hol_term*> third_intersection(8);
	array<hol_term*>* terms;
	unsigned int intersection_count;
	unsigned int* index_array;
	switch (first->type) {
	case hol_term_type::NOT:
		if (second->type != hol_term_type::NOT
		 || !intersect<BuiltInPredicates, ComputeIntersection>(first_intersection, first->unary.operand, second->unary.operand))
			return false;
		if (!ComputeIntersection) {
			return true;
		} else if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_intersection.length == 1) {
			if (first_intersection[0] == first->unary.operand) {
				dst[dst.length++] = first;
				first->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (first_intersection[0] == second->unary.operand) {
				dst[dst.length++] = second;
				second->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			}
		}
		for (hol_term* first_child : first_intersection) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = hol_term_type::NOT;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->unary.operand = first_child;
			first_child->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return true;

	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
		if (second->type != first->type
		 || !intersect<BuiltInPredicates, ComputeIntersection>(first_intersection, first->binary.left, second->binary.left))
		{
			return false;
		} else if (!intersect<BuiltInPredicates, ComputeIntersection>(second_intersection, first->binary.right, second->binary.right)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}

		if (!ComputeIntersection)
			return true;
		intersection_count = first_intersection.length * second_intersection.length;
		if (!dst.ensure_capacity(dst.length + intersection_count)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (intersection_count == 1) {
			if (first_intersection[0] == first->binary.left && second_intersection[0] == first->binary.right) {
				dst[dst.length++] = first;
				first->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (first_intersection[0] == second->binary.left && second_intersection[0] == second->binary.right) {
				dst[dst.length++] = second;
				second->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			}
		}
		for (hol_term* first_child : first_intersection) {
			for (hol_term* second_child : second_intersection) {
				if (!new_hol_term(dst[dst.length])) {
					for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				dst[dst.length]->type = first->type;
				dst[dst.length]->reference_count = 1;
				dst[dst.length]->binary.left = first_child;
				dst[dst.length]->binary.right = second_child;
				first_child->reference_count++;
				second_child->reference_count++;
				dst.length++;
			}
		}
		for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return true;

	case hol_term_type::BINARY_APPLICATION:
		if (second->type != first->type
		 || !intersect<BuiltInPredicates, ComputeIntersection>(first_intersection, first->ternary.first, second->ternary.first))
		{
			return false;
		} else if (!intersect<BuiltInPredicates, ComputeIntersection>(second_intersection, first->ternary.second, second->ternary.second)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (!intersect<BuiltInPredicates, ComputeIntersection>(third_intersection, first->ternary.third, second->ternary.third)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}

		if (!ComputeIntersection)
			return true;
		intersection_count = first_intersection.length * second_intersection.length * third_intersection.length;
		if (!dst.ensure_capacity(dst.length + intersection_count)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (intersection_count == 1) {
			if (first_intersection[0] == first->ternary.first && second_intersection[0] == first->ternary.second && third_intersection[0] == first->ternary.third) {
				dst[dst.length++] = first;
				first->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (first_intersection[0] == second->ternary.first && second_intersection[0] == second->ternary.second && third_intersection[0] == second->ternary.third) {
				dst[dst.length++] = second;
				second->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			}
		}
		for (hol_term* first_child : first_intersection) {
			for (hol_term* second_child : second_intersection) {
				for (hol_term* third_child : third_intersection) {
					if (!new_hol_term(dst[dst.length])) {
						for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					dst[dst.length]->type = first->type;
					dst[dst.length]->reference_count = 1;
					dst[dst.length]->ternary.first = first_child;
					dst[dst.length]->ternary.second = second_child;
					dst[dst.length]->ternary.third = third_child;
					first_child->reference_count++;
					second_child->reference_count++;
					third_child->reference_count++;
					dst.length++;
				}
			}
		}
		for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : third_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return true;

	case hol_term_type::IFF:
	case hol_term_type::AND:
	case hol_term_type::OR:
		if (second->type != first->type
		 || second->array.length != first->array.length)
			return false;
		if (ComputeIntersection) {
			terms = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * first->array.length);
			if (terms == nullptr) {
				fprintf(stderr, "intersect ERROR: Insufficient memory for `terms`.\n");
				return false;
			}
		}
		for (unsigned int i = 0; i < first->array.length; i++) {
			if (ComputeIntersection && !array_init(terms[i], 8)) {
				for (unsigned int j = 0; j < i; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms);
				return false;
			} else if (!intersect<BuiltInPredicates, ComputeIntersection>(ComputeIntersection ? terms[i] : first_intersection, first->array.operands[i], second->array.operands[i])) {
				if (ComputeIntersection) {
					for (unsigned int j = 0; j <= i; j++) {
						for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
						free(terms[j]);
					}
					free(terms);
				}
				return false;
			}
		}

		if (!ComputeIntersection)
			return true;
		intersection_count = terms[0].length;
		for (unsigned int i = 1; i < first->array.length; i++)
			intersection_count *= terms[i].length;
		if (intersection_count == 0 || !dst.ensure_capacity(dst.length + intersection_count)) {
			for (unsigned int j = 0; j < first->array.length; j++) {
				for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(terms[j]);
			}
			free(terms);
			return false;
		} else if (intersection_count == 1) {
			bool same_as_first = true;
			bool same_as_second = true;
			for (unsigned int i = 0; i < first->array.length; i++) {
				if (terms[i][0] != first->array.operands[i]) same_as_first = false;
				if (terms[i][0] != second->array.operands[i]) same_as_second = false;
			}
			if (same_as_first) {
				dst[dst.length++] = first;
				first->reference_count++;
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms);
				return true;
			} if (same_as_second) {
				dst[dst.length++] = second;
				second->reference_count++;
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms);
				return true;
			}
		}
		index_array = (unsigned int*) calloc(first->array.length, sizeof(unsigned int));
		if (index_array == nullptr) {
			fprintf(stderr, "intersect ERROR: Insufficient memory for `index_array`.\n");
			for (unsigned int j = 0; j < first->array.length; j++) {
				for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(terms[j]);
			}
			free(terms);
			return false;
		}
		while (true) {
			if (!new_hol_term(dst[dst.length])) {
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms); free(index_array);
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->array.length = first->array.length;
			dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first->array.length);
			if (dst[dst.length]->array.operands == nullptr) {
				for (unsigned int j = 0; j < first->array.length; j++) {
					for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
					free(terms[j]);
				}
				free(terms); free(index_array); free(dst[dst.length]);
				return false;
			}
			for (unsigned int i = 0; i < first->array.length; i++) {
				dst[dst.length]->array.operands[i] = terms[i][index_array[i]];
				dst[dst.length]->array.operands[i]->reference_count++;
			}
			dst.length++;

			/* increment `index_array` */
			bool remaining = false;
			for (unsigned int j = first->array.length; j > 0; j--) {
				index_array[j - 1]++;
				if (index_array[j - 1] == terms[j - 1].length) {
					index_array[j - 1] = 0;
				} else {
					remaining = true;
					break;
				}
			}
			if (!remaining) break;
		}
		for (unsigned int j = 0; j < first->array.length; j++) {
			for (hol_term* term : terms[j]) { free(*term); if (term->reference_count == 0) free(term); }
			free(terms[j]);
		}
		free(terms); free(index_array);
		return true;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (second->type != first->type || first->quantifier.variable != second->quantifier.variable
		 || !intersect<BuiltInPredicates, ComputeIntersection>(first_intersection, first->quantifier.operand, second->quantifier.operand))
			return false;
		if (!ComputeIntersection) {
			return true;
		} else if (!dst.ensure_capacity(dst.length + first_intersection.length)) {
			for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (first_intersection.length == 1) {
			if (first_intersection[0] == first->quantifier.operand) {
				dst[dst.length++] = first;
				first->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (first_intersection[0] == second->quantifier.operand) {
				dst[dst.length++] = second;
				second->reference_count++;
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			}
		}
		for (hol_term* first_child : first_intersection) {
			if (!new_hol_term(dst[dst.length])) {
				for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			dst[dst.length]->type = first->type;
			dst[dst.length]->reference_count = 1;
			dst[dst.length]->quantifier.variable = first->quantifier.variable;
			dst[dst.length]->quantifier.operand = first_child;
			first_child->reference_count++;
			dst.length++;
		}
		for (hol_term* term : first_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return true;

	case hol_term_type::INTEGER:
		if (second->type != hol_term_type::INTEGER
		 || second->integer != first->integer)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::STRING:
		if (second->type != hol_term_type::STRING
		 || second->str != first->str)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::UINT_LIST:
		if (second->type != hol_term_type::UINT_LIST
		 || second->uint_list != first->uint_list)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::CONSTANT:
		if (second->type != hol_term_type::CONSTANT
		 || second->constant != first->constant)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::VARIABLE:
		if (second->type != hol_term_type::VARIABLE
		 || second->variable != first->variable)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::PARAMETER:
		if (second->type != hol_term_type::PARAMETER
		 || second->parameter != first->parameter)
			return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		if (second->type != first->type) return false;
		if (!ComputeIntersection) return true;
		if (!dst.add(first)) return false;
		first->reference_count++;
		return true;
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_QUANTIFIER:
		/* we already handle this before the switch statement */
		break;
	}
	fprintf(stderr, "intersect ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<typename BuiltInPredicates>
inline bool has_intersection(hol_term* first, hol_term* second) {
	array<hol_term*> dummy(1);
#if !defined(NDEBUG)
	bool result = intersect<BuiltInPredicates, false>(dummy, first, second);
	if (dummy.length != 0)
		fprintf(stderr, "has_intersection WARNING: `intersect` returned non-empty array.\n");
	return result;
#else
	return intersect<BuiltInPredicates, false>(dummy, first, second);
#endif
}


/**
 * Code for tokenizing/lexing higher-order logic formulas in TPTP-like format.
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
	ARROW,
	IF_THEN,
	FOR_ALL,
	EXISTS,
	LAMBDA,
	EQUALS,

	IDENTIFIER,
	STRING,
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
	case tptp_token_type::ARROW:
		return print("->", stream);
	case tptp_token_type::IF_THEN:
		return print("=>", stream);
	case tptp_token_type::EQUALS:
		return print('=', stream);
	case tptp_token_type::FOR_ALL:
		return print('!', stream);
	case tptp_token_type::EXISTS:
		return print('?', stream);
	case tptp_token_type::LAMBDA:
		return print('^', stream);
	case tptp_token_type::SEMICOLON:
		return print(';', stream);
	case tptp_token_type::IDENTIFIER:
		return print("IDENTIFIER", stream);
	case tptp_token_type::STRING:
		return print("STRING", stream);
	}
	fprintf(stderr, "print ERROR: Unknown tptp_token_type.\n");
	return false;
}

enum class tptp_lexer_state {
	DEFAULT,
	IDENTIFIER,
	QUOTE
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
	case '^':
		return emit_token(tokens, start, start + 1, tptp_token_type::LAMBDA);
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
	 || next == '~' || next == '!' || next == '?' || next == '^'
	 || next == ';')
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
	} else if (next == '-') {
		fpos_t pos;
		fgetpos(input, &pos);
		next = fgetwc(input);
		if (next != '>') {
			read_error("Expected '>' after '-'", current);
			return false;
		} else {
			if (!emit_token(tokens, current, current + 2, tptp_token_type::ARROW)) return false;
		}
		current.column++;
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
			 || next == '^' || next == ';' || next == '-')
			{
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER)
				 || !tptp_lex_symbol(tokens, input, next, current))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear(); shift = {0};
			} else if (next == '"') {
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER))
					return false;
				state = tptp_lexer_state::QUOTE;
				tokens.clear(); shift = {0};
				start = current;
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

		case tptp_lexer_state::QUOTE:
			if (next == '\\') {
				/* this is an escape character */
				next = fgetwc(input);
				current.column++;
				new_line = (next == '\n');
				if (!append_to_token(token, next, shift)) return false;
			} else if (next == '"') {
				if (!emit_token(tokens, token, start, current, tptp_token_type::STRING))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear(); shift = {0};
			} else {
				if (!append_to_token(token, next, shift)) return false;
			}
			break;

		case tptp_lexer_state::DEFAULT:
			if (next == ',' || next == ':' || next == '(' || next == ')'
			 || next == '[' || next == ']' || next == '&' || next == '|'
			 || next == '~' || next == '!' || next == '?' || next == '='
			 || next == '^' || next == ';' || next == '-')
			{
				if (!tptp_lex_symbol(tokens, input, next, current))
					return false;
			} else if (next == '"') {
				state = tptp_lexer_state::QUOTE;
				tokens.clear(); shift = {0};
				start = current;
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
 * Recursive-descent parser for higher-order logic formulas in TPTP-like format.
 */

bool tptp_interpret_unary_term(
	const array<tptp_token>&,
	unsigned int&, hol_term&,
	hash_map<string, unsigned int>&,
	array_map<string, unsigned int>&);
bool tptp_interpret(
	const array<tptp_token>&,
	unsigned int&, hol_term&,
	hash_map<string, unsigned int>&,
	array_map<string, unsigned int>&);
bool tptp_interpret_type(
	const array<tptp_token>&,
	unsigned int&, hol_type&,
	hash_map<string, unsigned int>&,
	array_map<string, unsigned int>&);

bool tptp_interpret_argument_list(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	array<hol_term*>& terms)
{
	while (true) {
		if (!terms.ensure_capacity(terms.length + 1))
			return false;

		hol_term*& next_term = terms[terms.length];
		if (!new_hol_term(next_term)) return false;

		if (!tptp_interpret(tokens, index, *next_term, names, variables)) {
			free(next_term); return false;
		}
		terms.length++;

		if (index >= tokens.length || tokens[index].type != tptp_token_type::COMMA)
			return true;
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
		/*if (names.table.contains(tokens[index].text)) {
			fprintf(stderr, "WARNING at %d:%d: Variable '", tokens[index].start.line, tokens[index].start.column);
			print(tokens[index].text, stderr); print("' shadows previously declared identifier.\n", stderr);
		}*/
		if (variables.contains(tokens[index].text)) {
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

template<hol_term_type QuantifierType>
bool tptp_interpret_quantifier(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_term& term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	unsigned int old_variable_count = variables.size;
	if (!tptp_interpret_variable_list(tokens, index, names, variables)
	 || !expect_token(tokens, index, tptp_token_type::COLON, "colon for quantified term"))
		return false;
	index++;

	hol_term* operand;
	if (!new_hol_term(operand)) return false;
	operand->reference_count = 1;
	if (!tptp_interpret_unary_term(tokens, index, *operand, names, variables)) {
		free(operand); return false;
	}

	hol_term* inner = operand;
	for (unsigned int i = variables.size - 1; i > old_variable_count; i--) {
		hol_term* quantified;
		if (!new_hol_term(quantified)) {
			free(*inner); free(inner);
			return false;
		}
		quantified->quantifier.variable = variables.values[i];
		quantified->quantifier.operand = inner;
		quantified->type = QuantifierType;
		quantified->reference_count = 1;
		inner = quantified;
	}

	term.quantifier.variable = variables.values[old_variable_count];
	term.quantifier.operand = inner;
	term.type = QuantifierType;
	term.reference_count = 1;

	for (unsigned int i = old_variable_count; i < variables.size; i++)
		free(variables.keys[i]);
	variables.size = old_variable_count;
	return true;
}

bool tptp_interpret_unary_term(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_term& term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	if (index >= tokens.length) {
		fprintf(stderr, "ERROR: Unexpected end of input.\n");
		return false;

	} else if (tokens[index].type == tptp_token_type::NOT) {
		/* this is a negation of the form ~U */
		index++;
		hol_term* operand;
		if (!new_hol_term(operand)) return false;
		operand->reference_count = 1;
		if (!tptp_interpret_unary_term(tokens, index, *operand, names, variables)) {
			free(operand); return false;
		}

		term.unary.operand = operand;
		term.type = hol_term_type::NOT;
		term.reference_count = 1;

	} else if (tokens[index].type == tptp_token_type::LPAREN) {
		/* these are just grouping parenthesis of the form (F) */
		index++;
		if (!tptp_interpret(tokens, index, term, names, variables)) {
			return false;
		} if (!expect_token(tokens, index, tptp_token_type::RPAREN, "closing parenthesis")) {
			free(term); return false;
		}
		index++;

	} else if (tokens[index].type == tptp_token_type::FOR_ALL) {
		/* this is a universal quantifier of the form ![v_1,...,v_n]:U */
		index++;
		if (!tptp_interpret_quantifier<hol_term_type::FOR_ALL>(tokens, index, term, names, variables))
			return false;

	} else if (tokens[index].type == tptp_token_type::EXISTS) {
		/* this is an existential quantifier of the form ?[v_1,...,v_n]:U */
		index++;
		if (!tptp_interpret_quantifier<hol_term_type::EXISTS>(tokens, index, term, names, variables))
			return false;

	} else if (tokens[index].type == tptp_token_type::LAMBDA) {
		/* this is a lambda expression of the form ^[v_1,...,v_n]:U */
		index++;
		if (!tptp_interpret_quantifier<hol_term_type::LAMBDA>(tokens, index, term, names, variables))
			return false;

	} else if (tokens[index].type == tptp_token_type::STRING) {
		/* this is a string */
		term.type = hol_term_type::STRING;
		term.reference_count = 1;
		term.str = tokens[index].text;
		index++;

	} else if (tokens[index].type == tptp_token_type::IDENTIFIER) {
		int integer;
		if (tokens[index].text == "T") {
			/* this the constant true */
			term.type = hol_term_type::TRUE;
			term.reference_count = 1;
			index++;
		} else if (tokens[index].text == "F") {
			/* this the constant false */
			term.type = hol_term_type::FALSE;
			term.reference_count = 1;
			index++;
		} else if (parse_int(tokens[index].text, integer)) {
			/* this is an integer */
			term.integer = integer;
			term.type = hol_term_type::INTEGER;
			term.reference_count = 1;
			index++;
		} else {
			/* this is a constant or variable */
			bool contains;
			unsigned int variable = variables.get(tokens[index].text, contains);
			if (contains) {
				/* this argument is a variable */
				term.variable = variable;
				term.type = hol_term_type::VARIABLE;
				term.reference_count = 1;
			} else {
				if (!get_token(tokens[index].text, term.constant, names))
					return false;
				term.type = hol_term_type::CONSTANT;
				term.reference_count = 1;
			}
			index++;
		}

	} else {
		read_error("Unexpected symbol. Expected a unary term", tokens[index].start);
		return false;
	}

	while (index < tokens.length && tokens[index].type == tptp_token_type::LPAREN) {
		hol_term* left;
		if (!new_hol_term(left)) {
			free(term); return false;
		}
		move(term, *left);

		index++;
		array<hol_term*> terms = array<hol_term*>(4);
		if (!tptp_interpret_argument_list(tokens, index, names, variables, terms)) {
			free(*left); free(left);
			return false;
		} else if (!expect_token(tokens, index, tptp_token_type::RPAREN, "closing parenthesis for application")) {
			free(*left); free(left);
			for (hol_term* term : terms) {
				free(*term); if (term->reference_count == 0) free(term);
			}
			return false;
		}
		core::position lparen_position = tokens[index].start;
		index++;

		if (terms.length == 1) {
			term.type = hol_term_type::UNARY_APPLICATION;
			term.reference_count = 1;
			term.binary.left = left;
			term.binary.right = terms[0];
		} else if (terms.length == 2) {
			term.type = hol_term_type::BINARY_APPLICATION;
			term.reference_count = 1;
			term.ternary.first = left;
			term.ternary.second = terms[0];
			term.ternary.third = terms[1];
		} else {
			read_error("Application with arity greater than 2 is unsupported", lparen_position);
			free(*left); free(left);
			for (hol_term* term : terms) {
				free(*term); if (term->reference_count == 0) free(term);
			}
			return false;
		}
	}
	return true;
}

template<hol_term_type OperatorType>
inline bool tptp_interpret_binary_term(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_term& term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	hol_term* left)
{
	hol_term* right;
	if (!new_hol_term(right)) {
		free(*left); free(left); return false;
	} else if (!tptp_interpret_unary_term(tokens, index, *right, names, variables)) {
		free(*left); free(left); free(right); return false;
	}

	term.binary.left = left;
	term.binary.right = right;
	term.type = OperatorType;
	term.reference_count = 1;
	return true;
}

template<hol_term_type OperatorType> struct tptp_operator_type { };
template<> struct tptp_operator_type<hol_term_type::AND> { static constexpr tptp_token_type type = tptp_token_type::AND; };
template<> struct tptp_operator_type<hol_term_type::OR> { static constexpr tptp_token_type type = tptp_token_type::OR; };

template<hol_term_type OperatorType>
bool tptp_interpret_binary_sequence(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_term& term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables,
	hol_term* left)
{
	array<hol_term*>& operands = *((array<hol_term*>*) alloca(sizeof(array<hol_term*>)));
	if (!array_init(operands, 8)) {
		free(*left); free(left); return false;
	}
	operands[0] = left;
	operands.length = 1;

	while (true) {
		hol_term* next;
		if (!new_hol_term(next) || !operands.ensure_capacity(operands.length + 1)
		 || !tptp_interpret_unary_term(tokens, index, *next, names, variables))
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

	term.type = OperatorType;
	term.reference_count = 1;
	term.array.operands = operands.data;
	term.array.length = operands.length;
	return true;
}

bool tptp_interpret(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_term& term,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	hol_term* left;
	if (!new_hol_term(left)) {
		fprintf(stderr, "tptp_interpret ERROR: Out of memory.\n");
		return false;
	} if (!tptp_interpret_unary_term(tokens, index, *left, names, variables)) {
		free(left); return false;
	}

	if (index >= tokens.length) {
		move(*left, term); free(left);
		return true;
	} else if (tokens[index].type == tptp_token_type::AND) {
		index++;
		if (!tptp_interpret_binary_sequence<hol_term_type::AND>(tokens, index, term, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::OR) {
		index++;
		if (!tptp_interpret_binary_sequence<hol_term_type::OR>(tokens, index, term, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::IF_THEN) {
		index++;
		if (!tptp_interpret_binary_term<hol_term_type::IF_THEN>(tokens, index, term, names, variables, left))
			return false;

	} else if (tokens[index].type == tptp_token_type::EQUALS) {
		index++;
		if (!tptp_interpret_binary_term<hol_term_type::EQUALS>(tokens, index, term, names, variables, left))
			return false;

	} else {
		move(*left, term); free(left);
	}
	return true;
}

template<typename Stream>
inline bool parse(Stream& in, hol_term& term,
		hash_map<string, unsigned int>& names,
		position start = position(1, 1))
{
	array<tptp_token> tokens = array<tptp_token>(128);
	if (!tptp_lex(tokens, in, start)) {
		read_error("Unable to parse higher-order formula (lexical analysis failed)", start);
		free_tokens(tokens); return false;
	}

	unsigned int index = 0;
	array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
	if (!tptp_interpret(tokens, index, term, names, variables)) {
		read_error("Unable to parse higher-order formula", start);
		for (auto entry : variables) free(entry.key);
		free_tokens(tokens); return false;
	}
	free_tokens(tokens);
	return true;
}

bool tptp_interpret_unary_type(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_type& type,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	if (index >= tokens.length) {
		fprintf(stderr, "ERROR: Unexpected end of input.\n");
		return false;
	} else if (tokens[index].type == tptp_token_type::LPAREN) {
		index++;
		if (!tptp_interpret_type(tokens, index, type, names, variables)
		 || !expect_token(tokens, index, tptp_token_type::RPAREN, "closing parenthesis in type expression"))
			return false;
		index++;
	} else if (tokens[index].type == tptp_token_type::IDENTIFIER) {
		if (tokens[index].text == "o" || tokens[index].text == "𝝄") {
			index++;
			if (!init(type, hol_constant_type::BOOLEAN)) return false;
		} else if (tokens[index].text == "i" || tokens[index].text == "𝜾") {
			index++;
			if (!init(type, hol_constant_type::INDIVIDUAL)) return false;
		} else if (tokens[index].text == "*") {
			index++;
			if (!init(type, hol_type_kind::ANY)) return false;
		} else {
			read_error("Expected a type expression", tokens[index].start);
			return false;
		}
	}
	return true;
}

bool tptp_interpret_type(
	const array<tptp_token>& tokens,
	unsigned int& index, hol_type& type,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	array<hol_type> types(8);
	while (true) {
		if (!types.ensure_capacity(types.length + 1)
		 || !tptp_interpret_unary_type(tokens, index, types[types.length], names, variables))
		{
			for (hol_type& type : types) free(type);
			return false;
		}
		types.length++;

		if (index < tokens.length && tokens[index].type == tptp_token_type::ARROW)
			index++;
		else break;
	}

	if (types.length == 0) {
		read_error("Expected a type expression", tokens.last().end);
		return false;
	} else if (types.length == 1) {
		move(types[0], type);
		return true;
	}

	if (!init(type, types.last())) {
		for (hol_type& type : types) free(type);
		return false;
	}
	for (unsigned int i = types.length - 1; i > 0; i--) {
		hol_type new_type(types[i - 1], type);
		swap(new_type, type);
	}
	for (hol_type& type : types) free(type);
	return true;
}

inline bool tptp_interpret(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hol_term& term, hol_type& type,
	hash_map<string, unsigned int>& names,
	array_map<string, unsigned int>& variables)
{
	if (!tptp_interpret(tokens, index, term, names, variables))
		return false;
	if (!expect_token(tokens, index, tptp_token_type::COLON, "colon in typing statement")) {
		free(term); return false;
	}
	index++;
	if (!tptp_interpret_type(tokens, index, type, names, variables)) {
		free(term); return false;
	}
	return true;
}

#endif /* HIGHER_ORDER_LOGIC_H_ */
