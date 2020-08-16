#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <math/multiset.h>
#include <stdexcept>

#if !defined(NDEBUG)
#include <functional>
#endif

#include "array_view.h"
#include "set_reasoning.h"
#include "built_in_predicates.h"

using namespace core;


template<typename Proof>
inline void visit_node(const Proof& proof) { }

template<bool Negated, typename Term>
constexpr bool visit_unary_atom(const Term* term) { return true; }

template<bool Negated>
constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2) { return true; }

template<typename Proof>
constexpr bool visit_subset_axiom(const Proof& proof) { return true; }

constexpr bool visit_existential_intro() { return true; }
constexpr bool visit_negated_conjunction() { return true; }
constexpr bool visit_negated_universal_intro() { return true; }
constexpr bool visit_disjunction_intro() { return true; }
inline void on_subtract_changes() { }

enum class instance_type {
	ANY,
	CONSTANT,
	INTEGER,
	STRING
};

struct instance {
	instance_type type;
	union {
		unsigned int constant;
		int64_t integer;
		string* str;
	};

	static inline void swap(instance& first, instance& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(instance); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void move(const instance& src, instance& dst) {
		dst.type = src.type;
		switch (src.type) {
		case instance_type::ANY: return;
		case instance_type::CONSTANT: dst.constant = src.constant; return;
		case instance_type::INTEGER: dst.integer = src.integer; return;
		case instance_type::STRING: dst.str = src.str; return;
		}
		fprintf(stderr, "instance.move ERROR: Unrecognized instance_type.\n");
		exit(EXIT_FAILURE);
	}
};


/* forward declarations */

struct instantiation;
bool intersect(instantiation&, const instantiation&, const instantiation&);

struct any_instantiation {
	instantiation* excluded;
	uint_fast8_t excluded_count;
};

struct any_integer_instantiation {
	pair<int64_t, int64_t>* included;
	uint_fast8_t included_count;

	inline int64_t min() const {
		return included[0].key;
	}

	inline int64_t max() const {
		return included[included_count - 1].value;
	}
};

inline bool init(any_integer_instantiation& any_integer, const any_integer_instantiation& src)
{
	any_integer.included_count = src.included_count;
	any_integer.included = (pair<int64_t, int64_t>*) malloc(max((size_t) 1, sizeof(pair<int64_t, int64_t>) * src.included_count));
	if (any_integer.included == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	for (uint_fast8_t i = 0; i < any_integer.included_count; i++)
		any_integer.included[i] = src.included[i];
	return true;
}

enum class instantiation_type {
	ANY,
	ANY_INTEGER,
	CONSTANT,
	INTEGER,
	STRING
};

struct instantiation {
	instantiation_type type;
	union {
		any_instantiation any;
		any_integer_instantiation any_integer;
		unsigned int constant;
		int64_t integer;
		const string* str;
	};

	instantiation(const instantiation& src) {
		if (!init_helper(src)) throw std::bad_alloc();
	}

	inline bool operator = (const instantiation& src) {
		return init_helper(src);
	}

	static inline void move(const instantiation& src, instantiation& dst) {
		dst.type = src.type;
		switch (src.type) {
		case instantiation_type::ANY:
			dst.any.excluded = src.any.excluded;
			dst.any.excluded_count = src.any.excluded_count;
			return;
		case instantiation_type::ANY_INTEGER:
			dst.any_integer.included = src.any_integer.included;
			dst.any_integer.included_count = src.any_integer.included_count;
			return;
		case instantiation_type::CONSTANT:
			dst.constant = src.constant; return;
		case instantiation_type::INTEGER:
			dst.integer = src.integer; return;
		case instantiation_type::STRING:
			dst.str = src.str; return;
		}
		fprintf(stderr, "instantiation.move ERROR: Unrecognized `instantiation_type`.\n");
		exit(EXIT_FAILURE);
	}

	static inline void free(instantiation& value) {
		if (value.type == instantiation_type::ANY) {
			for (uint_fast8_t i = 0; i < value.any.excluded_count; i++)
				core::free(value.any.excluded[i]);
			core::free(value.any.excluded);
		} else if (value.type == instantiation_type::ANY_INTEGER) {
			core::free(value.any_integer.included);
		}
	}

private:
	inline bool init_helper(const instantiation& src) {
		type = src.type;
		switch (type) {
		case instantiation_type::ANY:
			any.excluded_count = src.any.excluded_count;
			any.excluded = (instantiation*) malloc(max((size_t) 1, sizeof(instantiation) * src.any.excluded_count));
			if (any.excluded == nullptr) {
				fprintf(stderr, "instantiation.init_helper ERROR: Out of memory.\n");
				return false;
			}
			for (uint_fast8_t i = 0; i < any.excluded_count; i++) {
				if (!init(any.excluded[i], src.any.excluded[i])) {
					for (unsigned int j = 0; j < i; j++) core::free(any.excluded[j]);
					core::free(any.excluded);
					return false;
				}
			}
			return true;
		case instantiation_type::ANY_INTEGER:
			return init(any_integer, src.any_integer);
		case instantiation_type::CONSTANT:
			constant = src.constant; return true;
		case instantiation_type::INTEGER:
			integer = src.integer; return true;
		case instantiation_type::STRING:
			str = src.str; return true;
		}
		fprintf(stderr, "instantiation.init_helper ERROR: Unrecognized `instantiation_type`.\n");
		return false;
	}

	friend bool init(instantiation&, const instantiation&);
	friend bool operator != (const instantiation&, const instantiation&);
	friend bool operator < (const instantiation&, const instantiation&);
};

inline bool init(instantiation& inst, instantiation_type type) {
	inst.type = type;
	if (type == instantiation_type::ANY) {
		inst.any.excluded = (instantiation*) malloc(1);
		if (inst.any.excluded == nullptr)
			return false;
		inst.any.excluded_count = 0;
	}
	return true;
}

inline bool init(instantiation& inst, const instantiation& src) {
	return inst.init_helper(src);
}

inline bool operator == (const any_instantiation& first, const any_instantiation& second) {
	if (first.excluded_count != second.excluded_count)
		return false;
	for (uint_fast8_t i = 0; i < first.excluded_count; i++)
		if (first.excluded[i] != second.excluded[i]) return false;
	return true;
}

inline bool operator == (const any_integer_instantiation& first, const any_integer_instantiation& second) {
	if (first.included_count != second.included_count)
		return false;
	for (uint_fast8_t i = 0; i < first.included_count; i++)
		if (first.included[i] != second.included[i]) return false;
	return true;
}

inline bool operator < (const any_instantiation& first, const any_instantiation& second) {
	if (first.excluded_count < second.excluded_count) return true;
	else if (first.excluded_count > second.excluded_count) return false;
	for (uint_fast8_t i = 0; i < first.excluded_count; i++) {
		if (first.excluded[i] < second.excluded[i]) return true;
		else if (second.excluded[i] < first.excluded[i]) return false;
	}
	return false;
}

inline bool operator < (const any_integer_instantiation& first, const any_integer_instantiation& second) {
	if (first.included_count < second.included_count) return true;
	else if (first.included_count > second.included_count) return false;
	for (uint_fast8_t i = 0; i < first.included_count; i++) {
		if (first.included[i].key < second.included[i].key) return true;
		else if (first.included[i].key > second.included[i].key) return false;
		else if (first.included[i].value < second.included[i].value) return true;
		else if (first.included[i].value > second.included[i].value) return false;
	}
	return false;
}

inline bool operator == (const instantiation& first, const instantiation& second) {
	if (first.type != second.type) return false;
	switch (first.type) {
	case instantiation_type::ANY:
		return first.any == second.any;
	case instantiation_type::ANY_INTEGER:
		return first.any_integer == second.any_integer;
	case instantiation_type::CONSTANT:
		return first.constant == second.constant;
	case instantiation_type::INTEGER:
		return first.integer == second.integer;
	case instantiation_type::STRING:
		return first.str == second.str || *first.str == *second.str;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized `instantiation_type`.\n");
	return false;
}

inline bool operator != (const instantiation& first, const instantiation& second) {
	return !(first == second);
}

inline bool operator < (const instantiation& first, const instantiation& second) {
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case instantiation_type::ANY:
		return first.any < second.any;
	case instantiation_type::ANY_INTEGER:
		return first.any_integer < second.any_integer;
	case instantiation_type::CONSTANT:
		return first.constant < second.constant;
	case instantiation_type::INTEGER:
		return first.integer < second.integer;
	case instantiation_type::STRING:
		return *first.str < *second.str;
	}
	fprintf(stderr, "operator < ERROR: Unrecognized `instantiation_type`.\n");
	return false;
}

inline bool subtract_with_any(instantiation& dst,
	const instantiation& first, const instantiation& second)
{
	dst.any.excluded = (instantiation*) malloc(sizeof(instantiation) * (first.any.excluded_count + 1));
	if (dst.any.excluded == nullptr) {
		fprintf(stderr, "subtract_with_any ERROR: Insufficient memory for `instantiation.any.excluded`.\n");
		return false;
	}

	for (unsigned int i = 0; i < first.any.excluded_count; i++) {
		if (!init(dst.any.excluded[i], first.any.excluded[i])) {
			for (unsigned int j = 0; j < i; j++) free(dst.any.excluded[j]);
			free(dst.any.excluded);
			return false;
		}
	}

	if (!init(dst.any.excluded[first.any.excluded_count], second)) {
		for (unsigned int j = 0; j < first.any.excluded_count; j++) free(dst.any.excluded[j]);
		free(dst.any.excluded);
		return false;
	}
	dst.any.excluded_count = first.any.excluded_count + 1;
	dst.type = instantiation_type::ANY;
	return true;
}

inline bool set_union(any_integer_instantiation& dst,
		const any_integer_instantiation& first, int64_t second)
{
	dst.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>) * (first.included_count + 1));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_union ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;

	uint_fast8_t i = 1;
	while (i < first.included_count) {
		if (second + 1 < first.included[i].key) {
			dst.included[dst.included_count].key = second;
			dst.included[dst.included_count++].value = second;
			break;
		} else if (second + 1 == first.included[i].key) {
			dst.included[dst.included_count].key = second;
			dst.included[dst.included_count++].value = first.included[i].value;
			i++; break;
		} else if (second <= first.included[i].value) {
			dst.included[dst.included_count].key = first.included[i].key;
			dst.included[dst.included_count++].value = first.included[i].value;
			i++; break;
		} else if (second - 1 == first.included[i].value) {
			if (i + 1 < first.included_count && first.included[i + 1].key == second + 1) {
				dst.included[dst.included_count].key = first.included[i].key;
				dst.included[dst.included_count++].value = first.included[i + 1].value;
				i += 2; break;
			} else {
				dst.included[dst.included_count].key = first.included[i].key;
				dst.included[dst.included_count++].value = first.included[i].value;
				i++; break;
			}
		} else {
			dst.included[dst.included_count].key = first.included[i].key;
			dst.included[dst.included_count++].value = first.included[i].value;
			i++;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].key = first.included[i].key;
		dst.included[dst.included_count++].value = first.included[i].value;
		i++;
	}
	return true;
}

inline bool set_union(any_integer_instantiation& dst,
		const any_integer_instantiation& first,
		const any_integer_instantiation& second)
{
	int64_t lower, upper;
	uint_fast8_t i, j;
	if (first.included[0].key < second.included[0].key) {
		lower = first.included[0].key;
		upper = first.included[0].value;
		i = 1; j = 0;
	} else {
		lower = second.included[0].key;
		upper = second.included[0].value;
		i = 0; j = 1;
	}

	dst.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_union ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	while (i < first.included_count && j < second.included_count) {
		if (first.included[i].key <= upper + 1) {
			upper = max(upper, first.included[i++].value);
		} else if (second.included[j].key <= upper + 1) {
			upper = max(upper, second.included[j++].value);
		} else {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = upper;
			if (first.included[i].key < second.included[j].key) {
				lower = first.included[i].key;
				upper = first.included[i++].value;
			} else {
				lower = second.included[j].key;
				upper = second.included[j++].value;
			}
		}
	}

	while (i < first.included_count) {
		if (first.included[i].key <= upper + 1) {
			upper = max(upper, first.included[i++].value);
		} else {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = upper;
			lower = first.included[i].key;
			upper = first.included[i++].value;
		}
	} while (j < second.included_count) {
		if (second.included[j].key <= upper + 1) {
			upper = max(upper, second.included[j++].value);
		} else {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = upper;
			lower = second.included[j].key;
			upper = second.included[j++].value;
		}
	}
	return true;
}

inline bool set_subtract(any_integer_instantiation& dst,
		const any_integer_instantiation& first,
		const any_integer_instantiation& second)
{
	dst.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_subtract ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0, j = 0;
	int64_t lower = first.included[0].key;
	while (i < first.included_count && j < second.included_count) {
		if (second.included[j].value < lower) {
			j++;
		} else if (second.included[j].key <= lower) {
			if (second.included[j].value < first.included[i].value) {
				lower = second.included[j].value + 1;
			} else if (i + 1 == first.included_count) {
				i++; break;
			} else {
				lower = first.included[++i].key;
			}
		} else if (second.included[j].key <= first.included[i].value) {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = second.included[j].key - 1;
			if (second.included[j].value < first.included[i].value) {
				lower = second.included[j++].value + 1;
			} else if (i + 1 == first.included_count) {
				i++; break;
			} else {
				lower = first.included[++i].key;
			}
		} else {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = first.included[i++].value;
			if (i == first.included_count) break;
			lower = first.included[i].key;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].key = lower;
		dst.included[dst.included_count++].value = first.included[i++].value;
		if (i < first.included_count)
			lower = first.included[i].key;
	}
	return true;
}

inline bool set_subtract(any_integer_instantiation& dst,
		const any_integer_instantiation& first, int64_t second)
{
	dst.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>) * (first.included_count + 1));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_subtract ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0;
	while (i < first.included_count) {
		if (first.included[i].value < second) {
			dst.included[dst.included_count].key = first.included[i].key;
			dst.included[dst.included_count++].value = first.included[i++].value;
		} else if (first.included[i].value == second) {
			if (first.included[i].key < first.included[i].value) {
				dst.included[dst.included_count].key = first.included[i].key;
				dst.included[dst.included_count++].value = second - 1;
			}
			i++; break;
		} else if (second > first.included[i].key) {
			dst.included[dst.included_count].key = first.included[i].key;
			dst.included[dst.included_count++].value = second - 1;
			dst.included[dst.included_count].key = second + 1;
			dst.included[dst.included_count++].value = first.included[i].value;
			i++; break;
		} else if (second == first.included[i].key) {
			if (first.included[i].key < first.included[i].value) {
				dst.included[dst.included_count].key = second + 1;
				dst.included[dst.included_count++].value = first.included[i].value;
			}
			i++; break;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].key = first.included[i].key;
		dst.included[dst.included_count++].value = first.included[i++].value;
	}
	return true;
}

inline bool set_intersect(any_integer_instantiation& dst,
		const any_integer_instantiation& first,
		const any_integer_instantiation& second)
{
	dst.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_union ERROR: Insufficient memory for `any_integer_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0, j = 0;
	int64_t lower = max(first.included[0].key, second.included[0].key);
	while (i < first.included_count && j < second.included_count) {
		if (first.included[i].value < lower) {
			i++;
		} else if (second.included[j].value < lower) {
			j++;
		} else {
			dst.included[dst.included_count].key = lower;
			dst.included[dst.included_count++].value = min(first.included[i++].value, second.included[j++].value);
			if (i == first.included_count || j == second.included_count)
				return true;
			lower = max(first.included[i].key, second.included[j].key);
		}
	}
	return true;
}

inline bool subtract(instantiation& dst,
		const instantiation& first, const instantiation& second)
{
	if (first.type == instantiation_type::ANY) {
		return subtract_with_any(dst, first, second);
	} else if (second.type == instantiation_type::ANY) {
		if (!init(dst, first)) return false;
		for (uint_fast8_t i = 0; i < second.any.excluded_count; i++) {
			instantiation& tmp = *((instantiation*) alloca(sizeof(instantiation)));
			if (!intersect(tmp, first, second.any.excluded[i])) {
				free(dst);
				return false;
			}

			/* compute the union of `dst` and `tmp` */
			if (dst.type == instantiation_type::ANY_INTEGER) {
				instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
				if (tmp.type == instantiation_type::ANY_INTEGER) {
					dummy.type = instantiation_type::ANY_INTEGER;
					if (!set_union(dummy.any_integer, dst.any_integer, tmp.any_integer)) return false;
				} else if (tmp.type == instantiation_type::INTEGER) {
					dummy.type = instantiation_type::ANY_INTEGER;
					if (!set_union(dummy.any_integer, dst.any_integer, tmp.integer)) return false;
				} else {
					fprintf(stderr, "subtract ERROR: Unclosed subtraction of `instantiation` objects.\n");
					return false;
				}
				free(dst); free(tmp);
				move(dummy, dst);
			} else if (tmp.type == instantiation_type::ANY_INTEGER) {
				instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
				if (dst.type == instantiation_type::INTEGER) {
					dummy.type = instantiation_type::ANY_INTEGER;
					if (!set_union(dummy.any_integer, tmp.any_integer, dst.integer)) return false;
				} else {
					fprintf(stderr, "subtract ERROR: Unclosed subtraction of `instantiation` objects.\n");
					return false;
				}
				free(dst); free(tmp);
				move(dummy, dst);
			} else {
				free(tmp);
			}
		}
	} else if (first.type == instantiation_type::ANY_INTEGER) {
		if (second.type == instantiation_type::ANY_INTEGER) {
			if (!set_subtract(dst.any_integer, first.any_integer, second.any_integer)) return false;
		} else if (second.type == instantiation_type::INTEGER) {
			if (!set_subtract(dst.any_integer, first.any_integer, second.integer)) return false;
		} else {
			return init(dst, first);
		}

		if (dst.any_integer.included_count == 1 && dst.any_integer.included[0].key == dst.any_integer.included[0].value) {
			int64_t value = dst.any_integer.included[0].key;
			free(dst.any_integer.included);
			dst.integer = value;
			dst.type = instantiation_type::INTEGER;
		} else {
			dst.type = instantiation_type::ANY_INTEGER;
		}
		return true;
	} else if (second.type == instantiation_type::ANY_INTEGER) {
		if (first.type == instantiation_type::INTEGER) {
			for (uint_fast8_t i = 0; i < second.any_integer.included_count; i++)
				if (first.integer >= second.any_integer.included[i].key && first.integer <= second.any_integer.included[i].value) return false;
			return init(dst, first);
		} else {
			return init(dst, first);
		}
	} else if (first.type != second.type) {
		return init(dst, first);
	}

	switch (first.type) {
	case instantiation_type::CONSTANT:
		if (first.constant == second.constant) return false;
		dst.type = instantiation_type::CONSTANT;
		dst.constant = first.constant;
		return true;
	case instantiation_type::INTEGER:
		if (first.integer == second.integer) return false;
		dst.type = instantiation_type::INTEGER;
		dst.integer = first.integer;
		return true;
	case instantiation_type::STRING:
		if (first.str == second.str || *first.str == *second.str) return false;
		dst.type = instantiation_type::STRING;
		dst.str = first.str;
		return true;
	case instantiation_type::ANY:
	case instantiation_type::ANY_INTEGER:
		break;
	}
	fprintf(stderr, "intersect ERROR: Unrecognized `instantiation_type`.\n");
	return false;
}

inline bool intersect_with_any(instantiation& dst,
	const instantiation& first, const instantiation& second)
{
	if (second.type == instantiation_type::ANY) {
		array<instantiation>& new_excluded = *((array<instantiation>*) alloca(sizeof(array<instantiation>)));
		if  (!array_init(new_excluded, first.any.excluded_count + second.any.excluded_count))
			return false;
		set_union(new_excluded.data, new_excluded.length, first.any.excluded, first.any.excluded_count, second.any.excluded, second.any.excluded_count);
		dst.type = instantiation_type::ANY;
		dst.any.excluded = new_excluded.data;
		dst.any.excluded_count = new_excluded.length;
		return true;
	} else if (second.type == instantiation_type::ANY_INTEGER) {
		if (!subtract(dst, second, first.any.excluded[0]))
			return false;
		if (dst.type == instantiation_type::INTEGER) {
			int64_t value = dst.integer;
			dst.any_integer.included = (pair<int64_t, int64_t>*) malloc(sizeof(pair<int64_t, int64_t>));
			if (dst.any_integer.included == nullptr) {
				fprintf(stderr, "intersect_with_any ERROR: Out of memory.\n");
				return false;
			}
			dst.any_integer.included[0] = {value, value};
			dst.any_integer.included_count = 1;
			dst.type = instantiation_type::ANY_INTEGER;
		}

		for (uint_fast8_t i = 1; i < first.any.excluded_count; i++) {
			instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
			if (!subtract(temp, second, first.any.excluded[i])) {
				free(dst);
				return false;
			}

			instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
			if (temp.type == instantiation_type::ANY_INTEGER) {
				if (!set_union(dummy.any_integer, dst.any_integer, temp.any_integer)) {
					free(dst); free(temp);
					return false;
				}
			} else if (temp.type == instantiation_type::INTEGER) {
				if (!set_union(dummy.any_integer, dst.any_integer, temp.integer)) {
					free(dst); free(temp);
					return false;
				}
			}
			dummy.type = instantiation_type::ANY_INTEGER;
			free(dst); free(temp);
			move(dummy, dst);
		}

		if (dst.any_integer.included_count == 1 && dst.any_integer.included[0].key == dst.any_integer.included[0].value) {
			int64_t value = dst.any_integer.included[0].value;
			free(dst.any_integer.included);
			dst.integer = value;
			dst.type = instantiation_type::INTEGER;
		}
		return true;
	}

	for (uint_fast8_t i = 0; i < first.any.excluded_count; i++) {
		if (intersect(dst, first.any.excluded[i], second)) {
			free(dst);
			return false;
		}
	}

	return init(dst, second);
}

inline bool intersect_with_any_integer(instantiation& dst,
	const instantiation& first, const instantiation& second)
{
	if (second.type == instantiation_type::ANY_INTEGER) {
		if (!set_intersect(dst.any_integer, first.any_integer, second.any_integer))
			return false;
		if (dst.any_integer.included_count == 1 && dst.any_integer.included[0].key == dst.any_integer.included[0].value) {
			int64_t value = dst.any_integer.included[0].value;
			free(dst.any_integer.included);
			dst.integer = value;
			dst.type = instantiation_type::INTEGER;
		} else {
			dst.type = instantiation_type::ANY_INTEGER;
		}
		return true;
	} else if (second.type == instantiation_type::INTEGER) {
		for (uint_fast8_t i = 0; i < first.any_integer.included_count; i++)
			if (second.integer >= first.any_integer.included[i].key && second.integer <= first.any_integer.included[i].value) return init(dst, second);
		return false;
	} else {
		return false;
	}
}

inline bool intersect(instantiation& dst,
	const instantiation& first, const instantiation& second)
{
	if (first.type == instantiation_type::ANY) {
		return intersect_with_any(dst, first, second);
	} else if (second.type == instantiation_type::ANY) {
		return intersect_with_any(dst, second, first);
	} else if (first.type == instantiation_type::ANY_INTEGER) {
		return intersect_with_any_integer(dst, first, second);
	} else if (second.type == instantiation_type::ANY_INTEGER) {
		return intersect_with_any_integer(dst, second, first);
	} else if (first.type != second.type) {
		return false;
	}

	switch (first.type) {
	case instantiation_type::CONSTANT:
		if (first.constant != second.constant) return false;
		dst.type = instantiation_type::CONSTANT;
		dst.constant = first.constant;
		return true;
	case instantiation_type::INTEGER:
		if (first.integer != second.integer) return false;
		dst.type = instantiation_type::INTEGER;
		dst.integer = first.integer;
		return true;
	case instantiation_type::STRING:
		if (first.str != second.str && *first.str != *second.str) return false;
		dst.type = instantiation_type::STRING;
		dst.str = first.str;
		return true;
	case instantiation_type::ANY:
	case instantiation_type::ANY_INTEGER:
		break;
	}
	fprintf(stderr, "intersect ERROR: Unrecognized `instantiation_type`.\n");
	return false;
}

struct pair_sorter { };

template<typename K, typename V>
inline bool less_than(const pair<K, V>& first, const pair<K, V>& second, const pair_sorter& sorter) {
	if (first.key < second.key) return true;
	else if (first.key > second.key) return false;
	else if (first.value < second.value) return true;
	else return false;
}

struct instantiation_tuple {
	instantiation* values;
	uint_fast8_t length;
	array<pair<uint_fast8_t, uint_fast8_t>> equal_indices;
	array<pair<uint_fast8_t, uint_fast8_t>> not_equal_indices;
	array<pair<uint_fast8_t, uint_fast8_t>> ge_indices;

	instantiation_tuple(const instantiation_tuple& src) :
			equal_indices(max((size_t) 1, src.equal_indices.length)),
			not_equal_indices(max((size_t) 1, src.not_equal_indices.length)),
			ge_indices(max((size_t) 1, src.ge_indices.length))
	{
		if (!init_helper(src)) throw std::bad_alloc();
	}

	~instantiation_tuple() { free_helper(); }

	inline bool operator = (const instantiation_tuple& src) {
		if (!array_init(equal_indices, max((size_t) 1, src.equal_indices.length))) {
			return false;
		} else if (!array_init(not_equal_indices, max((size_t) 1, src.not_equal_indices.length))) {
			core::free(equal_indices);
			return false;
		} else if (!array_init(ge_indices, max((size_t) 1, src.ge_indices.length))) {
			core::free(equal_indices);
			core::free(not_equal_indices);
			return false;
		} else if (!init_helper(src)) {
			core::free(equal_indices);
			core::free(not_equal_indices);
			core::free(ge_indices);
			return false;
		}
		return true;
	}

	bool change_value(uint_fast8_t index, instantiation& new_value) {
		if (values[index].type == instantiation_type::ANY || values[index].type == instantiation_type::ANY_INTEGER) {
			if (new_value.type == instantiation_type::ANY || new_value.type == instantiation_type::ANY_INTEGER) {
				uint_fast8_t root = length;
				for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
					if (entry.key == index || entry.value == index) {
						root = entry.key;
						break;
					}
				}
				if (root == length) root = index;
				array<uint_fast8_t> component(8);
				component[component.length++] = root;
				for (unsigned int i = 0; i < equal_indices.length; i++) {
					if (equal_indices[i].key == root) {
						if (!component.add(equal_indices[i].value))
							return false;
					}
				}
				for (uint_fast8_t other : component) {
					core::free(values[other]);
					if (!init(values[other], new_value))
						return false;
				}

				/* check the new value can satisfy the greater-than-or-equal-to graph */
				if (new_value.type == instantiation_type::ANY_INTEGER) {
					for (const pair<uint_fast8_t, uint_fast8_t>& pair : ge_indices) {
						if (pair.key == root || pair.value == root) {
							if (values[pair.key].any_integer.max() == values[pair.value].any_integer.min()) {
								instantiation& new_value = *((instantiation*) alloca(sizeof(instantiation)));
								new_value.type = instantiation_type::INTEGER;
								new_value.integer = values[pair.value].any_integer.min();
								change_value(index, new_value);
								return true;
							} else if (values[pair.key].any_integer.max() < values[pair.value].any_integer.min()) {
								return false;
							}
						}
					}
				}
			} else {
				uint_fast8_t root = length;
				for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
					if (entry.key == index || entry.value == index) {
						root = entry.key;
						break;
					}
				}
				if (root == length) root = index;
				array<uint_fast8_t> component(4);
				component[component.length++] = root;
				for (unsigned int i = 0; i < equal_indices.length; i++) {
					if (equal_indices[i].key == root) {
						if (!component.add(equal_indices[i].value))
							return false;
						equal_indices.remove(i--);
					}
				}
				for (uint_fast8_t other : component) {
					core::free(values[other]);
					if (!init(values[other], new_value))
						return false;
				}
				if (equal_indices.length > 1)
					insertion_sort(equal_indices, pair_sorter());

				for (unsigned int i = 0; i < not_equal_indices.length; i++) {
					uint_fast8_t other;
					if (component.contains(not_equal_indices[i].key))
						other = not_equal_indices[i].value;
					else if (component.contains(not_equal_indices[i].value))
						other = not_equal_indices[i].key;
					else continue;

					instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
					if (!subtract(dummy, values[other], new_value)
					 || !change_value(other, dummy))
					{
						core::free(dummy);
						return false;
					}
					core::free(dummy);
					not_equal_indices.remove(i--);
				}
				if (not_equal_indices.length > 1)
					insertion_sort(not_equal_indices, pair_sorter());

				if (values[index].type == instantiation_type::ANY_INTEGER) {
					for (unsigned int i = 0; i < ge_indices.length; i++) {
						if (ge_indices[i].key == index) {
							instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
							dummy.any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
							dummy.any_integer.included[0] = {INT64_MIN, new_value.integer};
							dummy.any_integer.included_count = 1;
							instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
							if (!intersect(temp, values[ge_indices[i].value], dummy)
							 || !change_value(ge_indices[i].value, temp))
								return false;
							core::free(temp);
							ge_indices.remove(i--);
						} else if (ge_indices[i].value == index) {
							instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
							dummy.any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
							dummy.any_integer.included[0] = {new_value.integer, INT64_MAX};
							dummy.any_integer.included_count = 1;
							instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
							if (!intersect(temp, values[ge_indices[i].key], dummy)
							 || !change_value(ge_indices[i].key, temp))
								return false;
							core::free(temp);
							ge_indices.remove(i--);
						}
					}
				}
			}
		} else {
			core::free(values[index]);
			return init(new_value, values[index]);
		}
		return true;
	}

	template<typename Term>
	bool unify_value(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (term->type == TermType::CONSTANT) {
			dummy[0].type = instantiation_type::CONSTANT;
			dummy[0].constant = term->constant;
			if (!intersect(dummy[1], values[index], dummy[0])) return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::INTEGER) {
			dummy[0].type = instantiation_type::INTEGER;
			dummy[0].integer = term->integer;
			if (!intersect(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::STRING) {
			dummy[0].type = instantiation_type::STRING;
			dummy[0].str = &term->str;
			if (!intersect(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::VARIABLE) {
			if (!intersect(dummy[1], values[index], values[term->variable - 1]))
				return false;

			if (dummy[1].type != instantiation_type::ANY && dummy[1].type != instantiation_type::ANY_INTEGER) {
				bool result = change_value(index, dummy[1]) && change_value(term->variable - 1, dummy[1]);
				core::free(dummy[1]); return result;
			}

			uint_fast8_t first_root = length;
			uint_fast8_t second_root = length;
			for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
				if (entry.key == index || entry.value == index) {
					first_root = entry.key;
					break;
				}
			} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
				if (entry.key == term->variable - 1 || entry.value == term->variable - 1) {
					second_root = entry.key;
					break;
				}
			}

			if (first_root == second_root) return true;

			/* check if a cycle is created in the greater-than-or-equal graph */
			array<bool> new_component(length);
			for (uint_fast8_t i = 0; i < length; i++)
				new_component[i] = false;
			if (dummy[1].type == instantiation_type::ANY_INTEGER) {
				array<bool> ancestors(length), descendants(length);
				for (uint_fast8_t i = 0; i < length; i++) {
					ancestors[i] = false;
					descendants[i] = false;
				}
				ancestors[first_root] = true;
				descendants[second_root] = true;
				bool changed = false;
				do {
					for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
						if (ancestors[edge.value] && !ancestors[edge.key]) {
							ancestors[edge.key] = true;
							changed = true;
						}
					}
				} while (changed);
				do {
					for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
						if (descendants[edge.key] && !descendants[edge.value]) {
							descendants[edge.value] = true;
							changed = true;
						}
					}
				} while (changed);

				for (uint_fast8_t i = 0; i < length; i++)
					new_component[i] = (ancestors[i] && descendants[i]);

				for (uint_fast8_t i = 0; i < length; i++) {
					ancestors[i] = false;
					descendants[i] = false;
				}
				ancestors[second_root] = true;
				descendants[first_root] = true;
				do {
					for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
						if (ancestors[edge.value] && !ancestors[edge.key]) {
							ancestors[edge.key] = true;
							changed = true;
						}
					}
				} while (changed);
				do {
					for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
						if (descendants[edge.key] && !descendants[edge.value]) {
							descendants[edge.value] = true;
							changed = true;
						}
					}
				} while (changed);

				for (uint_fast8_t i = 0; i < length; i++)
					new_component[i] = (new_component[i] || (ancestors[i] && descendants[i]));
			}
			new_component[first_root] = true;
			new_component[second_root] = true;

			uint_fast8_t new_root = length;
			for (uint_fast8_t i = 0; i < length; i++)
				if (new_component[i]) { new_root = i; break; }

			for (uint_fast8_t i = 0; i < length; i++) {
				if (!new_component[i] || i == index || i == term->variable - 1) continue;
				if (!intersect(dummy[0], dummy[1], values[i])) {
					core::free(dummy[1]);
					return false;
				}
				core::free(dummy[1]);
				core::move(dummy[0], dummy[1]);
			}

			for (uint_fast8_t i = 0; i < equal_indices.length; i++) {
				if (new_component[equal_indices[i].key])
					equal_indices[i].key = new_root;
				if (equal_indices[i].key == equal_indices[i].value)
					equal_indices.remove(i--);
			}
			insertion_sort(equal_indices, pair_sorter());

			for (uint_fast8_t i = 0; i < not_equal_indices.length; i++) {
				if (new_component[not_equal_indices[i].key])
					not_equal_indices[i].key = new_root;
				if (not_equal_indices[i].key == not_equal_indices[i].value) {
					core::free(dummy[1]);
					return false;
				}
			}
			insertion_sort(not_equal_indices, pair_sorter());
			unique(not_equal_indices);

			for (uint_fast8_t i = 0; i < length; i++) {
				if (!new_component[i]) continue;
				core::free(values[i]);
				if (!init(values[i], dummy[1])) {
					core::free(dummy[1]);
					return false;
				}
			} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
				if (entry.key == new_root) {
					core::free(values[entry.value]);
					if (!init(values[entry.value], dummy[1])) {
						core::free(dummy[1]);
						return false;
					}
				}
			}
			core::free(dummy[1]);
			return true;

		} else {
			return false;
		}
	}

	template<typename Term>
	bool antiunify_value(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (term->type == TermType::CONSTANT) {
			dummy[0].type = instantiation_type::CONSTANT;
			dummy[0].constant = term->constant;
			if (!subtract(dummy[1], values[index], dummy[0])) return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::INTEGER) {
			dummy[0].type = instantiation_type::INTEGER;
			dummy[0].integer = term->integer;
			if (!subtract(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::STRING) {
			dummy[0].type = instantiation_type::STRING;
			dummy[0].str = &term->str;
			if (!subtract(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::VARIABLE) {
			if (values[index].type == instantiation_type::ANY || values[index].type == instantiation_type::ANY_INTEGER) {
				if (values[term->variable - 1].type == instantiation_type::ANY || values[term->variable - 1].type == instantiation_type::ANY_INTEGER) {
					uint_fast8_t first_root = length;
					uint_fast8_t second_root = length;
					for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
						if (entry.key == index || entry.value == index) {
							first_root = entry.key;
							break;
						}
					} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
						if (entry.key == term->variable - 1 || entry.value == term->variable - 1) {
							second_root = entry.key;
							break;
						}
					}
					if (first_root == length) first_root = index;
					if (second_root == length) second_root = term->variable - 1;
					if (first_root == second_root) return false;

					if (!not_equal_indices.ensure_capacity(not_equal_indices.length + 1))
						return false;

					pair<uint_fast8_t, uint_fast8_t> new_pair = (first_root < second_root ? make_pair(first_root, second_root) : make_pair(second_root, first_root));
					uint_fast8_t index = 0;
					while (index < not_equal_indices.length && less_than(new_pair, not_equal_indices[index], pair_sorter()))
						index++;
					if (index < not_equal_indices.length && new_pair == not_equal_indices[index])
						return true;
					shift_right(not_equal_indices.data, not_equal_indices.length, index);
					not_equal_indices[index] = new_pair;
					not_equal_indices.length++;
					return true;
				} else {
					if (!subtract(dummy[1], values[index], values[term->variable - 1]))
						return false;
					bool result = change_value(index, dummy[1]);
					core::free(dummy[1]); return result;
				}
			} else {
				if (values[term->variable - 1].type == instantiation_type::ANY) {
					if (!subtract(dummy[1], values[term->variable - 1], values[index]))
						return false;
					bool result = change_value(term->variable - 1, dummy[1]);
					core::free(dummy[1]); return result;
				} else {
					return values[index] != values[term->variable - 1];
				}
			}
		} else {
			return false;
		}
	}

	inline bool unify_greater_than_or_equal(uint_fast8_t first, uint_fast8_t second)
	{
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (values[first].type != instantiation_type::ANY_INTEGER) {
			dummy[0].any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
			dummy[0].any_integer.included[0] = {INT64_MIN, INT64_MAX};
			dummy[0].any_integer.included_count = 1;
			if (!intersect(dummy[1], values[first], dummy[0]))
				return false;
			bool result = change_value(first, dummy[1]);
			core::free(dummy[1]);
			if (!result) return false;
		} if (values[second].type != instantiation_type::ANY_INTEGER) {
			dummy[0].any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
			dummy[0].any_integer.included[0] = {INT64_MIN, INT64_MAX};
			dummy[0].any_integer.included_count = 1;
			if (!intersect(dummy[1], values[second], dummy[0]))
				return false;
			bool result = change_value(second, dummy[1]);
			core::free(dummy[1]);
			if (!result) return false;
		}

		uint_fast8_t first_root = length;
		uint_fast8_t second_root = length;
		for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
			if (entry.key == first || entry.value == first) {
				first_root = entry.key;
				break;
			}
		} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
			if (entry.key == second || entry.value == second) {
				second_root = entry.key;
				break;
			}
		}

		if (first_root == second_root) return true;

		/* make sure there isn't already a path from `first_root` to `second_root` */
		array<bool> descendants(length);
		for (uint_fast8_t i = 0; i < length; i++)
			descendants[i] = false;
		descendants[first_root] = true;
		bool changed = false;
		do {
			for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
				if (descendants[edge.key] && !descendants[edge.value]) {
					descendants[edge.value] = true;
					changed = true;
				}
			}
		} while (changed);

		if (descendants[second_root])
			/* there already is a path from `first_root` to `second_root` */
			return true;

		array<bool> ancestors(length);
		for (uint_fast8_t i = 0; i < length; i++) {
			ancestors[i] = false;
			descendants[i] = false;
		}
		ancestors[first_root] = true;
		descendants[second_root] = true;
		do {
			for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
				if (ancestors[edge.value] && !ancestors[edge.key]) {
					ancestors[edge.key] = true;
					changed = true;
				}
			}
		} while (changed);
		do {
			for (const pair<uint_fast8_t, uint_fast8_t>& edge : ge_indices) {
				if (descendants[edge.key] && !descendants[edge.value]) {
					descendants[edge.value] = true;
					changed = true;
				}
			}
		} while (changed);

		array<bool> new_component(length);
		for (uint_fast8_t i = 0; i < length; i++)
			new_component[i] = (ancestors[i] && descendants[i]);

		uint_fast8_t new_root = length;
		for (uint_fast8_t i = 0; i < length; i++)
			if (new_component[i]) { new_root = i; break; }
		if (new_root != length) {
			instantiation& new_value = *((instantiation*) alloca(sizeof(instantiation)));
			if (!init(new_value, values[new_root])) return false;
			for (uint_fast8_t i = new_root + 1; i < length; i++) {
				if (!new_component[i]) continue;
				instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
				if (!intersect(temp, new_value, values[i])) {
					core::free(new_value);
					return false;
				}
				core::free(new_value);
				core::move(temp, new_value);
			}

			for (uint_fast8_t i = 0; i < equal_indices.length; i++) {
				if (new_component[equal_indices[i].key])
					equal_indices[i].key = new_root;
				if (equal_indices[i].key == equal_indices[i].value)
					equal_indices.remove(i--);
			}
			insertion_sort(equal_indices, pair_sorter());

			for (uint_fast8_t i = 0; i < not_equal_indices.length; i++) {
				if (new_component[not_equal_indices[i].key])
					not_equal_indices[i].key = new_root;
				if (not_equal_indices[i].key == not_equal_indices[i].value) {
					core::free(new_value);
					return false;
				}
			}
			insertion_sort(not_equal_indices, pair_sorter());
			unique(not_equal_indices);

			for (uint_fast8_t i = 0; i < length; i++) {
				if (!new_component[i]) continue;
				core::free(values[i]);
				if (!init(values[i], new_value)) {
					core::free(new_value);
					return false;
				}
			} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
				if (entry.key == new_root) {
					core::free(values[entry.value]);
					if (!init(values[entry.value], new_value)) {
						core::free(new_value);
						return false;
					}
				}
			}
			core::free(new_value);
		} else {
			/* add the actual edge to `ge_indices` */
			for (uint_fast8_t i = 0; i < ge_indices.length; i++) {
				if (ancestors.contains(ge_indices[i].key) && descendants.contains(ge_indices[i].value))
					ge_indices.remove(i--);
			}
			if (!ge_indices.ensure_capacity(ge_indices.length + 1))
				return false;
			ge_indices[ge_indices.length++] = {first_root, second_root};
			insertion_sort(ge_indices);
		}
		return true;
	}

	template<typename Term>
	bool unify_greater_than_or_equal(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (term->type == TermType::INTEGER) {
			dummy[0].any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
			dummy[0].any_integer.included[0] = {term->integer, INT64_MAX};
			dummy[0].any_integer.included_count = 1;
			if (!intersect(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::VARIABLE) {
			return unify_greater_than_or_equal(index, term->variable - 1);
		} else {
			return false;
		}
	}

	template<typename Term>
	bool unify_less_than_or_equal(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (term->type == TermType::INTEGER) {
			dummy[0].any_integer.included = (pair<int64_t, int64_t>*) alloca(sizeof(pair<int64_t, int64_t>));
			dummy[0].any_integer.included[0] = {INT64_MIN, term->integer};
			dummy[0].any_integer.included_count = 1;
			if (!intersect(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::VARIABLE) {
			return unify_greater_than_or_equal(term->variable - 1, index);
		} else {
			return false;
		}
	}

	static inline void free(instantiation_tuple& tuple) {
		tuple.free_helper();
		core::free(tuple.equal_indices);
		core::free(tuple.not_equal_indices);
		core::free(tuple.ge_indices);
	}

	static inline void move(const instantiation_tuple& src, instantiation_tuple& dst) {
		dst.values = src.values;
		dst.length = src.length;
		core::move(src.equal_indices, dst.equal_indices);
		core::move(src.not_equal_indices, dst.not_equal_indices);
		core::move(src.ge_indices, dst.ge_indices);
	}

	static inline void swap(instantiation_tuple& first, instantiation_tuple& second) {
		core::swap(first.values, second.values);
		core::swap(first.length, second.length);
		core::swap(first.equal_indices, second.equal_indices);
		core::swap(first.not_equal_indices, second.not_equal_indices);
		core::swap(first.ge_indices, second.ge_indices);
	}

private:
	inline bool init_helper(uint_fast8_t src_length) {
		length = src_length;
		values = (instantiation*) malloc(sizeof(instantiation) * length);
		if (values == nullptr) {
			fprintf(stderr, "instantiation_tuple.init_helper ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < length; i++) {
			if (!init(values[i], instantiation_type::ANY)) {
				for (unsigned int j = 0; j < i; j++) core::free(values[j]);
				core::free(values);
				return false;
			}
		}
		return true;
	}

	inline bool init_helper(const instantiation_tuple& src) {
		length = src.length;
		values = (instantiation*) malloc(sizeof(instantiation) * length);
		if (values == nullptr) {
			fprintf(stderr, "instantiation_tuple.init_helper ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < length; i++) {
			if (!init(values[i], src.values[i])) {
				for (unsigned int j = 0; j < i; j++) core::free(values[j]);
				core::free(values);
				return false;
			}
		}

		if (!equal_indices.append(src.equal_indices.data, src.equal_indices.length)
		 || !not_equal_indices.append(src.not_equal_indices.data, src.not_equal_indices.length)
		 || !ge_indices.append(src.ge_indices.data, src.ge_indices.length))
		{
			free_helper();
			return false;
		}
		return true;
	}

	inline void free_helper() {
		for (unsigned int i = 0; i < length; i++)
			core::free(values[i]);
		core::free(values);
	}

	friend bool init(instantiation_tuple&, uint_fast8_t);
	friend bool init(instantiation_tuple&, const instantiation_tuple&);
};

inline bool init(instantiation_tuple& new_tuple, uint_fast8_t length) {
	if (!array_init(new_tuple.equal_indices, 2)) {
		return false;
	} else if (!array_init(new_tuple.not_equal_indices, 2)) {
		free(new_tuple.equal_indices);
		return false;
	} else if (!array_init(new_tuple.ge_indices, 2)) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		return false;
	} else if (!new_tuple.init_helper(length)) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		free(new_tuple.ge_indices);
		return false;
	}
	return true;
}

inline bool init(instantiation_tuple& new_tuple, const instantiation_tuple& src) {
	if (!array_init(new_tuple.equal_indices, max((size_t) 1, src.equal_indices.length))) {
		return false;
	} else if (!array_init(new_tuple.not_equal_indices, max((size_t) 1, src.not_equal_indices.length))) {
		free(new_tuple.equal_indices);
		return false;
	} else if (!array_init(new_tuple.ge_indices, max((size_t) 1, src.ge_indices.length))) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		return false;
	} else if (!new_tuple.init_helper(src)) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		free(new_tuple.ge_indices);
		return false;
	}
	return true;
}

inline bool operator == (const instantiation_tuple& src, const instantiation_tuple& dst) {
	if (src.length != dst.length
	 || src.equal_indices.length != dst.equal_indices.length
	 || src.not_equal_indices.length != dst.not_equal_indices.length
	 || src.ge_indices.length != dst.ge_indices.length)
		return false;
	for (unsigned int i = 0; i < src.length; i++)
		if (src.values[i] != dst.values[i]) return false;
	for (unsigned int i = 0; i < src.equal_indices.length; i++)
		if (src.equal_indices[i] != dst.equal_indices[i]) return false;
	for (unsigned int i = 0; i < src.not_equal_indices.length; i++)
		if (src.not_equal_indices[i] != dst.not_equal_indices[i]) return false;
	for (unsigned int i = 0; i < src.ge_indices.length; i++)
		if (src.ge_indices[i] != dst.ge_indices[i]) return false;
	return true;
}

inline bool operator != (const instantiation_tuple& src, const instantiation_tuple& dst) {
	return !(src == dst);
}

inline bool operator < (const instantiation_tuple& src, const instantiation_tuple& dst) {
	if (src.length < dst.length) return true;
	else if (src.length > dst.length) return false;
	else if (src.equal_indices.length < dst.equal_indices.length) return true;
	else if (src.equal_indices.length > dst.equal_indices.length) return false;
	else if (src.not_equal_indices.length < dst.not_equal_indices.length) return true;
	else if (src.not_equal_indices.length > dst.not_equal_indices.length) return false;
	else if (src.ge_indices.length < dst.ge_indices.length) return true;
	else if (src.ge_indices.length > dst.ge_indices.length) return false;

	for (unsigned int i = 0; i < src.length; i++) {
		if (src.values[i] < dst.values[i]) return true;
		else if (dst.values[i] < src.values[i]) return false;
	} for (unsigned int i = 0; i < src.equal_indices.length; i++) {
		if (src.equal_indices[i].key < dst.equal_indices[i].key) return true;
		else if (dst.equal_indices[i].key < src.equal_indices[i].key) return true;
		else if (src.equal_indices[i].value < dst.equal_indices[i].value) return true;
		else if (dst.equal_indices[i].value < src.equal_indices[i].value) return true;
	} for (unsigned int i = 0; i < src.not_equal_indices.length; i++) {
		if (src.not_equal_indices[i].key < dst.not_equal_indices[i].key) return true;
		else if (dst.not_equal_indices[i].key < src.not_equal_indices[i].key) return true;
		else if (src.not_equal_indices[i].value < dst.not_equal_indices[i].value) return true;
		else if (dst.not_equal_indices[i].value < src.not_equal_indices[i].value) return true;
	} for (unsigned int i = 0; i < src.ge_indices.length; i++) {
		if (src.ge_indices[i].key < dst.ge_indices[i].key) return true;
		else if (dst.ge_indices[i].key < src.ge_indices[i].key) return true;
		else if (src.ge_indices[i].value < dst.ge_indices[i].value) return true;
		else if (dst.ge_indices[i].value < src.ge_indices[i].value) return true;
	}
	return false;
}

template<typename T>
unsigned int unique_and_cleanup(T* array, size_t length)
{
	unsigned int result = 0;
	for (unsigned int i = 1; i < length; i++) {
		if (array[result] != array[i])
			move(array[i], array[++result]);
		else free(array[i]);
	}
	return result + 1;
}

template<typename T>
void unique_and_cleanup(array<T>& a) {
	a.length = unique_and_cleanup(a.data, a.length);
}

struct relation {
	unsigned int predicate;
	unsigned int arg1; /* `0` here indicates the source vertex */
	unsigned int arg2; /* `0` here indicates the source vertex */

	relation() { }

	relation(unsigned int predicate, unsigned int arg1, unsigned int arg2) :
		predicate(predicate), arg1(arg1), arg2(arg2)
	{ }

	static inline bool is_empty(const relation& key) {
		return key.predicate == 0;
	}

	static inline unsigned int hash(const relation& key) {
		return default_hash(key.predicate) ^ default_hash(key.arg1) ^ default_hash(key.arg2);
	}

	static inline void move(const relation& src, relation& dst) {
		dst.predicate = src.predicate;
		dst.arg1 = src.arg1;
		dst.arg2 = src.arg2;
	}
};

inline bool operator == (const relation& first, const relation& second) {
	return first.predicate == second.predicate
		&& first.arg1 == second.arg1
		&& first.arg2 == second.arg2;
}

inline bool operator != (const relation& first, const relation& second) {
	return first.predicate != second.predicate
		|| first.arg1 != second.arg1
		|| first.arg2 != second.arg2;
}

template<typename ProofCalculus>
struct concept
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Term Term;

	array_map<Term, Proof*> types;
	array_map<Term, Proof*> negated_types;
	array_map<relation, Proof*> relations;
	array_map<relation, Proof*> negated_relations;
	array<Proof*> definitions;
	array<Proof*> existential_intro_nodes;

	array_map<unsigned int, Proof*> function_values;

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, const char* prefix, Printer&&... printer) const {
		for (auto entry : types) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_types) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : relations) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : negated_relations) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (Proof* definition : definitions) {
			if (!print(prefix, out) || !print(*definition->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		} for (auto entry : function_values) {
			if (!print(prefix, out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
		}
		return true;
	}

	bool check_axioms() const {
		bool success = true;
		for (auto entry : types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : negated_types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : relations) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : negated_relations) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		}
		return success;
	}

	static inline void move(const concept<ProofCalculus>& src, concept<ProofCalculus>& dst) {
		core::move(src.types, dst.types);
		core::move(src.negated_types, dst.negated_types);
		core::move(src.relations, dst.relations);
		core::move(src.negated_relations, dst.negated_relations);
		core::move(src.definitions, dst.definitions);
		core::move(src.existential_intro_nodes, dst.existential_intro_nodes);
		core::move(src.function_values, dst.function_values);
	}

	static inline void free(concept<ProofCalculus>& c) {
		/* we set the initial reference_count to 2 since `theory.free_proof`
		   will free definitions when their reference_count is 1 */
		core::free(*c.definitions[0]);
		for (auto entry : c.types) {
			core::free(entry.key);
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.negated_types) {
			core::free(entry.key);
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.relations) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : c.negated_relations) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (Proof* definition : c.definitions) {
			core::free(*definition);
			if (definition->reference_count == 0)
				core::free(definition);
		} for (auto entry : c.function_values) {
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}
		core::free(c.types);
		core::free(c.negated_types);
		core::free(c.relations);
		core::free(c.negated_relations);
		core::free(c.definitions);
		core::free(c.existential_intro_nodes);
		core::free(c.function_values);
	}
};

template<typename ProofCalculus>
inline bool init(concept<ProofCalculus>& c, unsigned int constant) {
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;

	if (!array_map_init(c.types, 8)) {
		return false;
	} else if (!array_map_init(c.negated_types, 8)) {
		free(c.types); return false;
	} else if (!array_map_init(c.relations, 8)) {
		free(c.negated_types);
		free(c.types); return false;
	} else if (!array_map_init(c.negated_relations, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); return false;
	} else if (!array_init(c.definitions, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		return false;
	} else if (!array_init(c.existential_intro_nodes, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); return false;
	} else if (!array_map_init(c.function_values, 8)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.existential_intro_nodes);
		return false;
	}

	Formula* constant_expr = Formula::new_constant(constant);
	if (constant_expr == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(c.existential_intro_nodes); return false;
	}
	Formula* axiom = Formula::new_equals(constant_expr,
			Formula::new_lambda(1, Formula::new_apply(constant_expr, &Term::template variables<1>::value)));
	if (constant_expr == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(*constant_expr); free(constant_expr);
		free(c.existential_intro_nodes); return false;
	}
	Term::template variables<1>::value.reference_count++;
	constant_expr->reference_count += 2 - 1;
	c.definitions[0] = ProofCalculus::new_axiom(axiom);
	if (c.definitions[0] == nullptr) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(*axiom); free(axiom);
		free(c.existential_intro_nodes); return false;
	}
	free(*axiom);
	/* we set the initial reference_count to 2 since `theory.free_proof`
	   will free definitions when their reference_count is 1 */
	c.definitions[0]->reference_count += 2;
	c.definitions.length++;
	return true;
}

template<typename Formula>
inline Formula* preprocess_formula(Formula* src) {
	Formula* first = remove_inverse(src);
	if (first == nullptr) return nullptr;

	Formula* second = same_to_equals(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = remove_exists(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = arg_of_to_arg(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = remove_tense(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = remove_object(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = normalize_set_operations(second);
	free(*second); if (second->reference_count == 0) free(second);
	return first;
}

/* this is useful in `theory.add_definition` where if any sets are created in
   `set_reasoning.get_subset_axiom`, we need the sets to have a specific size */
struct required_set_size {
	unsigned int set_size;
	required_set_size(unsigned int set_size) : set_size(set_size) { }
};

template<typename... Args>
inline bool compute_new_set_size(
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		required_set_size& required,
		Args&&... visitor)
{
	if (required.set_size == UINT_MAX) {
		if (!compute_new_set_size(out, min_set_size, max_set_size, std::forward<Args>(visitor)...))
			return false;
		required.set_size = out;
		return true;
	} else {
		return required.set_size >= min_set_size && required.set_size <= max_set_size
			&& compute_new_set_size(out, required.set_size, required.set_size, std::forward<Args>(visitor)...);
	}
}

template<typename BuiltInConstants, typename ProofCalculus, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus>& sets,
		const required_set_size& required,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline bool compute_new_set_size(unsigned int& out,
		unsigned int min_set_size, unsigned int max_set_size,
		const array<Proof*>& freeable_axioms,
		Args&&... visitor)
{
	return compute_new_set_size(out, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus>& sets,
		const array<typename ProofCalculus::Proof*>& freeable_axioms,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
}

template<typename ProofCalculus, typename Canonicalizer>
struct theory
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	template<unsigned int Value> using Constants = typename Term::template constants<Value>;
	template<unsigned int Value> using Variables = typename Term::template variables<Value>;

	unsigned int new_constant_offset;

	/* A map from `t` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]t` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]t`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]t` or `~[u/0]t`
	   are in the theory. */
	hash_map<Term, pair<array<unsigned int>, array<unsigned int>>> atoms;

	/* A map from `R` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]R` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]R`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]R` or `~[u/0]R`
	   are in the theory. */
	hash_map<relation, pair<array<unsigned int>, array<unsigned int>>> relations;

	concept<ProofCalculus>* ground_concepts;
	unsigned int ground_concept_capacity;
	unsigned int ground_axiom_count;

	hash_set<Proof*> observations;
	set_reasoning<built_in_predicates, ProofCalculus> sets;

	array<pair<Formula*, Proof*>> disjunction_intro_nodes;
	array<pair<Formula*, Proof*>> negated_conjunction_nodes;
	array<pair<Formula*, Proof*>> implication_intro_nodes;
	array<pair<Formula*, Proof*>> existential_intro_nodes;

	Proof* empty_set_axiom;
	Proof* maximality_axiom;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), atoms(64), relations(64),
			ground_concept_capacity(64), ground_axiom_count(0), observations(64),
			disjunction_intro_nodes(16), negated_conjunction_nodes(16),
			implication_intro_nodes(16), existential_intro_nodes(16)
	{
		ground_concepts = (concept<ProofCalculus>*) malloc(sizeof(concept<ProofCalculus>) * ground_concept_capacity);
		if (ground_concepts == NULL) {
			fprintf(stderr, "theory ERROR: Insufficient memory for `ground_concepts`.\n");
			exit(EXIT_FAILURE);
		}
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			ground_concepts[i].types.keys = NULL; /* this is used to indicate that this concept is uninitialized */

		Formula* empty_set_formula = Formula::new_for_all(1, Formula::new_equals(
			Formula::new_equals(Formula::new_atom((unsigned int) built_in_predicates::SIZE, &Variables<1>::value), Formula::new_int(0)),
			Formula::new_not(Formula::new_exists(2, Formula::new_apply(&Variables<1>::value, &Variables<2>::value)))
		));
		if (empty_set_formula == NULL) {
			core::free(ground_concepts);
			exit(EXIT_FAILURE);
		}
		Variables<1>::value.reference_count += 2;
		Variables<2>::value.reference_count++;
		empty_set_axiom = ProofCalculus::new_axiom(empty_set_formula);
		core::free(*empty_set_formula);
		if (empty_set_formula->reference_count == 0)
			core::free(empty_set_formula);
		if (empty_set_axiom == NULL) exit(EXIT_FAILURE);
		empty_set_axiom->reference_count++;

		/* add definition of maximality */
		Formula* left = Formula::new_exists(4, Formula::new_and(
				Formula::new_apply(Formula::new_apply(Term::new_constant((unsigned int) built_in_predicates::GREATEST), &Variables<2>::value), &Variables<4>::value),
				Formula::new_equals(Formula::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG1), &Variables<4>::value), &Variables<3>::value),
				Formula::new_equals(Formula::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG2), &Variables<4>::value), &Variables<1>::value)
			));
		if (left == nullptr) exit(EXIT_FAILURE);
		Variables<1>::value.reference_count++;
		Variables<2>::value.reference_count++;
		Variables<3>::value.reference_count++;
		Variables<4>::value.reference_count += 3;
		Formula* right = Formula::new_for_all(4, Formula::new_for_all(5, Formula::new_if_then(
				Formula::new_and(
					Formula::new_apply(&Variables<2>::value, &Variables<4>::value),
					Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG1), &Variables<4>::value), &Variables<1>::value),
					Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG2), &Variables<4>::value), &Variables<5>::value)
				),
				Formula::new_for_all(6, Formula::new_if_then(
					Formula::new_apply(&Variables<3>::value, &Variables<6>::value),
					Formula::new_for_all(7, Formula::new_for_all(8, Formula::new_if_then(
						Formula::new_and(
							Formula::new_apply(&Variables<2>::value, &Variables<7>::value),
							Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG1), &Variables<7>::value), &Variables<6>::value),
							Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG2), &Variables<7>::value), &Variables<8>::value)
						),
						Formula::new_apply(Term::new_constant((unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL), &Variables<5>::value, &Variables<8>::value)
					)))
				))
			)));
		if (right == nullptr) exit(EXIT_FAILURE);
		Variables<1>::value.reference_count++;
		Variables<2>::value.reference_count += 2;
		Variables<3>::value.reference_count++;
		Variables<4>::value.reference_count += 3;
		Variables<5>::value.reference_count += 2;
		Variables<6>::value.reference_count += 2;
		Variables<7>::value.reference_count += 3;
		Variables<8>::value.reference_count += 2;
		maximality_axiom = get_subset_axiom<true>(left, right, 3);
		free(*left); if (left->reference_count == 0) free(left);
		free(*right); if (right->reference_count == 0) free(right);
		if (maximality_axiom == NULL) exit(EXIT_FAILURE);
		maximality_axiom->reference_count++;
	}

	~theory() {
		for (auto entry : atoms) {
			core::free(entry.key);
			core::free(entry.value.key);
			core::free(entry.value.value);
		} for (auto entry : relations) {
			core::free(entry.value.key);
			core::free(entry.value.value);
		} for (Proof* proof : observations) {
			core::free(*proof);
			if (proof->reference_count == 0)
				core::free(proof);
		} for (auto entry : disjunction_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : negated_conjunction_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : implication_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		} for (auto entry : existential_intro_nodes) {
			core::free(*entry.key);
			if (entry.key->reference_count == 0)
				core::free(entry.key);
		}

		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;

			auto& c = ground_concepts[i];
			for (unsigned int i = 0; i < c.definitions.length; i++) {
				Proof* definition = c.definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					for (unsigned int j = i + 1; j < c.definitions.length; j++) {
						Proof* other_definition = c.definitions[j];
						if (other_definition->formula->binary.right->type != FormulaType::LAMBDA) continue;

						Proof* axiom = get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, 1);
						free(*axiom);
						if (axiom->reference_count == 1)
							sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand);

						axiom = get_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1);
						free(*axiom);
						if (axiom->reference_count == 1)
							sets.free_subset_axiom(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand);
					}
				}
			}

			core::free(c);
		}
		core::free(ground_concepts);

		core::free(*empty_set_axiom);
		if (empty_set_axiom->reference_count == 0)
			core::free(empty_set_axiom);
		core::free(*maximality_axiom);
		sets.free_subset_axiom(maximality_axiom);
	}

	unsigned int get_free_concept_id() {
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			if (ground_concepts[i].types.keys == NULL) return i + new_constant_offset;

		unsigned int new_capacity = ground_concept_capacity;
		expand_capacity(new_capacity, ground_concept_capacity + 1);

		if (!resize(ground_concepts, new_capacity))
			return EXIT_FAILURE;
		for (unsigned int i = ground_concept_capacity; i < new_capacity; i++)
			ground_concepts[i].types.keys = NULL; /* this is used to indicate that this concept is uninitialized */
		unsigned int old_capacity = ground_concept_capacity;
		ground_concept_capacity = new_capacity;
		return new_constant_offset + old_capacity;
	}

	inline bool try_init_concept(unsigned int id) {
		if (id < new_constant_offset) return true;
		concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.keys != NULL)
			return true;
		else return init(c, id);
	}

	void free_concept_id(unsigned int id) {
#if !defined(NDEBUG)
		if (id < new_constant_offset)
			fprintf(stderr, "theory.free_concept_id WARNING: The given `id` is less than `new_constant_offset`.\n");
#endif
		core::free(ground_concepts[id - new_constant_offset]);
		ground_concepts[id - new_constant_offset].types.keys = NULL;
	}

	void try_free_concept_id(unsigned int id) {
		if (id < new_constant_offset) return;
		const concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.size != 0 || c.negated_types.size != 0 || c.relations.size != 0
		 || c.negated_relations.size != 0 || c.definitions.length > 1 || c.definitions[0]->reference_count > 2
		 || c.existential_intro_nodes.length != 0 || c.function_values.size != 0)
			return;

		bool contains;
		unsigned int count = sets.symbols_in_formulas.counts.get(id, contains);
		if (contains && count > 0) return;

		if (sets.element_map.table.contains({&id, 1})) return;

		free_concept_id(id);
	}

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;
			if (!print("Concept ", out) || !print(i + new_constant_offset, out, std::forward<Printer>(printer)...) || !print(":\n", out)
			 || !ground_concepts[i].print_axioms(out, "  ", std::forward<Printer>(printer)...)) return false;
		}
		return sets.print_axioms(out, std::forward<Printer>(printer)...);
	}

	Proof* add_formula(Formula* formula, unsigned int& new_constant)
	{
		new_constant = 0;

		Formula* new_formula = preprocess_formula(formula);
		if (new_formula == NULL) return nullptr;

		array_map<unsigned int, unsigned int> variable_map(16);
		Formula* canonicalized = Canonicalizer::canonicalize(*new_formula, variable_map);
		free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
		if (canonicalized == NULL) return nullptr;

/* TODO: for debugging; delete this */
print_axioms(stderr);
print("canonicalized: ", stderr); print(*canonicalized, stderr); print('\n', stderr);
		Proof* new_proof = make_proof<false, true, true>(canonicalized, new_constant);
print_axioms(stderr);
if (new_proof != NULL) {
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put((unsigned int) built_in_predicates::UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof<built_in_predicates, typename ProofCalculus::ProofCanonicalizer>(*new_proof, expected_conclusion))
fprintf(stderr, "add_formula WARNING: `check_proof` failed.\n");
free(*expected_conclusion); if (expected_conclusion->reference_count == 0) free(expected_conclusion);
}
		core::free(*canonicalized);
		if (canonicalized->reference_count == 0)
			core::free(canonicalized);
		if (new_proof == NULL) {
			return nullptr;
		} else if (!observations.add(new_proof)) {
			free_proof(new_proof);
			return nullptr;
		}
		return new_proof;
	}

	void remove_formula(Proof* proof) {
		observations.remove(proof);
		free_proof(proof);
	}

	template<ProofType Type>
	bool check_disjunction_introductions(const array<pair<Formula*, Proof*>>& proof_steps) const
	{
		bool success = true;
		array<const Proof*> computed_steps(64);
		for (const Proof* proof : observations) {
			array<const Proof*> steps(16);
			if (!get_proof_steps<Type>(proof, steps)) {
				fprintf(stderr, "theory.check_proof_disjunctions ERROR: `get_proof_steps` failed.\n");
				return false;
			}

			for (const Proof* step : steps) {
				bool contains = false;
				for (const auto& entry : proof_steps) {
					if (*step == *entry.value) {
						contains = true;
						break;
					}
				}
				if (!contains) {
					print("theory.check_proof_disjunctions WARNING: Found step "
							"of a proof in `observations` that is not in `proof_steps`.\n", stderr);
					success = false;
				}
			}

			for (const Proof* step : steps) {
				if (!computed_steps.contains(step) && !computed_steps.add(step))
					return false;
			}
		}

		for (const auto& entry : proof_steps) {
			if (!computed_steps.contains(entry.value)) {
				print("theory.check_proof_disjunctions WARNING: `proof_steps` "
						"contains a step that does not belong to a proof in `observations`.\n", stderr);
				print("  Formula: ", stderr); print(*entry.key, stderr); print('\n', stderr);
				success = false;
			}
		}

		return success;
	}

	inline bool check_disjunction_introductions() const {
		return check_disjunction_introductions<ProofType::DISJUNCTION_INTRODUCTION>(disjunction_intro_nodes)
			&& check_disjunction_introductions<ProofType::EXISTENTIAL_INTRODUCTION>(existential_intro_nodes)
			&& check_disjunction_introductions<ProofType::IMPLICATION_INTRODUCTION>(implication_intro_nodes);
	}

	inline bool check_concept_axioms() const {
		bool success = true;
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;
			success &= ground_concepts[i].check_axioms();
		}
		return success;
	}

	inline bool are_elements_provable() const {
		bool success = true;
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == nullptr)
				continue;
			Formula* set_formula = sets.sets[i].size_axiom->formula->binary.left->binary.right;

			array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				quantifiers[quantifiers.length++] = set_formula;
				set_formula = set_formula->quantifier.operand;
			}

			for (unsigned int j = 0; j < sets.sets[i].element_count(); j++) {
				unsigned int* element = (unsigned int*) alloca(sizeof(unsigned int) * sets.sets[i].arity);
				const unsigned int* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * j);
				for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
					element[k] = element_src[k];
				array<instantiation_tuple> possible_values(1);
				if (!init(possible_values[0], sets.sets[i].arity))
					return false;
				possible_values.length = 1;
				for (uint_fast8_t k = 0; k < possible_values[0].length; k++) {
					free(possible_values[0].values[k]);
					possible_values[0].values[k].type = instantiation_type::CONSTANT;
					possible_values[0].values[k].constant = element[k];
				}

				sets.sets[i].remove_element_at(j);
				if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values)) {
					print("theory.are_elements_provable ERROR: The element ", stderr);
					if (sets.sets[i].arity == 1)
						print(element[0], stderr);
					else print(element, sets.sets[i].arity, stderr);
					print(" does not provably belong to set with ID ", stderr);
					print(i, stderr); print(".\n", stderr);
					success = false;
				}
				for (auto& element : possible_values) free(element);
				for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
					sets.sets[i].elements[sets.sets[i].elements.length++] = element[k];
			}
		}
		return success;
	}

	inline bool is_provably_not_a_set(unsigned int constant) const {
		if (constant < new_constant_offset || ground_concepts[constant - new_constant_offset].types.keys == NULL)
			return false;
		return ground_concepts[constant - new_constant_offset].types.size != 0
			|| ground_concepts[constant - new_constant_offset].function_values.size != 0;
	}

	inline bool is_provably_a_set(unsigned int constant) const {
		if (constant < new_constant_offset || ground_concepts[constant - new_constant_offset].types.keys == NULL)
			return false;
		const Formula* first_set_definition = ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.right->quantifier.operand;
		if (sets.set_ids.table.contains(*first_set_definition))
			return true;

		Term* atom = Term::new_apply(Term::new_constant(constant), &Variables<1>::value);
		if (atom == nullptr)
			return false;
		Variables<1>::value.reference_count++;

		bool contains;
		pair<array<unsigned int>, array<unsigned int>>& type_instances = atoms.get(*atom, contains);
		free(*atom); free(atom);
		if (contains && type_instances.key.length != 0)
			return true;
		return false;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_unary_atom(Term& atom, Proof* axiom, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (atom.binary.right->type != TermType::CONSTANT)
			fprintf(stderr, "add_unary_atom WARNING: The operand of this application is not constant.\n");
#endif
		unsigned int arg = atom.binary.right->constant;

		Formula* lifted_literal;
		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom);
				return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}

#if !defined(NDEBUG)
		if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_unary_atom WARNING: `ground_concepts` does not contain the concept %u.\n", arg);
#endif
		bool contains; unsigned int bucket;
		if (!atoms.check_size()) {
			free(*lifted_literal); free(lifted_literal);
			return false;
		}

		pair<array<unsigned int>, array<unsigned int>>& instance_pair = atoms.get(*lifted_literal, contains, bucket);
		if (!contains) {
			if (!array_init(instance_pair.key, 8)) {
				free(*lifted_literal); free(lifted_literal);
				return false;
			} else if (!array_init(instance_pair.value, 8)) {
				core::free(instance_pair.key);
				free(*lifted_literal); free(lifted_literal);
				return false;
			}
			atoms.table.keys[bucket] = *lifted_literal;
			atoms.table.size++;
		}

		array<unsigned int>& instances = (Negated ? instance_pair.value : instance_pair.key);
		array_map<Term, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);
		if (!instances.ensure_capacity(instances.length + 1)
		 || !ground_types.ensure_capacity(ground_types.size + 1))
		{
			free(*lifted_literal); free(lifted_literal);
			return false;
		}

		add_sorted<false>(instances, arg);
		ground_types.keys[ground_types.size] = *lifted_literal;
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		ground_axiom_count++;
		free(*lifted_literal); free(lifted_literal);

		if (!check_set_membership_after_addition<ResolveInconsistencies>(&atom, std::forward<Args>(visitor)...)) {
			remove_unary_atom<Negated>(atom, std::forward<Args>(visitor)...);
			return false;
		}
		return true;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_binary_atom(relation rel, Proof* axiom, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (rel.arg1 < new_constant_offset || ground_concepts[rel.arg1 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (rel.arg2 < new_constant_offset || ground_concepts[rel.arg2 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.add_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg2);
#endif

		bool contains; unsigned int bucket;
		if (!relations.check_size(relations.table.size + 4)) return false;
		pair<array<unsigned int>, array<unsigned int>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(predicate_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(predicate_instance_pair.value, 8)) {
				core::free(predicate_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {0, rel.arg1, rel.arg2};
			relations.table.size++;
		}

		pair<array<unsigned int>, array<unsigned int>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(arg1_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(arg1_instance_pair.value, 8)) {
				core::free(arg1_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, 0, rel.arg2};
			relations.table.size++;
		}

		pair<array<unsigned int>, array<unsigned int>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0}, contains, bucket);
		if (!contains) {
			if (!array_init(arg2_instance_pair.key, 8)) {
				return false;
			} else if (!array_init(arg2_instance_pair.value, 8)) {
				core::free(arg2_instance_pair.key); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, rel.arg1, 0};
			relations.table.size++;
		}

		array<unsigned int>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<unsigned int>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<unsigned int>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);
		if (!predicate_instances.ensure_capacity(predicate_instances.length + 1)
		 || !arg1_instances.ensure_capacity(arg1_instances.length + 1)
		 || !arg2_instances.ensure_capacity(arg2_instances.length + 1)
		 || !ground_arg1.ensure_capacity(ground_arg1.size + 3)
		 || !ground_arg2.ensure_capacity(ground_arg2.size + 3)) return false;

		if (rel.arg1 == rel.arg2) {
			pair<array<unsigned int>, array<unsigned int>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0}, contains, bucket);
			if (!contains) {
				if (!array_init(both_arg_instance_pair.key, 8)) {
					return false;
				} else if (!array_init(both_arg_instance_pair.value, 8)) {
					core::free(both_arg_instance_pair.key);
					return false;
				}
				relations.table.keys[bucket] = {rel.predicate, 0, 0};
				relations.table.size++;
			}

			array<unsigned int>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			if (!both_arg_instances.ensure_capacity(both_arg_instances.length + 1))
				return false;
			both_arg_instances[both_arg_instances.length++] = rel.arg1;
			insertion_sort(both_arg_instances);
		}

		add_sorted<false>(predicate_instances, rel.predicate);
		add_sorted<false>(arg1_instances, rel.arg1);
		add_sorted<false>(arg2_instances, rel.arg2);
		ground_arg1.keys[ground_arg1.size] = {rel.predicate, 0, rel.arg2};
		ground_arg1.values[ground_arg1.size++] = axiom;
		ground_arg2.keys[ground_arg2.size] = {rel.predicate, rel.arg1, 0};
		ground_arg2.values[ground_arg2.size++] = axiom;
		axiom->reference_count += 2;
		ground_axiom_count += 2;
		if (rel.arg1 == rel.arg2) {
			/* in this case, `ground_arg1` and `ground_arg2` are the same */
			ground_arg1.keys[ground_arg1.size] = {rel.predicate, 0, 0};
			ground_arg1.values[ground_arg1.size++] = axiom;
			axiom->reference_count++;
			ground_axiom_count++;
		}

		Formula* atom = Formula::new_atom(rel.predicate, Term::new_constant(rel.arg1), Term::new_constant(rel.arg2));
		if (atom == nullptr) return false;
		if (!check_set_membership_after_addition<ResolveInconsistencies>(atom, std::forward<Args>(visitor)...)) {
			remove_binary_atom<Negated>(rel, std::forward<Args>(visitor)...);
			free(*atom); free(atom);
			return false;
		}
		free(*atom); free(atom);
		return true;
	}

	template<bool Negated, typename... Args>
	inline bool remove_unary_atom(Term& atom, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (atom.binary.right->type != TermType::CONSTANT)
			fprintf(stderr, "remove_unary_atom WARNING: The operand of this application is not constant.\n");
#endif
		unsigned int arg = atom.binary.right->constant;

		Formula* lifted_literal;
		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				free(*lifted_atom); free(lifted_atom);
				return false;
			}
		} else {
			lifted_literal = lifted_atom;
		}

#if !defined(NDEBUG)
		if (!atoms.table.contains(*lifted_literal)) {
			print("theory.remove_unary_atom WARNING: `atoms` does not contain the key ", stderr);
			print(*lifted_literal, stderr); print(".\n", stderr);
		} if (arg < new_constant_offset || ground_concepts[arg - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", arg);
#endif

		array<unsigned int>& instances = (Negated ? atoms.get(*lifted_literal).value : atoms.get(*lifted_literal).key);
		array_map<Term, Proof*>& ground_types = (Negated ? ground_concepts[arg - new_constant_offset].negated_types : ground_concepts[arg - new_constant_offset].types);

		unsigned int index = instances.index_of(arg);
#if !defined(NDEBUG)
		if (index == instances.length)
			fprintf(stderr, "theory.remove_unary_atom WARNING: `instances` does not contain %u.\n", arg);
#endif
		shift_left(instances.data + index, instances.length - index - 1);
		instances.length--;

		index = ground_types.index_of(*lifted_literal);
#if !defined(NDEBUG)
		if (index == ground_types.size) {
			print("theory.remove_unary_atom WARNING: `ground_types` does not contain ", stderr);
			print(*lifted_literal, stderr); print(".\n", stderr);
		}
#endif
		Proof* axiom = ground_types.values[index];
		free(*axiom); if (axiom->reference_count == 0) free(axiom);
		free(ground_types.keys[index]);
		ground_types.remove_at(index);
		ground_axiom_count--;
		free(*lifted_literal); free(lifted_literal);
		return check_set_membership_after_subtraction(&atom, std::forward<Args>(visitor)...);
	}

	template<bool Negated, typename... Args>
	inline bool remove_binary_atom(relation rel, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (!relations.table.contains({0, rel.arg1, rel.arg2})
		 || !relations.table.contains({rel.predicate, 0, rel.arg2})
		 || !relations.table.contains({rel.predicate, rel.arg1, 0}))
			fprintf(stderr, "theory.remove_binary_atom WARNING: `relations` does not contain the necessary relations.\n");
		if (rel.arg1 < new_constant_offset || ground_concepts[rel.arg1 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg1);
		if (rel.arg2 < new_constant_offset || ground_concepts[rel.arg2 - new_constant_offset].types.keys == NULL)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `ground_concepts` does not contain the key %u.\n", rel.arg2);
#endif

		pair<array<unsigned int>, array<unsigned int>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2});
		pair<array<unsigned int>, array<unsigned int>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0});

		array<unsigned int>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<unsigned int>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<unsigned int>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);

		unsigned int index;
		if (rel.arg1 == rel.arg2) {
			pair<array<unsigned int>, array<unsigned int>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0});
			array<unsigned int>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			index = both_arg_instances.index_of(rel.arg1);
#if !defined(NDEBUG)
			if (index == both_arg_instances.length)
				fprintf(stderr, "theory.remove_binary_atom WARNING: `both_arg_instances` does not contain %u.\n", rel.arg1);
#endif
			shift_left(both_arg_instances.data + index, both_arg_instances.length - index - 1);
			both_arg_instances.length--;
		}

		index = predicate_instances.index_of(rel.predicate);
#if !defined(NDEBUG)
		if (index == predicate_instances.length)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `predicate_instances` does not contain %u.\n", rel.predicate);
#endif
		shift_left(predicate_instances.data + index, predicate_instances.length - index - 1);
		predicate_instances.length--;

		index = arg1_instances.index_of(rel.arg1);
#if !defined(NDEBUG)
		if (index == arg1_instances.length)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `arg1_instances` does not contain %u.\n", rel.arg1);
#endif
		shift_left(arg1_instances.data + index, arg1_instances.length - index - 1);
		arg1_instances.length--;

		index = arg2_instances.index_of(rel.arg2);
#if !defined(NDEBUG)
		if (index == arg2_instances.length)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `arg2_instances` does not contain %u.\n", rel.arg2);
#endif
		shift_left(arg2_instances.data + index, arg2_instances.length - index - 1);
		arg2_instances.length--;

		index = ground_arg1.index_of(relation(rel.predicate, 0, rel.arg2));
#if !defined(NDEBUG)
		if (index == ground_arg1.size)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `ground_arg1` does not contain the requested predicate.\n");
#endif
		Proof* axiom = ground_arg1.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		ground_arg1.remove_at(index);

		index = ground_arg2.index_of(relation(rel.predicate, rel.arg1, 0));
#if !defined(NDEBUG)
		if (index == ground_arg2.size)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `ground_arg2` does not contain the requested predicate.\n");
#endif
		axiom = ground_arg2.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		ground_arg2.remove_at(index);
		ground_axiom_count -= 2;

		if (rel.arg1 == rel.arg2) {
			/* in this case, `ground_arg1` and `ground_arg2` are the same */
			index = ground_arg1.index_of(relation(rel.predicate, 0, 0));
#if !defined(NDEBUG)
			if (index == ground_arg1.size)
				fprintf(stderr, "theory.remove_binary_atom WARNING: `ground_arg1` does not contain the requested predicate.\n");
#endif
			axiom = ground_arg1.values[index];
			core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
			ground_arg1.remove_at(index);
			ground_axiom_count--;
		}

		Formula* atom = Formula::new_atom(rel.predicate, Term::new_constant(rel.arg1), Term::new_constant(rel.arg2));
		if (atom == nullptr) return false;
		bool result = check_set_membership_after_subtraction(atom, std::forward<Args>(visitor)...);
		free(*atom); free(atom);
		return result;
	}

	unsigned int index_of(const array<pair<Formula*, Proof*>>& elements, Proof* proof) const {
		unsigned int index = elements.length;
		for (unsigned int i = 0; i < elements.length; i++)
			if (elements[i].value == proof) { index = i; break; }
#if !defined(NDEBUG)
		if (index == elements.length)
			fprintf(stderr, "theory.index_of WARNING: `elements` does not contain `proof`.\n");
#endif
		return index;
	}

	enum class change_type {
		UNARY_ATOM,
		NEGATED_UNARY_ATOM,
		BINARY_ATOM,
		NEGATED_BINARY_ATOM,
		SUBSET_AXIOM,
		SET_SIZE_AXIOM,
		DEFINITION,
		FUNCTION_VALUE,
		IMPLICATION_INTRO_NODE,
		NEGATED_CONJUNCTION_NODE,
		EXISTENTIAL_INTRO_NODE,
		DISJUNCTION_INTRO_NODE
	};

	struct change {
		change_type type;
		union {
			pair<Term, Proof*> unary_atom;
			pair<relation, Proof*> binary_atom;
			Proof* axiom;
			pair<Formula*, Proof*> intro_node;
		};

		static inline void free(change& c) {
			switch (c.type) {
			case change_type::UNARY_ATOM:
			case change_type::NEGATED_UNARY_ATOM:
				core::free(c.unary_atom.key);
				core::free(*c.unary_atom.value); if (c.unary_atom.value->reference_count == 0) core::free(c.unary_atom.value);
				return;
			case change_type::BINARY_ATOM:
			case change_type::NEGATED_BINARY_ATOM:
				core::free(*c.binary_atom.value); if (c.binary_atom.value->reference_count == 0) core::free(c.binary_atom.value);
				return;
			case change_type::SUBSET_AXIOM:
			case change_type::SET_SIZE_AXIOM:
			case change_type::DEFINITION:
			case change_type::FUNCTION_VALUE:
				core::free(*c.axiom); if (c.axiom->reference_count == 0) core::free(c.axiom);
				return;
			case change_type::IMPLICATION_INTRO_NODE:
			case change_type::NEGATED_CONJUNCTION_NODE:
			case change_type::EXISTENTIAL_INTRO_NODE:
			case change_type::DISJUNCTION_INTRO_NODE:
				core::free(*c.intro_node.key); if (c.intro_node.key->reference_count == 0) core::free(c.intro_node.key);
				core::free(*c.intro_node.value); if (c.intro_node.value->reference_count == 0) core::free(c.intro_node.value);
				return;
			}
			fprintf(stderr, "theory.change.free ERROR: Unrecognized `change_type`.\n");
			exit(EXIT_FAILURE);
		}
	};

	struct changes {
		array<change> list;

		changes() : list(8) { }

		~changes() {
			for (change& c : list)
				core::free(c);
		}

		inline bool add(change_type type, const pair<Term, Proof*>& unary_atom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length].unary_atom.key = unary_atom.key;
			list[list.length++].unary_atom.value = unary_atom.value;
			unary_atom.value->reference_count++;
			return true;
		}

		inline bool add(change_type type, const pair<relation, Proof*>& binary_atom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].binary_atom = binary_atom;
			binary_atom.value->reference_count++;
			return true;
		}

		inline bool add(change_type type, Proof* axiom) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].axiom = axiom;
			axiom->reference_count++;
			return true;
		}

		inline bool add(change_type type, const pair<Formula*, Proof*>& intro_node) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length++].intro_node = intro_node;
			intro_node.key->reference_count++;
			intro_node.value->reference_count++;
			return true;
		}
	};

	template<typename... Args>
	bool add_changes(theory::changes& changes, Args&&... visitor)
	{
		for (change& c : changes.list) {
			array_multiset<unsigned int> constants(8);
			switch (c.type) {
			case change_type::UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<false, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::NEGATED_UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<true, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<false, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<true, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::SUBSET_AXIOM:
				{
					unsigned int variable = c.axiom->formula->quantifier.variable;
					if (c.axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
						unsigned int predicate; Term const* arg1; Term const* arg2;
						Formula* left = c.axiom->formula->quantifier.operand->binary.left;
						if (is_atomic(*left, predicate, arg1, arg2)
						 && arg1->type == TermType::VARIABLE
						 && arg1->variable == variable && arg2 == NULL)
						{
							if (!try_init_concept(predicate)) return false;
						}
					}
					if (!sets.add_subset_axiom(c.axiom, std::forward<Args>(visitor)...)) return false;
				}
				continue;
			case change_type::SET_SIZE_AXIOM:
				{
					unsigned int set_id, arity = 1;
					Formula* set_formula = c.axiom->formula->binary.left->binary.right->quantifier.operand;
					while (set_formula->type == TermType::LAMBDA) {
						set_formula = set_formula->quantifier.operand;
						arity++;
					}
					if (!sets.get_set_id(set_formula, arity, set_id, std::forward<Args>(visitor)...)
					 || !sets.sets[set_id].set_size_axiom(c.axiom)) return false;
				}
				continue;
			case change_type::DEFINITION:
				if (!add_definition<false>(c.axiom, std::forward<Args>(visitor)...))
					return false;
				continue;
			case change_type::FUNCTION_VALUE:
				if (!add_function_value<false>(c.axiom, std::forward<Args>(visitor)...))
					return false;
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				if (!implication_intro_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				if (!negated_conjunction_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					if (!existential_intro_nodes.add(c.intro_node)) return false;

					Term* term;
					if (c.intro_node.value->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.value->operands[2]->term;
					else term = c.intro_node.value->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						if (!try_init_concept(term->constant)
						 || !ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.add(c.intro_node.value)) {
							unsigned int index = existential_intro_nodes.index_of(c.intro_node);
							existential_intro_nodes.remove(index);
							return false;
						}
					}
					c.intro_node.key->reference_count++;
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				if (!disjunction_intro_nodes.add(c.intro_node)) return false;
				c.intro_node.key->reference_count++;
				continue;
			}
			fprintf(stderr, "theory.add_changes ERROR: Unrecognized `change_type`.\n");
			return false;
		}
		return true;
	}

	template<typename... Args>
	bool subtract_changes(theory::changes& changes, Args&&... visitor)
	{
		on_subtract_changes(std::forward<Args>(visitor)...);
		array<Proof*> freeable_set_size_axioms(8);
		for (unsigned int i = changes.list.length; i > 0; i--) {
			change& c = changes.list[i - 1];
			array_multiset<unsigned int> constants(8);
			switch (c.type) {
			case change_type::UNARY_ATOM:
				remove_unary_atom<false>(c.unary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					try_free_concept_id(constants.counts.keys[i]);
				continue;
			case change_type::NEGATED_UNARY_ATOM:
				remove_unary_atom<true>(c.unary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					try_free_concept_id(constants.counts.keys[i]);
				continue;
			case change_type::BINARY_ATOM:
				remove_binary_atom<false>(c.binary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				try_free_concept_id(c.binary_atom.key.predicate);
				try_free_concept_id(c.binary_atom.key.arg1);
				try_free_concept_id(c.binary_atom.key.arg2);
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				remove_binary_atom<true>(c.binary_atom.key, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				try_free_concept_id(c.binary_atom.key.predicate);
				try_free_concept_id(c.binary_atom.key.arg1);
				try_free_concept_id(c.binary_atom.key.arg2);
				continue;
			case change_type::SUBSET_AXIOM:
				{
					unsigned int variable = c.axiom->formula->quantifier.variable;
					if (c.axiom->formula->quantifier.operand->type == FormulaType::IF_THEN) {
						unsigned int predicate; Term const* arg1; Term const* arg2;
						Formula* left = c.axiom->formula->quantifier.operand->binary.left;
						if (is_atomic(*left, predicate, arg1, arg2)
						 && arg1->type == TermType::VARIABLE
						 && arg1->variable == variable && arg2 == NULL)
						{
							try_free_concept_id(predicate);
						}
					}
					unsigned int child_count = c.axiom->children.length;
					if (changes.list.length == 1) child_count++;
					if (c.axiom->reference_count == child_count + 2) {
						sets.free_subset_axiom(c.axiom, freeable_set_size_axioms, std::forward<Args>(visitor)...);
					} else {
						free(*c.axiom);
					}
					continue;
				}
			case change_type::SET_SIZE_AXIOM:
				{
					unsigned int set_id, arity = 1;
					Formula* set_formula = c.axiom->formula->binary.left->binary.right->quantifier.operand;
					while (set_formula->type == TermType::LAMBDA) {
						set_formula = set_formula->quantifier.operand;
						arity++;
					}
					if (!sets.get_set_id(set_formula, arity, set_id)
					 || !freeable_set_size_axioms.add(c.axiom)) return false;
					unsigned int old_ref_count = c.axiom->reference_count;
					c.axiom->reference_count = 1;
					bool is_freeable = sets.is_freeable(set_id);
					c.axiom->reference_count = old_ref_count;
					if (is_freeable) {
						on_free_set(set_id, sets, std::forward<Args>(visitor)...);
						sets.free_set(set_formula, set_id);
					}
					continue;
				}
			case change_type::DEFINITION:
				remove_definition(c.axiom, freeable_set_size_axioms, std::forward<Args>(visitor)...);
				continue;
			case change_type::FUNCTION_VALUE:
				remove_function_value(c.axiom, std::forward<Args>(visitor)...);
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				implication_intro_nodes.remove(index_of(implication_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				negated_conjunction_nodes.remove(index_of(negated_conjunction_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					Term* term;
					if (c.intro_node.value->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.value->operands[2]->term;
					else term = c.intro_node.value->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						unsigned int index = ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.index_of(c.intro_node.value);
						ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.remove(index);
						try_free_concept_id(term->constant);
					}

					existential_intro_nodes.remove(index_of(existential_intro_nodes, c.intro_node.value));
					c.intro_node.key->reference_count--;
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				disjunction_intro_nodes.remove(index_of(disjunction_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				continue;
			}
			fprintf(stderr, "theory.subtract_changes ERROR: Unrecognized `change_type`.\n");
			return false;
		}
		return true;
	}

	template<typename... Visitor>
	bool get_theory_changes(
			Proof& proof, array<const Proof*>& discharged_axioms,
			array_map<const Proof*, unsigned int>& reference_counts,
			theory::changes& changes, Visitor&&... visitor) const
	{
		visit_node(proof, std::forward<Visitor>(visitor)...);
#if !defined(NDEBUG)
		bool contains;
		unsigned int reference_count = reference_counts.get(&proof, contains);
		if (!contains) fprintf(stderr, "theory.remove_proof WARNING: The given proof is not in the map `reference_counts`.\n");
#else
		unsigned int reference_count = reference_counts.get(&proof);
#endif

		unsigned int old_discharged_axiom_count = discharged_axioms.length;

		Formula* formula = nullptr;
		Term* predicate = nullptr;
		Term* arg1 = nullptr;
		Term* arg2 = nullptr;
		switch (proof.type) {
		case ProofType::AXIOM:
			if (discharged_axioms.contains(&proof)) break;
			if (is_atomic(*proof.formula, predicate, arg1, arg2)) {
				if (reference_count != 1) return true;
				if (arg2 == NULL) {
					if (!visit_unary_atom<false>(proof.formula, std::forward<Visitor>(visitor)...)
					 || !changes.add(change_type::UNARY_ATOM, {*proof.formula, &proof})) return false;
				} else {
					if (predicate->type != TermType::CONSTANT) return false;
					if (!visit_binary_atom<false>(predicate->constant, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
					 || !changes.add(change_type::BINARY_ATOM, {{predicate->constant, arg1->constant, arg2->constant}, &proof})) return false;
				}
			} else if (proof.formula->type == FormulaType::NOT) {
				if (is_atomic(*proof.formula->unary.operand, predicate, arg1, arg2)) {
					if (reference_count != 1) return true;
					if (arg2 == NULL) {
						if (!visit_unary_atom<true>(proof.formula->unary.operand, std::forward<Visitor>(visitor)...)
						 || !changes.add(change_type::NEGATED_UNARY_ATOM, {*proof.formula->unary.operand, &proof})) return false;
					} else {
						if (predicate->type != TermType::CONSTANT) return false;
						if (!visit_binary_atom<true>(predicate->constant, arg1->constant, arg2->constant, std::forward<Visitor>(visitor)...)
						 || !changes.add(change_type::NEGATED_BINARY_ATOM, {{predicate->constant, arg1->constant, arg2->constant}, &proof})) return false;
					}
				} else {
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected axiom.\n");
				}
			} else if (proof.formula->type == FormulaType::FOR_ALL
					&& proof.formula->quantifier.operand->type == FormulaType::IF_THEN)
			{
				if (!visit_subset_axiom(proof, std::forward<Visitor>(visitor)...))
					return false;
				if (reference_count != 1) return true;
				return changes.add(change_type::SUBSET_AXIOM, &proof);
			} else if (proof.formula->type == FormulaType::EQUALS) {
				bool atomic = is_atomic(*proof.formula->binary.left, predicate, arg1, arg2);
				if (atomic && predicate->type == TermType::CONSTANT
				 && predicate->constant == (unsigned int) built_in_predicates::SIZE
				 && arg1->type == TermType::LAMBDA && arg2 == NULL
				 && proof.formula->binary.right->type == TermType::INTEGER)
				{
					if (reference_count != 1) return true;
					return changes.add(change_type::SET_SIZE_AXIOM, &proof);
				} else if (proof.formula->binary.left->type == TermType::CONSTANT) {
					if (reference_count != 1) return true;
					/* make sure we keep track of decrementing the reference counts of bidirectional subset edges */
					if (proof.formula->binary.right->type == TermType::LAMBDA) {
						unsigned int concept_id = proof.formula->binary.left->constant;
						for (Proof* other_definition : ground_concepts[concept_id - new_constant_offset].definitions) {
							if (proof.formula->binary.right == other_definition->formula->binary.right
							 || other_definition->formula->binary.right->type != FormulaType::LAMBDA)
								continue;
							if (!reference_counts.ensure_capacity(reference_counts.size + 2))
								return false;

							Proof* axiom = sets.template get_existing_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, proof.formula->binary.right->quantifier.operand);
							unsigned int index = reference_counts.index_of(axiom);
							if (index == reference_counts.size) {
								reference_counts.keys[index] = axiom;
								reference_counts.values[index] = axiom->reference_count;
								reference_counts.size++;
							}
							reference_counts.values[index]--;

							axiom = sets.template get_existing_subset_axiom<false>(proof.formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand);
							index = reference_counts.index_of(axiom);
							if (index == reference_counts.size) {
								reference_counts.keys[index] = axiom;
								reference_counts.values[index] = axiom->reference_count;
								reference_counts.size++;
							}
							reference_counts.values[index]--;
						}
					}
					return changes.add(change_type::DEFINITION, &proof);
				} else if (proof.formula->binary.left->type == TermType::UNARY_APPLICATION) {
					if (reference_count != 1) return true;
					return changes.add(change_type::FUNCTION_VALUE, &proof);
				} else {
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected equals axiom.\n");
				}
			} else {
				fprintf(stderr, "get_theory_changes WARNING: Found unexpected axiom.\n");
			}
			break;

		case ProofType::IMPLICATION_INTRODUCTION:
			if (reference_count != 0) return true;
			formula = implication_intro_nodes[index_of(implication_intro_nodes, &proof)].key;
			if (!changes.add(change_type::IMPLICATION_INTRO_NODE, {formula, &proof})
			 || !discharged_axioms.add(proof.operands[1]))
				return false;
			break;

		case ProofType::PROOF_BY_CONTRADICTION:
			if (reference_count != 0) return true;
			discharged_axioms.add(proof.operands[1]);
			if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
			  && proof.operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
			  && proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				formula = negated_conjunction_nodes[index_of(negated_conjunction_nodes, &proof)].key;
				if (!visit_negated_conjunction(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::NEGATED_CONJUNCTION_NODE, {formula, &proof})) return false;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::UNIVERSAL_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				formula = existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key;
				if (!visit_negated_universal_intro(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, {formula, &proof})) return false;
			} else if (proof.operands[0]->type == ProofType::DISJUNCTION_ELIMINATION
					&& proof.operands[0]->operands[0] == proof.operands[1])
			{
				break;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::IMPLICATION_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				break;
			} else {
				fprintf(stderr, "get_theory_changes WARNING: Found unexpected proof by contradiction step.\n");
			}
			break;

		case ProofType::EXISTENTIAL_INTRODUCTION:
			if (reference_count != 0) return true;
			formula = existential_intro_nodes[index_of(existential_intro_nodes, &proof)].key;
			if (!visit_existential_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, {formula, &proof})) return false;
			break;

		case ProofType::DISJUNCTION_INTRODUCTION:
		case ProofType::DISJUNCTION_INTRODUCTION_LEFT:
		case ProofType::DISJUNCTION_INTRODUCTION_RIGHT:
			if (reference_count != 0) return true;
			formula = disjunction_intro_nodes[index_of(disjunction_intro_nodes, &proof)].key;
			if (!visit_disjunction_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::DISJUNCTION_INTRO_NODE, {formula, &proof})) return false;
			break;

		case ProofType::CONJUNCTION_INTRODUCTION:
		case ProofType::BETA_EQUIVALENCE:
		case ProofType::CONJUNCTION_ELIMINATION:
		case ProofType::CONJUNCTION_ELIMINATION_LEFT:
		case ProofType::CONJUNCTION_ELIMINATION_RIGHT:
		case ProofType::DISJUNCTION_ELIMINATION:
		case ProofType::IMPLICATION_ELIMINATION:
		case ProofType::BICONDITIONAL_INTRODUCTION:
		case ProofType::BICONDITIONAL_ELIMINATION_LEFT:
		case ProofType::BICONDITIONAL_ELIMINATION_RIGHT:
		case ProofType::NEGATION_ELIMINATION:
		case ProofType::FALSITY_ELIMINATION:
		case ProofType::UNIVERSAL_INTRODUCTION:
		case ProofType::UNIVERSAL_ELIMINATION:
		case ProofType::EXISTENTIAL_ELIMINATION:
		case ProofType::EQUALITY_ELIMINATION:
		case ProofType::COMPARISON_INTRODUCTION:
		case ProofType::PARAMETER:
		case ProofType::TERM_PARAMETER:
		case ProofType::ARRAY_PARAMETER:
		case ProofType::FORMULA_PARAMETER:
		case ProofType::COUNT:
			if (reference_count != 0) return true;
			break;
		}

		unsigned int operand_count;
		Proof* const* operands;
		proof.get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			if (!reference_counts.ensure_capacity(reference_counts.size + 1))
				return false;
			unsigned int index = reference_counts.index_of(operands[i]);
			if (index == reference_counts.size) {
				reference_counts.keys[index] = operands[i];
				reference_counts.values[index] = operands[i]->reference_count;
				reference_counts.size++;
			}
			reference_counts.values[index]--;
			if (!get_theory_changes(*operands[i], discharged_axioms, reference_counts, changes, std::forward<Visitor>(visitor)...))
				return false;
		}
		discharged_axioms.length = old_discharged_axiom_count;
		return true;
	}

	template<typename... Visitor>
	inline bool get_theory_changes(Proof& proof, theory::changes& changes, Visitor&&... visitor) const
	{
		array<const Proof*> discharged_axioms(16);
		array_map<const Proof*, unsigned int> reference_counts(32);
		reference_counts.keys[0] = &proof;
		reference_counts.values[0] = proof.reference_count - 1;
		reference_counts.size++;
		return get_theory_changes(proof, discharged_axioms, reference_counts, changes, std::forward<Visitor>(visitor)...);
	}

	static bool unify_atom(Term* term, Term* atom,
			array<Formula*>& quantifiers,
			array_map<Formula*, Term*>& unifications)
	{
		if (term->type == TermType::VARIABLE) {
			if (!unifications.ensure_capacity(unifications.size + 1))
				return false;
			Formula* quantifier = quantifiers[term->variable - 1];
			unsigned int index = unifications.index_of(quantifier);
			if (index == unifications.size) {
				unifications.values[unifications.size] = atom;
				unifications.keys[unifications.size++] = quantifier;
			} else {
				if (unifications.values[index] != atom && *unifications.values[index] != *atom)
					return false;
			}
			return true;
		} else if (term->type != atom->type) {
			return false;
		}

		switch (term->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
			return true;
		case FormulaType::CONSTANT:
			return term->constant == atom->constant;
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
			return term->variable == atom->variable;
		case FormulaType::PARAMETER:
			return term->parameter == atom->parameter;
		case FormulaType::INTEGER:
			return term->integer == atom->integer;
		case FormulaType::STRING:
			return term->str == atom->str;
		case FormulaType::UINT_LIST:
			return term->uint_list == atom->uint_list;
		case FormulaType::EQUALS:
		case FormulaType::UNARY_APPLICATION:
			return unify_atom(term->binary.left, atom->binary.left, quantifiers, unifications)
				&& unify_atom(term->binary.right, atom->binary.right, quantifiers, unifications);
		case FormulaType::BINARY_APPLICATION:
			return unify_atom(term->ternary.first, atom->ternary.first, quantifiers, unifications)
				&& unify_atom(term->ternary.second, atom->ternary.second, quantifiers, unifications)
				&& unify_atom(term->ternary.third, atom->ternary.third, quantifiers, unifications);
		case TermType::NOT:
			return unify_atom(term->unary.operand, atom->unary.operand, quantifiers, unifications);
		case TermType::IF_THEN:
		case TermType::AND:
		case TermType::OR:
		case TermType::IFF:
		case TermType::FOR_ALL:
		case TermType::EXISTS:
		case TermType::LAMBDA:
			return false;
		case TermType::ANY:
		case TermType::ANY_RIGHT:
		case TermType::ANY_RIGHT_ONLY:
		case TermType::ANY_ARRAY:
		case TermType::ANY_CONSTANT:
		case TermType::ANY_CONSTANT_EXCEPT:
		case TermType::ANY_QUANTIFIER:
			break;
		}
		fprintf(stderr, "theory.unify_atom ERROR: Unrecognized TermType.\n");
		return false;
	}

	static bool unify(Formula* first, Formula* second,
			array<Formula*>& quantifiers,
			array_map<Formula*, Term*>& unifications)
	{
		if (first->type == FormulaType::VARIABLE) {
			if (!unifications.ensure_capacity(unifications.size + 1))
				return false;
			bool contains;
			Term*& unification = unifications.get(quantifiers[first->variable - 1], contains);
			if (!contains) {
				unification = second;
				unifications.keys[unifications.size++] = quantifiers[first->variable - 1];
			} else {
				if (unification != second && *unification != *second)
					return false;
			}
			return true;
		} else if (first->type != second->type) {
			return false;
		}

		switch (first->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
			return true;
		case FormulaType::CONSTANT:
			return first->constant == second->constant;
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
			return first->variable == second->variable;
		case FormulaType::PARAMETER:
			return first->parameter == second->parameter;
		case FormulaType::INTEGER:
			return first->integer == second->integer;
		case FormulaType::STRING:
			return first->str == second->str;
		case FormulaType::UINT_LIST:
			return first->uint_list == second->uint_list;
		case FormulaType::IF_THEN:
		case FormulaType::EQUALS:
		case FormulaType::UNARY_APPLICATION:
			return unify(first->binary.left, second->binary.left, quantifiers, unifications)
				&& unify(first->binary.right, second->binary.right, quantifiers, unifications);
		case FormulaType::BINARY_APPLICATION:
			return unify(first->ternary.first, second->ternary.first, quantifiers, unifications)
				&& unify(first->ternary.second, second->ternary.second, quantifiers, unifications)
				&& unify(first->ternary.third, second->ternary.third, quantifiers, unifications);
		case FormulaType::NOT:
			return unify(first->unary.operand, second->unary.operand, quantifiers, unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			if (first->array.length != second->array.length) return false;
			for (unsigned int i = 0; i < first->array.length; i++)
				if (!unify(first->array.operands[i], second->array.operands[i], quantifiers, unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!quantifiers.add(first)
			 || !unify(first->quantifier.operand, second->quantifier.operand, quantifiers, unifications))
				return false;
			quantifiers.length--;
			return true;
		case FormulaType::ANY:
		case FormulaType::ANY_RIGHT:
		case FormulaType::ANY_RIGHT_ONLY:
		case FormulaType::ANY_ARRAY:
		case FormulaType::ANY_CONSTANT:
		case FormulaType::ANY_CONSTANT_EXCEPT:
		case FormulaType::ANY_QUANTIFIER:
			break;
		}
		fprintf(stderr, "theory.unify ERROR: Unrecognized FormulaType.\n");
		return false;
	}

	static bool unify(Formula* first, Formula* second,
			array<Formula*>& first_quantifiers,
			array_map<Formula*, Term*>& first_unifications,
			array<Formula*>& second_quantifiers,
			array_map<Formula*, Term*>& second_unifications)
	{
		if (first->type == FormulaType::VARIABLE) {
			if (!first_unifications.ensure_capacity(first_unifications.size + 1))
				return false;
			Formula* first_quantifier = first_quantifiers[first->variable - 1];
			unsigned int index = first_unifications.index_of(first_quantifier);
			if (index == first_unifications.size) {
				first_unifications.values[first_unifications.size] = second;
				first_unifications.keys[first_unifications.size++] = first_quantifier;
			} else {
				if (first_unifications.values[index] != second && *first_unifications.values[index] != *second)
					return false;
			}
			return true;
		} else if (second->type == FormulaType::VARIABLE) {
			if (!second_unifications.ensure_capacity(second_unifications.size + 1))
				return false;
			Formula* second_quantifier = second_quantifiers[second->variable - 1];
			unsigned int index = second_unifications.index_of(second_quantifier);
			if (index == second_unifications.size) {
				second_unifications.values[second_unifications.size] = first;
				second_unifications.keys[second_unifications.size++] = second_quantifier;
			} else {
				if (second_unifications.values[index] != first && *second_unifications.values[index] != *first)
					return false;
			}
			return true;
		} else if (first->type != second->type) {
			return false;
		}

		switch (first->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
			return true;
		case FormulaType::CONSTANT:
			return first->constant == second->constant;
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
			return first->variable == second->variable;
		case FormulaType::PARAMETER:
			return first->parameter == second->parameter;
		case FormulaType::INTEGER:
			return first->integer == second->integer;
		case FormulaType::STRING:
			return first->str == second->str;
		case FormulaType::UINT_LIST:
			return first->uint_list == second->uint_list;
		case FormulaType::IF_THEN:
		case FormulaType::EQUALS:
		case FormulaType::UNARY_APPLICATION:
			return unify(first->binary.left, second->binary.left, first_quantifiers, first_unifications, second_quantifiers, second_unifications)
				&& unify(first->binary.right, second->binary.right, first_quantifiers, first_unifications, second_quantifiers, second_unifications);
		case FormulaType::BINARY_APPLICATION:
			return unify(first->ternary.first, second->ternary.first, first_quantifiers, first_unifications, second_quantifiers, second_unifications)
				&& unify(first->ternary.second, second->ternary.second, first_quantifiers, first_unifications, second_quantifiers, second_unifications)
				&& unify(first->ternary.third, second->ternary.third, first_quantifiers, first_unifications, second_quantifiers, second_unifications);
		case FormulaType::NOT:
			return unify(first->unary.operand, second->unary.operand, first_quantifiers, first_unifications, second_quantifiers, second_unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			if (first->array.length != second->array.length) return false;
			for (unsigned int i = 0; i < first->array.length; i++)
				if (!unify(first->array.operands[i], second->array.operands[i], first_quantifiers, first_unifications, second_quantifiers, second_unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!first_quantifiers.add(first) || !second_quantifiers.add(second)
			 || !unify(first->quantifier.operand, second->quantifier.operand, first_quantifiers, first_unifications, second_quantifiers, second_unifications))
				return false;
			first_quantifiers.length--;
			second_quantifiers.length--;
			return true;
		case FormulaType::ANY:
		case FormulaType::ANY_RIGHT:
		case FormulaType::ANY_RIGHT_ONLY:
		case FormulaType::ANY_ARRAY:
		case FormulaType::ANY_CONSTANT:
		case FormulaType::ANY_CONSTANT_EXCEPT:
		case FormulaType::ANY_QUANTIFIER:
			break;
		}
		fprintf(stderr, "theory.unify ERROR: Unrecognized FormulaType.\n");
		return false;
	}

	template<bool Contradiction>
	bool is_provable_by_theorem_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values) const
	{
		unsigned int old_size = new_possible_values.length;
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == nullptr) continue;

			array<Formula*> second_quantifiers(sets.sets[i].arity);
			Formula* set_formula = sets.sets[i].size_axiom->formula->binary.left->binary.right;
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				second_quantifiers.add(set_formula);
				set_formula = set_formula->quantifier.operand;
			}
			if (Contradiction) {
				if (set_formula->type != FormulaType::NOT) continue;
				set_formula = set_formula->unary.operand;
			}
			array_map<Formula*, Term*> first_unifications(4);
			array_map<Formula*, Term*> second_unifications(4);
			if (!unify(formula, set_formula, quantifiers, first_unifications, second_quantifiers, second_unifications))
				continue;

			hash_set<tuple> provable_elements(16);
			if (!sets.get_provable_elements(i, provable_elements))
				return false;

			for (unsigned int j = 0; j < possible_values.length; j++) {
				const instantiation_tuple& values = possible_values[j];
				if (!new_possible_values.ensure_capacity(new_possible_values.length + provable_elements.size + 1)) {
					for (tuple& tup : provable_elements) free(tup);
					return false;
				}
				instantiation_tuple new_value(values);

				bool unifies = true;
				for (const auto& unification : first_unifications) {
					if (unification.value->type != TermType::VARIABLE
					 && !new_value.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
						unifies = false;
						break;
					}
				}
				if (!unifies) continue;

				for (const tuple& element : provable_elements) {
					instantiation_tuple& new_new_value = new_possible_values[new_possible_values.length];
					if (!init(new_new_value, new_value)) {
						for (tuple& tup : provable_elements) free(tup);
						return false;
					}

					unifies = true;
					for (const auto& unification : first_unifications) {
						if (unification.value->type != TermType::VARIABLE) continue;
						if (unification.value->variable > element.length) {
							if (unification.value->variable != unification.key->quantifier.variable) {
								unifies = false;
								break;
							} else {
								continue;
							}
						}
						Term* constant = Term::new_constant(element[unification.value->variable - 1]);
						if (constant == nullptr) {
							for (tuple& tup : provable_elements) free(tup);
							return false;
						}
						if (!new_new_value.unify_value(unification.key->quantifier.variable - 1, constant)) {
							free(*constant); free(constant);
							unifies = false;
							break;
						}
						free(*constant); free(constant);
					}
					if (!unifies) {
						free(new_new_value);
						continue;
					}

					for (const auto& unification : second_unifications) {
						if (unification.value->type != TermType::CONSTANT
						 || unification.value->constant != element[unification.key->quantifier.variable - 1])
						{
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						free(new_new_value);
						continue;
					}

					new_possible_values.length++;
				}
			}
			for (tuple& tup : provable_elements) free(tup);
		}
		if (new_possible_values.length > old_size) {
			insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
			array<instantiation_tuple> temp(new_possible_values.length);
			set_union(temp.data, temp.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
			swap(temp, new_possible_values);
			for (auto& element : temp) free(element);
		}
		return true;
	}

	bool for_all_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values) const
	{
		unsigned int arity = 0;
		while (formula->type == TermType::FOR_ALL) {
			formula = formula->quantifier.operand;
			arity++;
		}

		if (formula->type == TermType::IF_THEN) {
			Formula* antecedent = formula->binary.left;
			Formula* consequent = formula->binary.right;

			/* find sets that unify with `antecedent` */
			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axiom == nullptr || sets.sets[i].arity < arity) continue;

				array<Formula*> second_quantifiers(sets.sets[i].arity);
				Formula* set_formula = sets.sets[i].size_axiom->formula->binary.left->binary.right;
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					second_quantifiers.add(set_formula);
					set_formula = set_formula->quantifier.operand;
				}
				array_map<Formula*, Term*> first_unifications(4);
				array_map<Formula*, Term*> second_unifications(4);
				if (!unify(set_formula, antecedent, second_quantifiers, second_unifications, quantifiers, first_unifications))
					continue;

				hash_set<tuple> provable_elements(16);
				if (!sets.get_provable_elements(i, provable_elements))
					return false;
				if (provable_elements.size != sets.sets[i].set_size) {
					for (tuple& tup : provable_elements) free(tup);
					continue;
				}

				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& values = possible_values[i];
					array<instantiation_tuple> temp_possible_values(4);
					instantiation_tuple& new_values = temp_possible_values[0];
					if (!init(new_values, values)) {
						for (tuple& tup : provable_elements) free(tup);
						return false;
					}
					temp_possible_values.length = 1;

					bool unifies = true;
					for (const auto& unification : first_unifications) {
						if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						for (auto& element : temp_possible_values) free(element);
						continue;
					}

					/* check that `consequent` is provable for values in `provable_elements` */
					for (const tuple& element : provable_elements) {
						Term** src_variables = (Term**) alloca(sizeof(Term*) * (element.length + second_unifications.size) * 2);
						for (unsigned int i = 0; i < element.length; i++) {
							src_variables[i] = Term::new_variable(i + 1);
							if (src_variables[i] == nullptr) {
								for (unsigned int j = 0; j < i; j++) { free(*src_variables[j]); free(src_variables[j]); }
								for (auto& element : temp_possible_values) free(element);
								for (tuple& tup : provable_elements) free(tup);
								return false;
							}
						} for (unsigned int i = 0; i < second_unifications.size; i++) {
							src_variables[element.length + i] = Term::new_variable(second_unifications.keys[i]->quantifier.variable);
							if (src_variables[element.length + i] == nullptr) {
								for (unsigned int j = 0; j < element.length + i; j++) { free(*src_variables[j]); free(src_variables[j]); }
								for (auto& element : temp_possible_values) free(element);
								for (tuple& tup : provable_elements) free(tup);
								return false;
							}
						}

						Term** dst_constants = src_variables + element.length + second_unifications.size;
						for (unsigned int i = 0; i < element.length; i++) {
							dst_constants[i] = Term::new_constant(element[i]);
							if (dst_constants[i] == nullptr) {
								for (unsigned int j = 0; j < element.length; j++) { free(*src_variables[j]); free(src_variables[j]); }
								for (unsigned int j = 0; j < i; j++) { free(*dst_constants[j]); free(dst_constants[j]); }
								for (auto& element : temp_possible_values) free(element);
								for (tuple& tup : provable_elements) free(tup);
								return false;
							}
						} for (unsigned int i = 0; i < second_unifications.size; i++) {
							dst_constants[element.length + i] = second_unifications.values[i];
							dst_constants[element.length + i]->reference_count++;
						}

						Formula* substituted_consequent = substitute_all(consequent, src_variables, dst_constants, element.length);
						for (unsigned int j = 0; j < element.length + second_unifications.size; j++) { free(*src_variables[j]); free(src_variables[j]); }
						for (unsigned int j = 0; j < element.length + second_unifications.size; j++) { free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) free(dst_constants[j]); }

						if (!is_provable_without_abduction<false>(substituted_consequent, quantifiers, temp_possible_values)) {
							free(*substituted_consequent); if (substituted_consequent->reference_count == 0) free(substituted_consequent);
							break;
						}
						free(*substituted_consequent); if (substituted_consequent->reference_count == 0) free(substituted_consequent);
					}

					if (temp_possible_values.length != 0) {
						array<instantiation_tuple> temp(new_possible_values.length + temp_possible_values.length);
						set_union(temp, new_possible_values, temp_possible_values);
						swap(temp, new_possible_values);
						for (auto& element : temp_possible_values) free(element);
						for (auto& element : temp) free(element);
					}
				}
				for (tuple& tup : provable_elements) free(tup);
			}
		}
		return true;
	}

	template<bool Contradiction>
	bool exists_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values) const
	{
		Formula* quantified = formula->quantifier.operand;
		unsigned int variable = formula->quantifier.variable;

		array<instance> constants(ground_concept_capacity + 1);
		hash_set<int64_t> integers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != NULL) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length++].constant = new_constant_offset + i;

				for (const auto& entry : ground_concepts[i].function_values) {
					Term* constant = entry.value->formula->binary.right;
					if (constant->type == TermType::INTEGER) {
						if (!integers.add(constant->integer)) return false;
					} else if (constant->type == TermType::STRING) {
						bool contains = false;
						for (const string* str : strings)
							if (str == &constant->str || *str == constant->str) { contains = true; break; }
						if (!contains && !strings.add(&constant->str)) return false;
					}
				}
			}
		}
		if (!constants.ensure_capacity(constants.length + integers.size + strings.length))
			return false;
		for (int64_t integer : integers) {
			constants[constants.length].type = instance_type::INTEGER;
			constants[constants.length++].integer = integer;
		} for (string* str : strings) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length++].str = str;
		}

		if (!filter_constants(*this, quantified, variable, constants))
			return true;

		Term* var = Term::new_variable(variable);
		if (var == nullptr)
			return false;

		for (const instance& id : constants) {
			Term* constant = nullptr;
			if (id.type == instance_type::CONSTANT) {
				constant = Formula::new_constant(id.constant);
			} else if (id.type == instance_type::INTEGER) {
				constant = Formula::new_int(id.integer);
			} else if (id.type == instance_type::STRING) {
				constant = Formula::new_string(*id.str);
			}
			if (constant == nullptr) {
				free(*var); if (var->reference_count == 0) free(var);
				return false;
			}
			Formula* substituted = substitute(quantified, var, constant);
			if (substituted == NULL) {
				free(*var); if (var->reference_count == 0) free(var);
				free(*constant); if (constant->reference_count == 0) free(constant);
				return false;
			}

			array<instantiation_tuple> copy(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!init(copy[copy.length], possible_values[j])) {
					for (auto& element : copy) free(element);
					free(*var); if (var->reference_count == 0) free(var);
					return false;
				}
				copy.length++;
			}

			bool result = is_provable_without_abduction<Contradiction>(substituted, quantifiers, copy);
			free(*substituted); if (substituted->reference_count == 0) free(substituted);
			free(*constant); if (constant->reference_count == 0) free(constant);
			if (!result) {
				for (auto& element : copy) free(element);
				continue;
			}

			array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
			set_union(union_result, new_possible_values, copy);
			swap(union_result, new_possible_values);
			for (auto& element : union_result) free(element);
			for (auto& element : copy) free(element);
		}
		free(*var); if (var->reference_count == 0) free(var);
		return true;
	}

	template<bool Contradiction>
	bool is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			array<instantiation_tuple>& possible_values) const
	{
		/* consider if this formula provable by other extensional edges (i.e. universally-quantified theorems) */
		array<instantiation_tuple> new_possible_values(possible_values.length);
		if (!is_provable_by_theorem_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values)) {
			for (auto& element : new_possible_values) free(element);
			return false;
		}

		if (formula->type == FormulaType::UNARY_APPLICATION) {
			if (formula->binary.right->type != TermType::CONSTANT
			 && formula->binary.right->type != TermType::VARIABLE)
			{
				for (auto& element : new_possible_values) free(element);
				for (auto& element : possible_values) free(element);
				possible_values.clear(); return false;
			}

			for (const auto& atom : atoms) {
				array_map<Formula*, Term*> unifications(2);
				if (!unify_atom(formula->binary.left, atom.key.binary.left, quantifiers, unifications))
					continue;
				for (unsigned int i = 0; i < possible_values.length; i++) {
					const array<unsigned int>& constants = (Contradiction ? atom.value.value : atom.value.key);
					array<instantiation_tuple> temp(possible_values.length + constants.length);
					const instantiation_tuple& values = possible_values[i];
					instantiation_tuple new_values(values);
					/* make sure the `possible_values` unify with the variables in this term */
					bool unifies = true;
					for (const auto& unification : unifications) {
						if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						for (auto& element : temp) free(element);
						continue;
					}

					if (formula->binary.right->type == TermType::VARIABLE) {
						for (unsigned int constant : constants) {
							instantiation_tuple& current_new_values = temp[temp.length];
							if (!init(current_new_values, new_values)) {
								for (auto& element : new_possible_values) free(element);
								for (auto& element : temp) free(element);
								return false;
							}
							Term* term = Term::new_constant(constant);
							if (term == nullptr) {
								for (auto& element : new_possible_values) free(element);
								for (auto& element : temp) free(element);
								return false;
							}
							unifies = current_new_values.unify_value(formula->binary.right->variable - 1, term);
							free(*term); free(term);
							if (!unifies) {
								free(current_new_values);
								continue;
							}
							temp.length++;
						}
						if (temp.length > 0)
							insertion_sort(temp, default_sorter());
					} else if (constants.contains(formula->binary.right->constant)) {
						if (!init(temp[0], new_values)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						temp.length++;
					}
					if (temp.length != 0) {
						array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
						set_union(union_result, new_possible_values, temp);
						swap(new_possible_values, union_result);
						for (auto& element : temp) free(element);
						for (auto& element : union_result) free(element);
					}
				}
			}
			for (auto& element : possible_values) free(element);
			swap(new_possible_values, possible_values);

		} else if (formula->type == FormulaType::BINARY_APPLICATION) {
			if (formula->ternary.first->constant == (unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL) {
				array<instantiation_tuple> temp(possible_values.length);
				if (formula->ternary.second->type == TermType::VARIABLE) {
					for (unsigned int i = 0; i < possible_values.length; i++) {
						if (!init(temp[temp.length], possible_values[i])) {
							for (auto& element : new_possible_values) free(element);
							for (auto& element : possible_values) free(element);
							possible_values.clear(); return false;
						} else if (!Contradiction && !temp[temp.length].unify_greater_than_or_equal(formula->ternary.second->variable - 1, formula->ternary.third)) {
							free(temp[temp.length]); continue;
						} else if (Contradiction
								&& (!temp[temp.length].unify_less_than_or_equal(formula->ternary.second->variable - 1, formula->ternary.third)
								 || !temp[temp.length].antiunify_value(formula->ternary.second->variable - 1, formula->ternary.third)))
						{
							free(temp[temp.length]); continue;
						}
						temp.length++;
					}
				} else if (formula->ternary.third->type == TermType::VARIABLE) {
					for (unsigned int i = 0; i < possible_values.length; i++) {
						if (!init(temp[temp.length], possible_values[i])) {
							for (auto& element : new_possible_values) free(element);
							for (auto& element : possible_values) free(element);
							possible_values.clear(); return false;
						} else if (!Contradiction && !temp[temp.length].unify_less_than_or_equal(formula->ternary.third->variable - 1, formula->ternary.second)) {
							free(temp[temp.length]); continue;
						} else if (Contradiction
								&& (!temp[temp.length].unify_greater_than_or_equal(formula->ternary.third->variable - 1, formula->ternary.second)
								 || !temp[temp.length].antiunify_value(formula->ternary.third->variable - 1, formula->ternary.second)))
						{
							free(temp[temp.length]); continue;
						}
						temp.length++;
					}
				} else if (formula->ternary.second->type == TermType::INTEGER && formula->ternary.third->type == TermType::INTEGER) {
					if (!Contradiction && formula->ternary.second->integer >= formula->ternary.third->integer) {
						for (auto& element : new_possible_values) free(element);
						return true;
					} else if (Contradiction && formula->ternary.second->integer < formula->ternary.third->integer) {
						for (auto& element : new_possible_values) free(element);
						return true;
					} else {
						for (auto& element : new_possible_values) free(element);
						for (auto& element : possible_values) free(element);
						possible_values.clear(); return false;
					}
				} else {
					for (auto& element : new_possible_values) free(element);
					for (auto& element : possible_values) free(element);
					possible_values.clear(); return false;
				}
				if (temp.length != 0) {
					sort(temp);
					array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
					set_union(union_result, new_possible_values, temp);
					for (auto& element : temp) free(element);
					for (auto& element : new_possible_values) free(element);
					for (auto& element : possible_values) free(element);
					swap(union_result, possible_values);
				} else {
					for (auto& element : possible_values) free(element);
					swap(new_possible_values, possible_values);
				}
				return possible_values.length != 0;
			}

			unsigned int predicate;
			if (formula->ternary.first->type == TermType::VARIABLE)
				predicate = 0;
			else if (formula->ternary.first->type == TermType::CONSTANT)
				predicate = formula->ternary.first->constant;
			else {
				for (auto& element : new_possible_values) free(element);
				for (auto& element : possible_values) free(element);
				possible_values.clear(); return false;
			}

			unsigned int arg1;
			if (formula->ternary.second->type == TermType::VARIABLE)
				arg1 = 0;
			else if (formula->ternary.second->type == TermType::CONSTANT)
				arg1 = formula->ternary.second->constant;
			else {
				for (auto& element : new_possible_values) free(element);
				for (auto& element : possible_values) free(element);
				possible_values.clear(); return false;
			}

			unsigned int arg2;
			if (formula->ternary.third->type == TermType::VARIABLE)
				arg2 = 0;
			else if (formula->ternary.third->type == TermType::CONSTANT)
				arg2 = formula->ternary.third->constant;
			else {
				for (auto& element : new_possible_values) free(element);
				for (auto& element : possible_values) free(element);
				possible_values.clear(); return false;
			}

			array<instantiation_tuple> temp(possible_values.length);
			for (const auto& rel : relations) {
				if (rel.key.predicate != 0 && predicate != 0 && rel.key.predicate != predicate) continue;
				if (rel.key.arg1 != 0 && arg1 != 0 && rel.key.arg1 != arg1) continue;
				if (rel.key.arg2 != 0 && arg2 != 0 && rel.key.arg2 != arg2) continue;

				const array<unsigned int>& constants = (Contradiction ? rel.value.value : rel.value.key);
				for (unsigned int i = 0; i < possible_values.length; i++) {
					if (!temp.ensure_capacity(temp.length + constants.length + 1)) {
						for (auto& element : new_possible_values) free(element);
						for (auto& element : temp) free(element);
						return false;
					}
					instantiation_tuple values(possible_values[i]);
					if (predicate == 0 && rel.key.predicate != 0) {
						Term* temp = Term::new_constant(rel.key.predicate);
						bool unifies = values.unify_value(formula->ternary.first->variable - 1, temp);
						free(*temp); free(temp);
						if (!unifies) continue;
					} if (arg1 == 0 && rel.key.arg1 != 0) {
						Term* temp = Term::new_constant(rel.key.arg1);
						bool unifies = values.unify_value(formula->ternary.second->variable - 1, temp);
						free(*temp); free(temp);
						if (!unifies) continue;
					} if (arg2 == 0 && rel.key.arg2 != 0) {
						Term* temp = Term::new_constant(rel.key.arg2);
						bool unifies = values.unify_value(formula->ternary.third->variable - 1, temp);
						free(*temp); free(temp);
						if (!unifies) continue;
					}

					if (rel.key.predicate == 0) {
						if (predicate == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.first->variable - 1, term);
								free(*term); free(term);
								if (!unifies) {
									free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(predicate)) {
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								temp.length++;
							}
						}
					} else if (rel.key.arg1 == 0) {
						if (arg1 == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.second->variable - 1, term);
								free(*term); free(term);
								if (!unifies) {
									free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(arg1)) {
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								temp.length++;
							}
						}
					} else {
						if (arg2 == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!init(new_values, values)) {
									for (auto& element : temp) free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.first->variable - 1, term);
								free(*term); free(term);
								if (!unifies) {
									free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(arg2)) {
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									for (auto& element : temp) free(element);
									return false;
								}
								temp.length++;
							}
						}
					}
				}
			}
			if (temp.length != 0) {
				sort(temp); unique_and_cleanup(temp);
				array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
				set_union(union_result, new_possible_values, temp);
				for (auto& element : temp) free(element);
				for (auto& element : new_possible_values) free(element);
				for (auto& element : possible_values) free(element);
				swap(union_result, possible_values);
			} else {
				for (auto& element : possible_values) free(element);
				swap(new_possible_values, possible_values);
			}

		} else if (formula->type == FormulaType::AND) {
			if (Contradiction) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					array<instantiation_tuple> copy(possible_values.length);
					for (unsigned int j = 0; j < possible_values.length; j++) {
						if (!init(copy[copy.length], possible_values[j])) {
							for (auto& element : copy) free(element);
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						copy.length++;
					}
					if (!is_provable_without_abduction<true>(formula->array.operands[i], quantifiers, copy))
						continue;

					array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
					set_union(union_result, new_possible_values, copy);
					for (auto& element : new_possible_values) free(element);
					for (auto& element : copy) free(element);
					swap(union_result, new_possible_values);
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) free(element);
				return possible_values.length != 0;
			}

			for (unsigned int i = 0; i < formula->array.length; i++)
				if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values)) break;
			if (possible_values.length != 0) {
				array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) free(element);
				for (auto& element : new_possible_values) free(element);
				swap(possible_values, union_result);
			} else {
				swap(possible_values, new_possible_values);
			}

		} else if (formula->type == FormulaType::OR) {
			if (Contradiction) {
				for (unsigned int i = 0; i < formula->array.length; i++)
					if (!is_provable_without_abduction<true>(formula->array.operands[i], quantifiers, possible_values)) break;
				if (possible_values.length != 0) {
					array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) free(element);
					for (auto& element : new_possible_values) free(element);
					swap(possible_values, union_result);
				} else {
					swap(possible_values, new_possible_values);
				}
				return possible_values.length != 0;
			}

			for (unsigned int i = 0; i < formula->array.length; i++) {
				array<instantiation_tuple> copy(possible_values.length);
				for (unsigned int j = 0; j < possible_values.length; j++) {
					if (!init(copy[copy.length], possible_values[j])) {
						for (auto& element : copy) free(element);
						for (auto& element : new_possible_values) free(element);
						return false;
					}
					copy.length++;
				}
				if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, copy))
					continue;

				array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
				set_union(union_result, new_possible_values, copy);
				for (auto& element : new_possible_values) free(element);
				for (auto& element : copy) free(element);
				swap(union_result, new_possible_values);
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) free(element);

		} else if (formula->type == FormulaType::IF_THEN) {
			if (Contradiction) {
				if (is_provable_without_abduction<false>(formula->binary.left, quantifiers, possible_values))
					is_provable_without_abduction<true>(formula->binary.right, quantifiers, possible_values);
				if (possible_values.length != 0) {
					array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) free(element);
					for (auto& element : new_possible_values) free(element);
					swap(union_result, possible_values);
				} else {
					swap(possible_values, new_possible_values);
				}
				return possible_values.length != 0;
			}

			array<instantiation_tuple> first(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!init(first[first.length], possible_values[j])) {
					for (auto& element : first) free(element);
					for (auto& element : new_possible_values) free(element);
					return false;
				}
				first.length++;
			}
			if (!is_provable_without_abduction<true>(formula->binary.left, quantifiers, first))
				first.length = 0;

			array<instantiation_tuple> second(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!init(second[second.length], possible_values[j])) {
					for (auto& element : first) free(element);
					for (auto& element : second) free(element);
					for (auto& element : new_possible_values) free(element);
					return false;
				}
				second.length++;
			}
			if (!is_provable_without_abduction<false>(formula->binary.left, quantifiers, second))
				second.length = 0;

			array<instantiation_tuple> temp(max((size_t) 1, new_possible_values.length + first.length));
			set_union(temp, new_possible_values, first);
			for (auto& element : first) free(element);
			for (auto& element : new_possible_values) free(element);
			new_possible_values.length = 0;
			if (!new_possible_values.ensure_capacity(temp.length + second.length)) {
				for (auto& element : temp) free(element);
				for (auto& element : second) free(element);
				return false;
			}
			set_union(new_possible_values, temp, second);
			for (auto& element : temp) free(element);
			for (auto& element : second) free(element);
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) free(element);

		} else if (formula->type == FormulaType::NOT) {
			if (!is_provable_without_abduction<!Contradiction>(formula->unary.operand, quantifiers, possible_values)) {
				for (auto& element : new_possible_values) free(element);
				return false;
			}
			if (possible_values.length != 0) {
				array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) free(element);
				for (auto& element : new_possible_values) free(element);
				swap(union_result, possible_values);
			} else {
				swap(possible_values, new_possible_values);
			}
			return possible_values.length != 0;

		} else if (formula->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!exists_is_provable_without_abduction<true>(formula, quantifiers, possible_values, new_possible_values)) {
					for (auto& element : new_possible_values) free(element);
					return false;
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) free(element);
				return possible_values.length != 0;
			}

			if (!for_all_is_provable_without_abduction(formula, quantifiers, possible_values, new_possible_values)) {
				for (auto& element : new_possible_values) free(element);
				return false;
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) free(element);

		} else if (formula->type == FormulaType::EXISTS) {
			if (Contradiction) {
				/* TODO: implement this */
				fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
				return false;
			}

			if (!exists_is_provable_without_abduction<false>(formula, quantifiers, possible_values, new_possible_values)) {
				for (auto& element : new_possible_values) free(element);
				return false;
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) free(element);

		} else if (formula->type == FormulaType::EQUALS) {
			if (formula->binary.left->type == TermType::UNARY_APPLICATION
			 && formula->binary.left->binary.left->type == TermType::CONSTANT
			 && formula->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE)
			{
				Term* set_definition = formula->binary.left->binary.right;
				if (set_definition->type == TermType::LAMBDA) {
					/* size(^[x]:f(x))=n */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;

						if (formula->binary.right->type == TermType::INTEGER) {
							if ((!Contradiction && formula->binary.right->integer != sets.sets[i].set_size)
							 || (Contradiction && formula->binary.right->integer == sets.sets[i].set_size))
								continue;
						}

						array<Formula*> quantifiers(4);
						array_map<Formula*, Term*> unifications(4);
						if (!unify(set_definition, sets.sets[i].size_axiom->formula->binary.left->binary.right, quantifiers, unifications))
							continue;

						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = temp[temp.length];
							if (!init(new_values, values)) {
								for (auto& element : new_possible_values) free(element);
								for (auto& element : temp) free(element);
								return false;
							}

							if (formula->binary.right->type == TermType::VARIABLE) {
								if ((!Contradiction && !new_values.unify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right))
								 || (Contradiction && !new_values.antiunify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right)))
								{
									free(new_values);
									continue;
								}
							}

							/* make sure the `possible_values` unify with the variables in this term */
							bool unifies = true;
							for (const auto& unification : unifications) {
								if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
									unifies = false;
									break;
								}
							}
							if (!unifies) {
								free(new_values);
								continue;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							insertion_sort(temp, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) free(element);
							for (auto& element : temp) free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (set_definition->type == TermType::VARIABLE) {
					/* size(x)=n */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;
						Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type != TermType::UNARY_APPLICATION
						 || set_formula->binary.left->type != TermType::CONSTANT
						 || set_formula->binary.right->type != TermType::VARIABLE
						 || set_formula->binary.right->variable != 1)
							continue;

						if (formula->binary.right->type == TermType::INTEGER) {
							if ((!Contradiction && formula->binary.right->integer != sets.sets[i].set_size)
							 || (Contradiction && formula->binary.right->integer == sets.sets[i].set_size))
								continue;
						}

						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = temp[temp.length];
							if (!init(new_values, values)) {
								for (auto& element : new_possible_values) free(element);
								for (auto& element : temp) free(element);
								return false;
							}

							if (!new_values.unify_value(set_definition->variable - 1, set_formula->binary.left)) {
								free(new_values);
								continue;
							}

							if (formula->binary.right->type == TermType::VARIABLE) {
								if ((!Contradiction && !new_values.unify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right))
								 || (Contradiction && !new_values.antiunify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right)))
								{
									free(new_values);
									continue;
								}
							}
							temp.length++;
						}
						if (temp.length != 0) {
							insertion_sort(temp, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) free(element);
							for (auto& element : temp) free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (set_definition->type == TermType::CONSTANT && set_definition->constant >= new_constant_offset) {
					/* size(a)=n */
					Term* set_formula = Term::new_apply(set_definition, &Variables<1>::value);
					if (set_formula == nullptr)
						return false;
					set_definition->reference_count++;
					bool contains;
					unsigned int set_id = sets.set_ids.get(*set_formula, contains);
					if (!contains) {
						for (auto& element : possible_values) free(element);
						possible_values.clear(); return false;
					}

					if (formula->binary.right->type == TermType::INTEGER) {
						if ((!Contradiction && formula->binary.right->integer != sets.sets[set_id].set_size)
						 || (Contradiction && formula->binary.right->integer == sets.sets[set_id].set_size))
						{
							for (auto& element : possible_values) free(element);
							possible_values.clear(); return false;
						}
					}
					array<instantiation_tuple> temp(possible_values.length);
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const instantiation_tuple& values = possible_values[i];
						instantiation_tuple& new_values = temp[temp.length];
						if (!init(new_values, values)) {
							for (auto& element : new_possible_values) free(element);
							for (auto& element : temp) free(element);
							return false;
						}
						if (formula->binary.right->type == TermType::VARIABLE) {
							if ((!Contradiction && !new_values.unify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right))
							 || (Contradiction && !new_values.antiunify_value(formula->binary.right->variable - 1, sets.sets[i].size_axiom->formula->binary.right)))
							{
								free(new_values);
								for (auto& element : possible_values) free(element);
								possible_values.clear(); return false;
							}
						}
						temp.length++;
					}
					if (temp.length != 0) {
						insertion_sort(temp, default_sorter());
						array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
						set_union(union_result, new_possible_values, temp);
						for (auto& element : new_possible_values) free(element);
						for (auto& element : temp) free(element);
						swap(union_result, new_possible_values);
					}

				} else {
					fprintf(stderr, "theory.is_provable_without_abduction ERROR: Unsupported set size axiom.\n");
					return false;
				}

				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) free(element);
				return possible_values.length != 0;
			}

			/* A=A */
			array_map<Formula*, Term*> unifications(4);
			if (!Contradiction && unify(formula->binary.left, formula->binary.right, quantifiers, unifications)) {
				if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
					for (auto& element : new_possible_values) free(element);
					return false;
				}
				unsigned int old_size = new_possible_values.length;
				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& values = possible_values[i];
					instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
					if (!init(new_values, values)) {
						for (auto& element : new_possible_values) free(element);
						return false;
					}

					/* make sure the `possible_values` unify with the variables in this term */
					bool unifies = true;
					for (const auto& unification : unifications) {
						if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						free(new_values);
						continue;
					}
					new_possible_values.length++;
				}
				if (new_possible_values.length > old_size) {
					insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
					array<instantiation_tuple> union_result(new_possible_values.length);
					set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
					for (auto& element : new_possible_values) free(element);
					swap(union_result, new_possible_values);
				}
			}

			/* a=b */
			if (!Contradiction) {
				if (formula->binary.left->type == TermType::CONSTANT) {
					if (formula->binary.right->type == TermType::CONSTANT && formula->binary.right->constant == formula->binary.left->constant) {
						for (auto& element : possible_values) free(element);
						possible_values.clear(); return false;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!init(new_values, values)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}

							if (!new_values.antiunify_value(formula->binary.right->variable - 1, formula->binary.left)) {
								free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (formula->binary.left->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) free(element);
						return false;
					}
					unsigned int old_size = new_possible_values.length;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const instantiation_tuple& values = possible_values[i];
						instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
						if (!init(new_values, values)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}

						if (!new_values.antiunify_value(formula->binary.left->variable - 1, formula->binary.right)) {
							free(new_values);
							continue;
						}
						new_possible_values.length++;
					}
					if (new_possible_values.length > old_size) {
						insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
						array<instantiation_tuple> union_result(new_possible_values.length);
						set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
						for (auto& element : new_possible_values) free(element);
						swap(union_result, new_possible_values);
					}
				}
			}

			/* a=A or A=a (definition of a) */
			Term* left = formula->binary.left;
			Term* right = formula->binary.right;
			if (formula->binary.left->type == TermType::CONSTANT || formula->binary.right->type == TermType::CONSTANT) {
				if (left->type != TermType::CONSTANT)
					swap(left, right);

				if (Contradiction) {
					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						if (ground_concepts[i].types.keys == nullptr || i + new_constant_offset == left->constant)
							continue;
						for (Proof* definition : ground_concepts[i].definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;
							if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const instantiation_tuple& values = possible_values[i];
								instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									return false;
								}

								/* make sure the `possible_values` unify with the variables in this term */
								bool unifies = true;
								for (const auto& unification : unifications) {
									if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
										unifies = false;
										break;
									}
								}
								if (!unifies) {
									free(new_values);
									continue;
								}
								new_possible_values.length++;
							}
							if (new_possible_values.length > old_size) {
								insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
								array<instantiation_tuple> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
								for (auto& element : new_possible_values) free(element);
								swap(union_result, new_possible_values);
							}
						}
					}
				} else if (left->constant >= new_constant_offset) {
					for (Proof* definition : ground_concepts[left->constant - new_constant_offset].definitions) {
						array_map<Formula*, Term*> unifications(4);
						if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
							continue;
						if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!init(new_values, values)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}

							/* make sure the `possible_values` unify with the variables in this term */
							bool unifies = true;
							for (const auto& unification : unifications) {
								if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
									unifies = false;
									break;
								}
							}
							if (!unifies) {
								free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) free(element);
							swap(union_result, new_possible_values);
						}
					}
				}
			} else if (formula->binary.left->type == TermType::VARIABLE || formula->binary.right->type == TermType::VARIABLE) {
				if (formula->binary.left->type != TermType::VARIABLE)
					swap(left, right);

				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& values = possible_values[i];
					if (!Contradiction && values.values[left->variable - 1].type == instantiation_type::CONSTANT) {
						const concept<ProofCalculus>& c = ground_concepts[values.values[left->variable - 1].constant - new_constant_offset];
						if (!new_possible_values.ensure_capacity(new_possible_values.length + c.definitions.length)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (Proof* definition : c.definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;

							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!init(new_values, values)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}

							/* make sure the `possible_values` unify with the variables in this term */
							bool unifies = true;
							for (const auto& unification : unifications) {
								if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
									unifies = false;
									break;
								}
							}
							if (!unifies) {
								free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) free(element);
							swap(union_result, new_possible_values);
						}
						continue;
					} else if (!Contradiction && values.values[left->variable - 1].type != instantiation_type::ANY) {
						continue;
					}

					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						const concept<ProofCalculus>& c = ground_concepts[i];
						if (c.types.keys == nullptr)
							continue;
						unsigned int old_size = new_possible_values.length;
						if (!new_possible_values.ensure_capacity(new_possible_values.length + c.definitions.length)) {
							for (auto& element : new_possible_values) free(element);
							return false;
						}
						instantiation_tuple new_values(values);
						if ((!Contradiction && !new_values.unify_value(left->variable - 1, c.definitions[0]->formula->binary.left))
						 || (Contradiction && !new_values.antiunify_value(left->variable - 1, c.definitions[0]->formula->binary.left)))
							continue;

						for (Proof* definition : c.definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;

							instantiation_tuple& new_new_values = new_possible_values[new_possible_values.length];
							if (!init(new_new_values, new_values)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}

							/* make sure the `possible_values` unify with the variables in this term */
							bool unifies = true;
							for (const auto& unification : unifications) {
								if (!new_new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
									unifies = false;
									break;
								}
							}
							if (!unifies) {
								free(new_new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) free(element);
							swap(union_result, new_possible_values);
						}
					}
				}
			}

			/* a(b)=A or A=a(b) (function definition of a evaluated at b) */
			if ((left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT)
			 || (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT))
			{
				if (left->type != TermType::UNARY_APPLICATION || left->binary.left->type != TermType::CONSTANT)
					swap(left, right);
				
				if (left->binary.right->type == TermType::CONSTANT && left->binary.right->constant >= new_constant_offset) {
					bool contains;
					Proof* function_value = ground_concepts[left->binary.right->constant - new_constant_offset].function_values.get(left->binary.left->constant, contains);
					if (contains) {
						array_map<Formula*, Term*> unifications(4);
						if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
							if (Contradiction) {
								swap(new_possible_values, possible_values);
								for (auto& element : new_possible_values) free(element);
								return true;
							}
						} else {
							if (Contradiction) {
								for (auto& element : new_possible_values) free(element);
								for (auto& element : possible_values) free(element);
								possible_values.clear(); return false;
							}
							if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
								for (auto& element : new_possible_values) free(element);
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const instantiation_tuple& values = possible_values[i];
								instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
								if (!init(new_values, values)) {
									for (auto& element : new_possible_values) free(element);
									return false;
								}

								/* make sure the `possible_values` unify with the variables in this term */
								bool unifies = true;
								for (const auto& unification : unifications) {
									if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
										unifies = false;
										break;
									}
								}
								if (!unifies) {
									free(new_values);
									continue;
								}
								new_possible_values.length++;
							}
							if (new_possible_values.length > old_size) {
								insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
								array<instantiation_tuple> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
								for (auto& element : new_possible_values) free(element);
								swap(union_result, new_possible_values);
							}
						}
					}
				} else if (left->binary.right->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) free(element);
						return false;
					}
					unsigned int old_size = new_possible_values.length;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const instantiation_tuple& values = possible_values[i];
						if (!Contradiction && values.values[left->binary.right->variable - 1].type == instantiation_type::CONSTANT) {
							const concept<ProofCalculus>& c = ground_concepts[values.values[left->binary.right->variable - 1].constant - new_constant_offset];
							bool contains;
							Proof* function_value = c.function_values.get(left->binary.left->constant, contains);
							if (contains) {
								array_map<Formula*, Term*> unifications(4);
								if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
									if (Contradiction) {
										instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
										if (!init(new_values, values)) {
											for (auto& element : new_possible_values) free(element);
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									if (Contradiction) {
										fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
										return false;
									}
									instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
									if (!init(new_values, values)) {
										for (auto& element : new_possible_values) free(element);
										return false;
									}

									/* make sure the `possible_values` unify with the variables in this term */
									bool unifies = true;
									for (const auto& unification : unifications) {
										if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
											unifies = false;
											break;
										}
									}
									if (!unifies) {
										free(new_values);
										continue;
									}
									new_possible_values.length++;
								}
							}
							continue;
						} else if (!Contradiction && values.values[left->binary.right->variable - 1].type != instantiation_type::ANY) {
							continue;
						}

						for (unsigned int i = 0; i < ground_concept_capacity; i++) {
							const concept<ProofCalculus>& c = ground_concepts[i];
							if (c.types.keys == nullptr)
								continue;
							bool contains;
							Proof* function_value = c.function_values.get(left->binary.left->constant, contains);
							if (contains) {
								array_map<Formula*, Term*> unifications(4);
								if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
									for (auto& element : new_possible_values) free(element);
									return false;
								}
								if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
									if (Contradiction) {
										instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
										if (!init(new_values, values)) {
											for (auto& element : new_possible_values) free(element);
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									if (Contradiction) {
										fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
										return false;
									}
									instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
									if (!init(new_values, values)) {
										for (auto& element : new_possible_values) free(element);
										return false;
									}

									/* make sure the `possible_values` unify with the variables in this term */
									bool unifies = true;
									for (const auto& unification : unifications) {
										if (!new_values.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
											unifies = false;
											break;
										}
									}
									if (!unifies) {
										free(new_values);
										continue;
									}
									new_possible_values.length++;
								}
							}
						}
					}
					if (new_possible_values.length > old_size) {
						insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size);
						array<instantiation_tuple> union_result(new_possible_values.length);
						set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
						for (auto& element : new_possible_values) free(element);
						swap(union_result, new_possible_values);
					}
				}
			}

			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) free(element);

		} else {
			fprintf(stderr, "is_provable_without_abduction ERROR: Unsupported FormulaType.\n");
			return false;
		}
		return possible_values.length != 0;
	}

	static bool get_unifying_atoms(
			Formula* formula, Formula* atom,
			array<Formula*>& quantifiers,
			array<array_map<Formula*, Term*>>& unifications)
	{
		switch (formula->type) {
		case FormulaType::CONSTANT:
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
		case FormulaType::PARAMETER:
		case FormulaType::INTEGER:
		case FormulaType::STRING:
		case FormulaType::UINT_LIST:
		case FormulaType::TRUE:
		case FormulaType::FALSE:
		case FormulaType::EQUALS:
		case FormulaType::UNARY_APPLICATION:
		case FormulaType::BINARY_APPLICATION:
			if (!unifications.ensure_capacity(unifications.length + 1)
			 || !array_map_init(unifications[unifications.length], 2))
				return false;
			if (unify_atom(formula, atom, quantifiers, unifications[unifications.length]))  {
				unifications.length++;
			} else {
				free(unifications[unifications.length]);
			}
			return true;
		case FormulaType::NOT:
			if (!unifications.ensure_capacity(unifications.length + 2)
			 || !array_map_init(unifications[unifications.length], 2))
				return false;
			if (unify_atom(formula, atom, quantifiers, unifications[unifications.length]))  {
				unifications.length++;
			} else {
				free(unifications[unifications.length]);
			}

			if (!array_map_init(unifications[unifications.length], 2))
				return false;
			if (unify_atom(formula->unary.operand, atom, quantifiers, unifications[unifications.length]))  {
				unifications.length++;
			} else {
				free(unifications[unifications.length]);
			}
			return get_unifying_atoms(formula->unary.operand, atom, quantifiers, unifications);

		case FormulaType::IF_THEN:
			return get_unifying_atoms(formula->binary.left, atom, quantifiers, unifications)
				&& get_unifying_atoms(formula->binary.right, atom, quantifiers, unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			for (unsigned int i = 0; i < formula->array.length; i++)
				if (!get_unifying_atoms(formula->array.operands[i], atom, quantifiers, unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!quantifiers.add(formula)) return false;
			if (!get_unifying_atoms(formula->quantifier.operand, atom, quantifiers, unifications)) {
				quantifiers.length--;
				return false;
			} else {
				quantifiers.length--;
				return true;
			}
		case FormulaType::ANY:
		case FormulaType::ANY_RIGHT:
		case FormulaType::ANY_RIGHT_ONLY:
		case FormulaType::ANY_ARRAY:
		case FormulaType::ANY_CONSTANT:
		case FormulaType::ANY_CONSTANT_EXCEPT:
		case FormulaType::ANY_QUANTIFIER:
			break;
		}
		fprintf(stderr, "theory.get_unifying_atoms ERROR: Unrecognized FormulaType.\n");
		return false;
	}

	template<bool ResolveInconsistencies>
	bool check_set_membership(unsigned int set_id,
			array<instantiation_tuple>& possible_values,
			array<tuple>& new_elements,
			array<Formula*>& quantifiers)
	{
		/* check if any of the possible values of the quantified variables are newly provable members of this set */
		hash_set<tuple> provable_elements(16);
		if (!sets.get_provable_elements(set_id, provable_elements)) {
			for (tuple& tup : provable_elements) free(tup);
			for (auto& element : possible_values) free(element);
			return false;
		}
		unsigned int provable_set_size = provable_elements.size;
		for (unsigned int j = 0; j < possible_values.length; j++) {
			const instantiation_tuple& values = possible_values[j];
			tuple new_tuple;
			if (!init(new_tuple, values.length)) {
				for (tuple& tup : provable_elements) free(tup);
				for (auto& element : possible_values) free(element);
				return false;
			}
			bool all_constants = true;
			for (unsigned int k = 0; k < values.length; k++) {
				if (values.values[k].type != instantiation_type::CONSTANT) {
					all_constants = false;
					break;
				} else {
					new_tuple[k] = values.values[k].constant;
				}
			}
			if (!all_constants) {
				free(new_tuple);
				continue;
			}

			bool is_old = provable_elements.contains(new_tuple);
			free(new_tuple);
			if (is_old) {
				free(possible_values[j]);
				move(possible_values[possible_values.length - 1], possible_values[j]);
				possible_values.length--;
				j--;
			}
		}
		for (tuple& tup : provable_elements) free(tup);

		if (possible_values.length == 0)
			return true;

		/* all values in `possible_values` are newly provable */
		for (const auto& entry : sets.extensional_graph.vertices[set_id].parents) {
			unsigned int other_set = entry.key;
			Formula* other_formula = sets.sets[other_set].set_formula();

			array<instantiation_tuple> copy(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!init(copy[copy.length], possible_values[j])) {
					for (auto& element : possible_values) free(element);
					return false;
				}
				copy.length++;
			}

			if (is_provable_without_abduction<true>(other_formula, quantifiers, copy)) {
				/* we found a contradiction with a universally-quantified theorem */
				for (auto& element : possible_values) free(element);
				for (auto& element : copy) free(element);
				return false;
			}
			for (auto& element : copy) free(element);
		}

		/* make sure the set is big enough to contain the new elements */
		provable_set_size += possible_values.length;
		if (provable_set_size > sets.sets[set_id].set_size
		 && !sets.template get_size_axiom<ResolveInconsistencies>(set_id, provable_set_size))
		{
			for (auto& element : possible_values) free(element);
			return false;
		}

		/* add the newly provable elements to `new_elements` */
		if (!new_elements.ensure_capacity(new_elements.length + possible_values.length))
			return false;
		for (unsigned int j = 0; j < possible_values.length; j++) {
			if (!init(new_elements[new_elements.length], possible_values[0].length)) {
				for (auto& element : possible_values) free(element);
				return false;
			}

			bool all_constant = true;
			const instantiation_tuple& values = possible_values[j];
			for (unsigned int k = 0; k < values.length; k++) {
				if (values.values[k].type != instantiation_type::CONSTANT) {
					all_constant = false;
					free(new_elements[new_elements.length]);
					break;
				}
				new_elements[new_elements.length][k] = values.values[k].constant;
			}
			if (!all_constant) continue;
			new_elements.length++;
		}
		for (auto& element : possible_values) free(element);
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_set_membership_after_addition(Formula* new_atom, Args&&... visitor) {
		/* We first need to find all pairs of tuples and sets such that, with
		   the addition of `new_atom` as an axiom, the tuple must necessarily
		   belong to the set (and is not the case otherwise without
		   `new_atom`). To do this, find all sets that contain an atom that
		   unifies with `new_atom`. */
		array_map<unsigned int, array<tuple>> new_elements(4);
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == nullptr)
				continue;
			Formula* set_formula = sets.sets[i].size_axiom->formula->binary.left->binary.right;

			array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				quantifiers[quantifiers.length++] = set_formula;
				set_formula = set_formula->quantifier.operand;
			}
			array<array_map<Formula*, Term*>> unifications(4);
			if (!get_unifying_atoms(set_formula, new_atom, quantifiers, unifications)) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				}
				return false;
			}

			array<instantiation_tuple> possible_values(max((size_t) 1, 4 * unifications.length));
			for (const array_map<Formula*, Term*>& unification : unifications) {
				/* given `unification`, find the values of the quantified variables that make `set_formula` necessarily true */
				bool valid_unification = true;
				array_map<unsigned int, Term*> bound_variable_map(unification.size);
				array<instantiation_tuple> temp_possible_values(4);
				instantiation_tuple& values = temp_possible_values[0];
				if (!init(values, sets.sets[i].arity)) {
					for (auto& element : unifications) free(element);
					for (auto& element : possible_values) free(element);
					for (auto entry : new_elements) {
						for (tuple& tup : entry.value) free(tup);
						free(entry.value);
					}
					return false;
				}
				temp_possible_values.length++;
				for (const auto& entry : unification) {
					bound_variable_map.put(entry.key->quantifier.variable, entry.value);
					if (entry.key->quantifier.variable <= sets.sets[i].arity) {
						if (entry.value->type == TermType::CONSTANT) {
							free(values.values[entry.key->quantifier.variable - 1]);
							values.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
							values.values[entry.key->quantifier.variable - 1].constant = entry.value->constant;
						} else if (entry.value->type == TermType::INTEGER) {
							free(values.values[entry.key->quantifier.variable - 1]);
							values.values[entry.key->quantifier.variable - 1].type = instantiation_type::INTEGER;
							values.values[entry.key->quantifier.variable - 1].constant = entry.value->integer;
						} else if (entry.value->type == TermType::STRING) {
							free(values.values[entry.key->quantifier.variable - 1]);
							values.values[entry.key->quantifier.variable - 1].type = instantiation_type::STRING;
							values.values[entry.key->quantifier.variable - 1].str = &entry.value->str;
						} else {
							valid_unification = false;
							break;
						}
					}
				}
				if (!valid_unification) {
					for (auto& element : temp_possible_values) free(element);
					continue;
				}

				Formula* new_formula = substitute_quantified_variables(set_formula, bound_variable_map);
				if (new_formula == nullptr) {
					for (auto& element : unifications) free(element);
					for (auto& element : possible_values) free(element);
					for (auto& element : temp_possible_values) free(element);
					for (auto entry : new_elements) {
						for (tuple& tup : entry.value) free(tup);
						free(entry.value);
					}
					return false;
				}

				if (!is_provable_without_abduction<false>(new_formula, quantifiers, temp_possible_values)) {
					for (auto& element : temp_possible_values) free(element);
					free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
					continue;
				}
				free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);

				array<instantiation_tuple> union_result(possible_values.length + temp_possible_values.length);
				set_union(union_result, possible_values, temp_possible_values);
				for (auto& element : possible_values) free(element);
				for (auto& element : temp_possible_values) free(element);
				swap(union_result, possible_values);
			}
			for (auto& element : unifications) free(element);

			if (!new_elements.ensure_capacity(new_elements.size + 1)
			 || !array_init(new_elements.values[new_elements.size], 4)) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				}
				return false;
			}
			new_elements.keys[new_elements.size++] = i;
			if (!check_set_membership<ResolveInconsistencies>(i, possible_values, new_elements.values[new_elements.size - 1], quantifiers)) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				}
				return false;
			}
		}

		/* add the elements to the sets */
		array_map<unsigned int, array<Formula*>> new_unary_elements(4);
		if (new_atom->type == FormulaType::UNARY_APPLICATION && new_atom->binary.right->type == TermType::CONSTANT
		 && new_atom->binary.right->constant >= new_constant_offset)
		{
			if (!array_init(new_unary_elements.values[0], 8)) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				}
				return false;
			}
			new_unary_elements.values[0][0] = Formula::new_apply(new_atom->binary.left, &Variables<1>::value);
			if (new_unary_elements.values[0][0] == nullptr) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				}
				free(new_unary_elements.values[0]);
				return false;
			}
			new_atom->binary.left->reference_count++;
			Variables<1>::value.reference_count++;
			new_unary_elements.keys[0] = new_atom->binary.right->constant;
			new_unary_elements.values[0].length++;
			new_unary_elements.size = 1;
		}
		for (const auto& entry : new_elements) {
			for (const tuple& element : entry.value) {
				if (element.length == 1) {
					if (!new_unary_elements.ensure_capacity(new_unary_elements.size + 1)) {
						for (auto entry : new_elements) {
							for (tuple& tup : entry.value) free(tup);
							free(entry.value);
						} for (auto entry : new_unary_elements) free(entry.value);
						return false;
					}
					unsigned int index = new_unary_elements.index_of(element[0]);
					array<Formula*>& conjuncts = new_unary_elements.values[index];
					if (index == new_unary_elements.size) {
						if (!array_init(conjuncts, 8)) {
							for (auto entry : new_elements) {
								for (tuple& tup : entry.value) free(tup);
								free(entry.value);
							} for (auto entry : new_unary_elements) free(entry.value);
							return false;
						}
						new_unary_elements.keys[index] = element[0];
						new_unary_elements.size++;
					}

					Formula* set_formula = sets.sets[entry.key].set_formula();
					if (set_formula->type == FormulaType::AND) {
						if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
							for (auto entry : new_elements) {
								for (tuple& tup : entry.value) free(tup);
								free(entry.value);
							} for (auto entry : new_unary_elements) free(entry.value);
							return false;
						}
					} else if (!conjuncts.add(set_formula)) {
						for (auto entry : new_elements) {
							for (tuple& tup : entry.value) free(tup);
							free(entry.value);
						} for (auto entry : new_unary_elements) free(entry.value);
						return false;
					}

				} else if (!sets.sets[entry.key].add_element(element)) {
					for (auto entry : new_elements) {
						for (tuple& tup : entry.value) free(tup);
						free(entry.value);
					} for (auto entry : new_unary_elements) free(entry.value);
					return false;
				}
			}
		}

		array_map<unsigned int, Formula*> new_unary_subsets(max((size_t) 1, new_unary_elements.size));
		for (auto entry : new_unary_elements) {
			Formula* new_subset;
			if (entry.value.length == 1)
				new_subset = entry.value[0];
			else new_subset = Formula::new_and(make_array_view(entry.value.data, entry.value.length));
			if (new_subset == nullptr) {
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				} for (auto entry : new_unary_subsets) {
					free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
				} for (auto entry : new_unary_elements) free(entry.value);
				return false;
			}
			for (Formula* conjunct : entry.value)
				conjunct->reference_count++;
			new_unary_subsets.keys[new_unary_subsets.size] = entry.key;
			new_unary_subsets.values[new_unary_subsets.size++] = new_subset;
		}
		for (auto entry : new_unary_elements) free(entry.value);

		/* We need to call `move_element_to_subset` in the exact reverse order
		   in which we call `move_element_to_superset` in
		   `check_set_membership_after_subtraction`, so we sort
		   `new_unary_subsets`. Otherwise, the order depends on the order of
		   `sets.sets` which can change between here and the corresponding call
		   to `check_set_membership_after_subtraction`. */
		insertion_sort(new_unary_subsets.keys, new_unary_subsets.values, new_unary_subsets.size);
		for (unsigned int i = 0; i < new_unary_subsets.size; i++) {
			if (!move_element_to_subset<ResolveInconsistencies>(new_unary_subsets.keys[i], new_unary_subsets.values[i], std::forward<Args>(visitor)...)) {
				/* we have to undo the changes in reverse order */
				for (unsigned int j = i; j > 0; j--)
					move_element_to_superset(new_unary_subsets.keys[j - 1], new_unary_subsets.values[j - 1], std::forward<Args>(visitor)...);					
				for (auto entry : new_elements)
					for (const tuple& element : entry.value)
						if (element.length != 1) sets.sets[entry.key].remove_element_at(sets.sets[entry.key].index_of_element(element));
				for (auto entry : new_elements) {
					for (tuple& tup : entry.value) free(tup);
					free(entry.value);
				} for (auto entry : new_unary_elements) free(entry.value);
				return false;
			}
		}
		for (auto entry : new_elements) {
			for (tuple& tup : entry.value) free(tup);
			free(entry.value);
		} for (auto entry : new_unary_subsets) {
			free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
		}
		return true;
	}

	template<typename... Args>
	bool check_set_membership_after_subtraction(Formula* old_atom, Args&&... visitor)
	{
		bool changed = true;
		array_map<unsigned int, array<Formula*>> old_unary_conjuncts(4);
		array_map<unsigned int, unsigned int> removed_elements(8);
		while (changed) {
			changed = false;
			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axiom == nullptr)
					continue;
				Formula* set_formula = sets.sets[i].size_axiom->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}
				array<array_map<Formula*, Term*>> unifications(4);
				if (!get_unifying_atoms(set_formula, old_atom, quantifiers, unifications)) {
					for (auto entry : old_unary_conjuncts) free(entry.value);
					return false;
				}

				for (unsigned int j = 0; j < sets.sets[i].element_count(); j++) {
					unsigned int* element = (unsigned int*) alloca(sizeof(unsigned int) * sets.sets[i].arity);
					const unsigned int* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * j);
					for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
						element[k] = element_src[k];
					bool has_satisfying_unification = false;
					for (const array_map<Formula*, Term*>& unification : unifications) {
						bool matches = true;
						for (const auto& entry : unification) {
							if (entry.key->quantifier.variable > sets.sets[i].arity)
								continue;
							if (entry.value->type != TermType::CONSTANT
							 || element[entry.key->quantifier.variable - 1] != entry.value->constant)
							{
								matches = false;
								break;
							}
						}
						if (matches) {
							has_satisfying_unification = true;
							break;
						}
					}
					if (!has_satisfying_unification) continue;

					array<Formula*> conjuncts(4);
					if (set_formula->type == FormulaType::AND) {
						if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
							for (auto& element : unifications) free(element);
							for (auto entry : old_unary_conjuncts) free(entry.value);
							return false;
						}
					} else {
						conjuncts[conjuncts.length++] = set_formula;
					}

					sets.sets[i].remove_element_at(j);
					for (unsigned int l = 0; l < conjuncts.length; l++) {
						array<instantiation_tuple> temp_possible_values(1);
						instantiation_tuple& values = temp_possible_values[0];
						if (!init(values, sets.sets[i].arity)) {
							for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
								sets.sets[i].elements[sets.sets[i].elements.length++] = element[k];
							for (auto& element : unifications) free(element);
							for (auto entry : old_unary_conjuncts) free(entry.value);
							return false;
						}
						for (unsigned int k = 0; k < values.length; k++) {
							free(values.values[k]);
							values.values[k].type = instantiation_type::CONSTANT;
							values.values[k].constant = element[k];
						}
						temp_possible_values.length++;

						if (is_provable_without_abduction<false>(conjuncts[l], quantifiers, temp_possible_values)) {
							for (auto& element : temp_possible_values) free(element);
							conjuncts.remove(l--);
						}
					}

					/* `conjuncts` is the collection of sets that `element` is no longer a member of */
					if (conjuncts.length == 0) {
						for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
							sets.sets[i].elements[sets.sets[i].elements.length++] = element[k];
						continue;
					} else if (sets.sets[i].arity == 1) {
						if (!old_unary_conjuncts.ensure_capacity(old_unary_conjuncts.size + 1)
						 || !removed_elements.ensure_capacity(removed_elements.size + 1))
						{
							for (auto& element : unifications) free(element);
							for (auto entry : old_unary_conjuncts) free(entry.value);
							return false;
						}
						unsigned int index = old_unary_conjuncts.index_of(element[0]);
						if (index == old_unary_conjuncts.size) {
							if (!array_init(old_unary_conjuncts.values[index], max((size_t) 8, conjuncts.length))) {
								for (auto& element : unifications) free(element);
								for (auto entry : old_unary_conjuncts) free(entry.value);
								return false;
							}
							old_unary_conjuncts.keys[index] = element[0];
							old_unary_conjuncts.size++;
						} if (!old_unary_conjuncts.values[index].append(conjuncts.data, conjuncts.length)) {
							for (auto& element : unifications) free(element);
							for (auto entry : old_unary_conjuncts) free(entry.value);
							return false;
						}
						removed_elements.keys[removed_elements.size] = element[0];
						removed_elements.values[removed_elements.size++] = i;
					}
					changed = true;
				}
				for (auto& element : unifications) free(element);
			}
		}

		for (const auto& entry : removed_elements)
			sets.sets[entry.value].elements[sets.sets[entry.value].elements.length++] = entry.key;

		array_map<unsigned int, Formula*> old_unary_subsets(max((size_t) 1, old_unary_conjuncts.size));
		for (auto entry : old_unary_conjuncts) {
			Formula* old_subset;
			if (entry.value.length == 1)
				old_subset = entry.value[0];
			else old_subset = Formula::new_and(make_array_view(entry.value.data, entry.value.length));
			if (old_subset == nullptr) {
				for (auto entry : old_unary_subsets) {
					free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
				} for (auto entry : old_unary_conjuncts) free(entry.value);
				return false;
			}
			for (Formula* conjunct : entry.value)
				conjunct->reference_count++;
			old_unary_subsets.keys[old_unary_subsets.size] = entry.key;
			old_unary_subsets.values[old_unary_subsets.size++] = old_subset;
		}
		for (auto entry : old_unary_conjuncts) free(entry.value);

		/* We need to call `move_element_to_superset` in the exact reverse order
		   in which we call `move_element_to_subset` in
		   `check_set_membership_after_addition`, so we sort
		   `old_unary_subsets`. Otherwise, the order depends on the order of
		   `sets.sets` which can change between here and the corresponding
		   earlier call to `check_set_membership_after_addition`. */
		insertion_sort(old_unary_subsets.keys, old_unary_subsets.values, old_unary_subsets.size);
		for (unsigned int i = old_unary_subsets.size; i > 0; i--) {
			if (!move_element_to_superset(old_unary_subsets.keys[i - 1], old_unary_subsets.values[i - 1], std::forward<Args>(visitor)...)) {
				for (auto entry : old_unary_subsets) {
					free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
				}
				return false;
			}
		}
		for (auto entry : old_unary_subsets) {
			free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_new_set_membership(unsigned int set_id, Args&&... visitor)
	{
		Formula* set_formula = sets.sets[set_id].size_axiom->formula->binary.left->binary.right;

		array<Formula*> quantifiers(1 << (core::log2(sets.sets[set_id].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[set_id].arity; j++) {
			quantifiers[quantifiers.length++] = set_formula;
			set_formula = set_formula->quantifier.operand;
		}

		array<instantiation_tuple> possible_values(4);
		instantiation_tuple& values = possible_values[0];
		if (!init(values, sets.sets[set_id].arity))
			return false;
		possible_values.length++;

		if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values)) {
			for (auto& element : possible_values) free(element);
			return true;
		}

		array<tuple> new_elements(4);
		if (!check_set_membership<ResolveInconsistencies>(set_id, possible_values, new_elements, quantifiers)) {
			for (tuple& tup : new_elements) free(tup);
			return false;
		}

		/* add the elements to the sets */
		for (const tuple& element : new_elements) {
			if ((element.length == 1 && !move_element_to_subset<ResolveInconsistencies>(element[0], sets.sets[set_id].set_formula(), std::forward<Args>(visitor)...))
			 || (element.length != 1 && !sets.sets[set_id].add_element(element)))
			{
				for (tuple& tup : new_elements) free(tup);
				return false;
			}
		}
		for (tuple& tup : new_elements) free(tup);
		return true;
	}

	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_proof(Formula* canonicalized, unsigned int& new_constant, Args&&... args)
	{
		Term* predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<DefinitionsAllowed, Contradiction, ResolveInconsistencies>(canonicalized, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::NOT) {
			return make_proof<!Contradiction, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->unary.operand, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
					return NULL;

				Term* constant;
				Proof* exists_not_proof = make_exists_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (exists_not_proof == NULL) return NULL;
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				 && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				{
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				Proof* assumption = ProofCalculus::new_axiom(canonicalized);
				if (assumption == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				assumption->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_universal_elim(assumption, constant), exists_not_proof), assumption);
				free(*assumption); if (assumption->reference_count == 0) free(assumption);
				if (proof == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, std::forward<Args>(args)...);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				free(*constant); if (constant->reference_count == 0) free(constant);
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				free(*exists_not_proof); if (exists_not_proof->reference_count == 0) free(exists_not_proof);
				return proof;
			}

			Formula* operand = canonicalized->quantifier.operand;
			array_map<unsigned int, unsigned int> variable_map(8);
			variable_map.keys[0] = canonicalized->quantifier.variable;
			variable_map.values[0] = 1;
			variable_map.size = 1;
			while (operand->type == FormulaType::FOR_ALL) {
				if (!variable_map.ensure_capacity(variable_map.size + 1))
					return nullptr;
				variable_map.keys[variable_map.size] = operand->quantifier.variable;
				variable_map.values[variable_map.size] = variable_map.size + 1;
				variable_map.size++;
				operand = operand->quantifier.operand;
			}

			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == FormulaType::IF_THEN) {
				Proof* new_axiom;
				Formula* left = canonicalized->quantifier.operand->binary.left;
				Formula* right = canonicalized->quantifier.operand->binary.right;

				Formula* new_left = map_variables(left, variable_map);
				if (new_left == nullptr) return nullptr;

				Formula* new_right = map_variables(right, variable_map);
				if (new_right == nullptr) {
					free(*new_left); free(new_left);
					return nullptr;
				}

				Term const* predicate; Term const* arg1; Term const* arg2;
				bool atomic = is_atomic(*left, predicate, arg1, arg2);
				if (variable_map.size == 1 && atomic && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL
				 && predicate->type == TermType::CONSTANT
				 && predicate->constant == (unsigned int) built_in_predicates::UNKNOWN)
				{
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_proof ERROR: Definitions are not allowed in this context.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					/* this is a definition of a type */
					/* check the right-hand side is a valid definition */
					if (!valid_definition(right, variable)) {
						fprintf(stderr, "theory.make_proof ERROR: This is not a valid type definition.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					new_constant = get_free_concept_id();

					/* TODO: definitions of new concepts should be biconditionals */
					if (!try_init_concept(new_constant)) return NULL;

					new_axiom = get_subset_axiom<ResolveInconsistencies>(new_left, new_right, variable_map.size, std::forward<Args>(args)...);
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				} else {
					/* make sure there are no unknown predicates */
					if (contains_constant(*left, (unsigned int) built_in_predicates::UNKNOWN)
					 || contains_constant(*right, (unsigned int) built_in_predicates::UNKNOWN))
					{
						fprintf(stderr, "theory.make_proof ERROR: Universally-quantified statement has unknown constants.\n");
						free(*new_left); free(new_left);
						free(*new_right); free(new_right);
						return NULL;
					}

					/* this is a formula of form `![x]:(t(x) => f(x))` */
					new_axiom = get_subset_axiom<ResolveInconsistencies>(new_left, new_right, variable_map.size, std::forward<Args>(args)...);
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				}
				return new_axiom;
			}

		} else if (canonicalized->type == FormulaType::AND) {
			if (Contradiction) {
				if (!negated_conjunction_nodes.ensure_capacity(negated_conjunction_nodes.length + 1))
					return NULL;

				array<unsigned int> indices(canonicalized->array.length);
				for (unsigned int i = 0; i < canonicalized->array.length; i++)
					indices[i] = i + 1;
				indices.length = canonicalized->array.length;
				shuffle(indices);

				unsigned int old_index_count = indices.length;
				if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

				Formula* negated_conjunction = Formula::new_not(canonicalized);
				if (negated_conjunction == NULL) return NULL;
				canonicalized->reference_count++;

				for (unsigned int index : indices) {
					unsigned int index_minus_one = index - 1;
					Proof* operand = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index_minus_one], new_constant, std::forward<Args>(args)...);
					if (operand != NULL) {
						/* we found a disproof of a conjunct */
						Proof* axiom = ProofCalculus::new_axiom(canonicalized);
						if (axiom == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
							free_proof(operand, std::forward<Args>(args)...); return NULL;
						}
						axiom->reference_count++;
						Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
								ProofCalculus::new_conjunction_elim(axiom, make_array_view(&index_minus_one, 1)), operand), axiom);
						free(*axiom); if (axiom->reference_count == 0) free(axiom);
						if (proof == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
							free_proof(operand, std::forward<Args>(args)...); return NULL;
						}
						free(*operand); if (operand->reference_count == 0) free(operand);
						proof->reference_count++;
						/* record the formula with this negated conjunction node */
						negated_conjunction_nodes[negated_conjunction_nodes.length++] = { negated_conjunction, proof };
						return proof;
					}

					if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
				}

				/* we couldn't find a disproof of any of the conjuncts */
				finished_constants(canonicalized, old_index_count, std::forward<Args>(args)...);
				free(*negated_conjunction); if (negated_conjunction->reference_count == 0) free(negated_conjunction);
				return NULL;
			}

			Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
			if (operands == NULL) {
				fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
				return NULL;
			}
			for (unsigned int i = 0; i < canonicalized->array.length; i++) {
				if (new_constant == 0)
					operands[i] = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);
				else operands[i] = make_proof<false, false, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);

				if (operands[i] == NULL) {
					for (unsigned int j = i; j > 0; j--)
						/* undo the changes made by the recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
					free(operands); return NULL;
				}
			}
			Proof* conjunction = ProofCalculus::new_conjunction_intro(make_array_view(operands, canonicalized->array.length));
			if (conjunction == NULL) {
				for (unsigned int j = canonicalized->array.length; j > 0; j--)
					/* undo the changes made by the recursive calls to `make_proof` */
					free_proof(operands[j - 1], std::forward<Args>(args)...);
				free(operands); return NULL;
			}
			for (unsigned int j = 0; j < canonicalized->array.length; j++) {
				free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
			}
			free(operands);
			conjunction->reference_count++;
			return conjunction;

		} else if (canonicalized->type == FormulaType::OR) {
			if (Contradiction) {
				Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
				if (operands == NULL) {
					fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
					return NULL;
				}
				for (unsigned int i = 0; i < canonicalized->array.length; i++) {
					if (new_constant == 0)
						operands[i] = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);
					else operands[i] = make_proof<true, false, ResolveInconsistencies>(canonicalized->array.operands[i], new_constant, std::forward<Args>(args)...);

					if (operands[i] == NULL) {
						for (unsigned int j = i; j > 0; j--)
							/* undo the changes made by recursive calls to `make_proof` */
							free_proof(operands[j - 1], std::forward<Args>(args)...);
						free(operands); return NULL;
					} else if (operands[i]->type == ProofType::PROOF_BY_CONTRADICTION) {
						Proof* absurdity = operands[i]->operands[0];
						absurdity->reference_count++;
						free(*operands[i]); if (operands[i]->reference_count == 0) free(operands[i]);
						operands[i] = absurdity;
					} else {
						Proof* absurdity = ProofCalculus::new_negation_elim(
								ProofCalculus::new_axiom(canonicalized->array.operands[i]), operands[i]);
						if (absurdity == NULL) {
							for (unsigned int j = i; j > 0; j--)
								/* undo the changes made by recursive calls to `make_proof` */
								free_proof(operands[j - 1], std::forward<Args>(args)...);
							free(operands); return NULL;
						}
						absurdity->reference_count++;
						free(*operands[i]); if (operands[i]->reference_count == 0) free(operands[i]);
						operands[i] = absurdity;
					}
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
					free(operands); return NULL;
				}
				axiom->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(
						ProofCalculus::new_disjunction_elim(axiom, make_array_view(operands, canonicalized->array.length)), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], std::forward<Args>(args)...);
					free(operands); return NULL;
				}
				proof->reference_count++;
				for (unsigned int j = 0; j < canonicalized->array.length; j++) {
					free(*operands[j]); if (operands[j]->reference_count == 0) free(operands[j]);
				}
				free(operands);
				return proof;
			}

			if (!disjunction_intro_nodes.ensure_capacity(disjunction_intro_nodes.length + 1))
				return NULL;

			array<unsigned int> indices(canonicalized->array.length);
			for (unsigned int i = 0; i < canonicalized->array.length; i++)
				indices[i] = i + 1;
			indices.length = canonicalized->array.length;
			shuffle(indices);

			unsigned int old_index_count = indices.length;
			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

			for (unsigned int index : indices) {
				Proof* operand = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index - 1], new_constant, std::forward<Args>(args)...);
				if (operand != NULL) {
					/* we found a proof of a disjunct */
					array<Formula*> other_disjuncts(max(1u, canonicalized->array.length - 1));
					for (unsigned int i = 0; i < canonicalized->array.length; i++) {
						if (i == index - 1) continue;
						other_disjuncts[other_disjuncts.length++] = canonicalized->array.operands[i];
					}
					Formula* other_disjunction;
					if (other_disjuncts.length == 1) {
						other_disjunction = other_disjuncts[0];
						other_disjuncts[0]->reference_count++;
					} else {
						other_disjunction = Formula::new_or(make_array_view(other_disjuncts.data, other_disjuncts.length));
						if (other_disjunction == NULL) {
							/* undo changes made by the recursive call to `make_proof` */
							free_proof(operand, std::forward<Args>(args)...); return NULL;
						}
						for (Formula* other_disjunct : other_disjuncts)
							other_disjunct->reference_count++;
					}
					Proof* proof = ProofCalculus::new_disjunction_intro(operand, other_disjunction, index - 1);
					free(*other_disjunction); if (other_disjunction->reference_count == 0) free(other_disjunction);
					if (proof == NULL) {
						/* undo changes made by the recursive call to `make_proof` */
						free_proof(operand, std::forward<Args>(args)...); return NULL;
					}
					free(*operand); if (operand->reference_count == 0) free(operand);
					proof->reference_count++;
					/* record the formula with this disjunction intro node */
					disjunction_intro_nodes[disjunction_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					return proof;
				}

				if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
			}

			/* we couldn't find a proof of any of the disjuncts */
			finished_constants(canonicalized, old_index_count, std::forward<Args>(args)...);
			return NULL;

		} else if (canonicalized->type == FormulaType::IF_THEN) {
			if (Contradiction) {
				Proof* left = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, new_constant, std::forward<Args>(args)...);
				if (left == NULL) return NULL;

				Proof* right;
				if (new_constant == 0 && DefinitionsAllowed)
					right = make_proof<true, true, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);
				else right = make_proof<true, false, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);

				if (right == NULL) {
					free_proof(left, std::forward<Args>(args)...); return NULL;
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					free_proof(right, std::forward<Args>(args)...);
					free_proof(left, std::forward<Args>(args)...);
					return NULL;
				}
				axiom->reference_count++;

				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_implication_elim(axiom, left), right), axiom);
				free(*axiom); if (axiom->reference_count == 0) free(axiom);
				if (proof == NULL) {
					free_proof(right, std::forward<Args>(args)...);
					free_proof(left, std::forward<Args>(args)...);
					return NULL;
				}
				free(*left); if (left->reference_count == 0) free(left);
				free(*right); if (right->reference_count == 0) free(right);
				proof->reference_count++;
				return proof;
			}

			if (!implication_intro_nodes.ensure_capacity(implication_intro_nodes.length + 1))
				return NULL;

			array<unsigned int> indices(2);
			indices[0] = 1; indices[1] = 2;
			indices.length = 2;
			if (sample_uniform(2) == 1)
				swap(indices[0], indices[1]);
			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;
			for (unsigned int i = 0; i < 2; i++) {
				if (indices[i] == 1) {
					Proof* left = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, new_constant, std::forward<Args>(args)...);
					if (left == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* axiom = ProofCalculus::new_axiom(canonicalized->binary.left);
					if (axiom == NULL) {
						free_proof(left, std::forward<Args>(args)...);
						return NULL;
					}
					axiom->reference_count++;

					Proof* proof = ProofCalculus::new_implication_intro(
							ProofCalculus::new_falsity_elim(
								ProofCalculus::new_negation_elim(left, axiom),
								canonicalized->binary.right),
							axiom);
					free(*axiom); if (axiom->reference_count == 0) free(axiom);
					if (proof == NULL) {
						free_proof(left, std::forward<Args>(args)...);
						return NULL;
					}
					free(*left); if (left->reference_count == 0) free(left);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				} else {
					Proof* right = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.right, new_constant, std::forward<Args>(args)...);
					if (right == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* proof = ProofCalculus::new_implication_intro(right, ProofCalculus::new_axiom(canonicalized->binary.left));
					if (proof == NULL) {
						free_proof(right, std::forward<Args>(args)...);
						return NULL;
					}
					free(*right); if (right->reference_count == 0) free(right);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				}
			}

			finished_constants(canonicalized, 2, std::forward<Args>(args)...);
			return NULL;

		} else if (canonicalized->type == FormulaType::EXISTS) {
			if (Contradiction) {
				/* get empty set size axiom, forcing inconsistencies to be resolved */
				Formula* set_formula = canonicalized->quantifier.operand;
				Formula* lambda_formula = Formula::new_lambda(1, set_formula);
				if (lambda_formula == NULL) return NULL;
				set_formula->reference_count++;
				unsigned int set_id; bool is_set_new;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(set_formula, 1, 0, set_id, is_set_new, std::forward<Args>(args)...);
				if (set_size_axiom == NULL) {
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					return NULL;
				} else if (is_set_new && !check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					sets.try_free_set(set_formula, set_id, std::forward<Args>(args)...); return NULL;
				}

				Formula* beta_left = Formula::new_not(Formula::new_exists(1, Formula::new_apply(lambda_formula, &Variables<1>::value)));
				Formula* beta_right = Formula::new_not(Formula::new_exists(1, set_formula));
				if (beta_left == NULL || beta_right == NULL) {
					if (beta_left != NULL) { free(*beta_left); if (beta_left->reference_count == 0) free(beta_left); }
					if (beta_right != NULL) { free(*beta_right); if (beta_right->reference_count == 0) free(beta_right); }
					free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
					sets.try_free_set(set_formula, std::forward<Args>(args)...); return NULL;
				}
				Variables<1>::value.reference_count++;
				set_formula->reference_count++;
				Proof* proof = ProofCalculus::new_equality_elim(
						ProofCalculus::new_beta(beta_left, beta_right),
						ProofCalculus::new_equality_elim(ProofCalculus::new_universal_elim(empty_set_axiom, lambda_formula), set_size_axiom, make_repeated_array_view(0u, 1)),
						make_repeated_array_view(0u, 1));
				free(*beta_left); if (beta_left->reference_count == 0) free(beta_left);
				free(*beta_right); if (beta_right->reference_count == 0) free(beta_right);
				free(*lambda_formula); if (lambda_formula->reference_count == 0) free(lambda_formula);
				if (proof == NULL) { sets.try_free_set(set_formula, std::forward<Args>(args)...); return NULL; }
				lambda_formula->reference_count++;
				proof->reference_count++;
				return proof;
			} else {
				Term* variable = Formula::new_variable(canonicalized->quantifier.variable);
				if (variable == NULL) return NULL;

				Term* constant;
				Proof* operand = make_exists_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, constant, new_constant, std::forward<Args>(args)...);
				if (operand == NULL) {
					free(*variable); if (variable->reference_count == 0) free(variable);
					return NULL;
				}

				array<unsigned int> indices(8);
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1)
				 || (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				  && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				 || !compute_indices(*canonicalized->quantifier.operand, *variable, indices))
				{
					free(*variable); if (variable->reference_count == 0) free(variable);
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_proof(operand, std::forward<Args>(args)...);
					return NULL;
				}
				free(*variable); if (variable->reference_count == 0) free(variable);
				Proof* proof = ProofCalculus::new_existential_intro(operand, indices.data, indices.length, constant);
				if (proof == NULL) {
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_proof(operand, std::forward<Args>(args)...);
					return NULL;
				}
				free(*operand); if (operand->reference_count == 0) free(operand);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				free(*constant); if (constant->reference_count == 0) free(constant);
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				return proof;
			}

		} else if (canonicalized->type == FormulaType::EQUALS) {
			Formula* left = canonicalized->binary.left;
			Formula* right = canonicalized->binary.right;
			Term* arg1; Term* arg2;
			bool atomic = is_atomic(*left, predicate, arg1, arg2);
			if (atomic && predicate->type == TermType::CONSTANT
			 && predicate->constant == (unsigned int) built_in_predicates::SIZE
			 && arg1->type == TermType::LAMBDA && arg2 == NULL
			 && right->type == TermType::INTEGER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				unsigned int arity = 1;
				Term* operand = arg1->quantifier.operand;
				while (operand->type == TermType::LAMBDA) {
					operand = operand->quantifier.operand;
					arity++;
				}

				unsigned int set_id; bool is_set_new;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(arg1->quantifier.operand, arity, right->integer, set_id, is_set_new, std::forward<Args>(args)...);
				if (set_size_axiom == nullptr) {
					return nullptr;
				} else if (is_set_new && !check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
					sets.try_free_set(arg1->quantifier.operand, set_id, std::forward<Args>(args)...);
					return nullptr;
				}
				set_size_axiom->reference_count++;
				return set_size_axiom;

			} else if (atomic && predicate->type == TermType::CONSTANT
					&& predicate->constant == (unsigned int) built_in_predicates::SIZE
			 		&& arg1->type == TermType::CONSTANT && arg2 == NULL
			 		&& right->type == TermType::INTEGER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				/* this is a statement on the size of a set */
				Proof* definition = ground_concepts[arg1->constant - new_constant_offset].definitions[0];
#if !defined(NDEBUG)
				if (ground_concepts[arg1->constant - new_constant_offset].definitions.length == 0
				 || definition->formula->binary.right->type != FormulaType::LAMBDA
				 || definition->formula->binary.right->quantifier.operand->type != FormulaType::UNARY_APPLICATION
				 || definition->formula->binary.right->quantifier.operand->binary.left->type != TermType::CONSTANT
				 || definition->formula->binary.right->quantifier.operand->binary.left->constant != arg1->constant
				 || definition->formula->binary.right->quantifier.operand->binary.right->type != TermType::VARIABLE
				 || definition->formula->binary.right->quantifier.operand->binary.right->variable != definition->formula->binary.right->quantifier.variable)
					fprintf(stderr, "theory.make_proof ERROR: Expected a set definition of the form c=λx.c(x).\n");
#endif
				unsigned int set_id; bool is_set_new;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(definition->formula->binary.right->quantifier.operand, 1, right->integer, set_id, is_set_new, std::forward<Args>(args)...);
				Proof* proof = ProofCalculus::new_equality_elim(definition, set_size_axiom, make_repeated_array_view(3u, 1));
				if (proof == NULL)
					return NULL;
				proof->reference_count++;

				if (is_set_new && !check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
					free_proof(proof);
					return nullptr;
				}
				return proof;

			} else if (left == right || *left == *right) {
				Proof* proof = ProofCalculus::new_beta(left, left);
				if (proof == NULL) return NULL;
				proof->reference_count++;
				return proof;

			} else if (left->type == TermType::CONSTANT && right->type == TermType::CONSTANT) {
				if (left->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return NULL;
					} else {
						new_constant = right->constant;
						Proof* proof = ProofCalculus::new_beta(right, right);
						if (proof == NULL) return NULL;
						proof->reference_count++;
						return proof;
					}
				} else if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					new_constant = left->constant;
					Proof* proof = ProofCalculus::new_beta(left, left);
					if (proof == NULL) return NULL;
					proof->reference_count++;
					return proof;
				} else {
					if (Contradiction) {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return NULL;
					} else {
						/* this is impossible */
						return NULL;
					}
				}
			} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT) {
				if (right->type == TermType::CONSTANT && right->constant != (unsigned int) built_in_predicates::UNKNOWN)
					swap(left, right);

				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				} else {
					if (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT
					 && right->binary.right->type == TermType::CONSTANT
					 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
					{
						/* we disallow objects to be arguments of themselves */
						if (right->binary.right->constant == left->constant)
							return NULL;

						/* check if the other object has this function value */
						if (right->binary.right->constant >= new_constant_offset
						 && ground_concepts[right->binary.right->constant - new_constant_offset].function_values.contains((unsigned int) right->binary.left->constant))
							return NULL;
					}

					/* check that this constant could be a set */
					if (right->type == FormulaType::LAMBDA && is_provably_not_a_set(left->constant))
						return NULL;

					unsigned int min_variable = UINT_MAX;
					min_bound_variable(*right, min_variable);
					Formula* new_right;
					if (min_variable != UINT_MAX) {
						new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
					} else {
						new_right = right;
						right->reference_count++;
					}

					/* check if anything else has this definition */
					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						if (ground_concepts[i].types.keys == NULL || i == left->constant - new_constant_offset) continue;
						for (Proof* proof : ground_concepts[i].definitions) {
							if (proof->formula->binary.right == new_right || *proof->formula->binary.right == *new_right) {
								/* we found a different concept with this definition */
								free(*new_right); if (new_right->reference_count == 0) free(new_right);
								return NULL;
							}
						}
					}

					Formula* new_formula = Formula::new_equals(left, new_right);
					if (new_formula == NULL) {
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						return NULL;
					}
					left->reference_count++;

					Proof* new_proof = ProofCalculus::new_axiom(new_formula);
					free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
					if (new_proof == NULL)
						return NULL;
					new_proof->reference_count++;

					Proof* definition = add_definition<ResolveInconsistencies>(new_proof, std::forward<Args>(args)...);
					if (definition != new_proof) {
						free(*new_proof); free(new_proof);
					}
					return definition;
				}
			} else if ((left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT && left->binary.right->type == TermType::CONSTANT)
					|| (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT && right->binary.right->type == TermType::CONSTANT))
			{
				bool swap_order = (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT && right->binary.right->type == TermType::CONSTANT);
				if (swap_order) swap(left, right);

				if (is_provably_a_set(left->binary.right->constant))
					return NULL;

				/* check if anything else has this as a definition */
				for (unsigned int i = 0; i < ground_concept_capacity; i++) {
					if (ground_concepts[i].types.keys == NULL) continue;
					for (Proof* proof : ground_concepts[i].definitions) {
						if (proof->formula->binary.right == left || *proof->formula->binary.right == *left)
							return NULL;
					}
				}

				/* this is a function value definition */
				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				Formula* new_formula = Formula::new_equals(left, new_right);
				if (new_formula == NULL) {
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					return NULL;
				}
				left->reference_count++;

				Proof* new_proof = ProofCalculus::new_axiom(new_formula);
				free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
				if (new_proof == NULL)
					return NULL;
				new_proof->reference_count++;

				Proof* function_value = add_function_value<ResolveInconsistencies>(new_proof, std::forward<Args>(args)...);
				if (function_value != new_proof) {
					free(*new_proof); free(new_proof);
				} if (function_value == nullptr) {
					return nullptr;
				}

				if (swap_order) {
					Proof* swapped_proof = ProofCalculus::new_equality_elim(
							function_value, ProofCalculus::new_beta(right, right), make_repeated_array_view(2u, 1));
					if (swapped_proof == NULL) {
						free_proof(function_value);
						return NULL;
					}
					free(*function_value);
					swapped_proof->reference_count++;
					return swapped_proof;
				} else {
					return function_value;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline Proof* get_subset_axiom(Formula* antecedent, Formula* consequent, unsigned int arity, Args&&... visitor)
	{
		unsigned int antecedent_set, consequent_set;
		bool is_antecedent_new, is_consequent_new;
		Proof* axiom = sets.template get_subset_axiom<ResolveInconsistencies>(
				antecedent, consequent, arity,
				antecedent_set, consequent_set,
				is_antecedent_new, is_consequent_new,
				std::forward<Args>(visitor)...);
		if (axiom == nullptr) return nullptr;

		if ((is_antecedent_new && !check_new_set_membership<ResolveInconsistencies>(antecedent_set, std::forward<Args>(visitor)...))
		 || (is_consequent_new && !check_new_set_membership<ResolveInconsistencies>(consequent_set, std::forward<Args>(visitor)...)))
		{
			sets.free_subset_axiom(antecedent, consequent, std::forward<Args>(visitor)...);
			return nullptr;
		}
		return axiom;
	}

private:
	template<bool ResolveInconsistencies, typename... Args>
	Proof* add_definition(Proof* definition, Args&&... args)
	{
		Formula* constant = definition->formula->binary.left;
		Formula* new_definition = definition->formula->binary.right;

		if (!try_init_concept(constant->constant))
			return NULL;

		/* check if the axiom already exists */
		for (Proof* definition : ground_concepts[constant->constant - new_constant_offset].definitions) {
			if ((definition->formula->binary.left == constant || *definition->formula->binary.left == *constant)
			 && (definition->formula->binary.right == new_definition || *definition->formula->binary.right == *new_definition))
			{
				definition->reference_count++;
				return definition;
			}
		}

		if (new_definition->type == FormulaType::LAMBDA) {
			/* check if this constant defines any other sets, and indicate to the set reasoning module that they are the same set */
			bool contains;
			unsigned int set_size = UINT_MAX;
			unsigned int set_id = sets.set_ids.get(*new_definition->quantifier.operand, contains);
			if (contains) {
				set_size = sets.sets[set_id].set_size;
			} else if (ground_concepts[constant->constant - new_constant_offset].definitions.length != 0) {
				for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
					Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
					if (definition->formula->binary.right->type != FormulaType::LAMBDA)
						continue;
					set_id = sets.set_ids.get(*definition->formula->binary.right->quantifier.operand, contains);
					if (contains) {
						set_size = sets.sets[set_id].set_size;
						break;
					}
				}
			}

			for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					required_set_size set_size_computer(set_size);
					Proof* first_subset_axiom = get_subset_axiom<ResolveInconsistencies>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, set_size_computer, std::forward<Args>(args)...);
					if (first_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Proof* axiom = get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

							axiom = get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
						}
						return NULL;
					}
					first_subset_axiom->reference_count++;

					Proof* second_subset_axiom = get_subset_axiom<ResolveInconsistencies>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
					if (second_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						free(*first_subset_axiom);
						if (first_subset_axiom->reference_count == 1)
							sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Proof* axiom = get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

							axiom = get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
							free(*axiom);
							if (axiom->reference_count == 1)
								sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
						}
						return NULL;
					}
					second_subset_axiom->reference_count++;
				}
			}
		}

		if (!ground_concepts[constant->constant - new_constant_offset].definitions.add(definition)) {
			/* undo the changes we've made so far */
			on_subtract_changes(std::forward<Args>(args)...);
			for (unsigned int j = 0; j < ground_concepts[constant->constant - new_constant_offset].definitions.length; j++) {
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
				Proof* axiom = get_subset_axiom<false>(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == 1)
					sets.free_subset_axiom(new_definition->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

				axiom = get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == 1)
					sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, new_definition->quantifier.operand, std::forward<Args>(args)...);
			}
			return NULL;
		}
		definition->reference_count++;

		if (!check_set_membership_after_addition<ResolveInconsistencies>(definition->formula, std::forward<Args>(args)...)) {
			remove_definition(definition, std::forward<Args>(args)...);
			return nullptr;
		}
		return definition;
	}

	template<typename... Args>
	void remove_definition(Proof* definition, Args&&... args) {
		unsigned int concept_id = definition->formula->binary.left->constant;

		/* remove subset edges from other set definitions for `concept_id` */
		unsigned int index = ground_concepts[concept_id - new_constant_offset].definitions.length;
		for (unsigned int i = ground_concepts[concept_id - new_constant_offset].definitions.length; i > 0; i--) {
			Proof* other_definition = ground_concepts[concept_id - new_constant_offset].definitions[i - 1];
			if (definition->formula->binary.right == other_definition->formula->binary.right) {
				index = i - 1;
				continue;
			}
			if (other_definition->formula->binary.right->type != FormulaType::LAMBDA)
				continue;
			if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
				Proof* axiom = get_subset_axiom<false>(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					sets.free_subset_axiom(other_definition->formula->binary.right->quantifier.operand, definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);

				axiom = get_subset_axiom<false>(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, 1, std::forward<Args>(args)...);
				free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					sets.free_subset_axiom(definition->formula->binary.right->quantifier.operand, other_definition->formula->binary.right->quantifier.operand, std::forward<Args>(args)...);
			}
		}

#if !defined(NDEBUG)
		if (index == ground_concepts[concept_id - new_constant_offset].definitions.length)
			fprintf(stderr, "remove_definition WARNING: Unable to find definition to remove.\n");
#endif

		ground_concepts[concept_id - new_constant_offset].definitions.remove(index);
		try_free_concept_id(concept_id);
		check_set_membership_after_subtraction(definition->formula, std::forward<Args>(args)...);
		free(*definition); if (definition->reference_count == 0) free(definition);
	}

	template<bool ResolveInconsistencies, typename... Args>
	Proof* add_function_value(Proof* function_value_axiom, Args&&... args)
	{
		Formula* constant = function_value_axiom->formula->binary.left->binary.right;
		if (!try_init_concept(constant->constant))
			return NULL;

		/* first check if this function value is already defined */
		array_map<unsigned int, Proof*>& function_values = ground_concepts[constant->constant - new_constant_offset].function_values;
		if (!function_values.ensure_capacity(function_values.size + 1))
			return NULL;
		unsigned int index = function_values.index_of(function_value_axiom->formula->binary.left->binary.left->constant);
		if (index < function_values.size) {
			if ((function_value_axiom->formula->binary.right == function_values.values[index]->formula->binary.right)
			 || (*function_value_axiom->formula->binary.right == *function_values.values[index]->formula->binary.right))
			{
				/* this definition is already an axiom, so return it */
				function_values.values[index]->reference_count++;
				return function_values.values[index];
			} else {
				return NULL;
			}
		}

		function_values.keys[index] = function_value_axiom->formula->binary.left->binary.left->constant;
		function_values.values[index] = function_value_axiom;
		function_values.size++;
		function_value_axiom->reference_count++;

		if (!check_set_membership_after_addition<ResolveInconsistencies>(function_value_axiom->formula, std::forward<Args>(args)...)) {
			remove_function_value(function_value_axiom, std::forward<Args>(args)...);
			return nullptr;
		}
		return function_value_axiom;
	}

	template<typename... Args>
	void remove_function_value(Proof* function_value_axiom, Args&&... args) {
		unsigned int concept_id = function_value_axiom->formula->binary.left->binary.right->constant;

		array_map<unsigned int, Proof*>& function_values = ground_concepts[concept_id - new_constant_offset].function_values;
		unsigned int index = function_values.index_of(function_value_axiom->formula->binary.left->binary.left->constant);
#if !defined(NDEBUG)
		if (index == function_values.size)
			fprintf(stderr, "remove_function_value WARNING: Unable to find axiom to remove.\n");
#endif

		function_values.remove_at(index);
		try_free_concept_id(concept_id);
		check_set_membership_after_subtraction(function_value_axiom->formula, std::forward<Args>(args)...);
		free(*function_value_axiom); if (function_value_axiom->reference_count == 0) free(function_value_axiom);
	}

	/* NOTE: this function finds a constant that proves `quantified`, and not ?[x]:`quantified` */
	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_exists_proof(Formula* quantified, unsigned int variable, Term*& constant, unsigned int& new_constant, Args&&... args)
	{
		Formula* var = Formula::new_variable(variable);
		if (var == NULL) return NULL;

		array<instance> constants(ground_concept_capacity + 1);
		hash_set<int64_t> integers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != NULL) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length++].constant = new_constant_offset + i;

				for (const auto& entry : ground_concepts[i].function_values) {
					Term* constant = entry.value->formula->binary.right;
					if (constant->type == TermType::INTEGER) {
						if (!integers.add(constant->integer)) return NULL;
					} else if (constant->type == TermType::STRING) {
						bool contains = false;
						for (const string* str : strings)
							if (str == &constant->str || *str == constant->str) { contains = true; break; }
						if (!contains && !strings.add(&constant->str)) return NULL;
					}
				}
			}
		}
		constants[constants.length++].type = instance_type::ANY;
		if (!constants.ensure_capacity(constants.length + integers.size + strings.length))
			return NULL;
		for (int64_t integer : integers) {
			constants[constants.length].type = instance_type::INTEGER;
			constants[constants.length++].integer = integer;
		} for (string* str : strings) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length++].str = str;
		}
		shuffle(constants);

		unsigned int original_constant_count = constants.length;
		if (!filter_constants(*this, quantified, variable, constants, std::forward<Args>(args)...)) {
			free(*var); if (var->reference_count == 0) free(var);
			return NULL;
		}

		for (const instance& id : constants) {
			if (id.type == instance_type::ANY || (id.type == instance_type::CONSTANT && ground_concepts[id.constant - new_constant_offset].types.keys == nullptr)) {
				unsigned int constant_id;
				if (id.type == instance_type::ANY)
					constant_id = get_free_concept_id();
				else constant_id = id.constant;

				if (!try_init_concept(constant_id)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}

				constant = Formula::new_constant(constant_id);
				if (constant == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free_concept_id(constant_id); return NULL;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free(*constant); if (constant->reference_count == 0) free(constant);
					free_concept_id(constant_id); return NULL;
				}

				Proof* proof = make_proof<Contradiction, false, ResolveInconsistencies>(substituted, new_constant, std::forward<Args>(args)...);
				free(*substituted); if (substituted->reference_count == 0) free(substituted);
				if (proof != NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return proof;
				}
				free(*constant); if (constant->reference_count == 0) free(constant);
				if (ground_concepts[constant_id - new_constant_offset].types.keys != NULL)
					/* `make_proof` could have already freed the new concept */
					free_concept_id(constant_id);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
			} else {
				if (id.type == instance_type::CONSTANT) {
					constant = Formula::new_constant(id.constant);
				} else if (id.type == instance_type::INTEGER) {
					constant = Formula::new_int(id.integer);
				} else if (id.type == instance_type::STRING) {
					constant = Formula::new_string(*id.str);
				}
				if (constant == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					free(*constant); if (constant->reference_count == 0) free(constant);
					return NULL;
				}

				Proof* proof = make_proof<Contradiction, DefinitionsAllowed, ResolveInconsistencies>(substituted, new_constant, std::forward<Args>(args)...);
				free(*substituted); if (substituted->reference_count == 0) free(substituted);
				if (proof != NULL) {
					free(*var); if (var->reference_count == 0) free(var);
					return proof;
				}
				free(*constant); if (constant->reference_count == 0) free(constant);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					free(*var); if (var->reference_count == 0) free(var);
					return NULL;
				}
			}
		}

		finished_constants(quantified, original_constant_count, std::forward<Args>(args)...);
		free(*var); if (var->reference_count == 0) free(var);
		return NULL;
	}

	template<bool DefinitionsAllowed, bool Negated, bool ResolveInconsistencies, typename... Args>
	Proof* make_atom_proof(
			Term* atom,
			unsigned int& new_constant,
			Args&&... args)
	{
		if (atom->type == TermType::UNARY_APPLICATION && atom->binary.right->type == TermType::CONSTANT)
		{
			/* this is a unary formula */
			Term* arg1 = atom->binary.right;
			if (atom->binary.left->type != TermType::CONSTANT || atom->binary.left->constant != (unsigned int) built_in_predicates::UNKNOWN) {
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					/* make sure that `predicate` could be a set, and that `arg` could be not a set */
					if (!Negated && atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant))
						return NULL;

					/* this is a definition of an object */
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_atom_proof ERROR: Definitions are not allowed in this context.\n");
						return NULL;
					}
					new_constant = get_free_concept_id();

					Formula* formula = Negated ?
							Formula::new_not(Term::new_apply(atom->binary.left, Term::new_constant(new_constant))) :
							Term::new_apply(atom->binary.left, Term::new_constant(new_constant));
					if (formula == NULL) return NULL;
					atom->binary.left->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;

					if (!try_init_concept(new_constant)) {
						free(*new_axiom); free(new_axiom);
						return NULL;
					} if (!add_unary_atom<Negated, ResolveInconsistencies>(Negated ? *formula->unary.operand : *formula, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom);
						free_concept_id(new_constant); return NULL;
					}
					return new_axiom;
				} else {
					/* this is a formula of form `t(c)` */

					/* we do not allow sets to contain themselves */
					if (!Negated && *atom->binary.left == *arg1)
						return NULL;

					/* make sure that `predicate` could be a set, and that `arg` could be not a set */
					if (!Negated && ((atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant)) || is_provably_a_set(arg1->constant)))
						return NULL;

					/* first check if there is already an extensional edge that proves this formula (for this particular instantiation) */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;
						const Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type == FormulaType::NOT) {
							if (Negated) set_formula = set_formula->unary.operand;
							else continue;
						} else if (Negated) continue;
						if (set_formula->type != FormulaType::UNARY_APPLICATION)
							continue;

						Term* expected_constant_eliminator = nullptr;
						if (set_formula->binary.left->type == TermType::VARIABLE) {
							expected_constant_eliminator = atom->binary.left;
						} else if (set_formula->binary.left->type != TermType::CONSTANT || *set_formula->binary.left != *atom->binary.left) {
							continue;
						}

						if (set_formula->binary.right->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg1;
							else if (*expected_constant_eliminator != *arg1)
								continue;
						} else if (*set_formula->binary.right != *arg1) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							for (Proof* child : entry.value->children) {
								if (child->type == ProofType::UNIVERSAL_ELIMINATION
								 && child->operands[1]->type == ProofType::TERM_PARAMETER
								 && *child->operands[1]->term == *expected_constant_eliminator)
								{
									for (Proof* grandchild : child->children) {
										if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
											/* we found a proof of the atomic formula */
											grandchild->reference_count++;
											return grandchild;
										}
									}
								}
							}
						}
					}

					/* there is no extensional edge that proves this atomic formula,
					  so check if the axiom already exists, and if not, create a new axiom */
					Term* lifted_atom = Term::new_apply(atom->binary.left, &Variables<1>::value);
					if (lifted_atom == nullptr) return nullptr;
					atom->binary.left->reference_count++;
					Variables<1>::value.reference_count++;

					bool contains;
					concept<ProofCalculus>& c = ground_concepts[arg1->constant - new_constant_offset];
					Proof* axiom = (Negated ?
							c.negated_types.get(*lifted_atom, contains) :
							c.types.get(*lifted_atom, contains));
					if (contains) {
						axiom->reference_count++;
						free(*lifted_atom); free(lifted_atom);
						return axiom;
					}
					free(*lifted_atom); free(lifted_atom);
					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) return NULL;
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated, ResolveInconsistencies>(*atom, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_atom_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else if (atom->type == TermType::BINARY_APPLICATION
				&& atom->ternary.second->type == TermType::CONSTANT
				&& atom->ternary.third->type == TermType::CONSTANT)
		{
			/* this is a binary formula */
			Term* arg1 = atom->ternary.second;
			Term* arg2 = atom->ternary.third;
			if (atom->ternary.first->type != TermType::CONSTANT || atom->ternary.first->constant != (unsigned int) built_in_predicates::UNKNOWN)
			{
				if (arg1->constant == (unsigned int) built_in_predicates::UNKNOWN
				 || arg2->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
					return NULL;
				} else if (atom->ternary.first->type == TermType::CONSTANT && atom->ternary.first->constant == (unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL) {
					if (arg1->type != TermType::INTEGER && arg2->type != TermType::INTEGER)
						return nullptr;
					if (arg1->integer >= arg2->integer) {
						if (Negated) return nullptr;
						Proof* new_proof = ProofCalculus::new_comparison_introduction(atom);
						if (new_proof == nullptr) return nullptr;
						new_proof->reference_count++;
						return new_proof;
					} else {
						if (!Negated) return nullptr;
						Formula* negated = Formula::new_not(atom);
						if (negated == nullptr) return nullptr;
						atom->reference_count++;
						Proof* new_proof = ProofCalculus::new_comparison_introduction(negated);
						free(*negated); if (negated->reference_count == 0) free(negated);
						if (new_proof == nullptr) return nullptr;
						new_proof->reference_count++;
						return new_proof;
					}
				} else {
					/* this is a formula of form `r(c_1,c_2)` */
					if (atom->ternary.first->type != TermType::CONSTANT) {
						fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
						return nullptr;
					}

					/* first check if there is already an extensional edge that proves this formula (for this particular instantiation) */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axiom == NULL) continue;
						const Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type == FormulaType::NOT) {
							if (Negated) set_formula = set_formula->unary.operand;
							else continue;
						} else if (Negated) continue;
						if (set_formula->type != FormulaType::BINARY_APPLICATION)
							continue;

						Term* expected_constant_eliminator = nullptr;
						if (set_formula->ternary.first->type == TermType::VARIABLE) {
							expected_constant_eliminator = atom->ternary.first;
						} else if (set_formula->ternary.first->type != TermType::CONSTANT || *set_formula->ternary.first != *atom->ternary.first) {
							continue;
						}

						if (set_formula->ternary.second->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg1;
							else if (*expected_constant_eliminator != *arg1)
								continue;
						} else if (*set_formula->ternary.second != *arg1) {
							continue;
						}

						if (set_formula->ternary.third->type == TermType::VARIABLE) {
							if (expected_constant_eliminator == nullptr)
								expected_constant_eliminator = arg2;
							else if (*expected_constant_eliminator != *arg2)
								continue;
						} else if (*set_formula->ternary.third != *arg2) {
							continue;
						}

						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[i].children) {
							for (Proof* child : entry.value->children) {
								if (child->type == ProofType::UNIVERSAL_ELIMINATION
								 && child->operands[1]->type == ProofType::TERM_PARAMETER
								 && *child->operands[1]->term == *expected_constant_eliminator)
								{
									for (Proof* grandchild : child->children) {
										if (grandchild->type == ProofType::IMPLICATION_ELIMINATION) {
											/* we found a proof of the atomic formula */
											grandchild->reference_count++;
											return grandchild;
										}
									}
								}
							}
						}
					}

					/* there is no extensional edge that proves this atomic formula,
					   so check if the axiom already exists, and if not, create a new axiom */
					bool contains;
					concept<ProofCalculus>& c = ground_concepts[arg1->constant - new_constant_offset];
					Proof* axiom = (Negated ?
							c.negated_relations.get(relation(atom->ternary.first->constant, 0, arg2->constant), contains) :
							c.relations.get(relation(atom->ternary.first->constant, 0, arg2->constant), contains));
					if (contains) {
						axiom->reference_count++;
						return axiom;
					}
					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) return NULL;
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					free(*formula); if (formula->reference_count == 0) free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated, ResolveInconsistencies>({atom->ternary.first->constant, arg1->constant, arg2->constant}, new_axiom, std::forward<Args>(args)...)) {
						free(*new_axiom); free(new_axiom); return NULL;
					}
					return new_axiom;
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_atom_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else {
			fprintf(stderr, "theory.make_atom_proof ERROR: Unsupported formula type.\n");
			return NULL;
		}
	}

	template<typename... Args>
	void free_proof(Proof* proof, Args&&... args) {
		theory::changes changes;
		if (!get_theory_changes(*proof, changes, std::forward<Args>(args)...)) return;
		subtract_changes(changes, std::forward<Args>(args)...);
		free(*proof); if (proof->reference_count == 0) free(proof);
	}

	inline bool valid_definition(const Formula* right, unsigned int quantified_variable)
	{
		unsigned int predicate; Term const* arg1; Term const* arg2;
		if (is_atomic(*right, predicate, arg1, arg2)) {
			return predicate != (unsigned int) built_in_predicates::UNKNOWN
				&& arg1->type == TermType::VARIABLE
				&& arg1->variable == quantified_variable
				&& arg2 == NULL;
		} else if (right->type == FormulaType::AND) {
			for (unsigned int i = 0; i < right->array.length; i++)
				if (valid_definition(right->array.operands[i], quantified_variable)) return true;
			return false;
		} else {
			return false;
		}
	}

	template<bool Negated>
	bool make_conjunct_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		unsigned int predicate; Term* arg1; Term* arg2;
		if (is_atomic(*constraint, predicate, arg1, arg2)) {
			bool contains;
			if (arg2 == NULL) {
				/* `constraint` is a unary literal */
				Term* atom = Term::new_apply(Term::new_constant(predicate), &Variables<1>::value);
				if (atom == nullptr) return false;
				Variables<1>::value.reference_count++;
				proof = (Negated ? c.negated_types.get(*atom, contains) : c.types.get(*atom, contains));
				free(*atom); free(atom);
				if (!contains) proof = NULL;
				return true;
			} else {
				/* `constraint` is a binary literal */
				relation r = { predicate,
						(arg1->type == TermType::VARIABLE ? 0 : arg1->constant),
						(arg2->type == TermType::VARIABLE ? 0 : arg2->constant) };
				proof = (Negated ? c.negated_relations.get(r, contains) : c.relations.get(r, contains));
				if (!contains) proof = NULL;
				return true;
			}
		} else if (constraint->type == FormulaType::NOT) {
			return make_conjunct_proof<true>(constraint->unary.operand, c, proof);
		} else {
			return false;
		}
	}

	bool make_universal_elim_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		if (constraint->type == FormulaType::AND) {
			Proof** conjunct_proofs = (Proof**) malloc(sizeof(Proof*) * constraint->array.length);
			if (conjunct_proofs == NULL) {
				fprintf(stderr, "theory.make_universal_elim_proof ERROR: Out of memory.\n");
				return false;
			}
			for (unsigned int i = 0; i < constraint->array.length; i++) {
				Formula* conjunct = constraint->array.operands[i];
				if (!make_conjunct_proof<false>(conjunct, c, conjunct_proofs[i])) return false;
				if (conjunct_proofs[i] == NULL) {
					proof = NULL; free(conjunct_proofs);
					return true;
				}
			}
			proof = ProofCalculus::new_conjunction_intro(make_array_view(conjunct_proofs, constraint->array.length));
			free(conjunct_proofs);
			return proof != NULL;
		} else {
			return make_conjunct_proof<false>(constraint, c, proof);
		}
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool move_element_to_subset(unsigned int element, Formula* lifted_literal, Args&&... visitor)
	{
		Formula* new_set_formula; bool contains;
		unsigned int old_set_id = sets.element_map.get({&element, 1}, contains);
		Formula* old_set_formula = (contains ? sets.sets[old_set_id].set_formula() : NULL);
		if (old_set_formula == NULL) {
			/* there are no other atomic formulas for this concept */
			new_set_formula = lifted_literal;
			new_set_formula->reference_count++;
		} else {
			new_set_formula = Formula::new_and(old_set_formula, lifted_literal);
			old_set_formula->reference_count++;
			lifted_literal->reference_count++;
		}
		if (new_set_formula == NULL)
			return false;
		array_map<unsigned int, unsigned int> variable_map(16);
		variable_map.keys[variable_map.size] = 1;
		variable_map.values[variable_map.size++] = 1;
		Formula* canonicalized = Canonicalizer::canonicalize(*new_set_formula, variable_map);
		free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		if (canonicalized == NULL)
			return false;
		if (!sets.template move_element_to_set<ResolveInconsistencies>({&element, 1}, canonicalized, std::forward<Args>(visitor)...)) {
			free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
			return false;
		}

		/* make sure we can't currently prove the negation of this new axiom
		   via set containment (can we prove that the instance does not belong to any ancestor?) */
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		/*unsigned int new_set_id = sets.set_ids.get(*canonicalized);
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axiom == NULL) continue;
			Formula* set_formula = sets.sets[i].set_formula();
			Formula* negated_set_formula = Formula::new_not(set_formula);
			if (negated_set_formula == NULL) {
				if (old_set_formula == NULL) sets.remove_element({&element, 1});
				else sets.move_element_to_superset({&element, 1}, old_set_formula);
				return false;
			}
			set_formula->reference_count++;
			bool contradiction = is_subset(lifted_literal, negated_set_formula) && sets.sets[i].descendants.contains(new_set_id);
			free(*negated_set_formula); free(negated_set_formula);
			if (contradiction) {
				if (old_set_formula == NULL) sets.remove_element({&element, 1});
				else sets.move_element_to_superset({&element, 1}, old_set_formula);
				return false;
			}
		}*/

		return true;
	}

	bool subtract_conjuncts(
		Formula** first, unsigned int first_count,
		Formula** second, unsigned int second_count,
		array<Formula*>& difference)
	{
		unsigned int i = 0, j = 0;
		while (i < first_count && j < second_count)
		{
			if (*first[i] == *second[j]) {
				i++; j++;
			} else if (*first[i] < *second[j]) {
				if (!difference.add(first[i])) return false;
				i++;
			} else {
				j++;
			}
		}

		while (i < first_count) {
			if (!difference.add(first[i])) return false;
			i++;
		}
		return true;
	}

	template<typename... Args>
	bool move_element_to_superset(unsigned int element, Formula* lifted_literal, Args&&... visitor)
	{
		unsigned int set_id = sets.element_map.get({&element, 1});
		Formula* old_set_formula = sets.sets[set_id].set_formula();

		array_map<unsigned int, unsigned int> variable_map(max(16u, sets.sets[set_id].arity));
		for (unsigned int i = 0; i < sets.sets[set_id].arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized = Canonicalizer::canonicalize(*lifted_literal, variable_map);
		array<Formula*> difference(8);
		if (canonicalized == NULL) {
			return false;
		} else if (old_set_formula->type != FormulaType::AND) {
			if (canonicalized->type == FormulaType::AND) {
				if (!subtract_conjuncts(&old_set_formula, 1, canonicalized->array.operands, canonicalized->array.length, difference)) {
					free(*canonicalized); free(canonicalized); return false;
				}
			} else {
				if (!subtract_conjuncts(&old_set_formula, 1, &canonicalized, 1, difference)) {
					free(*canonicalized); free(canonicalized); return false;
				}
			}
		} else if (canonicalized->type == FormulaType::AND) {
			if (!subtract_conjuncts(old_set_formula->array.operands, old_set_formula->array.length, canonicalized->array.operands, canonicalized->array.length, difference)) {
				free(*canonicalized); free(canonicalized); return false;
			}
		} else {
			if (!subtract_conjuncts(old_set_formula->array.operands, old_set_formula->array.length, &canonicalized, 1, difference)) {
				free(*canonicalized); free(canonicalized); return false;
			}
		}
		free(*canonicalized); free(canonicalized);

		Formula* new_set_formula;
		if (difference.length == 0)
			new_set_formula = NULL;
		else if (difference.length == 1)
			new_set_formula = difference[0];
		else new_set_formula = Formula::new_and(make_array_view(difference.data, difference.length));
		if (difference.length != 0 && new_set_formula == NULL) {
			return false;
		}
		for (Formula* conjunct : difference)
			conjunct->reference_count++;

		bool success;
		if (new_set_formula == NULL) {
			sets.remove_element({&element, 1}, std::forward<Args>(visitor)...);
			success = true;
		} else {
			success = sets.move_element_to_superset({&element, 1}, new_set_formula, std::forward<Args>(visitor)...);
			free(*new_set_formula); if (new_set_formula->reference_count == 0) free(new_set_formula);
		}
		return success;
	}
};

template<typename Formula>
constexpr bool filter_operands(const Formula* formula, const array<unsigned int>& constants) { return true; }

inline unsigned int index_of_any(const array<instance>& constants) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::ANY)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_constant(const array<instance>& constants, unsigned int constant) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::CONSTANT && constants[i].constant == constant)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_integer(const array<instance>& constants, int64_t integer) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::INTEGER && constants[i].integer == integer)
			return i;
	}
	return constants.length;
}

inline unsigned int index_of_string(const array<instance>& constants, string* str) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::STRING && (constants[i].str == str || *constants[i].str == *str))
			return i;
	}
	return constants.length;
}

inline bool contains_any(const array<instance>& constants) {
	return index_of_any(constants) < constants.length;
}

inline bool contains_constant(const array<instance>& constants, unsigned int constant) {
	return index_of_constant(constants, constant) < constants.length;
}

inline bool contains_integer(const array<instance>& constants, int64_t integer) {
	return index_of_integer(constants, integer) < constants.length;
}

inline bool contains_string(const array<instance>& constants, string* str) {
	return index_of_string(constants, str) < constants.length;
}

template<typename ProofCalculus, typename Canonicalizer>
bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Term Term;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	if (formula->type == FormulaType::EQUALS) {
		Term* left = formula->binary.left;
		Term* right = formula->binary.right;
		if (left == right || *left == *right) {
			/* no-op */
		} else if ((left->type == TermType::VARIABLE && left->variable == variable)
				|| (right->type == TermType::VARIABLE && right->variable == variable))
		{
			if (right->type == TermType::VARIABLE && right->variable == variable)
				swap(left, right);
			if (right->type == TermType::CONSTANT) {
				unsigned int constant_id = right->constant;
				if (!contains_constant(constants, constant_id) && (!contains_any(constants) || T.ground_concepts[constant_id - T.new_constant_offset].types.keys != nullptr))
					return false;
				constants[0].type = instance_type::CONSTANT;
				constants[0].constant = right->constant;
				constants.length = 1;
			} else if (right->type == TermType::INTEGER) {
				int64_t integer = right->integer;
				if (!contains_integer(constants, integer) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::INTEGER;
				constants[0].integer = right->integer;
				constants.length = 1;
			} else if (right->type == TermType::STRING) {
				string* str = &right->str;
				if (!contains_string(constants, str) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::STRING;
				constants[0].str = &right->str;
				constants.length = 1;
			}
			/* TODO: handle the case where `right` is an integer or a string */
			else if (right->type == TermType::VARIABLE) {
				/* no-op */
			} else {
				if (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT
				 && right->binary.right->type == TermType::CONSTANT
				 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
				{
					/* we disallow objects to be arguments of themselves */
					unsigned int index = index_of_constant(constants, right->binary.right->constant);
					if (index < constants.length)
						constants.remove(index);

					/* check if the other object has this function value */
					if (right->binary.right->constant >= T.new_constant_offset
					 && T.ground_concepts[right->binary.right->constant - T.new_constant_offset].function_values.contains((unsigned int) right->binary.left->constant))
						return false;
				}

				for (unsigned int i = 0; right->type == FormulaType::LAMBDA && i < constants.length; i++) {
					/* check that this constant could be a set */
					if (constants[i].type == instance_type::ANY) {
						continue;
					} else if (constants[i].type != instance_type::CONSTANT) {
						constants.remove(i--);
					} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
						continue;
					} else if (T.is_provably_not_a_set(constants[i].constant)) {
						constants.remove(i--);
					}
				}
				if (constants.length == 0) return false;

				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				/* check if anything else has this definition */
				for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
					if (T.ground_concepts[i].types.keys == nullptr) continue;
					for (Proof* proof : T.ground_concepts[i].definitions) {
						if (proof->formula->binary.right == new_right || *proof->formula->binary.right == *new_right) {
							/* we found a different concept with this definition */
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							unsigned int constant_id = T.new_constant_offset + i;
							if (!contains_constant(constants, constant_id) && (!contains_any(constants) || T.ground_concepts[i].types.keys != nullptr))
								return false;
							constants[0].type = instance_type::CONSTANT;
							constants[0].constant = constant_id;
							constants.length = 1;
							return true;
						}
					}
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
			}
		} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT
				|| left->type == TermType::VARIABLE || right->type == TermType::VARIABLE)
		{
			if ((left->type != TermType::CONSTANT || left->constant == (unsigned int) built_in_predicates::UNKNOWN) && left->type != TermType::VARIABLE)
				swap(left, right);
			if (left->type == TermType::VARIABLE)
				return true;

			/* check that this constant could be a set */
			if (right->type == FormulaType::LAMBDA && T.is_provably_not_a_set(left->constant))
				return false;

			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY || constants[i].type != instance_type::CONSTANT) {
					continue;
				} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
					continue;
				}

				/* substitute `variable` in the right-hand side with `constants[i]` */
				Term* src_var = Term::new_variable(variable);
				if (src_var == nullptr) return false;
				Term* dst_term = Term::new_constant(constants[i].constant);
				if (dst_term == nullptr) {
					free(*src_var); free(src_var);
					return false;
				}
				Formula* new_right = substitute(right, src_var, dst_term);
				free(*src_var); free(src_var);
				free(*dst_term); if (dst_term->reference_count == 0) free(dst_term);
				if (new_right == nullptr)
					return false;

				if (*left == *new_right) {
					free(*new_right); if (new_right->reference_count == 0) free(new_right);
					return true;
				} else if (new_right->type != TermType::CONSTANT) {
					if (new_right->type == TermType::UNARY_APPLICATION && new_right->binary.left->type == TermType::CONSTANT
					 && new_right->binary.right->type == TermType::CONSTANT
					 && (new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					  || new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					  || new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
					{
						/* we disallow objects to be arguments of themselves */
						if (new_right->binary.right->constant == left->constant) {
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							constants.remove(i--); continue;
						}

						/* check if the other object has this function value */
						if (new_right->binary.right->constant >= T.new_constant_offset
						 && T.ground_concepts[new_right->binary.right->constant - T.new_constant_offset].function_values.contains((unsigned int) new_right->binary.left->constant))
						{
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							constants.remove(i--); continue;
						}
					}

					unsigned int min_variable = UINT_MAX;
					min_bound_variable(*new_right, min_variable);
					Formula* new_new_right;
					if (min_variable != UINT_MAX) {
						new_new_right = shift_bound_variables(new_right, -((int) (min_variable - 1)));
					} else {
						new_new_right = new_right;
						new_right->reference_count++;
					}

					/* check if anything else has this definition */
					bool found_conflicting_definition = false;
					for (unsigned int j = 0; j < T.ground_concept_capacity && !found_conflicting_definition; j++) {
						if (T.ground_concepts[j].types.keys == NULL || j == left->constant - T.new_constant_offset) continue;
						for (Proof* proof : T.ground_concepts[j].definitions) {
							if (proof->formula->binary.right == new_new_right || *proof->formula->binary.right == *new_new_right) {
								/* we found a different concept with this definition */
								found_conflicting_definition = true;
								break;
							}
						}
					}
					free(*new_new_right); if (new_new_right->reference_count == 0) free(new_new_right);
					if (found_conflicting_definition)
						constants.remove(i--);
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
			}
		} else if (left->type == TermType::UNARY_APPLICATION
				&& ((left->binary.left->type == TermType::CONSTANT && left->binary.right->type == TermType::VARIABLE && left->binary.right->variable == variable)
				 || (left->binary.right->type == TermType::CONSTANT && left->binary.left->type == TermType::VARIABLE && left->binary.left->variable == variable)))
		{
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT) {
					constants.remove(i--); continue;
				} else if (T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
					continue;
				}

				/* substitute `variable` in the left-hand side with `constants[i]` */
				Term* src_var = Term::new_variable(variable);
				if (src_var == nullptr) return false;
				Term* dst_term = Term::new_constant(constants[i].constant);
				if (dst_term == nullptr) {
					free(*src_var); free(src_var);
					return false;
				}
				Formula* new_left = substitute(left, src_var, dst_term);
				free(*src_var); free(src_var);
				free(*dst_term); if (dst_term->reference_count == 0) free(dst_term);
				if (new_left == nullptr)
					return false;

				/* check if anything else has this as a definition */
				bool found_conflicting_definition = false;
				for (unsigned int i = 0; i < T.ground_concept_capacity && !found_conflicting_definition; i++) {
					if (T.ground_concepts[i].types.keys == NULL) continue;
					for (Proof* proof : T.ground_concepts[i].definitions) {
						if (*proof->formula->binary.right == *new_left) {
							found_conflicting_definition = true;
							break;
						}
					}
				}
				if (found_conflicting_definition) {
					free(*new_left); if (new_left->reference_count == 0) free(new_left);
					constants.remove(i--); continue;
				}

				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*right, min_variable);
				Formula* new_right;
				if (min_variable != UINT_MAX) {
					new_right = shift_bound_variables(right, -((int) (min_variable - 1)));
				} else {
					new_right = right;
					right->reference_count++;
				}

				/* check if this function value is already defined */
				const array_map<unsigned int, Proof*>& function_values = T.ground_concepts[new_left->binary.right->constant - T.new_constant_offset].function_values;
				unsigned int index = function_values.index_of(new_left->binary.left->constant);
				if (index < function_values.size) {
					if ((new_right != function_values.values[index]->formula->binary.right)
					 && (*new_right != *function_values.values[index]->formula->binary.right))
					{
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						free(*new_left); if (new_left->reference_count == 0) free(new_left);
						constants.remove(i--); continue;
					}
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
			}
		}
	} else if (formula->type == FormulaType::UNARY_APPLICATION) {
		Term* left = formula->binary.left;
		Term* right = formula->binary.right;
		if (left->type == TermType::VARIABLE && left->variable == variable) {
			if (right->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, right->constant);
				if (index < constants.length)
					constants.remove(index);
			}
			/* make sure x could be a set in `x(y)` */
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT || T.is_provably_not_a_set(constants[i].constant)) {
					constants.remove(i);
					i--;
				}
			}
		} if (right->type == TermType::VARIABLE && right->variable == variable) {
			if (left->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, left->constant);
				if (index < constants.length)
					constants.remove(index);
			}
			/* make sure y could not be a set in `x(y)` */
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT || T.is_provably_a_set(constants[i].constant)) {
					constants.remove(i);
					i--;
				}
			}
		}
	} else if (formula->type == FormulaType::AND) {
		for (unsigned int i = 0; i < formula->array.length; i++)
			if (!filter_constants(T, formula->array.operands[i], variable, constants)) return false;
	} else if (formula->type == FormulaType::EXISTS || formula->type == FormulaType::FOR_ALL) {
		return filter_constants(T, formula->quantifier.operand, variable, constants);
	}
	return (constants.length > 0);
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant) {
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index) {
	return true;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count) { }

template<bool Negated, typename Formula>
bool intensional_element_of(
		unsigned int predicate,
		const typename Formula::Term* arg1,
		const typename Formula::Term* arg2,
		const Formula* set_formula,
		typename Formula::Term*& unifying_substitution)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;

	switch (set_formula->type) {
	case FormulaType::UNARY_APPLICATION:
		if (Negated || arg1 == NULL || arg2 != NULL) return false;
		if (set_formula->binary.left->type == TermType::VARIABLE && set_formula->binary.left->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = Term::new_constant(predicate);
				if (unifying_substitution == NULL) return false;
			} else if (unifying_substitution->type != TermType::CONSTANT || unifying_substitution->constant != predicate) {
				return false;
			}
		} else if (set_formula->binary.left->type != TermType::CONSTANT || set_formula->binary.left->constant != predicate) {
			return false;
		}

		if (set_formula->binary.right->type == TermType::VARIABLE && set_formula->binary.right->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg1;
				arg1->reference_count++;
			} else if (*unifying_substitution != *arg1) {
				return false;
			}
		} else if (*set_formula->binary.right != *arg1) {
			return false;
		}
		return true;
	case FormulaType::BINARY_APPLICATION:
		if (Negated || arg1 == NULL || arg2 == NULL) return false;
		if (set_formula->ternary.first->type == TermType::VARIABLE && set_formula->ternary.first->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = Term::new_constant(predicate);
				if (unifying_substitution == NULL) return false;
			} else if (unifying_substitution->type != TermType::CONSTANT || unifying_substitution->constant != predicate) {
				return false;
			}
		} else if (set_formula->ternary.first->type != TermType::CONSTANT || set_formula->ternary.first->constant != predicate) {
			return false;
		}

		if (set_formula->ternary.second->type == TermType::VARIABLE && set_formula->ternary.second->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg1;
				arg1->reference_count++;
			} else if (*unifying_substitution != *arg1) {
				return false;
			}
		} else if (*set_formula->ternary.second != *arg1) {
			return false;
		}

		if (set_formula->ternary.third->type == TermType::VARIABLE && set_formula->ternary.third->variable == 1) {
			if (unifying_substitution == NULL) {
				unifying_substitution = arg2;
				arg2->reference_count++;
			} else if (*unifying_substitution != *arg2) {
				return false;
			}
		} else if (*set_formula->ternary.third != *arg2) {
			return false;
		}
		return true;
	case FormulaType::NOT:
		return intensional_element_of<!Negated>(predicate, arg1, arg2, set_formula, unifying_substitution);
	case FormulaType::AND:
		for (unsigned int i = 0; i < set_formula->array.length; i++)
			if (!intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->array.operands[i], unifying_substitution)) return false;
		return true;
	case FormulaType::OR:
		for (unsigned int i = 0; i < set_formula->array.length; i++) {
			if (intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->array.operands[i], unifying_substitution))
				return true;
			if (unifying_substitution != NULL) {
				free(*unifying_substitution);
				if (unifying_substitution->reference_count == 0)
					free(unifying_substitution);
				unifying_substitution = NULL;
			}
		}
		return false;
	case FormulaType::IF_THEN:
		return intensional_element_of<Negated>(predicate, arg1, arg2, set_formula->binary.left, unifying_substitution)
			&& intensional_element_of<!Negated>(predicate, arg1, arg2, set_formula->binary.right, unifying_substitution);
	case FormulaType::TRUE:
		return !Negated;
	case FormulaType::FALSE:
		return Negated;
	case FormulaType::EQUALS:
	case FormulaType::IFF:
	case FormulaType::FOR_ALL:
	case FormulaType::EXISTS:
	case FormulaType::LAMBDA:
		fprintf(stderr, "intensional_element_of ERROR: Not implemented.\n");
		return false;
	case FormulaType::VARIABLE:
	case FormulaType::CONSTANT:
	case FormulaType::PARAMETER:
	case FormulaType::INTEGER:
		break;
	}
	fprintf(stderr, "intensional_element_of ERROR: Unrecognized formula type.\n");
	return false;
}

template<typename ProofCalculus, typename Canonicalizer>
typename ProofCalculus::Language* make_lifted_conjunction(unsigned int concept_id, const theory<ProofCalculus, Canonicalizer>& T)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;

	const concept<ProofCalculus>& c = T.ground_concepts.get(concept_id);
	array<Formula*> conjuncts(16);
	for (const auto& entry : c.types) {
		Formula* new_conjunct = Formula::new_apply(&entry.key, &Term::template variables<1>::value);
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		entry.key.reference_count++;
		Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.negated_types) {
		Formula* new_conjunct = Formula::new_not(Formula::new_apply(&entry.key, &Term::template variables<1>::value));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		entry.key.reference_count++;
		Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.relations) {
		Formula* new_conjunct = Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg2));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		if (entry.key.arg1 == 0) Term::template variables<1>::value.reference_count++;
		if (entry.key.arg2 == 0) Term::template variables<1>::value.reference_count++;
	} for (const auto& entry : c.negated_relations) {
		Formula* new_conjunct = Formula::new_not(Formula::new_atom(entry.key.predicate,
				(entry.key.arg1 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg1),
				(entry.key.arg2 == 0) ? &Term::template variables<1>::value : Formula::new_constant(entry.key.arg2)));
		if (new_conjunct == NULL || !conjuncts.add(new_conjunct)) { free_formulas(conjuncts); return NULL; }
		if (entry.key.arg1 == 0) Term::template variables<1>::value.reference_count++;
		if (entry.key.arg2 == 0) Term::template variables<1>::value.reference_count++;
	}
	Formula* conjunction = Formula::new_and(conjuncts);
	return conjunction;
}

template<typename Proof>
struct theory_sample {
	Proof** proofs;
	unsigned int count;
	double log_probability;
#if !defined(NDEBUG)
	/* a unique identifier for debugging */
	unsigned int id;
#endif

	static inline unsigned int hash(const theory_sample<Proof>& key) {
		unsigned int hash_value = 0;
		for (unsigned int i = 0; i < key.count; i++)
			hash_value ^= Proof::hash(*key.proofs[i]);
		return hash_value;
	}

	static inline bool is_empty(const theory_sample<Proof>& key) {
		return key.proofs == nullptr;
	}

	static inline void move(const theory_sample<Proof>& src, theory_sample<Proof>& dst) {
		dst.proofs = src.proofs;
		dst.count = src.count;
		dst.log_probability = src.log_probability;
#if !defined(NDEBUG)
		dst.id = src.id;
#endif
	}

	static inline void free(theory_sample<Proof>& sample) { sample.free(); }

private:
	inline void free() {
		for (unsigned int i = 0; i < count; i++) { core::free(*proofs[i]); if (proofs[i]->reference_count == 0) core::free(proofs[i]); }
		core::free(proofs);
	}
};

template<typename Proof>
bool init(theory_sample<Proof>& sample, const hash_set<Proof*>& proofs, double log_probability) {
	sample.log_probability = log_probability;
	sample.proofs = (Proof**) malloc(max((size_t) 1, sizeof(Proof*) * proofs.size));
	if (sample.proofs == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `theory_sample.proofs`.\n");
		return false;
	}

	sample.count = 0;
	for (Proof* proof : proofs) {
		sample.proofs[sample.count] = (Proof*) malloc(sizeof(Proof));
		if (sample.proofs[sample.count] == nullptr
		 || !Proof::clone(*proof, *sample.proofs[sample.count]))
		{
			if (sample.proofs[sample.count] != nullptr) free(sample.proofs[sample.count]);
			for (unsigned int j = 0; j < sample.count; j++)
				free(sample.proofs[j]);
			free(sample.proofs);
			return false;
		}
		sample.count++;
	}

	/* sort the proofs in canonical order */
	if (sample.count > 1)
		sort(sample.proofs, sample.count, pointer_sorter());
	return true;
}

struct null_collector {
	template<typename Proof>
	constexpr inline bool has_prior(const Proof* proof) const {
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept(const hash_set<Proof*>& sample, double proof_prior_diff) const {
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept_with_observation_changes(const hash_set<Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes) const
	{
		return true;
	}
};

template<typename Proof>
inline bool operator == (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr || first.count != second.count)
		return false;
	for (unsigned int i = 0; i < first.count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return false;
	}
	return true;
}

template<typename Proof>
inline bool operator != (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr || first.count != second.count)
		return true;
	for (unsigned int i = 0; i < first.count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return true;
	}
	return false;
}

template<typename Proof>
inline double log_probability(const theory_sample<Proof>& sample) {
	return sample.log_probability;
}

template<typename ProofCalculus, typename Canonicalizer, typename OnProofSampleFunction = no_op>
struct model_evidence_collector
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;

	hash_set<theory_sample<Proof>> samples;
	unsigned int observation_count;
	double current_log_probability;
	Proof* test_proof;
	OnProofSampleFunction on_new_proof_sample;

#if !defined(NDEBUG)
	std::function<double(void)> compute_current_log_probability;
#endif

	template<typename ProofPrior>
	model_evidence_collector(const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior, Proof* test_proof) :
		model_evidence_collector(T, proof_prior, test_proof, *this)
	{ }

	template<typename ProofPrior, typename TheorySampleCollector>
	model_evidence_collector(const theory<ProofCalculus, Canonicalizer>& T,
			ProofPrior& proof_prior, Proof* test_proof, TheorySampleCollector& sample_collector,
			OnProofSampleFunction on_new_proof_sample = no_op()) :
		samples(1024), observation_count(T.observations.size), test_proof(test_proof), on_new_proof_sample(on_new_proof_sample)
	{
#if !defined(NDEBUG)
		if (test_proof->type != ProofType::EXISTENTIAL_INTRODUCTION) {
			fprintf(stderr, "model_evidence_collector ERROR: `test_proof` is not an existential introduction.\n");
			throw std::runtime_error("`test_proof` is not an existential introduction.");
		}
#endif

		/* initialize `current_log_probability` */
		current_log_probability = log_probability(T.observations, proof_prior, sample_collector);
#if !defined(NDEBUG)
		compute_current_log_probability = [&]() { return log_probability(T.observations, proof_prior, sample_collector); };
#endif

		/* add the first sample */
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!init(new_sample, T.observations, current_log_probability))
			throw std::runtime_error("Failed to initialize first theory_sample.");
		on_new_proof_sample(test_proof, current_log_probability);

#if !defined(NDEBUG)
		new_sample.id = samples.size;
#endif
		unsigned int bucket = samples.index_of(new_sample);
		move(new_sample, samples.keys[bucket]);
		samples.size++;
	}

	~model_evidence_collector() { free(); }

	constexpr inline bool has_prior(const Proof* proof) const {
		return (proof != test_proof);
	}

	bool accept(const hash_set<typename ProofCalculus::Proof*>& sample, double proof_prior_diff)
	{
		typedef typename ProofCalculus::Proof Proof;

		current_log_probability += proof_prior_diff;
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!samples.check_size()
		 || !init(new_sample, sample, current_log_probability))
			return false;

		bool contains;
		unsigned int bucket = samples.index_of(new_sample, contains);
		if (contains) {
			/* we've already seen this sample before */
#if !defined(NDEBUG)
			if (fabs(new_sample.log_probability - samples.keys[bucket].log_probability) > 1.0e-9)
				fprintf(stderr, "model_evidence_collector WARNING: Found an"
						" old sample with a different computed log probability."
						" Old log probability: %lf, new log probability: %lf.\n",
						samples.keys[bucket].log_probability, new_sample.log_probability);
#endif
			core::free(new_sample);
			current_log_probability = samples.keys[bucket].log_probability;
			return true;
		}

		/* we've never seen this sample before */
#if !defined(NDEBUG)
		new_sample.id = samples.size;
		double expected_log_probability = compute_current_log_probability();
		if (fabs(expected_log_probability - new_sample.log_probability) > 1.0e-9) {
			fprintf(stderr, "model_evidence_collector WARNING: The computed"
					" log probability of the sample (%lf) differs from the expected log probability (%lf).\n",
					new_sample.log_probability, expected_log_probability);
		}
#endif
		on_new_proof_sample(test_proof, current_log_probability);
		move(new_sample, samples.keys[bucket]);
		samples.size++;
		return true;
	}

	inline bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes)
	{
		if (test_proof != nullptr) {
			for (unsigned int i = 0; i < observation_changes.length; i++) {
				if (observation_changes[i].key == test_proof) {
					test_proof = observation_changes[i].value;
					break;
				}
			}
		}

		return accept(sample, proof_prior_diff);
	}

	inline double total_log_probability() const {
		return logsumexp(samples);
	}

private:
	void free() {
		for (auto& sample : samples)
			core::free(sample);
	}
};

template<typename ProofCalculus, typename Canonicalizer, typename OnProofSampleFunction = no_op>
struct provability_collector
{
	typedef typename ProofCalculus::Proof Proof;

	model_evidence_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction> internal_collector;

	template<typename ProofPrior>
	provability_collector(const theory<ProofCalculus, Canonicalizer>& T,
			ProofPrior& proof_prior, Proof* test_proof,
			OnProofSampleFunction on_new_proof_sample = no_op()) :
		internal_collector(T, proof_prior, test_proof, *this, on_new_proof_sample)
	{ }

	inline bool has_prior(const Proof* proof) const {
		return (proof != internal_collector.test_proof);
	}

	inline bool accept(const hash_set<typename ProofCalculus::Proof*>& sample, double proof_prior_diff)
	{
		return internal_collector.accept(sample, proof_prior_diff);
	}

	bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			double proof_prior_diff, const array<pair<Proof*, Proof*>>& observation_changes)
	{
		return internal_collector.accept_with_observation_changes(sample, proof_prior_diff, observation_changes);
	}

	inline double total_log_probability() const {
		return internal_collector.total_log_probability();
	}
};

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction = no_op>
provability_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction> make_provability_collector(
		const theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofCalculus::Proof* test_proof,
		OnProofSampleFunction on_new_proof_sample = no_op())
{
	return provability_collector<ProofCalculus, Canonicalizer, OnProofSampleFunction>(T, proof_prior, test_proof, on_new_proof_sample);
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, unsigned int num_samples)
{
	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, nullptr);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);
	return collector.total_log_probability();
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability_of_observation(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples)
{
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	Proof* new_proof = T.add_formula(logical_form, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!add(new_proof, proof_prior, proof_axioms)) {
		T.remove_formula(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);

	subtract(collector.test_proof, proof_prior, proof_axioms);
	T.remove_formula(collector.test_proof);
	return collector.total_log_probability();
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability_of_truth(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples)
{
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	Proof* new_proof = T.add_formula(logical_form, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!add(new_proof, proof_prior, proof_axioms)) {
		T.remove_formula(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	provability_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, collector);

	subtract(collector.internal_collector.test_proof, proof_prior, proof_axioms);
	T.remove_formula(collector.internal_collector.test_proof);
	return collector.total_log_probability();
}

template<typename ProofCanonicalizer, typename OnProofSampleFunction>
struct lambda_proof_sample_delegate
{
	OnProofSampleFunction on_new_proof_sample;

	lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) : on_new_proof_sample(on_new_proof_sample) { }

	template<typename Proof>
	inline void operator() (const Proof* test_proof, double log_probability)
	{
		typedef typename Proof::FormulaType Formula;
		typedef typename Formula::Term Term;

		/* get the current value of the term used to introduce the existential quantifier */
		Formula* instantiated_formula = compute_proof_conclusion<built_in_predicates, ProofCanonicalizer>(*test_proof->operands[0]);
		if (test_proof->operands[1]->parameters.length == 0) {
			fprintf(stderr, "lambda_proof_sample_delegate.operator () ERROR: Formula does not depend on the lambda variable.\n");
			core::free(*instantiated_formula); if (instantiated_formula->reference_count == 0) core::free(instantiated_formula);
			return;
		}
		Term* current_term = get_term_at_index(*instantiated_formula, test_proof->operands[1]->parameters[0]);
		on_new_proof_sample(current_term, log_probability);
		free(*instantiated_formula); if (instantiated_formula->reference_count == 0) free(instantiated_formula);
	}
};

template<typename ProofCanonicalizer, typename OnProofSampleFunction>
lambda_proof_sample_delegate<ProofCanonicalizer, OnProofSampleFunction> make_lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) {
	return lambda_proof_sample_delegate<ProofCanonicalizer, OnProofSampleFunction>(on_new_proof_sample);
}

/* TODO: for debugging; delete this */
#include <core/random.h>

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction>
bool log_joint_probability_of_lambda(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		OnProofSampleFunction on_new_proof_sample)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

#if !defined(NDEBUG)
	typedef typename Formula::Type FormulaType;
	if (logical_form->type != FormulaType::LAMBDA)
		fprintf(stderr, "log_joint_probability_of_lambda WARNING: `logical_form` is not a lambda expression.\n");
#endif

	Formula* existential = Formula::new_exists(logical_form->quantifier.variable, logical_form->quantifier.operand);
	if (existential == nullptr)
		return false;
	existential->quantifier.operand->reference_count++;

	unsigned int new_constant;
extern const string_map_scribe* debug_terminal_printer;
T.print_axioms(stderr, *debug_terminal_printer);
	Proof* new_proof = T.add_formula(existential, new_constant);
	free(*existential); if (existential->reference_count == 0) free(existential);
	if (new_proof == nullptr) {
		return false;
	} else if (!proof_axioms.template add<false>(new_proof, proof_prior)) {
		T.remove_formula(new_proof);
		return false;
	}

	auto new_proof_sample_delegate = make_lambda_proof_sample_delegate<typename ProofCalculus::ProofCanonicalizer>(on_new_proof_sample);
	auto collector = make_provability_collector(T, proof_prior, new_proof, new_proof_sample_delegate);
	for (unsigned int t = 0; t < num_samples; t++)
{
fprintf(stderr, "DEBUG: t = %u\n", t);
proof_axioms.check_proof_axioms(T);
proof_axioms.check_universal_eliminations(T, collector);
T.check_concept_axioms();
T.check_disjunction_introductions();
T.are_elements_provable();
T.sets.are_elements_unique();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
T.sets.are_set_sizes_valid();
T.sets.check_set_ids();
if (!T.observations.contains(collector.internal_collector.test_proof))
	fprintf(stderr, "log_joint_probability_of_lambda WARNING: `provability_collector.internal_collector.test_proof` is not an observation in the theory.\n");
if (t == 100)
fprintf(stderr, "DEBUG: BREAKPOINT\n");
T.print_axioms(stderr, *debug_terminal_printer);
		do_mh_step(T, proof_prior, proof_axioms, collector);
if (t % 1000 == 0)
	fprintf(stdout, "(seed = %u)\n", get_seed());
}

	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, proof_prior);
	T.remove_formula(collector.internal_collector.test_proof);
	return true;
}

#endif /* THEORY_H_ */
