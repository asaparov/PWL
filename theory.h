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

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula) { return true; }

template<typename Theory, typename Formula>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable) { return true; }

enum class instance_type {
	ANY,
	CONSTANT,
	NUMBER,
	STRING
};

struct instance {
	instance_type type;
	union {
		unsigned int constant;
		hol_number number;
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
		case instance_type::NUMBER: dst.number = src.number; return;
		case instance_type::STRING: dst.str = src.str; return;
		}
		fprintf(stderr, "instance.move ERROR: Unrecognized instance_type.\n");
		exit(EXIT_FAILURE);
	}
};

inline bool operator == (const instance& first, const instance& second) {
	if (first.type != second.type) return false;
	switch (first.type) {
	case instance_type::ANY: return true;
	case instance_type::CONSTANT:
		return first.constant == second.constant;
	case instance_type::NUMBER:
		return first.number == second.number;
	case instance_type::STRING:
		return (first.str == second.str || *first.str == *second.str);
	}
	fprintf(stderr, "operator == ERROR: Unrecognized `instance_type`.\n");
	exit(EXIT_FAILURE);
}

inline instance instance_any() {
	instance i;
	i.type = instance_type::ANY;
	return i;
}

inline instance instance_constant(unsigned int constant) {
	instance i;
	i.type = instance_type::CONSTANT;
	i.constant = constant;
	return i;
}

inline instance instance_number(hol_number number) {
	instance i;
	i.type = instance_type::NUMBER;
	i.number = number;
	return i;
}

inline instance instance_number(int64_t integer, uint64_t decimal) {
	instance i;
	i.type = instance_type::NUMBER;
	i.number.integer = integer;
	i.number.decimal = decimal;
	return i;
}

inline instance instance_string(string* str) {
	instance i;
	i.type = instance_type::STRING;
	i.str = str;
	return i;
}


/* forward declarations */

struct instantiation;
bool intersect(instantiation&, const instantiation&, const instantiation&);

struct any_instantiation {
	instantiation* excluded;
	uint_fast8_t excluded_count;
};

template<typename T>
struct interval {
	T min, max;
	bool left_inclusive;
	bool right_inclusive;

	inline bool contains(const T& value) {
		return (min < value || (value == min && left_inclusive))
			&& (value < max || (value == max && right_inclusive));
	}
};

template<typename T>
inline bool operator == (const interval<T>& first, const interval<T>& second) {
	return (first.min == second.min) && (first.max == second.max)
		&& (first.left_inclusive == second.left_inclusive)
		&& (first.right_inclusive == second.right_inclusive);
}

template<typename T>
inline bool operator != (const interval<T>& first, const interval<T>& second) {
	return (first.min != second.min) || (first.max != second.max)
		|| (first.left_inclusive != second.left_inclusive)
		|| (first.right_inclusive != second.right_inclusive);
}

template<typename T>
inline int_fast8_t compare(const interval<T>& first, const interval<T>& second) {
	if (first.min < second.min) return 1;
	else if (second.min < first.min) return -1;
	else if (first.max < second.max) return 1;
	else if (second.max < first.max) return -1;

	if (first.left_inclusive) {
		if (!second.left_inclusive)
			return 1;
	} else {
		if (second.left_inclusive)
			return -1;
	}

	if (first.right_inclusive) {
		if (!second.right_inclusive)
			return 1;
	} else {
		if (second.right_inclusive)
			return -1;
	}
	return 0;
}

struct any_number_instantiation {
	interval<hol_number>* included;
	uint_fast8_t included_count;

	const interval<hol_number>& min_interval() const {
		return included[0];
	}

	const interval<hol_number>& max_interval() const {
		return included[included_count - 1];
	}
};

inline bool init(any_number_instantiation& any_number, const any_number_instantiation& src)
{
	any_number.included_count = src.included_count;
	any_number.included = (interval<hol_number>*) malloc(max((size_t) 1, sizeof(interval<hol_number>) * src.included_count));
	if (any_number.included == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	for (uint_fast8_t i = 0; i < any_number.included_count; i++)
		any_number.included[i] = src.included[i];
	return true;
}

enum class instantiation_type {
	ANY,
	ANY_NUMBER,
	CONSTANT,
	NUMBER,
	STRING
};

struct instantiation {
	instantiation_type type;
	union {
		any_instantiation any;
		any_number_instantiation any_number;
		unsigned int constant;
		hol_number number;
		string str;
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
		case instantiation_type::ANY_NUMBER:
			dst.any_number.included = src.any_number.included;
			dst.any_number.included_count = src.any_number.included_count;
			return;
		case instantiation_type::CONSTANT:
			dst.constant = src.constant; return;
		case instantiation_type::NUMBER:
			dst.number = src.number; return;
		case instantiation_type::STRING:
			core::move(src.str, dst.str); return;
		}
		fprintf(stderr, "instantiation.move ERROR: Unrecognized `instantiation_type`.\n");
		exit(EXIT_FAILURE);
	}

	static inline void free(instantiation& value) {
		if (value.type == instantiation_type::ANY) {
			for (uint_fast8_t i = 0; i < value.any.excluded_count; i++)
				core::free(value.any.excluded[i]);
			core::free(value.any.excluded);
		} else if (value.type == instantiation_type::ANY_NUMBER) {
			core::free(value.any_number.included);
		} else if (value.type == instantiation_type::STRING) {
			core::free(value.str);
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
		case instantiation_type::ANY_NUMBER:
			return init(any_number, src.any_number);
		case instantiation_type::CONSTANT:
			constant = src.constant; return true;
		case instantiation_type::NUMBER:
			number = src.number; return true;
		case instantiation_type::STRING:
			return init(str, src.str);
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

inline bool operator == (const any_number_instantiation& first, const any_number_instantiation& second) {
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

inline bool operator < (const any_number_instantiation& first, const any_number_instantiation& second) {
	if (first.included_count < second.included_count) return true;
	else if (first.included_count > second.included_count) return false;
	for (uint_fast8_t i = 0; i < first.included_count; i++) {
		int_fast8_t comparison = compare(first.included[i], second.included[i]);
		if (comparison < 0) return true;
		else if (comparison > 0) return false;
	}
	return false;
}

inline bool operator == (const instantiation& first, const instantiation& second) {
	if (first.type != second.type) return false;
	switch (first.type) {
	case instantiation_type::ANY:
		return first.any == second.any;
	case instantiation_type::ANY_NUMBER:
		return first.any_number == second.any_number;
	case instantiation_type::CONSTANT:
		return first.constant == second.constant;
	case instantiation_type::NUMBER:
		return first.number == second.number;
	case instantiation_type::STRING:
		return first.str == second.str;
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
	case instantiation_type::ANY_NUMBER:
		return first.any_number < second.any_number;
	case instantiation_type::CONSTANT:
		return first.constant < second.constant;
	case instantiation_type::NUMBER:
		return first.number < second.number;
	case instantiation_type::STRING:
		return first.str < second.str;
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

inline bool set_union(any_number_instantiation& dst,
		const any_number_instantiation& first, hol_number second)
{
	dst.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>) * (first.included_count + 1));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_union ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;

	uint_fast8_t i = 1;
	while (i < first.included_count) {
		if (second < first.included[i].min) {
			dst.included[dst.included_count].min = second;
			dst.included[dst.included_count].max = second;
			dst.included[dst.included_count].left_inclusive = true;
			dst.included[dst.included_count++].right_inclusive = true;
			break;
		} else if (second == first.included[i].min) {
			dst.included[dst.included_count].min = second;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = true;
			dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
			i++; break;
		} else if (second < first.included[i].max) {
			dst.included[dst.included_count].min = first.included[i].min;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
			dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
			i++; break;
		} else if (second == first.included[i].max) {
			if (i + 1 < first.included_count && first.included[i + 1].min == second) {
				dst.included[dst.included_count].min = first.included[i].min;
				dst.included[dst.included_count].max = first.included[i + 1].max;
				dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
				dst.included[dst.included_count++].right_inclusive = first.included[i + 1].right_inclusive;
				i += 2; break;
			} else {
				dst.included[dst.included_count].min = first.included[i].min;
				dst.included[dst.included_count].max = first.included[i].max;
				dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
				dst.included[dst.included_count++].right_inclusive = true;
				i++; break;
			}
		} else {
			dst.included[dst.included_count].min = first.included[i].min;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
			dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
			i++;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].min = first.included[i].min;
		dst.included[dst.included_count].max = first.included[i].max;
		dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
		dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
		i++;
	}
	return true;
}

/* returns true iff the left interval `first` is a subset of the left interval `second` */
inline bool left_interval_is_subset(
		hol_number first, bool first_inclusive,
		hol_number second, bool second_inclusive)
{
	return second < first || (first == second && (!first_inclusive || second_inclusive));
}

/* returns true iff the left interval `first` is a proper subset of the left interval `second` (i.e. there exist elements in second but not in first) */
inline bool left_interval_is_proper_subset(
		hol_number first, bool first_inclusive,
		hol_number second, bool second_inclusive)
{
	return second < first || (first == second && !first_inclusive && second_inclusive);
}

/* returns true iff the right interval `first` is a subset of the right interval `second` */
inline bool right_interval_is_subset(
		hol_number first, bool first_inclusive,
		hol_number second, bool second_inclusive)
{
	return first < second || (first == second && (!first_inclusive || second_inclusive));
}

/* returns true iff the right interval `first` is a proper subset of the right interval `second` (i.e. there exist elements in second but not in first) */
inline bool right_interval_is_proper_subset(
		hol_number first, bool first_inclusive,
		hol_number second, bool second_inclusive)
{
	return first < second || (first == second && !first_inclusive && second_inclusive);
}

/* returns true iff there exists a number *not* in the union of the left interval `left` and the right interval `right` */
inline bool has_gap(
		hol_number left, bool left_inclusive,
		hol_number right, bool right_inclusive)
{
	return left < right && (left == right && !left_inclusive && !right_inclusive);
}

/* returns true iff there exists a number in the intersection of the left interval `left` and the right interval `right` */
inline bool has_intersection(
		hol_number left, bool left_inclusive,
		hol_number right, bool right_inclusive)
{
	return right < left && (left == right && left_inclusive && right_inclusive);
}

inline bool set_union(any_number_instantiation& dst,
		const any_number_instantiation& first,
		const any_number_instantiation& second)
{
	hol_number lower, upper;
	bool left_inclusive, right_inclusive;
	uint_fast8_t i, j;
	if (left_interval_is_subset(second.included[0].min, second.included[0].left_inclusive, first.included[0].min, first.included[0].left_inclusive)) {
		lower = first.included[0].min;
		upper = first.included[0].max;
		left_inclusive = first.included[0].left_inclusive;
		right_inclusive = first.included[0].right_inclusive;
		i = 1; j = 0;
	} else {
		lower = second.included[0].min;
		upper = second.included[0].max;
		left_inclusive = second.included[0].left_inclusive;
		right_inclusive = second.included[0].right_inclusive;
		i = 0; j = 1;
	}

	dst.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_union ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	while (i < first.included_count && j < second.included_count) {
		if (!has_gap(first.included[i].min, first.included[i].left_inclusive, upper, right_inclusive)) {
			if (right_interval_is_subset(first.included[i].max, first.included[i].right_inclusive, upper, right_inclusive)) {
				i++;
			} else {
				upper = first.included[i].max;
				right_inclusive = first.included[i++].right_inclusive;
			}
		} else if (!has_gap(second.included[j].min, second.included[j].left_inclusive, upper, right_inclusive)) {
			if (right_interval_is_subset(second.included[j].max, second.included[j].right_inclusive, upper, right_inclusive)) {
				j++;
			} else {
				upper = second.included[j].max;
				right_inclusive = second.included[j++].right_inclusive;
			}
		} else {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].max = upper;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			dst.included[dst.included_count++].right_inclusive = right_inclusive;
			if (left_interval_is_subset(second.included[j].min, second.included[j].left_inclusive, first.included[i].min, first.included[i].left_inclusive)) {
				lower = first.included[i].min;
				upper = first.included[i].max;
				left_inclusive = first.included[i].left_inclusive;
				right_inclusive = first.included[i++].right_inclusive;
			} else {
				lower = second.included[j].min;
				upper = second.included[j].max;
				left_inclusive = second.included[j].left_inclusive;
				right_inclusive = second.included[j++].right_inclusive;
			}
		}
	}

	while (i < first.included_count) {
		if (!has_gap(first.included[i].min, first.included[i].left_inclusive, upper, right_inclusive)) {
			if (right_interval_is_subset(first.included[i].max, first.included[i].right_inclusive, upper, right_inclusive)) {
				i++;
			} else {
				upper = first.included[i].max;
				right_inclusive = first.included[i++].right_inclusive;
			}
		} else {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].max = upper;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			dst.included[dst.included_count++].right_inclusive = right_inclusive;
			lower = first.included[i].min;
			upper = first.included[i].max;
			left_inclusive = first.included[i].left_inclusive;
			right_inclusive = first.included[i++].right_inclusive;
		}
	} while (j < second.included_count) {
		if (!has_gap(second.included[j].min, second.included[j].left_inclusive, upper, right_inclusive)) {
			if (right_interval_is_subset(second.included[j].max, second.included[j].right_inclusive, upper, right_inclusive)) {
				j++;
			} else {
				upper = second.included[j].max;
				right_inclusive = second.included[j++].right_inclusive;
			}
		} else {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].max = upper;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			dst.included[dst.included_count++].right_inclusive = right_inclusive;
			lower = second.included[j].min;
			upper = second.included[j].max;
			left_inclusive = second.included[j].left_inclusive;
			right_inclusive = second.included[j++].right_inclusive;
		}
	}
	return true;
}

inline bool set_subtract(any_number_instantiation& dst,
		const any_number_instantiation& first,
		const any_number_instantiation& second)
{
	dst.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_subtract ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0, j = 0;
	hol_number lower = first.included[0].min;
	bool left_inclusive = first.included[0].left_inclusive;
	while (i < first.included_count && j < second.included_count) {
		if (has_gap(lower, left_inclusive, second.included[j].max, second.included[j].right_inclusive)) {
			j++;
		} else if (left_interval_is_subset(lower, left_inclusive, second.included[j].min, second.included[j].left_inclusive)) {
			if (right_interval_is_proper_subset(second.included[j].max, second.included[j].right_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
				lower = second.included[j].max;
				left_inclusive = !second.included[j++].right_inclusive;
			} else if (i + 1 == first.included_count) {
				i++; break;
			} else {
				lower = first.included[++i].min;
				left_inclusive = first.included[i].left_inclusive;
			}
		} else if (!has_gap(second.included[j].min, second.included[j].left_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].max = second.included[j].min;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			dst.included[dst.included_count++].right_inclusive = !second.included[j].left_inclusive;
			if (right_interval_is_proper_subset(second.included[j].max, second.included[j].right_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
				lower = second.included[j].max;
				left_inclusive = !second.included[j++].right_inclusive;
			} else if (i + 1 == first.included_count) {
				i++; break;
			} else {
				lower = first.included[++i].min;
				left_inclusive = first.included[i].left_inclusive;
			}
		} else {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			dst.included[dst.included_count++].right_inclusive = first.included[i++].right_inclusive;
			if (i == first.included_count) break;
			lower = first.included[i].min;
			left_inclusive = first.included[i].left_inclusive;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].min = lower;
		dst.included[dst.included_count].max = first.included[i].max;
		dst.included[dst.included_count].left_inclusive = left_inclusive;
		dst.included[dst.included_count++].right_inclusive = first.included[i++].right_inclusive;
		if (i < first.included_count) {
			lower = first.included[i].min;
			left_inclusive = first.included[i].left_inclusive;
		}
	}
	return true;
}

inline bool set_subtract(any_number_instantiation& dst,
		const any_number_instantiation& first, hol_number second)
{
	dst.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>) * (first.included_count + 1));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_subtract ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0;
	while (i < first.included_count) {
		if (first.included[i].max < second || (first.included[i].max == second && !first.included[i].right_inclusive)) {
			dst.included[dst.included_count].min = first.included[i].min;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
			dst.included[dst.included_count++].right_inclusive = first.included[i++].right_inclusive;
		} else if (first.included[i].max == second) {
			if (first.included[i].min < first.included[i].max) {
				dst.included[dst.included_count].min = first.included[i].min;
				dst.included[dst.included_count].max = first.included[i].max;
				dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
				dst.included[dst.included_count++].right_inclusive = false;
			}
			i++; break;
		} else if (first.included[i].min < second) {
			dst.included[dst.included_count].min = first.included[i].min;
			dst.included[dst.included_count].max = second;
			dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
			dst.included[dst.included_count++].right_inclusive = false;
			dst.included[dst.included_count].min = second;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = false;
			dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
			i++; break;
		} else if (second == first.included[i].min) {
			dst.included[dst.included_count].min = second;
			dst.included[dst.included_count].max = first.included[i].max;
			dst.included[dst.included_count].left_inclusive = false;
			dst.included[dst.included_count++].right_inclusive = first.included[i].right_inclusive;
			i++; break;
		}
	}

	while (i < first.included_count) {
		dst.included[dst.included_count].min = first.included[i].min;
		dst.included[dst.included_count].max = first.included[i].max;
		dst.included[dst.included_count].left_inclusive = first.included[i].left_inclusive;
		dst.included[dst.included_count++].right_inclusive = first.included[i++].right_inclusive;
	}
	return true;
}

inline bool set_intersect(any_number_instantiation& dst,
		const any_number_instantiation& first,
		const any_number_instantiation& second)
{
	dst.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>) * (first.included_count + second.included_count));
	if (dst.included == nullptr) {
		fprintf(stderr, "set_intersect ERROR: Insufficient memory for `any_number_instantiation.included`.\n");
		return false;
	}
	dst.included_count = 0;
	uint_fast8_t i = 0, j = 0;
	hol_number lower;
	bool left_inclusive;
	if (left_interval_is_subset(second.included[0].min, second.included[0].left_inclusive, first.included[0].min, first.included[0].left_inclusive)) {
		lower = second.included[0].min;
		left_inclusive = second.included[0].left_inclusive;
	} else {
		lower = first.included[0].min;
		left_inclusive = first.included[0].left_inclusive;
	}
	while (i < first.included_count && j < second.included_count) {
		if (!has_intersection(lower, left_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
			i++;
		} else if (!has_intersection(lower, left_inclusive, second.included[j].max, second.included[j].right_inclusive)) {
			j++;
		} else {
			dst.included[dst.included_count].min = lower;
			dst.included[dst.included_count].left_inclusive = left_inclusive;
			if (right_interval_is_subset(second.included[j].max, second.included[j].right_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
				dst.included[dst.included_count].max = second.included[j].max;
				dst.included[dst.included_count].right_inclusive = second.included[j].right_inclusive;
			} else {
				dst.included[dst.included_count].max = first.included[i].max;
				dst.included[dst.included_count].right_inclusive = first.included[i].right_inclusive;
			}
			dst.included_count++;
			i++; j++;

			if (i == first.included_count || j == second.included_count)
				return true;
			if (left_interval_is_subset(second.included[j].min, second.included[j].left_inclusive, first.included[i].min, first.included[i].left_inclusive)) {
				lower = second.included[j].min;
				left_inclusive = second.included[j].left_inclusive;
			} else {
				lower = first.included[i].min;
				left_inclusive = first.included[i].left_inclusive;
			}
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
			if (dst.type == instantiation_type::ANY_NUMBER) {
				instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
				if (tmp.type == instantiation_type::ANY_NUMBER) {
					dummy.type = instantiation_type::ANY_NUMBER;
					if (!set_union(dummy.any_number, dst.any_number, tmp.any_number)) return false;
				} else if (tmp.type == instantiation_type::NUMBER) {
					dummy.type = instantiation_type::ANY_NUMBER;
					if (!set_union(dummy.any_number, dst.any_number, tmp.number)) return false;
				} else {
					fprintf(stderr, "subtract ERROR: Unclosed subtraction of `instantiation` objects.\n");
					return false;
				}
				free(dst); free(tmp);
				move(dummy, dst);
			} else if (tmp.type == instantiation_type::ANY_NUMBER) {
				instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
				if (dst.type == instantiation_type::NUMBER) {
					dummy.type = instantiation_type::ANY_NUMBER;
					if (!set_union(dummy.any_number, tmp.any_number, dst.number)) return false;
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
	} else if (first.type == instantiation_type::ANY_NUMBER) {
		if (second.type == instantiation_type::ANY_NUMBER) {
			if (!set_subtract(dst.any_number, first.any_number, second.any_number)) return false;
		} else if (second.type == instantiation_type::NUMBER) {
			if (!set_subtract(dst.any_number, first.any_number, second.number)) return false;
		} else {
			return init(dst, first);
		}

		if (dst.any_number.included_count == 1 && dst.any_number.included[0].min == dst.any_number.included[0].max) {
			hol_number value = dst.any_number.included[0].min;
			free(dst.any_number.included);
			dst.number = value;
			dst.type = instantiation_type::NUMBER;
		} else {
			dst.type = instantiation_type::ANY_NUMBER;
		}
		return true;
	} else if (second.type == instantiation_type::ANY_NUMBER) {
		if (first.type == instantiation_type::NUMBER) {
			for (uint_fast8_t i = 0; i < second.any_number.included_count; i++)
				if (second.any_number.included[i].contains(first.number)) return false;
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
	case instantiation_type::NUMBER:
		if (first.number == second.number) return false;
		dst.type = instantiation_type::NUMBER;
		dst.number = first.number;
		return true;
	case instantiation_type::STRING:
		if (first.str == second.str) return false;
		dst.type = instantiation_type::STRING;
		return init(dst.str, first.str);
	case instantiation_type::ANY:
	case instantiation_type::ANY_NUMBER:
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
	} else if (second.type == instantiation_type::ANY_NUMBER) {
		if (!subtract(dst, second, first.any.excluded[0]))
			return false;
		if (dst.type == instantiation_type::NUMBER) {
			hol_number value = dst.number;
			dst.any_number.included = (interval<hol_number>*) malloc(sizeof(interval<hol_number>));
			if (dst.any_number.included == nullptr) {
				fprintf(stderr, "intersect_with_any ERROR: Out of memory.\n");
				return false;
			}
			dst.any_number.included[0].min = value;
			dst.any_number.included[0].max = value;
			dst.any_number.included[0].left_inclusive = true;
			dst.any_number.included[0].right_inclusive = true;
			dst.any_number.included_count = 1;
			dst.type = instantiation_type::ANY_NUMBER;
		}

		for (uint_fast8_t i = 1; i < first.any.excluded_count; i++) {
			instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
			if (!subtract(temp, second, first.any.excluded[i])) {
				free(dst);
				return false;
			}

			instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
			if (temp.type == instantiation_type::ANY_NUMBER) {
				if (!set_union(dummy.any_number, dst.any_number, temp.any_number)) {
					free(dst); free(temp);
					return false;
				}
			} else if (temp.type == instantiation_type::NUMBER) {
				if (!set_union(dummy.any_number, dst.any_number, temp.number)) {
					free(dst); free(temp);
					return false;
				}
			}
			dummy.type = instantiation_type::ANY_NUMBER;
			free(dst); free(temp);
			move(dummy, dst);
		}

		if (dst.any_number.included_count == 1 && dst.any_number.included[0].min == dst.any_number.included[0].max) {
			hol_number value = dst.any_number.included[0].max;
			free(dst.any_number.included);
			dst.number = value;
			dst.type = instantiation_type::NUMBER;
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

inline bool intersect_with_any_number(instantiation& dst,
	const instantiation& first, const instantiation& second)
{
	if (second.type == instantiation_type::ANY_NUMBER) {
		if (!set_intersect(dst.any_number, first.any_number, second.any_number))
			return false;
		if (dst.any_number.included_count == 1 && dst.any_number.included[0].min == dst.any_number.included[0].max) {
			hol_number value = dst.any_number.included[0].max;
			free(dst.any_number.included);
			dst.number = value;
			dst.type = instantiation_type::NUMBER;
		} else {
			dst.type = instantiation_type::ANY_NUMBER;
		}
		return true;
	} else if (second.type == instantiation_type::NUMBER) {
		for (uint_fast8_t i = 0; i < first.any_number.included_count; i++)
			if (first.any_number.included[i].contains(second.number)) return init(dst, second);
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
	} else if (first.type == instantiation_type::ANY_NUMBER) {
		return intersect_with_any_number(dst, first, second);
	} else if (second.type == instantiation_type::ANY_NUMBER) {
		return intersect_with_any_number(dst, second, first);
	} else if (first.type != second.type) {
		return false;
	}

	switch (first.type) {
	case instantiation_type::CONSTANT:
		if (first.constant != second.constant) return false;
		dst.type = instantiation_type::CONSTANT;
		dst.constant = first.constant;
		return true;
	case instantiation_type::NUMBER:
		if (first.number != second.number) return false;
		dst.type = instantiation_type::NUMBER;
		dst.number = first.number;
		return true;
	case instantiation_type::STRING:
		if (first.str != second.str) return false;
		dst.type = instantiation_type::STRING;
		return init(dst.str, first.str);
	case instantiation_type::ANY:
	case instantiation_type::ANY_NUMBER:
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
		if (values[index].type == instantiation_type::ANY || values[index].type == instantiation_type::ANY_NUMBER) {
			if (new_value.type == instantiation_type::ANY || new_value.type == instantiation_type::ANY_NUMBER) {
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
				if (new_value.type == instantiation_type::ANY_NUMBER) {
					for (const pair<uint_fast8_t, uint_fast8_t>& pair : ge_indices) {
						if (pair.key == root || pair.value == root) {
							const interval<hol_number>& max_interval = values[pair.key].any_number.max_interval();
							const interval<hol_number>& min_interval = values[pair.value].any_number.min_interval();
							if (max_interval.max == min_interval.min && max_interval.right_inclusive && min_interval.left_inclusive) {
								instantiation& new_value = *((instantiation*) alloca(sizeof(instantiation)));
								new_value.type = instantiation_type::NUMBER;
								new_value.number = max_interval.max;
								change_value(index, new_value);
								return true;
							} else if (has_intersection(min_interval.min, min_interval.left_inclusive, max_interval.max, max_interval.right_inclusive)) {
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

				if (values[index].type == instantiation_type::ANY_NUMBER) {
					for (unsigned int i = 0; i < ge_indices.length; i++) {
						if (ge_indices[i].key == index) {
							instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
							dummy.any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
							dummy.any_number.included[0].min = hol_number::min();
							dummy.any_number.included[0].max = new_value.number;
							dummy.any_number.included[0].left_inclusive = true;
							dummy.any_number.included[0].right_inclusive = true;
							dummy.any_number.included_count = 1;
							instantiation& temp = *((instantiation*) alloca(sizeof(instantiation)));
							if (!intersect(temp, values[ge_indices[i].value], dummy)
							 || !change_value(ge_indices[i].value, temp))
								return false;
							core::free(temp);
							ge_indices.remove(i--);
						} else if (ge_indices[i].value == index) {
							instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
							dummy.any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
							dummy.any_number.included[0].min = new_value.number;
							dummy.any_number.included[0].max = hol_number::max();
							dummy.any_number.included[0].left_inclusive = true;
							dummy.any_number.included[0].right_inclusive = true;
							dummy.any_number.included_count = 1;
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
			core::free(new_value);
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
		} else if (term->type == TermType::NUMBER) {
			dummy[0].type = instantiation_type::NUMBER;
			dummy[0].number = term->number;
			if (!intersect(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::STRING) {
			dummy[0].type = instantiation_type::STRING;
			if (!init(dummy[0].str, term->str)) {
				return false;
			} else if (!intersect(dummy[1], values[index], dummy[0])) {
				core::free(dummy[0]);
				return false;
			}
			bool result = change_value(index, dummy[1]);
			core::free(dummy[0]); core::free(dummy[1]);
			return result;
		} else if (term->type == TermType::VARIABLE) {
			if (!intersect(dummy[1], values[index], values[term->variable - 1]))
				return false;

			if (dummy[1].type != instantiation_type::ANY && dummy[1].type != instantiation_type::ANY_NUMBER) {
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
			if (dummy[1].type == instantiation_type::ANY_NUMBER) {
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
		} else if (term->type == TermType::NUMBER) {
			dummy[0].type = instantiation_type::NUMBER;
			dummy[0].number = term->number;
			if (!subtract(dummy[1], values[index], dummy[0]))
				return false;
			bool result = change_value(index, dummy[1]);
			core::free(dummy[1]); return result;
		} else if (term->type == TermType::STRING) {
			dummy[0].type = instantiation_type::STRING;
			if (!init(dummy[0].str, term->str)) {
				return false;
			} else if (!subtract(dummy[1], values[index], dummy[0])) {
				core::free(dummy[0]);
				return false;
			}
			bool result = change_value(index, dummy[1]);
			core::free(dummy[0]); core::free(dummy[1]);
			return result;
		} else if (term->type == TermType::VARIABLE) {
			if (values[index].type == instantiation_type::ANY || values[index].type == instantiation_type::ANY_NUMBER) {
				if (values[term->variable - 1].type == instantiation_type::ANY || values[term->variable - 1].type == instantiation_type::ANY_NUMBER) {
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
		if (values[first].type != instantiation_type::ANY_NUMBER) {
			dummy[0].type = instantiation_type::ANY_NUMBER;
			dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
			dummy[0].any_number.included[0].min = hol_number::min();
			dummy[0].any_number.included[0].max = hol_number::max();
			dummy[0].any_number.included[0].left_inclusive = true;
			dummy[0].any_number.included[0].right_inclusive = true;
			dummy[0].any_number.included_count = 1;
			if (!intersect(dummy[1], values[first], dummy[0]))
				return false;
			bool result = change_value(first, dummy[1]);
			core::free(dummy[1]);
			if (!result) return false;
		} if (values[second].type != instantiation_type::ANY_NUMBER) {
			dummy[0].type = instantiation_type::ANY_NUMBER;
			dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
			dummy[0].any_number.included[0].min = hol_number::min();
			dummy[0].any_number.included[0].max = hol_number::max();
			dummy[0].any_number.included[0].left_inclusive = true;
			dummy[0].any_number.included[0].right_inclusive = true;
			dummy[0].any_number.included_count = 1;
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

	bool unify_number(uint_fast8_t index) {
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		dummy[0].type = instantiation_type::ANY_NUMBER;
		dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
		dummy[0].any_number.included[0].min = hol_number::min();
		dummy[0].any_number.included[0].max = hol_number::max();
		dummy[0].any_number.included[0].left_inclusive = true;
		dummy[0].any_number.included[0].right_inclusive = true;
		dummy[0].any_number.included_count = 1;
		if (!intersect(dummy[1], values[index], dummy[0]))
			return false;
		bool result = change_value(index, dummy[1]);
		core::free(dummy[1]); return result;
	}

	bool unify_not_number(uint_fast8_t index) {
		if (values[index].type == instantiation_type::NUMBER || values[index].type == instantiation_type::ANY_NUMBER)
			return false;
		/* TODO: we are not currently supporting the case where `values[index].type` is `ANY` */
		return true;
	}

	template<typename Term>
	bool unify_greater_than_or_equal(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (term->type == TermType::NUMBER) {
			dummy[0].type = instantiation_type::ANY_NUMBER;
			dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
			dummy[0].any_number.included[0].min = term->number;
			dummy[0].any_number.included[0].max = hol_number::max();
			dummy[0].any_number.included[0].left_inclusive = true;
			dummy[0].any_number.included[0].right_inclusive = true;
			dummy[0].any_number.included_count = 1;
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
		if (term->type == TermType::NUMBER) {
			dummy[0].type = instantiation_type::ANY_NUMBER;
			dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
			dummy[0].any_number.included[0].min = hol_number::min();
			dummy[0].any_number.included[0].max = term->number;
			dummy[0].any_number.included[0].left_inclusive = true;
			dummy[0].any_number.included[0].right_inclusive = true;
			dummy[0].any_number.included_count = 1;
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

	static inline bool clone(
			const concept<ProofCalculus>& src, concept<ProofCalculus>& dst,
			array_map<const Proof*, Proof*>& proof_map)
	{
		if (!array_map_init(dst.types, src.types.capacity)) {
			return false;
		} else if (!array_map_init(dst.negated_types, src.negated_types.capacity)) {
			core::free(dst.types); return false;
		} else if (!array_map_init(dst.relations, src.relations.capacity)) {
			core::free(dst.types); core::free(dst.negated_types);
			return false;
		} else if (!array_map_init(dst.negated_relations, src.negated_relations.capacity)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); return false;
		} else if (!array_init(dst.definitions, src.definitions.capacity)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); core::free(dst.negated_relations);
			return false;
		} else if (!array_init(dst.existential_intro_nodes, src.existential_intro_nodes.capacity)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); core::free(dst.negated_relations);
			core::free(dst.definitions); return false;
		} else if (!array_map_init(dst.function_values, src.function_values.capacity)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); core::free(dst.negated_relations);
			core::free(dst.definitions); core::free(dst.existential_intro_nodes);
			return false;
		}

		if (!Proof::clone(src.definitions[0], dst.definitions[0], proof_map)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); core::free(dst.negated_relations);
			core::free(dst.definitions); core::free(dst.existential_intro_nodes);
			core::free(dst.function_values); return false;
		}
		dst.definitions[0]->reference_count++;
		dst.definitions.length++;

		for (unsigned int i = 0; i < src.types.size; i++) {
			if (!init(dst.types.keys[dst.types.size], src.types.keys[i])) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.types.values[i], dst.types.values[dst.types.size], proof_map)) {
				core::free(dst.types.keys[dst.types.size]);
				core::free(dst); return false;
			}
			dst.types.size++;
		} for (unsigned int i = 0; i < src.negated_types.size; i++) {
			if (!init(dst.negated_types.keys[dst.negated_types.size], src.negated_types.keys[i])) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.negated_types.values[i], dst.negated_types.values[dst.negated_types.size], proof_map)) {
				core::free(dst.negated_types.keys[dst.negated_types.size]);
				core::free(dst); return false;
			}
			dst.negated_types.size++;
		} for (unsigned int i = 0; i < src.relations.size; i++) {
			if (!Proof::clone(src.relations.values[i], dst.relations.values[dst.relations.size], proof_map)) {
				core::free(dst);
				return false;
			}
			dst.relations.keys[dst.relations.size++] = src.relations.keys[i];
		} for (unsigned int i = 0; i < src.negated_relations.size; i++) {
			if (!Proof::clone(src.negated_relations.values[i], dst.negated_relations.values[dst.negated_relations.size], proof_map)) {
				core::free(dst);
				return false;
			}
			dst.negated_relations.keys[dst.negated_relations.size++] = src.negated_relations.keys[i];
		} for (unsigned int i = 1; i < src.definitions.length; i++) {
			if (!Proof::clone(src.definitions[i], dst.definitions[dst.definitions.length], proof_map)) {
				core::free(dst);
				return false;
			}
			dst.definitions.length++;
		} for (unsigned int i = 0; i < src.existential_intro_nodes.length; i++) {
			/* the `theory` struct own the memory for these proofs */
			unsigned int index = proof_map.index_of(src.existential_intro_nodes[i]);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "concept.clone WARNING: `src.existential_intro_nodes[%u]` does not exist in `proof_map`.\n", i);
#endif
			dst.existential_intro_nodes[dst.existential_intro_nodes.length++] = proof_map.values[index];
		} for (unsigned int i = 0; i < src.function_values.size; i++) {
			if (!Proof::clone(src.function_values.values[i], dst.function_values.values[dst.function_values.size], proof_map)) {
				core::free(dst);
				return false;
			}
			dst.function_values.keys[dst.function_values.size++] = src.function_values.keys[i];
		}
		return true;
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

	inline bool init_first_definition(unsigned int constant, unsigned int arity)
	{
		Formula* constant_expr = Formula::new_constant(constant);
		if (constant_expr == nullptr)
			return false;
		Formula* set_formula = constant_expr;
		constant_expr->reference_count++;
		for (unsigned int i = 0; i < arity; i++) {
			Formula* temp = Formula::new_apply(set_formula, Formula::new_variable(i + 1));
			if (temp == nullptr) {
				core::free(*set_formula); core::free(set_formula);
				core::free(*constant_expr); core::free(constant_expr);
				return false;
			}
			set_formula = temp;
		} for (unsigned int i = arity; i > 0; i--) {
			Formula* temp = Formula::new_lambda(i, set_formula);
			if (temp == nullptr) {
				core::free(*set_formula); core::free(set_formula);
				core::free(*constant_expr); core::free(constant_expr);
				return false;
			}
			set_formula = temp;
		}
		Formula* axiom = Formula::new_equals(constant_expr, set_formula);
		if (constant_expr == nullptr) {
			core::free(*set_formula); core::free(set_formula);
			core::free(*constant_expr); core::free(constant_expr);
			return false;
		}
		definitions[0] = ProofCalculus::new_axiom(axiom);
		if (definitions[0] == nullptr) {
			core::free(*axiom); core::free(axiom);
			return false;
		}
		core::free(*axiom);
		/* we set the initial reference_count to 2 since `theory.free_proof`
		will free definitions when their reference_count is 1 */
		definitions[0]->reference_count += 2;
		return true;
	}
};

template<typename ProofCalculus>
inline bool init(concept<ProofCalculus>& c, unsigned int constant, unsigned int arity) {
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
	} else if (!c.init_first_definition(constant, arity)) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.function_values);
		free(c.existential_intro_nodes); return false;
	}
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

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline bool compute_new_set_size(
		unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		required_set_size& required,
		Args&&... visitor)
{
	if (required.set_size == UINT_MAX) {
		if (!compute_new_set_size(set_id, sets, out, min_set_size, max_set_size, std::forward<Args>(visitor)...))
			return false;
		required.set_size = out;
		return true;
	} else {
		return required.set_size >= min_set_size && required.set_size <= max_set_size
			&& compute_new_set_size(set_id, sets, out, required.set_size, required.set_size, std::forward<Args>(visitor)...);
	}
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const required_set_size& required,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size,
		const array<typename ProofCalculus::Proof*>& freeable_axioms,
		Args&&... visitor)
{
	return compute_new_set_size(set_id, sets, out, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const array<typename ProofCalculus::Proof*>& freeable_axioms,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
}

template<typename Formula>
struct set_changes {
	array<Formula*> old_set_axioms;
	array<Formula*> new_set_axioms;

	set_changes() : old_set_axioms(4), new_set_axioms(4) { }
	~set_changes() { free_helper(); }

	inline bool new_set(Formula* axiom) {
		for (unsigned int i = 0; i < old_set_axioms.length; i++) {
			if (old_set_axioms[i] == axiom || *old_set_axioms[i] == *axiom) {
				core::free(*old_set_axioms[i]);
				if (old_set_axioms[i]->reference_count == 0)
					core::free(old_set_axioms[i]);
				old_set_axioms.remove(i);
				return true;
			}
		}
		if (!new_set_axioms.add(axiom))
			return false;
		axiom->reference_count++;
		return true;
	}

	inline bool old_set(Formula* axiom) {
		for (unsigned int i = 0; i < new_set_axioms.length; i++) {
			if (new_set_axioms[i] == axiom || *new_set_axioms[i] == *axiom) {
				core::free(*new_set_axioms[i]);
				if (new_set_axioms[i]->reference_count == 0)
					core::free(new_set_axioms[i]);
				new_set_axioms.remove(i);
				return true;
			}
		}
		if (!old_set_axioms.add(axiom))
			return false;
		axiom->reference_count++;
		return true;
	}

	inline void clear() {
		free_all(old_set_axioms);
		free_all(new_set_axioms);
		old_set_axioms.clear();
		new_set_axioms.clear();
	}

	static inline void free(set_changes<Formula>& sets) {
		sets.free_helper();
		core::free(sets.new_set_axioms);
		core::free(sets.old_set_axioms);
	}

private:
	inline void free_helper() {
		free_all(old_set_axioms);
		free_all(new_set_axioms);
	}
};

template<typename Formula>
inline bool init(set_changes<Formula>& sets) {
	if (!array_init(sets.old_set_axioms, 4)) {
		return false;
	} else if (!array_init(sets.new_set_axioms, 4)) {
		core::free(sets.old_set_axioms);
		return false;
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline bool compute_new_set_size(
		unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		set_changes<typename ProofCalculus::Language>& set_diff,
		Args&&... visitor)
{
	return compute_new_set_size(set_id, sets, out, min_set_size, max_set_size, std::forward<Args>(visitor)...)
		&& set_diff.new_set(sets.sets[set_id].size_axioms[0]->formula);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		set_changes<typename ProofCalculus::Language>& set_diff,
		Args&&... visitor)
{
	on_free_set(set_id, sets, std::forward<Args>(visitor)...);
	set_diff.old_set(sets.sets[set_id].size_axioms[0]->formula);
}

template<typename Formula, typename... Args>
inline void on_subtract_changes(const set_changes<Formula>& set_diff, Args&&... visitor) {
	on_subtract_changes(std::forward<Args>(visitor)...);
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
	set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer> sets;

	array<pair<Formula*, Proof*>> disjunction_intro_nodes;
	array<pair<Formula*, Proof*>> negated_conjunction_nodes;
	array<pair<Formula*, Proof*>> implication_intro_nodes;
	array<pair<Formula*, Proof*>> existential_intro_nodes;

	Proof* empty_set_axiom;
	Proof* maximality_axiom;
	Proof* function_axiom;

	array<unsigned int> built_in_sets;

	Term* NAME_ATOM;

	theory(unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), atoms(64), relations(64),
			ground_concept_capacity(64), ground_axiom_count(0), observations(64),
			disjunction_intro_nodes(16), negated_conjunction_nodes(16),
			implication_intro_nodes(16), existential_intro_nodes(16), built_in_sets(4)
	{
		ground_concepts = (concept<ProofCalculus>*) malloc(sizeof(concept<ProofCalculus>) * ground_concept_capacity);
		if (ground_concepts == NULL) {
			fprintf(stderr, "theory ERROR: Insufficient memory for `ground_concepts`.\n");
			exit(EXIT_FAILURE);
		}
		for (unsigned int i = 0; i < ground_concept_capacity; i++)
			ground_concepts[i].types.keys = NULL; /* this is used to indicate that this concept is uninitialized */

		Formula* empty_set_formula = Formula::new_for_all(1, Formula::new_equals(
			Formula::new_equals(Formula::new_atom((unsigned int) built_in_predicates::SIZE, &Variables<1>::value), Formula::new_number(0, 0)),
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

		NAME_ATOM = Formula::new_apply(&Constants<(unsigned int) built_in_predicates::NAME>::value, &Variables<1>::value);
		if (NAME_ATOM == nullptr) exit(EXIT_FAILURE);
		Constants<(unsigned int) built_in_predicates::NAME>::value.reference_count++;
		Variables<1>::value.reference_count++;

		/* add definition of maximality */
		Formula* left = Formula::new_exists(4, Formula::new_and(
				Formula::new_apply(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::GREATEST>::value, &Variables<2>::value), &Variables<4>::value),
				Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG1>::value, &Variables<4>::value), &Variables<3>::value),
				Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG2>::value, &Variables<4>::value), &Variables<1>::value)
			));
		if (left == nullptr) exit(EXIT_FAILURE);
		Constants<(unsigned int) built_in_predicates::GREATEST>::value.reference_count++;
		Constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count++;
		Constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count++;
		Variables<1>::value.reference_count++;
		Variables<2>::value.reference_count++;
		Variables<3>::value.reference_count++;
		Variables<4>::value.reference_count += 3;
		Formula* right = Formula::new_and(
				Formula::new_apply(&Variables<3>::value, &Variables<1>::value),
				Formula::new_for_all(4, Formula::new_if_then(
					Formula::new_apply(Formula::new_apply(&Variables<2>::value, &Variables<1>::value), &Variables<4>::value),
					Formula::new_for_all(5, Formula::new_if_then(
						Formula::new_apply(&Variables<3>::value, &Variables<5>::value),
						Formula::new_for_all(6, Formula::new_if_then(
							Formula::new_apply(Formula::new_apply(&Variables<2>::value, &Variables<5>::value), &Variables<6>::value),
							Formula::new_and(
								Formula::new_if_then(
									Formula::new_apply(&Constants<(unsigned int) built_in_predicates::NUMBER>::value, &Variables<4>::value),
									Formula::new_apply(&Constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, &Variables<4>::value, &Variables<6>::value)
								),
								Formula::new_for_all(7, Formula::new_if_then(
									Formula::new_and(
										Formula::new_apply(&Constants<(unsigned int) built_in_predicates::MEASURE>::value, &Variables<4>::value),
										Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG1>::value, &Variables<4>::value), &Variables<7>::value)
									),
									Formula::new_for_all(8, Formula::new_if_then(
										Formula::new_and(
											Formula::new_apply(&Constants<(unsigned int) built_in_predicates::MEASURE>::value, &Variables<6>::value),
											Formula::new_exists(9, Formula::new_and(
												Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG2>::value, &Variables<4>::value), &Variables<9>::value),
												Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG2>::value, &Variables<6>::value), &Variables<9>::value)
											)),
											Formula::new_equals(Formula::new_apply(&Constants<(unsigned int) built_in_predicates::ARG1>::value, &Variables<6>::value), &Variables<8>::value)
										),
										Formula::new_apply(&Constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, &Variables<7>::value, &Variables<8>::value)
									))
								))
							)
						))
					))
				)));
		if (right == nullptr) exit(EXIT_FAILURE);
		Constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count += 2;
		Constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count += 2;
		Constants<(unsigned int) built_in_predicates::NUMBER>::value.reference_count++;
		Constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value.reference_count += 2;
		Constants<(unsigned int) built_in_predicates::MEASURE>::value.reference_count += 2;
		Variables<1>::value.reference_count += 2;
		Variables<2>::value.reference_count += 2;
		Variables<3>::value.reference_count += 2;
		Variables<4>::value.reference_count += 6;
		Variables<5>::value.reference_count += 2;
		Variables<6>::value.reference_count += 5;
		Variables<7>::value.reference_count += 2;
		Variables<8>::value.reference_count += 2;
		Variables<9>::value.reference_count += 2;
		bool is_antecedent_new, is_consequent_new;
		unsigned int antecedent_set, consequent_set;
		maximality_axiom = get_subset_axiom<true>(left, right, 3, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new);
		core::free(*left); if (left->reference_count == 0) core::free(left);
		core::free(*right); if (right->reference_count == 0) core::free(right);
		if (maximality_axiom == NULL || !built_in_sets.add(antecedent_set) || !built_in_sets.add(consequent_set))
			exit(EXIT_FAILURE);
		maximality_axiom->reference_count++;

		/* add definition of a function */
		/* TODO: generalize this to all predicates that are functions */
		left = Formula::new_and(
				Formula::new_apply(Formula::new_constant((unsigned int) built_in_predicates::AREA), &Variables<1>::value),
				Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG1), &Variables<1>::value), &Variables<2>::value));
		if (left == nullptr) exit(EXIT_FAILURE);
		Variables<1>::value.reference_count += 2;
		Variables<2>::value.reference_count++;
		right = Formula::new_for_all(3, Formula::new_if_then(
					Formula::new_and(
						Formula::new_apply(Formula::new_constant((unsigned int) built_in_predicates::AREA), &Variables<3>::value),
						Formula::new_equals(Term::new_apply(Term::new_constant((unsigned int) built_in_predicates::ARG1), &Variables<3>::value), &Variables<2>::value)
					),
					Formula::new_equals(&Variables<1>::value, &Variables<3>::value)
				));
		if (right == nullptr) exit(EXIT_FAILURE);
		Variables<1>::value.reference_count++;
		Variables<2>::value.reference_count++;
		Variables<3>::value.reference_count += 3;
		function_axiom = get_subset_axiom<true>(left, right, 2, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new);
		core::free(*left); if (left->reference_count == 0) core::free(left);
		core::free(*right); if (right->reference_count == 0) core::free(right);
		if (function_axiom == NULL || !built_in_sets.add(antecedent_set) || !built_in_sets.add(consequent_set))
			exit(EXIT_FAILURE);
		function_axiom->reference_count++;
	}

	~theory() { free_helper(); }

	inline void free_helper() {
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
					unsigned int arity = 1;
					Formula* set_formula = definition->formula->binary.right->quantifier.operand;
					while (set_formula->type == hol_term_type::LAMBDA) {
						set_formula = set_formula->quantifier.operand;
						arity++;
					}

					for (unsigned int j = i + 1; j < c.definitions.length; j++) {
						Proof* other_definition = c.definitions[j];
						if (other_definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
						Formula* other_set_formula = other_definition->formula->binary.right->quantifier.operand;
						while (other_set_formula->type == hol_term_type::LAMBDA)
							other_set_formula = other_set_formula->quantifier.operand;

						bool dummy; unsigned int antecedent_set, consequent_set;
						Proof* axiom = sets.template get_subset_axiom<false, false>(
								set_formula, other_set_formula, arity,
								antecedent_set, consequent_set, dummy, dummy);
						core::free(*axiom);
						if (axiom->reference_count == 1)
							sets.template free_subset_axiom<true>(set_formula, other_set_formula, arity);

						axiom = sets.template get_subset_axiom<false, false>(
								other_set_formula, set_formula, arity,
								consequent_set, antecedent_set, dummy, dummy);
						core::free(*axiom);
						if (axiom->reference_count == 1)
							sets.template free_subset_axiom<true>(other_set_formula, set_formula, arity);
					}
				}
			}

			core::free(c);
		}
		core::free(ground_concepts);

		core::free(*NAME_ATOM);
		if (NAME_ATOM->reference_count == 0)
			core::free(NAME_ATOM);
		core::free(*empty_set_axiom);
		if (empty_set_axiom->reference_count == 0)
			core::free(empty_set_axiom);
		core::free(*maximality_axiom);
		sets.free_subset_axiom(maximality_axiom);
		core::free(*function_axiom);
		sets.free_subset_axiom(function_axiom);
	}

	static inline void free(theory<ProofCalculus, Canonicalizer>& T) {
		T.free_helper();
		core::free(T.atoms);
		core::free(T.relations);
		core::free(T.observations);
		core::free(T.sets);
		core::free(T.disjunction_intro_nodes);
		core::free(T.negated_conjunction_nodes);
		core::free(T.implication_intro_nodes);
		core::free(T.existential_intro_nodes);
	}

	unsigned int get_free_concept_id(unsigned int start = 0) {
		for (unsigned int i = start; i < ground_concept_capacity; i++)
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
		else return ::init(c, id, 1);
	}

	inline bool try_init_concept(unsigned int id, unsigned int arity) {
		if (id < new_constant_offset) return true;
		concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.keys != NULL) {
			Formula* set_formula = c.definitions[0]->formula->binary.right;
			unsigned int prev_arity = 0;
			while (set_formula->type == FormulaType::LAMBDA) {
				set_formula = set_formula->quantifier.operand;
				prev_arity++;
			}
			if (prev_arity == arity)
				return true;

			/* try to change the arity of this set */
			Proof* old_definition = c.definitions[0];
			if (old_definition->reference_count != 2 || sets.set_ids.table.contains(*set_formula))
				return false;
			if (!c.init_first_definition(id, arity))
				return false;
			core::free(*old_definition);
			if (old_definition->reference_count == 1)
				core::free(old_definition);
			return true;
		}
		else return ::init(c, id, arity);
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

	template<typename Stream, typename... Printer>
	bool print_existential_introduction(const Proof* proof, Stream& out, Printer&&... printer) const {
		if (proof->type == ProofType::EXISTENTIAL_INTRODUCTION) {
			/* find the formula corresponding to this node */
			unsigned int index;
			for (index = 0; index < existential_intro_nodes.length; index++)
				if (existential_intro_nodes[index].value == proof) break;
			if (index == existential_intro_nodes.length) {
				fprintf(stderr, "print_existential_introduction WARNING: Found existential introduction in proof that is not an element of `theory.existential_intro_nodes`.\n");
				if (!print("  Existential introduction instantiated with ", out)
				 || !print(*proof->operands[2]->term, out, std::forward<Printer>(printer)...)
				 || !print(".\n", out))
					return false;
			} else {
				if (!print("  Existential introduction at index ", out)
				 || !print(index, out) || !print(": ", out)
				 || !print(*existential_intro_nodes[index].key, out, std::forward<Printer>(printer)...)
				 || !print(" instantiated with ", out)
				 || !print(*proof->operands[2]->term, out, std::forward<Printer>(printer)...)
				 || !print(".\n", out))
					return false;
			}
		}

		unsigned int operand_count;
		const Proof* const* operands;
		proof->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			if (!print_existential_introduction(operands[i], out, std::forward<Printer>(printer)...))
				return false;
		}
		return true;
	}

	template<typename Stream, typename... Printer>
	bool print_existential_introductions(Stream& out, Printer&&... printer) const {
		for (const Proof* observation : observations) {
			unsigned int index;
			for (index = 0; index < existential_intro_nodes.length; index++)
				if (existential_intro_nodes[index].value == observation) break;
			if (observation->type == ProofType::EXISTENTIAL_ELIMINATION && index == existential_intro_nodes.length) {
				fprintf(stderr, "print_existential_introduction WARNING: Found existential introduction in proof that is not an element of `theory.existential_intro_nodes`.\n");
				if (fprintf(out, "Proof at address 0x%lx:\n", (size_t) observation) < 0)
					return false;
			} else {
				if (!print("Proof of ", out)
				 || !print(*existential_intro_nodes[index].key, out, std::forward<Printer>(printer)...)
				 || fprintf(out, " at address 0x%lx:\n", (size_t) observation) < 0)
					return false;
			}

			if (!print_existential_introduction(observation, out, std::forward<Printer>(printer)...))
				return false;
		}
		return true;
	}

	inline bool get_concept_names(unsigned int constant, array<Term*>& name_terms) const {
		if (constant < new_constant_offset)
			return true;
		for (Proof* definition : ground_concepts[constant - new_constant_offset].definitions) {
			if (definition->formula->binary.right->type != TermType::UNARY_APPLICATION
			 || definition->formula->binary.right->binary.left->type != TermType::CONSTANT
			 || definition->formula->binary.right->binary.left->constant != (unsigned int) built_in_predicates::ARG1
			 || definition->formula->binary.right->binary.right->type != TermType::CONSTANT)
				continue;

			unsigned int event = definition->formula->binary.right->binary.right->constant;
			if (!ground_concepts[event - new_constant_offset].types.contains(*NAME_ATOM))
				continue;

			bool contains;
			Proof* function_value_axiom = ground_concepts[event - new_constant_offset].function_values.get((unsigned int) built_in_predicates::ARG2, contains);
			if (!contains || function_value_axiom->formula->binary.right->type != TermType::STRING)
				continue;
			if (!name_terms.add(function_value_axiom->formula->binary.right))
				return false;
		}
		return true;
	}

	inline bool get_extra_axioms(array<Formula*>& extra_axioms) const
	{
		for (unsigned int i = 2; i < sets.set_count + 1; i++) {
			if (built_in_sets.contains(i) || sets.sets[i].size_axioms.data == nullptr)
				continue;
			if (!extra_axioms.contains(sets.sets[i].size_axioms[0]->formula)
			 && !extra_axioms.add(sets.sets[i].size_axioms[0]->formula))
				return false;
		}
		return true;
	}

	static inline bool clone(const theory<ProofCalculus, Canonicalizer>& src, theory<ProofCalculus, Canonicalizer>& dst)
	{
		dst.new_constant_offset = src.new_constant_offset;
		if (!hash_map_init(dst.atoms, src.atoms.table.capacity)) {
			return false;
		} else if (!hash_map_init(dst.relations, src.relations.table.capacity)) {
			core::free(dst.atoms);
			return false;
		}
		dst.ground_concept_capacity = src.ground_concept_capacity;
		dst.ground_concepts = (concept<ProofCalculus>*) malloc(sizeof(concept<ProofCalculus>) * dst.ground_concept_capacity);
		if (dst.ground_concepts == nullptr) {
			fprintf(stderr, "theory.clone ERROR: Insufficient memory for `ground_concepts`.\n");
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		}
		for (unsigned int i = 0; i < dst.ground_concept_capacity; i++)
			dst.ground_concepts[i].types.keys = nullptr; /* this is used to indicate that this concept is uninitialized */
		dst.ground_axiom_count = 0;
		if (!hash_set_init(dst.observations, src.observations.capacity)) {
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		}
		array_map<const Proof*, Proof*> proof_map(64);
		proof_map.keys[0] = src.empty_set_axiom;
		proof_map.values[0] = src.empty_set_axiom;
		dst.empty_set_axiom = src.empty_set_axiom;
		dst.empty_set_axiom->reference_count++;
		proof_map.keys[1] = src.maximality_axiom;
		proof_map.values[1] = src.maximality_axiom;
		dst.maximality_axiom = src.maximality_axiom;
		dst.maximality_axiom->reference_count++;
		proof_map.keys[2] = src.function_axiom;
		proof_map.values[2] = src.function_axiom;
		dst.function_axiom = src.function_axiom;
		dst.function_axiom->reference_count++;
		dst.NAME_ATOM = src.NAME_ATOM;
		dst.NAME_ATOM->reference_count++;
		if (!set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>::clone(src.sets, dst.sets, proof_map)) {
			core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			core::free(*dst.empty_set_axiom); core::free(*dst.maximality_axiom);
			core::free(*dst.function_axiom); core::free(*dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.disjunction_intro_nodes, src.disjunction_intro_nodes.capacity)) {
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			core::free(*dst.empty_set_axiom); core::free(*dst.maximality_axiom);
			core::free(*dst.function_axiom); core::free(*dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.negated_conjunction_nodes, src.negated_conjunction_nodes.capacity)) {
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			core::free(*dst.empty_set_axiom); core::free(*dst.maximality_axiom);
			core::free(*dst.function_axiom); core::free(*dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.implication_intro_nodes, src.implication_intro_nodes.capacity)) {
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			core::free(*dst.empty_set_axiom); core::free(*dst.maximality_axiom);
			core::free(*dst.function_axiom); core::free(*dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.existential_intro_nodes, src.existential_intro_nodes.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			core::free(*dst.empty_set_axiom); core::free(*dst.maximality_axiom);
			core::free(*dst.function_axiom); core::free(*dst.NAME_ATOM);
			return false;
		}

		for (const auto& entry : src.atoms) {
			unsigned int bucket = dst.atoms.table.index_to_insert(entry.key);
			if (!array_init(dst.atoms.values[bucket].key, entry.value.key.capacity)) {
				core::free(dst); return false;
			} else if (!array_init(dst.atoms.values[bucket].value, entry.value.value.capacity)) {
				core::free(dst.atoms.values[bucket].key);
				core::free(dst); return false;
			} else if (!::init(dst.atoms.table.keys[bucket], entry.key)) {
				core::free(dst.atoms.values[bucket].key);
				core::free(dst.atoms.values[bucket].value);
				core::free(dst); return false;
			}
			dst.atoms.table.size++;
			for (unsigned int id : entry.value.key)
				dst.atoms.values[bucket].key[dst.atoms.values[bucket].key.length++] = id;
			for (unsigned int id : entry.value.value)
				dst.atoms.values[bucket].value[dst.atoms.values[bucket].value.length++] = id;
		} for (const auto& entry : src.relations) {
			unsigned int bucket = dst.relations.table.index_to_insert(entry.key);
			if (!array_init(dst.relations.values[bucket].key, entry.value.key.capacity)) {
				core::free(dst); return false;
			} else if (!array_init(dst.relations.values[bucket].value, entry.value.value.capacity)) {
				core::free(dst.relations.values[bucket].key);
				core::free(dst); return false;
			}
			dst.relations.table.keys[bucket] = entry.key;
			dst.relations.table.size++;
			for (unsigned int id : entry.value.key)
				dst.relations.values[bucket].key[dst.relations.values[bucket].key.length++] = id;
			for (unsigned int id : entry.value.value)
				dst.relations.values[bucket].value[dst.relations.values[bucket].value.length++] = id;
		} for (Proof* observation : src.observations) {
			Proof* new_observation;
			if (!Proof::clone(observation, new_observation, proof_map)) {
				core::free(dst);
				return false;
			}
			dst.observations.add(new_observation);
		} for (unsigned int i = 0; i < src.ground_concept_capacity; i++) {
			if (src.ground_concepts[i].types.keys == nullptr) continue;
			if (!concept<ProofCalculus>::clone(src.ground_concepts[i], dst.ground_concepts[i], proof_map)) {
				core::free(dst);
				return false;
			}
			dst.ground_axiom_count++;
		} for (auto& entry : src.disjunction_intro_nodes) {
			unsigned int index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			dst.disjunction_intro_nodes[dst.disjunction_intro_nodes.length].value = proof_map.values[index];
			dst.disjunction_intro_nodes[dst.disjunction_intro_nodes.length++].key = entry.key;
			entry.key->reference_count++;
		} for (auto& entry : src.negated_conjunction_nodes) {
			unsigned int index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			dst.negated_conjunction_nodes[dst.negated_conjunction_nodes.length].value = proof_map.values[index];
			dst.negated_conjunction_nodes[dst.negated_conjunction_nodes.length++].key = entry.key;
			entry.key->reference_count++;
		} for (auto& entry : src.implication_intro_nodes) {
			unsigned int index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			dst.implication_intro_nodes[dst.implication_intro_nodes.length].value = proof_map.values[index];
			dst.implication_intro_nodes[dst.implication_intro_nodes.length++].key = entry.key;
			entry.key->reference_count++;
		} for (auto& entry : src.existential_intro_nodes) {
			unsigned int index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			dst.existential_intro_nodes[dst.existential_intro_nodes.length].value = proof_map.values[index];
			dst.existential_intro_nodes[dst.existential_intro_nodes.length++].key = entry.key;
			entry.key->reference_count++;
		}

		for (unsigned int i = 0; i < dst.ground_concept_capacity; i++) {
			if (dst.ground_concepts[i].types.keys == nullptr) continue;
			auto& c = dst.ground_concepts[i];
			for (unsigned int i = 0; i < c.definitions.length; i++) {
				Proof* definition = c.definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					unsigned int arity = 1;
					Formula* set_formula = definition->formula->binary.right->quantifier.operand;
					while (set_formula->type == FormulaType::LAMBDA) {
						set_formula = set_formula->quantifier.operand;
						arity++;
					}

					for (unsigned int j = i + 1; j < c.definitions.length; j++) {
						Proof* other_definition = c.definitions[j];
						if (other_definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
						Formula* other_set_formula = other_definition->formula->binary.right->quantifier.operand;
						while (other_set_formula->type == FormulaType::LAMBDA)
							other_set_formula = other_set_formula->quantifier.operand;

						bool dummy; unsigned int antecedent_set, consequent_set;
						Proof* axiom = dst.sets.template get_subset_axiom<false, false>(
								set_formula, other_set_formula, arity,
								antecedent_set, consequent_set, dummy, dummy);
						axiom->reference_count++;

						axiom = dst.sets.template get_subset_axiom<false, false>(
								other_set_formula, set_formula, arity,
								consequent_set, antecedent_set, dummy, dummy);
						axiom->reference_count++;
					}
				}
			}
		}
		return true;
	}

	template<typename... Args>
	Proof* add_formula(Formula* formula, set_changes<Formula>& set_diff, unsigned int& new_constant, Args&&... args)
	{
		new_constant = 0;

		Formula* new_formula = preprocess_formula(formula);
		if (new_formula == NULL) return nullptr;

		array_map<unsigned int, unsigned int> variable_map(16);
		Formula* canonicalized = Canonicalizer::canonicalize(*new_formula, variable_map);
		core::free(*new_formula); if (new_formula->reference_count == 0) core::free(new_formula);
		if (canonicalized == NULL) return nullptr;

/* TODO: for debugging; delete this */
print_axioms(stderr);
print("canonicalized: ", stderr); print(*canonicalized, stderr); print('\n', stderr);
		Proof* new_proof = make_proof<false, true, true>(canonicalized, set_diff, new_constant, std::forward<Args>(args)...);
print_axioms(stderr);
if (new_proof != NULL) {
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put((unsigned int) built_in_predicates::UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof<built_in_predicates, typename ProofCalculus::ProofCanonicalizer>(*new_proof, expected_conclusion))
fprintf(stderr, "add_formula WARNING: `check_proof` failed.\n");
core::free(*expected_conclusion); if (expected_conclusion->reference_count == 0) core::free(expected_conclusion);
}
		core::free(*canonicalized);
		if (canonicalized->reference_count == 0)
			core::free(canonicalized);
		if (new_proof == NULL) {
			return nullptr;
		} else if (!observations.add(new_proof)) {
			free_proof(new_proof, set_diff, std::forward<Args>(args)...);
			return nullptr;
		}
		return new_proof;
	}

	template<bool FreeProof = true>
	void remove_formula(Proof* proof, set_changes<Formula>& set_diff) {
		theory::changes changes;
		observations.remove(proof);
		if (!get_theory_changes(*proof, changes)) return;
		subtract_changes(changes, set_diff);
		if (FreeProof) {
			core::free(*proof); if (proof->reference_count == 0) core::free(proof);
		}
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
			if (sets.sets[i].size_axioms.data == nullptr)
				continue;
			Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

			array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				quantifiers[quantifiers.length++] = set_formula;
				set_formula = set_formula->quantifier.operand;
			}

			for (unsigned int j = sets.sets[i].element_count(); j > 0; j--) {
				tuple_element* element = (tuple_element*) alloca(sizeof(tuple_element) * sets.sets[i].arity);
				const tuple_element* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * (j - 1));
				for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) {
					if (!::init(element[k], element_src[k])) {
						for (uint_fast8_t l = 0; l < k; l++) core::free(element[l]);
					}
				}
				array<instantiation_tuple> possible_values(1);
				if (!::init(possible_values[0], sets.sets[i].arity)) {
					for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) core::free(element[k]);
					return false;
				}
				possible_values.length = 1;
				for (uint_fast8_t k = 0; k < possible_values[0].length; k++) {
					core::free(possible_values[0].values[k]);
					if (element[k].type == tuple_element_type::CONSTANT) {
						possible_values[0].values[k].type = instantiation_type::CONSTANT;
						possible_values[0].values[k].constant = element[k].constant;
					} else if (element[k].type == tuple_element_type::NUMBER) {
						possible_values[0].values[k].type = instantiation_type::NUMBER;
						possible_values[0].values[k].number = element[k].number;
					} else if (element[k].type == tuple_element_type::STRING) {
						possible_values[0].values[k].type = instantiation_type::STRING;
						if (!core::init(possible_values[0].values[k].str, element[k].str)) {
							possible_values[0].values[k].type = instantiation_type::CONSTANT;
							for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) core::free(element[k]);
							core::free(possible_values[0]); return false;
						}
					}
				}

				sets.sets[i].remove_element_at(j - 1);
				default_prover prover(sets);
				if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values, prover)) {
					print("theory.are_elements_provable ERROR: The element ", stderr);
					if (sets.sets[i].arity == 1)
						print(element[0], stderr);
					else print(element, sets.sets[i].arity, stderr);
					print(" does not provably belong to set with ID ", stderr);
					print(i, stderr); print(".\n", stderr);
					success = false;
				}
				for (auto& element : possible_values) core::free(element);
				for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++)
					move(element[k], sets.sets[i].elements[sets.sets[i].elements.length++]);
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
		core::free(*atom); core::free(atom);
		if (contains && type_instances.key.length != 0)
			return true;
		return false;
	}

	Term* get_arg(unsigned int event_constant, unsigned int arg) const {
		bool contains;
		Proof* proof = ground_concepts[event_constant - new_constant_offset].function_values.get(arg, contains);
		if (contains) return proof->formula->binary.right;

		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;
			for (Proof* proof : ground_concepts[i].definitions) {
				if (proof->formula->binary.right->type == TermType::UNARY_APPLICATION
				 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
				 && proof->formula->binary.right->binary.left->constant == arg
				 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
				 && proof->formula->binary.right->binary.right->constant == event_constant)
				{
					return proof->formula->binary.left;
				}
			}
		}
		return nullptr;
	}

	inline bool is_name_event(unsigned int event_constant) const {
		if (event_constant < new_constant_offset) return false;
		for (const auto& entry : ground_concepts[event_constant - new_constant_offset].types) {
			if (entry.key.type == TermType::UNARY_APPLICATION
			 && entry.key.binary.left->type == TermType::CONSTANT && entry.key.binary.left->constant == (unsigned int) built_in_predicates::NAME
			 && entry.key.binary.right->type == TermType::VARIABLE && entry.key.binary.right->variable == 1)
			{
				return true;
			}
		}
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

		/* check to make sure the args are type-correct */
		if (atom.binary.left->type == TermType::CONSTANT && atom.binary.left->constant == (unsigned int) built_in_predicates::NAME) {
			/* arg1 must be a non-name, non-set constant and arg2 must be a string */
			Term* arg1 = get_arg(arg, (unsigned int) built_in_predicates::ARG1);
			if (arg1 != nullptr) {
				if (arg1->type != TermType::CONSTANT || is_name_event(arg1->constant) || is_provably_a_set(arg1->constant))
					return false;
			}
			Term* arg2 = get_arg(arg, (unsigned int) built_in_predicates::ARG2);
			if (arg2 != nullptr && arg2->type != TermType::STRING)
				return false;
			/* this event also cannot be the arg1 of a name event */
			for (const Proof* definition : ground_concepts[arg - new_constant_offset].definitions) {
				Formula* right = definition->formula->binary.right;
				if (right->type == TermType::UNARY_APPLICATION
				 && right->binary.left->type == TermType::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				 && right->binary.right->type == TermType::CONSTANT && is_name_event(right->binary.right->constant))
					return false;
			}
		} else {
			/* arg1 can be anything and arg2 must be a non-string */
			bool contains;
			Proof* proof = ground_concepts[arg - new_constant_offset].function_values.get((unsigned int) built_in_predicates::ARG2, contains);
			if (contains && proof->formula->binary.right->type == TermType::STRING)
				return false;
		}

		Formula* lifted_literal;
		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;
		if (Negated) {
			lifted_literal = Formula::new_not(lifted_atom);
			if (lifted_literal == NULL) {
				core::free(*lifted_atom); core::free(lifted_atom);
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
			core::free(*lifted_literal); core::free(lifted_literal);
			return false;
		}

		pair<array<unsigned int>, array<unsigned int>>& instance_pair = atoms.get(*lifted_literal, contains, bucket);
		if (!contains) {
			if (!array_init(instance_pair.key, 8)) {
				core::free(*lifted_literal); core::free(lifted_literal);
				return false;
			} else if (!array_init(instance_pair.value, 8)) {
				core::free(instance_pair.key);
				core::free(*lifted_literal); core::free(lifted_literal);
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
			core::free(*lifted_literal); core::free(lifted_literal);
			return false;
		}

		add_sorted<false>(instances, arg);
		ground_types.keys[ground_types.size] = *lifted_literal;
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		ground_axiom_count++;
		core::free(*lifted_literal); core::free(lifted_literal);

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
			core::free(*atom); core::free(atom);
			return false;
		}
		core::free(*atom); core::free(atom);
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
				core::free(*lifted_atom); core::free(lifted_atom);
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
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		core::free(ground_types.keys[index]);
		ground_types.remove_at(index);
		ground_axiom_count--;
		core::free(*lifted_literal); core::free(lifted_literal);
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
		core::free(*atom); core::free(atom);
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

		~changes() { free_helper(); }

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

		static inline void free(changes& changes) {
			changes.free_helper();
			core::free(changes.list);
		}

private:
		inline void free_helper() {
			for (change& c : list)
				core::free(c);
		}
	};

	static inline bool init(changes& changes) {
		return array_init(changes.list, 8);
	}

	template<typename... Args>
	bool add_changes(theory::changes& changes, set_changes<Formula>& set_diff, Args&&... visitor)
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
				if (!try_init_concept(c.binary_atom.key.predicate, 2)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<false, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) return false;
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate, 2)
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
							if (!try_init_concept(predicate, 1)) return false;
						}
					}
					if (!sets.add_subset_axiom(c.axiom, set_diff, std::forward<Args>(visitor)...)) return false;
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
					bool is_new;
					if (!sets.get_set_id(set_formula, arity, set_id, is_new, set_diff, std::forward<Args>(visitor)...)
					 || !sets.sets[set_id].set_size_axiom(c.axiom)
					 || (is_new && !check_new_set_membership<false>(set_id, std::forward<Args>(visitor)...)))
						return false;
				}
				continue;
			case change_type::DEFINITION:
				if (!add_definition<false>(c.axiom, set_diff, std::forward<Args>(visitor)...))
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
	bool subtract_changes(theory::changes& changes, set_changes<Formula>& set_diff, Args&&... visitor)
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
						sets.free_subset_axiom(c.axiom, set_diff, freeable_set_size_axioms, std::forward<Args>(visitor)...);
					} else {
						core::free(*c.axiom);
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
						check_old_set_membership(set_id, std::forward<Args>(visitor)...);
						on_free_set(set_id, sets, set_diff, std::forward<Args>(visitor)...);
						sets.free_set(set_id);
					}
					continue;
				}
			case change_type::DEFINITION:
				remove_definition(c.axiom, freeable_set_size_axioms, set_diff, std::forward<Args>(visitor)...);
				continue;
			case change_type::FUNCTION_VALUE:
				remove_function_value(c.axiom, std::forward<Args>(visitor)...);
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				implication_intro_nodes.remove(index_of(implication_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				on_undo_filter_operands(c.intro_node.key, std::forward<Args>(visitor)...);
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				negated_conjunction_nodes.remove(index_of(negated_conjunction_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				on_undo_filter_operands(c.intro_node.key, std::forward<Args>(visitor)...);
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
					on_undo_filter_constants(*this, c.intro_node.key->quantifier.operand, c.intro_node.key->quantifier.variable, std::forward<Args>(visitor)...);
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				disjunction_intro_nodes.remove(index_of(disjunction_intro_nodes, c.intro_node.value));
				c.intro_node.key->reference_count--;
				on_undo_filter_operands(c.intro_node.key, std::forward<Args>(visitor)...);
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
				 && proof.formula->binary.right->type == TermType::NUMBER)
				{
					if (reference_count != 1) return true;
					return changes.add(change_type::SET_SIZE_AXIOM, &proof);
				} else if (proof.formula->binary.left->type == TermType::CONSTANT) {
					if (reference_count != 1) return true;
					/* make sure we keep track of decrementing the reference counts of bidirectional subset edges */
					if (proof.formula->binary.right->type == TermType::LAMBDA) {
						unsigned int arity = 1;
						Formula* set_formula = proof.formula->binary.right->quantifier.operand;
						while (set_formula->type == FormulaType::LAMBDA) {
							set_formula = set_formula->quantifier.operand;
							arity++;
						}

						unsigned int concept_id = proof.formula->binary.left->constant;
						for (Proof* other_definition : ground_concepts[concept_id - new_constant_offset].definitions) {
							if (proof.formula->binary.right == other_definition->formula->binary.right
							 || other_definition->formula->binary.right->type != FormulaType::LAMBDA)
								continue;
							if (!reference_counts.ensure_capacity(reference_counts.size + 2))
								return false;
							Formula* other_set_formula = other_definition->formula->binary.right->quantifier.operand;
							while (other_set_formula->type == FormulaType::LAMBDA)
								other_set_formula = other_set_formula->quantifier.operand;

							Proof* axiom = sets.template get_existing_subset_axiom<false>(other_set_formula, set_formula, arity);
							unsigned int index = reference_counts.index_of(axiom);
							if (index == reference_counts.size) {
								reference_counts.keys[index] = axiom;
								reference_counts.values[index] = axiom->reference_count;
								reference_counts.size++;
							}
							reference_counts.values[index]--;

							axiom = sets.template get_existing_subset_axiom<false>(set_formula, other_set_formula, arity);
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
			} else if (&proof != empty_set_axiom
					&& &proof != maximality_axiom
					&& &proof != function_axiom)
			{
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
		case FormulaType::NUMBER:
			return term->number == atom->number;
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
		case FormulaType::NUMBER:
			return first->number == second->number;
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
		case FormulaType::NUMBER:
			return first->number == second->number;
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
		{
			if (!first_quantifiers.add(first) || !second_quantifiers.add(second))
				return false;
			bool result = unify(first->quantifier.operand, second->quantifier.operand, first_quantifiers, first_unifications, second_quantifiers, second_unifications);
			first_quantifiers.length--;
			second_quantifiers.length--;
			return result;
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
		fprintf(stderr, "theory.unify ERROR: Unrecognized FormulaType.\n");
		return false;
	}

	struct default_prover {
		const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets;

		default_prover(
				const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets) : sets(sets) { }

		inline bool get_provable_elements(unsigned int set_id, hash_set<tuple>& provable_elements) {
			return sets.get_provable_elements(set_id, provable_elements);
		}
	};

	struct set_membership_prover {
		const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets;
		const array_map<tuple, array<unsigned int>>& new_elements;

		set_membership_prover(
				const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets,
				const array_map<tuple, array<unsigned int>>& new_elements) : sets(sets), new_elements(new_elements) { }

		inline bool get_provable_elements(unsigned int set_id, hash_set<tuple>& provable_elements)
		{
			if (!sets.get_provable_elements(set_id, provable_elements))
				return false;

			for (const auto& entry : new_elements) {
				for (unsigned int other_set : entry.value) {
					if (sets.sets[set_id].descendants.contains(other_set)) {
						if (!provable_elements.check_size()) return false;

						bool contains;
						unsigned int index = provable_elements.index_of(entry.key, contains);
						if (!contains) {
							if (!::init(provable_elements.keys[index], entry.key))
								return false;
							provable_elements.size++;
						}
					}
				}
			}
			return true;
		}
	};

	template<bool Contradiction, typename Prover>
	bool is_provable_by_theorem_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values,
			Prover& prover) const
	{
		unsigned int old_size = new_possible_values.length;
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr) continue;

			array<Formula*> second_quantifiers(sets.sets[i].arity);
			Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				second_quantifiers.add(set_formula);
				set_formula = set_formula->quantifier.operand;
			}

			array_view<Formula*> set_formula_conjuncts(
					(set_formula->type == FormulaType::AND) ? set_formula->array.operands : &set_formula,
					(set_formula->type == FormulaType::AND) ? set_formula->array.length : 1);
			for (unsigned int l = 0; l < set_formula_conjuncts.length; l++) {
				Formula* set_formula_conjunct = set_formula_conjuncts.array[l];
				if (Contradiction) {
					if (set_formula_conjunct->type != FormulaType::NOT) continue;
					set_formula_conjunct = set_formula_conjunct->unary.operand;
				}
				array_map<Formula*, Term*> first_unifications(4);
				array_map<Formula*, Term*> second_unifications(4);
				if (!unify(formula, set_formula_conjunct, quantifiers, first_unifications, second_quantifiers, second_unifications))
					continue;

				hash_set<tuple> provable_elements(16);
				if (!prover.get_provable_elements(i, provable_elements))
					return false;

				for (unsigned int j = 0; j < possible_values.length; j++) {
					const instantiation_tuple& values = possible_values[j];
					if (!new_possible_values.ensure_capacity(new_possible_values.length + provable_elements.size + 1)) {
						for (tuple& tup : provable_elements) core::free(tup);
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
						if (!::init(new_new_value, new_value)) {
							for (tuple& tup : provable_elements) core::free(tup);
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
							Term* constant = nullptr;
							if (element[unification.value->variable - 1].type == tuple_element_type::CONSTANT)
								constant = Term::new_constant(element[unification.value->variable - 1].constant);
							else if (element[unification.value->variable - 1].type == tuple_element_type::NUMBER)
								constant = Term::new_number(element[unification.value->variable - 1].number);
							else if (element[unification.value->variable - 1].type == tuple_element_type::STRING)
								constant = Term::new_string(element[unification.value->variable - 1].str);
							if (constant == nullptr) {
								for (tuple& tup : provable_elements) core::free(tup);
								return false;
							}
							if (!new_new_value.unify_value(unification.key->quantifier.variable - 1, constant)) {
								core::free(*constant); core::free(constant);
								unifies = false;
								break;
							}
							core::free(*constant); core::free(constant);
						}
						if (!unifies) {
							core::free(new_new_value);
							continue;
						}

						for (const auto& unification : second_unifications) {
							if (element[unification.key->quantifier.variable - 1].type == tuple_element_type::CONSTANT) {
								if (unification.value->type != TermType::CONSTANT
								 || unification.value->constant != element[unification.key->quantifier.variable - 1].constant)
								{
									unifies = false;
									break;
								}
							} else if (element[unification.key->quantifier.variable - 1].type == tuple_element_type::NUMBER) {
								if (unification.value->type != TermType::NUMBER
								 || unification.value->number != element[unification.key->quantifier.variable - 1].number)
								{
									unifies = false;
									break;
								}
							} else if (element[unification.key->quantifier.variable - 1].type == tuple_element_type::STRING) {
								if (unification.value->type != TermType::STRING
								 || unification.value->str != element[unification.key->quantifier.variable - 1].str)
								{
									unifies = false;
									break;
								}
							}
						}
						if (!unifies) {
							core::free(new_new_value);
							continue;
						}

						new_possible_values.length++;
					}
				}
				for (tuple& tup : provable_elements) core::free(tup);
			}
		}
		if (new_possible_values.length > old_size) {
			insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
			array<instantiation_tuple> temp(new_possible_values.length);
			set_union(temp.data, temp.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
			swap(temp, new_possible_values);
			for (auto& element : temp) core::free(element);
		}
		return true;
	}

	template<typename Prover>
	bool for_all_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values,
			Prover& prover) const
	{
		array<unsigned int> new_variables(4);
		while (formula->type == TermType::FOR_ALL) {
			if (!new_variables.add(formula->quantifier.variable))
				return false;
			formula = formula->quantifier.operand;
		}

		if (formula->type == TermType::IF_THEN) {
			Formula* antecedent = formula->binary.left;
			Formula* consequent = formula->binary.right;

			/* find sets that unify with `antecedent` */
			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axioms.data == nullptr || sets.sets[i].arity < new_variables.length) continue;

				array<Formula*> second_quantifiers(sets.sets[i].arity);
				Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					second_quantifiers.add(set_formula);
					set_formula = set_formula->quantifier.operand;
				}
				array_map<Formula*, Term*> first_unifications(4);
				array_map<Formula*, Term*> second_unifications(4);
				if (!unify(set_formula, antecedent, second_quantifiers, second_unifications, quantifiers, first_unifications))
					continue;

				/* make sure the newly declared variable in `formula` is covered by a variable in `set_formula` */
				bool are_variables_mapped = true;
				for (unsigned int new_variable : new_variables) {
					bool is_variable_mapped = false;
					for (unsigned int i = 0; !is_variable_mapped && i < second_unifications.size; i++)
						if (second_unifications.values[i]->type == TermType::VARIABLE && second_unifications.values[i]->variable == new_variable)
							is_variable_mapped = true;
					if (!is_variable_mapped) {
						are_variables_mapped = false;
						break;
					}
				}
				if (!are_variables_mapped) continue;

				hash_set<tuple> provable_elements(16);
				if (!prover.get_provable_elements(i, provable_elements))
					return false;
				if (provable_elements.size != sets.sets[i].set_size) {
					for (tuple& tup : provable_elements) core::free(tup);
					continue;
				}

				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& values = possible_values[i];
					array<instantiation_tuple> temp_possible_values(4);
					instantiation_tuple& new_values = temp_possible_values[0];
					if (!::init(new_values, values)) {
						for (tuple& tup : provable_elements) core::free(tup);
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
						for (auto& element : temp_possible_values) core::free(element);
						continue;
					}

					/* check that `consequent` is provable for values in `provable_elements` */
					for (const tuple& element : provable_elements) {
						Term** src_variables = (Term**) alloca(sizeof(Term*) * second_unifications.size * 2);
						for (unsigned int i = 0; i < second_unifications.size; i++) {
							src_variables[i] = second_unifications.values[i];
							src_variables[i]->reference_count++;
						}

						Term** dst_constants = src_variables + second_unifications.size;
						for (unsigned int i = 0; i < second_unifications.size; i++) {
							const tuple_element& value = element[second_unifications.keys[i]->quantifier.variable - 1];
							if (value.type == tuple_element_type::CONSTANT)
								dst_constants[i] = Term::new_constant(value.constant);
							else if (value.type == tuple_element_type::NUMBER)
								dst_constants[i] = Term::new_number(value.number);
							else if (value.type == tuple_element_type::STRING)
								dst_constants[i] = Term::new_string(value.str);
							if (dst_constants[i] == nullptr) {
								for (unsigned int j = 0; j < second_unifications.size; j++) { core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]); }
								for (unsigned int j = 0; j < i; j++) { core::free(*dst_constants[j]); core::free(dst_constants[j]); }
								for (auto& element : temp_possible_values) core::free(element);
								for (tuple& tup : provable_elements) core::free(tup);
								return false;
							}
						}

						Formula* substituted_consequent = substitute_all(consequent, src_variables, dst_constants, element.length);
						for (unsigned int j = 0; j < second_unifications.size; j++) { core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]); }
						for (unsigned int j = 0; j < second_unifications.size; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }

						if (!is_provable_without_abduction<false>(substituted_consequent, quantifiers, temp_possible_values, prover)) {
							core::free(*substituted_consequent); if (substituted_consequent->reference_count == 0) core::free(substituted_consequent);
							break;
						}
						core::free(*substituted_consequent); if (substituted_consequent->reference_count == 0) core::free(substituted_consequent);
					}

					if (temp_possible_values.length != 0) {
						array<instantiation_tuple> temp(new_possible_values.length + temp_possible_values.length);
						set_union(temp, new_possible_values, temp_possible_values);
						swap(temp, new_possible_values);
						for (auto& element : temp_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
					}
				}
				for (tuple& tup : provable_elements) core::free(tup);
			}
		}
		return true;
	}

	template<bool Contradiction, typename Prover>
	bool exists_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<instantiation_tuple>& possible_values,
			array<instantiation_tuple>& new_possible_values,
			Prover& prover) const
	{
		Formula* quantified = formula->quantifier.operand;
		unsigned int variable = formula->quantifier.variable;

		array<instance> constants(ground_concept_capacity + 1);
		array<hol_number> numbers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != nullptr) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length++].constant = new_constant_offset + i;

				for (const auto& entry : ground_concepts[i].function_values) {
					Term* constant = entry.value->formula->binary.right;
					if (constant->type == TermType::NUMBER) {
						if (!numbers.contains(constant->number) && !numbers.add(constant->number)) return false;
					} else if (constant->type == TermType::STRING) {
						bool contains = false;
						for (const string* str : strings)
							if (str == &constant->str || *str == constant->str) { contains = true; break; }
						if (!contains && !strings.add(&constant->str)) return false;
					}
				}
			}
		} for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr) continue;
			hol_number number;
			number.integer = sets.sets[i].set_size;
			number.decimal = 0;
			if (!numbers.contains(number) && !numbers.add(number))
				return false;
		}
		if (!constants.ensure_capacity(constants.length + numbers.length + strings.length))
			return false;
		for (hol_number number : numbers) {
			constants[constants.length].type = instance_type::NUMBER;
			constants[constants.length++].number = number;
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
			} else if (id.type == instance_type::NUMBER) {
				constant = Formula::new_number(id.number);
			} else if (id.type == instance_type::STRING) {
				constant = Formula::new_string(*id.str);
			}
			if (constant == nullptr) {
				core::free(*var); if (var->reference_count == 0) core::free(var);
				return false;
			}
			Formula* substituted = substitute(quantified, var, constant);
			core::free(*constant); if (constant->reference_count == 0) core::free(constant);
			if (substituted == nullptr) {
				core::free(*var); if (var->reference_count == 0) core::free(var);
				return false;
			}

			array<instantiation_tuple> copy(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(copy[copy.length], possible_values[j])) {
					for (auto& element : copy) core::free(element);
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return false;
				}
				copy.length++;
			}

			bool result = is_provable_without_abduction<Contradiction>(substituted, quantifiers, copy, prover);
			core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
			if (!result) {
				for (auto& element : copy) core::free(element);
				continue;
			}

			array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
			set_union(union_result, new_possible_values, copy);
			swap(union_result, new_possible_values);
			for (auto& element : union_result) core::free(element);
			for (auto& element : copy) core::free(element);
		}
		core::free(*var); if (var->reference_count == 0) core::free(var);
		return true;
	}

	template<bool Contradiction, typename Prover>
	bool is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			array<instantiation_tuple>& possible_values,
			Prover& prover) const
	{
		/* consider if this formula provable by other extensional edges (i.e. universally-quantified theorems) */
		array<instantiation_tuple> new_possible_values(possible_values.length);
		if (!is_provable_by_theorem_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values, prover)) {
			for (auto& element : new_possible_values) core::free(element);
			return false;
		}

		if (formula->type == FormulaType::UNARY_APPLICATION) {
			if (formula->binary.left->type == TermType::CONSTANT && formula->binary.left->constant == (unsigned int) built_in_predicates::NUMBER) {
				if (Contradiction) {
					if (formula->binary.right->type == TermType::CONSTANT) {
						for (auto& element : new_possible_values) core::free(element);
						return true;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							if (!::init(temp[temp.length], possible_values[i])) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : possible_values) core::free(element);
								possible_values.clear(); return false;
							} else if (!temp[temp.length].unify_not_number(formula->binary.right->variable - 1)) {
								core::free(temp[temp.length]); continue;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
						for (auto& element : possible_values) core::free(element);
						swap(new_possible_values, possible_values);
						return (possible_values.length > 0);
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				} else {
					if (formula->binary.right->type == TermType::NUMBER) {
						for (auto& element : new_possible_values) core::free(element);
						return true;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							if (!::init(temp[temp.length], possible_values[i])) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : possible_values) core::free(element);
								possible_values.clear(); return false;
							} else if (!temp[temp.length].unify_number(formula->binary.right->variable - 1)) {
								core::free(temp[temp.length]); continue;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
						for (auto& element : possible_values) core::free(element);
						swap(new_possible_values, possible_values);
						return (possible_values.length > 0);
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				}
			}

			if (formula->binary.right->type != TermType::CONSTANT
			 && formula->binary.right->type != TermType::VARIABLE)
			{
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
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
						for (auto& element : temp) core::free(element);
						continue;
					}

					if (formula->binary.right->type == TermType::VARIABLE) {
						for (unsigned int constant : constants) {
							instantiation_tuple& current_new_values = temp[temp.length];
							if (!::init(current_new_values, new_values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								return false;
							}
							Term* term = Term::new_constant(constant);
							if (term == nullptr) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								return false;
							}
							unifies = current_new_values.unify_value(formula->binary.right->variable - 1, term);
							core::free(*term); core::free(term);
							if (!unifies) {
								core::free(current_new_values);
								continue;
							}
							temp.length++;
						}
						if (temp.length > 0)
							insertion_sort(temp, default_sorter());
					} else if (constants.contains(formula->binary.right->constant)) {
						if (!::init(temp[0], new_values)) {
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}
						temp.length++;
					}
					if (temp.length != 0) {
						array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
						set_union(union_result, new_possible_values, temp);
						swap(new_possible_values, union_result);
						for (auto& element : temp) core::free(element);
						for (auto& element : union_result) core::free(element);
					}
				}
			}

			/* this could also be a set membership declaration */
			Term* left_most = formula;
			array<Term*> args(4);
			while (left_most->type == TermType::UNARY_APPLICATION)
			{
				if (left_most->binary.right->type != TermType::CONSTANT && left_most->binary.right->type != TermType::NUMBER
				 && left_most->binary.right->type != TermType::STRING && left_most->binary.right->type != TermType::VARIABLE)
				{
					if (Contradiction) {
						for (auto& element : new_possible_values) core::free(element);
						return true;
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				}
				if (!args.add(left_most->binary.right)) {
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : possible_values) core::free(element);
					possible_values.clear(); return false;
				}
				left_most = left_most->binary.left;
			}
			reverse(args);
			if (left_most->type == TermType::VARIABLE || (left_most->type == TermType::CONSTANT && left_most->constant >= new_constant_offset)) {
				for (unsigned int i = 0; i < ground_concept_capacity; i++) {
					if (ground_concepts[i].types.keys == nullptr)
						continue;
					bool contains;
					Formula* set_formula = ground_concepts[i].definitions[0]->formula->binary.right->quantifier.operand;
					if (left_most->type == TermType::CONSTANT && left_most->constant != new_constant_offset + i)
						continue;
					unsigned int set_id = sets.set_ids.get(*set_formula, contains);
					if (!contains) continue;
					if (!Contradiction && sets.sets[set_id].arity != args.length) continue;

					hash_set<tuple> provable_elements(16);
					if (!prover.get_provable_elements(set_id, provable_elements)) {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
					if (Contradiction && sets.sets[set_id].arity != args.length) {
						/* since the arity of the elements don't match, this formula is trivially disprovable */
						for (auto& tup : provable_elements) core::free(tup);
						for (auto& element : new_possible_values) core::free(element);
						return true;
					}
					if (Contradiction && provable_elements.size != sets.sets[set_id].set_size) {
						/* we can only prove a contradiction by exclusion */
						for (auto& tup : provable_elements) core::free(tup);
						continue;
					}

					for (unsigned int j = 0; j < possible_values.length; j++) {
						instantiation_tuple new_values(possible_values[j]);
						if (left_most->type == TermType::VARIABLE && !new_values.unify_value(left_most->variable - 1, set_formula->binary.left))
							continue;

						array<instantiation_tuple> temp(max(1u, provable_elements.size));
						if (Contradiction) {
							bool is_provable_element = false;
							for (const tuple& provable_element : provable_elements) {
								bool unifies_with_element = true;
								for (unsigned int k = 0; k < args.length; k++) {
									if (args[k]->type == TermType::CONSTANT) {
										if (provable_element[k].type != tuple_element_type::CONSTANT
										 || provable_element[k].constant != args[k]->constant)
										{
											unifies_with_element = false;
											break;
										}
									} else if (args[k]->type == TermType::NUMBER) {
										if (provable_element[k].type != tuple_element_type::NUMBER
										 || provable_element[k].number != args[k]->number)
										{
											unifies_with_element = false;
											break;
										}
									} else if (args[k]->type == TermType::STRING) {
										if (provable_element[k].type != tuple_element_type::STRING
										 || provable_element[k].str != args[k]->str)
										{
											unifies_with_element = false;
											break;
										}
									} else {
										Term* constant = nullptr;
										if (provable_element[k].type == tuple_element_type::CONSTANT)
											constant = Term::new_constant(provable_element[k].constant);
										else if (provable_element[k].type == tuple_element_type::NUMBER)
											constant = Term::new_number(provable_element[k].number);
										else if (provable_element[k].type == tuple_element_type::STRING)
											constant = Term::new_string(provable_element[k].str);
										if (new_values.antiunify_value(args[k]->variable - 1, constant)) {
											core::free(*constant); core::free(constant);
											unifies_with_element = false; break;
										}
										core::free(*constant); core::free(constant);
									}
								}
								if (unifies_with_element) {
									is_provable_element = true;
									break;
								}
							}
							if (!is_provable_element) {
								if (!::init(temp[temp.length], new_values)) {
									for (auto& tup : provable_elements) core::free(tup);
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : possible_values) core::free(element);
									possible_values.clear(); return false;
								}
								temp.length++;
							}
						} else {
							for (const tuple& provable_element : provable_elements) {
								instantiation_tuple& new_new_values = temp[temp.length];
								if (!::init(new_new_values, new_values)) {
									for (auto& tup : provable_elements) core::free(tup);
									for (auto& element : temp) core::free(element);
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : possible_values) core::free(element);
									possible_values.clear(); return false;
								}

								bool unifies = true;
								for (unsigned int k = 0; k < args.length; k++) {
									if (args[k]->type == TermType::CONSTANT) {
										if (provable_element[k].type != tuple_element_type::CONSTANT
										 || provable_element[k].constant != args[k]->constant)
										{
											unifies = false;
											break;
										}
									} else if (args[k]->type == TermType::NUMBER) {
										if (provable_element[k].type != tuple_element_type::NUMBER
										 || provable_element[k].number != args[k]->number)
										{
											unifies = false;
											break;
										}
									} else if (args[k]->type == TermType::STRING) {
										if (provable_element[k].type != tuple_element_type::STRING
										 || provable_element[k].str != args[k]->str)
										{
											unifies = false;
											break;
										}
									} else {
										Term* constant = nullptr;
										if (provable_element[k].type == tuple_element_type::CONSTANT)
											constant = Term::new_constant(provable_element[k].constant);
										else if (provable_element[k].type == tuple_element_type::NUMBER)
											constant = Term::new_number(provable_element[k].number);
										else if (provable_element[k].type == tuple_element_type::STRING)
											constant = Term::new_string(provable_element[k].str);
										if (!new_values.unify_value(args[k]->variable - 1, constant)) {
											core::free(*constant); core::free(constant);
											unifies = false; break;
										}
										core::free(*constant); core::free(constant);
									}
								}
								if (unifies) {
									temp.length++;
								} else {
									core::free(new_new_values);
								}
							}
						}

						if (temp.length != 0) {
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
					}
					for (auto& tup : provable_elements) core::free(tup);
				}

				for (unsigned int k = 0; k < args.length; k++) {
					if (args[k]->type == TermType::VARIABLE) {
						if (Contradiction) {
							array<instantiation_tuple> temp(max((size_t) 1, possible_values.length));
							for (const instantiation_tuple& value : possible_values) {
								instantiation_tuple& new_value = temp[temp.length];
								if (!::init(new_value, value)) {
									for (auto& element : temp) core::free(element);
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : possible_values) core::free(element);
									possible_values.clear(); return false;
								} if (!new_value.unify_value(args[k]->variable - 1, left_most)) {
									core::free(new_value);
									continue;
								}
								temp.length++;
							}

							if (temp.length != 0) {
								array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
								set_union(union_result, new_possible_values, temp);
								swap(new_possible_values, union_result);
								for (auto& element : temp) core::free(element);
								for (auto& element : union_result) core::free(element);
							}
						} else {
							for (unsigned int i = 0; i < new_possible_values.length; i++) {
								if (!new_possible_values[i].antiunify_value(args[k]->variable - 1, left_most)) {
									core::free(new_possible_values[i]);
									new_possible_values.remove(i--);
								}
							}
						}
					}
				}
				if (left_most->type == TermType::VARIABLE) {
					for (unsigned int k = 0; k < args.length; k++) {
						if (args[k]->type == TermType::VARIABLE)
							/* avoid redundant computation since we already checked for this case earlier */
							continue;
						if (Contradiction) {
							array<instantiation_tuple> temp(max((size_t) 1, possible_values.length));
							for (const instantiation_tuple& value : possible_values) {
								instantiation_tuple& new_value = temp[temp.length];
								if (!::init(new_value, value)) {
									for (auto& element : temp) core::free(element);
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : possible_values) core::free(element);
									possible_values.clear(); return false;
								} if (!new_value.unify_value(left_most->variable - 1, args[k])) {
									core::free(new_value);
									continue;
								}
								temp.length++;
							}

							if (temp.length != 0) {
								array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
								set_union(union_result, new_possible_values, temp);
								swap(new_possible_values, union_result);
								for (auto& element : temp) core::free(element);
								for (auto& element : union_result) core::free(element);
							}
						} else {
							for (unsigned int i = 0; i < new_possible_values.length; i++) {
								if (!new_possible_values[i].antiunify_value(left_most->variable - 1, args[k])) {
									core::free(new_possible_values[i]);
									new_possible_values.remove(i--);
								}
							}
						}
					}
				}
			}

			for (auto& element : possible_values) core::free(element);
			swap(new_possible_values, possible_values);

		} else if (formula->type == FormulaType::BINARY_APPLICATION) {
			if (formula->ternary.first->constant == (unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL) {
				array<instantiation_tuple> temp(possible_values.length);
				if (formula->ternary.second->type == TermType::VARIABLE) {
					for (unsigned int i = 0; i < possible_values.length; i++) {
						if (!::init(temp[temp.length], possible_values[i])) {
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : possible_values) core::free(element);
							possible_values.clear(); return false;
						} else if (!Contradiction && !temp[temp.length].unify_greater_than_or_equal(formula->ternary.second->variable - 1, formula->ternary.third)) {
							core::free(temp[temp.length]); continue;
						} else if (Contradiction
								&& (!temp[temp.length].unify_less_than_or_equal(formula->ternary.second->variable - 1, formula->ternary.third)
								 || !temp[temp.length].antiunify_value(formula->ternary.second->variable - 1, formula->ternary.third)))
						{
							core::free(temp[temp.length]); continue;
						}
						temp.length++;
					}
				} else if (formula->ternary.third->type == TermType::VARIABLE) {
					for (unsigned int i = 0; i < possible_values.length; i++) {
						if (!::init(temp[temp.length], possible_values[i])) {
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : possible_values) core::free(element);
							possible_values.clear(); return false;
						} else if (!Contradiction && !temp[temp.length].unify_less_than_or_equal(formula->ternary.third->variable - 1, formula->ternary.second)) {
							core::free(temp[temp.length]); continue;
						} else if (Contradiction
								&& (!temp[temp.length].unify_greater_than_or_equal(formula->ternary.third->variable - 1, formula->ternary.second)
								 || !temp[temp.length].antiunify_value(formula->ternary.third->variable - 1, formula->ternary.second)))
						{
							core::free(temp[temp.length]); continue;
						}
						temp.length++;
					}
				} else if (formula->ternary.second->type == TermType::NUMBER && formula->ternary.third->type == TermType::NUMBER) {
					if (!Contradiction && formula->ternary.third->number <= formula->ternary.second->number) {
						for (auto& element : new_possible_values) core::free(element);
						return true;
					} else if (Contradiction && formula->ternary.second->number < formula->ternary.third->number) {
						for (auto& element : new_possible_values) core::free(element);
						return true;
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				} else {
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : possible_values) core::free(element);
					possible_values.clear(); return false;
				}
				if (temp.length != 0) {
					sort(temp);
					array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
					set_union(union_result, new_possible_values, temp);
					for (auto& element : temp) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : possible_values) core::free(element);
					swap(union_result, possible_values);
				} else {
					for (auto& element : possible_values) core::free(element);
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
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				possible_values.clear(); return false;
			}

			unsigned int arg1;
			if (formula->ternary.second->type == TermType::VARIABLE)
				arg1 = 0;
			else if (formula->ternary.second->type == TermType::CONSTANT)
				arg1 = formula->ternary.second->constant;
			else {
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				possible_values.clear(); return false;
			}

			unsigned int arg2;
			if (formula->ternary.third->type == TermType::VARIABLE)
				arg2 = 0;
			else if (formula->ternary.third->type == TermType::CONSTANT)
				arg2 = formula->ternary.third->constant;
			else {
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
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
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
						return false;
					}
					instantiation_tuple values(possible_values[i]);
					if (predicate == 0 && rel.key.predicate != 0) {
						Term* temp = Term::new_constant(rel.key.predicate);
						bool unifies = values.unify_value(formula->ternary.first->variable - 1, temp);
						core::free(*temp); core::free(temp);
						if (!unifies) continue;
					} if (arg1 == 0 && rel.key.arg1 != 0) {
						Term* temp = Term::new_constant(rel.key.arg1);
						bool unifies = values.unify_value(formula->ternary.second->variable - 1, temp);
						core::free(*temp); core::free(temp);
						if (!unifies) continue;
					} if (arg2 == 0 && rel.key.arg2 != 0) {
						Term* temp = Term::new_constant(rel.key.arg2);
						bool unifies = values.unify_value(formula->ternary.third->variable - 1, temp);
						core::free(*temp); core::free(temp);
						if (!unifies) continue;
					}

					if (rel.key.predicate == 0) {
						if (predicate == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.first->variable - 1, term);
								core::free(*term); core::free(term);
								if (!unifies) {
									core::free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(predicate)) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								temp.length++;
							}
						}
					} else if (rel.key.arg1 == 0) {
						if (arg1 == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.second->variable - 1, term);
								core::free(*term); core::free(term);
								if (!unifies) {
									core::free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(arg1)) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								temp.length++;
							}
						}
					} else {
						if (arg2 == 0) {
							for (unsigned int constant : constants) {
								instantiation_tuple& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : temp) core::free(element);
									return false;
								}
								Term* term = Term::new_constant(constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									return false;
								}
								bool unifies = new_values.unify_value(formula->ternary.first->variable - 1, term);
								core::free(*term); core::free(term);
								if (!unifies) {
									core::free(new_values);
									continue;
								}
								temp.length++;
							}
						} else {
							instantiation_tuple& new_values = temp[temp.length];
							if (constants.contains(arg2)) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
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
				for (auto& element : temp) core::free(element);
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				swap(union_result, possible_values);
			} else {
				for (auto& element : possible_values) core::free(element);
				swap(new_possible_values, possible_values);
			}

		} else if (formula->type == FormulaType::AND) {
			if (Contradiction) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					array<instantiation_tuple> copy(possible_values.length);
					for (unsigned int j = 0; j < possible_values.length; j++) {
						if (!::init(copy[copy.length], possible_values[j])) {
							for (auto& element : copy) core::free(element);
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}
						copy.length++;
					}
					if (!is_provable_without_abduction<true>(formula->array.operands[i], quantifiers, copy, prover))
						continue;

					array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
					set_union(union_result, new_possible_values, copy);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : copy) core::free(element);
					swap(union_result, new_possible_values);
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				return possible_values.length != 0;
			}

			for (unsigned int i = 0; i < formula->array.length; i++)
				if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
			if (possible_values.length != 0) {
				array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) core::free(element);
				for (auto& element : new_possible_values) core::free(element);
				swap(possible_values, union_result);
			} else {
				swap(possible_values, new_possible_values);
			}

		} else if (formula->type == FormulaType::OR) {
			if (Contradiction) {
				for (unsigned int i = 0; i < formula->array.length; i++)
					if (!is_provable_without_abduction<true>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
				if (possible_values.length != 0) {
					array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					swap(possible_values, union_result);
				} else {
					swap(possible_values, new_possible_values);
				}
				return possible_values.length != 0;
			}

			for (unsigned int i = 0; i < formula->array.length; i++) {
				array<instantiation_tuple> copy(possible_values.length);
				for (unsigned int j = 0; j < possible_values.length; j++) {
					if (!::init(copy[copy.length], possible_values[j])) {
						for (auto& element : copy) core::free(element);
						for (auto& element : new_possible_values) core::free(element);
						return false;
					}
					copy.length++;
				}
				if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, copy, prover))
					continue;

				array<instantiation_tuple> union_result(new_possible_values.length + copy.length);
				set_union(union_result, new_possible_values, copy);
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : copy) core::free(element);
				swap(union_result, new_possible_values);
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::IF_THEN) {
			if (Contradiction) {
				if (is_provable_without_abduction<false>(formula->binary.left, quantifiers, possible_values, prover))
					is_provable_without_abduction<true>(formula->binary.right, quantifiers, possible_values, prover);
				if (possible_values.length != 0) {
					array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					swap(union_result, possible_values);
				} else {
					swap(possible_values, new_possible_values);
				}
				return possible_values.length != 0;
			}

			array<instantiation_tuple> first(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(first[first.length], possible_values[j])) {
					for (auto& element : first) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					return false;
				}
				first.length++;
			}
			if (!is_provable_without_abduction<true>(formula->binary.left, quantifiers, first, prover))
				first.length = 0;

			array<instantiation_tuple> second(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(second[second.length], possible_values[j])) {
					for (auto& element : first) core::free(element);
					for (auto& element : second) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					return false;
				}
				second.length++;
			}
			if (!is_provable_without_abduction<false>(formula->binary.left, quantifiers, second, prover))
				second.length = 0;

			array<instantiation_tuple> temp(max((size_t) 1, new_possible_values.length + first.length));
			set_union(temp, new_possible_values, first);
			for (auto& element : first) core::free(element);
			for (auto& element : new_possible_values) core::free(element);
			new_possible_values.length = 0;
			if (!new_possible_values.ensure_capacity(temp.length + second.length)) {
				for (auto& element : temp) core::free(element);
				for (auto& element : second) core::free(element);
				return false;
			}
			set_union(new_possible_values, temp, second);
			for (auto& element : temp) core::free(element);
			for (auto& element : second) core::free(element);
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::NOT) {
			if (!is_provable_without_abduction<!Contradiction>(formula->unary.operand, quantifiers, possible_values, prover)) {
				for (auto& element : new_possible_values) core::free(element);
				return false;
			}
			if (possible_values.length != 0) {
				array<instantiation_tuple> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) core::free(element);
				for (auto& element : new_possible_values) core::free(element);
				swap(union_result, possible_values);
			} else {
				swap(possible_values, new_possible_values);
			}
			return possible_values.length != 0;

		} else if (formula->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!exists_is_provable_without_abduction<true>(formula, quantifiers, possible_values, new_possible_values, prover)) {
					for (auto& element : new_possible_values) core::free(element);
					return false;
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				return possible_values.length != 0;
			}

			if (!for_all_is_provable_without_abduction(formula, quantifiers, possible_values, new_possible_values, prover)) {
				for (auto& element : new_possible_values) core::free(element);
				return false;
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::EXISTS) {
			if (Contradiction) {
				/* TODO: implement this */
				//fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				possible_values.clear(); return false;
			}

			if (!exists_is_provable_without_abduction<false>(formula, quantifiers, possible_values, new_possible_values, prover)) {
				for (auto& element : new_possible_values) core::free(element);
				return false;
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::EQUALS) {
			if ((formula->binary.left->type == TermType::UNARY_APPLICATION
			  && formula->binary.left->binary.left->type == TermType::CONSTANT
			  && formula->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE)
			 || (formula->binary.right->type == TermType::UNARY_APPLICATION
			  && formula->binary.right->binary.left->type == TermType::CONSTANT
			  && formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::SIZE))
			{
				Formula* left = formula->binary.left;
				Formula* right = formula->binary.right;
				if (left->type != TermType::UNARY_APPLICATION
				 || left->binary.left->type != TermType::CONSTANT
				 || left->binary.left->constant != (unsigned int) built_in_predicates::SIZE)
				{
					swap(left, right);
				}

				Term* set_definition = left->binary.right;
				if (set_definition->type == TermType::LAMBDA) {
					/* size(^[x]:f(x))=n */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axioms.data == nullptr) continue;

						if (right->type == TermType::NUMBER) {
							if ((!Contradiction && right->number != sets.sets[i].set_size)
							 || (Contradiction && right->number == sets.sets[i].set_size))
								continue;
						}

						array<Formula*> quantifiers(4);
						array_map<Formula*, Term*> unifications(4);
						if (!unify(set_definition, sets.sets[i].size_axioms[0]->formula->binary.left->binary.right, quantifiers, unifications))
							continue;

						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = temp[temp.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								return false;
							}

							if (right->type == TermType::VARIABLE) {
								if ((!Contradiction && !new_values.unify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right))
								 || (Contradiction && !new_values.antiunify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right)))
								{
									core::free(new_values);
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
								core::free(new_values);
								continue;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							insertion_sort(temp, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (set_definition->type == TermType::VARIABLE) {
					/* size(x)=n */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axioms.data == nullptr) continue;
						Formula* set_formula = sets.sets[i].set_formula();
						if (set_formula->type != TermType::UNARY_APPLICATION
						 || set_formula->binary.left->type != TermType::CONSTANT
						 || set_formula->binary.right->type != TermType::VARIABLE
						 || set_formula->binary.right->variable != 1)
							continue;

						if (right->type == TermType::NUMBER) {
							if ((!Contradiction && right->number != sets.sets[i].set_size)
							 || (Contradiction && right->number == sets.sets[i].set_size))
								continue;
						}

						array<instantiation_tuple> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = temp[temp.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								return false;
							}

							if (!new_values.unify_value(set_definition->variable - 1, set_formula->binary.left)) {
								core::free(new_values);
								continue;
							}

							if (right->type == TermType::VARIABLE) {
								if ((!Contradiction && !new_values.unify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right))
								 || (Contradiction && !new_values.antiunify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right)))
								{
									core::free(new_values);
									continue;
								}
							}
							temp.length++;
						}
						if (temp.length != 0) {
							insertion_sort(temp, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
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
					core::free(*set_formula); core::free(set_formula);
					if (!contains) {
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}

					if (right->type == TermType::NUMBER) {
						if ((!Contradiction && right->number != sets.sets[set_id].set_size)
						 || (Contradiction && right->number == sets.sets[set_id].set_size))
						{
							for (auto& element : possible_values) core::free(element);
							possible_values.clear(); return false;
						}
					}
					array<instantiation_tuple> temp(possible_values.length);
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const instantiation_tuple& values = possible_values[i];
						instantiation_tuple& new_values = temp[temp.length];
						if (!::init(new_values, values)) {
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
							return false;
						}
						if (right->type == TermType::VARIABLE) {
							if ((!Contradiction && !new_values.unify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right))
							 || (Contradiction && !new_values.antiunify_value(right->variable - 1, sets.sets[i].size_axioms[0]->formula->binary.right)))
							{
								core::free(new_values);
								for (auto& element : possible_values) core::free(element);
								possible_values.clear(); return false;
							}
						}
						temp.length++;
					}
					if (temp.length != 0) {
						insertion_sort(temp, default_sorter());
						array<instantiation_tuple> union_result(new_possible_values.length + temp.length);
						set_union(union_result, new_possible_values, temp);
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
						swap(union_result, new_possible_values);
					}

				} else {
					fprintf(stderr, "theory.is_provable_without_abduction ERROR: Unsupported set size axiom.\n");
					return false;
				}

				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				return possible_values.length != 0;
			}

			/* A=A */
			array_map<Formula*, Term*> unifications(4);
			if (!Contradiction && unify(formula->binary.left, formula->binary.right, quantifiers, unifications)) {
				if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
					for (auto& element : new_possible_values) core::free(element);
					return false;
				}
				unsigned int old_size = new_possible_values.length;
				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& values = possible_values[i];
					instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
					if (!::init(new_values, values)) {
						for (auto& element : new_possible_values) core::free(element);
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
						core::free(new_values);
						continue;
					}
					new_possible_values.length++;
				}
				if (new_possible_values.length > old_size) {
					insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
					array<instantiation_tuple> union_result(new_possible_values.length);
					set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
					for (auto& element : new_possible_values) core::free(element);
					swap(union_result, new_possible_values);
				}
			}

			/* a=b */
			if (Contradiction) {
				if (formula->binary.left->type == TermType::CONSTANT) {
					if (formula->binary.right->type == TermType::CONSTANT && formula->binary.right->constant == formula->binary.left->constant) {
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								return false;
							}

							if (!new_values.antiunify_value(formula->binary.right->variable - 1, formula->binary.left)) {
								core::free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (formula->binary.left->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) core::free(element);
						return false;
					}
					unsigned int old_size = new_possible_values.length;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const instantiation_tuple& values = possible_values[i];
						instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
						if (!::init(new_values, values)) {
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}

						if (!new_values.antiunify_value(formula->binary.left->variable - 1, formula->binary.right)) {
							core::free(new_values);
							continue;
						}
						new_possible_values.length++;
					}
					if (new_possible_values.length > old_size) {
						insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
						array<instantiation_tuple> union_result(new_possible_values.length);
						set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
						for (auto& element : new_possible_values) core::free(element);
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
								for (auto& element : new_possible_values) core::free(element);
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const instantiation_tuple& values = possible_values[i];
								instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
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
									core::free(new_values);
									continue;
								}
								new_possible_values.length++;
							}
							if (new_possible_values.length > old_size) {
								insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
								array<instantiation_tuple> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
								for (auto& element : new_possible_values) core::free(element);
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
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const instantiation_tuple& values = possible_values[i];
							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
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
								core::free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
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
							for (auto& element : new_possible_values) core::free(element);
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (Proof* definition : c.definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;

							instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
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
								core::free(new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
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
							for (auto& element : new_possible_values) core::free(element);
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
							if (!::init(new_new_values, new_values)) {
								for (auto& element : new_possible_values) core::free(element);
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
								core::free(new_new_values);
								continue;
							}
							new_possible_values.length++;
						}
						if (new_possible_values.length > old_size) {
							insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
							array<instantiation_tuple> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
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
								for (auto& element : new_possible_values) core::free(element);
								return true;
							}
						} else {
							if (Contradiction) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : possible_values) core::free(element);
								possible_values.clear(); return false;
							}
							if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
								for (auto& element : new_possible_values) core::free(element);
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const instantiation_tuple& values = possible_values[i];
								instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
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
									core::free(new_values);
									continue;
								}
								new_possible_values.length++;
							}
							if (new_possible_values.length > old_size) {
								insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
								array<instantiation_tuple> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
								for (auto& element : new_possible_values) core::free(element);
								swap(union_result, new_possible_values);
							}
						}
					}
				} else if (left->binary.right->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) core::free(element);
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
										if (!::init(new_values, values)) {
											for (auto& element : new_possible_values) core::free(element);
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									if (Contradiction) {
										fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
										for (auto& element : new_possible_values) core::free(element);
										for (auto& element : possible_values) core::free(element);
										possible_values.clear(); return false;
									}
									instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
									if (!::init(new_values, values)) {
										for (auto& element : new_possible_values) core::free(element);
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
										core::free(new_values);
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
									for (auto& element : new_possible_values) core::free(element);
									return false;
								}
								if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
									if (Contradiction) {
										instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
										if (!::init(new_values, values)) {
											for (auto& element : new_possible_values) core::free(element);
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									if (Contradiction) {
										fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
										for (auto& element : new_possible_values) core::free(element);
										for (auto& element : possible_values) core::free(element);
										possible_values.clear(); return false;
									}
									instantiation_tuple& new_values = new_possible_values[new_possible_values.length];
									if (!::init(new_values, values)) {
										for (auto& element : new_possible_values) core::free(element);
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
										core::free(new_values);
										continue;
									}
									new_possible_values.length++;
								}
							}
						}
					}
					if (new_possible_values.length > old_size) {
						insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
						array<instantiation_tuple> union_result(new_possible_values.length);
						set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
						for (auto& element : new_possible_values) core::free(element);
						swap(union_result, new_possible_values);
					}
				}
			}

			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::TRUE) {
			if (Contradiction) {
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				return false;
			} else {
				for (auto& element : new_possible_values) core::free(element);
				return true;
			}
		} else if (formula->type == FormulaType::FALSE) {
			if (Contradiction) {
				for (auto& element : new_possible_values) core::free(element);
				return true;
			} else {
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				return false;
			}
		} else {
			fprintf(stderr, "theory.is_provable_without_abduction ERROR: Unsupported FormulaType.\n");
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
		case FormulaType::NUMBER:
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
				core::free(unifications[unifications.length]);
			}
			return true;
		case FormulaType::NOT:
			if (!unifications.ensure_capacity(unifications.length + 2)
			 || !array_map_init(unifications[unifications.length], 2))
				return false;
			if (unify_atom(formula, atom, quantifiers, unifications[unifications.length]))  {
				unifications.length++;
			} else {
				core::free(unifications[unifications.length]);
			}

			if (!array_map_init(unifications[unifications.length], 2))
				return false;
			if (unify_atom(formula->unary.operand, atom, quantifiers, unifications[unifications.length]))  {
				unifications.length++;
			} else {
				core::free(unifications[unifications.length]);
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

	template<bool ResolveInconsistencies, bool SubtractProvableElements>
	bool check_set_membership(unsigned int set_id,
			array<instantiation_tuple>& possible_values,
			array_map<tuple, array<unsigned int>>& new_elements,
			array<Formula*>& quantifiers,
			array<Formula*>& new_atoms,
			unsigned int new_atom_index)
	{
		/* check if any of the possible values of the quantified variables are newly provable members of this set */
		hash_set<tuple> provable_elements(16);
		set_membership_prover prover(sets, new_elements);
		if (!prover.get_provable_elements(set_id, provable_elements)) {
			for (tuple& tup : provable_elements) core::free(tup);
			return false;
		}
		unsigned int provable_set_size = provable_elements.size;
		for (unsigned int j = 0; SubtractProvableElements && j < possible_values.length; j++) {
			const instantiation_tuple& values = possible_values[j];
			tuple new_tuple;
			if (!::init(new_tuple, values.length)) {
				for (tuple& tup : provable_elements) core::free(tup);
				return false;
			}
			bool has_any = false;
			for (unsigned int k = 0; k < values.length; k++) {
				if (values.values[k].type == instantiation_type::CONSTANT) {
					new_tuple[k].type = tuple_element_type::CONSTANT;
					new_tuple[k].constant = values.values[k].constant;
				} else if (values.values[k].type == instantiation_type::NUMBER) {
					new_tuple[k].type = tuple_element_type::NUMBER;
					new_tuple[k].number = values.values[k].number;
				} else if (values.values[k].type == instantiation_type::STRING) {
					new_tuple[k].type = tuple_element_type::STRING;
					if (!core::init(new_tuple[k].str, values.values[k].str)) {
						for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
						for (tuple& tup : provable_elements) core::free(tup);
						core::free(new_tuple.elements); return false;
					}
				} else if (values.values[k].type == instantiation_type::ANY || values.values[k].type == instantiation_type::ANY_NUMBER) {
					for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
					core::free(new_tuple.elements);
					has_any = true;
					break;
				}
			}
			if (has_any)
				continue;

			bool is_old = provable_elements.contains(new_tuple);
			core::free(new_tuple);
			if (is_old) {
				core::free(possible_values[j]);
				move(possible_values[possible_values.length - 1], possible_values[j]);
				possible_values.length--;
				j--;
			}
		}
		for (tuple& tup : provable_elements) core::free(tup);

		if (possible_values.length == 0)
			return true;

		/* all values in `possible_values` are newly provable */
		for (const auto& entry : sets.extensional_graph.vertices[set_id].parents) {
			unsigned int other_set = entry.key;
			Formula* other_formula = sets.sets[other_set].set_formula();

			array<instantiation_tuple> copy(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(copy[copy.length], possible_values[j]))
					return false;
				copy.length++;
			}

			if (is_provable_without_abduction<true>(other_formula, quantifiers, copy, prover)) {
				/* we found a contradiction with a universally-quantified theorem */
				for (auto& element : copy) core::free(element);
				return false;
			}
			for (auto& element : copy) core::free(element);
		}

		/* make sure the set is big enough to contain the new elements */
		if (provable_set_size > sets.sets[set_id].set_size
		 && (!ResolveInconsistencies || !sets.template set_size_axiom<true>(set_id, provable_set_size)))
		{
			return false;
		}

		/* add the newly provable elements to `new_elements` */
		if (!new_elements.ensure_capacity(new_elements.size + possible_values.length))
			return false;
		for (unsigned int j = 0; j < possible_values.length; j++) {
			tuple& new_tuple = *((tuple*) alloca(sizeof(tuple)));
			if (!::init(new_tuple, possible_values[0].length))
				return false;

			const instantiation_tuple& values = possible_values[j];
			bool has_any = false;
			for (unsigned int k = 0; k < values.length; k++) {
				if (values.values[k].type == instantiation_type::CONSTANT) {
					new_tuple[k].type = tuple_element_type::CONSTANT;
					new_tuple[k].constant = values.values[k].constant;
				} else if (values.values[k].type == instantiation_type::NUMBER) {
					new_tuple[k].type = tuple_element_type::NUMBER;
					new_tuple[k].number = values.values[k].number;
				} else if (values.values[k].type == instantiation_type::STRING) {
					new_tuple[k].type = tuple_element_type::STRING;
					if (!core::init(new_tuple[k].str, values.values[k].str)) {
						for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
						core::free(new_tuple.elements); return false;
					}
				} else if (values.values[k].type == instantiation_type::ANY || values.values[k].type == instantiation_type::ANY_NUMBER) {
					for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
					core::free(new_tuple.elements);
					has_any = true;
					break;
				}
			}
			if (has_any) continue;

			unsigned int index = new_elements.index_of(new_tuple);
			if (index == new_elements.size) {
				if (!array_init(new_elements.values[index], 4)) {
					core::free(new_tuple);
					return false;
				}
				move(new_tuple, new_elements.keys[index]);
				new_elements.size++;
			} else {
				core::free(new_tuple);
			}
			tuple& tup = new_elements.keys[index];
			bool is_element_of_descendant = false;
			for (unsigned int other_set : new_elements.values[index]) {
				if (sets.sets[set_id].descendants.contains(other_set)) {
					is_element_of_descendant = true;
					break;
				}
			}
			if (is_element_of_descendant) {
				continue;
			} else {
				if (!new_elements.values[index].add(set_id))
					return false;
				for (unsigned int l = 0; l < new_elements.values[index].length; l++)
					if (new_elements.values[index][l] != set_id && sets.sets[set_id].descendants.contains(new_elements.values[index][l]))
						new_elements.values[index].remove(l--);
			}

			/* get the ancestors of this set that don't already provably contain `tup` */
			hash_set<unsigned int> visited(16);
			array<unsigned int> stack(8);
			array<Formula*> new_ancestor_terms(8);
			stack[stack.length++] = set_id;
			visited.add(set_id);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				Formula* set_formula = sets.sets[current].set_formula();
				if (set_formula->type == FormulaType::AND) {
					if (!new_ancestor_terms.append(set_formula->array.operands, set_formula->array.length))
						return false;
				} else if (!new_ancestor_terms.add(set_formula)) {
					return false;
				}

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent))
						return false;
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key))
						return false;
				}
			}

			Term** src_variables = (Term**) alloca(sizeof(Term*) * tup.length * 2);
			for (unsigned int i = 0; i < tup.length; i++) {
				src_variables[i] = Term::new_variable(i + 1);
				if (src_variables[i] == nullptr) {
					for (unsigned int j = 0; j < i; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					return false;
				}
			}
			Term** dst_constants = src_variables + tup.length;
			for (unsigned int i = 0; i < tup.length; i++) {
				if (tup[i].type == tuple_element_type::CONSTANT)
					dst_constants[i] = Term::new_constant(tup[i].constant);
				else if (tup[i].type == tuple_element_type::NUMBER)
					dst_constants[i] = Term::new_number(tup[i].number);
				else if (tup[i].type == tuple_element_type::STRING)
					dst_constants[i] = Term::new_string(tup[i].str);
				if (dst_constants[i] == nullptr) {
					for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < i; j++) { core::free(*dst_constants[j]); core::free(dst_constants[j]); }
					return false;
				}
			}

			if (!new_atoms.ensure_capacity(new_atoms.length + new_ancestor_terms.length)) {
				for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
				for (unsigned int j = 0; j < tup.length; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
				return false;
			}
			for (Formula* new_descendant_term : new_ancestor_terms) {
				Formula* substituted_atom = substitute_all(new_descendant_term, src_variables, dst_constants, tup.length);
				if (substituted_atom == nullptr) {
					for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < tup.length; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
					return false;
				}
				bool atom_already_exists = false;
				for (unsigned int k = new_atom_index; !atom_already_exists && k < new_atoms.length; k++)
					if (*new_atoms[k] == *substituted_atom) atom_already_exists = true;
				if (atom_already_exists) {
					core::free(*substituted_atom); if (substituted_atom->reference_count == 0) core::free(substituted_atom);
				} else {
					new_atoms[new_atoms.length++] = substituted_atom;
				}
			}

			for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
			for (unsigned int j = 0; j < tup.length; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_set_membership_after_addition(array<Formula*>& new_atoms, Args&&... visitor) {
		/* We first need to find all pairs of tuples and sets such that, with
		   the addition of `new_atom` as an axiom, the tuple must necessarily
		   belong to the set (and is not the case otherwise without
		   `new_atom`). To do this, find all sets that contain an atom that
		   unifies with `new_atom`. */
		unsigned int new_atom_index = 0;
		array_map<tuple, array<unsigned int>> new_elements(4);
		while (new_atom_index < new_atoms.length) {
			Formula* next_new_atom = new_atoms[new_atom_index++];

			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axioms.data == nullptr)
					continue;
				Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}
				array<array_map<Formula*, Term*>> unifications(4);
				if (!get_unifying_atoms(set_formula, next_new_atom, quantifiers, unifications)) {
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}

				array<instantiation_tuple> possible_values(max((size_t) 1, 4 * unifications.length));
				for (const array_map<Formula*, Term*>& unification : unifications) {
					/* given `unification`, find the values of the quantified variables that make `set_formula` necessarily true */
					bool valid_unification = true;
					instantiation_tuple& values = possible_values[possible_values.length];
					if (!::init(values, sets.sets[i].arity)) {
						for (auto& element : unifications) core::free(element);
						for (auto& element : possible_values) core::free(element);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						return false;
					}
					for (const auto& entry : unification) {
						if (entry.key->quantifier.variable <= sets.sets[i].arity) {
							if (entry.value->type == TermType::CONSTANT) {
								core::free(values.values[entry.key->quantifier.variable - 1]);
								values.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
								values.values[entry.key->quantifier.variable - 1].constant = entry.value->constant;
							} else if (entry.value->type == TermType::NUMBER) {
								core::free(values.values[entry.key->quantifier.variable - 1]);
								values.values[entry.key->quantifier.variable - 1].type = instantiation_type::NUMBER;
								values.values[entry.key->quantifier.variable - 1].number = entry.value->number;
							} else if (entry.value->type == TermType::STRING) {
								core::free(values.values[entry.key->quantifier.variable - 1]);
								values.values[entry.key->quantifier.variable - 1].type = instantiation_type::STRING;
								if (!core::init(values.values[entry.key->quantifier.variable - 1].str, entry.value->str)) {
									values.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
									for (auto& element : unifications) core::free(element);
									for (auto& element : possible_values) core::free(element);
									for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
									core::free(values); return false;
								}
							} else if (entry.value->type != TermType::VARIABLE) {
								valid_unification = false;
								break;
							}
						}
					}
					if (!valid_unification) {
						core::free(values);
						continue;
					}
					possible_values.length++;
				}
				for (auto& element : unifications) core::free(element);
				if (possible_values.length == 0) continue;

				/* remove possible_values that are subsets of another possible_value */
				for (unsigned int i = 0; i < possible_values.length; i++) {
					const instantiation_tuple& first = possible_values[i];
					for (unsigned int j = 0; j < possible_values.length; j++) {
						if (j == i) continue;
						const instantiation_tuple& second = possible_values[j];

						/* check if `first` is a subset of `second` */
						/* NOTE: this is not a general `subset` function for `instantiation_tuple` */
						bool is_subset = true;
						for (uint_fast8_t k = 0; k < first.length; k++) {
							if (second.values[k].type == instantiation_type::ANY) {
								continue;
							} else if (first.values[k].type == instantiation_type::ANY) {
								is_subset = false;
								break;
							} else if (first.values[k].type == instantiation_type::CONSTANT) {
								if (second.values[k].type != instantiation_type::CONSTANT || first.values[k].constant != second.values[k].constant) {
									is_subset = false;
									break;
								}
							} else if (first.values[k].type == instantiation_type::NUMBER) {
								if (second.values[k].type != instantiation_type::NUMBER || first.values[k].number != second.values[k].number) {
									is_subset = false;
									break;
								}
							} else if (first.values[k].type == instantiation_type::STRING) {
								if (second.values[k].type != instantiation_type::STRING || first.values[k].str != second.values[k].str) {
									is_subset = false;
									break;
								}
							}
						}
						if (is_subset) {
							core::free(possible_values[i]);
							possible_values.remove(i--);
							break;
						}
					}
				}

				set_membership_prover prover(sets, new_elements);
				array<instantiation_tuple> new_possible_values(max((size_t) 1, 4 * possible_values.length));
				for (const instantiation_tuple& possible_value : possible_values) {
					array<instantiation_tuple> temp_possible_values(4);
					instantiation_tuple& values = temp_possible_values[0];
					if (!::init(values, possible_value)) {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						return false;
					}
					temp_possible_values.length++;

					if (!is_provable_without_abduction<false>(set_formula, quantifiers, temp_possible_values, prover)) {
						for (auto& element : temp_possible_values) core::free(element);
						continue;
					}

					array<instantiation_tuple> union_result(new_possible_values.length + temp_possible_values.length);
					set_union(union_result, new_possible_values, temp_possible_values);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : temp_possible_values) core::free(element);
					swap(union_result, new_possible_values);
				}
				for (auto& element : possible_values) core::free(element);

				/* the addition could cause this formula to be provably false */
				hash_set<tuple> provable_elements(16);
				if (!sets.get_provable_elements(i, provable_elements)) {
					for (auto& element : new_possible_values) core::free(element);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}
				array<instantiation_tuple> temp_possible_values(max(1, provable_elements.size));
				for (const tuple& element : provable_elements) {
					instantiation_tuple& values = temp_possible_values[temp_possible_values.length];
					if (!::init(values, element.length)) {
						for (auto& element : temp_possible_values) core::free(element);
						for (tuple& tup : provable_elements) core::free(tup);
						for (auto& element : new_possible_values) core::free(element);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						return false;
					}
					for (uint_fast8_t k = 0; k < element.length; k++) {
						if (element[k].type == tuple_element_type::CONSTANT) {
							core::free(values.values[k]);
							values.values[k].type = instantiation_type::CONSTANT;
							values.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							core::free(values.values[k]);
							values.values[k].type = instantiation_type::NUMBER;
							values.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							core::free(values.values[k]);
							values.values[k].type = instantiation_type::STRING;
							if (!core::init(values.values[k].str, element[k].str)) {
								values.values[k].type = instantiation_type::CONSTANT;
								for (auto& element : temp_possible_values) core::free(element);
								for (tuple& tup : provable_elements) core::free(tup);
								for (auto& element : new_possible_values) core::free(element);
								for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
								values.values[k].type = instantiation_type::CONSTANT;
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;
				}
				for (tuple& tup : provable_elements) core::free(tup);
				if (temp_possible_values.length != 0 && is_provable_without_abduction<true>(set_formula, quantifiers, temp_possible_values, prover)) {
					/* we found a contradiction */
					for (auto& element : temp_possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}

				if (new_possible_values.length == 0) continue;
				if (!check_set_membership<ResolveInconsistencies, true>(i, new_possible_values, new_elements, quantifiers, new_atoms, new_atom_index + 1)) {
					for (auto& element : new_possible_values) core::free(element);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}
				for (auto& element : new_possible_values) core::free(element);
			}
		}

		/* get the old sets that contain the elements in `new_elements` */
		array<unsigned int>* old_elements = (array<unsigned int>*) malloc(max((size_t) 1, sizeof(array<unsigned int>) * new_elements.size));
		if (old_elements == nullptr) {
			fprintf(stderr, "check_set_membership_after_addition ERROR: Out of memory.\n");
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			return false;
		}
		for (unsigned int i = 0; i < new_elements.size; i++) {
			if (!array_init(old_elements[i], 4)) {
				for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
				for (unsigned int j = 0; j < i; j++) core::free(old_elements[j]);
				core::free(old_elements); return false;
			}
		}
		for (unsigned int i = 0; i < new_elements.size; i++) {
			array<unsigned int> stack(max((size_t) 8, new_elements.values[i].length));
			hash_set<unsigned int> visited(stack.capacity * RESIZE_THRESHOLD_INVERSE);
			for (unsigned int new_set : new_elements.values[i]) {
				stack[stack.length++] = new_set;
				visited.add(new_set);
			}
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				if (sets.sets[current].index_of_element(new_elements.keys[i]) < sets.sets[current].element_count()) {
					if (!old_elements[i].add(current)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
					continue;
				}

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				}
			}
		}

		/* add the elements to the sets */
		for (unsigned int i = 0; i < new_elements.size; i++) {
			const tuple& element = new_elements.keys[i];
			for (unsigned int old_set : old_elements[i])
				sets.sets[old_set].remove_element_at(sets.sets[old_set].index_of_element(element));
			for (unsigned int new_set : new_elements.values[i]) {
				if (!sets.sets[new_set].add_element(element)) {
					for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					core::free(old_elements); return false;
				}
			}
		}

		/* make sure there are no inconsistencies */
		for (unsigned int i = 0; i < new_elements.size; i++) {
			/* get the old ancestors */
			array<unsigned int> stack(max((size_t) 8, old_elements[i].length));
			hash_set<unsigned int> old_ancestors(max((size_t) 16, old_elements[i].length));
			for (unsigned int old_set : old_elements[i]) {
				stack[stack.length++] = old_set;
				old_ancestors.add(old_set);
			}
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(parent)) continue;
					if (!old_ancestors.add(parent) || !stack.add(parent)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(entry.key)) continue;
					if (!old_ancestors.add(entry.key) || !stack.add(entry.key)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				}
			}

			/* get the new ancestors */
			hash_set<unsigned int> new_ancestors(max((size_t) 16, new_elements.values[i].length));
			if (!stack.ensure_capacity(new_elements.values[i].length)) {
				for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
				for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
				core::free(old_elements); return false;
			}
			for (unsigned int new_set : new_elements.values[i]) {
				stack[stack.length++] = new_set;
				new_ancestors.add(new_set);
			}
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(parent) || new_ancestors.contains(parent)) continue;
					if (!new_ancestors.add(parent) || !stack.add(parent)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(entry.key) || new_ancestors.contains(entry.key)) continue;
					if (!new_ancestors.add(entry.key) || !stack.add(entry.key)) {
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				}
			}

			/* check the new ancestors of `new_set_id` to make sure the size of each
			   set is at least as large as the number of provable elements of each set */
			for (unsigned int new_ancestor : new_ancestors) {
				hash_set<tuple> provable_elements(16);
				if (!sets.get_provable_elements(new_ancestor, provable_elements)) {
					for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					core::free(old_elements); return false;
				}

				if (sets.sets[new_ancestor].set_size >= provable_elements.size) {
					for (tuple& tup : provable_elements) core::free(tup);
					continue;
				}

				while (true) {
					/* compute the upper bound on the size of this set; if the new size
						violates this bound, change the sizes of other sets to increase the bound */
					array<unsigned int> stack(8); bool graph_changed;
					if (!ResolveInconsistencies || !sets.increase_set_size(new_ancestor, provable_elements.size, stack, graph_changed)) {
						for (tuple& tup : provable_elements) core::free(tup);
						/* undo the set changes */
						for (unsigned int i = 0; i < new_elements.size; i++) {
							const tuple& element = new_elements.keys[i];
							for (unsigned int new_set : new_elements.values[i])
								sets.sets[new_set].remove_element_at(sets.sets[new_set].index_of_element(element));
							for (unsigned int old_set : old_elements[i])
								sets.sets[old_set].add_element(element);
						}
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
					if (sets.sets[new_ancestor].set_size >= provable_elements.size)
						break;
					if (!graph_changed) {
						for (tuple& tup : provable_elements) core::free(tup);
						/* undo the set changes */
						for (unsigned int i = 0; i < new_elements.size; i++) {
							const tuple& element = new_elements.keys[i];
							for (unsigned int new_set : new_elements.values[i])
								sets.sets[new_set].remove_element_at(sets.sets[new_set].index_of_element(element));
							for (unsigned int old_set : old_elements[i])
								sets.sets[old_set].add_element(element);
						}
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				}
				for (tuple& tup : provable_elements) core::free(tup);
			}
		}
		for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
		core::free(old_elements);
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline bool check_set_membership_after_addition(Formula* new_atom, Args&&... visitor) {
		array<Formula*> new_atoms(4);
		new_atoms[new_atoms.length++] = new_atom;
		new_atom->reference_count++;
		bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, std::forward<Args>(visitor)...);
		free_all(new_atoms);
		return result;
	}

	bool check_set_membership_of_element(unsigned int element_set_id,
			unsigned int formula_set_id, Formula* set_formula, unsigned int element_index,
			array<Formula*>& quantifiers, array<Formula*>& old_atoms, unsigned int old_atom_index)
	{
		tuple& element = *((tuple*) alloca(sizeof(tuple)));
		element.elements = (tuple_element*) alloca(sizeof(tuple_element) * sets.sets[element_set_id].arity);
		element.length = sets.sets[element_set_id].arity;
		const tuple_element* element_src = sets.sets[element_set_id].elements.data + (sets.sets[element_set_id].arity * element_index);
		for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) {
			if (!::init(element.elements[k], element_src[k])) {
				for (uint_fast8_t l = 0; l < k; l++) core::free(element[l]);
				return false;
			}
		}

		array<Formula*> conjuncts(4);
		if (set_formula->type == FormulaType::AND) {
			if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
				for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
				return false;
			}
		} else {
			conjuncts[conjuncts.length++] = set_formula;
		}

		sets.sets[element_set_id].remove_element_at(element_index);
		for (unsigned int l = 0; l < conjuncts.length; l++) {
			array<instantiation_tuple> temp_possible_values(1);
			instantiation_tuple& values = temp_possible_values[0];
			if (!::init(values, sets.sets[element_set_id].arity)) {
				for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++)
					move(element[k], sets.sets[element_set_id].elements[sets.sets[element_set_id].elements.length++]);
				return false;
			}
			for (unsigned int k = 0; k < values.length; k++) {
				core::free(values.values[k]);
				if (element[k].type == tuple_element_type::CONSTANT) {
					values.values[k].type = instantiation_type::CONSTANT;
					values.values[k].constant = element[k].constant;
				} else if (element[k].type == tuple_element_type::NUMBER) {
					values.values[k].type = instantiation_type::NUMBER;
					values.values[k].number = element[k].number;
				} else if (element[k].type == tuple_element_type::STRING) {
					values.values[k].type = instantiation_type::STRING;
					if (!core::init(values.values[k].str, element[k].str)) {
						values.values[k].type = instantiation_type::CONSTANT;
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++)
							move(element[k], sets.sets[element_set_id].elements[sets.sets[element_set_id].elements.length++]);
						core::free(values); return false;
					}
				}
			}
			temp_possible_values.length++;

			default_prover prover(sets);
			if (is_provable_without_abduction<false>(conjuncts[l], quantifiers, temp_possible_values, prover)) {
				for (auto& element : temp_possible_values) core::free(element);
				conjuncts.remove(l--);
			}
		}

		if (conjuncts.length == 0) {
			for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++)
				move(element[k], sets.sets[element_set_id].elements[sets.sets[element_set_id].elements.length++]);
		} else {
			/* this removed element could belong to an ancestor */
			/* first get the sets that currently contain `element` */
			array<unsigned int> stack(8);
			array<unsigned int> old_sets(8);
			hash_set<unsigned int> visited(16);
			stack[stack.length++] = 1;
			visited.add(1);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				if (sets.sets[current].index_of_element(element) < sets.sets[current].element_count()) {
					if (!old_sets.add(current)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
					continue;
				}
				bool has_containing_descendant = false;
				for (unsigned int old_set : old_sets) {
					if (sets.sets[current].descendants.contains(old_set)) {
						has_containing_descendant = true;
						break;
					}
				}
				if (has_containing_descendant) continue;

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				}
			}
			visited.clear();

			array<unsigned int> new_sets(8);
			stack[stack.length++] = element_set_id;
			visited.add(element_set_id);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				bool is_ancestor_of_set_containing_element = false;
				for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < new_sets.length; j++)
					if (sets.sets[current].descendants.contains(new_sets[j])) is_ancestor_of_set_containing_element = true;
				for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < old_sets.length; j++)
					if (sets.sets[current].descendants.contains(old_sets[j])) is_ancestor_of_set_containing_element = true;
				if (is_ancestor_of_set_containing_element) continue;

				Formula* current_set_formula = sets.sets[current].set_formula();
				if (current_set_formula != set_formula && sets.sets[current].arity == sets.sets[element_set_id].arity) {
					array<instantiation_tuple> temp_possible_values(1);
					instantiation_tuple& values = temp_possible_values[0];
					if (!::init(values, sets.sets[element_set_id].arity)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
					for (unsigned int k = 0; k < values.length; k++) {
						core::free(values.values[k]);
						if (element[k].type == tuple_element_type::CONSTANT) {
							values.values[k].type = instantiation_type::CONSTANT;
							values.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							values.values[k].type = instantiation_type::NUMBER;
							values.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							values.values[k].type = instantiation_type::STRING;
							if (!core::init(values.values[k].str, element[k].str)) {
								values.values[k].type = instantiation_type::CONSTANT;
								for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;

					default_prover prover(sets);
					if (is_provable_without_abduction<false>(current_set_formula, quantifiers, temp_possible_values, prover)) {
						for (auto& element : temp_possible_values) core::free(element);
						for (unsigned int j = 0; j < new_sets.length; j++) {
							if (sets.sets[new_sets[j]].descendants.contains(current))
								new_sets.remove(j--);
						}
						if (!new_sets.add(current)) {
							for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
							return false;
						}
						continue;
					}
				}

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				}
			}

			for (unsigned int new_set : new_sets) {
				if (!sets.sets[new_set].add_element(element)) {
					for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
					return false;
				}
			}

			/* terms in the set formulae of the ancestors of set `formula_set_id` could also be no longer provable */
			stack[stack.length++] = formula_set_id;
			visited.clear(); visited.add(formula_set_id);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				bool current_has_descendant_in_new_set = false;
				for (unsigned int j = 0; !current_has_descendant_in_new_set && j < new_sets.length; j++)
					current_has_descendant_in_new_set = sets.sets[current].descendants.contains(new_sets[j]);
				if (current_has_descendant_in_new_set) continue;

				if (current != formula_set_id) {
					Formula* set_formula = sets.sets[current].set_formula();
					if (set_formula->type == FormulaType::AND) {
						if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
							for (uint_fast8_t k = 0; k < sets.sets[formula_set_id].arity; k++) core::free(element[k]);
							return false;
						}
					} else if (!conjuncts.add(set_formula)) {
						for (uint_fast8_t k = 0; k < sets.sets[formula_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				}

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent)) {
						for (uint_fast8_t k = 0; k < sets.sets[formula_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key)) {
						for (uint_fast8_t k = 0; k < sets.sets[formula_set_id].arity; k++) core::free(element[k]);
						return false;
					}
				}
			}

			Term** src_variables = (Term**) alloca(sizeof(Term*) * sets.sets[formula_set_id].arity * 2);
			for (unsigned int k = 0; k < sets.sets[formula_set_id].arity; k++) {
				src_variables[k] = Term::new_variable(k + 1);
				if (src_variables[k] == nullptr) {
					for (unsigned int j = 0; j < k; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
					return false;
				}
			}
			Term** dst_constants = src_variables + sets.sets[formula_set_id].arity;
			for (unsigned int k = 0; k < sets.sets[formula_set_id].arity; k++) {
				if (element[k].type == tuple_element_type::CONSTANT)
					dst_constants[k] = Term::new_constant(element[k].constant);
				else if (element[k].type == tuple_element_type::NUMBER)
					dst_constants[k] = Term::new_number(element[k].number);
				else if (element[k].type == tuple_element_type::STRING)
					dst_constants[k] = Term::new_string(element[k].str);
				if (dst_constants[k] == nullptr) {
					for (unsigned int j = 0; j < sets.sets[formula_set_id].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < k; j++) { core::free(*dst_constants[j]); core::free(dst_constants[j]); }
					for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
					return false;
				}
			}
			for (Formula* conjunct : conjuncts) {
				Formula* substituted_conjunct = substitute_all(conjunct, src_variables, dst_constants, sets.sets[formula_set_id].arity);
				if (substituted_conjunct == nullptr
				 || !old_atoms.ensure_capacity(old_atoms.length + 1))
				{
					if (substituted_conjunct != nullptr) { core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct); }
					for (unsigned int j = 0; j < sets.sets[formula_set_id].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < sets.sets[formula_set_id].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
					for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
					return false;
				}
				/* make sure `substituted_conjunct` doesn't already exist in `old_atoms` to avoid redundant computation */
				bool in_old_atoms = false;
				for (unsigned int j = old_atom_index + 1; !in_old_atoms && j < old_atoms.length; j++)
					if (*old_atoms[j] == *substituted_conjunct) in_old_atoms = true;
				if (in_old_atoms) {
					core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct);
					continue;
				}
				old_atoms[old_atoms.length++] = substituted_conjunct;
			}
			for (unsigned int j = 0; j < sets.sets[formula_set_id].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
			for (unsigned int j = 0; j < sets.sets[formula_set_id].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
			for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
		}
		return true;
	}

	template<typename... Args>
	bool check_set_membership_after_subtraction(array<Formula*>& old_atoms, Args&&... visitor)
	{
		unsigned int old_atom_index = 0;
		array_map<unsigned int, tuple> old_elements(8);
		while (old_atom_index < old_atoms.length) {
			Formula* next_old_atom = old_atoms[old_atom_index++];

			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axioms.data == nullptr)
					continue;
				Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}
				array<array_map<Formula*, Term*>> unifications(4);
				if (!get_unifying_atoms(set_formula, next_old_atom, quantifiers, unifications)) {
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}

				for (unsigned int j = sets.sets[i].element_count(); j > 0; j--) {
					tuple& element = *((tuple*) alloca(sizeof(tuple)));
					element.elements = (tuple_element*) malloc(sizeof(tuple_element) * sets.sets[i].arity);
					if (element.elements == nullptr) {
						fprintf(stderr, "theory.check_set_membership_after_subtraction ERROR: Out of memory.\n");
						for (auto entry : old_elements) core::free(entry.value);
						return false;
					}
					element.length = sets.sets[i].arity;
					const tuple_element* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * (j - 1));
					for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) {
						if (!::init(element.elements[k], element_src[k])) {
							for (auto entry : old_elements) core::free(entry.value);
							for (uint_fast8_t l = 0; l < k; l++) core::free(element[l]);
							core::free(element.elements); return false;
						}
					}
					bool has_satisfying_unification = false;
					for (const array_map<Formula*, Term*>& unification : unifications) {
						bool matches = true;
						for (const auto& entry : unification) {
							if (entry.key->quantifier.variable > sets.sets[i].arity)
								continue;
							if (element[entry.key->quantifier.variable - 1].type == tuple_element_type::CONSTANT) {
								if (entry.value->type != TermType::CONSTANT
								 || element[entry.key->quantifier.variable - 1].constant != entry.value->constant)
								{
									matches = false;
									break;
								}
							} else if (element[entry.key->quantifier.variable - 1].type == tuple_element_type::NUMBER) {
								if (entry.value->type != TermType::NUMBER
								 || element[entry.key->quantifier.variable - 1].number != entry.value->number)
								{
									matches = false;
									break;
								}
							} else if (element[entry.key->quantifier.variable - 1].type == tuple_element_type::STRING) {
								if (entry.value->type != TermType::STRING
								 || element[entry.key->quantifier.variable - 1].str != entry.value->str)
								{
									matches = false;
									break;
								}
							}
						}
						if (matches) {
							has_satisfying_unification = true;
							break;
						}
					}
					if (!has_satisfying_unification) {
						core::free(element);
						continue;
					}

					array<Formula*> conjuncts(4);
					if (set_formula->type == FormulaType::AND) {
						if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
							for (auto entry : old_elements) core::free(entry.value);
							core::free(element); return false;
						}
					} else {
						conjuncts[conjuncts.length++] = set_formula;
					}

					sets.sets[i].remove_element_at(j - 1);
					if (!old_elements.ensure_capacity(old_elements.size + 1)) {
						for (auto entry : old_elements) core::free(entry.value);
						core::free(element); return false;
					}
					old_elements.keys[old_elements.size] = i;
					move(element, old_elements.values[old_elements.size]);
					old_elements.size++;

					/* terms in the set formulae of the ancestors of set `formula_set_id` could also be no longer provable */
					array<unsigned int> stack(8);
					hash_set<unsigned int> visited(16);
					stack[stack.length++] = i;
					visited.add(i);
					while (stack.length != 0) {
						unsigned int current = stack.pop();

						if (current != i) {
							Formula* set_formula = sets.sets[current].set_formula();
							if (set_formula->type == FormulaType::AND) {
								if (!conjuncts.append(set_formula->array.operands, set_formula->array.length)) {
									for (auto entry : old_elements) core::free(entry.value);
									return false;
								}
							} else if (!conjuncts.add(set_formula)) {
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
						}

						for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
							if (visited.contains(parent)) continue;
							if (!visited.add(parent) || !stack.add(parent)) {
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
						} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
							if (visited.contains(entry.key)) continue;
							if (!visited.add(entry.key) || !stack.add(entry.key)) {
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
						}
					}

					Term** src_variables = (Term**) alloca(sizeof(Term*) * sets.sets[i].arity * 2);
					for (unsigned int k = 0; k < sets.sets[i].arity; k++) {
						src_variables[k] = Term::new_variable(k + 1);
						if (src_variables[k] == nullptr) {
							for (unsigned int j = 0; j < k; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					}
					Term** dst_constants = src_variables + sets.sets[i].arity;
					for (unsigned int k = 0; k < sets.sets[i].arity; k++) {
						tuple_element& element_value = old_elements.values[old_elements.size - 1][k];
						if (element_value.type == tuple_element_type::CONSTANT)
							dst_constants[k] = Term::new_constant(element_value.constant);
						else if (element_value.type == tuple_element_type::NUMBER)
							dst_constants[k] = Term::new_number(element_value.number);
						else if (element_value.type == tuple_element_type::STRING)
							dst_constants[k] = Term::new_string(element_value.str);
						if (dst_constants[k] == nullptr) {
							for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
							for (unsigned int j = 0; j < k; j++) { core::free(*dst_constants[j]); core::free(dst_constants[j]); }
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					}
					for (Formula* conjunct : conjuncts) {
						Formula* substituted_conjunct = substitute_all(conjunct, src_variables, dst_constants, sets.sets[i].arity);
						if (substituted_conjunct == nullptr
						 || !old_atoms.ensure_capacity(old_atoms.length + 1))
						{
							if (substituted_conjunct != nullptr) { core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct); }
							for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
							for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						/* make sure `substituted_conjunct` doesn't already exist in `old_atoms` to avoid redundant computation */
						bool in_old_atoms = false;
						for (unsigned int j = old_atom_index + 1; !in_old_atoms && j < old_atoms.length; j++)
							if (*old_atoms[j] == *substituted_conjunct) in_old_atoms = true;
						if (in_old_atoms) {
							core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct);
							continue;
						}
						old_atoms[old_atoms.length++] = substituted_conjunct;
					}
					for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
				}
				for (auto& element : unifications) core::free(element);
			}
		}

		bool sets_changed = true;
		while (sets_changed) {
			sets_changed = false;
			for (unsigned int i = 0; i < old_elements.size; i++) {
				/* check if this removed element is still proveable */
				Formula* set_formula = sets.sets[old_elements.keys[i]].size_axioms[0]->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[old_elements.keys[i]].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[old_elements.keys[i]].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}

				array<instantiation_tuple> temp_possible_values(1);
				instantiation_tuple& values = temp_possible_values[0];
				if (!::init(values, sets.sets[old_elements.keys[i]].arity)) {
					for (uint_fast8_t k = 0; k < sets.sets[old_elements.keys[i]].arity; k++)
						move(old_elements.values[i][k], sets.sets[old_elements.keys[i]].elements[sets.sets[old_elements.keys[i]].elements.length++]);
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}
				for (unsigned int k = 0; k < values.length; k++) {
					core::free(values.values[k]);
					if (old_elements.values[i][k].type == tuple_element_type::CONSTANT) {
						values.values[k].type = instantiation_type::CONSTANT;
						values.values[k].constant = old_elements.values[i][k].constant;
					} else if (old_elements.values[i][k].type == tuple_element_type::NUMBER) {
						values.values[k].type = instantiation_type::NUMBER;
						values.values[k].number = old_elements.values[i][k].number;
					} else if (old_elements.values[i][k].type == tuple_element_type::STRING) {
						values.values[k].type = instantiation_type::STRING;
						if (!core::init(values.values[k].str, old_elements.values[i][k].str)) {
							values.values[k].type = instantiation_type::CONSTANT;
							for (uint_fast8_t k = 0; k < sets.sets[old_elements.keys[i]].arity; k++)
								move(old_elements.values[i][k], sets.sets[old_elements.keys[i]].elements[sets.sets[old_elements.keys[i]].elements.length++]);
							for (auto entry : old_elements) core::free(entry.value);
							core::free(values); return false;
						}
					}
				}
				temp_possible_values.length++;

				default_prover prover(sets);
				if (is_provable_without_abduction<false>(set_formula, quantifiers, temp_possible_values, prover)) {
					/* this element is still provably a member of this set */
					for (auto& element : temp_possible_values) core::free(element);
					if (!sets.sets[old_elements.keys[i]].add_element(old_elements.values[i])) {
						for (auto entry : old_elements) core::free(entry.value);
						return false;
					}
					core::free(old_elements.values[i]);
					old_elements.remove_at(i--);
					sets_changed = true;
					continue;
				}

				/* this removed element could belong to an ancestor */
				/* first get the sets that currently contain `element` */
				array<unsigned int> stack(8);
				array<unsigned int> old_sets(8);
				hash_set<unsigned int> visited(16);
				stack[stack.length++] = 1;
				visited.add(1);
				while (stack.length != 0) {
					unsigned int current = stack.pop();
					if (sets.sets[current].index_of_element(old_elements.values[i]) < sets.sets[current].element_count()) {
						if (!old_sets.add(current)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						continue;
					}
					bool has_containing_descendant = false;
					for (unsigned int old_set : old_sets) {
						if (sets.sets[current].descendants.contains(old_set)) {
							has_containing_descendant = true;
							break;
						}
					}
					if (has_containing_descendant) continue;

					for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
						if (visited.contains(parent)) continue;
						if (!visited.add(parent) || !stack.add(parent)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
						if (visited.contains(entry.key)) continue;
						if (!visited.add(entry.key) || !stack.add(entry.key)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					}
				}
				visited.clear();

				array<unsigned int> new_sets(8);
				stack[stack.length++] = old_elements.keys[i];
				visited.add(old_elements.keys[i]);
				while (stack.length != 0) {
					unsigned int current = stack.pop();
					bool is_ancestor_of_set_containing_element = false;
					for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < new_sets.length; j++)
						if (sets.sets[current].descendants.contains(new_sets[j])) is_ancestor_of_set_containing_element = true;
					for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < old_sets.length; j++)
						if (sets.sets[current].descendants.contains(old_sets[j])) is_ancestor_of_set_containing_element = true;
					if (is_ancestor_of_set_containing_element) continue;

					Formula* current_set_formula = sets.sets[current].set_formula();
					if (current_set_formula != set_formula && sets.sets[current].arity == sets.sets[old_elements.keys[i]].arity) {
						array<instantiation_tuple> temp_possible_values(1);
						instantiation_tuple& values = temp_possible_values[0];
						if (!::init(values, sets.sets[old_elements.keys[i]].arity)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						for (unsigned int k = 0; k < values.length; k++) {
							core::free(values.values[k]);
							if (old_elements.values[i][k].type == tuple_element_type::CONSTANT) {
								values.values[k].type = instantiation_type::CONSTANT;
								values.values[k].constant = old_elements.values[i][k].constant;
							} else if (old_elements.values[i][k].type == tuple_element_type::NUMBER) {
								values.values[k].type = instantiation_type::NUMBER;
								values.values[k].number = old_elements.values[i][k].number;
							} else if (old_elements.values[i][k].type == tuple_element_type::STRING) {
								values.values[k].type = instantiation_type::STRING;
								if (!core::init(values.values[k].str, old_elements.values[i][k].str)) {
									values.values[k].type = instantiation_type::CONSTANT;
									for (auto entry : old_elements) core::free(entry.value);
									core::free(values); return false;
								}
							}
						}
						temp_possible_values.length++;

						default_prover prover(sets);
						if (is_provable_without_abduction<false>(current_set_formula, quantifiers, temp_possible_values, prover)) {
							for (auto& element : temp_possible_values) core::free(element);
							for (unsigned int j = 0; j < new_sets.length; j++) {
								if (sets.sets[new_sets[j]].descendants.contains(current))
									new_sets.remove(j--);
							}
							if (!new_sets.add(current)) {
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
							continue;
						}
					}

					for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
						if (visited.contains(parent)) continue;
						if (!visited.add(parent) || !stack.add(parent)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
						if (visited.contains(entry.key)) continue;
						if (!visited.add(entry.key) || !stack.add(entry.key)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					}
				}

				for (unsigned int new_set : new_sets) {
					if (!sets.sets[new_set].add_element(old_elements.values[i])) {
						for (auto entry : old_elements) core::free(entry.value);
						return false;
					}
					/* this element could also be in `old_elements`, so remove it
					from `old_elements` to avoid the same element being added to
					the set multiple times */
					for (unsigned int j = i + 1; j < old_elements.size; j++) {
						if (old_elements.keys[j] == new_set && old_elements.values[j] == old_elements.values[i]) {
							core::free(old_elements.values[j]);
							old_elements.remove_at(j--);
						}
					}
				}
			}
		}
		for (auto entry : old_elements) core::free(entry.value);
		return true;
	}

	template<typename... Args>
	inline bool check_set_membership_after_subtraction(Formula* old_atom, Args&&... visitor)
	{
		array<Formula*> old_atoms(8);
		old_atoms[old_atoms.length++] = old_atom;
		old_atom->reference_count++;
		bool result = check_set_membership_after_subtraction(old_atoms, std::forward<Args>(visitor)...);
		free_all(old_atoms);
		return result;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_new_set_membership(unsigned int set_id, Args&&... visitor)
	{
		Formula* set_formula = sets.sets[set_id].size_axioms[0]->formula->binary.left->binary.right;

		array<Formula*> quantifiers(1 << (core::log2(sets.sets[set_id].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[set_id].arity; j++) {
			quantifiers[quantifiers.length++] = set_formula;
			set_formula = set_formula->quantifier.operand;
		}

		array<instantiation_tuple> possible_values(4);
		instantiation_tuple& values = possible_values[0];
		if (!::init(values, sets.sets[set_id].arity))
			return false;
		possible_values.length++;

		default_prover prover(sets);
		if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values, prover)) {
			for (auto& element : possible_values) core::free(element);
			return true;
		}

		array_map<tuple, array<unsigned int>> new_elements(1);
		array<Formula*> new_atoms(1);
		if (!check_set_membership<ResolveInconsistencies, true>(set_id, possible_values, new_elements, quantifiers, new_atoms, 0)) {
			for (auto& element : possible_values) core::free(element);
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			return false;
		}
		for (auto& element : possible_values) core::free(element);
		free_all(new_atoms);

		/* make sure we don't add elements that belong to descendants of this set */
		for (unsigned int descendant : sets.sets[set_id].descendants) {
			for (unsigned int i = 0; i < new_elements.size; i++) {
				if (sets.sets[descendant].index_of_element(new_elements.keys[i]) < sets.sets[descendant].element_count()) {
					core::free(new_elements.keys[i]);
					core::free(new_elements.values[i]);
					new_elements.remove_at(i--);
				}
			}
		}

		/* make sure the number of elements in `set_id` does not exceed its set_size */
		unsigned int new_element_count = 0;
		for (const auto& entry : new_elements) {
			if (!entry.value.contains(set_id)) continue;
			new_element_count++;
		}
		hash_set<tuple> provable_elements(16);
		if (!sets.get_provable_elements(set_id, provable_elements)
		 || provable_elements.size + new_element_count > sets.sets[set_id].set_size)
		{
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			for (tuple& tup : provable_elements) core::free(tup);
			return false;
		}
		for (tuple& tup : provable_elements) core::free(tup);

		/* add the elements to the sets */
		for (const auto& entry : new_elements) {
			if (!entry.value.contains(set_id)) continue;
			if (!sets.sets[set_id].add_element(entry.key)) {
				for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
				return false;
			}
		}
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
		return true;
	}

	template<typename... Args>
	bool check_old_set_membership(unsigned int set_id, Args&&... visitor)
	{
		Formula* set_formula = sets.sets[set_id].size_axioms[0]->formula->binary.left->binary.right;
		array<Formula*> quantifiers(1 << (core::log2(sets.sets[set_id].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[set_id].arity; j++) {
			quantifiers[quantifiers.length++] = set_formula;
			set_formula = set_formula->quantifier.operand;
		}

		for (unsigned int i = 0; i < sets.sets[set_id].element_count(); i++) {
			tuple_element* element = sets.sets[set_id].elements.data + (sets.sets[set_id].arity * i);

			/* first get the sets that currently contain `element` */
			array<unsigned int> stack(8);
			array<unsigned int> old_sets(8);
			hash_set<unsigned int> visited(16);
			stack[stack.length++] = 1;
			visited.add(1);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				if (current != set_id && sets.sets[current].index_of_element({element, sets.sets[set_id].arity}) < sets.sets[current].element_count()) {
					if (!old_sets.add(current))
						return false;
					continue;
				}
				bool has_containing_descendant = false;
				for (unsigned int old_set : old_sets) {
					if (sets.sets[current].descendants.contains(old_set)) {
						has_containing_descendant = true;
						break;
					}
				}
				if (has_containing_descendant) continue;

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent))
						return false;
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key))
						return false;
				}
			}
			visited.clear();

			array<unsigned int> new_sets(8);
			stack[stack.length++] = set_id;
			visited.add(set_id);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				bool is_ancestor_of_set_containing_element = false;
				for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < new_sets.length; j++)
					if (sets.sets[current].descendants.contains(new_sets[j])) is_ancestor_of_set_containing_element = true;
				for (unsigned int j = 0; !is_ancestor_of_set_containing_element && j < old_sets.length; j++)
					if (sets.sets[current].descendants.contains(old_sets[j])) is_ancestor_of_set_containing_element = true;
				if (is_ancestor_of_set_containing_element) continue;

				if (current != set_id && sets.sets[current].arity == sets.sets[set_id].arity) {
					array<instantiation_tuple> temp_possible_values(1);
					instantiation_tuple& values = temp_possible_values[0];
					if (!::init(values, sets.sets[set_id].arity))
						return false;
					for (unsigned int k = 0; k < values.length; k++) {
						core::free(values.values[k]);
						if (element[k].type == tuple_element_type::CONSTANT) {
							values.values[k].type = instantiation_type::CONSTANT;
							values.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							values.values[k].type = instantiation_type::NUMBER;
							values.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							values.values[k].type = instantiation_type::STRING;
							if (!core::init(values.values[k].str, element[k].str)) {
								values.values[k].type = instantiation_type::CONSTANT;
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;

					default_prover prover(sets);
					if (is_provable_without_abduction<false>(sets.sets[current].set_formula(), quantifiers, temp_possible_values, prover)) {
						for (auto& element : temp_possible_values) core::free(element);
						for (unsigned int j = 0; j < new_sets.length; j++) {
							if (sets.sets[new_sets[j]].descendants.contains(current))
								new_sets.remove(j--);
						}
						if (!new_sets.add(current))
							return false;
						continue;
					}
				}

				for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
					if (visited.contains(parent)) continue;
					if (!visited.add(parent) || !stack.add(parent))
						return false;
				} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
					if (visited.contains(entry.key)) continue;
					if (!visited.add(entry.key) || !stack.add(entry.key))
						return false;
				}
			}

			for (unsigned int new_set : new_sets) {
				if (!sets.sets[new_set].add_element({element, sets.sets[set_id].arity}))
					return false;
			}
			sets.sets[set_id].remove_element_at(i--);
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_new_subset_membership(unsigned int antecedent_set, unsigned int consequent_set, Args&&... visitor)
	{
		/* this new extensional edge may cause the elements in the strongly
		   connected component containing `antecedent_set` may be moveable into
		   a descendant of `consequent_set` */
		hash_set<unsigned int> consequent_component(8);
		if (!sets.get_strongly_connected_component(consequent_set, consequent_component))
			return false;
		if (consequent_component.contains(antecedent_set))
			return true;
		/* check if any of the elements in `consequent_component` is provably an element of `antecedent_set` */
		array<instantiation_tuple> possible_values(8);
		for (unsigned int set_id : consequent_component) {
			if (!possible_values.ensure_capacity(possible_values.length + sets.sets[set_id].element_count())) {
				for (auto& element : possible_values) core::free(element);
				return false;
			}
			for (unsigned int i = 0; i < sets.sets[set_id].element_count(); i++) {
				tuple_element* element = sets.sets[set_id].elements.data + (sets.sets[set_id].arity * i);
				tuple& current_element = *((tuple*) alloca(sizeof(tuple)));
				current_element.elements = element;
				current_element.length = sets.sets[set_id].arity;

				instantiation_tuple& values = possible_values[possible_values.length];
				if (!::init(values, sets.sets[set_id].arity)) {
					for (auto& element : possible_values) core::free(element);
					return false;
				}
				for (unsigned int k = 0; k < values.length; k++) {
					core::free(values.values[k]);
					if (current_element[k].type == tuple_element_type::CONSTANT) {
						values.values[k].type = instantiation_type::CONSTANT;
						values.values[k].constant = current_element[k].constant;
					} else if (current_element[k].type == tuple_element_type::NUMBER) {
						values.values[k].type = instantiation_type::NUMBER;
						values.values[k].number = current_element[k].number;
					} else if (current_element[k].type == tuple_element_type::STRING) {
						values.values[k].type = instantiation_type::STRING;
						if (!core::init(values.values[k].str, current_element[k].str)) {
							values.values[k].type = instantiation_type::CONSTANT;
							for (auto& element : possible_values) core::free(element);
							core::free(values); return false;
						}
					}
				}
				possible_values.length++;
			}
		}
		if (possible_values.length == 0) return true;

		Formula* antecedent_formula = sets.sets[antecedent_set].size_axioms[0]->formula->binary.left->binary.right;

		array<Formula*> quantifiers(1 << (core::log2(sets.sets[antecedent_set].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[antecedent_set].arity; j++) {
			quantifiers[quantifiers.length++] = antecedent_formula;
			antecedent_formula = antecedent_formula->quantifier.operand;
		}
		default_prover prover(sets);
		if (!is_provable_without_abduction<false>(antecedent_formula, quantifiers, possible_values, prover))
			return true;

		array<Formula*> new_atoms(4);
		array_map<tuple, array<unsigned int>> new_elements(4);
		if (!check_set_membership<ResolveInconsistencies, false>(consequent_set, possible_values, new_elements, quantifiers, new_atoms, 0)) {
			for (auto& element : possible_values) core::free(element);
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			free_all(new_atoms); return false;
		}
		for (auto& element : possible_values) core::free(element);
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
		bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, std::forward<Args>(visitor)...);
		free_all(new_atoms);
		return result;
	}

	template<typename... Args>
	bool check_old_subset_membership(unsigned int antecedent_set, unsigned int consequent_set, Args&&... visitor)
	{
		Formula* consequent_formula = sets.sets[consequent_set].size_axioms[0]->formula->binary.left->binary.right;
		array<Formula*> consequent_quantifiers(1 << (core::log2(sets.sets[consequent_set].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[consequent_set].arity; j++) {
			consequent_quantifiers[consequent_quantifiers.length++] = consequent_formula;
			consequent_formula = consequent_formula->quantifier.operand;
		}

		array<Formula*> old_atoms(8);
		for (unsigned int descendant : sets.sets[antecedent_set].descendants) {
			for (unsigned int i = sets.sets[descendant].element_count(); i > 0; i--) {
				const tuple_element* element = sets.sets[descendant].elements.data + (sets.sets[descendant].arity * (i - 1));
				Term** src_variables = (Term**) alloca(sizeof(Term*) * sets.sets[descendant].arity * 2);
				for (unsigned int i = 0; i < sets.sets[descendant].arity; i++) {
					src_variables[i] = Term::new_variable(i + 1);
					if (src_variables[i] == nullptr) {
						for (unsigned int j = 0; j < i; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
						free_all(old_atoms); return false;
					}
				}
				Term** dst_constants = src_variables + sets.sets[descendant].arity;
				for (unsigned int i = 0; i < sets.sets[descendant].arity; i++) {
					if (element[i].type == tuple_element_type::CONSTANT)
						dst_constants[i] = Term::new_constant(element[i].constant);
					else if (element[i].type == tuple_element_type::NUMBER)
						dst_constants[i] = Term::new_number(element[i].number);
					else if (element[i].type == tuple_element_type::STRING)
						dst_constants[i] = Term::new_string(element[i].str);
					if (dst_constants[i] == nullptr) {
						for (unsigned int j = 0; j < sets.sets[descendant].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
						for (unsigned int j = 0; j < i; j++) { core::free(*dst_constants[j]); core::free(dst_constants[j]); }
						free_all(old_atoms); return false;
					}
				}

				Formula* substituted_formula = substitute_all(consequent_formula, src_variables, dst_constants, sets.sets[descendant].arity);
				for (unsigned int j = 0; j < sets.sets[descendant].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
				for (unsigned int j = 0; j < sets.sets[descendant].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
				if (substituted_formula == nullptr) {
					free_all(old_atoms);
					return false;
				} else if (substituted_formula->type == FormulaType::AND) {
					if (!old_atoms.ensure_capacity(old_atoms.length + substituted_formula->array.length)) {
						core::free(*substituted_formula); if (substituted_formula->reference_count == 0) core::free(substituted_formula);
						free_all(old_atoms); return false;
					}
					for (unsigned int k = 0; k < substituted_formula->array.length; k++) {
						/* make sure `substituted_conjunct` doesn't already exist in `old_atoms` to avoid redundant computation */
						bool in_old_atoms = false;
						for (unsigned int j = 0; !in_old_atoms && j < old_atoms.length; j++)
							if (*old_atoms[j] == *substituted_formula->array.operands[k]) in_old_atoms = true;
						if (in_old_atoms) continue;
						old_atoms[old_atoms.length++] = substituted_formula->array.operands[k];
						substituted_formula->array.operands[k]->reference_count++;
					}
					core::free(*substituted_formula); if (substituted_formula->reference_count == 0) core::free(substituted_formula);
				} else {
					if (!old_atoms.ensure_capacity(old_atoms.length + 1)) {
						core::free(*substituted_formula); if (substituted_formula->reference_count == 0) core::free(substituted_formula);
						free_all(old_atoms); return false;
					}
					/* make sure `substituted_conjunct` doesn't already exist in `old_atoms` to avoid redundant computation */
					bool in_old_atoms = false;
					for (unsigned int j = 0; !in_old_atoms && j < old_atoms.length; j++)
						if (*old_atoms[j] == *substituted_formula) in_old_atoms = true;
					if (in_old_atoms) continue;
					old_atoms[old_atoms.length++] = substituted_formula;
				}
			}
		}

		bool result = check_set_membership_after_subtraction(old_atoms, std::forward<Args>(visitor)...);
		free_all(old_atoms);
		return result;
	}

	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_proof(Formula* canonicalized, set_changes<Formula>& set_diff, unsigned int& new_constant, Args&&... args)
	{
		Term* predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			return make_atom_proof<DefinitionsAllowed, Contradiction, ResolveInconsistencies>(canonicalized, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::NOT) {
			return make_proof<!Contradiction, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->unary.operand, set_diff, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
					return NULL;

				Term* constant;
				Proof* exists_not_proof = make_exists_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, set_diff, constant, new_constant, std::forward<Args>(args)...);
				if (exists_not_proof == NULL) return NULL;
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				 && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				{
					free_proof(exists_not_proof, set_diff, std::forward<Args>(args)...);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					return NULL;
				}
				Proof* assumption = ProofCalculus::new_axiom(canonicalized);
				if (assumption == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, set_diff, std::forward<Args>(args)...);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					return NULL;
				}
				assumption->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_universal_elim(assumption, constant), exists_not_proof), assumption);
				core::free(*assumption); if (assumption->reference_count == 0) core::free(assumption);
				if (proof == NULL) {
					/* we need to undo the changes made by `make_exists_proof` */
					free_proof(exists_not_proof, set_diff, std::forward<Args>(args)...);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					return NULL;
				}
				core::free(*exists_not_proof); if (exists_not_proof->reference_count == 0) core::free(exists_not_proof);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				core::free(*constant); if (constant->reference_count == 0) core::free(constant);
				existential_intro_nodes[existential_intro_nodes.length++] = { canonicalized, proof };
				canonicalized->reference_count++;
				core::free(*exists_not_proof); if (exists_not_proof->reference_count == 0) core::free(exists_not_proof);
				return proof;
			}

			unsigned int variable = canonicalized->quantifier.variable;
			Formula* new_canonicalized = shift_bound_variables(canonicalized, -((int) (variable - 1)));

			unsigned int arity = 1;
			Formula* operand = new_canonicalized->quantifier.operand;
			while (operand->type == FormulaType::FOR_ALL) {
				operand = new_canonicalized->quantifier.operand;
				arity++;
			}

			if (operand->type == FormulaType::IF_THEN) {
				Proof* new_axiom;
				Formula* left = operand->binary.left;
				Formula* right = operand->binary.right;

				/* check if `left` or `right` are of the form `c(x)` and that `c` can be a set */
				if (left->type == FormulaType::UNARY_APPLICATION
				 && left->binary.left->type == FormulaType::CONSTANT
				 && left->binary.right->type == FormulaType::VARIABLE)
				{
					if (is_provably_not_a_set(left->binary.left->constant)) {
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return nullptr;
					}
					/* make sure `c` is not the arg1 of a name event */
					if (left->binary.left->constant >= new_constant_offset) {
						for (Proof* proof : ground_concepts[left->binary.left->constant - new_constant_offset].definitions) {
							if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
							 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
							 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
							 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
							 && is_name_event(proof->formula->binary.right->binary.right->constant))
							{
								core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
								return nullptr;
							}
						}
					}
				}
				if (right->type == FormulaType::UNARY_APPLICATION
				 && right->binary.left->type == FormulaType::CONSTANT
				 && right->binary.right->type == FormulaType::VARIABLE)
				{
					if (is_provably_not_a_set(right->binary.left->constant)) {
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return nullptr;
					}
					/* make sure `c` is not the arg1 of a name event */
					if (right->binary.left->constant >= new_constant_offset) {
						for (Proof* proof : ground_concepts[right->binary.left->constant - new_constant_offset].definitions) {
							if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
							 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
							 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
							 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
							 && is_name_event(proof->formula->binary.right->binary.right->constant))
							{
								core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
								return nullptr;
							}
						}
					}
				}

				Term const* predicate; Term const* arg1; Term const* arg2;
				bool atomic = is_atomic(*left, predicate, arg1, arg2);
				if (arity == 1 && atomic && arg1->type == TermType::VARIABLE
				 && arg1->variable == variable && arg2 == NULL
				 && predicate->type == TermType::CONSTANT
				 && predicate->constant == (unsigned int) built_in_predicates::UNKNOWN)
				{
					if (!DefinitionsAllowed) {
						fprintf(stderr, "theory.make_proof ERROR: Definitions are not allowed in this context.\n");
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return NULL;
					}

					/* this is a definition of a type */
					/* check the right-hand side is a valid definition */
					if (!valid_definition(right, variable)) {
						fprintf(stderr, "theory.make_proof ERROR: This is not a valid type definition.\n");
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return NULL;
					}

					new_constant = get_free_concept_id();

					/* TODO: definitions of new concepts should be biconditionals */
					if (!try_init_concept(new_constant, arity)) return NULL;

					new_axiom = get_subset_axiom<ResolveInconsistencies>(left, right, arity, set_diff, std::forward<Args>(args)...);
					core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				} else {
					/* make sure there are no unknown predicates */
					if (contains_constant(*left, (unsigned int) built_in_predicates::UNKNOWN)
					 || contains_constant(*right, (unsigned int) built_in_predicates::UNKNOWN))
					{
						fprintf(stderr, "theory.make_proof ERROR: Universally-quantified statement has unknown constants.\n");
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return NULL;
					}

					/* this is a formula of form `![x]:(t(x) => f(x))` */
					new_axiom = get_subset_axiom<ResolveInconsistencies>(left, right, arity, set_diff, std::forward<Args>(args)...);
					core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
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
					Proof* operand = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index_minus_one], set_diff, new_constant, std::forward<Args>(args)...);
					if (operand != NULL) {
						/* we found a disproof of a conjunct */
						Proof* axiom = ProofCalculus::new_axiom(canonicalized);
						if (axiom == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							core::free(*negated_conjunction); if (negated_conjunction->reference_count == 0) core::free(negated_conjunction);
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						}
						axiom->reference_count++;
						Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
								ProofCalculus::new_conjunction_elim(axiom, make_array_view(&index_minus_one, 1)), operand), axiom);
						core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
						if (proof == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							core::free(*negated_conjunction); if (negated_conjunction->reference_count == 0) core::free(negated_conjunction);
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						}
						core::free(*operand); if (operand->reference_count == 0) core::free(operand);
						proof->reference_count++;
						/* record the formula with this negated conjunction node */
						negated_conjunction_nodes[negated_conjunction_nodes.length++] = { negated_conjunction, proof };
						return proof;
					}

					if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
				}

				/* we couldn't find a disproof of any of the conjuncts */
				finished_constants(canonicalized, old_index_count, std::forward<Args>(args)...);
				core::free(*negated_conjunction); if (negated_conjunction->reference_count == 0) core::free(negated_conjunction);
				return NULL;
			}

			Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
			if (operands == NULL) {
				fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
				return NULL;
			}
			for (unsigned int i = 0; i < canonicalized->array.length; i++) {
				if (new_constant == 0)
					operands[i] = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], set_diff, new_constant, std::forward<Args>(args)...);
				else operands[i] = make_proof<false, false, ResolveInconsistencies>(canonicalized->array.operands[i], set_diff, new_constant, std::forward<Args>(args)...);

				if (operands[i] == NULL) {
					for (unsigned int j = i; j > 0; j--)
						/* undo the changes made by the recursive calls to `make_proof` */
						free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
					core::free(operands); return NULL;
				}
			}
			Proof* conjunction = ProofCalculus::new_conjunction_intro(make_array_view(operands, canonicalized->array.length));
			if (conjunction == NULL) {
				for (unsigned int j = canonicalized->array.length; j > 0; j--)
					/* undo the changes made by the recursive calls to `make_proof` */
					free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
				core::free(operands); return NULL;
			}
			for (unsigned int j = 0; j < canonicalized->array.length; j++) {
				core::free(*operands[j]); if (operands[j]->reference_count == 0) core::free(operands[j]);
			}
			core::free(operands);
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
						operands[i] = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], set_diff, new_constant, std::forward<Args>(args)...);
					else operands[i] = make_proof<true, false, ResolveInconsistencies>(canonicalized->array.operands[i], set_diff, new_constant, std::forward<Args>(args)...);

					if (operands[i] == NULL) {
						for (unsigned int j = i; j > 0; j--)
							/* undo the changes made by recursive calls to `make_proof` */
							free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
						core::free(operands); return NULL;
					} else if (operands[i]->type == ProofType::PROOF_BY_CONTRADICTION) {
						Proof* absurdity = operands[i]->operands[0];
						absurdity->reference_count++;
						core::free(*operands[i]); if (operands[i]->reference_count == 0) core::free(operands[i]);
						operands[i] = absurdity;
					} else {
						Proof* absurdity = ProofCalculus::new_negation_elim(
								ProofCalculus::new_axiom(canonicalized->array.operands[i]), operands[i]);
						if (absurdity == NULL) {
							for (unsigned int j = i; j > 0; j--)
								/* undo the changes made by recursive calls to `make_proof` */
								free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
							core::free(operands); return NULL;
						}
						absurdity->reference_count++;
						core::free(*operands[i]); if (operands[i]->reference_count == 0) core::free(operands[i]);
						operands[i] = absurdity;
					}
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
					core::free(operands); return NULL;
				}
				axiom->reference_count++;
				Proof* proof = ProofCalculus::new_proof_by_contradiction(
						ProofCalculus::new_disjunction_elim(axiom, make_array_view(operands, canonicalized->array.length)), axiom);
				core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
				if (proof == NULL) {
					for (unsigned int j = canonicalized->array.length; j > 0; j--)
						/* undo the changes made by recursive calls to `make_proof` */
						free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
					core::free(operands); return NULL;
				}
				proof->reference_count++;
				for (unsigned int j = 0; j < canonicalized->array.length; j++) {
					core::free(*operands[j]); if (operands[j]->reference_count == 0) core::free(operands[j]);
				}
				core::free(operands);
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
				Proof* operand = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index - 1], set_diff, new_constant, std::forward<Args>(args)...);
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
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						}
						for (Formula* other_disjunct : other_disjuncts)
							other_disjunct->reference_count++;
					}
					Proof* proof = ProofCalculus::new_disjunction_intro(operand, other_disjunction, index - 1);
					core::free(*other_disjunction); if (other_disjunction->reference_count == 0) core::free(other_disjunction);
					if (proof == NULL) {
						/* undo changes made by the recursive call to `make_proof` */
						free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
					}
					core::free(*operand); if (operand->reference_count == 0) core::free(operand);
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
				Proof* left = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, set_diff, new_constant, std::forward<Args>(args)...);
				if (left == NULL) return NULL;

				Proof* right;
				if (new_constant == 0 && DefinitionsAllowed)
					right = make_proof<true, true, ResolveInconsistencies>(canonicalized->binary.right, set_diff, new_constant, std::forward<Args>(args)...);
				else right = make_proof<true, false, ResolveInconsistencies>(canonicalized->binary.right, set_diff, new_constant, std::forward<Args>(args)...);

				if (right == NULL) {
					free_proof(left, set_diff, std::forward<Args>(args)...); return NULL;
				}

				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == NULL) {
					free_proof(right, set_diff, std::forward<Args>(args)...);
					free_proof(left, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				axiom->reference_count++;

				Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
						ProofCalculus::new_implication_elim(axiom, left), right), axiom);
				core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
				if (proof == NULL) {
					free_proof(right, set_diff, std::forward<Args>(args)...);
					free_proof(left, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				core::free(*left); if (left->reference_count == 0) core::free(left);
				core::free(*right); if (right->reference_count == 0) core::free(right);
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
					Proof* left = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, set_diff, new_constant, std::forward<Args>(args)...);
					if (left == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* axiom = ProofCalculus::new_axiom(canonicalized->binary.left);
					if (axiom == NULL) {
						free_proof(left, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					axiom->reference_count++;

					Proof* proof = ProofCalculus::new_implication_intro(
							ProofCalculus::new_falsity_elim(
								ProofCalculus::new_negation_elim(left, axiom),
								canonicalized->binary.right),
							axiom);
					core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
					if (proof == NULL) {
						free_proof(left, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					core::free(*left); if (left->reference_count == 0) core::free(left);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length++] = { canonicalized, proof };
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				} else {
					Proof* right = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.right, set_diff, new_constant, std::forward<Args>(args)...);
					if (right == NULL) {
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					Proof* proof = ProofCalculus::new_implication_intro(right, ProofCalculus::new_axiom(canonicalized->binary.left));
					if (proof == NULL) {
						free_proof(right, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					core::free(*right); if (right->reference_count == 0) core::free(right);
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
				unsigned int min_variable = UINT_MAX;
				min_bound_variable(*canonicalized, min_variable);
				Formula* new_canonicalized;
				if (min_variable != UINT_MAX) {
					new_canonicalized = shift_bound_variables(canonicalized, -((int) (min_variable - 1)));
				} else {
					new_canonicalized = canonicalized;
					canonicalized->reference_count++;
				}

				unsigned int arity = 1;
				Formula* set_formula = new_canonicalized->quantifier.operand;
				while (set_formula->type == FormulaType::EXISTS) {
					set_formula = set_formula->quantifier.operand;
					arity++;
				}
				Formula* lambda_formula = set_formula;
				lambda_formula->reference_count++;
				for (unsigned int i = arity; i > 0; i--) {
					Formula* temp = Formula::new_lambda(i, lambda_formula);
					if (lambda_formula == NULL) {
						core::free(*lambda_formula); if (lambda_formula->reference_count == 0) core::free(lambda_formula);
						core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
						return NULL;
					}
					lambda_formula = temp;
				}
				core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
				unsigned int set_id; bool is_set_new;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(set_formula, arity, 0, set_id, is_set_new, set_diff, std::forward<Args>(args)...);
				if (set_size_axiom == NULL) {
					core::free(*lambda_formula); if (lambda_formula->reference_count == 0) core::free(lambda_formula);
					return NULL;
				} else if (is_set_new) {
					if (!check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
						core::free(*lambda_formula); if (lambda_formula->reference_count == 0) core::free(lambda_formula);
						sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
					}
				}

				Formula* beta_left = lambda_formula;
				for (unsigned int i = 0; i < arity; i++) {
					Formula* temp = Formula::new_apply(beta_left, Formula::new_variable(i + 1));
					if (temp == nullptr) {
						core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
						if (is_set_new) {
							check_old_set_membership(set_id, std::forward<Args>(args)...);
						}
						sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
					}
					beta_left = temp;
				} for (unsigned int i = arity; i > 0; i--) {
					Formula* temp = Formula::new_exists(i, beta_left);
					if (temp == nullptr) {
						core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
						if (is_set_new) {
							check_old_set_membership(set_id, std::forward<Args>(args)...);
						}
						sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
					}
					beta_left = temp;
				}
				Formula* temp = Formula::new_not(beta_left);
				if (temp == nullptr) {
					core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
					if (is_set_new) {
						check_old_set_membership(set_id, std::forward<Args>(args)...);
					}
					sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
				}
				beta_left = temp;

				Formula* beta_right = set_formula;
				set_formula->reference_count++;
				for (unsigned int i = arity; i > 0; i--) {
					Formula* temp = Formula::new_exists(i, beta_right);
					if (temp == nullptr) {
						core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
						core::free(*beta_right); if (beta_right->reference_count == 0) core::free(beta_right);
						if (is_set_new) {
							check_old_set_membership(set_id, std::forward<Args>(args)...);
						}
						sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
					}
					beta_right = temp;
				}
				temp = Formula::new_not(beta_right);
				if (temp == nullptr) {
					core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
					core::free(*beta_right); if (beta_right->reference_count == 0) core::free(beta_right);
					if (is_set_new) {
						check_old_set_membership(set_id, std::forward<Args>(args)...);
					}
					sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
				}
				beta_right = temp;

				Proof* proof = ProofCalculus::new_equality_elim(
						ProofCalculus::new_beta(beta_left, beta_right),
						ProofCalculus::new_equality_elim(ProofCalculus::new_universal_elim(empty_set_axiom, lambda_formula), set_size_axiom, make_repeated_array_view(0u, 1)),
						make_repeated_array_view(0u, 1));
				core::free(*beta_left); if (beta_left->reference_count == 0) core::free(beta_left);
				core::free(*beta_right); if (beta_right->reference_count == 0) core::free(beta_right);
				core::free(*lambda_formula); if (lambda_formula->reference_count == 0) core::free(lambda_formula);
				if (proof == NULL) {
					if (is_set_new) {
						check_old_set_membership(set_id, std::forward<Args>(args)...);
					}
					sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
				}
				lambda_formula->reference_count++;
				proof->reference_count++;
				return proof;
			} else {
				Term* variable = Formula::new_variable(canonicalized->quantifier.variable);
				if (variable == NULL) return NULL;

				Term* constant;
				Proof* operand = make_exists_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, canonicalized->quantifier.variable, set_diff, constant, new_constant, std::forward<Args>(args)...);
				if (operand == NULL) {
					core::free(*variable); if (variable->reference_count == 0) core::free(variable);
					return NULL;
				}

				array<unsigned int> indices(8);
				if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1)
				 || (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				  && !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1))
				 || !compute_indices(*canonicalized->quantifier.operand, *variable, indices))
				{
					core::free(*variable); if (variable->reference_count == 0) core::free(variable);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					free_proof(operand, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				core::free(*variable); if (variable->reference_count == 0) core::free(variable);
				Proof* proof = ProofCalculus::new_existential_intro(operand, indices.data, indices.length, constant);
				if (proof == NULL) {
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					free_proof(operand, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				core::free(*operand); if (operand->reference_count == 0) core::free(operand);
				proof->reference_count++;
				/* record the formula with this existential */
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset)
					ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.add(proof);
				core::free(*constant); if (constant->reference_count == 0) core::free(constant);
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
			 && right->type == TermType::NUMBER)
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
				if (right->number.decimal != 0)
					return nullptr;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(arg1->quantifier.operand, arity, right->number.integer, set_id, is_set_new, set_diff, std::forward<Args>(args)...);
				if (set_size_axiom == nullptr) {
					return nullptr;
				} else if (is_set_new) {
					if (!check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
						sets.try_free_set(set_id, std::forward<Args>(args)...);
						return nullptr;
					}
				}
				set_size_axiom->reference_count++;
				return set_size_axiom;

			} else if (atomic && predicate->type == TermType::CONSTANT
					&& predicate->constant == (unsigned int) built_in_predicates::SIZE
			 		&& arg1->type == TermType::CONSTANT && arg2 == NULL
			 		&& right->type == TermType::NUMBER)
			{
				if (Contradiction) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}

				/* make sure that `arg1->constant` can be a set */
				if (is_provably_not_a_set(arg1->constant))
					return nullptr;
				/* make sure `arg1->constant` is not the arg1 of a name event */
				if (arg1->constant >= new_constant_offset) {
					for (Proof* proof : ground_concepts[arg1->constant - new_constant_offset].definitions) {
						if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
						 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
						 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
						 && is_name_event(proof->formula->binary.right->binary.right->constant))
						{
							return nullptr;
						}
					}
				}

				/* this is a statement on the size of a set */
				if (right->number.decimal != 0)
					return nullptr;
				Proof* definition = ground_concepts[arg1->constant - new_constant_offset].definitions[0];
#if !defined(NDEBUG)
				if (ground_concepts[arg1->constant - new_constant_offset].definitions.length == 0
				 || definition->formula->binary.right->type != FormulaType::LAMBDA)
					fprintf(stderr, "theory.make_proof ERROR: Expected a set definition of the form c=λx.c(x).\n");
#endif
				unsigned int arity = 1;
				Formula* set_formula = definition->formula->binary.right->quantifier.operand;
				while (set_formula->type == FormulaType::LAMBDA) {
					set_formula = set_formula->quantifier.operand;
					arity++;
				}
				unsigned int set_id; bool is_set_new;
				Proof* set_size_axiom = sets.template get_size_axiom<ResolveInconsistencies>(set_formula, arity, right->number.integer, set_id, is_set_new, set_diff, std::forward<Args>(args)...);
				Proof* proof = ProofCalculus::new_equality_elim(definition, set_size_axiom, make_repeated_array_view(3u, 1));
				if (proof == NULL)
					return NULL;
				proof->reference_count++;

				if (is_set_new && !check_new_set_membership<ResolveInconsistencies>(set_id, std::forward<Args>(args)...)) {
					free_proof(proof, set_diff);
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
					if (right->type == FormulaType::LAMBDA) {
						if (is_provably_not_a_set(left->constant))
							return nullptr;
						/* make sure `left->constant` is not the arg1 of a name event */
						if (left->constant >= new_constant_offset) {
							for (Proof* proof : ground_concepts[left->constant - new_constant_offset].definitions) {
								if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
								 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
								 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
								 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
								 && is_name_event(proof->formula->binary.right->binary.right->constant))
								{
									return nullptr;
								}
							}
						}
					}

					/* check to make sure the args are type-correct */
					if (right->type == TermType::UNARY_APPLICATION
					 && right->binary.left->type == TermType::CONSTANT
					 && right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					 && right->binary.right->type == TermType::CONSTANT)
					{
						/* make sure both `right->binary.right->constant` and `left->constant` are not name event */
						if (is_name_event(left->constant) && is_name_event(right->binary.right->constant))
							return nullptr;
						/* make sure the arg1 of a name event is not a set */
						if (is_name_event(right->binary.right->constant) && is_provably_a_set(left->constant))
							return nullptr;
					}
					if (right->type == TermType::UNARY_APPLICATION
					 && right->binary.left->type == TermType::CONSTANT
					 && right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					 && right->binary.right->type == TermType::CONSTANT)
					{
						/* make sure the event is not a `name` since the arg2 of `name` events must be a string */
						if (is_name_event(right->binary.right->constant))
							return nullptr;
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

					/* check if anything else has this definition */
					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						if (ground_concepts[i].types.keys == NULL || i == left->constant - new_constant_offset) continue;
						for (Proof* proof : ground_concepts[i].definitions) {
							if (proof->formula->binary.right == new_right || *proof->formula->binary.right == *new_right) {
								/* we found a different concept with this definition */
								core::free(*new_right); if (new_right->reference_count == 0) core::free(new_right);
								return NULL;
							}
						}
					}

					Formula* new_formula = Formula::new_equals(left, new_right);
					if (new_formula == NULL) {
						core::free(*new_right); if (new_right->reference_count == 0) core::free(new_right);
						return NULL;
					}
					left->reference_count++;

					Proof* new_proof = ProofCalculus::new_axiom(new_formula);
					core::free(*new_formula); if (new_formula->reference_count == 0) core::free(new_formula);
					if (new_proof == NULL)
						return NULL;
					new_proof->reference_count++;

					Proof* definition = add_definition<ResolveInconsistencies>(new_proof, set_diff, std::forward<Args>(args)...);
					if (definition != new_proof) {
						core::free(*new_proof);
						core::free(new_proof);
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

				/* check to make sure the args are type-correct */
				if (left->binary.right->constant >= new_constant_offset) {
					unsigned int event = left->binary.right->constant;
					if (left->binary.left->constant == (unsigned int) built_in_predicates::ARG1) {
						/* the event must not be `name` since `right` cannot be a constant here */
						if (is_name_event(event)) return nullptr;
					} else if (left->binary.left->constant == (unsigned int) built_in_predicates::ARG2) {
						if (right->type == TermType::STRING) {
							/* the event must be `name` */
							if (ground_concepts[event - new_constant_offset].types.size > 1)
								return nullptr;
							if (ground_concepts[event - new_constant_offset].types.size == 1) {
								Term& key = ground_concepts[event - new_constant_offset].types.keys[0];
								if (key.type != TermType::UNARY_APPLICATION
								 || key.binary.left->type != TermType::CONSTANT
								 || key.binary.left->constant != (unsigned int) built_in_predicates::NAME
								 || key.binary.right->type != TermType::VARIABLE
								 || key.binary.right->variable != 1)
								{
									return nullptr;
								}
							}
						} else {
							/* the event must not be `name` */
							if (is_name_event(event)) return nullptr;
						}
					}
				}

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
					core::free(*new_right); if (new_right->reference_count == 0) core::free(new_right);
					return NULL;
				}
				left->reference_count++;

				Proof* new_proof = ProofCalculus::new_axiom(new_formula);
				core::free(*new_formula); if (new_formula->reference_count == 0) core::free(new_formula);
				if (new_proof == NULL)
					return NULL;
				new_proof->reference_count++;

				Proof* function_value = add_function_value<ResolveInconsistencies>(new_proof, std::forward<Args>(args)...);
				if (function_value != new_proof) {
					core::free(*new_proof);
					core::free(new_proof);
				} if (function_value == nullptr) {
					return nullptr;
				}

				if (swap_order) {
					Proof* swapped_proof = ProofCalculus::new_equality_elim(
							function_value, ProofCalculus::new_beta(right, right), make_repeated_array_view(2u, 1));
					if (swapped_proof == NULL) {
						free_proof(function_value, set_diff);
						return NULL;
					}
					core::free(*function_value);
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
	inline Proof* get_subset_axiom_with_required_set_size(
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int& antecedent_set, unsigned int& consequent_set,
			bool& is_antecedent_new, bool& is_consequent_new,
			unsigned int set_size, Args&&... visitor)
	{
		required_set_size set_size_enforcer(set_size);
		if (!sets.get_set_id(antecedent, arity, antecedent_set, is_antecedent_new, set_size_enforcer, std::forward<Args>(visitor)...))
			return nullptr;
		if (is_antecedent_new && !check_new_set_membership<ResolveInconsistencies>(antecedent_set, std::forward<Args>(visitor)...)) {
			sets.try_free_set(antecedent_set);
			return nullptr;
		}

		if (!sets.get_set_id(consequent, arity, consequent_set, is_consequent_new, set_size_enforcer, std::forward<Args>(visitor)...)) {
			if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		} if (is_consequent_new && !check_new_set_membership<ResolveInconsistencies>(consequent_set, std::forward<Args>(visitor)...)) {
			sets.try_free_set(consequent_set);
			if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		bool is_subset_new = !sets.extensional_graph.vertices[consequent_set].children.contains(antecedent_set)
						  && !sets.intensional_graph.vertices[consequent_set].children.contains(antecedent_set);

		bool dummy;
		Proof* axiom = sets.template get_subset_axiom<ResolveInconsistencies, false>(
				antecedent, consequent, arity,
				antecedent_set, consequent_set,
				dummy, dummy, set_size_enforcer,
				std::forward<Args>(visitor)...);
		if (axiom == nullptr) {
			if (is_consequent_new) {
				check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(consequent_set);
			} if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		if (is_subset_new && !check_new_subset_membership<ResolveInconsistencies>(antecedent_set, consequent_set, std::forward<Args>(visitor)...)) {
			sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
			if (is_consequent_new) {
				check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(consequent_set);
			} if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		if (consequent->type == TermType::UNARY_APPLICATION && consequent->binary.right->type == TermType::VARIABLE && consequent->binary.right->variable == 1) {
			/* adding this edge could increase the number of provable elements
			   in `consequent_set` to be its set size, which can cause `a(b)`
			   to be newly provable or disprovable */
			hash_set<tuple> provable_elements(16);
			if (!sets.get_provable_elements(consequent_set, provable_elements)
			 || (provable_elements.size == sets.sets[consequent_set].set_size && !check_set_membership_after_addition<ResolveInconsistencies>(consequent, std::forward<Args>(visitor)...)))
			{
				for (tuple& tup : provable_elements) core::free(tup);
				if (is_subset_new) {
					sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
					check_old_subset_membership(antecedent_set, consequent_set, std::forward<Args>(visitor)...);
				} if (is_consequent_new) {
					check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
					sets.try_free_set(consequent_set);
				} if (is_antecedent_new) {
					check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
					sets.try_free_set(antecedent_set);
				}
				return nullptr;
			}
			for (tuple& tup : provable_elements) core::free(tup);
		}
		return axiom;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline Proof* get_subset_axiom(
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int& antecedent_set, unsigned int& consequent_set,
			bool& is_antecedent_new, bool& is_consequent_new, Args&&... visitor)
	{
		if (!sets.get_set_id(antecedent, arity, antecedent_set, is_antecedent_new, std::forward<Args>(visitor)...))
			return nullptr;
		if (is_antecedent_new && !check_new_set_membership<ResolveInconsistencies>(antecedent_set, std::forward<Args>(visitor)...)) {
			sets.try_free_set(antecedent_set);
			return nullptr;
		}

		if (!sets.get_set_id(consequent, arity, consequent_set, is_consequent_new, std::forward<Args>(visitor)...)) {
			if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		} if (is_consequent_new && !check_new_set_membership<ResolveInconsistencies>(consequent_set, std::forward<Args>(visitor)...)) {
			sets.try_free_set(consequent_set);
			if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		bool is_subset_new = !sets.extensional_graph.vertices[consequent_set].children.contains(antecedent_set)
						  && !sets.intensional_graph.vertices[consequent_set].children.contains(antecedent_set);

		bool dummy;
		Proof* axiom = sets.template get_subset_axiom<ResolveInconsistencies, false>(
				antecedent, consequent, arity,
				antecedent_set, consequent_set,
				dummy, dummy, std::forward<Args>(visitor)...);
		if (axiom == nullptr) {
			if (is_consequent_new) {
				check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(consequent_set);
			} if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		if (is_subset_new && !check_new_subset_membership<ResolveInconsistencies>(antecedent_set, consequent_set, std::forward<Args>(visitor)...)) {
			sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
			if (is_consequent_new) {
				check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(consequent_set);
			} if (is_antecedent_new) {
				check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
				sets.try_free_set(antecedent_set);
			}
			return nullptr;
		}

		if (consequent->type == TermType::UNARY_APPLICATION && consequent->binary.right->type == TermType::VARIABLE && consequent->binary.right->variable == 1) {
			/* adding this edge could increase the number of provable elements
			   in `consequent_set` to be its set size, which can cause `a(b)`
			   to be newly provable or disprovable */
			hash_set<tuple> provable_elements(16);
			if (!sets.get_provable_elements(consequent_set, provable_elements)
			 || (provable_elements.size == sets.sets[consequent_set].set_size && !check_set_membership_after_addition<ResolveInconsistencies>(consequent, std::forward<Args>(visitor)...)))
			{
				for (tuple& tup : provable_elements) core::free(tup);
				if (is_subset_new) {
					sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
					check_old_subset_membership(antecedent_set, consequent_set, std::forward<Args>(visitor)...);
				} if (is_consequent_new) {
					check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
					sets.try_free_set(consequent_set);
				} if (is_antecedent_new) {
					check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
					sets.try_free_set(antecedent_set);
				}
				return nullptr;
			}
			for (tuple& tup : provable_elements) core::free(tup);
		}
		return axiom;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline Proof* get_subset_axiom(Formula* antecedent, Formula* consequent, unsigned int arity, Args&&... visitor)
	{
		unsigned int antecedent_set, consequent_set;
		bool is_antecedent_new, is_consequent_new;
		return get_subset_axiom<ResolveInconsistencies>(antecedent, consequent, arity, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(visitor)...);
	}

	template<typename... Args>
	inline void free_subset_axiom(Formula* antecedent, Formula* consequent, unsigned int arity, unsigned int antecedent_set, unsigned int consequent_set, Args&&... visitor)
	{
		consequent->reference_count++;
		sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
		if (consequent->type == TermType::UNARY_APPLICATION && consequent->binary.right->type == TermType::VARIABLE && consequent->binary.right->variable == 1) {
			/* removing this edge could decrease the number of provable elements
			   in `consequent_set` to no longer be its set size, which can
			   cause `a(b)` to be newly provable or disprovable */
			hash_set<tuple> provable_elements(16);
			if (!sets.get_provable_elements(consequent_set, provable_elements)) {
				for (tuple& tup : provable_elements) core::free(tup);
				core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
				return;
			}
			unsigned int new_set_size = provable_elements.size;
			if (new_set_size < sets.sets[consequent_set].set_size) {
				if (!sets.get_provable_elements(antecedent_set, provable_elements)) {
					for (tuple& tup : provable_elements) core::free(tup);
					core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
					return;
				}
				if (provable_elements.size == sets.sets[consequent_set].set_size
				 && !check_set_membership_after_subtraction(consequent, std::forward<Args>(visitor)...))
				{
					for (tuple& tup : provable_elements) core::free(tup);
					core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
					return;
				}
			}
			for (tuple& tup : provable_elements) core::free(tup);
		}
		core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
		check_old_subset_membership(antecedent_set, consequent_set, std::forward<Args>(visitor)...);
		if (sets.is_freeable(consequent_set, std::forward<Args>(visitor)...)) {
			on_free_set(consequent_set, sets, std::forward<Args>(visitor)...);
			sets.free_set(consequent_set);
		} if (sets.is_freeable(antecedent_set, std::forward<Args>(visitor)...)) {
			on_free_set(antecedent_set, sets, std::forward<Args>(visitor)...);
			sets.free_set(antecedent_set);
		}
	}

private:
	template<bool ResolveInconsistencies, typename... Args>
	Proof* add_definition(Proof* definition, Args&&... args)
	{
		Formula* constant = definition->formula->binary.left;
		Formula* new_definition = definition->formula->binary.right;

		unsigned int arity;
		Formula* new_set_formula;
		if (new_definition->type == FormulaType::LAMBDA) {
			arity = 1;
			new_set_formula = new_definition->quantifier.operand;
			while (new_set_formula->type == FormulaType::LAMBDA) {
				new_set_formula = new_set_formula->quantifier.operand;
				arity++;
			}
			if (!try_init_concept(constant->constant, arity))
				return NULL;
		} else {
			if (!try_init_concept(constant->constant))
				return NULL;
		}

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
			unsigned int set_id = sets.set_ids.get(*new_set_formula, contains);
			if (contains) {
				set_size = sets.sets[set_id].set_size;
			} else if (ground_concepts[constant->constant - new_constant_offset].definitions.length != 0) {
				for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
					Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
					if (definition->formula->binary.right->type != FormulaType::LAMBDA)
						continue;
					Formula* other_set_formula = definition->formula->binary.right->quantifier.operand;
					while (other_set_formula->type == FormulaType::LAMBDA)
						other_set_formula = other_set_formula->quantifier.operand;
					set_id = sets.set_ids.get(*other_set_formula, contains);
					if (contains) {
						set_size = sets.sets[set_id].set_size;
						break;
					}
				}
			}

			for (unsigned int i = 0; i < ground_concepts[constant->constant - new_constant_offset].definitions.length; i++) {
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[i];
				if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
					Formula* other_set_formula = definition->formula->binary.right->quantifier.operand;
					while (other_set_formula->type == FormulaType::LAMBDA)
						other_set_formula = other_set_formula->quantifier.operand;

					unsigned int antecedent_set, consequent_set;
					bool is_antecedent_new, is_consequent_new;
					Proof* first_subset_axiom = get_subset_axiom_with_required_set_size<ResolveInconsistencies>(
							new_set_formula, other_set_formula, arity, antecedent_set, consequent_set,
							is_antecedent_new, is_consequent_new, set_size, std::forward<Args>(args)...);
					if (first_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							if (definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
							Proof* axiom = get_subset_axiom<false>(new_set_formula, other_set_formula, arity,
									antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
							core::free(*axiom);
							if (axiom->reference_count == 1)
								free_subset_axiom(new_set_formula, other_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);

							axiom = get_subset_axiom<false>(other_set_formula, new_set_formula, arity,
									antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
							core::free(*axiom);
							if (axiom->reference_count == 1)
								free_subset_axiom(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);
						}
						return NULL;
					}
					first_subset_axiom->reference_count++;

					Proof* second_subset_axiom = get_subset_axiom<ResolveInconsistencies>(other_set_formula, new_set_formula, arity,
							antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
					if (second_subset_axiom == NULL) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						core::free(*first_subset_axiom);
						if (first_subset_axiom->reference_count == 1)
							free_subset_axiom(new_set_formula, other_set_formula, arity, consequent_set, antecedent_set, std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							if (definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
							Proof* axiom = get_subset_axiom<false>(new_set_formula, other_set_formula, arity,
									antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
							core::free(*axiom);
							if (axiom->reference_count == 1)
								free_subset_axiom(new_set_formula, other_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);

							axiom = get_subset_axiom<false>(other_set_formula, new_set_formula, arity,
									antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
							core::free(*axiom);
							if (axiom->reference_count == 1)
								free_subset_axiom(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);
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
				unsigned int antecedent_set, consequent_set;
				bool is_antecedent_new, is_consequent_new;
				Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
				if (definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
				Formula* other_set_formula = definition->formula->binary.right->quantifier.operand;
				while (other_set_formula->type == FormulaType::LAMBDA)
					other_set_formula = other_set_formula->quantifier.operand;
				Proof* axiom = get_subset_axiom<false>(new_set_formula, other_set_formula, arity,
						antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
				core::free(*axiom);
				if (axiom->reference_count == 1)
					free_subset_axiom(new_set_formula, other_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);

				axiom = get_subset_axiom<false>(other_set_formula, new_set_formula, arity,
						antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
				core::free(*axiom);
				if (axiom->reference_count == 1)
					free_subset_axiom(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);
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

	inline unsigned int has_unremoved_extensional_edges(unsigned int set_id, const Formula* antecedent_to_remove, const Formula* consequent_to_remove)
	{
		for (const auto& entry : sets.extensional_graph.vertices[set_id].parents) {
			for (Proof* existing_axiom : entry.value) {
				Formula* formula = existing_axiom->formula->quantifier.operand;
				while (formula->type == FormulaType::FOR_ALL)
					formula = formula->quantifier.operand;
				if ((antecedent_to_remove != formula->binary.left && *antecedent_to_remove != *formula->binary.left)
				 || (consequent_to_remove != formula->binary.right && *consequent_to_remove != *formula->binary.right))
					return true;
			}
		} for (const auto& entry : sets.extensional_graph.vertices[set_id].children) {
			for (Proof* existing_axiom : entry.value) {
				Formula* formula = existing_axiom->formula->quantifier.operand;
				while (formula->type == FormulaType::FOR_ALL)
					formula = formula->quantifier.operand;
				if ((antecedent_to_remove != formula->binary.left && *antecedent_to_remove != *formula->binary.left)
				 || (consequent_to_remove != formula->binary.right && *consequent_to_remove != *formula->binary.right))
					return true;
			}
		}
		return false;
	}

	inline bool has_unremoved_extensional_edges(unsigned int set_id)
	{
		return sets.extensional_graph.vertices[set_id].parents.size != 0
			|| sets.extensional_graph.vertices[set_id].children.size != 0;
	}

	inline bool has_unremoved_extensional_edges(unsigned int set_id, unsigned int set_to_remove)
	{
		for (const auto& entry : sets.extensional_graph.vertices[set_id].parents) {
			if (entry.key != set_to_remove) return true;
		} for (const auto& entry : sets.extensional_graph.vertices[set_id].children) {
			if (entry.key != set_to_remove) return true;
		}
		return false;
	}

	template<typename... Args>
	void remove_definition(Proof* definition, Args&&... args)
	{
		unsigned int concept_id = definition->formula->binary.left->constant;

		unsigned int index;
		for (index = 0; index < ground_concepts[concept_id - new_constant_offset].definitions.length; index++) {
			Proof* other_definition = ground_concepts[concept_id - new_constant_offset].definitions[index];
			if (definition->formula->binary.right == other_definition->formula->binary.right) break;
		}

#if !defined(NDEBUG)
		if (index == ground_concepts[concept_id - new_constant_offset].definitions.length)
			fprintf(stderr, "remove_definition WARNING: Unable to find definition to remove.\n");
#endif

		ground_concepts[concept_id - new_constant_offset].definitions.remove(index);
		check_set_membership_after_subtraction(definition->formula, std::forward<Args>(args)...);

		/* remove subset edges from other set definitions for `concept_id` */
		for (unsigned int i = ground_concepts[concept_id - new_constant_offset].definitions.length; i > 0; i--) {
			Proof* other_definition = ground_concepts[concept_id - new_constant_offset].definitions[i - 1];
			if (other_definition->formula->binary.right->type != FormulaType::LAMBDA)
				continue;
			unsigned int arity = 1;
			Formula* other_set_formula = other_definition->formula->binary.right->quantifier.operand;
			while (other_set_formula->type == FormulaType::LAMBDA) {
				other_set_formula = other_set_formula->quantifier.operand;
				arity++;
			}
			if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
				Formula* set_formula = definition->formula->binary.right->quantifier.operand;
				while (set_formula->type == FormulaType::LAMBDA)
					set_formula = set_formula->quantifier.operand;

				unsigned int antecedent_set, consequent_set;
				bool is_antecedent_new, is_consequent_new;
				Proof* axiom = get_subset_axiom<false>(
						other_set_formula, set_formula, arity, antecedent_set, consequent_set,
						is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
				core::free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					free_subset_axiom(other_set_formula, set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);

				axiom = get_subset_axiom<false>(
						set_formula, other_set_formula, arity, antecedent_set, consequent_set,
						is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
				core::free(*axiom);
				if (axiom->reference_count == axiom->children.length + 1)
					free_subset_axiom(set_formula, other_set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);
			}
		}

		try_free_concept_id(concept_id);
		core::free(*definition); if (definition->reference_count == 0) core::free(definition);
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
		core::free(*function_value_axiom); if (function_value_axiom->reference_count == 0) core::free(function_value_axiom);
	}

	/* NOTE: this function finds a constant that proves `quantified`, and not ?[x]:`quantified` */
	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_exists_proof(Formula* quantified, unsigned int variable, set_changes<Formula>& set_diff, Term*& constant, unsigned int& new_constant, Args&&... args)
	{
		Formula* var = Formula::new_variable(variable);
		if (var == nullptr) return nullptr;

		array<instance> constants(ground_concept_capacity + 1);
		array<hol_number> numbers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != nullptr) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length++].constant = new_constant_offset + i;

				for (const auto& entry : ground_concepts[i].function_values) {
					Term* constant = entry.value->formula->binary.right;
					if (constant->type == TermType::NUMBER) {
						if (!numbers.contains(constant->number) && !numbers.add(constant->number)) return nullptr;
					} else if (constant->type == TermType::STRING) {
						bool contains = false;
						for (const string* str : strings)
							if (str == &constant->str || *str == constant->str) { contains = true; break; }
						if (!contains && !strings.add(&constant->str)) return nullptr;
					}
				}
			}
		} for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr) continue;
			hol_number number;
			number.integer = sets.sets[i].set_size;
			number.decimal = 0;
			if (!numbers.contains(number) && !numbers.add(number))
				return nullptr;
		}
		constants[constants.length++].type = instance_type::ANY;
		if (!constants.ensure_capacity(constants.length + numbers.length + strings.length))
			return nullptr;
		for (hol_number number : numbers) {
			constants[constants.length].type = instance_type::NUMBER;
			constants[constants.length++].number = number;
		} for (string* str : strings) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length++].str = str;
		}
		shuffle(constants);

		unsigned int original_constant_count = constants.length;
		if (!filter_constants(*this, quantified, variable, constants, std::forward<Args>(args)...)) {
			core::free(*var); if (var->reference_count == 0) core::free(var);
			return nullptr;
		}

		for (const instance& id : constants) {
			if (id.type == instance_type::ANY || (id.type == instance_type::CONSTANT && ground_concepts[id.constant - new_constant_offset].types.keys == nullptr)) {
				unsigned int constant_id;
				if (id.type == instance_type::ANY)
					constant_id = get_free_concept_id();
				else constant_id = id.constant;

				if (!try_init_concept(constant_id)) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return nullptr;
				}

				constant = Formula::new_constant(constant_id);
				if (constant == nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					free_concept_id(constant_id); return nullptr;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					free_concept_id(constant_id); return nullptr;
				}

				Proof* proof = make_proof<Contradiction, false, ResolveInconsistencies>(substituted, set_diff, new_constant, std::forward<Args>(args)...);
				core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
				if (proof != nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return proof;
				}
				core::free(*constant); if (constant->reference_count == 0) core::free(constant);
				if (ground_concepts[constant_id - new_constant_offset].types.keys != nullptr)
					/* `make_proof` could have already freed the new concept */
					free_concept_id(constant_id);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return nullptr;
				}
			} else {
				if (id.type == instance_type::CONSTANT) {
					constant = Formula::new_constant(id.constant);
				} else if (id.type == instance_type::NUMBER) {
					constant = Formula::new_number(id.number);
				} else if (id.type == instance_type::STRING) {
					constant = Formula::new_string(*id.str);
				}
				if (constant == nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return nullptr;
				}
				Formula* substituted = substitute(quantified, var, constant);
				if (substituted == nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					core::free(*constant); if (constant->reference_count == 0) core::free(constant);
					return nullptr;
				}

				Proof* proof = make_proof<Contradiction, DefinitionsAllowed, ResolveInconsistencies>(substituted, set_diff, new_constant, std::forward<Args>(args)...);
				core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
				if (proof != nullptr) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return proof;
				}
				core::free(*constant); if (constant->reference_count == 0) core::free(constant);

				if (!inconsistent_constant(quantified, id, std::forward<Args>(args)...)) {
					core::free(*var); if (var->reference_count == 0) core::free(var);
					return nullptr;
				}
			}
		}

		finished_constants(quantified, original_constant_count, std::forward<Args>(args)...);
		core::free(*var); if (var->reference_count == 0) core::free(var);
		return nullptr;
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
					if (!Negated && ((atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant)) || is_provably_a_set(arg1->constant)))
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
					core::free(*formula); if (formula->reference_count == 0) core::free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;

					if (!try_init_concept(new_constant, 1)) {
						core::free(*new_axiom); core::free(new_axiom);
						return NULL;
					} if (!add_unary_atom<Negated, ResolveInconsistencies>(Negated ? *formula->unary.operand : *formula, new_axiom, std::forward<Args>(args)...)) {
						core::free(*new_axiom); core::free(new_axiom);
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
						if (sets.sets[i].size_axioms.data == nullptr) continue;
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
							for (Proof* edge_axiom : entry.value) {
								for (Proof* child : edge_axiom->children) {
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
						core::free(*lifted_atom); core::free(lifted_atom);
						return axiom;
					}
					core::free(*lifted_atom); core::free(lifted_atom);
					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) return NULL;
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					core::free(*formula); if (formula->reference_count == 0) core::free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated, ResolveInconsistencies>(*atom, new_axiom, std::forward<Args>(args)...)) {
						core::free(*new_axiom); core::free(new_axiom);
						return NULL;
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
					if (arg1->type != TermType::NUMBER || arg2->type != TermType::NUMBER)
						return nullptr;
					if (arg2->number <= arg1->number) {
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
						core::free(*negated); if (negated->reference_count == 0) core::free(negated);
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
						if (sets.sets[i].size_axioms.data == nullptr) continue;
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
							for (Proof* edge_axiom : entry.value) {
								for (Proof* child : edge_axiom->children) {
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
					core::free(*formula); if (formula->reference_count == 0) core::free(formula);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
					if (!add_binary_atom<Negated, ResolveInconsistencies>({atom->ternary.first->constant, arg1->constant, arg2->constant}, new_axiom, std::forward<Args>(args)...)) {
						core::free(*new_axiom); core::free(new_axiom);
						return NULL;
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
	void free_proof(Proof* proof, set_changes<Formula>& set_diff, Args&&... args) {
		theory::changes changes;
		if (!get_theory_changes(*proof, changes, std::forward<Args>(args)...)) return;
		subtract_changes(changes, set_diff, std::forward<Args>(args)...);
		core::free(*proof); if (proof->reference_count == 0) core::free(proof);
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
				core::free(*atom); core::free(atom);
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
			core::free(conjunct_proofs);
			return proof != NULL;
		} else {
			return make_conjunct_proof<false>(constraint, c, proof);
		}
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

inline unsigned int index_of_number(const array<instance>& constants, hol_number number) {
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::NUMBER && constants[i].number == number)
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

inline bool contains_number(const array<instance>& constants, hol_number number) {
	return index_of_number(constants, number) < constants.length;
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
			} else if (right->type == TermType::NUMBER) {
				hol_number number = right->number;
				if (!contains_number(constants, number) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::NUMBER;
				constants[0].number = right->number;
				constants.length = 1;
			} else if (right->type == TermType::STRING) {
				string* str = &right->str;
				if (!contains_string(constants, str) && !contains_any(constants))
					return false;
				constants[0].type = instance_type::STRING;
				constants[0].str = &right->str;
				constants.length = 1;
			} else if (right->type == TermType::VARIABLE) {
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
					/* make sure `constants[i].constant` is not the arg1 of a name event */
					else if (constants[i].constant >= T.new_constant_offset) {
						for (Proof* proof : T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions) {
							if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
							 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
							 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
							 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
							 && T.is_name_event(proof->formula->binary.right->binary.right->constant))
							{
								constants.remove(i--);
								break;
							}
						}
					}
				}

				/* check to make sure the args are type-correct */
				if (right->type == TermType::UNARY_APPLICATION
				 && right->binary.left->type == TermType::CONSTANT
				 && right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				 && right->binary.right->type == TermType::CONSTANT)
				{
					/* make sure both `right->binary.right->constant` and `left->constant` are not name event */
					if (T.is_name_event(right->binary.right->constant)) {
						for (unsigned int i = 0; i < constants.length; i++) {
							if (constants[i].type == instance_type::ANY || constants[i].type != instance_type::CONSTANT)
								continue;
							if (T.is_name_event(constants[i].constant))
								constants.remove(i--);
						}
					}
				}
				if (right->type == TermType::UNARY_APPLICATION
				 && right->binary.left->type == TermType::CONSTANT
				 && right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
				 && right->binary.right->type == TermType::CONSTANT)
				{
					/* make sure the event is not a `name` since the arg2 of `name` events must be a string */
					if (T.is_name_event(right->binary.right->constant))
						return false;
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
			if (right->type == FormulaType::LAMBDA) {
				if (T.is_provably_not_a_set(left->constant))
					return false;
				/* make sure `left->constant` is not the arg1 of a name event */
				if (left->constant >= T.new_constant_offset) {
					for (Proof* proof : T.ground_concepts[left->constant - T.new_constant_offset].definitions) {
						if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
						 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
						 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
						 && T.is_name_event(proof->formula->binary.right->binary.right->constant))
						{
							return false;
						}
					}
				}
			}

			bool left_is_name_event = T.is_name_event(left->constant);

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

					/* check to make sure the args are type-correct */
					if (new_right->type == TermType::UNARY_APPLICATION
					 && new_right->binary.left->type == TermType::CONSTANT
					 && new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					 && new_right->binary.right->type == TermType::CONSTANT)
					{
						/* make sure both `right->binary.right->constant` and `left->constant` are not name event */
						if (left_is_name_event && T.is_name_event(new_right->binary.right->constant)) {
							free(*new_right); if (new_right->reference_count == 0) free(new_right);
							constants.remove(i--); continue;
						}
					}
					if (new_right->type == TermType::UNARY_APPLICATION
					 && new_right->binary.left->type == TermType::CONSTANT
					 && new_right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					 && new_right->binary.right->type == TermType::CONSTANT)
					{
						/* make sure the event is not a `name` since the arg2 of `name` events must be a string */
						if (T.is_name_event(new_right->binary.right->constant)) {
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

				/* check to make sure the args are type-correct */
				if (new_left->binary.right->constant >= T.new_constant_offset) {
					unsigned int event = new_left->binary.right->constant;
					if (new_left->binary.left->constant == (unsigned int) built_in_predicates::ARG1) {
						/* the event must not be `name` since `right` cannot be a constant here */
						if (T.is_name_event(event)) {
							free(*new_left); if (new_left->reference_count == 0) free(new_left);
							constants.remove(i--); continue;
						}
					} else if (new_left->binary.left->constant == (unsigned int) built_in_predicates::ARG2) {
						if (right->type == TermType::STRING) {
							/* the event must be `name` */
							if (T.ground_concepts[event - T.new_constant_offset].types.size > 1) {
								free(*new_left); if (new_left->reference_count == 0) free(new_left);
								constants.remove(i--); continue;
							}
							if (T.ground_concepts[event - T.new_constant_offset].types.size == 1) {
								Term& key = T.ground_concepts[event - T.new_constant_offset].types.keys[0];
								if (key.type != TermType::UNARY_APPLICATION
								 || key.binary.left->type != TermType::CONSTANT
								 || key.binary.left->constant != (unsigned int) built_in_predicates::NAME
								 || key.binary.right->type != TermType::VARIABLE
								 || key.binary.right->variable != 1)
								{
									free(*new_left); if (new_left->reference_count == 0) free(new_left);
									constants.remove(i--); continue;
								}
							}
						} else {
							/* the event must not be `name` */
							if (T.is_name_event(event)) {
								free(*new_left); if (new_left->reference_count == 0) free(new_left);
								constants.remove(i--); continue;
							}
						}
					}
				}

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
				/* make sure `constants[i].constant` is not the arg1 of a name event */
				else if (constants[i].constant >= T.new_constant_offset) {
					for (Proof* proof : T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions) {
						if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
						 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
						 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
						 && T.is_name_event(proof->formula->binary.right->binary.right->constant))
						{
							constants.remove(i--);
							break;
						}
					}
				}
			}

			/* check to make sure the args are type-correct */
			if (right->type == TermType::CONSTANT) {
				/* if arg1 is a `name` constant, `left` cannot be `name` */
				Term* arg1 = T.get_arg(right->constant, (unsigned int) built_in_predicates::ARG1);
				if (arg1 != nullptr && arg1->type == TermType::CONSTANT && T.is_name_event(arg1->constant)) {
					unsigned int index = index_of_constant(constants, (unsigned int) built_in_predicates::NAME);
					if (index < constants.length) constants.remove(index);
				}

				/* if arg2 is a string, the event must be `name` */
				Term* arg2 = T.get_arg(right->constant, (unsigned int) built_in_predicates::ARG2);
				if (arg2 != nullptr) {
					if (arg2->type == TermType::STRING) {
						if (contains_constant(constants, (unsigned int) built_in_predicates::NAME)) {
							constants[0] = instance_constant((unsigned int) built_in_predicates::NAME);
							constants.length = 1;
						} else {
							return false;
						}
					} else if (arg2->type != TermType::VARIABLE) {
						/* arg1 is not a string, so `left` cannot be `name` */
						unsigned int index = index_of_constant(constants, (unsigned int) built_in_predicates::NAME);
						if (index < constants.length) constants.remove(index);
					}
				}

				/* if this event is the arg1 of a name event, this event cannot be a name event */
				for (const Proof* definition : T.ground_concepts[right->constant - T.new_constant_offset].definitions) {
					Formula* right = definition->formula->binary.right;
					if (right->type == TermType::UNARY_APPLICATION
					 && right->binary.left->type == TermType::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					 && right->binary.right->type == TermType::CONSTANT && T.is_name_event(right->binary.right->constant))
					{
						unsigned int index = index_of_constant(constants, (unsigned int) built_in_predicates::NAME);
						if (index < constants.length) constants.remove(index);
						break;
					}
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

			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) continue;
				if (constants[i].type != instance_type::CONSTANT) {
					constants.remove(i--);
					continue;
				}
				/* check to make sure the args are type-correct */
				unsigned int arg = constants[i].constant;
				if (left->constant == (unsigned int) built_in_predicates::NAME) {
					/* arg1 must be a non-name constant and arg2 must be a string */
					Term* arg1 = T.get_arg(arg, (unsigned int) built_in_predicates::ARG1);
					if (arg1 != nullptr) {
						if (arg1->type != TermType::CONSTANT || T.is_name_event(arg1->constant)) {
							constants.remove(i--);
							continue;
						}
					}
					Term* arg2 = T.get_arg(arg, (unsigned int) built_in_predicates::ARG2);
					if (arg2 != nullptr && arg2->type != TermType::STRING) {
						constants.remove(i--);
						continue;
					}
					/* this event also cannot be the arg1 of a name event */
					bool is_arg1_of_name_event = false;
					for (const Proof* definition : T.ground_concepts[arg - T.new_constant_offset].definitions) {
						Formula* right = definition->formula->binary.right;
						if (right->type == TermType::UNARY_APPLICATION
						 && right->binary.left->type == TermType::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && right->binary.right->type == TermType::CONSTANT && T.is_name_event(right->binary.right->constant))
						{
							is_arg1_of_name_event = true;
							break;
						}
					}
					if (is_arg1_of_name_event) {
						constants.remove(i--);
						continue;
					}
				} else {
					/* arg1 can be anything and arg2 must be a non-string */
					bool contains;
					Proof* proof = T.ground_concepts[arg - T.new_constant_offset].function_values.get((unsigned int) built_in_predicates::ARG2, contains);
					if (contains && proof->formula->binary.right->type == TermType::STRING) {
						constants.remove(i--);
						continue;
					}
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
	typedef typename Proof::FormulaType Formula;

	Proof** proofs;
	unsigned int proof_count;
	Formula** extra_axioms;
	unsigned int extra_axiom_count;
	double log_probability;
#if !defined(NDEBUG)
	/* a unique identifier for debugging */
	unsigned int id;
#endif

	static inline unsigned int hash(const theory_sample<Proof>& key) {
		unsigned int hash_value = 0;
		for (unsigned int i = 0; i < key.proof_count; i++)
			hash_value ^= Proof::hash(*key.proofs[i]);
		for (unsigned int i = 0; i < key.extra_axiom_count; i++)
			hash_value ^= Formula::hash(*key.extra_axioms[i]);
		return hash_value;
	}

	static inline bool is_empty(const theory_sample<Proof>& key) {
		return key.proofs == nullptr;
	}

	static inline void move(const theory_sample<Proof>& src, theory_sample<Proof>& dst) {
		dst.proofs = src.proofs;
		dst.proof_count = src.proof_count;
		dst.extra_axioms = src.extra_axioms;
		dst.extra_axiom_count = src.extra_axiom_count;
		dst.log_probability = src.log_probability;
#if !defined(NDEBUG)
		dst.id = src.id;
#endif
	}

	static inline void free(theory_sample<Proof>& sample) { sample.free(); }

private:
	inline void free() {
		for (unsigned int i = 0; i < proof_count; i++) { core::free(*proofs[i]); if (proofs[i]->reference_count == 0) core::free(proofs[i]); }
		for (unsigned int i = 0; i < extra_axiom_count; i++) { core::free(*extra_axioms[i]); if (extra_axioms[i]->reference_count == 0) core::free(extra_axioms[i]); }
		core::free(proofs); core::free(extra_axioms);
	}
};

template<typename Proof>
bool init(theory_sample<Proof>& sample, const hash_set<Proof*>& proofs,
		const array<typename Proof::FormulaType*>& extra_axioms, double log_probability)
{
	typedef typename Proof::FormulaType Formula;

	sample.log_probability = log_probability;
	sample.proofs = (Proof**) malloc(max((size_t) 1, sizeof(Proof*) * proofs.size));
	if (sample.proofs == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `theory_sample.proofs`.\n");
		return false;
	}
	sample.extra_axioms = (Formula**) malloc(max((size_t) 1, sizeof(Proof*) * extra_axioms.length));
	if (sample.extra_axioms == nullptr) {
		fprintf(stderr, "init ERROR: Insufficient memory for `theory_sample.extra_axioms`.\n");
		free(sample.proofs); return false;
	}

	sample.proof_count = 0;
	array_map<const Proof*, Proof*> proof_map(32);
	for (Proof* proof : proofs) {
		if (!Proof::clone(proof, sample.proofs[sample.proof_count], proof_map)) {
			for (unsigned int j = 0; j < sample.proof_count; j++) { free(*sample.proofs[j]); free(sample.proofs[j]); }
			free(sample.proofs);
			return false;
		}
		sample.proof_count++;
	}

	sample.extra_axiom_count = 0;
	for (Formula* extra_axiom : extra_axioms) {
		sample.extra_axioms[sample.extra_axiom_count] = extra_axiom;
		extra_axiom->reference_count++;
		sample.extra_axiom_count++;
	}

	/* sort the proofs and extra axioms in canonical order */
	if (sample.proof_count > 1)
		sort(sample.proofs, sample.proof_count, pointer_sorter());
	if (sample.extra_axiom_count > 1)
		sort(sample.extra_axioms, sample.extra_axiom_count, pointer_sorter());
	return true;
}

template<typename Proof>
inline bool operator == (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr
	 || first.proof_count != second.proof_count
	 || first.extra_axiom_count != second.extra_axiom_count)
		return false;
	for (unsigned int i = 0; i < first.proof_count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return false;
	} for (unsigned int i = 0; i < first.extra_axiom_count; i++) {
		if (*first.extra_axioms[i] != *second.extra_axioms[i])
			return false;
	}
	return true;
}

template<typename Proof>
inline bool operator != (const theory_sample<Proof>& first, const theory_sample<Proof>& second)
{
	if (first.proofs == nullptr
	 || first.proof_count != second.proof_count
	 || first.extra_axiom_count != second.extra_axiom_count)
		return true;
	for (unsigned int i = 0; i < first.proof_count; i++) {
		if (*first.proofs[i] != *second.proofs[i])
			return true;
	} for (unsigned int i = 0; i < first.extra_axiom_count; i++) {
		if (*first.extra_axioms[i] != *second.extra_axioms[i])
			return true;
	}
	return false;
}

template<typename Proof>
inline double log_probability(const theory_sample<Proof>& sample) {
	return sample.log_probability;
}

struct null_collector {
	template<typename Proof>
	constexpr inline bool has_prior(const Proof* proof) const {
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept(const hash_set<Proof*>& sample,
			const array<typename Proof::FormulaType*>& extra_axioms,
			double proof_prior_diff) const
	{
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept_with_observation_changes(const hash_set<Proof*>& sample,
			const array<typename Proof::FormulaType*>& extra_axioms, double proof_prior_diff,
			const array<pair<Proof*, Proof*>>& observation_changes) const
	{
		return true;
	}
};

template<typename ProofCalculus, typename Canonicalizer>
struct log_probability_collector
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;

	double current_log_probability;

#if !defined(NDEBUG)
	std::function<double(void)> compute_current_log_probability;
#endif

	template<typename ProofPrior>
	log_probability_collector(const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior)
	{
		/* initialize `current_log_probability` */
		array<Formula*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		null_collector collector;
		current_log_probability = log_probability(T.observations, extra_axioms, proof_prior, collector);
#if !defined(NDEBUG)
		compute_current_log_probability = [&]() {
			null_collector collector;
			array<Formula*> extra_axioms(16);
			T.get_extra_axioms(extra_axioms);
			return log_probability(T.observations, extra_axioms, proof_prior, collector);
		};
#endif
	}

	constexpr inline bool has_prior(const Proof* proof) const {
		return true;
	}

	bool accept(const hash_set<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff)
	{
		current_log_probability += proof_prior_diff;

#if !defined(NDEBUG)
		double expected_log_probability = compute_current_log_probability();
		if (fabs(expected_log_probability - current_log_probability) > 1.0e-8) {
			fprintf(stderr, "log_probability_collector WARNING: The computed"
					" log probability of the sample (%lf) differs from the expected log probability (%lf).\n",
					current_log_probability, expected_log_probability);
		}
#endif
		return true;
	}

	inline bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff,
			const array<pair<Proof*, Proof*>>& observation_changes)
	{
		return accept(sample, extra_axioms, proof_prior_diff);
	}
};

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
inline log_probability_collector<ProofCalculus, Canonicalizer> make_log_probability_collector(const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior)
{
	return log_probability_collector<ProofCalculus, Canonicalizer>(T, proof_prior);
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
		array<Formula*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
		current_log_probability = log_probability(T.observations, extra_axioms, proof_prior, sample_collector);
#if !defined(NDEBUG)
		compute_current_log_probability = [&]() {
			array<Formula*> extra_axioms(16);
			T.get_extra_axioms(extra_axioms);
			return log_probability(T.observations, extra_axioms, proof_prior, sample_collector);
		};
#endif

		/* add the first sample */
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!init(new_sample, T.observations, extra_axioms, current_log_probability))
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

	bool accept(const hash_set<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff)
	{
		typedef typename ProofCalculus::Proof Proof;

		current_log_probability += proof_prior_diff;
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!samples.check_size()
		 || !init(new_sample, sample, extra_axioms, current_log_probability))
			return false;

		bool contains;
		unsigned int bucket = samples.index_of(new_sample, contains);
		if (contains) {
			/* we've already seen this sample before */
#if !defined(NDEBUG)
			if (fabs(new_sample.log_probability - samples.keys[bucket].log_probability) > 1.0e-8)
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
		if (fabs(expected_log_probability - new_sample.log_probability) > 1.0e-8) {
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
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff,
			const array<pair<Proof*, Proof*>>& observation_changes)
	{
		if (test_proof != nullptr) {
			for (unsigned int i = 0; i < observation_changes.length; i++) {
				if (observation_changes[i].key == test_proof) {
					test_proof = observation_changes[i].value;
					break;
				}
			}
		}

		return accept(sample, extra_axioms, proof_prior_diff);
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

	inline bool accept(const hash_set<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff)
	{
		return internal_collector.accept(sample, extra_axioms, proof_prior_diff);
	}

	bool accept_with_observation_changes(const hash_set<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff,
			const array<pair<Proof*, Proof*>>& observation_changes)
	{
		return internal_collector.accept_with_observation_changes(sample, extra_axioms, proof_prior_diff, observation_changes);
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
	} else if (!proof_axioms.template add<false>(new_proof, proof_prior)) {
		T.template remove_formula<true>(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
{
fprintf(stderr, "DEBUG: t = %u\n", t);
proof_axioms.check_proof_axioms(T);
proof_axioms.check_universal_eliminations(T, collector);
T.check_concept_axioms();
T.check_disjunction_introductions();
T.are_elements_provable();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
T.sets.are_set_sizes_valid();
T.sets.check_set_ids();
bool print_debug = false;
if (print_debug) T.print_axioms(stderr, *debug_terminal_printer);
if (print_debug) T.print_existential_introductions(stderr, *debug_terminal_printer);
		do_mh_step(T, proof_prior, proof_axioms, collector);
}

	T.template remove_formula<false>(collector.test_proof);
	proof_axioms.template subtract<false>(collector.test_proof, proof_prior);
	free(*collector.test_proof);
	if (collector.test_proof->reference_count == 0)
		free(collector.test_proof);
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
	} else if (!proof_axioms.template add<false>(new_proof, proof_prior)) {
		T.template remove_formula<true>(new_proof);
		return -std::numeric_limits<double>::infinity();
	}

	provability_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
		do_mh_step(T, proof_prior, proof_axioms, collector);

	T.template remove_formula<false>(collector.internal_collector.test_proof);
	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, proof_prior);
	free(*collector.internal_collector.test_proof);
	if (collector.internal_collector.test_proof->reference_count == 0)
		free(collector.internal_collector.test_proof);
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

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction, typename... Args>
bool log_joint_probability_of_lambda(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		theory<ProofCalculus, Canonicalizer>& T_map,
		OnProofSampleFunction on_new_proof_sample, Args&&... add_formula_args)
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
	set_changes<Formula> set_diff;
extern const string_map_scribe* debug_terminal_printer;
T.print_axioms(stderr, *debug_terminal_printer);
	Proof* new_proof = T.add_formula(existential, set_diff, new_constant, std::forward<Args>(add_formula_args)...);
	free(*existential); if (existential->reference_count == 0) free(existential);
	if (new_proof == nullptr) {
		return false;
	} else if (!proof_axioms.template add<false>(new_proof, set_diff.new_set_axioms, proof_prior)) {
		T.template remove_formula<true>(new_proof, set_diff);
		return false;
	}
	set_diff.clear();

	if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_map)) {
		T.template remove_formula<false>(new_proof, set_diff);
		proof_axioms.template subtract<false>(new_proof, set_diff.old_set_axioms, proof_prior);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		return false;
	}

	auto new_proof_sample_delegate = make_lambda_proof_sample_delegate<typename ProofCalculus::ProofCanonicalizer>(on_new_proof_sample);
	auto collector = make_provability_collector(T, proof_prior, new_proof, new_proof_sample_delegate);
	double max_log_probability = collector.internal_collector.current_log_probability;
	for (unsigned int t = 0; t < num_samples; t++)
{
fprintf(stderr, "DEBUG: t = %u\n", t);
proof_axioms.check_proof_axioms(T);
proof_axioms.check_universal_eliminations(T, collector);
T.check_concept_axioms();
T.check_disjunction_introductions();
T.are_elements_provable();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
T.sets.are_set_sizes_valid();
T.sets.check_set_ids();
if (!T.observations.contains(collector.internal_collector.test_proof))
	fprintf(stderr, "log_joint_probability_of_lambda WARNING: `provability_collector.internal_collector.test_proof` is not an observation in the theory.\n");
bool print_debug = false;
if (print_debug) T.print_axioms(stderr, *debug_terminal_printer);
if (print_debug) T.print_existential_introductions(stderr, *debug_terminal_printer);
		do_mh_step(T, proof_prior, proof_axioms, collector, collector.internal_collector.test_proof);
		if (collector.internal_collector.current_log_probability > max_log_probability) {
			free(T_map);
			if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_map)) {
				T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
				proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
				free(*collector.internal_collector.test_proof);
				if (collector.internal_collector.test_proof->reference_count == 0)
					free(collector.internal_collector.test_proof);
				return false;
			}
			max_log_probability = collector.internal_collector.current_log_probability;
		}
if (t % 1000 == 0)
	fprintf(stdout, "(seed = %u)\n", get_seed());
}
	T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
	free(*collector.internal_collector.test_proof);
	if (collector.internal_collector.test_proof->reference_count == 0)
		free(collector.internal_collector.test_proof);
	return true;
}

#endif /* THEORY_H_ */
