#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>
#include <core/random.h>
#include <math/multiset.h>
#include <stdexcept>

#if !defined(NDEBUG)
#include <functional>
#endif

#include "array_view.h"
#include "set_reasoning.h"
#include "built_in_predicates.h"
#include "lf_utils.h"
#include "context.h"

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

template<typename Theory, typename Formula, typename Proof>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, const typename Formula::Term* term, unsigned int variable, const array_map<unsigned int, Proof*>& set_definitions) { return true; }


typedef uint_fast8_t instance_type_specifier;
enum class instance_type : instance_type_specifier {
	ANY,
	CONSTANT,
	NUMBER,
	STRING
};

/* TODO: for debugging; delete this */
#include <atomic>
std::atomic<unsigned long long> total_reasoning(0);
std::atomic<unsigned long long> consistency_checking_ms(0);
thread_local bool consistency_checking = false;
struct time_aggregator {
	std::atomic<unsigned long long>& milliseconds;
	bool& guard;
	bool old_guard;
	timer stopwatch;

	time_aggregator(std::atomic<unsigned long long>& milliseconds, bool& guard) :
			milliseconds(milliseconds), guard(guard), old_guard(guard)
	{
		guard = true;
	}
	~time_aggregator() {
		if (!old_guard) {
			milliseconds += stopwatch.milliseconds();
			guard = false;
		}
	}
};

struct instance {
	instance_type type;
	union {
		unsigned int constant;
		hol_number number;
		string* str;
	};
	unsigned int matching_types;
	unsigned int mismatching_types;

	static inline void swap(instance& first, instance& second) {
		char* first_data = (char*) &first;
		char* second_data = (char*) &second;
		for (unsigned int i = 0; i < sizeof(instance); i++)
			core::swap(first_data[i], second_data[i]);
	}

	static inline void move(const instance& src, instance& dst) {
		dst.type = src.type;
		dst.matching_types = src.matching_types;
		dst.mismatching_types = src.mismatching_types;
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

inline bool operator < (const instance& first, const instance& second) {
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case instance_type::ANY: return false;
	case instance_type::CONSTANT: return first.constant < second.constant;
	case instance_type::NUMBER: return first.number < second.number;
	case instance_type::STRING: return *first.str < *second.str;
	}
	fprintf(stderr, "operator < ERROR: Unrecognized `instance_type`.\n");
	exit(EXIT_FAILURE);
}

inline bool operator >= (const instance& first, const instance& second) {
	return !(first < second);
}

template<typename Stream>
inline bool read(instance& inst, Stream& in) {
	instance_type_specifier type;
	if (!read(type, in)
	 || !read(inst.matching_types, in)
	 || !read(inst.mismatching_types, in))
	{
		return false;
	}
	inst.type = (instance_type) type;

	switch (inst.type) {
	case instance_type::ANY:
		return true;
	case instance_type::CONSTANT:
		return read(inst.constant, in);
	case instance_type::NUMBER:
		return read(inst.number, in);
	case instance_type::STRING:
		fprintf(stderr, "read ERROR: Serialization/deserialization of `instance` with type `STRING` is unsupported.");
		return false;
	}
	fprintf(stderr, "read ERROR: Unrecognized `instance_type`.");
	return false;
}

template<typename Stream>
inline bool write(const instance& inst, Stream& out) {
	if (!write((instance_type_specifier) inst.type, out)
	 || !write(inst.matching_types, out)
	 || !write(inst.mismatching_types, out))
	{
		return false;
	}

	switch (inst.type) {
	case instance_type::ANY:
		return true;
	case instance_type::CONSTANT:
		return write(inst.constant, out);
	case instance_type::NUMBER:
		return write(inst.number, out);
	case instance_type::STRING:
		fprintf(stderr, "write ERROR: Serialization/deserialization of `instance` with type `STRING` is unsupported.");
		return false;
	}
	fprintf(stderr, "write ERROR: Unrecognized `instance_type`.");
	return false;
}

inline instance instance_any() {
	instance i;
	i.type = instance_type::ANY;
	i.matching_types = 0;
	i.mismatching_types = 0;
	return i;
}

inline instance instance_constant(unsigned int constant) {
	instance i;
	i.type = instance_type::CONSTANT;
	i.constant = constant;
	i.matching_types = 0;
	i.mismatching_types = 0;
	return i;
}

inline instance instance_number(hol_number number) {
	instance i;
	i.type = instance_type::NUMBER;
	i.number = number;
	i.matching_types = 0;
	i.mismatching_types = 0;
	return i;
}

inline instance instance_number(int64_t integer, uint64_t decimal) {
	instance i;
	i.type = instance_type::NUMBER;
	i.number.integer = integer;
	i.number.decimal = decimal;
	i.matching_types = 0;
	i.mismatching_types = 0;
	return i;
}

inline instance instance_string(string* str) {
	instance i;
	i.type = instance_type::STRING;
	i.str = str;
	i.matching_types = 0;
	i.mismatching_types = 0;
	return i;
}

inline size_t index_of_any(const array<instance>& constants) {
	for (size_t i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::ANY)
			return i;
	}
	return constants.length;
}

inline size_t index_of_constant(const array<instance>& constants, unsigned int constant) {
	for (size_t i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::CONSTANT && constants[i].constant == constant)
			return i;
	}
	return constants.length;
}

inline size_t index_of_number(const array<instance>& constants, hol_number number) {
	for (size_t i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::NUMBER && constants[i].number == number)
			return i;
	}
	return constants.length;
}

inline size_t index_of_string(const array<instance>& constants, string* str) {
	for (size_t i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::STRING && (constants[i].str == str || *constants[i].str == *str))
			return i;
	}
	return constants.length;
}

inline bool contains_any(const array<instance>& constants) {
	return index_of_any(constants) < constants.length;
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

template<typename T, typename Stream>
inline bool print(const interval<T>& i, Stream& out) {
	if (i.left_inclusive) {
		if (fputc('[', out) == EOF) return false;
	} else {
		if (fputc('(', out) == EOF) return false;
	}
	if (!print(i.min, out)
	 || fputc(',', out) == EOF
	 || !print(i.max, out))
		return false;
	if (i.left_inclusive) {
		if (fputc(']', out) == EOF) return false;
	} else {
		if (fputc(')', out) == EOF) return false;
	}
	return true;
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

	instantiation(instantiation_type type) {
		if (!init_helper(type)) throw std::bad_alloc();
	}

	instantiation(const instantiation& src) {
		if (!init_helper(src)) throw std::bad_alloc();
	}

	~instantiation() { free(*this); }

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
	inline bool init_helper(instantiation_type inst_type) {
		type = inst_type;
		if (type == instantiation_type::ANY) {
			any.excluded = (instantiation*) malloc(1);
			if (any.excluded == nullptr)
				return false;
			any.excluded_count = 0;
		}
		return true;
	}

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

	friend bool init(instantiation&, instantiation_type);
	friend bool init(instantiation&, const instantiation&);
	friend bool operator != (const instantiation&, const instantiation&);
	friend bool operator < (const instantiation&, const instantiation&);
};

inline bool init(instantiation& inst, instantiation_type type) {
	return inst.init_helper(type);
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

template<typename Stream, typename... Printer>
bool print(const instantiation& inst, Stream& out, Printer&&... printer) {
	switch (inst.type) {
	case instantiation_type::ANY:
		if (fputc('*', out) == EOF) return false;
		if (inst.any.excluded_count != 0
		 && fputc('\\', out) == EOF) return false;
		if (inst.any.excluded_count > 1
		 && fputc('(', out) == EOF) return false;
		for (uint_fast8_t i = 0; i < inst.any.excluded_count; i++) {
			if (i != 0 && fprintf(out, " ⋃ ") <= 0) return false;
			if (!print(inst.any.excluded[i], out, std::forward<Printer>(printer)...))
				return false;
		}
		if (inst.any.excluded_count > 1
		 && fputc(')', out) == EOF) return false;
	case instantiation_type::ANY_NUMBER:
		for (uint_fast8_t i = 0; i < inst.any_number.included_count; i++) {
			if (i != 0 && fprintf(out, " ⋃ ") <= 0) return false;
			if (!print(inst.any_number.included[i], out)) return false;
		}
		return true;
	case instantiation_type::CONSTANT:
		return print(inst.constant, out, std::forward<Printer>(printer)...);
	case instantiation_type::NUMBER:
		return print(inst.number, out);
	case instantiation_type::STRING:
		return fputc('"', out) != EOF
			&& print(inst.str, out)
			&& fputc('"', out) != EOF;
	}
	fprintf(stderr, "print < ERROR: Unrecognized `instantiation_type`.\n");
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

inline bool is_subset(
		const any_number_instantiation& first,
		const any_number_instantiation& second)
{
	uint_fast8_t i = 0, j = 0;
	hol_number lower;
	bool left_inclusive;
	if (first.included[0].min < second.included[0].min || (second.included[0].min == first.included[0].min && first.included[0].left_inclusive)) {
		lower = first.included[0].min;
		left_inclusive = first.included[0].left_inclusive;
	} else {
		return false;
	}
	while (i < first.included_count && j < second.included_count) {
		if (!has_intersection(lower, left_inclusive, first.included[i].max, first.included[i].right_inclusive)) {
			i++;
		} else if (!has_intersection(lower, left_inclusive, second.included[j].max, second.included[j].right_inclusive)) {
			j++;
		} else {
			if (first.included[i].max < second.included[j].max || (second.included[j].max == first.included[i].max && !first.included[i].right_inclusive))
				return false;
			i++; j++;

			if (i == first.included_count || j == second.included_count)
				return true;
			if (first.included[i].min < second.included[j].min || (second.included[j].min == first.included[i].min && first.included[i].left_inclusive)) {
				lower = first.included[i].min;
				left_inclusive = first.included[i].left_inclusive;
			} else {
				return false;
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
		if  (!array_init(new_excluded, max(1, first.any.excluded_count + second.any.excluded_count)))
			return false;
		set_union(new_excluded.data, new_excluded.length, first.any.excluded, first.any.excluded_count, second.any.excluded, second.any.excluded_count);
		dst.type = instantiation_type::ANY;
		dst.any.excluded = new_excluded.data;
		dst.any.excluded_count = (uint_fast8_t) new_excluded.length;
		return true;
	} else if (second.type == instantiation_type::ANY_NUMBER) {
		if (first.any.excluded_count == 0)
			return init(dst, second);
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

inline bool is_subset(const instantiation& first, const instantiation& second)
{
	if (first.type == instantiation_type::ANY || second.type == instantiation_type::ANY) {
		instantiation& dst = *((instantiation*) alloca(sizeof(instantiation)));
		intersect(dst, first, second);
		bool result = (dst == first);
		free(dst);
		return result;
	} else if (first.type == instantiation_type::ANY_NUMBER) {
		if (second.type == instantiation_type::ANY_NUMBER) {
			return is_subset(first.any_number, second.any_number);
		} else {
			return false;
		}
	} else if (second.type == instantiation_type::ANY_NUMBER) {
		if (first.type == instantiation_type::NUMBER) {
			for (uint_fast8_t i = 0; i < second.any_number.included_count; i++)
				if (second.any_number.included[i].contains(first.number)) return true;
			return false;
		} else {
			return false;
		}
	} else if (first.type != second.type) {
		return false;
	}

	switch (first.type) {
	case instantiation_type::CONSTANT:
		return (first.constant == second.constant);
	case instantiation_type::NUMBER:
		return (first.number == second.number);
	case instantiation_type::STRING:
		return (first.str == second.str);
	case instantiation_type::ANY:
	case instantiation_type::ANY_NUMBER:
		break;
	}
	fprintf(stderr, "is_subset ERROR: Unrecognized `instantiation_type`.\n");
	return false;
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
		if (!init_helper(src, src.length)) throw std::bad_alloc();
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
		} else if (!init_helper(src, src.length)) {
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

	inline bool unify_indices(uint_fast8_t index, uint_fast8_t other_index) {
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (!intersect(dummy[1], values[index], values[other_index]))
			return false;

		if (dummy[1].type != instantiation_type::ANY && dummy[1].type != instantiation_type::ANY_NUMBER) {
			bool result = change_value(index, dummy[1]) && change_value(other_index, dummy[1]);
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
			if (entry.key == other_index - 1 || entry.value == other_index) {
				second_root = entry.key;
				break;
			}
		}
		if (first_root == length) first_root = index;
		if (second_root == length) second_root = other_index;
		if (first_root == second_root) {
			core::free(dummy[1]);
			return true;
		}

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
			if (!new_component[i] || i == index || i == other_index) continue;
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
			return unify_indices(index, term->variable - 1);
		} else {
			return false;
		}
	}

	inline bool antiunify_indices(uint_fast8_t index, uint_fast8_t other_index) {
		instantiation* dummy = (instantiation*) alloca(2 * sizeof(instantiation));
		if (values[index].type == instantiation_type::ANY || values[index].type == instantiation_type::ANY_NUMBER) {
			if (values[other_index].type == instantiation_type::ANY || values[other_index].type == instantiation_type::ANY_NUMBER) {
				uint_fast8_t first_root = length;
				uint_fast8_t second_root = length;
				for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
					if (entry.key == index || entry.value == index) {
						first_root = entry.key;
						break;
					}
				} for (const pair<uint_fast8_t, uint_fast8_t>& entry : equal_indices) {
					if (entry.key == other_index || entry.value == other_index) {
						second_root = entry.key;
						break;
					}
				}
				if (first_root == length) first_root = index;
				if (second_root == length) second_root = other_index;
				if (first_root == second_root) return false;

				if (!not_equal_indices.ensure_capacity(not_equal_indices.length + 1))
					return false;

				pair<uint_fast8_t, uint_fast8_t> new_pair = (first_root < second_root ? make_pair(first_root, second_root) : make_pair(second_root, first_root));
				uint_fast8_t index = 0;
				while (index < not_equal_indices.length && less_than(new_pair, not_equal_indices[index], pair_sorter()))
					index++;
				if (index < not_equal_indices.length && new_pair == not_equal_indices[index])
					return true;
				shift_right(not_equal_indices.data, (unsigned int) not_equal_indices.length, index);
				not_equal_indices[index] = new_pair;
				not_equal_indices.length++;
				return true;
			} else {
				if (!subtract(dummy[1], values[index], values[other_index]))
					return false;
				bool result = change_value(index, dummy[1]);
				core::free(dummy[1]); return result;
			}
		} else {
			if (values[other_index].type == instantiation_type::ANY) {
				if (!subtract(dummy[1], values[other_index], values[index]))
					return false;
				bool result = change_value(other_index, dummy[1]);
				core::free(dummy[1]); return result;
			} else {
				return values[index] != values[other_index];
			}
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
			return antiunify_indices(index, term->variable - 1);
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
			dummy[0].any_number.included[0].left_inclusive = true;
			dummy[0].any_number.included[0].max = hol_number::max();
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
			dummy[0].any_number.included[0].left_inclusive = true;
			if (values[first].type == instantiation_type::ANY_NUMBER) {
				dummy[0].any_number.included[0].max = values[first].any_number.max_interval().max;
				dummy[0].any_number.included[0].right_inclusive = values[first].any_number.max_interval().right_inclusive;
			} else if (values[first].type == instantiation_type::NUMBER) {
				dummy[0].any_number.included[0].max = values[first].number;
				dummy[0].any_number.included[0].right_inclusive = true;
			} else {
				dummy[0].any_number.included[0].max = hol_number::max();
				dummy[0].any_number.included[0].right_inclusive = true;
			}
			dummy[0].any_number.included_count = 1;
			if (!intersect(dummy[1], values[second], dummy[0]))
				return false;
			bool result = change_value(second, dummy[1]);
			core::free(dummy[1]);
			if (!result) return false;
		}
		dummy[0].type = instantiation_type::ANY_NUMBER;
		dummy[0].any_number.included = (interval<hol_number>*) alloca(sizeof(interval<hol_number>));
		if (values[second].type == instantiation_type::ANY_NUMBER) {
			dummy[0].any_number.included[0].min = values[second].any_number.min_interval().min;
			dummy[0].any_number.included[0].left_inclusive = values[second].any_number.min_interval().left_inclusive;
		} else if (values[second].type == instantiation_type::NUMBER) {
			dummy[0].any_number.included[0].min = values[second].number;
			dummy[0].any_number.included[0].left_inclusive = true;
		}
		dummy[0].any_number.included[0].max = hol_number::max();
		dummy[0].any_number.included[0].right_inclusive = true;
		dummy[0].any_number.included_count = 1;
		if (!intersect(dummy[1], values[first], dummy[0]))
			return false;
		bool result = change_value(first, dummy[1]);
		core::free(dummy[1]);
		if (!result) return false;

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
		if (first_root == length) first_root = first;
		if (second_root == length) second_root = second;
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
		values = (instantiation*) malloc(max(1, sizeof(instantiation) * length));
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

	inline bool init_helper(const instantiation_tuple& src, uint_fast8_t new_length) {
		length = new_length;
		values = (instantiation*) malloc(max(1, sizeof(instantiation) * length));
		if (values == nullptr) {
			fprintf(stderr, "instantiation_tuple.init_helper ERROR: Out of memory.\n");
			return false;
		}
		for (uint_fast8_t i = 0; i < min(src.length, new_length); i++) {
			if (!init(values[i], src.values[i])) {
				for (unsigned int j = 0; j < i; j++) core::free(values[j]);
				core::free(values);
				return false;
			}
		} for (uint_fast8_t i = src.length; i < new_length; i++) {
			if (!init(values[i], instantiation_type::ANY)) {
				for (unsigned int j = 0; j < i; j++) core::free(values[j]);
				core::free(values);
				return false;
			}
		}

		for (const pair<uint_fast8_t, uint_fast8_t>& index : src.equal_indices) {
			if (index.key < new_length && index.value < new_length && !equal_indices.add(index)) {
				free_helper();
				return false;
			}
		} for (const pair<uint_fast8_t, uint_fast8_t>& index : src.not_equal_indices) {
			if (index.key < new_length && index.value < new_length && !not_equal_indices.add(index)) {
				free_helper();
				return false;
			}
		} for (const pair<uint_fast8_t, uint_fast8_t>& index : src.ge_indices) {
			if (index.key < new_length && index.value < new_length && !ge_indices.add(index)) {
				free_helper();
				return false;
			}
		}
		return true;
	}

	inline void free_helper() {
		for (unsigned int i = 0; i < length; i++)
			core::free(values[i]);
		core::free(values);
	}

	friend bool init(instantiation_tuple&, uint_fast8_t);
	friend bool init(instantiation_tuple&, const instantiation_tuple&, uint_fast8_t);
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

inline bool init(instantiation_tuple& new_tuple, const instantiation_tuple& src, uint_fast8_t new_length) {
	if (!array_init(new_tuple.equal_indices, max((size_t) 1, src.equal_indices.length))) {
		return false;
	} else if (!array_init(new_tuple.not_equal_indices, max((size_t) 1, src.not_equal_indices.length))) {
		free(new_tuple.equal_indices);
		return false;
	} else if (!array_init(new_tuple.ge_indices, max((size_t) 1, src.ge_indices.length))) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		return false;
	} else if (!new_tuple.init_helper(src, new_length)) {
		free(new_tuple.equal_indices);
		free(new_tuple.not_equal_indices);
		free(new_tuple.ge_indices);
		return false;
	}
	return true;
}

inline bool init(instantiation_tuple& new_tuple, const instantiation_tuple& src) {
	return init(new_tuple, src, src.length);
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

inline bool is_subset(const instantiation_tuple& first, const instantiation_tuple& second)
{
	if (first.length != second.length) return false;
	for (uint_fast8_t i = 0; i < first.length; i++) {
		if (!is_subset(first.values[i], second.values[i]))
			return false;
	}
	return is_subset(second.equal_indices.data, (unsigned int) second.equal_indices.length, first.equal_indices.data, (unsigned int) first.equal_indices.length)
		&& is_subset(second.not_equal_indices.data, (unsigned int) second.not_equal_indices.length, first.not_equal_indices.data, (unsigned int) first.not_equal_indices.length)
		&& is_subset(second.ge_indices.data, (unsigned int) second.ge_indices.length, first.ge_indices.data, (unsigned int) first.ge_indices.length);
}

template<typename K, typename V, typename Stream, typename... Printer>
inline bool print(const pair<K, V>& p, Stream& out, Printer&&... printer) {
	return fputc('(', out) != EOF
		&& !print(p.key, out, std::forward<Printer>(printer)...)
		&& fputc(',', out) != EOF
		&& !print(p.value, out, std::forward<Printer>(printer)...)
		&& fputc(')', out) != EOF;
}

template<typename Stream, typename... Printer>
bool print(const instantiation_tuple& tuple, Stream& out, Printer&&... printer) {
	if (tuple.length == 0) {
		if (fprintf(out, "()") <= 0)
			return false;
	} else if (tuple.length == 1) {
		if (!print(tuple.values[0], out, std::forward<Printer>(printer)...))
			return false;
	} else {
		for (uint_fast8_t i = 0; i < tuple.length; i++) {
			if (i == 0) {
				if (fputc('(', out) == EOF)
					return false;
			} else {
				if (fputc(',', out) == EOF)
					return false;
			}
			if (!print(tuple.values[i], out, std::forward<Printer>(printer)...))
				return false;
		}
		if (fputc(')', out) == EOF)
			return false;
	}
	if (tuple.equal_indices.length != 0) {
		if (fprintf(stderr, " equal_indices: ") <= 0
		 || !print(tuple.equal_indices, out)) return false;
	} if (tuple.not_equal_indices.length != 0) {
		if (fprintf(stderr, " not_equal_indices: ") <= 0
		 || !print(tuple.not_equal_indices, out)) return false;
	} if (tuple.ge_indices.length != 0) {
		if (fprintf(stderr, " ge_indices: ") <= 0
		 || !print(tuple.ge_indices, out)) return false;
	}
	return true;
}

struct axiom_assignment {
	instantiation_tuple assignment;
	array_map<uint_fast8_t, uint_fast8_t> variable_map;

	static inline void free(axiom_assignment& assignment) {
		core::free(assignment.assignment);
		core::free(assignment.variable_map);
	}
};

inline bool init(axiom_assignment& assignment, const axiom_assignment& src) {
	if (!init(assignment.assignment, src.assignment)) {
		return false;
	} else if (!array_map_init(assignment.variable_map, src.variable_map.capacity)) {
		free(assignment.assignment);
		return false;
	}
	for (unsigned int i = 0; i < src.variable_map.size; i++) {
		assignment.variable_map.keys[i] = src.variable_map.keys[i];
		assignment.variable_map.values[i] = src.variable_map.values[i];
	}
	assignment.variable_map.size = src.variable_map.size;
	return true;
}

inline bool init(axiom_assignment& assignment, const instantiation_tuple& src_assignment) {
	if (!init(assignment.assignment, src_assignment)) {
		return false;
	} else if (!array_map_init(assignment.variable_map, 4)) {
		free(assignment.assignment);
		return false;
	}
	return true;
}

inline bool operator == (const axiom_assignment& src, const axiom_assignment& dst) {
	if (src.variable_map.size != dst.variable_map.size
	 || src.assignment != dst.assignment)
		return false;
	for (unsigned int i = 0; i < src.variable_map.size; i++) {
		if (src.variable_map.keys[i] != dst.variable_map.keys[i]
		 || src.variable_map.values[i] != dst.variable_map.values[i])
		{
			return false;
		}
	}
	return true;
}

inline bool operator != (const axiom_assignment& src, const axiom_assignment& dst) {
	return !(src == dst);
}

inline bool operator < (const axiom_assignment& src, const axiom_assignment& dst) {
	if (src.assignment < dst.assignment) return true;
	else if (dst.assignment < src.assignment) return false;

	if (src.variable_map.size < dst.variable_map.size) return true;
	else if (src.variable_map.size > dst.variable_map.size) return false;
	for (unsigned int i = 0; i < src.variable_map.size; i++) {
		if (src.variable_map.keys[i] < dst.variable_map.keys[i]) return true;
		else if (src.variable_map.keys[i] > dst.variable_map.keys[i]) return false;
		else if (src.variable_map.values[i] < dst.variable_map.values[i]) return true;
		else if (dst.variable_map.values[i] < src.variable_map.values[i]) return false;
	}
	return false;
}

inline bool is_subset(const axiom_assignment& first, const axiom_assignment& second)
{
	for (const auto& entry : second.variable_map) {
		size_t index = first.variable_map.index_of(entry.key);
		if (index == first.variable_map.size
		 || entry.value != first.variable_map.values[index])
		{
			return false;
		}
	}
	return is_subset(first.assignment, second.assignment);
}

struct variable_assignment {
	instantiation_tuple assignment;
	array_map<unsigned int, axiom_assignment> axiom_assignments;
	void* matching_axiom;

	variable_assignment(const variable_assignment& src) :
			assignment(src.assignment),
			axiom_assignments(src.axiom_assignments.capacity)
	{
		if (!init_helper(src)) throw std::bad_alloc();
	}

	~variable_assignment() { free_helper(); }

	inline bool operator = (const variable_assignment& src) {
		if (!init(assignment, src.assignment)) {
			return false;
		} else if (!array_map_init(axiom_assignments, src.axiom_assignments.capacity)) {
			core::free(assignment);
			return false;
		} else if (!init_helper(src)) {
			core::free(assignment);
			core::free(axiom_assignments);
			return false;
		}
		return true;
	}

	template<typename Term>
	inline bool unify_value(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		if (!assignment.unify_value(index, term))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (term->type == TermType::VARIABLE) {
				uint_fast8_t term_variable = current_assignment.variable_map.get(term->variable, contains);
				if (!contains) continue;
				if (!current_assignment.assignment.unify_indices(other_variable - 1, term_variable - 1))
					return false;
			} else if (!current_assignment.assignment.unify_value(other_variable - 1, term)) {
				return false;
			}
		}
		return true;
	}

	template<typename Term>
	inline bool antiunify_value(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		if (!assignment.antiunify_value(index, term))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (term->type == TermType::VARIABLE) {
				uint_fast8_t term_variable = current_assignment.variable_map.get(term->variable, contains);
				if (!contains) continue;
				if (!current_assignment.assignment.antiunify_indices(other_variable - 1, term_variable - 1))
					return false;
			} else if (!current_assignment.assignment.antiunify_value(other_variable - 1, term)) {
				return false;
			}
		}
		return true;
	}

	inline bool unify_number(uint_fast8_t index) {
		if (!assignment.unify_number(index))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (!current_assignment.assignment.unify_number(other_variable - 1))
				return false;
		}
		return true;
	}

	inline bool unify_not_number(uint_fast8_t index) {
		if (!assignment.unify_not_number(index))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (!current_assignment.assignment.unify_not_number(other_variable - 1))
				return false;
		}
		return true;
	}

	template<typename Term>
	inline bool unify_greater_than_or_equal(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		if (!assignment.unify_greater_than_or_equal(index, term))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (term->type == TermType::VARIABLE) {
				uint_fast8_t term_variable = current_assignment.variable_map.get(term->variable, contains);
				if (!contains) continue;
				if (!current_assignment.assignment.unify_greater_than_or_equal(other_variable - 1, term_variable - 1))
					return false;
			} else if (!current_assignment.assignment.unify_greater_than_or_equal(other_variable - 1, term)) {
				return false;
			}
		}
		return true;
	}

	template<typename Term>
	inline bool unify_less_than_or_equal(uint_fast8_t index, const Term* term) {
		typedef typename Term::Type TermType;
		if (!assignment.unify_less_than_or_equal(index, term))
			return false;
		for (unsigned int i = 0; i < axiom_assignments.size; i++) {
			axiom_assignment& current_assignment = axiom_assignments.values[i];
			bool contains;
			uint_fast8_t other_variable = current_assignment.variable_map.get(index + 1, contains);
			if (!contains) continue;
			if (term->type == TermType::VARIABLE) {
				uint_fast8_t term_variable = current_assignment.variable_map.get(term->variable, contains);
				if (!contains) continue;
				if (!current_assignment.assignment.unify_greater_than_or_equal(term_variable - 1, other_variable - 1))
					return false;
			} else if (!current_assignment.assignment.unify_less_than_or_equal(other_variable - 1, term)) {
				return false;
			}
		}
		return true;
	}

	static inline void swap(variable_assignment& first, variable_assignment& second) {
		core::swap(first.assignment, second.assignment);
		core::swap(first.axiom_assignments, second.axiom_assignments);
	}

	static inline void move(const variable_assignment& src, variable_assignment& dst) {
		core::move(src.assignment, dst.assignment);
		core::move(src.axiom_assignments, dst.axiom_assignments);
	}

	static inline void free(variable_assignment& assignment) {
		assignment.free_helper();
		core::free(assignment.assignment);
		core::free(assignment.axiom_assignments);
	}

private:
	inline bool init_helper(const variable_assignment& src) {
		matching_axiom = src.matching_axiom;
		for (unsigned int i = 0; i < src.axiom_assignments.size; i++) {
			axiom_assignments.keys[i] = src.axiom_assignments.keys[i];
			if (!init(axiom_assignments.values[i], src.axiom_assignments.values[i])) {
				free_helper();
				return false;
			}
			axiom_assignments.size++;
		}
		return true;
	}

	inline void free_helper() {
		for (auto entry : axiom_assignments)
			core::free(entry.value);
	}

	friend bool init(variable_assignment&, const variable_assignment&, uint_fast8_t);
	friend bool init(variable_assignment&, const instantiation_tuple&, const array_map<uint_fast8_t, uint_fast8_t>&);
};

inline bool init(variable_assignment& assignment, uint_fast8_t length) {
	if (!init(assignment.assignment, length)) {
		return false;
	} else if (!array_map_init(assignment.axiom_assignments, 2)) {
		free(assignment.assignment);
		return false;
	}
	assignment.matching_axiom = nullptr;
	return true;
}

inline bool init(variable_assignment& assignment, const variable_assignment& src, uint_fast8_t new_length) {
	if (!init(assignment.assignment, src.assignment, new_length)) {
		return false;
	} else if (!array_map_init(assignment.axiom_assignments, src.axiom_assignments.capacity)) {
		free(assignment.assignment);
		return false;
	} else if (!assignment.init_helper(src)) {
		free(assignment.assignment);
		free(assignment.axiom_assignments);
		return false;
	}
	return true;
}

inline bool init(variable_assignment& assignment, const variable_assignment& src) {
	return init(assignment, src, src.assignment.length);
}

inline bool init(variable_assignment& assignment, const instantiation_tuple& src_assignment) {
	if (!init(assignment.assignment, src_assignment)) {
		return false;
	} else if (!array_map_init(assignment.axiom_assignments, 2)) {
		free(assignment.assignment);
		return false;
	}
	assignment.matching_axiom = nullptr;
	return true;
}

inline bool operator == (const variable_assignment& src, const variable_assignment& dst) {
	if (src.axiom_assignments.size != dst.axiom_assignments.size
	 || src.assignment != dst.assignment)
		return false;
	for (unsigned int i = 0; i < src.axiom_assignments.size; i++) {
		if (src.axiom_assignments.keys[i] != dst.axiom_assignments.keys[i]
		 || src.axiom_assignments.values[i] != dst.axiom_assignments.values[i])
		{
			return false;
		}
	}
	return true;
}

inline bool operator != (const variable_assignment& src, const variable_assignment& dst) {
	return !(src == dst);
}

inline bool operator < (const variable_assignment& src, const variable_assignment& dst) {
	if (src.assignment < dst.assignment) return true;
	else if (dst.assignment < src.assignment) return false;

	if (src.axiom_assignments.size < dst.axiom_assignments.size) return true;
	else if (src.axiom_assignments.size > dst.axiom_assignments.size) return false;
	for (unsigned int i = 0; i < src.axiom_assignments.size; i++) {
		if (src.axiom_assignments.keys[i] < dst.axiom_assignments.keys[i]) return true;
		else if (src.axiom_assignments.keys[i] > dst.axiom_assignments.keys[i]) return false;
		else if (src.axiom_assignments.values[i] < dst.axiom_assignments.values[i]) return true;
		else if (dst.axiom_assignments.values[i] < src.axiom_assignments.values[i]) return false;
	}
	return false;
}

inline bool is_subset(const variable_assignment& first, const variable_assignment& second)
{
	for (const auto& entry : second.axiom_assignments) {
		size_t index = first.axiom_assignments.index_of(entry.key);
		if (index == first.axiom_assignments.size
		 || !is_subset(entry.value, first.axiom_assignments.values[index]))
		{
			return false;
		}
	}
	return is_subset(first.assignment, second.assignment);
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

template<typename Stream>
inline bool read(relation& rel, Stream& in) {
	return read(rel.predicate, in)
		&& read(rel.arg1, in)
		&& read(rel.arg2, in);
}

template<typename Stream>
inline bool write(const relation& rel, Stream& out) {
	return write(rel.predicate, out)
		&& write(rel.arg1, out)
		&& write(rel.arg2, out);
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
			array_map<const Proof*, Proof*>& proof_map,
			hash_map<const Formula*, Formula*>& formula_map)
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

		if (!Proof::clone(src.definitions[0], dst.definitions[0], proof_map, formula_map)) {
			core::free(dst.types); core::free(dst.negated_types);
			core::free(dst.relations); core::free(dst.negated_relations);
			core::free(dst.definitions); core::free(dst.existential_intro_nodes);
			core::free(dst.function_values); return false;
		}
		dst.definitions[0]->reference_count++;
		dst.definitions.length++;

		for (unsigned int i = 0; i < src.types.size; i++) {
			if (!::clone(src.types.keys[i], dst.types.keys[dst.types.size], formula_map)) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.types.values[i], dst.types.values[dst.types.size], proof_map, formula_map)) {
				core::free(dst.types.keys[dst.types.size]);
				core::free(dst); return false;
			}
			dst.types.size++;
		} for (unsigned int i = 0; i < src.negated_types.size; i++) {
			if (!::clone(src.negated_types.keys[i], dst.negated_types.keys[dst.negated_types.size], formula_map)) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.negated_types.values[i], dst.negated_types.values[dst.negated_types.size], proof_map, formula_map)) {
				core::free(dst.negated_types.keys[dst.negated_types.size]);
				core::free(dst); return false;
			}
			dst.negated_types.size++;
		} for (unsigned int i = 0; i < src.relations.size; i++) {
			if (!Proof::clone(src.relations.values[i], dst.relations.values[dst.relations.size], proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
			dst.relations.keys[dst.relations.size++] = src.relations.keys[i];
		} for (unsigned int i = 0; i < src.negated_relations.size; i++) {
			if (!Proof::clone(src.negated_relations.values[i], dst.negated_relations.values[dst.negated_relations.size], proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
			dst.negated_relations.keys[dst.negated_relations.size++] = src.negated_relations.keys[i];
		} for (unsigned int i = 1; i < src.definitions.length; i++) {
			if (!Proof::clone(src.definitions[i], dst.definitions[dst.definitions.length], proof_map, formula_map)) {
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
			if (!Proof::clone(src.function_values.values[i], dst.function_values.values[dst.function_values.size], proof_map, formula_map)) {
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

template<typename ProofCalculus>
inline bool get_proof_map(const concept<ProofCalculus>& c,
		hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	for (const auto& entry : c.types) {
		if (!get_formula_map(entry.key, formula_map)
		 || !get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (const auto& entry : c.negated_types) {
		if (!get_formula_map(entry.key, formula_map)
		 || !get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (const auto& entry : c.relations) {
		if (!get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (const auto& entry : c.negated_relations) {
		if (!get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (Proof* proof : c.definitions) {
		if (!get_proof_map(proof, proof_map, formula_map))
			return false;
	} for (Proof* proof : c.existential_intro_nodes) {
		if (!get_proof_map(proof, proof_map, formula_map))
			return false;
	} for (const auto& entry : c.function_values) {
		if (!get_proof_map(entry.value, proof_map, formula_map))
			return false;
	}
	return true;
}

template<typename ProofCalculus, typename Stream>
bool read(concept<ProofCalculus>& c, Stream& in,
		typename ProofCalculus::Proof** proofs,
		typename ProofCalculus::Language** formulas)
{
	decltype(c.types.size) type_count;
	decltype(c.negated_types.size) negated_type_count;
	decltype(c.relations.size) relation_count;
	decltype(c.negated_relations.size) negated_relation_count;
	decltype(c.definitions.length) definition_count;
	decltype(c.existential_intro_nodes.length) existential_intro_node_count;
	decltype(c.function_values.size) function_value_count;

	if (!read(type_count, in)
	 || !read(negated_type_count, in)
	 || !read(relation_count, in)
	 || !read(negated_relation_count, in)
	 || !read(definition_count, in)
	 || !read(existential_intro_node_count, in)
	 || !read(function_value_count, in))
		return false;

	if (!array_map_init(c.types, ((size_t) 1) << (core::log2(type_count + 1) + 1))) {
		return false;
	} else if (!array_map_init(c.negated_types, ((size_t) 1) << (core::log2(negated_type_count + 1) + 1))) {
		free(c.types); return false;
	} else if (!array_map_init(c.relations, ((size_t) 1) << (core::log2(relation_count + 1) + 1))) {
		free(c.negated_types);
		free(c.types); return false;
	} else if (!array_map_init(c.negated_relations, ((size_t) 1) << (core::log2(negated_relation_count + 1) + 1))) {
		free(c.relations); free(c.negated_types);
		free(c.types); return false;
	} else if (!array_init(c.definitions, ((size_t) 1) << (core::log2(definition_count + 1) + 1))) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		return false;
	} else if (!array_init(c.existential_intro_nodes, ((size_t) 1) << (core::log2(existential_intro_node_count + 1) + 1))) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); return false;
	} else if (!array_map_init(c.function_values, ((size_t) 1) << (core::log2(function_value_count + 1) + 1))) {
		free(c.relations); free(c.negated_types);
		free(c.types); free(c.negated_relations);
		free(c.definitions); free(c.existential_intro_nodes);
		return false;
	}

	unsigned int index;
	for (unsigned int i = 0; i < type_count; i++) {
		if (!read(c.types.keys[i], in, formulas)) {
			free(c); return false;
		}
		c.types.keys[i].reference_count = 1;
		if (!read(index, in)) {
			free(c.types.keys[i]);
			free(c); return false;
		}
		c.types.values[i] = proofs[index];
		c.types.values[i]->reference_count++;
		c.types.size++;
	} for (unsigned int i = 0; i < negated_type_count; i++) {
		if (!read(c.negated_types.keys[i], in, formulas)) {
			free(c); return false;
		}
		c.negated_types.keys[i].reference_count = 1;
		if (!read(index, in)) {
			free(c.negated_types.keys[i]);
			free(c); return false;
		}
		c.negated_types.values[i] = proofs[index];
		c.negated_types.values[i]->reference_count++;
		c.negated_types.size++;
	} for (unsigned int i = 0; i < relation_count; i++) {
		if (!read(c.relations.keys[i], in)
		 || !read(index, in))
		{
			free(c); return false;
		}
		c.relations.values[i] = proofs[index];
		c.relations.values[i]->reference_count++;
		c.relations.size++;
	} for (unsigned int i = 0; i < negated_relation_count; i++) {
		if (!read(c.negated_relations.keys[i], in)
		 || !read(index, in))
		{
			free(c); return false;
		}
		c.negated_relations.values[i] = proofs[index];
		c.negated_relations.values[i]->reference_count++;
		c.negated_relations.size++;
	} for (unsigned int i = 0; i < definition_count; i++) {
		if (!read(index, in)) {
			free(c); return false;
		}
		c.definitions[i] = proofs[index];
		if (i == 0) {
			/* we set the initial reference_count to 2 since `theory.free_proof`
			   will free definitions when their reference_count is 1 */
			c.definitions[i]->reference_count += 2;
		} else {
			c.definitions[i]->reference_count++;
		}
		c.definitions.length++;
	} for (unsigned int i = 0; i < existential_intro_node_count; i++) {
		if (!read(index, in)) {
			free(c); return false;
		}
		c.existential_intro_nodes[i] = proofs[index];
		c.existential_intro_nodes.length++;
	} for (unsigned int i = 0; i < function_value_count; i++) {
		if (!read(c.function_values.keys[i], in)
		 || !read(index, in))
		{
			free(c); return false;
		}
		c.function_values.values[i] = proofs[index];
		c.function_values.values[i]->reference_count++;
		c.function_values.size++;
	}
	return true;
}

template<typename ProofCalculus, typename Stream>
bool write(const concept<ProofCalculus>& c, Stream& out,
		const hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		const hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	if (!write(c.types.size, out)
	 || !write(c.negated_types.size, out)
	 || !write(c.relations.size, out)
	 || !write(c.negated_relations.size, out)
	 || !write(c.definitions.length, out)
	 || !write(c.existential_intro_nodes.length, out)
	 || !write(c.function_values.size, out))
		return false;
	for (const auto& entry : c.types) {
		if (!write(entry.key, out, formula_map)
		 || !write(proof_map.get(entry.value), out))
			return false;
	} for (const auto& entry : c.negated_types) {
		if (!write(entry.key, out, formula_map)
		 || !write(proof_map.get(entry.value), out))
			return false;
	} for (const auto& entry : c.relations) {
		if (!write(entry.key, out)
		 || !write(proof_map.get(entry.value), out))
			return false;
	} for (const auto& entry : c.negated_relations) {
		if (!write(entry.key, out)
		 || !write(proof_map.get(entry.value), out))
			return false;
	} for (const Proof* proof : c.definitions) {
		if (!write(proof_map.get(proof), out))
			return false;
	} for (const Proof* proof : c.existential_intro_nodes) {
		if (!write(proof_map.get(proof), out))
			return false;
	} for (const auto& entry : c.function_values) {
		if (!write(entry.key, out)
		 || !write(proof_map.get(entry.value), out))
			return false;
	}
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
	if (first == nullptr) return nullptr;

	second = normalize_quantifiers_with_equality(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = normalize_optimizations(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = simplify_subsets(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = normalize_comparatives(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = simplify_properties(first);
	free(*first); if (first->reference_count == 0) free(first);
	if (second == nullptr) return nullptr;

	first = normalize_functions(second);
	free(*second); if (second->reference_count == 0) free(second);
	if (first == nullptr) return nullptr;

	second = normalize_multiple_universal_quantifiers(first);
	free(*first); if (first->reference_count == 0) free(first);
	return second;
}

/* this is useful in `theory.add_definition` where if any sets are created in
   `set_reasoning.get_subset_axiom`, we need the sets to have a specific size */
struct required_set_size {
	unsigned int min_set_size;
	unsigned int max_set_size;

	required_set_size(unsigned int min_set_size, unsigned int max_set_size) :
			min_set_size(min_set_size), max_set_size(max_set_size) { }
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
	unsigned int new_min_set_size = max(min_set_size, required.min_set_size);
	unsigned int new_max_set_size = min(max_set_size, required.max_set_size);
	if (new_min_set_size > new_max_set_size)
		return false;
	return compute_new_set_size(set_id, sets, out, new_min_set_size, new_max_set_size, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int min_set_size,
		unsigned int max_set_size,
		const required_set_size& required,
		Args&&... visitor)
{
	unsigned int new_min_set_size = max(min_set_size, required.min_set_size);
	unsigned int new_max_set_size = min(max_set_size, required.max_set_size);
	on_free_set(set_id, sets, new_min_set_size, new_max_set_size, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline bool on_new_size_axiom(
		Proof* new_size_axiom,
		const required_set_size& required,
		Args&&... visitor)
{
	return on_new_size_axiom(new_size_axiom, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const required_set_size& required,
		Args&&... visitor)
{
	on_old_size_axiom(old_size_axiom, std::forward<Args>(visitor)...);
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
		unsigned int min_set_size,
		unsigned int max_set_size,
		const array<typename ProofCalculus::Proof*>& freeable_axioms,
		Args&&... visitor)
{
	on_free_set(set_id, sets, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline bool on_new_size_axiom(
		Proof* new_size_axiom,
		const array<Proof*>& freeable_axioms,
		Args&&... visitor)
{
	return on_new_size_axiom(new_size_axiom, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const array<Proof*>& freeable_axioms,
		Args&&... visitor)
{
	on_old_size_axiom(old_size_axiom, std::forward<Args>(visitor)...);
}

template<typename Formula>
struct set_changes {
	array<Formula*> old_set_axioms;
	array<Formula*> new_set_axioms;

	set_changes() : old_set_axioms(4), new_set_axioms(4) { }
	~set_changes() { free_helper(); }

	inline bool new_set(Formula* axiom) {
		if (!new_set_axioms.add(axiom))
			return false;
		axiom->reference_count++;
		return true;
	}

	inline bool old_set(Formula* axiom) {
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
	return compute_new_set_size(set_id, sets, out, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... Args>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int min_set_size,
		unsigned int max_set_size,
		set_changes<typename ProofCalculus::Language>& set_diff,
		Args&&... visitor)
{
	on_free_set(set_id, sets, min_set_size, max_set_size, std::forward<Args>(visitor)...);
}

template<typename Proof, typename... Args>
inline bool on_new_size_axiom(
		Proof* new_size_axiom,
		set_changes<typename Proof::FormulaType>& set_diff,
		Args&&... visitor)
{
	return on_new_size_axiom(new_size_axiom, std::forward<Args>(visitor)...)
		&& set_diff.new_set(new_size_axiom->formula);
}

template<typename Proof, typename... Args>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		set_changes<typename Proof::FormulaType>& set_diff,
		Args&&... visitor)
{
	on_old_size_axiom(old_size_axiom, std::forward<Args>(visitor)...);
	set_diff.old_set(old_size_axiom->formula);
}

template<typename Formula, typename... Args>
inline void on_subtract_changes(const set_changes<Formula>& set_diff, Args&&... visitor) {
	on_subtract_changes(std::forward<Args>(visitor)...);
}


/* forward declarations */

template<typename ProofCalculus, typename Canonicalizer> struct theory;

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer>
bool is_impossible(typename ProofCalculus::Language*, const theory<ProofCalculus, Canonicalizer>&);


template<typename ProofCalculus, typename Canonicalizer>
struct theory
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef Canonicalizer FormulaCanonicalizer;

	template<unsigned int Value> using Constants = typename Term::template constants<Value>;
	template<unsigned int Value> using Variables = typename Term::template variables<Value>;

	unsigned int new_constant_offset;

	/* A map from `t` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]t` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]t`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]t` or `~[u/0]t`
	   are in the theory. */
	hash_map<Term, pair<array<instance>, array<instance>>> atoms;

	/* A map from `R` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]R` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]R`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]R` or `~[u/0]R`
	   are in the theory. */
	hash_map<relation, pair<array<instance>, array<instance>>> relations;

	concept<ProofCalculus>* ground_concepts;
	unsigned int ground_concept_capacity;
	unsigned int ground_axiom_count;
	context ctx;

	hash_map<Term, unsigned int> reverse_definitions;

	array_map<Term, Proof*> constant_types;
	array_map<Term, Proof*> constant_negated_types;

	array<Proof*> observations;
	set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer> sets;

	struct proof_node {
		Formula* formula;
		Proof* proof;
		array_map<unsigned int, Proof*> set_definitions;

		proof_node(const proof_node& src) : formula(src.formula), proof(src.proof), set_definitions(src.set_definitions.capacity) {
			for (const auto& entry : src.set_definitions) {
				set_definitions.keys[set_definitions.size] = entry.key;
				set_definitions.values[set_definitions.size++] = entry.value;
			}
		}

		static inline void move(const proof_node& src, proof_node& dst) {
			dst.formula = src.formula;
			dst.proof = src.proof;
			core::move(src.set_definitions, dst.set_definitions);
		}

		static inline void free(proof_node& node) {
			core::free(*node.formula);
			if (node.formula->reference_count == 0)
				core::free(node.formula);
			core::free(node.set_definitions);
		}
	};

	array<proof_node> disjunction_intro_nodes;
	array<proof_node> negated_conjunction_nodes;
	array<proof_node> implication_intro_nodes;
	array<proof_node> existential_intro_nodes;

	array<pair<Proof*, bool>> implication_axioms;

	Proof* empty_set_axiom;
	array<Proof*> built_in_axioms;
	array<unsigned int> built_in_sets;

	Term* NAME_ATOM;

	theory(const array<Formula*>& seed_axioms, unsigned int new_constant_offset) :
			new_constant_offset(new_constant_offset), atoms(64), relations(64),
			ground_concept_capacity(64), ground_axiom_count(0),
			reverse_definitions(256), constant_types(8),
			constant_negated_types(8), observations(32),
			disjunction_intro_nodes(16), negated_conjunction_nodes(16),
			implication_intro_nodes(16), existential_intro_nodes(16),
			implication_axioms(16), built_in_axioms(8), built_in_sets(8)
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

		/* add the seed axioms */
		if (!built_in_axioms.ensure_capacity(seed_axioms.length))
			exit(EXIT_FAILURE);
		for (Formula* seed_axiom : seed_axioms) {
			uint_fast8_t arity = 0;
			Formula* operand = seed_axiom;
			while (operand->type == FormulaType::FOR_ALL) {
				operand = operand->quantifier.operand;
				arity++;
			}
			if (arity == 0) {
				fprintf(stderr, "ERROR: Non-universally-quantified seed axioms are not supported.\n");
				exit(EXIT_FAILURE);
			} else {
				Formula* left; Formula* right;
				if (operand->type != FormulaType::IF_THEN) {
					left = Formula::new_true();
					left->reference_count--;
					right = operand->binary.right;
				} else {
					left = operand->binary.left;
					right = operand->binary.right;
				}

				bool is_antecedent_new = false;
				bool is_consequent_new = false;
				unsigned int antecedent_set, consequent_set;
				built_in_axioms[built_in_axioms.length] = get_subset_axiom<true>(left, right, arity, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new);
				if (built_in_axioms[built_in_axioms.length] == nullptr || !built_in_sets.add(antecedent_set) || !built_in_sets.add(consequent_set)) {
					print("Failed to add built-in axiom '", stderr); print(*seed_axiom, stderr, *debug_terminal_printer); print("'.\n", stderr);
					exit(EXIT_FAILURE);
				}
				built_in_axioms[built_in_axioms.length++]->reference_count++;
			}
		}
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
		} for (auto entry : reverse_definitions) {
			core::free(entry.key);
		} for (Proof* proof : observations) {
			core::free(*proof);
			if (proof->reference_count == 0)
				core::free(proof);
		} for (proof_node& node : disjunction_intro_nodes) {
			core::free(node);
		} for (proof_node& node : negated_conjunction_nodes) {
			core::free(node);
		} for (proof_node& node : implication_intro_nodes) {
			core::free(node);
		} for (proof_node& node : existential_intro_nodes) {
			core::free(node);
		} for (auto entry : implication_axioms) {
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

		for (auto entry : constant_types) {
			core::free(entry.key);
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		} for (auto entry : constant_negated_types) {
			core::free(entry.key);
			core::free(*entry.value);
			if (entry.value->reference_count == 0)
				core::free(entry.value);
		}

		core::free(*NAME_ATOM);
		if (NAME_ATOM->reference_count == 0)
			core::free(NAME_ATOM);
		core::free(*empty_set_axiom);
		if (empty_set_axiom->reference_count == 0)
			core::free(empty_set_axiom);
		for (Proof* axiom : built_in_axioms) {
			core::free(*axiom);
			sets.free_subset_axiom(axiom);
		}
	}

	static inline void free(theory<ProofCalculus, Canonicalizer>& T) {
		T.free_helper();
		core::free(T.atoms);
		core::free(T.relations);
		core::free(T.reverse_definitions);
		core::free(T.constant_types);
		core::free(T.constant_negated_types);
		core::free(T.observations);
		core::free(T.ctx);
		core::free(T.sets);
		core::free(T.disjunction_intro_nodes);
		core::free(T.negated_conjunction_nodes);
		core::free(T.implication_intro_nodes);
		core::free(T.existential_intro_nodes);
		core::free(T.implication_axioms);
		core::free(T.built_in_axioms);
		core::free(T.built_in_sets);
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

		/* remove all elements in sets that contain this concept */
		for (unsigned int i = 2; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr)
				continue;
			for (unsigned int j = sets.sets[i].element_count(); j > 0; j--) {
				const tuple_element* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * (j - 1));
				bool contains_constant = false;
				for (uint_fast8_t k = 0; k < sets.sets[i].arity && !contains_constant; k++)
					if (element_src[k].type == tuple_element_type::CONSTANT && element_src[k].constant == id) contains_constant = true;
				if (contains_constant)
					sets.remove_element_at(i, j - 1);
			}
		}
	}

	void try_free_concept_id(unsigned int id) {
		if (id < new_constant_offset) return;
		const concept<ProofCalculus>& c = ground_concepts[id - new_constant_offset];
		if (c.types.keys == nullptr) return;
		if (c.types.size != 0 || c.negated_types.size != 0 || c.relations.size != 0
		 || c.negated_relations.size != 0 || c.definitions.length > 1 || c.definitions[0]->reference_count > 2
		 || c.existential_intro_nodes.length != 0 || c.function_values.size != 0)
			return;

		bool contains;
		unsigned int count = sets.symbols_in_formulas.counts.get(id, contains);
		if (contains && count > 0) return;
		free_concept_id(id);
	}

	template<bool PrintProvableElements = false, typename Stream, typename... Printer>
	bool print_axioms(Stream&& out, Printer&&... printer) const {
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys == NULL) continue;
			if (!print("Object ", out) || !print(i + new_constant_offset, out, std::forward<Printer>(printer)...) || !print(":\n", out)
			 || !ground_concepts[i].print_axioms(out, "  ", std::forward<Printer>(printer)...)) return false;
		}
		if (constant_types.size > 0 || constant_negated_types.size > 0) {
			if (!print("Other axioms:\n", out)) return false;
			for (auto entry : constant_types) {
				if (!print("  ", out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
			} for (auto entry : constant_negated_types) {
				if (!print("  ", out) || !print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out)) return false;
			}
		}

		if (implication_axioms.length > 0) {
			print("Implication axioms:\n", out);
			for (const auto entry : implication_axioms) {
				if (!print("  ", out) || !print(*entry.key->formula, out, std::forward<Printer>(printer)...) || !print(", Antecedent satisfied: ", out)) return false;
				if (entry.value && !print("true.\n", out)) return false;
				else if (!entry.value && !print("false.\n", out)) return false;
			}
		}
		return sets.template print_axioms<PrintProvableElements>(out, std::forward<Printer>(printer)...);
	}

	template<typename Stream, typename... Printer>
	bool print_disjunction_introduction(const Proof* proof, Stream&& out, Printer&&... printer) const {
		if (proof->type == ProofType::EXISTENTIAL_INTRODUCTION) {
			/* find the formula corresponding to this node */
			unsigned int index;
			for (index = 0; index < existential_intro_nodes.length; index++)
				if (existential_intro_nodes[index].proof == proof) break;
			if (index == existential_intro_nodes.length) {
				fprintf(stderr, "print_disjunction_introduction WARNING: Found existential introduction in proof that is not an element of `theory.existential_intro_nodes`.\n");
				if (!print("  Existential introduction instantiated with ", out)
				 || !print(*proof->operands[2]->term, out, std::forward<Printer>(printer)...)
				 || !print(".\n", out))
					return false;
			} else {
				if (!print("  Existential introduction at index ", out)
				 || !print(index, out) || !print(": ", out)
				 || !print(*existential_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
				 || !print(" instantiated with ", out)
				 || !print(*proof->operands[2]->term, out, std::forward<Printer>(printer)...)
				 || !print(".\n", out))
					return false;
			}
		} else if (proof->type == ProofType::IMPLICATION_INTRODUCTION) {
			unsigned int index;
			for (index = 0; index < implication_intro_nodes.length; index++)
				if (implication_intro_nodes[index].proof == proof) break;
			if (index == implication_intro_nodes.length) {
				fprintf(stderr, "print_disjunction_introduction WARNING: Found implication introduction in proof that is not an element of `theory.implication_intro_nodes`.\n");
				if (!print("  Implication introduction", out)
				 || !print((proof->operands[0]->type == ProofType::FALSITY_ELIMINATION) ?
					" proved via the negation of the antecedent.\n" :
					" proved via the consequent.\n", out))
					return false;
			} else {
				if (!print("  Implication introduction at index ", out)
				 || !print(index, out) || !print(": ", out)
				 || !print(*implication_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
				 || !print((proof->operands[0]->type == ProofType::FALSITY_ELIMINATION) ?
					" proved via the negation of the antecedent.\n" :
					" proved via the consequent.\n", out))
					return false;
			}
		} else if (proof->type == ProofType::DISJUNCTION_INTRODUCTION) {
			unsigned int index;
			for (index = 0; index < disjunction_intro_nodes.length; index++)
				if (disjunction_intro_nodes[index].proof == proof) break;
			if (index == disjunction_intro_nodes.length) {
				fprintf(stderr, "print_disjunction_introduction WARNING: Found disjunction introduction in proof that is not an element of `theory.disjunction_intro_nodes`.\n");
				if (!print("  Disjunction introduction instantiated with disjunct at index ", out)
				 || !print(proof->operands[2]->parameter, out) || print(".\n", out))
					return false;
			} else {
				if (!print("  Disjunction introduction at index ", out)
				 || !print(index, out) || !print(": ", out)
				 || !print(*disjunction_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
				 || !print(" instantiated with disjunct at index ", out)
				 || !print(proof->operands[2]->parameter, out) || !print(".\n", out))
					return false;
			}
		} else if (proof->type == ProofType::PROOF_BY_CONTRADICTION
				&& proof->operands[0]->type == ProofType::NEGATION_ELIMINATION
				&& proof->operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
				&& proof->operands[0]->operands[0]->operands[0] == proof->operands[1])
		{
			unsigned int index;
			for (index = 0; index < negated_conjunction_nodes.length; index++)
				if (negated_conjunction_nodes[index].proof == proof) break;
			if (index == negated_conjunction_nodes.length) {
				fprintf(stderr, "print_disjunction_introduction WARNING: Found proof of negated conjunction in proof that is not an element of `theory.negated_conjunction_nodes`.\n");
				if (!print("  Negated conjunction instantiated with conjunct at indices ", out)
				 || !print(proof->operands[0]->operands[0]->operands[1]->parameters[0], out) || print(".\n", out))
					return false;
			} else {
				if (!print("  Negated conjunction at index ", out)
				 || !print(index, out) || !print(": ", out)
				 || !print(*negated_conjunction_nodes[index].formula, out, std::forward<Printer>(printer)...)
				 || !print(" instantiated with conjunct at indices ", out)
				 || !print(proof->operands[0]->operands[0]->operands[1]->parameters[0], out) || !print(".\n", out))
					return false;
			}
		}

		unsigned int operand_count;
		const Proof* const* operands;
		proof->get_subproofs(operands, operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			if (operands[i] == NULL) continue;
			if (!print_disjunction_introduction(operands[i], out, std::forward<Printer>(printer)...))
				return false;
		}
		return true;
	}

	template<typename Stream, typename... Printer>
	bool print_disjunction_introductions(Stream&& out, Printer&&... printer) const {
		for (const Proof* observation : observations) {
			if (observation->type == ProofType::EXISTENTIAL_INTRODUCTION) {
				unsigned int index;
				for (index = 0; index < existential_intro_nodes.length; index++)
					if (existential_intro_nodes[index].proof == observation) break;
				if (index == existential_intro_nodes.length) {
					fprintf(stderr, "print_disjunction_introductions WARNING: Found existential introduction in proof that is not an element of `theory.existential_intro_nodes`.\n");
					if (fprintf(out, "Proof at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				} else {
					if (!print("Proof of ", out)
					 || !print(*existential_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
					 || fprintf(out, " at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				}
			} else if (observation->type == ProofType::IMPLICATION_INTRODUCTION) {
				unsigned int index;
				for (index = 0; index < implication_intro_nodes.length; index++)
					if (implication_intro_nodes[index].proof == observation) break;
				if (index == implication_intro_nodes.length) {
					fprintf(stderr, "print_disjunction_introductions WARNING: Found implication introduction in proof that is not an element of `theory.implication_intro_nodes`.\n");
					if (fprintf(out, "Proof at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				} else {
					if (!print("Proof of ", out)
					 || !print(*implication_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
					 || fprintf(out, " at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				}
			} else if (observation->type == ProofType::DISJUNCTION_INTRODUCTION) {
				unsigned int index;
				for (index = 0; index < disjunction_intro_nodes.length; index++)
					if (disjunction_intro_nodes[index].proof == observation) break;
				if (index == disjunction_intro_nodes.length) {
					fprintf(stderr, "print_disjunction_introductions WARNING: Found disjunction introduction in proof that is not an element of `theory.disjunction_intro_nodes`.\n");
					if (fprintf(out, "Proof at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				} else {
					if (!print("Proof of ", out)
					 || !print(*disjunction_intro_nodes[index].formula, out, std::forward<Printer>(printer)...)
					 || fprintf(out, " at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				}
			} else if (observation->type == ProofType::PROOF_BY_CONTRADICTION
					&& observation->operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& observation->operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
					&& observation->operands[0]->operands[0]->operands[0] == observation->operands[1])
			{
				unsigned int index;
				for (index = 0; index < negated_conjunction_nodes.length; index++)
					if (negated_conjunction_nodes[index].proof == observation) break;
				if (index == negated_conjunction_nodes.length) {
					fprintf(stderr, "print_disjunction_introductions WARNING: Found proof of negated conjunction in proof that is not an element of `theory.negated_conjunction_nodes`.\n");
					if (fprintf(out, "Proof at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				} else {
					if (!print("Proof of ", out)
					 || !print(*negated_conjunction_nodes[index].formula, out, std::forward<Printer>(printer)...)
					 || fprintf(out, " at address 0x%lx:\n", (size_t) observation) < 0)
						return false;
				}
			}

			if (!print_disjunction_introduction(observation, out, std::forward<Printer>(printer)...))
				return false;
		}
		return true;
	}

	template<typename Stream, typename... Printer>
	bool print_proofs(Stream&& out, Printer&&... printer) const {
		for (const Proof* observation : observations) {
			if (!print<built_in_predicates, typename ProofCalculus::ProofCanonicalizer, ProofCalculus::Intuitionistic>(*observation, out, std::forward<Printer>(printer)...) || !print('\n', out))
				return false;
		}
		return true;
	}

	static inline Term* get_name_from_scope(Formula* scope, unsigned int variable)
	{
		if (scope->type != FormulaType::EXISTS
		 || scope->quantifier.operand->type != FormulaType::AND)
			return nullptr;
		Formula* operand = scope->quantifier.operand;
		bool is_name = false;
		bool has_arg1 = false;
		Term* arg2 = nullptr;
		for (unsigned int i = 0; i < operand->array.length; i++)
		{
			Formula* conjunct = operand->array.operands[i];
			if (!is_name && conjunct->type == FormulaType::UNARY_APPLICATION
			 && conjunct->binary.left->type == TermType::CONSTANT
			 && conjunct->binary.left->constant == (unsigned int) built_in_predicates::NAME
			 && conjunct->binary.right->type == TermType::VARIABLE
			 && conjunct->binary.right->variable == scope->quantifier.variable)
			{
				is_name = true;
			} else if (!has_arg1 && conjunct->type == FormulaType::EQUALS
					&& conjunct->binary.left->type == TermType::UNARY_APPLICATION
					&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
					&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
					&& conjunct->binary.left->binary.right->variable == scope->quantifier.variable
					&& conjunct->binary.right->type == TermType::VARIABLE
					&& conjunct->binary.right->variable == variable)
			{
				has_arg1 = true;
			} else if (arg2 == nullptr && conjunct->type == FormulaType::EQUALS
					&& conjunct->binary.left->type == TermType::UNARY_APPLICATION
					&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
					&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
					&& conjunct->binary.left->binary.right->variable == scope->quantifier.variable)
			{
				arg2 = conjunct->binary.right;
			}
		}

		if (is_name && has_arg1)
			return arg2;
		else return nullptr;
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

		/* check if a name is defined by a set that contains this object */
		hash_set<unsigned int> ancestors(16);
		for (const auto& entry : ground_concepts[constant - new_constant_offset].types) {
			bool contains;
			unsigned int set_id = sets.set_ids.get(entry.key, contains);
			if (!contains || sets.sets[set_id].arity != 1)
				continue;

			if (!sets.get_ancestors(set_id, ancestors))
				return false;
		}
		for (unsigned int ancestor : ancestors) {
			Formula* formula = sets.sets[ancestor].set_formula();
			if (formula->type == FormulaType::AND) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					Term* arg2 = get_name_from_scope(formula->array.operands[i], 1);
					if (arg2 == nullptr || arg2->type != TermType::STRING || name_terms.contains(arg2))
						continue;
					if (!name_terms.add(arg2))
						return false;
				}
			} else {
				Term* arg2 = get_name_from_scope(formula, 1);
				if (arg2 == nullptr || arg2->type != TermType::STRING || name_terms.contains(arg2))
					continue;
				if (!name_terms.add(arg2))
					return false;
			}
		}
		return true;
	}

	inline bool is_concept_name(unsigned int constant, const string& name) const
	{
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
			if (function_value_axiom->formula->binary.right->str == name)
				return true;
		}
		return false;
	}

	inline bool get_extra_axioms(array<Formula*>& extra_axioms) const
	{
		for (unsigned int i = 2; i < sets.set_count + 1; i++) {
			if (built_in_sets.contains(i) || sets.sets[i].size_axioms.data == nullptr)
				continue;
			for (Proof* size_axiom : sets.sets[i].size_axioms) {
				if (!extra_axioms.contains(size_axiom->formula)
				 && !extra_axioms.add(size_axiom->formula))
					return false;
			}
		}
		return true;
	}

	static inline bool clone(
			const theory<ProofCalculus, Canonicalizer>& src,
			theory<ProofCalculus, Canonicalizer>& dst,
			array_map<const Proof*, Proof*>& proof_map,
			hash_map<const Formula*, Formula*>& formula_map)
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
		dst.ground_axiom_count = src.ground_axiom_count;
		if (!array_init(dst.observations, src.observations.capacity)) {
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		}
		if (!Proof::clone(src.empty_set_axiom, dst.empty_set_axiom, proof_map, formula_map)) {
			core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		} if (!array_init(dst.built_in_axioms, src.built_in_axioms.capacity)) {
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		}
		for (unsigned int i = 0; i < src.built_in_axioms.length; i++) {
			if (!Proof::clone(src.built_in_axioms[i], dst.built_in_axioms[i], proof_map, formula_map)) {
				for (unsigned int j = 0; j < i; j++) {
					core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
				}
				core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
				core::free(dst.observations);
				core::free(dst.ground_concepts);
				core::free(dst.atoms); core::free(dst.relations);
				return false;
			}
			dst.built_in_axioms.length++;
		}
		if (!::clone(src.NAME_ATOM, dst.NAME_ATOM, formula_map)) {
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			return false;
		}
		if (!set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>::clone(src.sets, dst.sets, proof_map, formula_map)) {
			core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.disjunction_intro_nodes, src.disjunction_intro_nodes.capacity)) {
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.negated_conjunction_nodes, src.negated_conjunction_nodes.capacity)) {
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.implication_intro_nodes, src.implication_intro_nodes.capacity)) {
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.existential_intro_nodes, src.existential_intro_nodes.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.implication_axioms, src.implication_axioms.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_init(dst.built_in_sets, src.built_in_sets.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.existential_intro_nodes);
			core::free(dst.implication_axioms);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_map_init(dst.constant_types, src.constant_types.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.existential_intro_nodes);
			core::free(dst.implication_axioms);
			core::free(dst.built_in_sets);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!array_map_init(dst.constant_negated_types, src.constant_negated_types.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.existential_intro_nodes);
			core::free(dst.implication_axioms);
			core::free(dst.built_in_sets);
			core::free(dst.constant_types);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!hash_map_init(dst.reverse_definitions, src.reverse_definitions.table.capacity)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.existential_intro_nodes);
			core::free(dst.implication_axioms);
			core::free(dst.built_in_sets);
			core::free(dst.constant_types);
			core::free(dst.constant_negated_types);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		} else if (!context::clone(src.ctx, dst.ctx, formula_map)) {
			core::free(dst.implication_intro_nodes);
			core::free(dst.negated_conjunction_nodes);
			core::free(dst.disjunction_intro_nodes);
			core::free(dst.existential_intro_nodes);
			core::free(dst.implication_axioms);
			core::free(dst.built_in_sets);
			core::free(dst.constant_types);
			core::free(dst.constant_negated_types);
			core::free(dst.reverse_definitions);
			core::free(dst.sets); core::free(dst.observations);
			core::free(dst.ground_concepts);
			core::free(dst.atoms); core::free(dst.relations);
			for (unsigned int j = 0; j < dst.built_in_axioms.length; j++) {
				core::free(*dst.built_in_axioms[j]); if (dst.built_in_axioms[j]->reference_count == 0) core::free(dst.built_in_axioms[j]);
			} core::free(dst.built_in_axioms);
			core::free(*dst.empty_set_axiom); if (dst.empty_set_axiom->reference_count == 0) core::free(dst.empty_set_axiom);
			core::free(*dst.NAME_ATOM); if (dst.NAME_ATOM->reference_count == 0) core::free(dst.NAME_ATOM);
			return false;
		}


		for (const auto& entry : src.atoms) {
			unsigned int bucket = dst.atoms.table.index_to_insert(entry.key);
			if (!array_init(dst.atoms.values[bucket].key, entry.value.key.capacity)) {
				core::free(dst); return false;
			} else if (!array_init(dst.atoms.values[bucket].value, entry.value.value.capacity)) {
				core::free(dst.atoms.values[bucket].key);
				core::free(dst); return false;
			} else if (!::clone(entry.key, dst.atoms.table.keys[bucket])) {
				core::free(dst.atoms.values[bucket].key);
				core::free(dst.atoms.values[bucket].value);
				core::free(dst); return false;
			}
			dst.atoms.table.size++;
			for (const instance& id : entry.value.key)
				dst.atoms.values[bucket].key[dst.atoms.values[bucket].key.length++] = id;
			for (const instance& id : entry.value.value)
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
			for (const instance& id : entry.value.key)
				dst.relations.values[bucket].key[dst.relations.values[bucket].key.length++] = id;
			for (const instance& id : entry.value.value)
				dst.relations.values[bucket].value[dst.relations.values[bucket].value.length++] = id;
		} for (Proof* observation : src.observations) {
			Proof* new_observation;
			if (!Proof::clone(observation, new_observation, proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
			dst.observations.add(new_observation);
		} for (unsigned int i = 0; i < src.ground_concept_capacity; i++) {
			if (src.ground_concepts[i].types.keys == nullptr) continue;
			if (!concept<ProofCalculus>::clone(src.ground_concepts[i], dst.ground_concepts[i], proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
		} for (unsigned int i = 0; i < src.constant_types.size; i++) {
			if (!::clone(src.constant_types.keys[i], dst.constant_types.keys[dst.constant_types.size], formula_map)) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.constant_types.values[i], dst.constant_types.values[dst.constant_types.size], proof_map, formula_map)) {
				core::free(dst.constant_types.keys[dst.constant_types.size]);
				core::free(dst); return false;
			}
			dst.constant_types.size++;
		} for (unsigned int i = 0; i < src.constant_negated_types.size; i++) {
			if (!::clone(src.constant_negated_types.keys[i], dst.constant_negated_types.keys[dst.constant_negated_types.size], formula_map)) {
				core::free(dst); return false;
			} else if (!Proof::clone(src.constant_negated_types.values[i], dst.constant_negated_types.values[dst.constant_negated_types.size], proof_map, formula_map)) {
				core::free(dst.constant_negated_types.keys[dst.constant_negated_types.size]);
				core::free(dst); return false;
			}
			dst.constant_negated_types.size++;
		} for (const proof_node& src_node : src.disjunction_intro_nodes) {
			unsigned int index = proof_map.index_of(src_node.proof);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			proof_node& dst_node = dst.disjunction_intro_nodes[dst.disjunction_intro_nodes.length];
			dst_node.proof = proof_map.values[index];
			if (!::clone(src_node.formula, dst_node.formula, formula_map)) {
				core::free(dst);
				return false;
			} else if (!array_map_init(dst_node.set_definitions, src_node.set_definitions.capacity)) {
				core::free(*src_node.formula); if (src_node.formula->reference_count == 0) core::free(src_node.formula);
				core::free(dst); return false;
			}
			for (const auto& entry : src_node.set_definitions) {
				index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
				if (index == proof_map.size)
					fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
				dst_node.set_definitions.keys[dst_node.set_definitions.size] = entry.key;
				dst_node.set_definitions.values[dst_node.set_definitions.size++] = proof_map.values[index];
			}
			dst.disjunction_intro_nodes.length++;
		} for (const proof_node& src_node : src.negated_conjunction_nodes) {
			unsigned int index = proof_map.index_of(src_node.proof);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			proof_node& dst_node = dst.negated_conjunction_nodes[dst.negated_conjunction_nodes.length];
			dst_node.proof = proof_map.values[index];
			if (!::clone(src_node.formula, dst_node.formula, formula_map)) {
				core::free(dst);
				return false;
			} else if (!array_map_init(dst_node.set_definitions, src_node.set_definitions.capacity)) {
				core::free(*src_node.formula); if (src_node.formula->reference_count == 0) core::free(src_node.formula);
				core::free(dst); return false;
			}
			for (const auto& entry : src_node.set_definitions) {
				index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
				if (index == proof_map.size)
					fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
				dst_node.set_definitions.keys[dst_node.set_definitions.size] = entry.key;
				dst_node.set_definitions.values[dst_node.set_definitions.size++] = proof_map.values[index];
			}
			dst.negated_conjunction_nodes.length++;
		} for (const proof_node& src_node : src.implication_intro_nodes) {
			unsigned int index = proof_map.index_of(src_node.proof);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			proof_node& dst_node = dst.implication_intro_nodes[dst.implication_intro_nodes.length];
			dst_node.proof = proof_map.values[index];
			if (!::clone(src_node.formula, dst_node.formula, formula_map)) {
				core::free(dst);
				return false;
			} else if (!array_map_init(dst_node.set_definitions, src_node.set_definitions.capacity)) {
				core::free(*src_node.formula); if (src_node.formula->reference_count == 0) core::free(src_node.formula);
				core::free(dst); return false;
			}
			for (const auto& entry : src_node.set_definitions) {
				index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
				if (index == proof_map.size)
					fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
				dst_node.set_definitions.keys[dst_node.set_definitions.size] = entry.key;
				dst_node.set_definitions.values[dst_node.set_definitions.size++] = proof_map.values[index];
			}
			dst.implication_intro_nodes.length++;
		} for (const proof_node& src_node : src.existential_intro_nodes) {
			unsigned int index = proof_map.index_of(src_node.proof);
#if !defined(NDEBUG)
			if (index == proof_map.size)
				fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
			proof_node& dst_node = dst.existential_intro_nodes[dst.existential_intro_nodes.length];
			dst_node.proof = proof_map.values[index];
			if (!::clone(src_node.formula, dst_node.formula, formula_map)) {
				core::free(dst);
				return false;
			} else if (!array_map_init(dst_node.set_definitions, src_node.set_definitions.capacity)) {
				core::free(*src_node.formula); if (src_node.formula->reference_count == 0) core::free(src_node.formula);
				core::free(dst); return false;
			}
			for (const auto& entry : src_node.set_definitions) {
				index = proof_map.index_of(entry.value);
#if !defined(NDEBUG)
				if (index == proof_map.size)
					fprintf(stderr, "theory.clone WARNING: Given proof does not exist in `proof_map`.\n");
#endif
				dst_node.set_definitions.keys[dst_node.set_definitions.size] = entry.key;
				dst_node.set_definitions.values[dst_node.set_definitions.size++] = proof_map.values[index];
			}
			dst.existential_intro_nodes.length++;
		} for (const auto& entry : src.implication_axioms) {
			if (!Proof::clone(entry.key, dst.implication_axioms[dst.implication_axioms.length].key, proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
			dst.implication_axioms[dst.implication_axioms.length].value = entry.value;
			dst.implication_axioms.length++;
		}

		for (unsigned int concept = 0; concept < dst.ground_concept_capacity; concept++) {
			if (dst.ground_concepts[concept].types.keys == nullptr) continue;
			auto& c = dst.ground_concepts[concept];
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

				/* add this definition to the reverse definition map */
				if (i > 0) {
					unsigned int index = dst.reverse_definitions.table.index_to_insert(*definition->formula->binary.right);
					dst.reverse_definitions.table.keys[index] = *definition->formula->binary.right;
					dst.reverse_definitions.values[index] = dst.new_constant_offset + concept;
					dst.reverse_definitions.table.size++;
				}
			}
		}
		for (unsigned int built_in_set : src.built_in_sets)
			dst.built_in_sets[dst.built_in_sets.length++] = built_in_set;
		return true;
	}

	static inline bool clone(
			const theory<ProofCalculus, Canonicalizer>& src,
			theory<ProofCalculus, Canonicalizer>& dst,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		array_map<const Proof*, Proof*> proof_map(64);
		return clone(src, dst, proof_map, formula_map);
	}

	static inline bool is_empty(const theory<ProofCalculus, Canonicalizer>& T) {
		return (T.ground_concepts == nullptr);
	}

	static inline void set_empty(theory<ProofCalculus, Canonicalizer>& T) {
		T.ground_concepts = nullptr;
	}

	template<typename... Args>
	Proof* add_formula(Formula* formula, set_changes<Formula>& set_diff, unsigned int& new_constant, Args&&... args)
	{
		if (!observations.ensure_capacity(observations.length + 1)
		 || !ctx.referent_iterators.ensure_capacity(ctx.referent_iterators.length + 1)
		 || !ctx.bindings.ensure_capacity(ctx.bindings.length + 1)
		 || !::init(ctx.referent_iterators[ctx.referent_iterators.length], formula))
			return nullptr;
		ctx.referent_iterators.length++;
		referent_iterator& ref_iterator = ctx.referent_iterators.last();

		array_map<unsigned int, const hol_term*> bound_variables(8);
		if (!get_referents(formula, ref_iterator, bound_variables)
		 || !::init(ctx.bindings[ctx.bindings.length], ref_iterator.anaphora.size))
		{
			core::free(ref_iterator);
			ctx.referent_iterators.length--;
			return nullptr;
		}
		ctx.bindings.length++;
		anaphora_binding& binding = ctx.bindings.last();

		std::set<referent_iterator_state> queue;
		referent_iterator_state initial_state(&ref_iterator, 0.0);
		for (unsigned int j = 0; j < ref_iterator.anaphora.size; j++)
			initial_state.indices[j] = 0;
		queue.insert(initial_state);

		hol_term* resolved_formulas[2];
		double resolved_log_probabilities[2];
		unsigned int resolved_formula_count = 0;
		Proof* new_proof = nullptr;
/* TODO: for debugging; delete this */
//print("add_formula: Adding formula ", stderr); print(*formula, stderr, *debug_terminal_printer); print('\n', stderr);
		while (!queue.empty() && new_proof == nullptr)
		{
			auto last = queue.cend(); last--;
			referent_iterator_state state = *last;
			queue.erase(last);

			process_referent_iterator(state, queue, resolved_formulas, resolved_log_probabilities, resolved_formula_count);

			for (unsigned int i = 0; i < resolved_formula_count; i++) {
/* TODO: for debugging; delete this */
//print("add_formula: Adding resolved formula ", stderr); print(*resolved_formulas[i], stderr, *debug_terminal_printer); print('\n', stderr);
				new_proof = add_formula_helper(resolved_formulas[i], set_diff, new_constant, std::forward<Args>(args)...);
				if (new_proof != nullptr) {
					observations[observations.length++] = new_proof;

					/* record this anaphora binding in `ctx` */
					for (unsigned int j = 0; j < ref_iterator.anaphora.size; j++)
						binding.indices[j] = state.indices[j];
					break;
				}
			}
			for (unsigned int i = 0; i < resolved_formula_count; i++) {
				core::free(*resolved_formulas[i]);
				if (resolved_formulas[i]->reference_count == 0)
					core::free(resolved_formulas[i]);
			}
		}

		return new_proof;
	}

	template<typename... Args>
	Proof* add_formula_helper(Formula* formula, set_changes<Formula>& set_diff, unsigned int& new_constant, Args&&... args)
	{
		new_constant = 0;

		Formula* new_formula = preprocess_formula(formula);
		if (new_formula == nullptr) return nullptr;

		array_map<unsigned int, unsigned int> variable_map(16);
/* TODO: for debugging; delete this */
//print("add_formula_helper: Adding formula ", stderr); print(*new_formula, stderr, *debug_terminal_printer); print('\n', stderr);
		Formula* canonicalized = Canonicalizer::canonicalize(*new_formula, variable_map);
		core::free(*new_formula); if (new_formula->reference_count == 0) core::free(new_formula);
		if (canonicalized == nullptr) return nullptr;

/* TODO: for debugging; delete this */
//print_axioms(stderr);
//print("canonicalized: ", stderr); print(*canonicalized, stderr, *debug_terminal_printer); print('\n', stderr);
		array_map<unsigned int, Proof*> set_definitions(4);
		array_map<unsigned int, unsigned int> requested_set_sizes(4);
		Proof* new_proof = make_proof<false, true, true>(canonicalized, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
//print_axioms(stderr);
if (new_proof != nullptr) {
array_map<unsigned int, unsigned int> constant_map(1);
constant_map.put((unsigned int) built_in_predicates::UNKNOWN, new_constant);
Formula* expected_conclusion = relabel_constants(canonicalized, constant_map);
if (!check_proof<built_in_predicates, typename ProofCalculus::ProofCanonicalizer, ProofCalculus::Intuitionistic>(*new_proof, expected_conclusion)) {
check_proof<built_in_predicates, typename ProofCalculus::ProofCanonicalizer, ProofCalculus::Intuitionistic>(*new_proof, expected_conclusion);
print<built_in_predicates, typename ProofCalculus::ProofCanonicalizer, ProofCalculus::Intuitionistic>(*new_proof, stderr, *debug_terminal_printer);
fprintf(stderr, "add_formula WARNING: `check_proof` failed.\n");
}
core::free(*expected_conclusion); if (expected_conclusion->reference_count == 0) core::free(expected_conclusion);
}
		core::free(*canonicalized);
		if (canonicalized->reference_count == 0)
			core::free(canonicalized);
		return new_proof;
	}

	template<bool FreeProof = true>
	void remove_formula(Proof* proof, set_changes<Formula>& set_diff) {
		theory::changes changes;
		observations.remove(observations.index_of(proof));
		if (!get_theory_changes(*proof, changes)) return;
		subtract_changes(changes, set_diff);
		if (FreeProof) {
			core::free(*proof); if (proof->reference_count == 0) core::free(proof);
		}
	}

	template<ProofType Type>
	bool check_disjunction_introductions(const array<proof_node>& proof_steps) const
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
					if (*step == *entry.proof) {
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
			if (!computed_steps.contains(entry.proof)) {
				print("theory.check_proof_disjunctions WARNING: `proof_steps` "
						"contains a step that does not belong to a proof in `observations`.\n", stderr);
				print("  Formula: ", stderr); print(*entry.formula, stderr); print('\n', stderr);
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
		for (auto entry : constant_types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		} for (auto entry : constant_negated_types) {
			if (entry.value->reference_count <= 1) {
				print("concept.check_axioms WARNING: Found axiom '", stderr);
				print(*entry.value->formula, stderr);
				print("' with reference count less than 2.\n", stderr);
				success = false;
			}
		}
		return success;
	}

	template<typename... Printer>
	inline bool are_elements_provable(Printer&&... printer) {
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

			array<variable_assignment> possible_values(8);
			if (!::init(possible_values[0], sets.sets[i].arity))
				return false;
			possible_values.length = 1;
			default_prover prover(sets, implication_axioms);
			is_provable_without_abduction<false>(set_formula, quantifiers, possible_values, prover);

			/* count the number of concrete provable values in `possible_values` */
			for (const variable_assignment& values : possible_values) {
				/* if assignment is "concrete", check to make sure it's in the list of provable elements */
				tuple tup;
				::init(tup, values.assignment.length);
				bool is_concrete = true;
				for (uint_fast8_t k = 0; k < values.assignment.length; k++) {
					if (values.assignment.values[k].type == instantiation_type::CONSTANT) {
						tup.elements[k].type = tuple_element_type::CONSTANT;
						tup.elements[k].constant = values.assignment.values[k].constant;
					} else if (values.assignment.values[k].type == instantiation_type::NUMBER) {
						tup.elements[k].type = tuple_element_type::NUMBER;
						tup.elements[k].number = values.assignment.values[k].number;
					} else if (values.assignment.values[k].type == instantiation_type::STRING) {
						tup.elements[k].type = tuple_element_type::STRING;
						tup.elements[k].str = values.assignment.values[k].str;
					} else {
						is_concrete = false;
						break;
					}
				}
				if (is_concrete) {
					if (!sets.sets[i].provable_elements.contains(tup)) {
						print("theory.are_elements_provable ERROR: The set with ID ", stderr); print(i, stderr);
						print(" provably contains the element ", stderr); print(tup, stderr, *debug_terminal_printer);
						print(" but is not in `provable_elements`.\n", stderr);
						success = false;
					}
				}
				core::free(tup);
			}
			for (auto& element : possible_values) core::free(element);

			hash_set<unsigned int> ancestors(16);
			if (!sets.template get_ancestors<true>(i, ancestors))
				return false;
			for (unsigned int j = sets.sets[i].element_count(); j > 0; j--) {
				tuple_element* element = (tuple_element*) alloca(sizeof(tuple_element) * sets.sets[i].arity);
				const tuple_element* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * (j - 1));
				for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) {
					if (!::init(element[k], element_src[k])) {
						for (uint_fast8_t l = 0; l < k; l++) core::free(element[l]);
						return false;
					}
				}
				array<variable_assignment> possible_values(1);
				if (!::init(possible_values[0], sets.sets[i].arity)) {
					for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) core::free(element[k]);
					return false;
				}
				possible_values.length = 1;
				for (uint_fast8_t k = 0; k < possible_values[0].assignment.length; k++) {
					core::free(possible_values[0].assignment.values[k]);
					if (element[k].type == tuple_element_type::CONSTANT) {
						possible_values[0].assignment.values[k].type = instantiation_type::CONSTANT;
						possible_values[0].assignment.values[k].constant = element[k].constant;
					} else if (element[k].type == tuple_element_type::NUMBER) {
						possible_values[0].assignment.values[k].type = instantiation_type::NUMBER;
						possible_values[0].assignment.values[k].number = element[k].number;
					} else if (element[k].type == tuple_element_type::STRING) {
						possible_values[0].assignment.values[k].type = instantiation_type::STRING;
						if (!core::init(possible_values[0].assignment.values[k].str, element[k].str)) {
							possible_values[0].assignment.values[k].type = instantiation_type::CONSTANT;
							for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) core::free(element[k]);
							core::free(possible_values[0]); return false;
						}
					}
				}

				sets.remove_element_at(i, j - 1);
				default_prover prover(sets, implication_axioms);
				prover.h.set_ids[0] = i;
				prover.h.set_ids.length = 1;
				if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values, prover)) {
					print("theory.are_elements_provable ERROR: The element ", stderr);
					if (sets.sets[i].arity == 1)
						print(element[0], stderr, std::forward<Printer>(printer)...);
					else print(element, sets.sets[i].arity, stderr, std::forward<Printer>(printer)...);
					print(" does not provably belong to set with ID ", stderr);
					print(i, stderr); print(".\n", stderr);
					success = false;
				}
				for (auto& element : possible_values) core::free(element);
				sets.add_to_provable_elements(i, {element, sets.sets[i].arity}, ancestors);
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
		unsigned int arity = 1;
		const Formula* first_set_definition = ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.right->quantifier.operand;
		while (first_set_definition->type == FormulaType::LAMBDA) {
			first_set_definition = first_set_definition->quantifier.operand;
			arity++;
		}
		if (sets.set_ids.table.contains(*first_set_definition))
			return true;

		Term* atom = Term::new_constant(constant);
		for (unsigned int i = 0; i < arity; i++) {
			Term* temp = Formula::new_apply(atom, Formula::new_variable(i + 1));
			if (temp == nullptr) {
				core::free(*atom); core::free(atom);
				return false;
			}
			atom = temp;
		}

		bool contains;
		pair<array<instance>, array<instance>>& type_instances = atoms.get(*atom, contains);
		core::free(*atom); core::free(atom);
		if (contains && type_instances.key.length != 0)
			return true;
		return false;
	}

	template<unsigned int Arg>
	Term* get_arg(unsigned int event_constant) const {
		bool contains;
		Proof* proof = ground_concepts[event_constant - new_constant_offset].function_values.get(Arg, contains);
		if (contains) return proof->formula->binary.right;

		Term* definition = Term::new_apply(&Constants<Arg>::value, Term::new_constant(event_constant));
		if (definition == nullptr)
			return nullptr;
		Constants<Arg>::value.reference_count++;

		unsigned int constant = reverse_definitions.get(*definition, contains);
		core::free(*definition); core::free(definition);
		if (!contains) return nullptr;
		return ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.left;
	}

	inline bool is_name_event(unsigned int event_constant) const {
		if (event_constant < new_constant_offset
		 || ground_concepts[event_constant - new_constant_offset].types.keys == nullptr)
			return false;
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
		if (atom.binary.right->type != TermType::CONSTANT && atom.binary.right->type != TermType::NUMBER)
			fprintf(stderr, "add_unary_atom WARNING: The operand of this application is not a constant or number.\n");
#endif

		/* check to make sure the args are type-correct */
		if (atom.binary.right->type == TermType::CONSTANT) {
			unsigned int arg = atom.binary.right->constant;
			if (atom.binary.left->type == TermType::CONSTANT && atom.binary.left->constant == (unsigned int) built_in_predicates::NAME) {
				/* arg1 must be a non-name, non-set constant and arg2 must be a string */
				Term* arg1 = get_arg<(unsigned int) built_in_predicates::ARG1>(arg);
				if (arg1 != nullptr) {
					if (arg1->type != TermType::CONSTANT || is_name_event(arg1->constant) || is_provably_a_set(arg1->constant))
						return false;
				}
				Term* arg2 = get_arg<(unsigned int) built_in_predicates::ARG2>(arg);
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
		} else {
			if (atom.binary.left->type == TermType::CONSTANT && atom.binary.left->constant == (unsigned int) built_in_predicates::NAME)
				return false;
		}

		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr) {
			if (atom.binary.right->type == TermType::CONSTANT)
				try_free_concept_id(atom.binary.right->constant);
			return false;
		}
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;

		bool contains; unsigned int bucket;
		if (!atoms.check_size()) {
			core::free(*lifted_atom); core::free(lifted_atom);
			if (atom.binary.right->type == TermType::CONSTANT)
				try_free_concept_id(atom.binary.right->constant);
			return false;
		}

		pair<array<instance>, array<instance>>& instance_pair = atoms.get(*lifted_atom, contains, bucket);
		if (!contains) {
			if (!array_init(instance_pair.key, 8)) {
				core::free(*lifted_atom); core::free(lifted_atom);
				if (atom.binary.right->type == TermType::CONSTANT)
					try_free_concept_id(atom.binary.right->constant);
				return false;
			} else if (!array_init(instance_pair.value, 8)) {
				core::free(instance_pair.key);
				core::free(*lifted_atom); core::free(lifted_atom);
				if (atom.binary.right->type == TermType::CONSTANT)
					try_free_concept_id(atom.binary.right->constant);
				return false;
			}
			atoms.table.keys[bucket] = *lifted_atom;
			atoms.table.size++;
		}

		array<instance>& instances = (Negated ? instance_pair.value : instance_pair.key);
		array_map<Term, Proof*>& ground_types = (atom.binary.right->type == TermType::CONSTANT
				? (Negated ? ground_concepts[atom.binary.right->constant - new_constant_offset].negated_types : ground_concepts[atom.binary.right->constant - new_constant_offset].types)
				: (Negated ? constant_negated_types : constant_types));
		if (!instances.ensure_capacity(instances.length + 1)
		 || !ground_types.ensure_capacity(ground_types.size + 1))
		{
			core::free(*lifted_atom); core::free(lifted_atom);
			if (atom.binary.right->type == TermType::CONSTANT)
				try_free_concept_id(atom.binary.right->constant);
			return false;
		}

		if (atom.binary.right->type == TermType::CONSTANT)
			add_sorted<false>(instances, instance_constant(atom.binary.right->constant));
		else add_sorted<false>(instances, instance_number(atom.binary.right->number));
		ground_types.keys[ground_types.size] = (atom.binary.right->type == TermType::CONSTANT ? *lifted_atom : atom);
		ground_types.values[ground_types.size++] = axiom;
		axiom->reference_count++;
		if (atom.binary.right->type == TermType::CONSTANT)
			ground_axiom_count++;
		core::free(*lifted_atom); core::free(lifted_atom);

		if (!check_set_membership_after_addition<ResolveInconsistencies>(&atom, std::forward<Args>(visitor)...)) {
			remove_unary_atom<Negated>(atom, std::forward<Args>(visitor)...);
			return false;
		}
		return true;
	}

	template<bool Negated, bool ResolveInconsistencies, typename... Args>
	inline bool add_binary_atom(relation rel, Proof* axiom, Args&&... visitor)
	{
		try_init_concept(rel.arg1);
		try_init_concept(rel.arg2);

		bool contains; unsigned int bucket;
		if (!relations.check_size(relations.table.size + 4)) return false;
		pair<array<instance>, array<instance>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(predicate_instance_pair.key, 8)) {
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1);
				return false;
			} else if (!array_init(predicate_instance_pair.value, 8)) {
				core::free(predicate_instance_pair.key);
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1); return false;
			}
			relations.table.keys[bucket] = {0, rel.arg1, rel.arg2};
			relations.table.size++;
		}

		pair<array<instance>, array<instance>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2}, contains, bucket);
		if (!contains) {
			if (!array_init(arg1_instance_pair.key, 8)) {
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1);
				return false;
			} else if (!array_init(arg1_instance_pair.value, 8)) {
				core::free(arg1_instance_pair.key);
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, 0, rel.arg2};
			relations.table.size++;
		}

		pair<array<instance>, array<instance>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0}, contains, bucket);
		if (!contains) {
			if (!array_init(arg2_instance_pair.key, 8)) {
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1);
				return false;
			} else if (!array_init(arg2_instance_pair.value, 8)) {
				core::free(arg2_instance_pair.key);
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1); return false;
			}
			relations.table.keys[bucket] = {rel.predicate, rel.arg1, 0};
			relations.table.size++;
		}

		array<instance>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<instance>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<instance>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);
		if (!predicate_instances.ensure_capacity(predicate_instances.length + 1)
		 || !arg1_instances.ensure_capacity(arg1_instances.length + 1)
		 || !arg2_instances.ensure_capacity(arg2_instances.length + 1)
		 || !ground_arg1.ensure_capacity(ground_arg1.size + 3)
		 || !ground_arg2.ensure_capacity(ground_arg2.size + 3))
		{
			try_free_concept_id(rel.arg2);
			try_free_concept_id(rel.arg1);
			return false;
		}

		if (rel.arg1 == rel.arg2) {
			pair<array<instance>, array<instance>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0}, contains, bucket);
			if (!contains) {
				if (!array_init(both_arg_instance_pair.key, 8)) {
					try_free_concept_id(rel.arg2);
					try_free_concept_id(rel.arg1);
					return false;
				} else if (!array_init(both_arg_instance_pair.value, 8)) {
					core::free(both_arg_instance_pair.key);
					try_free_concept_id(rel.arg2);
					try_free_concept_id(rel.arg1);
					return false;
				}
				relations.table.keys[bucket] = {rel.predicate, 0, 0};
				relations.table.size++;
			}

			array<instance>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			if (!both_arg_instances.ensure_capacity(both_arg_instances.length + 1)) {
				try_free_concept_id(rel.arg2);
				try_free_concept_id(rel.arg1);
				return false;
			}
			both_arg_instances[both_arg_instances.length++] = instance_constant(rel.arg1);
			insertion_sort(both_arg_instances);
		}

		add_sorted<false>(predicate_instances, instance_constant(rel.predicate));
		add_sorted<false>(arg1_instances, instance_constant(rel.arg1));
		add_sorted<false>(arg2_instances, instance_constant(rel.arg2));
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
		if (atom == nullptr) {
			try_free_concept_id(rel.arg2);
			try_free_concept_id(rel.arg1);
			return false;
		}
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
		if (atom.binary.right->type != TermType::CONSTANT && atom.binary.right->type != TermType::NUMBER)
			fprintf(stderr, "remove_unary_atom WARNING: The operand of this application is not a constant or number.\n");
#endif

		Formula* lifted_atom = Term::new_apply(atom.binary.left, &Variables<1>::value);
		if (lifted_atom == nullptr)
			return false;
		atom.binary.left->reference_count++;
		Variables<1>::value.reference_count++;

#if !defined(NDEBUG)
		if (!atoms.table.contains(*lifted_atom)) {
			print("theory.remove_unary_atom WARNING: `atoms` does not contain the key ", stderr);
			print(*lifted_atom, stderr); print(".\n", stderr);
		} if (atom.binary.right->type == TermType::CONSTANT && (new_constant_offset > atom.binary.right->constant || ground_concepts[atom.binary.right->constant - new_constant_offset].types.keys == NULL))
			fprintf(stderr, "theory.remove_unary_atom WARNING: `ground_concepts` does not contain the key %u.\n", atom.binary.right->constant);
#endif

		array<instance>& instances = (Negated ? atoms.get(*lifted_atom).value : atoms.get(*lifted_atom).key);
		array_map<Term, Proof*>& ground_types = (atom.binary.right->type == TermType::CONSTANT
				? (Negated ? ground_concepts[atom.binary.right->constant - new_constant_offset].negated_types : ground_concepts[atom.binary.right->constant - new_constant_offset].types)
				: (Negated ? constant_negated_types : constant_types));

		unsigned int index = (atom.binary.right->type == TermType::CONSTANT ? index_of_constant(instances, atom.binary.right->constant) : index_of_number(instances, atom.binary.right->number));
#if !defined(NDEBUG)
		if (index == instances.length) {
			fprintf(stderr, "theory.remove_unary_atom WARNING: `instances` does not contain ");
			if (atom.binary.right->type == TermType::CONSTANT) {
				fprintf(stderr, "constant %u.\n", atom.binary.right->constant);
			} else {
				print("number ", stderr);
				print(atom.binary.right->number.integer, stderr);
				print(',', stderr);
				print(atom.binary.right->number.decimal, stderr);
				print(".\n", stderr);
			}
		}
#endif
		shift_left(instances.data + index, instances.length - index - 1);
		instances.length--;

		index = (atom.binary.right->type == TermType::CONSTANT ? ground_types.index_of(*lifted_atom) : ground_types.index_of(atom));
#if !defined(NDEBUG)
		if (index == ground_types.size) {
			print("theory.remove_unary_atom WARNING: `ground_types` does not contain ", stderr);
			print(*lifted_atom, stderr); print(".\n", stderr);
		}
#endif
		Proof* axiom = ground_types.values[index];
		core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		core::free(ground_types.keys[index]);
		ground_types.remove_at(index);
		if (atom.binary.right->type == TermType::CONSTANT)
			ground_axiom_count--;
		core::free(*lifted_atom); core::free(lifted_atom);

		return check_set_membership_after_subtraction(&atom, 0, std::forward<Args>(visitor)...);
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

		pair<array<instance>, array<instance>>& predicate_instance_pair = relations.get({0, rel.arg1, rel.arg2});
		pair<array<instance>, array<instance>>& arg1_instance_pair = relations.get({rel.predicate, 0, rel.arg2});
		pair<array<instance>, array<instance>>& arg2_instance_pair = relations.get({rel.predicate, rel.arg1, 0});

		array<instance>& predicate_instances = (Negated ? predicate_instance_pair.value : predicate_instance_pair.key);
		array<instance>& arg1_instances = (Negated ? arg1_instance_pair.value : arg1_instance_pair.key);
		array<instance>& arg2_instances = (Negated ? arg2_instance_pair.value : arg2_instance_pair.key);
		array_map<relation, Proof*>& ground_arg1 = (Negated ? ground_concepts[rel.arg1 - new_constant_offset].negated_relations : ground_concepts[rel.arg1 - new_constant_offset].relations);
		array_map<relation, Proof*>& ground_arg2 = (Negated ? ground_concepts[rel.arg2 - new_constant_offset].negated_relations : ground_concepts[rel.arg2 - new_constant_offset].relations);

		unsigned int index;
		if (rel.arg1 == rel.arg2) {
			pair<array<instance>, array<instance>>& both_arg_instance_pair = relations.get({rel.predicate, 0, 0});
			array<instance>& both_arg_instances = (Negated ? both_arg_instance_pair.value : both_arg_instance_pair.key);
			index = index_of_constant(both_arg_instances, rel.arg1);
#if !defined(NDEBUG)
			if (index == both_arg_instances.length)
				fprintf(stderr, "theory.remove_binary_atom WARNING: `both_arg_instances` does not contain %u.\n", rel.arg1);
#endif
			shift_left(both_arg_instances.data + index, both_arg_instances.length - index - 1);
			both_arg_instances.length--;
		}

		index = index_of_constant(predicate_instances, rel.predicate);
#if !defined(NDEBUG)
		if (index == predicate_instances.length)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `predicate_instances` does not contain %u.\n", rel.predicate);
#endif
		shift_left(predicate_instances.data + index, predicate_instances.length - index - 1);
		predicate_instances.length--;

		index = index_of_constant(arg1_instances, rel.arg1);
#if !defined(NDEBUG)
		if (index == arg1_instances.length)
			fprintf(stderr, "theory.remove_binary_atom WARNING: `arg1_instances` does not contain %u.\n", rel.arg1);
#endif
		shift_left(arg1_instances.data + index, arg1_instances.length - index - 1);
		arg1_instances.length--;

		index = index_of_constant(arg2_instances, rel.arg2);
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
		bool result = check_set_membership_after_subtraction(atom, 0, std::forward<Args>(visitor)...);
		core::free(*atom); core::free(atom);
		try_free_concept_id(rel.arg2);
		try_free_concept_id(rel.arg1);
		return result;
	}

	unsigned int index_of(const array<proof_node>& elements, Proof* proof) const {
		unsigned int index = elements.length;
		for (unsigned int i = 0; i < elements.length; i++)
			if (elements[i].proof == proof) { index = i; break; }
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
		IMPLICATION_AXIOM,
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
			proof_node intro_node;
		};
		required_set_size requested_set_size;
		required_set_size consequent_set_size;

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
			case change_type::IMPLICATION_AXIOM:
			case change_type::SET_SIZE_AXIOM:
			case change_type::DEFINITION:
			case change_type::FUNCTION_VALUE:
				core::free(*c.axiom); if (c.axiom->reference_count == 0) core::free(c.axiom);
				return;
			case change_type::IMPLICATION_INTRO_NODE:
			case change_type::NEGATED_CONJUNCTION_NODE:
			case change_type::EXISTENTIAL_INTRO_NODE:
			case change_type::DISJUNCTION_INTRO_NODE:
				core::free(*c.intro_node.formula); if (c.intro_node.formula->reference_count == 0) core::free(c.intro_node.formula);
				core::free(*c.intro_node.proof); if (c.intro_node.proof->reference_count == 0) core::free(c.intro_node.proof);
				core::free(c.intro_node.set_definitions);
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

		inline bool add(change_type type, Proof* axiom, required_set_size req_set_size) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length].requested_set_size = req_set_size;
			list[list.length++].axiom = axiom;
			axiom->reference_count++;
			return true;
		}

		inline bool add(change_type type, Proof* axiom,
						required_set_size antecedent_set_size,
						required_set_size consequent_set_size)
		{
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length].requested_set_size = antecedent_set_size;
			list[list.length].consequent_set_size = consequent_set_size;
			list[list.length++].axiom = axiom;
			axiom->reference_count++;
			return true;
		}

		inline bool add(change_type type, const proof_node& intro_node) {
			if (!list.ensure_capacity(list.length + 1)) return false;
			list[list.length].type = type;
			list[list.length].intro_node.formula = intro_node.formula;
			list[list.length].intro_node.proof = intro_node.proof;
			if (!array_map_init(list[list.length].intro_node.set_definitions, intro_node.set_definitions.capacity))
				return false;
			for (unsigned int i = 0; i < intro_node.set_definitions.size; i++) {
				list[list.length].intro_node.set_definitions.keys[i] = intro_node.set_definitions.keys[i];
				list[list.length].intro_node.set_definitions.values[i] = intro_node.set_definitions.values[i];
			}
			list[list.length++].intro_node.set_definitions.size = intro_node.set_definitions.size;
			intro_node.formula->reference_count++;
			intro_node.proof->reference_count++;
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
		for (unsigned int i = changes.list.length; i > 0; i--) {
			change& c = changes.list[i - 1];
			array_multiset<unsigned int> constants(8);
			switch (c.type) {
			case change_type::UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<false, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  UNARY_ATOM: ", i);
					print(c.unary_atom.key, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::NEGATED_UNARY_ATOM:
				if (!get_constants(c.unary_atom.key, constants)) return false;
				for (unsigned int i = 0; i < constants.counts.size; i++)
					if (!try_init_concept(constants.counts.keys[i])) return false;
				if (!add_unary_atom<true, false>(c.unary_atom.key, c.unary_atom.value, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  NEGATED_UNARY_ATOM: ~", i);
					print(c.unary_atom.key, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate, 2)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<false, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  BINARY_ATOM: ", i);
					print(c.binary_atom.key.predicate, stderr, *debug_terminal_printer); print('(', stderr);
					print(c.binary_atom.key.arg1, stderr, *debug_terminal_printer); print(',', stderr);
					print(c.binary_atom.key.arg2, stderr, *debug_terminal_printer); print(").\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::NEGATED_BINARY_ATOM:
				if (!try_init_concept(c.binary_atom.key.predicate, 2)
				 || !try_init_concept(c.binary_atom.key.arg1) || !try_init_concept(c.binary_atom.key.arg2)
				 || !add_binary_atom<true, false>(c.binary_atom.key, c.binary_atom.value, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  NEGATED_BINARY_ATOM: ~", i);
					print(c.binary_atom.key.predicate, stderr, *debug_terminal_printer); print('(', stderr);
					print(c.binary_atom.key.arg1, stderr, *debug_terminal_printer); print(',', stderr);
					print(c.binary_atom.key.arg2, stderr, *debug_terminal_printer); print(").\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::SUBSET_AXIOM:
				{
					unsigned int arity = 1;
					Formula* operand = c.axiom->formula->quantifier.operand;
					while (operand->type == FormulaType::FOR_ALL) {
						operand = operand->quantifier.operand;
						arity++;
					}
					Formula* left = operand->binary.left;
					Formula* right = operand->binary.right;

					if (!get_constants(*operand, constants)) return false;
					for (unsigned int i = 0; i < constants.counts.size; i++)
						if (!try_init_concept(constants.counts.keys[i])) return false;

					unsigned int antecedent_set, consequent_set;
					bool is_antecedent_new, is_consequent_new;
					Proof* axiom = get_subset_axiom_with_required_set_size<false, false>(
							c.axiom, left, right, arity, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new,
							c.requested_set_size, c.consequent_set_size, set_diff, std::forward<Args>(visitor)...);
					if (axiom != c.axiom) {
#ifndef NDEBUG
						fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  SUBSET_AXIOM: ", i);
						print(*c.axiom->formula, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
						print_axioms<true>(stdout, *debug_terminal_printer);
#endif
						return false;
					}
				}
				continue;
			case change_type::IMPLICATION_AXIOM:
				if (!implication_axioms.ensure_capacity(implication_axioms.length + 1))
					return false;
				implication_axioms[implication_axioms.length].key = c.axiom;
				implication_axioms[implication_axioms.length++].value = false;
				c.axiom->reference_count++;
				if (!check_new_implication_satisfaction<false>(implication_axioms.length - 1, set_diff, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  IMPLICATION_AXIOM: ", i);
					print(*c.axiom->formula, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					implication_axioms.length--;
					c.axiom->reference_count--;
					return false;
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
					if (!get_constants(*set_formula, constants)) return false;
					for (unsigned int i = 0; i < constants.counts.size; i++)
						if (!try_init_concept(constants.counts.keys[i])) return false;
					bool is_new;
					if (!sets.get_set_id(set_formula, arity, set_id, is_new, set_diff, std::forward<Args>(visitor)...)
					 || !sets.sets[set_id].set_size_axiom(c.axiom, set_diff, std::forward<Args>(visitor)...)
					 || (is_new && !check_new_set_membership<false>(set_id, std::forward<Args>(visitor)...))) {
#ifndef NDEBUG
						fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  SET_SIZE_AXIOM: ", i);
						print(*c.axiom->formula, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
						print_axioms<true>(stdout, *debug_terminal_printer);
#endif
						return false;
					}
				}
				continue;
			case change_type::DEFINITION:
				if (!add_definition<false>(c.axiom, c.requested_set_size, set_diff, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
					fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  DEFINITION: ", i);
					print(*c.axiom->formula, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
					print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::FUNCTION_VALUE:
				if (!add_function_value<false>(c.axiom, std::forward<Args>(visitor)...)) {
#ifndef NDEBUG
						fprintf(stderr, "theory.add_changes ERROR: Failed to add change %u;\n  FUNCTION_VALUE: ", i);
						print(*c.axiom->formula, stderr, *debug_terminal_printer); print(".\nCurrent theory:\n", stderr);
						print_axioms<true>(stdout, *debug_terminal_printer);
#endif
					return false;
				}
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				if (!implication_intro_nodes.ensure_capacity(implication_intro_nodes.length + 1))
					return false;
				implication_intro_nodes[implication_intro_nodes.length].formula = c.intro_node.formula;
				implication_intro_nodes[implication_intro_nodes.length].proof = c.intro_node.proof;
				if (!array_map_init(implication_intro_nodes[implication_intro_nodes.length].set_definitions, c.intro_node.set_definitions.capacity))
					return false;
				for (unsigned int i = 0; i < c.intro_node.set_definitions.size; i++) {
					implication_intro_nodes[implication_intro_nodes.length].set_definitions.keys[i] = c.intro_node.set_definitions.keys[i];
					implication_intro_nodes[implication_intro_nodes.length].set_definitions.values[i] = c.intro_node.set_definitions.values[i];
				}
				implication_intro_nodes[implication_intro_nodes.length++].set_definitions.size = c.intro_node.set_definitions.size;
				c.intro_node.formula->reference_count++;
				continue;
			case change_type::NEGATED_CONJUNCTION_NODE:
				if (!negated_conjunction_nodes.ensure_capacity(negated_conjunction_nodes.length + 1))
					return false;
				negated_conjunction_nodes[negated_conjunction_nodes.length].formula = c.intro_node.formula;
				negated_conjunction_nodes[negated_conjunction_nodes.length].proof = c.intro_node.proof;
				if (!array_map_init(negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions, c.intro_node.set_definitions.capacity))
					return false;
				for (unsigned int i = 0; i < c.intro_node.set_definitions.size; i++) {
					negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions.keys[i] = c.intro_node.set_definitions.keys[i];
					negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions.values[i] = c.intro_node.set_definitions.values[i];
				}
				negated_conjunction_nodes[negated_conjunction_nodes.length++].set_definitions.size = c.intro_node.set_definitions.size;
				c.intro_node.formula->reference_count++;
				continue;
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					if (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1))
						return false;
					existential_intro_nodes[existential_intro_nodes.length].formula = c.intro_node.formula;
					existential_intro_nodes[existential_intro_nodes.length].proof = c.intro_node.proof;
					if (!array_map_init(existential_intro_nodes[existential_intro_nodes.length].set_definitions, c.intro_node.set_definitions.capacity))
						return false;
					for (unsigned int i = 0; i < c.intro_node.set_definitions.size; i++) {
						existential_intro_nodes[existential_intro_nodes.length].set_definitions.keys[i] = c.intro_node.set_definitions.keys[i];
						existential_intro_nodes[existential_intro_nodes.length].set_definitions.values[i] = c.intro_node.set_definitions.values[i];
					}
					existential_intro_nodes[existential_intro_nodes.length++].set_definitions.size = c.intro_node.set_definitions.size;

					Term* term;
					if (c.intro_node.proof->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.proof->operands[2]->term;
					else term = c.intro_node.proof->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						if (!try_init_concept(term->constant)
						 || !ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.add(c.intro_node.proof))
						{
							core::free(existential_intro_nodes[--existential_intro_nodes.length]);
							return false;
						}
					}
					c.intro_node.formula->reference_count++;
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				if (!disjunction_intro_nodes.ensure_capacity(disjunction_intro_nodes.length + 1))
					return false;
				disjunction_intro_nodes[disjunction_intro_nodes.length].formula = c.intro_node.formula;
				disjunction_intro_nodes[disjunction_intro_nodes.length].proof = c.intro_node.proof;
				if (!array_map_init(disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions, c.intro_node.set_definitions.capacity))
					return false;
				for (unsigned int i = 0; i < c.intro_node.set_definitions.size; i++) {
					disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions.keys[i] = c.intro_node.set_definitions.keys[i];
					disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions.values[i] = c.intro_node.set_definitions.values[i];
				}
				disjunction_intro_nodes[disjunction_intro_nodes.length++].set_definitions.size = c.intro_node.set_definitions.size;
				c.intro_node.formula->reference_count++;
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
		for (unsigned int i = 0; i < changes.list.length; i++) {
			change& c = changes.list[i];
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
					if (!get_constants(*c.axiom->formula->quantifier.operand, constants)) return false;
					for (unsigned int i = 0; i < constants.counts.size; i++)
						try_free_concept_id(constants.counts.keys[i]);
					unsigned int child_count = c.axiom->children.length;
					if (changes.list.length == 1) child_count++;
					if (c.axiom->reference_count == child_count + 2) {
						Formula* operand = c.axiom->formula->quantifier.operand;
						unsigned int arity = 1;
						while (operand->type == FormulaType::FOR_ALL) {
							operand = operand->quantifier.operand;
							arity++;
						}
						Formula* antecedent = operand->binary.left;
						Formula* consequent = operand->binary.right;
						unsigned int antecedent_set = sets.get_existing_set_id(antecedent, arity);
						unsigned int consequent_set = sets.get_existing_set_id(consequent, arity);

						free_subset_axiom_with_required_set_size<false>(antecedent, consequent, arity, antecedent_set, consequent_set, c.requested_set_size, c.consequent_set_size, freeable_set_size_axioms, set_diff, std::forward<Args>(visitor)...);
					} else {
						core::free(*c.axiom);
					}
					continue;
				}
			case change_type::IMPLICATION_AXIOM:
				for (unsigned int i = 0; i < implication_axioms.length; i++) {
					if (implication_axioms[i].key == c.axiom) {
						check_old_implication_satisfaction(i, set_diff, std::forward<Args>(visitor)...);
						core::free(*implication_axioms[i].key);
						if (implication_axioms[i].key->reference_count == 0)
							core::free(implication_axioms[i].key);
						implication_axioms.remove(i);
						break;
					}
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
					if (!sets.get_set_id(set_formula, arity, set_id)
					 || !freeable_set_size_axioms.add(c.axiom)) return false;
					unsigned int old_ref_count = c.axiom->reference_count;
					c.axiom->reference_count = 1;
					bool is_freeable = sets.is_freeable(set_id);
					c.axiom->reference_count = old_ref_count;
					if (is_freeable) {
						check_old_set_membership(set_id, std::forward<Args>(visitor)...);
						for (Proof* size_axiom : sets.sets[set_id].size_axioms)
							on_old_size_axiom(size_axiom, set_diff, std::forward<Args>(visitor)...);
						unsigned int min_set_size; unsigned int max_set_size;
						sets.set_size_bounds(set_id, min_set_size, max_set_size);
						on_free_set(set_id, sets, min_set_size, max_set_size, set_diff, std::forward<Args>(visitor)...);
						sets.free_set(set_id);
					}
					continue;
				}
			case change_type::DEFINITION:
				remove_definition(c.axiom, c.requested_set_size, freeable_set_size_axioms, set_diff, std::forward<Args>(visitor)...);
				continue;
			case change_type::FUNCTION_VALUE:
				remove_function_value(c.axiom, std::forward<Args>(visitor)...);
				continue;
			case change_type::IMPLICATION_INTRO_NODE:
				{
					unsigned int index = index_of(implication_intro_nodes, c.intro_node.proof);
					implication_intro_nodes[index].formula->reference_count--;
					core::free(implication_intro_nodes[index].set_definitions);
					implication_intro_nodes.remove(index);
					on_undo_filter_operands(c.intro_node.formula, std::forward<Args>(visitor)...);
					continue;
				}
			case change_type::NEGATED_CONJUNCTION_NODE:
				{
					unsigned int index = index_of(negated_conjunction_nodes, c.intro_node.proof);
					negated_conjunction_nodes[index].formula->reference_count--;
					core::free(negated_conjunction_nodes[index].set_definitions);
					negated_conjunction_nodes.remove(index);
					on_undo_filter_operands(c.intro_node.formula, std::forward<Args>(visitor)...);
					continue;
				}
			case change_type::EXISTENTIAL_INTRO_NODE:
				{
					Term* term;
					if (c.intro_node.proof->type == ProofType::EXISTENTIAL_INTRODUCTION)
						term = c.intro_node.proof->operands[2]->term;
					else term = c.intro_node.proof->operands[0]->operands[0]->operands[1]->term;
					if (term->type == TermType::CONSTANT && term->constant >= new_constant_offset) {
						unsigned int index = ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.index_of(c.intro_node.proof);
						ground_concepts[term->constant - new_constant_offset].existential_intro_nodes.remove(index);
						try_free_concept_id(term->constant);
					}

					unsigned int index = index_of(existential_intro_nodes, c.intro_node.proof);
					existential_intro_nodes[index].formula->reference_count--;
					on_undo_filter_constants(*this, c.intro_node.formula->quantifier.operand, term, c.intro_node.formula->quantifier.variable, existential_intro_nodes[index].set_definitions, std::forward<Args>(visitor)...);
					core::free(existential_intro_nodes[index].set_definitions);
					existential_intro_nodes.remove(index);
					continue;
				}
			case change_type::DISJUNCTION_INTRO_NODE:
				{
					unsigned int index = index_of(disjunction_intro_nodes, c.intro_node.proof);
					disjunction_intro_nodes[index].formula->reference_count--;
					core::free(disjunction_intro_nodes[index].set_definitions);
					disjunction_intro_nodes.remove(index);
					on_undo_filter_operands(c.intro_node.formula, std::forward<Args>(visitor)...);
					continue;
				}
			}
			fprintf(stderr, "theory.subtract_changes ERROR: Unrecognized `change_type`.\n");
			return false;
		}
		return true;
	}

	inline bool is_set_size_proof(Proof* proof, unsigned int& constant, unsigned int& set_size) const
	{
		if (proof->type == ProofType::AXIOM && proof->formula->type == FormulaType::EQUALS
		 && proof->formula->binary.left->type == TermType::UNARY_APPLICATION
		 && proof->formula->binary.left->binary.left->type == TermType::CONSTANT
		 && proof->formula->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
		 && proof->formula->binary.left->binary.right->type == TermType::CONSTANT
		 && proof->formula->binary.right->type == TermType::NUMBER
		 && proof->formula->binary.right->number.integer >= 0
		 && proof->formula->binary.right->number.decimal == 0)
		{
			constant = proof->formula->binary.left->binary.right->constant;
			set_size = proof->formula->binary.right->number.integer;
			return true;
		} else if (proof->type == ProofType::EQUALITY_ELIMINATION && proof->operands[0]->type == ProofType::AXIOM
				&& proof->operands[0]->formula->type == FormulaType::EQUALS
				&& proof->operands[0]->formula->binary.left->type == TermType::CONSTANT
				&& proof->operands[0]->formula->binary.right->type == TermType::LAMBDA
				&& proof->operands[1]->type == ProofType::AXIOM
				&& proof->operands[1]->formula->type == FormulaType::EQUALS
				&& proof->operands[1]->formula->binary.left->type == TermType::UNARY_APPLICATION
				&& proof->operands[1]->formula->binary.left->binary.left->type == TermType::CONSTANT
				&& proof->operands[1]->formula->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
				&& proof->operands[1]->formula->binary.left->binary.right->type == TermType::LAMBDA
				&& proof->operands[1]->formula->binary.right->type == TermType::NUMBER
				&& proof->operands[1]->formula->binary.right->number.integer >= 0
				&& proof->operands[1]->formula->binary.right->number.decimal == 0)
		{
			constant = proof->operands[0]->formula->binary.left->constant;
			set_size = proof->operands[1]->formula->binary.right->number.integer;
			return true;
		}
		return false;
	}

	template<typename... Visitor>
	bool get_theory_changes(
			Proof& proof, array<const Proof*>& discharged_axioms,
			array_map<const Proof*, unsigned int>& reference_counts,
			array_map<unsigned int, unsigned int>& requested_set_sizes,
			theory::changes& changes, Visitor&&... visitor) const
	{
		visit_node(proof, std::forward<Visitor>(visitor)...);
#if !defined(NDEBUG)
		bool contains;
		unsigned int reference_count = reference_counts.get(&proof, contains);
		if (!contains) fprintf(stderr, "theory.get_theory_changes WARNING: The given proof is not in the map `reference_counts`.\n");
#else
		unsigned int reference_count = reference_counts.get(&proof);
#endif

		unsigned int old_discharged_axiom_count = discharged_axioms.length;
		unsigned int old_requested_set_size_count = requested_set_sizes.size;

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

				Formula* operand = proof.formula->quantifier.operand;
				while (operand->type == FormulaType::FOR_ALL)
					operand = operand->quantifier.operand;
				Formula* left = operand->binary.left;
				Formula* right = operand->binary.right;

				required_set_size antecedent_set_size = {0, UINT_MAX};
				required_set_size consequent_set_size = {0, UINT_MAX};
				if (left->type == FormulaType::UNARY_APPLICATION
				 && left->binary.left->type == FormulaType::CONSTANT
				 && left->binary.right->type == FormulaType::VARIABLE)
				{
					unsigned int index = requested_set_sizes.index_of(left->binary.left->constant);
					if (index < requested_set_sizes.size) {
						antecedent_set_size.min_set_size = requested_set_sizes.values[index];
						antecedent_set_size.max_set_size = requested_set_sizes.values[index];
					}
				}
				if (right->type == FormulaType::UNARY_APPLICATION
				 && right->binary.left->type == FormulaType::CONSTANT
				 && right->binary.right->type == FormulaType::VARIABLE)
				{
					unsigned int index = requested_set_sizes.index_of(right->binary.left->constant);
					if (index < requested_set_sizes.size) {
						consequent_set_size.min_set_size = requested_set_sizes.values[index];
						consequent_set_size.max_set_size = requested_set_sizes.values[index];
					}
				}

				return changes.add(change_type::SUBSET_AXIOM, &proof, antecedent_set_size, consequent_set_size);

			} else if (proof.formula->type == FormulaType::IF_THEN) {
				if (reference_count != 1) return true;
				return changes.add(change_type::IMPLICATION_AXIOM, &proof);
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

					required_set_size requested_set_size = {0, UINT_MAX};
					unsigned int index = requested_set_sizes.index_of(proof.formula->binary.left->constant);
					if (index < requested_set_sizes.size) {
						requested_set_size.min_set_size = requested_set_sizes.values[index];
						requested_set_size.max_set_size = requested_set_sizes.values[index];
					}
					return changes.add(change_type::DEFINITION, &proof, requested_set_size);
				} else if (proof.formula->binary.left->type == TermType::UNARY_APPLICATION) {
					if (reference_count != 1) return true;
					return changes.add(change_type::FUNCTION_VALUE, &proof);
				} else {
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected equals axiom.\n");
				}
			} else {
				bool is_built_in = (&proof == empty_set_axiom);
				for (unsigned int i = 0; !is_built_in && i < built_in_axioms.length; i++)
					is_built_in = (&proof == built_in_axioms[i]);
				if (!is_built_in)
					fprintf(stderr, "get_theory_changes WARNING: Found unexpected axiom.\n");
			}
			break;

		case ProofType::IMPLICATION_INTRODUCTION:
		{
			if (reference_count != 0) return true;
			unsigned int index = index_of(implication_intro_nodes, &proof);
			if (!changes.add(change_type::IMPLICATION_INTRO_NODE, implication_intro_nodes[index])
			 || !discharged_axioms.add(proof.operands[1]))
				return false;
			break;
		}

		case ProofType::PROOF_BY_CONTRADICTION:
			if (reference_count != 0) return true;
			discharged_axioms.add(proof.operands[1]);
			if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
			  && proof.operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
			  && proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				unsigned int index = index_of(negated_conjunction_nodes, &proof);
				if (!visit_negated_conjunction(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::NEGATED_CONJUNCTION_NODE, negated_conjunction_nodes[index])) return false;
			} else if (proof.operands[0]->type == ProofType::NEGATION_ELIMINATION
					&& proof.operands[0]->operands[0]->type == ProofType::UNIVERSAL_ELIMINATION
					&& proof.operands[0]->operands[0]->operands[0] == proof.operands[1])
			{
				unsigned int index = index_of(existential_intro_nodes, &proof);
				if (!visit_negated_universal_intro(std::forward<Visitor>(visitor)...)
				 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, existential_intro_nodes[index])) return false;
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
		{
			if (reference_count != 0) return true;
			unsigned int index = index_of(existential_intro_nodes, &proof);
			if (!visit_existential_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::EXISTENTIAL_INTRO_NODE, existential_intro_nodes[index])) return false;
			break;
		}

		case ProofType::DISJUNCTION_INTRODUCTION:
		case ProofType::DISJUNCTION_INTRODUCTION_LEFT:
		case ProofType::DISJUNCTION_INTRODUCTION_RIGHT:
		{
			if (reference_count != 0) return true;
			unsigned int index = index_of(disjunction_intro_nodes, &proof);
			if (!visit_disjunction_intro(std::forward<Visitor>(visitor)...)
			 || !changes.add(change_type::DISJUNCTION_INTRO_NODE, disjunction_intro_nodes[index])) return false;
			break;
		}

		case ProofType::CONJUNCTION_INTRODUCTION:
		{
			/* check if any of the conjuncts declare a fixed set size */
			for (unsigned int i = 0; i < proof.operand_array.length; i++) {
				Proof* operand = proof.operand_array[i];
				unsigned int constant; unsigned int set_size;
				if (is_set_size_proof(operand, constant, set_size))
					if (!requested_set_sizes.put(constant, set_size)) return false;
			}
		}

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
		case ProofType::INEQUALITY_INTRODUCTION:
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
		for (unsigned int i = operand_count; i > 0; i--) {
			if (operands[i - 1] == NULL) continue;
			if (!reference_counts.ensure_capacity(reference_counts.size + 1))
				return false;
			unsigned int index = reference_counts.index_of(operands[i - 1]);
			if (index == reference_counts.size) {
				reference_counts.keys[index] = operands[i - 1];
				reference_counts.values[index] = operands[i - 1]->reference_count;
				reference_counts.size++;
			}
			reference_counts.values[index]--;
			if (!get_theory_changes(*operands[i - 1], discharged_axioms, reference_counts, requested_set_sizes, changes, std::forward<Visitor>(visitor)...))
				return false;
		}
		discharged_axioms.length = old_discharged_axiom_count;
		requested_set_sizes.size = old_requested_set_size_count;
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
		array_map<unsigned int, unsigned int> requested_set_sizes(4);
		return get_theory_changes(proof, discharged_axioms, reference_counts, requested_set_sizes, changes, std::forward<Visitor>(visitor)...);
	}

	static bool unify_atom(Term* term, Term* atom,
			array<Formula*>& quantifiers,
			array_map<Formula*, Term*>& unifications)
	{
		if (term->type == TermType::VARIABLE) {
			if (atom->type != FormulaType::VARIABLE && atom->type != FormulaType::CONSTANT && atom->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

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
			const array<Formula*>& quantifiers,
			array_map<Formula*, Term*>& unifications,
			array_map<unsigned int, unsigned int>& variable_map)
	{
		if (first->type == FormulaType::VARIABLE) {
			if (second->type != FormulaType::VARIABLE && second->type != FormulaType::CONSTANT && second->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

			bool contains;
			unsigned int second_variable = variable_map.get(first->variable, contains);
			if (contains)
				return (second->type == FormulaType::VARIABLE && second->variable == second_variable);

			if (!unifications.ensure_capacity(unifications.size + 1))
				return false;
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
			return unify(first->binary.left, second->binary.left, quantifiers, unifications, variable_map)
				&& unify(first->binary.right, second->binary.right, quantifiers, unifications, variable_map);
		case FormulaType::BINARY_APPLICATION:
			return unify(first->ternary.first, second->ternary.first, quantifiers, unifications, variable_map)
				&& unify(first->ternary.second, second->ternary.second, quantifiers, unifications, variable_map)
				&& unify(first->ternary.third, second->ternary.third, quantifiers, unifications, variable_map);
		case FormulaType::NOT:
			return unify(first->unary.operand, second->unary.operand, quantifiers, unifications, variable_map);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			if (first->array.length != second->array.length) return false;
			for (unsigned int i = 0; i < first->array.length; i++)
				if (!unify(first->array.operands[i], second->array.operands[i], quantifiers, unifications, variable_map)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!variable_map.put(first->quantifier.variable, second->quantifier.variable)
			 || !unify(first->quantifier.operand, second->quantifier.operand, quantifiers, unifications, variable_map))
				return false;
			variable_map.size--;
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

	static inline bool unify(Formula* first, Formula* second,
			const array<Formula*>& quantifiers,
			array_map<Formula*, Term*>& unifications)
	{
		array_map<unsigned int, unsigned int> variable_map(4);
		return unify(first, second, quantifiers, unifications, variable_map);
	}

	static bool unify(Formula* first, Formula* second,
			const array<Formula*>& first_quantifiers,
			array_map<Formula*, Term*>& first_unifications,
			array_map<unsigned int, Term*>& second_unifications,
			array_map<unsigned int, unsigned int>& variable_map)
	{
		if (first->type == FormulaType::VARIABLE) {
			if (second->type != FormulaType::VARIABLE && second->type != FormulaType::CONSTANT && second->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

			bool contains;
			unsigned int second_variable = variable_map.get(first->variable, contains);
			if (contains)
				return (second->type == FormulaType::VARIABLE && second->variable == second_variable);

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
			if (first->type != FormulaType::VARIABLE && first->type != FormulaType::CONSTANT && first->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

			if (core::index_of(second->variable, variable_map.values, variable_map.size) < variable_map.size)
				return false;
			if (!second_unifications.ensure_capacity(second_unifications.size + 1))
				return false;
			unsigned int index = second_unifications.index_of(second->variable);
			if (index == second_unifications.size) {
				second_unifications.values[second_unifications.size] = first;
				second_unifications.keys[second_unifications.size++] = second->variable;
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
			return unify(first->binary.left, second->binary.left, first_quantifiers, first_unifications, second_unifications, variable_map)
				&& unify(first->binary.right, second->binary.right, first_quantifiers, first_unifications, second_unifications, variable_map);
		case FormulaType::BINARY_APPLICATION:
			return unify(first->ternary.first, second->ternary.first, first_quantifiers, first_unifications, second_unifications, variable_map)
				&& unify(first->ternary.second, second->ternary.second, first_quantifiers, first_unifications, second_unifications, variable_map)
				&& unify(first->ternary.third, second->ternary.third, first_quantifiers, first_unifications, second_unifications, variable_map);
		case FormulaType::NOT:
			return unify(first->unary.operand, second->unary.operand, first_quantifiers, first_unifications, second_unifications, variable_map);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			if (first->array.length != second->array.length) return false;
			for (unsigned int i = 0; i < first->array.length; i++)
				if (!unify(first->array.operands[i], second->array.operands[i], first_quantifiers, first_unifications, second_unifications, variable_map)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!variable_map.put(first->quantifier.variable, second->quantifier.variable)
			 || !unify(first->quantifier.operand, second->quantifier.operand, first_quantifiers, first_unifications, second_unifications, variable_map))
				return false;
			variable_map.size--;
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
			const array<Formula*>& first_quantifiers,
			array_map<Formula*, Term*>& first_unifications,
			array_map<unsigned int, Term*>& second_unifications)
	{
		array_map<unsigned int, unsigned int> variable_map(4);
		return unify(first, second, first_quantifiers, first_unifications, second_unifications, variable_map);
	}

	static bool unify(Formula* first, Formula* second,
			const array<Formula*>& first_quantifiers,
			array_map<Formula*, Term*>& first_unifications,
			const array<Formula*>& second_quantifiers,
			array_map<Formula*, Term*>& second_unifications,
			array_map<unsigned int, unsigned int>& variable_map)
	{
		if (first->type == FormulaType::VARIABLE) {
			if (second->type != FormulaType::VARIABLE && second->type != FormulaType::CONSTANT && second->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

			bool contains;
			unsigned int second_variable = variable_map.get(first->variable, contains);
			if (contains)
				return (second->type == FormulaType::VARIABLE && second->variable == second_variable);

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
			if (first->type != FormulaType::VARIABLE && first->type != FormulaType::CONSTANT && first->type != FormulaType::NUMBER)
				/* optimization: for now, we only consider unifications involving constants and numbers */
				return false;

			if (second->variable - 1 >= second_quantifiers.length)
				return false;
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
			return unify(first->binary.left, second->binary.left, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map)
				&& unify(first->binary.right, second->binary.right, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map);
		case FormulaType::BINARY_APPLICATION:
			return unify(first->ternary.first, second->ternary.first, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map)
				&& unify(first->ternary.second, second->ternary.second, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map)
				&& unify(first->ternary.third, second->ternary.third, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map);
		case FormulaType::NOT:
			return unify(first->unary.operand, second->unary.operand, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			if (first->array.length != second->array.length) return false;
			for (unsigned int i = 0; i < first->array.length; i++)
				if (!unify(first->array.operands[i], second->array.operands[i], first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			if (!variable_map.put(first->quantifier.variable, second->quantifier.variable)
			 || !unify(first->quantifier.operand, second->quantifier.operand, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map))
				return false;
			variable_map.size--;
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

	static inline bool unify(Formula* first, Formula* second,
			const array<Formula*>& first_quantifiers,
			array_map<Formula*, Term*>& first_unifications,
			const array<Formula*>& second_quantifiers,
			array_map<Formula*, Term*>& second_unifications)
	{
		array_map<unsigned int, unsigned int> variable_map(4);
		return unify(first, second, first_quantifiers, first_unifications, second_quantifiers, second_unifications, variable_map);
	}

	struct unification {
		Formula* subformula;
		array<Formula*> second_quantifiers;
		array_map<Formula*, Term*> first_unifications;
		array_map<Formula*, Term*> second_unifications;

		static inline void free(unification& u) {
			core::free(u.second_quantifiers);
			core::free(u.first_unifications);
			core::free(u.second_unifications);
		}
	};

	template<bool Polarity>
	static bool unify_subformula_helper(
			Formula* first, Formula* second,
			const array<Formula*>& first_quantifiers,
			array<Formula*>& second_quantifiers,
			array<unification>& unifications)
	{
		if (Polarity && unify(first, second,
				first_quantifiers, unifications[unifications.length].first_unifications,
				second_quantifiers, unifications[unifications.length].second_unifications))
		{
			if (!array_init(unifications[unifications.length].second_quantifiers, second_quantifiers.length)) {
				core::free(unifications[unifications.length].first_unifications);
				core::free(unifications[unifications.length].second_unifications);
				for (unification& u : unifications) core::free(u);
				return false;
			}
			for (unsigned int i = 0; i < second_quantifiers.length; i++) {
				array<Formula*>& arr = unifications[unifications.length].second_quantifiers;
				arr[arr.length++] = second_quantifiers[i];
			}
			unifications[unifications.length++].subformula = second;

			if (!unifications.ensure_capacity(unifications.length + 1)) {
				for (unification& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length].first_unifications, 4)) {
				for (unification& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length].second_unifications, 4)) {
				for (unification& u : unifications) core::free(u);
				core::free(unifications[unifications.length].first_unifications);
				return false;
			}
			return true;
		}

		switch (second->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
		case FormulaType::CONSTANT:
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
		case FormulaType::PARAMETER:
		case FormulaType::NUMBER:
		case FormulaType::STRING:
		case FormulaType::UINT_LIST:
		case FormulaType::UNARY_APPLICATION:
		case FormulaType::BINARY_APPLICATION:
			return true;
		case FormulaType::EQUALS:
			if (second->binary.right->type == hol_term_type::LAMBDA)
				return unify_subformula_helper<Polarity>(first, second->binary.right, first_quantifiers, second_quantifiers, unifications);
			return true;
		case FormulaType::IF_THEN:
			return unify_subformula_helper<!Polarity>(first, second->binary.left, first_quantifiers, second_quantifiers, unifications)
				&& unify_subformula_helper<Polarity>(first, second->binary.right, first_quantifiers, second_quantifiers, unifications);
		case FormulaType::NOT:
			return unify_subformula_helper<!Polarity>(first, second->unary.operand, first_quantifiers, second_quantifiers, unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			for (unsigned int i = 0; i < second->array.length; i++)
				if (!unify_subformula_helper<Polarity>(first, second->array.operands[i], first_quantifiers, second_quantifiers, unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
		{
			if (!second_quantifiers.add(second)) {
				for (unification& u : unifications) core::free(u);
				core::free(unifications[unifications.length].first_unifications);
				core::free(unifications[unifications.length].second_unifications);
				return false;
			}
			bool result = unify_subformula_helper<Polarity>(first, second->quantifier.operand, first_quantifiers, second_quantifiers, unifications);
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
		fprintf(stderr, "theory.unify_subformula_helper ERROR: Unrecognized FormulaType.\n");
		for (unification& u : unifications) core::free(u);
		core::free(unifications[unifications.length].first_unifications);
		core::free(unifications[unifications.length].second_unifications);
		return false;
	}

	template<bool Polarity>
	static inline bool unify_subformula(
			Formula* first, Formula* second,
			const array<Formula*>& first_quantifiers,
			array<Formula*>& second_quantifiers,
			array<unification>& unifications)
	{
		if (!array_map_init(unifications[unifications.length].first_unifications, 4)) {
			return false;
		} else if (!array_map_init(unifications[unifications.length].second_unifications, 4)) {
			core::free(unifications[unifications.length].first_unifications);
			return false;
		} else if (!unify_subformula_helper<Polarity>(first, second, first_quantifiers, second_quantifiers, unifications)) {
			return false;
		}
		core::free(unifications[unifications.length].first_unifications);
		core::free(unifications[unifications.length].second_unifications);
		return true;
	}

	template<bool Polarity>
	static bool unify_subformula_helper(
			Formula* first, Formula* second,
			array<Formula*>& quantifiers,
			array<array_map<Formula*, Term*>>& unifications)
	{
		if (Polarity && unify(second, first, quantifiers, unifications[unifications.length])) {
			unifications.length++;

			if (!unifications.ensure_capacity(unifications.length + 1)) {
				for (auto& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length], 4)) {
				for (auto& u : unifications) core::free(u);
				return false;
			}
			return true;
		}

		if (!Polarity && first->type == FormulaType::NOT && second->type != FormulaType::NOT
		 && unify(second, first->unary.operand, quantifiers, unifications[unifications.length]))
		{
			unifications.length++;

			if (!unifications.ensure_capacity(unifications.length + 1)) {
				for (auto& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length], 4)) {
				for (auto& u : unifications) core::free(u);
				return false;
			}
			return true;
		}

		switch (second->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
		case FormulaType::CONSTANT:
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
		case FormulaType::PARAMETER:
		case FormulaType::NUMBER:
		case FormulaType::STRING:
		case FormulaType::UINT_LIST:
		case FormulaType::UNARY_APPLICATION:
		case FormulaType::BINARY_APPLICATION:
			return true;
		case FormulaType::EQUALS:
			if (second->binary.right->type == hol_term_type::LAMBDA)
				return unify_subformula_helper<Polarity>(first, second->binary.right, quantifiers, unifications);
			else return true;
		case FormulaType::IF_THEN:
			return unify_subformula_helper<!Polarity>(first, second->binary.left, quantifiers, unifications)
				&& unify_subformula_helper<Polarity>(first, second->binary.right, quantifiers, unifications);
		case FormulaType::NOT:
			return unify_subformula_helper<!Polarity>(first, second->unary.operand, quantifiers, unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			for (unsigned int i = 0; i < second->array.length; i++)
				if (!unify_subformula_helper<Polarity>(first, second->array.operands[i], quantifiers, unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
		{
			if (!quantifiers.add(second)) {
				for (auto& u : unifications) core::free(u);
				core::free(unifications[unifications.length]);
				return false;
			}
			bool result = unify_subformula_helper<Polarity>(first, second->quantifier.operand, quantifiers, unifications);
			quantifiers.length--;
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
		fprintf(stderr, "theory.unify_subformula_helper ERROR: Unrecognized FormulaType.\n");
		for (auto& u : unifications) core::free(u);
		core::free(unifications[unifications.length]);
		return false;
	}

	template<bool Polarity>
	static inline bool unify_subformula(
			Formula* first, Formula* second,
			array<Formula*>& quantifiers,
			array<array_map<Formula*, Term*>>& unifications)
	{
		if (!array_map_init(unifications[unifications.length], 4)) {
			return false;
		} else if (!unify_subformula_helper<Polarity>(first, second, quantifiers, unifications)) {
			return false;
		}
		core::free(unifications[unifications.length]);
		return true;
	}

	static bool unify_subformula_helper(
			Formula* first, Formula* second,
			array<Formula*>& quantifiers,
			array<array_map<Formula*, Term*>>& unifications)
	{
		if (unify(second, first, quantifiers, unifications[unifications.length])) {
			unifications.length++;

			if (!unifications.ensure_capacity(unifications.length + 1)) {
				for (auto& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length], 4)) {
				for (auto& u : unifications) core::free(u);
				return false;
			}
			return true;
		}

		if (first->type == FormulaType::NOT && second->type != FormulaType::NOT
		 && unify(second, first->unary.operand, quantifiers, unifications[unifications.length]))
		{
			unifications.length++;

			if (!unifications.ensure_capacity(unifications.length + 1)) {
				for (auto& u : unifications) core::free(u);
				return false;
			} else if (!array_map_init(unifications[unifications.length], 4)) {
				for (auto& u : unifications) core::free(u);
				return false;
			}
			return true;
		}

		switch (second->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
		case FormulaType::CONSTANT:
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
		case FormulaType::PARAMETER:
		case FormulaType::NUMBER:
		case FormulaType::STRING:
		case FormulaType::UINT_LIST:
		case FormulaType::UNARY_APPLICATION:
		case FormulaType::BINARY_APPLICATION:
			return true;
		case FormulaType::EQUALS:
			if (second->binary.right->type == hol_term_type::LAMBDA)
				return unify_subformula_helper(first, second->binary.right, quantifiers, unifications);
			else return true;
		case FormulaType::IF_THEN:
			return unify_subformula_helper(first, second->binary.left, quantifiers, unifications)
				&& unify_subformula_helper(first, second->binary.right, quantifiers, unifications);
		case FormulaType::NOT:
			return unify_subformula_helper(first, second->unary.operand, quantifiers, unifications);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			for (unsigned int i = 0; i < second->array.length; i++)
				if (!unify_subformula_helper(first, second->array.operands[i], quantifiers, unifications)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
		{
			if (!quantifiers.add(second)) {
				for (auto& u : unifications) core::free(u);
				core::free(unifications[unifications.length]);
				return false;
			}
			bool result = unify_subformula_helper(first, second->quantifier.operand, quantifiers, unifications);
			quantifiers.length--;
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
		fprintf(stderr, "theory.unify_subformula_helper ERROR: Unrecognized FormulaType.\n");
		for (auto& u : unifications) core::free(u);
		core::free(unifications[unifications.length]);
		return false;
	}

	static inline bool unify_subformula(
			Formula* first, Formula* second,
			array<Formula*>& quantifiers,
			array<array_map<Formula*, Term*>>& unifications)
	{
		if (!array_map_init(unifications[unifications.length], 4)) {
			return false;
		} else if (!unify_subformula_helper(first, second, quantifiers, unifications)) {
			return false;
		}
		core::free(unifications[unifications.length]);
		return true;
	}

	struct hypothetical_reasoner {
		struct axiom {
			Formula* formula;
			bool negated;
			instantiation_tuple tuple;

			static inline void move(const axiom& src, axiom& dst) {
				dst.formula = src.formula;
				dst.negated = src.negated;
				core::move(src.tuple, dst.tuple);
			}

			static inline void free(axiom& ax) {
				core::free(ax.tuple);
			}
		};

		//array<unsigned int> indices;
		array<axiom> axioms;
		array<unsigned int> set_ids;

		hypothetical_reasoner() : /*indices(4),*/ axioms(4), set_ids(16) { }

		~hypothetical_reasoner() {
			for (axiom& ax : axioms)
				core::free(ax);
		}

		template<bool Negated>
		inline bool push_axiom(
				unsigned int set_id, Formula* formula,
				array<variable_assignment>& possible_values,
				unsigned int& pushed_axiom_count,
				const array<unification>& unifications)
		{
			if (!axioms.ensure_capacity(axioms.length + possible_values.length))
				return false;
			pushed_axiom_count = 0;
			for (variable_assignment& assignment : possible_values) {
				if (assignment.matching_axiom != nullptr) continue;

				/* make sure there is a valid unification */
				bool found_valid_unification = false;
				for (const unification& u : unifications) {
					bool valid = true;
					for (const auto& entry : u.first_unifications) {
						const instantiation& key = assignment.assignment.values[entry.key->quantifier.variable - 1];
						if (entry.value->type == TermType::CONSTANT) {
							instantiation temp(instantiation_type::CONSTANT);
							temp.constant = entry.value->constant;
							if (!is_subset(temp, key)) {
								valid = false;
								break;
							}
						} else if (entry.value->type == TermType::NUMBER) {
							instantiation temp(instantiation_type::NUMBER);
							temp.number = entry.value->number;
							if (!is_subset(temp, key)) {
								valid = false;
								break;
							}
						} else if (entry.value->type == TermType::STRING) {
							instantiation temp(instantiation_type::STRING);
							temp.str = entry.value->str;
							if (!is_subset(temp, key)) {
								valid = false;
								break;
							}
						}
					}
					if (valid) {
						found_valid_unification = true;
						break;
					}
				}
				if (!found_valid_unification) continue;

				/* only push an axiom if `formula` is not already an axiom */
				bool already_has_axiom = false;
				for (const axiom& other_axiom : axioms) {
					if ((other_axiom.formula == formula || *other_axiom.formula == *formula) && other_axiom.tuple == assignment.assignment) {
						already_has_axiom = true;
						assignment.matching_axiom = formula;
						break;
					}
				}
				if (already_has_axiom) continue;

				axioms[axioms.length].formula = formula;
				axioms[axioms.length].negated = Negated;
				if (!::init(axioms[axioms.length].tuple, assignment.assignment))
					return false;
				axioms.length++;
				pushed_axiom_count++;
			}
			return set_ids.add(set_id);

			// if (!indices.ensure_capacity(indices.length + 1))
			// 	return false;
			// array<axiom> new_axioms(possible_values.length);
			// for (instantiation_tuple& tuple : possible_values) {
			// 	new_axioms[new_axioms.length].formula = formula;
			// 	new_axioms[new_axioms.length].negated = Negated;
			// 	if (!init(new_axioms[new_axioms.length].tuple, tuple)) {
			// 		for (axiom& ax : new_axioms)
			// 			core::free(ax);
			// 		return false;
			// 	}
			// 	new_axioms.length++;
			// }

			// for (unsigned int i = 0; i < old_length; i++) {
			// 	array_map<Formula*, Term*> unifications(4);
			// 	array_map<unsigned int, Term*> second_unifications(4);
			// 	if (Negated != axioms[i].negated)
			// 		continue;
			// 	if (!unify(formula, axioms[i].formula, quantifiers, unifications, second_unifications))
			// 		continue;

			// 	if (unifications.size == 0) {
			// 		new_axioms_pushed = 0;
			// 		return true;
			// 	}

			// 	array<axiom> temp_axioms((unifications.size + axioms[i].tuple.equal_indices.size + axioms[i].tuple.not_equal_indices.size + axioms[i].tuple.ge_indices.size) * new_axioms.length);
			// 	for (axiom& new_axiom : new_axioms) {
			// 		/* subtract `new_axiom` with `axioms[i]` and add the result
			// 		   to `temp_axioms` to ensure that `axioms` are disjoint */
			// 		instantiation_tuple intersection(new_axiom.tuple);
			// 		bool has_intersection = true;
			// 		for (unsigned int j = 0; j < unifications.size; j++) {
			// 			if (j > 0) {
			// 				if (unifications.values[j]->type == TermType::VARIABLE) {
			// 					instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
			// 					if (!intersect(dummy, intersection.values[unifications.keys[j]->quantifier.variable - 1], axioms[i].tuple.values[unifications.values[j]->variable - 1])) {
			// 						has_intersection = false;
			// 						break;
			// 					} else if (!intersection.change_value(unifications.keys[j]->quantifier.variable - 1, dummy)) {
			// 						has_intersection = false;
			// 						core::free(dummy);
			// 						break;
			// 					}
			// 					core::free(dummy);
			// 				} else {
			// 					if (!intersection.unify_value(unifications.keys[j]->quantifier.variable - 1, unifications.values[j])) {
			// 						has_intersection = false;
			// 						break;
			// 					}
			// 				}
			// 			}

			// 			instantiation_tuple& new_tuple = temp_axioms[temp_axioms.length].tuple;
			// 			temp_axioms[temp_axioms.length].formula = new_axiom.formula;
			// 			temp_axioms[temp_axioms.length].negated = new_axiom.negated;
			// 			if (!init(new_tuple, intersection)) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			if (unifications.values[j]->type == TermType::VARIABLE) {
			// 				instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
			// 				if (!subtract(dummy, new_tuple.values[unifications.keys[j]->quantifier.variable - 1], axioms[i].tuple.values[unifications.values[j]->variable - 1])) {
			// 					core::free(new_tuple);
			// 					continue;
			// 				} else if (!new_tuple.change_value(unifications.keys[j]->quantifier.variable - 1, dummy)) {
			// 					core::free(dummy);
			// 					core::free(new_tuple);
			// 					continue;
			// 				}
			// 				core::free(dummy);
			// 			} else {
			// 				if (!new_tuple.antiunify_value(unifications.keys[j]->quantifier.variable - 1, unifications.values[j])) {
			// 					core::free(new_tuple);
			// 					continue;
			// 				}
			// 			}
			// 			temp_axioms.length++;
			// 		}
			// 		if (!has_intersection) continue;

			// 		if (unifications.values[unifications.size - 1]->type == TermType::VARIABLE) {
			// 			instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
			// 			if (!intersect(dummy, intersection.values[unifications.keys[unifications.size - 1]->quantifier.variable - 1], axioms[i].tuple.values[unifications.values[unifications.size - 1]->variable - 1])) {
			// 				continue;
			// 			} else if (!intersection.change_value(unifications.keys[unifications.size - 1]->quantifier.variable - 1, dummy)) {
			// 				core::free(dummy);
			// 				continue;
			// 			}
			// 			core::free(dummy);
			// 		} else {
			// 			if (!intersection.unify_value(unifications.keys[unifications.size - 1]->quantifier.variable - 1, unifications.values[unifications.size - 1]))
			// 				continue;
			// 		}

			// 		/* map the inter-variable constraints from `axioms[i].tuple` into the space of `new_axiom.tuple` */
			// 		array_map<uint_fast8_t, uint_fast8_t> variable_map(4);
			// 		for (const auto& unification : unifications) {
			// 			if (unification.value->type == TermType::VARIABLE
			// 			 && !variable_map.put(unification.value->variable - 1, unification.key->quantifier.variable - 1))
			// 			{
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 		}
			// 		array<pair<uint_fast8_t, uint_fast8_t>> mapped_equal_indices(max(1, axioms[i].tuple.equal_indices.size));
			// 		array<pair<uint_fast8_t, uint_fast8_t>> mapped_not_equal_indices(max(1, axioms[i].tuple.not_equal_indices.size));
			// 		array<pair<uint_fast8_t, uint_fast8_t>> mapped_ge_indices(max(1, axioms[i].tuple.ge_indices.size));
			// 		for (const pair<uint_fast8_t, uint_fast8_t>& entry : axioms[i].tuple.equal_indices) {
			// 			bool contains;
			// 			uint_fast8_t mapped_key = variable_map.get(entry.key, contains);
			// 			if (!contains) continue;
			// 			uint_fast8_t mapped_value = variable_map.get(entry.value, contains);
			// 			if (!contains) continue;
			// 			mapped_equal_indices[mapped_equal_indices.length++] = {mapped_key, mapped_value};
			// 		} for (const pair<uint_fast8_t, uint_fast8_t>& entry : axioms[i].tuple.not_equal_indices) {
			// 			bool contains;
			// 			uint_fast8_t mapped_key = variable_map.get(entry.key, contains);
			// 			if (!contains) continue;
			// 			uint_fast8_t mapped_value = variable_map.get(entry.value, contains);
			// 			if (!contains) continue;
			// 			mapped_not_equal_indices[mapped_not_equal_indices.length++] = {mapped_key, mapped_value};
			// 		} for (const pair<uint_fast8_t, uint_fast8_t>& entry : axioms[i].tuple.ge_indices) {
			// 			bool contains;
			// 			uint_fast8_t mapped_key = variable_map.get(entry.key, contains);
			// 			if (!contains) continue;
			// 			uint_fast8_t mapped_value = variable_map.get(entry.value, contains);
			// 			if (!contains) continue;
			// 			mapped_ge_indices[mapped_ge_indices.length++] = {mapped_key, mapped_value};
			// 		}

			// 		for (unsigned int j = 0; j < mapped_equal_indices.length; j++) {
			// 			if (j > 0) {
			// 				Term* var = Term::new_variable(mapped_equal_indices[j].value + 1);
			// 				if (var == nullptr) {
			// 					for (axiom& entry : new_axioms) core::free(entry);
			// 					for (axiom& entry : temp_axioms) core::free(entry);
			// 					return false;
			// 				}
			// 				if (!intersection.unify_value(mapped_equal_indices[j].key, var)) {
			// 					core::free(*var); core::free(var);
			// 					has_intersection = false;
			// 					break;
			// 				}
			// 				core::free(*var); core::free(var);
			// 			}

			// 			instantiation_tuple& new_tuple = temp_axioms[temp_axioms.length].tuple;
			// 			temp_axioms[temp_axioms.length].formula = new_axiom.formula;
			// 			temp_axioms[temp_axioms.length].negated = new_axiom.negated;
			// 			if (!init(new_tuple, intersection)) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			Term* var = Term::new_variable(mapped_equal_indices[j].value + 1);
			// 			if (!new_tuple.antiunify_value(mapped_equal_indices[j].key, var)) {
			// 				core::free(*var); core::free(var);
			// 				core::free(new_tuple);
			// 				continue;
			// 			}
			// 			core::free(*var); core::free(var);
			// 			temp_axioms.length++;
			// 		}
			// 		if (!has_intersection) continue;

			// 		if (mapped_equal_indices.length != 0) {
			// 			Term* var = Term::new_variable(mapped_equal_indices[mapped_equal_indices.length - 1].value + 1);
			// 			if (var == nullptr) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			if (!intersection.unify_value(mapped_equal_indices[mapped_equal_indices.length - 1].key, var)) {
			// 				core::free(*var); core::free(var);
			// 				continue;
			// 			}
			// 			core::free(*var); core::free(var);
			// 		}

			// 		for (unsigned int j = 0; j < mapped_not_equal_indices.length; j++) {
			// 			if (j > 0) {
			// 				Term* var = Term::new_variable(mapped_not_equal_indices[j].value + 1);
			// 				if (var == nullptr) {
			// 					for (axiom& entry : new_axioms) core::free(entry);
			// 					for (axiom& entry : temp_axioms) core::free(entry);
			// 					return false;
			// 				}
			// 				if (!intersection.antiunify_value(mapped_not_equal_indices[j].key, var)) {
			// 					core::free(*var); core::free(var);
			// 					has_intersection = false;
			// 					break;
			// 				}
			// 				core::free(*var); core::free(var);
			// 			}

			// 			instantiation_tuple& new_tuple = temp_axioms[temp_axioms.length].tuple;
			// 			temp_axioms[temp_axioms.length].formula = new_axiom.formula;
			// 			temp_axioms[temp_axioms.length].negated = new_axiom.negated;
			// 			if (!init(new_tuple, intersection)) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			Term* var = Term::new_variable(mapped_not_equal_indices[j].value + 1);
			// 			if (!new_tuple.unify_value(mapped_not_equal_indices[j].key, var)) {
			// 				core::free(*var); core::free(var);
			// 				core::free(new_tuple);
			// 				continue;
			// 			}
			// 			core::free(*var); core::free(var);
			// 			temp_axioms.length++;
			// 		}
			// 		if (!has_intersection) continue;

			// 		if (mapped_not_equal_indices.length != 0) {
			// 			Term* var = Term::new_variable(mapped_not_equal_indices[mapped_not_equal_indices.length - 1].value + 1);
			// 			if (var == nullptr) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			if (!intersection.antiunify_value(mapped_not_equal_indices[mapped_not_equal_indices.length - 1].key, var)) {
			// 				core::free(*var); core::free(var);
			// 				continue;
			// 			}
			// 			core::free(*var); core::free(var);
			// 		}

			// 		for (unsigned int j = 0; j < mapped_ge_indices.length; j++) {
			// 			if (j > 0) {
			// 				Term* var = Term::new_variable(mapped_ge_indices[j].value + 1);
			// 				if (var == nullptr) {
			// 					for (axiom& entry : new_axioms) core::free(entry);
			// 					for (axiom& entry : temp_axioms) core::free(entry);
			// 					return false;
			// 				}
			// 				if (!intersection.unify_greater_than_or_equal(mapped_ge_indices[j].key, var)) {
			// 					core::free(*var); core::free(var);
			// 					has_intersection = false;
			// 					break;
			// 				}
			// 				core::free(*var); core::free(var);
			// 			}

			// 			instantiation_tuple& new_tuple = temp_axioms[temp_axioms.length].tuple;
			// 			temp_axioms[temp_axioms.length].formula = new_axiom.formula;
			// 			temp_axioms[temp_axioms.length].negated = new_axiom.negated;
			// 			if (!init(new_tuple, intersection)) {
			// 				for (axiom& entry : new_axioms) core::free(entry);
			// 				for (axiom& entry : temp_axioms) core::free(entry);
			// 				return false;
			// 			}
			// 			Term* var = Term::new_variable(mapped_ge_indices[j].value + 1);
			// 			if (!new_tuple.unify_less_than_or_equal(mapped_ge_indices[j].key, var)
			// 			 || !new_tuple.antiunify_value(mapped_ge_indices[j].key, var))
			// 			{
			// 				core::free(*var); core::free(var);
			// 				core::free(new_tuple);
			// 				continue;
			// 			}
			// 			core::free(*var); core::free(var);
			// 			temp_axioms.length++;
			// 		}
			// 	}
			// 	for (axiom& entry : new_axioms)
			// 		core::free(entry);
			// 	core::swap(temp_axioms, new_axioms);
			// 	if (new_axioms.length == 0) {
			// 		new_axioms_pushed = 0;
			// 		return true;
			// 	}
			//}

			//if (!axioms.ensure_capacity(axioms.length + new_axioms.length)) {
			//	for (axiom& entry : new_axioms)
			//		core::free(entry);
			//	return false;
			//}
			//indices[indices.length++] = axioms.length;
			//for (axiom& new_axiom : new_axioms)
			//	core::move(new_axiom, axioms[axioms.length++]);
			//new_axioms_pushed = new_axioms.length;
			//return true;
		}

		inline void pop_axiom(unsigned int count) {
			for (unsigned int i = axioms.length - count; i < axioms.length; i++)
				core::free(axioms[i]);
			axioms.length -= count;
			set_ids.length--;
		}

		inline bool is_set_on_stack(unsigned int set_id) {
			return set_ids.contains(set_id);
		}

		inline bool get_axiom_instantiation(
				variable_assignment& assignment,
				const variable_assignment& src,
				unsigned int axiom_count)
		{
			for (unsigned int i = 0; i < src.axiom_assignments.size; i++) {
				if (src.axiom_assignments.keys[i] >= axioms.length - axiom_count) {
					if (!::init(assignment, src.axiom_assignments.values[i].assignment))
						return false;
					if (!assignment.axiom_assignments.ensure_capacity(src.axiom_assignments.size - 1)) {
						core::free(assignment);
						return false;
					}
					for (unsigned int j = 0; j < src.axiom_assignments.size; j++) {
						if (j == i) continue;
						assignment.axiom_assignments.keys[assignment.axiom_assignments.size] = src.axiom_assignments.keys[j];
						if (!::init(assignment.axiom_assignments.values[assignment.axiom_assignments.size], src.axiom_assignments.values[j])) {
							core::free(assignment);
							return false;
						}
						assignment.axiom_assignments.size++;
					}
					return true;
				}
			}
			return false;
		}

		template<bool Contradiction>
		bool is_provable_by_axiom_without_abduction(
				Formula* formula, array<Formula*>& quantifiers,
				const array<variable_assignment>& possible_values,
				array<variable_assignment>& new_possible_values) const
		{
			for (unsigned int i = 0; i < axioms.length; i++) {
				const axiom& ax = axioms[i];
				if (Contradiction != ax.negated) continue;

				array_map<Formula*, Term*> unifications(4);
				array_map<unsigned int, Term*> second_unifications(4);
				if (!unify(formula, ax.formula, quantifiers, unifications, second_unifications))
					continue;

				for (const variable_assignment& assignment : possible_values) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + 1))
						return false;

					/* map the variable constraints from the space of `ax` to that of each `(formula, possible_values)` */
					variable_assignment& new_assignment = new_possible_values[new_possible_values.length];
					if (!::init(new_assignment, assignment)) {
						return false;
					} else if (!new_assignment.axiom_assignments.ensure_capacity(new_assignment.axiom_assignments.size + 1)) {
						core::free(new_assignment);
						return false;
					}
					bool is_new = false;
					unsigned int index = new_assignment.axiom_assignments.index_of(i);
					if (index == new_assignment.axiom_assignments.size) {
						new_assignment.axiom_assignments.keys[index] = i;
						if (!::init(new_assignment.axiom_assignments.values[index], ax.tuple)) {
							core::free(new_assignment);
							return false;
						}
						new_assignment.axiom_assignments.size++;
						is_new = true;
					}
					axiom_assignment& new_axiom_assignment = new_assignment.axiom_assignments.values[index];

					bool unifies = true;
					for (const auto& unification : second_unifications) {
						if (!new_axiom_assignment.assignment.unify_value(unification.key - 1, unification.value)) {
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						core::free(new_assignment);
						continue;
					}

					for (const auto& unification : unifications) {
						if (unification.value->type == TermType::VARIABLE) {
							if (!new_axiom_assignment.variable_map.ensure_capacity(new_axiom_assignment.variable_map.size + 1)) {
								core::free(new_assignment);
								return false;
							}
							unsigned int first_other_variable = 0;
							unsigned int second_other_variable = 0;
							bool contains = false;
							for (unsigned int j = 0; j < new_axiom_assignment.variable_map.size; j++) {
								if (new_axiom_assignment.variable_map.keys[j] == unification.key->quantifier.variable) {
									if (new_axiom_assignment.variable_map.values[j] == unification.value->variable) {
										contains = true;
										break;
									} else {
										second_other_variable = new_axiom_assignment.variable_map.values[j];
									}
								} else if (new_axiom_assignment.variable_map.values[j] == unification.value->variable) {
									first_other_variable = new_axiom_assignment.variable_map.keys[j];
								}
							}
							if (contains) continue;

							new_axiom_assignment.variable_map.keys[new_axiom_assignment.variable_map.size] = unification.key->quantifier.variable;
							new_axiom_assignment.variable_map.values[new_axiom_assignment.variable_map.size++] = unification.value->variable;
							if (first_other_variable != 0 && !new_assignment.assignment.unify_indices(unification.key->quantifier.variable - 1, first_other_variable - 1)) {
								unifies = false;
								break;
							} if (second_other_variable != 0 && !new_axiom_assignment.assignment.unify_indices(unification.value->variable - 1, second_other_variable - 1)) {
								unifies = false;
								break;
							}
						} else {
							if (!new_assignment.assignment.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
								unifies = false;
								break;
							}
						}
					}
					if (!unifies) {
						core::free(new_assignment);
						continue;
					}

					if (is_new) {
						for (const auto& entry : new_axiom_assignment.variable_map) {
		 					instantiation& dummy = *((instantiation*) alloca(sizeof(instantiation)));
		 					if (!intersect(dummy, new_assignment.assignment.values[entry.key - 1], new_axiom_assignment.assignment.values[entry.value - 1])) {
		 						unifies = false;
		 						break;
		 					} else if (!new_assignment.assignment.change_value(entry.key - 1, dummy)
									|| !new_axiom_assignment.assignment.change_value(entry.value - 1, dummy))
							{
		 						unifies = false;
		 						core::free(dummy);
		 						break;
		 					}
		 					core::free(dummy);
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_axiom_assignment.assignment.equal_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.keys[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.keys[second_index] - 1;
							if (mapped_second < mapped_first) swap(mapped_first, mapped_second);
							if (!new_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_assignment.assignment.unify_indices(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_axiom_assignment.assignment.not_equal_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.keys[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.keys[second_index] - 1;
							if (mapped_second < mapped_first) swap(mapped_first, mapped_second);
							if (!new_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_assignment.assignment.antiunify_indices(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_axiom_assignment.assignment.ge_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.values, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.keys[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.keys[second_index] - 1;
							if (!new_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_assignment.assignment.unify_greater_than_or_equal(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_assignment.assignment.equal_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.values[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.values[second_index] - 1;
							if (mapped_second < mapped_first) swap(mapped_first, mapped_second);
							if (!new_axiom_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_axiom_assignment.assignment.unify_indices(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_assignment.assignment.not_equal_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.values[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.values[second_index] - 1;
							if (mapped_second < mapped_first) swap(mapped_first, mapped_second);
							if (!new_axiom_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_axiom_assignment.assignment.antiunify_indices(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}

						for (const pair<uint_fast8_t, uint_fast8_t>& entry : new_assignment.assignment.ge_indices) {
							unsigned int first_index = core::index_of(entry.key, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (first_index == new_axiom_assignment.variable_map.size) continue;
							unsigned int second_index = core::index_of(entry.value, new_axiom_assignment.variable_map.keys, new_axiom_assignment.variable_map.size);
							if (second_index == new_axiom_assignment.variable_map.size) continue;
							uint_fast8_t mapped_first = new_axiom_assignment.variable_map.values[first_index] - 1;
							uint_fast8_t mapped_second = new_axiom_assignment.variable_map.values[second_index] - 1;
							if (!new_axiom_assignment.assignment.equal_indices.contains(make_pair(mapped_first, mapped_second))
							 && !new_axiom_assignment.assignment.unify_greater_than_or_equal(mapped_first, mapped_second))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_assignment);
							continue;
						}
					}
					new_possible_values.length++;
				}
			}
			return true;
		}
	};

	struct default_prover {
		typedef const array<tuple>& ProvableElementArray;

		const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets;
		const array<pair<Proof*, bool>>& implication_axioms;
		hypothetical_reasoner h;
		unsigned int removed_set;

		default_prover(
				const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets,
				const array<pair<Proof*, bool>>& implication_axioms,
				unsigned int removed_set = 0) : sets(sets), implication_axioms(implication_axioms), removed_set(removed_set) { }

		inline const array<tuple>& get_provable_elements(unsigned int set_id) const {
			return sets.sets[set_id].provable_elements;
		}

		inline unsigned int get_provable_element_count(unsigned int set_id) const {
			return sets.sets[set_id].provable_elements.length;
		}

		inline void free_provable_elements(const array<tuple>& elements) const { }

		inline bool is_antecedent_satisfied(unsigned int implication_id) const {
			return implication_axioms[implication_id].value;
		}

		inline bool is_set_removed(unsigned int set) const {
			return (set == removed_set);
		}

		template<bool Negated>
		inline bool push_axiom(
				unsigned int set_id, Formula* formula,
				array<variable_assignment>& possible_values,
				unsigned int& pushed_axiom_count,
				const array<unification>& unifications)
		{
			return h.template push_axiom<Negated>(set_id, formula, possible_values, pushed_axiom_count, unifications);
		}

		inline void pop_axiom(unsigned int count) {
			h.pop_axiom(count);
		}

		inline bool is_set_on_stack(unsigned int set_id) {
			return h.is_set_on_stack(set_id);
		}

		inline bool get_axiom_instantiation(
				variable_assignment& assignment,
				const variable_assignment& src,
				unsigned int axiom_count)
		{
			return h.get_axiom_instantiation(assignment, src, axiom_count);
		}

		template<bool Contradiction>
		inline bool is_provable_by_axiom_without_abduction(
				Formula* formula, array<Formula*>& quantifiers,
				const array<variable_assignment>& possible_values,
				array<variable_assignment>& new_possible_values) const
		{
			return h.template is_provable_by_axiom_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values);
		}
	};

	template<typename T>
	struct array_with_extra_elements {
		struct iterator {
			const array_with_extra_elements<T>& a;
			unsigned int index;

			iterator(const array_with_extra_elements<T>& a, unsigned int index) : a(a), index(index) { }

			inline void operator ++ () {
				index++;
			}

			inline const T& operator * () {
				if (index < a.data_length)
					return a.data[index];
				else return a.extra_data[index - a.data_length];
			}

			inline bool operator != (const iterator& other) {
				return index != other.index;
			}
		};

		T* data;
		unsigned int data_length;
		T* extra_data;

		unsigned int length;

		array_with_extra_elements(T* data, unsigned int data_length, T* extra_data, unsigned int extra_data_length) :
				data(data), data_length(data_length), extra_data(extra_data), length(data_length + extra_data_length)
		{ }

		iterator begin() const { return iterator(*this, 0); }
		iterator end() const { return iterator(*this, length); }

		inline bool contains(const tuple& element) const
		{
			if (core::index_of(element, data, data_length) < data_length)
				return true;
			unsigned int extra_count = length - data_length;
			return (core::index_of(element, extra_data, extra_count) < extra_count);
		}
	};

	struct set_membership_prover {
		typedef array_with_extra_elements<tuple> ProvableElementArray;

		const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets;
		const array<pair<Proof*, bool>>& implication_axioms;
		const array_map<tuple, array<unsigned int>>& new_elements;
		const array<unsigned int>& new_antecedents;
		hypothetical_reasoner h;

		set_membership_prover(
				const set_reasoning<built_in_predicates, ProofCalculus, Canonicalizer>& sets,
				const array<pair<Proof*, bool>>& implication_axioms,
				const array_map<tuple, array<unsigned int>>& new_elements,
				const array<unsigned int>& new_antecedents) : sets(sets), implication_axioms(implication_axioms), new_elements(new_elements), new_antecedents(new_antecedents) { }

		inline array_with_extra_elements<tuple> get_provable_elements(unsigned int set_id) const
		{
			array<tuple>& extra_data = *((array<tuple>*) alloca(sizeof(array<tuple>)));
			array_init(extra_data, 4);
			for (const auto& entry : new_elements) {
				for (unsigned int other_set : entry.value) {
					if (sets.sets[set_id].descendants.contains(other_set)) {
						extra_data.ensure_capacity(extra_data.length + 1);

						if (!sets.sets[set_id].provable_elements.contains(entry.key)) {
							::init(extra_data[extra_data.length], entry.key);
							extra_data.length++;
						}
					}
				}
			}
			return array_with_extra_elements<tuple>(sets.sets[set_id].provable_elements.data, sets.sets[set_id].provable_elements.length, extra_data.data, extra_data.length);
		}

		inline unsigned int get_provable_element_count(unsigned int set_id) const
		{
			unsigned int count = sets.sets[set_id].provable_elements.length;
			for (const auto& entry : new_elements) {
				for (unsigned int other_set : entry.value) {
					if (sets.sets[set_id].descendants.contains(other_set)) {
						if (!sets.sets[set_id].provable_elements.contains(entry.key)) {
							count++;
						}
					}
				}
			}
			return count;
		}

		inline void free_provable_elements(array_with_extra_elements<tuple>& elements) const {
			unsigned int count = elements.length - elements.data_length;
			for (unsigned int i = 0; i < count; i++)
				core::free(elements.extra_data[i]);
			core::free(elements.extra_data);
		}

		inline bool is_antecedent_satisfied(unsigned int implication_id) const {
			return implication_axioms[implication_id].value || new_antecedents.contains(implication_id);
		}

		inline constexpr bool is_set_removed(unsigned int set) const { return false; }

		template<bool Negated>
		inline bool push_axiom(
				unsigned int set_id, Formula* formula,
				array<variable_assignment>& possible_values,
				unsigned int& pushed_axiom_count,
				const array<unification>& unifications)
		{
			return h.template push_axiom<Negated>(set_id, formula, possible_values, pushed_axiom_count, unifications);
		}

		inline void pop_axiom(unsigned int count) {
			h.pop_axiom(count);
		}

		inline bool is_set_on_stack(unsigned int set_id) {
			return h.is_set_on_stack(set_id);
		}

		inline bool get_axiom_instantiation(
				variable_assignment& assignment,
				const variable_assignment& src,
				unsigned int axiom_count)
		{
			return h.get_axiom_instantiation(assignment, src, axiom_count);
		}

		template<bool Contradiction>
		inline bool is_provable_by_axiom_without_abduction(
				Formula* formula, array<Formula*>& quantifiers,
				const array<variable_assignment>& possible_values,
				array<variable_assignment>& new_possible_values) const
		{
			return h.template is_provable_by_axiom_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values);
		}
	};

	template<bool Contradiction, class OnSuccessFunc>
	bool unify_consequents(Formula* formula, array<Formula*>& quantifiers, Formula* consequent, array<Formula*>& second_quantifiers, array_map<Formula*, unsigned int>& antecedents, OnSuccessFunc on_success) const
	{
		if (consequent->type == FormulaType::AND) {
			for (unsigned int i = 0; i < consequent->array.length; i++)
				if (!unify_consequents<Contradiction>(formula, quantifiers, consequent->array.operands[i], second_quantifiers, antecedents, on_success)) return false;
			return true;
		} else if (consequent->type == FormulaType::FOR_ALL) {
			if (!second_quantifiers.add(consequent)) return false;
			if (!unify_consequents<Contradiction>(formula, quantifiers, consequent->quantifier.operand, second_quantifiers, antecedents, on_success))
				return false;
			second_quantifiers.length--;
			return true;
		} else if (consequent->type == FormulaType::IF_THEN) {
			if (!antecedents.put(consequent->binary.left, second_quantifiers.length))
				return false;
			if (!unify_consequents<Contradiction>(formula, quantifiers, consequent->binary.right, second_quantifiers, antecedents, on_success))
				return false;
			antecedents.size--;
			return true;
		} else {
			if (Contradiction) {
				if (consequent->type != FormulaType::NOT) return true;
				consequent = consequent->unary.operand;
			}

			array_map<Formula*, Term*> first_unifications(4);
			array_map<Formula*, Term*> second_unifications(4);
			if (unify(formula, consequent, quantifiers, first_unifications, second_quantifiers, second_unifications))
			{
				if (!on_success(antecedents, first_unifications, second_unifications))
					return false;
			}
			return true;
		}
	}

	template<bool Contradiction, typename Prover>
	bool is_provable_by_theorem_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<variable_assignment>& possible_values,
			array<variable_assignment>& new_possible_values,
			Prover& prover) const
	{
		unsigned int old_size = new_possible_values.length;
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i)) continue;

			array<Formula*> second_quantifiers(sets.sets[i].arity);
			Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				second_quantifiers.add(set_formula);
				set_formula = set_formula->quantifier.operand;
			}

			auto on_unify = [this,i,&possible_values,&second_quantifiers,&new_possible_values,&prover](array_map<Formula*, unsigned int>& antecedents, array_map<Formula*, Term*>& first_unifications, array_map<Formula*, Term*>& second_unifications)
			{
				typename Prover::ProvableElementArray provable_elements = prover.get_provable_elements(i);
				if (provable_elements.length == 0) {
					prover.free_provable_elements(provable_elements);
					return true;
				}

				array<variable_assignment> temp_possible_values(3);
				variable_assignment& values = temp_possible_values[0];
				if (!::init(values, second_quantifiers.length)) {
					prover.free_provable_elements(provable_elements);
					return false;
				}

				for (const auto& entry : second_unifications) {
					if (entry.value->type == TermType::CONSTANT) {
						core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
						values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
						values.assignment.values[entry.key->quantifier.variable - 1].constant = entry.value->constant;
					} else if (entry.value->type == TermType::NUMBER) {
						core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
						values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::NUMBER;
						values.assignment.values[entry.key->quantifier.variable - 1].number = entry.value->number;
					} else if (entry.value->type == TermType::STRING) {
						core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
						values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::STRING;
						if (!core::init(values.assignment.values[entry.key->quantifier.variable - 1].str, entry.value->str)) {
							prover.free_provable_elements(provable_elements);
							values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
							core::free(values); return false;
						}
					}
				}
				temp_possible_values.length = 1;

				for (unsigned int j = 0; j < possible_values.length; j++) {
					variable_assignment& src_values = temp_possible_values[1];
					variable_assignment& values = temp_possible_values[2];
					if (!::init(src_values, possible_values[j])) {
						prover.free_provable_elements(provable_elements);
						for (auto& element : temp_possible_values) core::free(element);
						return false;
					} else if (!::init(values, temp_possible_values[0])) {
						core::free(src_values);
						prover.free_provable_elements(provable_elements);
						for (auto& element : temp_possible_values) core::free(element);
						return false;
					}
					temp_possible_values.length = 3;

					bool unifies = true;
					for (const auto& entry : first_unifications) {
						if (entry.value->type == TermType::VARIABLE) {
							Term* term;
							if (possible_values[j].assignment.values[entry.key->quantifier.variable - 1].type == instantiation_type::CONSTANT) {
								term = Term::new_constant(possible_values[j].assignment.values[entry.key->quantifier.variable - 1].constant);
							} else if (possible_values[j].assignment.values[entry.key->quantifier.variable - 1].type == instantiation_type::NUMBER) {
								term = Term::new_number(possible_values[j].assignment.values[entry.key->quantifier.variable - 1].number);
							} else if (possible_values[j].assignment.values[entry.key->quantifier.variable - 1].type == instantiation_type::STRING) {
								term = Term::new_string(possible_values[j].assignment.values[entry.key->quantifier.variable - 1].str);
							} else {
								continue;
							}

							if (term == nullptr) {
								prover.free_provable_elements(provable_elements);
								for (auto& element : temp_possible_values) core::free(element);
								return false;
							} else if (!values.unify_value(entry.value->variable - 1, term)) {
								unifies = false;
								core::free(*term); core::free(term);
								break;
							}
							core::free(*term); core::free(term);
						} else {
							if (!src_values.unify_value(entry.key->quantifier.variable - 1, entry.value)) {
								unifies = false;
								break;
							}
						}
					}
					if (!unifies) {
						core::free(values); core::free(src_values);
						temp_possible_values.length = 1;
						continue;
					}

					array<variable_assignment> new_temp_possible_values(provable_elements.length);
					for (const tuple& element : provable_elements) {
						variable_assignment& new_values = new_temp_possible_values[new_temp_possible_values.length];
						if (!::init(new_values, values)) {
							prover.free_provable_elements(provable_elements);
							for (auto& element : temp_possible_values) core::free(element);
							return false;
						}

						bool unifies = true;
						for (unsigned int k = 0; k < element.length; k++) {
							Term* constant = nullptr;
							if (element[k].type == tuple_element_type::CONSTANT)
								constant = Term::new_constant(element[k].constant);
							else if (element[k].type == tuple_element_type::NUMBER)
								constant = Term::new_number(element[k].number);
							else if (element[k].type == tuple_element_type::STRING)
								constant = Term::new_string(element[k].str);
							if (constant == nullptr) {
								prover.free_provable_elements(provable_elements);
								for (auto& element : temp_possible_values) core::free(element);
								for (auto& element : new_temp_possible_values) core::free(element);
								core::free(new_values); return false;
							}
							if (!new_values.unify_value(k, constant)) {
								core::free(*constant); core::free(constant);
								unifies = false;
								break;
							}
							core::free(*constant); core::free(constant);
						}
						if (!unifies) {
							core::free(new_values);
							continue;
						}
						new_temp_possible_values.length++;
					}

					if (new_temp_possible_values.length == 0) {
						core::free(values); core::free(src_values);
						temp_possible_values.length = 1;
						continue;
					}
					insertion_sort(new_temp_possible_values, default_sorter());
					unique_and_cleanup(new_temp_possible_values);

					/* first restrict the set of possible_values to only those that can prove the antecedents */
					array<Formula*> antecedent_quantifiers(second_quantifiers.capacity);
					for (const auto& entry : antecedents) {
						Formula* antecedent = entry.key;
						unsigned int num_quantifiers = entry.value;
						for (unsigned int i = 0; i < num_quantifiers; i++)
							antecedent_quantifiers[antecedent_quantifiers.length++] = second_quantifiers[i];
						if (!is_provable_without_abduction<false>(antecedent, antecedent_quantifiers, new_temp_possible_values, prover))
							break;
					}
					if (new_temp_possible_values.length == 0) {
						core::free(values); core::free(src_values);
						temp_possible_values.length = 1;
						continue;
					}

					/* we can prove all antecedents; now we can use this consequent to prove `formula` */
					if (!new_possible_values.ensure_capacity(new_possible_values.length + new_temp_possible_values.length + 1)) {
						prover.free_provable_elements(provable_elements);
						for (auto& element : temp_possible_values) core::free(element);
						for (auto& element : new_temp_possible_values) core::free(element);
						return false;
					}
					for (unsigned int k = 0; k < new_temp_possible_values.length; k++) {
						variable_assignment& new_values = new_possible_values[new_possible_values.length];
						if (!::init(new_values, src_values)) {
							prover.free_provable_elements(provable_elements);
							for (auto& element : temp_possible_values) core::free(element);
							for (auto& element : new_temp_possible_values) core::free(element);
							return false;
						}

						bool unifies = true;
						for (const auto& unification : first_unifications) {
							if (unification.value->type != TermType::VARIABLE)
								continue;

							Term* term;
							if (new_temp_possible_values[k].assignment.values[unification.value->variable - 1].type == instantiation_type::CONSTANT) {
								term = Term::new_constant(new_temp_possible_values[k].assignment.values[unification.value->variable - 1].constant);
							} else if (new_temp_possible_values[k].assignment.values[unification.value->variable - 1].type == instantiation_type::NUMBER) {
								term = Term::new_number(new_temp_possible_values[k].assignment.values[unification.value->variable - 1].number);
							} else if (new_temp_possible_values[k].assignment.values[unification.value->variable - 1].type == instantiation_type::STRING) {
								term = Term::new_string(new_temp_possible_values[k].assignment.values[unification.value->variable - 1].str);
							} else {
								continue;
							}

							if (term == nullptr) {
								prover.free_provable_elements(provable_elements);
								for (auto& element : temp_possible_values) core::free(element);
								for (auto& element : new_temp_possible_values) core::free(element);
								core::free(new_values); return false;
							} else if (!new_values.unify_value(unification.key->quantifier.variable - 1, term)) {
								unifies = false;
								core::free(*term); core::free(term);
								break;
							}
							core::free(*term); core::free(term);
						}
						if (!unifies) {
							core::free(new_values);
							continue;
						}
						new_possible_values.length++;
					}
					for (auto& element : new_temp_possible_values) core::free(element);

					core::free(values);
					core::free(src_values);
					temp_possible_values.length = 1;
				}
				prover.free_provable_elements(provable_elements);
				for (auto& element : temp_possible_values) core::free(element);
				return true;
			};

			array_map<Formula*, unsigned int> antecedents(8);
/*debug_counter++;
if (debug_flag) {
fprintf(stderr, "[DEBUG: %u] calling unify_consequents:\n", debug_counter);
print("  formula: ", stderr); print(*formula, stderr, *debug_terminal_printer); print('\n', stderr);
print("  set_formula: ", stderr); print(*set_formula, stderr, *debug_terminal_printer); print('\n', stderr);
}*/
			if (!unify_consequents<Contradiction>(formula, quantifiers, set_formula, second_quantifiers, antecedents, on_unify))
				return false;
		}
		for (unsigned int i = 0; i < implication_axioms.length; i++) {
			if (!prover.is_antecedent_satisfied(i)) continue;
			Proof* axiom = implication_axioms[i].key;
			Formula* consequent = axiom->formula->binary.right;

			array_view<Formula*> consequent_conjuncts(
					(consequent->type == FormulaType::AND) ? consequent->array.operands : &consequent,
					(consequent->type == FormulaType::AND) ? consequent->array.length : 1);
			for (unsigned int l = 0; l < consequent_conjuncts.length; l++) {
				Formula* consequent_conjunct = consequent_conjuncts.array[l];
				if (Contradiction) {
					if (consequent_conjunct->type != FormulaType::NOT) continue;
					consequent_conjunct = consequent_conjunct->unary.operand;
				}
				array_map<Formula*, Term*> first_unifications(4);
				array_map<Formula*, Term*> second_unifications(4);
				array<Formula*> second_quantifiers(1);
				if (!unify(formula, consequent_conjunct, quantifiers, first_unifications, second_quantifiers, second_unifications))
					continue;

				for (unsigned int j = 0; j < possible_values.length; j++) {
					const variable_assignment& values = possible_values[j];
					if (!new_possible_values.ensure_capacity(new_possible_values.length + 1))
						return false;
					variable_assignment& new_value = new_possible_values[new_possible_values.length];
					if (!::init(new_value, values))
						return false;

					bool unifies = true;
					for (const auto& unification : first_unifications) {
						if (unification.value->type != TermType::VARIABLE
						 && !new_value.unify_value(unification.key->quantifier.variable - 1, unification.value)) {
							unifies = false;
							break;
						}
					}
					if (!unifies) {
						core::free(new_value);
						continue;
					}

					new_possible_values.length++;
				}
			}
		}
		if (new_possible_values.length > old_size) {
			insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
			array<variable_assignment> temp(new_possible_values.length);
			set_union(temp.data, temp.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
			swap(temp, new_possible_values);
			for (auto& element : temp) core::free(element);
		}
		return true;
	}

	template<bool Contradiction, typename Prover>
	bool is_provable_by_exclusion_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			array<variable_assignment>& possible_values,
			array<variable_assignment>& new_possible_values,
			Prover& prover) const
	{
		/* for performance, we cutoff this search to a maximum depth of 1; TODO: how does this limit provability? */
		if (prover.h.set_ids.length > 1)
			return true;

		/* also for performance, we only consider contradiction proofs of existentially quantified expressions */
		if (!(formula->type == FormulaType::EXISTS || (formula->type == FormulaType::UNARY_APPLICATION && (formula->binary.left->type == FormulaType::VARIABLE || formula->binary.left->type == FormulaType::CONSTANT))))
			return true;

		unsigned int old_size = new_possible_values.length;
		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr
			 || prover.is_set_removed(i)
			 || prover.is_set_on_stack(i)) continue;

			typename Prover::ProvableElementArray provable_elements = prover.get_provable_elements(i);
			if (provable_elements.length != sets.sets[i].set_size) {
				prover.free_provable_elements(provable_elements);
				continue;
			}

			array<Formula*> second_quantifiers(sets.sets[i].arity);
			Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				second_quantifiers.add(set_formula);
				set_formula = set_formula->quantifier.operand;
			}

			array<unification> unifications(4);
			if (!unify_subformula<Contradiction>(formula, set_formula, quantifiers, second_quantifiers, unifications) || unifications.length == 0) {
				prover.free_provable_elements(provable_elements);
				for (unification& u : unifications) core::free(u);
				continue;
			}

			/* optimization: if `set_formula` contains a name scope as a
			   conjunct, then we avoid searching over this set if
			   `possible_values` only contains constants that don't provably
			   have the same name */
			/* TODO: i think this optimization is not correct? */
			/*if (unifications.length == 1 && unifications[0].first_unifications.size == 1
			 && unifications[0].first_unifications.values[0]->type == TermType::VARIABLE
			 && possible_values.length == 1)
			{
				unsigned int src_variable = unifications[0].first_unifications.keys[0]->quantifier.variable;
				unsigned int dst_variable = unifications[0].first_unifications.values[0]->variable;

				if (possible_values[0].assignment.values[src_variable - 1].type == instantiation_type::CONSTANT) {
					bool found_matching_name = false;
					unsigned int src_constant = possible_values[0].assignment.values[src_variable - 1].constant;
					if (set_formula->type == FormulaType::AND) {
						for (unsigned int i = 0; i < set_formula->array.length; i++) {
							const Term* arg1; const Term* arg2;
							if (!is_name_scope(set_formula->array.operands[i], arg1, arg2))
								continue;
							if (arg1->type == TermType::VARIABLE && arg1->variable == dst_variable
							 && arg2->type == TermType::STRING && is_concept_name(src_constant, arg2->str))
							{
								found_matching_name = true;
							}
							break;
						}
					} else {
						const Term* arg1; const Term* arg2;
						if (is_name_scope(set_formula, arg1, arg2)
						 && arg1->type == TermType::VARIABLE && arg1->variable == dst_variable
						 && arg2->type == TermType::STRING && is_concept_name(src_constant, arg2->str))
						{
							found_matching_name = true;
						}
					}
					if (!found_matching_name) {
						prover.free_provable_elements(provable_elements);
						for (unification& u : unifications) core::free(u);
						continue;
					}
				}
			}*/

			unsigned int pushed_axiom_count;
			if (!prover.template push_axiom<!Contradiction>(i, formula, possible_values, pushed_axiom_count, unifications)) {
				prover.free_provable_elements(provable_elements);
				for (unification& u : unifications) core::free(u);
				return false;
			} else if (pushed_axiom_count == 0) {
				prover.free_provable_elements(provable_elements);
				for (unification& u : unifications) core::free(u);
				continue;
			}
			for (unification& u : unifications) core::free(u);

			array<variable_assignment> temp_possible_values(4);
			if (!::init(temp_possible_values[0], sets.sets[i].arity)) {
				prover.free_provable_elements(provable_elements);
				return false;
			}
			temp_possible_values.length = 1;
			is_provable_without_abduction<false>(set_formula, second_quantifiers, temp_possible_values, prover);

			for (const variable_assignment& value : temp_possible_values) {
				bool proved_from_this_axiom = false;
				for (const auto& entry : value.axiom_assignments) {
					if (entry.key >= prover.h.axioms.length - pushed_axiom_count) {
						proved_from_this_axiom = true;
						break;
					}
				}
				if (!proved_from_this_axiom) continue;

				tuple new_tuple;
				if (!::init(new_tuple, sets.sets[i].arity)) {
					prover.free_provable_elements(provable_elements);
					for (auto& element : temp_possible_values) core::free(element);
					return false;
				}
				bool has_any = false;
				for (uint_fast8_t k = 0; k < value.assignment.length; k++) {
					if (value.assignment.values[k].type == instantiation_type::CONSTANT) {
						new_tuple[k].type = tuple_element_type::CONSTANT;
						new_tuple[k].constant = value.assignment.values[k].constant;
					} else if (value.assignment.values[k].type == instantiation_type::NUMBER) {
						new_tuple[k].type = tuple_element_type::NUMBER;
						new_tuple[k].number = value.assignment.values[k].number;
					} else if (value.assignment.values[k].type == instantiation_type::STRING) {
						new_tuple[k].type = tuple_element_type::STRING;
						if (!core::init(new_tuple[k].str, value.assignment.values[k].str)) {
							prover.free_provable_elements(provable_elements);
							for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
							core::free(new_tuple.elements);
							for (auto& element : temp_possible_values) core::free(element);
							return false;
						}
					} else if (value.assignment.values[k].type == instantiation_type::ANY || value.assignment.values[k].type == instantiation_type::ANY_NUMBER) {
						/* we assume the universe has sufficiently many concepts that `provable_elements` is not a superset of the universe */
						for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
						core::free(new_tuple.elements);
						has_any = true;
						break;
					}
				}

				if (has_any || !provable_elements.contains(new_tuple)) {
					if (!has_any) core::free(new_tuple);
					/* (when Contradiction is false) if `formula` is false, then we have a contradiction, so `formula` must be true;
					   (when Contradiction is true) if `formula` is true, then we have a contradiction, so `formula` must be false */
					if (!new_possible_values.ensure_capacity(new_possible_values.length + 1)
					 || !prover.get_axiom_instantiation(new_possible_values[new_possible_values.length], value, pushed_axiom_count))
					{
						prover.free_provable_elements(provable_elements);
						for (auto& element : temp_possible_values) core::free(element);
						return false;
					}
					new_possible_values.length++;
				} else {
					core::free(new_tuple);
				}
			}
			prover.pop_axiom(pushed_axiom_count);
			prover.free_provable_elements(provable_elements);
			for (auto& element : temp_possible_values) core::free(element);
		}
		if (new_possible_values.length > old_size) {
			insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
			array<variable_assignment> temp(new_possible_values.length);
			set_union(temp.data, temp.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
			swap(temp, new_possible_values);
			for (auto& element : temp) core::free(element);
		}
		return true;
	}

	template<typename Prover>
	bool for_all_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<variable_assignment>& possible_values,
			array<variable_assignment>& new_possible_values,
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
				if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i) || sets.sets[i].arity < new_variables.length) continue;

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

				typename Prover::ProvableElementArray provable_elements = prover.get_provable_elements(i);
				if (provable_elements.length != sets.sets[i].set_size) {
					prover.free_provable_elements(provable_elements);
					continue;
				}

				for (unsigned int i = 0; i < possible_values.length; i++) {
					const variable_assignment& values = possible_values[i];
					array<variable_assignment> temp_possible_values(4);
					variable_assignment& new_values = temp_possible_values[0];
					if (!::init(new_values, values)) {
						prover.free_provable_elements(provable_elements);
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
								prover.free_provable_elements(provable_elements);
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
						array<variable_assignment> temp(new_possible_values.length + temp_possible_values.length);
						set_union(temp, new_possible_values, temp_possible_values);
						swap(temp, new_possible_values);
						for (auto& element : temp_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
					}
				}
				prover.free_provable_elements(provable_elements);
			}
		}
		return true;
	}

	inline bool does_conjunct_define_existential(Term* left, Term* right,
			Term** src_variables, Term** dst_variables, unsigned int& substitution_count) const
	{
		if (right->type == TermType::CONSTANT || right->type == TermType::STRING || right->type == TermType::NUMBER || right->type == TermType::VARIABLE) {
			src_variables[substitution_count] = left;
			dst_variables[substitution_count] = right;
			for (unsigned int j = 0; j < substitution_count; j++) {
				if (dst_variables[j] == left || *dst_variables[j] == *left)
					dst_variables[j] = right;
			}
			substitution_count++;
			return true;
		} else {
			array_map<unsigned int, unsigned int> variable_map(16);
			Formula* canonicalized = Canonicalizer::canonicalize(*right, variable_map);

			/* check if any definitions match the right-hand side */
			if (canonicalized->type == TermType::LAMBDA && canonicalized->quantifier.operand->type == TermType::UNARY_APPLICATION
			 && canonicalized->quantifier.operand->binary.left->type == TermType::CONSTANT
			 && canonicalized->quantifier.operand->binary.left->constant >= new_constant_offset
			 && canonicalized->quantifier.operand->binary.right->type == TermType::VARIABLE
			 && canonicalized->quantifier.operand->binary.right->variable == canonicalized->quantifier.variable)
			{
				Proof* axiom = ground_concepts[canonicalized->quantifier.operand->binary.left->constant - new_constant_offset].definitions[0];
				src_variables[substitution_count] = left;
				dst_variables[substitution_count] = axiom->formula->binary.left;
				for (unsigned int j = 0; j < substitution_count; j++) {
					if (dst_variables[j] == left || *dst_variables[j] == *left)
						dst_variables[j] = axiom->formula->binary.left;
				}
				substitution_count++;
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return true;
			}

			bool contains;
			unsigned int constant = reverse_definitions.get(*canonicalized, contains);
			if (contains) {
				Proof* axiom = ground_concepts[constant - new_constant_offset].definitions[0];
				src_variables[substitution_count] = left;
				dst_variables[substitution_count] = axiom->formula->binary.left;
				for (unsigned int j = 0; j < substitution_count; j++) {
					if (dst_variables[j] == left || *dst_variables[j] == *left)
						dst_variables[j] = axiom->formula->binary.left;
				}
				substitution_count++;
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return true;
			}
			core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
		}
		return false;
	}

	template<typename Prover>
	inline bool get_elements_of_full_supersets(Term* subset_term,
			array<unsigned int>& quantified_variables,
			array_map<Term*, array<tuple>>& elements,
			Prover& prover, bool& found_empty_set) const
	{
		if (subset_term->type != TermType::UNARY_APPLICATION
		 || subset_term->binary.left->type != TermType::CONSTANT
		 || (new_constant_offset > subset_term->binary.left->constant)
		 || subset_term->binary.right->type != TermType::VARIABLE
		 || !quantified_variables.contains(subset_term->binary.right->variable))
			return true;
		for (const auto& entry : elements)
			if (entry.key->variable == subset_term->binary.right->variable) return true;
		const concept<ProofCalculus>& c = ground_concepts[subset_term->binary.left->constant - new_constant_offset];
		Formula* set_formula = c.definitions[0]->formula->binary.right->quantifier.operand;

		bool contains;
		unsigned int set_id = sets.set_ids.get(*set_formula, contains);
		if (!contains || sets.sets[set_id].arity != 1) return true;

		elements.keys[elements.size] = subset_term->binary.right;
		if (sets.sets[set_id].set_size == 0) {
			found_empty_set = true;
			return true;
		} else if (sets.sets[set_id].provable_elements.length != sets.sets[set_id].set_size) {
			return true;
		}
		array<tuple>& new_elements = elements.values[elements.size];
		if (!array_init(new_elements, max(1, sets.sets[set_id].provable_elements.length))) {
			return false;
		} for (const tuple& tup : sets.sets[set_id].provable_elements) {
			if (!::init(new_elements[new_elements.length], tup)) {
				for (tuple& tup : new_elements) core::free(tup);
				core::free(new_elements);
				return false;
			}
			new_elements.length++;
		}
		quantified_variables.remove(quantified_variables.index_of(elements.keys[elements.size]->variable));
		elements.size++;
		return true;
	}

	template<typename Prover>
	bool not_exists_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<variable_assignment>& possible_values,
			array<variable_assignment>& new_possible_values,
			Prover& prover) const
	{
		Formula* operand = formula;
		unsigned int old_quantifier_length = quantifiers.length;
		array<unsigned int> quantified_variables(4);
		while (operand->type == FormulaType::EXISTS) {
			if (!quantifiers.ensure_capacity(operand->quantifier.variable)
			 || !quantified_variables.add(operand->quantifier.variable))
				return false;
			while (operand->quantifier.variable - 1 >= quantifiers.length)
				quantifiers[quantifiers.length++] = nullptr;
			quantifiers[operand->quantifier.variable - 1] = operand;
			operand = operand->quantifier.operand;
		}

		/* look for terms of the form `x=A` in the operand which would fix the value of the quantified variable `x` */
		unsigned int substitution_count = 0;
		Term** src_variables = (Term**) alloca(sizeof(Term*) * (quantifiers.length - old_quantifier_length) * 2);
		Term** dst_variables = src_variables + (quantifiers.length - old_quantifier_length);
		array<Formula*> new_operands(operand->type == FormulaType::AND ? operand->array.length : 1);
		if (operand->type == FormulaType::AND) {
			for (unsigned int i = 0; i < operand->array.length; i++) {
				if (operand->array.operands[i]->type == TermType::EQUALS) {
					Term* left = operand->array.operands[i]->binary.left;
					Term* right = operand->array.operands[i]->binary.right;
					if (left->type == TermType::VARIABLE && quantified_variables.contains(left->variable)) {
						if (!does_conjunct_define_existential(left, right, src_variables, dst_variables, substitution_count)) {
							if (right->type == TermType::VARIABLE && quantified_variables.contains(right->variable)) {
								if (!does_conjunct_define_existential(right, left, src_variables, dst_variables, substitution_count))
									new_operands[new_operands.length++] = operand->array.operands[i];
								else quantified_variables.remove(quantified_variables.index_of(right->variable));
							} else {
								new_operands[new_operands.length++] = operand->array.operands[i];
							}
						} else {
							quantified_variables.remove(quantified_variables.index_of(left->variable));
						}
					} else if (right->type == TermType::VARIABLE && quantified_variables.contains(right->variable)) {
						if (!does_conjunct_define_existential(right, left, src_variables, dst_variables, substitution_count))
							new_operands[new_operands.length++] = operand->array.operands[i];
						else quantified_variables.remove(quantified_variables.index_of(right->variable));
					} else {
						new_operands[new_operands.length++] = operand->array.operands[i];
					}
				} else {
					new_operands[new_operands.length++] = operand->array.operands[i];
				}
			}
		} else {
			if (operand->type == TermType::EQUALS) {
				Term* left = operand->binary.left;
				Term* right = operand->binary.right;
				if (left->type == TermType::VARIABLE && quantified_variables.contains(left->variable)) {
					if (!does_conjunct_define_existential(left, right, src_variables, dst_variables, substitution_count)) {
						if (right->type == TermType::VARIABLE && quantified_variables.contains(right->variable)) {
							if (!does_conjunct_define_existential(right, left, src_variables, dst_variables, substitution_count))
								new_operands[new_operands.length++] = operand;
							else quantified_variables.remove(quantified_variables.index_of(right->variable));
						} else {
							new_operands[new_operands.length++] = operand;
						}
					} else {
						quantified_variables.remove(quantified_variables.index_of(left->variable));
					}
				} else if (right->type == TermType::VARIABLE && quantified_variables.contains(right->variable)) {
					if (!does_conjunct_define_existential(right, left, src_variables, dst_variables, substitution_count))
						new_operands[new_operands.length++] = operand;
					else quantified_variables.remove(quantified_variables.index_of(right->variable));
				} else {
					new_operands[new_operands.length++] = operand;
				}
			} else {
				new_operands[new_operands.length++] = operand;
			}
		}

		if (new_operands.length == 0) {
			/* all terms in `operand` are provably true */
			for (auto& element : new_possible_values) core::free(element);
			new_possible_values.clear();
			return true;
		}

		for (unsigned int i = 0; i < new_operands.length; i++) {
			Term* new_operand = substitute_all(new_operands[i], src_variables, dst_variables, substitution_count);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) {
					core::free(*new_operands[j]); if (new_operands[j]->reference_count == 0) core::free(new_operands[j]);
				}
				return false;
			}
			new_operands[i] = new_operand;
		}

		Term* new_operand;
		if (new_operands.length == 1) {
			new_operand = new_operands[0];
		} else {
			new_operand = Formula::new_and(make_array_view(new_operands.data, new_operands.length));
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < new_operands.length; j++) {
					core::free(*new_operands[j]); if (new_operands[j]->reference_count == 0) core::free(new_operands[j]);
				}
				return false;
			}
		}

		/* consider existential quantification over the elements in full sets */
		bool found_empty_set = false;
		array_map<Term*, array<tuple>> elements(max(1, quantified_variables.length));
		if (new_operand->type == FormulaType::AND) {
			for (unsigned int i = 0; i < new_operand->array.length; i++) {
				if (!get_elements_of_full_supersets(new_operand->array.operands[i], quantified_variables, elements, prover, found_empty_set)) {
					for (auto entry : elements) {
						for (tuple& tup : entry.value) core::free(tup);
						core::free(entry.value);
					}
					core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
					return false;
				} else if (found_empty_set) {
					for (auto entry : elements) {
						for (tuple& tup : entry.value) core::free(tup);
						core::free(entry.value);
					}
					core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
					for (auto& element : new_possible_values) core::free(element);
					new_possible_values.clear();
					if (!new_possible_values.ensure_capacity(possible_values.length)) return false;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						if (!::init(new_possible_values[i], possible_values[i])) return false;
						new_possible_values.length++;
					}
					return true;
				}
			}
		} else {
			if (!get_elements_of_full_supersets(new_operand, quantified_variables, elements, prover, found_empty_set)) {
				for (auto entry : elements) {
					for (tuple& tup : entry.value) core::free(tup);
					core::free(entry.value);
				}
				core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
				return false;
			} else if (found_empty_set) {
				for (auto entry : elements) {
					for (tuple& tup : entry.value) core::free(tup);
					core::free(entry.value);
				}
				core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
				for (auto& element : new_possible_values) core::free(element);
				new_possible_values.clear();
				if (!new_possible_values.ensure_capacity(possible_values.length)) return false;
				for (unsigned int i = 0; i < possible_values.length; i++) {
					if (!::init(new_possible_values[i], possible_values[i])) return false;
					new_possible_values.length++;
				}
				return true;
			}
		}
		if (quantified_variables.length == 0) {
			array<variable_assignment> temp_possible_values(possible_values.length);
			for (unsigned int i = 0; i < possible_values.length; i++) {
				if (!::init(temp_possible_values[i], possible_values[i]))
					return false;
				temp_possible_values.length++;
			}
			bool failure = false;
			apply_to_cartesian_product(elements.values, elements.size, [this,&quantifiers,&temp_possible_values,&prover,&elements,new_operand,&failure](const unsigned int* index_array)
			{
				if (elements.size == 0)
					return is_provable_without_abduction<true>(new_operand, quantifiers, temp_possible_values, prover);

				Term** dst_variables = (Term**) alloca(sizeof(Term*) * elements.size);
				for (unsigned int i = 0; i < elements.size; i++) {
					const tuple& value = elements.values[i][index_array[i]];
					if (value[0].type == tuple_element_type::CONSTANT) {
						dst_variables[i] = Formula::new_constant(value[0].constant);
					} else if (value[0].type == tuple_element_type::NUMBER) {
						dst_variables[i] = Formula::new_number(value[0].number);
					} else if (value[0].type == tuple_element_type::STRING) {
						dst_variables[i] = Formula::new_string(value[0].str);
					}
					if (dst_variables[i] == nullptr) {
						for (unsigned int j = 0; j < i; j++) {
							core::free(*dst_variables[j]); core::free(dst_variables[j]);
						}
						failure = true;
						return false;
					}
				}

				Formula* substituted_operand = substitute_all(new_operand, elements.keys, dst_variables, elements.size);
				for (unsigned int j = 0; j < elements.size; j++) {
					core::free(*dst_variables[j]);
					if (dst_variables[j]->reference_count == 0)
						core::free(dst_variables[j]);
				}
				if (substituted_operand == nullptr) {
					failure = true;
					return false;
				}
				if (!is_provable_without_abduction<true>(substituted_operand, quantifiers, temp_possible_values, prover)) {
					for (auto& element : temp_possible_values) core::free(element);
					core::free(*substituted_operand); if (substituted_operand->reference_count == 0) core::free(substituted_operand);
					return false;
				}
				core::free(*substituted_operand); if (substituted_operand->reference_count == 0) core::free(substituted_operand);
				return true;
			});
			if (failure) {
				for (auto entry : elements) {
					for (tuple& tup : entry.value) core::free(tup);
					core::free(entry.value);
				}
				core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
				return false;
			}
			if (temp_possible_values.length != 0) {
				array<variable_assignment> temp(new_possible_values.length + temp_possible_values.length);
				set_union(temp, new_possible_values, temp_possible_values);
				swap(temp, new_possible_values);
				for (auto& element : temp_possible_values) core::free(element);
				for (auto& element : temp) core::free(element);
			}
		}
		for (auto entry : elements) {
			for (tuple& tup : entry.value) core::free(tup);
			core::free(entry.value);
		}

		/* check if this is a definition of a set and a declaration of it's size */
		if (new_operand->type == FormulaType::AND) {
			unsigned int set_definition_index = new_operand->array.length;
			unsigned int set_size_index = new_operand->array.length;
			unsigned int element_index = new_operand->array.length;
			for (unsigned int i = 0; i < new_operand->array.length; i++) {
				hol_term* conjunct = new_operand->array.operands[i];
				if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::VARIABLE
				 && conjunct->binary.left->variable == formula->quantifier.variable
				 && conjunct->binary.right->type == TermType::LAMBDA)
				{
					set_definition_index = i;
				} else if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::UNARY_APPLICATION
						&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
						&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
						&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
						&& conjunct->binary.left->binary.right->variable == formula->quantifier.variable)
				{
					set_size_index = i;
				} else if (conjunct->type == TermType::UNARY_APPLICATION
						&& conjunct->binary.left->type == TermType::VARIABLE
						&& conjunct->binary.left->variable == formula->quantifier.variable
						&& conjunct->binary.right->type == TermType::VARIABLE)
				{
					element_index = i;
				}
			}

			if (set_definition_index != new_operand->array.length && set_size_index != new_operand->array.length) {
				hol_term* set_definition = new_operand->array.operands[set_definition_index];
				hol_term* set_size = new_operand->array.operands[set_size_index];
				unsigned int temp_quantifier_count = quantifiers.length;
				hol_term* inner_set_formula = set_definition->binary.right;
				while (inner_set_formula->type == FormulaType::LAMBDA) {
					if (!quantifiers.ensure_capacity(inner_set_formula->quantifier.variable)) {
						quantifiers.length = old_quantifier_length;
						core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
						return false;
					}
					while (inner_set_formula->quantifier.variable - 1 >= quantifiers.length)
						quantifiers[quantifiers.length++] = nullptr;
					quantifiers[inner_set_formula->quantifier.variable - 1] = inner_set_formula;
					inner_set_formula = inner_set_formula->quantifier.operand;
				}

				array<variable_assignment> temp_possible_values(possible_values.length);
				for (const variable_assignment& assignment : possible_values) {
					variable_assignment& new_value = temp_possible_values[temp_possible_values.length];
					if (!::init(new_value, assignment, quantifiers.length)) {
						quantifiers.length = old_quantifier_length;
						core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
						return false;
					}
					temp_possible_values.length++;
				}

				is_provable_without_abduction<false>(inner_set_formula, quantifiers, temp_possible_values, prover);
				quantifiers.length = temp_quantifier_count;

				array<variable_assignment> sub_possible_values(possible_values.length);
				array<array<variable_assignment>> provable_elements(possible_values.length);
				for (const variable_assignment& assignment : temp_possible_values) {
					/* find a matching variable assignment */
					unsigned int index;
					for (index = 0; index < sub_possible_values.length; index++) {
						bool match = true;
						for (uint_fast8_t k = 0; match && k < temp_quantifier_count; k++)
							if (assignment.assignment.values[k] != sub_possible_values[index].assignment.values[k]) match = false;
						if (match) break;
					}

					if (index == sub_possible_values.length) {
						if (!sub_possible_values.ensure_capacity(sub_possible_values.length + 1)
						 || !provable_elements.ensure_capacity(provable_elements.length + 1)
						 || !::init(sub_possible_values[sub_possible_values.length], assignment, temp_quantifier_count))
						{
							quantifiers.length = old_quantifier_length;
							core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
							for (auto& element : temp_possible_values) core::free(element);
							for (auto& element : sub_possible_values) core::free(element);
							for (array<variable_assignment>& elements : provable_elements) {
								for (auto& element : elements) core::free(element);
								core::free(elements);
							}
							return false;
						} else if (!array_init(provable_elements[provable_elements.length], 4)) {
							core::free(sub_possible_values[sub_possible_values.length]);
							quantifiers.length = old_quantifier_length;
							core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
							for (auto& element : temp_possible_values) core::free(element);
							for (auto& element : sub_possible_values) core::free(element);
							for (array<variable_assignment>& elements : provable_elements) {
								for (auto& element : elements) core::free(element);
								core::free(elements);
							}
							return false;
						}
						sub_possible_values.length++;
						provable_elements.length++;
					}

					if (!provable_elements[index].ensure_capacity(provable_elements[index].length + 1)
					 || !::init(provable_elements[index][provable_elements[index].length], assignment))
					{
						quantifiers.length = old_quantifier_length;
						core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
						for (auto& element : temp_possible_values) core::free(element);
						for (auto& element : sub_possible_values) core::free(element);
						for (array<variable_assignment>& elements : provable_elements) {
							for (auto& element : elements) core::free(element);
							core::free(elements);
						}
						return false;
					}
					provable_elements[index].length++;
				}
				for (auto& element : temp_possible_values) core::free(element);

				for (unsigned int i = 0; i < sub_possible_values.length; i++) {
					unsigned int provable_element_count = 0;
					for (const variable_assignment& element : provable_elements[i]) {
						bool has_any = false;
						for (uint_fast8_t k = temp_quantifier_count; !has_any && k < element.assignment.length; k++)
							if (element.assignment.values[k].type == instantiation_type::ANY || element.assignment.values[k].type == instantiation_type::ANY_NUMBER) has_any = true;
						if (has_any) {
							provable_element_count = UINT_MAX;
							break;
						} else {
							provable_element_count++;
						}
					}

					if (element_index != new_operand->array.length && provable_element_count != UINT_MAX) {
						hol_term* element = new_operand->array.operands[element_index];
						const instantiation& extra_provable_element = sub_possible_values[i].assignment.values[element->binary.right->variable - 1];

						/* check if `extra_provable_element` is already an element of `provable_elements[i]` */
						bool is_element = false;
						for (const variable_assignment& element : provable_elements[i]) {
							if (element.assignment.values[temp_quantifier_count] == extra_provable_element) {
								is_element = true;
								break;
							}
						}
						if (!is_element)
							provable_element_count++;
					}

					if (provable_element_count == UINT_MAX) {
						core::free(sub_possible_values[i]);
						for (auto& element : provable_elements[i]) core::free(element);
						core::free(provable_elements[i]);
						sub_possible_values.remove(i);
						provable_elements.remove(i--);
						continue;
					}
					if (set_size->binary.right->type == TermType::NUMBER) {
						if (provable_element_count <= set_size->binary.right->number.integer) {
							core::free(sub_possible_values[i]);
							for (auto& element : provable_elements[i]) core::free(element);
							core::free(provable_elements[i]);
							sub_possible_values.remove(i);
							provable_elements.remove(i--);
						}
					} else if (set_size->binary.right->type == TermType::VARIABLE) {
						Term* num = Term::new_number(provable_element_count, 0);
						if (num == nullptr) {
							quantifiers.length = old_quantifier_length;
							core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
							for (auto& element : sub_possible_values) core::free(element);
							for (array<variable_assignment>& elements : provable_elements) {
								for (auto& element : elements) core::free(element);
								core::free(elements);
							}
							return false;
						}
						if (!sub_possible_values[i].unify_greater_than_or_equal(set_size->binary.right->variable - 1, num)) {
							core::free(sub_possible_values[i]);
							for (auto& element : provable_elements[i]) core::free(element);
							core::free(provable_elements[i]);
							sub_possible_values.remove(i);
							provable_elements.remove(i--);
						}
						core::free(*num); core::free(num);
					}
				}
				if (sub_possible_values.length != 0) {
					insertion_sort(sub_possible_values, default_sorter());
					array<variable_assignment> temp(new_possible_values.length + sub_possible_values.length);
					set_union(temp, new_possible_values, sub_possible_values);
					swap(temp, new_possible_values);
					for (auto& element : sub_possible_values) core::free(element);
					for (auto& element : temp) core::free(element);
				}

				/* check if the internal set unifies with an existing set */
				for (unsigned int i = 1; i < sets.set_count + 1; i++) {
					if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i)) continue;

					array<Formula*> second_quantifiers(sets.sets[i].arity);
					Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

					array_map<Formula*, Term*> first_unifications(4);
					array_map<Formula*, Term*> second_unifications(4);
					if (!unify(set_definition->binary.right, set_formula, quantifiers, first_unifications, second_quantifiers, second_unifications))
						continue;

					array<variable_assignment> temp_possible_values(possible_values.length);
					for (unsigned int j = 0; j < possible_values.length; j++) {
						const variable_assignment& values = possible_values[j];
						variable_assignment& new_values = temp_possible_values[temp_possible_values.length];
						if (!::init(new_values, values)) {
							quantifiers.length = old_quantifier_length;
							core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
							for (auto& element : temp_possible_values) core::free(element);
							return false;
						}

						bool unifies = true;
						for (const auto& unification : first_unifications) {
							if (unification.key->quantifier.variable <= old_quantifier_length
							 && !new_values.unify_value(unification.key->quantifier.variable - 1, unification.value))
							{
								unifies = false;
								break;
							}
						}
						if (!unifies) {
							core::free(new_values);
							continue;
						}

						if (set_size->binary.right->type == TermType::NUMBER) {
							if (sets.sets[i].set_size == set_size->binary.right->number.integer) {
								core::free(new_values);
								continue;
							}
						} else if (set_size->binary.right->type == TermType::VARIABLE) {
							Term* num = Term::new_number(sets.sets[i].set_size, 0);
							if (num == nullptr) {
								core::free(new_values);
								quantifiers.length = old_quantifier_length;
								core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
								for (auto& element : temp_possible_values) core::free(element);
								return false;
							}
							if (!new_values.antiunify_value(set_size->binary.right->variable - 1, num)) {
								core::free(new_values);
								continue;
							}
							core::free(*num); core::free(num);
						}

						temp_possible_values.length++;
					}

					if (temp_possible_values.length != 0) {
						array<variable_assignment> temp(new_possible_values.length + temp_possible_values.length);
						set_union(temp, new_possible_values, temp_possible_values);
						swap(temp, new_possible_values);
						for (auto& element : temp_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
					}
				}
			}
		}

		for (unsigned int i = 1; i < sets.set_count + 1; i++) {
			if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i) || sets.sets[i].set_size != 0) continue;

			array<Formula*> second_quantifiers(sets.sets[i].arity);
			Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
			for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
				second_quantifiers.add(set_formula);
				set_formula = set_formula->quantifier.operand;
			}

			Formula* set_operand = set_formula;
			while (set_operand->type == FormulaType::EXISTS) {
				second_quantifiers.add(set_operand);
				set_operand = set_operand->quantifier.operand;
			}

			array_map<Formula*, Term*> first_unifications(4);
			array_map<Formula*, Term*> second_unifications(4);
			if (!unify(new_operand, set_operand, quantifiers, first_unifications, second_quantifiers, second_unifications))
				continue;

			array<variable_assignment> temp_possible_values(possible_values.length);
			for (unsigned int i = 0; i < possible_values.length; i++) {
				const variable_assignment& values = possible_values[i];
				variable_assignment& new_values = temp_possible_values[temp_possible_values.length];
				if (!::init(new_values, values)) {
					quantifiers.length = old_quantifier_length;
					core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
					return false;
				}

				bool unifies = true;
				for (const auto& unification : first_unifications) {
					if (unification.key->quantifier.variable <= old_quantifier_length
					 && !new_values.unify_value(unification.key->quantifier.variable - 1, unification.value))
					{
						unifies = false;
						break;
					}
				}
				if (!unifies) {
					core::free(new_values);
					continue;
				}
				temp_possible_values.length++;
			}

			if (temp_possible_values.length != 0) {
				array<variable_assignment> temp(new_possible_values.length + temp_possible_values.length);
				set_union(temp, new_possible_values, temp_possible_values);
				swap(temp, new_possible_values);
				for (auto& element : temp_possible_values) core::free(element);
				for (auto& element : temp) core::free(element);
			}
		}
		quantifiers.length = old_quantifier_length;
		core::free(*new_operand); if (new_operand->reference_count == 0) core::free(new_operand);
		return true;
	}

	template<bool Contradiction, typename Prover>
	bool exists_is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			const array<variable_assignment>& possible_values,
			array<variable_assignment>& new_possible_values,
			Prover& prover) const
	{
		Formula* operand = formula;
		unsigned int old_quantifier_length = quantifiers.length;
		uint_fast8_t old_length = possible_values[0].assignment.length;
		array<unsigned int> quantified_variables(4);
		do {
			if (!quantifiers.ensure_capacity(operand->quantifier.variable)
			 || !quantified_variables.add(operand->quantifier.variable))
				return false;
			while (operand->quantifier.variable - 1 >= quantifiers.length)
				quantifiers[quantifiers.length++] = nullptr;
			quantifiers[operand->quantifier.variable - 1] = operand;
			operand = operand->quantifier.operand;
		} while (operand->type == FormulaType::EXISTS);

		array<variable_assignment> temp_possible_values(possible_values.length);
		for (const variable_assignment& assignment : possible_values) {
			variable_assignment& new_value = temp_possible_values[temp_possible_values.length];
			if (!::init(new_value, assignment, quantifiers.length)) {
				quantifiers.length = old_quantifier_length;
				for (auto& element : temp_possible_values) core::free(element);
				return false;
			}
			temp_possible_values.length++;
		}

		if (!is_provable_without_abduction<Contradiction>(operand, quantifiers, temp_possible_values, prover)) {
			quantifiers.length = old_quantifier_length;
			for (auto& element : temp_possible_values) core::free(element);
			return true;
		}
		quantifiers.length = old_quantifier_length;

		unsigned int next_position = 0;
		for (unsigned int i = 0; i < temp_possible_values.length; i++) {
			variable_assignment& assignment = temp_possible_values[i];
			assignment.assignment.length = old_length;
			for (unsigned int j = 0; j < assignment.assignment.equal_indices.length; j++) {
				if (assignment.assignment.equal_indices[j].key >= old_length
				 || assignment.assignment.equal_indices[j].value >= old_length)
					assignment.assignment.equal_indices.remove(j--);
			} for (unsigned int j = 0; j < assignment.assignment.not_equal_indices.length; j++) {
				if (assignment.assignment.not_equal_indices[j].key >= old_length
				 || assignment.assignment.not_equal_indices[j].value >= old_length)
					assignment.assignment.not_equal_indices.remove(j--);
			} for (unsigned int j = 0; j < assignment.assignment.ge_indices.length; j++) {
				if (assignment.assignment.ge_indices[j].key >= old_length
				 || assignment.assignment.ge_indices[j].value >= old_length)
					assignment.assignment.ge_indices.remove(j--);
			}

			/* check if this new assignment is already a subset of an existing assignment */
			bool found_superset = false;
			for (unsigned int j = 0; j < next_position; j++) {
				if (is_subset(assignment, temp_possible_values[j])) {
					found_superset = true;
					break;
				}
			}
			if (found_superset) {
				core::free(assignment);
				continue;
			}

			/* check if an existing assignment is a subset of this new assignment */
			for (unsigned int j = 0; j < next_position; j++) {
				if (is_subset(temp_possible_values[j], assignment)) {
					core::free(temp_possible_values[j]);
					core::move(temp_possible_values[next_position - 1], temp_possible_values[j]);
					next_position--; j--;
				}
			}

			core::move(assignment, temp_possible_values[next_position++]);
		}
		temp_possible_values.length = next_position;

		array<variable_assignment> union_result(new_possible_values.length + temp_possible_values.length);
		set_union(union_result, new_possible_values, temp_possible_values);
		swap(union_result, new_possible_values);
		for (auto& element : union_result) core::free(element);
		for (auto& element : temp_possible_values) core::free(element);
		return true;

		/*Formula* quantified = formula->quantifier.operand;
		unsigned int variable = formula->quantifier.variable;

		array<instance> constants(ground_concept_capacity + 1);
		array<hol_number> numbers(64); array<string*> strings(64);
		for (unsigned int i = 0; i < ground_concept_capacity; i++) {
			if (ground_concepts[i].types.keys != nullptr) {
				constants[constants.length].type = instance_type::CONSTANT;
				constants[constants.length].matching_types = 0;
				constants[constants.length].mismatching_types = 0;
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
			if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i)) continue;
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
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].number = number;
		} for (string* str : strings) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].str = str;
		}

		if (!filter_constants_helper(*this, quantified, variable, constants))
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

			array<variable_assignment> copy(possible_values.length);
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

			array<variable_assignment> union_result(new_possible_values.length + copy.length);
			set_union(union_result, new_possible_values, copy);
			swap(union_result, new_possible_values);
			for (auto& element : union_result) core::free(element);
			for (auto& element : copy) core::free(element);
		}
		core::free(*var); if (var->reference_count == 0) core::free(var);
		return true;*/
	}

	inline bool is_name_scope(const Formula* formula, const Term*& arg1, const Term*& arg2) const {
		if (formula->type != FormulaType::EXISTS)
			return false;
		unsigned int variable = formula->quantifier.variable;
		const Formula* operand = formula->quantifier.operand;
		if (operand->type != FormulaType::AND)
			return false;
		bool found_name = false;
		arg1 = nullptr;
		arg2 = nullptr;
		for (unsigned int i = 0; i < operand->array.length; i++) {
			const Formula* conjunct = operand->array.operands[i];
			if (conjunct->type == FormulaType::UNARY_APPLICATION
			 && conjunct->binary.left->type == TermType::CONSTANT
			 && conjunct->binary.left->constant == (unsigned int) built_in_predicates::NAME
			 && conjunct->binary.right->type == TermType::VARIABLE
			 && conjunct->binary.right->variable == variable)
			{
				found_name = true;
			} else if (conjunct->type == FormulaType::EQUALS
					&& conjunct->binary.left->type == TermType::UNARY_APPLICATION
					&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
					&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
					&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
					&& conjunct->binary.left->binary.right->variable == variable)
			{
				arg1 = conjunct->binary.right;
			} else if (conjunct->type == FormulaType::EQUALS
					&& conjunct->binary.left->type == TermType::UNARY_APPLICATION
					&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
					&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
					&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
					&& conjunct->binary.left->binary.right->variable == variable
					&& conjunct->binary.right->type == TermType::STRING)
			{
				arg2 = conjunct->binary.right;
			}
		}
		return (found_name && arg1 != nullptr && arg2 != nullptr);
	}

	template<bool Contradiction, typename Prover>
	bool is_provable_without_abduction(
			Formula* formula, array<Formula*>& quantifiers,
			array<variable_assignment>& possible_values,
			Prover& prover) const
	{
		/* consider if this formula provable by other extensional edges (i.e. universally-quantified theorems) */
		array<variable_assignment> new_possible_values(possible_values.length);
		if (!is_provable_by_theorem_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values, prover)) {
			for (auto& element : new_possible_values) core::free(element);
			return false;
		} if (!is_provable_by_exclusion_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values, prover)) {
			for (auto& element : new_possible_values) core::free(element);
			return false;
		}
		if (!prover.template is_provable_by_axiom_without_abduction<Contradiction>(formula, quantifiers, possible_values, new_possible_values)) {
			for (auto& element : new_possible_values) core::free(element);
			for (variable_assignment& assignment : possible_values)
				if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
			return false;
		}

		if (formula->type == FormulaType::UNARY_APPLICATION) {
			if (formula->binary.left->type == TermType::CONSTANT && formula->binary.left->constant == (unsigned int) built_in_predicates::NUMBER) {
				if (Contradiction) {
					if (formula->binary.right->type == TermType::CONSTANT) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return true;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						array<variable_assignment> temp(possible_values.length);
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
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
						for (auto& element : possible_values) core::free(element);
						swap(new_possible_values, possible_values);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return (possible_values.length > 0);
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				} else {
					if (formula->binary.right->type == TermType::NUMBER) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return true;
					} else if (formula->binary.right->type == TermType::VARIABLE) {
						array<variable_assignment> temp(possible_values.length);
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
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
						for (auto& element : possible_values) core::free(element);
						swap(new_possible_values, possible_values);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return (possible_values.length > 0);
					} else {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : possible_values) core::free(element);
						possible_values.clear(); return false;
					}
				}
			}

			if (formula->binary.right->type != TermType::CONSTANT
			 && formula->binary.right->type != TermType::NUMBER
			 && formula->binary.right->type != TermType::VARIABLE)
			{
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				possible_values.clear(); return false;
			}

			if (!has_free_variables(*formula->binary.left)) {
				/* optimization: if the operator has no free variables, we can simply perform a hash lookup in `atoms` instead of looping over them */
				Term* atom_formula = Term::new_apply(formula->binary.left, &Variables<1>::value);
				if (atom_formula == nullptr) {
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : possible_values) core::free(element);
					possible_values.clear(); return false;
				}
				formula->binary.left->reference_count++;
				Variables<1>::value.reference_count++;

				bool contains;
				const pair<array<instance>, array<instance>>& value = atoms.get(*atom_formula, contains);
				core::free(*atom_formula); core::free(atom_formula);
				if (contains) {
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const array<instance>& constants = (Contradiction ? value.value : value.key);
						array<variable_assignment> temp(possible_values.length + constants.length);
						const variable_assignment& values = possible_values[i];
						if (formula->binary.right->type == TermType::VARIABLE) {
							for (const instance& constant : constants) {
								variable_assignment& current_new_values = temp[temp.length];
								if (!::init(current_new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								Term* term = nullptr;
								if (constant.type == instance_type::CONSTANT)
									term = Term::new_constant(constant.constant);
								else if (constant.type == instance_type::NUMBER)
									term = Term::new_number(constant.number);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								bool unifies = current_new_values.unify_value(formula->binary.right->variable - 1, term);
								core::free(*term); core::free(term);
								if (!unifies) {
									core::free(current_new_values);
									continue;
								}
								temp.length++;
							}
							if (temp.length > 0)
								insertion_sort(temp, default_sorter());
						} else if ((formula->binary.right->type == TermType::CONSTANT && index_of_constant(constants, formula->binary.right->constant) < constants.length)
								|| (formula->binary.right->type == TermType::NUMBER && index_of_number(constants, formula->binary.right->number) < constants.length)) {
							if (!::init(temp[0], values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
								return false;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
					}
				}

			} else {
				for (const auto& atom : atoms) {
					array_map<Formula*, Term*> unifications(2);
					if (!unify_atom(formula->binary.left, atom.key.binary.left, quantifiers, unifications))
						continue;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const array<instance>& constants = (Contradiction ? atom.value.value : atom.value.key);
						array<variable_assignment> temp(possible_values.length + constants.length);
						const variable_assignment& values = possible_values[i];
						variable_assignment new_values(values);
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
							for (const instance& constant : constants) {
								variable_assignment& current_new_values = temp[temp.length];
								if (!::init(current_new_values, new_values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								Term* term = nullptr;
								if (constant.type == instance_type::CONSTANT)
									term = Term::new_constant(constant.constant);
								else if (constant.type == instance_type::NUMBER)
									term = Term::new_number(constant.number);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
						} else if ((formula->binary.right->type == TermType::CONSTANT && index_of_constant(constants, formula->binary.right->constant) < constants.length)
								|| (formula->binary.right->type == TermType::NUMBER && index_of_number(constants, formula->binary.right->number) < constants.length)) {
							if (!::init(temp[0], new_values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
								return false;
							}
							temp.length++;
						}
						if (temp.length != 0) {
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
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
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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

					typename Prover::ProvableElementArray provable_elements = prover.get_provable_elements(set_id);
					if (Contradiction && sets.sets[set_id].arity != args.length) {
						/* since the arity of the elements don't match, this formula is trivially disprovable */
						prover.free_provable_elements(provable_elements);
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return true;
					}
					if (Contradiction) {
						/* we can only prove a contradiction by exclusion, which is
						   already handled by `is_provable_by_exclusion_without_abduction` */
						prover.free_provable_elements(provable_elements);
						continue;
					}

					for (unsigned int j = 0; j < possible_values.length; j++) {
						variable_assignment new_values(possible_values[j]);
						if (left_most->type == TermType::VARIABLE && !new_values.unify_value(left_most->variable - 1, set_formula->binary.left))
							continue;

						array<variable_assignment> temp(max(1u, provable_elements.length));
						for (const tuple& provable_element : provable_elements) {
							variable_assignment& new_new_values = temp[temp.length];
							if (!::init(new_new_values, new_values)) {
								prover.free_provable_elements(provable_elements);
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
									if (!new_new_values.unify_value(args[k]->variable - 1, constant)) {
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

						if (temp.length != 0) {
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							swap(new_possible_values, union_result);
							for (auto& element : temp) core::free(element);
							for (auto& element : union_result) core::free(element);
						}
					}
					prover.free_provable_elements(provable_elements);
				}

				for (unsigned int k = 0; k < args.length; k++) {
					if (args[k]->type == TermType::VARIABLE) {
						if (Contradiction) {
							array<variable_assignment> temp(max((size_t) 1, possible_values.length));
							for (const variable_assignment& value : possible_values) {
								variable_assignment& new_value = temp[temp.length];
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
								array<variable_assignment> union_result(new_possible_values.length + temp.length);
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
							array<variable_assignment> temp(max((size_t) 1, possible_values.length));
							for (const variable_assignment& value : possible_values) {
								variable_assignment& new_value = temp[temp.length];
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
								array<variable_assignment> union_result(new_possible_values.length + temp.length);
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
				array<variable_assignment> temp(possible_values.length);
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
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return true;
					} else if (Contradiction && formula->ternary.second->number < formula->ternary.third->number) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
					sort(temp, default_sorter());
					array<variable_assignment> union_result(new_possible_values.length + temp.length);
					set_union(union_result, new_possible_values, temp);
					for (auto& element : temp) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : possible_values) core::free(element);
					swap(union_result, possible_values);
				} else {
					for (auto& element : possible_values) core::free(element);
					swap(new_possible_values, possible_values);
				}
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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

			array<variable_assignment> temp(possible_values.length);
			for (const auto& rel : relations) {
				if (rel.key.predicate != 0 && predicate != 0 && rel.key.predicate != predicate) continue;
				if (rel.key.arg1 != 0 && arg1 != 0 && rel.key.arg1 != arg1) continue;
				if (rel.key.arg2 != 0 && arg2 != 0 && rel.key.arg2 != arg2) continue;

				const array<instance>& constants = (Contradiction ? rel.value.value : rel.value.key);
				for (unsigned int i = 0; i < possible_values.length; i++) {
					if (!temp.ensure_capacity(temp.length + constants.length + 1)) {
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return false;
					}
					variable_assignment values(possible_values[i]);
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
							for (const instance& constant : constants) {
								variable_assignment& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								Term* term = Term::new_constant(constant.constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							variable_assignment& new_values = temp[temp.length];
							if (index_of_constant(constants, predicate) < constants.length) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								temp.length++;
							}
						}
					} else if (rel.key.arg1 == 0) {
						if (arg1 == 0) {
							for (const instance& constant : constants) {
								variable_assignment& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								Term* term = Term::new_constant(constant.constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							variable_assignment& new_values = temp[temp.length];
							if (index_of_constant(constants, arg1) < constants.length) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								temp.length++;
							}
						}
					} else {
						if (arg2 == 0) {
							for (const instance& constant : constants) {
								variable_assignment& new_values = temp[temp.length];
								if (!::init(new_values, values)) {
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								Term* term = Term::new_constant(constant.constant);
								if (term == nullptr) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							variable_assignment& new_values = temp[temp.length];
							if (index_of_constant(constants, arg2) < constants.length) {
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (auto& element : temp) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								temp.length++;
							}
						}
					}
				}
			}
			if (temp.length != 0) {
				sort(temp, default_sorter()); unique_and_cleanup(temp);
				array<variable_assignment> union_result(new_possible_values.length + temp.length);
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
					array<variable_assignment> copy(possible_values.length);
					for (unsigned int j = 0; j < possible_values.length; j++) {
						if (!::init(copy[copy.length], possible_values[j])) {
							for (auto& element : copy) core::free(element);
							for (auto& element : new_possible_values) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						copy.length++;
					}
					if (!is_provable_without_abduction<true>(formula->array.operands[i], quantifiers, copy, prover))
						continue;

					array<variable_assignment> union_result(new_possible_values.length + copy.length);
					set_union(union_result, new_possible_values, copy);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : copy) core::free(element);
					swap(union_result, new_possible_values);
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			/* optimization: to more quickly narrow the set of `possible_values`, look for conjuncts of the form `argn(x)=y` where x is a variable with only constant values in `possible_values` */
			bool* visited = (bool*) alloca(sizeof(bool) * formula->array.length);
			for (unsigned int i = 0; i < formula->array.length; i++) visited[i] = false;

			for (unsigned int i = 0; i < formula->array.length; i++) {
				const Term* arg1; const Term* arg2;
				if (is_name_scope(formula->array.operands[i], arg1, arg2)) {
					visited[i] = true;
					if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
				} else if (formula->array.operands[i]->type == TermType::EQUALS
						&& formula->array.operands[i]->binary.left->type == TermType::UNARY_APPLICATION
						&& formula->array.operands[i]->binary.left->binary.left->type == TermType::CONSTANT
						&& formula->array.operands[i]->binary.left->binary.right->type == TermType::VARIABLE)
				{
					if (formula->array.operands[i]->binary.right->type == TermType::STRING
					 || formula->array.operands[i]->binary.right->type == TermType::NUMBER
					 || (formula->array.operands[i]->binary.right->type == TermType::CONSTANT
					  && (new_constant_offset > formula->array.operands[i]->binary.right->constant
					   || 16 > ground_concepts[formula->array.operands[i]->binary.right->constant - new_constant_offset].definitions.length)))
					{
						visited[i] = true;
						if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
					} else {
						bool constant_valued_variable = true;
						for (const variable_assignment& value : possible_values) {
							if (value.assignment.values[formula->array.operands[i]->binary.left->binary.right->variable - 1].type != instantiation_type::CONSTANT) {
								constant_valued_variable = false;
								break;
							}
						}
						if (constant_valued_variable) {
							visited[i] = true;
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						}
					}
				}
			}

			if (possible_values.length != 0) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					if (!visited[i] && formula->array.operands[i]->type == TermType::EQUALS
					 && formula->array.operands[i]->binary.right->type == TermType::VARIABLE)
					{
						bool constant_valued_variable = true;
						for (const variable_assignment& value : possible_values) {
							if (value.assignment.values[formula->array.operands[i]->binary.right->variable - 1].type != instantiation_type::CONSTANT) {
								constant_valued_variable = false;
								break;
							}
						}
						if (constant_valued_variable) {
							visited[i] = true;
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						}
					}
				}
			}

			/* optimization: to more quickly narrow the set of `possible_values`, look for conjuncts of the form `c(x)` where `c` is a constant */
			if (possible_values.length != 0) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					if (!visited[i] && formula->array.operands[i]->type == TermType::UNARY_APPLICATION && !(formula->array.operands[i]->binary.left->type == TermType::VARIABLE && formula->array.operands[i]->binary.right->type == TermType::VARIABLE)) {
						visited[i] = true;
						if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
					}
				}
			}

			if (possible_values.length != 0) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					if (!visited[i] && formula->array.operands[i]->type != TermType::UNARY_APPLICATION) {
						visited[i] = true;
						if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
					}
				}
			}

			if (possible_values.length != 0) {
				for (unsigned int i = 0; i < formula->array.length; i++) {
					if (!visited[i]) {
						if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
					}
				}
			}

			if (possible_values.length != 0) {
				array<variable_assignment> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) core::free(element);
				for (auto& element : new_possible_values) core::free(element);
				swap(possible_values, union_result);
			} else {
				swap(possible_values, new_possible_values);
			}

		} else if (formula->type == FormulaType::OR) {
			if (Contradiction) {
				/* optimization: to more quickly narrow the set of `possible_values`, look for conjuncts of the form `argn(x)=y` where x is a variable with only constant values in `possible_values` */
				bool* visited = (bool*) alloca(sizeof(bool) * formula->array.length);
				for (unsigned int i = 0; i < formula->array.length; i++) visited[i] = false;

				for (unsigned int i = 0; i < formula->array.length; i++) {
					const Term* arg1; const Term* arg2;
					if (is_name_scope(formula->array.operands[i], arg1, arg2)) {
						visited[i] = true;
						if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
					} else if (formula->array.operands[i]->type == TermType::EQUALS
							&& formula->array.operands[i]->binary.left->type == TermType::UNARY_APPLICATION
							&& formula->array.operands[i]->binary.left->binary.left->type == TermType::CONSTANT
							&& formula->array.operands[i]->binary.left->binary.right->type == TermType::VARIABLE)
					{
						if (formula->array.operands[i]->binary.right->type == TermType::STRING
						 || formula->array.operands[i]->binary.right->type == TermType::NUMBER
						 || (formula->array.operands[i]->binary.right->type == TermType::CONSTANT
						  && (new_constant_offset > formula->array.operands[i]->binary.right->constant
						   || 16 > ground_concepts[formula->array.operands[i]->binary.right->constant - new_constant_offset].definitions.length)))
						{
							visited[i] = true;
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						} else {
							bool constant_valued_variable = true;
							for (const variable_assignment& value : possible_values) {
								if (value.assignment.values[formula->array.operands[i]->binary.left->binary.right->variable - 1].type != instantiation_type::CONSTANT) {
									constant_valued_variable = false;
									break;
								}
							}
							if (constant_valued_variable) {
								visited[i] = true;
								if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
							}
						}
					}
				}

				if (possible_values.length != 0) {
					for (unsigned int i = 0; i < formula->array.length; i++) {
						if (!visited[i] && formula->array.operands[i]->type == TermType::EQUALS
						 && formula->array.operands[i]->binary.right->type == TermType::VARIABLE)
						{
							bool constant_valued_variable = true;
							for (const variable_assignment& value : possible_values) {
								if (value.assignment.values[formula->array.operands[i]->binary.right->variable - 1].type != instantiation_type::CONSTANT) {
									constant_valued_variable = false;
									break;
								}
							}
							if (constant_valued_variable) {
								visited[i] = true;
								if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
							}
						}
					}
				}

				/* optimization: to more quickly narrow the set of `possible_values`, look for conjuncts of the form `c(x)` where `c` is a constant */
				if (possible_values.length != 0) {
					for (unsigned int i = 0; i < formula->array.length; i++) {
						if (!visited[i] && formula->array.operands[i]->type == TermType::UNARY_APPLICATION && !(formula->array.operands[i]->binary.left->type == TermType::VARIABLE && formula->array.operands[i]->binary.right->type == TermType::VARIABLE)) {
							visited[i] = true;
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						}
					}
				}

				if (possible_values.length != 0) {
					for (unsigned int i = 0; i < formula->array.length; i++) {
						if (!visited[i] && formula->array.operands[i]->type != TermType::UNARY_APPLICATION) {
							visited[i] = true;
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						}
					}
				}

				if (possible_values.length != 0) {
					for (unsigned int i = 0; i < formula->array.length; i++) {
						if (!visited[i]) {
							if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, possible_values, prover)) break;
						}
					}
				}

				if (possible_values.length != 0) {
					array<variable_assignment> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					swap(possible_values, union_result);
				} else {
					swap(possible_values, new_possible_values);
				}
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			for (unsigned int i = 0; i < formula->array.length; i++) {
				array<variable_assignment> copy(possible_values.length);
				for (unsigned int j = 0; j < possible_values.length; j++) {
					if (!::init(copy[copy.length], possible_values[j])) {
						for (auto& element : copy) core::free(element);
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return false;
					}
					copy.length++;
				}
				if (!is_provable_without_abduction<false>(formula->array.operands[i], quantifiers, copy, prover))
					continue;

				array<variable_assignment> union_result(new_possible_values.length + copy.length);
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
					array<variable_assignment> union_result(possible_values.length + new_possible_values.length);
					set_union(union_result, possible_values, new_possible_values);
					for (auto& element : possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					swap(union_result, possible_values);
				} else {
					swap(possible_values, new_possible_values);
				}
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			array<variable_assignment> first(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(first[first.length], possible_values[j])) {
					for (auto& element : first) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}
				first.length++;
			}
			if (!is_provable_without_abduction<true>(formula->binary.left, quantifiers, first, prover))
				first.length = 0;

			array<variable_assignment> second(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(second[second.length], possible_values[j])) {
					for (auto& element : first) core::free(element);
					for (auto& element : second) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}
				second.length++;
			}
			if (!is_provable_without_abduction<false>(formula->binary.left, quantifiers, second, prover))
				second.length = 0;

			array<variable_assignment> temp(max((size_t) 1, new_possible_values.length + first.length));
			set_union(temp, new_possible_values, first);
			for (auto& element : first) core::free(element);
			for (auto& element : new_possible_values) core::free(element);
			new_possible_values.length = 0;
			if (!new_possible_values.ensure_capacity(temp.length + second.length)) {
				for (auto& element : temp) core::free(element);
				for (auto& element : second) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return false;
			}
			if (possible_values.length != 0) {
				array<variable_assignment> union_result(possible_values.length + new_possible_values.length);
				set_union(union_result, possible_values, new_possible_values);
				for (auto& element : possible_values) core::free(element);
				for (auto& element : new_possible_values) core::free(element);
				swap(union_result, possible_values);
			} else {
				swap(possible_values, new_possible_values);
			}
			for (variable_assignment& assignment : possible_values)
				if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
			return possible_values.length != 0;

		} else if (formula->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				if (!exists_is_provable_without_abduction<true>(formula, quantifiers, possible_values, new_possible_values, prover)) {
					for (auto& element : new_possible_values) core::free(element);
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			if (!for_all_is_provable_without_abduction(formula, quantifiers, possible_values, new_possible_values, prover)) {
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return false;
			}
			swap(possible_values, new_possible_values);
			for (auto& element : new_possible_values) core::free(element);

		} else if (formula->type == FormulaType::EXISTS) {
			if (Contradiction) {
				if (!not_exists_is_provable_without_abduction(formula, quantifiers, possible_values, new_possible_values, prover)) {
					for (auto& element : new_possible_values) core::free(element);
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}

				/* TODO: implement this */
				//fprintf(stderr, "theory.is_provable_without_abduction ERROR: Not implemented.\n");
				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			if (!exists_is_provable_without_abduction<false>(formula, quantifiers, possible_values, new_possible_values, prover)) {
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
						if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i)) continue;

						if (right->type == TermType::NUMBER) {
							if ((!Contradiction && right->number != sets.sets[i].set_size)
							 || (Contradiction && right->number == sets.sets[i].set_size))
								continue;
						}

						array<Formula*> quantifiers(4);
						array_map<Formula*, Term*> unifications(4);
						if (!unify(set_definition, sets.sets[i].size_axioms[0]->formula->binary.left->binary.right, quantifiers, unifications))
							continue;

						array<variable_assignment> temp(possible_values.length);
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const variable_assignment& values = possible_values[i];
							variable_assignment& new_values = temp[temp.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (set_definition->type == TermType::VARIABLE) {
					/* size(x)=n */
					for (unsigned int i = 1; i < sets.set_count + 1; i++) {
						if (sets.sets[i].size_axioms.data == nullptr || prover.is_set_removed(i)) continue;
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

						array<variable_assignment> temp(possible_values.length);
						for (unsigned int j = 0; j < possible_values.length; j++) {
							const variable_assignment& values = possible_values[j];
							variable_assignment& new_values = temp[temp.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (auto& element : temp) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length + temp.length);
							set_union(union_result, new_possible_values, temp);
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (set_definition->type == TermType::CONSTANT && set_definition->constant >= new_constant_offset) {
					/* size(a)=n */
					Term* set_formula = Term::new_apply(set_definition, &Variables<1>::value);
					if (set_formula == nullptr) {
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return false;
					}
					set_definition->reference_count++;
					Variables<1>::value.reference_count++;
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
					array<variable_assignment> temp(possible_values.length);
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const variable_assignment& values = possible_values[i];
						variable_assignment& new_values = temp[temp.length];
						if (!::init(new_values, values)) {
							for (auto& element : new_possible_values) core::free(element);
							for (auto& element : temp) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						if (right->type == TermType::VARIABLE) {
							if ((!Contradiction && !new_values.unify_value(right->variable - 1, sets.sets[set_id].size_axioms[0]->formula->binary.right))
							 || (Contradiction && !new_values.antiunify_value(right->variable - 1, sets.sets[set_id].size_axioms[0]->formula->binary.right)))
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
						array<variable_assignment> union_result(new_possible_values.length + temp.length);
						set_union(union_result, new_possible_values, temp);
						for (auto& element : new_possible_values) core::free(element);
						for (auto& element : temp) core::free(element);
						swap(union_result, new_possible_values);
					}

				} else {
					fprintf(stderr, "theory.is_provable_without_abduction ERROR: Unsupported set size axiom.\n");
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}

				swap(possible_values, new_possible_values);
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return possible_values.length != 0;
			}

			/* A=A */
			array_map<Formula*, Term*> unifications(4);
			if (!Contradiction && unify(formula->binary.left, formula->binary.right, quantifiers, unifications)) {
				if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
					for (auto& element : new_possible_values) core::free(element);
					for (variable_assignment& assignment : possible_values)
						if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
					return false;
				}
				unsigned int old_size = new_possible_values.length;
				for (unsigned int i = 0; i < possible_values.length; i++) {
					const variable_assignment& values = possible_values[i];
					variable_assignment& new_values = new_possible_values[new_possible_values.length];
					if (!::init(new_values, values)) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
					array<variable_assignment> union_result(new_possible_values.length);
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
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const variable_assignment& values = possible_values[i];
							variable_assignment& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				} else if (formula->binary.left->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return false;
					}
					unsigned int old_size = new_possible_values.length;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const variable_assignment& values = possible_values[i];
						variable_assignment& new_values = new_possible_values[new_possible_values.length];
						if (!::init(new_values, values)) {
							for (auto& element : new_possible_values) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
						array<variable_assignment> union_result(new_possible_values.length);
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
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const variable_assignment& values = possible_values[i];
								variable_assignment& new_values = new_possible_values[new_possible_values.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
								array<variable_assignment> union_result(new_possible_values.length);
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
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (unsigned int i = 0; i < possible_values.length; i++) {
							const variable_assignment& values = possible_values[i];
							variable_assignment& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length);
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
					const variable_assignment& values = possible_values[i];
					if (!Contradiction
					 && values.assignment.values[left->variable - 1].type == instantiation_type::CONSTANT
					 && values.assignment.values[left->variable - 1].constant >= new_constant_offset)
					{
						const concept<ProofCalculus>& c = ground_concepts[values.assignment.values[left->variable - 1].constant - new_constant_offset];
						if (!new_possible_values.ensure_capacity(new_possible_values.length + c.definitions.length)) {
							for (auto& element : new_possible_values) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						unsigned int old_size = new_possible_values.length;
						for (Proof* definition : c.definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;

							variable_assignment& new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_values, values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
							swap(union_result, new_possible_values);
						}
						continue;
					} else if (!Contradiction && values.assignment.values[left->variable - 1].type != instantiation_type::ANY) {
						continue;
					} else if (right->type == TermType::UNARY_APPLICATION
							&& right->binary.left->type == TermType::CONSTANT
							&& right->binary.right->type == TermType::VARIABLE
							&& values.assignment.values[right->binary.right->variable - 1].type != instantiation_type::ANY
							&& values.assignment.values[right->binary.right->variable - 1].type != instantiation_type::ANY_NUMBER)
					{
						Term* right_term = nullptr;
						if (values.assignment.values[right->binary.right->variable - 1].type == instantiation_type::CONSTANT)
							right_term = Term::new_apply(right->binary.left, Term::new_constant(values.assignment.values[right->binary.right->variable - 1].constant));
						else if (values.assignment.values[right->binary.right->variable - 1].type == instantiation_type::NUMBER)
							right_term = Term::new_apply(right->binary.left, Term::new_number(values.assignment.values[right->binary.right->variable - 1].number));
						else if (values.assignment.values[right->binary.right->variable - 1].type == instantiation_type::STRING)
							right_term = Term::new_apply(right->binary.left, Term::new_string(values.assignment.values[right->binary.right->variable - 1].str));
						if (right_term == nullptr) {
							for (auto& element : new_possible_values) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						right->binary.left->reference_count++;

						bool contains;
						unsigned int constant = reverse_definitions.get(*right_term, contains);
						core::free(*right_term); core::free(right_term);
						if (contains) {
							if (!new_possible_values.ensure_capacity(new_possible_values.length + 1)
							 || !::init(new_possible_values[new_possible_values.length], values))
							{
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
								return false;
							}
							variable_assignment& new_values = new_possible_values[new_possible_values.length];
							if ((!Contradiction && new_values.unify_value(left->variable - 1, ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.left))
							 || (Contradiction && new_values.antiunify_value(left->variable - 1, ground_concepts[constant - new_constant_offset].definitions[0]->formula->binary.left)))
							{
								new_possible_values.length++;
								array<variable_assignment> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, new_possible_values.length - 1, new_possible_values.data + new_possible_values.length - 1, 1);
								for (auto& element : new_possible_values) core::free(element);
								swap(union_result, new_possible_values);
							} else {
								core::free(new_values);
							}
						}
						continue;
					}

					for (unsigned int i = 0; i < ground_concept_capacity; i++) {
						const concept<ProofCalculus>& c = ground_concepts[i];
						if (c.types.keys == nullptr)
							continue;
						unsigned int old_size = new_possible_values.length;
						if (!new_possible_values.ensure_capacity(new_possible_values.length + c.definitions.length)) {
							for (auto& element : new_possible_values) core::free(element);
							for (variable_assignment& assignment : possible_values)
								if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
							return false;
						}
						variable_assignment new_values(values);
						if ((!Contradiction && !new_values.unify_value(left->variable - 1, c.definitions[0]->formula->binary.left))
						 || (Contradiction && !new_values.antiunify_value(left->variable - 1, c.definitions[0]->formula->binary.left)))
							continue;

						for (Proof* definition : c.definitions) {
							array_map<Formula*, Term*> unifications(4);
							if (!unify(right, definition->formula->binary.right, quantifiers, unifications))
								continue;

							variable_assignment& new_new_values = new_possible_values[new_possible_values.length];
							if (!::init(new_new_values, new_values)) {
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
							array<variable_assignment> union_result(new_possible_values.length);
							set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
							for (auto& element : new_possible_values) core::free(element);
							swap(union_result, new_possible_values);
						}
					}
				}
			}

			/* a(b)=A or A=a(b) (function definition of a evaluated at b) */
			if ((left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT && right->type != TermType::CONSTANT)
			 || (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT && left->type != TermType::CONSTANT))
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
								for (auto& element : new_possible_values) core::free(element);
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
								for (variable_assignment& assignment : possible_values)
									if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
								return false;
							}
							unsigned int old_size = new_possible_values.length;
							for (unsigned int i = 0; i < possible_values.length; i++) {
								const variable_assignment& values = possible_values[i];
								variable_assignment& new_values = new_possible_values[new_possible_values.length];
								if (!::init(new_values, values)) {
									for (auto& element : new_possible_values) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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
								array<variable_assignment> union_result(new_possible_values.length);
								set_union(union_result.data, union_result.length, new_possible_values.data, old_size, new_possible_values.data + old_size, new_possible_values.length - old_size);
								for (auto& element : new_possible_values) core::free(element);
								swap(union_result, new_possible_values);
							}
						}
					}
				} else if (left->binary.right->type == TermType::VARIABLE) {
					if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
						for (auto& element : new_possible_values) core::free(element);
						for (variable_assignment& assignment : possible_values)
							if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
						return false;
					}
					unsigned int old_size = new_possible_values.length;
					for (unsigned int i = 0; i < possible_values.length; i++) {
						const variable_assignment& values = possible_values[i];
						if (right->type == TermType::VARIABLE && values.assignment.values[right->variable - 1].type == instantiation_type::CONSTANT)
						{
							/* if `right` is a constant, this is a definition and not a function value */
							continue;
						} else if (values.assignment.values[left->binary.right->variable - 1].type == instantiation_type::CONSTANT
								&& values.assignment.values[left->binary.right->variable - 1].constant >= new_constant_offset)
						{
							const concept<ProofCalculus>& c = ground_concepts[values.assignment.values[left->binary.right->variable - 1].constant - new_constant_offset];
							bool contains;
							Proof* function_value = c.function_values.get(left->binary.left->constant, contains);
							if (contains) {
								array_map<Formula*, Term*> unifications(4);
								if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
									if (Contradiction) {
										variable_assignment& new_values = new_possible_values[new_possible_values.length];
										if (!::init(new_values, values)) {
											for (auto& element : new_possible_values) core::free(element);
											for (variable_assignment& assignment : possible_values)
												if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									variable_assignment& new_values = new_possible_values[new_possible_values.length];
									if (!::init(new_values, values)) {
										for (auto& element : new_possible_values) core::free(element);
										for (variable_assignment& assignment : possible_values)
											if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
										return false;
									}

									/* make sure the `possible_values` unify with the variables in this term */
									bool unifies = true;
									for (const auto& unification : unifications) {
										/* TODO: i dont think this is exactly correct when Contradiction is true */
										if ((!Contradiction && !new_values.unify_value(unification.key->quantifier.variable - 1, unification.value))
										 || (Contradiction && !new_values.antiunify_value(unification.key->quantifier.variable - 1, unification.value)))
										{
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
						} else if (!Contradiction && values.assignment.values[left->binary.right->variable - 1].type != instantiation_type::ANY) {
							continue;
						}

						for (unsigned int i = 0; i < ground_concept_capacity; i++) {
							const concept<ProofCalculus>& c = ground_concepts[i];
							if (c.types.keys == nullptr)
								continue;
							variable_assignment new_values(values);
							if (!new_values.unify_value(left->binary.right->variable - 1, c.definitions[0]->formula->binary.left))
								continue;

							bool contains;
							Proof* function_value = c.function_values.get(left->binary.left->constant, contains);
							if (contains) {
								array_map<Formula*, Term*> unifications(4);
								if (!new_possible_values.ensure_capacity(new_possible_values.length + possible_values.length)) {
									for (auto& element : new_possible_values) core::free(element);
									for (variable_assignment& assignment : possible_values)
										if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
									return false;
								}
								if (!unify(right, function_value->formula->binary.right, quantifiers, unifications)) {
									if (Contradiction) {
										variable_assignment& new_new_values = new_possible_values[new_possible_values.length];
										if (!::init(new_new_values, new_values)) {
											for (auto& element : new_possible_values) core::free(element);
											for (variable_assignment& assignment : possible_values)
												if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
											return false;
										}
										new_possible_values.length++;
									}
								} else {
									variable_assignment& new_new_values = new_possible_values[new_possible_values.length];
									if (!::init(new_new_values, new_values)) {
										for (auto& element : new_possible_values) core::free(element);
										for (variable_assignment& assignment : possible_values)
											if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
										return false;
									}

									/* make sure the `possible_values` unify with the variables in this term */
									bool unifies = true;
									for (const auto& unification : unifications) {
										/* TODO: i dont think this is exactly correct when Contradiction is true */
										if ((!Contradiction && !new_new_values.unify_value(unification.key->quantifier.variable - 1, unification.value))
										 || (Contradiction && !new_new_values.antiunify_value(unification.key->quantifier.variable - 1, unification.value)))
										{
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
							}
						}
					}
					if (new_possible_values.length > old_size) {
						insertion_sort(new_possible_values.data + old_size, new_possible_values.length - old_size, default_sorter());
						array<variable_assignment> union_result(new_possible_values.length);
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
				possible_values.clear(); return false;
			} else {
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return true;
			}
		} else if (formula->type == FormulaType::FALSE) {
			if (Contradiction) {
				for (auto& element : new_possible_values) core::free(element);
				for (variable_assignment& assignment : possible_values)
					if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
				return true;
			} else {
				for (auto& element : new_possible_values) core::free(element);
				for (auto& element : possible_values) core::free(element);
				possible_values.clear(); return false;
			}
		} else {
			fprintf(stderr, "theory.is_provable_without_abduction ERROR: Unsupported FormulaType.\n");
			for (variable_assignment& assignment : possible_values)
				if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
			return false;
		}
		for (variable_assignment& assignment : possible_values)
			if (assignment.matching_axiom == formula) assignment.matching_axiom = nullptr;
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

	template<bool ResolveInconsistencies, bool SubtractProvableElements, bool AddSelfToAtoms>
	bool check_set_membership(unsigned int set_id,
			array<variable_assignment>& possible_values,
			const array<unsigned int>& new_antecedents,
			array_map<tuple, array<unsigned int>>& new_elements,
			array<Formula*>& quantifiers,
			array<Formula*>& new_atoms,
			unsigned int new_atom_index,
			array<unsigned int>& visited_sets)
	{
		/* check if any of the possible values of the quantified variables are newly provable members of this set */
		set_membership_prover prover(sets, implication_axioms, new_elements, new_antecedents);
		prover.h.set_ids[0] = 0;
		prover.h.set_ids.length = 1;
		typename set_membership_prover::ProvableElementArray provable_elements = prover.get_provable_elements(set_id);
		unsigned int provable_set_size = provable_elements.length;
		for (unsigned int j = 0; SubtractProvableElements && j < possible_values.length; j++) {
			const instantiation_tuple& values = possible_values[j].assignment;
			tuple new_tuple;
			if (!::init(new_tuple, values.length)) {
				prover.free_provable_elements(provable_elements);
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
						prover.free_provable_elements(provable_elements);
						core::free(new_tuple.elements); return false;
					}
				} else if (values.values[k].type == instantiation_type::ANY || values.values[k].type == instantiation_type::ANY_NUMBER) {
					for (unsigned int l = 0; l < k; l++) core::free(new_tuple[l]);
					core::free(new_tuple.elements);
					has_any = true;
					break;
				}
			}
			if (has_any) {
				/* the number of provable elements is infinite */
				prover.free_provable_elements(provable_elements);
				return false;
			}

			bool is_old = provable_elements.contains(new_tuple);
			core::free(new_tuple);
			if (is_old) {
				core::free(possible_values[j]);
				move(possible_values[possible_values.length - 1], possible_values[j]);
				possible_values.length--;
				j--;
			}
		}
		prover.free_provable_elements(provable_elements);

		if (possible_values.length == 0)
			return true;

		/* all values in `possible_values` are newly provable */
		hash_set<unsigned int> visited(16);
		hash_set<unsigned int> extensional_ancestors(8);
		array<unsigned int> stack(8);
		stack[stack.length++] = set_id;
		visited.add(set_id);
		while (stack.length != 0) {
			unsigned int current = stack.pop();
			for (const auto& entry : sets.extensional_graph.vertices[current].parents)
				if (!extensional_ancestors.add(entry.key)) return false;

			for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent))
					return false;
			}
		}
		for (unsigned int extensional_ancestor : extensional_ancestors) {
			Formula* other_formula = sets.sets[extensional_ancestor].set_formula();

			array<variable_assignment> copy(possible_values.length);
			for (unsigned int j = 0; j < possible_values.length; j++) {
				if (!::init(copy[copy.length], possible_values[j]))
					return false;
				copy.length++;
			}

			prover.h.set_ids[0] = extensional_ancestor;
			prover.h.set_ids.length = 1;
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
			if (!::init(new_tuple, possible_values[0].assignment.length))
				return false;

			const instantiation_tuple& values = possible_values[j].assignment;
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
			array<pair<Formula*, bool>> new_ancestor_terms(8);
			stack[stack.length++] = set_id;
			visited.clear();
			visited.add(set_id);
			while (stack.length != 0) {
				unsigned int current = stack.pop();
				unsigned int provable_element_count = prover.get_provable_element_count(current);
				if (AddSelfToAtoms || current != set_id) {
					Formula* set_formula = sets.sets[current].set_formula();
					if (!new_ancestor_terms.add(make_pair(set_formula, provable_element_count == sets.sets[current].set_size))) {
						return false;
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

			for (pair<Formula*, bool> new_descendant_term : new_ancestor_terms) {
				Formula* substituted;
				if (new_descendant_term.value) {
					substituted = new_descendant_term.key;
					substituted->reference_count++;
				} else {
					substituted = substitute_all(new_descendant_term.key, src_variables, dst_constants, tup.length);
					if (substituted == nullptr) {
						for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
						for (unsigned int j = 0; j < tup.length; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
						return false;
					}
				}

				add_all_atoms<true>(new_atoms, substituted, new_atom_index);
				core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
			}

			for (unsigned int j = 0; j < tup.length; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
			for (unsigned int j = 0; j < tup.length; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
		}
		visited_sets.clear();
		return true;
	}

	bool check_antecedent_satisfaction(
			unsigned int implication_id,
			array<unsigned int>& new_antecedents,
			array_map<tuple, array<unsigned int>>& new_elements,
			array<Formula*>& new_atoms,
			unsigned int new_atom_index,
			array<unsigned int>& visited_sets)
	{
		/* check if the antecedent is already satisfied */
		if (implication_axioms[implication_id].value || new_antecedents.contains(implication_id))
			return true;

		Proof* axiom = implication_axioms[implication_id].key;
		Formula* consequent = axiom->formula->binary.right;
		array<Formula*> quantifiers(4);
		set_membership_prover prover(sets, implication_axioms, new_elements, new_antecedents);
		prover.h.set_ids[0] = 0;
		prover.h.set_ids.length = 1;
		array<variable_assignment> possible_values(1);
		if (!::init(possible_values[0], 0))
			return false;
		possible_values.length = 1;
		if (is_provable_without_abduction<true>(consequent, quantifiers, possible_values, prover)) {
			for (auto& element : possible_values) core::free(element);
			return false;
		}
		for (auto& element : possible_values) core::free(element);

		if (!new_antecedents.add(implication_id))
			return false;

		if (!add_all_atoms<true>(new_atoms, consequent, new_atom_index))
			return false;

		visited_sets.clear();
		return true;
	}

	template<bool Polarity>
	bool add_all_atoms(array<Formula*>& atoms, Formula* formula, unsigned int start_index) {
		if (formula->type == FormulaType::UNARY_APPLICATION || formula->type == FormulaType::BINARY_APPLICATION || formula->type == FormulaType::EQUALS)
		{
			Formula* new_atom;
			if (Polarity) {
				new_atom = formula;
			} else {
				new_atom = Formula::new_not(formula);
				if (new_atom == nullptr)
					return false;
			}
			formula->reference_count++;

			bool atom_already_exists = false;
			for (unsigned int k = start_index; !atom_already_exists && k < atoms.length; k++)
				if (*atoms[k] == *new_atom) atom_already_exists = true;
			if (!atom_already_exists) {
				if (!atoms.add(new_atom)) {
					core::free(*new_atom); if (new_atom->reference_count == 0) core::free(new_atom);
					return false;
				}
			} else {
				core::free(*new_atom); if (new_atom->reference_count == 0) core::free(new_atom);
			}
			return true;
		}

		switch (formula->type) {
		case FormulaType::TRUE:
		case FormulaType::FALSE:
		case FormulaType::CONSTANT:
		case FormulaType::VARIABLE:
		case FormulaType::VARIABLE_PREIMAGE:
		case FormulaType::PARAMETER:
		case FormulaType::NUMBER:
		case FormulaType::STRING:
		case FormulaType::UINT_LIST:
		case FormulaType::UNARY_APPLICATION:
		case FormulaType::BINARY_APPLICATION:
		case FormulaType::EQUALS:
			return true;
		case FormulaType::IF_THEN:
			return add_all_atoms<!Polarity>(atoms, formula->binary.left, start_index)
				&& add_all_atoms<Polarity>(atoms, formula->binary.right, start_index);
		case FormulaType::NOT:
			return add_all_atoms<!Polarity>(atoms, formula->unary.operand, start_index);
		case FormulaType::AND:
		case FormulaType::OR:
		case FormulaType::IFF:
			for (unsigned int i = 0; i < formula->array.length; i++)
				if (!add_all_atoms<Polarity>(atoms, formula->array.operands[i], start_index)) return false;
			return true;
		case FormulaType::FOR_ALL:
		case FormulaType::EXISTS:
		case FormulaType::LAMBDA:
			return add_all_atoms<Polarity>(atoms, formula->quantifier.operand, start_index);
		case FormulaType::ANY:
		case FormulaType::ANY_RIGHT:
		case FormulaType::ANY_RIGHT_ONLY:
		case FormulaType::ANY_ARRAY:
		case FormulaType::ANY_CONSTANT:
		case FormulaType::ANY_CONSTANT_EXCEPT:
		case FormulaType::ANY_QUANTIFIER:
			break;
		}
		fprintf(stderr, "theory.add_all_atoms ERROR: Unrecognized FormulaType.\n");
		return false;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_set_membership_after_addition(array<Formula*>& new_atoms, array<unsigned int>& new_antecedents, Args&&... visitor) {
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		/* We first need to find all pairs of tuples and sets such that, with
		   the addition of `new_atom` as an axiom, the tuple must necessarily
		   belong to the set (and is not the case otherwise without
		   `new_atom`). To do this, find all sets that contain an atom that
		   unifies with `new_atom`. */
		unsigned int new_atom_index = 0;
		array_map<tuple, array<unsigned int>> new_elements(4);
		array<unsigned int> visited_sets(16);
		while (new_atom_index < new_atoms.length)
		{
			/* We have to consider all atoms in `new_atoms` to determine which
			   sets to visit. Otherwise, if for example we have the following
			   next atoms `person(x1)` and `x1(bob)`, we would visit the set
			   `^[x1,x2]:x1(x2)` only once to check if any elements of the form
			   (person,x1) were provable elements, then the set would be added
			   to `visited` and we would be unable to check whether any
			   elements of the form (x1,bob) were provable. */
			array_map<unsigned int, array<variable_assignment>> possible_values(4);
			array<array<unsigned int>> unifying_conjuncts(4);
			array<unsigned int> unifying_antecedents(4);
			array<unsigned int> unifying_consequents(4);
			while (new_atom_index < new_atoms.length) {
				Formula* next_new_atom = new_atoms[new_atom_index++];

				for (unsigned int i = 1; i < sets.set_count + 1; i++) {
					if (sets.sets[i].size_axioms.data == nullptr || visited_sets.contains(i))
						continue;
					Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

					array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
					for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
						quantifiers[quantifiers.length++] = set_formula;
						set_formula = set_formula->quantifier.operand;
					}
					array<array_map<Formula*, Term*>> unifications(4);
					array<unsigned int> curr_unifying_conjuncts(4);
					if (set_formula->type == FormulaType::AND) {
						for (unsigned int i = 0; i < set_formula->array.length; i++) {
							unsigned int old_unification_count = unifications.length;
							if (!unify_subformula(next_new_atom, set_formula->array.operands[i], quantifiers, unifications)) {
								for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
								for (auto entry : possible_values) {
									for (auto& value : entry.value) core::free(value);
									core::free(entry.value);
								}
								for (auto& a : unifying_conjuncts) core::free(a);
								return false;
							}
							if (unifications.length != old_unification_count)
								curr_unifying_conjuncts.add(i);
						}
					} else {
						if (!unify_subformula(next_new_atom, set_formula, quantifiers, unifications)) {
							for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
							for (auto entry : possible_values) {
								for (auto& value : entry.value) core::free(value);
								core::free(entry.value);
							}
							for (auto& a : unifying_conjuncts) core::free(a);
							return false;
						}
						if (unifications.length != 0)
							curr_unifying_conjuncts[curr_unifying_conjuncts.length++] = 0;
					}
					if (unifications.length == 0)
						continue;

					if (!possible_values.ensure_capacity(possible_values.size + 1)
					 || !unifying_conjuncts.ensure_capacity(unifying_conjuncts.length + 1))
					{
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						return false;
					}
					unsigned int index = possible_values.index_of(i);
					if (index == possible_values.size) {
						possible_values.keys[index] = i;
						if (!array_init(possible_values.values[index], max((size_t) 1, unifications.length))) {
							for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
							for (auto entry : possible_values) {
								for (auto& value : entry.value) core::free(value);
								core::free(entry.value);
							}
							for (auto& a : unifying_conjuncts) core::free(a);
							return false;
						}
						array_init(unifying_conjuncts[index], curr_unifying_conjuncts.length);
						possible_values.size++;
						unifying_conjuncts.length++;
					} else if (!possible_values.values[index].ensure_capacity(possible_values.values[index].length + unifications.length)
							|| !unifying_conjuncts[index].ensure_capacity(unifying_conjuncts[index].length + curr_unifying_conjuncts.length))
					{
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						return false;
					}
					unifying_conjuncts[index].append(curr_unifying_conjuncts.data, curr_unifying_conjuncts.length);

					array<variable_assignment>& curr_possible_values = possible_values.values[index];
					for (const array_map<Formula*, Term*>& unification : unifications) {
						/* given `unification`, find the values of the quantified variables that make `set_formula` necessarily true */
						bool valid_unification = true;
						variable_assignment& values = curr_possible_values[curr_possible_values.length];
						if (!::init(values, sets.sets[i].arity)) {
							for (auto& element : unifications) core::free(element);
							for (auto entry : possible_values) {
								for (auto& value : entry.value) core::free(value);
								core::free(entry.value);
							}
							for (auto& a : unifying_conjuncts) core::free(a);
							for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
							return false;
						}
						for (const auto& entry : unification) {
							if (entry.value->type != TermType::CONSTANT && entry.value->type != TermType::NUMBER && entry.value->type != TermType::STRING && entry.value->type != TermType::VARIABLE) {
								valid_unification = false;
								break;
							} else if (entry.key->quantifier.variable <= sets.sets[i].arity) {
								if (entry.value->type == TermType::CONSTANT) {
									core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
									values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
									values.assignment.values[entry.key->quantifier.variable - 1].constant = entry.value->constant;
								} else if (entry.value->type == TermType::NUMBER) {
									core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
									values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::NUMBER;
									values.assignment.values[entry.key->quantifier.variable - 1].number = entry.value->number;
								} else if (entry.value->type == TermType::STRING) {
									core::free(values.assignment.values[entry.key->quantifier.variable - 1]);
									values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::STRING;
									if (!core::init(values.assignment.values[entry.key->quantifier.variable - 1].str, entry.value->str)) {
										values.assignment.values[entry.key->quantifier.variable - 1].type = instantiation_type::CONSTANT;
										for (auto& element : unifications) core::free(element);
										for (auto entry : possible_values) {
											for (auto& value : entry.value) core::free(value);
											core::free(entry.value);
										}
										for (auto& a : unifying_conjuncts) core::free(a);
										for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
										core::free(values); return false;
									}
								}
							}
						}
						if (!valid_unification) {
							core::free(values);
							continue;
						}
						curr_possible_values.length++;
					}
					for (auto& element : unifications) core::free(element);
				}

				for (unsigned int i = 0; i < implication_axioms.length; i++) {
					Proof* axiom = implication_axioms[i].key;
					Formula* antecedent = axiom->formula->binary.left;

					array<Formula*> quantifiers(4);
					array<array_map<Formula*, Term*>> unifications(4);
					if (!unify_subformula<true>(next_new_atom, antecedent, quantifiers, unifications)) {
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						return false;
					}

					bool has_valid_unification = false;
					for (const array_map<Formula*, Term*>& unification : unifications) {
						bool valid_unification = true;
						for (const auto& entry : unification) {
							if (entry.value->type != TermType::CONSTANT && entry.value->type != TermType::NUMBER && entry.value->type != TermType::STRING && entry.value->type != TermType::VARIABLE) {
								valid_unification = false;
								break;
							}
						}
						if (valid_unification) {
							has_valid_unification = true;
							break;
						}
					}
					for (auto& element : unifications) core::free(element);
					if (has_valid_unification)
						unifying_antecedents.add(i);
				}
				for (unsigned int i = 0; i < implication_axioms.length; i++) {
					if (!implication_axioms[i].value) continue;
					Proof* axiom = implication_axioms[i].key;
					Formula* consequent = axiom->formula->binary.right;

					array<Formula*> quantifiers(4);
					array<array_map<Formula*, Term*>> unifications(4);
					if (!unify_subformula<true>(next_new_atom, consequent, quantifiers, unifications)) {
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						return false;
					}

					bool has_valid_unification = false;
					for (const array_map<Formula*, Term*>& unification : unifications) {
						bool valid_unification = true;
						for (const auto& entry : unification) {
							if (entry.value->type != TermType::CONSTANT && entry.value->type != TermType::NUMBER && entry.value->type != TermType::STRING && entry.value->type != TermType::VARIABLE) {
								valid_unification = false;
								break;
							}
						}
						if (valid_unification) {
							has_valid_unification = true;
							break;
						}
					}
					for (auto& element : unifications) core::free(element);
					if (!has_valid_unification)
						unifying_consequents.add(i);
				}
			}

			for (unsigned int index = 0; index < possible_values.size; index++) {
				unsigned int i = possible_values.keys[index];
				array<variable_assignment>& curr_possible_values = possible_values.values[index];
				array<unsigned int>& curr_unifying_conjuncts = unifying_conjuncts[index];
				Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;
				array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}

				if (!visited_sets.add(i)) {
					for (auto entry : possible_values) {
						for (auto& value : entry.value) core::free(value);
						core::free(entry.value);
					}
					for (auto& a : unifying_conjuncts) core::free(a);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}

				/* remove possible_values that are subsets of another possible_value */
				for (unsigned int i = 0; i < curr_possible_values.length; i++) {
					const instantiation_tuple& first = curr_possible_values[i].assignment;
					for (unsigned int j = 0; j < curr_possible_values.length; j++) {
						if (j == i) continue;
						const instantiation_tuple& second = curr_possible_values[j].assignment;

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
							core::free(curr_possible_values[i]);
							curr_possible_values.remove(i--);
							break;
						}
					}
				}

				set_membership_prover prover(sets, implication_axioms, new_elements, new_antecedents);
				prover.h.set_ids[0] = i;
				prover.h.set_ids.length = 1;
				typename set_membership_prover::ProvableElementArray provable_elements = prover.get_provable_elements(i);
				if (provable_elements.length == sets.sets[i].set_size && set_formula->type == FormulaType::AND) {
					/* this set being full could make other expressions provable by exclusion */
					for (const variable_assignment& possible_value : curr_possible_values) {
						/* optimization: if the set formula has a name scope, only consider proof-by-exclusion if the `possible_value` has that name */
						bool skip = false;
						for (unsigned int k = 0; k < set_formula->array.length; k++) {
							const Term* arg1; const Term* arg2;
							if (is_name_scope(set_formula->array.operands[k], arg1, arg2)
							 && arg1->type == TermType::VARIABLE && arg2->type == TermType::STRING)
							{
								if (possible_value.assignment.values[arg1->variable - 1].type == instantiation_type::CONSTANT
								 && !is_concept_name(possible_value.assignment.values[arg1->variable - 1].constant, arg2->str))
								{
									skip = true;
									break;
								} else if (possible_value.assignment.values[arg1->variable - 1].type == instantiation_type::NUMBER
										|| possible_value.assignment.values[arg1->variable - 1].type == instantiation_type::STRING)
								{
									skip = true;
									break;
								}
							}
						}
						if (skip) continue;

						unsigned int substitution_count = 0;
						Term** src_variables = (Term**) malloc(sizeof(Term*) * sets.sets[i].arity * 2);
						if (src_variables == nullptr) {
							for (auto entry : possible_values) {
								for (auto& value : entry.value) core::free(value);
								core::free(entry.value);
							}
							for (auto& a : unifying_conjuncts) core::free(a);
							for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
							return false;
						}
						Term** dst_variables = src_variables + sets.sets[i].arity;
						for (uint_fast8_t k = 0; k < possible_value.assignment.length; k++) {
							if (possible_value.assignment.values[k].type == instantiation_type::CONSTANT) {
								src_variables[substitution_count] = Formula::new_variable(k + 1);
								dst_variables[substitution_count] = Formula::new_constant(possible_value.assignment.values[k].constant);
							} else if (possible_value.assignment.values[k].type == instantiation_type::NUMBER) {
								src_variables[substitution_count] = Formula::new_variable(k + 1);
								dst_variables[substitution_count] = Formula::new_number(possible_value.assignment.values[k].number);
							} else if (possible_value.assignment.values[k].type == instantiation_type::STRING) {
								src_variables[substitution_count] = Formula::new_variable(k + 1);
								dst_variables[substitution_count] = Formula::new_string(possible_value.assignment.values[k].str);
							} else {
								continue;
							}
							if (src_variables[substitution_count] == nullptr || dst_variables[substitution_count] == nullptr) {
								if (src_variables[substitution_count] != nullptr) { core::free(*src_variables[substitution_count]); core::free(src_variables[substitution_count]); }
								if (dst_variables[substitution_count] != nullptr) { core::free(*dst_variables[substitution_count]); core::free(dst_variables[substitution_count]); }
								for (unsigned int j = 0; j < substitution_count; j++) {
									core::free(*src_variables[j]); core::free(src_variables[j]);
									core::free(*dst_variables[j]); core::free(dst_variables[j]);
								}
								core::free(src_variables);
								for (auto entry : possible_values) {
									for (auto& value : entry.value) core::free(value);
									core::free(entry.value);
								}
								for (auto& a : unifying_conjuncts) core::free(a);
								for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
								return false;
							}
							substitution_count++;
						}

						for (unsigned int k = 0; k < set_formula->array.length; k++) {
							/* we only check conjuncts that don't unify with `next_new_atom` since
							   we're already checking those sets anyway in this loop (this is why we
							   require `set_formula` to be a conjunction) */
							if (curr_unifying_conjuncts.contains(k)) continue;
							Formula* conjunct = set_formula->array.operands[k];
							Formula* substituted;
							if (substitution_count == 0) {
								substituted = conjunct;
								conjunct->reference_count++;
							} else {
								substituted = substitute_all(conjunct, src_variables, dst_variables, substitution_count);
								if (substituted == nullptr) {
									for (unsigned int j = 0; j < substitution_count; j++) {
										core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]);
										core::free(*dst_variables[j]); if (dst_variables[j]->reference_count == 0) core::free(dst_variables[j]);
									}
									core::free(src_variables);
									for (auto entry : possible_values) {
										for (auto& value : entry.value) core::free(value);
										core::free(entry.value);
									}
									for (auto& a : unifying_conjuncts) core::free(a);
									for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
									return false;
								}
							}

							unsigned int old_atom_count = new_atoms.length;
							add_all_atoms<true>(new_atoms, substituted, 0);
							if (new_atoms.length != old_atom_count)
								visited_sets.clear();
							core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
						}
						for (unsigned int j = 0; j < substitution_count; j++) {
							core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]);
							core::free(*dst_variables[j]); if (dst_variables[j]->reference_count == 0) core::free(dst_variables[j]);
						}
						core::free(src_variables);
					}
				}

				array<variable_assignment> new_possible_values(max((size_t) 1, 4 * curr_possible_values.length));
				for (const variable_assignment& possible_value : curr_possible_values) {
					array<variable_assignment> temp_possible_values(4);
					variable_assignment& values = temp_possible_values[0];
					if (!::init(values, possible_value)) {
						prover.free_provable_elements(provable_elements);
						for (auto& element : new_possible_values) core::free(element);
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						return false;
					}
					temp_possible_values.length++;

					if (!is_provable_without_abduction<false>(set_formula, quantifiers, temp_possible_values, prover)) {
						for (auto& element : temp_possible_values) core::free(element);
						continue;
					}

					array<variable_assignment> union_result(new_possible_values.length + temp_possible_values.length);
					set_union(union_result, new_possible_values, temp_possible_values);
					for (auto& element : new_possible_values) core::free(element);
					for (auto& element : temp_possible_values) core::free(element);
					swap(union_result, new_possible_values);
				}

				/* the addition could cause this formula to be provably false */
				array<variable_assignment> temp_possible_values(max(1, provable_elements.length));
				for (const tuple& element : provable_elements) {
					variable_assignment& values = temp_possible_values[temp_possible_values.length];
					if (!::init(values, element.length)) {
						for (auto& element : temp_possible_values) core::free(element);
						prover.free_provable_elements(provable_elements);
						for (auto& element : new_possible_values) core::free(element);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						for (auto entry : possible_values) {
							for (auto& value : entry.value) core::free(value);
							core::free(entry.value);
						}
						for (auto& a : unifying_conjuncts) core::free(a);
						return false;
					}
					for (uint_fast8_t k = 0; k < element.length; k++) {
						if (element[k].type == tuple_element_type::CONSTANT) {
							core::free(values.assignment.values[k]);
							values.assignment.values[k].type = instantiation_type::CONSTANT;
							values.assignment.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							core::free(values.assignment.values[k]);
							values.assignment.values[k].type = instantiation_type::NUMBER;
							values.assignment.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							core::free(values.assignment.values[k]);
							values.assignment.values[k].type = instantiation_type::STRING;
							if (!core::init(values.assignment.values[k].str, element[k].str)) {
								values.assignment.values[k].type = instantiation_type::CONSTANT;
								for (auto& element : temp_possible_values) core::free(element);
								prover.free_provable_elements(provable_elements);
								for (auto entry : possible_values) {
									for (auto& value : entry.value) core::free(value);
									core::free(entry.value);
								}
								for (auto& a : unifying_conjuncts) core::free(a);
								for (auto& element : new_possible_values) core::free(element);
								for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
								values.assignment.values[k].type = instantiation_type::CONSTANT;
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;
				}
				prover.free_provable_elements(provable_elements);
				if (temp_possible_values.length != 0 && is_provable_without_abduction<true>(set_formula, quantifiers, temp_possible_values, prover)) {
					/* we found a contradiction */
					for (auto& element : temp_possible_values) core::free(element);
					for (auto& element : new_possible_values) core::free(element);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					for (auto entry : possible_values) {
						for (auto& value : entry.value) core::free(value);
						core::free(entry.value);
					}
					for (auto& a : unifying_conjuncts) core::free(a);
					return false;
				}

				if (new_possible_values.length == 0) continue;
				if (!check_set_membership<ResolveInconsistencies, true, true>(i, new_possible_values, new_antecedents, new_elements, quantifiers, new_atoms, new_atom_index + 1, visited_sets)) {
					for (auto& element : new_possible_values) core::free(element);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					for (auto entry : possible_values) {
						for (auto& value : entry.value) core::free(value);
						core::free(entry.value);
					}
					for (auto& a : unifying_conjuncts) core::free(a);
					return false;
				}
				for (auto& element : new_possible_values) core::free(element);
			}
			for (auto entry : possible_values) {
				for (auto& value : entry.value) core::free(value);
				core::free(entry.value);
			}
			for (auto& a : unifying_conjuncts) core::free(a);

			for (unsigned int i : unifying_antecedents) {
				Proof* axiom = implication_axioms[i].key;
				Formula* antecedent = axiom->formula->binary.left;

				set_membership_prover prover(sets, implication_axioms, new_elements, new_antecedents);
				prover.h.set_ids[0] = 0;
				prover.h.set_ids.length = 1;
				array<variable_assignment> temp_possible_values(1);
				if (!::init(temp_possible_values[0], 0)) {
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}
				temp_possible_values.length = 1;
				array<Formula*> quantifiers(4);
				if (!is_provable_without_abduction<false>(antecedent, quantifiers, temp_possible_values, prover)) {
					for (auto& element : temp_possible_values) core::free(element);
					continue;
				}
				for (auto& element : temp_possible_values) core::free(element);

				if (!check_antecedent_satisfaction(i, new_antecedents, new_elements, new_atoms, new_atom_index + 1, visited_sets)) {
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}
			}
			for (unsigned int i : unifying_consequents) {
				if (!implication_axioms[i].value) continue;
				Proof* axiom = implication_axioms[i].key;
				Formula* consequent = axiom->formula->binary.right;

				/* the addition could cause this formula to be provably false */
				set_membership_prover prover(sets, implication_axioms, new_elements, new_antecedents);
				prover.h.set_ids[0] = 0;
				prover.h.set_ids.length = 1;
				array<variable_assignment> temp_possible_values(1);
				if (!::init(temp_possible_values[0], 0)) {
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					return false;
				}
				temp_possible_values.length = 1;
				array<Formula*> quantifiers(4);
				if (is_provable_without_abduction<true>(consequent, quantifiers, temp_possible_values, prover)) {
					for (auto& element : temp_possible_values) core::free(element);
					return false;
				}
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
				sets.remove_element_at(old_set, sets.sets[old_set].index_of_element(element));
			for (unsigned int new_set : new_elements.values[i]) {
				if (!sets.add_element(new_set, element)) {
					for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
					for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
					core::free(old_elements); return false;
				}
			}
		}

		/* set the appropriate antecedents to be provable */
		for (unsigned int index : new_antecedents)
			implication_axioms[index].value = true;

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
				if (sets.sets[new_ancestor].set_size >= sets.sets[new_ancestor].provable_elements.length)
					continue;

				while (true) {
					/* compute the upper bound on the size of this set; if the new size
						violates this bound, change the sizes of other sets to increase the bound */
					array<unsigned int> stack(8); bool graph_changed;
					if (!ResolveInconsistencies || !sets.increase_set_size(new_ancestor, sets.sets[new_ancestor].provable_elements.length, stack, graph_changed)) {
						/* undo the set changes */
						for (unsigned int i = 0; i < new_elements.size; i++) {
							const tuple& element = new_elements.keys[i];
							for (unsigned int new_set : new_elements.values[i])
								sets.remove_element_at(new_set, sets.sets[new_set].index_of_element(element));
							for (unsigned int old_set : old_elements[i])
								sets.add_element(old_set, element);
						}
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
					if (sets.sets[new_ancestor].set_size >= sets.sets[new_ancestor].provable_elements.length)
						break;
					if (!graph_changed) {
						/* undo the set changes */
						for (unsigned int i = 0; i < new_elements.size; i++) {
							const tuple& element = new_elements.keys[i];
							for (unsigned int new_set : new_elements.values[i])
								sets.remove_element_at(new_set, sets.sets[new_set].index_of_element(element));
							for (unsigned int old_set : old_elements[i])
								sets.add_element(old_set, element);
						}
						for (unsigned int j = 0; j < new_elements.size; j++) core::free(old_elements[j]);
						for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
						core::free(old_elements); return false;
					}
				}
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
		if (!add_all_atoms<true>(new_atoms, new_atom, 0))
			return false;
		array<unsigned int> new_antecedents(4);
		bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, new_antecedents, std::forward<Args>(visitor)...);
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

		sets.remove_element_at(element_set_id, element_index);
		for (unsigned int l = 0; l < conjuncts.length; l++) {
			array<variable_assignment> temp_possible_values(1);
			variable_assignment& values = temp_possible_values[0];
			if (!::init(values, sets.sets[element_set_id].arity)) {
				sets.add_to_provable_elements(element_set_id, element);
				for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++)
					move(element[k], sets.sets[element_set_id].elements[sets.sets[element_set_id].elements.length++]);
				return false;
			}
			for (unsigned int k = 0; k < values.assignment.length; k++) {
				core::free(values.assignment.values[k]);
				if (element[k].type == tuple_element_type::CONSTANT) {
					values.assignment.values[k].type = instantiation_type::CONSTANT;
					values.assignment.values[k].constant = element[k].constant;
				} else if (element[k].type == tuple_element_type::NUMBER) {
					values.assignment.values[k].type = instantiation_type::NUMBER;
					values.assignment.values[k].number = element[k].number;
				} else if (element[k].type == tuple_element_type::STRING) {
					values.assignment.values[k].type = instantiation_type::STRING;
					if (!core::init(values.assignment.values[k].str, element[k].str)) {
						values.assignment.values[k].type = instantiation_type::CONSTANT;
						sets.add_to_provable_elements(element_set_id, element);
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++)
							move(element[k], sets.sets[element_set_id].elements[sets.sets[element_set_id].elements.length++]);
						core::free(values); return false;
					}
				}
			}
			temp_possible_values.length++;

			default_prover prover(sets, implication_axioms);
			prover.h.set_ids[0] = formula_set_id;
			prover.h.set_ids.length = 1;
			if (is_provable_without_abduction<false>(conjuncts[l], quantifiers, temp_possible_values, prover)) {
				for (auto& element : temp_possible_values) core::free(element);
				conjuncts.remove(l--);
			}
		}

		if (conjuncts.length == 0) {
			sets.add_to_provable_elements(element_set_id, element);
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
					array<variable_assignment> temp_possible_values(1);
					variable_assignment& values = temp_possible_values[0];
					if (!::init(values, sets.sets[element_set_id].arity)) {
						for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
						return false;
					}
					for (unsigned int k = 0; k < values.assignment.length; k++) {
						core::free(values.assignment.values[k]);
						if (element[k].type == tuple_element_type::CONSTANT) {
							values.assignment.values[k].type = instantiation_type::CONSTANT;
							values.assignment.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							values.assignment.values[k].type = instantiation_type::NUMBER;
							values.assignment.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							values.assignment.values[k].type = instantiation_type::STRING;
							if (!core::init(values.assignment.values[k].str, element[k].str)) {
								values.assignment.values[k].type = instantiation_type::CONSTANT;
								for (uint_fast8_t k = 0; k < sets.sets[element_set_id].arity; k++) core::free(element[k]);
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;

					default_prover prover(sets, implication_axioms);
					prover.h.set_ids[0] = current;
					prover.h.set_ids.length = 1;
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
				if (!sets.add_element(new_set, element)) {
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
	bool check_set_membership_after_subtraction(array<Formula*>& old_atoms, unsigned int removed_set, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		unsigned int old_atom_index = 0;
		array_map<unsigned int, tuple> old_elements(8);
		array<unsigned int> old_antecedents(4);
		while (old_atom_index < old_atoms.length) {
			Formula* next_old_atom = old_atoms[old_atom_index++];

			for (unsigned int i = 1; i < sets.set_count + 1; i++) {
				if (sets.sets[i].size_axioms.data == nullptr || i == removed_set)
					continue;
				Formula* set_formula = sets.sets[i].size_axioms[0]->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[i].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[i].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}
				array<array_map<Formula*, Term*>> unifications(4);
				if (!unify_subformula(next_old_atom, set_formula, quantifiers, unifications)) {
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}

				for (unsigned int j = sets.sets[i].element_count(); j > 0; j--) {
					tuple& element = *((tuple*) alloca(sizeof(tuple)));
					element.elements = (tuple_element*) malloc(sizeof(tuple_element) * sets.sets[i].arity);
					if (element.elements == nullptr) {
						fprintf(stderr, "theory.check_set_membership_after_subtraction ERROR: Out of memory.\n");
						for (auto& element : unifications) core::free(element);
						for (auto entry : old_elements) core::free(entry.value);
						return false;
					}
					element.length = sets.sets[i].arity;
					const tuple_element* element_src = sets.sets[i].elements.data + (sets.sets[i].arity * (j - 1));
					for (uint_fast8_t k = 0; k < sets.sets[i].arity; k++) {
						if (!::init(element.elements[k], element_src[k])) {
							for (auto& element : unifications) core::free(element);
							for (auto entry : old_elements) core::free(entry.value);
							for (uint_fast8_t l = 0; l < k; l++) core::free(element[l]);
							core::free(element.elements); return false;
						}
					}
					bool has_satisfying_unification = false;
					for (const array_map<Formula*, Term*>& unification : unifications) {
						bool matches = true;
						for (const auto& entry : unification) {
							if (entry.key->quantifier.variable > sets.sets[i].arity || entry.value->type == TermType::VARIABLE)
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

					sets.remove_element_at(i, j - 1);
					if (!old_elements.ensure_capacity(old_elements.size + 1)) {
						for (auto& element : unifications) core::free(element);
						for (auto entry : old_elements) core::free(entry.value);
						core::free(element); return false;
					}
					old_elements.keys[old_elements.size] = i;
					move(element, old_elements.values[old_elements.size]);
					old_elements.size++;

					/* this set is no longer full, so we may no longer be able to prove other things that relied on proof by excluding elements from this set */
					bool no_longer_full = (sets.sets[i].provable_elements.length + 1 == sets.sets[i].set_size);

					array_map<Formula*, bool> conjuncts(4);
					if (set_formula->type == FormulaType::AND) {
						for (unsigned int j = 0; j < set_formula->array.length; j++) {
							if (!conjuncts.put(set_formula->array.operands[j], no_longer_full)) {
								for (auto& element : unifications) core::free(element);
								for (auto entry : old_elements) core::free(entry.value);
								core::free(element); return false;
							}
						}
					} else {
						conjuncts.keys[0] = set_formula;
						conjuncts.values[0] = no_longer_full;
						conjuncts.size++;
					}

					/* terms in the set formulae of the ancestors of set `formula_set_id` could also be no longer provable */
					array<unsigned int> stack(8);
					hash_set<unsigned int> visited(16);
					stack[stack.length++] = i;
					visited.add(i);
					while (stack.length != 0) {
						unsigned int current = stack.pop();

						if (current != i) {
							Formula* set_formula = sets.sets[current].set_formula();
							no_longer_full = (sets.sets[current].provable_elements.length + 1 == sets.sets[current].set_size);
							if (set_formula->type == FormulaType::AND) {
								for (unsigned int j = 0; j < set_formula->array.length; j++) {
									if (!conjuncts.put(set_formula->array.operands[j], no_longer_full)) {
										for (auto& element : unifications) core::free(element);
										for (auto entry : old_elements) core::free(entry.value);
										return false;
									}
								}
							} else if (!conjuncts.put(set_formula, no_longer_full)) {
								for (auto& element : unifications) core::free(element);
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
						}

						for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
							if (visited.contains(parent)) continue;
							if (!visited.add(parent) || !stack.add(parent)) {
								for (auto& element : unifications) core::free(element);
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
						} for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
							if (visited.contains(entry.key)) continue;
							if (!visited.add(entry.key) || !stack.add(entry.key)) {
								for (auto& element : unifications) core::free(element);
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
							for (auto& element : unifications) core::free(element);
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
							for (auto& element : unifications) core::free(element);
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
					}
					for (auto conjunct : conjuncts) {
						Formula* substituted_conjunct;
						if (conjunct.value) {
							substituted_conjunct = conjunct.key;
							substituted_conjunct->reference_count++;
						} else {
							substituted_conjunct = substitute_all(conjunct.key, src_variables, dst_constants, sets.sets[i].arity);
						}
						if (substituted_conjunct == nullptr
						 || !add_all_atoms<true>(old_atoms, substituted_conjunct, old_atom_index + 1))
						{
							if (substituted_conjunct != nullptr) { core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct); }
							for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
							for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
							for (auto& element : unifications) core::free(element);
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						core::free(*substituted_conjunct); if (substituted_conjunct->reference_count == 0) core::free(substituted_conjunct);
					}
					for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*src_variables[j]); core::free(src_variables[j]); }
					for (unsigned int j = 0; j < sets.sets[i].arity; j++) { core::free(*dst_constants[j]); if (dst_constants[j]->reference_count == 0) core::free(dst_constants[j]); }
				}
				if (unifications.length == 0) continue;

				if (sets.sets[i].provable_elements.length == sets.sets[i].set_size && set_formula->type == FormulaType::AND) {
					/* this set being full could make other expressions provable by exclusion */
					/* TODO: we should only check conjuncts that don't unify
					   with `next_new_atom` since we're already checking those
					   sets anyway in this loop (this is why we require
					   `set_formula` to be a conjunction) */
					for (const array_map<Formula*, Term*>& unification : unifications) {
						unsigned int substitution_count = 0;
						Term** src_variables = (Term**) malloc(sizeof(Term*) * sets.sets[i].arity * 2);
						if (src_variables == nullptr) {
							for (auto& element : unifications) core::free(element);
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						Term** dst_variables = src_variables + sets.sets[i].arity;
						bool valid_unification = true;
						for (const auto& entry : unification) {
							if (entry.key->quantifier.variable > sets.sets[i].arity) {
								continue;
							} else if (entry.value->type == TermType::CONSTANT) {
								src_variables[substitution_count] = Formula::new_variable(entry.key->quantifier.variable);
								dst_variables[substitution_count] = Formula::new_constant(entry.value->constant);
							} else if (entry.value->type == TermType::NUMBER) {
								src_variables[substitution_count] = Formula::new_variable(entry.key->quantifier.variable);
								dst_variables[substitution_count] = Formula::new_number(entry.value->number);
							} else if (entry.value->type == TermType::STRING) {
								src_variables[substitution_count] = Formula::new_variable(entry.key->quantifier.variable);
								dst_variables[substitution_count] = Formula::new_string(entry.value->str);
							} else if (entry.value->type != TermType::VARIABLE) {
								valid_unification = false;
								break;
							} else {
								continue;
							}
							if (src_variables[substitution_count] == nullptr || dst_variables[substitution_count] == nullptr) {
								if (src_variables[substitution_count] != nullptr) { core::free(*src_variables[substitution_count]); core::free(src_variables[substitution_count]); }
								if (dst_variables[substitution_count] != nullptr) { core::free(*dst_variables[substitution_count]); core::free(dst_variables[substitution_count]); }
								for (unsigned int j = 0; j < substitution_count; j++) {
									core::free(*src_variables[j]); core::free(src_variables[j]);
									core::free(*dst_variables[j]); core::free(dst_variables[j]);
								}
								core::free(src_variables);
								for (auto& element : unifications) core::free(element);
								for (auto entry : old_elements) core::free(entry.value);
								return false;
							}
							substitution_count++;
						}
						if (!valid_unification) {
							for (unsigned int j = 0; j < substitution_count; j++) {
								core::free(*src_variables[j]); core::free(src_variables[j]);
								core::free(*dst_variables[j]); core::free(dst_variables[j]);
							}
							core::free(src_variables);
							continue;
						}

						for (unsigned int k = 0; k < set_formula->array.length; k++) {
							Formula* conjunct = set_formula->array.operands[k];
							Formula* substituted;
							if (substitution_count == 0) {
								substituted = conjunct;
								conjunct->reference_count++;
							} else {
								substituted = substitute_all(conjunct, src_variables, dst_variables, substitution_count);
								if (substituted == nullptr) {
									for (unsigned int j = 0; j < substitution_count; j++) {
										core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]);
										core::free(*dst_variables[j]); if (dst_variables[j]->reference_count == 0) core::free(dst_variables[j]);
									}
									core::free(src_variables);
									for (auto& element : unifications) core::free(element);
									for (auto entry : old_elements) core::free(entry.value);
									return false;
								}
							}

							add_all_atoms<true>(old_atoms, substituted, 0);
							core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
						}
						for (unsigned int j = 0; j < substitution_count; j++) {
							core::free(*src_variables[j]); if (src_variables[j]->reference_count == 0) core::free(src_variables[j]);
							core::free(*dst_variables[j]); if (dst_variables[j]->reference_count == 0) core::free(dst_variables[j]);
						}
						core::free(src_variables);
					}
				}
				for (auto& element : unifications) core::free(element);
			}
			for (unsigned int i = 0; i < implication_axioms.length; i++) {
				if (!implication_axioms[i].value) continue;
				Proof* axiom = implication_axioms[i].key;
				Formula* antecedent = axiom->formula->binary.left;
				Formula* consequent = axiom->formula->binary.right;

				array<Formula*> quantifiers(1);
				array<array_map<Formula*, Term*>> unifications(4);
				if (!unify_subformula<true>(next_old_atom, antecedent, quantifiers, unifications)) {
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}

				if (unifications.length == 0)
					continue;
				for (auto& element : unifications) core::free(element);

				if (!old_antecedents.add(i)) {
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}
				implication_axioms[i].value = false;

				if (!add_all_atoms<true>(old_atoms, antecedent, old_atom_index + 1)
				 || !add_all_atoms<true>(old_atoms, consequent, old_atom_index + 1))
				{
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}
			}
		}

		bool sets_changed = true;
		while (sets_changed) {
			sets_changed = false;
			for (unsigned int i = 0; i < old_elements.size; i++) {
				/* check if this removed element is still provable */
				Formula* set_formula = sets.sets[old_elements.keys[i]].size_axioms[0]->formula->binary.left->binary.right;

				array<Formula*> quantifiers(1 << (core::log2(sets.sets[old_elements.keys[i]].arity) + 1));
				for (unsigned int j = 0; j < sets.sets[old_elements.keys[i]].arity; j++) {
					quantifiers[quantifiers.length++] = set_formula;
					set_formula = set_formula->quantifier.operand;
				}

				array<variable_assignment> temp_possible_values(1);
				variable_assignment& values = temp_possible_values[0];
				if (!::init(values, sets.sets[old_elements.keys[i]].arity)) {
					sets.add_to_provable_elements(old_elements.keys[i], old_elements.values[i]);
					for (uint_fast8_t k = 0; k < sets.sets[old_elements.keys[i]].arity; k++)
						move(old_elements.values[i][k], sets.sets[old_elements.keys[i]].elements[sets.sets[old_elements.keys[i]].elements.length++]);
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}
				for (unsigned int k = 0; k < values.assignment.length; k++) {
					core::free(values.assignment.values[k]);
					if (old_elements.values[i][k].type == tuple_element_type::CONSTANT) {
						values.assignment.values[k].type = instantiation_type::CONSTANT;
						values.assignment.values[k].constant = old_elements.values[i][k].constant;
					} else if (old_elements.values[i][k].type == tuple_element_type::NUMBER) {
						values.assignment.values[k].type = instantiation_type::NUMBER;
						values.assignment.values[k].number = old_elements.values[i][k].number;
					} else if (old_elements.values[i][k].type == tuple_element_type::STRING) {
						values.assignment.values[k].type = instantiation_type::STRING;
						if (!core::init(values.assignment.values[k].str, old_elements.values[i][k].str)) {
							values.assignment.values[k].type = instantiation_type::CONSTANT;
							sets.add_to_provable_elements(old_elements.keys[i], old_elements.values[i]);
							for (uint_fast8_t k = 0; k < sets.sets[old_elements.keys[i]].arity; k++)
								move(old_elements.values[i][k], sets.sets[old_elements.keys[i]].elements[sets.sets[old_elements.keys[i]].elements.length++]);
							for (auto entry : old_elements) core::free(entry.value);
							core::free(values); return false;
						}
					}
				}
				temp_possible_values.length++;

				default_prover prover(sets, implication_axioms, removed_set);
				prover.h.set_ids[0] = old_elements.keys[i];
				prover.h.set_ids.length = 1;
				if (is_provable_without_abduction<false>(set_formula, quantifiers, temp_possible_values, prover)) {
					/* this element is still provably a member of this set */
					for (auto& element : temp_possible_values) core::free(element);
					if (!sets.add_element(old_elements.keys[i], old_elements.values[i])) {
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
						array<variable_assignment> temp_possible_values(1);
						variable_assignment& values = temp_possible_values[0];
						if (!::init(values, sets.sets[old_elements.keys[i]].arity)) {
							for (auto entry : old_elements) core::free(entry.value);
							return false;
						}
						for (unsigned int k = 0; k < values.assignment.length; k++) {
							core::free(values.assignment.values[k]);
							if (old_elements.values[i][k].type == tuple_element_type::CONSTANT) {
								values.assignment.values[k].type = instantiation_type::CONSTANT;
								values.assignment.values[k].constant = old_elements.values[i][k].constant;
							} else if (old_elements.values[i][k].type == tuple_element_type::NUMBER) {
								values.assignment.values[k].type = instantiation_type::NUMBER;
								values.assignment.values[k].number = old_elements.values[i][k].number;
							} else if (old_elements.values[i][k].type == tuple_element_type::STRING) {
								values.assignment.values[k].type = instantiation_type::STRING;
								if (!core::init(values.assignment.values[k].str, old_elements.values[i][k].str)) {
									values.assignment.values[k].type = instantiation_type::CONSTANT;
									for (auto entry : old_elements) core::free(entry.value);
									core::free(values); return false;
								}
							}
						}
						temp_possible_values.length++;

						default_prover prover(sets, implication_axioms, removed_set);
						prover.h.set_ids[0] = current;
						prover.h.set_ids.length = 1;
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
					if (!sets.add_element(new_set, old_elements.values[i])) {
						for (auto entry : old_elements) core::free(entry.value);
						return false;
					}
					/* this element could also be in `old_elements`, so remove it
					   from `old_elements` to avoid the same element being added to
					   the set multiple times */
					for (unsigned int j = 0; j < i; j++) {
						if (old_elements.keys[j] == new_set && old_elements.values[j] == old_elements.values[i]) {
							core::free(old_elements.values[j]);
							core::move(old_elements.keys[i-1],old_elements.keys[j]);
							core::move(old_elements.values[i-1],old_elements.values[j]);
							core::move(old_elements.keys[i],old_elements.keys[i-1]);
							core::move(old_elements.values[i],old_elements.values[i-1]);
							old_elements.size--;
							core::move(old_elements.keys[old_elements.size],old_elements.keys[i]);
							core::move(old_elements.values[old_elements.size],old_elements.values[i]);
							i--;
						}
					}
					for (unsigned int j = 0; j < old_elements.size; j++) {
						if (old_elements.keys[j] == new_set && old_elements.values[j] == old_elements.values[i]) {
							core::free(old_elements.values[j]);
							old_elements.remove_at(j--);
						}
					}
				}
			}
			for (unsigned int i = 0; i < old_antecedents.length; i++) {
				/* check if this antecedent is still provable */
				Proof* proof = implication_axioms[old_antecedents[i]].key;
				Formula* antecedent = proof->formula->binary.left;

				array<Formula*> quantifiers(1);
				array<variable_assignment> temp_possible_values(1);
				variable_assignment& values = temp_possible_values[0];
				if (!::init(values, 0)) {
					sets.add_to_provable_elements(old_elements.keys[i], old_elements.values[i]);
					for (uint_fast8_t k = 0; k < sets.sets[old_elements.keys[i]].arity; k++)
						move(old_elements.values[i][k], sets.sets[old_elements.keys[i]].elements[sets.sets[old_elements.keys[i]].elements.length++]);
					for (auto entry : old_elements) core::free(entry.value);
					return false;
				}
				temp_possible_values.length++;

				default_prover prover(sets, implication_axioms, removed_set);
				prover.h.set_ids[0] = 0;
				prover.h.set_ids.length = 1;
				if (is_provable_without_abduction<false>(antecedent, quantifiers, temp_possible_values, prover)) {
					/* this antecedent is still provable */
					for (auto& element : temp_possible_values) core::free(element);
					implication_axioms[old_antecedents[i]].value = true;
					old_antecedents.remove(i);
					sets_changed = true;
					continue;
				}
			}
		}
		for (auto entry : old_elements) core::free(entry.value);
		return true;
	}

	template<typename... Args>
	inline bool check_set_membership_after_subtraction(Formula* old_atom, unsigned int removed_set, Args&&... visitor)
	{
		array<Formula*> old_atoms(8);
		if (!add_all_atoms<true>(old_atoms, old_atom, 0))
			return false;
		if (old_atom->type == FormulaType::EQUALS) {
			Formula* swapped = Formula::new_equals(old_atom->binary.right, old_atom->binary.left);
			if (swapped == nullptr) {
				free_all(old_atoms);
				return false;
			}
			old_atom->binary.left->reference_count++;
			old_atom->binary.right->reference_count++;
			old_atoms[old_atoms.length++] = swapped;
		}
		bool result = check_set_membership_after_subtraction(old_atoms, removed_set, std::forward<Args>(visitor)...);
		free_all(old_atoms);
		return result;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline bool check_new_set_membership(
			unsigned int set_id, array<Formula*>& quantifiers,
			array<variable_assignment>& possible_values, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		array<unsigned int> new_antecedents(1);
		array_map<tuple, array<unsigned int>> new_elements(1);
		array<Formula*> new_atoms(1);
		array<unsigned int> visited_sets(1);
		if (!check_set_membership<ResolveInconsistencies, true, false>(set_id, possible_values, new_antecedents, new_elements, quantifiers, new_atoms, 0, visited_sets)) {
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			return false;
		}
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
		if (sets.sets[set_id].provable_elements.length + new_element_count > sets.sets[set_id].set_size) {
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			return false;
		}

		/* add the elements to the sets */
		hash_set<unsigned int> ancestors(16);
		if (!sets.template get_ancestors<true>(set_id, ancestors)) {
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			return false;
		}
		for (const auto& entry : new_elements) {
			if (!entry.value.contains(set_id)) continue;
			if (!sets.add_element(set_id, entry.key, ancestors)) {
				for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
				return false;
			}
		}
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline bool check_new_set_membership(unsigned int set_id, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		Formula* set_formula = sets.sets[set_id].size_axioms[0]->formula->binary.left->binary.right;

		array<Formula*> quantifiers(1 << (core::log2(sets.sets[set_id].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[set_id].arity; j++) {
			quantifiers[quantifiers.length++] = set_formula;
			set_formula = set_formula->quantifier.operand;
		}

		/* if this set is empty, negations of existentials could now be provable in other set formula */
		if (sets.sets[set_id].set_size == 0) {
			array<Formula*> new_atoms(4);
			if (!add_all_atoms<true>(new_atoms, set_formula, 0))
				return false;
			array<unsigned int> new_antecedents(4);
			bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, new_antecedents, std::forward<Args>(visitor)...);
			free_all(new_atoms);
			return result;
		}

		array<variable_assignment> possible_values(4);
		variable_assignment& values = possible_values[0];
		if (!::init(values, sets.sets[set_id].arity))
			return false;
		possible_values.length++;

		default_prover prover(sets, implication_axioms);
		prover.h.set_ids[0] = set_id;
		prover.h.set_ids.length = 1;
		if (!is_provable_without_abduction<false>(set_formula, quantifiers, possible_values, prover)) {
			for (auto& element : possible_values) core::free(element);
			return true;
		}

		if (!check_new_set_membership<ResolveInconsistencies>(set_id, quantifiers, possible_values, std::forward<Args>(visitor)...)) {
			for (auto& element : possible_values) core::free(element);
			return false;
		}
		for (auto& element : possible_values) core::free(element);
		return true;
	}

	template<typename... Args>
	bool check_old_set_membership(unsigned int set_id, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		Formula* set_formula = sets.sets[set_id].size_axioms[0]->formula->binary.left->binary.right;
		array<Formula*> quantifiers(1 << (core::log2(sets.sets[set_id].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[set_id].arity; j++) {
			quantifiers[quantifiers.length++] = set_formula;
			set_formula = set_formula->quantifier.operand;
		}

		/* if this set is empty, negations of existentials could now no longer be provable in other set formula */
		if (sets.sets[set_id].set_size == 0)
			return check_set_membership_after_subtraction(set_formula, set_id, std::forward<Args>(visitor)...);

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
					array<variable_assignment> temp_possible_values(1);
					variable_assignment& values = temp_possible_values[0];
					if (!::init(values, sets.sets[set_id].arity))
						return false;
					for (unsigned int k = 0; k < values.assignment.length; k++) {
						core::free(values.assignment.values[k]);
						if (element[k].type == tuple_element_type::CONSTANT) {
							values.assignment.values[k].type = instantiation_type::CONSTANT;
							values.assignment.values[k].constant = element[k].constant;
						} else if (element[k].type == tuple_element_type::NUMBER) {
							values.assignment.values[k].type = instantiation_type::NUMBER;
							values.assignment.values[k].number = element[k].number;
						} else if (element[k].type == tuple_element_type::STRING) {
							values.assignment.values[k].type = instantiation_type::STRING;
							if (!core::init(values.assignment.values[k].str, element[k].str)) {
								values.assignment.values[k].type = instantiation_type::CONSTANT;
								core::free(values); return false;
							}
						}
					}
					temp_possible_values.length++;

					default_prover prover(sets, implication_axioms);
					prover.h.set_ids[0] = current;
					prover.h.set_ids.length = 1;
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
				if (!sets.add_element(new_set, {element, sets.sets[set_id].arity}))
					return false;
			}
			sets.remove_element_at(set_id, i--);
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_new_subset_membership(unsigned int antecedent_set, unsigned int consequent_set, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		/* this new extensional edge may cause the elements in the strongly
		   connected component containing `antecedent_set` may be moveable into
		   a descendant of `consequent_set` */
		hash_set<unsigned int> antecedent_component(8);
		if (!sets.get_strongly_connected_component(antecedent_set, antecedent_component, make_pair(antecedent_set, consequent_set)))
			return false;
		if (antecedent_component.contains(consequent_set))
			return true;
		/* check if any of the elements in `antecedent_component` is provably an element of `antecedent_set` */
		array<variable_assignment> possible_values(8);
		for (unsigned int set_id : antecedent_component) {
			if (!possible_values.ensure_capacity(possible_values.length + sets.sets[set_id].provable_elements.length)) {
				for (auto& element : possible_values) core::free(element);
				return false;
			}
			for (unsigned int i = 0; i < sets.sets[set_id].provable_elements.length; i++) {
				tuple& current_element = sets.sets[set_id].provable_elements[i];

				variable_assignment& values = possible_values[possible_values.length];
				if (!::init(values, sets.sets[set_id].arity)) {
					for (auto& element : possible_values) core::free(element);
					return false;
				}
				for (unsigned int k = 0; k < values.assignment.length; k++) {
					core::free(values.assignment.values[k]);
					if (current_element[k].type == tuple_element_type::CONSTANT) {
						values.assignment.values[k].type = instantiation_type::CONSTANT;
						values.assignment.values[k].constant = current_element[k].constant;
					} else if (current_element[k].type == tuple_element_type::NUMBER) {
						values.assignment.values[k].type = instantiation_type::NUMBER;
						values.assignment.values[k].number = current_element[k].number;
					} else if (current_element[k].type == tuple_element_type::STRING) {
						values.assignment.values[k].type = instantiation_type::STRING;
						if (!core::init(values.assignment.values[k].str, current_element[k].str)) {
							values.assignment.values[k].type = instantiation_type::CONSTANT;
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
		default_prover prover(sets, implication_axioms);
		prover.h.set_ids[0] = antecedent_set;
		prover.h.set_ids.length = 1;
		if (!is_provable_without_abduction<false>(antecedent_formula, quantifiers, possible_values, prover))
			return true;

		array<Formula*> new_atoms(4);
		array_map<tuple, array<unsigned int>> new_elements(4);
		array<unsigned int> new_antecedents(4);
		array<unsigned int> visited_sets(1);
		if (!check_set_membership<ResolveInconsistencies, false, true>(consequent_set, possible_values, new_antecedents, new_elements, quantifiers, new_atoms, 0, visited_sets)) {
			for (auto& element : possible_values) core::free(element);
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			free_all(new_atoms); return false;
		}
		for (auto& element : possible_values) core::free(element);
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }

		bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, new_antecedents, std::forward<Args>(visitor)...);
		free_all(new_atoms);
		return result;
	}

	template<typename... Args>
	bool check_old_subset_membership(unsigned int antecedent_set, unsigned int consequent_set, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		Formula* consequent_formula = sets.sets[consequent_set].size_axioms[0]->formula->binary.left->binary.right;
		array<Formula*> consequent_quantifiers(1 << (core::log2(sets.sets[consequent_set].arity) + 1));
		for (unsigned int j = 0; j < sets.sets[consequent_set].arity; j++) {
			consequent_quantifiers[consequent_quantifiers.length++] = consequent_formula;
			consequent_formula = consequent_formula->quantifier.operand;
		}

		for (unsigned int i = 0; i < sets.sets[antecedent_set].element_count(); i++) {
			/* check each element of `antecedent_set` that is no longer in `provable_elements` of `consequent_set` */
			tuple_element* element_src = sets.sets[antecedent_set].elements.data + i * sets.sets[antecedent_set].arity;
			tuple element;
			element.elements = element_src;
			element.length = sets.sets[antecedent_set].arity;
			if (sets.sets[consequent_set].provable_elements.contains(element))
				continue;

			array<variable_assignment> temp_possible_values(1);
			variable_assignment& values = temp_possible_values[0];
			if (!::init(values, sets.sets[consequent_set].arity))
				return false;
			for (unsigned int k = 0; k < values.assignment.length; k++) {
				core::free(values.assignment.values[k]);
				if (element[k].type == tuple_element_type::CONSTANT) {
					values.assignment.values[k].type = instantiation_type::CONSTANT;
					values.assignment.values[k].constant = element[k].constant;
				} else if (element[k].type == tuple_element_type::NUMBER) {
					values.assignment.values[k].type = instantiation_type::NUMBER;
					values.assignment.values[k].number = element[k].number;
				} else if (element[k].type == tuple_element_type::STRING) {
					values.assignment.values[k].type = instantiation_type::STRING;
					if (!core::init(values.assignment.values[k].str, element[k].str)) {
						values.assignment.values[k].type = instantiation_type::CONSTANT;
						core::free(values); return false;
					}
				}
			}
			temp_possible_values.length++;

			default_prover prover(sets, implication_axioms);
			if (is_provable_without_abduction<false>(consequent_formula, consequent_quantifiers, temp_possible_values, prover)) {
				/* this element is still provably a member of the consequent set */
				for (auto& element : temp_possible_values) core::free(element);
				if (!sets.add_element(consequent_set, element))
					return false;
			}
		}

		array<Formula*> old_atoms(8);
		for (unsigned int descendant : sets.sets[antecedent_set].descendants) {
			for (unsigned int i = sets.sets[descendant].provable_elements.length; i > 0; i--) {
				const tuple& element = sets.sets[descendant].provable_elements[i - 1];
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
				} else {
					if (!add_all_atoms<true>(old_atoms, substituted_formula, 0)) {
						core::free(*substituted_formula); if (substituted_formula->reference_count == 0) core::free(substituted_formula);
						free_all(old_atoms); return false;
					}
					core::free(*substituted_formula); if (substituted_formula->reference_count == 0) core::free(substituted_formula);
				}
			}
		}

		bool result = check_set_membership_after_subtraction(old_atoms, 0, std::forward<Args>(visitor)...);
		free_all(old_atoms);
		return result;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool check_new_implication_satisfaction(unsigned int implication_id, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		Proof* axiom = implication_axioms[implication_id].key;
		Formula* antecedent = axiom->formula->binary.left;

		array<Formula*> quantifiers(1);
		default_prover prover(sets, implication_axioms);
		prover.h.set_ids[0] = 0;
		prover.h.set_ids.length = 1;
		array<variable_assignment> possible_values(1);
		if (!::init(possible_values[0], 0))
			return false;
		possible_values.length++;
		if (!is_provable_without_abduction<false>(antecedent, quantifiers, possible_values, prover))
			return true;
		for (auto& element : possible_values) core::free(element);

		array<Formula*> new_atoms(4);
		array<unsigned int> new_antecedents(4);
		array_map<tuple, array<unsigned int>> new_elements(4);
		array<unsigned int> visited_sets(1);
		if (!check_antecedent_satisfaction(implication_id, new_antecedents, new_elements, new_atoms, 0, visited_sets)) {
			for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
			free_all(new_atoms); return false;
		}
		for (auto entry : new_elements) { core::free(entry.key); core::free(entry.value); }
		bool result = check_set_membership_after_addition<ResolveInconsistencies>(new_atoms, new_antecedents, std::forward<Args>(visitor)...);
		free_all(new_atoms);
		return result;
	}

	template<typename... Args>
	bool check_old_implication_satisfaction(unsigned int implication_id, Args&&... visitor)
	{
time_aggregator profiler(consistency_checking_ms, consistency_checking);
		if (!implication_axioms[implication_id].value)
			return true;
		implication_axioms[implication_id].value = false;
		Proof* axiom = implication_axioms[implication_id].key;
		Formula* consequent = axiom->formula->binary.left;

		array<Formula*> old_atoms(8);
		if (!add_all_atoms<true>(old_atoms, consequent, 0))
			return false;

		bool result = check_set_membership_after_subtraction(old_atoms, 0, std::forward<Args>(visitor)...);
		free_all(old_atoms);
		return result;
	}

	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_proof(Formula* canonicalized, array_map<unsigned int, Proof*>& set_definitions, array_map<unsigned int, unsigned int>& requested_set_sizes, set_changes<Formula>& set_diff, unsigned int& new_constant, Args&&... args)
	{
		if (is_impossible<Contradiction>(canonicalized, *this, std::forward<Args>(args)...))
			return nullptr;

		Term* predicate; Term* arg1; Term* arg2;
		if (is_atomic(*canonicalized, predicate, arg1, arg2)) {
			if (predicate->type == TermType::CONSTANT && predicate->constant >= new_constant_offset) {
				/* find a set definition of the predicate constant */
				unsigned int index = set_definitions.last_index_of(predicate->constant);
				if (index != static_cast<unsigned int>(-1)) {
					Proof* set_definition = set_definitions.values[index];
					Formula* lambda_formula = set_definition->formula->binary.right;

					uint_fast8_t arity = 0;
					Formula* set_formula = lambda_formula;
					while (set_formula->type == FormulaType::LAMBDA) {
						set_formula = set_formula->quantifier.operand;
						arity++;
					}
					if (arg2 == nullptr && arity != 1) return nullptr;
					if (arg2 != nullptr && arity != 2) return nullptr;

					Term* src_variables[2];
					Term* dst_variables[2];
					src_variables[0] = &Variables<1>::value;
					dst_variables[0] = arg1;
					if (arg2 != nullptr) {
						src_variables[1] = &Variables<2>::value;
						dst_variables[1] = arg2;
					}

					Formula* substituted = substitute_all(set_formula, src_variables, dst_variables, arity);
					if (substituted == nullptr)
						return nullptr;

					array_map<unsigned int, Proof*> new_set_definitions(4);
					Proof* new_proof = make_proof<Contradiction, DefinitionsAllowed, ResolveInconsistencies>(substituted, new_set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
					if (new_proof == nullptr) {
						core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
						return nullptr;
					}

					Formula* beta_right = Formula::new_apply(lambda_formula, arg1);
					if (beta_right == nullptr) {
						free_proof(new_proof, set_diff, std::forward<Args>(args)...);
						core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
						return nullptr;
					}
					lambda_formula->reference_count++;
					arg1->reference_count++;
					if (arg2 != nullptr) {
						Formula* temp = Formula::new_apply(beta_right, arg2);
						if (temp == nullptr) {
							core::free(*beta_right); core::free(beta_right);
							free_proof(new_proof, set_diff, std::forward<Args>(args)...);
							core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
							return nullptr;
						}
						arg2->reference_count++;
						beta_right = temp;
					}

					Formula* shifted = shift_bound_variables(substituted, -((int) arity));
					core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
					if (shifted == nullptr) {
						core::free(*beta_right); core::free(beta_right);
						free_proof(new_proof, set_diff, std::forward<Args>(args)...);
						core::free(*substituted); if (substituted->reference_count == 0) core::free(substituted);
						return nullptr;
					}

					Proof* proof = ProofCalculus::new_equality_elim(
							set_definition,
							ProofCalculus::new_equality_elim(ProofCalculus::new_beta(shifted, beta_right), new_proof, make_repeated_array_view(0u, 1)),
							make_repeated_array_view(1u, 1));
					core::free(*shifted); if (shifted->reference_count == 0) core::free(shifted);
					core::free(*beta_right); if (beta_right->reference_count == 0) core::free(beta_right);
					if (proof == nullptr) {
						free_proof(new_proof, set_diff, std::forward<Args>(args)...);
						return nullptr;
					}
					core::free(*new_proof); if (new_proof->reference_count == 0) core::free(new_proof);
					proof->reference_count++;
					return proof;
				}
			}

			return make_atom_proof<DefinitionsAllowed, Contradiction, ResolveInconsistencies>(canonicalized, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::NOT) {
			return make_proof<!Contradiction, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->unary.operand, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);

		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			if (Contradiction) {
				Term* constant;
				Proof* exists_not_proof = make_exists_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand, set_definitions, requested_set_sizes, canonicalized->quantifier.variable, set_diff, constant, new_constant, std::forward<Args>(args)...);
				if (exists_not_proof == NULL) return NULL;
				if (constant->type == TermType::CONSTANT && constant->constant >= new_constant_offset
				 && (!existential_intro_nodes.ensure_capacity(existential_intro_nodes.length + 1)
				  || !ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.ensure_capacity(ground_concepts[constant->constant - new_constant_offset].existential_intro_nodes.length + 1)))
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
				existential_intro_nodes[existential_intro_nodes.length].formula = canonicalized;
				existential_intro_nodes[existential_intro_nodes.length].proof = proof;
				if (!array_map_init(existential_intro_nodes[existential_intro_nodes.length].set_definitions, set_definitions.capacity)) {
					free_proof(proof, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				canonicalized->reference_count++;
				for (unsigned int i = 0; i < set_definitions.size; i++) {
					existential_intro_nodes[existential_intro_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
					existential_intro_nodes[existential_intro_nodes.length].set_definitions.values[i] = set_definitions.values[i];
				}
				existential_intro_nodes[existential_intro_nodes.length++].set_definitions.size = set_definitions.size;
				core::free(*exists_not_proof); if (exists_not_proof->reference_count == 0) core::free(exists_not_proof);
				return proof;
			}

			unsigned int variable = canonicalized->quantifier.variable;
			Formula* new_canonicalized = shift_bound_variables(canonicalized, -((int) (variable - 1)));

			unsigned int arity = 1;
			Formula* operand = new_canonicalized->quantifier.operand;
			while (operand->type == FormulaType::FOR_ALL) {
				operand = operand->quantifier.operand;
				arity++;
			}

			if (operand->type == FormulaType::IF_THEN) {
				Proof* new_axiom;
				Formula* left = operand->binary.left;
				Formula* right = operand->binary.right;

				/* check if `left` or `right` are of the form `c(x)` and that `c` can be a set */
				required_set_size antecedent_set_size = {0, UINT_MAX};
				required_set_size consequent_set_size = {0, UINT_MAX};
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
					unsigned int index = requested_set_sizes.index_of(left->binary.left->constant);
					if (index < requested_set_sizes.size) {
						antecedent_set_size.min_set_size = requested_set_sizes.values[index];
						antecedent_set_size.max_set_size = requested_set_sizes.values[index];
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
					unsigned int index = requested_set_sizes.index_of(right->binary.left->constant);
					if (index < requested_set_sizes.size) {
						consequent_set_size.min_set_size = requested_set_sizes.values[index];
						consequent_set_size.max_set_size = requested_set_sizes.values[index];
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

					unsigned int antecedent_set, consequent_set;
					bool is_antecedent_new, is_consequent_new;
					new_axiom = get_subset_axiom_with_required_set_size<ResolveInconsistencies, false>(
							nullptr, left, right, arity, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new,
							antecedent_set_size, consequent_set_size, set_diff, std::forward<Args>(args)...);
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
					unsigned int antecedent_set, consequent_set;
					bool is_antecedent_new, is_consequent_new;
					new_axiom = get_subset_axiom_with_required_set_size<ResolveInconsistencies, false>(
							nullptr, left, right, arity, antecedent_set, consequent_set, is_antecedent_new, is_consequent_new,
							antecedent_set_size, consequent_set_size, set_diff, std::forward<Args>(args)...);
					core::free(*new_canonicalized); if (new_canonicalized->reference_count == 0) core::free(new_canonicalized);
					if (new_axiom == NULL) return NULL;
					new_axiom->reference_count++;
				}
				return new_axiom;
			}

		} else if (canonicalized->type == FormulaType::AND) {
			if (Contradiction) {
				array<unsigned int> indices(canonicalized->array.length);
				for (unsigned int i = 0; i < canonicalized->array.length; i++)
					indices[i] = i + 1;
				indices.length = canonicalized->array.length;
				shuffle(indices);

				if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

				unsigned int old_set_definition_count = set_definitions.size;
				for (unsigned int index : indices) {
					unsigned int index_minus_one = index - 1;
					Proof* operand = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index_minus_one], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
					if (operand != NULL) {
						if (!negated_conjunction_nodes.ensure_capacity(negated_conjunction_nodes.length + 1)) {
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
							return NULL;
						}

						/* we found a disproof of a conjunct */
						Proof* axiom = ProofCalculus::new_axiom(canonicalized);
						if (axiom == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						}
						axiom->reference_count++;
						Proof* proof = ProofCalculus::new_proof_by_contradiction(ProofCalculus::new_negation_elim(
								ProofCalculus::new_conjunction_elim(axiom, make_array_view(&index_minus_one, 1)), operand), axiom);
						core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
						if (proof == NULL) {
							/* undo the changes made by the recursive call to `make_proof` */
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						}
						proof->reference_count++;

						/* record the formula with this negated conjunction node */
						Formula* negated_conjunction = Formula::new_not(canonicalized);
						if (negated_conjunction == NULL) {
							core::free(*proof); core::free(proof);
							free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
							return NULL;
						}
						canonicalized->reference_count++;
						core::free(*operand); if (operand->reference_count == 0) core::free(operand);

						negated_conjunction_nodes[negated_conjunction_nodes.length].formula = negated_conjunction;
						negated_conjunction_nodes[negated_conjunction_nodes.length].proof = proof;
						if (!array_map_init(negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions, set_definitions.capacity)) {
							core::free(*proof); core::free(proof);
							core::free(*negated_conjunction); core::free(negated_conjunction);
							return NULL;
						}
						for (unsigned int i = 0; i < set_definitions.size; i++) {
							negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
							negated_conjunction_nodes[negated_conjunction_nodes.length].set_definitions.values[i] = set_definitions.values[i];
						}
						negated_conjunction_nodes[negated_conjunction_nodes.length++].set_definitions.size = set_definitions.size;
						return proof;
					}

					set_definitions.size = old_set_definition_count;
					if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
				}

				/* we couldn't find a disproof of any of the conjuncts */
				finished_operand_indices(canonicalized, std::forward<Args>(args)...);
				return NULL;
			}

			/* check if there are any set size statements in the conjunction */
			unsigned int new_requested_set_sizes = 0;
			for (unsigned int i = 1; i < canonicalized->array.length; i++) {
				Term* conjunct = canonicalized->array.operands[i];
				if (conjunct->type != TermType::EQUALS) continue;
				Term* left = conjunct->binary.left;
				Term* right = conjunct->binary.right;
				if (left->type != TermType::UNARY_APPLICATION)
					swap(left, right);
				if (left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT
				 && left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
				 && left->binary.right->type == TermType::CONSTANT && right->type == TermType::NUMBER)
				{
					if (!requested_set_sizes.ensure_capacity(requested_set_sizes.size + 1)
					 || right->number.integer < 0 || right->number.decimal != 0)
						return nullptr;
					unsigned int index = requested_set_sizes.index_of(left->binary.right->constant);
					if (index < requested_set_sizes.size) {
						if (requested_set_sizes.values[index] != right->number.integer)
							return nullptr;
					} else {
						requested_set_sizes.keys[index] = left->binary.right->constant;
						requested_set_sizes.values[index] = right->number.integer;
						requested_set_sizes.size++;
						new_requested_set_sizes++;
					}
				}
			}

			Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
			if (operands == NULL) {
				fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
				return NULL;
			}
			for (unsigned int i = 0; i < canonicalized->array.length; i++) {
				if (new_constant == 0)
					operands[i] = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				else operands[i] = make_proof<false, false, ResolveInconsistencies>(canonicalized->array.operands[i], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);

				if (operands[i] == NULL) {
					requested_set_sizes.size -= new_requested_set_sizes;
					for (unsigned int j = i; j > 0; j--)
						/* undo the changes made by the recursive calls to `make_proof` */
						free_proof(operands[j - 1], set_diff, std::forward<Args>(args)...);
					core::free(operands); return NULL;
				}
			}
			requested_set_sizes.size -= new_requested_set_sizes;
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
				/* check if there are any set size statements in the conjunction */
				unsigned int new_requested_set_sizes = 0;
				for (unsigned int i = 1; i < canonicalized->array.length; i++) {
					Term* conjunct = canonicalized->array.operands[i];
					if (conjunct->type != TermType::EQUALS) continue;
					Term* left = conjunct->binary.left;
					Term* right = conjunct->binary.right;
					if (left->type != TermType::UNARY_APPLICATION)
						swap(left, right);
					if (left->type == TermType::UNARY_APPLICATION && left->binary.left->type == TermType::CONSTANT
					 && left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
					 && left->binary.right->type == TermType::CONSTANT && right->type == TermType::NUMBER)
					{
						if (!requested_set_sizes.ensure_capacity(requested_set_sizes.size + 1)
						 || right->number.integer < 0 || right->number.decimal != 0)
							return nullptr;
						unsigned int index = requested_set_sizes.index_of(left->binary.right->constant);
						if (index < requested_set_sizes.size) {
							if (requested_set_sizes.values[index] != right->number.integer)
								return nullptr;
						} else {
							requested_set_sizes.keys[index] = left->binary.right->constant;
							requested_set_sizes.values[index] = right->number.integer;
							requested_set_sizes.size++;
							new_requested_set_sizes++;
						}
					}
				}

				Proof** operands = (Proof**) malloc(sizeof(Proof*) * canonicalized->array.length);
				if (operands == NULL) {
					fprintf(stderr, "theory.make_proof ERROR: Out of memory.\n");
					return NULL;
				}
				for (unsigned int i = 0; i < canonicalized->array.length; i++) {
					if (new_constant == 0)
						operands[i] = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[i], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
					else operands[i] = make_proof<true, false, ResolveInconsistencies>(canonicalized->array.operands[i], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);

					if (operands[i] == NULL) {
						requested_set_sizes.size -= new_requested_set_sizes;
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
							requested_set_sizes.size -= new_requested_set_sizes;
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
				requested_set_sizes.size -= new_requested_set_sizes;

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

			array<unsigned int> indices(canonicalized->array.length);
			for (unsigned int i = 0; i < canonicalized->array.length; i++)
				indices[i] = i + 1;
			indices.length = canonicalized->array.length;
			shuffle(indices);

			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;

			unsigned int old_set_definition_count = set_definitions.size;
			for (unsigned int index : indices) {
				Proof* operand = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->array.operands[index - 1], set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				if (operand != NULL) {
					if (!disjunction_intro_nodes.ensure_capacity(disjunction_intro_nodes.length + 1)) {
						free_proof(operand, set_diff, std::forward<Args>(args)...); return NULL;
						return NULL;
					}

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
					disjunction_intro_nodes[disjunction_intro_nodes.length].formula = canonicalized;
					disjunction_intro_nodes[disjunction_intro_nodes.length].proof = proof;
					if (!array_map_init(disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions, set_definitions.capacity)) {
						free_proof(proof, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					for (unsigned int i = 0; i < set_definitions.size; i++) {
						disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
						disjunction_intro_nodes[disjunction_intro_nodes.length].set_definitions.values[i] = set_definitions.values[i];
					}
					disjunction_intro_nodes[disjunction_intro_nodes.length++].set_definitions.size = set_definitions.size;
					canonicalized->reference_count++;
					return proof;
				}

				set_definitions.size = old_set_definition_count;
				if (!inconsistent_constant(canonicalized, index, std::forward<Args>(args)...)) return NULL;
			}

			/* we couldn't find a proof of any of the disjuncts */
			finished_operand_indices(canonicalized, std::forward<Args>(args)...);
			return NULL;

		} else if (canonicalized->type == FormulaType::IF_THEN) {
			if (Contradiction) {
				Proof* left = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				if (left == NULL) return NULL;

				Proof* right;
				if (new_constant == 0 && DefinitionsAllowed)
					right = make_proof<true, true, ResolveInconsistencies>(canonicalized->binary.right, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				else right = make_proof<true, false, ResolveInconsistencies>(canonicalized->binary.right, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);

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

			if (ProofCalculus::Intuitionistic) {
				for (unsigned int i = 0; i < implication_axioms.length; i++) {
					if (*implication_axioms[i].key->formula == *canonicalized) {
						implication_axioms[i].key->reference_count++;
						return implication_axioms[i].key;
					}
				}

				if (!implication_axioms.ensure_capacity(implication_axioms.length + 1))
					return nullptr;
				Proof* axiom = ProofCalculus::new_axiom(canonicalized);
				if (axiom == nullptr)
					return nullptr;
				implication_axioms[implication_axioms.length].key = axiom;
				implication_axioms[implication_axioms.length++].value = false;
				axiom->reference_count += 2;
				if (!check_new_implication_satisfaction<ResolveInconsistencies>(implication_axioms.length - 1, set_diff, std::forward<Args>(args)...)) {
					implication_axioms.length--;
					axiom->reference_count--;
					core::free(*axiom); core::free(axiom);
					return nullptr;
				}
				return axiom;
			}

			array<unsigned int> indices(2);
			indices[0] = 1; indices[1] = 2;
			indices.length = 2;
			if (sample_uniform(2) == 1)
				swap(indices[0], indices[1]);
			if (!filter_operands(canonicalized, indices, std::forward<Args>(args)...)) return NULL;
			unsigned int old_set_definition_count = set_definitions.size;
			for (unsigned int i = 0; i < indices.length; i++) {
				if (indices[i] == 1) {
					Proof* left = make_proof<true, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.left, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
					if (left == NULL) {
						set_definitions.size = old_set_definition_count;
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					if (!implication_intro_nodes.ensure_capacity(implication_intro_nodes.length + 1)) {
						free_proof(left, set_diff, std::forward<Args>(args)...);
						return NULL;
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
					implication_intro_nodes[implication_intro_nodes.length].formula = canonicalized;
					implication_intro_nodes[implication_intro_nodes.length].proof = proof;
					if (!array_map_init(implication_intro_nodes[implication_intro_nodes.length].set_definitions, set_definitions.capacity)) {
						free_proof(proof, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					for (unsigned int i = 0; i < set_definitions.size; i++) {
						implication_intro_nodes[implication_intro_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
						implication_intro_nodes[implication_intro_nodes.length].set_definitions.values[i] = set_definitions.values[i];
					}
					implication_intro_nodes[implication_intro_nodes.length++].set_definitions.size = set_definitions.size;
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				} else {
					Proof* right = make_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->binary.right, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
					if (right == NULL) {
						set_definitions.size = old_set_definition_count;
						if (!inconsistent_constant(canonicalized, indices[i], std::forward<Args>(args)...)) return NULL;
						continue;
					}

					if (!implication_intro_nodes.ensure_capacity(implication_intro_nodes.length + 1)) {
						free_proof(right, set_diff, std::forward<Args>(args)...);
						return NULL;
					}

					Proof* proof = ProofCalculus::new_implication_intro(right, ProofCalculus::new_axiom(canonicalized->binary.left));
					if (proof == NULL) {
						free_proof(right, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					core::free(*right); if (right->reference_count == 0) core::free(right);
					/* record the formula with this implication intro node */
					implication_intro_nodes[implication_intro_nodes.length].formula = canonicalized;
					implication_intro_nodes[implication_intro_nodes.length].proof = proof;
					if (!array_map_init(implication_intro_nodes[implication_intro_nodes.length].set_definitions, set_definitions.capacity)) {
						free_proof(proof, set_diff, std::forward<Args>(args)...);
						return NULL;
					}
					for (unsigned int i = 0; i < set_definitions.size; i++) {
						implication_intro_nodes[implication_intro_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
						implication_intro_nodes[implication_intro_nodes.length].set_definitions.values[i] = set_definitions.values[i];
					}
					implication_intro_nodes[implication_intro_nodes.length++].set_definitions.size = set_definitions.size;
					canonicalized->reference_count++;
					proof->reference_count++;
					return proof;
				}
			}

			finished_operand_indices(canonicalized, std::forward<Args>(args)...);
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
						sets.try_free_set(set_id, set_diff, std::forward<Args>(args)...); return NULL;
					}
				}

				Formula* beta_left = shift_bound_variables(lambda_formula, arity);
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
				if (proof == NULL) {
					if (is_set_new) {
						check_old_set_membership(set_id, std::forward<Args>(args)...);
					}
					sets.try_free_set(set_id, std::forward<Args>(args)...); return NULL;
				}
				proof->reference_count++;
				return proof;
			} else {
				Term* variable = Formula::new_variable(canonicalized->quantifier.variable);
				if (variable == NULL) return NULL;

				Term* constant;
				Proof* operand = make_exists_proof<false, DefinitionsAllowed, ResolveInconsistencies>(canonicalized->quantifier.operand,
						set_definitions, requested_set_sizes, canonicalized->quantifier.variable, set_diff, constant, new_constant, std::forward<Args>(args)...);
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
				existential_intro_nodes[existential_intro_nodes.length].formula = canonicalized;
				existential_intro_nodes[existential_intro_nodes.length].proof = proof;
				if (!array_map_init(existential_intro_nodes[existential_intro_nodes.length].set_definitions, set_definitions.capacity)) {
					free_proof(proof, set_diff, std::forward<Args>(args)...);
					return NULL;
				}
				for (unsigned int i = 0; i < set_definitions.size; i++) {
					existential_intro_nodes[existential_intro_nodes.length].set_definitions.keys[i] = set_definitions.keys[i];
					existential_intro_nodes[existential_intro_nodes.length].set_definitions.values[i] = set_definitions.values[i];
				}
				existential_intro_nodes[existential_intro_nodes.length++].set_definitions.size = set_definitions.size;
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
						sets.try_free_set(set_id, set_diff, std::forward<Args>(args)...);
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
				if (Contradiction)
					return nullptr;
				Proof* proof = ProofCalculus::new_beta(left, left);
				if (proof == NULL) return NULL;
				proof->reference_count++;
				return proof;

			} else if (left->type == TermType::CONSTANT && right->type == TermType::CONSTANT) {
				if (left->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return nullptr;
					} else {
						new_constant = right->constant;
						Proof* proof = ProofCalculus::new_beta(right, right);
						if (proof == nullptr) return nullptr;
						proof->reference_count++;
						return proof;
					}
				} else if (right->constant == (unsigned int) built_in_predicates::UNKNOWN) {
					new_constant = left->constant;
					Proof* proof = ProofCalculus::new_beta(left, left);
					if (proof == nullptr) return nullptr;
					proof->reference_count++;
					return proof;
				} else {
					if (Contradiction) {
						Term* negated = Term::new_not(canonicalized);
						if (negated == nullptr)
							return nullptr;
						canonicalized->reference_count++;
						Proof* proof = ProofCalculus::new_inequality_introduction(negated);
						core::free(*negated); if (negated->reference_count == 0) core::free(negated);
						if (proof == nullptr) return nullptr;
						proof->reference_count++;
						return proof;
					} else {
						/* this is impossible */
						return nullptr;
					}
				}
			} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT) {
				bool swap_order = (right->type == TermType::CONSTANT && right->constant != (unsigned int) built_in_predicates::UNKNOWN);
				if (swap_order) swap(left, right);

				if (Contradiction) {
					if (right->type == TermType::NUMBER || right->type == TermType::STRING) {
						Term* negated = Term::new_not(canonicalized);
						if (negated == nullptr)
							return nullptr;
						canonicalized->reference_count++;
						Proof* proof = ProofCalculus::new_inequality_introduction(negated);
						core::free(*negated); if (negated->reference_count == 0) core::free(negated);
						if (proof == nullptr) return nullptr;
						proof->reference_count++;
						return proof;
					} else {
						/* TODO: implement this */
						fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
						return NULL;
					}
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
						if (left->constant >= new_constant_offset && ground_concepts[left->constant - new_constant_offset].types.keys != nullptr) {
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
					if (new_right->type == TermType::LAMBDA && new_right->quantifier.operand->type == TermType::UNARY_APPLICATION
					 && new_right->quantifier.operand->binary.left->type == TermType::CONSTANT
					 && new_right->quantifier.operand->binary.left->constant >= new_constant_offset
					 && new_right->quantifier.operand->binary.right->type == TermType::VARIABLE
					 && new_right->quantifier.operand->binary.right->variable == new_right->quantifier.variable)
					{
						if (new_right->quantifier.operand->binary.left->constant != left->constant
						 && ground_concepts[new_right->quantifier.operand->binary.left->constant - new_constant_offset].types.keys != nullptr)
						{
							/* we found a different concept with this definition */
							core::free(*new_right); if (new_right->reference_count == 0) core::free(new_right);
							return NULL;
						}
					} else {
						bool contains;
						unsigned int constant = reverse_definitions.get(*new_right, contains);
						if (contains && constant != left->constant) {
							/* we found a different concept with this definition */
							core::free(*new_right); if (new_right->reference_count == 0) core::free(new_right);
							return NULL;
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

					required_set_size requested_set_size = {0, UINT_MAX};
					unsigned int index = requested_set_sizes.index_of(left->constant);
					if (index < requested_set_sizes.size) {
						requested_set_size.min_set_size = requested_set_sizes.values[index];
						requested_set_size.max_set_size = requested_set_sizes.values[index];
					}
					Proof* definition = add_definition<ResolveInconsistencies>(new_proof, requested_set_size, set_diff, std::forward<Args>(args)...);
					if (definition != new_proof) {
						core::free(*new_proof);
						core::free(new_proof);
					} if (definition == nullptr) {
						return nullptr;
					}

					if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
						if (!set_definitions.ensure_capacity(set_definitions.size + 1)) {
							free_proof(definition, set_diff);
							return NULL;
						}
						set_definitions.keys[set_definitions.size] = left->constant;
						set_definitions.values[set_definitions.size++] = definition;
					}

					if (swap_order) {
						Proof* swapped_proof = ProofCalculus::new_equality_elim(
								definition, ProofCalculus::new_beta(right, right), make_repeated_array_view(4u, 1));
						if (swapped_proof == NULL) {
							free_proof(definition, set_diff);
							return NULL;
						}
						core::free(*definition);
						swapped_proof->reference_count++;
						return swapped_proof;
					} else {
						return definition;
					}
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
				if (left->type == TermType::LAMBDA && left->quantifier.operand->type == TermType::UNARY_APPLICATION
				 && left->quantifier.operand->binary.left->type == TermType::CONSTANT
				 && left->quantifier.operand->binary.left->constant >= new_constant_offset
				 && left->quantifier.operand->binary.right->type == TermType::VARIABLE
				 && left->quantifier.operand->binary.right->variable == left->quantifier.variable)
				{
					if (ground_concepts[left->quantifier.operand->binary.left->constant - new_constant_offset].types.keys != nullptr)
						return NULL;
				} else {
					if (reverse_definitions.table.contains(*left))
						return NULL;
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
			} else if (Contradiction && (left->type == TermType::STRING || left->type == TermType::NUMBER) && (right->type == TermType::STRING || right->type == TermType::NUMBER)) {
				Term* negated = Term::new_not(canonicalized);
				if (negated == nullptr)
					return nullptr;
				canonicalized->reference_count++;
				Proof* proof = ProofCalculus::new_inequality_introduction(negated);
				core::free(*negated); if (negated->reference_count == 0) core::free(negated);
				if (proof == nullptr) return nullptr;
				proof->reference_count++;
				return proof;
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

	template<typename... Args>
	inline bool get_set_id_with_provable_elements(
			Formula* formula, unsigned int arity, unsigned int& set_id, bool& is_new,
			array<variable_assignment>& possible_values, required_set_size set_size, Args&&... visitor)
	{
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized = Canonicalizer::canonicalize(*formula, variable_map);
		if (canonicalized == nullptr)
			return false;

		bool contains; unsigned int bucket;
		if (!sets.set_ids.check_size()) return false;
		set_id = sets.set_ids.get(*canonicalized, contains, bucket);
		if (!contains) {
			/* get the minimum provable set size for this set */
			array<Formula*> quantifiers(1 << (core::log2(arity) + 1));
			for (unsigned int j = 0; j < arity; j++) {
				quantifiers[quantifiers.length] = (Formula*) alloca(sizeof(Formula));
				quantifiers[quantifiers.length]->type = FormulaType::LAMBDA;
				quantifiers[quantifiers.length]->reference_count = 1;
				quantifiers[quantifiers.length]->quantifier.variable = j + 1;
				quantifiers[quantifiers.length]->quantifier.operand = nullptr;
				quantifiers.length++;
			}

			variable_assignment& values = possible_values[0];
			if (!::init(values, arity)) {
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return false;
			}
			possible_values.length++;

			default_prover prover(sets, implication_axioms);
			prover.h.set_ids.length = 0;
			is_provable_without_abduction<false>(formula, quantifiers, possible_values, prover);
			if (possible_values.length > set_size.max_set_size) {
				for (auto& element : possible_values) core::free(element);
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return false;
			}
			set_size.min_set_size = max(set_size.min_set_size, possible_values.length);

			/* initialize the actual set */
			array_multiset<unsigned int> symbols(16);
			if (!get_constants(*canonicalized, symbols)) {
				for (auto& element : possible_values) core::free(element);
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return false;
			}
			sets.symbols_in_formulas.add(symbols);

			if (!::init(sets.set_ids.table.keys[bucket], *canonicalized)) {
				for (auto& element : possible_values) core::free(element);
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				return false;
			} else if (!sets.new_set(canonicalized, arity, set_id, set_size, std::forward<Args>(visitor)...)) {
				for (auto& element : possible_values) core::free(element);
				core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
				core::free(sets.set_ids.table.keys[bucket]);
				core::set_empty(sets.set_ids.table.keys[bucket]);
				return false;
			}
			sets.set_ids.values[bucket] = set_id;
			sets.set_ids.table.size++;
			is_new = true;

			/* if this set is empty, negations of existentials could now be provable in other set formula */
			if (sets.sets[set_id].set_size == 0 && !check_set_membership_after_addition<false>(formula, std::forward<Args>(visitor)...)) {
				sets.try_free_set(set_id);
				return false;
			}
		} else {
			is_new = false;
		}
		core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
		return true;
	}

	template<bool ResolveInconsistencies, bool SameSize, typename... Args>
	inline Proof* get_subset_axiom_with_required_set_size(Proof* axiom,
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int& antecedent_set, unsigned int& consequent_set,
			bool& is_antecedent_new, bool& is_consequent_new,
			required_set_size antecedent_set_size,
			required_set_size consequent_set_size, Args&&... visitor)
	{
		array<variable_assignment> possible_values(4);
		if (!get_set_id_with_provable_elements(antecedent, arity, antecedent_set, is_antecedent_new, possible_values, antecedent_set_size, std::forward<Args>(visitor)...))
			return nullptr;
		if (is_antecedent_new) {
			Formula* set_formula = sets.sets[antecedent_set].size_axioms[0]->formula->binary.left->binary.right;
			array<Formula*> quantifiers(1 << (core::log2(sets.sets[antecedent_set].arity) + 1));
			for (unsigned int j = 0; j < sets.sets[antecedent_set].arity; j++) {
				quantifiers[quantifiers.length++] = set_formula;
				set_formula = set_formula->quantifier.operand;
			}
			if (!check_new_set_membership<ResolveInconsistencies>(antecedent_set, quantifiers, possible_values, std::forward<Args>(visitor)...)) {
				for (auto& element : possible_values) core::free(element);
				sets.try_free_set(antecedent_set);
				return nullptr;
			}
		}
		for (auto& element : possible_values) core::free(element);

		if (SameSize) {
			consequent_set_size.min_set_size = sets.sets[antecedent_set].set_size;
			consequent_set_size.max_set_size = sets.sets[antecedent_set].set_size;
		} else {
			consequent_set_size.min_set_size = max(consequent_set_size.min_set_size, sets.sets[antecedent_set].set_size);
		}
		if (consequent_set_size.min_set_size > consequent_set_size.max_set_size
		 || !sets.get_set_id(consequent, arity, consequent_set, is_consequent_new, consequent_set_size, std::forward<Args>(visitor)...))
		{
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

		if (axiom == nullptr) {
			bool dummy;
			axiom = sets.template get_subset_axiom<ResolveInconsistencies, false>(
					antecedent, consequent, arity,
					antecedent_set, consequent_set,
					dummy, dummy, std::forward<Args>(visitor)...);
		} else {
			if (!sets.add_existing_subset_axiom(axiom, antecedent, consequent, antecedent_set, consequent_set))
				axiom = nullptr;
		}
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
			if (sets.sets[consequent_set].provable_elements.length == sets.sets[consequent_set].set_size && !check_set_membership_after_addition<ResolveInconsistencies>(consequent, std::forward<Args>(visitor)...))
			{
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

		if (is_subset_new) {
			if (!check_new_subset_membership<ResolveInconsistencies>(antecedent_set, consequent_set, std::forward<Args>(visitor)...)) {
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
				if (sets.sets[consequent_set].provable_elements.length == sets.sets[consequent_set].set_size && !check_set_membership_after_addition<ResolveInconsistencies>(consequent, std::forward<Args>(visitor)...))
				{
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
			}
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

	template<bool SameSize, typename... Args>
	inline void free_subset_axiom_with_required_set_size(
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int antecedent_set, unsigned int consequent_set,
			required_set_size antecedent_set_size, required_set_size consequent_set_size,
			Args&&... visitor)
	{
		consequent->reference_count++;
		sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
		if (consequent->type == TermType::UNARY_APPLICATION && consequent->binary.right->type == TermType::VARIABLE && consequent->binary.right->variable == 1) {
			/* removing this edge could decrease the number of provable elements
			   in `consequent_set` to no longer be its set size, which can
			   cause `a(b)` to be newly provable or disprovable */
			unsigned int new_set_size = sets.sets[consequent_set].provable_elements.length;
			if (new_set_size < sets.sets[consequent_set].set_size) {
				unsigned int provable_element_count = sets.count_union_of_provable_elements(consequent_set, antecedent_set);
				if (provable_element_count == sets.sets[consequent_set].set_size
				 && !check_set_membership_after_subtraction(consequent, 0, std::forward<Args>(visitor)...))
				{
					core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
					return;
				}
			}
		}
		core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
		check_old_subset_membership(antecedent_set, consequent_set, std::forward<Args>(visitor)...);
		if (sets.is_freeable(consequent_set, std::forward<Args>(visitor)...)) {
			if (SameSize) {
				consequent_set_size.min_set_size = sets.sets[antecedent_set].set_size;
				consequent_set_size.max_set_size = sets.sets[antecedent_set].set_size;
			} else {
				consequent_set_size.min_set_size = max(consequent_set_size.min_set_size, sets.sets[antecedent_set].set_size);
			}
			for (Proof* size_axiom : sets.sets[consequent_set].size_axioms)
				on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
			unsigned int min_set_size; unsigned int max_set_size;
			sets.set_size_bounds(consequent_set, min_set_size, max_set_size);
			on_free_set(consequent_set, sets, min_set_size, max_set_size, consequent_set_size, std::forward<Args>(visitor)...);
			check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
			sets.free_set(consequent_set);
		} if (sets.is_freeable(antecedent_set, std::forward<Args>(visitor)...)) {
			for (Proof* size_axiom : sets.sets[antecedent_set].size_axioms)
				on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
			unsigned int min_set_size; unsigned int max_set_size;
			sets.set_size_bounds(antecedent_set, min_set_size, max_set_size);
			on_free_set(antecedent_set, sets, min_set_size, max_set_size, antecedent_set_size, std::forward<Args>(visitor)...);
			check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
			sets.free_set(antecedent_set);
		}
	}

	template<typename... Args>
	inline void free_subset_axiom(
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int antecedent_set, unsigned int consequent_set,
			Args&&... visitor)
	{
		consequent->reference_count++;
		sets.template free_subset_axiom<false>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
		if (consequent->type == TermType::UNARY_APPLICATION && consequent->binary.right->type == TermType::VARIABLE && consequent->binary.right->variable == 1) {
			/* removing this edge could decrease the number of provable elements
			   in `consequent_set` to no longer be its set size, which can
			   cause `a(b)` to be newly provable or disprovable */
			unsigned int new_set_size = sets.sets[consequent_set].provable_elements.length;
			if (new_set_size < sets.sets[consequent_set].set_size) {
				unsigned int provable_element_count = sets.count_union_of_provable_elements(consequent_set, antecedent_set);
				if (provable_element_count == sets.sets[consequent_set].set_size
				 && !check_set_membership_after_subtraction(consequent, 0, std::forward<Args>(visitor)...))
				{
					core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
					return;
				}
			}
		}
		core::free(*consequent); if (consequent->reference_count == 0) core::free(consequent);
		check_old_subset_membership(antecedent_set, consequent_set, std::forward<Args>(visitor)...);
		if (sets.is_freeable(consequent_set, std::forward<Args>(visitor)...)) {
			for (Proof* size_axiom : sets.sets[consequent_set].size_axioms)
				on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
			unsigned int min_set_size; unsigned int max_set_size;
			sets.set_size_bounds(consequent_set, min_set_size, max_set_size);
			on_free_set(consequent_set, sets, min_set_size, max_set_size, std::forward<Args>(visitor)...);
			check_old_set_membership(consequent_set, std::forward<Args>(visitor)...);
			sets.free_set(consequent_set);
		} if (sets.is_freeable(antecedent_set, std::forward<Args>(visitor)...)) {
			for (Proof* size_axiom : sets.sets[antecedent_set].size_axioms)
				on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
			unsigned int min_set_size; unsigned int max_set_size;
			sets.set_size_bounds(antecedent_set, min_set_size, max_set_size);
			on_free_set(antecedent_set, sets, min_set_size, max_set_size, std::forward<Args>(visitor)...);
			check_old_set_membership(antecedent_set, std::forward<Args>(visitor)...);
			sets.free_set(antecedent_set);
		}
	}

private:
	template<bool ResolveInconsistencies, typename... Args>
	Proof* add_definition(Proof* definition, required_set_size requested_set_size, Args&&... args)
	{
		if (!reverse_definitions.check_size())
			return nullptr;

		Formula* constant = definition->formula->binary.left;
		Formula* new_definition = definition->formula->binary.right;

		unsigned int arity = 1;
		Formula* new_set_formula = nullptr;
		if (new_definition->type == FormulaType::LAMBDA) {
			new_set_formula = new_definition->quantifier.operand;
			while (new_set_formula->type == FormulaType::LAMBDA) {
				new_set_formula = new_set_formula->quantifier.operand;
				arity++;
			}
			if (!try_init_concept(constant->constant, arity))
				return nullptr;
		} else {
			if (!try_init_concept(constant->constant))
				return nullptr;
		}

		/* check if the axiom already exists */
		if (constant->constant < new_constant_offset) {
			print("ERROR: Attempted to add definition for symbol '", stderr); print(*constant, stderr, *debug_terminal_printer); print("'.\n", stderr);
			exit(EXIT_FAILURE);
		}
		for (Proof* definition : ground_concepts[constant->constant - new_constant_offset].definitions) {
			if ((definition->formula->binary.left == constant || *definition->formula->binary.left == *constant)
			 && (definition->formula->binary.right == new_definition || *definition->formula->binary.right == *new_definition))
			{
				definition->reference_count++;
				return definition;
			}
		}

		/* make sure all the concepts referenced from the right-hand side are initialized */
		array<unsigned int> constants(4);
		constants[constants.length++] = constant->constant;
		if (!get_constants(*new_definition, constants, new_constant_offset)) {
			try_free_concept_id(constant->constant);
			return nullptr;
		}
		for (unsigned int i = 0; i < constants.length; i++) {
			if (!try_init_concept(constants[i])) {
				for (unsigned int j = i; j > 0; j--)
					try_free_concept_id(constants[j - 1]);
				return nullptr;
			}
		}

		if (new_definition->type == FormulaType::LAMBDA) {
			/* check if this constant defines any other sets, and indicate to the set reasoning module that they are the same set */
			bool contains;
			unsigned int set_id = sets.set_ids.get(*new_set_formula, contains);
			if (contains) {
				requested_set_size.min_set_size = sets.sets[set_id].set_size;
				requested_set_size.max_set_size = sets.sets[set_id].set_size;
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
						if (sets.sets[set_id].set_size < requested_set_size.min_set_size
						 || sets.sets[set_id].set_size > requested_set_size.max_set_size)
						{
							for (unsigned int j = constants.length; j > 0; j--)
								try_free_concept_id(constants[j - 1]);
							return nullptr;
						}
						requested_set_size.min_set_size = sets.sets[set_id].set_size;
						requested_set_size.max_set_size = sets.sets[set_id].set_size;
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
					Proof* first_subset_axiom = get_subset_axiom_with_required_set_size<ResolveInconsistencies, true>(
							nullptr, new_set_formula, other_set_formula, arity, antecedent_set, consequent_set,
							is_antecedent_new, is_consequent_new, requested_set_size, requested_set_size, std::forward<Args>(args)...);
					if (first_subset_axiom == nullptr) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Formula* other_set_formula = definition->formula->binary.right;
							if (other_set_formula->type != FormulaType::LAMBDA) continue;
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
								free_subset_axiom_with_required_set_size<true>(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, requested_set_size, requested_set_size, std::forward<Args>(args)...);
						}
						for (unsigned int j = constants.length; j > 0; j--)
							try_free_concept_id(constants[j - 1]);
						return nullptr;
					}
					first_subset_axiom->reference_count++;

					Proof* second_subset_axiom = get_subset_axiom<ResolveInconsistencies>(other_set_formula, new_set_formula, arity,
							antecedent_set, consequent_set, is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
					if (second_subset_axiom == nullptr) {
						/* undo the changes we've made so far */
						on_subtract_changes(std::forward<Args>(args)...);
						core::free(*first_subset_axiom);
						if (first_subset_axiom->reference_count == 1)
							free_subset_axiom_with_required_set_size<true>(new_set_formula, other_set_formula, arity, antecedent_set, consequent_set, requested_set_size, requested_set_size, std::forward<Args>(args)...);
						for (unsigned int j = 0; j < i; j++) {
							Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
							Formula* other_set_formula = definition->formula->binary.right;
							if (other_set_formula->type != FormulaType::LAMBDA) continue;
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
								free_subset_axiom_with_required_set_size<true>(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, requested_set_size, requested_set_size, std::forward<Args>(args)...);
						}
						for (unsigned int j = constants.length; j > 0; j--)
							try_free_concept_id(constants[j - 1]);
						return nullptr;
					}
					second_subset_axiom->reference_count++;
				}
			}
		}

		if (!ground_concepts[constant->constant - new_constant_offset].definitions.add(definition)) {
			/* undo the changes we've made so far */
			on_subtract_changes(std::forward<Args>(args)...);
			if (new_definition->type == FormulaType::LAMBDA) {
				for (unsigned int j = 0; j < ground_concepts[constant->constant - new_constant_offset].definitions.length; j++) {
					unsigned int antecedent_set, consequent_set;
					bool is_antecedent_new, is_consequent_new;
					Proof* definition = ground_concepts[constant->constant - new_constant_offset].definitions[j];
					if (definition->formula->binary.right->type != FormulaType::LAMBDA) continue;
					Formula* other_set_formula = definition->formula->binary.right->quantifier.operand;
					if (other_set_formula->type != FormulaType::LAMBDA) continue;
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
						free_subset_axiom_with_required_set_size<true>(other_set_formula, new_set_formula, arity, antecedent_set, consequent_set, requested_set_size, requested_set_size, std::forward<Args>(args)...);
				}
			}
			for (unsigned int j = constants.length; j > 0; j--)
				try_free_concept_id(constants[j - 1]);
			return nullptr;
		}
		definition->reference_count++;

		unsigned int index = reverse_definitions.table.index_to_insert(*new_definition);
		reverse_definitions.table.keys[index] = *new_definition;
		reverse_definitions.values[index] = constant->constant;
		reverse_definitions.table.size++;

		if (!check_set_membership_after_addition<ResolveInconsistencies>(definition->formula, std::forward<Args>(args)...)) {
			remove_definition(definition, requested_set_size, std::forward<Args>(args)...);
			return nullptr;
		} else if (definition->formula->type == FormulaType::EQUALS) {
			Formula* swapped = Formula::new_equals(definition->formula->binary.right, definition->formula->binary.left);
			if (swapped == nullptr) {
				remove_definition(definition, requested_set_size, std::forward<Args>(args)...);
				return nullptr;
			}
			definition->formula->binary.left->reference_count++;
			definition->formula->binary.right->reference_count++;
			if (!check_set_membership_after_addition<ResolveInconsistencies>(swapped, std::forward<Args>(args)...)) {
				core::free(*swapped); core::free(swapped);
				remove_definition(definition, requested_set_size, std::forward<Args>(args)...);
				return nullptr;
			}
			core::free(*swapped); core::free(swapped);
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
	void remove_definition(Proof* definition, required_set_size requested_set_size, Args&&... args)
	{
		unsigned int concept_id = definition->formula->binary.left->constant;

		unsigned int index;
		for (index = 1; index < ground_concepts[concept_id - new_constant_offset].definitions.length; index++) {
			Proof* other_definition = ground_concepts[concept_id - new_constant_offset].definitions[index];
			if (definition->formula->binary.right == other_definition->formula->binary.right) break;
		}

#if !defined(NDEBUG)
		if (index == ground_concepts[concept_id - new_constant_offset].definitions.length)
			fprintf(stderr, "remove_definition WARNING: Unable to find definition to remove.\n");
#endif

		ground_concepts[concept_id - new_constant_offset].definitions.remove(index);
		check_set_membership_after_subtraction(definition->formula, 0, std::forward<Args>(args)...);

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
				if (axiom->reference_count == 1)
					free_subset_axiom(other_set_formula, set_formula, arity, antecedent_set, consequent_set, std::forward<Args>(args)...);

				axiom = get_subset_axiom<false>(
						set_formula, other_set_formula, arity, antecedent_set, consequent_set,
						is_antecedent_new, is_consequent_new, std::forward<Args>(args)...);
				core::free(*axiom);
				if (axiom->reference_count == 1)
					free_subset_axiom_with_required_set_size<true>(set_formula, other_set_formula, arity, antecedent_set, consequent_set, requested_set_size, requested_set_size, std::forward<Args>(args)...);
			}
		}

		index = reverse_definitions.table.index_of(*definition->formula->binary.right);
		core::free(reverse_definitions.table.keys[index]);
		reverse_definitions.remove_at(index);

		/* make sure all the concepts referenced from the right-hand side are freed if necessary */
		array<unsigned int> constants(4);
		constants[constants.length++] = concept_id;
		get_constants(*definition->formula->binary.right, constants, new_constant_offset);
		for (unsigned int j = constants.length; j > 0; j--)
			try_free_concept_id(constants[j - 1]);
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
		check_set_membership_after_subtraction(function_value_axiom->formula, 0, std::forward<Args>(args)...);
		core::free(*function_value_axiom); if (function_value_axiom->reference_count == 0) core::free(function_value_axiom);
	}

	/* NOTE: this function finds a constant that proves `quantified`, and not ?[x]:`quantified` */
	template<bool Contradiction, bool DefinitionsAllowed, bool ResolveInconsistencies, typename... Args>
	Proof* make_exists_proof(Formula* quantified,
			array_map<unsigned int, Proof*>& set_definitions,
			array_map<unsigned int, unsigned int>& requested_set_sizes,
			unsigned int variable, set_changes<Formula>& set_diff,
			Term*& constant, unsigned int& new_constant, Args&&... args)
	{
		Formula* var = Formula::new_variable(variable);
		if (var == nullptr) return nullptr;

		array<instance> constants(ground_concept_capacity + constant_types.size + constant_negated_types.size + 1);
		if (!get_constants(*this, quantified, variable, constants, set_definitions, std::forward<Args>(args)...)) {
			core::free(*var); if (var->reference_count == 0) core::free(var);
			return nullptr;
		}

		unsigned int old_set_definition_count = set_definitions.size;
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

				Proof* proof = make_proof<Contradiction, false, ResolveInconsistencies>(substituted, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				set_definitions.size = old_set_definition_count;
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

				Proof* proof = make_proof<Contradiction, DefinitionsAllowed, ResolveInconsistencies>(substituted, set_definitions, requested_set_sizes, set_diff, new_constant, std::forward<Args>(args)...);
				set_definitions.size = old_set_definition_count;
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

		finished_constants(quantified, std::forward<Args>(args)...);
		core::free(*var); if (var->reference_count == 0) core::free(var);
		return nullptr;
	}

	template<bool DefinitionsAllowed, bool Negated, bool ResolveInconsistencies, typename... Args>
	Proof* make_atom_proof(
			Term* atom,
			unsigned int& new_constant,
			Args&&... args)
	{
		if (atom->type == TermType::UNARY_APPLICATION && (atom->binary.right->type == TermType::CONSTANT || atom->binary.right->type == TermType::NUMBER))
		{
			/* this is a unary formula */
			Term* arg1 = atom->binary.right;
			if (atom->binary.left->type != TermType::CONSTANT || atom->binary.left->constant != (unsigned int) built_in_predicates::UNKNOWN) {
				if (arg1->type == TermType::CONSTANT && arg1->constant == (unsigned int) built_in_predicates::UNKNOWN) {
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

					array_multiset<unsigned int> constants(8);
					if (!get_constants(*formula, constants)) {
						core::free(*new_axiom); core::free(new_axiom);
						return nullptr;
					}
					for (unsigned int i = 0; i < constants.counts.size; i++) {
						if (!try_init_concept(constants.counts.keys[i])) {
							core::free(*new_axiom); core::free(new_axiom);
							return nullptr;
						}
					}

					if (!add_unary_atom<Negated, ResolveInconsistencies>(Negated ? *formula->unary.operand : *formula, new_axiom, std::forward<Args>(args)...)) {
						core::free(*new_axiom); core::free(new_axiom);
						for (unsigned int i = 0; i < constants.counts.size; i++)
							try_free_concept_id(constants.counts.keys[i]);
						return NULL;
					}
					return new_axiom;
				} else {
					/* this is a formula of form `t(c)` */

					/* we do not allow sets to contain themselves */
					if (!Negated && *atom->binary.left == *arg1)
						return NULL;

					/* make sure that `predicate` could be a set, and that `arg` could be not a set */
					if (!Negated && ((atom->binary.left->type == TermType::CONSTANT && is_provably_not_a_set(atom->binary.left->constant)) || (arg1->type == TermType::CONSTANT && is_provably_a_set(arg1->constant))))
						return NULL;

					/* first check if there is already an extensional edge that proves this formula (for this particular instantiation) */
					Formula* set_formula = (Negated
							? Formula::new_not(Formula::new_apply(atom->binary.left, &Variables<1>::value))
							: Formula::new_apply(atom->binary.left, &Variables<1>::value));
					if (set_formula == nullptr)
						return nullptr;
					atom->binary.left->reference_count++;
					Variables<1>::value.reference_count++;
					bool contains;
					unsigned int set_id = sets.set_ids.get(*set_formula, contains);
					core::free(*set_formula); core::free(set_formula);
					if (contains) {
						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[set_id].children) {
							for (Proof* edge_axiom : entry.value) {
								for (Proof* child : edge_axiom->children) {
									if (child->type == ProofType::UNIVERSAL_ELIMINATION
									 && child->operands[1]->type == ProofType::TERM_PARAMETER
									 && *child->operands[1]->term == *arg1)
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
					set_formula = (Negated
							? Formula::new_not(Formula::new_apply(&Variables<1>::value, arg1))
							: Formula::new_apply(&Variables<1>::value, arg1));
					if (set_formula == nullptr)
						return nullptr;
					arg1->reference_count++;
					Variables<1>::value.reference_count++;
					set_id = sets.set_ids.get(*set_formula, contains);
					core::free(*set_formula); core::free(set_formula);
					if (contains) {
						/* the atomic formula unifies with `set_formula` of this set */
						for (auto entry : sets.extensional_graph.vertices[set_id].children) {
							for (Proof* edge_axiom : entry.value) {
								for (Proof* child : edge_axiom->children) {
									if (child->type == ProofType::UNIVERSAL_ELIMINATION
									 && child->operands[1]->type == ProofType::TERM_PARAMETER
									 && *child->operands[1]->term == *atom->binary.left)
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

					/* make sure we can't disprove this atom */
					array<Formula*> quantifiers(1);
					array<variable_assignment> temp_possible_values(1);
					if (!::init(temp_possible_values[0], 0))
						return nullptr;
					temp_possible_values.length = 1;
					default_prover prover(sets, implication_axioms);
					prover.h.set_ids[0] = 0;
					prover.h.set_ids.length = 1;
					if (is_provable_without_abduction<!Negated>(atom, quantifiers, temp_possible_values, prover)) {
						for (auto& element : temp_possible_values) core::free(element);
						return nullptr;
					}

					array_multiset<unsigned int> constants(8);
					if (!get_constants(*atom, constants)) return nullptr;
					for (unsigned int i = 0; i < constants.counts.size; i++)
						if (!try_init_concept(constants.counts.keys[i])) return nullptr;

					/* there is no extensional edge that proves this atomic formula,
					   so check if the axiom already exists, and if not, create a new axiom */
					if (arg1->type == TermType::CONSTANT) {
						Term* lifted_atom = Term::new_apply(atom->binary.left, &Variables<1>::value);
						if (lifted_atom == nullptr) return nullptr;
						atom->binary.left->reference_count++;
						Variables<1>::value.reference_count++;

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
					} else {
						Proof* axiom = (Negated ?
								constant_negated_types.get(*atom, contains) :
								constant_types.get(*atom, contains));
						if (contains) {
							axiom->reference_count++;
							return axiom;
						}
					}

					Formula* formula = Negated ? Formula::new_not(atom) : atom;
					if (formula == NULL) {
						for (unsigned int i = 0; i < constants.counts.size; i++)
							try_free_concept_id(constants.counts.keys[i]);
						return NULL;
					}
					atom->reference_count++;
					Proof* new_axiom = ProofCalculus::new_axiom(formula);
					core::free(*formula); if (formula->reference_count == 0) core::free(formula);
					if (new_axiom == NULL) {
						for (unsigned int i = 0; i < constants.counts.size; i++)
							try_free_concept_id(constants.counts.keys[i]);
						return NULL;
					}
					new_axiom->reference_count++;
					if (!add_unary_atom<Negated, ResolveInconsistencies>(*atom, new_axiom, std::forward<Args>(args)...)) {
						core::free(*new_axiom); core::free(new_axiom);
						for (unsigned int i = 0; i < constants.counts.size; i++)
							try_free_concept_id(constants.counts.keys[i]);
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

template<typename ProofCalculus, typename Canonicalizer>
bool get_proof_map(const theory<ProofCalculus, Canonicalizer>& T,
		hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	if (!get_proof_map(T.empty_set_axiom, proof_map, formula_map)
	 || !get_formula_map(T.NAME_ATOM, formula_map)
	 || !get_proof_map(T.sets, proof_map, formula_map)
	 || !get_formula_map(T.ctx, formula_map))
	{
		return false;
	}

	for (const Proof* axiom : T.built_in_axioms)
		if (!get_proof_map(axiom, proof_map, formula_map)) return false;

	for (const auto& entry : T.atoms) {
		if (!get_formula_map(entry.key, formula_map))
			return false;
	} for (const Proof* observation : T.observations) {
		if (!get_proof_map(observation, proof_map, formula_map))
			return false;
	} for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;
		if (!get_proof_map(T.ground_concepts[i], proof_map, formula_map))
			return false;
	} for (const auto& entry : T.constant_types) {
		if (!get_formula_map(entry.key, formula_map)
		 || !get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.constant_negated_types) {
		if (!get_formula_map(entry.key, formula_map)
		 || !get_proof_map(entry.value, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.disjunction_intro_nodes) {
		if (!get_formula_map(entry.formula, formula_map)
		 || !get_proof_map(entry.proof, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.negated_conjunction_nodes) {
		if (!get_formula_map(entry.formula, formula_map)
		 || !get_proof_map(entry.proof, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.implication_intro_nodes) {
		if (!get_formula_map(entry.formula, formula_map)
		 || !get_proof_map(entry.proof, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.existential_intro_nodes) {
		if (!get_formula_map(entry.formula, formula_map)
		 || !get_proof_map(entry.proof, proof_map, formula_map))
			return false;
	} for (const auto& entry : T.implication_axioms) {
		if (!get_proof_map(entry.key, proof_map, formula_map))
			return false;
	}
	return true;
}

template<typename ProofCalculus, typename Canonicalizer, typename Stream>
bool read(theory<ProofCalculus, Canonicalizer>& T, Stream& in,
		typename ProofCalculus::Proof** proofs,
		typename ProofCalculus::Language** formulas)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;

	unsigned int empty_set_axiom_index;
	unsigned int NAME_ATOM_index;
	decltype(T.atoms.table.size) atom_count;
	decltype(T.built_in_axioms.length) built_in_axiom_count;
	if (!read(T.new_constant_offset, in)
	 || !read(T.ground_concept_capacity, in)
	 || !read(T.ground_axiom_count, in)
	 || !read(built_in_axiom_count, in)
	 || !read(empty_set_axiom_index, in)
	 || !read(NAME_ATOM_index, in)
	 || !read(T.sets, in, proofs, formulas)
	 || !read(T.ctx, in, formulas))
	{
		return false;
	} else if (!read(atom_count, in)) {
		free(T.sets);
		return false;
	}

	T.empty_set_axiom = proofs[empty_set_axiom_index];
	T.empty_set_axiom->reference_count++;
	T.NAME_ATOM = formulas[NAME_ATOM_index];
	T.NAME_ATOM->reference_count++;

	if (!array_init(T.built_in_axioms, ((size_t) 1) << (core::log2(built_in_axiom_count == 0 ? 1 : built_in_axiom_count) + 1))) {
		free(*T.empty_set_axiom);
		free(*T.NAME_ATOM);
		free(T.sets); return false;
	}
	for (unsigned int i = 0; i < built_in_axiom_count; i++) {
		unsigned int index;
		if (!read(index, in)) {
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom);
			free(*T.NAME_ATOM);
			free(T.sets); return false;
		}
		T.built_in_axioms[T.built_in_axioms.length++] = proofs[index];
		proofs[index]->reference_count++;
	}

	if (!hash_map_init(T.atoms, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (atom_count == 0 ? 1 : atom_count)) + 1))) {
		free(*T.empty_set_axiom);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.NAME_ATOM);
		free(T.sets); return false;
	}
	Term& term = *((Term*) alloca(sizeof(Term)));
	for (unsigned int i = 0; i < atom_count; i++) {
		if (!read(term, in, formulas)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); return false;
		}
		term.reference_count = 1;
		unsigned int bucket = T.atoms.table.index_to_insert(term);
		if (!read(T.atoms.values[bucket].key, in)) {
			free(term);
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); return false;
		} else if (!read(T.atoms.values[bucket].value, in)) {
			free(term); free(T.atoms.values[bucket].key);
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); return false;
		}
		move(term, T.atoms.table.keys[bucket]);
		T.atoms.table.size++;
	}

	if (!read(atom_count, in)
	 || !hash_map_init(T.relations, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (atom_count == 0 ? 1 : atom_count)) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		}
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); return false;
	}
	for (unsigned int i = 0; i < atom_count; i++) {
		relation rel;
		if (!read(rel, in)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			return false;
		}
		unsigned int bucket = T.relations.table.index_to_insert(rel);
		if (!read(T.relations.values[bucket].key, in)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			return false;
		} else if (!read(T.relations.values[bucket].value, in)) {
			free(T.relations.values[bucket].key);
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			return false;
		}
		move(rel, T.relations.table.keys[bucket]);
		T.relations.table.size++;
	}

	decltype(T.observations.size) observation_count;
	if (!read(observation_count, in)
	 || !hash_set_init(T.observations, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (observation_count == 0 ? 1 : observation_count)) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		}
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		return false;
	}
	unsigned int index;
	for (unsigned int i = 0; i < observation_count; i++) {
		if (!read(index, in)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); return false;
		}
		Proof* observation = proofs[index];
		unsigned int bucket = T.observations.index_to_insert(observation);
		T.observations.keys[bucket] = observation;
		observation->reference_count++;
		T.observations.size++;
	}

	unsigned int ground_concept_count;
	T.ground_concepts = (concept<ProofCalculus>*) calloc(T.ground_concept_capacity, sizeof(concept<ProofCalculus>));
	if (T.ground_concepts == nullptr || !read(ground_concept_count, in)) {
		if (T.ground_concepts != nullptr) free(T.ground_concepts);
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); return false;
	}
	for (unsigned int i = 0; i < ground_concept_count; i++) {
		if (!read(index, in)
		 || !read(T.ground_concepts[index], in, proofs, formulas))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < i; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			return false;
		}
	}

	decltype(T.constant_types.size) constant_type_count;
	if (!read(constant_type_count, in)
	 || !array_map_init(T.constant_types, ((size_t) 1) << (core::log2(constant_type_count == 0 ? 1 : constant_type_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		return false;
	}
	for (unsigned int i = 0; i < constant_type_count; i++) {
		if (!read(T.constant_types.keys[i], in, formulas)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); return false;
		}
		T.constant_types.keys[i].reference_count = 1;
		if (!read(index, in)) {
			free(T.constant_types.keys[i]);
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); return false;
		}
		T.constant_types.values[i] = proofs[index];
		proofs[index]->reference_count++;
		T.constant_types.size++;
	}

	if (!read(constant_type_count, in)
	 || !array_map_init(T.constant_negated_types, ((size_t) 1) << (core::log2(constant_type_count == 0 ? 1 : constant_type_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); return false;
	}
	for (unsigned int i = 0; i < constant_type_count; i++) {
		if (!read(T.constant_negated_types.keys[i], in, formulas)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			return false;
		}
		T.constant_negated_types.keys[i].reference_count = 1;
		if (!read(index, in)) {
			free(T.constant_negated_types.keys[i]);
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			return false;
		}
		T.constant_negated_types.values[i] = proofs[index];
		proofs[index]->reference_count++;
		T.constant_negated_types.size++;
	}

	decltype(T.disjunction_intro_nodes.length) proof_count;
	if (!read(proof_count, in)
	 || !array_init(T.disjunction_intro_nodes, ((size_t) 1) << (core::log2(proof_count == 0 ? 1 : proof_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		return false;
	}
	for (unsigned int i = 0; i < proof_count; i++) {
		unsigned int formula_index;
		unsigned int proof_index;
		decltype(T.disjunction_intro_nodes[i].set_definitions.size) set_definition_count;
		if (!read(formula_index, in)
		 || !read(proof_index, in)
		 || !read(set_definition_count, in)
		 || !array_map_init(T.disjunction_intro_nodes[i].set_definitions, ((size_t) 1) << (core::log2(set_definition_count == 0 ? 1 : set_definition_count) + 1)))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.disjunction_intro_nodes);
			return false;
		}
		for (unsigned int j = 0; j < set_definition_count; j++) {
			unsigned int key;
			unsigned int index;
			if (!read(key, in)
			 || !read(index, in))
			{
				for (auto entry : T.atoms) {
					free(entry.key);
					free(entry.value.key);
					free(entry.value.value);
				} for (auto entry : T.relations) {
					free(entry.value.key);
					free(entry.value.value);
				} for (Proof* observation : T.observations)
					free(*observation);
				for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
					if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
				for (auto entry : T.constant_types) {
					free(entry.key);
					free(*entry.value);
				} for (auto entry : T.constant_negated_types) {
					free(entry.key);
					free(*entry.value);
				}
				for (auto& node : T.disjunction_intro_nodes) free(node);
				for (Proof* axiom : T.built_in_axioms) free(*axiom);
				free(T.built_in_axioms);
				free(*T.empty_set_axiom); free(*T.NAME_ATOM);
				free(T.sets); free(T.atoms); free(T.relations);
				free(T.observations); free(T.ground_concepts);
				free(T.constant_types); free(T.constant_negated_types);
				free(T.disjunction_intro_nodes);
				return false;
			}
			T.disjunction_intro_nodes[i].set_definitions.keys[j] = key;
			T.disjunction_intro_nodes[i].set_definitions.values[j] = proofs[index];
		}
		T.disjunction_intro_nodes[i].set_definitions.size = set_definition_count;
		T.disjunction_intro_nodes[i].formula = formulas[formula_index];
		T.disjunction_intro_nodes[i].proof = proofs[proof_index];
		T.disjunction_intro_nodes[i].formula->reference_count++;
		T.disjunction_intro_nodes.length++;
	}

	if (!read(proof_count, in)
	 || !array_init(T.negated_conjunction_nodes, ((size_t) 1) << (core::log2(proof_count == 0 ? 1 : proof_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.disjunction_intro_nodes);
		return false;
	}
	for (unsigned int i = 0; i < proof_count; i++) {
		unsigned int formula_index;
		unsigned int proof_index;
		decltype(T.negated_conjunction_nodes[i].set_definitions.size) set_definition_count;
		if (!read(formula_index, in)
		 || !read(proof_index, in)
		 || !read(set_definition_count, in)
		 || !array_map_init(T.negated_conjunction_nodes[i].set_definitions, ((size_t) 1) << (core::log2(set_definition_count == 0 ? 1 : set_definition_count) + 1)))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (auto& node : T.negated_conjunction_nodes) free(node);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.disjunction_intro_nodes);
			free(T.negated_conjunction_nodes);
			return false;
		}
		for (unsigned int j = 0; j < set_definition_count; j++) {
			unsigned int key;
			unsigned int index;
			if (!read(key, in)
			 || !read(index, in))
			{
				for (auto entry : T.atoms) {
					free(entry.key);
					free(entry.value.key);
					free(entry.value.value);
				} for (auto entry : T.relations) {
					free(entry.value.key);
					free(entry.value.value);
				} for (Proof* observation : T.observations)
					free(*observation);
				for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
					if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
				for (auto entry : T.constant_types) {
					free(entry.key);
					free(*entry.value);
				} for (auto entry : T.constant_negated_types) {
					free(entry.key);
					free(*entry.value);
				}
				for (auto& node : T.disjunction_intro_nodes) free(node);
				for (auto& node : T.negated_conjunction_nodes) free(node);
				for (Proof* axiom : T.built_in_axioms) free(*axiom);
				free(T.built_in_axioms);
				free(*T.empty_set_axiom); free(*T.NAME_ATOM);
				free(T.sets); free(T.atoms); free(T.relations);
				free(T.observations); free(T.ground_concepts);
				free(T.constant_types); free(T.constant_negated_types);
				free(T.disjunction_intro_nodes);
				free(T.negated_conjunction_nodes);
				return false;
			}
			T.negated_conjunction_nodes[i].set_definitions.keys[j] = key;
			T.negated_conjunction_nodes[i].set_definitions.values[j] = proofs[index];
		}
		T.negated_conjunction_nodes[i].set_definitions.size = set_definition_count;
		T.negated_conjunction_nodes[i].formula = formulas[formula_index];
		T.negated_conjunction_nodes[i].proof = proofs[proof_index];
		T.negated_conjunction_nodes[i].formula->reference_count++;
		T.negated_conjunction_nodes.length++;
	}

	if (!read(proof_count, in)
	 || !array_init(T.implication_intro_nodes, ((size_t) 1) << (core::log2(proof_count == 0 ? 1 : proof_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (auto& node : T.negated_conjunction_nodes) free(node);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.disjunction_intro_nodes);
		free(T.negated_conjunction_nodes);
		return false;
	}
	for (unsigned int i = 0; i < proof_count; i++) {
		unsigned int formula_index;
		unsigned int proof_index;
		decltype(T.implication_intro_nodes[i].set_definitions.size) set_definition_count;
		if (!read(formula_index, in)
		 || !read(proof_index, in)
		 || !read(set_definition_count, in)
		 || !array_map_init(T.implication_intro_nodes[i].set_definitions, ((size_t) 1) << (core::log2(set_definition_count == 0 ? 1 : set_definition_count) + 1)))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (auto& node : T.negated_conjunction_nodes) free(node);
			for (auto& node : T.implication_intro_nodes) free(node);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.disjunction_intro_nodes);
			free(T.negated_conjunction_nodes);
			free(T.implication_intro_nodes);
			return false;
		}
		for (unsigned int j = 0; j < set_definition_count; j++) {
			unsigned int key;
			unsigned int index;
			if (!read(key, in)
			 || !read(index, in))
			{
				for (auto entry : T.atoms) {
					free(entry.key);
					free(entry.value.key);
					free(entry.value.value);
				} for (auto entry : T.relations) {
					free(entry.value.key);
					free(entry.value.value);
				} for (Proof* observation : T.observations)
					free(*observation);
				for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
					if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
				for (auto entry : T.constant_types) {
					free(entry.key);
					free(*entry.value);
				} for (auto entry : T.constant_negated_types) {
					free(entry.key);
					free(*entry.value);
				}
				for (auto& node : T.disjunction_intro_nodes) free(node);
				for (auto& node : T.negated_conjunction_nodes) free(node);
				for (auto& node : T.implication_intro_nodes) free(node);
				for (Proof* axiom : T.built_in_axioms) free(*axiom);
				free(T.built_in_axioms);
				free(*T.empty_set_axiom); free(*T.NAME_ATOM);
				free(T.sets); free(T.atoms); free(T.relations);
				free(T.observations); free(T.ground_concepts);
				free(T.constant_types); free(T.constant_negated_types);
				free(T.disjunction_intro_nodes);
				free(T.negated_conjunction_nodes);
				free(T.implication_intro_nodes);
				return false;
			}
			T.implication_intro_nodes[i].set_definitions.keys[j] = key;
			T.implication_intro_nodes[i].set_definitions.values[j] = proofs[index];
		}
		T.implication_intro_nodes[i].set_definitions.size = set_definition_count;
		T.implication_intro_nodes[i].formula = formulas[formula_index];
		T.implication_intro_nodes[i].proof = proofs[proof_index];
		T.implication_intro_nodes[i].formula->reference_count++;
		T.implication_intro_nodes.length++;
	}

	if (!read(proof_count, in)
	 || !array_init(T.existential_intro_nodes, ((size_t) 1) << (core::log2(proof_count == 0 ? 1 : proof_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (auto& node : T.negated_conjunction_nodes) free(node);
		for (auto& node : T.implication_intro_nodes) free(node);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.disjunction_intro_nodes);
		free(T.negated_conjunction_nodes);
		free(T.implication_intro_nodes);
		return false;
	}
	for (unsigned int i = 0; i < proof_count; i++) {
		unsigned int formula_index;
		unsigned int proof_index;
		decltype(T.existential_intro_nodes[i].set_definitions.size) set_definition_count;
		if (!read(formula_index, in)
		 || !read(proof_index, in)
		 || !read(set_definition_count, in)
		 || !array_map_init(T.existential_intro_nodes[i].set_definitions, ((size_t) 1) << (core::log2(set_definition_count == 0 ? 1 : set_definition_count) + 1)))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (auto& node : T.negated_conjunction_nodes) free(node);
			for (auto& node : T.implication_intro_nodes) free(node);
			for (auto& node : T.existential_intro_nodes) free(node);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.disjunction_intro_nodes);
			free(T.negated_conjunction_nodes);
			free(T.implication_intro_nodes);
			free(T.existential_intro_nodes);
			return false;
		}
		for (unsigned int j = 0; j < set_definition_count; j++) {
			unsigned int key;
			unsigned int index;
			if (!read(key, in)
			 || !read(index, in))
			{
				for (auto entry : T.atoms) {
					free(entry.key);
					free(entry.value.key);
					free(entry.value.value);
				} for (auto entry : T.relations) {
					free(entry.value.key);
					free(entry.value.value);
				} for (Proof* observation : T.observations)
					free(*observation);
				for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
					if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
				for (auto entry : T.constant_types) {
					free(entry.key);
					free(*entry.value);
				} for (auto entry : T.constant_negated_types) {
					free(entry.key);
					free(*entry.value);
				}
				for (auto& node : T.disjunction_intro_nodes) free(node);
				for (auto& node : T.negated_conjunction_nodes) free(node);
				for (auto& node : T.implication_intro_nodes) free(node);
				for (auto& node : T.existential_intro_nodes) free(node);
				for (Proof* axiom : T.built_in_axioms) free(*axiom);
				free(T.built_in_axioms);
				free(*T.empty_set_axiom); free(*T.NAME_ATOM);
				free(T.sets); free(T.atoms); free(T.relations);
				free(T.observations); free(T.ground_concepts);
				free(T.constant_types); free(T.constant_negated_types);
				free(T.disjunction_intro_nodes);
				free(T.negated_conjunction_nodes);
				free(T.implication_intro_nodes);
				free(T.existential_intro_nodes);
				return false;
			}
			T.existential_intro_nodes[i].set_definitions.keys[j] = key;
			T.existential_intro_nodes[i].set_definitions.values[j] = proofs[index];
		}
		T.existential_intro_nodes[i].set_definitions.size = set_definition_count;
		T.existential_intro_nodes[i].formula = formulas[formula_index];
		T.existential_intro_nodes[i].proof = proofs[proof_index];
		T.existential_intro_nodes[i].formula->reference_count++;
		T.existential_intro_nodes.length++;
	}

	decltype(T.implication_axioms.length) implication_axiom_count;
	if (!read(implication_axiom_count, in)
	 || !array_init(T.implication_axioms, ((size_t) 1) << (core::log2(implication_axiom_count == 0 ? 1 : implication_axiom_count) + 1)))
	{
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (auto& node : T.negated_conjunction_nodes) free(node);
		for (auto& node : T.implication_intro_nodes) free(node);
		for (auto& node : T.existential_intro_nodes) free(node);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.disjunction_intro_nodes);
		free(T.negated_conjunction_nodes);
		free(T.implication_intro_nodes);
		free(T.existential_intro_nodes);
		return false;
	}
	for (unsigned int i = 0; i < implication_axiom_count; i++) {
		if (!read(index, in)
		 || !read(T.implication_axioms[i].value, in))
		{
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (auto& node : T.negated_conjunction_nodes) free(node);
			for (auto& node : T.implication_intro_nodes) free(node);
			for (auto& node : T.existential_intro_nodes) free(node);
			for (auto& element : T.implication_axioms)
				free(*element.key);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.disjunction_intro_nodes);
			free(T.negated_conjunction_nodes);
			free(T.implication_intro_nodes);
			free(T.existential_intro_nodes);
			free(T.implication_axioms); return false;
		}
		T.implication_axioms[i].key = proofs[index];
		T.implication_axioms[i].key->reference_count++;
		T.implication_axioms.length++;
	}

	if (!read(T.built_in_sets, in)) {
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (auto& node : T.negated_conjunction_nodes) free(node);
		for (auto& node : T.implication_intro_nodes) free(node);
		for (auto& node : T.existential_intro_nodes) free(node);
		for (auto& element : T.implication_axioms)
			free(*element.key);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.disjunction_intro_nodes);
		free(T.negated_conjunction_nodes);
		free(T.implication_intro_nodes);
		free(T.existential_intro_nodes);
		free(T.implication_axioms); return false;
	}

	if (!hash_map_init(T.reverse_definitions, 512)) {
		for (auto entry : T.atoms) {
			free(entry.key);
			free(entry.value.key);
			free(entry.value.value);
		} for (auto entry : T.relations) {
			free(entry.value.key);
			free(entry.value.value);
		} for (Proof* observation : T.observations)
			free(*observation);
		for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
			if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
		for (auto entry : T.constant_types) {
			free(entry.key);
			free(*entry.value);
		} for (auto entry : T.constant_negated_types) {
			free(entry.key);
			free(*entry.value);
		}
		for (auto& node : T.disjunction_intro_nodes) free(node);
		for (auto& node : T.negated_conjunction_nodes) free(node);
		for (auto& node : T.implication_intro_nodes) free(node);
		for (auto& node : T.existential_intro_nodes) free(node);
		for (auto& element : T.implication_axioms)
			free(*element.key);
		for (Proof* axiom : T.built_in_axioms) free(*axiom);
		free(T.built_in_axioms);
		free(*T.empty_set_axiom); free(*T.NAME_ATOM);
		free(T.sets); free(T.atoms); free(T.relations);
		free(T.observations); free(T.ground_concepts);
		free(T.constant_types); free(T.constant_negated_types);
		free(T.built_in_sets);
		free(T.disjunction_intro_nodes);
		free(T.negated_conjunction_nodes);
		free(T.implication_intro_nodes);
		free(T.existential_intro_nodes);
		free(T.implication_axioms); return false;
	}
	for (unsigned int concept = 0; concept < T.ground_concept_capacity; concept++) {
		if (T.ground_concepts[concept].types.keys == nullptr) continue;
		auto& c = T.ground_concepts[concept];
		if (!T.reverse_definitions.check_size(T.reverse_definitions.table.size + c.definitions.length)) {
			for (auto entry : T.atoms) {
				free(entry.key);
				free(entry.value.key);
				free(entry.value.value);
			} for (auto entry : T.relations) {
				free(entry.value.key);
				free(entry.value.value);
			} for (Proof* observation : T.observations)
				free(*observation);
			for (unsigned int j = 0; j < T.ground_concept_capacity; j++)
				if (T.ground_concepts[j].types.keys != nullptr) free(T.ground_concepts[j]);
			for (auto entry : T.constant_types) {
				free(entry.key);
				free(*entry.value);
			} for (auto entry : T.constant_negated_types) {
				free(entry.key);
				free(*entry.value);
			}
			for (auto& node : T.disjunction_intro_nodes) free(node);
			for (auto& node : T.negated_conjunction_nodes) free(node);
			for (auto& node : T.implication_intro_nodes) free(node);
			for (auto& node : T.existential_intro_nodes) free(node);
			for (auto& element : T.implication_axioms)
				free(*element.key);
			for (auto entry : T.reverse_definitions)
				core::free(entry.key);
			for (Proof* axiom : T.built_in_axioms) free(*axiom);
			free(T.built_in_axioms);
			free(*T.empty_set_axiom); free(*T.NAME_ATOM);
			free(T.sets); free(T.atoms); free(T.relations);
			free(T.observations); free(T.ground_concepts);
			free(T.constant_types); free(T.constant_negated_types);
			free(T.built_in_sets); free(T.reverse_definitions);
			free(T.disjunction_intro_nodes);
			free(T.negated_conjunction_nodes);
			free(T.implication_intro_nodes);
			free(T.existential_intro_nodes);
			free(T.implication_axioms); return false;
		}
		for (unsigned int i = 0; i < c.definitions.length; i++) {
			/* add this definition to the reverse definition map */
			Proof* definition = c.definitions[i];
			if (i > 0) {
				unsigned int index = T.reverse_definitions.table.index_to_insert(*definition->formula->binary.right);
				T.reverse_definitions.table.keys[index] = *definition->formula->binary.right;
				T.reverse_definitions.values[index] = T.new_constant_offset + concept;
				T.reverse_definitions.table.size++;
			}
		}
	}

	/* we need to increment the reference counters for the proofs of set equivalence of sets defined as the same constant */
	for (unsigned int k = 0; k < T.ground_concept_capacity; k++) {
		if (T.ground_concepts[k].types.keys == nullptr) continue;

		auto& c = T.ground_concepts[k];
		for (unsigned int i = 0; i < c.definitions.length; i++) {
			Proof* definition = c.definitions[i];
			if (definition->formula->binary.right->type == FormulaType::LAMBDA) {
				uint_fast8_t arity = 1;
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

					Proof* axiom = T.sets.template get_existing_subset_axiom<false>(other_set_formula, set_formula, arity);
					axiom->reference_count++;

					axiom = T.sets.template get_existing_subset_axiom<false>(set_formula, other_set_formula, arity);
					axiom->reference_count++;
				}
			}
		}
	}
	return true;
}

template<typename ProofCalculus, typename Canonicalizer, typename Stream>
bool write(const theory<ProofCalculus, Canonicalizer>& T, Stream& out,
			const hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
			const hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	if (!write(T.new_constant_offset, out)
	 || !write(T.ground_concept_capacity, out)
	 || !write(T.ground_axiom_count, out)
	 || !write(T.built_in_axioms.length, out)
	 || !write(proof_map.get(T.empty_set_axiom), out)
	 || !write(formula_map.get(T.NAME_ATOM), out)
	 || !write(T.sets, out, proof_map, formula_map)
	 || !write(T.atoms.table.size, out))
	{
		return false;
	}

	for (const Proof* axiom : T.built_in_axioms)
		if (!write(proof_map.get(axiom), out)) return false;

	for (const auto& entry : T.atoms) {
		if (!write(entry.key, out, formula_map)
		 || !write(entry.value.key, out)
		 || !write(entry.value.value, out))
			return false;
	}
	if (!write(T.relations.table.size, out))
		return false;
	for (const auto& entry : T.relations) {
		if (!write(entry.key, out)
		 || !write(entry.value.key, out)
		 || !write(entry.value.value, out))
			return false;
	}

	if (!write(T.observations.size, out))
		return false;
	for (Proof* observation : T.observations) {
		if (!write(proof_map.get(observation), out))
			return false;
	}
	unsigned int ground_concept_count = 0;
	for (unsigned int i = 0; i < T.ground_concept_capacity; i++)
		if (T.ground_concepts[i].types.keys != nullptr) ground_concept_count++;
	if (!write(ground_concept_count, out))
		return false;
	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys == nullptr) continue;
		if (!write(i, out)
		 || !write(T.ground_concepts[i], out, proof_map, formula_map))
			return false;
	}
	if (!write(T.constant_types.size, out))
		return false;
	for (unsigned int i = 0; i < T.constant_types.size; i++) {
		if (!write(T.constant_types.keys[i], out, formula_map)
		 || !write(proof_map.get(T.constant_types.values[i]), out))
			return false;
	}
	if (!write(T.constant_negated_types.size, out))
		return false;
	for (unsigned int i = 0; i < T.constant_negated_types.size; i++) {
		if (!write(T.constant_negated_types.keys[i], out, formula_map)
		 || !write(proof_map.get(T.constant_negated_types.values[i]), out))
			return false;
	}
	if (!write(T.disjunction_intro_nodes.length, out))
		return false;
	for (const auto& entry : T.disjunction_intro_nodes) {
		if (!write(formula_map.get(entry.formula), out)
		 || !write(proof_map.get(entry.proof), out)
		 || !write(entry.set_definitions.size, out))
			return false;
		for (const auto& set_definition : entry.set_definitions) {
			if (!write(set_definition.key, out)
			 || !write(proof_map.get(set_definition.value), out))
				return false;
		}
	}
	if (!write(T.negated_conjunction_nodes.length, out))
		return false;
	for (const auto& entry : T.negated_conjunction_nodes) {
		if (!write(formula_map.get(entry.formula), out)
		 || !write(proof_map.get(entry.proof), out)
		 || !write(entry.set_definitions.size, out))
			return false;
		for (const auto& set_definition : entry.set_definitions) {
			if (!write(set_definition.key, out)
			 || !write(proof_map.get(set_definition.value), out))
				return false;
		}
	}
	if (!write(T.implication_intro_nodes.length, out))
		return false;
	for (const auto& entry : T.implication_intro_nodes) {
		if (!write(formula_map.get(entry.formula), out)
		 || !write(proof_map.get(entry.proof), out)
		 || !write(entry.set_definitions.size, out))
			return false;
		for (const auto& set_definition : entry.set_definitions) {
			if (!write(set_definition.key, out)
			 || !write(proof_map.get(set_definition.value), out))
				return false;
		}
	}
	if (!write(T.existential_intro_nodes.length, out))
		return false;
	for (const auto& entry : T.existential_intro_nodes) {
		if (!write(formula_map.get(entry.formula), out)
		 || !write(proof_map.get(entry.proof), out)
		 || !write(entry.set_definitions.size, out))
			return false;
		for (const auto& set_definition : entry.set_definitions) {
			if (!write(set_definition.key, out)
			 || !write(proof_map.get(set_definition.value), out))
				return false;
		}
	}
	if (!write(T.implication_axioms.length, out))
		return false;
	for (const auto& entry : T.implication_axioms) {
		if (!write(proof_map.get(entry.key), out)
		 || !write(entry.value, out))
			return false;
	}
	return write(T.built_in_sets, out);
}

template<typename ProofCalculus, typename Canonicalizer, typename Stream, typename PriorState>
inline bool read(theory<ProofCalculus, Canonicalizer>& T, Stream& in, PriorState& prior_state)
{
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::Language Formula;

	size_t formula_count, proof_count;
	if (!read(formula_count, in)
	 || !read(proof_count, in))
		return false;

	Formula** formulas = (Formula**) malloc(max(1, sizeof(Formula*) * formula_count));
	if (formulas == nullptr) {
		fprintf(stderr, "read ERROR: Insufficient memory for `formulas` array.\n");
		return false;
	}
	for (size_t i = 0; i < formula_count; i++) {
		formulas[i] = Formula::new_constant(0);
		if (formulas[i] == nullptr) {
			fprintf(stderr, "read ERROR: Insufficient memory for `formulas` arrays.\n");
			for (size_t j = 0; j < i; j++) free(formulas[j]);
			free(formulas); return false;
		}
	}

	Proof** proofs;
	if (!init_proof_array(proofs, proof_count)) {
		for (size_t j = 0; j < formula_count; j++) free(formulas[j]);
		free(formulas); return false;
	}

	/* read the formulas and proofs arrays */
	for (size_t i = 0; i < formula_count; i++) {
		if (!read(*formulas[i], in, formulas)) {
			for (size_t j = 0; j < formula_count; j++) {
				free(*formulas[j]); if (formulas[j]->reference_count == 0) free(formulas[j]);
			} for (size_t j = 0; j < proof_count; j++)
				free(*proofs[j]);
			free(formulas); free(proofs);
			return false;
		}
	} for (size_t i = 0; i < proof_count; i++) {
		if (!read(*proofs[i], in, proofs, formulas)) {
			for (size_t j = 0; j < formula_count; j++) {
				free(*formulas[j]); if (formulas[j]->reference_count == 0) free(formulas[j]);
			} for (size_t j = 0; j < proof_count; j++) {
				free(*proofs[j]); if (proofs[j]->reference_count == 0) free(proofs[j]);
			}
			free(formulas); free(proofs);
			return false;
		}
	}

	/* read the actual theory */
	if (!read(T, in, proofs, formulas)) {
		for (size_t j = 0; j < formula_count; j++) {
			free(*formulas[j]); if (formulas[j]->reference_count == 0) free(formulas[j]);
		} for (size_t j = 0; j < proof_count; j++) {
			free(*proofs[j]); if (proofs[j]->reference_count == 0) free(proofs[j]);
		}
		free(formulas); free(proofs);
		return false;
	} else if (!PriorState::read(prior_state, in, formulas)) {
		free(T);
		for (size_t j = 0; j < formula_count; j++) {
			free(*formulas[j]); if (formulas[j]->reference_count == 0) free(formulas[j]);
		} for (size_t j = 0; j < proof_count; j++) {
			free(*proofs[j]); if (proofs[j]->reference_count == 0) free(proofs[j]);
		}
		free(formulas); free(proofs);
		return false;
	}

	/* cleanup */
	for (size_t j = 0; j < formula_count; j++) {
#if !defined(NDEBUG)
		if (formulas[j]->reference_count == 1) {
			print("read WARNING: Found an orphan deserialized formula `", stderr);
			print(*formulas[j], stderr); print("`.\n", stderr);
			free(*formulas[j]); free(formulas[j]);
			continue;
		}
#endif
		free(*formulas[j]);
	} for (size_t j = 0; j < proof_count; j++) {
#if !defined(NDEBUG)
		if (proofs[j]->reference_count == 1) {
			print("read WARNING: Found an orphan deserialized proof:\n", stderr);
			print<built_in_predicates, typename ProofCalculus::ProofCanonicalizer, ProofCalculus::Intuitionistic>(*proofs[j], stderr);
			print('\n', stderr);
			free(*proofs[j]); free(proofs[j]);
			continue;
		}
#endif
		free(*proofs[j]);
	}
	free(formulas); free(proofs);
	return true;
}

template<typename MapType>
inline typename MapType::key_type* invert_values(const MapType& map) {
	typename MapType::key_type* inverse =
			(typename MapType::key_type*) calloc(size(map) + 1, sizeof(typename MapType::key_type));
	if (inverse == NULL) {
		fprintf(stderr, "invert ERROR: Unable to invert map. Out of memory.\n");
		return NULL;
	}
	for (const auto& entry : map)
		inverse[entry.value] = entry.key;
	return inverse;
}

template<typename ProofCalculus, typename Canonicalizer, typename Stream, typename PriorState>
inline bool write(const theory<ProofCalculus, Canonicalizer>& T, Stream& out, const PriorState& prior_state)
{
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::Language Formula;

	hash_map<const Proof*, unsigned int> proof_map(1024);
	hash_map<const Formula*, unsigned int> formula_map(2048);
	if (!get_proof_map(T, proof_map, formula_map)
	 || !PriorState::get_formula_map(prior_state, formula_map)
	 || !write((size_t) formula_map.table.size, out)
	 || !write((size_t) proof_map.table.size, out))
		return false;

	/* write the formulas array */
	const Formula** formulas = invert_values(formula_map);
	for (size_t i = 0; i < formula_map.table.size; i++) {
		if (!write(*formulas[i], out, formula_map)) {
			free(formulas);
			return false;
		}
	}
	free(formulas);

	/* write the proofs array */
	const Proof** proofs = invert_values(proof_map);
	for (size_t i = 0; i < proof_map.table.size; i++) {
		if (!write(*proofs[i], out, proof_map, formula_map)) {
			free(proofs);
			return false;
		}
	}
	free(proofs);

	/* write the actual theory */
	return write(T, out, proof_map, formula_map)
		&& PriorState::write(prior_state, out, formula_map);
}

template<typename Formula>
constexpr bool filter_operands(const Formula* formula, const array<unsigned int>& constants) { return true; }

template<typename ProofCalculus, typename Canonicalizer>
inline bool is_name_event(
		unsigned int constant, unsigned int variable_constant,
		const theory<ProofCalculus, Canonicalizer>& T,
		const array<pair<unsigned int, bool>>& new_name_events)
{
	for (const pair<unsigned int, bool>& entry : new_name_events) {
		if (entry.key == 0 && variable_constant == constant)
			return true;
		if (entry.key != constant) continue;
		if (!entry.value) return true;
		if (variable_constant == (unsigned int) built_in_predicates::NAME)
			return true;
	}
	return T.is_name_event(constant);
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool is_name_event(unsigned int constant,
		const theory<ProofCalculus, Canonicalizer>& T,
		const array<pair<unsigned int, bool>>& new_name_events)
{
	for (const pair<unsigned int, bool>& entry : new_name_events) {
		if (entry.key != constant) continue;
		if (!entry.value) return true;
	}
	return T.is_name_event(constant);
}

template<bool Hypothetical, typename ProofCalculus, typename Canonicalizer>
bool filter_constants_helper(
		const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions,
		array<pair<unsigned int, bool>>& new_name_events,
		array<const typename ProofCalculus::Language*>& arg1_of,
		array<const typename ProofCalculus::Language*>& arg2_of)
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
			if (Hypothetical)
				return constants.length != 0;
			if (right->type == TermType::VARIABLE && right->variable == variable)
				swap(left, right);
			if (right->type == TermType::CONSTANT) {
				unsigned int constant_id = right->constant;
				unsigned int index = index_of_constant(constants, constant_id);
				if (index == constants.length && (!contains_any(constants) || T.ground_concepts[constant_id - T.new_constant_offset].types.keys != nullptr))
					return false;
				constants[0] = instance_constant(constant_id);
				constants.length = 1;
			} else if (right->type == TermType::NUMBER) {
				hol_number number = right->number;
				unsigned int index = index_of_number(constants, number);
				if (index == constants.length && !contains_any(constants))
					return false;
				constants[0] = instance_number(number);
				constants.length = 1;
			} else if (right->type == TermType::STRING) {
				string* str = &right->str;
				unsigned int index = index_of_string(constants, str);
				if (index == constants.length && !contains_any(constants))
					return false;
				constants[0] = instance_string(str);
				constants.length = 1;
			} else if (right->type == TermType::VARIABLE) {
				/* no-op */
			} else {
				if (right->type == TermType::UNARY_APPLICATION && right->binary.left->type == TermType::CONSTANT
				 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
				{
					if (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1)
						arg1_of.add(right->binary.right);
					else if (right->binary.left->constant == (unsigned int) built_in_predicates::ARG2)
						arg2_of.add(right->binary.right);

					if (right->binary.right->type == TermType::CONSTANT) {
						/* we disallow objects to be arguments of themselves */
						unsigned int index = index_of_constant(constants, right->binary.right->constant);
						if (index < constants.length)
							constants.remove(index);

						/* check if the other object has this function value */
						if (right->binary.right->constant >= T.new_constant_offset
						 && T.ground_concepts[right->binary.right->constant - T.new_constant_offset].function_values.contains((unsigned int) right->binary.left->constant))
							return false;
					}
				}

				unsigned int arity = 0;
				Formula* set_definition = right;
				while (set_definition->type == FormulaType::LAMBDA) {
					set_definition = set_definition->quantifier.operand;
					arity++;
				}
				for (unsigned int i = 0; arity != 0 && i < constants.length; i++) {
					/* check that this constant could be a set */
					if (constants[i].type == instance_type::ANY) {
						continue;
					} else if (constants[i].type != instance_type::CONSTANT
							|| T.is_provably_not_a_set(constants[i].constant))
					{
						constants.remove(i--);
					} else if (constants[i].constant >= T.new_constant_offset) {
						concept<ProofCalculus>& c = T.ground_concepts[constants[i].constant - T.new_constant_offset];
						if (c.types.keys == nullptr)
							continue;
						Formula* other_set_formula = c.definitions[0]->formula->binary.right;
						unsigned int prev_arity = 0;
						while (other_set_formula->type == FormulaType::LAMBDA) {
							other_set_formula = other_set_formula->quantifier.operand;
							prev_arity++;
						}

						/* check that the arity matches, or can be changed */
						if (arity != prev_arity && (c.definitions[0]->reference_count != 2 || T.sets.set_ids.table.contains(*other_set_formula))) {
							constants.remove(i--);
							continue;
						}

						/* make sure `constants[i].constant` is not the arg1 of a name event */
						for (Proof* proof : c.definitions) {
							if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
							 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
							 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
							 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
							 && is_name_event(proof->formula->binary.right->binary.right->constant, constants[i].constant, T, new_name_events))
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
					if (is_name_event(right->binary.right->constant, T, new_name_events)) {
						for (unsigned int i = 0; i < constants.length; i++) {
							if (constants[i].type == instance_type::ANY || constants[i].type != instance_type::CONSTANT)
								continue;
							if (is_name_event(constants[i].constant, T, new_name_events))
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
					if (is_name_event(right->binary.right->constant, T, new_name_events))
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
				if (new_right->type == TermType::LAMBDA && new_right->quantifier.operand->type == TermType::UNARY_APPLICATION
				 && new_right->quantifier.operand->binary.left->type == TermType::CONSTANT
				 && new_right->quantifier.operand->binary.left->constant >= T.new_constant_offset
				 && new_right->quantifier.operand->binary.right->type == TermType::VARIABLE
				 && new_right->quantifier.operand->binary.right->variable == new_right->quantifier.variable)
				{
					if (T.ground_concepts[new_right->quantifier.operand->binary.left->constant - T.new_constant_offset].types.keys != nullptr) {
						/* we found a different concept with this definition */
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						unsigned int index = index_of_constant(constants, new_right->quantifier.operand->binary.left->constant);
						if (index == constants.length)
							return false;
						move(constants[index], constants[0]);
						constants.length = 1;
						return true;
					}
				} else {
					bool contains;
					unsigned int constant_id = T.reverse_definitions.get(*new_right, contains);
					if (contains) {
						/* we found a different concept with this definition */
						free(*new_right); if (new_right->reference_count == 0) free(new_right);
						unsigned int index = index_of_constant(constants, constant_id);
						if (index == constants.length)
							return false;
						move(constants[index], constants[0]);
						constants.length = 1;
						return true;
					}
				}
				free(*new_right); if (new_right->reference_count == 0) free(new_right);

				/* discourage sampling constants that are already defined as other sets */
				for (unsigned int i = 0; arity != 0 && i < constants.length; i++) {
					if (constants[i].type == instance_type::ANY) {
						constants[i].matching_types++;
						continue;
					} else if (constants[i].constant < T.new_constant_offset)
						continue;
					concept<ProofCalculus>& c = T.ground_concepts[constants[i].constant - T.new_constant_offset];
					for (unsigned int j = 1; j < c.definitions.length; j++) {
						if (c.definitions[j]->formula->binary.right->type != FormulaType::LAMBDA) continue;
						constants[i].mismatching_types++;
					}
				}

				/* check if this is a definition of a set of named objects */
				const Term* arg1; const Term* arg2;
				if (arity == 1 && T.is_name_scope(set_definition, arg1, arg2) && arg1->type == TermType::VARIABLE && arg1->variable == right->quantifier.variable) {
					for (unsigned int i = 0; i < constants.length; i++) {
						/* decrease the probability of sampling a set of objects with different names */
						if (constants[i].type != instance_type::CONSTANT || constants[i].constant < T.new_constant_offset
						 || T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr)
							continue;
						for (unsigned int j = 0; j < T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions.length; j++) {
							Proof* definition = T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions[j];
							Formula* other_right = definition->formula->binary.right;
							if (other_right->type != FormulaType::LAMBDA) continue;
							const Term* other_arg2;
							if (T.is_name_scope(other_right->quantifier.operand, arg1, other_arg2) && arg1->type == TermType::VARIABLE && arg1->variable == other_right->quantifier.variable) {
								if (*arg2 == *other_arg2)
									constants[i].matching_types += 2;
								else constants[i].mismatching_types += 2;
							}
						}
					}
				}
			}
		} else if (left->type == TermType::CONSTANT || right->type == TermType::CONSTANT
				|| left->type == TermType::VARIABLE || right->type == TermType::VARIABLE)
		{
			if (Hypothetical)
				return constants.length != 0;
			if ((left->type != TermType::CONSTANT || left->constant == (unsigned int) built_in_predicates::UNKNOWN) && left->type != TermType::VARIABLE)
				swap(left, right);
			if (left->type == TermType::VARIABLE)
				return true;

			/* check that this constant could be a set */
			if (right->type == FormulaType::LAMBDA) {
				if (T.is_provably_not_a_set(left->constant))
					return false;

				/* make sure `left->constant` is not the arg1 of a name event */
				if (left->constant >= T.new_constant_offset && T.ground_concepts[left->constant - T.new_constant_offset].types.keys != nullptr) {
					for (Proof* proof : T.ground_concepts[left->constant - T.new_constant_offset].definitions) {
						if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
						 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
						 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
						 && is_name_event(proof->formula->binary.right->binary.right->constant, T, new_name_events))
						{
							return false;
						}
					}
				}
			}

			bool left_is_name_event = is_name_event(left->constant, T, new_name_events);

			for (unsigned int i = 0; i < constants.length; i++) {
				if ((constants[i].type == instance_type::NUMBER || constants[i].type == instance_type::STRING)
				 && right->type == TermType::UNARY_APPLICATION && right->binary.right->type == TermType::VARIABLE
				 && right->binary.right->variable == variable && right->binary.left->type == TermType::CONSTANT
				 && (right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG2
				  || right->binary.left->constant == (unsigned int) built_in_predicates::ARG3))
				{
					constants.remove(i--);
					continue;
				} else if (constants[i].type == instance_type::ANY || constants[i].type != instance_type::CONSTANT) {
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
						if (left_is_name_event && is_name_event(new_right->binary.right->constant, constants[i].constant, T, new_name_events)) {
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
						if (is_name_event(new_right->binary.right->constant, constants[i].constant, T, new_name_events)) {
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
					if (new_new_right->type == TermType::LAMBDA && new_new_right->quantifier.operand->type == TermType::UNARY_APPLICATION
					 && new_new_right->quantifier.operand->binary.left->type == TermType::CONSTANT
					 && new_new_right->quantifier.operand->binary.left->constant >= T.new_constant_offset
					 && new_new_right->quantifier.operand->binary.right->type == TermType::VARIABLE
					 && new_new_right->quantifier.operand->binary.right->variable == new_new_right->quantifier.variable)
					{
						if (T.ground_concepts[new_new_right->quantifier.operand->binary.left->constant - T.new_constant_offset].types.keys != nullptr
						 && new_new_right->quantifier.operand->binary.left->constant != left->constant)
						{
							/* we found a different concept with this definition */
							found_conflicting_definition = true;
						}
					} else {
						bool contains;
						unsigned int constant_id = T.reverse_definitions.get(*new_new_right, contains);
						if (contains && constant_id != left->constant) {
							/* we found a different concept with this definition */
							found_conflicting_definition = true;
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
			if (Hypothetical)
				return constants.length != 0;
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
						if (is_name_event(event, constants[i].constant, T, new_name_events)) {
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
							if (is_name_event(event, constants[i].constant, T, new_name_events)) {
								free(*new_left); if (new_left->reference_count == 0) free(new_left);
								constants.remove(i--); continue;
							}
						}
					}
				}

				/* check if anything else has this as a definition */
				bool found_conflicting_definition = false;
				if (new_left->type == TermType::LAMBDA && new_left->quantifier.operand->type == TermType::UNARY_APPLICATION
				 && new_left->quantifier.operand->binary.left->type == TermType::CONSTANT
				 && new_left->quantifier.operand->binary.left->constant >= T.new_constant_offset
				 && new_left->quantifier.operand->binary.right->type == TermType::VARIABLE
				 && new_left->quantifier.operand->binary.right->variable == new_left->quantifier.variable)
				{
					if (T.ground_concepts[new_left->quantifier.operand->binary.left->constant - T.new_constant_offset].types.keys != nullptr)
						found_conflicting_definition = true;
				} else {
					if (T.reverse_definitions.table.contains(*new_left))
						found_conflicting_definition = true;
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
	} else if (formula->type == FormulaType::UNARY_APPLICATION || (formula->type == FormulaType::NOT && formula->unary.operand->type == FormulaType::UNARY_APPLICATION)) {
		bool negated = false;
		if (formula->type == FormulaType::NOT) {
			formula = formula->unary.operand;
			negated = true;
		}

		Term* left = formula->binary.left;
		Term* right = formula->binary.right;
		if (left->type == TermType::VARIABLE && left->variable == variable) {
			if (right->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, right->constant);
				if (index < constants.length)
					constants.remove(index);

				/* count the number of matching and mismatching types */
				if (!Hypothetical && right->constant >= T.new_constant_offset) {
					const concept<ProofCalculus>& c = T.ground_concepts[right->constant - T.new_constant_offset];
					for (unsigned int i = 0; i < constants.length; i++) {
						if (constants[i].type != instance_type::CONSTANT) continue;
						for (const auto& entry : c.types) {
							if (!negated && entry.key.binary.left->type == TermType::CONSTANT && entry.key.binary.left->constant == constants[i].constant) {
								constants[i].matching_types++;
							} else if (constants[i].constant >= T.new_constant_offset) {
								constants[i].mismatching_types++;
							}
						} for (const auto& entry : c.negated_types) {
							if (negated && entry.key.binary.left->type == TermType::CONSTANT && entry.key.binary.left->constant == constants[i].constant) {
								constants[i].matching_types++;
							} else if (constants[i].constant >= T.new_constant_offset) {
								constants[i].mismatching_types++;
							}
						}
					}
				}

				/* if this is a subset statement, increment `matching_types` if the set does indeed provably contain the element */
				for (unsigned int i = 0; !Hypothetical && i < constants.length; i++) {
					if (constants[i].type != instance_type::CONSTANT || constants[i].constant < T.new_constant_offset)
						continue;
					const concept<ProofCalculus>& c = T.ground_concepts[constants[i].constant - T.new_constant_offset];
					if (c.types.keys == nullptr)
						continue;
					Formula* set_formula = c.definitions[0]->formula->binary.right->quantifier.operand;

					bool contains;
					unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
					if (!contains || T.sets.sets[set_id].arity != 1) continue;

					tuple_element& element = *((tuple_element*) alloca(sizeof(tuple_element)));
					element.type = tuple_element_type::CONSTANT;
					element.constant = right->constant;
					const tuple tup = {&element, 1};
					if (T.sets.sets[set_id].provable_elements.contains(tup)) {
						constants[i].matching_types++;
					} else if (T.sets.sets[set_id].provable_elements.length == T.sets.sets[set_id].set_size) {
						/* this set is full so this `right` cannot be an element of this set */
						constants.remove(i--);
					}
				}

				/* the `true` flag here indicates that `right->constant` is a name
				   event only if the variable is instantiated as the symbol `NAME` */
				if (!negated && !new_name_events.add(make_pair(right->constant, true)))
					return false;
			} else {
				for (unsigned int i = 0; i < constants.length; i++) {
					if (constants[i].type == instance_type::CONSTANT && constants[i].constant < T.new_constant_offset) {
						if (right->type == TermType::NUMBER) {
							if (constants[i].constant != (unsigned int) built_in_predicates::NUMBER)
								constants.remove(i--);
						} else {
							constants.remove(i--);
						}
					}
				}
			}

			/* make sure x could be a set in `x(y)` */
			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) {
					continue;
				} else if (constants[i].type != instance_type::CONSTANT || T.is_provably_not_a_set(constants[i].constant)) {
					constants.remove(i);
					i--;
				} else if (Hypothetical) {
					continue;
				}
				/* make sure `constants[i].constant` is not the arg1 of a name event */
				else if (!negated && constants[i].constant >= T.new_constant_offset
					  && T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys != nullptr)
				{
					for (Proof* proof : T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions) {
						if (proof->formula->binary.right->type == FormulaType::UNARY_APPLICATION
						 && proof->formula->binary.right->binary.left->type == TermType::CONSTANT
						 && proof->formula->binary.right->binary.left->constant == (unsigned int) built_in_predicates::ARG1
						 && proof->formula->binary.right->binary.right->type == TermType::CONSTANT
						 && is_name_event(proof->formula->binary.right->binary.right->constant, constants[i].constant, T, new_name_events))
						{
							constants.remove(i--);
							break;
						}
					}
				}
			}

			/* check to make sure the args are type-correct */
			if (!Hypothetical && !negated && right->type == TermType::CONSTANT
			 && right->constant >= T.new_constant_offset
			 && T.ground_concepts[right->constant - T.new_constant_offset].types.keys != nullptr)
			{
				/* if arg1 is a `name` constant, `left` cannot be `name` */
				Term* arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(right->constant);
				if (arg1 != nullptr && arg1->type == TermType::CONSTANT && is_name_event(arg1->constant, T, new_name_events)) {
					unsigned int index = index_of_constant(constants, (unsigned int) built_in_predicates::NAME);
					if (index < constants.length) constants.remove(index);
				}

				/* if arg2 is a string, the event must be `name` */
				Term* arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(right->constant);
				if (arg2 != nullptr) {
					if (arg2->type == TermType::STRING) {
						/* this could be `name` or any set since, the set could be a collection of name events */
						for (unsigned int i = 0; i < constants.length; i++) {
							if (constants[i].type == instance_type::ANY) {
								continue;
							} else if (constants[i].type != instance_type::CONSTANT || (constants[i].constant < T.new_constant_offset && constants[i].constant != (unsigned int) built_in_predicates::NAME)) {
								constants.remove(i--);
								continue;
							}
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
					 && right->binary.right->type == TermType::CONSTANT && is_name_event(right->binary.right->constant, T, new_name_events))
					{
						unsigned int index = index_of_constant(constants, (unsigned int) built_in_predicates::NAME);
						if (index < constants.length) constants.remove(index);
						break;
					}
				}
			}
		} else if (!negated && left->type == TermType::CONSTANT && left->constant == (unsigned int) built_in_predicates::NAME) {
			/* the `false` flag here indicates that `left->constant` is
			   a name event regardless of the value of `variable` */
			if (right->type == TermType::CONSTANT) {
				if (!new_name_events.add(make_pair(right->constant, false)))
					return false;
			} else if (right->type == TermType::VARIABLE && right->variable == variable) {
				if (!new_name_events.add(make_pair(0u, false)))
					return false;
			}

			/* if `right` is in `arg1_of`, then we must exclude all constants that are name events or provably sets */
			bool is_arg1_of = false;
			for (unsigned int i = 0; !is_arg1_of && i < arg1_of.length; i++)
				if (*arg1_of[i] == *right) is_arg1_of = true;
			if (is_arg1_of) {
				for (unsigned int i = 0; i < constants.length; i++) {
					if (constants[i].type == instance_type::ANY) continue;
					else if (constants[i].type != instance_type::CONSTANT || is_name_event(constants[i].constant, T, new_name_events) || T.is_provably_a_set(constants[i].constant)) {
						constants.remove(i--);
						continue;
					}
				}
			}
		} else if (!negated && left->type == TermType::CONSTANT && right->type == TermType::VARIABLE && right->variable == variable) {
			unsigned int index = set_definitions.last_index_of(left->constant);
			if (index != static_cast<unsigned int>(-1)) {
				Proof* set_definition = set_definitions.values[index];
				Formula* lambda_formula = set_definition->formula->binary.right;

				uint_fast8_t arity = 0;
				Formula* set_formula = lambda_formula;
				while (set_formula->type == FormulaType::LAMBDA) {
					set_formula = set_formula->quantifier.operand;
					arity++;
				}
				if (arity != 1) return false;

				if (!filter_constants_helper<Hypothetical>(T, set_formula, lambda_formula->quantifier.variable, constants, set_definitions, new_name_events, arg1_of, arg2_of))
					return false;
			}
		}
		if (right->type == TermType::VARIABLE && right->variable == variable) {
			if (left->type == TermType::CONSTANT) {
				/* disallow statements of the form `x(x)` */
				unsigned int index = index_of_constant(constants, left->constant);
				if (index < constants.length)
					constants.remove(index);
			}
			/* make sure y could not be a set in `x(y)` */
			for (unsigned int i = 0; !Hypothetical && i < constants.length; i++) {
				if (constants[i].type != instance_type::ANY && constants[i].type != instance_type::CONSTANT) {
					if (left->type == TermType::CONSTANT && T.new_constant_offset > left->constant && left->constant != (unsigned int) built_in_predicates::NUMBER)
						constants.remove(i--);
					continue;
				} else if (constants[i].type == instance_type::CONSTANT && T.is_provably_a_set(constants[i].constant)) {
					constants.remove(i--);
				}
			}

			if (!Hypothetical && !is_ambiguous(*left) && !(left->type == TermType::CONSTANT && left->constant >= T.new_constant_offset)) {
				/* count the number of matching and mismatching types */
				for (unsigned int i = 0; i < constants.length; i++) {
					if (constants[i].type != instance_type::CONSTANT) continue;
					const concept<ProofCalculus>& c = T.ground_concepts[constants[i].constant - T.new_constant_offset];
					if (c.types.keys == nullptr) continue;
					for (const auto& entry : c.types) {
						if (!negated && *entry.key.binary.left == *left) {
							constants[i].matching_types++;
						} else if (constants[i].constant >= T.new_constant_offset) {
							constants[i].mismatching_types++;
						}
					} for (const auto& entry : c.negated_types) {
						if (negated && *entry.key.binary.left == *left) {
							constants[i].matching_types++;
						} else if (constants[i].constant >= T.new_constant_offset) {
							constants[i].mismatching_types++;
						}
					}
				}
			}

			/* if this is a subset statement, increment `matching_types` if the set does indeed provably contain the element */
			if (left->type == TermType::CONSTANT && left->constant >= T.new_constant_offset) {
				const concept<ProofCalculus>& c = T.ground_concepts[left->constant - T.new_constant_offset];
				Formula* set_formula = c.definitions[0]->formula->binary.right->quantifier.operand;

				bool contains;
				unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
				if (contains && T.sets.sets[set_id].arity == 1) {
					for (unsigned int i = 0; i < constants.length; i++) {
						if (constants[i].type == instance_type::ANY) {
							if (T.sets.sets[set_id].provable_elements.length == T.sets.sets[set_id].set_size) {
								/* the set is full, so a new constant cannot be a member of the set */
								constants.remove(i--);
							}
							continue;
						}

						tuple_element& element = *((tuple_element*) alloca(sizeof(tuple_element)));
						if (constants[i].type == instance_type::CONSTANT) {
							element.type = tuple_element_type::CONSTANT;
							element.constant = constants[i].constant;
						} else if (constants[i].type == instance_type::NUMBER) {
							element.type = tuple_element_type::NUMBER;
							element.number = constants[i].number;
						} else if (constants[i].type == instance_type::STRING) {
							element.type = tuple_element_type::STRING;
							element.str = *constants[i].str;
						}
						const tuple tup = {&element, 1};
						if (T.sets.sets[set_id].provable_elements.contains(tup)) {
							constants[i].matching_types++;
						} else if (T.sets.sets[set_id].provable_elements.length == T.sets.sets[set_id].set_size) {
							/* the set is full, so `constants[i]` cannot be a member of the set */
							constants.remove(i--);
						}
						free(element);
					}
				}
			}

			for (unsigned int i = 0; i < constants.length; i++) {
				if (constants[i].type == instance_type::ANY) continue;
				if (constants[i].type != instance_type::CONSTANT && constants[i].type != instance_type::NUMBER) {
					constants.remove(i--);
					continue;
				}
				if (Hypothetical || negated || constants[i].type != instance_type::CONSTANT)
					continue;

				/* check to make sure the args are type-correct */
				unsigned int arg = constants[i].constant;
				if (arg < T.new_constant_offset || T.ground_concepts[arg - T.new_constant_offset].types.keys == nullptr)
					continue;
				if (left->type == TermType::CONSTANT && left->constant == (unsigned int) built_in_predicates::NAME) {
					/* arg1 must be a non-name constant and arg2 must be a string */
					Term* arg1 = T.template get_arg<(unsigned int) built_in_predicates::ARG1>(arg);
					if (arg1 != nullptr) {
						if (arg1->type != TermType::CONSTANT || is_name_event(arg1->constant, constants[i].constant, T, new_name_events)) {
							constants.remove(i--);
							continue;
						}
					}
					Term* arg2 = T.template get_arg<(unsigned int) built_in_predicates::ARG2>(arg);
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
						 && right->binary.right->type == TermType::CONSTANT && is_name_event(right->binary.right->constant, constants[i].constant, T, new_name_events))
						{
							is_arg1_of_name_event = true;
							break;
						}
					}
					if (is_arg1_of_name_event) {
						constants.remove(i--);
						continue;
					}
				} else if (left->type == TermType::CONSTANT && !set_definitions.contains(left->constant)) { /* ignore the possibility where we have `c(x)` where `c` is a constant with a set definition */
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
			if (!filter_constants_helper<Hypothetical>(T, formula->array.operands[i], variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)) return false;
	} else if (formula->type == FormulaType::OR) {
		for (unsigned int i = 0; i < formula->array.length; i++)
			if (!filter_constants_helper<true>(T, formula->array.operands[i], variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)) return false;
	} else if (formula->type == FormulaType::EXISTS) {
		size_t old_arg1_of_length = arg1_of.length;
		size_t old_arg2_of_length = arg2_of.length;
		if (!filter_constants_helper<Hypothetical>(T, formula->quantifier.operand, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)) return false;
		for (size_t i = old_arg1_of_length; i < arg1_of.length; i++) {
			if (arg1_of[i]->type == TermType::VARIABLE && arg1_of[i]->variable == formula->quantifier.variable)
				arg1_of.remove(i--);
		} for (size_t i = old_arg2_of_length; i < arg2_of.length; i++) {
			if (arg2_of[i]->type == TermType::VARIABLE && arg2_of[i]->variable == formula->quantifier.variable)
				arg2_of.remove(i--);
		}
		if (!Hypothetical) {
			if (formula->quantifier.operand->type == FormulaType::AND) {
				/* optimization: if this is a name scope, and it declares the name of
				   `variable`, then increase the probability of sampling a constant
				   with the same name */
				bool is_name_scope = false;
				const Term* arg1 = nullptr;
				const Term* arg2 = nullptr;
				for (unsigned int i = 0; i < formula->quantifier.operand->array.length; i++) {
					const Formula* conjunct = formula->quantifier.operand->array.operands[i];
					if (conjunct->type == TermType::UNARY_APPLICATION && conjunct->binary.left->type == TermType::CONSTANT
					 && conjunct->binary.left->constant == (unsigned int) built_in_predicates::NAME)
					{
						is_name_scope = true;
					} else if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::UNARY_APPLICATION
							&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
							&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
							&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
							&& conjunct->binary.left->binary.right->variable == formula->quantifier.variable)
					{
						arg1 = conjunct->binary.right;
					} else if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::UNARY_APPLICATION
							&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
							&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
							&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
							&& conjunct->binary.left->binary.right->variable == formula->quantifier.variable)
					{
						arg2 = conjunct->binary.right;
					}
				}
				if (is_name_scope && arg1 != nullptr && arg1->type == TermType::VARIABLE && arg1->variable == variable && arg2 != nullptr) {
					unsigned int any_index = constants.length;
					bool found_concept_with_name = false;
					for (unsigned int i = 0; i < constants.length; i++) {
						if (constants[i].type == instance_type::ANY) {
							any_index = i;
							continue;
						} else if (constants[i].type != instance_type::CONSTANT) {
							constants.remove(i--);
							continue;
						} else if (constants[i].constant < T.new_constant_offset || T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr) {
							continue;
						}
						array<Term*> names(4);
						if (!T.get_concept_names(constants[i].constant, names))
							return false;
						for (Term* name : names) {
							if (*name == *arg2) {
								found_concept_with_name = true;
								constants[i].matching_types += 2;
							} else {
								constants[i].mismatching_types += 2;
							}
						}
					}
					if (!found_concept_with_name && any_index != constants.length)
						constants[any_index].matching_types += 2;
				}
				/* optimization: if this quantifier defines a set with a fixed
				   size and members of that set, then make sure the set is
				   large enough to fit all elements */
				Formula* set_definition = nullptr;
				uint64_t set_size = UINT_MAX;
				array<Term*> members(4);
				for (unsigned int i = 0; i < formula->quantifier.operand->array.length; i++) {
					Formula* conjunct = formula->quantifier.operand->array.operands[i];
					if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::VARIABLE
					 && conjunct->binary.left->variable == formula->quantifier.variable
					 && conjunct->binary.right->type == TermType::LAMBDA
					 && conjunct->binary.right->quantifier.operand->type != TermType::LAMBDA)
					{
						set_definition = conjunct->binary.right;
					} else if (conjunct->type == TermType::EQUALS && conjunct->binary.left->type == TermType::UNARY_APPLICATION
							&& conjunct->binary.left->binary.left->type == TermType::CONSTANT
							&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
							&& conjunct->binary.left->binary.right->type == TermType::VARIABLE
							&& conjunct->binary.left->binary.right->variable == formula->quantifier.variable
							&& conjunct->binary.right->type == TermType::NUMBER)
					{
						if (conjunct->binary.right->number.decimal != 0
						 || conjunct->binary.right->number.integer < 0)
						{
							constants.clear();
							return false;
						}
						set_size = (uint64_t) conjunct->binary.right->number.integer;
					} else if (conjunct->type == TermType::UNARY_APPLICATION && conjunct->binary.left->type == TermType::VARIABLE
							&& conjunct->binary.left->variable == formula->quantifier.variable
							&& (conjunct->binary.right->type == TermType::CONSTANT || conjunct->binary.right->type == TermType::NUMBER || conjunct->binary.right->type == TermType::STRING))
					{
						bool contains = false;
						for (unsigned int i = 0; !contains && i < members.length; i++)
							if (*members[i] == *conjunct->binary.right) contains = true;
						if (!contains) members.add(conjunct->binary.right);
					}
				}
				if (set_definition != nullptr && set_size != UINT_MAX) {
					array<unsigned int> free_variables(4);
					get_free_variables(*set_definition, free_variables);
					if (free_variables.length == 1 && free_variables[0] == variable) {
						array<Formula*> quantifiers(set_definition->quantifier.variable);
						quantifiers[variable - 1] = hol_term::new_exists(variable, &HOL_ANY);
						HOL_ANY.reference_count++;
						quantifiers[set_definition->quantifier.variable - 1] = set_definition;
						quantifiers.length = set_definition->quantifier.variable;
						array<variable_assignment> provable_values(1);
						if (!init(provable_values[0], set_definition->quantifier.variable))
							return false;
						provable_values.length = 1;
						typename theory<ProofCalculus, Canonicalizer>::default_prover prover(T.sets, T.implication_axioms);
						prover.h.set_ids[0] = 0;
						prover.h.set_ids.length = 1;
						T.template is_provable_without_abduction<false>(set_definition->quantifier.operand, quantifiers, provable_values, prover);
						free(*quantifiers[variable - 1]); free(quantifiers[variable - 1]);

						for (unsigned int i = 0; i < constants.length; i++) {
							if (constants[i].type == instance_type::ANY) continue;

							/* count the number of provable elements of the set */
							unsigned int provable_set_size = members.length;
							for (unsigned int j = 0; provable_set_size <= set_size && j < provable_values.length; j++) {
								const instantiation& variable_value = provable_values[j].assignment.values[variable - 1];
								if (constants[i].type == instance_type::CONSTANT) {
									if (variable_value.type != instantiation_type::CONSTANT || variable_value.constant != constants[i].constant)
										continue;
								} else if (constants[i].type == instance_type::NUMBER) {
									if (variable_value.type != instantiation_type::NUMBER || variable_value.number != constants[i].number)
										continue;
								} else if (constants[i].type == instance_type::STRING) {
									if (variable_value.type != instantiation_type::STRING || variable_value.str != *constants[i].str)
										continue;
								}

								/* check if the provable element is already in `members` */
								const instantiation& provable_member = provable_values[j].assignment.values[set_definition->quantifier.variable - 1];
								bool is_member = false;
								for (const Term* member : members) {
									if (provable_member.type == instantiation_type::CONSTANT && member->type == TermType::CONSTANT && provable_member.constant == member->constant) {
										is_member = true; break;
									} else if (provable_member.type == instantiation_type::NUMBER && member->type == TermType::NUMBER && provable_member.number == member->number) {
										is_member = true; break;
									} else if (provable_member.type == instantiation_type::STRING && member->type == TermType::STRING && provable_member.str == member->str) {
										is_member = true; break;
									}
								}
								if (is_member) continue;
								provable_set_size++;
							}
							if (provable_set_size > set_size) {
								/* the set is provably too large */
								constants.remove(i--);
							}
						}
						for (auto& element : provable_values) core::free(element);
					}
				}
			}
			/* optimization: check if this quantification is relativized by a
			   provably full set; if so, make sure this quantifier can be
			   instantiated by a member of that set */
			array_view<Formula* const> set_formula_conjuncts(
					(formula->quantifier.operand->type == FormulaType::AND) ? formula->quantifier.operand->array.operands : &formula->quantifier.operand,
					(formula->quantifier.operand->type == FormulaType::AND) ? formula->quantifier.operand->array.length : 1);
			for (const Formula* conjunct : set_formula_conjuncts) {
				if (conjunct->type == TermType::UNARY_APPLICATION && conjunct->binary.right->type == TermType::VARIABLE && conjunct->binary.right->variable == formula->quantifier.variable) {
					if (conjunct->binary.left->type == TermType::CONSTANT && conjunct->binary.left->constant >= T.new_constant_offset) {
						/* check if `conjunct->binary.left->constant` is a full set */
						Formula* set_formula = T.ground_concepts[conjunct->binary.left->constant - T.new_constant_offset].definitions[0]->formula->binary.right;
						while (set_formula->type == FormulaType::LAMBDA)
							set_formula = set_formula->quantifier.operand;
						bool contains;
						unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
						if (!contains) continue;
						if (T.sets.sets[set_id].arity != 1) return false;
						if (T.sets.sets[set_id].provable_elements.length != T.sets.sets[set_id].set_size)
							continue;

						array<instance> other_constants(max(T.sets.sets[set_id].provable_elements.length, 1));
						for (const tuple& element : T.sets.sets[set_id].provable_elements) {
							switch (element.elements[0].type) {
							case tuple_element_type::CONSTANT:
								other_constants[other_constants.length].type = instance_type::CONSTANT;
								other_constants[other_constants.length].constant = element.elements[0].constant;
								other_constants[other_constants.length].matching_types = 0;
								other_constants[other_constants.length++].mismatching_types = 0;
								break;
							case tuple_element_type::NUMBER:
								other_constants[other_constants.length].type = instance_type::NUMBER;
								other_constants[other_constants.length].number = element.elements[0].number;
								other_constants[other_constants.length].matching_types = 0;
								other_constants[other_constants.length++].mismatching_types = 0;
								break;
							case tuple_element_type::STRING:
								other_constants[other_constants.length].type = instance_type::STRING;
								other_constants[other_constants.length].str = &element.elements[0].str;
								other_constants[other_constants.length].matching_types = 0;
								other_constants[other_constants.length++].mismatching_types = 0;
								break;
							}
						}

						if (!filter_constants_helper<Hypothetical>(T, formula->quantifier.operand, formula->quantifier.variable, other_constants, set_definitions, new_name_events, arg1_of, arg2_of) || other_constants.length == 0)
							return false;
						arg1_of.length = old_arg1_of_length;
						arg2_of.length = old_arg2_of_length;

					} else if (conjunct->binary.left->type == TermType::VARIABLE && conjunct->binary.left->variable == variable) {
						for (unsigned int i = 0; i < constants.length; i++) {
							if (constants[i].type == instance_type::ANY) {
								continue;
							} else if (constants[i].type != instance_type::CONSTANT) {
								constants.remove(i--);
								continue;
							} else if (constants[i].constant < T.new_constant_offset
									|| T.ground_concepts[constants[i].constant - T.new_constant_offset].types.keys == nullptr)
							{
								continue;
							}

							/* check if `constants[i].constant` is a full set */
							Formula* set_formula = T.ground_concepts[constants[i].constant - T.new_constant_offset].definitions[0]->formula->binary.right;
							while (set_formula->type == FormulaType::LAMBDA)
								set_formula = set_formula->quantifier.operand;
							bool contains;
							unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
							if (!contains) continue;
							if (T.sets.sets[set_id].arity != 1) return false;
							if (T.sets.sets[set_id].provable_elements.length != T.sets.sets[set_id].set_size)
								continue;

							array<instance> other_constants(max(T.sets.sets[set_id].provable_elements.length, 1));
							for (const tuple& element : T.sets.sets[set_id].provable_elements) {
								switch (element.elements[0].type) {
								case tuple_element_type::CONSTANT:
									other_constants[other_constants.length].type = instance_type::CONSTANT;
									other_constants[other_constants.length].constant = element.elements[0].constant;
									other_constants[other_constants.length].matching_types = 0;
									other_constants[other_constants.length++].mismatching_types = 0;
									break;
								case tuple_element_type::NUMBER:
									other_constants[other_constants.length].type = instance_type::NUMBER;
									other_constants[other_constants.length].number = element.elements[0].number;
									other_constants[other_constants.length].matching_types = 0;
									other_constants[other_constants.length++].mismatching_types = 0;
									break;
								case tuple_element_type::STRING:
									other_constants[other_constants.length].type = instance_type::STRING;
									other_constants[other_constants.length].str = &element.elements[0].str;
									other_constants[other_constants.length].matching_types = 0;
									other_constants[other_constants.length++].mismatching_types = 0;
									break;
								}
							}

							if (!filter_constants_helper<Hypothetical>(T, formula->quantifier.operand, formula->quantifier.variable, other_constants, set_definitions, new_name_events, arg1_of, arg2_of) || other_constants.length == 0)
								constants.remove(i--);
							arg1_of.length = old_arg1_of_length;
							arg2_of.length = old_arg2_of_length;
						}
					}
				}
			}
		}
	} else if (formula->type == FormulaType::FOR_ALL) {
		size_t old_arg1_of_length = arg1_of.length;
		size_t old_arg2_of_length = arg2_of.length;
		if (!filter_constants_helper<Hypothetical>(T, formula->quantifier.operand, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)) return false;
		for (size_t i = old_arg1_of_length; i < arg1_of.length; i++) {
			if (arg1_of[i]->type == TermType::VARIABLE && arg1_of[i]->variable == formula->quantifier.variable)
				arg1_of.remove(i--);
		} for (size_t i = old_arg2_of_length; i < arg2_of.length; i++) {
			if (arg2_of[i]->type == TermType::VARIABLE && arg2_of[i]->variable == formula->quantifier.variable)
				arg2_of.remove(i--);
		}
		return true;
	} else if (formula->type == FormulaType::IF_THEN) {
		return filter_constants_helper<true>(T, formula->binary.left, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)
			&& filter_constants_helper<true>(T, formula->binary.right, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of);
	}
	return (constants.length > 0);
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool get_possible_constants(const theory<ProofCalculus, Canonicalizer>& T, array<instance>& constants)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

	array<hol_number> numbers(64); array<string*> strings(64);
	for (unsigned int i = 0; i < T.ground_concept_capacity; i++) {
		if (T.ground_concepts[i].types.keys != nullptr) {
			constants[constants.length].type = instance_type::CONSTANT;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].constant = T.new_constant_offset + i;

			for (const auto& entry : T.ground_concepts[i].function_values) {
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
	} for (unsigned int i = 1; i < T.sets.set_count + 1; i++) {
		if (T.sets.sets[i].size_axioms.data == nullptr) continue;
		hol_number number;
		number.integer = T.sets.sets[i].set_size;
		number.decimal = 0;
		if (!numbers.contains(number) && !numbers.add(number))
			return false;
	}
	constants[constants.length].matching_types = 0;
	constants[constants.length].mismatching_types = 0;
	constants[constants.length++].type = instance_type::ANY;
	if (!constants.ensure_capacity(constants.length + numbers.length + strings.length + T.constant_types.size + T.constant_negated_types.size))
		return false;
	for (hol_number number : numbers) {
		constants[constants.length].type = instance_type::NUMBER;
		constants[constants.length].matching_types = 0;
		constants[constants.length].mismatching_types = 0;
		constants[constants.length++].number = number;
	} for (string* str : strings) {
		constants[constants.length].type = instance_type::STRING;
		constants[constants.length].matching_types = 0;
		constants[constants.length].mismatching_types = 0;
		constants[constants.length++].str = str;
	} for (const auto& entry : T.constant_types) {
		if (entry.key.binary.right->type == TermType::NUMBER) {
			constants[constants.length].type = instance_type::NUMBER;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].number = entry.key.binary.right->number;
		} else if (entry.key.binary.right->type == TermType::STRING) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].str = &entry.key.binary.right->str;
		}
	} for (const auto& entry : T.constant_negated_types) {
		if (entry.key.binary.right->type == TermType::NUMBER) {
			constants[constants.length].type = instance_type::NUMBER;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].number = entry.key.binary.right->number;
		} else if (entry.key.binary.right->type == TermType::STRING) {
			constants[constants.length].type = instance_type::STRING;
			constants[constants.length].matching_types = 0;
			constants[constants.length].mismatching_types = 0;
			constants[constants.length++].str = &entry.key.binary.right->str;
		}
	}
	return true;
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool get_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions)
{
	typedef typename ProofCalculus::Language Formula;

	if (!get_possible_constants(T, constants))
		return false;

	array<pair<unsigned int, bool>> new_name_events(8);
	array<const Formula*> arg1_of(4);
	array<const Formula*> arg2_of(4);
	if (!filter_constants_helper<false>(T, formula, variable, constants, set_definitions, new_name_events, arg1_of, arg2_of))
		return false;
	shuffle(constants);

	unsigned int to_swap = constants.length;
	for (unsigned int i = 0; i < constants.length; i++) {
		if (constants[i].type == instance_type::ANY) {
			if (to_swap == constants.length) to_swap = i;
		} else if (constants[i].matching_types >= 4) {
			to_swap = i;
		}
	}
	if (to_swap < constants.length)
		swap(constants[to_swap], constants[0]);
	return true;
}

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer>
bool is_impossible(
		typename ProofCalculus::Language* formula,
		const theory<ProofCalculus, Canonicalizer>& T)
{
	typedef typename ProofCalculus::Language Formula;

	array<Formula*> quantifiers(1);
	array<variable_assignment> temp_possible_values(1);
	if (!init(temp_possible_values[0], 0))
		return true;
	temp_possible_values.length = 1;
	typename theory<ProofCalculus, Canonicalizer>::default_prover prover(T.sets, T.implication_axioms);
	prover.h.set_ids[0] = 0;
	prover.h.set_ids.length = 1;
	bool result = T.template is_provable_without_abduction<!Contradiction>(formula, quantifiers, temp_possible_values, prover);
	for (auto& element : temp_possible_values) free(element);
	return result;
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
inline void finished_constants(const Formula* formula) { }

template<typename Formula>
inline void finished_operand_indices(const Formula* formula) { }

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
bool init(theory_sample<Proof>& sample, const array<Proof*>& proofs,
		const array<typename Proof::FormulaType*>& extra_axioms, double log_probability)
{
	typedef typename Proof::FormulaType Formula;

	sample.log_probability = log_probability;
	sample.proofs = (Proof**) malloc(max((size_t) 1, sizeof(Proof*) * proofs.length));
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
	hash_map<const Formula*, Formula*> formula_map(64);
	for (Proof* proof : proofs) {
		if (!Proof::clone(proof, sample.proofs[sample.proof_count], proof_map, formula_map)) {
			for (unsigned int j = 0; j < sample.proof_count; j++) { free(*sample.proofs[j]); free(sample.proofs[j]); }
			free(sample.proofs); free(sample.extra_axioms);
			return false;
		}
		sample.proof_count++;
	}

	sample.extra_axiom_count = 0;
	for (Formula* extra_axiom : extra_axioms) {
		if (!clone(extra_axiom, sample.extra_axioms[sample.extra_axiom_count], formula_map)) {
			for (unsigned int j = 0; j < sample.proof_count; j++) { free(*sample.proofs[j]); free(sample.proofs[j]); }
			for (unsigned int j = 0; j < sample.extra_axiom_count; j++) {
				free(*sample.extra_axioms[j]); if (sample.extra_axioms[j]->reference_count == 0) free(sample.extra_axioms[j]);
			}
			free(sample.proofs); free(sample.extra_axioms);
			return false;
		}
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
	constexpr inline bool accept(const array<Proof*>& sample,
			const array<typename Proof::FormulaType*>& extra_axioms,
			double proof_prior_diff) const
	{
		return true;
	}

	template<typename Proof>
	constexpr inline bool accept_with_observation_changes(const array<Proof*>& sample,
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
	Proof* test_proof;

#if !defined(NDEBUG)
	std::function<double(void)> compute_current_log_probability;
#endif

	template<typename ProofPrior>
	log_probability_collector(const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior, Proof* test_proof = nullptr) : test_proof(test_proof)
	{
		/* initialize `current_log_probability` */
		auto compute_log_probability = [&]() {
			null_collector collector;
			array<Formula*> extra_axioms(16);
			T.get_extra_axioms(extra_axioms);
			double value = log_probability(T.observations, extra_axioms, proof_prior, collector);
//fprintf(stderr, "log probability of theory: %lf\n", value);
//T.print_axioms(stderr, *debug_terminal_printer);
//T.print_disjunction_introductions(stderr, *debug_terminal_printer);
			return value;
		};
#if !defined(NDEBUG)
		compute_current_log_probability = compute_log_probability;
#endif
		current_log_probability = compute_log_probability();
	}

	constexpr inline bool has_prior(const Proof* proof) const {
		return true;
	}

	bool accept(const array<typename ProofCalculus::Proof*>& sample,
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

	inline bool accept_with_observation_changes(const array<typename ProofCalculus::Proof*>& sample,
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
};

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
inline log_probability_collector<ProofCalculus, Canonicalizer> make_log_probability_collector(
		const theory<ProofCalculus, Canonicalizer>& T, ProofPrior& proof_prior,
		typename ProofCalculus::Proof* test_proof = nullptr)
{
	return log_probability_collector<ProofCalculus, Canonicalizer>(T, proof_prior, test_proof);
}

template<typename ProofCalculus, typename Canonicalizer, typename OnProofSampleFunction = no_op>
struct model_evidence_collector
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename ProofCalculus::ProofType ProofType;

	const theory<ProofCalculus, Canonicalizer>& T;
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
	model_evidence_collector(const theory<ProofCalculus, Canonicalizer>& T_src,
			ProofPrior& proof_prior, Proof* test_proof, TheorySampleCollector& sample_collector,
			OnProofSampleFunction on_new_proof_sample = no_op()) :
		T(T_src), samples(1024), observation_count(T.observations.length), test_proof(test_proof), on_new_proof_sample(on_new_proof_sample)
	{
		/* initialize `current_log_probability` */
		array<Formula*> extra_axioms(16);
		T.get_extra_axioms(extra_axioms);
#if !defined(NDEBUG)
		compute_current_log_probability = [&]() {
			array<Formula*> extra_axioms(16);
			T.get_extra_axioms(extra_axioms);
			double value = log_probability(T.observations, extra_axioms, proof_prior, sample_collector);
/*fprintf(stderr, "log probability of theory: %lf\n", value);
T.print_axioms(stderr, *debug_terminal_printer);
T.print_disjunction_introductions(stderr, *debug_terminal_printer);*/
			return value;
		};
		current_log_probability = compute_current_log_probability();
#else
		current_log_probability = log_probability(T.observations, extra_axioms, proof_prior, sample_collector);
#endif

		/* add the first sample */
		theory_sample<Proof>& new_sample = *((theory_sample<Proof>*) alloca(sizeof(theory_sample<Proof>)));
		if (!init(new_sample, T.observations, extra_axioms, current_log_probability))
			throw std::runtime_error("Failed to initialize first theory_sample.");
		on_new_proof_sample(T, test_proof, current_log_probability);

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

	bool accept(const array<typename ProofCalculus::Proof*>& sample,
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
		on_new_proof_sample(T, test_proof, current_log_probability);
		move(new_sample, samples.keys[bucket]);
		samples.size++;
		return true;
	}

	inline bool accept_with_observation_changes(const array<typename ProofCalculus::Proof*>& sample,
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

	inline bool accept(const array<typename ProofCalculus::Proof*>& sample,
			const array<typename ProofCalculus::Language*>& extra_axioms, double proof_prior_diff)
	{
		return internal_collector.accept(sample, extra_axioms, proof_prior_diff);
	}

	bool accept_with_observation_changes(const array<typename ProofCalculus::Proof*>& sample,
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
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	set_changes<Formula> set_diff;
	Proof* new_proof = T.add_formula(logical_form, set_diff, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!proof_axioms.template add<false>(new_proof, set_diff.new_set_axioms, proof_prior)) {
		T.template remove_formula<true>(new_proof, set_diff);
		return -std::numeric_limits<double>::infinity();
	}
	set_diff.clear();

	model_evidence_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	for (unsigned int t = 0; t < num_samples; t++)
{
/*fprintf(stderr, "DEBUG: t = %u\n", t);
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
if (print_debug) T.print_disjunction_introductions(stderr, *debug_terminal_printer);*/
		do_mh_step(T, proof_prior, proof_axioms, collector);
}

	T.template remove_formula<false>(collector.test_proof, set_diff);
	proof_axioms.template subtract<false>(collector.test_proof, set_diff.old_set_axioms, proof_prior);
	free(*collector.test_proof);
	if (collector.test_proof->reference_count == 0)
		free(collector.test_proof);
	return collector.total_log_probability();
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior>
double log_joint_probability_of_truth(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form,
		unsigned int num_samples, unsigned int num_restarts, unsigned int burn_in,
		theory<ProofCalculus, Canonicalizer>& T_MAP,
		typename ProofCalculus::Proof*& proof_MAP)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	set_changes<Formula> set_diff;
	Proof* new_proof = T.add_formula(logical_form, set_diff, new_constant);
	if (new_proof == nullptr) {
		return -std::numeric_limits<double>::infinity();
	} else if (!proof_axioms.template add<false>(new_proof, set_diff.new_set_axioms, proof_prior)) {
		T.template remove_formula<true>(new_proof, set_diff);
		return -std::numeric_limits<double>::infinity();
	}
	set_diff.clear();

	array_map<const Proof*, Proof*> proof_map(64);
	hash_map<const Formula*, Formula*> formula_map(128);
	if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_MAP, proof_map, formula_map)) {
		T.template remove_formula<false>(new_proof, set_diff);
		proof_axioms.template subtract<false>(new_proof, set_diff.old_set_axioms, proof_prior);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		return -std::numeric_limits<double>::infinity();
	}
	proof_MAP = proof_map.get(new_proof);

	provability_collector<ProofCalculus, Canonicalizer> collector(T, proof_prior, new_proof);
	double max_log_probability = collector.internal_collector.current_log_probability;
	for (unsigned int i = 0; i < num_restarts; i++) {
		for (unsigned int t = 0; t < num_samples; t++)
		{
/*fprintf(stderr, "DEBUG: i = %u, t = %u\n", i, t);
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
if (print_debug) T.template print_axioms<true>(stderr, *debug_terminal_printer);
if (print_debug) T.print_disjunction_introductions(stderr, *debug_terminal_printer);*/
			do_mh_step(T, proof_prior, proof_axioms, collector);
			if (collector.internal_collector.current_log_probability > max_log_probability) {
				free(T_MAP); proof_map.clear(); formula_map.clear();
				if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_MAP, proof_map, formula_map)) {
					T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
					proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
					free(*collector.internal_collector.test_proof);
					if (collector.internal_collector.test_proof->reference_count == 0)
						free(collector.internal_collector.test_proof);
					return -std::numeric_limits<double>::infinity();
				}
				proof_MAP = proof_map.get(collector.internal_collector.test_proof);
				max_log_probability = collector.internal_collector.current_log_probability;
			}
		}

		if (i + 1 < num_restarts) {
			for (unsigned int t = 0; t < burn_in; t++)
				do_exploratory_mh_step(T, proof_prior, proof_axioms, collector);
		}
	}

	T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
	free(*collector.internal_collector.test_proof);
	if (collector.internal_collector.test_proof->reference_count == 0)
		free(collector.internal_collector.test_proof);
	return collector.total_log_probability();
}

template<typename OnProofSampleFunction>
struct lambda_proof_sample_delegate
{
	OnProofSampleFunction on_new_proof_sample;

	lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) : on_new_proof_sample(on_new_proof_sample) { }

	template<typename ProofCalculus, typename Canonicalizer>
	inline void operator() (const theory<ProofCalculus, Canonicalizer>& T,
			const typename ProofCalculus::Proof* test_proof, double log_probability)
	{
		typedef typename ProofCalculus::Language Formula;
		typedef typename Formula::Term Term;

		/* get the current value of the term used to introduce the existential quantifier */
		Term* current_term = test_proof->operands[2]->term;
		on_new_proof_sample(T, current_term, log_probability);
	}
};

template<typename OnProofSampleFunction>
lambda_proof_sample_delegate<OnProofSampleFunction> make_lambda_proof_sample_delegate(OnProofSampleFunction on_new_proof_sample) {
	return lambda_proof_sample_delegate<OnProofSampleFunction>(on_new_proof_sample);
}

template<typename Theory, typename ProofAxioms, typename Collector>
bool check_consistency(
	Theory& T, const ProofAxioms& proof_axioms, const Collector& collector)
{
	bool success = proof_axioms.check_proof_axioms(T);
	success &= proof_axioms.check_universal_eliminations(T, collector);
	success &= T.check_concept_axioms();
	success &= T.check_disjunction_introductions();
	success &= T.are_elements_provable(*debug_terminal_printer);
	success &= T.sets.check_freeable_sets();
	success &= T.sets.are_descendants_valid();
	success &= T.sets.are_set_sizes_valid();
	success &= T.sets.check_set_ids();
	success &= T.sets.are_provable_elements_valid();
	if (!observations_has_test_proof(T, collector)) {
		fprintf(stderr, "WARNING: `collector.test_proof` is not an observation in the theory.\n");
		success = false;
	}
	return success;
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction, typename... Args>
bool log_joint_probability_of_lambda(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		theory<ProofCalculus, Canonicalizer>& T_MAP,
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
//T.print_axioms(stderr, *debug_terminal_printer);
	Proof* new_proof = T.add_formula(existential, set_diff, new_constant, std::forward<Args>(add_formula_args)...);
	free(*existential); if (existential->reference_count == 0) free(existential);
	if (new_proof == nullptr) {
		return false;
	} else if (!proof_axioms.template add<false>(new_proof, set_diff.new_set_axioms, proof_prior)) {
		T.template remove_formula<true>(new_proof, set_diff);
		return false;
	}
	set_diff.clear();

	hash_map<const Formula*, Formula*> formula_map(128);
	if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_MAP, formula_map)) {
		T.template remove_formula<false>(new_proof, set_diff);
		proof_axioms.template subtract<false>(new_proof, set_diff.old_set_axioms, proof_prior);
		free(*new_proof); if (new_proof->reference_count == 0) free(new_proof);
		return false;
	}

	auto new_proof_sample_delegate = make_lambda_proof_sample_delegate(on_new_proof_sample);
	auto collector = make_provability_collector(T, proof_prior, new_proof, new_proof_sample_delegate);
	double max_log_probability = collector.internal_collector.current_log_probability;
	for (unsigned int t = 0; t < num_samples; t++) {
/*fprintf(stderr, "DEBUG: t = %u\n", t);
T.print_axioms(stderr, *debug_terminal_printer);
T.print_disjunction_introductions(stderr, *debug_terminal_printer);
if (!check_consistency(T, proof_axioms, collector)) exit(0);
if (t == 252)
	debug_flag = true;
else debug_flag = false;*/
		do_mh_step(T, proof_prior, proof_axioms, collector, collector.internal_collector.test_proof, t < num_samples / 4 ? 1.0 : 0.1);
		if (collector.internal_collector.current_log_probability > max_log_probability) {
			free(T_MAP); formula_map.clear();
			if (!theory<ProofCalculus, Canonicalizer>::clone(T, T_MAP, formula_map)) {
				T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
				proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
				free(*collector.internal_collector.test_proof);
				if (collector.internal_collector.test_proof->reference_count == 0)
					free(collector.internal_collector.test_proof);
				return false;
			}
			max_log_probability = collector.internal_collector.current_log_probability;
		}
	}
	T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
	free(*collector.internal_collector.test_proof);
	if (collector.internal_collector.test_proof->reference_count == 0)
		free(collector.internal_collector.test_proof);
print("Best theory while answering question:\n", stdout);
T_MAP.print_axioms(stdout, *debug_terminal_printer);
fprintf(stdout, "Theory log probability: %lf\n", max_log_probability);
	return true;
}

template<typename Term, typename OnProofSampleFunction>
struct proof_sample_delegate
{
	Term* term;
	OnProofSampleFunction on_new_proof_sample;

	proof_sample_delegate(Term* term, OnProofSampleFunction on_new_proof_sample) : term(term), on_new_proof_sample(on_new_proof_sample) { }

	template<typename ProofCalculus, typename Canonicalizer>
	inline void operator() (const theory<ProofCalculus, Canonicalizer>& T,
			const typename ProofCalculus::Proof* test_proof, double log_probability)
	{
		on_new_proof_sample(T, term, log_probability);
	}
};

template<typename Term, typename OnProofSampleFunction>
proof_sample_delegate<Term, OnProofSampleFunction> make_proof_sample_delegate(Term* term, OnProofSampleFunction on_new_proof_sample) {
	return proof_sample_delegate<Term, OnProofSampleFunction>(term, on_new_proof_sample);
}

struct proof_disjunction_nodes {
	array<instance> expected_constants;
	array<unsigned int> expected_operand_indices;
	unsigned int constant_position;
	unsigned int operand_position;

	proof_disjunction_nodes() : expected_constants(8), expected_operand_indices(4), constant_position(0), operand_position(0) { }

	inline void reset() {
		constant_position = 0;
		operand_position = 0;
	}

	inline void clear() {
		expected_constants.clear();
		expected_operand_indices.clear();
		reset();
	}

	static inline void free(proof_disjunction_nodes& initializer) {
		core::free(initializer.expected_constants);
		core::free(initializer.expected_operand_indices);
	}
};

inline bool init(proof_disjunction_nodes& nodes) {
	if (!array_init(nodes.expected_constants, 8)) {
		return false;
	} else if (!array_init(nodes.expected_operand_indices, 8)) {
		free(nodes.expected_constants);
		return false;
	}
	nodes.constant_position = 0;
	nodes.operand_position = 0;
	return true;
}

template<typename Formula>
inline bool get_proof_disjunction_nodes(nd_step<Formula>* proof, proof_disjunction_nodes& nodes)
{
	typedef typename Formula::TermType TermType;

	switch (proof->type) {
	case nd_step_type::IMPLICATION_INTRODUCTION:
		if (proof->operands[0]->type == nd_step_type::FALSITY_ELIMINATION
		 && proof->operands[0]->operands[0]->type == nd_step_type::NEGATION_ELIMINATION
		 && proof->operands[0]->operands[0]->operands[1] == proof->operands[1])
		{
			if (!nodes.expected_operand_indices.add(1)) return false;
		} else {
			if (!nodes.expected_operand_indices.add(0)) return false;
		}
		break;
	case nd_step_type::PROOF_BY_CONTRADICTION:
		if (proof->operands[0]->type == nd_step_type::NEGATION_ELIMINATION
		 && proof->operands[0]->operands[0]->type == nd_step_type::CONJUNCTION_ELIMINATION
		 && proof->operands[0]->operands[0]->operands[0] == proof->operands[1])
		{
			unsigned int index = proof->operands[0]->operands[0]->operands[1]->parameters[0] + 1;
			if (!nodes.expected_operand_indices.add(index)) return false;
		}
		break;

	case nd_step_type::EXISTENTIAL_INTRODUCTION:
		if (proof->operands[2]->term->type == TermType::CONSTANT) {
			if (!nodes.expected_constants.add(instance_constant(proof->operands[2]->term->constant)))
				return false;
		} else if (proof->operands[2]->term->type == TermType::NUMBER) {
			if (!nodes.expected_constants.add(instance_number(proof->operands[2]->term->number)))
				return false;
		} else if (proof->operands[2]->term->type == TermType::STRING) {
			if (!nodes.expected_constants.add(instance_string(&proof->operands[2]->term->str)))
				return false;
		}
		break;

	case nd_step_type::DISJUNCTION_INTRODUCTION:
		if (!nodes.expected_operand_indices.add(proof->operands[2]->parameter + 1))
			return false;
		break;

	case nd_step_type::AXIOM:
	case nd_step_type::DISJUNCTION_INTRODUCTION_LEFT:
	case nd_step_type::DISJUNCTION_INTRODUCTION_RIGHT:
	case nd_step_type::CONJUNCTION_INTRODUCTION:
	case nd_step_type::BETA_EQUIVALENCE:
	case nd_step_type::CONJUNCTION_ELIMINATION:
	case nd_step_type::CONJUNCTION_ELIMINATION_LEFT:
	case nd_step_type::CONJUNCTION_ELIMINATION_RIGHT:
	case nd_step_type::DISJUNCTION_ELIMINATION:
	case nd_step_type::IMPLICATION_ELIMINATION:
	case nd_step_type::BICONDITIONAL_INTRODUCTION:
	case nd_step_type::BICONDITIONAL_ELIMINATION_LEFT:
	case nd_step_type::BICONDITIONAL_ELIMINATION_RIGHT:
	case nd_step_type::NEGATION_ELIMINATION:
	case nd_step_type::FALSITY_ELIMINATION:
	case nd_step_type::UNIVERSAL_INTRODUCTION:
	case nd_step_type::UNIVERSAL_ELIMINATION:
	case nd_step_type::EXISTENTIAL_ELIMINATION:
	case nd_step_type::EQUALITY_ELIMINATION:
	case nd_step_type::COMPARISON_INTRODUCTION:
	case nd_step_type::INEQUALITY_INTRODUCTION:
	case nd_step_type::PARAMETER:
	case nd_step_type::TERM_PARAMETER:
	case nd_step_type::ARRAY_PARAMETER:
	case nd_step_type::FORMULA_PARAMETER:
	case nd_step_type::COUNT:
		break;
	}

	unsigned int operand_count;
	nd_step<Formula>* const* operands;
	proof->get_subproofs(operands, operand_count);
	for (unsigned int i = 0; i < operand_count; i++) {
		if (operands[i] == nullptr) continue;
		if (!get_proof_disjunction_nodes(operands[i], nodes))
			return false;
	}
	return true;
}

struct cached_proof_sampler {
	proof_disjunction_nodes prev_proof;
	bool is_prev_proof;

	cached_proof_sampler() : is_prev_proof(true) { }

	inline void reset() {
		prev_proof.reset();
		is_prev_proof = true;
	}

	inline void clear() {
		prev_proof.clear();
		is_prev_proof = true;
	}
};

template<typename Proof>
inline void visit_node(const Proof& proof, const cached_proof_sampler& visitor) { }

template<bool Negated, typename Term> constexpr bool visit_unary_atom(const Term* term, const cached_proof_sampler& visitor) { return true; }
template<bool Negated> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const cached_proof_sampler& visitor) { return true; }
template<typename Proof> constexpr bool visit_subset_axiom(const Proof& proof, const cached_proof_sampler& visitor) { return true; }
constexpr bool visit_existential_intro(const cached_proof_sampler& visitor) { return true; }
constexpr bool visit_negated_universal_intro(const cached_proof_sampler& visitor) { return true; }
constexpr bool visit_negated_conjunction(const cached_proof_sampler& visitor) { return true; }
constexpr bool visit_disjunction_intro(const cached_proof_sampler& visitor) { return true; }
inline void on_subtract_changes(cached_proof_sampler& visitor) { }

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const cached_proof_sampler& visitor) { return true; }

template<typename Theory, typename Formula, typename Proof>
constexpr bool on_undo_filter_constants(
		const Theory& T, const Formula* quantified, const typename Formula::Term* term, unsigned int variable,
		const array_map<unsigned int, Proof*>& set_definitions, const cached_proof_sampler& visitor)
{ return true; }

template<typename ProofCalculus, typename Canonicalizer>
inline bool get_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		const array_map<unsigned int, typename ProofCalculus::Proof*>& set_definitions,
		cached_proof_sampler& sampler)
{
	if (!get_constants(T, formula, variable, constants, set_definitions))
		return false;

	if (sampler.is_prev_proof) {
		if (sampler.prev_proof.constant_position == sampler.prev_proof.expected_constants.length)
			/* this is currently possible when the new proof has more existential
			   introductions than the old proof, for example if the proof of the
			   expression `a(b)` chose a different set definition of `a` than in
			   the old proof */
			return false;

		unsigned int index;
		const instance& key = sampler.prev_proof.expected_constants[sampler.prev_proof.constant_position];
		for (index = 0; index < constants.length; index++) {
			if (key.type == instance_type::CONSTANT && key.constant >= T.new_constant_offset
			 && (T.ground_concepts[key.constant - T.new_constant_offset].types.keys == nullptr
			  || key.constant - T.new_constant_offset >= T.ground_concept_capacity)
			 && constants[index].type == instance_type::ANY)
			{
				break;
			} else if (key == constants[index]) {
				break;
			}
		}
		if (index == constants.length) {
			sampler.is_prev_proof = false;
		} else {
			swap(constants[0], constants[index]);
			sampler.prev_proof.constant_position++;
		}
	}
	return true;
}

template<typename Formula>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, cached_proof_sampler& sampler)
{
	if (!filter_operands(formula, indices))
		return false;

	if (sampler.is_prev_proof) {
		if (sampler.prev_proof.operand_position == sampler.prev_proof.expected_operand_indices.length) {
			fprintf(stderr, "filter_operands ERROR: `cached_proof_sampler.prev_proof` has no further `expected_operand_indices`.\n");
			exit(EXIT_FAILURE);
		}

		unsigned int index = indices.index_of(sampler.prev_proof.expected_operand_indices[sampler.prev_proof.operand_position]);
		if (index == indices.length) {
			sampler.is_prev_proof = false;
		} else {
			swap(indices[0], indices[index]);
			sampler.prev_proof.operand_position++;
		}
	}
	return true;
}

template<bool Contradiction, typename ProofCalculus, typename Canonicalizer>
bool is_impossible(
		typename ProofCalculus::Language* formula,
		const theory<ProofCalculus, Canonicalizer>& T,
		cached_proof_sampler& sampler)
{
	return is_impossible<Contradiction>(formula, T);
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, const instance& constant, cached_proof_sampler& sampler) {
	return true;
}

template<typename Formula>
inline bool inconsistent_constant(const Formula* formula, unsigned int index, cached_proof_sampler& sampler) {
	sampler.is_prev_proof = false;
	return false;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, cached_proof_sampler& sampler) {
	sampler.prev_proof.constant_position--;
}

template<typename Formula>
inline void finished_operand_indices(const Formula* formula, cached_proof_sampler& sampler) {
	sampler.prev_proof.operand_position--;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int min_set_size, unsigned int max_set_size,
		cached_proof_sampler& sampler)
{ }

template<typename Proof>
constexpr bool on_new_size_axiom(
		Proof* new_size_axiom,
		const cached_proof_sampler& visitor)
{
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(
		Proof* old_size_axiom,
		const cached_proof_sampler& visitor)
{ }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(
		unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		cached_proof_sampler& sampler)
{
	return compute_new_set_size(set_id, sets, out, min_set_size, max_set_size);
}

template<bool InitWithPrevProof, typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction, typename... Args>
bool log_joint_probability_of_lambda_by_linear_search_helper(
		theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		cached_proof_sampler& sampler, OnProofSampleFunction on_new_proof_sample)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	set_changes<Formula> set_diff;
	Proof* new_proof;
	if (InitWithPrevProof)
		new_proof = T.add_formula(logical_form, set_diff, new_constant, sampler);
	else new_proof = T.add_formula(logical_form, set_diff, new_constant);
	sampler.reset();
	if (new_proof == nullptr) {
		return false;
	} else if (!proof_axioms.template add<false>(new_proof, set_diff.new_set_axioms, proof_prior)) {
		T.template remove_formula<true>(new_proof, set_diff);
		return false;
	}
	set_diff.clear();

	if (!InitWithPrevProof) {
		sampler.clear();
		get_proof_disjunction_nodes(new_proof, sampler.prev_proof);
	}

	auto collector = make_provability_collector(T, proof_prior, new_proof, on_new_proof_sample);
	for (unsigned int t = 0; t < num_samples; t++) {
		//if (debug_flag) T.template print_axioms<true>(stderr, *debug_terminal_printer);
		//if (debug_flag) T.print_disjunction_introductions(stderr, *debug_terminal_printer);
		do_mh_step(T, proof_prior, proof_axioms, collector, collector.internal_collector.test_proof, 1.0);
	}
	//if (debug_flag) T.template print_axioms<true>(stderr, *debug_terminal_printer);
	//if (debug_flag) T.print_disjunction_introductions(stderr, *debug_terminal_printer);
	T.template remove_formula<false>(collector.internal_collector.test_proof, set_diff);
	proof_axioms.template subtract<false>(collector.internal_collector.test_proof, set_diff.old_set_axioms, proof_prior);
	free(*collector.internal_collector.test_proof);
	if (collector.internal_collector.test_proof->reference_count == 0)
		free(collector.internal_collector.test_proof);
	return true;
}

template<typename ProofCalculus, typename Canonicalizer, typename ProofPrior, typename OnProofSampleFunction>
bool log_joint_probability_of_lambda_by_linear_search(
		const theory<ProofCalculus, Canonicalizer>& T,
		ProofPrior& proof_prior, typename ProofPrior::PriorState& proof_axioms,
		typename ProofCalculus::Language* logical_form, unsigned int num_samples,
		OnProofSampleFunction on_new_proof_sample)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Term Term;
	typedef theory<ProofCalculus, Canonicalizer> Theory;
	typedef typename ProofPrior::PriorState PriorStateType;

	Formula* preprocessed = preprocess_formula(logical_form);
	if (preprocessed == nullptr)
		return false;

	array_map<unsigned int, unsigned int> variable_map(16);
	Formula* canonicalized = Canonicalizer::canonicalize(*preprocessed, variable_map);
	core::free(*preprocessed); if (preprocessed->reference_count == 0) core::free(preprocessed);
	if (canonicalized == nullptr)
		return false;

#if !defined(NDEBUG)
	typedef typename Formula::Type FormulaType;
	if (canonicalized->type != FormulaType::LAMBDA)
		fprintf(stderr, "log_joint_probability_of_lambda_by_linear_search WARNING: `canonicalized` is not a lambda expression.\n");
#endif

	array<instance> constants(T.ground_concept_capacity + T.constant_types.size + T.constant_negated_types.size + 1);
	if (!get_possible_constants(T, constants)) {
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		return false;
	}

	array<pair<unsigned int, bool>> new_name_events(8);
	array<const Formula*> arg1_of(4);
	array<const Formula*> arg2_of(4);
	const array_map<unsigned int, Proof*> set_definitions(4);
	if (!filter_constants_helper<false>(T, canonicalized->quantifier.operand, canonicalized->quantifier.variable, constants, set_definitions, new_name_events, arg1_of, arg2_of)) {
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		return false;
	}

	const std::minstd_rand prng_engine = core::engine;

	Term* var = Term::new_variable(canonicalized->quantifier.variable);
	if (var == nullptr) {
		free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
		return false;
	}
	cached_proof_sampler prev_proof;
	for (unsigned int i = 0; i < constants.length; i++)
	{
timer stopwatch;
		/* copy the theory */
		Theory& T_copy = *((Theory*) alloca(sizeof(Theory)));
		PriorStateType& proof_axioms_copy = *((PriorStateType*) alloca(sizeof(PriorStateType)));
		hash_map<const hol_term*, hol_term*> formula_map(128);
		Theory::clone(T, T_copy, formula_map);
		PriorStateType::clone(proof_axioms, proof_axioms_copy, formula_map);

		Term* constant = nullptr;
		unsigned int constant_id = 0;
		if (constants[i].type == instance_type::ANY) {
			constant_id = T_copy.get_free_concept_id();
			constant_id = T_copy.get_free_concept_id(constant_id + 100);
			if (!T_copy.try_init_concept(constant_id)) {
				free(*var); free(var);
				free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
				free(T_copy); free(proof_axioms_copy);
				return false;
			}

			constant = Formula::new_constant(constant_id);
			if (constant == nullptr) {
				free(*var); free(var);
				free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
				T_copy.free_concept_id(constant_id);
				free(T_copy); free(proof_axioms_copy);
				return false;
			}

		} else {
			if (constants[i].type == instance_type::CONSTANT) {
				constant = Term::new_constant(constants[i].constant);
			} else if (constants[i].type == instance_type::NUMBER) {
				constant = Term::new_number(constants[i].number);
			} else if (constants[i].type == instance_type::STRING) {
				constant = Term::new_string(*constants[i].str);
			}
			if (constant == nullptr) {
				free(*var); free(var);
				free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
				free(T_copy); free(proof_axioms_copy);
				return false;
			}
		}

		Formula* substituted = substitute(canonicalized->quantifier.operand, var, constant);
		if (substituted == nullptr) {
			free(*var); free(var);
			free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
			free(T_copy); free(proof_axioms_copy);
			return false;
		}

		core::engine = prng_engine;

		auto new_proof_sample_delegate = make_proof_sample_delegate(constant, on_new_proof_sample);
		if (prev_proof.prev_proof.expected_constants.length == 0 && prev_proof.prev_proof.expected_operand_indices.length == 0)
			log_joint_probability_of_lambda_by_linear_search_helper<false>(T_copy, proof_prior, proof_axioms_copy, substituted, num_samples, prev_proof, new_proof_sample_delegate);
		else log_joint_probability_of_lambda_by_linear_search_helper<true>(T_copy, proof_prior, proof_axioms_copy, substituted, num_samples, prev_proof, new_proof_sample_delegate);
total_reasoning += stopwatch.milliseconds();
fprintf(stderr, "consistency checking time: %llums, total reasoning time: %llums\n", consistency_checking_ms.load(), total_reasoning.load());
		free(*constant); if (constant->reference_count == 0) free(constant);
		free(*substituted); if (substituted->reference_count == 0) free(substituted);
		if (constants[i].type == instance_type::ANY && T_copy.ground_concepts[constant_id - T_copy.new_constant_offset].types.keys != nullptr)
			T_copy.try_free_concept_id(constant_id);

		free(T_copy); free(proof_axioms_copy);
	}
	free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
	return true;
}

#endif /* THEORY_H_ */
