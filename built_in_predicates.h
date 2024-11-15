#ifndef BUILT_IN_PREDICATES_H_
#define BUILT_IN_PREDICATES_H_

enum class built_in_predicates : unsigned int {
	ZERO = 0,
	UNKNOWN,
	ARG1,
	ARG2,
	ARG3,
	ARG1_OF,
	ARG2_OF,
	ARG3_OF,
	SIZE,
	INVERSE,
	HAS,
	EXIST,
	PRESENT,
	PRESENT_PROGRESSIVE,
	PRESENT_PERFECT,
	PRESENT_PERFECT_PROGRESSIVE,
	PAST,
	PAST_PROGRESSIVE,
	PAST_PERFECT,
	PAST_PERFECT_PROGRESSIVE,
	FUTURE,
	FUTURE_PROGRESSIVE,
	FUTURE_PERFECT,
	FUTURE_PERFECT_PROGRESSIVE,
	EMPTY,
	EMPTY_REF,
	WIDE_SCOPE,
	GAP,
	TRACE,
	GREATER,
	GREATEST,
	GREATER_THAN_OR_EQUAL,
	LEAST,
	NUMBER,
	LOCATE,
	NEGATIVE,
	HIGH_DEGREE,
	LOW_DEGREE,
	CARDINALITY,
	VALUE,
	LAMBDA_PREDICATE,
	PROPERTY,

	EQUALS, /* this is only used to refer to set definitions of the form `A=^[x]:f(x)` */
	SUBSET,
	MAXIMAL_SUBSET,
	SAME,
	OBJECT,
	NAME,
	NAMED_ENTITY, /* this is only used in parsing */
	MEASURE,
	REF,
	PLURAL_REF,
	ANIMATE,
	CAPABLE_OF,

	COUNT
};

/* WARNING: The below should preserve the order of the entries in the enum. */

unsigned int ARGS[] = {
	(unsigned int) built_in_predicates::ARG1,
	(unsigned int) built_in_predicates::ARG2,
	(unsigned int) built_in_predicates::ARG3
};

unsigned int ARGS_OF[] = {
	(unsigned int) built_in_predicates::ARG1_OF,
	(unsigned int) built_in_predicates::ARG2_OF,
	(unsigned int) built_in_predicates::ARG3_OF
};

template<built_in_predicates Arg>
struct invert_arg { };

template<> struct invert_arg<built_in_predicates::ARG1> { static constexpr built_in_predicates value = built_in_predicates::ARG1_OF; };
template<> struct invert_arg<built_in_predicates::ARG2> { static constexpr built_in_predicates value = built_in_predicates::ARG2_OF; };
template<> struct invert_arg<built_in_predicates::ARG3> { static constexpr built_in_predicates value = built_in_predicates::ARG3_OF; };
template<> struct invert_arg<built_in_predicates::ARG1_OF> { static constexpr built_in_predicates value = built_in_predicates::ARG1; };
template<> struct invert_arg<built_in_predicates::ARG2_OF> { static constexpr built_in_predicates value = built_in_predicates::ARG2; };
template<> struct invert_arg<built_in_predicates::ARG3_OF> { static constexpr built_in_predicates value = built_in_predicates::ARG3; };

unsigned int PAST_OR_PRESENT[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PAST
};

unsigned int PRESENT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE
};

unsigned int FUTURE_PREDICATES[] = {
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int PERFECT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int NON_PERFECT_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE
};

unsigned int PROGRESSIVE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int NON_PROGRESSIVE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT
};

unsigned int PRESENT_PROGRESSIVE_PREDICATE[] = {
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE
};

unsigned int PRESENT_PREDICATE[] = {
	(unsigned int) built_in_predicates::PRESENT
};

unsigned int PAST_PROGRESSIVE_PREDICATE[] = {
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE
};

unsigned int PAST_PREDICATE[] = {
	(unsigned int) built_in_predicates::PAST
};

unsigned int TENSE_PREDICATES[] = {
	(unsigned int) built_in_predicates::PRESENT,
	(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PRESENT_PERFECT,
	(unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST,
	(unsigned int) built_in_predicates::PAST_PROGRESSIVE,
	(unsigned int) built_in_predicates::PAST_PERFECT,
	(unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE,
	(unsigned int) built_in_predicates::FUTURE_PROGRESSIVE,
	(unsigned int) built_in_predicates::FUTURE_PERFECT,
	(unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE
};

unsigned int ONLY_ASPECT_PREDICATES[2][3] = {
	{(unsigned int) built_in_predicates::PRESENT_PROGRESSIVE, (unsigned int) built_in_predicates::PRESENT_PERFECT, (unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE},
	{(unsigned int) built_in_predicates::PAST_PROGRESSIVE, (unsigned int) built_in_predicates::PAST_PERFECT, (unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE}
};

unsigned int ASPECT_PREDICATES[2][4] = {
	{(unsigned int) built_in_predicates::PRESENT, (unsigned int) built_in_predicates::PRESENT_PROGRESSIVE, (unsigned int) built_in_predicates::PRESENT_PERFECT, (unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE},
	{(unsigned int) built_in_predicates::PAST, (unsigned int) built_in_predicates::PAST_PROGRESSIVE, (unsigned int) built_in_predicates::PAST_PERFECT, (unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE}
};

inline bool is_tense_predicate(unsigned int predicate) {
	return index_of(predicate, TENSE_PREDICATES, array_length(TENSE_PREDICATES)) < array_length(TENSE_PREDICATES);
}

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("<ZERO>", (unsigned int) built_in_predicates::ZERO)
		&& names.put("unknown", (unsigned int) built_in_predicates::UNKNOWN)
		&& names.put("arg1", (unsigned int) built_in_predicates::ARG1)
		&& names.put("arg2", (unsigned int) built_in_predicates::ARG2)
		&& names.put("arg3", (unsigned int) built_in_predicates::ARG3)
		&& names.put("arg1_of", (unsigned int) built_in_predicates::ARG1_OF)
		&& names.put("arg2_of", (unsigned int) built_in_predicates::ARG2_OF)
		&& names.put("arg3_of", (unsigned int) built_in_predicates::ARG3_OF)
		&& names.put("inverse", (unsigned int) built_in_predicates::INVERSE)
		&& names.put("has", (unsigned int) built_in_predicates::HAS)
		&& names.put("exist", (unsigned int) built_in_predicates::EXIST)
		&& names.put("size", (unsigned int) built_in_predicates::SIZE)
		&& names.put("present", (unsigned int) built_in_predicates::PRESENT)
		&& names.put("present_progressive", (unsigned int) built_in_predicates::PRESENT_PROGRESSIVE)
		&& names.put("present_perfect", (unsigned int) built_in_predicates::PRESENT_PERFECT)
		&& names.put("present_perfect_progressive", (unsigned int) built_in_predicates::PRESENT_PERFECT_PROGRESSIVE)
		&& names.put("past", (unsigned int) built_in_predicates::PAST)
		&& names.put("past_progressive", (unsigned int) built_in_predicates::PAST_PROGRESSIVE)
		&& names.put("past_perfect", (unsigned int) built_in_predicates::PAST_PERFECT)
		&& names.put("past_perfect_progressive", (unsigned int) built_in_predicates::PAST_PERFECT_PROGRESSIVE)
		&& names.put("future", (unsigned int) built_in_predicates::FUTURE)
		&& names.put("future_progressive", (unsigned int) built_in_predicates::FUTURE_PROGRESSIVE)
		&& names.put("future_perfect", (unsigned int) built_in_predicates::FUTURE_PERFECT)
		&& names.put("future_perfect_progressive", (unsigned int) built_in_predicates::FUTURE_PERFECT_PROGRESSIVE)
		&& names.put("empty", (unsigned int) built_in_predicates::EMPTY)
		&& names.put("empty_ref", (unsigned int) built_in_predicates::EMPTY_REF)
		&& names.put("W", (unsigned int) built_in_predicates::WIDE_SCOPE)
		&& names.put("GAP", (unsigned int) built_in_predicates::GAP)
		&& names.put("TRACE", (unsigned int) built_in_predicates::TRACE)
		&& names.put("greater", (unsigned int) built_in_predicates::GREATER)
		&& names.put("greatest", (unsigned int) built_in_predicates::GREATEST)
		&& names.put("≥", (unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL)
		&& names.put("least", (unsigned int) built_in_predicates::LEAST)
		&& names.put("=", (unsigned int) built_in_predicates::EQUALS)
		&& names.put("subset", (unsigned int) built_in_predicates::SUBSET)
		&& names.put("maximal_subset", (unsigned int) built_in_predicates::MAXIMAL_SUBSET)
		&& names.put("name", (unsigned int) built_in_predicates::NAME)
		&& names.put("same", (unsigned int) built_in_predicates::SAME)
		&& names.put("object", (unsigned int) built_in_predicates::OBJECT)
		&& names.put("number", (unsigned int) built_in_predicates::NUMBER)
		&& names.put("locate", (unsigned int) built_in_predicates::LOCATE)
		&& names.put("negative", (unsigned int) built_in_predicates::NEGATIVE)
		&& names.put("high_degree", (unsigned int) built_in_predicates::HIGH_DEGREE)
		&& names.put("low_degree", (unsigned int) built_in_predicates::LOW_DEGREE)
		&& names.put("cardinality", (unsigned int) built_in_predicates::CARDINALITY)
		&& names.put("value", (unsigned int) built_in_predicates::VALUE)
		&& names.put("lambda_predicate", (unsigned int) built_in_predicates::LAMBDA_PREDICATE)
		&& names.put("property", (unsigned int) built_in_predicates::PROPERTY)
		&& names.put("named_entity", (unsigned int) built_in_predicates::NAMED_ENTITY)
		&& names.put("measure", (unsigned int) built_in_predicates::MEASURE)
		&& names.put("ref", (unsigned int) built_in_predicates::REF)
		&& names.put("plural_ref", (unsigned int) built_in_predicates::PLURAL_REF)
		&& names.put("animate", (unsigned int) built_in_predicates::ANIMATE)
		&& names.put("capable_of", (unsigned int) built_in_predicates::CAPABLE_OF);
}

struct no_op {
	template<typename... Args>
	inline constexpr bool operator() (Args&&... args) const { return true; }
};

thread_local hol_term HOL_ZERO(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::ZERO);
thread_local hol_term HOL_EMPTY(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::EMPTY);
thread_local hol_term HOL_UNKNOWN(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::UNKNOWN);
thread_local hol_term HOL_GAP(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::GAP);
thread_local hol_term HOL_NAME(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::NAME);
thread_local hol_term HOL_EXIST(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::EXIST);
thread_local hol_term HOL_SAME(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::SAME);

#endif /* BUILT_IN_PREDICATES_H_ */
