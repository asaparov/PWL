#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include "array_view.h"
#include "morphology_en.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

thread_local hol_term HOL_ZERO(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::ZERO);
thread_local hol_term HOL_EMPTY(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::EMPTY);
thread_local hol_term HOL_UNKNOWN(hol_term_type::CONSTANT, (unsigned int) built_in_predicates::UNKNOWN);

enum class grammatical_conjunction : uint_fast8_t {
	NONE = 0,
	THAT,
	IF,
	WHETHER,
	BECAUSE,
	FOR,
	ANY
};

enum class auxiliary_flag : uint_fast8_t {
	NONE = 0,
	REQ_AUX,
	REQ_NO_AUX,
	AUX,

	ANY,
	NO_REQ_AUX, /* NONE, REQ_NO_AUX, or AUX */
	NO_REQ_NO_AUX, /* NONE, REQ_AUX, or AUX */
	NONE_OR_AUX, /* NONE or AUX */
	NONE_OR_REQ_AUX /* NONE or REQ_AUX */
};

enum class grammatical_flag : uint_fast8_t {
	IS_ADJUNCT,
	NULLABLE_SUBJECT,
	SUBORDINATE,
	PREPOSITION,
	PARTICLE,
	NEGATIVE,
	ADV,
	TION,
	LY,
	GENITIVE,
	COMMA,
	COUNT
};

enum class grammatical_flag_value : uint_fast8_t {
	FALSE = 0,
	TRUE,
	ANY
};

enum class correlator : uint_fast8_t {
	NONE = 0,
	BOTH,
	EITHER,
	NEITHER,
	ANY
};

enum class coordination : uint_fast8_t {
	NONE = 0,
	AND,
	OR,
	NOR,

	ANY,
	NOT_NONE
};

template<typename Stream>
static inline bool print_comma(bool& first, Stream& out) {
	if (first) {
		first = false;
		return true;
	} else {
		return print(',', out);
	}
}

template<typename Stream>
static inline bool print_feature(const char* feature_name, bool& first, grammatical_num number, Stream& out) {
	switch (number) {
	case grammatical_num::SINGULAR:
		return print_comma(first, out) && core::print(feature_name, out) && core::print(":sg", out);
	case grammatical_num::PLURAL:
		return print_comma(first, out) && core::print(feature_name, out) && core::print(":pl", out);
	case grammatical_num::ANY:
		return print_comma(first, out) && core::print(feature_name, out) && core::print(":*", out);
	case grammatical_num::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_num.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, grammatical_comparison comp, Stream& out) {
	switch (comp) {
	case grammatical_comparison::COMPARATIVE:
		return print_comma(first, out) && core::print("comp", out);
	case grammatical_comparison::SUPERLATIVE:
		return print_comma(first, out) && core::print("sup", out);
	case grammatical_comparison::ANY:
		return print_comma(first, out) && core::print("{comp,sup}", out);
	case grammatical_comparison::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_comparison.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, grammatical_conjunction number, Stream& out) {
	switch (number) {
	case grammatical_conjunction::THAT:
		return print_comma(first, out) && core::print("cnj:that", out);
	case grammatical_conjunction::IF:
		return print_comma(first, out) && core::print("cnj:if", out);
	case grammatical_conjunction::WHETHER:
		return print_comma(first, out) && core::print("cnj:whether", out);
	case grammatical_conjunction::BECAUSE:
		return print_comma(first, out) && core::print("cnj:because", out);
	case grammatical_conjunction::FOR:
		return print_comma(first, out) && core::print("cnj:for", out);
	case grammatical_conjunction::ANY:
		return print_comma(first, out) && core::print("cnj:*", out);
	case grammatical_conjunction::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_conjunction.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, auxiliary_flag aux, Stream& out) {
	switch (aux) {
	case auxiliary_flag::REQ_AUX:
		return print_comma(first, out) && core::print("req_aux", out);
	case auxiliary_flag::REQ_NO_AUX:
		return print_comma(first, out) && core::print("req_no_aux", out);
	case auxiliary_flag::AUX:
		return print_comma(first, out) && core::print("aux", out);
	case auxiliary_flag::ANY:
		return print_comma(first, out) && core::print("aux:*", out);
	case auxiliary_flag::NO_REQ_AUX:
		return print_comma(first, out) && core::print("~req_aux", out);
	case auxiliary_flag::NO_REQ_NO_AUX:
		return print_comma(first, out) && core::print("~req_no_aux", out);
	case auxiliary_flag::NONE_OR_AUX:
		return print_comma(first, out) && core::print("(aux)", out);
	case auxiliary_flag::NONE_OR_REQ_AUX:
		return print_comma(first, out) && core::print("(req_aux)", out);
	case auxiliary_flag::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized auxiliary_flag.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, grammatical_mood inf, Stream& out) {
	switch (inf) {
	case grammatical_mood::INDICATIVE:
		return true;
	case grammatical_mood::PAST_PARTICIPLE:
		return print_comma(first, out) && core::print("past_ptc", out);
	case grammatical_mood::PRESENT_PARTICIPLE:
		return print_comma(first, out) && core::print("pres_ptc", out);
	case grammatical_mood::SUBJUNCTIVE:
		return print_comma(first, out) && core::print("subjunctive", out);
	case grammatical_mood::BARE_INFINITIVE:
		return print_comma(first, out) && core::print("bare_inf", out);
	case grammatical_mood::TO_INFINITIVE:
		return print_comma(first, out) && core::print("to_inf", out);
	case grammatical_mood::ANY:
		return print_comma(first, out) && core::print("mood:*", out);
	case grammatical_mood::NOT_TO_INFINITIVE:
		return print_comma(first, out) && core::print("~to_inf", out);
	case grammatical_mood::NOT_SUBJUNCTIVE:
		return print_comma(first, out) && core::print("~subjunctive", out);
	case grammatical_mood::NOT_TO_INF_OR_SUBJ:
		return print_comma(first, out) && core::print("~{to_inf,subj}", out);
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_mood.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature_name(grammatical_flag flag, Stream& out) {
	switch (flag) {
	case grammatical_flag::IS_ADJUNCT:
		return core::print("adj", out);
	case grammatical_flag::NULLABLE_SUBJECT:
		return core::print("nullable_subj", out);
	case grammatical_flag::SUBORDINATE:
		return core::print("subord", out);
	case grammatical_flag::PREPOSITION:
		return core::print("prep", out);
	case grammatical_flag::PARTICLE:
		return core::print("particle", out);
	case grammatical_flag::NEGATIVE:
		return core::print("negative", out);
	case grammatical_flag::ADV:
		return core::print("adv", out);
	case grammatical_flag::TION:
		return core::print("tion", out);
	case grammatical_flag::LY:
		return core::print("ly", out);
	case grammatical_flag::GENITIVE:
		return core::print("gen", out);
	case grammatical_flag::COMMA:
		return core::print("comma", out);
	case grammatical_flag::COUNT: break;
	}
	fprintf(stderr, "print_feature_name ERROR: Unrecognized grammatical_flag.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, grammatical_flag flag, grammatical_flag_value value, Stream& out) {
	switch (value) {
	case grammatical_flag_value::TRUE:
		return print_comma(first, out) && print_feature_name(flag, out);
	case grammatical_flag_value::ANY:
		return print_comma(first, out) && print_feature_name(flag, out) && print('*', out);
	case grammatical_flag_value::FALSE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_flag_value.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(const char* prefix, bool& first, correlator corr, Stream& out) {
	switch (corr) {
	case correlator::BOTH:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":both", out);
	case correlator::EITHER:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":either", out);
	case correlator::NEITHER:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":neither", out);
	case correlator::ANY:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":*", out);
	case correlator::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized correlator.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(const char* prefix, bool& first, coordination coord, Stream& out) {
	switch (coord) {
	case coordination::AND:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":and", out);
	case coordination::OR:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":or", out);
	case coordination::NOR:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":nor", out);
	case coordination::ANY:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":*", out);
	case coordination::NOT_NONE:
		return print_comma(first, out) && core::print(prefix, out) && core::print(":~none", out);
	case coordination::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized coordination.\n");
	return false;
}

template<typename GrammaticalFlagType>
static inline bool intersect(GrammaticalFlagType& out, GrammaticalFlagType first, GrammaticalFlagType second)
{
	if (first == GrammaticalFlagType::ANY) {
		out = second;
		return true;
	} else if (second == GrammaticalFlagType::ANY) {
		out = first;
		return true;
	} else if (first != second) {
		return false;
	} else {
		out = first;
		return true;
	}
}

static inline bool intersect(coordination& out, coordination first, coordination second)
{
	if (first == coordination::ANY) {
		out = second;
		return true;
	} else if (second == coordination::ANY) {
		out = first;
		return true;
	} else if (first == coordination::NOT_NONE) {
		if (second == coordination::NONE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == coordination::NOT_NONE) {
		if (first == coordination::NONE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == second) {
		out = first;
		return true;
	} else {
		return false;
	}
}

static inline bool has_intersection(coordination first, coordination second) {
	coordination dummy;
	return intersect(dummy, first, second);
}

static inline bool intersect(auxiliary_flag& out, auxiliary_flag first, auxiliary_flag second)
{
	if (first == auxiliary_flag::ANY) {
		out = second;
		return true;
	} else if (second == auxiliary_flag::ANY) {
		out = first;
		return true;
	} else if (first == auxiliary_flag::NO_REQ_AUX) {
		if (second == auxiliary_flag::REQ_AUX) {
			return false;
		} else if (second == auxiliary_flag::NO_REQ_NO_AUX) {
			out = auxiliary_flag::NONE_OR_AUX;
			return true;
		} else if (second == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = auxiliary_flag::NONE;
			return true;
		} else {
			out = second;
			return true;
		}
	} else if (second == auxiliary_flag::NO_REQ_AUX) {
		if (first == auxiliary_flag::REQ_AUX) {
			return false;
		} else if (first == auxiliary_flag::NO_REQ_NO_AUX) {
			out = auxiliary_flag::NONE_OR_AUX;
			return true;
		} else if (first == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = auxiliary_flag::NONE;
			return true;
		} else {
			out = first;
			return true;
		}
	} else if (first == auxiliary_flag::NO_REQ_NO_AUX) {
		if (second == auxiliary_flag::REQ_NO_AUX) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == auxiliary_flag::NO_REQ_NO_AUX) {
		if (first == auxiliary_flag::REQ_NO_AUX) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == auxiliary_flag::NONE_OR_AUX) {
		if (second == auxiliary_flag::NONE || second == auxiliary_flag::AUX || second == auxiliary_flag::NONE_OR_AUX) {
			out = second;
			return true;
		} else if (second == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = auxiliary_flag::NONE;
			return true;
		} else {
			return false;
		}
	} else if (second == auxiliary_flag::NONE_OR_AUX) {
		if (first == auxiliary_flag::NONE || first == auxiliary_flag::AUX) {
			out = first;
			return true;
		} else if (first == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = auxiliary_flag::NONE;
			return true;
		} else {
			return false;
		}
	} else if (first == auxiliary_flag::NONE_OR_REQ_AUX) {
		if (second == auxiliary_flag::NONE || second == auxiliary_flag::REQ_AUX || second == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = second;
			return true;
		} else {
			return false;
		}
	} else if (second == auxiliary_flag::NONE_OR_REQ_AUX) {
		if (first == auxiliary_flag::NONE || first == auxiliary_flag::REQ_AUX || first == auxiliary_flag::NONE_OR_REQ_AUX) {
			out = first;
			return true;
		} else {
			return false;
		}
	} else if (first != second) {
		return false;
	} else {
		out = first;
		return true;
	}
}

static inline bool has_intersection(auxiliary_flag first, auxiliary_flag second) {
	auxiliary_flag dummy;
	return intersect(dummy, first, second);
}

static inline bool intersect(grammatical_mood& out, grammatical_mood first, grammatical_mood second)
{
	if (first == grammatical_mood::ANY) {
		out = second;
		return true;
	} else if (second == grammatical_mood::ANY) {
		out = first;
		return true;
	} else if (first == grammatical_mood::NOT_TO_INFINITIVE) {
		if (second == grammatical_mood::TO_INFINITIVE) {
			return false;
		} else if (second == grammatical_mood::NOT_SUBJUNCTIVE) {
			out = grammatical_mood::NOT_TO_INF_OR_SUBJ;
		} else {
			out = second;
		}
		return true;
	} else if (second == grammatical_mood::NOT_TO_INFINITIVE) {
		if (first == grammatical_mood::TO_INFINITIVE) {
			return false;
		} else if (first == grammatical_mood::NOT_SUBJUNCTIVE) {
			out = grammatical_mood::NOT_TO_INF_OR_SUBJ;
		} else {
			out = first;
		}
		return true;
	} else if (first == grammatical_mood::NOT_SUBJUNCTIVE) {
		if (second == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_mood::NOT_SUBJUNCTIVE) {
		if (first == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == grammatical_mood::NOT_TO_INF_OR_SUBJ) {
		if (second == grammatical_mood::TO_INFINITIVE || second == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_mood::NOT_TO_INF_OR_SUBJ) {
		if (first == grammatical_mood::TO_INFINITIVE || first == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first != second) {
		return false;
	} else {
		out = first;
		return true;
	}
}

static inline bool has_intersection(grammatical_mood first, grammatical_mood second) {
	grammatical_mood dummy;
	return intersect(dummy, first, second);
}

struct grammatical_flags {
	grammatical_num index_number;
	grammatical_num concord_number;
	grammatical_comparison comp;
	grammatical_conjunction cnj;
	grammatical_flag_value flags[(uint_fast8_t) grammatical_flag::COUNT] = { grammatical_flag_value::FALSE };
	correlator corr;
	correlator correlated_by;
	coordination coord;
	auxiliary_flag aux;
	grammatical_mood mood;
	bool aux_or_subjunctive_or_inf_or_to_inf;
	bool is_first_token_capital;

	grammatical_flags() :
			index_number(grammatical_num::NONE), concord_number(grammatical_num::NONE),
			comp(grammatical_comparison::NONE), cnj(grammatical_conjunction::NONE),
			corr(correlator::NONE), correlated_by(correlator::NONE), coord(coordination::NONE),
			aux(auxiliary_flag::NONE), mood(grammatical_mood::INDICATIVE),
			aux_or_subjunctive_or_inf_or_to_inf(false), is_first_token_capital(false) { }

	grammatical_flags(const grammatical_flags& src) {
		init(src);
	}

	inline void operator = (const grammatical_flags& src) {
		init(src);
	}

	static inline unsigned int hash(const grammatical_flags& key) {
		return default_hash(key.index_number)
			 ^ (3 * default_hash(key.concord_number))
			 ^ (7 * default_hash(key.comp))
			 ^ (31 * default_hash(key.cnj))
			 ^ (53 * default_hash(key.corr))
			 ^ (97 * default_hash(key.correlated_by))
			 ^ (193 * default_hash(key.coord))
			 ^ (389 * default_hash(key.flags, (uint_fast8_t) grammatical_flag::COUNT))
			 ^ (769 * default_hash(key.aux))
			 ^ (1543 * default_hash(key.mood))
			 ^ (3079 * default_hash(key.aux_or_subjunctive_or_inf_or_to_inf))
			 ^ (6151 * default_hash(key.is_first_token_capital));
	}

	static inline void swap(grammatical_flags& first, grammatical_flags& second) {
		core::swap(first.index_number, second.index_number);
		core::swap(first.concord_number, second.concord_number);
		core::swap(first.comp, second.comp);
		core::swap(first.cnj, second.cnj);
		core::swap(first.corr, second.corr);
		core::swap(first.correlated_by, second.correlated_by);
		core::swap(first.coord, second.coord);
		core::swap(first.aux, second.aux);
		core::swap(first.mood, second.mood);
		core::swap(first.aux_or_subjunctive_or_inf_or_to_inf, second.aux_or_subjunctive_or_inf_or_to_inf);
		core::swap(first.is_first_token_capital, second.is_first_token_capital);
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			core::swap(first.flags[i], second.flags[i]);
	}

	inline bool intersect_aux_or_subjunctive_or_inf_or_to_inf(bool first, bool second) {
		if (first && second) {
			if (aux == auxiliary_flag::AUX || mood == grammatical_mood::SUBJUNCTIVE || mood == grammatical_mood::BARE_INFINITIVE || mood == grammatical_mood::TO_INFINITIVE)
				aux_or_subjunctive_or_inf_or_to_inf = false;
			if (!has_intersection(aux, auxiliary_flag::AUX) && !has_intersection(mood, grammatical_mood::SUBJUNCTIVE)
			 && !has_intersection(mood, grammatical_mood::BARE_INFINITIVE) && !has_intersection(mood, grammatical_mood::TO_INFINITIVE))
				return false;
			aux_or_subjunctive_or_inf_or_to_inf = true;
		} else {
			aux_or_subjunctive_or_inf_or_to_inf = false;
		}
		return true;
	}

private:
	inline void init(const grammatical_flags& src) {
		index_number = src.index_number;
		concord_number = src.concord_number;
		comp = src.comp;
		cnj = src.cnj;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			flags[i] = src.flags[i];
		corr = src.corr;
		correlated_by = src.correlated_by;
		coord = src.coord;
		aux = src.aux;
		mood = src.mood;
		aux_or_subjunctive_or_inf_or_to_inf = src.aux_or_subjunctive_or_inf_or_to_inf;
		is_first_token_capital = src.is_first_token_capital;
	}
};

inline bool intersect(grammatical_flags& dst,
		const grammatical_flags& first,
		const grammatical_flags& second)
{
	if (!intersect(dst.index_number, first.index_number, second.index_number)
	 || !intersect(dst.concord_number, first.concord_number, second.concord_number)
	 || !intersect(dst.comp, first.comp, second.comp)
	 || !intersect(dst.cnj, first.cnj, second.cnj)
	 || !intersect(dst.corr, first.corr, second.corr)
	 || !intersect(dst.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(dst.coord, first.coord, second.coord)
	 || !intersect(dst.aux, first.aux, second.aux)
	 || !intersect(dst.mood, first.mood, second.mood))
		return false;
	dst.is_first_token_capital = (first.is_first_token_capital || second.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(dst.flags[i], first.flags[i], second.flags[i])) return false;
	return dst.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.aux_or_subjunctive_or_inf_or_to_inf, second.aux_or_subjunctive_or_inf_or_to_inf);
}

template<typename K, typename V>
struct static_pair {
	K key;
	V value;
};

template<typename Formula>
struct flagged_logical_form
{
	grammatical_flags flags;
	Formula* root;

	flagged_logical_form() { }

	flagged_logical_form(Formula& src) : root(&src) {
		src.reference_count++;
	}

	flagged_logical_form(const flagged_logical_form<Formula>& src) : flags(src.flags), root(src.root) {
		src.root->reference_count++;
	}

	~flagged_logical_form() { free_helper(); }

	inline void operator = (const flagged_logical_form<Formula>& src) {
		flags = src.flags;
		root = src.root;
		root->reference_count++;
	}

	static inline bool is_empty(const flagged_logical_form<Formula>& src) {
		return src.root == nullptr;
	}

	static inline unsigned int hash(const flagged_logical_form<Formula>& src) {
		return hasher<Formula>::hash(*src.root) ^ grammatical_flags::hash(src.flags);
	}

	static inline void move(
			const flagged_logical_form<Formula>& src,
			flagged_logical_form<Formula>& dst)
	{
		dst.root = src.root;
		dst.flags = src.flags;
	}

	static inline void swap(
			flagged_logical_form<Formula>& first,
			flagged_logical_form<Formula>& second)
	{
		core::swap(first.root, second.root);
		core::swap(first.flags, second.flags);
	}

	static inline void free(flagged_logical_form<Formula>& lf) {
		lf.free_helper();
	}

	inline void recompute_hash() { }

	enum class feature {
		EMPTY = 0,
		CONSTANT,
		PREDICATE,
		PREDICATE_ONLY,
		SET_DEFINITION,
		LEFT_ARG
	};

	enum class function_type {
		EMPTY = 0,
		IDENTITY,
		SELECT_RIGHT_CONJUNCT,
		REMOVE_RIGHT_CONJUNCT,
		SELECT_LEFT_CONJUNCT,
		SELECT_LEFT_CONJUNCT_AND_NEGATION,
		REMOVE_LEFT_CONJUNCT,
		REMOVE_LEFT_CONJUNCT_AND_NEGATION,
		REMOVE_SECOND_LEFT_CONJUNCT,
		SELECT_LEFT_CONJUNCT_IN_SET,
		REMOVE_LEFT_CONJUNCT_IN_SET,
		SELECT_RIGHT_CONJUNCT_IN_SET,
		REMOVE_RIGHT_CONJUNCT_IN_SET,
		SELECT_RIGHT_SUBSET_IN_SET,
		REQUIRE_NO_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE_OWN,
		REQUIRE_LEFT_PREDICATE_EXIST,
		REQUIRE_NO_LEFT_PREDICATE_EXIST,
		REMOVE_INVERSE,
		SELECT_RIGHT_ARG1_WITHOUT_HEAD,
		SELECT_RIGHT_ARG2_WITHOUT_HEAD,
		SELECT_RIGHT_ARG3_WITHOUT_HEAD,
		SELECT_RIGHT_ARG1_OF_WITHOUT_HEAD,
		SELECT_RIGHT_ARG2_OF_WITHOUT_HEAD,
		SELECT_RIGHT_ARG3_OF_WITHOUT_HEAD,
		SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE,
		SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE,
		SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE,
		SELECT_SINGLETON_ARG1_IN_SET_WITHOUT_HEAD_PREDICATIVE,
		SELECT_SINGLETON_ARG2_IN_SET_WITHOUT_HEAD_PREDICATIVE,
		SELECT_SINGLETON_ARG3_IN_SET_WITHOUT_HEAD_PREDICATIVE,
		REMOVE_FUTURE,
		REQUIRE_NO_FUTURE,
		REMOVE_PERFECT,
		REQUIRE_NO_PERFECT,
		REMOVE_PROGRESSIVE,
		REQUIRE_NO_PROGRESSIVE,
		REMOVE_NOT,
		REQUIRE_NO_EMPTY_REF,
		PREDICATE_ONLY,
		PREDICATE,
		PREDICATE_AND_TENSE,
		EMPTY_AND_TENSE,
		REQUIRE_PREDICATIVE_UNIVERSAL,
		REQUIRE_PREDICATIVE_EXISTENTIAL,
		REMOVE_PREDICATIVE_NOT,
		SELECT_PREDICATE_IN_SET,
		MARK_WIDE_SCOPE,
		REQUIRE_WIDE_SCOPE,
		REMOVE_WIDE_SCOPE,
		REQUIRE_NARROW_SCOPE,
		REQUIRE_CONSTANT_IN_SET,
		REQUIRE_NO_CONSTANT_IN_SET,
		REQUIRE_SINGLETON,
		SELECT_SECOND_LEFT_SET_CONJUNCT,
		SELECT_SECOND_LEFT_SET_CONJUNCT_ROOT,
		REMOVE_SECOND_LEFT_SET_CONJUNCT,
		SIZE,
		SET_SIZE,
		REQUIRE_LEFT_ARG1,
		REQUIRE_LAMBDA,
		REQUIRE_NO_LAMBDA,
		REMOVE_LEFT_PREDICATE,
		FACTOR,
		FACTOR_PREDICATIVE,
		SELECT_LEFT_PREDICATE_AND_TENSE,
		SET_PREDICATE_EMPTY,
		REQUIRE_CONJUNCTION,
		REQUIRE_BINARY_CONJUNCTION,
		REQUIRE_DISJUNCTION,
		REQUIRE_NEGATIVE_CONJUNCTION,
		SELECT_LEFT_OPERAND,
		REMOVE_LEFT_OPERAND,
		REPLACE_PREDICATIVE_SUBSET_WITH_EQUALITY,

		/* functions that modify grammatical features */
		ADD_SINGULAR,
		ADD_PLURAL,
		REQUIRE_SINGULAR,
		REQUIRE_PLURAL,
		TRY_REMOVE_NUMBER,
		ADD_CONCORD_SINGULAR,
		ADD_CONCORD_PLURAL,
		ADD_THAT,
		REMOVE_THAT,
		REQUIRE_NO_THAT,
		ADD_WHETHER,
		REMOVE_WHETHER,
		ADD_IF,
		REMOVE_IF,
		ADD_BECAUSE,
		REMOVE_BECAUSE,
		ADD_FOR,
		REMOVE_FOR,
		REQUIRE_NO_CONJUNCTION,
		ADD_IS_ADJUNCT,
		TRY_REMOVE_IS_ADJUNCT,
		REQUIRE_NOT_ADJUNCT,
		ADD_NULLABLE_SUBJECT,
		REMOVE_NULLABLE_SUBJECT,
		TRY_REMOVE_NULLABLE_SUBJECT,
		ADD_SUBORDINATE,
		REMOVE_SUBORDINATE,
		TRY_REMOVE_SUBORDINATE,
		REQUIRE_NO_SUBORDINATE,
		ADD_PREPOSITION,
		REQUIRE_PREPOSITION,
		REQUIRE_NO_PREPOSITION,
		ADD_PARTICLE,
		ADD_AUX,
		TRY_REMOVE_AUX,
		ADD_REQ_AUX,
		TRY_ADD_REQ_AUX,
		REQUIRE_NO_REQ_AUX,
		TRY_REMOVE_REQ_AUX,
		ADD_REQ_NO_AUX,
		REQUIRE_NO_REQ_NO_AUX,
		ADD_INFINITIVE,
		ADD_TO_INFINITIVE,
		REMOVE_TO_INFINITIVE,
		REQUIRE_TO_INFINITIVE,
		REQUIRE_NO_TO_INFINITIVE,
		ADD_SUBJUNCTIVE,
		REQUIRE_NO_SUBJUNCTIVE,
		REQUIRE_AUX_OR_SUBJUNCTIVE_OR_INFINITIVE_OR_TO_INFINITIVE,
		ADD_BOTH,
		ADD_EITHER,
		ADD_NEITHER,
		REMOVE_BOTH,
		REMOVE_EITHER,
		REMOVE_NEITHER,
		TRY_REMOVE_CORRELATOR,
		REQUIRE_NO_CORRELATOR,
		ADD_CORRELATED_BY_BOTH,
		ADD_CORRELATED_BY_EITHER,
		ADD_CORRELATED_BY_NEITHER,
		TRY_REMOVE_CORRELATED,
		REQUIRE_NOT_CORRELATED,
		ADD_PAST_PARTICIPLE,
		ADD_PRESENT_PARTICIPLE,
		TRY_REMOVE_PARTICIPLE,
		REQUIRE_PAST_PARTICIPLE,
		REQUIRE_PRESENT_PARTICIPLE,
		ADD_NEGATIVE,
		REQUIRE_NEGATIVE,
		ADD_ADV,
		REMOVE_ADV,
		TRY_REMOVE_ADV,
		REQUIRE_NO_ADV,
		ADD_TION,
		ADD_LY,
		ADD_GENITIVE,
		TRY_REMOVE_GENITIVE,
		REQUIRE_NO_GENITIVE,
		ADD_COMMA,
		REMOVE_COMMA,
		REQUIRE_NO_COMMA,
		ADD_AND,
		ADD_OR,
		ADD_NOR,
		REMOVE_AND,
		REMOVE_OR,
		REMOVE_NOR,
		REMOVE_COORD
	};

	static const static_pair<feature, const char*> FEATURE_NAMES[];
	static const static_pair<function_type, const char*> FUNCTION_NAMES[];

	struct function {
		function_type type;

		constexpr function(const function_type& type) : type(type) { }

		static inline unsigned int hash(const function& f) {
			return default_hash(f.type);
		}

		static inline bool is_empty(const function& f) {
			return f.type == function_type::EMPTY;
		}

		static inline void set_empty(function& f) {
			f.type = function_type::EMPTY;
		}

		inline bool operator == (const function& other) const {
			return type == other.type;
		}

		inline bool operator != (const function& other) const {
			return type != other.type;
		}

		inline bool operator < (const function& other) const {
			return type < other.type;
		}
	};

	template<typename Stream>
	static bool print(const feature& f, Stream& out) {
		for (unsigned int i = 0; i < array_length(FEATURE_NAMES); i++) {
			if (FEATURE_NAMES[i].key == f)
				return core::print(FEATURE_NAMES[i].value, out);
		}
		fprintf(stderr, "flagged_logical_form.print ERROR: Unrecognized semantic feature.\n");
		return false;
	}

	static bool interpret(feature& f, const string& name) {
		for (unsigned int i = 0; i < array_length(FEATURE_NAMES); i++) {
			if (name == FEATURE_NAMES[i].value) {
				f = FEATURE_NAMES[i].key;
				return true;
			}
		}
		core::print("flagged_logical_form.parse ERROR: Unrecognized semantic feature name '", stderr);
		core::print(name, stderr); core::print("'.\n", stderr);
		return false;
	}

	template<typename Stream>
	static bool print(const function& f, Stream& out) {
		for (unsigned int i = 0; i < array_length(FUNCTION_NAMES); i++) {
			if (FUNCTION_NAMES[i].key == f.type)
				return core::print(FUNCTION_NAMES[i].value, out);
		}
		fprintf(stderr, "flagged_logical_form.print ERROR: Unrecognized semantic transformation function.\n");
		return false;
	}

	static bool interpret(function& f, const string& name) {
		for (unsigned int i = 0; i < array_length(FUNCTION_NAMES); i++) {
			if (name == FUNCTION_NAMES[i].value) {
				f = FUNCTION_NAMES[i].key;
				return true;
			}
		}
		core::print("flagged_logical_form.parse ERROR: Unrecognized semantic transformation function name '", stderr);
		core::print(name, stderr); core::print("'.\n", stderr);
		return false;
	}

	static bool is_feature_pruneable(feature f)
	{
		switch (f) {
		case feature::CONSTANT:
		case feature::PREDICATE:
		case feature::PREDICATE_ONLY:
		case feature::SET_DEFINITION:
		case feature::LEFT_ARG:
			return true;
		case feature::EMPTY: break;
		}
		fprintf(stderr, "flagged_logical_form.is_feature_pruneable ERROR: Unrecognized semantic feature.\n");
		exit(EXIT_FAILURE);
	}

	static constexpr function default_function() {
		return function(function_type::EMPTY);
	}

private:
	inline void free_helper() {
		core::free(*root);
		if (root->reference_count == 0)
			core::free(root);
	}
};

template<typename Formula>
const static_pair<typename flagged_logical_form<Formula>::feature, const char*> flagged_logical_form<Formula>::FEATURE_NAMES[] = {
	{feature::EMPTY, "empty"},
	{feature::CONSTANT, "constant"},
	{feature::PREDICATE, "predicate"},
	{feature::PREDICATE_ONLY, "predicate_only"},
	{feature::SET_DEFINITION, "set_definition"},
	{feature::LEFT_ARG, "left_arg"}
};

template<typename Formula>
const static_pair<typename flagged_logical_form<Formula>::function_type, const char*> flagged_logical_form<Formula>::FUNCTION_NAMES[] = {
	{function_type::EMPTY, "empty"},
	{function_type::IDENTITY, "identity"},
	{function_type::SELECT_RIGHT_CONJUNCT, "select_right_conjunct"},
	{function_type::REMOVE_RIGHT_CONJUNCT, "remove_right_conjunct"},
	{function_type::SELECT_LEFT_CONJUNCT, "select_left_conjunct"},
	{function_type::SELECT_LEFT_CONJUNCT_AND_NEGATION, "select_left_conjunct_and_negation"},
	{function_type::REMOVE_LEFT_CONJUNCT, "remove_left_conjunct"},
	{function_type::REMOVE_LEFT_CONJUNCT_AND_NEGATION, "remove_left_conjunct_and_negation"},
	{function_type::REMOVE_SECOND_LEFT_CONJUNCT, "remove_second_left_conjunct"},
	{function_type::SELECT_LEFT_CONJUNCT_IN_SET, "select_left_conjunct_in_set"},
	{function_type::REMOVE_LEFT_CONJUNCT_IN_SET, "remove_left_conjunct_in_set"},
	{function_type::SELECT_RIGHT_CONJUNCT_IN_SET, "select_right_conjunct_in_set"},
	{function_type::REMOVE_RIGHT_CONJUNCT_IN_SET, "remove_right_conjunct_in_set"},
	{function_type::SELECT_RIGHT_SUBSET_IN_SET, "select_right_subset_in_set"},
	{function_type::REQUIRE_NO_INVERSE, "require_no_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE, "require_left_predicate_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN, "require_left_predicate_inverse_own"},
	{function_type::REQUIRE_LEFT_PREDICATE_EXIST, "require_left_predicate_exist"},
	{function_type::REQUIRE_NO_LEFT_PREDICATE_EXIST, "require_no_left_predicate_exist"},
	{function_type::REMOVE_INVERSE, "remove_inverse"},
	{function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD, "select_right_arg1_without_head"},
	{function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD, "select_right_arg2_without_head"},
	{function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD, "select_right_arg3_without_head"},
	{function_type::SELECT_RIGHT_ARG1_OF_WITHOUT_HEAD, "select_right_arg1_of_without_head"},
	{function_type::SELECT_RIGHT_ARG2_OF_WITHOUT_HEAD, "select_right_arg2_of_without_head"},
	{function_type::SELECT_RIGHT_ARG3_OF_WITHOUT_HEAD, "select_right_arg3_of_without_head"},
	{function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE, "select_right_arg1_without_head_predicative"},
	{function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE, "select_right_arg2_without_head_predicative"},
	{function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE, "select_right_arg3_without_head_predicative"},
	{function_type::SELECT_SINGLETON_ARG1_IN_SET_WITHOUT_HEAD_PREDICATIVE, "select_singleton_arg1_in_set_without_head_predicative"},
	{function_type::SELECT_SINGLETON_ARG2_IN_SET_WITHOUT_HEAD_PREDICATIVE, "select_singleton_arg2_in_set_without_head_predicative"},
	{function_type::SELECT_SINGLETON_ARG3_IN_SET_WITHOUT_HEAD_PREDICATIVE, "select_singleton_arg3_in_set_without_head_predicative"},
	{function_type::REMOVE_FUTURE, "remove_future"},
	{function_type::REQUIRE_NO_FUTURE, "require_no_future"},
	{function_type::REMOVE_PERFECT, "remove_perfect"},
	{function_type::REQUIRE_NO_PERFECT, "require_no_perfect"},
	{function_type::REMOVE_PROGRESSIVE, "remove_progressive"},
	{function_type::REQUIRE_NO_PROGRESSIVE, "require_no_progressive"},
	{function_type::REMOVE_NOT, "remove_not"},
	{function_type::REQUIRE_NO_EMPTY_REF, "require_no_empty_ref"},
	{function_type::PREDICATE_ONLY, "predicate_only"},
	{function_type::PREDICATE, "predicate"},
	{function_type::PREDICATE_AND_TENSE, "predicate_and_tense"},
	{function_type::EMPTY_AND_TENSE, "empty_and_tense"},
	{function_type::REQUIRE_PREDICATIVE_UNIVERSAL, "require_predicative_universal"},
	{function_type::REQUIRE_PREDICATIVE_EXISTENTIAL, "require_predicative_existential"},
	{function_type::REMOVE_PREDICATIVE_NOT, "remove_predicative_not"},
	{function_type::SELECT_PREDICATE_IN_SET, "select_predicate_in_set"},
	{function_type::MARK_WIDE_SCOPE, "mark_wide_scope"},
	{function_type::REQUIRE_WIDE_SCOPE, "require_wide_scope"},
	{function_type::REQUIRE_NARROW_SCOPE, "require_narrow_scope"},
	{function_type::REMOVE_WIDE_SCOPE, "remove_wide_scope"},
	{function_type::REQUIRE_CONSTANT_IN_SET, "require_constant_in_set"},
	{function_type::REQUIRE_NO_CONSTANT_IN_SET, "require_no_constant_in_set"},
	{function_type::REQUIRE_SINGLETON, "require_singleton"},
	{function_type::SELECT_SECOND_LEFT_SET_CONJUNCT, "select_second_left_set_conjunct"},
	{function_type::SELECT_SECOND_LEFT_SET_CONJUNCT_ROOT, "select_second_left_set_conjunct_root"},
	{function_type::REMOVE_SECOND_LEFT_SET_CONJUNCT, "remove_second_left_set_conjunct"},
	{function_type::SIZE, "size"},
	{function_type::SET_SIZE, "set_size"},
	{function_type::REQUIRE_LEFT_ARG1, "require_left_arg1"},
	{function_type::REQUIRE_LAMBDA, "require_lambda"},
	{function_type::REQUIRE_NO_LAMBDA, "require_no_lambda"},
	{function_type::REMOVE_LEFT_PREDICATE, "remove_left_predicate"},
	{function_type::FACTOR, "factor"},
	{function_type::FACTOR_PREDICATIVE, "factor_predicative"},
	{function_type::SELECT_LEFT_PREDICATE_AND_TENSE, "select_left_predicate_and_tense"},
	{function_type::SET_PREDICATE_EMPTY, "set_predicate_empty"},
	{function_type::REQUIRE_CONJUNCTION, "require_conjunction"},
	{function_type::REQUIRE_BINARY_CONJUNCTION, "require_binary_conjunction"},
	{function_type::REQUIRE_DISJUNCTION, "require_disjunction"},
	{function_type::REQUIRE_NEGATIVE_CONJUNCTION, "require_negative_conjunction"},
	{function_type::SELECT_LEFT_OPERAND, "select_left_operand"},
	{function_type::REMOVE_LEFT_OPERAND, "remove_left_operand"},
	{function_type::REPLACE_PREDICATIVE_SUBSET_WITH_EQUALITY, "replace_predicative_subset_with_equality"},
	{function_type::ADD_SINGULAR, "add_singular"},
	{function_type::ADD_PLURAL, "add_plural"},
	{function_type::REQUIRE_SINGULAR, "require_singular"},
	{function_type::REQUIRE_PLURAL, "require_plural"},
	{function_type::TRY_REMOVE_NUMBER, "try_remove_number"},
	{function_type::ADD_CONCORD_SINGULAR, "add_concord_singular"},
	{function_type::ADD_CONCORD_PLURAL, "add_concord_plural"},
	{function_type::ADD_THAT, "add_that"},
	{function_type::REMOVE_THAT, "remove_that"},
	{function_type::REQUIRE_NO_THAT, "require_no_that"},
	{function_type::ADD_WHETHER, "add_whether"},
	{function_type::REMOVE_WHETHER, "remove_whether"},
	{function_type::ADD_IF, "add_if"},
	{function_type::REMOVE_IF, "remove_if"},
	{function_type::ADD_BECAUSE, "add_because"},
	{function_type::REMOVE_BECAUSE, "remove_because"},
	{function_type::ADD_FOR, "add_for"},
	{function_type::REMOVE_FOR, "remove_for"},
	{function_type::REQUIRE_NO_CONJUNCTION, "require_no_conjunction"},
	{function_type::ADD_IS_ADJUNCT, "add_is_adjunct"},
	{function_type::TRY_REMOVE_IS_ADJUNCT, "try_remove_is_adjunct"},
	{function_type::REQUIRE_NOT_ADJUNCT, "require_not_adjunct"},
	{function_type::ADD_NULLABLE_SUBJECT, "add_nullable_subject"},
	{function_type::REMOVE_NULLABLE_SUBJECT, "remove_nullable_subject"},
	{function_type::TRY_REMOVE_NULLABLE_SUBJECT, "try_remove_nullable_subject"},
	{function_type::ADD_SUBORDINATE, "add_subordinate"},
	{function_type::REMOVE_SUBORDINATE, "remove_subordinate"},
	{function_type::TRY_REMOVE_SUBORDINATE, "try_remove_subordinate"},
	{function_type::REQUIRE_NO_SUBORDINATE, "require_no_subordinate"},
	{function_type::ADD_PREPOSITION, "add_preposition"},
	{function_type::REQUIRE_PREPOSITION, "require_preposition"},
	{function_type::REQUIRE_NO_PREPOSITION, "require_no_preposition"},
	{function_type::ADD_PARTICLE, "add_particle"},
	{function_type::ADD_AUX, "add_aux"},
	{function_type::TRY_REMOVE_AUX, "try_remove_aux"},
	{function_type::ADD_REQ_AUX, "add_req_aux"},
	{function_type::TRY_ADD_REQ_AUX, "try_add_req_aux"},
	{function_type::REQUIRE_NO_REQ_AUX, "require_no_req_aux"},
	{function_type::TRY_REMOVE_REQ_AUX, "try_remove_req_aux"},
	{function_type::ADD_REQ_NO_AUX, "add_req_no_aux"},
	{function_type::REQUIRE_NO_REQ_NO_AUX, "require_no_req_no_aux"},
	{function_type::ADD_INFINITIVE, "add_infinitive"},
	{function_type::ADD_TO_INFINITIVE, "add_to_infinitive"},
	{function_type::REMOVE_TO_INFINITIVE, "remove_to_infinitive"},
	{function_type::REQUIRE_TO_INFINITIVE, "require_to_infinitive"},
	{function_type::REQUIRE_NO_TO_INFINITIVE, "require_no_to_infinitive"},
	{function_type::ADD_SUBJUNCTIVE, "add_subjunctive"},
	{function_type::REQUIRE_NO_SUBJUNCTIVE, "require_no_subjunctive"},
	{function_type::REQUIRE_AUX_OR_SUBJUNCTIVE_OR_INFINITIVE_OR_TO_INFINITIVE, "require_aux_or_subjunctive_or_infinitive_or_to_infinitive"},
	{function_type::ADD_BOTH, "add_both"},
	{function_type::ADD_EITHER, "add_either"},
	{function_type::ADD_NEITHER, "add_neither"},
	{function_type::REMOVE_BOTH, "remove_both"},
	{function_type::REMOVE_EITHER, "remove_either"},
	{function_type::REMOVE_NEITHER, "remove_neither"},
	{function_type::TRY_REMOVE_CORRELATOR, "try_remove_correlator"},
	{function_type::REQUIRE_NO_CORRELATOR, "require_no_correlator"},
	{function_type::ADD_CORRELATED_BY_BOTH, "add_correlated_by_both"},
	{function_type::ADD_CORRELATED_BY_EITHER, "add_correlated_by_either"},
	{function_type::ADD_CORRELATED_BY_NEITHER, "add_correlated_by_neither"},
	{function_type::TRY_REMOVE_CORRELATED, "try_remove_correlated"},
	{function_type::REQUIRE_NOT_CORRELATED, "require_not_correlated"},
	{function_type::ADD_PAST_PARTICIPLE, "add_past_participle"},
	{function_type::ADD_PRESENT_PARTICIPLE, "add_present_participle"},
	{function_type::TRY_REMOVE_PARTICIPLE, "try_remove_participle"},
	{function_type::REQUIRE_PAST_PARTICIPLE, "require_past_participle"},
	{function_type::REQUIRE_PRESENT_PARTICIPLE, "require_present_participle"},
	{function_type::ADD_NEGATIVE, "add_negative"},
	{function_type::REQUIRE_NEGATIVE, "require_negative"},
	{function_type::ADD_ADV, "add_adv"},
	{function_type::REMOVE_ADV, "remove_adv"},
	{function_type::TRY_REMOVE_ADV, "try_remove_adv"},
	{function_type::REQUIRE_NO_ADV, "require_no_adv"},
	{function_type::ADD_TION, "add_tion"},
	{function_type::ADD_LY, "add_ly"},
	{function_type::ADD_GENITIVE, "add_genitive"},
	{function_type::TRY_REMOVE_GENITIVE, "try_remove_genitive"},
	{function_type::REQUIRE_NO_GENITIVE, "require_no_genitive"},
	{function_type::ADD_COMMA, "add_comma"},
	{function_type::REMOVE_COMMA, "remove_comma"},
	{function_type::REQUIRE_NO_COMMA, "require_no_comma"},
	{function_type::ADD_AND, "add_and"},
	{function_type::ADD_OR, "add_or"},
	{function_type::ADD_NOR, "add_nor"},
	{function_type::REMOVE_AND, "remove_and"},
	{function_type::REMOVE_OR, "remove_or"},
	{function_type::REMOVE_NOR, "remove_nor"},
	{function_type::REMOVE_COORD, "remove_coord"}
};

inline void initialize_any(grammatical_flags& flags) {
	flags.concord_number = grammatical_num::ANY;
	flags.index_number = grammatical_num::ANY;
	flags.comp = grammatical_comparison::NONE;
	flags.cnj = grammatical_conjunction::ANY;
	flags.corr = correlator::ANY;
	flags.correlated_by = correlator::ANY;
	flags.coord = coordination::ANY;
	flags.aux = auxiliary_flag::ANY;
	flags.mood = grammatical_mood::ANY;
	flags.aux_or_subjunctive_or_inf_or_to_inf = false;
	flags.is_first_token_capital = false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		flags.flags[i] = grammatical_flag_value::ANY;
}

template<typename Formula>
inline bool initialize_any(flagged_logical_form<Formula>& lf) {
	lf.root = (Formula*) malloc(sizeof(Formula));
	if (lf.root == nullptr) {
		fprintf(stderr, "initialize_any ERROR: Out of memory.\n");
		return false;
	}
	initialize_any(*lf.root);
	initialize_any(lf.flags);
	return true;
}

inline bool operator == (
		const grammatical_flags& first,
		const grammatical_flags& second)
{
	if (first.index_number != second.index_number
	 || first.concord_number != second.concord_number
	 || first.comp != second.comp
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by
	 || first.coord != second.coord
	 || first.aux != second.aux
	 || first.mood != second.mood
	 || first.aux_or_subjunctive_or_inf_or_to_inf != second.aux_or_subjunctive_or_inf_or_to_inf
	 || first.is_first_token_capital != second.is_first_token_capital)
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (first.flags[i] != second.flags[i]) return false;
	return true;
}

inline bool is_subset(
		const grammatical_flags& first,
		const grammatical_flags& second)
{
	grammatical_flags dst;
	if (!intersect(dst, first, second))
		return false;
	dst.is_first_token_capital = first.is_first_token_capital;
	return (dst == first);
}

template<typename Formula>
inline bool operator == (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	/* `first` may be uninitialized, but `second` must be initialized */
	if (first.root == nullptr)
		return false;
	return (first.flags == second.flags)
		&& ((first.root == second.root) || (*first.root == *second.root));
}

inline bool operator != (
		const grammatical_flags& first,
		const grammatical_flags& second)
{
	if (first.index_number != second.index_number
	 || first.concord_number != second.concord_number
	 || first.comp != second.comp
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by
	 || first.coord != second.coord
	 || first.aux != second.aux
	 || first.mood != second.mood
	 || first.aux_or_subjunctive_or_inf_or_to_inf != second.aux_or_subjunctive_or_inf_or_to_inf
	 || first.is_first_token_capital != second.is_first_token_capital)
		return true;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (first.flags[i] != second.flags[i]) return true;
	return false;
}

template<typename Formula>
inline bool operator != (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	return (first.flags != second.flags)
		|| ((first.root != second.root) && (*first.root != *second.root));
}

template<typename Formula>
inline bool equivalent(
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{

	Formula* first_logical_form = relabel_variables(first.root);
	if (first_logical_form == nullptr)
		return false;

	Formula* second_logical_form = relabel_variables(second.root);
	if (second_logical_form == nullptr) {
		free(*first_logical_form); if (first_logical_form->reference_count == 0) free(first_logical_form);
		return false;
	}

/* TODO: for debugging; remove this */
//print("first_logical_form:  ", stderr); print(*first_logical_form, stderr, *debug_terminal_printer); print('\n', stderr);
//print("second_logical_form: ", stderr); print(*second_logical_form, stderr, *debug_terminal_printer); print('\n', stderr);
	bool result = (first_logical_form == second_logical_form || *first_logical_form == *second_logical_form);
	free(*first_logical_form); if (first_logical_form->reference_count == 0) free(first_logical_form);
	free(*second_logical_form); if (second_logical_form->reference_count == 0) free(second_logical_form);
	return result;
}

template<typename Stream>
inline bool print(const grammatical_flags& flags, Stream& out)
{
	static grammatical_flags ZERO_FLAGS;
	if (flags == ZERO_FLAGS) return true;

	bool first = true;
	if (!print('[', out)
	 || !print_feature("index", first, flags.index_number, out)
	 || !print_feature("concord", first, flags.concord_number, out)
	 || !print_feature(first, flags.comp, out)
	 || !print_feature(first, flags.cnj, out)
	 || !print_feature("corr:", first, flags.corr, out)
	 || !print_feature("corr_by:", first, flags.correlated_by, out)
	 || !print_feature("coord", first, flags.coord, out)
	 || !print_feature(first, flags.aux, out)
	 || !print_feature(first, flags.mood, out))
		return false;

	if (flags.aux_or_subjunctive_or_inf_or_to_inf) {
		if ((!first && !print(',', out)) || !print("aux_or_subjunctive_or_inf_or_to_inf", out)) return false;
	}

	if (flags.is_first_token_capital) {
		if ((!first && !print(',', out)) || !print("cap", out)) return false;
	}

	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!print_feature(first, (grammatical_flag) i, flags.flags[i], out)) return false;
	return print(']', out);
}

template<typename Formula, typename Stream, typename... Printer>
inline bool print(
		const flagged_logical_form<Formula>& lf,
		Stream& out, Printer&&... printer)
{
	return print(*lf.root, out, std::forward<Printer>(printer)...) && print(lf.flags, out);
}

/* TODO: take into account the semantic prior during parsing */
template<bool Complete>
constexpr double log_probability(const flagged_logical_form<hol_term>& lf) {
	return 0.0;
}

inline bool parse_flag(
		char* str, unsigned int length,
		grammatical_flags& flags,
		hash_map<string, unsigned int>& names,
		position current)
{
	if (compare_strings("index:sg", str, length)) {
		flags.index_number = grammatical_num::SINGULAR;
	} else if (compare_strings("index:pl", str, length)) {
		flags.index_number = grammatical_num::PLURAL;
	} else if (compare_strings("index:*", str, length)) {
		flags.index_number = grammatical_num::ANY;
	} else if (compare_strings("concord:sg", str, length)) {
		flags.concord_number = grammatical_num::SINGULAR;
	} else if (compare_strings("concord:pl", str, length)) {
		flags.concord_number = grammatical_num::PLURAL;
	} else if (compare_strings("concord:*", str, length)) {
		flags.concord_number = grammatical_num::ANY;
	} else {
		read_error("Unrecognized grammatical flag label", current);
		return false;
	}
	return true;
}

bool parse(char* str, unsigned int length,
		grammatical_flags& flags,
		hash_map<string, unsigned int>& names,
		position start = position(1, 1))
{
	position current = start;
	unsigned int token_start = 0;
	for (unsigned int i = 0; i < length; i++) {
		if (str[i] == ',') {
			if (!parse_flag(str + token_start, i - token_start, flags, names, current))
				return false;
			token_start = i + 1;
		}

		if (str[i] == '\r' && i + 1 < length && str[i + 1] == '\n')
			i++;
		if (str[i] == '\r' || str[i] == '\n') {
			current.line++;
			current.column = 1;
		} else {
			current.column++;
		}
	}

	return parse_flag(str + token_start, length - token_start, flags, names, current);
}

template<typename Formula>
inline bool parse(char* str, unsigned int length,
		flagged_logical_form<Formula>& logical_form,
		hash_map<string, unsigned int>& names,
		position start = position(1, 1))
{
	static grammatical_flags ZERO_FLAGS;
	logical_form.flags = ZERO_FLAGS;

	/* first check if there are any grammatical flag labels */
	if (str[length - 1] == ']') {
		/* find the corresponding '[' */
		unsigned int index = last_index_of('[', str, length - 1);
		if (index != static_cast<unsigned int>(-1)) {
			/* appropriately shift the `start` position */
			position current = start;
			for (unsigned int i = 0; i < index + 1; i++) {
				if (str[i] == '\r' && i + 1 < length && str[i + 1] == '\n')
					i++;
				if (str[i] == '\r' || str[i] == '\n') {
					current.line++;
					current.column = 1;
				} else {
					current.column++;
				}
			}

			if (!parse(str + index + 1, length - index - 2, logical_form.flags, names, current))
				return false;
			length = index;
		}
	}

	logical_form.root = (Formula*) malloc(sizeof(Formula));
	if (logical_form.root == nullptr) {
		fprintf(stderr, "parse ERROR: Insufficient memory for `logical_form.root`.\n");
		return false;
	}

	memory_stream& in = *((memory_stream*) alloca(sizeof(memory_stream)));
	in.buffer = str;
	in.length = length;
	in.position = 0; in.shift = {0};
	if (!parse(in, *logical_form.root, names, start)) {
		free(logical_form.root);
		return false;
	}
	return true;
}

typedef sequence_distribution<token_distribution<double>> terminal_prior_type;
template<typename Formula> using hdp_grammar_type =
	hdp_grammar<rule_list_prior<terminal_prior<terminal_prior_type>, flagged_logical_form<Formula>>, flagged_logical_form<Formula>>;

template<typename Formula, typename Stream>
void print_nonterminal_hdps(hdp_grammar_type<Formula>& G, Stream& out,
		const string_map_scribe& terminal_printer, const string_map_scribe& nonterminal_printer)
{
	auto printers = make_pair<const string_map_scribe&, const string_map_scribe&>(terminal_printer, nonterminal_printer);
	for (unsigned int i = 0; i < G.nonterminals.length; i++) {
		if (G.nonterminals[i].rule_distribution.observations.sum == 0) continue;
		print(G.nonterminals[i].rule_distribution.type, out); print(' ', out);
		print(G.nonterminals[i].name, out); fprintf(out, " (%u) HDP: ", G.nonterminals[i].id);
		print(G.nonterminals[i].rule_distribution.sampler, out, terminal_printer, printers); print('\n', out);
		print(G.nonterminals[i].rule_distribution.h.alpha, G.nonterminals[i].rule_distribution.feature_count + 1, out); print('\n', out);
	}
}

template<typename Formula>
inline bool init(sequence& seq, const sentence<rooted_syntax_node<flagged_logical_form<Formula>>>& src)
{
	seq.tokens = (unsigned int*) malloc(sizeof(unsigned int) * src.length);
	if (seq.tokens == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for `sequence.tokens`.\n");
		return false;
	}
	for (unsigned int i = 0; i < src.length; i++)
		seq.tokens[i] = src.tokens[i].id;
	seq.length = src.length;
	return true;
}

struct number_parser_en
{
	static constexpr static_pair<int, const char*> ONES_NAMES[] = {
		{1, "one"},
		{2, "two"},
		{3, "three"},
		{4, "four"},
		{5, "five"},
		{6, "six"},
		{7, "seven"},
		{8, "eight"},
		{9, "nine"}
	};

	static constexpr static_pair<int, const char*> TEENS_NAMES[] = {
		{10, "ten"},
		{11, "eleven"},
		{12, "twelve"},
		{13, "thirteen"},
		{14, "fourteen"},
		{15, "fifteen"},
		{16, "sixteen"},
		{17, "seventeen"},
		{18, "eighteen"},
		{19, "nineteen"}
	};

	static constexpr static_pair<int, const char*> TENS_NAMES[] = {
		{20, "twenty"},
		{30, "thirty"},
		{40, "fourty"},
		{50, "fifty"},
		{60, "sixty"},
		{70, "seventy"},
		{80, "eighty"},
		{90, "ninety"}
	};

	static constexpr static_pair<int, const char*> LARGE_NUMBER_NAMES[] = {
		{1000000000, "billion"},
		{1000000, "million"},
		{1000, "thousand"}
	};

	unsigned int ones_name_ids[array_length(ONES_NAMES)];
	unsigned int teens_name_ids[array_length(TEENS_NAMES)];
	unsigned int tens_name_ids[array_length(TENS_NAMES)];
	unsigned int large_name_ids[array_length(LARGE_NUMBER_NAMES)];
	unsigned int zero_id, hundred_id, hyphen_id;

	number_parser_en(hash_map<string, unsigned int>& names) {
		for (unsigned int i = 0; i < array_length(ONES_NAMES); i++) {
			if (!get_token(ONES_NAMES[i].value, ones_name_ids[i], names))
				exit(EXIT_FAILURE);
		} for (unsigned int i = 0; i < array_length(TEENS_NAMES); i++) {
			if (!get_token(TEENS_NAMES[i].value, teens_name_ids[i], names))
				exit(EXIT_FAILURE);
		} for (unsigned int i = 0; i < array_length(TENS_NAMES); i++) {
			if (!get_token(TENS_NAMES[i].value, tens_name_ids[i], names))
				exit(EXIT_FAILURE);
		} for (unsigned int i = 0; i < array_length(LARGE_NUMBER_NAMES); i++) {
			if (!get_token(LARGE_NUMBER_NAMES[i].value, large_name_ids[i], names))
				exit(EXIT_FAILURE);
		}
		if (!get_token("zero", zero_id, names)
		 || !get_token("hundred", hundred_id, names)
		 || !get_token("-", hyphen_id, names))
			exit(EXIT_FAILURE);
	}
};

constexpr static_pair<int, const char*> number_parser_en::ONES_NAMES[];
constexpr static_pair<int, const char*> number_parser_en::TEENS_NAMES[];
constexpr static_pair<int, const char*> number_parser_en::TENS_NAMES[];
constexpr static_pair<int, const char*> number_parser_en::LARGE_NUMBER_NAMES[];

template<typename Formula>
bool get_unrecognized_terminals(
	array<unsigned int>& terminal_indices,
	const syntax_node<flagged_logical_form<Formula>>& syntax,
	const flagged_logical_form<Formula>& logical_form,
	unsigned int& left)
{
	const rule<flagged_logical_form<Formula>>& r = syntax.right;
	if (r.is_terminal()) {
		if (is_ambiguous(*logical_form.root)) {
			for (unsigned int i = 0; i < r.t.length; i++)
				if (!terminal_indices.add(left + i)) return false;
		}
		left += r.t.length;
		return true;
	}

	for (unsigned int i = 0; i < r.nt.length; i++) {
		flagged_logical_form<Formula>& transformed = *((flagged_logical_form<Formula>*) alloca(sizeof(flagged_logical_form<Formula>)));
		if (syntax.children[i] == nullptr) continue;
		if (!apply(r.nt.transformations[i], logical_form, transformed)) {
			return false;
		} else if (!get_unrecognized_terminals(terminal_indices, *syntax.children[i], transformed, left)) {
			free(transformed);
			return false;
		}
		free(transformed);
	}
	return true;
}

struct ambiguous_to_unknown_converter { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::ANY_CONSTANT_EXCEPT>::type* = nullptr>
inline hol_term* apply(hol_term* src, ambiguous_to_unknown_converter& converter) {
	HOL_UNKNOWN.reference_count++;
	return &HOL_UNKNOWN;
}

inline hol_term* ambiguous_to_unknown(hol_term* src)
{
	ambiguous_to_unknown_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct same_to_equals_converter { };

inline bool process_event_conjunct(
		unsigned int variable, hol_term* conjunct,
		hol_term*& arg1, hol_term*& arg2, unsigned int& predicate)
{
	if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
			&& conjunct->binary.left->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.left->binary.right->variable == variable)
	{
		arg1 = conjunct->binary.right;
	} else if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
			&& conjunct->binary.left->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.left->binary.right->variable == variable)
	{
		arg2 = conjunct->binary.right;
	} else if (conjunct->type == hol_term_type::UNARY_APPLICATION && conjunct->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.right->type == hol_term_type::VARIABLE && conjunct->binary.right->variable == variable)
	{
		bool is_tense = false;
		for (unsigned int i = 0; !is_tense && i < array_length(TENSE_PREDICATES); i++)
			if (TENSE_PREDICATES[i] == conjunct->binary.left->constant) is_tense = true;
		if (is_tense) return true;
		if (predicate != 0) return false;
		predicate = conjunct->binary.left->constant;
	} else {
		return false;
	}
	return true;
}

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, same_to_equals_converter& converter) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND || operand->array.length != 4)
		return default_apply<Type>(src, converter);

	hol_term* arg1 = nullptr;
	hol_term* arg2 = nullptr;
	unsigned int predicate = 0;
	for (unsigned int i = 0; i < operand->array.length; i++) {
		if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, predicate))
			return default_apply<Type>(src, converter);
	}

	if (arg1 == nullptr || arg2 == nullptr || predicate != (unsigned int) built_in_predicates::SAME)
		return default_apply<Type>(src, converter);

	hol_term* new_term = hol_term::new_equals(arg1, arg2);
	if (new_term == nullptr)
		return nullptr;
	arg1->reference_count++;
	arg2->reference_count++;
	return new_term;
}

inline hol_term* same_to_equals(hol_term* src)
{
	same_to_equals_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct exists_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, exists_remover& converter) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND || operand->array.length != 3)
		return default_apply<Type>(src, converter);

	hol_term* arg1 = nullptr;
	hol_term* arg2 = nullptr;
	unsigned int predicate = 0;
	for (unsigned int i = 0; i < operand->array.length; i++) {
		if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, predicate))
			return default_apply<Type>(src, converter);
	}

	if (arg1 == nullptr || arg2 != nullptr || predicate != (unsigned int) built_in_predicates::EXIST)
		return default_apply<Type>(src, converter);

	HOL_TRUE.reference_count++;
	return &HOL_TRUE;
}

inline hol_term* remove_exists(hol_term* src)
{
	exists_remover converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct arg_of_to_arg_converter { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EQUALS>::type* = nullptr>
inline hol_term* apply(hol_term* src, arg_of_to_arg_converter& converter) {
	if (src->binary.left->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.left->binary.left->type == hol_term_type::CONSTANT)
	{
		unsigned int index = index_of(src->binary.left->binary.left->constant, ARGS_OF, array_length(ARGS_OF));
		if (index < array_length(ARGS_OF)) {
			hol_term* left = apply(src->binary.left->binary.right, converter);
			if (left == nullptr) return nullptr;
			hol_term* right = apply(src->binary.right, converter);
			if (right == nullptr) {
				if (left != src->binary.left->binary.right) { free(*left); if (left->reference_count == 0) free(left); }
				return nullptr;
			}

			hol_term* new_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ARGS[index]), right), left);
			if (new_term == nullptr) {
				if (left != src->binary.left->binary.right) { free(*left); if (left->reference_count == 0) free(left); }
				if (right != src->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				return nullptr;
			}
			if (left == src->binary.left->binary.right) left->reference_count++;
			if (right == src->binary.right) right->reference_count++;
			return new_term;
		}
	} if (src->binary.right->type == hol_term_type::UNARY_APPLICATION
	   && src->binary.right->binary.left->type == hol_term_type::CONSTANT)
	{
		unsigned int index = index_of(src->binary.right->binary.left->constant, ARGS_OF, array_length(ARGS_OF));
		if (index < array_length(ARGS_OF)) {
			hol_term* right = apply(src->binary.right->binary.right, converter);
			if (right == nullptr) return nullptr;
			hol_term* left = apply(src->binary.left, converter);
			if (left == nullptr) {
				if (right != src->binary.right->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				return nullptr;
			}

			hol_term* new_term = hol_term::new_equals(right, hol_term::new_apply(hol_term::new_constant(ARGS[index]), left));
			if (new_term == nullptr) {
				if (right != src->binary.right->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				if (left != src->binary.left) { free(*left); if (left->reference_count == 0) free(left); }
				return nullptr;
			}
			if (right == src->binary.right->binary.right) right->reference_count++;
			if (left == src->binary.left) left->reference_count++;
			return new_term;
		}
	}
	return default_apply<Type>(src, converter);
}

inline hol_term* arg_of_to_arg(hol_term* src)
{
	arg_of_to_arg_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct tense_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, tense_remover& remover) {
	if (src->binary.left->type == hol_term_type::CONSTANT)
	{
		unsigned int index = index_of(src->binary.left->constant, TENSE_PREDICATES, array_length(TENSE_PREDICATES));
		if (index < array_length(TENSE_PREDICATES)) {
			HOL_TRUE.reference_count++;
			return &HOL_TRUE;
		}
	}
	return default_apply<Type>(src, remover);
}

inline hol_term* remove_tense(hol_term* src)
{
	tense_remover remover;
	hol_term* dst = apply(src, remover);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct object_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, object_remover& remover) {
	if (src->binary.left->type == hol_term_type::CONSTANT
	 && src->binary.left->constant == (unsigned int) built_in_predicates::OBJECT)
	{
		HOL_TRUE.reference_count++;
		return &HOL_TRUE;
	}
	return default_apply<Type>(src, remover);
}

inline hol_term* remove_object(hol_term* src)
{
	object_remover remover;
	hol_term* dst = apply(src, remover);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

template<typename Formula>
struct hdp_parser
{
	typedef flagged_logical_form<Formula> logical_form_type;

	morphology_en morph;
	hdp_grammar_type<Formula> G;
	const string** reverse_name_map;
	unsigned int name_count;

	number_parser_en number_parser;

	hdp_parser(unsigned int unknown_id,
			hash_map<string, unsigned int>& names,
			const char* morphology_filepath,
			const char* grammar_filepath) : reverse_name_map(nullptr), name_count(0), number_parser(names)
	{
		printf("Loading morphology data..."); fflush(stdout);
		if (!morph.initialize(names)
		 || !morphology_read(morph, names, morphology_filepath))
		{
			fflush(stdout); fprintf(stderr, "\nERROR: Unable to initialize morphology model.\n");
			exit(EXIT_FAILURE);
		}
		printf(" done.\n");

		printf("Loading grammar..."); fflush(stdout);
		if (!read_grammar(G, names, grammar_filepath)) {
			fflush(stdout); fprintf(stderr, "\nERROR: Unable to read grammar at '%s'.\n", grammar_filepath);
			exit(EXIT_FAILURE);
		}
		printf(" done.\n");
	}

	~hdp_parser() {
		if (reverse_name_map != nullptr)
			free(reverse_name_map);
	}

	bool invert_name_map(hash_map<string, unsigned int>& names) {
		if (!init_capitalization_map(morph, names))
			return false;
		if (reverse_name_map != NULL) free(reverse_name_map);
		reverse_name_map = invert(names);
		if (reverse_name_map == NULL) return false;
		name_count = names.table.size + 1;
		return true;
	}

	inline const string& map_to_string(unsigned int id) const {
		return *reverse_name_map[id];
	}

	bool train(
			const array<array_map<sentence<rooted_syntax_node<flagged_logical_form<Formula>>>, flagged_logical_form<Formula>>>& data,
			hash_map<string, unsigned int>& names,
			unsigned int iteration_count)
	{
		if (!init_capitalization_map(morph, names))
			return false;
		reverse_name_map = invert(names);
		const string** nonterminal_name_map = invert(G.nonterminal_names);
		if (reverse_name_map == NULL || nonterminal_name_map == NULL) {
			if (reverse_name_map != NULL) free(reverse_name_map);
			reverse_name_map = NULL;
			return false;
		}
		string_map_scribe terminal_printer = { reverse_name_map, names.table.size + 1 };
		string_map_scribe nonterminal_printer = { nonterminal_name_map, G.nonterminal_names.table.size + 1 };
/* TODO: for debugging; remove this */
debug_terminal_printer = &terminal_printer;
debug_nonterminal_printer = &nonterminal_printer;
detect_duplicate_logical_forms = true;
		/* construct the initial derivation trees (running the parser with an empty grammar) */
		rooted_syntax_node<logical_form_type>** syntax = (rooted_syntax_node<logical_form_type>**)
				calloc(data.length, sizeof(rooted_syntax_node<logical_form_type>*));
		unsigned int* order = (unsigned int*) malloc(sizeof(unsigned int) * data.length);
		if (syntax == NULL || order == NULL) {
			fprintf(stderr, "hdp_parser.train ERROR: Out of memory.\n");
			cleanup(data, nonterminal_name_map, syntax, order);
			return false;
		}
		for (unsigned int i = 0; i < data.length; i++) {
			syntax[i] = (rooted_syntax_node<logical_form_type>*) calloc(data[i].size, sizeof(rooted_syntax_node<logical_form_type>));
			if (syntax[i] == NULL) {
				fprintf(stderr, "hdp_parser.train ERROR: Out of memory.\n");
				cleanup(data, nonterminal_name_map, syntax, order);
				return false;
			}
		}
		for (unsigned int i = 0; i < data.length; i++) order[i] = i;
		shuffle(order, (unsigned int) data.length);
		for (unsigned int i = 0; i < data.length; i++) {
			unsigned int id = order[i];
			for (unsigned int j = 0; j < data[id].size; j++) {
				/* TODO: add a discourse model instead of treating each sentence in each passage as independent */
				const logical_form_type& logical_form = data[id].values[j];
				if (!is_empty(data[id].keys[j].derivation)) {
					/* this sentence has a fully specified derivation tree */
					logical_form_type any(HOL_ANY);
					initialize_any(any.flags);
					if (!is_parseable(*data[id].keys[j].derivation.tree, logical_form, G, morph, any, nonterminal_printer, terminal_printer, *this, data[id].keys[j].derivation.root))
					{
						fprintf(stderr, "sample ERROR: Derivation for example %u, sentence %u is not parseable: '", id, j);
						print(data[id].keys[j], stderr, terminal_printer); print("'\n", stderr);
						print(logical_form, stderr, terminal_printer); print("\n", stderr);
						cleanup(data, nonterminal_name_map, syntax, order);
						return false;
					}
					syntax[id][j] = data[id].keys[j].derivation;
				} else {
					sequence& seq = *((sequence*) alloca(sizeof(sequence)));
					if (!init(seq, data[id].keys[j])) {
						cleanup(data, nonterminal_name_map, syntax, order);
						return false;
					}
					auto sentence = tokenized_sentence<logical_form_type>(seq);
					free(seq);
					syntax[id][j].tree = (syntax_node<logical_form_type>*) malloc(sizeof(syntax_node<logical_form_type>));
					/* NOTE: sample can set syntax[id] to null */
					if (syntax[id][j].tree == NULL || !sample(syntax[id][j].tree, G, logical_form, sentence, morph, *this, syntax[id][j].root) || syntax[id][j].tree == NULL)
					{
						fprintf(stderr, "sample ERROR: Unable to sample derivation for example %u, sentence %u: '", id, j);
						print(data[id].keys[j], stderr, terminal_printer); print("'\n", stderr);
						print(logical_form, stderr, terminal_printer); print("\n", stderr);
						cleanup(data, nonterminal_name_map, syntax, order);
						return false;
					}

					print(logical_form, stdout, terminal_printer); print('\n', stdout);
					print(*syntax[id][j].tree, stdout, nonterminal_printer, terminal_printer, syntax[id][j].root); print("\n\n", stdout);
					fflush(stdout);
				}

				if (!add_tree(syntax[id][j].root, *syntax[id][j].tree, logical_form, G)) {
					cleanup(data, nonterminal_name_map, syntax, order);
					return false;
				}
			}
		}

		/* perform MCMC */
		for (unsigned int t = 0; t < iteration_count; t++) {
			shuffle(order, (unsigned int) data.length);
			for (unsigned int i = 0; i < data.length; i++) {
				unsigned int id = order[i];
				for (unsigned int j = 0; j < data[id].size; j++) {
					if (!is_empty(data[id].keys[j].derivation))
						/* do not resample training examples labeled with derivation trees */
						continue;
					const logical_form_type& logical_form = data[id].values[j];
					sequence& seq = *((sequence*) alloca(sizeof(sequence)));
					if (!init(seq, data[id].keys[j])) {
						cleanup(data, nonterminal_name_map, syntax, order);
						return false;
					}
					auto sentence = tokenized_sentence<logical_form_type>(seq);
					free(seq);
					resample(syntax[id][j].tree, G, logical_form, sentence, morph, *this, syntax[id][j].root);
				}
			}
			sample_grammar(G);
			fprintf(stdout, "Unnormalized log posterior probability: %lf\n",
					log_probability(G, syntax, data, *this));

			if (t % 1 == 0) {
				fprintf(stdout, "[iteration %u]\n", t);
				print_nonterminal_hdps(G, stdout, terminal_printer, nonterminal_printer);
				fprintf(stdout, "(seed = %u)\n", get_seed());
				fflush(stdout);
			}
		}

		/* cleanup */
		cleanup(data, nonterminal_name_map, syntax, order);
		return true;
	}

	template<unsigned int K, typename TheoryType>
	bool parse(const sentence<rooted_syntax_node<flagged_logical_form<Formula>>>& s,
			Formula** logical_forms, double* log_probabilities, unsigned int& parse_count,
			const TheoryType& T, array<sentence_token>& unrecognized)
	{
		static_assert(K > 0, "`K` must be at least 1.");
#if !defined(NDEBUG)
		if (reverse_name_map == NULL) {
			fprintf(stderr, "hdp_parser.parse ERROR: `hdp_parser.invert_name_map` must be called before `hdp_parser.parse`.\n");
			return false;
		}
#endif

		logical_form_type logical_form(HOL_ANY);
		syntax_node<logical_form_type>* parsed_syntax =
				(syntax_node<logical_form_type>*) alloca(K * sizeof(syntax_node<logical_form_type>));
		logical_form_type* logical_form_output =
				(logical_form_type*) alloca(K * sizeof(logical_form_type));
		sequence& seq = *((sequence*) alloca(sizeof(sequence)));
		if (!init(seq, s)) return false;
		auto sentence = tokenized_sentence<logical_form_type>(seq);
		free(seq);

/* TODO: for debugging; remove this */
const string** nonterminal_name_map = invert(G.nonterminal_names);
string_map_scribe terminal_printer = { reverse_name_map, name_count };
string_map_scribe nonterminal_printer = { nonterminal_name_map, G.nonterminal_names.table.size + 1 };
debug_terminal_printer = &terminal_printer;
debug_nonterminal_printer = &nonterminal_printer;
debug_flag = true;
		if (!::parse<true, false, K>(parsed_syntax, parse_count,
				logical_form, logical_form_output, G, sentence, morph, *this))
{
/* TODO: for debugging; remove this */
free(nonterminal_name_map);
			return false;
}

		array<unsigned int> ambiguous_terminal_indices(8);
		for (unsigned int i = 0; i < parse_count; i++) {
/* TODO: for debugging; remove this */
print(logical_form_output[i], stderr, terminal_printer); print('\n', stdout);
print(parsed_syntax[i], stderr, nonterminal_printer, terminal_printer, logical_form_output[i]); print("\n\n", stdout);
			double log_likelihood = log_probability(G, parsed_syntax[i], logical_form_output[i], *this);
			/* TODO: compute this prior */
			double log_prior = 0.0; //log_probability<true>(T, logical_form_output[i]);
			log_probabilities[i] = log_likelihood + log_prior;

			if (i == 0) {
				unsigned int left = 0;
				if (!get_unrecognized_terminals(ambiguous_terminal_indices, parsed_syntax[i], logical_form_output[i], left)
				 || !unrecognized.ensure_capacity(ambiguous_terminal_indices.length))
				{
					for (unsigned int j = 0; j < i; j++) {
						free(*logical_forms[j]); free(logical_forms[j]);
					} for (unsigned int j = 0; j < parse_count; j++) {
						free(parsed_syntax[j]);
						free(logical_form_output[j]);
					}
/* TODO: for debugging; remove this */
free(nonterminal_name_map);
					return false;
				}
			}

			logical_forms[i] = ambiguous_to_unknown(logical_form_output[i].root);
			if (logical_forms[i] == NULL) {
				fprintf(stderr, "hdp_parser.parse ERROR: Out of memory.\n");
				for (unsigned int j = 0; j < i; j++) {
					free(*logical_forms[j]); free(logical_forms[j]);
				} for (unsigned int j = 0; j < parse_count; j++) {
					free(parsed_syntax[j]);
					free(logical_form_output[j]);
				}
/* TODO: for debugging; remove this */
free(nonterminal_name_map);
				return false;
			}
		}
/* TODO: for debugging; remove this */
free(nonterminal_name_map);

		for (unsigned int index : ambiguous_terminal_indices)
			unrecognized[unrecognized.length++] = { s.tokens[index] };

		for (unsigned int j = 0; j < parse_count; j++) {
			free(parsed_syntax[j]);
			free(logical_form_output[j]);
		}
		return true;
	}

	bool add_definition(const sentence<rooted_syntax_node<flagged_logical_form<Formula>>>& s, Formula* definition, unsigned int new_constant)
	{
		/* first replace `unknown` in `definition` with `new_constant` */
		hol_term* constant = hol_term::new_constant(new_constant);
		if (constant == nullptr)
			return false;
		hol_term* new_definition = substitute(definition, &HOL_UNKNOWN, constant);
		free(*constant); if (constant->reference_count == 0) free(constant);
		if (new_definition == nullptr)
			return false;

		const logical_form_type logical_form(*new_definition);
		free(*new_definition); if (new_definition->reference_count == 0) free(new_definition);
		sequence& seq = *((sequence*) alloca(sizeof(sequence)));
		if (!init(seq, s)) return false;
		auto sentence = tokenized_sentence<logical_form_type>(seq);
		free(seq);
		syntax_node<logical_form_type>* syntax = (syntax_node<logical_form_type>*) malloc(sizeof(syntax_node<logical_form_type>));
		/* NOTE: sample can set syntax.tree to null */
		if (syntax == NULL || !sample(syntax, G, logical_form, sentence, morph, *this) || syntax == NULL)
		{
			fprintf(stderr, "hdp_parser.add_definition ERROR: Unable to sample derivation for new sentence.\n");
			return false;
		}

		const string** nonterminal_name_map = invert(G.nonterminal_names);
		if (nonterminal_name_map == nullptr) {
			free(*syntax); free(syntax);
			return false;
		}
		string_map_scribe terminal_printer = { reverse_name_map, name_count };
		string_map_scribe nonterminal_printer = { nonterminal_name_map, G.nonterminal_names.table.size + 1 };
		auto printer = get_printer(terminal_printer);
		print(logical_form, stdout, printer); print('\n', stdout);
		print(*syntax, stdout, nonterminal_printer, printer, 1); print("\n\n", stdout);
		fflush(stdout);
		free(nonterminal_name_map);

		if (!add_tree(1, *syntax, logical_form, G)) {
			free(*syntax); free(syntax);
			return false;
		}
		free(*syntax); free(syntax);

		/* TODO: do some MCMC iterations */
		return true;
	}

	template<typename Printer>
	constexpr Printer& get_printer(Printer& constant_printer) const {
		return constant_printer;
	}

private:
	void cleanup(
			const array<array_map<sentence<rooted_syntax_node<flagged_logical_form<Formula>>>, flagged_logical_form<Formula>>>& data,
			const string** nonterminal_name_map, rooted_syntax_node<logical_form_type>** syntax, unsigned int* order)
	{
		if (syntax != NULL) {
			for (unsigned int k = 0; k < data.length; k++) {
				if (syntax[k] == NULL) continue;
				for (unsigned int l = 0; l < data[k].size; l++)
					if (!is_empty(syntax[k][l])) free(syntax[k][l]);
				free(syntax[k]);
			}
			free(syntax);
		}
		if (order != NULL) free(order);
		if (reverse_name_map != NULL) free(reverse_name_map);
		if (nonterminal_name_map != NULL) free(nonterminal_name_map);
		reverse_name_map = NULL;
	}
};

/* Computes the log joint probability of the grammar and given derivations */
template<typename Formula>
double log_probability(
	hdp_grammar_type<Formula>& G,
	const rooted_syntax_node<flagged_logical_form<Formula>>* const* syntax,
	const array<array_map<sentence<rooted_syntax_node<flagged_logical_form<Formula>>>, flagged_logical_form<Formula>>>& data,
	const hdp_parser<Formula>& parser)
{
	typedef flagged_logical_form<Formula> logical_form_type;

	double score = 0.0;
	for (unsigned int i = 0; i < G.nonterminals.length; i++)
		score += log_probability(G.nonterminals[i].rule_distribution);
	for (unsigned int i = 0; i < data.length; i++) {
		for (unsigned int j = 0; j < data[i].size; j++) {
			const logical_form_type& logical_form = data[i].values[j];
			score += log_probability(G, *syntax[i][j].tree, logical_form, parser, syntax[i][j].root);
		}
	}
	return score;
}

template<typename Formula>
inline bool is_ambiguous(const flagged_logical_form<Formula>& exp) {
	return is_ambiguous(*exp.root);
}

template<typename BuiltInPredicates>
struct hol_non_head_constants {
	hol_term* constants[18];

	hol_non_head_constants() {
		constants[0] = hol_term::new_constant((unsigned int) BuiltInPredicates::UNKNOWN);
		constants[1] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG1);
		constants[2] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG2);
		constants[3] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG3);
		constants[4] = hol_term::new_constant((unsigned int) BuiltInPredicates::SIZE);
		constants[5] = hol_term::new_constant((unsigned int) BuiltInPredicates::WIDE_SCOPE);
		constants[6] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT);
		constants[7] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PROGRESSIVE);
		constants[8] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PERFECT);
		constants[9] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PERFECT_PROGRESSIVE);
		constants[10] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST);
		constants[11] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PROGRESSIVE);
		constants[12] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PERFECT);
		constants[13] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PERFECT_PROGRESSIVE);
		constants[14] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE);
		constants[15] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PROGRESSIVE);
		constants[16] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PERFECT);
		constants[17] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PERFECT_PROGRESSIVE);
	}

	~hol_non_head_constants() {
		for (unsigned int i = 0; i < array_length(constants); i++) {
			free(*constants[i]);
			if (constants[i]->reference_count == 0)
				free(constants[i]);
		}
	}

	static inline hol_non_head_constants<BuiltInPredicates>& get() {
		static hol_non_head_constants<BuiltInPredicates> constants;
		return constants;
	}

	static inline hol_term** get_terms() {
		return get().constants;
	}

	static inline void increment_terms() {
		for (unsigned int i = 0; i < array_length(get().constants); i++)
			get().constants[i]->reference_count++;
	}

	static inline unsigned int count() {
		return array_length(get().constants);
	}
};

template<typename BuiltInPredicates>
inline bool is_built_in(unsigned int constant) {
	return constant == (unsigned int) BuiltInPredicates::UNKNOWN
		|| constant == (unsigned int) BuiltInPredicates::ARG1
		|| constant == (unsigned int) BuiltInPredicates::ARG2
		|| constant == (unsigned int) BuiltInPredicates::ARG3
		|| constant == (unsigned int) BuiltInPredicates::ARG1_OF
		|| constant == (unsigned int) BuiltInPredicates::ARG2_OF
		|| constant == (unsigned int) BuiltInPredicates::ARG3_OF
		|| constant == (unsigned int) BuiltInPredicates::SIZE
		|| constant == (unsigned int) BuiltInPredicates::WIDE_SCOPE
		|| constant == (unsigned int) BuiltInPredicates::PRESENT
		|| constant == (unsigned int) BuiltInPredicates::PRESENT_PROGRESSIVE
		|| constant == (unsigned int) BuiltInPredicates::PRESENT_PERFECT
		|| constant == (unsigned int) BuiltInPredicates::PRESENT_PERFECT_PROGRESSIVE
		|| constant == (unsigned int) BuiltInPredicates::PAST
		|| constant == (unsigned int) BuiltInPredicates::PAST_PROGRESSIVE
		|| constant == (unsigned int) BuiltInPredicates::PAST_PERFECT
		|| constant == (unsigned int) BuiltInPredicates::PAST_PERFECT_PROGRESSIVE
		|| constant == (unsigned int) BuiltInPredicates::FUTURE
		|| constant == (unsigned int) BuiltInPredicates::FUTURE_PROGRESSIVE
		|| constant == (unsigned int) BuiltInPredicates::FUTURE_PERFECT
		|| constant == (unsigned int) BuiltInPredicates::FUTURE_PERFECT_PROGRESSIVE;
}

template<typename BuiltInPredicates>
inline hol_term* get_predicate_of_literal(
		hol_term* src, unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.right->type == hol_term_type::VARIABLE
	 && src->binary.right->variable == scope_variable
	 && src->binary.left->type == hol_term_type::CONSTANT
	 && !is_built_in<BuiltInPredicates>(src->binary.left->constant))
	{
		return src;
	}
	return nullptr;
}

inline hol_term* get_variable_definition_of_literal(
		hol_term* src, unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.left->type == hol_term_type::VARIABLE
	 && src->binary.left->variable == scope_variable)
	{
		return src;
	}
	return nullptr;
}

enum class head_position {
	NONE = 0,
	LEFT,
	RIGHT,
	ANY
};

struct head_index {
	head_position position;
	unsigned int index;

	inline bool compare(int conjunct_index) const {
		if (position == head_position::LEFT && conjunct_index >= 0) {
			return index == (unsigned int) conjunct_index;
		} else if (position == head_position::RIGHT && conjunct_index < 0) {
			return index == (unsigned int) (-conjunct_index) - 1;
		}
		return false;
	}
};

template<typename BuiltInPredicates>
inline hol_term* find_predicate(unsigned int head_variable, hol_term* conjunction, head_index& predicate_index)
{
	hol_term* expected_predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, hol_non_head_constants<BuiltInPredicates>::get_terms(), hol_non_head_constants<BuiltInPredicates>::count()),
			hol_term::new_variable(head_variable));
	if (expected_predicate == nullptr)
		return nullptr;
	hol_non_head_constants<built_in_predicates>::increment_terms();

	hol_term* predicate = nullptr;
	predicate_index.position = head_position::NONE;
	if (conjunction->type == hol_term_type::ANY_ARRAY) {
		for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.left.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::LEFT, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.left.operands[i], expected_predicate)) {
				predicate_index = {head_position::LEFT, i};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.right.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.right.operands[i], expected_predicate)) {
				predicate_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.any.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::ANY, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.any.operands[i], expected_predicate)) {
				predicate_index = {head_position::ANY, i};
			}
		}
	} else if (conjunction->type == hol_term_type::AND) {
		predicate_index = {head_position::NONE, 0};
		for (unsigned int i = 0; i < conjunction->array.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->array.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::LEFT, i};
				break;
			} else if (predicate == nullptr && has_intersection<BuiltInPredicates>(conjunction->array.operands[i], expected_predicate)) {
				predicate_index = {head_position::LEFT, i};
			}
		}
	} else {
		/* this conjunct could be singleton */
		predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction, head_variable);
		if (predicate != nullptr) {
			predicate_index = {head_position::LEFT, 0};
		} else if (predicate == nullptr && has_intersection<BuiltInPredicates>(conjunction, expected_predicate)) {
			predicate_index = {head_position::LEFT, 0};
		}
	}
	free(*expected_predicate); free(expected_predicate);
	return predicate;
}

template<typename BuiltInPredicates>
inline hol_term* find_variable_definition(unsigned int head_variable, hol_term* conjunction, head_index& definition_index)
{
	hol_term* expected_definition = hol_term::new_equals(
			hol_term::new_variable(head_variable), &HOL_ANY);
	if (expected_definition == nullptr)
		return nullptr;
	HOL_ANY.reference_count++;

	hol_term* variable_definition = nullptr;
	definition_index.position = head_position::NONE;
	if (conjunction->type == hol_term_type::ANY_ARRAY) {
		for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.left.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::LEFT, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.left.operands[i], expected_definition)) {
				definition_index = {head_position::LEFT, i};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.right.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.right.operands[i], expected_definition)) {
				definition_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.any.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::ANY, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.any.operands[i], expected_definition)) {
				definition_index = {head_position::ANY, i};
			}
		}
	} else if (conjunction->type == hol_term_type::AND) {
		definition_index = {head_position::NONE, 0};
		for (unsigned int i = 0; i < conjunction->array.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->array.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::LEFT, i};
				break;
			} else if (variable_definition == nullptr && has_intersection<BuiltInPredicates>(conjunction->array.operands[i], expected_definition)) {
				definition_index = {head_position::LEFT, i};
			}
		}
	} else {
		/* this conjunct could be singleton */
		variable_definition = get_variable_definition_of_literal(conjunction, head_variable);
		if (variable_definition != nullptr) {
			definition_index = {head_position::LEFT, 0};
		} else if (variable_definition == nullptr && has_intersection<BuiltInPredicates>(conjunction, expected_definition)) {
			definition_index = {head_position::LEFT, 0};
		}
	}
	free(*expected_definition); free(expected_definition);
	return variable_definition;
}

template<typename BuiltInPredicates, bool AllowVariableDefinitions = false, bool IncludeNegation = true>
inline void find_head(
		hol_term* src, hol_term*& head,
		head_index& predicate_index)
{
	head = nullptr;
	if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
		if (src->any.included == nullptr) {
			head = src;
		} else {
			find_head<BuiltInPredicates>(src->any.included, head, predicate_index);
			if (head != nullptr)
				head = src;
		}
	} else if (src->type == hol_term_type::EXISTS) {
		head = src;
		unsigned int head_variable = src->quantifier.variable;

		/* make sure this scope has no term `arg1(*)=x` or `arg2(*)=x` where `x` is the scope variable */
		hol_term wildcard;
		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_constant((unsigned int) BuiltInPredicates::ARG1), &HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		bool not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) { head = nullptr; return; }

		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_constant((unsigned int) BuiltInPredicates::ARG2), &HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) { head = nullptr; return; }

		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_constant((unsigned int) BuiltInPredicates::ARG3), &HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) { head = nullptr; return; }

		/* find the predicate */
		find_predicate<BuiltInPredicates>(head_variable, head->quantifier.operand, predicate_index);
		if (predicate_index.position == head_position::NONE) {
			if (AllowVariableDefinitions)
				find_variable_definition<BuiltInPredicates>(head_variable, head->quantifier.operand, predicate_index);

			if (!AllowVariableDefinitions || predicate_index.position == head_position::NONE) {
				head = nullptr;
				return;
			}
		}
	} else if (IncludeNegation && src->type == hol_term_type::NOT) {
		find_head<BuiltInPredicates>(src->unary.operand, head, predicate_index);
		if (head != nullptr)
			head = src;
	}
}

template<typename BuiltInPredicates>
inline void find_head_or_universal(
		hol_term* src, hol_term*& head,
		head_index& predicate_index)
{
	if ((src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) && src->any.included != nullptr) {
		find_head_or_universal<BuiltInPredicates>(src->any.included, head, predicate_index);
		if (head != nullptr) return;
	}

	hol_term* universal = hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, &HOL_ANY);
	if (universal == nullptr)
		return;
	HOL_ANY.reference_count++;

	if (has_intersection<BuiltInPredicates>(src, universal)) {
		free(*universal); free(universal);
		head = src; return;
	}
	free(*universal); free(universal);

	find_head<BuiltInPredicates>(src, head, predicate_index);
}

template<typename BuiltInPredicates, unsigned int Constant>
inline void find_head_or_unary_application(
		hol_term* src, hol_term*& head,
		head_index& predicate_index)
{
	if (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == Constant) {
		head = src;
	} else {
		find_head<BuiltInPredicates>(src, head, predicate_index);
	}
}

template<typename BuiltInPredicates>
struct predicative_head_finder {
	unsigned int lambda_variable;

	predicative_head_finder(unsigned int lambda_variable) : lambda_variable(lambda_variable) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		apply<true>(src, head, predicate_index);
	}

	inline bool operator == (const predicative_head_finder<BuiltInPredicates>& other) const {
		return lambda_variable == other.lambda_variable;
	}

	template<typename T>
	constexpr inline bool operator == (const T& other) const {
		return false;
	}

private:
	template<bool First>
	void apply(
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
			if (src->any.included == nullptr) {
				head = src;
			} else {
				apply<First>(src->any.included, head, predicate_index);
				if (head != nullptr)
					head = src;
			}
		} else if (src->type == hol_term_type::EXISTS) {
			unsigned int set_variable = src->quantifier.variable;
			hol_term* operand = src->quantifier.operand;
			hol_term* right = nullptr;
			if (operand->type == hol_term_type::ANY) {
				head = nullptr;
				return;
			} else if (operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::ANY_ARRAY || operand->any_array.oper == hol_term_type::AND)) {
				if (operand->any_array.right.length == 0) {
					head = nullptr;
					return;
				}
				right = operand->any_array.right.operands[operand->any_array.right.length - 1];
			} else if (operand->type == hol_term_type::AND) {
				if (operand->array.length == 1) {
					head = nullptr;
					return;
				}
				right = operand->array.operands[operand->array.length - 1];
			} else {
				head = nullptr;
				return;
			}

			while (right->type == hol_term_type::NOT)
				right = right->unary.operand;

			bool right_is_any = false;
			if (right->type == hol_term_type::ANY_RIGHT && right->any.included != nullptr) {
				right = right->any.included;
				right_is_any = true;
			} if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
				right = right->binary.right;
			} if (right->type == hol_term_type::ANY_RIGHT && right->any.included != nullptr) {
				right = right->any.included;
				right_is_any = true;
			}

			while (right->type == hol_term_type::NOT)
				right = right->unary.operand;

			hol_term* inner_right;
			unsigned int quantified_variable = 0;
			predicate_index = {head_position::LEFT, 0};
			if (right->type == hol_term_type::AND && right->array.length == 2 && right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
			 && right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[0]->binary.left->variable == set_variable
			 && right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
			{
				head = src;
				return;
			} else if (right->type == hol_term_type::FOR_ALL && right->quantifier.operand->type == hol_term_type::IF_THEN
					&& right->quantifier.operand->binary.left->type == hol_term_type::UNARY_APPLICATION && right->quantifier.operand->binary.left->binary.left->type == hol_term_type::VARIABLE
					&& right->quantifier.operand->binary.left->binary.left->variable == set_variable)
			{
				quantified_variable = right->quantifier.variable;
				inner_right = right->quantifier.operand->binary.right;
			} else if (right->type == hol_term_type::EXISTS && right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.length == 2
					&& right->quantifier.operand->array.operands[0]->type == hol_term_type::UNARY_APPLICATION && right->quantifier.operand->array.operands[0]->binary.left->type == hol_term_type::VARIABLE
					&& right->quantifier.operand->array.operands[0]->binary.left->variable == set_variable)
			{
				quantified_variable = right->quantifier.variable;
				inner_right = right->quantifier.operand->array.operands[1];
			} else if (right_is_any) {
				inner_right = right;
			} else {
				head = nullptr;
				return;
			}

			if (inner_right->type == hol_term_type::ANY_RIGHT && inner_right->any.included != nullptr)
				inner_right = inner_right->any.included;
			else if (inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
				inner_right = inner_right->binary.right;

			if (inner_right->type == hol_term_type::UNARY_APPLICATION
			 && inner_right->binary.left->type == hol_term_type::VARIABLE && inner_right->binary.left->variable == lambda_variable
			 && inner_right->binary.right->type == hol_term_type::VARIABLE && (quantified_variable == 0 || inner_right->binary.right->variable == quantified_variable))
			{
				head = src;
			} else {
				head = nullptr;
				return;
			}
		} else if (src->type == hol_term_type::NOT && First) {
			apply<false>(src->unary.operand, head, predicate_index);
			if (head != nullptr)
				head = src;
		} else {
			head = nullptr;
		}
	}
};

template<typename FindHead>
struct array_finder {
	FindHead& find_head;

	array_finder(FindHead& find_head) : find_head(find_head) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
			if (src->any.included == nullptr) {
				head = src;
			} else {
				head = nullptr;
			}
		} else if (src->type == hol_term_type::AND || src->type == hol_term_type::OR) {
			find_head(src->array.operands[0], head, predicate_index);
			if (head != nullptr)
				head = src;
		} else if (src->type == hol_term_type::ANY_ARRAY && (src->any_array.oper == hol_term_type::ANY_ARRAY || src->any_array.oper == hol_term_type::AND || src->any_array.oper == hol_term_type::OR)) {
			find_head(src->any_array.all, head, predicate_index);
			if (head != nullptr)
				head = src;
		} else {
			/* check if the "array" is a singleton */
			find_head(src, head, predicate_index);
		}
	}

	inline bool operator == (const array_finder<FindHead>& other) const {
		return find_head == other.find_head;
	}

	template<typename T>
	constexpr inline bool operator == (const T& other) const {
		return false;
	}
};

template<typename FindHead>
inline array_finder<FindHead> make_array_finder(FindHead& find_head) {
	return array_finder<FindHead>(find_head);
}

inline void find_root(hol_term* src, hol_term*& head, head_index& predicate_index) {
	head = src;
}

inline void find_zero(hol_term* src, hol_term*& head, head_index& predicate_index) {
	if (src == &HOL_ZERO) {
		head = src;
	} else {
		head = nullptr;
	}
}

struct node_finder {
	hol_term* node;

	node_finder(hol_term* node) : node(node) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		if (src == node) {
			head = src;
		} else {
			head = nullptr;
		}
	}

	inline bool operator == (const node_finder& other) const {
		return node == other.node;
	}

	template<typename T>
	constexpr inline bool operator == (const T& other) const {
		return false;
	}
};

bool prune_independent_siblings(
		array<unsigned int>& keep_variables,
		array<hol_term*>& siblings)
{
	array<unsigned int>* sibling_variables = (array<unsigned int>*) malloc(max((size_t) 1, sizeof(array<unsigned int>) * siblings.length));
	for (unsigned int i = 0; i < siblings.length; i++) {
		if (!array_init(sibling_variables[i], 8)) {
			for (unsigned int j = 0; j < i; j++) free(sibling_variables[j]);
			free(sibling_variables);
			return false;
		}
		get_free_variables(*siblings[i], sibling_variables[i]);
		if (sibling_variables[i].length > 1)
			insertion_sort(sibling_variables[i]);
	}

	unsigned int old_arg_variable_count;
	do {
		old_arg_variable_count = keep_variables.length;
		for (unsigned int i = 0; i < siblings.length; i++) {
			if (has_intersection(sibling_variables[i], keep_variables)) {
				array<unsigned int> union_variables(sibling_variables[i].length + keep_variables.length);
				set_union(union_variables, keep_variables, sibling_variables[i]);
				swap(union_variables, keep_variables);
			}
		}
	} while (old_arg_variable_count != keep_variables.length);

	/* remove terms from `siblings` that are unneeded */
	for (unsigned int i = 0; i < siblings.length; i++) {
		if (!has_intersection(sibling_variables[i], keep_variables)) {
			free(sibling_variables[i]);
			move(sibling_variables[siblings.length - 1], sibling_variables[i]);
			siblings.remove(i--);
		} else {
			free(sibling_variables[i]);
		}
	}
	free(sibling_variables);
	return true;
}

template<bool InArray, typename MakeDstFunction>
inline hol_term* try_make_dst(hol_term* head,
		unsigned int& max_variable,
		head_index predicate_index,
		bool& remove_wide_scope_marker,
		array<hol_term*>& siblings,
		MakeDstFunction make_dst)
{
	unsigned int head_variable = 0;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		if (predicate_index.position == head_position::NONE) {
			if (head->any.included != nullptr)
				max_bound_variable(*head->any.included, max_variable);
			head_variable = max_variable + 1;
		} else {
			head_variable = head->any.included->quantifier.variable;
		}
	} else {
		hol_term* inner = head;
		while (inner->type == hol_term_type::NOT)
			inner = inner->unary.operand;
		if (inner->type == hol_term_type::EXISTS)
			head_variable = inner->quantifier.variable;
	}
	return make_dst(head, head_variable, predicate_index, InArray, remove_wide_scope_marker, siblings);
}

template<typename TryFindHeadFunction, typename MakeDstFunction>
inline bool apply_head_to_array(
		hol_array_term& src,
		hol_term**& dst,
		unsigned int& max_variable,
		bool& remove_wide_scope_marker,
		array<hol_term*>& siblings,
		TryFindHeadFunction try_find_head,
		MakeDstFunction make_dst)
{
	dst = nullptr;
	hol_term** heads = (hol_term**) malloc(sizeof(hol_term*) * src.length);
	if (heads == nullptr) {
		fprintf(stderr, "apply_head_to_array ERROR: Insufficient memory for `heads`.\n");
		return false;
	}
	head_index* predicate_indices = (head_index*) calloc(src.length, sizeof(head_index));
	if (predicate_indices == nullptr) {
		fprintf(stderr, "apply_head_to_array ERROR: Insufficient memory for `predicate_indices`.\n");
		free(heads); return false;
	}
	bool success = true;
	for (unsigned int i = 0; i < src.length; i++) {
		try_find_head(src.operands[i], heads[i], predicate_indices[i]);
		if (heads[i] == nullptr) {
			success = false;
			break;
		}
	}

	if (!success) {
		free(heads);
		free(predicate_indices);
		return true;
	}

	dst = (hol_term**) calloc(src.length, sizeof(hol_term*));
	if (dst == nullptr) {
		fprintf(stderr, "apply_head_to_array ERROR: Insufficient memory for `dst`.\n");
		free(heads); free(predicate_indices); return false;
	}
	for (unsigned int i = 0; i < src.length; i++) {
		dst[i] = try_make_dst<true>(heads[i], max_variable, predicate_indices[i], remove_wide_scope_marker, siblings, make_dst);
		if (dst[i] == nullptr) {
			success = false;
			break;
		}
	}

	free(heads);
	free(predicate_indices);
	if (!success) {
		for (unsigned int i = 0; i < src.length; i++) {
			if (dst[i] == nullptr) break;
			free(*dst[i]); if (dst[i]->reference_count == 0) free(dst[i]);
		}
		free(dst);
		dst = nullptr;
		return false;
	}
	return true;
}

template<typename TryFindHeadFunction, typename MakeDstFunction, typename ApplyFunction>
inline bool apply_head(
		hol_term* src, hol_term*& dst,
		array<unsigned int>& dst_variables,
		array<hol_term*>& siblings,
		unsigned int max_variable,
		bool& removed_quantifier,
		bool& remove_negations,
		bool& remove_wide_scope_marker,
		TryFindHeadFunction try_find_head,
		MakeDstFunction make_dst,
		ApplyFunction apply)
{
	if (!apply(src)) return false;

	/* check if the current scope is the head */
	hol_term* head;
	head_index predicate_index = {head_position::NONE, 0};
	try_find_head(src, head, predicate_index);
	if (head != nullptr) {
		dst = try_make_dst<false>(head, max_variable, predicate_index, remove_wide_scope_marker, siblings, make_dst);
		if (dst == nullptr) return false;

		get_free_variables(*dst, dst_variables);
		if (dst_variables.length > 1) insertion_sort(dst_variables);
		if (!prune_independent_siblings(dst_variables, siblings)) {
			free(*dst); if (dst->reference_count == 0) free(dst);
			return false;
		}
		removed_quantifier = false;
		return true;
	}

	if (src->type == hol_term_type::AND || src->type == hol_term_type::OR) {
		hol_term** new_dst;
		if (!apply_head_to_array(src->array, new_dst, max_variable, remove_wide_scope_marker, siblings, try_find_head, make_dst))
			return false;
		if (new_dst != nullptr) {
			for (unsigned int i = 0; i < src->array.length; i++)
				get_free_variables(*new_dst[i], dst_variables);
			if (dst_variables.length > 1) insertion_sort(dst_variables);
			if (!prune_independent_siblings(dst_variables, siblings)) {
				for (unsigned int i = 0; i < src->array.length; i++) {
					free(*new_dst[i]); if (new_dst[i]->reference_count == 0) free(new_dst[i]);
				}
				free(new_dst);
				return false;
			}

			removed_quantifier = false;
			if (src->type == hol_term_type::AND) {
				dst = hol_term::new_and(make_array_view(new_dst, src->array.length));
			} else {
				dst = hol_term::new_or(make_array_view(new_dst, src->array.length));
			}
			if (dst == nullptr) {
				for (unsigned int i = 0; i < src->array.length; i++) {
					free(*new_dst[i]); if (new_dst[i]->reference_count == 0) free(new_dst[i]);
				}
				free(new_dst);
				return false;
			}
			free(new_dst);
			return true;
		}
	} else if (src->type == hol_term_type::ANY_ARRAY && (src->any_array.oper == hol_term_type::ANY_ARRAY || src->any_array.oper == hol_term_type::AND || src->any_array.oper == hol_term_type::OR)) {
		hol_term* all;
		try_find_head(src->any_array.all, all, predicate_index);
		if (all != nullptr) {
			hol_term* new_all = try_make_dst<true>(all, max_variable, predicate_index, remove_wide_scope_marker, siblings, make_dst);
			if (new_all == nullptr) return false;

			hol_term** new_left; hol_term** new_right; hol_term** new_any;
			if (!apply_head_to_array(src->any_array.left, new_left, max_variable, remove_wide_scope_marker, siblings, try_find_head, make_dst)) {
				free(*new_all); if (new_all->reference_count == 0) free(new_all);
				return false;
			}
			if (new_left != nullptr) {
				if (!apply_head_to_array(src->any_array.right, new_right, max_variable, remove_wide_scope_marker, siblings, try_find_head, make_dst)) {
					free(*new_all); if (new_all->reference_count == 0) free(new_all);
					for (unsigned int i = 0; i < src->any_array.left.length; i++) {
						free(*new_left[i]); if (new_left[i]->reference_count == 0) free(new_left[i]);
					}
					free(new_left);
					return false;
				}

				if (new_right != nullptr) {
					if (!apply_head_to_array(src->any_array.any, new_any, max_variable, remove_wide_scope_marker, siblings, try_find_head, make_dst)) {
						free(*new_all); if (new_all->reference_count == 0) free(new_all);
						for (unsigned int i = 0; i < src->any_array.left.length; i++) {
							free(*new_left[i]); if (new_left[i]->reference_count == 0) free(new_left[i]);
						} for (unsigned int i = 0; i < src->any_array.right.length; i++) {
							free(*new_right[i]); if (new_right[i]->reference_count == 0) free(new_right[i]);
						}
						free(new_left); free(new_right);
						return false;
					}

					if (new_any != nullptr) {
						get_free_variables(*new_all, dst_variables);
						for (unsigned int i = 0; i < src->any_array.left.length; i++)
							get_free_variables(*new_left[i], dst_variables);
						for (unsigned int i = 0; i < src->any_array.right.length; i++)
							get_free_variables(*new_right[i], dst_variables);
						for (unsigned int i = 0; i < src->any_array.any.length; i++)
							get_free_variables(*new_any[i], dst_variables);
						if (dst_variables.length > 1) insertion_sort(dst_variables);
						if (!prune_independent_siblings(dst_variables, siblings)) {
							free(*new_all); if (new_all->reference_count == 0) free(new_all);
							for (unsigned int i = 0; i < src->any_array.left.length; i++) {
								free(*new_left[i]); if (new_left[i]->reference_count == 0) free(new_left[i]);
							} for (unsigned int i = 0; i < src->any_array.right.length; i++) {
								free(*new_right[i]); if (new_right[i]->reference_count == 0) free(new_right[i]);
							} for (unsigned int i = 0; i < src->any_array.any.length; i++) {
								free(*new_any[i]); if (new_any[i]->reference_count == 0) free(new_any[i]);
							}
							free(new_left); free(new_right); free(new_any);
							return false;
						}

						removed_quantifier = false;
						dst = hol_term::new_any_array(src->any_array.oper, new_all,
								make_array_view(new_any, src->any_array.any.length),
								make_array_view(new_left, src->any_array.left.length),
								make_array_view(new_right, src->any_array.right.length));
						if (dst == nullptr) {
							free(*new_all); if (new_all->reference_count == 0) free(new_all);
							for (unsigned int i = 0; i < src->any_array.left.length; i++) {
								free(*new_left[i]); if (new_left[i]->reference_count == 0) free(new_left[i]);
							} for (unsigned int i = 0; i < src->any_array.right.length; i++) {
								free(*new_right[i]); if (new_right[i]->reference_count == 0) free(new_right[i]);
							} for (unsigned int i = 0; i < src->any_array.any.length; i++) {
								free(*new_any[i]); if (new_any[i]->reference_count == 0) free(new_any[i]);
							}
							free(new_left); free(new_right); free(new_any);
							return false;
						}
						free(new_left); free(new_right); free(new_any);
						return true;
					}

					for (unsigned int i = 0; i < src->any_array.right.length; i++) {
						free(*new_right[i]); if (new_right[i]->reference_count == 0) free(new_right[i]);
					}
					free(new_right);
				}

				for (unsigned int i = 0; i < src->any_array.left.length; i++) {
					free(*new_left[i]); if (new_left[i]->reference_count == 0) free(new_left[i]);
				}
				free(new_left);
			}

			free(*new_all); if (new_all->reference_count == 0) free(new_all);
		}
	}

	/* we didn't find the head, so keep looking */
	/* NOTE: this function should mirror the semantics of
	   `apply_head`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	hol_term* new_formula;
	switch (src->type) {
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::EQUALS:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IFF:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
	case hol_term_type::ANY_RIGHT_ONLY:
		dst = nullptr;
		removed_quantifier = false;
		remove_negations = false;
		return true;

	case hol_term_type::NOT:
		if (!apply_head(src->unary.operand, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (remove_negations) {
				dst = new_formula;
			} else if (new_formula == src->unary.operand) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				dst = hol_term::new_not(new_formula);
				if (dst == nullptr) {
					free(*new_formula);
					if (new_formula->reference_count == 0)
						free(new_formula);
					return false;
				}
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::AND:
	case hol_term_type::OR:
		if (!siblings.ensure_capacity(siblings.length + src->array.length - 1))
			return false;
		for (unsigned int i = 0; i + 1 < src->array.length; i++)
			siblings[siblings.length++] = src->array.operands[i];
		if (!apply_head(src->array.operands[src->array.length - 1], new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		remove_negations = false;
		if (new_formula != nullptr) {
			hol_term** new_dst = (hol_term**) malloc(sizeof(hol_term*) * src->array.length);
			if (new_dst == nullptr) {
				free(*new_formula);
				if (new_formula->reference_count == 0)
					free(new_formula);
				return false;
			}
			new_dst[src->array.length - 1] = new_formula;

			unsigned int new_start = src->array.length - 1;
			unsigned int new_length = 1;
			for (unsigned int i = src->array.length - 1; i > 0; i--) {
				if (siblings.contains(src->array.operands[i - 1])) {
					new_dst[--new_start] = src->array.operands[i - 1];
					src->array.operands[i - 1]->reference_count++;
					new_length++;
				} else {
					if (new_length == 1) {
						dst = new_dst[new_start];
					} else if (src->type == hol_term_type::AND) {
						dst = hol_term::new_and(make_array_view(new_dst + new_start, new_length));
					} else {
						dst = hol_term::new_or(make_array_view(new_dst + new_start, new_length));
					}
					if (dst == nullptr) {
						for (unsigned int j = 0; j < new_length; j++) {
							free(*new_dst[new_start + j]);
							if (new_dst[new_start + j]->reference_count == 0)
								free(new_dst[new_start + j]);
						}
						free(new_dst); return false;
					}
					new_length = 1;
					new_dst[new_start] = dst;
				}
			}

			if (new_length == 1) {
				dst = new_dst[new_start];
			} else if (new_length == src->array.length && new_formula == src->array.operands[src->array.length - 1]) {
				for (unsigned int i = 0; i < new_length; i++) {
					free(*new_dst[i]);
					if (new_dst[i]->reference_count == 0)
						free(new_dst[i]);
				}
				dst = src;
				dst->reference_count++;
			} else if (src->type == hol_term_type::AND) {
				dst = hol_term::new_and(make_array_view(new_dst + new_start, new_length));
			} else {
				dst = hol_term::new_or(make_array_view(new_dst + new_start, new_length));
			}
			if (dst == nullptr) {
				for (unsigned int i = 0; i < new_length; i++) {
					free(*new_dst[i]);
					if (new_dst[i]->reference_count == 0)
						free(new_dst[i]);
				}
				free(new_dst);
				return false;
			}
			free(new_dst);
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!apply_head(src->quantifier.operand, new_formula, dst_variables, siblings, max(max_variable, src->quantifier.variable), removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (!dst_variables.contains(src->quantifier.variable)) {
				removed_quantifier = true;
				remove_negations = true;
				dst = new_formula;
			} else {
				remove_negations = false;
				if (new_formula == src->quantifier.operand) {
					free(*new_formula);
					dst = src;
					dst->reference_count++;
				} else {
					if (src->type == hol_term_type::FOR_ALL)
						dst = hol_term::new_for_all(src->quantifier.variable_type, src->quantifier.variable, new_formula);
					else if (src->type == hol_term_type::EXISTS)
						dst = hol_term::new_exists(src->quantifier.variable_type, src->quantifier.variable, new_formula);
					else if (src->type == hol_term_type::LAMBDA)
						dst = hol_term::new_lambda(src->quantifier.variable_type, src->quantifier.variable, new_formula);
					if (dst == nullptr) {
						free(*new_formula);
						if (new_formula->reference_count == 0)
							free(new_formula);
						return false;
					}
				}
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::IF_THEN:
	case hol_term_type::UNARY_APPLICATION:
		if (src->type == hol_term_type::IF_THEN) {
			if (!siblings.ensure_capacity(siblings.length + 1))
				return false;
			siblings[siblings.length++] = src->binary.left;
		}
		if (!apply_head(src->binary.right, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		remove_negations = false;
		if (new_formula != nullptr) {
			if (remove_wide_scope_marker && src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
				dst = new_formula;
			} else if (src->type == hol_term_type::IF_THEN && !siblings.contains(src->binary.left)) {
				dst = new_formula;
			} else if (new_formula == src->binary.right) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				if (src->type == hol_term_type::IF_THEN)
					dst = hol_term::new_if_then(src->binary.left, new_formula);
				else dst = hol_term::new_apply(src->binary.left, new_formula);
				if (dst == nullptr) {
					free(*new_formula);
					if (new_formula->reference_count == 0)
						free(new_formula);
					return false;
				}
				src->binary.left->reference_count++;
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		if (src->any.included != nullptr) {
			if (!apply_head(src->any.included, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
				return false;
			remove_negations = false;
			if (new_formula != nullptr) {
				if (new_formula == src->any.included) {
					free(*new_formula);
					dst = src;
					dst->reference_count++;
				} else if (new_formula->type == hol_term_type::ANY && new_formula->any.excluded_tree_count == 0) {
					dst = new_formula;
				} else if (new_formula->type == hol_term_type::ANY_RIGHT && new_formula->any.excluded_tree_count == 0) {
					if (src->type == hol_term_type::ANY_RIGHT) {
						dst = new_formula;
					} else {
						dst = hol_term::new_any(new_formula->any.included);
						if (dst == nullptr) {
							free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
							return false;
						}
					}
				} else if (new_formula->type == hol_term_type::ANY || new_formula->type == hol_term_type::ANY_RIGHT) {
					array<hol_term*> excluded(max(1u, new_formula->any.excluded_tree_count + src->any.excluded_tree_count));
					for (unsigned int i = 0; i < new_formula->any.excluded_tree_count; i++)
						excluded[excluded.length++] = new_formula->any.excluded_trees[i];
					for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
						excluded[excluded.length++] = src->any.excluded_trees[i];
					if (src->type == hol_term_type::ANY || new_formula->type == hol_term_type::ANY)
						dst = hol_term::new_any(new_formula->any.included, excluded.data, excluded.length);
					else dst = hol_term::new_any_right(new_formula->any.included, excluded.data, excluded.length);
					if (dst == nullptr) {
						free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
						return false;
					}
					if (new_formula->any.included != nullptr)
						new_formula->any.included->reference_count++;
					for (hol_term* tree : excluded)
						tree->reference_count++;
					free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
				} else {
					if (src->type == hol_term_type::ANY)
						dst = hol_term::new_any(new_formula, src->any.excluded_trees, src->any.excluded_tree_count);
					else dst = hol_term::new_any_right(new_formula, src->any.excluded_trees, src->any.excluded_tree_count);
					if (dst == nullptr) {
						free(*new_formula); if (new_formula->reference_count == 0) free(new_formula);
						return false;
					}
					for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
						src->any.excluded_trees[i]->reference_count++;
				}
			} else {
				dst = nullptr;
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::ANY_ARRAY:
		if (src->any_array.right.length == 0) {
			dst = nullptr;
			return true;
		} if (!apply_head(src->any_array.right.operands[src->any_array.right.length - 1], new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		remove_negations = false;
		if (new_formula != nullptr) {
			if (new_formula == src->any_array.right.operands[src->any_array.right.length - 1]) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				dst = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), new_formula));
				if (dst == nullptr) {
					free(*new_formula);
					if (new_formula->reference_count == 0)
						free(new_formula);
					return false;
				}
				dst->any_array.all->reference_count++;
				for (unsigned int i = 0; i < dst->any_array.any.length; i++)
					dst->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 0; i < dst->any_array.left.length; i++)
					dst->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i + 1 < dst->any_array.right.length; i++)
					dst->any_array.right.operands[i]->reference_count++;
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::ANY_QUANTIFIER:
		if (!apply_head(src->any_quantifier.operand, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_negations, remove_wide_scope_marker, try_find_head, make_dst, apply))
			return false;
		remove_negations = false;
		if (new_formula != nullptr) {
			if (new_formula == src->any_quantifier.operand) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				dst = hol_term::new_any_quantifier(src->any_quantifier.quantifier, new_formula);
				if (dst == nullptr) {
					free(*new_formula);
					if (new_formula->reference_count == 0)
						free(new_formula);
					return false;
				}
			}
		} else {
			dst = nullptr;
		}
		return true;
	}
	fprintf(stderr, "apply_head ERROR: Unrecognized hol_term_type.\n");
	return false;
}

enum class any_node_position {
	LEFT,
	RIGHT,
	NONE
};

template<any_node_position AnyNodePosition, bool AnyRightOnly>
struct head_substituter {
	bool found_wide_scope;
	bool could_have_wide_scope;
	bool first_any;
	bool has_any;
	bool parent_is_negation;
	const hol_term* src;
	hol_term* dst;
	unsigned int last_declared_variable;

	head_substituter(bool found_wide_scope, bool could_have_wide_scope, bool first_any, const hol_term* src, hol_term* dst) :
			found_wide_scope(found_wide_scope), could_have_wide_scope(could_have_wide_scope), first_any(first_any), has_any(false), parent_is_negation(false), src(src), dst(dst), last_declared_variable(0) { }
};

template<bool AnyRightOnly>
inline hol_term* wrap_any_right(hol_term* operand, unsigned int prev_declared_variable, bool exclude_wide_scope)
{
	if (prev_declared_variable != 0) {
		unsigned int excluded_tree_count = 3 + (exclude_wide_scope ? 1 : 0);
		hol_term* excluded_trees[4];
		excluded_trees[0] = hol_term::new_any(hol_term::new_for_all(prev_declared_variable, &HOL_ANY));
		excluded_trees[1] = hol_term::new_any(hol_term::new_exists(prev_declared_variable, &HOL_ANY));
		excluded_trees[2] = hol_term::new_any(hol_term::new_lambda(prev_declared_variable, &HOL_ANY));
		if (exclude_wide_scope)
			excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
		if (excluded_trees[0] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[1] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[2] != nullptr) HOL_ANY.reference_count++;
		if (exclude_wide_scope && excluded_trees[3] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || (exclude_wide_scope && excluded_trees[3] == nullptr)) {
			if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
			if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
			if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
			return nullptr;
		}

		hol_term* term;
		if (AnyRightOnly)
			term = hol_term::new_any_right_only(operand, excluded_trees, excluded_tree_count);
		else term = hol_term::new_any_right(operand, excluded_trees, excluded_tree_count);
		if (term == nullptr) {
			free(*excluded_trees[0]); free(excluded_trees[0]);
			free(*excluded_trees[1]); free(excluded_trees[1]);
			free(*excluded_trees[2]); free(excluded_trees[2]);
			free(*excluded_trees[3]); free(excluded_trees[3]);
			return nullptr;
		}
		return term;
	} else if (exclude_wide_scope) {
		hol_term* excluded = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
		if (excluded == nullptr) return nullptr;
		HOL_ANY.reference_count++;
		hol_term* term;
		if (AnyRightOnly)
			term = hol_term::new_any_right_only(operand, &excluded, 1);
		else term = hol_term::new_any_right(operand, &excluded, 1);
		if (term == nullptr) {
			free(*excluded); free(excluded);
			return nullptr;
		}
		return term;
	} else {
		if (AnyRightOnly)
			return hol_term::new_any_right_only(operand);
		else return hol_term::new_any_right(operand);
	}
}

template<hol_term_type Type, any_node_position AnyNodePosition, bool AnyRightOnly>
inline hol_term* apply(hol_term* src, head_substituter<AnyNodePosition, AnyRightOnly>& substituter)
{
	if (src == substituter.src) {
		if (AnyNodePosition == any_node_position::LEFT) {
			if (substituter.dst->type == hol_term_type::ANY || substituter.dst->type == hol_term_type::ANY_RIGHT
			 || (substituter.dst->type == hol_term_type::UNARY_APPLICATION && substituter.dst->binary.left->type == hol_term_type::CONSTANT && substituter.dst->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE))
				substituter.could_have_wide_scope = true;
			if (substituter.has_any || substituter.dst->type == hol_term_type::ANY || substituter.dst->type == hol_term_type::ANY_RIGHT) {
				if (substituter.src != substituter.dst)
					substituter.dst->reference_count++;
				return substituter.dst;
			} else {
				hol_term* new_term = wrap_any_right<AnyRightOnly>(substituter.dst, substituter.last_declared_variable, substituter.found_wide_scope || (!substituter.could_have_wide_scope && substituter.first_any));
				substituter.first_any = false;
				if (new_term == nullptr)
					return nullptr;
				substituter.dst->reference_count++;
				return new_term;
			}
		} else if (AnyNodePosition == any_node_position::RIGHT && (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT)) {
			hol_term* new_term = wrap_any_right<AnyRightOnly>(substituter.dst, substituter.last_declared_variable, substituter.found_wide_scope || (!substituter.could_have_wide_scope && substituter.first_any));
			substituter.first_any = false;
			if (new_term == nullptr)
				return nullptr;
			substituter.dst->reference_count++;
			return new_term;
		}
		if (substituter.src != substituter.dst)
			substituter.dst->reference_count++;
		return substituter.dst;
	}

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	bool changed; hol_term* new_term = nullptr;
	hol_term* first; hol_term* second; hol_term* third;
	unsigned int prev_declared_variable = 0;
	bool found_wide_scope = false, first_any = false, has_any = false, parent_is_negation = false;
	switch (Type) {
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		if ((AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && src->binary.right->type != hol_term_type::ANY && src->binary.right->type != hol_term_type::ANY_RIGHT)
		 || (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::UNARY_APPLICATION))
		{
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.found_wide_scope = (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
		substituter.parent_is_negation = false;
		first = src->binary.left;
		second = apply(src->binary.right, substituter);
		if (second == nullptr)
			return nullptr;

		if (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
			substituter.could_have_wide_scope = true;

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && src->binary.right->type != hol_term_type::ANY && src->binary.right->type != hol_term_type::ANY_RIGHT && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(second, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (second != src->binary.right) { free(*second); if (second->reference_count == 0) free(second); }
				return nullptr;
			}
			if (second == src->binary.right) second->reference_count++;
			second = term;
		}

		if (second == src->binary.right) {
			new_term = src;
		} else {
			if (!new_hol_term(new_term)) {
				if (second != src->binary.right) {
					free(*second); if (second->reference_count == 0) free(second);
				}
				return nullptr;
			}
			new_term->binary.left = first;
			new_term->binary.right = second;
			first->reference_count++;
			if (second == src->binary.right) second->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
		}

		if (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::UNARY_APPLICATION && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::BINARY_APPLICATION:
		if (AnyNodePosition == any_node_position::RIGHT || AnyNodePosition == any_node_position::LEFT) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
			substituter.parent_is_negation = false;
		}
		first = src->ternary.first;
		second = src->ternary.second;
		third = apply(src->ternary.third, substituter);
		if (third == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT) {
			hol_term* term = wrap_any_right<AnyRightOnly>(third, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (third != src->unary.operand) { free(*third); if (third->reference_count == 0) free(third); }
				return nullptr;
			}
			if (third == src->unary.operand) third->reference_count++;
			third = term;
		}

		if (third == src->ternary.third) {
			new_term = src;
		} else {
			if (!new_hol_term(new_term)) {
				if (third != src->ternary.third) {
					free(*third); if (third->reference_count == 0) free(third);
				}
				return nullptr;
			}
			new_term->ternary.first = first;
			new_term->ternary.second = second;
			new_term->ternary.third = third;
			first->reference_count++;
			second->reference_count++;
			if (third == src->ternary.third) third->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
		}

		if (AnyNodePosition == any_node_position::LEFT && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		changed = false;
		substituter.has_any = false;
		substituter.parent_is_negation = false;
		first = apply(src->array.operands[src->array.length - 1], substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->array.operands[src->array.length - 1])
			changed = true;

		if (!changed) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->array.operands[src->array.length - 1]) {
					free(*first);
					if (first->reference_count == 0)
						free(first);
				}
				return nullptr;
			}
			new_term->type = Type;
			new_term->reference_count = 1;
			if (first->type == Type && first->reference_count == 1) {
				/* the new child is the same type as the parent, so merge them */
				new_term->array.operands = (hol_term**) realloc(first->array.operands, sizeof(hol_term*) * (src->array.length + first->array.length - 1));
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first); if (first->reference_count == 0) free(first);
					}
					return nullptr;
				}
				for (unsigned int i = first->array.length; i > 0; i--)
					new_term->array.operands[i + src->array.length - 2] = new_term->array.operands[i - 1];
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.length = src->array.length + first->array.length - 1;
				free(first);
			} else {
				new_term->array.length = src->array.length;
				new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * new_term->array.length);
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					free(new_term);
					return nullptr;
				}
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.operands[src->array.length - 1] = first;
				if (first == src->array.operands[src->array.length - 1])
					first->reference_count++;
			}
			return new_term;
		}

	case hol_term_type::ANY_ARRAY:
		changed = false;
		substituter.has_any = false;
		substituter.parent_is_negation = false;
		if (src->any_array.right.length == 0)
			first = apply(src->any_array.all, substituter);
		else first = apply(src->any_array.right.operands[src->any_array.right.length - 1], substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != (src->any_array.right.length == 0 ? src->any_array.all : src->any_array.right.operands[src->any_array.right.length - 1])) {
			changed = true;
		}

		if (!changed) {
			new_term = src;
		} else {
			if (src->any_array.right.length == 0) {
				new_term = hol_term::new_any_array(src->any_array.oper, first,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_array_view(src->any_array.right.operands, src->any_array.right.length));
			} else {
				new_term = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), first));
			}
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
			new_term->any_array.all->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.any.length; i++)
				new_term->any_array.any.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.left.length; i++)
				new_term->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.right.length; i++)
				new_term->any_array.right.operands[i]->reference_count++;
			first->reference_count--;
		}
		return new_term;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if ((AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA) || AnyNodePosition == any_node_position::LEFT) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		if (!(AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA))
			substituter.last_declared_variable = src->quantifier.variable;
		substituter.parent_is_negation = false;
		first = apply(src->quantifier.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA) {
			hol_term* term = wrap_any_right<AnyRightOnly>(first, src->quantifier.variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (first != src->quantifier.operand) { free(*first); if (first->reference_count == 0) free(first); }
				return nullptr;
			}
			if (first == src->quantifier.operand) first->reference_count++;
			first = term;
		}

		if (first == src->quantifier.operand) {
			new_term = src;
		} else {
			if (src->type == hol_term_type::FOR_ALL)
				new_term = hol_term::new_for_all(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::EXISTS)
				new_term = hol_term::new_exists(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::LAMBDA)
				new_term = hol_term::new_lambda(src->quantifier.variable_type, src->quantifier.variable, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && !parent_is_negation && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::ANY_QUANTIFIER:
		if (AnyNodePosition == any_node_position::LEFT && src->any_quantifier.quantifier != hol_quantifier_type::EXISTS) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.parent_is_negation = false;
		first = apply(src->any_quantifier.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (first == src->any_quantifier.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_any_quantifier(src->any_quantifier.quantifier, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && src->any_quantifier.quantifier != hol_quantifier_type::EXISTS && !parent_is_negation && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::NOT:
		if (AnyNodePosition == any_node_position::RIGHT || AnyNodePosition == any_node_position::LEFT) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.parent_is_negation = true;
		first = apply(src->unary.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT) {
			hol_term* term = wrap_any_right<AnyRightOnly>(first, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (first != src->unary.operand) { free(*first); if (first->reference_count == 0) free(first); }
				return nullptr;
			}
			if (first == src->unary.operand) first->reference_count++;
			first = term;
		}

		if (first == src->unary.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_not(first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope));
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_RIGHT_ONLY:
		if (src->any.included == nullptr)
			return nullptr;
		substituter.has_any = true;
		substituter.parent_is_negation = false;
		first = apply(src->any.included, substituter);
		if (first == nullptr)
			return nullptr;
		substituter.could_have_wide_scope = true;

		if (first == src->any.included) {
			new_term = src;
		} else {
			if (AnyNodePosition == any_node_position::LEFT && first->type == src->type) {
				array<hol_term*> excluded(max(1u, first->any.excluded_tree_count + src->any.excluded_tree_count));
				for (unsigned int i = 0; i < first->any.excluded_tree_count; i++)
					excluded[excluded.length++] = first->any.excluded_trees[i];
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
					excluded[excluded.length++] = src->any.excluded_trees[i];
				if (src->type == hol_term_type::ANY)
					new_term = hol_term::new_any(first->any.included, excluded.data, excluded.length);
				else new_term = hol_term::new_any_right(first->any.included, excluded.data, excluded.length);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				if (first->any.included != nullptr)
					first->any.included->reference_count++;
				for (hol_term* tree : excluded)
					tree->reference_count++;
				free(*first); if (first->reference_count == 0) free(first);
			} else if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
				if (src->type == hol_term_type::ANY)
					new_term = hol_term::new_any(first, src->any.excluded_trees, src->any.excluded_tree_count);
				else new_term = hol_term::new_any_right(first, src->any.excluded_trees, src->any.excluded_tree_count);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
					src->any.excluded_trees[i]->reference_count++;
			} else {
				new_term = hol_term::new_any_right_only(first, src->any.excluded_trees, src->any.excluded_tree_count);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
					src->any.excluded_trees[i]->reference_count++;
			}
		}
		return new_term;

	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
		return default_apply<Type>(src, substituter);
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type when substituting head.\n");
	return NULL;
}

template<any_node_position AnyNodePosition, bool AnyRightOnly = true>
inline hol_term* substitute_head(hol_term* src,
		const hol_term* src_term, hol_term* dst_term,
		bool could_have_wide_scope = false)
{
	head_substituter<AnyNodePosition, AnyRightOnly> substituter(false, could_have_wide_scope, true, src_term, dst_term);
	hol_term* dst = apply(src, substituter);
	if (dst == nullptr)
		return nullptr;

	if (AnyNodePosition == any_node_position::RIGHT && dst->type == hol_term_type::ANY_RIGHT_ONLY) {
		hol_term* term = hol_term::new_any_right(dst->any.included, dst->any.excluded_trees, dst->any.excluded_tree_count);
		if (term->any.included != nullptr)
			term->any.included->reference_count++;
		for (unsigned int i = 0; i < term->any.excluded_tree_count; i++)
			term->any.excluded_trees[i]->reference_count++;
		if (dst != src) { free(*dst); if (dst->reference_count == 0) free(dst); }
		if (term == nullptr)
			return nullptr;
		dst = term;
	} else if (AnyNodePosition == any_node_position::RIGHT && dst->type != hol_term_type::LAMBDA && dst->type != hol_term_type::ANY && dst->type != hol_term_type::ANY_RIGHT) {
		hol_term* term = hol_term::new_any_right(dst);
		if (term == nullptr) {
			if (dst != src) { free(*dst); if (dst->reference_count == 0) free(dst); }
			return nullptr;
		}
		if (dst == src) dst->reference_count++;
		dst = term;
	}

	if (dst == src)
		dst->reference_count++;
	return dst;
}

template<typename TryFindHeadFunction>
struct any_node_remover {
	TryFindHeadFunction try_find_head;
};

template<hol_term_type Type, typename TryFindHeadFunction>
inline hol_term* apply(hol_term* src, any_node_remover<TryFindHeadFunction>& remover)
{
	if (src->type != hol_term_type::ANY_RIGHT_ONLY) {
		hol_term* head; head_index predicate_index;
		remover.try_find_head(src, head, predicate_index);
		if (head != nullptr)
			return src;
	}

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	bool changed; hol_term* new_term = nullptr;
	hol_term* first; hol_term* second; hol_term* third;
	switch (Type) {
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		first = src->binary.left;
		second = apply(src->binary.right, remover);
		if (second == nullptr)
			return nullptr;

		if (second == src->binary.right) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (second != src->binary.right) {
					free(*second); if (second->reference_count == 0) free(second);
				}
				return nullptr;
			}
			new_term->binary.left = first;
			new_term->binary.right = second;
			first->reference_count++;
			if (second == src->binary.right) second->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
			return new_term;
		}

	case hol_term_type::BINARY_APPLICATION:
		first = src->ternary.first;
		second = src->ternary.second;
		third = apply(src->ternary.third, remover);
		if (third == nullptr)
			return nullptr;

		if (third == src->ternary.third) {
			new_term = src;
		} else {
			if (!new_hol_term(new_term)) {
				if (third != src->ternary.third) {
					free(*third); if (third->reference_count == 0) free(third);
				}
				return nullptr;
			}
			new_term->ternary.first = first;
			new_term->ternary.second = second;
			new_term->ternary.third = third;
			first->reference_count++;
			second->reference_count++;
			if (third == src->ternary.third) third->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
		}
		return new_term;

	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		changed = false;
		first = apply(src->array.operands[src->array.length - 1], remover);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->array.operands[src->array.length - 1])
			changed = true;

		if (!changed) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->array.operands[src->array.length - 1]) {
					free(*first);
					if (first->reference_count == 0)
						free(first);
				}
				return nullptr;
			}
			new_term->type = Type;
			new_term->reference_count = 1;
			if (first->type == Type && first->reference_count == 1) {
				/* the new child is the same type as the parent, so merge them */
				new_term->array.operands = (hol_term**) realloc(first->array.operands, sizeof(hol_term*) * (src->array.length + first->array.length - 1));
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					return nullptr;
				}
				for (unsigned int i = first->array.length; i > 0; i--)
					new_term->array.operands[i + src->array.length - 2] = new_term->array.operands[i - 1];
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.length = src->array.length + first->array.length - 1;
				free(first);
			} else {
				new_term->array.length = src->array.length;
				new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * new_term->array.length);
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					free(new_term);
					return nullptr;
				}
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.operands[src->array.length - 1] = first;
				if (first == src->array.operands[src->array.length - 1])
					first->reference_count++;
			}
			return new_term;
		}

	case hol_term_type::ANY_ARRAY:
		changed = false;
		if (src->any_array.right.length == 0)
			first = apply(src->any_array.all, remover);
		else first = apply(src->any_array.right.operands[src->any_array.right.length - 1], remover);
		if (first == nullptr) {
			return nullptr;
		} else if (first != (src->any_array.right.length == 0 ? src->any_array.all : src->any_array.right.operands[src->any_array.right.length - 1])) {
			changed = true;
		}

		if (!changed) {
			new_term = src;
		} else {
			if (src->any_array.right.length == 0) {
				new_term = hol_term::new_any_array(src->any_array.oper, first,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_array_view(src->any_array.right.operands, src->any_array.right.length));
			} else {
				new_term = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), first));
			}
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
			new_term->any_array.all->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.any.length; i++)
				new_term->any_array.any.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.left.length; i++)
				new_term->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i + 1 < new_term->any_array.right.length; i++)
				new_term->any_array.right.operands[i]->reference_count++;
		}
		return new_term;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		first = apply(src->quantifier.operand, remover);
		if (first == nullptr)
			return nullptr;

		if (first == src->quantifier.operand) {
			new_term = src;
		} else {
			if (src->type == hol_term_type::FOR_ALL)
				new_term = hol_term::new_for_all(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::EXISTS)
				new_term = hol_term::new_exists(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::LAMBDA)
				new_term = hol_term::new_lambda(src->quantifier.variable_type, src->quantifier.variable, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}
		return new_term;

	case hol_term_type::ANY_QUANTIFIER:
		first = apply(src->any_quantifier.operand, remover);
		if (first == nullptr)
			return nullptr;

		if (first == src->any_quantifier.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_any_quantifier(src->any_quantifier.quantifier, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}
		return new_term;

	case hol_term_type::NOT:
		first = apply(src->unary.operand, remover);
		if (first == nullptr)
			return nullptr;

		if (first == src->unary.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_not(first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}
		return new_term;

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
		if (src->any.included == nullptr)
			return nullptr;

		first = apply(src->any.included, remover);
		if (first == nullptr)
			return nullptr;

		if (first == src->any.included) {
			new_term = src;
		} else {
			if (src->type == hol_term_type::ANY)
				new_term = hol_term::new_any(first, src->any.excluded_trees, src->any.excluded_tree_count);
			else new_term = hol_term::new_any_right(first, src->any.excluded_trees, src->any.excluded_tree_count);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
			for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
				src->any.excluded_trees[i]->reference_count++;
		}
		return new_term;

	case hol_term_type::ANY_RIGHT_ONLY:
		first = apply(src->any.included, remover);
		if (first == src->any.included)
			first->reference_count++;
		return first;

	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
		return default_apply<Type>(src, remover);
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type when substituting head.\n");
	return NULL;
}

template<typename TryFindHeadFunction>
inline hol_term* remove_any_nodes(hol_term* src, TryFindHeadFunction try_find_head)
{
	if (src->type == hol_term_type::ANY_RIGHT_ONLY)
		src = src->any.included;

	any_node_remover<TryFindHeadFunction> remover = {try_find_head};
	hol_term* dst = apply(src, remover);

	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct apply_head_inverter {
	array<hol_term*> outer;
	bool has_bound_variables;
	unsigned int max_variable;

	apply_head_inverter() : outer(8), has_bound_variables(false), max_variable(0) { }

	inline void operator() (hol_term* node) {
		/* NOTE: this function should mirror the semantics of
		  `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
		  `find_head` in `hdp_parser.h` */
		if (node->type == hol_term_type::IF_THEN || node->type == hol_term_type::EQUALS || node->type == hol_term_type::UNARY_APPLICATION) {
			has_bound_variables |= max_bound_variable(*node->binary.left, max_variable);
		} else if (node->type == hol_term_type::AND || node->type == hol_term_type::OR) {
			for (unsigned int i = 0; i + 1 < node->array.length; i++)
				has_bound_variables |= max_bound_variable(*node->array.operands[i], max_variable);
		} else if (node->type == hol_term_type::FOR_ALL || node->type == hol_term_type::EXISTS || node->type == hol_term_type::LAMBDA) {
			max_variable = max(max_variable, node->quantifier.variable);
		}
		outer.add(node);
	}
};

template<typename TryFindHeadFunction, typename Function>
hol_term* find_head(hol_term* term, head_index& predicate_index, TryFindHeadFunction& try_find_head, Function& apply)
{
	apply(term);

	hol_term* head;
	try_find_head(term, head, predicate_index);
	if (head != nullptr) return head;

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`,
	   `apply_arg`, and `find_head` in `hdp_parser.h` */
	switch (term->type) {
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IFF:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
		return nullptr;

	case hol_term_type::NOT:
		return find_head(term->unary.operand, predicate_index, try_find_head, apply);

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return find_head(term->quantifier.operand, predicate_index, try_find_head, apply);

	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return find_head(term->binary.right, predicate_index, try_find_head, apply);

	case hol_term_type::AND:
	case hol_term_type::OR:
		return find_head(term->array.operands[term->array.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_RIGHT_ONLY:
		if (term->any.included == nullptr) return nullptr;
		return find_head(term->any.included, predicate_index, try_find_head, apply);

	case hol_term_type::ANY_ARRAY:
		if (term->any_array.right.length == 0)
			return find_head(term->any_array.all, predicate_index, try_find_head, apply);
		return find_head(term->any_array.right.operands[term->any_array.right.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY_QUANTIFIER:
		return find_head(term->any_quantifier.operand, predicate_index, try_find_head, apply);
	}
	fprintf(stderr, "find_head ERROR: Unrecognied hol_term_type.\n");
	return nullptr;
}

struct no_op {
	template<typename... Args>
	inline constexpr bool operator() (Args&&... args) const { return true; }
};

template<int_fast8_t ConjunctIndex, bool SelectNegation = false>
inline bool select_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					hol_term* head_var = hol_term::new_variable(head_variable);
					if (head_var == nullptr) return (hol_term*) nullptr;
					constexpr unsigned int additional_excluded_tree_count = 4;
					unsigned int excluded_tree_count = hol_non_head_constants<built_in_predicates>::count() + additional_excluded_tree_count;
					hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term*) * excluded_tree_count);
					excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
					excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
					excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
					excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
					if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
						if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
						if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
						if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
						free(*head_var); free(head_var);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++)
						excluded_trees[i + additional_excluded_tree_count] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
					free(*head_var);

					hol_term* predicate_term = hol_term::new_apply(
							hol_term::new_any(nullptr, excluded_trees, hol_non_head_constants<built_in_predicates>::count() + additional_excluded_tree_count),
							hol_term::new_variable(head_variable));
					if (predicate_term == nullptr) {
						free(*excluded_trees[0]); free(excluded_trees[0]);
						free(*excluded_trees[1]); free(excluded_trees[1]);
						free(*excluded_trees[2]); free(excluded_trees[2]);
						free(*excluded_trees[3]); free(excluded_trees[3]);
						return (hol_term*) nullptr;
					}
					hol_non_head_constants<built_in_predicates>::increment_terms();

					excluded_trees[additional_excluded_tree_count] = predicate_term;
					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(
						predicate_term,
						hol_term::new_any(nullptr, excluded_trees, additional_excluded_tree_count + 1)
					));
					if (dst == nullptr) {
						free(*predicate_term); free(predicate_term);
						return (hol_term*) nullptr;
					}
					predicate_term->reference_count += 2 - 1;
					for (unsigned int i = 0; i < 4; i++)
						excluded_trees[i]->reference_count++;

					if (negated && SelectNegation) {
						hol_term* new_dst = hol_term::new_not(dst);
						if (new_dst->reference_count == 0) {
							free(*dst); if (dst->reference_count == 0) free(dst);
							return (hol_term*) nullptr;
						}
						dst = new_dst;
					}

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "select_conjunct WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count == 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					if (is_array) {
						dst = hol_term::new_any_right(intersection[0]);
					} else if (SelectNegation && (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT)) {
						/* in case the head is negated */
						dst = hol_term::new_any_right(hol_term::new_any_array(hol_term_type::ANY_ARRAY, hol_term::new_any_right(intersection[0]),
								make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
					} else {
						dst = hol_term::new_any_right(hol_term::new_any_array(hol_term_type::ANY_ARRAY, intersection[0],
								make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
					}
					if (dst == nullptr) {
						free(*intersection[0]); if (intersection[0]->reference_count == 0) free(intersection[0]);
						return (hol_term*) nullptr;
					}
					return dst;
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
#if !defined(NDEBUG)
					if (predicate_index.compare(ConjunctIndex)) {
						fprintf(stderr, "select_conjunct WARNING: The index of the conjunct and the index of the predicate are the same.\n");
						return (hol_term*) nullptr;
					}
#endif
					hol_term* operand = head->quantifier.operand;
					hol_term* conjunct;
					if (ConjunctIndex >= 0) {
						if (ConjunctIndex < operand->any_array.left.length)
							conjunct = operand->any_array.left.operands[ConjunctIndex];
						else return (hol_term*) nullptr;
					} else if (ConjunctIndex < 0) {
						unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						if (index < operand->any_array.right.length)
							conjunct = operand->any_array.right.operands[operand->any_array.right.length - index - 1];
						else return (hol_term*) nullptr;
					}

					hol_term* predicate;
					if (predicate_index.position == head_position::LEFT) {
						predicate = operand->any_array.left.operands[predicate_index.index];
					} else if (predicate_index.position == head_position::RIGHT) {
						predicate = operand->any_array.right.operands[operand->any_array.right.length - predicate_index.index - 1];
					} else {
						predicate = operand->any_array.any.operands[predicate_index.index];
					}
					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(predicate, conjunct));
					if (dst != nullptr) {
						predicate->reference_count++;
						conjunct->reference_count++;
					}
					if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
						hol_term* new_dst = hol_term::new_any_right(dst);
						if (new_dst == nullptr) {
							free(*dst); free(dst);
							return (hol_term*) nullptr;
						}
						dst = new_dst;
					}

					if (negated && SelectNegation) {
						hol_term* new_dst = hol_term::new_not(dst);
						if (new_dst->reference_count == 0) {
							free(*dst); if (dst->reference_count == 0) free(dst);
							return (hol_term*) nullptr;
						}
						dst = new_dst;
					}
					return dst;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
						fprintf(stderr, "select_conjunct WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					unsigned int conjunct_index = (ConjunctIndex < 0) ? (operand->array.length + ConjunctIndex) : ConjunctIndex;
#if !defined(NDEBUG)
					if (conjunct_index == predicate_index.index) {
						fprintf(stderr, "select_conjunct WARNING: The index of the conjunct and the index of the predicate are the same.\n");
						return (hol_term*) nullptr;
					}
#endif
					hol_term* conjunct = operand->array.operands[conjunct_index];
					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(operand->array.operands[predicate_index.index], conjunct));
					if (dst == nullptr)
						return (hol_term*) nullptr;
					operand->array.operands[predicate_index.index]->reference_count++;
					conjunct->reference_count++;
					if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
						hol_term* new_dst = hol_term::new_any_right(dst);
						if (new_dst == nullptr) {
							free(*dst); free(dst);
							return (hol_term*) nullptr;
						}
						dst = new_dst;
					}

					if (negated && SelectNegation) {
						hol_term* new_dst = hol_term::new_not(dst);
						if (new_dst->reference_count == 0) {
							free(*dst); if (dst->reference_count == 0) free(dst);
							return (hol_term*) nullptr;
						}
						dst = new_dst;
					}
					return dst;
				}
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex, bool RemoveNegation = false>
inline bool remove_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					hol_term* head_var = hol_term::new_variable(head_variable);
					if (head_var == nullptr) return (hol_term*) nullptr;
					constexpr unsigned int excluded_tree_count = 4;
					hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term) * (excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()));
					excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
					excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
					excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
					excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
					if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
						if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
						if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
						if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
						free(*head_var); free(head_var);
						return (hol_term*) nullptr;
					}
					free(*head_var);

					for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++) {
						excluded_trees[excluded_tree_count + i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
						excluded_trees[excluded_tree_count + i]->reference_count++;
					}

					hol_term* expected_predicate = hol_term::new_apply(
								hol_term::new_any(nullptr, excluded_trees, excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()),
								hol_term::new_variable(head_variable));
					if (expected_predicate == nullptr) {
						for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
							free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
						}
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++)
						excluded_trees[i]->reference_count++;

					dst = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND,
							hol_term::new_any(nullptr, excluded_trees, excluded_tree_count), make_array_view(&expected_predicate, 1),
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
					if (dst == nullptr) {
						free(*expected_predicate); free(expected_predicate);
						for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
							free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
						}
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < excluded_tree_count; i++)
						excluded_trees[i]->reference_count++;
					for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
						free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
					}

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "remove_conjunct WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count == 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunction = nullptr;
					if (ConjunctIndex >= 0) {
						if (ConjunctIndex < operand->any_array.left.length) {
							conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
									make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
									make_excluded_array_view(operand->any_array.left.operands, operand->any_array.left.length, ConjunctIndex),
									make_array_view(operand->any_array.right.operands, operand->any_array.right.length));
						} else {
							conjunction = operand;
						}
					} else if (ConjunctIndex < 0) {
						if (-ConjunctIndex - 1 < operand->any_array.right.length) {
							unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
							index = operand->any_array.right.length - index - 1;
							conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
									make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
									make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
									make_excluded_array_view(operand->any_array.right.operands, operand->any_array.right.length, index));
						} else {
							conjunction = operand;
						}
					}
					if (conjunction == nullptr)
						return (hol_term*) nullptr;
					if (conjunction == operand) {
						operand->reference_count++;
					} else {
						conjunction->any_array.all->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
							conjunction->any_array.any.operands[i]->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
							conjunction->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
							conjunction->any_array.right.operands[i]->reference_count++;
					}

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
						return (hol_term*) nullptr;
					}
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::AND) {
					hol_term* operand = head->quantifier.operand;
					unsigned int conjunct_index = (ConjunctIndex < 0) ? (operand->array.length + ConjunctIndex) : ConjunctIndex;
					if (operand->array.length == 2) {
						dst = hol_term::new_exists(head_variable, (conjunct_index == 0 ? operand->array.operands[1] : operand->array.operands[0]));
					} else {
						dst = hol_term::new_exists(head_variable, hol_term::new_and(make_excluded_array_view(operand->array.operands, operand->array.length, conjunct_index)));
					}
					if (dst != nullptr) {
						for (unsigned int i = 0; i < operand->array.length; i++)
							if (i != conjunct_index) operand->array.operands[i]->reference_count++;
					}
				} else {
					return (hol_term*) nullptr;
				}

				if (negated && !RemoveNegation) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst, old_head->any.excluded_trees, old_head->any.excluded_tree_count);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < old_head->any.excluded_tree_count; i++)
						old_head->any.excluded_trees[i]->reference_count++;
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex, bool IsDstRoot = false>
inline bool select_set_conjunct(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
		src = src->quantifier.operand;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				if (head->type == hol_term_type::NOT)
					head = head->unary.operand;

				hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
				if (excluded_quantifier == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
				if (expected_conjunct == nullptr) {
					free(*excluded_quantifier); free(excluded_quantifier);
					return (hol_term*) nullptr;
				}

				hol_term* expected_head;
				if (ConjunctIndex >= 0) {
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
							make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(expected_conjunct, ConjunctIndex + 1), make_array_view((hol_term**) nullptr, 0)));
					if (expected_head == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return (hol_term*) nullptr;
					}
					expected_conjunct->reference_count += ConjunctIndex + 1;
				} else {
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(expected_conjunct, (unsigned int) (-ConjunctIndex))));
					if (expected_head == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return (hol_term*) nullptr;
					}
					expected_conjunct->reference_count += (unsigned int) (-ConjunctIndex);
				}

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "select_set_conjunct ERROR: Expected intersection to be unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				hol_term* conjunct;
				hol_term* operand = intersection[0]->quantifier.operand;
				if (operand->type == hol_term_type::ANY_ARRAY) {
					if (ConjunctIndex >= 0) {
						conjunct = operand->any_array.left.operands[ConjunctIndex];
					} else {
						unsigned int index = (unsigned int) (operand->any_array.right.length + ConjunctIndex);
						conjunct = operand->any_array.right.operands[index];
					}
				} else {
					if (ConjunctIndex >= 0) {
						conjunct = operand->array.operands[ConjunctIndex];
					} else {
						unsigned int index = (unsigned int) (operand->array.length + ConjunctIndex);
						conjunct = operand->array.operands[index];
					}
				}

				hol_term* dst = hol_term::new_exists(head_variable, conjunct);
				if (dst == nullptr) {
					free_all(intersection);
					return (hol_term*) nullptr;
				}
				conjunct->reference_count++;
				free_all(intersection);

				if (IsDstRoot) {
					/* check that there are no undeclared variables */
					array<unsigned int> free_variables(4);
					if (!get_free_variables(*dst, free_variables) || free_variables.length != 0) {
						free(*dst); free(dst);
						return (hol_term*) nullptr;
					}
				}

				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool remove_set_conjunct(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				if (head->type == hol_term_type::NOT)
					head = head->unary.operand;

				hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
				if (excluded_quantifier == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
				if (expected_conjunct == nullptr) {
					free(*excluded_quantifier); free(excluded_quantifier);
					return (hol_term*) nullptr;
				}

				hol_term* expected_head;
				if (ConjunctIndex >= 0) {
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
							make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(expected_conjunct, ConjunctIndex + 1), make_array_view((hol_term**) nullptr, 0)));
					if (expected_head == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return (hol_term*) nullptr;
					}
					expected_conjunct->reference_count += ConjunctIndex + 1;
				} else {
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(expected_conjunct, (unsigned int) (-ConjunctIndex))));
					if (expected_head == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return (hol_term*) nullptr;
					}
					expected_conjunct->reference_count += (unsigned int) (-ConjunctIndex);
				}

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "remove_set_conjunct ERROR: Expected intersection to be unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				hol_term* new_conjunction;
				hol_term* operand = intersection[0]->quantifier.operand;
				if (operand->type == hol_term_type::ANY_ARRAY) {
					if (ConjunctIndex >= 0) {
						new_conjunction = hol_term::new_any_array(hol_term_type::AND,
								operand->any_array.all, make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
								make_excluded_array_view(operand->any_array.left.operands, operand->any_array.left.length, ConjunctIndex),
								make_array_view(operand->any_array.right.operands, operand->any_array.right.length));
					} else {
						unsigned int index = (unsigned int) (operand->any_array.right.length + ConjunctIndex);
						new_conjunction = hol_term::new_any_array(hol_term_type::AND,
								operand->any_array.all, make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
								make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
								make_excluded_array_view(operand->any_array.right.operands, operand->any_array.right.length, index));
					}
					if (new_conjunction == nullptr) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					new_conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.any.length; i++)
						new_conjunction->any_array.any.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.left.length; i++)
						new_conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.right.length; i++)
						new_conjunction->any_array.right.operands[i]->reference_count++;
				} else {
					if (ConjunctIndex >= 0) {
						new_conjunction = hol_term::new_and(make_excluded_array_view(operand->array.operands, operand->array.length, ConjunctIndex));
					} else {
						unsigned int index = (unsigned int) (operand->array.length + ConjunctIndex);
						new_conjunction = hol_term::new_and(make_excluded_array_view(operand->array.operands, operand->array.length, index));
					}
					if (new_conjunction == nullptr) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < new_conjunction->array.length; i++)
						new_conjunction->array.operands[i]->reference_count++;
				}
				free_all(intersection);

				hol_term* dst = hol_term::new_exists(head_variable, new_conjunction);
				if (dst == nullptr) {
					free(*new_conjunction); free(new_conjunction);
					return (hol_term*) nullptr;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool select_conjunct_in_set(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
		src = src->quantifier.operand;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				unsigned int set_variable, element_variable = 0;
				hol_term* left = nullptr;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT)
				{
					set_variable = ++max_variable;
					element_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS)
						fprintf(stderr, "select_conjunct_in_set WARNING: Expected an existential quantification.\n");
#endif
					set_variable = head->quantifier.variable;
					hol_term* operand = head->quantifier.operand;

					hol_term* right = nullptr;
					if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
						element_variable = ++max_variable;
					} else if (operand->type == hol_term_type::AND && operand->array.length != 0) {
						left = operand->array.operands[0];
						right = operand->array.operands[operand->array.length - 1];
					} else if (operand->type == hol_term_type::ANY_ARRAY) {
						if (operand->any_array.left.length != 0)
							left = operand->any_array.left.operands[0];
						else left = operand->any_array.all;
						if (operand->any_array.right.length != 0)
							right = operand->any_array.right.operands[operand->any_array.right.length - 1];
						else right = operand->any_array.all;
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* set_definition = nullptr;
					if (left != nullptr) {
						if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
							if (left->any.included != nullptr)
								set_definition = left->any.included;
						} else if (left->type == hol_term_type::EQUALS) {
							set_definition = left->binary.right;
						} else if (left->type == hol_term_type::BINARY_APPLICATION) {
							set_definition = left->ternary.third;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (set_definition != nullptr) {
						if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
							/* no-op */
						} else if (set_definition->type == hol_term_type::LAMBDA) {
							element_variable = set_definition->quantifier.variable;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (element_variable == 0 && right != nullptr) {
						/* try to get the element variable from the right conjunct, since we couldn't get it from the left conjunct */
						if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
							right = right->any.included;
						if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
							right = right->binary.right;
						if (right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) {
							element_variable = ++max_variable;
						} else if (right->type == hol_term_type::EXISTS || right->type == hol_term_type::FOR_ALL) {
							element_variable = right->quantifier.variable;
						} else if (right->type == hol_term_type::AND && right->array.length == 2 && right->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
								&& right->array.operands[1]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[1]->binary.left->variable == set_variable
								&& right->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = right->array.operands[1]->binary.right->variable;
						} else {
							return (hol_term*) nullptr;
						}
					}
				}

				hol_term* head_var = hol_term::new_variable(element_variable);
				if (head_var == nullptr) return (hol_term*) nullptr;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return (hol_term*) nullptr;
				}
				free(*head_var);

				hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
				if (conjunct == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;

				hol_term* conjunction = nullptr;
				if (ConjunctIndex >= 0) {
					conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct, make_array_view((hol_term**) nullptr, 0),
							make_repeated_array_view(conjunct, ConjunctIndex + 1),
							make_array_view((hol_term**) nullptr, 0));
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct, make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0),
							make_repeated_array_view(conjunct, index + 1));
				}
				if (conjunction == nullptr) {
					free(*conjunct); free(conjunct);
					return (hol_term*) nullptr;
				}
				conjunction->any_array.all->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
					conjunction->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
					conjunction->any_array.right.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
					conjunction->any_array.any.operands[i]->reference_count++;
				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);

				hol_term* set_definition;
				if (left == nullptr || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
					set_definition = hol_term::new_any_right(hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::EQUALS) {
					set_definition = hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = hol_term::new_apply(&HOL_ANY, hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
					HOL_ANY.reference_count++;
				}

				hol_term* dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition, 1), make_array_view((hol_term**) nullptr, 0)));
				if (dst == nullptr) {
					free(*set_definition); free(set_definition);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, dst, head);
				free(*dst); if (dst->reference_count == 0) free(dst);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "select_conjunct_in_set ERROR: Expected intersection to be unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				hol_term* operand = intersection[0]->quantifier.operand;
				hol_term* new_left;
				if (operand->type == hol_term_type::AND) {
					new_left = operand->array.operands[0];
				} else {
					new_left = operand->any_array.left.operands[0];
				}

				hol_term* inner_operand;
				if (new_left->type == hol_term_type::ANY || new_left->type == hol_term_type::ANY_RIGHT) {
					inner_operand = new_left->any.included->quantifier.operand;
				} else if (new_left->type == hol_term_type::EQUALS) {
					inner_operand = new_left->binary.right->quantifier.operand;
				} else {
					inner_operand = new_left->ternary.third->quantifier.operand;
				}

				/* find the predicate in `inner_operand` */
				head_index inner_predicate_index;
				hol_term* predicate = find_predicate<built_in_predicates>(element_variable, inner_operand, inner_predicate_index);
				if (inner_predicate_index.position == head_position::NONE) {
					/* the "predicate" could be an equals term */
					predicate = find_variable_definition<built_in_predicates>(element_variable, inner_operand, inner_predicate_index);
					if (inner_predicate_index.position == head_position::NONE) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
				}

				hol_term* new_conjunct;
				if (inner_operand->type == hol_term_type::AND) {
					unsigned int index = (unsigned int) ((ConjunctIndex >= 0) ? ConjunctIndex : (inner_operand->array.length + ConjunctIndex));
					new_conjunct = inner_operand->array.operands[index];
					if (predicate == nullptr)
						predicate = inner_operand->array.operands[inner_predicate_index.index];
				} else if (inner_operand->type == hol_term_type::ANY_ARRAY) {
					if (ConjunctIndex >= 0) {
						new_conjunct = inner_operand->any_array.left.operands[ConjunctIndex];
					} else {
						unsigned int index = (unsigned int) (inner_operand->any_array.right.length + ConjunctIndex);
						new_conjunct = inner_operand->any_array.right.operands[index];
					}
					if (predicate == nullptr) {
						if (inner_predicate_index.position == head_position::LEFT)
							predicate = inner_operand->any_array.left.operands[inner_predicate_index.index];
						else if (inner_predicate_index.position == head_position::RIGHT)
							predicate = inner_operand->any_array.right.operands[inner_operand->any_array.right.length - inner_predicate_index.index - 1];
						else if (inner_predicate_index.position == head_position::ANY)
							predicate = inner_operand->any_array.any.operands[inner_predicate_index.index];
					}
				} else if (ConjunctIndex == 0 || ConjunctIndex == -1) {
					new_conjunct = inner_operand;
					if (predicate == nullptr && inner_predicate_index.position == head_position::LEFT && inner_predicate_index.index == 0)
						predicate = inner_operand;
				} else {
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				/* map any `^[x]` to `?[x]` where x is the element_variable */
				hol_term* temp = change_quantifier(new_conjunct, element_variable, hol_term_type::LAMBDA, hol_term_type::EXISTS);
				if (temp == nullptr) {
					free_all(intersection);
					return (hol_term*) nullptr;
				}
				new_conjunct = temp;

				dst = hol_term::new_exists(element_variable, hol_term::new_and(predicate, new_conjunct));
				if (dst == nullptr) {
					free(*new_conjunct); if (new_conjunct->reference_count == 0) free(new_conjunct);
					free_all(intersection);
					return (hol_term*) nullptr;
				}
				predicate->reference_count++;
				free_all(intersection);

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool remove_conjunct_in_set(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[lambda_variable,&max_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* left = nullptr;
				unsigned int set_variable, element_variable = 0;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT)
				{
					set_variable = ++max_variable;
					element_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS)
						fprintf(stderr, "remove_conjunct_in_set WARNING: Expected an existential quantification.\n");
#endif
					set_variable = head->quantifier.variable;
					hol_term* operand = head->quantifier.operand;

					if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
						element_variable = ++max_variable;
					} else if (operand->type == hol_term_type::AND && operand->array.length != 0) {
						left = operand->array.operands[0];
					} else if (operand->type == hol_term_type::ANY_ARRAY) {
						if (operand->any_array.left.length != 0) {
							left = operand->any_array.left.operands[0];
						} if (operand->any_array.right.length != 0) {
							hol_term* right = operand->any_array.right.operands[operand->any_array.right.length - 1];

							bool has_any_right = false;
							if (right->type == hol_term_type::ANY_RIGHT && right->any.included != nullptr) {
								has_any_right = true;
								right = right->any.included;
							}

							if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT
							 && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
								right = right->binary.right;
							while (right->type == hol_term_type::NOT)
								right = right->unary.operand;
							if (right->type == hol_term_type::FOR_ALL || right->type == hol_term_type::EXISTS) {
								element_variable = right->quantifier.variable;
							} else if (right->type == hol_term_type::AND && right->array.length == 2
									&& right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
									&& right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE
									&& right->array.operands[0]->binary.left->variable == set_variable
									&& right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
							{
								element_variable = right->array.operands[0]->binary.right->variable;
							} else if (has_any_right && right->type == hol_term_type::UNARY_APPLICATION
									&& right->binary.left->type == hol_term_type::VARIABLE
									&& right->binary.left->variable == lambda_variable
									&& right->binary.right->type == hol_term_type::VARIABLE)
							{
								element_variable = right->binary.right->variable;
							} else {
								return (hol_term*) nullptr;
							}
						}
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* set_definition = nullptr;
					if (left != nullptr) {
						if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
							if (element_variable == 0)
								element_variable = ++max_variable;
						} else if (left->type == hol_term_type::EQUALS) {
							set_definition = left->binary.right;
						} else if (left->type == hol_term_type::BINARY_APPLICATION) {
							set_definition = left->ternary.third;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (set_definition != nullptr) {
						if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
							element_variable = ++max_variable;
						} else if (set_definition->type == hol_term_type::LAMBDA) {
							element_variable = set_definition->quantifier.variable;
						} else {
							return (hol_term*) nullptr;
						}
					}
				}

				hol_term* head_var = hol_term::new_variable(element_variable);
				if (head_var == nullptr) return (hol_term*) nullptr;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return (hol_term*) nullptr;
				}
				free(*head_var);

				hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
				if (conjunct == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;

				hol_term* conjunction = nullptr;
				if (ConjunctIndex >= 0) {
					conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct, make_array_view((hol_term**) nullptr, 0),
							make_repeated_array_view(conjunct, ConjunctIndex + 1),
							make_array_view((hol_term**) nullptr, 0));
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct, make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0),
							make_repeated_array_view(conjunct, index + 1));
				}
				if (conjunction == nullptr) {
					free(*conjunct); free(conjunct);
					return (hol_term*) nullptr;
				}
				conjunction->any_array.all->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
					conjunction->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
					conjunction->any_array.right.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
					conjunction->any_array.any.operands[i]->reference_count++;
				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);

				hol_term* set_definition;
				if (left == nullptr || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
					set_definition = hol_term::new_any_right(hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::EQUALS) {
					set_definition = hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = hol_term::new_apply(&HOL_ANY, hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, conjunction));
					if (set_definition == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
					HOL_ANY.reference_count++;
				}

				hol_term* dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition, 1), make_array_view((hol_term**) nullptr, 0)));
				if (dst == nullptr) {
					free(*set_definition); free(set_definition);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, dst, head);
				free(*dst); if (dst->reference_count == 0) free(dst);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "remove_conjunct_in_set ERROR: Expected intersection to be unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				hol_term* operand = intersection[0]->quantifier.operand;
				if (operand->type == hol_term_type::AND) {
					left = operand->array.operands[0];
				} else {
					left = operand->any_array.left.operands[0];
				}

				hol_term* inner_operand;
				if (left->type == hol_term_type::ANY_RIGHT) {
					inner_operand = left->any.included->quantifier.operand;
				} else if (left->type == hol_term_type::EQUALS) {
					inner_operand = left->binary.right->quantifier.operand;
				} else {
					inner_operand = left->ternary.third->quantifier.operand;
				}

				hol_term* new_set_definition;
				if (inner_operand->type == hol_term_type::AND) {
					unsigned int index = (unsigned int) ((ConjunctIndex >= 0) ? ConjunctIndex : (inner_operand->array.length + ConjunctIndex));
					if (inner_operand->array.length == 2) {
						new_set_definition = hol_term::new_lambda(element_variable,
								(index == 0 ? inner_operand->array.operands[1] : inner_operand->array.operands[0]));
					} else {
						new_set_definition = hol_term::new_lambda(element_variable, hol_term::new_and(
								make_excluded_array_view(inner_operand->array.operands, inner_operand->array.length, index)));
					}
					if (new_set_definition == nullptr) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < inner_operand->array.length; i++)
						if (i != index) inner_operand->array.operands[i]->reference_count++;
				} else if (inner_operand->type == hol_term_type::ANY_ARRAY) {
					if (ConjunctIndex >= 0) {
						new_set_definition = hol_term::new_lambda(element_variable, hol_term::new_any_array(
								hol_term_type::AND, inner_operand->any_array.all, make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
								make_excluded_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length, ConjunctIndex),
								make_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length)));
					} else {
						new_set_definition = hol_term::new_lambda(element_variable, hol_term::new_any_array(
								hol_term_type::AND, inner_operand->any_array.all, make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
								make_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length),
								make_excluded_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length, (unsigned int) (inner_operand->any_array.right.length + ConjunctIndex))));
					}
					if (new_set_definition == nullptr) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					hol_term* new_conjunction = new_set_definition->quantifier.operand;
					new_conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.any.length; i++)
						new_conjunction->any_array.any.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.left.length; i++)
						new_conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.right.length; i++)
						new_conjunction->any_array.right.operands[i]->reference_count++;
				} else {
					new_set_definition = hol_term::new_lambda(element_variable, hol_term::new_true());
				}

				hol_term* new_left;
				if (left->type == hol_term_type::ANY_RIGHT) {
					new_left = hol_term::new_any_right(new_set_definition, left->any.excluded_trees, left->any.excluded_tree_count);
					if (new_left == nullptr) {
						free(*new_set_definition); free(new_set_definition);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < left->any.excluded_tree_count; i++)
						left->any.excluded_trees[i]->reference_count++;
				} else if (left->type == hol_term_type::EQUALS) {
					new_left = hol_term::new_equals(left->binary.left, new_set_definition);
					if (new_left == nullptr) {
						free(*new_set_definition); free(new_set_definition);
						return (hol_term*) nullptr;
					}
					left->binary.left->reference_count++;
				} else {
					new_left = hol_term::new_apply(left->ternary.first, left->ternary.second, new_set_definition);
					if (new_left == nullptr) {
						free(*new_set_definition); free(new_set_definition);
						return (hol_term*) nullptr;
					}
					left->ternary.first->reference_count++;
					left->ternary.second->reference_count++;
				}

				if (operand->type == hol_term_type::AND) {
					dst = hol_term::new_exists(set_variable, hol_term::new_and(
							make_prepended_array_view(new_left, make_array_view(operand->array.operands + 1, operand->array.length - 1))));
					if (dst == nullptr) {
						free_all(intersection);
						free(*new_left); free(new_left);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 1; i < operand->array.length; i++)
						operand->array.operands[i]->reference_count++;
				} else {
					dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND,
							operand->any_array.all, make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
							make_prepended_array_view(new_left, make_array_view(operand->any_array.left.operands + 1, operand->any_array.left.length - 1)),
							make_array_view(operand->any_array.right.operands, operand->any_array.right.length)));
					if (dst == nullptr) {
						free_all(intersection);
						free(*new_left); free(new_left);
						return (hol_term*) nullptr;
					}
					operand->any_array.all->reference_count++;
					for (unsigned int i = 0; i < operand->any_array.any.length; i++)
						operand->any_array.any.operands[i]->reference_count++;
					for (unsigned int i = 1; i < operand->any_array.left.length; i++)
						operand->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < operand->any_array.right.length; i++)
						operand->any_array.right.operands[i]->reference_count++;
				}
				free_all(intersection);

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool select_subset_in_set(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array_map<unsigned int, hol_term*> scopes(8);
	auto gather_scope = [&scopes](hol_term* term) {
		if (term->type == hol_term_type::ANY || term->type == hol_term_type::ANY_RIGHT) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = 0;
			scopes.values[scopes.size++] = term;
		} else if (term->type == hol_term_type::FOR_ALL || term->type == hol_term_type::EXISTS || term->type == hol_term_type::LAMBDA) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = term->quantifier.variable;
			scopes.values[scopes.size++] = term;
		}
		return true;
	};

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[lambda_variable,&max_variable,&scopes](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				if (head->type == hol_term_type::NOT)
					head = head->unary.operand;

				bool has_new_set = false;
				unsigned int element_variable = 0, new_set_variable;
				hol_term* left = nullptr;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT)
				{
					element_variable = ++max_variable;
					new_set_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS)
						fprintf(stderr, "select_subset_in_set WARNING: Expected an existential quantification.\n");
#endif
					hol_term* operand = head->quantifier.operand;

					if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
						element_variable = ++max_variable;
						new_set_variable = ++max_variable;
					} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
						if (operand->any_array.left.length == 0) {
							left = operand->any_array.all;
						} else {
							left = operand->any_array.left.operands[0];
						}
					} else if (operand->type == hol_term_type::AND) {
						left = operand->array.operands[0];
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* set_definition = nullptr;
					if (left != nullptr) {
						/* make sure `left` can be of the form `subset(*,*)` */
						hol_term* subset_term = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SUBSET), &HOL_ANY, &HOL_ANY);
						if (subset_term == nullptr)
							return (hol_term*) nullptr;
						if (!has_intersection<built_in_predicates>(left, subset_term)) {
							free(*subset_term); free(subset_term);
							return (hol_term*) nullptr;
						}
						free(*subset_term); free(subset_term);

						if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
							if (left->any.included != nullptr) {
								set_definition = left->any.included;
							} else {
								new_set_variable = ++max_variable;
							}
						} else if (left->type == hol_term_type::BINARY_APPLICATION) {
							set_definition = left->ternary.third;
						} else {
							return (hol_term*) nullptr;
						}
					}

					hol_term* inner_operand = nullptr;
					if (set_definition != nullptr) {
						if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
							new_set_variable = ++max_variable;
						} else if (set_definition->type == hol_term_type::ANY_QUANTIFIER && has_intersection(set_definition->any_quantifier.quantifier, hol_quantifier_type::LAMBDA)) {
							new_set_variable = ++max_variable;
							inner_operand = set_definition->any_quantifier.operand;
						} else if (set_definition->type == hol_term_type::LAMBDA) {
							element_variable = set_definition->quantifier.variable;
							inner_operand = set_definition->quantifier.operand;
						} else {
							return (hol_term*) nullptr;
						}
					}

					hol_term* subset_term = nullptr;
					if (inner_operand != nullptr) {
						if (inner_operand->type == hol_term_type::ANY || inner_operand->type == hol_term_type::ANY_RIGHT) {
							new_set_variable = ++max_variable;
						} else if (inner_operand->type == hol_term_type::ANY_ARRAY) {
							if (ConjunctIndex >= 0) {
								if (ConjunctIndex < inner_operand->any_array.left.length)
									subset_term = inner_operand->any_array.left.operands[ConjunctIndex];
								else subset_term = inner_operand->any_array.all;
							} else {
								if (-ConjunctIndex - 1 < inner_operand->any_array.right.length)
									subset_term = inner_operand->any_array.right.operands[inner_operand->any_array.right.length + ConjunctIndex];
								else subset_term = inner_operand->any_array.all;
							}
						} else if (inner_operand->type == hol_term_type::AND) {
							if (ConjunctIndex >= 0 && (unsigned int) ConjunctIndex >= inner_operand->array.length)
								return (hol_term*) nullptr;
							else if (ConjunctIndex < 0 && (unsigned int) (-ConjunctIndex) > inner_operand->array.length)
								return (hol_term*) nullptr;
							unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (unsigned int) (inner_operand->array.length + ConjunctIndex);
							subset_term = inner_operand->array.operands[index];
						} else if (ConjunctIndex == 0 || ConjunctIndex == -1) {
							subset_term = inner_operand;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (subset_term == nullptr || subset_term->type == hol_term_type::ANY || subset_term->type == hol_term_type::ANY_RIGHT) {
						new_set_variable = ++max_variable;
					} else if (subset_term->type == hol_term_type::UNARY_APPLICATION) {
						if (subset_term->binary.left->type == hol_term_type::VARIABLE) {
							has_new_set = true;
							new_set_variable = subset_term->binary.left->variable;
						} else {
							new_set_variable = ++max_variable;
						}
					} else {
						return (hol_term*) nullptr;
					}

					if (element_variable == 0) {
						/* try to get the variable from the right conjunct in the set scope */
						hol_term* right = nullptr;
						if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
							if (operand->any_array.right.length == 0) {
								right = operand->any_array.all;
							} else {
								right = operand->any_array.right.operands[operand->any_array.right.length - 1];
							}
						} else if (operand->type == hol_term_type::AND) {
							right = operand->array.operands[operand->array.length - 1];
						}

						if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
							right = right->any.included;

						hol_term* lambda_application = nullptr;
						if (right->type == hol_term_type::FOR_ALL && right->quantifier.operand->type == hol_term_type::IF_THEN) {
							lambda_application = right->quantifier.operand->binary.right;
						} else if (right->type == hol_term_type::EXISTS && right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.length == 2) {
							lambda_application = right->quantifier.operand->array.operands[1];
						} else if (right->type == hol_term_type::AND && right->array.length == 2) {
							lambda_application = right->array.operands[1];
						} else {
							lambda_application = right;
						}

						if ((lambda_application->type == hol_term_type::ANY || lambda_application->type == hol_term_type::ANY_RIGHT) && lambda_application->any.included != nullptr)
							lambda_application = lambda_application->any.included;

						if (lambda_application->type == hol_term_type::UNARY_APPLICATION
						 && lambda_application->binary.left->type == hol_term_type::VARIABLE
						 && lambda_application->binary.left->variable == lambda_variable
						 && lambda_application->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = lambda_application->binary.right->variable;
						} else {
							return (hol_term*) nullptr;
						}
					}
				}

				if (has_new_set) {
					/* the new set is already defined, so we don't need to define it here */
					hol_term* new_set_scope = scopes.get(new_set_variable)->quantifier.operand;

					/* get the element variable of this set */
					hol_term* left = nullptr;
					unsigned int new_element_variable = 0;
					if (new_set_scope->type == hol_term_type::ANY || new_set_scope->type == hol_term_type::ANY_RIGHT) {
						new_element_variable = element_variable;
					} else if (new_set_scope->type == hol_term_type::ANY_ARRAY && new_set_scope->any_array.oper == hol_term_type::AND) {
						if (new_set_scope->any_array.left.length == 0)
							left = new_set_scope->any_array.all;
						else left = new_set_scope->any_array.left.operands[0];
					} else if (new_set_scope->type == hol_term_type::AND) {
						left = new_set_scope->array.operands[0];
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* new_set_definition = nullptr;
					if (left != nullptr) {
						if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
							new_element_variable = element_variable;
						} else if (left->type == hol_term_type::EQUALS) {
							new_set_definition = left->binary.right;
						} else if (left->type == hol_term_type::BINARY_APPLICATION) {
							new_set_definition = left->ternary.third;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (new_set_definition != nullptr) {
						if (new_set_definition->type == hol_term_type::ANY || new_set_definition->type == hol_term_type::ANY_RIGHT
						 || (new_set_definition->type == hol_term_type::ANY_QUANTIFIER && has_intersection(new_set_definition->any_quantifier.quantifier, hol_quantifier_type::LAMBDA)))
						{
							new_element_variable = element_variable;
						} else if (new_set_definition->type == hol_term_type::LAMBDA) {
							new_element_variable = new_set_definition->quantifier.variable;
						} else {
							return (hol_term*) nullptr;
						}
					}

					hol_term* new_element_var = hol_term::new_variable(new_element_variable);
					if (new_element_var == nullptr)
						return (hol_term*) nullptr;
					hol_term* dst = hol_term::new_exists(new_element_variable, hol_term::new_and(
							hol_term::new_apply(hol_term::new_variable(new_set_variable), new_element_var),
							hol_term::new_apply(hol_term::new_variable(lambda_variable), new_element_var)));
					if (dst == nullptr) {
						free(*new_element_var); free(new_element_var);
						return (hol_term*) nullptr;
					}
					new_element_var->reference_count += 2 - 1;
					return dst;
				} else {
					unsigned int new_element_variable = element_variable;
					hol_term* new_element_var = hol_term::new_variable(new_element_variable);
					if (new_element_var == nullptr)
						return (hol_term*) nullptr;
					hol_term* expected_right = hol_term::new_exists(new_element_variable, hol_term::new_and(
							hol_term::new_apply(hol_term::new_variable(new_set_variable), new_element_var),
							hol_term::new_apply(hol_term::new_variable(lambda_variable), new_element_var)));
					if (expected_right == nullptr) {
						free(*new_element_var); free(new_element_var);
						return (hol_term*) nullptr;
					}
					new_element_var->reference_count += 2 - 1;

					hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(new_set_variable, &HOL_ANY));
					if (excluded_quantifier == nullptr) {
						free(*expected_right); free(expected_right);
						return (hol_term*) nullptr;
					}

					hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
					if (expected_conjunct == nullptr) {
						free(*expected_right); free(expected_right);
						free(*excluded_quantifier); free(excluded_quantifier);
						return (hol_term*) nullptr;
					}

					hol_term* dst = hol_term::new_exists(new_set_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&expected_right, 1)));
					if (dst == nullptr) {
						free(*expected_right); free(expected_right);
						free(*expected_conjunct); free(expected_conjunct);
						return (hol_term*) nullptr;
					}
					return dst;
				}
			}, gather_scope);
}

template<int_fast8_t PredicateIndex>
inline bool require_predicate(
		hol_term* src, hol_term*& dst,
		hol_term* predicate)
{
	array_map<unsigned int, hol_term*> predicates(4);

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	bool result = apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[predicate, &predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if (!predicates.ensure_capacity(predicates.size + 1))
					return (hol_term*) nullptr;
				hol_term* current_predicate;
				unsigned int index = predicates.index_of(head_variable);
				if (index == predicates.size) {
					array_map<unsigned int, unsigned int> variable_map(1);
					variable_map.put(0, head_variable);
					predicates.keys[index] = head_variable;
					predicates.values[index] = map_variables(predicate, variable_map);
					if (predicates.values[index] == nullptr)
						return (hol_term*) nullptr;
					predicates.size++;
				}
				current_predicate = predicates.values[index];

				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					hol_term* head_var = hol_term::new_variable(head_variable);
					if (head_var == nullptr) return (hol_term*) nullptr;
					constexpr unsigned int excluded_tree_count = 4;
					hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term) * (excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()));
					excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
					excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
					excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
					excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
					if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
						if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
						if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
						if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
						free(*head_var); free(head_var);
						return (hol_term*) nullptr;
					}
					free(*head_var);

					for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++) {
						excluded_trees[excluded_tree_count + i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
						excluded_trees[excluded_tree_count + i]->reference_count++;
					}

					hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
					if (conjunct == nullptr) {
						for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
							free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
						}
						return (hol_term*) nullptr;
					}

					hol_term* conjunction;
					if (!new_hol_term(conjunction)) return (hol_term*) nullptr;
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = conjunct;
					if (PredicateIndex >= 0 && PredicateIndex != INT_FAST8_MAX) {
						conjunction->any_array.left.length = PredicateIndex + 1;
						conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * (PredicateIndex + 1));
						if (conjunction->any_array.left.operands == nullptr) {
							free(conjunction); free(*conjunct); free(conjunct);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < PredicateIndex; i++)
							conjunction->any_array.left.operands[i] = conjunct;
						conjunction->any_array.left.operands[PredicateIndex] = current_predicate;
						conjunction->any_array.right.length = 0;
						conjunction->any_array.right.operands = nullptr;
						conjunction->any_array.any.length = 0;
						conjunction->any_array.any.operands = nullptr;
						conjunct->reference_count += 1 + PredicateIndex;
					} else if (PredicateIndex < 0) {
						unsigned int index = (unsigned int) (-PredicateIndex) - 1;
						conjunction->any_array.right.length = index + 1;
						conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * (index + 1));
						if (conjunction->any_array.right.operands == nullptr) {
							free(conjunction); free(*conjunct); free(conjunct);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 1; i < index + 1; i++)
							conjunction->any_array.right.operands[i] = conjunct;
						conjunction->any_array.right.operands[0] = current_predicate;
						conjunction->any_array.left.length = 0;
						conjunction->any_array.left.operands = nullptr;
						conjunction->any_array.any.length = 0;
						conjunction->any_array.any.operands = nullptr;
						conjunct->reference_count += 1 + index;
					} else if (PredicateIndex == INT_FAST8_MAX) {
						conjunction->any_array.left.length = 0;
						conjunction->any_array.left.operands = nullptr;
						conjunction->any_array.right.length = 0;
						conjunction->any_array.right.operands = nullptr;
						conjunction->any_array.any.length = 1;
						conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * 1);
						if (conjunction->any_array.any.operands == nullptr) {
							free(conjunction); free(*conjunct); free(conjunct);
							return (hol_term*) nullptr;
						}
						conjunction->any_array.any.operands[0] = current_predicate;
						conjunct->reference_count += 1;
					} else {
						fprintf(stderr, "require_predicate ERROR: Unsupported valueof `PredicateIndex`.\n");
						free(conjunction); free(*conjunct); free(conjunct);
						return (hol_term*) nullptr;
					}
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
					current_predicate->reference_count++;

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "require_predicate ERROR: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count == 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;
					if (PredicateIndex != INT_FAST8_MAX && predicate_index.compare(PredicateIndex))
						return (hol_term*) nullptr;

					hol_term* predicate = nullptr;
					if (predicate_index.position == head_position::LEFT) {
						predicate = operand->any_array.left.operands[predicate_index.index];
					} else if (predicate_index.position == head_position::RIGHT) {
						predicate = operand->any_array.right.operands[operand->any_array.right.length - predicate_index.index - 1];
					} else {
						predicate = operand->any_array.any.operands[predicate_index.index];
					}
					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, predicate, current_predicate);
					if (intersection.length > 1) {
						fprintf(stderr, "require_predicate ERROR: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return (hol_term*) nullptr;
					}
					if (intersection.length == 0)
						return (hol_term*) nullptr;

					if (intersection[0] == predicate) {
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						dst = head;
						head->reference_count++;
					} else {
						hol_term* conjunction;
						if (!new_hol_term(conjunction)) {
							for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
							return (hol_term*) nullptr;
						}
						conjunction->type = hol_term_type::ANY_ARRAY;
						conjunction->reference_count = 1;
						conjunction->any_array.oper = hol_term_type::AND;
						conjunction->any_array.all = operand->any_array.all;
						conjunction->any_array.all->reference_count++;
						conjunction->any_array.left.length = 0;
						conjunction->any_array.left.operands = nullptr;
						conjunction->any_array.right.length = 0;
						conjunction->any_array.right.operands = nullptr;
						conjunction->any_array.any.length = 0;
						conjunction->any_array.any.operands = nullptr;

						conjunction->any_array.left.length = operand->any_array.left.length;
						if (conjunction->any_array.left.length == 0) {
							conjunction->any_array.left.operands = nullptr;
						} else {
							conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.left.length);
							if (conjunction->any_array.left.operands == nullptr) {
								for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunction); free(conjunction);
								return (hol_term*) nullptr;
							}
							for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
								if (predicate_index.position == head_position::LEFT && predicate_index.index == i) {
									conjunction->any_array.left.operands[i] = intersection[0];
								} else {
									conjunction->any_array.left.operands[i] = operand->any_array.left.operands[i];
								}
								conjunction->any_array.left.operands[i]->reference_count++;
							}
						}
						conjunction->any_array.right.length = operand->any_array.right.length;
						if (conjunction->any_array.right.length == 0) {
							conjunction->any_array.right.operands = nullptr;
						} else {
							conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.right.length);
							if (conjunction->any_array.right.operands == nullptr) {
								for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunction); free(conjunction);
								return (hol_term*) nullptr;
							}
							for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
								if (predicate_index.position == head_position::RIGHT && predicate_index.index == conjunction->any_array.right.length - i - 1) {
									conjunction->any_array.right.operands[i] = intersection[0];
								} else {
									conjunction->any_array.right.operands[i] = operand->any_array.right.operands[i];
								}
								conjunction->any_array.right.operands[i]->reference_count++;
							}
						}
						conjunction->any_array.any.length = operand->any_array.any.length;
						if (conjunction->any_array.any.length == 0) {
							conjunction->any_array.any.operands = nullptr;
						} else {
							conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.any.length);
							if (conjunction->any_array.any.operands == nullptr) {
								for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunction); free(conjunction);
								return (hol_term*) nullptr;
							}
							for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
								if (predicate_index.position == head_position::ANY && predicate_index.index == i) {
									conjunction->any_array.any.operands[i] = intersection[0];
								} else {
									conjunction->any_array.any.operands[i] = operand->any_array.any.operands[i];
								}
								conjunction->any_array.any.operands[i]->reference_count++;
							}
						}
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }

						dst = hol_term::new_exists(head->quantifier.variable, conjunction);
						if (dst == nullptr) {
							free(*conjunction); free(conjunction);
							return (hol_term*) nullptr;
						}
					}
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS)
						fprintf(stderr, "require_predicate WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					if (PredicateIndex != INT_FAST8_MAX) {
						int expected_predicate_index = PredicateIndex;
						if (expected_predicate_index < 0)
							expected_predicate_index += operand->array.length;
						if (!predicate_index.compare(expected_predicate_index))
							return (hol_term*) nullptr;
					}

					if (operand->type != hol_term_type::AND && predicate_index.index != 0)
						return (hol_term*) nullptr;

					hol_term* predicate = (operand->type == hol_term_type::AND ? operand->array.operands[predicate_index.index] : operand);
					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, predicate, current_predicate);
					if (intersection.length > 1) {
						fprintf(stderr, "require_predicate ERROR: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return (hol_term*) nullptr;
					}
					if (intersection.length == 0)
						return (hol_term*) nullptr;

					if (intersection[0] == predicate) {
						free_all(intersection);
						dst = head;
						head->reference_count++;
					} else if (operand->type == hol_term_type::AND) {
						hol_term* conjunction;
						if (!new_hol_term(conjunction)) {
							free_all(intersection);
							return (hol_term*) nullptr;
						}
						conjunction->type = hol_term_type::AND;
						conjunction->reference_count = 1;
						conjunction->array.length = operand->array.length;
						conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * operand->array.length);
						if (conjunction->array.operands == nullptr) {
							fprintf(stderr, "require_predicate ERROR: Out of memory.\n");
							free_all(intersection); return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < operand->array.length; i++) {
							if (i == predicate_index.index)
								conjunction->array.operands[i] = intersection[0];
							else
								conjunction->array.operands[i] = operand->array.operands[i];
							conjunction->array.operands[i]->reference_count++;
						}
						free_all(intersection);

						dst = hol_term::new_exists(head->quantifier.variable, conjunction);
						if (dst == nullptr) {
							free(*conjunction); free(conjunction);
							return (hol_term*) nullptr;
						}
					} else {
						dst = hol_term::new_exists(head->quantifier.variable, intersection[0]);
						if (dst == nullptr) {
							free_all(intersection);
							return (hol_term*) nullptr;
						}
					}
				}

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst, old_head->any.excluded_trees, old_head->any.excluded_tree_count);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < old_head->any.excluded_tree_count; i++)
						old_head->any.excluded_trees[i]->reference_count++;
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
	for (auto entry : predicates) {
		free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
	}
	return result;
}

inline bool require_no_inverse(hol_term* src, hol_term*& dst)
{
	hol_term** excluded_constants = (hol_term**) malloc(sizeof(hol_term*) * (hol_non_head_constants<built_in_predicates>::count() + 1));
	if (excluded_constants == nullptr) {
		fprintf(stderr, "require_no_inverse ERROR: Out of memory.\n");
		return false;
	}
	for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++)
		excluded_constants[i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
	excluded_constants[hol_non_head_constants<built_in_predicates>::count()] = hol_term::new_constant((unsigned int) built_in_predicates::INVERSE);
	if (excluded_constants[hol_non_head_constants<built_in_predicates>::count()] == nullptr) {
		free(excluded_constants);
		return false;
	}

	hol_term* predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, excluded_constants, hol_non_head_constants<built_in_predicates>::count() + 1),
			hol_term::new_variable(0));
	if (predicate == nullptr) {
		free(*excluded_constants[hol_non_head_constants<built_in_predicates>::count()]);
		if (excluded_constants[hol_non_head_constants<built_in_predicates>::count()]->reference_count == 0)
			free(excluded_constants[hol_non_head_constants<built_in_predicates>::count()]);
		free(excluded_constants);
		return false;
	}
	hol_non_head_constants<built_in_predicates>::increment_terms();
	free(excluded_constants);

	bool result = require_predicate<INT_FAST8_MAX>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool require_left_predicate_inverse(hol_term* src, hol_term*& dst)
{
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_apply(
				hol_term::new_constant((unsigned int) built_in_predicates::INVERSE),
				hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count())),
			hol_term::new_variable(0));
	if (predicate == nullptr) return false;
	hol_non_head_constants<built_in_predicates>::increment_terms();

	bool result = require_predicate<0>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool require_left_predicate_inverse_own(hol_term* src, hol_term*& dst)
{
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_apply(
				hol_term::new_constant((unsigned int) built_in_predicates::INVERSE),
				hol_term::new_constant((unsigned int) built_in_predicates::OWN)),
			hol_term::new_variable(0));
	if (predicate == nullptr) return false;

	bool result = require_predicate<0>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool require_left_predicate_exist(hol_term* src, hol_term*& dst)
{
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_constant((unsigned int) built_in_predicates::EXIST),
			hol_term::new_variable(0));
	if (predicate == nullptr) return false;

	bool result = require_predicate<0>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool require_no_left_predicate_exist(hol_term* src, hol_term*& dst)
{
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_any_constant_except((unsigned int) built_in_predicates::EXIST),
			hol_term::new_variable(0));
	if (predicate == nullptr) return false;

	bool result = require_predicate<0>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool remove_inverse(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					dst = head;
					dst->reference_count++;
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;

					hol_term* predicate;
					if (predicate_index.position == head_position::LEFT) {
						predicate = operand->any_array.left.operands[predicate_index.index];
					} else if (predicate_index.position == head_position::RIGHT) {
						predicate = operand->any_array.right.operands[operand->any_array.right.length - predicate_index.index - 1];
					} else {
						predicate = operand->any_array.any.operands[predicate_index.index];
					}

					hol_term* expected_predicate = hol_term::new_apply(
							hol_term::new_apply(
								hol_term::new_constant((unsigned int) built_in_predicates::INVERSE),
								&HOL_ANY),
							hol_term::new_variable(head_variable));
					if (expected_predicate == nullptr) return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					if (!has_intersection<built_in_predicates>(expected_predicate, predicate)) {
						free(*expected_predicate); free(expected_predicate);
						return (hol_term*) nullptr;
					}
					free(*expected_predicate); free(expected_predicate);

					hol_term* new_predicate = nullptr;
					if (predicate->type == hol_term_type::UNARY_APPLICATION) {
						if (predicate->binary.left->type == hol_term_type::UNARY_APPLICATION) {
							new_predicate = hol_term::new_apply(predicate->binary.left->binary.right, predicate->binary.right);
							if (new_predicate == nullptr) return (hol_term*) nullptr;
							new_predicate->binary.left->reference_count++;
							new_predicate->binary.right->reference_count++;
						} else if (predicate->binary.left->type == hol_term_type::ANY) {
							new_predicate = predicate;
							new_predicate->reference_count++;
						}
					} else if (predicate->type == hol_term_type::ANY) {
						new_predicate = predicate;
						new_predicate->reference_count++;
					}

					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						free(*new_predicate); free(new_predicate);
						return (hol_term*) nullptr;
					}
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = operand->any_array.all;
					conjunction->any_array.all->reference_count++;
					conjunction->any_array.left.length = 0;
					conjunction->any_array.left.operands = nullptr;
					conjunction->any_array.right.length = 0;
					conjunction->any_array.right.operands = nullptr;
					conjunction->any_array.any.length = 0;
					conjunction->any_array.any.operands = nullptr;

					conjunction->any_array.left.length = operand->any_array.left.length;
					if (conjunction->any_array.left.length == 0) {
						conjunction->any_array.left.operands = nullptr;
					} else {
						conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.left.length);
						if (conjunction->any_array.left.operands == nullptr) {
							free(*conjunction); free(conjunction);
							free(*new_predicate); free(new_predicate);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
							if (predicate_index.position == head_position::LEFT && predicate_index.index == i) {
								conjunction->any_array.left.operands[i] = new_predicate;
							} else {
								conjunction->any_array.left.operands[i] = operand->any_array.left.operands[i];
							}
							conjunction->any_array.left.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.right.length = operand->any_array.right.length;
					if (conjunction->any_array.right.length == 0) {
						conjunction->any_array.right.operands = nullptr;
					} else {
						conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.right.length);
						if (conjunction->any_array.right.operands == nullptr) {
							free(*conjunction); free(conjunction);
							free(*new_predicate); free(new_predicate);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
							if (predicate_index.position == head_position::RIGHT && predicate_index.index == conjunction->any_array.right.length - i - 1) {
								conjunction->any_array.right.operands[i] = new_predicate;
							} else {
								conjunction->any_array.right.operands[i] = operand->any_array.right.operands[i];
							}
							conjunction->any_array.right.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.any.length = operand->any_array.any.length;
					if (conjunction->any_array.any.length == 0) {
						conjunction->any_array.any.operands = nullptr;
					} else {
						conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.any.length);
						if (conjunction->any_array.any.operands == nullptr) {
							free(*conjunction); free(conjunction);
							free(*new_predicate); free(new_predicate);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
							if (predicate_index.position == head_position::ANY && predicate_index.index == i) {
								conjunction->any_array.any.operands[i] = new_predicate;
							} else {
								conjunction->any_array.any.operands[i] = operand->any_array.any.operands[i];
							}
							conjunction->any_array.any.operands[i]->reference_count++;
						}
					}
					free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}

				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
						fprintf(stderr, "remove_inverse WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* expected_predicate = hol_term::new_apply(
							hol_term::new_apply(
								hol_term::new_constant((unsigned int) built_in_predicates::INVERSE),
								&HOL_ANY),
							hol_term::new_variable(head_variable));
					if (expected_predicate == nullptr) return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					hol_term* operand = head->quantifier.operand;
					hol_term* predicate = operand->array.operands[predicate_index.index];
					if (!has_intersection<built_in_predicates>(expected_predicate, predicate)) {
						free(*expected_predicate); free(expected_predicate);
						return (hol_term*) nullptr;
					}
					free(*expected_predicate); free(expected_predicate);

					hol_term* new_predicate = nullptr;
					if (predicate->type == hol_term_type::UNARY_APPLICATION) {
						if (predicate->binary.left->type == hol_term_type::UNARY_APPLICATION) {
							new_predicate = hol_term::new_apply(predicate->binary.left->binary.right, predicate->binary.right);
							if (new_predicate == nullptr) return (hol_term*) nullptr;
							new_predicate->binary.left->reference_count++;
							new_predicate->binary.right->reference_count++;
						} else if (predicate->binary.left->type == hol_term_type::ANY) {
							new_predicate = predicate;
							new_predicate->reference_count++;
						}
					} else if (predicate->type == hol_term_type::ANY) {
						new_predicate = predicate;
						new_predicate->reference_count++;
					}

					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						free(*new_predicate); free(new_predicate);
						return (hol_term*) nullptr;
					}
					conjunction->type = hol_term_type::AND;
					conjunction->reference_count = 1;
					conjunction->array.length = operand->array.length;
					conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * operand->array.length);
					if (conjunction->array.operands == nullptr) {
						fprintf(stderr, "remove_inverse ERROR: Out of memory.\n");
						free(*new_predicate); free(new_predicate); free(conjunction);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < operand->array.length; i++) {
						if (i == predicate_index.index) {
							conjunction->array.operands[i] = new_predicate;
						} else {
							conjunction->array.operands[i] = operand->array.operands[i];
							conjunction->array.operands[i]->reference_count++;
						}
					}

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				}

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline hol_term* make_selected_arg_without_head(
		unsigned int negation_count, unsigned int element_variable,
		unsigned int set_variable, unsigned int lambda_variable)
{
	hol_term* element_var = hol_term::new_variable(element_variable);
	if (element_var == nullptr)
		return nullptr;
	hol_term* lambda_var = hol_term::new_variable(lambda_variable);
	if (lambda_var == nullptr) {
		free(*element_var); free(element_var);
		return nullptr;
	}

	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
	if (excluded_quantifier == nullptr) {
		free(*element_var); free(element_var);
		free(*lambda_var); free(lambda_var);
		return nullptr;
	}
	HOL_ANY.reference_count++;

	/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
	hol_term* quantified_term;
	if (negation_count > 0)
		quantified_term = hol_term::new_any_right(hol_term::new_apply(lambda_var, element_var));
	else quantified_term = hol_term::new_any_right(hol_term::new_apply(lambda_var, element_var), &excluded_quantifier, 1);
	if (quantified_term == nullptr) {
		free(*element_var); free(element_var);
		free(*lambda_var); free(lambda_var);
		free(*excluded_quantifier); free(excluded_quantifier);
		return nullptr;
	}

	for (unsigned int i = 0; i < negation_count; i++) {
		hol_term* temp = hol_term::new_not(quantified_term);
		if (temp == nullptr) {
			free(*quantified_term); free(quantified_term);
			if (negation_count > 0) { free(*excluded_quantifier); free(excluded_quantifier); }
			return nullptr;
		}
		quantified_term = temp;
	}

	hol_term* right;
	if (negation_count > 0)
		right = hol_term::new_any_right(quantified_term, &excluded_quantifier, 1);
	else right = quantified_term;
	if (right == nullptr) {
		free(*quantified_term); free(quantified_term);
		if (negation_count > 0) { free(*excluded_quantifier); free(excluded_quantifier); }
	}

	hol_term* dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&right, 1)));
	if (dst == nullptr) {
		free(*quantified_term); free(quantified_term);
		return nullptr;
	}
	HOL_ANY.reference_count++;
	return dst;
}

struct scope_info {
	hol_term* scope;
	bool narrow_scope;
	unsigned int negations;
};

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant>
inline bool select_arg_without_head_predicative(
		hol_term* src, hol_term*& dst)
{
	unsigned int lambda_variable = 0;
	max_bound_variable(*src, lambda_variable);
	lambda_variable++;

	bool narrow_scope = false;
	unsigned int parent_negations = 0;
	hol_term* next_universal_scope = nullptr;
	hol_term* parent = nullptr;
	array_map<unsigned int, scope_info> scopes(8);
	auto gather_scope = [&parent_negations,&scopes,&narrow_scope,&next_universal_scope,&parent](hol_term* term) {
		if (term->type == hol_term_type::ANY || term->type == hol_term_type::ANY_RIGHT) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = 0;
			scopes.values[scopes.size++] = {term, narrow_scope, parent_negations};
		} else if (term->type == hol_term_type::FOR_ALL || term->type == hol_term_type::EXISTS || term->type == hol_term_type::LAMBDA) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = term->quantifier.variable;
			scopes.values[scopes.size++] = {term, narrow_scope, parent_negations};
		} else if (term->type == hol_term_type::UNARY_APPLICATION && term->binary.left->type == hol_term_type::CONSTANT && term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
			narrow_scope = true;
			if (term->binary.right->type == hol_term_type::FOR_ALL)
				next_universal_scope = term->binary.right;
		}
		if (term->type == hol_term_type::NOT)
			parent_negations++;
		else parent_negations = 0;
		parent = term;
		return true;
	};

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	bool result = apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[&scopes,lambda_variable,&narrow_scope,&next_universal_scope,&parent](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				unsigned int negation_count = 0;
				while (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negation_count++;
				}

				hol_term* dst = nullptr;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					unsigned int set_variable = lambda_variable + 1;
					unsigned int element_variable = lambda_variable + 2;
					hol_term* new_term = make_selected_arg_without_head(negation_count, element_variable, set_variable, lambda_variable);
					if (new_term == nullptr) return (hol_term*) nullptr;
					if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
						dst = hol_term::new_any_right(new_term);
						if (dst == nullptr) {
							free(*new_term); if (new_term->reference_count == 0) free(new_term);
							return (hol_term*) nullptr;
						}
					} else {
						dst = new_term;
					}
					remove_wide_scope_marker = true;
				} else if (head->type == hol_term_type::EXISTS) {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunct;
					if (operand->type == hol_term_type::ANY_ARRAY) {
						if (ConjunctIndex >= 0) {
							if (ConjunctIndex < operand->any_array.left.length)
								conjunct = operand->any_array.left.operands[ConjunctIndex];
							else conjunct = operand->any_array.all;
						} else if (ConjunctIndex < 0) {
							unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
							if (index < operand->any_array.right.length)
								conjunct = operand->any_array.right.operands[operand->any_array.right.length - index - 1];
							else conjunct = operand->any_array.all;
						}
					} else if (operand->type == hol_term_type::AND) {
						int conjunct_index = ConjunctIndex;
						if (ConjunctIndex < 0)
							conjunct_index += operand->array.length;
						conjunct = operand->array.operands[conjunct_index];
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* expected_left_side = hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(head->quantifier.variable));
					if (expected_left_side == nullptr)
						return (hol_term*) nullptr;
					hol_term* expected_conjunct = hol_term::new_any_right(hol_term::new_equals(expected_left_side, &HOL_ANY));
					if (expected_left_side == nullptr) {
						free(*expected_left_side); free(expected_left_side);
						return (hol_term*) nullptr;
					}
					expected_left_side->reference_count++;
					HOL_ANY.reference_count++;
					if (!has_intersection<built_in_predicates>(expected_conjunct, conjunct)) {
						free(*expected_conjunct); free(expected_conjunct);
						free(*expected_left_side); free(expected_left_side);
						return (hol_term*) nullptr;
					}
					free(*expected_conjunct); free(expected_conjunct);

					if (conjunct->type == hol_term_type::ANY || conjunct->type == hol_term_type::ANY_RIGHT) {
						free(*expected_left_side); free(expected_left_side);
						unsigned int set_variable = lambda_variable + 1;
						unsigned int element_variable = lambda_variable + 2;
						dst = make_selected_arg_without_head(negation_count, element_variable, set_variable, lambda_variable);
						if (dst == nullptr) return (hol_term*) nullptr;
						if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
							hol_term* new_dst = hol_term::new_any_right(dst);
							if (dst == nullptr) {
								free(*dst); if (dst->reference_count == 0) free(dst);
								return (hol_term*) nullptr;
							}
							dst = new_dst;
						}
						remove_wide_scope_marker = true;
					} else if (conjunct->type == hol_term_type::EQUALS) {
						if (!has_intersection<built_in_predicates>(expected_left_side, conjunct->binary.left)) {
							free(*expected_left_side); free(expected_left_side);
							return (hol_term*) nullptr;
						}
						free(*expected_left_side); free(expected_left_side);

						hol_term* arg = conjunct->binary.right;
						if (arg->type == hol_term_type::ANY) {
							unsigned int set_variable = lambda_variable + 1;
							unsigned int element_variable = lambda_variable + 2;
							dst = make_selected_arg_without_head(negation_count, element_variable, set_variable, lambda_variable);
							if (dst == nullptr) return (hol_term*) nullptr;
							if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
								hol_term* new_dst = hol_term::new_any_right(dst);
								if (dst == nullptr) {
									free(*dst); if (dst->reference_count == 0) free(dst);
									return (hol_term*) nullptr;
								}
								dst = new_dst;
							}
							remove_wide_scope_marker = true;
						} else if (arg->type == hol_term_type::CONSTANT || arg->type == hol_term_type::ANY_CONSTANT || arg->type == hol_term_type::ANY_CONSTANT_EXCEPT) {
							unsigned int set_variable = lambda_variable + 1;
							unsigned int element_variable = lambda_variable + 2;
							hol_term* element_var = hol_term::new_variable(element_variable);
							if (element_var == nullptr)
								return (hol_term*) nullptr;
							hol_term* set_var = hol_term::new_variable(set_variable);
							if (set_var == nullptr) {
								free(*element_var); free(element_var);
								return (hol_term*) nullptr;
							}
							hol_term* quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
									hol_term::new_apply(set_var, element_var),
									hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)
								));
							if (quantified_term == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;

							for (unsigned int i = 0; i < negation_count; i++) {
								hol_term* temp = hol_term::new_not(quantified_term);
								if (temp == nullptr) {
									free(*quantified_term); free(quantified_term);
									return (hol_term*) nullptr;
								}
								quantified_term = temp;
							}

							hol_term* set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, hol_term::new_equals(element_var, arg)));
							if (set_definition == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
							set_var->reference_count++;
							arg->reference_count++;
							dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term));
							if (dst == nullptr) {
								free(*set_definition); free(set_definition);
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							HOL_ANY.reference_count++;
						} else if (arg->type == hol_term_type::VARIABLE) {
							unsigned int element_variable = arg->variable;
							unsigned int index = scopes.index_of(element_variable);
							scope_info& scope = scopes.values[index];
							hol_term* operand = scope.scope->quantifier.operand;

							bool wide_scope_must_be_before_quantifier = false;
							bool wide_scope_must_be_after_quantifier = false;
							bool wide_scope_can_be_before_quantifier = false;
							if (narrow_scope && !scope.narrow_scope) {
								wide_scope_must_be_after_quantifier = true;
								wide_scope_can_be_before_quantifier = false;
							} else if (narrow_scope && scope.narrow_scope) {
								wide_scope_must_be_before_quantifier = true;
							} else {
								wide_scope_can_be_before_quantifier = false;
								for (unsigned int i = 0; i < index; i++) {
									if (scopes.keys[i] == 0) {
										wide_scope_can_be_before_quantifier = true;
										break;
									}
								}
							}

							if (scope.scope->type == hol_term_type::FOR_ALL && operand->type == hol_term_type::IF_THEN) {
								/* check if the scope is a relativized quantification over a set */
								hol_term* antecedent = operand->binary.left;
								if (antecedent->type == hol_term_type::UNARY_APPLICATION && antecedent->binary.right->type == hol_term_type::VARIABLE
								 && antecedent->binary.right->variable == element_variable && antecedent->binary.left->type == hol_term_type::VARIABLE)
								{
									hol_term* new_head = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (new_head == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;

									unsigned int index = scopes.index_of(antecedent->binary.left->variable);
									scope_info& set_scope = scopes.values[index];

									auto check_wide_scope = [&remove_wide_scope_marker](hol_term* term) {
										remove_wide_scope_marker |= (term->type == hol_term_type::UNARY_APPLICATION
											&& term->binary.left->type == hol_term_type::CONSTANT
											&& term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
										return true;
									};
									array<hol_term*> child_siblings(8); array<unsigned int> dst_variables(8);
									bool removed_quantifier, dummy = false, remove_negations;
									bool result = apply_head(set_scope.scope, dst, dst_variables, child_siblings, 0, removed_quantifier, remove_negations, dummy, node_finder(is_array ? parent : head),
									[new_head](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings) { return new_head; }, check_wide_scope);
									if (!result || dst == nullptr) {
										free(*new_head); free(new_head);
										return (hol_term*) nullptr;
									}
									for (hol_term* sibling : child_siblings) {
										unsigned int index = siblings.index_of(sibling);
										if (index < siblings.length) siblings.remove(index);
									}
								} else {
									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr)
										return (hol_term*) nullptr;
									hol_term* set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, antecedent));
									if (set_definition == nullptr) {
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									antecedent->reference_count++;
									hol_term* quantified_term;
									if (wide_scope_must_be_after_quantifier) {
										quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
										if (quantified_term == nullptr) {
											free(*set_definition); free(set_definition);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;
										remove_wide_scope_marker = true;

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition); free(set_definition);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									} else if (wide_scope_must_be_before_quantifier || !wide_scope_can_be_before_quantifier) {
										quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)));
										if (quantified_term == nullptr) {
											free(*set_definition); free(set_definition);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;

										if (scope.scope == next_universal_scope) {
											hol_term* temp = hol_term::new_apply(
													hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), quantified_term);
											remove_wide_scope_marker = true;
											if (temp == nullptr) {
												free(*set_definition); free(set_definition);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition); free(set_definition);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									} else {
										quantified_term = hol_term::new_any_right(hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
										if (quantified_term == nullptr) {
											free(*set_definition); free(set_definition);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition); free(set_definition);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									}
									dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term));
									if (dst == nullptr) {
										free(*set_definition); free(set_definition);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									HOL_ANY.reference_count++;
								}
							} else if (scope.scope->type == hol_term_type::EXISTS && operand->type == hol_term_type::AND) {
								hol_term* first_conjunct = operand->array.operands[0];
								if (first_conjunct->type == hol_term_type::UNARY_APPLICATION && first_conjunct->binary.right->type == hol_term_type::VARIABLE
								 && first_conjunct->binary.right->variable == element_variable && first_conjunct->binary.left->type == hol_term_type::VARIABLE)
								{
									/* check if the scope is a relativized quantification over a set */
									hol_term* new_head = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (new_head == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;

									unsigned int index = scopes.index_of(first_conjunct->binary.left->variable);
									scope_info& set_scope = scopes.values[index];

									auto check_wide_scope = [&remove_wide_scope_marker](hol_term* term) {
										remove_wide_scope_marker |= (term->type == hol_term_type::UNARY_APPLICATION
											&& term->binary.left->type == hol_term_type::CONSTANT
											&& term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
										return true;
									};
									array<hol_term*> child_siblings(8); array<unsigned int> dst_variables(8);
									bool removed_quantifier, dummy = false, remove_negations;
									bool result = apply_head(set_scope.scope, dst, dst_variables, child_siblings, 0, removed_quantifier, remove_negations, dummy, node_finder(is_array ? parent : head),
									[new_head](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings) { return new_head; }, check_wide_scope);
									if (!result || dst == nullptr) {
										free(*new_head); free(new_head);
										return (hol_term*) nullptr;
									}
									for (hol_term* sibling : child_siblings) {
										unsigned int index = siblings.index_of(sibling);
										if (index < siblings.length) siblings.remove(index);
									}
								} else {
									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr)
										return (hol_term*) nullptr;
									hol_term* set_definition;
									if (operand->array.length == 2) {
										set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, operand->array.operands[0]));
									} else {
										set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable,
												hol_term::new_and(make_array_view(operand->array.operands, operand->array.length - 1))));
									}
									if (set_definition == nullptr) {
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									for (unsigned int i = 0; i + 1 < operand->array.length; i++)
										operand->array.operands[i]->reference_count++;
									hol_term* quantified_term;
									if (wide_scope_must_be_after_quantifier) {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									} else if (wide_scope_must_be_before_quantifier || !wide_scope_can_be_before_quantifier) {
										quantified_term = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE),
												hol_term::new_exists(element_variable, hol_term::new_and(
													hol_term::new_apply(set_var, arg),
													hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									} else {
										quantified_term = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									}
									if (quantified_term == nullptr) {
										free(*set_definition); free(set_definition);
										return (hol_term*) nullptr;
									}
									set_var->reference_count++;
									arg->reference_count += 2;

									for (unsigned int i = 0; i < scope.negations; i++) {
										hol_term* temp = hol_term::new_not(quantified_term);
										if (temp == nullptr) {
											free(*set_definition); free(set_definition);
											free(*quantified_term); free(quantified_term);
											return (hol_term*) nullptr;
										}
										quantified_term = temp;
									}

									dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term));
									if (dst == nullptr) {
										free(*set_definition); free(set_definition);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									remove_wide_scope_marker = true;
								}
							} else if (scope.scope->type == hol_term_type::EXISTS && operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::AND || operand->any_array.oper == hol_term_type::ANY_ARRAY)) {
								hol_term* first_conjunct = (operand->any_array.left.length == 0 ? nullptr : operand->any_array.left.operands[0]);
								/* check if the scope is a relativized quantification over a set */
								if (first_conjunct != nullptr && first_conjunct->type == hol_term_type::UNARY_APPLICATION && first_conjunct->binary.right->type == hol_term_type::VARIABLE
								 && first_conjunct->binary.right->variable == element_variable && first_conjunct->binary.left->type == hol_term_type::VARIABLE)
								{
									hol_term* new_head = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (new_head == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;

									unsigned int index = scopes.index_of(first_conjunct->binary.left->variable);
									scope_info& set_scope = scopes.values[index];

									auto check_wide_scope = [&remove_wide_scope_marker](hol_term* term) {
										remove_wide_scope_marker |= (term->type == hol_term_type::UNARY_APPLICATION
											&& term->binary.left->type == hol_term_type::CONSTANT
											&& term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
										return true;
									};
									array<hol_term*> child_siblings(8); array<unsigned int> dst_variables(8);
									bool removed_quantifier, dummy = false, remove_negations;
									bool result = apply_head(set_scope.scope, dst, dst_variables, child_siblings, 0, removed_quantifier, remove_negations, dummy, node_finder(is_array ? parent : head),
									[new_head](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings) { return new_head; }, check_wide_scope);
									if (!result || dst == nullptr) {
										free(*new_head); free(new_head);
										return (hol_term*) nullptr;
									}
									for (hol_term* sibling : child_siblings) {
										unsigned int index = siblings.index_of(sibling);
										if (index < siblings.length) siblings.remove(index);
									}
								} else {
									hol_term* set_definition = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
											make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
											make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
											make_array_view(operand->any_array.right.operands, operand->any_array.right.length - 1));
									if (set_definition == nullptr)
										return (hol_term*) nullptr;
									operand->any_array.all->reference_count++;
									for (unsigned int i = 0; i < set_definition->any_array.any.length; i++)
										set_definition->any_array.any.operands[i]->reference_count++;
									for (unsigned int i = 0; i < set_definition->any_array.left.length; i++)
										set_definition->any_array.left.operands[i]->reference_count++;
									for (unsigned int i = 0; i < set_definition->any_array.right.length; i++)
										set_definition->any_array.right.operands[i]->reference_count++;

									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr) {
										free(*set_definition); free(set_definition);
										return (hol_term*) nullptr;
									}
									hol_term* set_definition_term = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, set_definition));
									if (set_definition_term == nullptr) {
										free(*set_definition); free(set_definition);
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									hol_term* quantified_term;
									if (wide_scope_must_be_after_quantifier) {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
										if (quantified_term == nullptr) {
											free(*set_definition_term); free(set_definition_term);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition_term); free(set_definition_term);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									} else if (wide_scope_can_be_before_quantifier) {
										quantified_term = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
										if (quantified_term == nullptr) {
											free(*set_definition_term); free(set_definition_term);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition_term); free(set_definition_term);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									} else {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)));
										if (quantified_term == nullptr) {
											free(*set_definition_term); free(set_definition_term);
											return (hol_term*) nullptr;
										}
										set_var->reference_count++;
										arg->reference_count += 2;

										for (unsigned int i = 0; i < scope.negations; i++) {
											hol_term* temp = hol_term::new_not(quantified_term);
											if (temp == nullptr) {
												free(*set_definition_term); free(set_definition_term);
												free(*quantified_term); free(quantified_term);
												return (hol_term*) nullptr;
											}
											quantified_term = temp;
										}
									}
									dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
											make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition_term, 1), make_array_view(&quantified_term, 1)));
									if (dst == nullptr) {
										free(*set_definition_term); free(set_definition_term);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									HOL_ANY.reference_count++;
									remove_wide_scope_marker = true;
								}
							} else if (scope.scope->type == hol_term_type::LAMBDA) {
								dst = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
								if (dst == nullptr)
									return (hol_term*) nullptr;
								arg->reference_count++;
							}
						}
					} else if (conjunct->type == hol_term_type::EXISTS) {
						hol_term* inner_operand = conjunct->quantifier.operand;
						hol_term* equals_term = hol_term::new_equals(expected_left_side, hol_term::new_variable(conjunct->quantifier.variable));
						if (equals_term == nullptr) {
							free(*expected_left_side); free(expected_left_side);
							return (hol_term*) nullptr;
						}
						hol_term* expected_operand = hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&equals_term, 1));
						if (expected_operand == nullptr) {
							free(*equals_term); free(equals_term);
							return (hol_term*) nullptr;
						}
						equals_term->reference_count++;
						HOL_ANY.reference_count++;
						if (!has_intersection<built_in_predicates>(expected_operand, inner_operand)) {
							free(*expected_operand); free(expected_operand);
							free(*equals_term); free(equals_term);
							return (hol_term*) nullptr;
						}
						free(*expected_operand); free(expected_operand);

						unsigned int set_variable = lambda_variable + 1;
						unsigned int element_variable = conjunct->quantifier.variable;
						hol_term* element_var = hol_term::new_variable(conjunct->quantifier.variable);
						if (element_var == nullptr) {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
							return (hol_term*) nullptr;
						}
						hol_term* set_var = hol_term::new_variable(set_variable);
						if (set_var == nullptr) {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
							free(*element_var); free(element_var);
							return (hol_term*) nullptr;
						}
						if (inner_operand->type == hol_term_type::ANY) {
							array<hol_term*> diff(1);
							subtract_any<built_in_predicates>(diff, inner_operand, equals_term);
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
#if !defined(NDEBUG)
							if (diff.length > 1)
								fprintf(stderr, "select_arg_without_head_predicative WARNING: `diff` has more than one logical form.\n");
#endif

							hol_term* quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)));
							if (quantified_term == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								free_all(diff); return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;

							for (unsigned int i = 0; i < negation_count; i++) {
								hol_term* temp = hol_term::new_not(quantified_term);
								if (temp == nullptr) {
									free(*quantified_term); free(quantified_term);
									free_all(diff); return (hol_term*) nullptr;
								}
								quantified_term = temp;
							}

							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, diff[0])), quantified_term));
							for (hol_term* term : diff) { free(*term); if (term->reference_count == 0) free(term); }
							if (dst == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							set_var->reference_count++;
						} else if (inner_operand->type == hol_term_type::ANY_ARRAY) {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
							hol_term* new_inner_operand = hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all,
									make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
									make_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length),
									make_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length - 1));
							if (new_inner_operand == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							new_inner_operand->any_array.all->reference_count++;
							for (unsigned int i = 0; i < new_inner_operand->any_array.any.length; i++)
								new_inner_operand->any_array.any.operands[i]->reference_count++;
							for (unsigned int i = 0; i < new_inner_operand->any_array.left.length; i++)
								new_inner_operand->any_array.left.operands[i]->reference_count++;
							for (unsigned int i = 0; i < new_inner_operand->any_array.right.length; i++)
								new_inner_operand->any_array.right.operands[i]->reference_count++;

							hol_term* quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)));
							if (quantified_term == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;

							for (unsigned int i = 0; i < negation_count; i++) {
								hol_term* temp = hol_term::new_not(quantified_term);
								if (temp == nullptr) {
									free(*quantified_term); free(quantified_term);
									return (hol_term*) nullptr;
								}
								quantified_term = temp;
							}

							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, new_inner_operand)), quantified_term));
							if (dst == nullptr) {
								free(*quantified_term); free(quantified_term);
								free(*new_inner_operand); free(new_inner_operand);
								return (hol_term*) nullptr;
							}
							set_var->reference_count++;
						} else if (inner_operand->type == hol_term_type::AND) {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
							hol_term* set_definition;
							if (inner_operand->array.length == 2) {
								set_definition = inner_operand->array.operands[0];
							} else {
								set_definition = hol_term::new_and(make_array_view(inner_operand->array.operands, inner_operand->array.length - 1));
							}

							hol_term* quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)));
							if (quantified_term == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;

							for (unsigned int i = 0; i < negation_count; i++) {
								hol_term* temp = hol_term::new_not(quantified_term);
								if (temp == nullptr) {
									free(*quantified_term); free(quantified_term);
									return (hol_term*) nullptr;
								}
								quantified_term = temp;
							}

							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, set_definition)), quantified_term));
							if (dst == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							set_var->reference_count++;
							for (unsigned int i = 0; i < inner_operand->array.length - 1; i++)
								inner_operand->array.operands[i]->reference_count++;
						} else {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
						}
					} else {
						free(*expected_left_side); free(expected_left_side);
						return (hol_term*) nullptr;
					}
				}
				return dst;
			}, gather_scope) && dst != nullptr;

	if (!result) return false;
	hol_term* new_dst = hol_term::new_lambda(lambda_variable, dst);
	if (new_dst == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		return false;
	}
	dst = new_dst;
	return true;
}

template<unsigned int ArgConstant>
inline bool select_singleton_arg_in_set_without_head_predicative(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
		src = src->quantifier.operand;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				unsigned int set_variable, element_variable = 0;
				hol_term* left = nullptr;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT)
				{
					set_variable = ++max_variable;
					element_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS)
						fprintf(stderr, "select_singleton_arg_in_set_without_head_predicative WARNING: Expected an existential quantification.\n");
#endif
					set_variable = head->quantifier.variable;
					hol_term* operand = head->quantifier.operand;

					hol_term* right = nullptr;
					if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
						element_variable = ++max_variable;
					} else if (operand->type == hol_term_type::AND && operand->array.length != 0) {
						left = operand->array.operands[0];
						right = operand->array.operands[operand->array.length - 1];
					} else if (operand->type == hol_term_type::ANY_ARRAY) {
						if (operand->any_array.left.length != 0)
							left = operand->any_array.left.operands[0];
						else left = operand->any_array.all;
						if (operand->any_array.right.length != 0)
							right = operand->any_array.right.operands[operand->any_array.right.length - 1];
						else right = operand->any_array.all;
					} else {
						return (hol_term*) nullptr;
					}

					hol_term* set_definition = nullptr;
					if (left != nullptr) {
						if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
							if (left->any.included != nullptr)
								set_definition = left->any.included;
						} else if (left->type == hol_term_type::EQUALS) {
							set_definition = left->binary.right;
						} else if (left->type == hol_term_type::BINARY_APPLICATION) {
							set_definition = left->ternary.third;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (set_definition != nullptr) {
						if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
							/* no-op */
						} else if (set_definition->type == hol_term_type::LAMBDA) {
							element_variable = set_definition->quantifier.variable;
						} else {
							return (hol_term*) nullptr;
						}
					}

					if (element_variable == 0 && right != nullptr) {
						/* try to get the element variable from the right conjunct, since we couldn't get it from the left conjunct */
						if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
							right = right->any.included;
						if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
							right = right->binary.right;
						if (right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) {
							element_variable = ++max_variable;
						} else if (right->type == hol_term_type::EXISTS || right->type == hol_term_type::FOR_ALL) {
							element_variable = right->quantifier.variable;
						} else if (right->type == hol_term_type::AND && right->array.length == 2 && right->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
								&& right->array.operands[1]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[1]->binary.left->variable == set_variable
								&& right->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = right->array.operands[1]->binary.right->variable;
						} else {
							return (hol_term*) nullptr;
						}
					}
				}

				hol_term* head_var = hol_term::new_variable(element_variable);
				if (head_var == nullptr) return (hol_term*) nullptr;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return (hol_term*) nullptr;
				}
				free(*head_var);

				hol_term* singleton = hol_term::new_any_right(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), &HOL_ANY), head_var), excluded_trees, excluded_tree_count);
				if (singleton == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;
				head_var->reference_count++;

				hol_term* set_definition;
				if (left == nullptr || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
					set_definition = hol_term::new_any_right(hol_term::new_lambda(element_variable, singleton));
					if (set_definition == nullptr) {
						free(*singleton); free(singleton);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::EQUALS) {
					set_definition = hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, singleton));
					if (set_definition == nullptr) {
						free(*singleton); free(singleton);
						return (hol_term*) nullptr;
					}
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = hol_term::new_apply(&HOL_ANY, hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, singleton));
					if (set_definition == nullptr) {
						free(*singleton); free(singleton);
						return (hol_term*) nullptr;
					}
					HOL_ANY.reference_count++;
				}

				hol_term* lambda_var = hol_term::new_variable(lambda_variable);
				if (lambda_var == nullptr) {
					free(*set_definition); free(set_definition);
					return (hol_term*) nullptr;
				}

				hol_term* set_var = hol_term::new_variable(set_variable);
				if (set_var == nullptr) {
					free(*set_definition); free(set_definition);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}

				hol_term* element_var = hol_term::new_variable(element_variable);
				if (element_var == nullptr) {
					free(*set_definition); free(set_definition);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}

				hol_term* dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, hol_term::new_exists(element_variable, hol_term::new_and(
						hol_term::new_apply(set_var, element_var), hol_term::new_apply(lambda_var, element_var)))));
				if (dst == nullptr) {
					free(*set_definition); free(set_definition);
					free(*element_var); free(element_var);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}
				element_var->reference_count += 2 - 1;
				set_var->reference_count++;
				lambda_var->reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, dst, head);
				free(*dst); if (dst->reference_count == 0) free(dst);
				if (intersection.length == 0) {
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "select_singleton_arg_in_set_without_head_predicative ERROR: Expected intersection to be unique.\n");
					free_all(intersection);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}

				hol_term* operand = intersection[0]->quantifier.operand;
				hol_term* new_left;
				if (operand->type == hol_term_type::AND) {
					new_left = operand->array.operands[0];
				} else {
					new_left = operand->any_array.left.operands[0];
				}

				hol_term* inner_operand;
				if (new_left->type == hol_term_type::ANY || new_left->type == hol_term_type::ANY_RIGHT) {
					inner_operand = new_left->any.included->quantifier.operand;
				} else if (new_left->type == hol_term_type::EQUALS) {
					inner_operand = new_left->binary.right->quantifier.operand;
				} else {
					inner_operand = new_left->ternary.third->quantifier.operand;
				}

				hol_term* second_inner_operand;
				if (inner_operand->type == hol_term_type::EXISTS) {
					second_inner_operand = inner_operand->quantifier.operand;
				} else if (inner_operand->type == hol_term_type::ANY_RIGHT) {
					free_all(intersection);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);

					unsigned int set_variable = lambda_variable + 1;
					unsigned int element_variable = lambda_variable + 2;
					hol_term* new_term = make_selected_arg_without_head(0, element_variable, set_variable, lambda_variable);
					if (new_term == nullptr) return (hol_term*) nullptr;
					remove_wide_scope_marker = true;
					return new_term;
				} else {
					free_all(intersection);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}

				hol_term* last_conjunct;
				if (second_inner_operand->type == hol_term_type::AND) {
					last_conjunct = second_inner_operand->array.operands[second_inner_operand->array.length - 1];
				} else if (second_inner_operand->type == hol_term_type::ANY_ARRAY && second_inner_operand->any_array.oper == hol_term_type::AND) {
					if (second_inner_operand->any_array.right.length != 0) {
						last_conjunct = second_inner_operand->any_array.right.operands[second_inner_operand->any_array.right.length - 1];
					} else {
						last_conjunct = second_inner_operand->any_array.all;
					}
				} else {
					free_all(intersection);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}

				if (last_conjunct->type != hol_term_type::EQUALS
				 || last_conjunct->binary.right->type != hol_term_type::VARIABLE
				 || last_conjunct->binary.right->variable != element_variable
				 || last_conjunct->binary.left->type != hol_term_type::UNARY_APPLICATION
				 || last_conjunct->binary.left->binary.left->type != hol_term_type::CONSTANT
				 || last_conjunct->binary.left->binary.right->type != hol_term_type::VARIABLE
				 || last_conjunct->binary.left->binary.right->variable != inner_operand->quantifier.variable)
				{
					free_all(intersection);
					free(*set_var); free(set_var);
					free(*lambda_var); free(lambda_var);
					return (hol_term*) nullptr;
				}
				unsigned int new_element_variable = inner_operand->quantifier.variable;

				hol_term* new_second_inner_operand;
				if (second_inner_operand->type == hol_term_type::AND) {
					if (second_inner_operand->array.length == 2) {
						new_second_inner_operand = second_inner_operand->array.operands[0];
						new_second_inner_operand->reference_count++;
					} else {
						new_second_inner_operand = hol_term::new_and(make_array_view(second_inner_operand->array.operands, second_inner_operand->array.length - 1));
						if (new_second_inner_operand == nullptr) {
							free_all(intersection);
							free(*set_var); free(set_var);
							free(*lambda_var); free(lambda_var);
							return (hol_term*) nullptr;
						}
					}
				} else {
					if (second_inner_operand->any_array.right.length != 0) {
						new_second_inner_operand = hol_term::new_any_array(second_inner_operand->any_array.oper, second_inner_operand->any_array.all,
								make_array_view(second_inner_operand->any_array.any.operands, second_inner_operand->any_array.any.length),
								make_array_view(second_inner_operand->any_array.left.operands, second_inner_operand->any_array.left.length),
								make_array_view(second_inner_operand->any_array.right.operands, second_inner_operand->any_array.right.length - 1));
						if (new_second_inner_operand == nullptr) {
							free_all(intersection);
							free(*set_var); free(set_var);
							free(*lambda_var); free(lambda_var);
							return(hol_term*) nullptr;
						}
						new_second_inner_operand->any_array.all->reference_count++;
						for (unsigned int i = 0; i < new_second_inner_operand->any_array.any.length; i++)
							new_second_inner_operand->any_array.any.operands[i]->reference_count++;
						for (unsigned int i = 0; i < new_second_inner_operand->any_array.left.length; i++)
							new_second_inner_operand->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < new_second_inner_operand->any_array.right.length; i++)
							new_second_inner_operand->any_array.right.operands[i]->reference_count++;
					} else {
						new_second_inner_operand = second_inner_operand;
						second_inner_operand->reference_count++;
					}
				}
				free_all(intersection);

				hol_term* new_element_var = hol_term::new_variable(new_element_variable);
				if (new_element_var == nullptr) {
					free(*new_second_inner_operand); if (new_second_inner_operand->reference_count == 0) free(new_second_inner_operand);
					free(*set_var); if (set_var->reference_count == 0) free(set_var);
					free(*lambda_var); if (lambda_var->reference_count == 0) free(lambda_var);
					return (hol_term*) nullptr;
				}

				dst = hol_term::new_exists(set_variable, hol_term::new_and(
						hol_term::new_equals(set_var, hol_term::new_lambda(new_element_variable, new_second_inner_operand)),
						hol_term::new_exists(new_element_variable, hol_term::new_and(
								hol_term::new_apply(set_var, new_element_var),
								hol_term::new_apply(lambda_var, new_element_var)
						))
					));
				if (dst == nullptr) {
					free(*new_second_inner_operand); if (new_second_inner_operand->reference_count == 0) free(new_second_inner_operand);
					free(*set_var); if (set_var->reference_count == 0) free(set_var);
					free(*lambda_var); if (lambda_var->reference_count == 0) free(lambda_var);
					return (hol_term*) nullptr;
				}
				new_element_var->reference_count += 2 - 1;
				set_var->reference_count += 2 - 1;

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline hol_term* select_arg_without_head_process_operand(hol_term* head, hol_term* operand)
{
	if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
		HOL_ANY.reference_count++;
		return &HOL_ANY;
	} else if (operand->type == hol_term_type::ANY_QUANTIFIER) {
		/* TODO: this isn't really the result of applying this function to `src` */
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		unsigned int head_variable = ++max_variable;
		hol_term* dst = hol_term::new_exists(head_variable, &HOL_ANY);
		if (dst == nullptr) return nullptr;
		HOL_ANY.reference_count++;
		return dst;
	} else if (operand->type == hol_term_type::EXISTS) {
		hol_term* inner_operand = operand->quantifier.operand;
		unsigned int head_variable = operand->quantifier.variable;

		hol_term* dst;
		if (inner_operand->type == hol_term_type::ANY_ARRAY) {
			dst = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all,
					make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
					make_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length),
					make_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length - 1)));
			if (dst == nullptr) {
				return nullptr;					
			}
			inner_operand->any_array.all->reference_count++;
			for (unsigned int i = 0; i < inner_operand->any_array.any.length; i++)
				inner_operand->any_array.any.operands[i]->reference_count++;
			for (unsigned int i = 0; i < inner_operand->any_array.left.length; i++)
				inner_operand->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i < inner_operand->any_array.right.length - 1; i++)
				inner_operand->any_array.right.operands[i]->reference_count++;
		} else if (inner_operand->type == hol_term_type::AND) {
			if (inner_operand->array.length == 2) {
				dst = hol_term::new_exists(head_variable, inner_operand->array.operands[0]);
			} else {
				dst = hol_term::new_exists(head_variable, hol_term::new_and(make_array_view(inner_operand->array.operands, inner_operand->array.length - 1)));
			}
			if (dst == nullptr) {
				return nullptr;					
			}
			for (unsigned int i = 0; i < inner_operand->array.length - 1; i++)
				inner_operand->array.operands[i]->reference_count++;
		} else {
			return nullptr;
		}
		return dst;
	} else {
		return nullptr;
	}
}

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant, bool InvertArg>
inline bool select_arg_without_head(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates, true>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				if (head->type == hol_term_type::NOT)
					head = head->unary.operand;

				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					/* TODO: this isn't really the result of applying this function to `src` */
					unsigned int max_variable = 0;
					max_bound_variable(*head, max_variable);
					unsigned int head_variable = ++max_variable;
					hol_term* dst = hol_term::new_exists(head_variable, &HOL_ANY);
					if (dst == nullptr) return (hol_term*) nullptr;
					HOL_ANY.reference_count++;
					return dst;
				} else {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunct = nullptr;
					if (operand->type == hol_term_type::ANY_ARRAY) {
						if (ConjunctIndex >= 0) {
							if (ConjunctIndex < operand->any_array.left.length)
								conjunct = operand->any_array.left.operands[ConjunctIndex];
							else conjunct = operand->any_array.all;
						} else if (ConjunctIndex < 0) {
							unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
							if (index < operand->any_array.right.length)
								conjunct = operand->any_array.right.operands[operand->any_array.right.length - index - 1];
							else conjunct = operand->any_array.all;
						}
					} else if (operand->type == hol_term_type::AND) {
						int conjunct_index = ConjunctIndex;
						if (ConjunctIndex < 0)
							conjunct_index += operand->array.length;
						conjunct = operand->array.operands[conjunct_index];
					}
#if !defined(NDEBUG)
					else fprintf(stderr, "select_arg_without_head ERROR: Expected the head to be an existentially-quantified conjunction.\n");
#endif

					hol_term* expected_conjunct;
					if (InvertArg) {
						expected_conjunct = hol_term::new_any_right(hol_term::new_equals(
								hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(head->quantifier.variable)),
								&HOL_ANY));
					} else {
						expected_conjunct = hol_term::new_any_right(hol_term::new_equals(
								hol_term::new_apply(hol_term::new_constant(ArgConstant), &HOL_ANY),
								hol_term::new_variable(head->quantifier.variable)));
					}
					if (expected_conjunct == nullptr)
						return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, expected_conjunct, conjunct);
					free(*expected_conjunct); if (expected_conjunct->reference_count == 0) free(expected_conjunct);
					if (intersection.length == 0) {
						return (hol_term*) nullptr;
					} else if (intersection.length != 1) {
						fprintf(stderr, "select_arg_without_head ERROR: Intersection is not unique.\n");
						free_all(intersection);
						return (hol_term*) nullptr;
					}

					if (intersection[0]->type == hol_term_type::ANY || intersection[0]->type == hol_term_type::ANY_RIGHT) {
						free_all(intersection);
						HOL_ANY.reference_count++;
						return &HOL_ANY;

					} else if (intersection[0]->type == hol_term_type::ANY_ARRAY) {
						hol_term* all = select_arg_without_head_process_operand(head, intersection[0]->any_array.all);

						array<hol_term*> any(max(1u, intersection[0]->any_array.any.length));
						for (unsigned int i = 0; i < intersection[0]->any_array.any.length; i++) {
							any[i] = select_arg_without_head_process_operand(head, intersection[0]->any_array.any.operands[i]);
							if (any[i] == nullptr) {
								free(*all); if (all->reference_count == 0) free(all);
								free_all(intersection); free_all(any);
								return (hol_term*) nullptr;
							}
							any.length++;
						}

						array<hol_term*> left(max(1u, intersection[0]->any_array.left.length));
						for (unsigned int i = 0; i < intersection[0]->any_array.left.length; i++) {
							left[i] = select_arg_without_head_process_operand(head, intersection[0]->any_array.left.operands[i]);
							if (left[i] == nullptr) {
								free(*all); if (all->reference_count == 0) free(all);
								free_all(intersection); free_all(any); free_all(left);
								return (hol_term*) nullptr;
							}
							left.length++;
						}

						array<hol_term*> right(max(1u, intersection[0]->any_array.right.length));
						for (unsigned int i = 0; i < intersection[0]->any_array.right.length; i++) {
							right[i] = select_arg_without_head_process_operand(head, intersection[0]->any_array.right.operands[i]);
							if (right[i] == nullptr) {
								free(*all); if (all->reference_count == 0) free(all);
								free_all(intersection); free_all(any); free_all(left); free_all(right);
								return (hol_term*) nullptr;
							}
							right.length++;
						}
						free_all(intersection);

						hol_term* dst = hol_term::new_any_array(intersection[0]->any_array.oper, all,
								make_array_view(any.data, any.length), make_array_view(left.data, left.length), make_array_view(right.data, right.length));
						if (dst == nullptr) {
							free(*all); if (all->reference_count == 0) free(all);
							free_all(any); free_all(left); free_all(right);
							return (hol_term*) nullptr;
						}
						return dst;

					} else if (intersection[0]->type == hol_term_type::AND || intersection[0]->type == hol_term_type::OR) {
						array<hol_term*> new_operands(intersection[0]->array.length);
						for (unsigned int i = 0; i < intersection[0]->array.length; i++) {
							new_operands[i] = select_arg_without_head_process_operand(head, intersection[0]->array.operands[i]);
							if (new_operands[i] == nullptr) {
								free_all(intersection); free_all(new_operands);
								return (hol_term*) nullptr;
							}
							new_operands.length++;
						}
						hol_term_type type = intersection[0]->type;
						free_all(intersection);

						hol_term* dst;
						if (type == hol_term_type::AND)
							dst = hol_term::new_and(make_array_view(new_operands.data, new_operands.length));
						else dst = hol_term::new_or(make_array_view(new_operands.data, new_operands.length));
						if (dst == nullptr) {
							free_all(new_operands);
							return (hol_term*) nullptr;
						}
						return dst;

					} else {
						hol_term* dst = select_arg_without_head_process_operand(head, intersection[0]);
						free_all(intersection);
						return dst;
					}
				}
			}, no_op()) && dst != nullptr;
}

template<size_t N>
inline unsigned int map_tense_constant(
		unsigned int constant,
		const unsigned int (&input_constants)[N],
		const unsigned int (&output_constants)[N])
{
	for (unsigned int i = 0; i < N; i++)
		if (input_constants[i] == constant) return output_constants[i];
	fprintf(stderr, "map_tense_constant ERROR: Unexpected input constant.\n");
	return 0;
}

template<bool UniqueOutput, size_t N>
inline bool map_tense_conjunct(
		array<hol_term*>& out,
		hol_term* src_conjunct, hol_term* input_conjunct,
		const unsigned int (&input_constants)[N],
		const unsigned int (&output_constants)[N])
{
	array<hol_term*> intersections(4);
	intersect<built_in_predicates>(intersections, src_conjunct, input_conjunct);
	if (UniqueOutput && intersections.length > 1) {
		fprintf(stderr, "map_tense_conjunct WARNING: Expected intersection size to be 1.\n");
		for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	} else if (intersections.length == 0) {
		return true;
	}

	if (out.ensure_capacity(out.length + intersections.length)) {
		for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	for (hol_term* intersection : intersections) {
		if (intersection->binary.left->type == hol_term_type::CONSTANT) {
			unsigned int new_constant = map_tense_constant(intersection->binary.left->constant, input_constants, output_constants);
			out[out.length] = hol_term::new_apply(hol_term::new_constant(new_constant), intersection->binary.right);
		} else {
#if !defined(NDEBUG)
			if (intersection->binary.left->type != hol_term_type::ANY_CONSTANT)
				fprintf(stderr, "map_tense_conjunct WARNING: Unexpected hol_term_type.\n");
#endif

			array<unsigned int> new_constants(intersection->binary.left->any_constant.length);
			for (unsigned int i = 0; i < intersection->binary.left->any_constant.length; i++)
				new_constants[new_constants.length++] = map_tense_constant(intersection->binary.left->any_constant.constants[i], input_constants, output_constants);
			if (new_constants.length > 1) {
				insertion_sort(new_constants);
				unique(new_constants);
			}

			if (new_constants.length == 1) {
				out[out.length] = hol_term::new_apply(hol_term::new_constant(new_constants[0]), intersection->binary.right);
			} else {
				out[out.length] = hol_term::new_apply(hol_term::new_any_constant(make_array_view(new_constants.data, new_constants.length)), intersection->binary.right);
			}
		}

		if (out[out.length] == nullptr) {
			for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		out.length++;
	}
	for (hol_term* term : intersections) { free(*term); if (term->reference_count == 0) free(term); }
	return true;
}

template<bool UniqueOutput, size_t N>
inline bool map_tense_predicate(array<hol_term*>& out, hol_term* head,
		unsigned int head_variable, head_index predicate_index,
		const unsigned int(&input_constants)[N], const unsigned int (&output_constants)[N])
{
	hol_term* input_conjunct = hol_term::new_apply(hol_term::new_any_constant(make_array_view(input_constants, N)), hol_term::new_variable(head_variable));
	if (input_conjunct == nullptr) return false;

	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
	 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
	 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
	{
		if (!out.ensure_capacity(out.length + 1))
			return false;
		out[out.length] = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, hol_term::new_any(nullptr, &input_conjunct, 1),
				make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
		if (out[out.length] == nullptr) {
			free(*input_conjunct); free(input_conjunct);
			return false;
		}
		out.length++;
	} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
		hol_term* operand = head->quantifier.operand;
		if ((operand->any_array.oper != hol_term_type::AND && operand->any_array.oper != hol_term_type::ANY_ARRAY)
		 || (predicate_index.position == head_position::RIGHT && predicate_index.index == 0))
		{
			free(*input_conjunct); free(input_conjunct);
			return false;
		}

		array<hol_term*> all_differences(4);
		subtract<built_in_predicates>(all_differences, operand->any_array.all, input_conjunct);
		if (UniqueOutput && all_differences.length > 1) {
			fprintf(stderr, "map_tense_predicate ERROR: Set difference is not unique.\n");
			for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
			free(*input_conjunct); free(input_conjunct);
			return false;
		} else if (all_differences.length == 0) {
			free(*input_conjunct); free(input_conjunct);
			return false;
		}

		unsigned int new_left_length = operand->any_array.left.length;
		if (predicate_index.position == head_position::LEFT && predicate_index.index == new_left_length - 1)
			new_left_length++;
		unsigned int new_any_length = operand->any_array.any.length;
		if (predicate_index.position == head_position::ANY && predicate_index.index == new_any_length - 1)
			new_any_length++;

		array<hol_term*>* left_differences = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * new_left_length));
		array<hol_term*>* right_differences = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * operand->any_array.right.length));
		array<hol_term*>* any_differences = (array<hol_term*>*) malloc(max((size_t) 1, sizeof(array<hol_term*>) * new_any_length));
		if (left_differences == nullptr || right_differences == nullptr || any_differences == nullptr) {
			if (left_differences != nullptr) free(left_differences);
			if (right_differences != nullptr) free(right_differences);
			for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
			free(*input_conjunct); free(input_conjunct);
			return false;
		}

		for (unsigned int i = 0; i < new_left_length; i++) {
			if (!array_init(left_differences[i], 4)) {
				for (unsigned int j = 0; j < i; j++) free(left_differences[j]);
				free(left_differences); free(right_differences); free(any_differences);
				for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*input_conjunct); free(input_conjunct);
				return false;
			}
		} for (unsigned int i = 0; i < operand->any_array.right.length; i++) {
			if (!array_init(right_differences[i], 4)) {
				for (unsigned int j = 0; j < new_left_length; j++) free(left_differences[j]);
				for (unsigned int j = 0; j < i; j++) free(right_differences[j]);
				free(left_differences); free(right_differences); free(any_differences);
				for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*input_conjunct); free(input_conjunct);
				return false;
			}
		} for (unsigned int i = 0; i < new_any_length; i++) {
			if (!array_init(any_differences[i], 4)) {
				for (unsigned int j = 0; j < new_left_length; j++) free(left_differences[j]);
				for (unsigned int j = 0; j < operand->any_array.right.length; j++) free(right_differences[j]);
				for (unsigned int j = 0; j < i; j++) free(any_differences[j]);
				free(left_differences); free(right_differences); free(any_differences);
				for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
				free(*input_conjunct); free(input_conjunct);
				return false;
			}
		}

		auto cleanup = [operand,new_left_length,new_any_length,left_differences,right_differences,any_differences,input_conjunct,&all_differences]() {
			for (unsigned int j = 0; j < new_left_length; j++) {
				for (hol_term* term : left_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(left_differences[j]);
			} for (unsigned int j = 0; j < operand->any_array.right.length; j++) {
				for (hol_term* term : right_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(right_differences[j]);
			} for (unsigned int j = 0; j < new_any_length; j++) {
				for (hol_term* term : any_differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(any_differences[j]);
			}
			free(left_differences); free(right_differences); free(any_differences);
			for (hol_term* term : all_differences) { free(*term); if (term->reference_count == 0) free(term); }
			free(*input_conjunct); free(input_conjunct);
		};

		for (unsigned int i = 0; i < operand->any_array.left.length; i++) {
			if (predicate_index.position == head_position::LEFT && i == predicate_index.index + 1) {
				map_tense_conjunct<UniqueOutput>(left_differences[i], operand->any_array.left.operands[i], input_conjunct, input_constants, output_constants);
			} else {
				subtract<built_in_predicates>(left_differences[i], operand->any_array.left.operands[i], input_conjunct);
				if (UniqueOutput && left_differences[i].length > 0) {
					fprintf(stderr, "map_tense_predicate ERROR: Set difference is not unique.\n");
					cleanup(); return false;
				}
			}
			if (left_differences[i].length == 0) {
				cleanup();
				return false;
			}
		}
		if (predicate_index.position == head_position::LEFT && operand->any_array.left.length == predicate_index.index + 1) {
			map_tense_conjunct<UniqueOutput>(left_differences[operand->any_array.left.length], operand->any_array.all, input_conjunct, input_constants, output_constants);
			if (left_differences[operand->any_array.left.length].length == 0) {
				cleanup();
				return false;
			}
		}

		for (unsigned int i = 0; i < operand->any_array.right.length; i++) {
			if (predicate_index.position == head_position::RIGHT && i == predicate_index.index - 1) {
				map_tense_conjunct<UniqueOutput>(right_differences[i], operand->any_array.right.operands[i], input_conjunct, input_constants, output_constants);
			} else {
				subtract<built_in_predicates>(right_differences[i], operand->any_array.right.operands[i], input_conjunct);
				if (UniqueOutput && right_differences[i].length > 0) {
					fprintf(stderr, "map_tense_predicate ERROR: Set difference is not unique.\n");
					cleanup(); return false;
				}
			}
			if (right_differences[i].length == 0) {
				cleanup();
				return false;
			}
		}

		for (unsigned int i = 0; i < operand->any_array.any.length; i++) {
			if (predicate_index.position == head_position::ANY && i == predicate_index.index + 1) {
				map_tense_conjunct<UniqueOutput>(any_differences[i], operand->any_array.any.operands[i], input_conjunct, input_constants, output_constants);
			} else {
				subtract<built_in_predicates>(any_differences[i], operand->any_array.any.operands[i], input_conjunct);
				if (UniqueOutput && any_differences[i].length > 0) {
					fprintf(stderr, "map_tense_predicate ERROR: Set difference is not unique.\n");
					cleanup(); return false;
				}
			}
			if (any_differences[i].length == 0) {
				cleanup();
				return false;
			}
		}
		if (predicate_index.position == head_position::LEFT && operand->any_array.any.length == predicate_index.index + 1) {
			map_tense_conjunct<UniqueOutput>(any_differences[operand->any_array.any.length], operand->any_array.all, input_conjunct, input_constants, output_constants);
			if (any_differences[operand->any_array.any.length].length == 0) {
				cleanup();
				return false;
			}
		}

		for (hol_term* all : all_differences) {
			bool result = apply_to_cartesian_product(left_differences, new_left_length, [&out,head,all,left_differences,right_differences,any_differences,operand,new_left_length,new_any_length](const unsigned int* left_index_array) {
				return apply_to_cartesian_product(right_differences, operand->any_array.right.length, [&out,head,all,left_differences,right_differences,any_differences,operand,new_left_length,new_any_length,left_index_array](const unsigned int* right_index_array) {
					return apply_to_cartesian_product(any_differences, new_any_length, [&out,head,all,left_differences,right_differences,any_differences,operand,new_left_length,new_any_length,left_index_array,right_index_array](const unsigned int* any_index_array) {
						if (!out.ensure_capacity(out.length + 1)) return false;
						out[out.length] = hol_term::new_exists(head->quantifier.variable, hol_term::new_any_array(hol_term_type::AND, all,
								lookup_table_array_view<array<hol_term*>, hol_term*>(any_differences, any_index_array, new_any_length),
								lookup_table_array_view<array<hol_term*>, hol_term*>(left_differences, left_index_array, new_left_length),
								lookup_table_array_view<array<hol_term*>, hol_term*>(right_differences, right_index_array, operand->any_array.right.length)));
						if (out[out.length] == nullptr)
							return false;

						out[out.length]->any_array.all->reference_count++;
						for (unsigned int i = 0; i < out[out.length]->any_array.left.length; i++)
							out[out.length]->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < out[out.length]->any_array.right.length; i++)
							out[out.length]->any_array.right.operands[i]->reference_count++;
						for (unsigned int i = 0; i < out[out.length]->any_array.any.length; i++)
							out[out.length]->any_array.any.operands[i]->reference_count++;
						out.length++;
						return true;
					});
				});
			});
			if (!result) {
				cleanup();
				return false;
			}
		}
		cleanup();
	} else {
#if !defined(NDEBUG)
		if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
			fprintf(stderr, "map_tense_predicate WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
		hol_term* operand = head->quantifier.operand;
		if (predicate_index.index + 1 == operand->array.length) {
			free(*input_conjunct); free(input_conjunct);
			return false;
		}

		array<hol_term*>* differences = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * operand->array.length);
		if (differences == nullptr) {
			free(*input_conjunct); free(input_conjunct);
			return false;
		}

		for (unsigned int i = 0; i < operand->array.length; i++) {
			if (!array_init(differences[i], 4)) {
				for (unsigned int j = 0; j < i; j++) free(differences[j]);
				free(differences); free(*input_conjunct); free(input_conjunct);
				return false;
			}
		}

		auto cleanup = [operand,differences,input_conjunct]() {
			for (unsigned int j = 0; j < operand->array.length; j++) {
				for (hol_term* term : differences[j]) { free(*term); if (term->reference_count == 0) free(term); }
				free(differences[j]);
			}
			free(differences); free(*input_conjunct); free(input_conjunct);
		};

		for (unsigned int i = 0; i < operand->array.length; i++) {
			if (i == predicate_index.index + 1) {
				map_tense_conjunct<UniqueOutput>(differences[i], operand->array.operands[i], input_conjunct, input_constants, output_constants);
			} else {
				subtract<built_in_predicates>(differences[i], operand->array.operands[i], input_conjunct);
				if (UniqueOutput && differences[i].length > 1) {
					fprintf(stderr, "map_tense_predicate ERROR: Set difference is not unique.\n");
					cleanup(); return false;
				}
			}
			if (differences[i].length == 0) {
				cleanup();
				return false;
			}
		}

		bool result = apply_to_cartesian_product(differences, operand->array.length, [&out,head,differences,operand](const unsigned int* index_array) {
			if (!out.ensure_capacity(out.length + 1)) return false;
			out[out.length] = hol_term::new_exists(head->quantifier.variable, hol_term::new_and(lookup_table_array_view<array<hol_term*>, hol_term*>(differences, index_array, operand->array.length)));
			if (out[out.length] == nullptr) return false;

			for (unsigned int i = 0; i < operand->array.length; i++)
				out[out.length]->array.operands[i]->reference_count++;
			out.length++;
			return true;
		});

		cleanup();
		if (!result) return false;
	}
	return true;
}

template<size_t N>
inline bool apply_tense_predicate(
		hol_term* src, hol_term*& dst,
		const unsigned int (&input_tense_predicates)[N],
		const unsigned int (&output_tense_predicates)[N])
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[input_tense_predicates,output_tense_predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				array<hol_term*> out(1);
				if (!map_tense_predicate<true>(out, head, head_variable, predicate_index, input_tense_predicates, output_tense_predicates) || out.length == 0)
					return (hol_term*) nullptr;
				hol_term* dst = out[0];

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<size_t N>
inline bool require_no_tense_predicates(
		hol_term* src, hol_term*& dst,
		const unsigned int (&tense_predicates)[N])
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[tense_predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* input_conjunct = hol_term::new_apply(hol_term::new_any_constant(make_array_view(tense_predicates, N)), hol_term::new_variable(head_variable));
				if (input_conjunct == nullptr)
					return (hol_term*) nullptr;

				hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, hol_term::new_any(nullptr, &input_conjunct, 1),
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*input_conjunct); free(input_conjunct);
					return (hol_term*) nullptr;
				}

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length > 1) {
					fprintf(stderr, "require_no_tense_predicates ERROR: Intersection is not unique.\n");
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return (hol_term*) nullptr;
				}

				hol_term* dst = intersection[0];
				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst, old_head->any.excluded_trees, old_head->any.excluded_tree_count);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < old_head->any.excluded_tree_count; i++)
						old_head->any.excluded_trees[i]->reference_count++;
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool remove_perfect(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				array<hol_term*> out(1);
				if (!map_tense_predicate<true>(out, head, head_variable, predicate_index, PERFECT_PREDICATES, NON_PERFECT_PREDICATES) || out.length == 0)
					return (hol_term*) nullptr;
				hol_term* dst = out[0];

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool remove_not(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if (head->type == hol_term_type::NOT) {
					head->reference_count++;
					return head;
				} else {
					return (hol_term*) nullptr;
				}
			}, no_op()) && dst != nullptr;
}

inline bool require_no_empty_ref(hol_term* src, hol_term*& dst)
{
	hol_term* hol_empty_ptr = &HOL_EMPTY;
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, &hol_empty_ptr, 1),
			hol_term::new_variable(0));
	if (predicate == nullptr) 
		return false;
	HOL_EMPTY.reference_count++;

	bool result = require_predicate<INT_FAST8_MAX>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool predicate_only(hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_apply(&HOL_ANY, hol_term::new_variable(head_variable)));
				if (expected_head == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, head, expected_head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length > 1) {
					fprintf(stderr, "predicate_only ERROR: Expected intersection size to be 1.\n");
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return (hol_term*) nullptr;
				} else if (intersection.length == 0) {
					return (hol_term*) nullptr;
				}

				hol_term* dst = intersection[0]->quantifier.operand->binary.left;
				dst->reference_count++;
				for (unsigned int i = 0; i < intersection.length; i++) {
					hol_term* term = intersection[i];
					free(*term); if (term->reference_count == 0) free(term);
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool predicate(hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					dst = &HOL_ANY;
					HOL_ANY.reference_count++;

				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || (head->quantifier.operand->type != hol_term_type::AND && head->quantifier.operand->type != hol_term_type::ANY_ARRAY))
						fprintf(stderr, "predicate WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif

					hol_term* operand = head->quantifier.operand;
					hol_term* predicate = nullptr;
					if (head->quantifier.operand->type == hol_term_type::AND) {
						predicate = operand->array.operands[predicate_index.index];
					} else if (predicate_index.position == head_position::LEFT) {
						predicate = operand->any_array.left.operands[predicate_index.index];
					} else if (predicate_index.position == head_position::RIGHT) {
						predicate = operand->any_array.right.operands[operand->any_array.right.length - predicate_index.index - 1];
					} else {
						predicate = operand->any_array.any.operands[predicate_index.index];
					}

					hol_term* expected_predicate = hol_term::new_apply(&HOL_ANY, hol_term::new_variable(head->quantifier.variable));
					if (expected_predicate == nullptr)
						return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, predicate, expected_predicate);
					free(*expected_predicate); if (expected_predicate->reference_count == 0) free(expected_predicate);
					if (intersection.length > 1) {
						fprintf(stderr, "predicate WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return (hol_term*) nullptr;
					}
					if (intersection.length == 0)
						return (hol_term*) nullptr;

					dst = intersection[0]->binary.left;
					dst->reference_count++;
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool predicate_and_tense(
		hol_term* head, hol_term*& dst,
		hol_term* head_variable,
		hol_term* expected_predicate)
{
	hol_term* expected_head = hol_term::new_exists(head_variable->variable, hol_term::new_and(expected_predicate,
			hol_term::new_apply(hol_term::new_any_constant(make_array_view(PAST_OR_PRESENT, array_length(PAST_OR_PRESENT))), head_variable)));
	if (expected_head == nullptr)
		return false;
	expected_predicate->reference_count++;
	head_variable->reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, head, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (intersection.length > 1) {
		fprintf(stderr, "predicate_and_tense ERROR: Expected intersection size to be 1.\n");
		free_all(intersection); return false;
	} else if (intersection.length == 0) {
		return false;
	}
	dst = intersection[0];
	return true;
}

inline bool predicate_and_tense(
		hol_term* src, hol_term*& dst)
{
	head_index predicate_index; no_op apply;
	auto find_array_head = make_array_finder(find_head<built_in_predicates>);
	hol_term* head = find_head(src, predicate_index, find_array_head, apply);
	if (head == nullptr)
		return false;

	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else if (head->type == hol_term_type::ANY_ARRAY) {
		if (head->any_array.all->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.all->quantifier.variable;
		} else if (head->any_array.right.length != 0 && head->any_array.right.operands[head->any_array.right.length - 1]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.right.operands[head->any_array.right.length - 1]->quantifier.variable;
		} else {
			head_variable = 0;
		}
	} else if (head->type == hol_term_type::AND || head->type == hol_term_type::OR) {
		return false;
	} else {
		head_variable = 0;
	}

	if (head_variable == 0) {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* src_head;
	if (head->type == hol_term_type::ANY_ARRAY) {
		/* make sure this ANY_ARRAY can be a singleton */
		if (head->any_array.left.length > 1 || head->any_array.right.length > 1 || head->any_array.any.length > 1)
			return false;

		array<hol_term*> intersection(2);
		if (head->any_array.left.length != 0) {
			intersection[0] = head->any_array.left.operands[0];
			intersection[0]->reference_count++;
			intersection.length++;
		} if (head->any_array.right.length != 0) {
			if (intersection.length == 0) {
				intersection[0] = head->any_array.right.operands[0];
				intersection[0]->reference_count++;
				intersection.length++;
			} else {
				array<hol_term*> new_intersection(2);
				for (hol_term* term : intersection)
					intersect<built_in_predicates>(new_intersection, term, head->any_array.right.operands[0]);
				free_all(intersection);
				swap(intersection, new_intersection);
				if (intersection.length == 0)
					return false;
			}
		} if (head->any_array.any.length != 0) {
			if (intersection.length == 0) {
				intersection[0] = head->any_array.any.operands[0];
				intersection[0]->reference_count++;
				intersection.length++;
			} else {
				array<hol_term*> new_intersection(2);
				for (hol_term* term : intersection)
					intersect<built_in_predicates>(new_intersection, term, head->any_array.any.operands[0]);
				free_all(intersection);
				swap(intersection, new_intersection);
				if (intersection.length == 0)
					return false;
			}
		}

		if (intersection.length == 0) {
			src_head = head->any_array.all;
			src_head->reference_count++;
		} else if (intersection.length != 1) {
			fprintf(stderr, "predicate_and_tense ERROR: Intersection is not unique.\n");
			free_all(intersection); return false;
		} else {
			src_head = intersection[0];
		}
	} else {
		src_head = head;
		src_head->reference_count++;
	}

	hol_term* var_term = hol_term::new_variable(head_variable);
	if (var_term == nullptr) {
		free(*src_head); if (src_head->reference_count == 0) free(src_head);
		return false;
	}

	hol_term* expected_predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), var_term);
	if (expected_predicate == nullptr) {
		free(*src_head); if (src_head->reference_count == 0) free(src_head);
		free(*var_term); free(var_term); return false;
	}
	var_term->reference_count++;
	hol_non_head_constants<built_in_predicates>::increment_terms();

	bool result = predicate_and_tense(src_head, dst, var_term, expected_predicate);
	free(*src_head); if (src_head->reference_count == 0) free(src_head);
	free(*var_term); if (var_term->reference_count == 0) free(var_term);
	free(*expected_predicate); if (expected_predicate->reference_count == 0) free(expected_predicate);
	return result;
}

inline bool empty_and_tense(
		hol_term* src, hol_term*& dst)
{
	head_index predicate_index; no_op apply;
	hol_term* head = find_head(src, predicate_index, find_head<built_in_predicates>, apply);
	if (head == nullptr)
		return false;

	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* var_term = hol_term::new_variable(head_variable);
	if (var_term == nullptr)
		return false;

	hol_term* expected_predicate = hol_term::new_apply(&HOL_EMPTY, var_term);
	if (expected_predicate == nullptr) {
		free(*var_term); free(var_term);
		return false;
	}
	HOL_EMPTY.reference_count++;
	var_term->reference_count++;

	bool result = predicate_and_tense(head, dst, var_term, expected_predicate);
	free(*var_term); if (var_term->reference_count == 0) free(var_term);
	free(*expected_predicate); if (expected_predicate->reference_count == 0) free(expected_predicate);
	return result;
}

template<hol_term_type QuantifierType>
inline bool require_predicative_quantifier(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* dst = nullptr;
				unsigned int set_variable, element_variable;
				bool must_have_wide_scope_marker_before_quantifier = false;
				bool could_have_wide_scope_marker_before_quantifier = false;
				bool must_have_wide_scope_marker_before_lambda_term = false;
				bool could_have_wide_scope_marker_before_lambda_term = false;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					set_variable = head_variable;
					element_variable = head_variable + 1;
					could_have_wide_scope_marker_before_lambda_term = true;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "require_predicative_quantifier ERROR: Expected existential quantification of set.\n");
						return (hol_term*) nullptr;
					}
#endif

					hol_term* operand = head->quantifier.operand;
					set_variable = head->quantifier.variable;
					if (operand->type == hol_term_type::ANY) {
						element_variable = ++max_variable;
						could_have_wide_scope_marker_before_lambda_term = true;
					} else {
						hol_term* last;
						if (operand->type == hol_term_type::ANY_ARRAY) {
							if (operand->any_array.right.length == 0) {
								fprintf(stderr, "require_predicative_quantifier ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
								return (hol_term*) nullptr;
							}
							last = operand->any_array.right.operands[operand->any_array.right.length - 1];
						} else if (operand->type == hol_term_type::AND) {
							last = operand->array.operands[operand->array.length - 1];
						} else {
							last = operand;
						}

						bool last_is_any = false;
						if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
							last = last->any.included;
							last_is_any = true;
							could_have_wide_scope_marker_before_quantifier = true;
						} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
							last = last->binary.right;
							must_have_wide_scope_marker_before_quantifier = true;
							could_have_wide_scope_marker_before_lambda_term = false;
						}

						while (last->type == hol_term_type::NOT)
							last = operand->unary.operand;

						if (last->type == hol_term_type::ANY_RIGHT) {
							element_variable = ++max_variable;
							could_have_wide_scope_marker_before_quantifier = true;
						} else if (last->type != QuantifierType && !last_is_any) {
							return (hol_term*) nullptr;
						} else if (last_is_any && last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::VARIABLE
								&& last->binary.left->variable == lambda_variable && last->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = last->binary.right->variable;
							if (last_is_any)
								could_have_wide_scope_marker_before_lambda_term = true;
						} else {
							dst = head;
							head->reference_count++;

							hol_term* inner_right = nullptr;
							if (QuantifierType == hol_term_type::EXISTS && last->quantifier.operand->type == hol_term_type::AND)
								inner_right = last->quantifier.operand->array.operands[1];
							else if (QuantifierType == hol_term_type::FOR_ALL && last->quantifier.operand->type == hol_term_type::IF_THEN)
								inner_right = last->quantifier.operand->binary.right;

							if (inner_right != nullptr) {
								if (inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
									must_have_wide_scope_marker_before_lambda_term = true;
								else if (inner_right->type == hol_term_type::ANY || inner_right->type == hol_term_type::ANY_RIGHT)
									could_have_wide_scope_marker_before_lambda_term = true;
							}
						}
					}
				}

				if (dst == nullptr) {
					hol_term* element_var = hol_term::new_variable(element_variable);
					if (element_var == nullptr)
						return (hol_term*) nullptr;

					hol_term* lambda_apply_term = hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var);
					if (lambda_apply_term == nullptr) {
						free(*element_var); free(element_var);
						return (hol_term*) nullptr;
					}

					if (must_have_wide_scope_marker_before_lambda_term) {
						hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), lambda_apply_term);
						if (temp == nullptr) {
							free(*lambda_apply_term); free(lambda_apply_term);
							return (hol_term*) nullptr;
						}
						lambda_apply_term = temp;
					} else if (could_have_wide_scope_marker_before_lambda_term) {
						hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
						if (excluded_quantifier == nullptr) {
							free(*lambda_apply_term); free(lambda_apply_term);
							return (hol_term*) nullptr;
						}
						HOL_ANY.reference_count++;

						hol_term* temp = hol_term::new_any_right(lambda_apply_term, &excluded_quantifier, 1);
						if (temp == nullptr) {
							free(*lambda_apply_term); free(lambda_apply_term);
							free(*excluded_quantifier); free(excluded_quantifier);
							return (hol_term*) nullptr;
						}
						lambda_apply_term = temp;
					}

					hol_term* quantifier;
					if (QuantifierType == hol_term_type::FOR_ALL) {
						/* check if we need to add a WIDE_SCOPE */
						if (!must_have_wide_scope_marker_before_quantifier && can_have_free_variables(*head)) {
							hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								return (hol_term*) nullptr;
							}

							quantifier = hol_term::new_any_right(hol_term::new_for_all(element_variable, hol_term::new_if_then(
									hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term)),
									&excluded_quantifier, 1);
							if (quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								free(*excluded_quantifier); free(excluded_quantifier);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
						} else {
							quantifier = hol_term::new_apply(
									hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE),
									hol_term::new_for_all(element_variable, hol_term::new_if_then(
										hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term)));
							if (quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
						}
					} else if (QuantifierType == hol_term_type::EXISTS) {
						if (could_have_wide_scope_marker_before_quantifier && !must_have_wide_scope_marker_before_quantifier) {
							hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								return (hol_term*) nullptr;
							}

							quantifier = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(
									hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term)),
									&excluded_quantifier, 1);
							if (quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								free(*excluded_quantifier); free(excluded_quantifier);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
						} else {
							quantifier = hol_term::new_exists(element_variable, hol_term::new_and(
									hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term));
							if (quantifier == nullptr) {
								free(*lambda_apply_term); free(lambda_apply_term);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
						}
					} else {
						fprintf(stderr, "require_predicative_quantifier ERROR: Unsupported `QuantifierType`.\n");
						free(*lambda_apply_term); free(lambda_apply_term);
						return (hol_term*) nullptr;
					}

					if (must_have_wide_scope_marker_before_quantifier) {
						hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), quantifier);
						if (temp == nullptr) {
							free(*quantifier); free(quantifier);
							return (hol_term*) nullptr;
						}
						quantifier = temp;
					}

					hol_term* new_head;
					if (QuantifierType == hol_term_type::FOR_ALL) {
						hol_term* excluded_set_definition = hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, hol_term::new_equals(element_var, &HOL_ANY)));
						if (excluded_set_definition == nullptr) {
							free(*quantifier); free(quantifier);
							return (hol_term*) nullptr;
						}
						element_var->reference_count++;
						HOL_ANY.reference_count++;
						hol_term* left = hol_term::new_any(nullptr, &excluded_set_definition, 1);
						if (left == nullptr) {
							free(*quantifier); free(quantifier);
							free(*excluded_set_definition); free(excluded_set_definition);
							return (hol_term*) nullptr;
						}
						new_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
								make_array_view((hol_term**) nullptr, 0), make_array_view(&left, 1), make_array_view(&quantifier, 1)));
						if (new_head == nullptr) {
							free(*quantifier); free(quantifier);
							free(*left); free(left);
							return (hol_term*) nullptr;
						}
					} else {
						new_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
								make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&quantifier, 1)));
						if (new_head == nullptr) {
							free(*quantifier); free(quantifier);
							return (hol_term*) nullptr;
						}
					}
					HOL_ANY.reference_count++;

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, new_head, head);
					free(*new_head); if (new_head->reference_count == 0) free(new_head);
					if (intersection.length == 0) {
						return (hol_term*) nullptr;
					} else if (intersection.length > 1) {
						fprintf(stderr, "require_predicative_quantifier ERROR: Intersection is not unique.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						return (hol_term*) nullptr;
					}
					dst = intersection[0];
					for (unsigned int i = 1; i < intersection.length; i++) {
						hol_term* term = intersection[i];
						free(*term); if (term->reference_count == 0) free(term);
					}
				}

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool remove_predicative_not(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				bool could_have_not = false;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					head->reference_count++;
					return head;
				}

#if !defined(NDEBUG)
				if (head->type != hol_term_type::EXISTS) {
					fprintf(stderr, "remove_predicative_not ERROR: Expected existential quantification of set.\n");
					return (hol_term*) nullptr;
				}
#endif

				hol_term* operand = head->quantifier.operand;
				if (operand->type == hol_term_type::ANY) {
					head->reference_count++;
					return head;
				}

				hol_term* last;
				if (operand->type == hol_term_type::ANY_ARRAY) {
					if (operand->any_array.right.length == 0) {
						fprintf(stderr, "remove_predicative_not ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
						return (hol_term*) nullptr;
					}
					last = operand->any_array.right.operands[operand->any_array.right.length - 1];
				} else if (operand->type == hol_term_type::AND) {
					last = operand->array.operands[operand->array.length - 1];
				} else {
					last = operand;
				}

				if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
					last = last->any.included;
					could_have_not = true;
				} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					last = last->binary.right;
				}

				if (last->type == hol_term_type::NOT)
					return substitute_head<any_node_position::NONE>(head, last, last->unary.operand);

				if (!could_have_not) {
					return (hol_term*) nullptr;
				} else {
					head->reference_count++;
					return head;
				}
			}, no_op()) && dst != nullptr;
}

template<bool UniqueOutput, bool StorePredicates, bool EqualityOnly, bool AllowSubset>
inline bool select_predicate_in_set(array<hol_term*>& out, hol_term* head,
		hol_term* predicate, unsigned int max_variable, unsigned int lambda_variable)
{
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	hol_term* left = nullptr;
	unsigned int set_variable, element_variable = 0;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		set_variable = ++max_variable;
		element_variable = ++max_variable;
	} else {
#if !defined(NDEBUG)
		if (head->type != hol_term_type::EXISTS) {
			fprintf(stderr, "select_predicate_in_set ERROR: Expected existential quantification of set.\n");
			return false;
		}
#endif

		hol_term* operand = head->quantifier.operand;
		set_variable = head->quantifier.variable;
		if (operand->type == hol_term_type::ANY) {
			element_variable = ++max_variable;
		} else {
			hol_term* last;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (operand->any_array.right.length == 0) {
					fprintf(stderr, "select_predicate_in_set ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
					return false;
				}
				last = operand->any_array.right.operands[operand->any_array.right.length - 1];
				if (operand->any_array.left.length != 0)
					left = operand->any_array.left.operands[0];
			} else if (operand->type == hol_term_type::AND) {
				last = operand->array.operands[operand->array.length - 1];
				left = operand->array.operands[0];
			} else {
				last = operand;
				left = operand;
			}

			if (left != nullptr) {
				hol_term* set_definition = nullptr;
				if ((left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) && left->any.included != nullptr) {
					set_definition = left->any.included;
				} else if (left->type == hol_term_type::EQUALS) {
					set_definition = left->binary.right;
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = left->ternary.third;
				}

				if (set_definition != nullptr && set_definition->type == hol_term_type::LAMBDA)
					element_variable = set_definition->quantifier.variable;
			}

			if (element_variable == 0) {
				bool last_is_any = false;
				if (last->type == hol_term_type::ANY && last->any.included != nullptr) {
					last = last->any.included;
					last_is_any = true;
				} else if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
					last = last->any.included;
					last_is_any = true;
				} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					last = last->binary.right;
				}

				while (last->type == hol_term_type::NOT)
					last = operand->unary.operand;

				if (last->type == hol_term_type::ANY) {
					element_variable = ++max_variable;
				} else if (last_is_any && last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::VARIABLE
						&& last->binary.left->variable == lambda_variable && last->binary.right->type == hol_term_type::VARIABLE)
				{
					element_variable = last->binary.right->variable;
				} else if (last->type == hol_term_type::FOR_ALL || last->type == hol_term_type::EXISTS) {
					element_variable = last->quantifier.variable;
				} else {
					return false;
				}
			}
		}
	}

	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
	if (excluded_quantifier == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* excluded_lambda = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
	if (excluded_lambda == nullptr) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}
	HOL_ANY.reference_count++;

	hol_term* excluded_quantifiers[2];
	excluded_quantifiers[0] = excluded_quantifier;
	excluded_quantifiers[1] = excluded_lambda;
	hol_term* expected_predicate = hol_term::new_any(nullptr, excluded_quantifiers, array_length(excluded_quantifiers));
	if (expected_predicate == nullptr) {
		free(*excluded_lambda); free(excluded_lambda);
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}
	excluded_quantifier->reference_count++;

	array<hol_term*> predicates(2);
	intersect<built_in_predicates>(predicates, expected_predicate, predicate);
	free(*expected_predicate); if (expected_predicate->reference_count == 0) free(expected_predicate);
	if (predicates.length == 0) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	} else if (predicates.length != 1) {
		fprintf(stderr, "select_predicate_in_set ERROR: Not implemented.\n");
		free_all(predicates);
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}

	if (left == nullptr || left->type != hol_term_type::BINARY_APPLICATION) {
		hol_term* first_expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, hol_term::new_equals(hol_term::new_variable(element_variable), predicates[0]))),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
		if (first_expected_head == nullptr) {
			free_all(predicates);
			free(*excluded_quantifier); free(excluded_quantifier);
			return false;
		}
		predicates[0]->reference_count++;
		excluded_quantifier->reference_count++;

		intersect<built_in_predicates>(out, first_expected_head, head);
		free(*first_expected_head); if (first_expected_head->reference_count == 0) free(first_expected_head);
		if (UniqueOutput && out.length > 1) {
			fprintf(stderr, "select_predicate_in_set ERROR: Intersection is not unique.\n");
			for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}

		for (unsigned int i = 0; StorePredicates && i < out.length; i++) {
			hol_term* term = out[i]->quantifier.operand->array.operands[0]->binary.right->quantifier.operand->binary.right;
			term->reference_count++;
			free(*out[i]); if (out[i]->reference_count == 0) free(out[i]);
			out[i] = term;
		}

		if ((UniqueOutput && out.length == 1) || EqualityOnly) {
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return true;
		}
	}

	if (EqualityOnly) {
		free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
		free_all(predicates); return false;
	}

	array<hol_term*> second_expected_heads(2);
	if (UniqueOutput && (left == nullptr || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT)) {
		hol_term* excluded_set_definition = hol_term::new_equals(hol_term::new_variable(set_variable),
				hol_term::new_lambda(element_variable, hol_term::new_equals(hol_term::new_variable(element_variable), predicates[0])));
		if (excluded_set_definition == nullptr) {
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}
		predicates[0]->reference_count++;

		second_expected_heads[second_expected_heads.length] = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_any_right(hol_term::new_lambda(element_variable, hol_term::new_apply(predicates[0], hol_term::new_variable(element_variable))), &excluded_set_definition, 1),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
		if (second_expected_heads[second_expected_heads.length] == nullptr) {
			free(*excluded_set_definition); free(excluded_set_definition);
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}
		predicates[0]->reference_count++;
		excluded_quantifier->reference_count++;
		second_expected_heads.length++;
	} if (!UniqueOutput || left->type == hol_term_type::EQUALS) {
		second_expected_heads[second_expected_heads.length] = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, hol_term::new_apply(predicates[0], hol_term::new_variable(element_variable)))),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
		if (second_expected_heads[second_expected_heads.length] == nullptr) {
			free_all(second_expected_heads);
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}
		predicates[0]->reference_count++;
		excluded_quantifier->reference_count++;
		second_expected_heads.length++;
	} if (AllowSubset && (!UniqueOutput || left->type == hol_term_type::BINARY_APPLICATION)) {
		hol_term* set_var = hol_term::new_variable(set_variable);
		if (set_var == nullptr) {
			free_all(second_expected_heads);
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}

		hol_term* element_var = hol_term::new_variable(element_variable);
		if (element_var == nullptr) {
			free_all(second_expected_heads);
			free(*set_var); free(set_var);
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}

		second_expected_heads[second_expected_heads.length] = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SUBSET), set_var, hol_term::new_lambda(element_variable, hol_term::new_apply(predicates[0], element_var))),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
		if (second_expected_heads[second_expected_heads.length] == nullptr) {
			free_all(second_expected_heads);
			free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
			free_all(predicates); return false;
		}
		predicates[0]->reference_count++;
		excluded_quantifier->reference_count++;
		second_expected_heads.length++;
	}
	free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
	free_all(predicates);

	unsigned int old_out_length = out.length;
	for (hol_term* second_expected_head : second_expected_heads) {
		intersect<built_in_predicates>(out, second_expected_head, head);
		if (UniqueOutput && out.length > 1) {
			fprintf(stderr, "select_predicate_in_set ERROR: Intersection is not unique.\n");
			free_all(second_expected_heads);
			for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	}
	free_all(second_expected_heads);

	for (unsigned int i = old_out_length; StorePredicates && i < out.length; i++) {
		hol_term* term;
		hol_term* left = out[i]->quantifier.operand->array.operands[0];
		if (left->type == hol_term_type::ANY_RIGHT) {
			term = left->any.included->quantifier.operand->binary.left;
		} else if (left->type == hol_term_type::EQUALS) {
			term = left->binary.right->quantifier.operand->binary.left;
		} else {
			term = left->ternary.third->quantifier.operand->binary.left;
		}
		term->reference_count++;
		free(*out[i]); if (out[i]->reference_count == 0) free(out[i]);
		out[i] = term;
	}
	return true;
}

inline bool select_predicate_in_set(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				array<hol_term*> out(2);
				if (!select_predicate_in_set<true, true, false, true>(out, head, &HOL_ANY, max_variable, lambda_variable) || out.length == 0)
					return (hol_term*) nullptr;
				remove_wide_scope_marker = true;
				return out[0];
			}, no_op()) && dst != nullptr;
}

inline bool mark_wide_scope(hol_term* src, hol_term*& dst)
{
	head_index predicate_index; no_op apply;
	hol_term* head = find_head(src, predicate_index, find_head_or_universal<built_in_predicates>, apply);
	if (head == nullptr)
		return false;

	hol_term* new_head;
	if (head->type == hol_term_type::FOR_ALL) {
		new_head = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), head);
		if (new_head == nullptr) return false;
		head->reference_count++;
	} else {
		new_head = head;
		head->reference_count++;
	}

	hol_term* new_term = substitute_head<any_node_position::NONE>(src, head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	if (new_term == nullptr)
		return false;

	hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
	if (any_wide_scope == nullptr) {
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		return false;
	}
	HOL_ANY.reference_count++;

	hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
	if (duplicate_wide_scopes == nullptr) {
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		free(*any_wide_scope); free(any_wide_scope);
		return false;
	}

	hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
	if (excluded_universal == nullptr) {
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		return false;
	}
	any_wide_scope->reference_count++;

	hol_term* excluded_trees[2];
	excluded_trees[0] = excluded_universal;
	excluded_trees[1] = duplicate_wide_scopes;
	hol_term* unique_wide_scope = hol_term::new_any_right(nullptr, excluded_trees, array_length(excluded_trees));
	if (unique_wide_scope == nullptr) {
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		free(*excluded_universal); free(excluded_universal);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		return false;
	}
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, new_term, unique_wide_scope);
	free(*new_term); if (new_term->reference_count == 0) free(new_term);
	free(*unique_wide_scope); if (unique_wide_scope->reference_count == 0) free(unique_wide_scope);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length > 1) {
		fprintf(stderr, "mark_wide_scope ERROR: Intersection is not unique.\n");
		free_all(intersection);
		return false;
	}
	dst = intersection[0];
	return true;
}

template<bool WideScope>
inline hol_term* do_require_narrow_or_wide_scope(hol_term* head,
		unsigned int lambda_variable, unsigned int max_variable, bool can_have_wide_scope)
{
	hol_term* old_head = head;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	hol_term* dst = nullptr;
	unsigned int set_variable, element_variable = 0;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		set_variable = ++max_variable;
		element_variable = ++max_variable;
	} else {
#if !defined(NDEBUG)
		if (head->type != hol_term_type::EXISTS) {
			fprintf(stderr, "do_require_narrow_or_wide_scope ERROR: Expected existential quantification of set.\n");
			return (hol_term*) nullptr;
		}
#endif

		hol_term* operand = head->quantifier.operand;
		set_variable = head->quantifier.variable;
		if (operand->type == hol_term_type::ANY) {
			element_variable = ++max_variable;
		} else {
			hol_term* last;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (operand->any_array.right.length == 0) {
					fprintf(stderr, "do_require_narrow_or_wide_scope ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
					return (hol_term*) nullptr;
				}
				last = operand->any_array.right.operands[operand->any_array.right.length - 1];
			} else if (operand->type == hol_term_type::AND) {
				last = operand->array.operands[operand->array.length - 1];
			} else {
				last = operand;
			}

			if (last->type == hol_term_type::ANY && last->any.included != nullptr) {
				last = last->any.included;
			} else if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
				last = last->any.included;
			} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
				last = last->binary.right;
				if (WideScope) {
					head->reference_count++;
					return head;
				} else {
					return (hol_term*) nullptr;
				}
			} else if (WideScope && !can_have_wide_scope) {
				/* it's not possible to have wide scope */
				return (hol_term*) nullptr;
			}

			while (last->type == hol_term_type::NOT)
				last = operand->unary.operand;

			if (last->type == hol_term_type::ANY || last->type == hol_term_type::ANY_RIGHT) {
				element_variable = ++max_variable;
			} else if (last->type == hol_term_type::EXISTS || last->type == hol_term_type::FOR_ALL) {

				hol_term* inner_right = nullptr;
				if (last->type == hol_term_type::EXISTS && last->quantifier.operand->type == hol_term_type::AND)
					inner_right = last->quantifier.operand->array.operands[1];
				else if (last->type == hol_term_type::FOR_ALL && last->quantifier.operand->type == hol_term_type::IF_THEN)
					inner_right = last->quantifier.operand->binary.right;

				if (inner_right != nullptr) {
					if (inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
						if (WideScope) {
							return (hol_term*) nullptr;
						} else {
							head->reference_count++;
							return head;
						}
					}
				}

				dst = head;
				head->reference_count++;
			}
		}
	}

	if (dst == nullptr) {
		hol_term* element_var = hol_term::new_variable(element_variable);
		if (element_var == nullptr)
			return (hol_term*) nullptr;

		hol_term* lambda_apply_term = hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var);
		if (lambda_apply_term == nullptr) {
			free(*element_var); free(element_var);
			return (hol_term*) nullptr;
		}

		if (!WideScope) {
			hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), lambda_apply_term);
			if (temp == nullptr) {
				free(*lambda_apply_term); free(lambda_apply_term);
				return (hol_term*) nullptr;
			}
			lambda_apply_term = temp;
		}

		hol_term* quantifier = hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL_OR_EXISTS,
				hol_term::new_if_then(hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term));
		if (quantifier == nullptr) {
			free(*lambda_apply_term); free(lambda_apply_term);
			return (hol_term*) nullptr;
		}
		element_var->reference_count++;

		if (WideScope) {
			hol_term* temp = hol_term::new_any_right(quantifier);
			if (temp == nullptr) {
				free(*quantifier); free(quantifier);
				return (hol_term*) nullptr;
			}
			quantifier = temp;
		}

		hol_term* new_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
				make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&quantifier, 1)));
		if (new_head == nullptr) {
			free(*quantifier); free(quantifier);
			return (hol_term*) nullptr;
		}
		HOL_ANY.reference_count++;

		array<hol_term*> intersection(2);
		intersect<built_in_predicates>(intersection, new_head, head);
		free(*new_head); if (new_head->reference_count == 0) free(new_head);
		if (intersection.length == 0) {
			return (hol_term*) nullptr;
		} else if (intersection.length > 1) {
			fprintf(stderr, "do_require_narrow_or_wide_scope ERROR: Intersection is not unique.\n");
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return (hol_term*) nullptr;
		}
		dst = intersection[0];
		for (unsigned int i = 1; i < intersection.length; i++) {
			hol_term* term = intersection[i];
			free(*term); if (term->reference_count == 0) free(term);
		}
	}

	if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
		hol_term* new_dst = hol_term::new_any_right(dst);
		if (new_dst == nullptr) {
			free(*dst); if (dst->reference_count == 0) free(dst);
			return (hol_term*) nullptr;
		}
		dst = new_dst;
	}
	return dst;
}

template<bool WideScope>
bool require_narrow_or_wide_scope(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	bool can_have_wide_scope = false;
	auto check_wide_scope = [&can_have_wide_scope](hol_term* term) {
		if (term->type == hol_term_type::ANY || term->type == hol_term_type::ANY_RIGHT) {
			can_have_wide_scope = true;
		} else if (term->type == hol_term_type::UNARY_APPLICATION && term->binary.left->type == hol_term_type::CONSTANT && term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
			can_have_wide_scope = true;
		}
		return true;
	};

	head_index predicate_index;
	auto head_finder = predicative_head_finder<built_in_predicates>(lambda_variable);
	hol_term* head = find_head(src, predicate_index, head_finder, check_wide_scope);
	if (head == nullptr)
		return false;

	hol_term* new_head = do_require_narrow_or_wide_scope<WideScope>(head, lambda_variable, max_variable, can_have_wide_scope);
	if (new_head == nullptr)
		return false;

	dst = substitute_head<any_node_position::NONE>(src, head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	return (dst != nullptr);
}

bool remove_wide_scope(hol_term* src, hol_term*& dst)
{
	head_index predicate_index; no_op apply;
	hol_term* head = find_head(src, predicate_index, find_head_or_unary_application<built_in_predicates, (unsigned int) built_in_predicates::WIDE_SCOPE>, apply);
	if (head == nullptr) return false;

	hol_term* new_head;
	if (head->type == hol_term_type::UNARY_APPLICATION && head->binary.left->type == hol_term_type::CONSTANT && head->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
		new_head = head->binary.right;
	} else {
		new_head = head;
	}

	hol_term* new_term = substitute_head<any_node_position::NONE>(src, head, new_head);
	if (new_term == nullptr)
		return false;

	hol_term* excluded = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
	if (excluded == nullptr) {
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		return false;
	}
	HOL_ANY.reference_count++;

	array<hol_term*> difference(2);
	subtract<built_in_predicates>(difference, new_term, excluded);
	free(*new_term); if (new_term->reference_count == 0) free(new_term);
	free(*excluded); if (excluded->reference_count == 0) free(excluded);
	if (difference.length == 0) {
		return false;
	} else if (difference.length != 1) {
		fprintf(stderr, "remove_wide_scope ERROR: Set difference is not unique.\n");
		return false;
	}
	dst = difference[0];
	return true;
}

inline bool require_constant_in_set(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* any_constant = hol_term::new_any_constant_except();
				if (any_constant == nullptr)
					return (hol_term*) nullptr;
				array<hol_term*> out(2);
				if (!select_predicate_in_set<true, false, true, true>(out, head, any_constant, max_variable, lambda_variable) || out.length == 0) {
					free(*any_constant); free(any_constant);
					return (hol_term*) nullptr;
				}
				free(*any_constant); if (any_constant->reference_count == 0) free(any_constant);
				return out[0];
			}, no_op()) && dst != nullptr;
}

inline bool require_no_constant_in_set(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				unsigned int set_variable, element_variable;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					set_variable = ++max_variable;
					element_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "require_no_constant_in_set ERROR: Expected existential quantification of set.\n");
						return (hol_term*) nullptr;
					}
#endif

					hol_term* operand = head->quantifier.operand;
					set_variable = head->quantifier.variable;
					if (operand->type == hol_term_type::ANY) {
						element_variable = ++max_variable;
					} else {
						hol_term* last;
						if (operand->type == hol_term_type::ANY_ARRAY) {
							if (operand->any_array.right.length == 0) {
								fprintf(stderr, "require_no_constant_in_set ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
								return (hol_term*) nullptr;
							}
							last = operand->any_array.right.operands[operand->any_array.right.length - 1];
						} else if (operand->type == hol_term_type::AND) {
							last = operand->array.operands[operand->array.length - 1];
						} else {
							last = operand;
						}

						bool last_is_any = false;
						if (last->type == hol_term_type::ANY && last->any.included != nullptr) {
							last = last->any.included;
							last_is_any = true;
						} else if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
							last = last->any.included;
							last_is_any = true;
						} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
							last = last->binary.right;
						}

						while (last->type == hol_term_type::NOT)
							last = operand->unary.operand;

						if (last->type == hol_term_type::ANY) {
							element_variable = ++max_variable;
						} else if (last_is_any && last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::VARIABLE
								&& last->binary.left->variable == lambda_variable && last->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = last->binary.right->variable;
						} else if (last->type == hol_term_type::FOR_ALL || last->type == hol_term_type::EXISTS) {
							element_variable = last->quantifier.variable;
						} else {
							return (hol_term*) nullptr;
						}
					}
				}

				hol_term* excluded_trees[3];
				excluded_trees[0] = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
				excluded_trees[1] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
				excluded_trees[2] = hol_term::new_equals(hol_term::new_variable(element_variable), &HOL_ANY);
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					return (hol_term*) nullptr;
				}

				hol_term* set_definition = hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable,
						hol_term::new_any_array(hol_term_type::AND, hol_term::new_any(nullptr, excluded_trees, array_length(excluded_trees)),
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0))));
				if (set_definition == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					return (hol_term*) nullptr;
				}

				hol_term* expected_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND,
						hol_term::new_any(nullptr, &excluded_trees[0], 1), make_array_view((hol_term**) nullptr, 0),
						make_array_view(&set_definition, 1), make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*set_definition); free(set_definition);
					return (hol_term*) nullptr;
				}
				excluded_trees[0]->reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length > 1) {
					fprintf(stderr, "require_no_constant_in_set ERROR: Intersection is not unique.\n");
					free_all(intersection); return (hol_term*) nullptr;
				}

				hol_term* dst = intersection[0];
				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_term = hol_term::new_any_right(dst, old_head->any.excluded_trees, old_head->any.excluded_tree_count);
					if (new_term == nullptr) {
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					for (unsigned int j = 0; j < old_head->any.excluded_tree_count; j++)
						old_head->any.excluded_trees[j]->reference_count++;
					dst = new_term;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool require_singleton(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				unsigned int set_variable;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					set_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "require_singleton ERROR: Expected existential quantification of set.\n");
						return (hol_term*) nullptr;
					}
#endif
					set_variable = head->quantifier.variable;
				}

				hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
				if (excluded_quantifier == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				hol_term* conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
				if (conjunct == nullptr) {
					free(*excluded_quantifier); free(excluded_quantifier);
					return (hol_term*) nullptr;
				}

				hol_term* set_size_conjunct = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SIZE), hol_term::new_variable(set_variable)), hol_term::new_int(1));
				if (set_size_conjunct == nullptr) {
					free(*conjunct); free(conjunct);
					return (hol_term*) nullptr;
				}

				hol_term* expected_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
						make_array_view((hol_term**) nullptr, 0),
						make_appended_array_view(make_repeated_array_view(conjunct, 1), set_size_conjunct),
						make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*set_size_conjunct); free(set_size_conjunct);
					free(*conjunct); free(conjunct);
					return (hol_term*) nullptr;
				}
				conjunct->reference_count += 2 - 1;

				array<hol_term*> out(2);
				intersect<built_in_predicates>(out, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (out.length == 0) {
					return (hol_term*) nullptr;
				} else if (out.length > 1) {
					fprintf(stderr, "require_singleton ERROR: Intersection is not unique.\n");
					for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
					return (hol_term*) nullptr;
				}

				hol_term* dst = out[0];
				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool size(hol_term* src, hol_term*& dst, grammatical_flags& dst_flags)
{
	unsigned int head_variable;
	if (src->type == hol_term_type::EXISTS) {
		head_variable = src->quantifier.variable;
	} else if ((src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) && src->any.included != nullptr && src->any.included->type == hol_term_type::EXISTS) {
		src = src->any.included;
		head_variable = src->quantifier.variable;
	} else {
		head_variable = 1;
	}

	hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_equals(hol_term::new_apply(
			hol_term::new_constant((unsigned int) built_in_predicates::SIZE), hol_term::new_variable(head_variable)), &HOL_ANY));
	if (expected_head == nullptr)
		return (hol_term*) nullptr;
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, src, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (intersection.length == 0) {
		return (hol_term*) nullptr;
	} else if (intersection.length > 1) {
		fprintf(stderr, "size ERROR: Intersection is not unique.\n");
		free_all(intersection);
		return (hol_term*) nullptr;
	}

	dst = intersection[0]->quantifier.operand->binary.right;
	if (dst->type == hol_term_type::INTEGER) {
		if (dst->integer == 1) {
			if (!intersect(dst_flags.index_number, grammatical_num::SINGULAR, dst_flags.index_number)) {
				free_all(intersection);
				return false;
			}
		} else {
			if (!intersect(dst_flags.index_number, grammatical_num::PLURAL, dst_flags.index_number)) {
				free_all(intersection);
				return false;
			}
		}
	}

	dst->reference_count++;
	free_all(intersection);
	return true;
}

inline bool set_size(hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				unsigned int set_variable;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					set_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "set_size ERROR: Expected existential quantification of set.\n");
						return (hol_term*) nullptr;
					}
#endif
					set_variable = head->quantifier.variable;
				}

				hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
				if (excluded_quantifier == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
				if (expected_conjunct == nullptr) {
					free(*excluded_quantifier); free(excluded_quantifier);
					return (hol_term*) nullptr;
				}

				hol_term* expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
						expected_conjunct,
						hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SIZE), hol_term::new_variable(set_variable)), expected_conjunct),
						expected_conjunct));
				if (expected_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return (hol_term*) nullptr;
				}
				expected_conjunct->reference_count += 3 - 1;

				array<hol_term*> out(2);
				intersect<built_in_predicates>(out, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (out.length == 0) {
					return (hol_term*) nullptr;
				} else if (out.length != 1) {
					fprintf(stderr, "set_size ERROR: Intersection is not unique.\n");
					for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_quantifier); if (excluded_quantifier->reference_count == 0) free(excluded_quantifier);
					return (hol_term*) nullptr;
				}

				hol_term* dst = out[0]->quantifier.operand->array.operands[1]->binary.right;
				dst->reference_count++;
				free_all(out);
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool remove_predicate(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* old_head = head;
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				bool negated = false;
				if (head->type == hol_term_type::NOT) {
					head = head->unary.operand;
					negated = true;
				}

				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
				 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
				{
					hol_term* head_var = hol_term::new_variable(head_variable);
					if (head_var == nullptr) return (hol_term*) nullptr;
					constexpr unsigned int excluded_tree_count = 4;
					hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term) * (excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()));
					excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
					excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
					excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
					excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
					if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
						if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
						if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
						if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
						free(*head_var); free(head_var);
						return (hol_term*) nullptr;
					}
					free(*head_var);

					for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++) {
						excluded_trees[excluded_tree_count + i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
						excluded_trees[excluded_tree_count + i]->reference_count++;
					}

					dst = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND,
							hol_term::new_any(nullptr, excluded_trees, excluded_tree_count), make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
					if (dst == nullptr) {
						for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
							free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
						}
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < excluded_tree_count; i++)
						excluded_trees[i]->reference_count++;
					for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
						free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
					}

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "remove_predicate WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count == 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;

					if (ConjunctIndex >= 0 && (predicate_index.position != head_position::LEFT || predicate_index.index != ConjunctIndex))
						return (hol_term*) nullptr;
					if (ConjunctIndex < 0 && (predicate_index.position != head_position::RIGHT || predicate_index.index != (unsigned int) (-ConjunctIndex) - 1))
						return (hol_term*) nullptr;

					hol_term* conjunction = nullptr;
					if (ConjunctIndex >= 0) {
						if (ConjunctIndex < operand->any_array.left.length) {
							conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
									make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
									make_excluded_array_view(operand->any_array.left.operands, operand->any_array.left.length, ConjunctIndex),
									make_array_view(operand->any_array.right.operands, operand->any_array.right.length));
						} else {
							conjunction = operand;
						}
					} else if (ConjunctIndex < 0) {
						if (-ConjunctIndex - 1 < operand->any_array.right.length) {
							unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
							index = operand->any_array.right.length - index - 1;
							conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
									make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
									make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
									make_excluded_array_view(operand->any_array.right.operands, operand->any_array.right.length, index));
						} else {
							conjunction = operand;
						}
					}
					if (conjunction == nullptr)
						return (hol_term*) nullptr;
					if (conjunction == operand) {
						operand->reference_count++;
					} else {
						conjunction->any_array.all->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
							conjunction->any_array.any.operands[i]->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
							conjunction->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
							conjunction->any_array.right.operands[i]->reference_count++;
					}

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
						return (hol_term*) nullptr;
					}
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
						fprintf(stderr, "remove_predicate WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					unsigned int conjunct_index = (ConjunctIndex < 0) ? (operand->array.length + ConjunctIndex) : ConjunctIndex;
					if (predicate_index.index != conjunct_index)
						return (hol_term*) nullptr;
					if (operand->array.length == 2) {
						dst = hol_term::new_exists(head_variable, (conjunct_index == 0 ? operand->array.operands[1] : operand->array.operands[0]));
					} else {
						dst = hol_term::new_exists(head_variable, hol_term::new_and(make_excluded_array_view(operand->array.operands, operand->array.length, conjunct_index)));
					}
					if (dst != nullptr) {
						for (unsigned int i = 0; i < operand->array.length; i++)
							if (i != conjunct_index) operand->array.operands[i]->reference_count++;
					}
				}

				if (negated) {
					hol_term* new_dst = hol_term::new_not(dst);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					dst = new_dst;
				}

				if (old_head->type == hol_term_type::ANY || old_head->type == hol_term_type::ANY_RIGHT) {
					hol_term* new_dst = hol_term::new_any_right(dst, old_head->any.excluded_trees, old_head->any.excluded_tree_count);
					if (new_dst == nullptr) {
						free(*dst); if (dst->reference_count == 0) free(dst);
						return (hol_term*) nullptr;
					}
					for (unsigned int i = 0; i < old_head->any.excluded_tree_count; i++)
						old_head->any.excluded_trees[i]->reference_count++;
					dst = new_dst;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<typename FindHeadFunction>
inline bool factor(
		hol_term* src, hol_term*& dst,
		FindHeadFunction find_head_function)
{
	hol_term* parent = nullptr;
	hol_term* current = nullptr;
	auto apply = [&parent,&current](hol_term* term) {
		parent = current;
		current = term;
	};

	head_index predicate_index;
	auto find_array_head = make_array_finder(find_head_function);
	hol_term* head = find_head(src, predicate_index, find_array_head, apply);
	if (head == nullptr)
		return false;

	hol_term* new_head;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		hol_term* excluded_array = hol_term::new_any_array(
				hol_term_type::ANY_ARRAY, &HOL_ANY, make_array_view((hol_term**) nullptr, 0),
				make_repeated_array_view(&HOL_ANY, 2), make_array_view((hol_term**) nullptr, 0));
		if (excluded_array == nullptr)
			return false;
		HOL_ANY.reference_count += 3;

		array<hol_term*> difference(2);
		subtract<built_in_predicates>(difference, head, excluded_array);
		free(*excluded_array); if (excluded_array->reference_count == 0) free(excluded_array);
		if (difference.length == 0) {
			return false;
		} else if (difference.length != 1) {
			fprintf(stderr, "factor ERROR: Set difference is not unique.\n");
			free_all(difference);
			return false;
		}
		new_head = difference[0];

	} else if (head->type == hol_term_type::AND || head->type == hol_term_type::OR) {
		array<hol_term*> intersection(2);
		intersection[intersection.length] = head->array.operands[0];
		intersection[intersection.length++]->reference_count++;
		for (unsigned int i = 1; i < head->array.length; i++) {
			array<hol_term*> new_intersection(4);
			for (hol_term* prev : intersection)
				intersect<built_in_predicates>(new_intersection, prev, head->array.operands[i]);
			free_all(intersection);
			swap(new_intersection, intersection);
		}

		if (intersection.length == 0) {
			return false;
		} else if (intersection.length != 1) {
			fprintf(stderr, "factor ERROR: Set intersection is not unique.\n");
			free_all(intersection);
			return false;
		}
		new_head = intersection[0];

	} else if (head->type == hol_term_type::ANY_ARRAY && (head->any_array.oper == hol_term_type::ANY_ARRAY || head->any_array.oper == hol_term_type::AND || head->any_array.oper == hol_term_type::OR)) {
		array<hol_term*> intersection(2);
		intersection[intersection.length] = head->any_array.all;
		intersection[intersection.length++]->reference_count++;
		for (unsigned int i = 0; i < head->any_array.left.length; i++) {
			array<hol_term*> new_intersection(4);
			for (hol_term* prev : intersection)
				intersect<built_in_predicates>(new_intersection, prev, head->any_array.left.operands[i]);
			free_all(intersection);
			swap(new_intersection, intersection);
		} for (unsigned int i = 0; i < head->any_array.right.length; i++) {
			array<hol_term*> new_intersection(4);
			for (hol_term* prev : intersection)
				intersect<built_in_predicates>(new_intersection, prev, head->any_array.right.operands[i]);
			free_all(intersection);
			swap(new_intersection, intersection);
		} for (unsigned int i = 0; i < head->any_array.any.length; i++) {
			array<hol_term*> new_intersection(4);
			for (hol_term* prev : intersection)
				intersect<built_in_predicates>(new_intersection, prev, head->any_array.any.operands[i]);
			free_all(intersection);
			swap(new_intersection, intersection);
		}

		if (intersection.length == 0) {
			return false;
		} else if (intersection.length != 1) {
			fprintf(stderr, "factor ERROR: Set intersection is not unique.\n");
			free_all(intersection);
			return false;
		}
		new_head = intersection[0];
	} else {
		new_head = head;
		new_head->reference_count++;
	}

	if ((new_head->type == hol_term_type::ANY || new_head->type == hol_term_type::ANY_RIGHT)
	 && parent != nullptr && (parent->type == hol_term_type::ANY || parent->type == hol_term_type::ANY_RIGHT))
	{
		hol_term_type new_type = (new_head->type == hol_term_type::ANY || parent->type == hol_term_type::ANY) ? hol_term_type::ANY : hol_term_type::ANY_RIGHT;
		if (parent->any.excluded_tree_count == 0 && new_type == new_head->type) {
			/* the merger of `parent` and `new_head` is just `new_head` (ignoring `included`) */
			dst = substitute_head<any_node_position::NONE>(src, parent, new_head);
			free(*new_head); if (new_head->reference_count == 0) free(new_head);
		} else {
			array<hol_term*> excluded(max(1u, parent->any.excluded_tree_count + new_head->any.excluded_tree_count));
			for (unsigned int i = 0; i < parent->any.excluded_tree_count; i++)
				excluded[excluded.length++] = parent->any.excluded_trees[i];
			for (unsigned int i = 0; i < new_head->any.excluded_tree_count; i++)
				excluded[excluded.length++] = new_head->any.excluded_trees[i];
			hol_term* temp;
			if (new_type == hol_term_type::ANY)
				temp = hol_term::new_any(new_head->any.included, excluded.data, excluded.length);
			else temp = hol_term::new_any_right(new_head->any.included, excluded.data, excluded.length);
			if (new_head->any.included != nullptr)
				new_head->any.included->reference_count++;
			for (hol_term* term : excluded)
				term->reference_count++;
			free(*new_head); if (new_head->reference_count == 0) free(new_head);
			dst = substitute_head<any_node_position::NONE>(src, parent, temp);
			free(*temp); if (temp->reference_count == 0) free(temp);
		}
	} else {
		dst = substitute_head<any_node_position::NONE>(src, head, new_head);
		free(*new_head); if (new_head->reference_count == 0) free(new_head);
	}
	return (dst != nullptr);
}

inline bool factor(
		hol_term* src, hol_term*& dst)
{
	return factor(src, dst, find_head<built_in_predicates>);
}

inline bool factor_predicative(
		hol_term* src, hol_term*& dst)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	return factor(src, dst, predicative_head_finder<built_in_predicates>(lambda_variable));
}

template<int_fast8_t ConjunctIndex>
inline bool select_predicate_and_tense(
		hol_term* src, hol_term*& dst)
{
	static_assert(ConjunctIndex != -1, "select_predicate_and_tense ERROR: `ConjunctIndex` cannot be -1.");

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, find_head<built_in_predicates>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && predicate_index.position != head_position::NONE)
					head = head->any.included;

				hol_term* head_var = hol_term::new_variable(head_variable);
				if (head_var == nullptr) return (hol_term*) nullptr;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term) * (excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()));
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return (hol_term*) nullptr;
				}
				free(*head_var);

				for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++) {
					excluded_trees[excluded_tree_count + i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
					excluded_trees[excluded_tree_count + i]->reference_count++;
				}

				hol_term* expected_predicate = hol_term::new_apply(
							hol_term::new_any(nullptr, excluded_trees, excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()),
							hol_term::new_variable(head_variable));
				if (expected_predicate == nullptr) {
					for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
						free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
					}
					return (hol_term*) nullptr;
				}

				hol_term* tense_conjunct = hol_term::new_apply(hol_term::new_any_constant(make_array_view(TENSE_PREDICATES, array_length(TENSE_PREDICATES))), hol_term::new_variable(head_variable));
				if (tense_conjunct == nullptr) {
					free(*expected_predicate); free(expected_predicate);
					return (hol_term*) nullptr;
				}

				hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
				if (conjunct == nullptr) {
					free(*expected_predicate); free(expected_predicate);
					free(*tense_conjunct); free(tense_conjunct);
					return (hol_term*) nullptr;
				}
				for (unsigned int i = 0; i < excluded_tree_count; i++)
					excluded_trees[i]->reference_count++;

				hol_term* expected_head;
				if (ConjunctIndex >= 0) {
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, conjunct, make_array_view((hol_term**) nullptr, 0),
							make_appended_array_view(make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), expected_predicate), tense_conjunct),
							make_array_view((hol_term**) nullptr, 0)));
					if (expected_head == nullptr) {
						free(*expected_predicate); free(expected_predicate);
						free(*tense_conjunct); free(tense_conjunct);
						free(*conjunct); free(conjunct);
						return (hol_term*) nullptr;
					}
					conjunct->reference_count += ConjunctIndex;
				} else {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 2;
					expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0),
							make_prepended_array_view(expected_predicate, make_prepended_array_view(tense_conjunct, make_repeated_array_view(conjunct, index)))));
					if (expected_head == nullptr) {
						free(*expected_predicate); free(expected_predicate);
						free(*tense_conjunct); free(tense_conjunct);
						free(*conjunct); free(conjunct);
						return (hol_term*) nullptr;
					}
					conjunct->reference_count += index;
				}

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, expected_head, head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length > 1) {
					fprintf(stderr, "select_predicate_and_tense ERROR: Intersection is not unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}

				hol_term* left = nullptr; hol_term* right = nullptr;
				hol_term* operand = intersection[0]->quantifier.operand;
				if (operand->type == hol_term_type::AND) {
					unsigned index = (ConjunctIndex >= 0) ? ConjunctIndex : (operand->array.length + ConjunctIndex);
					left = operand->array.operands[index];
					right = operand->array.operands[index + 1];
				} else if (operand->type == hol_term_type::ANY_ARRAY) {
					if (ConjunctIndex >= 0) {
						left = operand->any_array.left.operands[ConjunctIndex];
						right = operand->any_array.left.operands[ConjunctIndex + 1];
					} else {
						left = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex];
						right = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex + 1];
					}
				}

				hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(left, right));
				if (dst == nullptr) {
					free_all(intersection);
					return (hol_term*) nullptr;
				}
				left->reference_count++;
				right->reference_count++;
				free_all(intersection);
				return dst;
			}, no_op()) && dst != nullptr;
}

inline hol_term* apply_predicate(
		hol_term* head,
		unsigned int head_variable,
		head_index predicate_index,
		hol_term* src_predicate,
		hol_term* dst_predicate)
{
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT
	 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)
	 || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_RIGHT))
	{
		hol_term* head_var = hol_term::new_variable(head_variable);
		if (head_var == nullptr) return (hol_term*) nullptr;
		constexpr unsigned int excluded_tree_count = 4;
		hol_term* excluded_trees[excluded_tree_count];
		excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
		excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
		excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
		excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
		if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
		if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
		if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
		if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
		if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
			if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
			if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
			if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
			free(*head_var); free(head_var);
			return (hol_term*) nullptr;
		}
		free(*head_var);

		hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
		if (conjunct == nullptr) {
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
			}
			return (hol_term*) nullptr;
		}

		hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
				make_array_view(&dst_predicate, 1), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
		if (expected_head == nullptr) {
			free(*conjunct); free(conjunct);
			return (hol_term*) nullptr;
		}
		dst_predicate->reference_count++;

		array<hol_term*> intersection(2);
		intersect<built_in_predicates>(intersection, head, expected_head);
		free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
		if (intersection.length == 0) {
			return nullptr;
		} else if (intersection.length != 1) {
			fprintf(stderr, "apply_predicate ERROR: Intersection is not unique.\n");
			free_all(intersection); return nullptr;
		}
		return intersection[0];

	} else if (head->type == hol_term_type::EXISTS) {
		hol_term* new_head;
		hol_term* operand = head->quantifier.operand;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			if (predicate_index.position == head_position::ANY) {
				if (src_predicate != &HOL_ANY && !has_intersection<built_in_predicates>(src_predicate, operand->any_array.any.operands[predicate_index.index]))
					return nullptr;
				new_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
						make_replaced_array_view(operand->any_array.any.operands, operand->any_array.any.length, dst_predicate, predicate_index.index),
						make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
						make_array_view(operand->any_array.right.operands, operand->any_array.right.length)));
			} else if (predicate_index.position == head_position::LEFT) {
				if (src_predicate != &HOL_ANY && !has_intersection<built_in_predicates>(src_predicate, operand->any_array.left.operands[predicate_index.index]))
					return nullptr;
				new_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
						make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
						make_replaced_array_view(operand->any_array.left.operands, operand->any_array.left.length, dst_predicate, predicate_index.index),
						make_array_view(operand->any_array.right.operands, operand->any_array.right.length)));
			} else if (predicate_index.position == head_position::RIGHT) {
				unsigned int index = operand->any_array.right.length - predicate_index.index - 1;
				if (src_predicate != &HOL_ANY && !has_intersection<built_in_predicates>(src_predicate, operand->any_array.right.operands[index]))
					return nullptr;
				new_head = hol_term::new_exists(head_variable, hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
						make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
						make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
						make_replaced_array_view(operand->any_array.right.operands, operand->any_array.right.length, dst_predicate, index)));
			} else {
				return nullptr;
			}
			if (new_head == nullptr)
				return nullptr;
			hol_term* new_operand = new_head->quantifier.operand;
			new_operand->any_array.all->reference_count++;
			for (unsigned int i = 0; i < new_operand->any_array.left.length; i++)
				new_operand->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_operand->any_array.right.length; i++)
				new_operand->any_array.right.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_operand->any_array.any.length; i++)
				new_operand->any_array.any.operands[i]->reference_count++;
			return new_head;
		} else if (operand->type == hol_term_type::AND) {
			if (predicate_index.position != head_position::LEFT)
				return nullptr;
			if (src_predicate != &HOL_ANY && !has_intersection<built_in_predicates>(src_predicate, operand->array.operands[predicate_index.index]))
				return nullptr;
			hol_term* new_head = hol_term::new_exists(head_variable, hol_term::new_and(make_replaced_array_view(operand->array.operands, operand->array.length, dst_predicate, predicate_index.index)));
			if (new_head == nullptr) return nullptr;
			hol_term* new_operand = new_head->quantifier.operand;
			for (unsigned int i = 0; i < new_operand->array.length; i++)
				new_operand->array.operands[i]->reference_count++;
			return new_head;
		} else {
			if (predicate_index.position != head_position::LEFT || predicate_index.index != 0)
				return nullptr;
			if (src_predicate != &HOL_ANY && !has_intersection<built_in_predicates>(src_predicate, operand))
				return nullptr;
			hol_term* new_head = hol_term::new_exists(head_variable, dst_predicate);
			if (new_head == nullptr) return nullptr;
			dst_predicate->reference_count++;
			return new_head;
		}
	} else {
		return nullptr;
	}
}

template<typename ApplyHead>
inline hol_term* apply_array(hol_term*& head, hol_term* parent, ApplyHead apply_head_function)
{
	bool head_is_array = false;
	if (parent != nullptr && (parent->type == hol_term_type::AND || parent->type == hol_term_type::OR)) {
		hol_term* head; head_index predicate_index;
		find_head<built_in_predicates>(parent->array.operands[0], head, predicate_index);
		if (head != nullptr)
			head_is_array = true;
	} else if (parent != nullptr && parent->type == hol_term_type::ANY_ARRAY && (parent->any_array.oper == hol_term_type::ANY_ARRAY || parent->any_array.oper == hol_term_type::AND || parent->any_array.oper == hol_term_type::OR)) {
		if (head == parent->any_array.all) {
			head_is_array = true;
		} else {
			hol_term* head; head_index predicate_index;
			find_head<built_in_predicates>(parent->any_array.all, head, predicate_index);
			if (head != nullptr)
				head_is_array = true;
		}
	}

	hol_term* new_head;
	if (head_is_array) {
		if (parent->type == hol_term_type::AND || parent->type == hol_term_type::OR) {
			array<hol_term*> new_heads(parent->array.length);
			for (unsigned int i = 0; i < parent->array.length; i++) {
				new_heads[i] = apply_head_function(parent->array.operands[i]);
				if (new_heads[i] == nullptr) {
					free_all(new_heads);
					return nullptr;
				}
				new_heads.length++;
			}
			if (parent->type == hol_term_type::AND)
				new_head = hol_term::new_and(make_array_view(new_heads.data, new_heads.length));
			else new_head = hol_term::new_or(make_array_view(new_heads.data, new_heads.length));
			if (new_head == nullptr) {
				free_all(new_heads);
				return nullptr;
			}
		} else {
			hol_term* new_all = apply_head_function(parent->any_array.all);
			if (new_all == nullptr)
				return nullptr;
			array<hol_term*> new_left(max(1u, parent->any_array.left.length));
			array<hol_term*> new_right(max(1u, parent->any_array.right.length));
			array<hol_term*> new_any(max(1u, parent->any_array.any.length));
			for (unsigned int i = 0; i < parent->any_array.left.length; i++) {
				new_left[i] = apply_head_function(parent->any_array.left.operands[i]);
				if (new_left[i] == nullptr) {
					free_all(new_left);
					free(*new_all); if (new_all->reference_count == 0) free(new_all);
					return nullptr;
				}
				new_left.length++;
			} for (unsigned int i = 0; i < parent->any_array.right.length; i++) {
				new_right[i] = apply_head_function(parent->any_array.right.operands[i]);
				if (new_right[i] == nullptr) {
					free_all(new_left); free_all(new_right);
					free(*new_all); if (new_all->reference_count == 0) free(new_all);
					return nullptr;
				}
				new_right.length++;
			} for (unsigned int i = 0; i < parent->any_array.any.length; i++) {
				new_any[i] = apply_head_function(parent->any_array.any.operands[i]);
				if (new_any[i] == nullptr) {
					free_all(new_left); free_all(new_right); free_all(new_any);
					free(*new_all); if (new_all->reference_count == 0) free(new_all);
					return nullptr;
				}
				new_any.length++;
			}
			new_head = hol_term::new_any_array(parent->any_array.oper, new_all,
					make_array_view(new_any.data, new_any.length),
					make_array_view(new_left.data, new_left.length),
					make_array_view(new_right.data, new_right.length));
			if (new_head == nullptr) {
				free_all(new_left); free_all(new_right); free_all(new_any);
				free(*new_all); if (new_all->reference_count == 0) free(new_all);
				return nullptr;
			}
		}
		head = parent;
	} else {
		new_head = apply_head_function(head);
	}
	return new_head;
}

inline hol_term* do_set_predicate_empty(hol_term* head)
{
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* new_predicate = hol_term::new_apply(&HOL_EMPTY, hol_term::new_variable(head_variable));
	if (new_predicate == nullptr)
		return nullptr;
	HOL_EMPTY.reference_count++;

	head_index predicate_index;
	find_predicate<built_in_predicates>(head_variable, head->quantifier.operand, predicate_index);
	hol_term* new_head = apply_predicate(head, head_variable, predicate_index, &HOL_ANY, new_predicate);
	free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);
	return new_head;
}

inline bool set_predicate_empty(
		hol_term* src, hol_term*& dst)
{
	hol_term* parent = nullptr;
	hol_term* current = nullptr;
	auto apply = [&parent,&current](hol_term* term) {
		parent = current;
		current = term;
	};

	head_index predicate_index;
	hol_term* head = find_head(src, predicate_index, find_head<built_in_predicates>, apply);
	if (head == nullptr)
		return false;

	hol_term* new_head = apply_array(head, parent, do_set_predicate_empty);
	if (new_head == nullptr)
		return false;

	dst = substitute_head<any_node_position::NONE>(src, head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	return (dst != nullptr);
}

template<hol_term_type Type, bool Negative, int_fast8_t Length = -1>
inline bool require_array(
		hol_term* src, hol_term*& dst)
{
	static_assert(Length == -1 || Length >= 2, "require_array ERROR: `Length` must be -1 or at least 2.");
	static_assert(Type == hol_term_type::AND || Type == hol_term_type::OR, "require_array ERROR: Unsupported value for `Type`.");

	head_index predicate_index; no_op apply;
	auto find_array_head = make_array_finder(find_head<built_in_predicates>);
	hol_term* head = find_head(src, predicate_index, find_array_head, apply);
	if (head == nullptr)
		return false;

	hol_term* expected_conjunct = (Negative ? hol_term::new_not(&HOL_ANY) : &HOL_ANY);
	if (expected_conjunct == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* expected_head;
	if (Length == -1) {
		expected_head = hol_term::new_any_array(Type, expected_conjunct,
				make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0));
		if (expected_head == nullptr) {
			free(*expected_conjunct); if (Negative && expected_conjunct->reference_count == 0) free(expected_conjunct);
			return false;
		}
	} else if (Type == hol_term_type::AND) {
		expected_head = hol_term::new_and(make_repeated_array_view(expected_conjunct, Length));
		if (expected_head == nullptr) {
			free(*expected_conjunct); if (Negative && expected_conjunct->reference_count == 0) free(expected_conjunct);
			return false;
		}
		expected_conjunct->reference_count += Length - 1;
	} else if (Type == hol_term_type::OR) {
		expected_head = hol_term::new_or(make_repeated_array_view(expected_conjunct, Length));
		if (expected_head == nullptr) {
			free(*expected_conjunct); if (Negative && expected_conjunct->reference_count == 0) free(expected_conjunct);
			return false;
		}
		expected_conjunct->reference_count += Length - 1;
	}

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, head, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "require_array ERROR: Set intersection is not unique.\n");
		free_all(intersection);
		return false;
	}
	hol_term* new_head = intersection[0];

	dst = substitute_head<any_node_position::NONE>(src, head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	return (dst != nullptr);
}

template<int_fast8_t ConjunctIndex>
inline bool select_operand(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, make_array_finder(find_root),
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					hol_term* excluded_array = hol_term::new_any_array(
							hol_term_type::ANY_ARRAY, &HOL_ANY, make_array_view((hol_term**) nullptr, 0),
							make_repeated_array_view(&HOL_ANY, 2), make_array_view((hol_term**) nullptr, 0));
					if (excluded_array == nullptr)
						return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					array<hol_term*> difference(2);
					intersect<built_in_predicates>(difference, head, excluded_array);
					free(*excluded_array); if (excluded_array->reference_count == 0) free(excluded_array);
					if (difference.length == 0) {
						return (hol_term*) nullptr;
					} else if (difference.length != 1) {
						fprintf(stderr, "select_operand ERROR: Set difference is not unique.\n");
						free_all(difference);
						return (hol_term*) nullptr;
					}
					dst = difference[0];

				} else if (head->type == hol_term_type::AND || head->type == hol_term_type::OR) {
					if (ConjunctIndex >= 0 && ConjunctIndex >= head->array.length)
						return (hol_term*) nullptr;
					if (ConjunctIndex < 0 && -ConjunctIndex > head->array.length)
						return (hol_term*) nullptr;

					unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (ConjunctIndex + head->array.length);
					dst = head->array.operands[index];
					dst->reference_count++;

				} else if (head->type == hol_term_type::ANY_ARRAY && (head->any_array.oper == hol_term_type::ANY_ARRAY || head->any_array.oper == hol_term_type::AND || head->any_array.oper == hol_term_type::OR)) {
					if (ConjunctIndex >= 0) {
						if (ConjunctIndex >= head->any_array.left.length) {
							dst = head->any_array.all;
							dst->reference_count++;
						} else {
							dst = head->any_array.left.operands[ConjunctIndex];
							dst->reference_count++;
						}
					} else {
						if (-ConjunctIndex > head->any_array.right.length) {
							dst = head->any_array.all;
							dst->reference_count++;
						} else {
							unsigned int index = ConjunctIndex + head->any_array.right.length;
							dst = head->any_array.right.operands[index];
							dst->reference_count++;
						}
					}

				} else {
					if (ConjunctIndex != 0 && ConjunctIndex != -1)
						return (hol_term*) nullptr;
					dst = head;
					dst->reference_count++;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool remove_operand(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, make_array_finder(find_root),
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					hol_term* expected_head = hol_term::new_any_array(
							hol_term_type::ANY_ARRAY, &HOL_ANY, make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0));
					if (expected_head == nullptr)
						return (hol_term*) nullptr;
					HOL_ANY.reference_count++;

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, head, expected_head);
					free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
					if (intersection.length == 0) {
						return (hol_term*) nullptr;
					} else if (intersection.length != 1) {
						fprintf(stderr, "remove_operand ERROR: Set intersection is not unique.\n");
						free_all(intersection);
						return (hol_term*) nullptr;
					}
					dst = intersection[0];

				} else if (head->type == hol_term_type::AND || head->type == hol_term_type::OR) {
					if (ConjunctIndex >= 0 && ConjunctIndex >= head->array.length)
						return (hol_term*) nullptr;
					if (ConjunctIndex < 0 && -ConjunctIndex > head->array.length)
						return (hol_term*) nullptr;

					unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (ConjunctIndex + head->array.length);
					if (head->array.length == 2)
						dst = (index == 0 ? head->array.operands[1] : head->array.operands[0]);
					else if (head->type == hol_term_type::AND)
						dst = hol_term::new_and(make_excluded_array_view(head->array.operands, head->array.length, index));
					else dst = hol_term::new_or(make_excluded_array_view(head->array.operands, head->array.length, index));

					if (dst == nullptr)
						return (hol_term*) nullptr;
					for (unsigned int i = 0; i < head->array.length; i++)
						if (i != index) head->array.operands[i]->reference_count++;

				} else if (head->type == hol_term_type::ANY_ARRAY && (head->any_array.oper == hol_term_type::ANY_ARRAY || head->any_array.oper == hol_term_type::AND || head->any_array.oper == hol_term_type::OR)) {
					if (ConjunctIndex >= 0) {
						if (ConjunctIndex >= head->any_array.left.length) {
							dst = head;
							dst->reference_count++;
						} else {
							dst = hol_term::new_any_array(head->any_array.oper, head->any_array.all,
									make_array_view(head->any_array.any.operands, head->any_array.any.length),
									make_excluded_array_view(head->any_array.left.operands, head->any_array.left.length, ConjunctIndex),
									make_array_view(head->any_array.right.operands, head->any_array.right.length));
							if (dst == nullptr)
								return (hol_term*) nullptr;
							dst->any_array.all->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.left.length; i++)
								dst->any_array.left.operands[i]->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.right.length; i++)
								dst->any_array.right.operands[i]->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.any.length; i++)
								dst->any_array.any.operands[i]->reference_count++;

							if (dst->any_array.left.length <= 1 && dst->any_array.right.length <= 1 && dst->any_array.any.length <= 1) {
								if (dst->any_array.right.length == 1 && is_subset<built_in_predicates>(dst, dst->any_array.right.operands[0])) {
									hol_term* new_dst = dst->any_array.right.operands[0];
									new_dst->reference_count++;
									free(*dst); free(dst);
									dst = new_dst;
								} else if (dst->any_array.left.length == 1 && is_subset<built_in_predicates>(dst, dst->any_array.left.operands[0])) {
									hol_term* new_dst = dst->any_array.left.operands[0];
									new_dst->reference_count++;
									free(*dst); free(dst);
									dst = new_dst;
								} else if (dst->any_array.any.length == 1 && is_subset<built_in_predicates>(dst, dst->any_array.any.operands[0])) {
									hol_term* new_dst = dst->any_array.any.operands[0];
									new_dst->reference_count++;
									free(*dst); free(dst);
									dst = new_dst;
								} else if (dst->any_array.right.length == 0 && dst->any_array.left.length == 0 && dst->any_array.any.length == 0
										&& is_subset<built_in_predicates>(dst, dst->any_array.all))
								{
									hol_term* new_dst = dst->any_array.all;
									new_dst->reference_count++;
									free(*dst); free(dst);
									dst = new_dst;
								}
							}
						}
					} else {
						if (-ConjunctIndex > head->any_array.right.length) {
							dst = head;
							dst->reference_count++;
						} else {
							unsigned int index = ConjunctIndex + head->any_array.right.length;
							dst = hol_term::new_any_array(head->any_array.oper, head->any_array.all,
									make_array_view(head->any_array.any.operands, head->any_array.any.length),
									make_array_view(head->any_array.left.operands, head->any_array.left.length),
									make_excluded_array_view(head->any_array.right.operands, head->any_array.right.length, index));
							if (dst == nullptr)
								return (hol_term*) nullptr;
							dst->any_array.all->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.left.length; i++)
								dst->any_array.left.operands[i]->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.right.length; i++)
								dst->any_array.right.operands[i]->reference_count++;
							for (unsigned int i = 0; i < dst->any_array.any.length; i++)
								dst->any_array.any.operands[i]->reference_count++;
						}
					}

				} else {
					return (hol_term*) nullptr;
				}
				return dst;
			}, no_op()) && dst != nullptr;
}

bool apply_to_predicative_set_function(hol_term* head,
		array<hol_term*>& dst, unsigned int& max_variable,
		unsigned int src_predicate, unsigned int dst_predicate)
{
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int set_variable;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		set_variable = ++max_variable;
	} else {
#if !defined(NDEBUG)
		if (head->type != hol_term_type::EXISTS) {
			fprintf(stderr, "apply_to_predicative_set_function ERROR: Expected existential quantification of set.\n");
			return false;
		}
#endif
		set_variable = head->quantifier.variable;
	}

	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
	if (excluded_quantifier == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* expected_head;
	if (src_predicate == (unsigned int) built_in_predicates::EQUALS) {
		expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_any(nullptr, &excluded_quantifier, 1)),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
	} else {
		expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
				hol_term::new_apply(hol_term::new_constant(src_predicate), hol_term::new_variable(set_variable), hol_term::new_any(nullptr, &excluded_quantifier, 1)),
				hol_term::new_any(nullptr, &excluded_quantifier, 1)));
	}
	if (expected_head == nullptr) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}
	excluded_quantifier->reference_count += 2 - 1;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, expected_head, head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (intersection.length == 0)
		return false;

	if (!dst.ensure_capacity(dst.length + intersection.length)) {
		free_all(intersection);
		return false;
	}

	for (hol_term* term : intersection) {
		hol_term* set_var; hol_term* set_definition;
		if (src_predicate == (unsigned int) built_in_predicates::EQUALS) {
			set_var = term->quantifier.operand->array.operands[0]->binary.left;
			set_definition = term->quantifier.operand->array.operands[0]->binary.right;
		} else {
			set_var = term->quantifier.operand->array.operands[0]->ternary.second;
			set_definition = term->quantifier.operand->array.operands[0]->ternary.third;
		}
		hol_term* right = term->quantifier.operand->array.operands[1];
		if (dst_predicate == (unsigned int) built_in_predicates::EQUALS) {
			dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_and(
					hol_term::new_equals(set_var, set_definition), right));
		} else {
			dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_and(
					hol_term::new_apply(hol_term::new_constant(dst_predicate), set_var, set_definition), right));
		}
		if (dst[dst.length] == nullptr) {
			free_all(intersection); free_all(dst);
			return false;
		}
		set_var->reference_count++;
		set_definition->reference_count++;
		right->reference_count++;
		dst.length++;
	}
	free_all(intersection);
	return true;
}

inline bool apply_to_predicative_set_function(
		hol_term* src, hol_term*& dst,
		unsigned int src_predicate,
		unsigned int dst_predicate)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	return apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker, predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable,src_predicate,dst_predicate](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				array<hol_term*> dst(2);
				if (!apply_to_predicative_set_function(head, dst, max_variable, src_predicate, dst_predicate))
					return (hol_term*) nullptr;
				if (dst.length != 1) {
					fprintf(stderr, "apply_to_predicative_set_function ERROR: Intersection is not unique.\n");
					free_all(dst); return (hol_term*) nullptr;
				}
				return dst[0];
			}, no_op()) && dst != nullptr;
}

template<typename Formula>
inline bool add_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE
	 && src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags.flags[(uint_fast8_t) flag] = grammatical_flag_value::TRUE;
	return true;
}

template<typename Formula>
inline bool remove_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE
	 && src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool try_remove_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	dst = src;
	dst.flags.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool require_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE
	 && src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags.flags[(uint_fast8_t) flag] = grammatical_flag_value::TRUE;
	return true;
}

template<typename Formula>
inline bool require_no_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE
	 && src.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool add_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.flags.cnj != grammatical_conjunction::NONE && src.flags.cnj != grammatical_conjunction::ANY)
		return false;
	dst = src;
	dst.flags.cnj = cnj;
	return true;
}

template<typename Formula>
inline bool remove_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.flags.cnj != cnj && src.flags.cnj != grammatical_conjunction::ANY)
		return false;
	dst = src;
	dst.flags.cnj = grammatical_conjunction::NONE;
	return true;
}

template<typename Formula>
inline bool require_no_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.flags.cnj == cnj)
		return false;
	dst = src;
	dst.flags.cnj = src.flags.cnj;
	return true;
}

template<typename Formula>
inline bool add_correlator(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator corr)
{
	if (src.flags.corr != correlator::NONE && src.flags.corr != correlator::ANY)
		return false;
	dst = src;
	dst.flags.corr = corr;
	return true;
}

template<typename Formula>
inline bool remove_correlator(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator corr)
{
	if (src.flags.corr != corr && src.flags.corr != correlator::ANY)
		return false;
	dst = src;
	dst.flags.corr = correlator::NONE;
	return true;
}

template<typename Formula>
inline bool add_correlated_by(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator correlated_by)
{
	if (src.flags.correlated_by != correlator::NONE && src.flags.correlated_by != correlator::ANY)
		return false;
	dst = src;
	dst.flags.correlated_by = correlated_by;
	return true;
}

template<typename Formula>
inline bool add_coordination(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		coordination coord)
{
	if (src.flags.coord != coordination::NONE && src.flags.coord != coordination::ANY)
		return false;
	dst = src;
	dst.flags.coord = coord;
	return true;
}

template<typename Formula>
inline bool remove_coordination(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		coordination coord)
{
	if (!has_intersection(src.flags.coord, coord))
		return false;
	dst = src;
	dst.flags.coord = coordination::NONE;
	return true;
}

template<typename Formula>
inline bool apply_auxiliary(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		const auxiliary_flag expected_aux,
		auxiliary_flag dst_aux)
{
	if (!has_intersection(src.flags.aux, expected_aux))
		return false;
	dst = src;
	dst.flags.aux = dst_aux;
	return true;
}

template<typename Formula>
inline bool require_auxiliary(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		auxiliary_flag aux)
{
	auxiliary_flag dst_aux;
	if (!intersect(dst_aux, src.flags.aux, aux))
		return false;
	dst = src;
	dst.flags.aux = dst_aux;
	return true;
}

template<typename Formula>
inline bool add_mood(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_mood mood)
{
	if (!has_intersection(grammatical_mood::INDICATIVE, src.flags.mood))
		return false;
	dst = src;
	dst.flags.mood = mood;
	return true;
}

template<typename Formula>
inline bool remove_mood(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_mood mood)
{
	if (!has_intersection(mood, src.flags.mood))
		return false;
	dst = src;
	dst.flags.mood = grammatical_mood::INDICATIVE;
	return true;
}

template<typename Formula>
inline bool try_remove_mood(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst)
{
	dst = src;
	dst.flags.mood = grammatical_mood::INDICATIVE;
	return true;
}

template<typename Formula>
inline bool require_mood(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_mood mood)
{
	if (!has_intersection(mood, src.flags.mood))
		return false;
	dst = src;
	dst.flags.mood = mood;
	return true;
}

template<typename Formula>
inline bool require_no_to_infinitive(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst)
{
	grammatical_mood intersection;
	if (!intersect(intersection, src.flags.mood, grammatical_mood::NOT_TO_INFINITIVE))
		return false;
	dst = src;
	dst.flags.mood = intersection;
	return true;
}

template<typename Formula>
inline bool require_no_subjunctive(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst)
{
	grammatical_mood intersection;
	if (!intersect(intersection, src.flags.mood, grammatical_mood::NOT_SUBJUNCTIVE))
		return false;
	dst = src;
	dst.flags.mood = intersection;
	return true;
}

template<typename Formula>
inline bool require_aux_or_subjunctive_or_infinitive_or_to_infinitive(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst)
{
	if (src.flags.aux == auxiliary_flag::AUX || src.flags.mood == grammatical_mood::SUBJUNCTIVE
	 || src.flags.mood == grammatical_mood::BARE_INFINITIVE || src.flags.mood == grammatical_mood::TO_INFINITIVE) {
		dst = src;
		return true;
	} else if (!has_intersection(src.flags.aux, auxiliary_flag::AUX) && !has_intersection(src.flags.mood, grammatical_mood::SUBJUNCTIVE)
			&& !has_intersection(src.flags.mood, grammatical_mood::TO_INFINITIVE) && !has_intersection(src.flags.mood, grammatical_mood::BARE_INFINITIVE))
	{
		return false;
	} else {
		dst = src;
		dst.flags.aux_or_subjunctive_or_inf_or_to_inf = true;
		return true;
	}
}


/* forward declarations for semantic feature functions */

template<int_fast8_t ConjunctIndex>
bool set_arg(hol_term*, hol_term*&, unsigned int);


template<typename Formula>
bool apply(typename flagged_logical_form<Formula>::function function,
		const flagged_logical_form<Formula>& src, flagged_logical_form<Formula>& dst)
{
	typedef typename flagged_logical_form<Formula>::function_type function_type;

	uint_fast8_t i;
	switch (function.type) {
	case function_type::EMPTY:
		dst.root = Formula::new_false();
		dst.flags.index_number = grammatical_num::NONE;
		dst.flags.concord_number = grammatical_num::NONE;
		dst.flags.comp = grammatical_comparison::NONE;
		dst.flags.cnj = grammatical_conjunction::NONE;
		dst.flags.corr = correlator::NONE;
		dst.flags.correlated_by = correlator::NONE;
		dst.flags.coord = coordination::NONE;
		dst.flags.aux = auxiliary_flag::NONE;
		dst.flags.mood = grammatical_mood::INDICATIVE;
		dst.flags.aux_or_subjunctive_or_inf_or_to_inf = false;
		dst.flags.is_first_token_capital = false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			dst.flags.flags[i] = grammatical_flag_value::FALSE;
		return dst.root != nullptr;
	case function_type::IDENTITY:
		dst = src;
		return true;
	case function_type::SELECT_RIGHT_CONJUNCT:
		dst.flags = src.flags;
		return select_conjunct<-1>(src.root, dst.root);
	case function_type::REMOVE_RIGHT_CONJUNCT:
		dst.flags = src.flags;
		return remove_conjunct<-1>(src.root, dst.root);
	case function_type::SELECT_LEFT_CONJUNCT:
		dst.flags = src.flags;
		return select_conjunct<0>(src.root, dst.root);
	case function_type::SELECT_LEFT_CONJUNCT_AND_NEGATION:
		dst.flags = src.flags;
		return select_conjunct<0, true>(src.root, dst.root);
	case function_type::REMOVE_LEFT_CONJUNCT:
		dst.flags = src.flags;
		return remove_conjunct<0>(src.root, dst.root);
	case function_type::REMOVE_LEFT_CONJUNCT_AND_NEGATION:
		dst.flags = src.flags;
		return remove_conjunct<0, true>(src.root, dst.root);
	case function_type::REMOVE_SECOND_LEFT_CONJUNCT:
		dst.flags = src.flags;
		return remove_conjunct<1>(src.root, dst.root);
	case function_type::SELECT_SECOND_LEFT_SET_CONJUNCT:
		dst.flags = src.flags;
		return select_set_conjunct<1>(src.root, dst.root);
	case function_type::SELECT_SECOND_LEFT_SET_CONJUNCT_ROOT:
		dst.flags = src.flags;
		return select_set_conjunct<1, true>(src.root, dst.root);
	case function_type::REMOVE_SECOND_LEFT_SET_CONJUNCT:
		dst.flags = src.flags;
		return remove_set_conjunct<1>(src.root, dst.root);
	case function_type::SELECT_LEFT_CONJUNCT_IN_SET:
		dst.flags = src.flags;
		return select_conjunct_in_set<0>(src.root, dst.root);
	case function_type::REMOVE_LEFT_CONJUNCT_IN_SET:
		dst.flags = src.flags;
		return remove_conjunct_in_set<0>(src.root, dst.root);
	case function_type::SELECT_RIGHT_CONJUNCT_IN_SET:
		dst.flags = src.flags;
		return select_conjunct_in_set<-1>(src.root, dst.root);
	case function_type::REMOVE_RIGHT_CONJUNCT_IN_SET:
		dst.flags = src.flags;
		return remove_conjunct_in_set<-1>(src.root, dst.root);
	case function_type::SELECT_RIGHT_SUBSET_IN_SET:
		dst.flags = src.flags;
		return select_subset_in_set<-1>(src.root, dst.root);
	case function_type::REQUIRE_NO_INVERSE:
		dst.flags = src.flags;
		return require_no_inverse(src.root, dst.root);
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE:
		dst.flags = src.flags;
		return require_left_predicate_inverse(src.root, dst.root);
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN:
		dst.flags = src.flags;
		return require_left_predicate_inverse_own(src.root, dst.root);
	case function_type::REQUIRE_LEFT_PREDICATE_EXIST:
		dst.flags = src.flags;
		return require_left_predicate_exist(src.root, dst.root);
	case function_type::REQUIRE_NO_LEFT_PREDICATE_EXIST:
		dst.flags = src.flags;
		return require_no_left_predicate_exist(src.root, dst.root);
	case function_type::REMOVE_INVERSE:
		dst.flags = src.flags;
		return remove_inverse(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG1, false>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG2, false>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG3, false>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG1_OF_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG1_OF, true>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG2_OF_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG2_OF, true>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG3_OF_WITHOUT_HEAD:
		dst.flags = src.flags;
		return select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG3_OF, true>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG1>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG2>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG3>(src.root, dst.root);
	case function_type::SELECT_SINGLETON_ARG1_IN_SET_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_singleton_arg_in_set_without_head_predicative<(unsigned int) built_in_predicates::ARG1>(src.root, dst.root);
	case function_type::SELECT_SINGLETON_ARG2_IN_SET_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_singleton_arg_in_set_without_head_predicative<(unsigned int) built_in_predicates::ARG2>(src.root, dst.root);
	case function_type::SELECT_SINGLETON_ARG3_IN_SET_WITHOUT_HEAD_PREDICATIVE:
		dst.flags = src.flags;
		return select_singleton_arg_in_set_without_head_predicative<(unsigned int) built_in_predicates::ARG3>(src.root, dst.root);
	case function_type::REMOVE_FUTURE:
		dst.flags = src.flags;
		return apply_tense_predicate(src.root, dst.root, FUTURE_PREDICATES, PRESENT_PREDICATES);
	case function_type::REQUIRE_NO_FUTURE:
		dst.flags = src.flags;
		return require_no_tense_predicates(src.root, dst.root, FUTURE_PREDICATES);
	case function_type::REMOVE_PERFECT:
		dst.flags = src.flags;
		return apply_tense_predicate(src.root, dst.root, PERFECT_PREDICATES, NON_PERFECT_PREDICATES);
	case function_type::REQUIRE_NO_PERFECT:
		dst.flags = src.flags;
		return require_no_tense_predicates(src.root, dst.root, PERFECT_PREDICATES);
	case function_type::REMOVE_PROGRESSIVE:
		dst.flags = src.flags;
		return apply_tense_predicate(src.root, dst.root, PROGRESSIVE_PREDICATES, NON_PROGRESSIVE_PREDICATES);
	case function_type::REQUIRE_NO_PROGRESSIVE:
		dst.flags = src.flags;
		return require_no_tense_predicates(src.root, dst.root, PROGRESSIVE_PREDICATES);
	case function_type::REMOVE_NOT:
		dst.flags = src.flags;
		return remove_not(src.root, dst.root);
	case function_type::REQUIRE_NO_EMPTY_REF:
		dst.flags = src.flags;
		return require_no_empty_ref(src.root, dst.root);
	case function_type::REQUIRE_PREDICATIVE_UNIVERSAL:
		dst.flags = src.flags;
		return require_predicative_quantifier<hol_term_type::FOR_ALL>(src.root, dst.root);
	case function_type::REQUIRE_PREDICATIVE_EXISTENTIAL:
		dst.flags = src.flags;
		return require_predicative_quantifier<hol_term_type::EXISTS>(src.root, dst.root);
	case function_type::REMOVE_PREDICATIVE_NOT:
		dst.flags = src.flags;
		return remove_predicative_not(src.root, dst.root);
	case function_type::PREDICATE_ONLY:
		dst.flags = src.flags;
		return predicate_only(src.root, dst.root);
	case function_type::PREDICATE:
		dst.flags = src.flags;
		return predicate(src.root, dst.root);
	case function_type::PREDICATE_AND_TENSE:
		dst.flags = src.flags;
		return predicate_and_tense(src.root, dst.root);
	case function_type::EMPTY_AND_TENSE:
		dst.flags = src.flags;
		return empty_and_tense(src.root, dst.root);
	case function_type::SELECT_PREDICATE_IN_SET:
		dst.flags = src.flags;
		return select_predicate_in_set(src.root, dst.root);
	case function_type::MARK_WIDE_SCOPE:
		dst.flags = src.flags;
		return mark_wide_scope(src.root, dst.root);
	case function_type::REQUIRE_WIDE_SCOPE:
		dst.flags = src.flags;
		return require_narrow_or_wide_scope<true>(src.root, dst.root);
	case function_type::REQUIRE_NARROW_SCOPE:
		dst.flags = src.flags;
		return require_narrow_or_wide_scope<false>(src.root, dst.root);
	case function_type::REMOVE_WIDE_SCOPE:
		dst.flags = src.flags;
		return remove_wide_scope(src.root, dst.root);
	case function_type::REQUIRE_CONSTANT_IN_SET:
		dst.flags = src.flags;
		return require_constant_in_set(src.root, dst.root);
	case function_type::REQUIRE_NO_CONSTANT_IN_SET:
		dst.flags = src.flags;
		return require_no_constant_in_set(src.root, dst.root);
	case function_type::REQUIRE_SINGLETON:
		dst.flags = src.flags;
		return require_singleton(src.root, dst.root);
	case function_type::SIZE:
		dst.flags = src.flags;
		return size(src.root, dst.root, dst.flags);
	case function_type::SET_SIZE:
		dst.flags = src.flags;
		return set_size(src.root, dst.root);
	case function_type::REQUIRE_LEFT_ARG1:
		dst.flags = src.flags;
		return set_arg<0>(src.root, dst.root, (unsigned int) built_in_predicates::ARG1);
	case function_type::REQUIRE_LAMBDA:
		dst.flags = src.flags;
		return require_lambda(src.root, dst.root);
	case function_type::REQUIRE_NO_LAMBDA:
		dst.flags = src.flags;
		return require_no_lambda(src.root, dst.root);
	case function_type::REMOVE_LEFT_PREDICATE:
		dst.flags = src.flags;
		return remove_predicate<0>(src.root, dst.root);
	case function_type::FACTOR:
		dst.flags = src.flags;
		return factor(src.root, dst.root);
	case function_type::FACTOR_PREDICATIVE:
		dst.flags = src.flags;
		return factor_predicative(src.root, dst.root);
	case function_type::SELECT_LEFT_PREDICATE_AND_TENSE:
		dst.flags = src.flags;
		return select_predicate_and_tense<0>(src.root, dst.root);
	case function_type::SET_PREDICATE_EMPTY:
		dst.flags = src.flags;
		return set_predicate_empty(src.root, dst.root);
	case function_type::REQUIRE_CONJUNCTION:
		dst.flags = src.flags;
		return require_array<hol_term_type::AND, false>(src.root, dst.root);
	case function_type::REQUIRE_BINARY_CONJUNCTION:
		dst.flags = src.flags;
		return require_array<hol_term_type::AND, false, 2>(src.root, dst.root);
	case function_type::REQUIRE_DISJUNCTION:
		dst.flags = src.flags;
		return require_array<hol_term_type::OR, false>(src.root, dst.root);
	case function_type::REQUIRE_NEGATIVE_CONJUNCTION:
		dst.flags = src.flags;
		return require_array<hol_term_type::AND, true>(src.root, dst.root);
	case function_type::SELECT_LEFT_OPERAND:
		dst.flags = src.flags;
		return select_operand<0>(src.root, dst.root);
	case function_type::REMOVE_LEFT_OPERAND:
		dst.flags = src.flags;
		return remove_operand<0>(src.root, dst.root);
	case function_type::REPLACE_PREDICATIVE_SUBSET_WITH_EQUALITY:
		dst.flags = src.flags;
		return apply_to_predicative_set_function(src.root, dst.root, (unsigned int) built_in_predicates::SUBSET, (unsigned int) built_in_predicates::EQUALS);
	case function_type::ADD_SINGULAR:
		if (src.flags.index_number != grammatical_num::NONE && src.flags.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.index_number = grammatical_num::SINGULAR;
		return true;
	case function_type::ADD_PLURAL:
		if (src.flags.index_number != grammatical_num::NONE && src.flags.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.index_number = grammatical_num::PLURAL;
		return true;
	case function_type::REQUIRE_SINGULAR:
		if (src.flags.index_number != grammatical_num::SINGULAR && src.flags.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.index_number = grammatical_num::SINGULAR;
		return true;
	case function_type::REQUIRE_PLURAL:
		if (src.flags.index_number != grammatical_num::PLURAL && src.flags.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.index_number = grammatical_num::PLURAL;
		return true;
	case function_type::TRY_REMOVE_NUMBER:
		dst = src;
		dst.flags.index_number = grammatical_num::NONE;
		return true;
	case function_type::ADD_CONCORD_SINGULAR:
		if (src.flags.concord_number != grammatical_num::NONE && src.flags.concord_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.concord_number = grammatical_num::SINGULAR;
		return true;
	case function_type::ADD_CONCORD_PLURAL:
		if (src.flags.concord_number != grammatical_num::NONE && src.flags.concord_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.flags.concord_number = grammatical_num::PLURAL;
		return true;
	case function_type::ADD_THAT:
		return add_conjunction(src, dst, grammatical_conjunction::THAT);
	case function_type::REMOVE_THAT:
		return remove_conjunction(src, dst, grammatical_conjunction::THAT);
	case function_type::REQUIRE_NO_THAT:
		return require_no_conjunction(src, dst, grammatical_conjunction::THAT);
	case function_type::ADD_WHETHER:
		return add_conjunction(src, dst, grammatical_conjunction::WHETHER);
	case function_type::REMOVE_WHETHER:
		return remove_conjunction(src, dst, grammatical_conjunction::WHETHER);
	case function_type::ADD_IF:
		return add_conjunction(src, dst, grammatical_conjunction::IF);
	case function_type::REMOVE_IF:
		return remove_conjunction(src, dst, grammatical_conjunction::IF);
	case function_type::ADD_BECAUSE:
		return add_conjunction(src, dst, grammatical_conjunction::BECAUSE);
	case function_type::REMOVE_BECAUSE:
		return remove_conjunction(src, dst, grammatical_conjunction::BECAUSE);
	case function_type::ADD_FOR:
		return add_conjunction(src, dst, grammatical_conjunction::FOR);
	case function_type::REMOVE_FOR:
		return remove_conjunction(src, dst, grammatical_conjunction::FOR);
	case function_type::REQUIRE_NO_CONJUNCTION:
		if (src.flags.cnj != grammatical_conjunction::NONE && src.flags.cnj != grammatical_conjunction::ANY)
			return false;
		dst = src;
		dst.flags.cnj = grammatical_conjunction::NONE;
		return true;
	case function_type::ADD_IS_ADJUNCT:
		return add_flag(src, dst, grammatical_flag::IS_ADJUNCT);
	case function_type::TRY_REMOVE_IS_ADJUNCT:
		return try_remove_flag(src, dst, grammatical_flag::IS_ADJUNCT);
	case function_type::REQUIRE_NOT_ADJUNCT:
		return require_no_flag(src, dst, grammatical_flag::IS_ADJUNCT);
	case function_type::ADD_NULLABLE_SUBJECT:
		return add_flag(src, dst, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::REMOVE_NULLABLE_SUBJECT:
		return remove_flag(src, dst, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::TRY_REMOVE_NULLABLE_SUBJECT:
		return try_remove_flag(src, dst, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::ADD_SUBORDINATE:
		return add_flag(src, dst, grammatical_flag::SUBORDINATE);
	case function_type::REMOVE_SUBORDINATE:
		return remove_flag(src, dst, grammatical_flag::SUBORDINATE);
	case function_type::TRY_REMOVE_SUBORDINATE:
		return try_remove_flag(src, dst, grammatical_flag::SUBORDINATE);
	case function_type::REQUIRE_NO_SUBORDINATE:
		return require_no_flag(src, dst, grammatical_flag::SUBORDINATE);
	case function_type::ADD_PREPOSITION:
		return add_flag(src, dst, grammatical_flag::PREPOSITION);
	case function_type::REQUIRE_PREPOSITION:
		return require_flag(src, dst, grammatical_flag::PREPOSITION);
	case function_type::REQUIRE_NO_PREPOSITION:
		return require_no_flag(src, dst, grammatical_flag::PREPOSITION);
	case function_type::ADD_PARTICLE:
		return add_flag(src, dst, grammatical_flag::PARTICLE);
	case function_type::ADD_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NONE_OR_REQ_AUX, auxiliary_flag::AUX);
	case function_type::TRY_REMOVE_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NO_REQ_AUX, auxiliary_flag::NONE);
	case function_type::ADD_REQ_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NONE_OR_AUX, auxiliary_flag::REQ_AUX);
	case function_type::TRY_ADD_REQ_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NO_REQ_NO_AUX, auxiliary_flag::REQ_AUX);
	case function_type::REQUIRE_NO_REQ_AUX:
		return require_auxiliary(src, dst, auxiliary_flag::NO_REQ_AUX);
	case function_type::TRY_REMOVE_REQ_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NO_REQ_NO_AUX, auxiliary_flag::NONE);
	case function_type::ADD_REQ_NO_AUX:
		return apply_auxiliary(src, dst, auxiliary_flag::NONE_OR_AUX, auxiliary_flag::REQ_NO_AUX);
	case function_type::REQUIRE_NO_REQ_NO_AUX:
		return require_auxiliary(src, dst, auxiliary_flag::NO_REQ_NO_AUX);
	case function_type::ADD_INFINITIVE:
		return add_mood(src, dst, grammatical_mood::BARE_INFINITIVE);
	case function_type::ADD_TO_INFINITIVE:
		return add_mood(src, dst, grammatical_mood::TO_INFINITIVE);
	case function_type::REMOVE_TO_INFINITIVE:
		return remove_mood(src, dst, grammatical_mood::TO_INFINITIVE);
	case function_type::REQUIRE_TO_INFINITIVE:
		return require_mood(src, dst, grammatical_mood::TO_INFINITIVE);
	case function_type::REQUIRE_NO_TO_INFINITIVE:
		return require_no_to_infinitive(src, dst);
	case function_type::ADD_SUBJUNCTIVE:
		return add_mood(src, dst, grammatical_mood::SUBJUNCTIVE);
	case function_type::REQUIRE_NO_SUBJUNCTIVE:
		return require_no_subjunctive(src, dst);
	case function_type::REQUIRE_AUX_OR_SUBJUNCTIVE_OR_INFINITIVE_OR_TO_INFINITIVE:
		return require_aux_or_subjunctive_or_infinitive_or_to_infinitive(src, dst);
	case function_type::ADD_BOTH:
		return add_correlator(src, dst, correlator::BOTH);
	case function_type::ADD_EITHER:
		return add_correlator(src, dst, correlator::EITHER);
	case function_type::ADD_NEITHER:
		return add_correlator(src, dst, correlator::NEITHER);
	case function_type::REMOVE_BOTH:
		return remove_correlator(src, dst, correlator::BOTH);
	case function_type::REMOVE_EITHER:
		return remove_correlator(src, dst, correlator::EITHER);
	case function_type::REMOVE_NEITHER:
		return remove_correlator(src, dst, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATOR:
		dst = src;
		dst.flags.corr = correlator::NONE;
		return true;
	case function_type::REQUIRE_NO_CORRELATOR:
		if (src.flags.corr != correlator::NONE && src.flags.corr != correlator::ANY)
			return false;
		dst = src;
		dst.flags.corr = correlator::NONE;
		return true;
	case function_type::ADD_CORRELATED_BY_BOTH:
		return add_correlated_by(src, dst, correlator::BOTH);
	case function_type::ADD_CORRELATED_BY_EITHER:
		return add_correlated_by(src, dst, correlator::EITHER);
	case function_type::ADD_CORRELATED_BY_NEITHER:
		return add_correlated_by(src, dst, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATED:
		dst = src;
		dst.flags.correlated_by = correlator::NONE;
		return true;
	case function_type::REQUIRE_NOT_CORRELATED:
		if (src.flags.correlated_by != correlator::NONE && src.flags.correlated_by != correlator::ANY)
			return false;
		dst = src;
		dst.flags.correlated_by = correlator::NONE;
		return true;
	case function_type::ADD_PAST_PARTICIPLE:
		return add_mood(src, dst, grammatical_mood::PAST_PARTICIPLE);
	case function_type::ADD_PRESENT_PARTICIPLE:
		return add_mood(src, dst, grammatical_mood::PRESENT_PARTICIPLE);
	case function_type::TRY_REMOVE_PARTICIPLE:
		return try_remove_mood(src, dst);
	case function_type::REQUIRE_PAST_PARTICIPLE:
		return require_mood(src, dst, grammatical_mood::PAST_PARTICIPLE);
	case function_type::REQUIRE_PRESENT_PARTICIPLE:
		return require_mood(src, dst, grammatical_mood::PRESENT_PARTICIPLE);
	case function_type::ADD_NEGATIVE:
		return add_flag(src, dst, grammatical_flag::NEGATIVE);
	case function_type::REQUIRE_NEGATIVE:
		return require_flag(src, dst, grammatical_flag::NEGATIVE);
	case function_type::ADD_ADV:
		return add_flag(src, dst, grammatical_flag::ADV);
	case function_type::REMOVE_ADV:
		return remove_flag(src, dst, grammatical_flag::ADV);
	case function_type::TRY_REMOVE_ADV:
		return try_remove_flag(src, dst, grammatical_flag::ADV);
	case function_type::REQUIRE_NO_ADV:
		return require_no_flag(src, dst, grammatical_flag::ADV);
	case function_type::ADD_TION:
		return add_flag(src, dst, grammatical_flag::TION);
	case function_type::ADD_LY:
		return add_flag(src, dst, grammatical_flag::LY);
	case function_type::ADD_GENITIVE:
		return add_flag(src, dst, grammatical_flag::GENITIVE);
	case function_type::TRY_REMOVE_GENITIVE:
		return try_remove_flag(src, dst, grammatical_flag::GENITIVE);
	case function_type::REQUIRE_NO_GENITIVE:
		return require_no_flag(src, dst, grammatical_flag::GENITIVE);
	case function_type::ADD_COMMA:
		return add_flag(src, dst, grammatical_flag::COMMA);
	case function_type::REMOVE_COMMA:
		return remove_flag(src, dst, grammatical_flag::COMMA);
	case function_type::REQUIRE_NO_COMMA:
		return require_no_flag(src, dst, grammatical_flag::COMMA);
	case function_type::ADD_AND:
		return add_coordination(src, dst, coordination::AND);
	case function_type::ADD_OR:
		return add_coordination(src, dst, coordination::OR);
	case function_type::ADD_NOR:
		return add_coordination(src, dst, coordination::NOR);
	case function_type::REMOVE_AND:
		return remove_coordination(src, dst, coordination::AND);
	case function_type::REMOVE_OR:
		return remove_coordination(src, dst, coordination::OR);
	case function_type::REMOVE_NOR:
		return remove_coordination(src, dst, coordination::NOR);
	case function_type::REMOVE_COORD:
		return remove_coordination(src, dst, coordination::NOT_NONE);
	}
	fprintf(stderr, "apply ERROR: Unrecognized transformation function.\n");
	return false;
}

template<
	typename FindFirstHeadFunction,
	typename FindSecondHeadFunction,
	typename OnRemapFunction>
inline hol_term* remap_to_invert_apply_head(
		hol_term* first, hol_term* second,
		FindFirstHeadFunction find_first_head,
		FindSecondHeadFunction find_second_head,
		OnRemapFunction on_remap_variables,
		unsigned int& max_variable,
		bool& first_head_is_array,
		bool& second_head_is_array)
{
	apply_head_inverter first_inverter; head_index first_predicate_index;
	apply_head_inverter second_inverter; head_index second_predicate_index;
/* TODO: for debugging; remove this */
//print("first:  ", stderr); print(*first, stderr, *debug_terminal_printer); print('\n', stderr);
//print("second: ", stderr); print(*second, stderr, *debug_terminal_printer); print('\n', stderr);
	hol_term* first_head = find_head(first, first_predicate_index, find_first_head, first_inverter);
	hol_term* second_head = find_head(second, second_predicate_index, find_second_head, second_inverter);
	if (first_head == nullptr || second_head == nullptr)
		return nullptr;

/* TODO: for debugging; remove this */
//print("first_head:  ", stderr); print(*first_head, stderr, *debug_terminal_printer); print('\n', stderr);
//print("second_head: ", stderr); print(*second_head, stderr, *debug_terminal_printer); print('\n', stderr);
	max_variable = first_inverter.max_variable;
	max_bound_variable(*first_head, max_variable);

	array<hol_term*> second_scopes(8);
	if (!get_scopes(*second_head, second_scopes))
		return nullptr;

	array_map<const hol_term*, unsigned int> second_variable_map(8);
	for (unsigned int i = first_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = first_inverter.outer[i - 1];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;

		bool has_variable = false;
		for (const hol_term* scope : second_scopes) {
			if (scope->type != hol_term_type::FOR_ALL && scope->type != hol_term_type::EXISTS && scope->type != hol_term_type::LAMBDA)
				continue;
			if (scope->quantifier.variable != node->quantifier.variable)
				continue;

			if (!has_variable) {
				++max_variable;
				has_variable = true;
			}
			if (!second_variable_map.put(scope, max_variable))
				return nullptr;
		}
	}

	for (unsigned int i = second_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = second_inverter.outer[i - 1];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;
		if (node->quantifier.variable > max_variable)
			continue;

		if (!second_variable_map.put(node, ++max_variable))
			return nullptr;
	}

	unsigned int first_head_length = 0, second_head_length = 0, first_head_min_length = 0, second_head_min_length = 0;
	if (first_inverter.outer.length > 1) {
		hol_term* first_parent = first_inverter.outer[first_inverter.outer.length - 2];
		if (first_parent->type == hol_term_type::AND || first_parent->type == hol_term_type::OR) {
			hol_term* head; head_index predicate_index;
			find_first_head(first_parent->array.operands[0], head, predicate_index);
			if (head != nullptr) {
				first_head_is_array = true;
				first_head_length = first_parent->array.length;
			}
		} else if (first_parent->type == hol_term_type::ANY_ARRAY && (first_parent->any_array.oper == hol_term_type::ANY_ARRAY || first_parent->any_array.oper == hol_term_type::AND || first_parent->any_array.oper == hol_term_type::OR)) {
			if (first_head == first_parent->any_array.all) {
				first_head_is_array = true;
				first_head_min_length = max(first_parent->any_array.any.length, max(first_parent->any_array.left.length, first_parent->any_array.right.length));
			} else {
				hol_term* head; head_index predicate_index;
				find_first_head(first_parent->any_array.all, head, predicate_index);
				if (head != nullptr) {
					first_head_is_array = true;
					first_head_min_length = max(first_parent->any_array.any.length, max(first_parent->any_array.left.length, first_parent->any_array.right.length));
				}
			}
		}
	} if (second_inverter.outer.length > 1) {
		hol_term* second_parent = second_inverter.outer[second_inverter.outer.length - 2];
		if (second_parent->type == hol_term_type::AND || second_parent->type == hol_term_type::OR) {
			hol_term* head; head_index predicate_index;
			find_second_head(second_parent->array.operands[0], head, predicate_index);
			if (head != nullptr) {
				second_head_is_array = true;
				second_head_length = second_parent->array.length;
			}
		} else if (second_parent->type == hol_term_type::ANY_ARRAY && (second_parent->any_array.oper == hol_term_type::ANY_ARRAY || second_parent->any_array.oper == hol_term_type::AND || second_parent->any_array.oper == hol_term_type::OR)) {
			if (second_head == second_parent->any_array.all) {
				second_head_is_array = true;
				second_head_min_length = max(second_parent->any_array.any.length, max(second_parent->any_array.left.length, second_parent->any_array.right.length));
			} else {
				hol_term* head; head_index predicate_index;
				find_second_head(second_parent->any_array.all, head, predicate_index);
				if (head != nullptr) {
					second_head_is_array = true;
					second_head_min_length = max(second_parent->any_array.any.length, max(second_parent->any_array.left.length, second_parent->any_array.right.length));
				}
			}
		}
	}
	if ((first_head_length != 0 && second_head_length != 0 && first_head_length != second_head_length)
	 || (first_head_length != 0 && first_head_length < second_head_min_length)
	 || (second_head_length != 0 && second_head_length < first_head_min_length))
		return nullptr;

	if (!on_remap_variables(first_head, second_head, first_inverter, second_inverter, second_variable_map, max_variable, first_head_is_array, second_head_is_array))
		return nullptr;

	return map_variables(second, second_variable_map);
}

template<typename InvertSecondFunction,
	typename FindFirstHeadFunction,
	typename FindSecondHeadFunction,
	typename OnRemapFunction>
inline bool invert_apply_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second,
		FindFirstHeadFunction find_first_head,
		FindSecondHeadFunction find_second_head,
		OnRemapFunction on_remap_variables,
		InvertSecondFunction invert_second_head)
{
	unsigned int max_variable;
	bool first_head_is_array = false, second_head_is_array = false;
	hol_term* remapped_second = remap_to_invert_apply_head(first, second, find_first_head, find_second_head, on_remap_variables, max_variable, first_head_is_array, second_head_is_array);
	if (remapped_second == nullptr)
		return false;

	apply_head_inverter first_inverter; head_index first_predicate_index;
	apply_head_inverter second_inverter; head_index second_predicate_index;
	hol_term* first_head = find_head(first, first_predicate_index, find_first_head, first_inverter);
	hol_term* second_head = find_head(remapped_second, second_predicate_index, find_second_head, second_inverter);
	if (first_head == nullptr || second_head == nullptr) {
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		return false;
	}

	hol_term* conjunct = nullptr;
	array<hol_term*> second_heads(8);
	array<array<hol_term*>> second_outer(8);
	bool any_right_only = true;
	bool could_have_wide_scope = false;
	if (second_inverter.outer.length > 1) {
		hol_term* parent = second_inverter.outer[second_inverter.outer.length - 2];
		second_inverter.outer.length--;
		/* check if the conjunction/disjunction is part of the head */
		if (parent->type == hol_term_type::AND || parent->type == hol_term_type::OR) {
			array<hol_term*>* new_heads = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->array.length);
			array<hol_term*>* new_outer = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->array.length);
			if (new_heads == nullptr || new_outer == nullptr) {
				fprintf(stderr, "invert_apply_head ERROR: Out of memory.\n");
				if (new_heads != nullptr) free(new_heads);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
				return false;
			}
			for (unsigned int i = 0; i < parent->array.length; i++) {
				if (!array_init(new_heads[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < i; j++) free(new_heads[j]);
					free(new_heads); free(new_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->array.length; i++) {
				if (!array_init(new_outer[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->array.length; j++) free(new_heads[j]);
					for (unsigned int j = 0; j < i; j++) free(new_outer[j]);
					free(new_heads); free(new_outer);
					return false;
				}
			}
			/* make sure this is an array of heads, rather than just an array */
			hol_term* temp_second_head; head_index temp_predicate_index;
			find_second_head(parent->array.operands[0], temp_second_head, temp_predicate_index);
			if (temp_second_head == nullptr) {
				for (unsigned int j = 0; j < parent->array.length; j++) { free(new_heads[j]); free(new_outer[j]); }
				free(new_heads); free(new_outer); new_heads = nullptr;
			}
			for (unsigned int i = 0; temp_second_head != nullptr && i < parent->array.length; i++) {
				hol_term* current_first_head = first_head;
				if (first_head_is_array) {
					hol_term* first_parent = first_inverter.outer[first_inverter.outer.length - 2];
					if (first_parent->type == hol_term_type::AND || first_parent->type == hol_term_type::OR) {
						current_first_head = first_parent->array.operands[i];
					} else if (first_parent->type == hol_term_type::ANY_ARRAY && (first_parent->any_array.oper == hol_term_type::ANY_ARRAY || first_parent->any_array.oper == hol_term_type::AND || first_parent->any_array.oper == hol_term_type::OR)) {
						if (i < first_parent->any_array.left.length) {
							current_first_head = first_parent->any_array.left.operands[i];
						} else if (parent->array.length - i - 1 < first_parent->any_array.right.length) {
							current_first_head = first_parent->any_array.right.operands[first_parent->any_array.right.length - parent->array.length + i];
						} else {
							current_first_head = first_parent->any_array.all;
						}
					}
				}
				if (!invert_second_head(new_heads[i], new_outer[i], current_first_head, parent->array.operands[i], first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, true)) {
					for (unsigned int j = 0; j < i; j++) { free_all(new_heads[j]); free_all(new_outer[j]); }
					for (unsigned int j = 0; j < parent->array.length; j++) { free(new_heads[j]); free(new_outer[j]); }
					free(new_heads); free(new_outer);
					new_heads = nullptr; break;
				}
			}

			if (new_heads != nullptr) {
				bool success = apply_to_cartesian_product(new_heads, parent->array.length, [parent,new_heads,new_outer,&second_heads,&second_outer](const unsigned int* index_array) {
					if (!second_heads.ensure_capacity(second_heads.length + 1)
					 || !second_outer.ensure_capacity(second_outer.length + 1)
					 || !array_init(second_outer[second_outer.length], 4))
						return false;
					second_outer.length++;

					array<hol_term*>& outer_intersections = second_outer.last();
					outer_intersections[0] = &HOL_ANY;
					outer_intersections.length = 1;
					HOL_ANY.reference_count++;
					for (unsigned int i = 0; i < parent->array.length; i++) {
						array<hol_term*> new_outer_terms(4);
						for (hol_term* outer : outer_intersections)
							intersect<built_in_predicates>(new_outer_terms, outer, new_outer[i][index_array[i]]);
						free_all(outer_intersections);
						swap(new_outer_terms, outer_intersections);
					}

					hol_term* new_head;
					if (parent->type == hol_term_type::AND)
						new_head = hol_term::new_and(lookup_table_array_view<array<hol_term*>, hol_term*>(new_heads, index_array, parent->array.length));
					else new_head = hol_term::new_or(lookup_table_array_view<array<hol_term*>, hol_term*>(new_heads, index_array, parent->array.length));
					if (new_head == nullptr)
						return false;
					for (unsigned int i = 0; i < parent->array.length; i++)
						new_head->array.operands[i]->reference_count++;
					second_heads[second_heads.length++] = new_head;
					return true;
				});
				for (unsigned int j = 0; j < parent->array.length; j++) { free_all(new_heads[j]); free(new_heads[j]); free_all(new_outer[j]); free(new_outer[j]); }
				free(new_heads); free(new_outer);
				if (!success) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
					free_all(second_heads); return false;
				}
				second_head = parent;
			}
		} else if (parent->type == hol_term_type::ANY_ARRAY && (parent->any_array.oper == hol_term_type::ANY_ARRAY || parent->any_array.oper == hol_term_type::AND || parent->any_array.oper == hol_term_type::OR)) {
			array<hol_term*> new_all(8); array<hol_term*> new_all_outer(8);
			array<hol_term*>* new_left_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.left.length);
			array<hol_term*>* new_left_outer = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.left.length);
			array<hol_term*>* new_right_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.right.length);
			array<hol_term*>* new_right_outer = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.right.length);
			array<hol_term*>* new_any_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.any.length);
			array<hol_term*>* new_any_outer = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->any_array.any.length);
			if (new_left_array == nullptr || new_left_outer == nullptr || new_right_array == nullptr || new_right_outer == nullptr || new_any_array == nullptr || new_any_outer == nullptr) {
				fprintf(stderr, "invert_apply_head ERROR: Out of memory.\n");
				if (new_left_array != nullptr) free(new_left_array);
				if (new_left_outer != nullptr) free(new_left_outer);
				if (new_right_array != nullptr) free(new_right_array);
				if (new_right_outer != nullptr) free(new_right_outer);
				if (new_any_array != nullptr) free(new_any_array);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
				return false;
			}
			for (unsigned int i = 0; i < parent->any_array.left.length; i++) {
				if (!array_init(new_left_array[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < i; j++) free(new_left_array[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->any_array.left.length; i++) {
				if (!array_init(new_left_outer[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) free(new_left_array[j]);
					for (unsigned int j = 0; j < i; j++) free(new_left_outer[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->any_array.right.length; i++) {
				if (!array_init(new_right_array[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free(new_left_array[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < i; j++) free(new_right_array[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->any_array.right.length; i++) {
				if (!array_init(new_right_outer[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free(new_left_array[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) free(new_right_array[j]);
					for (unsigned int j = 0; j < i; j++) free(new_right_outer[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->any_array.any.length; i++) {
				if (!array_init(new_any_array[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free(new_left_array[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free(new_right_array[j]); free(new_right_outer[j]); }
					for (unsigned int j = 0; j < i; j++) free(new_any_array[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			} for (unsigned int i = 0; i < parent->any_array.any.length; i++) {
				if (!array_init(new_any_outer[i], 8)) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free(new_left_array[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free(new_right_array[j]); free(new_right_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.any.length; j++) free(new_any_array[j]);
					for (unsigned int j = 0; j < i; j++) free(new_any_outer[j]);
					free(new_left_array); free(new_right_array); free(new_any_array);
					free(new_left_outer); free(new_right_outer); free(new_any_outer);
					return false;
				}
			}

			if (!invert_second_head(new_all, new_all_outer, first_head, parent->any_array.all, first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, true)) {
				for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free_all(new_left_array[j]); free(new_left_array[j]); free_all(new_left_outer[j]); free(new_left_outer[j]); }
				for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free_all(new_right_array[j]); free(new_right_array[j]); free_all(new_right_outer[j]); free(new_right_outer[j]); }
				for (unsigned int j = 0; j < parent->any_array.any.length; j++) { free_all(new_any_array[j]); free(new_any_array[j]); free_all(new_any_outer[j]); free(new_any_outer[j]); }
				free(new_left_array); free(new_right_array); free(new_any_array); free(new_left_outer); free(new_right_outer); free(new_any_outer);
				new_left_array = nullptr;
			}
			for (unsigned int i = 0; new_left_array != nullptr && i < parent->any_array.left.length; i++) {
				hol_term* current_first_head = first_head;
				if (first_head_is_array) {
					hol_term* first_parent = first_inverter.outer[first_inverter.outer.length - 2];
					if (first_parent->type == hol_term_type::AND || first_parent->type == hol_term_type::OR) {
						current_first_head = first_parent->array.operands[i];
					} else if (first_parent->type == hol_term_type::ANY_ARRAY && (first_parent->any_array.oper == hol_term_type::ANY_ARRAY || first_parent->any_array.oper == hol_term_type::AND || first_parent->any_array.oper == hol_term_type::OR)) {
						if (i < first_parent->any_array.left.length) {
							current_first_head = first_parent->any_array.left.operands[i];
						} else {
							current_first_head = first_parent->any_array.all;
						}
					}
				}
				if (!invert_second_head(new_left_array[i], new_left_outer[i], current_first_head, parent->any_array.left.operands[i], first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, true)) {
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free_all(new_left_array[j]); free(new_left_array[j]); free_all(new_left_outer[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free_all(new_right_array[j]); free(new_right_array[j]); free_all(new_right_outer[j]); free(new_right_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.any.length; j++) { free_all(new_any_array[j]); free(new_any_array[j]); free_all(new_any_outer[j]); free(new_any_outer[j]); }
					free(new_left_array); free(new_right_array); free(new_any_array); free(new_left_outer); free(new_right_outer); free(new_any_outer);
					free_all(new_all); free_all(new_all_outer); new_left_array = nullptr;
					break;
				}
			} for (unsigned int i = 0; new_left_array != nullptr && i < parent->any_array.right.length; i++) {
				hol_term* current_first_head = first_head;
				if (first_head_is_array) {
					hol_term* first_parent = first_inverter.outer[first_inverter.outer.length - 2];
					if (first_parent->type == hol_term_type::AND || first_parent->type == hol_term_type::OR) {
						current_first_head = first_parent->array.operands[first_parent->array.length - parent->any_array.right.length + i];
					} else if (first_parent->type == hol_term_type::ANY_ARRAY && (first_parent->any_array.oper == hol_term_type::ANY_ARRAY || first_parent->any_array.oper == hol_term_type::AND || first_parent->any_array.oper == hol_term_type::OR)) {
						if (first_parent->any_array.right.length - i - 1 < parent->any_array.right.length) {
							current_first_head = first_parent->any_array.right.operands[parent->any_array.right.length - first_parent->any_array.right.length + i];
						} else {
							current_first_head = first_parent->any_array.all;
						}
					}
				}
				if (!invert_second_head(new_right_array[i], new_right_outer[i], current_first_head, parent->any_array.right.operands[i], first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, true)) {
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free_all(new_left_array[j]); free(new_left_array[j]); free_all(new_left_outer[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free_all(new_right_array[j]); free(new_right_array[j]); free_all(new_right_outer[j]); free(new_right_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.any.length; j++) { free_all(new_any_array[j]); free(new_any_array[j]); free_all(new_any_outer[j]); free(new_any_outer[j]); }
					free(new_left_array); free(new_right_array); free(new_any_array); free(new_left_outer); free(new_right_outer); free(new_any_outer);
					free_all(new_all); free_all(new_all_outer); new_left_array = nullptr;
					break;
				}
			} for (unsigned int i = 0; new_left_array != nullptr && i < parent->any_array.any.length; i++) {
				hol_term* current_first_head = first_head;
				if (first_head_is_array) {
					hol_term* first_parent = first_inverter.outer[first_inverter.outer.length - 2];
					if (first_parent->type == hol_term_type::ANY_ARRAY && (first_parent->any_array.oper == hol_term_type::ANY_ARRAY || first_parent->any_array.oper == hol_term_type::AND || first_parent->any_array.oper == hol_term_type::OR)) {
						current_first_head = first_parent->any_array.all;
					}
				}
				if (!invert_second_head(new_any_array[i], new_any_outer[i], current_first_head, parent->any_array.any.operands[i], first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, true)) {
					for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free_all(new_left_array[j]); free(new_left_array[j]); free_all(new_left_outer[j]); free(new_left_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free_all(new_right_array[j]); free(new_right_array[j]); free_all(new_right_outer[j]); free(new_right_outer[j]); }
					for (unsigned int j = 0; j < parent->any_array.any.length; j++) { free_all(new_any_array[j]); free(new_any_array[j]); free_all(new_any_outer[j]); free(new_any_outer[j]); }
					free(new_left_array); free(new_right_array); free(new_any_array); free(new_left_outer); free(new_right_outer); free(new_any_outer);
					free_all(new_all); free_all(new_all_outer); new_left_array = nullptr;
					break;
				}
			}
			if (new_left_array != nullptr) {
				bool success = true;
				for (unsigned int i = 0; i < new_all.length; i++) {
					hol_term* all = new_all[i];
					hol_term* all_outer = new_all_outer[i];
					success = apply_to_cartesian_product(new_left_array, parent->any_array.left.length, [all,all_outer,new_left_array,new_left_outer,new_right_array,new_right_outer,new_any_array,new_any_outer,&second_heads,&second_outer,parent](const unsigned int* left_index_array) {
						return apply_to_cartesian_product(new_right_array, parent->any_array.right.length, [all,all_outer,new_left_array,new_left_outer,new_right_array,new_right_outer,new_any_array,new_any_outer,&second_heads,&second_outer,parent,left_index_array](const unsigned int* right_index_array) {
							return apply_to_cartesian_product(new_any_array, parent->any_array.any.length, [all,all_outer,new_left_array,new_left_outer,new_right_array,new_right_outer,new_any_array,new_any_outer,&second_heads,&second_outer,parent,left_index_array,right_index_array](const unsigned int* any_index_array) {
								if (!second_heads.ensure_capacity(second_heads.length + 1)
								 || !second_outer.ensure_capacity(second_outer.length + 1)
								 || !array_init(second_outer[second_outer.length], 4))
									return false;
								second_outer.length++;

								array<hol_term*>& outer_terms = second_outer.last();
								outer_terms[0] = all_outer;
								outer_terms.length = 1;
								all_outer->reference_count++;
								for (unsigned int i = 0; i < parent->any_array.left.length; i++) {
									array<hol_term*> new_outer_terms(4);
									for (hol_term* outer : outer_terms)
										intersect<built_in_predicates>(new_outer_terms, outer, new_left_outer[i][left_index_array[i]]);
									free_all(outer_terms);
									swap(new_outer_terms, outer_terms);
								} for (unsigned int i = 0; i < parent->any_array.right.length; i++) {
									array<hol_term*> new_outer_terms(4);
									for (hol_term* outer : outer_terms)
										intersect<built_in_predicates>(new_outer_terms, outer, new_right_outer[i][right_index_array[i]]);
									free_all(outer_terms);
									swap(new_outer_terms, outer_terms);
								} for (unsigned int i = 0; i < parent->any_array.any.length; i++) {
									array<hol_term*> new_outer_terms(4);
									for (hol_term* outer : outer_terms)
										intersect<built_in_predicates>(new_outer_terms, outer, new_any_outer[i][any_index_array[i]]);
									free_all(outer_terms);
									swap(new_outer_terms, outer_terms);
								}

								hol_term* new_head = hol_term::new_any_array(parent->any_array.oper, all,
										lookup_table_array_view<array<hol_term*>, hol_term*>(new_any_array, any_index_array, parent->any_array.any.length),
										lookup_table_array_view<array<hol_term*>, hol_term*>(new_left_array, left_index_array, parent->any_array.left.length),
										lookup_table_array_view<array<hol_term*>, hol_term*>(new_right_array, right_index_array, parent->any_array.right.length));
								if (new_head == nullptr)
									return false;
								new_head->any_array.all->reference_count++;
								for (unsigned int i = 0; i < new_head->any_array.left.length; i++)
									new_head->any_array.left.operands[i]->reference_count++;
								for (unsigned int i = 0; i < new_head->any_array.right.length; i++)
									new_head->any_array.right.operands[i]->reference_count++;
								for (unsigned int i = 0; i < new_head->any_array.any.length; i++)
									new_head->any_array.any.operands[i]->reference_count++;
								second_heads[second_heads.length++] = new_head;
								return true;
							});
						});
					});
					if (!success) break;
				}
				for (unsigned int j = 0; j < parent->any_array.left.length; j++) { free_all(new_left_array[j]); free(new_left_array[j]); free_all(new_left_outer[j]); free(new_left_outer[j]); }
				for (unsigned int j = 0; j < parent->any_array.right.length; j++) { free_all(new_right_array[j]); free(new_right_array[j]); free_all(new_right_outer[j]); free(new_right_outer[j]); }
				for (unsigned int j = 0; j < parent->any_array.any.length; j++) { free_all(new_any_array[j]); free(new_any_array[j]); free_all(new_any_outer[j]); free(new_any_outer[j]); }
				free(new_left_array); free(new_right_array); free(new_any_array); free(new_left_outer); free(new_right_outer); free(new_any_outer); free_all(new_all); free_all(new_all_outer);
				if (!success) {
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
					free_all(second_heads); return false;
				}
				second_head = parent;
			}
		}
	}

	if (second_heads.length == 0) {
		array<hol_term*> outer(8);
		if (!invert_second_head(second_heads, outer, first_head, second_head, first_inverter, second_inverter, first_predicate_index, second_predicate_index, conjunct, max_variable, any_right_only, could_have_wide_scope, false)
		 || !second_outer.ensure_capacity(outer.length))
		{
			free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
			if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
			free_all(second_heads); free_all(outer); return false;
		}
		for (unsigned int i = 0; i < outer.length; i++) {
			if (!array_init(second_outer[i], 1)) {
				for (unsigned int j = 0; j < i; j++) free(second_outer[j]);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
				free_all(second_heads); free_all(outer); return false;
			}
			second_outer[i][0] = outer[i];
			second_outer[i].length = 1;
			second_outer.length++;
		}
	}
	if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }

	if (first_head_is_array) {
		hol_term* parent = first_inverter.outer[first_inverter.outer.length - 2];
		first_head = parent;
	}

	hol_term* first_outer = substitute_head<any_node_position::RIGHT>(first, first_head, &HOL_ZERO, true);
	if (first_outer == nullptr) {
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		return false;
	}

	array<hol_term*> inverted_logical_forms(8);
	for (unsigned int i = 0; i < second_heads.length; i++) {
		array<hol_term*> new_heads(8);
/* TODO: for debugging; remove this */
//fprintf(stderr, "first_head:       "); print(*first_head, stderr, *debug_terminal_printer); print('\n', stderr);
//fprintf(stderr, "second_heads[%u]: ", i); print(*second_heads[i], stderr, *debug_terminal_printer); print('\n', stderr);
		intersect<built_in_predicates>(new_heads, first_head, second_heads[i]);
		if (new_heads.length == 0)
			continue;

		for (unsigned int j = 0; j < second_outer[i].length; j++) {
			hol_term* new_second_outer;
			if (any_right_only)
				new_second_outer = substitute_head<any_node_position::LEFT, true>(remapped_second, second_head, second_outer[i][j], could_have_wide_scope);
			else new_second_outer = substitute_head<any_node_position::LEFT, false>(remapped_second, second_head, second_outer[i][j], could_have_wide_scope);
			if (new_second_outer == nullptr) {
				free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
				free_all(second_heads); free_all(new_heads); free_all(inverted_logical_forms); return false;
			}

			array<hol_term*> new_outer(8);
/* TODO: for debugging; remove this */
//fprintf(stderr, "first_outer:      "); print(*first_outer, stderr, *debug_terminal_printer); print('\n', stderr);
//fprintf(stderr, "new_second_outer: "); print(*new_second_outer, stderr, *debug_terminal_printer); print('\n', stderr);
			intersect<built_in_predicates>(new_outer, first_outer, new_second_outer);
			free(*new_second_outer); if (new_second_outer->reference_count == 0) free(new_second_outer);

			unsigned int index = 0;
			for (unsigned int k = 0; k < new_outer.length; k++) {
				hol_term* cleaned_outer = remove_any_nodes(new_outer[k], find_zero);
				if (cleaned_outer == nullptr) {
					free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
					free_all(second_heads); free_all(new_heads); free_all(inverted_logical_forms);
					for (unsigned int l = 0; l < index; l++) { free(*new_outer[l]); if (new_outer[l]->reference_count == 0) free(new_outer[l]); }
					for (unsigned int l = k; l < new_outer.length; l++) { free(*new_outer[l]); if (new_outer[l]->reference_count == 0) free(new_outer[l]); }
					return false;
				}
				free(*new_outer[k]); if (new_outer[k]->reference_count == 0) free(new_outer[k]);

				/* check if we already have `cleaned_outer` */
				bool is_duplicate = false;
				for (unsigned int l = 0; l < index; l++) {
					if (*cleaned_outer == *new_outer[l]) {
						is_duplicate = true;
						break;
					}
				}
				if (is_duplicate) {
					free(*cleaned_outer); if (cleaned_outer->reference_count == 0) free(cleaned_outer);
					new_outer[k] = new_outer[new_outer.length - 1];
					k--; new_outer.length--;
				} else {
					new_outer[index++] = cleaned_outer;
				}
			}
			new_outer.length = index;

			if (!inverted_logical_forms.ensure_capacity(inverted_logical_forms.length + new_outer.length * new_heads.length)) {
				free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
				free_all(second_heads); free_all(new_heads); free_all(new_outer); free_all(inverted_logical_forms);
				return false;
			}
			for (hol_term* outer : new_outer) {
				for (hol_term* head : new_heads) {
					inverted_logical_forms[inverted_logical_forms.length] = substitute_head<any_node_position::NONE>(outer, &HOL_ZERO, head);
					if (inverted_logical_forms[inverted_logical_forms.length] == nullptr) {
						free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
						free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
						for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
						free_all(second_heads); free_all(new_heads); free_all(new_outer); free_all(inverted_logical_forms);
						return false;
					}
					inverted_logical_forms.length++;
				}
			}
			free_all(new_outer);
		}
		free_all(new_heads);
	}
	free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
	free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
	for (array<hol_term*>& outer : second_outer) { free_all(outer); free(outer); }
	free_all(second_heads);

	if (inverted_logical_forms.length == 0)
		return false;

#if !defined(NDEBUG)
	bool is_disjoint = true;
	for (unsigned int i = 1; is_disjoint && i < inverted_logical_forms.length; i++) {
		for (unsigned int j = 0; j < i; j++) {
			if (has_intersection<built_in_predicates>(inverted_logical_forms[i], inverted_logical_forms[j])) {
				fprintf(stderr, "invert_apply_head WARNING: `inverted_logical_forms` is not pairwise-disjoint.\n");
/* TODO: for debugging; delete this */
has_intersection<built_in_predicates>(inverted_logical_forms[i], inverted_logical_forms[j]);
				is_disjoint = false; break;
			}
		}
	}
#endif

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * inverted_logical_forms.length);
	if (inverse == nullptr) {
		free_all(inverted_logical_forms);
		return false;
	}
	for (unsigned int i = 0; i < inverted_logical_forms.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = inverted_logical_forms[i];
/* TODO: for debugging; remove this */
//fprintf(stderr, "inverse[%u]: ", i); print(inverse[i], stderr, *debug_terminal_printer); print('\n', stderr);
	}
	inverse_count = inverted_logical_forms.length;
	return true;
}

inline bool get_target_variable(
		unsigned int src_variable, unsigned int& target_var,
		const variable_set& target_var_set, unsigned int& max_variable)
{
	target_var = 0;
	switch (target_var_set.type) {
	case variable_set_type::SINGLETON:
		target_var = target_var_set.variable;
		break;
	case variable_set_type::ANY:
		if (index_of(src_variable, target_var_set.any.array, target_var_set.any.length) < target_var_set.any.length) {
			target_var = src_variable;
			break;
		}
		for (unsigned int i = 0; i < target_var_set.any.length; i++) {
			if (target_var_set.any.array[i] >= max_variable) {
				target_var = target_var_set.any.array[i];
				max_variable = target_var;
				break;
			}
		}
		break;
	case variable_set_type::ANY_EXCEPT:
		if (index_of(src_variable, target_var_set.any.array, target_var_set.any.length) == target_var_set.any.length) {
			target_var = src_variable;
			break;
		}
		target_var = max_variable;
		for (unsigned int i = 0; i < target_var_set.any.length; i++)
			target_var = max(target_var, target_var_set.any.array[i]);
		target_var++;
		max_variable = target_var;
		break;
	}
	return (target_var != 0);
}

bool remap_scopes(array_map<unsigned int, variable_set>& free_variable_map,
		apply_head_inverter& first_head_inverter, apply_head_inverter& second_head_inverter,
		array_map<const hol_term*, unsigned int>& second_variable_map, unsigned int& max_variable,
		bool first_head_is_array, bool second_head_is_array)
{
	unsigned int prev_first_inverter_index = first_head_inverter.outer.length;
	unsigned int prev_second_inverter_index = second_head_inverter.outer.length;
	while (free_variable_map.size != 0) {
		unsigned int index = 0;

		/* find the next inner-most scope in `second` whose variable is in `free_variable_map.keys` */
		unsigned int second_inverter_index;
		for (second_inverter_index = prev_second_inverter_index; second_inverter_index > 0; second_inverter_index--) {
			const hol_term* node = second_head_inverter.outer[second_inverter_index - 1];
			if (node->type == hol_term_type::FOR_ALL || node->type == hol_term_type::EXISTS || node->type == hol_term_type::LAMBDA) {
				index = free_variable_map.index_of(node->quantifier.variable);
				if (index < free_variable_map.size) break;
				prev_second_inverter_index = second_inverter_index;
			}
		}

		unsigned int second_variable = free_variable_map.keys[index];
		if (free_variable_map.values[index].type != variable_set_type::SINGLETON) {
			free(free_variable_map.values[index]);
			free_variable_map.remove_at(index);
			prev_second_inverter_index = second_inverter_index;
			continue;
		}
		unsigned int first_variable = free_variable_map.values[index].variable;
		free_variable_map.remove_at(index);

		index = second_variable_map.index_of(second_head_inverter.outer[second_inverter_index - 1]);
		if (index < second_variable_map.size) {
			second_variable_map.values[index] = first_variable;
			if (first_variable == second_variable)
				second_variable_map.remove_at(index);
		}

		/* find the scopes where `first_variable` and `second_variable` are declared */
		unsigned int first_inverter_index;
		for (first_inverter_index = prev_first_inverter_index; first_inverter_index > 0; first_inverter_index--) {
			const hol_term* node = first_head_inverter.outer[first_inverter_index - 1];
			if (node->type == hol_term_type::FOR_ALL || node->type == hol_term_type::EXISTS || node->type == hol_term_type::LAMBDA) {
				if (node->quantifier.variable == first_variable) break;
				prev_first_inverter_index = first_inverter_index;
			}
		}

		/* make sure the two scopes can intersect */
		hol_term* first_scope = substitute_head<any_node_position::NONE>(
				first_head_inverter.outer[first_inverter_index - 1], first_head_inverter.outer[prev_first_inverter_index - 1], &HOL_ZERO);
		hol_term* second_scope = substitute_head<any_node_position::NONE>(
				second_head_inverter.outer[second_inverter_index - 1], second_head_inverter.outer[prev_second_inverter_index - 1], &HOL_ZERO);
		if (first_scope == nullptr || second_scope == nullptr) {
			if (first_scope == nullptr) { free(*first_scope); free(first_scope); }
			return false;
		}

		array<pair<hol_term*, variable_map>> intersection(2);
		intersect<built_in_predicates, true, true>(intersection, first_scope, second_scope);
		if (intersection.length == 0) {
			free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
			free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
			return false;
		} else if (intersection.length != 1) {
			fprintf(stderr, "remap_scopes ERROR: Intersection is not unique.\n");
			free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
			free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
			free_all(intersection); return false;
		}

		const variable_map& var_map = intersection[0].value;
		for (const auto& entry : var_map.scope_map) {
			const variable_set& set = entry.value;
			unsigned int src_variable = entry.key.src->quantifier.variable;
			unsigned int target_var;
			if (!get_target_variable(src_variable, target_var, set, max_variable)
			 || !second_variable_map.ensure_capacity(second_variable_map.size + 1))
			{
				free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
				free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
				free_all(intersection); return false;
			}

			unsigned int index = second_variable_map.index_of(entry.key.src);
			if (index < second_variable_map.size) {
				second_variable_map.values[index] = target_var;
				if (src_variable == target_var)
					second_variable_map.remove_at(index);
			} else if (src_variable != target_var) {
				second_variable_map.keys[second_variable_map.size] = entry.key.src;
				second_variable_map.values[second_variable_map.size] = target_var;
				second_variable_map.size++;
			}
		} for (const auto& entry : var_map.free_variables) {
			const variable_set& set = entry.value;
			unsigned int src_variable = entry.key;
			if (!free_variable_map.ensure_capacity(free_variable_map.size + 1)) {
				free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
				free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
				free_all(intersection); return false;
			}
			unsigned int index = free_variable_map.index_of(src_variable);
			if (index < free_variable_map.size) {
				variable_set& new_set = *((variable_set*) alloca(sizeof(variable_map)));
				if (!intersect(new_set, set, free_variable_map.values[index])) {
					free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
					free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
					free_all(intersection); return false;
				}
				swap(free_variable_map.values[index], new_set);
				free(new_set);
			} else {
				free_variable_map.keys[free_variable_map.size] = src_variable;
				if (!init(free_variable_map.values[free_variable_map.size], set)) {
					free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
					free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
					free_all(intersection); return false;
				}
				free_variable_map.size++;
			}
		}
		free(*first_scope); if (first_scope->reference_count == 0) free(first_scope);
		free(*second_scope); if (second_scope->reference_count == 0) free(second_scope);
		free_all(intersection);

		prev_first_inverter_index = first_inverter_index;
		prev_second_inverter_index = second_inverter_index;

		if (first_head_is_array && prev_first_inverter_index == first_head_inverter.outer.length)
			prev_first_inverter_index--;
		if (second_head_is_array && prev_second_inverter_index == second_head_inverter.outer.length)
			prev_second_inverter_index--;
	}
	return true;
}

bool intersect(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		grammatical_flags flags,
		hol_term* first, hol_term* second)
{
	array<hol_term*> intersection(8);
	if (!intersect<built_in_predicates>(intersection, first, second))
		return false;
	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * intersection.length);
	if (inverse == nullptr) {
		fprintf(stderr, "intersect ERROR: Insufficient memory for `inverse`.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	for (unsigned int i = 0; i < intersection.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = intersection[i];
	}
	inverse_count = intersection.length;
	return true;
}

bool intersect_with_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		grammatical_flags flags,
		hol_term* first, hol_term* second)
{
	head_index first_predicate_index, second_predicate_index; no_op apply;
	auto find_array_head = make_array_finder(find_head<built_in_predicates>);
	hol_term* first_head = find_head(first, first_predicate_index, find_array_head, apply);
	hol_term* second_head = find_head(second, second_predicate_index, find_array_head, apply);
	if (first_head == nullptr || second_head == nullptr)
		return false;

	hol_term* first_outer = substitute_head<any_node_position::NONE>(first, first_head, &HOL_ZERO);
	if (first_outer == nullptr)
		return false;
	HOL_ZERO.reference_count++;

	hol_term* second_outer = substitute_head<any_node_position::NONE>(second, second_head, &HOL_ZERO);
	if (second_outer == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		return false;
	}
	HOL_ZERO.reference_count++;

	/* intersect the heads and the outer terms respectively */
	array<hol_term*> new_heads(4);
	intersect<built_in_predicates>(new_heads, first_head, second_head);
	if (new_heads.length == 0) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
		return false;
	}

	array<hol_term*> new_outer(4);
	intersect<built_in_predicates>(new_outer, first_outer, second_outer);
	free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
	free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
	if (new_outer.length == 0) {
		free_all(new_heads);
		return false;
	}

	array<hol_term*> intersection(new_heads.length * new_outer.length);
	for (hol_term* outer : new_outer) {
		for (hol_term* head : new_heads) {
			hol_term* new_term = substitute_head<any_node_position::NONE>(outer, &HOL_ZERO, head);
			if (new_term == nullptr) {
				free_all(new_heads); free_all(new_outer); free_all(intersection);
				return false;
			}
			intersection[intersection.length++] = new_term;
		}
	}
	free_all(new_heads);
	free_all(new_outer);

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * intersection.length);
	if (inverse == nullptr) {
		fprintf(stderr, "intersect_with_head ERROR: Out of memory.\n");
		free_all(intersection); return false;
	}
	for (unsigned int i = 0; i < intersection.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = intersection[i];
	}
	inverse_count = intersection.length;
	return true;
}

template<int_fast8_t ConjunctIndex, bool SelectNegation = false>
inline bool invert_select_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	auto on_remap_variables = [](hol_term* first_head, hol_term* second_head,
			apply_head_inverter& first_head_inverter, apply_head_inverter& second_head_inverter,
			array_map<const hol_term*, unsigned int>& second_variable_map, unsigned int max_variable,
			bool first_head_is_array, bool second_head_is_array)
	{
		if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr && first_head->any.included->type == hol_term_type::EXISTS)
			first_head = first_head->any.included;
		if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr && second_head->any.included->type == hol_term_type::EXISTS)
			second_head = second_head->any.included;

		while (first_head->type == hol_term_type::NOT)
			first_head = first_head->unary.operand;
		while (second_head->type == hol_term_type::NOT)
			second_head = second_head->unary.operand;

		hol_term* expected_head;
		if (ConjunctIndex >= 0) {
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(&HOL_ANY, ConjunctIndex + 1), make_array_view((hol_term**) nullptr, 0)));
			if (expected_head == nullptr)
				return false;
			HOL_ANY.reference_count += 2 + ConjunctIndex;
		} else {
			unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(&HOL_ANY, index + 1)));
			if (expected_head == nullptr)
				return false;
			HOL_ANY.reference_count += 2 + index;
		}

		array<hol_term*> intersection(2);
		intersect<built_in_predicates>(intersection, first_head, expected_head);
		free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
		if (intersection.length == 0) {
			return false;
		} else if (intersection.length != 1) {
			fprintf(stderr, "invert_select_conjunct ERROR: Intersection is not unique.\n");
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}

		hol_term* operand;
		if (intersection[0]->type == hol_term_type::ANY_QUANTIFIER)
			operand = intersection[0]->any_quantifier.operand;
		else operand = intersection[0]->quantifier.operand;

		hol_term* conjunct;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			if (ConjunctIndex >= 0) {
				conjunct = operand->any_array.left.operands[ConjunctIndex];
			} else {
				conjunct = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex];
			}
		} else {
			unsigned int index;
			if (ConjunctIndex >= 0) {
				index = ConjunctIndex;
			} else {
				index = operand->array.length + ConjunctIndex;
			}
			conjunct = operand->array.operands[index];
		}

		/* now get the conjunct from `second_head` */
		expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_and(&HOL_ANY, &HOL_ANY));
		if (expected_head == nullptr) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		HOL_ANY.reference_count += 2;

		array<hol_term*> second_intersection(2);
		intersect<built_in_predicates>(second_intersection, second_head, expected_head);
		free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
		if (second_intersection.length == 0) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		} else if (second_intersection.length != 1) {
			fprintf(stderr, "invert_select_conjunct ERROR: Intersection is not unique.\n");
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}

		if (second_intersection[0]->type == hol_term_type::ANY_QUANTIFIER)
			operand = second_intersection[0]->any_quantifier.operand;
		else operand = second_intersection[0]->quantifier.operand;

		hol_term* second_conjunct = operand->array.operands[1];
		array<pair<hol_term*, variable_map>> conjunct_intersection(2);
		intersect<built_in_predicates>(conjunct_intersection, second_conjunct, conjunct);
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
		if (conjunct_intersection.length == 0) {
			return false;
		} else if (conjunct_intersection.length != 1) {
			fprintf(stderr, "invert_select_conjunct ERROR: Intersection module variable relabeling is not unique.\n");
			free_all(conjunct_intersection); return false;
		}

		/* make sure the variables in `second_head` map to the correct variables in `first_head` */
		const variable_map& var_map = conjunct_intersection[0].value;
		for (const auto& entry : var_map.scope_map) {
			if (!second_variable_map.contains(entry.key.src))
				continue;

			unsigned int target_var = 0;
			if (!get_target_variable(entry.key.src->quantifier.variable, target_var, entry.value, max_variable)) {
				free_all(conjunct_intersection);
				return false;
			}

			if (!second_variable_map.ensure_capacity(second_variable_map.size + 1)) {
				free_all(conjunct_intersection);
				return false;
			}
			unsigned int index = second_variable_map.index_of(entry.key.src);
			if (target_var == second_variable_map.values[index]) {
				/* `second_variable_map` should map `var` to itself */
				second_variable_map.remove_at(index);
			} else {
				/* `second_variable_map` should map `var` to `target_var` */
				second_variable_map.values[index] = target_var;
				if (index == second_variable_map.size) {
					second_variable_map.keys[index] = entry.key.src;
					second_variable_map.size++;
				}
			}
		}

		array_map<unsigned int, variable_set> free_variable_map(8);
		for (const auto& entry : var_map.free_variables) {
			const variable_set& set = entry.value;
			unsigned int src_variable = entry.key;
			if (!free_variable_map.ensure_capacity(free_variable_map.size + 1)) {
				free_all(conjunct_intersection);
				return false;
			}
			unsigned int index = free_variable_map.index_of(src_variable);
			if (index < free_variable_map.size) {
				variable_set& new_set = *((variable_set*) alloca(sizeof(variable_map)));
				if (!intersect(new_set, set, free_variable_map.values[index])) {
					free_all(conjunct_intersection);
					return false;
				}
				swap(free_variable_map.values[index], new_set);
				free(new_set);
			} else {
				free_variable_map.keys[free_variable_map.size] = src_variable;
				if (!init(free_variable_map.values[free_variable_map.size], set)) {
					free_all(conjunct_intersection);
					return false;
				}
				free_variable_map.size++;
			}
		}
		free_all(conjunct_intersection);

		if (first_head_inverter.outer.last() != first_head && !first_head_inverter.outer.add(first_head)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		} if (second_head_inverter.outer.last() != second_head && !second_head_inverter.outer.add(second_head)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		}

		if (!remap_scopes(free_variable_map, first_head_inverter, second_head_inverter, second_variable_map, max_variable, first_head_is_array, second_head_is_array)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		}
		return true;
	};

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, on_remap_variables,
		[second](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			any_right_only = false;
			hol_term* old_second_head = second_head;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_predicate_index.position != head_position::NONE)
				second_head = second_head->any.included;

			unsigned int first_min_negations = 0;
			unsigned int second_negations = 0;
			while (first_head->type == hol_term_type::NOT) {
				first_head = first_head->unary.operand;
				first_min_negations++;
			} while (second_head->type == hol_term_type::NOT) {
				second_head = second_head->unary.operand;
				second_negations++;
			}

			/* if `SelectNegation` is true, `first_negations` should be equal to `second_negations`,
			   otherwise `first_negations` should be equal to either `second_negations + 1` or `0` */
			unsigned int expected_first_negations = (SelectNegation ? second_negations : second_negations + 1);
			bool could_be_zero_negations = (!SelectNegation && second_negations == 0);
			if (first_min_negations > expected_first_negations)
				return false;

			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;
				if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT || (first_head->type == hol_term_type::EXISTS
				 && (first_head->quantifier.operand->type == hol_term_type::ANY || first_head->quantifier.operand->type == hol_term_type::ANY_RIGHT || first_head->quantifier.operand->type == hol_term_type::ANY_ARRAY)))
				{
					hol_term* head_var = hol_term::new_variable(second_head->quantifier.variable);
					if (head_var == nullptr) return false;
					constexpr unsigned int excluded_tree_count = 5;
					hol_term* excluded_trees[excluded_tree_count];
					excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
					excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
					excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
					excluded_trees[3] = hol_term::new_any(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY));
					excluded_trees[4] = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
					if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
					if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[4] != nullptr) { HOL_ANY.reference_count++; }
					if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr || excluded_trees[4] == nullptr) {
						if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
						if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
						if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
						if (excluded_trees[3] != nullptr) { free(*excluded_trees[3]); free(excluded_trees[3]); }
						free(*head_var); free(head_var);
						return false;
					}
					free(*head_var);

					hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
					if (conjunct == nullptr) {
						free(*excluded_trees[0]); free(excluded_trees[0]);
						free(*excluded_trees[1]); free(excluded_trees[1]);
						free(*excluded_trees[2]); free(excluded_trees[2]);
						free(*excluded_trees[3]); free(excluded_trees[3]);
						free(*excluded_trees[4]); free(excluded_trees[4]);
						return false;
					}

					hol_term* conjunction = nullptr;
					if (ConjunctIndex >= 0) {
						conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct,
								make_array_view(&second_head_operand->array.operands[0], 1),
								make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), second_head_operand->array.operands[1]),
								make_array_view((hol_term**) nullptr, 0));
					} else if (ConjunctIndex < 0) {
						unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						conjunction = hol_term::new_any_array(hol_term_type::AND, conjunct,
								make_array_view(&second_head_operand->array.operands[0], 1),
								make_array_view((hol_term**) nullptr, 0),
								make_prepended_array_view(second_head_operand->array.operands[1], make_repeated_array_view(conjunct, index)));
					}
					if (conjunction == nullptr) {
						free(*conjunct); free(conjunct);
						return false;
					}
					conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
						conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
						conjunction->any_array.right.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
						conjunction->any_array.any.operands[i]->reference_count++;
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);

					hol_term* new_second_head = hol_term::new_exists(second_head->quantifier.variable, conjunction);
					if (new_second_head == nullptr) {
						free(*conjunction); free(conjunction);
						return false;
					}

					unsigned int dst_start = 0, dst_end = 0;
					if (could_be_zero_negations) {
						dst_start = dst.length;
						intersect<built_in_predicates>(dst, new_second_head, first_head);
						dst_end = dst.length;
					}

					if (expected_first_negations > 0) {
						for (unsigned int i = 0; i < expected_first_negations; i++) {
							hol_term* temp = hol_term::new_not(new_second_head);
							if (temp == nullptr) {
								free(*new_second_head); free(new_second_head);
								return false;
							}
							new_second_head = temp;
						}
						intersect<built_in_predicates>(dst, new_second_head, first_head);
					}
					free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);

					for (hol_term* new_head : dst) {
						if (!can_have_free_variables(*new_head))
							any_right_only = true;
					}

					if (!dst_outer.ensure_capacity(dst.length)) {
						free_all(dst);
						return false;
					}
					for (unsigned int i = 0; i < dst.length; i++)
						dst_outer[i] = &HOL_ZERO;
					dst_outer.length = dst.length;
					HOL_ZERO.reference_count += dst.length;

					if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
						/* gather the other quantifiers outside the head in `second` */
						array<hol_term*> excluded_trees(first_head->any.excluded_tree_count + 4);
						for (unsigned int i = 0; i < first_head->any.excluded_tree_count; i++) {
							excluded_trees[i] = first_head->any.excluded_trees[i];
							excluded_trees[i]->reference_count++;
						}
						excluded_trees.length = first_head->any.excluded_tree_count;
						excluded_trees[excluded_trees.length] = hol_term::new_any_right(hol_term::new_apply(
								hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
						if (excluded_trees[excluded_trees.length] == nullptr) {
							free_all(excluded_trees); free_all(dst); free_all(dst_outer);
							return false;
						}
						HOL_ANY.reference_count++;
						excluded_trees.length++;
						for (unsigned int i = second_inverter.outer.length - 1; i > 0; i--) {
							const hol_term* node = second_inverter.outer[i - 1];
							if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
								continue;
							if (!excluded_trees.ensure_capacity(excluded_trees.length + 1)) {
								free_all(excluded_trees); free_all(dst); free_all(dst_outer);
								return false;
							}

							if (node->type == hol_term_type::FOR_ALL) {
								excluded_trees[excluded_trees.length] = hol_term::new_any(hol_term::new_for_all(node->quantifier.variable, &HOL_ANY));
							} else if (node->type == hol_term_type::EXISTS) {
								excluded_trees[excluded_trees.length] = hol_term::new_any(hol_term::new_exists(node->quantifier.variable, &HOL_ANY));
							} else {
								excluded_trees[excluded_trees.length] = hol_term::new_any(hol_term::new_lambda(node->quantifier.variable, &HOL_ANY));
							}
							if (excluded_trees[excluded_trees.length] == nullptr) {
								free_all(excluded_trees); free_all(dst); free_all(dst_outer);
								return false;
							}
							HOL_ANY.reference_count++;
							excluded_trees.length++;
						}

						/* add excluded tree for when the head is not negated */
						if (!excluded_trees.ensure_capacity(excluded_trees.length + 1)) {
							free_all(excluded_trees); free_all(dst); free_all(dst_outer);
							return false;
						}
						excluded_trees[excluded_trees.length] = hol_term::new_any(hol_term::new_not(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY)));
						if (excluded_trees[excluded_trees.length] == nullptr) {
							free_all(excluded_trees); free_all(dst); free_all(dst_outer);
							return false;
						}
						HOL_ANY.reference_count++;
						excluded_trees.length++;

						for (unsigned int i = 0; i < dst.length; i++) {
							if (dst[i]->type == hol_term_type::ANY || dst[i]->type == hol_term_type::ANY_RIGHT || !can_have_free_variables(*dst[i])) continue;
							if ((!is_array && second != second_head) || (is_array && second_inverter.outer.length != 1)) {
								if (i < dst_start || i >= dst_end) continue;
								hol_term* new_outer = hol_term::new_any_right(dst_outer[i], excluded_trees.data + excluded_trees.length - 1, 1);
								if (new_outer == nullptr) {
									free_all(excluded_trees); free_all(dst); free_all(dst_outer);
									return false;
								}
								dst_outer[i] = new_outer;
								excluded_trees[excluded_trees.length - 1]->reference_count++;
							} else {
								unsigned int excluded_count = (i >= dst_start && i < dst_end) ? excluded_trees.length : (excluded_trees.length - 1);
								hol_term* new_outer = hol_term::new_any_right(dst_outer[i], excluded_trees.data, excluded_count);
								if (new_outer == nullptr) {
									free_all(excluded_trees); free_all(dst); free_all(dst_outer);
									return false;
								}
								dst_outer[i] = new_outer;
								for (unsigned int i = 0; i < excluded_count; i++)
									excluded_trees[i]->reference_count++;
							}
						}
						free_all(excluded_trees);
					} else {
						hol_term* excluded_negation = hol_term::new_any(hol_term::new_not(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY)));
						if (excluded_negation == nullptr) {
							free_all(dst); free_all(dst_outer);
							return false;
						}
						HOL_ANY.reference_count++;

						for (unsigned int i = dst_start; i < dst_end; i++) {
							hol_term* new_outer = hol_term::new_any_right(dst_outer[i], &excluded_negation, 1);
							if (new_outer == nullptr) {
								free(*excluded_negation); if (excluded_negation->reference_count == 0) free(excluded_negation);
								free_all(dst); free_all(dst_outer); return false;
							}
							dst_outer[i] = new_outer;
							excluded_negation->reference_count++;
						}
						free(*excluded_negation); if (excluded_negation->reference_count == 0) free(excluded_negation);
					}
					return (dst.length > 0);
				} else if (first_head->type == hol_term_type::EXISTS && first_head->quantifier.operand->type == hol_term_type::AND) {
					hol_term* first_head_operand = first_head->quantifier.operand;

					if ((!could_be_zero_negations && first_min_negations != expected_first_negations)
					 || (could_be_zero_negations && (first_min_negations != expected_first_negations && first_min_negations != 0)))
						return false;

					int conjunct_index = ConjunctIndex;
					if (ConjunctIndex < 0) conjunct_index += first_head_operand->array.length;
#if !defined(NDEBUG)
					if (conjunct_index == (int) first_predicate_index.index)
						fprintf(stderr, "invert_select_conjunct ERROR: `conjunct_index` and `first_predicate_index` are the same.\n");
#endif

					array<hol_term*> predicates(8);
					array<hol_term*> conjuncts(8);
					intersect<built_in_predicates>(predicates, first_head_operand->array.operands[first_predicate_index.index], second_head_operand->array.operands[0]);
					intersect<built_in_predicates>(conjuncts, first_head_operand->array.operands[conjunct_index], second_head_operand->array.operands[1]);
					unsigned int intersection_count = predicates.length * conjuncts.length;
					if (!dst.ensure_capacity(dst.length + intersection_count) || !dst_outer.ensure_capacity(dst_outer.length + intersection_count)) {
						for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					for (hol_term* conjunct : conjuncts) {
						static unsigned int ARGS[] = { (unsigned int) built_in_predicates::ARG1, (unsigned int) built_in_predicates::ARG2, (unsigned int) built_in_predicates::ARG3 };
						hol_term* arg_trees[array_length(ARGS)];
						unsigned int arg_tree_count = 0;
						for (unsigned int i = 0; i < array_length(ARGS); i++) {
							arg_trees[arg_tree_count] = hol_term::new_any(hol_term::new_apply(hol_term::new_constant(ARGS[i]), hol_term::new_variable(second_head->quantifier.variable)));
							if (arg_trees[arg_tree_count] == nullptr) {
								for (unsigned int j = 0; j < i; j++) { free(*arg_trees[j]); free(arg_trees[j]); }
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								return false;
							}
							if (has_intersection<built_in_predicates>(arg_trees[arg_tree_count], conjunct)) {
								arg_tree_count++;
							} else {
								free(*arg_trees[arg_tree_count]); free(arg_trees[arg_tree_count]);
							}
						}
						for (hol_term* predicate : predicates) {
							hol_term* conjunction;
							if (!new_hol_term(conjunction)) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
								return false;
							}
							conjunction->type = hol_term_type::AND;
							conjunction->reference_count = 1;
							conjunction->array.length = first_head_operand->array.length;
							conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first_head_operand->array.length);
							if (conjunction->array.operands == nullptr) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
								free(conjunction); return false;
							}
							bool skip = false;
							for (unsigned int i = 0; i < first_head_operand->array.length; i++) {
								if (i == first_predicate_index.index) {
									conjunction->array.operands[i] = predicate;
								} else if (i == (unsigned int) conjunct_index) {
									conjunction->array.operands[i] = conjunct;
								} else {
									for (unsigned int k = 0; k < arg_tree_count; k++) {
										if (has_intersection<built_in_predicates>(arg_trees[k], first_head_operand->array.operands[i])) {
											for (unsigned int j = 0; j < i; j++) {
												free(*conjunction->array.operands[i]);
												if (conjunction->array.operands[i]->reference_count == 0)
													free(conjunction->array.operands[i]);
											}
											free(conjunction); skip = true; break;
										}
									}
									if (skip) break;
									conjunction->array.operands[i] = first_head_operand->array.operands[i];
								}
								conjunction->array.operands[i]->reference_count++;
							}
							if (skip) continue;

							hol_term* new_head = hol_term::new_exists(first_head->quantifier.variable, conjunction);
							if (new_head == nullptr) {
								free(*conjunction); free(conjunction);
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
								free_all(dst); free_all(dst_outer); return false;
							}

							/* check if there `new_head` can have any free variables; if not, then set `any_right_only` to true */
							if (!can_have_free_variables(*new_head))
								any_right_only = true;

							if (could_be_zero_negations) {
								dst[dst.length++] = new_head;
								new_head->reference_count++;
							}
							if (expected_first_negations > 0) {
								for (unsigned int i = 0; i < expected_first_negations; i++) {
									hol_term* temp = hol_term::new_not(new_head);
									if (temp == nullptr) {
										for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
										for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
										free(*new_head); if (new_head->reference_count == 0) free(new_head);
										free_all(dst); free_all(dst_outer); return false;
									}
									new_head = temp;
								}
								dst[dst.length++] = new_head;
							}

							for (unsigned int j = 0; j < dst.length; j++) {
								if (old_second_head->type == hol_term_type::ANY || old_second_head->type == hol_term_type::ANY_RIGHT) {
									dst_outer[dst_outer.length] = hol_term::new_any_right(&HOL_ZERO, old_second_head->any.excluded_trees, old_second_head->any.excluded_tree_count);
									if (dst_outer[dst_outer.length] == nullptr) {
										for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
										for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
										free_all(dst); free_all(dst_outer); return false;
									}
									for (unsigned int i = 0; i < old_second_head->any.excluded_tree_count; i++)
										old_second_head->any.excluded_trees[i]->reference_count++;
								} else {
									dst_outer[dst_outer.length] = &HOL_ZERO;
								}
								HOL_ZERO.reference_count++;
								dst_outer.length++;
							}
						}
						for (unsigned int k = 0; k < arg_tree_count; k++) { free(*arg_trees[k]); free(arg_trees[k]); }
					}
					for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
					return (dst.length > 0);
				} else {
					fprintf(stderr, "invert_select_conjunct ERROR: Unexpected type of `first_head`.\n");
					return false;
				}
			} else {
				return false;
			}
		});
}

template<int_fast8_t ConjunctIndex, bool RemoveNegation = false>
inline bool invert_remove_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	auto on_remap_variables = [](hol_term* first_head, hol_term* second_head,
			apply_head_inverter& first_head_inverter, apply_head_inverter& second_head_inverter,
			array_map<const hol_term*, unsigned int>& second_variable_map, unsigned int max_variable,
			bool first_head_is_array, bool second_head_is_array)
	{
		if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
			second_head = second_head->any.included;

		while (first_head->type == hol_term_type::NOT)
			first_head = first_head->unary.operand;
		while (second_head->type == hol_term_type::NOT)
			second_head = second_head->unary.operand;

		if (second_head->type != hol_term_type::EXISTS)
			return false;
		if (first_head->type != hol_term_type::EXISTS)
			return true;

		hol_term* first_operand = first_head->quantifier.operand;
		hol_term* new_conjunction;
		if (first_operand->type == hol_term_type::ANY || first_operand->type == hol_term_type::ANY_RIGHT) {
			return true;
		} else if (first_operand->type == hol_term_type::ANY_ARRAY && first_operand->any_array.oper == hol_term_type::AND) {
			if (ConjunctIndex >= 0) {
				if (ConjunctIndex < first_operand->any_array.left.length) {
					new_conjunction = hol_term::new_any_array(hol_term_type::AND,
							first_operand->any_array.all, make_array_view(first_operand->any_array.any.operands, first_operand->any_array.any.length),
							make_excluded_array_view(first_operand->any_array.left.operands, first_operand->any_array.left.length, ConjunctIndex),
							make_array_view(first_operand->any_array.right.operands, first_operand->any_array.right.length));
					if (new_conjunction == nullptr)
						return false;
					new_conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.left.length; i++)
						new_conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.right.length; i++)
						new_conjunction->any_array.right.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.any.length; i++)
						new_conjunction->any_array.any.operands[i]->reference_count++;
				} else {
					new_conjunction = first_operand;
					first_operand->reference_count++;
				}
			} else {
				if (-ConjunctIndex <= first_operand->any_array.right.length) {
					new_conjunction = hol_term::new_any_array(hol_term_type::AND,
							first_operand->any_array.all, make_array_view(first_operand->any_array.any.operands, first_operand->any_array.any.length),
							make_array_view(first_operand->any_array.left.operands, first_operand->any_array.left.length),
							make_excluded_array_view(first_operand->any_array.right.operands, first_operand->any_array.right.length, first_operand->any_array.right.length + ConjunctIndex));
					if (new_conjunction == nullptr)
						return false;
					new_conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.left.length; i++)
						new_conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.right.length; i++)
						new_conjunction->any_array.right.operands[i]->reference_count++;
					for (unsigned int i = 0; i < new_conjunction->any_array.any.length; i++)
						new_conjunction->any_array.any.operands[i]->reference_count++;
				} else {
					new_conjunction = first_operand;
					first_operand->reference_count++;
				}
			}
		} else if (first_operand->type == hol_term_type::AND) {
			unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (unsigned int) (first_operand->array.length + ConjunctIndex);
			if (first_operand->array.length == 2) {
				new_conjunction = (index == 0) ? first_operand->array.operands[1] : first_operand->array.operands[0];
				new_conjunction->reference_count++;
			} else {
				new_conjunction = hol_term::new_and(make_excluded_array_view(first_operand->array.operands, first_operand->array.length, index));
				if (new_conjunction == nullptr)
					return false;
				for (unsigned int i = 0; i < new_conjunction->array.length; i++)
					new_conjunction->array.operands[i]->reference_count++;
			}
		} else {
			return false;
		}

		hol_term* new_first_head = hol_term::new_exists(first_head->quantifier.variable, new_conjunction);
		if (new_first_head == nullptr) {
			free(*new_conjunction); if (new_conjunction->reference_count == 0) free(new_conjunction);
			return false;
		}

		array<pair<hol_term*, variable_map>> intersection(2);
		intersect<built_in_predicates>(intersection, new_first_head, second_head);
		free(*new_first_head); if (new_first_head->reference_count == 0) free(new_first_head);
		if (intersection.length == 0) {
			return false;
		} else if (intersection.length != 1) {
			fprintf(stderr, "invert_remove_conjunct ERROR: Intersection module variable relabeling is not unique.\n");
			free_all(intersection); return false;
		}

		/* make sure the variables in `second_head` map to the correct variables in `first_head` */
		const variable_map& var_map = intersection[0].value;
		for (const auto& entry : var_map.scope_map) {
			if (!second_variable_map.contains(entry.key.src))
				continue;

			unsigned int target_var = 0;
			if (!get_target_variable(entry.key.src->quantifier.variable, target_var, entry.value, max_variable)) {
				free_all(intersection);
				return false;
			}

			if (!second_variable_map.ensure_capacity(second_variable_map.size + 1)) {
				free_all(intersection);
				return false;
			}
			unsigned int index = second_variable_map.index_of(entry.key.src);
			if (target_var == second_variable_map.values[index]) {
				/* `second_variable_map` should map `var` to itself */
				second_variable_map.remove_at(index);
			} else {
				/* `second_variable_map` should map `var` to `target_var` */
				second_variable_map.values[index] = target_var;
				if (index == second_variable_map.size) {
					second_variable_map.keys[index] = entry.key.src;
					second_variable_map.size++;
				}
			}
		}

		array_map<unsigned int, variable_set> free_variable_map(8);
		for (const auto& entry : var_map.free_variables) {
			const variable_set& set = entry.value;
			unsigned int src_variable = entry.key;
			if (!free_variable_map.ensure_capacity(free_variable_map.size + 1)) {
				free_all(intersection);
				return false;
			}
			unsigned int index = free_variable_map.index_of(src_variable);
			if (index < free_variable_map.size) {
				variable_set& new_set = *((variable_set*) alloca(sizeof(variable_map)));
				if (!intersect(new_set, set, free_variable_map.values[index])) {
					free_all(intersection);
					return false;
				}
				swap(free_variable_map.values[index], new_set);
				free(new_set);
			} else {
				free_variable_map.keys[free_variable_map.size] = src_variable;
				if (!init(free_variable_map.values[free_variable_map.size], set)) {
					free_all(intersection);
					return false;
				}
				free_variable_map.size++;
			}
		}
		free_all(intersection);

		if (first_head_inverter.outer.last() != first_head && !first_head_inverter.outer.add(first_head)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		} if (second_head_inverter.outer.last() != second_head && !second_head_inverter.outer.add(second_head)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		}

		if (!remap_scopes(free_variable_map, first_head_inverter, second_head_inverter, second_variable_map, max_variable, first_head_is_array, second_head_is_array)) {
			for (auto entry : second_variable_map) free(entry.value);
			return false;
		}
		return true;
	};

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, on_remap_variables,
		[second](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
	{
		hol_term* old_second_head = second_head;
		if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
			second_head = second_head->any.included;

		unsigned int first_min_negations = 0;
		unsigned int second_negations = 0;
		while (first_head->type == hol_term_type::NOT) {
			first_head = first_head->unary.operand;
			first_min_negations++;
		} while (second_head->type == hol_term_type::NOT) {
			second_head = second_head->unary.operand;
			second_negations++;
		}

		/* if `RemoveNegation` is true, `first_negations` should be equal to `second_negations + 1`,
			otherwise `first_negations` should be equal to `second_negations` */
		unsigned int expected_first_negations = (RemoveNegation ? second_negations + 1 : second_negations);
		if (first_min_negations > expected_first_negations)
			return false;

		if (second_head->type == hol_term_type::EXISTS) {
			hol_term* second_head_operand = second_head->quantifier.operand;
			hol_term* head_var = hol_term::new_variable(second_head->quantifier.variable);
			if (head_var == nullptr) return false;

			hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
			if (any_wide_scope == nullptr) {
				free(*head_var); free(head_var);
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
			if (duplicate_wide_scopes == nullptr) {
				free(*head_var); free(head_var);
				free(*any_wide_scope); free(any_wide_scope);
				return false;
			}

			hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
			if (excluded_universal == nullptr) {
				free(*head_var); free(head_var);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				return false;
			}
			any_wide_scope->reference_count++;

			static unsigned int ARGS[] = { (unsigned int) built_in_predicates::ARG1, (unsigned int) built_in_predicates::ARG2, (unsigned int) built_in_predicates::ARG3 };
			unsigned int excluded_tree_count = 7;
			hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term*) * (excluded_tree_count + array_length(ARGS)));
			if (excluded_trees == nullptr) {
				free(*head_var); free(head_var);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				free(*excluded_universal); free(excluded_universal);
				return false;
			}
			excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
			excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
			excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
			excluded_trees[3] = hol_term::new_any(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY));
			excluded_trees[4] = hol_term::new_apply(
					hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), head_var);
			excluded_trees[5] = duplicate_wide_scopes;
			excluded_trees[6] = excluded_universal;
			if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
			if (excluded_trees[4] != nullptr) { hol_non_head_constants<built_in_predicates>::increment_terms(); head_var->reference_count++; }
			if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr || excluded_trees[4] == nullptr) {
				if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
				if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
				if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
				if (excluded_trees[3] != nullptr) { free(*excluded_trees[3]); free(excluded_trees[3]); }
				free(*head_var); free(head_var);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				free(*excluded_universal); free(excluded_universal);
				return false;
			}
			free(*head_var);

			for (unsigned int i = 0; i < array_length(ARGS); i++) {
				excluded_trees[excluded_tree_count] = hol_term::new_any(hol_term::new_apply(hol_term::new_constant(ARGS[i]), hol_term::new_variable(second_head->quantifier.variable)));
				if (excluded_trees[excluded_tree_count] == nullptr) {
					for (unsigned int j = 0; j < excluded_tree_count; j++) { free(*excluded_trees[j]); free(excluded_trees[j]); }
					return false;
				}
				if (has_intersection<built_in_predicates>(excluded_trees[excluded_tree_count], second_head_operand)) {
					excluded_tree_count++;
				} else {
					free(*excluded_trees[excluded_tree_count]); free(excluded_trees[excluded_tree_count]);
				}
			}

			hol_term** second_head_operands;
			unsigned int second_head_length;
			if (second_head_operand->type == hol_term_type::AND) {
				second_head_operands = second_head_operand->array.operands;
				second_head_length = second_head_operand->array.length;
			} else if (second_head_operand->type != hol_term_type::ANY && second_head_operand->type != hol_term_type::ANY_RIGHT && second_head_operand->type != hol_term_type::ANY_ARRAY) {
				second_head_operands = &second_head_operand;
				second_head_length = 1;
			} else {
				for (unsigned int j = 0; j < excluded_tree_count; j++) { free(*excluded_trees[j]); free(excluded_trees[j]); }
				return false;
			}

			hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
			if (conjunct == nullptr) {
				for (unsigned int j = 0; j < excluded_tree_count; j++) { free(*excluded_trees[j]); free(excluded_trees[j]); }
				return false;
			}

			hol_term* conjunction;
			if (!new_hol_term(conjunction)) {
				free(*conjunct); free(conjunct);
				return false;
			}
			conjunction->type = hol_term_type::AND;
			conjunction->reference_count = 1;
			conjunction->array.length = second_head_length + 1;
			conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * (second_head_length + 1));
			if (conjunction->array.operands == nullptr) {
				free(*conjunct); free(conjunct);
				free(conjunction); return false;
			}
			int conjunct_index = ConjunctIndex;
			if (ConjunctIndex < 0) conjunct_index += second_head_length + 1;
			for (unsigned int i = 0; i < conjunction->array.length; i++) {
				if (i == (unsigned int) conjunct_index) {
					conjunction->array.operands[i] = conjunct;
				} else {
					unsigned int second_head_index = (i > (unsigned) conjunct_index) ? (i - 1) : i;
					conjunction->array.operands[i] = second_head_operands[second_head_index];
					conjunction->array.operands[i]->reference_count++;
				}
			}

			dst[dst.length] = hol_term::new_exists(second_head->quantifier.variable, conjunction);
			if (dst[dst.length] == nullptr) {
				free(*conjunction); free(conjunction);
				return false;
			}
			dst.length++;

			for (unsigned int i = 0; i < expected_first_negations; i++) {
				hol_term* temp = hol_term::new_not(dst.last());
				if (temp == nullptr) {
					free_all(dst);
					return false;
				}
				dst.last() = temp;
			}

			if (old_second_head->type == hol_term_type::ANY || old_second_head->type == hol_term_type::ANY_RIGHT) {
				dst_outer[dst_outer.length] = hol_term::new_any_right(&HOL_ZERO, old_second_head->any.excluded_trees, old_second_head->any.excluded_tree_count);
				if (dst_outer[dst_outer.length] == nullptr) {
					free_all(dst);
					return false;
				}
				for (unsigned int i = 0; i < old_second_head->any.excluded_tree_count; i++)
					old_second_head->any.excluded_trees[i]->reference_count++;
			} else if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT)
					&& (second == second_head || (second->type == hol_term_type::UNARY_APPLICATION && second->binary.left->type == hol_term_type::CONSTANT && second->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE && second->binary.right == second_head)))
			{
				any_right_only = false;
				dst_outer[dst_outer.length] = &HOL_ZERO;
			} else {
				dst_outer[dst_outer.length] = &HOL_ZERO;
			}
			HOL_ZERO.reference_count++;
			dst_outer.length++;
			return true;
		} else {
			return false;
		}
	});
}

template<int_fast8_t ConjunctIndex, bool IsSecondRoot = false>
inline bool invert_select_set_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*first, max_variable);

	unsigned int lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		lambda_variable = first->quantifier.variable;
	} else if (first->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), IsSecondRoot ? find_root : find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;

#if !defined(NDEBUG)
			if (first_head->type != hol_term_type::EXISTS || second_head->type != hol_term_type::EXISTS)
				fprintf(stderr, "invert_select_set_conjunct WARNING: Expected an existential quantification.\n");
#endif

			unsigned int set_variable = first_head->quantifier.variable;

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr)
				return false;
			HOL_ANY.reference_count++;

			hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (expected_conjunct == nullptr) {
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			if (ConjunctIndex >= 0) {
				dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
						make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(expected_conjunct, ConjunctIndex), second_head->quantifier.operand), make_array_view((hol_term**) nullptr, 0)));
				if (dst[dst.length] == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				second_head->quantifier.operand->reference_count++;
				expected_conjunct->reference_count += ConjunctIndex;
			} else {
				dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(second_head->quantifier.operand, make_repeated_array_view(expected_conjunct, (unsigned int) (-ConjunctIndex - 1)))));
				if (dst[dst.length] == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				second_head->quantifier.operand->reference_count++;
				expected_conjunct->reference_count += (unsigned int) (-ConjunctIndex - 1);
			}
			dst.length++;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;

			hol_term* right = (ConjunctIndex == -1) ? second_head->quantifier.operand : expected_conjunct;
			hol_term* wide_scope = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY);
			if (wide_scope == nullptr) {
				free_all(dst); free_all(dst_outer);
				return false;
			}
			if (has_intersection<built_in_predicates>(wide_scope, right))
				could_have_wide_scope = true;
			free(*wide_scope); free(wide_scope);

			return true;
		});
	if (!result) return false;

	for (unsigned int i = 0; i < inverse_count; i++) {
		if (inverse[i].root->type != hol_term_type::LAMBDA) {
			hol_term* new_inverse = hol_term::new_lambda(lambda_variable, inverse[i].root);
			if (new_inverse == nullptr) {
				for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
				free(inverse); return false;
			}
			inverse[i].root = new_inverse;
		}
	}
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_remove_set_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*first, max_variable);

	unsigned int first_lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		first_lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else if (first->type == hol_term_type::ANY) {
		first_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	unsigned int second_lambda_variable = 0;
	if (second->type == hol_term_type::LAMBDA) {
		second_lambda_variable = second->quantifier.variable;
		second = second->quantifier.operand;
	} else if (second->type == hol_term_type::ANY) {
		second_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(first_lambda_variable), predicative_head_finder<built_in_predicates>(second_lambda_variable), no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
				second_head = second_head->any.included;

#if !defined(NDEBUG)
			if (first_head->type != hol_term_type::EXISTS || second_head->type != hol_term_type::EXISTS)
				fprintf(stderr, "invert_remove_set_conjunct WARNING: Expected an existential quantification.\n");
#endif

			unsigned int set_variable = first_head->quantifier.variable;

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr)
				return false;
			HOL_ANY.reference_count++;

			hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (expected_conjunct == nullptr) {
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			hol_term* conjunction;
			hol_term* operand = second_head->quantifier.operand;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (ConjunctIndex >= 0) {
					conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all, make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
							make_included_array_view(operand->any_array.left.operands, operand->any_array.left.length, expected_conjunct, ConjunctIndex), make_array_view(operand->any_array.right.operands, operand->any_array.right.length));
				} else {
					unsigned int index = (unsigned int) (operand->any_array.right.length + ConjunctIndex);
					conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all, make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
							make_array_view(operand->any_array.left.operands, operand->any_array.left.length), make_included_array_view(operand->any_array.right.operands, operand->any_array.right.length, expected_conjunct, index));
				}
				if (conjunction == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				conjunction->any_array.all->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
					conjunction->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
					conjunction->any_array.right.operands[i]->reference_count++;
				for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
					conjunction->any_array.any.operands[i]->reference_count++;
			} else {
				unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (unsigned int) (ConjunctIndex + operand->array.length);
				conjunction = hol_term::new_and(make_included_array_view(operand->array.operands, operand->array.length, expected_conjunct, index));
				if (conjunction == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				for (unsigned int i = 0; i < operand->array.length; i++)
					operand->array.operands[i]->reference_count++;
			}

			dst[dst.length] = hol_term::new_exists(set_variable, conjunction);
			if (dst[dst.length] == nullptr) {
				free(*conjunction); free(conjunction);
				return false;
			}
			dst.length++;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
	if (!result) return false;

	for (unsigned int i = 0; i < inverse_count; i++) {
		if (inverse[i].root->type != hol_term_type::LAMBDA) {
			hol_term* new_inverse = hol_term::new_lambda(first_lambda_variable, inverse[i].root);
			if (new_inverse == nullptr) {
				for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
				free(inverse); return false;
			}
			inverse[i].root = new_inverse;
		}
	}
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_select_conjunct_in_set(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*first, max_variable);

	unsigned int lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		lambda_variable = first->quantifier.variable;
	} else if (first->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS || second_head->quantifier.operand->type != hol_term_type::AND)
				fprintf(stderr, "invert_select_conjunct_in_set WARNING: Expected an existentially-quantified conjunction.\n");
#endif

			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;

			unsigned int set_variable;
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				set_variable = ++max_variable;
			} else {
#if !defined(NDEBUG)
				if (first_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_select_conjunct_in_set WARNING: Expected an existential quantification.\n");
#endif
				set_variable = first_head->quantifier.variable;
			}

			unsigned int element_variable = second_head->quantifier.variable;
			hol_term* head_var = hol_term::new_variable(element_variable);
			if (head_var == nullptr) return false;
			constexpr unsigned int excluded_tree_count = 5;
			hol_term* excluded_trees[excluded_tree_count];
			excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
			excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
			excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
			excluded_trees[3] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
			excluded_trees[4] = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
			if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
			if (excluded_trees[4] != nullptr) { HOL_ANY.reference_count++; }
			if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr || excluded_trees[4] == nullptr) {
				if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
				if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
				if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
				if (excluded_trees[3] != nullptr) { free(*excluded_trees[3]); free(excluded_trees[3]); }
				free(*head_var); free(head_var);
				return false;
			}
			free(*head_var);

			hol_term* expected_conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
			if (expected_conjunct == nullptr) {
				free(*excluded_trees[0]); free(excluded_trees[0]);
				free(*excluded_trees[1]); free(excluded_trees[1]);
				free(*excluded_trees[2]); free(excluded_trees[2]);
				free(*excluded_trees[3]); free(excluded_trees[3]);
				free(*excluded_trees[4]); free(excluded_trees[4]);
				return false;
			}

			hol_term* conjunction = nullptr;
			if (ConjunctIndex >= 0) {
				conjunction = hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
						make_array_view(&second_head->quantifier.operand->array.operands[0], 1),
						make_appended_array_view(make_repeated_array_view(expected_conjunct, ConjunctIndex), second_head->quantifier.operand->array.operands[1]),
						make_array_view((hol_term**) nullptr, 0));
			} else if (ConjunctIndex < 0) {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				conjunction = hol_term::new_any_array(hol_term_type::AND, expected_conjunct,
						make_array_view(&second_head->quantifier.operand->array.operands[0], 1),
						make_array_view((hol_term**) nullptr, 0),
						make_prepended_array_view(second_head->quantifier.operand->array.operands[1], make_repeated_array_view(expected_conjunct, index)));
			}
			if (conjunction == nullptr) {
				free(*expected_conjunct); free(expected_conjunct);
				return false;
			}
			conjunction->any_array.all->reference_count++;
			for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
				conjunction->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
				conjunction->any_array.right.operands[i]->reference_count++;
			for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
				conjunction->any_array.any.operands[i]->reference_count++;
			free(*expected_conjunct); if (expected_conjunct->reference_count == 0) free(expected_conjunct);

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr) {
				free(*conjunction); free(conjunction);
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* set_definition = hol_term::new_any_right(hol_term::new_lambda(element_variable, conjunction), &excluded_quantifier, 1);
			if (set_definition == nullptr) {
				free(*conjunction); free(conjunction);
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			hol_term* outer_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (outer_conjunct == nullptr) {
				free(*set_definition); free(set_definition);
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}
			excluded_quantifier->reference_count++;

			hol_term* new_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, outer_conjunct,
					make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition, 1), make_array_view((hol_term**) nullptr, 0)));
			if (new_head == nullptr) {
				free(*set_definition); free(set_definition);
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}
			could_have_wide_scope = true;

			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
	if (!result) return false;

	for (unsigned int i = 0; i < inverse_count; i++) {
		if (inverse[i].root->type != hol_term_type::LAMBDA) {
			hol_term* new_inverse = hol_term::new_lambda(lambda_variable, inverse[i].root);
			if (new_inverse == nullptr) {
				for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
				free(inverse); return false;
			}
			inverse[i].root = new_inverse;
		}
	}
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_remove_conjunct_in_set(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*first, max_variable);

	unsigned int first_lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		first_lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else if (first->type == hol_term_type::ANY) {
		first_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	unsigned int second_lambda_variable = 0;
	if (second->type == hol_term_type::LAMBDA) {
		second_lambda_variable = second->quantifier.variable;
		second = second->quantifier.operand;
	} else if (second->type == hol_term_type::ANY) {
		second_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(first_lambda_variable), predicative_head_finder<built_in_predicates>(second_lambda_variable), no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
				second_head = second_head->any.included;

#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS)
				fprintf(stderr, "invert_remove_conjunct_in_set WARNING: Expected an existential quantification.\n");
#endif

			hol_term* operand = second_head->quantifier.operand;
			unsigned int set_variable = second_head->quantifier.variable;

			hol_term* left;
			if (operand->type == hol_term_type::AND) {
				left = operand->array.operands[0];
			} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.left.length != 0) {
				left = operand->any_array.left.operands[0];
			} else {
				return false;
			}

			hol_term* inner_operand;
			unsigned int element_variable;
			if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
				inner_operand = left->any.included->quantifier.operand;
				element_variable = left->any.included->quantifier.variable;
			} else if (left->type == hol_term_type::EQUALS) {
				inner_operand = left->binary.right->quantifier.operand;
				element_variable = left->binary.right->quantifier.variable;
			} else {
				inner_operand = left->ternary.third->quantifier.operand;
				element_variable = left->ternary.third->quantifier.variable;
			}

			hol_term* head_var = hol_term::new_variable(element_variable);
			if (head_var == nullptr) return false;
			unsigned int excluded_tree_count = (inner_operand->type == hol_term_type::TRUE) ? 5 : 4;
			hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term*) * excluded_tree_count);
			excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
			excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
			excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
			excluded_trees[3] = hol_term::new_any(hol_term::new_lambda(element_variable, &HOL_ANY));
			if (inner_operand->type == hol_term_type::TRUE)
				excluded_trees[4] = hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(&HOL_ANY, 2), make_array_view((hol_term**) nullptr, 0));
			if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
			if (inner_operand->type == hol_term_type::TRUE && excluded_trees[4] != nullptr) { HOL_ANY.reference_count += 3; }
			if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr
			 || (inner_operand->type == hol_term_type::TRUE && excluded_trees[4] == nullptr))
			{
				if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
				if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
				if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
				if (excluded_trees[3] != nullptr) { free(*excluded_trees[3]); free(excluded_trees[3]); }
				free(*head_var); free(head_var);
				return false;
			}
			free(*head_var);

			hol_term* expected_conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
			if (expected_conjunct == nullptr) {
				free(*excluded_trees[0]); free(excluded_trees[0]);
				free(*excluded_trees[1]); free(excluded_trees[1]);
				free(*excluded_trees[2]); free(excluded_trees[2]);
				free(*excluded_trees[3]); free(excluded_trees[3]);
				if (inner_operand->type == hol_term_type::TRUE) {
					free(*excluded_trees[4]); free(excluded_trees[4]);
				}
				return false;
			}

			hol_term* conjunction = nullptr;
			if (inner_operand->type == hol_term_type::AND) {
				unsigned int index = (unsigned int) ((ConjunctIndex < 0) ? (ConjunctIndex + inner_operand->array.length + 1) : ConjunctIndex);
				conjunction = hol_term::new_and(make_included_array_view(inner_operand->array.operands, inner_operand->array.length, expected_conjunct, index));
				if (conjunction == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				for (unsigned int i = 0; i < inner_operand->array.length; i++)
					inner_operand->array.operands[i]->reference_count++;
			} else if (inner_operand->type == hol_term_type::ANY_ARRAY) {
				/* TODO: we can avoid constructing `expected_conjunct` in this case */
				if (ConjunctIndex >= 0 && (unsigned int) ConjunctIndex < inner_operand->any_array.left.length) {
					conjunction = hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all,
							make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
							make_included_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length, expected_conjunct, ConjunctIndex),
							make_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length));
					if (conjunction == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return false;
					}
					conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
						conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
						conjunction->any_array.right.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
						conjunction->any_array.any.operands[i]->reference_count++;
					free(*expected_conjunct); if (expected_conjunct->reference_count == 0) free(expected_conjunct);
				} else if (ConjunctIndex < 0 && -ConjunctIndex < inner_operand->any_array.right.length + 1) {
					unsigned int index = (unsigned int) (inner_operand->any_array.right.length + ConjunctIndex + 1);
					conjunction = hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all,
							make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
							make_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length),
							make_included_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length, expected_conjunct, index));
					if (conjunction == nullptr) {
						free(*expected_conjunct); free(expected_conjunct);
						return false;
					}
					conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
						conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
						conjunction->any_array.right.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
						conjunction->any_array.any.operands[i]->reference_count++;
					free(*expected_conjunct); if (expected_conjunct->reference_count == 0) free(expected_conjunct);
				} else {
					free(*expected_conjunct); free(expected_conjunct);
					conjunction = inner_operand;
					conjunction->reference_count++;
				}
			} else if (inner_operand->type == hol_term_type::TRUE) {
				conjunction = expected_conjunct;
			} else {
				unsigned int index = (unsigned int) ((ConjunctIndex < 0) ? (ConjunctIndex + inner_operand->array.length + 1) : ConjunctIndex);
				if (index == 0)
					conjunction = hol_term::new_and(expected_conjunct, inner_operand);
				else conjunction = hol_term::new_and(inner_operand, expected_conjunct);
				if (conjunction == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				inner_operand->reference_count++;
			}

			hol_term* set_definition;
			if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
				set_definition = hol_term::new_any_right(hol_term::new_lambda(element_variable, conjunction), left->any.excluded_trees, left->any.excluded_tree_count);
				if (set_definition == nullptr) {
					free(*conjunction); free(conjunction);
					return false;
				}
				for (unsigned int i = 0; i < left->any.excluded_tree_count; i++)
					left->any.excluded_trees[i]->reference_count++;
			} else if (left->type == hol_term_type::EQUALS) {
				set_definition = hol_term::new_equals(left->binary.left, hol_term::new_lambda(element_variable, conjunction));
				if (set_definition == nullptr) {
					free(*conjunction); free(conjunction);
					return false;
				}
				left->binary.left->reference_count++;
			} else {
				set_definition = hol_term::new_apply(left->ternary.first, left->ternary.second, hol_term::new_lambda(element_variable, conjunction));
				if (set_definition == nullptr) {
					free(*conjunction); free(conjunction);
					return false;
				}
				left->ternary.first->reference_count++;
				left->ternary.second->reference_count++;
			}

			if (operand->type == hol_term_type::AND) {
				dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_and(make_prepended_array_view(set_definition,
						make_array_view(operand->array.operands + 1, operand->array.length - 1))));
				if (dst[dst.length] == nullptr) {
					free(*set_definition); free(set_definition);
					return false;
				}
				for (unsigned int i = 1; i < operand->array.length; i++)
					operand->array.operands[i]->reference_count++;
			} else {
				dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
						make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
						make_prepended_array_view(set_definition, make_array_view(operand->any_array.left.operands + 1, operand->any_array.left.length - 1)),
						make_array_view(operand->any_array.right.operands, operand->any_array.right.length)));
				if (dst[dst.length] == nullptr) {
					free(*set_definition); free(set_definition);
					return false;
				}
				operand->any_array.all->reference_count++;
				for (unsigned int i = 0; i < operand->any_array.any.length; i++)
					operand->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 1; i < operand->any_array.left.length; i++)
					operand->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < operand->any_array.right.length; i++)
					operand->any_array.right.operands[i]->reference_count++;
			}
			dst.length++;

			dst_outer[dst_outer.length] = hol_term::new_any_right(&HOL_ZERO);
			if (dst_outer[dst_outer.length] == nullptr) {
				free_all(dst);
				return false;
			}
			HOL_ZERO.reference_count++;
			dst_outer.length++;
			return true;
		});
	if (!result) return false;

	for (unsigned int i = 0; i < inverse_count; i++) {
		hol_term* new_inverse = hol_term::new_lambda(first_lambda_variable, inverse[i].root);
		if (new_inverse == nullptr) {
			for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
			free(inverse); return false;
		}
		inverse[i].root = new_inverse;
	}
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_select_subset_in_set(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*first, max_variable);

	unsigned int first_lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		first_lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else if (first->type == hol_term_type::ANY) {
		first_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	unsigned int second_lambda_variable = 0;
	if (second->type == hol_term_type::LAMBDA) {
		second_lambda_variable = second->quantifier.variable;
		second = second->quantifier.operand;
	} else if (second->type == hol_term_type::ANY) {
		second_lambda_variable = ++max_variable;
	} else {
		return false;
	}

	auto on_remap_variables = [](hol_term* first_head, hol_term* second_head,
			apply_head_inverter& first_head_inverter, apply_head_inverter& second_head_inverter,
			array_map<const hol_term*, unsigned int>& second_variable_map, unsigned int& max_variable,
			bool first_head_is_array, bool second_head_is_array)
	{
		if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
			first_head = first_head->any.included;
		if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
			second_head = second_head->any.included;

#if !defined(NDEBUG)
		if (second_head->type != hol_term_type::EXISTS || (second_head->quantifier.operand->type != hol_term_type::ANY_ARRAY && second_head->quantifier.operand->type != hol_term_type::AND))
			fprintf(stderr, "invert_select_subset_in_set WARNING: Expected an existentially-quantified conjunction.\n");
#endif

		unsigned int first_superset_variable;
		if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
			first_superset_variable = 0;
		} else {
#if !defined(NDEBUG)
			if (first_head->type != hol_term_type::EXISTS)
				fprintf(stderr, "invert_select_subset_in_set WARNING: Expected an existential quantification.\n");
#endif
			hol_term* left = nullptr;
			hol_term* operand = first_head->quantifier.operand;
			if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
				first_superset_variable = 0;
			} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
				if (operand->any_array.left.length == 0) {
					left = operand->any_array.all;
				} else {
					left = operand->any_array.left.operands[0];
				}
			} else if (operand->type == hol_term_type::AND) {
				left = operand->array.operands[0];
			} else {
				return false;
			}

			hol_term* set_definition = nullptr;
			if (left != nullptr) {
				if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
					first_superset_variable = 0;
				} else if (left->type == hol_term_type::EQUALS) {
					set_definition = left->binary.right;
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = left->ternary.third;
				} else {
					return false;
				}
			}

			hol_term* inner_operand = nullptr;
			if (set_definition != nullptr) {
				if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
					first_superset_variable = 0;
				} else if (set_definition->type == hol_term_type::ANY_QUANTIFIER && has_intersection(set_definition->any_quantifier.quantifier, hol_quantifier_type::LAMBDA)) {
					inner_operand = set_definition->any_quantifier.operand;
				} else if (set_definition->type == hol_term_type::LAMBDA) {
					inner_operand = set_definition->any_quantifier.operand;
				} else {
					return false;
				}
			}

			hol_term* inner_conjunct = nullptr;
			if (inner_operand != nullptr) {
				if (inner_operand->type == hol_term_type::ANY || inner_operand->type == hol_term_type::ANY_RIGHT) {
					first_superset_variable = 0;
				} else if (inner_operand->type == hol_term_type::ANY_ARRAY && inner_operand->any_array.oper == hol_term_type::AND) {
					if (ConjunctIndex >= 0){
						if (ConjunctIndex < inner_operand->any_array.left.length)
							inner_conjunct = inner_operand->any_array.left.operands[ConjunctIndex];
						else inner_conjunct = inner_operand->any_array.all;
					} else {
						if (-ConjunctIndex <= inner_operand->any_array.right.length)
							inner_conjunct = inner_operand->any_array.right.operands[inner_operand->any_array.right.length + ConjunctIndex];
						else inner_conjunct = inner_operand->any_array.all;
					}
				} else if (inner_operand->type == hol_term_type::AND) {
					if (ConjunctIndex >= 0)
						inner_conjunct = inner_operand->array.operands[ConjunctIndex];
					else inner_conjunct = inner_operand->array.operands[inner_operand->array.length + ConjunctIndex];
				} else {
					inner_conjunct = inner_operand;
				}
			}

			if (inner_conjunct != nullptr) {
				if (inner_conjunct->type == hol_term_type::ANY || inner_conjunct->type == hol_term_type::ANY_RIGHT) {
					first_superset_variable = 0;
				} else if (inner_conjunct->type == hol_term_type::UNARY_APPLICATION) {
					hol_term* superset_var = inner_conjunct->binary.left;
					if (superset_var->type == hol_term_type::VARIABLE) {
						first_superset_variable = superset_var->variable;
					} else {
						first_superset_variable = 0;
					}
				} else {
					return false;
				}
			}
		}

		if (first_superset_variable != 0) {
			/* make sure `second_superset_variable` maps to `first_superset_variable` */
			for (unsigned int i = 0; i < second_variable_map.size; i++) {
				if (second_variable_map.keys[i] == second_head) {
					second_variable_map.values[i] = first_superset_variable;
					if (second_head->quantifier.variable == first_superset_variable)
						second_variable_map.remove_at(i);
					break;
				}
			}
		}
		return true;
	};

	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(first_lambda_variable), predicative_head_finder<built_in_predicates>(second_lambda_variable), on_remap_variables,
		[first_lambda_variable,second_lambda_variable](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
				second_head = second_head->any.included;

#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS || (second_head->quantifier.operand->type != hol_term_type::ANY_ARRAY && second_head->quantifier.operand->type != hol_term_type::AND))
				fprintf(stderr, "invert_select_subset_in_set WARNING: Expected an existentially-quantified conjunction.\n");
#endif

			hol_term* left = nullptr;
			unsigned int first_set_variable, first_element_variable;
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				first_set_variable = ++max_variable;
				first_element_variable = 0;
			} else {
#if !defined(NDEBUG)
				if (first_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_select_subset_in_set WARNING: Expected an existential quantification.\n");
#endif
				first_set_variable = first_head->quantifier.variable;
				hol_term* operand = first_head->quantifier.operand;

				if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
					first_element_variable = 0;
				} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
					if (operand->any_array.left.length == 0) {
						left = operand->any_array.all;
					} else {
						left = operand->any_array.left.operands[0];
					}
				} else if (operand->type == hol_term_type::AND) {
					left = operand->array.operands[0];
				} else {
					return false;
				}

				hol_term* set_definition = nullptr;
				if (left != nullptr) {
					if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
						first_element_variable = 0;
					} else if (left->type == hol_term_type::BINARY_APPLICATION) {
						set_definition = left->ternary.third;
					} else {
						return false;
					}
				}

				if (set_definition != nullptr) {
					if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
						first_element_variable = 0;
					} else if (set_definition->type == hol_term_type::ANY_QUANTIFIER && has_intersection(set_definition->any_quantifier.quantifier, hol_quantifier_type::LAMBDA)) {
						first_element_variable = 0;
					} else if (set_definition->type == hol_term_type::LAMBDA) {
						first_element_variable = set_definition->quantifier.variable;
					} else {
						return false;
					}
				}

				if (first_element_variable == 0) {
					/* try to get the variable from the right conjunct in the set scope */
					hol_term* right = nullptr;
					if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
						if (operand->any_array.right.length == 0) {
							right = operand->any_array.all;
						} else {
							right = operand->any_array.right.operands[operand->any_array.right.length - 1];
						}
					} else if (operand->type == hol_term_type::AND) {
						right = operand->array.operands[operand->array.length - 1];
					}

					if (right != nullptr) {
						if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
							right = right->any.included;
						if (right->type == hol_term_type::NOT)
							right = right->unary.operand;

						if (right->type == hol_term_type::FOR_ALL && right->quantifier.operand->type == hol_term_type::IF_THEN
							&& right->quantifier.operand->binary.right->type == hol_term_type::UNARY_APPLICATION
							&& right->quantifier.operand->binary.right->binary.left->type == hol_term_type::VARIABLE
							&& right->quantifier.operand->binary.right->binary.left->variable == first_lambda_variable
							&& right->quantifier.operand->binary.right->binary.right->type == hol_term_type::VARIABLE)
						{
							first_element_variable = right->quantifier.operand->binary.right->binary.right->variable;
						} else if (right->type == hol_term_type::EXISTS && right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.length == 2
								&& right->quantifier.operand->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
								&& right->quantifier.operand->array.operands[1]->binary.left->type == hol_term_type::VARIABLE
								&& right->quantifier.operand->array.operands[1]->binary.left->variable == first_lambda_variable
								&& right->quantifier.operand->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
						{
							first_element_variable = right->quantifier.operand->array.operands[1]->binary.right->variable;
						} else if (right->type == hol_term_type::AND && right->array.length == 2
								&& right->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
								&& right->array.operands[1]->binary.left->type == hol_term_type::VARIABLE
								&& right->array.operands[1]->binary.left->variable == first_lambda_variable
								&& right->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
						{
							first_element_variable = right->array.operands[1]->binary.right->variable;
						} else if (right->type == hol_term_type::UNARY_APPLICATION
								&& right->binary.left->type == hol_term_type::VARIABLE
								&& right->binary.left->variable == first_lambda_variable
								&& right->binary.right->type == hol_term_type::VARIABLE)
						{
							first_element_variable = right->array.operands[1]->binary.right->variable;
						} else {
							return false;
						}
					}
				}
			}

			if (first_element_variable == 0) {
				/* we couldn't find the first element variable from `first_head`, so get it from `second_head` */
				hol_term* operand = second_head->quantifier.operand;

				hol_term* right;
				if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
					right = nullptr;
				} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
					if (operand->any_array.right.length == 0) {
						right = operand->any_array.all;
					} else {
						right = operand->any_array.right.operands[operand->any_array.right.length - 1];
					}
				} else if (operand->type == hol_term_type::AND) {
					right = operand->array.operands[operand->array.length - 1];
				} else {
					return false;
				}

				if (right != nullptr) {
					if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
						right = right->any.included;
					if (right->type == hol_term_type::NOT)
						right = right->unary.operand;

					if (right->type == hol_term_type::FOR_ALL && right->quantifier.operand->type == hol_term_type::IF_THEN
						&& right->quantifier.operand->binary.right->type == hol_term_type::UNARY_APPLICATION
						&& right->quantifier.operand->binary.right->binary.left->type == hol_term_type::VARIABLE
						&& right->quantifier.operand->binary.right->binary.left->variable == second_lambda_variable
						&& right->quantifier.operand->binary.right->binary.right->type == hol_term_type::VARIABLE)
					{
						first_element_variable = right->quantifier.operand->binary.right->binary.right->variable;
					} else if (right->type == hol_term_type::EXISTS && right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.length == 2
							&& right->quantifier.operand->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
							&& right->quantifier.operand->array.operands[1]->binary.left->type == hol_term_type::VARIABLE
							&& right->quantifier.operand->array.operands[1]->binary.left->variable == second_lambda_variable
							&& right->quantifier.operand->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
					{
						first_element_variable = right->quantifier.operand->array.operands[1]->binary.right->variable;
					} else if (right->type == hol_term_type::AND && right->array.length == 2
							&& right->array.operands[1]->type == hol_term_type::UNARY_APPLICATION
							&& right->array.operands[1]->binary.left->type == hol_term_type::VARIABLE
							&& right->array.operands[1]->binary.left->variable == second_lambda_variable
							&& right->array.operands[1]->binary.right->type == hol_term_type::VARIABLE)
					{
						first_element_variable = right->array.operands[1]->binary.right->variable;
					} else if (right->type == hol_term_type::UNARY_APPLICATION
							&& right->binary.left->type == hol_term_type::VARIABLE
							&& right->binary.left->variable == second_lambda_variable
							&& right->binary.right->type == hol_term_type::VARIABLE)
					{
						first_element_variable = right->array.operands[1]->binary.right->variable;
					} else {
						return false;
					}
				}
			}

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(first_set_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr)
				return false;
			HOL_ANY.reference_count++;

			hol_term* wildcard = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (wildcard == nullptr) {
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			unsigned int second_set_variable = second_head->quantifier.variable;
			hol_term* subset_term = hol_term::new_apply(hol_term::new_variable(second_set_variable), hol_term::new_variable(first_element_variable));
			if (subset_term == nullptr) {
				free(*wildcard); free(wildcard);
				return false;
			}

			hol_term* new_inner_operand;
			if (ConjunctIndex >= 0) {
				new_inner_operand = hol_term::new_any_array(hol_term_type::AND, wildcard, make_array_view((hol_term**) nullptr, 0),
						make_appended_array_view(make_repeated_array_view(wildcard, ConjunctIndex), subset_term),
						make_array_view((hol_term**) nullptr, 0));
				if (new_inner_operand == nullptr) {
					free(*wildcard); free(wildcard);
					free(*subset_term); free(subset_term);
					return false;
				}
				wildcard->reference_count += 1 + ConjunctIndex;
			} else {
				new_inner_operand = hol_term::new_any_array(hol_term_type::AND, wildcard,
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0),
						make_prepended_array_view(subset_term, make_repeated_array_view(wildcard, -ConjunctIndex - 1)));
				if (new_inner_operand == nullptr) {
					free(*wildcard); free(wildcard);
					free(*subset_term); free(subset_term);
					return false;
				}
				wildcard->reference_count += -ConjunctIndex;
			}

			hol_term* new_left = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SUBSET),
					hol_term::new_variable(first_set_variable), hol_term::new_lambda(first_element_variable, new_inner_operand));
			if (new_left == nullptr) {
				free(*new_inner_operand); free(new_inner_operand);
				free(*wildcard); if (wildcard->reference_count == 0) free(wildcard);
				return false;
			}

			hol_term* new_head = hol_term::new_exists(first_set_variable, hol_term::new_any_array(hol_term_type::AND,
					wildcard, make_array_view((hol_term**) nullptr, 0), make_array_view(&new_left, 1), make_array_view((hol_term**) nullptr, 0)));
			if (new_head == nullptr) {
				free(*new_left); free(new_left);
				free(*wildcard); if (wildcard->reference_count == 0) free(wildcard);
				return false;
			}
			HOL_ANY.reference_count++;
			dst[dst.length++] = new_head;

			hol_term* new_right = hol_term::new_any_right_only(&HOL_ZERO);
			if (new_right == nullptr)
				return false;
			HOL_ZERO.reference_count++;

			hol_term* second_operand = second_head->quantifier.operand;
			if (second_operand->type == hol_term_type::ANY_ARRAY && second_operand->any_array.oper == hol_term_type::AND) {
				dst_outer[dst_outer.length] = hol_term::new_any_array(second_operand->any_array.oper, second_operand->any_array.all,
						make_array_view(second_operand->any_array.any.operands, second_operand->any_array.any.length),
						make_array_view(second_operand->any_array.left.operands, second_operand->any_array.left.length),
						make_appended_array_view(make_array_view(second_operand->any_array.right.operands, second_operand->any_array.right.length - 1), new_right));
				if (dst_outer[dst_outer.length] == nullptr) {
					free(*new_right); free(new_right);
					return false;
				}
				dst_outer.length++;
				second_operand->any_array.all->reference_count++;
				for (unsigned int i = 0; i < second_operand->any_array.any.length; i++)
					second_operand->any_array.any.operands[i]->reference_count++;
				for (unsigned int i = 0; i < second_operand->any_array.left.length; i++)
					second_operand->any_array.left.operands[i]->reference_count++;
				for (unsigned int i = 0; i < second_operand->any_array.right.length - 1; i++)
					second_operand->any_array.right.operands[i]->reference_count++;
			} else {
				dst_outer[dst_outer.length] = hol_term::new_exists(second_set_variable, hol_term::new_and(make_appended_array_view(make_array_view(second_operand->array.operands, second_operand->array.length - 1), new_right)));
				if (dst_outer[dst_outer.length] == nullptr) {
					free(*new_right); free(new_right);
					return false;
				}
				dst_outer.length++;
				for (unsigned int i = 0; i < second_operand->array.length - 1; i++)
					second_operand->array.operands[i]->reference_count++;
			}
			return true;
		});
	if (!result) return false;

	for (unsigned int i = 0; i < inverse_count; i++) {
		hol_term* new_inverse = hol_term::new_lambda(first_lambda_variable, inverse[i].root);
		if (new_inverse == nullptr) {
			for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
			free(inverse); return false;
		}
		inverse[i].root = new_inverse;
	}
	return true;
}


inline bool invert_remove_inverse(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;

				unsigned int second_index = 0;
				if (second_predicate_index.position == head_position::LEFT)
					second_index = second_predicate_index.index;
				else if (second_predicate_index.position == head_position::RIGHT)
					second_index = second_head_operand->array.length - second_predicate_index.index - 1;
				hol_term* second_predicate = second_head_operand->array.operands[second_index];
				array<hol_term*> new_predicates(8);
				hol_term* new_predicate = hol_term::new_apply(
					hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::INVERSE), &HOL_ANY),
					hol_term::new_variable(second_head->quantifier.variable));
				if (new_predicate == nullptr) return false;
				HOL_ANY.reference_count++;
				intersect<built_in_predicates>(new_predicates, second_predicate, new_predicate);
				free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);

#if !defined(NDEBUG)
				if (new_predicates.length == 0)
					fprintf(stderr, "invert_remove_inverse WARNING: `new_predicates` is empty.\n");
#endif

				for (hol_term* new_predicate : new_predicates) {
					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					conjunction->type = hol_term_type::AND;
					conjunction->reference_count = 1;
					conjunction->array.length = second_head_operand->array.length;
					conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * second_head_operand->array.length);
					if (conjunction->array.operands == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(conjunction); return false;
					}
					for (unsigned int i = 0; i < conjunction->array.length; i++) {
						if (i == second_index) {
							conjunction->array.operands[i] = new_predicate;
						} else {
							conjunction->array.operands[i] = second_head_operand->array.operands[i];
						}
						conjunction->array.operands[i]->reference_count++;
					}

					hol_term* new_second_head = hol_term::new_exists(second_head->quantifier.variable, conjunction);
					if (new_second_head == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunction); free(conjunction);
						return false;
					}
					intersect<built_in_predicates>(dst, new_second_head, first_head);
					free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);

					if (!dst_outer.ensure_capacity(dst.length)) {
						free_all(dst); free_all(new_predicates);
						return false;
					}
					for (unsigned int i = 0; i < dst.length; i++)
						dst_outer[i] = &HOL_ZERO;
					dst_outer.length = dst.length;
					HOL_ZERO.reference_count += dst.length;
				}
				free_all(new_predicates);
				return (dst.length > 0);
			} else if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
				hol_term* second_head_operand = second_head->quantifier.operand;

				hol_term* second_predicate = nullptr;
				if (second_predicate_index.position == head_position::LEFT)
					second_predicate = second_head_operand->any_array.left.operands[second_predicate_index.index];
				else if (second_predicate_index.position == head_position::RIGHT)
					second_predicate = second_head_operand->any_array.right.operands[second_head_operand->any_array.right.length - second_predicate_index.index - 1];
				else if (second_predicate_index.position == head_position::ANY)
					second_predicate = second_head_operand->any_array.any.operands[second_predicate_index.index];

				array<hol_term*> new_predicates(8);
				hol_term* new_predicate = hol_term::new_apply(
					hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::INVERSE), &HOL_ANY),
					hol_term::new_variable(second_head->quantifier.variable));
				if (new_predicate == nullptr) return false;
				HOL_ANY.reference_count++;
				intersect<built_in_predicates>(new_predicates, second_predicate, new_predicate);
				free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);

#if !defined(NDEBUG)
				if (new_predicates.length == 0)
					fprintf(stderr, "invert_remove_inverse WARNING: `new_predicates` is empty.\n");
#endif

				for (hol_term* new_predicate : new_predicates) {
					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = second_head_operand->any_array.all;
					conjunction->any_array.all->reference_count++;
					conjunction->any_array.left.length = 0;
					conjunction->any_array.left.operands = nullptr;
					conjunction->any_array.right.length = 0;
					conjunction->any_array.right.operands = nullptr;
					conjunction->any_array.any.length = 0;
					conjunction->any_array.any.operands = nullptr;

					conjunction->any_array.left.length = second_head_operand->any_array.left.length;
					if (conjunction->any_array.left.length == 0) {
						conjunction->any_array.left.operands = nullptr;
					} else {
						conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.left.length);
						if (conjunction->any_array.left.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
							if (second_predicate_index.position == head_position::LEFT && second_predicate_index.index == i) {
								conjunction->any_array.left.operands[i] = new_predicate;
							} else {
								conjunction->any_array.left.operands[i] = second_head_operand->any_array.left.operands[i];
							}
							conjunction->any_array.left.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.right.length = second_head_operand->any_array.right.length;
					if (conjunction->any_array.right.length == 0) {
						conjunction->any_array.right.operands = nullptr;
					} else {
						conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.right.length);
						if (conjunction->any_array.right.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
							if (second_predicate_index.position == head_position::RIGHT && second_predicate_index.index == conjunction->any_array.right.length - i - 1) {
								conjunction->any_array.right.operands[i] = new_predicate;
							} else {
								conjunction->any_array.right.operands[i] = second_head_operand->any_array.right.operands[i];
							}
							conjunction->any_array.right.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.any.length = second_head_operand->any_array.any.length;
					if (conjunction->any_array.any.length == 0) {
						conjunction->any_array.any.operands = nullptr;
					} else {
						conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.any.length);
						if (conjunction->any_array.any.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
							if (second_predicate_index.position == head_position::ANY && second_predicate_index.index == i) {
								conjunction->any_array.any.operands[i] = new_predicate;
							} else {
								conjunction->any_array.any.operands[i] = second_head_operand->any_array.any.operands[i];
							}
							conjunction->any_array.any.operands[i]->reference_count++;
						}
					}

					hol_term* new_second_head = hol_term::new_exists(second_head->quantifier.variable, conjunction);
					if (new_second_head == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunction); free(conjunction);
						return false;
					}
					intersect<built_in_predicates>(dst, new_second_head, first_head);
					free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);

					if (!dst_outer.ensure_capacity(dst.length)) {
						free_all(dst); free_all(new_predicates);
						return false;
					}
					for (unsigned int i = 0; i < dst.length; i++)
						dst_outer[i] = &HOL_ZERO;
					dst_outer.length = dst.length;
					HOL_ZERO.reference_count += dst.length;
				}
				free_all(new_predicates);
				return (dst.length > 0);
			} else {
				fprintf(stderr, "invert_remove_inverse ERROR: Expected `second_head` to be an existentially quantified conjunction.\n");
				return false;
			}
		});
}

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant>
inline bool invert_select_arg_without_head_predicative(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	if (second->type != hol_term_type::LAMBDA)
		return false;
	unsigned int lambda_variable = second->quantifier.variable;
	second = second->quantifier.operand;

	auto on_remap_variables = [lambda_variable](hol_term* first_head, hol_term* second_head,
			apply_head_inverter& first_head_inverter, apply_head_inverter& second_head_inverter,
			array_map<const hol_term*, unsigned int>& second_variable_map, unsigned int& max_variable,
			bool first_head_is_array, bool second_head_is_array)
	{
		if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr && first_head->any.included->type == hol_term_type::EXISTS)
			first_head = first_head->any.included;
		if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
			second_head = second_head->any.included;

		while (first_head->type == hol_term_type::NOT)
			first_head = first_head->unary.operand;

		hol_term* expected_arg = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), &HOL_ANY), &HOL_ANY);
		if (expected_arg == nullptr)
			return false;
		HOL_ANY.reference_count += 2;

		hol_term* expected_head;
		if (ConjunctIndex >= 0) {
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(&HOL_ANY, ConjunctIndex), expected_arg), make_array_view((hol_term**) nullptr, 0)));
			if (expected_head == nullptr) {
				free(*expected_arg); free(expected_arg);
				return false;
			}
			HOL_ANY.reference_count += 1 + ConjunctIndex;
		} else {
			unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(expected_arg, make_repeated_array_view(&HOL_ANY, index))));
			if (expected_head == nullptr) {
				free(*expected_arg); free(expected_arg);
				return false;
			}
			HOL_ANY.reference_count += 1 + index;
		}

		array<hol_term*> intersection(2);
		unsigned int first_element_variable;
		intersect<built_in_predicates>(intersection, first_head, expected_head);
		bool element_wide_scope = (intersection.length != 0);
		if (intersection.length != 0) {
			free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);

			hol_term* operand;
			if (intersection[0]->type == hol_term_type::ANY_QUANTIFIER) {
				operand = intersection[0]->any_quantifier.operand;
			} else {
				operand = intersection[0]->quantifier.operand;
			}

			hol_term* conjunct;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (ConjunctIndex >= 0) {
					conjunct = operand->any_array.left.operands[ConjunctIndex];
				} else {
					conjunct = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex];
				}
			} else {
				unsigned int index;
				if (ConjunctIndex >= 0) {
					index = ConjunctIndex;
				} else {
					index = operand->array.length + ConjunctIndex;
				}
				conjunct = operand->array.operands[index];
			}

			if (conjunct->binary.right->type == hol_term_type::ANY || conjunct->binary.right->type == hol_term_type::ANY_RIGHT) {
				/* the "element variable" is not defined in `first_head` so we don't have to worry about variable agreement */
				first_element_variable = 0;
			} else if (conjunct->binary.right->type == hol_term_type::VARIABLE) {
				first_element_variable = conjunct->binary.right->variable;
			} else if (conjunct->binary.right->type == hol_term_type::CONSTANT) {
				/* `first` has no element variable, so we don't need to remap variables */
				first_element_variable = 0;
			} else {
				/* this should be a variable */
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }

		} else {
			expected_arg->reference_count++;
			free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);

			/* consider the case where the left conjunct in the head scope is an existential quantification */
			hol_term* expected_conjunct = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&expected_arg, 1)));
			if (expected_conjunct == nullptr) {
				free(*expected_arg); free(expected_arg);
				return false;
			}
			HOL_ANY.reference_count++;

			if (ConjunctIndex >= 0) {
				expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(&HOL_ANY, ConjunctIndex), expected_conjunct), make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				HOL_ANY.reference_count += 1 + ConjunctIndex;
			} else {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(expected_conjunct, make_repeated_array_view(&HOL_ANY, index))));
				if (expected_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				HOL_ANY.reference_count += 1 + index;
			}

			intersect<built_in_predicates>(intersection, first_head, expected_head);
			free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
			if (intersection.length == 0)
				return false;

			hol_term* operand;
			if (intersection[0]->type == hol_term_type::ANY_QUANTIFIER) {
				operand = intersection[0]->any_quantifier.operand;
			} else {
				operand = intersection[0]->quantifier.operand;
			}

			hol_term* conjunct;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (ConjunctIndex >= 0) {
					conjunct = operand->any_array.left.operands[ConjunctIndex];
				} else {
					conjunct = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex];
				}
			} else {
				unsigned int index;
				if (ConjunctIndex >= 0) {
					index = ConjunctIndex;
				} else {
					index = operand->array.length + ConjunctIndex;
				}
				conjunct = operand->array.operands[index];
			}

			hol_term* inner_operand;
			if (conjunct->type == hol_term_type::ANY_QUANTIFIER)
				inner_operand = conjunct->any_quantifier.operand;
			else inner_operand = conjunct->quantifier.operand;

			hol_term* inner_conjunct;
			if (inner_operand->type == hol_term_type::ANY_ARRAY)
				inner_conjunct = inner_operand->any_array.right.operands[inner_operand->any_array.right.length - 1];
			else inner_conjunct = inner_operand->array.operands[inner_operand->array.length - 1];

			if (inner_conjunct->binary.right->type == hol_term_type::ANY) {
				/* the "element variable" is not defined in `first_head` so we don't have to worry about variable agreement */
				first_element_variable = 0;
			} else if (inner_conjunct->binary.right->type == hol_term_type::VARIABLE) {
				first_element_variable = inner_conjunct->binary.right->variable;
			} else {
				/* this should be a variable */
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		}

		hol_term* right; hol_term* inner_right;
		hol_term* set_definition = nullptr;
		unsigned int second_element_variable;
		if (second_head->type == hol_term_type::EXISTS) {
			unsigned int set_variable = second_head->quantifier.variable;
			hol_term* operand = second_head->quantifier.operand;

			hol_term* left = nullptr;
			if (operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::ANY_ARRAY || operand->any_array.oper == hol_term_type::AND)) {
				left = (operand->any_array.left.length == 0 ? operand->any_array.all : operand->any_array.left.operands[0]);
				right = (operand->any_array.right.length == 0 ? operand->any_array.all : operand->any_array.right.operands[operand->any_array.right.length - 1]);
			} else if (operand->type == hol_term_type::AND) {
				left = operand->array.operands[0];
				right = operand->array.operands[operand->array.length - 1];
			} else if (operand->type == hol_term_type::ANY) {
				second_element_variable = 0;
				left = nullptr;
				right = nullptr;
			} else {
				return false;
			}

			if (left != nullptr) {
				if (left->type == hol_term_type::EQUALS) {
					set_definition = left->binary.right;
				} else if (left->type == hol_term_type::BINARY_APPLICATION) {
					set_definition = left->ternary.third;
				} else if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
					set_definition = left->any.included;
				} else {
					set_definition = nullptr;
				}
			}

			if (right != nullptr) {
				while (right->type == hol_term_type::NOT)
					right = right->unary.operand;
				if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
					right = right->any.included;
				if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
					right = right->binary.right;
				if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
					right = right->any.included;

				if (right->type == hol_term_type::FOR_ALL) {
					second_element_variable = right->quantifier.variable;
					inner_right = right->quantifier.operand->binary.right;
				} else if (right->type == hol_term_type::EXISTS) {
					second_element_variable = right->quantifier.variable;
					inner_right = right->quantifier.operand->array.operands[right->quantifier.operand->array.length - 1];
				} else if (right->type == hol_term_type::AND && right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
						&& right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[0]->binary.left->variable == set_variable
						&& right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
				{
					second_element_variable = right->array.operands[0]->binary.right->variable;
					inner_right = right->array.operands[right->array.length - 1];
				} else if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::VARIABLE
						&& right->binary.left->variable == lambda_variable && right->binary.right->type == hol_term_type::VARIABLE)
				{
					second_element_variable = right->binary.right->variable;
					inner_right = nullptr;
				} else {
					return false;
				}
			}
		} else if (second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT
				|| (second_head->type == hol_term_type::ANY_QUANTIFIER && has_intersection(hol_quantifier_type::EXISTS, second_head->any_quantifier.quantifier))) {
			second_element_variable = 0;
		} else {
			return false;
		}

		if (first_element_variable != 0 && second_element_variable != 0) {
			if (!element_wide_scope) {
				for (unsigned int i = 0; i < second_variable_map.size; i++) {
					if (second_variable_map.keys[i]->quantifier.variable == second_element_variable) {
						second_variable_map.values[i] = first_element_variable;
						if (first_element_variable == second_element_variable)
							second_variable_map.remove_at(i);
						break;
					}
				}
				return true;
			}

			if (right == nullptr || inner_right == nullptr)
				return false;

			array_map<unsigned int, variable_set> free_variable_map(4);
			free_variable_map.keys[0] = second_element_variable;
			free_variable_map.values[0].type = variable_set_type::SINGLETON;
			free_variable_map.values[0].variable = first_element_variable;
			free_variable_map.size = 1;

			if (first_head_inverter.outer.last() != first_head && !first_head_inverter.outer.add(first_head))
				return false;
			if (second_head_inverter.outer.last() != second_head && !second_head_inverter.outer.add(second_head))
				return false;

			if (!second_head_inverter.outer.add(right)
			 || !second_head_inverter.outer.add(inner_right))
				return false;

			if (remap_scopes(free_variable_map, first_head_inverter, second_head_inverter, second_variable_map, max_variable, first_head_is_array, second_head_is_array)) {
				/* check to see what `right` was mapped to, and make sure `set_definition` is mapped to the same variable */
				if (set_definition != nullptr && set_definition->type == hol_term_type::LAMBDA) {
					bool contains;
					unsigned int new_var = second_variable_map.get(right, contains);
					if (!contains) new_var = right->quantifier.variable;
					unsigned int index = second_variable_map.index_of(set_definition);
					if (index < second_variable_map.size) {
						second_variable_map.values[index] = new_var;
						if (set_definition->quantifier.variable == new_var)
							second_variable_map.remove_at(index);
					} else if (set_definition->quantifier.variable != new_var) {
						if (!second_variable_map.put(set_definition, new_var)) return false;
					}
				}
			} else {
				/* try removing the set */
				for (auto entry : second_variable_map) free(entry.value);

				for (unsigned int i = 0; i < second_variable_map.size; i++) {
					if (second_variable_map.keys[i]->quantifier.variable == second_element_variable) {
						second_variable_map.values[i] = first_element_variable;
						if (first_element_variable == second_element_variable)
							second_variable_map.remove_at(i);
						break;
					}
				}

				hol_term* second_set_definition;
				if (set_definition == nullptr || set_definition->type != hol_term_type::LAMBDA) {
					return true;
				} else {
					second_set_definition = set_definition->quantifier.operand;
				}

				/* find the relativization of the set in `first` */
				hol_term* first_scope = nullptr;
				for (unsigned int i = first_head_inverter.outer.length - 1; i > 0; i--) {
					first_scope = first_head_inverter.outer[i - 1];
					if ((first_scope->type == hol_term_type::FOR_ALL || first_scope->type == hol_term_type::EXISTS || first_scope->type == hol_term_type::LAMBDA) && first_scope->quantifier.variable == first_element_variable)
						break;
				}

				hol_term* first_set_definition;
				if (first_scope->type == hol_term_type::FOR_ALL && first_scope->quantifier.operand->type == hol_term_type::IF_THEN) {
					first_set_definition = first_scope->quantifier.operand->binary.left;
					first_set_definition->reference_count++;
				} else if (first_scope->type == hol_term_type::EXISTS && first_scope->quantifier.operand->type == hol_term_type::AND) {
					if (first_scope->quantifier.operand->array.length > 2) {
						first_set_definition = hol_term::new_and(make_array_view(first_scope->quantifier.operand->array.operands, first_scope->quantifier.operand->array.length - 1));
						if (first_set_definition == nullptr) return false;
						for (unsigned int i = 0; i < first_set_definition->array.length; i++)
							first_set_definition->array.operands[i]->reference_count++;
					} else {
						first_set_definition = first_scope->quantifier.operand->array.operands[0];
						first_set_definition->reference_count++;
					}
				} else {
					return true;
				}

				array<pair<hol_term*, variable_map>> intersection(2);
				intersect<built_in_predicates, true, true>(intersection, first_set_definition, second_set_definition);
				if (intersection.length == 0) {
					free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
					return false;
				} else if (intersection.length != 1) {
					fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: Intersection is not unique.\n");
					free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
					free_all(intersection); return false;
				}

				const variable_map& var_map = intersection[0].value;
				for (const auto& entry : var_map.scope_map) {
					const variable_set& set = entry.value;
					unsigned int src_variable = entry.key.src->quantifier.variable;
					unsigned int target_var;
					if (!get_target_variable(src_variable, target_var, set, max_variable)
					 || !second_variable_map.ensure_capacity(second_variable_map.size + 1))
					{
						free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
						free_all(intersection); return false;
					}

					unsigned int index = second_variable_map.index_of(entry.key.src);
					if (index < second_variable_map.size) {
						second_variable_map.values[index] = target_var;
						if (src_variable == target_var)
							second_variable_map.remove_at(index);
					} else if (src_variable != target_var) {
						second_variable_map.keys[second_variable_map.size] = entry.key.src;
						second_variable_map.values[second_variable_map.size] = target_var;
						second_variable_map.size++;
					}
				} for (const auto& entry : var_map.free_variables) {
					const variable_set& set = entry.value;
					unsigned int src_variable = entry.key;
					if (src_variable == second_element_variable)
						continue;
					if (!free_variable_map.ensure_capacity(free_variable_map.size + 1)) {
						free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
						free_all(intersection); return false;
					}
					unsigned int index = free_variable_map.index_of(src_variable);
					if (index < free_variable_map.size) {
						variable_set& new_set = *((variable_set*) alloca(sizeof(variable_map)));
						if (!intersect(new_set, set, free_variable_map.values[index])) {
							free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
							free_all(intersection); return false;
						}
						swap(free_variable_map.values[index], new_set);
						free(new_set);
					} else {
						free_variable_map.keys[free_variable_map.size] = src_variable;
						if (!init(free_variable_map.values[free_variable_map.size], set)) {
							free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
							free_all(intersection); return false;
						}
						free_variable_map.size++;
					}
				}
				free(*first_set_definition); if (first_set_definition->reference_count == 0) free(first_set_definition);
				free_all(intersection);

				if (!remap_scopes(free_variable_map, first_head_inverter, second_head_inverter, second_variable_map, max_variable, first_head_is_array, second_head_is_array)) {
					for (auto entry : second_variable_map) free(entry.value);
					return false;
				}
			}
		}
		return true;
	};

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates, false, true>, predicative_head_finder<built_in_predicates>(lambda_variable), on_remap_variables,
		[lambda_variable](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_predicate_index.position != head_position::NONE)
				first_head = first_head->any.included;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
				second_head = second_head->any.included;

			while (first_head->type == hol_term_type::NOT)
				first_head = first_head->unary.operand;

			unsigned int predicate_variable;
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				predicate_variable = max_variable + 1;
			} else {
#if !defined(NDEBUG)
				if (first_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_select_arg_without_head_predicative WARNING: Expected `first_head` to be an existential quantification.\n");
#endif
				predicate_variable = first_head->quantifier.variable;
			}

			if (second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) {
				if (!dst.ensure_capacity(dst.length + 2) || !dst_outer.ensure_capacity(dst_outer.length + 2))
					return false;
				hol_term* any_head_quantifier = hol_term::new_any(hol_term::new_exists(predicate_variable, &HOL_ANY));
				if (any_head_quantifier == nullptr)
					return false;
				hol_term* arg = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(predicate_variable)),
						hol_term::new_any(nullptr, &any_head_quantifier, 1));
				if (arg == nullptr) {
					free(*any_head_quantifier); free(any_head_quantifier);
					return false;
				}

				if (ConjunctIndex >= 0) {
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(&HOL_ANY, ConjunctIndex), arg), make_array_view((hol_term**) nullptr, 0)));
					if (dst[dst.length] == nullptr)
						return false;
					HOL_ANY.reference_count += 1 + ConjunctIndex;
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg, make_repeated_array_view(&HOL_ANY, index))));
					if (dst[dst.length] == nullptr)
						return false;
					HOL_ANY.reference_count += 1 + index;
				} else {
					fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: Unsupported value of `ConjunctIndex`.\n");
					return false;
				}
				dst.length++;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;

				unsigned int conjunct_variable = max_variable + 2;
				hol_term* conjunct = hol_term::new_exists(conjunct_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&arg, 1)));
				if (conjunct == nullptr)
					return false;
				arg->reference_count++;
				HOL_ANY.reference_count++;
				if (ConjunctIndex >= 0) {
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(&HOL_ANY, ConjunctIndex), conjunct), make_array_view((hol_term**) nullptr, 0)));
					HOL_ANY.reference_count += 1 + ConjunctIndex;
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(conjunct, make_repeated_array_view(&HOL_ANY, index))));
					HOL_ANY.reference_count += 1 + index;
				} else {
					fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: Unsupported value of `ConjunctIndex`.\n");
					return false;
				}
				if (dst[dst.length] == nullptr) {
					free(*conjunct); free(conjunct);
					return false;
				}
				dst.length++;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
				return true;
			} else if (second_head->type == hol_term_type::EXISTS) {
				unsigned int set_variable = second_head->quantifier.variable;
				hol_term* operand = second_head->quantifier.operand;

				hol_term* left = nullptr;
				hol_term* right = nullptr;
				hol_term* set_definition = nullptr;
				bool can_keep_set_variable = true;
				bool can_remove_set_variable = true;
				bool must_be_simple_set_def = false;
				if (operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::ANY_ARRAY || operand->any_array.oper == hol_term_type::AND)) {
					left = (operand->any_array.left.length == 0 ? operand->any_array.all : operand->any_array.left.operands[0]);
					right = (operand->any_array.right.length == 0 ? operand->any_array.all : operand->any_array.right.operands[operand->any_array.right.length - 1]);
					if (left->type == hol_term_type::EQUALS && left->binary.left->type == hol_term_type::VARIABLE && left->binary.left->variable == set_variable && left->binary.right->type == hol_term_type::LAMBDA) {
						must_be_simple_set_def = true;
						set_definition = left->binary.right;
					} else if (left->type == hol_term_type::BINARY_APPLICATION && left->ternary.second->type == hol_term_type::VARIABLE && left->ternary.second->variable == set_variable && left->ternary.third->type == hol_term_type::LAMBDA) {
						must_be_simple_set_def = false;
						can_remove_set_variable = false;
						set_definition = left->ternary.third;
					}
				} else if (operand->type == hol_term_type::AND) {
					left = operand->array.operands[0];
					right = operand->array.operands[operand->array.length - 1];
					if (left->type == hol_term_type::EQUALS && left->binary.left->type == hol_term_type::VARIABLE && left->binary.left->variable == set_variable && left->binary.right->type == hol_term_type::LAMBDA) {
						must_be_simple_set_def = true;
						set_definition = left->binary.right;
					} else if (left->type == hol_term_type::BINARY_APPLICATION && left->ternary.second->type == hol_term_type::VARIABLE && left->ternary.second->variable == set_variable && left->ternary.third->type == hol_term_type::LAMBDA) {
						must_be_simple_set_def = false;
						can_remove_set_variable = false;
						set_definition = left->ternary.third;
					}
					if (operand->array.length == 2 && must_be_simple_set_def)
						can_keep_set_variable = false;
					if (operand->array.length > 2)
						can_remove_set_variable = false;
				}

				if (right->type == hol_term_type::ANY) {
					fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: `ANY` type is unsupported for the right conjunct.");
					return false;
				}

				if (!dst.ensure_capacity(dst.length + 2) || !dst_outer.ensure_capacity(dst_outer.length + 2))
					return false;

				hol_term* head_var = hol_term::new_variable(predicate_variable);
				if (head_var == nullptr) return false;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_exists(predicate_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return false;
				}
				free(*head_var);

				hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
				if (conjunct == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return false;
				}
				HOL_ANY.reference_count++;

				unsigned int negation_count = 0;
				while (right->type == hol_term_type::NOT) {
					right = right->unary.operand;
					negation_count++;
				}

				hol_term* old_right = right;
				if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
					right = right->any.included;

				hol_term* element_var;
				unsigned int element_variable;
				bool wide_scope_marker = false;
				bool element_narrow_scope = false;
				if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					right = right->binary.right;
					element_narrow_scope = true;
					wide_scope_marker = true;
				}
				if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
					right = right->any.included;
				if (right->type == hol_term_type::FOR_ALL) {
					if (right->quantifier.operand->type == hol_term_type::IF_THEN && right->quantifier.operand->binary.right->type == hol_term_type::UNARY_APPLICATION
					 && right->quantifier.operand->binary.right->binary.left->type == hol_term_type::CONSTANT && right->quantifier.operand->binary.right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
						wide_scope_marker = true;
					element_variable = right->quantifier.variable;
					element_var = hol_term::new_variable(element_variable);
					if (element_var == nullptr) {
						free(*conjunct); free(conjunct);
						return false;
					}
				} else if (right->type == hol_term_type::EXISTS) {
					if (right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.operands[right->quantifier.operand->array.length - 1]->type == hol_term_type::UNARY_APPLICATION
					 && right->quantifier.operand->array.operands[right->quantifier.operand->array.length - 1]->binary.left->type == hol_term_type::CONSTANT
					 && right->quantifier.operand->array.operands[right->quantifier.operand->array.length - 1]->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
						wide_scope_marker = true;
					if (right->quantifier.operand->type == hol_term_type::ANY_ARRAY && (right->quantifier.operand->any_array.oper == hol_term_type::ANY_ARRAY || right->quantifier.operand->any_array.oper == hol_term_type::AND)
					 && right->quantifier.operand->any_array.right.length > 0 && right->quantifier.operand->any_array.right.operands[right->quantifier.operand->any_array.right.length - 1]->type == hol_term_type::UNARY_APPLICATION
					 && right->quantifier.operand->any_array.right.operands[right->quantifier.operand->any_array.right.length - 1]->binary.left->type == hol_term_type::CONSTANT
					 && right->quantifier.operand->any_array.right.operands[right->quantifier.operand->any_array.right.length - 1]->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
						wide_scope_marker = true;
					element_variable = right->quantifier.variable;
					element_var = hol_term::new_variable(element_variable);
					if (element_var == nullptr) {
						free(*conjunct); free(conjunct);
						return false;
					}
				} else if (right->type == hol_term_type::AND && right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
						&& right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[0]->binary.left->variable == set_variable
						&& right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
				{
					if (right->array.operands[right->array.length - 1]->type == hol_term_type::UNARY_APPLICATION && right->array.operands[right->array.length - 1]->binary.left->type == hol_term_type::CONSTANT && right->array.operands[right->array.length - 1]->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
						wide_scope_marker = true;
					element_var = right->array.operands[0]->binary.right;
					element_var->reference_count++;
					element_variable = element_var->variable;
				} else if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::VARIABLE
						&& right->binary.left->variable == lambda_variable && right->binary.right->type == hol_term_type::VARIABLE)
				{
					element_var = right->binary.right;
					element_var->reference_count++;
					element_variable = element_var->variable;
				} else {
					free(*conjunct); free(conjunct);
					return false;
				}

				if (must_be_simple_set_def && can_remove_set_variable
				 && set_definition->quantifier.operand->type == hol_term_type::EQUALS
				 && set_definition->quantifier.operand->binary.left->type == hol_term_type::VARIABLE
				 && set_definition->quantifier.operand->binary.left->variable == set_definition->quantifier.variable
				 && (set_definition->quantifier.operand->binary.right->type == hol_term_type::CONSTANT
				  || set_definition->quantifier.operand->binary.right->type == hol_term_type::ANY_CONSTANT
				  || set_definition->quantifier.operand->binary.right->type == hol_term_type::ANY_CONSTANT_EXCEPT
				  || set_definition->quantifier.operand->binary.right->type == hol_term_type::ANY)
				 && right->type != hol_term_type::FOR_ALL && right->type != hol_term_type::AND)
				{
					/* the set contains a single constant */
					hol_term* constant = set_definition->quantifier.operand->binary.right;
					hol_term* arg_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), head_var), constant);
					if (arg_term == nullptr) {
						free(*conjunct); free(conjunct);
						free(*element_var); if (element_var->reference_count == 0) free(element_var);
						return false;
					}
					head_var->reference_count++;
					constant->reference_count++;

					array<hol_term*> arg_terms(2);
					intersect<built_in_predicates>(arg_terms, arg_term, conjunct);
					free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
					if (!dst.ensure_capacity(dst.length + arg_terms.length) || !dst_outer.ensure_capacity(dst_outer.length + arg_terms.length)) {
						for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*element_var); if (element_var->reference_count == 0) free(element_var);
						return false;
					}
					for (hol_term* arg_term : arg_terms) {
						if (ConjunctIndex >= 0) {
							dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
									make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0)));
							conjunct->reference_count += 1 + ConjunctIndex;
						} else if (ConjunctIndex < 0) {
							unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
							dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
									make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index))));
							conjunct->reference_count += 1 + index;
						} else {
							fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: Unsupported value of `ConjunctIndex`.\n");
							for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*element_var); if (element_var->reference_count == 0) free(element_var);
							return false;
						}
						if (dst[dst.length] == nullptr) {
							for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*element_var); if (element_var->reference_count == 0) free(element_var);
							return false;
						}
						arg_term->reference_count++;
						dst.length++;

						for (unsigned int i = 0; i < negation_count; i++) {
							hol_term* temp = hol_term::new_not(dst.last());
							if (temp == nullptr) {
								for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*element_var); if (element_var->reference_count == 0) free(element_var);
								return false;
							}
							dst.last() = temp;
						}

						dst_outer[dst_outer.length++] = &HOL_ZERO;
						HOL_ZERO.reference_count++;
					}
					for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
					if (set_definition->quantifier.operand->binary.right->type != hol_term_type::ANY) {
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*element_var); if (element_var->reference_count == 0) free(element_var);
						return true;
					}
				}

				hol_term* set_var = hol_term::new_variable(set_variable);
				if (set_var == nullptr) {
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
					free(*element_var); if (element_var->reference_count == 0) free(element_var);
					return false;
				}

				hol_term* arg_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), head_var), element_var);
				if (arg_term == nullptr) {
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
					free(*element_var); if (element_var->reference_count == 0) free(element_var);
					free(*set_var); free(set_var);
					return false;
				}
				head_var->reference_count++;

				if (can_remove_set_variable) {
					hol_term* expected_left = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, &HOL_ANY));
					if (expected_left == nullptr) {
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); free(arg_term);
						free(*set_var); free(set_var);
						return false;
					}
					HOL_ANY.reference_count++;
					set_var->reference_count++;

					array<hol_term*> left_intersections(4);
					intersect<built_in_predicates>(left_intersections, expected_left, left);
					free(*expected_left); free(expected_left);

					expected_left = hol_term::new_apply(&HOL_ANY, set_var, hol_term::new_lambda(element_variable, &HOL_ANY));
					if (expected_left == nullptr) {
						free_all(left_intersections);
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); free(arg_term);
						free(*set_var); free(set_var);
						return false;
					}
					HOL_ANY.reference_count += 2;
					set_var->reference_count++;

					intersect<built_in_predicates>(left_intersections, expected_left, left);
					free(*expected_left); free(expected_left);

					if (!dst.ensure_capacity(dst.length + left_intersections.length + 1) || !dst_outer.ensure_capacity(dst_outer.length + left_intersections.length + 1)) {
						for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						return false;
					}
					for (hol_term* left_intersection : left_intersections) {
						hol_term* set_definition;
						if (left_intersection->type == hol_term_type::EQUALS)
							set_definition = left_intersection->binary.right->quantifier.operand;
						else set_definition = left_intersection->ternary.third->quantifier.operand;

						if (right->type == hol_term_type::FOR_ALL || right->type == hol_term_type::EXISTS) {
							/* this case was originally only for FOR_ALL, since we always moved universal quantifiers
							   outside the head scope, and we would keep existential quantifiers inside the head scope
							   if possible; but now we move both quantifiers outside the head scope */
							hol_term* excluded_quantifier;
							if (right->type == hol_term_type::FOR_ALL)
								excluded_quantifier = hol_term::new_any(hol_term::new_for_all(element_variable, &HOL_ANY));
							else excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							HOL_ANY.reference_count++;

							if (ConjunctIndex >= 0) {
								dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0)));
								if (dst[dst.length] == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									free(*excluded_quantifier); free(excluded_quantifier);
									return false;
								}
								arg_term->reference_count++;
								conjunct->reference_count += 1 + ConjunctIndex;
							} else if (ConjunctIndex < 0) {
								unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
								dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index))));
								if (dst[dst.length] == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									free(*excluded_quantifier); free(excluded_quantifier);
									return false;
								}
								arg_term->reference_count++;
								conjunct->reference_count += 1 + index;
							}
							dst.length++;

							if (is_array) {
								hol_term* quantified_term;
								if (right->type == hol_term_type::FOR_ALL) {
									quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(set_definition, hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1)));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									set_definition->reference_count++;
									HOL_ZERO.reference_count++;
								} else if (set_definition->type == hol_term_type::AND) {
									quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(
											make_array_view(set_definition->array.operands, set_definition->array.length),
											hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1))));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									for (unsigned int i = 0; i < set_definition->array.length; i++)
										set_definition->array.operands[i]->reference_count++;
									HOL_ZERO.reference_count++;
								} else {
									quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1)));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									set_definition->reference_count++;
									HOL_ZERO.reference_count++;
								}

								for (unsigned int i = 0; i < negation_count; i++) {
									hol_term* temp = hol_term::new_not(quantified_term);
									if (temp == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*quantified_term); free(quantified_term);
										free_all(dst); free_all(dst_outer); return false;
									}
									quantified_term = temp;
								}

								if (right->type == hol_term_type::FOR_ALL) {
									dst_outer[dst_outer.length] = hol_term::new_apply(
											hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), quantified_term);
									if (dst_outer[dst_outer.length] == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*quantified_term); free(quantified_term);
										free_all(dst); free_all(dst_outer); return false;
									}
								} else {
									dst_outer[dst_outer.length] = quantified_term;
								}
								dst_outer.length++;
							} else {
								hol_term* quantified_term;
								if (right->type == hol_term_type::FOR_ALL) {
									quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(set_definition, hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1)));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									set_definition->reference_count++;
									HOL_ZERO.reference_count++;
								} else if (set_definition->type == hol_term_type::AND) {
									quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(
											make_array_view(set_definition->array.operands, set_definition->array.length),
											hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1))));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									for (unsigned int i = 0; i < set_definition->array.length; i++)
										set_definition->array.operands[i]->reference_count++;
									HOL_ZERO.reference_count++;
								} else {
									quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, hol_term::new_any_right(&HOL_ZERO, &excluded_quantifier, 1)));
									if (quantified_term == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_quantifier); free(excluded_quantifier);
										free_all(dst); free_all(dst_outer); return false;
									}
									set_definition->reference_count++;
									HOL_ZERO.reference_count++;
								}

								for (unsigned int i = 0; i < negation_count; i++) {
									hol_term* temp = hol_term::new_not(quantified_term);
									if (temp == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*quantified_term); free(quantified_term);
										free_all(dst); free_all(dst_outer); return false;
									}
									quantified_term = temp;
								}

								if (right->type == hol_term_type::FOR_ALL) {
									hol_term* excluded_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
									if (excluded_wide_scope == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*quantified_term); free(quantified_term);
										free_all(dst); free_all(dst_outer); return false;
									}
									HOL_ANY.reference_count++;

									dst_outer[dst_outer.length] = hol_term::new_any_right_only(hol_term::new_apply(
											hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE),
											hol_term::new_any_right_only(quantified_term, &excluded_wide_scope, 1)));
									if (dst_outer[dst_outer.length] == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										free(*excluded_wide_scope); free(excluded_wide_scope);
										free(*quantified_term); free(quantified_term);
										free_all(dst); free_all(dst_outer); return false;
									}
								} else {
									dst_outer[dst_outer.length] = quantified_term;
								}
								dst_outer.length++;
							}
							could_have_wide_scope = true;
						} else {
							for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							return false;
						}
						// else if (wide_scope_marker && !element_narrow_scope) {
						// 	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
						// 	if (excluded_quantifier == nullptr) {
						// 		for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 		free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 		free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 		free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 		return false;
						// 	}
						// 	HOL_ANY.reference_count++;

						// 	hol_term* narrow_scope;
						// 	if (ConjunctIndex >= 0) {
						// 		narrow_scope = hol_term::new_any_right_only(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right_only(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
						// 				make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0))))), &excluded_quantifier, 1);
						// 		if (narrow_scope == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			free(*excluded_quantifier); free(excluded_quantifier);
						// 			return false;
						// 		}
						// 		arg_term->reference_count++;
						// 		conjunct->reference_count += 1 + ConjunctIndex;
						// 	} else if (ConjunctIndex < 0) {
						// 		unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						// 		narrow_scope = hol_term::new_any_right_only(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right_only(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
						// 				make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index)))))), &excluded_quantifier, 1);
						// 		if (narrow_scope == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			free(*excluded_quantifier); free(excluded_quantifier);
						// 			return false;
						// 		}
						// 		arg_term->reference_count++;
						// 		conjunct->reference_count += 1 + index;
						// 	}

						// 	if (set_definition->type == hol_term_type::AND) {
						// 		dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(make_array_view(set_definition->array.operands, set_definition->array.length), narrow_scope)));
						// 		if (dst[dst.length] == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*narrow_scope); free(narrow_scope);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			return false;
						// 		}
						// 		for (unsigned int i = 0; i < set_definition->array.length; i++)
						// 			set_definition->array.operands[i]->reference_count++;
						// 	} else if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::AND) {
						// 		dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_any_array(hol_term_type::AND, set_definition->any_array.all,
						// 				make_array_view(set_definition->any_array.any.operands, set_definition->any_array.any.length),
						// 				make_array_view(set_definition->any_array.left.operands, set_definition->any_array.left.length),
						// 				make_appended_array_view(make_array_view(set_definition->any_array.right.operands, set_definition->any_array.right.length), narrow_scope)));
						// 		if (dst[dst.length] == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*narrow_scope); free(narrow_scope);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			return false;
						// 		}
						// 		set_definition->any_array.all->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.left.length; i++)
						// 			set_definition->any_array.left.operands[i]->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.right.length; i++)
						// 			set_definition->any_array.right.operands[i]->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.any.length; i++)
						// 			set_definition->any_array.any.operands[i]->reference_count++;
						// 	} else {
						// 		if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::ANY_ARRAY)
						// 			fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: The set definition has `ANY_ARRAY` type `oper` is not `AND`.\n");
						// 		dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, narrow_scope));
						// 		if (dst[dst.length] == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*narrow_scope); free(narrow_scope);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			return false;
						// 		}
						// 		set_definition->reference_count++;
						// 	}
						// 	dst.length++;

						// 	for (unsigned int i = 0; i < negation_count; i++) {
						// 		hol_term* temp = hol_term::new_not(dst.last());
						// 		if (temp == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			free_all(dst); free_all(dst_outer); return false;
						// 		}
						// 		dst.last() = temp;
						// 	}

						// 	dst_outer[dst_outer.length++] = &HOL_ZERO;
						// 	HOL_ZERO.reference_count++;
						// 	could_have_wide_scope = true;
						// } else {
						// 	hol_term* new_arg_term = hol_term::new_equals(
						// 			hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable_preimage(predicate_variable)),
						// 			hol_term::new_variable_preimage(element_variable));
						// 	if (new_arg_term == nullptr) {
						// 		for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 		free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 		free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 		free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 		return false;
						// 	}

						// 	hol_term* new_definition;
						// 	if (set_definition->type == hol_term_type::AND) {
						// 		new_definition = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(make_array_view(set_definition->array.operands, set_definition->array.length), new_arg_term)));
						// 		if (new_definition == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			free(*new_arg_term); free(new_arg_term);
						// 			return false;
						// 		}
						// 		for (unsigned int i = 0; i < set_definition->array.length; i++)
						// 			set_definition->array.operands[i]->reference_count++;
						// 	} else if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::AND) {
						// 		new_definition = hol_term::new_exists(element_variable, hol_term::new_any_array(hol_term_type::AND, set_definition->any_array.all,
						// 				make_array_view(set_definition->any_array.any.operands, set_definition->any_array.any.length),
						// 				make_array_view(set_definition->any_array.left.operands, set_definition->any_array.left.length),
						// 				make_appended_array_view(make_array_view(set_definition->any_array.right.operands, set_definition->any_array.right.length), new_arg_term)));
						// 		if (new_definition == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			free(*new_arg_term); free(new_arg_term);
						// 			return false;
						// 		}
						// 		set_definition->any_array.all->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.left.length; i++)
						// 			set_definition->any_array.left.operands[i]->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.right.length; i++)
						// 			set_definition->any_array.right.operands[i]->reference_count++;
						// 		for (unsigned int i = 0; i < set_definition->any_array.any.length; i++)
						// 			set_definition->any_array.any.operands[i]->reference_count++;
						// 	} else {
						// 		if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::ANY_ARRAY)
						// 			fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: The set definition has `ANY_ARRAY` type `oper` is not `AND`.\n");
						// 		new_definition = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, new_arg_term));
						// 		if (new_definition == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			return false;
						// 		}
						// 		set_definition->reference_count++;
						// 	}

						// 	array<pair<hol_term*, variable_map>> new_definitions(2);
						// 	intersect<built_in_predicates, true>(new_definitions, new_definition, conjunct);
						// 	free(*new_definition); if (new_definition->reference_count == 0) free(new_definition);
						// 	if (!dst.ensure_capacity(dst.length + new_definitions.length) || !dst_outer.ensure_capacity(dst_outer.length + new_definitions.length)) {
						// 		for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 		free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 		free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 		free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 		return false;
						// 	}
						// 	for (pair<hol_term*, variable_map>& pair : new_definitions) {
						// 		hol_term* new_definition = pair.key;

						// 		/* check if we need to relabel variables */
						// 		array_map<const hol_term*, unsigned int> variable_map(max((size_t) 1, pair.value.scope_map.size));
						// 		for (const auto& entry : pair.value.scope_map) {
						// 			const variable_set& set = entry.value;
						// 			unsigned int src_variable = entry.key.dst->quantifier.variable;
						// 			unsigned int target_var;
						// 			if (!get_target_variable(src_variable, target_var, set, max_variable)) {
						// 				for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 				free_all(new_definitions);
						// 				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 				free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 				return false;
						// 			}

						// 			if (target_var != src_variable)
						// 				variable_map.put(entry.key.dst, target_var);
						// 		}

						// 		hol_term* temp = map_variables<true>(new_definition, variable_map);
						// 		if (temp == nullptr) {
						// 			for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 			free_all(new_definitions);
						// 			free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 			free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 			free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 			return false;
						// 		}
						// 		free(*new_definition); if (new_definition->reference_count == 0) free(new_definition);
						// 		pair.key = temp;
						// 		new_definition = temp;

						// 		if (ConjunctIndex >= 0) {
						// 			dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
						// 					make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), new_definition), make_array_view((hol_term**) nullptr, 0)));
						// 			if (dst[dst.length] == nullptr) {
						// 				for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 				free_all(new_definitions);
						// 				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 				free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 				return false;
						// 			}
						// 			conjunct->reference_count += 1 + ConjunctIndex;
						// 		} else if (ConjunctIndex < 0) {
						// 			unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						// 			dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
						// 					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(new_definition, make_repeated_array_view(conjunct, index))));
						// 			if (dst[dst.length] == nullptr) {
						// 				for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 				free_all(new_definitions);
						// 				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 				free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 				return false;
						// 			}
						// 			conjunct->reference_count += 1 + index;
						// 		}
						// 		new_definition->reference_count++;
						// 		dst.length++;

						// 		for (unsigned int i = 0; i < negation_count; i++) {
						// 			hol_term* temp = hol_term::new_not(dst.last());
						// 			if (temp == nullptr) {
						// 				for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						// 				free_all(new_definitions);
						// 				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						// 				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						// 				free(*set_var); if (set_var->reference_count == 0) free(set_var);
						// 				return false;
						// 			}
						// 			dst.last() = temp;
						// 		}

						// 		dst_outer[dst_outer.length++] = &HOL_ZERO;
						// 		HOL_ZERO.reference_count++;
						// 		if (ConjunctIndex != -1)
						// 			could_have_wide_scope = true;
						// 	}
						// 	free_all(new_definitions);
						// }
					}
					for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
				}

				if (can_keep_set_variable && (right->type == hol_term_type::EXISTS || right->type == hol_term_type::FOR_ALL)) {
					hol_term* new_second_head;
					if (ConjunctIndex >= 0) {
						new_second_head = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
								make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0)));
						if (new_second_head == nullptr) {
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							return false;
						}
						arg_term->reference_count++;
						conjunct->reference_count += 1 + ConjunctIndex;
					} else if (ConjunctIndex < 0) {
						unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						new_second_head = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
								make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index))));
						if (new_second_head == nullptr) {
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							return false;
						}
						arg_term->reference_count++;
						conjunct->reference_count += 1 + index;
					} else {
						fprintf(stderr, "invert_select_arg_without_head_predicative ERROR: Unsupported value of `ConjunctIndex`.\n");
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						return false;
					}

					array<hol_term*> new_heads(2);
					intersect<built_in_predicates>(new_heads, new_second_head, first_head);
					free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);
					if (new_heads.length == 0 || !dst.ensure_capacity(dst.length + new_heads.length) || !dst_outer.ensure_capacity(dst_outer.length + new_heads.length)) {
						free_all(new_heads);
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						return false;
					}
					for (unsigned int i = 0; i < new_heads.length; i++) {
						hol_term* new_head = new_heads[i];
						if (wide_scope_marker && element_narrow_scope && right->type != hol_term_type::FOR_ALL) {
							hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right_only(new_head));
							if (temp == nullptr) {
								free_all(new_heads);
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							new_head = temp;
							new_heads[i] = temp;
						}
						dst[dst.length++] = new_head;
						new_head->reference_count++;

						hol_term* new_right_conjunct;
						hol_term* excluded_quantifier;
						if (right->type == hol_term_type::EXISTS) {
							excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								free_all(new_heads);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								return false;
							}
							HOL_ANY.reference_count++;
							new_right_conjunct = hol_term::new_any_right_only(hol_term::new_exists(element_variable, hol_term::new_and(hol_term::new_apply(set_var, element_var), hol_term::new_any_right_only(&HOL_ZERO, &excluded_quantifier, 1))));
						} else {
							excluded_quantifier = hol_term::new_any(hol_term::new_for_all(element_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								free_all(new_heads);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								return false;
							}
							HOL_ANY.reference_count++;
							if (wide_scope_marker) {
								new_right_conjunct = hol_term::new_any_right_only(hol_term::new_for_all(element_variable,
										hol_term::new_if_then(hol_term::new_apply(set_var, element_var), hol_term::new_any_right_only(&HOL_ZERO, &excluded_quantifier, 1))));
							} else {
								new_right_conjunct = hol_term::new_any_right_only(hol_term::new_apply(
										hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE),
										hol_term::new_any_right_only(hol_term::new_for_all(element_variable,
											hol_term::new_if_then(hol_term::new_apply(set_var, element_var), hol_term::new_any_right_only(&HOL_ZERO, &excluded_quantifier, 1))))));
							}
						}
						if (new_right_conjunct == nullptr) {
							free_all(new_heads);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*excluded_quantifier); free(excluded_quantifier);
							return false;
						}
						HOL_ZERO.reference_count++;
						element_var->reference_count++;
						set_var->reference_count++;

						if ((old_right->type == hol_term_type::ANY || old_right->type == hol_term_type::ANY_RIGHT) && !can_have_free_variables(*new_head))
							dst_outer[dst_outer.length] = substitute_head<any_node_position::NONE>(second_head, old_right, new_right_conjunct);
						else dst_outer[dst_outer.length] = substitute_head<any_node_position::NONE>(second_head, right, new_right_conjunct);
						free(*new_right_conjunct); if (new_right_conjunct->reference_count == 0) free(new_right_conjunct);
						if (dst_outer[dst_outer.length] == nullptr) {
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							return false;
						}
						dst_outer.length++;
						could_have_wide_scope = true;
					}
					free_all(new_heads);
				}
				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
				free(*set_var); if (set_var->reference_count == 0) free(set_var);
				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);

				return (dst.length > 0);
			} else {
				return false;
			}
		});
}

template<unsigned int ArgConstant, bool InvertArg>
hol_term* invert_select_arg_without_head_process_operand(hol_term* operand, unsigned int new_head_variable)
{
	if (operand->type != hol_term_type::EXISTS) {
		fprintf(stderr, "invert_select_arg_without_head_process_operand ERROR: Expected `operand` to be an existentially-quantified expression.\n");
		return nullptr;
	}

	hol_term* new_inner_right;
	if (InvertArg) {
		new_inner_right = hol_term::new_equals(
				hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(new_head_variable)),
				hol_term::new_variable(operand->quantifier.variable));
	} else {
		new_inner_right = hol_term::new_equals(
				hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(operand->quantifier.variable)),
				hol_term::new_variable(new_head_variable));
	}
	if (new_inner_right == nullptr)
		return nullptr;

	hol_term* new_conjunct;
	hol_term* inner_operand = operand->quantifier.operand;
	if (inner_operand->type == hol_term_type::ANY_ARRAY && inner_operand->any_array.oper == hol_term_type::AND) {
		new_conjunct = hol_term::new_exists(inner_operand->quantifier.variable, hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all,
				make_array_view(inner_operand->any_array.any.operands, inner_operand->any_array.any.length),
				make_array_view(inner_operand->any_array.left.operands, inner_operand->any_array.left.length),
				make_appended_array_view(make_array_view(inner_operand->any_array.right.operands, inner_operand->any_array.right.length), new_inner_right)));
		if (new_conjunct == nullptr) {
			free(*new_inner_right); free(new_inner_right);
			return nullptr;
		}
		inner_operand->any_array.all->reference_count++;
		for (unsigned int i = 0; i < inner_operand->any_array.any.length; i++)
			inner_operand->any_array.any.operands[i]->reference_count++;
		for (unsigned int i = 0; i < inner_operand->any_array.left.length; i++)
			inner_operand->any_array.left.operands[i]->reference_count++;
		for (unsigned int i = 0; i < inner_operand->any_array.right.length; i++)
			inner_operand->any_array.right.operands[i]->reference_count++;
	} else if (inner_operand->type == hol_term_type::AND) {
		new_conjunct = hol_term::new_exists(operand->quantifier.variable, hol_term::new_and(
				make_appended_array_view(make_array_view(inner_operand->array.operands, inner_operand->array.length), new_inner_right)));
		if (new_conjunct == nullptr) {
			free(*new_inner_right); free(new_inner_right);
			return nullptr;
		}
		for (unsigned int i = 0; i < inner_operand->array.length; i++)
			inner_operand->array.operands[i]->reference_count++;
	} else {
		new_conjunct = hol_term::new_exists(operand->quantifier.variable, hol_term::new_and(inner_operand, new_inner_right));
		if (new_conjunct == nullptr) {
			free(*new_inner_right); free(new_inner_right);
			return nullptr;
		}
		inner_operand->reference_count++;
	}
	return new_conjunct;
}

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant, bool InvertArg>
inline bool invert_select_arg_without_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates, true>, make_array_finder(find_head<built_in_predicates>), no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_predicate_index.position != head_position::NONE)
				first_head = first_head->any.included;

			unsigned int new_head_variable;
			if (first_head->type == hol_term_type::EXISTS) {
				new_head_variable = first_head->quantifier.variable;
			} else {
				new_head_variable = ++max_variable;
			}

			hol_term* new_conjunct;
			if (second_head->type == hol_term_type::ANY_ARRAY) {
				hol_term* all = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head->any_array.all, new_head_variable);
				if (all == nullptr) return false;

				array<hol_term*> left(second_head->array.length);
				for (unsigned int i = 0; i < second_head->array.length; i++) {
					left[i] = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head->any_array.left.operands[i], new_head_variable);
					if (left[i] == nullptr) {
						free(*all); if (all->reference_count == 0) free(all);
						free_all(left); return false;
					}
					left.length++;
				}

				array<hol_term*> right(second_head->array.length);
				for (unsigned int i = 0; i < second_head->array.length; i++) {
					right[i] = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head->any_array.right.operands[i], new_head_variable);
					if (right[i] == nullptr) {
						free(*all); if (all->reference_count == 0) free(all);
						free_all(left); free_all(right); return false;
					}
					right.length++;
				}

				array<hol_term*> any(second_head->array.length);
				for (unsigned int i = 0; i < second_head->array.length; i++) {
					any[i] = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head->any_array.any.operands[i], new_head_variable);
					if (any[i] == nullptr) {
						free(*all); if (all->reference_count == 0) free(all);
						free_all(left); free_all(right); free_all(any); return false;
					}
					any.length++;
				}

				new_conjunct = hol_term::new_any_array(second_head->any_array.oper, all,
						make_array_view(any.data, any.length), make_array_view(left.data, left.length), make_array_view(right.data, right.length));
				if (new_conjunct == nullptr) {
					free(*all); if (all->reference_count == 0) free(all);
					free_all(left); free_all(right); free_all(any); return false;
				}

			} else if (second_head->type == hol_term_type::AND || second_head->type == hol_term_type::OR) {
				array<hol_term*> new_operands(second_head->array.length);
				for (unsigned int i = 0; i < second_head->array.length; i++) {
					new_operands[i] = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head->array.operands[i], new_head_variable);
					if (new_operands[i] == nullptr) {
						free_all(new_operands);
						return false;
					}
					new_operands.length++;
				}

				if (second_head->type == hol_term_type::AND)
					new_conjunct = hol_term::new_and(make_array_view(new_operands.data, new_operands.length));
				else new_conjunct = hol_term::new_or(make_array_view(new_operands.data, new_operands.length));
				if (new_conjunct == nullptr) {
					free_all(new_operands);
					return false;
				}

			} else {
				new_conjunct = invert_select_arg_without_head_process_operand<ArgConstant, InvertArg>(second_head, new_head_variable);
				if (new_conjunct == nullptr) return false;
			}

			hol_term* new_head;
			if (ConjunctIndex >= 0) {
				new_head = hol_term::new_exists(new_head_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(&HOL_ANY, (unsigned int) ConjunctIndex), new_conjunct),
						make_array_view((hol_term**) nullptr, 0)));
				if (new_head == nullptr) {
					free(*new_conjunct); free(new_conjunct);
					return false;
				}
				HOL_ANY.reference_count += ConjunctIndex + 1;
			} else {
				unsigned int index = -(ConjunctIndex + 1);
				new_head = hol_term::new_exists(new_head_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0),
						make_prepended_array_view(new_conjunct, make_repeated_array_view(&HOL_ANY, (unsigned int) index))));
				if (new_head == nullptr) {
					free(*new_conjunct); free(new_conjunct);
					return false;
				}
				HOL_ANY.reference_count += index + 1;
			}

			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

template<size_t N>
inline bool invert_apply_tense_predicate(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second,
		const unsigned int (&input_tense_predicates)[N],
		const unsigned int (&output_tense_predicates)[N])
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[input_tense_predicates,output_tense_predicates](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			unsigned int predicate_variable;
			if (second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) {
				predicate_variable = max_variable + 1;
			} else {
#if !defined(NDEBUG)
				if (second_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_apply_tense_predicate WARNING: Expected `second_head` to be an existential quantification.\n");
#endif
				predicate_variable = second_head->quantifier.variable;
			}

			if (!map_tense_predicate<false>(dst, second_head, predicate_variable, second_predicate_index, output_tense_predicates, input_tense_predicates))
				return false;
			if (!dst_outer.ensure_capacity(dst.length)) {
				free_all(dst);
				return false;
			}
			for (unsigned int i = 0; i < dst.length; i++)
				dst_outer[i] = &HOL_ZERO;
			dst_outer.length = dst.length;
			HOL_ZERO.reference_count += dst.length;
			return (dst.length > 0);
		});
}

inline bool invert_remove_perfect(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			unsigned int predicate_variable;
			if (second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) {
				predicate_variable = max_variable + 1;
			} else {
#if !defined(NDEBUG)
				if (second_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_remove_perfect WARNING: Expected `second_head` to be an existential quantification.\n");
#endif
				predicate_variable = second_head->quantifier.variable;
			}

			if (!map_tense_predicate<false>(dst, second_head, predicate_variable, second_predicate_index, NON_PERFECT_PREDICATES, PERFECT_PREDICATES))
				return false;
			if (!dst_outer.ensure_capacity(dst.length)) {
				free_all(dst);
				return false;
			}
			for (unsigned int i = 0; i < dst.length; i++)
				dst_outer[i] = &HOL_ZERO;
			dst_outer.length = dst.length;
			HOL_ZERO.reference_count += dst.length;
			return (dst.length > 0);
		});
}

inline bool invert_remove_not(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			if (!dst.ensure_capacity(dst.length + 1))
				return false;

			dst[dst.length] = hol_term::new_not(second_head);
			if (dst[dst.length] == nullptr)
				return false;
			second_head->reference_count++;
			dst.length++;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

inline bool invert_remove_predicative_not(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	if (second->type != hol_term_type::LAMBDA)
		return false;
	unsigned int lambda_variable = second->quantifier.variable;
	second = second->quantifier.operand;

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), predicative_head_finder<built_in_predicates>(lambda_variable), no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = second_head->any.included;
			if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
				second_head = second_head->any.included;

#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS) {
				fprintf(stderr, "invert_remove_predicative_not ERROR: Expected existential quantification of set.\n");
				return false;
			}
#endif

			hol_term* last;
			hol_term* operand = second_head->quantifier.operand;
			if (operand->type == hol_term_type::ANY_ARRAY) {
				if (operand->any_array.right.length == 0)
					return false;
				last = operand->any_array.right.operands[operand->any_array.right.length - 1];
			} else if (operand->type == hol_term_type::AND) {
				last = operand->array.operands[operand->array.length - 1];
			} else {
				last = operand;
			}

			if (last->type == hol_term_type::ANY_RIGHT && last->any.included != nullptr) {
				bool has_wide_scope = false;
				hol_term* inner_last = last->any.included;
				while (inner_last->type == hol_term_type::NOT)
					inner_last = inner_last->unary.operand;
				if (inner_last->type == hol_term_type::UNARY_APPLICATION && inner_last->binary.left->type == hol_term_type::CONSTANT && inner_last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
					has_wide_scope = true;
				while (inner_last->type == hol_term_type::NOT)
					inner_last = inner_last->unary.operand;

				if (inner_last->type != hol_term_type::UNARY_APPLICATION && !has_wide_scope) {
					/* make sure the tree correctly excludes everything in `last->any.excluded_trees` */
					array<hol_term*> new_last(4);
					new_last[new_last.length++] = inner_last;
					inner_last->reference_count++;
					for (unsigned int i = 0; i < last->any.excluded_tree_count; i++) {
						array<hol_term*> temp(4);
						for (hol_term* term : new_last)
							subtract<built_in_predicates>(temp, term, last->any.excluded_trees[i]);
						free_all(new_last);
						swap(temp, new_last);
					}

					/* try adding a negated wide scope */
					if (!dst.ensure_capacity(dst.length + new_last.length + 2) || !dst_outer.ensure_capacity(dst.length + new_last.length + 2)) {
						free_all(new_last);
						return false;
					}
					for (hol_term* term : new_last) {
						hol_term* negated = hol_term::new_not(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), term));
						if (negated == nullptr) return false;
						term->reference_count++;
						hol_term* new_head = substitute_head<any_node_position::NONE>(second_head, last, negated);
						free(*negated); if (negated->reference_count == 0) free(negated);
						if (new_head == nullptr) {
							free_all(new_last);
							return false;
						}
						dst[dst.length++] = new_head;
						dst_outer[dst_outer.length++] = &HOL_ZERO;
						HOL_ZERO.reference_count++;
					}
					free_all(new_last);
				}

				last = last->any.included;
			} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
				/* try adding the negation here */
				hol_term* negated = hol_term::new_not(last);
				if (negated == nullptr) return false;
				last->reference_count++;
				hol_term* new_head = substitute_head<any_node_position::NONE>(second_head, last, negated);
				free(*negated); if (negated->reference_count == 0) free(negated);
				if (new_head == nullptr)
					return false;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;

				last = last->binary.right;
			}

			if (last->type == hol_term_type::NOT) {
				/* try adding the negation here */
				hol_term* negated = hol_term::new_not(last);
				if (negated == nullptr) return false;
				last->reference_count++;
				hol_term* new_head = substitute_head<any_node_position::NONE>(second_head, last, negated);
				free(*negated); if (negated->reference_count == 0) free(negated);
				if (new_head == nullptr)
					return false;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
				return true;
			}

			if (last->type == hol_term_type::UNARY_APPLICATION)
				return true;

			/* try adding the negation here */
			hol_term* negated = hol_term::new_not(last);
			if (negated == nullptr) return false;
			last->reference_count++;
			hol_term* new_head = substitute_head<any_node_position::NONE>(second_head, last, negated);
			free(*negated); if (negated->reference_count == 0) free(negated);
			if (new_head == nullptr)
				return false;
			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

inline bool invert_predicate_only(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	/* we expect the head of `first` to be at the root */
	hol_term* head; head_index predicate_index;
	find_head<built_in_predicates>(first, head, predicate_index);
	if (head == nullptr) return false;

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_root, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			array<hol_term*> new_heads(4);
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				unsigned int predicate_variable = max_variable + 1;
				new_heads[new_heads.length] = hol_term::new_exists(predicate_variable, hol_term::new_apply(second_head, hol_term::new_variable(predicate_variable)));
				if (new_heads[new_heads.length] == nullptr)
					return false;
				second_head->reference_count++;
				new_heads.length++;

			} else if (first_head->type == hol_term_type::EXISTS && (first_head->quantifier.operand->type == hol_term_type::AND || first_head->quantifier.operand->type == hol_term_type::ANY_ARRAY)) {
				hol_term* new_head = hol_term::new_exists(first_head->quantifier.variable, hol_term::new_apply(second_head, hol_term::new_variable(first_head->quantifier.variable)));
				if (new_head == nullptr) return false;
				second_head->reference_count++;
				intersect<built_in_predicates>(new_heads, first_head, new_head);
				free(*new_head); if (new_head->reference_count == 0) free(new_head);

			} else {
				fprintf(stderr, "invert_predicate_only ERROR: Expected `first_head` to be an existentially quantified conjunction.\n");
				return false;
			}

			if (!dst.ensure_capacity(dst.length + new_heads.length) || !dst_outer.ensure_capacity(dst_outer.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}

			hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
			if (any_wide_scope == nullptr) {
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
			if (duplicate_wide_scopes == nullptr) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				free(*any_wide_scope); free(any_wide_scope);
				return false;
			}

			hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
			if (excluded_universal == nullptr) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				return false;
			}
			any_wide_scope->reference_count++;

			for (hol_term* new_head : new_heads) {
				dst[dst.length++] = new_head;
				new_head->reference_count++;

				hol_term* excluded_trees[2];
				excluded_trees[0] = excluded_universal;
				excluded_trees[1] = duplicate_wide_scopes;
				dst_outer[dst_outer.length] = hol_term::new_any_right(&HOL_ZERO, excluded_trees, array_length(excluded_trees));
				if (dst_outer[dst_outer.length] == nullptr) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
					free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
					return false;
				}
				excluded_universal->reference_count++;
				duplicate_wide_scopes->reference_count++;
				HOL_ZERO.reference_count++;
				dst_outer.length++;
			}
			for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
			free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
			free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
			return (dst.length > 0);
		});
}

inline bool invert_predicate(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	/* we expect the head of `first` to be at the root */
	hol_term* head; head_index predicate_index;
	find_head<built_in_predicates>(first, head, predicate_index);
	if (head == nullptr) return false;

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_root, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			array<hol_term*> new_heads(4);
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				unsigned int predicate_variable = max_variable + 1;
				hol_term* head_var = hol_term::new_variable(predicate_variable);
				if (head_var == nullptr) return false;
				constexpr unsigned int excluded_tree_count = 4;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_exists(predicate_variable, &HOL_ANY));
				if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
				if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
				if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
					if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
					if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
					if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
					free(*head_var); free(head_var);
					return false;
				}
				free(*head_var);

				hol_term* head_conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
				if (head_conjunct == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return false;
				}

				hol_term* predicate = hol_term::new_apply(second_head, head_var);
				if (predicate == nullptr) {
					free(*head_conjunct); free(head_conjunct);
					return false;
				}
				second_head->reference_count++;
				head_var->reference_count++;

				array<hol_term*> new_predicates(2);
				intersect<built_in_predicates>(new_predicates, predicate, head_conjunct);
				free(*predicate); if (predicate->reference_count == 0) free(predicate);
				if (!new_heads.ensure_capacity(new_heads.length + new_predicates.length)) {
					free_all(new_predicates);
					free(*head_conjunct); free(head_conjunct);
					return false;
				}
				for (hol_term* new_predicate : new_predicates) {
					new_heads[new_heads.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, head_conjunct,
							make_array_view(&new_predicate, 1), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
					if (new_heads[new_heads.length] == nullptr) {
						free_all(new_heads); free_all(new_predicates);
						free(*head_conjunct); free(head_conjunct);
						return false;
					}
					head_conjunct->reference_count++;
					new_predicate->reference_count++;
					new_heads.length++;
				}
				free_all(new_predicates);
				free(*head_conjunct); if (head_conjunct->reference_count == 0) free(head_conjunct);

			} else if (first_head->type == hol_term_type::EXISTS && first_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* first_head_operand = first_head->quantifier.operand;

				unsigned int first_index = 0;
				if (first_predicate_index.position == head_position::LEFT)
					first_index = first_predicate_index.index;
				else if (first_predicate_index.position == head_position::RIGHT)
					first_index = first_head_operand->array.length - first_predicate_index.index - 1;
				hol_term* first_predicate = first_head_operand->array.operands[first_index];
				array<hol_term*> new_predicates(8);
				hol_term* new_predicate = hol_term::new_apply(second_head, hol_term::new_variable(first_head->quantifier.variable));
				if (new_predicate == nullptr) return false;
				second_head->reference_count++;
				intersect<built_in_predicates>(new_predicates, first_predicate, new_predicate);
				free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);

#if !defined(NDEBUG)
				if (new_predicates.length == 0)
					fprintf(stderr, "invert_predicate WARNING: `new_predicates` is empty.\n");
#endif

				if (!new_heads.ensure_capacity(new_heads.length + new_predicates.length)) {
					for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				for (hol_term* new_predicate : new_predicates) {
					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					conjunction->type = hol_term_type::AND;
					conjunction->reference_count = 1;
					conjunction->array.length = first_head_operand->array.length;
					conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first_head_operand->array.length);
					if (conjunction->array.operands == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(conjunction); return false;
					}
					for (unsigned int i = 0; i < conjunction->array.length; i++) {
						if (i == first_index) {
							conjunction->array.operands[i] = new_predicate;
						} else {
							conjunction->array.operands[i] = first_head_operand->array.operands[i];
						}
						conjunction->array.operands[i]->reference_count++;
					}

					new_heads[new_heads.length] = hol_term::new_exists(first_head->quantifier.variable, conjunction);
					if (new_heads[new_heads.length] == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunction); free(conjunction);
						return false;
					}
					new_heads.length++;
				}
				for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }

			} else if (first_head->type == hol_term_type::EXISTS && first_head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
				hol_term* first_head_operand = first_head->quantifier.operand;

				hol_term* first_predicate = nullptr;
				if (first_predicate_index.position == head_position::LEFT)
					first_predicate = first_head_operand->any_array.left.operands[first_predicate_index.index];
				else if (first_predicate_index.position == head_position::RIGHT)
					first_predicate = first_head_operand->any_array.right.operands[first_head_operand->any_array.right.length - first_predicate_index.index - 1];
				else if (first_predicate_index.position == head_position::ANY)
					first_predicate = first_head_operand->any_array.any.operands[first_predicate_index.index];

				array<hol_term*> new_predicates(8);
				hol_term* new_predicate = hol_term::new_apply(second_head, hol_term::new_variable(first_head->quantifier.variable));
				if (new_predicate == nullptr) return false;
				second_head->reference_count++;
				intersect<built_in_predicates>(new_predicates, first_predicate, new_predicate);
				free(*new_predicate); if (new_predicate->reference_count == 0) free(new_predicate);

#if !defined(NDEBUG)
				if (new_predicates.length == 0)
					fprintf(stderr, "invert_predicate WARNING: `new_predicates` is empty.\n");
#endif

				if (!new_heads.ensure_capacity(new_heads.length + new_predicates.length)) {
					for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}
				for (hol_term* new_predicate : new_predicates) {
					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = first_head_operand->any_array.all;
					conjunction->any_array.all->reference_count++;
					conjunction->any_array.left.length = 0;
					conjunction->any_array.left.operands = nullptr;
					conjunction->any_array.right.length = 0;
					conjunction->any_array.right.operands = nullptr;
					conjunction->any_array.any.length = 0;
					conjunction->any_array.any.operands = nullptr;

					conjunction->any_array.left.length = first_head_operand->any_array.left.length;
					if (conjunction->any_array.left.length == 0) {
						conjunction->any_array.left.operands = nullptr;
					} else {
						conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.left.length);
						if (conjunction->any_array.left.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
							if (first_predicate_index.position == head_position::LEFT && first_predicate_index.index == i) {
								conjunction->any_array.left.operands[i] = new_predicate;
							} else {
								conjunction->any_array.left.operands[i] = first_head_operand->any_array.left.operands[i];
							}
							conjunction->any_array.left.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.right.length = first_head_operand->any_array.right.length;
					if (conjunction->any_array.right.length == 0) {
						conjunction->any_array.right.operands = nullptr;
					} else {
						conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.right.length);
						if (conjunction->any_array.right.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
							if (first_predicate_index.position == head_position::RIGHT && first_predicate_index.index == conjunction->any_array.right.length - i - 1) {
								conjunction->any_array.right.operands[i] = new_predicate;
							} else {
								conjunction->any_array.right.operands[i] = first_head_operand->any_array.right.operands[i];
							}
							conjunction->any_array.right.operands[i]->reference_count++;
						}
					}
					conjunction->any_array.any.length = first_head_operand->any_array.any.length;
					if (conjunction->any_array.any.length == 0) {
						conjunction->any_array.any.operands = nullptr;
					} else {
						conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * conjunction->any_array.any.length);
						if (conjunction->any_array.any.operands == nullptr) {
							for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunction); free(conjunction);
							return false;
						}
						for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
							if (first_predicate_index.position == head_position::ANY && first_predicate_index.index == i) {
								conjunction->any_array.any.operands[i] = new_predicate;
							} else {
								conjunction->any_array.any.operands[i] = first_head_operand->any_array.any.operands[i];
							}
							conjunction->any_array.any.operands[i]->reference_count++;
						}
					}

					new_heads[new_heads.length] = hol_term::new_exists(first_head->quantifier.variable, conjunction);
					if (new_heads[new_heads.length] == nullptr) {
						for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunction); free(conjunction);
						return false;
					}
					new_heads.length++;
				}
				for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
			} else {
				fprintf(stderr, "invert_predicate ERROR: Expected `first_head` to be an existentially quantified conjunction.\n");
				return false;
			}

			if (!dst.ensure_capacity(dst.length + new_heads.length) || !dst_outer.ensure_capacity(dst_outer.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}

			hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
			if (any_wide_scope == nullptr) {
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
			if (duplicate_wide_scopes == nullptr) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				free(*any_wide_scope); free(any_wide_scope);
				return false;
			}

			hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
			if (excluded_universal == nullptr) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				return false;
			}
			any_wide_scope->reference_count++;

			for (hol_term* new_head : new_heads) {
				dst[dst.length++] = new_head;
				new_head->reference_count++;

				hol_term* excluded_trees[2];
				excluded_trees[0] = excluded_universal;
				excluded_trees[1] = duplicate_wide_scopes;
				dst_outer[dst_outer.length] = hol_term::new_any_right(&HOL_ZERO, excluded_trees, array_length(excluded_trees));
				if (dst_outer[dst_outer.length] == nullptr) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
					free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
					return false;
				}
				excluded_universal->reference_count++;
				duplicate_wide_scopes->reference_count++;
				HOL_ZERO.reference_count++;
				dst_outer.length++;
			}
			for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
			free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
			free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
			return (dst.length > 0);
		});
}

inline bool invert_predicate_and_tense(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	head_index predicate_index; no_op apply;
	auto find_array_head = make_array_finder(find_head<built_in_predicates>);
	hol_term* first_head = find_head(first, predicate_index, find_array_head, apply);
	if (first_head == nullptr)
		return false;

	if (second->type != hol_term_type::EXISTS)
		return false;
	unsigned int predicate_variable = second->quantifier.variable;

	hol_term* first_outer = substitute_head<any_node_position::NONE>(first, first_head, &HOL_ZERO);
	if (first_outer == nullptr)
		return false;

	hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
	if (any_wide_scope == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		return false;
	}
	HOL_ANY.reference_count++;

	hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
	if (duplicate_wide_scopes == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		free(*any_wide_scope); free(any_wide_scope); return false;
	}

	hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
	if (excluded_universal == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes); return false;
	}
	any_wide_scope->reference_count++;

	hol_term* excluded_negation = hol_term::new_not(hol_term::new_exists(predicate_variable, &HOL_ANY));
	if (excluded_negation == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		free(*excluded_universal); free(excluded_universal); return false;
	}
	HOL_ANY.reference_count++;

	hol_term* excluded_trees[3];
	excluded_trees[0] = excluded_universal;
	excluded_trees[1] = duplicate_wide_scopes;
	excluded_trees[2] = excluded_negation;
	hol_term* new_second_outer = hol_term::new_any_right(&HOL_ZERO, excluded_trees, array_length(excluded_trees));
	if (new_second_outer == nullptr) {
		free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		free(*excluded_universal); free(excluded_universal);
		free(*excluded_negation); free(excluded_negation);
		return false;
	}
	HOL_ZERO.reference_count++;

	array<hol_term*> new_outer(4);
	intersect<built_in_predicates>(new_outer, first_outer, new_second_outer);
	free(*first_outer); if (first_outer->reference_count == 0) free(first_outer);
	free(*new_second_outer); if (new_second_outer->reference_count == 0) free(new_second_outer);
	if (new_outer.length == 0)
		return false;

	array<hol_term*> inverted_logical_forms(new_outer.length);
	for (hol_term* outer : new_outer) {
		hol_term* new_term = substitute_head<any_node_position::NONE>(outer, &HOL_ZERO, second);
		if (new_term == nullptr) {
			free_all(new_outer); free_all(inverted_logical_forms);
			return false;
		}
		inverted_logical_forms[inverted_logical_forms.length++] = new_term;
	}
	free_all(new_outer);

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * inverted_logical_forms.length);
	for (unsigned int i = 0; i < inverted_logical_forms.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = inverted_logical_forms[i];
	}
	inverse_count = inverted_logical_forms.length;
	return true;
}

inline bool invert_select_predicate_in_set(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else if (first->type == hol_term_type::ANY) {
		unsigned int max_variable = 0;
		max_bound_variable(*first, max_variable);
		lambda_variable = ++max_variable;
	} else {
		return false;
	}
	if (!invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), find_root, no_op(),
		[lambda_variable](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			array<hol_term*> new_heads(2);
			select_predicate_in_set<false, false, false, false>(new_heads, first_head, second_head, max_variable, lambda_variable);
			if (new_heads.length == 0)
				return false;

			for (hol_term* new_head : new_heads) {
				if (!dst.ensure_capacity(dst.length + 2) || !dst_outer.ensure_capacity(dst_outer.length + 2)) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					return false;
				}

				hol_term* right = new_head->quantifier.operand->array.operands[new_head->quantifier.operand->array.length - 1];
				if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length++] = &HOL_ZERO;
					HOL_ZERO.reference_count++;
					continue;
				}

				hol_term* operand;
				if (right->type == hol_term_type::ANY_QUANTIFIER) {
					operand = right->any_quantifier.operand;
				} else if (right->type == hol_term_type::FOR_ALL || right->type == hol_term_type::EXISTS) {
					operand = right->quantifier.operand;
				} else {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length++] = &HOL_ZERO;
					HOL_ZERO.reference_count++;
					continue;
				}

				hol_term* inner_right;
				if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.right.length != 0) {
					inner_right = operand->any_array.right.operands[operand->any_array.right.length - 1];
				} else if (operand->type == hol_term_type::AND) {
					inner_right = operand->array.operands[operand->array.length - 1];
				} else if (operand->type == hol_term_type::IF_THEN) {
					inner_right = operand->binary.right;
				} else {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length++] = &HOL_ZERO;
					HOL_ZERO.reference_count++;
					continue;
				}

				if (inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length++] = &HOL_ZERO;
					HOL_ZERO.reference_count++;
					continue;
				} else {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length++] = &HOL_ZERO;
					HOL_ZERO.reference_count++;

					dst[dst.length++] = new_head;
					new_head->reference_count++;
					dst_outer[dst_outer.length] = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ZERO);
					if (dst_outer[dst_outer.length] == nullptr) {
						free_all(new_heads);
						return false;
					}
					HOL_ZERO.reference_count++;
					dst_outer.length++;
				}
			}
			for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
			return true;
		})) return false;
	for (unsigned int i = 0; i < inverse_count; i++) {
		hol_term* new_inverse = hol_term::new_lambda(lambda_variable, inverse[i].root);
		if (new_inverse == nullptr) {
			for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
			free(inverse); return false;
		}
		inverse[i].root = new_inverse;
	}
	return true;
}

inline bool invert_mark_wide_scope(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	bool found_universal_quantifier = false;
	auto second_apply = [&found_universal_quantifier](hol_term* term) {
		found_universal_quantifier |= (term->type == hol_term_type::FOR_ALL || (term->type == hol_term_type::ANY_QUANTIFIER && term->any_quantifier.quantifier == hol_quantifier_type::FOR_ALL));
	};

	head_index second_predicate_index;
	hol_term* second_head = find_head(second, second_predicate_index, find_head_or_unary_application<built_in_predicates, (unsigned int) built_in_predicates::WIDE_SCOPE>, second_apply);
	if (second_head == nullptr || found_universal_quantifier)
		return false;

	while (second_head->type == hol_term_type::NOT)
		second_head = second_head->unary.operand;

	hol_term* new_head;
	if (second_head->type == hol_term_type::UNARY_APPLICATION && second_head->binary.left->type == hol_term_type::CONSTANT && second_head->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
		if (second_head->binary.right->type != hol_term_type::FOR_ALL)
			return false;
		new_head = second_head->binary.right;
		second_head->binary.right->reference_count++;
	} else if (second_head->type == hol_term_type::EXISTS) {
		new_head = second_head;
		second_head->reference_count++;
	} else {
		return false;
	}

	hol_term* new_term = substitute_head<any_node_position::NONE>(second, second_head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	if (new_term == nullptr)
		return false;

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * 1);
	if (inverse == nullptr) {
		fprintf(stderr, "invert_mark_wide_scope ERROR: Out of memory.\n");
		return false;
	}
	inverse[0].flags = flags;
	inverse[0].root = new_term;
	inverse_count = 1;
	return true;
}

template<bool WideScope>
inline bool invert_require_narrow_or_wide_scope(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	unsigned int lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		lambda_variable = first->quantifier.variable;
	} else if (first->type == hol_term_type::ANY) {
		max_bound_variable(*first, max_variable);
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	head_index predicate_index; no_op apply;
	predicative_head_finder<built_in_predicates> try_find_head(lambda_variable);
	hol_term* second_head = find_head(second, predicate_index, try_find_head, apply);
	if (second_head == nullptr)
		return false;

	if ((second_head->type == hol_term_type::ANY || second_head->type == hol_term_type::ANY_RIGHT) && second_head->any.included != nullptr)
		second_head = second_head->any.included;

#if !defined(NDEBUG)
	if (second_head->type != hol_term_type::EXISTS)
		fprintf(stderr, "invert_require_narrow_or_wide_scope WARNING: Expected `second_head` to be an existential quantification.\n");
#endif

	hol_term* operand = second_head->quantifier.operand;

	hol_term* last;
	if (operand->type == hol_term_type::ANY_ARRAY) {
		if (operand->any_array.right.length == 0) {
			fprintf(stderr, "invert_require_narrow_or_wide_scope ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
			return (hol_term*) nullptr;
		}
		last = operand->any_array.right.operands[operand->any_array.right.length - 1];
	} else if (operand->type == hol_term_type::AND) {
		last = operand->array.operands[operand->array.length - 1];
	} else {
		last = operand;
	}

	bool could_have_negation = false;
	unsigned int negation_count = 0;
	if ((last->type == hol_term_type::ANY || last->type == hol_term_type::ANY_RIGHT) && last->any.included != nullptr) {
		last = last->any.included;
		could_have_negation = true;
	} if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
		last = last->binary.right;
		if (WideScope) {
			return intersect(inverse, inverse_count, flags, first, second);
		} else {
			return false;
		}
	}

	while (last->type == hol_term_type::NOT) {
		last = operand->unary.operand;
		negation_count++;
	}

	array<hol_term*> inverted_logical_forms(8);
	hol_term* expected_second;
	if (WideScope) {
		hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
		if (any_wide_scope == nullptr)
			return false;
		HOL_ANY.reference_count++;

		hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
		if (duplicate_wide_scopes == nullptr) {
			free(*any_wide_scope); free(any_wide_scope);
			return false;
		}

		hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
		if (excluded_universal == nullptr) {
			free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
			return false;
		}
		any_wide_scope->reference_count++;

		hol_term* expected_last;
		if ((last->type == hol_term_type::ANY || last->type == hol_term_type::ANY_RIGHT) && last->any.included != nullptr
		 && last->any.included->type == hol_term_type::UNARY_APPLICATION && last->any.included->binary.left->type == hol_term_type::VARIABLE
		 && last->any.included->binary.left->variable == lambda_variable)
		{
			expected_last = last;
			last->reference_count++;
			could_have_negation = true;
		} else {
			hol_term* inner_last;
			if (last->type == hol_term_type::AND && last->array.length == 2) {
				inner_last = last->array.operands[1];
			} else if (last->type == hol_term_type::FOR_ALL && last->quantifier.operand->type == hol_term_type::IF_THEN) {
				inner_last = last->quantifier.operand->binary.right;
			} else if (last->type == hol_term_type::EXISTS && last->quantifier.operand->type == hol_term_type::AND && last->quantifier.operand->array.length == 2) {
				inner_last = last->quantifier.operand->array.operands[1];
			} else if (last->type == hol_term_type::ANY_QUANTIFIER && last->any_quantifier.operand->type == hol_term_type::IF_THEN) {
				inner_last = last->quantifier.operand->binary.right;
			} else if (last->type == hol_term_type::ANY_QUANTIFIER && last->any_quantifier.operand->type == hol_term_type::AND && last->quantifier.operand->array.length == 2) {
				inner_last = last->quantifier.operand->array.operands[1];
			} else {
				free(*excluded_universal); free(excluded_universal);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				return false;
			}

			if ((inner_last->type == hol_term_type::ANY || inner_last->type == hol_term_type::ANY_RIGHT) && inner_last->any.included != nullptr) {
				expected_last = substitute_head<any_node_position::NONE>(last, inner_last, inner_last->any.included);
			} else {
				expected_last = last;
				last->reference_count++;
			}
		}

		if (could_have_negation) {
			hol_term* temp = hol_term::new_any_right(expected_last);
			if (temp == nullptr) {
				free(*excluded_universal); free(excluded_universal);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				free(*expected_last); if (expected_last->reference_count == 0) free(expected_last);
				return false;
			}
			expected_last = temp;
		} for (unsigned int i = 0; i < negation_count; i++) {
			hol_term* temp = hol_term::new_not(expected_last);
			if (temp == nullptr) {
				free(*excluded_universal); free(excluded_universal);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				free(*expected_last); if (expected_last->reference_count == 0) free(expected_last);
				return false;
			}
			expected_last = temp;
		}

		hol_term* excluded_trees[2];
		excluded_trees[0] = excluded_universal;
		excluded_trees[1] = duplicate_wide_scopes;
		expected_second = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), expected_last), excluded_trees, 2);
		if (expected_second == nullptr) {
			free(*excluded_trees[0]); free(excluded_trees[0]);
			free(*excluded_trees[1]); free(excluded_trees[1]);
			free(*expected_last); if (expected_last->reference_count == 0) free(expected_last);
			return false;
		}
	} else {
		hol_term* wide_scope = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY);
		if (wide_scope == nullptr) return false;
		expected_second = hol_term::new_any_right(nullptr, &wide_scope, 1);
		if (expected_second == nullptr) {
			free(*wide_scope); free(wide_scope);
			return false;
		}
	}

	array<hol_term*> intersection(4);
	intersect<built_in_predicates>(intersection, second, expected_second);
	free(*expected_second); if (expected_second->reference_count == 0) free(expected_second);

	for (hol_term* new_second : intersection)
		intersect<built_in_predicates>(inverted_logical_forms, first, new_second);
	free_all(intersection);

	if (inverted_logical_forms.length == 0)
		return false;

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * inverted_logical_forms.length);
	for (unsigned int i = 0; i < inverted_logical_forms.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = inverted_logical_forms[i];
	}
	inverse_count = inverted_logical_forms.length;
	return true;
}

inline bool invert_remove_wide_scope(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	hol_term* excluding_term = nullptr;
	unsigned int excluded_tree_index;
	auto apply = [&excluding_term,&excluded_tree_index](hol_term* term) {
		if (term->type == hol_term_type::ANY || term->type == hol_term_type::ANY_RIGHT) {
			for (unsigned int i = 0; i < term->any.excluded_tree_count; i++) {
				hol_term* excluded = term->any.excluded_trees[i];
				if ((excluded->type == hol_term_type::ANY || excluded->type == hol_term_type::ANY_RIGHT) && excluded->any.included != nullptr && excluded->any.included->type == hol_term_type::UNARY_APPLICATION
				 && excluded->any.included->binary.left->type == hol_term_type::CONSTANT && excluded->any.included->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE
				 && (excluded->any.included->binary.right->type == hol_term_type::ANY || excluded->any.included->binary.right->type == hol_term_type::ANY_RIGHT)
				 && excluded->any.included->binary.right->any.included == nullptr)
				{
					excluding_term = term;
					excluded_tree_index = i;
					break;
				}
			}
		}
	};

	head_index predicate_index;
	auto find_array_head = make_array_finder(find_head<built_in_predicates>);
	hol_term* second_head = find_head(second, predicate_index, find_array_head, apply);
	if (second_head == nullptr)
		return false;

	hol_term* new_second;
	if (excluding_term != nullptr) {
		hol_term** new_excluded_trees;
		if (excluding_term->any.excluded_tree_count == 1) {
			new_excluded_trees = nullptr;
		} else {
			new_excluded_trees = (hol_term**) malloc(sizeof(hol_term*) * (excluding_term->any.excluded_tree_count - 1));
			if (new_excluded_trees == nullptr) {
				fprintf(stderr, "invert_remove_wide_scope ERROR: Out of memory.\n");
				return false;
			}
			unsigned int index = 0;
			for (unsigned int i = 0; i < excluding_term->any.excluded_tree_count; i++)
				if (i != excluded_tree_index) new_excluded_trees[index++] = excluding_term->any.excluded_trees[i];
		}

		hol_term* new_term;
		if (!new_hol_term(new_term)) {
			free(new_excluded_trees);
			return false;
		}
		new_term->type = excluding_term->type;
		new_term->reference_count = 1;
		new_term->any.included = excluding_term->any.included;
		new_term->any.included->reference_count++;
		new_term->any.excluded_trees = new_excluded_trees;
		new_term->any.excluded_tree_count = excluding_term->any.excluded_tree_count - 1;
		for (unsigned int i = 0; i < new_term->any.excluded_tree_count; i++)
			new_excluded_trees[i]->reference_count++;

		new_second = substitute_head<any_node_position::NONE>(second, excluding_term, new_term);
		free(*new_term); if (new_term->reference_count == 0) free(new_term);
		if (new_second == nullptr) {
			free(*new_term); free(new_term);
			return false;
		}
	} else {
		new_second = second;
		second->reference_count++;
	}

	hol_term* second_outer = substitute_head<any_node_position::LEFT>(new_second, second_head, &HOL_ZERO, true);
	free(*new_second); if (new_second->reference_count == 0) free(new_second);
	if (second_outer == nullptr)
		return false;

	hol_term* any_wide_scope = hol_term::new_any(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
	if (any_wide_scope == nullptr) {
		free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
		return false;
	}
	HOL_ANY.reference_count++;

	hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
	if (duplicate_wide_scopes == nullptr) {
		free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
		free(*any_wide_scope); free(any_wide_scope);
		return false;
	}

	hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
	if (excluded_universal == nullptr) {
		free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		return false;
	}
	any_wide_scope->reference_count++;

	hol_term* excluded_trees[2];
	excluded_trees[0] = excluded_universal;
	excluded_trees[1] = duplicate_wide_scopes;
	hol_term* wide_scope = hol_term::new_any_right(hol_term::new_apply(
				hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE),
				hol_term::new_any_right(&HOL_ZERO)),
			excluded_trees, 2);
	if (wide_scope == nullptr) {
		free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
		free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
		free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
		return false;
	}
	HOL_ZERO.reference_count++;

	array<hol_term*> intersection(4);
	intersect<built_in_predicates>(intersection, second_outer, wide_scope);
	free(*second_outer); if (second_outer->reference_count == 0) free(second_outer);
	free(*wide_scope); if (wide_scope->reference_count == 0) free(wide_scope);

	/* also consider the possibility that there is no universal quantifier */
	hol_term* any_right_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, &HOL_ANY));
	if (any_right_universal == nullptr) {
		free_all(intersection);
		return false;
	}

	hol_term* no_wide_scope = hol_term::new_any_right(&HOL_ZERO, &any_right_universal, 1);
	if (no_wide_scope == nullptr) {
		free(*any_right_universal); if (any_right_universal->reference_count == 0) free(any_right_universal);
		free_all(intersection); return false;
	}
	HOL_ZERO.reference_count++;

	intersect<built_in_predicates>(intersection, second_outer, no_wide_scope);
	free(*no_wide_scope); if (no_wide_scope->reference_count == 0) free(no_wide_scope);
	if (intersection.length == 0)
		return false;

	array<hol_term*> inverted_logical_forms(intersection.length);
	for (hol_term* term : intersection) {
		hol_term* new_outer = remove_any_nodes(term, find_zero);
		if (new_outer == nullptr) {
			free_all(inverted_logical_forms); free_all(intersection);
			return false;
		}

		hol_term* new_term = substitute_head<any_node_position::NONE>(new_outer, &HOL_ZERO, second_head);
		free(*new_outer); if (new_outer->reference_count == 0) free(new_outer);
		if (new_term == nullptr) {
			free_all(inverted_logical_forms); free_all(intersection);
			return false;
		}
		inverted_logical_forms[inverted_logical_forms.length++] = new_term;
	}
	free_all(intersection);

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * inverted_logical_forms.length);
	for (unsigned int i = 0; i < inverted_logical_forms.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = inverted_logical_forms[i];
	}
	inverse_count = inverted_logical_forms.length;
	return true;
}

inline bool invert_size(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int head_variable;
	if (first->type == hol_term_type::EXISTS) {
		head_variable = first->quantifier.variable;
	} else if ((first->type == hol_term_type::ANY || first->type == hol_term_type::ANY_RIGHT) && first->any.included != nullptr && first->any.included->type == hol_term_type::EXISTS) {
		first = first->any.included;
		head_variable = first->quantifier.variable;
	} else {
		head_variable = 1;
	}

	if (second->type == hol_term_type::INTEGER) {
		if (second->integer == 1) {
			if (!intersect(flags.index_number, grammatical_num::SINGULAR, flags.index_number))
				return false;
		} else {
			if (!intersect(flags.index_number, grammatical_num::PLURAL, flags.index_number))
				return false;
		}
	}

	hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SIZE), hol_term::new_variable(head_variable)), second));
	if (expected_head == nullptr)
		return (hol_term*) nullptr;
	second->reference_count++;

	array<hol_term*> inverted_logical_forms(2);
	intersect<built_in_predicates>(inverted_logical_forms, first, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (inverted_logical_forms.length == 0)
		return (hol_term*) nullptr;

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * inverted_logical_forms.length);
	for (unsigned int i = 0; i < inverted_logical_forms.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = inverted_logical_forms[i];
	}
	inverse_count = inverted_logical_forms.length;
	return true;
}

inline bool invert_set_size(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else if (first->type == hol_term_type::ANY) {
		unsigned int max_variable = 0;
		max_bound_variable(*first, max_variable);
		lambda_variable = ++max_variable;
	} else {
		return false;
	}
	bool result = invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), find_root, no_op(),
		[lambda_variable](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			if ((first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) && first_head->any.included != nullptr)
				first_head = first_head->any.included;

			hol_term* left = nullptr;
			unsigned int set_variable, element_variable;
			if (first_head->type == hol_term_type::ANY || first_head->type == hol_term_type::ANY_RIGHT) {
				set_variable = ++max_variable;
				element_variable = ++max_variable;
			} else {
#if !defined(NDEBUG)
				if (first_head->type != hol_term_type::EXISTS) {
					fprintf(stderr, "invert_set_size ERROR: Expected existential quantification of set.\n");
					return false;
				}
#endif
				set_variable = first_head->quantifier.variable;
				hol_term* operand = first_head->quantifier.operand;

				if (operand->type == hol_term_type::ANY || operand->type == hol_term_type::ANY_RIGHT) {
					element_variable = 0;
				} else if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
					if (operand->any_array.left.length == 0) {
						left = operand->any_array.all;
					} else {
						left = operand->any_array.left.operands[0];
					}
				} else if (operand->type == hol_term_type::AND) {
					left = operand->array.operands[0];
				} else {
					return false;
				}

				hol_term* set_definition = nullptr;
				if (left != nullptr) {
					if (left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
						element_variable = 0;
					} else if (left->type == hol_term_type::EQUALS) {
						set_definition = left->binary.right;
					} else if (left->type == hol_term_type::BINARY_APPLICATION) {
						set_definition = left->ternary.third;
					} else {
						return false;
					}
				}

				if (set_definition != nullptr) {
					if (set_definition->type == hol_term_type::ANY || set_definition->type == hol_term_type::ANY_RIGHT) {
						element_variable = 0;
					} else if (set_definition->type == hol_term_type::ANY_QUANTIFIER && has_intersection(set_definition->any_quantifier.quantifier, hol_quantifier_type::LAMBDA)) {
						element_variable = 0;
					} else if (set_definition->type == hol_term_type::LAMBDA) {
						element_variable = set_definition->quantifier.variable;
					} else {
						return false;
					}
				}

				if (element_variable == 0) {
					/* try to get the variable from the right conjunct in the set scope */
					hol_term* right = nullptr;
					if (operand->type == hol_term_type::ANY_ARRAY && operand->any_array.oper == hol_term_type::AND) {
						if (operand->any_array.right.length == 0) {
							right = operand->any_array.all;
						} else {
							right = operand->any_array.right.operands[operand->any_array.right.length - 1];
						}
					} else if (operand->type == hol_term_type::AND) {
						right = operand->array.operands[operand->array.length - 1];
					}

					if (right != nullptr) {
						if ((right->type == hol_term_type::ANY || right->type == hol_term_type::ANY_RIGHT) && right->any.included != nullptr)
							right = right->any.included;
						if (right->type == hol_term_type::NOT)
							right = right->unary.operand;

						hol_term* lambda_application = nullptr;
						if (right->type == hol_term_type::FOR_ALL && right->quantifier.operand->type == hol_term_type::IF_THEN) {
							lambda_application = right->quantifier.operand->binary.right;
						} else if (right->type == hol_term_type::EXISTS && right->quantifier.operand->type == hol_term_type::AND && right->quantifier.operand->array.length == 2) {
							lambda_application = right->quantifier.operand->array.operands[1];
						} else if (right->type == hol_term_type::AND && right->array.length == 2) {
							lambda_application = right->array.operands[1];
						} else {
							lambda_application = right;
						}

						if ((lambda_application->type == hol_term_type::ANY || lambda_application->type == hol_term_type::ANY_RIGHT) && lambda_application->any.included != nullptr)
							lambda_application = lambda_application->any.included;

						if (lambda_application->type == hol_term_type::UNARY_APPLICATION
						 && lambda_application->binary.left->type == hol_term_type::VARIABLE
						 && lambda_application->binary.left->variable == lambda_variable
						 && lambda_application->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = lambda_application->binary.right->variable;
						} else {
							return false;
						}
					}
				}
			}

			hol_term* set_var = hol_term::new_variable(set_variable);
			if (set_var == nullptr)
				return false;

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(set_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr) {
				free(*set_var); free(set_var);
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (expected_conjunct == nullptr) {
				free(*set_var); free(set_var);
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			array<hol_term*> expected_lefts(2);
			if (left == nullptr || left->type == hol_term_type::EQUALS || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
				expected_lefts[expected_lefts.length] = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, hol_term::new_true()));
				if (expected_lefts[expected_lefts.length] == nullptr) {
					free(*set_var); free(set_var);
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				set_var->reference_count++;
				expected_lefts.length++;
			} if (left == nullptr || left->type == hol_term_type::BINARY_APPLICATION || left->type == hol_term_type::ANY || left->type == hol_term_type::ANY_RIGHT) {
				expected_lefts[expected_lefts.length] = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SUBSET), set_var, hol_term::new_lambda(element_variable, hol_term::new_true()));
				if (expected_lefts[expected_lefts.length] == nullptr) {
					free_all(expected_lefts);
					free(*set_var); if (set_var->reference_count == 0) free(set_var);
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				set_var->reference_count++;
				expected_lefts.length++;
			}

			if (!dst.ensure_capacity(dst.length + expected_lefts.length) || !dst_outer.ensure_capacity(dst_outer.length + expected_lefts.length)) {
				free_all(expected_lefts);
				free(*set_var); if (set_var->reference_count == 0) free(set_var);
				free(*expected_conjunct); free(expected_conjunct);
				return false;
			}
			for (hol_term* expected_left : expected_lefts) {
				dst[dst.length] = hol_term::new_exists(set_variable, hol_term::new_and(
						expected_left,
						hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::SIZE), set_var), second_head),
						expected_conjunct));
				if (dst[dst.length] == nullptr) {
					free_all(dst); free_all(expected_lefts);
					free(*set_var); if (set_var->reference_count == 0) free(set_var);
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				expected_left->reference_count++;
				set_var->reference_count++;
				expected_conjunct->reference_count++;
				second_head->reference_count++;
				dst.length++;

				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			}
			free_all(expected_lefts);
			free(*set_var); if (set_var->reference_count == 0) free(set_var);
			free(*expected_conjunct); if (expected_conjunct->reference_count == 0) free(expected_conjunct);
			return true;
		});
	if (!result) return false;
	for (unsigned int i = 0; i < inverse_count; i++) {
		hol_term* new_inverse = hol_term::new_lambda(lambda_variable, inverse[i].root);
		if (new_inverse == nullptr) {
			for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
			free(inverse); return false;
		}
		inverse[i].root = new_inverse;
	}
	return true;
}

template<typename FindHeadFunction>
inline bool invert_factor(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second,
		FindHeadFunction find_head_function)
{
	head_index first_predicate_index, second_predicate_index; no_op apply;
	auto find_array_head = make_array_finder(find_head_function);
	hol_term* first_head = find_head(first, first_predicate_index, find_array_head, apply);
	hol_term* second_head = find_head(second, second_predicate_index, find_head_function, apply);
	if (first_head == nullptr || second_head == nullptr)
		return false;

	if (first_head->type != hol_term_type::ANY && first_head->type != hol_term_type::ANY_RIGHT
	 && first_head->type != hol_term_type::ANY_ARRAY && first_head->type != hol_term_type::AND
	 && first_head->type != hol_term_type::OR)
	{
		return intersect(inverse, inverse_count, flags, first, second);
	}

	hol_term* new_second_head = hol_term::new_any_array(hol_term_type::ANY_ARRAY, second_head,
			make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0));
	if (new_second_head == nullptr)
		return false;
	second_head->reference_count++;

	array<hol_term*> new_heads(4);
	intersect<built_in_predicates>(new_heads, first_head, new_second_head);
	free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);
	if (new_heads.length == 0)
		return false;

	array<hol_term*> intersection(new_heads.length);
	for (hol_term* new_head : new_heads) {
		intersection[intersection.length] = substitute_head<any_node_position::NONE>(second, second_head, new_head);
		if (intersection[intersection.length] == nullptr) {
			free_all(new_heads); free_all(intersection);
			return false;
		}
		intersection.length++;
	}
	free_all(new_heads);

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * intersection.length);
	if (inverse == nullptr) {
		free_all(intersection);
		return false;
	}
	for (unsigned int i = 0; i < intersection.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = intersection[i];
	}
	inverse_count = intersection.length;
	return true;
}

inline bool invert_factor(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_factor(inverse, inverse_count, flags, first, second, find_head<built_in_predicates>);
}

inline bool invert_factor_predicative(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	unsigned int max_variable = 0;
	max_bound_variable(*second, max_variable);

	unsigned int first_lambda_variable = 0;
	if (first->type == hol_term_type::LAMBDA) {
		first_lambda_variable = first->quantifier.variable;
		first = first->quantifier.operand;
	} else {
		return false;
	}

	unsigned int second_lambda_variable = 0;
	if (second->type == hol_term_type::LAMBDA) {
		second_lambda_variable = second->quantifier.variable;
		second = second->quantifier.operand;
	} else {
		return false;
	}

	bool result = invert_factor(inverse, inverse_count, flags, first, second, predicative_head_finder<built_in_predicates>(second_lambda_variable));
	if (!result) return false;
	for (unsigned int i = 0; i < inverse_count; i++) {
		hol_term* new_inverse = hol_term::new_lambda(first_lambda_variable, inverse[i].root);
		if (new_inverse == nullptr) {
			for (unsigned int j = 0; j < inverse_count; j++) free(inverse[j]);
			free(inverse); return false;
		}
		inverse[i].root = new_inverse;
	}
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_select_conjunct_without_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS)
				fprintf(stderr, "invert_select_conjunct_without_head WARNING: Expected an extential quantification.\n");
#endif

			hol_term* operand = second_head->quantifier.operand;
			unsigned int head_variable = second_head->quantifier.variable;

			hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
			if (excluded_quantifier == nullptr)
				return false;
			HOL_ANY.reference_count++;

			hol_term* expected_conjunct = hol_term::new_any(nullptr, &excluded_quantifier, 1);
			if (expected_conjunct == nullptr) {
				free(*excluded_quantifier); free(excluded_quantifier);
				return false;
			}

			hol_term* new_head;
			if (ConjunctIndex >= 0) {
				new_head = hol_term::new_any_array(hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
						make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), operand), make_array_view((hol_term**) nullptr, 0));
				if (new_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				operand->reference_count++;
				expected_conjunct->reference_count += ConjunctIndex;
			} else {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				new_head = hol_term::new_any_array(hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
						make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(operand, make_repeated_array_view(expected_conjunct, index)));
				if (new_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				operand->reference_count++;
				expected_conjunct->reference_count += index;
			}

			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

template<int_fast8_t ConjunctIndex, uint_fast8_t ConjunctCount>
inline bool invert_select_conjuncts_without_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	static_assert(ConjunctCount != 0, "invert_select_conjuncts_without_head ERROR: `ConjunctCount` must be non-zero.");

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
#if !defined(NDEBUG)
			if (second_head->type != hol_term_type::EXISTS || (ConjunctCount > 1 && (second_head->quantifier.operand->type != hol_term_type::AND || second_head->quantifier.operand->array.length != ConjunctCount))) {
				fprintf(stderr, "invert_select_conjuncts_without_head ERROR: Expected an existentially-quantified conjunction of length %u.\n", ConjunctCount);
				return false;
			}
#endif

			hol_term* head_var = hol_term::new_variable(second_head->quantifier.variable);
			if (head_var == nullptr) return false;
			constexpr unsigned int excluded_tree_count = 4;
			hol_term* excluded_trees[excluded_tree_count];
			excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
			excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
			excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
			excluded_trees[3] = hol_term::new_any(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY));
			if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
			if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
			if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
				if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
				if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
				if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
				free(*head_var); free(head_var);
				return false;
			}
			free(*head_var);

			hol_term* expected_conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
			if (expected_conjunct == nullptr) {
				free(*excluded_trees[0]); free(excluded_trees[0]);
				free(*excluded_trees[1]); free(excluded_trees[1]);
				free(*excluded_trees[2]); free(excluded_trees[2]);
				free(*excluded_trees[3]); free(excluded_trees[3]);
				return false;
			}

			hol_term* new_head;
			if (ConjunctIndex >= 0) {
				if (ConjunctCount == 1) {
					new_head = hol_term::new_exists(second_head->quantifier.variable, hol_term::new_any_array(
							hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
							make_appended_array_view(make_repeated_array_view(expected_conjunct, ConjunctIndex), second_head->quantifier.operand),
							make_array_view((hol_term**) nullptr, 0)));
				} else {
					new_head = hol_term::new_exists(second_head->quantifier.variable, hol_term::new_any_array(
							hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
							make_concat_array_view(make_repeated_array_view(expected_conjunct, ConjunctIndex), make_array_view(second_head->quantifier.operand->array.operands, second_head->quantifier.operand->array.length)),
							make_array_view((hol_term**) nullptr, 0)));
				}
				if (new_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				expected_conjunct->reference_count += ConjunctIndex;
			} else {
				unsigned int index = (unsigned int) (-ConjunctIndex) - ConjunctCount;
				if (ConjunctCount == 1) {
					new_head = hol_term::new_exists(second_head->quantifier.variable, hol_term::new_any_array(
							hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0),
							make_prepended_array_view(second_head->quantifier.operand, make_repeated_array_view(expected_conjunct, index))));
				} else {
					new_head = hol_term::new_exists(second_head->quantifier.variable, hol_term::new_any_array(
							hol_term_type::AND, expected_conjunct, make_array_view((hol_term**) nullptr, 0),
							make_array_view((hol_term**) nullptr, 0),
							make_concat_array_view(make_array_view(second_head->quantifier.operand->array.operands, second_head->quantifier.operand->array.length), make_repeated_array_view(expected_conjunct, index))));
				}
				if (new_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				expected_conjunct->reference_count += index;
			}
			if (ConjunctCount == 1) {
				second_head->quantifier.operand->reference_count++;
			} else {
				for (unsigned int i = 0; i < second_head->quantifier.operand->array.length; i++)
					second_head->quantifier.operand->array.operands[i]->reference_count++;
			}
			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

inline hol_term* do_invert_set_predicate_empty(hol_term* second_head)
{
	unsigned int head_variable;
	if (second_head->type == hol_term_type::EXISTS) {
		head_variable = second_head->quantifier.variable;
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*second_head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* head_var = hol_term::new_variable(head_variable);
	if (head_var == nullptr) return (hol_term*) nullptr;
	constexpr unsigned int excluded_tree_count = 4;
	hol_term** excluded_trees = (hol_term**) alloca(sizeof(hol_term) * (excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()));
	excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
	excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
	excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
	excluded_trees[3] = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
	if (excluded_trees[0] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
	if (excluded_trees[1] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
	if (excluded_trees[2] != nullptr) { HOL_ANY.reference_count++; head_var->reference_count++; }
	if (excluded_trees[3] != nullptr) { HOL_ANY.reference_count++; }
	if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || excluded_trees[3] == nullptr) {
		if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
		if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
		if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
		free(*head_var); free(head_var);
		return nullptr;
	}
	free(*head_var);

	for (unsigned int i = 0; i < hol_non_head_constants<built_in_predicates>::count(); i++) {
		excluded_trees[excluded_tree_count + i] = hol_non_head_constants<built_in_predicates>::get_terms()[i];
		excluded_trees[excluded_tree_count + i]->reference_count++;
	}

	hol_term* expected_predicate = hol_term::new_apply(
				hol_term::new_any(nullptr, excluded_trees, excluded_tree_count + hol_non_head_constants<built_in_predicates>::count()), head_var);
	if (expected_predicate == nullptr) {
		for (unsigned int i = 0; i < excluded_tree_count + hol_non_head_constants<built_in_predicates>::count(); i++) {
			free(*excluded_trees[i]); if (excluded_trees[i]->reference_count == 0) free(excluded_trees[i]);
		}
		return nullptr;
	}
	head_var->reference_count++;

	hol_term* empty_predicate = hol_term::new_apply(&HOL_EMPTY, head_var);
	if (empty_predicate == nullptr) {
		free(*expected_predicate); free(expected_predicate);
		return nullptr;
	}
	HOL_EMPTY.reference_count++;
	head_var->reference_count++;

	head_index second_predicate_index;
	find_predicate<built_in_predicates>(head_variable, second_head->quantifier.operand, second_predicate_index);
	hol_term* new_head = apply_predicate(second_head, head_variable, second_predicate_index, empty_predicate, expected_predicate);
	free(*empty_predicate); if (empty_predicate->reference_count == 0) free(empty_predicate);
	free(*expected_predicate); if (expected_predicate->reference_count == 0) free(expected_predicate);
	return new_head;
}

inline bool invert_set_predicate_empty(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	hol_term* second_parent = nullptr;
	hol_term* current = nullptr;
	auto apply = [&second_parent,&current](hol_term* term) {
		second_parent = current;
		current = term;
	};

	head_index second_predicate_index;
	hol_term* second_head = find_head(second, second_predicate_index, find_head<built_in_predicates>, apply);
	if (second_head == nullptr)
		return false;

	hol_term* new_head = apply_array(second_head, second_parent, do_invert_set_predicate_empty);
	if (new_head == nullptr)
		return false;

	hol_term* new_second = substitute_head<any_node_position::NONE>(second, second_head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	if (new_second == nullptr)
		return false;

	array<hol_term*> intersection(4);
	intersect<built_in_predicates>(intersection, first, new_second);
	free(*new_second); if (new_second->reference_count == 0) free(new_second);
	if (intersection.length == 0)
		return false;

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * intersection.length);
	if (inverse == nullptr) {
		free_all(intersection);
		return false;
	}
	for (unsigned int i = 0; i < intersection.length; i++) {
		inverse[i].flags = flags;
		inverse[i].root = intersection[i];
	}
	inverse_count = intersection.length;
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_select_operand(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		make_array_finder(find_root), find_root, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array) {
			hol_term* new_head;
			if (ConjunctIndex >= 0) {
				new_head = hol_term::new_any_array(hol_term_type::ANY_ARRAY, &HOL_ANY, make_array_view((hol_term**) nullptr, 0),
						make_appended_array_view(make_repeated_array_view(&HOL_ANY, ConjunctIndex), second_head), make_array_view((hol_term**) nullptr, 0));
				if (new_head == nullptr)
					return false;
				HOL_ANY.reference_count += 1 + ConjunctIndex;
			} else {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				new_head = hol_term::new_any_array(hol_term_type::ANY_ARRAY, &HOL_ANY, make_array_view((hol_term**) nullptr, 0),
						make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(second_head, make_repeated_array_view(&HOL_ANY, index)));
				if (new_head == nullptr)
					return false;
				HOL_ANY.reference_count += 1 + index;
			}
			second_head->reference_count++;
			dst[dst.length++] = new_head;
			dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count++;
			return true;
		});
}

template<int_fast8_t ConjunctIndex>
inline bool invert_remove_operand(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		make_array_finder(find_root), find_root, no_op(),
		[](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			hol_term* hol_any_ptr = &HOL_ANY;
			if (second_head->type == hol_term_type::ANY_ARRAY) {
				hol_term* new_head;
				if (ConjunctIndex >= 0) {
					if (ConjunctIndex >= second_head->any_array.left.length) {
						new_head = second_head;
						new_head->reference_count++;
					} else {
						new_head = hol_term::new_any_array(second_head->any_array.oper, second_head->any_array.all,
								make_array_view(second_head->any_array.any.operands, second_head->any_array.any.length),
								make_included_array_view(second_head->any_array.left.operands, second_head->any_array.left.length, hol_any_ptr, ConjunctIndex),
								make_array_view(second_head->any_array.right.operands, second_head->any_array.right.length));
						if (new_head == nullptr)
							return false;
						HOL_ANY.reference_count++;
						second_head->any_array.all->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.any.length; i++)
							second_head->any_array.any.operands[i]->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.left.length; i++)
							second_head->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.right.length; i++)
							second_head->any_array.right.operands[i]->reference_count++;
					}
				} else {
					if (-ConjunctIndex > second_head->any_array.right.length) {
						new_head = second_head;
						new_head->reference_count++;
					} else {
						unsigned int index = second_head->any_array.right.length + ConjunctIndex;
						new_head = hol_term::new_any_array(second_head->any_array.oper, second_head->any_array.all,
								make_array_view(second_head->any_array.any.operands, second_head->any_array.any.length),
								make_array_view(second_head->any_array.left.operands, second_head->any_array.left.length),
								make_included_array_view(second_head->any_array.right.operands, second_head->any_array.right.length, hol_any_ptr, index));
						if (new_head == nullptr)
							return false;
						HOL_ANY.reference_count++;
						second_head->any_array.all->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.any.length; i++)
							second_head->any_array.any.operands[i]->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.left.length; i++)
							second_head->any_array.left.operands[i]->reference_count++;
						for (unsigned int i = 0; i < second_head->any_array.right.length; i++)
							second_head->any_array.right.operands[i]->reference_count++;
					}
				}
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			} else if (second_head->type == hol_term_type::AND) {
				unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (second_head->array.length + ConjunctIndex);
				hol_term* new_head = hol_term::new_and(make_included_array_view(second_head->array.operands, second_head->array.length, hol_any_ptr, index));
				if (new_head == nullptr)
					return false;
				for (unsigned int i = 0; i < new_head->array.length; i++)
					new_head->array.operands[i]->reference_count++;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			} else if (second_head->type == hol_term_type::OR) {
				unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (second_head->array.length + ConjunctIndex);
				hol_term* new_head = hol_term::new_or(make_included_array_view(second_head->array.operands, second_head->array.length, hol_any_ptr, index));
				if (new_head == nullptr)
					return false;
				for (unsigned int i = 0; i < new_head->array.length; i++)
					new_head->array.operands[i]->reference_count++;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			} else {
				if (ConjunctIndex < -2 || ConjunctIndex > 1) return false;
				unsigned int index = (ConjunctIndex >= 0) ? ConjunctIndex : (2 + ConjunctIndex);
				hol_term* new_head = (index == 0) ? hol_term::new_and(hol_any_ptr, second_head) : hol_term::new_and(second_head, hol_any_ptr);
				if (new_head == nullptr)
					return false;
				HOL_ANY.reference_count++;
				second_head->reference_count++;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;

				new_head = (index == 0) ? hol_term::new_or(hol_any_ptr, second_head) : hol_term::new_or(second_head, hol_any_ptr);
				if (new_head == nullptr)
					return false;
				HOL_ANY.reference_count++;
				second_head->reference_count++;
				dst[dst.length++] = new_head;
				dst_outer[dst_outer.length++] = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			}
			return true;
		});
}

inline bool invert_apply_to_predicative_set_function(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second,
		unsigned int src_predicate,
		unsigned int dst_predicate)
{
	unsigned int max_variable = 0;
	max_bound_variable(*second, max_variable);

	unsigned int lambda_variable = 0;
	if (second->type == hol_term_type::LAMBDA) {
		lambda_variable = second->quantifier.variable;
		second = second->quantifier.operand;
	} else if (second->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		predicative_head_finder<built_in_predicates>(lambda_variable), predicative_head_finder<built_in_predicates>(lambda_variable), no_op(),
		[src_predicate,dst_predicate](array<hol_term*>& dst, array<hol_term*>& dst_outer, hol_term* first_head, hol_term* second_head, const apply_head_inverter& first_inverter, const apply_head_inverter& second_inverter, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable, bool& any_right_only, bool& could_have_wide_scope, bool is_array)
		{
			if (!apply_to_predicative_set_function(second_head, dst, max_variable, dst_predicate, src_predicate))
				return false;
			if (!dst_outer.ensure_capacity(dst.length)) {
				free_all(dst);
				return false;
			}
			for (unsigned int i = 0; i < dst.length; i++)
				dst_outer[dst_outer.length++] = &HOL_ZERO;
			HOL_ZERO.reference_count += dst.length;
			return true;
		});
}

template<typename Formula>
inline bool invert_add_flag(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	grammatical_flags flags;
	if ((second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE && second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_flag(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	grammatical_flags flags;
	if ((second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_try_remove_flag(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	grammatical_flags flags;
	if ((second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::ANY))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_require_flag(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	grammatical_flags flags;
	if ((second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE && second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_require_no_flag(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	grammatical_flags flags;
	if ((second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_conjunction(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	grammatical_flags flags;
	if ((second.flags.cnj != cnj && second.flags.cnj != grammatical_conjunction::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_conjunction(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	grammatical_flags flags;
	if ((second.flags.cnj != grammatical_conjunction::NONE && second.flags.cnj != grammatical_conjunction::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, cnj))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_require_no_conjunction(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	grammatical_flags flags;
	if (second.flags.cnj == cnj
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::ANY))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_correlator(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator corr)
{
	grammatical_flags flags;
	if ((second.flags.corr != corr && second.flags.corr != correlator::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_correlator(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator corr)
{
	grammatical_flags flags;
	if ((second.flags.corr != correlator::NONE && second.flags.corr != correlator::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.corr, first.flags.corr, corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_correlated_by(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator correlated_by)
{
	grammatical_flags flags;
	if ((second.flags.correlated_by != correlated_by && second.flags.correlated_by != correlator::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_coordination(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		coordination coord)
{
	grammatical_flags flags;
	if ((second.flags.coord != coord && second.flags.coord != coordination::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.coord, first.flags.coord, coordination::NONE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_coordination(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		coordination coord)
{
	grammatical_flags flags;
	if ((second.flags.coord != coordination::NONE && second.flags.coord != coordination::ANY)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_apply_auxiliary(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		const auxiliary_flag expected_aux,
		auxiliary_flag dst_aux)
{
	grammatical_flags flags;
	if (!has_intersection(dst_aux, second.flags.aux)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.aux, expected_aux, first.flags.aux))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_mood(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_mood mood)
{
	grammatical_flags flags;
	if (!has_intersection(second.flags.mood, mood)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::INDICATIVE))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_mood(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_mood mood)
{
	grammatical_flags flags;
	if (!has_intersection(second.flags.mood, grammatical_mood::INDICATIVE)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::ANY))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
inline bool invert_try_remove_mood(
		flagged_logical_form<Formula>*& inverse,
		unsigned int& inverse_count,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	grammatical_flags flags;
	if (!has_intersection(second.flags.mood, grammatical_mood::INDICATIVE)
	 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
	 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
	 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::ANY))
		return false;
	flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
		&& intersect(inverse, inverse_count, flags, first.root, second.root);
}

template<typename Formula>
bool init_array(flagged_logical_form<Formula>*& array, unsigned int length)
{
	array = (flagged_logical_form<Formula>*) malloc(sizeof(flagged_logical_form<Formula>) * length);
	if (array == nullptr) {
		fprintf(stderr, "init_array ERROR: Out of memory.\n");
		return false;
	}
	return true;
}

template<typename Formula>
bool invert(
	flagged_logical_form<Formula>*& inverse,
	unsigned int& inverse_count,
	typename flagged_logical_form<Formula>::function function,
	const flagged_logical_form<Formula>& first,
	const flagged_logical_form<Formula>& second)
{
	typedef typename flagged_logical_form<Formula>::function_type function_type;

	uint_fast8_t i;
	grammatical_flags flags;
	switch (function.type) {
	case function_type::EMPTY:
		inverse_count = 1; if (!init_array(inverse, inverse_count)) return false;
		*inverse = first;
		return true;
	case function_type::IDENTITY:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NO_INVERSE:
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE:
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN:
	case function_type::REQUIRE_LEFT_PREDICATE_EXIST:
	case function_type::REQUIRE_NO_LEFT_PREDICATE_EXIST:
	case function_type::REQUIRE_NO_FUTURE:
	case function_type::REQUIRE_NO_PERFECT:
	case function_type::REQUIRE_NO_PROGRESSIVE:
	case function_type::REQUIRE_NO_EMPTY_REF:
	case function_type::REQUIRE_NO_REQ_AUX:
	case function_type::REQUIRE_NO_REQ_NO_AUX:
	case function_type::REQUIRE_TO_INFINITIVE:
	case function_type::REQUIRE_NO_TO_INFINITIVE:
	case function_type::REQUIRE_NO_SUBJUNCTIVE:
	case function_type::REQUIRE_AUX_OR_SUBJUNCTIVE_OR_INFINITIVE_OR_TO_INFINITIVE:
	case function_type::REQUIRE_PAST_PARTICIPLE:
	case function_type::REQUIRE_PRESENT_PARTICIPLE:
	case function_type::REQUIRE_PREDICATIVE_UNIVERSAL:
	case function_type::REQUIRE_PREDICATIVE_EXISTENTIAL:
	case function_type::REQUIRE_CONSTANT_IN_SET:
	case function_type::REQUIRE_NO_CONSTANT_IN_SET:
	case function_type::REQUIRE_SINGLETON:
	case function_type::REQUIRE_LEFT_ARG1:
	case function_type::REQUIRE_LAMBDA:
	case function_type::REQUIRE_NO_LAMBDA:
		/* the forward application already ensures that `second` satisfies this requirement */
		if (!intersect(flags, first.flags, second.flags)) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_CONJUNCTION:
	case function_type::REQUIRE_BINARY_CONJUNCTION:
	case function_type::REQUIRE_DISJUNCTION:
	case function_type::REQUIRE_NEGATIVE_CONJUNCTION:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return intersect_with_head(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_CONJUNCT_AND_NEGATION:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<0, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_LEFT_CONJUNCT:
	case function_type::REMOVE_LEFT_PREDICATE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_LEFT_CONJUNCT_AND_NEGATION:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<0, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_SECOND_LEFT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_SECOND_LEFT_SET_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_set_conjunct<1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_SECOND_LEFT_SET_CONJUNCT_ROOT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_set_conjunct<1, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_SECOND_LEFT_SET_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_set_conjunct<1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_CONJUNCT_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct_in_set<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_LEFT_CONJUNCT_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct_in_set<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_CONJUNCT_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct_in_set<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_RIGHT_CONJUNCT_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct_in_set<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_SUBSET_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_subset_in_set<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_INVERSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_inverse(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG1, false>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG2, false>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG3, false>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG1_OF_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG1_OF, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG2_OF_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG2_OF, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG3_OF_WITHOUT_HEAD:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head<-1, (unsigned int) built_in_predicates::ARG3_OF, true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG2>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative<-1, (unsigned int) built_in_predicates::ARG3>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_FUTURE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_apply_tense_predicate(inverse, inverse_count, flags, first.root, second.root, FUTURE_PREDICATES, PRESENT_PREDICATES);
	case function_type::REMOVE_PERFECT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_apply_tense_predicate(inverse, inverse_count, flags, first.root, second.root, PERFECT_PREDICATES, NON_PERFECT_PREDICATES);
	case function_type::REMOVE_PROGRESSIVE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_apply_tense_predicate(inverse, inverse_count, flags, first.root, second.root, PROGRESSIVE_PREDICATES, NON_PROGRESSIVE_PREDICATES);
	case function_type::REMOVE_NOT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_not(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_PREDICATIVE_NOT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_predicative_not(inverse, inverse_count, flags, first.root, second.root);
	case function_type::PREDICATE_ONLY:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate_only(inverse, inverse_count, flags, first.root, second.root);
	case function_type::PREDICATE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate(inverse, inverse_count, flags, first.root, second.root);
	case function_type::PREDICATE_AND_TENSE:
	case function_type::EMPTY_AND_TENSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate_and_tense(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_PREDICATE_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_predicate_in_set(inverse, inverse_count, flags, first.root, second.root);
	case function_type::MARK_WIDE_SCOPE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_mark_wide_scope(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_WIDE_SCOPE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_require_narrow_or_wide_scope<true>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NARROW_SCOPE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_require_narrow_or_wide_scope<false>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_WIDE_SCOPE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_wide_scope(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SIZE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_size(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SET_SIZE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_set_size(inverse, inverse_count, flags, first.root, second.root);
	case function_type::FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_factor(inverse, inverse_count, flags, first.root, second.root);
	case function_type::FACTOR_PREDICATIVE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_factor_predicative(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_PREDICATE_AND_TENSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjuncts_without_head<0, 2>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SET_PREDICATE_EMPTY:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_set_predicate_empty(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_OPERAND:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_operand<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_LEFT_OPERAND:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_operand<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REPLACE_PREDICATIVE_SUBSET_WITH_EQUALITY:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_apply_to_predicative_set_function(inverse, inverse_count, flags, first.root, second.root,
				(unsigned int) built_in_predicates::SUBSET, (unsigned int) built_in_predicates::EQUALS);
	case function_type::ADD_SINGULAR:
		if ((second.flags.index_number != grammatical_num::SINGULAR && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_PLURAL:
		if ((second.flags.index_number != grammatical_num::PLURAL && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_SINGULAR:
		if ((second.flags.index_number != grammatical_num::SINGULAR && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::SINGULAR)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_PLURAL:
		if ((second.flags.index_number != grammatical_num::PLURAL && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::PLURAL)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::TRY_REMOVE_NUMBER:
		if ((second.flags.index_number != grammatical_num::NONE && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		flags.index_number = first.flags.index_number;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_CONCORD_SINGULAR:
		if ((second.flags.concord_number != grammatical_num::SINGULAR && second.flags.concord_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, grammatical_num::NONE)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_CONCORD_PLURAL:
		if ((second.flags.concord_number != grammatical_num::PLURAL && second.flags.concord_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, grammatical_num::NONE)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_THAT:
		return invert_add_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::THAT);
	case function_type::REMOVE_THAT:
		return invert_remove_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::THAT);
	case function_type::REQUIRE_NO_THAT:
		return invert_require_no_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::THAT);
	case function_type::ADD_WHETHER:
		return invert_add_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::WHETHER);
	case function_type::REMOVE_WHETHER:
		return invert_remove_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::WHETHER);
	case function_type::ADD_IF:
		return invert_add_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::IF);
	case function_type::REMOVE_IF:
		return invert_remove_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::IF);
	case function_type::ADD_BECAUSE:
		return invert_add_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::BECAUSE);
	case function_type::REMOVE_BECAUSE:
		return invert_remove_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::BECAUSE);
	case function_type::ADD_FOR:
		return invert_add_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::FOR);
	case function_type::REMOVE_FOR:
		return invert_remove_conjunction(inverse, inverse_count, first, second, grammatical_conjunction::FOR);
	case function_type::REQUIRE_NO_CONJUNCTION:
		if ((second.flags.cnj != grammatical_conjunction::NONE && second.flags.cnj != grammatical_conjunction::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_IS_ADJUNCT:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::TRY_REMOVE_IS_ADJUNCT:
		return invert_try_remove_flag(inverse, inverse_count, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::REQUIRE_NOT_ADJUNCT:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::ADD_NULLABLE_SUBJECT:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::REMOVE_NULLABLE_SUBJECT:
		return invert_remove_flag(inverse, inverse_count, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::TRY_REMOVE_NULLABLE_SUBJECT:
		return invert_try_remove_flag(inverse, inverse_count, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::ADD_SUBORDINATE:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::SUBORDINATE);
	case function_type::REMOVE_SUBORDINATE:
		return invert_remove_flag(inverse, inverse_count, first, second, grammatical_flag::SUBORDINATE);
	case function_type::TRY_REMOVE_SUBORDINATE:
		return invert_try_remove_flag(inverse, inverse_count, first, second, grammatical_flag::SUBORDINATE);
	case function_type::REQUIRE_NO_SUBORDINATE:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::SUBORDINATE);
	case function_type::ADD_PREPOSITION:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::PREPOSITION);
	case function_type::REQUIRE_PREPOSITION:
		return invert_require_flag(inverse, inverse_count, first, second, grammatical_flag::PREPOSITION);
	case function_type::REQUIRE_NO_PREPOSITION:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::PREPOSITION);
	case function_type::ADD_PARTICLE:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::PARTICLE);
	case function_type::ADD_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NONE_OR_REQ_AUX, auxiliary_flag::AUX);
	case function_type::TRY_REMOVE_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NO_REQ_AUX, auxiliary_flag::NONE);
	case function_type::ADD_REQ_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NONE_OR_AUX, auxiliary_flag::REQ_AUX);
	case function_type::TRY_ADD_REQ_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NO_REQ_NO_AUX, auxiliary_flag::REQ_AUX);
	case function_type::TRY_REMOVE_REQ_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NO_REQ_NO_AUX, auxiliary_flag::NONE);
	case function_type::ADD_REQ_NO_AUX:
		return invert_apply_auxiliary(inverse, inverse_count, first, second, auxiliary_flag::NONE_OR_AUX, auxiliary_flag::REQ_NO_AUX);
	case function_type::ADD_INFINITIVE:
		return invert_add_mood(inverse, inverse_count, first, second, grammatical_mood::BARE_INFINITIVE);
	case function_type::ADD_TO_INFINITIVE:
		return invert_add_mood(inverse, inverse_count, first, second, grammatical_mood::TO_INFINITIVE);
	case function_type::REMOVE_TO_INFINITIVE:
		return invert_remove_mood(inverse, inverse_count, first, second, grammatical_mood::TO_INFINITIVE);
	case function_type::ADD_SUBJUNCTIVE:
		return invert_add_mood(inverse, inverse_count, first, second, grammatical_mood::SUBJUNCTIVE);
	case function_type::ADD_BOTH:
		return invert_add_correlator(inverse, inverse_count, first, second, correlator::BOTH);
	case function_type::ADD_EITHER:
		return invert_add_correlator(inverse, inverse_count, first, second, correlator::EITHER);
	case function_type::ADD_NEITHER:
		return invert_add_correlator(inverse, inverse_count, first, second, correlator::NEITHER);
	case function_type::REMOVE_BOTH:
		return invert_remove_correlator(inverse, inverse_count, first, second, correlator::BOTH);
	case function_type::REMOVE_EITHER:
		return invert_remove_correlator(inverse, inverse_count, first, second, correlator::EITHER);
	case function_type::REMOVE_NEITHER:
		return invert_remove_correlator(inverse, inverse_count, first, second, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATOR:
		if ((second.flags.corr != correlator::NONE && second.flags.corr != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.corr, first.flags.corr, correlator::ANY))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NO_CORRELATOR:
		if ((second.flags.corr != correlator::NONE && second.flags.corr != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_CORRELATED_BY_BOTH:
		return invert_add_correlated_by(inverse, inverse_count, first, second, correlator::BOTH);
	case function_type::ADD_CORRELATED_BY_EITHER:
		return invert_add_correlated_by(inverse, inverse_count, first, second, correlator::EITHER);
	case function_type::ADD_CORRELATED_BY_NEITHER:
		return invert_add_correlated_by(inverse, inverse_count, first, second, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATED:
		if ((second.flags.correlated_by != correlator::NONE && second.flags.correlated_by != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::ANY)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NOT_CORRELATED:
		if ((second.flags.correlated_by != correlator::NONE && second.flags.correlated_by != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.comp, first.flags.comp, second.flags.comp)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE)
		 || !intersect(flags.coord, first.flags.coord, second.flags.coord)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.is_first_token_capital = (first.flags.is_first_token_capital || second.flags.is_first_token_capital);
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_PAST_PARTICIPLE:
		return invert_add_mood(inverse, inverse_count, first, second, grammatical_mood::PAST_PARTICIPLE);
	case function_type::ADD_PRESENT_PARTICIPLE:
		return invert_add_mood(inverse, inverse_count, first, second, grammatical_mood::PRESENT_PARTICIPLE);
	case function_type::TRY_REMOVE_PARTICIPLE:
		return invert_try_remove_mood(inverse, inverse_count, first, second);
	case function_type::ADD_NEGATIVE:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::NEGATIVE);
	case function_type::REQUIRE_NEGATIVE:
		return invert_require_flag(inverse, inverse_count, first, second, grammatical_flag::NEGATIVE);
	case function_type::ADD_ADV:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::ADV);
	case function_type::REMOVE_ADV:
		return invert_remove_flag(inverse, inverse_count, first, second, grammatical_flag::ADV);
	case function_type::TRY_REMOVE_ADV:
		return invert_try_remove_flag(inverse, inverse_count, first, second, grammatical_flag::ADV);
	case function_type::REQUIRE_NO_ADV:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::ADV);
	case function_type::ADD_TION:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::TION);
	case function_type::ADD_LY:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::LY);
	case function_type::ADD_GENITIVE:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::GENITIVE);
	case function_type::TRY_REMOVE_GENITIVE:
		return invert_try_remove_flag(inverse, inverse_count, first, second, grammatical_flag::GENITIVE);
	case function_type::REQUIRE_NO_GENITIVE:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::GENITIVE);
	case function_type::ADD_COMMA:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::COMMA);
	case function_type::REMOVE_COMMA:
		return invert_remove_flag(inverse, inverse_count, first, second, grammatical_flag::COMMA);
	case function_type::REQUIRE_NO_COMMA:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::COMMA);
	case function_type::ADD_AND:
		return invert_add_coordination(inverse, inverse_count, first, second, coordination::AND);
	case function_type::ADD_OR:
		return invert_add_coordination(inverse, inverse_count, first, second, coordination::OR);
	case function_type::ADD_NOR:
		return invert_add_coordination(inverse, inverse_count, first, second, coordination::NOR);
	case function_type::REMOVE_AND:
		return invert_remove_coordination(inverse, inverse_count, first, second, coordination::AND);
	case function_type::REMOVE_OR:
		return invert_remove_coordination(inverse, inverse_count, first, second, coordination::OR);
	case function_type::REMOVE_NOR:
		return invert_remove_coordination(inverse, inverse_count, first, second, coordination::NOR);
	case function_type::REMOVE_COORD:
		return invert_remove_coordination(inverse, inverse_count, first, second, coordination::NOT_NONE);
	}
	fprintf(stderr, "invert ERROR: Unrecognized transformation function.\n");
	return false;
}

template<typename Formula>
inline bool is_subset(
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	if (!is_subset(first.flags, second.flags))
		return false;

/* TODO: for debugging; remove this */
//print("first.root:  ", stderr); print(*first.root, stderr, *debug_terminal_printer); print('\n', stderr);
//print("second.root: ", stderr); print(*second.root, stderr, *debug_terminal_printer); print('\n', stderr);
	array<pair<hol_term*, variable_map>> intersection(2);
	intersect<built_in_predicates>(intersection, first.root, second.root);
	if (intersection.length == 0)
		return false;
	free_all(intersection);
	return true;
}

template<typename Formula>
static void is_separable(
		const transformation<flagged_logical_form<Formula>>* functions,
		unsigned int rule_length, bool* separable)
{
	for (unsigned int i = 0; i < rule_length; i++)
		separable[i] = false;
}

inline bool copy_array(
		const unsigned int* src, unsigned int src_length,
		unsigned int*& dst, unsigned int& dst_length)
{
	if (src_length == 0) {
		dst_length = 0;
		return true;
	}
	dst = (unsigned int*) malloc(sizeof(unsigned int) * src_length);
	if (dst == nullptr) {
		fprintf(stderr, "copy_array ERROR: Out of memory.\n");
		return false;
	}
	memcpy(dst, src, sizeof(unsigned int) * src_length);
	dst_length = src_length;
	return true;
}

bool get_constant(hol_term* src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	hol_term* term = hol_term::new_any_constant_except();
	if (term == nullptr) return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "get_constant ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	if (intersection[0]->type == hol_term_type::CONSTANT) {
		value = intersection[0]->constant;
		excluded_count = 0;
	} else if (intersection[0]->type == hol_term_type::ANY_CONSTANT) {
		value = UNION_NODE;
		if (!copy_array(intersection[0]->any_constant.constants, intersection[0]->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	} else {
#if !defined(NDEBUG)
		if (intersection[0]->type != hol_term_type::ANY_CONSTANT_EXCEPT)
			fprintf(stderr, "get_constant WARNING: Unexpected formula type in interseciton.\n");
#endif
		value = IMPLICIT_NODE;
		if (intersection[0]->any_constant.length > 0 && !copy_array(intersection[0]->any_constant.constants, intersection[0]->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
		excluded_count = intersection[0]->any_constant.length;
	}
	for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
	return true;
}

bool get_predicate(hol_term* src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS,
			hol_term::new_and(hol_term::new_apply(hol_term::new_any_constant_except(), &HOL_ANY), &HOL_ANY));
	if (term == nullptr) return false;
	HOL_ANY.reference_count += 2;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "get_predicate ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	hol_term* operand;
	if (intersection[0]->type == hol_term_type::EXISTS) {
		operand = intersection[0]->quantifier.operand;
	} else {
		operand = intersection[0]->any_quantifier.operand;
	}

	hol_term* predicate;
	if (operand->type == hol_term_type::AND) {
		predicate = operand->array.operands[0]->binary.left;
	} else if (operand->type == hol_term_type::ANY_ARRAY) {
		predicate = operand->any_array.left.operands[0]->binary.left;
	} else {
		predicate = operand->binary.left;
	}

	if (predicate->type == hol_term_type::CONSTANT) {
		value = predicate->constant;
		excluded_count = 0;
	} else if (predicate->type == hol_term_type::ANY_CONSTANT) {
		value = UNION_NODE;
		if (!copy_array(predicate->any_constant.constants, predicate->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	} else {
#if !defined(NDEBUG)
		if (predicate->type != hol_term_type::ANY_CONSTANT_EXCEPT)
			fprintf(stderr, "get_predicate WARNING: Unexpected formula type in interseciton.\n");
#endif
		value = IMPLICIT_NODE;
		if (!copy_array(predicate->any_constant.constants, predicate->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	}
	for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
	return true;
}

bool get_predicate_only(hol_term* src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	hol_term* head = src;
	if ((src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) && src->any.included != nullptr && src->any.included->type == hol_term_type::EXISTS)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else if (head->type == hol_term_type::ANY_ARRAY) {
		if (head->any_array.any.length > 1 || head->any_array.left.length > 1 || head->any_array.right.length > 1)
			return false;
		if (head->any_array.left.length == 1 && head->any_array.left.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.left.operands[0]->quantifier.variable;
		} else if (head->any_array.right.length == 1 && head->any_array.right.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.right.operands[0]->quantifier.variable;
		} else if (head->any_array.any.length == 1 && head->any_array.any.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.any.operands[0]->quantifier.variable;
		} else if (head->any_array.all->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.all->quantifier.variable;
		} else {
			unsigned int max_variable = 0;
			max_bound_variable(*head, max_variable);
			head_variable = ++max_variable;
		}
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* term = hol_term::new_exists(head_variable,
			hol_term::new_apply(hol_term::new_any_constant_except(), hol_term::new_variable(head_variable)));
	if (term == nullptr)
		return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "get_predicate_only ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	hol_term* operand;
	if (intersection[0]->type == hol_term_type::EXISTS) {
		operand = intersection[0]->quantifier.operand;
	} else {
		operand = intersection[0]->any_quantifier.operand;
	}

	hol_term* predicate = operand->binary.left;
	if (predicate->type == hol_term_type::CONSTANT) {
		value = predicate->constant;
		excluded_count = 0;
	} else if (predicate->type == hol_term_type::ANY_CONSTANT) {
		value = UNION_NODE;
		if (!copy_array(predicate->any_constant.constants, predicate->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	} else {
#if !defined(NDEBUG)
		if (predicate->type != hol_term_type::ANY_CONSTANT_EXCEPT)
			fprintf(stderr, "get_predicate_only WARNING: Unexpected formula type in interseciton.\n");
#endif
		value = IMPLICIT_NODE;
		if (!copy_array(predicate->any_constant.constants, predicate->any_constant.length, excluded, excluded_count)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}
	}
	for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
	return true;
}

bool get_set_definition(hol_term* src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else {
		return false;
	}

	head_index predicate_index; no_op apply;
	auto head_finder = predicative_head_finder<built_in_predicates>(lambda_variable);
	hol_term* head = find_head(src, predicate_index, head_finder, apply);
	if (head == nullptr)
		return false;

	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int set_variable;
	if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		value = IMPLICIT_NODE;
		excluded_count = 0;
		return true;
	} else {
#if !defined(NDEBUG)
		if (head->type != hol_term_type::EXISTS) {
			fprintf(stderr, "get_set_definition ERROR: Expected existential quantification of set.\n");
			return false;
		}
#endif
		set_variable = head->quantifier.variable;
	}

	hol_term* expected_left = hol_term::new_equals(hol_term::new_variable(set_variable), &HOL_ANY);
	if (expected_left == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* expected_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
			make_array_view((hol_term**) nullptr, 0), make_array_view(&expected_left, 1), make_array_view((hol_term**) nullptr, 0)));
	if (expected_head == nullptr) {
		free(*expected_left); free(expected_left);
		return (hol_term*) nullptr;
	}
	HOL_ANY.reference_count++;

	bool is_complement = false;
	array<unsigned int> allowed_set_definitions(4);
	if (has_intersection<built_in_predicates>(head, expected_head))
		allowed_set_definitions[allowed_set_definitions.length++] = (unsigned int) built_in_predicates::EQUALS;
	free(*expected_head); free(expected_head);

	expected_left = hol_term::new_apply(hol_term::new_any_constant_except(), hol_term::new_variable(set_variable), &HOL_ANY);
	if (expected_left == nullptr)
		return false;
	HOL_ANY.reference_count++;

	expected_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
			make_array_view((hol_term**) nullptr, 0), make_array_view(&expected_left, 1), make_array_view((hol_term**) nullptr, 0)));
	if (expected_head == nullptr) {
		free(*expected_left); free(expected_left);
		return (hol_term*) nullptr;
	}
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(4);
	intersect<built_in_predicates>(intersection, head, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	for (hol_term* new_head : intersection) {
		hol_term* predicate;
		hol_term* operand = new_head->quantifier.operand;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			predicate = operand->any_array.left.operands[0]->ternary.first;
		} else {
			predicate = operand->array.operands[0]->ternary.first;
		}

		if (predicate->type == hol_term_type::CONSTANT) {
			if (is_complement) {
				unsigned int index = allowed_set_definitions.index_of(predicate->constant);
				if (index < allowed_set_definitions.length) allowed_set_definitions.remove(index);
			} else {
				if (!allowed_set_definitions.contains(predicate->constant)) {
					if (!allowed_set_definitions.add(predicate->constant)) {
						free_all(intersection);
						return false;
					}
					insertion_sort(allowed_set_definitions);
				}
			}
		} else if (predicate->type == hol_term_type::ANY_CONSTANT) {
			if (is_complement) {
				set_subtract(allowed_set_definitions.data, allowed_set_definitions.length, predicate->any_constant.constants, predicate->any_constant.length);
			} else {
				array<unsigned int> new_allowed_set_definitions(allowed_set_definitions.length + predicate->any_constant.length);
				set_union(new_allowed_set_definitions.data, new_allowed_set_definitions.length,
						allowed_set_definitions.data, allowed_set_definitions.length,
						predicate->any_constant.constants, predicate->any_constant.length);
				swap(new_allowed_set_definitions, allowed_set_definitions);
			}
		} else if (predicate->type == hol_term_type::ANY_CONSTANT_EXCEPT) {
			if (is_complement) {
				set_intersect(allowed_set_definitions, predicate->any_constant.constants, predicate->any_constant.length);
			} else {
				array<unsigned int> new_allowed_set_definitions(max(1u, predicate->any_constant.length));
				set_subtract(new_allowed_set_definitions.data, new_allowed_set_definitions.length,
						predicate->any_constant.constants, predicate->any_constant.length,
						allowed_set_definitions.data, allowed_set_definitions.length);
				swap(new_allowed_set_definitions, allowed_set_definitions);
				is_complement = true;
			}
		}
	}
	free_all(intersection);

	if (allowed_set_definitions.length == 0) {
		if (!is_complement)
			return false;
		value = IMPLICIT_NODE;
		excluded_count = 0;
		return true;
	} else if (allowed_set_definitions.length == 1 && !is_complement) {
		value = allowed_set_definitions[0];
		excluded_count = 0;
		return true;
	} else {
		if (is_complement) value = IMPLICIT_NODE;
		else value = UNION_NODE;
		return copy_array(allowed_set_definitions.data, allowed_set_definitions.length, excluded, excluded_count);
	}
}

template<typename BuiltInPredicates>
struct arg_finder {
	unsigned int head_variable;

	arg_finder(unsigned int head_variable) : head_variable(head_variable) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		if (src->type == hol_term_type::EQUALS
		 && src->binary.left->type == hol_term_type::UNARY_APPLICATION
		 && src->binary.left->binary.right->type == hol_term_type::VARIABLE
		 && src->binary.left->binary.right->variable == head_variable
		 && (src->binary.left->binary.left->type == hol_term_type::CONSTANT
		  || src->binary.left->binary.left->type == hol_term_type::ANY_CONSTANT))
		{
			head = src->binary.left->binary.left;
		} else {
			head = nullptr;
		}
	}
};

template<int_fast8_t ConjunctIndex>
bool get_arg(hol_term* src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	head_index predicate_index; no_op apply;
	hol_term* head = find_head(src, predicate_index, find_head<built_in_predicates>, apply);
	if (head == nullptr)
		return false;

	unsigned int head_variable;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;
	while (head->type == hol_term_type::NOT)
		head = head->unary.operand;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
		value = UNION_NODE;
		excluded = (unsigned int*) malloc(sizeof(unsigned int) * 4);
		if (excluded == nullptr) {
			fprintf(stderr, "get_arg ERROR: Out of memory.\n");
			return false;
		}
		excluded[0] = (unsigned int) built_in_predicates::UNKNOWN;
		excluded[1] = (unsigned int) built_in_predicates::ARG1;
		excluded[2] = (unsigned int) built_in_predicates::ARG2;
		excluded[3] = (unsigned int) built_in_predicates::ARG3;
		excluded_count = 4;
		return true;
	} else {
		return false;
	}

	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
	if (excluded_quantifier == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* excluded_tree = hol_term::new_any(nullptr, &excluded_quantifier, 1);
	if (excluded_tree == nullptr) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}

	hol_term* expected_conjunct = hol_term::new_any_right(hol_term::new_equals(
			hol_term::new_apply(hol_term::new_any_constant(
				(unsigned int) built_in_predicates::ARG1,
				(unsigned int) built_in_predicates::ARG2,
				(unsigned int) built_in_predicates::ARG3),
			hol_term::new_variable(head_variable)), excluded_tree));
	if (expected_conjunct == nullptr) {
		free(*excluded_tree); free(excluded_tree);
		return false;
	}
	excluded_tree->reference_count++;

	hol_term* expected_head;
	if (ConjunctIndex >= 0) {
		expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(
				hol_term_type::AND, excluded_tree, make_array_view((hol_term**) nullptr, 0),
				make_appended_array_view(make_repeated_array_view(excluded_tree, ConjunctIndex), expected_conjunct),
				make_array_view((hol_term**) nullptr, 0)));
		if (expected_head == nullptr) {
			free(*expected_conjunct); free(expected_conjunct);
			free(*excluded_tree); free(excluded_tree);
			return false;
		}
		excluded_tree->reference_count += ConjunctIndex;
	} else {
		unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
		expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(
				hol_term_type::AND, excluded_tree, make_array_view((hol_term**) nullptr, 0),
				make_array_view((hol_term**) nullptr, 0),
				make_prepended_array_view(expected_conjunct, make_repeated_array_view(excluded_tree, index))));
		if (expected_head == nullptr) {
			free(*expected_conjunct); free(expected_conjunct);
			free(*excluded_tree); free(excluded_tree);
			return false;
		}
		excluded_tree->reference_count += index;
	}

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, head, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (intersection.length == 0) {
		value = (unsigned int) built_in_predicates::UNKNOWN;
		return true;
	} else if (intersection.length != 1) {
		fprintf(stderr, "get_arg ERROR: Intersection is not unique.\n");
		free_all(intersection); return false;
	}

	bool is_subset = (intersection[0] == head || *intersection[0] == *head);

	hol_term* conjunct;
	hol_term* operand = intersection[0]->quantifier.operand;
	if (operand->type == hol_term_type::AND) {
		if (ConjunctIndex >= 0)
			conjunct = operand->array.operands[ConjunctIndex];
		else conjunct = operand->array.operands[operand->array.length + ConjunctIndex];
	} else {
		if (ConjunctIndex >= 0)
			conjunct = operand->any_array.left.operands[ConjunctIndex];
		else conjunct = operand->any_array.right.operands[operand->any_array.right.length + ConjunctIndex];
	}

	arg_finder<built_in_predicates> find_arg(head_variable);
	hol_term* arg = find_head(conjunct, predicate_index, find_arg, apply);

	if (arg->type == hol_term_type::CONSTANT) {
		if (!is_subset) {
			value = UNION_NODE;
			excluded = (unsigned int*) malloc(sizeof(unsigned int) * 2);
			if (excluded == nullptr) {
				fprintf(stderr, "get_arg ERROR: Out of memory.\n");
				free_all(intersection); return false;
			}
			excluded[0] = (unsigned int) built_in_predicates::UNKNOWN;
			excluded[1] = arg->constant;
			excluded_count = 2;
		} else {
			value = arg->constant;
			excluded_count = 0;
		}
	} else if (arg->type == hol_term_type::ANY_CONSTANT) {
		if (!is_subset) {
			value = UNION_NODE;
			excluded = (unsigned int*) malloc(sizeof(unsigned int) * (arg->any_constant.length + 1));
			if (excluded == nullptr) {
				fprintf(stderr, "get_arg ERROR: Out of memory.\n");
				free_all(intersection); return false;
			}
			excluded[0] = (unsigned int) built_in_predicates::UNKNOWN;
			for (unsigned int i = 0; i < arg->any_constant.length; i++)
				excluded[1 + i] = arg->any_constant.constants[i];
			excluded_count = arg->any_constant.length + 1;
		} else {
			value = UNION_NODE;
			if (!copy_array(arg->any_constant.constants, arg->any_constant.length, excluded, excluded_count)) {
				free_all(intersection);
				return false;
			}
		}
	}
	free_all(intersection);
	return true;
}

template<typename Formula>
bool get_feature(
		typename flagged_logical_form<Formula>::feature feature,
		const flagged_logical_form<Formula>& src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	switch (feature) {
	case feature_type::CONSTANT:
		return get_constant(src.root, value, excluded, excluded_count);
	case feature_type::PREDICATE:
		return get_predicate(src.root, value, excluded, excluded_count);
	case feature_type::PREDICATE_ONLY:
		return get_predicate_only(src.root, value, excluded, excluded_count);
	case feature_type::SET_DEFINITION:
		return get_set_definition(src.root, value, excluded, excluded_count);
	case feature_type::LEFT_ARG:
		return get_arg<0>(src.root, value, excluded, excluded_count);
	case feature_type::EMPTY: break;
	}
	fprintf(stderr, "get_feature ERROR: Unrecognized semantic feature.\n");
	return false;
}

bool set_constant(hol_term* src, hol_term*& dst, unsigned int value)
{
	hol_term* term = hol_term::new_constant(value);
	if (term == nullptr) return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "set_constant ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

bool set_predicate(hol_term* src, hol_term*& dst, unsigned int value)
{
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS,
			hol_term::new_and(hol_term::new_apply(hol_term::new_constant(value), &HOL_ANY), &HOL_ANY));
	if (term == nullptr) return false;
	HOL_ANY.reference_count += 2;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "set_predicate ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

bool set_predicate_only(hol_term* src, hol_term*& dst, unsigned int value)
{
	hol_term* head = src;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else if (head->type == hol_term_type::ANY_ARRAY) {
		if (head->any_array.any.length > 1 || head->any_array.left.length > 1 || head->any_array.right.length > 1)
			return false;
		if (head->any_array.left.length == 1 && head->any_array.left.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.left.operands[0]->quantifier.variable;
		} else if (head->any_array.right.length == 1 && head->any_array.right.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.right.operands[0]->quantifier.variable;
		} else if (head->any_array.any.length == 1 && head->any_array.any.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.any.operands[0]->quantifier.variable;
		} else if (head->any_array.all->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.all->quantifier.variable;
		} else {
			unsigned int max_variable = 0;
			max_bound_variable(*head, max_variable);
			head_variable = ++max_variable;
		}
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* term = hol_term::new_exists(head_variable,
			hol_term::new_apply(hol_term::new_constant(value), hol_term::new_variable(head_variable)));
	if (term == nullptr)
		return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "set_predicate_only ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

bool set_set_definition(hol_term* src, hol_term*& dst, unsigned int value)
{
	unsigned int max_variable = 0;
	max_bound_variable(*src, max_variable);

	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
	} else if (src->type == hol_term_type::ANY) {
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false, remove_negations;
	bool result = apply_head(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_negations, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable),
			[&max_variable, value](hol_term* head, unsigned int head_variable, head_index predicate_index, bool is_array, bool& remove_wide_scope_marker, array<hol_term*>& siblings)
			{
				if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
					head = head->any.included;

				unsigned int set_variable;
				if (head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) {
					set_variable = ++max_variable;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "set_set_definition ERROR: Expected existential quantification of set.\n");
						return (hol_term*) nullptr;
					}
#endif
					set_variable = head->quantifier.variable;
				}

				hol_term* expected_left;
				if (value == (unsigned int) built_in_predicates::EQUALS) {
					expected_left = hol_term::new_equals(hol_term::new_variable(set_variable), &HOL_ANY);
				} else {
					expected_left = hol_term::new_apply(hol_term::new_constant(value), hol_term::new_variable(set_variable), &HOL_ANY);
				}
				if (expected_left == nullptr)
					return (hol_term*) nullptr;
				HOL_ANY.reference_count++;

				hol_term* expected_head = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view(&expected_left, 1), make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*expected_left); free(expected_left);
					return (hol_term*) nullptr;
				}
				HOL_ANY.reference_count++;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, head, expected_head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length == 0) {
					return (hol_term*) nullptr;
				} else if (intersection.length != 1) {
					fprintf(stderr, "set_set_definition ERROR: Intersection is not unique.\n");
					free_all(intersection);
					return (hol_term*) nullptr;
				}
				return intersection[0];
			}, no_op());

	if (!result) return false;
	else return dst;
}

template<int_fast8_t ConjunctIndex>
bool set_arg(hol_term* src, hol_term*& dst, unsigned int value)
{
	head_index predicate_index; no_op apply;
	hol_term* head = find_head(src, predicate_index, find_head<built_in_predicates>, apply);
	if (head == nullptr)
		return false;

	unsigned int head_variable;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;
	while (head->type == hol_term_type::NOT)
		head = head->unary.operand;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*src, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(head_variable, &HOL_ANY));
	if (excluded_quantifier == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* excluded_tree = hol_term::new_any(nullptr, &excluded_quantifier, 1);
	if (excluded_tree == nullptr) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}

	hol_term* expected_conjunct;
	if (value == (unsigned int) built_in_predicates::UNKNOWN) {
		expected_conjunct = hol_term::new_any_right(hol_term::new_equals(
				hol_term::new_apply(hol_term::new_any_constant(
					(unsigned int) built_in_predicates::ARG1,
					(unsigned int) built_in_predicates::ARG2,
					(unsigned int) built_in_predicates::ARG3),
				hol_term::new_variable(head_variable)), excluded_tree));
	} else {
		expected_conjunct = hol_term::new_any_right(hol_term::new_equals(
				hol_term::new_apply(hol_term::new_constant(value),
				hol_term::new_variable(head_variable)), excluded_tree));
	}
	if (expected_conjunct == nullptr) {
		free(*excluded_tree); free(excluded_tree);
		return false;
	}
	excluded_tree->reference_count++;

	hol_term* expected_head;
	if (ConjunctIndex >= 0) {
		expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(
				hol_term_type::AND, excluded_tree, make_array_view((hol_term**) nullptr, 0),
				make_appended_array_view(make_repeated_array_view(excluded_tree, ConjunctIndex), expected_conjunct),
				make_array_view((hol_term**) nullptr, 0)));
		if (expected_head == nullptr) {
			free(*expected_conjunct); free(expected_conjunct);
			free(*excluded_tree); free(excluded_tree);
			return false;
		}
		excluded_tree->reference_count += ConjunctIndex;
	} else {
		unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
		expected_head = hol_term::new_exists(head_variable, hol_term::new_any_array(
				hol_term_type::AND, excluded_tree, make_array_view((hol_term**) nullptr, 0),
				make_array_view((hol_term**) nullptr, 0),
				make_prepended_array_view(expected_conjunct, make_repeated_array_view(excluded_tree, index))));
		if (expected_head == nullptr) {
			free(*expected_conjunct); free(expected_conjunct);
			free(*excluded_tree); free(excluded_tree);
			return false;
		}
		excluded_tree->reference_count += index;
	}

	array<hol_term*> result(2);
	if (value == (unsigned int) built_in_predicates::UNKNOWN)
		subtract<built_in_predicates>(result, head, expected_head);
	else intersect<built_in_predicates>(result, head, expected_head);
	free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
	if (result.length == 0) {
		return false;
	} else if (result.length != 1) {
		fprintf(stderr, "set_arg ERROR: Set operation result is not unique.\n");
		free_all(result); return false;
	}

	dst = substitute_head<any_node_position::NONE>(src, head, result[0]);
	free_all(result);
	return (dst != nullptr);
}

bool require_lambda(hol_term* src, hol_term*& dst)
{
	unsigned int lambda_variable = 0;
	if (src->type == hol_term_type::LAMBDA) {
		lambda_variable = src->quantifier.variable;
		src = src->quantifier.operand;
	} else if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
		unsigned int max_variable = 0;
		max_bound_variable(*src, max_variable);
		lambda_variable = ++max_variable;
	} else {
		return false;
	}

	hol_term* expected = hol_term::new_lambda(lambda_variable, &HOL_ANY);
	if (expected == nullptr)
		return false;
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, src, expected);
	free(*expected); if (expected->reference_count == 0) free(expected);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "require_lambda ERROR: Intersection is not unique.\n");
		return false;
	}

	dst = intersection[0];
	return true;
}

bool require_no_lambda(hol_term* src, hol_term*& dst)
{
	hol_term* excluded_lambda = hol_term::new_any_quantifier(hol_quantifier_type::LAMBDA, &HOL_ANY);
	if (excluded_lambda == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* expected = hol_term::new_any(nullptr, &excluded_lambda, 1);
	if (expected == nullptr) {
		free(*excluded_lambda); free(excluded_lambda);
		return false;
	}

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, src, expected);
	free(*expected); if (expected->reference_count == 0) free(expected);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "require_no_lambda ERROR: Intersection is not unique.\n");
		return false;
	}

	dst = intersection[0];
	return true;
}

bool set_tense(hol_term* src, hol_term*& dst, unsigned int value)
{
	hol_term* head = src;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else if (head->type == hol_term_type::ANY_ARRAY) {
		if (head->any_array.any.length > 1 || head->any_array.left.length > 1 || head->any_array.right.length > 1)
			return false;
		if (head->any_array.left.length == 1 && head->any_array.left.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.left.operands[0]->quantifier.variable;
		} else if (head->any_array.right.length == 1 && head->any_array.right.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.right.operands[0]->quantifier.variable;
		} else if (head->any_array.any.length == 1 && head->any_array.any.operands[0]->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.any.operands[0]->quantifier.variable;
		} else if (head->any_array.all->type == hol_term_type::EXISTS) {
			head_variable = head->any_array.all->quantifier.variable;
		} else {
			unsigned int max_variable = 0;
			max_bound_variable(*head, max_variable);
			head_variable = ++max_variable;
		}
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* excluded_quantifier = hol_term::new_any_quantifier(hol_quantifier_type::ANY, &HOL_ANY);
	if (excluded_quantifier == nullptr)
		return false;
	HOL_ANY.reference_count++;

	hol_term* variable = hol_term::new_any(nullptr, &excluded_quantifier, 1);
	if (variable == nullptr) {
		free(*excluded_quantifier); free(excluded_quantifier);
		return false;
	}

	hol_term* term = hol_term::new_exists(head_variable, hol_term::new_and(&HOL_ANY, hol_term::new_apply(hol_term::new_constant(value), variable)));
	if (term == nullptr) {
		free(*variable); free(variable);
		return false;
	}
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "set_tense ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

template<typename Formula>
bool set_feature(
		typename flagged_logical_form<Formula>::feature feature,
		flagged_logical_form<Formula>& exp, unsigned int value)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	hol_term* new_logical_form;
	switch (feature) {
	case feature_type::CONSTANT:
		if (!set_constant(exp.root, new_logical_form, value))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::PREDICATE:
		if (!set_predicate(exp.root, new_logical_form, value))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::PREDICATE_ONLY:
		if (!set_predicate_only(exp.root, new_logical_form, value))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::SET_DEFINITION:
		if (!set_set_definition(exp.root, new_logical_form, value))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::LEFT_ARG:
		if (!set_arg<0>(exp.root, new_logical_form, value))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::EMPTY: break;
	}
	fprintf(stderr, "set_feature ERROR: Unrecognized semantic feature.\n");
	exit(EXIT_FAILURE);
}

bool exclude_constants(hol_term* src, hol_term*& dst,
		const unsigned int* values, unsigned int count)
{
	hol_term* term = hol_term::new_any_constant_except(make_array_view(values, count));
	if (term == nullptr) return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "exclude_constants ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

bool exclude_predicates(hol_term* src, hol_term*& dst,
		const unsigned int* values, unsigned int count)
{
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS,
			hol_term::new_and(hol_term::new_apply(hol_term::new_any_constant_except(make_array_view(values, count)), &HOL_ANY), &HOL_ANY));
	if (term == nullptr) return false;
	HOL_ANY.reference_count += 2;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "exclude_predicates ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

bool exclude_predicates_only(hol_term* src, hol_term*& dst,
		const unsigned int* values, unsigned int count)
{
	hol_term* head = src;
	if ((head->type == hol_term_type::ANY || head->type == hol_term_type::ANY_RIGHT) && head->any.included != nullptr)
		head = head->any.included;

	unsigned int head_variable;
	if (head->type == hol_term_type::EXISTS) {
		head_variable = head->quantifier.variable;
	} else {
		unsigned int max_variable = 0;
		max_bound_variable(*head, max_variable);
		head_variable = ++max_variable;
	}

	hol_term* term = hol_term::new_exists(head_variable,
			hol_term::new_apply(hol_term::new_any_constant_except(make_array_view(values, count)), hol_term::new_variable(head_variable)));
	if (term == nullptr)
		return false;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, term, src);
	free(*term); if (term->reference_count == 0) free(term);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length != 1) {
		fprintf(stderr, "exclude_predicates_only ERROR: Expected intersection size to be 1.\n");
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}
	dst = intersection[0];
	return true;
}

template<typename Formula>
bool exclude_features(typename flagged_logical_form<Formula>::feature feature,
		flagged_logical_form<Formula>& exp, const unsigned int* values, unsigned int count)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	hol_term* new_logical_form;
	switch (feature) {
	case feature_type::CONSTANT:
		if (!exclude_constants(exp.root, new_logical_form, values, count))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::PREDICATE:
		if (!exclude_predicates(exp.root, new_logical_form, values, count))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::PREDICATE_ONLY:
		if (!exclude_predicates_only(exp.root, new_logical_form, values, count))
			return false;
		free(*exp.root); if (exp.root->reference_count == 0) free(exp.root);
		exp.root = new_logical_form;
		return true;
	case feature_type::SET_DEFINITION:
	case feature_type::LEFT_ARG:
		return false;
	case feature_type::EMPTY: break;
	}
	fprintf(stderr, "exclude_features ERROR: Unrecognized semantic feature.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula>
inline bool any_number(const flagged_logical_form<Formula>& src) {
	return any_number(*src.root);
}

template<typename Formula>
inline bool get_number(const flagged_logical_form<Formula>& src, int& value) {
	return get_number(*src.root, value);
}

template<typename Formula>
inline bool set_number(flagged_logical_form<Formula>& exp,
		const flagged_logical_form<Formula>& set, int value)
{
	exp.root = (Formula*) malloc(sizeof(Formula));
	if (exp.root == nullptr) return false;
	exp.flags = set.flags;
	if (!set_number(*exp.root, *set.root, value)) {
		free(exp.root);
		return false;
	}
	return true;
}

template<typename Formula>
inline bool any_string(const flagged_logical_form<Formula>& src) {
	return any_uint_list(*src.root);
}

template<typename Formula>
inline bool get_string(const flagged_logical_form<Formula>& src, sequence& value) {
	return get_uint_list(*src.root, value);
}

template<typename Formula>
inline bool set_string(flagged_logical_form<Formula>& exp,
		const flagged_logical_form<Formula>& set, const sequence& value)
{
	exp.root = (Formula*) malloc(sizeof(Formula));
	if (exp.root == nullptr) return false;
	exp.flags = set.flags;
	if (!set_uint_list(*exp.root, *set.root, value)) {
		free(exp.root);
		return false;
	}
	return true;
}

bool parse_written_hundreds(const number_parser_en& parser,
		unsigned int* terminals, unsigned int& length,
		unsigned int& index, int& integer)
{
	integer = 0;
	for (unsigned int i = 0; i < array_length(parser.ones_name_ids); i++) {
		if (terminals[index] == parser.ones_name_ids[i]) {
			integer = parser.ONES_NAMES[i].key;
			break;
		}
	}
	if (integer != 0) {
		index++;
		if (index == length)
			return true;
	}

	bool is_teen = false;
	if (integer == 0) {
		for (unsigned int i = 0; i < array_length(parser.teens_name_ids); i++) {
			if (terminals[index] == parser.teens_name_ids[i]) {
				integer = parser.TEENS_NAMES[i].key;
				is_teen = true;
				break;
			}
		}
		if (integer != 0) {
			index++;
			if (index == length)
				return true;
		}
	}

	if (terminals[index] == parser.hundred_id) {
		if (integer == 0)
			integer = 100;
		else integer *= 100;
		index++;
		if (index == length)
			return true;
		is_teen = false;
	} else if (integer != 0) {
		/* there is only a ones/teens digit */
		return true;
	}

	int tens = 0;
	for (unsigned int i = 0; i < array_length(parser.tens_name_ids); i++) {
		if (terminals[index] == parser.tens_name_ids[i]) {
			tens = parser.TENS_NAMES[i].key;
			break;
		}
	}
	if (tens != 0) {
		integer += tens;
		index++;
		if (index == length)
			return true;

		/* the next token could be a hyphen */
		if (terminals[index] == parser.hyphen_id) {
			index++;
			if (index == length)
				return false; /* dangling hyphen */
		}
	}

	if (tens == 0) {
		int teens = 0;
		for (unsigned int i = 0; i < array_length(parser.teens_name_ids); i++) {
			if (terminals[index] == parser.teens_name_ids[i]) {
				teens = parser.TEENS_NAMES[i].key;
				is_teen = true;
				break;
			}
		}
		if (teens != 0) {
			integer += teens;
			index++;
			if (index == length)
				return true;
		}
	}

	if (!is_teen) {
		int ones = 0;
		for (unsigned int i = 0; i < array_length(parser.ones_name_ids); i++) {
			if (terminals[index] == parser.ones_name_ids[i]) {
				ones = parser.ONES_NAMES[i].key;
				break;
			}
		}
		if (ones != 0) {
			integer += ones;
			index++;
			if (index == length)
				return true;
		}
	}

	return integer != 0;
}

template<typename Formula>
bool parse_written_number(const hdp_parser<Formula>& parser,
		unsigned int* terminals, unsigned int length, int& integer)
{
	if (length == 0) {
		return false;
	} else if (length == 1 && terminals[0] == parser.number_parser.zero_id) {
		integer = 0;
		return true;
	}

	unsigned int index = 0;
	unsigned int large_number_index = 0;
	integer = 0;
	while (true) {
		int leading_hundreds;
		parse_written_hundreds(parser.number_parser, terminals, length, index, leading_hundreds);
		if (index == length) {
			integer += leading_hundreds;
			return true;
		}

		/* find the next large number that matches */
		while (large_number_index < array_length(parser.number_parser.large_name_ids) && terminals[index] != parser.number_parser.large_name_ids[large_number_index])
			large_number_index++;

		if (large_number_index < array_length(parser.number_parser.large_name_ids)) {
			if (leading_hundreds != 0)
				integer += leading_hundreds * parser.number_parser.LARGE_NUMBER_NAMES[large_number_index].key;	
			else integer += parser.number_parser.LARGE_NUMBER_NAMES[large_number_index].key;
		} else {
			/* no large number matches */
			return false;
		}
	}
}


/**
 * Code to perform morphological parsing and generation.
 */


template<typename PartOfSpeechType, typename Formula>
constexpr bool morphology_is_valid(const morphology_en& morph,
		const sequence& terminal, PartOfSpeechType pos,
		const flagged_logical_form<Formula>& logical_form)
{
	return true;
}

template<bool First, typename Formula, typename PartOfSpeechType, typename EmitRootFunction>
bool morphology_parse(
		const morphology_en& morphology_parser, const sequence& words, PartOfSpeechType pos,
		const flagged_logical_form<Formula>& logical_form, EmitRootFunction emit_root)
{
	if (First) {
		/* try to decapitalize the word */
		unsigned int decapitalized_word;
		if (morphology_parser.decapitalize(words[0], decapitalized_word)) {
			sequence decapitalized_words(NULL, 0); decapitalized_words = words;
			decapitalized_words[0] = decapitalized_word;
			bool result = morphology_parse<false>(morphology_parser, decapitalized_words, pos, logical_form, emit_root);
			free(decapitalized_words);
			if (!result) return false;
		}
	}

	if (pos == POS_VERB) {
		bool contains;
		const array<inflected_verb>& forms = morphology_parser.inflected_verbs.get(words, contains);
		if (!contains) return true;

		flagged_logical_form<Formula> marked_logical_form = logical_form;
		marked_logical_form.flags.is_first_token_capital = (First ? morphology_parser.is_capitalized(words[0]) : false);
		for (const inflected_verb& form : forms) {
			if (!has_intersection(grammatical_person::THIRD, form.person)
			 || !intersect(marked_logical_form.flags.mood, logical_form.flags.mood, form.mood)
			 || !intersect(marked_logical_form.flags.index_number, logical_form.flags.index_number, form.number))
				continue;

			hol_term* new_logical_form;
			switch (form.tense) {
			case grammatical_tense::PRESENT:
				if (!set_tense(marked_logical_form.root, new_logical_form, (unsigned int) built_in_predicates::PRESENT))
					continue;
				free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
				marked_logical_form.root = new_logical_form;
				break;
			case grammatical_tense::PAST:
				if (!set_tense(marked_logical_form.root, new_logical_form, (unsigned int) built_in_predicates::PAST))
					continue;
				free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
				marked_logical_form.root = new_logical_form;
				break;
			case grammatical_tense::ANY:
				break;
			}

			if (!emit_root(form.root, marked_logical_form))
				return false;
			free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
			marked_logical_form.root = logical_form.root;
			logical_form.root->reference_count++;
		}
		return true;
	} else if (pos == POS_NOUN) {
		bool contains;
		const array<inflected_noun>& forms = morphology_parser.inflected_nouns.get(words, contains);
		if (!contains) return true;

		flagged_logical_form<Formula> marked_logical_form = logical_form;
		marked_logical_form.flags.is_first_token_capital = (First ? morphology_parser.is_capitalized(words[0]) : false);
		for (const inflected_noun& form : forms) {
			if (!intersect(marked_logical_form.flags.index_number, logical_form.flags.index_number, form.number))
				continue;

			if (!emit_root(form.root, marked_logical_form))
				return false;
			free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
			marked_logical_form.root = logical_form.root;
			logical_form.root->reference_count++;
		}
		return true;
	} else if (pos == POS_ADJECTIVE) {
		/* check if its an adverb formed from an adjective and '-ly' */
		flagged_logical_form<Formula> marked_logical_form = logical_form;
		marked_logical_form.flags.is_first_token_capital = (First ? morphology_parser.is_capitalized(words[0]) : false);
		if (intersect(marked_logical_form.flags.flags[(unsigned int) grammatical_flag::LY], logical_form.flags.flags[(unsigned int) grammatical_flag::LY], grammatical_flag_value::TRUE)) {
			bool contains;
			const array<inflected_adverb>& forms = morphology_parser.inflected_adverbs.get(words, contains);
			if (contains) {
				for (const inflected_adverb& form : forms) {
					if (!intersect(marked_logical_form.flags.comp, logical_form.flags.comp, form.comp))
						continue;

					const adverb_root& root = morphology_parser.adverbs.get(form.root);
					if (root.adj_root.length == 0)
						continue;

					if (!emit_root(root.adj_root, marked_logical_form))
						return false;
					free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
					marked_logical_form.root = logical_form.root;
					marked_logical_form.flags.flags[(unsigned int) grammatical_flag::LY] = grammatical_flag_value::TRUE;
					logical_form.root->reference_count++;
				}
			}
		}

		if (intersect(marked_logical_form.flags.flags[(unsigned int) grammatical_flag::LY], logical_form.flags.flags[(unsigned int) grammatical_flag::LY], grammatical_flag_value::FALSE)) {
			bool contains;
			const array<inflected_adjective>& forms = morphology_parser.inflected_adjectives.get(words, contains);
			if (contains) {
				for (const inflected_adjective& form : forms) {
					/* add grammatical flags for comparative, superlative, adverb formation using '-ly' */
					if (!intersect(marked_logical_form.flags.comp, logical_form.flags.comp, form.comp))
						continue;

					if (!emit_root(form.root, marked_logical_form))
						return false;
					free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
					marked_logical_form.root = logical_form.root;
					marked_logical_form.flags.flags[(unsigned int) grammatical_flag::LY] = grammatical_flag_value::FALSE;
					logical_form.root->reference_count++;
				}
			}
		}
		return true;
	} else if (pos == POS_ADVERB) {
		/* TODO: implement this */
		return false;
		//return emit_root(words, logical_form);
	} else {
		flagged_logical_form<Formula> marked_logical_form = logical_form;
		marked_logical_form.flags.is_first_token_capital = (First ? morphology_parser.is_capitalized(words[0]) : false);
		return emit_root(words, marked_logical_form);
	}
}

#endif /* HDP_PARSER_H_ */
