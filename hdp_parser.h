#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include "array_view.h"
#include "morphology_en.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

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
	grammatical_conjunction cnj;
	grammatical_flag_value flags[(uint_fast8_t) grammatical_flag::COUNT] = { grammatical_flag_value::FALSE };
	correlator corr;
	correlator correlated_by;
	auxiliary_flag aux;
	grammatical_mood mood;
	bool aux_or_subjunctive_or_inf_or_to_inf;

	grammatical_flags() :
			index_number(grammatical_num::NONE), concord_number(grammatical_num::NONE),
			cnj(grammatical_conjunction::NONE), corr(correlator::NONE),
			correlated_by(correlator::NONE), aux(auxiliary_flag::NONE),
			mood(grammatical_mood::INDICATIVE),
			aux_or_subjunctive_or_inf_or_to_inf(false) { }

	grammatical_flags(const grammatical_flags& src) {
		init(src);
	}

	inline void operator = (const grammatical_flags& src) {
		init(src);
	}

	static inline unsigned int hash(const grammatical_flags& key) {
		return default_hash(key.index_number)
			 ^ (3 * default_hash(key.concord_number))
			 ^ (7 * default_hash(key.cnj))
			 ^ (31 * default_hash(key.corr))
			 ^ (53 * default_hash(key.correlated_by))
			 ^ (97 * default_hash(key.flags, (uint_fast8_t) grammatical_flag::COUNT))
			 ^ (193 * default_hash(key.aux))
			 ^ (389 * default_hash(key.mood))
			 ^ (769 * default_hash(key.aux_or_subjunctive_or_inf_or_to_inf));
	}

	static inline void swap(grammatical_flags& first, grammatical_flags& second) {
		core::swap(first.index_number, second.index_number);
		core::swap(first.concord_number, second.concord_number);
		core::swap(first.cnj, second.cnj);
		core::swap(first.corr, second.corr);
		core::swap(first.correlated_by, second.correlated_by);
		core::swap(first.aux, second.aux);
		core::swap(first.mood, second.mood);
		core::swap(first.aux_or_subjunctive_or_inf_or_to_inf, second.aux_or_subjunctive_or_inf_or_to_inf);
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
		cnj = src.cnj;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			flags[i] = src.flags[i];
		corr = src.corr;
		correlated_by = src.correlated_by;
		aux = src.aux;
		mood = src.mood;
		aux_or_subjunctive_or_inf_or_to_inf = src.aux_or_subjunctive_or_inf_or_to_inf;
	}
};

inline bool intersect(grammatical_flags& dst,
		const grammatical_flags& first,
		const grammatical_flags& second)
{
	if (!intersect(dst.index_number, first.index_number, second.index_number)
	 || !intersect(dst.concord_number, first.concord_number, second.concord_number)
	 || !intersect(dst.cnj, first.cnj, second.cnj)
	 || !intersect(dst.corr, first.corr, second.corr)
	 || !intersect(dst.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(dst.aux, first.aux, second.aux)
	 || !intersect(dst.mood, first.mood, second.mood))
		return false;
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
		PREDICATE
	};

	enum class function_type {
		EMPTY = 0,
		IDENTITY,
		SELECT_RIGHT_CONJUNCT,
		REMOVE_RIGHT_CONJUNCT,
		SELECT_LEFT_CONJUNCT,
		REMOVE_LEFT_CONJUNCT,
		REQUIRE_NO_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE_OWN,
		REMOVE_INVERSE,
		SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,
		SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,
		SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,
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
		REQUIRE_PREDICATIVE_UNIVERSAL,
		REQUIRE_PREDICATIVE_EXISTENTIAL,
		SELECT_PREDICATE_IN_SET,
		MARK_WIDE_SCOPE,
		REQUIRE_WIDE_SCOPE,

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
		ADD_LY
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
	{feature::PREDICATE, "predicate"}
};

template<typename Formula>
const static_pair<typename flagged_logical_form<Formula>::function_type, const char*> flagged_logical_form<Formula>::FUNCTION_NAMES[] = {
	{function_type::EMPTY, "empty"},
	{function_type::IDENTITY, "identity"},
	{function_type::SELECT_RIGHT_CONJUNCT, "select_right_conjunct"},
	{function_type::REMOVE_RIGHT_CONJUNCT, "remove_right_conjunct"},
	{function_type::SELECT_LEFT_CONJUNCT, "select_left_conjunct"},
	{function_type::REMOVE_LEFT_CONJUNCT, "remove_left_conjunct"},
	{function_type::REQUIRE_NO_INVERSE, "require_no_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE, "require_left_predicate_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN, "require_left_predicate_inverse_own"},
	{function_type::REMOVE_INVERSE, "remove_inverse"},
	{function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg1_without_head_predicative_and_factor"},
	{function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg2_without_head_predicative_and_factor"},
	{function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg3_without_head_predicative_and_factor"},
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
	{function_type::REQUIRE_PREDICATIVE_UNIVERSAL, "require_predicative_universal"},
	{function_type::REQUIRE_PREDICATIVE_EXISTENTIAL, "require_predicative_existential"},
	{function_type::SELECT_PREDICATE_IN_SET, "select_predicate_in_set"},
	{function_type::MARK_WIDE_SCOPE, "mark_wide_scope"},
	{function_type::REQUIRE_WIDE_SCOPE, "require_wide_scope"},
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
	{function_type::ADD_LY, "add_ly"}
};

inline void initialize_any(grammatical_flags& flags) {
	flags.concord_number = grammatical_num::ANY;
	flags.index_number = grammatical_num::ANY;
	flags.cnj = grammatical_conjunction::ANY;
	flags.corr = correlator::ANY;
	flags.correlated_by = correlator::ANY;
	flags.aux = auxiliary_flag::ANY;
	flags.mood = grammatical_mood::ANY;
	flags.aux_or_subjunctive_or_inf_or_to_inf = false;
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
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by
	 || first.aux != second.aux
	 || first.mood != second.mood
	 || first.aux_or_subjunctive_or_inf_or_to_inf != second.aux_or_subjunctive_or_inf_or_to_inf)
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
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by
	 || first.aux != second.aux
	 || first.mood != second.mood
	 || first.aux_or_subjunctive_or_inf_or_to_inf != second.aux_or_subjunctive_or_inf_or_to_inf)
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

print("first_logical_form:  ", stderr); print(*first_logical_form, stderr, *debug_terminal_printer); print('\n', stderr);
print("second_logical_form: ", stderr); print(*second_logical_form, stderr, *debug_terminal_printer); print('\n', stderr);
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
	 || !print_feature(first, flags.cnj, out)
	 || !print_feature("corr:", first, flags.corr, out)
	 || !print_feature("corr_by:", first, flags.correlated_by, out)
	 || !print_feature(first, flags.aux, out)
	 || !print_feature(first, flags.mood, out))
		return false;

	if (flags.aux_or_subjunctive_or_inf_or_to_inf) {
		if ((!first && !print(',', out)) || !print("aux_or_subjunctive_or_inf_or_to_inf", out)) return false;
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

/* Computes the log joint probability of the grammar and given derivations */
template<typename Formula>
double log_probability(
	hdp_grammar_type<Formula>& G,
	const syntax_node<flagged_logical_form<Formula>>* const* const* syntax,
	const array<array_map<sentence<syntax_node<flagged_logical_form<Formula>>>, Formula>>& data,
	const string** reverse_name_map)
{
	typedef flagged_logical_form<Formula> logical_form_type;

	double score = 0.0;
	for (unsigned int i = 0; i < G.nonterminals.length; i++)
		score += log_probability(G.nonterminals[i].rule_distribution);
	for (unsigned int i = 0; i < data.length; i++) {
		for (unsigned int j = 0; j < data[i].size; j++) {
			logical_form_type logical_form = logical_form_type(data[i].values[j]);
			score += log_probability(G, *syntax[i][j], logical_form, reverse_name_map);
		}
	}
	return score;
}

template<typename Formula>
inline bool init(sequence& seq, const sentence<syntax_node<flagged_logical_form<Formula>>>& src)
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

template<typename Formula>
struct hdp_parser
{
	typedef flagged_logical_form<Formula> logical_form_type;

	morphology_en morph;
	hdp_grammar_type<Formula> G;
	const string** reverse_name_map;

	hdp_parser(unsigned int unknown_id,
			hash_map<string, unsigned int>& names,
			const char* morphology_filepath,
			const char* grammar_filepath) : reverse_name_map(nullptr)
	{
		if (!morph.initialize(names)
		 || !morphology_read(morph, names, morphology_filepath))
		{
			fprintf(stderr, "ERROR: Unable to initialize morphology model.\n");
			exit(EXIT_FAILURE);
		} else if (!read_grammar(G, names, grammar_filepath)) {
			fprintf(stderr, "ERROR: Unable to read grammar at '%s'.\n", grammar_filepath);
			exit(EXIT_FAILURE);
		}
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
		return true;
	}

	bool train(
			const array<array_map<sentence<syntax_node<flagged_logical_form<Formula>>>, hol_term>>& data,
			hash_map<string, unsigned int>& names,
			unsigned int iteration_count)
	{
		if (!init_capitalization_map(morph, names))
			return false;
		const string** reverse_name_map = invert(names);
		const string** nonterminal_name_map = invert(G.nonterminal_names);
		if (reverse_name_map == NULL || nonterminal_name_map == NULL) {
			if (reverse_name_map != NULL) free(reverse_name_map);
			return false;
		}
		string_map_scribe terminal_printer = { reverse_name_map, names.table.size + 1 };
		string_map_scribe nonterminal_printer = { nonterminal_name_map, G.nonterminal_names.table.size + 1 };
/* TODO: for debugging; remove this */
debug_terminal_printer = &terminal_printer;
		/* construct the initial derivation trees (running the parser with an empty grammar) */
		syntax_node<logical_form_type>*** syntax = (syntax_node<logical_form_type>***)
				calloc(data.length, sizeof(syntax_node<logical_form_type>**));
		unsigned int* order = (unsigned int*) malloc(sizeof(unsigned int) * data.length);
		if (syntax == NULL || order == NULL) {
			fprintf(stderr, "hdp_parser.train ERROR: Out of memory.\n");
			cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
			return false;
		}
		for (unsigned int i = 0; i < data.length; i++) {
			syntax[i] = (syntax_node<logical_form_type>**) calloc(data[i].size, sizeof(syntax_node<logical_form_type>*));
			if (syntax[i] == NULL) {
				fprintf(stderr, "hdp_parser.train ERROR: Out of memory.\n");
				cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
				return false;
			}
		}
		for (unsigned int i = 0; i < data.length; i++) order[i] = i;
		shuffle(order, (unsigned int) data.length);
		for (unsigned int i = 0; i < data.length; i++) {
			unsigned int id = order[i];
			for (unsigned int j = 0; j < data[id].size; j++) {
				/* TODO: add a discourse model instead of treating each sentence in each passage as independent */
				logical_form_type logical_form = logical_form_type(data[id].values[j]);
				if (data[id].keys[j].derivation != nullptr) {
					/* this sentence has a fully specified derivation tree */
					logical_form_type any(HOL_ANY);
					if (!is_parseable(*data[id].keys[j].derivation, logical_form, G, morph, any, nonterminal_printer, terminal_printer, terminal_printer.map))
					{
						fprintf(stderr, "sample ERROR: Derivation for example %u, sentence %u is not parseable: '", id, j);
						print(data[id].keys[j], stderr, terminal_printer); print("'\n", stderr);
						print(logical_form, stderr, terminal_printer); print("\n", stderr);
						cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
						return false;
					}
					syntax[id][j] = data[id].keys[j].derivation;
					syntax[id][j]->reference_count++;
				} else {
					sequence& seq = *((sequence*) alloca(sizeof(sequence)));
					if (!init(seq, data[id].keys[j])) {
						cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
						return false;
					}
					auto sentence = tokenized_sentence<logical_form_type>(seq);
					free(seq);
					syntax[id][j] = (syntax_node<logical_form_type>*) malloc(sizeof(syntax_node<logical_form_type>));
					/* NOTE: sample can set syntax[id] to null */
					if (syntax[id][j] == NULL || !sample(syntax[id][j], G, logical_form, sentence, morph, reverse_name_map) || syntax[id][j] == NULL)
					{
						fprintf(stderr, "sample ERROR: Unable to sample derivation for example %u, sentence %u: '", id, j);
						print(data[id].keys[j], stderr, terminal_printer); print("'\n", stderr);
						print(logical_form, stderr, terminal_printer); print("\n", stderr);
						cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
						return false;
					}

					print(logical_form, stdout, terminal_printer); print('\n', stdout);
					print(*syntax[id][j], stdout, nonterminal_printer, terminal_printer); print("\n\n", stdout);
					fflush(stdout);
				}

				if (!add_tree(1, *syntax[id][j], logical_form, G)) {
					cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
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
					if (data[id].keys[j].derivation != nullptr)
						/* do not resample training examples labeled with derivation trees */
						continue;
					logical_form_type logical_form = logical_form_type(data[id].values[j]);
					sequence& seq = *((sequence*) alloca(sizeof(sequence)));
					if (!init(seq, data[id].keys[j])) {
						cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
						return false;
					}
					auto sentence = tokenized_sentence<logical_form_type>(seq);
					free(seq);
					resample(syntax[id][j], G, logical_form, sentence, morph, reverse_name_map);
				}
			}
			sample_grammar(G);
			fprintf(stdout, "Unnormalized log posterior probability: %lf\n",
					log_probability(G, syntax, data, reverse_name_map));

			if (t % 1 == 0) {
				fprintf(stdout, "[iteration %u]\n", t);
				print_nonterminal_hdps(G, stdout, terminal_printer, nonterminal_printer);
				fprintf(stdout, "(seed = %u)\n", get_seed());
				fflush(stdout);
			}
		}

		/* cleanup */
		cleanup(data, reverse_name_map, nonterminal_name_map, syntax, order);
		return true;
	}

	template<unsigned int K, typename TheoryType>
	bool parse(const sentence<syntax_node<flagged_logical_form<Formula>>>& s,
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

		logical_form_type logical_form(HOL_ANY); /* TODO: initialize the features */
		syntax_node<logical_form_type>* parsed_syntax =
				(syntax_node<logical_form_type>*) alloca(K * sizeof(syntax_node<logical_form_type>));
		logical_form_type* logical_form_output =
				(logical_form_type*) alloca(K * sizeof(logical_form_type));
		sequence& seq = *((sequence*) alloca(sizeof(sequence)));
		if (!init(seq, s)) return false;
		auto sentence = tokenized_sentence<logical_form_type>(seq);
		free(seq);

		if (!::parse<false, false, K>(parsed_syntax, parse_count, logical_form,
				logical_form_output, G, sentence, morph, reverse_name_map))
			return false;

		for (unsigned int i = 0; i < parse_count; i++) {
			double log_likelihood = log_probability(G, parsed_syntax[i], logical_form_output[i], reverse_name_map);
			/* TODO: compute this prior */
			double log_prior = 0.0; //log_probability<true>(T, logical_form_output[i]);
			log_probabilities[i] = log_likelihood + log_prior;

			logical_forms[i] = (Formula*) malloc(sizeof(Formula));
			if (logical_forms[i] == NULL) {
				fprintf(stderr, "hdp_parser.parse ERROR: Out of memory.\n");
				for (unsigned int j = 0; j < i; j++) {
					free(*logical_forms[j]); free(logical_forms[j]);
				} for (unsigned int j = 0; j < parse_count; j++) {
					free(parsed_syntax[j]);
					free(logical_form_output[j]);
				}
				return false;
			} else if (!init(*logical_forms[i], *logical_form_output[i].root)) {
				free(logical_forms[i]);
				for (unsigned int j = 0; j < i; j++) {
					free(*logical_forms[j]); free(logical_forms[j]);
				} for (unsigned int j = 0; j < parse_count; j++) {
					free(parsed_syntax[j]);
					free(logical_form_output[j]);
				}
				return false;
			}
		}

		for (unsigned int j = 0; j < parse_count; j++) {
			free(parsed_syntax[j]);
			free(logical_form_output[j]);
		}
		return true;
	}

	bool add_definition(const sentence<syntax_node<flagged_logical_form<Formula>>>& s, const Formula* definition, unsigned int new_constant)
	{
	}

	template<typename Printer>
	constexpr Printer& get_printer(Printer& constant_printer) const {
		return constant_printer;
	}

private:
	void cleanup(
			const array<array_map<sentence<syntax_node<flagged_logical_form<Formula>>>, hol_term>>& data,
			const string** reverse_name_map, const string** nonterminal_name_map,
			syntax_node<logical_form_type>*** syntax, unsigned int* order)
	{
		if (syntax != NULL) {
			for (unsigned int k = 0; k < data.length; k++) {
				if (syntax[k] == NULL) continue;
				for (unsigned int l = 0; l < data[k].size; l++) {
					if (syntax[k][l] != NULL) {
						free(*syntax[k][l]);
						if (syntax[k][l]->reference_count == 0)
							free(syntax[k][l]);
					}
				}
				free(syntax[k]);
			}
			free(syntax);
		}
		if (order != NULL) free(order);
		if (reverse_name_map != NULL) free(reverse_name_map);
		if (nonterminal_name_map != NULL) free(nonterminal_name_map);
	}
};

template<typename Formula>
inline bool is_ambiguous(const flagged_logical_form<Formula>& exp) {
	return is_ambiguous(*exp.root);
}

template<typename BuiltInPredicates>
struct hol_non_head_constants {
	hol_term* constants[17];

	hol_non_head_constants() {
		constants[0] = hol_term::new_constant((unsigned int) BuiltInPredicates::UNKNOWN);
		constants[1] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG1);
		constants[2] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG2);
		constants[3] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG3);
		constants[4] = hol_term::new_constant((unsigned int) BuiltInPredicates::SIZE);
		constants[5] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT);
		constants[6] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PROGRESSIVE);
		constants[7] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PERFECT);
		constants[8] = hol_term::new_constant((unsigned int) BuiltInPredicates::PRESENT_PERFECT_PROGRESSIVE);
		constants[9] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST);
		constants[10] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PROGRESSIVE);
		constants[11] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PERFECT);
		constants[12] = hol_term::new_constant((unsigned int) BuiltInPredicates::PAST_PERFECT_PROGRESSIVE);
		constants[13] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE);
		constants[14] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PROGRESSIVE);
		constants[15] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PERFECT);
		constants[16] = hol_term::new_constant((unsigned int) BuiltInPredicates::FUTURE_PERFECT_PROGRESSIVE);
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
		|| constant == (unsigned int) BuiltInPredicates::SIZE
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
		const hol_term* src, unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.right->type == hol_term_type::VARIABLE
	 && src->binary.right->variable == scope_variable
	 && src->binary.left->type == hol_term_type::CONSTANT
	 && !is_built_in<BuiltInPredicates>(src->binary.left->constant))
	{
		return src->binary.left;
	}
	return 0;
}

enum class head_position {
	LEFT,
	RIGHT,
	ANY,
	NONE
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
inline void find_head(
		hol_term* src, hol_term*& head,
		head_index& predicate_index,
		bool& negated)
{
	head = nullptr;
	if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
		if (src->any.included == nullptr) {
			head = src;
			negated = false;
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

		if (head->quantifier.operand->type == hol_term_type::ANY)
			negated = false;

		/* find the predicate */
		hol_term* expected_predicate = hol_term::new_apply(
				hol_term::new_any(nullptr, hol_non_head_constants<BuiltInPredicates>::get_terms(), hol_non_head_constants<BuiltInPredicates>::count()),
				hol_term::new_variable(head_variable));
		if (expected_predicate == nullptr) {
			head = nullptr;
			return;
		}
		hol_non_head_constants<built_in_predicates>::increment_terms();
		hol_term* operand = head->quantifier.operand;

		hol_term* predicate = nullptr;
		predicate_index.position = head_position::NONE;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			for (unsigned int i = 0; i < operand->any_array.left.length; i++) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(operand->any_array.left.operands[i], head_variable);
				if (predicate != nullptr) {
					predicate_index = {head_position::LEFT, i};
					break;
				} else if (has_intersection<BuiltInPredicates>(operand->any_array.left.operands[i], expected_predicate)) {
					predicate_index = {head_position::LEFT, i};
				}
			} for (unsigned int i = 0; i < operand->any_array.right.length; i++) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(operand->any_array.right.operands[i], head_variable);
				if (predicate != nullptr) {
					predicate_index = {head_position::RIGHT, operand->any_array.right.length - i - 1};
					break;
				} else if (has_intersection<BuiltInPredicates>(operand->any_array.right.operands[i], expected_predicate)) {
					predicate_index = {head_position::RIGHT, operand->any_array.right.length - i - 1};
				}
			} for (unsigned int i = 0; i < operand->any_array.any.length; i++) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(operand->any_array.any.operands[i], head_variable);
				if (predicate != nullptr) {
					predicate_index = {head_position::ANY, i};
					break;
				} else if (has_intersection<BuiltInPredicates>(operand->any_array.any.operands[i], expected_predicate)) {
					predicate_index = {head_position::ANY, i};
				}
			}
		} else if (operand->type == hol_term_type::AND) {
			predicate_index = {head_position::NONE, 0};
			for (unsigned int i = 0; i < operand->array.length; i++) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(operand->array.operands[i], head_variable);
				if (predicate != nullptr) {
					predicate_index = {head_position::LEFT, i};
					break;
				} else if (predicate == nullptr && has_intersection<BuiltInPredicates>(operand->array.operands[i], expected_predicate)) {
					predicate_index = {head_position::LEFT, i};
				}
			}
		} else {
			free(*expected_predicate); free(expected_predicate);
			return;
		}
		free(*expected_predicate); free(expected_predicate);

		if (predicate_index.position == head_position::NONE) {
			head = nullptr;
			return;
		}
		negated = false;
	} else if (src->type == hol_term_type::NOT) {
		find_head<BuiltInPredicates>(src->unary.operand, head, predicate_index, negated);
		negated = true;
	}
}

template<typename BuiltInPredicates>
inline void find_head_or_universal(
		hol_term* src, hol_term*& head,
		head_index& predicate_index,
		bool& negated)
{
	hol_term* universal = hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, &HOL_ANY);
	if (universal == nullptr)
		return;
	HOL_ANY.reference_count++;

	if (has_intersection<BuiltInPredicates>(src, universal)) {
		free(*universal); free(universal);
		head = src;
		return;
	}
	free(*universal); free(universal);

	find_head<BuiltInPredicates>(src, head, predicate_index, negated);
}

template<unsigned int Constant>
inline void find_unary_application(
		hol_term* src, hol_term*& head,
		head_index& predicate_index,
		bool& negated)
{
	if (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == Constant) {
		head = src;
	} else {
		head = nullptr;
	}
}

template<typename BuiltInPredicates>
struct predicative_head_finder {
	unsigned int lambda_variable;

	predicative_head_finder(unsigned int lambda_variable) : lambda_variable(lambda_variable) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index,
			bool& negated)
	{
		apply<true>(src, head, predicate_index, negated);
	}

private:
	template<bool First>
	void apply(
			hol_term* src, hol_term*& head,
			head_index& predicate_index,
			bool& negated)
	{
		if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
			if (src->any.included == nullptr) {
				head = src;
				negated = false;
			} else {
				head = nullptr;
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
			}

			bool right_is_any = false;
			if (right->type == hol_term_type::ANY_RIGHT && right->any.included != nullptr) {
				right = right->any.included;
				right_is_any = true;
			}
			else if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
				right = right->binary.right;

			hol_term* inner_right;
			unsigned int quantified_variable = 0;
			if (right->type == hol_term_type::AND && right->array.length == 2 && right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
			 && right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[0]->binary.left->variable == set_variable
			 && right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
			{
				head = src;
				negated = false;
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
			negated = false;
		} else if (src->type == hol_term_type::NOT && First) {
			apply<false>(src->unary.operand, head, predicate_index, negated);
			negated = true;
		} else {
			head = nullptr;
		}
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
	}

	unsigned int old_arg_variable_count = keep_variables.length;
	do {
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

template<bool Factor, typename TryFindHeadFunction,
	typename MakeApplyHeadContext, typename MakeDstFunction, typename ApplyFunction>
inline bool apply_head(
		hol_term* src, hol_term*& dst,
		array<unsigned int>& dst_variables,
		array<hol_term*>& siblings,
		unsigned int max_variable,
		bool& removed_quantifier,
		bool& remove_wide_scope_marker,
		TryFindHeadFunction try_find_head,
		MakeApplyHeadContext make_context,
		MakeDstFunction make_dst,
		ApplyFunction apply)
{
	if (!apply(src)) return false;

	/* check if the current scope is the head */
	hol_term* head; head_index predicate_index; bool negated;
	try_find_head(src, head, predicate_index, negated);
	if (head != nullptr) {
		unsigned int head_variable = 0;
		if (head->type == hol_term_type::ANY) {
			if (head->any.included != nullptr)
				max_bound_variable(*head->any.included, max_variable);
			head_variable = max_variable + 1;
		} else if (head->type == hol_term_type::EXISTS) {
			head_variable = head->quantifier.variable;
		}
		dst = make_dst(head, head_variable, predicate_index, negated, remove_wide_scope_marker);
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
		if (!siblings.ensure_capacity(siblings.length + src->array.length)) return false;
		try_find_head(src->array.operands[src->array.length - 1], head, predicate_index, negated);

		if (head != nullptr) {
			unsigned int head_variable = 0;
			if (head->type == hol_term_type::ANY) {
				if (head->any.included != nullptr)
					max_bound_variable(*head->any.included, max_variable);
				head_variable = max_variable + 1;
			} else if (head->type == hol_term_type::EXISTS) {
				head_variable = head->quantifier.variable;
			}
			dst = make_dst(head, head_variable, predicate_index, negated, remove_wide_scope_marker);
			if (dst == nullptr) return false;

			get_free_variables(*dst, dst_variables);
			if (dst_variables.length > 1) insertion_sort(dst_variables);

			auto context = make_context(head);

			/* check if the other conjuncts are also part of the head scope (i.e. the head scope is a conjunction) */
			hol_term** new_dst = (hol_term**) malloc(sizeof(hol_term*) * src->array.length);
			if (new_dst == nullptr) {
				free(*dst); if (dst->reference_count == 0) free(dst);
				return false;
			}
			new_dst[src->array.length - 1] = dst;
			array<unsigned int> new_variables(1 << (core::log2(dst_variables.length == 0 ? 1 : dst_variables.length) + 1));
			new_variables.append(dst_variables.data, dst_variables.length);
			for (unsigned int i = 0; i + 1 < src->array.length; i++) {
				try_find_head(src->array.operands[i], head, predicate_index, negated);
				auto new_context = make_context(head);
				if (head == nullptr || !compatible_context(new_context, context)) {
					for (unsigned int j = 0; !Factor && j < i; j++) {
						free(*new_dst[j]);
						if (new_dst[j]->reference_count == 0)
							free(new_dst[j]);
					}
					head = nullptr;
					break;
				}

				if (!Factor) {
					if (head->type == hol_term_type::ANY) {
						head_variable = max_variable + 1;
					} else if (head->type == hol_term_type::EXISTS) {
						head_variable = head->quantifier.variable;
					}
					new_dst[i] = make_dst(head, head_variable, predicate_index, negated, remove_wide_scope_marker);
					if (new_dst[i] == nullptr) {
						for (unsigned int j = 0; j < i; j++) {
							free(*new_dst[j]);
							if (new_dst[j]->reference_count == 0)
								free(new_dst[j]);
						}
						head = nullptr;
						break;
					}

					array<unsigned int> temp_variables(8);
					get_free_variables(*dst, temp_variables);
					temp_variables.remove(temp_variables.index_of(head_variable));
					if (temp_variables.length > 1) insertion_sort(temp_variables);

					array<unsigned int> union_variables(new_variables.length + temp_variables.length);
					set_union(union_variables, new_variables, temp_variables);
					swap(union_variables, new_variables);
				}
			}

			if (head == nullptr) {
				/* only the last conjunct is the head */
				for (unsigned int i = 0; i + 1 < src->array.length; i++)
					siblings[siblings.length++] = src->array.operands[i];

				if (!prune_independent_siblings(dst_variables, siblings)) {
					free(*dst); if (dst->reference_count == 0) free(dst);
					free(new_dst); return false;
				}

				unsigned int new_start = src->array.length - 1;
				unsigned int new_length = 0;
				for (unsigned int i = src->array.length - 1; i > 0; i--) {
					if (siblings.contains(src->array.operands[i - 1])) {
						new_dst[--new_start] = src->array.operands[i];
						src->array.operands[i]->reference_count++;
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
					removed_quantifier = false;
					return true;
				} else if (new_length == src->array.length && new_dst[src->array.length - 1] == src->array.operands[src->array.length - 1]) {
					for (unsigned int i = 0; i < new_length; i++) {
						free(*new_dst[i]);
						if (new_dst[i]->reference_count == 0)
							free(new_dst[i]);
					}
					free(new_dst);
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
				removed_quantifier = false;
				return true;
			} else {
				/* all conjuncts are part of the head */
				if (!Factor) {
					swap(new_variables, dst_variables);
					if (src->type == hol_term_type::AND)
						dst = hol_term::new_and(make_array_view(new_dst, src->array.length));
					else dst = hol_term::new_or(make_array_view(new_dst, src->array.length));
					if (dst == nullptr) {
						for (unsigned int i = 0; i < src->array.length; i++) {
							free(*new_dst[i]);
							if (new_dst[i]->reference_count == 0)
								free(new_dst[i]);
						}
						free(new_dst);
						return false;
					}
				}
				free(new_dst);

				if (!prune_independent_siblings(dst_variables, siblings)) {
					free(*dst); if (dst->reference_count == 0) free(dst);
					return false;
				}
				removed_quantifier = false;
				return true;
			}
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
		dst = nullptr;
		removed_quantifier = false;
		return true;

	case hol_term_type::NOT:
		if (!apply_head<Factor>(src->unary.operand, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (new_formula == src->unary.operand) {
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
		if (!apply_head<Factor>(src->array.operands[src->array.length - 1], new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
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
			unsigned int new_length = 0;
			for (unsigned int i = src->array.length - 1; i > 0; i--) {
				if (siblings.contains(src->array.operands[i - 1])) {
					new_dst[--new_start] = src->array.operands[i];
					src->array.operands[i]->reference_count++;
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
				return true;
			} else if (new_length == src->array.length && new_formula == src->array.operands[src->array.length - 1]) {
				for (unsigned int i = 0; i < new_length; i++) {
					free(*new_dst[i]);
					if (new_dst[i]->reference_count == 0)
						free(new_dst[i]);
				}
				free(new_dst);
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
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!apply_head<Factor>(src->quantifier.operand, new_formula, dst_variables, siblings, max(max_variable, src->quantifier.variable), removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (!dst_variables.contains(src->quantifier.variable)) {
				removed_quantifier = true;
				dst = new_formula;
			} else {
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
		if (!apply_head<Factor>(src->binary.right, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
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
			if (!apply_head<Factor>(src->any.included, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
				return false;
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
				} else {
					if (src->type == hol_term_type::ANY)
						dst = hol_term::new_any(new_formula, src->any.excluded_trees, src->any.excluded_tree_count);
					else dst = hol_term::new_any_right(new_formula, src->any.excluded_trees, src->any.excluded_tree_count);
					if (dst == nullptr) {
						free(*new_formula);
						if (new_formula->reference_count == 0)
							free(new_formula);
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
		} if (!apply_head<Factor>(src->any_array.right.operands[src->any_array.right.length - 1], new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
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
				for (unsigned int i = 0; i < dst->any_array.right.length; i++)
					dst->any_array.right.operands[i]->reference_count++;
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::ANY_QUANTIFIER:
		if (!apply_head<Factor>(src->any_quantifier.operand, new_formula, dst_variables, siblings, max_variable, removed_quantifier, remove_wide_scope_marker, try_find_head, make_context, make_dst, apply))
			return false;
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

template<int_fast8_t ConjunctIndex>
hol_term* get_conjunct_context(hol_term* head) {
	if (head == nullptr) {
		return nullptr;
	} else if (head->type == hol_term_type::ANY) {
		return &HOL_ANY;
	} else if (head->type == hol_term_type::EXISTS) {
		hol_term* operand = head->quantifier.operand;
		if (operand->type == hol_term_type::ANY)
			return &HOL_ANY;

		/* find the requested conjunct */
		if (operand->type == hol_term_type::ANY_ARRAY) {
			if (ConjunctIndex >= 0) {
				if (ConjunctIndex < operand->any_array.left.length)
					return operand->any_array.left.operands[ConjunctIndex];
				else return nullptr;
			} else if (ConjunctIndex < 0) {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				if (index < operand->any_array.right.length)
					return operand->any_array.right.operands[operand->any_array.right.length - index - 1];
				else return nullptr;
			}
		} else if (operand->type == hol_term_type::AND) {
			int conjunct_index = ConjunctIndex;
			if (ConjunctIndex < 0)
				conjunct_index = ConjunctIndex + operand->array.length;
			if (conjunct_index < 0 || conjunct_index >= (int) operand->array.length)
				return nullptr;
			return operand->array.operands[conjunct_index];
		}
	}
#if !defined(NDEBUG)
	else fprintf(stderr, "get_conjunct_context ERROR: Unexpected hol_term_type.\n");
#endif
	return nullptr;
}

inline bool compatible_context(hol_term* first, hol_term* second) {
	return first == second || *first == *second;
}

struct no_op {
	template<typename... Args>
	inline constexpr bool operator() (Args&&... args) { return true; }
};

template<int_fast8_t ConjunctIndex>
inline bool select_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, get_conjunct_context<ConjunctIndex>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
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
					dst = hol_term::new_any_right(intersection[0]);
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
					if (dst != nullptr) {
						operand->array.operands[predicate_index.index]->reference_count++;
						conjunct->reference_count++;
					}
					return dst;
				}
			}, no_op()) && dst != nullptr;
}

template<int_fast8_t ConjunctIndex>
inline bool remove_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, get_conjunct_context<ConjunctIndex>,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
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

					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_any(
						hol_term::new_apply(
								hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()),
								hol_term::new_variable(head_variable)),
						excluded_trees, excluded_tree_count));
					if (dst == nullptr) {
						free(*excluded_trees[0]); free(excluded_trees[0]);
						free(*excluded_trees[1]); free(excluded_trees[1]);
						free(*excluded_trees[2]); free(excluded_trees[2]);
						free(*excluded_trees[3]); free(excluded_trees[3]);
						return (hol_term*) nullptr;
					}
					hol_non_head_constants<built_in_predicates>::increment_terms();

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
						conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
								make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
								make_excluded_array_view(operand->any_array.left.operands, operand->any_array.left.length, ConjunctIndex),
								make_array_view(operand->any_array.right.operands, operand->any_array.right.length));
					} else if (ConjunctIndex < 0) {
						unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
						index = operand->any_array.right.length - index - 1;
						conjunction = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all,
								make_array_view(operand->any_array.any.operands, operand->any_array.any.length),
								make_array_view(operand->any_array.left.operands, operand->any_array.left.length),
								make_excluded_array_view(operand->any_array.right.operands, operand->any_array.right.length, index));
					}
#if !defined(NDEBUG)
					else fprintf(stderr, "remove_conjunct WARNING: Unsupported value of `ConjunctIndex`.\n");
#endif
					if (conjunction == nullptr)
						return (hol_term*) nullptr;

					conjunction->any_array.all->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.any.length; i++)
						conjunction->any_array.any.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.left.length; i++)
						conjunction->any_array.left.operands[i]->reference_count++;
					for (unsigned int i = 0; i < conjunction->any_array.right.length; i++)
						conjunction->any_array.right.operands[i]->reference_count++;

					dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
						fprintf(stderr, "remove_conjunct WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					unsigned int conjunct_index = (ConjunctIndex < 0) ? (operand->array.length + ConjunctIndex) : ConjunctIndex;
					dst = hol_term::new_exists(head_variable, hol_term::new_and(make_excluded_array_view(operand->array.operands, operand->array.length, conjunct_index)));
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

struct empty_context { };
constexpr bool compatible_context(const empty_context& first, const empty_context& second) { return true; }
constexpr empty_context make_empty_context(const hol_term* head) { return empty_context(); }

template<int_fast8_t PredicateIndex>
inline bool require_predicate(
		hol_term* src, hol_term*& dst,
		hol_term* predicate)
{
	array_map<unsigned int, hol_term*> predicates(4);

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	bool result = apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[predicate, &predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
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

				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
					hol_term* conjunction;
					if (!new_hol_term(conjunction)) return (hol_term*) nullptr;
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = &HOL_ANY;
					if (PredicateIndex >= 0 && PredicateIndex != INT_FAST8_MAX) {
						conjunction->any_array.left.length = PredicateIndex + 1;
						conjunction->any_array.left.operands = (hol_term**) malloc(sizeof(hol_term*) * (PredicateIndex + 1));
						if (conjunction->any_array.left.operands == nullptr) {
							free(conjunction);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < PredicateIndex; i++)
							conjunction->any_array.left.operands[i] = &HOL_ANY;
						conjunction->any_array.left.operands[PredicateIndex] = current_predicate;
						conjunction->any_array.right.length = 0;
						conjunction->any_array.right.operands = nullptr;
						conjunction->any_array.any.length = 0;
						conjunction->any_array.any.operands = nullptr;
						HOL_ANY.reference_count += 1 + PredicateIndex;
					} else if (PredicateIndex < 0) {
						unsigned int index = (unsigned int) (-PredicateIndex) - 1;
						conjunction->any_array.right.length = index + 1;
						conjunction->any_array.right.operands = (hol_term**) malloc(sizeof(hol_term*) * (index + 1));
						if (conjunction->any_array.right.operands == nullptr) {
							free(conjunction);
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 1; i < index + 1; i++)
							conjunction->any_array.right.operands[i] = &HOL_ANY;
						conjunction->any_array.right.operands[0] = current_predicate;
						conjunction->any_array.left.length = 0;
						conjunction->any_array.left.operands = nullptr;
						conjunction->any_array.any.length = 0;
						conjunction->any_array.any.operands = nullptr;
						HOL_ANY.reference_count += 1 + index;
					} else if (PredicateIndex == INT_FAST8_MAX) {
						conjunction->any_array.left.length = 0;
						conjunction->any_array.left.operands = nullptr;
						conjunction->any_array.right.length = 0;
						conjunction->any_array.right.operands = nullptr;
						conjunction->any_array.any.length = 1;
						conjunction->any_array.any.operands = (hol_term**) malloc(sizeof(hol_term*) * 1);
						if (conjunction->any_array.any.operands == nullptr) {
							free(conjunction);
							return (hol_term*) nullptr;
						}
						conjunction->any_array.any.operands[0] = current_predicate;
						HOL_ANY.reference_count += 1;
					} else {
						fprintf(stderr, "require_predicate ERROR: Unsupported valueof `PredicateIndex`.\n");
						free(conjunction); return (hol_term*) nullptr;
					}
					current_predicate->reference_count++;

					hol_term* dst = hol_term::new_exists(head_variable, conjunction);
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
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
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

					hol_term* predicate = operand->array.operands[predicate_index.index];
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
						conjunction->type = hol_term_type::AND;
						conjunction->reference_count = 1;
						conjunction->array.length = operand->array.length;
						conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * operand->array.length);
						if (conjunction->array.operands == nullptr) {
							fprintf(stderr, "require_predicate ERROR: Out of memory.\n");
							for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < operand->array.length; i++) {
							if (i == predicate_index.index)
								conjunction->array.operands[i] = intersection[0];
							else
								conjunction->array.operands[i] = operand->array.operands[i];
							conjunction->array.operands[i]->reference_count++;
						}
						for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }

						dst = hol_term::new_exists(head->quantifier.variable, conjunction);
						if (dst == nullptr) {
							free(*conjunction); free(conjunction);
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

inline bool remove_inverse(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY))
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

struct arg_context {
	unsigned int head_variable;
	hol_term* operand;
};

template<int_fast8_t ConjunctIndex>
arg_context get_arg_context(hol_term* head) {
	if (head == nullptr) {
		return {0, nullptr};
	} else if (head->type == hol_term_type::ANY) {
		return {0, &HOL_ANY};
	} else if (head->type == hol_term_type::EXISTS) {
		hol_term* operand = head->quantifier.operand;
		if (operand->type == hol_term_type::ANY)
			return {0, &HOL_ANY};

		/* find the requested conjunct */
		hol_term* conjunct = nullptr;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			if (ConjunctIndex >= 0) {
				if (ConjunctIndex < operand->any_array.left.length)
					conjunct = operand->any_array.left.operands[ConjunctIndex];
				else return {0, nullptr};
			} else if (ConjunctIndex < 0) {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				if (index < operand->any_array.right.length)
					conjunct = operand->any_array.right.operands[operand->any_array.right.length - index - 1];
				else return {0, nullptr};
			}
		} else if (operand->type == hol_term_type::AND) {
			int conjunct_index = ConjunctIndex;
			if (ConjunctIndex < 0)
				conjunct_index = ConjunctIndex + operand->array.length;
			if (conjunct_index < 0 || conjunct_index >= (int) operand->array.length)
				return {0, nullptr};
			conjunct = operand->array.operands[conjunct_index];
		}

		return {head->quantifier.variable, conjunct};
	}
#if !defined(NDEBUG)
	else fprintf(stderr, "select_arg_context ERROR: Unexpected hol_term_type.\n");
#endif
	return {0, nullptr};
}

inline bool compatible_context(const arg_context& first, const arg_context& second) {
	if (first.operand == nullptr || second.operand == nullptr)
		return false;
	hol_term src_variable, dst_variable;
	src_variable.type = hol_term_type::VARIABLE;
	src_variable.variable = first.head_variable;
	dst_variable.type = hol_term_type::VARIABLE;
	dst_variable.variable = second.head_variable;
	hol_term* substituted = substitute(first.operand, &src_variable, &dst_variable);
	if (substituted == nullptr) return false;
	bool result = (*substituted == *second.operand);
	free(*substituted);
	if (substituted->reference_count == 0)
		free(substituted);
	return result;
}

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant>
inline bool select_arg_without_head_predicative_and_factor(
		hol_term* src, hol_term*& dst)
{
	unsigned int lambda_variable = 0;
	max_bound_variable(*src, lambda_variable);
	lambda_variable++;

	bool narrow_scope = false;
	array_map<unsigned int, pair<hol_term*, bool>> scopes(8);
	auto gather_scope = [&scopes,&narrow_scope](hol_term* term) {
		if (term->type == hol_term_type::FOR_ALL || term->type == hol_term_type::EXISTS || term->type == hol_term_type::LAMBDA) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = term->quantifier.variable;
			scopes.values[scopes.size++] = {term, narrow_scope};
		} else if (term->type == hol_term_type::UNARY_APPLICATION && term->binary.left->type == hol_term_type::CONSTANT && term->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
			narrow_scope = true;
		}
		return true;
	};

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	bool result = apply_head<true>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, get_arg_context<ConjunctIndex>,
			[&scopes,lambda_variable,&narrow_scope](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst = nullptr;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY))
				{
					unsigned int set_variable = lambda_variable + 1;
					unsigned int element_variable = lambda_variable + 2;
					/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
					hol_term* quantified_term;
					if (negated)
						quantified_term = hol_term::new_any_right(hol_term::new_not(hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)))));
					else quantified_term = hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)));
					if (quantified_term == nullptr)
						return (hol_term*) nullptr;
					dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&quantified_term, 1))));
					if (dst == nullptr) {
						free(*quantified_term); free(quantified_term);
						return (hol_term*) nullptr;
					}
					HOL_ANY.reference_count++;
					remove_wide_scope_marker = true;
				} else if (head->type == hol_term_type::EXISTS) {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunct;
					if (operand->type == hol_term_type::ANY_ARRAY) {
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
#if !defined(NDEBUG)
						else fprintf(stderr, "select_arg_without_head_predicative_and_factor ERROR: `operand` has type `ANY_ARRAY`, so `ConjunctIndex` should be either 0 or -1.\n");
#endif
					} else if (operand->type == hol_term_type::AND) {
						int conjunct_index = ConjunctIndex;
						if (ConjunctIndex < 0)
							conjunct_index += operand->array.length;
						conjunct = operand->array.operands[conjunct_index];
					}
#if !defined(NDEBUG)
					else fprintf(stderr, "select_arg_without_head_predicative_and_factor ERROR: Expected the head to be an existentially-quantified conjunction.\n");
#endif

					hol_term* expected_left_side = hol_term::new_apply(hol_term::new_constant(ArgConstant), hol_term::new_variable(head->quantifier.variable));
					if (expected_left_side == nullptr)
						return (hol_term*) nullptr;
					hol_term* expected_conjunct = hol_term::new_any(hol_term::new_equals(expected_left_side, &HOL_ANY));
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

					if (conjunct->type == hol_term_type::ANY) {
						free(*expected_left_side); free(expected_left_side);
						unsigned int set_variable = lambda_variable + 1;
						unsigned int element_variable = lambda_variable + 2;
						/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
						hol_term* quantified_term;
						if (negated)
							quantified_term = hol_term::new_any_right(hol_term::new_not(hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)))));
						else quantified_term = hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)));
						if (quantified_term == nullptr)
							return (hol_term*) nullptr;
						dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&quantified_term, 1))));
						if (dst == nullptr) {
							free(*quantified_term); free(quantified_term);
							return (hol_term*) nullptr;
						}
						HOL_ANY.reference_count++;
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
							/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
							hol_term* quantified_term;
							if (negated)
								quantified_term = hol_term::new_any_right(hol_term::new_not(hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)))));
							else quantified_term = hol_term::new_any_right(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)));
							if (quantified_term == nullptr)
								return (hol_term*) nullptr;
							dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&quantified_term, 1))));
							if (dst == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							HOL_ANY.reference_count++;
							remove_wide_scope_marker = true;
						} else if (arg->type == hol_term_type::CONSTANT) {
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
							hol_term* quantified_term;
							if (negated)
								quantified_term = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)
									)));
							else quantified_term = hol_term::new_any_right(hol_term::new_not(hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)
									))));
							if (quantified_term == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;
							hol_term* set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, hol_term::new_equals(element_var, hol_term::new_constant(arg->constant))));
							if (set_definition == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
							set_var->reference_count++;
							dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
									make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition, 1), make_array_view(&quantified_term, 1))));
							if (dst == nullptr) {
								free(*set_definition); free(set_definition);
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							HOL_ANY.reference_count++;
							remove_wide_scope_marker = true;
						} else if (arg->type == hol_term_type::VARIABLE) {
							unsigned int element_variable = arg->variable;
							pair<hol_term*, bool> scope = scopes.get(element_variable);
							hol_term* operand = scope.key->quantifier.operand;
							if (scope.key->type == hol_term_type::FOR_ALL && operand->type == hol_term_type::IF_THEN) {
								hol_term* antecedent = operand->binary.left;
								/* check if the scope is a relativized quantification over a set */
								if (antecedent->type == hol_term_type::UNARY_APPLICATION && antecedent->binary.right->type == hol_term_type::VARIABLE
								 && antecedent->binary.right->variable == element_variable && antecedent->binary.left->type == hol_term_type::VARIABLE)
								{
									dst = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (dst == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;
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
									if (!scope.value && narrow_scope) {
										quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
										remove_wide_scope_marker = true;
									} else if (narrow_scope) {
										quantified_term = hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)));
									} else {
										quantified_term = hol_term::new_any_right(hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									}
									if (quantified_term == nullptr) {
										free(*set_definition); free(set_definition);
										return (hol_term*) nullptr;
									}
									set_var->reference_count++;
									arg->reference_count += 2;
									if (narrow_scope) {
										dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term));
									} else {
										dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term)));
										remove_wide_scope_marker = true;
									}
									if (dst == nullptr) {
										free(*set_definition); free(set_definition);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									HOL_ANY.reference_count++;
								}
							} else if (scope.key->type == hol_term_type::EXISTS && operand->type == hol_term_type::AND) {
								hol_term* first_conjunct = operand->array.operands[0];
								/* check if the scope is a relativized quantification over a set */
								if (first_conjunct->type == hol_term_type::UNARY_APPLICATION && first_conjunct->binary.right->type == hol_term_type::VARIABLE
								 && first_conjunct->binary.right->variable == element_variable && first_conjunct->binary.left->type == hol_term_type::VARIABLE)
								{
									dst = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (dst == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;
								} else {
									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr)
										return (hol_term*) nullptr;
									hol_term* set_definition = hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, hol_term::new_and(make_array_view(operand->array.operands, operand->array.length - 1))));
									if (set_definition == nullptr) {
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									for (unsigned int i = 0; i + 1 < operand->array.length; i++)
										operand->array.operands[i]->reference_count++;
									hol_term* quantified_term;
									if (!scope.value && narrow_scope) {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									} else if (narrow_scope) {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)));
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
									if (narrow_scope) {
										dst = hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term));
									} else {
										dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_and(set_definition, quantified_term)));
									}
									if (dst == nullptr) {
										free(*set_definition); free(set_definition);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									HOL_ANY.reference_count++;
									remove_wide_scope_marker = true;
								}
							} else if (scope.key->type == hol_term_type::EXISTS && operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::AND || operand->any_array.oper == hol_term_type::ANY_ARRAY)) {
								hol_term* first_conjunct = (operand->any_array.left.length == 0 ? nullptr : operand->any_array.left.operands[0]);
								/* check if the scope is a relativized quantification over a set */
								if (first_conjunct != nullptr && first_conjunct->type == hol_term_type::UNARY_APPLICATION && first_conjunct->binary.right->type == hol_term_type::VARIABLE
								 && first_conjunct->binary.right->variable == element_variable && first_conjunct->binary.left->type == hol_term_type::VARIABLE)
								{
									dst = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (dst == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;
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
									if (!scope.value && narrow_scope) {
										quantified_term = hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									} else {
										quantified_term = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg))));
									}
									if (quantified_term == nullptr) {
										free(*set_definition_term); free(set_definition_term);
										return (hol_term*) nullptr;
									}
									set_var->reference_count++;
									arg->reference_count += 2;
									if (!scope.value && narrow_scope)
										dst = hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
												make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition_term, 1), make_array_view(&quantified_term, 1)));
									else dst = hol_term::new_any_right(hol_term::new_exists(set_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
												make_array_view((hol_term**) nullptr, 0), make_array_view(&set_definition_term, 1), make_array_view(&quantified_term, 1))));
									if (dst == nullptr) {
										free(*set_definition_term); free(set_definition_term);
										free(*quantified_term); free(quantified_term);
										return (hol_term*) nullptr;
									}
									HOL_ANY.reference_count++;
									remove_wide_scope_marker = true;
								}
							} else if (scope.key->type == hol_term_type::LAMBDA) {
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
								fprintf(stderr, "select_arg_without_head_predicative_and_factor WARNING: `diff` has more than one logical form.\n");
#endif

							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, diff[0])),
									hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)))
								));
							for (hol_term* term : diff) { free(*term); if (term->reference_count == 0) free(term); }
							if (dst == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;
							set_var->reference_count += 2 - 1;
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

							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, new_inner_operand)),
									hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)))
								));
							if (dst == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								free(*new_inner_operand); free(new_inner_operand);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;
							set_var->reference_count += 2 - 1;
						} else if (inner_operand->type == hol_term_type::AND) {
							free(*equals_term); if (equals_term->reference_count == 0) free(equals_term);
							hol_term* set_definition;
							if (inner_operand->array.length - 1 == 1) {
								set_definition = inner_operand->array.operands[0];
							} else {
								set_definition = hol_term::new_and(make_array_view(inner_operand->array.operands, inner_operand->array.length - 1));
							}
							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, set_definition)),
									hol_term::new_exists(element_variable, hol_term::new_and(
										hol_term::new_apply(set_var, element_var),
										hol_term::new_apply(hol_term::new_variable(lambda_variable), element_var)))
								));
							if (dst == nullptr) {
								free(*element_var); free(element_var);
								free(*set_var); free(set_var);
								return (hol_term*) nullptr;
							}
							element_var->reference_count += 2 - 1;
							set_var->reference_count += 2 - 1;
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

	if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[input_tense_predicates,output_tense_predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[tense_predicates](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
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

inline bool remove_perfect(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				if (!negated)
					return (hol_term*) nullptr;

				hol_term* dst = head;
				dst->reference_count++;
				return dst;
			}, no_op()) && dst != nullptr;
}

inline bool require_no_empty_ref(hol_term* src, hol_term*& dst)
{
	hol_term* empty_ref = hol_term::new_constant((unsigned int) built_in_predicates::EMPTY_REF);
	if (empty_ref == nullptr)
		return false;

	hol_term* predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, &empty_ref, 1),
			hol_term::new_variable(0));
	if (predicate == nullptr) {
		free(*empty_ref); free(empty_ref);
		return false;
	}

	bool result = require_predicate<INT_FAST8_MAX>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool predicate_only(hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
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

inline bool predicate_and_tense(hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* var_term = hol_term::new_variable(head_variable);
				if (var_term == nullptr)
					return (hol_term*) nullptr;

				hol_term* expected_head = hol_term::new_exists(head_variable, hol_term::new_and(
					hol_term::new_apply(
						hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), var_term),
					hol_term::new_apply(hol_term::new_any_constant(make_array_view(PAST_OR_PRESENT, array_length(PAST_OR_PRESENT))), var_term)));
				if (expected_head == nullptr) {
					free(*var_term); free(var_term);
					return (hol_term*) nullptr;
				}
				hol_non_head_constants<built_in_predicates>::increment_terms();
				var_term->reference_count += 2 - 1;

				array<hol_term*> intersection(2);
				intersect<built_in_predicates>(intersection, head, expected_head);
				free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
				if (intersection.length > 1) {
					fprintf(stderr, "predicate_and_tense ERROR: Expected intersection size to be 1.\n");
					for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
					return (hol_term*) nullptr;
				} else if (intersection.length == 0) {
					return (hol_term*) nullptr;
				}
				remove_wide_scope_marker = true;
				return intersection[0];
			}, no_op()) && dst != nullptr;
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable), make_empty_context,
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst = nullptr;
				unsigned int set_variable, element_variable;
				bool wide_scope_marker_before_quantifier = false;
				bool wide_scope_marker_before_lambda_term = false;
				if (head->type == hol_term_type::ANY) {
					set_variable = head_variable;
					element_variable = head_variable + 1;
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
						} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
							last = last->binary.right;
							wide_scope_marker_before_quantifier = true;
						}

						while (last->type == hol_term_type::NOT)
							last = operand->unary.operand;

						if (last->type == hol_term_type::ANY_RIGHT) {
							element_variable = ++max_variable;
						} else if (last->type != QuantifierType && !last_is_any) {
							return (hol_term*) nullptr;
						} else if (last_is_any && last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::VARIABLE
								&& last->binary.left->variable == lambda_variable && last->binary.right->type == hol_term_type::VARIABLE)
						{
							element_variable = last->binary.right->variable;
						} else {
							dst = head;
							head->reference_count++;

							hol_term* inner_right = nullptr;
							if (QuantifierType == hol_term_type::EXISTS && last->quantifier.operand->type == hol_term_type::AND)
								inner_right = last->quantifier.operand->array.operands[1];
							else if (QuantifierType == hol_term_type::FOR_ALL && last->quantifier.operand->type == hol_term_type::IF_THEN)
								inner_right = last->quantifier.operand->binary.right;

							if (inner_right != nullptr && inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
								wide_scope_marker_before_lambda_term = true;
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

					if (wide_scope_marker_before_lambda_term) {
						hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), lambda_apply_term);
						if (temp == nullptr) {
							free(*lambda_apply_term); free(lambda_apply_term);
							return (hol_term*) nullptr;
						}
						lambda_apply_term = temp;
					}

					hol_term* quantifier;
					if (QuantifierType == hol_term_type::FOR_ALL) {
						quantifier = hol_term::new_for_all(element_variable, hol_term::new_if_then(
								hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term));
					} else if (QuantifierType == hol_term_type::EXISTS) {
						quantifier = hol_term::new_exists(element_variable, hol_term::new_and(
								hol_term::new_apply(hol_term::new_variable(set_variable), element_var), lambda_apply_term));
					} else {
						fprintf(stderr, "require_predicative_quantifier ERROR: Unsupported `QuantifierType`.\n");
						free(*lambda_apply_term); free(lambda_apply_term);
						return (hol_term*) nullptr;
					}
					if (quantifier == nullptr) {
						free(*lambda_apply_term); free(lambda_apply_term);
						return (hol_term*) nullptr;
					}
					element_var->reference_count++;

					if (wide_scope_marker_before_quantifier) {
						hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), quantifier);
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

template<bool UniqueOutput, bool StorePredicates>
inline bool select_predicate_in_set(array<hol_term*>& out, hol_term* head, hol_term* predicate, unsigned int max_variable, unsigned int lambda_variable)
{
	unsigned int set_variable, element_variable;
	if (head->type == hol_term_type::ANY) {
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
				return false;
			}
		}
	}

	hol_term* first_expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
			hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, hol_term::new_equals(hol_term::new_variable(element_variable), predicate))),
			&HOL_ANY));
	if (first_expected_head == nullptr)
		return false;
	predicate->reference_count++;
	HOL_ANY.reference_count++;

	intersect<built_in_predicates>(out, first_expected_head, head);
	free(*first_expected_head); if (first_expected_head->reference_count == 0) free(first_expected_head);
	if (UniqueOutput && out.length > 1) {
		fprintf(stderr, "select_predicate_in_set ERROR: Intersection is not unique.\n");
		for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	for (unsigned int i = 0; StorePredicates && i < out.length; i++) {
		hol_term* term = out[i]->quantifier.operand->array.operands[0]->binary.right->quantifier.operand->binary.right;
		term->reference_count++;
		free(*out[i]); if (out[i]->reference_count == 0) free(out[i]);
		out[i] = term;
	}

	if (UniqueOutput && out.length == 1)
		return true;

	hol_term* second_expected_head = hol_term::new_exists(set_variable, hol_term::new_and(
			hol_term::new_equals(hol_term::new_variable(set_variable), hol_term::new_lambda(element_variable, hol_term::new_apply(predicate, hol_term::new_variable(element_variable)))),
			&HOL_ANY));
	if (second_expected_head == nullptr)
		return false;
	predicate->reference_count++;
	HOL_ANY.reference_count++;

	unsigned int old_out_length = out.length;
	intersect<built_in_predicates>(out, second_expected_head, head);
	free(*second_expected_head); if (second_expected_head->reference_count == 0) free(second_expected_head);
	if (UniqueOutput && out.length > 1) {
		fprintf(stderr, "select_predicate_in_set ERROR: Intersection is not unique.\n");
		for (hol_term* term : out) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	for (unsigned int i = old_out_length; StorePredicates && i < out.length; i++) {
		hol_term* term = out[i]->quantifier.operand->array.operands[0]->binary.right->quantifier.operand->binary.left;
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable), make_empty_context,
			[max_variable, lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				array<hol_term*> out(2);
				if (!select_predicate_in_set<true, true>(out, head, &HOL_ANY, max_variable, lambda_variable) || out.length == 0)
					return (hol_term*) nullptr;
				remove_wide_scope_marker = true;
				return out[0];
			}, no_op()) && dst != nullptr;
}

inline bool mark_wide_scope(hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool removed_quantifier, remove_wide_scope_marker = false;
	bool result = apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			find_head_or_universal<built_in_predicates>, make_empty_context,
			[](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* new_head = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), head);
				if (new_head == nullptr) return (hol_term*) nullptr;
				head->reference_count++;
				if (head->type == hol_term_type::ANY) {
					hol_term* any = hol_term::new_any(new_head);
					if (any == nullptr) {
						free(*new_head); free(new_head);
						return (hol_term*) nullptr;
					}
					new_head = any;
				} else if (head->type == hol_term_type::ANY_RIGHT) {
					hol_term* any = hol_term::new_any_right(new_head);
					if (any == nullptr) {
						free(*new_head); free(new_head);
						return (hol_term*) nullptr;
					}
					new_head = any;
				}
				return new_head;
			}, no_op()) && dst != nullptr;
	if (!result) return false;

	hol_term* any_wide_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
	if (any_wide_scope == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		return false;
	}
	HOL_ANY.reference_count++;

	hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
	if (duplicate_wide_scopes == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		free(*any_wide_scope); free(any_wide_scope);
		return false;
	}

	hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
	if (excluded_universal == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		return false;
	}
	any_wide_scope->reference_count++;

	hol_term* excluded_trees[2];
	excluded_trees[0] = excluded_universal;
	excluded_trees[1] = duplicate_wide_scopes;
	hol_term* unique_wide_scope = hol_term::new_any_right(
			hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY),
			excluded_trees, array_length(excluded_trees));
	if (unique_wide_scope == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		free(*excluded_universal); free(excluded_universal);
		free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
		return false;
	}
	HOL_ANY.reference_count++;

	array<hol_term*> intersection(2);
	intersect<built_in_predicates>(intersection, dst, unique_wide_scope);
	free(*dst); if (dst->reference_count == 0) free(dst);
	free(*unique_wide_scope); if (unique_wide_scope->reference_count == 0) free(unique_wide_scope);
	if (intersection.length == 0) {
		return false;
	} else if (intersection.length > 1) {
		fprintf(stderr, "mark_wide_scope ERROR: Intersection is not unique.\n");
		return false;
	}
	dst = intersection[0];
	return true;
}

template<bool WideScope>
inline bool require_narrow_or_wide_scope(hol_term* src, hol_term*& dst)
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
	bool removed_quantifier, remove_wide_scope_marker = false;
	return apply_head<false>(src, dst, dst_variables, siblings, 0, removed_quantifier, remove_wide_scope_marker,
			predicative_head_finder<built_in_predicates>(lambda_variable), make_empty_context,
			[&max_variable,lambda_variable](hol_term* head, unsigned int head_variable, head_index predicate_index, bool negated, bool& remove_wide_scope_marker)
			{
				hol_term* dst = nullptr;
				unsigned int set_variable, element_variable;
				if (head->type == hol_term_type::ANY) {
					set_variable = head_variable;
					element_variable = head_variable + 1;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS) {
						fprintf(stderr, "require_narrow_or_wide_scope ERROR: Expected existential quantification of set.\n");
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
								fprintf(stderr, "require_narrow_or_wide_scope ERROR: Expected an existentially quantified conjunction with a right-most operand.\n");
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
						} else if (last->type == hol_term_type::UNARY_APPLICATION && last->binary.left->type == hol_term_type::CONSTANT && last->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
							last = last->binary.right;
							if (WideScope) {
								head->reference_count++;
								return head;
							} else {
								return (hol_term*) nullptr;
							}
						}

						while (last->type == hol_term_type::NOT)
							last = operand->unary.operand;

						if (last->type == hol_term_type::ANY) {
							element_variable = ++max_variable;
						} else if (last->type == hol_term_type::EXISTS || last->type == hol_term_type::FOR_ALL) {
							dst = head;
							head->reference_count++;

							hol_term* inner_right = nullptr;
							if (last->type == hol_term_type::EXISTS && last->quantifier.operand->type == hol_term_type::AND)
								inner_right = last->quantifier.operand->array.operands[1];
							else if (last->type == hol_term_type::FOR_ALL && last->quantifier.operand->type == hol_term_type::IF_THEN)
								inner_right = last->quantifier.operand->binary.right;

							if (inner_right != nullptr && inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
								if (WideScope) {
									return (hol_term*) nullptr;
								} else {
									head->reference_count++;
									return head;
								}
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
						fprintf(stderr, "require_narrow_or_wide_scope ERROR: Intersection is not unique.\n");
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
		dst.flags.cnj = grammatical_conjunction::NONE;
		dst.flags.corr = correlator::NONE;
		dst.flags.correlated_by = correlator::NONE;
		dst.flags.aux = auxiliary_flag::NONE;
		dst.flags.mood = grammatical_mood::INDICATIVE;
		dst.flags.aux_or_subjunctive_or_inf_or_to_inf = false;
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
	case function_type::REMOVE_LEFT_CONJUNCT:
		dst.flags = src.flags;
		return remove_conjunct<0>(src.root, dst.root);
	case function_type::REQUIRE_NO_INVERSE:
		dst.flags = src.flags;
		return require_no_inverse(src.root, dst.root);
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE:
		dst.flags = src.flags;
		return require_left_predicate_inverse(src.root, dst.root);
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN:
		dst.flags = src.flags;
		return require_left_predicate_inverse_own(src.root, dst.root);
	case function_type::REMOVE_INVERSE:
		dst.flags = src.flags;
		return remove_inverse(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		dst.flags = src.flags;
		return select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG1>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		dst.flags = src.flags;
		return select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG2>(src.root, dst.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		dst.flags = src.flags;
		return select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG3>(src.root, dst.root);
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
	case function_type::PREDICATE_ONLY:
		dst.flags = src.flags;
		return predicate_only(src.root, dst.root);
	case function_type::PREDICATE:
		dst.flags = src.flags;
		return predicate(src.root, dst.root);
	case function_type::PREDICATE_AND_TENSE:
		dst.flags = src.flags;
		return predicate_and_tense(src.root, dst.root);
	case function_type::SELECT_PREDICATE_IN_SET:
		dst.flags = src.flags;
		return select_predicate_in_set(src.root, dst.root);
	case function_type::MARK_WIDE_SCOPE:
		dst.flags = src.flags;
		return mark_wide_scope(src.root, dst.root);
	case function_type::REQUIRE_WIDE_SCOPE:
		dst.flags = src.flags;
		return require_narrow_or_wide_scope<true>(src.root, dst.root);
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
	}
	fprintf(stderr, "apply ERROR: Unrecognized transformation function.\n");
	return false;
}

enum class any_node_position {
	LEFT,
	RIGHT,
	NONE
};

template<any_node_position AnyNodePosition>
struct head_substituter {
	bool found_wide_scope;
	const hol_term* src;
	hol_term* dst;
	unsigned int last_declared_variable;

	head_substituter(bool found_wide_scope, const hol_term* src, hol_term* dst) :
			found_wide_scope(found_wide_scope), src(src), dst(dst), last_declared_variable(0) { }
};

inline hol_term* wrap_any_right(hol_term* operand, unsigned int prev_declared_variable, bool found_wide_scope)
{
	if (prev_declared_variable != 0) {
		unsigned int excluded_tree_count = 3 + (found_wide_scope ? 1 : 0);
		hol_term* excluded_trees[4];
		excluded_trees[0] = hol_term::new_any(hol_term::new_for_all(prev_declared_variable, &HOL_ANY));
		excluded_trees[1] = hol_term::new_any(hol_term::new_exists(prev_declared_variable, &HOL_ANY));
		excluded_trees[2] = hol_term::new_any(hol_term::new_lambda(prev_declared_variable, &HOL_ANY));
		if (found_wide_scope)
			excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
		if (excluded_trees[0] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[1] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[2] != nullptr) HOL_ANY.reference_count++;
		if (found_wide_scope && excluded_trees[3] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || (found_wide_scope && excluded_trees[3] == nullptr)) {
			if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
			if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
			if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
			return nullptr;
		}

		hol_term* term = hol_term::new_any_right(operand, excluded_trees, excluded_tree_count);
		if (term == nullptr) {
			free(*excluded_trees[0]); free(excluded_trees[0]);
			free(*excluded_trees[1]); free(excluded_trees[1]);
			free(*excluded_trees[2]); free(excluded_trees[2]);
			free(*excluded_trees[3]); free(excluded_trees[3]);
			return nullptr;
		}
		return term;
	} else if (found_wide_scope) {
		hol_term* excluded = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
		if (excluded == nullptr) return nullptr;
		HOL_ANY.reference_count++;
		hol_term* term = hol_term::new_any_right(operand, &excluded, 1);
		if (term == nullptr) {
			free(*excluded); free(excluded);
			return nullptr;
		}
		return term;
	} else {
		return hol_term::new_any_right(operand);
	}
}

template<hol_term_type Type, any_node_position AnyNodePosition>
inline hol_term* apply(hol_term* src, head_substituter<AnyNodePosition>& substituter)
{
	if (src == substituter.src) {
		if (AnyNodePosition == any_node_position::LEFT) {
			hol_term* new_term = wrap_any_right(substituter.dst, substituter.last_declared_variable, substituter.found_wide_scope);
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
	bool found_wide_scope;
	unsigned int prev_declared_variable;
	switch (Type) {
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && src->binary.right->type != hol_term_type::ANY && src->binary.right->type != hol_term_type::ANY_RIGHT) {
			found_wide_scope = substituter.found_wide_scope;
			prev_declared_variable = substituter.last_declared_variable;
			substituter.found_wide_scope = false;
			substituter.last_declared_variable = 0;
		}
		substituter.found_wide_scope = (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
		first = src->binary.left;
		second = apply(src->binary.right, substituter);
		if (second == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && src->binary.right->type != hol_term_type::ANY && src->binary.right->type != hol_term_type::ANY_RIGHT) {
			hol_term* term = wrap_any_right(second, prev_declared_variable, found_wide_scope);
			if (term == nullptr) {
				if (second != src->binary.right) { free(*second); if (second->reference_count == 0) free(second); }
				return nullptr;
			}
			if (second == src->binary.right) second->reference_count++;
			second = term;
		}

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
		if (AnyNodePosition == any_node_position::RIGHT || AnyNodePosition == any_node_position::LEFT) {
			found_wide_scope = substituter.found_wide_scope;
			prev_declared_variable = substituter.last_declared_variable;
			substituter.found_wide_scope = false;
			substituter.last_declared_variable = 0;
		}
		first = src->ternary.first;
		second = src->ternary.second;
		third = apply(src->ternary.third, substituter);
		if (third == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT) {
			hol_term* term = wrap_any_right(third, prev_declared_variable, found_wide_scope);
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

		if (AnyNodePosition == any_node_position::LEFT) {
			hol_term* term = wrap_any_right(new_term, prev_declared_variable, found_wide_scope);
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
			if (first->type == Type && first != src->array.operands[src->array.length - 1]) {
				/* the new child is the same type as the parent, so merge them */
				new_term = (hol_term*) realloc(first->array.operands, sizeof(hol_term*) * (src->array.length + first->array.length - 1));
				if (new_term == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					return nullptr;
				}
				for (unsigned int i = new_term->array.length; i > 0; i--)
					new_term->array.operands[i + src->array.length - 1] = new_term->array.operands[i];
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.operands[src->array.length - 1] = first;
				new_term->array.length = src->array.length + first->array.length - 1;
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
		if (src->any_array.right.length == 0)
			return nullptr;
		changed = false;
		first = apply(src->any_array.right.operands[src->any_array.right.length - 1], substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->any_array.right.operands[src->any_array.right.length - 1]) {
			changed = true;
		}

		if (!changed) {
			new_term = src;
		} else {
			new_term = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
					make_array_view(src->any_array.any.operands, src->any_array.any.length),
					make_array_view(src->any_array.left.operands, src->any_array.left.length),
					make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), first));
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
		if ((AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA) || (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::FOR_ALL)) {
			found_wide_scope = substituter.found_wide_scope;
			prev_declared_variable = substituter.last_declared_variable;
			substituter.found_wide_scope = false;
			substituter.last_declared_variable = 0;
		}
		if (!(AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA))
			substituter.last_declared_variable = src->quantifier.variable;
		first = apply(src->quantifier.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA) {
			hol_term* term = wrap_any_right(first, src->quantifier.variable, found_wide_scope);
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

		if (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::FOR_ALL) {
			hol_term* term = wrap_any_right(new_term, prev_declared_variable, found_wide_scope);
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
			found_wide_scope = substituter.found_wide_scope;
			prev_declared_variable = substituter.last_declared_variable;
			substituter.found_wide_scope = false;
			substituter.last_declared_variable = 0;
		}
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

		if (AnyNodePosition == any_node_position::LEFT && src->any_quantifier.quantifier != hol_quantifier_type::EXISTS) {
			hol_term* term = wrap_any_right(new_term, prev_declared_variable, found_wide_scope);
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
			found_wide_scope = substituter.found_wide_scope;
			prev_declared_variable = substituter.last_declared_variable;
			substituter.found_wide_scope = false;
			substituter.last_declared_variable = 0;
		}
		first = apply(src->unary.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT) {
			hol_term* term = wrap_any_right(first, prev_declared_variable, found_wide_scope);
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

		if (AnyNodePosition == any_node_position::LEFT) {
			hol_term* term = wrap_any_right(new_term, prev_declared_variable, found_wide_scope);
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
		if (src->any.included == nullptr)
			return nullptr;
		first = apply(src->any.included, substituter);
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

template<any_node_position AnyNodePosition>
inline hol_term* substitute_head(hol_term* src,
		const hol_term* src_term, hol_term* dst_term)
{
	head_substituter<AnyNodePosition> substituter(false, src_term, dst_term);
	hol_term* dst = apply(src, substituter);

	if (AnyNodePosition == any_node_position::RIGHT && dst->type != hol_term_type::LAMBDA && dst->type != hol_term_type::ANY_RIGHT && dst->type != hol_term_type::ANY) {
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
	hol_term* head; head_index predicate_index; bool negated;
	remover.try_find_head(src, head, predicate_index, negated);
	if (head != nullptr) {
		if (src->type == hol_term_type::ANY_RIGHT && src->any.included == head) {
			hol_term* term = src->any.included;
			term->reference_count++;
			return term;
		} else {
			return src;
		}
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
			if (first->type == Type && first != src->array.operands[src->array.length - 1]) {
				/* the new child is the same type as the parent, so merge them */
				new_term = (hol_term*) realloc(first->array.operands, sizeof(hol_term*) * (src->array.length + first->array.length - 1));
				if (new_term == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					return nullptr;
				}
				for (unsigned int i = new_term->array.length; i > 0; i--)
					new_term->array.operands[i + src->array.length - 1] = new_term->array.operands[i];
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.operands[src->array.length - 1] = first;
				new_term->array.length = src->array.length + first->array.length - 1;
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
		if (src->any_array.right.length == 0)
			return nullptr;
		changed = false;
		first = apply(src->any_array.right.operands[src->any_array.right.length - 1], remover);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->any_array.right.operands[src->any_array.right.length - 1]) {
			changed = true;
		}

		if (!changed) {
			new_term = src;
		} else {
			new_term = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
					make_array_view(src->any_array.any.operands, src->any_array.any.length),
					make_array_view(src->any_array.left.operands, src->any_array.left.length),
					make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), first));
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
	if (src->type != hol_term_type::LAMBDA && !(src->type == hol_term_type::ANY_RIGHT && src->any.excluded_tree_count > 0)) {
		if (src->type == hol_term_type::ANY_RIGHT) {
			if (src->any.included != nullptr)
				src = src->any.included;
#if !defined(NDEBUG)
			else fprintf(stderr, "remove_any_nodes ERROR: Unexpected null `included` field in `src`.\n");
#endif
		}
	}

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
		}
		outer.add(node);
	}
};

template<typename TryFindHeadFunction, typename Function>
hol_term* find_head(hol_term* term, head_index& predicate_index, TryFindHeadFunction& try_find_head, Function& apply)
{
	apply(term);

	hol_term* head; bool negated;
	try_find_head(term, head, predicate_index, negated);
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
		if (term->any.included == nullptr) return nullptr;
		return find_head(term->any.included, predicate_index, try_find_head, apply);

	case hol_term_type::ANY_ARRAY:
		if (term->any_array.right.length == 0) return nullptr;
		return find_head(term->any_array.right.operands[term->any_array.right.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY_QUANTIFIER:
		return find_head(term->any_quantifier.operand, predicate_index, try_find_head, apply);
	}
	fprintf(stderr, "find_head ERROR: Unrecognied hol_term_type.\n");
	return nullptr;
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
		unsigned int& max_variable)
{
	apply_head_inverter first_inverter; head_index first_predicate_index;
	apply_head_inverter second_inverter; head_index second_predicate_index;
print("first:  ", stderr); print(*first, stderr, *debug_terminal_printer); print('\n', stderr);
print("second: ", stderr); print(*second, stderr, *debug_terminal_printer); print('\n', stderr);
	hol_term* first_head = find_head(first, first_predicate_index, find_first_head, first_inverter);
	hol_term* second_head = find_head(second, second_predicate_index, find_second_head, second_inverter);
	if (first_head == nullptr || second_head == nullptr)
		return nullptr;

print("first_head:  ", stderr); print(*first_head, stderr, *debug_terminal_printer); print('\n', stderr);
print("second_head: ", stderr); print(*second_head, stderr, *debug_terminal_printer); print('\n', stderr);
	max_variable = first_inverter.max_variable;
	max_bound_variable(*first_head, max_variable);

	array_map<unsigned int, unsigned int> second_variable_map(8);
	for (unsigned int i = first_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = first_inverter.outer[i - 1];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;
		if (!second_variable_map.put(node->quantifier.variable, ++max_variable))
			return nullptr;
	}

	for (unsigned int i = second_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = second_inverter.outer[i - 1];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;
		if (node->quantifier.variable <= max_variable) {
			if (!second_variable_map.put(node->quantifier.variable, ++max_variable))
				return nullptr;
		}
	}

	if (!on_remap_variables(first_head, second_head, second_variable_map, max_variable))
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
	hol_term* remapped_second = remap_to_invert_apply_head(first, second, find_first_head, find_second_head, on_remap_variables, max_variable);
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

	hol_term* remapped_first = substitute_head<any_node_position::RIGHT>(first, first_head, first_head);
	if (remapped_first == nullptr) {
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		return false;
	}

	hol_term* conjunct = nullptr;
	array<hol_term*> new_second_heads(8);
	if (!invert_second_head(new_second_heads, first_head, second_head, first_predicate_index, second_predicate_index, conjunct, max_variable)) {
		free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
		return false;
	}

	if (second_inverter.outer.length > 1) {
		hol_term* parent = second_inverter.outer[second_inverter.outer.length - 2];
		if (parent->type == hol_term_type::AND || parent->type == hol_term_type::OR) {
			/* check if the conjunction/disjunction is part of the head */
			array<hol_term*>* new_second_head_array = (array<hol_term*>*) malloc(sizeof(array<hol_term*>) * parent->array.length);
			if (new_second_head_array == nullptr) {
				fprintf(stderr, "invert_apply_head ERROR: Out of memory.\n");
				free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
				free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
				if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
				for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			for (unsigned int i = 0; i < parent->array.length; i++) {
				if (!array_init(new_second_head_array[i], 8)) {
					free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
					for (unsigned int j = 0; j < i; j++) free(new_second_head_array[j]);
					free(new_second_head_array);
					return false;
				}
			}
			for (unsigned int i = 0; i + 1 < parent->array.length; i++) {
				if (!invert_second_head(new_second_head_array[i], first_head, parent->array.operands[i], first_predicate_index, second_predicate_index, conjunct, max_variable)) {
					for (unsigned int j = 0; j < i; j++) {
						for (hol_term* term : new_second_head_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
					}
					free(new_second_head_array);
					new_second_head_array = nullptr;
					break;
				}
			}

			if (new_second_head_array != nullptr) {
				swap(new_second_head_array[parent->array.length - 1], new_second_heads);
				unsigned int intersection_count = new_second_head_array[0].length;
				for (unsigned int j = 1; j < parent->array.length; j++)
					intersection_count *= new_second_head_array[j].length;
				unsigned int* index_array = (unsigned int*) calloc(parent->array.length, sizeof(unsigned int));
				if (!new_second_heads.ensure_capacity(intersection_count) || index_array == nullptr) {
					free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
					free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
					if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
					for (unsigned int j = 0; j < parent->array.length; j++) {
						for (hol_term* term : new_second_head_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
					}
					free(new_second_head_array);
					return false;
				}
				while (true) {
					hol_term*& new_head = new_second_heads[new_second_heads.length];
					if (!new_hol_term(new_head)) {
						free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
						free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
						if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
						for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
						for (unsigned int j = 0; j < parent->array.length; j++) {
							for (hol_term* term : new_second_head_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
						}
						free(new_second_head_array); free(index_array);
						return false;
					}
					new_head->type = parent->type;
					new_head->reference_count = 1;
					new_head->array.length = parent->array.length;
					new_head->array.operands = (hol_term**) malloc(sizeof(hol_term*) * parent->array.length);
					if (new_head->array.operands == nullptr) {
						free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
						free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
						if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }
						for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
						for (unsigned int j = 0; j < parent->array.length; j++) {
							for (hol_term* term : new_second_head_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
						}
						free(new_second_head_array); free(index_array); free(new_head);
						return false;
					}
					for (unsigned int i = 0; i < parent->array.length; i++) {
						new_head->array.operands[i] = new_second_head_array[i][index_array[i]];
						new_head->array.operands[i]->reference_count++;
					}
					new_second_heads.length++;

					/* increment `index_array` */
					bool remaining = false;
					for (unsigned int j = parent->array.length; j > 0; j--) {
						index_array[j - 1]++;
						if (index_array[j - 1] == new_second_head_array[j - 1].length) {
							index_array[j - 1] = 0;
						} else {
							remaining = true;
							break;
						}
					}
					if (!remaining) break;
				}
				for (unsigned int j = 0; j < parent->array.length; j++) {
					for (hol_term* term : new_second_head_array[j]) { free(*term); if (term->reference_count == 0) free(term); }
				}
				free(new_second_head_array); free(index_array);
			}
		}
	}
	if (conjunct != nullptr) { free(*conjunct); if (conjunct->reference_count == 0) free(conjunct); }

	array<hol_term*> new_seconds(new_second_heads.length);
	for (hol_term* new_second_head : new_second_heads) {
		new_seconds[new_seconds.length] = substitute_head<any_node_position::LEFT>(remapped_second, second_head, new_second_head);
		if (new_seconds[new_seconds.length] == nullptr) {
			for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : new_seconds) { free(*term); if (term->reference_count == 0) free(term); }
			free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
			free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
			return false;
		}
		new_seconds.length++;
	}
	for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }
	free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);

	array<hol_term*> intersection(16);
	for (hol_term* new_second : new_seconds)
{
print("new_second:     ", stderr); print(*new_second, stderr, *debug_terminal_printer); print('\n', stderr);
print("remapped_first: ", stderr); print(*remapped_first, stderr, *debug_terminal_printer); print('\n', stderr);
		intersect<built_in_predicates>(intersection, remapped_first, new_second);
}
	for (hol_term* term : new_seconds) { free(*term); if (term->reference_count == 0) free(term); }
	free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);

	if (intersection.length == 0)
		return false;

	inverse = (flagged_logical_form<hol_term>*) malloc(sizeof(flagged_logical_form<hol_term>) * intersection.length);
	if (inverse == nullptr) {
		for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		return false;
	}

	for (unsigned int i = 0; i < intersection.length; i++) {
		inverse[i].flags = flags;
fprintf(stderr, "intersection[%u]: ", i); print(*intersection[i], stderr, *debug_terminal_printer); print('\n', stderr);
		inverse[i].root = remove_any_nodes(intersection[i], find_first_head);
		if (inverse[i].root == nullptr) {
			for (unsigned int j = 0; j < i; j++) {
				free(*inverse[j].root); if (inverse[j].root->reference_count == 0) free(inverse[j].root);
			} for (unsigned int j = i; j < intersection.length; j++) {
				free(*intersection[j]); if (intersection[j]->reference_count == 0) free(intersection[j]);
			}
			free(inverse);
			return false;
		}
fprintf(stderr, "inverse[%u]:      ", i); print(inverse[i], stderr, *debug_terminal_printer); print('\n', stderr);
		free(*intersection[i]); if (intersection[i]->reference_count == 0) free(intersection[i]);
	}
	inverse_count = intersection.length;
	return true;
}

template<int_fast8_t ConjunctIndex>
inline bool invert_select_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	auto on_remap_variables = [](hol_term* first_head, hol_term* second_head, array_map<unsigned int, unsigned int>& second_variable_map, unsigned int max_variable) {
		hol_term* expected_head;
		hol_term* hol_any_ptr = &HOL_ANY;
		if (ConjunctIndex >= 0) {
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(hol_any_ptr, ConjunctIndex + 1), make_array_view((hol_term**) nullptr, 0)));
			if (expected_head == nullptr)
				return false;
			HOL_ANY.reference_count += 2 + ConjunctIndex;
		} else {
			unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_repeated_array_view(hol_any_ptr, index + 1)));
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
		array<unsigned int> second_variables(8);
		if (!get_variables(*second_conjunct, second_variables)) {
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
			for (hol_term* term : second_intersection) { free(*term); if (term->reference_count == 0) free(term); }
			return false;
		}

		array<pair<hol_term*, array_map<unsigned int, variable_set>>> conjunct_intersection(2);
print("second_conjunct: ", stderr); print(*second_conjunct, stderr, *debug_terminal_printer); print('\n', stderr);
print("conjunct:        ", stderr); print(*conjunct, stderr, *debug_terminal_printer); print('\n', stderr);
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
		const array_map<unsigned int, variable_set>& variable_map = conjunct_intersection[0].value;
		for (unsigned int var : second_variables) {
			unsigned int index = variable_map.index_of(var);
			if (index == variable_map.size) continue;

			unsigned int target_var = 0;
			const variable_set& target_var_set = variable_map.values[index];
			switch (target_var_set.type) {
			case variable_set_type::SINGLETON:
				target_var = target_var_set.variable;
				break;
			case variable_set_type::ANY:
				for (unsigned int i = 0; i < target_var_set.any.length; i++) {
					if (target_var_set.any.array[i] >= max_variable) {
						target_var = target_var_set.any.array[i];
						max_variable = target_var + 1;
						break;
					}
				}
				break;
			case variable_set_type::ANY_EXCEPT:
				target_var = max_variable;
				for (unsigned int i = 0; i < target_var_set.any.length; i++)
					target_var = max(target_var, target_var_set.any.array[i]);
				target_var++;
				max_variable = target_var + 1;
				break;
			}
			if (target_var == 0) {
				free_all(conjunct_intersection);
				return false;
			}

			if (!second_variable_map.ensure_capacity(second_variable_map.size + 1)) {
				free_all(conjunct_intersection);
				return false;
			}
			index = second_variable_map.index_of(var);
			if (target_var == var) {
				/* `second_variable_map` should map `var` to itself */
				if (index < second_variable_map.size)
					second_variable_map.remove_at(index);
			} else {
				/* `second_variable_map` should map `var` to `target_var` */
				second_variable_map.values[index] = target_var;
				if (index == second_variable_map.size) {
					second_variable_map.keys[index] = var;
					second_variable_map.size++;
				}
			}
		}
		free_all(conjunct_intersection);
		return true;
	};

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, on_remap_variables,
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;
				if (first_head->type == hol_term_type::ANY || (first_head->type == hol_term_type::EXISTS
				 && (first_head->quantifier.operand->type == hol_term_type::ANY || first_head->quantifier.operand->type == hol_term_type::ANY_ARRAY)))
				{
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

					hol_term* conjunct = hol_term::new_any(nullptr, excluded_trees, excluded_tree_count);
					if (conjunct == nullptr) {
						free(*excluded_trees[0]); free(excluded_trees[0]);
						free(*excluded_trees[1]); free(excluded_trees[1]);
						free(*excluded_trees[2]); free(excluded_trees[2]);
						free(*excluded_trees[3]); free(excluded_trees[3]);
						return false;
					}
					HOL_ANY.reference_count++;

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
					intersect<built_in_predicates>(dst, new_second_head, first_head);
					free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);
					return (dst.length > 0);
				} else if (first_head->type == hol_term_type::EXISTS && first_head->quantifier.operand->type == hol_term_type::AND) {
					hol_term* first_head_operand = first_head->quantifier.operand;

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
					if (!dst.ensure_capacity(dst.length + intersection_count)) {
						for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					for (hol_term* predicate : predicates) {
						for (hol_term* conjunct : conjuncts) {
							hol_term* conjunction;
							if (!new_hol_term(conjunction)) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								return false;
							}
							conjunction->type = hol_term_type::AND;
							conjunction->reference_count = 1;
							conjunction->array.length = first_head_operand->array.length;
							conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first_head_operand->array.length);
							if (conjunction->array.operands == nullptr) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								free(conjunction); return false;
							}
							for (unsigned int i = 0; i < first_head_operand->array.length; i++) {
								if (i == first_predicate_index.index) {
									conjunction->array.operands[i] = predicate;
								} else if (i == (unsigned int) conjunct_index) {
									conjunction->array.operands[i] = conjunct;
								} else {
									conjunction->array.operands[i] = first_head_operand->array.operands[i];
								}
								conjunction->array.operands[i]->reference_count++;
							}

							dst[dst.length] = hol_term::new_exists(first_head->quantifier.variable, conjunction);
							if (dst[dst.length] == nullptr) {
								free(*conjunction); free(conjunction);
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								return false;
							}
							dst.length++;
						}
					}
					for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
					for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
					return (dst.length > 0);
				} else {
					fprintf(stderr, "invert_select_conjunct ERROR: Unexpected type of `first_head`.\n");
					return false;
				}
			} else {
				fprintf(stderr, "invert_select_conjunct ERROR: Expected `second_head` to be an existentially quantified conjunction.\n");
				return false;
			}
		});
}

template<int_fast8_t ConjunctIndex>
inline bool invert_remove_conjunct(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;
				hol_term* head_var = hol_term::new_variable(second_head->quantifier.variable);
				if (head_var == nullptr) return false;
				constexpr unsigned int excluded_tree_count = 5;
				hol_term* excluded_trees[excluded_tree_count];
				excluded_trees[0] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), head_var));
				excluded_trees[1] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), head_var));
				excluded_trees[2] = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), head_var));
				excluded_trees[3] = hol_term::new_any(hol_term::new_exists(second_head->quantifier.variable, &HOL_ANY));
				excluded_trees[4] = hol_term::new_apply(
						hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), head_var);
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

				hol_term* conjunction;
				if (!new_hol_term(conjunction)) {
					free(*conjunct); free(conjunct);
					return false;
				}
				conjunction->type = hol_term_type::AND;
				conjunction->reference_count = 1;
				conjunction->array.length = second_head_operand->array.length + 1;
				conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * (second_head_operand->array.length + 1));
				if (conjunction->array.operands == nullptr) {
					free(*conjunct); free(conjunct);
					free(conjunction); return false;
				}
				int conjunct_index = ConjunctIndex;
				if (ConjunctIndex < 0) conjunct_index += second_head_operand->array.length + 1;
				for (unsigned int i = 0; i < conjunction->array.length; i++) {
					if (i == (unsigned int) conjunct_index) {
						conjunction->array.operands[i] = conjunct;
					} else {
						unsigned int second_head_index = (i > (unsigned) conjunct_index) ? (i - 1) : i;
						conjunction->array.operands[i] = second_head_operand->array.operands[second_head_index];
						conjunction->array.operands[i]->reference_count++;
					}
				}

				if (!dst.ensure_capacity(dst.length + 1)) {
					free(*conjunction); free(conjunction);
					return false;
				}
				dst[dst.length] = hol_term::new_exists(second_head->quantifier.variable, conjunction);
				if (dst[dst.length] == nullptr) {
					free(dst[dst.length]);
					free(*conjunction); free(conjunction);
					return false;
				}
				dst.length++;
				return true;
			} else {
				fprintf(stderr, "invert_remove_conjunct ERROR: Expected `second_head` to be an existentially quantified conjunction.\n");
				return false;
			}
		});
}

inline bool invert_remove_inverse(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
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
				}
				for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
				return (dst.length > 0);
			} else if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
				hol_term* second_head_operand = second_head->quantifier.operand;

				hol_term* second_predicate;
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
				}
				for (hol_term* term : new_predicates) { free(*term); if (term->reference_count == 0) free(term); }
				return (dst.length > 0);
			} else {
				fprintf(stderr, "invert_remove_inverse ERROR: Expected `second_head` to be an existentially quantified conjunction.\n");
				return false;
			}
		});
}

template<int_fast8_t ConjunctIndex, unsigned int ArgConstant>
inline bool invert_select_arg_without_head_predicative_and_factor(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	if (second->type != hol_term_type::LAMBDA)
		return false;
	unsigned int lambda_variable = second->quantifier.variable;
	second = second->quantifier.operand;

	auto on_remap_variables = [](hol_term* first_head, hol_term* second_head, array_map<unsigned int, unsigned int>& second_variable_map, unsigned int max_variable) {
		hol_term* expected_arg = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), &HOL_ANY), &HOL_ANY);
		if (expected_arg == nullptr)
			return false;
		HOL_ANY.reference_count += 2;

		hol_term* expected_head;
		hol_term* hol_any_ptr = &HOL_ANY;
		if (ConjunctIndex >= 0) {
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(hol_any_ptr, ConjunctIndex), expected_arg), make_array_view((hol_term**) nullptr, 0)));
			if (expected_head == nullptr) {
				free(*expected_arg); free(expected_arg);
				return false;
			}
			HOL_ANY.reference_count += 1 + ConjunctIndex;
		} else {
			unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
			expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
					make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(expected_arg, make_repeated_array_view(hol_any_ptr, index))));
			if (expected_head == nullptr) {
				free(*expected_arg); free(expected_arg);
				return false;
			}
			HOL_ANY.reference_count += 1 + index;
		}

		array<hol_term*> intersection(2);
		unsigned int first_element_variable;
		intersect<built_in_predicates>(intersection, first_head, expected_head);
		if (intersection.length != 0) {
			free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);

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

			if (conjunct->binary.right->type == hol_term_type::ANY) {
				/* the "element variable" is not defined in `first_head` so we don't have to worry about variable agreement */
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (conjunct->binary.right->type == hol_term_type::VARIABLE) {
				first_element_variable = conjunct->binary.right->variable;
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

			hol_term* hol_any_ptr = &HOL_ANY;
			if (ConjunctIndex >= 0) {
				expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(hol_any_ptr, ConjunctIndex), expected_conjunct), make_array_view((hol_term**) nullptr, 0)));
				if (expected_head == nullptr) {
					free(*expected_conjunct); free(expected_conjunct);
					return false;
				}
				HOL_ANY.reference_count += 1 + ConjunctIndex;
			} else {
				unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
				expected_head = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
						make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(expected_conjunct, make_repeated_array_view(hol_any_ptr, index))));
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
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return true;
			} else if (inner_conjunct->binary.right->type == hol_term_type::VARIABLE) {
				first_element_variable = inner_conjunct->binary.right->variable;
			} else {
				/* this should be a variable */
				for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			for (hol_term* term : intersection) { free(*term); if (term->reference_count == 0) free(term); }
		}

		unsigned int second_element_variable;
		if (second_head->type == hol_term_type::EXISTS) {
			unsigned int set_variable = second_head->quantifier.variable;
			hol_term* operand = second_head->quantifier.operand;

			hol_term* right;
			if (operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::ANY_ARRAY || operand->any_array.oper == hol_term_type::AND)) {
				right = (operand->any_array.right.length == 0 ? operand->any_array.all : operand->any_array.right.operands[operand->any_array.right.length - 1]);
			} else if (operand->type == hol_term_type::AND) {
				right = operand->array.operands[operand->array.length - 1];
			} else {
				return (operand->type == hol_term_type::ANY);
			}

			if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
				right = right->binary.right;
			if (right->type == hol_term_type::FOR_ALL) {
				second_element_variable = right->quantifier.variable;
			} else if (right->type == hol_term_type::EXISTS) {
				second_element_variable = right->quantifier.variable;
			} else if (right->type == hol_term_type::AND && right->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
					&& right->array.operands[0]->binary.left->type == hol_term_type::VARIABLE && right->array.operands[0]->binary.left->variable == set_variable
					&& right->array.operands[0]->binary.right->type == hol_term_type::VARIABLE)
			{
				second_element_variable = right->array.operands[0]->binary.right->variable;
			}
		} else {
			return (second_head->type == hol_term_type::ANY || (second_head->type == hol_term_type::ANY_QUANTIFIER && has_intersection(hol_quantifier_type::EXISTS, second_head->any_quantifier.quantifier)));
		}

		for (unsigned int i = 0; i < second_variable_map.size; i++) {
			if (second_variable_map.keys[i] == second_element_variable) {
				second_variable_map.values[i] = first_element_variable;
				if (first_element_variable == second_element_variable)
					second_variable_map.remove_at(i);
				break;
			}
		}
		return true;
	};

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, predicative_head_finder<built_in_predicates>(lambda_variable), on_remap_variables,
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable)
		{
			unsigned int predicate_variable;
			if (first_head->type == hol_term_type::ANY) {
				predicate_variable = max_variable + 1;
			} else {
#if !defined(NDEBUG)
				if (first_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor WARNING: Expected `first_head` to be an existential quantification.\n");
#endif
				predicate_variable = first_head->quantifier.variable;
			}

			if (second_head->type == hol_term_type::ANY) {
				if (!dst.ensure_capacity(dst.length + 2))
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

				hol_term* hol_any_ptr = &HOL_ANY;
				if (ConjunctIndex >= 0) {
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(hol_any_ptr, ConjunctIndex), arg), make_array_view((hol_term**) nullptr, 0)));
					if (dst[dst.length] == nullptr)
						return false;
					HOL_ANY.reference_count += 1 + ConjunctIndex;
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg, make_repeated_array_view(hol_any_ptr, index))));
					if (dst[dst.length] == nullptr)
						return false;
					HOL_ANY.reference_count += 1 + index;
				} else {
					fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: Unsupported value of `ConjunctIndex`.\n");
					return false;
				}
				dst.length++;

				unsigned int conjunct_variable = max_variable + 2;
				hol_term* conjunct = hol_term::new_exists(conjunct_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_array_view(&arg, 1)));
				if (conjunct == nullptr)
					return false;
				arg->reference_count++;
				HOL_ANY.reference_count++;
				if (ConjunctIndex >= 0) {
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(hol_any_ptr, ConjunctIndex), conjunct), make_array_view((hol_term**) nullptr, 0)));
					HOL_ANY.reference_count += 1 + ConjunctIndex;
				} else if (ConjunctIndex < 0) {
					unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
					dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, &HOL_ANY,
							make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(conjunct, make_repeated_array_view(hol_any_ptr, index))));
					HOL_ANY.reference_count += 1 + index;
				} else {
					fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: Unsupported value of `ConjunctIndex`.\n");
					return false;
				}
				if (dst[dst.length] == nullptr) {
					free(*conjunct); free(conjunct);
					return false;
				}
				dst.length++;
				return true;
			} else {
#if !defined(NDEBUG)
				if (second_head->type != hol_term_type::EXISTS)
					fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: Expected `second_head` to be an existential quantification.\n");
#endif
				unsigned int set_variable = second_head->quantifier.variable;
				hol_term* operand = second_head->quantifier.operand;

				hol_term* left = nullptr;
				hol_term* right = nullptr;
				bool can_keep_set_variable = true;
				bool can_remove_set_variable = true;
				bool must_be_simple_set_def;
				if (operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::ANY_ARRAY || operand->any_array.oper == hol_term_type::AND)) {
					left = (operand->any_array.left.length == 0 ? operand->any_array.all : operand->any_array.left.operands[0]);
					right = (operand->any_array.right.length == 0 ? operand->any_array.all : operand->any_array.right.operands[operand->any_array.right.length - 1]);
					must_be_simple_set_def = (left->type == hol_term_type::EQUALS && left->binary.left->type == hol_term_type::VARIABLE && left->binary.left->variable == set_variable && left->binary.right->type == hol_term_type::LAMBDA);
				} else if (operand->type == hol_term_type::AND) {
					left = operand->array.operands[0];
					right = operand->array.operands[operand->array.length - 1];
					must_be_simple_set_def = (left->type == hol_term_type::EQUALS && left->binary.left->type == hol_term_type::VARIABLE && left->binary.left->variable == set_variable && left->binary.right->type == hol_term_type::LAMBDA);
					if (operand->array.length == 2 && must_be_simple_set_def)
						can_keep_set_variable = false;
					if (operand->array.length > 2)
						can_remove_set_variable = false;
				}

				if (right->type == hol_term_type::ANY) {
					fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: `ANY` type is unsupported for the right conjunct.");
					return false;
				}

				if (!dst.ensure_capacity(dst.length + 2))
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

				hol_term* element_var;
				unsigned int element_variable;
				bool wide_scope_marker = false;
				bool element_narrow_scope = false;
				if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					right = right->binary.right;
					element_narrow_scope = true;
					wide_scope_marker = true;
				}
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
				}

				if (must_be_simple_set_def && left->binary.right->quantifier.operand->type == hol_term_type::EQUALS
				 && left->binary.right->quantifier.operand->binary.left->type == hol_term_type::VARIABLE
				 && left->binary.right->quantifier.operand->binary.left->variable == left->binary.right->quantifier.variable
				 && (left->binary.right->quantifier.operand->binary.right->type == hol_term_type::CONSTANT
				  || left->binary.right->quantifier.operand->binary.right->type == hol_term_type::ANY_CONSTANT
				  || left->binary.right->quantifier.operand->binary.right->type == hol_term_type::ANY)
				 && right->type != hol_term_type::FOR_ALL && right->type != hol_term_type::AND)
				{
					/* the set contains a single constant */
					hol_term* constant = left->binary.right->quantifier.operand->binary.right;
					hol_term* arg_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), head_var), constant);
					if (arg_term == nullptr) {
						free(*conjunct); free(conjunct);
						free(*element_var); free(element_var);
						return false;
					}
					head_var->reference_count++;
					constant->reference_count++;

					array<hol_term*> arg_terms(2);
					intersect<built_in_predicates>(arg_terms, arg_term, conjunct);
					free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
					for (hol_term* arg_term : arg_terms) {
						if (!dst.ensure_capacity(dst.length + 1)) {
							for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); free(conjunct);
							free(*element_var); free(element_var);
							return false;
						}
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
							fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: Unsupported value of `ConjunctIndex`.\n");
							for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); free(conjunct);
							free(*element_var); free(element_var);
							return false;
						}
						if (dst[dst.length] == nullptr) {
							for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
							free(*conjunct); free(conjunct);
							free(*element_var); free(element_var);
							return false;
						}
						arg_term->reference_count++;
						dst.length++;
					}
					for (hol_term* term : arg_terms) { free(*term); if (term->reference_count == 0) free(term); }
					if (left->binary.right->quantifier.operand->binary.right->type != hol_term_type::ANY) {
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*element_var); free(element_var);
						return true;
					}
				}

				hol_term* set_var = hol_term::new_variable(set_variable);
				if (set_var == nullptr) {
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
					free(*element_var); free(element_var);
					return false;
				}

				hol_term* arg_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ArgConstant), head_var), element_var);
				if (arg_term == nullptr) {
					free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
					free(*element_var); free(element_var);
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
					if (!dst.ensure_capacity(dst.length + left_intersections.length + 1)) {
						for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						return false;
					}
					for (hol_term* left_intersection : left_intersections) {
						hol_term* set_definition = left_intersection->binary.right->quantifier.operand;

						if (right->type == hol_term_type::FOR_ALL) {
							hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_for_all(element_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							HOL_ANY.reference_count++;

							hol_term* narrow_scope;
							if (ConjunctIndex >= 0) {
								narrow_scope = hol_term::new_any_right(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0))), &excluded_quantifier, 1);
								if (narrow_scope == nullptr) {
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
								narrow_scope = hol_term::new_any_right(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index)))), &excluded_quantifier, 1);
								if (narrow_scope == nullptr) {
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

							dst[dst.length] = hol_term::new_for_all(element_variable, hol_term::new_if_then(set_definition, narrow_scope));
							if (dst[dst.length] == nullptr) {
								for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*narrow_scope); free(narrow_scope);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							set_definition->reference_count++;
							dst.length++;
						} else if (wide_scope_marker && !element_narrow_scope) {
							hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
							if (excluded_quantifier == nullptr) {
								for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							HOL_ANY.reference_count++;

							hol_term* narrow_scope;
							if (ConjunctIndex >= 0) {
								narrow_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), arg_term), make_array_view((hol_term**) nullptr, 0))))), &excluded_quantifier, 1);
								if (narrow_scope == nullptr) {
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
								narrow_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right(hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
										make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(arg_term, make_repeated_array_view(conjunct, index)))))), &excluded_quantifier, 1);
								if (narrow_scope == nullptr) {
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

							if (set_definition->type == hol_term_type::AND) {
								dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(make_array_view(set_definition->array.operands, set_definition->array.length), narrow_scope)));
								if (dst[dst.length] == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*narrow_scope); free(narrow_scope);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								for (unsigned int i = 0; i < set_definition->array.length; i++)
									set_definition->array.operands[i]->reference_count++;
							} else if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::AND) {
								dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_any_array(hol_term_type::AND, set_definition->any_array.all,
										make_array_view(set_definition->any_array.any.operands, set_definition->any_array.any.length),
										make_array_view(set_definition->any_array.left.operands, set_definition->any_array.left.length),
										make_appended_array_view(make_array_view(set_definition->any_array.right.operands, set_definition->any_array.right.length), narrow_scope)));
								if (dst[dst.length] == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*narrow_scope); free(narrow_scope);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								set_definition->any_array.all->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.left.length; i++)
									set_definition->any_array.left.operands[i]->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.right.length; i++)
									set_definition->any_array.right.operands[i]->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.any.length; i++)
									set_definition->any_array.any.operands[i]->reference_count++;
							} else {
								if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::ANY_ARRAY)
									fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: The set definition has `ANY_ARRAY` type `oper` is not `AND`.\n");
								dst[dst.length] = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, narrow_scope));
								if (dst[dst.length] == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*narrow_scope); free(narrow_scope);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								set_definition->reference_count++;
							}
							dst.length++;
						} else {
							hol_term* new_definition;
							if (set_definition->type == hol_term_type::AND) {
								new_definition = hol_term::new_exists(element_variable, hol_term::new_and(make_appended_array_view(make_array_view(set_definition->array.operands, set_definition->array.length), arg_term)));
								if (new_definition == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								arg_term->reference_count++;
								for (unsigned int i = 0; i < set_definition->array.length; i++)
									set_definition->array.operands[i]->reference_count++;
							} else if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::AND) {
								new_definition = hol_term::new_exists(element_variable, hol_term::new_any_array(hol_term_type::AND, set_definition->any_array.all,
										make_array_view(set_definition->any_array.any.operands, set_definition->any_array.any.length),
										make_array_view(set_definition->any_array.left.operands, set_definition->any_array.left.length),
										make_appended_array_view(make_array_view(set_definition->any_array.right.operands, set_definition->any_array.right.length), arg_term)));
								if (new_definition == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								arg_term->reference_count++;
								set_definition->any_array.all->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.left.length; i++)
									set_definition->any_array.left.operands[i]->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.right.length; i++)
									set_definition->any_array.right.operands[i]->reference_count++;
								for (unsigned int i = 0; i < set_definition->any_array.any.length; i++)
									set_definition->any_array.any.operands[i]->reference_count++;
							} else {
								if (set_definition->type == hol_term_type::ANY_ARRAY && set_definition->any_array.oper == hol_term_type::ANY_ARRAY)
									fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: The set definition has `ANY_ARRAY` type `oper` is not `AND`.\n");
								new_definition = hol_term::new_exists(element_variable, hol_term::new_and(set_definition, arg_term));
								if (new_definition == nullptr) {
									for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
									free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
									free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
									free(*set_var); if (set_var->reference_count == 0) free(set_var);
									return false;
								}
								arg_term->reference_count++;
								set_definition->reference_count++;
							}

							array<hol_term*> new_definitions(2);
							intersect<built_in_predicates>(new_definitions, new_definition, conjunct);
							free(*new_definition); if (new_definition->reference_count == 0) free(new_definition);
							if (!dst.ensure_capacity(dst.length + new_definitions.length)) {
								for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
								free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
								free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
								free(*set_var); if (set_var->reference_count == 0) free(set_var);
								return false;
							}
							for (hol_term* new_definition : new_definitions) {
								if (ConjunctIndex >= 0) {
									dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
											make_array_view((hol_term**) nullptr, 0), make_appended_array_view(make_repeated_array_view(conjunct, ConjunctIndex), new_definition), make_array_view((hol_term**) nullptr, 0)));
									if (dst[dst.length] == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : new_definitions) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										return false;
									}
									conjunct->reference_count += 1 + ConjunctIndex;
								} else if (ConjunctIndex < 0) {
									unsigned int index = (unsigned int) (-ConjunctIndex) - 1;
									dst[dst.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, conjunct,
											make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0), make_prepended_array_view(new_definition, make_repeated_array_view(conjunct, index))));
									if (dst[dst.length] == nullptr) {
										for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
										for (hol_term* term : new_definitions) { free(*term); if (term->reference_count == 0) free(term); }
										free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
										free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
										free(*set_var); if (set_var->reference_count == 0) free(set_var);
										return false;
									}
									conjunct->reference_count += 1 + index;
								}
								new_definition->reference_count++;
								dst.length++;
							}
							for (hol_term* term : new_definitions) { free(*term); if (term->reference_count == 0) free(term); }
						}
					}
					for (hol_term* term : left_intersections) { free(*term); if (term->reference_count == 0) free(term); }
				}

				if (can_keep_set_variable) {
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
						fprintf(stderr, "invert_select_arg_without_head_predicative_and_factor ERROR: Unsupported value of `ConjunctIndex`.\n");
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						return false;
					}

					if (wide_scope_marker && element_narrow_scope) {
						hol_term* temp = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), hol_term::new_any_right(new_second_head));
						if (temp == nullptr) {
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							free(*new_second_head); free(new_second_head);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							return false;
						}
						new_second_head = temp;
					}

					hol_term* new_right_conjunct;
					hol_term* excluded_quantifier;
					if (right->type == hol_term_type::EXISTS) {
						hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_exists(element_variable, &HOL_ANY));
						if (excluded_quantifier == nullptr) {
							free(*new_second_head); free(new_second_head);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							return false;
						}
						HOL_ANY.reference_count++;
						new_right_conjunct = hol_term::new_any_right(hol_term::new_exists(element_variable, hol_term::new_and(hol_term::new_apply(set_var, element_var), hol_term::new_any_right(new_second_head, &excluded_quantifier, 1))));
					} else if (right->type == hol_term_type::FOR_ALL) {
						hol_term* excluded_quantifier = hol_term::new_any(hol_term::new_for_all(element_variable, &HOL_ANY));
						if (excluded_quantifier == nullptr) {
							free(*new_second_head); free(new_second_head);
							free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
							free(*set_var); if (set_var->reference_count == 0) free(set_var);
							free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
							return false;
						}
						HOL_ANY.reference_count++;
						new_right_conjunct = hol_term::new_any_right(hol_term::new_for_all(element_variable, hol_term::new_if_then(hol_term::new_apply(set_var, element_var), hol_term::new_any_right(new_second_head, &excluded_quantifier, 1))));
					} else {
						excluded_quantifier = nullptr;
						new_right_conjunct = hol_term::new_and(hol_term::new_apply(set_var, element_var), hol_term::new_any_right(new_second_head));
					}
					if (new_right_conjunct == nullptr) {
						free(*new_second_head); free(new_second_head);
						free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						if (excluded_quantifier != nullptr) { free(*excluded_quantifier); free(excluded_quantifier); }
						return false;
					}
					element_var->reference_count++;
					set_var->reference_count++;

					dst[dst.length] = substitute_head<any_node_position::NONE>(second_head, right, new_right_conjunct);
					if (dst[dst.length] == nullptr) {
						free(*new_right_conjunct); free(new_right_conjunct);
						free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
						free(*set_var); if (set_var->reference_count == 0) free(set_var);
						free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);
						return false;
					}
					dst.length++;
				}
				free(*arg_term); if (arg_term->reference_count == 0) free(arg_term);
				free(*set_var); if (set_var->reference_count == 0) free(set_var);
				free(*conjunct); if (conjunct->reference_count == 0) free(conjunct);

				return (dst.length > 0);
			}
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
		[input_tense_predicates,output_tense_predicates](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable)
		{
			unsigned int predicate_variable;
			if (second_head->type == hol_term_type::ANY) {
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
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable)
		{
			unsigned int predicate_variable;
			if (second_head->type == hol_term_type::ANY) {
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
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable)
		{
			if (!dst.ensure_capacity(dst.length + 1))
				return false;

			dst[dst.length] = hol_term::new_not(second_head);
			if (dst[dst.length] == nullptr)
				return false;
			second_head->reference_count++;
			dst.length++;
			return true;
		});
}

inline void find_root(hol_term* src, hol_term*& head, head_index& predicate_index, bool& negated) {
	head = src;
}

inline bool invert_predicate_only(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second)
{
	/* we expect the head of `first` to be at the root */
	hol_term* head; head_index predicate_index; bool negated;
	find_head<built_in_predicates>(first, head, predicate_index, negated);
	if (head == nullptr) return false;

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_root, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			array<hol_term*> new_heads(4);
			if (first_head->type == hol_term_type::ANY) {
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

			if (!dst.ensure_capacity(dst.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}

			hol_term* any_wide_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
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
				hol_term* excluded_trees[2];
				excluded_trees[0] = excluded_universal;
				excluded_trees[1] = duplicate_wide_scopes;
				dst[dst.length] = hol_term::new_any_right(
						hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), new_head),
						excluded_trees, array_length(excluded_trees));
				if (dst[dst.length] == nullptr) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
					free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
					return false;
				}
				excluded_universal->reference_count++;
				duplicate_wide_scopes->reference_count++;
				new_head->reference_count++;
				dst.length++;
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
	hol_term* head; head_index predicate_index; bool negated;
	find_head<built_in_predicates>(first, head, predicate_index, negated);
	if (head == nullptr) return false;

	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_root, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			array<hol_term*> new_heads(4);
			if (first_head->type == hol_term_type::ANY) {
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

				hol_term* head_conjunct = hol_term::new_any(
					hol_term::new_apply(
							hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), head_var),
					excluded_trees, excluded_tree_count);
				if (head_conjunct == nullptr) {
					free(*excluded_trees[0]); free(excluded_trees[0]);
					free(*excluded_trees[1]); free(excluded_trees[1]);
					free(*excluded_trees[2]); free(excluded_trees[2]);
					free(*excluded_trees[3]); free(excluded_trees[3]);
					return false;
				}
				hol_non_head_constants<built_in_predicates>::increment_terms();

				hol_term* predicate = hol_term::new_apply(second_head, head_var);
				if (predicate == nullptr) {
					free(*head_conjunct); free(head_conjunct);
					return false;
				}
				head_var->reference_count++;

				new_heads[new_heads.length] = hol_term::new_exists(predicate_variable, hol_term::new_any_array(hol_term_type::AND, head_conjunct,
						make_array_view(&predicate, 1), make_array_view((hol_term**) nullptr, 0), make_array_view((hol_term**) nullptr, 0)));
				if (new_heads[new_heads.length] == nullptr) {
					free(*predicate); free(predicate);
					free(*head_conjunct); free(head_conjunct);
					return false;
				}
				new_heads.length++;

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

				hol_term* first_predicate;
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

			if (!dst.ensure_capacity(dst.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}

			hol_term* any_wide_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
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
				hol_term* excluded_trees[2];
				excluded_trees[0] = excluded_universal;
				excluded_trees[1] = duplicate_wide_scopes;
				dst[dst.length] = hol_term::new_any_right(
						hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), new_head),
						excluded_trees, array_length(excluded_trees));
				if (dst[dst.length] == nullptr) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
					free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
					return false;
				}
				excluded_universal->reference_count++;
				duplicate_wide_scopes->reference_count++;
				new_head->reference_count++;
				dst.length++;
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
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_root, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			unsigned int predicate_variable;
			if (first_head->type == hol_term_type::ANY) {
				predicate_variable = max_variable + 1;
			} else if (first_head->type == hol_term_type::EXISTS) {
				predicate_variable = first_head->quantifier.variable;
			} else {
				fprintf(stderr, "invert_predicate_and_tense ERROR: Expected `first_head` to be an existentially quantified conjunction.\n");
				return false;
			}

			hol_term* var_term = hol_term::new_variable(predicate_variable);
			if (var_term == nullptr)
				return false;

			hol_term* expected_head = hol_term::new_exists(predicate_variable, hol_term::new_and(
				hol_term::new_apply(
					hol_term::new_any(nullptr, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()), var_term),
				hol_term::new_apply(hol_term::new_any_constant(make_array_view(PAST_OR_PRESENT, array_length(PAST_OR_PRESENT))), var_term)));
			if (expected_head == nullptr) {
				free(*var_term); free(var_term);
				return false;
			}
			hol_non_head_constants<built_in_predicates>::increment_terms();
			var_term->reference_count += 2 - 1;

			hol_term* any_wide_scope = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), &HOL_ANY));
			if (any_wide_scope == nullptr) {
				free(*expected_head); free(expected_head);
				return false;
			}
			HOL_ANY.reference_count++;

			hol_term* duplicate_wide_scopes = hol_term::new_any_right(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), any_wide_scope));
			if (duplicate_wide_scopes == nullptr) {
				free(*expected_head); free(expected_head);
				free(*any_wide_scope); free(any_wide_scope);
				return false;
			}

			hol_term* excluded_universal = hol_term::new_any_right(hol_term::new_any_quantifier(hol_quantifier_type::FOR_ALL, any_wide_scope));
			if (excluded_universal == nullptr) {
				free(*expected_head); free(expected_head);
				free(*duplicate_wide_scopes); free(duplicate_wide_scopes);
				return false;
			}
			any_wide_scope->reference_count++;

			array<hol_term*> new_heads(2);
			intersect<built_in_predicates>(new_heads, expected_head, second_head);
			free(*expected_head); if (expected_head->reference_count == 0) free(expected_head);
			if (!dst.ensure_capacity(dst.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}

			for (hol_term* new_head : new_heads) {
				hol_term* excluded_trees[2];
				excluded_trees[0] = excluded_universal;
				excluded_trees[1] = duplicate_wide_scopes;
				dst[dst.length] = hol_term::new_any_right(
						hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), new_head),
						excluded_trees, array_length(excluded_trees));
				if (dst[dst.length] == nullptr) {
					for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
					free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
					free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
					return false;
				}
				excluded_universal->reference_count++;
				duplicate_wide_scopes->reference_count++;
				new_head->reference_count++;
				dst.length++;
			}
			for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
			free(*excluded_universal); if (excluded_universal->reference_count == 0) free(excluded_universal);
			free(*duplicate_wide_scopes); if (duplicate_wide_scopes->reference_count == 0) free(duplicate_wide_scopes);
			return (dst.length > 0);
		});
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
		[lambda_variable](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			array<hol_term*> new_heads(2);
			select_predicate_in_set<false, false>(new_heads, first_head, second_head, max_variable, lambda_variable);
			if (new_heads.length == 0)
				return false;

			if (!dst.ensure_capacity(dst.length + new_heads.length)) {
				for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
				return false;
			}
			for (hol_term* new_head : new_heads) {
				hol_term* right = new_head->quantifier.operand->array.operands[new_head->quantifier.operand->array.length - 1];
				if (right->type == hol_term_type::UNARY_APPLICATION && right->binary.left->type == hol_term_type::CONSTANT && right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
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
					continue;
				}

				if (inner_right->type == hol_term_type::UNARY_APPLICATION && inner_right->binary.left->type == hol_term_type::CONSTANT && inner_right->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
					dst[dst.length++] = new_head;
					new_head->reference_count++;
					continue;
				} else {
					dst[dst.length] = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::WIDE_SCOPE), new_head);
					if (dst[dst.length] == nullptr) {
						for (hol_term* term : new_heads) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					new_head->reference_count++;
					dst.length++;
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
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head_or_universal<built_in_predicates>, find_unary_application<(unsigned int) built_in_predicates::WIDE_SCOPE>, no_op(),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, head_index first_predicate_index, head_index second_predicate_index, hol_term*& conjunct, unsigned int& max_variable) {
			if (second_head->type == hol_term_type::UNARY_APPLICATION && second_head->binary.left->type == hol_term_type::CONSTANT && second_head->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
				if (!dst.ensure_capacity(dst.length + 1)) return false;
				dst[dst.length++] = second_head->binary.right;
				second_head->binary.right->reference_count++;
				return true;
			} else {
				return false;
			}
		});
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
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::ANY))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
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
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
		return false;
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
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, cnj))
		return false;
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
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::ANY))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
		return false;
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
	 || !intersect(flags.corr, first.flags.corr, corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
	 || !intersect(flags.aux, expected_aux, first.flags.aux))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::INDICATIVE))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::ANY))
		return false;
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
	 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
	 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
	 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
	 || !intersect(flags.mood, first.flags.mood, grammatical_mood::ANY))
		return false;
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
	case function_type::REQUIRE_WIDE_SCOPE:
		/* the forward application already ensures that `second` satisfies this requirement */
		if (!intersect(flags, first.flags, second.flags)) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_LEFT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_LEFT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<0>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_INVERSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_inverse(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG2>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG3>(inverse, inverse_count, flags, first.root, second.root);
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
	case function_type::PREDICATE_ONLY:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate_only(inverse, inverse_count, flags, first.root, second.root);
	case function_type::PREDICATE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate(inverse, inverse_count, flags, first.root, second.root);
	case function_type::PREDICATE_AND_TENSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_predicate_and_tense(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_PREDICATE_IN_SET:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_predicate_in_set(inverse, inverse_count, flags, first.root, second.root);
	case function_type::MARK_WIDE_SCOPE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_mark_wide_scope(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_SINGULAR:
		if ((second.flags.index_number != grammatical_num::SINGULAR && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_PLURAL:
		if ((second.flags.index_number != grammatical_num::PLURAL && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_SINGULAR:
		if ((second.flags.index_number != grammatical_num::SINGULAR && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::SINGULAR)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_PLURAL:
		if ((second.flags.index_number != grammatical_num::PLURAL && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::PLURAL)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::TRY_REMOVE_NUMBER:
		if ((second.flags.index_number != grammatical_num::NONE && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		flags.index_number = first.flags.index_number;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_CONCORD_SINGULAR:
		if ((second.flags.concord_number != grammatical_num::SINGULAR && second.flags.concord_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, grammatical_num::NONE)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_CONCORD_PLURAL:
		if ((second.flags.concord_number != grammatical_num::PLURAL && second.flags.concord_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, grammatical_num::NONE)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
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
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
			return false;
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
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.corr, first.flags.corr, correlator::ANY))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NO_CORRELATOR:
		if ((second.flags.corr != correlator::NONE && second.flags.corr != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood)
		 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
			return false;
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
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::ANY)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return flags.intersect_aux_or_subjunctive_or_inf_or_to_inf(first.flags.aux_or_subjunctive_or_inf_or_to_inf, second.flags.aux_or_subjunctive_or_inf_or_to_inf)
			&& intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NOT_CORRELATED:
		if ((second.flags.correlated_by != correlator::NONE && second.flags.correlated_by != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE)
		 || !intersect(flags.aux, first.flags.aux, second.flags.aux)
		 || !intersect(flags.mood, first.flags.mood, second.flags.mood))
			return false;
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

	Formula* first_logical_form = relabel_variables(first.root);
	if (first_logical_form == nullptr)
		return false;

	Formula* second_logical_form = relabel_variables(second.root);
	if (second_logical_form == nullptr) {
		free(*first_logical_form); if (first_logical_form->reference_count == 0) free(first_logical_form);
		return false;
	}

	bool result = is_subset<built_in_predicates>(first_logical_form, second_logical_form);
	free(*first_logical_form); if (first_logical_form->reference_count == 0) free(first_logical_form);
	free(*second_logical_form); if (second_logical_form->reference_count == 0) free(second_logical_form);
	return result;
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
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_and(hol_term::new_apply(hol_term::new_any_constant_except(), &HOL_ANY), &HOL_ANY));
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

	hol_term* predicate = intersection[0]->any_quantifier.operand->array.operands[0]->binary.left;
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
			fprintf(stderr, "get_constant WARNING: Unexpected formula type in interseciton.\n");
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
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_and(hol_term::new_apply(hol_term::new_constant(value), &HOL_ANY), &HOL_ANY));
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

bool set_tense(hol_term* src, hol_term*& dst, unsigned int value)
{
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_and(&HOL_ANY, hol_term::new_apply(hol_term::new_constant(value), &HOL_ANY)));
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
	hol_term* term = hol_term::new_any_quantifier(hol_quantifier_type::EXISTS, hol_term::new_and(hol_term::new_apply(hol_term::new_any_constant_except(make_array_view(values, count)), &HOL_ANY), &HOL_ANY));
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
	return set_number(*exp.root, *set.root, value);
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
	return set_uint_list(*exp.root, *set.root, value);
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
			return result;
		}
	}

	if (pos == POS_VERB) {
		if (words.length != 1)
			return false;

		bool contains;
		const array<inflected_verb>& forms = morphology_parser.inflected_verbs.get(words[0], contains);
		if (!contains) return true;

		sequence root(NULL, 0); root = words;
		flagged_logical_form<Formula> marked_logical_form = logical_form;
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

			root[0] = form.root;
			if (!emit_root(root, words, marked_logical_form)) {
				free(marked_logical_form);
				free(root); return false;
			}
			free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
			marked_logical_form.root = logical_form.root;
			logical_form.root->reference_count++;
		}
		free(root);
		return true;
	} else if (pos == POS_NOUN) {
		if (words.length != 1)
			return false;

		bool contains;
		const array<inflected_noun>& forms = morphology_parser.inflected_nouns.get(words[0], contains);
		if (!contains) return true;

		sequence root(NULL, 0); root = words;
		flagged_logical_form<Formula> marked_logical_form = logical_form;
		for (const inflected_noun& form : forms) {
			if (!intersect(marked_logical_form.flags.index_number, logical_form.flags.index_number, form.number))
				continue;

			root[0] = form.root;
			if (!emit_root(root, words, marked_logical_form)) {
				free(marked_logical_form);
				free(root); return false;
			}
			free(*marked_logical_form.root); if (marked_logical_form.root->reference_count == 0) free(marked_logical_form.root);
			marked_logical_form.root = logical_form.root;
			logical_form.root->reference_count++;
		}
		free(root);
		return true;
	} else if (pos == POS_ADJECTIVE) {
		/* TODO: implement this */
		return emit_root(words, {nullptr, 0}, logical_form);
	} else if (pos == POS_ADVERB) {
		/* TODO: implement this */
		return emit_root(words, {nullptr, 0}, logical_form);
	} else {
		return emit_root(words, {nullptr, 0}, logical_form);
	}
}

#endif /* HDP_PARSER_H_ */
