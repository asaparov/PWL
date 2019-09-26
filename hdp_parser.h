#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include "array_view.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

enum class grammatical_num : uint_fast8_t {
	SINGULAR,
	PLURAL,
	ANY,
	NONE
};

enum class grammatical_conjunction : uint_fast8_t {
	NONE = 0,
	THAT,
	IF,
	WHETHER,
	BECAUSE,
	FOR,
	ANY
};

enum class grammatical_flag : uint_fast8_t {
	IS_ADJUNCT,
	NULLABLE_SUBJECT,
	SUBORDINATE,
	PREPOSITION,
	PARTICLE,
	REQ_AUX,
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
static inline bool print_feature(const char* feature_name, bool& first, grammatical_num number, Stream& out) {
	switch (number) {
	case grammatical_num::SINGULAR:
		first = false; return core::print(feature_name, out) && core::print(":sg", out);
	case grammatical_num::PLURAL:
		first = false; return core::print(feature_name, out) && core::print(":pl", out);
	case grammatical_num::ANY:
		first = false; return core::print(feature_name, out) && core::print(":*", out);
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
		first = false; return core::print("cnj:that", out);
	case grammatical_conjunction::IF:
		first = false; return core::print("cnj:if", out);
	case grammatical_conjunction::WHETHER:
		first = false; return core::print("cnj:whether", out);
	case grammatical_conjunction::BECAUSE:
		first = false; return core::print("cnj:because", out);
	case grammatical_conjunction::FOR:
		first = false; return core::print("cnj:for", out);
	case grammatical_conjunction::ANY:
		first = false; return core::print("cnj:*", out);
	case grammatical_conjunction::NONE:
		return true;
	}
	fprintf(stderr, "print_feature ERROR: Unrecognized grammatical_conjunction.\n");
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
	case grammatical_flag::REQ_AUX:
		return core::print("requires_aux", out);
	case grammatical_flag::COUNT: break;
	}
	fprintf(stderr, "print_feature_name ERROR: Unrecognized grammatical_flag.\n");
	return false;
}

template<typename Stream>
static inline bool print_feature(bool& first, grammatical_flag flag, grammatical_flag_value value, Stream& out) {
	switch (value) {
	case grammatical_flag_value::TRUE:
		first = false; return print_feature_name(flag, out);
	case grammatical_flag_value::ANY:
		first = false; return print_feature_name(flag, out) && print('*', out);
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
		first = false; return core::print(prefix, out) && core::print(":both", out);
	case correlator::EITHER:
		first = false; return core::print(prefix, out) && core::print(":either", out);
	case correlator::NEITHER:
		first = false; return core::print(prefix, out) && core::print(":neither", out);
	case correlator::ANY:
		first = false; return core::print(prefix, out) && core::print(":*", out);
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

struct grammatical_flags {
	grammatical_num index_number;
	grammatical_num concord_number;
	grammatical_conjunction cnj;
	grammatical_flag_value flags[(uint_fast8_t) grammatical_flag::COUNT] = { grammatical_flag_value::FALSE };
	correlator corr;
	correlator correlated_by;

	grammatical_flags() :
			index_number(grammatical_num::NONE), concord_number(grammatical_num::NONE),
			cnj(grammatical_conjunction::NONE), corr(correlator::NONE), correlated_by(correlator::NONE) { }

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
			 ^ (127 * default_hash(key.correlated_by))
			 ^ (8191 * default_hash(key.flags, (uint_fast8_t) grammatical_flag::COUNT));
	}

	static inline void swap(grammatical_flags& first, grammatical_flags& second) {
		core::swap(first.index_number, second.index_number);
		core::swap(first.concord_number, second.concord_number);
		core::swap(first.cnj, second.cnj);
		core::swap(first.corr, second.corr);
		core::swap(first.correlated_by, second.correlated_by);
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			core::swap(first.flags[i], second.flags[i]);
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
	 || !intersect(dst.correlated_by, first.correlated_by, second.correlated_by))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(dst.flags[i], first.flags[i], second.flags[i])) return false;
	return true;
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
		CONSTANT
	};

	enum class function_type {
		EMPTY = 0,
		IDENTITY,
		SELECT_RIGHT_CONJUNCT,
		REMOVE_RIGHT_CONJUNCT,
		REQUIRE_NO_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE,
		REQUIRE_LEFT_PREDICATE_INVERSE_OWN,
		REMOVE_INVERSE,
		SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,
		SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,
		SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR,

		/* functions that modify grammatical features */
		ADD_SINGULAR,
		ADD_PLURAL,
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
		ADD_REQ_AUX,
		REMOVE_REQ_AUX,
		REQUIRE_NO_REQ_AUX,
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
		REQUIRE_NOT_CORRELATED
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
	{feature::CONSTANT, "constant"}
};

template<typename Formula>
const static_pair<typename flagged_logical_form<Formula>::function_type, const char*> flagged_logical_form<Formula>::FUNCTION_NAMES[] = {
	{function_type::EMPTY, "empty"},
	{function_type::IDENTITY, "identity"},
	{function_type::SELECT_RIGHT_CONJUNCT, "select_right_conjunct"},
	{function_type::REMOVE_RIGHT_CONJUNCT, "remove_right_conjunct"},
	{function_type::REQUIRE_NO_INVERSE, "require_no_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE, "require_left_predicate_inverse"},
	{function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN, "require_left_predicate_inverse_own"},
	{function_type::REMOVE_INVERSE, "remove_inverse"},
	{function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg1_without_head_predicative_and_factor"},
	{function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg2_without_head_predicative_and_factor"},
	{function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR, "select_right_arg3_without_head_predicative_and_factor"},
	{function_type::ADD_SINGULAR, "add_singular"},
	{function_type::ADD_PLURAL, "add_plural"},
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
	{function_type::ADD_REQ_AUX, "add_req_aux"},
	{function_type::REMOVE_REQ_AUX, "remove_req_aux"},
	{function_type::REQUIRE_NO_REQ_AUX, "require_no_req_aux"},
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
	{function_type::REQUIRE_NOT_CORRELATED, "require_not_correlated"}
};

inline void initialize_any(grammatical_flags& flags) {
	flags.concord_number = grammatical_num::ANY;
	flags.index_number = grammatical_num::ANY;
	flags.cnj = grammatical_conjunction::ANY;
	flags.corr = correlator::ANY;
	flags.correlated_by = correlator::ANY;
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
	 || first.correlated_by != second.correlated_by)
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (first.flags[i] != second.flags[i]) return false;
	return true;
}

template<typename Formula>
inline bool operator == (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
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
	 || first.correlated_by != second.correlated_by)
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

template<typename Stream>
inline bool print(const grammatical_flags& flags, Stream& out)
{
	bool first = true;
	if (!print('[', out)
	 || (first && !print(',', out)) || !print_feature("index", first, flags.index_number, out)
	 || (first && !print(',', out)) || !print_feature("concord", first, flags.concord_number, out)
	 || (first && !print(',', out)) || !print_feature(first, flags.cnj, out)
	 || (first && !print(',', out)) || !print_feature("corr:", first, flags.corr, out)
	 || (first && !print(',', out)) || !print_feature("corr_by:", first, flags.correlated_by, out))
		return false;

	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if ((first && !print(',', out)) || !print_feature(first, (grammatical_flag) i, flags.flags[i], out)) return false;
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
	const array<array_map<sentence, Formula>>& data,
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

inline bool init(sequence& seq, const sentence& src)
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

	morphology morph;
	hdp_grammar_type<Formula> G;
	const string** reverse_name_map;

	hdp_parser(unsigned int unknown_id,
			hash_map<string, unsigned int>& names,
			const char* agid_filepath,
			const char* uncountable_filepath,
			const char* grammar_filepath) : reverse_name_map(nullptr)
	{
		if (!morph.initialize(names)
		 || !morphology_read(morph, names, agid_filepath, uncountable_filepath))
		{
			fprintf(stderr, "ERROR: Unable to initialize morphology model.\n");
			exit(EXIT_FAILURE);
		} else if (!read_grammar(G, names, grammar_filepath)) {
			fprintf(stderr, "ERROR: Unable to read grammar at '%s'.\n", grammar_filepath);
			exit(EXIT_FAILURE);
		}
	}

	bool invert_name_map(const hash_map<string, unsigned int>& names) {
		if (reverse_name_map != NULL) free(reverse_name_map);
		reverse_name_map = invert(names);
		if (reverse_name_map == NULL) return false;
		return true;
	}

	bool train(
			const array<array_map<sentence, hol_term>>& data,
			hash_map<string, unsigned int>& names,
			unsigned int iteration_count)
	{
		const string** reverse_name_map = invert(names);
		const string** nonterminal_name_map = invert(G.nonterminal_names);
		if (reverse_name_map == NULL || nonterminal_name_map == NULL) {
			if (reverse_name_map != NULL) free(reverse_name_map);
			return false;
		}
		string_map_scribe terminal_printer = { reverse_name_map, names.table.size + 1 };
		string_map_scribe nonterminal_printer = { nonterminal_name_map, G.nonterminal_names.table.size + 1 };

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
	bool parse(const sentence& s, Formula** logical_forms,
			double* log_probabilities, unsigned int& parse_count,
			const TheoryType& T, array<sentence_token>& unrecognized)
	{
		static_assert(K > 0, "`K` must be at least 1.");
#if !defined(NDEBUG)
		if (reverse_name_map == NULL) {
			fprintf(stderr, "hdp_parser.parse ERROR: `hdp_parser.invert_name_map` must be called before `hdp_parser.parse`.\n");
			return false;
		}
#endif

		logical_form_type logical_form; /* TODO: initialize the features */
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

	bool add_definition(const sentence& s, const Formula* definition, unsigned int new_constant)
	{
	}

	template<typename Printer>
	constexpr Printer& get_printer(Printer& constant_printer) const {
		return constant_printer;
	}

private:
	void cleanup(const array<array_map<sentence, hol_term>>& data,
			const string** reverse_name_map, const string** nonterminal_name_map,
			syntax_node<logical_form_type>*** syntax, unsigned int* order)
	{
		if (syntax != NULL) {
			for (unsigned int k = 0; k < data.length; k++) {
				if (syntax[k] == NULL) continue;
				for (unsigned int l = 0; l < data[k].size; l++)
					if (syntax[k][l] != NULL) free(syntax[k][l]);
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
	hol_term* constants[5];

	hol_non_head_constants() {
		constants[0] = hol_term::new_constant((unsigned int) BuiltInPredicates::UNKNOWN);
		constants[1] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG1);
		constants[2] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG2);
		constants[3] = hol_term::new_constant((unsigned int) BuiltInPredicates::ARG3);
		constants[4] = hol_term::new_constant((unsigned int) BuiltInPredicates::SIZE);
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
		|| constant == (unsigned int) BuiltInPredicates::SIZE;
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
	} else if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.right->type == hol_term_type::VARIABLE
	 && src->binary.right->variable == scope_variable
	 && src->binary.left->type == hol_term_type::ANY)
	{
		return src->binary.left;
	}
	return 0;
}

template<typename BuiltInPredicates>
inline void find_head(
		hol_term* src, hol_term*& head,
		int& predicate_index,
		bool& negated)
{
	if (src->type == hol_term_type::ANY) {
		if (src->any.included == nullptr) {
			head = src;
			negated = false;
		} else {
			find_head<BuiltInPredicates>(src->any.included, head, predicate_index, negated);
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
		if (not_head) return;

		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_constant((unsigned int) BuiltInPredicates::ARG2), &HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) return;

		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_constant((unsigned int) BuiltInPredicates::ARG3), &HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) return;

		if (head->quantifier.operand->type == hol_term_type::ANY)
			negated = false;

		/* find the predicate */
		wildcard.any.included = hol_term::new_apply(
				hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<BuiltInPredicates>::get_terms(), hol_non_head_constants<BuiltInPredicates>::count()),
				hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		hol_term* operand = head->quantifier.operand;

		hol_term* predicate = nullptr;
		if (operand->type == hol_term_type::ANY_ARRAY) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(src->any_array.left, head_variable);
			if (predicate != nullptr) {
				predicate_index = 0;
			} else if (has_intersection<BuiltInPredicates>(src->any_array.left, &wildcard)) {
				predicate_index = 0;
			}

			predicate = get_predicate_of_literal<BuiltInPredicates>(src->any_array.right, head_variable);
			if (predicate != nullptr) {
				predicate_index = -1;
			} else if (has_intersection<BuiltInPredicates>(src->any_array.right, &wildcard)) {
				predicate_index = -1;
			}

			if (src->any_array.any != nullptr) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(src->any_array.any, head_variable);
				if (predicate != nullptr) {
					predicate_index = INT_MAX;
				} else if (has_intersection<BuiltInPredicates>(src->any_array.any, &wildcard)) {
					predicate_index = INT_MAX;
				}
			}
		} else if (operand->type == hol_term_type::AND) {
			predicate_index = operand->array.length;
			for (unsigned int i = 0; i < operand->array.length; i++) {
				predicate = get_predicate_of_literal<BuiltInPredicates>(src->array.operands[i], head_variable);
				if (predicate != nullptr) {
					predicate_index = i;
					break;
				} else if (has_intersection<BuiltInPredicates>(src->array.operands[i], &wildcard)) {
					predicate_index = i;
				}
			}
		} else {
			return;
		}

		if (predicate == nullptr) {
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
struct predicative_head_finder {
	unsigned int lambda_variable;

	predicative_head_finder(unsigned int lambda_variable) : lambda_variable(lambda_variable) { }

	void operator () (
			hol_term* src, hol_term*& head,
			int& predicate_index,
			bool& negated)
	{
		if (src->type == hol_term_type::ANY) {
			if (src->any.included == nullptr) {
				head = src;
				negated = false;
			} else {
				operator () (src->any.included, head, predicate_index, negated);
			}
		} else if (src->type == hol_term_type::UNARY_APPLICATION) {
			hol_term* lambda_var = hol_term::new_variable(lambda_variable);
			if (lambda_var == nullptr) return;
			if (!has_intersection<BuiltInPredicates>(lambda_var, src))
				return;

			head = src;
			negated = false;
		} else if (src->type == hol_term_type::NOT) {
			operator () (src->unary.operand, head, predicate_index, negated);
			negated = true;
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

template<bool Factor, typename MakeApplyHeadContext, typename MakeDstFunction, typename ApplyFunction>
inline bool apply_head(
		hol_term* src, hol_term*& dst,
		array<unsigned int>& dst_variables,
		array<hol_term*>& siblings,
		unsigned int max_variable,
		MakeApplyHeadContext make_context,
		MakeDstFunction make_dst,
		ApplyFunction apply)
{
	if (!apply(src)) return false;

	/* check if the current scope is the head */
	hol_term* head; int predicate_index; bool negated;
	find_head<built_in_predicates>(src, head, predicate_index, negated);
	if (head != nullptr) {
		unsigned int head_variable = 0;
		if (head->type == hol_term_type::ANY) {
			max_bound_variable(*head->any.included, max_variable);
			head_variable = max_variable + 1;
		} else if (head->type == hol_term_type::EXISTS) {
			head_variable = head->quantifier.variable;
		}
		dst = make_dst(head, head_variable, predicate_index, negated);
		if (dst == nullptr) return false;

		get_free_variables(*dst, dst_variables);
		dst_variables.remove(dst_variables.index_of(head_variable));
		if (dst_variables.length > 1) insertion_sort(dst_variables);
		return true;
	}

	if (src->type == hol_term_type::AND || src->type == hol_term_type::OR) {
		if (!siblings.ensure_capacity(siblings.length + src->array.length)) return false;
		find_head<built_in_predicates>(src->array.operands[src->array.length - 1], head, predicate_index, negated);

		if (head != nullptr) {
			unsigned int head_variable = 0;
			if (head->type == hol_term_type::ANY) {
				max_bound_variable(*head->any.included, max_variable);
				head_variable = max_variable + 1;
			} else if (head->type == hol_term_type::EXISTS) {
				head_variable = head->quantifier.variable;
			}
			dst = make_dst(head, head_variable, predicate_index, negated);
			if (dst == nullptr) return false;

			get_free_variables(*dst, dst_variables);
			dst_variables.remove(dst_variables.index_of(head_variable));
			if (dst_variables.length > 1) insertion_sort(dst_variables);

			auto context = make_context(head);

			/* check if the other conjuncts are also part of the head scope (i.e. the head scope is a conjunction) */
			hol_term** new_dst = (hol_term**) malloc(sizeof(hol_term*) * src->array.length);
			if (new_dst == nullptr) {
				free(*dst);
				if (dst->reference_count == 0)
					free(dst);
				return false;
			}
			new_dst[src->array.length - 1] = dst;
			array<unsigned int> new_variables(1 << (core::log2(dst_variables.length == 0 ? 1 : dst_variables.length) + 1));
			new_variables.append(dst_variables.data, dst_variables.length);
			for (unsigned int i = 0; i + 1 < src->array.length; i++) {
				find_head<built_in_predicates>(src->array.operands[i], head, predicate_index, negated);
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
					new_dst[i] = make_dst(head, head_variable, predicate_index, negated);
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
					free(*dst);
					if (dst->reference_count == 0)
						free(dst);
					free(new_dst);
					return false;
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
							if (new_dst[new_start]->type == hol_term_type::ANY && new_dst[new_start]->any.excluded_subtree_count == 0 && new_dst[new_start]->any.excluded_tree_count == 0) {
								dst = new_dst[new_start];
							} else {
								dst = hol_term::new_any(new_dst[new_start], nullptr, 0, nullptr, 0);
							}
						} else if (src->type == hol_term_type::AND) {
							dst = hol_term::new_any(hol_term::new_and(make_array_view(new_dst + new_start, new_length)), nullptr, 0, nullptr, 0);
						} else {
							dst = hol_term::new_any(hol_term::new_or(make_array_view(new_dst + new_start, new_length)), nullptr, 0, nullptr, 0);
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
					free(*dst);
					if (dst->reference_count == 0)
						free(dst);
					return false;
				}
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
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IFF:
		dst = nullptr;
		return true;

	case hol_term_type::NOT:
		if (!apply_head<Factor>(src->unary.operand, new_formula, dst_variables, siblings, max_variable, make_context, make_dst, apply))
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
		if (!apply_head<Factor>(src->array.operands[src->array.length - 1], new_formula, dst_variables, siblings, max_variable, make_context, make_dst, apply))
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

			for (unsigned int i = 0; i + 1 < src->array.length; i++)
				siblings[siblings.length++] = src->array.operands[i];

			if (!prune_independent_siblings(dst_variables, siblings)) {
				free(*new_formula);
				if (new_formula->reference_count == 0)
					free(new_formula);
				free(new_dst);
				return false;
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
						if (new_dst[new_start]->type == hol_term_type::ANY && new_dst[new_start]->any.excluded_subtree_count == 0 && new_dst[new_start]->any.excluded_tree_count == 0) {
							dst = new_dst[new_start];
						} else {
							dst = hol_term::new_any(new_dst[new_start], nullptr, 0, nullptr, 0);
						}
					} else if (src->type == hol_term_type::AND) {
						dst = hol_term::new_any(hol_term::new_and(make_array_view(new_dst + new_start, new_length)), nullptr, 0, nullptr, 0);
					} else {
						dst = hol_term::new_any(hol_term::new_or(make_array_view(new_dst + new_start, new_length)), nullptr, 0, nullptr, 0);
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
		}
		return true;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!apply_head<Factor>(src->quantifier.operand, new_formula, dst_variables, siblings, max(max_variable, src->quantifier.variable), make_context, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (!dst_variables.contains(src->quantifier.variable)) {
				if (new_formula->type == hol_term_type::ANY && new_formula->any.excluded_subtree_count == 0 && new_formula->any.excluded_tree_count == 0) {
					dst = new_formula;
				} else {
					dst = hol_term::new_any(new_formula, nullptr, 0, nullptr, 0);
				}
			} else {
				if (new_formula == src->quantifier.operand) {
					free(*new_formula);
					dst = src;
					dst->reference_count++;
				} else {
					if (src->type == hol_term_type::FOR_ALL)
						dst = hol_term::new_for_all(src->quantifier.variable, new_formula);
					else if (src->type == hol_term_type::EXISTS)
						dst = hol_term::new_exists(src->quantifier.variable, new_formula);
					else if (src->type == hol_term_type::LAMBDA)
						dst = hol_term::new_lambda(src->quantifier.variable, new_formula);
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
		if (!apply_head<Factor>(src->binary.right, new_formula, dst_variables, siblings, max_variable, make_context, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (new_formula == src->binary.right) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				dst = hol_term::new_if_then(src->binary.left, new_formula);
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
		if (src->any.included != nullptr) {
			if (!apply_head<Factor>(src->any.included, new_formula, dst_variables, siblings, max_variable, make_context, make_dst, apply))
				return false;
			if (new_formula != nullptr) {
				if (new_formula == src->any.included) {
					free(*new_formula);
					dst = src;
					dst->reference_count++;
				} else {
					dst = hol_term::new_any(new_formula, src->any.excluded_subtrees, src->any.excluded_subtree_count, src->any.excluded_trees, src->any.excluded_tree_count);
					if (dst == nullptr) {
						free(*new_formula);
						if (new_formula->reference_count == 0)
							free(new_formula);
						return false;
					}
					for (unsigned int i = 0; i < src->any.excluded_subtree_count; i++)
						src->any.excluded_subtrees[i]->reference_count++;
					for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
						src->any.excluded_trees[i]->reference_count++;
				}
			}
		} else {
			dst = nullptr;
		}
		return true;

	case hol_term_type::ANY_ARRAY:
		if (!apply_head<Factor>(src->any_array.right, new_formula, dst_variables, siblings, max_variable, make_context, make_dst, apply))
			return false;
		if (new_formula != nullptr) {
			if (new_formula == src->any_array.right) {
				free(*new_formula);
				dst = src;
				dst->reference_count++;
			} else {
				if (!new_hol_term(dst)) {
					free(*new_formula);
					if (new_formula->reference_count == 0)
						free(new_formula);
					return false;
				}
				dst->any_array.oper = src->any_array.oper;
				dst->any_array.all = src->any_array.all;
				dst->any_array.any = src->any_array.any;
				dst->any_array.left = src->any_array.left;
				dst->any_array.right = new_formula;
				dst->any_array.all->reference_count++;
				dst->any_array.left->reference_count++;
				if (dst->any_array.any != nullptr)
					dst->any_array.any->reference_count++;
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
			if (ConjunctIndex == 0)
				return operand->any_array.left;
			else if (ConjunctIndex == -1)
				return operand->any_array.right;
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
	else fprintf(stderr, "select_conjunct ERROR: Unexpected hol_term_type.\n");
#endif
	return nullptr;
}

inline bool compatible_context(hol_term* first, hol_term* second) {
	return first == second || *first == *second;
}

inline constexpr bool no_op(const hol_term* src) { return true; }

template<int_fast8_t ConjunctIndex>
inline bool select_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	return apply_head<false>(src, dst, dst_variables, siblings, 0, get_conjunct_context<ConjunctIndex>,
			[](hol_term* head, unsigned int head_variable, int predicate_index, bool negated)
			{
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
					hol_term* excluded_subtrees[3];
					excluded_subtrees[0] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), hol_term::new_variable(head_variable));
					excluded_subtrees[1] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), hol_term::new_variable(head_variable));
					excluded_subtrees[2] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), hol_term::new_variable(head_variable));
					if (excluded_subtrees[0] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[1] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[2] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[0] == nullptr || excluded_subtrees[1] == nullptr || excluded_subtrees[2] == nullptr) {
						if (excluded_subtrees[0] != nullptr) { free(*excluded_subtrees[0]); free(excluded_subtrees[0]); }
						if (excluded_subtrees[1] != nullptr) { free(*excluded_subtrees[1]); free(excluded_subtrees[1]); }
						return (hol_term*) nullptr;
					}

					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(
						hol_term::new_any(hol_term::new_apply(
								hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()),
								hol_term::new_variable(head_variable)),
							excluded_subtrees, 3, nullptr, 0),
						hol_term::new_any(nullptr, excluded_subtrees, 3, nullptr, 0)
					));
					if (dst == nullptr) {
						free(*excluded_subtrees[0]); free(excluded_subtrees[0]);
						free(*excluded_subtrees[1]); free(excluded_subtrees[1]);
						free(*excluded_subtrees[2]); free(excluded_subtrees[2]);
						return (hol_term*) nullptr;
					}
					excluded_subtrees[0]->reference_count += 2 - 1;
					excluded_subtrees[1]->reference_count += 2 - 1;
					excluded_subtrees[2]->reference_count += 2 - 1;

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "select_conjunct WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						free(*dst); if (dst->reference_count != 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count != 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
					return dst;
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
#if !defined(NDEBUG)
					if (ConjunctIndex == predicate_index)
						fprintf(stderr, "select_conjunct WARNING: The index of the conjunct and the index of the predicate are the same.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					hol_term* conjunct;
					if (ConjunctIndex == 0)
						conjunct = operand->any_array.left;
					else if (ConjunctIndex == -1)
						conjunct = operand->any_array.right;
					hol_term* predicate;
					if (predicate_index == 0)
						predicate = operand->any_array.left;
					else if (predicate_index == -1)
						predicate = operand->any_array.right;
					else if (predicate_index == INT_MAX)
						predicate = operand->any_array.any;
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
					hol_term* conjunct = operand->array.operands[(ConjunctIndex < 0) ? (operand->array.length + ConjunctIndex) : ConjunctIndex];
					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(operand->array.operands[predicate_index], conjunct));
					if (dst != nullptr) {
						operand->array.operands[predicate_index]->reference_count++;
						conjunct->reference_count++;
					}
					return dst;
				}
			}, no_op);
}

template<int_fast8_t ConjunctIndex>
inline bool remove_conjunct(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	return apply_head<false>(src, dst, dst_variables, siblings, 0, get_conjunct_context<ConjunctIndex>,
			[](hol_term* head, unsigned int head_variable, int predicate_index, bool negated)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY)) {
					hol_term* excluded_subtrees[3];
					excluded_subtrees[0] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), hol_term::new_variable(head_variable));
					excluded_subtrees[1] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), hol_term::new_variable(head_variable));
					excluded_subtrees[2] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), hol_term::new_variable(head_variable));
					if (excluded_subtrees[0] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[1] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[2] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[0] == nullptr || excluded_subtrees[1] == nullptr || excluded_subtrees[2] == nullptr) {
						if (excluded_subtrees[0] != nullptr) { free(*excluded_subtrees[0]); free(excluded_subtrees[0]); }
						if (excluded_subtrees[1] != nullptr) { free(*excluded_subtrees[1]); free(excluded_subtrees[1]); }
						return (hol_term*) nullptr;
					}

					hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_any(
						hol_term::new_apply(
								hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()),
								hol_term::new_variable(head_variable)),
						excluded_subtrees, 3, nullptr, 0));
					if (dst == nullptr) {
						free(*excluded_subtrees[0]); free(excluded_subtrees[0]);
						free(*excluded_subtrees[1]); free(excluded_subtrees[1]);
						free(*excluded_subtrees[2]); free(excluded_subtrees[2]);
						return (hol_term*) nullptr;
					}

					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, dst, head);
					if (intersection.length > 1) {
						fprintf(stderr, "remove_conjunct WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						free(*dst); if (dst->reference_count != 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count != 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunction;
					if (!new_hol_term(conjunction))
						return (hol_term*) nullptr;
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = operand->any_array.all;
					conjunction->any_array.all->reference_count++;
					if (ConjunctIndex == 0) {
						conjunction->any_array.left = operand->any_array.all;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.right = operand->any_array.right;
						conjunction->any_array.right->reference_count++;
						conjunction->any_array.any = operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					} else if (ConjunctIndex == -1) {
						conjunction->any_array.right = operand->any_array.all;
						conjunction->any_array.right->reference_count++;
						conjunction->any_array.left = operand->any_array.left;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.any = operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					}
#if !defined(NDEBUG)
					else fprintf(stderr, "remove_conjunct WARNING: Unsupported value of `ConjunctIndex`.\n");
#endif

					hol_term* dst = hol_term::new_exists(head_variable, conjunction);
					if (dst == nullptr) {
						free(*conjunction); free(conjunction);
						return (hol_term*) nullptr;
					}
					return dst;
				} else {
#if !defined(NDEBUG)
					if (head->type != hol_term_type::EXISTS || head->quantifier.operand->type != hol_term_type::AND)
						fprintf(stderr, "remove_conjunct WARNING: Expected `head` to be an existentially quantified conjunction.\n");
#endif
					hol_term* operand = head->quantifier.operand;
					unsigned int conjunct_index = (ConjunctIndex < 0) ? (head->array.length + ConjunctIndex) : ConjunctIndex;
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
					return new_dst;
				} else {
					return dst;
				}
			}, no_op);
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
	bool result = apply_head<false>(src, dst, dst_variables, siblings, 0, make_empty_context,
			[predicate, &predicates](hol_term* head, unsigned int head_variable, int predicate_index, bool negated)
			{
				if (predicates.ensure_capacity(predicates.size + 1))
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
					if (PredicateIndex == 0) {
						conjunction->any_array.left = current_predicate;
						conjunction->any_array.right = &HOL_ANY;
						conjunction->any_array.any = nullptr;
						HOL_ANY.reference_count += 2;
					} else if (PredicateIndex == -1) {
						conjunction->any_array.right = current_predicate;
						conjunction->any_array.left = &HOL_ANY;
						conjunction->any_array.any = nullptr;
						HOL_ANY.reference_count += 2;
					} else if (PredicateIndex == INT_FAST8_MAX) {
						conjunction->any_array.right = &HOL_ANY;
						conjunction->any_array.left = &HOL_ANY;
						conjunction->any_array.any = current_predicate;
						HOL_ANY.reference_count += 3;
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
						fprintf(stderr, "require_predicate WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						free(*dst); if (dst->reference_count != 0) free(dst);
						return (hol_term*) nullptr;
					}
					free(*dst); if (dst->reference_count != 0) free(dst);
					if (intersection.length == 0)
						return (hol_term*) nullptr;
					dst = intersection[0];
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;
					if (PredicateIndex != INT_FAST8_MAX && PredicateIndex != predicate_index)
						return (hol_term*) nullptr;

					hol_term* predicate = nullptr;
					if (predicate_index == 0)
						predicate = operand->any_array.left;
					else if (predicate_index == -1)
						predicate = operand->any_array.right;
					else if (predicate_index == INT_MAX)
						predicate = operand->any_array.any;
					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, predicate, current_predicate);
					if (intersection.length > 1) {
						fprintf(stderr, "require_predicate WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						return (hol_term*) nullptr;
					}
					if (intersection.length == 0)
						return (hol_term*) nullptr;

					if (intersection[0] == predicate) {
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						dst = head;
						head->reference_count++;
					} else {
						hol_term* conjunction;
						if (!new_hol_term(conjunction)) {
							for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
							return (hol_term*) nullptr;
						}
						conjunction->type = hol_term_type::ANY_ARRAY;
						conjunction->reference_count = 1;
						conjunction->any_array.oper = hol_term_type::AND;
						conjunction->any_array.all = operand->any_array.all;
						conjunction->any_array.all->reference_count++;
						if (predicate_index == 0) {
							conjunction->any_array.left = intersection[0];
							conjunction->any_array.right = operand->any_array.right;
							conjunction->any_array.right->reference_count++;
							conjunction->any_array.any = operand->any_array.any;
							if (conjunction->any_array.any != nullptr)
								conjunction->any_array.any->reference_count++;
						} else if (predicate_index == -1) {
							conjunction->any_array.right = intersection[0];
							conjunction->any_array.left = operand->any_array.left;
							conjunction->any_array.left->reference_count++;
							conjunction->any_array.any = operand->any_array.any;
							if (conjunction->any_array.any != nullptr)
								conjunction->any_array.any->reference_count++;
						} else if (predicate_index == INT_MAX) {
							conjunction->any_array.any = intersection[0];
							conjunction->any_array.left = operand->any_array.left;
							conjunction->any_array.left->reference_count++;
							conjunction->any_array.right = operand->any_array.right;
							conjunction->any_array.right->reference_count++;
						}

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
						if (predicate_index != expected_predicate_index)
							return (hol_term*) nullptr;
					}

					hol_term* predicate = operand->array.operands[predicate_index];
					array<hol_term*> intersection(2);
					intersect<built_in_predicates>(intersection, predicate, current_predicate);
					if (intersection.length > 1) {
						fprintf(stderr, "require_predicate WARNING: Expected intersection size to be 1.\n");
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						return (hol_term*) nullptr;
					}
					if (intersection.length == 0)
						return (hol_term*) nullptr;

					if (intersection[0] == predicate) {
						for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
						dst = head;
						head->reference_count++;
					} else {
						hol_term* conjunction;
						if (!new_hol_term(conjunction)) {
							for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
							return (hol_term*) nullptr;
						}
						conjunction->type = hol_term_type::AND;
						conjunction->reference_count = 1;
						conjunction->array.length = operand->array.length;
						conjunction->array.operands = (hol_term**) malloc(sizeof(hol_term*) * operand->array.length);
						if (conjunction->array.operands == nullptr) {
							fprintf(stderr, "require_predicate ERROR: Out of memory.\n");
							for (hol_term* term : intersection) { free(*term); if (term->reference_count != 0) free(term); }
							return (hol_term*) nullptr;
						}
						for (unsigned int i = 0; i < operand->array.length; i++) {
							if ((int) i == predicate_index)
								conjunction->array.operands[i] = intersection[0];
							else
								conjunction->array.operands[i] = operand->array.operands[i];
							conjunction->array.operands[i]->reference_count++;
						}

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
					return new_dst;
				} else {
					return dst;
				}
			}, no_op);
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
	excluded_constants[hol_non_head_constants<built_in_predicates>::count() + 1] = hol_term::new_constant((unsigned int) built_in_predicates::INVERSE);
	if (excluded_constants[hol_non_head_constants<built_in_predicates>::count() + 1] == nullptr) {
		free(excluded_constants);
		return false;
	}

	hol_term* predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, nullptr, 0, excluded_constants, hol_non_head_constants<built_in_predicates>::count() + 1),
			hol_term::new_variable(0));
	if (predicate == nullptr) {
		free(*excluded_constants[hol_non_head_constants<built_in_predicates>::count() + 1]);
		if (excluded_constants[hol_non_head_constants<built_in_predicates>::count() + 1]->reference_count == 0)
			free(excluded_constants[hol_non_head_constants<built_in_predicates>::count() + 1]);
		free(excluded_constants);
		return false;
	}

	bool result = require_predicate<INT_FAST8_MAX>(src, dst, predicate);
	free(*predicate); if (predicate->reference_count == 0) free(predicate);
	return result;
}

inline bool require_left_predicate_inverse(hol_term* src, hol_term*& dst)
{
	hol_term* predicate = hol_term::new_apply(
			hol_term::new_apply(
				hol_term::new_constant((unsigned int) built_in_predicates::INVERSE),
				hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count())),
			hol_term::new_variable(0));
	if (predicate == nullptr) return false;

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
	return apply_head<false>(src, dst, dst_variables, siblings, 0, make_empty_context,
			[](hol_term* head, unsigned int head_variable, int predicate_index, bool negated)
			{
				hol_term* dst;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY))
				{
					dst = head;
					dst->reference_count++;
				} else if (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY_ARRAY) {
					hol_term* operand = head->quantifier.operand;

					hol_term* predicate;
					if (predicate_index == 0)
						predicate = operand->any_array.left;
					else if (predicate_index == -1)
						predicate = operand->any_array.right;
					else if (predicate_index == INT_MAX)
						predicate = operand->any_array.any;

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
					if (predicate_index == 0) {
						conjunction->any_array.left = new_predicate;
						conjunction->any_array.right = operand->any_array.right;
						conjunction->any_array.right->reference_count++;
						conjunction->any_array.any = operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					} else if (predicate_index == -1) {
						conjunction->any_array.right = new_predicate;
						conjunction->any_array.left = operand->any_array.left;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.any = operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					} else if (predicate_index == INT_MAX) {
						conjunction->any_array.any = new_predicate;
						conjunction->any_array.left = operand->any_array.left;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.right = operand->any_array.right;
						conjunction->any_array.right->reference_count++;
					}

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
					hol_term* predicate = operand->array.operands[predicate_index];
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
						if ((int) i == predicate_index) {
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
					return new_dst;
				} else {
					return dst;
				}
			}, no_op);
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
			if (ConjunctIndex == 0)
				conjunct = operand->any_array.left;
			else if (ConjunctIndex == -1)
				conjunct = operand->any_array.right;
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
	else fprintf(stderr, "select_conjunct ERROR: Unexpected hol_term_type.\n");
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

	array_map<unsigned int, hol_term*> scopes(8);
	auto gather_scope = [&scopes](hol_term* term) {
		if (term->type == hol_term_type::FOR_ALL || term->type == hol_term_type::EXISTS || term->type == hol_term_type::LAMBDA) {
			if (!scopes.ensure_capacity(scopes.size + 1))
				return false;
			scopes.keys[scopes.size] = term->quantifier.variable;
			scopes.values[scopes.size++] = term;
		}
		return true;
	};

	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	bool result = apply_head<true>(src, dst, dst_variables, siblings, 0, get_arg_context<ConjunctIndex>,
			[&scopes,lambda_variable](hol_term* head, unsigned int head_variable, int predicate_index, bool negated)
			{
				hol_term* dst = nullptr;
				if (head->type == hol_term_type::ANY || (head->type == hol_term_type::EXISTS && head->quantifier.operand->type == hol_term_type::ANY))
				{
					unsigned int set_variable = lambda_variable + 1;
					unsigned int element_variable = lambda_variable + 2;
					/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
					if (negated)
						dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_not(hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0))));
					else dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0)));
					if (dst == nullptr) return (hol_term*) nullptr;
					HOL_ANY.reference_count++;
				} else if (head->type == hol_term_type::EXISTS) {
					hol_term* operand = head->quantifier.operand;

					hol_term* conjunct;
					if (operand->type == hol_term_type::ANY_ARRAY) {
						if (ConjunctIndex == 0)
							conjunct = operand->any_array.left;
						else if (ConjunctIndex == -1)
							conjunct = operand->any_array.right;
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
					hol_term* expected_conjunct = hol_term::new_any(hol_term::new_equals(expected_left_side, &HOL_ANY), nullptr, 0, nullptr, 0);
					if (expected_left_side == nullptr) {
						free(*expected_left_side); free(expected_left_side);
						return (hol_term*) nullptr;
					}
					expected_left_side->reference_count++;
					if (!has_intersection<built_in_predicates>(expected_conjunct, conjunct)) {
						free(*expected_conjunct); free(expected_conjunct);
						free(*expected_left_side); free(expected_left_side);
						return (hol_term*) nullptr;
					}
					free(*expected_conjunct); free(expected_conjunct);

					if (conjunct->type == hol_term_type::ANY) {
						unsigned int set_variable = lambda_variable + 1;
						unsigned int element_variable = lambda_variable + 2;
						/* TODO: this set isn't really the same as the result of applying the function to `ANY` */
						if (negated)
							dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_not(hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0))));
						else dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0)));
						if (dst == nullptr) return (hol_term*) nullptr;
						HOL_ANY.reference_count++;
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
							if (negated)
								dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_not(hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0))));
							else dst = hol_term::new_exists(set_variable, hol_term::new_and(&HOL_ANY, hol_term::new_any(hol_term::new_apply(hol_term::new_variable(lambda_variable), hol_term::new_variable(element_variable)), nullptr, 0, nullptr, 0)));
							if (dst == nullptr) return (hol_term*) nullptr;
							HOL_ANY.reference_count++;
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
							if (negated) {
								hol_term* new_quantified_term = hol_term::new_not(quantified_term);
								if (new_quantified_term == nullptr) {
									free(*quantified_term); free(quantified_term);
									return (hol_term*) nullptr;
								}
								quantified_term = new_quantified_term;
							}
							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, hol_term::new_equals(element_var, hol_term::new_constant(arg->constant)))), quantified_term));
							if (dst == nullptr) {
								free(*quantified_term); free(quantified_term);
								return (hol_term*) nullptr;
							}
							element_var->reference_count++;
							set_var->reference_count++;
						} else if (arg->type == hol_term_type::VARIABLE) {
							unsigned int element_variable = arg->variable;
							hol_term* scope = scopes.get(element_variable);
							hol_term* operand = scope->quantifier.operand;
							if (scope->type == hol_term_type::FOR_ALL && operand->type == hol_term_type::IF_THEN) {
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
									dst = hol_term::new_exists(set_variable, hol_term::new_and(
											hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, antecedent)),
											hol_term::new_for_all(element_variable, hol_term::new_if_then(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)))
										));
									if (dst == nullptr) {
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									antecedent->reference_count++;
									set_var->reference_count += 2;
									arg->reference_count += 2;
								}
							} else if (scope->type == hol_term_type::EXISTS && operand->type == hol_term_type::AND) {
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
									hol_term* set_definition;
									if (operand->array.length == 2) {
										set_definition = first_conjunct;
										set_definition->reference_count++;
									} else {
										set_definition = hol_term::new_and(make_array_view(operand->array.operands, operand->array.length - 1));
										if (set_definition == nullptr)
											return (hol_term*) nullptr;
										for (unsigned int i = 0; i < set_definition->array.length; i++)
											set_definition->array.operands[i]->reference_count++;
									}

									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr) {
										free(*set_definition); if (set_definition->reference_count == 0) free(set_definition);
										return (hol_term*) nullptr;
									}
									dst = hol_term::new_exists(set_variable, hol_term::new_and(
											hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, set_definition)),
											hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)))
										));
									if (dst == nullptr) {
										free(*set_definition); if (set_definition->reference_count == 0) free(set_definition);
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									set_var->reference_count += 2;
									arg->reference_count += 2;
								}
							} else if (scope->type == hol_term_type::EXISTS && operand->type == hol_term_type::ANY_ARRAY && (operand->any_array.oper == hol_term_type::AND || operand->any_array.oper == hol_term_type::ANY_ARRAY)) {
								hol_term* first_conjunct = operand->any_array.left;
								/* check if the scope is a relativized quantification over a set */
								if (first_conjunct->type == hol_term_type::UNARY_APPLICATION && first_conjunct->binary.right->type == hol_term_type::VARIABLE
								 && first_conjunct->binary.right->variable == element_variable && first_conjunct->binary.left->type == hol_term_type::VARIABLE)
								{
									dst = hol_term::new_apply(hol_term::new_variable(lambda_variable), arg);
									if (dst == nullptr)
										return (hol_term*) nullptr;
									arg->reference_count++;
								} else {
									hol_term* set_definition = hol_term::new_any_array(hol_term_type::AND, operand->any_array.all, operand->any_array.any, operand->any_array.left, operand->any_array.all);
									if (set_definition == nullptr)
										return (hol_term*) nullptr;
									operand->any_array.all->reference_count += 2;
									operand->any_array.left->reference_count++;
									if (operand->any_array.any != nullptr)
										operand->any_array.any->reference_count++;

									unsigned int set_variable = lambda_variable + 1;
									hol_term* set_var = hol_term::new_variable(set_variable);
									if (set_var == nullptr) {
										free(*set_definition); if (set_definition->reference_count == 0) free(set_definition);
										return (hol_term*) nullptr;
									}
									dst = hol_term::new_exists(set_variable, hol_term::new_and(
											hol_term::new_equals(set_var, hol_term::new_lambda(element_variable, set_definition)),
											hol_term::new_exists(element_variable, hol_term::new_and(
												hol_term::new_apply(set_var, arg),
												hol_term::new_apply(hol_term::new_variable(lambda_variable), arg)))
										));
									if (dst == nullptr) {
										free(*set_definition); if (set_definition->reference_count == 0) free(set_definition);
										free(*set_var); free(set_var);
										return (hol_term*) nullptr;
									}
									set_var->reference_count += 2;
									arg->reference_count += 2;
								}
							} else if (scope->type == hol_term_type::LAMBDA) {
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
						hol_term* expected_operand = hol_term::new_any_array(hol_term_type::AND, &HOL_ANY, nullptr, &HOL_ANY, equals_term);
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
							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable,
										hol_term::new_any_array(hol_term_type::AND, inner_operand->any_array.all, inner_operand->any_array.left, inner_operand->any_array.any, inner_operand->any_array.all))),
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
							inner_operand->any_array.all->reference_count += 2;
							inner_operand->any_array.left->reference_count++;
							if (inner_operand->any_array.any != nullptr)
								inner_operand->any_array.any->reference_count++;
						} else if (inner_operand->type == hol_term_type::AND) {
							dst = hol_term::new_exists(set_variable, hol_term::new_and(
									hol_term::new_equals(set_var, hol_term::new_lambda(element_variable,
										hol_term::new_and(make_array_view(inner_operand->array.operands, inner_operand->array.length - 1)))),
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
						}
					} else {
						return (hol_term*) nullptr;
					}
				}
				return dst;
			}, gather_scope);

	if (!result) return false;
	hol_term* new_dst = hol_term::new_lambda(lambda_variable, dst);
	if (new_dst == nullptr) {
		free(*dst); if (dst->reference_count == 0) free(dst);
		return false;
	}
	dst = new_dst;
	return true;
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
	case function_type::ADD_REQ_AUX:
		return add_flag(src, dst, grammatical_flag::REQ_AUX);
	case function_type::REMOVE_REQ_AUX:
		return remove_flag(src, dst, grammatical_flag::REQ_AUX);
	case function_type::REQUIRE_NO_REQ_AUX:
		return require_no_flag(src, dst, grammatical_flag::REQ_AUX);
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
	}
	fprintf(stderr, "apply ERROR: Unrecognized transformation function.\n");
	return false;
}

struct head_substituter {
	const hol_term* src;
	hol_term* dst;
};

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, const head_substituter& substituter)
{
	if (src == substituter.src) {
		substituter.dst->reference_count++;
		return substituter.dst;
	}

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	bool changed; hol_term* new_term;
	hol_term* first; hol_term* second; hol_term* third;
	switch (Type) {
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		first = src->binary.left;
		second = apply(src->binary.right, substituter);
		if (second == nullptr) {
			return nullptr;
		} else if (second == src->binary.right) {
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
		third = apply(src->ternary.third, substituter);
		if (third == nullptr) {
			return nullptr;
		} else if (third == src->ternary.third) {
			return src;
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
			return new_term;
		}

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
		changed = false;
		first = apply(src->any_array.right, substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->any_array.right) {
			changed = true;
		}

		if (!changed) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->any_array.right) { free(*first); if (first->reference_count == 0) free(first); }
				return nullptr;
			}
			new_term->type = Type;
			new_term->reference_count = 1;
			new_term->any_array.oper = src->any_array.oper;
			new_term->any_array.all = src->any_array.all;
			new_term->any_array.any = src->any_array.any;
			new_term->any_array.left = src->any_array.left;
			new_term->any_array.right = first;
			new_term->any_array.all->reference_count++;
			new_term->any_array.left->reference_count++;
			if (new_term->any_array.any != nullptr)
				new_term->any_array.any->reference_count++;
			if (first == src->any_array.right)
				first->reference_count++;
			return new_term;
		}

	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::NOT:
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY:
		return default_apply<Type>(src, substituter);
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type when substituting head.\n");
	return NULL;
}

inline hol_term* substitute_head(hol_term* src,
		const hol_term* src_term, hol_term* dst_term)
{
	const head_substituter substituter = {src_term, dst_term};
	hol_term* dst = apply(src, substituter);
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
		if (node->type == hol_term_type::IF_THEN) {
			has_bound_variables |= max_bound_variable(*node->binary.left, max_variable);
		} else if (node->type == hol_term_type::AND || node->type == hol_term_type::OR) {
			for (unsigned int i = 0; i + 1 < node->array.length; i++)
				has_bound_variables |= max_bound_variable(*node->array.operands[i], max_variable);
		}
		outer.add(node);
	}
};

template<typename TryFindHeadFunction, typename Function>
hol_term* find_head(hol_term* term, int& predicate_index, TryFindHeadFunction try_find_head, Function apply)
{
	apply(term);

	hol_term* head; bool negated;
	try_find_head(term, head, predicate_index, negated);
	if (head != nullptr) return head;

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	switch (term->type) {
	case hol_term_type::VARIABLE:
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IFF:
		return nullptr;

	case hol_term_type::NOT:
		return find_head(term->unary.operand, predicate_index, try_find_head, apply);

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return find_head(term->quantifier.operand, predicate_index, try_find_head, apply);

	case hol_term_type::IF_THEN:
		return find_head(term->binary.right, predicate_index, try_find_head, apply);

	case hol_term_type::AND:
	case hol_term_type::OR:
		return find_head(term->array.operands[term->array.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY:
		return find_head(term->any.included, predicate_index, try_find_head, apply);

	case hol_term_type::ANY_ARRAY:
		return find_head(term->any_array.right, predicate_index, try_find_head, apply);
	}
	fprintf(stderr, "find_head ERROR: Unrecognied hol_term_type.\n");
	return nullptr;
}

template<
	typename FindFirstHeadFunction,
	typename FindSecondHeadFunction>
inline bool remap_to_invert_apply_head(
		hol_term* first, hol_term* second,
		hol_term*& remapped_first,
		hol_term*& remapped_second,
		FindFirstHeadFunction find_first_head,
		FindSecondHeadFunction find_second_head)
{
	apply_head_inverter first_inverter; int first_predicate_index;
	apply_head_inverter second_inverter; int second_predicate_index;
	hol_term* first_head = find_head(first, first_predicate_index, find_first_head, first_inverter);
	hol_term* second_head = find_head(second, second_predicate_index, find_second_head, second_inverter);
	if (first_head == nullptr || second_head == nullptr)
		return false;

	unsigned int max_variable = second_inverter.max_variable;
	bool has_bound_variables = second_inverter.has_bound_variables;
	has_bound_variables |= max_bound_variable(*first_head, max_variable);
	has_bound_variables |= max_bound_variable(*second_head, max_variable);

	array_map<unsigned int, unsigned int> first_variable_map(8);
	for (unsigned int i = first_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = first_inverter.outer[i];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;
		if (has_bound_variables && node->quantifier.variable <= max_variable) {
			if (!first_variable_map.put(node->quantifier.variable, max_variable + 1))
				return false;
			max_variable++;
		} else {
			has_bound_variables = true;
		}
	}

	has_bound_variables |= first_inverter.has_bound_variables;
	max_variable = max(max_variable, first_inverter.max_variable);
	array_map<unsigned int, unsigned int> second_variable_map(8);
	for (unsigned int i = second_inverter.outer.length - 1; i > 0; i--) {
		const hol_term* node = second_inverter.outer[i];
		if (node->type != hol_term_type::FOR_ALL && node->type != hol_term_type::EXISTS && node->type != hol_term_type::LAMBDA)
			continue;
		if (has_bound_variables && node->quantifier.variable <= max_variable) {
			if (!second_variable_map.put(node->quantifier.variable, max_variable + 1))
				return false;
			max_variable++;
		} else {
			has_bound_variables = true;
		}
	}

	remapped_first = map_variables(remapped_first, first_variable_map);
	remapped_second = map_variables(remapped_second, second_variable_map);
	if (remapped_first == nullptr || remapped_second) {
		if (remapped_first != nullptr) {
			free(*remapped_first);
			if (remapped_first->reference_count == 0)
				free(remapped_first);
		}
		return false;
	}
	return true;
}

template<typename InvertSecondFunction,
	typename FindFirstHeadFunction,
	typename FindSecondHeadFunction>
inline bool invert_apply_head(
		flagged_logical_form<hol_term>*& inverse,
		unsigned int& inverse_count,
		const grammatical_flags& flags,
		hol_term* first, hol_term* second,
		FindFirstHeadFunction find_first_head,
		FindSecondHeadFunction find_second_head,
		InvertSecondFunction invert_second_head)
{
	hol_term* remapped_first;
	hol_term* remapped_second;
	if (!remap_to_invert_apply_head(first, second, remapped_first, remapped_second, find_first_head, find_second_head))
		return false;

	apply_head_inverter first_inverter; int first_predicate_index;
	apply_head_inverter second_inverter; int second_predicate_index;
	hol_term* first_head = find_head(remapped_first, first_predicate_index, find_first_head, first_inverter);
	hol_term* second_head = find_head(remapped_second, second_predicate_index, find_second_head, second_inverter);
	if (first_head == nullptr || second_head == nullptr) {
		free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		return false;
	}

	hol_term* conjunct;
	array<hol_term*> new_second_heads(8);
	if (!invert_second_head(new_second_heads, first_head, second_head, first_predicate_index, second_predicate_index, conjunct)) {
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
				if (!invert_second_head(new_second_head_array[i], first_head, parent->array.operands[i], first_predicate_index, second_predicate_index, conjunct)) {
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
		new_seconds[new_seconds.length] = substitute_head(remapped_second, second_head, new_second_head);
		free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
		free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);
		if (new_seconds[new_seconds.length] == nullptr) {
			free(*remapped_first); if (remapped_first->reference_count == 0) free(remapped_first);
			free(*remapped_second); if (remapped_second->reference_count == 0) free(remapped_second);
			return false;
		}
		new_seconds.length++;
	}
	for (hol_term* term : new_second_heads) { free(*term); if (term->reference_count == 0) free(term); }

	array<hol_term*> intersection(16);
	for (hol_term* new_second : new_seconds)
		intersect<built_in_predicates>(intersection, remapped_first, new_second);
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
		inverse[i].root = intersection[i];
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
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, find_head<built_in_predicates>,
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, int first_predicate_index, int second_predicate_index, hol_term*& conjunct) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;
				if (first_head->type == hol_term_type::ANY || (first_head->type == hol_term_type::EXISTS
				 && (first_head->quantifier.operand->type == hol_term_type::ANY || first_head->quantifier.operand->type == hol_term_type::ANY_ARRAY)))
				{
					hol_term* excluded_subtrees[4];
					excluded_subtrees[0] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
					excluded_subtrees[1] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
					excluded_subtrees[2] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
					excluded_subtrees[3] = hol_term::new_apply(
							hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()),
							hol_term::new_variable(second_head->quantifier.variable));
					if (excluded_subtrees[0] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[1] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[2] != nullptr) HOL_ANY.reference_count++;
					if (excluded_subtrees[0] == nullptr || excluded_subtrees[1] == nullptr || excluded_subtrees[2] == nullptr || excluded_subtrees[3] == nullptr) {
						if (excluded_subtrees[0] != nullptr) { free(*excluded_subtrees[0]); free(excluded_subtrees[0]); }
						if (excluded_subtrees[1] != nullptr) { free(*excluded_subtrees[1]); free(excluded_subtrees[1]); }
						if (excluded_subtrees[2] != nullptr) { free(*excluded_subtrees[2]); free(excluded_subtrees[2]); }
						return false;
					}

					hol_term* conjunct = hol_term::new_any(&HOL_ANY, excluded_subtrees, 4, nullptr, 0);
					if (conjunct == nullptr) {
						free(*excluded_subtrees[0]); free(excluded_subtrees[0]);
						free(*excluded_subtrees[1]); free(excluded_subtrees[1]);
						free(*excluded_subtrees[2]); free(excluded_subtrees[2]);
						free(*excluded_subtrees[3]); free(excluded_subtrees[3]);
						return false;
					}
					HOL_ANY.reference_count++;

					hol_term* conjunction;
					if (!new_hol_term(conjunction)) {
						free(*conjunct); free(conjunct);
						return false;
					}
					conjunction->type = hol_term_type::ANY_ARRAY;
					conjunction->reference_count = 1;
					conjunction->any_array.oper = hol_term_type::AND;
					conjunction->any_array.all = conjunct;
					conjunction->any_array.any = second_head_operand->array.operands[0];
					conjunction->any_array.any->reference_count++;
					if (ConjunctIndex == 0) {
						conjunction->any_array.left = second_head_operand->array.operands[1];
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.right = conjunct;
						conjunction->any_array.right->reference_count++;
					} else if (ConjunctIndex == -1) {
						conjunction->any_array.right = second_head_operand->array.operands[1];
						conjunction->any_array.right->reference_count++;
						conjunction->any_array.left = conjunct;
						conjunction->any_array.left->reference_count++;
					} else {
						fprintf(stderr, "invert_select_conjunct ERROR: Unsupported value for `ConjunctIndex`.\n");
						free(*conjunct); free(conjunct);
						free(conjunction);
						return false;
					}

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
					if (conjunct_index == first_predicate_index)
						fprintf(stderr, "invert_select_conjunct ERROR: `conjunct_index` and `first_predicate_index` are the same.\n");
#endif

					array<hol_term*> predicates(8);
					array<hol_term*> conjuncts(8);
					intersect<built_in_predicates>(predicates, first_head_operand->array.operands[first_predicate_index], second_head_operand->array.operands[0]);
					intersect<built_in_predicates>(conjuncts, first_head_operand->array.operands[conjunct_index], second_head_operand->array.operands[1]);
					unsigned int intersection_count = predicates.length * conjuncts.length;
					if (!dst.ensure_capacity(dst.length + intersection_count)) {
						for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
						for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
						return false;
					}
					for (hol_term* predicate : predicates) {
						for (hol_term* conjunct : conjuncts) {
							if (!new_hol_term(dst[dst.length])) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								return false;
							}
							dst[dst.length]->type = hol_term_type::AND;
							dst[dst.length]->reference_count = 1;
							dst[dst.length]->array.length = first_head_operand->array.length;
							dst[dst.length]->array.operands = (hol_term**) malloc(sizeof(hol_term*) * first_head_operand->array.length);
							if (dst[dst.length]->array.operands == nullptr) {
								for (hol_term* term : predicates) { free(*term); if (term->reference_count == 0) free(term); }
								for (hol_term* term : conjuncts) { free(*term); if (term->reference_count == 0) free(term); }
								free(dst[dst.length]); return false;
							}
							for (unsigned int i = 0; i < first_head_operand->array.length; i++) {
								if ((int) i == first_predicate_index) {
									dst[dst.length]->array.operands[i] = predicate;
								} else if (i == (unsigned int) conjunct_index) {
									dst[dst.length]->array.operands[i] = conjunct;
								} else {
									dst[dst.length]->array.operands[i] = first_head_operand->array.operands[i];
								}
								dst[dst.length]->array.operands[i]->reference_count++;
							}
							dst.length++;
						}
					}
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
		find_head<built_in_predicates>, find_head<built_in_predicates>,
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, int first_predicate_index, int second_predicate_index, hol_term*& conjunct) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;
				hol_term* excluded_subtrees[4];
				excluded_subtrees[0] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG1), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
				excluded_subtrees[1] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG2), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
				excluded_subtrees[2] = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::ARG3), &HOL_ANY), hol_term::new_variable(second_head->quantifier.variable));
				excluded_subtrees[3] = hol_term::new_apply(
						hol_term::new_any(nullptr, nullptr, 0, hol_non_head_constants<built_in_predicates>::get_terms(), hol_non_head_constants<built_in_predicates>::count()),
						hol_term::new_variable(second_head->quantifier.variable));
				if (excluded_subtrees[0] != nullptr) HOL_ANY.reference_count++;
				if (excluded_subtrees[1] != nullptr) HOL_ANY.reference_count++;
				if (excluded_subtrees[2] != nullptr) HOL_ANY.reference_count++;
				if (excluded_subtrees[0] == nullptr || excluded_subtrees[1] == nullptr || excluded_subtrees[2] == nullptr || excluded_subtrees[3] == nullptr) {
					if (excluded_subtrees[0] != nullptr) { free(*excluded_subtrees[0]); free(excluded_subtrees[0]); }
					if (excluded_subtrees[1] != nullptr) { free(*excluded_subtrees[1]); free(excluded_subtrees[1]); }
					if (excluded_subtrees[2] != nullptr) { free(*excluded_subtrees[2]); free(excluded_subtrees[2]); }
					return false;
				}

				hol_term* conjunct = hol_term::new_any(&HOL_ANY, excluded_subtrees, 4, nullptr, 0);
				if (conjunct == nullptr) {
					free(*excluded_subtrees[0]); free(excluded_subtrees[0]);
					free(*excluded_subtrees[1]); free(excluded_subtrees[1]);
					free(*excluded_subtrees[2]); free(excluded_subtrees[2]);
					free(*excluded_subtrees[3]); free(excluded_subtrees[3]);
					return false;
				}
				HOL_ANY.reference_count++;

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
						conjunction->array.operands[i] = second_head_operand->array.operands[i];
						conjunction->array.operands[i]->reference_count++;
					}
				}

				hol_term* new_second_head = hol_term::new_exists(second_head->quantifier.variable, conjunction);
				if (new_second_head == nullptr) {
					free(*conjunction); free(conjunction);
					return false;
				}
				intersect<built_in_predicates>(dst, new_second_head, first_head);
				free(*new_second_head); if (new_second_head->reference_count == 0) free(new_second_head);
				return (dst.length > 0);
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
		find_head<built_in_predicates>, find_head<built_in_predicates>,
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, int first_predicate_index, int second_predicate_index, hol_term*& conjunct) {
			if (second_head->type == hol_term_type::EXISTS && second_head->quantifier.operand->type == hol_term_type::AND) {
				hol_term* second_head_operand = second_head->quantifier.operand;

				if (second_predicate_index < 0)
					second_predicate_index += second_head_operand->array.length;
				hol_term* second_predicate = second_head_operand->array.operands[second_predicate_index];
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
						if (i == (unsigned int) second_predicate_index) {
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
				if (second_predicate_index == 0)
					second_predicate = second_head_operand->any_array.left;
				else if (second_predicate_index == -1)
					second_predicate = second_head_operand->any_array.right;
				else if (second_predicate_index == INT_MAX)
					second_predicate = second_head_operand->any_array.any;
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
					if (second_predicate_index == 0) {
						conjunction->any_array.left = new_predicate;
						conjunction->any_array.right = second_head_operand->any_array.right;
						conjunction->any_array.right->reference_count++;
						conjunction->any_array.any = second_head_operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					} else if (second_predicate_index == -1) {
						conjunction->any_array.right = new_predicate;
						conjunction->any_array.left = second_head_operand->any_array.left;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.any = second_head_operand->any_array.any;
						if (conjunction->any_array.any != nullptr)
							conjunction->any_array.any->reference_count++;
					} else if (second_predicate_index == INT_MAX) {
						conjunction->any_array.any = new_predicate;
						conjunction->any_array.left = second_head_operand->any_array.left;
						conjunction->any_array.left->reference_count++;
						conjunction->any_array.right = second_head_operand->any_array.right;
						conjunction->any_array.right->reference_count++;
					}
					new_predicate->reference_count++;

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
				fprintf(stderr, "invert_remove_conjunct ERROR: Expected `second_head` to be an existentially quantified conjunction.\n");
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
	return invert_apply_head(inverse, inverse_count, flags, first, second,
		find_head<built_in_predicates>, predicative_head_finder<built_in_predicates>(lambda_variable),
		[](array<hol_term*>& dst, hol_term* first_head, hol_term* second_head, int first_predicate_index, int second_predicate_index, hol_term*& conjunct) {
			/* TODO: implement this */
			return false;
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
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::ANY))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.flags[(uint_fast8_t) flag], first.flags.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	}
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.cnj, first.flags.cnj, cnj))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::ANY))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
	return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	flagged_logical_form<Formula> exp;
	switch (function.type) {
	case function_type::EMPTY:
		inverse_count = 1; if (!init_array(inverse, inverse_count)) return false;
		*inverse = first;
		return true;
	case function_type::IDENTITY:
	case function_type::REQUIRE_NO_INVERSE:
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE:
	case function_type::REQUIRE_LEFT_PREDICATE_INVERSE_OWN:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_RIGHT_CONJUNCT:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_conjunct<-1>(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REMOVE_INVERSE:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_remove_inverse(inverse, inverse_count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG1_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG1>(inverse, inverse_Count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG2_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG2>(inverse, inverse_Count, flags, first.root, second.root);
	case function_type::SELECT_RIGHT_ARG3_WITHOUT_HEAD_PREDICATIVE_AND_FACTOR:
		if (!intersect(flags, first.flags, second.flags)) return false;
		return invert_select_arg_without_head_predicative_and_factor<-1, (unsigned int) built_in_predicates::ARG3>(inverse, inverse_Count, flags, first.root, second.root);
	case function_type::ADD_SINGULAR:
		if ((second.flags.index_number != grammatical_num::SINGULAR && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::ADD_PLURAL:
		if ((second.flags.index_number != grammatical_num::PLURAL && second.flags.index_number != grammatical_num::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, grammatical_num::NONE)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
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
		 || !intersect(flags.cnj, first.flags.cnj, grammatical_conjunction::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
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
	case function_type::ADD_REQ_AUX:
		return invert_add_flag(inverse, inverse_count, first, second, grammatical_flag::REQ_AUX);
	case function_type::REMOVE_REQ_AUX:
		return invert_remove_flag(inverse, inverse_count, first, second, grammatical_flag::REQ_AUX);
	case function_type::REQUIRE_NO_REQ_AUX:
		return invert_require_no_flag(inverse, inverse_count, first, second, grammatical_flag::REQ_AUX);
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
		 || !intersect(flags.corr, first.flags.corr, correlator::ANY))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NO_CORRELATOR:
		if ((second.flags.corr != correlator::NONE && second.flags.corr != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, second.flags.correlated_by)
		 || !intersect(flags.corr, first.flags.corr, correlator::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
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
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::ANY))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	case function_type::REQUIRE_NOT_CORRELATED:
		if ((second.flags.correlated_by != correlator::NONE && second.flags.correlated_by != correlator::ANY)
		 || !intersect(flags.index_number, first.flags.index_number, second.flags.index_number)
		 || !intersect(flags.concord_number, first.flags.concord_number, second.flags.concord_number)
		 || !intersect(flags.cnj, first.flags.cnj, second.flags.cnj)
		 || !intersect(flags.corr, first.flags.corr, second.flags.corr)
		 || !intersect(flags.correlated_by, first.flags.correlated_by, correlator::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(flags.flags[i], first.flags.flags[i], second.flags.flags[i])) return false;
		return intersect(inverse, inverse_count, flags, first.root, second.root);
	}
	fprintf(stderr, "invert ERROR: Unrecognized transformation function.\n");
	return false;
}

template<typename Formula>
static void is_separable(
		const transformation<flagged_logical_form<Formula>>* functions,
		unsigned int rule_length, bool* separable)
{
	for (unsigned int i = 0; i < rule_length; i++)
		separable[i] = false;
}

template<typename Formula>
bool get_feature(
		typename flagged_logical_form<Formula>::feature feature,
		const flagged_logical_form<Formula>& src, unsigned int& value,
		unsigned int*& excluded, unsigned int& excluded_count)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	switch (feature) {
	case feature_type::CONSTANT: /* TODO: implement this */
	case feature_type::EMPTY: break;
	}
	fprintf(stderr, "get_feature ERROR: Unrecognized semantic feature.\n");
	return false;
}

template<typename Formula>
bool set_feature(
		typename flagged_logical_form<Formula>::feature feature,
		flagged_logical_form<Formula>& exp, unsigned int value)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	switch (feature) {
	case feature_type::CONSTANT: /* TODO: implement this */
	case feature_type::EMPTY: break;
	}
	fprintf(stderr, "set_feature ERROR: Unrecognized semantic feature.\n");
	exit(EXIT_FAILURE);
}

template<typename Formula>
bool exclude_features(typename flagged_logical_form<Formula>::feature feature,
		flagged_logical_form<Formula>& exp, const unsigned int* values, unsigned int count)
{
	typedef typename flagged_logical_form<Formula>::feature feature_type;
	switch (feature) {
	case feature_type::CONSTANT: /* TODO: implement this */
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
constexpr bool morphology_is_valid(const morphology& morph,
		const sequence& terminal, PartOfSpeechType pos,
		const flagged_logical_form<Formula>& logical_form)
{
	return true;
}

template<typename Formula, typename PartOfSpeechType, typename EmitRootFunction>
bool morphology_parse(
		const morphology& morph, const sequence& words, PartOfSpeechType pos,
		const flagged_logical_form<Formula>& logical_form, EmitRootFunction emit_root)
{
	return emit_root(words, logical_form);
}

#endif /* HDP_PARSER_H_ */
