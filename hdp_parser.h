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

template<typename K, typename V>
struct static_pair {
	K key;
	V value;
};

template<typename Formula>
struct flagged_logical_form
{
	Formula* root;
	grammatical_num index_number;
	grammatical_num concord_number;
	grammatical_conjunction cnj;
	grammatical_flag_value flags[(uint_fast8_t) grammatical_flag::COUNT] = {0};
	correlator corr;
	correlator correlated_by;

	flagged_logical_form() :
			index_number(grammatical_num::NONE), concord_number(grammatical_num::NONE),
			cnj(grammatical_conjunction::NONE), corr(correlator::NONE), correlated_by(correlator::NONE) { }

	flagged_logical_form(Formula& src) : root(&src),
			index_number(grammatical_num::NONE), concord_number(grammatical_num::NONE),
			cnj(grammatical_conjunction::NONE), corr(correlator::NONE), correlated_by(correlator::NONE)
	{
		src.reference_count++;
	}

	~flagged_logical_form() { free_helper(); }

	inline void copy_flags_from(const flagged_logical_form<Formula>& src) {
		index_number = src.index_number;
		concord_number = src.concord_number;
		cnj = src.cnj;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			flags[i] = src.flags[i];
		corr = src.corr;
		correlated_by = src.correlated_by;
	}

	inline void operator = (const flagged_logical_form<Formula>& src) {
		copy_flags_from(src);
		root = src.root;
		root->reference_count++;
	}

	static inline bool is_empty(const flagged_logical_form<Formula>& src) {
		return src.root == nullptr;
	}

	static inline unsigned int hash(const flagged_logical_form<Formula>& src) {
		return hasher<Formula>::hash(*src.root)
			 ^ (3 * default_hash(src.index_number))
			 ^ (7 * default_hash(src.concord_number))
			 ^ (31 * default_hash(src.cnj))
			 ^ (127 * default_hash(src.corr))
			 ^ (8191 * default_hash(src.correlated_by))
			 ^ (131071 * default_hash(src.flags, (uint_fast8_t) grammatical_flag::COUNT));
	}

	static inline void move(
			const flagged_logical_form<Formula>& src,
			flagged_logical_form<Formula>& dst)
	{
		dst.root = src.root;
		dst.copy_flags_from(src);
	}

	static inline void swap(
			flagged_logical_form<Formula>& first,
			flagged_logical_form<Formula>& second)
	{
		core::swap(first.root, second.root);
		core::swap(first.index_number, second.index_number);
		core::swap(first.concord_number, second.concord_number);
		core::swap(first.cnj, second.cnj);
		core::swap(first.corr, second.corr);
		core::swap(first.correlated_by, second.correlated_by);
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			core::swap(first.flags[i], second.flags[i]);
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
		SELECT_ARG1_KEEP_HEAD,
		REMOVE_ARG1,

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
	{function_type::SELECT_ARG1_KEEP_HEAD, "select_arg1_keep_head"},
	{function_type::REMOVE_ARG1, "remove_arg1"},
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

template<typename Formula>
inline bool initialize_any(flagged_logical_form<Formula>& lf) {
	lf.root = (Formula*) malloc(sizeof(Formula));
	if (lf.root == nullptr) {
		fprintf(stderr, "initialize_any ERROR: Out of memory.\n");
		return false;
	}
	initialize_any(*lf.root);
	lf.concord_number = grammatical_num::ANY;
	lf.index_number = grammatical_num::ANY;
	lf.cnj = grammatical_conjunction::ANY;
	lf.corr = correlator::ANY;
	lf.correlated_by = correlator::ANY;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		lf.flags[i] = grammatical_flag_value::ANY;
	return true;
}

template<typename Formula>
inline bool operator == (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	if (first.index_number != second.index_number
	 || first.concord_number != second.concord_number
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by)
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (first.flags[i] != second.flags[i]) return false;
	return (first.root == second.root) || (*first.root == *second.root);
}

template<typename Formula>
inline bool operator != (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	if (first.index_number != second.index_number
	 || first.concord_number != second.concord_number
	 || first.cnj != second.cnj
	 || first.corr != second.corr
	 || first.correlated_by != second.correlated_by)
		return true;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (first.flags[i] != second.flags[i]) return true;
	return (first.root != second.root) && (*first.root != *second.root);
}

template<typename Formula, typename Stream, typename... Printer>
inline bool print(
		const flagged_logical_form<Formula>& lf,
		Stream& out, Printer&&... printer)
{
	bool first = true;
	if (!print(*lf.root, out, std::forward<Printer>(printer)...) || !print('[', out)
	 || (first && !print(',', out)) || !print_feature("index", first, lf.index_number, out)
	 || (first && !print(',', out)) || !print_feature("concord", first, lf.concord_number, out)
	 || (first && !print(',', out)) || !print_feature(first, lf.cnj, out)
	 || (first && !print(',', out)) || !print_feature("corr:", first, lf.corr, out)
	 || (first && !print(',', out)) || !print_feature("corr_by:", first, lf.correlated_by, out))
		return false;

	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if ((first && !print(',', out)) || !print_feature(first, (grammatical_flag) i, lf.flags[i], out)) return false;
	return print(']', out);
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

inline bool is_built_in(unsigned int constant) {
	return constant != (unsigned int) built_in_predicates::UNKNOWN
		&& constant != (unsigned int) built_in_predicates::ARG1
		&& constant != (unsigned int) built_in_predicates::ARG2
		&& constant != (unsigned int) built_in_predicates::ARG3
		&& constant != (unsigned int) built_in_predicates::SIZE;
}

inline unsigned int get_predicate_of_literal(
		const hol_term* src,
		unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.right->type == hol_term_type::VARIABLE
	 && src->binary.right->variable == scope_variable
	 && src->binary.left->type == hol_term_type::CONSTANT
	 && is_built_in(src->binary.left->constant))
	{
		return src->binary.left->constant;
	}
	return 0;
}

inline unsigned int get_predicate_of_scope(
		const hol_term* src,
		unsigned int scope_variable,
		unsigned int& index)
{
	if (src->type == hol_term_type::AND) {
		for (unsigned int i = 0; i < src->array.length; i++) {
			unsigned int predicate = get_predicate_of_literal(src->array.operands[i], scope_variable);
			if (predicate != 0) {
				index = i;
				return predicate;
			}
		}
	} else {
		unsigned int predicate = get_predicate_of_literal(src, scope_variable);
		if (predicate != 0) {
			index = 0;
			return predicate;
		}
	}
	return 0;
}

template<unsigned int Arg>
struct arg_finder {
	unsigned int head_variable;
	hol_term* arg;
};

template<hol_term_type Type, unsigned int Arg>
inline bool visit(const hol_term& term, arg_finder<Arg>& visitor) {
	if (Type == hol_term_type::EQUALS) {
		if (term.binary.left->type == hol_term_type::UNARY_APPLICATION
		 && term.binary.left->binary.left->type == hol_term_type::CONSTANT
		 && term.binary.left->binary.left->constant == Arg
		 && term.binary.left->binary.right->type == hol_term_type::VARIABLE
		 && term.binary.left->binary.right->variable == visitor.head_variable)
		{
			visitor.arg = term.binary.right;
			return false;
		}
	}
	return true;
}

template<unsigned int Arg>
inline hol_term* get_arg(const hol_term& src, unsigned int head_variable) {
	arg_finder<Arg> visitor;
	visitor.head_variable = head_variable;
	visitor.arg = nullptr;
	visit(src, visitor);
	return visitor.arg;
}

template<unsigned int Arg>
inline void find_arg(
		const hol_term* src,
		hol_term*& head,
		unsigned int& head_variable,
		unsigned int& predicate_index,
		unsigned int& arg_index,
		hol_term*& arg,
		bool& negated)
{
	arg = nullptr;
	if (src->type == hol_term_type::EXISTS) {
		unsigned int predicate = get_predicate_of_scope(src->quantifier.operand, src->quantifier.variable, predicate_index);
		if (predicate != (unsigned int) built_in_predicates::CAPABLE_OF)
		{
			/* we found the head */
			head = src->quantifier.operand;
			head_variable = src->quantifier.variable;

			if (head->type != hol_term_type::AND)
				return;

			/* now find the term containing `Arg(x)` where `x` is the head variable */
			for (unsigned int i = 0; i < head->array.length; i++) {
				arg = get_arg<Arg>(*head->array.operands[i], head_variable);
				if (arg != nullptr) {
					arg_index = i;
					break;
				}
			}
			negated = false;
		}
	} else if (src->type == hol_term_type::NOT) {
		find_arg<Arg>(src->unary.operand, head, head_variable, predicate_index, arg_index, arg, negated);
		negated = true;
	}
}

bool prune_independent_siblings(
		array<unsigned int>& arg_variables,
		array<hol_term*>& siblings)
{
	array<unsigned int>* sibling_variables = (array<unsigned int>*) malloc(max((size_t) 1, sizeof(array<unsigned int>) * siblings.length));
	for (unsigned int i = 0; i < siblings.length; i++) {
		if (!array_init(sibling_variables[i], 8)) {
			for (unsigned int j = 0; j < i; j++) free(sibling_variables[j]);
			free(sibling_variables);
			return false;
		}
		get_variables(*siblings[i], sibling_variables[i]);
	}

	unsigned int old_arg_variable_count = arg_variables.length;
	do {
		for (unsigned int i = 0; i < siblings.length; i++) {
			if (has_intersection(sibling_variables[i], arg_variables)) {
				array<unsigned int> union_variables(sibling_variables[i].length + arg_variables.length);
				set_union(union_variables, arg_variables, sibling_variables[i]);
				swap(union_variables, arg_variables);
			}
		}
	} while (old_arg_variable_count != arg_variables.length);

	/* remove terms from `siblings` that are unneeded */
	for (unsigned int i = 0; i < siblings.length; i++) {
		if (!has_intersection(sibling_variables[i], arg_variables)) {
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

template<unsigned int Arg, typename MakeDstFunction>
inline bool apply_arg(
		const hol_term* src, hol_term*& dst,
		array<unsigned int>& dst_variables,
		array<hol_term*>& siblings,
		MakeDstFunction make_dst)
{
	/* check if the current scope is the head */
	hol_term* head; hol_term* arg; bool negated;
	unsigned int head_variable, predicate_index, arg_index;
	find_arg<Arg>(src, head, head_variable, predicate_index, arg_index, arg, negated);
	if (arg != nullptr) {
		dst = make_dst(head, head_variable, predicate_index, arg_index, negated);
		if (dst == nullptr) return false;

		get_variables(*dst, dst_variables);
		dst_variables.remove(dst_variables.index_of(head_variable));
		if (dst_variables.length > 1) insertion_sort(dst_variables);
		return true;
	}

	if (src->type == hol_term_type::AND || src->type == hol_term_type::OR) {
		if (!siblings.ensure_capacity(siblings.length + src->array.length)) return false;
		find_arg<Arg>(src->array.operands[src->array.length - 1], head, head_variable, predicate_index, arg_index, arg, negated);

		if (arg != nullptr) {
			dst = make_dst(head, head_variable, predicate_index, arg_index, negated);
			if (dst == nullptr) return false;

			get_variables(*dst, dst_variables);
			dst_variables.remove(dst_variables.index_of(head_variable));
			if (dst_variables.length > 1) insertion_sort(dst_variables);

			/* check if the other conjuncts are also part of the head scope (i.e. the head scope is a conjunction) */
			hol_term* new_arg;
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
				find_arg<Arg>(src->array.operands[i], head, head_variable, predicate_index, arg_index, new_arg, negated);
				if (new_arg == nullptr || *new_arg != *arg) {
					for (unsigned int j = 0; j < i; j++) {
						free(*new_dst[j]);
						if (new_dst[j]->reference_count == 0)
							free(new_dst[j]);
					}
					new_arg = nullptr;
					break;
				}

				new_dst[i] = make_dst(head, head_variable, predicate_index, arg_index, negated);
				if (new_dst[i] == nullptr) {
					for (unsigned int j = 0; j < i; j++) {
						free(*new_dst[j]);
						if (new_dst[j]->reference_count == 0)
							free(new_dst[j]);
					}
					new_arg = nullptr;
					break;
				}

				array<unsigned int> temp_variables(8);
				get_variables(*dst, temp_variables);
				temp_variables.remove(temp_variables.index_of(head_variable));
				if (temp_variables.length > 1) insertion_sort(temp_variables);

				array<unsigned int> union_variables(new_variables.length + temp_variables.length);
				set_union(union_variables, new_variables, temp_variables);
				swap(union_variables, new_variables);
			}

			if (new_arg == nullptr) {
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

				unsigned int new_length = 0;
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					if (siblings.contains(src->array.operands[i])) {
						new_dst[new_length++] = src->array.operands[i];
						src->array.operands[i]->reference_count++;
					}
				}
				new_dst[new_length++] = new_dst[src->array.length - 1];
				if (new_length == 1) return true;
				if (src->type == hol_term_type::AND)
					dst = hol_term::new_and(make_array_view(new_dst, new_length));
				else dst = hol_term::new_or(make_array_view(new_dst, new_length));
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

	case hol_term_type::ANY:
		break;

	case hol_term_type::NOT:
		if (!apply_arg<Arg>(src->unary.operand, new_formula, dst_variables, siblings, make_dst))
			return false;
		if (new_formula != nullptr) {
			dst = hol_term::new_not(new_formula);
			if (dst == nullptr) {
				free(*new_formula);
				if (new_formula->reference_count == 0)
					free(new_formula);
				return false;
			}
		}
		return true;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!apply_arg<Arg>(src->unary.operand, new_formula, dst_variables, siblings, make_dst))
			return false;
		if (new_formula != nullptr) {
			if (!dst_variables.contains(src->quantifier.variable)) {
				dst = new_formula;
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
		return true;

	case hol_term_type::IF_THEN:
		if (!apply_arg<Arg>(src->unary.operand, new_formula, dst_variables, siblings, make_dst))
			return false;
		if (new_formula != nullptr) {
			dst = hol_term::new_if_then(src->binary.left, new_formula);
			if (dst == nullptr) {
				free(*new_formula);
				if (new_formula->reference_count == 0)
					free(new_formula);
				return false;
			}
		}
		return true;

	case hol_term_type::AND:
	case hol_term_type::OR:
		/* we already handle this case above */
		break;
	}
	fprintf(stderr, "apply_arg ERROR: Unrecognized hol_term_type.\n");
	return false;
}

template<unsigned int Arg>
inline bool select_arg_keep_head(
		const hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	return apply_arg<Arg>(src, dst, dst_variables, siblings,
			[](hol_term* head, unsigned int head_variable, unsigned int predicate_index, unsigned int arg_index, bool negated) {
				hol_term* predicate = head->array.operands[predicate_index];
				hol_term* arg = head->array.operands[arg_index];
				hol_term* dst = hol_term::new_exists(head_variable, hol_term::new_and(predicate, arg));
				if (dst != nullptr) {
					predicate->reference_count++;
					arg->reference_count++;
				}
				return dst;
			});
}

template<unsigned int Arg>
inline bool remove_arg(
		const hol_term* src, hol_term*& dst)
{
	array<hol_term*> siblings(8);
	array<unsigned int> dst_variables(8);
	return apply_arg<Arg>(src, dst, dst_variables, siblings,
			[](hol_term* head, unsigned int head_variable, unsigned int predicate_index, unsigned int arg_index, bool negated) {
				hol_term* dst;
				if (negated)
					dst = hol_term::new_not(hol_term::new_exists(head_variable, hol_term::new_and(make_excluded_array_view(head->array.operands, head->array.length, arg_index))));
				else dst = hol_term::new_exists(head_variable, hol_term::new_and(make_excluded_array_view(head->array.operands, head->array.length, arg_index)));
				if (dst != nullptr) {
					for (unsigned int i = 0; i < head->array.length; i++)
						if (i != arg_index) head->array.operands[i]->reference_count++;
				}
				return dst;
			});
}

template<typename Formula>
inline bool add_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE
	 && src.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags[(uint_fast8_t) flag] = grammatical_flag_value::TRUE;
	return true;
}

template<typename Formula>
inline bool remove_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE
	 && src.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool try_remove_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	dst = src;
	dst.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool require_no_flag(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_flag flag)
{
	if (src.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE
	 && src.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
		return false;
	dst = src;
	dst.flags[(uint_fast8_t) flag] = grammatical_flag_value::FALSE;
	return true;
}

template<typename Formula>
inline bool add_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.cnj != grammatical_conjunction::NONE && src.cnj != grammatical_conjunction::ANY)
		return false;
	dst = src;
	dst.cnj = cnj;
	return true;
}

template<typename Formula>
inline bool remove_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.cnj != cnj && src.cnj != grammatical_conjunction::ANY)
		return false;
	dst = src;
	dst.cnj = grammatical_conjunction::NONE;
	return true;
}

template<typename Formula>
inline bool require_no_conjunction(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		grammatical_conjunction cnj)
{
	if (src.cnj == cnj)
		return false;
	dst = src;
	dst.cnj = src.cnj;
	return true;
}

template<typename Formula>
inline bool add_correlator(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator corr)
{
	if (src.corr != correlator::NONE && src.corr != correlator::ANY)
		return false;
	dst = src;
	dst.corr = corr;
	return true;
}

template<typename Formula>
inline bool remove_correlator(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator corr)
{
	if (src.corr != corr && src.corr != correlator::ANY)
		return false;
	dst = src;
	dst.corr = correlator::NONE;
	return true;
}

template<typename Formula>
inline bool add_correlated_by(
		const flagged_logical_form<Formula>& src,
		flagged_logical_form<Formula>& dst,
		correlator correlated_by)
{
	if (src.correlated_by != correlator::NONE && src.correlated_by != correlator::ANY)
		return false;
	dst = src;
	dst.correlated_by = correlated_by;
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
		dst.index_number = grammatical_num::NONE;
		dst.concord_number = grammatical_num::NONE;
		dst.cnj = grammatical_conjunction::NONE;
		dst.corr = correlator::NONE;
		dst.correlated_by = correlator::NONE;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			dst.flags[i] = grammatical_flag_value::FALSE;
		return dst.root != nullptr;
	case function_type::IDENTITY:
		dst = src;
		return true;
	case function_type::SELECT_ARG1_KEEP_HEAD:
		dst.copy_flags_from(src);
		return select_arg_keep_head<(unsigned int) built_in_predicates::ARG1>(src.root, dst.root);
	case function_type::REMOVE_ARG1:
		dst.copy_flags_from(src);
		return remove_arg<(unsigned int) built_in_predicates::ARG1>(src.root, dst.root);
	case function_type::ADD_SINGULAR:
		if (src.index_number != grammatical_num::NONE && src.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.index_number = grammatical_num::SINGULAR;
		return true;
	case function_type::ADD_PLURAL:
		if (src.index_number != grammatical_num::NONE && src.index_number != grammatical_num::ANY)
			return false;
		dst = src;
		dst.index_number = grammatical_num::PLURAL;
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
		if (src.cnj != grammatical_conjunction::NONE && src.cnj != grammatical_conjunction::ANY)
			return false;
		dst = src;
		dst.cnj = grammatical_conjunction::NONE;
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
		dst.corr = correlator::NONE;
		return true;
	case function_type::REQUIRE_NO_CORRELATOR:
		if (src.corr != correlator::NONE && src.corr != correlator::ANY)
			return false;
		dst = src;
		dst.corr = correlator::NONE;
		return true;
	case function_type::ADD_CORRELATED_BY_BOTH:
		return add_correlated_by(src, dst, correlator::BOTH);
	case function_type::ADD_CORRELATED_BY_EITHER:
		return add_correlated_by(src, dst, correlator::EITHER);
	case function_type::ADD_CORRELATED_BY_NEITHER:
		return add_correlated_by(src, dst, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATED:
		dst = src;
		dst.correlated_by = correlator::NONE;
		return true;
	case function_type::REQUIRE_NOT_CORRELATED:
		if (src.correlated_by != correlator::NONE && src.correlated_by != correlator::ANY)
			return false;
		dst = src;
		dst.correlated_by = correlator::NONE;
		return true;
	}
	fprintf(stderr, "apply ERROR: Unrecognized transformation function.\n");
	return false;
}

template<typename Formula>
inline bool invert_add_flag(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	if ((second.flags[(uint_fast8_t) flag] != grammatical_flag_value::TRUE && second.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.flags[(uint_fast8_t) flag], first.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	}
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_flag(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	if ((second.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.flags[(uint_fast8_t) flag], first.flags[(uint_fast8_t) flag], grammatical_flag_value::TRUE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	}
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_try_remove_flag(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	if ((second.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.flags[(uint_fast8_t) flag], first.flags[(uint_fast8_t) flag], grammatical_flag_value::ANY))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	}
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_require_no_flag(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_flag flag)
{
	if ((second.flags[(uint_fast8_t) flag] != grammatical_flag_value::FALSE && second.flags[(uint_fast8_t) flag] != grammatical_flag_value::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.flags[(uint_fast8_t) flag], first.flags[(uint_fast8_t) flag], grammatical_flag_value::FALSE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++) {
		if ((grammatical_flag) i == flag) continue;
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	}
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_conjunction(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	if ((second.cnj != cnj && second.cnj != grammatical_conjunction::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.cnj, first.cnj, grammatical_conjunction::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_conjunction(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	if ((second.cnj != grammatical_conjunction::NONE && second.cnj != grammatical_conjunction::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.cnj, first.cnj, cnj))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_require_no_conjunction(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		grammatical_conjunction cnj)
{
	if (second.cnj == cnj
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.cnj, first.cnj, grammatical_conjunction::ANY))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_correlator(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator corr)
{
	if ((second.corr != corr && second.corr != correlator::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.corr, first.corr, correlator::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_remove_correlator(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator corr)
{
	if ((second.corr != correlator::NONE && second.corr != correlator::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.corr, first.corr, corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
	 || !intersect(inverse.cnj, first.cnj, second.cnj))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert_add_correlated_by(
		flagged_logical_form<Formula>& inverse,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second,
		correlator correlated_by)
{
	if ((second.correlated_by != correlated_by && second.correlated_by != correlator::ANY)
	 || !intersect(inverse.index_number, first.index_number, second.index_number)
	 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
	 || !intersect(inverse.cnj, first.cnj, second.cnj)
	 || !intersect(inverse.corr, first.corr, second.corr)
	 || !intersect(inverse.correlated_by, first.correlated_by, correlator::NONE))
		return false;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
	return intersect(inverse.root, first.root, second.root);
}

template<typename Formula>
inline bool invert(
		flagged_logical_form<Formula>& inverse,
		typename flagged_logical_form<Formula>::function function,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	typedef typename flagged_logical_form<Formula>::function_type function_type;

	uint_fast8_t i;
	flagged_logical_form<Formula> exp;
	switch (function.type) {
	case function_type::EMPTY:
		inverse = first;
		return true;
	case function_type::IDENTITY:
		if (!intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::ADD_SINGULAR:
		if ((second.index_number != grammatical_num::SINGULAR && second.index_number != grammatical_num::ANY)
		 || !intersect(inverse.index_number, first.index_number, grammatical_num::NONE)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::ADD_PLURAL:
		if ((second.index_number != grammatical_num::PLURAL && second.index_number != grammatical_num::ANY)
		 || !intersect(inverse.index_number, first.index_number, grammatical_num::NONE)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::ADD_THAT:
		return invert_add_conjunction(inverse, first, second, grammatical_conjunction::THAT);
	case function_type::REMOVE_THAT:
		return invert_remove_conjunction(inverse, first, second, grammatical_conjunction::THAT);
	case function_type::REQUIRE_NO_THAT:
		return invert_require_no_conjunction(inverse, first, second, grammatical_conjunction::THAT);
	case function_type::ADD_WHETHER:
		return invert_add_conjunction(inverse, first, second, grammatical_conjunction::WHETHER);
	case function_type::REMOVE_WHETHER:
		return invert_remove_conjunction(inverse, first, second, grammatical_conjunction::WHETHER);
	case function_type::ADD_IF:
		return invert_add_conjunction(inverse, first, second, grammatical_conjunction::IF);
	case function_type::REMOVE_IF:
		return invert_remove_conjunction(inverse, first, second, grammatical_conjunction::IF);
	case function_type::ADD_BECAUSE:
		return invert_add_conjunction(inverse, first, second, grammatical_conjunction::BECAUSE);
	case function_type::REMOVE_BECAUSE:
		return invert_remove_conjunction(inverse, first, second, grammatical_conjunction::BECAUSE);
	case function_type::ADD_FOR:
		return invert_add_conjunction(inverse, first, second, grammatical_conjunction::FOR);
	case function_type::REMOVE_FOR:
		return invert_remove_conjunction(inverse, first, second, grammatical_conjunction::FOR);
	case function_type::REQUIRE_NO_CONJUNCTION:
		if ((second.cnj != grammatical_conjunction::NONE && second.cnj != grammatical_conjunction::ANY)
		 || !intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
		 || !intersect(inverse.cnj, first.cnj, grammatical_conjunction::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::ADD_IS_ADJUNCT:
		return invert_add_flag(inverse, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::TRY_REMOVE_IS_ADJUNCT:
		return invert_try_remove_flag(inverse, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::REQUIRE_NOT_ADJUNCT:
		return invert_require_no_flag(inverse, first, second, grammatical_flag::IS_ADJUNCT);
	case function_type::ADD_NULLABLE_SUBJECT:
		return invert_add_flag(inverse, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::REMOVE_NULLABLE_SUBJECT:
		return invert_remove_flag(inverse, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::TRY_REMOVE_NULLABLE_SUBJECT:
		return invert_try_remove_flag(inverse, first, second, grammatical_flag::NULLABLE_SUBJECT);
	case function_type::ADD_SUBORDINATE:
		return invert_add_flag(inverse, first, second, grammatical_flag::SUBORDINATE);
	case function_type::REMOVE_SUBORDINATE:
		return invert_remove_flag(inverse, first, second, grammatical_flag::SUBORDINATE);
	case function_type::TRY_REMOVE_SUBORDINATE:
		return invert_try_remove_flag(inverse, first, second, grammatical_flag::SUBORDINATE);
	case function_type::REQUIRE_NO_SUBORDINATE:
		return invert_require_no_flag(inverse, first, second, grammatical_flag::SUBORDINATE);
	case function_type::ADD_BOTH:
		return invert_add_correlator(inverse, first, second, correlator::BOTH);
	case function_type::ADD_EITHER:
		return invert_add_correlator(inverse, first, second, correlator::EITHER);
	case function_type::ADD_NEITHER:
		return invert_add_correlator(inverse, first, second, correlator::NEITHER);
	case function_type::REMOVE_BOTH:
		return invert_remove_correlator(inverse, first, second, correlator::BOTH);
	case function_type::REMOVE_EITHER:
		return invert_remove_correlator(inverse, first, second, correlator::EITHER);
	case function_type::REMOVE_NEITHER:
		return invert_remove_correlator(inverse, first, second, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATOR:
		if ((second.corr != correlator::NONE && second.corr != correlator::ANY)
		 || !intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
		 || !intersect(inverse.corr, first.corr, correlator::ANY))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::REQUIRE_NO_CORRELATOR:
		if ((second.corr != correlator::NONE && second.corr != correlator::ANY)
		 || !intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.correlated_by, first.correlated_by, second.correlated_by)
		 || !intersect(inverse.corr, first.corr, correlator::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::ADD_CORRELATED_BY_BOTH:
		return invert_add_correlated_by(inverse, first, second, correlator::BOTH);
	case function_type::ADD_CORRELATED_BY_EITHER:
		return invert_add_correlated_by(inverse, first, second, correlator::EITHER);
	case function_type::ADD_CORRELATED_BY_NEITHER:
		return invert_add_correlated_by(inverse, first, second, correlator::NEITHER);
	case function_type::TRY_REMOVE_CORRELATED:
		if ((second.correlated_by != correlator::NONE && second.correlated_by != correlator::ANY)
		 || !intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, correlator::ANY))
			return false;
		for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	case function_type::REQUIRE_NOT_CORRELATED:
		if ((second.correlated_by != correlator::NONE && second.correlated_by != correlator::ANY)
		 || !intersect(inverse.index_number, first.index_number, second.index_number)
		 || !intersect(inverse.concord_number, first.concord_number, second.concord_number)
		 || !intersect(inverse.cnj, first.cnj, second.cnj)
		 || !intersect(inverse.corr, first.corr, second.corr)
		 || !intersect(inverse.correlated_by, first.correlated_by, correlator::NONE))
			return false;
		for (i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
			if (!intersect(inverse.flags[i], first.flags[i], second.flags[i])) return false;
		return intersect(inverse.root, first.root, second.root);
	}
	fprintf(stderr, "invert ERROR: Unrecognized transformation function.\n");
	return false;
}

template<typename Formula>
bool invert(
	flagged_logical_form<Formula>*& inverse,
	unsigned int& inverse_count,
	typename flagged_logical_form<Formula>::function function,
	const flagged_logical_form<Formula>& first,
	const flagged_logical_form<Formula>& second)
{
	inverse_count = 1;
	inverse = (flagged_logical_form<Formula>*) malloc(sizeof(flagged_logical_form<Formula>));
	if (inverse == NULL) {
		fprintf(stderr, "invert ERROR: Out of memory.\n");
		return false;
	} else if (!invert(*inverse, function, first, second)) {
		free(inverse);
		return false;
	}
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
	exp.index_number = set.index_number;
	exp.concord_number = set.concord_number;
	exp.cnj = set.cnj;
	exp.corr = set.corr;
	exp.correlated_by = set.correlated_by;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		exp.flags[i] = set.flags[i];
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
	exp.index_number = set.index_number;
	exp.concord_number = set.concord_number;
	exp.cnj = set.cnj;
	exp.corr = set.corr;
	exp.correlated_by = set.correlated_by;
	for (uint_fast8_t i = 0; i < (uint_fast8_t) grammatical_flag::COUNT; i++)
		exp.flags[i] = set.flags[i];
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
