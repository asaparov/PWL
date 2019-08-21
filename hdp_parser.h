#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

template<typename K, typename V>
struct static_pair {
	K key;
	V value;
};

template<typename Formula>
struct flagged_logical_form {
	Formula* root;
	/* TODO: add features */

	flagged_logical_form() { }

	flagged_logical_form(Formula& src) : root(&src) {
		src.reference_count++;
	}

	~flagged_logical_form() { free_helper(); }

	inline void operator = (const flagged_logical_form<Formula>& src) {
		root = src.root;
		root->reference_count++;
	}

	static inline bool is_empty(const flagged_logical_form<Formula>& src) {
		return src.root == nullptr;
	}

	static inline unsigned int hash(const flagged_logical_form<Formula>& src) {
		return hasher<Formula>::hash(*src.root);
	}

	static inline void move(
			const flagged_logical_form<Formula>& src,
			flagged_logical_form<Formula>& dst)
	{
		dst.root = src.root;
	}

	static inline void swap(
			flagged_logical_form<Formula>& first,
			flagged_logical_form<Formula>& second)
	{
		core::swap(first.root, second.root);
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
		IDENTITY
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
	{function_type::IDENTITY, "identity"}
};

template<typename Formula>
inline bool initialize_any(flagged_logical_form<Formula>& lf) {
	lf.root = (Formula*) malloc(sizeof(Formula));
	if (lf.root == nullptr) {
		fprintf(stderr, "initialize_any ERROR: Out of memory.\n");
		return false;
	}
	initialize_any(*lf.root);
	return true;
}

template<typename Formula>
inline bool operator == (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	return (first.root == second.root)
		|| (*first.root == *second.root);
}

template<typename Formula>
inline bool operator != (
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	return (first.root != second.root)
		&& (*first.root != *second.root);
}

template<typename Formula, typename Stream, typename... Printer>
inline bool print(
		const flagged_logical_form<Formula>& lf,
		Stream& out, Printer&&... printer)
{
	/* TODO: print the features too */
	return print(*lf.root, out, std::forward<Printer>(printer)...);
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

template<typename Formula>
bool apply(typename flagged_logical_form<Formula>::function function,
		const flagged_logical_form<Formula>& src, flagged_logical_form<Formula>& dst)
{
	typedef typename flagged_logical_form<Formula>::function_type function_type;
	switch (function.type) {
	case function_type::EMPTY:
		dst.root = Formula::new_false();
		return dst.root != nullptr;
	case function_type::IDENTITY:
		dst.root = src.root;
		dst.root->reference_count++;
		return true;
	}
	fprintf(stderr, "apply ERROR: Unrecognized transformation function.\n");
	return false;
}

template<typename Formula>
inline bool invert(
		flagged_logical_form<Formula>& inverse,
		typename flagged_logical_form<Formula>::function function,
		const flagged_logical_form<Formula>& first,
		const flagged_logical_form<Formula>& second)
{
	typedef typename flagged_logical_form<Formula>::function_type function_type;
	flagged_logical_form<Formula> exp;
	switch (function.type) {
	case function_type::EMPTY:
		inverse.root = first.root;
		inverse.root->reference_count++;
		return true;
	case function_type::IDENTITY:
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
	/* TODO: copy features from `set` to `exp` */
	exp.root = (Formula*) malloc(sizeof(Formula));
	if (exp.root == nullptr) return false;
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
	/* TODO: copy features from `set` to `exp` */
	exp.root = (Formula*) malloc(sizeof(Formula));
	if (exp.root == nullptr) return false;
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
