#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

template<typename Formula>
struct logical_form_with_features {
	Formula root;
	/* TODO: add features */

	logical_form_with_features(const Formula& src) : root(src) { }
};

typedef sequence_distribution<token_distribution<double>> terminal_prior_type;
template<typename Formula> using hdp_grammar_type =
	hdp_grammar<rule_list_prior<terminal_prior<terminal_prior_type>, logical_form_with_features<Formula>>, logical_form_with_features<Formula>>;

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
	const syntax_node<logical_form_with_features<Formula>>* const* const* syntax,
	const array<array_map<sentence, Formula>>& data,
	const string** reverse_name_map)
{
	typedef logical_form_with_features<Formula> logical_form_type;

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

inline sequence to_sequence(const sentence& src)
{
	sequence.
}

template<typename Formula>
struct hdp_parser
{
	typedef logical_form_with_features<Formula> logical_form_type;

	hdp_grammar_type<Formula> G;

	hdp_parser(unsigned int unknown_id,
			hash_map<string, unsigned int>& names,
			const char* grammar_filepath)
	{
		if (!read_grammar(G, names, grammar_filepath)) {
			fprintf(stderr, "ERROR: Unable to read grammar at '%s'.\n", grammar_filepath);
			exit(EXIT_FAILURE);
		}
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
			if (syntax != NULL) free(syntax);
			free(reverse_name_map); free(nonterminal_name_map);
			return false;
		}
		for (unsigned int i = 0; i < data.length; i++) {
			syntax[i] = (syntax_node<logical_form_type>**) calloc(data[i].size, sizeof(syntax_node<logical_form_type>*));
			if (syntax[i] == NULL) {
				fprintf(stderr, "hdp_parser.train ERROR: Out of memory.\n");
				for (unsigned int k = 0; k < i; k++) free(syntax[k]);
				free(syntax); free(order);
				free(reverse_name_map); free(nonterminal_name_map);
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
				sequence seq = to_sequence(data[id].keys[j]);
				auto sentence = tokenized_sentence<logical_form_type>(sequence(data[id].keys[j].tokens, data[id].keys[j].length));
				syntax[id][j] = (syntax_node<logical_form_type>*) malloc(sizeof(syntax_node<logical_form_type>));
				/* NOTE: sample can set syntax[id] to null */
				if (syntax[id][j] == NULL || !sample(syntax[id][j], G, logical_form, sentence, reverse_name_map) || syntax[id][j] == NULL)
				{
					fprintf(stderr, "sample ERROR: Unable to sample derivation for example %u, sentence %u: '", id, j);
					print(data[id].keys[j], stderr, terminal_printer); print("'\n", stderr);
					print(logical_form, stderr, terminal_printer); print("\n", stderr);
					for (unsigned int k = 0; k < data.length; k++) {
						for (unsigned int l = 0; l < data[k].size; l++) free(syntax[k][l]);
						free(syntax[k]);
					}
					free(syntax); free(order);
					free(reverse_name_map); free(nonterminal_name_map);
					return false;
				}

				print(logical_form, stdout, terminal_printer); print('\n', stdout);
				print(*syntax[id][j], stdout, nonterminal_printer, terminal_printer); print("\n\n", stdout);
				fflush(stdout);

				if (!add_tree(1, *syntax[id][j], logical_form, G)) {
					for (unsigned int k = 0; k < data.length; k++) {
						for (unsigned int l = 0; l < data[k].size; l++) free(syntax[k][l]);
						free(syntax[k]);
					}
					free(syntax); free(order);
					free(reverse_name_map); free(nonterminal_name_map);
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
					auto sentence = tokenized_sentence<logical_form_type>(sequence(data[id].keys[j].tokens, data[id].keys[j].length));
					resample(syntax[id][j], G, logical_form, sentence, reverse_name_map);
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
		for (unsigned int k = 0; k < data.length; k++) {
			for (unsigned int l = 0; l < data[k].size; l++) free(syntax[k][l]);
			free(syntax[k]);
		}
		free(syntax); free(order);
		free(reverse_name_map); free(nonterminal_name_map);
		return true;
	}

	template<unsigned int K, typename TheoryType>
	bool parse(const sentence& s, Formula** logical_forms,
			double* log_probabilities, unsigned int& parse_count,
			const TheoryType& T, array<token>& unrecognized) const
	{
		static_assert(K > 0, "`K` must be at least 1.");

		logical_form_type logical_form; /* TODO: initialize the features */
		syntax_node<logical_form_type>* parsed_syntax =
				(syntax_node<logical_form_type>*) alloca(K * sizeof(syntax_node<logical_form_type>));
		logical_form_type* logical_form_output =
				(logical_form_type*) alloca(K * sizeof(logical_form_type));
		auto sentence = tokenized_sentence<logical_form_type>(sequence(s.tokens, s.length));

		if (!parse<false, false, K>(parsed_syntax, parse_count, logical_form,
				logical_form_output, G, sentence, terminal_printer.map))
			return false;

		for (unsigned int i = 0; i < parse_count; i++) {
			double log_likelihood = log_probability(G, parsed_syntax[i], logical_form_output[i], terminal_printer.map);
			double log_prior = log_probability<true>(T, logical_form_output[i]);
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
			} else if (!init(*logical_forms[i], logical_form_output[i].root)) {
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
};

#endif /* HDP_PARSER_H_ */
