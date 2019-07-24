#ifndef HDP_PARSER_H_
#define HDP_PARSER_H_

#include "higher_order_logic.h"
#include <grammar/parser.h>
#include <grammar/hdp_grammar_io.h>

template<typename Formula>
struct logical_form_with_features {
	Formula root;
	/* TODO: add features */
};

typedef sequence_distribution<token_distribution<double>> terminal_prior_type;
template<typename Formula> using hdp_grammar_type = hdp_grammar<rule_list_prior<terminal_prior<terminal_prior_type>, logical_form_with_features<Formula>>, logical_form_with_features<Formula>>;

template<typename Formula>
struct hdp_parser
{
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

	template<unsigned int K, typename TheoryType>
	bool parse(const sentence& s, Formula** logical_forms,
			double* log_probabilities, unsigned int& parse_count,
			const TheoryType& T, array<token>& unrecognized) const
	{
		static_assert(K > 0, "`K` must be at least 1.");

		logical_form_with_features<Formula> logical_form; /* TODO: initialize the features */
		syntax_node<logical_form_with_features<Formula>>* parsed_syntax =
				(syntax_node<logical_form_with_features<Formula>>*) alloca(K * sizeof(syntax_node<logical_form_with_features<Formula>>));
		logical_form_with_features<Formula>* logical_form_output =
				(logical_form_with_features<Formula>*) alloca(K * sizeof(logical_form_with_features<Formula>));
		auto sentence = tokenized_sentence<logical_form_with_features<Formula>>(s);

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
