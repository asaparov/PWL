#ifndef CONSOLE_H_
#define CONSOLE_H_

#include <core/lex.h>
#include "console_utils.h"
#include "higher_order_logic.h"
#include "natural_deduction.h"
#include "natural_deduction_mh.h"
#include "theory_prior.h"
#include "theory.h"
#include "executive.h"

template<typename Stream, typename Parser>
void run_console(
		Stream& input, const char* prompt, Parser& parser, hash_map<string, unsigned int>& names,
		array<array_map<typename Parser::SentenceType, typename Parser::logical_form_type>>& training_set)
{
	array<char> line = array<char>(256);
	while (true) {
		if (prompt) {
			printf("%s", prompt);
			fflush(stdout);
		}

		line.clear();
		if (!read_line(line, input)) {
			break;
		} else {
			if (!line.ensure_capacity(line.length + 1))
				break;
			line[line.length++] = '\0';
			parse_sentence(parser, line.data, names, training_set);
		}
	}
}

/*template<typename Stream, typename Parser>
inline void run_console(Stream& input, const char* prompt,
		Parser& parser, hash_map<string, unsigned int>& names)
{
	array<array_map<typename Parser::SentenceType, typename Parser::logical_form_type>> dummy_training_set(1);
	run_console(input, prompt, parser, names, dummy_training_set);
}*/

inline void print_help(unsigned int max_parse_count)
{
	printf(CONSOLE_BOLD "Available commands:\n" CONSOLE_RESET);
	printf(CONSOLE_BOLD "read" CONSOLE_RESET " [sentence without surrounding quotes]       Find up to %u logical forms that maximize the likelihood.\n", max_parse_count);
	printf(CONSOLE_BOLD "rerank" CONSOLE_RESET "                                           For each parsed logical form, compute its prior, and re-rank them according to their posterior.\n");
	printf(CONSOLE_BOLD "add" CONSOLE_RESET " [logical form 0-based index]                 Add the selected logical form as an observation to the theory. (default is the first logical form)\n");
	printf(CONSOLE_BOLD "mcmc" CONSOLE_RESET " [iterations]                                Perform MCMC for the specified number of iterations.\n");
	printf(CONSOLE_BOLD "print_theory" CONSOLE_RESET "                                     Print the current theory.\n");
	printf(CONSOLE_BOLD "print_proofs" CONSOLE_RESET "                                     Print the proof of each observation in the theory.\n");
	printf(CONSOLE_BOLD "answer" CONSOLE_RESET " [logical form 0-based index] [iterations] Try to answer the given question, using the specified number of MCMC iterations.\n");
	printf(CONSOLE_BOLD "research" CONSOLE_RESET " [question without surrounding quotes]   Try to answer the given question, searching the web for more information as needed.\n");
	printf(CONSOLE_BOLD "generate" CONSOLE_RESET " [logical form 0-based index]            Generate sentences from the selected logical form. (default is the first logical form)\n");
	printf(CONSOLE_BOLD "examples" CONSOLE_RESET "                                         Suggest some interesting examples.\n");
}

inline void print_examples()
{
	printf("Try the following to see how knowledge in the theory can be exploited to resolve syntactic ambiguity (PP-attachment):\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read Sally caught a butterfly with a net.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " rerank\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read No butterfly has a net.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " add 0\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read Sally caught a butterfly with a net.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " rerank\n\n");
	printf("Another example where the PP attaches to the noun with a \"softer\" constraint:\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read Sally caught a butterfly with a spot.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " rerank\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read The butterfly has a spot.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " add 0\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " read Sally caught a butterfly with a spot.\n");
	printf(" " CONSOLE_BOLD ">" CONSOLE_RESET " rerank\n\n");
}

template<typename Stream>
bool read_line(array<char>& line, Stream& input)
{
	while (true) {
		/* if `fgets` does not read a full line, then `line` will end with '\0' without a preceding newline */
		line[line.capacity - 1] = '1';
		if (fgets(line.data + line.length, (line.capacity - line.length) * sizeof(char), input) == nullptr)
			return false;

		if (line[line.capacity - 1] == '\0' && line[line.capacity - 2] != '\n') {
			/* we did not read a full line */
			line.length = line.capacity - 1;
			if (!line.ensure_capacity(line.capacity + 1)) {
				fprintf(stderr, "read_line ERROR: Out of memory.\n");
				return false;
			}
		} else {
			line.length += strlen(line.data + line.length) - 1;
			break;
		}
	}

	return true;
}

template<typename Stream, typename Parser>
void run_console(
		Stream& input, const char* prompt, Parser& parser,
		array<hol_term*>& seed_axioms,
		hash_map<string, unsigned int>& names)
{
	if (!parser.invert_name_map(names)) {
		fprintf(stderr, "ERROR: `run_console.invert_name_map` failed.\n");
		return;
	}

	theory<natural_deduction<hol_term, false>, polymorphic_canonicalizer<true, false, built_in_predicates>> T(seed_axioms, 1000000000);
	extern unsigned int constant_offset;
	constant_offset = T.new_constant_offset;
	auto constant_prior = make_simple_constant_distribution(
			iid_uniform_distribution<unsigned int>(100), chinese_restaurant_process<unsigned int>(1.0, 0.0),
			make_dirichlet_process(1.0e-12, make_dirichlet_process(1000.0, make_iid_uniform_distribution<hol_term>(10000))));
	auto theory_element_prior = make_simple_hol_term_distribution<built_in_predicates>(
						constant_prior, geometric_distribution(0.0001), very_light_tail_distribution(-40.0),
						0.0199999, 0.01, 0.0000001, 0.17, 0.1, 0.1, 0.01, 0.57, 0.01, 0.01,
						0.1099999, 0.01, 0.0000001, 0.1999999, 0.26, 0.01, 0.01, 0.0000001, 0.2, 0.2,
						0.999999998, 0.000000001, 0.000000001, 0.3, 0.4, 0.2, 0.4, -2000.0);
	auto axiom_prior = make_dirichlet_process(1.0e-1, theory_element_prior);
	auto conjunction_introduction_prior = uniform_subset_distribution<const nd_step<hol_term>*>(0.8);
	auto conjunction_elimination_prior = make_levy_process(poisson_distribution(2.0), poisson_distribution(1.0));
	auto universal_introduction_prior = unif_distribution<unsigned int>();
	auto universal_elimination_prior = chinese_restaurant_process<hol_term>(1.0, 0.0);
	auto term_indices_prior = make_levy_process(poisson_distribution(4.0), poisson_distribution(1.5));
	auto proof_prior = make_canonicalized_proof_prior(axiom_prior, conjunction_introduction_prior, conjunction_elimination_prior,
			universal_introduction_prior, universal_elimination_prior, term_indices_prior, poisson_distribution(20.0), 0.00001);

	typedef decltype(T) Theory;
	typedef decltype(proof_prior) ProofPrior;
	typedef typename ProofPrior::PriorState PriorStateType;

	PriorStateType proof_axioms;

	constexpr unsigned int max_parse_count = 4;
	unsigned int parse_count = 0;
	hol_term* logical_forms[max_parse_count];
	double log_likelihoods[max_parse_count];

	putc('\n', stdout);
	print_help(max_parse_count);

	array<char> line = array<char>(256);
	while (true) {
		if (prompt) {
			printf(CONSOLE_BOLD CONSOLE_BLUE "%s" CONSOLE_RESET, prompt);
			fflush(stdout);
		}

		line.clear();
		if (!read_line(line, input)) {
			break;
		} else {
			if (!line.ensure_capacity(line.length + 1))
				break;
			line[line.length] = '\0';

			unsigned int index = line.index_of(' ');

			/* find the command */
			if (compare_strings("print_theory", line.data, index)) {
				T.print_axioms(stdout, parser.get_printer());

				null_collector collector;
				array<hol_term*> extra_axioms(16);
				T.get_extra_axioms(extra_axioms);
				double value = log_probability(T.observations, extra_axioms, proof_prior, collector);
				print("Log probability of theory and proofs: ", stdout); print(value, stdout); print('\n', stdout);

			} else if (compare_strings("print_proofs", line.data, index)) {
				T.print_proofs(stdout, parser.get_printer());

				null_collector collector;
				array<hol_term*> extra_axioms(16);
				T.get_extra_axioms(extra_axioms);
				double value = log_probability(T.observations, extra_axioms, proof_prior, collector);
				print("Log probability of theory and proofs: ", stdout); print(value, stdout); print('\n', stdout);

			} else if (compare_strings("read", line.data, index)) {
				while (isspace(line[index])) index++;
				if (index == line.length) {
					printf("ERROR: Missing input sentence.\n");
					continue;
				}

				for (unsigned int i = 0; i < parse_count; i++) {
					free(*logical_forms[i]); if (logical_forms[i]->reference_count == 0) free(logical_forms[i]);
				}
				if (!parse_sentence(parser, line.data + index, names, logical_forms, log_likelihoods, parse_count)) {
					printf("Unable to parse sentence \"%s\"", line.data + index);
					parse_count = 0;
				}

			} else if (compare_strings("research", line.data, index)) {
				while (isspace(line[index])) index++;
				if (index == line.length) {
					printf("ERROR: Missing input sentence.\n");
					continue;
				}

				array<string> answers(4);
				hash_set<unsigned int> seed_entities(2);
				const in_memory_article_store<typename Parser::DerivationType> corpus;
				if (answer_question<true>(answers, line.data + index, 1000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
					print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
				}
				for (string& str : answers)
					free(str);

			} else if (compare_strings("answer", line.data, index)) {
				if (parse_count == 0) {
					printf("ERROR: There are no logical forms to answer. Use 'read' to parse a sentence into logical forms.\n");
					continue;
				}

				unsigned int logical_form_index = 0;
				while (isspace(line[index])) index++;
				if (index < line.length) {
					char* end_ptr;
					logical_form_index = strtoul(line.data + index, &end_ptr, 0);
					if (isspace(*end_ptr)) {
						index = ((size_t) end_ptr - (size_t) line.data) / sizeof(char);
					} else if (*end_ptr != '\0') {
						printf("ERROR: Invalid logical form index.\n");
						continue;
					}
				}
				if (logical_form_index >= parse_count) {
					printf("ERROR: Logical form index is out of bounds.\n");
					continue;
				}

				unsigned int mcmc_iterations = 400;
				while (isspace(line[index])) index++;
				if (index < line.length) {
					char* end_ptr;
					mcmc_iterations = strtoul(line.data + index, &end_ptr, 0);
					if (*end_ptr != '\0') {
						printf("ERROR: Invalid number of MCMC iterations.\n");
						continue;
					}
				}


				hash_set<unsigned int> seed_entities(2);
				array_map<string, double> answers(16);
				if (answer_question<false>(answers, logical_forms[logical_form_index], mcmc_iterations, parser.get_printer(), T, proof_prior, proof_axioms)) {
					print("Answers:\n", stdout);
					sort(answers.values, answers.keys, answers.size, default_sorter());
					for (unsigned int i = answers.size; i > 0; i--) {
						print("  ", stdout); print(answers.keys[i - 1], stdout);
						print(" with log probability ", stdout); print(answers.values[i - 1], stdout);
						print('\n', stdout);
					}
				}
				for (auto entry : answers)
					free(entry.key);

			} else if (compare_strings("rerank", line.data, index)) {
				if (parse_count == 0) {
					printf("ERROR: There are no logical forms to rerank. Use 'read' to parse a sentence into logical forms.\n");
					continue;
				}

				double log_priors[max_parse_count];
				double log_posteriors[max_parse_count];
				for (unsigned int i = 0; i < parse_count; i++) {
					log_priors[i] = log_joint_probability_of_observation(T, proof_prior, proof_axioms, logical_forms[i], 400);
					log_posteriors[i] = log_likelihoods[i] + log_priors[i];
				}
				unsigned int sorted_indices[max_parse_count];
				for (unsigned int i = 0; i < parse_count; i++)
					sorted_indices[i] = i;
				if (parse_count > 1) {
					insertion_sort(log_posteriors, sorted_indices, parse_count);
					reverse(log_posteriors, parse_count);
					reverse(sorted_indices, parse_count);
				}

				hol_term* temp_logical_forms[max_parse_count];
				double temp_log_likelihoods[max_parse_count];
				double temp_log_priors[max_parse_count];
				for (unsigned int i = 0; i < parse_count; i++) {
					temp_logical_forms[i] = logical_forms[i];
					temp_log_likelihoods[i] = log_likelihoods[i];
					temp_log_priors[i] = log_priors[i];
				} for (unsigned int i = 0; i < parse_count; i++) {
					logical_forms[i] = temp_logical_forms[sorted_indices[i]];
					log_likelihoods[i] = temp_log_likelihoods[sorted_indices[i]];
					log_priors[i] = temp_log_priors[sorted_indices[i]];
				}

				for (unsigned int i = 0; i < parse_count; i++) {
					print(CONSOLE_BOLD "[", stdout); print(i, stdout); print("] " CONSOLE_RESET, stdout);
					print(*logical_forms[i], stdout, parser.get_printer());
					print(" with log likelihood ", stdout); print(log_likelihoods[i], stdout);
					print(" + log prior ", stdout); print(log_priors[i], stdout);
					print(" = log posterior ", stdout); print(log_posteriors[i], stdout);
					print('\n', stdout);
				}

			} else if (compare_strings("add", line.data, index)) {
				if (parse_count == 0) {
					printf("ERROR: There are no logical forms to add. Use 'read' to parse a sentence into logical forms.\n");
					continue;
				}

				unsigned int logical_form_index = 0;
				while (isspace(line[index])) index++;
				if (index < line.length) {
					char* end_ptr;
					logical_form_index = strtoul(line.data + index, &end_ptr, 0);
					if (*end_ptr != '\0') {
						printf("ERROR: Invalid logical form index.\n");
						continue;
					}
				}
				if (logical_form_index >= parse_count) {
					printf("ERROR: Logical form index is out of bounds.\n");
					continue;
				}

				set_changes<hol_term> set_diff;
				unsigned int new_constant;
				auto* new_proof = T.add_formula(logical_forms[logical_form_index], set_diff, new_constant);
				for (unsigned int i = 0; new_proof == nullptr && i < 100; i++) {
					set_diff.clear();
					null_collector collector;
					for (unsigned int t = 0; t < 10; t++)
						do_exploratory_mh_step(T, proof_prior, proof_axioms, collector);
					new_proof = T.add_formula(logical_forms[logical_form_index], set_diff, new_constant);
				}
				if (new_proof != nullptr && !proof_axioms.add(new_proof, set_diff.new_set_axioms, proof_prior)) {
					T.remove_formula(new_proof, set_diff);
					new_proof = nullptr;
				}
				if (new_proof == nullptr) {
					print("ERROR: Unable to add logical form to theory.\n", stdout);
					print("  Logical form: ", stdout); print(*logical_forms[logical_form_index], stdout, parser.get_printer()); print("\n", stdout);
					continue;
				}
				print("Successfully added logical form as an observation to the theory. Use 'print_theory' to inspect the current theory. Note that no additional MCMC iterations were performed.\n", stdout);
				print("  Logical form: ", stdout); print(*logical_forms[logical_form_index], stdout, parser.get_printer()); print("\n", stdout);

				for (unsigned int i = 0; i < parse_count; i++) {
					free(*logical_forms[i]); if (logical_forms[i]->reference_count == 0) free(logical_forms[i]);
				}
				parse_count = 0;

            } else if (compare_strings("mcmc", line.data, index)) {
				while (isspace(line[index])) index++;
				if (index == line.length) {
					printf("Missing number of MCMC iterations argument.\n");
					continue;
				}

                char* end_ptr;
                unsigned int mcmc_iterations = strtoul(line.data + index, &end_ptr, 0);
                if (*end_ptr != '\0') {
                    printf("ERROR: Invalid number of MCMC iterations.\n");
                    continue;
                }

				Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
				PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
				hash_map<const hol_term*, hol_term*> formula_map(128);
				Theory::clone(T, T_MAP, formula_map);
				new (&proof_axioms_MAP) PriorStateType(proof_axioms, formula_map);
				auto collector = make_log_probability_collector(T, proof_prior);
				double max_log_probability = collector.current_log_probability;
				for (unsigned int j = 0; j < mcmc_iterations; j++) {
					do_mh_step(T, proof_prior, proof_axioms, collector);
					if (collector.current_log_probability > max_log_probability) {
						free(T_MAP); proof_axioms_MAP.~PriorStateType(); formula_map.clear();
						Theory::clone(T, T_MAP, formula_map);
						new (&proof_axioms_MAP) PriorStateType(proof_axioms, formula_map);
						max_log_probability = collector.current_log_probability;
					}
				}
				printf("Maximum a posteriori theory:\n");
				T_MAP.print_axioms(stderr, parser.get_printer());
				free(T_MAP); proof_axioms_MAP.~PriorStateType();

			} else if (compare_strings("generate", line.data, index)) {
				if (parse_count == 0) {
					printf("ERROR: There are no logical forms from which to generate. Use 'read' to parse a sentence into logical forms.\n");
					continue;
				}

				unsigned int logical_form_index = 0;
				while (isspace(line[index])) index++;
				if (index < line.length) {
					char* end_ptr;
					logical_form_index = strtoul(line.data + index, &end_ptr, 0);
					if (*end_ptr != '\0') {
						printf("ERROR: Invalid logical form index.\n");
						continue;
					}
				}
				if (logical_form_index >= parse_count) {
					printf("ERROR: Logical form index is out of bounds.\n");
					continue;
				}

				unsigned int generated_derivation_count;
				constexpr unsigned int max_generated_derivation_count = 4;
				double log_likelihoods[max_generated_derivation_count];
				syntax_node<typename Parser::logical_form_type>* generated_derivations =
						(syntax_node<typename Parser::logical_form_type>*) alloca(sizeof(syntax_node<typename Parser::logical_form_type>) * max_generated_derivation_count);
				if (!parser.template generate<max_generated_derivation_count>(generated_derivations, log_likelihoods, generated_derivation_count, logical_forms[logical_form_index], names) || generated_derivation_count == 0) {
					printf("ERROR: Failed to generate derivation.\n");
					continue;
				}

				const string** nonterminal_name_map = invert(parser.G.nonterminal_names);
				string_map_scribe nonterminal_printer = { nonterminal_name_map, parser.G.nonterminal_names.table.size + 1 };
				printf("Generated %u derivations:\n", generated_derivation_count);
				for (unsigned int i = 0; i < generated_derivation_count; i++) {
					/* compute the yield of the derivation */
					sequence new_sentence = sequence(NULL, 0);
					if (!parser.yield(generated_derivations[i], logical_forms[0], new_sentence))
						break;
					print('"', stdout); print(new_sentence, stdout, parser.get_printer());
					print("\" with log likelihood ", stdout); print(log_likelihoods[i], stdout);
					print(" and derivation tree:\n", stdout);
					print(generated_derivations[i], stdout, nonterminal_printer, parser.get_printer());
					print('\n', stdout);
					free(new_sentence);
				}
				free(nonterminal_name_map);
				for (unsigned int i = 0; i < generated_derivation_count; i++)
					free(generated_derivations[i]);

			} else if (compare_strings("help", line.data, index)) {
				print_help(max_parse_count);

			} else if (compare_strings("examples", line.data, index)) {
				print_examples();

			} else {
				line[index] = '\0';
				printf("%s: Command not found.\n", line.data);
			}
		}
	}

	for (unsigned int i = 0; i < parse_count; i++) {
		free(*logical_forms[i]); if (logical_forms[i]->reference_count == 0) free(logical_forms[i]);
	}
}

#endif /* CONSOLE_H_ */
