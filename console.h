#ifndef CONSOLE_H_
#define CONSOLE_H_

#include <core/lex.h>

template<typename Stream, typename Parser>
void run_console(
		Stream& input, const char* prompt, Parser& parser, hash_map<string, unsigned int>& names,
		array<array_map<typename Parser::SentenceType, flagged_logical_form<hol_term>>>& training_set)
{
	array<char> line = array<char>(256);
	while (true) {
		if (prompt) {
			printf("%s", prompt);
			fflush(stdout);
		}

		line.clear();
		int read = read_line(line, input);
		if (read == 0) {
			break;
		} else if (read > 0) {
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
	array<array_map<typename Parser::SentenceType, flagged_logical_form<hol_term>>> dummy_training_set(1);
	run_console(input, prompt, parser, names, dummy_training_set);
}*/

template<typename Stream, typename Parser>
void run_console(
		Stream& input, const char* prompt, Parser& parser, hash_map<string, unsigned int>& names)
{
	theory<natural_deduction<hol_term, true>, polymorphic_canonicalizer<true, false, built_in_predicates>> T(1000000000);
	constant_offset = T.new_constant_offset;
	auto constant_prior = make_simple_constant_distribution(
			iid_uniform_distribution<unsigned int>(100), chinese_restaurant_process<unsigned int>(1.0, 0.0),
			make_dirichlet_process(1.0e-12, make_dirichlet_process(1000.0, make_iid_uniform_distribution<hol_term>(10000))));
	auto theory_element_prior = make_simple_hol_term_distribution<built_in_predicates>(
						constant_prior, geometric_distribution(0.02), very_light_tail_distribution(-40.0),
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

	array<char> line = array<char>(256);
	while (true) {
		if (prompt) {
			printf("%s", prompt);
			fflush(stdout);
		}

		line.clear();
		int read = read_line(line, input);
		if (read == 0) {
			break;
		} else if (read > 0) {
			if (!line.ensure_capacity(line.length + 1))
				break;
			line[line.length++] = '\0';

			unsigned int index = line.index_of(' ');

			/* find the command */
			if (compare_strings("print_theory", line.data, index)) {
				T.print_axioms(stdout, *debug_terminal_printer);

			} else if (compare_strings("print_proofs", line.data, index)) {
				T.print_disjunction_introductions(stdout, *debug_terminal_printer);

			} else if (compare_strings("read", line.data, index)) {
				if (index == line.length) {
					printf("ERROR: Missing input sentence.\n");
					continue;
				}

				for (unsigned int i = 0; i < parse_count; i++) {
					free(*logical_forms[i]); if (logical_forms[i]->reference_count == 0) free(logical_forms[i]);
				}
				if (!parse_sentence(parser, line.data + index + 1, names, logical_forms, log_likelihoods, parse_count)) {
					printf("Failed to parse sentence \"%s\"", line.data + index + 1);
					parse_count = 0;
				}

			} else if (compare_strings("rerank", line.data, index)) {
				if (parse_count == 0) {
					printf("ERROR: There are no logical forms to rerank. Use 'read' to parse a sentence into logical forms.\n");
					continue;
				}

				double log_priors[max_parse_count];
				double log_posteriors[max_parse_count];
				for (unsigned int i = 0; i < parse_count; i++) {
					log_priors[i] = log_joint_probability_of_observation(T, proof_prior, proof_axioms, logical_forms[i], 1000);
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
				} else if (index == line.length) {
					printf("ERROR: Missing logical form index argument.\n");
					continue;
				}

				char* end_ptr;
				unsigned int logical_form_index = strtoul(line.data + index + 1, &end_ptr, 0);
				if (*end_ptr != '\0') {
					printf("ERROR: Invalid logical form index.\n");
					continue;
				} else if (logical_form_index >= parse_count) {
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
					print("read_sentence ERROR: Unable to add logical form to theory.\n", stdout);
					print("  Logical form: ", stdout); print(*logical_forms[logical_form_index], stdout, parser.get_printer()); print("\n", stdout);
					continue;
				}
				print("Successfully added logical form as an observation to the theory. Use 'print_theory' to inspect the current theory. Note that no additional MCMC iterations were performed.", stdout);
				print("  Logical form: ", stdout); print(*logical_forms[logical_form_index], stdout, parser.get_printer()); print("\n", stdout);

            } else if (compare_strings("mcmc", line.data, index)) {
				if (index == line.length) {
					printf("Missing number of MCMC iterations argument.\n");
					continue;
				}

                char* end_ptr;
                unsigned int mcmc_iterations = strtoul(line.data + index + 1, &end_ptr, 0);
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
				T_MAP.print_axioms(stderr, *debug_terminal_printer);
				free(T_MAP); proof_axioms_MAP.~PriorStateType();

			} else if (compare_strings("help", line.data, index)) {
				printf("Available commands:\n");
				printf("read <sentence without surrounding quotes>      Find up to %u logical forms that maximize the likelihood.\n", max_parse_count);
				printf("rerank                                          For each parsed logical form, compute its prior, and re-rank them according to their posterior.\n");
                printf("add <logical form 0-based index>                Add the selected logical form as an observation to the theory.\n");
				printf("mcmc <iterations>                               Perform MCMC for the specified number of iterations.\n");
				printf("print_theory                                    Print the current theory.\n");
				printf("print_proofs                                    Print the proof of each observation in the theory.\n");

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
