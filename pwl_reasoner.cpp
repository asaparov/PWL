#include <core/utility.h>
const thread_local core::string_map_scribe* debug_terminal_printer = nullptr;
bool debug_flag = false;
unsigned int debug_counter = 0;

#include "theory_prior.h"
#include "built_in_predicates.h"
#include "natural_deduction.h"
#include "natural_deduction_mh.h"
#include "theory.h"
#include "executive.h"

#include <locale.h>
#include <cstdlib>


unsigned int constant_offset = 0;

template<typename Stream>
bool print_special_string(unsigned int key, Stream& out) {
	return print('c', out) && print_subscript(key - constant_offset, out);
}

template<typename Stream>
bool read_terms(
	array<hol_term*>& terms, Stream& in,
	hash_map<string, unsigned int>& names)
{
	array<tptp_token> tokens = array<tptp_token>(512);
	if (!tptp_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis failed.\n");
		free_tokens(tokens); return false;
	}

	unsigned int index = 0;
	while (index < tokens.length) {
		array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
		hol_term* term = (hol_term*) malloc(sizeof(hol_term));
		if (term == NULL) {
			fprintf(stderr, "read_terms ERROR: Out of memory.\n");
			free_tokens(tokens); return false;
		} else if (!tptp_interpret(tokens, index, *term, names, variables)) {
			fprintf(stderr, "ERROR: Unable to parse higher-order term.\n");
			for (auto entry : variables) free(entry.key);
			free(term); free_tokens(tokens); return false;
		} else if (!expect_token(tokens, index, tptp_token_type::SEMICOLON, "semicolon at end of higher-order term") || !terms.add(term)) {
			free(*term); free(term);
			free_tokens(tokens); return false;
		}
		index++;

		if (variables.size != 0)
			fprintf(stderr, "WARNING: Variable map is not empty.\n");
	}
	free_tokens(tokens);
	return true;
}

template<typename Theory, typename Collector>
inline bool observations_has_test_proof(const Theory& T, const Collector& collector) {
	return true;
}

template<typename Theory, typename ProofCalculus, typename Canonicalizer>
inline bool observations_has_test_proof(const Theory& T, const log_probability_collector<ProofCalculus, Canonicalizer>& collector) {
	return T.observations.contains(collector.test_proof);
}

template<typename Theory, typename ProofAxioms, typename Collector>
bool check_consistency(
	Theory& T, const ProofAxioms& proof_axioms, const Collector& collector,
	unsigned int i, unsigned int j, unsigned int t, const char* phase)
{
	bool success = proof_axioms.check_proof_axioms(T);
	success &= proof_axioms.check_universal_eliminations(T, collector);
	success &= T.check_concept_axioms();
	success &= T.check_disjunction_introductions();
	success &= T.are_elements_provable(*debug_terminal_printer);
	success &= T.sets.check_freeable_sets();
	success &= T.sets.are_descendants_valid();
	success &= T.sets.are_set_sizes_valid();
	success &= T.sets.check_set_ids();
	success &= T.sets.check_symbols_in_formulas();
	success &= T.sets.are_provable_elements_valid();
	if (T.ctx.referent_iterators.length != T.observations.length) {
		fprintf(stderr, "WARNING: `T.ctx.referent_iterators.length` is not equal to `T.observations.length`.\n");
		success = false;
	}
	if (!observations_has_test_proof(T, collector)) {
		fprintf(stderr, "WARNING: `collector.test_proof` is not an observation in the theory.\n");
		success = false;
	}
	if (!success) {
		T.template print_axioms<true>(stderr, *debug_terminal_printer);
		T.print_disjunction_introductions(stderr, *debug_terminal_printer);
		fprintf(stderr, "ERROR: Theory consistency check failed at iteration: i = %u, j = %u, t = %u", i, j, t);
		fprintf(stderr, " (phase: %s)\n", phase);
		return false;
	}
	return true;
}

inline void print_usage()
{
	fprintf(stdout, "Usage: pwl_reasoner <file with logical forms> [r] [n] [h] [N] [M]\n"
					"  r = Number of retry attempts after failing to add LF (default: 100)\n"
					"  n = Number of high-temp MCMC steps for each retry attempt (defualt: 10)\n"
					"  h = Number of reheating phases after adding each LF (default: 4)\n"
					"  N = Number of MCMC iterations per reheating phase (default: 500)\n"
					"  M = Number of high-temp MCMC iterations between each reheating phase (default: 40)\n"
		);
	fflush(stdout);
}

int main(int argc, const char** argv)
{
#if defined(_WIN32)
	/* set the terminal to correctly display UTF-8 characters */
	SetConsoleOutputCP(CP_UTF8);
#endif
	setlocale(LC_ALL, "en_US.UTF-8");
	log_cache<double>::instance().ensure_size(1024);
set_seed(1356941742);
	fprintf(stdout, "(seed = %u)\n", get_seed());

	/* parse command-line arguments */
	if (argc < 2) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	unsigned int retry_attempts = 100;
	unsigned int iterations_per_retry = 10;
	unsigned int num_reheating_phases = 4;
	unsigned int iterations_per_reheat = 800;
	unsigned int reheating_iterations = 200;

	/* parse the optional command-line arguments */
	if (argc > 2 && !parse_uint(string(argv[2]), retry_attempts)) {
		print_usage();
		exit(EXIT_FAILURE);
	} if (argc > 3 && !parse_uint(string(argv[3]), iterations_per_retry)) {
		print_usage();
		exit(EXIT_FAILURE);
	} if (argc > 4 && !parse_uint(string(argv[4]), num_reheating_phases)) {
		print_usage();
		exit(EXIT_FAILURE);
	} if (argc > 5 && !parse_uint(string(argv[5]), iterations_per_reheat)) {
		print_usage();
		exit(EXIT_FAILURE);
	} if (argc > 6 && !parse_uint(string(argv[6]), reheating_iterations)) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	hash_map<string, unsigned int> names(256);
	if (!add_constants_to_string_map(names))
		return EXIT_FAILURE;

	/* read the seed axioms */
	array<hol_term*> seed_axioms(8);
	const char* axioms_filename = "seed_axioms.txt";
	FILE* in = fopen(axioms_filename, "rb");
	if (in == nullptr) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", axioms_filename);
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	} else if (!read_terms(seed_axioms, in, names)) {
		fprintf(stderr, "ERROR: Failed to parse logical forms in '%s'.\n", axioms_filename);
		fclose(in); free_all(seed_axioms);
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}
	printf("Finished reading seed axioms.\n"); fflush(stdout);
	fclose(in);

	/* read the input logical forms */
	array<hol_term*> lfs(8);
	const char* input_filename = argv[1];
	in = fopen(input_filename, "rb");
	if (in == nullptr) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", input_filename);
		free_all(seed_axioms);
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	} else if (!read_terms(lfs, in, names)) {
		fprintf(stderr, "ERROR: Failed to parse logical forms in '%s'.\n", input_filename);
		fclose(in); free_all(seed_axioms); free_all(lfs);
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}
	printf("Finished reading input logical forms.\n"); fflush(stdout);
	fclose(in);

	const string** name_map = invert(names);
	if (name_map == nullptr) {
		free_all(lfs);
		free_all(seed_axioms);
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}
	string_map_scribe printer = {name_map, names.table.size + 1};
	debug_terminal_printer = &printer;

	/* initialize the theory */
	theory<natural_deduction<hol_term, false>, polymorphic_canonicalizer<true, false, built_in_predicates>> T(seed_axioms, 1000000000);
	printf("Finished constructing initial theory.\n"); fflush(stdout);

	constant_offset = T.new_constant_offset;
	auto constant_prior = make_simple_constant_distribution(
			iid_uniform_distribution<unsigned int>(10000), chinese_restaurant_process<unsigned int>(1.0, 0.0),
			make_dirichlet_process(1.0e-1, make_dirichlet_process(1000.0, make_iid_uniform_distribution<hol_term>(10000))));
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

	typedef decltype(proof_prior) ProofPrior;
	typedef typename ProofPrior::PriorState PriorStateType;
	typedef decltype(T) Theory;

	PriorStateType proof_axioms;

	for (unsigned int i = 0; i + 1 < lfs.length; i++) {
		/* add the logical form to the theory */
		hol_term* lf = lfs[i];
		set_changes<hol_term> set_diff;
		unsigned int new_constant;
		auto* new_proof = T.add_formula(lf, set_diff, new_constant);
		for (unsigned int j = 0; new_proof == nullptr && j < retry_attempts; j++) {
			set_diff.clear();
			null_collector collector;
			for (unsigned int t = 0; t < iterations_per_retry; t++) {
//fprintf(stderr, "(add_formula) j = %u, t = %u\n", j, t);
//if (!check_consistency(T, proof_axioms, collector, i, j, t, "add_formula")) exit(EXIT_FAILURE);
/*T.template print_axioms<true>(stdout, *debug_terminal_printer);
T.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout);
if (i == 13 && j == 76 && t == 2) {
fprintf(stderr, "DEBUG\n");
debug_flag = true;
}*/
				do_exploratory_mh_step(T, proof_prior, proof_axioms, collector);
			}
			new_proof = T.add_formula(lf, set_diff, new_constant);
		}
		array<hol_term*> new_set_axioms(8);
		for (hol_term* formula : set_diff.new_set_axioms) {
			bool is_old_formula = false;
			for (hol_term* old_formula : set_diff.old_set_axioms) {
				if (*formula == *old_formula) {
					is_old_formula = true;
					break;
				}
			}
			if (!is_old_formula)
				new_set_axioms.add(formula);
		}
		if (new_proof != nullptr && !proof_axioms.add(new_proof, new_set_axioms, proof_prior)) {
			T.remove_formula(new_proof, set_diff);
			new_proof = nullptr;
		}
		if (new_proof == nullptr) {
			print("ERROR: Unable to add logical form to theory.\n", stdout);
			print("  Logical form: ", stdout); print(*lf, stdout, printer); print("\n", stdout);
			continue;
		}
		print("Successfully added logical form as an observation to the theory.\n", stdout);
		print("  Logical form: ", stdout); print(*lf, stdout, printer); print("\n", stdout);

		print("Running some iterations of MCMC to optimize theory...\n", stdout);
		Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
		PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
		hash_map<const hol_term*, hol_term*> formula_map(128);
		Theory::clone(T, T_MAP, formula_map);
		PriorStateType::clone(proof_axioms, proof_axioms_MAP, formula_map);
		auto collector = make_log_probability_collector(T, proof_prior, new_proof);
		double max_log_probability = collector.current_log_probability;
		for (unsigned int j = 0; j < num_reheating_phases; j++) {
			for (unsigned int t = 0; t < iterations_per_reheat; t++) {
//fprintf(stderr, "i = %u, j = %u, t = %u\n", i, j, t);
//if (!check_consistency(T, proof_axioms, collector, i, j, t, "intermediate MCMC")) exit(EXIT_FAILURE);
/*T.template print_axioms<true>(stdout, *debug_terminal_printer);
T.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout);
if (i == 10 && j == 0 && t == 4) {
fprintf(stderr, "DEBUG\n");
debug_flag = true;
} else {
debug_flag = false;
}*/
				do_mh_step(T, proof_prior, proof_axioms, collector, collector.test_proof, (t < 40 ? 1.0 : 0.01));

				if (collector.current_log_probability > max_log_probability + 1.0e-8) {
					free(T_MAP); free(proof_axioms_MAP); formula_map.clear();
					Theory::clone(T, T_MAP, formula_map);
					PriorStateType::clone(proof_axioms, proof_axioms_MAP, formula_map);
					max_log_probability = collector.current_log_probability;
				}
			}

			if (j + 1 < num_reheating_phases) {
				for (unsigned int t = 0; t < reheating_iterations; t++) {
//fprintf(stderr, "i = %u, j = %u, t = %u\n", i, j, t);
//if (!check_consistency(T, proof_axioms, collector, i, j, t, "exploratory")) exit(EXIT_FAILURE);
//T.template print_axioms<true>(stdout, *debug_terminal_printer);
//T.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout);
					do_exploratory_mh_step(T, proof_prior, proof_axioms, collector, collector.test_proof, 1.0);
				}
			}
		}
		free(T); free(proof_axioms); formula_map.clear();
		Theory::clone(T_MAP, T, formula_map);
		PriorStateType::clone(proof_axioms_MAP, proof_axioms, formula_map);
		print("Best theory so far:\n", stdout);
		T_MAP.template print_axioms<true>(stdout, *debug_terminal_printer);
		print("Theory log probability: ", stdout); print(max_log_probability, stdout); print("\n", stdout);
		print("Proof of newly-added logical form in the best theory:\n", stdout);
		print<built_in_predicates, identity_canonicalizer, false>(*T.observations.last(), stdout, *debug_terminal_printer);
		print('\n', stdout); fflush(stdout);
		free(T_MAP); free(proof_axioms_MAP);
	}

	print("Finished reading declarative sentences. Attempting to answer question:\n", stdout);
	print("  Logical form: ", stdout); print(*lfs.last(), stdout, printer); print("\n", stdout);

	if (lfs.last()->type != hol_term_type::LAMBDA) {
		typedef typename Theory::Proof Proof;
		Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
		Theory::set_empty(T_MAP);
		Proof* proof_MAP;
		double log_probability_true = log_joint_probability_of_truth(T, proof_prior, proof_axioms, lfs.last(), iterations_per_reheat, num_reheating_phases, reheating_iterations, T_MAP, proof_MAP);
		if (!Theory::is_empty(T_MAP)) {
			print("Highest probability theory after testing whether query is true:\n", stdout);
			T_MAP.print_axioms<true>(stdout, *debug_terminal_printer);
			print("Proof of newly-added logical form in the best theory:\n", stdout);
			print<built_in_predicates, polymorphic_canonicalizer<true, false, built_in_predicates>, false>(*T_MAP.observations.last(), stdout, *debug_terminal_printer);
			print('\n', stdout); fflush(stdout);
			free(T_MAP); Theory::set_empty(T_MAP);
		} else {
			print("Unable to find theory where query is true.\n\n", stdout);
		}

		hol_term* negation;
		if (lfs.last()->type == hol_term_type::NOT) {
			negation = lfs.last()->unary.operand;
			negation->reference_count++;
		} else {
			negation = hol_term::new_not(lfs.last());
			lfs.last()->reference_count++;
		}
		double log_probability_false = log_joint_probability_of_truth(T, proof_prior, proof_axioms, negation, iterations_per_reheat, num_reheating_phases, reheating_iterations, T_MAP, proof_MAP);
		free(*negation); if (negation->reference_count == 0) free(negation);
		if (!Theory::is_empty(T_MAP)) {
			print("Highest probability theory after testing whether query is false:\n", stdout);
			T_MAP.print_axioms<true>(stdout, *debug_terminal_printer);
			print("Proof of newly-added logical form in the best theory:\n", stdout);
			print<built_in_predicates, polymorphic_canonicalizer<true, false, built_in_predicates>, false>(*T_MAP.observations.last(), stdout, *debug_terminal_printer);
			print('\n', stdout); fflush(stdout);
			free(T_MAP);
		} else {
			print("Unable to find theory where query is false.\n\n", stdout);
		}

		double probabilities[] = {log_probability_true, log_probability_false};
		normalize_exp(probabilities, array_length(probabilities));
		print("Answer:\n", stdout);
		print("  true : ", stdout);
		print(probabilities[0], stdout); print('\n', stdout);
		print("  false : ", stdout);
		print(probabilities[1], stdout); print('\n', stdout);
		fflush(stdout);

	} else {
		array_map<string, double> answers(8);
		/* TODO: add support for multiple reheating phases in `answer_question`  */
		answer_question<false>(answers, lfs.last(), iterations_per_reheat * num_reheating_phases, printer, T, proof_prior, proof_axioms);
		sort(answers.values, answers.keys, answers.size, default_sorter());
		//print("Theory after attempting to answer question:\n", stdout);
		//T.template print_axioms<true>(stdout, *debug_terminal_printer); print('\n', stdout).
		reverse(answers.keys, answers.size);
		reverse(answers.values, answers.size);
		normalize_exp(answers.values, answers.size);
		print("Answers: (with estimated probabilities)\n", stdout);
		for (unsigned int i = 0; i < answers.size; i++) {
			print("  ", stdout); print(answers.keys[i], stdout);
			print(" : ", stdout); print(answers.values[i], stdout);
			print('\n', stdout);
		}
		fflush(stdout);
		for (auto entry : answers) core::free(entry.key);
	}

	free(name_map);
	free_all(lfs);
	free_all(seed_axioms);
	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
