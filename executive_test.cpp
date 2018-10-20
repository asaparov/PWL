#include "executive.h"
#include "fake_parser.h"

struct chinese_restaurant_process
{
	double alpha, log_alpha;
	array_multiset<unsigned int> tables;
};

template<typename PredicateDistribution, typename ConstantDistribution>
struct simple_fol_formula_distribution
{
	double log_ground_literal_probability;
	double log_universal_probability;

	double log_negation_probability;
	double log_positive_probability;

	double log_unary_probability;
	double log_binary_probability;

	double log_antecedent_continue_probability;
	double log_antecedent_stop_probability;
	double log_consequent_continue_probability;
	double log_consequent_stop_probability;

	PredicateDistribution predicate_distribution;
	ConstantDistribution constant_distribution;
};

template<typename PredicateDistribution, typename ConstantDistribution>
double log_probability_atom(const fol_atom& atom,
		const simple_fol_formula_distribution<PredicateDistribution, ConstantDistribution>& prior,
		typename PredicateDistribution::ObservationCollection& predicates,
		typename ConstantDistribution::ObservationCollection& constants)
{
	predicates.add(formula->atom.predicate);
	if (formula->atom.arg2.type == fol_term_type::NONE) {
		constants.add(formula->atom.arg2.constant);
		return prior.log_unary_probability;
	} else {
		constants.add(formula->atom.arg1.constant);
		constants.add(formula->atom.arg2.constant);
		return prior.log_binary_probability;
	}
}

template<typename PredicateDistribution, typename ConstantDistribution>
double log_probability_literal(const fol_formula* literal,
		const simple_fol_formula_distribution<PredicateDistribution, ConstantDistribution>& prior,
		typename PredicateDistribution::ObservationCollection& predicates,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (literal->type == fol_formula_type::ATOM) {
		return prior.log_positive_probability + log_probability_atom(literal->atom, prior, predicates, constants);
	} else if (literal->type == fol_formula_type::NOT && literal->unary.operand->type == fol_formula_type::ATOM) {
		return prior.log_negation_probability + log_probability_atom(literal->unary.operand->atom, prior, predicates, constants);
	} else {
		return std::numeric_limits<double>::infinity();
	}
}

template<typename PredicateDistribution, typename ConstantDistribution>
double log_probability(const fol_formula* formula,
		const simple_fol_formula_distribution<PredicateDistribution, ConstantDistribution>& prior,
		typename PredicateDistribution::ObservationCollection& predicates,
		typename ConstantDistribution::ObservationCollection& constants)
{
	double value;
	const fol_formula* antecedent;
	const fol_formula* consequent;
	const fol_formula* operand;
	switch (formula->type)
	{
	case fol_formula_type::ATOM:
	case fol_formula_type::NOT:
		return prior.log_ground_literal_probability + log_probability_literal(formula, prior, predicates, constants);
	case fol_formula_type::FOR_ALL:
		if (formula->quantifier.operand->type == fol_formula_type::IF_THEN) {
			antecedent = formula->quantifier.operand->binary.left;
			consequent = formula->quantifier.operand->binary.right;
			if (antecedent->type == fol_formula_type::AND) {
				value = (antecedent->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < antecedent->array.length; i++)
					value += log_probability_literal(antecedent->array.operands[i], prior, predicates, constants);
			} else {
				value = prior.log_antecedent_stop_probability + log_probability_literal(antecedent, prior, predicates, constants);
			}

			if (consequent->type == fol_formula_type::AND) {
				value = (consequent->array.length - 1) * prior.log_consequent_continue_probability + prior.log_consequent_stop_probability;
				for (unsigned int i = 0; i < consequent->array.length; i++)
					value += log_probability_literal(consequent->array.operands[i], prior, predicates, constants);
			} else {
				value = prior.log_consequent_stop_probability + log_probability_literal(consequent, prior, predicates, constants);
			}
			return value;
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	case fol_formula_type::AND:
	case fol_formula_type::OR:
	case fol_formula_type::IF_THEN:
	case fol_formula_type::IFF:
	case fol_formula_type::EXISTS:
	case fol_formula_type::TRUE:
	case fol_formula_type::FALSE:
		return -std::numeric_limits<double>::infinity();
	}
	fprintf(stderr, "log_probability ERROR: Unrecognized fol_formula_type.\n");
	exit(EXIT_FAILURE);
}

template<template<typename> class BaseDistribution, typename T>
struct identically_sampled_distribution {
	BaseDistribution<T> base_distribution;

	typedef array<T> ObservationCollection;
};

template<template<typename> class Collection,
	template<typename> class BaseDistribution,
	template<typename> class ObservationCollection, typename T,
	typename NewBaseDistributionParameters,
	typename OldBaseDistributionParameters>
double log_probability_ratio(const Collection<T>& observations,
		const identically_sampled_distribution<BaseDistribution, T>& prior,
		const ObservationCollection<T>& old_observations,
		const ObservationCollection<T>& new_observations,
		const NewBaseDistributionParameters& new_base_distribution_parameters,
		const OldBaseDistributionParameters& old_base_distribution_parameters)
{
	double value = 0.0;
	for (const T& observation : new_observations)
		value += log_probability(observation, prior.base_distribution, new_base_distribution_parameters);
	for (const T& observation : old_observations)
		value -= log_probability(observation, prior.base_distribution, old_base_distribution_parameters);
	return value;
}

template<template<typename> class BaseDistribution, typename T>
struct dirichlet_process
{
	double alpha, log_alpha;
	BaseDistribution<T> base_distribution;
	hash_multiset<T> tables;

	dirichlet_process(double alpha, const BaseDistribution<T>& base_distribution) :
			alpha(alpha), log_alpha(log(alpha)), base_distribution(base_distribution), tables(64)
	{ }

	inline bool add(const T& customer) {
		return tables.add(customer);
	}

	inline void remove(const T& customer) {
		bool contains; unsigned int bucket;
		unsigned int& count = tables.counts.get(customer, contains, bucket);
#if !defined(NDEBUG)
		if (!contains || count == 0) {
			fprintf(stderr, "dirichlet_process.remove WARNING: Given customer is not seated.\n");
			return;
		}
#endif
		if (count == 1)
			tables.counts.remove_at(bucket);
		else count--;
		tables.sum--;
	}
};

template<template<typename> class BaseDistribution, typename T>
inline double log_probability_ratio_old_cluster(
		const T& observation, unsigned int frequency,
		const dirichlet_process<BaseDistribution, T>& prior,
		typename BaseDistribution<T>::ObservationCollection& old_clusters)
{
	unsigned int count = prior.tables.counts.get(observation);
	double value = lgamma(count - frequency) - lgamma(count);
	if (count == frequency) {
		value -= prior.log_alpha;
		old_clusters.add(observation);
	}
	return value;
}

template<template<typename> class BaseDistribution, typename T>
inline double log_probability_ratio_new_cluster(
		const T& observation, unsigned int frequency,
		const dirichlet_process<BaseDistribution, T>& prior,
		typename BaseDistribution<T>::ObservationCollection& new_clusters)
{
	bool contains;
	unsigned int count = prior.tables.counts.get(observation, contains);
	if (contains) {
		return lgamma(count + frequency) - lgamma(count);
	} else {
		new_clusters.add(observation);
		return prior.log_alpha + lgamma(frequency);
	}
}

template<template<typename> class BaseDistribution, typename T,
		typename BaseDistributionOldParameter,
		typename BaseDistributionNewParameter>
double log_probability_ratio(
		const array_multiset<T>& old_observations,
		const array_multiset<T>& new_observations,
		const dirichlet_process<BaseDistribution, T>& prior,
		BaseDistributionOldParameter& base_distribution_old_parameter,
		BaseDistributionNewParameter& base_distribution_new_parameter)
{
	typedef typename BaseDistribution<T>::ObservationCollection Clusters;

	unsigned int i = 0, j = 0; double value = 0.0;
	Clusters old_clusters, new_clusters;
	while (i < old_observations.counts.size && j < new_observations.counts.size)
	{
		if (old_observations.counts.keys[i] == new_observations.counts.keys[j]) {
			unsigned int count = prior.tables.counts.get(old_observations.counts.keys[i]);
			int diff = (int) new_observations.counts.values[j] - (int) old_observations.counts.values[i];
			value += lgamma(count + diff) - lgamma(count);
			i++; j++;
		} else if (old_observations.counts.keys[i] < new_observations.counts.keys[j]) {
			value += log_probability_ratio_old_cluster(
					old_observations.counts.keys[i], old_observations.counts.values[i], prior, old_clusters);
			i++;
		} else {
			value += log_probability_ratio_new_cluster(
					new_observations.counts.keys[j], new_observations.counts.values[j], prior, new_clusters);
			j++;
		}
	}

	while (i < old_observations.counts.size) {
		value += log_probability_ratio_old_cluster(
				old_observations.counts.keys[i], old_observations.counts.values[i], prior, old_clusters);
		i++;
	} while (j < new_observations.counts.size) {
		value += log_probability_ratio_new_cluster(
				new_observations.counts.keys[j], new_observations.counts.values[j], prior, new_clusters);
		j++;
	}

	value += lgamma(prior.alpha + prior.tables.sum) - lgamma(prior.alpha + prior.tables.sum - old_observations.sum + new_observations.sum);
	return value + log_probability_ratio(prior.tables.counts.table, prior.base_distribution,
			old_clusters, new_clusters, base_distribution_old_parameter, base_distribution_new_parameter);
}

const string* get_name(const hash_map<string, unsigned int>& names, unsigned int id)
{
	for (const auto& entry : names)
		if (entry.value == id) return &entry.key;
	return NULL;
}

int main(int argc, const char** argv)
{
	FILE* in = fopen("simple_cat_articles.txt", "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open file for reading.\n");
		return EXIT_FAILURE;
	}

	array<article_token> tokens = array<article_token>(256);
	if (!article_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis of article failed.\n");
		fclose(in); free_tokens(tokens); return EXIT_FAILURE;
	}
	fclose(in);

	unsigned int index = 0;
	in_memory_article_store corpus;
	fake_parser parser = fake_parser(PREDICATE_UNKNOWN);
	hash_map<string, unsigned int> names(256);
	while (index < tokens.length) {
		unsigned int article_name = 0;
		article& new_article = *((article*) alloca(sizeof(article)));
		if (!corpus.articles.check_size() || !article_interpret(tokens, index, new_article, article_name, parser.table, names)) {
			const string* article_name_str = get_name(names, article_name);
			if (article_name_str == NULL) {
				fprintf(stderr, "ERROR: Unable to parse article %u.\n", corpus.articles.table.size + 1);
			} else {
				print("ERROR: Unable to parse article with title '", stderr); print(*article_name_str, stderr); print("'.\n", stderr);
			}
			for (auto entry : names) free(entry.key);
			free_tokens(tokens); return EXIT_FAILURE;
		}

		bool contains; unsigned int bucket;
		article& value = corpus.articles.get(article_name, contains, bucket);
		if (contains) {
			const string* article_name_str = get_name(names, article_name);
			print("ERROR: Article with title '", stderr); print(*article_name_str, stderr); print("' already exists.\n", stderr);
			for (auto entry : names) free(entry.key);
			free_tokens(tokens); free(new_article); return EXIT_FAILURE;
		}
		move(new_article, value);
		corpus.articles.table.keys[bucket] = article_name;
		corpus.articles.table.size++;
	}
	free_tokens(tokens);

	/* read the articles */
	double proof_stop_probability = 0.2;
	double log_proof_stop_probability = log(proof_stop_probability);
	double log_proof_continue_probability = log(1.0 - proof_stop_probability);
	theory<fol_formula, natural_deduction<fol_formula, true>> T;
	/*read_article(names.get("bob"), corpus, parser, T,
		log_proof_stop_probability, log_proof_continue_probability,
		TheoryPrior& theory_prior, AxiomPrior& axiom_prior,
		UniversalIntroductionPrior& universal_introduction_prior,
		UniversalEliminationPrior& universal_elimination_prior);*/

	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
