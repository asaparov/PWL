#include "executive.h"
#include "fake_parser.h"

template<typename ProofCalculus>
double log_probability_ratio(theory<fol_formula, ProofCalculus>& T) {

}

template<typename T>
struct uniform_distribution { };

template<template<typename, typename> class MapType, typename T>
inline double log_probability(unsigned int size,
		const uniform_distribution<T>& prior)
{
	if (!log_cache<double>::instance().ensure_size(size)) exit(EXIT_FAILURE);
	return -log_cache<double>::instance().get(size);
}

template<template<typename, typename> class MapType, typename T>
double log_probability(
		const MapType<T, unsigned int>& multiset,
		const uniform_distribution<T>& prior)
{
	return log_probability(multiset.size, prior);
}

template<template<typename, typename> class MapType, typename T>
double log_probability(
		const MapType<T, unsigned int>& multiset,
		const array<T>& new_clusters,
		const uniform_distribution<T>& prior)
{
	return log_probability(multiset.size + new_clusters.length, prior);
}

template<typename Formula>
struct uniform_axiom_distribution {
	struct counter {
		unsigned int count;

		counter() : count(0) { }

		inline bool add(const Formula* observation) {
			count++; return true;
		}

		inline void remove(const Formula* observation) {
			count--;
		}
	};

	typedef counter ObservationCollection;
};

template<template<typename> class ContainerType, typename Formula, typename ProofCalculus, bool Negated, unsigned int Arity>
double log_probability_ratio(const ContainerType<Formula*>& axioms,
		const uniform_axiom_distribution<Formula>& prior,
		const typename uniform_axiom_distribution<Formula>::counter& new_axioms,
		const typename uniform_axiom_distribution<Formula>::counter& old_axioms)
{
	
}

template<typename BaseDistribution, typename T>
struct dirichlet_process
{
	double alpha, log_alpha;
	BaseDistribution base_distribution;
	hash_multiset<T> tables;

	dirichlet_process(double alpha, const BaseDistribution& base_distribution) :
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

template<typename BaseDistribution, typename T>
inline double log_probability_ratio_old_cluster(
		const T& observation, unsigned int frequency,
		const dirichlet_process<BaseDistribution, T>& prior,
		typename BaseDistribution::ObservationCollection& old_clusters)
{
	unsigned int count = prior.tables.counts.get(observation);
	double value = lgamma(count - frequency) - lgamma(count);
	if (count == frequency) {
		value -= prior.log_alpha;
		old_clusters.add(observation);
	}
	return value;
}

template<typename BaseDistribution, typename T>
inline double log_probability_ratio_new_cluster(
		const T& observation, unsigned int frequency,
		const dirichlet_process<BaseDistribution, T>& prior,
		typename BaseDistribution::ObservationCollection& new_clusters)
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

template<typename BaseDistribution, typename T>
double log_probability_ratio(
		const array_multiset<T>& new_observations,
		const array_multiset<T>& old_observations,
		const dirichlet_process<BaseDistribution, T>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

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
	return value + log_probability_ratio(prior.tables.counts.table, new_clusters, old_clusters, prior.base_distribution);
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
