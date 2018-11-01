#include "executive.h"
#include "fake_parser.h"

template<typename T>
struct chinese_restaurant_process
{
	double alpha, log_alpha;
	hash_multiset<T> tables;

	chinese_restaurant_process(double alpha) :
		alpha(alpha), log_alpha(log(alpha)), tables(64)
	{ }

	chinese_restaurant_process(const chinese_restaurant_process<T>& src) :
		alpha(src.alpha), log_alpha(log(alpha)), tables(src.tables.counts.table.capacity)
	{
		if (!core::copy(src.tables, tables))
			exit(EXIT_FAILURE);
	}

	~chinese_restaurant_process() {
		core::free(tables);
	}

	inline bool add(const T& customer) {
		return tables.add(customer);
	}

	template<typename BaseDistribution>
	inline bool add(const T& customer, BaseDistribution& base_distribution) {
		if (!tables.counts.check_size()) return false;

		bool contains; unsigned int bucket;
		unsigned int& count = tables.counts.get(customer, contains, bucket);
		if (!contains) {
			if (!base_distribution.add(customer)) return false;
			tables.counts.table.keys[bucket] = customer;
			tables.counts.table.size++;
			count = 1;
		} else {
			count++;
		}
		return true;
	}

	inline void remove(const T& customer) {
		bool contains; unsigned int bucket;
		unsigned int& count = tables.counts.get(customer, contains, bucket);
#if !defined(NDEBUG)
		if (!contains || count == 0) {
			fprintf(stderr, "chinese_restaurant_process.remove WARNING: Given customer is not seated.\n");
			return;
		}
#endif
		if (count == 1)
			tables.counts.remove_at(bucket);
		else count--;
		tables.sum--;
	}

	template<typename BaseDistribution>
	inline void remove(const T& customer, BaseDistribution& base_distribution) {
		bool contains; unsigned int bucket;
		unsigned int& count = tables.counts.get(customer, contains, bucket);
#if !defined(NDEBUG)
		if (!contains || count == 0) {
			fprintf(stderr, "chinese_restaurant_process.remove WARNING: Given customer is not seated.\n");
			return;
		}
#endif
		if (count == 1) {
			tables.counts.remove_at(bucket);
			base_distribution.remove(customer);
		}
		else count--;
		tables.sum--;
	}
};

template<typename T, template<typename> class Collection>
inline double log_probability_ratio_old_cluster(
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters)
{
	unsigned int count = prior.tables.counts.get(observation);
	double value = lgamma(count - frequency) - lgamma(count);
	if (count == frequency) {
		value -= prior.log_alpha;
		old_clusters.add(observation);
	}
	return value;
}

template<typename T, template<typename> class Collection>
inline double log_probability_ratio_new_cluster(
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& new_clusters)
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

template<typename T, template<typename> class Collection>
double log_probability_ratio(
		const array_multiset<T>& old_observations,
		const array_multiset<T>& new_observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters, Collection<T>& new_clusters)
{
	unsigned int i = 0, j = 0; double value = 0.0;
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

	return value + lgamma(prior.alpha + prior.tables.sum)
		 - lgamma(prior.alpha + prior.tables.sum - old_observations.sum + new_observations.sum);
}

template<typename T>
struct dummy_collection {
	constexpr bool add(const T& i) { return true; }
};

template<typename T, template<typename> class Collection>
inline double log_probability_ratio(
		const array_multiset<T>& old_observations,
		const array_multiset<T>& new_observations,
		const chinese_restaurant_process<T>& prior)
{
	dummy_collection<T> old_clusters, new_clusters;
	return log_probability_ratio(old_observations, new_observations, prior, old_clusters, new_clusters);
}

template<typename T, template<typename> class Collection>
inline double log_probability_ratio_old_cluster(
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters)
{
	unsigned int count = prior.tables.counts.get(observation);
	double value = lgamma(count - frequency) - lgamma(count);
	if (count == frequency) {
		value -= prior.log_alpha;
		old_clusters.add(observation);
	}
	return value;
}

template<typename T, template<typename> class Collection>
inline double log_probability_ratio_new_cluster(
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& new_clusters)
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

template<typename T, template<typename> class Collection>
double log_probability_ratio(
		const array_multiset<T>& old_observations,
		const array_multiset<T>& new_observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters, Collection<T>& new_clusters)
{
	unsigned int i = 0, j = 0; double value = 0.0;
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

	return value + lgamma(prior.alpha + prior.tables.sum)
		 - lgamma(prior.alpha + prior.tables.sum - old_observations.sum + new_observations.sum);
}

template<typename BaseDistribution>
struct dirichlet_process
{
	typedef typename BaseDistribution::ObservationType ObservationType;

	chinese_restaurant_process<ObservationType> restaurant;
	BaseDistribution base_distribution;

	dirichlet_process(double alpha, const BaseDistribution& base_distribution) :
		restaurant(alpha), base_distribution(base_distribution)
	{ }

	dirichlet_process(const dirichlet_process<BaseDistribution>& src) :
		restaurant(src.restaurant), base_distribution(src.base_distribution)
	{ }

	inline bool add(const ObservationType& customer) {
		return restaurant.add(customer, base_distribution);
	}

	inline void remove(const ObservationType& customer) {
		restaurant.remove(customer, base_distribution);
	}
};

template<typename BaseDistribution>
inline dirichlet_process<BaseDistribution> make_dirichlet_process(
		double alpha, const BaseDistribution& base_distribution)
{
	return dirichlet_process<BaseDistribution>(alpha, base_distribution);
}

template<typename BaseDistribution, typename T,
		typename BaseDistributionOldParameter,
		typename BaseDistributionNewParameter>
double log_probability_ratio(
		const array_multiset<T>& old_observations,
		const array_multiset<T>& new_observations,
		const dirichlet_process<BaseDistribution>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters old_clusters, new_clusters;
	return log_probability_ratio(old_observations, new_observations, prior.restaurant, old_clusters, new_clusters)
		 + log_probability_ratio(prior.restaurant.tables.counts.table, prior.base_distribution, old_clusters, new_clusters);
}

template<typename ConstantDistribution, typename PredicateDistribution>
struct simple_constant_distribution
{
	typedef constant_array ObservationCollection;

	struct constant_array {
		array<unsigned int> constants;
		array<unsigned int> predicates;

		constant_array() : constants(8), predicates(8) { }

		inline bool add_constant(unsigned int constant) {
			return constants.add(constant);
		}

		inline bool add_predicate(unsigned int predicate) {
			return predicates.add(predicate);
		}
	};

	ConstantDistribution constant_distribution;
	PredicateDistribution predicate_distribution;

	simple_constant_distribution(
			const ConstantDistribution& constant_distribution,
			const PredicateDistribution& predicate_distribution) :
		constant_distribution(constant_distribution),
		predicate_distribution(predicate_distribution)
	{ }

	simple_constant_distribution(const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& src) :
		constant_distribution(src.constant_distribution), predicate_distribution(src.predicate_distribution)
	{ }

	bool add_constant(unsigned int constant) {
		return constant_distribution.add(constant);
	}

	bool add_predicate(unsigned int predicate) {
		return predicate_distribution.add(predicate);
	}

	void remove_constant(unsigned int constant) {
		constant_distribution.remove(constant);
	}

	void remove_predicate(unsigned int predicate) {
		predicate_distribution.remove(predicate);
	}
};

template<typename ConstantDistribution, typename PredicateDistribution>
inline simple_constant_distribution<ConstantDistribution, PredicateDistribution>
make_simple_constant_distribution(
		const ConstantDistribution& constant_distribution,
		const PredicateDistribution& predicate_distribution)
{
	return simple_constant_distribution<ConstantDistribution, PredicateDistribution>(constant_distribution, predicate_distribution);
}

template<typename ConstantCollection, typename ConstantDistribution, typename PredicateDistribution>
inline double log_probability_ratio(
		const ConstantCollection& old_constants, const ConstantCollection& new_constants,
		const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& prior)
{
	return log_probability_ratio(old_constants.constants, new_constants.constants, prior.constant_distribution)
		 + log_probability_ratio(old_constants.predicates, new_constants.predicates, prior.predicate_distribution);
}

template<typename T>
struct default_array {
	array<T> a;

	default_array() : a(16) { }

	inline bool add(const T& item) {
		return a.add(item);
	}

	inline auto begin() -> decltype(a.begin()) {
		return a.begin();
	}

	inline auto end() -> decltype(a.end()) {
		return a.end();
	}
};

template<typename T>
inline unsigned int size(const default_array<T>& array) {
	return array.a.length;
}

template<typename ConstantDistribution>
struct constant_adder {
	ConstantDistribution& constant_distribution;
};

template<typename ConstantDistribution>
inline bool visit_constant(unsigned int constant, const constant_adder<ConstantDistribution>& visitor) {
	return visitor.constant_distribution.add_constant(constant);
}

template<typename ConstantDistribution>
inline bool visit_predicate(unsigned int predicate, const constant_adder<ConstantDistribution>& visitor) {
	return visitor.constant_distribution.add_predicate(constant);
}

template<typename ConstantDistribution> constexpr bool visit_variable(unsigned int variable, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_parameter(unsigned int parameter, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_true(const fol_formula& formula, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_false(const fol_formula& formula, const constant_adder<ConstantDistribution>& visitor) { return true; }

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula& formula, const parameter_collector& visitor) { return true; }

template<typename ConstantDistribution>
struct constant_remover {
	ConstantDistribution& constant_distribution;
};

template<typename ConstantDistribution>
inline bool visit_constant(unsigned int constant, const constant_remover<ConstantDistribution>& visitor) {
	visitor.constant_distribution.remove_constant(constant);
	return true;
}

template<typename ConstantDistribution>
inline bool visit_predicate(unsigned int predicate, const constant_remover<ConstantDistribution>& visitor) {
	visitor.constant_distribution.remove_predicate(constant);
	return true;
}

template<typename ConstantDistribution> constexpr bool visit_variable(unsigned int variable, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_parameter(unsigned int parameter, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_true(const fol_formula& formula, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_false(const fol_formula& formula, const constant_remover<ConstantDistribution>& visitor) { return true; }

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula& formula, const parameter_collector& visitor) { return true; }

template<typename ConstantDistribution>
struct simple_fol_formula_distribution
{
	typedef fol_formula* ObservationType;
	typedef default_array<fol_formula*> ObservationCollection;

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

	ConstantDistribution constant_distribution;

	simple_fol_formula_distribution(const ConstantDistribution& constant_distribution,
			double ground_literal_probability, double negation_probability,
			double unary_probability, double antecedent_stop_probability,
			double consequent_stop_probability) :
		log_ground_literal_probability(log(ground_literal_probability)),
		log_universal_probability(log(1.0 - ground_literal_probability)),
		log_negation_probability(log(negation_probability)),
		log_positive_probability(log(1.0 - negation_probability)),
		log_unary_probability(log(unary_probability)),
		log_binary_probability(log(1.0 - unary_probability)),
		log_antecedent_continue_probability(log(1.0 - antecedent_stop_probability)),
		log_antecedent_stop_probability(log(antecedent_stop_probability)),
		log_consequent_continue_probability(log(1.0 - consequent_stop_probability)),
		log_consequent_stop_probability(log(consequent_stop_probability)),
		constant_distribution(constant_distribution)
	{ }

	simple_fol_formula_distribution(const simple_fol_formula_distribution<ConstantDistribution>& src) :
		log_ground_literal_probability(src.log_ground_literal_probability),
		log_universal_probability(src.log_universal_probability),
		log_negation_probability(src.log_negation_probability),
		log_positive_probability(src.log_positive_probability),
		log_unary_probability(src.log_unary_probability),
		log_binary_probability(src.log_binary_probability),
		log_antecedent_continue_probability(src.log_antecedent_continue_probability),
		log_antecedent_stop_probability(src.log_antecedent_stop_probability),
		log_consequent_continue_probability(src.log_consequent_continue_probability),
		log_consequent_stop_probability(src.log_consequent_stop_probability),
		constant_distribution(src.constant_distribution)
	{ }

	inline bool add(fol_formula* formula) {
		constant_adder<ConstantDistribution> adder = { constant_distribution };
		return visit(*formula, adder);
	}

	inline bool remove(fol_formula* formula) {
		constant_remover<ConstantDistribution> remover = { constant_distribution };
		return visit(*formula, remover);
	}
};

template<typename ConstantDistribution>
inline simple_fol_formula_distribution<ConstantDistribution> make_simple_fol_formula_distribution(
		const ConstantDistribution& constant_distribution,
		double ground_literal_probability, double negation_probability,
		double unary_probability, double antecedent_stop_probability,
		double consequent_stop_probability)
{
	return simple_fol_formula_distribution<ConstantDistribution>(
			constant_distribution, ground_literal_probability,
			negation_probability, unary_probability,
			antecedent_stop_probability, consequent_stop_probability);
}

template<typename ConstantDistribution>
double log_probability_atom(const fol_atom& atom,
		const simple_fol_formula_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	constants.add_predicate(atom.predicate);
	if (atom.arg2.type == fol_term_type::NONE) {
		constants.add_constant(atom.arg2.constant);
		return prior.log_unary_probability;
	} else {
		constants.add_constant(atom.arg1.constant);
		constants.add_constant(atom.arg2.constant);
		return prior.log_binary_probability;
	}
}

template<typename ConstantDistribution>
double log_probability_literal(const fol_formula* literal,
		const simple_fol_formula_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (literal->type == fol_formula_type::ATOM) {
		return prior.log_positive_probability + log_probability_atom(literal->atom, prior, constants);
	} else if (literal->type == fol_formula_type::NOT && literal->unary.operand->type == fol_formula_type::ATOM) {
		return prior.log_negation_probability + log_probability_atom(literal->unary.operand->atom, prior, constants);
	} else {
		return std::numeric_limits<double>::infinity();
	}
}

template<typename ConstantDistribution>
double log_probability_helper(const fol_formula* formula,
		const simple_fol_formula_distribution<ConstantDistribution>& prior,
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
		return prior.log_ground_literal_probability + log_probability_literal(formula, prior, constants);
	case fol_formula_type::FOR_ALL:
		if (formula->quantifier.operand->type == fol_formula_type::IF_THEN) {
			antecedent = formula->quantifier.operand->binary.left;
			consequent = formula->quantifier.operand->binary.right;
			if (antecedent->type == fol_formula_type::AND) {
				value = (antecedent->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < antecedent->array.length; i++)
					value += log_probability_literal(antecedent->array.operands[i], prior, constants);
			} else {
				value = prior.log_antecedent_stop_probability + log_probability_literal(antecedent, prior, constants);
			}

			if (consequent->type == fol_formula_type::AND) {
				value = (consequent->array.length - 1) * prior.log_consequent_continue_probability + prior.log_consequent_stop_probability;
				for (unsigned int i = 0; i < consequent->array.length; i++)
					value += log_probability_literal(consequent->array.operands[i], prior, constants);
			} else {
				value = prior.log_consequent_stop_probability + log_probability_literal(consequent, prior, constants);
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

template<typename ConstantDistribution,
	template<typename> class ObservationCollection,
	template<typename> class Collection>
double log_probability_ratio(
		const Collection<fol_formula*>& observations,
		const simple_fol_formula_distribution<ConstantDistribution>& prior,
		const Collection<fol_formula*>& old_clusters,
		const Collection<fol_formula*>& new_clusters)
{
	double value = 0.0;
	typename ConstantDistribution::ObservationCollection old_constants, new_constants;
	for (fol_formula* formula : new_clusters)
		value += log_probability_helper(formula, prior, new_constants);
	for (fol_formula* formula : old_clusters)
		value -= log_probability_helper(formula, prior, old_constants);
	return value + log_probability_ratio(old_constants, new_constants, prior.constant_distribution);
}

template<typename T>
struct uniform_subset_distribution {
	double log_stop_probability;
	double log_continue_probability;

	typedef default_array<T> ObservationCollection;

	uniform_subset_distribution(double stop_probability) :
		log_stop_probability(log(stop_probability)),
		log_continue_probability(log(1.0 - stop_probability))
	{ }

	uniform_subset_distribution(const uniform_subset_distribution<T>& src) :
		log_stop_probability(src.log_stop_probability),
		log_continue_probability(src.log_continue_probability)
	{ }
};

template<template<typename> class Collection, typename T>
double log_probability(
		const pair<Collection<T>, Collection<T>>& observation,
		const uniform_subset_distribution<T>& prior)
{
	log_cache<double>::instance().ensure_size(size(observation.value));
	return -log_cache<double>::instance().get(size(observation.value)) * size(observation.key)
		 + size(observation) * prior.log_continue_probability + prior.log_stop_probability;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
double log_probability(
		const ObservationCollection<pair<Collection<T>, Collection<T>>>& observations,
		const uniform_subset_distribution<T>& prior)
{
	double value = 0.0;
	for (const pair<Collection<T>, Collection<T>>& observation : observations)
		value += log_probability(observation, prior);
	return value;
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
	theory<fol_formula, natural_deduction<fol_formula, true>> T;
	auto constant_prior = make_simple_constant_distribution(
			chinese_restaurant_process<unsigned int>(1.0), chinese_restaurant_process<unsigned int>(1.0));
	auto theory_element_prior = make_simple_fol_formula_distribution(constant_prior, 0.5, 0.3, 0.4, 0.2, 0.4);
	auto axiom_prior = make_dirichlet_process(1.0, theory_element_prior);
	auto conjunction_prior = uniform_subset_distribution<nd_step<fol_formula, true>*>(0.5);
	auto universal_introduction_prior = uniform_subset_distribution<unsigned int >(0.5);
	auto universal_elimination_prior = chinese_restaurant_process<fol_term>(1.0);
	auto proof_prior = make_canonicalized_proof_prior(0.1, axiom_prior,
			conjunction_prior, universal_introduction_prior, universal_elimination_prior);
	read_article(names.get("bob"), corpus, parser, T, proof_prior);

	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
