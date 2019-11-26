#include "higher_order_logic.h"
#include "executive.h"
#include "hdp_parser.h"

#include <boost/math/special_functions/gamma.hpp>

struct poisson_distribution {
	typedef unsigned int ObservationType;

	double lambda, log_lambda;

	poisson_distribution(double lambda) : lambda(lambda), log_lambda(log(lambda)) { }
	poisson_distribution(const poisson_distribution& src) : lambda(src.lambda), log_lambda(src.log_lambda) { }
};

inline double log_probability(unsigned int k, const poisson_distribution& prior) {
	return k * prior.log_lambda - prior.lambda - lgamma(k + 1);
}

struct geometric_distribution {
	typedef unsigned int ObservationType;

	double p, log_p, log_one_minus_p;

	geometric_distribution(double p) : p(p), log_p(log(p)), log_one_minus_p(log(1.0 - p)) { }
	geometric_distribution(const geometric_distribution& src) : p(src.p), log_p(src.log_p), log_one_minus_p(src.log_one_minus_p) { }
};

inline double log_probability(unsigned int k, const geometric_distribution& prior) {
	return k * prior.log_one_minus_p + prior.log_p;
}

unsigned int sample(const geometric_distribution& prior, unsigned int min, unsigned int max) {
	/* we use inverse transform sampling */
	double u = sample_uniform<double>();
	if (max == UINT_MAX)
		return floor((min * prior.log_one_minus_p + log(1 - u)) / prior.log_one_minus_p);
	return floor(log(pow(1 - prior.p, min) * (1 - u) + pow(1 - prior.p, max + 1) * u) / prior.log_one_minus_p);
}

template<typename T>
struct default_array {
	array<T> a;

	default_array() : a(16) { }

	inline bool add(const T& item) {
		return a.add(item);
	}

	inline auto begin() const -> decltype(a.begin()) {
		return a.begin();
	}

	inline auto end() const -> decltype(a.end()) {
		return a.end();
	}
};

template<typename T>
inline unsigned int size(const default_array<T>& array) {
	return array.a.length;
}

template<typename T>
struct default_array_multiset {
	array_multiset<T, false> a;

	default_array_multiset() : a(16) { }

	inline bool add(const T& item) {
		return a.add(item);
	}
};

template<typename T>
inline unsigned int size(const default_array_multiset<T>& multiset) {
	return multiset.a.counts.size;
}

template<typename T>
inline const T& get_key(const default_array_multiset<T>& multiset, unsigned int i) {
	return multiset.a.counts.keys[i];
}

template<typename T>
inline unsigned int get_value(const default_array_multiset<T>& multiset, unsigned int i) {
	return multiset.a.counts.values[i];
}

template<typename T>
unsigned int sum(const default_array_multiset<T>& multiset) {
	return multiset.a.sum;
}

template<typename T, bool AutomaticallyFree>
inline unsigned int size(const array_multiset<T, AutomaticallyFree>& multiset) {
	return multiset.counts.size;
}

template<typename T, bool AutomaticallyFree>
inline const T& get_key(const array_multiset<T, AutomaticallyFree>& multiset, unsigned int i) {
	return multiset.counts.keys[i];
}

template<typename T, bool AutomaticallyFree>
inline unsigned int get_value(const array_multiset<T, AutomaticallyFree>& multiset, unsigned int i) {
	return multiset.counts.values[i];
}

template<typename T, bool AutomaticallyFree>
unsigned int sum(const array_multiset<T, AutomaticallyFree>& multiset) {
	return multiset.sum;
}

template<typename T>
struct chinese_restaurant_process
{
	typedef default_array_multiset<T> ObservationCollection;

	double alpha, log_alpha;

	chinese_restaurant_process(double alpha) :
		alpha(alpha), log_alpha(log(alpha))
	{ }

	chinese_restaurant_process(const chinese_restaurant_process<T>& src) :
		alpha(src.alpha), log_alpha(log(alpha))
	{ }

	~chinese_restaurant_process() { }
};

template<typename MultisetType,
	template<typename> class Collection, typename T>
double log_probability(
		const MultisetType& observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& clusters)
{
	double value = 0.0;
	for (unsigned int i = 0; i < size(observations); i++) {
		clusters.add(get_key(observations, i));
		value += prior.log_alpha + lgamma(get_value(observations, i));
	}
	return value + lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
}

template<typename MultisetType,
	template<typename> class Collection, typename T>
inline double log_probability_ratio(
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters, Collection<T>& new_clusters)
{
	return log_probability(new_observations, prior, new_clusters) - log_probability(old_observations, prior, old_clusters);
}

template<typename T>
struct dummy_collection {
	constexpr bool add(const T& i) const { return true; }
};

template<typename MultisetType, typename T>
inline double log_probability_ratio(
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior)
{
	dummy_collection<T> old_clusters, new_clusters;
	return log_probability_ratio(old_observations, new_observations, prior, old_clusters, new_clusters);
}

template<bool AutomaticallyFree, typename T, template<typename> class Collection>
inline double log_probability_ratio_old_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters)
{
#if !defined(NDEBUG)
	if (!tables.counts.table.contains(observation))
		fprintf(stderr, "log_probability_ratio_old_cluster WARNING: This chinese_restaurant_process does not contain the given observation.\n");
#endif
	unsigned int count = tables.counts.get(observation);
	if (count == frequency) {
		old_clusters.add(observation);
		return -lgamma(count) - prior.log_alpha;
	} else {
		return lgamma(count - frequency) - lgamma(count);
	}
}

template<bool AutomaticallyFree, typename T, template<typename> class Collection>
inline double log_probability_ratio_new_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& new_clusters)
{
	bool contains;
	unsigned int count = tables.counts.get(observation, contains);
#if !defined(NDEBUG)
	if (contains && count == 0)
		fprintf(stderr, "log_probability_ratio_new_cluster WARNING: The hash_multiset has an observation with zero count.\n");
#endif
	if (contains) {
		return lgamma(count + frequency) - lgamma(count);
	} else {
		new_clusters.add(observation);
		return prior.log_alpha + lgamma(frequency);
	}
}

template<bool AutomaticallyFree, typename MultisetType,
	template<typename> class Collection, typename T>
double log_probability_ratio(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters, Collection<T>& new_clusters)
{
	unsigned int i = 0, j = 0; double value = 0.0;
	while (i < size(old_observations) && j < size(new_observations))
	{
		if (get_key(old_observations, i) == get_key(new_observations, j)) {
#if !defined(NDEBUG)
			bool contains;
			const T& key = get_key(old_observations, i);
			unsigned int count = tables.counts.get(key, contains);
			if (!contains)
				print("log_probability_ratio WARNING: The given observation is not in the chinese_restaurant_process.\n", stderr);
#else
			unsigned int count = tables.counts.get(get_key(old_observations, i));
#endif
			int diff = (int) get_value(new_observations, j) - (int) get_value(old_observations, i);
			if (count + diff == 0) {
				/* this cluster is being removed */
				old_clusters.add(get_key(old_observations, i));
				value += -lgamma(count) - prior.log_alpha;
			} else {
				value += lgamma(count + diff) - lgamma(count);
			}
			i++; j++;
		} else if (get_key(old_observations, i) < get_key(new_observations, j)) {
			value += log_probability_ratio_old_cluster(tables,
					get_key(old_observations, i), get_value(old_observations, i), prior, old_clusters);
			i++;
		} else {
			value += log_probability_ratio_new_cluster(tables,
					get_key(new_observations, j), get_value(new_observations, j), prior, new_clusters);
			j++;
		}
	}

	while (i < size(old_observations)) {
		value += log_probability_ratio_old_cluster(tables,
				get_key(old_observations, i), get_value(old_observations, i), prior, old_clusters);
		i++;
	} while (j < size(new_observations)) {
		value += log_probability_ratio_new_cluster(tables,
				get_key(new_observations, j), get_value(new_observations, j), prior, new_clusters);
		j++;
	}

	return value + lgamma(prior.alpha + tables.sum)
		 - lgamma(prior.alpha + tables.sum - sum(old_observations) + sum(new_observations));
}

template<typename MultisetType, typename T>
inline double log_probability_ratio(
		const hash_multiset<T>& tables,
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior)
{
	dummy_collection<T> old_clusters, new_clusters;
	return log_probability_ratio(tables, old_observations, new_observations, prior, old_clusters, new_clusters);
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
};

template<typename BaseDistribution>
inline dirichlet_process<BaseDistribution> make_dirichlet_process(
		double alpha, const BaseDistribution& base_distribution)
{
	return dirichlet_process<BaseDistribution>(alpha, base_distribution);
}

template<typename MultisetType, typename T, bool AutomaticallyFree, typename BaseDistribution>
double log_probability_ratio(
		const MultisetType& restaurant_tables,
		const array_multiset<T, AutomaticallyFree>& old_observations,
		const array_multiset<T, AutomaticallyFree>& new_observations,
		const dirichlet_process<BaseDistribution>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters old_clusters, new_clusters;
	return log_probability_ratio(restaurant_tables, old_observations, new_observations, prior.restaurant, old_clusters, new_clusters)
		 + log_probability_ratio(restaurant_tables.counts.table, prior.base_distribution, old_clusters, new_clusters);
}

template<typename ConstantDistribution, typename PredicateDistribution>
struct simple_constant_distribution
{
	struct constant_array {
		array_multiset<unsigned int, false> constants;
		array_multiset<unsigned int, false> predicates;

		constant_array() : constants(8), predicates(8) { }

		inline bool add_constant(unsigned int constant) {
			return constants.add(constant);
		}

		inline bool add_predicate(unsigned int predicate) {
			return predicates.add(predicate);
		}
	};

	typedef constant_array ObservationCollection;

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

template<typename ConstantDistribution>
struct simple_hol_term_distribution
{
	typedef hol_term* ObservationType;
	typedef default_array<hol_term*> ObservationCollection;

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

	simple_hol_term_distribution(const ConstantDistribution& constant_distribution,
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

	simple_hol_term_distribution(const simple_hol_term_distribution<ConstantDistribution>& src) :
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
};

template<typename ConstantDistribution>
inline simple_hol_term_distribution<ConstantDistribution> make_simple_hol_term_distribution(
		const ConstantDistribution& constant_distribution,
		double ground_literal_probability, double negation_probability,
		double unary_probability, double antecedent_stop_probability,
		double consequent_stop_probability)
{
	return simple_hol_term_distribution<ConstantDistribution>(
			constant_distribution, ground_literal_probability,
			negation_probability, unary_probability,
			antecedent_stop_probability, consequent_stop_probability);
}

template<bool Quantified, typename ConstantDistribution>
double log_probability_atom(const hol_term* function, const hol_term* arg1,
		const simple_hol_term_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type != hol_term_type::CONSTANT)
		return std::numeric_limits<double>::infinity();

	constants.add_predicate(function->constant);
	if (Quantified && arg1->type == hol_term_type::VARIABLE) {
		return prior.log_unary_probability;
	} else if (!Quantified && arg1->type == hol_term_type::CONSTANT) {
		constants.add_constant(arg1->constant);
		return prior.log_unary_probability;
	} else {
		return std::numeric_limits<double>::infinity();
	}	
}

template<bool Quantified, typename ConstantDistribution>
double log_probability_atom(
		const hol_term* function, const hol_term* arg1, const hol_term* arg2,
		const simple_hol_term_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type != hol_term_type::CONSTANT)
		return std::numeric_limits<double>::infinity();

	constants.add_predicate(function->constant);
	if (arg1->type == hol_term_type::VARIABLE) {
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else if (arg2->type == hol_term_type::CONSTANT) {
			constants.add_constant(arg2->constant);
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else {
			return std::numeric_limits<double>::infinity();
		}
	} else if (arg1->type == hol_term_type::CONSTANT) {
		constants.add_constant(arg1->constant);
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else { /* both `arg1` and `arg2` can't be constants */
			return std::numeric_limits<double>::infinity();
		}
	} else {
		return std::numeric_limits<double>::infinity();
	}
}

template<bool Quantified, typename ConstantDistribution>
double log_probability_literal(const hol_term* literal,
		const simple_hol_term_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (literal->type == hol_term_type::UNARY_APPLICATION) {
		return prior.log_positive_probability + log_probability_atom<Quantified>(literal->binary.left, literal->binary.right, prior, constants);
	} else if (literal->type == hol_term_type::BINARY_APPLICATION) {
		return prior.log_positive_probability + log_probability_atom<Quantified>(literal->ternary.first, literal->ternary.second, literal->ternary.third, prior, constants);
	} else if (literal->type == hol_term_type::NOT && literal->unary.operand->type == hol_term_type::UNARY_APPLICATION) {
		return prior.log_negation_probability + log_probability_atom<Quantified>(literal->unary.operand->binary.left, literal->unary.operand->binary.right, prior, constants);
	} else if (literal->type == hol_term_type::NOT && literal->unary.operand->type == hol_term_type::BINARY_APPLICATION) {
		return prior.log_negation_probability + log_probability_atom<Quantified>(literal->unary.operand->ternary.first, literal->unary.operand->ternary.second, literal->unary.operand->ternary.third, prior, constants);
	} else {
		return std::numeric_limits<double>::infinity();
	}
}

template<typename ConstantDistribution>
double log_probability_helper(const hol_term* term,
		const simple_hol_term_distribution<ConstantDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	double value;
	const hol_term* antecedent;
	const hol_term* consequent;
	switch (term->type) {
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::NOT:
		return prior.log_ground_literal_probability + log_probability_literal<false>(term, prior, constants);
	case hol_term_type::FOR_ALL:
		if (term->quantifier.operand->type == hol_term_type::IF_THEN) {
			antecedent = term->quantifier.operand->binary.left;
			consequent = term->quantifier.operand->binary.right;
			if (antecedent->type == hol_term_type::AND) {
				value = (antecedent->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < antecedent->array.length; i++)
					value += log_probability_literal<true>(antecedent->array.operands[i], prior, constants);
			} else {
				value = prior.log_antecedent_stop_probability + log_probability_literal<true>(antecedent, prior, constants);
			}

			if (consequent->type == hol_term_type::AND) {
				value = (consequent->array.length - 1) * prior.log_consequent_continue_probability + prior.log_consequent_stop_probability;
				for (unsigned int i = 0; i < consequent->array.length; i++)
					value += log_probability_literal<true>(consequent->array.operands[i], prior, constants);
			} else {
				value = prior.log_consequent_stop_probability + log_probability_literal<true>(consequent, prior, constants);
			}
			return value;
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IF_THEN:
	case hol_term_type::IFF:
	case hol_term_type::EQUALS:
	case hol_term_type::EXISTS:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::LAMBDA:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::INTEGER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
	case hol_term_type::ANY_QUANTIFIER:
	case hol_term_type::VARIABLE_PREIMAGE:
		return -std::numeric_limits<double>::infinity();
	}
	fprintf(stderr, "log_probability ERROR: Unrecognized hol_term_type.\n");
	exit(EXIT_FAILURE);
}

template<typename ConstantDistribution,
	template<typename> class ObservationCollection,
	template<typename> class Collection>
double log_probability_ratio(
		const ObservationCollection<hol_term*>& observations,
		const simple_hol_term_distribution<ConstantDistribution>& prior,
		const Collection<hol_term*>& old_clusters,
		const Collection<hol_term*>& new_clusters)
{
	double value = 0.0;
	typename ConstantDistribution::ObservationCollection old_constants, new_constants;
	for (hol_term* formula : new_clusters)
		value += log_probability_helper(formula, prior, new_constants);
	for (hol_term* formula : old_clusters)
		value -= log_probability_helper(formula, prior, old_constants);
	return value + log_probability_ratio(old_constants, new_constants, prior.constant_distribution);
}

template<typename T>
struct uniform_subset_distribution {
	double log_stop_probability;
	double log_continue_probability;

	typedef default_array<pair<array_view<const T>, array_view<const T>>> ObservationCollection;

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
		const pair<Collection<const T>, Collection<const T>>& observation,
		const uniform_subset_distribution<T>& prior)
{
	if (observation.key.size() < 2)
		return -std::numeric_limits<double>::infinity();
	log_cache<double>::instance().ensure_size(observation.value.size() + 1);
	return -log_cache<double>::instance().get(observation.value.size()) * observation.key.size()
		 + (observation.key.size() - 2) * prior.log_continue_probability + prior.log_stop_probability;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
double log_probability(
		const ObservationCollection<pair<Collection<const T>, Collection<const T>>>& observations,
		const uniform_subset_distribution<T>& prior)
{
	double value = 0.0;
	for (const pair<Collection<const T>, Collection<const T>>& observation : observations)
		value += log_probability(observation, prior);
	return value;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
inline double log_probability_ratio(
		const ObservationCollection<pair<Collection<const T>, Collection<const T>>>& old_observations,
		const ObservationCollection<pair<Collection<const T>, Collection<const T>>>& new_observations,
		const uniform_subset_distribution<T>& prior)
{
	return log_probability(new_observations, prior) - log_probability(old_observations, prior);
}

template<typename T>
struct unif_distribution {
	typedef default_array<pair<T, array<T>>> ObservationCollection;
};

template<template<typename> class Collection, typename T>
double log_probability(
		const pair<T, Collection<T>>& observation,
		const unif_distribution<T>& prior)
{
	log_cache<double>::instance().ensure_size(size(observation.value) + 1);
	return -log_cache<double>::instance().get(size(observation.value));
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
double log_probability(
		const ObservationCollection<pair<T, Collection<T>>>& observations,
		const unif_distribution<T>& prior)
{
	double value = 0.0;
	for (const pair<T, Collection<T>>& observation : observations)
		value += log_probability(observation, prior);
	return value;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
inline double log_probability_ratio(
		const ObservationCollection<pair<T, Collection<T>>>& old_observations,
		const ObservationCollection<pair<T, Collection<T>>>& new_observations,
		const unif_distribution<T>& prior)
{
	return log_probability(new_observations, prior) - log_probability(old_observations, prior);
}

template<typename JumpDistribution, typename LengthDistribution>
struct levy_process {
	typedef typename JumpDistribution::ObservationType ValueType;
	typedef array_view<ValueType> ObservationType;
	typedef default_array<ObservationType> ObservationCollection;

	JumpDistribution jump_distribution;
	LengthDistribution length_distribution;

	levy_process(
			const JumpDistribution& jump_distribution,
			const LengthDistribution& length_distribution) :
		jump_distribution(jump_distribution),
		length_distribution(length_distribution)
	{ }
};

template<typename JumpDistribution, typename LengthDistribution>
inline levy_process<JumpDistribution, LengthDistribution> make_levy_process(
		const JumpDistribution& jump_distribution,
		const LengthDistribution& length_distribution)
{
	return levy_process<JumpDistribution, LengthDistribution>(
		jump_distribution, length_distribution);
}

template<template<typename> class Collection,
	typename JumpDistribution, typename LengthDistribution>
double log_probability(
		const Collection<typename JumpDistribution::ObservationType>& observation,
		const levy_process<JumpDistribution, LengthDistribution>& prior)
{
	double value = log_probability(observation.size() - 1, prior.length_distribution);
	typename JumpDistribution::ObservationType prev = 0;
	for (const auto& entry : observation) {
		value += log_probability(entry - prev, prior.jump_distribution);
		prev = entry;
	}
	return value;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection,
	typename JumpDistribution, typename LengthDistribution>
double log_probability(
		const ObservationCollection<Collection<typename JumpDistribution::ObservationType>>& observations,
		const levy_process<JumpDistribution, LengthDistribution>& prior)
{
	double value = 0.0;
	for (const auto& observation : observations)
		value += log_probability(observation, prior);
	return value;
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection,
	typename JumpDistribution, typename LengthDistribution>
inline double log_probability_ratio(
		const ObservationCollection<Collection<typename JumpDistribution::ObservationType>>& old_observations,
		const ObservationCollection<Collection<typename JumpDistribution::ObservationType>>& new_observations,
		const levy_process<JumpDistribution, LengthDistribution>& prior)
{
	return log_probability(new_observations, prior) - log_probability(old_observations, prior);
}

const string* get_name(const hash_map<string, unsigned int>& names, unsigned int id)
{
	for (const auto& entry : names)
		if (entry.value == id) return &entry.key;
	return NULL;
}

template<typename Formula, typename ProofCalculus, typename Canonicalizer>
inline bool contains_subset_axiom(
		const theory<Formula, ProofCalculus, Canonicalizer>& T,
		const Formula* antecedent, const Formula* consequent)
{
	bool contains;
	unsigned int antecedent_set = T.sets.set_ids.get(*antecedent, contains);
	if (!contains) return false;
	unsigned int consequent_set = T.sets.set_ids.get(*consequent, contains);
	if (!contains) return false;

	return T.sets.extensional_graph.vertices[antecedent_set].parents.contains(consequent_set)
		|| T.sets.intensional_graph.vertices[antecedent_set].parents.contains(consequent_set);
}

template<typename Formula, typename ProofCalculus, typename Canonicalizer>
bool contains_axiom(const theory<Formula, ProofCalculus, Canonicalizer>& T, const Formula* formula)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::Term Term;

	if (formula == NULL) {
		fprintf(stderr, "contains_axiom WARNING: `formula` is null.\n");
		return false;
	}

	bool contains;
	unsigned int predicate = 0; Term const* arg1 = nullptr; Term const* arg2 = nullptr;
	if (formula->type == FormulaType::FOR_ALL) {
		if (formula->quantifier.operand->type == FormulaType::IF_THEN)
			return contains_subset_axiom(T, formula->quantifier.operand->binary.left, formula->quantifier.operand->binary.right);
		else return false;
	} else if (is_atomic(*formula, predicate, arg1, arg2)) {
		if (arg2 == NULL) {
			/* `formula` is an atom of form `f(a)` */
			const pair<array<unsigned int>, array<unsigned int>>& types = T.types.get(predicate, contains);
			if (!contains) return false;
			return types.key.contains(arg1->constant);
		} else {
			/* `formula` is an atom of form `f(a,b)` */
			relation rel = { 0, arg1->constant, arg2->constant };
			const pair<array<unsigned int>, array<unsigned int>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return relations.key.contains(predicate);
		}
	} else if (formula->type == FormulaType::NOT && is_atomic(*formula->unary.operand)) {
		if (arg2 == NULL) {
			/* `formula` is an atom of form `~f(a)` */
			const pair<array<unsigned int>, array<unsigned int>>& types = T.types.get(predicate, contains);
			if (!contains) return false;
			return types.value.contains(arg1->constant);
		} else {
			/* `formula` is an atom of form `~f(a,b)` */
			relation rel = { 0, arg1->constant, arg2->constant };
			const pair<array<unsigned int>, array<unsigned int>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return relations.value.contains(predicate);
		}
	} else {
		return false;
	}
}

unsigned int constant_offset = 0;

template<typename Stream>
bool print_special_string(unsigned int key, Stream& out) {
	return print('c', out) && print_subscript(key - constant_offset, out);
}

int main(int argc, const char** argv)
{
	setlocale(LC_ALL, "en_US.UTF-8");
	log_cache<double>::instance().ensure_size(1024);

	hash_map<string, unsigned int> names(256);
	if (!add_constants_to_string_map(names)) {
		return EXIT_FAILURE;
	}

	/* construct the parser */
	hdp_parser<hol_term> parser = hdp_parser<hol_term>(
			(unsigned int) built_in_predicates::UNKNOWN,
			names, "english.morph", "english.gram");

	/* read the seed training set of sentences labeled with logical forms */
	FILE* in = fopen("seed_training_set.txt", "rb");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open file for reading.\n");
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

	array<article_token> tokens = array<article_token>(256);
	if (!article_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis of training data failed.\n");
		for (auto entry : names) free(entry.key);
		fclose(in); free_tokens(tokens); return EXIT_FAILURE;
	}
	fclose(in);

	unsigned int index = 0;
	typedef sentence<syntax_node<flagged_logical_form<hol_term>>> sentence_type;
	array<array_map<sentence_type, hol_term>> seed_training_set(64);
	while (index < tokens.length) {
		if (!seed_training_set.ensure_capacity(seed_training_set.length + 1)
		 || !array_map_init(seed_training_set[seed_training_set.length], 4))
		{
			for (array_map<sentence_type, hol_term>& paragraph : seed_training_set) {
				for (auto entry : paragraph) { free(entry.key); free(entry.value); }
				free(paragraph);
			}
			free_tokens(tokens);
			for (auto entry : names) free(entry.key);
			return EXIT_FAILURE;
		}
		seed_training_set.length++;
		if (!article_interpret(tokens, index, seed_training_set.last(), names, parser.G.nonterminal_names)) {
			fprintf(stderr, "ERROR: Failed to parse training data.\n");
			for (array_map<sentence_type, hol_term>& paragraph : seed_training_set) {
				for (auto entry : paragraph) { free(entry.key); free(entry.value); }
				free(paragraph);
			}
			free_tokens(tokens);
			for (auto entry : names) free(entry.key);
			return EXIT_FAILURE;
		}
	}
	free_tokens(tokens); tokens.clear();

	/* train the parser */
	if (!parser.train(seed_training_set, names, 10)) {
		for (auto entry : names) free(entry.key);
		for (array_map<sentence_type, hol_term>& paragraph : seed_training_set) {
			for (auto entry : paragraph) { free(entry.key); free(entry.value); }
			free(paragraph);
		}
		return EXIT_FAILURE;
	}
	for (array_map<sentence_type, hol_term>& paragraph : seed_training_set) {
		for (auto entry : paragraph) { free(entry.key); free(entry.value); }
		free(paragraph);
	}

fprintf(stderr, "Parser constructed and trained. Exiting...\n");
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;

	/* read the articles */
	in = fopen("simple_set_reasoning_articles.txt", "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open file for reading.\n");
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

	if (!article_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis of article failed.\n");
		for (auto entry : names) free(entry.key);
		fclose(in); free_tokens(tokens); return EXIT_FAILURE;
	}
	fclose(in);

	index = 0;
	typedef article<syntax_node<flagged_logical_form<hol_term>>> article_type;
	typedef in_memory_article_store<syntax_node<flagged_logical_form<hol_term>>> article_store_type;
	article_store_type corpus;
	while (index < tokens.length) {
		unsigned int article_name = 0;
		article_type& new_article = *((article_type*) alloca(sizeof(article_type)));
		if (!corpus.articles.check_size() || !article_interpret(tokens, index, new_article, article_name, names, parser.G.nonterminal_names)) {
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
		article_type& value = corpus.articles.get(article_name, contains, bucket);
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
	theory<hol_term, natural_deduction<hol_term>, standard_canonicalizer<true, false>> T(names.table.size + 1);
	constant_offset = T.new_constant_offset;
	auto constant_prior = make_simple_constant_distribution(
			chinese_restaurant_process<unsigned int>(1.0), chinese_restaurant_process<unsigned int>(1.0));
	auto theory_element_prior = make_simple_hol_term_distribution(constant_prior, 0.01, 0.3, 0.4, 0.2, 0.4);
	auto axiom_prior = make_dirichlet_process(1.0e-3, theory_element_prior);
	auto conjunction_prior = uniform_subset_distribution<const nd_step<hol_term>*>(0.1);
	auto universal_introduction_prior = unif_distribution<unsigned int>();
	auto universal_elimination_prior = chinese_restaurant_process<hol_term>(1.0);
	auto term_indices_prior = make_levy_process(poisson_distribution(1.0), poisson_distribution(1.0));
	auto proof_prior = make_canonicalized_proof_prior(axiom_prior, conjunction_prior,
			universal_introduction_prior, universal_elimination_prior, term_indices_prior, poisson_distribution(5.0));
	auto theory_prior = make_theory_prior(proof_prior, geometric_distribution(0.2));
	const string** reverse_name_map = invert(names);
	string_map_scribe printer = { reverse_name_map, names.table.size + 1 };

	bool success = true;
	success &= read_article(names.get("Nemo"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Dory"), corpus, parser, T, theory_prior, printer);
	//success &= read_article(names.get("red"), corpus, parser, T, theory_prior, printer);
	//success &= read_article(names.get("blue"), corpus, parser, T, theory_prior, printer);
	//success &= read_article(names.get("red_or_blue"), corpus, parser, T, theory_prior, printer);
	//success &= read_article(names.get("red_and_blue"), corpus, parser, T, theory_prior, printer);
	/*success &= read_article(names.get("Bob"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Kate"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Sam"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Byron"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Alex"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Lee"), corpus, parser, T, theory_prior, printer);
	success &= read_article(names.get("Amy"), corpus, parser, T, theory_prior, printer);*/

	array_map<hol_term*, unsigned int> tracked_logical_forms(2);
	/*tracked_logical_forms.put(
		hol_term::new_for_all(1,
			hol_term::new_if_then(
				hol_term::new_atom(parser.symbol_map.get(names.get("cat")), hol_term::new_variable(1)),
				hol_term::new_atom(parser.symbol_map.get(names.get("mammal")), hol_term::new_variable(1))
			)
		), 0);
	tracked_logical_forms.put(
		hol_term::new_for_all(1,
			hol_term::new_if_then(
				hol_term::new_atom(parser.symbol_map.get(names.get("mammal")), hol_term::new_variable(1)),
				hol_term::new_atom(parser.symbol_map.get(names.get("cat")), hol_term::new_variable(1))
			)
		), 0);*/

	unsigned int iterations = (success ? 120000 : 0);
	timer stopwatch;
	auto scribe = parser.get_printer(printer);
array_multiset<unsigned int> set_size_distribution(16);
	for (unsigned int t = 0; t < iterations; t++) {
		//printf("[%u]\n", t);
T.check_proof_axioms();
T.sets.are_elements_unique();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
		//T.print_axioms(stdout, scribe);
		//print('\n', stdout); fflush(stdout);
		if (stopwatch.milliseconds() > 1000) {
			print("[iteration ", stdout); print(t, stdout); print("]\n", stdout);
			for (const auto& entry : tracked_logical_forms) {
				print("p(", stdout); print(*entry.key, stdout, scribe); print(" axiom) ≈ ", stdout); print((double) entry.value / t, stdout); print('\n', stdout);
			}
for (const auto& entry : set_size_distribution.counts)
fprintf(stderr, "%u %lf\n", entry.key, (double) entry.value / set_size_distribution.sum);
			stopwatch.start();
		}
if (t == 7)
fprintf(stderr, "DEBUG: BREAKPOINT\n");
		do_mh_step(T, theory_prior);

		for (auto entry : tracked_logical_forms)
			if (contains_axiom(T, entry.key)) entry.value++;
hol_term* set_formula = hol_term::new_atom(names.get("fish"), hol_term::new_variable(1));
set_size_distribution.add(T.sets.sets[T.sets.set_ids.get(*set_formula)].set_size);
free(*set_formula); free(set_formula);
	}

FILE* histogram_file = fopen("histogram.txt", "w");
for (const auto& entry : set_size_distribution.counts)
fprintf(histogram_file, "%u %lf\n", entry.key, (double) entry.value / set_size_distribution.sum);
fflush(histogram_file); fclose(histogram_file);

	print("[iteration ", stdout); print(iterations, stdout); print("]\n", stdout);
	for (const auto& entry : tracked_logical_forms) {
		print("p(", stdout); print(*entry.key, stdout, scribe); print(" axiom) ≈ ", stdout); print((double) entry.value / iterations, stdout); print('\n', stdout);
	}
for (const auto& entry : set_size_distribution.counts)
fprintf(stderr, "%u %lf\n", entry.key, (double) entry.value / set_size_distribution.sum);

	for (auto entry : tracked_logical_forms) {
		free(*entry.key); if (entry.key->reference_count == 0) free(entry.key);
	}
	free(reverse_name_map);
	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
