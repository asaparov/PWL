#include "higher_order_logic.h"
#include "hdp_parser.h"
#include "executive.h"

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
struct dummy_array {
	dummy_array() { }

	inline constexpr bool add(const T& item) const { return true; }
};

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
struct default_hash_multiset {
	hash_multiset<T, false> a;

	default_hash_multiset() : a(16) { }

	inline bool add(const T& item) {
		return a.add(item);
	}

	inline bool add(const default_array_multiset<T>& changes) {
		return a.add(changes.a);
	}

	inline void subtract(const default_array_multiset<T>& changes) {
		a.template subtract<true>(changes.a);
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
	typedef default_hash_multiset<T> PriorState;
	typedef default_array_multiset<T> PriorStateChanges;

	double alpha, log_alpha;

	chinese_restaurant_process(double alpha) :
		alpha(alpha), log_alpha(log(alpha))
	{ }

	chinese_restaurant_process(const chinese_restaurant_process<T>& src) :
		alpha(src.alpha), log_alpha(log(alpha))
	{ }

	~chinese_restaurant_process() { }
};

template<typename MultisetType, typename T>
double log_probability(
		const MultisetType& observations,
		const chinese_restaurant_process<T>& prior)
{
	double value = 0.0;
	for (unsigned int i = 0; i < size(observations); i++)
		value += prior.log_alpha + lgamma(get_value(observations, i));
	return value + lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
}

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

template<typename T, bool AutomaticallyFree, typename MultisetType>
inline double log_probability_ratio(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior)
{
	dummy_array<T> dummy;
	return log_probability_ratio(tables, old_observations, new_observations, prior, dummy, dummy);
}

template<typename T, typename MultisetType>
inline double log_probability_ratio(
		const default_hash_multiset<T>& tables,
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior)
{
	return log_probability_ratio(tables.a, old_observations, new_observations, prior);
}

template<typename BaseDistribution>
struct dirichlet_process
{
	struct prior_state {
		typename BaseDistribution::PriorState base_prior_state;

		inline bool add(const typename BaseDistribution::PriorStateChanges& changes) {
			return base_prior_state.add(changes);
		}

		inline bool add(
				const typename BaseDistribution::ObservationType& observation,
				const dirichlet_process<BaseDistribution>& prior)
		{
			return base_prior_state.add(observation, prior.base_distribution);
		}

		inline void subtract(const typename BaseDistribution::PriorStateChanges& changes) {
			base_prior_state.subtract(changes);
		}

		inline void subtract(
				const typename BaseDistribution::ObservationType& observation,
				const dirichlet_process<BaseDistribution>& prior)
		{
			base_prior_state.subtract(observation, prior.base_distribution);
		}
	};

	typedef typename BaseDistribution::ObservationType ObservationType;
	typedef prior_state PriorState;
	typedef typename BaseDistribution::PriorStateChanges PriorStateChanges;

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

template<typename MultisetType, typename BaseDistribution>
double log_probability(
		const MultisetType& restaurant_tables,
		const dirichlet_process<BaseDistribution>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters clusters;
	return log_probability(restaurant_tables, prior.restaurant, clusters)
		 + log_probability(clusters, prior.base_distribution);
}

template<typename MultisetType, typename T, bool AutomaticallyFree, typename BaseDistribution>
double log_probability_ratio(
		const MultisetType& restaurant_tables,
		const array_multiset<T, AutomaticallyFree>& old_observations,
		const array_multiset<T, AutomaticallyFree>& new_observations,
		const dirichlet_process<BaseDistribution>& prior,
		const typename dirichlet_process<BaseDistribution>::prior_state& prior_state,
		typename BaseDistribution::PriorStateChanges& old_prior_changes,
		typename BaseDistribution::PriorStateChanges& new_prior_changes)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters old_clusters, new_clusters;
	return log_probability_ratio(restaurant_tables, old_observations, new_observations, prior.restaurant, old_clusters, new_clusters)
		 + log_probability_ratio(old_clusters, new_clusters, prior.base_distribution, prior_state.base_prior_state, old_prior_changes, new_prior_changes);
}

template<typename ConstantDistribution, typename PredicateDistribution>
struct simple_constant_distribution
{
	struct prior_state_changes {
		typename ConstantDistribution::PriorStateChanges constants;
		typename PredicateDistribution::PriorStateChanges predicates;

		inline bool add_constant(unsigned int constant) {
			return constants.add(constant);
		}

		inline bool add_predicate(unsigned int predicate) {
			return predicates.add(predicate);
		}
	};

	struct prior_state {
		typename ConstantDistribution::PriorState constants;
		typename PredicateDistribution::PriorState predicates;

		inline bool add(const prior_state_changes& changes) {
			if (!constants.add(changes.constants)) {
				return false;
			} else if (!predicates.add(changes.predicates)) {
				constants.subtract(changes.constants);
				return false;
			}
			return true;
		}

		inline bool add(const prior_state_changes& changes,
				const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& prior)
		{
			return add(changes);
		}

		inline void subtract(const prior_state_changes& changes) {
			constants.subtract(changes.constants);
			predicates.subtract(changes.predicates);
		}

		inline void subtract(const prior_state_changes& changes,
				const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& prior)
		{
			subtract(changes);
		}
	};

	typedef prior_state PriorState;
	typedef prior_state_changes PriorStateChanges;
	typedef prior_state_changes ObservationCollection;

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
inline double log_probability(const ConstantCollection& constants,
		const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& prior)
{
	return log_probability(constants.constants, prior.constant_distribution)
		 + log_probability(constants.predicates, prior.predicate_distribution);
}

template<typename ObservationCollection, typename ConstantCollection,
	typename ConstantDistribution, typename PredicateDistribution>
inline double log_probability_ratio(const ObservationCollection& existing_constants,
		const ConstantCollection& old_constants, const ConstantCollection& new_constants,
		const simple_constant_distribution<ConstantDistribution, PredicateDistribution>& prior)
{
	return log_probability_ratio(existing_constants.constants, old_constants.constants, new_constants.constants, prior.constant_distribution)
		 + log_probability_ratio(existing_constants.predicates, old_constants.predicates, new_constants.predicates, prior.predicate_distribution);
}

template<typename ConstantDistribution, typename SetSizeDistribution>
struct simple_hol_term_distribution
{
	struct prior_state {
		typename ConstantDistribution::PriorState constant_prior_state;

		inline bool add(const typename ConstantDistribution::PriorStateChanges& changes) {
			return constant_prior_state.add(changes);
		}

		inline bool add(const hol_term* observation,
				const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior)
		{
			typename ConstantDistribution::ObservationCollection constants;
			log_probability_helper(observation, prior, constants);
			return constant_prior_state.add(constants, prior.constant_distribution);
		}

		inline void subtract(const typename ConstantDistribution::PriorStateChanges& changes) {
			constant_prior_state.subtract(changes);
		}

		inline void subtract(const hol_term* observation,
				const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior)
		{
			typename ConstantDistribution::ObservationCollection constants;
			log_probability_helper(observation, prior, constants);
			constant_prior_state.subtract(constants, prior.constant_distribution);
		}
	};

	typedef hol_term* ObservationType;
	typedef default_array<hol_term*> ObservationCollection;
	typedef prior_state PriorState;
	typedef typename ConstantDistribution::PriorStateChanges PriorStateChanges;

	double log_ground_literal_probability;
	double log_universal_probability;
	double log_existential_probability;
	double log_disjunction_probability;
	double log_set_size_axiom_probability;

	double log_negation_probability;
	double log_positive_probability;

	double log_unary_probability;
	double log_binary_probability;

	double log_antecedent_continue_probability;
	double log_antecedent_stop_probability;
	double log_consequent_continue_probability;
	double log_consequent_stop_probability;

	ConstantDistribution constant_distribution;
	SetSizeDistribution set_size_distribution;

	simple_hol_term_distribution(
			const ConstantDistribution& constant_distribution,
			const SetSizeDistribution& set_size_distribution,
			double ground_literal_probability, double universal_probability,
			double existential_probability, double disjunction_probability,
			double set_size_axiom_probability, double negation_probability,
			double unary_probability, double antecedent_stop_probability,
			double consequent_stop_probability) :
		log_ground_literal_probability(log(ground_literal_probability)),
		log_universal_probability(log(universal_probability)),
		log_existential_probability(log(existential_probability)),
		log_disjunction_probability(log(disjunction_probability)),
		log_set_size_axiom_probability(log(set_size_axiom_probability)),
		log_negation_probability(log(negation_probability)),
		log_positive_probability(log(1.0 - negation_probability)),
		log_unary_probability(log(unary_probability)),
		log_binary_probability(log(1.0 - unary_probability)),
		log_antecedent_continue_probability(log(1.0 - antecedent_stop_probability)),
		log_antecedent_stop_probability(log(antecedent_stop_probability)),
		log_consequent_continue_probability(log(1.0 - consequent_stop_probability)),
		log_consequent_stop_probability(log(consequent_stop_probability)),
		constant_distribution(constant_distribution),
		set_size_distribution(set_size_distribution)
	{
		if (fabs(ground_literal_probability + universal_probability + existential_probability + disjunction_probability + set_size_axiom_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `ground_literal_probability + universal_probability + existential_probability + disjunction_probability + set_size_axiom_probability` is not 1.\n");
	}

	simple_hol_term_distribution(const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& src) :
		log_ground_literal_probability(src.log_ground_literal_probability),
		log_universal_probability(src.log_universal_probability),
		log_existential_probability(src.log_existential_probability),
		log_disjunction_probability(src.log_disjunction_probability),
		log_set_size_axiom_probability(src.log_set_size_axiom_probability),
		log_negation_probability(src.log_negation_probability),
		log_positive_probability(src.log_positive_probability),
		log_unary_probability(src.log_unary_probability),
		log_binary_probability(src.log_binary_probability),
		log_antecedent_continue_probability(src.log_antecedent_continue_probability),
		log_antecedent_stop_probability(src.log_antecedent_stop_probability),
		log_consequent_continue_probability(src.log_consequent_continue_probability),
		log_consequent_stop_probability(src.log_consequent_stop_probability),
		constant_distribution(src.constant_distribution),
		set_size_distribution(src.set_size_distribution)
	{ }
};

template<typename ConstantDistribution, typename SetSizeDistribution>
inline simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution> make_simple_hol_term_distribution(
		const ConstantDistribution& constant_distribution,
		const SetSizeDistribution set_size_distribution,
		double ground_literal_probability, double universal_probability,
		double existential_probability, double disjunction_probability,
		double set_size_axiom_probability, double negation_probability,
		double unary_probability, double antecedent_stop_probability,
		double consequent_stop_probability)
{
	return simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>(
			constant_distribution, set_size_distribution,
			ground_literal_probability, universal_probability,
			existential_probability, disjunction_probability,
			set_size_axiom_probability, negation_probability,
			unary_probability, antecedent_stop_probability,
			consequent_stop_probability);
}

template<bool Quantified, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_atom(const hol_term* function, const hol_term* arg1,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type != hol_term_type::CONSTANT)
		return -std::numeric_limits<double>::infinity();

	constants.add_predicate(function->constant);
	if (Quantified && arg1->type == hol_term_type::VARIABLE) {
		return prior.log_unary_probability;
	} else if (!Quantified && arg1->type == hol_term_type::CONSTANT) {
		constants.add_constant(arg1->constant);
		return prior.log_unary_probability;
	} else {
		return -std::numeric_limits<double>::infinity();
	}	
}

template<bool Quantified, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_atom(
		const hol_term* function, const hol_term* arg1, const hol_term* arg2,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type != hol_term_type::CONSTANT)
		return -std::numeric_limits<double>::infinity();

	constants.add_predicate(function->constant);
	if (arg1->type == hol_term_type::VARIABLE) {
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else if (arg2->type == hol_term_type::CONSTANT) {
			constants.add_constant(arg2->constant);
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	} else if (arg1->type == hol_term_type::CONSTANT) {
		constants.add_constant(arg1->constant);
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(3);
		} else { /* both `arg1` and `arg2` can't be constants */
			return -std::numeric_limits<double>::infinity();
		}
	} else {
		return -std::numeric_limits<double>::infinity();
	}
}

template<bool Quantified, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_literal(const hol_term* literal,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior,
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
		return -std::numeric_limits<double>::infinity();
	}
}

inline bool is_literal(const hol_term* term) {
	if (term->type == hol_term_type::UNARY_APPLICATION || term->type == hol_term_type::BINARY_APPLICATION)
		return true;
	if (term->type == hol_term_type::NOT && is_literal(term->unary.operand))
		return true;
	return false;
}

template<bool Quantified = false, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_helper(const hol_term* term,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (is_literal(term))
		return prior.log_ground_literal_probability + log_probability_literal<Quantified>(term, prior, constants);

	double value;
	const hol_term* antecedent;
	const hol_term* consequent;
	switch (term->type) {
	case hol_term_type::FOR_ALL:
		value = prior.log_universal_probability;
		if (term->quantifier.operand->type == hol_term_type::IF_THEN) {
			antecedent = term->quantifier.operand->binary.left;
			consequent = term->quantifier.operand->binary.right;
			if (antecedent->type == hol_term_type::AND) {
				value = (antecedent->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < antecedent->array.length; i++)
					value += log_probability_helper<true>(antecedent->array.operands[i], prior, constants);
			} else {
				value = prior.log_antecedent_stop_probability + log_probability_helper<true>(antecedent, prior, constants);
			}

			if (consequent->type == hol_term_type::AND) {
				value = (consequent->array.length - 1) * prior.log_consequent_continue_probability + prior.log_consequent_stop_probability;
				for (unsigned int i = 0; i < consequent->array.length; i++)
					value += log_probability_helper<true>(consequent->array.operands[i], prior, constants);
			} else {
				value = prior.log_consequent_stop_probability + log_probability_helper<true>(consequent, prior, constants);
			}
			return value;
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	case hol_term_type::EQUALS:
		if (term->binary.left->type == hol_term_type::UNARY_APPLICATION
		 && term->binary.left->binary.left->type == hol_term_type::CONSTANT
		 && term->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
		 && term->binary.left->binary.right->type == hol_term_type::LAMBDA
		 && term->binary.right->type == hol_term_type::INTEGER)
		{
			value = prior.log_set_size_axiom_probability;
			hol_term* operand = term->binary.left->binary.right->quantifier.operand;
			if (operand->type == hol_term_type::AND) {
				value = (operand->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < operand->array.length; i++)
					value += log_probability_helper<true>(operand->array.operands[i], prior, constants);
			} else {
				value = prior.log_antecedent_stop_probability + log_probability_helper<true>(operand, prior, constants);
			}
			value += log_probability(term->binary.right->integer, prior.set_size_distribution);
			return value;
		} else {
			// TODO: implement this (we also have definition axioms for constants)
			return prior.log_set_size_axiom_probability;
			/*fprintf(stderr, "log_probability_helper WARNING: Found equality term that is not a set size axiom: ");
			print(*term, stderr); print('\n', stderr);
			return -std::numeric_limits<double>::infinity();*/
		}
	case hol_term_type::EXISTS:
		value = prior.log_existential_probability;
		if (term->quantifier.operand->type == hol_term_type::AND) {
			hol_term* operand = term->quantifier.operand;
			value += (operand->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
			for (unsigned int i = 0; i < operand->array.length; i++)
				value += log_probability_helper<true>(operand->array.operands[i], prior, constants);
		} else {
			value += prior.log_antecedent_stop_probability + log_probability_helper<true>(term->quantifier.operand, prior, constants);
		}
		return value;
	case hol_term_type::OR:
		value = prior.log_disjunction_probability;
		value += (term->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
		for (unsigned int i = 0; i < term->array.length; i++)
			value += log_probability_helper<true>(term->array.operands[i], prior, constants);
		return value;
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::NOT:
	case hol_term_type::AND:
	case hol_term_type::IF_THEN:
	case hol_term_type::IFF:
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
	case hol_term_type::ANY_RIGHT_ONLY:
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
	typename SetSizeDistribution,
	template<typename> class Collection>
double log_probability(const Collection<hol_term*>& clusters,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior)
{
	double value = 0.0;
	typename ConstantDistribution::ObservationCollection constants;
	for (hol_term* formula : clusters)
		value += log_probability_helper(formula, prior, constants);
	return value + log_probability(constants, prior.constant_distribution);
}

template<typename ConstantDistribution,
	typename SetSizeDistribution,
	template<typename> class Collection>
double log_probability_ratio(
		const Collection<hol_term*>& old_clusters,
		const Collection<hol_term*>& new_clusters,
		const simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>& prior,
		const typename simple_hol_term_distribution<ConstantDistribution, SetSizeDistribution>::prior_state& prior_state,
		typename ConstantDistribution::PriorStateChanges& old_prior_changes,
		typename ConstantDistribution::PriorStateChanges& new_prior_changes)
{
	double value = 0.0;
	for (hol_term* formula : new_clusters)
		value += log_probability_helper(formula, prior, new_prior_changes);
	for (hol_term* formula : old_clusters)
		value -= log_probability_helper(formula, prior, old_prior_changes);
	return value + log_probability_ratio(prior_state.constant_prior_state, old_prior_changes, new_prior_changes, prior.constant_distribution);
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

template<typename ProofCalculus, typename Canonicalizer>
inline bool contains_subset_axiom(
		const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* antecedent,
		const typename ProofCalculus::Language* consequent)
{
	bool contains;
	unsigned int antecedent_set = T.sets.set_ids.get(*antecedent, contains);
	if (!contains) return false;
	unsigned int consequent_set = T.sets.set_ids.get(*consequent, contains);
	if (!contains) return false;

	return T.sets.extensional_graph.vertices[antecedent_set].parents.contains(consequent_set)
		|| T.sets.intensional_graph.vertices[antecedent_set].parents.contains(consequent_set);
}

template<typename ProofCalculus, typename Canonicalizer>
bool contains_axiom(
		const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula)
{
	typedef typename ProofCalculus::Language Formula;
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
			Term* atom = Term::new_apply(Term::new_constant(predicate), &Term::template variables<0>::value);
			if (atom == nullptr) return false;
			Term::template variables<0>::value.reference_count++;
			const pair<array<unsigned int>, array<unsigned int>>& types = T.atoms.get(*atom, contains);
			free(*atom); free(atom);
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
			Term* atom = Term::new_apply(Term::new_constant(predicate), &Term::template variables<0>::value);
			if (atom == nullptr) return false;
			Term::template variables<0>::value.reference_count++;
			const pair<array<unsigned int>, array<unsigned int>>& types = T.atoms.get(*atom, contains);
			free(*atom); free(atom);
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

template<typename Stream>
unsigned int read_line(array<char>& line, Stream& input)
{
	unsigned int bytes_read = 0;
	while (true) {
		int width;
		wint_t next = fgetwc(input);
		if (!line.ensure_capacity(line.length + MB_CUR_MAX)) {
			fprintf(stderr, "read_line ERROR: Out of memory.\n");
			return 0;
		}
		switch (next) {
		case WEOF:
			return bytes_read;

		case '\n':
			return bytes_read + 1;

		default:
#if defined(_WIN32)
			wctomb_s(&width, line.data + line.length, (line.capacity - line.length) * sizeof(char), next);
#else
			width = wctomb(line.data + line.length, next);
#endif
			if (width == -1)
				return 0;
			line.length += width;
			bytes_read += width;
		}
	}

	return bytes_read;
}

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

template<typename Stream, typename Parser>
inline void run_console(Stream& input, const char* prompt,
		Parser& parser, hash_map<string, unsigned int>& names)
{
	array<array_map<typename Parser::SentenceType, flagged_logical_form<hol_term>>> dummy_training_set(1);
	run_console(input, prompt, parser, names, dummy_training_set);
}

template<bool AllConstantsDistinct, bool PolymorphicEquality, typename BuiltInPredicates>
struct polymorphic_canonicalizer
{
	template<bool Quiet = false>
	static inline hol_term* canonicalize(const hol_term& src)
	{
		equals_arg_types<simple_type> types(16);
		array_map<unsigned int, hol_type<simple_type>> constant_types(8);
		array_map<unsigned int, hol_type<simple_type>> variable_types(8);
		array_map<unsigned int, hol_type<simple_type>> parameter_types(8);

		constant_types.keys[0] = (unsigned int) BuiltInPredicates::ARG1;
		constant_types.keys[1] = (unsigned int) BuiltInPredicates::ARG2;
		constant_types.keys[2] = (unsigned int) BuiltInPredicates::ARG3;
		if (!init(constant_types.values[0], 1, hol_type<simple_type>(2, hol_type<simple_type>(hol_type<simple_type>(1), hol_type<simple_type>(2))))) {
			return nullptr;
		} else if (!init(constant_types.values[1], 1, hol_type<simple_type>(2, hol_type<simple_type>(hol_type<simple_type>(1), hol_type<simple_type>(2))))) {
			free(constant_types.values[0]);
			return nullptr;
		} else if (!init(constant_types.values[2], 1, hol_type<simple_type>(2, hol_type<simple_type>(hol_type<simple_type>(1), hol_type<simple_type>(2))))) {
			free(constant_types.values[0]);
			free(constant_types.values[1]);
			return nullptr;
		}
		constant_types.size += 3;

		bool success = compute_type<PolymorphicEquality, Quiet>(src, types, constant_types, variable_types, parameter_types);
		for (unsigned int j = 0; j < constant_types.size; j++) free(constant_types.values[j]);
		for (unsigned int j = 0; j < variable_types.size; j++) free(variable_types.values[j]);
		for (unsigned int j = 0; j < parameter_types.size; j++) free(parameter_types.values[j]);
		if (!success) return nullptr;

		array_map<unsigned int, unsigned int> variable_map(16);
		hol_scope& scope = *((hol_scope*) alloca(sizeof(hol_scope)));
		if (!canonicalize_scope<AllConstantsDistinct>(src, scope, variable_map, types))
			return nullptr;
		hol_term* canonicalized = scope_to_term(scope);
		free(scope);
		return canonicalized;
	}
};


unsigned int constant_offset = 0;

template<typename Stream>
bool print_special_string(unsigned int key, Stream& out) {
	return print('c', out) && print_subscript(key - constant_offset, out);
}

int main(int argc, const char** argv)
{
	setlocale(LC_ALL, "en_US.UTF-8");
	log_cache<double>::instance().ensure_size(1024);
set_seed(1356941742);
	fprintf(stdout, "(seed = %u)\n", get_seed());

	hash_map<string, unsigned int> names(256);
	if (!add_constants_to_string_map(names)) {
		return EXIT_FAILURE;
	}

	/* construct the parser */
	hdp_parser<hol_term> parser = hdp_parser<hol_term>(
			(unsigned int) built_in_predicates::UNKNOWN,
			names, "english.morph.short", "english.gram");

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
	typedef sentence<rooted_syntax_node<flagged_logical_form<hol_term>>> sentence_type;
	array<array_map<sentence_type, flagged_logical_form<hol_term>>> seed_training_set(64);
	while (index < tokens.length) {
		if (!seed_training_set.ensure_capacity(seed_training_set.length + 1)
		 || !array_map_init(seed_training_set[seed_training_set.length], 4))
		{
			for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
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
			for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
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
		for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
			for (auto entry : paragraph) { free(entry.key); free(entry.value); }
			free(paragraph);
		}
		return EXIT_FAILURE;
	}

	/* set the named entities in the seed training set to be "known", so we don't go looking for their definitions later */
	hash_set<unsigned int> seed_entities(64);
	for (const array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
		for (const auto& entry : paragraph) {
			array<string> named_entities(16);
			if (!get_named_entities(*entry.value.root, named_entities)
			 || !seed_entities.check_size(seed_entities.size + named_entities.length))
			{
				for (auto entry : names) free(entry.key);
				for (string& entity : named_entities) free(entity);
				for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
					for (auto entry : paragraph) { free(entry.key); free(entry.value); }
					free(paragraph);
				}
				return EXIT_FAILURE;
			}

			/* if the logical form is a string, add it as a known entity */
			if (entry.value.root->type == hol_term_type::STRING) {
				if (!init(named_entities[named_entities.length], entry.value.root->str)) {
					for (auto entry : names) free(entry.key);
					for (string& entity : named_entities) free(entity);
					for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
						for (auto entry : paragraph) { free(entry.key); free(entry.value); }
						free(paragraph);
					}
					return EXIT_FAILURE;
				}
				named_entities.length++;
			}

			for (const string& named_entity : named_entities) {
				unsigned int id;
				if (!get_token(named_entity, id, names)
				 || !seed_entities.add(id))
				{
					for (auto entry : names) free(entry.key);
					for (string& entity : named_entities) free(entity);
					for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
						for (auto entry : paragraph) { free(entry.key); free(entry.value); }
						free(paragraph);
					}
					return EXIT_FAILURE;
				}
			}
			for (string& entity : named_entities) free(entity);
		}
	}

/*run_console(stdin, "\nEnter sentence to parse: ", parser, names, seed_training_set);
for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
	for (auto entry : paragraph) { free(entry.key); free(entry.value); }
	free(paragraph);
}
for (auto entry : names) free(entry.key);
if (seed_training_set.length > 0)
return EXIT_SUCCESS;*/

	for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
		for (auto entry : paragraph) { free(entry.key); free(entry.value); }
		free(paragraph);
	}

	/* read the articles */
	in = fopen("geoquery.txt", "r");
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
	typedef article<rooted_syntax_node<flagged_logical_form<hol_term>>> article_type;
	typedef in_memory_article_store<rooted_syntax_node<flagged_logical_form<hol_term>>> article_store_type;
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

	if (!parser.invert_name_map(names)) {
		fprintf(stderr, "ERROR: `hdp_parser.invert_name_map` failed.\n");
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

	/* read the articles */
	theory<natural_deduction<hol_term>, polymorphic_canonicalizer<true, false, built_in_predicates>> T(1000000000);
	constant_offset = T.new_constant_offset;
	auto constant_prior = make_simple_constant_distribution(
			chinese_restaurant_process<unsigned int>(1.0), chinese_restaurant_process<unsigned int>(1.0));
	auto theory_element_prior = make_simple_hol_term_distribution(constant_prior, geometric_distribution(0.2), 0.01, 0.45, 0.35, 0.1, 0.09, 0.3, 0.4, 0.2, 0.4);
	auto axiom_prior = make_dirichlet_process(1.0e-3, theory_element_prior);
	auto conjunction_prior = uniform_subset_distribution<const nd_step<hol_term>*>(0.1);
	auto universal_introduction_prior = unif_distribution<unsigned int>();
	auto universal_elimination_prior = chinese_restaurant_process<hol_term>(1.0);
	auto term_indices_prior = make_levy_process(poisson_distribution(1.0), poisson_distribution(1.0));
	auto proof_prior = make_canonicalized_proof_prior(axiom_prior, conjunction_prior,
			universal_introduction_prior, universal_elimination_prior, term_indices_prior, poisson_distribution(5.0));
	decltype(proof_prior)::PriorState proof_axioms;
	if (!parser.invert_name_map(names)) {
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

	/*read_article(names.get("Des Moines"), corpus, parser, T, names, seed_entities, proof_prior);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;*/

	read_sentence(corpus, parser, "Louisiana is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "Arkansas is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	/*read_sentence(corpus, parser, "Oklahoma is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "New Mexico is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);*/
	read_sentence(corpus, parser, "There are 2 states that border Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of Louisiana is 104.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of Arkansas is 103.", T, names, seed_entities, proof_prior, proof_axioms);
	/*read_sentence(corpus, parser, "The area of Oklahoma is 200.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of New Mexico is 1082.", T, names, seed_entities, proof_prior, proof_axioms);*/
	/*read_sentence(corpus, parser, "The largest state bordering Texas is Arkansas.", T, names, seed_entities, proof_prior, proof_axioms);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;*/

	array<string> answers(4);
	/*if (answer_question<true>(answers, "Pittsburgh is in what state?", 10000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "Des Moines is located in what state?", 10000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "The population of Arizona is what?", 10000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	}*/ if (answer_question<true>(answers, "What is the largest state bordering Texas?", 10000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	}
for (string& str : answers) free(str);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;

	array_map<hol_term*, unsigned int> tracked_logical_forms(2);
	/*tracked_logical_forms.put(
		hol_term::new_for_all(1,
			hol_term::new_if_then(
				hol_term::new_atom(parser.symbol_map.get(names.get("cat")), &hol_term::variables<1>::value),
				hol_term::new_atom(parser.symbol_map.get(names.get("mammal")), &hol_term::variables<1>::value)
			)
		), 0);
	hol_term::variables<1>::value.reference_count += 2;
	tracked_logical_forms.put(
		hol_term::new_for_all(1,
			hol_term::new_if_then(
				hol_term::new_atom(parser.symbol_map.get(names.get("mammal")), &hol_term::variables<1>::value),
				hol_term::new_atom(parser.symbol_map.get(names.get("cat")), &hol_term::variables<1>::value)
			)
		), 0);
	hol_term::variables<1>::value.reference_count += 2;*/

	unsigned int iterations = 120000;
	timer stopwatch;
	auto scribe = parser.get_printer();
array_multiset<unsigned int> set_size_distribution(16);
	for (unsigned int t = 0; t < iterations; t++) {
proof_axioms.check_proof_axioms(T);
T.check_disjunction_introductions();
T.sets.are_elements_unique();
T.sets.check_freeable_sets();
T.sets.are_descendants_valid();
T.sets.are_set_sizes_valid();
T.sets.check_set_ids();
		if (stopwatch.milliseconds() > 1000) {
			print("[iteration ", stdout); print(t, stdout); print("]\n", stdout);
			for (const auto& entry : tracked_logical_forms) {
				print("p(", stdout); print(*entry.key, stdout, scribe); print(" axiom) ≈ ", stdout); print((double) entry.value / t, stdout); print('\n', stdout);
			}
for (const auto& entry : set_size_distribution.counts)
fprintf(stderr, "%u %lf\n", entry.key, (double) entry.value / set_size_distribution.sum);
			stopwatch.start();
		}
/*if (t == 21)
fprintf(stderr, "DEBUG: BREAKPOINT\n");*/
		null_collector collector;
		do_mh_step(T, proof_prior, proof_axioms, collector);

		for (auto entry : tracked_logical_forms)
			if (contains_axiom(T, entry.key)) entry.value++;
hol_term* set_formula = hol_term::new_atom(names.get("fish"), &hol_term::variables<1>::value);
hol_term::variables<1>::value.reference_count++;
bool contains;
unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
if (contains) set_size_distribution.add(T.sets.sets[set_id].set_size);
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
	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
