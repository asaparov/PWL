#include "higher_order_logic.h"
#include "hdp_parser.h"
#include "executive.h"
#include "ruletaker.h"

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

struct very_light_tail_distribution {
	typedef unsigned int ObservationType;

	double lambda, log_lambda, log_normalization;

	very_light_tail_distribution(double lambda) : lambda(lambda), log_lambda(log(lambda)), log_normalization(0.0) {
		for (unsigned int k = 0; ; k++) {
			double old_log_normalization = log_normalization;
			log_normalization = logsumexp(log_normalization, (k * k) * log_lambda);
			if (log_normalization == old_log_normalization)
				break;
		}
	}

	very_light_tail_distribution(const very_light_tail_distribution& src) : lambda(src.lambda), log_lambda(src.log_lambda), log_normalization(src.log_normalization) { }
};

inline double log_probability(unsigned int k, const very_light_tail_distribution& prior) {
	return prior.log_normalization + (k * k) * prior.log_lambda;
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

	default_array(const default_array<T>& src) = delete;

	inline bool add(const T& item) {
		return a.add(item);
	}

	inline auto begin() -> decltype(a.begin()) {
		return a.begin();
	}

	inline auto end() -> decltype(a.end()) {
		return a.end();
	}

	inline auto begin() const -> decltype(a.begin()) {
		return a.begin();
	}

	inline auto end() const -> decltype(a.end()) {
		return a.end();
	}

	static inline void move(const default_array<T>& src, default_array<T>& dst) {
		core::move(src.a, dst.a);
	}

	static inline void free(default_array<T>& array) {
		core::free(array.a);
	}
};

template<typename T>
inline bool init(default_array<T>& array) {
	return array_init(array.a, 16);
}

template<typename T>
inline size_t size(const default_array<T>& array) {
	return array.a.length;
}

template<typename T, typename Stream>
inline bool print(const default_array<T>& array, Stream& out) {
	return print(array.a, out);
}

template<typename T>
struct default_array_multiset {
	array_multiset<T, true> a;

	default_array_multiset() : a(16) { }

	default_array_multiset(const default_array_multiset<T>& src) = delete;

	inline bool add(const T& item) {
		return a.add(item);
	}

	static inline void move(const default_array_multiset<T>& src, default_array_multiset<T>& dst) {
		core::move(src.a, dst.a);
	}

	static inline void free(default_array_multiset<T>& multiset) {
		core::free(multiset.a);
	}
};

template<typename T>
inline bool init(default_array_multiset<T>& multiset) {
	return init(multiset.a, 16);
}

template<typename T>
struct default_hash_multiset {
	hash_multiset<T, true> a;

	default_hash_multiset() : a(16) { }

	default_hash_multiset(const default_hash_multiset<T>& src) : a(src.a.counts.table.capacity)
	{
		for (const auto& entry : src.a.counts)
			a.counts.put(entry.key, entry.value);
		a.sum = src.a.sum;
	}

	inline bool add(const T& item) {
		return a.add(item);
	}

	inline bool add(const default_array_multiset<T>& changes) {
		return a.add(changes.a);
	}

	inline void subtract(const default_array_multiset<T>& changes) {
		a.template subtract<true>(changes.a);
	}

	static inline void move(const default_hash_multiset<T>& src, default_hash_multiset<T>& dst) {
		core::move(src.a, dst.a);
	}

	static inline void free(default_hash_multiset<T>& multiset) {
		core::free(multiset.a);
	}
};

template<typename T>
inline bool init(default_hash_multiset<T>& multiset) {
	return init(multiset.a, 16);
}

template<typename T>
inline bool init(default_hash_multiset<T>& multiset, const default_hash_multiset<T>& src) {
	if (!init(multiset.a, src.a.counts.table.capacity))
		return false;
	for (const auto& entry : src.a.counts)
		multiset.a.counts.put(entry.key, entry.value);
	multiset.a.sum = src.a.sum;
	return true;
}

template<typename T>
inline unsigned int size(const default_array_multiset<T>& multiset) {
	return multiset.a.counts.size;
}

template<typename T>
inline unsigned int size(const default_hash_multiset<T>& multiset) {
	return multiset.a.counts.table.size;
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

struct empty_prior_state {
	template<typename T>
	inline constexpr bool add(const T& changes) const { return true; }

	template<typename T>
	inline void subtract(const T& changes) const { }

	static inline void move(const empty_prior_state& src, const empty_prior_state& dst) { }
	static inline void free(const empty_prior_state& state) { }
};

constexpr bool init(const empty_prior_state& new_state) { return true; }
constexpr unsigned int size(const empty_prior_state& new_state) { return 0; }

template<typename T>
struct iid_uniform_distribution
{
	typedef empty_prior_state PriorState;
	typedef default_array<T> PriorStateChanges;
	typedef T ObservationType;
	typedef default_array<T> ObservationCollection;

	double log_n;

	iid_uniform_distribution(unsigned int n) : log_n(log(n)) { }
	iid_uniform_distribution(const iid_uniform_distribution<T>& src) : log_n(src.log_n) { }
	~iid_uniform_distribution() { }

private:
	iid_uniform_distribution() { }
	template<typename A> friend iid_uniform_distribution<A> make_iid_uniform_distribution_from_log_n(double);
};

template<typename T>
inline iid_uniform_distribution<T> make_iid_uniform_distribution(unsigned int n) {
	return iid_uniform_distribution<T>(n);
}

template<typename T>
inline iid_uniform_distribution<T> make_iid_uniform_distribution_from_log_n(double log_n) {
	iid_uniform_distribution<T> dist;
	dist.log_n = log_n;
	return dist;
}

template<typename T, typename ObservationCollection>
inline double log_probability(
		const ObservationCollection& observations,
		const iid_uniform_distribution<T>& prior)
{
	double value = -(size(observations) * prior.log_n);
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of `iid_uniform_distribution`: %lf.\n", value);
	print("  Observations: ", stderr); print(observations, stderr); print('\n', stderr);
#endif
	return value;
}

template<typename T, typename ObservationCollection>
inline double log_probability_ratio(
		const empty_prior_state& dummy,
		const ObservationCollection& old_observations,
		const ObservationCollection& new_observations,
		const iid_uniform_distribution<T>& prior)
{
	return (size(old_observations) * prior.log_n) - (size(new_observations) * prior.log_n);
}

template<typename T, typename ObservationCollection>
inline double log_probability_ratio(
		const empty_prior_state& dummy,
		const ObservationCollection& old_observations,
		const ObservationCollection& new_observations,
		const iid_uniform_distribution<T>& prior,
		const default_array<T>& old_prior_changes,
		const default_array<T>& new_prior_changes)
{
	return (size(old_observations) * prior.log_n) - (size(new_observations) * prior.log_n);
}

template<typename T>
struct chinese_restaurant_process
{
	typedef default_array_multiset<T> ObservationCollection;
	typedef default_hash_multiset<T> PriorState;
	typedef default_array_multiset<T> PriorStateChanges;

	double alpha, log_alpha;
	double d, log_d, lgamma_one_minus_d;

	chinese_restaurant_process(double alpha, double d) :
		alpha(alpha), log_alpha(log(alpha)), d(d), log_d(log(d)), lgamma_one_minus_d(lgamma(1.0 - d))
	{ }

	chinese_restaurant_process(const chinese_restaurant_process<T>& src) :
		alpha(src.alpha), log_alpha(log(alpha)), d(src.d), log_d(src.log_d), lgamma_one_minus_d(src.lgamma_one_minus_d)
	{ }

	~chinese_restaurant_process() { }
};

template<typename MultisetType, typename T>
double log_probability(
		const MultisetType& observations,
		const chinese_restaurant_process<T>& prior)
{
	double value = 0.0;
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of `chinese_restaurant_process`:\n");
#endif
	for (unsigned int i = 0; i < size(observations); i++) {
		double current_value = lgamma(get_value(observations, i) - prior.d) - prior.lgamma_one_minus_d;
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "  Observation "); print(get_key(observations, i), stderr);
		fprintf(stderr, " has count %u and log probability: %lf.\n", get_value(observations, i), current_value);
#endif
		value += current_value;
	}
	double normalization = lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
	if (prior.d == 0.0)
		normalization += prior.log_alpha * size(observations);
	else normalization += prior.log_d * size(observations) + lgamma(prior.alpha/prior.d + size(observations)) - lgamma(prior.alpha/prior.d);
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Prior: %lf.\n", normalization);
#endif
	value += normalization;
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total: %lf.\n", value);
#endif
	return value;
}

template<typename Stream>
inline bool print(const hol_term* term, Stream& stream) {
	return print(*term, stream);
}

template<typename MultisetType,
	template<typename> class Collection, typename T>
double log_probability(
		const MultisetType& observations,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& clusters)
{
	double value = 0.0;
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of `chinese_restaurant_process`:\n");
#endif
	for (unsigned int i = 0; i < size(observations); i++) {
		clusters.add(get_key(observations, i));
		double current_value = lgamma(get_value(observations, i) - prior.d) - prior.lgamma_one_minus_d;
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "  Observation "); print(get_key(observations, i), stderr);
		fprintf(stderr, " has count %u and log probability: %lf.\n", get_value(observations, i), current_value);
#endif
		value += current_value;
	}
	double normalization = lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
	if (prior.d == 0.0)
		normalization += prior.log_alpha * size(observations);
	else normalization += prior.log_d * size(observations) + lgamma(prior.alpha/prior.d + size(observations)) - lgamma(prior.alpha/prior.d);
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Prior: %lf.\n", normalization);
#endif
	value += normalization;
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total: %lf.\n", value);
#endif
	return value;
}

template<bool AutomaticallyFree, typename T, template<typename> class Collection>
inline double log_probability_ratio_old_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& old_clusters,
		unsigned int& old_cluster_count)
{
#if !defined(NDEBUG)
	if (!tables.counts.table.contains(observation))
		fprintf(stderr, "log_probability_ratio_old_cluster WARNING: This chinese_restaurant_process does not contain the given observation.\n");
#endif
	unsigned int count = tables.counts.get(observation);
	if (count == frequency) {
		old_clusters.add(observation);
		old_cluster_count++;
		return -lgamma(count - prior.d) - (prior.d == 0 ? prior.log_alpha : prior.log_d) + prior.lgamma_one_minus_d;
	} else {
		return lgamma(count - frequency - prior.d) - lgamma(count - prior.d);
	}
}

template<bool AutomaticallyFree, typename T, template<typename> class Collection>
inline double log_probability_ratio_new_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection<T>& new_clusters,
		unsigned int& new_cluster_count)
{
	bool contains;
	unsigned int count = tables.counts.get(observation, contains);
#if !defined(NDEBUG)
	if (contains && count == 0)
		fprintf(stderr, "log_probability_ratio_new_cluster WARNING: The hash_multiset has an observation with zero count.\n");
#endif
	if (contains) {
		return lgamma(count + frequency - prior.d) - lgamma(count - prior.d);
	} else {
		new_clusters.add(observation);
		new_cluster_count++;
		return (prior.d == 0 ? prior.log_alpha : prior.log_d) + lgamma(frequency - prior.d) - prior.lgamma_one_minus_d;
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
	unsigned int old_cluster_count = 0, new_cluster_count = 0;
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
				old_cluster_count++;
				value += -lgamma(count - prior.d) - (prior.d == 0 ? prior.log_alpha : prior.log_d) + prior.lgamma_one_minus_d;
			} else {
				value += lgamma(count + diff - prior.d) - lgamma(count - prior.d);
			}
			i++; j++;
		} else if (get_key(old_observations, i) < get_key(new_observations, j)) {
			value += log_probability_ratio_old_cluster(tables,
					get_key(old_observations, i), get_value(old_observations, i), prior, old_clusters, old_cluster_count);
			i++;
		} else {
			value += log_probability_ratio_new_cluster(tables,
					get_key(new_observations, j), get_value(new_observations, j), prior, new_clusters, new_cluster_count);
			j++;
		}
	}

	while (i < size(old_observations)) {
		value += log_probability_ratio_old_cluster(tables,
				get_key(old_observations, i), get_value(old_observations, i), prior, old_clusters, old_cluster_count);
		i++;
	} while (j < size(new_observations)) {
		value += log_probability_ratio_new_cluster(tables,
				get_key(new_observations, j), get_value(new_observations, j), prior, new_clusters, new_cluster_count);
		j++;
	}

	value += lgamma(prior.alpha + tables.sum) - lgamma(prior.alpha + (tables.sum - sum(old_observations) + sum(new_observations)));
	if (prior.d != 0.0)
		value += -lgamma(prior.alpha/prior.d + tables.counts.table.size) + lgamma(prior.alpha/prior.d + (tables.counts.table.size - old_cluster_count + new_cluster_count));
	return value;
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

		prior_state() { }
		prior_state(const prior_state& src) : base_prior_state(src.base_prior_state) { }

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
		restaurant(alpha, 0.0), base_distribution(base_distribution)
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

template<typename T, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
inline void free_all(const default_array<T>& a) { }

template<typename T, typename std::enable_if<!std::is_pointer<T>::value>::type* = nullptr>
inline void free_all(default_array<T>& a) {
	for (T& element : a)
		free(element);
}

template<typename MultisetType, typename BaseDistribution>
double log_probability(
		const MultisetType& restaurant_tables,
		const dirichlet_process<BaseDistribution>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters clusters;
	double value = log_probability(restaurant_tables, prior.restaurant, clusters);
	value += log_probability(clusters, prior.base_distribution);
	free_all(clusters);
	return value;
}

template<typename MultisetType, typename BaseDistribution>
double log_probability(
		const MultisetType& restaurant_tables,
		const array<typename BaseDistribution::ObservationType>& extra_observations,
		const dirichlet_process<BaseDistribution>& prior)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters clusters;
	double value = log_probability(restaurant_tables, prior.restaurant, clusters);
	for (const typename BaseDistribution::ObservationType& extra_observation : extra_observations) {
		bool contains = false;
		for (const typename BaseDistribution::ObservationType& existing_observation : clusters) {
			if (existing_observation == extra_observation || *existing_observation == *extra_observation) {
				contains = true;
				break;
			}
		}
		if (!contains)
			clusters.add(extra_observation);
	}
	value += log_probability(clusters, prior.base_distribution);
	free_all(clusters);
	return value;
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
	double value = log_probability_ratio(restaurant_tables, old_observations, new_observations, prior.restaurant, old_clusters, new_clusters);
	value += log_probability_ratio(prior_state.base_prior_state, old_clusters, new_clusters, prior.base_distribution, old_prior_changes, new_prior_changes);
	free_all(old_clusters);
	free_all(new_clusters);
	return value;
}

template<typename MultisetType, typename T, bool AutomaticallyFree, typename BaseDistribution>
double log_probability_ratio(
		const MultisetType& restaurant_tables,
		const array_multiset<T, AutomaticallyFree>& old_observations,
		const array_multiset<T, AutomaticallyFree>& new_observations,
		const array<T>& old_extra_observations,
		const array<T>& new_extra_observations,
		const dirichlet_process<BaseDistribution>& prior,
		const typename dirichlet_process<BaseDistribution>::prior_state& prior_state,
		typename BaseDistribution::PriorStateChanges& old_prior_changes,
		typename BaseDistribution::PriorStateChanges& new_prior_changes)
{
	typedef typename BaseDistribution::ObservationCollection Clusters;

	Clusters old_clusters, new_clusters;
	double value = log_probability_ratio(restaurant_tables, old_observations, new_observations, prior.restaurant, old_clusters, new_clusters);
	for (const T& extra_observation : old_extra_observations) {
		bool contains = false;
		for (const T& existing_observation : old_clusters) {
			if (existing_observation == extra_observation || *existing_observation == *extra_observation) {
				contains = true;
				break;
			}
		}
		if (!contains)
			old_clusters.add(extra_observation);
	} for (const T& extra_observation : new_extra_observations) {
		bool contains = false;
		for (const T& existing_observation : new_clusters) {
			if (existing_observation == extra_observation || *existing_observation == *extra_observation) {
				contains = true;
				break;
			}
		}
		if (!contains)
			new_clusters.add(extra_observation);
	}
	value += log_probability_ratio(prior_state.base_prior_state, old_clusters, new_clusters, prior.base_distribution, old_prior_changes, new_prior_changes);
	free_all(old_clusters);
	free_all(new_clusters);
	return value;
}

template<typename ConstantDistribution, typename PredicateDistribution, typename TypeDistribution>
struct simple_constant_distribution
{
	struct prior_state_changes {
		typename ConstantDistribution::PriorStateChanges constants;
		typename PredicateDistribution::PriorStateChanges predicates;
		array_map<unsigned int, array_multiset<hol_term>> types;

		prior_state_changes() : types(8) { }
		~prior_state_changes() {
			for (auto entry : types) core::free(entry.value);
		}

		inline bool add_constant(unsigned int constant) {
			return constants.add(constant);
		}

		inline bool add_predicate(unsigned int predicate) {
			return predicates.add(predicate);
		}

		inline bool add_type(const hol_term* type, unsigned int constant) {
			if (!types.ensure_capacity(types.size + 1))
				return false;
			unsigned int index = types.index_of(constant);
			if (index == types.size) {
				if (!init(types.values[index], 4)) {
					return false;
				} else if (!types.values[index].add(*type)) {
					free(types.values[index]);
					return false;
				}
				types.keys[index] = constant;
				types.size++;
				insertion_sort(types.keys, types.values, types.size, default_sorter());
				return true;
			} else {
				return types.values[index].add(*type);
			}
		}
	};

	struct prior_state {
		typename ConstantDistribution::PriorState constants;
		typename PredicateDistribution::PriorState predicates;
		hash_map<unsigned int, hash_multiset<hol_term>> types;

		prior_state() : types(16) { }

		prior_state(const prior_state& src) : constants(src.constants), predicates(src.predicates), types(src.types.table.capacity) {
			for (const auto& entry : src.types) {
				bool contains; unsigned int bucket;
				types.get(entry.key, contains, bucket);
				if (!init(types.values[bucket], entry.value.counts.table.capacity))
					exit(EXIT_FAILURE);
				types.values[bucket].counts.put_all(entry.value.counts);
				types.values[bucket].sum = entry.value.sum;
				types.table.keys[bucket] = entry.key;
				types.table.size++;
			}
		}

		~prior_state() {
			for (auto entry : types)
				free(entry.value);
		}

		inline bool add(const prior_state_changes& changes) {
			if (!constants.add(changes.constants)) {
				return false;
			} else if (!predicates.add(changes.predicates)) {
				constants.subtract(changes.constants);
				return false;
			}

			if (!types.check_size(types.table.size + changes.types.size)) {
				constants.subtract(changes.constants);
				predicates.subtract(changes.predicates);
				return false;
			}
			for (unsigned int i = 0; i < changes.types.size; i++) {
				bool contains; unsigned int index;
				hash_multiset<hol_term>& value = types.get(changes.types.keys[i], contains, index);
				if (!contains) {
					if (!init(types.values[index], 4)) {
						constants.subtract(changes.constants);
						predicates.subtract(changes.predicates);
						for (unsigned int j = 0; j < i; j++)
							types.get(changes.types.keys[j]).subtract(changes.types.values[j]);
						return false;
					}
					types.table.keys[index] = changes.types.keys[i];
					types.table.size++;
				}
				if (!value.add(changes.types.values[i])) {
					constants.subtract(changes.constants);
					predicates.subtract(changes.predicates);
					for (unsigned int j = 0; j < i; j++)
						types.get(changes.types.keys[j]).subtract(changes.types.values[j]);
					return false;
				}
			}
			return true;
		}

		inline bool add(const prior_state_changes& changes,
				const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& prior)
		{
			return add(changes);
		}

		inline void subtract(const prior_state_changes& changes) {
			constants.subtract(changes.constants);
			predicates.subtract(changes.predicates);
			for (unsigned int i = 0; i < changes.types.size; i++) {
				bool contains; unsigned int index;
				hash_multiset<hol_term>& value = types.get(changes.types.keys[i], contains, index);
#if !defined(NDEBUG)
				if (!contains)
					fprintf(stderr, "simple_constant_distribution.prior_state.subtract WARNING: The given key does not exist in `types`.\n");
#endif
				value.subtract<true>(changes.types.values[i]);
				if (value.sum == 0) {
					free(value);
					types.remove_at(index);
				}
			}
		}

		inline void subtract(const prior_state_changes& changes,
				const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& prior)
		{
			subtract(changes);
		}
	};

	typedef prior_state PriorState;
	typedef prior_state_changes PriorStateChanges;
	typedef prior_state_changes ObservationCollection;

	ConstantDistribution constant_distribution;
	PredicateDistribution predicate_distribution;
	TypeDistribution type_distribution;

	simple_constant_distribution(
			const ConstantDistribution& constant_distribution,
			const PredicateDistribution& predicate_distribution,
			const TypeDistribution& type_distribution) :
		constant_distribution(constant_distribution),
		predicate_distribution(predicate_distribution),
		type_distribution(type_distribution)
	{ }

	simple_constant_distribution(const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& src) :
		constant_distribution(src.constant_distribution), predicate_distribution(src.predicate_distribution), type_distribution(src.type_distribution)
	{ }
};

template<typename ConstantDistribution, typename PredicateDistribution, typename TypeDistribution>
inline simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>
make_simple_constant_distribution(
		const ConstantDistribution& constant_distribution,
		const PredicateDistribution& predicate_distribution,
		const TypeDistribution& type_distribution)
{
	return simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>(constant_distribution, predicate_distribution, type_distribution);
}

template<typename ConstantCollection, typename ConstantDistribution, typename PredicateDistribution, typename TypeDistribution>
inline double log_probability(const ConstantCollection& constants,
		const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& prior)
{
	double value = log_probability(constants.constants, prior.constant_distribution)
				 + log_probability(constants.predicates, prior.predicate_distribution);
	for (const auto& entry : constants.types) {
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "log_probability of type distribution for constant %u:\n", entry.key);
#endif
		value += log_probability(entry.value, prior.type_distribution);
	}
	return value;
}

template<typename ObservationCollection, typename ConstantCollection,
	typename ConstantDistribution, typename PredicateDistribution, typename TypeDistribution>
inline double log_probability_ratio(const ObservationCollection& existing_constants,
		const ConstantCollection& old_constants, const ConstantCollection& new_constants,
		const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& prior)
{
	double value = log_probability_ratio(existing_constants.constants, old_constants.constants, new_constants.constants, prior.constant_distribution)
				 + log_probability_ratio(existing_constants.predicates, old_constants.predicates, new_constants.predicates, prior.predicate_distribution);

	static typename TypeDistribution::PriorState type_prior_state;
	static hash_multiset<hol_term> empty_type_prior_state(1);
	static array_multiset<hol_term> empty_type_prior_state_changes(1);

	unsigned int i = 0, j = 0; bool contains;
	default_array<hol_term> old_prior_changes, new_prior_changes;
	while (i < old_constants.types.size && j < new_constants.types.size) {
		if (old_constants.types.keys[i] == new_constants.types.keys[j]) {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
			if (contains) {
				value += log_probability_ratio(existing, old_constants.types.values[i], new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			} else {
				value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			}
			i++; j++;
		} else if (old_constants.types.keys[i] < new_constants.types.keys[j]) {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
			if (contains) {
				value += log_probability_ratio(existing, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			} else {
				value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			}
			i++;
		} else {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(new_constants.types.keys[j], contains);
			if (contains) {
				value += log_probability_ratio(existing, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			} else {
				value += log_probability_ratio(empty_type_prior_state, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
			}
			j++;
		}
	} while (i < old_constants.types.size) {
		const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
		if (contains) {
			value += log_probability_ratio(existing, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
		} else {
			value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
		}
		i++;
	} while (j < new_constants.types.size) {
		const hash_multiset<hol_term>& existing = existing_constants.types.get(new_constants.types.keys[j], contains);
		if (contains) {
			value += log_probability_ratio(existing, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
		} else {
			value += log_probability_ratio(empty_type_prior_state, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution, type_prior_state, old_prior_changes, new_prior_changes);
		}
		j++;
	}
	return value;
}

template<typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
struct simple_hol_term_distribution
{
	struct prior_state_changes {
		typename ConstantDistribution::PriorStateChanges constants;
		array_map<unsigned int, unsigned int> arg1_map;
		array_map<unsigned int, string*> arg2_string_map;

		prior_state_changes() : arg1_map(4), arg2_string_map(4) { }
	};

	struct prior_state {
		typename ConstantDistribution::PriorState constant_prior_state;
		array_map<unsigned int, unsigned int> arg1_map;
		array_map<unsigned int, string*> arg2_string_map;

		prior_state() : arg1_map(16), arg2_string_map(16) { }

		prior_state(const prior_state& src) : constant_prior_state(src.constant_prior_state), arg1_map(src.arg1_map.capacity), arg2_string_map(src.arg2_string_map.capacity) {
			for (unsigned int i = 0; i < src.arg1_map.size; i++) {
				arg1_map.keys[i] = src.arg1_map.keys[i];
				arg1_map.values[i] = src.arg1_map.values[i];
				arg1_map.size++;
			} for (unsigned int i = 0; i < src.arg2_string_map.size; i++) {
				arg2_string_map.keys[i] = src.arg2_string_map.keys[i];
				arg2_string_map.values[i] = src.arg2_string_map.values[i];
				arg2_string_map.size++;
			}
		}

		inline bool add(const prior_state_changes& changes) {
			for (const auto& entry : changes.arg1_map) {
#if !defined(NDEBUG)
				if (arg1_map.contains(entry.key))
					fprintf(stderr, "simple_hol_term_distribution.prior_state.add ERROR: `arg1_map` already contains the key %u.\n", entry.key);
#endif
				if (!arg1_map.put(entry.key, entry.value))
					return false;
			} for (const auto& entry : changes.arg2_string_map) {
#if !defined(NDEBUG)
				if (arg2_string_map.contains(entry.key))
					fprintf(stderr, "simple_hol_term_distribution.prior_state.add ERROR: `arg2_string_map` already contains the key %u.\n", entry.key);
#endif
				if (!arg2_string_map.put(entry.key, entry.value))
					return false;
			}
			return constant_prior_state.add(changes.constants);
		}

		inline bool add(const hol_term* observation,
				const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior)
		{
			prior_state_changes changes;
			log_probability_helper(observation, prior, changes);
			return add(changes);
		}

		inline void subtract(const prior_state_changes& changes) {
			constant_prior_state.subtract(changes.constants);
			for (const auto& entry : changes.arg1_map) {
				unsigned int index = arg1_map.index_of(entry.key);
#if !defined(NDEBUG)
				if (index == arg1_map.size)
					fprintf(stderr, "simple_hol_term_distribution.prior_state.subtract ERROR: `arg1_map` doesn't contain the key %u.\n", entry.key);
#endif
				arg1_map.remove_at(index);
			} for (const auto& entry : changes.arg2_string_map) {
				unsigned int index = arg2_string_map.index_of(entry.key);
#if !defined(NDEBUG)
				if (index == arg2_string_map.size)
					fprintf(stderr, "simple_hol_term_distribution.prior_state.subtract ERROR: `arg2_string_map` doesn't contain the key %u.\n", entry.key);
#endif
				arg2_string_map.remove_at(index);
			}
		}

		inline void subtract(const hol_term* observation,
				const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior)
		{
			prior_state_changes changes;
			log_probability_helper(observation, prior, changes);
			return subtract(changes);
		}
	};

	typedef hol_term* ObservationType;
	typedef default_array<hol_term*> ObservationCollection;
	typedef prior_state PriorState;
	typedef prior_state_changes PriorStateChanges;

	double log_ground_literal_probability;
	double log_negated_expression_probability;
	double log_universal_probability;
	double log_existential_probability;
	double log_conjunction_probability;
	double log_disjunction_probability;
	double log_set_size_axiom_probability;
	double log_arg1_probability;
	double log_arg2_probability;

	double log_quantified_ground_literal_probability;
	double log_quantified_negated_expression_probability;
	double log_quantified_universal_probability;
	double log_quantified_existential_probability;
	double log_quantified_conjunction_probability;
	double log_quantified_disjunction_probability;
	double log_quantified_set_size_axiom_probability;
	double log_quantified_arg1_probability;
	double log_quantified_arg2_probability;

	double log_arg_constant_probability;
	double log_arg_number_probability;
	double log_arg_string_probability;

	double log_negation_probability;
	double log_positive_probability;

	double log_unary_probability;
	double log_binary_probability;

	double log_antecedent_continue_probability;
	double log_antecedent_stop_probability;
	double log_consequent_continue_probability;
	double log_consequent_stop_probability;

	very_light_tail_distribution name_count_distribution;

	ConstantDistribution constant_distribution;
	SetSizeDistribution set_size_distribution;

	simple_hol_term_distribution(
			const ConstantDistribution& constant_distribution,
			const SetSizeDistribution& set_size_distribution,
			double ground_literal_probability, double negated_expression_probability,
			double universal_probability, double existential_probability,
			double conjunction_probability, double disjunction_probability,
			double set_size_axiom_probability, double arg1_probability,
			double arg2_probability,
			double quantified_ground_literal_probability, double quantified_negated_expression_probability,
			double quantified_universal_probability, double quantified_existential_probability,
			double quantified_conjunction_probability, double quantified_disjunction_probability,
			double quantified_set_size_axiom_probability, double quantified_arg1_probability,
			double quantified_arg2_probability,
			double arg_constant_probability, double arg_number_probability,
			double arg_string_probability, double negation_probability,
			double unary_probability, double antecedent_stop_probability,
			double consequent_stop_probability, double name_count_parameter) :
		log_ground_literal_probability(log(ground_literal_probability)),
		log_negated_expression_probability(log(negated_expression_probability)),
		log_universal_probability(log(universal_probability)),
		log_existential_probability(log(existential_probability)),
		log_conjunction_probability(log(conjunction_probability)),
		log_disjunction_probability(log(disjunction_probability)),
		log_set_size_axiom_probability(log(set_size_axiom_probability)),
		log_arg1_probability(log(arg1_probability)),
		log_arg2_probability(log(arg2_probability)),
		log_quantified_ground_literal_probability(log(quantified_ground_literal_probability)),
		log_quantified_negated_expression_probability(log(quantified_negated_expression_probability)),
		log_quantified_universal_probability(log(quantified_universal_probability)),
		log_quantified_existential_probability(log(quantified_existential_probability)),
		log_quantified_conjunction_probability(log(quantified_conjunction_probability)),
		log_quantified_disjunction_probability(log(quantified_disjunction_probability)),
		log_quantified_set_size_axiom_probability(log(quantified_set_size_axiom_probability)),
		log_quantified_arg1_probability(log(quantified_arg1_probability)),
		log_quantified_arg2_probability(log(quantified_arg2_probability)),
		log_arg_constant_probability(log(arg_constant_probability)),
		log_arg_number_probability(log(arg_number_probability)),
		log_arg_string_probability(log(arg_string_probability)),
		log_negation_probability(log(negation_probability)),
		log_positive_probability(log(1.0 - negation_probability)),
		log_unary_probability(log(unary_probability)),
		log_binary_probability(log(1.0 - unary_probability)),
		log_antecedent_continue_probability(log(1.0 - antecedent_stop_probability)),
		log_antecedent_stop_probability(log(antecedent_stop_probability)),
		log_consequent_continue_probability(log(1.0 - consequent_stop_probability)),
		log_consequent_stop_probability(log(consequent_stop_probability)),
		name_count_distribution(name_count_parameter),
		constant_distribution(constant_distribution),
		set_size_distribution(set_size_distribution)
	{
		if (fabs(ground_literal_probability + negated_expression_probability + universal_probability + existential_probability + conjunction_probability + disjunction_probability + set_size_axiom_probability + arg1_probability + arg2_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `ground_literal_probability + negated_expression_probability + universal_probability + existential_probability + conjunction_probability + disjunction_probability + set_size_axiom_probability + arg1_probability + arg2_probability` is not 1.\n");
		if (fabs(quantified_ground_literal_probability + quantified_negated_expression_probability + quantified_universal_probability + quantified_existential_probability + quantified_conjunction_probability + quantified_disjunction_probability + quantified_set_size_axiom_probability + quantified_arg1_probability + quantified_arg2_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `quantified_ground_literal_probability + quantified_negated_expression_probability + quantified_universal_probability + quantified_existential_probability + quantified_conjunction_probability + quantified_disjunction_probability + quantified_set_size_axiom_probability + quantified_arg1_probability + quantified_arg2_probability` is not 1.\n");
		if (fabs(arg_constant_probability + arg_number_probability + arg_string_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `arg_constant_probability + arg_number_probability + arg_string_probability` is not 1.\n");
	}

	simple_hol_term_distribution(const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& src) :
		log_ground_literal_probability(src.log_ground_literal_probability),
		log_negated_expression_probability(src.log_negated_expression_probability),
		log_universal_probability(src.log_universal_probability),
		log_existential_probability(src.log_existential_probability),
		log_conjunction_probability(src.log_conjunction_probability),
		log_disjunction_probability(src.log_disjunction_probability),
		log_set_size_axiom_probability(src.log_set_size_axiom_probability),
		log_arg1_probability(src.log_arg1_probability),
		log_arg2_probability(src.log_arg2_probability),
		log_quantified_ground_literal_probability(src.log_quantified_ground_literal_probability),
		log_quantified_negated_expression_probability(src.log_quantified_negated_expression_probability),
		log_quantified_universal_probability(src.log_quantified_universal_probability),
		log_quantified_existential_probability(src.log_quantified_existential_probability),
		log_quantified_conjunction_probability(src.log_quantified_conjunction_probability),
		log_quantified_disjunction_probability(src.log_quantified_disjunction_probability),
		log_quantified_set_size_axiom_probability(src.log_quantified_set_size_axiom_probability),
		log_quantified_arg1_probability(src.log_quantified_arg1_probability),
		log_quantified_arg2_probability(src.log_quantified_arg2_probability),
		log_arg_constant_probability(src.log_arg_constant_probability),
		log_arg_number_probability(src.log_arg_number_probability),
		log_arg_string_probability(src.log_arg_string_probability),
		log_negation_probability(src.log_negation_probability),
		log_positive_probability(src.log_positive_probability),
		log_unary_probability(src.log_unary_probability),
		log_binary_probability(src.log_binary_probability),
		log_antecedent_continue_probability(src.log_antecedent_continue_probability),
		log_antecedent_stop_probability(src.log_antecedent_stop_probability),
		log_consequent_continue_probability(src.log_consequent_continue_probability),
		log_consequent_stop_probability(src.log_consequent_stop_probability),
		name_count_distribution(src.name_count_distribution),
		constant_distribution(src.constant_distribution),
		set_size_distribution(src.set_size_distribution)
	{ }
};

template<typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
inline simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution> make_simple_hol_term_distribution(
		const ConstantDistribution& constant_distribution,
		const SetSizeDistribution set_size_distribution,
		double ground_literal_probability, double negated_expression_probability,
		double universal_probability, double existential_probability,
		double conjunction_probability, double disjunction_probability,
		double set_size_axiom_probability, double arg1_probability,
		double arg2_probability,
		double quantified_ground_literal_probability, double quantified_negated_expression_probability,
		double quantified_universal_probability, double quantified_existential_probability,
		double quantified_conjunction_probability, double quantified_disjunction_probability,
		double quantified_set_size_axiom_probability, double quantified_arg1_probability,
		double quantified_arg2_probability,
		double arg_constant_probability, double arg_number_probability,
		double arg_string_probability, double negation_probability,
		double unary_probability, double antecedent_stop_probability,
		double consequent_stop_probability, double name_count_parameter)
{
	return simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>(
			constant_distribution, set_size_distribution,
			ground_literal_probability, negated_expression_probability,
			universal_probability, existential_probability,
			conjunction_probability, disjunction_probability,
			set_size_axiom_probability, arg1_probability,
			arg2_probability,
			quantified_ground_literal_probability, quantified_negated_expression_probability,
			quantified_universal_probability, quantified_existential_probability,
			quantified_conjunction_probability, quantified_disjunction_probability,
			quantified_set_size_axiom_probability, quantified_arg1_probability,
			quantified_arg2_probability,
			arg_constant_probability, arg_number_probability,
			arg_string_probability, negation_probability,
			unary_probability, antecedent_stop_probability,
			consequent_stop_probability, name_count_parameter);
}

template<bool Quantified, bool IsRoot, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_atom(const hol_term* function, const hol_term* arg1,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type == hol_term_type::UNARY_APPLICATION) {
		double value = 0.0;
		if (Quantified && arg1->type == hol_term_type::VARIABLE) {
			value = prior.log_binary_probability;
		} else if (!Quantified && arg1->type == hol_term_type::CONSTANT) {
			constants.add_constant(arg1->constant);
			value = prior.log_binary_probability;
		} else {
			return -std::numeric_limits<double>::infinity();
		}

		/* TODO: make sure this is a proper prior (for higher arity
		   applications, the arguments can be a mix of constants and
		   variables) */
		const hol_term* left = function;
		while (left->type == hol_term_type::UNARY_APPLICATION) {
			if (Quantified && left->binary.right->type == hol_term_type::VARIABLE) {
				value += prior.log_binary_probability;
			} else if (left->binary.right->type == hol_term_type::CONSTANT) {
				constants.add_predicate(left->binary.right->constant);
				value += prior.log_binary_probability;
			} else {
				return -std::numeric_limits<double>::infinity();
			}
			left = left->binary.left;
		}
		if (left->type == hol_term_type::CONSTANT) {
			constants.add_predicate(left->constant);
		} else if (!Quantified) {
			return -std::numeric_limits<double>::infinity();
		}
		return value;
	}

	if (function->type == hol_term_type::CONSTANT) {
		constants.add_predicate(function->constant);
	} else if (!Quantified) {
		return -std::numeric_limits<double>::infinity();
	}
	if (Quantified && arg1->type == hol_term_type::VARIABLE) {
		return prior.log_unary_probability;
	} else if (!Quantified && arg1->type == hol_term_type::CONSTANT) {
		if (IsRoot)
			constants.add_type(function, arg1->constant);
		else constants.add_constant(arg1->constant);
		return prior.log_unary_probability;
	} else {
		return -std::numeric_limits<double>::infinity();
	}
}

template<bool Quantified, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_atom(
		const hol_term* function, const hol_term* arg1, const hol_term* arg2,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior,
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

template<bool Quantified, bool IsRoot, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_literal(const hol_term* literal,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (literal->type == hol_term_type::UNARY_APPLICATION) {
		return prior.log_positive_probability + log_probability_atom<Quantified, IsRoot>(literal->binary.left, literal->binary.right, prior, constants);
	} else if (literal->type == hol_term_type::BINARY_APPLICATION) {
		return prior.log_positive_probability + log_probability_atom<Quantified>(literal->ternary.first, literal->ternary.second, literal->ternary.third, prior, constants);
	} else if (literal->type == hol_term_type::NOT && literal->unary.operand->type == hol_term_type::UNARY_APPLICATION) {
		return prior.log_negation_probability + log_probability_atom<Quantified, false>(literal->unary.operand->binary.left, literal->unary.operand->binary.right, prior, constants);
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

template<bool Quantified = false, bool IsRoot = true, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution>
double log_probability_helper(const hol_term* term,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>::prior_state_changes& changes)
{
	if (is_literal(term))
		return (Quantified ? prior.log_quantified_ground_literal_probability : prior.log_ground_literal_probability) + log_probability_literal<Quantified, IsRoot>(term, prior, changes.constants);

	double value;
	const hol_term* antecedent;
	const hol_term* consequent;
	switch (term->type) {
	case hol_term_type::FOR_ALL:
		value = (Quantified ? prior.log_quantified_universal_probability : prior.log_universal_probability);
		if (term->quantifier.operand->type == hol_term_type::IF_THEN) {
			antecedent = term->quantifier.operand->binary.left;
			consequent = term->quantifier.operand->binary.right;
			if (antecedent->type == hol_term_type::AND) {
				value += (antecedent->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < antecedent->array.length; i++)
					value += log_probability_helper<true, false>(antecedent->array.operands[i], prior, changes);
			} else {
				value += prior.log_antecedent_stop_probability + log_probability_helper<true, false>(antecedent, prior, changes);
			}

			if (consequent->type == hol_term_type::AND) {
				value += (consequent->array.length - 1) * prior.log_consequent_continue_probability + prior.log_consequent_stop_probability;
				for (unsigned int i = 0; i < consequent->array.length; i++)
					value += log_probability_helper<true, false>(consequent->array.operands[i], prior, changes);
			} else {
				value += prior.log_consequent_stop_probability + log_probability_helper<true, false>(consequent, prior, changes);
			}
			return value;
		} else {
			// TODO: implement this (this is currently only possible if `term` is a prespecified theorem)
			return 0.0;
		}
	case hol_term_type::EQUALS:
		if (term->binary.left->type == hol_term_type::UNARY_APPLICATION
		 && term->binary.left->binary.left->type == hol_term_type::CONSTANT
		 && term->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
		 && term->binary.left->binary.right->type == hol_term_type::LAMBDA
		 && term->binary.right->type == hol_term_type::NUMBER)
		{
			value = (Quantified ? prior.log_quantified_set_size_axiom_probability : prior.log_set_size_axiom_probability);
			hol_term* operand = term->binary.left->binary.right->quantifier.operand;
			while (operand->type == hol_term_type::LAMBDA)
				operand = operand->quantifier.operand;
			if (operand->type == hol_term_type::AND) {
				value += (operand->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
				for (unsigned int i = 0; i < operand->array.length; i++)
					value += log_probability_helper<true, false>(operand->array.operands[i], prior, changes);
			} else {
				value += prior.log_antecedent_stop_probability + log_probability_helper<true, false>(operand, prior, changes);
			}
			value += log_probability(term->binary.right->number.integer, prior.set_size_distribution);
			return value;
		} else if ((term->binary.right->type == hol_term_type::UNARY_APPLICATION
				 && term->binary.right->binary.left->type == hol_term_type::CONSTANT
				 && term->binary.right->binary.left->constant == (unsigned int) BuiltInPredicates::ARG1)
				|| (term->binary.left->type == hol_term_type::UNARY_APPLICATION
				 && term->binary.left->binary.left->type == hol_term_type::CONSTANT
				 && term->binary.left->binary.left->constant == (unsigned int) BuiltInPredicates::ARG1))
		{
			value = (Quantified ? prior.log_quantified_arg1_probability : prior.log_arg1_probability);
			hol_term* left = term->binary.left;
			hol_term* right = term->binary.right;
			if (left->type != hol_term_type::UNARY_APPLICATION)
				swap(left, right);
			if (left->binary.right->type == hol_term_type::CONSTANT)
				changes.constants.add_constant(left->binary.right->constant);
			if (right->type == hol_term_type::CONSTANT) {
				value += prior.log_arg_constant_probability;
				changes.constants.add_constant(right->constant);
				if (left->binary.right->type == hol_term_type::CONSTANT)
					changes.arg1_map.put(left->binary.right->constant, right->constant);
			} else if (right->type == hol_term_type::NUMBER) {
				value += prior.log_arg_number_probability;
				/* TODO: implement this */
				value += -7.0f;
			} else if (right->type == hol_term_type::STRING) {
				value += prior.log_arg_string_probability;
				/* TODO: implement this */
				value += -7.0f;
			}
			return value;
		} else if ((term->binary.right->type == hol_term_type::UNARY_APPLICATION
				 && term->binary.right->binary.left->type == hol_term_type::CONSTANT
				 && term->binary.right->binary.left->constant == (unsigned int) BuiltInPredicates::ARG2)
				|| (term->binary.left->type == hol_term_type::UNARY_APPLICATION
				 && term->binary.left->binary.left->type == hol_term_type::CONSTANT
				 && term->binary.left->binary.left->constant == (unsigned int) BuiltInPredicates::ARG2))
		{
			value = (Quantified ? prior.log_quantified_arg2_probability : prior.log_arg2_probability);
			hol_term* left = term->binary.left;
			hol_term* right = term->binary.right;
			if (left->type != hol_term_type::UNARY_APPLICATION)
				swap(left, right);
			changes.constants.add_constant(left->binary.right->constant);
			if (right->type == hol_term_type::CONSTANT) {
				value += prior.log_arg_constant_probability;
				changes.constants.add_constant(right->constant);
			} else if (right->type == hol_term_type::NUMBER) {
				value += prior.log_arg_number_probability;
				/* TODO: implement this */
				value += -7.0f;
			} else if (right->type == hol_term_type::STRING) {
				value += prior.log_arg_string_probability;
				changes.arg2_string_map.put(left->binary.right->constant, &right->str);
			}
			return value;
		} else {
			// TODO: implement this (we also have definition axioms for constants)
			return (Quantified ? prior.log_quantified_set_size_axiom_probability : prior.log_set_size_axiom_probability);
			/*fprintf(stderr, "log_probability_helper WARNING: Found equality term that is not a set size axiom: ");
			print(*term, stderr); print('\n', stderr);
			return -std::numeric_limits<double>::infinity();*/
		}
	case hol_term_type::EXISTS:
		value = (Quantified ? prior.log_quantified_existential_probability : prior.log_existential_probability);
		if (term->quantifier.operand->type == hol_term_type::AND) {
			hol_term* operand = term->quantifier.operand;
			value += (operand->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
			for (unsigned int i = 0; i < operand->array.length; i++)
				value += log_probability_helper<true, false>(operand->array.operands[i], prior, changes);
		} else {
			value += prior.log_antecedent_stop_probability + log_probability_helper<true, false>(term->quantifier.operand, prior, changes);
		}
		return value;
	case hol_term_type::AND:
		value = (Quantified ? prior.log_quantified_conjunction_probability : prior.log_conjunction_probability);
		value += (term->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
		for (unsigned int i = 0; i < term->array.length; i++)
			value += log_probability_helper<Quantified, false>(term->array.operands[i], prior, changes);
		return value;
	case hol_term_type::OR:
		value = (Quantified ? prior.log_quantified_disjunction_probability : prior.log_disjunction_probability);
		value += (term->array.length - 1) * prior.log_antecedent_continue_probability + prior.log_antecedent_stop_probability;
		for (unsigned int i = 0; i < term->array.length; i++)
			value += log_probability_helper<Quantified, false>(term->array.operands[i], prior, changes);
		return value;
	case hol_term_type::NOT:
		value = (Quantified ? prior.log_quantified_negated_expression_probability : prior.log_negated_expression_probability);
		value += log_probability_helper<Quantified, false>(term->unary.operand, prior, changes);
		return value;
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IF_THEN:
	case hol_term_type::IFF:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::LAMBDA:
	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::PARAMETER:
	case hol_term_type::NUMBER:
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

template<typename BuiltInPredicates,
	typename ConstantDistribution,
	typename SetSizeDistribution,
	template<typename> class Collection>
double log_probability(const Collection<hol_term*>& clusters,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior)
{
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of `simple_hol_term_distribution`:\n");
#endif
	double value = 0.0;
	typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>::prior_state_changes changes;
	array_map<unsigned int, unsigned int> arg1_map(8);
	array_map<unsigned int, string*> arg2_string_map(8);
	for (hol_term* formula : clusters) {
		double current_value = log_probability_helper(formula, prior, changes);
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "  Log probability of ");
		print(*formula, stderr);
		fprintf(stderr, " is %lf.\n", current_value);
#endif
		value += current_value;
	}
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total log likelihood: %lf.\n", value);
#endif

	/* get the name events */
	array_map<unsigned int, array<const string*>> name_map(max(1, changes.arg2_string_map.size));
	for (const auto& entry : changes.arg2_string_map) {
		bool contains;
		unsigned int named_entity = changes.arg1_map.get(entry.key, contains);
		if (!contains) continue;

		unsigned int index = name_map.index_of(named_entity);
		if (index == name_map.size) {
			name_map.keys[index] = named_entity;
			if (!array_init(name_map.values[index], 4)) {
				for (auto entry : name_map) free(entry.value);
				return false;
			}
			name_map.size++;
		} if (!name_map.values[index].add(entry.value)) {
			for (auto entry : name_map) free(entry.value);
			return false;
		}
	}

#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of name entities:\n");
#endif
	double name_prior = 0.0;
	for (const auto& entry : name_map) {
		double current_value = log_probability(entry.value.length, prior.name_count_distribution) - log_probability(0, prior.name_count_distribution);
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "  Concept %u has names ", entry.key);
		print(entry.value, stderr, pointer_scribe());
		fprintf(stderr, " with log probability %lf.\n", current_value);
#endif
		name_prior += current_value;
	}
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total: %lf\n", name_prior);
#endif
	for (auto entry : name_map) free(entry.value);
	value += name_prior;

	return value + log_probability(changes.constants, prior.constant_distribution);
}

template<typename BuiltInPredicates,
	typename ConstantDistribution,
	typename SetSizeDistribution,
	template<typename> class Collection>
double log_probability_ratio(
		const typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>::prior_state& prior_state,
		const Collection<hol_term*>& old_clusters,
		const Collection<hol_term*>& new_clusters,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>& prior,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>::prior_state_changes& old_prior_changes,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution>::prior_state_changes& new_prior_changes)
{
	double value = 0.0;
	for (hol_term* formula : new_clusters)
		value += log_probability_helper(formula, prior, new_prior_changes);
	for (hol_term* formula : old_clusters)
		value -= log_probability_helper(formula, prior, old_prior_changes);

	/* get the name events */
	array_map<unsigned int, array<const string*>> old_name_map(max(1, old_prior_changes.arg2_string_map.size));
	for (const auto& entry : old_prior_changes.arg2_string_map) {
		bool contains;
		unsigned int named_entity = prior_state.arg1_map.get(entry.key, contains);
		if (!contains) continue;

		unsigned int index = old_name_map.index_of(named_entity);
		if (index == old_name_map.size) {
			old_name_map.keys[index] = named_entity;
			if (!array_init(old_name_map.values[index], 4)) {
				for (auto entry : old_name_map) free(entry.value);
				return false;
			}
			old_name_map.size++;
		} if (!old_name_map.values[index].add(entry.value)) {
			for (auto entry : old_name_map) free(entry.value);
			return false;
		}
	} for (const auto& entry : old_prior_changes.arg1_map) {
		bool contains;
		if (old_prior_changes.arg2_string_map.contains(entry.key))
			continue;
		const string* name = prior_state.arg2_string_map.get(entry.key, contains);
		if (!contains) continue;

		unsigned int named_entity = entry.value;
		unsigned int index = old_name_map.index_of(named_entity);
		if (index == old_name_map.size) {
			old_name_map.keys[index] = named_entity;
			if (!array_init(old_name_map.values[index], 4)) {
				for (auto entry : old_name_map) free(entry.value);
				return false;
			}
			old_name_map.size++;
		} if (!old_name_map.values[index].add(name)) {
			for (auto entry : old_name_map) free(entry.value);
			return false;
		}
	}

	array_map<unsigned int, array<const string*>> new_name_map(max(1, new_prior_changes.arg2_string_map.size));
	for (const auto& entry : new_prior_changes.arg2_string_map) {
		bool contains;
		unsigned int named_entity = new_prior_changes.arg1_map.get(entry.key, contains);
		if (!contains) {
			if (old_prior_changes.arg1_map.contains(entry.key))
				continue;
			named_entity = prior_state.arg1_map.get(entry.key, contains);
			if (!contains) continue;
		}

		unsigned int index = new_name_map.index_of(named_entity);
		if (index == new_name_map.size) {
			new_name_map.keys[index] = named_entity;
			if (!array_init(new_name_map.values[index], 4)) {
				for (auto entry : old_name_map) free(entry.value);
				for (auto entry : new_name_map) free(entry.value);
				return false;
			}
			new_name_map.size++;
		} if (!new_name_map.values[index].add(entry.value)) {
			for (auto entry : old_name_map) free(entry.value);
			for (auto entry : new_name_map) free(entry.value);
			return false;
		}
	} for (const auto& entry : new_prior_changes.arg1_map) {
		bool contains;
		if (old_prior_changes.arg2_string_map.contains(entry.key))
			continue;
		const string* name = prior_state.arg2_string_map.get(entry.key, contains);
		if (!contains) continue;

		unsigned int named_entity = entry.value;
		unsigned int index = new_name_map.index_of(named_entity);
		if (index == new_name_map.size) {
			new_name_map.keys[index] = named_entity;
			if (!array_init(new_name_map.values[index], 4)) {
				for (auto entry : old_name_map) free(entry.value);
				for (auto entry : new_name_map) free(entry.value);
				return false;
			}
			new_name_map.size++;
		} if (!new_name_map.values[index].add(name)) {
			for (auto entry : old_name_map) free(entry.value);
			for (auto entry : new_name_map) free(entry.value);
			return false;
		}
	}

	if (old_name_map.size > 1) insertion_sort(old_name_map.keys, old_name_map.values, old_name_map.size, default_sorter());
	if (new_name_map.size > 1) insertion_sort(new_name_map.keys, new_name_map.values, new_name_map.size, default_sorter());

	/* get the full set of name events before the changes (TODO: we can actually keep `name_map` in `prior_state`) */
	array_map<unsigned int, array<const string*>> name_map(max(1, prior_state.arg2_string_map.size));
	for (const auto& entry : prior_state.arg2_string_map) {
		bool contains;
		unsigned int named_entity = prior_state.arg1_map.get(entry.key, contains);
		if (!contains) continue;

		unsigned int index = name_map.index_of(named_entity);
		if (index == name_map.size) {
			name_map.keys[index] = named_entity;
			if (!array_init(name_map.values[index], 4)) {
				for (auto entry : name_map) free(entry.value);
				for (auto entry : old_name_map) free(entry.value);
				for (auto entry : new_name_map) free(entry.value);
				return false;
			}
			name_map.size++;
		} if (!name_map.values[index].add(entry.value)) {
			for (auto entry : name_map) free(entry.value);
			for (auto entry : old_name_map) free(entry.value);
			for (auto entry : new_name_map) free(entry.value);
			return false;
		}
	}

	bool contains;
	double name_prior = 0.0;
	unsigned int i = 0, j = 0;
	while (i < old_name_map.size && j < new_name_map.size) {
		if (old_name_map.keys[i] == new_name_map.keys[j]) {
			array<const string*>& names = name_map.get(old_name_map.keys[i], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count + ((ssize_t) new_name_map.values[j].length - old_name_map.values[i].length), prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			i++; j++;
		} else if (old_name_map.keys[i] < new_name_map.keys[j]) {
			array<const string*>& names = name_map.get(old_name_map.keys[i], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count - ((ssize_t) old_name_map.values[i].length), prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			i++;
		} else {
			array<const string*>& names = name_map.get(new_name_map.keys[j], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count + new_name_map.values[j].length, prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			j++;
		}
	} while (i < old_name_map.size) {
		array<const string*>& names = name_map.get(old_name_map.keys[i], contains);
		unsigned int old_count = (contains ? names.length : 0);
		name_prior += log_probability(old_count - ((ssize_t) old_name_map.values[i].length), prior.name_count_distribution)
					- log_probability(old_count, prior.name_count_distribution);
		i++;
	} while (j < new_name_map.size) {
		array<const string*>& names = name_map.get(new_name_map.keys[j], contains);
		unsigned int old_count = (contains ? names.length : 0);
		name_prior += log_probability(old_count + new_name_map.values[j].length, prior.name_count_distribution)
					- log_probability(old_count, prior.name_count_distribution);
		j++;
	}
	for (auto entry : old_name_map) free(entry.value);
	for (auto entry : new_name_map) free(entry.value);
	for (auto entry : name_map) free(entry.value);
	value += name_prior;

	return value + log_probability_ratio(prior_state.constant_prior_state, old_prior_changes.constants, new_prior_changes.constants, prior.constant_distribution);
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
	static inline hol_term* canonicalize(const hol_term& src, array_map<unsigned int, unsigned int>& variable_map)
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

		hol_scope& scope = *((hol_scope*) alloca(sizeof(hol_scope)));
		if (!canonicalize_scope<AllConstantsDistinct>(src, scope, variable_map, types))
			return nullptr;
		hol_term* canonicalized = scope_to_term(scope);
		free(scope);
		return canonicalized;
	}

	template<bool Quiet = false>
	static inline hol_term* canonicalize(const hol_term& src)
	{
		array_map<unsigned int, unsigned int> variable_map(16);
		return canonicalize<Quiet>(src, variable_map);
	}
};

struct theory_initializer {
	array<instance> expected_constants;
	unsigned int constant_position;

	theory_initializer(unsigned int initial_capacity) : expected_constants(initial_capacity), constant_position(0) { }
};

template<typename Proof>
inline void visit_node(const Proof& proof, const theory_initializer& visitor) { }

template<bool Negated, typename Term> constexpr bool visit_unary_atom(const Term* term, const theory_initializer& visitor) { return true; }
template<bool Negated> constexpr bool visit_binary_atom(unsigned int predicate, unsigned int arg1, unsigned int arg2, const theory_initializer& visitor) { return true; }
template<typename Proof> constexpr bool visit_subset_axiom(const Proof& proof, const theory_initializer& visitor) { return true; }
constexpr bool visit_existential_intro(const theory_initializer& visitor) { return true; }
constexpr bool visit_negated_universal_intro(const theory_initializer& visitor) { return true; }
constexpr bool visit_negated_conjunction(const theory_initializer& visitor) { return true; }
constexpr bool visit_disjunction_intro(const theory_initializer& visitor) { return true; }

inline void on_subtract_changes(const theory_initializer& visitor) { }

template<typename Formula>
constexpr bool on_undo_filter_operands(const Formula* formula, const theory_initializer& visitor) { return true; }

template<typename Theory, typename Formula>
constexpr bool on_undo_filter_constants(const Theory& T, const Formula* quantified, unsigned int variable, const theory_initializer& visitor) { return true; }

template<typename Formula>
inline bool filter_operands(const Formula* formula, array<unsigned int>& indices, theory_initializer& initializer)
{
	return filter_operands(formula, indices);
}

template<typename ProofCalculus, typename Canonicalizer>
inline bool filter_constants(const theory<ProofCalculus, Canonicalizer>& T,
		const typename ProofCalculus::Language* formula,
		unsigned int variable, array<instance>& constants,
		theory_initializer& initializer)
{
	if (!filter_constants_helper(T, formula, variable, constants))
		return false;

	if (initializer.constant_position == initializer.expected_constants.length) {
		fprintf(stderr, "filter_constants ERROR: `theory_initializer` has no further `expected_constants`.\n");
		exit(EXIT_FAILURE);
	} else if (!constants.contains(initializer.expected_constants[initializer.constant_position])) {
		fprintf(stderr, "filter_constants ERROR: The next constant in `theory_initializer.expected_constants` is unavailable.\n");
		exit(EXIT_FAILURE);
	}
	constants[0] = initializer.expected_constants[initializer.constant_position++];
	constants.length = 1;
	return true;
}

template<typename Formula>
constexpr bool inconsistent_constant(const Formula* formula, unsigned int index, theory_initializer& initializer) { return true; }

template<typename Formula>
inline bool inconsistent_constant(const Formula* formula, const instance& constant, theory_initializer& initializer) {
	initializer.constant_position--;
	return true;
}

template<typename Formula>
inline void finished_constants(const Formula* formula, unsigned int original_constant_count, theory_initializer& initializer) { }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(
		unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out,
		unsigned int min_set_size,
		unsigned int max_set_size,
		theory_initializer& initializer)
{
	return compute_new_set_size(set_id, sets, out, min_set_size, max_set_size);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		theory_initializer& initializer)
{
	on_free_set(set_id, sets);
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
// to avoid breakpoints being moved due to eliminated code
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
			iid_uniform_distribution<unsigned int>(100), chinese_restaurant_process<unsigned int>(1.0, 0.0), make_dirichlet_process(1.0e-12, make_iid_uniform_distribution_from_log_n<hol_term>(24.0)));
	auto theory_element_prior = make_simple_hol_term_distribution<built_in_predicates>(constant_prior, geometric_distribution(0.2),
			0.0199999, 0.01, 0.0000001, 0.17, 0.1, 0.1, 0.2, 0.2, 0.2,
			0.1099999, 0.01, 0.0000001, 0.27, 0.1099999, 0.1, 0.0000001, 0.2, 0.2,
			0.999999998, 0.000000001, 0.000000001, 0.3, 0.4, 0.2, 0.4, 0.00000000001);
	auto axiom_prior = make_dirichlet_process(1.0e-1, theory_element_prior);
	auto conjunction_introduction_prior = uniform_subset_distribution<const nd_step<hol_term>*>(0.5);
	auto conjunction_elimination_prior = make_levy_process(poisson_distribution(2.0), poisson_distribution(1.0));
	auto universal_introduction_prior = unif_distribution<unsigned int>();
	auto universal_elimination_prior = chinese_restaurant_process<hol_term>(1.0, 0.0);
	auto term_indices_prior = make_levy_process(poisson_distribution(9.0), poisson_distribution(2.0));
	auto proof_prior = make_canonicalized_proof_prior(axiom_prior, conjunction_introduction_prior, conjunction_elimination_prior,
			universal_introduction_prior, universal_elimination_prior, term_indices_prior, poisson_distribution(20.0), 0.5);
	decltype(proof_prior)::PriorState proof_axioms;
	if (!parser.invert_name_map(names)) {
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

/* run RuleTaker experiments */
run_ruletaker_experiments(corpus, parser, T, proof_axioms, proof_prior, names, seed_entities, "ruletaker/birds-electricity/test.jsonl", 1);
for (auto entry : names) free(entry.key);
// to avoid breakpoints being moved due to eliminated code
if (seed_training_set.length > 0)
return EXIT_SUCCESS;

	/*read_sentence(corpus, parser, "A butterfly has a spot.", T, names, seed_entities, proof_prior, proof_axioms);

	parse_sentence_with_prior(parser, "Sally caught a butterfly with a spot.", T, names, proof_prior, proof_axioms);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;*/

	/*read_article(names.get("Des Moines"), corpus, parser, T, names, seed_entities, proof_prior);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;*/

/*theory_initializer initializer(29);
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 2));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());

initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_constant(1000000000 + 1));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 5));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());

initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_constant(1000000000 + 1));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 8));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());

initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_constant(1000000000 + 1));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 11));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_any());

initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 1));

initializer.expected_constants.add(instance_number(52069, 0));
initializer.expected_constants.add(instance_constant(1000000000 + 2));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 3));

initializer.expected_constants.add(instance_number(53179, 0));
initializer.expected_constants.add(instance_constant(1000000000 + 5));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 6));

initializer.expected_constants.add(instance_number(69899, 0));
initializer.expected_constants.add(instance_constant(1000000000 + 8));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 9));

initializer.expected_constants.add(instance_number(121590, 0));
initializer.expected_constants.add(instance_constant(1000000000 + 11));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 12));

initializer.expected_constants.add(instance_constant(1000000000 + 11));
initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_constant(1000000000 + 14));
initializer.expected_constants.add(instance_constant(1000000000 + 11));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 1));*/
/*initializer.expected_constants.add(instance_constant(1000000000 + 2));
initializer.expected_constants.add(instance_constant(1000000000 + 0));
initializer.expected_constants.add(instance_constant(1000000000 + 8));
initializer.expected_constants.add(instance_constant(1000000000 + 2));
initializer.expected_constants.add(instance_any());
initializer.expected_constants.add(instance_constant(1000000000 + 1));*/

	read_sentence(corpus, parser, "Louisiana is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "Arkansas is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "Oklahoma is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "New Mexico is a state that borders Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "There are 4 states that border Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of Louisiana is 52069 square miles.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of Arkansas is 53179 square miles.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of Oklahoma is 69899 square miles.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The area of New Mexico is 121590 square miles.", T, names, seed_entities, proof_prior, proof_axioms);

	/*read_sentence(corpus, parser, "The Red River is a river in Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The Canadian River is a river in Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The Rio Grande is a river in Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The Pecos River is a river in Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "The Washita River is a river in Texas.", T, names, seed_entities, proof_prior, proof_axioms);
	read_sentence(corpus, parser, "There are 5 rivers in Texas.", T, names, seed_entities, proof_prior, proof_axioms);*/

typedef decltype(T) TheoryType;
TheoryType& T_map = *((TheoryType*) alloca(sizeof(TheoryType)));
hash_map<const hol_term*, hol_term*> formula_map(128);
TheoryType::clone(T, T_map, formula_map);
auto collector = make_log_probability_collector(T, proof_prior);
double max_log_probability = collector.current_log_probability;
for (unsigned int t = 0; t < 200 /*20000*/; t++) {
	bool print_debug = false;
	if (print_debug) T.print_axioms(stderr, *debug_terminal_printer);
	if (print_debug) T.print_disjunction_introductions(stderr, *debug_terminal_printer);
	do_mh_step(T, proof_prior, proof_axioms, collector);
	if (collector.current_log_probability > max_log_probability) {
		free(T_map); formula_map.clear();
		TheoryType::clone(T, T_map, formula_map);
		max_log_probability = collector.current_log_probability;
	}
}
T_map.print_axioms(stderr, *debug_terminal_printer);
free(T_map);

	array<string> answers(4);
	/*if (answer_question<true>(answers, "Pittsburgh is in what state?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "Des Moines is located in what state?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "The population of Arizona is what?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	}*/ if (answer_question<true>(answers, "What is the largest state bordering Texas?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} /*if (answer_question<true>(answers, "What is the state with the highest population?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "Which is the longest river in the USA?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "San Antonio is in what state?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	} if (answer_question<true>(answers, "What are all the rivers in Texas?", 100000, corpus, parser, T, names, seed_entities, proof_prior, proof_axioms)) {
		print("Answers: ", stdout); print(answers, stdout); print('\n', stdout);
	}*/
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
