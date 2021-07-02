#ifndef THEORY_PRIOR_H_
#define THEORY_PRIOR_H_

#include <core/random.h>
#include <math/log.h>

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
	double min_cdf = 1.0 - pow(1 - prior.p, min);
	double max_cdf = (max == UINT_MAX ? 1.0 : (1.0 - pow(1 - prior.p, max + 1)));
	double u = sample_uniform<double>() * (max_cdf - min_cdf) + min_cdf;
	return floor(log(1.0 - u) / prior.log_one_minus_p);
}

struct very_light_tail_distribution {
	typedef unsigned int ObservationType;

	double log_lambda, log_normalization;

	very_light_tail_distribution(double log_lambda) : log_lambda(log_lambda), log_normalization(0.0) {
		for (unsigned int k = 0; ; k++) {
			double old_log_normalization = log_normalization;
			log_normalization = logsumexp(log_normalization, (k * k) * log_lambda);
			if (log_normalization == old_log_normalization)
				break;
		}
	}

	very_light_tail_distribution(const very_light_tail_distribution& src) : log_lambda(src.log_lambda), log_normalization(src.log_normalization) { }
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
inline void free_all(default_array_multiset<T>& a) { }

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
inline bool clone(const default_hash_multiset<T>& src, default_hash_multiset<T>& dst) {
	return init(dst, src);
}

template<typename T, typename Stream>
bool read(default_hash_multiset<T>& state, Stream& in) {
	decltype(state.a.counts.table.size) histogram_size;
	if (!read(histogram_size, in)
	 || !init(state.a, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (histogram_size == 0 ? 1 : histogram_size)) + 1)))
		return false;
	T& key = *((T*) alloca(sizeof(T)));
	for (unsigned int i = 0; i < histogram_size; i++) {
		if (!read(key, in)) {
			free(state.a);
			return false;
		}
		unsigned int bucket = state.a.counts.table.index_to_insert(key);
		if (!read(state.a.counts.values[bucket], in)) {
			free(key);
			free(state.a);
			return false;
		}
		move(key, state.a.counts.table.keys[bucket]);
		state.a.counts.table.size++;
		state.a.sum += state.a.counts.values[bucket];
	}
	return true;
}

template<typename T, typename Stream>
bool write(const default_hash_multiset<T>& state, Stream& out) {
	if (!write(state.a.counts.table.size, out))
		return false;
	for (const auto& entry : state.a.counts) {
		if (!write(entry.key, out)
		 || !write(entry.value, out))
			return false;
	}
	return true;
}

template<typename T>
inline void free_all(default_hash_multiset<T>& a) { }

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
constexpr bool clone(const empty_prior_state& src, const empty_prior_state& dst) { return true; }

template<typename Stream, typename... Reader>
constexpr bool read(const empty_prior_state& state, Stream& in, Reader&&... reader) { return true; }

template<typename Stream, typename... Writer>
constexpr bool write(const empty_prior_state& state, Stream& out, Writer&&... writer) { return true; }

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
	double normalization;
	if (prior.alpha <= 1.0e11)
		normalization = lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
	else normalization = -(double) sum(observations) * log(prior.alpha + sum(observations));
	if (prior.d == 0.0) {
		normalization += prior.log_alpha * size(observations);
	} else {
		normalization += prior.log_d * size(observations);
		double ratio = prior.alpha/prior.d;
		if (ratio <= 1.0e11)
			normalization += lgamma(ratio + size(observations)) - lgamma(ratio);
		else normalization += (double) size(observations) * log(ratio + size(observations));
	}
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
	double normalization;
	if (prior.alpha <= 1.0e11)
		normalization = lgamma(prior.alpha) - lgamma(prior.alpha + sum(observations));
	else normalization = -(double) sum(observations) * log(prior.alpha + sum(observations));
	if (prior.d == 0.0) {
		normalization += prior.log_alpha * size(observations);
	} else {
		normalization += prior.log_d * size(observations);
		double ratio = prior.alpha/prior.d;
		if (ratio <= 1.0e11)
			normalization += lgamma(ratio + size(observations)) - lgamma(ratio);
		else normalization += (double) size(observations) * log(ratio + size(observations));
	}
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Prior: %lf.\n", normalization);
#endif
	value += normalization;
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total: %lf.\n", value);
#endif
	return value;
}

template<bool AutomaticallyFree, typename T, typename Collection>
inline double log_probability_ratio_old_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection& old_clusters,
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

template<bool AutomaticallyFree, typename T, typename Collection>
inline double log_probability_ratio_new_cluster(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const T& observation, unsigned int frequency,
		const chinese_restaurant_process<T>& prior,
		Collection& new_clusters,
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

template<bool AutomaticallyFree, typename MultisetType, typename Collection, typename T>
double log_probability_ratio(
		const hash_multiset<T, AutomaticallyFree>& tables,
		const MultisetType& old_observations,
		const MultisetType& new_observations,
		const chinese_restaurant_process<T>& prior,
		Collection& old_clusters, Collection& new_clusters)
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

	if (prior.alpha + tables.sum <= 1.0e11) {
		value += lgamma(prior.alpha + tables.sum) - lgamma(prior.alpha + (tables.sum - sum(old_observations) + sum(new_observations)));
	} else {
		if (sum(new_observations) < sum(old_observations))
			value += -((double) sum(new_observations) - sum(old_observations)) * log(prior.alpha + tables.sum);
		else value += -((double) sum(new_observations) - sum(old_observations)) * log(prior.alpha + (tables.sum - sum(old_observations) + sum(new_observations)));
	}
	if (prior.d != 0.0) {
		if (prior.alpha/prior.d + tables.counts.table.size <= 1.0e11) {
			value += -lgamma(prior.alpha/prior.d + tables.counts.table.size) + lgamma(prior.alpha/prior.d + (tables.counts.table.size - old_cluster_count + new_cluster_count));
		} else {
			if (new_cluster_count < old_cluster_count)
				value += ((double) new_cluster_count - old_cluster_count) * log(prior.alpha/prior.d + tables.counts.table.size);
			else value += ((double) new_cluster_count - old_cluster_count) * log(prior.alpha/prior.d + (tables.counts.table.size - old_cluster_count + new_cluster_count));
		}
	}
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

		template<typename... Args>
		prior_state(const prior_state& src, Args&&... args) : base_prior_state(src.base_prior_state, std::forward<Args>(args)...) { }

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

		template<typename... Args>
		static inline bool get_formula_map(const prior_state& state, Args&&... args) {
			return BaseDistribution::PriorState::get_formula_map(state.base_prior_state, std::forward<Args>(args)...);
		}

		template<typename Stream, typename... Reader>
		static inline bool read(prior_state& state,
				Stream& in, Reader&&... reader)
		{
			return BaseDistribution::PriorState::read(state.base_prior_state, in, std::forward<Reader>(reader)...);
		}

		template<typename Stream, typename... Writer>
		static inline bool write(
				const prior_state& state,
				Stream& out, Writer&&... writer)
		{
			return BaseDistribution::PriorState::write(state.base_prior_state, out, std::forward<Writer>(writer)...);
		}

		template<typename... Cloner>
		static inline bool clone(const prior_state& src, prior_state& dst, Cloner&&... cloner)
		{
			return BaseDistribution::PriorState::clone(src.base_prior_state, dst.base_prior_state, std::forward<Cloner>(cloner)...);
		}

		static inline void free(prior_state& state) {
			core::free(state.base_prior_state);
		}
	};

	typedef typename BaseDistribution::ObservationType ObservationType;
	typedef default_array_multiset<ObservationType> ObservationCollection;
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
inline double log_probability(
		const MultisetType& restaurant_tables,
		const dirichlet_process<BaseDistribution>& prior,
		typename BaseDistribution::ObservationCollection& clusters)
{
	return log_probability(restaurant_tables, prior.restaurant, clusters);
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
bool debug_flag = false;
if (debug_flag) {
	print("new_clusters:\n", stderr);
	for (hol_term* formula : new_clusters) {
		print("  ", stderr); print(*formula, stderr); print('\n', stderr);
	}
	print("old_clusters:\n", stderr);
	for (hol_term* formula : old_clusters) {
		print("  ", stderr); print(*formula, stderr); print('\n', stderr);
	}
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
		hash_multiset<hol_term> root_types;

		prior_state() : types(16), root_types(16) { }

		prior_state(const prior_state& src) : constants(src.constants), predicates(src.predicates), types(src.types.table.capacity), root_types(src.root_types.counts.table.capacity) {
			if (!init_helper(src))
				exit(EXIT_FAILURE);
		}

		~prior_state() { free_helper(); }

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
			array_multiset<hol_term> new_root_types(8);
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

				if (!value.counts.check_size(value.counts.table.size + changes.types.values[i].counts.size)
				 || !new_root_types.counts.ensure_capacity(new_root_types.counts.size + changes.types.values[i].counts.size))
				{
					constants.subtract(changes.constants);
					predicates.subtract(changes.predicates);
					for (unsigned int j = 0; j < i; j++)
						types.get(changes.types.keys[j]).subtract(changes.types.values[j]);
					return false;
				}
				for (const auto& entry : changes.types.values[i].counts) {
					unsigned int& count = value.counts.get(entry.key, contains, index);
					if (contains) {
						count += entry.value;
					} else {
						value.counts.values[index] = entry.value;
						value.counts.table.keys[index] = entry.key;
						value.counts.table.size++;
						new_root_types.add(entry.key);
					}
					value.sum += entry.value;
				}
			}

			if (!root_types.add(new_root_types)) {
				constants.subtract(changes.constants);
				predicates.subtract(changes.predicates);
				for (unsigned int j = 0; j < changes.types.size; j++)
					types.get(changes.types.keys[j]).subtract(changes.types.values[j]);
				return false;
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

				const array_multiset<hol_term>& items = changes.types.values[i];
				for (unsigned int j = 0; j < items.counts.size; j++) {
					bool contains; unsigned int bucket;
					unsigned int& count = value.counts.get(items.counts.keys[j], contains, bucket);
#if !defined(NDEBUG)
					if (!contains || count < items.counts.values[j]) {
						fprintf(stderr, "simple_constant_distribution.prior_state.subtract "
								"WARNING: Attempted to remove more items from a bin than it contains.\n");
						count = 0;
					} else count -= items.counts.values[j];
#else
					count -= items.counts.values[j];
#endif
					if (count == 0) {
						value.counts.remove_at(bucket);
						bucket = root_types.counts.table.index_of(items.counts.keys[j]);
						if (root_types.counts.values[bucket] == 1) {
							root_types.counts.remove_at(bucket);
						} else {
							root_types.counts.values[bucket]--;
						}
						root_types.sum--;
					}
				}
				value.sum -= items.sum;
				if (value.sum == 0) {
					core::free(value);
					types.remove_at(index);
				}
			}
		}

		inline void subtract(const prior_state_changes& changes,
				const simple_constant_distribution<ConstantDistribution, PredicateDistribution, TypeDistribution>& prior)
		{
			subtract(changes);
		}

		static inline bool get_formula_map(const prior_state& state,
				hash_map<const hol_term*, unsigned int>& formula_map)
		{
			for (const auto& entry : state.types) {
				for (const auto& inner_entry : entry.value.counts) {
					if (!::get_formula_map(inner_entry.key, formula_map))
						return false;
				}
			} for (const auto& entry : state.root_types.counts) {
				if (!::get_formula_map(entry.key, formula_map))
					return false;
			}
			return true;
		}

		template<typename Stream>
		static inline bool read(prior_state& state,
				Stream& in, hol_term** terms)
		{
			if (!::read(state.constants, in)) {
				return false;
			} else if (!::read(state.predicates, in)) {
				core::free(state.constants);
				return false;
			}

			decltype(state.types.table.size) type_count;
			if (!core::read(type_count, in)
			 || !hash_map_init(state.types, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (type_count == 0 ? 1 : type_count)) + 1)))
			{
				core::free(state.constants);
				core::free(state.predicates);
				return false;
			}
			hol_term& term = *((hol_term*) alloca(sizeof(hol_term)));
			for (unsigned int i = 0; i < type_count; i++) {
				unsigned int key;
				unsigned int type_histogram_size;
				if (!core::read(key, in)
				 || !core::read(type_histogram_size, in))
				{
					state.free_helper();
					core::free(state.constants);
					core::free(state.predicates);
					core::free(state.types);
					return false;
				}
				unsigned int bucket = state.types.table.index_to_insert(key);
				hash_multiset<hol_term>& type_histogram = state.types.values[bucket];
				if (!init(type_histogram, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (type_histogram_size == 0 ? 1 : type_histogram_size)) + 1))) {
					state.free_helper();
					core::free(state.constants);
					core::free(state.predicates);
					core::free(state.types);
					return false;
				}
				state.types.table.keys[bucket] = key;
				state.types.table.size++;
				for (unsigned int j = 0; j < type_histogram_size; j++) {
					if (!::read(term, in, terms)) {
						state.free_helper();
						core::free(state.constants);
						core::free(state.predicates);
						core::free(state.types);
						return false;
					}
					term.reference_count = 1;
					bucket = type_histogram.counts.table.index_to_insert(term);
					if (!core::read(type_histogram.counts.values[bucket], in)) {
						core::free(term);
						state.free_helper();
						core::free(state.constants);
						core::free(state.predicates);
						core::free(state.types);
						return false;
					}
					move(term, type_histogram.counts.table.keys[bucket]);
					type_histogram.counts.table.size++;
					type_histogram.sum += type_histogram.counts.values[bucket];
				}
			}

			decltype(state.root_types.counts.table.size) root_type_count;
			if (!core::read(root_type_count, in)
			 || !init(state.root_types, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (root_type_count == 0 ? 1 : root_type_count)) + 1)))
			{
				state.free_helper();
				core::free(state.constants);
				core::free(state.predicates);
				core::free(state.types);
				return false;
			}
			for (unsigned int i = 0; i < root_type_count; i++) {
				if (!::read(term, in, terms)) {
					core::free(state);
					return false;
				}
				term.reference_count = 1;
				unsigned int bucket = state.root_types.counts.table.index_to_insert(term);
				if (!core::read(state.root_types.counts.values[bucket], in)) {
					core::free(term);
					core::free(state);
					return false;
				}
				move(term, state.root_types.counts.table.keys[bucket]);
				state.root_types.counts.table.size++;
				state.root_types.sum += state.root_types.counts.values[bucket];
			}

			return true;
		}

		template<typename Stream>
		static inline bool write(const prior_state& state, Stream& out,
				const hash_map<const hol_term*, unsigned int>& formula_map)
		{
			if (!::write(state.constants, out)
			 || !::write(state.predicates, out)
			 || !core::write(state.types.table.size, out))
				return false;
			for (const auto& entry : state.types) {
				if (!core::write(entry.key, out)
				 || !core::write((unsigned int) entry.value.counts.table.size, out))
					return false;
				for (const auto& inner_entry : entry.value.counts) {
					if (!::write(inner_entry.key, out, formula_map)
					 || !core::write(inner_entry.value, out))
						return false;
				}
			}

			if (!core::write(state.root_types.counts.table.size, out))
				return false;
			for (const auto& entry : state.root_types.counts) {
				if (!::write(entry.key, out, formula_map)
				 || !core::write(entry.value, out))
					return false;
			}
			return true;
		}

		static bool clone(const prior_state& src, prior_state& dst) {
			if (!::clone(src.constants, dst.constants)) {
				return false;
			} else if (!::clone(src.predicates, dst.predicates)) {
				core::free(dst.constants);
				return false;
			} else if (!hash_map_init(dst.types, src.types.table.capacity)) {
				core::free(dst.constants);
				core::free(dst.predicates);
				return false;
			} else if (!init(dst.root_types, src.root_types.counts.table.capacity)) {
				core::free(dst.constants);
				core::free(dst.predicates);
				core::free(dst.types);
				return false;
			} else if (!dst.init_helper(src)) {
				core::free(dst.constants);
				core::free(dst.predicates);
				core::free(dst.types);
				core::free(dst.root_types);
				return false;
			}
			return true;
		}

		static inline void free(prior_state& state) {
			state.free_helper();
			core::free(state.constants);
			core::free(state.predicates);
			core::free(state.types);
			core::free(state.root_types);
		}

		inline bool init_helper(const prior_state& src) {
			for (const auto& entry : src.types) {
				bool contains; unsigned int bucket;
				types.get(entry.key, contains, bucket);
				if (!init(types.values[bucket], entry.value.counts.table.capacity))
					return false;
				types.values[bucket].counts.put_all(entry.value.counts);
				types.values[bucket].sum = entry.value.sum;
				types.table.keys[bucket] = entry.key;
				types.table.size++;
			}
			root_types.counts.put_all(src.root_types.counts);
			root_types.sum = src.root_types.sum;
			return true;
		}

		inline void free_helper() {
			for (auto entry : types)
				core::free(entry.value);
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

	typedef typename decltype(prior.type_distribution.base_distribution)::ObservationCollection Clusters;
	Clusters clusters;
	for (const auto& entry : constants.types) {
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "log_probability of type distribution for constant %u:\n", entry.key);
#endif
		value += log_probability(entry.value, prior.type_distribution, clusters);
	}
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "log_probability of root type distribution:\n");
#endif
	value += log_probability(clusters, prior.type_distribution.base_distribution);
	free_all(clusters);
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
	array_multiset<hol_term> old_clusters(8);
	array_multiset<hol_term> new_clusters(8);
	while (i < old_constants.types.size && j < new_constants.types.size) {
		if (old_constants.types.keys[i] == new_constants.types.keys[j]) {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
			if (contains) {
				value += log_probability_ratio(existing, old_constants.types.values[i], new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
			} else {
				value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
			}
			i++; j++;
		} else if (old_constants.types.keys[i] < new_constants.types.keys[j]) {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
			if (contains) {
				value += log_probability_ratio(existing, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution.restaurant, old_clusters, new_clusters);
			} else {
				value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution.restaurant, old_clusters, new_clusters);
			}
			i++;
		} else {
			const hash_multiset<hol_term>& existing = existing_constants.types.get(new_constants.types.keys[j], contains);
			if (contains) {
				value += log_probability_ratio(existing, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
			} else {
				value += log_probability_ratio(empty_type_prior_state, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
			}
			j++;
		}
	} while (i < old_constants.types.size) {
		const hash_multiset<hol_term>& existing = existing_constants.types.get(old_constants.types.keys[i], contains);
		if (contains) {
			value += log_probability_ratio(existing, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution.restaurant, old_clusters, new_clusters);
		} else {
			value += log_probability_ratio(empty_type_prior_state, old_constants.types.values[i], empty_type_prior_state_changes, prior.type_distribution.restaurant, old_clusters, new_clusters);
		}
		i++;
	} while (j < new_constants.types.size) {
		const hash_multiset<hol_term>& existing = existing_constants.types.get(new_constants.types.keys[j], contains);
		if (contains) {
			value += log_probability_ratio(existing, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
		} else {
			value += log_probability_ratio(empty_type_prior_state, empty_type_prior_state_changes, new_constants.types.values[j], prior.type_distribution.restaurant, old_clusters, new_clusters);
		}
		j++;
	}

	default_array<hol_term> old_prior_changes, new_prior_changes;
	value += log_probability_ratio(existing_constants.root_types, old_clusters, new_clusters, prior.type_distribution.base_distribution, type_prior_state.base_prior_state, old_prior_changes, new_prior_changes);
	free_all(old_prior_changes);
	free_all(new_prior_changes);
	return value;
}

template<typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
struct simple_hol_term_distribution
{
	struct prior_state_changes {
		typename ConstantDistribution::PriorStateChanges constants;
		array_map<unsigned int, unsigned int> arg1_map;
		array_map<unsigned int, hol_term*> arg2_string_map;
		array_map<unsigned int, array<hol_term>> definitions;

		prior_state_changes() : arg1_map(4), arg2_string_map(4), definitions(4) { }

		~prior_state_changes() {
			for (auto entry : definitions) {
				for (hol_term& term : entry.value)
					core::free(term);
				core::free(entry.value);
			}
		}

		bool add_definition(unsigned int constant, const hol_term* definition) {
			if (definition->type == hol_term_type::LAMBDA
			 && definition->quantifier.operand->type == hol_term_type::UNARY_APPLICATION
			 && definition->quantifier.operand->binary.left->type == hol_term_type::CONSTANT
			 && definition->quantifier.operand->binary.left->constant == constant
			 && definition->quantifier.operand->binary.right->type == hol_term_type::VARIABLE
			 && definition->quantifier.operand->binary.right->variable == definition->quantifier.variable)
				return true;
			if (!definitions.ensure_capacity(definitions.size + 1))
				return false;
			unsigned int index = definitions.index_of(constant);
			if (index == definitions.size) {
				if (!array_init(definitions.values[index], 4))
					return false;
				definitions.keys[index] = constant;
				definitions.size++;
				if (!init(definitions.values[index][0], *definition))
					return false;
				definitions.values[index].length = 1;
				return true;
			} else {
#if !defined(NDEBUG)
				if (definitions.values[index].contains(*definition))
					fprintf(stderr, "simple_hol_term_distribution.prior_state_changes.add_definition ERROR: `definitions[%u]` already contains the given formula.\n", index);
#endif
				return definitions.values[index].add(*definition);
			}
		}
	};

	struct prior_state {
		typename ConstantDistribution::PriorState constant_prior_state;
		array_map<unsigned int, unsigned int> arg1_map;
		array_map<unsigned int, hol_term*> arg2_string_map;
		array_map<unsigned int, array<hol_term>> definitions;

		prior_state() : arg1_map(16), arg2_string_map(16), definitions(16) { }

		prior_state(const prior_state& src, const hash_map<const hol_term*, hol_term*>& formula_map) :
			constant_prior_state(src.constant_prior_state), arg1_map(src.arg1_map.capacity), arg2_string_map(src.arg2_string_map.capacity), definitions(src.definitions.capacity)
		{
			if (!init_helper(src, formula_map))
				throw std::bad_alloc();
		}

		~prior_state() { free_helper(); }

		inline bool add(const prior_state_changes& changes) {
			for (const auto& entry : changes.arg1_map) {
				if (!arg1_map.ensure_capacity(arg1_map.size + 2))
					return false;
				unsigned int index = arg1_map.index_of(entry.key);
#if !defined(NDEBUG)
				if (index < arg1_map.size)
					fprintf(stderr, "simple_hol_term_distribution.prior_state.add ERROR: `arg1_map` already contains the key %u.\n", entry.key);
#endif
				arg1_map.keys[index] = entry.key;
				arg1_map.values[index] = entry.value;
				arg1_map.size++;
			} for (const auto& entry : changes.arg2_string_map) {
				if (!arg2_string_map.ensure_capacity(arg2_string_map.size + 2))
					return false;
				unsigned int index = arg2_string_map.index_of(entry.key);
#if !defined(NDEBUG)
				if (index < arg2_string_map.size)
					fprintf(stderr, "simple_hol_term_distribution.prior_state.add ERROR: `arg2_string_map` already contains the key %u.\n", entry.key);
#endif
				arg2_string_map.keys[index] = entry.key;
				arg2_string_map.values[index] = entry.value;
				arg2_string_map.size++;
			} for (const auto& entry : changes.definitions) {
				if (!definitions.ensure_capacity(definitions.size + 1))
					return false;
				unsigned int index = definitions.index_of(entry.key);
				if (index == definitions.size) {
					if (!array_init(definitions.values[index], entry.value.capacity))
						return false;
					for (unsigned int i = 0; i < entry.value.length; i++)
						definitions.values[index][i] = entry.value[i];
					definitions.values[index].length = entry.value.length;
					definitions.keys[index] = entry.key;
					definitions.size++;
				} else {
#if !defined(NDEBUG)
					for (const hol_term& term : entry.value) {
						if (definitions.values[index].contains(term))
							fprintf(stderr, "simple_hol_term_distribution.prior_state.add ERROR: `definitions.values[%u]` already contains the given formula.\n", index);
					}
#endif
					if (!definitions.values[index].ensure_capacity(definitions.values[index].length + entry.value.length))
						return false;
					for (unsigned int i = 0; i < entry.value.length; i++) {
						definitions.values[index][definitions.values[index].length] = entry.value[i];
						definitions.values[index].length++;
					}
				}
			}
			return constant_prior_state.add(changes.constants);
		}

		inline bool add(const hol_term* observation,
				const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior)
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
			} for (const auto& entry : changes.definitions) {
				unsigned int index = definitions.index_of(entry.key);
#if !defined(NDEBUG)
				if (index == definitions.size)
					fprintf(stderr, "simple_hol_term_distribution.prior_state.subtract ERROR: `definitions` doesn't contain the given formula.\n");
#endif
				for (const hol_term& term : entry.value) {
					unsigned int term_index = definitions.values[index].index_of(term);
#if !defined(NDEBUG)
					if (term_index == definitions.values[index].length)
						fprintf(stderr, "simple_hol_term_distribution.prior_state.subtract ERROR: `definitions.values[%u]` doesn't contain the given formula.\n", index);
#endif
					core::free(definitions.values[index][term_index]);
					definitions.values[index].remove(term_index);
				}
				if (definitions.values[index].length == 0) {
					core::free(definitions.values[index]);
					definitions.remove_at(index);
				}
			}
		}

		inline void subtract(const hol_term* observation,
				const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior)
		{
			prior_state_changes changes;
			log_probability_helper(observation, prior, changes);
			return subtract(changes);
		}

		static inline bool get_formula_map(const prior_state& state,
				hash_map<const hol_term*, unsigned int>& formula_map)
		{
			for (const auto& entry : state.arg2_string_map) {
				if (!::get_formula_map(entry.value, formula_map))
					return false;
			} for (const auto& entry : state.definitions) {
				for (const hol_term& definition : entry.value) {
					if (!::get_formula_map(definition, formula_map))
						return false;
				}
			}
			return ConstantDistribution::PriorState::get_formula_map(state.constant_prior_state, formula_map);
		}

		template<typename Stream>
		static inline bool read(prior_state& state,
				Stream& in, hol_term** terms)
		{
			decltype(state.arg1_map.size) arg1_map_size;
			if (!core::read(arg1_map_size, in)
			 || !array_map_init(state.arg1_map, ((size_t) 1) << (core::log2(arg1_map_size == 0 ? 1 : arg1_map_size) + 1)))
				return false;
			for (unsigned int i = 0; i < arg1_map_size; i++) {
				if (!core::read(state.arg1_map.keys[i], in)
				 || !core::read(state.arg1_map.values[i], in))
				{
					core::free(state.arg1_map);
					return false;
				}
				state.arg1_map.size++;
			}

			decltype(state.arg2_string_map.size) arg2_string_map_size;
			if (!core::read(arg2_string_map_size, in)
			 || !array_map_init(state.arg2_string_map, ((size_t) 1) << (core::log2(arg2_string_map_size == 0 ? 1 : arg2_string_map_size) + 1)))
			{
				core::free(state.arg1_map);
				return false;
			}
			for (unsigned int i = 0; i < arg2_string_map_size; i++) {
				unsigned int index;
				if (!core::read(state.arg2_string_map.keys[i], in)
				 || !core::read(index, in))
				{
					core::free(state.arg1_map);
					core::free(state.arg2_string_map);
					return false;
				}
				state.arg2_string_map.values[i] = terms[index];
				state.arg2_string_map.size++;
			}

			decltype(state.definitions.size) definition_count;
			if (!core::read(definition_count, in)
			 || !array_map_init(state.definitions, ((size_t) 1) << (core::log2(definition_count == 0 ? 1 : definition_count) + 1)))
			{
				core::free(state.arg1_map);
				core::free(state.arg2_string_map);
				return false;
			}
			for (unsigned int i = 0; i < definition_count; i++) {
				size_t definition_array_length;
				if (!core::read(state.definitions.keys[i], in)
				 || !core::read(definition_array_length, in)
				 || !array_init(state.definitions.values[i], ((size_t) 1) << (core::log2(definition_array_length == 0 ? 1 : definition_array_length) + 1)))
				{
					state.free_helper();
					core::free(state.arg1_map);
					core::free(state.arg2_string_map);
					core::free(state.definitions);
					return false;
				}
				state.definitions.size++;
				for (size_t j = 0; j < definition_array_length; j++) {
					if (!::read(state.definitions.values[i][j], in, terms)) {
						state.free_helper();
						core::free(state.arg1_map);
						core::free(state.arg2_string_map);
						core::free(state.definitions);
						return false;
					}
					state.definitions.values[i][j].reference_count = 1;
					state.definitions.values[i].length++;
				}
			}

			if (!ConstantDistribution::PriorState::read(state.constant_prior_state, in, terms)) {
				state.free_helper();
				core::free(state.arg1_map);
				core::free(state.arg2_string_map);
				core::free(state.definitions);
				return false;
			}
			return true;
		}

		template<typename Stream>
		static inline bool write(const prior_state& state, Stream& out,
				const hash_map<const hol_term*, unsigned int>& formula_map)
		{
			if (!core::write(state.arg1_map.size, out))
				return false;
			for (const auto& entry : state.arg1_map) {
				if (!core::write(entry.key, out)
				 || !core::write(entry.value, out))
					return false;
			}

			if (!core::write(state.arg2_string_map.size, out))
				return false;
			for (const auto& entry : state.arg2_string_map) {
				if (!core::write(entry.key, out)
				 || !core::write(formula_map.get(entry.value), out))
					return false;
			}

			if (!core::write(state.definitions.size, out))
				return false;
			for (const auto& entry : state.definitions) {
				if (!core::write(entry.key, out)
				 || !core::write((size_t) entry.value.length, out))
					return false;
				for (const hol_term& term : entry.value)
					if (!::write(term, out, formula_map)) return false;
			}

			return ConstantDistribution::PriorState::write(state.constant_prior_state, out, formula_map);
		}

		static bool clone(const prior_state& src, prior_state& dst,
				const hash_map<const hol_term*, hol_term*>& formula_map)
		{
			if (!ConstantDistribution::PriorState::clone(src.constant_prior_state, dst.constant_prior_state)) {
				return false;
			} else if (!array_map_init(dst.arg1_map, src.arg1_map.capacity)) {
				core::free(dst.constant_prior_state);
				return false;
			} else if (!array_map_init(dst.arg2_string_map, src.arg2_string_map.capacity)) {
				core::free(dst.constant_prior_state);
				core::free(dst.arg1_map);
				return false;
			} else if (!array_map_init(dst.definitions, src.definitions.capacity)) {
				core::free(dst.constant_prior_state);
				core::free(dst.arg1_map);
				core::free(dst.arg2_string_map);
				return false;
			} else if (!dst.init_helper(src, formula_map)) {
				core::free(dst.constant_prior_state);
				core::free(dst.arg1_map);
				core::free(dst.arg2_string_map);
				core::free(dst.definitions);
				return false;
			}
			return true;
		}

		static inline void free(prior_state& state) {
			state.free_helper();
			core::free(state.constant_prior_state);
			core::free(state.arg1_map);
			core::free(state.arg2_string_map);
			core::free(state.definitions);
		}

		inline bool init_helper(const prior_state& src,
				const hash_map<const hol_term*, hol_term*>& formula_map)
		{
			for (unsigned int i = 0; i < src.arg1_map.size; i++) {
				arg1_map.keys[i] = src.arg1_map.keys[i];
				arg1_map.values[i] = src.arg1_map.values[i];
				arg1_map.size++;
			} for (unsigned int i = 0; i < src.arg2_string_map.size; i++) {
				arg2_string_map.keys[i] = src.arg2_string_map.keys[i];
				arg2_string_map.values[i] = formula_map.get(src.arg2_string_map.values[i]);
				arg2_string_map.size++;
			} for (unsigned int i = 0; i < src.definitions.size; i++) {
				if (!array_init(definitions.values[i], src.definitions.values[i].capacity)) {
					for (unsigned int j = 0; j < i; j++) core::free(definitions.values[j]);
					return false;
				}
				definitions.keys[i] = src.definitions.keys[i];
				for (unsigned int j = 0; j < src.definitions.values[i].length; j++)
					definitions.values[i][j] = src.definitions.values[i][j];
				definitions.values[i].length = src.definitions.values[i].length;
				definitions.size++;
			}
			return true;
		}

		inline void free_helper() {
			for (auto entry : definitions) {
				for (hol_term& term : entry.value)
					core::free(term);
				core::free(entry.value);
			}
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
	double log_implication_probability;
	double log_set_size_axiom_probability;
	double log_arg1_probability;
	double log_arg2_probability;

	double log_quantified_ground_literal_probability;
	double log_quantified_negated_expression_probability;
	double log_quantified_universal_probability;
	double log_quantified_existential_probability;
	double log_quantified_conjunction_probability;
	double log_quantified_disjunction_probability;
	double log_quantified_implication_probability;
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
	DefinitionCountDistribution definition_count_distribution;

	simple_hol_term_distribution(
			const ConstantDistribution& constant_distribution,
			const SetSizeDistribution& set_size_distribution,
			const DefinitionCountDistribution& definition_count_distribution,
			double ground_literal_probability, double negated_expression_probability,
			double universal_probability, double existential_probability,
			double conjunction_probability, double disjunction_probability,
			double implication_probability, double set_size_axiom_probability,
			double arg1_probability, double arg2_probability,
			double quantified_ground_literal_probability, double quantified_negated_expression_probability,
			double quantified_universal_probability, double quantified_existential_probability,
			double quantified_conjunction_probability, double quantified_disjunction_probability,
			double quantified_implication_probability, double quantified_set_size_axiom_probability,
			double quantified_arg1_probability, double quantified_arg2_probability,
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
		log_implication_probability(log(implication_probability)),
		log_set_size_axiom_probability(log(set_size_axiom_probability)),
		log_arg1_probability(log(arg1_probability)),
		log_arg2_probability(log(arg2_probability)),
		log_quantified_ground_literal_probability(log(quantified_ground_literal_probability)),
		log_quantified_negated_expression_probability(log(quantified_negated_expression_probability)),
		log_quantified_universal_probability(log(quantified_universal_probability)),
		log_quantified_existential_probability(log(quantified_existential_probability)),
		log_quantified_conjunction_probability(log(quantified_conjunction_probability)),
		log_quantified_disjunction_probability(log(quantified_disjunction_probability)),
		log_quantified_implication_probability(log(quantified_implication_probability)),
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
		set_size_distribution(set_size_distribution),
		definition_count_distribution(definition_count_distribution)
	{
		if (fabs(ground_literal_probability + negated_expression_probability + universal_probability + existential_probability + conjunction_probability + disjunction_probability + implication_probability + set_size_axiom_probability + arg1_probability + arg2_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `ground_literal_probability + negated_expression_probability + universal_probability + existential_probability + conjunction_probability + disjunction_probability + implication_probability + set_size_axiom_probability + arg1_probability + arg2_probability` is not 1.\n");
		if (fabs(quantified_ground_literal_probability + quantified_negated_expression_probability + quantified_universal_probability + quantified_existential_probability + quantified_conjunction_probability + quantified_disjunction_probability + quantified_implication_probability + quantified_set_size_axiom_probability + quantified_arg1_probability + quantified_arg2_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `quantified_ground_literal_probability + quantified_negated_expression_probability + quantified_universal_probability + quantified_existential_probability + quantified_conjunction_probability + quantified_disjunction_probability + quantified_implication_probability + quantified_set_size_axiom_probability + quantified_arg1_probability + quantified_arg2_probability` is not 1.\n");
		if (fabs(arg_constant_probability + arg_number_probability + arg_string_probability - 1.0) > 1.0e-12)
			fprintf(stderr, "simple_hol_term_distribution WARNING: `arg_constant_probability + arg_number_probability + arg_string_probability` is not 1.\n");
	}

	simple_hol_term_distribution(const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& src) :
		log_ground_literal_probability(src.log_ground_literal_probability),
		log_negated_expression_probability(src.log_negated_expression_probability),
		log_universal_probability(src.log_universal_probability),
		log_existential_probability(src.log_existential_probability),
		log_conjunction_probability(src.log_conjunction_probability),
		log_disjunction_probability(src.log_disjunction_probability),
		log_implication_probability(src.log_implication_probability),
		log_set_size_axiom_probability(src.log_set_size_axiom_probability),
		log_arg1_probability(src.log_arg1_probability),
		log_arg2_probability(src.log_arg2_probability),
		log_quantified_ground_literal_probability(src.log_quantified_ground_literal_probability),
		log_quantified_negated_expression_probability(src.log_quantified_negated_expression_probability),
		log_quantified_universal_probability(src.log_quantified_universal_probability),
		log_quantified_existential_probability(src.log_quantified_existential_probability),
		log_quantified_conjunction_probability(src.log_quantified_conjunction_probability),
		log_quantified_disjunction_probability(src.log_quantified_disjunction_probability),
		log_quantified_implication_probability(src.log_quantified_implication_probability),
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
		set_size_distribution(src.set_size_distribution),
		definition_count_distribution(src.definition_count_distribution)
	{ }
};

template<typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
inline simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution> make_simple_hol_term_distribution(
		const ConstantDistribution& constant_distribution,
		const SetSizeDistribution set_size_distribution,
		const DefinitionCountDistribution& definition_count_distribution,
		double ground_literal_probability, double negated_expression_probability,
		double universal_probability, double existential_probability,
		double conjunction_probability, double disjunction_probability,
		double implication_probability, double set_size_axiom_probability,
		double arg1_probability, double arg2_probability,
		double quantified_ground_literal_probability, double quantified_negated_expression_probability,
		double quantified_universal_probability, double quantified_existential_probability,
		double quantified_conjunction_probability, double quantified_disjunction_probability,
		double quantified_implication_probability, double quantified_set_size_axiom_probability,
		double quantified_arg1_probability, double quantified_arg2_probability,
		double arg_constant_probability, double arg_number_probability,
		double arg_string_probability, double negation_probability,
		double unary_probability, double antecedent_stop_probability,
		double consequent_stop_probability, double name_count_parameter)
{
	return simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>(
			constant_distribution, set_size_distribution, definition_count_distribution,
			ground_literal_probability, negated_expression_probability,
			universal_probability, existential_probability,
			conjunction_probability, disjunction_probability,
			implication_probability, set_size_axiom_probability,
			arg1_probability, arg2_probability,
			quantified_ground_literal_probability, quantified_negated_expression_probability,
			quantified_universal_probability, quantified_existential_probability,
			quantified_conjunction_probability, quantified_disjunction_probability,
			quantified_implication_probability, quantified_set_size_axiom_probability,
			quantified_arg1_probability, quantified_arg2_probability,
			arg_constant_probability, arg_number_probability,
			arg_string_probability, negation_probability,
			unary_probability, antecedent_stop_probability,
			consequent_stop_probability, name_count_parameter);
}

template<bool Quantified, bool IsRoot, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
double log_probability_atom(const hol_term* function, const hol_term* arg1,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior,
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
	} else if (!Quantified && arg1->type == hol_term_type::NUMBER) {
		/* TODO: implement this */
		return prior.log_arg_number_probability - 7.0f;
	} else {
		return -std::numeric_limits<double>::infinity();
	}
}

template<bool Quantified, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
double log_probability_atom(
		const hol_term* function, const hol_term* arg1, const hol_term* arg2,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior,
		typename ConstantDistribution::ObservationCollection& constants)
{
	if (function->type != hol_term_type::CONSTANT)
		return -std::numeric_limits<double>::infinity();

	constants.add_predicate(function->constant);
	if (arg1->type == hol_term_type::VARIABLE) {
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(5);
		} else if (arg2->type == hol_term_type::CONSTANT) {
			constants.add_constant(arg2->constant);
			return prior.log_binary_probability - log_cache<double>::instance().get(5);
		} else if (arg2->type == hol_term_type::NUMBER) {
			/* TODO: add a proper prior distribution for this case */
			return prior.log_binary_probability - log_cache<double>::instance().get(5);
		} else {
			return -std::numeric_limits<double>::infinity();
		}
	} else if (arg1->type == hol_term_type::CONSTANT) {
		constants.add_constant(arg1->constant);
		if (arg2->type == hol_term_type::VARIABLE) {
			return prior.log_binary_probability - log_cache<double>::instance().get(5);
		} else if (arg2->type == hol_term_type::NUMBER) {
			/* TODO: add a proper prior distribution for this case */
			return prior.log_binary_probability - log_cache<double>::instance().get(5);
		} else { /* both `arg1` and `arg2` can't be constants */
			return -std::numeric_limits<double>::infinity();
		}
	} else {
		return -std::numeric_limits<double>::infinity();
	}
}

template<bool Quantified, bool IsRoot, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
double log_probability_literal(const hol_term* literal,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior,
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

template<bool Quantified = false, bool IsRoot = true, typename BuiltInPredicates, typename ConstantDistribution, typename SetSizeDistribution, typename DefinitionCountDistribution>
double log_probability_helper(const hol_term* term,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>::prior_state_changes& changes)
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
			if (left->binary.right->type == hol_term_type::CONSTANT)
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
				if (left->binary.right->type == hol_term_type::CONSTANT)
					changes.arg2_string_map.put(left->binary.right->constant, right);
			}
			return value;
		} else if (term->binary.left->type == hol_term_type::CONSTANT || term->binary.right->type == hol_term_type::CONSTANT) {
			if (IsRoot) {
				hol_term* left = term->binary.left;
				hol_term* right = term->binary.right;
				if (left->type != hol_term_type::CONSTANT)
					swap(left, right);
				changes.add_definition(left->constant, right);
			}
			return (Quantified ? prior.log_quantified_set_size_axiom_probability : prior.log_set_size_axiom_probability);
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
	case hol_term_type::IF_THEN:
		value = (Quantified ? prior.log_quantified_implication_probability : prior.log_implication_probability);
		value += log_probability_helper<Quantified, false>(term->binary.left, prior, changes);
		value += log_probability_helper<Quantified, false>(term->binary.right, prior, changes);
		return value;
	case hol_term_type::UNARY_APPLICATION:
	case hol_term_type::BINARY_APPLICATION:
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
	typename DefinitionCountDistribution,
	template<typename> class Collection>
double log_probability(const Collection<hol_term*>& clusters,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior)
{
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of `simple_hol_term_distribution`:\n");
#endif
	double value = 0.0;
	typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>::prior_state_changes changes;
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
	array_map<unsigned int, array<const hol_term*>> name_map(max(1, changes.arg2_string_map.size));
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

#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "log_probability of definitions:\n");
#endif
	double definition_prior = 0.0;
	for (const auto& entry : changes.definitions) {
		double current_value = log_probability(entry.value.length, prior.definition_count_distribution) - log_probability(0, prior.definition_count_distribution);
#if defined(DEBUG_LOG_PROBABILITY)
		fprintf(stderr, "  Concept %u has definitions ", entry.key);
		print(entry.value, stderr);
		fprintf(stderr, " with log probability %lf.\n", current_value);
#endif
		definition_prior += current_value;
	}
#if defined(DEBUG_LOG_PROBABILITY)
	fprintf(stderr, "  Total: %lf\n", definition_prior);
#endif
	value += definition_prior;

	return value + log_probability(changes.constants, prior.constant_distribution);
}

template<typename BuiltInPredicates,
	typename ConstantDistribution,
	typename SetSizeDistribution,
	typename DefinitionCountDistribution,
	template<typename> class Collection>
double log_probability_ratio(
		const typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>::prior_state& prior_state,
		const Collection<hol_term*>& old_clusters,
		const Collection<hol_term*>& new_clusters,
		const simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>& prior,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>::prior_state_changes& old_prior_changes,
		typename simple_hol_term_distribution<BuiltInPredicates, ConstantDistribution, SetSizeDistribution, DefinitionCountDistribution>::prior_state_changes& new_prior_changes)
{
	double value = 0.0;
	for (hol_term* formula : new_clusters)
		value += log_probability_helper(formula, prior, new_prior_changes);
	for (hol_term* formula : old_clusters)
		value -= log_probability_helper(formula, prior, old_prior_changes);

	/* get the name events */
	array_map<unsigned int, array<const hol_term*>> old_name_map(max(1, old_prior_changes.arg2_string_map.size));
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
		const hol_term* name = prior_state.arg2_string_map.get(entry.key, contains);
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

	array_map<unsigned int, array<const hol_term*>> new_name_map(max(1, new_prior_changes.arg2_string_map.size));
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
		const hol_term* name = prior_state.arg2_string_map.get(entry.key, contains);
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
	array_map<unsigned int, array<const hol_term*>> name_map(max(1, prior_state.arg2_string_map.size));
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
			array<const hol_term*>& names = name_map.get(old_name_map.keys[i], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count + ((ssize_t) new_name_map.values[j].length - old_name_map.values[i].length), prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			i++; j++;
		} else if (old_name_map.keys[i] < new_name_map.keys[j]) {
			array<const hol_term*>& names = name_map.get(old_name_map.keys[i], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count - ((ssize_t) old_name_map.values[i].length), prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			i++;
		} else {
			array<const hol_term*>& names = name_map.get(new_name_map.keys[j], contains);
			unsigned int old_count = (contains ? names.length : 0);
			name_prior += log_probability(old_count + new_name_map.values[j].length, prior.name_count_distribution)
						- log_probability(old_count, prior.name_count_distribution);
			j++;
		}
	} while (i < old_name_map.size) {
		array<const hol_term*>& names = name_map.get(old_name_map.keys[i], contains);
		unsigned int old_count = (contains ? names.length : 0);
		name_prior += log_probability(old_count - ((ssize_t) old_name_map.values[i].length), prior.name_count_distribution)
					- log_probability(old_count, prior.name_count_distribution);
		i++;
	} while (j < new_name_map.size) {
		array<const hol_term*>& names = name_map.get(new_name_map.keys[j], contains);
		unsigned int old_count = (contains ? names.length : 0);
		name_prior += log_probability(old_count + new_name_map.values[j].length, prior.name_count_distribution)
					- log_probability(old_count, prior.name_count_distribution);
		j++;
	}
	for (auto entry : old_name_map) free(entry.value);
	for (auto entry : new_name_map) free(entry.value);
	for (auto entry : name_map) free(entry.value);
	value += name_prior;

	double definition_prior = 0.0;
	for (const auto& entry : new_prior_changes.definitions) {
		const array<hol_term>& definitions = prior_state.definitions.get(entry.key, contains);
		unsigned int old_count = (contains ? definitions.length : 0);
		unsigned int index = old_prior_changes.definitions.index_of(entry.key);
		if (index < old_prior_changes.definitions.size) {
			definition_prior += log_probability(old_count + entry.value.length - old_prior_changes.definitions.values[index].length, prior.definition_count_distribution)
							  - log_probability(old_count, prior.definition_count_distribution);
		} else {
			definition_prior += log_probability(old_count + entry.value.length, prior.definition_count_distribution)
							  - log_probability(old_count, prior.definition_count_distribution);
		}
	} for (const auto& entry : old_prior_changes.definitions) {
		if (new_prior_changes.definitions.contains(entry.key)) continue;
		unsigned int old_count = prior_state.definitions.get(entry.key).length;
		definition_prior += log_probability(old_count - entry.value.length, prior.definition_count_distribution)
						  - log_probability(old_count, prior.definition_count_distribution);
	}
	value += definition_prior;

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

#endif /* THEORY_PRIOR_H_ */
