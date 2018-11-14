#include "executive.h"
#include "fake_parser.h"

struct poisson_distribution {
	double lambda, log_lambda;

	poisson_distribution(double lambda) : lambda(lambda), log_lambda(log(lambda)) { }
	poisson_distribution(const poisson_distribution& src) : lambda(src.lambda), log_lambda(src.log_lambda) { }
};

inline double log_probability(unsigned int k, const poisson_distribution& prior) {
	return k * prior.log_lambda - prior.lambda - lgamma(k + 1);
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
			unsigned int count = tables.counts.get(get_key(old_observations, i));
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
struct constant_adder {
	ConstantDistribution& constant_distribution;
};

template<typename ConstantDistribution>
inline bool visit_constant(unsigned int constant, const constant_adder<ConstantDistribution>& visitor) {
	return visitor.constant_distribution.add_constant(constant);
}

template<typename ConstantDistribution>
inline bool visit_predicate(unsigned int predicate, const constant_adder<ConstantDistribution>& visitor) {
	return visitor.constant_distribution.add_predicate(predicate);
}

template<typename ConstantDistribution> constexpr bool visit_variable(unsigned int variable, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_parameter(unsigned int parameter, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_true(const fol_formula& formula, const constant_adder<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_false(const fol_formula& formula, const constant_adder<ConstantDistribution>& visitor) { return true; }

template<fol_formula_type Operator, typename ConstantDistribution>
constexpr bool visit_operator(const fol_formula& formula, const constant_adder<ConstantDistribution>& visitor) { return true; }

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
	visitor.constant_distribution.remove_predicate(predicate);
	return true;
}

template<typename ConstantDistribution> constexpr bool visit_variable(unsigned int variable, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_parameter(unsigned int parameter, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_true(const fol_formula& formula, const constant_remover<ConstantDistribution>& visitor) { return true; }
template<typename ConstantDistribution> constexpr bool visit_false(const fol_formula& formula, const constant_remover<ConstantDistribution>& visitor) { return true; }

template<fol_formula_type Operator, typename ConstantDistribution>
constexpr bool visit_operator(const fol_formula& formula, const constant_remover<ConstantDistribution>& visitor) { return true; }

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
		constants.add_constant(atom.arg1.constant);
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
		const ObservationCollection<fol_formula*>& observations,
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
	if (size(observation.key) < 2)
		return -std::numeric_limits<double>::infinity();
	log_cache<double>::instance().ensure_size(size(observation.value));
	return -log_cache<double>::instance().get(size(observation.value)) * size(observation.key)
		 + (size(observation.key) - 2) * prior.log_continue_probability + prior.log_stop_probability;
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
struct uniform_distribution {
	typedef default_array<pair<T, array<T>>> ObservationCollection;
};

template<template<typename> class Collection, typename T>
double log_probability(
		const pair<T, Collection<T>>& observation,
		const uniform_distribution<T>& prior)
{
	log_cache<double>::instance().ensure_size(size(observation.value));
	return -log_cache<double>::instance().get(size(observation.value));
}

template<template<typename> class ObservationCollection,
	template<typename> class Collection, typename T>
double log_probability(
		const ObservationCollection<pair<T, Collection<T>>>& observations,
		const uniform_distribution<T>& prior)
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
		const uniform_distribution<T>& prior)
{
	return log_probability(new_observations, prior) - log_probability(old_observations, prior);
}

const string* get_name(const hash_map<string, unsigned int>& names, unsigned int id)
{
	for (const auto& entry : names)
		if (entry.value == id) return &entry.key;
	return NULL;
}

template<typename Formula, typename ProofCalculus>
bool contains_axiom(const theory<Formula, ProofCalculus>& T, const Formula* formula)
{
	typedef typename Formula::Type FormulaType;
	typedef typename Formula::TermType TermType;
	typedef typename ProofCalculus::Proof Proof;

	if (formula == NULL) {
		fprintf(stderr, "contains_axiom WARNING: `formula` is null.\n");
		return false;
	}

	bool contains;
	if (formula->type == FormulaType::FOR_ALL) {
		for (const Proof* axiom : T.universal_quantifications)
			if (*axiom->formula == *formula) return true;
		return false;
	} else if (formula->type == FormulaType::ATOM) {
		if (formula->atom.arg2.type == TermType::NONE) {
			/* `formula` is an atom of form `f(a)` */
			const pair<array<unsigned int>, array<unsigned int>>& types = T.types.get(formula->atom.predicate, contains);
			if (!contains) return false;
			return types.key.contains(formula->atom.arg1.constant);
		} else {
			/* `formula` is an atom of form `f(a,b)` */
			relation rel = { 0, formula->atom.arg1.constant, formula->atom.arg2.constant };
			const pair<array<unsigned int>, array<unsigned int>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return relations.key.contains(formula->atom.predicate);
		}
	} else if (formula->type == FormulaType::NOT && formula->unary.operand->type == FormulaType::ATOM) {
		formula = formula->unary.operand;
		if (formula->atom.arg2.type == TermType::NONE) {
			/* `formula` is an atom of form `~f(a)` */
			const pair<array<unsigned int>, array<unsigned int>>& types = T.types.get(formula->atom.predicate, contains);
			if (!contains) return false;
			return types.value.contains(formula->atom.arg1.constant);
		} else {
			/* `formula` is an atom of form `~f(a,b)` */
			relation rel = { 0, formula->atom.arg1.constant, formula->atom.arg2.constant };
			const pair<array<unsigned int>, array<unsigned int>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return relations.value.contains(formula->atom.predicate);
		}
	} else {
		return false;
	}
}

int main(int argc, const char** argv)
{
	FILE* in = fopen("flat_ontology_articles.txt", "r");
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

	hash_map<string, unsigned int> names(256);
	if (!add_constants_to_string_map(names)) {
		free_tokens(tokens); return EXIT_FAILURE;
	}

	unsigned int index = 0;
	in_memory_article_store corpus;
	fake_parser parser = fake_parser(PREDICATE_UNKNOWN);
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
	theory<fol_formula, natural_deduction<fol_formula, true>> T(names.table.size + 1);
	auto constant_prior = make_simple_constant_distribution(
			chinese_restaurant_process<unsigned int>(1.0), chinese_restaurant_process<unsigned int>(1.0));
	auto theory_element_prior = make_simple_fol_formula_distribution(constant_prior, 0.01, 0.3, 0.4, 0.2, 0.4);
	auto axiom_prior = make_dirichlet_process(1.0e-4, theory_element_prior);
	auto conjunction_prior = uniform_subset_distribution<const nd_step<fol_formula, true>*>(0.1);
	auto universal_introduction_prior = uniform_distribution<unsigned int>();
	auto universal_elimination_prior = chinese_restaurant_process<fol_term>(1.0);
	auto proof_prior = make_canonicalized_proof_prior(axiom_prior, conjunction_prior,
			universal_introduction_prior, universal_elimination_prior, poisson_distribution(1.0));
	const string** reverse_name_map = invert(names);
	string_map_scribe printer = { reverse_name_map, names.table.size + 1 };
	read_article(names.get("Bob"), corpus, parser, T, proof_prior, printer);
	read_article(names.get("Kate"), corpus, parser, T, proof_prior, printer);
	read_article(names.get("Sam"), corpus, parser, T, proof_prior, printer);
	read_article(names.get("Byron"), corpus, parser, T, proof_prior, printer);
	read_article(names.get("Alex"), corpus, parser, T, proof_prior, printer);

	fol_formula* all_cats_are_mammals = fol_formula::new_for_all(1,
			fol_formula::new_if_then(
				fol_formula::new_atom(parser.symbol_map.get(names.get("cat")), fol_formula::new_variable(1)),
				fol_formula::new_atom(parser.symbol_map.get(names.get("mammal")), fol_formula::new_variable(1))
			)
		);
	fol_formula* all_mammals_are_cats = fol_formula::new_for_all(1,
			fol_formula::new_if_then(
				fol_formula::new_atom(parser.symbol_map.get(names.get("mammal")), fol_formula::new_variable(1)),
				fol_formula::new_atom(parser.symbol_map.get(names.get("cat")), fol_formula::new_variable(1))
			)
		);

	unsigned int all_cats_are_mammals_count = 0;
	unsigned int all_mammals_are_cats_count = 0;
	unsigned int neither_count = 0;
	constexpr unsigned int iterations = 10000000;
	timer stopwatch;
	for (unsigned int t = 0; t < iterations; t++) {
		//T.print_axioms(stdout, parser.get_printer(printer));
		//print('\n', stdout); fflush(stdout);
		if (stopwatch.milliseconds() > 1000) {
			fprintf(stderr, "all_cats_are_mammals: %lf\n", (double) all_cats_are_mammals_count / t);
			fprintf(stderr, "all_mammals_are_cats: %lf\n", (double) all_mammals_are_cats_count / t);
			stopwatch.start();
		}
		do_mh_step(T, proof_prior);

		bool has_all_cats_are_mammals = contains_axiom(T, all_cats_are_mammals);
		bool has_all_mammals_are_cats = contains_axiom(T, all_mammals_are_cats);
		if (has_all_cats_are_mammals) all_cats_are_mammals_count++;
		if (has_all_mammals_are_cats) all_mammals_are_cats_count++;
		if (!has_all_cats_are_mammals && !has_all_mammals_are_cats) neither_count++;
	}

	fprintf(stderr, "all_cats_are_mammals: %lf\n", (double) all_cats_are_mammals_count / iterations);
	fprintf(stderr, "all_mammals_are_cats: %lf\n", (double) all_mammals_are_cats_count / iterations);

	free(*all_cats_are_mammals); free(all_cats_are_mammals);
	free(*all_mammals_are_cats); free(all_mammals_are_cats);
	free(reverse_name_map);
	for (auto entry : names) free(entry.key);
	return EXIT_SUCCESS;
}
