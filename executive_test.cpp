#include "higher_order_logic.h"
#include "hdp_parser.h"
#include "executive.h"
#include "ruletaker.h"
#include "geoquery.h"
#include "console.h"
#include "theory_prior.h"

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
	typedef typename Formula::TermType TermType;

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
			const pair<array<instance>, array<instance>>& types = T.atoms.get(*atom, contains);
			free(*atom); free(atom);
			if (!contains) return false;
			if (arg1->type == TermType::CONSTANT)
				return index_of_constant(types.key, arg1->constant) < types.key.length;
			else return index_of_number(types.key, arg1->number) < types.key.length;
		} else {
			/* `formula` is an atom of form `f(a,b)` */
			relation rel = { 0, arg1->constant, arg2->constant };
			const pair<array<instance>, array<instance>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return index_of_constant(relations.key, predicate) < relations.key.length;
		}
	} else if (formula->type == FormulaType::NOT && is_atomic(*formula->unary.operand)) {
		if (arg2 == NULL) {
			/* `formula` is an atom of form `~f(a)` */
			Term* atom = Term::new_apply(Term::new_constant(predicate), &Term::template variables<0>::value);
			if (atom == nullptr) return false;
			Term::template variables<0>::value.reference_count++;
			const pair<array<instance>, array<instance>>& types = T.atoms.get(*atom, contains);
			free(*atom); free(atom);
			if (!contains) return false;
			if (arg1->type == TermType::CONSTANT)
				return index_of_constant(types.value, arg1->constant) < types.value.length;
			else return index_of_number(types.value, arg1->number) < types.value.length;
		} else {
			/* `formula` is an atom of form `~f(a,b)` */
			relation rel = { 0, arg1->constant, arg2->constant };
			const pair<array<instance>, array<instance>>& relations = T.relations.get(rel, contains);
			if (!contains) return false;
			return index_of_constant(relations.value, predicate) < relations.value.length;
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

//run_console(stdin, "\nEnter command: ", parser, names);

//run_console(stdin, "\nEnter sentence to parse: ", parser, names, seed_training_set);
/*for (array_map<sentence_type, flagged_logical_form<hol_term>>& paragraph : seed_training_set) {
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

	/* construct the theory */
	typedef theory<natural_deduction<hol_term, false>, polymorphic_canonicalizer<true, false, built_in_predicates>> Theory;
	Theory T(1000000000);
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
	typedef decltype(proof_prior)::PriorState PriorStateType;
	PriorStateType proof_axioms;
	if (!parser.invert_name_map(names)) {
		for (auto entry : names) free(entry.key);
		return EXIT_FAILURE;
	}

/* run RuleTaker experiments */
/*run_ruletaker_experiments_single_threaded(corpus, parser, T, proof_axioms, proof_prior, names, seed_entities, "proofwriter/OWA/birds-electricity/meta-test.jsonl", "ruletaker_results.txt");
for (auto entry : names) free(entry.key);
// to avoid breakpoints being moved due to eliminated code
if (seed_training_set.length > 0)
return EXIT_SUCCESS;*/

Theory& T_copy = *((Theory*) alloca(sizeof(Theory)));
PriorStateType& proof_axioms_copy = *((PriorStateType*) alloca(sizeof(PriorStateType)));
hash_map<const hol_term*, hol_term*> formula_map(128);
Theory::clone(T, T_copy, formula_map);
PriorStateType::clone(proof_axioms, proof_axioms_copy, formula_map);

/* read the GeoBase sentences */
/*in = fopen("geobase_simple.txt", "rb");
if (in == nullptr) {
	fprintf(stderr, "ERROR: Unable to open file for reading.\n");
	for (auto entry : names) free(entry.key);
	return EXIT_FAILURE;
}
array<char> line(1024);
for (unsigned int counter = 0; ; counter++) {
	line.clear();
	int num_read = read_line(line, in);
	if (num_read == 0) {
		break;
	} else if (num_read > 0) {
		if (!line.ensure_capacity(line.length + 1))
			break;
		line[line.length++] = '\0';
		if (!read_sentence(corpus, parser, line.data, T_copy, names, seed_entities, proof_prior, proof_axioms_copy)) {
			fprintf(stderr, "ERROR: Unable to read sentence '%s'.\n", line.data);
			for (auto entry : names) free(entry.key);
			return EXIT_FAILURE;
		}

		Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
		PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
		formula_map.clear();
		Theory::clone(T_copy, T_MAP, formula_map);
		PriorStateType::clone(proof_axioms_copy, proof_axioms_MAP, formula_map);
		auto collector = make_log_probability_collector(T_copy, proof_prior);
		double max_log_probability = collector.current_log_probability;
		for (unsigned int j = 0; j < 4; j++) {
			for (unsigned int t = 0; t < 100; t++) {
				fprintf(stderr, "j = %u, t = %u\n", j, t);
				bool print_debug = false;
				if (print_debug) T_copy.template print_axioms<true>(stdout, *debug_terminal_printer);
				if (print_debug) { T_copy.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout); }
				do_mh_step(T_copy, proof_prior, proof_axioms_copy, collector);

				if (collector.current_log_probability > max_log_probability) {
					free(T_MAP); free(proof_axioms_MAP); formula_map.clear();
					Theory::clone(T_copy, T_MAP, formula_map);
					PriorStateType::clone(proof_axioms_copy, proof_axioms_MAP, formula_map);
					max_log_probability = collector.current_log_probability;
				}
			}

			if (j + 1 < 4) {
				for (unsigned int t = 0; t < 20; t++)
					do_exploratory_mh_step(T_copy, proof_prior, proof_axioms_copy, collector);
			}
		}
		free(T_copy); free(proof_axioms_copy); formula_map.clear();
		Theory::clone(T_MAP, T_copy, formula_map);
		PriorStateType::clone(proof_axioms_MAP, proof_axioms_copy, formula_map);
		T_MAP.template print_axioms<true>(stdout, *debug_terminal_printer); fflush(stdout);
		free(T_MAP); free(proof_axioms_MAP);
	}

	FILE* theory_stream = (FILE*) fopen("geoquery.theory", "wb");
	write_random_state(theory_stream);
	write(T_copy, theory_stream, proof_axioms_copy);
	fclose(theory_stream);
}
fclose(in);*/

free(T_copy); free(proof_axioms_copy);
FILE* theory_stream = (FILE*) fopen("geoquery.theory", "rb");
read_random_state(theory_stream);
read(T_copy, theory_stream, proof_axioms_copy);
fclose(theory_stream);

/* run GeoQuery experiments */
run_geoquery_experiments_single_threaded(corpus, parser, T_copy, proof_axioms_copy, proof_prior, names, seed_entities, "geoquery-rephrased.jsonl", "geoquery_results.txt");
free(T_copy); free(proof_axioms_copy);
for (auto entry : names) free(entry.key);
// to avoid breakpoints being moved due to eliminated code
if (seed_training_set.length > 0)
return EXIT_SUCCESS;

	/*read_sentence(corpus, parser, "A butterfly has a spot.", T, names, seed_entities, proof_prior, proof_axioms);

	parse_sentence_with_prior(parser, "Sally caught a butterfly with a spot.", T, names, proof_prior, proof_axioms);
for (auto entry : names) free(entry.key);
return EXIT_SUCCESS;*/

	/* read the articles */
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
hash_map<const hol_term*, hol_term*> term_map(128);
TheoryType::clone(T, T_map, term_map);
auto collector = make_log_probability_collector(T, proof_prior);
double max_log_probability = collector.current_log_probability;
for (unsigned int t = 0; t < 200 /*20000*/; t++) {
	bool print_debug = false;
	if (print_debug) T.print_axioms(stderr, *debug_terminal_printer);
	if (print_debug) T.print_disjunction_introductions(stderr, *debug_terminal_printer);
	do_mh_step(T, proof_prior, proof_axioms, collector);
	if (collector.current_log_probability > max_log_probability) {
		free(T_map); term_map.clear();
		TheoryType::clone(T, T_map, term_map);
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
