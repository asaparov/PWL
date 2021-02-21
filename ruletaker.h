/**
 * ruletaker.h - Code to read RukeTaker data and attempt to answer the questions.
 *
 *  Created on: Dec 8, 2020
 *      Author: asaparov
 */

#ifndef RULETAKER_H_
#define RULETAKER_H_

#include <stdio.h>

#include "json.h"

enum class ruletaker_key {
	THEORY,
	QUESTIONS,
	QUESTION,
	ANSWER,
	OTHER
};

enum class ruletaker_reader_state {
	START,
	ROOT_OBJECT,
	QUESTIONS,
	QUESTION,
	TRIPLES,
	TRIPLE,
	PROOFS_WITH_INTERMEDIATES,
	INTERMEDIATES,
	INTERMEDIATE
};

enum class ruletaker_label : uint_fast8_t {
	FALSE,
	TRUE,
	UNKNOWN
};

struct ruletaker_reader {
	ruletaker_reader_state state;
	ruletaker_key current_key;
	char* context;
	pair<char*, ruletaker_label> next_question;
	array<pair<string, ruletaker_label>> questions;

	ruletaker_reader() : state(ruletaker_reader_state::START), context(nullptr), next_question(nullptr, ruletaker_label::UNKNOWN), questions(44) { }

	~ruletaker_reader() {
		if (context != nullptr)
			free(context);
		if (next_question.key != nullptr)
			free(next_question.key);
		for (pair<string, ruletaker_label>& question : questions)
			free(question.key);
	}
};

inline bool begin_object(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::START) {
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	} else if (reader.current_key == ruletaker_key::QUESTIONS) {
		reader.state = ruletaker_reader_state::QUESTIONS;
	} else if (reader.state == ruletaker_reader_state::QUESTIONS) {
		reader.state = ruletaker_reader_state::QUESTION;
	} else if (reader.state == ruletaker_reader_state::QUESTION) {
		reader.state = ruletaker_reader_state::PROOFS_WITH_INTERMEDIATES;
	} else if (reader.state == ruletaker_reader_state::PROOFS_WITH_INTERMEDIATES) {
		reader.state = ruletaker_reader_state::INTERMEDIATES;
	} else if (reader.state == ruletaker_reader_state::INTERMEDIATES) {
		reader.state = ruletaker_reader_state::INTERMEDIATE;
	} else if (reader.state == ruletaker_reader_state::ROOT_OBJECT) {
		reader.state = ruletaker_reader_state::TRIPLES;
	} else if (reader.state == ruletaker_reader_state::TRIPLES) {
		reader.state = ruletaker_reader_state::TRIPLE;
	}
	return true;
}

inline bool end_object(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::ROOT_OBJECT) {
		reader.state = ruletaker_reader_state::START;
	} else if (reader.state == ruletaker_reader_state::QUESTION) {
		if (!reader.questions.ensure_capacity(reader.questions.length + 1))
			return false;
		reader.questions[reader.questions.length].key.data = reader.next_question.key;
		reader.questions[reader.questions.length].key.length = strlen(reader.next_question.key);
		reader.questions[reader.questions.length++].value = reader.next_question.value;
		reader.next_question.key = nullptr;
		reader.state = ruletaker_reader_state::QUESTIONS;
	} else if (reader.state == ruletaker_reader_state::QUESTIONS) {
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	} else if (reader.state == ruletaker_reader_state::PROOFS_WITH_INTERMEDIATES) {
		reader.state = ruletaker_reader_state::QUESTION;
	} else if (reader.state == ruletaker_reader_state::INTERMEDIATES) {
		reader.state = ruletaker_reader_state::PROOFS_WITH_INTERMEDIATES;
	} else if (reader.state == ruletaker_reader_state::INTERMEDIATE) {
		reader.state = ruletaker_reader_state::INTERMEDIATES;
	} else if (reader.state == ruletaker_reader_state::TRIPLES) {
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	} else if (reader.state == ruletaker_reader_state::TRIPLE) {
		reader.state = ruletaker_reader_state::TRIPLES;
	}
	return true;
}

inline bool emit_key(const array<char>& key_name, const position& pos, ruletaker_reader& reader) {
	if (compare_strings(key_name, "theory")) {
		reader.current_key = ruletaker_key::THEORY;
	} else if (compare_strings(key_name, "questions")) {
		reader.current_key = ruletaker_key::QUESTIONS;
	} else if (compare_strings(key_name, "question")) {
		reader.current_key = ruletaker_key::QUESTION;
	} else if (compare_strings(key_name, "answer")) {
		reader.current_key = ruletaker_key::ANSWER;
	} else {
		reader.current_key = ruletaker_key::OTHER;
	}
	return true;
}

inline constexpr bool begin_list(const position& pos, ruletaker_reader& reader) { return true; }
inline constexpr bool end_list(const position& pos, ruletaker_reader& reader) { return true; }

inline bool emit_string(const array<char>& str, const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::THEORY) {
		if (reader.context != nullptr) {
			read_error("Found duplicate `theory` entry", pos);
			return false;
		}
		reader.context = (char*) malloc(sizeof(char) * (str.length + 1));
		if (reader.context == nullptr) {
			fprintf(stderr, "emit_string ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str.length; i++)
			reader.context[i] = str[i];
		reader.context[str.length] = '\0';
	} else if (reader.current_key == ruletaker_key::QUESTION) {
		if (reader.next_question.key != nullptr) {
			read_error("Found duplicate `question` entry", pos);
			return false;
		}
		reader.next_question.key = (char*) malloc(sizeof(char) * (str.length + 1));
		if (reader.next_question.key == nullptr) {
			fprintf(stderr, "emit_string ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str.length; i++)
			reader.next_question.key[i] = str[i];
		reader.next_question.key[str.length] = '\0';
	} else if (reader.current_key == ruletaker_key::ANSWER) {
		if (compare_strings(str, "Unknown")) {
			reader.next_question.value = ruletaker_label::UNKNOWN;
		} else {
			read_error("Label must be either `true`, `false`, or \"Unknown\"", pos);
			return false;
		}
	}
	return true;
}

inline bool emit_true(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::ANSWER)
		reader.next_question.value = ruletaker_label::TRUE;
	return true;
}

inline bool emit_false(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::ANSWER)
		reader.next_question.value = ruletaker_label::FALSE;
	return true;
}

constexpr bool emit_null(const position& pos, ruletaker_reader& reader) { return true; }
constexpr bool emit_number(double value, const position& pos, ruletaker_reader& reader) { return true; }

template<typename Stream>
struct line_reader {
	Stream& in;
	bool eof;

	line_reader(Stream& in) : in(in), eof(false) { }
};

template<typename Stream>
int fgetc(line_reader<Stream>& lr) {
	int c = fgetc(lr.in);
	if (c == -1)
		lr.eof = true;
	else if (c == '\n')
		return -1;
	return c;
}

template<typename Stream>
int ungetc(int c, line_reader<Stream>& lr) {
	return ungetc(c, lr.in);
}

template<typename ProcessQuestions>
bool read_ruletaker_data(const char* filename, ProcessQuestions process_questions)
{
	FILE* in = (FILE*) fopen(filename, "rb");
	if (in == nullptr) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	line_reader<FILE*> lr(in);
	position current(1, 1);
	while (!lr.eof) {
		/* check if this is an empty line */
		int c = fgetc(lr);
		if (c == -1) continue;
		ungetc(c, lr);

		ruletaker_reader reader;
		if (!json_parse(lr, reader, current)) {
			fclose(in);
			return false;
		}
		current.line++;
		current.column = 1;

		if (!process_questions(reader.context, reader.questions)) {
			fclose(in);
			return false;
		}
	}

	fclose(in);
	return true;
}

#include <atomic>

enum class ruletaker_work_item_type {
	READ_CONTEXT,
	ANSWER_QUESTION
};

template<typename Theory, typename PriorStateType>
struct ruletaker_context_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	std::minstd_rand prng_engine;
	char* context;
	array<pair<string, ruletaker_label>> questions;

	static inline void free(ruletaker_context_item<Theory, PriorStateType>& item) {
		core::free(item.context);
		for (auto& entry : item.questions)
			core::free(entry.key);
		core::free(item.questions);
		if (!core::is_empty(item.T)) {
			core::free(item.T);
			item.proof_axioms.~PriorStateType();
		}
	}
};

template<typename Theory, typename PriorStateType>
struct ruletaker_question_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	unsigned int question_id;
	string question;
	ruletaker_label label;

	static inline void free(ruletaker_question_item& item) {
		core::free(item.question);
		if (!core::is_empty(item.T)) {
			core::free(item.T);
			item.proof_axioms.~PriorStateType();
		}
	}
};

template<typename BuiltInPredicates>
inline void find_head_or_not(
		hol_term* src, hol_term*& head,
		head_index& predicate_index)
{
	if (src->type == hol_term_type::NOT) {
		head = src;
		return;
	}
	find_head<BuiltInPredicates>(src, head, predicate_index);
}

inline bool negate_head(
		hol_term* src, hol_term*& dst)
{
	array<hol_term*> scopes(8);
	auto gather_scopes = [&scopes](hol_term* term) {
		if (term->type == hol_term_type::EXISTS || term->type == hol_term_type::FOR_ALL || term->type == hol_term_type::LAMBDA)
			return scopes.add(term);
		return true;
	};

	head_index predicate_index;
	hol_term* head = find_head(src, predicate_index, find_head_or_not<built_in_predicates>, gather_scopes);
	if (head == nullptr)
		return false;

	hol_term* new_head;
	if (head->type == hol_term_type::NOT) {
		new_head = head->unary.operand;
		new_head->reference_count++;
	} else {
		/* check if the head is a `same` event */
		bool is_same = true;
		unsigned int head_index = scopes.length - 1;
		do {
			is_same = false;
			if (head->type == hol_term_type::EXISTS) {
				hol_term* operand = head->quantifier.operand;
				if (operand->type == hol_term_type::AND) {
					for (unsigned int i = 0; !is_same && i < operand->array.length; i++) {
						hol_term* conjunct = operand->array.operands[i];
						if (conjunct->type == hol_term_type::UNARY_APPLICATION
						 && (*conjunct->binary.left == hol_term::constants<(unsigned int) built_in_predicates::SAME>::value
						  || *conjunct->binary.left == hol_term::constants<(unsigned int) built_in_predicates::NAME>::value)
						 && conjunct->binary.right->type == hol_term_type::VARIABLE && conjunct->binary.right->variable == head->quantifier.variable)
						{
							is_same = true;
						} else if (get_scope<built_in_predicates, (unsigned int) built_in_predicates::NAME>(conjunct, head->quantifier.variable) != nullptr) {
							is_same = true;
						}
					}
				} else {
					if (operand->type == hol_term_type::UNARY_APPLICATION
					 && (*operand->binary.left == hol_term::constants<(unsigned int) built_in_predicates::SAME>::value
					  || *operand->binary.left == hol_term::constants<(unsigned int) built_in_predicates::NAME>::value)
					 && operand->binary.right->type == hol_term_type::VARIABLE && operand->binary.right->variable == head->quantifier.variable)
					{
						is_same = true;
					} else if (get_scope<built_in_predicates, (unsigned int) built_in_predicates::NAME>(operand, head->quantifier.variable) != nullptr) {
						is_same = true;
					}
				}

				if (is_same) {
					if (head_index > 0) {
						/* negate the scope immediately preceding `head` */
						head_index--;
						head = scopes[head_index];
					} else {
						fprintf(stderr, "negate_head ERROR: Unable to find appropriate scope to negate.\n");
						return false;
					}
				}
			}
		} while (is_same);

		new_head = hol_term::new_not(head);
		if (new_head == nullptr) return false;
		head->reference_count++;
	}

	dst = substitute_head<any_node_position::NONE>(src, head, new_head);
	free(*new_head); if (new_head->reference_count == 0) free(new_head);
	return (dst != nullptr);
}

template<typename Proof>
struct dummy_collector {
	const Proof* test_proof;

	dummy_collector(const Proof* test_proof) : test_proof(test_proof) { }

	constexpr inline bool has_prior(const Proof* proof) const {
		return (proof != test_proof);
	}
};

template<typename Proof>
dummy_collector<Proof> make_dummy_collector(const Proof* test_proof) {
	return dummy_collector<Proof>(test_proof);
}

template<typename Theory, typename ProofPrior>
inline void print_theory(const Theory& T, const typename Theory::Proof* test_proof, ProofPrior& proof_prior)
{
	typedef typename Theory::Formula Formula;

	array<Formula*> extra_axioms(16);
	T.get_extra_axioms(extra_axioms);
	auto collector = make_dummy_collector(test_proof);
	double value = log_probability(T.observations, extra_axioms, proof_prior, collector);
	fprintf(stderr, "log probability of theory: %lf\n", value);
	T.print_axioms(stderr, *debug_terminal_printer);
	T.print_disjunction_introductions(stderr, *debug_terminal_printer);
}

struct question_result {
	unsigned int context_id;
	unsigned int question_id;
	double log_probability_diff;
	ruletaker_label true_label;

	static inline void swap(question_result& first, question_result& second) {
		core::swap(first.context_id, second.context_id);
		core::swap(first.question_id, second.question_id);
		core::swap(first.log_probability_diff, second.log_probability_diff);
		core::swap(first.true_label, second.true_label);
	}
};

inline bool operator < (const question_result& first, const question_result& second) {
	if (first.context_id < second.context_id) return true;
	else if (first.context_id > second.context_id) return false;
	else return first.question_id < second.question_id;
}

inline bool operator == (const question_result& first, const question_result& second) {
	return first.context_id == second.context_id
		&& first.question_id == second.question_id;
}

constexpr unsigned int MAX_CONTEXT_COUNT = 140;
constexpr unsigned int MAX_QUESTION_COUNT = 5270;
constexpr double PREDICT_UNKNOWN_THRESHOLD = 2000.0;


template<typename ProofCalculus, typename Canonicalizer>
inline bool is_formula_possible(
		theory<ProofCalculus, Canonicalizer>& T,
		typename ProofCalculus::Language* logical_form)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	unsigned int new_constant;
	set_changes<Formula> set_diff;
	Proof* new_proof = T.add_formula(logical_form, set_diff, new_constant);
	if (new_proof == nullptr)
		return false;
	T.template remove_formula<true>(new_proof, set_diff);
	return true;
}

template<typename Theory, typename PriorStateType, typename ProofPrior>
bool resample_observations(Theory& T, PriorStateType& proof_axioms, ProofPrior& proof_prior)
{
	typedef typename Theory::Proof Proof;
	typedef typename Theory::ProofType ProofType;

	array<pair<hol_term*, Proof*>> observations(T.observations.size);
	for (Proof* proof : T.observations) {
		/* get the logical form conclusion of this proof */
		hol_term* formula;
		if (proof->type == ProofType::AXIOM) {
			formula = proof->formula;
		} else if (proof->type == ProofType::EXISTENTIAL_INTRODUCTION) {
			unsigned int index;
			for (index = 0; index < T.existential_intro_nodes.length; index++)
				if (T.existential_intro_nodes[index].value == proof) break;
			if (index == T.existential_intro_nodes.length) {
				fprintf(stderr, "resample_observations ERROR: Unable to find proof in `theory.existential_intro_nodes`.\n");
				return false;
			}
			formula = T.existential_intro_nodes[index].key;
		} else if (proof->type == ProofType::IMPLICATION_INTRODUCTION) {
			unsigned int index;
			for (index = 0; index < T.implication_intro_nodes.length; index++)
				if (T.implication_intro_nodes[index].value == proof) break;
			if (index == T.implication_intro_nodes.length) {
				fprintf(stderr, "resample_observations ERROR: Unable to find proof in `theory.implication_intro_nodes`.\n");
				return false;
			}
			formula = T.implication_intro_nodes[index].key;
		} else if (proof->type == ProofType::DISJUNCTION_INTRODUCTION) {
			unsigned int index;
			for (index = 0; index < T.disjunction_intro_nodes.length; index++)
				if (T.disjunction_intro_nodes[index].value == proof) break;
			if (index == T.disjunction_intro_nodes.length) {
				fprintf(stderr, "resample_observations ERROR: Unable to find proof in `theory.disjunction_intro_nodes`.\n");
				return false;
			}
			formula = T.disjunction_intro_nodes[index].key;
		} else if (proof->type == ProofType::PROOF_BY_CONTRADICTION
				&& proof->operands[0]->type == ProofType::NEGATION_ELIMINATION
				&& proof->operands[0]->operands[0]->type == ProofType::CONJUNCTION_ELIMINATION
				&& proof->operands[0]->operands[0]->operands[0] == proof->operands[1])
		{
			unsigned int index;
			for (index = 0; index < T.negated_conjunction_nodes.length; index++)
				if (T.negated_conjunction_nodes[index].value == proof) break;
			if (index == T.negated_conjunction_nodes.length) {
				fprintf(stderr, "resample_observations ERROR: Unable to find proof in `theory.negated_conjunction_nodes`.\n");
				return false;
			}
			formula = T.negated_conjunction_nodes[index].key;
		} else {
			fprintf(stderr, "resample_observations ERROR: Unsupported `ProofType`.\n");
			return false;
		}

		observations[observations.length++] = {formula, proof};
	}

	for (const pair<hol_term*, Proof*>& entry : observations) {
		entry.key->reference_count++;
		set_changes<hol_term> set_diff;
		T.template remove_formula<false>(entry.value, set_diff);
		proof_axioms.template subtract(entry.value, set_diff.old_set_axioms, proof_prior);
		free(*entry.value); if (entry.value->reference_count == 0) free(entry.value);
	}

	for (unsigned int i = observations.length - 1; i > 0; i--) {
		unsigned int next = sample_uniform(i + 1);
		if (next != i) {
			core::swap(observations[next].key, observations[i].key);
			core::swap(observations[next].value, observations[i].value);
		}
	}
	for (unsigned int i = 0; i < observations.length; i++) {
		unsigned int new_constant;
		set_changes<hol_term> set_diff;
		Proof* new_proof = T.add_formula(observations[i].key, set_diff, new_constant);
		while (new_proof == nullptr) {
			null_collector collector;
			set_diff.clear();
			for (unsigned int t = 0; t < 10; t++)
				do_exploratory_mh_step(T, proof_prior, proof_axioms, collector);
			new_proof = T.add_formula(observations[i].key, set_diff, new_constant);
		}
		if (new_proof == nullptr) {
			for (auto entry : observations) { free(*entry.key); if (entry.key->reference_count == 0) free(entry.key); }
			return false;
		} else if (!proof_axioms.template add(new_proof, set_diff.new_set_axioms, proof_prior)) {
			T.template remove_formula<true>(new_proof, set_diff);
			for (auto entry : observations) { free(*entry.key); if (entry.key->reference_count == 0) free(entry.key); }
			return false;
		}
		observations[i].value = new_proof;

		Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
		PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
		hash_map<const hol_term*, hol_term*> formula_map(128);
		Theory::clone(T, T_MAP, formula_map);
		new (&proof_axioms_MAP) PriorStateType(proof_axioms, formula_map);
		auto collector = make_log_probability_collector(T, proof_prior);
		double max_log_probability = collector.current_log_probability;
		for (unsigned int j = 0; j < 4; j++) {
			for (unsigned int t = 0; t < 400; t++) {
				do_mh_step(T, proof_prior, proof_axioms, collector);
				if (collector.current_log_probability > max_log_probability) {
					free(T_MAP); proof_axioms_MAP.~PriorStateType(); formula_map.clear();
					Theory::clone(T, T_MAP, formula_map);
					new (&proof_axioms_MAP) PriorStateType(proof_axioms, formula_map);
					max_log_probability = collector.current_log_probability;
				}
			}

			if (j + 1 < 4) {
				for (unsigned int t = 0; t < 20; t++)
					do_exploratory_mh_step(T, proof_prior, proof_axioms, collector);
			}
		}
		free(T); proof_axioms.~PriorStateType(); formula_map.clear();
		Theory::clone(T_MAP, T, formula_map);
		new (&proof_axioms) PriorStateType(proof_axioms_MAP, formula_map);
		T_MAP.print_axioms(stderr, *debug_terminal_printer);
		free(T_MAP); proof_axioms_MAP.~PriorStateType();
	}

	for (auto entry : observations) { free(*entry.key); if (entry.key->reference_count == 0) free(entry.key); }
	return true;
}

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
void do_ruletaker_experiments(bool& status,
		ruletaker_context_item<Theory, PriorStateType>* context_queue,
		ruletaker_question_item<Theory, PriorStateType>* question_queue,
		unsigned int& context_queue_start,
		unsigned int& question_queue_start,
		unsigned int& context_queue_length,
		unsigned int& question_queue_length,
		std::mutex& work_queue_lock,
		std::condition_variable& work_queue_cv,
		const std::minstd_rand prng_engine,
		ArticleSource& corpus, const Parser& parser_src,
		ProofPrior& proof_prior,
		const hash_map<string, unsigned int>& names_src,
		hash_set<unsigned int>& seed_entities,
		std::mutex& results_lock,
		array<question_result>& results,
		array<pair<unsigned int, string>>& unparseable_context,
		std::atomic_uint& total,
		std::atomic_uint& num_threads_reading_context,
		std::atomic_uint& num_threads_running)
{
	num_threads_running++;
	Parser& parser = *((Parser*) alloca(sizeof(Parser)));
	if (!init(parser, parser_src)) {
		status = false;
		num_threads_running--;
		work_queue_cv.notify_all();
		return;
	}
	hash_map<string, unsigned int> names(names_src.table.capacity);
	for (const auto& entry : names_src) {
		unsigned int index = names.table.index_to_insert(entry.key);
		if (!init(names.table.keys[index], entry.key)) {
			status = false;
			num_threads_running--;
			work_queue_cv.notify_all();
			for (auto entry : names) free(entry.key);
			free(parser); return;
		}
		names.values[index] = entry.value;
		names.table.size++;
	}
	if (!parser.invert_name_map(names)) {
		status = false;
		num_threads_running--;
		work_queue_cv.notify_all();
		for (auto entry : names) free(entry.key);
		free(parser); return;
	}

	while (status)
	{
		std::unique_lock<std::mutex> lock(work_queue_lock);
		while (status && context_queue_length == context_queue_start && question_queue_length == question_queue_start && num_threads_reading_context != 0)
			work_queue_cv.wait(lock);
		if (!status || (num_threads_reading_context == 0 && context_queue_start == context_queue_length && question_queue_start == question_queue_length)) {
			num_threads_running--;
			for (auto entry : names) free(entry.key);
			free(parser); return;
		}

		if (question_queue_start < question_queue_length) {
			ruletaker_question_item<Theory, PriorStateType>& job = question_queue[question_queue_start++];
			lock.unlock();
if (job.question_id < 23 - 1)
{
total++;
free(job);
continue;
}

			/* for reproducibility, reset the PRNG state */
			core::engine = context_queue[job.context_id].prng_engine;

			/* first clone the theory from the appropriate work item */
			Theory& T_copy = *((Theory*) alloca(sizeof(Theory)));
			hash_map<const hol_term*, hol_term*> formula_map(128);
			if (!Theory::clone(job.T, T_copy, formula_map)) {
				status = false;
				num_threads_running--;
				work_queue_cv.notify_all();
				free(job);
				for (auto entry : names) free(entry.key);
				free(parser); return;
			}
			PriorStateType proof_axioms_copy(job.proof_axioms, formula_map);

			unsigned int parse_count;
			constexpr unsigned int max_parse_count = 2;
			hol_term* logical_forms[max_parse_count];
			double log_probabilities[max_parse_count];
			if (parse_sentence(parser, job.question.data, names, logical_forms, log_probabilities, parse_count))
			{
				typedef typename Theory::Proof Proof;
				Theory& T_MAP_true = *((Theory*) alloca(sizeof(Theory)));
				Proof* proof_MAP_true; Proof* proof_MAP_false;
				double log_probability_true = log_joint_probability_of_truth(job.T, proof_prior, job.proof_axioms, logical_forms[0], 100, 4, 20, T_MAP_true, proof_MAP_true);
				for (unsigned int j = 0; isinf(log_probability_true) && j < 400; j++) {
					null_collector collector;
					for (unsigned int t = 0; t < 10; t++)
						do_exploratory_mh_step(job.T, proof_prior, job.proof_axioms, collector);
					log_probability_true = log_joint_probability_of_truth(job.T, proof_prior, job.proof_axioms, logical_forms[0], 100, 4, 20, T_MAP_true, proof_MAP_true);
				}

				if (!isinf(log_probability_true)) {
					hol_term* negated;
					if (!negate_head(logical_forms[0], negated) || negated == nullptr) {
						free_logical_forms(logical_forms, parse_count);
						status = false;
						num_threads_running--;
						work_queue_cv.notify_all();
						free(job); free(T_copy); free(T_MAP_true);
						total++;
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}

					/* for reproducibility, reset the PRNG state */
					core::engine = context_queue[job.context_id].prng_engine;

					Theory& T_MAP_false = *((Theory*) alloca(sizeof(Theory)));
T_copy.print_axioms(stderr, *debug_terminal_printer);
T_copy.print_disjunction_introductions(stderr, *debug_terminal_printer);
					double log_probability_false = log_joint_probability_of_truth(T_copy, proof_prior, proof_axioms_copy, negated, 100, 4, 20, T_MAP_false, proof_MAP_false);
					for (unsigned int j = 0; isinf(log_probability_false) && j < 400; j++) {
						null_collector collector;
						for (unsigned int t = 0; t < 10; t++)
							do_exploratory_mh_step(T_copy, proof_prior, proof_axioms_copy, collector);
						log_probability_false = log_joint_probability_of_truth(T_copy, proof_prior, proof_axioms_copy, negated, 100, 4, 20, T_MAP_false, proof_MAP_false);
					}
					free(*negated); if (negated->reference_count == 0) free(negated);

					if (fabs(log_probability_true - log_probability_false) < PREDICT_UNKNOWN_THRESHOLD) {
						if (job.label != ruletaker_label::UNKNOWN) {
							if (!isinf(log_probability_true)) print_theory(T_MAP_true, proof_MAP_true, proof_prior);
							if (!isinf(log_probability_false)) print_theory(T_MAP_false, proof_MAP_false, proof_prior);
						}
					} else if (log_probability_true > log_probability_false) {
						if (job.label != ruletaker_label::TRUE) {
							if (!isinf(log_probability_true)) print_theory(T_MAP_true, proof_MAP_true, proof_prior);
							if (!isinf(log_probability_false)) print_theory(T_MAP_false, proof_MAP_false, proof_prior);
						}
					} else if (log_probability_false > log_probability_true) {
						if (job.label != ruletaker_label::FALSE) {
							if (!isinf(log_probability_true)) print_theory(T_MAP_true, proof_MAP_true, proof_prior);
							if (!isinf(log_probability_false)) print_theory(T_MAP_false, proof_MAP_false, proof_prior);
						}
					}

					results_lock.lock();
					results.add({job.context_id, job.question_id, log_probability_true - log_probability_false, job.label});
					results_lock.unlock();
					if (!isinf(log_probability_true)) free(T_MAP_true);
					if (!isinf(log_probability_false)) free(T_MAP_false);
				} else {
					if (job.label != ruletaker_label::FALSE) {
						if (!isinf(log_probability_true)) print_theory(T_MAP_true, proof_MAP_true, proof_prior);
					}

					results_lock.lock();
					results.add({job.context_id, job.question_id, -std::numeric_limits<double>::infinity(), job.label});
					results_lock.unlock();
				}
				total++;
				free_logical_forms(logical_forms, parse_count);
			} else {
				total++;
			}
			free(T_copy);
			free(job);

		} else {
			num_threads_reading_context++;
			ruletaker_context_item<Theory, PriorStateType>& job = context_queue[context_queue_start++];
			lock.unlock();
if (job.context_id != 63 - 1) { // != 6 - 1) { //< 10 - 1 || job.context_id >= 139 - 1) {
total += job.questions.length;
num_threads_reading_context--;
free(job);
continue;
}

			/* for reproducibility, reset the PRNG state */
			core::engine = prng_engine;

			/* read the context sentences */
			unsigned int i = 0;
			unsigned int start = 0;
			for (; job.context[i] != '\0'; i++) {
				if (job.context[i] == '.') {
					const char old_next = job.context[i + 1];
					job.context[i + 1] = '\0';
					if (!read_sentence(corpus, parser, job.context + start, job.T, names, seed_entities, proof_prior, job.proof_axioms, 10, UINT_MAX)) {
						std::unique_lock<std::mutex> lock(results_lock);
						if (!unparseable_context.ensure_capacity(unparseable_context.length + 1)
						 || !init(unparseable_context[unparseable_context.length].value, job.context + start))
						{
							job.context[i + 1] = old_next;
							status = false;
							num_threads_running--;
							num_threads_reading_context--;
							work_queue_cv.notify_all();
							free(job);
							for (auto entry : names) free(entry.key);
							free(parser); return;
						}
						unparseable_context[unparseable_context.length++].key = job.context_id;
						job.context[i + 1] = old_next;
						break;
					}
					job.context[i + 1] = old_next;
					start = i + 1;
					while (isspace(job.context[start])) start++;

					Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
					PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
					hash_map<const hol_term*, hol_term*> formula_map(128);
					Theory::clone(job.T, T_MAP, formula_map);
					new (&proof_axioms_MAP) PriorStateType(job.proof_axioms, formula_map);
					auto collector = make_log_probability_collector(job.T, proof_prior);
					double max_log_probability = collector.current_log_probability;
					for (unsigned int j = 0; j < 4; j++) {
						for (unsigned int t = 0; t < 100; t++) {
							bool print_debug = false;
							if (print_debug) job.T.template print_axioms<true>(stderr, *debug_terminal_printer);
							if (print_debug) job.T.print_disjunction_introductions(stderr, *debug_terminal_printer);
							do_mh_step(job.T, proof_prior, job.proof_axioms, collector);
							if (collector.current_log_probability > max_log_probability) {
								free(T_MAP); proof_axioms_MAP.~PriorStateType(); formula_map.clear();
								Theory::clone(job.T, T_MAP, formula_map);
								new (&proof_axioms_MAP) PriorStateType(job.proof_axioms, formula_map);
								max_log_probability = collector.current_log_probability;
							}
						}

						if (j + 1 < 4) {
							for (unsigned int t = 0; t < 20; t++)
								do_exploratory_mh_step(job.T, proof_prior, job.proof_axioms, collector);
						}
					}
					free(job.T); job.proof_axioms.~PriorStateType(); formula_map.clear();
					Theory::clone(T_MAP, job.T, formula_map);
					new (&job.proof_axioms) PriorStateType(proof_axioms_MAP, formula_map);
					T_MAP.print_axioms(stderr, *debug_terminal_printer);
					free(T_MAP); proof_axioms_MAP.~PriorStateType();
				}
			}

			if (job.context[i] == '\0') {
				/* if we successfully read the context, enqueue the jobs for reading/answering the associated questions */
				job.prng_engine = core::engine;
				std::unique_lock<std::mutex> lock(work_queue_lock);
				if (question_queue_length + job.questions.length > MAX_QUESTION_COUNT) {
					fprintf(stderr, "do_ruletaker_experiments ERROR: Requested question queue length exceeds `MAX_QUESTION_COUNT`.\n");
					status = false;
					num_threads_running--;
					num_threads_reading_context--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				for (unsigned int j = 0; j < job.questions.length; j++) {
					ruletaker_question_item<Theory, PriorStateType>& new_question = question_queue[question_queue_length];
					set_empty(new_question.T);
					new_question.context_id = job.context_id;
					new_question.question_id = j;
					if (!init(new_question.question, job.questions[j].key.length + 1)) {
						status = false;
						num_threads_running--;
						num_threads_reading_context--;
						work_queue_cv.notify_all();
						free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
					for (unsigned int k = 0; k < job.questions[j].key.length; k++)
						new_question.question[k] = job.questions[j].key[k];
					new_question.question[job.questions[j].key.length] = '\0';
					new_question.question.length = job.questions[j].key.length;
					new_question.label = job.questions[j].value;
					hash_map<const hol_term*, hol_term*> formula_map(128);
					if (!Theory::clone(job.T, new_question.T, formula_map)) {
						status = false;
						num_threads_running--;
						num_threads_reading_context--;
						work_queue_cv.notify_all();
						set_empty(new_question.T);
						free(new_question); free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					} else if (new (&new_question.proof_axioms) PriorStateType(job.proof_axioms, formula_map) == nullptr) {
						status = false;
						num_threads_running--;
						num_threads_reading_context--;
						work_queue_cv.notify_all();
						free(new_question.T);
						set_empty(new_question.T);
						free(new_question); free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
					question_queue_length++;
					work_queue_cv.notify_one();
				}
			} else {
				total += job.questions.length;
			}
			num_threads_reading_context--;
			free(job);
		}
	}

	num_threads_running--;
	for (auto entry : names) free(entry.key);
	free(parser);
}

inline char label_to_char(ruletaker_label label) {
	switch (label) {
	case ruletaker_label::UNKNOWN: return '?';
	case ruletaker_label::TRUE: return 'T';
	case ruletaker_label::FALSE: return 'F';
	}
	fprintf(stderr, "label_to_char ERROR: Unrecognized `ruletaker_label`.\n");
	exit(EXIT_FAILURE);
}

inline void print_ruletaker_results(
		const std::atomic_uint& total,
		array<question_result>& results,
		array<pair<unsigned int, string>>& unparseable_context,
		std::mutex& results_lock,
		const char* output_filepath)
{
	std::unique_lock<std::mutex> lock(results_lock);
	insertion_sort(results);
	FILE* out = open_file(output_filepath, "w");
	if (out == nullptr) {
		fprintf(stderr, "ERROR: Unable to open `%s` for writing.\n", output_filepath);
		return;
	}
	array<pair<unsigned int, unsigned int>> incorrect(64);
	for (const question_result& result : results) {
		if (result.question_id + 1 < 10)
			fprintf(out, "[%u, %u]  ", result.context_id + 1, result.question_id + 1);
		else fprintf(out, "[%u, %u] ", result.context_id + 1, result.question_id + 1);
		fprintf(out, "%c %lf", label_to_char(result.true_label), result.log_probability_diff);
		bool correct = true;
		if (fabs(result.log_probability_diff) < PREDICT_UNKNOWN_THRESHOLD) {
			if (result.true_label != ruletaker_label::UNKNOWN)
				correct = false;
		} else if (result.log_probability_diff > 0.0) {
			if (result.true_label != ruletaker_label::TRUE)
				correct = false;
		} else {
			if (result.true_label != ruletaker_label::FALSE)
				correct = false;
		}

		if (correct) {
			fputc('\n', out);
		} else {
			incorrect.add(make_pair(result.context_id, result.question_id));
			fprintf(out, " *\n");
		}
	}
	fprintf(out,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %lu\n"
			"  Incorrect questions: %lu\n",
			total.load(), results.length,
			incorrect.length);
	if (unparseable_context.length != 0) {
		fprintf(out, "Failed to parse following context sentences:\n");
		insertion_sort(unparseable_context, pair_sorter());
		for (const auto& entry : unparseable_context) {
			fprintf(out, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, out); print("\"\n", out);
		}
	}
	fclose(out);

	fprintf(stderr,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %lu\n"
			"  Incorrect questions: %lu\n",
			total.load(), results.length,
			incorrect.length);
	if (incorrect.length != 0) {
		fprintf(stderr, "Incorrect questions:\n");
		insertion_sort(incorrect, pair_sorter());
		for (const auto& entry : incorrect)
			fprintf(stderr, "  Context ID %u, Question ID %u\n", entry.key + 1, entry.value + 1);
	} if (unparseable_context.length != 0) {
		fprintf(stderr, "Failed to parse following context sentences:\n");
		insertion_sort(unparseable_context, pair_sorter());
		for (const auto& entry : unparseable_context) {
			fprintf(stderr, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, stderr); print("\"\n", stderr);
		}
	}
}

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
bool run_ruletaker_experiments(
		ArticleSource& corpus, Parser& parser,
		Theory& T, PriorStateType& proof_axioms,
		ProofPrior& proof_prior,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& seed_entities,
		const char* data_filepath,
		const char* results_filepath,
		unsigned int thread_count)
{
	bool status = true;
	ruletaker_context_item<Theory, PriorStateType>* context_queue = (ruletaker_context_item<Theory, PriorStateType>*) malloc(sizeof(ruletaker_context_item<Theory, PriorStateType>) * MAX_CONTEXT_COUNT);
	if (context_queue == nullptr) return false;
	ruletaker_question_item<Theory, PriorStateType>* question_queue = (ruletaker_question_item<Theory, PriorStateType>*) malloc(sizeof(ruletaker_question_item<Theory, PriorStateType>) * MAX_QUESTION_COUNT);
	if (question_queue == nullptr) {
		free(context_queue);
		return false;
	}
	unsigned int context_queue_start = 0;
	unsigned int question_queue_start = 0;
	unsigned int context_queue_length = 0;
	unsigned int question_queue_length = 0;
	std::mutex work_queue_lock;
	std::condition_variable work_queue_cv;
	std::minstd_rand prng_engine = core::engine;
	std::mutex results_lock;
	array<question_result> results(64);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(1);
	std::atomic_uint num_threads_running(0);

	std::thread* workers = new std::thread[thread_count];
	for (unsigned int i = 0; i < thread_count; i++) {
		workers[i] = std::thread(
				do_ruletaker_experiments<ArticleSource, Parser, Theory, PriorStateType, ProofPrior>,
				std::ref(status), context_queue, question_queue,
				std::ref(context_queue_start), std::ref(question_queue_start),
				std::ref(context_queue_length), std::ref(question_queue_length),
				std::ref(work_queue_lock), std::ref(work_queue_cv),
				prng_engine, std::ref(corpus),
				std::ref(parser), std::ref(proof_prior),
				std::ref(names), std::ref(seed_entities),
				std::ref(results_lock), std::ref(results),
				std::ref(unparseable_context), std::ref(total),
				std::ref(num_threads_reading_context),
				std::ref(num_threads_running));
	}

	unsigned int context_id = 0;
	auto process_ruletaker_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms](char* context, array<pair<string, ruletaker_label>>& questions)
	{
		if (context_queue_length + 1 > MAX_CONTEXT_COUNT) {
			fprintf(stderr, "run_ruletaker_experiments ERROR: Requested context queue length exceeds `MAX_CONTEXT_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		ruletaker_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_ruletaker_experiments ERROR: Out of memory.\n");
			return false;
		} else if (!array_init(new_context.questions, questions.length)) {
			free(new_context.context);
			return false;
		}
		unsigned int i;
		for (i = 0; context[i] != '\0'; i++)
			new_context.context[i] = context[i];
		new_context.context[i] = '\0';
		for (const auto& entry : questions) {
			if (!init(new_context.questions[new_context.questions.length].key, entry.key)) {
				free(new_context);
				return false;
			}
			new_context.questions[new_context.questions.length].value = entry.value;
			new_context.questions.length++;
		}
		hash_map<const hol_term*, hol_term*> formula_map(128);
		if (!Theory::clone(T, new_context.T, formula_map)) {
			set_empty(new_context.T);
			free(new_context);
			return false;
		} else if (new (&new_context.proof_axioms) PriorStateType(proof_axioms, formula_map) == nullptr) {
			free(new_context.T); set_empty(new_context.T);
			free(new_context); return false;
		}
		context_queue_length++;
		work_queue_cv.notify_one();
		return true;
	};

	if (!read_ruletaker_data(data_filepath, process_ruletaker_questions))
		status = false;
	num_threads_reading_context--;

	timer stopwatch;
	while (status) {
		if (num_threads_reading_context == 0 && context_queue_start == context_queue_length && question_queue_start == question_queue_length)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_ruletaker_results(total, results, unparseable_context, results_lock, results_filepath);
			stopwatch.start();
		}
	}

	work_queue_cv.notify_all();
	while (status) {
		if (num_threads_running == 0)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_ruletaker_results(total, results, unparseable_context, results_lock, results_filepath);
			stopwatch.start();
		}
	}
	for (unsigned int i = 0; i < thread_count; i++) {
		if (!workers[i].joinable()) continue;
		try {
			workers[i].join();
		} catch (...) { }
	}
	print_ruletaker_results(total, results, unparseable_context, results_lock, results_filepath);
	delete[] workers;
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
bool run_ruletaker_experiments_single_threaded(
		ArticleSource& corpus, Parser& parser,
		Theory& T, PriorStateType& proof_axioms,
		ProofPrior& proof_prior,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& seed_entities,
		const char* data_filepath,
		const char* results_filepath)
{
	bool status = true;
	ruletaker_context_item<Theory, PriorStateType>* context_queue = (ruletaker_context_item<Theory, PriorStateType>*) malloc(sizeof(ruletaker_context_item<Theory, PriorStateType>) * MAX_CONTEXT_COUNT);
	if (context_queue == nullptr) return false;
	ruletaker_question_item<Theory, PriorStateType>* question_queue = (ruletaker_question_item<Theory, PriorStateType>*) malloc(sizeof(ruletaker_question_item<Theory, PriorStateType>) * MAX_QUESTION_COUNT);
	if (question_queue == nullptr) {
		free(context_queue);
		return false;
	}
	unsigned int context_queue_start = 0;
	unsigned int question_queue_start = 0;
	unsigned int context_queue_length = 0;
	unsigned int question_queue_length = 0;
	std::mutex work_queue_lock;
	std::condition_variable work_queue_cv;
	std::minstd_rand prng_engine = core::engine;
	std::mutex results_lock;
	array<question_result> results(64);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(0);
	std::atomic_uint num_threads_running(0);

	unsigned int context_id = 0;
	auto process_ruletaker_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms](char* context, array<pair<string, ruletaker_label>>& questions)
	{
		if (context_queue_length + 1 > MAX_CONTEXT_COUNT) {
			fprintf(stderr, "run_ruletaker_experiments_single_threaded ERROR: Requested context queue length exceeds `MAX_CONTEXT_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		ruletaker_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_ruletaker_experiments_single_threaded ERROR: Out of memory.\n");
			return false;
		} else if (!array_init(new_context.questions, questions.length)) {
			free(new_context.context);
			return false;
		}
		unsigned int i;
		for (i = 0; context[i] != '\0'; i++)
			new_context.context[i] = context[i];
		new_context.context[i] = '\0';
		for (const auto& entry : questions) {
			if (!init(new_context.questions[new_context.questions.length].key, entry.key)) {
				free(new_context);
				return false;
			}
			new_context.questions[new_context.questions.length].value = entry.value;
			new_context.questions.length++;
		}
		hash_map<const hol_term*, hol_term*> formula_map(128);
		if (!Theory::clone(T, new_context.T, formula_map)) {
			set_empty(new_context.T);
			free(new_context);
			return false;
		} else if (new (&new_context.proof_axioms) PriorStateType(proof_axioms, formula_map) == nullptr) {
			free(new_context.T); set_empty(new_context.T);
			free(new_context); return false;
		}
		context_queue_length++;
		work_queue_cv.notify_one();
		return true;
	};

	if (!read_ruletaker_data(data_filepath, process_ruletaker_questions))
		status = false;

	do_ruletaker_experiments(status, context_queue, question_queue,
			context_queue_start, question_queue_start, context_queue_length,
			question_queue_length, work_queue_lock, work_queue_cv, prng_engine,
			corpus, parser, proof_prior, names, seed_entities,
			results_lock, results, unparseable_context, total,
			num_threads_reading_context, num_threads_running);

	print_ruletaker_results(total, results, unparseable_context, results_lock, results_filepath);
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

#endif /* RULETAKER_H_ */
