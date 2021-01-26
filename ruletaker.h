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
	CONTEXT,
	QUESTIONS,
	TEXT,
	LABEL,
	OTHER
};

enum class ruletaker_reader_state {
	START,
	ROOT_OBJECT,
	QUESTIONS,
	QUESTION
};

struct ruletaker_reader {
	ruletaker_reader_state state;
	ruletaker_key current_key;
	char* context;
	pair<char*, bool> next_question;
	array<pair<string, bool>> questions;

	ruletaker_reader() : state(ruletaker_reader_state::START), context(nullptr), next_question(nullptr, false), questions(44) { }

	~ruletaker_reader() {
		if (context != nullptr)
			free(context);
		if (next_question.key != nullptr)
			free(next_question.key);
		for (pair<string, bool>& question : questions)
			free(question.key);
	}
};

inline bool begin_object(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::START) {
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	} else if (reader.state == ruletaker_reader_state::QUESTIONS) {
		reader.state = ruletaker_reader_state::QUESTION;
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
	}
	return true;
}

inline bool emit_key(const array<char>& key_name, const position& pos, ruletaker_reader& reader) {
	if (compare_strings(key_name, "context")) {
		reader.current_key = ruletaker_key::CONTEXT;
	} else if (compare_strings(key_name, "questions")) {
		reader.current_key = ruletaker_key::QUESTIONS;
	} else if (compare_strings(key_name, "text")) {
		reader.current_key = ruletaker_key::TEXT;
	} else if (compare_strings(key_name, "label")) {
		reader.current_key = ruletaker_key::LABEL;
	} else {
		reader.current_key = ruletaker_key::OTHER;
	}
	return true;
}

inline bool begin_list(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::QUESTIONS)
		reader.state = ruletaker_reader_state::QUESTIONS;
	return true;
}

inline bool end_list(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::QUESTIONS)
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	return true;
}

inline bool emit_string(const array<char>& str, const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::CONTEXT) {
		if (reader.context != nullptr) {
			read_error("Found duplicate `context` entry", pos);
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
	} else if (reader.current_key == ruletaker_key::TEXT) {
		if (reader.next_question.key != nullptr) {
			read_error("Found duplicate `text` entry for question", pos);
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
	}
	return true;
}

inline bool emit_true(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::LABEL)
		reader.next_question.value = true;
	return true;
}

inline bool emit_false(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::LABEL)
		reader.next_question.value = false;
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
	array<pair<string, bool>> questions;

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
	bool label;

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
		if (head->type == hol_term_type::EXISTS) {
			hol_term* operand = head->quantifier.operand;
			bool is_same = false;
			if (operand->type == hol_term_type::AND) {
				for (unsigned int i = 0; !is_same && i < operand->array.length; i++) {
					hol_term* conjunct = operand->array.operands[i];
					if (conjunct->type == hol_term_type::UNARY_APPLICATION && *conjunct->binary.left == hol_term::constants<(unsigned int) built_in_predicates::SAME>::value
					 && conjunct->binary.right->type == hol_term_type::VARIABLE && conjunct->binary.right->variable == head->quantifier.variable)
					{
						is_same = true;
					}
				}
			} else {
				if (operand->type == hol_term_type::UNARY_APPLICATION && *operand->binary.left == hol_term::constants<(unsigned int) built_in_predicates::SAME>::value
				 && operand->binary.right->type == hol_term_type::VARIABLE && operand->binary.right->variable == head->quantifier.variable)
				{
					is_same = true;
				}
			}

			if (is_same && scopes.length > 1) {
				/* negate the scope immediately preceding `head` */
				head = scopes[scopes.length - 2];
			}
		}

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

constexpr unsigned int MAX_CONTEXT_COUNT = 140;
constexpr unsigned int MAX_QUESTION_COUNT = 5270;

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
void do_ruletaker_experiments(
		bool& status, bool& running,
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
		std::mutex& incorrect_question_ids_lock,
		array<pair<unsigned int, unsigned int>>& incorrect,
		array<pair<unsigned int, unsigned int>>& half_correct,
		array<pair<unsigned int, string>>& unparseable_context,
		std::atomic_uint& total,
		std::atomic_uint& answered,
		std::atomic_uint& num_threads_reading_context)
{
	Parser& parser = *((Parser*) alloca(sizeof(Parser)));
	if (!init(parser, parser_src)) {
		status = false;
		work_queue_cv.notify_all();
		return;
	}
	hash_map<string, unsigned int> names(names_src.table.capacity);
	for (const auto& entry : names_src) {
		unsigned int index = names.table.index_to_insert(entry.key);
		if (!init(names.table.keys[index], entry.key)) {
			status = false;
			work_queue_cv.notify_all();
			for (auto entry : names) free(entry.key);
			free(parser); return;
		}
		names.values[index] = entry.value;
		names.table.size++;
	}
	if (!parser.invert_name_map(names)) {
		status = false;
		work_queue_cv.notify_all();
		for (auto entry : names) free(entry.key);
		free(parser); return;
	}

	while (status && running)
	{
		std::unique_lock<std::mutex> lock(work_queue_lock);
		while (status && running && context_queue_length == context_queue_start && question_queue_length == question_queue_start)
			work_queue_cv.wait(lock);
		if (!status || !running) {
			for (auto entry : names) free(entry.key);
			free(parser); return;
		}

		if (question_queue_start < question_queue_length) {
			ruletaker_question_item<Theory, PriorStateType>& job = question_queue[question_queue_start++];
			lock.unlock();
if (job.question_id < 18 - 1)
{
total++;
continue;
}

			/* for reproducibility, reset the PRNG state */
			core::engine = context_queue[job.context_id].prng_engine;

			/* first clone the theory from the appropriate work item */
			Theory& T_copy = *((Theory*) alloca(sizeof(Theory)));
			hash_map<const hol_term*, hol_term*> formula_map(128);
			if (!Theory::clone(job.T, T_copy, formula_map)) {
				status = false;
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
				double log_probability_true = log_joint_probability_of_truth(job.T, proof_prior, job.proof_axioms, logical_forms[0], 10000, T_MAP_true, proof_MAP_true);
				for (unsigned int j = 0; isinf(log_probability_true) && j < 1000; j++) {
					null_collector collector;
					for (unsigned int t = 0; t < 10; t++)
						do_exploratory_mh_step(job.T, proof_prior, job.proof_axioms, collector);
					log_probability_true = log_joint_probability_of_truth(job.T, proof_prior, job.proof_axioms, logical_forms[0], 10000, T_MAP_true, proof_MAP_true);
				}

				hol_term* negated;
				if (!negate_head(logical_forms[0], negated) || negated == nullptr) {
					free_logical_forms(logical_forms, parse_count);
					status = false;
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
				double log_probability_false = log_joint_probability_of_truth(T_copy, proof_prior, proof_axioms_copy, negated, 10000, T_MAP_false, proof_MAP_false);
				for (unsigned int j = 0; isinf(log_probability_false) && j < 1000; j++) {
					null_collector collector;
					for (unsigned int t = 0; t < 10; t++)
						do_exploratory_mh_step(T_copy, proof_prior, proof_axioms_copy, collector);
					log_probability_false = log_joint_probability_of_truth(T_copy, proof_prior, proof_axioms_copy, negated, 10000, T_MAP_false, proof_MAP_false);
				}
				free(*negated); if (negated->reference_count == 0) free(negated);
				if (log_probability_true > log_probability_false) {
					if (!job.label) {
						print_theory(T_MAP_true, proof_MAP_true, proof_prior);
						print_theory(T_MAP_false, proof_MAP_false, proof_prior);
 						std::unique_lock<std::mutex> lock(incorrect_question_ids_lock);
						incorrect.add(make_pair(job.context_id, job.question_id));
					}
				} else if (log_probability_false > log_probability_true) {
					if (job.label) {
						print_theory(T_MAP_true, proof_MAP_true, proof_prior);
						print_theory(T_MAP_false, proof_MAP_false, proof_prior);
						std::unique_lock<std::mutex> lock(incorrect_question_ids_lock);
						incorrect.add(make_pair(job.context_id, job.question_id));
					}
				} else {
					std::unique_lock<std::mutex> lock(incorrect_question_ids_lock);
					half_correct.add(make_pair(job.context_id, job.question_id));
				}
				free(T_MAP_true);
				free(T_MAP_false);
				answered++;
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
if (job.context_id != 139 - 1) { // != 6 - 1) { //< 10 - 1 || job.context_id >= 139 - 1) {
total += job.questions.length;
num_threads_reading_context--;
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
						std::unique_lock<std::mutex> lock(incorrect_question_ids_lock);
						if (!unparseable_context.ensure_capacity(unparseable_context.length + 1)
						 || !init(unparseable_context[unparseable_context.length].value, job.context + start))
						{
							job.context[i + 1] = old_next;
							status = false;
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
					hash_map<const hol_term*, hol_term*> formula_map(128);
					Theory::clone(job.T, T_MAP, formula_map);
					auto collector = make_log_probability_collector(job.T, proof_prior);
					double max_log_probability = collector.current_log_probability;
					for (unsigned int t = 0; t < 4000; t++) {
if (i >= 329) {
fprintf(stderr, "DEBUG: t = %u\n", t);
job.proof_axioms.check_proof_axioms(job.T);
job.proof_axioms.check_universal_eliminations(job.T, collector);
job.T.check_concept_axioms();
job.T.check_disjunction_introductions();
job.T.are_elements_provable();
job.T.sets.check_freeable_sets();
job.T.sets.are_descendants_valid();
job.T.sets.are_set_sizes_valid();
job.T.sets.check_set_ids();
job.T.print_axioms(stderr, *debug_terminal_printer);
job.T.print_disjunction_introductions(stderr, *debug_terminal_printer);
}
						bool print_debug = false;
						if (print_debug) job.T.print_axioms(stderr, *debug_terminal_printer);
						if (print_debug) job.T.print_disjunction_introductions(stderr, *debug_terminal_printer);
						do_mh_step(job.T, proof_prior, job.proof_axioms, collector);
						if (collector.current_log_probability > max_log_probability) {
							free(T_MAP); formula_map.clear();
							Theory::clone(job.T, T_MAP, formula_map);
							max_log_probability = collector.current_log_probability;
						}
					}
					T_MAP.print_axioms(stderr, *debug_terminal_printer);
					free(T_MAP);
				}
			}

			if (job.context[i] == '\0') {
				/* if we successfully read the context, enqueue the jobs for reading/answering the associated questions */
				job.prng_engine = core::engine;
				std::unique_lock<std::mutex> lock(work_queue_lock);
				if (question_queue_length + job.questions.length > MAX_QUESTION_COUNT) {
					fprintf(stderr, "do_ruletaker_experiments ERROR: Requested question queue length exceeds `MAX_QUESTION_COUNT`.\n");
					status = false;
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
						num_threads_reading_context--;
						work_queue_cv.notify_all();
						set_empty(new_question.T);
						free(new_question); free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					} else if (new (&new_question.proof_axioms) PriorStateType(job.proof_axioms, formula_map) == nullptr) {
						status = false;
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

	for (auto entry : names) free(entry.key);
	free(parser);
}

inline void print_ruletaker_results(
		const std::atomic_uint& total, const std::atomic_uint& answered,
		array<pair<unsigned int, unsigned int>>& incorrect,
		array<pair<unsigned int, unsigned int>>& half_correct,
		array<pair<unsigned int, string>>& unparseable_context,
		std::mutex& incorrect_question_ids_lock)
{
	fprintf(stderr,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %u\n"
			"  Incorrect questions: %lu\n"
			"  Half-correct questions: %lu\n",
			total.load(), answered.load(),
			incorrect.length, half_correct.length);
	std::unique_lock<std::mutex> lock(incorrect_question_ids_lock);
	if (incorrect.length != 0) {
		fprintf(stderr, "Incorrect questions:\n");
		insertion_sort(incorrect, pair_sorter());
		for (const auto& entry : incorrect)
			fprintf(stderr, "  Context ID %u, Question ID %u\n", entry.key + 1, entry.value + 1);
	} if (half_correct.length != 0) {
		fprintf(stderr, "Half-correct questions:\n");
		insertion_sort(half_correct, pair_sorter());
		for (const auto& entry : half_correct)
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
		unsigned int thread_count)
{
	bool status = true; bool running = true;
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
	std::mutex incorrect_question_ids_lock;
	array<pair<unsigned int, unsigned int>> incorrect(64);
	array<pair<unsigned int, unsigned int>> half_correct(4);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint answered(0);
	std::atomic_uint num_threads_reading_context(0);

	std::thread* workers = new std::thread[thread_count];
	for (unsigned int i = 0; i < thread_count; i++) {
		workers[i] = std::thread(
				do_ruletaker_experiments<ArticleSource, Parser, Theory, PriorStateType, ProofPrior>,
				std::ref(status), std::ref(running), context_queue, question_queue,
				std::ref(context_queue_start), std::ref(question_queue_start),
				std::ref(context_queue_length), std::ref(question_queue_length),
				std::ref(work_queue_lock), std::ref(work_queue_cv),
				prng_engine, std::ref(corpus),
				std::ref(parser), std::ref(proof_prior),
				std::ref(names), std::ref(seed_entities),
				std::ref(incorrect_question_ids_lock),
				std::ref(incorrect), std::ref(half_correct),
				std::ref(unparseable_context),
				std::ref(total), std::ref(answered),
				std::ref(num_threads_reading_context));
	}

	unsigned int context_id = 0;
	auto process_ruletaker_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms](char* context, array<pair<string, bool>>& questions)
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

	timer stopwatch;
	while (status) {
		if (num_threads_reading_context == 0 && context_queue_start == context_queue_length && question_queue_start == question_queue_length)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_ruletaker_results(total, answered, incorrect, half_correct, unparseable_context, incorrect_question_ids_lock);
			stopwatch.start();
		}
	}

	running = false;
	work_queue_cv.notify_all();
	for (unsigned int i = 0; i < thread_count; i++) {
		if (!workers[i].joinable()) continue;
		try {
			workers[i].join();
		} catch (...) { }
	}
	print_ruletaker_results(total, answered, incorrect, half_correct, unparseable_context, incorrect_question_ids_lock);
	delete[] workers;
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

#endif /* RULETAKER_H_ */
