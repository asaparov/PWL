/**
 * geoquery.h - Code to read GeoQuery data and attempt to answer the questions.
 *
 *  Created on: Mar 22, 2021
 *      Author: asaparov
 */

#ifndef GEOQUERY_H_
#define GEOQUERY_H_

#include <core/utility.h>
#include <stdio.h>
#include <atomic>

using namespace core;

#include <atomic>

enum class geoquery_work_item_type {
	READ_CONTEXT,
	ANSWER_QUESTION
};

template<typename Theory, typename PriorStateType>
struct geoquery_context_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	std::minstd_rand prng_engine;
	char* context;
	array<pair<string, string>> questions;

	static inline void free(geoquery_context_item<Theory, PriorStateType>& item) {
		core::free(item.context);
		for (auto& entry : item.questions) {
			core::free(entry.key);
			core::free(entry.value);
		}
		core::free(item.questions);
		if (!core::is_empty(item.T)) {
			core::free(item.T);
			item.proof_axioms.~PriorStateType();
		}
	}
};

template<typename Theory, typename PriorStateType>
struct geoquery_question_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	unsigned int question_id;
	string question;
	string label;

	static inline void free(geoquery_question_item& item) {
		core::free(item.question);
		core::free(item.label);
		if (!core::is_empty(item.T)) {
			core::free(item.T);
			item.proof_axioms.~PriorStateType();
		}
	}
};

struct geoquery_question_result {
	unsigned int context_id;
	unsigned int question_id;
	string answer;
	string label;

	static inline void swap(geoquery_question_result& first, geoquery_question_result& second) {
		core::swap(first.context_id, second.context_id);
		core::swap(first.question_id, second.question_id);
		core::swap(first.answer, second.answer);
		core::swap(first.label, second.label);
	}
};

inline bool operator < (const geoquery_question_result& first, const geoquery_question_result& second) {
	if (first.context_id < second.context_id) return true;
	else if (first.context_id > second.context_id) return false;
	else return first.question_id < second.question_id;
}

inline bool operator == (const geoquery_question_result& first, const geoquery_question_result& second) {
	return first.context_id == second.context_id
		&& first.question_id == second.question_id;
}

constexpr unsigned int MAX_GEOQUERY_QUESTION_COUNT = 280;

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
void do_geoquery_experiments(bool& status,
		geoquery_context_item<Theory, PriorStateType>* context_queue,
		geoquery_question_item<Theory, PriorStateType>* question_queue,
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
		array<geoquery_question_result>& results,
		array<pair<unsigned int, string>>& unparseable_questions,
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
			geoquery_question_item<Theory, PriorStateType>& job = question_queue[question_queue_start++];
			lock.unlock();

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
				/* TODO: add question answering logic here */

				results_lock.lock();
				results.ensure_capacity(results.length + 1);
				results[results.length].context_id = job.context_id;
				results[results.length].question_id = job.question_id;
				if (!init(results[results.length].answer, "")) {
					free_logical_forms(logical_forms, parse_count);
					status = false;
					num_threads_running--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				} else if (!init(results[results.length].label, job.label)) {
					free(results[results.length].answer);
					free_logical_forms(logical_forms, parse_count);
					status = false;
					num_threads_running--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				results.length++;
				results_lock.unlock();

				total++;
				free_logical_forms(logical_forms, parse_count);
			} else {
				results_lock.lock();
				if (!unparseable_questions.ensure_capacity(unparseable_questions.length + 1)
				 || !init(unparseable_questions[unparseable_questions.length].value, job.question))
				{
					status = false;
					num_threads_running--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				unparseable_questions[unparseable_questions.length++].key = job.context_id;
				results_lock.unlock();

				total++;
			}
			free(T_copy);
			free(job);

		} else {
			num_threads_reading_context++;
			geoquery_context_item<Theory, PriorStateType>& job = context_queue[context_queue_start++];
			lock.unlock();

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
							if (print_debug) job.T.template print_axioms<true>(stdout, *debug_terminal_printer);
							if (print_debug) { job.T.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout); }
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
					T_MAP.print_axioms(stdout, *debug_terminal_printer); fflush(stdout);
					free(T_MAP); proof_axioms_MAP.~PriorStateType();
				}
			}

			if (job.context[i] == '\0') {
				/* if we successfully read the context, enqueue the jobs for reading/answering the associated questions */
				job.prng_engine = core::engine;
				std::unique_lock<std::mutex> lock(work_queue_lock);
				if (question_queue_length + job.questions.length > MAX_GEOQUERY_QUESTION_COUNT) {
					fprintf(stderr, "do_geoquery_experiments ERROR: Requested question queue length exceeds `MAX_GEOQUERY_QUESTION_COUNT`.\n");
					status = false;
					num_threads_running--;
					num_threads_reading_context--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				for (unsigned int j = 0; j < job.questions.length; j++) {
					geoquery_question_item<Theory, PriorStateType>& new_question = question_queue[question_queue_length];
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
					if (!init(new_question.label, job.questions[j].value)) {
						status = false;
						num_threads_running--;
						num_threads_reading_context--;
						work_queue_cv.notify_all();
						free(new_question.question); free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
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

inline void print_geoquery_results(
		const std::atomic_uint& total,
		array<geoquery_question_result>& results,
		array<pair<unsigned int, string>>& unparseable_questions,
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
	for (const geoquery_question_result& result : results) {
		fprintf(out, "[%u]\n", result.context_id + 1);
		fprintf(out, "Predicted: "); print(result.answer, out); print('\n', out);
		fprintf(out, "    Label: "); print(result.label, out); print("\n\n", out);
	}
	fprintf(out,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %lu\n",
			total.load(), results.length);
	if (unparseable_context.length != 0) {
		fprintf(out, "Failed to parse following context sentences:\n");
		insertion_sort(unparseable_context, pair_sorter());
		for (const auto& entry : unparseable_context) {
			fprintf(out, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, out); print("\"\n", out);
		}
	} if (unparseable_questions.length != 0) {
		fprintf(out, "Failed to parse query sentences:\n");
		insertion_sort(unparseable_questions, pair_sorter());
		for (const auto& entry : unparseable_questions) {
			fprintf(out, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, out); print("\"\n", out);
		}
	}
	fclose(out);

	fprintf(stdout,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %lu\n",
			total.load(), results.length);
	if (unparseable_context.length != 0) {
		fprintf(stdout, "Failed to parse following context sentences:\n");
		insertion_sort(unparseable_context, pair_sorter());
		for (const auto& entry : unparseable_context) {
			fprintf(stdout, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, stdout); print("\"\n", stdout);
		}
	} if (unparseable_questions.length != 0) {
		fprintf(stdout, "Failed to parse query sentences:\n");
		insertion_sort(unparseable_questions, pair_sorter());
		for (const auto& entry : unparseable_questions) {
			fprintf(stdout, "  Context ID %u: \"", entry.key + 1);
			print(entry.value, stdout); print("\"\n", stdout);
		}
	}
	fflush(stdout);
}

template<typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
bool run_geoquery_experiments(
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
	geoquery_context_item<Theory, PriorStateType>* context_queue = (geoquery_context_item<Theory, PriorStateType>*) malloc(sizeof(geoquery_context_item<Theory, PriorStateType>) * MAX_GEOQUERY_QUESTION_COUNT);
	if (context_queue == nullptr) return false;
	geoquery_question_item<Theory, PriorStateType>* question_queue = (geoquery_question_item<Theory, PriorStateType>*) malloc(sizeof(geoquery_question_item<Theory, PriorStateType>) * MAX_GEOQUERY_QUESTION_COUNT);
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
	array<geoquery_question_result> results(64);
	array<pair<unsigned int, string>> unparseable_questions(4);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(1);
	std::atomic_uint num_threads_running(0);

	std::thread* workers = new std::thread[thread_count];
	for (unsigned int i = 0; i < thread_count; i++) {
		workers[i] = std::thread(
				do_geoquery_experiments<ArticleSource, Parser, Theory, PriorStateType, ProofPrior>,
				std::ref(status), context_queue, question_queue,
				std::ref(context_queue_start), std::ref(question_queue_start),
				std::ref(context_queue_length), std::ref(question_queue_length),
				std::ref(work_queue_lock), std::ref(work_queue_cv),
				prng_engine, std::ref(corpus),
				std::ref(parser), std::ref(proof_prior),
				std::ref(names), std::ref(seed_entities),
				std::ref(results_lock), std::ref(results),
				std::ref(unparseable_questions),
				std::ref(unparseable_context), std::ref(total),
				std::ref(num_threads_reading_context),
				std::ref(num_threads_running));
	}

	unsigned int context_id = 0;
	auto process_geoquery_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms](char* context, array<pair<string, string>>& questions)
	{
		if (context_queue_length + 1 > MAX_GEOQUERY_QUESTION_COUNT) {
			fprintf(stderr, "run_geoquery_experiments ERROR: Requested context queue length exceeds `MAX_GEOQUERY_QUESTION_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		geoquery_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_geoquery_experiments ERROR: Out of memory.\n");
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
			} else if (!init(new_context.questions[new_context.questions.length].value, entry.value)) {
				free(new_context.questions[new_context.questions.length].key);
				free(new_context); return false;
			}
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

	if (!read_ruletaker_data<string>(data_filepath, process_geoquery_questions))
		status = false;
	num_threads_reading_context--;

	timer stopwatch;
	while (status) {
		if (num_threads_reading_context == 0 && context_queue_start == context_queue_length && question_queue_start == question_queue_length)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_geoquery_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath);
			stopwatch.start();
		}
	}

	work_queue_cv.notify_all();
	while (status) {
		if (num_threads_running == 0)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_geoquery_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath);
			stopwatch.start();
		}
	}
	for (unsigned int i = 0; i < thread_count; i++) {
		if (!workers[i].joinable()) continue;
		try {
			workers[i].join();
		} catch (...) { }
	}
	print_geoquery_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath);
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
bool run_geoquery_experiments_single_threaded(
		ArticleSource& corpus, Parser& parser,
		Theory& T, PriorStateType& proof_axioms,
		ProofPrior& proof_prior,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& seed_entities,
		const char* data_filepath,
		const char* results_filepath)
{
	bool status = true;
	geoquery_context_item<Theory, PriorStateType>* context_queue = (geoquery_context_item<Theory, PriorStateType>*) malloc(sizeof(geoquery_context_item<Theory, PriorStateType>) * MAX_GEOQUERY_QUESTION_COUNT);
	if (context_queue == nullptr) return false;
	geoquery_question_item<Theory, PriorStateType>* question_queue = (geoquery_question_item<Theory, PriorStateType>*) malloc(sizeof(geoquery_question_item<Theory, PriorStateType>) * MAX_GEOQUERY_QUESTION_COUNT);
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
	array<geoquery_question_result> results(64);
	array<pair<unsigned int, string>> unparseable_questions(4);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(0);
	std::atomic_uint num_threads_running(0);

	unsigned int context_id = 0;
	auto process_geoquery_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms](char* context, array<pair<string, string>>& questions)
	{
		if (context_queue_length + 1 > MAX_GEOQUERY_QUESTION_COUNT) {
			fprintf(stderr, "run_geoquery_experiments_single_threaded ERROR: Requested context queue length exceeds `MAX_GEOQUERY_QUESTION_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		geoquery_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_geoquery_experiments_single_threaded ERROR: Out of memory.\n");
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
			} else if (!init(new_context.questions[new_context.questions.length].value, entry.value)) {
				free(new_context.questions[new_context.questions.length].key);
				free(new_context); return false;
			}
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

	if (!read_ruletaker_data<string>(data_filepath, process_geoquery_questions))
		status = false;

	do_geoquery_experiments(status, context_queue, question_queue,
			context_queue_start, question_queue_start, context_queue_length,
			question_queue_length, work_queue_lock, work_queue_cv, prng_engine,
			corpus, parser, proof_prior, names, seed_entities,
			results_lock, results, unparseable_questions, unparseable_context,
			total, num_threads_reading_context, num_threads_running);

	print_geoquery_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath);
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

#endif /* GEOQUERY_H_ */