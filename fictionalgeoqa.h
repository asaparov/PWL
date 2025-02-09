/**
 * fictionalgeoqa.h - Code to read FictionalGeoQA data and attempt to answer the questions.
 *
 *  Created on: Mar 22, 2021
 *      Author: asaparov
 */

#ifndef FICTIONALGEOQA_H_
#define FICTIONALGEOQA_H_

#include <core/utility.h>
#include <stdio.h>
#include <atomic>

using namespace core;

#include <atomic>

enum class fictionalgeo_work_item_type {
	READ_CONTEXT,
	ANSWER_QUESTION
};

template<typename Theory, typename PriorStateType>
struct fictionalgeo_context_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	std::minstd_rand prng_engine;
	char* context;
	array<pair<string, string>> questions;

	static inline void free(fictionalgeo_context_item<Theory, PriorStateType>& item) {
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
struct fictionalgeo_question_item
{
	Theory T;
	PriorStateType proof_axioms;
	unsigned int context_id;
	unsigned int question_id;
	string question;
	string label;

	static inline void free(fictionalgeo_question_item& item) {
		core::free(item.question);
		core::free(item.label);
		if (!core::is_empty(item.T)) {
			core::free(item.T);
			item.proof_axioms.~PriorStateType();
		}
	}
};

struct fictionalgeo_question_result {
	unsigned int context_id;
	unsigned int question_id;
	string answer;
	string label;

	static inline void swap(fictionalgeo_question_result& first, fictionalgeo_question_result& second) {
		core::swap(first.context_id, second.context_id);
		core::swap(first.question_id, second.question_id);
		core::swap(first.answer, second.answer);
		core::swap(first.label, second.label);
	}
};

inline bool operator < (const fictionalgeo_question_result& first, const fictionalgeo_question_result& second) {
	if (first.context_id < second.context_id) return true;
	else if (first.context_id > second.context_id) return false;
	else return first.question_id < second.question_id;
}

inline bool operator == (const fictionalgeo_question_result& first, const fictionalgeo_question_result& second) {
	return first.context_id == second.context_id
		&& first.question_id == second.question_id;
}

inline bool get_answer(string& out, const array<string>& answers) {
	if (answers.length == 0) {
		return init(out, "");
	} else if (answers.length == 1) {
		unsigned int start = 0;
		unsigned int length = answers[0].length;
		if (answers[0][0] == '{' && answers[0][answers[0].length - 1] == '}') {
			start = 1;
			length -= 2;
			if (length >= 4 && answers[0][answers[0].length - 2] == '.'
			 && answers[0][answers[0].length - 3] == '.'
			 && answers[0][answers[0].length - 4] == '.'
			 && answers[0][answers[0].length - 5] == ',')
			{
				length -= 4;
			} else if (length == 3 && answers[0][answers[0].length - 2] == '.'
					&& answers[0][answers[0].length - 3] == '.'
					&& answers[0][answers[0].length - 4] == '.')
			{
				length -= 3;
			}
		}
		return init(out, answers[0].data + start, length);
	} else {
		unsigned int length = answers[0].length;
		for (unsigned int i = 1; i < answers.length; i++)
			length += 2 + answers[i].length;

		if (!init(out, length))
			return false;
		out.length = 0;
		for (unsigned int i = 0; i < answers.length; i++) {
			if (i != 0) {
				out[out.length++] = ',';
				out[out.length++] = ' ';
			}
			for (unsigned int j = 0; j < answers[i].length; j++)
				out[out.length++] = answers[i][j];
		}
		return true;
	}
}

constexpr unsigned int MAX_FICTIONALGEO_QUESTION_COUNT = 2000;

#if defined(SANITIZE_ADDRESS)
/* TODO: for memory debugging; delete this */
#include <sanitizer/lsan_interface.h>
#endif

template<bool LinearSearch, bool ParseOnly, typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
void do_fictionalgeo_experiments(bool& status,
		fictionalgeo_context_item<Theory, PriorStateType>* context_queue,
		fictionalgeo_question_item<Theory, PriorStateType>* question_queue,
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
		const array<string>& geobase,
		std::mutex& results_lock,
		array<fictionalgeo_question_result>& results,
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
			fictionalgeo_question_item<Theory, PriorStateType>& job = question_queue[question_queue_start++];
			lock.unlock();

			/* for reproducibility, reset the PRNG state */
			core::engine = context_queue[job.context_id].prng_engine;

			unsigned int parse_count;
			constexpr unsigned int max_parse_count = 3;
			hol_term* logical_forms[max_parse_count];
			double log_probabilities[max_parse_count];
			if (ParseOnly) {
				if (!parse_sentence(parser, job.question.data, names, logical_forms, log_probabilities, parse_count))
					parse_count = 0;
				results_lock.lock();
				results.ensure_capacity(results.length + 1);
				unsigned int index;
				for (index = 0; index < results.length; index++)
					if (results[index].context_id == job.context_id) break;
				if (index == results.length) {
					results[index].context_id = job.context_id;
					results[index].question_id = 0;
					if (!init(results[index].answer, "[]")) {
						free_logical_forms(logical_forms, parse_count);
						status = false;
						num_threads_running--;
						work_queue_cv.notify_all();
						free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					} else if (!init(results[index].label, 1)) {
						free(results[index].answer);
						free_logical_forms(logical_forms, parse_count);
						status = false;
						num_threads_running--;
						work_queue_cv.notify_all();
						free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
					results.length++;
				}

				memory_stream stream(results[index].answer.length + 1024);
				write(results[index].answer.data, stream, results[index].answer.length - 1);
				if (stream.position != 1)
					fputc(',', stream);
				fputc('[', stream);
				for (unsigned int i = 0; i < parse_count; i++) {
					hol_term* preprocessed = preprocess_formula(logical_forms[i]);
					if (preprocessed == nullptr) {
						free_logical_forms(logical_forms, parse_count);
						status = false;
						num_threads_running--;
						work_queue_cv.notify_all();
						free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
					array_map<unsigned int, unsigned int> variable_map(16);
					hol_term* canonicalized = Theory::FormulaCanonicalizer::canonicalize(*preprocessed, variable_map);
					core::free(*preprocessed); if (preprocessed->reference_count == 0) core::free(preprocessed);
					if (canonicalized == nullptr) {
						free_logical_forms(logical_forms, parse_count);
						status = false;
						num_threads_running--;
						work_queue_cv.notify_all();
						free(job);
						for (auto entry : names) free(entry.key);
						free(parser); return;
					}
					if (i != 0) fputc(',', stream);
					fputc('"', stream);
					memory_stream temp_stream(1024);
					print<hol_term_syntax::TPTP>(*canonicalized, temp_stream, parser.terminal_printer);
					for (unsigned int i = 0; i < temp_stream.position; i++) {
						if (temp_stream.buffer[i] == '"')
							fputc('\\', stream);
						fputc(temp_stream.buffer[i], stream);
					}
					fputc('"', stream);
					free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
				}
				fputc(']', stream);
				fputc(']', stream);
				swap(stream.buffer, results[index].answer.data);
				results[index].answer.length = stream.position;

				printf("Answer for question %u: ", job.context_id + 1);
				print(results[index].answer, stdout); printf("\n");
				results_lock.unlock();

				total++;
				free_logical_forms(logical_forms, parse_count);

			} else if (parse_sentence(parser, job.question.data, names, logical_forms, log_probabilities, parse_count)) {
				/* try to answer the question */
				array<string> best_answers(4);
				double best_answer_probability = -std::numeric_limits<double>::infinity();
#if defined(SANITIZE_ADDRESS)
/* TODO: for memory debugging; delete this */
__lsan_do_leak_check();
#endif

				double first_question_probability = -std::numeric_limits<double>::infinity();
				for (unsigned int i = 0; i < parse_count && log_probabilities[i] + 1.0 > first_question_probability; i++) {
					array<string> answers(4);
					double answer_probability = -std::numeric_limits<double>::infinity();
					if (!answer_question<LinearSearch>(answers, logical_forms[i], LinearSearch ? 40 : 400, parser.get_printer(), job.T, proof_prior, job.proof_axioms, answer_probability) || answers.length == 0)
						continue;
					bool confident = true;
					for (unsigned int j = 0; j < answers.length; j++) {
						if (answers[j] == UNKNOWN_CONCEPT_NAME) {
							confident = false;
							break;
						}
					}
					if (confident) {
						first_question_probability = max(first_question_probability, log_probabilities[i]);
						if (answer_probability > best_answer_probability) {
							for (string& str : best_answers) free(str);
							best_answers.clear();
							for (const string& answer : answers)
								best_answers.add(answer);
							best_answer_probability = answer_probability;
						}
					}
					for (string& str : answers) free(str);
				}
				if (best_answers.length == 0) {
					best_answers[0] = "<failed to answer question>";
					best_answers.length = 1;
				}
#if defined(SANITIZE_ADDRESS)
/* TODO: for memory debugging; delete this */
__lsan_do_leak_check();
#endif

				results_lock.lock();
				results.ensure_capacity(results.length + 1);
				results[results.length].context_id = job.context_id;
				results[results.length].question_id = job.question_id;
				if (!get_answer(results[results.length].answer, best_answers)) {
					for (string& str : best_answers) free(str);
					free_logical_forms(logical_forms, parse_count);
					status = false;
					num_threads_running--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				} else if (!init(results[results.length].label, job.label)) {
					for (string& str : best_answers) free(str);
					free(results[results.length].answer);
					free_logical_forms(logical_forms, parse_count);
					status = false;
					num_threads_running--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				printf("Answer for question %u: ", job.context_id + 1);
				print(results[results.length].answer, stdout); printf("\n");
				results.length++;
				results_lock.unlock();
				for (string& str : best_answers) free(str);

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
			free(job);

		} else {
			num_threads_reading_context++;
			fictionalgeo_context_item<Theory, PriorStateType>& job = context_queue[context_queue_start++];
			lock.unlock();
/*if (job.context_id != 8 - 1) {
total += job.questions.length;
num_threads_reading_context--;
free(job);
continue;
}*/

			/* for reproducibility, reset the PRNG state */
			core::engine = prng_engine;

			/* parse the list of line numbers */
			// unsigned int i = 0;
			// unsigned int start = 0;
			// unsigned int first_line = UINT_MAX;
			// array<pair<unsigned int, unsigned int>> line_numbers(8);
			// bool error = false;
			// for (; !error; i++) {
			// 	char old_char;
			// 	switch (job.context[i]) {
			// 	case '-':
			// 		/* parse the first line */
			// 		if (first_line != UINT_MAX) {
			// 			fprintf(stderr, "WARNING: Found a line number range with more than one hyphen.\n");
			// 		} else {
			// 			job.context[i] = '\0';
			// 			if (!parse_uint(string(job.context + start, i - start), first_line)) {
			// 				fprintf(stderr, "ERROR: Failed to parse line number.\n");
			// 				job.context[i] = '-'; error = true; break;
			// 			}
			// 			job.context[i] = '-';
			// 			start = i + 1;
			// 		}
			// 		break;

			// 	case ',':
			// 	case '\0':
			// 		if (!line_numbers.ensure_capacity(line_numbers.length + 1)) {
			// 			error = true;
			// 			break;
			// 		}
			// 		old_char = job.context[i];
			// 		job.context[i] = '\0';
			// 		if (!parse_uint(string(job.context + start, i - start), line_numbers[line_numbers.length].value)) {
			// 			fprintf(stderr, "ERROR: Failed to parse line number.\n");
			// 			job.context[i] = ','; error = true; break;
			// 		}
			// 		job.context[i] = old_char;
			// 		if (first_line == UINT_MAX) {
			// 			line_numbers[line_numbers.length].key = line_numbers[line_numbers.length].value;
			// 			line_numbers.length++;
			// 		} else {
			// 			line_numbers[line_numbers.length].key = first_line;
			// 			line_numbers.length++;
			// 			first_line = UINT_MAX;
			// 		}
			// 		start = i + 1;
			// 		break;
			// 	}
			// 	if (job.context[i] == '\0')
			// 		break;
			// }

			/* read the context sentences */
			unsigned int i = 0;
			unsigned int start = 0;
			array<string> context_sentences(16);
			bool error = false;
			for (; job.context[i] != '\0'; i++) {
				if (job.context[i] == '.' && !isdigit(job.context[i + 1])) {
					if (!context_sentences.ensure_capacity(context_sentences.length + 1)
					 || !init(context_sentences[context_sentences.length], job.context + start, i - start + 2))
					{
						error = true;
						break;
					}
					context_sentences[context_sentences.length][i - start + 1] = '\0';
					context_sentences[context_sentences.length].length = i - start + 1;
					context_sentences.length++;
					start = i + 1;
					while (isspace(job.context[start])) start++;
				}
			}

			if (!error) {
				/*char filename[256];
				snprintf(filename, 256, "fictionalgeo_theories/%u.th", job.context_id);*/

				unsigned int result_index;
				if (ParseOnly) {
					results_lock.lock();
					results.ensure_capacity(results.length + 1);
					for (result_index = 0; result_index < results.length; result_index++)
						if (results[result_index].context_id == job.context_id) break;
					if (result_index == results.length) {
						results[result_index].context_id = job.context_id;
						results[result_index].question_id = 0;
						if (!init(results[result_index].answer, "[]")) {
							status = false;
							num_threads_running--;
							work_queue_cv.notify_all();
							free(job);
							for (auto entry : names) free(entry.key);
							free(parser); return;
						} else if (!init(results[result_index].label, 1)) {
							free(results[result_index].answer);
							status = false;
							num_threads_running--;
							work_queue_cv.notify_all();
							free(job);
							for (auto entry : names) free(entry.key);
							free(parser); return;
						}
						results.length++;
					}
					results_lock.unlock();
				}

				/* read the context sentences */
				/*free(job.T); free(job.proof_axioms);
				FILE* theory_stream = (FILE*) fopen(filename, "rb");
				read_random_state(theory_stream);
				read(job.T, theory_stream, job.proof_axioms);
				fclose(theory_stream);
debug_terminal_printer = &parser.terminal_printer;
job.T.template print_axioms<true>(stdout, *debug_terminal_printer);*/
				unsigned int sentence_counter = 0;
				for (unsigned int i = 0; i < context_sentences.length; i++) {
//if (sentence_counter < 7) { sentence_counter++; continue; }

					unsigned int parse_count;
					constexpr unsigned int max_parse_count = 2;
					hol_term* logical_forms[max_parse_count];
					double log_probabilities[max_parse_count];
					if (ParseOnly) {
						if (!parse_sentence(parser, context_sentences[i].data, names, logical_forms, log_probabilities, parse_count))
							parse_count = 0;
						results_lock.lock();
						memory_stream stream(results[result_index].answer.length + 1024);
						write(results[result_index].answer.data, stream, results[result_index].answer.length - 1);
						if (stream.position != 1)
							fputc(',', stream);
						fputc('[', stream);
						for (unsigned int i = 0; i < parse_count; i++) {
							hol_term* preprocessed = preprocess_formula(logical_forms[i]);
							if (preprocessed == nullptr) {
								free_logical_forms(logical_forms, parse_count);
								status = false;
								num_threads_running--;
								work_queue_cv.notify_all();
								free(job);
								for (auto entry : names) free(entry.key);
								free(parser); return;
							}
							array_map<unsigned int, unsigned int> variable_map(16);
							hol_term* canonicalized = Theory::FormulaCanonicalizer::canonicalize(*preprocessed, variable_map);
							core::free(*preprocessed); if (preprocessed->reference_count == 0) core::free(preprocessed);
							if (canonicalized == nullptr) {
								free_logical_forms(logical_forms, parse_count);
								status = false;
								num_threads_running--;
								work_queue_cv.notify_all();
								free(job);
								for (auto entry : names) free(entry.key);
								free(parser); return;
							}
							if (i != 0) fputc(',', stream);
							fputc('"', stream);
							memory_stream temp_stream(1024);
							print<hol_term_syntax::TPTP>(*canonicalized, temp_stream, parser.terminal_printer);
							for (unsigned int i = 0; i < temp_stream.position; i++) {
								if (temp_stream.buffer[i] == '"')
									fputc('\\', stream);
								fputc(temp_stream.buffer[i], stream);
							}
							fputc('"', stream);
							free(*canonicalized); if (canonicalized->reference_count == 0) free(canonicalized);
						}
						fputc(']', stream);
						fputc(']', stream);
						swap(stream.buffer, results[result_index].answer.data);
						results[result_index].answer.length = stream.position;

						printf("Answer for question %u: ", job.context_id + 1);
						print(results[result_index].answer, stdout); printf("\n");
						results_lock.unlock();

						free_logical_forms(logical_forms, parse_count);
						sentence_counter++;
						continue;
					}

					// TODO: this is kind of a hacky way to get the new proof
					hash_set<nd_step<hol_term>*> old_proofs(job.T.observations.capacity);
					for (nd_step<hol_term>* term : job.T.observations)
						old_proofs.add(term);
					printf("Sentence index: %u\n", sentence_counter);
					if (!read_sentence(corpus, parser, context_sentences[i].data, job.T, names, seed_entities, proof_prior, job.proof_axioms)) {
						std::unique_lock<std::mutex> lock(results_lock);
						if (!unparseable_context.ensure_capacity(unparseable_context.length + 1)
						 || !init(unparseable_context[unparseable_context.length].value, context_sentences[i]))
						{
							status = false;
							num_threads_running--;
							num_threads_reading_context--;
							work_queue_cv.notify_all();
							free(job);
							for (auto entry : names) free(entry.key);
							free(parser); return;
						}
						unparseable_context[unparseable_context.length++].key = job.context_id;
						error = true;
						break;
					}

					nd_step<hol_term>* new_proof = nullptr;
					for (nd_step<hol_term>* proof : job.T.observations) {
						if (!old_proofs.contains(proof)) {
							new_proof = proof;
							break;
						}
					}

#if defined(SANITIZE_ADDRESS)
// TODO: for memory debugging; delete this
__lsan_do_leak_check();
#endif
					Theory& T_MAP = *((Theory*) alloca(sizeof(Theory)));
					PriorStateType& proof_axioms_MAP = *((PriorStateType*) alloca(sizeof(PriorStateType)));
					hash_map<const hol_term*, hol_term*> formula_map(128);
					Theory::clone(job.T, T_MAP, formula_map);
					PriorStateType::clone(job.proof_axioms, proof_axioms_MAP, formula_map);
					auto collector = make_log_probability_collector(job.T, proof_prior, new_proof);
					double max_log_probability = collector.current_log_probability;
timer stopwatch;
					for (unsigned int j = 0; j < 4; j++) {
						for (unsigned int t = 0; t < 150; t++) {
							fprintf(stderr, "j = %u, t = %u\n", j, t);
/*job.proof_axioms.check_proof_axioms(job.T);
job.proof_axioms.check_universal_eliminations(job.T, collector);
job.T.check_concept_axioms();
job.T.check_disjunction_introductions();
job.T.are_elements_provable();
job.T.sets.check_freeable_sets();
job.T.sets.are_descendants_valid();
job.T.sets.are_set_sizes_valid();
job.T.sets.check_set_ids();
if (!job.T.observations.contains(collector.test_proof))
fprintf(stderr, "WARNING: `log_probability_collector.test_proof` is not an observation in the theory.\n");*/
							bool print_debug = false;
							if (print_debug) job.T.template print_axioms<true>(stdout, *debug_terminal_printer);
							if (print_debug) { job.T.print_disjunction_introductions(stdout, *debug_terminal_printer); fflush(stdout); }
							do_mh_step(job.T, proof_prior, job.proof_axioms, collector, collector.test_proof, (t < 40 ? 1.0 : 0.01));

							if (collector.current_log_probability > max_log_probability + 1.0e-8) {
								free(T_MAP); free(proof_axioms_MAP); formula_map.clear();
								Theory::clone(job.T, T_MAP, formula_map);
								PriorStateType::clone(job.proof_axioms, proof_axioms_MAP, formula_map);
								max_log_probability = collector.current_log_probability;
							}
						}

						if (j + 1 < 4) {
							for (unsigned int t = 0; t < 20; t++)
								do_exploratory_mh_step(job.T, proof_prior, job.proof_axioms, collector, collector.test_proof, 1.0);
						}
					}
					free(job.T); free(job.proof_axioms); formula_map.clear();
					Theory::clone(T_MAP, job.T, formula_map);
					PriorStateType::clone(proof_axioms_MAP, job.proof_axioms, formula_map);
					T_MAP.template print_axioms<true>(stdout, *debug_terminal_printer); fflush(stdout);
					free(T_MAP); free(proof_axioms_MAP);
total_reasoning += stopwatch.milliseconds();
fprintf(stderr, "consistency checking time: %llums, total reasoning time: %llums\n", consistency_checking_ms.load(), total_reasoning.load());

					/*FILE* theory_stream = (FILE*) fopen(filename, "wb");
					write_random_state(theory_stream);
					write(job.T, theory_stream, job.proof_axioms);
					fclose(theory_stream);*/
					sentence_counter++;
#if defined(SANITIZE_ADDRESS)
// TODO: for memory debugging; delete this
__lsan_do_leak_check();
#endif
				}
			}

			if (!error) {
				/* if we successfully read the context, enqueue the jobs for reading/answering the associated questions */
				job.prng_engine = core::engine;
				std::unique_lock<std::mutex> lock(work_queue_lock);
				if (question_queue_length + job.questions.length > MAX_FICTIONALGEO_QUESTION_COUNT) {
					fprintf(stderr, "do_fictionalgeo_experiments ERROR: Requested question queue length exceeds `MAX_FICTIONALGEO_QUESTION_COUNT`.\n");
					status = false;
					num_threads_running--;
					num_threads_reading_context--;
					work_queue_cv.notify_all();
					free(job);
					for (auto entry : names) free(entry.key);
					free(parser); return;
				}
				for (unsigned int j = 0; j < job.questions.length; j++) {
					fictionalgeo_question_item<Theory, PriorStateType>& new_question = question_queue[question_queue_length];
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

inline void print_fictionalgeo_results(
		const std::atomic_uint& total,
		array<fictionalgeo_question_result>& results,
		array<pair<unsigned int, string>>& unparseable_questions,
		array<pair<unsigned int, string>>& unparseable_context,
		std::mutex& results_lock,
		const char* output_filepath,
		unsigned int total_question_count)
{
	std::unique_lock<std::mutex> lock(results_lock);
	insertion_sort(results);
	FILE* out = open_file(output_filepath, "w");
	if (out == nullptr) {
		fprintf(stderr, "ERROR: Unable to open `%s` for writing.\n", output_filepath);
		return;
	}
	unsigned int i = 0;
	for (unsigned int j = 0; j < total_question_count; j++) {
		if (i < results.length && results[i].context_id == j) {
			print(results[i].answer, out);
			print('\n', out);
			i++;
		} else {
			print("<no answer>\n", out);
		}
	}
	fclose(out);

	fprintf(stdout,
			"Results so far:\n"
			"  Total questions: %u\n"
			"  Answered questions: %zu\n",
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

template<bool LinearSearch, bool ParseOnly = false, typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
bool run_fictionalgeoqa_experiments(
		ArticleSource& corpus, Parser& parser,
		Theory& T, PriorStateType& proof_axioms,
		ProofPrior& proof_prior,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& seed_entities,
		const array<string>& geobase,
		const char* data_filepath,
		const char* results_filepath,
		unsigned int thread_count)
{
	bool status = true;
	fictionalgeo_context_item<Theory, PriorStateType>* context_queue = (fictionalgeo_context_item<Theory, PriorStateType>*) malloc(sizeof(fictionalgeo_context_item<Theory, PriorStateType>) * MAX_FICTIONALGEO_QUESTION_COUNT);
	if (context_queue == nullptr) return false;
	fictionalgeo_question_item<Theory, PriorStateType>* question_queue = (fictionalgeo_question_item<Theory, PriorStateType>*) malloc(sizeof(fictionalgeo_question_item<Theory, PriorStateType>) * MAX_FICTIONALGEO_QUESTION_COUNT);
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
	array<fictionalgeo_question_result> results(64);
	array<pair<unsigned int, string>> unparseable_questions(4);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(1);
	std::atomic_uint num_threads_running(0);

	std::thread* workers = new std::thread[thread_count];
	for (unsigned int i = 0; i < thread_count; i++) {
		workers[i] = std::thread(
				do_fictionalgeo_experiments<LinearSearch, ParseOnly, ArticleSource, Parser, Theory, PriorStateType, ProofPrior>,
				std::ref(status), context_queue, question_queue,
				std::ref(context_queue_start), std::ref(question_queue_start),
				std::ref(context_queue_length), std::ref(question_queue_length),
				std::ref(work_queue_lock), std::ref(work_queue_cv),
				prng_engine, std::ref(corpus),
				std::ref(parser), std::ref(proof_prior),
				std::ref(names), std::ref(seed_entities),
				std::ref(geobase), std::ref(results_lock),
				std::ref(results), std::ref(unparseable_questions),
				std::ref(unparseable_context), std::ref(total),
				std::ref(num_threads_reading_context),
				std::ref(num_threads_running));
	}

	unsigned int context_id = 0;
	unsigned int total_question_count = 0;
	auto process_fictionalgeo_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms,&total_question_count](char* context, array<pair<string, string>>& questions)
	{
		if (context_queue_length + 1 > MAX_FICTIONALGEO_QUESTION_COUNT) {
			fprintf(stderr, "run_fictionalgeoqa_experiments ERROR: Requested context queue length exceeds `MAX_FICTIONALGEO_QUESTION_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		fictionalgeo_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_fictionalgeoqa_experiments ERROR: Out of memory.\n");
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
		total_question_count++;
		work_queue_cv.notify_one();
		return true;
	};

	if (!read_ruletaker_data<string>(data_filepath, process_fictionalgeo_questions))
		status = false;
	num_threads_reading_context--;

	timer stopwatch;
	while (status) {
		if (num_threads_reading_context == 0 && context_queue_start == context_queue_length && question_queue_start == question_queue_length)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_fictionalgeo_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath, total_question_count);
			stopwatch.start();
		}
	}

	work_queue_cv.notify_all();
	while (status) {
		if (num_threads_running == 0)
			break;

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		if (stopwatch.milliseconds() > 1000) {
			print_fictionalgeo_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath, total_question_count);
			stopwatch.start();
		}
	}
	for (unsigned int i = 0; i < thread_count; i++) {
		if (!workers[i].joinable()) continue;
		try {
			workers[i].join();
		} catch (...) { }
	}
	print_fictionalgeo_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath, total_question_count);
	delete[] workers;
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

template<bool LinearSearch, bool ParseOnly = false, typename ArticleSource, typename Parser, typename Theory, typename PriorStateType, typename ProofPrior>
bool run_fictionalgeoqa_experiments_single_threaded(
		ArticleSource& corpus, Parser& parser,
		Theory& T, PriorStateType& proof_axioms,
		ProofPrior& proof_prior,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& seed_entities,
		const array<string>& geobase,
		const char* data_filepath,
		const char* results_filepath)
{
	bool status = true;
	fictionalgeo_context_item<Theory, PriorStateType>* context_queue = (fictionalgeo_context_item<Theory, PriorStateType>*) malloc(sizeof(fictionalgeo_context_item<Theory, PriorStateType>) * MAX_FICTIONALGEO_QUESTION_COUNT);
	if (context_queue == nullptr) return false;
	fictionalgeo_question_item<Theory, PriorStateType>* question_queue = (fictionalgeo_question_item<Theory, PriorStateType>*) malloc(sizeof(fictionalgeo_question_item<Theory, PriorStateType>) * MAX_FICTIONALGEO_QUESTION_COUNT);
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
	array<fictionalgeo_question_result> results(64);
	array<pair<unsigned int, string>> unparseable_questions(4);
	array<pair<unsigned int, string>> unparseable_context(4);
	std::atomic_uint total(0);
	std::atomic_uint num_threads_reading_context(0);
	std::atomic_uint num_threads_running(0);

	unsigned int context_id = 0;
	unsigned int total_question_count = 0;
	auto process_fictionalgeo_questions = [context_queue,&context_queue_length,&work_queue_lock,&work_queue_cv,&context_id,&T,&proof_axioms,&total_question_count](char* context, array<pair<string, string>>& questions)
	{
		if (context_queue_length + 1 > MAX_FICTIONALGEO_QUESTION_COUNT) {
			fprintf(stderr, "run_fictionalgeoqa_experiments_single_threaded ERROR: Requested context queue length exceeds `MAX_FICTIONALGEO_QUESTION_COUNT`.\n");
			return false;
		}

		std::unique_lock<std::mutex> lock(work_queue_lock);
		fictionalgeo_context_item<Theory, PriorStateType>& new_context = context_queue[context_queue_length];
		set_empty(new_context.T);
		new_context.context_id = context_id++;
		new_context.context = (char*) malloc(sizeof(char) * (strlen(context) + 1));
		if (new_context.context == nullptr) {
			fprintf(stderr, "run_fictionalgeoqa_experiments_single_threaded ERROR: Out of memory.\n");
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
		total_question_count++;
		work_queue_cv.notify_one();
		return true;
	};

	if (!read_ruletaker_data<string>(data_filepath, process_fictionalgeo_questions))
		status = false;

	do_fictionalgeo_experiments<LinearSearch, ParseOnly>(status, context_queue, question_queue,
			context_queue_start, question_queue_start, context_queue_length,
			question_queue_length, work_queue_lock, work_queue_cv, prng_engine,
			corpus, parser, proof_prior, names, seed_entities, geobase,
			results_lock, results, unparseable_questions, unparseable_context,
			total, num_threads_reading_context, num_threads_running);

	print_fictionalgeo_results(total, results, unparseable_questions, unparseable_context, results_lock, results_filepath, total_question_count);
	for (unsigned int i = context_queue_start; i < context_queue_length; i++)
		free(context_queue[i]);
	for (unsigned int i = question_queue_start; i < question_queue_length; i++)
		free(question_queue[i]);
	free(context_queue);
	free(question_queue);
	return status;
}

#endif /* FICTIONALGEOQA_H_ */
