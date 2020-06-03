#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include "article.h"
#include "natural_deduction_mh.h"

constexpr double PERPLEXITY_THRESHOLD = 0.01;

template<typename Formula>
inline void free_logical_forms(Formula** logical_forms, unsigned int count)
{
	for (unsigned int i = 0; i < count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
}

inline bool concatenate(
		sentence_token* tokens, unsigned int token_count,
		string& out, const string** reverse_string_map)
{
	unsigned int length = reverse_string_map[tokens[0].id]->length;
	for (unsigned int i = 1; i < token_count; i++)
		length += 1 + reverse_string_map[tokens[i].id]->length;

	if (!core::resize(out.data, length)) return false;
	out.length = length;

	unsigned int index = 0;
	for (unsigned int i = 0; i < token_count; i++) {
		if (i > 0) { out[index] = ' '; index++; }
		for (unsigned int j = 0; j < reverse_string_map[tokens[i].id]->length; j++) {
			out[index] = reverse_string_map[tokens[i].id]->data[j];
			index++;
		}
	}
	return true;
}

template<typename ArticleSource, typename Parser,
	typename Formula, typename Canonicalizer, typename TheoryPrior, typename Printer>
bool read_sentence(
		const ArticleSource& articles, Parser& parser, const typename Parser::SentenceType& s,
		theory<Formula, natural_deduction<Formula>, Canonicalizer>& T,
		unsigned int article_name, const hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& visited_articles, TheoryPrior& theory_prior, Printer& printer)
{
	unsigned int parse_count, new_constant;
	Formula* logical_forms[2];
	double log_probabilities[2];
	while (true) {
		/* attempt to parse the sentence */
		array<array<sentence_token>> unrecognized(16);
		print("Reading sentence: '", stdout); print(s, stdout, printer); print("'\n", stdout);
		if (!parser.template parse<2>(s, logical_forms, log_probabilities, parse_count, T, unrecognized)) {
			print("read_sentence ERROR: Unable to parse sentence '", stderr); print(s, stderr, printer); print("'.\n", stderr);
			return false;
		}

		/* concatenate the unrecognized tokens */
		array<unsigned int> unrecognized_concatenated(max((size_t) 1, unrecognized.length));
		for (unsigned int i = 0; i < unrecognized.length; i++) {
			string concatenated(16);
			if (!concatenate(unrecognized[i].data, unrecognized[i].length, concatenated, parser.reverse_name_map)) {
				for (array<sentence_token>& tokens : unrecognized) free(tokens);
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
			bool contains;
			unsigned int concatenated_id = names.get(concatenated, contains);
			if (contains)
				unrecognized_concatenated[unrecognized_concatenated.length++] = concatenated_id;
		}
		for (array<sentence_token>& tokens : unrecognized) free(tokens);

		/* get unrecognized named entities */
		array<string> named_entities(16);
		if (parse_count > 0 && !get_named_entities(*logical_forms[0], named_entities)) {
			free_logical_forms(logical_forms, parse_count);
			return false;
		}
		for (const string& named_entity : named_entities) {
			bool contains;
			unsigned int named_entity_id = names.get(named_entity, contains);
			if (contains
			 && (named_entity_id == article_name || !visited_articles.contains(named_entity_id))
			 && !unrecognized_concatenated.add(named_entity_id))
			{
				for (string& entity : named_entities) free(entity);
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
		}
		for (string& entity : named_entities) free(entity);

		/* read the article on the unrecognized word */
		if (unrecognized_concatenated.length == 0) {
			break;
		} else if (parse_count > 0 && unrecognized_concatenated.length == 1 && unrecognized_concatenated[0] == article_name) {
			/* this could be a definition so try adding it to the theory */
			bool success = (T.add_formula(logical_forms[0], new_constant) != nullptr)
						&& parser.add_definition(s, logical_forms[0], new_constant);
			if (!success) {
				print("read_sentence ERROR: Unable to add definition to theory.\n", stderr);
				print("  Sentence:     '", stderr); print(s, stderr, printer); print("'\n", stderr);
				print("  Logical form: ", stderr); print(*logical_forms[0], stderr, printer); print("\n", stderr);
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
			free_logical_forms(logical_forms, parse_count);

			//for (unsigned int t = 0; t < 10; t++)
			//	if (!do_mh_step(T, theory_prior)) return false;
			return true;
		}

		/* find an article in order to learn about the unrecognized word */
		unsigned int next_article = unrecognized_concatenated[0];
		if (next_article == article_name) {
			if (unrecognized_concatenated.length == 1) {
				fprintf(stderr, "read_sentence ERROR: Unable to parse definitional sentence.\n");
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
			next_article = unrecognized_concatenated[1];
		}
		read_article(next_article, articles, parser, T, names, visited_articles, theory_prior, printer);
		free_logical_forms(logical_forms, parse_count);
	}

	if (parse_count == 0) {
		fprintf(stderr, "read_sentence ERROR: Given sentence has no valid parses.\n");
		return false;
	} else if (parse_count > 1 && (log_probabilities[0] - log_probabilities[1]) / s.length < PERPLEXITY_THRESHOLD) {
		/* this parse is too ambiguous */
		free_logical_forms(logical_forms, parse_count);
		return true;
	}

	/* add the most probable logical form to the theory */
	if (T.add_formula(logical_forms[0], new_constant) == nullptr) {
		print("read_sentence ERROR: Unable to add logical form to theory.\n", stderr);
		print("  Sentence:     '", stderr); print(s, stderr, printer); print("'\n", stderr);
		print("  Logical form: ", stderr); print(*logical_forms[0], stderr, printer); print("\n", stderr);
		free_logical_forms(logical_forms, parse_count);
		return false;
	}

	//for (unsigned int t = 0; t < 10; t++)
	//	if (!do_mh_step(T, theory_prior)) return false;

	free_logical_forms(logical_forms, parse_count);
	return true;
}

template<typename ArticleSource, typename Parser, typename Formula,
	typename Canonicalizer, typename TheoryPrior, typename Printer>
bool read_article(
		unsigned int article_name, const ArticleSource& articles, Parser& parser,
		theory<Formula, natural_deduction<Formula>, Canonicalizer>& T,
		const hash_map<string, unsigned int>& names, hash_set<unsigned int>& visited_articles,
		TheoryPrior& theory_prior, Printer& printer)
{
	print("Reading article: '", stdout); print(article_name, stdout, printer); print("'\n", stdout);

	bool article_exists;
	const auto& doc = articles.get(article_name, article_exists);
	if (!article_exists) {
		print("read_article ERROR: No such article '", stderr); print(article_name, stderr, printer); print("'.\n", stderr);
		return false;
	} else if (!visited_articles.add(article_name)) {
		return false;
	}

	for (unsigned int i = 0; i < doc.sentence_count; i++) {
		if (!read_sentence(articles, parser, doc.sentences[i], T, article_name, names, visited_articles, theory_prior, printer))
			return false;
	}
	return true;
}

template<typename Parser, size_t ParseCount>
inline bool parse_sentence(Parser& parser,
		const typename Parser::SentenceType& sentence,
		hash_map<string, unsigned int>& names,
		hol_term* (&logical_forms)[ParseCount],
		double (&log_probabilities)[ParseCount],
		unsigned int& parse_count)
{
	array<array<sentence_token>> unrecognized(4);
	if (parser.invert_name_map(names)) {
		if (parser.template parse<ParseCount>(sentence, logical_forms, log_probabilities, parse_count, nullptr, unrecognized)) {
			for (array<sentence_token>& tokens : unrecognized) free(tokens);
			return true;
		} else {
			fprintf(stderr, "ERROR: Parsing failed.\n");
			return false;
		}
	} else {
		fprintf(stderr, "ERROR: `invert_name_map` failed.\n");
		return false;
	}
}

template<typename Parser>
inline bool parse_sentence(Parser& parser, const char* input_sentence, hash_map<string, unsigned int>& names)
{
	typename Parser::SentenceType sentence;
	if (!tokenize(input_sentence, sentence, names))
		return false;

	constexpr unsigned int max_parse_count = 4;
	hol_term* logical_forms[max_parse_count];
	double log_probabilities[max_parse_count];
	unsigned int parse_count;
	if (!parse_sentence(parser, sentence, names, logical_forms, log_probabilities, parse_count))
		return false;
	for (unsigned int i = 0; i < parse_count; i++) {
		string_map_scribe terminal_printer = { parser.reverse_name_map, names.table.size + 1 };
		print(*logical_forms[i], stderr, terminal_printer); print(" with log probability ", stderr); print(log_probabilities[i], stderr); print('\n', stderr);
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
	return true;
}

template<typename Parser, typename Formula>
inline bool parse_sentence(
		Parser& parser, const char* input_sentence, hash_map<string, unsigned int>& names,
		array<array_map<typename Parser::SentenceType, Formula>>& training_set)
{
	typename Parser::SentenceType sentence;
	if (!tokenize(input_sentence, sentence, names))
		return false;

	constexpr unsigned int max_parse_count = 4;
	hol_term* logical_forms[max_parse_count];
	double log_probabilities[max_parse_count];
	unsigned int parse_count = 0;
	bool result = parse_sentence(parser, sentence, names, logical_forms, log_probabilities, parse_count);

	string_map_scribe terminal_printer = { parser.reverse_name_map, names.table.size + 1 };
extern string_map_scribe* debug_terminal_printer;
debug_terminal_printer = &terminal_printer;
	for (const auto& paragraph : training_set) {
		bool found_training_sentence = false;
		for (const auto& entry : paragraph) {
			sequence expected_sentence(nullptr, 0);
			if (is_empty(entry.key.derivation)) {
				expected_sentence.tokens = (unsigned int*) malloc(sizeof(unsigned int) * entry.key.length);
				if (expected_sentence.tokens == nullptr) {
					fprintf(stderr, "parse_sentence ERROR: Out of memory.\n");
					result = false; found_training_sentence = true;
					break;
				}
				for (unsigned int i = 0; i < entry.key.length; i++)
					expected_sentence.tokens[i] = entry.key.tokens[i].id;
				expected_sentence.length = entry.key.length;
			} else {
				if (!yield(parser.G, *entry.key.derivation.tree, entry.value, expected_sentence, terminal_printer, entry.key.derivation.root))
					continue;
			}

			if (expected_sentence.length != sentence.length) {
				free(expected_sentence);
				continue;
			}
			bool are_sentences_identical = true;
			for (unsigned int i = 0; i < sentence.length && are_sentences_identical; i++)
				if (sentence.tokens[i].id != expected_sentence.tokens[i]) are_sentences_identical = false;
			if (!are_sentences_identical) {
				free(expected_sentence);
				continue;
			}
			free(expected_sentence);

			found_training_sentence = true;
			if (parse_count == 0) {
				fprintf(stderr, "parse_sentence WARNING: Unable to parse sentence '%s' despite being in the training data.\n", input_sentence);
				break;
			} if (*logical_forms[0] != *entry.value.root) {
				fprintf(stderr, "parse_sentence WARNING: The parsed logical form does not match the label logical form in the training data:\n");
				fprintf(stderr, "  Sentence: '%s'\n", input_sentence);
				print("  Parsed logical form:   ", stderr); print(*logical_forms[0], stderr, terminal_printer); print('\n', stderr);
				print("  Expected logical form: ", stderr); print(*entry.value.root, stderr, terminal_printer); print('\n', stderr);
				if (!is_empty(entry.key.derivation)) {
					double expected_log_probability = log_probability(parser.G, *entry.key.derivation.tree, entry.value, parser, entry.key.derivation.root);
					print("  with log probability ", stderr); print(expected_log_probability, stderr); print('\n', stderr);
				}
			}
			break;
		}
		if (found_training_sentence) break;
	}

	for (unsigned int i = 0; i < parse_count; i++) {
		print(*logical_forms[i], stderr, terminal_printer); print(" with log probability ", stderr); print(log_probabilities[i], stderr); print('\n', stderr);
	}

	for (unsigned int i = 0; i < parse_count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
	free(sentence);
	return result;
}

template<bool LookupUnknownWords, typename ArticleSource,
	typename Formula, typename ProofCalculus, typename Canonicalizer,
	typename TheoryPrior, typename Parser, typename Printer>
inline bool answer_question(array<string>& answers,
		const char* input_question, unsigned int num_samples,
		const ArticleSource& articles, Parser& parser,
		theory<Formula, ProofCalculus, Canonicalizer>& T,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& visited_articles,
		TheoryPrior& theory_prior, Printer& printer)
{
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;
	typedef typename ProofCalculus::Proof Proof;

	typename Parser::SentenceType sentence;
	if (!tokenize(input_question, sentence, names))
		return false;

	unsigned int parse_count;
	constexpr unsigned int max_parse_count = 2;
	Formula* logical_forms[max_parse_count];
	double log_parse_probabilities[max_parse_count];
	while (true) {
		if (!parser.invert_name_map(names)) {
			fprintf(stderr, "ERROR: `invert_name_map` failed.\n");
			free(sentence); return false;
		}

		/* attempt to parse the sentence */
		array<array<sentence_token>> unrecognized(4);
		print("Reading question '", stdout); print(sentence, stdout, printer); print("'\n", stdout);
		if (!parser.template parse<max_parse_count>(sentence, logical_forms, log_parse_probabilities, parse_count, nullptr, unrecognized)) {
			fprintf(stderr, "ERROR: Unable to parse question.\n");
			free(sentence); return false;
		}

		if (LookupUnknownWords)
		{
			/* concatenate the unrecognized tokens */
			array<unsigned int> unrecognized_concatenated(max((size_t) 1, unrecognized.length));
			for (unsigned int i = 0; i < unrecognized.length; i++) {
				string concatenated(16);
				if (!concatenate(unrecognized[i].data, unrecognized[i].length, concatenated, parser.reverse_name_map)) {
					for (array<sentence_token>& tokens : unrecognized) free(tokens);
					free_logical_forms(logical_forms, parse_count);
					free(sentence); return false;
				}
				bool contains;
				unsigned int concatenated_id = names.get(concatenated, contains);
				if (contains)
					unrecognized_concatenated[unrecognized_concatenated.length++] = concatenated_id;
			}
			for (array<sentence_token>& tokens : unrecognized) free(tokens);

			/* get unrecognized named entities */
			array<string> named_entities(16);
			if (parse_count > 0 && !get_named_entities(*logical_forms[0], named_entities)) {
				free_logical_forms(logical_forms, parse_count);
				free(sentence); return false;
			}
			for (const string& named_entity : named_entities) {
				bool contains;
				unsigned int named_entity_id = names.get(named_entity, contains);
				if (contains && !visited_articles.contains(named_entity_id)
				&& !unrecognized_concatenated.add(named_entity_id))
				{
					for (string& entity : named_entities) free(entity);
					free_logical_forms(logical_forms, parse_count);
					free(sentence); return false;
				}
			}
			for (string& entity : named_entities) free(entity);

			if (unrecognized_concatenated.length == 0)
				break;
			free_logical_forms(logical_forms, parse_count);

			/* find an article in order to learn about the unrecognized word */
			unsigned int next_article = unrecognized_concatenated[0];
			read_article(next_article, articles, parser, T, names, visited_articles, theory_prior, printer);
		} else {
			for (array<sentence_token>& tokens : unrecognized) free(tokens);
			break;
		}
	}
	free(sentence);

	array_map<Term, double> log_probabilities(1024);
	if (!log_joint_probability_of_lambda(T, theory_prior, logical_forms[0], num_samples, log_probabilities)) {
		fprintf(stderr, "ERROR: Failed to answer question.\n");
		free_logical_forms(logical_forms, parse_count);
		return false;
	}
	free_logical_forms(logical_forms, parse_count);

	/* get the name of each term */
	array_map<string, double> temp_answers(8);
	for (const auto& entry : log_probabilities) {
		const Term* term = &entry.key;
		if (term->type == TermType::STRING) {
			if (!temp_answers.ensure_capacity(temp_answers.size + 1)) {
				for (auto pair : log_probabilities) free(pair.key);
				for (auto pair : temp_answers) free(pair.key);
				return false;
			}
			unsigned int index = temp_answers.index_of(term->str);
			if (index < temp_answers.size) {
				temp_answers.values[index] = logsumexp(temp_answers.values[index], entry.value);
			} else {
				if (!init(temp_answers.keys[index], term->str)) {
					for (auto pair : log_probabilities) free(pair.key);
					for (auto pair : temp_answers) free(pair.key);
					return false;
				}
				temp_answers.values[index] = logsumexp(temp_answers.values[index], entry.value);
				temp_answers.size++;
			}
		} else if (term->type == TermType::INTEGER) {
			int length = snprintf(NULL, 0, "%d", term->integer);
			if (!temp_answers.ensure_capacity(temp_answers.size + 1) || length < 0) {
				for (auto pair : log_probabilities) free(pair.key);
				for (auto pair : temp_answers) free(pair.key);
				return false;
			}

			string new_name(length + 1);
			snprintf(new_name.data, length + 1, "%d", term->integer);
			new_name.length = length;

			unsigned int index = temp_answers.index_of(new_name);
			if (index < temp_answers.size) {
				temp_answers.values[index] = logsumexp(temp_answers.values[index], entry.value);
			} else {
				if (!init(temp_answers.keys[index], new_name)) {
					for (auto pair : log_probabilities) free(pair.key);
					for (auto pair : temp_answers) free(pair.key);
					return false;
				}
				temp_answers.values[index] = logsumexp(temp_answers.values[index], entry.value);
				temp_answers.size++;
			}
		} else if (term->type == TermType::CONSTANT) {
			bool contains;
			for (Proof* definition : T.ground_concepts[term->constant - T.new_constant_offset].definitions) {
				if (definition->formula->binary.right->type != TermType::UNARY_APPLICATION
				 || definition->formula->binary.right->binary.left->type != TermType::CONSTANT
				 || definition->formula->binary.right->binary.left->constant != (unsigned int) built_in_predicates::ARG1
				 || definition->formula->binary.right->binary.right->type != TermType::CONSTANT)
					continue;

				unsigned int event = definition->formula->binary.right->binary.right->constant;
				if (!T.ground_concepts[event - T.new_constant_offset].types.contains((unsigned int) built_in_predicates::NAME))
					continue;

				Proof* function_value_axiom = T.ground_concepts[event - T.new_constant_offset].function_values.get((unsigned int) built_in_predicates::ARG2, contains);
				if (!contains || function_value_axiom->formula->binary.right->type != TermType::STRING)
					continue;

				unsigned int index = temp_answers.index_of(function_value_axiom->formula->binary.right->str);
				if (index < temp_answers.size) {
					temp_answers.values[index] = logsumexp(temp_answers.values[index], entry.value);
				} else {
					if (!init(temp_answers.keys[index], function_value_axiom->formula->binary.right->str)) {
						for (auto pair : log_probabilities) free(pair.key);
						for (auto pair : temp_answers) free(pair.key);
						return false;
					}
					temp_answers.values[index] = entry.value;
					temp_answers.size++;
				}
			}
		} else {
			fprintf(stderr, "ERROR: Unable to convert semantic answer into text.\n");
			for (auto pair : log_probabilities) free(pair.key);
			for (auto pair : temp_answers) free(pair.key);
			return false;
		}
	}
	for (auto pair : log_probabilities) free(pair.key);

	/* keep only the answers with highest probability */
	double max_log_probability = -std::numeric_limits<double>::infinity();
	for (const auto& pair : temp_answers) {
		if (pair.value > max_log_probability) {
			for (string& str : answers) free(str);
			answers.clear();
			if (!init(answers[0], pair.key)) {
				for (auto pair : temp_answers) free(pair.key);
				for (string& str : answers) free(str);
				return false;
			}
			answers.length = 1;
			max_log_probability = pair.value;
		} else if (pair.value == max_log_probability) {
			if (!answers.ensure_capacity(answers.length + 1)
			 || !init(answers[answers.length], pair.key))
			{
				for (auto pair : temp_answers) free(pair.key);
				for (string& str : answers) free(str);
				return false;
			}
			answers.length++;
		}
	}
	for (auto pair : temp_answers) free(pair.key);
	return true;
}

#endif /* EXECUTIVE_H_ */
