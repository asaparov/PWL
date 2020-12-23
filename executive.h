#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include "article.h"
#include "natural_deduction_mh.h"

constexpr double PERPLEXITY_THRESHOLD = 0.0; //0.01;
constexpr double SUFFICIENT_KNOWLEDGE_THRESHOLD = 8.0;

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
	typename Formula, typename Canonicalizer, typename TheoryPrior, typename... Args>
bool read_sentence(
		const ArticleSource& articles, Parser& parser, const typename Parser::SentenceType& s,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		unsigned int article_name, hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& visited_articles, TheoryPrior& theory_prior,
		typename TheoryPrior::PriorState& proof_axioms,
		unsigned int mcmc_iterations_per_retry = 100,
		unsigned int max_retries = 0, Args&&... add_formula_args)
{
	null_collector collector;
	unsigned int parse_count, new_constant;
	Formula* logical_forms[2];
	double log_probabilities[2];
	while (true) {
		/* attempt to parse the sentence */
		array<array<sentence_token>> unrecognized(16);
		print("Reading sentence: '", stdout); print(s, stdout, parser.get_printer()); print("'\n", stdout);
		if (!parser.template parse<2>(s, logical_forms, log_probabilities, parse_count, T, unrecognized, names)) {
			print("read_sentence ERROR: Unable to parse sentence '", stderr); print(s, stderr, parser.get_printer()); print("'.\n", stderr);
			return false;
		}

		/* concatenate the unrecognized tokens */
		array<unsigned int> unrecognized_concatenated(max((size_t) 1, unrecognized.length));
		for (unsigned int i = 0; i < unrecognized.length; i++) {
			string concatenated(16);
			if (!concatenate(unrecognized[i].data, unrecognized[i].length, concatenated, parser.reverse_string_map())) {
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

		/* preemptively remove tokens for which an article doesn't exist in the corpus */
		for (unsigned int i = 0; i < unrecognized_concatenated.length; i++) {
			if (!articles.contains(unrecognized_concatenated[i])) {
				if (!visited_articles.add(unrecognized_concatenated[i])) {
					free_logical_forms(logical_forms, parse_count);
					return false;
				}
				unrecognized_concatenated.remove(i--);
			}
		}

		/* read the article on the unrecognized word */
		if (unrecognized_concatenated.length == 0) {
			break;
		} else if (parse_count > 0 && unrecognized_concatenated.length == 1 && unrecognized_concatenated[0] == article_name) {
			/* this could be a definition so try adding it to the theory */
			set_changes<Formula> set_diff;
			auto* new_proof = T.add_formula(logical_forms[0], set_diff, new_constant, std::forward<Args>(add_formula_args)...);
			for (unsigned int i = 0; new_proof == nullptr && i < max_retries; i++) {
				set_diff.clear();
				for (unsigned int t = 0; t < mcmc_iterations_per_retry; t++)
					do_exploratory_mh_step(T, theory_prior, proof_axioms, collector);
				new_proof = T.add_formula(logical_forms[0], set_diff, new_constant, std::forward<Args>(add_formula_args)...);
			}
			if (new_proof != nullptr) {
				if (proof_axioms.add(new_proof, set_diff.new_set_axioms, theory_prior)) {
					if (!parser.add_definition(s, logical_forms[0], new_constant, names)) {
						proof_axioms.subtract(new_proof, set_diff.new_set_axioms, theory_prior);
						T.remove_formula(new_proof, set_diff);
						new_proof = nullptr;
/* TODO: for debugging; delete this */
exit(EXIT_FAILURE);
					}
				} else {
					T.remove_formula(new_proof, set_diff);
					new_proof = nullptr;
				}
			}
			if (new_proof == nullptr) {
				print("read_sentence ERROR: Unable to add definition to theory.\n", stderr);
				print("  Sentence:     '", stderr); print(s, stderr, parser.get_printer()); print("'\n", stderr);
				print("  Logical form: ", stderr); print(*logical_forms[0], stderr, parser.get_printer()); print("\n", stderr);
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
			free_logical_forms(logical_forms, parse_count);

			//for (unsigned int t = 0; t < 10; t++)
			//	do_mh_step(T, theory_prior, proof_axioms, collector);
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
		read_article(next_article, articles, parser, T, names, visited_articles, theory_prior, proof_axioms, std::forward<Args>(add_formula_args)...);
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
	set_changes<Formula> set_diff;
	auto* new_proof = T.add_formula(logical_forms[0], set_diff, new_constant, std::forward<Args>(add_formula_args)...);
	for (unsigned int i = 0; new_proof == nullptr && i < max_retries; i++) {
		set_diff.clear();
auto collector = make_log_probability_collector(T, theory_prior);
		for (unsigned int t = 0; t < mcmc_iterations_per_retry; t++)
			do_exploratory_mh_step(T, theory_prior, proof_axioms, collector);
		new_proof = T.add_formula(logical_forms[0], set_diff, new_constant, std::forward<Args>(add_formula_args)...);
	}
	if (new_proof != nullptr && !proof_axioms.add(new_proof, set_diff.new_set_axioms, theory_prior)) {
		T.remove_formula(new_proof, set_diff);
		new_proof = nullptr;
	}
	if (new_proof == nullptr) {
		print("read_sentence ERROR: Unable to add logical form to theory.\n", stderr);
		print("  Sentence:     '", stderr); print(s, stderr, parser.get_printer()); print("'\n", stderr);
		print("  Logical form: ", stderr); print(*logical_forms[0], stderr, parser.get_printer()); print("\n", stderr);
		free_logical_forms(logical_forms, parse_count);
		return false;
	}

	//for (unsigned int t = 0; t < 10; t++)
	//	do_mh_step(T, theory_prior, proof_axioms, collector);

	free_logical_forms(logical_forms, parse_count);
	return true;
}

template<typename ArticleSource, typename Parser,
	typename Formula, typename Canonicalizer, typename TheoryPrior, typename... Args>
bool read_article(
		unsigned int article_name, const ArticleSource& articles, Parser& parser,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		hash_map<string, unsigned int>& names, hash_set<unsigned int>& visited_articles,
		TheoryPrior& theory_prior, typename TheoryPrior::PriorState& proof_axioms,
		unsigned int mcmc_iterations_per_retry = 100,
		unsigned int max_retries = 0, Args&&... add_formula_args)
{
	print("Reading article: '", stdout); print(article_name, stdout, parser.get_printer()); print("'\n", stdout);

	bool article_exists;
	const auto& doc = articles.get(article_name, article_exists);
	if (!visited_articles.add(article_name)) {
		return false;
	} else if (!article_exists) {
		print("read_article ERROR: No such article '", stderr); print(article_name, stderr, parser.get_printer()); print("'.\n", stderr);
		return false;
	}

	for (unsigned int i = 0; i < doc.sentence_count; i++) {
		if (!read_sentence(articles, parser, doc.sentences[i], T, article_name, names, visited_articles, theory_prior, proof_axioms, mcmc_iterations_per_retry, max_retries, std::forward<Args>(add_formula_args)...))
			return false;
	}
	return true;
}

template<typename ArticleSource, typename Parser,
	typename Formula, typename Canonicalizer, typename TheoryPrior, typename... Args>
inline bool read_sentence(
		const ArticleSource& articles, Parser& parser, const char* input_sentence,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		hash_map<string, unsigned int>& names, hash_set<unsigned int>& visited_articles,
		TheoryPrior& theory_prior, typename TheoryPrior::PriorState& proof_axioms,
		unsigned int mcmc_iterations_per_retry = 100,
		unsigned int max_retries = 0, Args&&... add_formula_args)
{
	typename Parser::SentenceType sentence;
	if (!tokenize(input_sentence, sentence, names)
	 || !parser.invert_name_map(names))
		return false;

	bool result = read_sentence(articles, parser, sentence, T, UINT_MAX, names, visited_articles, theory_prior, proof_axioms, mcmc_iterations_per_retry, max_retries, std::forward<Args>(add_formula_args)...);
	free(sentence);
	return result;
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
		if (parser.template parse<ParseCount>(sentence, logical_forms, log_probabilities, parse_count, nullptr, unrecognized, names)) {
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

template<typename Parser, size_t ParseCount>
inline bool parse_sentence(Parser& parser,
		const char* input_sentence,
		hash_map<string, unsigned int>& names,
		hol_term* (&logical_forms)[ParseCount],
		double (&log_probabilities)[ParseCount],
		unsigned int& parse_count)
{
	typename Parser::SentenceType sentence;
	if (!tokenize(input_sentence, sentence, names))
		return false;

	if (!parse_sentence(parser, sentence, names, logical_forms, log_probabilities, parse_count)) {
		free(sentence);
		return false;
	}
	free(sentence);
	return true;
}

template<typename Parser, typename ProofCalculus, typename Canonicalizer, typename TheoryPrior>
inline bool parse_sentence_with_prior(Parser& parser, const char* input_sentence,
		theory<ProofCalculus, Canonicalizer>& T, hash_map<string, unsigned int>& names,
		TheoryPrior& theory_prior, typename TheoryPrior::PriorState& proof_axioms)
{
	constexpr unsigned int max_parse_count = 6;
	hol_term* logical_forms[max_parse_count];
	double log_likelihoods[max_parse_count];
	unsigned int parse_count;
	if (!parse_sentence(parser, input_sentence, names, logical_forms, log_likelihoods, parse_count))
		return false;

	double log_priors[max_parse_count];
	double log_posteriors[max_parse_count];
	for (unsigned int i = 0; i < parse_count; i++) {
		log_priors[i] = log_joint_probability_of_observation(T, theory_prior, proof_axioms, logical_forms[i], 100000);
		log_posteriors[i] = log_likelihoods[i] + log_priors[i];
	}
	unsigned int sorted_indices[max_parse_count];
	for (unsigned int i = 0; i < parse_count; i++)
		sorted_indices[i] = i;
	if (parse_count > 1) {
		insertion_sort(log_posteriors, sorted_indices, parse_count);
		reverse(log_posteriors, parse_count);
		reverse(sorted_indices, parse_count);
	}
	for (unsigned int i = 0; i < parse_count; i++) {
		print(*logical_forms[sorted_indices[i]], stderr, parser.get_printer());
		print(" with log likelihood ", stderr); print(log_likelihoods[sorted_indices[i]], stderr);
		print(" + log prior ", stderr); print(log_priors[sorted_indices[i]], stderr);
		print(" = log posterior ", stderr); print(log_posteriors[i], stderr);
		print('\n', stderr);
		free(*logical_forms[sorted_indices[i]]);
		if (logical_forms[sorted_indices[i]]->reference_count == 0)
			free(logical_forms[sorted_indices[i]]);
	}
	return true;
}

template<typename Parser>
inline bool parse_sentence(Parser& parser, const char* input_sentence, hash_map<string, unsigned int>& names)
{
	constexpr unsigned int max_parse_count = 4;
	hol_term* logical_forms[max_parse_count];
	double log_probabilities[max_parse_count];
	unsigned int parse_count;
	if (!parse_sentence(parser, input_sentence, names, logical_forms, log_probabilities, parse_count))
		return false;
	for (unsigned int i = 0; i < parse_count; i++) {
		print(*logical_forms[i], stderr, parser.get_printer()); print(" with log probability ", stderr); print(log_probabilities[i], stderr); print('\n', stderr);
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

extern const string_map_scribe* debug_terminal_printer;
debug_terminal_printer = &parser.get_printer();
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
				if (!yield(parser.G, *entry.key.derivation.tree, entry.value, expected_sentence, parser.get_printer(), entry.key.derivation.root))
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
			}

			if (!equivalent(logical_forms[0], entry.value.root)) {
				fprintf(stderr, "parse_sentence WARNING: The parsed logical form does not match the label logical form in the training data:\n");
				fprintf(stderr, "  Sentence: '%s'\n", input_sentence);
				print("  Parsed logical form:   ", stderr); print(*logical_forms[0], stderr, parser.get_printer()); print('\n', stderr);
				print("  Expected logical form: ", stderr); print(*entry.value.root, stderr, parser.get_printer()); print('\n', stderr);
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
		print(*logical_forms[i], stderr, parser.get_printer()); print(" with log probability ", stderr); print(log_probabilities[i], stderr); print('\n', stderr);
	}

	if (parse_count != 0) {
		unsigned int generated_derivation_count;
		constexpr unsigned int max_generated_derivation_count = 12;
		double log_likelihoods[max_generated_derivation_count];
		syntax_node<typename Parser::logical_form_type>* generated_derivations =
				(syntax_node<typename Parser::logical_form_type>*) alloca(sizeof(syntax_node<typename Parser::logical_form_type>) * max_generated_derivation_count);
		if (!parser.template generate<max_generated_derivation_count>(generated_derivations, log_likelihoods, generated_derivation_count, logical_forms[0], names) || generated_derivation_count == 0) {
			fprintf(stderr, "parse_sentence ERROR: Failed to generate derivation.\n");
			for (unsigned int i = 0; i < parse_count; i++) {
				free(*logical_forms[i]);
				if (logical_forms[i]->reference_count == 0)
					free(logical_forms[i]);
			}
			free(sentence);
			return false;
		}

		const string** nonterminal_name_map = invert(parser.G.nonterminal_names);
		string_map_scribe nonterminal_printer = { nonterminal_name_map, parser.G.nonterminal_names.table.size + 1 };
		fprintf(stderr, "Generated %u derivations:\n", generated_derivation_count);
		for (unsigned int i = 0; i < generated_derivation_count; i++) {
			/* compute the yield of the derivation */
			sequence new_sentence = sequence(NULL, 0);
			if (!parser.yield(generated_derivations[i], logical_forms[0], new_sentence)) {
				for (unsigned int i = 0; i < parse_count; i++) {
					free(*logical_forms[i]);
					if (logical_forms[i]->reference_count == 0)
						free(logical_forms[i]);
				} for (unsigned int i = 0; i < generated_derivation_count; i++)
					free(generated_derivations[i]);
				free(sentence); free(nonterminal_name_map);
				return false;
			}
			print('"', stderr); print(new_sentence, stderr, parser.get_printer());
			print("\" with log likelihood ", stderr); print(log_likelihoods[i], stderr);
			print(" and derivation tree:\n", stderr);
			print(generated_derivations[i], stderr, nonterminal_printer, parser.get_printer());
			print('\n', stderr);
			free(new_sentence);
		}
		free(nonterminal_name_map);
		for (unsigned int i = 0; i < generated_derivation_count; i++)
			free(generated_derivations[i]);
	}

	for (unsigned int i = 0; i < parse_count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
	free(sentence);
	return result;
}

#include "network.h"
#include <fcntl.h>
#include <openssl/ssl.h>
#include <openssl/err.h>
#include <openssl/conf.h>

inline ssize_t read(SSL* ssl, void* buf, size_t nbytes) {
	return SSL_read(ssl, buf, nbytes);
}

template<unsigned int BufferSize, typename Stream>
struct buffered_stream {
	Stream& underlying_stream;
	unsigned long long timeout_ms;
	array<char> buffer;
	size_t position;
	bool eof;

	buffered_stream(Stream& underlying_stream, unsigned long long timeout_ms) :
			underlying_stream(underlying_stream), timeout_ms(timeout_ms), buffer(4096), position(0), eof(false) { }

	bool fill_buffer() {
		if (!buffer.ensure_capacity(buffer.length + BufferSize))
			return false;
		timer stopwatch;
		while (true) {
			int bytes_read = read(underlying_stream, buffer.data + buffer.length, BufferSize);
			if (bytes_read == 0) {
				eof = true;
				return true;
			} else if (bytes_read > 0) {
				buffer.length += bytes_read;
				return true;
			} else if (errno != EAGAIN && errno != EWOULDBLOCK) {
				fprintf(stderr, "buffered_stream.fill_buffer ERROR: `read` failed.\n");
				return false;
			}
			if (stopwatch.milliseconds() > timeout_ms)
				return false;
		}
	}
};

template<unsigned int BufferSize, typename Stream>
int fgetc(buffered_stream<BufferSize, Stream>& in) {
	if (in.position == in.buffer.length && (!in.fill_buffer() || in.eof))
		return EOF;
	return in.buffer[in.position++];
}

template<typename Stream>
bool read_line(Stream& in, array<char>& line)
{
	while (true) {
		int next = fgetc(in);
		if (next == EOF) {
			return false;
		} else if (next == '\r') {
			next = fgetc(in);
			if (next == '\n') {
				break;
			} else {
				if (!line.add('\r') || !line.add(next)) return false;
			}
		} else {
			if (!line.add(next)) return false;
		}
	}
	return true;
}

struct throttler
{
	array_map<string, unsigned long long> next_request_times;

	throttler() : next_request_times(16) { }

	~throttler() {
		for (auto entry : next_request_times) free(entry.key);
	}

	inline void prune_old_request_times()
	{
		unsigned long long current_time = milliseconds();
		for (unsigned int i = 0; i < next_request_times.size; i++) {
			if (current_time >= next_request_times.values[i]) {
				free(next_request_times.keys[i]);
				next_request_times.remove_at(i--);
			}
		}
	}

	inline unsigned long long get_next_request_time(const string& hostname) const
	{
		unsigned int index = next_request_times.index_of(hostname);
		if (index < next_request_times.size)
			return next_request_times.values[index];
		else return 0;
	}

	inline bool set_next_request_time(const string& hostname, unsigned long long next_request_time) {
		return next_request_times.put(hostname, next_request_time);
	}
};

/* TODO: make this thread-safe */
throttler GLOBAL_THROTTLER;
constexpr unsigned long long THROTTLE_DURATION_MILLISECONDS = 400;

template<typename Stream>
bool parse_http_status(Stream& in,
		string& version, unsigned int& status, string& reason)
{
	/* read the first line */
	array<char> line(64);
	if (!read_line(in, line)) return false;

	/* first try tokenizing the line */
	unsigned int index = line.index_of(' ');
	if (index == line.length) {
		fprintf(stderr, "parse_http_status ERROR: HTTP response status is malformed.\n");
		return false;
	} else if (!init(version, line.data, index)) {
		return false;
	}

	/* make sure the version starts with 'HTTP/' */
	static string VERSION_PREFIX("HTTP/");
	if (!compare_strings(VERSION_PREFIX, version.data, min(version.length, VERSION_PREFIX.length))) {
		fprintf(stderr, "parse_http_status ERROR: HTTP response version string is malformed.\n");
		return false;
	}

	unsigned int second_index = line.index_of(' ', index + 1);
	if (second_index == line.length) {
		if (!init(reason, 0)) {
			free(version);
			return false;
		}
	} else {
		if (!init(reason, line.data + second_index + 1, line.length - second_index - 1)) {
			free(version);
			return false;
		}
	}

	if (!parse_uint(string(line.data + index + 1, second_index - index - 1), status)) {
		fprintf(stderr, "parse_http_status ERROR: Unable to parse HTTP response code.\n");
		free(version); free(reason); return false;
	} else if (status < 100 || status > 999) {
		fprintf(stderr, "parse_http_status ERROR: Illegal HTTP response code.\n");
		free(version); free(reason); return false;
	}
	return true;
}

template<typename Stream, typename FilterResponseHeader>
bool parse_http_response(Stream& in, array<char>& payload,
		FilterResponseHeader filter_response_header)
{
	unsigned int status;
	string& version = *((string*) alloca(sizeof(string)));
	string& reason = *((string*) alloca(sizeof(string)));

	if (!parse_http_status(in, version, status, reason))
		return false;
	/* TODO: handle status 100 (CONTINUE) */
	free(version); free(reason);

	/* parse the headers */
	array_map<string, string> headers(16);
	while (true) {
		array<char> line(64);
		if (!read_line(in, line)) return false;
		if (line.length == 0) break;

		if (!headers.ensure_capacity(headers.size + 1)) {
			for (auto entry : headers) { free(entry.key); free(entry.value); }
			return false;
		}

		unsigned int index = line.index_of(':');
		for (unsigned int i = 0; i < index; i++)
			line[i] = tolower(line[i]);

		unsigned int second_start = index + 1;
		while (second_start < line.length && isspace(line[second_start])) second_start++;
		while (second_start < line.length && isspace(line[line.length - 1])) line.length--;

		if (!init(headers.keys[headers.size], line.data, index)) {
			for (auto entry : headers) { free(entry.key); free(entry.value); }
			return false;
		} else if (!init(headers.values[headers.size], line.data + second_start, line.length - second_start)) {
			free(headers.keys[headers.size]);
			for (auto entry : headers) { free(entry.key); free(entry.value); }
			return false;
		}
		headers.size++;
	}

	if (!filter_response_header(status, headers)) {
		for (auto entry : headers) { free(entry.key); free(entry.value); }
		return true;
	}

	/* check if the transfer encoding is chunked */
	static const string TRANSFER_ENCODING("transfer-encoding");
	bool contains, chunked = false;
	string& transfer_encoding = headers.get(TRANSFER_ENCODING, contains);
	if (contains) {
		for (unsigned int i = 0; i < transfer_encoding.length; i++)
			transfer_encoding[i] = tolower(transfer_encoding[i]);
		static const string CHUNKED("chunked");
		if (transfer_encoding == CHUNKED)
			chunked = true;
	}

	/* get the content length, if provided */
	unsigned int content_length;
	if (chunked) {
		content_length = UINT_MAX;
	} else if (status >= 100 && status < 200) {
		/* TODO: also check if `status` is NO_CONTENT or NOT_MODIFIED, or the request was "HEAD" */
		content_length = 0;
	} else {
		static const string CONTENT_LENGTH("content-length");
		string& content_length_str = headers.get(CONTENT_LENGTH, contains);
		if (contains) {
			if (!parse_uint(content_length_str, content_length)) {
				fprintf(stderr, "parse_http_response WARNING: Invalid content-length value.\n");
				content_length = UINT_MAX;
			}
		} else {
			content_length = UINT_MAX;
		}
	}
	for (auto entry : headers) { free(entry.key); free(entry.value); }

	/* read the payload */
	if (chunked) {
		unsigned int chunk_bytes_remaining = UINT_MAX;
		while (true) {
			if (chunk_bytes_remaining == 0) {
				/* discard the CLRF at the end of the chunk */
				fgetc(in); fgetc(in);
				chunk_bytes_remaining = UINT_MAX;
			} if (chunk_bytes_remaining == UINT_MAX) {
				/* read the next chunk size */
				array<char> line(64);
				if (!read_line(in, line)) return false;
				unsigned int index = line.index_of(';');
				line.length = index;

				if (!parse_uint(line, chunk_bytes_remaining, 16)) {
					fprintf(stderr, "parse_http_response WARNING: Invalid chunk length.\n");
					return false;
				}
				if (chunk_bytes_remaining == 0)
					break;
			}

			/* read the chunk */
			bool eof = false;
			for (; chunk_bytes_remaining > 0; chunk_bytes_remaining--) {
				int next = fgetc(in);
				if (next == EOF) {
					fprintf(stderr, "parse_http_response WARNING: Payload is smaller than chunk length.\n");
					eof = true; break;
				}
				if (!payload.add(next)) return false;
			}
			if (eof) break;
		}
	} else {
		for (unsigned int i = 0; i < content_length; i++) {
			int next = fgetc(in);
			if (next == EOF) {
				fprintf(stderr, "parse_http_response WARNING: Payload is smaller than content-length.\n");
				break;
			}
			if (!payload.add(next)) return false;
		}
	}
	return true;
}

struct ssl_context_provider {
	SSL_CTX* ctx;
	SSL_CONF_CTX* conf;

	ssl_context_provider() {
		ERR_load_crypto_strings();
		SSL_load_error_strings();

		conf = SSL_CONF_CTX_new();
		SSL_CONF_CTX_set_flags(conf, SSL_CONF_FLAG_CLIENT);

		ctx = SSL_CTX_new(TLS_client_method());
		//SSL_CTX_clear_mode(ctx, SSL_MODE_AUTO_RETRY);
		SSL_CONF_CTX_set_ssl_ctx(conf, ctx);
		if (!SSL_CONF_CTX_finish(conf)) ERR_print_errors_fp(stderr);
		SSL_CTX_set_min_proto_version(ctx, TLS1_2_VERSION);
		//SSL_CTX_set_verify(ctx, SSL_VERIFY_NONE, nullptr);
        //if (SSL_CTX_set_default_verify_file(ctx) <= 0
		// || SSL_CTX_set_default_verify_dir(ctx) <= 0)
		//	ERR_print_errors_fp(stderr);
		//SSL_CTX_set_session_cache_mode(ctx, SSL_SESS_CACHE_CLIENT | SSL_SESS_CACHE_NO_INTERNAL_STORE);
	}

	~ssl_context_provider() {
		SSL_CTX_free(ctx);
		SSL_CONF_CTX_free(conf);
	}
};

ssl_context_provider GLOBAL_SSL_PROVIDER;

void print_openssl_error(SSL* ssl, int ret, bool& can_shutdown)
{
	char msg[1024];
	can_shutdown = true;
	switch (SSL_get_error(ssl, ret)) {
	case SSL_ERROR_ZERO_RETURN: fprintf(stderr, "TLS/SSL peer has closed the connection.\n"); break;
	case SSL_ERROR_WANT_READ: fprintf(stderr, "Read operation did not complete.\n"); break;
	case SSL_ERROR_WANT_WRITE: fprintf(stderr, "Write operation did not complete.\n"); break;
	case SSL_ERROR_WANT_CONNECT: fprintf(stderr, "Connect operation did not complete.\n"); break;
	case SSL_ERROR_WANT_ACCEPT: fprintf(stderr, "Accept operation did not complete.\n"); break;
	case SSL_ERROR_WANT_X509_LOOKUP: fprintf(stderr, "X509 lookup operation did not complete.\n"); break;
	case SSL_ERROR_WANT_ASYNC: fprintf(stderr, "Asynchronous engine is still processing data.\n"); break;
	case SSL_ERROR_WANT_ASYNC_JOB: fprintf(stderr, "No asynchronous jobs are available.\n"); break;
	case SSL_ERROR_WANT_CLIENT_HELLO_CB: fprintf(stderr, "Client callback operation did not complete.\n"); break;
	case SSL_ERROR_SYSCALL: fprintf(stderr, "I/O error.\n"); can_shutdown = false; break;
	case SSL_ERROR_SSL:
		ERR_error_string_n(ERR_get_error(), msg, sizeof(msg));
		fprintf(stderr, "%s %s %s %s.\n", msg, ERR_lib_error_string(0), ERR_func_error_string(0), ERR_reason_error_string(0));
		can_shutdown = false; break;
	default: fprintf(stderr, "Unknown error.\n"); break;
	}
}

template<bool UseSSL, typename FilterResponseHeader>
bool get_http_page(
		const char* hostname, const char* query, const char* port,
		unsigned long long timeout_ms, array<char>& response,
		FilterResponseHeader filter_response_header)
{
	static constexpr int MAX_REQUEST_LEN = 1024;
	static char REQUEST_TEMPLATE[] =
			"GET %s HTTP/1.1\r\n"
			"Host: %s\r\n"
			"Accept: text/html,application/xhtml+xml,application/xml,text/plain\r\n"
			"Accept-Encoding: identity\r\n"
			"User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.108 Safari/537.36\r\n"
			"Connection: keep-alive\r\n"
			"DNT: 1\r\n\r\n";
	int request_length = snprintf(NULL, 0, REQUEST_TEMPLATE, query, hostname);
	if (request_length >= MAX_REQUEST_LEN) {
		fprintf(stderr, "get_http_page ERROR: Request length is at least `MAX_REQUEST_LENGTH`.\n");
		return false;
	}

	char* request = (char*) malloc(sizeof(char) * request_length + 1);
	if (request == nullptr) {
		fprintf(stderr, "get_http_page ERROR: Out of memory.\n");
		return false;
	}
	snprintf(request, request_length + 1, REQUEST_TEMPLATE, query, hostname);

	bool success = true;
	auto process_connection = [&response,&success,hostname,request,request_length,filter_response_header,timeout_ms](socket_type& connection)
	{
		/* make the underlying socket non-blocking */
		int flags;
		while ((flags = fcntl(connection.handle, F_GETFL)) == -1 && errno == EINTR) { }
		if (flags == -1) {
			fprintf(stderr, "get_http_page ERROR: Failed to make socket non-blocking; %s.\n", strerror(errno));
			shutdown(connection.handle, 2);
			success = false; return false;
		}

		int rv;
		while ((rv = fcntl(connection.handle, F_SETFL, flags | O_NONBLOCK)) == -1 && errno == EINTR) { }
		if (rv != 0) {
			fprintf(stderr, "get_http_page ERROR: Failed to make socket non-blocking; %s.\n", strerror(errno));
			shutdown(connection.handle, 2);
			success = false; return false;
		}

		SSL* ssl;
		if (UseSSL) {
			ssl = SSL_new(GLOBAL_SSL_PROVIDER.ctx);
			if (ssl == nullptr) {
				fprintf(stderr, "get_http_page ERROR: SSL_new failed; ");
				char msg[1024];
				ERR_error_string_n(ERR_get_error(), msg, sizeof(msg));
				fprintf(stderr, "%s %s %s %s.\n", msg, ERR_lib_error_string(0), ERR_func_error_string(0), ERR_reason_error_string(0));
				shutdown(connection.handle, 2);
				success = false; return false;
			}
			if (!SSL_set_tlsext_host_name(ssl, hostname)) {
				fprintf(stderr, "get_http_page ERROR: SSL_set_tlsext_host_name failed; ");
				ERR_print_errors_fp(stderr);
				SSL_free(ssl); shutdown(connection.handle, 2);
				success = false; return false;
			}
			int ret = SSL_set_fd(ssl, connection.handle);
			if (ret != 1) {
				fprintf(stderr, "get_http_page ERROR: SSL_set_fd failed; ");

				bool can_shutdown;
				print_openssl_error(ssl, ret, can_shutdown);
				if (can_shutdown) SSL_shutdown(ssl);
				SSL_free(ssl); shutdown(connection.handle, 2);
				success = false; return false;
			}
			while (true) {
				ret = SSL_connect(ssl);
				if (ret == 1)
					break;
				int error = SSL_get_error(ssl, ret);
				if (error != SSL_ERROR_WANT_READ && error != SSL_ERROR_WANT_WRITE) {
					fprintf(stderr, "get_http_page ERROR: SSL_connect failed; ");

					bool can_shutdown;
					print_openssl_error(ssl, ret, can_shutdown);
					if (can_shutdown) SSL_shutdown(ssl);
					SSL_free(ssl); shutdown(connection.handle, 2);
					success = false; return false;
				}
			}
		}

		/* send HTTP request */
		unsigned int total_written = 0;
		while (total_written < (unsigned int) request_length) {
			int written;
			if (UseSSL)
				written = SSL_write(ssl, request + total_written, request_length - total_written);
			else written = write(connection.handle, request + total_written, request_length - total_written);
			if (written == -1) {
				fprintf(stderr, "get_http_page ERROR: Failed to send HTTP request.\n");
				if (UseSSL) { SSL_shutdown(ssl); SSL_free(ssl); }
				shutdown(connection.handle, 2); success = false; return false;
			}
			total_written += written;
		}

		/* read the response */
		if (UseSSL) {
			buffered_stream<BUFSIZ, SSL*> in(ssl, timeout_ms);
			success = parse_http_response(in, response, filter_response_header);
			SSL_shutdown(ssl); SSL_free(ssl);
		} else {
			buffered_stream<BUFSIZ, decltype(connection.handle)> in(connection.handle, timeout_ms);
			success = parse_http_response(in, response, filter_response_header);
		}
		shutdown(connection.handle, 2);
		return success;
	};

	unsigned long long current_time = milliseconds();
	unsigned long long next_request_time = GLOBAL_THROTTLER.get_next_request_time(string(hostname));
	if (next_request_time > current_time)
		std::this_thread::sleep_for(std::chrono::milliseconds(next_request_time - current_time));
	if (!GLOBAL_THROTTLER.set_next_request_time(string(hostname), milliseconds() + (THROTTLE_DURATION_MILLISECONDS / 2) + sample_uniform(THROTTLE_DURATION_MILLISECONDS))
	 || !run_client(hostname, port, process_connection))
		success = false;
	free(request);
	return success;
}

template<typename ProcessResultFunction>
bool search_google(
		const string* query, unsigned int query_length,
		unsigned long long timeout_ms,
		unsigned int start, bool& has_next,
		ProcessResultFunction process_result)
{
	memory_stream query_stream(1024);
	if (!print("/search?q=%22", query_stream)) return false;
	if (!print(query[0], query_stream)) return false;
	for (unsigned int i = 1; i < query_length; i++) {
		if (!print("+", query_stream)
		 || !print(query[i], query_stream)) return false;
	}
	if (start == 0) {
		if (!print("%22", query_stream)) return false;
	} else {
		if (fprintf(query_stream, "%22&start=%u", start) <= 0) return false;
	}
	if (!write('\0', query_stream)) return false;

	array<char> response(4096);
	if (!get_http_page<true>("www.google.com", query_stream.buffer, "443", timeout_ms, response, [](unsigned int status, const array_map<string, string>& headers) { return true; })) {
		print("search_google ERROR: Unable to retrieve webpage at '", stderr);
		print(query_stream.buffer, stderr); print("'.\n", stderr);
		return false;
	}

	/* find all links beginning with '<a href="' */
	array<string> results(10);
	unsigned int i; has_next = false;
	static string URL_PREFIX = "<div class=\"r\"><a href=\"";
	static string NEXT_PAGE = "Next</span></a>";
	for (i = 0; i < response.length; i++) {
		if (compare_strings(URL_PREFIX, response.data + i, min((size_t) URL_PREFIX.length, response.length - i))) {
			/* find a closing quote */
			unsigned int j;
			for (j = i + URL_PREFIX.length; j < response.length; j++)
				if (response[j] == '"') break;

			if (j == response.length) {
				fprintf(stderr, "search_google ERROR: Unexpected end of HTML response.\n");
				for (string& result : results) free(result);
				return false;
			} else if (!results.ensure_capacity(results.length + 1)
					|| !init(results[results.length], response.data + i + URL_PREFIX.length, j - i - URL_PREFIX.length))
			{
				for (string& result : results) free(result);
				return false;
			}
			results.length++;
			i = j;
		} else if (compare_strings(NEXT_PAGE, response.data + i, min((size_t) NEXT_PAGE.length, response.length - i))) {
			has_next = true;
			break;
		}
	}

	/* find all links beginning with '<div class="*"><a href=""' */
	static string DIV_CLASS_PREFIX = "<div class=\"";
	static string A_HREF_PREFIX = "\"><a href=\"";
	static string LOGO_STRING = "logo";
	static string SEARCH_URL_PREFIX = "/search?";
	for (i = 0; i < response.length; i++) {
		if (compare_strings(DIV_CLASS_PREFIX, response.data + i, min((size_t) DIV_CLASS_PREFIX.length, response.length - i))) {
			/* find a closing quote */
			unsigned int j;
			for (j = i + DIV_CLASS_PREFIX.length; j < response.length; j++)
				if (response[j] == '"') break;

			/* make sure the div class is not "logo" */
			if (compare_strings(LOGO_STRING, response.data + i + DIV_CLASS_PREFIX.length, j - i - DIV_CLASS_PREFIX.length)
			 || !compare_strings(A_HREF_PREFIX, response.data + j, min((size_t) A_HREF_PREFIX.length, response.length - j)))
			{
				i = j;
				continue;
			}
			j += A_HREF_PREFIX.length;

			/* make sure the url does not begin with "/search?" */
			if (compare_strings(SEARCH_URL_PREFIX, response.data + j, min((size_t) SEARCH_URL_PREFIX.length, response.length - j))) {
				i = j;
				continue;
			}

			/* find a closing quote */
			unsigned int url_start = j;
			for (; j < response.length; j++)
				if (response[j] == '"') break;

			if (j == response.length) {
				fprintf(stderr, "search_google ERROR: Unexpected end of HTML response.\n");
				for (string& result : results) free(result);
				return false;
			} else if (!results.ensure_capacity(results.length + 1)
					|| !init(results[results.length], response.data + url_start, j - url_start))
			{
				for (string& result : results) free(result);
				return false;
			}
			results.length++;
			i = j;
		}
	}

	for (const string& result : results)
		if (!process_result(result)) { has_next = false; break; }
	for (string& result : results) free(result);
	return true;
}

template<typename ProcessResultFunction>
bool search_google(
		const string* query, unsigned int query_length,
		unsigned long long timeout_ms,
		ProcessResultFunction process_result)
{
	bool has_next = true;
	for (unsigned int page = 0; has_next; page++) {
		print("Retrieving Google search results with query: \"", stdout);
		print(query[0], stdout);
		for (unsigned int i = 1; i < query_length; i++) {
			print(' ', stdout);
			print(query[i], stdout);
		}
		print('"', stdout);
		if (page > 0) {
			print(" (page: ", stdout); print(page + 1, stdout); print(")\n", stdout);
		} else {
			print('\n', stdout);
		}

		if (!search_google(query, query_length, timeout_ms, page * 10, has_next, process_result)) return false;
	}
	return true;
}

bool parse_path(const string& url, string& path, string& query, string& fragment, unsigned int start = 0)
{
	unsigned int position = start;
	unsigned int question_index = position + index_of('?', url.data + position, url.length - position);
	unsigned int hash_index = position + index_of('#', url.data + position, url.length - position);
	if (question_index < hash_index) {
		if (!init(path, url.data + position, question_index - position))
			return false;
		position = question_index + 1;
		if (!init(query, url.data + position, hash_index - position)) {
			free(path);
			return false;
		}
		position = min(hash_index + 1, url.length);
		if (!init(fragment, url.data + position, url.length - position)) {
			free(path); free(query);
			return false;
		}
	} else {
		if (!init(path, url.data + position, hash_index - position)) {
			return false;
		} else if (!init(query, 0)) {
			free(path);
			return false;
		}
		position = min(hash_index + 1, url.length);
		if (!init(fragment, url.data + position, url.length - position)) {
			free(path); free(query);
			return false;
		}
	}
	return true;
}

bool parse_url(const string& url,
		string& scheme, string& userinfo, string& hostname,
		string& port, string& path, string& query, string& fragment)
{
	unsigned int position = url.index_of(':');
	if (position == url.length
	 || !init(scheme, url.data, position)) return false;

	if (url.length > position + 2
	 && url[position + 1] == '/'
	 && url[position + 2] == '/')
	{
		position += 3;
		unsigned int at_index = position + index_of('@', url.data + position, url.length - position);
		if (at_index == url.length) {
			if (!init(userinfo, 0)) { free(scheme); return false; }
		} else {
			if (!init(userinfo, url.data + position, at_index - position)) { free(scheme); return false; }
			position = at_index + 1;
		}

		unsigned int colon_index = position + index_of(':', url.data + position, url.length - position);
		unsigned int slash_index = position + index_of('/', url.data + position, url.length - position);
		if (colon_index < slash_index) {
			/* there is a port number */
			if (!init(hostname, url.data + position, colon_index - position)) { free(scheme); free(userinfo); return false; }
			position = colon_index + 1;
			if (!init(port, url.data + position, slash_index - position)) { free(scheme); free(userinfo); free(hostname); return false; }
		} else {
			/* there is no port number */
			if (!init(hostname, url.data + position, slash_index - position)) { free(scheme); free(userinfo); return false; }
			if (!init(port, 0)) { free(scheme); free(userinfo); free(hostname); return false; }
		}
		position = slash_index;

	} else {
		if (!init(userinfo, 0) || !init(hostname, 0) || !init(port, 0)) {
			free(scheme);
			return false;
		}
		position++;
	}

	if (!parse_path(url, path, query, fragment, position)) {
		free(scheme); free(userinfo);
		free(hostname); free(port); return false;
	}
	return true;
}

bool get_website(const string& address, array<char>& response, unsigned long long timeout_ms = 5000)
{
	string& scheme = *((string*) alloca(sizeof(string)));
	string& userinfo = *((string*) alloca(sizeof(string)));
	string& hostname = *((string*) alloca(sizeof(string)));
	string& port = *((string*) alloca(sizeof(string)));
	string& path = *((string*) alloca(sizeof(string)));
	string& query = *((string*) alloca(sizeof(string)));
	string& fragment = *((string*) alloca(sizeof(string)));
	if (!parse_url(address, scheme, userinfo, hostname, port, path, query, fragment)) {
		print("find_answer_in_website ERROR: Unable to parse URL '\n", stderr);
		print(address, stderr); print("'.\n", stderr); return true;
	}
	free(query); free(fragment);

	if ((scheme != "http" && scheme != "https") || userinfo.length != 0) {
		free(scheme); free(userinfo); free(hostname);
		free(port); free(path); return true;
	}
	bool use_ssl = (scheme == "https");

	if (port.length == 0) {
		free(port);
		if (!init(port, use_ssl ? "443" : "80")) {
			free(scheme); free(userinfo);
			free(hostname); free(path);
			return false;
		}
	}

	hostname += " "; hostname[--hostname.length] = '\0';
	path += " "; path[--path.length] = '\0';
	port += " "; port[--port.length] = '\0';

	static const string CONTENT_TYPE("content-type");
	static const string LOCATION("location");
	static const string TEXT_HTML("text/html");
	static const string TEXT_PLAIN("text/plain");
	string& redirect = *((string*) alloca(sizeof(string)));
	redirect.data = nullptr; redirect.length = 0;
	auto filter_response_header = [&redirect](unsigned int status, const array_map<string, string>& headers) {
		bool contains;
		if (status == 400) {
			return false;
		} else if (status == 301 || status == 302) {
			const string& location = headers.get(LOCATION, contains);
			if (!contains || location.length == 0) return false;
			return init(redirect, location);
		}
		const string& content_type = headers.get(CONTENT_TYPE, contains);
		if (!contains) return true;
		unsigned int end = content_type.index_of(';');
		return compare_strings(TEXT_HTML, content_type.data, end)
			|| compare_strings(TEXT_PLAIN, content_type.data, end);
	};

	if ((use_ssl && !get_http_page<true>(hostname.data, path.data, port.data, timeout_ms, response, filter_response_header))
	 || (!use_ssl && !get_http_page<false>(hostname.data, path.data, port.data, timeout_ms, response, filter_response_header)))
	{
		print("find_answer_in_website ERROR: Unable to retrieve webpage at '", stderr);
		print(address, stderr); print("'.\n", stderr);
		free(scheme); free(userinfo); free(hostname);
		free(path); free(port); return true;
	}

	while (redirect.length != 0) {
		response.clear();
		if (redirect[0] == '/') {
			/* the URL is relative */
			free(path);
			if (!parse_path(redirect, path, query, fragment, 0)) {
				free(scheme); free(userinfo);
				free(hostname); free(port); return false;
			}
			free(query); free(fragment);
			free(redirect); redirect.length = 0;

			path += " "; path[--path.length] = '\0';

			if ((use_ssl && !get_http_page<true>(hostname.data, path.data, port.data, timeout_ms, response, filter_response_header))
			 || (!use_ssl && !get_http_page<false>(hostname.data, path.data, port.data, timeout_ms, response, filter_response_header)))
			{
				print("find_answer_in_website ERROR: Unable to retrieve webpage at '", stderr);
				print(address, stderr); print("'.\n", stderr);
				free(scheme); free(userinfo); free(hostname);
				free(path); free(port); return true;
			}
		} else {
			if (!get_website(redirect, response, timeout_ms)) {
				free(scheme); free(userinfo); free(hostname);
				free(path); free(port); return false;
			}
			free(redirect); redirect.length = 0;
		}
	}
	free(scheme); free(userinfo);
	free(hostname); free(path); free(port);
	return true;
}

inline unsigned int eat_whitespace(const char* input, unsigned int length, unsigned int position)
{
	while (position < length && isspace(input[position] == ' ')) position++;
	return position;
}

inline unsigned int get_next_token(const char* input, unsigned int length, unsigned int position)
{
	while (position < length || isspace(input[position] == ' ') || input[position] == '<') position++;
	return position;
}

enum class html_lexer_state {
	DEFAULT,
	TOKEN,
	TAG
};

struct html_lexer_token {
	html_lexer_state type;
	const char* text;
	unsigned int length;
};

template<typename EmitSentenceFunction>
bool find_answer_in_website(const string& address,
		unsigned long long timeout_ms,
		const sequence& left_question,
		const sequence& right_question,
		hash_map<string, unsigned int>& names,
		EmitSentenceFunction emit_sentence)
{
	print("Searching for sentence in website '", stdout); print(address, stdout); print("'.\n", stdout);

	array<char> response(4096);
	if (!get_website(address, response, timeout_ms))
		return false;

	unsigned int start = 0;
	html_lexer_state state = html_lexer_state::DEFAULT;
	array<html_lexer_token> tokens(4096);
	hash_set<unsigned int> fake_periods(64); /* a list of periods that do not demarcate the end of a sentence */
	for (unsigned int position = 0; position < response.length; position++) {
		switch (response[position]) {
		case ' ':
		case '\t':
		case '\r':
		case '\n':
			if (state == html_lexer_state::TOKEN) {
				/* emit the token from response[start:position] */
				if (!tokens.add({html_lexer_state::TOKEN, response.data + start, position - start})) return false;
				state = html_lexer_state::DEFAULT;
			}
			break;

		case '<':
			if (state == html_lexer_state::TOKEN) {
				/* emit the token from response[start:position] */
				if (!tokens.add({html_lexer_state::TOKEN, response.data + start, position - start})) return false;
				state = html_lexer_state::TAG;
				start = position;
			} else if (state == html_lexer_state::DEFAULT) {
				state = html_lexer_state::TAG;
				start = position;
			}
			break;

		case '>':
			if (state == html_lexer_state::TAG) {
				char* text = response.data + start;
				unsigned int length = position - start + 1;

				unsigned int start = 1;

				/* get the tag name */
				unsigned int name_end = start;
				while (start < length && !isspace(text[name_end]) && text[name_end] != '>') {
					text[name_end] = tolower(text[name_end]);
					name_end++;
				}
				const char* tag_name = text + start;
				unsigned int tag_name_length = name_end - start;

				/* ignore all <a>, </a>, <span>, and </span> tags */
				if (!(tag_name_length == 1 && tag_name[0] == 'a') && !(tag_name_length == 2 && tag_name[0] == '/' && tag_name[1] == 'a')
				 && !(tag_name_length == 4 && tag_name[0] == 's' && tag_name[1] == 'p' && tag_name[2] == 'a' && tag_name[3] == 'n')
				 && !(tag_name_length == 5 && tag_name[0] == '/' && tag_name[1] == 's' && tag_name[2] == 'p' && tag_name[3] == 'a' && tag_name[4] == 'n'))
				{
					/* emit the tag from response[start:position] */
					if (!tokens.add({html_lexer_state::TAG, text, length})) return false;
				}
				state = html_lexer_state::DEFAULT;
			}
			break;

		case '.':
		case '?':
		case '!':
		case ',':
		case ':':
		case '-':
		case '(':
		case ')':
		case '%':
		case '"':
		case '\'':
		case '{':
		case '}':
		case '[':
		case ']':
			if (state == html_lexer_state::TOKEN) {
				if (position + 1 < response.length && response[position] == '.' && isdigit(response[position - 1]) && isdigit(response[position + 1]))
					if (!fake_periods.add(tokens.length + 1)) return false;
				/* emit the token from response[start:position] */
				if (!tokens.add({html_lexer_state::TOKEN, response.data + start, position - start})
				 || !tokens.add({html_lexer_state::TOKEN, response.data + position, 1})) return false;
				state = html_lexer_state::DEFAULT;
			} else if (state == html_lexer_state::DEFAULT) {
				if (!tokens.add({html_lexer_state::TOKEN, response.data + position, 1})) return false;
				state = html_lexer_state::DEFAULT;
			}
			break;

		default:
			if (state == html_lexer_state::DEFAULT) {
				state = html_lexer_state::TOKEN;
				start = position;
			}
			break;
		}
	}

	if (state == html_lexer_state::TOKEN) {
		if (!tokens.add({html_lexer_state::TOKEN, response.data + start, (unsigned int) (response.length - start)})) return false;
	}

	if (tokens.length < left_question.length + right_question.length)
		return true;

	/* find all matches of `left_question` */
	array<unsigned int> left_match_indices(4);
	for (unsigned int i = 0; i < tokens.length - left_question.length; i++) {
		/* make sure this is the beginning of a sentence */
		if (tokens[i].type != html_lexer_state::TOKEN)
			continue;
		if (i > 0 && tokens[i - 1].type == html_lexer_state::TOKEN
		 && (tokens[i - 1].text[0] != '.' || fake_periods.contains(i - 1))
		 && tokens[i - 1].text[0] != '?' && tokens[i - 1].text[0] != '!')
			continue;

		bool match = true;
		for (unsigned int j = 0; j < left_question.length; j++) {
			bool contains;
			string key(tokens[i + j].text, tokens[i + j].length);
			unsigned int token_id = names.get(key, contains);
			if (!contains || token_id != left_question[j]) {
				match = false;
				break;
			}
		}

		if (match && !left_match_indices.add(i))
			return false;
	}

	/* find all matches of `right_question` */
	array<unsigned int> right_match_indices(4);
	for (unsigned int i = 1; i < tokens.length - right_question.length; i++) {
		bool match = true;
		for (unsigned int j = 0; j < right_question.length; j++) {
			bool contains;
			string key(tokens[i + j].text, tokens[i + j].length);
			unsigned int token_id = names.get(key, contains);
			if (!contains || token_id != right_question[j]) {
				match = false;
				break;
			}
		}

		/* make sure this is the end of a sentence */
		if (tokens[i + right_question.length - 1].type != html_lexer_state::TOKEN)
			continue;
		if (tokens[i + right_question.length].type == html_lexer_state::TOKEN
		 && (tokens[i + right_question.length].text[0] != '.' || fake_periods.contains(i + right_question.length)))
			continue;

		if (match && !right_match_indices.add(i))
			return false;
	}

	unsigned int i = 0;
	if (!tokens.ensure_capacity(tokens.length + 1)) return false;
	for (unsigned int left_index : left_match_indices) {
		/* find the smallest `right_match_index[i]` larger than `left_index` */
		while (i < right_match_indices.length && right_match_indices[i] <= left_index + left_question.length) i++;
		if (i == right_match_indices.length) {
			/* there is no such `right_match_index[i]` */
			break;
		}

		/* make sure there are no boundaries between `left_index` and `right_match_indices[i]` */
		bool is_sentence_valid = true;
		for (unsigned int j = left_index + left_question.length; j < right_match_indices[i]; j++) {
			if (tokens[j].type == html_lexer_state::TAG || (tokens[j].type == html_lexer_state::TOKEN
			 && (tokens[j].text[0] == '.' || tokens[j].text[0] == '?' || tokens[j].text[0] == '!')))
			{
				is_sentence_valid = false;
				break;
			}
		}

		if (!is_sentence_valid) continue;

		/* we found a valid matching sentence, so construct a `sentence` struct */
		html_lexer_token old_token = tokens[right_match_indices[i] + right_question.length];
		tokens[right_match_indices[i] + right_question.length].text = ".";
		tokens[right_match_indices[i] + right_question.length].length = 1;
		if (!emit_sentence(tokens.data + left_index, right_match_indices[i] + right_question.length - left_index + 1))
			return false;
		tokens[right_match_indices[i] + right_question.length] = old_token;
	}
	return true;
}

template<typename ArticleSource, typename Parser,
	typename Formula, typename Canonicalizer, typename TheoryPrior>
bool read_sentence(
		const ArticleSource& articles, Parser& parser,
		const html_lexer_token* tokens, unsigned int length,
		theory<natural_deduction<Formula>, Canonicalizer>& T,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& visited_articles, TheoryPrior& theory_prior,
		typename TheoryPrior::PriorState& proof_axioms,
		unsigned int mcmc_iterations_per_retry = 100,
		unsigned int max_retries = 0)
{
	typename Parser::SentenceType sentence;
	if (!tokenize(tokens, length, sentence, names)) {
		return false;
	} else if (!parser.invert_name_map(names)) {
		free(sentence);
		return false;
	}
	bool result = read_sentence(articles, parser, sentence, T, UINT_MAX, names, visited_articles, theory_prior, proof_axioms, mcmc_iterations_per_retry, max_retries);
	free(sentence);
	return result;
}

struct string_array_sorter { };

inline bool less_than(
		const array<string>& first,
		const array<string>& second,
		const string_array_sorter& sorter)
{
	unsigned int i = 0;
	while (i < first.length && i < second.length) {
		if (first[i] < second[i]) return true;
		else if (second[i] < first[i]) return false;
		i++;
	}
	if (i < second.length) return true;
	else if (i < first.length) return false;
	return false;
}

#include <grammar/grammar.h>

template<bool LookupUnknownWords, typename ArticleSource,
	typename ProofCalculus, typename Canonicalizer,
	typename TheoryPrior, typename Parser, typename... Args>
inline bool answer_question(array<string>& answers,
		const char* input_question, unsigned int num_samples,
		const ArticleSource& articles, Parser& parser,
		theory<ProofCalculus, Canonicalizer>& T,
		hash_map<string, unsigned int>& names,
		hash_set<unsigned int>& visited_articles,
		TheoryPrior& theory_prior,
		typename TheoryPrior::PriorState& proof_axioms,
		Args&&... add_formula_args)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename Formula::Term Term;
	typedef typename Formula::TermType TermType;

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
		print("Reading question '", stdout); print(sentence, stdout, parser.get_printer()); print("'\n", stdout);
		if (!parser.template parse<max_parse_count>(sentence, logical_forms, log_parse_probabilities, parse_count, nullptr, unrecognized, names)) {
			fprintf(stderr, "ERROR: Unable to parse question.\n");
			free(sentence); return false;
		}

		if (LookupUnknownWords)
		{
			/* concatenate the unrecognized tokens */
			array<unsigned int> unrecognized_concatenated(max((size_t) 1, unrecognized.length));
			for (unsigned int i = 0; i < unrecognized.length; i++) {
				string concatenated(16);
				if (!concatenate(unrecognized[i].data, unrecognized[i].length, concatenated, parser.reverse_string_map())) {
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
			read_article(next_article, articles, parser, T, names, visited_articles, theory_prior, proof_axioms, std::forward<Args>(add_formula_args)...);
		} else {
			for (array<sentence_token>& tokens : unrecognized) free(tokens);
			break;
		}
	}
	free(sentence);

	constexpr const char* UNKNOWN_CONCEPT_NAME = "<unknown concept>";

	array_map<string, double> temp_answers(8);
	auto on_new_proof_sample = [&T, &temp_answers, &parser, UNKNOWN_CONCEPT_NAME](const Term* term, double log_probability)
	{
		/* get the name of the term */
		if (term->type == TermType::STRING) {
			if (!temp_answers.ensure_capacity(temp_answers.size + 1))
				return;
			unsigned int index = temp_answers.index_of(term->str);
			if (index < temp_answers.size) {
				temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
			} else {
				if (!init(temp_answers.keys[index], term->str))
					return;
				temp_answers.values[index] = log_probability;
				temp_answers.size++;
			}
		} else if (term->type == TermType::NUMBER) {
			int length;
			if (term->number.decimal == 0)
				length = snprintf(NULL, 0, "%" PRId64, term->number.integer);
			else length = snprintf(NULL, 0, "%" PRId64 ".%" PRIu64, term->number.integer, term->number.decimal);
			if (!temp_answers.ensure_capacity(temp_answers.size + 1) || length < 0)
				return;

			string new_name(length + 1);
			if (term->number.decimal == 0)
				snprintf(new_name.data, length + 1, "%" PRId64, term->number.integer);
			else snprintf(new_name.data, length + 1, "%" PRId64 ".%" PRIu64, term->number.integer, term->number.decimal);
			new_name.length = length;

			unsigned int index = temp_answers.index_of(new_name);
			if (index < temp_answers.size) {
				temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
			} else {
				if (!init(temp_answers.keys[index], new_name))
					return;
				temp_answers.values[index] = log_probability;
				temp_answers.size++;
			}
		} else if (term->type == TermType::CONSTANT) {
			/* check if the constant is named */
			bool named_constant_or_set;
			if (T.new_constant_offset > term->constant) {
				const string_map_scribe& printer = parser.get_printer();
				if (printer.length > term->constant) {
					named_constant_or_set = true;
					if (!temp_answers.ensure_capacity(temp_answers.size + 1))
						return;
					unsigned int index = temp_answers.index_of(*printer.map[term->constant]);
					if (index < temp_answers.size) {
						temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
					} else {
						if (!init(temp_answers.keys[index], *printer.map[term->constant]))
							return;
						temp_answers.values[index] = log_probability;
						temp_answers.size++;
					}
				} else {
					named_constant_or_set = false;
				}
			} else {
				array<Term*> name_terms(2);
				if (!T.get_concept_names(term->constant, name_terms))
					return;
				named_constant_or_set = (name_terms.length != 0);
				for (Term* name_term : name_terms) {
print(term->constant, stderr); print(": \"", stderr);
print(name_term->str, stderr);
print("\", log probability: ", stderr); print(log_probability, stderr); print('\n', stderr);
					if (!temp_answers.ensure_capacity(temp_answers.size + 1))
						return;
					unsigned int index = temp_answers.index_of(name_term->str);
					if (index < temp_answers.size) {
						temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
					} else {
						if (!init(temp_answers.keys[index], name_term->str))
							return;
						temp_answers.values[index] = log_probability;
						temp_answers.size++;
					}
				}
			}

			/* check if the constant is a set */
			Term* set_formula = Term::new_apply(Term::new_constant(term->constant), Term::new_variable(1));
			if (set_formula == nullptr) return;
			bool contains;
			unsigned int set_id = T.sets.set_ids.get(*set_formula, contains);
			free(*set_formula); free(set_formula);
			if (contains) {
				hash_set<tuple> provable_elements(16);
				if (!T.sets.get_provable_elements(set_id, provable_elements))
					return;
				array<array<string>> element_names(provable_elements.size + 1);
				bool has_unnamed_elements = (provable_elements.size < T.sets.sets[set_id].set_size);
				named_constant_or_set = true;
				for (const tuple& tup : provable_elements) {
					array<string>& current_element_names = element_names[element_names.length];
					if (!array_init(current_element_names, max(1, tup.length))) {
						for (array<string>& name_array : element_names) {
							for (string& str : name_array) free(str);
							free(name_array);
						} for (tuple& tup : provable_elements) free(tup);
						return;
					}
					element_names.length++;
					for (unsigned int i = 0; i < tup.length; i++) {
						int length;
						string& next_name = current_element_names[current_element_names.length];
						switch (tup[i].type) {
						case tuple_element_type::NUMBER:
							if (tup[i].number.decimal == 0)
								length = snprintf(NULL, 0, "%" PRId64, tup[i].number.integer);
							else length = snprintf(NULL, 0, "%" PRId64 ".%" PRIu64, tup[i].number.integer, tup[i].number.decimal);
							if (!init(next_name, length)) {
								for (array<string>& name_array : element_names) {
									for (string& str : name_array) free(str);
									free(name_array);
								} for (tuple& tup : provable_elements) free(tup);
								return;
							}
							current_element_names.length++;
							if (tup[i].number.decimal == 0)
								snprintf(next_name.data, length + 1, "%" PRId64, tup[i].number.integer);
							else snprintf(next_name.data, length + 1, "%" PRId64 ".%" PRIu64, tup[i].number.integer, tup[i].number.decimal);
							next_name.length = length;
							break;
						case tuple_element_type::STRING:
							if (!init(next_name, tup[i].str)) {
								for (array<string>& name_array : element_names) {
									for (string& str : name_array) free(str);
									free(name_array);
								} for (tuple& tup : provable_elements) free(tup);
								return;
							}
							current_element_names.length++;
							break;
						case tuple_element_type::CONSTANT:
							if (tup[i].constant < T.new_constant_offset) {
								const string_map_scribe& printer = parser.get_printer();
								if (tup[i].constant < printer.length) {
									if (!init(next_name, *printer.map[tup[i].constant])) {
										for (array<string>& name_array : element_names) {
											for (string& str : name_array) free(str);
											free(name_array);
										} for (tuple& tup : provable_elements) free(tup);
										return;
									}
									current_element_names.length++;
								} else if (tup.length == 1) {
									if (!init(next_name, UNKNOWN_CONCEPT_NAME)) {
										for (array<string>& name_array : element_names) {
											for (string& str : name_array) free(str);
											free(name_array);
										} for (tuple& tup : provable_elements) free(tup);
										return;
									}
									current_element_names.length++;
								} else {
									has_unnamed_elements = true;
								}
							} else {
								array<Term*> name_terms(2);
								if (!T.get_concept_names(tup[i].constant, name_terms))
									return;
								if (name_terms.length == 0) {
									if (tup.length == 1) {
										if (!init(next_name, UNKNOWN_CONCEPT_NAME)) {
											for (array<string>& name_array : element_names) {
												for (string& str : name_array) free(str);
												free(name_array);
											} for (tuple& tup : provable_elements) free(tup);
											return;
										}
										current_element_names.length++;
									} else {
										has_unnamed_elements = true;
									}
								} else {
									insertion_sort(name_terms, pointer_sorter());
									unsigned int name_length = (name_terms.length - 1);
									for (Term* name_term : name_terms)
										name_length += name_term->str.length;
									if (!init(next_name, name_length)) {
										for (array<string>& name_array : element_names) {
											for (string& str : name_array) free(str);
											free(name_array);
										} for (tuple& tup : provable_elements) free(tup);
										return;
									}
									current_element_names.length++;
									name_length = 0;
									for (unsigned int i = 0; i < name_terms.length; i++) {
										if (i != 0) next_name[name_length++] = '/';
										for (unsigned int j = 0; j < name_terms[i]->str.length; j++)
											next_name[name_length++] = name_terms[i]->str[j];
									}
								}
							}
							break;
						}
					}
					if (current_element_names.length == 0) {
						free(current_element_names);
						element_names.length--;
					}
				}
				for (tuple& tup : provable_elements) free(tup);

				string& new_name = *((string*) alloca(sizeof(string)));
				if (element_names.length == 0) {
					if (!init(new_name, (has_unnamed_elements ? "{...}" : "{}"))) {
						for (array<string>& name_array : element_names) {
							for (string& str : name_array) free(str);
							free(name_array);
						}
						return;
					}
				} else {
					for (array<string>& name_array : element_names) {
						if (name_array.length > 1)
							insertion_sort(name_array);
					}
					insertion_sort(element_names, string_array_sorter());

					unsigned int string_length = 2 + (element_names.length - 1);
					for (const array<string>& name_array : element_names) {
						if (name_array.length > 1)
							string_length += 2 + name_array.length - 1;
						for (const string& str : name_array)
							string_length += str.length;
					}
					if (has_unnamed_elements)
						string_length += 4;

					if (!init(new_name, string_length)) {
						for (array<string>& name_array : element_names) {
							for (string& str : name_array) free(str);
							free(name_array);
						}
						return;
					}
					string_length = 0;
					new_name[string_length++] = '{';
					for (unsigned int i = 0; i < element_names.length; i++) {
						if (i != 0) new_name[string_length++] = ',';
						const array<string>& name_array = element_names[i];
						if (name_array.length > 1)
							new_name[string_length++] = '(';
						for (unsigned int j = 0; j < name_array.length; j++) {
							if (j != 0) new_name[string_length++] = ',';
							const string& src = name_array[j];
							for (unsigned int k = 0; k < src.length; k++)
								new_name[string_length++] = src[k];
						}
						if (name_array.length > 1)
							new_name[string_length++] = ')';
					}
					if (has_unnamed_elements) {
						new_name[string_length++] = ',';
						new_name[string_length++] = '.';
						new_name[string_length++] = '.';
						new_name[string_length++] = '.';
					}
					new_name[string_length++] = '}';
					for (array<string>& name_array : element_names) {
						for (string& str : name_array) free(str);
						free(name_array);
					}
				}

				if (!temp_answers.ensure_capacity(temp_answers.size + 1))
					return;
				unsigned int index = temp_answers.index_of(new_name);
				if (index < temp_answers.size) {
					free(new_name);
					temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
				} else {
					move(new_name, temp_answers.keys[index]);
					temp_answers.values[index] = log_probability;
					temp_answers.size++;
				}
			}

			if (!named_constant_or_set) {
print(term->constant, stderr); print(": <unnamed>, log probability: ", stderr);
print(log_probability, stderr); print('\n', stderr);
T.print_axioms(stderr); print('\n', stderr);
				if (!temp_answers.ensure_capacity(temp_answers.size + 1))
					return;
				unsigned int index = temp_answers.index_of(UNKNOWN_CONCEPT_NAME);
				if (index < temp_answers.size) {
					temp_answers.values[index] = logsumexp(temp_answers.values[index], log_probability);
				} else {
					if (!init(temp_answers.keys[index], UNKNOWN_CONCEPT_NAME))
						return;
					temp_answers.values[index] = log_probability;
					temp_answers.size++;
				}
			}
		} else {
			fprintf(stderr, "ERROR: Unable to convert semantic answer into text.\n");
		}
print("Totals so far:\n", stderr);
for (const auto& entry : temp_answers) {
print('"', stderr); print(entry.key, stderr); print("\": ", stderr);
print(entry.value, stderr); print('\n', stderr);
}
	};

/* TODO: for debugging; delete this */
extern const string_map_scribe* debug_terminal_printer;
debug_terminal_printer = &parser.get_printer();
	theory<ProofCalculus, Canonicalizer>& T_map = *((theory<ProofCalculus, Canonicalizer>*) alloca(sizeof(theory<ProofCalculus, Canonicalizer>)));
	if (!log_joint_probability_of_lambda(T, theory_prior, proof_axioms, logical_forms[0], num_samples, T_map, on_new_proof_sample, std::forward<Args>(add_formula_args)...)) {
		fprintf(stderr, "ERROR: Failed to answer question.\n");
		free_logical_forms(logical_forms, parse_count);
		for (auto entry : temp_answers) free(entry.key);
		return false;
	}
T_map.print_axioms(stderr, *debug_terminal_printer);
	free(T_map);

	/* check if we are confident enough in the most probable answer to stop looking for more information */
	if (temp_answers.size > 1)
		sort(temp_answers.values, temp_answers.keys, temp_answers.size, default_sorter());
	if ((temp_answers.size > 1 && temp_answers.values[temp_answers.size - 1] - temp_answers.values[temp_answers.size - 2] < SUFFICIENT_KNOWLEDGE_THRESHOLD)
	 || (temp_answers.size != 0 && temp_answers.keys[temp_answers.size - 1] == UNKNOWN_CONCEPT_NAME))
	{
		for (auto pair : temp_answers) free(pair.key);
		temp_answers.clear();

		/* TODO: iterate over possible explanations of the logical form */

		/* generate a search query */
		unsigned int generated_derivation_count;
		double generated_log_likelihood;
		syntax_node<typename Parser::logical_form_type>& generated_derivation =
				*((syntax_node<typename Parser::logical_form_type>*) alloca(sizeof(syntax_node<typename Parser::logical_form_type>)));
		/* TODO: disable the production rules that govern wh-movement */
		if (!parser.template generate<1>(&generated_derivation, &generated_log_likelihood, generated_derivation_count, logical_forms[0], names) || generated_derivation_count == 0) {
			fprintf(stderr, "ERROR: Failed to generate search query derivation.\n");
			free_logical_forms(logical_forms, parse_count);
			/* TODO: re-enable the production rules that govern wh-movement */
			return false;
		}
		/* TODO: re-enable the production rules that govern wh-movement */

		/* compute the yield of the derivation */
		sequence search_query = sequence(NULL, 0);
		if (!parser.yield_search_query(generated_derivation, logical_forms[0], search_query)) {
			free_logical_forms(logical_forms, parse_count);
			free(generated_derivation); return false;
		}
		free(generated_derivation);

		/* remove the punctuation at the end of the sentence; if the asterisk
		   appears at the beginning or the end, delete it */
		search_query.length--;

		/* find the left and right portions of the query */
		unsigned int index = index_of(parser.ASTERISK_ID, search_query.tokens, search_query.length);
		sequence left_query(nullptr, 0), right_query(nullptr, 0);
		if (index == 0) {
			shift_left(search_query.tokens, --search_query.length);
			right_query.tokens = search_query.tokens;
			right_query.length = search_query.length;
		} else if (index == search_query.length - 1) {
			search_query.length--;
			left_query.tokens = search_query.tokens;
			left_query.length = search_query.length;
		} else {
			left_query.tokens = search_query.tokens;
			left_query.length = index;
			right_query.tokens = search_query.tokens + index + 1;
			right_query.length = search_query.length - index - 1;
		}

		bool failure = false;
		auto process_matched_sentence = [&failure, &articles, &parser, &T, &names, &visited_articles, &theory_prior, &proof_axioms, &logical_forms, num_samples, on_new_proof_sample, &temp_answers](const html_lexer_token* tokens, unsigned int length) {
			if (read_sentence(articles, parser, tokens, length, T, names, visited_articles, theory_prior, proof_axioms)) {
				theory<ProofCalculus, Canonicalizer>& T_map = *((theory<ProofCalculus, Canonicalizer>*) alloca(sizeof(theory<ProofCalculus, Canonicalizer>)));
				if (!log_joint_probability_of_lambda(T, theory_prior, proof_axioms, logical_forms[0], num_samples, T_map, on_new_proof_sample)) {
					fprintf(stderr, "ERROR: Failed to answer question.\n");
					failure = true; return false;
				}
T_map.print_axioms(stderr, *debug_terminal_printer);
				free(T_map);

				if (temp_answers.size > 1) {
					sort(temp_answers.values, temp_answers.keys, temp_answers.size, default_sorter());
					if (temp_answers.values[temp_answers.size - 1] - temp_answers.values[temp_answers.size - 2] >= SUFFICIENT_KNOWLEDGE_THRESHOLD)
						return false; /* break the Google search */
				}
				return true;
			}
			return true;
		};
		constexpr unsigned long long TIMEOUT_MS = 5000;
		auto process_search_result = [&left_query,&right_query,&names,process_matched_sentence](const string& result) {
			return find_answer_in_website(result, TIMEOUT_MS, left_query, right_query, names, process_matched_sentence);
		};

		memory_stream out(32);
		string* query = (string*) malloc(sizeof(string) * search_query.length);
		for (unsigned int i = 0; i < search_query.length; i++) {
			if (!print(search_query[i], out, parser.get_printer())
			 || !init(query[i], out.buffer, out.position))
			{
				free_logical_forms(logical_forms, parse_count);
				for (unsigned int j = 0; j < i; j++) free(query[j]);
				free(query); free(search_query); return false;
			}
			out.position = 0;
			out.shift = {0};
		}

		if (!search_google(query, search_query.length, TIMEOUT_MS, process_search_result)) {
			free_logical_forms(logical_forms, parse_count);
			for (unsigned int j = 0; j < search_query.length; j++) free(query[j]);
			free(query); free(search_query); return false;
		}
		free_logical_forms(logical_forms, parse_count);
		for (unsigned int j = 0; j < search_query.length; j++) free(query[j]);
		free(query); free(search_query);

		if (failure) return false;

	} else {
		free_logical_forms(logical_forms, parse_count);
	}

	/* keep only the answers with highest probability */
	if (temp_answers.size == 0)
		return true;
	for (unsigned int i = temp_answers.size; i > 0; i--) {
		if (!answers.ensure_capacity(answers.length + 1)
		 || !init(answers[answers.length], temp_answers.keys[i - 1]))
		{
			for (auto pair : temp_answers) free(pair.key);
			for (string& str : answers) free(str);
			return false;
		}
		answers.length++;
		if (i > 1 && temp_answers.values[i - 2] < temp_answers.values[i - 1])
			break;
	}
	for (auto pair : temp_answers) free(pair.key);
	return true;
}

#endif /* EXECUTIVE_H_ */
