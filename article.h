#ifndef ARTICLE_H_
#define ARTICLE_H_

#include <core/lex.h>
#include <cstdint>

#include "first_order_logic.h"

using namespace core;

struct token {
	unsigned int id;

	static inline void move(const token& src, token& dst) {
		dst.id = src.id;
	}

	static inline unsigned int hash(const token& key) {
		return default_hash(key.id);
	}

	static inline bool is_empty(const token& key) {
		return key.id == 0;
	}
};

inline bool operator == (const token& first, const token& second) {
	return first.id == second.id;
}

inline bool operator != (const token& first, const token& second) {
	return first.id != second.id;
}

struct sentence {
	token* tokens;
	unsigned int length;

	static inline bool is_empty(const sentence& key) {
		return key.tokens == NULL;
	}

	static inline unsigned int hash(const sentence& key) {
		return default_hash(key.tokens, key.length);
	}

	static inline void move(const sentence& src, sentence& dst) {
		dst.tokens = src.tokens;
		dst.length = src.length;
	}

	static inline void free(sentence& s) {
		core::free(s.tokens);
	}
};

inline bool init(sentence& dst, const sentence& src) {
	dst.tokens = (token*) malloc(max((size_t) 1, sizeof(token) * src.length));
	if (dst.tokens == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for sentence.tokens.\n");
		return false;
	}
	for (unsigned int i = 0; i < src.length; i++)
		dst.tokens[i] = src.tokens[i];
	dst.length = src.length;
	return true;
}

inline bool operator == (const sentence& first, const sentence& second) {
	if (first.length != second.length) return false;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.tokens[i] != second.tokens[i]) return false;
	return true;
}

struct article {
	sentence* sentences;
	unsigned int sentence_count;

	~article() { free(); }

	static inline void move(const article& src, article& dst) {
		dst.sentences = src.sentences;
		dst.sentence_count = src.sentence_count;
	}

	static inline void free(article& a) { a.free(); }

private:
	inline void free() {
		for (unsigned int i = 0; i < sentence_count; i++)
			core::free(sentences[i]);
		core::free(sentences);
	}
};

struct in_memory_article_store {
	hash_map<unsigned int, article> articles;

	in_memory_article_store() : articles(64) { }

	~in_memory_article_store() {
		for (auto entry : articles)
			free(entry.value);
	}

	const article& get(unsigned int article_id, bool& contains) const {
		return articles.get(article_id, contains);
	}
};


/**
 * Code for tokenizing/lexing the articles.
 */

enum class article_token_type {
	PERIOD,
	COLON,
	NEWLINE,
	TOKEN,
	FORMULA
};

typedef lexical_token<article_token_type> article_token;

template<typename Stream>
inline bool print(article_token_type type, Stream& stream) {
	switch (type) {
	case article_token_type::PERIOD:
		return print('.', stream);
	case article_token_type::COLON:
		return print(':', stream);
	case article_token_type::NEWLINE:
		return print("NEWLINE", stream);
	case article_token_type::TOKEN:
		return print("TOKEN", stream);
	case article_token_type::FORMULA:
		return print("FORMULA", stream);
	}
	fprintf(stderr, "print ERROR: Unknown article_token_type.\n");
	return false;
}

enum class article_lexer_state {
	DEFAULT,
	TOKEN,
	FORMULA
};

bool article_emit_symbol(array<article_token>& tokens, const position& start, char symbol) {
	switch (symbol) {
	case '.':
		return emit_token(tokens, start, start + 1, article_token_type::PERIOD);
	case ':':
		return emit_token(tokens, start, start + 1, article_token_type::COLON);
	default:
		fprintf(stderr, "article_emit_symbol ERROR: Unexpected symbol.\n");
		return false;
	}
}

template<typename Stream>
bool article_lex(array<article_token>& tokens, Stream& input) {
	position start = position(1, 1);
	position current = position(1, 1);
	article_lexer_state state = article_lexer_state::DEFAULT;
	array<char> token = array<char>(1024);

	std::mbstate_t shift = {0};
	wint_t next = fgetwc(input);
	bool new_line = false;
	while (next != WEOF) {
		switch (state) {
		case article_lexer_state::TOKEN:
			if (next == '.' || next == ':') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN)
				 || !article_emit_symbol(tokens, current, next))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
			} else if (next == '{') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN))
					return false;
				state = article_lexer_state::FORMULA;
				token.clear(); shift = {0};
				start = current;
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
				new_line = (next == '\n');
				if (new_line && !emit_token(tokens, current, current + 1, article_token_type::NEWLINE))
					return false;
			} else {
				if (!append_to_token(token, next, shift)) return false;
			}
			break;

		case article_lexer_state::FORMULA:
			if (next == '}') {
				if (!emit_token(tokens, token, start, current, article_token_type::FORMULA))
					return false;
				state = article_lexer_state::DEFAULT;
			} else {
				if (!append_to_token(token, next, shift)) return false;
				new_line = (next == '\n');
			}
			break;

		case article_lexer_state::DEFAULT:
			if (next == '.' || next == ':') {
				if (!article_emit_symbol(tokens, current, next))
					return false;
			} else if (next == '{') {
				state = article_lexer_state::FORMULA;
				token.clear(); shift = {0};
				start = current;
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
				if (new_line && !emit_token(tokens, current, current + 1, article_token_type::NEWLINE))
					return false;
			} else {
				if (!append_to_token(token, next, shift)) return false;
				state = article_lexer_state::TOKEN;
				start = current;
			}
			break;
		}

		if (new_line) {
			current.line++;
			current.column = 1;
			new_line = false;
		} else current.column++;
		next = fgetwc(input);
	}

	if (state == article_lexer_state::TOKEN)
		return emit_token(tokens, token, start, current, article_token_type::TOKEN);
	return true;
}


/**
 * Recursive-descent parser for articles.
 */

template<typename TokenType>
inline bool concatenate(lexical_token<TokenType>* tokens, unsigned int token_count, string& out)
{
	unsigned int length = tokens[0].text.length;
	for (unsigned int i = 1; i < token_count; i++)
		length += 1 + tokens[i].text.length;

	if (!core::resize(out.data, length)) return false;
	out.length = length;

	unsigned int index = 0;
	for (unsigned int i = 0; i < token_count; i++) {
		if (i > 0) { out[index] = ' '; index++; }
		for (unsigned int j = 0; j < tokens[i].text.length; j++) {
			out[index] = tokens[i].text[j];
			index++;
		}
	}
	return true;
}

bool article_interpret_sentence(
		const array<article_token>& tokens,
		unsigned int& index, sentence& out,
		array_map<unsigned int, unsigned int>& labels,
		hash_map<string, unsigned int>& names)
{
	array<token>& sentence_tokens = *((array<token>*) alloca(sizeof(array<token>)));
	if (!array_init(sentence_tokens, 16)) return false;
	while (true) {
		unsigned int token_id;
		if (!expect_token(tokens, index, article_token_type::TOKEN, "sentence token")
		 || !get_token(tokens[index].text, token_id, names) || !sentence_tokens.add({token_id}))
		{
			free(sentence_tokens); return false;
		}
		index++;

		if (index < tokens.length && tokens[index].type == article_token_type::COLON) {
			index++;
			unsigned int token_label;
			if (!expect_token(tokens, index, article_token_type::TOKEN, "token label")
			 || !get_token(tokens[index].text, token_label, names) || !labels.put(token_id, token_label))
			{
				free(sentence_tokens); return false;
			}
			index++;
		}

		if (index >= tokens.length) {
			read_error("Expected a token or period at end of sentence", tokens[index].start);
			free(sentence_tokens); return false;
		} else if (tokens[index].type == article_token_type::PERIOD) {
			index++; break;
		} else if (tokens[index].type != article_token_type::TOKEN) {
			read_error("Expected a token or period at end of sentence", tokens[index].start);
			free(sentence_tokens); return false;
		}
	}

	move(sentence_tokens.data, out.tokens);
	out.length = sentence_tokens.length;
	return true;
}

struct sentence_label {
	fol_formula* logical_form;
	array_map<unsigned int, unsigned int> labels;

	static inline void move(const sentence_label& src, sentence_label& dst) {
		dst.logical_form = src.logical_form;
		core::move(src.labels, dst.labels);
	}

	static inline void free(sentence_label& label) {
		if (label.logical_form->reference_count > 0) {
			core::free(*label.logical_form);
			if (label.logical_form->reference_count == 0)
				core::free(label.logical_form);
		} else {
			core::free(label.logical_form);
		}
		core::free(label.labels);
	}
};

inline bool init(sentence_label& label) {
	label.logical_form = (fol_formula*) malloc(sizeof(fol_formula));
	if (label.logical_form == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for sentence_label.logical_form.\n");
		return false;
	} else if (!array_map_init(label.labels, 8)) {
		free(label.logical_form);
		return false;
	}
	label.logical_form->reference_count = 0;
	return true;
}

inline bool operator != (const sentence_label& first, const sentence_label& second) {
	if (*first.logical_form != *second.logical_form)
		return true;
	if (first.labels.size != second.labels.size)
		return true;
	for (unsigned int i = 0; i < first.labels.size; i++)
		if (first.labels.keys[i] != second.labels.keys[i]
		 || first.labels.values[i] != second.labels.values[i]) return true;
	return false;
}

bool article_interpret(
		const array<article_token>& tokens,
		unsigned int& index,
		article& out, unsigned int& article_name,
		hash_map<sentence, sentence_label>& logical_forms,
		hash_map<string, unsigned int>& names)
{
	unsigned int start = index;
	if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
		return false;
	index++;
	while (tokens[index].type != article_token_type::COLON) {
		if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
			return false;
		index++;
	}

	string name = string(16);
	if (!concatenate(tokens.data + start, index - start, name)) {
		fprintf(stderr, "article_interpret ERROR: Unable to construct string for article name.\n");
		return false;
	} else if (!get_token(name, article_name, names)) {
		return false;
	}
	index++;

	array<sentence>& sentences = *((array<sentence>*) alloca(sizeof(array<sentence>)));
	if (!array_init(sentences, 8))
		return false;
	while (true) {
		sentence_label& new_label = *((sentence_label*) alloca(sizeof(sentence_label)));
		if (!init(new_label)) {
			for (sentence& s : sentences) free(s);
			free(sentences); return false;
		} else if (!sentences.ensure_capacity(sentences.length + 1)
				|| !article_interpret_sentence(tokens, index, sentences[sentences.length], new_label.labels, names))
		{
			for (sentence& s : sentences) free(s);
			free(sentences); free(new_label); return false;
		}
		sentences.length++;

		if (index < tokens.length && tokens[index].type == article_token_type::FORMULA) {
			/* parse the formula */
			array<tptp_token> formula_tokens = array<tptp_token>(128);
			memory_stream& in = *((memory_stream*) alloca(sizeof(memory_stream)));
			in.buffer = tokens[index].text.data;
			in.length = tokens[index].text.length;
			in.position = 0; in.shift = {0};
			if (!tptp_lex(formula_tokens, in, tokens[index].start + 1)
			 || !logical_forms.check_size())
			{
				read_error("Unable to parse first-order formula (lexical analysis failed)", tokens[index].start);
				for (sentence& s : sentences) free(s);
				free(sentences); free(new_label); free_tokens(formula_tokens);
				return false;
			}

			unsigned int formula_index = 0;
			array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
			if (!tptp_interpret(formula_tokens, formula_index, *new_label.logical_form, names, variables)) {
				read_error("Unable to parse first-order formula", tokens[index].start);
				for (auto entry : variables) free(entry.key);
				for (sentence& s : sentences) free(s);
				free(sentences); free(new_label); free_tokens(formula_tokens);
				return false;
			}
			free_tokens(formula_tokens);

			bool contains; unsigned int bucket;
			sentence_label& value = logical_forms.get(sentences.last(), contains, bucket);
			if (!contains) {
				if (!init(logical_forms.table.keys[bucket], sentences.last())) {
					for (sentence& s : sentences) free(s);
					free(sentences); free(new_label);
					free_tokens(formula_tokens);
					return false;
				}
				move(new_label, value);
				logical_forms.table.size++;
			} else if (value != new_label) {
				read_error("A different logical form was previously mapped to this sentence", tokens[index].start);
				for (sentence& s : sentences) free(s);
				free(sentences); free(new_label);
				free_tokens(formula_tokens);
				return false;
			}
			index++;
		}

		if (index >= tokens.length) {
			break;
		} else if (tokens[index].type == article_token_type::NEWLINE) {
			index++; break;
		}
	}

	move(sentences.data, out.sentences);
	out.sentence_count = sentences.length;
	return true;
}

bool articles_interpret(
		const array<article_token>& tokens,
		in_memory_article_store& articles,
		hash_map<sentence, sentence_label>& logical_forms,
		hash_map<string, unsigned int>& names)
{
	unsigned int index = 0;
	while (index < tokens.length) {
		if (!articles.articles.check_size())
			return false;

		position article_pos = tokens[index].start;

		unsigned int article_name;
		article& new_article = *((article*) alloca(sizeof(article)));
		if (!article_interpret(tokens, index, new_article, article_name, logical_forms, names))
			return false;

		bool contains; unsigned int bucket;
		article& value = articles.articles.get(article_name, contains, bucket);
		if (!contains) {
			move(new_article, value);
			articles.articles.table.keys[bucket] = article_name;
			articles.articles.table.size++;
		} else {
			read_error("Article name already exists", article_pos);
			free(value); return false;
		}
	}
	return true;
}

#endif /* ARTICLE_H_ */
