#ifndef ARTICLE_H_
#define ARTICLE_H_

#include <core/lex.h>
#include <cstdint>

using namespace core;

struct token {
	unsigned int id;
};

struct sentence {
	token* tokens;
	unsigned int length;

	static inline void free() {
		free(tokens);
	}
};

struct article {
	sentence* sentences;
	unsigned int sentence_count;

	~article() { free(); }

	static inline void free(article& a) { a.free(); }

private:
	inline void free() {
		for (unsigned int i = 0; i < sentence_count; i++)
			core::free(sentences[i]);
		core::free(sentences);
	}
};

struct in_memory_article_store {
	hash_map<unsigned int, article>& articles;

	~in_memory_article_store() {
		for (auto entry : articles)
			free(entry.value);
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

inline bool append_to_token(
	array<char>& token, wint_t next, std::mbstate_t& shift)
{
	if (!token.ensure_capacity(token.length + MB_CUR_MAX))
		return false;
	size_t written = wcrtomb(token.data + token.length, next, &shift);
	if (written == static_cast<std::size_t>(-1))
		return false;
	token.length += written;
	return true;
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
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
				new_line = (next == '\n');
				if (new_line && !emit_token(tokens, start, start + 1, article_token_type::NEWLINE))
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
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
				if (new_line && !emit_token(tokens, start, start + 1, article_token_type::NEWLINE))
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

	if (!init(out, length)) return false;

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
		hash_map<string, unsigned int>& names)
{
	array<token>& sentence_tokens = *((array<token>*) alloca(sizeof(array<tokens>)));
	if (!array_init(sentence_tokens, 16)) return false;
	while (true) {
		unsigned int token_id;
		if (!expect_token(tokens, index, article_token_type::TOKEN, "sentence token")
		 || !get_token(tokens[index], token_id, names) || !sentence_tokens.add({token_id}))
		{
			free(sentence_tokens); return false;
		}
		index++;

		if (index >= tokens.length) {
			read_error("Expected a token or period at end of sentence", tokens[index].start);
			free(sentence_tokens); return false;
		} else if (tokens[index].type == article_token_type::PERIOD) {
			index++; break;
		} else if (tokens[index].type == article_token_type::TOKEN) {
			read_error("Expected a token or period at end of sentence", tokens[index].start);
			free(sentence_tokens); return false;
		}
	}

	move(sentence_tokens.data, out.tokens);
	move(sentence_tokens.length, out.length);
	return true;
}

bool article_interpret(
		const array<article_token>& tokens,
		unsigned int& index,
		article& out, unsigned int& article_name,
		hash_map<sentence, fol_formula*> logical_forms,
		hash_map<string, unsigned int>& names)
{
	unsigned int start = index;
	if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
		return false;
	index++;
	while (tokens[index].type != article_token_type::COLON) {
		if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
			return false;
		index+;
	}

	string name;
	if (!concatenate(tokens.data + start, index - start, name)) {
		fprintf(stderr, "article_interpret ERROR: Unable to construct string for article name.\n");
		return false;
	} else if (!get_token(name, article_name, names)) {
		free(name);
		return false;
	}
	free(name);

	array<sentence>& sentences = *((array<sentence>*) alloca(sizeof(array<sentence>)));
	if (!array_init(sentences, 8))
		return false;
	while (true) {
		if (!sentences.ensure_capacity(sentences.length + 1)
		 || !article_interpret_sentence(tokens, index, sentences[sentences.length], names))
		{
			for (sentence& s : sentences) free(s);
			free(sentences); return false;
		}
		sentences.length++;

		if (index >= tokens.length) {
			break;
		} else if (tokens[index].type == article_token_type::NEWLINE) {
			index++; break;
		} else if (tokens[index].type == article_token_type::FORMULA) {
			/* parse the formula */
			array<tptp_token> formula_tokens = array<tptp_token>(128);
			memory_stream& in = *((memory_stream*) alloca(sizeof(memory_stream)));
			in.buffer = tokens[index].text.data;
			in.length = tokens[index].text.length;
			in.position = 0; in.shift = {0};
			if (!tptp_lex(formula_tokens, in)
			 || !logical_forms.check_size())
			{
				read_error("Unable to parse first-order formula (lexical analysis failed)", tokens[index].start);
				for (sentence& s : sentences) free(s);
				free(sentences); free_tokens(formula_tokens);
				return false;
			}

			unsigned int formula_index = 0;
			array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
			fol_formula* formula = (fol_formula*) malloc(sizeof(fol_formula));
			if (formula == NULL) {
				fprintf(stderr, "article_interpret ERROR: Insufficient memory for fol_formula.\n");
				for (sentence& s : sentences) free(s);
				free(sentences); free_tokens(formula_tokens);
				return false;
			} else if (!tptp_interpret(formula_tokens, formula_index, *formula, names, variables)) {
				read_error("Unable to parse first-order formula", tokens[index].start);
				for (auto entry : variables) free(entry.key);
				for (sentence& s : sentences) free(s);
				free(sentences); free(formula); free_tokens(formula_tokens);
				return false;
			}
			free_tokens(formula_tokens);

			bool contains; unsigned int bucket;
			fol_formula*& value = logical_forms.get(sentences.last(), contains, bucket);
			if (!contains) {
				value = formula;
				logical_forms.table.keys[bucket] = sentences.last();
				logical_forms.table.size++;
			} else if (*value != *formula) {
				read_error("A different logical form was previously mapped to this sentence", tokens[index].start);
				for (sentence& s : sentences) free(s);
				free(sentences); free(*formula); free(formula);
				free_tokens(formula_tokens);
				return false;
			}
			index++;
		}
	}

	move(sentences.data, article.sentences);
	move(sentences.length, article.sentence_count);
	return true;
}

bool articles_interpret(
		const array<article_token>& tokens,
		in_memory_article_store& articles,
		hash_map<sentence, fol_formula*> logical_forms,
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
