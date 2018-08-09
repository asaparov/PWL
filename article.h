#ifndef ARTICLE_H_
#define ARTICLE_H_

#include <core/lex.h>

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

	~article() {
		for (unsigned int i = 0; i < sentence_count; i++)
			free(sentences[i]);
		free(sentences);
	}
};


/**
 * Code for tokenizing/lexing the articles.
 */

enum class article_token_type {
	LBRACKET,
	RBRACKET,
	PERIOD,
	COLON,
	IDENTIFIER
};

typedef lexical_token<article_token_type> article_token;

template<typename Stream>
inline bool print(article_token_type type, Stream& stream) {
	switch (type) {
	case article_token_type::LBRACKET:
		return print('[', stream);
	case article_token_type::RBRACKET:
		return print(']', stream);
	case article_token_type::PERIOD:
		return print('.', stream);
	case article_token_type::COLON:
		return print(':', stream);
	case article_token_type::IDENTIFIER:
		return print("IDENTIFIER", stream);
	}
	fprintf(stderr, "print ERROR: Unknown article_token_type.\n");
	return false;
}

enum class article_lexer_state {
	DEFAULT,
	IDENTIFIER,
};

bool article_emit_symbol(array<article_token>& tokens, const position& start, char symbol) {
	switch (symbol) {
	case '[':
		return emit_token(tokens, start, start + 1, article_token_type::LBRACKET);
	case ']':
		return emit_token(tokens, start, start + 1, article_token_type::RBRACKET);
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
		case article_lexer_state::IDENTIFIER:
			if (next == '[' || next == ']' || next == '.' || next == ':')
			{
				if (!emit_token(tokens, token, start, current, article_token_type::IDENTIFIER)
				 || !article_emit_symbol(tokens, current, next))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, article_token_type::IDENTIFIER))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
				new_line = (next == '\n');
			} else {
				if (!append_to_token(token, next, shift)) return false;
			}
			break;

		case article_lexer_state::DEFAULT:
			if (next == '[' || next == ']' || next == '.' || next == ':')
			{
				if (!article_emit_symbol(tokens, current, next))
					return false;
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
			} else {
				if (!append_to_token(token, next, shift)) return false;
				state = article_lexer_state::IDENTIFIER;
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

	if (state == article_lexer_state::IDENTIFIER)
		return emit_token(tokens, token, start, current, article_token_type::IDENTIFIER);
	return true;
}


/**
 * Recursive-descent parser for articles.
 */

bool article_interpret(
	const array<article_token>& tokens, article& out,
	hash_map<sentence, fol_formula*> logical_forms)
{
	unsigned int index = 0;
	while (true) {
		
	}
}

#endif /* ARTICLE_H_ */
