#ifndef ARTICLE_H_
#define ARTICLE_H_

#include <core/lex.h>
#include <cstdint>

using namespace core;

struct sentence_token {
	unsigned int id;

	static inline void move(const sentence_token& src, sentence_token& dst) {
		dst.id = src.id;
	}

	static inline unsigned int hash(const sentence_token& key) {
		return default_hash(key.id);
	}

	static inline bool is_empty(const sentence_token& key) {
		return key.id == 0;
	}
};

inline bool operator == (const sentence_token& first, const sentence_token& second) {
	return first.id == second.id;
}

inline bool operator != (const sentence_token& first, const sentence_token& second) {
	return first.id != second.id;
}

template<typename Stream, typename... Printer>
inline bool print(const sentence_token& t, Stream& out, Printer&&... printer) {
	return print(t.id, out, std::forward<Printer>(printer)...);
}

template<typename Derivation>
struct sentence {
	sentence_token* tokens;
	unsigned int length;
	Derivation derivation;

	static inline bool is_empty(const sentence& key) {
		return key.tokens == nullptr;
	}

	static inline unsigned int hash(const sentence& key) {
		if (!core::is_empty(key.derivation))
			return hasher<Derivation>::hash(key.derivation);
		return default_hash(key.tokens, key.length);
	}

	static inline void move(const sentence& src, sentence& dst) {
		dst.tokens = src.tokens;
		dst.length = src.length;
		dst.derivation = src.derivation;
	}

	static inline void free(sentence& s) {
		if (s.tokens != nullptr)
			core::free(s.tokens);
		if (!core::is_empty(s.derivation))
			core::free(s.derivation);
	}
};

template<typename Derivation>
inline bool init(sentence<Derivation>& dst, const sentence<Derivation>& src) {
	if (src.tokens == nullptr) {
		dst.tokens = nullptr;
	} else {
		dst.tokens = (sentence_token*) malloc(sizeof(sentence_token) * max((size_t) 1, src.length));
		if (dst.tokens == NULL) {
			fprintf(stderr, "init ERROR: Insufficient memory for sentence.tokens.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.length; i++)
			dst.tokens[i] = src.tokens[i];
		dst.length = src.length;
	}

	dst.derivation = src.derivation;
	return true;
}

template<typename Derivation>
inline bool operator == (const sentence<Derivation>& first, const sentence<Derivation>& second) {
	if (first.tokens == nullptr) {
		if (second.tokens != nullptr) return false;
	} else {
		if (second.tokens == nullptr) {
			return false;
		} else {
			if (first.length != second.length) return false;
			for (unsigned int i = 0; i < first.length; i++)
				if (first.tokens[i] != second.tokens[i]) return false;
		}
	}

	if (core::is_empty(first.derivation)) {
		if (core::is_empty(second.derivation)) {
			return true;
		} else {
			return false;
		}
	} else {
		if (core::is_empty(second.derivation)) {
			return false;
		} else {
			return first.derivation == second.derivation;
		}
	}
}

template<typename Derivation, typename Stream, typename... Printer>
bool print(const sentence<Derivation>& s, Stream&& out, Printer&&... printer) {
	if (!core::is_empty(s.derivation)) {
		return print("<derivation tree printing not implemented>", out);
		//return print(s.derivation, out, std::forward<Printer>(printer)...);
	}

	if (s.length == 0) return true;
	if (!print(s.tokens[0], out, std::forward<Printer>(printer)...)) return false;
	for (unsigned int i = 1; i < s.length; i++) {
		if (!print(' ', out) || !print(s.tokens[i], out, std::forward<Printer>(printer)...)) return false;
	}
	return true;
}

template<typename Derivation>
struct article {
	sentence<Derivation>* sentences;
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

template<typename Derivation>
struct in_memory_article_store {
	hash_map<unsigned int, article<Derivation>> articles;

	in_memory_article_store() : articles(64) { }

	~in_memory_article_store() {
		for (auto entry : articles)
			free(entry.value);
	}

	const article<Derivation>& get(unsigned int article_id, bool& contains) const {
		return articles.get(article_id, contains);
	}

	inline bool contains(unsigned int article_id) const {
		return articles.table.contains(article_id);
	}
};


/**
 * Code for tokenizing/lexing the articles.
 */

enum class article_token_type {
	PERIOD,
	QUESTION,
	EXCLAMATION,
	COMMA,
	COLON,
	SEMICOLON,
	HYPHEN,
	LPAREN,
	RPAREN,
	LSQUARE_BRACE,
	RSQUARE_BRACE,
	LANGLE_BRACKET,
	RANGLE_BRACKET,
	SINGLE_QUOTE,
	DOUBLE_QUOTE,
	TOKEN,
	FORMULA
};

typedef lexical_token<article_token_type> article_token;

template<typename Stream>
inline bool print(article_token_type type, Stream& stream) {
	switch (type) {
	case article_token_type::PERIOD:
		return print('.', stream);
	case article_token_type::QUESTION:
		return print('?', stream);
	case article_token_type::EXCLAMATION:
		return print('!', stream);
	case article_token_type::COMMA:
		return print(',', stream);
	case article_token_type::COLON:
		return print(':', stream);
	case article_token_type::SEMICOLON:
		return print(';', stream);
	case article_token_type::HYPHEN:
		return print('-', stream);
	case article_token_type::LPAREN:
		return print('(', stream);
	case article_token_type::RPAREN:
		return print(')', stream);
	case article_token_type::LSQUARE_BRACE:
		return print('[', stream);
	case article_token_type::RSQUARE_BRACE:
		return print(']', stream);
	case article_token_type::LANGLE_BRACKET:
		return print('<', stream);
	case article_token_type::RANGLE_BRACKET:
		return print('>', stream);
	case article_token_type::SINGLE_QUOTE:
		return print('\'', stream);
	case article_token_type::DOUBLE_QUOTE:
		return print('"', stream);
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
	FORMULA,
	COMMENT
};

bool article_emit_symbol(array<article_token>& tokens, const position& start, char symbol) {
	switch (symbol) {
	case '.': return emit_token(tokens, start, start + 1, article_token_type::PERIOD);
	case '?': return emit_token(tokens, start, start + 1, article_token_type::QUESTION);
	case '!': return emit_token(tokens, start, start + 1, article_token_type::EXCLAMATION);
	case ',': return emit_token(tokens, start, start + 1, article_token_type::COMMA);
	case '[': return emit_token(tokens, start, start + 1, article_token_type::LSQUARE_BRACE);
	case ']': return emit_token(tokens, start, start + 1, article_token_type::RSQUARE_BRACE);
	case '<': return emit_token(tokens, start, start + 1, article_token_type::LANGLE_BRACKET);
	case '>': return emit_token(tokens, start, start + 1, article_token_type::RANGLE_BRACKET);
	case '\'': return emit_token(tokens, start, start + 1, article_token_type::SINGLE_QUOTE);
	case '"': return emit_token(tokens, start, start + 1, article_token_type::DOUBLE_QUOTE);
	case ':': return emit_token(tokens, start, start + 1, article_token_type::COLON);
	case ';': return emit_token(tokens, start, start + 1, article_token_type::SEMICOLON);
	case '-': return emit_token(tokens, start, start + 1, article_token_type::HYPHEN);
	case '(': return emit_token(tokens, start, start + 1, article_token_type::LPAREN);
	case ')': return emit_token(tokens, start, start + 1, article_token_type::RPAREN);
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

	mbstate_t shift = {0};
	buffered_stream<MB_LEN_MAX, Stream> wrapper(input);
	char32_t next = fgetc32(wrapper);
	bool new_line = false;
	while (next != static_cast<char32_t>(-1)) {
		switch (state) {
		case article_lexer_state::TOKEN:
			if (next == '.' || next == '?' || next == '!' || next == ',' || next == ':' || next == ';' || next == '-'
			 || next == '(' || next == ')' || next == '[' || next == ']' || next == '<' || next == '>' || next == '\'' || next == '"')
			{
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
			} else if (next == '#') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN))
					return false;
				state = article_lexer_state::COMMENT;
				token.clear(); shift = {0};
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, article_token_type::TOKEN))
					return false;
				state = article_lexer_state::DEFAULT;
				token.clear(); shift = {0};
				new_line = (next == '\n');
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

		case article_lexer_state::COMMENT:
			if (next == '\n') {
				state = article_lexer_state::DEFAULT;
				new_line = true;
			}
			break;

		case article_lexer_state::DEFAULT:
			if (next == '.' || next == '?' || next == '!' || next == ',' || next == ':' || next == ';' || next == '-'
			 || next == '(' || next == ')' || next == '[' || next == ']' || next == '<' || next == '>' || next == '\'' || next == '"')
			{
				if (!article_emit_symbol(tokens, current, next))
					return false;
			} else if (next == '{') {
				state = article_lexer_state::FORMULA;
				token.clear(); shift = {0};
				start = current;
			} else if (next == '#') {
				state = article_lexer_state::COMMENT;
				token.clear(); shift = {0};
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
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
		next = fgetc32(wrapper);
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

inline bool article_interpret_sentence_token(
		const array<article_token>& tokens,
		unsigned int& index,
		array<sentence_token>& sentence_tokens,
		hash_map<string, unsigned int>& names,
		unsigned int& token_id)
{
	if (index >= tokens.length) {
		read_error("Unexpected end of input. Expected a sentence token", tokens.last().end);
		return false;
	} else if (tokens[index].type == article_token_type::TOKEN) {
		if (!get_token(tokens[index].text, token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::PERIOD) {
		if (!get_token(".", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::QUESTION) {
		if (!get_token("?", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::EXCLAMATION) {
		if (!get_token("!", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::COMMA) {
		if (!get_token(",", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::HYPHEN) {
		if (!get_token("-", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::LPAREN) {
		if (!get_token("(", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::RPAREN) {
		if (!get_token(")", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::LSQUARE_BRACE) {
		if (!get_token("[", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::RSQUARE_BRACE) {
		if (!get_token("]", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::LANGLE_BRACKET) {
		if (!get_token("<", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::RANGLE_BRACKET) {
		if (!get_token(">", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::SINGLE_QUOTE) {
		if (!get_token("'", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else if (tokens[index].type == article_token_type::DOUBLE_QUOTE) {
		if (!get_token("\"", token_id, names) || !sentence_tokens.add({token_id}))
			return false;
	} else {
		read_error("Expected a sentence token", tokens[index].start);
		return false;
	}
	index++;
	return true;
}

template<typename Derivation>
inline bool interpret_derivation_tree(
		const array<article_token>& tokens,
		unsigned int& index, sentence<Derivation>& out,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	memory_stream& in = *((memory_stream*) alloca(sizeof(memory_stream)));
	in.buffer = tokens[index].text.data;
	in.length = tokens[index].text.length;
	in.position = 0;
	Derivation derivation;
	if (!parse(in, derivation, names, nonterminal_names, tokens[index].start + 1))
		return false;
	index++;

	out.derivation = derivation;
	out.tokens = nullptr;
	out.length = 0;
	return true;
}

template<typename Derivation>
bool article_interpret_sentence(
		const array<article_token>& tokens,
		unsigned int& index, sentence<Derivation>& out,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	if (index < tokens.length && tokens[index].type == article_token_type::FORMULA) {
		/* this is a derivation tree, so we parse it */
		return interpret_derivation_tree(tokens, index, out, names, nonterminal_names);
	}

	array<sentence_token>& sentence_tokens = *((array<sentence_token>*) alloca(sizeof(array<sentence_token>)));
	if (!array_init(sentence_tokens, 16)) return false;
	while (true) {
		unsigned int token_id;
		if (!article_interpret_sentence_token(tokens, index, sentence_tokens, names, token_id)) {
			free(sentence_tokens);
			return false;
		}

		if (index < tokens.length && tokens[index].type == article_token_type::COLON) {
			index++;
			if (!expect_token(tokens, index, article_token_type::TOKEN, "token label")) {
				free(sentence_tokens); return false;
			}
			/* ignore the token label */
			index++;
		}

		if (index >= tokens.length || (tokens[index].type != article_token_type::TOKEN && tokens[index].type != article_token_type::HYPHEN
		 && tokens[index].type != article_token_type::LPAREN && tokens[index].type != article_token_type::RPAREN
		 && tokens[index].type != article_token_type::LSQUARE_BRACE && tokens[index].type != article_token_type::RSQUARE_BRACE
		 && tokens[index].type != article_token_type::LANGLE_BRACKET && tokens[index].type != article_token_type::RANGLE_BRACKET
		 && tokens[index].type != article_token_type::SINGLE_QUOTE && tokens[index].type != article_token_type::DOUBLE_QUOTE
		 && tokens[index].type != article_token_type::PERIOD && tokens[index].type != article_token_type::COMMA
		 && tokens[index].type != article_token_type::QUESTION && tokens[index].type != article_token_type::EXCLAMATION))
			break;
	}

	move(sentence_tokens.data, out.tokens);
	out.length = sentence_tokens.length;
	set_empty(out.derivation);
	return true;
}

template<typename Derivation>
bool article_interpret_sentence(
		const array<article_token>& tokens,
		unsigned int& index, sentence<Derivation>& out,
		array_map<unsigned int, unsigned int>& labels,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	if (index < tokens.length && tokens[index].type == article_token_type::FORMULA) {
		/* this is a derivation tree, so we parse it */
		return interpret_derivation_tree(tokens, index, out, names, nonterminal_names);
	}

	array<sentence_token>& sentence_tokens = *((array<sentence_token>*) alloca(sizeof(array<sentence_token>)));
	if (!array_init(sentence_tokens, 16)) return false;
	while (true) {
		unsigned int token_id;
		if (!article_interpret_sentence_token(tokens, index, sentence_tokens, names, token_id)) {
			free(sentence_tokens);
			return false;
		}

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

		if (index >= tokens.length || (tokens[index].type != article_token_type::TOKEN && tokens[index].type != article_token_type::PERIOD
		 && tokens[index].type != article_token_type::QUESTION && tokens[index].type != article_token_type::EXCLAMATION))
			break;
	}

	move(sentence_tokens.data, out.tokens);
	out.length = sentence_tokens.length;
	set_empty(out.derivation);
	return true;
}

template<typename Formula>
struct sentence_label {
	Formula* logical_form;
	array_map<unsigned int, unsigned int> labels;

	static inline void move(const sentence_label<Formula>& src, sentence_label<Formula>& dst) {
		dst.logical_form = src.logical_form;
		core::move(src.labels, dst.labels);
	}

	static inline void free(sentence_label<Formula>& label) {
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

template<typename Formula>
inline bool init(sentence_label<Formula>& label) {
	label.logical_form = (Formula*) malloc(sizeof(Formula));
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

template<typename Formula>
inline bool operator != (const sentence_label<Formula>& first, const sentence_label<Formula>& second) {
	if (*first.logical_form != *second.logical_form)
		return true;
	if (first.labels.size != second.labels.size)
		return true;
	for (unsigned int i = 0; i < first.labels.size; i++)
		if (first.labels.keys[i] != second.labels.keys[i]
		 || first.labels.values[i] != second.labels.values[i]) return true;
	return false;
}

/* NOTE: we assume that `index < tokens.length` */
template<typename Derivation>
bool article_interpret(
		const array<article_token>& tokens,
		unsigned int& index, article<Derivation>& out,
		unsigned int& article_name,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	unsigned int start = index;
	while (tokens[index].type != article_token_type::COLON) {
		if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
			return false;
		index++;
	}

	if (start == index) {
		/* the name is empty so generate a unique name */
		if (!names.check_size()) return false;
		for (unsigned int article_counter = article_name; ; article_counter++) {
			int length = snprintf(NULL, 0, "%u", article_counter);
			if (length < 0) return false;

			string new_name = string(length + 1);
			snprintf(new_name.data, length + 1, "%u", article_counter);
			new_name.length = length;

			bool contains; unsigned int bucket;
			article_name = names.get(new_name, contains, bucket);
			if (!contains) {
				names.table.keys[bucket] = new_name;
				names.table.size++;
				names.values[bucket] = names.table.size;
				break;
			}
		}
	} else {
		string name = string(16);
		if (!concatenate(tokens.data + start, index - start, name)) {
			fprintf(stderr, "article_interpret ERROR: Unable to construct string for article name.\n");
			return false;
		} else if (!get_token(name, article_name, names)) {
			return false;
		}
	}
	index++;

	array<sentence<Derivation>>& sentences = *((array<sentence<Derivation>>*) alloca(sizeof(array<sentence<Derivation>>)));
	if (!array_init(sentences, 8))
		return false;
	while (true) {
		if (!sentences.ensure_capacity(sentences.length + 1)
		 || !article_interpret_sentence(tokens, index, sentences[sentences.length], names, nonterminal_names))
		{
			for (sentence<Derivation>& s : sentences) free(s);
			free(sentences); return false;
		}
		sentences.length++;

		if (index < tokens.length && tokens[index].type == article_token_type::FORMULA) {
			/* skip the logical form */
			index++;
		}

		if (index >= tokens.length) {
			break;
		} else if (tokens[index].type == article_token_type::SEMICOLON) {
			index++; break;
		}
	}

	move(sentences.data, out.sentences);
	out.sentence_count = sentences.length;
	return true;
}

/* NOTE: we assume that `index < tokens.length` */
template<typename Derivation, typename Formula>
bool article_interpret(
		const array<article_token>& tokens,
		unsigned int& index, article<Derivation>& out,
		unsigned int& article_name,
		hash_map<sentence<Derivation>, sentence_label<Formula>>& logical_forms,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	unsigned int start = index;
	while (tokens[index].type != article_token_type::COLON) {
		if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
			return false;
		index++;
	}

	if (start == index) {
		/* the name is empty so generate a unique name */
		if (!names.check_size()) return false;
		for (unsigned int article_counter = article_name; ; article_counter++) {
			int length = snprintf(NULL, 0, "%u", article_counter);
			if (length < 0) return false;

			string new_name = string(length + 1);
			snprintf(new_name.data, length + 1, "%u", article_counter);
			new_name.length = length;

			bool contains; unsigned int bucket;
			article_name = names.get(new_name, contains, bucket);
			if (!contains) {
				names.table.keys[bucket] = new_name;
				names.table.size++;
				names.values[bucket] = names.table.size;
				break;
			}
		}
	} else {
		string name = string(16);
		if (!concatenate(tokens.data + start, index - start, name)) {
			fprintf(stderr, "article_interpret ERROR: Unable to construct string for article name.\n");
			return false;
		} else if (!get_token(name, article_name, names)) {
			return false;
		}
	}
	index++;

	array<sentence<Derivation>>& sentences = *((array<sentence<Derivation>>*) alloca(sizeof(array<sentence<Derivation>>)));
	if (!array_init(sentences, 8))
		return false;
	while (true) {
		sentence_label<Formula>& new_label = *((sentence_label<Formula>*) alloca(sizeof(sentence_label<Formula>)));
		if (!init(new_label)) {
			for (sentence<Derivation>& s : sentences) free(s);
			free(sentences); return false;
		} else if (!sentences.ensure_capacity(sentences.length + 1)
				|| !article_interpret_sentence(tokens, index, sentences[sentences.length], new_label.labels, names, nonterminal_names))
		{
			for (sentence<Derivation>& s : sentences) free(s);
			free(sentences); free(new_label); return false;
		}
		sentences.length++;

		if (index < tokens.length && tokens[index].type == article_token_type::FORMULA) {
			/* parse the formula */
			if (!logical_forms.check_size()
			 || !parse(tokens[index].text.data, tokens[index].text.length, *new_label.logical_form, names, tokens[index].start + 1)) {
				for (sentence<Derivation>& s : sentences) free(s);
				free(sentences); free(new_label);
				return false;
			}

			bool contains; unsigned int bucket;
			sentence_label<Formula>& value = logical_forms.get(sentences.last(), contains, bucket);
			if (!contains) {
				if (!init(logical_forms.table.keys[bucket], sentences.last())) {
					for (sentence<Derivation>& s : sentences) free(s);
					free(sentences); free(new_label);
					return false;
				}
				move(new_label, value);
				logical_forms.table.size++;
			} else if (value != new_label) {
				read_error("A different logical form was previously mapped to this sentence", tokens[index].start);
				for (sentence<Derivation>& s : sentences) free(s);
				free(sentences); free(new_label);
				return false;
			}
			index++;
		}

		if (index >= tokens.length) {
			break;
		} else if (tokens[index].type == article_token_type::SEMICOLON) {
			index++; break;
		}
	}

	move(sentences.data, out.sentences);
	out.sentence_count = sentences.length;
	return true;
}

/* NOTE: we assume that `index < tokens.length` */
template<typename Derivation, typename Formula>
bool article_interpret(
		const array<article_token>& tokens,
		unsigned int& index,
		array_map<sentence<Derivation>, Formula>& logical_forms,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	unsigned int start = index;
	while (tokens[index].type != article_token_type::COLON) {
		if (!expect_token(tokens, index, article_token_type::TOKEN, "article name"))
			return false;
		index++;
	}

	unsigned int article_name;
	if (start == index) {
		/* the name is empty so generate a unique name */
		if (!names.check_size()) return false;
		for (unsigned int article_counter = index; ; article_counter++) {
			int length = snprintf(NULL, 0, "%u", article_counter);
			if (length < 0) return false;

			string new_name = string(length + 1);
			snprintf(new_name.data, length + 1, "%u", article_counter);
			new_name.length = length;

			bool contains; unsigned int bucket;
			article_name = names.get(new_name, contains, bucket);
			if (!contains) {
				names.table.keys[bucket] = new_name;
				names.table.size++;
				names.values[bucket] = names.table.size;
				article_name = names.table.size;
				break;
			}
		}
	} else {
		string name = string(16);
		if (!concatenate(tokens.data + start, index - start, name)) {
			fprintf(stderr, "article_interpret ERROR: Unable to construct string for article name.\n");
			return false;
		} else if (!get_token(name, article_name, names)) {
			return false;
		}
	}
	index++;

	while (true) {
		array_map<unsigned int, unsigned int> dummy_labels(1);
		if (!logical_forms.ensure_capacity(logical_forms.size + 1)
		 || !article_interpret_sentence(tokens, index, logical_forms.keys[logical_forms.size], dummy_labels, names, nonterminal_names))
		{
			return false;
		}

		if (!expect_token(tokens, index, article_token_type::FORMULA, "logical form label")) {
			free(logical_forms.keys[logical_forms.size]);
			return false;
		}

		/* parse the formula */
		if (!parse(tokens[index].text.data, tokens[index].text.length, logical_forms.values[logical_forms.size], names, tokens[index].start + 1)) {
			free(logical_forms.keys[logical_forms.size]);
			return false;
		}
		logical_forms.size++;
		index++;

		if (index >= tokens.length) {
			break;
		} else if (tokens[index].type == article_token_type::SEMICOLON) {
			index++; break;
		}
	}

	return true;
}

template<typename Derivation, typename Formula>
bool articles_interpret(
		const array<article_token>& tokens,
		in_memory_article_store<Derivation>& articles,
		hash_map<sentence<Derivation>, sentence_label<Formula>>& logical_forms,
		hash_map<string, unsigned int>& names,
		const hash_map<string, unsigned int>& nonterminal_names)
{
	unsigned int index = 0;
	while (index < tokens.length) {
		if (!articles.articles.check_size())
			return false;

		position article_pos = tokens[index].start;

		unsigned int article_name = articles.articles.table.size + 1;
		article<Derivation>& new_article = *((article<Derivation>*) alloca(sizeof(article<Derivation>)));
		if (!article_interpret(tokens, index, new_article, article_name, logical_forms, names, nonterminal_names))
			return false;

		bool contains; unsigned int bucket;
		article<Derivation>& value = articles.articles.get(article_name, contains, bucket);
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

/* a utility function for tokenizing individual sentences represented as C strings */
template<typename Derivation>
inline bool tokenize(
		const char* input, sentence<Derivation>& out,
		hash_map<string, unsigned int>& names)
{
	memory_stream in(input, strlen(input));

	array<article_token> tokens(256);
	if (!article_lex(tokens, in)) {
		free_tokens(tokens);
		return false;
	}

	if (tokens.length == 0) {
		fprintf(stderr, "ERROR: Input is empty.\n");
		return false;
	}

	unsigned int index = 0;
	hash_map<string, unsigned int> dummy(1);
	if (!article_interpret_sentence(tokens, index, out, names, dummy)) {
		free_tokens(tokens);
		return false;
	}
	free_tokens(tokens);
	return true;
}

/* a utility function for tokenizing individual sentences represented as a collection of tokens */
template<typename TokenType, typename Derivation>
inline bool tokenize(
		const TokenType* input, unsigned int input_length,
		sentence<Derivation>& out,
		hash_map<string, unsigned int>& names)
{
	array<sentence_token>& sentence_tokens = *((array<sentence_token>*) alloca(sizeof(array<sentence_token>)));
	if (!array_init(sentence_tokens, 16)) return false;
	for (unsigned int i = 0; i < input_length; i++) {
		unsigned int token_id;
		if (!get_token(string(input[i].text, input[i].length), token_id, names) || !sentence_tokens.add({token_id})) {
			free(sentence_tokens);
			return false;
		}
	}

	move(sentence_tokens.data, out.tokens);
	out.length = sentence_tokens.length;
	set_empty(out.derivation);
	return true;
}

#endif /* ARTICLE_H_ */
