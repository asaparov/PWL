/**
 * first_order_logic.h
 *
 *  Created on: Jun 26, 2018
 *      Author: asaparov
 */

#ifndef FIRST_ORDER_LOGIC_H_
#define FIRST_ORDER_LOGIC_H_

#include <core/lex.h>

using namespace core;


/* forward declarations */

struct fol_expression;

enum class fol_term_type {
	VARIABLE,
	CONSTANT
};

struct fol_term {
	fol_term_type type;
	union {
		unsigned int variable;
		unsigned int constant;
	};
};

enum class fol_expression_type {
	ATOM,

	AND,
	OR,
	IF_THEN,
	IFF,
	NOT,

	FOR_ALL,
	EXISTS
};

struct fol_atom {
	unsigned int predicate;
	fol_term arg1;
	fol_term arg2;
};

struct fol_and_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_or_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_if_then_expression {
	fol_expression* antecedent;
	fol_expression* consequent;
};

struct fol_iff_expression {
	fol_expression* left;
	fol_expression* right;
};

struct fol_not_expression {
	fol_expression* operand;
};

struct fol_for_all_expression {
	unsigned int variable;
	fol_expression* operand;
};

struct fol_exists_expression {
	unsigned int variable;
	fol_expression* operand;
};

struct fol_expression {
	fol_expression_type type;
	unsigned int reference_count;
	union {
		fol_atom atom;
		fol_and_expression conj;
		fol_or_expression disj;
		fol_if_then_expression if_then;
		fol_iff_expression iff;
		fol_for_all_expression for_all;
		fol_exists_expression exists;
	};
};


/**
 * Code for tokenizing/lexing first-order logic formulas in TPTP-like format.
 */

enum class tptp_token_type {
	LBRACKET,
	RBRACKET,
	LPAREN,
	RPAREN,
	COMMA,
	COLON,

	AND,
	OR,
	NOT,
	IF_THEN,
	IFF,
	FOR_ALL,
	EXISTS,

	IDENTIFIER
};

typedef lexical_token<tptp_token_type> tptp_token;

template<typename Stream>
inline bool print(tptp_token_type type, Stream& stream) {
	switch (type) {
	case tptp_token_type::LBRACKET:
		return print('[', stream);
	case tptp_token_type::RBRACKET:
		return print(']', stream);
	case tptp_token_type::LPAREN:
		return print('(', stream);
	case tptp_token_type::RPAREN:
		return print(')', stream);
	case tptp_token_type::COMMA:
		return print(',', stream);
	case tptp_token_type::COLON:
		return print(':', stream);
	case tptp_token_type::AND:
		return print('&', stream);
	case tptp_token_type::OR:
		return print('|', stream);
	case tptp_token_type::NOT:
		return print('~', stream);
	case tptp_token_type::IF_THEN:
		return print("=>", stream);
	case tptp_token_type::IFF:
		return print("<=>", stream);
	case tptp_token_type::FOR_ALL:
		return print('!', stream);
	case tptp_token_type::EXISTS:
		return print('?', stream);
	case tptp_token_type::IDENTIFIER:
		return print("IDENTIFIER", stream);
	}
	fprintf(stderr, "print ERROR: Unknown tptp_token_type.\n");
	return false;
}

enum class tptp_lexer_state {
	DEFAULT,
	IDENTIFIER,
}

bool datalog_lex(array<tptp_token>& tokens, FILE* input) {
	position start = position(1, 1);
	position current = position(1, 1);
	tptp_lexer_state state = tptp_lexer_state::DEFAULT;
	array<char> token = array<char>(1024);

	int prev = 0, next = fgetc(input);
	bool new_line = false;
	while (next != -1) {
		switch (state) {
		case tptp_lexer_state::IDENTIFIER:
			if (next == ',' || next == ':' || next == '(' || next == ')'
			 || next == '[' || next == ']' || next == '&' || next == '|'
			 || next == '~')
			{
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER)
				 || !tptp_emit_symbol(tokens, current, next))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear();
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER))
					return false;
				state = tptp_lexer_state::DEFAULT;
				token.clear();
				new_line = (next == '\n');
			} else {
				if (!token.add(next)) return false;
			}
			break;

		case tptp_lexer_state::DEFAULT:
			if (next == ',' || next == ':' || next == '(' || next == ')'
			 || next == '[' || next == ']' || next == '&' || next == '|'
			 || next == '~')
			{
				if (!tptp_emit_symbol(tokens, current, next))
					return false;
			} else if (next == '\\') {
				next = fgetc(input);
				if (next != '+') {
					read_error("Expected '+'", current);
					return false;
				} if (!emit_token(tokens, current, current + 1, DATALOG_TOKEN_SLASH_PLUS))
					return false;
				current.column++;
			} else if (next == ' ' || next == '\t' || next == '\n' || next == '\r') {
				new_line = (next == '\n');
			} else {
				if (!token.add(next)) return false;
				state = tptp_lexer_state::IDENTIFIER;
				start = current;
			}
			break;
		}

		if (new_line) {
			current.line++;
			current.column = 1;
			new_line = false;
		} else current.column++;
		prev = next;
		next = fgetc(input);
	}

	if (state == tptp_lexer_state::IDENTIFIER) {
		return emit_token(tokens, token, start, current, tptp_lexer_state::IDENTIFIER);
	}
	return true;
}

#endif /* FIRST_ORDER_LOGIC_H_ */
