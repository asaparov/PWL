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

struct fol_formula;

enum class fol_term_type {
	VARIABLE,
	CONSTANT,
	NONE
};

struct fol_term {
	fol_term_type type;
	union {
		unsigned int variable;
		unsigned int constant;
	};
};

enum class fol_formula_type {
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

	static inline void free(fol_and_formula& formula) { }
};

struct fol_and_formula {
	fol_formula* left;
	fol_formula* right;

	static inline void free(fol_and_formula& formula) {
		core::free(*formula.left);
		if (formula.left->reference_count == 0)
			core::free(formula.left);

		core::free(*formula.right);
		if (formula.right->reference_count == 0)
			core::free(formula.right);
	}
};

struct fol_or_formula {
	fol_formula* left;
	fol_formula* right;

	static inline void free(fol_or_formula& formula) {
		core::free(*formula.left);
		if (formula.left->reference_count == 0)
			core::free(formula.left);

		core::free(*formula.right);
		if (formula.right->reference_count == 0)
			core::free(formula.right);
	}
};

struct fol_if_then_formula {
	fol_formula* antecedent;
	fol_formula* consequent;

	static inline void free(fol_if_then_formula& formula) {
		core::free(*formula.antecedent);
		if (formula.antecedent->reference_count == 0)
			core::free(formula.antecedent);

		core::free(*formula.consequent);
		if (formula.consequent->reference_count == 0)
			core::free(formula.consequent);
	}
};

struct fol_iff_formula {
	fol_formula* left;
	fol_formula* right;

	static inline void free(fol_iff_formula& formula) {
		core::free(*formula.left);
		if (formula.left->reference_count == 0)
			core::free(formula.left);

		core::free(*formula.right);
		if (formula.right->reference_count == 0)
			core::free(formula.right);
	}
};

struct fol_not_formula {
	fol_formula* operand;

	static inline void free(fol_not_formula& formula) {
		core::free(*formula.operand);
		if (formula.operand->reference_count == 0)
			core::free(formula.operand);
	}
};

struct fol_for_all_formula {
	unsigned int variable;
	fol_formula* operand;

	static inline void free(fol_for_all_formula& formula) {
		core::free(*formula.operand);
		if (formula.operand->reference_count == 0)
			core::free(formula.operand);
	}
};

struct fol_exists_formula {
	unsigned int variable;
	fol_formula* operand;

	static inline void free(fol_exists_formula& formula) {
		core::free(*formula.operand);
		if (formula.operand->reference_count == 0)
			core::free(formula.operand);
	}
};

struct fol_formula {
	fol_formula_type type;
	unsigned int reference_count;
	union {
		fol_atom atom;
		fol_and_formula conj;
		fol_or_formula disj;
		fol_if_then_formula if_then;
		fol_iff_formula iff;
		fol_for_all_formula for_all;
		fol_exists_formula exists;
		fol_formula* negation;
	};

	static inline void free(fol_formula& formula) {
		formula.reference_count--;
		if (formula.reference_count == 0) {
			switch (formula.type) {
			case fol_formula_type::ATOM:
				core::free(formula.atom); break;
			case fol_formula_type::AND:
				core::free(formula.conj); break;
			case fol_formula_type::OR:
				core::free(formula.disj); break;
			case fol_formula_type::IF_THEN:
				core::free(formula.if_then); break;
			case fol_formula_type::IFF:
				core::free(formula.iff); break;
			case fol_formula_type::FOR_ALL:
				core::free(formula.for_all); break;
			case fol_formula_type::EXISTS:
				core::free(formula.exists); break;
			case fol_formula_type::NOT:
				core::free(formula.negation); break;
			}
			fprintf(stderr, "fol_formula.free ERROR: Unrecognized fol_formula_type.\n");
		}
	}
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
};

bool tptp_emit_symbol(array<tptp_token>& tokens, const position& start, char symbol) {
	switch (symbol) {
	case ',':
		return emit_token(tokens, start, start + 1, tptp_token_type::COMMA);
	case ':':
		return emit_token(tokens, start, start + 1, tptp_token_type::COLON);
	case '(':
		return emit_token(tokens, start, start + 1, tptp_token_type::LPAREN);
	case ')':
		return emit_token(tokens, start, start + 1, tptp_token_type::RPAREN);
	case '[':
		return emit_token(tokens, start, start + 1, tptp_token_type::LBRACKET);
	case ']':
		return emit_token(tokens, start, start + 1, tptp_token_type::RBRACKET);
	case '&':
		return emit_token(tokens, start, start + 1, tptp_token_type::AND);
	case '|':
		return emit_token(tokens, start, start + 1, tptp_token_type::OR);
	case '~':
		return emit_token(tokens, start, start + 1, tptp_token_type::NOT);
	case '!':
		return emit_token(tokens, start, start + 1, tptp_token_type::FOR_ALL);
	case '?':
		return emit_token(tokens, start, start + 1, tptp_token_type::EXISTS);
	default:
		fprintf(stderr, "tptp_emit_symbol ERROR: Unexpected symbol.\n");
		return false;
	}
}

inline bool tptp_lex_symbol(array<tptp_token>& tokens, FILE* input, int next, position& current)
{
	if (next == ',' || next == ':' || next == '(' || next == ')'
	 || next == '[' || next == ']' || next == '&' || next == '|'
	 || next == '~' || next == '!' || next == '?')
	{
		return tptp_emit_symbol(tokens, current, next);
	} else if (next == '=') {
		next = fgetc(input);
		if (next != '>') {
			read_error("Expected '>' after '='", current);
			return false;
		} if (!emit_token(tokens, current, current + 2, tptp_token_type::IF_THEN))
			return false;
		current.column++;
	} else if (next == '<') {
		next = fgetc(input);
		if (next != '=') {
			read_error("Expected '=' after '<'", current);
			return false;
		}
		next = fgetc(input);
		if (next != '>') {
			read_error("Expected '>' after '='", current);
			return false;
		} if (!emit_token(tokens, current, current + 3, tptp_token_type::IF_THEN))
			return false;
		current.column += 2;
	} else {
		fprintf(stderr, "tptp_lex_symbol ERROR: Unrecognized symbol.\n");
		return false;
	}
	return true;
}

bool tptp_lex(array<tptp_token>& tokens, FILE* input) {
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
			 || next == '~' || next == '!' || next == '?' || next == '='
			 || next == '<')
			{
				if (!emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER)
				 || !tptp_lex_symbol(tokens, input, next, current))
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
			 || next == '~' || next == '!' || next == '?' || next == '='
			 || next == '<')
			{
				if (!tptp_lex_symbol(tokens, input, next, current))
					return false;
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

	if (state == tptp_lexer_state::IDENTIFIER)
		return emit_token(tokens, token, start, current, tptp_token_type::IDENTIFIER);
	return true;
}


/**
 * Recursive-descent parser for first-order logic formulas in TPTP-like format.
 */

bool tptp_interpret_argument_list(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hash_map<string, unsigned int>& names,
	hash_map<string, unsigned int>& variables,
	array<fol_term>& terms)
{
	if (!expect_token(tokens, index, tptp_token_type::LPAREN,
			"opening parenthesis for list of arguments in atomic formula"))
		return false;
	index++;

	while (true) {
		if (terms.ensure_capacity(terms.length + 1)
		 || !expect_token(tokens, index, tptp_token_type::IDENTIFIER,
				"identifier in list of arguments in atomic formula"))
			return false;

		bool contains;
		fol_term& next_term = terms[terms.length];
		unsigned int variable = variables.get(tokens[index].text, contains);
		if (contains) {
			/* this argument is a variable */
			next_term.variable = variable;
			next_term.type = fol_term_type::VARIABLE;
		} else {
			if (!get_token(tokens[index].text, next_term.constant, names))
				return false;
			next_term.type = fol_term_type::CONSTANT;
		}
		terms.length++;
		index++;

		if (index >= tokens.length) {
			read_error("Unexpected end of input", tokens.last().end);
			return false;
		} else if (tokens[index].type == tptp_token_type::RPAREN) {
			index++;
			return true;
		} else if (tokens[index].type != tptp_token_type::COMMA) {
			read_error("Unexpected symbol. Expected a comma", tokens[index].start);
			return false;
		}
		index++;
	}
}

bool tptp_interpret_variable_list(
	const array<tptp_token>& tokens,
	unsigned int& index,
	hash_map<string, unsigned int>& names,
	hash_map<string, unsigned int>& variables,
	array<unsigned int>& var_list)
{
	if (!expect_token(tokens, index, tptp_token_type::LBRACKET, "left bracket for list of quantified variables"))
		return false;
	index++;

	while (true) {
		if (!expect_token(tokens, index, tptp_token_type::IDENTIFIER, "variable in list of quantified variables"))
			return false;
		unsigned int var;
		if (names.table.contains(tokens[index].text)) {
			fprintf(stderr, "WARNING at %d:%d: Variable '", tokens[index].start.line, tokens[index].start.column);
			print(tokens[index].text, stderr); print("' shadows previously declared identifier.\n", stderr);
		} if (variables.table.contains(tokens[index].text)) {
			read_error("Variable redeclared", tokens[index].start);
		} if (!get_token(tokens[index].text, var, variables) || var_list.add(var)) {
			return false;
		}
		index++;

		if (index >= tokens.length) {
			read_error("Unexpected end of input", tokens.last().end);
			return false;
		} else if (tokens[index].type == tptp_token_type::RBRACKET) {
			index++;
			return true;
		} else if (tokens[index].type != tptp_token_type::COMMA) {
			read_error("Unexpected symbol. Expected a comma", tokens[index].start);
			return false;
		}
		index++;
	}
}

bool tptp_interpret_unary_formula(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names,
	hash_map<string, unsigned int>& variables)
{
	if (index >= tokens.length) {
		fprintf(stderr, "ERROR: Unexpected end of input.\n");
		return false;

	} else if (tokens[index].type == tptp_token_type::NOT) {
		/* this is a negation of the form ~U */
		index++;
		fol_formula* operand = (fol_formula*) malloc(sizeof(fol_formula));
		operand->reference_count = 1;
		if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
			free(operand); return false;
		}

		formula.type = fol_formula_type::NOT;
		formula.negation = operand;

	} else if (tokens[index].type == tptp_token_type::LPAREN) {
		/* these are just grouping parenthesis of the form (U) */
		index++;
		if (!tptp_interpret_unary_formula(tokens, index, formula, names, variables)) {
			return false;
		} if (!expect_token(tokens, index, tptp_token_type::RPAREN, "closing parenthesis")) {
			free(formula); return false;
		}
		index++;

	} else if (tokens[index].type == tptp_token_type::FOR_ALL) {
		/* this is a universal quantifier of the form ![v_1,...,v_n]:U */
		index++;
		array<unsigned int> var_list = array<unsigned int>(4);
		if (!tptp_interpret_variable_list(tokens, index, names, variables, var_list)
		 || !expect_token(tokens, index, tptp_token_type::COLON, "colon for univerally-quantified formula"))
			return false;
		index++;

		fol_formula* operand = (fol_formula*) malloc(sizeof(fol_formula));
		operand->reference_count = 1;
		if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
			free(operand); return false;
		}

		fol_formula* inner = operand;
		for (unsigned int i = var_list.length - 1; i > 0; i--) {
			fol_formula* quantified = (fol_formula*) malloc(sizeof(fol_formula));
			if (quantified == NULL) {
				fprintf(stderr, "tptp_interpret_unary_formula ERROR: Out of memory.\n");
				free(*inner); free(inner);
			}
			quantified->for_all.variable = var_list[i];
			quantified->for_all.operand = inner;
			quantified->type = fol_formula_type::FOR_ALL;
			quantified->reference_count = 1;
			inner = quantified;
		}

		formula.for_all.variable = var_list[0];
		formula.for_all.operand = inner;
		formula.type = fol_formula_type::FOR_ALL;
		formula.reference_count = 1;

	} else if (tokens[index].type == tptp_token_type::EXISTS) {
		/* this is an existential quantifier of the form ?[v_1,...,v_n]:U */
		index++;
		array<unsigned int> var_list = array<unsigned int>(4);
		if (!tptp_interpret_variable_list(tokens, index, names, variables, var_list)
		 || !expect_token(tokens, index, tptp_token_type::COLON, "colon for existentially-quantified formula"))
			return false;
		index++;

		fol_formula* operand = (fol_formula*) malloc(sizeof(fol_formula));
		operand->reference_count = 1;
		if (!tptp_interpret_unary_formula(tokens, index, *operand, names, variables)) {
			free(operand); return false;
		}

		fol_formula* inner = operand;
		for (unsigned int i = var_list.length - 1; i > 0; i--) {
			fol_formula* quantified = (fol_formula*) malloc(sizeof(fol_formula));
			if (quantified == NULL) {
				fprintf(stderr, "tptp_interpret_unary_formula ERROR: Out of memory.\n");
				free(*inner); free(inner);
			}
			quantified->exists.variable = var_list[i];
			quantified->exists.operand = inner;
			quantified->type = fol_formula_type::EXISTS;
			quantified->reference_count = 1;
			inner = quantified;
		}

		formula.exists.variable = var_list[0];
		formula.exists.operand = inner;
		formula.type = fol_formula_type::EXISTS;
		formula.reference_count = 1;

	} else if (tokens[index].type == tptp_token_type::IDENTIFIER) {
		/* this is an atomic formula of the form P(T_1,...,T_n) */
		unsigned int predicate;
		if (!get_token(tokens[index].text, predicate, names))
			return false;
		index++;

		array<fol_term> terms = array<fol_term>(2);
		if (!tptp_interpret_argument_list(tokens, index, names, variables, terms))
			return false;
		if (terms.length == 0) {
			formula.atom.arg1.type = fol_term_type::NONE;
			formula.atom.arg2.type = fol_term_type::NONE;
		} else if (terms.length == 1) {
			formula.atom.arg1 = terms[0];
			formula.atom.arg2.type = fol_term_type::NONE;
		} else if (terms.length == 2) {
			formula.atom.arg1 = terms[0];
			formula.atom.arg2 = terms[1];
		} else {
			read_error("Atomic formulas with arity greater than 2 are not supported", tokens[index - 1].end);
			return false;
		}
		formula.atom.predicate = predicate;
		formula.type = fol_formula_type::ATOM;

	} else {
		read_error("Unexpected symbol. Expected a unary formula", tokens[index].start);
		return false;
	}
	return true;
}

bool tptp_interpret(
	const array<tptp_token>& tokens,
	unsigned int& index, fol_formula& formula,
	hash_map<string, unsigned int>& names)
{
	fol_formula left;

}

#endif /* FIRST_ORDER_LOGIC_H_ */
