#ifndef FAKE_PARSER_H_
#define FAKE_PARSER_H_

#include "article.h"
#include "first_order_logic.h"
#include "theory.h"

struct constant_relabeler {
	const array<unsigned int>& old_constants;
	unsigned int new_constant;
};

inline bool clone_constant(unsigned int src_constant, unsigned int& dst_constant, constant_relabeler& relabeler) {
	if (relabeler.old_constants.contains(src_constant))
		dst_constant = relabeler.new_constant;
	else return clone_constant(src_constant, dst_constant);
	return true;
}

inline bool clone_predicate(unsigned int src_predicate, unsigned int& dst_predicate, constant_relabeler& relabeler) {
	return clone_predicate(src_predicate, dst_predicate);
}

inline bool clone_variable(unsigned int src_variable, unsigned int& dst_variable, constant_relabeler& relabeler) {
	return clone_variable(src_variable, dst_variable);
}

inline fol_formula* relabel_constants(const fol_formula* src,
		const array<unsigned int>& old_constants, unsigned int new_constant)
{
	fol_formula* dst;
	if (!new_fol_formula(dst)) return NULL;

	constant_relabeler relabeler = {old_constants, new_constant};
	if (!clone(*src, *dst, relabeler)) {
		free(dst); return NULL;
	}
	return dst;
}

struct constant_subtracter {
	array<unsigned int>& constants;
};

constexpr bool visit_predicate(unsigned int predicate, const constant_subtracter& subtracter) { return true; }
constexpr bool visit_variable(unsigned int variable, const constant_subtracter& subtracter) { return true; }

template<fol_formula_type Operator>
constexpr bool visit_operator(const fol_formula_type& formula, const constant_subtracter& subtracter) { return true; }

inline bool visit_constant(unsigned int constant, const constant_subtracter& subtracter) {
	unsigned int index = subtracter.constants.index_of(constant);
	if (index < subtracter.constants.length)
		subtracter.constants.remove(index);
	return true;
}

struct fake_parser {
	hash_map<sentence, fol_formula*> table;
	array<token> unknown_tokens;
	array<unsigned int> unknown_concepts;
	unsigned int unknown_id;

	fake_parser(unsigned int unknown_id) : table(64),
			unknown_tokens(16), unknown_concepts(16), unknown_id(unknown_id) { }

	~fake_parser() {
		for (auto entry : table) {
			free(entry.key);
			free(*entry.value);
			if (entry.value->reference_count == 0)
				free(entry.value);
		}
	}

	bool parse(const sentence& s, fol_formula** logical_forms,
			double* log_probabilities, unsigned int& parse_count,
			const theory& T, array<token>& unrecognized) const
	{
		for (unsigned int i = 0; i < s.length; i++) {
			if (unknown_tokens.contains(s.tokens[i])
			 && !unrecognized.add(s.tokens[i]))
				return false;
		}

		bool contains;
		const fol_formula* parse = table.get(s, contains);
		if (!contains) {
			parse_count = 0;
			return true;
		}

		fol_formula* out = relabel_constants(parse, unknown_concepts, unknown_id);
		if (out == NULL) {
			fprintf(stderr, "fake_parser.parse ERROR: Failed to relabel unknown constants.\n");
			return false;
		}
		return true;
	}

	bool is_definition_of(fol_formula* logical_form,
			unsigned int word, const theory& T, fol_formula*& definition) const
	{
		if (logical_form->type == fol_formula_type::ATOM) {
			if (logical_form->atom.predicate == PREDICATE_TYPE
			 && logical_form->atom.arg1.type == fol_term_type::CONSTANT
			 && logical_form->atom.arg1.constant == PREDICATE_UNKNOWN)
			{
				definition = logical_form;
				definition->reference_count++;
				return true;
			}
		} else if (logical_form->type == fol_formula_type::FOR_ALL) {
			unsigned int variable = logical_form->quantifier.variable;
			if (logical_form->quantifier.operand->type == fol_formula_type::IF_THEN) {
				const fol_formula* left = logical_form->quantifier.operand->binary.left;
				if (left->type == fol_formula_type::ATOM
				 && left->atom.arg1.type == fol_term_type::VARIABLE
				 && left->atom.arg1.variable == variable
				 && left->atom.arg2.type == fol_term_type::CONSTANT
				 && left->atom.arg2.constant == PREDICATE_UNKNOWN)
				{
					definition = logical_form;
					definition->reference_count++;
					return true;
				}
			}
		}
		return false;
	}

	bool add_definition(const sentence& s, const fol_formula* definition)
	{
		for (unsigned int i = 0; i < s.length; i++) {
			unsigned int index = unknown_tokens.index_of(s.tokens[i]);
			if (index < unknown_tokens.length)
				unknown_tokens.remove(index);
		}

		/* TODO: check if `definition` has any constants in `unknown_concepts`; if so, remove them from `unknown_concepts` */
		constant_subtracter subtracter = {unknown_concepts};
		return visit(*definition, subtracter);
	}
};

#endif /* FAKE_PARSER_H_ */
