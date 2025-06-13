#ifndef LF_UTILS_H_
#define LF_UTILS_H_

#include "higher_order_logic.h"

enum class head_position : uint_fast8_t {
	NONE = 0,
	LEFT,
	RIGHT,
	ANY,
	ALL
};

struct head_index {
	head_position position;
	unsigned int index;

	inline bool compare(int conjunct_index) const {
		if (position == head_position::LEFT && conjunct_index >= 0) {
			return index == (unsigned int) conjunct_index;
		} else if (position == head_position::RIGHT && conjunct_index < 0) {
			return index == (unsigned int) (-conjunct_index) - 1;
		}
		return false;
	}
};

template<typename BuiltInPredicates>
struct hol_tense_constants {
	hol_term* constants;

	hol_tense_constants() {
		constants = hol_term::new_any_constant(make_array_view(TENSE_PREDICATES, array_length(TENSE_PREDICATES)));
	}

	~hol_tense_constants() {
		free(*constants); if (constants->reference_count == 0) free(constants);
	}

	static inline hol_tense_constants<BuiltInPredicates>& get() {
		static thread_local hol_tense_constants<BuiltInPredicates> constants;
		return constants;
	}
};

template<typename BuiltInPredicates>
struct hol_non_head_constants {
	static constexpr unsigned int CONSTANTS[]  =  {
			(unsigned int) BuiltInPredicates::UNKNOWN,
			(unsigned int) BuiltInPredicates::ARG1,
			(unsigned int) BuiltInPredicates::ARG2,
			(unsigned int) BuiltInPredicates::ARG3,
			(unsigned int) BuiltInPredicates::ARG1_OF,
			(unsigned int) BuiltInPredicates::ARG2_OF,
			(unsigned int) BuiltInPredicates::ARG3_OF,
			(unsigned int) BuiltInPredicates::SIZE,
			(unsigned int) BuiltInPredicates::PRESENT,
			(unsigned int) BuiltInPredicates::PRESENT_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::PRESENT_PERFECT,
			(unsigned int) BuiltInPredicates::PRESENT_PERFECT_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::PAST,
			(unsigned int) BuiltInPredicates::PAST_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::PAST_PERFECT,
			(unsigned int) BuiltInPredicates::PAST_PERFECT_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::FUTURE,
			(unsigned int) BuiltInPredicates::FUTURE_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::FUTURE_PERFECT,
			(unsigned int) BuiltInPredicates::FUTURE_PERFECT_PROGRESSIVE,
			(unsigned int) BuiltInPredicates::WIDE_SCOPE
		};
	static constexpr unsigned int CONSTANT_COUNT = (unsigned int) array_length(CONSTANTS);

	hol_term* constants;

	hol_non_head_constants() {
		constants = hol_term::new_any(hol_term::new_any_constant(make_array_view(CONSTANTS, CONSTANT_COUNT)));
	}

	~hol_non_head_constants() {
		free(*constants);
		if (constants->reference_count == 0)
			free(constants);
	}

	static inline hol_non_head_constants<BuiltInPredicates>& get() {
		static thread_local hol_non_head_constants<BuiltInPredicates> constants;
		return constants;
	}

	static inline hol_term** get_terms() {
		return &get().constants;
	}

	static inline void increment_terms() {
		get().constants->reference_count++;
	}

	static inline constexpr unsigned int count() {
		return 1;
	}
};

template<typename BuiltInPredicates>
constexpr unsigned int hol_non_head_constants<BuiltInPredicates>::CONSTANTS[];

template<typename BuiltInPredicates>
constexpr unsigned int hol_non_head_constants<BuiltInPredicates>::CONSTANT_COUNT;

template<typename BuiltInPredicates>
inline bool is_built_in(unsigned int constant) {
	return index_of(constant, hol_non_head_constants<BuiltInPredicates>::CONSTANTS, hol_non_head_constants<BuiltInPredicates>::CONSTANT_COUNT) < hol_non_head_constants<BuiltInPredicates>::CONSTANT_COUNT;
}

template<typename BuiltInPredicates>
inline hol_term* get_predicate_of_literal(
		hol_term* src, unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.right->type == hol_term_type::VARIABLE
	 && src->binary.right->variable == scope_variable
	 && src->binary.left->type == hol_term_type::CONSTANT
	 && !is_built_in<BuiltInPredicates>(src->binary.left->constant))
	{
		return src;
	}
	return nullptr;
}

inline hol_term* get_variable_definition_of_literal(
		hol_term* src, unsigned int scope_variable)
{
	if (src->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.left->type == hol_term_type::VARIABLE
	 && src->binary.left->variable == scope_variable)
	{
		return src;
	}
	return nullptr;
}

template<typename BuiltInPredicates>
inline hol_term* find_predicate(unsigned int head_variable, hol_term* conjunction, head_index& predicate_index)
{
	predicate_index.position = head_position::NONE;
	hol_term* expected_predicate = hol_term::new_apply(
			hol_term::new_any(nullptr, hol_non_head_constants<BuiltInPredicates>::get_terms(), hol_non_head_constants<BuiltInPredicates>::count()),
			hol_term::new_variable(head_variable));
	if (expected_predicate == nullptr)
		return nullptr;
	hol_non_head_constants<built_in_predicates>::increment_terms();

	hol_term* predicate = nullptr;
	head_index first_predicate_index = predicate_index;
	if (conjunction->type == hol_term_type::ANY_ARRAY) {
		for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.any.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::ANY, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.any.operands[i], expected_predicate) && first_predicate_index.position == head_position::NONE) {
				first_predicate_index = {head_position::ANY, i};
			}
		} for (unsigned int i = 0; predicate == nullptr && i < conjunction->any_array.left.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.left.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::LEFT, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.left.operands[i], expected_predicate) && first_predicate_index.position == head_position::NONE) {
				first_predicate_index = {head_position::LEFT, i};
			}
		} for (unsigned int i = 0; predicate == nullptr && i < conjunction->any_array.right.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->any_array.right.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.right.operands[i], expected_predicate) && first_predicate_index.position == head_position::NONE) {
				first_predicate_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
			}
		} if (predicate == nullptr && has_intersection<BuiltInPredicates>(conjunction->any_array.all, expected_predicate) && first_predicate_index.position == head_position::NONE) {
			first_predicate_index = {head_position::ALL, 0};
		}
	} else if (conjunction->type == hol_term_type::AND) {
		predicate_index = {head_position::NONE, 0};
		for (unsigned int i = 0; i < conjunction->array.length; i++) {
			predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction->array.operands[i], head_variable);
			if (predicate != nullptr) {
				predicate_index = {head_position::LEFT, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->array.operands[i], expected_predicate) && first_predicate_index.position == head_position::NONE) {
				first_predicate_index = {head_position::LEFT, i};
			}
		}
	} else {
		/* this conjunct could be singleton */
		predicate = get_predicate_of_literal<BuiltInPredicates>(conjunction, head_variable);
		if (predicate != nullptr) {
			predicate_index = {head_position::LEFT, 0};
		} else if (has_intersection<BuiltInPredicates>(conjunction, expected_predicate) && first_predicate_index.position == head_position::NONE) {
			first_predicate_index = {head_position::LEFT, 0};
		}
	}
	free(*expected_predicate); free(expected_predicate);
	if (predicate == nullptr)
		predicate_index = first_predicate_index;
	return predicate;
}

template<typename BuiltInPredicates>
inline hol_term* find_variable_definition(unsigned int head_variable, hol_term* conjunction, head_index& definition_index)
{
	hol_term* expected_definition = hol_term::new_equals(
			hol_term::new_variable(head_variable), &HOL_ANY);
	if (expected_definition == nullptr)
		return nullptr;
	HOL_ANY.reference_count++;

	hol_term* variable_definition = nullptr;
	definition_index.position = head_position::NONE;
	if (conjunction->type == hol_term_type::ANY_ARRAY) {
		for (unsigned int i = 0; i < conjunction->any_array.left.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.left.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::LEFT, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.left.operands[i], expected_definition)) {
				definition_index = {head_position::LEFT, i};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.right.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.right.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.right.operands[i], expected_definition)) {
				definition_index = {head_position::RIGHT, conjunction->any_array.right.length - i - 1};
			}
		} for (unsigned int i = 0; i < conjunction->any_array.any.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->any_array.any.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::ANY, i};
				break;
			} else if (has_intersection<BuiltInPredicates>(conjunction->any_array.any.operands[i], expected_definition)) {
				definition_index = {head_position::ANY, i};
			}
		}
	} else if (conjunction->type == hol_term_type::AND) {
		definition_index = {head_position::NONE, 0};
		for (unsigned int i = 0; i < conjunction->array.length; i++) {
			variable_definition = get_variable_definition_of_literal(conjunction->array.operands[i], head_variable);
			if (variable_definition != nullptr) {
				definition_index = {head_position::LEFT, i};
				break;
			} else if (variable_definition == nullptr && has_intersection<BuiltInPredicates>(conjunction->array.operands[i], expected_definition)) {
				definition_index = {head_position::LEFT, i};
			}
		}
	} else {
		/* this conjunct could be singleton */
		variable_definition = get_variable_definition_of_literal(conjunction, head_variable);
		if (variable_definition != nullptr) {
			definition_index = {head_position::LEFT, 0};
		} else if (variable_definition == nullptr && has_intersection<BuiltInPredicates>(conjunction, expected_definition)) {
			definition_index = {head_position::LEFT, 0};
		}
	}
	free(*expected_definition); free(expected_definition);
	return variable_definition;
}

template<typename BuiltInPredicates, bool AllowVariableDefinitions = false, bool IncludeNegation = true, bool IncludeAny = true>
inline void find_head(
		hol_term* src, hol_term*& head,
		head_index& predicate_index)
{
	head = nullptr;
	if (IncludeAny && (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT)) {
		if (src->any.included == nullptr) {
			head = src;
		} else {
			find_head<BuiltInPredicates>(src->any.included, head, predicate_index);
			if (head != nullptr)
				head = src;
		}
	} else if (src->type == hol_term_type::EXISTS) {
		head = src;
		unsigned int head_variable = src->quantifier.variable;

		/* make sure this scope has no term `arg1(*)=x` or `arg2(*)=x` where `x` is the scope variable */
		hol_term wildcard;
		wildcard.any.included = hol_term::new_equals(hol_term::new_apply(
				hol_term::new_any_constant(
					(unsigned int) BuiltInPredicates::ARG1, (unsigned int) BuiltInPredicates::ARG2, (unsigned int) BuiltInPredicates::ARG3,
					(unsigned int) BuiltInPredicates::ARG1_OF, (unsigned int) BuiltInPredicates::ARG2_OF, (unsigned int) BuiltInPredicates::ARG3_OF),
				&HOL_ANY), hol_term::new_variable(head_variable));
		if (wildcard.any.included == nullptr) return;
		HOL_ANY.reference_count++;
		bool not_head = has_intersection<BuiltInPredicates>(head->quantifier.operand, &wildcard);
		free(*wildcard.any.included); free(wildcard.any.included);
		wildcard.any.included = nullptr;
		if (not_head) { head = nullptr; return; }

		/* find the predicate */
		find_predicate<BuiltInPredicates>(head_variable, head->quantifier.operand, predicate_index);
		if (predicate_index.position == head_position::NONE) {
			if (AllowVariableDefinitions)
				find_variable_definition<BuiltInPredicates>(head_variable, head->quantifier.operand, predicate_index);

			if (!AllowVariableDefinitions || predicate_index.position == head_position::NONE) {
				head = nullptr;
				return;
			}
		}
	} else if (IncludeNegation && src->type == hol_term_type::NOT) {
		find_head<BuiltInPredicates, AllowVariableDefinitions, IncludeNegation, false>(src->unary.operand, head, predicate_index);
		if (head != nullptr)
			head = src;
	}
}

template<typename TryFindHeadFunction, typename Function>
hol_term* find_head(hol_term* term, head_index& predicate_index, TryFindHeadFunction& try_find_head, Function& apply)
{
	apply(term);

	hol_term* head;
	try_find_head(term, head, predicate_index);
	if (head != nullptr) return head;

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`,
	   `apply_arg`, and `find_head` in `hdp_parser.h` */
	switch (term->type) {
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::NUMBER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::BINARY_APPLICATION:
	case hol_term_type::IFF:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
		return nullptr;

	case hol_term_type::NOT:
		return find_head(term->unary.operand, predicate_index, try_find_head, apply);

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		return find_head(term->quantifier.operand, predicate_index, try_find_head, apply);

	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		return find_head(term->binary.right, predicate_index, try_find_head, apply);

	case hol_term_type::AND:
	case hol_term_type::OR:
		return find_head(term->array.operands[term->array.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_RIGHT_ONLY:
		if (term->any.included == nullptr) return nullptr;
		return find_head(term->any.included, predicate_index, try_find_head, apply);

	case hol_term_type::ANY_ARRAY:
		if (term->any_array.right.length == 0)
			return find_head(term->any_array.all, predicate_index, try_find_head, apply);
		return find_head(term->any_array.right.operands[term->any_array.right.length - 1], predicate_index, try_find_head, apply);

	case hol_term_type::ANY_QUANTIFIER:
		return find_head(term->any_quantifier.operand, predicate_index, try_find_head, apply);
	}
	fprintf(stderr, "find_head ERROR: Unrecognied hol_term_type.\n");
	return nullptr;
}

enum class any_node_position {
	LEFT,
	RIGHT,
	NONE
};

template<any_node_position AnyNodePosition, bool AnyRightOnly>
struct head_substituter {
	bool found_wide_scope;
	bool could_have_wide_scope;
	bool first_any;
	bool has_any;
	bool parent_is_negation;
	const hol_term* src;
	hol_term* dst;
	hol_term* last_quantifier;
	unsigned int last_declared_variable;

	head_substituter(bool found_wide_scope, bool could_have_wide_scope, bool first_any, const hol_term* src, hol_term* dst) :
			found_wide_scope(found_wide_scope), could_have_wide_scope(could_have_wide_scope), first_any(first_any), has_any(false), parent_is_negation(false), src(src), dst(dst), last_quantifier(nullptr), last_declared_variable(0) { }
};

template<bool AnyRightOnly>
inline hol_term* wrap_any_right(hol_term* operand, unsigned int prev_declared_variable, bool exclude_wide_scope, hol_term* last_quantifier)
{
	if (prev_declared_variable != 0) {
		unsigned int excluded_tree_count = 3 + ((exclude_wide_scope || last_quantifier != nullptr) ? 1 : 0);
		hol_term* excluded_trees[4];
		excluded_trees[0] = hol_term::new_any(hol_term::new_for_all(prev_declared_variable, &HOL_ANY));
		excluded_trees[1] = hol_term::new_any(hol_term::new_exists(prev_declared_variable, &HOL_ANY));
		excluded_trees[2] = hol_term::new_any(hol_term::new_lambda(prev_declared_variable, &HOL_ANY));
		if (exclude_wide_scope)
			excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, &HOL_ANY));
		else if (last_quantifier != nullptr) {
			if (last_quantifier->type == hol_term_type::FOR_ALL) {
				excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_for_all(last_quantifier->quantifier.variable, &HOL_ANY))));
			} else if (last_quantifier->type == hol_term_type::EXISTS) {
				excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_exists(last_quantifier->quantifier.variable, &HOL_ANY))));
			} else {
				excluded_trees[3] = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_lambda(last_quantifier->quantifier.variable, &HOL_ANY))));
			}
		}
		if (excluded_trees[0] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[1] != nullptr) HOL_ANY.reference_count++;
		if (excluded_trees[2] != nullptr) HOL_ANY.reference_count++;
		if ((exclude_wide_scope || last_quantifier != nullptr) && excluded_trees[3] != nullptr) {
			hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value.reference_count++;
			HOL_ANY.reference_count++;
		}
		if (excluded_trees[0] == nullptr || excluded_trees[1] == nullptr || excluded_trees[2] == nullptr || ((exclude_wide_scope || last_quantifier != nullptr) && excluded_trees[3] == nullptr)) {
			if (excluded_trees[0] != nullptr) { free(*excluded_trees[0]); free(excluded_trees[0]); }
			if (excluded_trees[1] != nullptr) { free(*excluded_trees[1]); free(excluded_trees[1]); }
			if (excluded_trees[2] != nullptr) { free(*excluded_trees[2]); free(excluded_trees[2]); }
			return nullptr;
		}

		hol_term* term;
		if (AnyRightOnly)
			term = hol_term::new_any_right_only(operand, excluded_trees, excluded_tree_count);
		else term = hol_term::new_any_right(operand, excluded_trees, excluded_tree_count);
		if (term == nullptr) {
			for (unsigned int i = 0; i < excluded_tree_count; i++) {
				free(*excluded_trees[i]); free(excluded_trees[i]);
			}
			return nullptr;
		}
		return term;
	} else if (exclude_wide_scope) {
		hol_term* excluded = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, &HOL_ANY));
		if (excluded == nullptr) return nullptr;
		hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value.reference_count++;
		HOL_ANY.reference_count++;
		hol_term* term;
		if (AnyRightOnly)
			term = hol_term::new_any_right_only(operand, &excluded, 1);
		else term = hol_term::new_any_right(operand, &excluded, 1);
		if (term == nullptr) {
			free(*excluded); free(excluded);
			return nullptr;
		}
		return term;
	} else if (last_quantifier != nullptr) {
		hol_term* excluded;
		if (last_quantifier->type == hol_term_type::FOR_ALL)
			excluded = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_for_all(last_quantifier->quantifier.variable, &HOL_ANY))));
		else if (last_quantifier->type == hol_term_type::EXISTS)
			excluded = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_exists(last_quantifier->quantifier.variable, &HOL_ANY))));
		else
			excluded = hol_term::new_any_right(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value, hol_term::new_any_right(hol_term::new_lambda(last_quantifier->quantifier.variable, &HOL_ANY))));
		if (excluded == nullptr) return nullptr;
		hol_term::constants<(unsigned int) built_in_predicates::WIDE_SCOPE>::value.reference_count++;
		HOL_ANY.reference_count++;
		hol_term* term;
		if (AnyRightOnly)
			term = hol_term::new_any_right_only(operand, &excluded, 1);
		else term = hol_term::new_any_right(operand, &excluded, 1);
		if (term == nullptr) {
			free(*excluded); free(excluded);
			return nullptr;
		}
		return term;
	} else {
		if (AnyRightOnly)
			return hol_term::new_any_right_only(operand);
		else return hol_term::new_any_right(operand);
	}
}

template<hol_term_type Type, any_node_position AnyNodePosition, bool AnyRightOnly>
inline hol_term* apply(hol_term* src, head_substituter<AnyNodePosition, AnyRightOnly>& substituter)
{
	if (src == substituter.src) {
		if (AnyNodePosition == any_node_position::LEFT) {
			if (substituter.dst->type == hol_term_type::UNARY_APPLICATION && substituter.dst->binary.left->type == hol_term_type::CONSTANT && substituter.dst->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE)
				substituter.could_have_wide_scope = true;
			if (substituter.dst->type == hol_term_type::ANY || substituter.dst->type == hol_term_type::ANY_RIGHT || substituter.dst->type == hol_term_type::ANY_RIGHT_ONLY) {
				/* check if wide scope is excluded */
				bool wide_scope_excluded = false;
				for (unsigned int i = 0; i < substituter.dst->any.excluded_tree_count; i++) {
					hol_term* tree = substituter.dst->any.excluded_trees[i];
					if (tree->type == hol_term_type::ANY_RIGHT && tree->any.included != nullptr
						&& tree->any.included->type == hol_term_type::UNARY_APPLICATION
						&& tree->any.included->binary.left->type == hol_term_type::CONSTANT
						&& tree->any.included->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE
						&& *tree->any.included->binary.right == HOL_ANY)
					{
						wide_scope_excluded = true;
						break;
					}
				}

				if (!wide_scope_excluded)
					substituter.could_have_wide_scope = true;
			}
			if (substituter.has_any || substituter.dst->type == hol_term_type::ANY || substituter.dst->type == hol_term_type::ANY_RIGHT) {
				if (substituter.src != substituter.dst)
					substituter.dst->reference_count++;
				if (substituter.dst->type == hol_term_type::EXISTS || substituter.dst->type == hol_term_type::FOR_ALL || substituter.dst->type == hol_term_type::LAMBDA)
					substituter.last_quantifier = substituter.dst;
				return substituter.dst;
			} else {
				hol_term* new_term = wrap_any_right<AnyRightOnly>(substituter.dst, substituter.last_declared_variable, substituter.found_wide_scope || (!substituter.could_have_wide_scope && substituter.first_any), substituter.last_quantifier);
				substituter.first_any = false;
				if (new_term == nullptr)
					return nullptr;
				substituter.dst->reference_count++;
				if (substituter.dst->type == hol_term_type::EXISTS || substituter.dst->type == hol_term_type::FOR_ALL || substituter.dst->type == hol_term_type::LAMBDA)
					substituter.last_quantifier = substituter.dst;
				return new_term;
			}
		} else if (AnyNodePosition == any_node_position::RIGHT && (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT)) {
			hol_term* new_term = wrap_any_right<false>(substituter.dst, substituter.last_declared_variable, substituter.found_wide_scope || (!substituter.could_have_wide_scope && substituter.first_any), substituter.last_quantifier);
			substituter.first_any = false;
			if (new_term == nullptr)
				return nullptr;
			substituter.dst->reference_count++;
			if (substituter.dst->type == hol_term_type::EXISTS || substituter.dst->type == hol_term_type::FOR_ALL || substituter.dst->type == hol_term_type::LAMBDA)
				substituter.last_quantifier = substituter.dst;

			array<hol_term*> excluded(max(1u, new_term->any.excluded_tree_count + src->any.excluded_tree_count));
			for (unsigned int i = 0; i < new_term->any.excluded_tree_count; i++) {
				/* make sure this tree isn't a subset of any tree in `src->any.excluded_trees` */
				bool irreducible = true;
				for (unsigned int j = 0; j < src->any.excluded_tree_count; j++) {
					if (is_subset<built_in_predicates>(new_term->any.excluded_trees[i], src->any.excluded_trees[j])) {
						irreducible = false;
						break;
					}
				}
				if (irreducible)
					excluded[excluded.length++] = new_term->any.excluded_trees[i];
			}
			unsigned int old_excluded_tree_count = excluded.length;
			for (unsigned int i = 0; i < src->any.excluded_tree_count; i++) {
				/* make sure this tree isn't a subset of any tree in `excluded` so far */
				bool irreducible = true;
				for (unsigned int j = 0; j < old_excluded_tree_count; j++) {
					if (is_subset<built_in_predicates>(src->any.excluded_trees[i], excluded[j])) {
						irreducible = false;
						break;
					}
				}
				if (irreducible)
					excluded[excluded.length++] = src->any.excluded_trees[i];
			}

			hol_term* temp;
			if (new_term->type == hol_term_type::ANY_RIGHT) {
				temp = hol_term::new_any_right(new_term->any.included, excluded.data, excluded.length);
			} else if (new_term->type == hol_term_type::ANY_RIGHT_ONLY) {
				temp = hol_term::new_any_right_only(new_term->any.included, excluded.data, excluded.length);
			} else {
				temp = hol_term::new_any(new_term->any.included, excluded.data, excluded.length);
			}
			if (temp == nullptr) {
				free(*new_term); if (new_term->reference_count == 0) free(new_term);
				return nullptr;
			}
			new_term->any.included->reference_count++;
			for (hol_term* tree : excluded)
				tree->reference_count++;
			free(*new_term); if (new_term->reference_count == 0) free(new_term);
			return temp;
		}
		if (substituter.src != substituter.dst)
			substituter.dst->reference_count++;
		if (substituter.dst->type == hol_term_type::EXISTS || substituter.dst->type == hol_term_type::FOR_ALL || substituter.dst->type == hol_term_type::LAMBDA)
			substituter.last_quantifier = substituter.dst;
		return substituter.dst;
	}

	/* NOTE: this function should mirror the semantics of
	   `apply_head_conjunct`, `apply<head_substituter>`, `apply_arg`, and
	   `find_head` in `hdp_parser.h` */
	bool changed; hol_term* new_term = nullptr;
	hol_term* first; hol_term* second; hol_term* third;
	unsigned int prev_declared_variable = 0;
	bool found_wide_scope = false, first_any = false, has_any = false, parent_is_negation = false;
	switch (Type) {
	case hol_term_type::IF_THEN:
	case hol_term_type::EQUALS:
	case hol_term_type::UNARY_APPLICATION:
		if ((AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && src->binary.right->type != hol_term_type::ANY && src->binary.right->type != hol_term_type::ANY_RIGHT)
		 || (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::UNARY_APPLICATION))
		{
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.found_wide_scope = (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE);
		substituter.parent_is_negation = false;
		first = src->binary.left;
		second = apply(src->binary.right, substituter);
		if (second == nullptr)
			return nullptr;

		if (src->type == hol_term_type::UNARY_APPLICATION && src->binary.left->type == hol_term_type::CONSTANT && src->binary.left->constant == (unsigned int) built_in_predicates::WIDE_SCOPE) {
			substituter.could_have_wide_scope = true;
			substituter.last_quantifier = nullptr;
		}

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::IF_THEN && second->type != hol_term_type::ANY && second->type != hol_term_type::ANY_RIGHT && second->type != hol_term_type::ANY_RIGHT_ONLY && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(second, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (second != src->binary.right) { free(*second); if (second->reference_count == 0) free(second); }
				return nullptr;
			}
			if (second == src->binary.right) second->reference_count++;
			second = term;
		}

		if (second == src->binary.right) {
			new_term = src;
		} else {
			if (!new_hol_term(new_term)) {
				if (second != src->binary.right) {
					free(*second); if (second->reference_count == 0) free(second);
				}
				return nullptr;
			}
			new_term->binary.left = first;
			new_term->binary.right = second;
			first->reference_count++;
			if (second == src->binary.right) second->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
		}

		if (AnyNodePosition == any_node_position::LEFT && src->type == hol_term_type::UNARY_APPLICATION && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::BINARY_APPLICATION:
		if (AnyNodePosition == any_node_position::RIGHT || AnyNodePosition == any_node_position::LEFT) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
			substituter.parent_is_negation = false;
		}
		first = src->ternary.first;
		second = src->ternary.second;
		third = apply(src->ternary.third, substituter);
		if (third == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT && third->type != hol_term_type::ANY && third->type != hol_term_type::ANY_RIGHT && third->type != hol_term_type::ANY_RIGHT_ONLY) {
			hol_term* term = wrap_any_right<AnyRightOnly>(third, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (third != src->unary.operand) { free(*third); if (third->reference_count == 0) free(third); }
				return nullptr;
			}
			if (third == src->unary.operand) third->reference_count++;
			third = term;
		}

		if (third == src->ternary.third) {
			new_term = src;
		} else {
			if (!new_hol_term(new_term)) {
				if (third != src->ternary.third) {
					free(*third); if (third->reference_count == 0) free(third);
				}
				return nullptr;
			}
			new_term->ternary.first = first;
			new_term->ternary.second = second;
			new_term->ternary.third = third;
			first->reference_count++;
			second->reference_count++;
			if (third == src->ternary.third) third->reference_count++;
			new_term->type = Type;
			new_term->reference_count = 1;
		}

		if (AnyNodePosition == any_node_position::LEFT && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		changed = false;
		substituter.has_any = false;
		substituter.parent_is_negation = false;
		first = apply(src->array.operands[src->array.length - 1], substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != src->array.operands[src->array.length - 1])
			changed = true;

		if (!changed) {
			return src;
		} else {
			if (!new_hol_term(new_term)) {
				if (first != src->array.operands[src->array.length - 1]) {
					free(*first);
					if (first->reference_count == 0)
						free(first);
				}
				return nullptr;
			}
			new_term->type = Type;
			new_term->reference_count = 1;
			if (first->type == Type && first->reference_count == 1) {
				/* the new child is the same type as the parent, so merge them */
				new_term->array.operands = (hol_term**) realloc(first->array.operands, sizeof(hol_term*) * (src->array.length + first->array.length - 1));
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first); if (first->reference_count == 0) free(first);
					}
					return nullptr;
				}
				for (unsigned int i = first->array.length; i > 0; i--)
					new_term->array.operands[i + src->array.length - 2] = new_term->array.operands[i - 1];
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.length = src->array.length + first->array.length - 1;
				free(first);
			} else {
				new_term->array.length = src->array.length;
				new_term->array.operands = (hol_term**) malloc(sizeof(hol_term*) * new_term->array.length);
				if (new_term->array.operands == nullptr) {
					if (first != src->array.operands[src->array.length - 1]) {
						free(*first);
						if (first->reference_count == 0)
							free(first);
					}
					free(new_term);
					return nullptr;
				}
				for (unsigned int i = 0; i + 1 < src->array.length; i++) {
					new_term->array.operands[i] = src->array.operands[i];
					new_term->array.operands[i]->reference_count++;
				}
				new_term->array.operands[src->array.length - 1] = first;
				if (first == src->array.operands[src->array.length - 1])
					first->reference_count++;
			}
			return new_term;
		}

	case hol_term_type::ANY_ARRAY:
		changed = false;
		substituter.has_any = false;
		substituter.parent_is_negation = false;
		if (src->any_array.right.length == 0)
			first = apply(src->any_array.all, substituter);
		else first = apply(src->any_array.right.operands[src->any_array.right.length - 1], substituter);
		if (first == nullptr) {
			return nullptr;
		} else if (first != (src->any_array.right.length == 0 ? src->any_array.all : src->any_array.right.operands[src->any_array.right.length - 1])) {
			changed = true;
		}

		if (!changed) {
			new_term = src;
		} else {
			if (src->any_array.right.length == 0) {
				new_term = hol_term::new_any_array(src->any_array.oper, first,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_array_view(src->any_array.right.operands, src->any_array.right.length));
			} else {
				new_term = hol_term::new_any_array(src->any_array.oper, src->any_array.all,
						make_array_view(src->any_array.any.operands, src->any_array.any.length),
						make_array_view(src->any_array.left.operands, src->any_array.left.length),
						make_appended_array_view(make_array_view(src->any_array.right.operands, src->any_array.right.length - 1), first));
			}
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
			new_term->any_array.all->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.any.length; i++)
				new_term->any_array.any.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.left.length; i++)
				new_term->any_array.left.operands[i]->reference_count++;
			for (unsigned int i = 0; i < new_term->any_array.right.length; i++)
				new_term->any_array.right.operands[i]->reference_count++;
			first->reference_count--;
		}
		return new_term;

	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if ((AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA)
		 || (AnyNodePosition == any_node_position::LEFT && src->type != hol_term_type::LAMBDA))
		{
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.last_declared_variable = src->quantifier.variable;
		substituter.parent_is_negation = false;
		first = apply(src->quantifier.operand, substituter);
		if (first == nullptr)
			return nullptr;
		if (substituter.last_quantifier == nullptr && (first->type == hol_term_type::FOR_ALL || first->type == hol_term_type::EXISTS || first->type == hol_term_type::LAMBDA))
			substituter.last_quantifier = first;

		if (AnyNodePosition == any_node_position::RIGHT && src->type == hol_term_type::LAMBDA && first->type != hol_term_type::ANY && first->type != hol_term_type::ANY_RIGHT && first->type != hol_term_type::ANY_RIGHT_ONLY) {
			hol_term* term = wrap_any_right<AnyRightOnly>(first, src->quantifier.variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (first != src->quantifier.operand) { free(*first); if (first->reference_count == 0) free(first); }
				return nullptr;
			}
			if (first == src->quantifier.operand) first->reference_count++;
			first = term;
		}

		if (first == src->quantifier.operand) {
			new_term = src;
		} else {
			if (src->type == hol_term_type::FOR_ALL)
				new_term = hol_term::new_for_all(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::EXISTS)
				new_term = hol_term::new_exists(src->quantifier.variable_type, src->quantifier.variable, first);
			else if (src->type == hol_term_type::LAMBDA)
				new_term = hol_term::new_lambda(src->quantifier.variable_type, src->quantifier.variable, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && src->type != hol_term_type::LAMBDA && !parent_is_negation && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::ANY_QUANTIFIER:
		if (AnyNodePosition == any_node_position::LEFT && src->any_quantifier.quantifier != hol_quantifier_type::EXISTS) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.parent_is_negation = false;
		first = apply(src->any_quantifier.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (first == src->any_quantifier.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_any_quantifier(src->any_quantifier.quantifier, first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && src->any_quantifier.quantifier != hol_quantifier_type::EXISTS && !parent_is_negation && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::NOT:
		if (AnyNodePosition == any_node_position::RIGHT || AnyNodePosition == any_node_position::LEFT) {
			prev_declared_variable = substituter.last_declared_variable;
			found_wide_scope = substituter.found_wide_scope;
			first_any = substituter.first_any;
			has_any = substituter.has_any;
			parent_is_negation = substituter.parent_is_negation;
			substituter.last_declared_variable = 0;
			substituter.found_wide_scope = false;
			substituter.first_any = false;
			substituter.has_any = false;
		}
		substituter.parent_is_negation = true;
		first = apply(src->unary.operand, substituter);
		if (first == nullptr)
			return nullptr;

		if (AnyNodePosition == any_node_position::RIGHT && first->type != hol_term_type::ANY && first->type != hol_term_type::ANY_RIGHT && first->type != hol_term_type::ANY_RIGHT_ONLY) {
			hol_term* term = wrap_any_right<AnyRightOnly>(first, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (first != src->unary.operand) { free(*first); if (first->reference_count == 0) free(first); }
				return nullptr;
			}
			if (first == src->unary.operand) first->reference_count++;
			first = term;
		}

		if (first == src->unary.operand) {
			new_term = src;
		} else {
			new_term = hol_term::new_not(first);
			if (new_term == nullptr) {
				free(*first); if (first->reference_count == 0) free(first);
				return nullptr;
			}
		}

		if (AnyNodePosition == any_node_position::LEFT && !has_any) {
			hol_term* term = wrap_any_right<AnyRightOnly>(new_term, prev_declared_variable, found_wide_scope || (first_any && !substituter.could_have_wide_scope), substituter.last_quantifier);
			if (term == nullptr) {
				if (new_term != src) { free(*new_term); if (new_term->reference_count == 0) free(new_term); }
				return nullptr;
			}
			if (new_term == src) new_term->reference_count++;
			new_term = term;
		}
		return new_term;

	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_RIGHT_ONLY:
		if (src->any.included == nullptr)
			return nullptr;
		substituter.has_any = true;
		substituter.parent_is_negation = false;
		first = apply(src->any.included, substituter);
		if (first == nullptr)
			return nullptr;
		substituter.could_have_wide_scope = true;
		substituter.last_quantifier = nullptr;

		if (first == src->any.included) {
			new_term = src;
		} else {
			if (first->type == src->type) {
				array<hol_term*> excluded(max(1u, first->any.excluded_tree_count + src->any.excluded_tree_count));
				for (unsigned int i = 0; i < first->any.excluded_tree_count; i++) {
					/* make sure this tree isn't a subset of any tree in `src->any.excluded_trees` */
					bool irreducible = true;
					for (unsigned int j = 0; j < src->any.excluded_tree_count; j++) {
						if (is_subset<built_in_predicates>(first->any.excluded_trees[i], src->any.excluded_trees[j])) {
							irreducible = false;
							break;
						}
					}
					if (irreducible)
						excluded[excluded.length++] = first->any.excluded_trees[i];
				}
				unsigned int old_excluded_tree_count = excluded.length;
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++) {
					/* make sure this tree isn't a subset of any tree in `excluded` so far */
					bool irreducible = true;
					for (unsigned int j = 0; j < old_excluded_tree_count; j++) {
						if (is_subset<built_in_predicates>(src->any.excluded_trees[i], excluded[j])) {
							irreducible = false;
							break;
						}
					}
					if (irreducible)
						excluded[excluded.length++] = src->any.excluded_trees[i];
				}
				if (src->type == hol_term_type::ANY)
					new_term = hol_term::new_any(first->any.included, excluded.data, excluded.length);
				else new_term = hol_term::new_any_right(first->any.included, excluded.data, excluded.length);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				if (first->any.included != nullptr)
					first->any.included->reference_count++;
				for (hol_term* tree : excluded)
					tree->reference_count++;
				free(*first); if (first->reference_count == 0) free(first);
			} else if (src->type == hol_term_type::ANY || src->type == hol_term_type::ANY_RIGHT) {
				if (src->type == hol_term_type::ANY)
					new_term = hol_term::new_any(first, src->any.excluded_trees, src->any.excluded_tree_count);
				else new_term = hol_term::new_any_right(first, src->any.excluded_trees, src->any.excluded_tree_count);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
					src->any.excluded_trees[i]->reference_count++;
			} else {
				new_term = hol_term::new_any_right_only(first, src->any.excluded_trees, src->any.excluded_tree_count);
				if (new_term == nullptr) {
					free(*first); if (first->reference_count == 0) free(first);
					return nullptr;
				}
				for (unsigned int i = 0; i < src->any.excluded_tree_count; i++)
					src->any.excluded_trees[i]->reference_count++;
			}
		}
		return new_term;

	case hol_term_type::CONSTANT:
	case hol_term_type::VARIABLE:
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::PARAMETER:
	case hol_term_type::NUMBER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
		return default_apply<Type>(src, substituter);
	}
	fprintf(stderr, "apply ERROR: Unrecognized hol_term_type when substituting head.\n");
	return NULL;
}

template<any_node_position AnyNodePosition, bool AnyRightOnly = true>
inline hol_term* substitute_head(hol_term* src,
		const hol_term* src_term, hol_term* dst_term,
		bool could_have_wide_scope = false)
{
	head_substituter<AnyNodePosition, AnyRightOnly> substituter(false, could_have_wide_scope, true, src_term, dst_term);
	hol_term* dst = apply(src, substituter);
	if (dst == nullptr)
		return nullptr;

	if (AnyNodePosition == any_node_position::RIGHT && dst->type == hol_term_type::ANY_RIGHT_ONLY) {
		hol_term* term = hol_term::new_any_right(dst->any.included, dst->any.excluded_trees, dst->any.excluded_tree_count);
		if (term->any.included != nullptr)
			term->any.included->reference_count++;
		for (unsigned int i = 0; i < term->any.excluded_tree_count; i++)
			term->any.excluded_trees[i]->reference_count++;
		if (dst != src) { free(*dst); if (dst->reference_count == 0) free(dst); }
		if (term == nullptr)
			return nullptr;
		dst = term;
	} else if (AnyNodePosition == any_node_position::RIGHT && dst->type != hol_term_type::LAMBDA && dst->type != hol_term_type::ANY && dst->type != hol_term_type::ANY_RIGHT) {
		hol_term* term = hol_term::new_any_right(dst);
		if (term == nullptr) {
			if (dst != src) { free(*dst); if (dst->reference_count == 0) free(dst); }
			return nullptr;
		}
		if (dst == src) dst->reference_count++;
		dst = term;
	}

	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct ambiguous_to_unknown_converter { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::ANY_CONSTANT_EXCEPT>::type* = nullptr>
inline hol_term* apply(hol_term* src, ambiguous_to_unknown_converter& converter) {
	HOL_UNKNOWN.reference_count++;
	return &HOL_UNKNOWN;
}

inline hol_term* ambiguous_to_unknown(hol_term* src)
{
	ambiguous_to_unknown_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct same_to_equals_converter { };

inline bool process_event_conjunct(
		unsigned int variable, hol_term* conjunct,
		hol_term*& arg1, hol_term*& arg2,
		hol_term*& arg1_of, hol_term*& arg2_of,
		hol_term*& predicate)
{
	if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
			&& conjunct->binary.left->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.left->binary.right->variable == variable)
	{
		arg1 = conjunct->binary.right;
	} else if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
			&& conjunct->binary.left->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.left->binary.right->variable == variable)
	{
		arg2 = conjunct->binary.right;
	} else if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1_OF
			&& conjunct->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.right->variable == variable)
	{
		arg1_of = conjunct->binary.left->binary.right;
	} else if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
			&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2_OF
			&& conjunct->binary.right->type == hol_term_type::VARIABLE
			&& conjunct->binary.right->variable == variable)
	{
		arg2_of = conjunct->binary.left->binary.right;
	} else if (conjunct->type == hol_term_type::UNARY_APPLICATION
			&& conjunct->binary.right->type == hol_term_type::VARIABLE && conjunct->binary.right->variable == variable)
	{
		if (conjunct->binary.left->type == hol_term_type::CONSTANT && is_tense_predicate(conjunct->binary.left->constant))
			return true;
		if (predicate != nullptr) return false;
		predicate = conjunct->binary.left;
	} else if (conjunct->type == hol_term_type::TRUE) {
		return true;
	} else {
		return false;
	}
	return true;
}

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, same_to_equals_converter& converter) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND || (operand->array.length != 3 && operand->array.length != 4))
		return default_apply<Type>(src, converter);

	hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
	hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
	hol_term* predicate = nullptr;
	for (unsigned int i = 0; i < operand->array.length; i++) {
		if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate))
			return default_apply<Type>(src, converter);
	}

	if (arg1 == nullptr || arg2 == nullptr || arg1_of != nullptr || arg2_of != nullptr
	 || predicate->type != hol_term_type::CONSTANT || predicate->constant != (unsigned int) built_in_predicates::SAME)
		return default_apply<Type>(src, converter);

	hol_term* new_term = hol_term::new_equals(arg1, arg2);
	if (new_term == nullptr)
		return nullptr;
	arg1->reference_count++;
	arg2->reference_count++;
	return new_term;
}

inline hol_term* same_to_equals(hol_term* src)
{
	same_to_equals_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct exists_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, exists_remover& converter) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND || operand->array.length != 3)
		return default_apply<Type>(src, converter);

	hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
	hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
	hol_term* predicate = nullptr;
	for (unsigned int i = 0; i < operand->array.length; i++) {
		if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate))
			return default_apply<Type>(src, converter);
	}

	if (arg1 == nullptr || arg2 != nullptr || arg1_of != nullptr || arg2_of != nullptr
	 || predicate->type != hol_term_type::CONSTANT || predicate->constant != (unsigned int) built_in_predicates::EXIST)
		return default_apply<Type>(src, converter);

	HOL_TRUE.reference_count++;
	return &HOL_TRUE;
}

inline hol_term* remove_exists(hol_term* src)
{
	exists_remover converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct arg_of_to_arg_converter { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EQUALS>::type* = nullptr>
inline hol_term* apply(hol_term* src, arg_of_to_arg_converter& converter) {
	if (src->binary.left->type == hol_term_type::UNARY_APPLICATION
	 && src->binary.left->binary.left->type == hol_term_type::CONSTANT)
	{
		unsigned int index = index_of(src->binary.left->binary.left->constant, ARGS_OF, array_length(ARGS_OF));
		if (index < array_length(ARGS_OF)) {
			hol_term* left = apply(src->binary.left->binary.right, converter);
			if (left == nullptr) return nullptr;
			hol_term* right = apply(src->binary.right, converter);
			if (right == nullptr) {
				if (left != src->binary.left->binary.right) { free(*left); if (left->reference_count == 0) free(left); }
				return nullptr;
			}

			hol_term* new_term = hol_term::new_equals(hol_term::new_apply(hol_term::new_constant(ARGS[index]), right), left);
			if (new_term == nullptr) {
				if (left != src->binary.left->binary.right) { free(*left); if (left->reference_count == 0) free(left); }
				if (right != src->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				return nullptr;
			}
			if (left == src->binary.left->binary.right) left->reference_count++;
			if (right == src->binary.right) right->reference_count++;
			return new_term;
		}
	} if (src->binary.right->type == hol_term_type::UNARY_APPLICATION
	   && src->binary.right->binary.left->type == hol_term_type::CONSTANT)
	{
		unsigned int index = index_of(src->binary.right->binary.left->constant, ARGS_OF, array_length(ARGS_OF));
		if (index < array_length(ARGS_OF)) {
			hol_term* right = apply(src->binary.right->binary.right, converter);
			if (right == nullptr) return nullptr;
			hol_term* left = apply(src->binary.left, converter);
			if (left == nullptr) {
				if (right != src->binary.right->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				return nullptr;
			}

			hol_term* new_term = hol_term::new_equals(right, hol_term::new_apply(hol_term::new_constant(ARGS[index]), left));
			if (new_term == nullptr) {
				if (right != src->binary.right->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				if (left != src->binary.left) { free(*left); if (left->reference_count == 0) free(left); }
				return nullptr;
			}
			if (right == src->binary.right->binary.right) right->reference_count++;
			if (left == src->binary.left) left->reference_count++;
			return new_term;
		}
	}
	return default_apply<Type>(src, converter);
}

inline hol_term* arg_of_to_arg(hol_term* src)
{
	arg_of_to_arg_converter converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct tense_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, tense_remover& remover) {
	if (src->binary.left->type == hol_term_type::CONSTANT)
	{
		if (is_tense_predicate(src->binary.left->constant)) {
			HOL_TRUE.reference_count++;
			return &HOL_TRUE;
		}
	}
	return default_apply<Type>(src, remover);
}

inline hol_term* remove_tense(hol_term* src)
{
	tense_remover remover;
	hol_term* dst = apply(src, remover);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct object_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, object_remover& remover) {
	if (src->binary.left->type == hol_term_type::CONSTANT
	 && (src->binary.left->constant == (unsigned int) built_in_predicates::OBJECT
	  || src->binary.left->constant == (unsigned int) built_in_predicates::ANIMATE))
	{
		HOL_TRUE.reference_count++;
		return &HOL_TRUE;
	}
	return default_apply<Type>(src, remover);
}

inline hol_term* remove_object(hol_term* src)
{
	object_remover remover;
	hol_term* dst = apply(src, remover);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct inverse_remover { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, inverse_remover& converter) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND)
		return default_apply<Type>(src, converter);

	pair<hol_term*, unsigned int> arg1 = {nullptr, 0};
	pair<hol_term*, unsigned int> arg1_of = {nullptr, 0};
	pair<hol_term*, unsigned int> arg2 = {nullptr, 0};
	pair<hol_term*, unsigned int> arg2_of = {nullptr, 0};
	pair<hol_term*, unsigned int> predicate = {nullptr, 0};
	for (unsigned int i = 0; i < operand->array.length; i++) {
		bool found_arg1 = (arg1.key != nullptr);
		bool found_arg2 = (arg2.key != nullptr);
		bool found_arg1_of = (arg1_of.key != nullptr);
		bool found_arg2_of = (arg2_of.key != nullptr);
		bool found_predicate = (predicate.key != 0);
		if (process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1.key, arg2.key, arg1_of.key, arg2_of.key, predicate.key)) {
			if (!found_arg1 && arg1.key != nullptr) arg1.value = i;
			if (!found_arg2 && arg2.key != nullptr) arg2.value = i;
			if (!found_arg1_of && arg1_of.key != nullptr) arg1_of.value = i;
			if (!found_arg2_of && arg2_of.key != nullptr) arg2_of.value = i;
			if (!found_predicate && predicate.key != nullptr) predicate.value = i;
		}
	}

	if (predicate.key == nullptr || predicate.key->type != hol_term_type::UNARY_APPLICATION
	 || predicate.key->binary.left->type != hol_term_type::CONSTANT
	 || predicate.key->binary.left->constant != (unsigned int) built_in_predicates::INVERSE)
		return default_apply<Type>(src, converter);

	hol_term** new_operands = (hol_term**) malloc(sizeof(hol_term*) * operand->array.length);
	if (new_operands == nullptr) {
		fprintf(stderr, "apply ERROR: Insufficient memory for `new_operands` array when removing inverse events.\n");
		return nullptr;
	}
	for (unsigned int i = 0; i < operand->array.length; i++) {
		hol_term* new_operand;
		if (i == predicate.value) {
			new_operand = hol_term::new_apply(predicate.key->binary.right, operand->array.operands[i]->binary.right);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
				free(new_operands); return nullptr;
			}
			predicate.key->binary.right->reference_count++;
			operand->array.operands[i]->binary.right->reference_count++;
		} else if (i == arg1.value) {
			new_operand = hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value, operand->array.operands[i]->binary.left->binary.right), arg1.key);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
				free(new_operands); return nullptr;
			}
			hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count++;
			operand->array.operands[i]->binary.left->binary.right->reference_count++;
			arg1.key->reference_count++;
		} else if (i == arg2.value) {
			new_operand = hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value, operand->array.operands[i]->binary.left->binary.right), arg2.key);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
				free(new_operands); return nullptr;
			}
			hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count++;
			operand->array.operands[i]->binary.left->binary.right->reference_count++;
			arg2.key->reference_count++;
		} else if (i == arg1_of.value) {
			new_operand = hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2_OF>::value, arg1_of.key), operand->array.operands[i]->binary.right);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
				free(new_operands); return nullptr;
			}
			hol_term::constants<(unsigned int) built_in_predicates::ARG2_OF>::value.reference_count++;
			arg1_of.key->reference_count++;
			operand->array.operands[i]->binary.right->reference_count++;
		} else if (i == arg2_of.value) {
			new_operand = hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1_OF>::value, arg2_of.key), operand->array.operands[i]->binary.right);
			if (new_operand == nullptr) {
				for (unsigned int j = 0; j < i; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
				free(new_operands); return nullptr;
			}
			hol_term::constants<(unsigned int) built_in_predicates::ARG1_OF>::value.reference_count++;
			arg2_of.key->reference_count++;
			operand->array.operands[i]->binary.right->reference_count++;
		} else {
			new_operand = apply(operand->array.operands[i], converter);
			if (new_operand == operand->array.operands[i])
				new_operand->reference_count++;
		}
		new_operands[i] = new_operand;
	}

	hol_term* new_term = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(make_array_view(new_operands, operand->array.length)));
	if (new_term == nullptr) {
		for (unsigned int j = 0; j < operand->array.length; j++) { free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]); }
		free(new_operands); return nullptr;
	}
	free(new_operands);
	return new_term;
}

inline hol_term* remove_inverse(hol_term* src)
{
	inverse_remover converter;
	hol_term* dst = apply(src, converter);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct named_entity_collector {
	array<string>& named_entities;
};

template<hol_term_type Type>
inline bool visit(const hol_term& term, named_entity_collector& visitor) {
	if (Type == hol_term_type::EXISTS) {
		hol_term* operand = term.quantifier.operand;
		if (operand->type != hol_term_type::AND)
			return true;

		hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
		hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
		hol_term* predicate = nullptr;
		for (unsigned int i = 0; i < operand->array.length; i++) {
			if (!process_event_conjunct(term.quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate))
				return true;
		}

		if (arg2 == nullptr || predicate->type != hol_term_type::CONSTANT
		 || predicate->constant != (unsigned int) built_in_predicates::NAME || arg2->type != hol_term_type::STRING)
			return true;

		if (!visitor.named_entities.ensure_capacity(visitor.named_entities.length + 1)
		 || !init(visitor.named_entities[visitor.named_entities.length], arg2->str))
			return false;
		visitor.named_entities.length++;
	}
	return true;
}

inline bool get_named_entities(const hol_term& src, array<string>& named_entities) {
	named_entity_collector visitor = {named_entities};
	return visit(src, visitor);
}

struct set_operation_normalizer {
	array<unsigned int> variable_stack;
	array_map<unsigned int, hol_term*> maximal_subsets;

	set_operation_normalizer() : variable_stack(8), maximal_subsets(4) { }
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, set_operation_normalizer& normalizer) {
	if (src->binary.left->type == hol_term_type::CONSTANT
	 && src->binary.left->constant == (unsigned int) built_in_predicates::SIZE
	 && src->binary.right->type == hol_term_type::CONSTANT)
	{
		unsigned int max_variable = normalizer.variable_stack[0];
		for (unsigned int i = 1; i < normalizer.variable_stack.length; i++)
			max_variable = max(max_variable, normalizer.variable_stack[i]);
		hol_term* new_term = hol_term::new_apply(src->binary.left, hol_term::new_lambda(max_variable + 1, hol_term::new_apply(src->binary.left, hol_term::new_variable(max_variable + 1))));
		if (new_term == nullptr)
			return nullptr;
		src->binary.left->reference_count += 2;
		return new_term;
	}
	return default_apply<Type>(src, normalizer);
}

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::BINARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, set_operation_normalizer& normalizer) {
	if (src->ternary.first->type == hol_term_type::CONSTANT
	 && src->ternary.first->constant == (unsigned int) built_in_predicates::SUBSET)
	{
		hol_term* antecedent;
		hol_term* consequent;

		unsigned int variable;
		if (src->ternary.second->type == hol_term_type::LAMBDA) {
			variable = src->ternary.second->quantifier.variable;

			/* check that the variable is not used in the consequent */
			if (src->ternary.third->type == hol_term_type::LAMBDA
			 && src->ternary.third->quantifier.variable != variable)
			{
				array<unsigned int> variables(8);
				if (!get_bound_variables(*src->ternary.third->quantifier.operand, variables))
					return nullptr;
				if (variables.contains(variable)) {
					/* we can't use the antecedent variable, so check if we can use the consequent variable */
					variable = src->ternary.third->quantifier.variable;
					if (!get_bound_variables(*src->ternary.second->quantifier.operand, variables))
						return nullptr;
					if (variables.contains(variable)) {
						/* we can't use either variable, so let's use a new one */
						variable = normalizer.variable_stack[0];
						for (unsigned int i = 1; i < normalizer.variable_stack.length; i++)
							variable = max(variable, normalizer.variable_stack[i]);
						for (unsigned int i = 0; i < variables.length; i++)
							variable = max(variable, variables[i]);
					}
				}
			}
		} else if (src->ternary.third->type == hol_term_type::LAMBDA) {
			variable = src->ternary.third->quantifier.variable;
		} else {
			/* neither antecedent or consequent are lambda expressions, so just use a new variable */
			variable = normalizer.variable_stack[0];
			for (unsigned int i = 1; i < normalizer.variable_stack.length; i++)
				variable = max(variable, normalizer.variable_stack[i]);
		}

		hol_term* dst_var = hol_term::new_variable(variable);
		if (dst_var == nullptr)
			return nullptr;

		if (src->ternary.second->type == hol_term_type::LAMBDA) {
			antecedent = src->ternary.second->quantifier.operand;
			if (src->ternary.second->quantifier.variable != variable) {
				hol_term* src_var = hol_term::new_variable(src->ternary.second->quantifier.variable);
				if (src_var == nullptr) {
					free(*dst_var); free(dst_var);
					return nullptr;
				}
				antecedent = substitute(antecedent, src_var, dst_var);
				free(*src_var); free(src_var);
				if (antecedent == nullptr) {
					free(*dst_var); free(dst_var);
					return nullptr;
				}
			} else {
				antecedent->reference_count++;
			}
		} else {
			antecedent = hol_term::new_apply(src->ternary.second, dst_var);
			src->ternary.second->reference_count++;
			dst_var->reference_count++;
		}

		if (src->ternary.third->type == hol_term_type::LAMBDA) {
			consequent = src->ternary.third->quantifier.operand;
			if (src->ternary.third->quantifier.variable != variable) {
				hol_term* src_var = hol_term::new_variable(src->ternary.third->quantifier.variable);
				if (src_var == nullptr) {
					free(*dst_var); free(dst_var);
					return nullptr;
				}
				consequent = substitute(consequent, src_var, dst_var);
				free(*src_var); free(src_var);
				if (consequent == nullptr) {
					free(*dst_var); free(dst_var);
					return nullptr;
				}
			} else {
				consequent->reference_count++;
			}
		} else {
			consequent = hol_term::new_apply(src->ternary.third, dst_var);
			src->ternary.third->reference_count++;
			dst_var->reference_count++;
		}

		free(*dst_var); if (dst_var->reference_count == 0) free(dst_var);
		hol_term* new_term = hol_term::new_for_all(variable, hol_term::new_if_then(antecedent, consequent));
		if (new_term == nullptr) {
			free(*antecedent); if (antecedent->reference_count == 0) free(antecedent);
			free(*consequent); if (consequent->reference_count == 0) free(consequent);
			return nullptr;
		}
		return new_term;
	}
	return default_apply<Type>(src, normalizer);
}

struct dependent_scope_finder {
	unsigned int variable;

	dependent_scope_finder(unsigned int variable) : variable(variable) { }

	inline void operator () (
			hol_term* src, hol_term*& head,
			head_index& predicate_index)
	{
		if (src->type == hol_term_type::EXISTS || src->type == hol_term_type::FOR_ALL) {
			hol_term* operand = src->quantifier.operand;
			while (true) {
				if (operand->type == hol_term_type::AND || operand->type == hol_term_type::OR || operand->type == hol_term_type::IFF) {
					if (operand->array.operands[operand->array.length - 1]->type == hol_term_type::EQUALS
					 || operand->array.operands[operand->array.length - 1]->type == hol_term_type::UNARY_APPLICATION
					 || operand->array.operands[operand->array.length - 1]->type == hol_term_type::BINARY_APPLICATION
					 || operand->array.operands[operand->array.length - 1]->type == hol_term_type::CONSTANT
					 || operand->array.operands[operand->array.length - 1]->type == hol_term_type::VARIABLE)
					{
						head = src;
						return;
					}
					for (unsigned int i = 0; i + 1 < operand->array.length; i++) {
						if (is_variable_free(*operand->array.operands[i], variable)) {
							head = src;
							return;
						}
					}
					operand = operand->array.operands[operand->array.length - 1];
				} else if (operand->type == hol_term_type::IF_THEN) {
					if (operand->binary.right->type == hol_term_type::EQUALS
					 || operand->binary.right->type == hol_term_type::UNARY_APPLICATION
					 || operand->binary.right->type == hol_term_type::BINARY_APPLICATION
					 || operand->binary.right->type == hol_term_type::CONSTANT
					 || operand->binary.right->type == hol_term_type::VARIABLE)
					{
						head = src;
						return;
					}
					if (is_variable_free(*operand->binary.left, variable)) {
						head = src;
						return;
					}
					operand = operand->binary.right;
				} else {
					break;
				}
			}

			head = nullptr;
			return;
		} else if (src->type == hol_term_type::AND || src->type == hol_term_type::OR || src->type == hol_term_type::IFF || src->type == hol_term_type::IF_THEN) {
			head = nullptr;
			return;
		}
		head = src;
	}
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, set_operation_normalizer& normalizer)
{
	normalizer.variable_stack.add(src->quantifier.variable);

	hol_term* dst = nullptr;
	if (Type == hol_term_type::EXISTS && src->quantifier.operand->type == hol_term_type::AND) {
		hol_term* operand = src->quantifier.operand;
		hol_term* left = operand->array.operands[0];
		if (left->type == hol_term_type::BINARY_APPLICATION && left->ternary.first->type == hol_term_type::CONSTANT
		 && left->ternary.first->constant == (unsigned int) built_in_predicates::MAXIMAL_SUBSET
		 && left->ternary.second->type == hol_term_type::VARIABLE && left->ternary.second->variable == src->quantifier.variable
		 && left->ternary.third->type == hol_term_type::LAMBDA && left->ternary.third->quantifier.operand->type != hol_term_type::LAMBDA)
		{
			hol_term* set_definition = left->ternary.third->quantifier.operand;
			array<hol_term*> new_operands(operand->array.length);
			unsigned int index = normalizer.maximal_subsets.size;
			normalizer.maximal_subsets.put(src->quantifier.variable, nullptr);
			unsigned int operand_with_universal_index = operand->array.length;
			for (unsigned int i = 1; i < operand->array.length; i++) {
				hol_term* new_operand = apply(operand->array.operands[i], normalizer);
				if (new_operand == nullptr) {
					free_all(new_operands);
					return nullptr;
				} else if (new_operand == operand->array.operands[i]) {
					new_operand->reference_count++;
				}
				if (operand_with_universal_index == operand->array.length && normalizer.maximal_subsets.values[index] != nullptr)
					operand_with_universal_index = new_operands.length;
				new_operands[new_operands.length++] = new_operand;
			}
			hol_term* inner_operand = normalizer.maximal_subsets.values[index];
			normalizer.maximal_subsets.size--;

			unsigned int max_variable = 0;
			max_bound_variable(*left->ternary.third, max_variable);
			max_bound_variable(*inner_operand, max_variable);
			max_bound_variable(*new_operands[operand_with_universal_index], max_variable);
			for (unsigned int variable : normalizer.variable_stack)
				max_variable = max(max_variable, variable);

			unsigned int new_variable = ++max_variable;
			hol_term* dst_variable = hol_term::new_variable(new_variable);
			if (dst_variable == nullptr) {
				free_all(new_operands);
				return nullptr;
			}
			hol_term* src_variable = hol_term::new_variable(left->ternary.third->quantifier.variable);
			if (src_variable == nullptr) {
				free(*dst_variable); free(dst_variable);
				free_all(new_operands); return nullptr;
			}
			hol_term* new_set_definition = substitute<hol_term_type::VARIABLE>(set_definition, src_variable, dst_variable);
			free(*src_variable); free(src_variable);
			if (new_set_definition == nullptr) {
				free(*dst_variable); free(dst_variable);
				free_all(new_operands); return nullptr;
			}

			src_variable = hol_term::new_variable(inner_operand->quantifier.variable);
			if (src_variable == nullptr) {
				free(*new_set_definition); if (new_set_definition->reference_count == 0) free(new_set_definition);
				free(*dst_variable); free(dst_variable);
				free_all(new_operands); return nullptr;
			}
			hol_term* relativized = (inner_operand->type == hol_term_type::EXISTS ? inner_operand->quantifier.operand->array.operands[1] : inner_operand->quantifier.operand->binary.right);
			hol_term* new_inner_operand = substitute<hol_term_type::VARIABLE>(relativized, src_variable, dst_variable);
			free(*src_variable); free(src_variable);
			free(*dst_variable); if (dst_variable->reference_count == 0) free(dst_variable);
			if (new_inner_operand == nullptr) {
				free(*new_set_definition); if (new_set_definition->reference_count == 0) free(new_set_definition);
				free_all(new_operands); return nullptr;
			}

			normalizer.variable_stack.add(new_variable);
			hol_term* dst_set_definition = apply(new_set_definition, normalizer);
			if (dst_set_definition == nullptr) {
				free(*new_set_definition); if (new_set_definition->reference_count == 0) free(new_set_definition);
				free(*new_inner_operand); if (new_inner_operand->reference_count == 0) free(new_inner_operand);
				free_all(new_operands); return nullptr;
			} else if (dst_set_definition != new_set_definition) {
				free(*new_set_definition); if (new_set_definition->reference_count == 0) free(new_set_definition);
			}
			hol_term* dst_inner_operand = apply(new_inner_operand, normalizer);
			if (dst_inner_operand == nullptr) {
				free(*dst_set_definition); if (dst_set_definition->reference_count == 0) free(dst_set_definition);
				free(*new_inner_operand); if (new_inner_operand->reference_count == 0) free(new_inner_operand);
				free_all(new_operands); return nullptr;
			} else if (dst_inner_operand != new_inner_operand) {
				free(*new_inner_operand); if (new_inner_operand->reference_count == 0) free(new_inner_operand);
			}
			normalizer.variable_stack.length--;

			/* find the right-most scope in `dst_inner_operand` that depends on the universally-quantified variable */
			head_index predicate_index; no_op apply;
			dependent_scope_finder head_finder(new_variable);
			hol_term* inner_operand_head = find_head(dst_inner_operand, predicate_index, head_finder, apply);

			hol_term* new_left;
			if (dst_set_definition->type == hol_term_type::AND) {
				if (inner_operand_head->type == hol_term_type::AND) {
					new_left = hol_term::new_equals(hol_term::new_variable(src->quantifier.variable),
							hol_term::new_lambda(new_variable, hol_term::new_and(make_concat_array_view(
								make_array_view(dst_set_definition->array.operands, dst_set_definition->array.length),
								make_array_view(inner_operand_head->array.operands, inner_operand_head->array.length)))));
				} else {
					new_left = hol_term::new_equals(hol_term::new_variable(src->quantifier.variable),
							hol_term::new_lambda(new_variable, hol_term::new_and(make_appended_array_view(
								make_array_view(dst_set_definition->array.operands, dst_set_definition->array.length), inner_operand_head))));
				}
			} else {
				if (inner_operand_head->type == hol_term_type::AND) {
					new_left = hol_term::new_equals(hol_term::new_variable(src->quantifier.variable),
							hol_term::new_lambda(new_variable, hol_term::new_and(make_prepended_array_view(
								dst_set_definition, make_array_view(inner_operand_head->array.operands, inner_operand_head->array.length)))));
				} else {
					new_left = hol_term::new_equals(hol_term::new_variable(src->quantifier.variable),
							hol_term::new_lambda(new_variable, hol_term::new_and(dst_set_definition, inner_operand_head)));
				}
			}
			if (new_left == nullptr) {
				free(*dst_set_definition); if (dst_set_definition->reference_count == 0) free(dst_set_definition);
				free(*dst_inner_operand); if (dst_inner_operand->reference_count == 0) free(dst_inner_operand);
				free_all(new_operands); return nullptr;
			}
			if (dst_set_definition->type == hol_term_type::AND) {
				for (unsigned int i = 0; i < dst_set_definition->array.length; i++)
					dst_set_definition->array.operands[i]->reference_count++;
				free(*dst_set_definition); if (dst_set_definition->reference_count == 0) free(dst_set_definition);
			} if (inner_operand_head->type == hol_term_type::AND) {
				for (unsigned int i = 0; i < inner_operand_head->array.length; i++)
					inner_operand_head->array.operands[i]->reference_count++;
			} else {
				inner_operand_head->reference_count++;
			}

			hol_term* outer = new_operands[operand_with_universal_index];
			shift_left(new_operands.data + operand_with_universal_index, new_operands.length - operand_with_universal_index - 1);
			new_operands.length--;

			dst = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(
					make_prepended_array_view(new_left, make_array_view(new_operands.data, new_operands.length))));
			if (dst == nullptr) {
				free(*dst_inner_operand); if (dst_inner_operand->reference_count == 0) free(dst_inner_operand);
				free(*outer); if (outer->reference_count == 0) free(outer);
				free(*new_left); free(new_left);
				free_all(new_operands); return nullptr;
			}

			hol_term* temp = substitute_head<any_node_position::NONE>(dst_inner_operand, inner_operand_head, dst);
			free(*dst_inner_operand); if (dst_inner_operand->reference_count == 0) free(dst_inner_operand);
			free(*dst); if (dst->reference_count == 0) free(dst);
			if (temp == nullptr) {
				free(*outer); if (outer->reference_count == 0) free(outer);
				return nullptr;
			}
			dst = temp;

			temp = substitute<hol_term_type::CONSTANT>(outer, &HOL_ZERO, dst);
			free(*outer); if (outer->reference_count == 0) free(outer);
			free(*dst); if (dst->reference_count == 0) free(dst);
			if (temp == nullptr)
				return nullptr;
			dst = temp;
		} else if (src->quantifier.operand->array.length == 2
				&& src->quantifier.operand->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
				&& src->quantifier.operand->array.operands[0]->binary.left->type == hol_term_type::VARIABLE
				&& src->quantifier.operand->array.operands[0]->binary.right->type == hol_term_type::VARIABLE
				&& src->quantifier.operand->array.operands[0]->binary.right->variable == src->quantifier.variable)
		{
			unsigned int index = normalizer.maximal_subsets.index_of(src->quantifier.operand->array.operands[0]->binary.left->variable);
			if (index != normalizer.maximal_subsets.size) {
				normalizer.maximal_subsets.values[index] = src;
				dst = &HOL_ZERO;
				HOL_ZERO.reference_count++;
			}
		}
	} else if (Type == hol_term_type::FOR_ALL && src->quantifier.operand->type == hol_term_type::IF_THEN
			&& src->quantifier.operand->binary.left->type == hol_term_type::UNARY_APPLICATION
			&& src->quantifier.operand->binary.left->binary.left->type == hol_term_type::VARIABLE
			&& src->quantifier.operand->binary.left->binary.right->type == hol_term_type::VARIABLE
			&& src->quantifier.operand->binary.left->binary.right->variable == src->quantifier.variable)
	{
		unsigned int index = normalizer.maximal_subsets.index_of(src->quantifier.operand->binary.left->binary.left->variable);
		if (index != normalizer.maximal_subsets.size) {
			normalizer.maximal_subsets.values[index] = src;
			dst = &HOL_ZERO;
			HOL_ZERO.reference_count++;
		}
	}

	if (dst == nullptr)
		dst = default_apply<Type>(src, normalizer);

#if !defined(NDEBUG)
	if (normalizer.variable_stack.length == 0) {
		fprintf(stderr, "apply ERROR: Quantified term is not well-formed while normalizing set operations.\n");
		return dst;
	}
#endif
	normalizer.variable_stack.length--;
	return dst;
}

inline hol_term* normalize_set_operations(hol_term* src)
{
	set_operation_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct conjunction_iterator {
	array<pair<hol_term*, unsigned int>> stack;

	conjunction_iterator(hol_term* conjunction) : stack(4) {
		stack[0].key = conjunction;
		stack[0].value = 0;
		stack.length = 1;
	}

	hol_term* next() {
		pair<hol_term*, unsigned int>& current = stack[stack.length - 1];
		if (current.value < current.key->array.length) {
			hol_term* next_item = current.key->array.operands[current.value];
			if (next_item->type == hol_term_type::AND) {
				if (!stack.add(make_pair(next_item, 0u)))
					return nullptr;
				return next();
			} else {
				current.value++;
				return next_item;
			}
		}

		stack.length--;
		if (stack.length == 0)
			return nullptr;
		stack[stack.length - 1].value++;
		return next();
	}
};

struct equality_quantifier_normalizer { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS>::type* = nullptr>
inline hol_term* apply(hol_term* src, equality_quantifier_normalizer& normalizer) {
	hol_term* operand = src->quantifier.operand;
	if (operand->type != hol_term_type::AND) {
		hol_term* other = nullptr;
		if (operand->type == hol_term_type::EQUALS) {
			if (operand->binary.left->type == hol_term_type::VARIABLE && operand->binary.left->variable == src->quantifier.variable
			 && (operand->binary.right->type == hol_term_type::CONSTANT || operand->binary.right->type == hol_term_type::VARIABLE))
			{
				other = operand->binary.right;
			} else if (operand->binary.right->type == hol_term_type::VARIABLE && operand->binary.right->variable == src->quantifier.variable
					&& (operand->binary.left->type == hol_term_type::CONSTANT || operand->binary.left->type == hol_term_type::VARIABLE))
			{
				other = operand->binary.left;
			}
		}
		if (other != nullptr) {
			HOL_TRUE.reference_count++;
			return &HOL_TRUE;
		} else {
			return default_apply<Type>(src, normalizer);
		}
	} else {
		hol_term* var = nullptr;
		hol_term* substitution = nullptr;
		array<hol_term*> new_conjuncts(8);
		conjunction_iterator iterator(operand);
		hol_term* conjunct = iterator.next();
		while (conjunct != nullptr) {
			hol_term* other = nullptr;
			if (substitution == nullptr && conjunct->type == hol_term_type::EQUALS) {
				if (conjunct->binary.left->type == hol_term_type::VARIABLE && conjunct->binary.left->variable == src->quantifier.variable
				 && (conjunct->binary.right->type == hol_term_type::CONSTANT || conjunct->binary.right->type == hol_term_type::VARIABLE))
				{
					var = conjunct->binary.left;
					other = conjunct->binary.right;
				} else if (conjunct->binary.right->type == hol_term_type::VARIABLE && conjunct->binary.right->variable == src->quantifier.variable
						&& (conjunct->binary.left->type == hol_term_type::CONSTANT || conjunct->binary.left->type == hol_term_type::VARIABLE))
				{
					var = conjunct->binary.right;
					other = conjunct->binary.left;
				}
			}
			if (other == nullptr) {
				if (!new_conjuncts.ensure_capacity(new_conjuncts.length + 1)) {
					free_all(new_conjuncts);
					return nullptr;
				}
				new_conjuncts[new_conjuncts.length] = apply(conjunct, normalizer);
				if (new_conjuncts[new_conjuncts.length] == conjunct)
					new_conjuncts[new_conjuncts.length]->reference_count++;
				new_conjuncts.length++;
			} else {
				if (other->type == hol_term_type::VARIABLE && other->variable == src->quantifier.variable) {
					conjunct = iterator.next();
					continue;
				}
				substitution = other;
			}
			conjunct = iterator.next();
		}

		if (new_conjuncts.length == 0) {
			HOL_TRUE.reference_count++;
			return &HOL_TRUE;
		}

		if (substitution == nullptr) {
			bool same_as_src = (new_conjuncts.length == operand->array.length);
			for (unsigned int i = 0; same_as_src && i < new_conjuncts.length; i++)
				if (new_conjuncts[i] != operand->array.operands[i]) same_as_src = false;
			if (same_as_src) {
				for (unsigned int i = 0; i < new_conjuncts.length; i++)
					free(*new_conjuncts[i]);
				return src;
			} else {
				hol_term* dst;
				if (new_conjuncts.length == 1) {
					dst = hol_term::new_exists(src->quantifier.variable, new_conjuncts[0]);
				} else {
					dst = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(make_array_view(new_conjuncts.data, new_conjuncts.length)));
				}
				if (dst == nullptr) {
					free_all(new_conjuncts);
					return nullptr;
				}
				return dst;
			}
		}

		for (unsigned int i = 0; i < new_conjuncts.length; i++) {
			hol_term* substituted = substitute(new_conjuncts[i], var, substitution);
			if (substituted == nullptr) {
				free_all(new_conjuncts);
				return nullptr;
			}
			free(*new_conjuncts[i]); if (new_conjuncts[i]->reference_count == 0) free(new_conjuncts[i]);
			new_conjuncts[i] = substituted;
		}

		hol_term* dst;
		if (new_conjuncts.length == 1) {
			dst = new_conjuncts[0];
		} else {
			dst = hol_term::new_and(make_array_view(new_conjuncts.data, new_conjuncts.length));
		}
		if (dst == nullptr) {
			free_all(new_conjuncts);
			return nullptr;
		}
		return dst;
	}
}

inline hol_term* normalize_quantifiers_with_equality(hol_term* src)
{
	equality_quantifier_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct negative_remover {
	unsigned int target_variable;
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, negative_remover& normalizer)
{
	hol_term* operand = src->quantifier.operand;
	hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
	hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
	hol_term* predicate = nullptr;
	hol_term* dst = nullptr;
	if (operand->type == hol_term_type::AND) {
		unsigned int predicate_index = operand->array.length;
		for (unsigned int i = 0; i < operand->array.length; i++) {
			if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate)) {
				predicate = nullptr;
				break;
			}
			if (predicate != nullptr && predicate_index == operand->array.length)
				predicate_index = i;
		}

		if (predicate != nullptr && predicate->type == hol_term_type::UNARY_APPLICATION
		 && predicate->binary.left->type == hol_term_type::CONSTANT
		 && predicate->binary.left->constant == (unsigned int) built_in_predicates::NEGATIVE
		 && arg2 != nullptr && arg2->type == hol_term_type::VARIABLE
		 && arg2->variable == normalizer.target_variable)
		{
			hol_term* new_predicate = hol_term::new_apply(predicate->binary.right, operand->array.operands[predicate_index]->binary.right);
			if (new_predicate == nullptr)
				return nullptr;
			predicate->binary.right->reference_count++;
			operand->array.operands[predicate_index]->binary.right->reference_count++;

			if (Type == hol_term_type::EXISTS)
				dst = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(
						make_replaced_array_view(make_array_view(operand->array.operands, operand->array.length), new_predicate, predicate_index)));
			else if (Type == hol_term_type::FOR_ALL)
				dst = hol_term::new_for_all(src->quantifier.variable, hol_term::new_and(
						make_replaced_array_view(make_array_view(operand->array.operands, operand->array.length), new_predicate, predicate_index)));
			else if (Type == hol_term_type::LAMBDA)
				dst = hol_term::new_lambda(src->quantifier.variable, hol_term::new_and(
						make_replaced_array_view(make_array_view(operand->array.operands, operand->array.length), new_predicate, predicate_index)));
			if (dst == nullptr) {
				free(*new_predicate); free(new_predicate);
				return nullptr;
			}
			for (unsigned int i = 0; i < operand->array.length; i++) {
				if (i == predicate_index) continue;
				operand->array.operands[i]->reference_count++;
			}
		}
	}

	if (dst == nullptr)
		dst = default_apply<Type>(src, normalizer);

	return dst;
}

inline hol_term* remove_negative(hol_term* src, unsigned int target_variable)
{
	negative_remover normalizer = {target_variable};
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct optimization_normalizer {
	array_map<unsigned int, hol_term*> scopes;
	array_map<unsigned int, hol_term*> negate_functions;

	optimization_normalizer() : scopes(8), negate_functions(4) { }
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, optimization_normalizer& normalizer) {
	if (!normalizer.scopes.put(src->quantifier.variable, src))
		return nullptr;

	hol_term* dst = nullptr;
	if (Type == hol_term_type::EXISTS) {
		hol_term* operand = src->quantifier.operand;

		hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
		hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
		hol_term* predicate = nullptr;
		unsigned int predicate_index = operand->array.length;
		if (operand->type == hol_term_type::AND) {
			for (unsigned int i = 0; i < operand->array.length; i++) {
				if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate)) {
					predicate = nullptr;
					break;
				}
				if (predicate_index == operand->array.length && predicate != nullptr
				 && predicate->type == hol_term_type::UNARY_APPLICATION
				 && predicate->binary.left->type == hol_term_type::CONSTANT
				 && predicate->binary.left->constant == (unsigned int) built_in_predicates::GREATEST
				 && predicate->binary.right->type == hol_term_type::VARIABLE)
				{
					predicate_index = i;
				}
			}
		}

		if (predicate_index != operand->array.length) {
			unsigned int function_variable = predicate->binary.right->variable;
			hol_term* function_definition = normalizer.scopes.get(function_variable);
			hol_term* left;
			if (function_definition->quantifier.operand->type == hol_term_type::AND) {
				left = function_definition->quantifier.operand->array.operands[0];
			} else {
				left = function_definition->quantifier.operand;
			}

			if (left->type == hol_term_type::EQUALS && left->binary.left->type == hol_term_type::VARIABLE
			 && left->binary.left->variable == function_variable
			 && left->binary.right->type == hol_term_type::LAMBDA
			 && left->binary.right->quantifier.operand->type == hol_term_type::LAMBDA
			 && left->binary.right->quantifier.operand->quantifier.operand->type == hol_term_type::EXISTS)
			{
				hol_term* function = left->binary.right->quantifier.operand->quantifier.operand;
				hol_term* new_function = remove_negative(function, left->binary.right->quantifier.operand->quantifier.variable);
				if (new_function == function) {
					function->reference_count--;
				} else {
					/* we found a `greatest(f)` scope where `f` is a `negative` function */
					if (!normalizer.negate_functions.put(function_variable, new_function)) {
						free(*new_function); free(new_function);
						return nullptr;
					}

					array<hol_term*> operands(operand->array.length);
					for (unsigned int i = 0; i < operand->array.length; i++) {
						if (i == predicate_index) {
							operands[i] = hol_term::new_apply(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::LEAST>::value, predicate->binary.right), operand->array.operands[predicate_index]->binary.right);
							if (operands[i] == nullptr) {
								free_all(operands);
								return nullptr;
							}
							hol_term::constants<(unsigned int) built_in_predicates::LEAST>::value.reference_count++;
							predicate->binary.right->reference_count++;
							operand->array.operands[predicate_index]->binary.right->reference_count++;
						} else {
							operands[i] = apply(operand->array.operands[i], normalizer);
							if (operands[i] == nullptr) {
								free_all(operands);
								return nullptr;
							} else if (operands[i] == operand->array.operands[i]) {
								operands[i]->reference_count++;
							}
						}
						operands.length++;
					}

					dst = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(make_array_view(operands.data, operands.length)));
					if (dst == nullptr) {
						free_all(operands);
						return nullptr;
					}
				}
			}
		}
	}

	if (dst == nullptr) {
		dst = default_apply<Type>(src, normalizer);
		if (dst == nullptr) return nullptr;
		else if (dst == src) dst->reference_count++;
	}

	unsigned int index = normalizer.negate_functions.index_of(dst->quantifier.variable);
	if (index != normalizer.negate_functions.size) {
		hol_term* new_function = normalizer.negate_functions.values[index];
		normalizer.negate_functions.remove_at(index);

		hol_term* left;
		if (dst->quantifier.operand->type == hol_term_type::AND) {
			left = dst->quantifier.operand->array.operands[0];
		} else {
			left = dst->quantifier.operand;
		}

		hol_term* new_left = hol_term::new_equals(left->binary.left, hol_term::new_lambda(left->binary.right->quantifier.variable, hol_term::new_lambda(left->binary.right->quantifier.operand->quantifier.variable, new_function)));
		if (new_left == nullptr) {
			free(*dst); if (dst->reference_count == 0) free(dst);
			free(*new_function); if (new_function->reference_count == 0) free(new_function);
			return nullptr;
		}
		left->binary.left->reference_count++;

		hol_term* new_dst;
		if (dst->quantifier.operand->type == hol_term_type::AND) {
			new_dst = hol_term::new_exists(dst->quantifier.variable, hol_term::new_and(make_replaced_array_view(
					make_array_view(dst->quantifier.operand->array.operands, dst->quantifier.operand->array.length), new_left, 0)));
			if (new_dst == nullptr) {
				free(*dst); if (dst->reference_count == 0) free(dst);
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
				return nullptr;
			}
			for (unsigned int i = 1; i < new_dst->quantifier.operand->array.length; i++)
				new_dst->quantifier.operand->array.operands[i]->reference_count++;
		} else {
			new_dst = hol_term::new_exists(dst->quantifier.variable, new_left);
			if (new_dst == nullptr) {
				free(*dst); if (dst->reference_count == 0) free(dst);
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
				return nullptr;
			}
		}

		free(*dst); if (dst->reference_count == 0) free(dst);
		dst = new_dst;
	}

	normalizer.scopes.size--;
	if (dst == src) dst->reference_count--;
	return dst;
}

inline hol_term* normalize_optimizations(hol_term* src)
{
	optimization_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct subset_substituter {
	unsigned int variable;
	hol_term* set;
	hol_term* lambda_var;

	subset_substituter(unsigned int variable, hol_term* set) : variable(variable), set(set) {
		lambda_var = hol_term::new_variable(set->quantifier.variable);
		if (lambda_var == nullptr) throw std::bad_alloc();
	}

	~subset_substituter() {
		free(*lambda_var); free(lambda_var);
	}
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::UNARY_APPLICATION>::type* = nullptr>
inline hol_term* apply(hol_term* src, subset_substituter& normalizer)
{
	if (src->binary.left->type == hol_term_type::VARIABLE && src->binary.left->variable == normalizer.variable) {
		return substitute(normalizer.set->quantifier.operand, normalizer.lambda_var, src->binary.right);
	} else {
		return default_apply<Type>(src, normalizer);
	}
}

inline hol_term* substitute_subsets(hol_term* src, unsigned int variable, hol_term* set)
{
	subset_substituter normalizer(variable, set);
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct subset_simplifier {
	array<unsigned int> irreducible_variables;

	subset_simplifier() : irreducible_variables(4) { }
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::UNARY_APPLICATION || Type == hol_term_type::VARIABLE>::type* = nullptr>
inline hol_term* apply(hol_term* src, subset_simplifier& normalizer)
{
	if (Type == hol_term_type::EXISTS) {
		if (src->quantifier.operand->type != hol_term_type::AND)
			return default_apply<Type>(src, normalizer);

		unsigned int operand_count = src->quantifier.operand->array.length;
		unsigned int set_definition_index = operand_count;
		array<hol_term*> new_operands(operand_count);
		for (unsigned int i = 0; i < operand_count; i++) {
			hol_term* operand = src->quantifier.operand->array.operands[i];
			if (operand->type == hol_term_type::EQUALS && operand->binary.left->type == hol_term_type::VARIABLE && operand->binary.left->variable == src->quantifier.variable
			 && operand->binary.right->type == hol_term_type::LAMBDA && operand->binary.right->quantifier.operand->type != hol_term_type::LAMBDA)
			{
				set_definition_index = i;
				hol_term* right = apply(operand->binary.right, normalizer);
				if (right == operand->binary.right) {
					new_operands[i] = operand;
				} else {
					new_operands[i] = hol_term::new_equals(operand->binary.left, right);
					if (new_operands[i] == nullptr) {
						free(*right); free(right);
						new_operands[i] = nullptr;
					} else {
						operand->binary.left->reference_count++;
					}
				}
			} else {
				new_operands[i] = apply(operand, normalizer);
			}

			if (new_operands[i] == nullptr) {
				for (unsigned int j = 0; j < i; j++) {
					if (new_operands[j] != src->quantifier.operand->array.operands[j]) {
						free(*new_operands[j]); free(new_operands[j]);
					}
				}
				return nullptr;
			}
		}

		if (set_definition_index < operand_count) {
			unsigned int index = normalizer.irreducible_variables.index_of(src->quantifier.variable);
			if (index < normalizer.irreducible_variables.length) {
				normalizer.irreducible_variables.remove(index);
			} else {
				for (unsigned int i = 0; i < operand_count; i++)
					if (new_operands[i] == src->quantifier.operand->array.operands[i]) new_operands[i]->reference_count++;
				for (unsigned int i = 0; i < operand_count; i++) {
					if (i == set_definition_index) continue;
					hol_term* new_operand = substitute_subsets(new_operands[i], src->quantifier.variable, new_operands[set_definition_index]->binary.right);
					if (new_operand == nullptr) {
						for (unsigned int j = 0; j < operand_count; j++) {
							free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]);
						}
						return nullptr;
					}
					free(*new_operands[i]); if (new_operands[i]->reference_count == 0) free(new_operands[i]);
					new_operands[i] = new_operand;
				}

				if (operand_count == 2) {
					free(*new_operands[set_definition_index]); if (new_operands[set_definition_index]->reference_count == 0) free(new_operands[set_definition_index]);
					hol_term* dst = new_operands[1 - set_definition_index];
					if (dst == src->quantifier.operand->array.operands[1 - set_definition_index])
						dst->reference_count--;
					return dst;
				} else {
					hol_term* dst = hol_term::new_and(make_excluded_array_view(new_operands.data, operand_count, set_definition_index));
					if (dst == nullptr) {
						for (unsigned int j = 0; j < operand_count; j++) {
							free(*new_operands[j]); if (new_operands[j]->reference_count == 0) free(new_operands[j]);
						}
						return nullptr;
					}
					free(*new_operands[set_definition_index]); if (new_operands[set_definition_index]->reference_count == 0) free(new_operands[set_definition_index]);
					return dst;
				}
			}
		} else {
			unsigned int index = normalizer.irreducible_variables.index_of(src->quantifier.variable);
			if (index != normalizer.irreducible_variables.length)
				normalizer.irreducible_variables.remove(index);
		}

		bool same_as_src = true;
		for (unsigned int i = 0; same_as_src && i < operand_count; i++)
			if (new_operands[i] != src->quantifier.operand->array.operands[i]) same_as_src = false;
		if (same_as_src)
			return src;

		hol_term* dst = hol_term::new_exists(src->quantifier.variable, hol_term::new_and(make_array_view(new_operands.data, operand_count)));
		if (dst == nullptr) {
			for (unsigned int j = 0; j < operand_count; j++) {
				if (new_operands[j] != src->quantifier.operand->array.operands[j]) {
					free(*new_operands[j]); free(new_operands[j]);
				}
			}
			return nullptr;
		}
		for (unsigned int i = 0; i < operand_count; i++) {
			if (new_operands[i] == src->quantifier.operand->array.operands[i])
				new_operands[i]->reference_count++;
		}
		return dst;

	} else if (Type == hol_term_type::VARIABLE) {
		if (!normalizer.irreducible_variables.contains(src->variable) && !normalizer.irreducible_variables.add(src->variable))
			return nullptr;
		return default_apply<Type>(src, normalizer);

	} else if (Type == hol_term_type::UNARY_APPLICATION) {
		if (src->binary.left->type == hol_term_type::VARIABLE) {
			hol_term* right = apply(src->binary.right, normalizer);
			if (right == nullptr) {
				return nullptr;
			} else if (right == src->binary.right) {
				return src;
			}

			hol_term* dst = hol_term::new_apply(src->binary.right, right);
			if (dst == nullptr) {
				free(*right); free(right);
				return nullptr;
			}
			src->binary.right->reference_count++;
			return dst;
		} else {
			return default_apply<Type>(src, normalizer);
		}
	}

	/* unreachable */
	return nullptr;
}

inline hol_term* simplify_subsets(hol_term* src)
{
	subset_simplifier normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct comparative_normalizer {
	array_map<unsigned int, hol_term*> scopes;
	array<unsigned int> remove_scopes;

	comparative_normalizer() : scopes(8), remove_scopes(4) { }
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, comparative_normalizer& normalizer)
{
	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA) {
		if (!normalizer.scopes.put(src->quantifier.variable, src))
			return nullptr;
	}

	hol_term* dst = nullptr;
	if (Type == hol_term_type::EXISTS) {
		/* check if this is a `greater` scope */
		hol_term* operand = src->quantifier.operand;
		if (operand->type == hol_term_type::AND) {
			hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
			hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
			hol_term* predicate = nullptr;
			for (unsigned int i = 0; i < operand->array.length; i++) {
				if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate)) {
					predicate = nullptr;
					break;
				}
			}

			if (predicate != nullptr && predicate->type == hol_term_type::UNARY_APPLICATION
			 && predicate->binary.left->type == hol_term_type::CONSTANT
			 && predicate->binary.left->constant == (unsigned int) built_in_predicates::GREATER
			 && arg1 != nullptr && arg2 != nullptr)
			{
				hol_term* function = predicate->binary.right;
				bool inverted = false;
				if (function->type == hol_term_type::UNARY_APPLICATION && function->binary.left->type == hol_term_type::CONSTANT
				 && function->binary.left->constant == (unsigned int) built_in_predicates::NEGATIVE)
				{
					function = function->binary.right;
					inverted = true;
				}

				unsigned int max_variable = 0;
				max_bound_variable(*src, max_variable);
				for (const auto& entry : normalizer.scopes)
					max_variable = max(max_variable, entry.key);

				/* check if the arg2 is a `measure` instance */
				if (arg2->type == hol_term_type::VARIABLE) {
					hol_term* other_scope = normalizer.scopes.get(arg2->variable);
					hol_term* other_operand = other_scope->quantifier.operand;
					hol_term* other_arg1 = nullptr; hol_term* other_arg1_of = nullptr;
					hol_term* other_arg2 = nullptr; hol_term* other_arg2_of = nullptr;
					hol_term* other_predicate = nullptr;
					unsigned int predicate_index = -1, arg1_index = -1, arg2_index = -1;
					if (other_operand->type == hol_term_type::AND) {
						predicate_index = other_operand->array.length;
						arg1_index = other_operand->array.length;
						arg2_index = other_operand->array.length;
						for (unsigned int i = 0; i < other_operand->array.length; i++) {
							if (other_operand->array.operands[i]->type == hol_term_type::EXISTS
							 || other_operand->array.operands[i]->type == hol_term_type::FOR_ALL
							 || other_operand->array.operands[i]->type == hol_term_type::NOT)
								continue;
							if (!process_event_conjunct(other_scope->quantifier.variable, other_operand->array.operands[i], other_arg1, other_arg2, other_arg1_of, other_arg2_of, other_predicate)) {
								other_predicate = nullptr;
								break;
							}
							if (other_predicate != nullptr && predicate_index == other_operand->array.length)
								predicate_index = i;
							else if (other_arg1 != nullptr && arg1_index == other_operand->array.length)
								arg1_index = i;
							else if (other_arg2 != nullptr && arg2_index == other_operand->array.length)
								arg2_index = i;
						}
					}

					if (other_predicate != nullptr && other_predicate->type == hol_term_type::CONSTANT
					 && other_predicate->constant == (unsigned int) built_in_predicates::MEASURE
					 && other_arg1 != nullptr && other_arg2 != nullptr)
					{
						if (!normalizer.remove_scopes.contains(other_scope->quantifier.variable)
						 && !normalizer.remove_scopes.add(other_scope->quantifier.variable))
							return (hol_term*) nullptr;

						/* arg2 is a `measure` instance */
						unsigned int value_variable = ++max_variable;
						hol_term* value_var = hol_term::new_variable(value_variable);
						if (value_var == nullptr)
							return (hol_term*) nullptr;

						unsigned int measure_variable = other_scope->quantifier.variable;
						hol_term* measure_var = other_operand->array.operands[predicate_index]->binary.right;
						unsigned int function_variable = ++max_variable;
						hol_term* function_var = hol_term::new_variable(function_variable);
						if (function_var == nullptr) {
							free(*value_var); free(value_var);
							return (hol_term*) nullptr;
						}

						dst = hol_term::new_exists(value_variable, hol_term::new_and(
								hol_term::new_exists(measure_variable, hol_term::new_and(
									other_operand->array.operands[predicate_index],
									hol_term::new_equals(other_operand->array.operands[arg1_index]->binary.left, value_var),
									other_operand->array.operands[arg2_index],
									hol_term::new_exists(function_variable, hol_term::new_and(
										hol_term::new_apply(function, function_var),
										hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value, function_var), arg1),
										hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value, function_var), measure_var)
									))
								)),
								(inverted
									? hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, other_arg1, value_var)
									: hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, value_var, other_arg1)
								)
							));
						if (dst == nullptr) {
							free(*value_var); free(value_var);
							free(*function_var); free(function_var);
							return (hol_term*) nullptr;
						}
						value_var->reference_count += 3 - 1;
						function_var->reference_count += 3 - 1;
						other_operand->array.operands[predicate_index]->reference_count++;
						other_operand->array.operands[arg1_index]->binary.left->reference_count++;
						other_operand->array.operands[arg2_index]->reference_count++;
						function->reference_count++;
						hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count++;
						arg1->reference_count++;
						hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count++;
						measure_var->reference_count++;
						hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value.reference_count++;
						other_arg1->reference_count++;
					}
				}

				if (dst == nullptr) {
					if (function->type == hol_term_type::CONSTANT && function->constant == (unsigned int) built_in_predicates::VALUE) {
						dst = (inverted
								? hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, arg2, arg1)
								: hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, arg1, arg2));
						if (dst == nullptr)
							return (hol_term*) nullptr;
						hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value.reference_count++;
						arg1->reference_count++;
						arg2->reference_count++;
					} else {
						unsigned int arg1_value_variable = ++max_variable;
						hol_term* arg1_value_var = hol_term::new_variable(arg1_value_variable);
						if (arg1_value_var == nullptr)
							return (hol_term*) nullptr;

						unsigned int arg2_value_variable = ++max_variable;
						hol_term* arg2_value_var = hol_term::new_variable(arg2_value_variable);
						if (arg2_value_var == nullptr) {
							free(*arg1_value_var); free(arg1_value_var);
							return (hol_term*) nullptr;
						}

						unsigned int function_variable = ++max_variable;
						hol_term* function_var = hol_term::new_variable(function_variable);
						if (function_var == nullptr) {
							free(*arg1_value_var); free(arg1_value_var);
							free(*arg2_value_var); free(arg2_value_var);
							return (hol_term*) nullptr;
						}

						dst = hol_term::new_exists(arg1_value_variable, hol_term::new_and(
								hol_term::new_exists(function_variable, hol_term::new_and(
									hol_term::new_apply(function, function_var),
									hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value, function_var), arg1),
									hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value, function_var), arg1_value_var)
								)),
								hol_term::new_exists(arg2_value_variable, hol_term::new_and(
									hol_term::new_exists(function_variable, hol_term::new_and(
										hol_term::new_apply(function, function_var),
										hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value, function_var), arg2),
										hol_term::new_equals(hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value, function_var), arg2_value_var)
									)),
									(inverted
										? hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, arg2_value_var, arg1_value_var)
										: hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value, arg1_value_var, arg2_value_var)
									)
								))
							));
						if (function_var == nullptr) {
							free(*arg1_value_var); free(arg1_value_var);
							free(*arg2_value_var); free(arg2_value_var);
							free(*function_var); free(function_var);
							return (hol_term*) nullptr;
						}
						arg1_value_var->reference_count += 3 - 1;
						arg2_value_var->reference_count += 3 - 1;
						function_var->reference_count += 6 - 1;
						function->reference_count += 2;
						hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count += 2;
						arg1->reference_count++;
						hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count += 2;
						hol_term::constants<(unsigned int) built_in_predicates::GREATER_THAN_OR_EQUAL>::value.reference_count++;
						arg2->reference_count++;
					}
				}
			}
		}
	}

	if (dst == nullptr) {
		if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA) {
			hol_term* operand = src->quantifier.operand;
			array<hol_term*> new_operands(operand->type == hol_term_type::AND ? operand->array.length : 1);
			bool is_same = true;
			if (operand->type == hol_term_type::AND) {
				for (unsigned int i = 0; i < operand->array.length; i++) {
					new_operands[i] = apply(operand->array.operands[i], normalizer);
					if (new_operands[i] == nullptr) {
						free_all(new_operands);
						return nullptr;
					} else if (new_operands[i] == operand->array.operands[i]) {
						new_operands[i]->reference_count++;
					} else {
						is_same = false;
					}
					new_operands.length++;
				}
			} else {
				new_operands[0] = apply(operand, normalizer);
				if (new_operands[0] == nullptr) {
					return nullptr;
				} else if (new_operands[0] == operand) {
					new_operands[0]->reference_count++;
				} else {
					is_same = false;
				}
				new_operands.length = 1;
			}

			unsigned int index = normalizer.remove_scopes.index_of(src->quantifier.variable);
			if (index != normalizer.remove_scopes.length) {
				normalizer.remove_scopes.remove(index);

				for (unsigned int i = 0; i < new_operands.length; i++) {
					if (is_variable_free(*new_operands[i], src->quantifier.variable)) {
						free(*new_operands[i]); if (new_operands[i]->reference_count == 0) free(new_operands[i]);
						new_operands.remove(i--);
					}
				}

				if (new_operands.length == 1) {
					dst = new_operands[0];
				} else {
					dst = hol_term::new_and(make_array_view(new_operands.data, new_operands.length));
				}
				if (dst == nullptr) {
					free_all(new_operands);
					return nullptr;
				}

			} else {
				if (is_same) {
					for (hol_term* new_operand : new_operands)
						new_operand->reference_count--;
					dst = src;
				} else {
					if (Type == hol_term_type::EXISTS)
						dst = hol_term::new_exists(src->quantifier.variable,
								(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
					else if (Type == hol_term_type::FOR_ALL)
						dst = hol_term::new_for_all(src->quantifier.variable,
								(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
					else if (Type == hol_term_type::LAMBDA)
						dst = hol_term::new_lambda(src->quantifier.variable,
								(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
					if (dst == nullptr) {
						free_all(new_operands);
						return nullptr;
					}
				}
			}
		} else {
			dst = default_apply<Type>(src, normalizer);
		}
	}

	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
		normalizer.scopes.size--;

	return dst;
}

inline hol_term* normalize_comparatives(hol_term* src)
{
	comparative_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct variable_name_universalizer {
	unsigned int max_variable;
	array_map<char, hol_term*> new_variables;
	array_map<unsigned int, hol_term*> variable_map;

	variable_name_universalizer(unsigned int max_variable) : max_variable(max_variable), new_variables(4), variable_map(4) { }

	~variable_name_universalizer() {
		for (auto entry : new_variables) {
			free(*entry.value);
			if (entry.value->reference_count == 0)
				free(entry.value);
		} for (auto entry : variable_map) {
			free(*entry.value);
			if (entry.value->reference_count == 0)
				free(entry.value);
		}
	}
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::IF_THEN || Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, variable_name_universalizer& normalizer)
{
	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA) {
		hol_term* operand = src->quantifier.operand;
		if (Type == hol_term_type::EXISTS && operand->type == hol_term_type::AND && operand->array.length == 3) {
			hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
			hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
			hol_term* predicate = nullptr;
			for (unsigned int i = 0; i < operand->array.length; i++) {
				if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate)) {
					predicate = nullptr;
					break;
				}
			}

			if (predicate != nullptr && arg1_of != nullptr && arg2 != nullptr
			 && predicate->type == hol_term_type::CONSTANT
			 && predicate->constant == (unsigned int) built_in_predicates::NAME
			 && arg1_of->type == hol_term_type::VARIABLE
			 && arg2->type == hol_term_type::STRING
			 && arg2->str.length == 1 && isupper(arg2->str[0]))
			{
				if (!normalizer.new_variables.ensure_capacity(normalizer.new_variables.size + 1))
					return nullptr;
				unsigned int index = normalizer.new_variables.index_of(arg2->str[0]);
				hol_term* new_variable;
				if (index != normalizer.new_variables.size) {
					new_variable = normalizer.new_variables.values[index];
				} else {
					/* we found a new "variable name" */
					new_variable = hol_term::new_variable(++normalizer.max_variable);
					if (new_variable == nullptr)
						return nullptr;
					normalizer.new_variables.keys[normalizer.new_variables.size] = arg2->str[0];
					normalizer.new_variables.values[normalizer.new_variables.size++] = new_variable;
				}

				if (!normalizer.variable_map.ensure_capacity(normalizer.variable_map.size + 1))
					return nullptr;
				index = normalizer.variable_map.index_of(arg1_of->variable);
				if (index == normalizer.variable_map.size) {
					normalizer.variable_map.keys[normalizer.variable_map.size] = arg1_of->variable;
					normalizer.variable_map.values[normalizer.variable_map.size++] = new_variable;
					new_variable->reference_count++;
				}

				HOL_TRUE.reference_count++;
				return &HOL_TRUE;
			}
		}

		/* check if we need to remove this quantifier */
		hol_term* new_operand = apply(src->quantifier.operand, normalizer);
		if (new_operand == nullptr)
			return nullptr;

		unsigned int index = normalizer.variable_map.index_of(src->quantifier.variable);
		if (index != normalizer.variable_map.size) {
			hol_term* src_var = hol_term::new_variable(src->quantifier.variable);
			if (src_var == nullptr) {
				if (new_operand != src->quantifier.operand) { free(*new_operand); free(new_operand); }
				return nullptr;
			}
			hol_term* new_variable = normalizer.variable_map.values[index];
			hol_term* substituted = substitute<hol_term_type::VARIABLE>(new_operand, src_var, new_variable);
			free(*src_var); free(src_var);
			if (new_operand != src->quantifier.operand) {
				free(*new_operand); if (new_operand->reference_count == 0) free(new_operand);
			}
			normalizer.variable_map.remove_at(index);
			free(*new_variable); if (new_variable->reference_count == 0) free(new_variable);
			return substituted;
		} else {
			if (new_operand == src->quantifier.operand)
				return src;
			hol_term* dst = nullptr;
			if (Type == hol_term_type::EXISTS)
				dst = hol_term::new_exists(src->quantifier.variable, new_operand);
			else if (Type == hol_term_type::FOR_ALL)
				dst = hol_term::new_for_all(src->quantifier.variable, new_operand);
			else if (Type == hol_term_type::LAMBDA)
				dst = hol_term::new_lambda(src->quantifier.variable, new_operand);
			if (dst == nullptr) {
				if (new_operand != src->quantifier.operand) { free(*new_operand); free(new_operand); }
				return nullptr;
			}
			if (new_operand == src->quantifier.operand) new_operand->reference_count++;
			return dst;
		}

	} else if (Type == hol_term_type::IF_THEN) {
		hol_term* left = apply(src->binary.left, normalizer);
		if (left == nullptr) return nullptr;
		hol_term* right = apply(src->binary.right, normalizer);
		if (right == nullptr) {
			if (left != src->binary.left) { free(*left); free(left); }
			return nullptr;
		}

		if (normalizer.new_variables.size == 0) {
			if (left == src->binary.left && right == src->binary.right)
				return src;
			hol_term* dst = hol_term::new_if_then(left, right);
			if (dst == nullptr) {
				if (left != src->binary.left) { free(*left); free(left); }
				if (right != src->binary.right) { free(*right); free(right); }
				return nullptr;
			}
			if (left == src->binary.left) left->reference_count++;
			if (right == src->binary.right) right->reference_count++;
			return dst;
		} else {
			hol_term* operand = hol_term::new_if_then(left, right);
			if (operand == nullptr) {
				if (left != src->binary.left) { free(*left); free(left); }
				if (right != src->binary.right) { free(*right); free(right); }
				return nullptr;
			}
			if (left == src->binary.left) left->reference_count++;
			if (right == src->binary.right) right->reference_count++;

			for (unsigned int i = 0; i < normalizer.new_variables.size; i++) {
				hol_term* new_operand = hol_term::new_for_all(normalizer.new_variables.values[i]->variable, operand);
				if (new_operand == nullptr) {
					free(*operand); free(operand);
					return nullptr;
				}
				operand = new_operand;

				free(*normalizer.new_variables.values[i]);
				if (normalizer.new_variables.values[i]->reference_count == 0)
					free(normalizer.new_variables.values[i]);
			}
			normalizer.new_variables.size = 0;
			return operand;
		}
	}

	/* unreachable */
	return nullptr;
}

inline hol_term* universalize_variable_names(hol_term* src)
{
	unsigned int max_variable = 0;
	if (!max_bound_variable(*src, max_variable))
		return nullptr;

	variable_name_universalizer normalizer(max_variable);
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

bool add_flattened_conjuncts(array<hol_term*>& conjuncts, hol_term* new_conjunct) {
	if (new_conjunct->type == hol_term_type::AND) {
		for (unsigned int i = 0; i < new_conjunct->array.length; i++) {
			if (!add_flattened_conjuncts(conjuncts, new_conjunct->array.operands[i]))
				return false;
		}
		return true;
	} else {
		return conjuncts.add(new_conjunct);
	}
}

struct property_simplifier {
	array_map<unsigned int, hol_term*> scopes;

	property_simplifier() : scopes(8) { }
};

bool is_event(hol_term* scope, const array_map<unsigned int, hol_term*>& scopes)
{
	/* check if `variable` has any arguments */
	hol_term* any_arg_term = hol_term::new_any(hol_term::new_equals(hol_term::new_apply(
			hol_term::new_any_constant(
				(unsigned int) built_in_predicates::ARG1,
				(unsigned int) built_in_predicates::ARG2,
				(unsigned int) built_in_predicates::ARG3),
			hol_term::new_variable(scope->quantifier.variable)), &HOL_ANY));
	if (any_arg_term == nullptr)
		return false;
	HOL_ANY.reference_count++;

	bool has_arg_term = has_intersection<built_in_predicates>(scope, any_arg_term);
	free(*any_arg_term); free(any_arg_term);
	if (has_arg_term) return true;

	/* check if `scope->quantifier.variable` is an element of a set of events */
	hol_term* operand = scope->quantifier.operand;
	unsigned int set_variable = 0;
	if (operand->type == hol_term_type::AND) {
		for (unsigned int i = 0; i < operand->array.length; i++) {
			hol_term* conjunct = operand->array.operands[i];
			if (conjunct->type == hol_term_type::UNARY_APPLICATION
			 && conjunct->binary.left->type == hol_term_type::VARIABLE
			 && conjunct->binary.left->variable != scope->quantifier.variable
			 && conjunct->binary.right->type == hol_term_type::VARIABLE
			 && conjunct->binary.right->variable == scope->quantifier.variable)
			{
				set_variable = conjunct->binary.left->variable;
				break;
			}
		}
	} else if (operand->type == hol_term_type::UNARY_APPLICATION
			&& operand->binary.left->type == hol_term_type::VARIABLE
			&& operand->binary.left->variable != scope->quantifier.variable
			&& operand->binary.right->type == hol_term_type::VARIABLE
			&& operand->binary.right->variable == scope->quantifier.variable)
	{
		set_variable = operand->binary.left->variable;
	}

	if (set_variable == 0)
		return false;
	hol_term* set_scope = scopes.get(set_variable);
	hol_term* set_operand = set_scope->quantifier.operand;
	hol_term* set_definition = nullptr;
	if (set_operand->type == hol_term_type::AND) {
		for (unsigned int i = 0; i < set_operand->array.length; i++) {
			hol_term* conjunct = set_operand->array.operands[i];
			if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::VARIABLE
			 && conjunct->binary.left->variable == set_variable && conjunct->binary.right->type == hol_term_type::LAMBDA
			 && conjunct->binary.right->quantifier.operand->type != hol_term_type::LAMBDA)
			{
				set_definition = conjunct->binary.right;
				break;
			}
		}
	} else if (set_operand->type == hol_term_type::EQUALS && set_operand->binary.left->type == hol_term_type::VARIABLE
			&& set_operand->binary.left->variable == set_variable && set_operand->binary.right->type == hol_term_type::LAMBDA
			&& set_operand->binary.right->quantifier.operand->type != hol_term_type::LAMBDA)
	{
		set_definition = set_operand->binary.right;
	}

	if (set_definition == nullptr)
		return false;

	return is_event(set_definition, scopes);
}

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, property_simplifier& normalizer)
{
	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA) {
		if (!normalizer.scopes.put(src->quantifier.variable, src))
			return nullptr;
	}

	if (Type == hol_term_type::EXISTS && src->quantifier.operand->type == hol_term_type::AND) {
		hol_term* operand = src->quantifier.operand;
		hol_term* arg1 = nullptr; hol_term* arg1_of = nullptr;
		hol_term* arg2 = nullptr; hol_term* arg2_of = nullptr;
		hol_term* predicate = nullptr;
		for (unsigned int i = 0; i < operand->array.length; i++) {
			if (!process_event_conjunct(src->quantifier.variable, operand->array.operands[i], arg1, arg2, arg1_of, arg2_of, predicate)) {
				predicate = nullptr;
				break;
			}
		}

		if (predicate != nullptr && arg1 != nullptr && arg2 != nullptr) {
			if (predicate->type == hol_term_type::CONSTANT && predicate->constant == (unsigned int) built_in_predicates::CARDINALITY) {
				hol_term* dst = hol_term::new_equals(
							hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::SIZE>::value, arg1), arg2);
				if (dst == nullptr)
					return nullptr;
				hol_term::constants<(unsigned int) built_in_predicates::SIZE>::value.reference_count++;
				arg1->reference_count++;
				arg2->reference_count++;
				if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
					normalizer.scopes.size--;
				return dst;
			} else if (predicate->type == hol_term_type::CONSTANT && predicate->constant == (unsigned int) built_in_predicates::VALUE) {
				hol_term* dst = hol_term::new_equals(
							hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value, arg1), arg2);
				if (dst == nullptr)
					return nullptr;
				hol_term::constants<(unsigned int) built_in_predicates::ARG2>::value.reference_count++;
				arg1->reference_count++;
				arg2->reference_count++;
				if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
					normalizer.scopes.size--;
				return dst;
			} else if (predicate->type == hol_term_type::CONSTANT && predicate->constant == (unsigned int) built_in_predicates::HAS
					&& arg2->type == hol_term_type::VARIABLE && is_event(normalizer.scopes.get(arg2->variable), normalizer.scopes))
			{
				hol_term* dst = hol_term::new_equals(
							hol_term::new_apply(&hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value, arg2), arg1);
				if (dst == nullptr)
					return nullptr;
				hol_term::constants<(unsigned int) built_in_predicates::ARG1>::value.reference_count++;
				arg1->reference_count++;
				arg2->reference_count++;
				if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
					normalizer.scopes.size--;
				return dst;
			}
		}

		/* check if this scope declares a set of function values */
		array<hol_term*> new_operands(operand->array.length);
		bool is_same = true;
		for (unsigned int i = 0; i < operand->array.length; i++) {
			hol_term* new_operand = apply(operand->array.operands[i], normalizer);
			if (new_operand == nullptr) {
				free_all(new_operands);
				return nullptr;
			} else if (new_operand == operand->array.operands[i]) {
				new_operand->reference_count++;
			} else {
				is_same = false;
			}

			if (!add_flattened_conjuncts(new_operands, new_operand)) {
				free_all(new_operands);
				return nullptr;
			}
		}

		hol_term* function = nullptr;
		hol_term* argument = nullptr;
		for (hol_term* conjunct : new_operands) {
			if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::VARIABLE
			 && conjunct->binary.left->variable == src->quantifier.variable
			 && conjunct->binary.right->type == hol_term_type::LAMBDA
			 && conjunct->binary.right->quantifier.operand->type == hol_term_type::EQUALS
			 && conjunct->binary.right->quantifier.operand->binary.left->type == hol_term_type::UNARY_APPLICATION
			 && conjunct->binary.right->quantifier.operand->binary.right->type == hol_term_type::VARIABLE
			 && conjunct->binary.right->quantifier.operand->binary.right->variable == conjunct->binary.right->quantifier.variable)
			{
				function = conjunct->binary.right->quantifier.operand->binary.left;
			} else if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::UNARY_APPLICATION
					&& conjunct->binary.left->binary.left->type == hol_term_type::CONSTANT
					&& conjunct->binary.left->binary.left->constant == (unsigned int) built_in_predicates::SIZE
					&& conjunct->binary.left->binary.right->type == hol_term_type::VARIABLE
					&& conjunct->binary.left->binary.right->variable == src->quantifier.variable
					&& conjunct->binary.right->type == hol_term_type::NUMBER)
			{
				if (conjunct->binary.right->number.integer != 1 || conjunct->binary.right->number.decimal != 0) {
					free_all(new_operands);
					HOL_FALSE.reference_count++;
					if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
						normalizer.scopes.size--;
					return &HOL_FALSE;
				}
			} else if (conjunct->type == hol_term_type::UNARY_APPLICATION
					&& conjunct->binary.left->type == hol_term_type::VARIABLE
					&& conjunct->binary.left->variable == src->quantifier.variable)
			{
				argument = conjunct->binary.right;
			} else if (conjunct->type == hol_term_type::TRUE) {
				continue;
			} else {
				function = nullptr;
				break;
			}
		}

		if (function != nullptr && argument != nullptr) {
			hol_term* dst = hol_term::new_equals(function, argument);
			if (dst == nullptr) {
				free_all(new_operands);
				return nullptr;
			}
			function->reference_count++;
			argument->reference_count++;
			free_all(new_operands);
			if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
				normalizer.scopes.size--;
			return dst;
		}

		hol_term* dst;
		if (is_same) {
			for (hol_term* new_operand : new_operands)
				new_operand->reference_count--;
			dst = src;
		} else {
			if (Type == hol_term_type::EXISTS)
				dst = hol_term::new_exists(src->quantifier.variable,
						(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
			else if (Type == hol_term_type::FOR_ALL)
				dst = hol_term::new_for_all(src->quantifier.variable,
						(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
			else if (Type == hol_term_type::LAMBDA)
				dst = hol_term::new_lambda(src->quantifier.variable,
						(new_operands.length == 1 ? new_operands[0] : hol_term::new_and(make_array_view(new_operands.data, new_operands.length))));
			if (dst == nullptr) {
				free_all(new_operands);
				return nullptr;
			}
		}
		if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
			normalizer.scopes.size--;
		return dst;
	}

	hol_term* dst = default_apply<Type>(src, normalizer);
	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA)
		normalizer.scopes.size--;
	return dst;
}

inline hol_term* simplify_properties(hol_term* src)
{
	property_simplifier normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct function_normalizer { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::LAMBDA>::type* = nullptr>
inline hol_term* apply(hol_term* src, function_normalizer& normalizer)
{
	hol_term* dst = nullptr;
	if (src->quantifier.operand->type == hol_term_type::LAMBDA && src->quantifier.operand->quantifier.operand->type == hol_term_type::EXISTS
	 && src->quantifier.operand->quantifier.operand->quantifier.operand->type == hol_term_type::AND)
	{
		hol_term* function_definition = src->quantifier.operand->quantifier.operand;
		hol_term* function_operand = function_definition->quantifier.operand;
		hol_term* set_definition = nullptr;
		hol_term* inner_function = nullptr;
		for (unsigned int i = 0; i < function_operand->array.length; i++) {
			hol_term* conjunct = function_operand->array.operands[i];
			if (conjunct->type == hol_term_type::EQUALS && conjunct->binary.left->type == hol_term_type::VARIABLE
			 && conjunct->binary.left->variable == function_definition->quantifier.variable
			 && conjunct->binary.right->type == hol_term_type::LAMBDA
			 && conjunct->binary.right->quantifier.operand->type != hol_term_type::LAMBDA)
			{
				set_definition = conjunct->binary.right;
			} else if (conjunct->type == hol_term_type::EXISTS && conjunct->quantifier.operand->type == hol_term_type::AND
					&& conjunct->quantifier.operand->array.operands[0]->type == hol_term_type::UNARY_APPLICATION
					&& conjunct->quantifier.operand->array.operands[0]->binary.left->type == hol_term_type::VARIABLE
					&& conjunct->quantifier.operand->array.operands[0]->binary.left->variable == function_definition->quantifier.variable
					&& conjunct->quantifier.operand->array.operands[0]->binary.right->type == hol_term_type::VARIABLE
					&& conjunct->quantifier.operand->array.operands[0]->binary.right->variable == conjunct->quantifier.variable)
			{
				inner_function = conjunct;
			}
		}

		if (set_definition != nullptr && inner_function != nullptr) {
			hol_term* dst_var = hol_term::new_variable(function_definition->quantifier.variable);
			if (dst_var == nullptr)
				return nullptr;

			hol_term* src_var = hol_term::new_variable(set_definition->quantifier.variable);
			if (src_var == nullptr) {
				free(*dst_var); free(dst_var);
				return nullptr;
			}

			hol_term* new_set_definition = substitute<hol_term_type::VARIABLE>(set_definition->quantifier.operand, src_var, dst_var);
			free(*src_var); free(src_var);
			if (new_set_definition == nullptr) {
				free(*dst_var); free(dst_var);
				return nullptr;
			}

			array<hol_term*> new_operands((new_set_definition->type == hol_term_type::AND ? new_set_definition->array.length : 1) + inner_function->quantifier.operand->array.length - 1);
			if (new_set_definition->type == hol_term_type::AND) {
				for (unsigned int i = 0; i < new_set_definition->array.length; i++) {
					new_operands[i] = new_set_definition->array.operands[i];
					new_operands[i]->reference_count++;
				}
				new_operands.length = new_set_definition->array.length;
				free(*new_set_definition); if (new_set_definition->reference_count == 0) free(new_set_definition);
			} else {
				new_operands[0] = new_set_definition;
				new_operands.length++;
			}

			src_var = hol_term::new_variable(inner_function->quantifier.variable);
			if (src_var == nullptr) {
				free_all(new_operands);
				free(*dst_var); free(dst_var);
				return nullptr;
			}

			for (unsigned int i = 1; i < inner_function->quantifier.operand->array.length; i++) {
				hol_term* new_operand = substitute<hol_term_type::VARIABLE>(inner_function->quantifier.operand->array.operands[i], src_var, dst_var);
				if (new_operand == nullptr) {
					free_all(new_operands);
					free(*src_var); free(src_var);
					free(*dst_var); free(dst_var);
					return nullptr;
				}
				new_operands[new_operands.length++] = new_operand;
			}
			free(*src_var); free(src_var);
			free(*dst_var); if (dst_var->reference_count == 0) free(dst_var);

			dst = hol_term::new_lambda(src->quantifier.variable, hol_term::new_lambda(src->quantifier.operand->quantifier.variable,
					hol_term::new_exists(function_definition->quantifier.variable, hol_term::new_and(make_array_view(new_operands.data, new_operands.length)))));
			if (dst == nullptr) {
				free_all(new_operands);
				return nullptr;
			}
		}
	}

	if (dst == nullptr)
		dst = default_apply<Type>(src, normalizer);
	return dst;
}

inline hol_term* normalize_functions(hol_term* src)
{
	function_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct multiple_universal_quantifier_normalizer { };

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::FOR_ALL>::type* = nullptr>
inline hol_term* apply(hol_term* src, multiple_universal_quantifier_normalizer& normalizer)
{
	array<hol_term*> quantifiers(4);
	quantifiers[0] = src;
	quantifiers.length++;
	hol_term* operand = src->quantifier.operand;
	while (operand->type == hol_term_type::FOR_ALL) {
		if (!quantifiers.add(operand))
			return nullptr;
		operand = operand->quantifier.operand;
	}

	if (operand->type == hol_term_type::IF_THEN) {
		hol_term* left = apply(operand->binary.left, normalizer);
		if (left == nullptr)
			return nullptr;
		hol_term* right = apply(operand->binary.right, normalizer);
		if (right == nullptr) {
			if (left != operand->binary.left) { free(*left); free(left); }
			return nullptr;
		}

		array<unsigned int> left_free_variables(8);
		array<unsigned int> right_free_variables(8);
		if (!get_free_variables(*left, left_free_variables)
		 || !get_free_variables(*right, right_free_variables))
		{
			if (left != operand->binary.left) { free(*left); free(left); }
			if (right != operand->binary.left) { free(*right); free(right); }
			return nullptr;
		}

		/* check if there are any univerally-quantified variables used in `left` but not `right` */
		array<hol_term*> new_outer_quantifiers(4);
		array<hol_term*> new_left_quantifiers(4);
		array<hol_term*> new_right_quantifiers(4);
		for (hol_term* quantifier : quantifiers) {
			bool in_left = left_free_variables.contains(quantifier->quantifier.variable);
			bool in_right = right_free_variables.contains(quantifier->quantifier.variable);
			if (in_left && !in_right) {
				if (!new_left_quantifiers.add(quantifier)) {
					if (left != operand->binary.left) { free(*left); free(left); }
					if (right != operand->binary.left) { free(*right); free(right); }
					return nullptr;
				}
			} else if (!in_left && in_right) {
				if (!new_right_quantifiers.add(quantifier)) {
					if (left != operand->binary.left) { free(*left); free(left); }
					if (right != operand->binary.left) { free(*right); free(right); }
					return nullptr;
				}
			} else {
				if (!new_outer_quantifiers.add(quantifier)) {
					if (left != operand->binary.left) { free(*left); free(left); }
					if (right != operand->binary.left) { free(*right); free(right); }
					return nullptr;
				}
			}
		}

		if (new_outer_quantifiers.length == quantifiers.length) {
			if (left == operand->binary.left && right == operand->binary.right)
				return src;
			hol_term* dst = hol_term::new_if_then(left, right);
			if (dst == nullptr) {
				if (left != operand->binary.left) { free(*left); free(left); }
				if (right != operand->binary.left) { free(*right); free(right); }
				return nullptr;
			}
			if (left == operand->binary.left) left->reference_count++;
			if (right == operand->binary.right) right->reference_count++;
			for (unsigned int i = quantifiers.length; i > 0; i--) {
				hol_term* quantifier = quantifiers[i - 1];
				hol_term* temp = hol_term::new_for_all(quantifier->quantifier.variable, dst);
				if (temp == nullptr) {
					free(*dst); if (dst->reference_count == 0) free(dst);
					return nullptr;
				}
				dst = temp;
			}
			return dst;
		}
		if (left == operand->binary.left) left->reference_count++;
		if (right == operand->binary.right) right->reference_count++;

		hol_term* new_left = left;
		for (unsigned int i = new_left_quantifiers.length; i > 0; i--) {
			hol_term* quantifier = new_left_quantifiers[i - 1];
			hol_term* temp = hol_term::new_exists(quantifier->quantifier.variable, new_left);
			if (temp == nullptr) {
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
				free(*right); if (right->reference_count == 0) free(right);
				return nullptr;
			}
			new_left = temp;
		}

		hol_term* new_right = right;
		for (unsigned int i = new_right_quantifiers.length; i > 0; i--) {
			hol_term* quantifier = new_right_quantifiers[i - 1];
			hol_term* temp = hol_term::new_for_all(quantifier->quantifier.variable, new_right);
			if (temp == nullptr) {
				free(*new_left); if (new_left->reference_count == 0) free(new_left);
				free(*new_right); if (new_right->reference_count == 0) free(new_right);
				return nullptr;
			}
			new_right = temp;
		}

		hol_term* dst = hol_term::new_if_then(new_left, new_right);
		if (dst == nullptr) {
			free(*new_left); if (new_left->reference_count == 0) free(new_left);
			free(*new_right); if (new_right->reference_count == 0) free(new_right);
			return nullptr;
		}
		for (unsigned int i = new_outer_quantifiers.length; i > 0; i--) {
			hol_term* quantifier = new_outer_quantifiers[i - 1];
			hol_term* temp = hol_term::new_for_all(quantifier->quantifier.variable, dst);
			if (temp == nullptr) {
				free(*dst); if (dst->reference_count == 0) free(dst);
				return nullptr;
			}
			dst = temp;
		}
		return dst;

	} else {
		return default_apply<Type>(src, normalizer);
	}
}

inline hol_term* normalize_multiple_universal_quantifiers(hol_term* src)
{
	multiple_universal_quantifier_normalizer normalizer;
	hol_term* dst = apply(src, normalizer);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct spurious_quantifier_remover {
	array_map<unsigned int, bool> is_quantifier_spurious;

	spurious_quantifier_remover() : is_quantifier_spurious(8) { }
};

template<hol_term_type Type, typename std::enable_if<Type == hol_term_type::FOR_ALL || Type == hol_term_type::EXISTS || Type == hol_term_type::VARIABLE>::type* = nullptr>
inline hol_term* apply(hol_term* src, spurious_quantifier_remover& remover)
{
	if (Type == hol_term_type::VARIABLE) {
		bool contains;
		bool& is_spurious = remover.is_quantifier_spurious.get(src->variable, contains);
		if (contains) is_spurious = false;
		return src;
	}

	if (!remover.is_quantifier_spurious.put(src->quantifier.variable, true))
		return nullptr;
	hol_term* operand = apply(src->quantifier.operand, remover);
	if (operand == nullptr)
		return nullptr;

	unsigned int index = remover.is_quantifier_spurious.index_of(src->quantifier.variable);
	if (index == remover.is_quantifier_spurious.size) {
		print("remove_spurious_quantifiers WARNING: Duplicate quantifier detected in logical form.\n", stderr);
		print("  ", stderr); print(*src, stderr, *debug_terminal_printer); print('\n', stderr);
		if (operand == src->quantifier.operand)
			operand->reference_count++;
		return operand;
	}

	bool is_spurious = remover.is_quantifier_spurious.values[index];
	remover.is_quantifier_spurious.remove_at(index);
	if (is_spurious) {
		print("WARNING: Logical form contains an unused quantified variable:\n", stderr);
		print("  ", stderr); print(*src, stderr, *debug_terminal_printer); print('\n', stderr);
		if (operand == src->quantifier.operand)
			operand->reference_count++;
		return operand;
	} else {
		if (operand == src->quantifier.operand)
			return src;
		if (Type == hol_term_type::EXISTS) {
			hol_term* new_lf = hol_term::new_exists(src->quantifier.variable, operand);
			if (new_lf == nullptr)
				return nullptr;
			operand->reference_count++;
			return new_lf;
		} else {
			hol_term* new_lf = hol_term::new_for_all(src->quantifier.variable, operand);
			if (new_lf == nullptr)
				return nullptr;
			operand->reference_count++;
			return new_lf;
		}
	}
}

inline hol_term* remove_spurious_quantifiers(hol_term* src)
{
	spurious_quantifier_remover remover;
	hol_term* dst = apply(src, remover);
	if (dst == src)
		dst->reference_count++;
	return dst;
}

enum class binding_state : uint_fast8_t {
	NOT_BOUND = 0,
	EXISTENTIAL,
	UNIVERSAL
};

struct anaphora_binder {
	const array_map<const hol_term*, const hol_term*>& assignments;
	array_map<unsigned int, hol_term*> variable_map;
	unsigned int variable_offset;

	pair<binding_state, bool>* binding_states;
	bool polarity;

	anaphora_binder(const array_map<const hol_term*, const hol_term*>& assignments, unsigned int variable_offset) :
			assignments(assignments), variable_map(8), variable_offset(variable_offset), polarity(true)
	{
		binding_states = (pair<binding_state, bool>*) calloc(assignments.size, sizeof(pair<binding_state, bool>));
		if (binding_states == nullptr)
			throw std::bad_alloc();
	}

	~anaphora_binder() {
		free(binding_states);
	}
};

template<hol_term_type Type>
inline hol_term* apply(hol_term* src, anaphora_binder& binder)
{
	pair<bool, bool>* old_binding_states = (pair<bool, bool>*) malloc(sizeof(pair<bool, bool>) * binder.assignments.size);
	if (old_binding_states == nullptr)
		return nullptr;
	for (unsigned int i = 0; i < binder.assignments.size; i++) {
		old_binding_states[i].key = (binder.binding_states[i].key != binding_state::NOT_BOUND);
		old_binding_states[i].value = binder.binding_states[i].value;
	}

	hol_term* result = nullptr;
	if (Type == hol_term_type::EXISTS || Type == hol_term_type::FOR_ALL || Type == hol_term_type::LAMBDA) {
		unsigned int index = index_of(src, binder.assignments.keys, binder.assignments.size);
		if (index < binder.assignments.size) {
			/* this is an anaphora scope */
			unsigned int var_index = 0;
			while (var_index < index) {
				if (binder.assignments.values[var_index] == binder.assignments.values[index])
					break;
				var_index++;
			}
			hol_term* new_variable = hol_term::new_variable(binder.variable_offset + var_index);
			if (new_variable == nullptr) {
				free(old_binding_states);
				return nullptr;
			} else if (!binder.variable_map.put(src->quantifier.variable, new_variable)) {
				free(old_binding_states);
				free(*new_variable); free(new_variable);
				return nullptr;
			}

			if (src->quantifier.operand->type != hol_term_type::AND) {
				fprintf(stderr, "apply ERROR: Expected anaphora scope operand to be a conjunction.\n");
				free(old_binding_states);
				return nullptr;
			}

			hol_term** new_operands = (hol_term**) malloc(sizeof(hol_term*) * src->quantifier.operand->array.length);
			unsigned int new_operand_count = 0;
			for (unsigned int i = 0; i < src->quantifier.operand->array.length; i++) {
				hol_term* operand = src->quantifier.operand->array.operands[i];
				if (operand->type == hol_term_type::UNARY_APPLICATION
				 && operand->binary.left->type == hol_term_type::CONSTANT)
				{
					if (operand->binary.left->constant == (unsigned int) built_in_predicates::PLURAL_REF
					 || operand->binary.left->constant == (unsigned int) built_in_predicates::REF)
					{
						continue;
					} else if (operand->binary.left->constant == (unsigned int) built_in_predicates::LOC_REF) {
						hol_term* new_term = hol_term::new_apply(hol_term::new_constant((unsigned int) built_in_predicates::LOCATION), new_variable);
						if (new_term == nullptr) {
							for (unsigned int i = 0; i < new_operand_count; i++) {
								free(*new_operands[i]);
								if (new_operands[i]->reference_count == 0)
									free(new_operands[i]);
							}
							free(old_binding_states);
							return nullptr;
						}
						new_variable->reference_count++;
						new_operands[new_operand_count++] = new_term;
						continue;
					}
				}

				new_operands[new_operand_count] = apply(src->quantifier.operand->array.operands[i], binder);
				if (new_operands[new_operand_count] == nullptr) {
					for (unsigned int i = 0; i < new_operand_count; i++) {
						free(*new_operands[i]);
						if (new_operands[i]->reference_count == 0)
							free(new_operands[i]);
					}
					free(old_binding_states);
					return nullptr;
				}
				if (new_operands[new_operand_count] == src->quantifier.operand->array.operands[i])
					new_operands[new_operand_count]->reference_count++;
				new_operand_count++;
			}

			if (new_operand_count == 1) {
				result = new_operands[0];
				free(new_operands);
			} else {
				if (!new_hol_term(result)) {
					for (unsigned int i = 0; i < new_operand_count; i++) {
						free(*new_operands[i]);
						if (new_operands[i]->reference_count == 0)
							free(new_operands[i]);
					}
					free(old_binding_states);
					return nullptr;
				}
				result->type = hol_term_type::AND;
				result->reference_count = 1;
				result->array.operands = new_operands;
				result->array.length = new_operand_count;
			}

			binder.variable_map.size--;
			free(*new_variable); if (new_variable->reference_count == 0) free(new_variable);
			binder.binding_states[index].value = true;
		} else {
			index = index_of(src, binder.assignments.values, binder.assignments.size);
			if (index < binder.assignments.size) {
				/* this is a referent scope */
				unsigned int var_index = 0;
				while (var_index < index) {
					if (binder.assignments.values[var_index] == binder.assignments.values[index])
						break;
					var_index++;
				}
				hol_term* new_variable = hol_term::new_variable(binder.variable_offset + var_index);
				if (new_variable == nullptr) {
					free(old_binding_states);
					return nullptr;
				} else if (!binder.variable_map.put(src->quantifier.variable, new_variable)) {
					free(old_binding_states);
					free(*new_variable); free(new_variable);
					return nullptr;
				}

				result = apply(src->quantifier.operand, binder);
				if (result == nullptr) {
					free(old_binding_states);
					return nullptr;
				}

				binder.variable_map.size--;
				free(*new_variable); if (new_variable->reference_count == 0) free(new_variable);
				if (Type == hol_term_type::EXISTS) {
					binder.binding_states[index].key = (binder.polarity ? binding_state::EXISTENTIAL : binding_state::UNIVERSAL);
				} else {
					binder.binding_states[index].key = (binder.polarity ? binding_state::UNIVERSAL : binding_state::EXISTENTIAL);
				}
			}
		}
	} else if (Type == hol_term_type::VARIABLE) {
		bool contains;
		hol_term* new_variable = binder.variable_map.get(src->variable, contains);
		if (contains) {
			result = new_variable;
			new_variable->reference_count++;
		}
	} else if (Type == hol_term_type::NOT) {
		binder.polarity = !binder.polarity;
		result = default_apply<Type>(src, binder);
		if (result == nullptr) {
			free(old_binding_states);
			return nullptr;
		}
		binder.polarity = !binder.polarity;
	} else if (Type == hol_term_type::IF_THEN) {
		binder.polarity = !binder.polarity;
		hol_term* left = apply(src->binary.left, binder);
		if (left == nullptr) {
			free(old_binding_states);
			return nullptr;
		}
		binder.polarity = !binder.polarity;

		hol_term* right = apply(src->binary.right, binder);
		if (right == nullptr) {
			free(old_binding_states);
			if (left != src->binary.left) { free(*left); if (left->reference_count == 0) free(left); }
			return nullptr;
		}

		if (left == src->binary.left && right == src->binary.right) {
			result = src;
		} else {
			result = hol_term::new_if_then(left, right);
			if (result == nullptr) {
				free(old_binding_states);
				if (left != src->binary.left) { free(*left); if (left->reference_count == 0) free(left); }
				if (right != src->binary.right) { free(*right); if (right->reference_count == 0) free(right); }
				return nullptr;
			}
			if (left == src->binary.left) left->reference_count++;
			if (right == src->binary.right) right->reference_count++;
		}
	}

	if (result == nullptr) {
		result = default_apply<Type>(src, binder);
		if (result == nullptr) {
			free(old_binding_states);
			return nullptr;
		}
	}

	/* check if this is the most recent common ancestor scope of an anaphora and its referent */
	for (unsigned int i = 0; i < binder.assignments.size; i++) {
		if (!old_binding_states[i].key && !old_binding_states[i].value
		 && binder.binding_states[i].key != binding_state::NOT_BOUND
		 && binder.binding_states[i].value)
		{
			unsigned int var_index = 0;
			while (var_index < i) {
				if (binder.assignments.values[var_index] == binder.assignments.values[i])
					break;
				var_index++;
			}

			hol_term* new_scope;
			if (binder.binding_states[i].key == binding_state::EXISTENTIAL) {
				new_scope = hol_term::new_exists(binder.variable_offset + var_index, result);
			} else {
				new_scope = hol_term::new_for_all(binder.variable_offset + var_index, result);
			}
			if (new_scope == nullptr) {
				free(old_binding_states);
				if (result != src) { free(*result); if (result->reference_count == 0) free(result); }
				return nullptr;
			}
			if (result == src) result->reference_count++;
			result = new_scope;

			binder.binding_states[i].key = binding_state::NOT_BOUND;
			binder.binding_states[i].value = false;
		}
	}
	free(old_binding_states);
	return result;
}

inline hol_term* bind_anaphora(hol_term* src,
		const array_map<const hol_term*, const hol_term*>& anaphora_assignments,
		unsigned int new_variable_offset)
{
	anaphora_binder binder(anaphora_assignments, new_variable_offset);
	hol_term* dst = apply(src, binder);
	if (dst == nullptr) {
		for (auto entry : binder.variable_map) {
			free(*entry.value);
			if (entry.value->reference_count == 0)
				free(entry.value);
		}
	}
	if (dst == src)
		dst->reference_count++;
	return dst;
}

struct referent_iterator
{
	hol_term* logical_form;
	unsigned int max_variable;

	/* an array map of referents and their parent scopes */
	array_map<const hol_term*, const hol_term*> referents;

	array_map<const hol_term*, unsigned int> anaphora;

	static inline void free(referent_iterator& iterator) {
		core::free(iterator.anaphora);
		core::free(iterator.referents);
		core::free(*iterator.logical_form);
		if (iterator.logical_form->reference_count == 0)
			core::free(iterator.logical_form);
	}

	static inline bool clone(
			const referent_iterator& src, referent_iterator& dst,
			hash_map<const hol_term*, hol_term*>& formula_map)
	{
		dst.max_variable = src.max_variable;
		if (!array_map_init(dst.referents, src.referents.capacity)) {
			return false;
		} else if (!array_map_init(dst.anaphora, src.anaphora.capacity)) {
			core::free(dst.referents);
			return false;
		} else if (!::clone(src.logical_form, dst.logical_form, formula_map)) {
			core::free(dst.referents);
			core::free(dst.anaphora);
			return false;
		}

		for (const auto& entry : src.referents) {
			if (!dst.referents.put(entry.key, entry.value)) {
				core::free(dst);
				return false;
			}
		} for (const auto& entry : src.anaphora) {
			if (!dst.anaphora.put(entry.key, entry.value)) {
				core::free(dst);
				return false;
			}
		}
		return true;
	}
};

inline bool init(referent_iterator& iterator, hol_term* logical_form) {
	iterator.max_variable = 0;
	if (!array_map_init(iterator.referents, 16)) {
		return false;
	} else if (!array_map_init(iterator.anaphora, 8)) {
		free(iterator.referents);
		return false;
	} else if (!max_bound_variable(*logical_form, iterator.max_variable)) {
		free(iterator.anaphora);
		free(iterator.referents);
		return false;
	}
	iterator.logical_form = logical_form;
	logical_form->reference_count++;
	return true;
}

inline bool get_formula_map(const referent_iterator& iterator, hash_map<const hol_term*, unsigned int>& formula_map)
{
	return get_formula_map(iterator.logical_form, formula_map);
}

template<typename Stream>
bool read(referent_iterator& iterator, Stream& in, hol_term** formulas)
{
	size_t referent_count, anaphora_count;
	unsigned int lf_index;
	if (!read(iterator.max_variable, in)
	 || !read(referent_count, in)
	 || !read(anaphora_count, in)
	 || !read(lf_index, in)
	 || !array_map_init(iterator.referents, ((size_t) 1) << (core::log2(referent_count == 0 ? 1 : referent_count) + 1)))
	{
		return false;
	} else if (!array_map_init(iterator.anaphora, ((size_t) 1) << (core::log2(anaphora_count == 0 ? 1 : anaphora_count) + 1))) {
		core::free(iterator.referents);
		return false;
	}

	iterator.logical_form = formulas[lf_index];
	iterator.logical_form->reference_count++;

	for (size_t i = 0; i < referent_count; i++) {
		unsigned int key_index;
		unsigned int value_index;
		if (!read(key_index, in)
		 || !read(value_index, in))
		{
			core::free(iterator);
			return false;
		}
		iterator.referents.keys[iterator.referents.size] = formulas[key_index];
		iterator.referents.values[iterator.referents.size++] = formulas[value_index];
	} for (size_t i = 0; i < anaphora_count; i++) {
		unsigned int key_index;
		unsigned int value;
		if (!read(key_index, in)
		 || !read(value, in))
		{
			core::free(iterator);
			return false;
		}
		iterator.anaphora.keys[iterator.anaphora.size] = formulas[key_index];
		iterator.anaphora.values[iterator.anaphora.size++] = value;
	}
	return true;
}

template<typename Stream>
bool write(const referent_iterator& iterator, Stream& out,
		const hash_map<const hol_term*, unsigned int>& formula_map)
{
	if (!write(iterator.max_variable, out)
	 || !write(iterator.referents.size, out)
	 || !write(iterator.anaphora.size, out)
	 || !write(formula_map.get(iterator.logical_form), out))
	{
		return false;
	}

	for (size_t i = 0; i < iterator.referents.size; i++) {
		if (!write(formula_map.get(iterator.referents.keys[i]), out)
		 || !write(formula_map.get(iterator.referents.values[i]), out))
		{
			return false;
		}
	} for (size_t i = 0; i < iterator.anaphora.size; i++) {
		if (!write(formula_map.get(iterator.anaphora.keys[i]), out)
		 || !write(iterator.anaphora.values[i], out))
		{
			return false;
		}
	}
	return true;
}

struct referent_iterator_state {
	const referent_iterator* iterator;
	unsigned int* indices;
	double log_probability;

	referent_iterator_state(const referent_iterator* iterator, double log_probability) :
			iterator(iterator), log_probability(log_probability)
	{
		indices = (unsigned int*) malloc(sizeof(unsigned int) * iterator->anaphora.size);
		if (indices == nullptr)
			throw std::bad_alloc();
	}

	referent_iterator_state(const referent_iterator_state& other) :
			iterator(other.iterator), log_probability(other.log_probability)
	{
		indices = (unsigned int*) malloc(sizeof(unsigned int) * iterator->anaphora.size);
		if (indices == nullptr)
			throw std::bad_alloc();
		for (unsigned int i = 0; i < iterator->anaphora.size; i++)
			indices[i] = other.indices[i];
	}

	~referent_iterator_state() {
		free(indices);
	}

	referent_iterator_state& operator = (referent_iterator_state const& other) = delete;
};

inline bool operator == (const referent_iterator_state& first, const referent_iterator_state& second) {
	if (first.iterator != second.iterator)
		return false;
	for (unsigned int i = 0; i < first.iterator->anaphora.size; i++) {
		if (first.indices[i] != second.indices[i])
			return false;
	}
	return true;
}

inline bool operator < (const referent_iterator_state& first, const referent_iterator_state& second) {
	return first.log_probability < second.log_probability;
}

inline bool is_anaphora_binding_legal(
		const referent_iterator& iterator,
		const unsigned int* binding)
{
	/* check that the anaphora don't point referents that share the same parent */
	for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
		bool contains;
		const hol_term* anaphora = iterator.anaphora.keys[i];
		unsigned int anaphora_start = iterator.anaphora.values[i];
		const hol_term* anaphora_parent = iterator.referents.get(anaphora, contains);
		if (!contains) continue;

		const hol_term* referent_parent = iterator.referents.values[anaphora_start - binding[i] - 1];
		if (anaphora_parent == referent_parent)
			return false;
	}

	/* check that no anaphora refers to another anaphora */
	for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
		unsigned int anaphora_start = iterator.anaphora.values[i];
		const hol_term* referent = iterator.referents.keys[anaphora_start - binding[i] - 1];
		if (iterator.anaphora.contains(referent))
			return false;
	}

	/* if the referent is universal, make sure it is an ancestor of the anaphora */
	for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
		const hol_term* anaphora = iterator.anaphora.keys[i];
		unsigned int anaphora_start = iterator.anaphora.values[i];
		const hol_term* referent = iterator.referents.keys[anaphora_start - binding[i] - 1];
		if (referent->type == hol_term_type::FOR_ALL && !contains_term_ptr(*referent, anaphora))
			return false;
	}
	return true;
}

inline hol_term* bind_anaphora(
	const referent_iterator& iterator,
	const unsigned int* binding)
{
	array_map<const hol_term*, const hol_term*> assignment(iterator.anaphora.size);
	for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
		unsigned int anaphora_start = iterator.anaphora.values[i];
		assignment.keys[i] = iterator.anaphora.keys[i];
		assignment.values[i] = iterator.referents.keys[anaphora_start - binding[i] - 1];
	}
	assignment.size = iterator.anaphora.size;

	return bind_anaphora(iterator.logical_form, assignment, iterator.max_variable + 1);
}

bool process_referent_iterator(
		const referent_iterator_state& state,
		std::set<referent_iterator_state>& queue,
		hol_term** resolved_logical_forms,
		double* resolved_log_probabilities,
		unsigned int& resolved_logical_form_count)
{
	if (is_anaphora_binding_legal(*state.iterator, state.indices)) {
		/* consider the case where there are no anaphora */
		if (state.iterator->anaphora.size == 0) {
			resolved_logical_forms[resolved_logical_form_count] = state.iterator->logical_form;
			resolved_log_probabilities[resolved_logical_form_count++] = state.log_probability;
			state.iterator->logical_form->reference_count++;
			return true;
		}

		/* resolve the logical form with the current anaphora assignment */
		hol_term* resolved_logical_form = bind_anaphora(*state.iterator, state.indices);
		if (resolved_logical_form != nullptr) {
			resolved_logical_forms[resolved_logical_form_count] = resolved_logical_form;
			resolved_log_probabilities[resolved_logical_form_count++] = state.log_probability;
		}
	}

	for (unsigned int i = 0; i < state.iterator->anaphora.size; i++) {
		unsigned int anaphora_start = state.iterator->anaphora.values[i];
		if (state.indices[i] + 1 >= anaphora_start)
			/* there are no more referents for this anaphora to bind to */
			continue;

		referent_iterator_state new_state(state.iterator, state.log_probability - 1.0);
		for (unsigned int j = 0; j < state.iterator->anaphora.size; j++)
			new_state.indices[j] = state.indices[j];
		new_state.indices[i]++;

		queue.insert(new_state);
	}
	return true;
}

bool get_referents(const hol_term* logical_form,
		referent_iterator& iterator,
		array_map<unsigned int, const hol_term*>& bound_variables)
{
	if (logical_form->type == hol_term_type::EXISTS || logical_form->type == hol_term_type::FOR_ALL) {
		if (!iterator.referents.put(logical_form, nullptr))
			return false;
	}

	switch (logical_form->type) {
	case hol_term_type::VARIABLE:
	case hol_term_type::CONSTANT:
	case hol_term_type::PARAMETER:
	case hol_term_type::NUMBER:
	case hol_term_type::STRING:
	case hol_term_type::UINT_LIST:
	case hol_term_type::TRUE:
	case hol_term_type::FALSE:
		return true;
	case hol_term_type::EQUALS:
		if (logical_form->binary.left->type == hol_term_type::UNARY_APPLICATION
		 && logical_form->binary.left->binary.left->type == hol_term_type::CONSTANT
		 && (logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1
		  || logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2
		  || logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG3)
		 && logical_form->binary.left->binary.right->type == hol_term_type::VARIABLE)
		{
			const hol_term* scope = bound_variables.get(logical_form->binary.left->binary.right->variable);
			size_t index = iterator.referents.index_of(scope);
			if (index < iterator.referents.size) {
				shift_left(iterator.referents.keys + index, (unsigned int) (iterator.referents.size - index - 1));
				shift_left(iterator.referents.values + index, (unsigned int) (iterator.referents.size - index - 1));
				iterator.referents.size--;
				for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
					if (iterator.anaphora.values[i] > index) {
						iterator.anaphora.values[i]--;
						if (iterator.anaphora.values[i] == 0)
							iterator.anaphora.remove_at(i--);
					}
				}
			}

			if (logical_form->binary.right->type == hol_term_type::VARIABLE) {
				/* `logical_form->binary.left->binary.right->variable` is the parent scope of `logical_form->binary.right->variable` */
				const hol_term* child_scope = bound_variables.get(logical_form->binary.right->variable);
				index = iterator.referents.index_of(child_scope);
				if (index < iterator.referents.size)
					iterator.referents.values[index] = scope;
			}
		} else if (logical_form->binary.left->type == hol_term_type::UNARY_APPLICATION
				&& logical_form->binary.left->binary.left->type == hol_term_type::CONSTANT
				&& (logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG1_OF
				 || logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG2_OF
				 || logical_form->binary.left->binary.left->constant == (unsigned int) built_in_predicates::ARG3_OF)
				&& logical_form->binary.right->type == hol_term_type::VARIABLE)
		{
			const hol_term* scope = bound_variables.get(logical_form->binary.right->variable);
			size_t index = iterator.referents.index_of(scope);
			if (index < iterator.referents.size) {
				shift_left(iterator.referents.keys + index, (unsigned int) (iterator.referents.size - index - 1));
				shift_left(iterator.referents.values + index, (unsigned int) (iterator.referents.size - index - 1));
				iterator.referents.size--;
				for (unsigned int i = 0; i < iterator.anaphora.size; i++) {
					if (iterator.anaphora.values[i] > index) {
						iterator.anaphora.values[i]--;
						if (iterator.anaphora.values[i] == 0)
							iterator.anaphora.remove_at(i--);
					}
				}
			}

			if (logical_form->binary.left->binary.right->type == hol_term_type::VARIABLE) {
				/* `logical_form->binary.right->variable` is the parent scope of `logical_form->binary.left->binary.right->variable` */
				const hol_term* child_scope = bound_variables.get(logical_form->binary.left->binary.right->variable);
				index = iterator.referents.index_of(child_scope);
				if (index < iterator.referents.size)
					iterator.referents.values[index] = scope;
			}
		} else {
			return get_referents(logical_form->binary.left, iterator, bound_variables)
				&& get_referents(logical_form->binary.right, iterator, bound_variables);
		}
		return true;
	case hol_term_type::NOT:
		if (!get_referents(logical_form->unary.operand, iterator, bound_variables))
			return false;
		return true;
	case hol_term_type::AND:
	case hol_term_type::OR:
	case hol_term_type::IFF:
		for (unsigned int i = 0; i < logical_form->array.length; i++) {
			if (!get_referents(logical_form->array.operands[i], iterator, bound_variables))
				return false;
		}
		return true;
	case hol_term_type::FOR_ALL:
	case hol_term_type::EXISTS:
	case hol_term_type::LAMBDA:
		if (!bound_variables.put(logical_form->quantifier.variable, logical_form)
		 || !get_referents(logical_form->quantifier.operand, iterator, bound_variables))
			return false;
		bound_variables.size--;
		return true;
	case hol_term_type::IF_THEN:
	case hol_term_type::UNARY_APPLICATION:
		if (logical_form->type == hol_term_type::UNARY_APPLICATION
		 && logical_form->binary.left->type == hol_term_type::CONSTANT
		 && (logical_form->binary.left->constant == (unsigned int) built_in_predicates::REF
		  || logical_form->binary.left->constant == (unsigned int) built_in_predicates::PLURAL_REF
		  || logical_form->binary.left->constant == (unsigned int) built_in_predicates::LOC_REF)
		 && logical_form->binary.right->type == hol_term_type::VARIABLE)
		{
			if (iterator.referents.size > 1
			 && !iterator.anaphora.put(bound_variables.get(logical_form->binary.right->variable), (unsigned int) iterator.referents.size - 1))
				return false;

			/* do not allow anaphora scopes to be referents */
			iterator.referents.remove(bound_variables.get(logical_form->binary.right->variable));
		} else {
			return get_referents(logical_form->binary.left, iterator, bound_variables)
				&& get_referents(logical_form->binary.right, iterator, bound_variables);
		}
		return true;
	case hol_term_type::BINARY_APPLICATION:
		if (!get_referents(logical_form->ternary.first, iterator, bound_variables)
		 || !get_referents(logical_form->ternary.second, iterator, bound_variables)
		 || !get_referents(logical_form->ternary.third, iterator, bound_variables))
			return false;
		return true;
	case hol_term_type::VARIABLE_PREIMAGE:
	case hol_term_type::ANY:
	case hol_term_type::ANY_RIGHT:
	case hol_term_type::ANY_ARRAY:
	case hol_term_type::ANY_QUANTIFIER:
	case hol_term_type::ANY_CONSTANT:
	case hol_term_type::ANY_CONSTANT_EXCEPT:
	case hol_term_type::ANY_RIGHT_ONLY:
		return true;
	}
	fprintf(stderr, "get_referents ERROR: Unrecognized hol_term_type.\n");
	return false;
}

/* TODO: for debugging; delete this */
#include "console_utils.h"

template<unsigned int K>
inline bool resolve_coreference(hol_term** logical_forms, double* log_probabilities, unsigned int& parse_count)
{
	if (parse_count == 0)
		return true;
	referent_iterator* iterators = (referent_iterator*) malloc(sizeof(referent_iterator) * parse_count);
	if (iterators == nullptr)
		return false;
	for (unsigned int i = 0; i < parse_count; i++) {
		if (!init(iterators[i], logical_forms[i])) {
			for (unsigned int j = 0; j < i; j++)
				core::free(iterators[j]);
			core::free(iterators);
			return false;
		}

		array_map<unsigned int, const hol_term*> bound_variables(8);
		if (!get_referents(logical_forms[i], iterators[i], bound_variables)) {
			for (unsigned int j = 0; j < i + 1; j++)
				core::free(iterators[j]);
			core::free(iterators);
			return false;
		}
	}

	std::set<referent_iterator_state> queue;
	for (unsigned int i = 0; i < parse_count; i++) {
		referent_iterator_state initial_state(&iterators[i], log_probabilities[i]);
		for (unsigned int j = 0; j < iterators[i].anaphora.size; j++)
			initial_state.indices[j] = 0;
		queue.insert(initial_state);
	}

	hol_term* resolved_formulas[K];
	double resolved_log_probabilities[K];
	unsigned int resolved_formula_count = 0;
	while (!queue.empty()) {
		auto last = queue.cend(); last--;
		referent_iterator_state state = *last;
		queue.erase(last);

		process_referent_iterator(state, queue, resolved_formulas, resolved_log_probabilities, resolved_formula_count);
		if (resolved_formula_count == K)
			break;
	}

/* TODO: for debugging; delete this */
bool has_anaphora = false;
for (unsigned int i = 0; !has_anaphora && i < parse_count; i++)
has_anaphora = (iterators[i].anaphora.size != 0);
if (has_anaphora) {
if (resolved_formula_count == 0) {
	print(CONSOLE_BOLD "\nCoreference resolution failed.\n" CONSOLE_RESET, stdout);
} else {
	print(CONSOLE_BOLD "\nCoreference resolution results:\n" CONSOLE_RESET, stdout);
	for (unsigned int i = 0; i < resolved_formula_count; i++) {
		print(CONSOLE_BOLD "[", stdout); print(i, stdout); print("] " CONSOLE_RESET, stdout);
		print(*resolved_formulas[i], stdout, *debug_terminal_printer); print('\n', stdout);
	}
}
}

	for (unsigned int i = 0; i < parse_count; i++) {
		core::free(iterators[i]);
		core::free(*logical_forms[i]); if (logical_forms[i]->reference_count == 0) core::free(logical_forms[i]);
	} for (unsigned int i = 0; i < resolved_formula_count; i++) {
		logical_forms[i] = resolved_formulas[i];
		log_probabilities[i] = resolved_log_probabilities[i];
	}
	core::free(iterators);
	parse_count = resolved_formula_count;
	return true;
}

#endif /* LF_UTILS_H_ */
