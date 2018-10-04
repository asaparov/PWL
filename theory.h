#ifndef THEORY_H_
#define THEORY_H_

#include <core/array.h>

#include "first_order_logic.h"

using namespace core;


enum built_in_predicates : unsigned int {
	PREDICATE_UNKNOWN = 1,
	PREDICATE_COUNT
};

inline bool add_constants_to_string_map(hash_map<string, unsigned int>& names)
{
	return names.put("unknown", PREDICATE_UNKNOWN);
}

struct relation {
	unsigned int type;
	unsigned int arg1; /* `0` here indicates the source vertex */
	unsigned int arg2; /* `0` here indicates the source vertex */
};

template<typename ProofCalculus>
struct concept
{
	array_map<unsigned int, Proof*> types;
	array_map<unsigned int, Proof*> negated_types;
	array_map<relation, Proof*> relations;
	array_map<relation, Proof*> negated_relations;
};

template<typename Formula, typename ProofCalculus>
struct theory
{
	typedef typename Formula::Type FormulaType;

	/* A map from `x` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `x(y_i)` and for any `z_i` there is an axiom in the theory
	   `~x(z_i)`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `x(u)` or `~x(u)`
	   are in the theory. */
	hash_map<unsigned int, pair<array<unsigned int>, array<unsigned int>>> types;

	/* A map from `R` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]R` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]R`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]R` or `~[u/0]R`
	   are in the theory. */
	hash_map<relation, pair<array<unsigned int>, array<unsigned int>>> relations;

	hash_map<unsigned int, concept> ground_concepts;

	array<Proof*> universal_quantications;

	hash_set<Proof*> observations;

	bool add_formula(Formula* formula)
	{
		Formula* canonicalized = canonicalize(*formula);
		if (canonicalized == NULL) {
			return false;
		}

		for (unsigned int i = 0; i < observations.length; i++) {
			if (*observations[i] == *canonicalized) {
				/* this formula has already been added to the theory */
				free(*canonicalized);
				if (canonicalized->reference_count == 0)
					free(canonicalized);
				return true;
			}
		}

		if (!observations.add(canonicalized)) {
			free(*canonicalized);
			if (canonicalized->reference_count == 0)
				free(canonicalized);
		}

		Proof* proof = make_proof(canonicalized);
		if (proof == NULL)
			return false;
		else proof->reference_count++;

		if (!proofs.add(proof)) {
			free(*proof);
			if (proof->reference_count == 0)
				free(proof);
			free(*canonicalized);
			if (canonicalized->reference_count == 0)
				free(canonicalized);
			observations.length--;
			return false;
		}

return check_proof(*proof, canonicalized);
		return true;
	}

private:
	Proof* make_proof(Formula* canonicalized)
	{
		if (canonicalized->type == FormulaType::ATOM) {
			if (canonicalized->atom.arg1.type == fol_term_type::CONSTANT)
			{
				if (canonicalized->atom.predicate != PREDICATE_UNKNOWN) {
					if (canonicalized->atom.arg1.constant == PREDICATE_UNKNOWN) {
						/* this is a definition of an object */
						canonicalized->atom.arg1.constant = PREDICATE_COUNT + definitions.length;
						return ProofCalculus::new_axiom(canonicalized);
					} else {
						/* this is a formula of form `t(c)` */
						return ProofCalculus::new_axiom(canonicalized);
					}
				} else {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				}
			}
		} else if (canonicalized->type == FormulaType::FOR_ALL) {
			unsigned int variable = canonicalized->quantifier.variable;
			if (canonicalized->quantifier.operand->type == FormulaType::IF_THEN) {
				const fol_formula* left = canonicalized->quantifier.operand->binary.left;
				if (left->type == FormulaType::ATOM
				 && left->atom.arg1.type == fol_term_type::VARIABLE
				 && left->atom.arg1.variable == variable
				 && left->atom.arg2.type == fol_term_type::NONE)
				{
					if (left->atom.predicate == PREDICATE_UNKNOWN) {
						/* this is a definition of a type */
						fol_formula* right = canonicalized->quantifier.operand.binary.right;

						/* check the right-hand side is a valid definition */
						if (!valid_definition(right, variable)) {
							fprintf(stderr, "theory.make_proof ERROR: This is not a valid type definition.\n");
							return NULL;
						}

						unsigned int new_type = PREDICATE_COUNT + definitions.length;
						fol_formula* definition = Formula::new_for_all(variable, Formula::new_iff(
								Formula::new_atom(new_type, Formula::new_variable(variable)), right));
						if (definition == NULL) return NULL;
						right->reference_count++;

						Proof* proof = ProofCalculus::new_universal_intro(
							ProofCalculus::new_implication_intro(
								ProofCalculus::new_biconditional_elim_left(
									ProofCalculus::new_universal_elim(ProofCalculus::new_axiom(definition), Formula::new_parameter(1)),
									ProofCalculus::new_axiom(left)),
								ProofCalculus::new_axiom(left)), 1);
						if (proof == NULL) {
							free(*definition); free(definition);
						}
						return proof;
					} else {
						/* this is a formula of form `![x]:(t(x) => f(x))` */
						return ProofCalculus::new_axiom(canonicalized);
					}
				} else if (left->type == FormulaType::ATOM) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				} else if (left->type == FormulaType::AND) {
					/* TODO: implement this */
					fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
					return NULL;
				} else {
					fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
					return NULL;
				}
			}
		} else if (canonicalized->type == FormulaType::AND) {
			return ProofCalculus::new_conjunction_intro(
					make_proof(canonicalized->binary.left),
					make_proof(canonicalized->binary.right));
		} else if (canonicalized->type == FormulaType::EXISTS) {
			/* TODO: implement this */
			fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
			return NULL;
		} else {
			fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
			return NULL;
		}
	}

	inline bool valid_definition(const Formula* right,
			unsigned int quantified_variable)
	{
		if (right->type == FormulaType::ATOM) {
			return right->atom.predicate != PREDICATE_UNKNOWN
				&& right->atom.arg1.type == fol_term_type::VARIABLE
				&& right->atom.arg1.variable == quantified_variable
				&& right->atom.arg2.type == NONE);
		} else if (right->type == FormulaType::AND) {
			return valid_definition(right->binary.left, quantified_variable)
				|| valid_definition(right->binary.right, quantified_variable);
		} else {
			return false;
		}
	}
};

#endif /* THEORY_H_ */
