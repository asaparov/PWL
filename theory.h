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
	typedef typename ProofCalculus::Proof Proof;

	/* A map from `x` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `x(y_i)` and for any `z_i` there is an axiom in the theory
	   `~x(z_i)`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `x(u)` or `~x(u)`
	   are in the theory. */
	hash_map<unsigned int, pair<array<unsigned int>, array<unsigned int>>> types;
	unsigned int type_count;

	/* A map from `R` to two lists of constants `{y_1, ..., y_n}` and
	   `{z_1, ..., z_m}` such that for any `y_i`, there is an axiom in the
	   theory `[y_i/0]R` and for any `z_i` there is an axiom in the theory
	   `~[z_i/0]R`. Note that this map is exhaustive, and there are no
	   other constants `u` such that the axiom `[u/0]R` or `~[u/0]R`
	   are in the theory. */
	hash_map<relation, pair<array<unsigned int>, array<unsigned int>>> relations;
	unsigned int relation_count;

	hash_map<unsigned int, concept<ProofCalculus>> ground_concepts;

	array<Proof*> universal_quantifications;

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
			return make_atom_proof(canonicalized, canonicalized);
		} else if (canonicalized->type == FormulaType::NOT) {
			if (canonicalized->unary.operand->type == FormulaType::ATOM) {
				return make_atom_proof(canonicalized, canonicalized->unary.operand);
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
						/* TODO: check that this is not implied by an existing univerally-quantified axiom */
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
		}
		fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
		return NULL;
	}

	Proof* make_atom_proof(Formula* canonicalized, Formula* atom)
	{
		if (atom->atom.arg1.type == fol_term_type::CONSTANT
		 && atom->atom.arg2.type == fol_term_type::NONE)
		{
			/* this is a unary formula */
			if (atom->atom.predicate != PREDICATE_UNKNOWN) {
				if (atom->atom.arg1.constant == PREDICATE_UNKNOWN) {
					/* this is a definition of an object */
					atom->atom.arg1.constant = PREDICATE_COUNT + definitions.length;
					return ProofCalculus::new_axiom(canonicalized);
				} else {
					/* this is a formula of form `t(c)` */

					/* check that this is not implied by an existing univerally-quantified axiom */
					// for (unsigned int k = 0; k < universal_quantifications.length; k++) {
					// 	Proof* axiom = universal_quantifications[k];
					// 	Formula* antecedent = axiom->formula->quantifier.operand->binary.left;
					// 	Formula* consequent = axiom->formula->quantifier.operand->binary.right;
					// 	unsigned int i = 0;
					// 	if (consequent->type == FormulaType::AND) {
					// 		for (; i < consequent->array.length; i++) {
					// 			fol_term dst;
					// 			if (unify(*consequent->operands[i], *canonicalized, Formula::new_variable(axiom->formula->quantifier.variable), dst))
					// 				break;
					// 		}
					// 		if (i == consequent->array.length) continue;
					// 	} else if (*consequent != *canonicalized) {
					// 		i = consequent->array.length;
					// 		continue;
					// 	}

					// 	/* check if the antecedent is satisfied by `c` */
					// 	Proof* proof;
					// 	if (!make_universal_elim_proof(consequent, ground_concepts.get(atom->atom.arg1.constant), proof))
					// 		return false;
					// 	if (proof != NULL) {
					// 		Proof* new_proof;
					// 		proof->reference_count++;
					// 		if (i == consequent->array.length) {
					// 			new_proof = ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, atom->atom.arg1), proof);
					// 		} else {
					// 			new_proof = ProofCalculus::new_conjunction_elim(
					// 					ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, atom->atom.arg1), proof),
					// 					array_view(&i, 1));
					// 		}
					// 		if (new_proof == NULL) {
					// 			free(*proof); if (proof->reference_count == 0) free(proof);
					// 			return false;
					// 		}
					// 		return new_proof;
					// 	}
					// }

					/* no existing universally-quantified axiom implies this observation */
					return ProofCalculus::new_axiom(canonicalized);
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
		} else if (atom->atom.arg1.type == fol_term_type::CONSTANT
				&& atom->atom.arg2.type == fol_term_type::CONSTANT)
		{
			/* this is a binary formula */
			if (atom->atom.predicate != PREDICATE_UNKNOWN)
			{
				if (atom->atom.arg1.constant == PREDICATE_UNKNOWN
				 || atom->atom.arg2.constant == PREDICATE_UNKNOWN) {
					fprintf(stderr, "theory.make_proof ERROR: Unsupported formula type.\n");
					return NULL;
				} else {
					/* this is a formula of form `r(c_1,c_2)` */

					/* check that this is not implied by an existing univerally-quantified axiom */
					// for (unsigned int k = 0; k < universal_quantifications.length; k++) {
					// 	Proof* axiom = universal_quantifications[k];
					// 	Formula* antecedent = axiom->formula->quantifier.operand->binary.left;
					// 	Formula* consequent = axiom->formula->quantifier.operand->binary.right;
					// 	fol_term unifying_term;
					// 	unsigned int i = 0;
					// 	if (consequent->type == FormulaType::AND) {
					// 		for (; i < consequent->array.length; i++) {
					// 			if (unify(*consequent->operands[i], *canonicalized, Formula::new_variable(axiom->formula->quantifier.variable), unifying_term))
					// 				break;
					// 		}
					// 		if (i == consequent->array.length) continue;
					// 	} else if (*consequent != *canonicalized) {
					// 		i = consequent->array.length;
					// 		continue;
					// 	}

					// 	/* check if the antecedent is satisfied by `c_1` or `c_2` (whichever unifies with the consequent) */
					// 	Proof* proof;
					// 	if (!make_universal_elim_proof(consequent, ground_concepts.get(unifying_term.constant), proof))
					// 		return false;
					// 	if (proof != NULL) {
					// 		Proof* new_proof;
					// 		proof->reference_count++;
					// 		if (i == consequent->array.length) {
					// 			new_proof = ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, canonlicalized->atom.arg1), proof);
					// 		} else {
					// 			new_proof = ProofCalculus::new_conjunction_elim(
					// 					ProofCalculus::new_implication_elim(ProofCalculus::new_universal_elim(axiom, canonlicalized->atom.arg1), proof),
					// 					array_view(&i, 1));
					// 		}
					// 		if (new_proof == NULL) {
					// 			free(*proof); if (proof->reference_count == 0) free(proof);
					// 			return false;
					// 		}
					// 		return new_proof;
					// 	}
					// }

					/* no existing universally-quantified axiom implies this observation */
					return ProofCalculus::new_axiom(canonicalized);
				}
			} else {
				/* TODO: implement this */
				fprintf(stderr, "theory.make_proof ERROR: Not implemented.\n");
				return NULL;
			}
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

	template<bool Negated = false>
	bool make_atom_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		if (constraint->type == FormulaType::ATOM) {
			bool contains;
			if (constraint->atom.arg2.type == TermType::NONE) {
				/* `constraint` is a unary literal */
				proof = (Negated ? c.negated_types.get(constraint->atom.predicate, contains) : c.types.contains(constraint->atom.predicate, contains));
				if (!contains) proof = NULL;
				return true;
			} else {
				/* `constraint` is a binary literal */
				relation r = { constraint->atom.predicate,
						(constraint->atom.arg1.type == TermType::VARIABLE ? 0 : constraint->atom.arg1.constant),
						(constraint->atom.arg2.type == TermType::VARIABLE ? 0 : constraint->atom.arg2.constant) };
				proof = (Negated ? c.negated_relations.get(r, contains) : c.relations.get(r, contains));
				if (!contains) proof = NULL;
				return true;
			}
		} else if (constraint->type == FormulaType::NOT) {
			return satisfies_atom<true>(constraint->unary.operand, c);
		} else {
			return false;
		}
	}

	bool make_universal_elim_proof(Formula* constraint, const concept<ProofCalculus>& c, Proof*& proof) const
	{
		if (constraint->type == FormulaType::AND) {
			Proof** conjunct_proofs = (Proof**) malloc(sizeof(Proof*) * constraint->array.length);
			if (conjunct_proofs == NULL) {
				fprintf(stderr, "theory.make_universal_elim_proof ERROR: Out of memory.\n");
				return false;
			}
			for (unsigned int i = 0; i < constraint->array.length; i++) {
				Formula* conjunct = constraint->array.operands[i];
				if (!make_atom_proof(conjunct, c, conjunct_proofs[i])) return false;
				if (conjunct_proofs[i] == NULL) { proof = NULL; return true; }
			}
			proof = ProofCalculus::new_conjunction_intro(array_view(conjunct_proofs, constraint->array.length));
			free(conjunct_proofs);
			return proof != NULL;
		} else {
			return make_atom_proof(constraint, c, proof);
		}
	}
};

#endif /* THEORY_H_ */
