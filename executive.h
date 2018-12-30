#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include "article.h"
#include "natural_deduction_mh.h"

constexpr double PERPLEXITY_THRESHOLD = 10.0;

inline void free_logical_forms(fol_formula** logical_forms, unsigned int count)
{
	for (unsigned int i = 0; i < count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
}

template<typename ArticleSource, typename Parser,
	typename Canonicalizer, typename ProofPrior, typename Printer>
bool read_sentence(
		const ArticleSource& articles, Parser& parser, const sentence& s,
		theory<fol_formula, natural_deduction<fol_formula>, Canonicalizer>& T,
		unsigned int article_name, ProofPrior& proof_prior, Printer& printer)
{
	unsigned int parse_count, new_constant;
	fol_formula* logical_forms[2];
	double log_probabilities[2];
	while (true) {
		/* attempt to parse the sentence */
		array<token> unrecognized = array<token>(16);
		if (!parser.template parse<2>(s, logical_forms, log_probabilities, parse_count, T, unrecognized)) {
			print("read_sentence ERROR: Unable to parse sentence '", stderr); print(s, stderr, printer); print("'.\n", stderr);
			return false;
		}

		/* read the article on the unrecognized word */
		if (unrecognized.length == 0) {
			break;
		} else if (parse_count > 0 && unrecognized.length == 1 && unrecognized[0].id == article_name) {
			/* this could be a definition so try adding it to the theory */
			new_constant = 0;
			bool success = T.add_formula(logical_forms[0], new_constant)
						&& parser.add_definition(s, logical_forms[0], new_constant);
			free_logical_forms(logical_forms, parse_count);
			return success;
		}

		/* find an article in order to learn about the unrecognized word */
		unsigned int next_article = unrecognized[0].id;
		if (next_article == article_name) {
			if (unrecognized.length == 1) {
				fprintf(stderr, "read_sentence ERROR: Unable to parse definitional sentence.\n");
				free_logical_forms(logical_forms, parse_count);
				return false;
			}
			next_article = unrecognized[1].id;
		}
		if (!read_article(next_article, articles, parser, T, proof_prior, printer)) {
			free_logical_forms(logical_forms, parse_count);
			return false;
		}
		free_logical_forms(logical_forms, parse_count);
	}

	if (parse_count == 0) {
		fprintf(stderr, "read_sentence ERROR: Given sentence has no valid parses.\n");
		return false;
	} else if (parse_count > 1 && (log_probabilities[0] - log_probabilities[1]) / s.length < PERPLEXITY_THRESHOLD) {
		/* this parse is too ambiguous */
		free_logical_forms(logical_forms, parse_count);
		return true;
	}

	/* add the most probable logical form to the theory */
	if (!T.add_formula(logical_forms[0], new_constant)) {
		free_logical_forms(logical_forms, parse_count);
		return false;
	}

// TODO: uncomment this; and also add it after the earlier call to `add_formula` too
//	for (unsigned int t = 0; t < 10; t++)
//		if (!do_mh_step(T, proof_prior)) return false;

	free_logical_forms(logical_forms, parse_count);
	return true;
}

template<typename ArticleSource, typename Parser,
	typename Canonicalizer, typename ProofPrior, typename Printer>
bool read_article(
		unsigned int article_name, const ArticleSource& articles, Parser& parser,
		theory<fol_formula, natural_deduction<fol_formula>, Canonicalizer>& T,
		ProofPrior& proof_prior, Printer& printer)
{
	bool article_exists;
	const article& doc = articles.get(article_name, article_exists);
	if (!article_exists) {
		print("read_article ERROR: No such article '", stderr); print(article_name, stderr, printer); print("'.\n", stderr);
		return false;
	}

	for (unsigned int i = 0; i < doc.sentence_count; i++) {
		if (!read_sentence(articles, parser, doc.sentences[i], T, article_name, proof_prior, printer))
			return false;
	}
	return true;
}

#endif /* EXECUTIVE_H_ */
