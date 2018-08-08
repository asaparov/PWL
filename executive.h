#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include "article.h"
#include "theory.h"

constexpr double PERPLEXITY_THRESHOLD = 10.0;

template<typename ArticleSource, typename Parser>
bool read_sentence(const ArticleSource& articles,
		Parser& parser, const sentence& s, theory& T)
{
	/* check if the sentence has any unrecognized words */
	while (true) {
		token unrecognized;
		if (!parse.get_unrecognized_word(s, unrecognized));
			break;

		/* read the article on the unrecognized word */
		if (!read_article(unrecognized.id, articles, parser, T))
			return false;
	}

	unsigned int parse_count;
	fol_formula* logical_forms[2];
	double log_probabilities[2];
	if (!parser.parse<2>(s, logical_forms, log_probabilities, parse_count)) return false;

	if (parse_count == 0) {
		fprintf(stderr, "read_sentence ERROR: Given sentence has no valid parses.\n");
		return false;
	} else if (parse_count > 1 && (log_probabilities[0] - log_probabilities[1]) / s.length < PERPLEXITY_THRESHOLD) {
		/* this parse is too ambiguous */
		return true;
	}

	/* add the most probable logical form to the theory */
	if (!T.add_formula(logical_forms[0]))
		return false;

	for (unsigned int i = 0; i < parse_count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i].reference_count == 0)
			free(logical_forms[i]);
	}
}

template<typename ArticleSource, typename Parser>
bool read_article(unsigned int article_name,
		const ArticleSource& articles,
		Parser& parser, theory& T)
{
	const article doc = articles.get(article_name);
	for (unsigned int i = 0; i < doc.sentence_count; i++) {
		if (!read_sentence(articles, parser, doc.sentences[i], T))
			return false;
	}
	return true;
}

#endif /* EXECUTIVE_H_ */
