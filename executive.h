#ifndef EXECUTIVE_H_
#define EXECUTIVE_H_

#include "article.h"
#include "theory.h"

constexpr double PERPLEXITY_THRESHOLD = 10.0;

inline void free_logical_forms(fol_formula** logical_forms, unsigned int count)
{
	for (unsigned int i = 0; i < count; i++) {
		free(*logical_forms[i]);
		if (logical_forms[i]->reference_count == 0)
			free(logical_forms[i]);
	}
}

template<typename ArticleSource, typename Parser>
bool read_sentence(const ArticleSource& articles,
		Parser& parser, const sentence& s, theory& T,
		unsigned int article_name)
{
	unsigned int parse_count;
	fol_formula* logical_forms[2];
	double log_probabilities[2];
	while (true) {
		/* attempt to parse the sentence */
		array<token> unrecognized = array<token>(16);
		if (!parser.parse<2>(s, logical_forms, log_probabilities, parse_count, T, unrecognized)) return false;

		/* read the article on the unrecognized word */
		fol_formula* definition;
		if (unrecognized.length == 0) {
			break;
		} else if (parse_count > 0 && parser.is_definition_of(logical_forms[0], article_name, T, definition)) {
			/* this is a definition */
			free_logical_forms(logical_forms, parse_count);
			bool success = T.add_definition(definition) && parser.add_definition(s, definition);
			free(*definition); if (definition->reference_count == 0) free(definition);
			return success;
		}

		/* find an article in order to learn about the unrecognized word */
		if (!read_article(unrecognized[0].id, articles, parser, T)) {
			free_logical_forms(logical_forms, parse_count);
			return false;
		}
		free_logical_forms(logical_forms, parse_count);
		continue;
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
	if (!T.add_formula(logical_forms[0])) {
		free_logical_forms(logical_forms, parse_count);
		return false;
	}

	free_logical_forms(logical_forms, parse_count);
	return true;
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
