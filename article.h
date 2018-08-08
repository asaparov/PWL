#ifndef ARTICLE_H_
#define ARTICLE_H_

struct token {
	unsigned int id;
};

struct sentence {
	token* tokens;
	unsigned int length;

	static inline void free() {
		free(tokens);
	}
};

struct article {
	sentence* sentences;
	unsigned int sentence_count;

	~article() {
		for (unsigned int i = 0; i < sentence_count; i++)
			free(sentences[i]);
		free(sentences);
	}
};

#endif /* ARTICLE_H_ */
