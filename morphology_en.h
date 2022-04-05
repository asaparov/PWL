#ifndef MORPHOLOGY_H_
#define MORPHOLOGY_H_

#include <core/map.h>

enum part_of_speech : uint_fast8_t {
	POS_NOUN = 1,
	POS_VERB = 2,
	POS_ADJECTIVE = 3,
	POS_ADVERB = 4,
	POS_OTHER = 5
};

template<typename T>
inline bool merge_arrays(
		T*& dst, unsigned int& dst_length,
		const T* src, unsigned int src_length)
{
	if (dst_length + src_length > 0) {
		T* new_dst = (T*) malloc(sizeof(T) * (dst_length + src_length));
		if (new_dst == nullptr) {
			fprintf(stderr, "merge_arrays ERROR: Out of memory.\n");
			return false;
		}
		unsigned int new_length = 0;
		set_union(new_dst, new_length, dst, dst_length, src, src_length);
		for (unsigned int i = 0; i < dst_length; i++)
			free(dst[i]);
		free(dst);
		dst = new_dst;
		dst_length = new_length;
	}
	return true;
}

inline void try_free(sequence& seq) {
	free(seq);
}

inline void try_free(pair<sequence, sequence>& entry) {
	if (entry.key.tokens != nullptr) free(entry.key);
	if (entry.value.tokens != nullptr) free(entry.value);
}

template<typename T>
unsigned int do_unique(T* array, size_t length)
{
	unsigned int result = 0;
	for (unsigned int i = 1; i < length; i++) {
		if (array[result] != array[i])
			move(array[i], array[++result]);
		else try_free(array[i]);
	}
	return result + 1;
}

inline bool init_string_id_array(
		sequence*& dst, unsigned int& dst_length,
		const string* strs, unsigned int str_count,
		hash_map<string, unsigned int>& names)
{
	dst_length = str_count;
	if (str_count == 0) {
		dst = nullptr;
	} else {
		dst = (sequence*) malloc(sizeof(sequence) * str_count);
		if (dst == nullptr) {
			fprintf(stderr, "init_string_id_array ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str_count; i++) {
			array<unsigned int> tokens(2);
			if (!tokenize(strs[i].data, strs[i].length, tokens, names)) {
				for (unsigned int j = 0; j < i; j++)
					free(dst[j]);
				free(dst); return false;
			}

			sequence src(tokens.data, (unsigned int) tokens.length);
			dst[i] = src;
		}

		/* sort the IDs */
		insertion_sort(dst, dst_length, default_sorter());
		dst_length = do_unique(dst, dst_length);
	}
	return true;
}

enum class properness {
	IMPROPER,
	PROPER,
	BOTH
};

enum class countability {
	COUNTABLE,
	UNCOUNTABLE,
	BOTH,
	PLURAL_UNATTESTED,
	PLURAL_UNKNOWN,
	PLURAL_ONLY
};

struct noun_root {
	countability count;
	properness is_proper;
	sequence* plural;
	unsigned int plural_count;

	inline bool merge(const noun_root& src) {
		/* first merge countabilities */
		switch (count) {
		case countability::COUNTABLE:
			if (src.count == countability::UNCOUNTABLE || src.count == countability::BOTH)
				count = countability::BOTH;
			break;
		case countability::UNCOUNTABLE:
			if (src.count == countability::COUNTABLE) {
				count = countability::BOTH;
			} else if (src.count == countability::PLURAL_ONLY) {
				fprintf(stderr, "noun_root.merge ERROR: PLURAL_ONLY can only be merged with PLURAL_ONLY or COUNTABLE.\n");
				return false;
			} else {
				count = src.count;
			}
			break;
		case countability::BOTH:
			if (src.count == countability::PLURAL_ONLY) {
				fprintf(stderr, "noun_root.merge ERROR: PLURAL_ONLY can only be merged with PLURAL_ONLY or COUNTABLE.\n");
				return false;
			}
			break;
		case countability::PLURAL_UNATTESTED:
			if (src.count == countability::PLURAL_ONLY) {
				fprintf(stderr, "noun_root.merge ERROR: PLURAL_ONLY can only be merged with PLURAL_ONLY or COUNTABLE.\n");
				return false;
			} else if (src.count != countability::UNCOUNTABLE && src.count != countability::PLURAL_UNKNOWN) {
				count = src.count;
			}
			break;
		case countability::PLURAL_UNKNOWN:
			if (src.count == countability::PLURAL_ONLY) {
				fprintf(stderr, "noun_root.merge ERROR: PLURAL_ONLY can only be merged with PLURAL_ONLY or COUNTABLE.\n");
				return false;
			} else if (src.count != countability::UNCOUNTABLE) {
				count = src.count;
			}
			break;
		case countability::PLURAL_ONLY:
			if (src.count == countability::COUNTABLE) {
				count = countability::COUNTABLE;
			} else if (src.count != countability::PLURAL_ONLY) {
				fprintf(stderr, "noun_root.merge ERROR: PLURAL_ONLY can only be merged with PLURAL_ONLY or COUNTABLE.\n");
				return false;
			}
			break;
		}

		/* first merge properness */
		switch (is_proper) {
		case properness::BOTH:
			break;
		case properness::IMPROPER:
			if (src.is_proper == properness::PROPER)
				is_proper = properness::BOTH;
			break;
		case properness::PROPER:
			if (src.is_proper == properness::IMPROPER)
				is_proper = properness::BOTH;
			break;
		}

		return merge_arrays(plural, plural_count, src.plural, src.plural_count);
	}

	static inline void move(const noun_root& src, noun_root& dst) {
		dst.count = src.count;
		dst.is_proper = src.is_proper;
		core::move(src.plural, dst.plural);
		dst.plural_count = src.plural_count;
	}

	static inline void free(noun_root& root) {
		if (root.plural_count != 0) {
			for (unsigned int i = 0; i < root.plural_count; i++)
				core::free(root.plural[i]);
			core::free(root.plural);
		}
	}
};

inline bool init(noun_root& root, bool is_proper, countability count,
		const string* plural_forms, unsigned int plural_form_count,
		hash_map<string, unsigned int>& names)
{
	root.count = count;
	root.is_proper = (is_proper ? properness::PROPER : properness::IMPROPER);
	return init_string_id_array(root.plural, root.plural_count, plural_forms, plural_form_count, names);
}

inline bool init(noun_root& dst, const noun_root& src) {
	dst.count = src.count;
	dst.is_proper = src.is_proper;
	dst.plural_count = src.plural_count;
	if (src.plural_count == 0) {
		dst.plural = nullptr;
	} else {
		dst.plural = (sequence*) malloc(sizeof(sequence) * src.plural_count);
		if (dst.plural == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `noun_root.plural`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.plural_count; i++) {
			if (!init(dst.plural[i], src.plural[i])) {
				for (unsigned int j = 0; j < i; j++) free(dst.plural[j]);
				free(dst.plural); return false;
			}
		}
	}
	return true;
}

enum class comparability {
	COMPARABLE,
	INCOMPARABLE,
	UNKNOWN,
	COMPARABLE_ONLY
};

struct adjective_root {
	comparability comp;
	pair<sequence, sequence>* inflected_forms; /* comparatives and superlatives, respectively */
	unsigned int inflected_form_count;
	sequence adj_root;

	inline bool merge(const adjective_root& src) {
		/* first merge comparabilities */
		switch (comp) {
		case comparability::COMPARABLE:
			break;
		case comparability::INCOMPARABLE:
			if (src.comp == comparability::COMPARABLE || src.comp == comparability::UNKNOWN)
				comp = src.comp;
			break;
		case comparability::UNKNOWN:
			comp = src.comp;
			break;
		case comparability::COMPARABLE_ONLY:
			if (src.comp == comparability::COMPARABLE) {
				comp = src.comp;
			} else if (src.comp == comparability::INCOMPARABLE) {
				fprintf(stderr, "adjective_root.merge ERROR: INCOMPARABLE cannot be merged with COMPARABLE_ONLY.\n");
				return false;
			}
			break;
		}

		return merge_arrays(inflected_forms, inflected_form_count, src.inflected_forms, src.inflected_form_count);
	}

	static inline void move(const adjective_root& src, adjective_root& dst) {
		dst.comp = src.comp;
		core::move(src.inflected_forms, dst.inflected_forms);
		dst.inflected_form_count = src.inflected_form_count;
		core::move(src.adj_root, dst.adj_root);
	}

	static inline void free(adjective_root& root) {
		if (root.inflected_form_count != 0) {
			for (unsigned int i = 0; i < root.inflected_form_count; i++) {
				if (root.inflected_forms[i].key.tokens != nullptr)
					core::free(root.inflected_forms[i].key);
				if (root.inflected_forms[i].value.tokens != nullptr)
					core::free(root.inflected_forms[i].value);
			}
			core::free(root.inflected_forms);
		}
		if (root.adj_root.length != 0)
			core::free(root.adj_root);
	}
};

typedef adjective_root adverb_root;

inline bool init(adjective_root& root, comparability comp,
		const pair<string, string>* inflected_forms,
		unsigned int inflected_form_count,
		const sequence& adj_root,
		hash_map<string, unsigned int>& names)
{
	root.comp = comp;
	if (adj_root.length == 0) {
		root.adj_root.tokens = nullptr;
		root.adj_root.length = 0;
	} else if (!init(root.adj_root, adj_root)) {
		return false;
	}

	root.inflected_form_count = inflected_form_count;
	if (inflected_form_count == 0) {
		root.inflected_forms = nullptr;
	} else {
		root.inflected_forms = (pair<sequence, sequence>*) malloc(sizeof(pair<sequence, sequence>) * inflected_form_count);
		if (root.inflected_forms == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `adjective_root.inflected_forms`.\n");
			if (root.adj_root.length != 0) free(root.adj_root);
			return false;
		}
		for (unsigned int i = 0; i < inflected_form_count; i++) {
			if (inflected_forms[i].key.length == 0) {
				root.inflected_forms[i].key = {nullptr, 0};
			} else {
				array<unsigned int> tokens(2);
				if (!tokenize(inflected_forms[i].key.data, inflected_forms[i].key.length, tokens, names)) {
					for (unsigned int j = 0; j < i; j++) {
						if (root.inflected_forms[i].key.tokens != nullptr)
							free(root.inflected_forms[i].key);
						if (root.inflected_forms[i].value.tokens != nullptr)
							free(root.inflected_forms[i].value);
					}
					free(root.inflected_forms);
					if (root.adj_root.length != 0) free(root.adj_root);
					return false;
				}
				sequence src(tokens.data, (unsigned int) tokens.length);
				root.inflected_forms[i].key = src;
			}

			if (inflected_forms[i].value.length == 0) {
				root.inflected_forms[i].value = {nullptr, 0};
			} else {
				array<unsigned int> tokens(2);
				if (!tokenize(inflected_forms[i].value.data, inflected_forms[i].value.length, tokens, names)) {
					for (unsigned int j = 0; j < i; j++) {
						if (root.inflected_forms[i].key.tokens != nullptr)
							free(root.inflected_forms[i].key);
						if (root.inflected_forms[i].value.tokens != nullptr)
							free(root.inflected_forms[i].value);
					}
					free(root.inflected_forms);
					if (root.adj_root.length != 0) free(root.adj_root);
					return false;
				}
				sequence src(tokens.data, (unsigned int) tokens.length);
				root.inflected_forms[i].value = src;
			}
		}

		/* remove any entries where both the comparative and superlative are empty */
		for (unsigned int i = 0; i < root.inflected_form_count; i++) {
			if (root.inflected_forms[i].key.tokens == nullptr && root.inflected_forms[i].value.tokens == nullptr) {
				move(root.inflected_forms[root.inflected_form_count - 1], root.inflected_forms[i]);
				root.inflected_form_count--;
				i--;
			}
		}

		if (root.inflected_form_count == 0) {
			free(root.inflected_forms);
			root.inflected_forms = nullptr;
		} else {
			/* sort the inflected forms */
			insertion_sort(root.inflected_forms, root.inflected_form_count, default_sorter());
			root.inflected_form_count = do_unique(root.inflected_forms, root.inflected_form_count);
		}
	}
	return true;
}

bool init(adjective_root& dst, const adjective_root& src) {
	dst.comp = src.comp;
	dst.inflected_form_count = src.inflected_form_count;
	if (src.inflected_form_count != 0) {
		dst.inflected_forms = (pair<sequence, sequence>*) malloc(sizeof(pair<sequence, sequence>) * src.inflected_form_count);
		if (dst.inflected_forms == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `adjective_root.inflected_form_count`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.inflected_form_count; i++) {
			if (!init(dst.inflected_forms[i].key, src.inflected_forms[i].key)) {
				for (unsigned int j = 0; j < i; j++) { free(dst.inflected_forms[j].key); free(dst.inflected_forms[j].value); }
				free(dst.inflected_forms); return false;
			} else if (!init(dst.inflected_forms[i].value, src.inflected_forms[i].value)) {
				free(dst.inflected_forms[i].key);
				for (unsigned int j = 0; j < i; j++) { free(dst.inflected_forms[j].key); free(dst.inflected_forms[j].value); }
				free(dst.inflected_forms); return false;
			}
		}
	} if (src.adj_root.length == 0) {
		dst.adj_root.tokens = nullptr;
		dst.adj_root.length = 0;
	} else if (!init(dst.adj_root, src.adj_root)) {
		for (unsigned int j = 0; j < src.inflected_form_count; j++) { free(dst.inflected_forms[j].key); free(dst.inflected_forms[j].value); }
		free(dst.inflected_forms); return false;
	}
	return true;
}

struct verb_root {
	sequence* present_3sg;
	unsigned int present_3sg_count;
	sequence* present_participle;
	unsigned int present_participle_count;
	sequence* simple_past;
	unsigned int simple_past_count;
	sequence* past_participle;
	unsigned int past_participle_count;

	inline bool merge(const verb_root& src) {
		return merge_arrays(present_3sg, present_3sg_count, src.present_3sg, src.present_3sg_count)
			&& merge_arrays(present_participle, present_participle_count, src.present_participle, src.present_participle_count)
			&& merge_arrays(simple_past, simple_past_count, src.simple_past, src.simple_past_count)
			&& merge_arrays(past_participle, past_participle_count, src.past_participle, src.past_participle_count);
	}

	static inline void move(const verb_root& src, verb_root& dst) {
		core::move(src.present_3sg, dst.present_3sg);
		dst.present_3sg_count = src.present_3sg_count;
		core::move(src.present_participle, dst.present_participle);
		dst.present_participle_count = src.present_participle_count;
		core::move(src.simple_past, dst.simple_past);
		dst.simple_past_count = src.simple_past_count;
		core::move(src.past_participle, dst.past_participle);
		dst.past_participle_count = src.past_participle_count;
	}

	static inline void free(verb_root& root) {
		if (root.present_3sg_count != 0) {
			for (unsigned int i = 0; i < root.present_3sg_count; i++)
				core::free(root.present_3sg[i]);
			core::free(root.present_3sg);
		} if (root.present_participle_count != 0) {
			for (unsigned int i = 0; i < root.present_participle_count; i++)
				core::free(root.present_participle[i]);
			core::free(root.present_participle);
		} if (root.simple_past_count != 0) {
			for (unsigned int i = 0; i < root.simple_past_count; i++)
				core::free(root.simple_past[i]);
			core::free(root.simple_past);
		} if (root.past_participle_count != 0) {
			for (unsigned int i = 0; i < root.past_participle_count; i++)
				core::free(root.past_participle[i]);
			core::free(root.past_participle);
		}
	}
};

inline bool init(verb_root& root,
		const string* present_3sg,
		unsigned int present_3sg_count,
		const string* present_participle,
		unsigned int present_participle_count,
		const string* simple_past,
		unsigned int simple_past_count,
		const string* past_participle,
		unsigned int past_participle_count,
		hash_map<string, unsigned int>& names)
{
	if (!init_string_id_array(root.present_3sg, root.present_3sg_count, present_3sg, present_3sg_count, names)) {
		return false;
	} else if (!init_string_id_array(root.present_participle, root.present_participle_count, present_participle, present_participle_count, names)) {
		free(root.present_3sg);
		return false;
	} else if (!init_string_id_array(root.simple_past, root.simple_past_count, simple_past, simple_past_count, names)) {
		free(root.present_3sg);
		free(root.present_participle);
		return false;
	} else if (!init_string_id_array(root.past_participle, root.past_participle_count, past_participle, past_participle_count, names)) {
		free(root.present_3sg);
		free(root.present_participle);
		free(root.simple_past);
		return false;
	}
	return true;
}

inline bool init(verb_root& dst, const verb_root& src) {
	dst.present_3sg_count = src.present_3sg_count;
	dst.present_participle_count = src.present_participle_count;
	dst.simple_past_count = src.simple_past_count;
	dst.past_participle_count = src.past_participle_count;
	if (src.present_3sg_count == 0) {
		dst.present_3sg = nullptr;
	} else {
		dst.present_3sg = (sequence*) malloc(sizeof(sequence) * src.present_3sg_count);
		if (dst.present_3sg == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `verb_root.present_3sg`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.present_3sg_count; i++) {
			if (!init(dst.present_3sg[i], src.present_3sg[i])) {
				for (unsigned int j = 0; j < i; j++) free(dst.present_3sg[j]);
				free(dst.present_3sg); return false;
			}
		}
	}
	if (src.present_participle_count == 0) {
		dst.present_participle = nullptr;
	} else {
		dst.present_participle = (sequence*) malloc(sizeof(sequence) * src.present_participle_count);
		if (dst.present_participle == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `verb_root.present_participle`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.present_participle_count; i++) {
			if (!init(dst.present_participle[i], src.present_participle[i])) {
				if (dst.present_3sg_count != 0) {
					for (unsigned int j = 0; j < dst.present_3sg_count; j++) free(dst.present_3sg[j]);
					free(dst.present_3sg);
				}
				for (unsigned int j = 0; j < i; j++) free(dst.present_participle[j]);
				free(dst.present_participle); return false;
			}
		}
	}
	if (src.simple_past_count == 0) {
		dst.simple_past = nullptr;
	} else {
		dst.simple_past = (sequence*) malloc(sizeof(sequence) * src.simple_past_count);
		if (dst.simple_past == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `verb_root.simple_past`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.simple_past_count; i++) {
			if (!init(dst.simple_past[i], src.simple_past[i])) {
				if (dst.present_3sg_count != 0) {
					for (unsigned int j = 0; j < dst.present_3sg_count; j++) free(dst.present_3sg[j]);
					free(dst.present_3sg);
				} if (dst.present_participle_count != 0) {
					for (unsigned int j = 0; j < dst.present_participle_count; j++) free(dst.present_participle[j]);
					free(dst.present_participle);
				}
				for (unsigned int j = 0; j < i; j++) free(dst.simple_past[j]);
				free(dst.simple_past); return false;
			}
		}
	}
	if (src.past_participle_count == 0) {
		dst.past_participle = nullptr;
	} else {
		dst.past_participle = (sequence*) malloc(sizeof(sequence) * src.past_participle_count);
		if (dst.past_participle == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `verb_root.past_participle`.\n");
			return false;
		}
		for (unsigned int i = 0; i < src.past_participle_count; i++) {
			if (!init(dst.past_participle[i], src.past_participle[i])) {
				if (dst.present_3sg_count != 0) {
					for (unsigned int j = 0; j < dst.present_3sg_count; j++) free(dst.present_3sg[j]);
					free(dst.present_3sg);
				} if (dst.present_participle_count != 0) {
					for (unsigned int j = 0; j < dst.present_participle_count; j++) free(dst.present_participle[j]);
					free(dst.present_participle);
				} if (dst.simple_past_count != 0) {
					for (unsigned int j = 0; j < dst.simple_past_count; j++) free(dst.simple_past[j]);
					free(dst.simple_past);
				}
				for (unsigned int j = 0; j < i; j++) free(dst.past_participle[j]);
				free(dst.past_participle); return false;
			}
		}
	}
	return true;
}

enum class grammatical_person : uint_fast8_t {
	FIRST,
	SECOND,
	THIRD,

	ANY,
	FIRST_OR_SECOND,
	FIRST_OR_THIRD
};

enum class grammatical_num : uint_fast8_t {
	NONE = 0,
	SINGULAR,
	PLURAL,

	ANY,
	SINGULAR_OR_NONE,
	PLURAL_OR_NONE
};

enum class grammatical_mood : uint_fast8_t {
	INDICATIVE = 0,
	PRESENT_PARTICIPLE,
	PAST_PARTICIPLE,
	BARE_INFINITIVE,
	TO_INFINITIVE,
	SUBJUNCTIVE,

	ANY,
	NOT_TO_INFINITIVE,
	NOT_SUBJUNCTIVE,
	NOT_TO_INF_OR_SUBJ
};

enum class grammatical_comparison : uint_fast8_t {
	NONE,
	COMPARATIVE,
	SUPERLATIVE,

	ANY,
	NOT_COMPARATIVE
};

enum class grammatical_tense : uint_fast8_t {
	PRESENT,
	PAST,
	ANY
};

inline bool intersect(grammatical_person& out, grammatical_person first, grammatical_person second) {
	if (first == grammatical_person::ANY) {
		out = second;
		return true;
	} else if (second == grammatical_person::ANY) {
		out = first;
		return true;
	} else if (first == grammatical_person::FIRST_OR_SECOND) {
		if (second == grammatical_person::THIRD) {
			return false;
		} else if (second == grammatical_person::FIRST_OR_THIRD) {
			out = grammatical_person::FIRST;
		} else {
			out = second;
		}
		return true;
	} else if (second == grammatical_person::FIRST_OR_SECOND) {
		if (first == grammatical_person::THIRD) {
			return false;
		} else if (first == grammatical_person::FIRST_OR_THIRD) {
			out = grammatical_person::FIRST;
		} else {
			out = first;
		}
		return true;
	} else if (first == grammatical_person::FIRST_OR_THIRD) {
		if (second == grammatical_person::SECOND) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_person::FIRST_OR_THIRD) {
		if (first == grammatical_person::SECOND) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == second) {
		return true;
	} else {
		return false;
	}
}

inline bool has_intersection(grammatical_person first, grammatical_person second) {
	grammatical_person dummy;
	return intersect(dummy, first, second);
}

static inline bool intersect(grammatical_num& out, grammatical_num first, grammatical_num second)
{
	if (first == grammatical_num::ANY) {
		out = second;
		return true;
	} else if (second == grammatical_num::ANY) {
		out = first;
		return true;
	} else if (first == grammatical_num::SINGULAR_OR_NONE) {
		if (second == grammatical_num::PLURAL) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_num::SINGULAR_OR_NONE) {
		if (first == grammatical_num::PLURAL) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == grammatical_num::PLURAL_OR_NONE) {
		if (second == grammatical_num::SINGULAR) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_num::PLURAL_OR_NONE) {
		if (first == grammatical_num::SINGULAR) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == second) {
		out = first;
		return true;
	} else {
		return false;
	}
}

static inline bool has_intersection(grammatical_num first, grammatical_num second) {
	grammatical_num dummy;
	return intersect(dummy, first, second);
}

inline bool intersect(grammatical_mood& out, grammatical_mood first, grammatical_mood second)
{
	if (first == grammatical_mood::ANY) {
		out = second;
		return true;
	} else if (second == grammatical_mood::ANY) {
		out = first;
		return true;
	} else if (first == grammatical_mood::NOT_TO_INFINITIVE) {
		if (second == grammatical_mood::TO_INFINITIVE) {
			return false;
		} else if (second == grammatical_mood::NOT_SUBJUNCTIVE) {
			out = grammatical_mood::NOT_TO_INF_OR_SUBJ;
		} else {
			out = second;
		}
		return true;
	} else if (second == grammatical_mood::NOT_TO_INFINITIVE) {
		if (first == grammatical_mood::TO_INFINITIVE) {
			return false;
		} else if (first == grammatical_mood::NOT_SUBJUNCTIVE) {
			out = grammatical_mood::NOT_TO_INF_OR_SUBJ;
		} else {
			out = first;
		}
		return true;
	} else if (first == grammatical_mood::NOT_SUBJUNCTIVE) {
		if (second == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_mood::NOT_SUBJUNCTIVE) {
		if (first == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == grammatical_mood::NOT_TO_INF_OR_SUBJ) {
		if (second == grammatical_mood::TO_INFINITIVE || second == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_mood::NOT_TO_INF_OR_SUBJ) {
		if (first == grammatical_mood::TO_INFINITIVE || first == grammatical_mood::SUBJUNCTIVE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first != second) {
		return false;
	} else {
		out = first;
		return true;
	}
}

inline bool has_intersection(grammatical_mood first, grammatical_mood second) {
	grammatical_mood dummy;
	return intersect(dummy, first, second);
}

inline bool intersect(grammatical_comparison& out, grammatical_comparison first, grammatical_comparison second)
{
	if (first == grammatical_comparison::ANY) {
		out = second;
	} else if (second == grammatical_comparison::ANY) {
		out = first;
	} else if (first == grammatical_comparison::NOT_COMPARATIVE) {
		if (second == grammatical_comparison::COMPARATIVE) {
			return false;
		} else {
			out = second;
			return true;
		}
	} else if (second == grammatical_comparison::NOT_COMPARATIVE) {
		if (first == grammatical_comparison::COMPARATIVE) {
			return false;
		} else {
			out = first;
			return true;
		}
	} else if (first == second) {
		out = first;
	} else {
		return false;
	}
	return true;
}

inline bool has_intersection(grammatical_comparison first, grammatical_comparison second) {
	grammatical_comparison dummy;
	return intersect(dummy, first, second);
}

inline bool intersect(grammatical_tense& out, grammatical_tense first, grammatical_tense second)
{
	if (first == grammatical_tense::ANY) {
		out = second;
	} else if (second == grammatical_tense::ANY) {
		out = first;
	} else if (first == second) {
		out = first;
	} else {
		return false;
	}
	return true;
}

inline bool has_intersection(grammatical_tense first, grammatical_tense second) {
	grammatical_tense dummy;
	return intersect(dummy, first, second);
}

struct inflected_noun {
	sequence root;
	properness is_proper;
	grammatical_num number;

	static inline void free(inflected_noun& key) {
		core::free(key.root);
	}
};

inline bool init(inflected_noun& dst, const inflected_noun& src) {
	if (!init(dst.root, src.root))
		return false;
	dst.is_proper = src.is_proper;
	dst.number = src.number;
	return true;
}

inline bool operator == (const inflected_noun& first, const inflected_noun& second) {
	return first.root == second.root
		&& first.number == second.number;
}

struct inflected_adjective {
	sequence root;
	grammatical_comparison comp;

	static inline void free(inflected_adjective& key) {
		core::free(key.root);
	}
};

inline bool init(inflected_adjective& dst, const inflected_adjective& src) {
	if (!init(dst.root, src.root))
		return false;
	dst.comp = src.comp;
	return true;
}

inline bool operator == (const inflected_adjective& first, const inflected_adjective& second) {
	return first.root == second.root
		&& first.comp == second.comp;
}

typedef inflected_adjective inflected_adverb;

struct inflected_verb {
	sequence root;
	grammatical_person person;
	grammatical_num number;
	grammatical_mood mood;
	grammatical_tense tense;

	static inline void free(inflected_verb& key) {
		core::free(key.root);
	}
};

inline bool init(inflected_verb& dst, const inflected_verb& src) {
	if (!init(dst.root, src.root))
		return false;
	dst.person = src.person;
	dst.number = src.number;
	dst.mood = src.mood;
	dst.tense = src.tense;
	return true;
}

inline bool operator == (const inflected_verb& first, const inflected_verb& second) {
	return first.root == second.root
		&& first.person == second.person
		&& first.number == second.number
		&& first.mood == second.mood
		&& first.tense == second.tense;
}

struct morphology_en {
	hash_map<sequence, noun_root> nouns;
	hash_map<sequence, adjective_root> adjectives;
	hash_map<sequence, adverb_root> adverbs;
	hash_map<sequence, verb_root> verbs;

	hash_map<sequence, array<inflected_noun>> inflected_nouns;
	hash_map<sequence, array<inflected_adjective>> inflected_adjectives;
	hash_map<sequence, array<inflected_adverb>> inflected_adverbs;
	hash_map<sequence, array<inflected_verb>> inflected_verbs;

	hash_map<unsigned int, unsigned int> capitalization_map;
	hash_map<unsigned int, unsigned int> decapitalization_map;

	hash_map<sequence, array<sequence>> adjective_adverb_map;

	unsigned int MORE_COMPARATIVE_ID;
	unsigned int MOST_SUPERLATIVE_ID;
	unsigned int FURTHER_COMPARATIVE_ID;
	unsigned int FURTHEST_SUPERLATIVE_ID;

	sequence BE_SEQ, AM_SEQ, ARE_SEQ, IS_SEQ, WAS_SEQ, WERE_SEQ, BEING_SEQ, BEEN_SEQ;

	static string MORE_COMPARATIVE_STRING;
	static string MOST_SUPERLATIVE_STRING;
	static string FURTHER_COMPARATIVE_STRING;
	static string FURTHEST_SUPERLATIVE_STRING;

	morphology_en() :
		nouns(1024), adjectives(1024), adverbs(1024), verbs(1024), inflected_nouns(2048),
		inflected_adjectives(2048), inflected_adverbs(2048), inflected_verbs(2048),
		capitalization_map(2048), decapitalization_map(2048), adjective_adverb_map(1024),
		BE_SEQ(nullptr, 0), AM_SEQ(nullptr, 0), ARE_SEQ(nullptr, 0), IS_SEQ(nullptr, 0),
		WAS_SEQ(nullptr, 0), WERE_SEQ(nullptr, 0), BEING_SEQ(nullptr, 0), BEEN_SEQ(nullptr, 0)
	{ }

	~morphology_en() { free_helper(); }

	static inline void free(morphology_en& morph) {
		morph.free_helper();
		core::free(morph.nouns);
		core::free(morph.adjectives);
		core::free(morph.adverbs);
		core::free(morph.verbs);
		core::free(morph.inflected_nouns);
		core::free(morph.inflected_adjectives);
		core::free(morph.inflected_adverbs);
		core::free(morph.inflected_verbs);
		core::free(morph.capitalization_map);
		core::free(morph.decapitalization_map);
		core::free(morph.adjective_adverb_map);
	}

	inline bool initialize(hash_map<string, unsigned int>& names) {
		/* add highly irregular verbs */
		unsigned int be, am, are, is, was, were, being, been;
		if (!get_token("be", be, names)
		 || !get_token("am", am, names)
		 || !get_token("are", are, names)
		 || !get_token("is", is, names)
		 || !get_token("was", was, names)
		 || !get_token("were", were, names)
		 || !get_token("being", being, names)
		 || !get_token("been", been, names)
		 || !init(BE_SEQ, sequence(&be, 1))
		 || !init(AM_SEQ, sequence(&am, 1))
		 || !init(ARE_SEQ, sequence(&are, 1))
		 || !init(IS_SEQ, sequence(&is, 1))
		 || !init(WAS_SEQ, sequence(&was, 1))
		 || !init(WERE_SEQ, sequence(&were, 1))
		 || !init(BEING_SEQ, sequence(&being, 1))
		 || !init(BEEN_SEQ, sequence(&been, 1)))
			return false;

		/* add the infinitive form */
		if (!add_inflected_form(inflected_verbs, BE_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::BARE_INFINITIVE, grammatical_tense::ANY}))
			return false;

		/* add the subjunctive forms */
		if (!add_inflected_form(inflected_verbs, BE_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::SUBJUNCTIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, WERE_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::SUBJUNCTIVE, grammatical_tense::PAST}))
			return false;

		/* add the indicative verb forms */
		if (!add_inflected_form(inflected_verbs, AM_SEQ, {BE_SEQ, grammatical_person::FIRST, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, ARE_SEQ, {BE_SEQ, grammatical_person::SECOND, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, IS_SEQ, {BE_SEQ, grammatical_person::THIRD, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, WAS_SEQ, {BE_SEQ, grammatical_person::FIRST_OR_THIRD, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PAST})
		 || !add_inflected_form(inflected_verbs, WERE_SEQ, {BE_SEQ, grammatical_person::SECOND, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PAST})
		 || !add_inflected_form(inflected_verbs, ARE_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::PLURAL, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, WERE_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::PLURAL, grammatical_mood::INDICATIVE, grammatical_tense::PAST}))
			return false;

		/* add the participle forms */
		if (!add_inflected_form(inflected_verbs, BEING_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PRESENT_PARTICIPLE, grammatical_tense::ANY})
		 || !add_inflected_form(inflected_verbs, BEEN_SEQ, {BE_SEQ, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PAST_PARTICIPLE, grammatical_tense::ANY}))
			return false;

		return true;
	}

	inline bool inflect_verb(const inflected_verb& verb, array<sequence>& inflections) const
	{
		if (verb.root == BE_SEQ) {
			/* handle the highly irregular verb "be" */
			if (has_intersection(verb.mood, grammatical_mood::BARE_INFINITIVE) || has_intersection(verb.mood, grammatical_mood::SUBJUNCTIVE)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], BE_SEQ))
					return false;
				inflections.length++;
			}

			if ((has_intersection(verb.mood, grammatical_mood::SUBJUNCTIVE) && has_intersection(verb.tense, grammatical_tense::PAST))
			 || (has_intersection(verb.person, grammatical_person::SECOND) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PAST))
			 || (has_intersection(verb.number, grammatical_num::PLURAL) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PAST)))
			{
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], WERE_SEQ))
					return false;
				inflections.length++;
			}

			if (has_intersection(verb.person, grammatical_person::FIRST) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], AM_SEQ))
					return false;
				inflections.length++;
			}

			if ((has_intersection(verb.person, grammatical_person::SECOND) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT))
			 || (has_intersection(verb.number, grammatical_num::PLURAL) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT)))
			{
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], ARE_SEQ))
					return false;
				inflections.length++;
			}

			if (has_intersection(verb.person, grammatical_person::THIRD) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], IS_SEQ))
					return false;
				inflections.length++;
			}

			if (has_intersection(verb.person, grammatical_person::FIRST_OR_THIRD) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PAST)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], WAS_SEQ))
					return false;
				inflections.length++;
			}

			if (has_intersection(verb.mood, grammatical_mood::PRESENT_PARTICIPLE)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], BEING_SEQ))
					return false;
				inflections.length++;
			}

			if (has_intersection(verb.mood, grammatical_mood::PAST_PARTICIPLE)) {
				if (!inflections.ensure_capacity(inflections.length + 1)
				 || !init(inflections[inflections.length], BEEN_SEQ))
					return false;
				inflections.length++;
			}
			return true;
		}

		/* handle (mostly) regular verbs */
		bool contains;
		const verb_root& root = verbs.get(verb.root, contains);
		if (!contains) return false;
		if (has_intersection(verb.mood, grammatical_mood::BARE_INFINITIVE) || has_intersection(verb.mood, grammatical_mood::SUBJUNCTIVE)
		 || (has_intersection(verb.person, grammatical_person::FIRST_OR_SECOND) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT))
		 || (has_intersection(verb.number, grammatical_num::PLURAL) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT)))
		{
			if (!inflections.ensure_capacity(inflections.length + 1)
			 || !init(inflections[inflections.length], verb.root))
				return false;
			inflections.length++;
		}

		if (has_intersection(verb.person, grammatical_person::THIRD) && has_intersection(verb.number, grammatical_num::SINGULAR) && has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PRESENT))
		{
			if (!inflections.ensure_capacity(inflections.length + root.present_3sg_count))
				return false;
			for (unsigned int i = 0; i < root.present_3sg_count; i++) {
				if (!init(inflections[inflections.length], root.present_3sg[i]))
					return false;
				inflections.length++;
			}
		}

		if (has_intersection(verb.mood, grammatical_mood::INDICATIVE) && has_intersection(verb.tense, grammatical_tense::PAST))
		{
			if (!inflections.ensure_capacity(inflections.length + root.simple_past_count))
				return false;
			for (unsigned int i = 0; i < root.simple_past_count; i++) {
				if (!init(inflections[inflections.length], root.simple_past[i]))
					return false;
				inflections.length++;
			}
		}

		if (has_intersection(verb.mood, grammatical_mood::PRESENT_PARTICIPLE))
		{
			if (!inflections.ensure_capacity(inflections.length + root.present_participle_count))
				return false;
			for (unsigned int i = 0; i < root.present_participle_count; i++) {
				if (!init(inflections[inflections.length], root.present_participle[i]))
					return false;
				inflections.length++;
			}
		}

		if (has_intersection(verb.mood, grammatical_mood::PRESENT_PARTICIPLE))
		{
			if (!inflections.ensure_capacity(inflections.length + root.past_participle_count))
				return false;
			for (unsigned int i = 0; i < root.past_participle_count; i++) {
				if (!init(inflections[inflections.length], root.past_participle[i]))
					return false;
				inflections.length++;
			}
		}
		return true;
	}

	inline bool inflect_noun(const inflected_noun& noun, array<sequence>& inflections) const
	{
		bool contains;
		const noun_root& root = nouns.get(noun.root, contains);
		if (!contains) return false;

		if (has_intersection(noun.number, grammatical_num::SINGULAR)) {
			if (!inflections.ensure_capacity(inflections.length + 1)
			 || !init(inflections[inflections.length], noun.root))
				return false;
			inflections.length++;
		}

		if (has_intersection(noun.number, grammatical_num::PLURAL)) {
			if (!inflections.ensure_capacity(inflections.length + root.plural_count))
				return false;
			for (unsigned int i = 0; i < root.plural_count; i++) {
				if (!init(inflections[inflections.length], root.plural[i]))
					return false;
				inflections.length++;
			}
		}
		return true;
	}

	inline bool inflect_adjective(const inflected_adjective& adjective, array<sequence>& inflections) const
	{
		bool contains;
		const adjective_root& root = adjectives.get(adjective.root, contains);
		if (!contains) return false;

		if (has_intersection(adjective.comp, grammatical_comparison::NONE)) {
			if (!inflections.ensure_capacity(inflections.length + 1)
			 || !init(inflections[inflections.length], adjective.root))
				return false;
			inflections.length++;
		}

		if (has_intersection(adjective.comp, grammatical_comparison::COMPARATIVE)) {
			if (!inflections.ensure_capacity(inflections.length + root.inflected_form_count))
				return false;
			for (unsigned int i = 0; i < root.inflected_form_count; i++) {
				if (root.inflected_forms[i].key.tokens == nullptr) continue;
				if (!init(inflections[inflections.length], root.inflected_forms[i].key))
					return false;
				inflections.length++;
			}
		}

		if (has_intersection(adjective.comp, grammatical_comparison::SUPERLATIVE)) {
			if (!inflections.ensure_capacity(inflections.length + root.inflected_form_count))
				return false;
			for (unsigned int i = 0; i < root.inflected_form_count; i++) {
				if (root.inflected_forms[i].value.tokens == nullptr) continue;
				if (!init(inflections[inflections.length], root.inflected_forms[i].value))
					return false;
				inflections.length++;
			}
		}
		return true;
	}

	inline bool inflect_adverb(const inflected_adverb& adverb, array<sequence>& inflections) const
	{
		bool contains;
		const adverb_root& root = adverbs.get(adverb.root, contains);
		if (!contains) return false;

		if (has_intersection(adverb.comp, grammatical_comparison::NONE)) {
			if (!inflections.ensure_capacity(inflections.length + 1)
			 || !init(inflections[inflections.length], adverb.root))
				return false;
			inflections.length++;
		}

		if (has_intersection(adverb.comp, grammatical_comparison::COMPARATIVE)) {
			if (!inflections.ensure_capacity(inflections.length + root.inflected_form_count))
				return false;
			for (unsigned int i = 0; i < root.inflected_form_count; i++) {
				if (root.inflected_forms[i].key.tokens == nullptr) continue;
				if (!init(inflections[inflections.length], root.inflected_forms[i].key))
					return false;
				inflections.length++;
			}
		}

		if (has_intersection(adverb.comp, grammatical_comparison::SUPERLATIVE)) {
			if (!inflections.ensure_capacity(inflections.length + root.inflected_form_count))
				return false;
			for (unsigned int i = 0; i < root.inflected_form_count; i++) {
				if (root.inflected_forms[i].value.tokens == nullptr) continue;
				if (!init(inflections[inflections.length], root.inflected_forms[i].value))
					return false;
				inflections.length++;
			}
		}
		return true;
	}

	inline bool add_capitalized_form(unsigned int decapitalized, unsigned int capitalized) {
		return capitalization_map.put(decapitalized, capitalized)
			&& decapitalization_map.put(capitalized, decapitalized);
	}

	inline bool capitalize(unsigned int word_id, unsigned int& capitalized_word_id) const {
		bool contains;
		capitalized_word_id = capitalization_map.get(word_id, contains);
		return contains;
	}

	inline bool decapitalize(unsigned int word_id, unsigned int& decapitalized_word_id) const {
		bool contains;
		decapitalized_word_id = decapitalization_map.get(word_id, contains);
		return contains;
	}

	inline bool is_capitalized(unsigned int word_id) const {
		return decapitalization_map.table.contains(word_id);
	}

	inline bool is_decapitalized(unsigned int word_id) const {
		return capitalization_map.table.contains(word_id);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_noun_root(const sequence& root_id, noun_root& root) {
		return add_root(nouns, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adjective_root(const sequence& root_id, adjective_root& root) {
		return add_root(adjectives, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adverb_root(const sequence& root_id, adverb_root& root) {
		return add_root(adverbs, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_verb_root(const sequence& root_id, verb_root& root) {
		return add_root(verbs, root_id, root);
	}

	bool add_all_inflected_forms() {
		for (auto entry : nouns)
			if (!add_inflected_forms(inflected_nouns, entry.key, entry.value)) return false;
		for (auto entry : adjectives)
			if (!add_inflected_forms(inflected_adjectives, entry.key, entry.value)) return false;
		for (auto entry : adverbs) {
			if (entry.value.adj_root.length != 0) {
				if (!adjective_adverb_map.check_size()) return false;
				bool contains; unsigned int bucket;
				array<sequence>& adj_roots = adjective_adverb_map.get(entry.value.adj_root, contains, bucket);
				if (!contains) {
					if (!array_init(adj_roots, 2)) {
						return false;
					} else if (!init(adjective_adverb_map.table.keys[bucket], entry.value.adj_root)) {
						core::free(adj_roots);
						return false;
					}
					adjective_adverb_map.table.size++;
				}
				if (!adj_roots.ensure_capacity(adj_roots.length + 1)
				 || !init(adj_roots[adj_roots.length], entry.key))
					return false;
				adj_roots.length++;
			}

			if (!add_inflected_forms(inflected_adverbs, entry.key, entry.value)) return false;
		} for (auto entry : verbs)
			if (!add_inflected_forms(inflected_verbs, entry.key, entry.value)) return false;
		return true;
	}

private:
	template<typename T>
	static inline bool add_root(
			hash_map<sequence, T>& root_map,
			const sequence& root_id, T& root)
	{
		if (!root_map.check_size()) {
			core::free(root);
			return false;
		}

		bool contains; unsigned int bucket;
		T& value = root_map.get(root_id, contains, bucket);
		if (!contains) {
			move(root, value);
			root_map.table.keys[bucket] = root_id;
			root_map.table.size++;
		} else {
			if (!value.merge(root)) {
				core::free(root);
				return false;
			}
			core::free(root);
		}
		return true;
	}

	template<typename T>
	static inline bool add_inflected_form(
			hash_map<sequence, array<T>>& inflected_form_map,
			const sequence& new_inflected_form_id, T new_inflected_form)
	{
		bool contains; unsigned int bucket;
		array<T>& inflected_forms = inflected_form_map.get(new_inflected_form_id, contains, bucket);
		if (!contains) {
			if (!array_init(inflected_forms, 1))
				return false;
			inflected_form_map.table.keys[bucket] = new_inflected_form_id;
			inflected_form_map.table.size++;
		}

		if (!inflected_forms.contains(new_inflected_form))
			return inflected_forms.add(new_inflected_form);
		else return true;
	}

	static inline bool add_inflected_forms(
			hash_map<sequence, array<inflected_noun>>& inflected_form_map,
			const sequence& root_id, const noun_root& root)
	{
		if (!inflected_form_map.check_size(inflected_form_map.table.size + root.plural_count + 1))
			return false;

		/* first add the singular form */
		if (root.count != countability::PLURAL_ONLY) {
			if (!add_inflected_form(inflected_form_map, root_id, {root_id, root.is_proper, grammatical_num::SINGULAR}))
				return false;
		}

		/* add the plural forms */
		for (unsigned int i = 0; i < root.plural_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.plural[i], {root_id, root.is_proper, grammatical_num::PLURAL}))
				return false;
		}
		return true;
	}

	static inline bool add_inflected_forms(
			hash_map<sequence, array<inflected_adjective>>& inflected_form_map,
			const sequence& root_id, const adjective_root& root)
	{
		if (!inflected_form_map.check_size(inflected_form_map.table.size + root.inflected_form_count * 2 + 1))
			return false;

		/* first add the base form */
		if (root.comp != comparability::COMPARABLE_ONLY) {
			if (!add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_comparison::NONE}))
				return false;
		}

		/* add the comparative and superlative forms */
		for (unsigned int i = 0; i < root.inflected_form_count; i++) {
			if (root.inflected_forms[i].key.tokens != nullptr) {
				if (!add_inflected_form(inflected_form_map, root.inflected_forms[i].key, {root_id, grammatical_comparison::COMPARATIVE}))
					return false;
			}

			if (root.inflected_forms[i].value.tokens != nullptr) {
				if (!add_inflected_form(inflected_form_map, root.inflected_forms[i].value, {root_id, grammatical_comparison::SUPERLATIVE}))
					return false;
			}
		}
		return true;
	}

	static inline bool add_inflected_forms(
			hash_map<sequence, array<inflected_verb>>& inflected_form_map,
			const sequence& root_id, const verb_root& root)
	{
		if (!inflected_form_map.check_size(inflected_form_map.table.size + root.present_3sg_count + root.simple_past_count + root.present_participle_count + root.past_participle_count + 3))
			return false;

		/* add the infinitive form */
		if (!add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::BARE_INFINITIVE, grammatical_tense::ANY}))
			return false;

		/* add the subjunctive form */
		if (!add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::SUBJUNCTIVE, grammatical_tense::ANY}))
			return false;

		/* add the indicative verb forms */
		for (unsigned int i = 0; i < root.present_3sg_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.present_3sg[i], {root_id, grammatical_person::THIRD, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT}))
				return false;
		}
		if (!add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_person::FIRST_OR_SECOND, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_person::ANY, grammatical_num::PLURAL, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT}))
			return false;
		for (unsigned int i = 0; i < root.simple_past_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.simple_past[i], {root_id, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::INDICATIVE, grammatical_tense::PAST}))
				return false;
		}

		/* add the participle forms */
		for (unsigned int i = 0; i < root.present_participle_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.present_participle[i], {root_id, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PRESENT_PARTICIPLE, grammatical_tense::ANY}))
				return false;
		} for (unsigned int i = 0; i < root.past_participle_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.past_participle[i], {root_id, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PAST_PARTICIPLE, grammatical_tense::ANY}))
				return false;
		}
		return true;
	}

	inline void free_helper() {
		for (auto entry : nouns) {
			core::free(entry.key);
			core::free(entry.value);
		} for (auto entry : adjectives) {
			core::free(entry.key);
			core::free(entry.value);
		} for (auto entry : adverbs) {
			core::free(entry.key);
			core::free(entry.value);
		} for (auto entry : verbs) {
			core::free(entry.key);
			core::free(entry.value);
		} for (auto entry : inflected_nouns) {
			for (auto& element : entry.value)
				core::free(element);
			core::free(entry.value);
			core::free(entry.key);
		} for (auto entry : inflected_adjectives) {
			for (auto& element : entry.value)
				core::free(element);
			core::free(entry.value);
			core::free(entry.key);
		} for (auto entry : inflected_adverbs) {
			for (auto& element : entry.value)
				core::free(element);
			core::free(entry.value);
			core::free(entry.key);
		} for (auto entry : inflected_verbs) {
			for (auto& element : entry.value)
				core::free(element);
			core::free(entry.value);
			core::free(entry.key);
		} for (auto entry : adjective_adverb_map) {
			for (auto& element : entry.value) core::free(element);
			core::free(entry.value);
			core::free(entry.key);
		}
		if (BE_SEQ.tokens != nullptr) core::free(BE_SEQ);
		if (AM_SEQ.tokens != nullptr) core::free(AM_SEQ);
		if (ARE_SEQ.tokens != nullptr) core::free(ARE_SEQ);
		if (IS_SEQ.tokens != nullptr) core::free(IS_SEQ);
		if (WAS_SEQ.tokens != nullptr) core::free(WAS_SEQ);
		if (WERE_SEQ.tokens != nullptr) core::free(WERE_SEQ);
		if (BEING_SEQ.tokens != nullptr) core::free(BEING_SEQ);
		if (BEEN_SEQ.tokens != nullptr) core::free(BEEN_SEQ);
	}
};

string morphology_en::MORE_COMPARATIVE_STRING = "_more";
string morphology_en::MOST_SUPERLATIVE_STRING = "_most";
string morphology_en::FURTHER_COMPARATIVE_STRING = "_further";
string morphology_en::FURTHEST_SUPERLATIVE_STRING = "_furthest";

bool init(morphology_en& dst, const morphology_en& src) {
	if (!hash_map_init(dst.nouns, src.nouns.table.capacity)) {
		return false;
	} else if (!hash_map_init(dst.adjectives, src.adjectives.table.capacity)) {
		free(dst.nouns);
		return false;
	} else if (!hash_map_init(dst.adverbs, src.adverbs.table.capacity)) {
		free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.verbs, src.verbs.table.capacity)) {
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.inflected_nouns, src.inflected_nouns.table.capacity)) {
		free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.inflected_adjectives, src.inflected_adjectives.table.capacity)) {
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.inflected_adverbs, src.inflected_adverbs.table.capacity)) {
		free(dst.inflected_adjectives);
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.inflected_verbs, src.inflected_verbs.table.capacity)) {
		free(dst.inflected_adverbs); free(dst.inflected_adjectives);
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.capitalization_map, src.capitalization_map.table.capacity)) {
		free(dst.inflected_verbs);
		free(dst.inflected_adverbs); free(dst.inflected_adjectives);
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.decapitalization_map, src.decapitalization_map.table.capacity)) {
		free(dst.capitalization_map); free(dst.inflected_verbs);
		free(dst.inflected_adverbs); free(dst.inflected_adjectives);
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	} else if (!hash_map_init(dst.adjective_adverb_map, src.adjective_adverb_map.table.capacity)) {
		free(dst.decapitalization_map);
		free(dst.capitalization_map); free(dst.inflected_verbs);
		free(dst.inflected_adverbs); free(dst.inflected_adjectives);
		free(dst.inflected_nouns); free(dst.verbs);
		free(dst.adverbs); free(dst.adjectives);
		free(dst.nouns); return false;
	}
	dst.BE_SEQ.tokens = nullptr;
	dst.AM_SEQ.tokens = nullptr;
	dst.ARE_SEQ.tokens = nullptr;
	dst.IS_SEQ.tokens = nullptr;
	dst.WAS_SEQ.tokens = nullptr;
	dst.WERE_SEQ.tokens = nullptr;
	dst.BEING_SEQ.tokens = nullptr;
	dst.BEEN_SEQ.tokens = nullptr;

	dst.MORE_COMPARATIVE_ID = src.MORE_COMPARATIVE_ID;
	dst.MOST_SUPERLATIVE_ID = src.MOST_SUPERLATIVE_ID;
	dst.FURTHER_COMPARATIVE_ID = src.FURTHER_COMPARATIVE_ID;
	dst.FURTHEST_SUPERLATIVE_ID = src.FURTHEST_SUPERLATIVE_ID;
	for (const auto& entry : src.nouns) {
		unsigned int index = dst.nouns.table.index_to_insert(entry.key);
		if (!init(dst.nouns.values[index], entry.value)) {
			free(dst);
			return false;
		} else if (!init(dst.nouns.table.keys[index], entry.key)) {
			free(dst.nouns.values[index]);
			free(dst);
			return false;
		}
		dst.nouns.table.size++;
	} for (const auto& entry : src.adjectives) {
		unsigned int index = dst.adjectives.table.index_to_insert(entry.key);
		if (!init(dst.adjectives.values[index], entry.value)) {
			free(dst);
			return false;
		} else if (!init(dst.adjectives.table.keys[index], entry.key)) {
			free(dst.adjectives.values[index]);
			free(dst);
			return false;
		}
		dst.adjectives.table.size++;
	} for (const auto& entry : src.adverbs) {
		unsigned int index = dst.adverbs.table.index_to_insert(entry.key);
		if (!init(dst.adverbs.values[index], entry.value)) {
			free(dst);
			return false;
		} else if (!init(dst.adverbs.table.keys[index], entry.key)) {
			free(dst.adverbs.values[index]);
			free(dst);
			return false;
		}
		dst.adverbs.table.size++;
	} for (const auto& entry : src.verbs) {
		unsigned int index = dst.verbs.table.index_to_insert(entry.key);
		if (!init(dst.verbs.values[index], entry.value)) {
			free(dst);
			return false;
		} else if (!init(dst.verbs.table.keys[index], entry.key)) {
			free(dst.verbs.values[index]);
			free(dst);
			return false;
		}
		dst.verbs.table.size++;
	} for (const auto& entry : src.inflected_nouns) {
		unsigned int index = dst.inflected_nouns.table.index_to_insert(entry.key);
		array<inflected_noun>& new_array = dst.inflected_nouns.values[index];
		if (!array_init(new_array, entry.value.length)) {
			free(dst);
			return false;
		}
		for (const inflected_noun& noun : entry.value) {
			if (!init(new_array[new_array.length], noun)) {
				for (inflected_noun& element : new_array) free(element);
				free(new_array); free(dst);
				return false;
			}
			new_array.length++;
		}
		if (!init(dst.inflected_nouns.table.keys[index], entry.key)) {
			for (inflected_noun& element : new_array) free(element);
			free(new_array); free(dst);
			return false;
		}
		dst.inflected_nouns.table.size++;
	} for (const auto& entry : src.inflected_adjectives) {
		unsigned int index = dst.inflected_adjectives.table.index_to_insert(entry.key);
		array<inflected_adjective>& new_array = dst.inflected_adjectives.values[index];
		if (!array_init(new_array, entry.value.length)) {
			free(dst);
			return false;
		}
		for (const inflected_adjective& adj : entry.value) {
			if (!init(new_array[new_array.length], adj)) {
				for (inflected_adjective& element : new_array) free(element);
				free(new_array); free(dst);
				return false;
			}
			new_array.length++;
		}
		if (!init(dst.inflected_adjectives.table.keys[index], entry.key)) {
			for (inflected_adjective& element : new_array) free(element);
			free(new_array); free(dst);
			return false;
		}
		dst.inflected_adjectives.table.size++;
	} for (const auto& entry : src.inflected_adverbs) {
		unsigned int index = dst.inflected_adverbs.table.index_to_insert(entry.key);
		array<inflected_adverb>& new_array = dst.inflected_adverbs.values[index];
		if (!array_init(new_array, entry.value.length)) {
			free(dst);
			return false;
		}
		for (const inflected_adverb& adv : entry.value) {
			if (!init(new_array[new_array.length], adv)) {
				for (inflected_adverb& element : new_array) free(element);
				free(new_array); free(dst);
				return false;
			}
			new_array.length++;
		}
		if (!init(dst.inflected_adverbs.table.keys[index], entry.key)) {
			for (inflected_adverb& element : new_array) free(element);
			free(new_array); free(dst);
			return false;
		}
		dst.inflected_adverbs.table.size++;
	} for (const auto& entry : src.inflected_verbs) {
		unsigned int index = dst.inflected_verbs.table.index_to_insert(entry.key);
		array<inflected_verb>& new_array = dst.inflected_verbs.values[index];
		if (!array_init(new_array, entry.value.length)) {
			free(dst);
			return false;
		}
		for (const inflected_verb& verb : entry.value) {
			if (!init(new_array[new_array.length], verb)) {
				for (inflected_verb& element : new_array) free(element);
				free(new_array); free(dst);
				return false;
			}
			new_array.length++;
		}
		if (!init(dst.inflected_verbs.table.keys[index], entry.key)) {
			for (inflected_verb& element : new_array) free(element);
			free(new_array); free(dst);
			return false;
		}
		dst.inflected_verbs.table.size++;
	} for (const auto& entry : src.adjective_adverb_map) {
		unsigned int index = dst.adjective_adverb_map.table.index_to_insert(entry.key);
		array<sequence>& new_array = dst.adjective_adverb_map.values[index];
		if (!array_init(new_array, entry.value.length)) {
			free(dst);
			return false;
		}
		for (const sequence& seq : entry.value) {
			if (!init(new_array[new_array.length], seq)) {
				for (sequence& element : new_array) free(element);
				free(new_array); free(dst);
				return false;
			}
			new_array.length++;
		}
		if (!init(dst.adjective_adverb_map.table.keys[index], entry.key)) {
			for (sequence& element : new_array) free(element);
			free(new_array); free(dst);
			return false;
		}
		dst.adjective_adverb_map.table.size++;
	}

	dst.capitalization_map.put_all(src.capitalization_map);
	dst.decapitalization_map.put_all(src.decapitalization_map);

	if ((src.BE_SEQ.tokens != nullptr && !init(dst.BE_SEQ, src.BE_SEQ))
	 || (src.AM_SEQ.tokens != nullptr && !init(dst.AM_SEQ, src.AM_SEQ))
	 || (src.ARE_SEQ.tokens != nullptr && !init(dst.ARE_SEQ, src.ARE_SEQ))
	 || (src.IS_SEQ.tokens != nullptr && !init(dst.IS_SEQ, src.IS_SEQ))
	 || (src.WAS_SEQ.tokens != nullptr && !init(dst.WAS_SEQ, src.WAS_SEQ))
	 || (src.WERE_SEQ.tokens != nullptr && !init(dst.WERE_SEQ, src.WERE_SEQ))
	 || (src.BEING_SEQ.tokens != nullptr && !init(dst.BEING_SEQ, src.BEING_SEQ))
	 || (src.BEEN_SEQ.tokens != nullptr && !init(dst.BEEN_SEQ, src.BEEN_SEQ)))
	{
		free(dst);
		return false;
	}
	return true;
}

enum class morphology_state {
	DEFAULT,
	ENTRY
};

inline bool starts_with(const string& src, const char* pattern) {
	for (unsigned int i = 0; i < src.length; i++) {
		if (pattern[i] == '\0') {
			return true;
		} else if (src[i] != pattern[i]) {
			return false;
		}
	}
	return (pattern[src.length] == '\0');
}

inline bool string_compare(
		const char* str, unsigned int str_length,
		const char* pattern, unsigned int pattern_length)
{
	if (str_length < pattern_length) return false;
	for (unsigned int i = 0; i < pattern_length; i++)
		if (str[i] != pattern[i]) return false;
	return true;
}

inline bool get_parameter(
		const string& src, const string& left_qual_name, const string& right_qual_name,
		unsigned int& index, const char*& qualifier, unsigned int& qualifier_length)
{
#if !defined(NDEBUG)
	if (right_qual_name.length == 0)
		fprintf(stderr, "get_parameter WARNING: `right_qual_name` is empty.\n");
#endif

	if (src.length < left_qual_name.length + right_qual_name.length)
		return false;
	if (!string_compare(src.data, src.length, left_qual_name.data, left_qual_name.length))
		return false;
	unsigned int right_index = index_of(right_qual_name[0], src.data + left_qual_name.length, src.length - left_qual_name.length) + left_qual_name.length;

	if (src.length - right_index < right_qual_name.length) return false;
	if (right_index == left_qual_name.length) {
		index = 1;
	} else if (!parse_uint(string(src.data + left_qual_name.length, right_index - left_qual_name.length), index)) {
		return false;
	}
	index--;

	if (!string_compare(src.data + right_index, src.length - right_index, right_qual_name.data, right_qual_name.length))
		return false;
	qualifier = src.data + right_index + right_qual_name.length;
	qualifier_length = src.length - right_index - right_qual_name.length;
	return true;
}

inline bool is_vowel(char c) {
	return (c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u');
}

inline bool inflect_verb(array<string>& dst, const string& first) {
	if (dst.length == 0) {
		if (!init(dst[0], first.data, first.length)) return false;
		dst.length = 1;
	} else if (dst[0].data != nullptr && dst[0].length == 0) {
		string inflected(first.data, first.length);
		swap(dst[0], inflected);
	}
	return true;
}

inline bool inflect_verb(array<string>& dst, const string& first, const char* second) {
	if (dst.length == 0) {
		if (!init(dst[0], first.data, first.length)) return false;
		dst[0] += second;
		dst.length = 1;
	} else if (dst[0].data != nullptr && dst[0].length == 0) {
		string inflected(first.data, first.length);
		inflected += second;
		swap(dst[0], inflected);
	}
	return true;
}

inline bool inflect_verb(array<string>& dst, const string& first, const string& second, const char* third) {
	if (dst.length == 0) {
		if (!init(dst[0], first.data, first.length)) return false;
		dst[0] += second;
		dst[0] += third;
		dst.length = 1;
	} else if (dst[0].data != nullptr && dst[0].length == 0) {
		string inflected(first.data, first.length);
		inflected += second;
		inflected += third;
		swap(dst[0], inflected);
	}
	return true;
}

inline bool inflect_verb(array<string>& dst, unsigned int index, const char* first, unsigned int first_length)
{
	if (index < dst.length) {
		if (dst[index].data == nullptr)
			return true;
		string inflected(first, first_length);
		swap(inflected, dst[index]);
	} else {
		if (!dst.ensure_capacity(index + 1)) return false;
		while (dst.length < index) {
			if (!init(dst[dst.length], "")) return false;
			dst.length++;
		}
		if (!init(dst[index], first, first_length)) return false;
		dst.length++;
	}
	return true;
}

struct wikt_entry {
	string root;
	array<string> entries;
	position pos;

	wikt_entry(position initial_pos) : root(""), entries(16), pos(initial_pos) { }

	~wikt_entry() { free_helper(); }

	static inline void free(wikt_entry& entry) {
		entry.free_helper();
		core::free(entry.entries);
		core::free(entry.root);
	}

private:
	inline void free_helper() {
		for (string& str : entries) core::free(str);
	}
};

inline bool init(wikt_entry& entry, const wikt_entry& src) {
	if (!init(entry.root, src.root)) {
		return false;
	} else if (!array_init(entry.entries, src.entries.length)) {
		free(entry.root);
		return false;
	}
	for (unsigned int i = 0; i < src.entries.length; i++) {
		if (!init(entry.entries[i], src.entries[i])) {
			core::free(entry);
			return false;
		}
		entry.entries.length++;
	}
	entry.pos = src.pos;
	return true;
}

template<bool IsAdverb = false>
inline bool emit_adj_entry(const sequence& root_ids,
		morphology_en& m, const wikt_entry& entry,
		const sequence& adj_root,
		hash_map<string, unsigned int>& names)
{
	if (entry.entries.length == 1) {
		/* the comparative is formed with "more" and the superlative with "most" */
		pair<string, string>& inflected_forms = *((pair<string, string>*) alloca(sizeof(pair<string, string>)));
		inflected_forms.key = morphology_en::MORE_COMPARATIVE_STRING;
		inflected_forms.value = morphology_en::MOST_SUPERLATIVE_STRING;

		adjective_root& new_root = *((adjective_root*) alloca(sizeof(adjective_root)));
		if (!init(new_root, comparability::COMPARABLE, &inflected_forms, 1, adj_root, names)) {
			free(inflected_forms);
			return false;
		}
		free(inflected_forms);
		if (IsAdverb)
			return m.add_adverb_root(root_ids, new_root);
		else return m.add_adjective_root(root_ids, new_root);
	} else {
		array<pair<string, string>> inflected_forms(entry.entries.length - 1);
		array_map<unsigned int, string> superlative_forms(entry.entries.length - 1);
		comparability comp = comparability::COMPARABLE;
		for (unsigned int i = 1; i < entry.entries.length; i++) {
			if (entry.entries[i] == "more" && entry.root != "many" && entry.root != "most") {
				if (!init(inflected_forms[inflected_forms.length].key, morphology_en::MORE_COMPARATIVE_STRING)) {
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				} else if (!init(inflected_forms[inflected_forms.length].value, morphology_en::MOST_SUPERLATIVE_STRING)) {
					free(inflected_forms[inflected_forms.length].key);
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				}
				inflected_forms.length++;
			} else if (entry.entries[i] == "further" && entry.root != "far") {
				if (!init(inflected_forms[inflected_forms.length].key, morphology_en::FURTHER_COMPARATIVE_STRING)) {
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				} else if (!init(inflected_forms[inflected_forms.length].value, morphology_en::FURTHEST_SUPERLATIVE_STRING)) {
					free(inflected_forms[inflected_forms.length].key);
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				}
				inflected_forms.length++;
			} else if (entry.entries[i] == "er") {
				string stem(entry.root.data, entry.root.length);
				if (entry.root.length > 0 && entry.root[entry.root.length - 1] == 'e') {
					stem.length--;
				} else if (entry.root.length > 2 && entry.root[entry.root.length - 2] == 'e' && entry.root[entry.root.length - 1] == 'y'
						&& (entry.root.length < 3 || !is_vowel(entry.root[entry.root.length - 3])))
				{
					stem.length--;
					stem[stem.length - 1] = 'i';
				} else if (entry.root.length > 1 && entry.root[entry.root.length - 1] == 'y'
						&& (entry.root.length < 2 || !is_vowel(entry.root[entry.root.length - 2])))
				{
					stem[stem.length - 1] = 'i';
				}

				if (!init(inflected_forms[inflected_forms.length].key, stem)) {
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				} else if (!init(inflected_forms[inflected_forms.length].value, stem)) {
					free(inflected_forms[inflected_forms.length].key);
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				}

				inflected_forms[inflected_forms.length].key += "er";
				inflected_forms[inflected_forms.length++].value += "est";
			} else if (entry.entries[i] == "-") {
				comp = comparability::INCOMPARABLE;
			} else if (entry.entries[i] == "?") {
				comp = comparability::UNKNOWN;
			} else if (entry.entries[i] == "+") {
				comp = comparability::COMPARABLE_ONLY;
			} else {
				unsigned int comparative_index;
				const char* superlative; unsigned int superlative_length;
				if (get_parameter(entry.entries[i], "sup", "=", comparative_index, superlative, superlative_length)) {
					unsigned int index = superlative_forms.index_of(comparative_index);
					if (index < superlative_forms.size)
						free(superlative_forms.values[index]);
					superlative_forms.keys[index] = comparative_index;
					if (!init(superlative_forms.values[index], superlative, superlative_length)) {
						if (index < superlative_forms.size) {
							move(superlative_forms.values[superlative_forms.size - 1], superlative_forms.values[index]);
							superlative_forms.size--;
						}
						for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
						for (auto entry : superlative_forms) free(entry.value);
						return false;
					}
					if (index == superlative_forms.size)
						superlative_forms.size++;
				} else {
					if (entry.entries[i].index_of('=') < entry.entries[i].length || entry.entries[i].index_of('[') < entry.entries[i].length) {
						for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
						for (auto entry : superlative_forms) free(entry.value);
						return true;
					}

					if (!init(inflected_forms[inflected_forms.length].key, entry.entries[i])) {
						for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
						for (auto entry : superlative_forms) free(entry.value);
						return false;
					}

					if (entry.entries[i].length > 2 && entry.entries[i][entry.entries[i].length - 2] == 'e' && entry.entries[i][entry.entries[i].length - 1] == 'r') {
						string superlative(entry.entries[i].data, entry.entries[i].length);
						superlative.length -= 2;
						superlative += "est";
						if (!init(inflected_forms[inflected_forms.length].value, superlative)) {
							free(inflected_forms[inflected_forms.length].key);
							for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
							for (auto entry : superlative_forms) free(entry.value);
							return false;
						}
					} else if (!init(inflected_forms[inflected_forms.length].value, "")) {
						free(inflected_forms[inflected_forms.length].key);
						for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
						for (auto entry : superlative_forms) free(entry.value);
						return false;
					}
					inflected_forms.length++;
				}
			}
		}

		for (auto entry : superlative_forms) {
			/* if `inflected_forms.length` is smaller than the index, expand
				`inflected_forms` with empty entries until its sufficiently large */
			if (!inflected_forms.ensure_capacity(entry.key + 1)) {
				for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
				for (auto entry : superlative_forms) free(entry.value);
				return false;
			}
			while (entry.key >= inflected_forms.length) {
				if (!init(inflected_forms[inflected_forms.length].key, "")) {
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				} else if (!init(inflected_forms[inflected_forms.length].value, "")) {
					free(inflected_forms[inflected_forms.length].key);
					for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
					for (auto entry : superlative_forms) free(entry.value);
					return false;
				}
				inflected_forms.length++;
			}
			swap(inflected_forms[entry.key].value, entry.value);
		}
		for (auto entry : superlative_forms) free(entry.value);

		adjective_root& new_root = *((adjective_root*) alloca(sizeof(adjective_root)));
		if (!init(new_root, comp, inflected_forms.data, inflected_forms.length, adj_root, names))
			return false;
		for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
		if (IsAdverb)
			return m.add_adverb_root(root_ids, new_root);
		else return m.add_adjective_root(root_ids, new_root);
	}
}

inline bool emit_adv_entry(
		morphology_en& m, const wikt_entry& entry,
		hash_map<string, unsigned int>& names)
{
	array<unsigned int> tokenization(2);
	if (!tokenize(entry.root.data, entry.root.length, tokenization, names))
		return false;
	sequence root_ids(tokenization.data, (unsigned int) tokenization.length);

	/* check if we can form this adverb from an adjective */
	bool contains;
	m.adjectives.get(root_ids, contains);
	if (contains) {
		/* the adverb and adjective forms are the same */
		return emit_adj_entry<true>(root_ids, m, entry, root_ids, names);
	}

	if (entry.root.length > 2 && entry.root[entry.root.length - 2] == 'l' && entry.root[entry.root.length - 1] == 'y') {
		string substring(entry.root.data, entry.root.length - 2);
		array<unsigned int> new_tokenization(2);
		if (!tokenize(substring.data, substring.length, new_tokenization, names))
			return false;
		sequence new_root_ids(new_tokenization.data, (unsigned int) new_tokenization.length);
		m.adjectives.get(new_root_ids, contains);
		if (contains)
			/* we have an adverb of the form '-ly' (e.g. "quickly", "slowly") */
			return emit_adj_entry<true>(root_ids, m, entry, new_root_ids, names);
	}

	if (entry.root.length > 2 && entry.root[entry.root.length - 2] == 'l' && entry.root[entry.root.length - 1] == 'y') {
		string substring(entry.root.data, entry.root.length - 1);
		substring += "e";

		array<unsigned int> new_tokenization(2);
		if (!tokenize(substring.data, substring.length, new_tokenization, names))
			return false;
		sequence new_root_ids(new_tokenization.data, (unsigned int) new_tokenization.length);
		m.adjectives.get(new_root_ids, contains);
		if (contains)
			/* we have an adverb of the form '-ly' (e.g. "possibly", "forcibly") */
			return emit_adj_entry<true>(root_ids, m, entry, new_root_ids, names);
	}

	if (entry.root.length > 6
	 && entry.root[entry.root.length - 6] == 'i' && entry.root[entry.root.length - 5] == 'c'
	 && entry.root[entry.root.length - 4] == 'a' && entry.root[entry.root.length - 3] == 'l'
	 && entry.root[entry.root.length - 2] == 'l' && entry.root[entry.root.length - 1] == 'y')
	{
		string substring(entry.root.data, entry.root.length - 4);
		array<unsigned int> new_tokenization(2);
		if (!tokenize(substring.data, substring.length, new_tokenization, names))
			return false;
		sequence new_root_ids(new_tokenization.data, (unsigned int) new_tokenization.length);
		m.adjectives.get(new_root_ids, contains);
		if (contains)
			/* we have an adverb of the form '-ically' */
			return emit_adj_entry<true>(root_ids, m, entry, new_root_ids, names);
	}

	if (entry.root.length > 3
	 && entry.root[entry.root.length - 3] == 'i' && entry.root[entry.root.length - 2] == 'l'
	 && entry.root[entry.root.length - 1] == 'y')
	{
		string substring(entry.root.data, entry.root.length - 3);
		substring += "y";

		array<unsigned int> new_tokenization(2);
		if (!tokenize(substring.data, substring.length, new_tokenization, names))
			return false;
		sequence new_root_ids(new_tokenization.data, (unsigned int) new_tokenization.length);
		m.adjectives.get(new_root_ids, contains);
		if (contains)
			/* we have an adverb of the form '-ily' */
			return emit_adj_entry<true>(root_ids, m, entry, new_root_ids, names);
	}

	/* the adverb is not formed from an adjective */
	return emit_adj_entry<true>(root_ids, m, entry, sequence(nullptr, 0), names);
}

inline bool emit_entry(
		morphology_en& m, const wikt_entry& entry,
		hash_map<string, unsigned int>& names,
		array<wikt_entry>& adv_entries)
{
	array<unsigned int> tokenization(2);
	if (!tokenize(entry.root.data, entry.root.length, tokenization, names))
		return false;
	sequence root_ids(tokenization.data, (unsigned int) tokenization.length);

	if (entry.entries.length == 0) {
		read_error("Found an entry with no elements", entry.pos);
		return false;
	} else if (entry.entries[0] == "n") {
		if (entry.entries.length == 1) {
			/* the noun is countable and has a simple plural form "-s" */
			string plural(entry.root.data, entry.root.length);
			plural += "s";

			noun_root new_root;
			if (!init(new_root, false, countability::COUNTABLE, &plural, 1, names))
				return false;
			return m.add_noun_root(root_ids, new_root);
		} else {
			static constexpr const char* rare_str = "rare";
			static unsigned int rare_str_length = (unsigned int) strlen(rare_str);
			static constexpr const char* by_suppletion_str = "by [[suppletion]]";
			static unsigned int by_suppletion_str_length = (unsigned int) strlen(by_suppletion_str);
			static constexpr const char* nonstandard_str = "nonstandard";
			static unsigned int nonstandard_str_length = (unsigned int) strlen(nonstandard_str);
			static constexpr const char* archaic_str = "archaic";
			static unsigned int archaic_str_length = (unsigned int) strlen(archaic_str);

			array<string> plural_forms(entry.entries.length - 1);
			countability count = countability::COUNTABLE;
			for (unsigned int i = 1; i < entry.entries.length; i++) {
				if (entry.entries[i] == "s") {
					if (!init(plural_forms[plural_forms.length], entry.root)) {
						for (string& str : plural_forms) free(str);
						return false;
					}
					plural_forms[plural_forms.length++] += "s";
				} else if (entry.entries[i] == "es") {
					if (!init(plural_forms[plural_forms.length], entry.root)) {
						for (string& str : plural_forms) free(str);
						return false;
					}
					plural_forms[plural_forms.length++] += "es";
				} else if (entry.entries[i] == "-") {
					count = countability::UNCOUNTABLE;
				} else if (entry.entries[i] == "~") {
					count = countability::BOTH;
				} else if (entry.entries[i] == "!") {
					count = countability::PLURAL_UNATTESTED;
				} else if (entry.entries[i] == "?") {
					count = countability::PLURAL_UNKNOWN;
				} else if (starts_with(entry.entries[i], "head=")) {
					/* ignore these for now */
					for (string& str : plural_forms) free(str);
					return true;
				} else {
					unsigned int qualifier_index; const char* qualifier; unsigned int qualifier_length;
					if (get_parameter(entry.entries[i], "pl", "qual=", qualifier_index, qualifier, qualifier_length)) {
						if (!string_compare(qualifier, qualifier_length, rare_str, rare_str_length)
						 && !string_compare(qualifier, qualifier_length, by_suppletion_str, by_suppletion_str_length)
						 && !string_compare(qualifier, qualifier_length, nonstandard_str, nonstandard_str_length)
						 && !string_compare(qualifier, qualifier_length, archaic_str, archaic_str_length))
						{
							for (string& str : plural_forms) free(str);
							return true;
						}
						continue;
					} else if (entry.entries[i].index_of('=') < entry.entries[i].length || entry.entries[i].index_of('[') < entry.entries[i].length) {
						for (string& str : plural_forms) free(str);
						return true;
					}
					plural_forms[plural_forms.length++] = entry.entries[i];
				}
			}

			noun_root new_root;
			if (!init(new_root, false, count, plural_forms.data, (unsigned int) plural_forms.length, names))
				return false;
			for (string& str : plural_forms) free(str);
			return m.add_noun_root(root_ids, new_root);
		}
	} else if (entry.entries[0] == "pr") {
		/* this is a proper noun */
		array<string> plural_forms(max((size_t) 1, entry.entries.length - 1));
		countability count = countability::UNCOUNTABLE;
		for (unsigned int i = 1; i < entry.entries.length; i++) {
			if (entry.entries[i] == "s") {
				if (!init(plural_forms[plural_forms.length], entry.root)) {
					for (string& str : plural_forms) free(str);
					return false;
				}
				plural_forms[plural_forms.length++] += "s";
			} else if (entry.entries[i] == "es") {
				if (!init(plural_forms[plural_forms.length], entry.root)) {
					for (string& str : plural_forms) free(str);
					return false;
				}
				plural_forms[plural_forms.length++] += "es";
			} else if (entry.entries[i] == "-") {
				count = countability::UNCOUNTABLE;
			} else if (entry.entries[i] == "~") {
				count = countability::BOTH;
			} else if (starts_with(entry.entries[i], "head=")) {
				/* ignore these for now */
				continue;
			} else {
				if (entry.entries[i].index_of('=') < entry.entries[i].length || entry.entries[i].index_of('[') < entry.entries[i].length)
					continue;
				plural_forms[plural_forms.length++] = entry.entries[i];
			}
		}

		noun_root new_root;
		if (!init(new_root, true, count, plural_forms.data, (unsigned int) plural_forms.length, names))
			return false;
		for (string& str : plural_forms) free(str);
		return m.add_noun_root(root_ids, new_root);
	} else if (entry.entries[0] == "pl") {
		/* the noun is plural only; we ignore these for now */
		return true;
	} else if (entry.entries[0] == "adv") {
		if (!adv_entries.ensure_capacity(adv_entries.length + 1)
		 || !init(adv_entries[adv_entries.length], entry))
		{
			return false;
		}
		adv_entries.length++;
		return true;
	} else if (entry.entries[0] == "adj") {
		return emit_adj_entry(root_ids, m, entry, sequence(nullptr, 0), names);
	} else if (entry.entries[0] == "v") {
		if (entry.entries.length == 1) {
			/* the noun is regular, i.e. the present 3rd person singular is
			   formed by adding "-s", the present participle by adding "-ing",
			   and the past by adding "-ed" */
			string present_3sg(entry.root.data, entry.root.length);
			present_3sg += "s";
			string present_participle(entry.root.data, entry.root.length);
			present_participle += "ing";
			string past(entry.root.data, entry.root.length);
			past += "ed";

			verb_root new_root;
			if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
				return false;
			return m.add_verb_root(root_ids, new_root);
		}

		static constexpr const char* obsolete_str = "obsolete";
		static unsigned int obsolete_str_length = (unsigned int) strlen(obsolete_str);

		/* first process the flags */
		array<string> args(entry.entries.length - 1);
		array<string> pres_3sg(5);
		array<string> pres_ptc(5);
		array<string> past(5);
		array<string> past_ptc(5);
		auto free_strings = [&]() {
			for (string& str : pres_3sg) { if (str.data != nullptr) free(str); }
			for (string& str : pres_ptc) { if (str.data != nullptr) free(str); }
			for (string& str : past) { if (str.data != nullptr) free(str); }
			for (string& str : past_ptc) { if (str.data != nullptr) free(str); }
			for (string& str : args) { free(str); }
		};
		for (unsigned int i = 1; i < entry.entries.length; i++) {
			unsigned int qualifier_index; const char* qualifier; unsigned int qualifier_length;
			if (entry.entries[i].index_of('[') < entry.entries[i].length) { free_strings(); return true; }
			if (get_parameter(entry.entries[i], "head", "=", qualifier_index, qualifier, qualifier_length)) { free_strings(); return true; }

			if (get_parameter(entry.entries[i], "pres_3sg", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(pres_3sg, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry.entries[i], "pres_ptc", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(pres_ptc, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry.entries[i], "past_ptc", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(past_ptc, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry.entries[i], "past", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(past, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry.entries[i], "pres_3sg", "_qual=", qualifier_index, qualifier, qualifier_length)) {
				if (string_compare(qualifier, qualifier_length, obsolete_str, obsolete_str_length)) {
					if (qualifier_index < pres_3sg.length) {
						if (pres_3sg[qualifier_index].data != nullptr) {
							free(pres_3sg[qualifier_index]); pres_3sg[qualifier_index].data = nullptr;
						}
					} else {
						if (!pres_3sg.ensure_capacity(qualifier_index + 1)) { free_strings(); return false; }
						while (pres_3sg.length < qualifier_index) {
							if (!init(pres_3sg[pres_3sg.length], "")) { free_strings(); return false; }
							pres_3sg.length++;
						}
						pres_3sg[qualifier_index].data = nullptr;
						pres_3sg.length++;
					}
				}
				continue;
			} else if (get_parameter(entry.entries[i], "pres_ptc", "_qual=", qualifier_index, qualifier, qualifier_length)) {
				if (string_compare(qualifier, qualifier_length, obsolete_str, obsolete_str_length)) {
					if (qualifier_index < pres_ptc.length) {
						if (pres_ptc[qualifier_index].data != nullptr) {
							free(pres_ptc[qualifier_index]); pres_ptc[qualifier_index].data = nullptr;
						}
					} else {
						if (!pres_ptc.ensure_capacity(qualifier_index + 1)) { free_strings(); return false; }
						while (pres_ptc.length < qualifier_index) {
							if (!init(pres_ptc[pres_ptc.length], "")) { free_strings(); return false; }
							pres_ptc.length++;
						}
						pres_ptc[qualifier_index].data = nullptr;
						pres_ptc.length++;
					}
				}
				continue;
			} else if (get_parameter(entry.entries[i], "past_ptc", "_qual=", qualifier_index, qualifier, qualifier_length)) {
				if (string_compare(qualifier, qualifier_length, obsolete_str, obsolete_str_length)) {
					if (qualifier_index < past_ptc.length) {
						if (past_ptc[qualifier_index].data != nullptr) {
							free(past_ptc[qualifier_index]); past_ptc[qualifier_index].data = nullptr;
						}
					} else {
						if (!past_ptc.ensure_capacity(qualifier_index + 1)) { free_strings(); return false; }
						while (past_ptc.length < qualifier_index) {
							if (!init(past_ptc[past_ptc.length], "")) { free_strings(); return false; }
							past_ptc.length++;
						}
						past_ptc[qualifier_index].data = nullptr;
						past_ptc.length++;
					}
				}
				continue;
			} else if (get_parameter(entry.entries[i], "past", "_qual=", qualifier_index, qualifier, qualifier_length)) {
				if (string_compare(qualifier, qualifier_length, obsolete_str, obsolete_str_length)) {
					if (qualifier_index < past.length) {
						if (past[qualifier_index].data != nullptr) {
							free(past[qualifier_index]); past[qualifier_index].data = nullptr;
						}
					} else {
						if (!past.ensure_capacity(qualifier_index + 1)) { free_strings(); return false; }
						while (past.length < qualifier_index) {
							if (!init(past[past.length], "")) { free_strings(); return false; }
							past.length++;
						}
						past[qualifier_index].data = nullptr;
						past.length++;
					}
				}
				continue;
			}

#if !defined(NDEBUG)
			unsigned int eq_index = entry.entries[i].index_of('=');
			if (eq_index < entry.entries[i].length) {
				fprintf(stderr, "WARNING at %u:%u: Unrecognized parameter in entry '", entry.pos.line, entry.pos.column);
				print(entry.entries[i], stderr); print("'.\n", stderr);
				continue;
			}
#endif

			if (!init(args[args.length], entry.entries[i])) {
				free_strings(); return false;
			}
			args.length++;
		}

		bool first_empty = (args.length < 1 || args[0].length == 0);
		bool second_empty = (args.length < 2 || args[1].length == 0);
		bool third_empty = (args.length < 3 || args[2].length == 0);
		if (!first_empty && second_empty && third_empty) {
			if (args[0] == "es") {
				if (!inflect_verb(pres_3sg, entry.root, "es")
				 || !inflect_verb(pres_ptc, entry.root, "ing")
				 || !inflect_verb(past, entry.root, "ed"))
				{ free_strings(); return false; }
			} else if (args[0] == "ies") {
				if (entry.root.length == 0 || entry.root[entry.root.length - 1] != 'y') {
					fprintf(stderr, "ERROR at %u:%u: Verb root '", entry.pos.line, entry.pos.column);
					print(entry.root, stderr); fprintf(stderr, "' does not end in 'y'.\n");
					return false;
				}

				string stem(entry.root.data, entry.root.length - 1);
				if (!inflect_verb(pres_3sg, stem, "ies")
				 || !inflect_verb(pres_ptc, stem, "ying")
				 || !inflect_verb(past, stem, "ied"))
				{ free_strings(); return false; }
			} else if (args[0] == "d") {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, entry.root, "ing")
				 || !inflect_verb(past, entry.root, "d"))
				{ free_strings(); return false; }
			} else {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "ed"))
				{ free_strings(); return false; }
			}
		} else if (third_empty) {
			if (args.length > 1 && args[1] == "es") {
				if (!inflect_verb(pres_3sg, args[0], "es")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "ed"))
				{ free_strings(); return false; }
			} else if (args.length > 1 && args[1] == "ies") {
				if (!inflect_verb(pres_3sg, args[0], "ies")
				 || !inflect_verb(pres_ptc, args[0], "ying")
				 || !inflect_verb(past, args[0], "ied"))
				{ free_strings(); return false; }
			} else if (args.length > 1 && (args[1] == "ing" || args[1] == "ed")) {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "ed"))
				{ free_strings(); return false; }
			} else if (args.length > 1 && args[1] == "d") {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "d"))
				{ free_strings(); return false; }
			} else {
				if (first_empty) {
					if (!inflect_verb(pres_3sg, entry.root, "s")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, args[0])) { free_strings(); return false; }
				}
				if (second_empty) {
					if (!inflect_verb(pres_ptc, entry.root, "ing")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_ptc, args[1])) { free_strings(); return false; }
				}
				if (!inflect_verb(past, entry.root, "ed")) { free_strings(); return false; }
			}
		} else {
			if (args[2] == "es") {
				if (!inflect_verb(pres_3sg, args[0], args[1], "es")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing")
				 || !inflect_verb(past, args[0], args[1], "ed"))
				{ free_strings(); return false; }
			} else if (args[2] == "ing") {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing"))
				{ free_strings(); return false; }

				if (args[1] == "y") {
					if (!inflect_verb(past, entry.root, "d")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(past, args[0], args[1], "ed")) { free_strings(); return false; }
				}
			} else if (args[2] == "ed") {
				if (args[1] == "i") {
					if (!inflect_verb(pres_3sg, args[0], args[1], "es")
					 || !inflect_verb(pres_ptc, entry.root, "ing"))
					{ free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, entry.root, "s")
					 || !inflect_verb(pres_ptc, args[0], args[1], "ing"))
					{ free_strings(); return false; }
				}
				if (!inflect_verb(past, args[0], args[1], "ed")) { free_strings(); return false; }
			} else if (args[2] == "d") {
				if (!inflect_verb(pres_3sg, entry.root, "s")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing")
				 || !inflect_verb(past, args[0], args[1], "d"))
				{ free_strings(); return false; }
			} else {
				if (first_empty) {
					if (!inflect_verb(pres_3sg, entry.root, "s")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, args[0])) { free_strings(); return false; }
				}
				if (second_empty) {
					if (!inflect_verb(pres_ptc, entry.root, "ing")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_ptc, args[1])) { free_strings(); return false; }
				}
				if (!inflect_verb(past, args[2])) { free_strings(); return false; }
			}
		}

		if (past_ptc.length == 0 || (past_ptc[0].data != nullptr && past_ptc[0].length == 0)) {
			if (args.length > 3 && args[3].length != 0) {
				if (past_ptc.length != 0) {
					string inflected(args[3].data, args[3].length);
					swap(past_ptc[0], inflected);
				} else {
					if (!init(past_ptc[0], args[3].data, args[3].length)) { free_strings(); return false; }
					if (past_ptc.length == 0) past_ptc.length++;
				}
			} else if (past_ptc.length == 0) {
				for (const string& str : past) {
					if (str.data == nullptr) continue;
					if (!init(past_ptc[past_ptc.length], str)) { free_strings(); return false; }
					past_ptc.length++;
				}
			}
		}

		/* remove empty inflected forms */
		for (unsigned int i = 0; i < pres_3sg.length; i++) {
			if (pres_3sg[i].data == nullptr) { pres_3sg.remove(i--); continue; }
			if (pres_3sg[i] == "-") { free(pres_3sg[i]); pres_3sg.remove(i--); continue; }
		} for (unsigned int i = 0; i < pres_ptc.length; i++) {
			if (pres_ptc[i].data == nullptr) { pres_ptc.remove(i--); continue; }
			if (pres_ptc[i] == "-") { free(pres_ptc[i]); pres_ptc.remove(i--); continue; }
			if (pres_ptc[i] == "-ing") { free_strings(); return true; }
		} for (unsigned int i = 0; i < past.length; i++) {
			if (past[i].data == nullptr) { past.remove(i--); continue; }
			if (past[i] == "-") { free(past[i]); past.remove(i--); continue; }
		} for (unsigned int i = 0; i < past_ptc.length; i++) {
			if (past_ptc[i].data == nullptr) { past_ptc.remove(i--); continue; }
			if (past_ptc[i] == "-") { free(past_ptc[i]); past_ptc.remove(i--); continue; }
		}

		if (pres_3sg.length == 0 || pres_ptc.length == 0 || past.length == 0 || past_ptc.length == 0) {
			free_strings();
			return true;
		}

		verb_root new_root;
		if (!init(new_root,
				pres_3sg.data, (unsigned int) pres_3sg.length,
				pres_ptc.data, (unsigned int) pres_ptc.length,
				past.data, (unsigned int) past.length,
				past_ptc.data, (unsigned int) past_ptc.length, names))
		{
			free_strings();
			return false;
		}
		free_strings();
		return m.add_verb_root(root_ids, new_root);
	} else {
		read_error("Unrecognized part of speech", entry.pos);
		return false;
	}
	return true;
}

template<typename Stream>
bool morphology_read(morphology_en& m,
		hash_map<string, unsigned int>& names,
		Stream& input)
{
	if (!get_token(morphology_en::MORE_COMPARATIVE_STRING, m.MORE_COMPARATIVE_ID, names)
	 || !get_token(morphology_en::MOST_SUPERLATIVE_STRING, m.MOST_SUPERLATIVE_ID, names)
	 || !get_token(morphology_en::FURTHER_COMPARATIVE_STRING, m.FURTHER_COMPARATIVE_ID, names)
	 || !get_token(morphology_en::FURTHEST_SUPERLATIVE_STRING, m.FURTHEST_SUPERLATIVE_ID, names))
		return false;

	position start(1, 1);
	morphology_state state = morphology_state::DEFAULT;
	array<char> token = array<char>(1024);
	wikt_entry current_entry(start);
	array<wikt_entry> adv_entries(1024);

	mbstate_t shift = {0};
	buffered_stream<MB_LEN_MAX, Stream> wrapper(input);
	char32_t next = fgetc32(wrapper);
	bool new_line = false;
	while (next != static_cast<char32_t>(-1)) {
		switch (state) {
		case morphology_state::DEFAULT:
			if (next == '\t') {
				string new_root(token.data, token.length);
				swap(current_entry.root, new_root);
				state = morphology_state::ENTRY;
				token.clear(); shift = {0};
			} else if (next == '\n') {
				fprintf(stderr, "WARNING: Found root with no entries on line %u.\n", current_entry.pos.line);
				state = morphology_state::DEFAULT;
				new_line = true;
				token.clear(); shift = {0};
			} else {
				if (!append_to_token(token, next, shift)) {
					for (auto& entry : adv_entries) free(entry);
					return false;
				}
			}
			break;

		case morphology_state::ENTRY:
			if (next == '\t' || next == '\n') {
				if (!current_entry.entries.ensure_capacity(current_entry.entries.length + 1)
				 || !init(current_entry.entries[current_entry.entries.length++], token.data, token.length)
				 || !emit_entry(m, current_entry, names, adv_entries))
				{
					for (auto& entry : adv_entries) free(entry);
					return false;
				}
				for (string& str : current_entry.entries) { free(str); }
				current_entry.entries.clear();
				token.clear(); shift = {0};
				if (next == '\n') {
					state = morphology_state::DEFAULT;
					new_line = true;
				}
			} else if (next == '\\') {
				if (!current_entry.entries.ensure_capacity(current_entry.entries.length + 1)
				 || !init(current_entry.entries[current_entry.entries.length++], token.data, token.length))
				{
					for (auto& entry : adv_entries) free(entry);
					return false;
				}
				token.clear(); shift = {0};
			} else {
				if (!append_to_token(token, next, shift)) {
					for (auto& entry : adv_entries) free(entry);
					return false;
				}
			}
			break;
		}

		if (new_line) {
			current_entry.pos.line++;
			current_entry.pos.column = 1;
			new_line = false;
		} else current_entry.pos.column++;
		next = fgetc32(wrapper);
	}

	if (!feof(input)) {
		perror("Error occurred while reading morphology file");
		for (auto& entry : adv_entries) free(entry);
		return false;
	} else if (state != morphology_state::DEFAULT) {
		read_error("Unexpected end of input", current_entry.pos);
		for (auto& entry : adv_entries) free(entry);
		return false;
	}

	for (wikt_entry& entry : adv_entries) {
		if (!emit_adv_entry(m, entry, names)) {
			for (auto& entry : adv_entries) free(entry);
			return false;
		}
	}
	for (auto& entry : adv_entries) free(entry);
	return m.add_all_inflected_forms();
}

inline bool morphology_read(morphology_en& m,
		hash_map<string, unsigned int>& names,
		const char* filepath)
{
	FILE* in = fopen(filepath, "r");
	if (in == nullptr) {
		fprintf(stderr, "morphology_read ERROR: Unable to open '%s' for reading.\n", filepath);
		return false;
	} else if (!morphology_read(m, names, in)) {
		fclose(in); return false;
	}
	fclose(in);
	return true;
}

inline bool init_capitalization_map(
		morphology_en& m,
		hash_map<string, unsigned int>& names)
{
	array<pair<string, unsigned int>> to_capitalize(1024);
	for (const auto& entry : names) {
		if (entry.key.length == 0 || !islower(entry.key[0]))
			continue;
		if (!to_capitalize.ensure_capacity(to_capitalize.length + 1)
		 || !init(to_capitalize[to_capitalize.length].key, entry.key))
		{
			for (pair<string, unsigned int>& element : to_capitalize) free(element.key);
			return false;
		}
		to_capitalize[to_capitalize.length++].value = entry.value;
	}

	for (pair<string, unsigned int>& element : to_capitalize) {
		unsigned int capitalized_id;
		element.key[0] = toupper(element.key[0]);
		if (!get_token(element.key, capitalized_id, names)
		 || !m.add_capitalized_form(element.value, capitalized_id))
		{
			for (pair<string, unsigned int>& element : to_capitalize) free(element.key);
			return false;
		}
	}

	for (pair<string, unsigned int>& element : to_capitalize) free(element.key);
	return true;
}

#endif /* MORPHOLOGY_H_ */
