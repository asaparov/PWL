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
		free(dst);
		dst = new_dst;
		dst_length = new_length;
	}
	return true;
}

inline bool init_string_id_array(
		unsigned int*& dst, unsigned int& dst_length,
		const string* strs, unsigned int str_count,
		hash_map<string, unsigned int>& names)
{
	dst_length = str_count;
	if (str_count == 0) {
		dst = nullptr;
	} else {
		dst = (unsigned int*) malloc(sizeof(unsigned int) * str_count);
		if (dst == nullptr) {
			fprintf(stderr, "init_string_id_array ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str_count; i++) {
			if (!get_token(strs[i], dst[i], names)) {
				free(dst);
				return false;
			}
		}

		/* sort the IDs */
		insertion_sort(dst, dst_length);
		dst_length = unique(dst, dst_length);
	}
	return true;
}

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
	unsigned int* plural;
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

		return merge_arrays(plural, plural_count, src.plural, src.plural_count);
	}

	static inline void move(const noun_root& src, noun_root& dst) {
		dst.count = src.count;
		dst.plural = src.plural;
		dst.plural_count = src.plural_count;
	}

	static inline void free(noun_root& root) {
		if (root.plural_count != 0)
			core::free(root.plural);
	}
};

inline bool init(noun_root& root, countability count,
		const string* plural_forms, unsigned int plural_form_count,
		hash_map<string, unsigned int>& names)
{
	root.count = count;
	return init_string_id_array(root.plural, root.plural_count, plural_forms, plural_form_count, names);
}

enum class comparability {
	COMPARABLE,
	INCOMPARABLE,
	UNKNOWN,
	COMPARABLE_ONLY
};

struct adjective_root {
	comparability comp;
	pair<unsigned int, unsigned int>* inflected_forms; /* comparatives and superlatives, respectively */
	unsigned int inflected_form_count;

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
		dst.inflected_forms = src.inflected_forms;
		dst.inflected_form_count = src.inflected_form_count;
	}

	static inline void free(adjective_root& root) {
		if (root.inflected_form_count != 0)
			core::free(root.inflected_forms);
	}
};

typedef adjective_root adverb_root;

inline bool init(adjective_root& root, comparability comp,
		const pair<string, string>* inflected_forms,
		unsigned int inflected_form_count,
		hash_map<string, unsigned int>& names)
{
	root.comp = comp;
	root.inflected_form_count = inflected_form_count;
	if (inflected_form_count == 0) {
		root.inflected_forms = nullptr;
	} else {
		root.inflected_forms = (pair<unsigned int, unsigned int>*) malloc(sizeof(pair<unsigned int, unsigned int>) * inflected_form_count);
		if (root.inflected_forms == nullptr) {
			fprintf(stderr, "init ERROR: Insufficient memory for `adjective_root.inflected_forms`.\n");
			return false;
		}
		for (unsigned int i = 0; i < inflected_form_count; i++) {
			if (inflected_forms[i].key.length == 0) {
				root.inflected_forms[i].key = 0;
			} else if (!get_token(inflected_forms[i].key, root.inflected_forms[i].key, names)) {
				free(root.inflected_forms);
				return false;
			}
			if (inflected_forms[i].value.length == 0) {
				root.inflected_forms[i].value = 0;
			} else if (!get_token(inflected_forms[i].value, root.inflected_forms[i].value, names)) {
				free(root.inflected_forms);
				return false;
			}
		}

		/* remove any entries where both the comparative and superlative are empty */
		for (unsigned int i = 0; i < root.inflected_form_count; i++) {
			if (root.inflected_forms[i].key == 0 && root.inflected_forms[i].value == 0) {
				root.inflected_forms[i] = root.inflected_forms[root.inflected_form_count - 1];
				root.inflected_form_count--;
				i--;
			}
		}

		if (root.inflected_form_count == 0) {
			free(root.inflected_forms);
			root.inflected_forms = nullptr;
		} else {
			/* sort the inflected forms */
			insertion_sort(root.inflected_forms, root.inflected_form_count);
			root.inflected_form_count = unique(root.inflected_forms, root.inflected_form_count);
		}
	}
	return true;
}

struct verb_root {
	unsigned int* present_3sg;
	unsigned int present_3sg_count;
	unsigned int* present_participle;
	unsigned int present_participle_count;
	unsigned int* simple_past;
	unsigned int simple_past_count;
	unsigned int* past_participle;
	unsigned int past_participle_count;

	inline bool merge(const verb_root& src) {
		return merge_arrays(present_3sg, present_3sg_count, src.present_3sg, src.present_3sg_count)
			&& merge_arrays(present_participle, present_participle_count, src.present_participle, src.present_participle_count)
			&& merge_arrays(simple_past, simple_past_count, src.simple_past, src.simple_past_count)
			&& merge_arrays(past_participle, past_participle_count, src.past_participle, src.past_participle_count);
	}

	static inline void move(const verb_root& src, verb_root& dst) {
		dst.present_3sg = src.present_3sg;
		dst.present_3sg_count = src.present_3sg_count;
		dst.present_participle = src.present_participle;
		dst.present_participle_count = src.present_participle_count;
		dst.simple_past = src.simple_past;
		dst.simple_past_count = src.simple_past_count;
		dst.past_participle = src.past_participle;
		dst.past_participle_count = src.past_participle_count;
	}

	static inline void free(verb_root& root) {
		if (root.present_3sg_count != 0)
			core::free(root.present_3sg);
		if (root.present_participle_count != 0)
			core::free(root.present_participle);
		if (root.simple_past_count != 0)
			core::free(root.simple_past);
		if (root.past_participle_count != 0)
			core::free(root.past_participle);
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

	ANY
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
	SUPERLATIVE
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

struct inflected_noun {
	unsigned int root;
	grammatical_num number;
};

struct inflected_adjective {
	unsigned int root;
	grammatical_comparison comp;
};

typedef inflected_adjective inflected_adverb;

struct inflected_verb {
	unsigned int root;
	grammatical_person person;
	grammatical_num number;
	grammatical_mood mood;
	grammatical_tense tense;
};

struct morphology_en {
	hash_map<unsigned int, noun_root> nouns;
	hash_map<unsigned int, adjective_root> adjectives;
	hash_map<unsigned int, adverb_root> adverbs;
	hash_map<unsigned int, verb_root> verbs;

	hash_map<unsigned int, array<inflected_noun>> inflected_nouns;
	hash_map<unsigned int, array<inflected_adjective>> inflected_adjectives;
	hash_map<unsigned int, array<inflected_adverb>> inflected_adverbs;
	hash_map<unsigned int, array<inflected_verb>> inflected_verbs;

	hash_map<unsigned int, unsigned int> decapitalization_map;

	unsigned int MORE_COMPARATIVE_ID;
	unsigned int MOST_SUPERLATIVE_ID;
	unsigned int FURTHER_COMPARATIVE_ID;
	unsigned int FURTHEST_SUPERLATIVE_ID;

	static string MORE_COMPARATIVE_STRING;
	static string MOST_SUPERLATIVE_STRING;
	static string FURTHER_COMPARATIVE_STRING;
	static string FURTHEST_SUPERLATIVE_STRING;

	morphology_en() : nouns(1024), adjectives(1024), adverbs(1024), verbs(1024),
		inflected_nouns(2048), inflected_adjectives(2048), inflected_adverbs(2048),
		inflected_verbs(2048), decapitalization_map(2048)
	{ }

	~morphology_en() {
		for (auto entry : nouns)
			free(entry.value);
		for (auto entry : adjectives)
			free(entry.value);
		for (auto entry : adverbs)
			free(entry.value);
		for (auto entry : verbs)
			free(entry.value);
		for (auto entry : inflected_nouns)
			free(entry.value);
		for (auto entry : inflected_adjectives)
			free(entry.value);
		for (auto entry : inflected_adverbs)
			free(entry.value);
		for (auto entry : inflected_verbs)
			free(entry.value);
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
		 || !get_token("been", been, names))
			return false;

		/* add the infinitive form */
		if (!add_inflected_form(inflected_verbs, be, {be, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::BARE_INFINITIVE, grammatical_tense::ANY}))
			return false;

		/* add the subjunctive forms */
		if (!add_inflected_form(inflected_verbs, be, {be, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::SUBJUNCTIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, were, {be, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::SUBJUNCTIVE, grammatical_tense::PAST}))
			return false;

		/* add the indicative verb forms */
		if (!add_inflected_form(inflected_verbs, am, {be, grammatical_person::FIRST, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, are, {be, grammatical_person::SECOND, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, is, {be, grammatical_person::THIRD, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, was, {be, grammatical_person::FIRST_OR_THIRD, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PAST})
		 || !add_inflected_form(inflected_verbs, were, {be, grammatical_person::SECOND, grammatical_num::SINGULAR, grammatical_mood::INDICATIVE, grammatical_tense::PAST})
		 || !add_inflected_form(inflected_verbs, are, {be, grammatical_person::ANY, grammatical_num::PLURAL, grammatical_mood::INDICATIVE, grammatical_tense::PRESENT})
		 || !add_inflected_form(inflected_verbs, were, {be, grammatical_person::ANY, grammatical_num::PLURAL, grammatical_mood::INDICATIVE, grammatical_tense::PAST}))
			return false;

		/* add the participle forms */
		if (!add_inflected_form(inflected_verbs, being, {be, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PRESENT_PARTICIPLE, grammatical_tense::ANY})
		 || !add_inflected_form(inflected_verbs, been, {be, grammatical_person::ANY, grammatical_num::ANY, grammatical_mood::PAST_PARTICIPLE, grammatical_tense::ANY}))
			return false;

		return true;
	}

	inline bool add_capitalized_form(unsigned int decapitalized, unsigned int capitalized) {
		return decapitalization_map.put(capitalized, decapitalized);
	}

	inline bool decapitalize(unsigned int word_id, unsigned int& decapitalized_word_id) const {
		bool contains;
		decapitalized_word_id = decapitalization_map.get(word_id, contains);
		return contains;
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_noun_root(unsigned int root_id, noun_root& root) {
		return add_inflected_forms(inflected_nouns, root_id, root)
			&& add_root(nouns, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adjective_root(unsigned int root_id, adjective_root& root) {
		return add_inflected_forms(inflected_adjectives, root_id, root)
			&& add_root(adjectives, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adverb_root(unsigned int root_id, adverb_root& root) {
		return add_inflected_forms(inflected_adverbs, root_id, root)
			&& add_root(adverbs, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_verb_root(unsigned int root_id, verb_root& root) {
		return add_inflected_forms(inflected_verbs, root_id, root)
			&& add_root(verbs, root_id, root);
	}

private:
	template<typename T>
	static inline bool add_root(
			hash_map<unsigned int, T>& root_map,
			unsigned int root_id, T& root)
	{
		if (!root_map.check_size()) {
			free(root);
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
				free(root);
				return false;
			}
			free(root);
		}
		return true;
	}

	template<typename T>
	static inline bool add_inflected_form(
			hash_map<unsigned int, array<T>>& inflected_form_map,
			unsigned int new_inflected_form_id, T new_inflected_form)
	{
		bool contains; unsigned int bucket;
		array<T>& inflected_forms = inflected_form_map.get(new_inflected_form_id, contains, bucket);
		if (!contains) {
			if (!array_init(inflected_forms, 1))
				return false;
			inflected_form_map.table.keys[bucket] = new_inflected_form_id;
			inflected_form_map.table.size++;
		}
		return inflected_forms.add(new_inflected_form);
	}

	static inline bool add_inflected_forms(
			hash_map<unsigned int, array<inflected_noun>>& inflected_form_map,
			unsigned int root_id, const noun_root& root)
	{
		if (!inflected_form_map.check_size(inflected_form_map.table.size + root.plural_count + 1))
			return false;

		/* first add the singular form */
		if (root.count != countability::PLURAL_ONLY) {
			if (!add_inflected_form(inflected_form_map, root_id, {root_id, grammatical_num::SINGULAR}))
				return false;
		}

		/* add the plural forms */
		for (unsigned int i = 0; i < root.plural_count; i++) {
			if (!add_inflected_form(inflected_form_map, root.plural[i], {root_id, grammatical_num::PLURAL}))
				return false;
		}
		return true;
	}

	static inline bool add_inflected_forms(
			hash_map<unsigned int, array<inflected_adjective>>& inflected_form_map,
			unsigned int root_id, const adjective_root& root)
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
			if (root.inflected_forms[i].key != 0) {
				if (!add_inflected_form(inflected_form_map, root.inflected_forms[i].key, {root_id, grammatical_comparison::COMPARATIVE}))
					return false;
			}

			if (root.inflected_forms[i].value != 0) {
				if (!add_inflected_form(inflected_form_map, root.inflected_forms[i].value, {root_id, grammatical_comparison::SUPERLATIVE}))
					return false;
			}
		}
		return true;
	}

	static inline bool add_inflected_forms(
			hash_map<unsigned int, array<inflected_verb>>& inflected_form_map,
			unsigned int root_id, const verb_root& root)
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
};

string morphology_en::MORE_COMPARATIVE_STRING = "_more";
string morphology_en::MOST_SUPERLATIVE_STRING = "_most";
string morphology_en::FURTHER_COMPARATIVE_STRING = "_further";
string morphology_en::FURTHEST_SUPERLATIVE_STRING = "_furthest";

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

inline bool emit_entry(morphology_en& m,
		const string& root, unsigned int root_id,
		const array<string>& entry,
		hash_map<string, unsigned int>& names,
		position current)
{
	if (entry.length == 0) {
		read_error("Found an entry with no elements", current);
		return false;
	} else if (entry[0] == "n") {
		if (entry.length == 1) {
			/* the noun is countable and has a simple plural form "-s" */
			string plural(root.data, root.length);
			plural += "s";

			noun_root new_root;
			if (!init(new_root, countability::COUNTABLE, &plural, 1, names))
				return false;
			return m.add_noun_root(root_id, new_root);
		} else {
			array<string> plural_forms(entry.length - 1);
			countability count = countability::COUNTABLE;
			for (unsigned int i = 1; i < entry.length; i++) {
				if (entry[i] == "s") {
					if (!init(plural_forms[plural_forms.length], root)) {
						for (string& str : plural_forms) free(str);
						return false;
					}
					plural_forms[plural_forms.length++] += "s";
				} else if (entry[i] == "es") {
					if (!init(plural_forms[plural_forms.length], root)) {
						for (string& str : plural_forms) free(str);
						return false;
					}
					plural_forms[plural_forms.length++] += "es";
				} else if (entry[i] == "-") {
					count = countability::UNCOUNTABLE;
				} else if (entry[i] == "~") {
					count = countability::BOTH;
				} else if (entry[i] == "!") {
					count = countability::PLURAL_UNATTESTED;
				} else if (entry[i] == "?") {
					count = countability::PLURAL_UNKNOWN;
				} else if (starts_with(entry[i], "head=")) {
					/* ignore these for now */
					for (string& str : plural_forms) free(str);
					return true;
				} else {
					if (entry[i].index_of('=') < entry[i].length || entry[i].index_of('[') < entry[i].length) {
						for (string& str : plural_forms) free(str);
						return true;
					}
					plural_forms[plural_forms.length++] = entry[i];
				}
			}

			noun_root new_root;
			if (!init(new_root, count, plural_forms.data, plural_forms.length, names))
				return false;
			for (string& str : plural_forms) free(str);
			return m.add_noun_root(root_id, new_root);
		}
	} else if (entry[0] == "pl") {
		/* the noun is plural only; we ignore these for now */
		return true;
	} else if (entry[0] == "adj" || entry[0] == "adv") {
		if (entry.length == 1) {
			/* the comparative is formed with "more" and the superlative with "most" */
			pair<string, string>& inflected_forms = *((pair<string, string>*) alloca(sizeof(pair<string, string>)));
			inflected_forms.key = morphology_en::MORE_COMPARATIVE_STRING;
			inflected_forms.value = morphology_en::MOST_SUPERLATIVE_STRING;

			adjective_root new_root;
			if (!init(new_root, comparability::COMPARABLE, &inflected_forms, 1, names)) {
				free(inflected_forms);
				return false;
			}
			free(inflected_forms);
			if (entry[0] == "adj")
				return m.add_adjective_root(root_id, new_root);
			else return m.add_adverb_root(root_id, new_root);
		} else {
			array<pair<string, string>> inflected_forms(entry.length - 1);
			array_map<unsigned int, string> superlative_forms(entry.length - 1);
			comparability comp = comparability::COMPARABLE;
			for (unsigned int i = 1; i < entry.length; i++) {
				if (entry[i] == "more" && root != "many" && root != "most") {
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
				} else if (entry[i] == "further" && root != "far") {
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
				} else if (entry[i] == "er") {
					string stem(root.data, root.length);
					if (root.length > 0 && root[root.length - 1] == 'e') {
						stem.length--;
					} else if (root.length > 2 && root[root.length - 2] == 'e' && root[root.length - 1] == 'y'
							&& (root.length < 3 || !is_vowel(root[root.length - 3])))
					{
						stem.length--;
						stem[stem.length - 1] = 'i';
					} else if (root.length > 1 && root[root.length - 1] == 'y'
							&& (root.length < 2 || !is_vowel(root[root.length - 2])))
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
				} else if (entry[i] == "-") {
					comp = comparability::INCOMPARABLE;
				} else if (entry[i] == "?") {
					comp = comparability::UNKNOWN;
				} else if (entry[i] == "+") {
					comp = comparability::COMPARABLE_ONLY;
				} else {
					unsigned int comparative_index;
					const char* superlative; unsigned int superlative_length;
					if (get_parameter(entry[i], "sup", "=", comparative_index, superlative, superlative_length)) {
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
						if (entry[i].index_of('=') < entry[i].length || entry[i].index_of('[') < entry[i].length) {
							for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
							for (auto entry : superlative_forms) free(entry.value);
							return true;
						}

						if (!init(inflected_forms[inflected_forms.length].key, entry[i])) {
							for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
							for (auto entry : superlative_forms) free(entry.value);
							return false;
						}

						if (entry[i].length > 2 && entry[i][entry[i].length - 2] == 'e' && entry[i][entry[i].length - 1] == 'r') {
							string superlative(entry[i].data, entry[i].length);
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

			adjective_root new_root;
			if (!init(new_root, comp, inflected_forms.data, inflected_forms.length, names))
				return false;
			for (auto& pair : inflected_forms) { free(pair.key); free(pair.value); }
			if (entry[0] == "adj")
				return m.add_adjective_root(root_id, new_root);
			else return m.add_adverb_root(root_id, new_root);
		}
	} else if (entry[0] == "v") {
		if (entry.length == 1) {
			/* the noun is regular, i.e. the present 3rd person singular is
			   formed by adding "-s", the present participle by adding "-ing",
			   and the past by adding "-ed" */
			string present_3sg(root.data, root.length);
			present_3sg += "s";
			string present_participle(root.data, root.length);
			present_participle += "ing";
			string past(root.data, root.length);
			past += "ed";

			verb_root new_root;
			if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
				return false;
			return m.add_verb_root(root_id, new_root);
		}

		static constexpr const char* obsolete_str = "obsolete";
		static unsigned int obsolete_str_length = strlen(obsolete_str);

		/* first process the flags */
		array<string> args(entry.length - 1);
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
		for (unsigned int i = 1; i < entry.length; i++) {
			unsigned int qualifier_index; const char* qualifier; unsigned int qualifier_length;
			if (entry[i].index_of('[') < entry[i].length) { free_strings(); return true; }
			if (get_parameter(entry[i], "head", "=", qualifier_index, qualifier, qualifier_length)) { free_strings(); return true; }

			if (get_parameter(entry[i], "pres_3sg", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(pres_3sg, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry[i], "pres_ptc", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(pres_ptc, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry[i], "past_ptc", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(past_ptc, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry[i], "past", "=", qualifier_index, qualifier, qualifier_length)) {
				if (!inflect_verb(past, qualifier_index, qualifier, qualifier_index)) { free_strings(); return false; }
				continue;
			} else if (get_parameter(entry[i], "pres_3sg", "_qual=", qualifier_index, qualifier, qualifier_length)) {
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
			} else if (get_parameter(entry[i], "pres_ptc", "_qual=", qualifier_index, qualifier, qualifier_length)) {
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
			} else if (get_parameter(entry[i], "past_ptc", "_qual=", qualifier_index, qualifier, qualifier_length)) {
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
			} else if (get_parameter(entry[i], "past", "_qual=", qualifier_index, qualifier, qualifier_length)) {
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
			unsigned int eq_index = entry[i].index_of('=');
			if (eq_index < entry[i].length) {
				fprintf(stderr, "WARNING at %u:%u: Unrecognized parameter in entry '", current.line, current.column);
				print(entry[i], stderr); print("'.\n", stderr);
				continue;
			}
#endif

			if (!init(args[args.length], entry[i])) {
				free_strings(); return false;
			}
			args.length++;
		}

		bool first_empty = (args.length < 1 || args[0].length == 0);
		bool second_empty = (args.length < 2 || args[1].length == 0);
		bool third_empty = (args.length < 3 || args[2].length == 0);
		if (!first_empty && second_empty && third_empty) {
			if (args[0] == "es") {
				if (!inflect_verb(pres_3sg, root, "es")
				 || !inflect_verb(pres_ptc, root, "ing")
				 || !inflect_verb(past, root, "ed"))
				{ free_strings(); return false; }
			} else if (args[0] == "ies") {
				if (root.length == 0 || root[root.length - 1] != 'y') {
					fprintf(stderr, "ERROR at %u:%u: Verb root '", current.line, current.column);
					print(root, stderr); fprintf(stderr, "' does not end in 'y'.\n");
					return false;
				}

				string stem(root.data, root.length - 1);
				if (!inflect_verb(pres_3sg, stem, "ies")
				 || !inflect_verb(pres_ptc, stem, "ying")
				 || !inflect_verb(past, stem, "ied"))
				{ free_strings(); return false; }
			} else if (args[0] == "d") {
				if (!inflect_verb(pres_3sg, root, "s")
				 || !inflect_verb(pres_ptc, root, "ing")
				 || !inflect_verb(past, root, "d"))
				{ free_strings(); return false; }
			} else {
				if (!inflect_verb(pres_3sg, root, "s")
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
				if (!inflect_verb(pres_3sg, root, "s")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "ed"))
				{ free_strings(); return false; }
			} else if (args.length > 1 && args[1] == "d") {
				if (!inflect_verb(pres_3sg, root, "s")
				 || !inflect_verb(pres_ptc, args[0], "ing")
				 || !inflect_verb(past, args[0], "d"))
				{ free_strings(); return false; }
			} else {
				if (first_empty) {
					if (!inflect_verb(pres_3sg, root, "s")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, args[0])) { free_strings(); return false; }
				}
				if (second_empty) {
					if (!inflect_verb(pres_ptc, root, "ing")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_ptc, args[1])) { free_strings(); return false; }
				}
				if (!inflect_verb(past, root, "ed")) { free_strings(); return false; }
			}
		} else {
			if (args[2] == "es") {
				if (!inflect_verb(pres_3sg, args[0], args[1], "es")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing")
				 || !inflect_verb(past, args[0], args[1], "ed"))
				{ free_strings(); return false; }
			} else if (args[2] == "ing") {
				if (!inflect_verb(pres_3sg, root, "s")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing"))
				{ free_strings(); return false; }

				if (args[1] == "y") {
					if (!inflect_verb(past, root, "d")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(past, args[0], args[1], "ed")) { free_strings(); return false; }
				}
			} else if (args[2] == "ed") {
				if (args[1] == "i") {
					if (!inflect_verb(pres_3sg, args[0], args[1], "es")
					 || !inflect_verb(pres_ptc, root, "ing"))
					{ free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, root, "s")
					 || !inflect_verb(pres_ptc, args[0], args[1], "ing"))
					{ free_strings(); return false; }
				}
				if (!inflect_verb(past, args[0], args[1], "ed")) { free_strings(); return false; }
			} else if (args[2] == "d") {
				if (!inflect_verb(pres_3sg, root, "s")
				 || !inflect_verb(pres_ptc, args[0], args[1], "ing")
				 || !inflect_verb(past, args[0], args[1], "d"))
				{ free_strings(); return false; }
			} else {
				if (first_empty) {
					if (!inflect_verb(pres_3sg, root, "s")) { free_strings(); return false; }
				} else {
					if (!inflect_verb(pres_3sg, args[0])) { free_strings(); return false; }
				}
				if (second_empty) {
					if (!inflect_verb(pres_ptc, root, "ing")) { free_strings(); return false; }
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
		if (!init(new_root, pres_3sg.data, pres_3sg.length, pres_ptc.data, pres_ptc.length, past.data, past.length, past_ptc.data, past_ptc.length, names)) {
			free_strings();
			return false;
		}
		free_strings();
		return m.add_verb_root(root_id, new_root);
	} else {
		read_error("Unrecognized part of speech", current);
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
	position current = start;
	morphology_state state = morphology_state::DEFAULT;
	array<char> token = array<char>(1024);

	string current_root("");
	unsigned int current_root_id = 0;
	array<string> current_entry(16);

	std::mbstate_t shift = {0};
	wint_t next = fgetwc(input);
	bool new_line = false;
	while (next != WEOF) {
		switch (state) {
		case morphology_state::DEFAULT:
			if (next == '\t') {
				string new_root(token.data, token.length);
				swap(current_root, new_root);
				if (!get_token(current_root, current_root_id, names)) {
					for (string& str : current_entry) { free(str); } return false;
				}
				state = morphology_state::ENTRY;
				token.clear(); shift = {0};
			} else if (next == '\n') {
				fprintf(stderr, "WARNING: Found root with no entries on line %u.\n", current.line);
				state = morphology_state::DEFAULT;
				new_line = true;
				token.clear(); shift = {0};
			} else {
				if (!append_to_token(token, next, shift)) {
					for (string& str : current_entry) { free(str); } return false;
				}
			}
			break;

		case morphology_state::ENTRY:
			if (next == '\t' || next == '\n') {
				if (!current_entry.ensure_capacity(current_entry.length + 1)
				 || !init(current_entry[current_entry.length++], token.data, token.length))
				{
					for (string& str : current_entry) { free(str); } return false;
				} else if (!emit_entry(m, current_root, current_root_id, current_entry, names, current)) {
					for (string& str : current_entry) { free(str); } return false;
				}
				for (string& str : current_entry) { free(str); }
				current_entry.clear();
				token.clear(); shift = {0};
				if (next == '\n') {
					state = morphology_state::DEFAULT;
					new_line = true;
				}
			} else if (next == '\\') {
				if (!current_entry.ensure_capacity(current_entry.length + 1)
				 || !init(current_entry[current_entry.length++], token.data, token.length))
				{
					for (string& str : current_entry) { free(str); } return false;
				}
				token.clear(); shift = {0};
			} else {
				if (!append_to_token(token, next, shift)) {
					for (string& str : current_entry) { free(str); } return false;
				}
			}
			break;
		}

		if (new_line) {
			current.line++;
			current.column = 1;
			new_line = false;
		} else current.column++;
		next = fgetwc(input);
	}

	for (string& str : current_entry) { free(str); }
	if (!feof(input)) {
		perror("Error occurred while reading morphology file");
		return false;
	} else if (state != morphology_state::DEFAULT) {
		read_error("Unexpected end of input", current);
		return false;
	}
	return true;
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