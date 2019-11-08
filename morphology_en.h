#ifndef MORPHOLOGY_H_
#define MORPHOLOGY_H_

#include <core/map.h>

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

struct morphology_en {
	hash_map<unsigned int, noun_root> nouns;
	hash_map<unsigned int, adjective_root> adjectives;
	hash_map<unsigned int, adverb_root> adverbs;
	hash_map<unsigned int, verb_root> verbs;

	unsigned int MORE_COMPARATIVE_ID;
	unsigned int MOST_SUPERLATIVE_ID;
	unsigned int FURTHER_COMPARATIVE_ID;
	unsigned int FURTHEST_SUPERLATIVE_ID;

	static string MORE_COMPARATIVE_STRING;
	static string MOST_SUPERLATIVE_STRING;
	static string FURTHER_COMPARATIVE_STRING;
	static string FURTHEST_SUPERLATIVE_STRING;

	morphology_en() : nouns(1024), adjectives(1024), adverbs(1024), verbs(1024) { }

	~morphology_en() {
		for (auto entry : nouns)
			free(entry.value);
		for (auto entry : adjectives)
			free(entry.value);
		for (auto entry : adverbs)
			free(entry.value);
		for (auto entry : verbs)
			free(entry.value);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_noun_root(unsigned int root_id, noun_root& root) {
		return add_root(nouns, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adjective_root(unsigned int root_id, adjective_root& root) {
		return add_root(adjectives, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_adverb_root(unsigned int root_id, adverb_root& root) {
		return add_root(adverbs, root_id, root);
	}

	/* NOTE: this function takes ownership of the memory of `root` */
	bool add_verb_root(unsigned int root_id, verb_root& root) {
		return add_root(verbs, root_id, root);
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

inline bool qualified_plural(
		const string& src, unsigned int& plural_index,
		const char*& qualifier, unsigned int& qualifier_length)
{
	if (src.length < 8) return false;
	if (src[0] != 'p' || src[1] != 'l') return false;
	unsigned int q_index = index_of('q', src.data + 2, src.length - 2) + 2;
	if (src.length - q_index < 5) return false;
	if (!parse_uint(string(src.data + 2, q_index - 2), plural_index))
		return false;
	if (src[q_index + 1] != 'u' || src[q_index + 2] != 'a' || src[q_index + 3] != 'l' || src[q_index + 4] != '=')
		return false;
	qualifier = src.data + q_index + 5;
	qualifier_length = src.length - q_index - 5;
	return true;
}

inline bool specified_superlative(
		const string& src, unsigned int& comparative_index,
		const char*& superlative, unsigned int& superlative_length)
{
	if (src.length < 4) return false;
	if (src[0] != 's' || src[1] != 'u' || src[2] != 'p') return false;
	unsigned int eq_index = index_of('=', src.data + 3, src.length - 3) + 3;
	if (eq_index == src.length) return false;
	if (eq_index == 3) {
		comparative_index = 1;
	} else if (!parse_uint(string(src.data + 3, eq_index - 3), comparative_index)) {
		return false;
	}
	comparative_index--;

	superlative = src.data + eq_index + 1;
	superlative_length = src.length - eq_index - 1;
	return true;
}

inline bool is_vowel(char c) {
	return (c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u');
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
					if (specified_superlative(entry[i], comparative_index, superlative, superlative_length)) {
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

		/* first process the flags */
		array<string> args(entry.length - 1);
		for (unsigned int i = 1; i < entry.length; i++) {
			if (entry[i].index_of('[') < entry[i].length) {
				for (string& str : args) free(str);
				return true;
			}

			unsigned int eq_index = entry[i].index_of('=');
			/* TODO: continue from here */
		}

		if (entry.length == 2) {
			if (entry[1] == "es") {
				string present_3sg(root.data, root.length);
				present_3sg += "es";
				string present_participle(root.data, root.length);
				present_participle += "ing";
				string past(root.data, root.length);
				past += "ed";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[1] == "ies") {
				if (root.length == 0 || root[root.length - 1] != 'y') {
					fprintf(stderr, "ERROR at %u:%u: Verb root '", current.line, current.column);
					print(root, stderr); fprintf(stderr, "' does not end in 'y'.\n");
					return false;
				}

				string present_3sg(root.data, root.length - 1);
				present_3sg += "ies";
				string present_participle(root.data, root.length - 1);
				present_participle += "ying";
				string past(root.data, root.length - 1);
				past += "ied";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[1] == "d") {
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
			} else {
				string present_3sg(root.data, root.length);
				present_3sg += "s";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += "ed";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			}
		} else if (entry.length == 3) {
			if (entry[2] == "es") {
				string present_3sg(entry[1].data, entry[1].length);
				present_3sg += "es";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += "ed";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[2] == "ies") {
				string present_3sg(entry[1].data, entry[1].length);
				present_3sg += "ies";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += "ying";
				string past(entry[1].data, entry[1].length);
				past += "ied";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[2] == "ing" || entry[2] == "ed") {
				string present_3sg(root.data, root.length);
				present_3sg += "s";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += "ed";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[2] == "d") {
				string present_3sg(root.data, root.length);
				present_3sg += "s";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += "d";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else {
				string past(root.data, root.length);
				past += "ed";

				verb_root new_root;
				if (!init(new_root, entry.data + 1, 1, entry.data + 2, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			}
		} else {
			if (entry[3] == "es") {
				string present_3sg(entry[1].data, entry[1].length);
				present_3sg += entry[2];
				present_3sg += "es";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += entry[2];
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += entry[2];
				past += "ed";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else if (entry[3] == "ing") {
				string present_3sg(root.data, root.length);
				present_3sg += "s";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += entry[2];
				present_participle += "ing";
				if (entry[2] == "y") {
					string past(root.data, root.length);
					past += "d";

					verb_root new_root;
					if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
						return false;
					return m.add_verb_root(root_id, new_root);
				} else {
					string past(entry[1].data, entry[1].length);
					past += entry[2];
					past += "ed";

					verb_root new_root;
					if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
						return false;
					return m.add_verb_root(root_id, new_root);
				}
			} else if (entry[3] == "ed") {
				string past(entry[1].data, entry[1].length);
				past += entry[2];
				past += "ed";

				if (entry[2] == "i") {
					string present_3sg(entry[1].data, entry[1].length);
					present_3sg += entry[2];
					present_3sg += "es";
					string present_participle(root.data, root.length);
					present_participle += "ing";

					verb_root new_root;
					if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
						return false;
					return m.add_verb_root(root_id, new_root);
				} else {
					string present_3sg(root.data, root.length);
					present_3sg += "s";
					string present_participle(entry[1].data, entry[1].length);
					present_participle += entry[2];
					present_participle += "ing";

					verb_root new_root;
					if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
						return false;
					return m.add_verb_root(root_id, new_root);
				}
			} else if (entry[3] == "d") {
				string present_3sg(root.data, root.length);
				present_3sg += "s";
				string present_participle(entry[1].data, entry[1].length);
				present_participle += entry[2];
				present_participle += "ing";
				string past(entry[1].data, entry[1].length);
				past += entry[2];
				past += "d";

				verb_root new_root;
				if (!init(new_root, &present_3sg, 1, &present_participle, 1, &past, 1, &past, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			} else {
				verb_root new_root;
				if (!init(new_root, entry.data + 1, 1, entry.data + 2, 1, entry.data + 3, 1, entry.data + 3, 1, names))
					return false;
				return m.add_verb_root(root_id, new_root);
			}
		}
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
				if (!get_token(new_root, current_root_id, names)) {
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

	for (string& str : current_entry) { free(str); } return false;
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

#endif /* MORPHOLOGY_H_ */
