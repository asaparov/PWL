#include "xml.h"

#include <locale.h>

enum class wikt_xml_element {
	TITLE,
	TEXT,
	OTHER
};

struct wikt_xml_reader {
	wikt_xml_element element;
	string current_word;

	wikt_xml_reader() : element(wikt_xml_element::OTHER), current_word("") { }
};

unsigned int substring(
		const char* str, unsigned int str_length,
		const char* pattern, unsigned int pattern_length)
{
#if !defined(NDEBUG)
	if (pattern_length == 0)
		fprintf(stderr, "substring WARNING: `pattern` is empty.\n");
#endif

	for (unsigned int i = 0; i + pattern_length <= str_length; i++) {
		bool match = true;
		for (unsigned int j = 0; j < pattern_length; j++) {
			if (str[i + j] != pattern[j]) {
				match = false;
				break;
			}
		}

		if (match) return i;
	}

	return str_length;
}

inline bool tokenize(array<string>& tokens, const char* str, unsigned int str_length) {
	if (str[0] == '}') {
		return true;
	} else if (str[0] != '|') {
		fprintf(stderr, "tokenize ERROR: Expected either '}' or '|'.\n");
		return false;
	}

	unsigned int start = 1;
	bool found_end = false;
	for (unsigned int i = start; i < str_length; i++) {
		if (str[i] == '|') {
			if (!tokens.add(string(str + start, i - start))) return false;
			start = i + 1;
		} else if (str[i] == '}') {
			if (!tokens.add(string(str + start, i - start))) return false;
			start = i + 1;
			found_end = true;
			break;
		}
	}

	if (!found_end && !tokens.add(string(str + start, str_length - start)))
		return false;
	return true;
}

inline bool emit_text(const array<char>& text, const position& start, const position& end, wikt_xml_reader& reader) {
	if (reader.element == wikt_xml_element::TITLE) {
		string new_word(text.data, text.length);
		swap(new_word, reader.current_word);
	} else if (reader.element == wikt_xml_element::TEXT) {
		if (reader.current_word.index_of(' ') < reader.current_word.length
		 || reader.current_word.index_of(':') < reader.current_word.length)
			return true;

		static string ENGLISH_HEADER = "==English==";
		unsigned int english_section = substring(text.data, text.length, ENGLISH_HEADER.data, ENGLISH_HEADER.length);
		if (english_section == text.length)
			return true;

		bool has_subsections = false;
		static string SUBSECTION_PATTERNS[] = { "{{en-noun", "{{en-plural noun", "{{en-verb", "{{en-adjective", "{{en-adj", "{{en-adverb", "{{en-adv", "{{head" };
		static string POS_NAMES[] = { "n", "pl", "v", "adj", "adj", "adv", "adv", "" };
		for (unsigned int i = english_section + ENGLISH_HEADER.length; i < text.length; i++) {
			/* check if we reached the end of the English section */
			if (i + 3 < text.length && text[i] != '=' && text[i + 1] == '=' && text[i + 2] == '=' && text[i + 3] != '=')
				break;

			/* check if any of the patterns match at this position */
			for (unsigned int k = 0; k < array_length(SUBSECTION_PATTERNS); k++) {
				if (i + SUBSECTION_PATTERNS[k].length >= text.length) continue;
				bool match = true;
				for (unsigned int j = 0; j < SUBSECTION_PATTERNS[k].length; j++) {
					if (text[i + j] != SUBSECTION_PATTERNS[k][j]) {
						match = false;
						break;
					}
				}

				if (match) {
					if (k == array_length(SUBSECTION_PATTERNS) - 1) {
						if (text[i + SUBSECTION_PATTERNS[k].length] != '|')
							break;

						/* we matched the "head" pattern */
						/* parse the tokens separated by `|` until `}` */
						array<string> tokens(4);
						unsigned int offset = i + SUBSECTION_PATTERNS[k].length;
						if (!tokenize(tokens, text.data + offset, text.length - offset))
							return false;

						if (tokens.length >= 2 && tokens[0] == "en") {
							if (tokens[1] == "numeral" || tokens[1] == "abbreviation"
							 || tokens[1] == "noun form" || tokens[1] == "verb form"
							 || tokens[1] == "phrase" || tokens[1] == "determiner"
							 || tokens[1] == "proper noun" || tokens[1] == "suffix"
							 || tokens[1] == "initialism" || tokens[1] == "misspelling"
							 || tokens[1] == "acronym" || tokens[1] == "pronoun"
							 || tokens[1] == "contraction" || tokens[1] == "comparative adjective"
							 || tokens[1] == "superlative adjective" || tokens[1] == "interjection"
							 || tokens[1] == "particle" || tokens[1] == "infix"
							 || tokens[1] == "contractions" || tokens[1] == "symbol"
							 || tokens[1] == "article" || tokens[1] == "prefix"
							 || tokens[1] == "comparative adverb" || tokens[1] == "superlative adverb"
							 || tokens[1] == "suffixes" || tokens[1] == "preposition"
							 || tokens[1] == "misspellings" || tokens[1] == "number"
							 || tokens[1] == "prepositional phrase" || tokens[1] == "proper nouns"
							 || tokens[1] == "conjunction" || tokens[1] == "prefixes"
							 || tokens[1] == "postposition" || tokens[1] == "obsolete verb form"
							 || tokens[1] == "abbreviations" || tokens[1] == "proper noun plural form"
							 || tokens[1] == "suffix forms" || tokens[1] == "interfix"
							 || tokens[1] == "adjective form" || tokens[1] == "proper noun form"
							 || tokens[1] == "noun forms" || tokens[1] == "affix"
							 || tokens[1] == "numeral symbol" || tokens[1] == "pronominal adverbs"
							 || tokens[1] == "verb forms" || tokens[1] == "noun plural form"
							 || tokens[1] == " abbreviation" || tokens[1] == "suffix form"
							 || tokens[1] == "verb")
							{
								for (string& token : tokens) free(token);
								break;
							}
						}

						const char* pos_name = nullptr;
						if (tokens.length >= 2 && tokens[0] == "en") {
							if (tokens[1] == "adjective" || tokens[1] == "adjectives") {
								if (tokens.length == 2) {
									pos_name = "adj";
								} else if (tokens.length == 3 && tokens[2].length >= 4 && tokens[2][0] == 'h' && tokens[2][1] == 'e' && tokens[2][2] == 'a' && tokens[2][3] == 'd') {
									pos_name = "adj";
								}
							}
						}

						if (pos_name == nullptr) {
							for (string& token : tokens) free(token);
							break;
						}

						if (!has_subsections) {
							print(reader.current_word, stdout);
							has_subsections = true;
						}
						print('\t', stdout); print(pos_name, stdout);
						for (unsigned int i = 2; i < tokens.length; i++) {
							print('\\', stdout); print(tokens[i], stdout);
						}
						for (string& token : tokens) free(token);
						break;
					}

					if (!has_subsections) {
						print(reader.current_word, stdout);
						has_subsections = true;
					}

					/* parse the tokens separated by `|` until `}` */
					array<string> tokens(4);
					unsigned int offset = i + SUBSECTION_PATTERNS[k].length;
					if (!tokenize(tokens, text.data + offset, text.length - offset))
						return false;

					print('\t', stdout); print(POS_NAMES[k], stdout);
					for (const string& token : tokens) {
						print('\\', stdout); print(token, stdout);
					}
					for (string& token : tokens) free(token);
					break;
				}
			}
		}

		if (has_subsections)
			print('\n', stdout);
	}
	return true;
}

inline bool emit_start_element(const xml_element& element, position start, position end, wikt_xml_reader& reader) {
	if (element.name == "title") {
		reader.element = wikt_xml_element::TITLE;
	} else if (element.name == "text") {
		reader.element = wikt_xml_element::TEXT;
	}
	return true;
}

inline bool emit_end_element(const xml_element& element, position start, position end, wikt_xml_reader& reader) {
	reader.element = wikt_xml_element::OTHER;
	return true;
}

#include "morphology_en.h"

inline bool test_morphology() {
	hash_map<string, unsigned int> names(64);
	morphology_en m;
	bool result = morphology_read(m, names, "english.morph");
	for (auto entry : names) free(entry.key);
	fprintf(stderr, "Read %u noun roots.\n", m.nouns.table.size);
	fprintf(stderr, "Read %u adjective roots.\n", m.adjectives.table.size);
	fprintf(stderr, "Read %u adverb roots.\n", m.adverbs.table.size);
	fprintf(stderr, "Read %u verb roots.\n", m.verbs.table.size);
	return result;
}

int main(int argc, const char** argv) {
	setlocale(LC_ALL, "en_US.UTF-8");
test_morphology();
return EXIT_SUCCESS;
	FILE* in = fopen("/home/asaparov/Desktop/enwiktionary-20191101-pages-articles.xml", "rb");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open XML file for reading.\n");
		return false;
	}

	wikt_xml_reader reader;
	if (!xml_parse<true>(in, reader))
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
