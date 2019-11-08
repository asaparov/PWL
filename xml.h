#ifndef XML_H_
#define XML_H_

#include <core/array.h>
#include <core/lex.h>

using namespace core;

enum class xml_state {
	DEFAULT,
	START_TAG_NAME,
	START_TAG,
	ATTRIBUTE_NAME,
	ATTRIBUTE_VALUE_SINGLE,
	ATTRIBUTE_VALUE_DOUBLE,
	END_TAG_NAME,
	END_TAG,
	REFERENCE
};

inline bool is_whitespace(wint_t c) {
	return (c == ' ' || c == '\t' || c == '\r' || c == '\n');
}

inline bool is_whitespace_only(const array<char>& string) {
	for (char c : string)
		if (!is_whitespace(c)) return false;
	return true;
}

inline bool is_valid_name_start(wint_t c) {
	return (c == ':' || (c >= 'A' && c <= 'Z') || c == '_' || (c >= 'a' && c <= 'z')
		|| (c >= 0x20 && c <= 0xD6) || (c >= 0xD8 && c <= 0xF6)
		|| (c >= 0xF8 && c <= 0x2FF) || (c >= 0x370 && c <= 0x37D)
		|| (c >= 0x37F && c <= 0x1FFF) || (c >= 0x200C && c <= 0x200D)
		|| (c >= 0x2070 && c <= 0x218F) || (c >= 0x2C00 && c <= 0x2FEF)
		|| (c >= 0x3001 && c <= 0xD7FF) || (c >= 0xF900 && c <= 0xFDCF)
		|| (c >= 0xFDF0 && c <= 0xFFFD) || (c >= 0x10000 && c <= 0xEFFFF));
}

inline bool is_valid_name_char(wint_t c) {
	return (is_valid_name_start(c) || c == '-' || c == '.' || (c >= '0' && c <= '9')
		|| c == 0xB7 || (c >= 0x0300 && c <= 0x036F) || (c >= 0x203F && c <= 0x2040));
}

bool expand_reference(array<char>& reference_name, array<char>& text, std::mbstate_t& shift, position current) {
	if (compare_strings(reference_name, "lt")) {
		return append_to_token(text, '<', shift);
	} else if (compare_strings(reference_name, "gt")) {
		return append_to_token(text, '>', shift);
	} else if (compare_strings(reference_name, "amp")) {
		return append_to_token(text, '&', shift);
	} else if (compare_strings(reference_name, "quot")) {
		return append_to_token(text, '"', shift);
	} else if (compare_strings(reference_name, "apos")) {
		return append_to_token(text, '\'', shift);
	} else if (reference_name.length != 0 && reference_name[0] == '#') {
		unsigned int codepoint;
		reference_name.data++;
		reference_name.length--;
		if (!parse_uint(reference_name, codepoint, 16)) {
			read_error("Unable to interpret Unicode codepoint in reference", current);
			reference_name.data--;
			reference_name.length++;
			return false;
		}
		reference_name.data--;
		reference_name.length++;
		return append_to_token(text, (wint_t) codepoint, shift);
	} else {
		read_error("Unrecognized reference name", current);
		return false;
	}
}

struct xml_element {
	string name;
	array_map<string, string> attributes;
	bool preserve_space;

	static inline void free(xml_element& element) {
		core::free(element.name);
		for (auto entry : element.attributes) {
			core::free(entry.key);
			core::free(entry.value);
		}
		core::free(element.attributes);
	}
};

inline bool init(xml_element& element, const array<char>& name, bool preserve_space) {
	if (!init(element.name, name.data, name.length)) {
		return false;
	} else if (!array_map_init(element.attributes, 2)) {
		free(element.name);
		return false;
	}
	element.preserve_space = preserve_space;
	return true;
}

template<bool Validate = false, typename Stream, typename Reader>
bool xml_parse(Stream& input, Reader& reader) {
	position start(1, 1);
	position current = start;
	xml_state state = xml_state::DEFAULT;
	array<char> token = array<char>(1024);
	array<char> reference_token = array<char>(16);
	array<char> attr_token = array<char>(16);

	xml_state prev_state;
	array<xml_element> element_stack(16);
	bool preserve_space = false;
	std::mbstate_t shift = {0};
	std::mbstate_t reference_shift = {0};
	std::mbstate_t attr_shift = {0};
	wint_t next = fgetwc(input);
	bool new_line = false;
	while (next != WEOF) {
		switch (state) {
		case xml_state::DEFAULT:
			if (next == '<') {
				if (token.length != 0 && (preserve_space || !is_whitespace_only(token)) && !emit_text(token, start, current, reader)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				wint_t peek = fgetwc(input);
				start = current;
				if (peek == '/') {
					/* found an end tag */
					state = xml_state::END_TAG_NAME;
					current.column++;
				} else {
					/* found a start tag */
					state = xml_state::START_TAG_NAME;
					ungetwc(peek, input);
				}
				token.clear(); shift = {0};
			} else if (next == '&') {
				prev_state = state;
				state = xml_state::REFERENCE;
			} else if (next == '\r') {
				wint_t peek = fgetwc(input);
				if (peek != '\n')
					ungetwc(peek, input);
				new_line = true;
				if (!append_to_token(token, '\n', shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			} else if (next == '\n') {
				new_line = true;
				if (!append_to_token(token, '\n', shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			} else {
				if (!append_to_token(token, next, shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			}
			break;

		case xml_state::REFERENCE:
			if (next == ';') {
				if (!expand_reference(reference_token, token, shift, current)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				state = prev_state;
				reference_token.clear(); reference_shift = {0};
			} else {
				if (Validate && ((reference_token.length == 0 && !is_valid_name_start(next)) || (reference_token.length != 0 && !is_valid_name_char(next)))) {
					read_error("Invalid character in reference name", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(reference_token, next, reference_shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			}
			break;

		case xml_state::START_TAG_NAME:
			if (is_whitespace(next)) {
				if (!element_stack.ensure_capacity(element_stack.length + 1)
				|| !init(element_stack[element_stack.length], token, preserve_space)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				element_stack.length++;
				state = xml_state::START_TAG;
				token.clear(); shift = {0};
			} else if (next == '>') {
				if (!element_stack.ensure_capacity(element_stack.length + 1)
				|| !init(element_stack[element_stack.length], token, preserve_space)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				element_stack.length++;
				if (!emit_start_element(element_stack.last(), start, current, reader)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				start = current;
				state = xml_state::DEFAULT;
				token.clear(); shift = {0};
			} else {
				if (Validate && ((token.length == 0 && !is_valid_name_start(next)) || (token.length != 0 && !is_valid_name_char(next)))) {
					read_error("Invalid character in tag name", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(token, next, shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			}
			break;

		case xml_state::END_TAG_NAME:
			if (is_whitespace(next)) {
				/* make sure the end tag name matches the name of the corresponding start tag */
				if (element_stack.length == 0) {
					read_error("End tag with no corresponding begin tag", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!compare_strings(element_stack.last().name, token.data, token.length)) {
					fprintf(stderr, "ERROR at %d:%d: End tag with unexpected name. Expected '", current.line, current.column);
					print(element_stack.last().name, stderr); print("'.\n", stderr);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				}
				state = xml_state::END_TAG;
			} else if (next == '>') {
				if (!emit_end_element(element_stack.last(), start, current, reader)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				free(element_stack.last());
				element_stack.length--;
				if (element_stack.length != 0)
					preserve_space = element_stack.last().preserve_space;
				start = current;
				state = xml_state::DEFAULT;
				token.clear(); shift = {0};
			} else {
				if (Validate && ((token.length == 0 && !is_valid_name_start(next)) || (token.length != 0 && !is_valid_name_char(next)))) {
					read_error("Invalid character in tag name", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(token, next, shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			}
			break;

		case xml_state::START_TAG:
			if (is_whitespace(next)) {
				break;
			} else if (next == '/') {
				wint_t peek = fgetwc(input);
				if (peek != '>') {
					/* this is a tag attribute */
					ungetwc(peek, input);
					attr_token.clear(); attr_shift = {0};
					if (Validate && !is_valid_name_start(next)) {
						read_error("Invalid character in attribute name", current);
						for (xml_element& element : element_stack) { free(element); }
						return false;
					} else if (!append_to_token(attr_token, next, attr_shift)) {
						for (xml_element& element : element_stack) { free(element); } return false;
					}
					state = xml_state::ATTRIBUTE_NAME;
				} else {
					/* this is an empty tag */
					current.column++;
					if (!emit_start_element(element_stack.last(), start, current, reader) || !emit_end_element(element_stack.last(), start, current, reader)) {
						for (xml_element& element : element_stack) { free(element); } return false;
					}
					free(element_stack.last());
					element_stack.length--;
					if (element_stack.length != 0)
						preserve_space = element_stack.last().preserve_space;
					start = current;
					state = xml_state::DEFAULT;
					token.clear(); shift = {0};
				}
			} else if (next == '>') {
				if (!emit_start_element(element_stack.last(), start, current, reader)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				start = current;
				state = xml_state::DEFAULT;
				token.clear(); shift = {0};
			} else {
				/* this is a tag attribute */
				attr_token.clear(); attr_shift = {0};
				if (Validate && !is_valid_name_start(next)) {
					read_error("Invalid character in attribute name", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(attr_token, next, attr_shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				state = xml_state::ATTRIBUTE_NAME;
			}
			break;

		case xml_state::END_TAG:
			if (is_whitespace(next)) {
				break;
			} else if (next == '>') {
				if (!emit_end_element(element_stack.last(), start, current, reader)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
				free(element_stack.last());
				element_stack.length--;
				if (element_stack.length != 0)
					preserve_space = element_stack.last().preserve_space;
				start = current;
				state = xml_state::DEFAULT;
				token.clear(); shift = {0};
			} else {
				read_error("Invalid character in end tag", current);
				for (xml_element& element : element_stack) { free(element); }
				return false;
			}
			break;

		case xml_state::ATTRIBUTE_NAME:
			if (next == '=') {
				wint_t peek = fgetwc(input);
				current.column++;
				if (peek == '\'') {
					state = xml_state::ATTRIBUTE_VALUE_SINGLE;
				} else if (peek == '"') {
					state = xml_state::ATTRIBUTE_VALUE_DOUBLE;
				} else {
					read_error("Expected a single or double quote before attribute value", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				}
			} else {
				if (Validate && !is_valid_name_char(next)) {
					read_error("Invalid character in attribute name", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(attr_token, next, attr_shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			}
			break;

		case xml_state::ATTRIBUTE_VALUE_SINGLE:
		case xml_state::ATTRIBUTE_VALUE_DOUBLE:
			if ((state == xml_state::ATTRIBUTE_VALUE_SINGLE && next == '\'')
			 || (state == xml_state::ATTRIBUTE_VALUE_DOUBLE && next == '"'))
			{
				/* add a new attribute */
				if (Validate) {
					/* check that the attribute name is unique */
					for (const auto& entry : element_stack.last().attributes) {
						if (compare_strings(entry.key, attr_token.data, attr_token.length)) {
							read_error("Duplicate attribute name in this tag", current);
							for (xml_element& element : element_stack) { free(element); }
							return false;
						}
					}
				}
				if (!element_stack.last().attributes.ensure_capacity(element_stack.last().attributes.size + 1)
				 || !init(element_stack.last().attributes.keys[element_stack.last().attributes.size], attr_token.data, attr_token.length))
				{
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!init(element_stack.last().attributes.values[element_stack.last().attributes.size], token.data, token.length)) {
					free(element_stack.last().attributes.keys[element_stack.last().attributes.size]);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				}
				element_stack.last().attributes.size++;
				state = xml_state::START_TAG;

				if (compare_strings(attr_token, "xml:space")) {
					if (compare_strings(token, "preserve")) {
						element_stack.last().preserve_space = true;
						preserve_space = true;
					} else if (compare_strings(token, "default")) {
						element_stack.last().preserve_space = false;
						preserve_space = false;
					} else if (Validate) {
						read_error("Illegal value for attribute 'xml:space'", current);
						for (xml_element& element : element_stack) { free(element); }
						return false;
					}
				}
			} else if (next == '&') {
				prev_state = state;
				state = xml_state::REFERENCE;
			} else if (next == '\r') {
				wint_t peek = fgetwc(input);
				if (peek != '\n')
					ungetwc(peek, input);
				new_line = true;
				if (!append_to_token(token, '\n', shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			} else if (next == '\n') {
				new_line = true;
				if (!append_to_token(token, '\n', shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
				}
			} else {
				if (Validate && next == '<') {
					read_error("Invalid character in attribute value", current);
					for (xml_element& element : element_stack) { free(element); }
					return false;
				} else if (!append_to_token(token, next, shift)) {
					for (xml_element& element : element_stack) { free(element); } return false;
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

	for (xml_element& element : element_stack) { free(element); }
	if (!feof(input)) {
		perror("Error occurred while reading XML file");
		return false;
	} else if (Validate && state != xml_state::DEFAULT) {
		read_error("Unexpected end of input", current);
		return false;
	} else if (Validate && element_stack.length != 0) {
		read_error("Some tags were not closed", current);
		return false;
	}
	return true;
}

#endif /* XML_H_ */
