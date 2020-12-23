/**
 * ruletaker.h - Code to read RukeTaker data and attempt to answer the questions.
 *
 *  Created on: Dec 8, 2020
 *      Author: asaparov
 */

#ifndef RULETAKER_H_
#define RULETAKER_H_

#include <stdio.h>

#include "json.h"

enum class ruletaker_key {
	CONTEXT,
	QUESTIONS,
	TEXT,
	LABEL,
	OTHER
};

enum class ruletaker_reader_state {
	START,
	ROOT_OBJECT,
	QUESTIONS,
	QUESTION
};

struct ruletaker_reader {
	ruletaker_reader_state state;
	ruletaker_key current_key;
	char* context;
	pair<char*, bool> next_question;
	array<pair<string, bool>> questions;

	ruletaker_reader() : state(ruletaker_reader_state::START), context(nullptr), next_question(nullptr, false), questions(44) { }

	~ruletaker_reader() {
		if (context != nullptr)
			free(context);
		if (next_question.key != nullptr)
			free(next_question.key);
		for (pair<string, bool>& question : questions)
			free(question.key);
	}
};

inline bool begin_object(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::START) {
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	} else if (reader.state == ruletaker_reader_state::QUESTIONS) {
		reader.state = ruletaker_reader_state::QUESTION;
	}
	return true;
}

inline bool end_object(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::ROOT_OBJECT) {
		reader.state = ruletaker_reader_state::START;
	} else if (reader.state == ruletaker_reader_state::QUESTION) {
		if (!reader.questions.ensure_capacity(reader.questions.length + 1))
			return false;
		reader.questions[reader.questions.length].key.data = reader.next_question.key;
		reader.questions[reader.questions.length].key.length = strlen(reader.next_question.key);
		reader.questions[reader.questions.length++].value = reader.next_question.value;
		reader.next_question.key = nullptr;
		reader.state = ruletaker_reader_state::QUESTIONS;
	}
	return true;
}

inline bool emit_key(const array<char>& key_name, const position& pos, ruletaker_reader& reader) {
	if (compare_strings(key_name, "context")) {
		reader.current_key = ruletaker_key::CONTEXT;
	} else if (compare_strings(key_name, "questions")) {
		reader.current_key = ruletaker_key::QUESTIONS;
	} else if (compare_strings(key_name, "text")) {
		reader.current_key = ruletaker_key::TEXT;
	} else if (compare_strings(key_name, "label")) {
		reader.current_key = ruletaker_key::LABEL;
	} else {
		reader.current_key = ruletaker_key::OTHER;
	}
	return true;
}

inline bool begin_list(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::QUESTIONS)
		reader.state = ruletaker_reader_state::QUESTIONS;
	return true;
}

inline bool end_list(const position& pos, ruletaker_reader& reader) {
	if (reader.state == ruletaker_reader_state::QUESTIONS)
		reader.state = ruletaker_reader_state::ROOT_OBJECT;
	return true;
}

inline bool emit_string(const array<char>& str, const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::CONTEXT) {
		if (reader.context != nullptr) {
			read_error("Found duplicate `context` entry", pos);
			return false;
		}
		reader.context = (char*) malloc(sizeof(char) * (str.length + 1));
		if (reader.context == nullptr) {
			fprintf(stderr, "emit_string ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str.length; i++)
			reader.context[i] = str[i];
		reader.context[str.length] = '\0';
	} else if (reader.current_key == ruletaker_key::TEXT) {
		if (reader.next_question.key != nullptr) {
			read_error("Found duplicate `text` entry for question", pos);
			return false;
		}
		reader.next_question.key = (char*) malloc(sizeof(char) * (str.length + 1));
		if (reader.next_question.key == nullptr) {
			fprintf(stderr, "emit_string ERROR: Out of memory.\n");
			return false;
		}
		for (unsigned int i = 0; i < str.length; i++)
			reader.next_question.key[i] = str[i];
		reader.next_question.key[str.length] = '\0';
	}
	return true;
}

inline bool emit_true(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::LABEL)
		reader.next_question.value = true;
	return true;
}

inline bool emit_false(const position& pos, ruletaker_reader& reader) {
	if (reader.current_key == ruletaker_key::LABEL)
		reader.next_question.value = false;
	return true;
}

constexpr bool emit_null(const position& pos, ruletaker_reader& reader) { return true; }
constexpr bool emit_number(double value, const position& pos, ruletaker_reader& reader) { return true; }

template<typename Stream>
struct line_reader {
	Stream& in;
	bool eof;

	line_reader(Stream& in) : in(in), eof(false) { }
};

template<typename Stream>
int fgetc(line_reader<Stream>& lr) {
	int c = fgetc(lr.in);
	if (c == -1)
		lr.eof = true;
	else if (c == '\n')
		return -1;
	return c;
}

template<typename ProcessQuestions>
bool read_ruletaker_data(const char* filename, ProcessQuestions process_questions)
{
	FILE* in = (FILE*) fopen(filename, "rb");
	if (in == nullptr) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	line_reader<FILE*> lr(in);
	position current(1, 1);
	while (!lr.eof) {
		ruletaker_reader reader;
		if (!json_parse(lr, reader, current)) {
			fclose(in);
			return false;
		}
		current.line++;
		current.column = 1;

		if (!process_questions(reader.context, reader.questions)) {
			fclose(in);
			return false;
		}
	}

	fclose(in);
	return true;
}

#endif /* RULETAKER_H_ */
