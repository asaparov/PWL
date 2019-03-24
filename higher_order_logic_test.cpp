#include "higher_order_logic.h"

#include <locale.h>

using namespace core;

template<typename Stream>
bool read_terms(
	array<hol_term*>& terms, Stream& in,
	hash_map<string, unsigned int>& names)
{
	array<tptp_token> tokens = array<tptp_token>(512);
	if (!tptp_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis failed.\n");
		free_tokens(tokens); return false;
	}

	unsigned int index = 0;
	while (index < tokens.length) {
		array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
		hol_term* term = (hol_term*) malloc(sizeof(hol_term));
		if (term == NULL) {
			fprintf(stderr, "read_terms ERROR: Out of memory.\n");
			free_tokens(tokens); return false;
		} else if (!tptp_interpret(tokens, index, *term, names, variables)) {
			fprintf(stderr, "ERROR: Unable to parse higher-order term.\n");
			for (auto entry : variables) free(entry.key);
			free(term); free_tokens(tokens); return false;
		} else if (!expect_token(tokens, index, tptp_token_type::SEMICOLON, "semicolon at end of higher-order term") || !terms.add(term)) {
			free(*term); free(term);
			free_tokens(tokens); return false;
		}
		index++;

		if (variables.size != 0)
			fprintf(stderr, "WARNING: Variable map is not empty.\n");
	}
	free_tokens(tokens);
	return true;
}

struct type_statements {
	array_map<hol_term, hol_type> types;

	static inline void free(type_statements& example) {
		for (auto entry : example.types) {
			core::free(entry.key);
			core::free(entry.value);
		}
		core::free(example.types);
	}
};

inline bool init(type_statements& t) {
	return array_map_init(t.types, 8);
}

template<typename Stream>
bool read_types(
	array<type_statements>& examples, Stream& in,
	hash_map<string, unsigned int>& names)
{
	array<tptp_token> tokens = array<tptp_token>(512);
	if (!tptp_lex(tokens, in)) {
		fprintf(stderr, "ERROR: Lexical analysis failed.\n");
		free_tokens(tokens); return false;
	}

	unsigned int index = 0;
	while (index < tokens.length) {
		if (!examples.ensure_capacity(examples.length + 1)
		 || !init(examples[examples.length]))
		{
			free_tokens(tokens);
			return false;
		}
		type_statements& example = examples[examples.length];
		examples.length++;

		while (true) {
			array_map<string, unsigned int> variables = array_map<string, unsigned int>(16);
			if (!example.types.ensure_capacity(example.types.size + 1)) {
				free_tokens(tokens);
				return false;
			}

			hol_term& term = example.types.keys[example.types.size];
			hol_type& type = example.types.values[example.types.size];
			if (!tptp_interpret(tokens, index, term, type, names, variables)) {
				fprintf(stderr, "ERROR: Unable to parse higher-order type statement.\n");
				for (auto entry : variables) free(entry.key);
				free_tokens(tokens); return false;
			}
			example.types.size++;

			if (variables.size != 0) {
				fprintf(stderr, "WARNING: Variable map is not empty.\n");
				for (auto entry : variables) free(entry.key);
			}

			if (index >= tokens.length) {
				read_error("Unexpected end of input", tokens.last().end);
				free_tokens(tokens);
				return false;
			} else if (tokens[index].type == tptp_token_type::SEMICOLON) {
				index++; break;
			} else if (tokens[index].type != tptp_token_type::COMMA) {
				read_error("Expected a comma separating type statements", tokens[index].start);
				free_tokens(tokens);
				return false;
			}
			index++;
		}
	}
	free_tokens(tokens);
	return true;
}

void cleanup(
	hash_map<string, unsigned int>& names,
	array<hol_term*>& terms)
{
	for (auto entry : names) free(entry.key);
	for (hol_term* term : terms) {
		free(*term); free(term);
	}
}

void cleanup(
	hash_map<string, unsigned int>& names,
	array<type_statements>& examples)
{
	for (auto entry : names) free(entry.key);
	free_elements(examples);
}

bool hol_test_simple(const char* filename = "hol_logical_forms.txt")
{
	setlocale(LC_CTYPE, "en_US.UTF-8");
	FILE* in = open_file(filename, "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	array<hol_term*> terms = array<hol_term*>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_terms(terms, in, names)) {
		fclose(in); cleanup(names, terms); return false;
	}
	fclose(in);

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (const hol_term* term : terms) {
		print(*term, stderr, printer); print('\n', stderr);
	}
	cleanup(names, terms); free(name_ids);
	print("hol_test_simple: Test complete.\n", stdout); fflush(stdout);
	return true;
}

bool hol_test_canonicalization(const char* filename = "hol_canonicalization_test.txt")
{
	setlocale(LC_CTYPE, "en_US.UTF-8");
	FILE* in = open_file(filename, "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	array<hol_term*> terms = array<hol_term*>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_terms(terms, in, names)) {
		fclose(in); cleanup(names, terms); return false;
	}
	fclose(in);

	if (terms.length % 2 != 0) {
		fprintf(stderr, "ERROR: The canonicalization test requires an even number of input terms.\n");
		cleanup(names, terms); return false;
	}

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (unsigned int i = 0; i < terms.length; i += 2) {
		hol_term* canonicalized = canonicalize(*terms[i], standard_canonicalizer<true, false>());
		if (canonicalized == NULL) {
			fprintf(stderr, "ERROR: Unable to canonicalize example %u.\n", i / 2);
			continue;
		} else if (*canonicalized != *terms[i + 1]) {
			fprintf(stderr, "ERROR: The canonicalized form of example %u does not match the expected output.\n", i / 2);
			print("  Original form: ", stderr); print(*terms[i], stderr, printer); print('\n', stderr);
			print("  Canonicalized: ", stderr); print(*canonicalized, stderr, printer); print('\n', stderr);
			print("  Expected form: ", stderr); print(*terms[i + 1], stderr, printer); print('\n', stderr);
		}
		free(*canonicalized);
		if (canonicalized->reference_count == 0)
			free(canonicalized);
	}
	cleanup(names, terms); free(name_ids);
	print("hol_test_canonicalization: Test complete.\n", stdout); fflush(stdout);
	return true;
}

bool hol_test_compute_type(const char* filename = "hol_type_check_test.txt")
{
	setlocale(LC_CTYPE, "en_US.UTF-8");
	FILE* in = open_file(filename, "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	array<type_statements> examples = array<type_statements>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_types(examples, in, names)) {
		fclose(in); cleanup(names, examples); return false;
	}
	fclose(in);

	unsigned int not_well_typed;
	if (!get_token("not_well_typed", not_well_typed, names)) {
		cleanup(names, examples); return false;
	}

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (unsigned int i = 0; i < examples.length; i++)
	{
		bool success;
		const type_statements& example = examples[i];
		if (example.types.size == 0) {
			fprintf(stderr, "ERROR: Test example %u is empty.\n", i);
			continue;
		}
		const hol_term& term = example.types.keys[0];
		const hol_type& term_type = example.types.values[0];

		bool expect_success = true;
		for (const auto& entry : example.types) {
			if (entry.key.type == hol_term_type::CONSTANT && entry.key.constant == not_well_typed) {
				expect_success = false;
				break;
			}
		}

		type_map types(16);
		array_map<unsigned int, hol_type> constant_types(8);
		array_map<unsigned int, hol_type> variable_types(8);
		array_map<unsigned int, hol_type> parameter_types(8);

		/*array<hol_type> type_variables(8);
		if (!init(type_variables[0], hol_type_kind::ANY)) {
			cleanup(names, examples); free(name_ids);
			return false;
		}
		type_variables.length++;

		print(CONSOLE_BOLD "Example ", stdout); print(i, stdout); print(":" CONSOLE_RESET "\n", stdout);
		hol_type type(0);
		success = compute_type<false>(term, types, type, constant_types, variable_types, parameter_types, type_variables);

		if (success) {
			print("Term type: ", stdout); print(type, stdout); print('\n', stdout);

			if (constant_types.size > 0) {
				print("Constant types:\n", stdout);
				for (auto entry : constant_types) {
					print("  ", stdout); print(entry.key, stdout, printer); print(" : ", stdout);
					print(entry.value, stdout); print('\n', stdout);
				}
			} if (type_variables.length > 0) {
				print("where:\n", stdout);
				for (unsigned int i = 0; i < type_variables.length; i++) {
					print("  ", stdout); print_variable(i, stdout); print(" = ", stdout);
					print(type_variables[i], stdout); print("\n", stdout);
				}
			}
		} else {
			fprintf(stderr, "ERROR: On example %u, 'compute_type' returned false.\n", i);
		}

		free_elements(type_variables);
		for (unsigned int j = 0; j < constant_types.size; j++) free(constant_types.values[j]);
		for (unsigned int j = 0; j < variable_types.size; j++) free(variable_types.values[j]);
		for (unsigned int j = 0; j < parameter_types.size; j++) free(parameter_types.values[j]);*/

		types.clear(); constant_types.clear();
		variable_types.clear(); parameter_types.clear();
		if (!expect_success)
			print("Expecting 'compute_type' to output error on this example:\n", stdout);
		success = compute_type<false>(term, types, constant_types, variable_types, parameter_types);
		if (success) {
			if (!expect_success)
				fprintf(stderr, "ERROR: On example %u, expected 'compute_type' to return false, but it returned true.\n", i);
			if (types.types.get(&term) != term_type) {
				fprintf(stderr, "ERROR: On example %u, computed term type differs from expected term type.\n", i);
				print("  Computed type: ", stderr); print(types.types.get(&term), stderr); print('\n', stderr);
				print("  Expected type: ", stderr); print(term_type, stderr); print('\n', stderr);
			} /*else {
				print("Term type: ", stdout); print(type, stdout); print('\n', stdout);	
			}*/

			/*if (constant_types.size > 0) {
				print("Flattened constant types:\n", stdout);
				for (auto entry : constant_types) {
					print("  ", stdout); print(entry.key, stdout, printer); print(" : ", stdout);
					print(entry.value, stdout); print('\n', stdout);
				}
			}*/

			for (unsigned int j = 1; expect_success && j < example.types.size; j++) {
				const hol_term& term = example.types.keys[j];
				const hol_type& expected_type = example.types.values[j];
				if (term.type == hol_term_type::CONSTANT) {
					bool contains;
					const hol_type& computed_type = constant_types.get(term.constant, contains);
					if (!contains) {
						fprintf(stderr, "ERROR: On example %u, 'compute_type' did not compute the type of '", i);
						print(term, stderr, printer); print("'.\n", stderr);
					} else if (computed_type != expected_type) {
						fprintf(stderr, "ERROR: On example %u, the computed type of '", i);
						print(term, stderr, printer); print("' does not match the expected type.\n", stderr);
						print("  Computed type: ", stderr); print(computed_type, stderr); print('\n', stderr);
						print("  Expected type: ", stderr); print(expected_type, stderr); print('\n', stderr);
					}
				}
			}
		} else {
			if (expect_success)
				fprintf(stderr, "ERROR: On example %u, 'compute_type' returned false.\n", i);
		}

		for (unsigned int j = 0; j < constant_types.size; j++) free(constant_types.values[j]);
		for (unsigned int j = 0; j < variable_types.size; j++) free(variable_types.values[j]);
		for (unsigned int j = 0; j < parameter_types.size; j++) free(parameter_types.values[j]);
		//print('\n', stdout);
	}
	cleanup(names, examples); free(name_ids);
	print("hol_test_compute_type: Test complete.\n", stdout); fflush(stdout);
	return true;
}

int main(int argc, const char** argv)
{
	hol_test_compute_type();
}
