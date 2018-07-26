#include "first_order_logic.h"

#include <locale.h>

using namespace core;

template<typename Stream>
bool read_formulas(
	array<fol_formula*>& formulas, Stream& in,
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
		fol_formula* formula = (fol_formula*) malloc(sizeof(fol_formula));
		if (formula == NULL) {
			fprintf(stderr, "read_formulas ERROR: Out of memory.\n");
			free_tokens(tokens); return false;
		} else if (!tptp_interpret(tokens, index, *formula, names, variables)) {
			fprintf(stderr, "ERROR: Unable to parse first-order formula.\n");
			free(formula); free_tokens(tokens); return false;
		} else if (!formulas.add(formula)) {
			free(*formula); free(formula);
			free_tokens(tokens); return false;
		}

		if (variables.size != 0)
			fprintf(stderr, "WARNING: Variable map is not empty.\n");
	}
	free_tokens(tokens);
	return true;
}

void cleanup(
	hash_map<string, unsigned int>& names,
	array<fol_formula*>& formulas)
{
	for (auto entry : names) free(entry.key);
	for (fol_formula* formula : formulas) {
		free(*formula); free(formula);
	}
}

bool fol_test_simple(const char* filename = "logical_forms.txt")
{
	setlocale(LC_CTYPE, "en_US.UTF-8");
	FILE* in = open_file(filename, "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	array<fol_formula*> formulas = array<fol_formula*>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_formulas(formulas, in, names)) {
		fclose(in); cleanup(names, formulas); return false;
	}
	fclose(in);

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (const fol_formula* formula : formulas) {
		print(*formula, stderr, printer); print('\n', stderr);
	}
	cleanup(names, formulas); free(name_ids);
	return true;
}

bool fol_test_canonicalization(const char* filename = "canonicalization_test.txt")
{
	setlocale(LC_CTYPE, "en_US.UTF-8");
	FILE* in = open_file(filename, "r");
	if (in == NULL) {
		fprintf(stderr, "ERROR: Unable to open '%s' for reading.\n", filename);
		return false;
	}

	array<fol_formula*> formulas = array<fol_formula*>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_formulas(formulas, in, names)) {
		fclose(in); cleanup(names, formulas); return false;
	}
	fclose(in);

	if (formulas.length % 2 != 0) {
		fprintf(stderr, "ERROR: The canonicalization test requires an even number of input formulas.\n");
		cleanup(names, formulas); return false;
	}

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (unsigned int i = 0; i < formulas.length; i += 2) {
		fol_formula* canonicalized = canonicalize(*formulas[i]);
		if (canonicalized == NULL) {
			fprintf(stderr, "ERROR: Unable to canonicalize example %u.\n", i / 2);
			continue;
		} else if (*canonicalized != *formulas[i + 1]) {
			fprintf(stderr, "ERROR: The canonicalized form of example %u does not match the expected output.\n", i / 2);
			print("  Original form: ", stderr); print(*formulas[i], stderr, printer); print('\n', stderr);
			print("  Canonicalized: ", stderr); print(*canonicalized, stderr, printer); print('\n', stderr);
			print("  Expected form: ", stderr); print(*formulas[i + 1], stderr, printer); print('\n', stderr);
		}
		free(*canonicalized);
		if (canonicalized->reference_count == 0)
			free(canonicalized);
	}
	cleanup(names, formulas); free(name_ids);
	return true;
}

int main(int argc, const char** argv)
{
	fol_test_canonicalization();
}
