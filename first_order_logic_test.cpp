#include "first_order_logic.h"

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

int main(int argc, const char** argv)
{
	const char* formulas_src = "has_color(frog,green)\n"
		"![X]:(frog(X)=>green(X))\n"
		"\t![X]: (  frog( X)  =>\t green (  X\t) )\n"
		"![X,Y]:((cat(X) & dog(Y))=> ~like(X,Y)) | ?[Z]:((((~frog(Z)))))\n";
	memory_stream in = memory_stream(formulas_src, strlen(formulas_src));

	array<fol_formula*> formulas = array<fol_formula*>(16);
	hash_map<string, unsigned int> names = hash_map<string, unsigned int>(1024);
	if (!read_formulas(formulas, in, names)) {
		cleanup(names, formulas); return EXIT_FAILURE;
	}

	const string** name_ids = invert(names);
	string_map_scribe printer = { name_ids, names.table.size + 1 };
	for (const fol_formula* formula : formulas) {
		print(*formula, stderr, printer); print('\n', stderr);
	}
	cleanup(names, formulas); free(name_ids);
	return EXIT_SUCCESS;
}
