#include "higher_order_logic.h"
#include "set_reasoning.h"

hol_term* canonicalize(hol_term* src) {
	standard_canonicalizer<true, false> canonicalizer;
	hol_term* canonicalized = canonicalize(*src, canonicalizer);
	free(*src); if (src->reference_count == 0) free(src);
	if (canonicalized == NULL) return NULL;
	return canonicalized;
}

int main(int argc, const char** argv)
{
	set_reasoning<hol_term> sets;
	constexpr unsigned int RED = 1, BLUE = 2, FLUFFY = 3, CAT = 4;
	hol_term* cats = hol_term::new_atom(CAT, hol_term::new_variable(1));
	hol_term* red = hol_term::new_atom(RED, hol_term::new_variable(1));
	hol_term* blue = hol_term::new_atom(BLUE, hol_term::new_variable(1));
	hol_term* fluffy = hol_term::new_atom(FLUFFY, hol_term::new_variable(1));
	hol_term* all = hol_term::new_true();
	red->reference_count += 3;
	blue->reference_count += 3;
	cats->reference_count += 5;
	fluffy->reference_count += 1;
	hol_term* red_cats = canonicalize(hol_term::new_and(red, cats));
	hol_term* blue_cats = canonicalize(hol_term::new_and(blue, cats));
	hol_term* red_blue_cats = canonicalize(hol_term::new_and(red, blue, cats));
	hol_term* red_blue_fluffy_cats = canonicalize(hol_term::new_and(red, blue, fluffy, cats));
	hol_term* no_cats = hol_term::new_not(cats);

	hash_map<string, unsigned int> names(256);
	names.put("cat", CAT);
	names.put("red", RED);
	names.put("blue", BLUE);
	names.put("fluffy", FLUFFY);

	const string** reverse_name_map = invert(names);
	string_map_scribe printer = { reverse_name_map, names.table.size + 1 };

	sets.set_size(cats, 6);
	sets.set_size(no_cats, 9);
	sets.set_size(red, 10);
	sets.set_size(blue, 7);
	sets.set_size(red_cats, 4);
	sets.set_size(blue_cats, 2);
	sets.set_size(red_blue_cats, 0);
	sets.set_size(red_blue_fluffy_cats, 0);
	sets.set_size(all, 15);

	unsigned int* clique = NULL; unsigned int clique_count; unsigned int ancestor_of_clique;
	find_largest_disjoint_clique_with_set<INT_MIN>(sets, sets.set_ids.get(*blue_cats), clique, clique_count, ancestor_of_clique);
	print("Largest pairwise-disjoint family of subsets: {", stdout);
	if (clique_count > 0) {
		print(*sets.sets[clique[0]].set_formula, stdout, printer);
		for (unsigned int i = 1; i < clique_count; i++) {
			print(", ", stdout);
			print(*sets.sets[clique[i]].set_formula, stdout, printer);
		}
	}
	print("}.\n", stdout);
	if (clique != NULL) {
		print("  with ancestor ", stdout); print(*sets.sets[ancestor_of_clique].set_formula, stdout, printer); print('\n', stdout);
	}
	free(clique);

	free(*red_blue_cats); free(*red_blue_fluffy_cats);
	free(*red_cats); free(*blue_cats);
	free(*cats); free(*red); free(*blue); free(*fluffy);
	free(*no_cats); free(*all);
	if (cats->reference_count == 0) free(cats);
	if (red->reference_count == 0) free(red);
	if (blue->reference_count == 0) free(blue);
	if (fluffy->reference_count == 0) free(fluffy);
	if (red_cats->reference_count == 0) free(red_cats);
	if (blue_cats->reference_count == 0) free(blue_cats);
	if (red_blue_cats->reference_count == 0) free(red_blue_cats);
	if (red_blue_fluffy_cats->reference_count == 0) free(red_blue_fluffy_cats);
	if (no_cats->reference_count == 0) free(no_cats);
	if (all->reference_count == 0) free(all);
	free(reverse_name_map);
	for (auto entry : names) free(entry.key);
}
