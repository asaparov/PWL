#ifndef FAKE_PARSER_H_
#define FAKE_PARSER_H_

#include "article.h"

template<typename Formula, typename Printer>
struct fake_parser_printer {
	Printer& constant_printer;
	const hash_map<unsigned int, unsigned int>& reverse_map;

	fake_parser_printer(Printer& constant_printer, const hash_map<unsigned int, unsigned int>& reverse_symbol_map) :
			constant_printer(constant_printer), reverse_map(reverse_symbol_map) { }
};

template<typename Formula>
struct fake_parser {
	hash_map<sentence, sentence_label<Formula>> table;
	hash_map<sentence_token, unsigned int> learned_tokens;
	hash_map<unsigned int, unsigned int> symbol_map;
	hash_map<unsigned int, unsigned int> reverse_symbol_map;
	unsigned int unknown_id;

	fake_parser(unsigned int unknown_id) : table(64), learned_tokens(32), symbol_map(32), reverse_symbol_map(32), unknown_id(unknown_id) { }

	~fake_parser() {
		for (auto entry : table) {
			free(entry.key);
			free(entry.value);
		}
	}

	template<unsigned int K, typename TheoryType>
	bool parse(const sentence& s, Formula** logical_forms,
			double* log_probabilities, unsigned int& parse_count,
			const TheoryType& T, array<sentence_token>& unrecognized) const
	{
		static_assert(K > 0, "`K` must be at least 1.");

		bool contains;
		const sentence_label<Formula>& parse = table.get(s, contains);
		if (!contains) {
			parse_count = 0;
			return true;
		}

		array_map<unsigned int, unsigned int> constant_map(parse.labels.size);
		for (const auto& entry : parse.labels) {
			unsigned int constant = learned_tokens.get({entry.key}, contains);
			if (contains) {
				constant_map.put({entry.value}, constant);
			} else {
				if (!unrecognized.add({entry.key})) return false;
				constant_map.put({entry.value}, unknown_id);
			}
		}

		Formula* out = relabel_constants(parse.logical_form, constant_map);
		if (out == NULL) {
			fprintf(stderr, "fake_parser.parse ERROR: Failed to relabel unknown constants.\n");
			return false;
		}

		logical_forms[0] = out;
		log_probabilities[0] = 0.0;
		parse_count = 1;
		return true;
	}

	bool add_definition(const sentence& s, const Formula* definition, unsigned int new_constant)
	{
		bool contains;
		sentence_label<Formula>& label = table.get(s, contains);
		if (!contains || !learned_tokens.check_size(learned_tokens.table.size + label.labels.size))
			return false;

		unsigned int bucket;
		for (const auto& entry : label.labels) {
			unsigned int& token = learned_tokens.get({entry.key}, contains, bucket);
			if (!contains) {
				learned_tokens.table.keys[bucket] = {entry.key};
				learned_tokens.table.size++;
				token = new_constant;
			}
		}

		if (label.labels.size > 0)
			return symbol_map.put(label.labels.values[0], new_constant) && reverse_symbol_map.put(new_constant, label.labels.values[0]);
		return true;
	}

	template<typename Printer>
	inline fake_parser_printer<Formula, Printer> get_printer(Printer& constant_printer) const {
		return fake_parser_printer<Formula, Printer>(constant_printer, reverse_symbol_map);
	}
};

template<typename Stream, typename Formula, typename Printer>
inline bool print(unsigned int constant, Stream& out, const fake_parser_printer<Formula, Printer>& printer) {
	bool contains;
	unsigned int value = printer.reverse_map.get(constant, contains);
	if (contains) constant = value;
	return print(constant, out, printer.constant_printer);
}

#endif /* FAKE_PARSER_H_ */
