#ifndef FAKE_PARSER_H_
#define FAKE_PARSER_H_

#include "article.h"

struct constant_relabeler {
	const array_map<unsigned int, unsigned int>& map;
};

inline bool clone_constant(unsigned int src_constant, unsigned int& dst_constant, constant_relabeler& relabeler) {
	bool contains;
	const unsigned int& dst = relabeler.map.get(src_constant, contains);
	if (contains) {
		dst_constant = dst;
		return true;
	} else {
		dst_constant = src_constant;
		return true;
	}
}

inline bool clone_predicate(unsigned int src_predicate, unsigned int& dst_predicate, constant_relabeler& relabeler) {
	bool contains;
	const unsigned int& dst = relabeler.map.get(src_predicate, contains);
	if (contains) {
		dst_predicate = dst;
		return true;
	} else {
		dst_predicate = src_predicate;
		return true;
	}
}

inline bool clone_variable(unsigned int src_variable, unsigned int& dst_variable, constant_relabeler& relabeler) {
	return clone_variable(src_variable, dst_variable);
}

inline bool clone_parameter(unsigned int src_parameter, unsigned int& dst_parameter, constant_relabeler& relabeler) {
	return clone_parameter(src_parameter, dst_parameter);
}

template<typename Formula>
inline Formula* relabel_constants(const Formula* src,
		const array_map<unsigned int, unsigned int>& constant_map)
{
	Formula* dst = (Formula*) malloc(sizeof(Formula));
	if (dst == NULL) {
		fprintf(stderr, "relabel_constants ERROR: Out of memory.\n");
		return NULL;
	}

	constant_relabeler relabeler = {constant_map};
	if (!clone(*src, *dst, relabeler)) {
		free(dst); return NULL;
	}
	return dst;
}

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
	hash_map<token, unsigned int> learned_tokens;
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
			const TheoryType& T, array<token>& unrecognized) const
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
