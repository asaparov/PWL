#ifndef FAKE_PARSER_H_
#define FAKE_PARSER_H_

#include "article.h"
#include "first_order_logic.h"

struct fake_parser {
	hash_map<sentence, fol_formula*> table;

	fake_parser() : table(64) { }
};

#endif /* FAKE_PARSER_H_ */
