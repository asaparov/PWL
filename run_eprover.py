import sys
import os
import json
import subprocess
from fol import *

def first_orderize(formula):
	def func(f):
		if type(f) == FOLFuncApplication and type(f.function) == FOLVariable:
			return FOLFuncApplication('element_of', f.args + [f.function])
		elif type(f) == FOLEquals and type(f.right) == FOLForAll:
			variables = []
			operand = f.right
			while type(operand) == FOLForAll:
				variables.append(operand.variable)
				operand = operand.operand
			new_operand = FOLIff(FOLFuncApplication('element_of', [FOLVariable(v) for v in variables] + [f.left]), operand.apply(func))
			for variable in reversed(variables):
				new_operand = FOLForAll(variable, new_operand)
			return new_operand
		else:
			return f.apply(func)

	return formula.apply(func)



if len(sys.argv) < 3:
	print("Missing arguments")
	print("Usage: python run_eprover.py [path to E theorem prover directory] [output file]")
	sys.exit(1)

test_file = open('fictionalgeoqa_parse_outputs.txt', 'r')
output_file = open(sys.argv[2], 'w')
line_number = 1
for line in test_file:
	#if line_number < 446 or line_number > 446:
	#	line_number += 1
	#	continue
	if line.startswith('<no answer>'):
		output_file.write('\n')
		line_number += 1
		continue
	example = json.loads('[' + line + ']')[0]
	print(f"LINE NUMBER IS {line_number}")
	for sentence in example:
		for logical_form in sentence:
			print(logical_form)
			if logical_form[0] == '^':
				logical_form[0] = '?'
			formula = do_parse_fol_from_tptp(logical_form.replace('^','!').replace(')(',','))
			print(fol_to_tptp(first_orderize(formula)))
	line_number += 1

output_file.close()
