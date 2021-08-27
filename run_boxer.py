import sys
import os
import json
import subprocess
from fol import *

def normalize_names_and_numbers(formula, numbers):
	def func(f):
		nonlocal numbers
		if type(f) != FOLFuncApplication:
			return f.apply(func)
		function = f.function

		if len(function) > 3 and function[-3].isupper() and function[-2].isdigit() and function[-1].isdigit():
			function = function[:-3]

		if (function.startswith('orgnam') or function.startswith('geonam') or function.startswith('namnam')):
			if len(f.args) != 1:
				raise Exception("normalize_names_and_numbers ERROR: Name declaration must be unary.")
			return FOLFuncApplication('value', [f.args[0], FOLConstant(function[6:])])
		elif function.startswith('c') and function.endswith('number'):
			if len(f.args) != 1:
				raise Exception("normalize_names_and_numbers ERROR: Number literal declaration must be unary.")
			literal = function[len('c'):-len('number')]
			numbers.append(literal)
			return FOLFuncApplication('value', [f.args[0], FOLConstant('num' + literal)])
		elif function.startswith('n1') and function[2:].isdigit():
			if len(f.args) != 1:
				raise Exception("normalize_names_and_numbers ERROR: Number literal declaration must be unary.")
			numbers.append(function[2:])
			return FOLFuncApplication('value', [f.args[0], FOLConstant('num' + function[2:])])
		else:
			if function == 'n1cities':
				function = 'n1city'
			elif function.startswith('n1') and function.endswith('s'):
				function = function[:-1]
			return FOLFuncApplication(function, f.args)

	return formula.apply(func)

def remove_lambda_variable(formula):
	new_variable = max_variable(formula) + 1

	counter = 0
	remove_quantifier = 0
	def remove_lambda_var_apply(f):
		nonlocal counter
		nonlocal remove_quantifier
		if (type(f) == FOLNot and type(f.operand) == FOLExists and type(f.operand.operand) == FOLAnd and len(f.operand.operand.operands) == 2
		 and type(f.operand.operand.operands[0]) == FOLFuncApplication and type(f.operand.operand.operands[1]) == FOLNot
		 and f.operand.operand.operands[0].function == "n12thing" and len(f.operand.operand.operands[0].args) == 1
		 and type(f.operand.operand.operands[0].args[0]) == FOLVariable
		 and f.operand.operand.operands[0].args[0].variable == f.operand.variable):
			counter += 1
			return substitute(f.operand.operand.operands[1].operand, FOLVariable(f.operand.variable), FOLVariable(new_variable))
		elif type(f) == FOLAnd:
			for i in range(len(f.operands)):
				if type(f.operands[i]) == FOLFuncApplication and (f.operands[i].function == "n1what" or f.operands[i].function == "r1what") and len(f.operands[i].args) == 1 and type(f.operands[i].args[0]) == FOLVariable:
					counter += 1
					remove_quantifier = f.operands[i].args[0].variable
					new_operands = [x for j,x in enumerate(f.operands) if j != i]
					if len(new_operands) == 1:
						return new_operands[0]
					else:
						return FOLAnd(new_operands)

		result = f.apply(remove_lambda_var_apply)
		if type(result) == FOLExists and result.variable == remove_quantifier:
			remove_quantifier = 0
			return substitute(result.operand, FOLVariable(result.variable), FOLVariable(new_variable))
		else:
			return result

	new_formula = formula.apply(remove_lambda_var_apply)
	if type(new_formula) == FOLExists and new_formula.variable == remove_quantifier:
		remove_quantifier = 0
		new_formula = substitute(new_formula.operand, FOLVariable(new_formula.variable), FOLVariable(new_variable))

	if counter == 0:
		raise Exception("remove_lambda_variable ERROR: Unable to find lambda variable.")
	elif counter != 1:
		raise Exception("remove_lambda_variable ERROR: Found multiple candidates for lambda variable.")
	else:
		return FOLExists(new_variable + 1,
			FOLExists(new_variable,
				FOLAnd([
					FOLOr([
						FOLFuncApplication('value', [FOLVariable(new_variable), FOLVariable(new_variable + 1)]),
						FOLExists(new_variable + 2, FOLAnd([
							FOLFuncApplication('card', [FOLVariable(new_variable), FOLVariable(new_variable + 2)]),
							FOLFuncApplication('value', [FOLVariable(new_variable + 2), FOLVariable(new_variable + 1)])
						]))
					]),
					new_formula
				])
			)
		)



if len(sys.argv) < 4:
	print("Missing arguments")
	print("Usage: python run_boxer.py [path to C&C/Boxer directory] [path to E theorem prover directory] [output file]")
	sys.exit(1)

test_file = open('fictionalgeoqa.jsonl', 'r')
output_file = open(sys.argv[3], 'w')
line_number = 1
for line in test_file:
	#if line_number < 446 or line_number > 446:
	#	line_number += 1
	#	continue
	example = json.loads(line)
	input_file = open(sys.argv[1] + 'boxer_sentences.txt', 'w')
	input_file.write(example["theory"])

	for question in example["questions"].values():
		input_file.write(' ')
		input_file.write(question['question'].replace('Which','What').replace('which','what'))
	input_file.close()

	subprocess.call(['bin/t', 'a', '--input', 'boxer_sentences.txt', '--output', 'boxer_sentences.tok'], cwd=sys.argv[1])
	subprocess.call(['bin/candc', '--input', 'boxer_sentences.tok', '--output', 'boxer_sentences.ccg', '--models', 'models/boxer', '--candc-printer', 'boxer'], cwd=sys.argv[1])
	subprocess.call(['bin/boxer', 'a', '--input', 'boxer_sentences.ccg', '--output', 'boxer_sentences.out', '--resolve', '--semantics', 'fol'], cwd=sys.argv[1])

	# parse the first-order logic output from Boxer and convert it into TPTP notation
	fol_file = open(sys.argv[1] + 'boxer_sentences.out', 'r')
	formulas = []
	numbers = []
	while True:
		line = fol_file.readline().strip()
		if line == '':
			break
		if not line.startswith('fol('):
			continue
		index = line.find(',', len('fol('))
		formulas.append(normalize_names_and_numbers(do_parse_fol_from_prolog(line[(index + 1):-2]), numbers))
	fol_file.close()

	os.remove(sys.argv[1] + 'boxer_sentences.txt')
	os.remove(sys.argv[1] + 'boxer_sentences.tok')
	os.remove(sys.argv[1] + 'boxer_sentences.ccg')
	os.remove(sys.argv[1] + 'boxer_sentences.out')

	# run the E theorem prover to answer the question
	question_count = len(example["questions"].values())
	skip = False
	for i in range(question_count):
		try:
			formulas[-i-1] = remove_lambda_variable(formulas[-i-1])
		except:
			print(f'ERROR at {line_number}: `remove_lambda_variable` failed with formula {fol_to_tptp(formulas[-i-1])}')
			skip = True
			break
	if skip:
		output_file.write('\n')
		output_file.flush()
		line_number += 1
		continue
	for i in range(question_count):
		prover_file = open(sys.argv[2] + 'boxer_sentences.p', 'w')
		prover_file.write('fof(a0,axiom,![X1]:![L2]:(?[L1]:(n1length(L1) & r1of(L1,X1) & n1kilometer(L1) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2)))) => ?[X4]:(a1longer(X4) & r1than(X4,L2) & r1Theme(X4,X1)))).\n')
		prover_file.write('fof(a1,axiom,![X1]:![L2]:(?[L1]:(n1length(L1) & r1of(L1,X1) & n1kilometer(L1) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1)))) => ?[X4]:(a1shorter(X4) & r1than(X4,L2) & r1Theme(X4,X1)))).\n')
		prover_file.write('fof(a2,axiom,![X1]:![L2]:(?[L1]:(n1length(L1) & r1of(L1,X1) & n1meter(L1) & n1meter(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2)))) => ?[X4]:(a1longer(X4) & r1than(X4,L2) & r1Theme(X4,X1)))).\n')
		prover_file.write('fof(a3,axiom,![X1]:![L2]:(?[L1]:(n1length(L1) & r1of(L1,X1) & n1meter(L1) & n1meter(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1)))) => ?[X4]:(a1shorter(X4) & r1than(X4,L2) & r1Theme(X4,X1)))).\n')
		prover_file.write('fof(a4,axiom,![L1]:![L2]:(?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2))) => ?[X4]:(a1greater(X4) & r1than(X4,L2) & r1Theme(X4,L1)))).\n')
		prover_file.write('fof(a5,axiom,![L1]:![L2]:(?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1))) => ?[X4]:(a1smaller(X4) & r1than(X4,L2) & r1Theme(X4,L1)))).\n')
		prover_file.write('fof(a6,axiom,a1less=a1smaller).\n')
		prover_file.write('fof(a7,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1length(L1) & r1of(L1,X1) & n1kilometer(L1) & n1length(L2) & r1of(L2,X2) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2))))) => ?[M]:(a1longest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a8,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1length(L1) & r1of(L1,X1) & n1kilometer(L1) & n1length(L2) & r1of(L2,X2) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1))))) => ?[M]:(a1shortest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a9,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1length(L1) & r1of(L1,X1) & n1meter(L1) & n1length(L2) & r1of(L2,X2) & n1meter(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2))))) => ?[M]:(a1longest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a10,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1length(L1) & r1of(L1,X1) & n1meter(L1) & n1length(L2) & r1of(L2,X2) & n1meter(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1))))) => ?[M]:(a1shortest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a11,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1population(L1) & r1of(L1,X1) & n1population(L2) & r1of(L2,X2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2))))) => ?[M]:(a1largest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a12,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1population(L1) & r1of(L1,X1) & n1population(L2) & r1of(L2,X2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1))))) => ?[M]:(a1smallest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a13,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1area(L1) & r1of(L1,X1) & n1kilometer(L1) & n1area(L2) & r1of(L2,X2) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V1,V2))))) => ?[M]:(a1largest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a14,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(n1area(L1) & r1of(L1,X1) & n1kilometer(L1) & n1area(L2) & r1of(L2,X2) & n1kilometer(L2) & ?[C1,C2]:(card(L1,C1) & card(L2,C2) & ?[V1,V2]:(value(C1,V1) & value(C2,V2) & greater(V2,V1))))) => ?[M]:(a1smallest(M) & r1Theme(M,X1)))).\n')
		prover_file.write('fof(a15,axiom,a1largest=a1biggest).\n')
		prover_file.write('fof(a16,axiom,![X]:![Y]:(r1in(Y,X) => ?[H]:(v1have(H) & r1Actor(H,X) & r1Theme(H,Y)))).\n')
		prover_file.write('fof(a17,axiom,v1run=v1flow).\n')
		axiom_count = 0
		for a in numbers:
			for b in numbers:
				if float(a) > float(b):
					prover_file.write(f'fof(g{axiom_count},axiom,greater(n{a},n{b})).\n')
					axiom_count += 1
				elif a != b:
					prover_file.write(f'fof(g{axiom_count},axiom,~greater(n{a},n{b})).\n')
					axiom_count += 1
		for j in range(len(formulas) - question_count):
			prover_file.write(f'fof(c{j},axiom,({fol_to_tptp(formulas[j])})).\n')
		prover_file.write(f'fof(q,question,({fol_to_tptp(formulas[-i-1])})).\n')
		prover_file.close()

		result = subprocess.run(['PROVER/eprover', '-s', 'boxer_sentences.p', '--soft-cpu-limit=45'], cwd=sys.argv[2], stdout=subprocess.PIPE, universal_newlines=True)
		answers = []
		for line in result.stdout.split('\n'):
			if line.startswith('# SZS answers Tuple [['):
				index = line.find(']')
				candidates = line[len('# SZS answers Tuple [['):index].split(',')
				for candidate in candidates:
					candidate = candidate.strip()
					if candidate.startswith('esk'):
						continue
					if candidate.startswith('num'):
						answers.append(candidate[3:])
					else:
						answers.append(candidate)
		output_file.write(', '.join(answers) + '\n')
		output_file.flush()
		#os.remove(sys.argv[2] + 'boxer_sentences.p')

	line_number += 1

output_file.close()
