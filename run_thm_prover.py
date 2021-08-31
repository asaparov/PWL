import sys
import os
import json
import subprocess
import threading
from fol import *

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def number_symbol(num):
	return 'num' + str(num).replace('.','_')

def first_orderize(formula, numbers, names):
	def func(f):
		nonlocal numbers
		if type(f) == FOLFuncApplication and type(f.function) == FOLVariable:
			return FOLFuncApplication('element_of', [func(arg) for arg in f.args] + [f.function])
		elif type(f) == FOLFuncApplication and f.function == '≥':
			return FOLFuncApplication('greater', [func(arg) for arg in f.args])
		elif type(f) == FOLEquals and type(f.right) == FOLForAll:
			variables = []
			operand = f.right
			while type(operand) == FOLForAll:
				variables.append(operand.variable)
				operand = operand.operand
			new_operand = FOLIff(FOLFuncApplication('element_of', [FOLVariable(v) for v in variables] + [func(f.left)]), func(operand))
			for variable in reversed(variables):
				new_operand = FOLForAll(variable, new_operand)
			return new_operand
		elif type(f) == FOLConstant and isfloat(f.constant):
			numbers.add(f.constant)
			return FOLConstant(number_symbol(f.constant))
		elif type(f) == FOLConstant and f.constant[0] == '"':
			names.add(f.constant[1:-1])
			return f
		else:
			return f.apply(func)

	return func(formula)

def run_thm_prover(formulas, presuppositions, question_count, question_lf, numbers, use_vampire, candidate):
	if candidate != None:
		print("attempting proof with candidate " + candidate)
	filename = f'input_sentences.{threading.current_thread().ident}.p'
	prover_file = open(sys.argv[1] + filename, 'w')
	prover_file.write(f'fof(a0,axiom,size(empty_set)={number_symbol(0)} & ~?[X]:element_of(X,empty_set) & ![S]:(~?[X]:element_of(X,S) => (S=empty_set))).\n')
	prover_file.write('fof(a1,axiom,![S,X]:(~element_of(X,S) => ?[U]:(![Y]:(element_of(Y,U) <=> (Y=X | element_of(Y,S))) & size(U)=successor(size(S))))).\n')
	if has_greatest or has_least:
		prover_file.write('fof(a2,axiom,![X]:((square(X) & kilometer(X)) => (X=unit_sqkm))).\n')
		prover_file.write('fof(a3,axiom,![X]:((~square(X) & kilometer(X)) => (X=unit_km))).\n')
		prover_file.write('fof(a4,axiom,![X]:(meter(X) => (X=unit_m))).\n')
		prover_file.write('fof(a5,axiom,![X,Y,W,V]:((measure(X) & measure(Y) & ?[U]:(arg2(X)=U & arg2(Y)=U) & arg1(X)=W & arg1(Y)=V & greater(W,V)) => greater(X,Y))).\n')
		prover_file.write('fof(a6,axiom,![X,Y,W,V]:((measure(X) & measure(Y) & ?[U]:(arg2(X)=U & arg2(Y)=U) & arg1(X)=W & arg1(Y)=V & (W=V)) => (X=Y))).\n')
	if has_greatest:
		prover_file.write('fof(a7,axiom,![X,F,S]:(?[G]:(greatest(F,G) & (arg1(G)=S) & (arg2(G)=X)) <=> (element_of(X,S) & ![V]:(element_of(X,V,F) => ![Y]:(element_of(Y,S) => ![W]:(element_of(Y,W,F) => (greater(V,W) | V=W))))))).\n')
	if has_least:
		prover_file.write('fof(a8,axiom,![X,F,S]:(?[G]:(least(F,G)    & (arg1(G)=S) & (arg2(G)=X)) <=> (element_of(X,S) & ![V]:(element_of(X,V,F) => ![Y]:(element_of(Y,S) => ![W]:(element_of(Y,W,F) => (greater(W,V) | V=W))))))).\n')
	successor_count = 10
	for n in range(successor_count):
		prover_file.write(f'fof(s{n},axiom,successor({number_symbol(n)})={number_symbol(n+1)}).\n')
	axiom_count = 0
	for a in numbers:
		for b in numbers:
			if float(a) > float(b):
				prover_file.write(f'fof(g{axiom_count},axiom,greater({number_symbol(a)},{number_symbol(b)})).\n')
				axiom_count += 1
			elif a != b:
				prover_file.write(f'fof(g{axiom_count},axiom,~greater({number_symbol(a)},{number_symbol(b)})).\n')
				axiom_count += 1
				prover_file.write(f'fof(g{axiom_count},axiom,~({number_symbol(a)}={number_symbol(b)})).\n')
				axiom_count += 1
	for j in range(len(formulas) - question_count):
		prover_file.write(f'fof(c{j},axiom,({fol_to_tptp(formulas[j])})).\n')
	for j in range(len(presuppositions)):
		prover_file.write(f'fof(p{j},axiom,({fol_to_tptp(presuppositions[j])})).\n')
	if candidate != None:
		if isfloat(candidate[0]):
			question_lf = substitute(question_lf.operand, FOLVariable(question_lf.variable), FOLConstant(number_symbol(candidate)))
		else:
			if type(question_lf.operand) == FOLAnd:
				new_operands = question_lf.operand.operands
			else:
				new_operands = [question_lf.operand]
			new_var = max_variable(question_lf) + 1
			question_lf = FOLExists(question_lf.variable, FOLAnd([
					FOLExists(new_var, FOLAnd([
						FOLFuncApplication("name", [FOLVariable(new_var)]),
						FOLEquals(FOLFuncApplication("arg1", [FOLVariable(new_var)]), FOLVariable(question_lf.variable)),
						FOLEquals(FOLFuncApplication("arg2", [FOLVariable(new_var)]), FOLConstant(f'"{candidate}"'))
					]))
				] + new_operands))
	prover_file.write(f'fof(q,{"question" if candidate == None else "conjecture"},({fol_to_tptp(question_lf)})).\n')
	prover_file.close()

	if candidate == None:
		if use_vampire:
			result = subprocess.run(['./vampire_rel__', '--mode', 'casc', '-t', '300', '-av', 'off', '-qa', 'answer_literal', filename], cwd=sys.argv[1], stdout=subprocess.PIPE, universal_newlines=True)
		else:
			result = subprocess.run(['PROVER/eprover', '-s', filename, '--soft-cpu-limit=60', '--answers'], cwd=sys.argv[1], stdout=subprocess.PIPE, universal_newlines=True)
	else:
		if use_vampire:
			result = subprocess.run(['./vampire_rel__', '--mode', 'casc', '-t', '30', filename], cwd=sys.argv[1], stdout=subprocess.PIPE, universal_newlines=True)
		else:
			result = subprocess.run(['PROVER/eprover', '-s', filename, '--soft-cpu-limit=10'], cwd=sys.argv[1], stdout=subprocess.PIPE, universal_newlines=True)
	answers = []
	found_proof = False
	for line in result.stdout.split('\n'):
		if candidate == None and (line.startswith('# SZS answers Tuple [[') or line.startswith('% SZS answers Tuple [[')):
			index = line.find(']')
			answer_candidates = line[len('# SZS answers Tuple [['):index].split(',')
			print(answer_candidates)
			for answer_candidate in answer_candidates:
				answer_candidate = answer_candidate.strip()
				if answer_candidate.startswith('esk') or answer_candidate.startswith('sK'):
					continue
				if answer_candidate.startswith('num'):
					answers.append(answer_candidate[3:].replace('_','.'))
				else:
					answers.append(answer_candidate)
		elif candidate != None and (line.startswith('# SZS status Theorem') or line.startswith('% SZS status Theorem')):
			found_proof = True
			print('found proof')
			break
	os.remove(sys.argv[1] + filename)
	if candidate != None:
		return found_proof
	return answers


if len(sys.argv) < 3:
	print("Missing arguments")
	print("Usage: python run_thm_prover.py [path to theorem prover directory] [output file] [--use-vampire] [--conjectures]")
	sys.exit(1)

use_vampire = False
conjectures = False
for i in range(3, len(sys.argv)):
	if sys.argv[i] == '--use-vampire':
		use_vampire = True
	elif sys.argv[i] == '--conjectures':
		conjectures = True

test_file = open('fictionalgeoqa_parse_outputs.txt', 'r')
output_file = open(sys.argv[2], 'w')
line_number = 1
for line in test_file:
	if line_number < 166 or line_number > 500:
		output_file.write('\n')
		line_number += 1
		continue
	if line.startswith('<no answer>'):
		output_file.write('\n')
		line_number += 1
		continue
	has_greatest = 'greatest' in line
	has_least = 'least' in line
	example = json.loads('[' + line + ']')[0]

	print(f"LINE NUMBER: {line_number}")
	formulas = []
	numbers = set()
	names = set()
	skip = False
	for sentence in example:
		if len(sentence) == 0:
			skip = True
			break
		logical_form = sentence[0]
		if logical_form[0] == '^':
			logical_form = '?' + logical_form[1:]
		formula = do_parse_fol_from_tptp(logical_form.replace('^','!').replace(')(',','))
		formulas.append(first_orderize(formula, numbers, names))
	if skip:
		output_file.write('\n')
		line_number += 1
		continue

	# run the theorem prover to answer the question
	question_count = 1
	presuppositions = []
	for i in range(question_count):
		question = formulas[-i-1]
		presuppositions.append(question)
		if not conjectures:
			new_var = max_variable(question) + 1
			if type(question.operand) == FOLAnd:
				new_operands = question.operand.operands
			else:
				new_operands = [question.operand]
			new_question = FOLExists(new_var,
					FOLExists(question.variable,
						FOLAnd([FOLOr([
							FOLEquals(FOLVariable(question.variable), FOLVariable(new_var)),
							FOLFuncApplication("element_of", [FOLVariable(new_var), FOLVariable(question.variable)]),
							FOLExists(new_var+1, FOLAnd([
								FOLFuncApplication("name", [FOLVariable(new_var+1)]),
								FOLEquals(FOLFuncApplication("arg1",[FOLVariable(new_var+1)]), FOLVariable(question.variable)),
								FOLEquals(FOLFuncApplication("arg2",[FOLVariable(new_var+1)]), FOLVariable(new_var))
							])
						)])] + new_operands)
					)
				)
			formulas[-i-1] = new_question
	for i in range(question_count):
		answers = []
		if conjectures:
			for number in numbers:
				if run_thm_prover(formulas, presuppositions, question_count, formulas[-i-1], numbers, use_vampire, number):
					answers.append(number)
			for name in names:
				question_lf = formulas[-i-1]
				if run_thm_prover(formulas, presuppositions, question_count, question_lf, numbers, use_vampire, name):
					answers.append(name)
				new_var = max_variable(question_lf) + 1
				if type(question_lf.operand) == FOLAnd:
					new_operands = question_lf.operand.operands
				else:
					new_operands = [question_lf.operand]
				question_lf = FOLExists(new_var, FOLExists(question_lf.variable, FOLAnd([
						FOLFuncApplication("element_of", [FOLVariable(new_var), FOLVariable(question_lf.variable)])
					] + new_operands)))
				if run_thm_prover(formulas, presuppositions, question_count, question_lf, numbers, use_vampire, name):
					answers.append(name)
		else:
			answers = run_thm_prover(formulas, presuppositions, question_count, formulas[-i-1], numbers, use_vampire, None)
		output_file.write(', '.join(answers) + '\n')
		output_file.flush()

	line_number += 1

output_file.close()
