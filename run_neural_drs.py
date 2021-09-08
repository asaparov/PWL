import sys
import os
import glob
import json
import subprocess
import threading
import re
from drs import *

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def normalize_names_and_numbers(formula, numbers, names):
	def func(f):
		nonlocal numbers
		if type(f) != FOLFuncApplication:
			return f.apply(func)

		if f.function == 'Name' and f.args[1].constant != '"?"':
			names.add(f.args[1].constant)
			return FOLFuncApplication('name', f.args)
		elif f.function == 'Quantity' or f.function == 'YearOfCentury':
			if type(f.args[1]) == FOLConstant:
				literal = f.args[1].constant[1:-1]
				if isfloat(literal):
					numbers.add(literal)
				return FOLFuncApplication('quantity', [f.args[0], FOLConstant('num' + literal)])
			else:
				return FOLFuncApplication('quantity', f.args)
		elif f.function.isdigit():
			numbers.add(f.function)
			return FOLFuncApplication('quantity', [f.args[0], FOLConstant('num' + f.function)])
		elif f.function == 'LEQ':
			return FOLOr([FOLEquals(f.args[0], f.args[1]), FOLFuncApplication('greater_than', [f.args[1], f.args[0]])])
		elif f.function == 'LES':
			return FOLFuncApplication('greater_than', [f.args[1], f.args[0]])
		else:
			function = f.function.replace('-','')
			if function == 'quantity':
				function = 'quantity_event'
			elif function == 'name':
				function = 'name_event'
			function = function[0].lower() + function[1:]
			return FOLFuncApplication(function, f.args)

	return func(formula)

def remove_lambda_variable(formula):
	new_variable = max_variable(formula) + 1

	counter = 0
	remove_quantifier = 0
	def remove_lambda_var_apply(f):
		nonlocal counter
		nonlocal remove_quantifier
		if type(f) == FOLAnd:
			new_operands = []
			for operand in f.operands:
				new_operand = remove_lambda_var_apply(operand)
				if type(new_operand) == FOLFuncApplication and new_operand.function in ['name', 'value'] and new_operand.args[1].constant == '"?"':
					counter += 1
					remove_quantifier = new_operand.args[0].variable
				elif type(new_operand) == FOLFuncApplication and new_operand.function == 'entity' and new_operand.args[0].variable == remove_quantifier:
					continue
				elif type(new_operand) == FOLFuncApplication and new_operand.function == 'quantity' and type(new_operand.args[1]) == FOLConstant and new_operand.args[1].constant == 'num?':
					counter += 1
					remove_quantifier = new_operand.args[0].variable
				else:
					new_operands.append(new_operand)
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

	new_formula = remove_lambda_var_apply(formula)

	if counter == 0:
		raise Exception("remove_lambda_variable ERROR: Unable to find lambda variable.")
	elif counter != 1:
		raise Exception("remove_lambda_variable ERROR: Found multiple candidates for lambda variable.")
	else:
		return FOLExists(new_variable + 1,
			FOLExists(new_variable,
				FOLAnd([
					FOLOr([
						FOLFuncApplication('name', [FOLVariable(new_variable), FOLVariable(new_variable + 1)]),
						FOLFuncApplication('value', [FOLVariable(new_variable), FOLVariable(new_variable + 1)]),
						FOLFuncApplication('quantity', [FOLVariable(new_variable), FOLVariable(new_variable + 1)])
					]),
					new_formula
				])
			)
		)

def run_thm_prover(formulas, question_count, question_lf, numbers, use_vampire, candidate):
	if candidate != None:
		print("attempting proof with candidate " + candidate)
	filename = f'input_sentences.{threading.current_thread().ident}.p'
	prover_file = open(sys.argv[2] + filename, 'w')
	prover_file.write('fof(a0,axiom,![A]:![B]:(?[X]:(be(X) & theme(X,A) & coTheme(X,B)) => (A=B))).\n')
	prover_file.write('fof(a1,axiom,![X1]:![L2]:(?[L1]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & unit(L2,"kilometer") & ?[V1,V2]:(quantity(L1,V1) & quantity(L2,V2) & greater_than(V1,V2))) => ?[X4]:(longer(X4) & value(X4,L2) & theme(X4,X1)))).\n')
	prover_file.write('fof(a2,axiom,![X1]:![L2]:(?[L1]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & unit(L2,"kilometer") & ?[V1,V2]:(quantity(L1,V1) & quantity(L2,V2) & greater_than(V2,V1))) => ?[X4]:(shorter(X4) & value(X4,L2) & theme(X4,X1)))).\n')
	prover_file.write('fof(a3,axiom,![X1]:![L2]:(?[L1]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & unit(L2,"meter") & ?[V1,V2]:(quantity(L1,V1) & quantity(L2,V2) & greater_than(V1,V2))) => ?[X4]:(longer(X4) & value(X4,L2) & theme(X4,X1)))).\n')
	prover_file.write('fof(a4,axiom,![X1]:![L2]:(?[L1]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & unit(L2,"meter") & ?[V1,V2]:(quantity(L1,V1) & quantity(L2,V2) & greater_than(V2,V1))) => ?[X4]:(shorter(X4) & value(X4,L2) & theme(X4,X1)))).\n')
	prover_file.write('fof(a5,axiom,![L1]:![L2]:(?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V1,V2)) => ?[X4]:(greater(X4) & value(X4,L2) & theme(X4,L1)))).\n')
	prover_file.write('fof(a6,axiom,![L1]:![L2]:(?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V2,V1)) => ?[X4]:(smaller(X4) & value(X4,L2) & theme(X4,L1)))).\n')
	prover_file.write('fof(a7,axiom,less=smaller).\n')
	prover_file.write('fof(a8,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & length(L2) & partOf(L2,X2) & unit(L2,"kilometer") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V1,V2)))) => ?[M]:(longest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a9,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(length(L1) & partOf(L1,X1) & unit(L1,"kilometer") & length(L2) & partOf(L2,X2) & unit(L2,"kilometer") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V2,V1)))) => ?[M]:(shortest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a10,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(length(L1) & partOf(L1,X1) & unit(L1,"meter") & length(L2) & partOf(L2,X2) & unit(L2,"meter") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V1,V2)))) => ?[M]:(longest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a11,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(length(L1) & partOf(L1,X1) & unit(L1,"meter") & length(L2) & partOf(L2,X2) & unit(L2,"meter") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V2,V1)))) => ?[M]:(shortest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a12,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(population(L1) & partOf(L1,X1) & population(L2) & partOf(L2,X2) & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V1,V2)))) => ?[M]:(largest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a13,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(population(L1) & partOf(L1,X1) & population(L2) & partOf(L2,X2) & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V2,V1)))) => ?[M]:(smallest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a14,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(area(L1) & partOf(L1,X1) & unit(L1,"kilometer") & area(L2) & partOf(L2,X2) & unit(L2,"kilometer") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V1,V2)))) => ?[M]:(largest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a15,axiom,![X1]:(![X2]:(~(X1=X2) => ?[L1]:?[L2]:(area(L1) & partOf(L1,X1) & unit(L1,"kilometer") & area(L2) & partOf(L2,X2) & unit(L2,"kilometer") & ?[V1,V2]:(value(L1,V1) & value(L2,V2) & greater_than(V2,V1)))) => ?[M]:(smallest(M) & theme(M,X1)))).\n')
	prover_file.write('fof(a16,axiom,largest=biggest).\n')
	prover_file.write('fof(a17,axiom,![X]:![Y]:(location(Y,X) => ?[H]:(have(H) & pivot(H,X) & theme(H,Y)))).\n')
	prover_file.write('fof(a18,axiom,run=flow).\n')
	prover_file.write('fof(a19,axiom,![L]:![X]:((length(L) & partOf(L,X)) <=> (long(L) & theme(L,X)))).\n')
	axiom_count = 0
	for a in numbers:
		for b in numbers:
			if float(a) > float(b):
				prover_file.write(f'fof(g{axiom_count},axiom,greater_than(num{a},num{b})).\n')
				axiom_count += 1
			elif a != b:
				prover_file.write(f'fof(g{axiom_count},axiom,~greater_than(num{a},num{b})).\n')
				axiom_count += 1
				prover_file.write(f'fof(g{axiom_count},axiom,~(num{a}=num{b})).\n')
				axiom_count += 1
	for j in range(len(formulas) - question_count):
		prover_file.write(f'fof(c{j},axiom,({fol_to_tptp(formulas[j])})).\n')
	if candidate != None:
		if isfloat(candidate):
			question_lf = substitute(question_lf.operand, FOLVariable(question_lf.variable), FOLConstant('num' + candidate))
		else:
			question_lf = substitute(question_lf.operand, FOLVariable(question_lf.variable), FOLConstant(candidate))
	prover_file.write(f'fof(q,{"question" if candidate == None else "conjecture"},({fol_to_tptp(question_lf)})).\n')
	prover_file.close()

	if candidate == None:
		if use_vampire:
			result = subprocess.run(['./vampire_rel__', '--mode', 'casc', '-t', '300', '-av', 'off', '-qa', 'answer_literal', filename], cwd=sys.argv[2], stdout=subprocess.PIPE, universal_newlines=True)
		else:
			result = subprocess.run(['PROVER/eprover', '-s', filename, '--soft-cpu-limit=60', '--answers'], cwd=sys.argv[2], stdout=subprocess.PIPE, universal_newlines=True)
	else:
		if use_vampire:
			result = subprocess.run(['./vampire_rel__', '--mode', 'casc', '-t', '30', filename], cwd=sys.argv[2], stdout=subprocess.PIPE, universal_newlines=True)
		else:
			result = subprocess.run(['PROVER/eprover', '-s', filename, '--soft-cpu-limit=10'], cwd=sys.argv[2], stdout=subprocess.PIPE, universal_newlines=True)
	answers = []
	found_proof = False
	for line in result.stdout.split('\n'):
		if candidate == None and (line.startswith('# SZS answers Tuple [[') or line.startswith('% SZS answers Tuple [[')):
			index = line.find(']')
			candidates = line[len('# SZS answers Tuple [['):index].split(',')
			print(candidates)
			for answer_candidate in candidates:
				answer_candidate = answer_candidate.strip()
				if answer_candidate.startswith('esk') or answer_candidate.startswith('sK'):
					continue
				if answer_candidate.startswith('num'):
					answers.append(answer_candidate[3:])
				else:
					answers.append(answer_candidate)
		elif candidate != None and (line.startswith('# SZS status Theorem') or line.startswith('% SZS status Theorem')):
			found_proof = True
			print('found proof')
			break
	os.remove(sys.argv[2] + filename)
	if candidate != None:
		return found_proof
	return answers

def count_missing_names(first, second):
	if len(first.presuppositions) != len(second.presuppositions) or first.referents != second.referents:
		raise Exception('count_missing_names ERROR: Input DRS have different structures.')
	n = 0
	for i in range(len(first.presuppositions)):
		n += count_missing_names(first.presuppositions[i], second.presuppositions[i])
	i, j = 0, 0
	while i < len(first.conditions) and j < len(second.conditions):
		if type(second.conditions[j]) == DRSFuncApplication and second.conditions[j].function == 'name' and first.conditions[i] != second.conditions[j]:
			n += 1
			j += 1
		elif type(second.conditions[j]) == DRSNegation:
			if first.conditions[i] != DRSNegation:
				raise Exception('count_missing_names ERROR: Input DRS have different structures.')
			n += count_missing_names(first.conditions[i].operand, second.conditions[j].operand)
			i += 1
			j += 1
		elif type(second.conditions[j]) == DRSImplication:
			if first.conditions[i] != DRSImplication:
				raise Exception('count_missing_names ERROR: Input DRS have different structures.')
			n += count_missing_names(first.conditions[i].antecedent, second.conditions[j].antecedent)
			n += count_missing_names(first.conditions[i].consequent, second.conditions[j].consequent)
			i += 1
			j += 1
		elif type(second.conditions[j]) in [DRSFuncApplication, DRSEquals]:
			if first.conditions[i] != second.conditions[j]:
				raise Exception('count_missing_names ERROR: Input DRS have different structures.')
	if i < len(first.conditions):
		raise Exception('count_missing_names ERROR: Input DRS have different structures.')
	while j < len(second.conditions):
		if type(second.conditions[j]) == DRSFuncApplication and second.conditions[j].function == 'name':
			n += 1
			j += 1
		else:
			raise Exception('count_missing_names ERROR: Input DRS have different structures.')
	return n

def substitute_strings(drs, map):
	def visit(exp):
		if type(exp) == DRSString:
			if exp.string in map:
				exp.string = map[exp.string]
		else:
			exp.visit(visit)
	visit(drs)

fictional_cities = ["Wilgalway", "Elfincaster", "Ballaterma", "Emprizia", "Veronizzia", "Baritolloti", "Brescitabia", "Dorbincester", "Gilkarney", "Waterghal", "Lomtrieste", "Albucapua", "Messinitina", "Lomberona", "Geoturin", "Koruhashi", "Dogayashi", "Kyoukashino", "Agarikoshi", "Kagenegawa", "Yotsuyamashi", "Sennouhama", "Mbalam", "Marba Adama", "Alharwassa", "Sharhamane", "Saare Daewa", "Gondarbawa"]
fictional_rivers = ["Merasardu River", "Merafagole River", "Mbalam", "Kolufori River", "River Giffeleney", "River Wulstershire", "River Elsuir", "Begliomento", "Terravipacco", "Pernatisone", "Chetershka", "Kolonolga", "Voronolga", "Brasyugan", "Karbankaya"]
fictional_states = ["Merafagole", "Wulstershire", "Kangoyaken", "Gunmaishyu", "Ibarakishyo", "Toyusuma", "Senkuoka", "Grappulia", "Cartabitan", "Regnobenoa", "Bascilitina", "Umbriazzo", "Lordanidia", "Estmolise", "Baermerick", "Eldmunster", "Efanangole", "Bolurofi", "Kolufori", "Timbuqt"]
fictional_provinces = ["Voronolga", "Abdorostan", "Fordgorod", "Galininograd", "Puotorsk", "Getarovo", "Bripetrsk"]
fictional_countries = ["Efanangole", "Gyoshoru", "Dogoreoku", "Catardinia", "Gaelgiland", "Mofubali", "Bolurofi", "Bievorsk"]
real_cities = ["Boston", "Cardiff", "Naples", "Turin", "Lyon", "Calgary", "Montreal", "Monterey", "Tyre", "Birmingham", "Stockholm", "Riga", "Portland", "Wales", "Albany", "Mumbai", "Brisbane", "Tonga", "Damascus", "Fukuoka", "Vancouver", "Lagos", "Chicago", "Buenos Aires", "Palermo", "Nairobi", "Burkina Faso", "Kathmandu"]
real_rivers = ["Mississippi River", "Minnesota River", "Chicago", "Alabama River", "River Shannon", "River Barrow", "River Till", "Reno", "Po", "Danube", "Ob", "Oder", "Volga", "Congo", "Yukon"]
real_states = ["Minnesota", "Barrow", "Queensland", "Tasmania", "Victoria", "Gujarat", "Goa", "Kerala", "Rajasthan", "Bavaria", "Hesse", "Saxony", "Bremen", "Hidalgo", "Oaxaca", "Puebla", "Mexico", "Georgia", "Alabama", "Sonora"]
real_provinces = ["Volga", "Nara", "Ontario", "Alberta", "Manitoba", "Fujian", "Chiba"]
real_countries = ["Mexico", "Korea", "Malaysia", "Lithuania", "Switzerland", "Bolivia", "Zimbabwe", "Poland"]

str_map = {}
reverse_map = {}
for i in range(len(fictional_cities)):
	str_map[re.escape(fictional_cities[i])] = real_cities[i]
	reverse_map['"' + real_cities[i].lower().replace(' ', '~') + '"'] = '"' + fictional_cities[i].lower().replace(' ', '~') + '"'
for i in range(len(fictional_rivers)):
	str_map[re.escape(fictional_rivers[i])] = real_rivers[i]
	reverse_map['"' + real_rivers[i].lower().replace(' ', '~') + '"'] = '"' + fictional_rivers[i].lower().replace(' ', '~') + '"'
for i in range(len(fictional_states)):
	str_map[re.escape(fictional_states[i])] = real_states[i]
	reverse_map['"' + real_states[i].lower().replace(' ', '~') + '"'] = '"' + fictional_states[i].lower().replace(' ', '~') + '"'
for i in range(len(fictional_provinces)):
	str_map[re.escape(fictional_provinces[i])] = real_provinces[i]
	reverse_map['"' + real_provinces[i].lower().replace(' ', '~') + '"'] = '"' + fictional_provinces[i].lower().replace(' ', '~') + '"'
for i in range(len(fictional_countries)):
	str_map[re.escape(fictional_countries[i])] = real_countries[i]
	reverse_map['"' + real_countries[i].lower().replace(' ', '~') + '"'] = '"' + fictional_countries[i].lower().replace(' ', '~') + '"'
pattern = re.compile("|".join(str_map.keys()))


if len(sys.argv) < 4:
	print("Missing arguments")
	print("Usage: python run_neural_drs.py [Neural DRS directory] [path to theorem prover directory] [output file] [--use-vampire] [--conjectures]")
	sys.exit(1)

use_vampire = False
use_marian = True
conjectures = False
for i in range(4, len(sys.argv)):
	if sys.argv[i] == '--use-vampire':
		use_vampire = True
	elif sys.argv[i] == '--conjectures':
		conjectures = True

parser_env = os.environ.copy()
new_python_path = sys.argv[1] + "/DRS_parsing/:" + sys.argv[1] + "/DRS_parsing/evaluation/"
if "PYTHONPATH" in parser_env:
	parser_env["PYTHONPATH"] = new_python_path + ":" + parser_env["PYTHONPATH"]
else:
	parser_env["PYTHONPATH"] = new_python_path

test_file = open('fictionalgeoqa.jsonl', 'r')
output_file = open(sys.argv[3], 'w')
line_number = 1
for line in test_file:
	if line_number < 141 or line_number > 150:
		line_number += 1
		continue
	example = json.loads(line)
	input_filename = f'input_sentences.{threading.current_thread().ident}'
	input_file = open(sys.argv[1] + input_filename + '.txt', 'w')
	theory = example["theory"].split('. ')
	for sentence in theory:
		input_file.write(sentence)
		if sentence[-1] != '.' and sentence[-1] != '?' and sentence[-1] != '!':
			input_file.write('.')
		input_file.write('\n')

	for question in example["questions"].values():
		input_file.write(question['question'])
	input_file.close()

	if not use_marian:
		other_filename = f'other_sentences.{threading.current_thread().ident}'
		other_file = open(sys.argv[1] + other_filename + '.txt', 'w')
		sub_counts = []
		for sentence in theory:
			# replace all fictional place names with names that exist in the decoder's target vocabulary
			new_sentence, n = pattern.subn(lambda m: str_map[re.escape(m.group(0))], sentence)
			sub_counts.append(n)
			other_file.write(new_sentence)
			if new_sentence[-1] != '.' and new_sentence[-1] != '?' and new_sentence[-1] != '!':
				other_file.write('.')
			other_file.write('\n')
		for question in example["questions"].values():
			question_text = question['question']
			question_text, n = pattern.subn(lambda m: str_map[re.escape(m.group(0))], question_text)
			sub_counts.append(n)
			other_file.write(question_text)
		other_file.close()

	if use_marian:
		absolute_path = os.path.abspath(sys.argv[1])
		subprocess.call(['./src/marian_scripts/extract_ling_features.sh', absolute_path + '/' + input_filename + '.txt'], env=parser_env, cwd=sys.argv[1])
		clem_file = open(sys.argv[1] + input_filename + '.clem', "w")
		subprocess.call(['python', 'src/merge_tags.py', '-f', input_filename + '.txt.feat.lem', '--char_exts', '.feat.lem'], env=parser_env, cwd=sys.argv[1], stdout=clem_file)
		clem_file.close()
		subprocess.call(['./src/marian_scripts/parse_raw_text.sh', 'config/marian/best_gold_silver.sh', 'models/marian/best_gold_silver.npz', input_filename + '.out', input_filename + '.txt', input_filename + '.clem'], env=parser_env, cwd=sys.argv[1])
		drs_file = open(sys.argv[1] + input_filename + '.out.res', 'r')
	else:
		subprocess.call(['./src/allennlp_scripts/parse.sh', input_filename + '.txt', 'models/allennlp/bert_char_1enc.tar.gz', 'vocabs/allennlp/tgt_bert_char_1enc.txt'], env=parser_env, cwd=sys.argv[1])
		subprocess.call(['./src/allennlp_scripts/parse.sh', other_filename + '.txt', 'models/allennlp/bert_char_1enc.tar.gz', 'vocabs/allennlp/tgt_bert_char_1enc.txt'], env=parser_env, cwd=sys.argv[1])
		drs_file = open(sys.argv[1] + input_filename + '.txt.drs.out', 'r')
	drs_structs = []
	while True:
		drs_struct, eof = parse_drs(drs_file)
		if eof:
			break
		drs_structs.append(drs_struct)
	drs_file.close()

	if not use_marian:
		drs_file = open(sys.argv[1] + other_filename + '.txt.drs.out', 'r')
		other_drs_structs = []
		while True:
			drs_struct, eof = parse_drs(drs_file)
			if eof:
				break
			other_drs_structs.append(drs_struct)
		drs_file.close()

		for i in range(len(drs_structs)):
			print('count_missing_names with i = ' + str(i))
			n = count_missing_names(drs_structs[i], other_drs_structs[i])
			if n != sub_counts[i]:
				raise Exception(f'Unexpected number of missing place names. (expected {sub_counts[i]}, found {n})')
			substitute_strings(other_drs_structs[i], reverse_map)
			drs_structs[i] = other_drs_structs[i]

	for filename in glob.glob(sys.argv[1] + input_filename + '*'):
		os.remove(filename)

	formulas = []
	numbers = set()
	names = set()
	skip = False
	for drs in drs_structs:
		if drs == None:
			skip = True
			break
		resolve_presuppositions(drs)
		remove_tense_info(drs)
		remove_roles(drs)
		formula = drs_to_fol(drs)
		if formula == True:
			continue
		elif formula == False:
			skip = True
			break
		formulas.append(normalize_names_and_numbers(formula, numbers, names))
		print(fol_to_tptp(formulas[-1]) + '\n')
	if skip:
		output_file.write('\n')
		output_file.flush()
		line_number += 1
		continue

	# run the theorem prover to answer the question
	question_count = len(example["questions"].values())
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
		answers = []
		if conjectures:
			for number in numbers:
				if run_thm_prover(formulas, question_count, formulas[-i-1], numbers, use_vampire, number):
					answers.append(number)
			for name in names:
				if run_thm_prover(formulas, question_count, formulas[-i-1], numbers, use_vampire, name):
					answers.append(name)
		else:
			answers = run_thm_prover(formulas, question_count, formulas[-i-1], numbers, use_vampire, None)
		output_file.write(', '.join(answers) + '\n')
		output_file.flush()

	line_number += 1

output_file.close()
