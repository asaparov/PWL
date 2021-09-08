import re
from fol import *

def line_split(line):
    return re.findall(r'[^"\s]\S*|".+?"', line)

class DRSCondition(object):
	def __init__(self):
		return

class DRSTerm(object):
	def __init__(self):
		return

	def visit(self, visit):
		return

class DRSReference(DRSTerm):
	def __init__(self, referent_id):
		self.referent_id = referent_id

	def __eq__(self, other):
		return type(other) == DRSReference and self.referent_id == other.referent_id

	def __ne__(self, other):
		return type(other) != DRSReference or self.referent_id != other.referent_id

class DRSNumber(DRSTerm):
	def __init__(self, number):
		self.number = number

	def __eq__(self, other):
		return type(other) == DRSNumber and self.number == other.number

	def __ne__(self, other):
		return type(other) != DRSNumber or self.number != other.number

class DRSString(DRSTerm):
	def __init__(self, string):
		self.string = string

	def __eq__(self, other):
		return type(other) == DRSString and self.string == other.string

	def __ne__(self, other):
		return type(other) != DRSString or self.string != other.string

class DRSFuncApplication(DRSCondition):
	def __init__(self, function, args):
		self.function = function
		self.args = args

	def visit(self, visit):
		for arg in self.args:
			visit(arg)

	def substitute(self, map):
		for i in range(len(self.args)):
			if type(self.args[i]) == DRSReference and self.args[i].referent_id in map:
				self.args[i] = map[self.args[i].referent_id]

class DRSNegation(DRSCondition):
	def __init__(self, operand):
		self.operand = operand

	def visit(self, visit):
		visit(self.operand)

	def substitute(self, map):
		self.operand.substitute(map)

class DRSImplication(DRSCondition):
	def __init__(self, antecedent, consequent):
		self.antecedent = antecedent
		self.consequent = consequent

	def visit(self, visit):
		visit(self.antecedent)
		visit(self.consequent)

	def substitute(self, map):
		self.antecedent.substitute(map)
		self.consequent.substitute(map)

class DRSEquals(DRSCondition):
	def __init__(self, first, second):
		self.first = first
		self.second = second

	def visit(self, visit):
		visit(self.first)
		visit(self.second)

	def substitute(self, map):
		if type(self.first) == DRSReference and self.first.referent_id in map:
			self.first = map[self.first.referent_id]
		if type(self.second) == DRSReference and self.second.referent_id in map:
			self.second = map[self.second.referent_id]

class DRS(object):
	def __init__(self, referents, conditions, presuppositions, parent):
		self.referents = referents
		self.conditions = conditions
		self.presuppositions = presuppositions
		self.parent = parent # (inverse) subordination

	def visit(self, visit):
		for presupposition in self.presuppositions:
			visit(presupposition)
		for condition in self.conditions:
			visit(condition)

	def empty(self):
		return len(self.conditions) == 0 and len(self.presuppositions) == 0

	def substitute(self, map):
		for condition in self.conditions:
			condition.substitute(map)
		for presupposition in self.presuppositions:
			presupposition.substitute(map)

def get_free_variables(drs):
	declared_referents = []
	free_variables = set()
	def visit(exp):
		nonlocal declared_referents
		if type(exp) == DRS:
			old_len = len(declared_referents)
			declared_referents.extend(get_declared_variables(exp))

			for presupposition in exp.presuppositions:
				presupposition.visit(visit)
			for condition in exp.conditions:
				visit(condition)

			declared_referents = declared_referents[:old_len]
		elif type(exp) == DRSReference:
			if exp.referent_id not in declared_referents:
				free_variables.add(exp.referent_id)
		elif type(exp) == DRSImplication:
			old_len = len(declared_referents)
			declared_referents.extend(get_declared_variables(exp.antecedent))
			exp.visit(visit)
			declared_referents = declared_referents[:old_len]
		else:
			exp.visit(visit)
	visit(drs)
	return free_variables

def is_var_declared(drs, var):
	if var in drs.referents:
		return True
	for presupposition in drs.presuppositions:
		if is_var_declared(presupposition, var):
			return True
	return False

def get_declared_variables(drs):
	variables = list(drs.referents)
	for presupposition in drs.presuppositions:
		variables.extend(get_declared_variables(presupposition))
	return variables

def has_variable(drs, variables):
	found_var = False
	def visit(exp):
		nonlocal found_var
		if type(exp) == DRSReference and exp.referent_id in variables:
			found_var = True
		else:
			exp.visit(visit)
	visit(drs)
	return found_var

def parse_drs_term(string, ref_map):
	if string[0] == '"':
		return DRSString(string[1:-1])
	elif string[0].isdigit():
		return DRSNumber(string)
	else:
		if string not in ref_map:
			ref_map[string] = len(ref_map) + 1
		return DRSReference(ref_map[string])

def parse_drs(input_file):
	drs_map = {}
	ref_map = {}
	eof = False
	while True:
		line = input_file.readline()
		if line == '':
			eof = True
			break
		line = line.strip()
		print(line)
		if line == '':
			break

		tokens = line_split(line)
		if tokens[0] not in drs_map:
			drs_map[tokens[0]] = DRS([], [], [], None)
		curr_drs = drs_map[tokens[0]]

		if tokens[1] == 'REF':
			if tokens[2] not in ref_map:
				ref_map[tokens[2]] = len(ref_map) + 1
			curr_drs.referents.append(ref_map[tokens[2]])
		elif tokens[1] == 'EQU':
			arg1 = parse_drs_term(tokens[2], ref_map)
			arg2 = parse_drs_term(tokens[3], ref_map)
			curr_drs.conditions.append(DRSEquals(arg1, arg2))
		elif tokens[1] == 'NEQ':
			arg1 = parse_drs_term(tokens[2], ref_map)
			arg2 = parse_drs_term(tokens[3], ref_map)
			curr_drs.conditions.append(DRSNegation(DRS([], [DRSEquals(arg1, arg2)], [], None)))
		elif tokens[1] == 'PRESUPPOSITION':
			if tokens[2] not in drs_map:
				drs_map[tokens[2]] = DRS([], [], [], None)
			drs_map[tokens[2]].presuppositions.append(curr_drs)
			curr_drs.parent = drs_map[tokens[2]]
		elif tokens[1] == 'NOT':
			if tokens[2] not in drs_map:
				drs_map[tokens[2]] = DRS([], [], [], None)
			curr_drs.conditions.append(DRSNegation(drs_map[tokens[2]]))
			drs_map[tokens[2]].parent = curr_drs
		elif tokens[1] == 'IMP':
			if tokens[2] not in drs_map:
				drs_map[tokens[2]] = DRS([], [], [], None)
			if tokens[3] not in drs_map:
				drs_map[tokens[3]] = DRS([], [], [], None)
			curr_drs.conditions.append(DRSImplication(drs_map[tokens[2]], drs_map[tokens[3]]))
			drs_map[tokens[2]].parent = curr_drs
			drs_map[tokens[3]].parent = drs_map[tokens[2]]
		elif tokens[1] == 'CONDITION':
			if tokens[2] not in drs_map:
				drs_map[tokens[2]] = DRS([], [], [], None)
			curr_drs.conditions.append(DRSImplication(drs_map[tokens[2]], None))
			drs_map[tokens[2]].parent = curr_drs
		elif tokens[1] == 'CONSEQUENCE':
			if tokens[2] not in drs_map:
				drs_map[tokens[2]] = DRS([], [], [], None)
			for condition in curr_drs.parent.conditions:
				if type(condition) == DRSImplication and condition.antecedent == curr_drs:
					condition.consequent = drs_map[tokens[2]]
					drs_map[tokens[2]].parent = curr_drs
					break
		elif tokens[1].isupper() and tokens[1] not in ['APX', 'LEQ', 'LES', 'TPR', 'TAB', 'SZP', 'SZN', 'SXP', 'SXN', 'STI', 'STO', 'SY1', 'SY2', 'SXY']:
			# we throw an exception to make sure to properly handle all clauses
			raise Exception(f'Found unhandled clause {tokens[1]}.')
		else:
			args = []
			for token in tokens[2:]:
				args.append(parse_drs_term(token, ref_map))
			curr_drs.conditions.append(DRSFuncApplication(tokens[1], args))

	root_drs = []
	free_drs = []
	for drs_name, drs in drs_map.items():
		if drs.parent == None:
			free_variables = list(get_free_variables(drs))
			if len(free_variables) == 0:
				root_drs.append(drs)
			else:
				free_drs.append((drs_name, drs, free_variables))
	while len(free_drs) != 0:
		drs_name, drs, free_variables = free_drs.pop()
		while len(free_variables) != 0:
			# find the root DRS that declares this referent
			free_variable = free_variables.pop()
			found_declaration = False
			for i in range(len(root_drs)):
				if is_var_declared(root_drs[i], free_variable):
					drs.presuppositions.append(root_drs[i])
					root_drs[i].parent = drs
					del root_drs[i]
					found_declaration = True
					break
			if not found_declaration:
				inv_ref_map = {v:k for k, v in ref_map.items()}
				print(f'ERROR: Found free variable {free_variable} ({inv_ref_map[free_variable]}) when processing DRS {drs_name}\n')
				return None, eof
			free_variables = get_free_variables(drs)
		root_drs.append(drs)
	if len(root_drs) == 0:
		return None, eof
	elif len(root_drs) != 1:
		print('ERROR: Found multiple DRS roots.\n')
		return None, eof
	return root_drs[0], eof

def unify_term(first, second, ref_map):
	if type(first) == DRSReference:
		if first.referent_id in ref_map:
			return ref_map[first.referent_id] == second
		else:
			ref_map[first.referent_id] = second
			return True
	else:
		return first == second

def unify(first, second, ref_map):
	if type(first) != type(second):
		return False
	elif type(first) == DRSFuncApplication:
		if first.function != second.function or len(first.args) != len(second.args):
			return False
		for i in range(len(first.args)):
			if not unify_term(first.args[i], second.args[i], ref_map):
				return False
		return True
	elif type(first) == DRSEquals:
		return unify_term(first.first, second.first, ref_map) and unify_term(first.second, second.second, ref_map)
	elif type(first) == DRSNegation or type(first) == DRSImplication:
		return False
	else:
		raise Exception('unify ERROR: Unrecognized DRS condition type.')

def try_bind(src, dst):
	ref_map = {}
	for src_condition in src.conditions:
		# find a condition in `dst` that unifies with `src_condition`
		found_unification = False
		for dst_condition in dst.conditions:
			if unify(src_condition, dst_condition, ref_map):
				found_unification = True
				break
		if not found_unification:
			return False

	for referent in src.referents:
		if referent not in ref_map:
			dst.referents.append(referent)
	dst.substitute(ref_map)
	src.referents.clear()
	src.conditions.clear()
	return True

def try_accomodate(src, dst):
	# make sure no referent becomes free
	free_variables = get_free_variables(src)
	for var in free_variables:
		ancestor = dst
		is_free = True
		while ancestor != None:
			if is_var_declared(ancestor, var):
				is_free = False
				break
			ancestor = ancestor.parent
		if is_free:
			return False

	dst.referents.extend(src.referents)
	dst.conditions.extend(src.conditions)
	src.referents.clear()
	src.conditions.clear()
	return True

def resolve_presupposition(drs):
	# first try to bind this with the closest DRS in its projection line
	candidate = drs.parent
	projection_line = []
	while candidate != None:
		if try_bind(drs, candidate):
			# binding is successful
			return
		projection_line.append(candidate)
		candidate = candidate.parent

	# next try accomodating with the highest DRS in its projection line
	for candidate in reversed(projection_line):
		if try_accomodate(drs, candidate):
			# accomodation is successful
			return

def resolve_presuppositions(drs):
	for condition in drs.conditions:
		if type(condition) == DRSNegation:
			resolve_presuppositions(condition.operand)
		elif type(condition) == DRSImplication:
			resolve_presuppositions(condition.antecedent)
			resolve_presuppositions(condition.consequent)
		elif type(condition) not in [DRSFuncApplication, DRSEquals]: 
			raise Exception('Unrecognized DRS condition.')
	index = 0
	while index < len(drs.presuppositions):
		resolve_presuppositions(drs.presuppositions[index])

		if len(drs.presuppositions[index].presuppositions) == 0:
			resolve_presupposition(drs.presuppositions[index])
			if drs.presuppositions[index].empty():
				del drs.presuppositions[index]
				index -= 1
		index += 1

def get_tense_variables(drs, variables, tense_variables):
	for condition in drs.conditions:
		if type(condition) == DRSFuncApplication:
			if condition.function == 'Time' and type(condition.args[1]) == DRSReference and condition.args[1].referent_id in variables:
				tense_variables.add(condition.args[1].referent_id)
			elif condition.function == 'time' and type(condition.args[0]) == DRSReference and condition.args[0].referent_id in variables:
				tense_variables.add(condition.args[0].referent_id)
		elif type(condition) == DRSNegation:
			get_tense_variables(condition.operand, variables, tense_variables)
		elif type(condition) == DRSImplication:
			get_tense_variables(condition.antecedent, variables, tense_variables)
			get_tense_variables(condition.consequent, variables, tense_variables)
		elif type(condition) != DRSEquals:
			raise Exception('get_tense_variables ERROR: Unrecognized DRS condition type.')
	for presupposition in drs.presuppositions:
		get_tense_variables(presupposition, variables, tense_variables)

def remove_conditions(drs, variables_to_remove):
	index = 0
	while index < len(drs.conditions):
		condition = drs.conditions[index]
		if type(condition) in [DRSFuncApplication, DRSEquals]:
			if has_variable(condition, variables_to_remove):
				del drs.conditions[index]
				index -= 1
		elif type(condition) == DRSImplication:
			remove_conditions(condition.antecedent, variables_to_remove)
			remove_conditions(condition.consequent, variables_to_remove)
		elif type(condition) == DRSNegation:
			remove_conditions(condition.operand, variables_to_remove)
		else:
			raise Exception('remove_conditions ERROR: Unrecognized DRS condition type.')
		index += 1

	index = 0
	while index < len(drs.presuppositions):
		remove_conditions(drs.presuppositions[index], variables_to_remove)
		if drs.presuppositions[index].empty():
			del drs.presuppositions[index]
			index -= 1
		index += 1

def remove_tense_info(drs):
	for condition in drs.conditions:
		if type(condition) == DRSNegation:
			remove_tense_info(condition.operand)
		elif type(condition) == DRSImplication:
			remove_tense_info(condition.antecedent)
			remove_tense_info(condition.consequent)
		elif type(condition) not in [DRSFuncApplication, DRSEquals]:
			raise Exception('remove_tense_info ERROR: Unrecognized DRS condition type.')
	for presupposition in drs.presuppositions:
		remove_tense_info(presupposition)

	tense_variables = set()
	get_tense_variables(drs, drs.referents, tense_variables)
	remove_conditions(drs, tense_variables)
	drs.referents = [ref for ref in drs.referents if ref not in tense_variables]

def remove_roles(drs):
	def visit(exp):
		if type(exp) == DRSFuncApplication:
			index = 0
			while index < len(exp.args):
				arg = exp.args[index]
				if type(arg) == DRSString and len(arg.string) > 2 and arg.string[0].islower() and arg.string[1] == '.' and arg.string[2:].isdigit():
					del exp.args[index]
					index -= 1
				index += 1
		else:
			exp.visit(visit)
	visit(drs)

def drs_term_to_fol(term):
	if type(term) == DRSReference:
		return FOLVariable(term.referent_id)
	elif type(term) == DRSNumber:
		return FOLNumber(term.number)
	elif type(term) == DRSString:
		return FOLConstant('"' + term.string + '"')
	else:
		raise Exception('drs_to_fol ERROR: Unrecognized DRS term type.')

def drs_condition_to_fol(condition):
	if type(condition) == DRSFuncApplication:
		return FOLFuncApplication(condition.function, [drs_term_to_fol(arg) for arg in condition.args])
	elif type(condition) == DRSNegation:
		operand = drs_to_fol(condition.operand)
		if operand == True:
			return False
		elif operand == False:
			return True
		return FOLNot(operand)
	elif type(condition) == DRSImplication:
		formula = FOLIfThen(drs_conditions_to_fol(condition.antecedent), drs_to_fol(condition.consequent))
		if formula.consequent == True:
			return True
		elif formula.consequent == False:
			if formula.antecedent == True:
				return False
			elif formula.antecedent == False:
				return True
			else:
				formula = FOLNot(formula)
		if formula.antecedent == True:
			formula = formula.consequent
		elif formula.antecedent == False:
			return True
		for referent in condition.antecedent.referents:
			formula = FOLForAll(referent, formula)
		return formula
	elif type(condition) == DRSEquals:
		return FOLEquals(drs_term_to_fol(condition.first), drs_term_to_fol(condition.second))
	else:
		raise Exception('drs_to_fol ERROR: Unrecognized DRS condition type.')

def drs_conditions_to_fol(drs):
	conjuncts = []
	for condition in drs.conditions:
		conjunct = drs_condition_to_fol(condition)
		if conjunct == True:
			continue
		elif conjunct == False:
			return False
		conjuncts.append(conjunct)
	if len(conjuncts) == 0:
		return True
	elif len(conjuncts) == 1:
		return conjuncts[0]
	else:
		return FOLAnd(conjuncts)

def drs_to_fol(drs):
	formula = drs_conditions_to_fol(drs)
	if formula == True:
		return True
	elif formula == False:
		return False
	for referent in drs.referents:
		formula = FOLExists(referent, formula)
	return formula
