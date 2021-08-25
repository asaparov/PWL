class FOLFormula(object):
	def __init__(self):
		return

class FOLAnd(FOLFormula):
	def __init__(self, conjuncts):
		self.operands = conjuncts

	def __eq__(self, other):
		if type(other) != FOLAnd:
			return False
		if len(self.operands) != len(other.operands):
			return False
		for i in range(len(self.operands)):
			if self.operands[i] != other.operands[i]:
				return False
		return True

	def __ne__(self, other):
		if type(other) != FOLAnd:
			return True
		if len(self.operands) != len(other.operands):
			return True
		for i in range(len(self.operands)):
			if self.operands[i] != other.operands[i]:
				return True
		return False

	def apply(self, func):
		operands = []
		for operand in self.operands:
			operands.append(func(operand))
		return FOLAnd(operands)

	def visit(self, visit):
		for operand in self.operands:
			visit(operand)

class FOLOr(FOLFormula):
	def __init__(self, disjuncts):
		self.operands = disjuncts

	def __eq__(self, other):
		if type(other) != FOLOr:
			return False
		if len(self.operands) != len(other.operands):
			return False
		for i in range(len(self.operands)):
			if self.operands[i] != other.operands[i]:
				return False
		return True

	def __ne__(self, other):
		if type(other) != FOLOr:
			return True
		if len(self.operands) != len(other.operands):
			return True
		for i in range(len(self.operands)):
			if self.operands[i] != other.operands[i]:
				return True
		return False

	def apply(self, func):
		operands = []
		for operand in self.operands:
			operands.append(func(operand))
		return FOLOr(operands)

	def visit(self, visit):
		for operand in self.operands:
			visit(operand)

class FOLNot(FOLFormula):
	def __init__(self, operand):
		self.operand = operand

	def __eq__(self, other):
		if type(other) != FOLNot:
			return False
		return self.operand == other.operand

	def __ne__(self, other):
		if type(other) != FOLNot:
			return True
		return self.operand != other.operand

	def apply(self, func):
		return FOLNot(func(self.operand))

	def visit(self, visit):
		visit(self.operand)

class FOLIfThen(FOLFormula):
	def __init__(self, antecedent, consequent):
		self.antecedent = antecedent
		self.consequent = consequent

	def __eq__(self, other):
		if type(other) != FOLIfThen:
			return False
		return self.antecedent == other.antecedent and self.consequent == other.consequent

	def __ne__(self, other):
		if type(other) != FOLIfThen:
			return True
		return self.antecedent != other.antecedent or self.consequent != other.consequent

	def apply(self, func):
		return FOLIfThen(func(self.antecedent), func(self.consequent))

	def visit(self, visit):
		visit(self.antecedent)
		visit(self.consequent)

class FOLEquals(FOLFormula):
	def __init__(self, left, right):
		self.left = left
		self.right = right

	def __eq__(self, other):
		if type(other) != FOLEquals:
			return False
		return self.left == other.left and self.right == other.right

	def __ne__(self, other):
		if type(other) != FOLEquals:
			return True
		return self.left != other.left or self.right != other.right

	def apply(self, func):
		return FOLEquals(func(self.left), func(self.right))

	def visit(self, visit):
		visit(self.left)
		visit(self.right)

class FOLForAll(FOLFormula):
	def __init__(self, variable, operand):
		self.variable = variable
		self.operand = operand

	def __eq__(self, other):
		if type(other) != FOLForAll:
			return False
		return self.variable == other.variable and self.operand == other.operand

	def __ne__(self, other):
		if type(other) != FOLForAll:
			return True
		return self.variable != other.variable or self.operand != other.operand

	def apply(self, func):
		return FOLForAll(self.variable, func(self.operand))

	def visit(self, visit):
		visit(self.operand)

class FOLExists(FOLFormula):
	def __init__(self, variable, operand):
		self.variable = variable
		self.operand = operand

	def __eq__(self, other):
		if type(other) != FOLExists:
			return False
		return self.variable == other.variable and self.operand == other.operand

	def __ne__(self, other):
		if type(other) != FOLExists:
			return True
		return self.variable != other.variable or self.operand != other.operand

	def apply(self, func):
		return FOLExists(self.variable, func(self.operand))

	def visit(self, visit):
		visit(self.operand)

class FOLFuncApplication(FOLFormula):
	def __init__(self, function, args):
		self.function = function
		self.args = args

	def __eq__(self, other):
		if type(other) != FOLFuncApplication:
			return False
		if len(self.args) != len(other.args):
			return False
		if self.function != other.function:
			return False
		for i in range(len(self.args)):
			if self.args[i] != other.args[i]:
				return False
		return True

	def __ne__(self, other):
		if type(other) != FOLFuncApplication:
			return True
		if len(self.args) != len(other.args):
			return True
		if self.function != other.function:
			return True
		for i in range(len(self.args)):
			if self.args[i] != other.args[i]:
				return True
		return False

	def apply(self, func):
		return FOLFuncApplication(self.function, [func(arg) for arg in self.args])

	def visit(self, visit):
		for arg in self.args:
			visit(arg)

class FOLTerm(object):
	def __init__(self):
		return

class FOLVariable(FOLTerm):
	def __init__(self, variable):
		self.variable = variable

	def __eq__(self, other):
		if type(other) != FOLVariable:
			return True
		return self.variable == other.variable

	def __ne__(self, other):
		if type(other) != FOLVariable:
			return False
		return self.variable != other.variable

	def apply(self, func):
		return FOLVariable(self.variable)

	def visit(self, visit):
		return

class FOLConstant(FOLTerm):
	def __init__(self, constant):
		self.constant = constant

	def __eq__(self, other):
		if type(other) != FOLConstant:
			return True
		return self.constant == other.constant

	def __ne__(self, other):
		if type(other) != FOLConstant:
			return False
		return self.constant != other.constant

	def apply(self, func):
		return FOLConstant(self.constant)

	def visit(self, visit):
		return

class FOLNumber(FOLTerm):
	def __init__(self, number):
		self.number = number

	def __eq__(self, other):
		if type(other) != FOLNumber:
			return True
		return self.number == other.number

	def __ne__(self, other):
		if type(other) != FOLNumber:
			return False
		return self.number != other.number

	def apply(self, func):
		return FOLNumber(self.number)

	def visit(self, visit):
		return

def substitute(formula, src, dst):
	def apply_substitute(f):
		if f == src:
			return dst
		else:
			return f.apply(apply_substitute)
	return formula.apply(apply_substitute)

def max_variable(formula):
	max_var = 0
	def max_variable_visit(f):
		nonlocal max_var
		if type(f) == FOLVariable or type(f) == FOLForAll or type(f) == FOLExists:
			max_var = max(max_var, f.variable)
		f.visit(max_variable_visit)
	formula.visit(max_variable_visit)
	return max_var

def parse_fol_term_from_prolog(string, pos, var_map):
	end = pos
	while string[end].isdigit() or string[end].isalpha():
		end += 1
	if string[pos:end] in var_map:
		term = FOLVariable(var_map[string[pos:end]])
	else:
		term = FOLConstant(string[pos:end])
	return term, end

def parse_fol_term_array_from_prolog(string, pos, var_map):
	operands = []
	while True:
		(term, pos) = parse_fol_term_from_prolog(string, pos, var_map)
		operands.append(term)
		if string[pos] == ')':
			pos += 1
			break
		elif string[pos] != ',':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis or comma in argument array.")
		pos += 1
	return operands, pos

def parse_fol_formula_array_from_prolog(string, pos, var_map):
	operands = []
	while True:
		(formula, pos) = parse_fol_from_prolog(string, pos, var_map)
		operands.append(formula)
		if string[pos] == ')':
			pos += 1
			break
		elif string[pos] != ',':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis or comma in conjunction/disjunction.")
		pos += 1
	return operands, pos

def parse_fol_quantifier_from_prolog(string, pos, var_map):
	index = string.find(',', pos)
	var_name = string[pos:index]
	if var_name in var_map:
		raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Variable '{var_name}' already declared in this scope.")
	variable = len(var_map) + 1
	var_map[var_name] = variable
	(formula, pos) = parse_fol_from_prolog(string, index + 1, var_map)
	del var_map[var_name]
	if string[pos] != ')':
		raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis for operand of quantification.")
	return variable, formula, pos + 1

def parse_fol_from_prolog(string, pos, var_map):
	if string.startswith('some(', pos):
		(variable, formula, pos) = parse_fol_quantifier_from_prolog(string, pos + len('some('), var_map)
		return FOLExists(variable, formula), pos
	elif string.startswith('all(', pos):
		(variable, formula, pos) = parse_fol_quantifier_from_prolog(string, pos + len('all('), var_map)
		return FOLExists(variable, formula), pos
	elif string.startswith('not(', pos):
		(formula, pos) = parse_fol_from_prolog(string, pos + len('not('), var_map)
		if string[pos] != ')':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis for operand of negation.")
		return FOLNot(formula), pos + 1
	elif string.startswith('and(', pos):
		(operands, pos) = parse_fol_formula_array_from_prolog(string, pos + len('and('), var_map)
		return FOLAnd(operands), pos
	elif string.startswith('or(', pos):
		(operands, pos) = parse_fol_formula_array_from_prolog(string, pos + len('or('), var_map)
		return FOLOr(operands), pos
	elif string.startswith('eq(', pos):
		(arg1, pos) = parse_fol_term_from_prolog(string, pos + len('eq('), var_map)
		if string[pos] != ',':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a comma in equality declaration.")
		(arg2, pos) = parse_fol_term_from_prolog(string, pos + 1, var_map)
		if string[pos] != ')':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis in equality declaration.")
		return FOLEquals(arg1, arg2), pos + 1
	elif string.startswith('imp(', pos):
		(arg1, pos) = parse_fol_from_prolog(string, pos + len('imp('), var_map)
		if string[pos] != ',':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a comma in implication declaration.")
		(arg2, pos) = parse_fol_from_prolog(string, pos + 1, var_map)
		if string[pos] != ')':
			raise Exception(f"parse_fol_from_prolog ERROR at {pos+1}: Expected a closing parenthesis in implication declaration.")
		return FOLIfThen(arg1, arg2), pos + 1
	else:
		index = string.find('(', pos)
		function = string[pos:index]
		(operands, pos) = parse_fol_term_array_from_prolog(string, index + 1, var_map)
		return FOLFuncApplication(function, operands), pos

def do_parse_fol_from_prolog(string):
	var_map = {}
	formula, pos = parse_fol_from_prolog(string, 0, var_map)
	return formula

def fol_term_to_tptp(term):
	if type(term) == FOLVariable:
		return f'X{term.variable}'
	elif type(term) == FOLConstant:
		return term.constant
	elif type(term) == FOLNumber:
		return term.number
	else:
		raise Exception("fol_term_to_tptp ERROR: Unrecognized term type.")

def fol_to_tptp(formula):
	if type(formula) == FOLAnd:
		return '(' + ' & '.join([fol_to_tptp(operand) for operand in formula.operands]) + ')'
	elif type(formula) == FOLOr:
		return '(' + ' | '.join([fol_to_tptp(operand) for operand in formula.operands]) + ')'
	elif type(formula) == FOLNot:
		return '~' + fol_to_tptp(formula.operand)
	elif type(formula) == FOLIfThen:
		return '(' + fol_to_tptp(formula.antecedent) + ' => ' + fol_to_tptp(formula.consequent) + ')'
	elif type(formula) == FOLEquals:
		return '(' + fol_term_to_tptp(formula.left) + '=' + fol_term_to_tptp(formula.right) + ')'
	elif type(formula) == FOLForAll:
		return f'![X{formula.variable}]:(' + fol_to_tptp(formula.operand) + ')'
	elif type(formula) == FOLExists:
		return f'?[X{formula.variable}]:(' + fol_to_tptp(formula.operand) + ')'
	elif type(formula) == FOLFuncApplication:
		return formula.function + '(' + ','.join([fol_term_to_tptp(arg) for arg in formula.args]) + ')'
	else:
		raise Exception(f"fol_to_tptp ERROR: Unrecognized formula type {type(formula)}.")

