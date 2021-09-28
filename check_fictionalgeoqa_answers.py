import json
import sys

def test_candidate(candidate, actual_answers):
	candidates = [candidate]
	if not candidate.startswith('river') and not candidate.endswith('river'):
		candidates.append(candidate + ' river')
		candidates.append('river ' + candidate)
	for c in candidates:
		try:
			actual_answers.remove(c)
			return True
		except ValueError:
			pass
	return False

if len(sys.argv) < 2:
	print("check_fictionalgeoqa_answers [answers file]")
	sys.exit(1)

try:
	actual_file = open('fictionalgeoqa.jsonl')
	predicted_file = open(sys.argv[1])
except IOError as e:
	print(e)
	sys.exit(1)

total = 0
correct = 0
correct_per_flag = {}
total_per_flag = {}
correct_per_flag_count = {}
total_per_flag_count = {}
while True:
	actual_line = actual_file.readline()
	predicted_line = predicted_file.readline()

	if len(actual_line) == 0:
		if len(predicted_line) != 0:
			print("WARNING: There are more predictions than questions.")
		break

	predicted_answer = predicted_line
	while True:
		old_predicted_answer = predicted_answer
		predicted_answer = predicted_line.strip()
		if len(predicted_answer) > 0 and predicted_answer[0] == '<':
			predicted_answer = predicted_answer[1:]
		if len(predicted_answer) > 0 and predicted_answer[-1] == '>':
			predicted_answer = predicted_answer[:-1]
		if len(predicted_answer) > 0 and predicted_answer[0] == '(':
			predicted_answer = predicted_answer[1:]
		if len(predicted_answer) > 0 and predicted_answer[-1] == ')':
			predicted_answer = predicted_answer[:-1]
		if len(predicted_answer) > 0 and predicted_answer[0] == '[':
			predicted_answer = predicted_answer[1:]
		if len(predicted_answer) > 0 and predicted_answer[-1] == ']':
			predicted_answer = predicted_answer[:-1]
		if predicted_answer == old_predicted_answer:
			break

	predicted_answer = predicted_answer.lower()
	if predicted_answer == 'no answer' or predicted_answer == 'none' or predicted_answer == 'nothing' or predicted_answer == 'empty':
		predicted_answer = ''

	example = json.loads(actual_line)
	for question in example["questions"].values():
		actual_answer = question["answer"]
		if actual_answer == "":
			actual_answers = []
		else:
			actual_answers = [x.lower().strip() for x in actual_answer.split(",")]

		start = 0
		next_index = 0
		while next_index < len(predicted_answer):
			index = predicted_answer.find(',', next_index)
			if index != -1:
				next_index = index + 1
			else:
				index = predicted_answer.find('.', next_index)
				if index != -1:
					next_index = index + 1
				else:
					index = predicted_answer.find('and', next_index)
					if index != -1:
						next_index = index + 3
					else:
						index = len(predicted_answer)
						next_index = index

			# check if predicted_answer[start:index] is a correct answer
			candidate = predicted_answer[start:index].strip()
			if test_candidate(candidate, actual_answers):
				start = next_index
				continue

			if candidate.startswith('the ') and test_candidate(candidate[4:], actual_answers):
				start = next_index
				continue
			elif candidate.startswith('a ') and test_candidate(candidate[2:], actual_answers):
				start = next_index
				continue

			if "answer_templates" in question:
				answer_templates = question["answer_templates"]
				found_matching_template = False
				for answer_template in answer_templates:
					ans_index = answer_template.find('{0}')
					if candidate.startswith(answer_template[:ans_index].lower()) and candidate.endswith(answer_template[(ans_index + 3):].lower()) and test_candidate(candidate[ans_index:-(len(answer_template) - ans_index - 3)].lower(), actual_answers):
						found_matching_template = True
						break
					if not found_matching_template and answer_template[-1] == '.':
						answer_template = answer_template[:-1]
						if candidate.startswith(answer_template[:ans_index].lower()) and candidate.endswith(answer_template[(ans_index + 3):].lower()) and test_candidate(candidate[ans_index:-(len(answer_template) - ans_index - 3)].lower(), actual_answers):
							found_matching_template = True
							break
				if found_matching_template:
					start = next_index
					continue

		if len(actual_answers) != 0:
			# not all correct answers were predicted
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' INCORRECT")
		elif start != len(predicted_answer):
			# there are more predicted answers than actual answers
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' INCORRECT")
		else:
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' CORRECT")

		is_correct = (len(actual_answers) == 0 and start == len(predicted_answer))
		if is_correct:
			correct += 1
		total += 1

		flag_count = 0
		if "flags" in example:
			flag_count = len(example["flags"])
			for flag in example["flags"]:
				if flag not in correct_per_flag:
					correct_per_flag[flag] = 0
					total_per_flag[flag] = 0
				if is_correct:
					correct_per_flag[flag] += 1
				total_per_flag[flag] += 1
		if flag_count not in correct_per_flag_count:
			correct_per_flag_count[flag_count] = 0
			total_per_flag_count[flag_count] = 0
		if is_correct:
			correct_per_flag_count[flag_count] += 1
		total_per_flag_count[flag_count] += 1

print("Accuracy: " + str(correct) + "/" + str(total) + " = " + str(float(correct)/total))
for (flag, total_flag) in total_per_flag.items():
	correct_flag = correct_per_flag[flag]
	print("[" + flag + "] " + str(correct_flag) + "/" + str(total_flag) + " = " + str(float(correct_flag)/total_flag))
	print("[all except " + flag + "] " + str(correct - correct_flag) + "/" + str(total - total_flag) + " = " + str(float(correct - correct_flag)/(total - total_flag)))
for (flag_count, total_flag) in sorted(total_per_flag_count.items()):
	correct_flag = correct_per_flag_count[flag_count]
	print("[" + str(flag_count) + " flags] " + str(correct_flag) + "/" + str(total_flag) + " = " + str(float(correct_flag)/(total_flag)))
