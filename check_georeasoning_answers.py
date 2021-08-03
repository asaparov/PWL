import json
import sys

def test_candidate(candidate, actual_answers):
	try:
		actual_answers.remove(candidate)
		return True
	except ValueError:
		return False

if len(sys.argv) < 2:
	print("check_georeasoning_answers [answers file]")
	sys.exit(1)

try:
	actual_file = open('georeasoning.jsonl')
	predicted_file = open(sys.argv[1])
except IOError as e:
	print(e)
	sys.exit(1)

total = 0
correct = 0
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

		if len(actual_answers) != 0:
			# not all correct answers were predicted
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' INCORRECT")
		elif start != len(predicted_answer):
			# there are more predicted answers than actual answers
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' INCORRECT")
		else:
			print("[" + str(total + 1) + "] Predicted: '" + predicted_line.strip() + "', Actual: '" + actual_answer + "' CORRECT")

		if len(actual_answers) == 0 and start == len(predicted_answer):
			correct += 1
		total += 1

print("Accuracy: " + str(correct) + "/" + str(total) + " = " + str(correct/total))
