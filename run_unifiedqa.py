from transformers import AutoTokenizer, T5ForConditionalGeneration
import json

model_name = "allenai/unifiedqa-t5-large" # you can specify the model size here
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = T5ForConditionalGeneration.from_pretrained(model_name)

def run_model(input_string, **generator_args):
    input_ids = tokenizer.encode(input_string, return_tensors="pt")
    res = model.generate(input_ids, **generator_args)
    return tokenizer.batch_decode(res, skip_special_tokens=True)

test_file = open('georeasoning.jsonl')
line_number = 1
for line in test_file:
	if line_number < 191 or line_number > 210:
		line_number += 1
		continue
	example = json.loads(line)
	theory = example["theory"]
	for question in example["questions"].values():
		predicted = run_model(question["question"] + " \\\\n " + theory)[0]
		print(predicted + " (correct answer: " + question["answer"] + ")")
	line_number += 1
test_file.close()
