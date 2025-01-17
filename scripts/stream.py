import json
import re
import sys
from datetime import datetime, timezone
import uuid


def parse_and_emit_section(text):
    sections = re.split(
        r'(Successfully added logical form as an observation to the theory\.|Finished reading declarative sentences\. Attempting to answer question:)', text)

    line_id_map = {}
    buffer = ""

    for i, section in enumerate(sections):
        if section.strip() in [
            "Successfully added logical form as an observation to the theory.",
            "Finished reading declarative sentences. Attempting to answer question:"
        ]:
            if buffer.strip():
                process_section(buffer)
            buffer = section
        else:
            buffer += section

    if buffer.strip():
        process_section(buffer)


def process_section(current_section):
    lines = current_section.strip().split('\n')
    theory_log_prob = None

    for line in lines:
        if "Theory log probability:" in line:
            theory_log_prob = line.split("Theory log probability:")[1].strip()
            break

    for line_idx, line in enumerate(lines):
        if "Logical form:" in line:
            logical_form_id = str(uuid.uuid4())
            logical_form = line.split("Logical form:")[1].strip()
            logical_form_obj = {
                "id": logical_form_id,
                "type": "logical_form",
                "axiom": logical_form,
                "description": current_section[:100].strip(),
                "theory_log_probability": theory_log_prob,
                "created_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
            }

            # Create proof object before printing anything
            proof_obj = None
            proof_start = current_section.find(
                'Proof of newly-added logical form')
            if proof_start != -1:
                proof_section = current_section[proof_start:]
                proof_id = str(uuid.uuid4())
                steps = []

                proof_lines = proof_section.split('\n')
                for proof_line_idx, proof_line in enumerate(proof_lines):
                    if proof_line.strip() and '[' in proof_line:
                        match = re.match(
                            r'\s*\[(\d+)\]\s*(.*?)\s*by\s*(.*?)$', proof_line.strip())
                        if match:
                            step_id = str(uuid.uuid4())
                            step_num = int(match.group(1))
                            formula = match.group(2)
                            justification = "by " + match.group(3)
                            theory_objects = sorted(
                                set(re.findall(r'c[₀₁₂₃₄₅₆₇₈₉]+', formula)))
                            steps.append({
                                "id": step_id,
                                "step": step_num,
                                "formula": formula,
                                "justification": justification,
                                "theory_objects": theory_objects,
                            })

                if steps:
                    proof_obj = {
                        "id": proof_id,
                        "type": "proof",
                        "steps": steps,
                        "created_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
                    }
                    logical_form_obj["proof_id"] = proof_id

            # Print objects only after both are created and linked
            print(json.dumps(logical_form_obj, ensure_ascii=False))
            if proof_obj:
                print(json.dumps(proof_obj, ensure_ascii=False))
            break

    answer_match = re.search(
        r'Answer:\s*(\w+)\s*:\s*([\d.]+)\s*\w+\s*:\s*([\d.]+)', current_section)
    if answer_match:
        answer_id = str(uuid.uuid4())
        answer_obj = {
            "id": answer_id,
            "type": "answer",
            "result": answer_match.group(1),
            "true_rate": answer_match.group(2),
            "false_rate": answer_match.group(3),
            "created_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
        }
        print(json.dumps(answer_obj, ensure_ascii=False))


def main():
    buffer = ""
    try:
        while True:
            line = sys.stdin.readline()
            if not line:
                break
            buffer += line
            if "Successfully added logical form as an observation to the theory." in line or \
               "Finished reading declarative sentences. Attempting to answer question:" in line:
                parse_and_emit_section(buffer)
                buffer = ""

        if buffer.strip():
            parse_and_emit_section(buffer)
    except OSError as e:
        print(json.dumps({
            "error": str(e),
            "message": "Error reading input stream",
            "created_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
        }, ensure_ascii=False))

    stream_end_id = str(uuid.uuid4())
    stream_end_obj = {
        "id": stream_end_id,
        "type": "stream_end",
        "message": "Reasoning session completed",
        "created_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
    }
    print(json.dumps(stream_end_obj, ensure_ascii=False))


if __name__ == "__main__":
    main()
