#!/bin/bash
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_file> <id>"
  exit 1
fi

input_file="$1"
id="$2"

./pwl_reasoner_dbg "$input_file" | tee >(python3 ./scripts/stream.py >> "/tmp/centaur/${id}.jsonl")
