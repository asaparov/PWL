# Probabilistic Worldbuilding from Language

This code depends on OpenSSL as well as the following repositories: [core](https://github.com/asaparov/core), [math](https://github.com/asaparov/math), [hdp](https://github.com/asaparov/hdp), and [grammar](https://github.com/asaparov/grammar).

To use this code, download the repository and dependencies and run `make executive_test` and then `./executive_test`. In order to add paths to the `make` command to search for additional include files, run `make executive_test CPPFLAGS+='-I[1st additional path] -I[2nd additional path]'`.

### Console

To run the code in **_console_** mode, where the user can input custom sentences and inspect the learned theory and proofs, run `./executive_test console`.

## Experiments

### ProofWriter question-answering

Download the [ProofWriter data](https://allenai.org/data/proofwriter) and run
```bash
./executive_test proofwriter --data=[ProofWriter data filepath] --out=[output predicted answers filepath]
```
`ProofWriter data filepath` points to `OWA/birds-electricity/meta-test.jsonl` within the ProofWriter data.

### FictionalGeoQA question-answering

Download the [FictionalGeoQA data](https://github.com/asaparov/fictionalgeoqa) and run
```bash
./executive_test fictionalgeoqa --data=[fictionalgeoqa.jsonl path] --out=[output predicted answers filepath]
```
