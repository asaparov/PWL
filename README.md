# Probabilistic Worldbuilding from Language

If you use this data or code in your research, please cite:

```bibtex
@article{10.1162/tacl_a_00463,
  author = {Abulhair Saparov and Tom M. Mitchell},
  title = {Towards General Natural Language Understanding with Probabilistic Worldbuilding},
  journal = {Transactions of the Association for Computational Linguistics},
  volume = {10},
  pages = {325-342},
  year = {2022},
  month = {04},
  issn = {2307-387X},
  doi = {10.1162/tacl_a_00463},
  url = {https://doi.org/10.1162/tacl\_a\_00463}
}
```

This code depends on OpenSSL as well as the following repositories: [core](https://github.com/asaparov/core), [math](https://github.com/asaparov/math), [hdp](https://github.com/asaparov/hdp), and [grammar](https://github.com/asaparov/grammar).

To use this code:
1. Download the dependencies into a single directory.
2. Download this repository and run `make executive_test CPPFLAGS+="-I[deps_directory]"` where `deps_directory` is the folder containing the dependency folders `core`, `math`, `hdp`, and `grammar`.
3. Run `./executive_test`.

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
