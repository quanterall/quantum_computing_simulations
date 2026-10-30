# Credits

## qx — Craig Richards

The quantum demos on quanterall.com began on **qx**, an Elixir quantum computing simulator by Craig Richards:
https://github.com/richarc/qx — licensed under the Apache License 2.0. We started from version 0.5.1 (released
7 March 2026), and earlier versions of our README and demo page called it "qx_sim 0.5.1".

The simulator that runs now, `QxSim` in `quanterall_quantum/lib/qx_demo.ex`, is a separate implementation: a
state vector kept as an Erlang tuple of complex pairs, gates applied in one pass over it, sampling by binary
search. On 29 October 2026 we compared our code with every qx release, by a published method — `tools/qx_compare.py`,
result in `site/data/qx-comparison.json`: (1) every qx release depends on Nx, and our service does not depend on
qx; (2) none of 19 distinctive names of ours occurs anywhere in qx's history (207 commits, lib/ and examples/);
(3) cut into 8-token sequences (comments, docs and strings removed), at most 2 % of any of our files occurs in any
qx release, while consecutive qx releases share 92–100 % (the positive control: the measure does find shared code).
We found no shared code. The idea of a quantum
simulator in Elixir, and the start of this project, came from qx, and we credit it for that.

If qx code is ever copied into this repository, it keeps its Apache-2.0 licence and copyright notice, and this file
says which files.

## Stim, PyMatching, Qiskit Aer — used to check our work

The gate-level noise study on quanterall.com checks our stabilizer simulator, QxStab, against **Stim**
(Craig Gidney, https://github.com/quantumlib/Stim, Apache-2.0): our circuits are written in Stim's format, Stim
derives their detector error models and code distances, and samples the same circuits for comparison. Both
samples are decoded by **PyMatching** (Oscar Higgott, https://github.com/oscarhiggott/PyMatching, Apache-2.0). The
state-vector benchmark compares with **Qiskit Aer** (https://github.com/Qiskit/qiskit-aer, Apache-2.0). None of
their code is in this repository; they run in `tools/` only (`.venv-tools`, see `tools/stab_study.py`). QxStab
follows the published methods: the tableau of Aaronson and Gottesman (Phys. Rev. A 70, 052328, 2004) and the
Pauli-frame sampling described for Stim (Gidney, Quantum 5, 497, 2021).
