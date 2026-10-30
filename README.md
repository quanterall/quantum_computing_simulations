# Quantum computing simulations

The code behind [quanterall.com](https://quanterall.com): quantum algorithms, quantum error correction and photonic
schemes from the literature, reproduced in simulation on the Erlang VM (BEAM) and in Python, each result checked
against an independent tool. This is simulation, not quantum hardware.

## What is here

| Path | What it is |
|---|---|
| `quanterall_quantum/` | The live service: Elixir and Erlang on the BEAM, HTTP API on [Nova](https://github.com/novaframework/nova). `QxSim` (state vectors on Erlang tuples), `QxStab` (stabilizer tableau and a bit-packed Pauli-frame sampler with gate-level noise), `QxDecode` (a union-find decoder with erasures), `QxFock` (linear optics in Fock space), `QxCV` (one optical mode: GKP and cat states), `QxProcesses` (entanglement and teleportation as linked processes). See its README for the endpoints. |
| `qsim/` | The Python engine: the quantum module of the TechUPGRADE digital-twin platform's simulation engine (`engine.py` state vectors, `fock.py` Fock space, `cv.py` continuous variables), with a FastAPI app and tests. |
| `tools/` | The studies, checks and benchmarks: against Stim, PyMatching, Qiskit Aer and Perceval; the decoder comparison and load test; the photonic and erasure checks. Each script says how it is run. |
| `site/` | The website and the data and charts it publishes (`site/data/`). |
| `docs/` | The photonic-schemes plan and the optics conventions shared by the Elixir and Python code. |
| `NOTES.md` | Errata found in the sources while reproducing them, with how each was checked. |
| `CREDITS.md` | Where the work began and the tools it is checked against. |
| `deploy/` | An example deployment (nginx and the services). |

## Running it

The BEAM service:

```sh
docker build -t quanterall-quantum quanterall_quantum
docker run --rm -p 8080:8080 quanterall-quantum
curl 'http://localhost:8080/quantum/bell?shots=1000&seed=1'
quanterall_quantum/test_all_quantum_demos.sh http://localhost:8080     # every endpoint, with checks
```

The Python engine (tests run in the build):

```sh
docker build --target test qsim                 # the tests
docker build --target run -t quanterall-qsim qsim && docker run --rm -p 8000:8000 quanterall-qsim
```

The site, the service and nginx together: `docker compose -f deploy/docker-compose.example.yml up -d`, then
`http://localhost:8100`.

The tools need a Python environment with `numpy scipy stim pymatching qiskit qiskit-aer perceval-quandela`
(the published results used stim 1.16, pymatching 2.4, qiskit 2.6.1, qiskit-aer 0.18.2, perceval-quandela 1.3.1).
Scripts that drive the BEAM run their Elixir part inside the service's container.

## How the results are checked

Every published number comes with the check behind it: an independent implementation (Stim, PyMatching, Qiskit
Aer, Perceval, or the other of our two engines), a fixed seed, the data next to the chart, and a power check — the
same check with a deliberate error must fail. Where a source turned out to be wrong, `NOTES.md` says how and how
it was found.

The code was written with an AI coding assistant (Claude); the checks above are how its results are verified.

## License

Apache License 2.0 — see `LICENSE` and `NOTICE`. Copyright 2026 Quanterall Ltd.
