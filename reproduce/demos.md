# Reproducing the quantum demos

On quanterall.com ([#demos](https://quanterall.com/#demos)): ten circuits (Bell states, teleportation, Grover, a small
quantum classifier, Shor's algorithm for 15, error-correcting codes and others) simulated live on the Erlang virtual
machine by the service in [quanterall_quantum/](../quanterall_quantum). This page runs the same service on your
machine and checks every endpoint.

## What you need

- Docker, git, curl and bash; about 2 GB of disk for the service's image (Elixir 1.18 and its dependencies, at the
  versions in [mix.lock](../quanterall_quantum/mix.lock)).
- Any Linux or macOS machine.

## Steps

```sh
git clone https://github.com/quanterall/quantum_computing_simulations.git
cd quantum_computing_simulations
git checkout 3c31766
docker build -t quanterall-quantum quanterall_quantum
docker run -d --rm --name qcs-demo -p 127.0.0.1:8080:8080 quanterall-quantum
curl 'http://127.0.0.1:8080/quantum/bell?shots=1000&seed=1'
bash quanterall_quantum/test_all_quantum_demos.sh http://127.0.0.1:8080
docker stop qcs-demo
```

`git checkout 3c31766` **reproduces the published run**: the commit this page was tested at. To **run the latest**
version, leave it out; the results must still be the published ones.

The image builds in about 2 minutes. The service starts in a few seconds (`curl http://127.0.0.1:8080/health` answers
when it is ready).

## What you should see

The `curl` line returns JSON: the counts of the 1,000 shots of a Bell pair (about half `00`, half `11`, never `01` or
`10`) and the circuit in OpenQASM. The same seed gives the same counts every time.

The test script checks every endpoint, its answer's shape and its physics (for example that the Bell pair only gives
`00` and `11`, that Shor's algorithm finds 3 × 5, that each error-correcting code corrects its error):

```
Quanterall Quantum Demo Suite — Integration Test
Target: http://127.0.0.1:8080
────────────────────────────────────────────────

Algorithms, protocols and QML:
  ✓ Bell State
  ✓ Teleportation
  ✓ Grover Search
  ...
Results: 39 passed, 0 failed
```

Run it with **bash**: under `sh` the counts at the end come out as 0.

## How to compare with the published results

- Every check passes: `39 passed, 0 failed`.
- The same request with the same `seed` gives the same counts as quanterall.com's service, for example
  `curl 'https://quanterall.com/quantum/bell?shots=1024&seed=42'` against your
  `curl 'http://127.0.0.1:8080/quantum/bell?shots=1024&seed=42'`. The page's "Reproduce this" box lists more.

## Tested

From the repository on 15 November 2026, aarch64: the image built in 1 minute 51 seconds, 39 of 39 checks passed in
11 seconds.
