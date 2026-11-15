# Reproducing the benchmarks against Qiskit Aer and Stim

On quanterall.com ([#benchmarks](https://quanterall.com/#benchmarks)): two timing tables.

- **State vector:** random circuits of 8 to 20 qubits (n layers of Rx on every qubit and CNOTs on neighbours,
  1,024 shots), on QxSim (our state-vector simulator on the Erlang virtual machine), on our NumPy engine
  ([qsim/](../qsim)) and on Qiskit Aer.
- **Stabilizer sampling:** 100,000 shots of repetition and surface-code memories, on QxStab (our stabilizer simulator
  on the BEAM, 20 schedulers) and on Stim, on one core and over 20 processes.

The published data: [site/data/benchmarks.csv](../site/data/benchmarks.csv) and `benchmarks.json`.

## What you need

The same as for the noise study ([noise.md](noise.md)): Docker, git, the service image and the tools image. The
timings mean most on a machine like the published one (20 cores; the BEAM given 1 GB of memory, as in production).

## Steps

The benchmarks come out of the same pipeline as the noise study. From a clone of the repository, after
`git checkout 3c31766` to **reproduce the published run** (leave it out to **run the latest** version):

```sh
docker build -t quanterall-quantum quanterall_quantum
docker build -t qcs-tools -f tools/Dockerfile .
mkdir -p bench out
for job in study sv stab; do
  docker run --rm --memory=1g -v "$PWD/tools":/tools:ro -v "$PWD/bench":/bench -w /app quanterall-quantum \
    mix run --no-start /tools/beam_runs.exs $job /bench
done
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e OPENBLAS_NUM_THREADS=1 -v "$PWD":/repo:ro \
  -v "$PWD/bench":/bench -v "$PWD/out":/out -w /repo qcs-tools python tools/stab_study.py /bench /out
```

- `sv` times QxSim from 8 to 19 qubits. It stops after a run over 60 seconds, and does not try 20 qubits, which do
  not fit the service's 1 GB; the table says so.
- `stab` times QxStab on the five memories; `stab_study.py` times Stim (one core, then 20 worker processes), the
  NumPy engine and Qiskit Aer on the same circuits, and writes `out/benchmarks.csv` and `out/benchmarks.json`.
- `OPENBLAS_NUM_THREADS=1` keeps NumPy on one thread, as in the published run.
- The whole pipeline takes about 7 minutes after the images are built.

## What you should see

`stab_study.py` prints one line per state-vector size, for example

```
{'qubits': 18, 'gates': 630, 'runs': 3, 'beam_median_ms': 17546.4, ..., 'numpy_median_ms': 229.98, ..., 'aer_median_ms': 138.58, ...}
{'qubits': 20, 'gates': 780, 'runs': 3, 'beam_median_ms': None, ..., 'beam_note': 'out of memory in the 1 GB container', ...}
done: stim 1.16.0, pymatching 2.4.0
```

## How to compare with the published results

These are timings: they depend on the machine and on what else runs on it. Compare the shape, not the milliseconds:

- **State vector** (published): up to 7 qubits QxSim is faster than Qiskit Aer (whose fixed cost dominates there);
  from 8 qubits it is the slowest of the three, and its gap to our NumPy engine grows with the size, from 1.8 times
  at 8 qubits to about 33 times at 18 and 19. Aer is the fastest from 17 qubits up. The BEAM does not reach 20
  qubits in 1 GB.
- **Stabilizer sampling** (published, QxStab on 20 schedulers against Stim on one core): level on the repetition
  code (14.0 ms each), ahead at surface-code distance 5 (20.4 against 22.6 ms), 1.3 to 2.7 times behind on the other
  surface codes; against Stim over 20 processes, 27 times behind at distance 9.
- Your machine's numbers will differ, sometimes by more than a factor of 2 for the BEAM (it shares the cores with
  everything else). Where two programs are within a factor of about 1.5 of each other, they can change places from
  run to run; the larger differences are what should hold.
- The published table is in `site/data/benchmarks.csv`; `diff` it against `out/benchmarks.csv` to see your
  machine's numbers beside them. The columns `qubits` and `gates` are equal; only the times differ.

## Tested

From the repository on 15 November 2026, aarch64 (Cortex-X925 and Cortex-A725): the pipeline ran in about 7 minutes;
the circuits (qubits, gates) equal to the published ones. On a machine busy with other
work, the BEAM's state-vector times were up to 2.7 times the published ones, and two close pairs changed places (at
10 qubits QxSim and Aer, at 17 Aer and NumPy); QxStab was 1.4 to 3.3 times Stim's one-core time (published: 0.9 to
2.7), so its lead at distance 5 did not hold in this run. The large gaps held: QxSim far behind NumPy and Aer from
12 qubits, QxStab far behind Stim over 20 processes.
