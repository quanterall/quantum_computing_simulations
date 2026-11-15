# Reproducing the surface code under circuit noise

On quanterall.com ([#noise](https://quanterall.com/#noise)): surface-code memories of distance 3, 5 and 7 under
circuit-level noise (p from 0.1 % to 1.2 %), sampled by QxStab, our stabilizer simulator on the Erlang virtual
machine, and by Stim; both decoded by PyMatching; the logical error per round compared. The published data:
[site/data/surface-circuit-noise.csv](../site/data/surface-circuit-noise.csv) and its chart.

## What you need

- Docker, git, and about 4 GB of disk for two images: the service ([quanterall_quantum/](../quanterall_quantum),
  Elixir, with QxStab) and the tools ([tools/Dockerfile](../tools/Dockerfile): Stim 1.16.0, PyMatching 2.4.0, Qiskit
  2.6.1, Qiskit Aer 0.18.2 and NumPy, at the versions of the published results).
- Any Linux machine. The published run was on aarch64.

## Steps

```sh
git clone https://github.com/quanterall/quantum_computing_simulations.git
cd quantum_computing_simulations
docker build -t quanterall-quantum quanterall_quantum
docker build -t qcs-tools -f tools/Dockerfile .
mkdir -p bench out
# 1. on the BEAM: the noisy memories sampled by QxStab (study), state-vector timings (sv), stabilizer timings (stab)
for job in study sv stab; do
  docker run --rm --memory=1g -v "$PWD/tools":/tools:ro -v "$PWD/bench":/bench -w /app quanterall-quantum \
    mix run --no-start /tools/beam_runs.exs $job /bench
done
# 2. the same circuits sampled by Stim, both samples decoded by PyMatching, and the benchmarks' other programs
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/repo:ro -v "$PWD/bench":/bench -v "$PWD/out":/out \
  -w /repo qcs-tools python tools/stab_study.py /bench /out
```

- The image builds take about 2 minutes (the service) and 5 (the tools: on aarch64, Stim and PyMatching compile from
  source).
- Step 1 takes about 4.5 minutes; step 2 about 2. `--memory=1g` gives the BEAM the published run's limit (it matters
  only for the benchmarks' timings).
- The same pipeline writes the benchmarks' data too ([benchmarks.md](benchmarks.md)).

## What you should see

Step 1 prints one line per memory (`surface_d3_p0.001: 100000 shots in ... ms`, 21 of them) and the timings. Step 2
writes to `out/`: `surface-circuit-noise.csv` and `.svg`, `benchmarks.csv` and `.json`, and ends with
`done: stim 1.16.0, pymatching 2.4.0`. On the way it prints the statistical comparison of the two samplers'
detection events (`detector_max_abs_z`, ...).

## How to compare with the published results

```sh
diff site/data/surface-circuit-noise.csv out/surface-circuit-noise.csv
```

- **Every count must be equal:** QxStab's and Stim's failures, logical error rates, per-round rates and detection
  rates, on all 21 rows. Both samplers are seeded, and decoding is deterministic. Only the two timing columns
  (`beam_ms`, `stim_ms`) differ.
- Stim's samples are the same for the same seed, Stim version and machine architecture; on another architecture its
  counts may differ within their statistical errors (the `_err` columns), while QxStab's stay equal.
- The chart `out/surface-circuit-noise.svg` is drawn from the counts and is then identical too.

## Exact tests

The page also reports exact tests of QxStab against Stim ([tools/stim_exact.py](../tools/stim_exact.py)): exact
two-sided binomial tests of the two samplers' counts, detector by detector. After step 1 (8 seconds):

```sh
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/repo:ro -v "$PWD/bench":/bench -v "$PWD/out":/out \
  -w /repo qcs-tools python tools/stim_exact.py /bench /out
diff site/data/stim-agreement.json out/stim-agreement.json
```

The result, `out/stim-agreement.json`, must be identical to the published one.

## Tested

From the repository on 15 November 2026, aarch64 (Cortex-X925 and Cortex-A725), with the images above: all 21 rows
equal to the published ones in every column but the two timings; the chart and the exact tests identical; steps 1 and 2 in 4.5 and 2
minutes.
