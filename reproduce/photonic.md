# Reproducing the photonic schemes

On quanterall.com ([#photonic](https://quanterall.com/#photonic)): photonic schemes from the literature, simulated
and checked. The Hong–Ou–Mandel dip, the KLM gates, fusion, growing a cluster chain, measurement-based computation,
detectors and loss, the 2.7 % loss figure, erasure decoding, and one optical mode (squeezing, GKP and cat states).
The published data are the files in [site/data/](../site/data) that each step below names.

The simulators are written twice, from one conventions file ([docs/optics_conventions.md](../docs/optics_conventions.md)):
in Python ([qsim/qsim/](../qsim/qsim): `fock.py`, `cv.py`) and in Elixir on the BEAM
([quanterall_quantum/](../quanterall_quantum): QxFock, QxCV). The checks run the Python side; the cross-checks compare
the two.

## What you need

- Docker and git; the tools image ([tools/Dockerfile](../tools/Dockerfile), Python with NumPy, SciPy, Qiskit,
  Perceval, Stim and PyMatching at the versions of the published results).
- For the cross-checks and the erasure study, which also drive the BEAM service through `docker`: the service image
  and a Python environment on the host (Python 3.12) from [tools/requirements.txt](../tools/requirements.txt).

```sh
git clone https://github.com/quanterall/quantum_computing_simulations.git
cd quantum_computing_simulations
git checkout 3c31766
docker build -t qcs-tools -f tools/Dockerfile .          # about 5 minutes
mkdir -p out
```

`git checkout 3c31766` **reproduces the published run**: the commit this page was tested at. To **run the latest**
version, leave it out; the results must still be the published ones.

## 1. The checks (Python, in the tools image)

```sh
for f in ns_check klm_checks xcheck_linear fusion_checks cluster_growth mbqc loss_checks fbqc_loss cv_checks; do
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/repo:ro -v "$PWD/out":/out -w /repo qcs-tools \
    python tools/optics/$f.py /out
done
```

| script | what it checks | output | time |
| --- | --- | --- | --- |
| `ns_check.py` | KLM's NS-gate unitary (Nature 409, 46, Eq. 9): probability 1/4, and which reading of the review's angles works | printed | 1 s |
| `klm_checks.py` | dual-rail gates against the qubit engine, the HOM dip, the NS gate and CZ, with checks that must fail when a parameter is off | `klm-checks.json`, `hom-dip.csv` | 1 s |
| `xcheck_linear.py` | the Fock module's two routes (creation operators, Ryser permanents) against each other and against Perceval | `fock-vs-perceval.json` | 2 s |
| `fusion_checks.py` | type-I and type-II fusion (Browne and Rudolph 2005) photon by photon | `fusion-checks.json` | 1 s |
| `cluster_growth.py` | growing a linear cluster by fusion: Bell pairs per added qubit | `cluster-growth.json` | 19 s |
| `mbqc.py` | measurement-based computation on small clusters, against Qiskit | `mbqc-checks.json` | 1 s |
| `loss_checks.py` | detectors of efficiency η, and a lost photon as an erasure | `loss-checks.json`, `heralded-vs-efficiency.csv` | 1 s |
| `fbqc_loss.py` | the review's 2.7 % loss figure traced to its source | printed | 1 s |
| `cv_checks.py` | one mode: squeezing (Eqs. 22–26), GKP states (Eqs. 38–39), cat states, the heralded comb | `cv-checks.json`, `gkp-*.csv`, `cv-wigner.json`, `heralded-comb*.csv` | 2.5 min |

**Compare:** `for f in out/*; do cmp "$f" "site/data/$(basename "$f")" && echo "same $f"; done`. Every file must be
identical, except `cv-checks.json` and `gkp-readout.csv`, which differ only in the last digits of floating-point
sums (about 1e-14 relative). One value, `p_exact` at σ = 0.1 in `gkp-readout.csv`, is about 1e-12 and differs in
its third digit: it is a numerical integral at the noise floor of the integrator (the closed form there is 4.6e-13),
and the page does not use it.

**Fig. 5 of the review.** `cv_checks.py` also reads the review's Fig. 5 (Romero and Milburn,
[arXiv:2404.03367](https://arxiv.org/abs/2404.03367), version 2) from the vector curve in its PDF, when the PDF is
given as a second argument. Without it, those two checks (`fig5_reading`, `fig5_digitised`) are skipped and
`cv-checks.json` has 21 entries instead of 23. With it:

```sh
curl -sL -o 2404.03367v2.pdf https://arxiv.org/pdf/2404.03367v2
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/repo:ro -v "$PWD/out":/out -w /repo qcs-tools \
  python tools/optics/cv_checks.py /out /repo/2404.03367v2.pdf
```

`out/gkp-fig5-digitised.json` must then be identical to the published file, and `cv-checks.json` has all 23 entries.

## 2. The cross-checks, Python against the BEAM

```sh
docker build -t quanterall-quantum quanterall_quantum
docker run -d --rm --name qcs-check quanterall-quantum
python3 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt    # on aarch64 Stim compiles: minutes
QX_CONTAINER=qcs-check .venv/bin/python tools/optics/xcheck_beam.py
QX_CONTAINER=qcs-check .venv/bin/python tools/optics/xcheck_cv.py out
```

`QX_CONTAINER` names the running service container the scripts use (`docker cp` and `docker exec`).

- `xcheck_beam.py` sends 20 random interferometers to QxFock and compares every output amplitude with `fock.py`:
  it prints `{"cases":20,"max_diff":1.5e-15}` (any value near 1e-15 is the rounding of the two programs).
- `xcheck_cv.py` does the same for QxCV against `cv.py` and writes `out/cv-xcheck.json`; it equals the published
  file up to the last digits (about 1e-15).

## 3. Erasure decoding

With the service container from step 2 still running:

```sh
mkdir -p cl work
docker run --rm -v "$PWD/tools":/tools:ro -v "$PWD/cl":/cl -w /app quanterall-quantum \
  mix run --no-start /tools/erasure/circuit_runs.exs /cl 20000
QX_CONTAINER=qcs-check .venv/bin/python tools/erasure/circuit.py cl out
QX_CONTAINER=qcs-check .venv/bin/python tools/erasure/sweep.py work out
docker stop qcs-check
```

- `circuit_runs.exs`: surface-code memories with heralded erasure, sampled by QxStab on the BEAM.
- `circuit.py`: those samples decoded by QxDecode (our union-find decoder, on the BEAM) with the heralds, and by
  PyMatching with and without them; `out/erasure-circuit.csv`. About 3.5 minutes.
- `sweep.py`: the code-capacity sweeps, QxDecode against PyMatching on the same samples; `out/erasure-pure.csv` and
  `out/erasure-mixed.csv`. About 3 minutes.

**Compare:** `cmp` the three files with `site/data/`: they must be identical. Every sample is seeded, on the BEAM
and in Python, and both decoders are deterministic.

## Tested

From the repository on 15 November 2026, aarch64 (Cortex-X925 and Cortex-A725): every check script ran, with the
outputs above; the Fig. 5 checks (with the arXiv PDF, version 2) identical to the published
ones; the cross-checks agreed to 1.6e-15; the three erasure files identical to the published ones (100 rows), in
about 7 minutes.
