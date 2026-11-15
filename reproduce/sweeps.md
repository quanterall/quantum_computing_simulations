# Reproducing the repetition-code sweep

On quanterall.com ([#sweeps](https://quanterall.com/#sweeps)): the logical error rate of repetition codes of distance
3, 5, 7 and 9 against the physical error rate, 22 values of p for each, with a seeded Monte Carlo beside the exact
values. The published data: [site/data/repetition-sweep.csv](../site/data/repetition-sweep.csv) (88 rows).

## What was computed, and what you can recompute

The published sweep ([tools/repetition_sweep.py](../tools/repetition_sweep.py)) ran on the simulation engine we built
for the EU project TechUPGRADE, which is not public. Each row has:

- **logical_exact:** every one of the 2^d bit-flip patterns sent through the code's syndrome circuit on that engine,
  weighted by its probability, and decoded;
- **logical_closed_form:** the binomial sum over k > d/2 of C(d, k) p^k (1 − p)^(d − k);
- **the Monte Carlo columns:** one seeded binomial draw of 100,000 shots per point.

[tools/repetition_sweep_check.py](../tools/repetition_sweep_check.py) recomputes every column with public code only
(Python and NumPy): the exhaustive sum over the 2^d patterns with the minimum-weight (majority) decoder, the closed
form, and the Monte Carlo draws with the same seed in the same order. It compares each with the published CSV.

## What you need

Docker (or Python 3.9+ with NumPy) and git. Nothing else; it takes seconds.

## Steps

```sh
git clone https://github.com/quanterall/quantum_computing_simulations.git
cd quantum_computing_simulations
git checkout 3c31766
docker run --rm -v "$PWD":/w -w /w python:3.12-slim \
  sh -c 'pip install -q numpy==1.26.4 && python tools/repetition_sweep_check.py site/data/repetition-sweep.csv'
```

`git checkout 3c31766` **reproduces the published run**: the commit this page was tested at. To **run the latest**
version, leave it out; the results must still be the published ones.

(With Python and NumPy installed: `python3 tools/repetition_sweep_check.py site/data/repetition-sweep.csv`.)

## What you should see

```
88 rows (distances 3, 5, 7, 9, 22 values of p each)
logical_exact:       largest relative difference from the published value 3.9e-11
logical_closed_form: largest relative difference from the published value 3.9e-11
Monte Carlo failures: 88 of 88 equal to the published counts
PASS: the published sweep is reproduced
```

- The relative differences of about 4e-11 are the rounding of the published values, which carry 11 significant
  digits.
- The Monte Carlo counts are equal exactly: the same seed, the same draws. NumPy versions with the same `Generator`
  stream (1.17 and later) give the same counts.

## The chart

`tools/repetition_sweep.py` also drew [site/data/repetition-sweep.svg](../site/data/repetition-sweep.svg) from the same
rows; it needs the TechUPGRADE engine for its first column and is kept as the record of the published run.

## Tested

From a fresh clone on 15 November 2026, in `python:3.12-slim` with NumPy 1.26.4: PASS, all 88 rows.
