"""The repetition-code sweep of quanterall.com (#sweeps), reproduced with public code only, and compared with the
published data (site/data/repetition-sweep.csv).

The published sweep (tools/repetition_sweep.py) ran on the TechUPGRADE simulation engine, which is not public. This
script recomputes every column independently:
  - logical_exact: every one of the 2^d bit-flip patterns of a distance-d repetition code, weighted by its probability
    and decoded by the minimum-weight correction (for independent flips on a repetition code: the majority vote);
  - logical_closed_form: sum over k > d/2 of C(d, k) p^k (1 - p)^(d - k);
  - the Monte Carlo columns: the same seed and draws as the published run (numpy's default_rng(20261029), one binomial
    draw of 100,000 shots per point, in the same order).
It prints the largest differences from the published CSV and exits 1 if any column disagrees.
    python3 tools/repetition_sweep_check.py site/data/repetition-sweep.csv
(needs Python 3.9+ and NumPy; the published run used NumPy 1.26)
"""
import csv
import itertools
import math
import sys

import numpy as np

SEED, SHOTS, DISTANCES = 20261029, 100_000, (3, 5, 7, 9)


def logical_exact(d, p):
    """sum over all 2^d flip patterns of P(pattern), for the patterns the decoder gets wrong"""
    fail = 0.0
    for flips in itertools.product((0, 1), repeat=d):
        w = sum(flips)
        if w > d // 2:                     # the majority vote (the minimum-weight correction) decodes to the wrong word
            fail += p ** w * (1 - p) ** (d - w)
    return fail


def closed_form(d, p):
    return sum(math.comb(d, k) * p ** k * (1 - p) ** (d - k) for k in range(d // 2 + 1, d + 1))


def main(path):
    published = list(csv.DictReader(open(path)))
    ps = np.geomspace(0.001, 0.5, 22)
    rng = np.random.default_rng(SEED)
    rows = []
    for d in DISTANCES:
        for p in ps:
            ex, cf = logical_exact(d, p), closed_form(d, p)
            fails = int(rng.binomial(SHOTS, ex))
            rows.append((d, 2 * d - 1, p, ex, cf, fails))
    if len(rows) != len(published):
        print(f"FAIL: {len(rows)} rows here, {len(published)} published")
        return 1
    worst_ex = worst_cf = 0.0
    mc_diff = 0
    for (d, q, p, ex, cf, fails), pub in zip(rows, published):
        assert int(pub["distance"]) == d and int(pub["qubits"]) == q and abs(float(pub["p"]) / p - 1) < 1e-5
        # the published values carry 11 significant digits
        worst_ex = max(worst_ex, abs(ex - float(pub["logical_exact"])) / max(ex, 1e-300))
        worst_cf = max(worst_cf, abs(cf - float(pub["logical_closed_form"])) / max(cf, 1e-300))
        mc_diff += fails != int(pub["mc_failures"])
    print(f"{len(rows)} rows (distances {', '.join(map(str, DISTANCES))}, 22 values of p each)")
    print(f"logical_exact:       largest relative difference from the published value {worst_ex:.1e}")
    print(f"logical_closed_form: largest relative difference from the published value {worst_cf:.1e}")
    print(f"Monte Carlo failures: {len(rows) - mc_diff} of {len(rows)} equal to the published counts")
    ok = worst_ex < 1e-9 and worst_cf < 1e-9 and mc_diff == 0
    print("PASS: the published sweep is reproduced" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "site/data/repetition-sweep.csv"))
