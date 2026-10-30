"""Our stabilizer sampler against Stim with exact tests (quanterall.com, Results · gate-level noise).

For each detector at every one of the study's 21 points, and each pair of detectors (joint firing) at d = 5,
p = 0.6 %: the counts from our sample (k1) and from Stim's (k2), equal shot numbers. Under agreement k1 given
k1 + k2 is Binomial(k1 + k2, 1/2): an exact two-sided binomial test, valid for small counts. The tests are not
independent (detectors share errors), so the share of small p-values is a guide, not a calibrated test; the
smallest p-value is also given Bonferroni-corrected (times the number of tests).

    ~/quanterall/.venv-tools/bin/python tools/stim_exact.py BENCH_DIR site/data
"""
import json
import sys
from pathlib import Path

import numpy as np
import stim
from scipy.stats import binomtest

bench, out = Path(sys.argv[1]), Path(sys.argv[2])


def ours(name, c, shots):
    raw = np.fromfile(bench / f"{name}.bin", dtype=np.uint8).reshape(c.num_detectors + c.num_observables, -1)
    return np.unpackbits(raw, axis=1, bitorder="little")[:, :shots].T.astype(bool)[:, :c.num_detectors]


def exact(k1, k2):
    return np.array([binomtest(int(a), int(a + b), 0.5).pvalue if a + b > 0 else 1.0 for a, b in zip(k1, k2)])


def ours_obs(name, c, shots):
    raw = np.fromfile(bench / f"{name}.bin", dtype=np.uint8).reshape(c.num_detectors + c.num_observables, -1)
    return np.unpackbits(raw, axis=1, bitorder="little")[:, :shots].T.astype(bool)[:, c.num_detectors:]


def summarise(ps):
    ps = np.asarray(ps)
    return {"tests": int(ps.size), "share_p_below_0.0027": float((ps < 0.0027).mean()), "expected_share": 0.0027,
            "smallest_p": float(ps.min()), "smallest_p_bonferroni": float(min(1.0, ps.min() * ps.size))}


single, pairs, det_obs, logical = [], [], [], []
power = None
for run in json.load(open(bench / "study.json")):
    text = open(bench / f"{run['name']}.stim").read()
    c = stim.Circuit(text)
    a, ao = ours(run["name"], c, run["shots"]), ours_obs(run["name"], c, run["shots"])
    b, bo = c.compile_detector_sampler(seed=run["seed"]).sample(run["shots"], separate_observables=True)
    single.append(exact(a.sum(0), b.sum(0)))
    # each detector together with the logical observable: the correlation that sets the logical error rate
    det_obs.append(exact((a & ao).sum(0), (b & bo).sum(0)))
    if run["d"] == 5 and run["p"] == 0.006:
        iu = np.triu_indices(c.num_detectors, 1)
        joint = lambda d: sum(((d[k:k + 5000, :, None] & d[k:k + 5000, None, :])[:, iu[0], iu[1]]).sum(0) for k in range(0, len(d), 5000))
        pairs = exact(joint(a), joint(b))
        # power: the same tests between Stim at p = 0.6 % and Stim at p = 0.65 % must find the difference
        c2 = stim.Circuit(text.replace("(0.006)", "(0.0065)"))
        b2, _ = c2.compile_detector_sampler(seed=99).sample(run["shots"], separate_observables=True)
        power = summarise(exact(b2.sum(0), b.sum(0)))
single, det_obs = np.concatenate(single), np.concatenate(det_obs)

# the logical error rates themselves, decoded, from the study's table (tools/stab_study.py)
import csv  # noqa: E402
for row in csv.DictReader(open(out / "surface-circuit-noise.csv")):
    k1, k2 = int(row["beam_failures"]), int(row["stim_failures"])
    logical.append({"d": int(row["distance"]), "p": float(row["p"]), "beam_failures": k1, "stim_failures": k2,
                    "p_value": binomtest(k1, k1 + k2, 0.5).pvalue if k1 + k2 else 1.0})

res = {"detectors": summarise(single), "detector_and_observable": summarise(det_obs), "pairs": summarise(pairs),
       "logical_error_rates": {"tests": len(logical), "smallest_p": min(r["p_value"] for r in logical),
                               "smallest_p_bonferroni": min(1.0, min(r["p_value"] for r in logical) * len(logical)),
                               "rows": logical},
       "power_check": {"what": "Stim at p = 0.65 % against Stim at p = 0.6 %, d = 5, per detector: the tests must find this", **power},
       "note": ("exact two-sided binomial tests of our count against Stim's, conditional on their sum; the tests are "
                "correlated (detectors share error mechanisms), so the shares are a guide rather than a calibrated test")}
json.dump(res, open(out / "stim-agreement.json", "w"), indent=1)
print(json.dumps({k: (v if k != "logical_error_rates" else {kk: vv for kk, vv in v.items() if kk != "rows"}) for k, v in res.items()}, indent=1))
