"""Our union-find decoder (lib/qx_decode.ex, run by tools/decode_runs.exs) against PyMatching, on the same shots.

For each point of the threshold study (d = 3, 5, 7; p = 0.1, 0.4, 0.6, 0.8 %): the logical error rate of each
decoder on the same 100,000 syndromes (our sampler's), and McNemar's exact test on the shots where they
disagree; time per syndrome on one core (PyMatching in one thread, decode_batch; ours in one BEAM process), and
our throughput on all 20 cores. The decoding graph of both comes from Stim's detector error model.

    ~/quanterall/.venv-tools/bin/python tools/decoder_compare.py BENCH_DIR UF_DIR site/data
"""
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pymatching
import stim
from scipy.stats import binomtest

bench, ufdir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
uf = {r["name"]: r for r in json.load(open(ufdir / "uf.json"))}
rows = []
for name, r in uf.items():
    c = stim.Circuit(open(bench / f"{name}.stim").read())
    D, S = c.num_detectors, r["shots"]
    raw = np.fromfile(bench / f"{name}.bin", dtype=np.uint8).reshape(D + 1, -1)
    bits = np.unpackbits(raw, axis=1, bitorder="little")[:, :S].T.astype(bool)
    det, obs = bits[:, :D], bits[:, D]
    m = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
    m.decode_batch(det[:2000])                                              # warm-up
    t0 = time.perf_counter()
    pm = m.decode_batch(det)[:, 0].astype(bool)
    t_pm = time.perf_counter() - t0
    ours = np.frombuffer(open(ufdir / f"{name}.pred", "rb").read(), dtype=np.uint8).astype(bool)
    pm_fail, uf_fail = pm != obs, ours != obs
    only_uf, only_pm = int(np.sum(uf_fail & ~pm_fail)), int(np.sum(pm_fail & ~uf_fail))
    p_mcnemar = binomtest(only_uf, only_uf + only_pm, 0.5).pvalue if only_uf + only_pm else 1.0
    pl, ul = pm_fail.mean(), uf_fail.mean()
    rows.append({"distance": r["d"], "p": r["p"], "shots": S,
                 "pymatching_logical": round(float(pl), 6), "pymatching_err": round(math.sqrt(pl * (1 - pl) / S), 6),
                 "uf_logical": round(float(ul), 6), "uf_err": round(math.sqrt(ul * (1 - ul) / S), 6),
                 "uf_over_pymatching": round(float(ul / pl), 3) if pl else None,
                 "only_uf_fails": only_uf, "only_pymatching_fails": only_pm, "mcnemar_p": p_mcnemar,
                 "pymatching_us_per_syndrome_1core": round(t_pm / S * 1e6, 3),
                 "uf_us_per_syndrome_1process": round(r["uf_one_process_per_syndrome_us"], 3),
                 "uf_syndromes_per_s_all_cores": round(S / r["uf_all_cores_s"]),
                 "pymatching_syndromes_per_s_1core": round(S / t_pm)})
    print(rows[-1])
rows.sort(key=lambda x: (x["distance"], x["p"]))
with open(out / "decoder-vs-pymatching.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
json.dump({"rows": rows, "versions": f"pymatching {pymatching.__version__}, stim {stim.__version__}",
           "note": "same 100,000 syndromes per point (sampled by QxStab); McNemar's exact test on the shots where the two decoders disagree; PyMatching single-threaded (decode_batch), ours in one BEAM process for the per-syndrome time and on 20 schedulers for the throughput"},
          open(out / "decoder-vs-pymatching.json", "w"), indent=1)
