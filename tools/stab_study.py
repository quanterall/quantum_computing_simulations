"""The gate-level noise study and the benchmarks published on quanterall.com (#noise, #benchmarks).

  1. On the BEAM (tools/beam_runs.exs): surface-code memories sampled by QxStab, our stabilizer simulator;
     state-vector timings of QxSim; stabilizer timings. Their files land in BENCH (default /tmp/claude-bench).
  2. Here: the same circuits (Stim's text, written by QxStab.to_stim) sampled by Stim; both samples decoded by
     PyMatching on the detector error model Stim derives from the circuit; the logical error per round; Stim's
     sampling time; the state-vector circuits on our NumPy engine and on Qiskit Aer.

    ~/quanterall/.venv-tools/bin/python tools/stab_study.py BENCH_DIR SITE_DATA_DIR
    ~/quanterall/.venv-tools/bin/python tools/stab_study.py --chart-only SITE_DATA_DIR    (redraw from the CSV)
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

CHART_ONLY = sys.argv[1] == "--chart-only"
bench, out = (None, Path(sys.argv[2])) if CHART_ONLY else (Path(sys.argv[1]), Path(sys.argv[2]))
out.mkdir(parents=True, exist_ok=True)
VERS = f"stim {stim.__version__}, pymatching {pymatching.__version__}"


def ours(name, c, shots):
    """our packed sample: detectors then observables, ceil(shots/8) bytes each, bit k = shot k"""
    raw = np.fromfile(bench / f"{name}.bin", dtype=np.uint8).reshape(c.num_detectors + c.num_observables, -1)
    bits = np.unpackbits(raw, axis=1, bitorder="little")[:, :shots].T.astype(bool)
    return bits[:, :c.num_detectors], bits[:, c.num_detectors:]


def per_round(p_total, rounds):
    return 1 - (1 - p_total) ** (1 / rounds) if p_total < 0.5 else float("nan")


# ── 1. the noise study ────────────────────────────────────────────────────────────────────────────────────
rows = []
for run in ([] if CHART_ONLY else json.load(open(bench / "study.json"))):
    c = stim.Circuit(open(bench / f"{run['name']}.stim").read())
    m = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
    d_ours, o_ours = ours(run["name"], c, run["shots"])
    t0 = time.perf_counter()
    d_stim, o_stim = c.compile_detector_sampler(seed=run["seed"]).sample(run["shots"], separate_observables=True)
    t_stim = time.perf_counter() - t0
    row = {"distance": run["d"], "p": run["p"], "rounds": run["rounds"], "qubits": run["qubits"], "shots": run["shots"],
           "graphlike_distance": len(c.shortest_graphlike_error())}
    for who, (dd, oo) in (("beam", (d_ours, o_ours)), ("stim", (d_stim, o_stim))):
        fails = int(np.sum(np.any(m.decode_batch(dd) != oo, axis=1)))
        pl = fails / run["shots"]
        row[f"{who}_failures"] = fails
        row[f"{who}_logical"] = pl
        row[f"{who}_logical_err"] = math.sqrt(max(pl * (1 - pl), 1e-12) / run["shots"])
        row[f"{who}_per_round"] = per_round(pl, run["rounds"])
        row[f"{who}_detection_rate"] = float(dd.mean())
    row["beam_ms"], row["stim_ms"] = run["beam_ms"], t_stim * 1000
    rows.append(row)
    print(f"d={run['d']} p={run['p']}: logical {row['beam_logical']:.5f} (ours) / {row['stim_logical']:.5f} (Stim)")

if CHART_ONLY:
    rows = [{k: (float(v) if k not in ("distance", "rounds", "qubits", "shots", "graphlike_distance", "beam_failures", "stim_failures") else int(v))
             for k, v in row.items()} for row in csv.DictReader(open(out / "surface-circuit-noise.csv"))]
else:
    with open(out / "surface-circuit-noise.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

# the chart: logical error per round against p, log-log, ours (lines, filled) and Stim (rings)
W, H, L, R, T, B = 760, 420, 74, 20, 20, 50
PX, PY = (0.001, 0.012), (1e-6, 0.3)
x = lambda p: L + math.log10(p / PX[0]) / math.log10(PX[1] / PX[0]) * (W - L - R)
y = lambda v: T + (math.log10(PY[1]) - math.log10(min(max(v, PY[0]), PY[1]))) / (math.log10(PY[1]) - math.log10(PY[0])) * (H - T - B)
col = {3: "#6ab0ff", 5: "#a8ff3e", 7: "#f472b6"}
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Surface code under circuit-level noise: '
       'logical error per round against the physical error rate, d = 3, 5, 7, sampled on the BEAM and in Stim" '
       'font-family="IBM Plex Mono, monospace" font-size="11">', f'<rect width="{W}" height="{H}" fill="#0f1726"/>']
for e in range(-6, 0):
    yy = y(10.0 ** e)
    svg.append(f'<line x1="{L}" y1="{yy:.1f}" x2="{W - R}" y2="{yy:.1f}" stroke="#1e2a40"/><text x="{L - 8}" y="{yy + 4:.1f}" fill="#7a879c" text-anchor="end">1e{e}</text>')
for p in (0.001, 0.002, 0.004, 0.008, 0.012):
    xx = x(p)
    svg.append(f'<line x1="{xx:.1f}" y1="{T}" x2="{xx:.1f}" y2="{H - B}" stroke="#1e2a40"/><text x="{xx:.1f}" y="{H - B + 16}" fill="#7a879c" text-anchor="middle">{p * 100:g}%</text>')
svg.append(f'<text x="{(L + W - R) / 2}" y="{H - 12}" fill="#9aa6bb" text-anchor="middle">physical error rate p (every gate, reset and measurement)</text>')
svg.append(f'<text x="16" y="{(T + H - B) / 2}" fill="#9aa6bb" text-anchor="middle" transform="rotate(-90 16 {(T + H - B) / 2})">logical error per round</text>')
lx, ly = W - R - 300, H - B - 96                                   # the legend, lower right, below the curves
for i, d in enumerate((3, 5, 7)):
    pts = [(r["p"], r["beam_per_round"], r["stim_per_round"], r["beam_failures"]) for r in rows if r["distance"] == d]
    shown = [(p, b, s) for p, b, s, f in pts if f >= 10]          # fewer than 10 failures: the estimate is not shown
    if shown:
        svg.append(f'<polyline points="{" ".join(f"{x(p):.1f},{y(b):.1f}" for p, b, _ in shown)}" fill="none" stroke="{col[d]}" stroke-width="2"/>')
    for p, b, s in shown:
        svg.append(f'<circle cx="{x(p):.1f}" cy="{y(b):.1f}" r="2.4" fill="{col[d]}"/><circle cx="{x(p):.1f}" cy="{y(s):.1f}" r="4.5" fill="none" stroke="{col[d]}"/>')
    svg.append(f'<line x1="{lx}" y1="{ly + i * 16}" x2="{lx + 20}" y2="{ly + i * 16}" stroke="{col[d]}" stroke-width="2"/>'
               f'<text x="{lx + 26}" y="{ly + 4 + i * 16}" fill="#c8d4e8">d = {d} ({2 * d * d - 1} qubits, {d} rounds)</text>')
svg.append(f'<text x="{lx}" y="{ly + 4 + 3 * 16 + 2}" fill="#7a879c">● our BEAM sampler   ○ Stim, same circuit</text>')
svg.append(f'<text x="{lx}" y="{ly + 4 + 4 * 16}" fill="#7a879c">both decoded by PyMatching, 100,000 shots;</text>')
svg.append(f'<text x="{lx}" y="{ly + 4 + 5 * 16}" fill="#7a879c">points with fewer than 10 failures not drawn</text>')
svg.append("</svg>")
(out / "surface-circuit-noise.svg").write_text("\n".join(svg))
if CHART_ONLY:
    sys.exit(0)

# ── 1b. a stricter check of the sampler: every detector's firing rate, and every pair's joint rate ─────────
# The detector error model (and so the decoding) comes from Stim; what is compared here is our sampling alone.
def zscores(a, b):
    """z of the difference of two rates estimated from the same number of shots"""
    pa, pb = a.mean(0), b.mean(0)
    var = (pa * (1 - pa) + pb * (1 - pb)) / a.shape[0]
    return (pa - pb) / np.sqrt(np.maximum(var, 1e-12))


strict = {"single": [], "pairs": []}
if not CHART_ONLY:
    for run in json.load(open(bench / "study.json")):
        c = stim.Circuit(open(bench / f"{run['name']}.stim").read())
        d_ours, _ = ours(run["name"], c, run["shots"])
        d_stim, _ = c.compile_detector_sampler(seed=run["seed"]).sample(run["shots"], separate_observables=True)
        strict["single"].extend(zscores(d_ours, d_stim).tolist())
        if run["d"] == 5 and run["p"] == 0.006:                          # every pair of 120 detectors: 7,140 pairs
            iu = np.triu_indices(c.num_detectors, 1)
            jo = lambda d: (d[:, :, None] & d[:, None, :])[:, iu[0], iu[1]]
            a = np.concatenate([jo(d_ours[k:k + 5000]) for k in range(0, run["shots"], 5000)])
            b = np.concatenate([jo(d_stim[k:k + 5000]) for k in range(0, run["shots"], 5000)])
            strict["pairs"] = zscores(a, b).tolist()
    zs, zp = np.abs(np.array(strict["single"])), np.abs(np.array(strict["pairs"]))
    validation = {
        "detectors_compared": int(zs.size), "detector_max_abs_z": float(zs.max()),
        "detector_share_over_3": float((zs > 3).mean()), "expected_share_over_3": 0.0027,
        "pairs_compared": int(zp.size), "pair_max_abs_z": float(zp.max()), "pair_share_over_3": float((zp > 3).mean()),
        "note": "z of (our rate - Stim's rate) per detector over all 21 points, and per pair of detectors (joint firing) for d=5, p=0.6 %; under agreement |z| > 3 has probability 0.27 %; the largest of N values grows like sqrt(2 ln N)",
    }
    print("strict (normal approximation; the published exact tests are tools/stim_exact.py):", validation)

# ── 2. the stabilizer benchmark: the same circuit, 100,000 shots ──────────────────────────────────────────
from multiprocessing import Pool  # noqa: E402


def _stim_part(args):
    text, shots, seed = args
    stim.Circuit(text).compile_detector_sampler(seed=seed).sample(shots, separate_observables=True)
    return shots


def med_min(f, reps, warm=True):
    if warm:
        f()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        f()
        ts.append((time.perf_counter() - t0) * 1000)
    return float(np.median(ts)), float(min(ts))


stab = []
with Pool(20) as pool:
    pool.map(_stim_part, [("H 0\nM 0", 1, k) for k in range(20)])          # start the workers before timing
    for run in json.load(open(bench / "stab.json")):
        text = open(bench / f"{run['name']}.stim").read()
        c = stim.Circuit(text)
        one = med_min(lambda: c.compile_detector_sampler(seed=1).sample(run["shots"], separate_observables=True), 5)
        twenty = med_min(lambda: pool.map(_stim_part, [(text, run["shots"] // 20, k) for k in range(20)]), 5)
        stab.append({"circuit": run["name"], "qubits": run["qubits"], "detectors": c.num_detectors, "shots": run["shots"],
                     "beam_median_ms": round(run["median_ms"], 1), "beam_min_ms": round(run["min_ms"], 1),
                     "stim_1core_median_ms": round(one[0], 1), "stim_1core_min_ms": round(one[1], 1),
                     "stim_20proc_median_ms": round(twenty[0], 1), "stim_20proc_min_ms": round(twenty[1], 1)})
        print(stab[-1])

# ── 3. the state-vector benchmark: the same random-layer circuits, 1,024 shots ────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "qsim"))
from qsim.engine import Circuit  # noqa: E402
from qiskit import QuantumCircuit, transpile  # noqa: E402
import qiskit_aer  # noqa: E402
from qiskit_aer import AerSimulator  # noqa: E402

aer = AerSimulator(method="statevector")
sv = []
beam_sv = {r["n"]: r for r in json.load(open(bench / "sv.json"))}
for n in range(4, 21):
    def layers(c, rx, cx):
        for layer in range(n):
            for q in range(n):
                rx(c, q, 0.1 * (q + 1) + 0.01 * layer)
            for q in range(0, n - 1, 2):
                cx(c, q, q + 1)
            for q in range(1, n - 1, 2):
                cx(c, q, q + 1)
        return c
    reps = 3 if n >= 17 else 5
    npy = med_min(lambda: layers(Circuit(n), lambda c, q, a: c.rx(q, a), lambda c, a, b: c.cx(a, b)).sample(1024, np.random.default_rng(1)), reps)
    qc = layers(QuantumCircuit(n), lambda c, q, a: c.rx(a, q), lambda c, a, b: c.cx(a, b))
    qc.measure_all()
    tq = transpile(qc, aer)                                                  # not timed: Qiskit's compilation
    ae = med_min(lambda: aer.run(tq, shots=1024, seed_simulator=1).result(), reps)
    b = beam_sv.get(n)
    sv.append({"qubits": n, "gates": n * n + n * (n - 1), "runs": reps,
               "beam_median_ms": round(b["median_ms"], 1) if b else None, "beam_min_ms": round(b["min_ms"], 1) if b else None,
               "beam_note": "" if b else "out of memory in the 1 GB container",
               "numpy_median_ms": round(npy[0], 2), "numpy_min_ms": round(npy[1], 2),
               "aer_median_ms": round(ae[0], 2), "aer_min_ms": round(ae[1], 2)})
    print(sv[-1])

with open(out / "benchmarks.csv", "w", newline="") as f:
    f.write("# state vector: random layers (Rx on every qubit, CNOTs on even then odd pairs) x n layers, 1,024 shots;"
            " one warm-up, then 5 runs (3 from 17 qubits); median and minimum in ms\n")
    f.write("# timed: QxSim building the state and sampling inside the VM; the NumPy engine the same; Aer's run, its transpile excluded\n")
    w = csv.DictWriter(f, fieldnames=list(sv[0]))
    w.writeheader()
    w.writerows(sv)
    f.write("\n# stabilizer: memory circuits (circuit-level noise p = 0.001, d rounds), 100,000 shots; one warm-up, 5 runs\n")
    f.write("# timed: QxStab reference tableau + sampling in parallel batches on 20 schedulers; Stim compile + sample, on one core,"
            " and in 20 worker processes of 5,000 shots each (pool already started)\n")
    w = csv.DictWriter(f, fieldnames=list(stab[0]))
    w.writeheader()
    w.writerows(stab)
json.dump({"sv": sv, "stab": stab, "versions": VERS + f", qiskit-aer {qiskit_aer.__version__} (statevector)",
           "date": time.strftime("%Y-%m-%d"),
           "where": "one machine, 20 Arm (aarch64) cores; QxSim and QxStab in the quanterall_quantum container (1 GB memory limit, 20 schedulers), Python on the host"},
          open(out / "benchmarks.json", "w"), indent=1)
print("done:", VERS)
