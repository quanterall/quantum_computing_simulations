"""Parity: the NumPy back end against the BEAM one, demo by demo, as distributions (total variation
distance between the two sampled histograms; identical circuits give sampling noise only, ~0.01-0.03
at 10,000 shots).   python3 parity.py BEAM_BASE QSIM_BASE"""
import json, sys, urllib.request

beam, qsim = sys.argv[1], sys.argv[2]
SHOTS = 10000
CASES = [(f"/quantum/{d}", None) for d in ("bell", "teleport", "grover", "qml", "shor15")] + \
        [(f"/quantum/qec/{d}", e) for d in ("bitflip", "phaseflip", "shor9", "steane7", "surface3") for e in ("none", "1")]


def get(base, path, eq):
    q = f"?shots={SHOTS}" + (f"&error_qubit={eq}" if eq else "")
    with urllib.request.urlopen(base + path + q, timeout=60) as r:
        return json.loads(r.read())


def tvd(a, b):
    ka = set(a) | set(b); na, nb = sum(a.values()), sum(b.values())
    return 0.5 * sum(abs(a.get(k, 0) / na - b.get(k, 0) / nb) for k in ka)


worst = 0.0
for path, eq in CASES:
    A, B = get(beam, path, eq), get(qsim, path, eq)
    fields = ["counts"] if "counts" in A else ["counts_before", "syndrome_counts", "corrected_counts"]
    d = max(tvd(A[f], B[f]) for f in fields)
    worst = max(worst, d)
    extra = f" fidelity {A['logical_fidelity']} / {B['logical_fidelity']}" if "logical_fidelity" in A else ""
    print(f"{path:24s} err={eq or '-':4s} TVD {d:.4f}  time {A['execution_time_ms']} / {B['execution_time_ms']} ms{extra}")
print(f"worst TVD {worst:.4f}", "OK" if worst < 0.05 else "MISMATCH")
