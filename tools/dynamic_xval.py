"""QxSim's measured circuits against Qiskit Aer's dynamic circuits (quanterall.com, Results).

The error-correcting codes on QxSim measure ancillas mid-circuit, reset them, and correct from the syndrome;
teleportation corrects Bob's qubit from Alice's bits. Here each is rebuilt in Qiskit as a dynamic circuit —
real mid-circuit measurements, resets, and corrections conditioned on the measured bits (for teleportation,
measured bits, not the deferred-measurement form QxSim uses) — and run on Aer. For every code and every
injected error, the full data-readout distributions (with the error left in, and corrected) are compared:
total variation distance, and a chi-square test of homogeneity.

    ~/quanterall/.venv-tools/bin/python tools/dynamic_xval.py http://<quanterall_quantum>:8080 site/data BENCH_DIR
"""
import json
import math
import sys
import urllib.request
from pathlib import Path

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit.classical import expr
from qiskit_aer import AerSimulator
from scipy.stats import chi2_contingency

BASE, OUT = sys.argv[1], Path(sys.argv[2])
SHOTS = 4000
aer = AerSimulator()


def beam(path):
    return json.load(urllib.request.urlopen(f"{BASE}/quantum/{path}"))


def counts_of(qc, reg_index, n_regs, seed=7):
    """the counts of one classical register (Qiskit prints registers last-added first, separated by spaces)"""
    raw = aer.run(qc, shots=SHOTS, seed_simulator=seed).result().get_counts()
    out = {}
    for key, c in raw.items():
        part = key.split(" ")[n_regs - 1 - reg_index]
        out[part] = out.get(part, 0) + c
    return out


def compare(a, b):
    keys = sorted(set(a) | set(b))
    va, vb = np.array([a.get(k, 0) for k in keys]), np.array([b.get(k, 0) for k in keys])
    tvd = 0.5 * np.abs(va / va.sum() - vb / vb.sum()).sum()
    if len(keys) == 1:
        return float(tvd), 1.0
    _, p, _, _ = chi2_contingency(np.vstack([va, vb]))
    return float(tvd), float(p)


def check(qc, anc, reg, kind, qs, reset=True):
    """one stabilizer check on the ancilla, measured into reg's next bit, the ancilla reset (as QxSim)"""
    if kind == "z":
        for q in qs:
            qc.cx(q, anc)
    else:
        qc.h(anc)
        for q in qs:
            qc.cx(anc, q)
        qc.h(anc)
    qc.measure(anc, reg)
    if reset:
        qc.reset(anc)


def pauli(qc, g, q):
    {"x": qc.x, "y": qc.y, "z": qc.z}[g](q)


# ── the codes, as in QxSim's qec_cycle ────────────────────────────────────────────────────────────────────
REP = {(1, 0): 0, (1, 1): 1, (0, 1): 2}
SURF_X, SURF_Z = [[0, 1, 3, 4], [4, 5, 7, 8], [1, 2], [6, 7]], [[1, 2, 4, 5], [3, 4, 6, 7], [0, 3], [5, 8]]
STEANE = [[0, 2, 4, 6], [1, 2, 5, 6], [3, 4, 5, 6]]


def build(code, err, correct, bug=False):
    """the dynamic circuit; returns (circuit, index of the data register, number of registers)"""
    k = {"phaseflip": 3, "shor9": 9, "steane7": 7, "surface3": 9}[code]
    q = QuantumRegister(k + 1, "q")
    anc = k
    regs = []

    def creg(name, size):
        r = ClassicalRegister(size, name)
        regs.append(r)
        return r

    if code == "phaseflip":
        s = creg("s", 2)
        qc = QuantumCircuit(q, *regs)
        qc.cx(0, 1); qc.cx(0, 2); qc.h(0); qc.h(1); qc.h(2)
        if err is not None:
            qc.z(err)
        check(qc, anc, s[0], "x", [0, 1]); check(qc, anc, s[1], "x", [1, 2])
        if correct:
            for (a, b), qq in REP.items():
                with qc.if_test((s, a + 2 * b)):
                    qc.z(qq)
        qc.h(0); qc.h(1); qc.h(2)
    elif code == "shor9":
        zb = [creg(f"zb{i}", 2) for i in range(3)]
        xs = creg("xs", 2)
        qc = QuantumCircuit(q, *regs)
        qc.cx(0, 3); qc.cx(0, 6); qc.h(0); qc.h(3); qc.h(6)
        for a, b in ((0, 1), (0, 2), (3, 4), (3, 5), (6, 7), (6, 8)):
            qc.cx(a, b)
        if err is not None:
            qc.y(err)
        for blk in range(3):
            check(qc, anc, zb[blk][0], "z", [3 * blk, 3 * blk + 1]); check(qc, anc, zb[blk][1], "z", [3 * blk + 1, 3 * blk + 2])
        check(qc, anc, xs[0], "x", list(range(0, 6))); check(qc, anc, xs[1], "x", list(range(3, 9)))
        if correct:
            for blk in range(3):
                for (a, b), qq in REP.items():
                    with qc.if_test((zb[blk], a + 2 * b)):
                        qc.x(3 * blk + qq)
            for (a, b), blk in REP.items():
                with qc.if_test((xs, a + 2 * b)):
                    qc.z(3 * blk)
        for i in range(9):
            qc.h(i)
    elif code == "steane7":
        zr, xr = creg("zr", 3), creg("xr", 3)
        qc = QuantumCircuit(q, *regs)
        qc.h(0); qc.h(1); qc.h(3)
        for a, b in ((3, 4), (3, 5), (3, 6), (1, 2), (1, 5), (1, 6), (0, 2), (0, 4), (0, 6)):
            qc.cx(a, b)
        if err is not None:
            qc.y(err)
        for i, g in enumerate(STEANE):
            check(qc, anc, zr[i], "z", g)
        for i, g in enumerate(STEANE):
            check(qc, anc, xr[i], "x", g)
        if correct:
            for v in range(1, 8):
                with qc.if_test((zr, v)):
                    qc.x(v % 7 if bug else v - 1)                    # the bug: X one position off
                with qc.if_test((xr, v)):
                    qc.z(v - 1)
    else:  # surface3: a reference round, the error, a second round; the syndrome is their XOR
        rx, rz, sx, sz = creg("rx", 4), creg("rz", 4), creg("sx", 4), creg("sz", 4)
        qc = QuantumCircuit(q, *regs)
        for i, g in enumerate(SURF_X):
            check(qc, anc, rx[i], "x", g)
        for i, g in enumerate(SURF_Z):
            check(qc, anc, rz[i], "z", g)
        if err is not None:
            qc.y(err)
        for i, g in enumerate(SURF_X):
            check(qc, anc, sx[i], "x", g)
        for i, g in enumerate(SURF_Z):
            check(qc, anc, sz[i], "z", g)
        if correct:
            sig = lambda cs, qq: sum(1 << i for i, c in enumerate(cs) if qq in c)
            xfix = {}
            zfix = {}
            for qq in range(8, -1, -1):
                xfix[sig(SURF_Z, qq)] = qq
                zfix[sig(SURF_X, qq)] = qq
            for v, qq in xfix.items():
                if v:
                    with qc.if_test(expr.equal(expr.bit_xor(sz, rz), v)):
                        qc.x(qq)
            for v, qq in zfix.items():
                if v:
                    with qc.if_test(expr.equal(expr.bit_xor(sx, rx), v)):
                        qc.z(qq)
    c = ClassicalRegister(k, "c")
    qc.add_register(c)
    for i in range(k):
        qc.measure(i, c[i])
    return qc, len(regs), len(regs) + 1


rows = []
for code, k in (("phaseflip", 3), ("shor9", 9), ("steane7", 7), ("surface3", 9)):
    for e in [None] + list(range(k)):
        # every comparison its own seeds, on both sides: otherwise identical circuits repeat one sample
        sd = 1000 * len(rows) + 17
        b = beam(f"qec/{code}?shots={SHOTS}&error_qubit={'none' if e is None else e}&seed={sd}")
        for stage, correct, ours in (("with the error", False, b["data_before"]), ("corrected", True, b["data_after"])):
            qc, idx, nregs = build(code, e, correct)
            ref = counts_of(qc, idx, nregs, seed=sd + (1 if correct else 2))
            tvd, p = compare(ours, ref)
            rows.append({"code": code, "error_qubit": e, "stage": stage, "outcomes_ours": len(ours), "outcomes_aer": len(ref),
                         "tvd": round(tvd, 4), "chi2_p": round(p, 4)})
        print(code, e, rows[-2]["tvd"], rows[-1]["tvd"], rows[-1]["chi2_p"])

# teleportation: QxSim's deferred-measurement circuit against Aer with Alice's bits measured and fed forward
def teleport_aer(th, ph, seed, wrong=False):
    q, m, bob = QuantumRegister(3, "q"), ClassicalRegister(2, "m"), ClassicalRegister(1, "bob")
    qc = QuantumCircuit(q, m, bob)
    qc.ry(th, 0); qc.rz(ph, 0); qc.h(1); qc.cx(1, 2); qc.cx(0, 1); qc.h(0)
    qc.measure(0, m[0]); qc.measure(1, m[1])
    with qc.if_test((m[0] if wrong else m[1], 1)):                  # the wrong version swaps which bit drives X and Z
        qc.x(2)
    with qc.if_test((m[1] if wrong else m[0], 1)):
        qc.z(2)
    qc.rz(-ph, 2); qc.ry(-th, 2)
    qc.measure(2, bob[0])
    raw = aer.run(qc, shots=SHOTS, seed_simulator=seed).result().get_counts()
    ref = {}
    for key, cc in raw.items():
        bb, mm = key.split(" ")
        ref[bb + mm] = ref.get(bb + mm, 0) + cc
    return ref


for seed in range(1, 6):
    b = beam(f"teleport?shots={SHOTS}&seed={seed}")
    th, ph = math.radians(b["state"]["theta_deg"]), math.radians(b["state"]["phi_deg"])
    ref = teleport_aer(th, ph, seed)
    tvd, p = compare(b["counts"], ref)
    rows.append({"code": "teleport", "error_qubit": None, "stage": f"random ψ, seed {seed}", "outcomes_ours": len(b["counts"]),
                 "outcomes_aer": len(ref), "tvd": round(tvd, 4), "chi2_p": round(p, 4)})
    print("teleport", seed, tvd, p)

for row in rows:
    row["deterministic"] = row["outcomes_ours"] == 1 and row["outcomes_aer"] == 1
det = [r for r in rows if r["deterministic"]]
rnd = [r for r in rows if not r["deterministic"]]
ps = [r["chi2_p"] for r in rnd]

# power: each check must catch a known bug
power = []
b = beam(f"teleport?shots={SHOTS}&seed=1")
th, ph = math.radians(b["state"]["theta_deg"]), math.radians(b["state"]["phi_deg"])
t, p = compare(b["counts"], teleport_aer(th, ph, 1, wrong=True))
power.append({"bug": "teleportation, Aer with X and Z driven by the wrong bits", "tvd": round(t, 4), "chi2_p": p})
b = beam(f"qec/steane7?shots={SHOTS}&error_qubit=3&seed=11")
qc, idx, nregs = build("steane7", 3, True, bug=True)
t, p = compare(b["data_after"], counts_of(qc, idx, nregs))
power.append({"bug": "Steane, Aer's decoder off by one position (error on qubit 3)", "tvd": round(t, 4), "chi2_p": p})
from qiskit.quantum_info import Statevector  # noqa: E402
case = json.load(open(Path(sys.argv[3]) / "xval.json"))[-1] if len(sys.argv) > 3 else None
if case:
    qn = QuantumCircuit(case["n"])
    first = True
    for g in case["gates"]:
        k = g[0]
        if k in ("rx", "ry", "rz"):
            getattr(qn, k)(g[2] + (0.01 if first else 0), g[1]); first = False
        elif k == "cp":
            qn.cp(g[3], g[1], g[2])
        elif k in ("cx", "cz", "swap"):
            getattr(qn, k)(g[1], g[2])
        else:
            getattr(qn, k)(g[1])
    ours = np.array([complex(x, y) for x, y in case["amps"]])
    fid = abs(np.vdot(Statevector(qn).data, ours)) ** 2
    power.append({"bug": f"state vectors: one rotation angle changed by 0.01 rad ({case['n']} qubits)", "fidelity_gap": float(1 - fid)})

summary = {"comparisons": len(rows), "shots_each": SHOTS,
           "deterministic": {"count": len(det), "all_equal": all(r["tvd"] == 0 for r in det)},
           "random": {"count": len(rnd), "min_chi2_p": min(ps), "max_tvd": max(r["tvd"] for r in rnd),
                      "p_histogram_tenths": np.histogram(ps, bins=10, range=(0, 1))[0].tolist()},
           "power_checks": power,
           "note": "QxSim's data readouts (with the error left in, and corrected) against Qiskit Aer dynamic circuits "
                   "with mid-circuit measurement, reset and classically conditioned corrections; teleportation with "
                   "Alice's bits measured and fed forward. Deterministic comparisons (one outcome on both sides) are "
                   "checked for equality; for the others the chi-square p-values should be roughly uniform. Many "
                   "comparisons share a circuit after correction, so they are not independent.",
           "rows": rows}
json.dump(summary, open(OUT / "dynamic-vs-aer.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))
