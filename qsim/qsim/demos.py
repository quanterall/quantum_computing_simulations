"""The ten Quanterall demos on the NumPy engine.

Gate for gate the circuits of the BEAM demo suite (quanterall_quantum/lib/qx_demo.ex, version 3.0), with
its response fields, so the demo page reads them unchanged: five quantum-AI circuits and five
quantum-error-correction (QEC) demos. The QEC post-processing is the original's too: a majority vote over
each measured string, and the share of all-zero strings after it as the "logical fidelity". These are
teaching circuits, simplified as the original says (no full encoders, no noise model).
"""
from __future__ import annotations

import copy
import time
from datetime import datetime, timezone
from math import gcd

import numpy as np

from .engine import Circuit, X, Y, Z, run

VERSION = "4.0"
ENGINE = "numpy state-vector (qsim)"


# ── responses ─────────────────────────────────────────────────────────────────────────────────────────────

def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def vega(counts: dict[str, int]) -> dict:
    return {
        "mark": "bar",
        "data": {"values": [{"state": s, "count": c} for s, c in counts.items()]},
        "encoding": {
            "x": {"field": "state", "type": "nominal", "title": "State"},
            "y": {"field": "count", "type": "quantitative", "title": "Counts"},
        },
    }


def resp(name, qubits, shots, ms, counts, description, qasm) -> dict:
    return {
        "status": "ok", "name": name, "qubits": qubits, "shots": shots,
        "execution_time_ms": round(ms, 2),
        "counts": counts,
        "visualisation_hint": "histogram",
        "plots": {"counts": vega(counts)},
        "circuit_qasm": qasm,
        "metadata": {"description": description, "version": VERSION, "engine": ENGINE, "date": _now()},
    }


def majority_decode(counts: dict[str, int]) -> dict[str, int]:
    out: dict[str, int] = {}
    for bits, c in counts.items():
        w = len(bits)
        key = "1" * w if bits.count("1") > w // 2 else "0" * w
        out[key] = out.get(key, 0) + c
    return out


def fidelity(corrected: dict[str, int]) -> float:
    total = sum(corrected.values())
    if not total:
        return 0.0
    w = len(next(iter(corrected)))
    return round(corrected.get("0" * w, 0) / total, 4)


def qec_resp(name, qubits, shots, ms, before, syndrome, error_type, error_qubit, qasm) -> dict:
    corrected = majority_decode(before)
    return {
        "status": "ok", "name": name, "qubits": qubits, "shots": shots,
        "execution_time_ms": round(ms, 2),
        "error_type": error_type, "error_qubit": error_qubit,
        "counts_before": before, "syndrome_counts": syndrome, "corrected_counts": corrected,
        "logical_fidelity": fidelity(corrected),
        "visualisation_hint": "qec",
        "plots": {"counts_before": vega(before), "corrected_counts": vega(corrected), "syndrome": vega(syndrome)},
        "circuit_qasm": qasm,
        "metadata": {"description": f"{name} QEC demo", "code_distance": 3, "error_correction": "majority_vote",
                     "version": VERSION, "engine": ENGINE, "date": _now()},
    }


def parse_error_qubit(eq, n: int):
    """'none', or a qubit index 0..n-1 (anything else: no error), as the original"""
    if eq is None or eq == "none":
        return None
    try:
        q = int(str(eq).strip())
    except ValueError:
        return None
    return q if 0 <= q < n else None


def inject(c, n: int, eq, gate: str):
    q = parse_error_qubit(eq, n)
    if q is not None:
        c.x(q) if gate == "x" else c.z(q)
    return c


# ── quantum-AI demos ─────────────────────────────────────────────────────────────────────────────────────

def bell(shots):
    n = 2
    counts, ms = run(n, shots, lambda c: c.h(0).cx(0, 1))
    return resp("bell", n, shots, ms, counts,
                "Bell state |Φ⁺⟩ = (|00⟩+|11⟩)/√2 — maximally entangled pair",
                "OPENQASM 3.0;\nqubit[2] q; bit[2] c;\nh q[0];\ncx q[0],q[1];\nc[0]=measure q[0]; c[1]=measure q[1];\n")


def teleport(shots):
    """A random state ψ = Rz(φ)·Ry(θ)|0⟩ on Alice's q0, a Bell pair q1–q2, Alice's Bell measurement, Bob's
    corrections X and Z as controlled gates (deferred measurement); then Bob undoes the preparation on q2, so q2
    reads 0 in every shot if ψ arrived intact. The same circuit as quanterall.com's BEAM simulator."""
    n = 3
    rng = np.random.default_rng()
    theta, phi = float(np.arccos(1 - 2 * rng.random())), float(2 * np.pi * rng.random())
    counts, ms = run(n, shots, lambda c: c.ry(0, theta).rz(0, phi).h(1).cx(1, 2).cx(0, 1).h(0)
                     .cx(1, 2).cz(0, 2).rz(2, -phi).ry(2, -theta))
    td, pd = round(np.degrees(theta), 2), round(np.degrees(phi), 2)
    out = resp("teleport", n, shots, ms, counts,
               f"Quantum teleportation of a random state (θ = {td}°, φ = {pd}°): Bob undoes the preparation, so q2 reads 0 if ψ arrived",
               f"OPENQASM 3.0;\nqubit[3] q; bit[3] c;\nry({theta:.6f}) q[0]; rz({phi:.6f}) q[0];   // ψ, random\nh q[1]; cx q[1],q[2];        // the Bell pair\ncx q[0],q[1]; h q[0];        // Alice's Bell measurement\ncx q[1],q[2]; cz q[0],q[2];  // Bob's corrections X and Z (deferred measurement)\nrz({-phi:.6f}) q[2]; ry({-theta:.6f}) q[2];   // Bob undoes the preparation\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2];   // q[2] = 0 every time\n")
    out["state"] = {"theta_deg": td, "phi_deg": pd}
    out["bob_zero_fraction"] = round(sum(m for k, m in counts.items() if k[0] == "0") / shots, 4)
    return out


def grover(shots):
    n = 2
    counts, ms = run(n, shots, lambda c: c.h(0).h(1).cz(0, 1).h(0).h(1).z(0).z(1).cz(0, 1).h(0).h(1))
    return resp("grover", n, shots, ms, counts,
                "Grover's search — oracle marks |11⟩, amplitude amplification boosts it",
                "OPENQASM 3.0;\nqubit[2] q; bit[2] c;\nh q[0]; h q[1];\ncz q[0],q[1];\nh q[0]; h q[1]; z q[0]; z q[1];\ncz q[0],q[1]; h q[0]; h q[1];\nc[0]=measure q[0]; c[1]=measure q[1];\n")


def qml(shots):
    n = 4
    counts, ms = run(n, shots, lambda c: c.rx(0, 0.5).rx(1, 0.3).cx(0, 2).cx(1, 3).rx(2, 0.7).rx(3, 0.7).cx(2, 3))
    return resp("qml", n, shots, ms, counts,
                "Variational quantum classifier — Rx rotations + CNOT entanglement layer",
                "OPENQASM 3.0;\nqubit[4] q; bit[4] c;\nrx(0.5) q[0]; rx(0.3) q[1];\ncx q[0],q[2]; cx q[1],q[3];\nrx(0.7) q[2]; rx(0.7) q[3]; cx q[2],q[3];\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2]; c[3]=measure q[3];\n")


def shor15(shots):
    """Shor's algorithm for N = 15, a = 7: counting q0..q2, work q3..q6 = |1⟩; controlled ×7 (U) from q0, ×4 (U²)
    from q1, U⁴ = identity; the inverse QFT; each measured k/8 = s/r gives the period 4 and gcd(7² ± 1, 15) = 3, 5.
    The same circuit as quanterall.com's BEAM simulator."""
    n = 7

    def build(c):
        c.h(0).h(1).h(2).x(3)
        c.cswap(0, 3, 4).cswap(0, 4, 5).cswap(0, 5, 6).cx(0, 3).cx(0, 4).cx(0, 5).cx(0, 6)
        c.cswap(1, 4, 6).cswap(1, 3, 5)
        c.swap(0, 2).h(0).cp(0, 1, -np.pi / 2).h(1).cp(0, 2, -np.pi / 4).cp(1, 2, -np.pi / 2).h(2)

    full, ms = run(n, shots, build)
    counting: dict[str, int] = {}
    for bits, m in full.items():                         # the counting register only: the last three characters
        counting[bits[-3:]] = counting.get(bits[-3:], 0) + m
    out = resp("shor15", 3, shots, ms, counting,
               "Shor's algorithm, N = 15, a = 7 — controlled modular exponentiation, inverse QFT, continued fractions: 15 = 3 × 5",
               "OPENQASM 3.0;\n// Shor N=15, a=7: counting q[0..2], work q[3..6]\nqubit[7] q; bit[3] c;\nh q[0]; h q[1]; h q[2]; x q[3];\n// controlled ×7 mod 15 (U) from q[0]\ncswap q[0],q[3],q[4]; cswap q[0],q[4],q[5]; cswap q[0],q[5],q[6];\ncx q[0],q[3]; cx q[0],q[4]; cx q[0],q[5]; cx q[0],q[6];\n// controlled ×4 mod 15 (U²) from q[1]; U⁴ = identity\ncswap q[1],q[4],q[6]; cswap q[1],q[3],q[5];\n// inverse QFT on q[0..2]\nswap q[0],q[2];\nh q[0];\ncp(-pi/2) q[0],q[1]; h q[1];\ncp(-pi/4) q[0],q[2]; cp(-pi/2) q[1],q[2]; h q[2];\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2];\n")
    out["classical"] = _shor_classical(counting)
    out["qubits_total"] = n
    return out


def _shor_classical(counting):
    """each measured k: the phase k/8, the period it gives, and the factors (or why it gives none)"""
    rows = []
    for bits, m in sorted(counting.items()):
        k = int(bits, 2)
        if k == 0:
            r, note = None, "phase 0: no information, run again"
        else:
            r0 = 8 // gcd(k, 8)                                       # k/8 = s/r in lowest terms
            r = next(x for x in (r0, 2 * r0, 4 * r0) if pow(7, x, 15) == 1)
            note = f"period {r}" if r == r0 else f"{k}/8 = 1/{r0}: 7^{r0} ≠ 1 mod 15, so the period is a multiple: {r}"
        factors = ([gcd(7 ** (r // 2) - 1, 15), gcd(7 ** (r // 2) + 1, 15)]
                   if r is not None and r % 2 == 0 and pow(7, r // 2, 15) != 14 else [])
        rows.append({"measured": bits, "k": k, "count": m, "phase": f"{k}/8", "period": r, "factors": factors, "note": note})
    return rows


# ── QEC demos ────────────────────────────────────────────────────────────────────────────────────────────

def bitflip(shots, eq="none"):
    n, ns = 3, 5
    before, ms1 = run(n, shots, lambda c: inject(c, n, eq, "x"))
    syn, ms2 = run(ns, shots, lambda c: inject(c, ns, eq, "x").cx(0, 3).cx(1, 3).cx(1, 4).cx(2, 4))
    return qec_resp("bitflip", n, shots, ms1 + ms2, before, syn, "bit_flip", eq,
                    "OPENQASM 3.0;\n// 3-qubit bit-flip repetition code\nqubit[5] q; bit[2] s;\ncx q[0],q[3]; cx q[1],q[3];\ncx q[1],q[4]; cx q[2],q[4];\ns[0]=measure q[3]; s[1]=measure q[4];\n")


# ── full correction cycles: phase-flip, Shor-9, Steane-7 ─────────────────────────────────────────────────
# encode logical |0⟩ → inject the error → measure every stabiliser check with one ancilla, reset after each
# check → look the syndrome up in the decoder → apply the correction → read the logical qubit out in the code's
# basis. The same circuit without the correction gives "before". With noiseless gates a single Pauli error's
# syndrome is certain, so it is measured once on the state vector and holds for every shot. The same cycle
# as quanterall.com's BEAM simulator (QxDemo.qec_cycle).

def _parity(bits: str) -> str:
    return str(bits.count("1") % 2)


def _majority(bits: str) -> str:
    return "1" if 2 * bits.count("1") > len(bits) else "0"


def _check(c: Circuit, anc: int, kind: str, qs) -> int:
    """one stabiliser check on the ancilla: Z-type copies the parity in; X-type runs it in the Hadamard basis"""
    if kind == "z":
        for q in qs:
            c.cx(q, anc)
    else:
        c.h(anc)
        for q in qs:
            c.cx(anc, q)
        c.h(anc)
    return c.measure_reset(anc)


def _cycle(name, k, error, encode, checks, decode, readout, logical, keep_data, qasm, shots, eq, error_type,
           reference=False):
    t0 = time.perf_counter()
    n, anc = k + 1, k
    err = parse_error_qubit(eq, k)
    c = encode(Circuit(n))
    ref = [_check(c, anc, kind, qs) for kind, qs in checks] if reference else [0] * len(checks)   # a code made by measurement
    if err is not None:
        c.gate({"x": X, "y": Y, "z": Z}[error], err)
    syn = [a ^ b for a, b in zip((_check(c, anc, kind, qs) for kind, qs in checks), ref)]
    fixes = decode(syn)
    before_c = copy.deepcopy(c)
    for g, q in fixes:
        c.gate({"x": X, "z": Z}[g], q)
    before = readout(before_c).sample_low(k, shots)
    fixed = readout(c).sample_low(k, shots)

    def logical_counts(counts):
        out: dict[str, int] = {}
        for b, m in counts.items():
            out[logical(b)] = out.get(logical(b), 0) + m
        return out

    lb, la = logical_counts(before), logical_counts(fixed)
    cb, ca = (before, fixed) if keep_data else (lb, la)
    s = "".join(map(str, syn))
    return {
        "status": "ok", "name": name, "qubits": n, "shots": shots,
        "execution_time_ms": round((time.perf_counter() - t0) * 1000, 2),
        "error_type": error_type, "error_qubit": eq,
        "counts_before": cb, "syndrome_counts": {s: shots}, "corrected_counts": ca,
        "logical_before": lb, "logical_after": la,
        "logical_fidelity": round(la.get("0", 0) / shots, 4),
        "syndrome": s, "correction": [f"{g} q[{q}]" for g, q in fixes],
        "visualisation_hint": "qec",
        "plots": {"counts_before": vega(cb), "corrected_counts": vega(ca), "syndrome": vega({s: shots})},
        "circuit_qasm": qasm,
        "metadata": {"description": f"{name}: encode |0_L⟩, inject {error}, measure the stabiliser checks with a "
                                    "reused ancilla, decode the syndrome, correct, read the logical qubit",
                     "code_distance": 3, "error_correction": "syndrome decoding (lookup table)",
                     "readout": "data qubits in the X basis" if keep_data else "logical qubit: 0 or 1",
                     "version": VERSION, "engine": ENGINE, "date": _now()},
    }


_REP = {(1, 0): 0, (1, 1): 1, (0, 1): 2}                   # a repetition code's two checks → the qubit at fault


def _hs(c, qs):
    for q in qs:
        c.h(q)
    return c


def phaseflip(shots, eq="none"):
    """|0_L⟩ = |+++⟩: a repetition code in the X basis. Checks X0X1 and X1X2 find a Z error; Z corrects it."""
    return _cycle(
        "phaseflip", 3, "z",
        encode=lambda c: _hs(c.cx(0, 1).cx(0, 2), [0, 1, 2]),
        checks=[("x", [0, 1]), ("x", [1, 2])],
        decode=lambda s: [("z", _REP[tuple(s)])] if tuple(s) in _REP else [],
        readout=lambda c: _hs(c, [0, 1, 2]), logical=_majority, keep_data=True,
        qasm="OPENQASM 3.0;\n// 3-qubit phase-flip code: encode |0_L⟩ = |+++⟩, checks X0X1, X1X2 on a reused ancilla q[3]\nqubit[4] q; bit[2] s; bit[3] c;\ncx q[0],q[1]; cx q[0],q[2]; h q[0]; h q[1]; h q[2];\n// error: z q[e]\nh q[3]; cx q[3],q[0]; cx q[3],q[1]; h q[3]; s[0]=measure q[3]; reset q[3];\nh q[3]; cx q[3],q[1]; cx q[3],q[2]; h q[3]; s[1]=measure q[3]; reset q[3];\n// correction: z on the qubit the syndrome names (10→0, 11→1, 01→2)\nh q[0]; h q[1]; h q[2]; c = measure q[0:2];\n",
        shots=shots, eq=eq, error_type="phase_flip")


def _shor9_decode(s):
    fixes = [("x", 3 * blk + _REP[(s[2 * blk], s[2 * blk + 1])]) for blk in range(3) if (s[2 * blk], s[2 * blk + 1]) in _REP]
    if (s[6], s[7]) in _REP:
        fixes.append(("z", 3 * _REP[(s[6], s[7])]))
    return fixes


def shor9(shots, eq="none"):
    """Shor's code: a phase-flip code over three bit-flip codes; the error is Y, both at once"""
    return _cycle(
        "shor9", 9, "y",
        encode=lambda c: _hs(c.cx(0, 3).cx(0, 6), [0, 3, 6]).cx(0, 1).cx(0, 2).cx(3, 4).cx(3, 5).cx(6, 7).cx(6, 8),
        checks=[("z", [0, 1]), ("z", [1, 2]), ("z", [3, 4]), ("z", [4, 5]), ("z", [6, 7]), ("z", [7, 8]),
                ("x", [0, 1, 2, 3, 4, 5]), ("x", [3, 4, 5, 6, 7, 8])],
        decode=_shor9_decode, readout=lambda c: _hs(c, range(9)), logical=_parity, keep_data=False,
        qasm="OPENQASM 3.0;\n// Shor 9-qubit code: encode |0_L⟩, six Z checks and two X checks on a reused ancilla q[9]\nqubit[10] q; bit[8] s; bit[9] c;\ncx q[0],q[3]; cx q[0],q[6]; h q[0]; h q[3]; h q[6];\ncx q[0],q[1]; cx q[0],q[2]; cx q[3],q[4]; cx q[3],q[5]; cx q[6],q[7]; cx q[6],q[8];\n// error: y q[e]\n// Z checks (bit flips), e.g. Z0Z1:\ncx q[0],q[9]; cx q[1],q[9]; s[0]=measure q[9]; reset q[9];\n// ... Z1Z2, Z3Z4, Z4Z5, Z6Z7, Z7Z8 alike\n// X checks (phase flips), X0..X5 and X3..X8:\nh q[9]; cx q[9],q[0]; cx q[9],q[1]; cx q[9],q[2]; cx q[9],q[3]; cx q[9],q[4]; cx q[9],q[5]; h q[9]; s[6]=measure q[9]; reset q[9];\nh q[9]; cx q[9],q[3]; cx q[9],q[4]; cx q[9],q[5]; cx q[9],q[6]; cx q[9],q[7]; cx q[9],q[8]; h q[9]; s[7]=measure q[9]; reset q[9];\n// correction: x on the flipped qubit of a block, z on the first qubit of the flipped block\nh q[0:8]; c = measure q[0:8];   // logical bit = parity of c\n",
        shots=shots, eq=eq, error_type="bit_phase")


_STEANE = [[0, 2, 4, 6], [1, 2, 5, 6], [3, 4, 5, 6]]       # positions 1,3,5,7 · 2,3,6,7 · 4,5,6,7


def _steane_decode(s):
    zpos, xpos = s[0] + 2 * s[1] + 4 * s[2], s[3] + 2 * s[4] + 4 * s[5]
    return ([("x", zpos - 1)] if zpos else []) + ([("z", xpos - 1)] if xpos else [])


def steane7(shots, eq="none"):
    """Steane's [[7,1,3]] code from the Hamming code: each syndrome, in binary, is the error's position"""
    return _cycle(
        "steane7", 7, "y",
        encode=lambda c: _hs(c, [0, 1, 3]).cx(3, 4).cx(3, 5).cx(3, 6).cx(1, 2).cx(1, 5).cx(1, 6).cx(0, 2).cx(0, 4).cx(0, 6),
        checks=[("z", g) for g in _STEANE] + [("x", g) for g in _STEANE],
        decode=_steane_decode, readout=lambda c: c, logical=_parity, keep_data=False,
        qasm="OPENQASM 3.0;\n// Steane [[7,1,3]] code: encode |0_L⟩, three Z and three X checks on a reused ancilla q[7]\nqubit[8] q; bit[6] s; bit[7] c;\nh q[0]; h q[1]; h q[3];\ncx q[3],q[4]; cx q[3],q[5]; cx q[3],q[6]; cx q[1],q[2]; cx q[1],q[5]; cx q[1],q[6]; cx q[0],q[2]; cx q[0],q[4]; cx q[0],q[6];\n// error: y q[e]\n// Z checks on {0,2,4,6}, {1,2,5,6}, {3,4,5,6}, e.g.:\ncx q[0],q[7]; cx q[2],q[7]; cx q[4],q[7]; cx q[6],q[7]; s[0]=measure q[7]; reset q[7];\n// X checks on the same supports, e.g.:\nh q[7]; cx q[7],q[0]; cx q[7],q[2]; cx q[7],q[4]; cx q[7],q[6]; h q[7]; s[3]=measure q[7]; reset q[7];\n// correction: x at position s[0]+2s[1]+4s[2], z at position s[3]+2s[4]+4s[5]\nc = measure q[0:6];   // logical bit = parity of c\n",
        shots=shots, eq=eq, error_type="bit_phase")


_SURF_X = [[0, 1, 3, 4], [4, 5, 7, 8], [1, 2], [6, 7]]
_SURF_Z = [[1, 2, 4, 5], [3, 4, 6, 7], [0, 3], [5, 8]]


def _surface3_decode(s):
    sig = lambda cs, q: tuple(int(q in c) for c in cs)
    xfix = {sig(_SURF_Z, q): q for q in range(8, -1, -1)}        # Z checks see X errors (lowest qubit kept)
    zfix = {sig(_SURF_X, q): q for q in range(8, -1, -1)}        # X checks see Z errors
    xb, zb = tuple(s[:4]), tuple(s[4:])
    return ([("x", xfix[zb])] if any(zb) else []) + ([("z", zfix[xb])] if any(xb) else [])


def surface3(shots, eq="none"):
    """The rotated surface code of distance 3: data q0..q8 on a 3×3 grid, eight stabilizers on one reused ancilla
    q9; logical |0⟩ made by measurement (a reference round), a Y error, a lookup decoder, the logical Z0Z1Z2.
    The same cycle as quanterall.com's BEAM simulator."""
    return _cycle(
        "surface3", 9, "y",
        encode=lambda c: c,
        checks=[("x", g) for g in _SURF_X] + [("z", g) for g in _SURF_Z],
        decode=_surface3_decode, readout=lambda c: c, logical=lambda b: _parity(b[-3:]), keep_data=False,
        qasm="OPENQASM 3.0;\n// Rotated surface code, d = 3: data q[0..8] (3×3, row by row), one reused ancilla q[9]\nqubit[10] q; bit[16] s; bit[9] c;\n// X checks {0,1,3,4} {4,5,7,8} {1,2} {6,7}; Z checks {1,2,4,5} {3,4,6,7} {0,3} {5,8}\n// round 1 (reference) then the error y q[e], round 2; syndrome = round 2 xor round 1; lookup correction\nc = measure q[0:8];   // logical bit = c[0] xor c[1] xor c[2]\n",
        shots=shots, eq=eq, error_type="bit_phase", reference=True)


AI = {"bell": bell, "teleport": teleport, "grover": grover, "qml": qml, "shor15": shor15}
QEC = {"bitflip": bitflip, "phaseflip": phaseflip, "shor9": shor9, "steane7": steane7, "surface3": surface3}
