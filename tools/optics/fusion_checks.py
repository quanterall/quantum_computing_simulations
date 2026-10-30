"""Phase 3a: type-I and type-II fusion (Browne & Rudolph, PRL 95, 010501, 2005) at the photon level.
Modes per photon: (H, V) of its spatial port. A polarising beamsplitter (PBS) transmits H and reflects V: the
V modes swap ports. A 45° polarisation rotation is a balanced beamsplitter between a port's H and V modes.
    ~/quanterall/.venv-tools/bin/python tools/optics/fusion_checks.py site/data"""
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim.fock import FockState  # noqa: E402

OUT = Path(sys.argv[1])
res = {}
Q = [(1, 0), (0, 1)]                                  # a polarisation qubit: |0⟩ = H, |1⟩ = V


def two_photon(psi):
    """ψ over |q1 q2⟩ (q1 the high bit) → photon 1 in port 1 (modes 0 H, 1 V), photon 2 in port 2 (modes 2 H, 3 V)"""
    st = FockState(4, {})
    for i, (a, b) in enumerate(itertools.product(range(2), repeat=2)):
        if abs(psi[i]) > 0:
            st = st + FockState.basis(Q[a] + Q[b]).scale(psi[i])
    return st


def pbs(st):
    """V modes swap ports: (1H, 1V, 2H, 2V) → (1H, 2V, 2H, 1V)"""
    M = np.zeros((4, 4), dtype=complex)
    M[0, 0] = M[2, 2] = 1
    M[3, 1] = M[1, 3] = 1
    return st.apply(M)


# ── type I: PBS, 45° on port 2, a polarisation-resolving detector on port 2; success = one photon there ──────
def type1(psi):
    st = pbs(two_photon(psi)).bs(2, 3, math.pi / 4)
    outs = {}
    for name, counts in (("+", (1, 0)), ("−", (0, 1))):
        h = st.project([2, 3], list(counts))
        outs[name] = np.array([h.amps.get((1, 0), 0), h.amps.get((0, 1), 0)])        # port 1's photon: the fused qubit
    fail = sum(st.project([2, 3], list(c)).norm2() for c in ((0, 0), (2, 0), (1, 1), (0, 2)))
    return outs, fail


# Kraus operators, column by column from the basis inputs
K = {"+": np.zeros((2, 4), dtype=complex), "−": np.zeros((2, 4), dtype=complex)}
for col in range(4):
    e = np.zeros(4); e[col] = 1
    outs, _ = type1(e)
    for s in K:
        K[s][:, col] = outs[s]
ideal = {"+": np.array([[1, 0, 0, 0], [0, 0, 0, 1]]) / math.sqrt(2) / math.sqrt(2) * math.sqrt(2),
         "−": np.array([[1, 0, 0, 0], [0, 0, 0, -1]]) / math.sqrt(2)}
ideal["+"] = np.array([[1, 0, 0, 0], [0, 0, 0, 1]]) / math.sqrt(2)


def same_up_to_phase(A, B):
    a, b = A.flatten(), B.flatten()
    ph = np.vdot(b, a) / max(abs(np.vdot(b, a)), 1e-30)
    return float(np.abs(a - ph * b).max())


res["type1"] = {
    "kraus_plus_vs_|0⟩⟨00|+|1⟩⟨11| over √2": same_up_to_phase(K["+"], ideal["+"]),
    "kraus_minus_vs_|0⟩⟨00|−|1⟩⟨11| over √2": same_up_to_phase(K["−"], ideal["−"]),
}
rng = np.random.default_rng(5)
succ, failz = [], 0.0
for _ in range(200):
    psi = rng.normal(size=4) + 1j * rng.normal(size=4); psi /= np.linalg.norm(psi)
    outs, fail = type1(psi)
    p = sum(float(np.vdot(o, o).real) for o in outs.values())
    # success probability equals the input's weight on HH and VV
    succ.append(abs(p - (abs(psi[0]) ** 2 + abs(psi[3]) ** 2)))
res["type1"]["success_equals_HH_plus_VV_weight_max_diff"] = max(succ)
# on a product of |+⟩ states (cluster qubits), success is 1/2
plus = np.ones(4) / 2
o, f = type1(plus)
res["type1"]["success_on_plus_plus"] = sum(float(np.vdot(x, x).real) for x in o.values())
# failure: zero or two photons on port 2 ⇒ HV or VH — both qubits measured in Z (their values known from the count)
st = pbs(two_photon(np.array([0, 1, 0, 0]))).bs(2, 3, math.pi / 4)                    # |0 1⟩ = H V
res["type1"]["HV_fails_with_port2_empty"] = st.project([2, 3], [0, 0]).norm2()
st = pbs(two_photon(np.array([0, 0, 1, 0]))).bs(2, 3, math.pi / 4)                    # |1 0⟩ = V H
res["type1"]["VH_fails_with_two_photons_on_port2"] = sum(st.project([2, 3], list(c)).norm2() for c in ((2, 0), (1, 1), (0, 2)))
print("type I:", json.dumps(res["type1"], ensure_ascii=False))


# ── type II: 45° on both inputs, PBS, 45° on both outputs, polarisation-resolving detectors on both ports;
#    success = one photon on each port; the measured pair projects onto a Bell state (XX and ZZ known) ─────────
def type2_povm():
    els = {}
    for col in range(4):
        e = np.zeros(4); e[col] = 1
        st = two_photon(e).bs(0, 1, math.pi / 4).bs(2, 3, math.pi / 4)
        st = pbs(st).bs(0, 1, math.pi / 4).bs(2, 3, math.pi / 4)
        for c1 in ((1, 0), (0, 1)):
            for c2 in ((1, 0), (0, 1)):
                amp = st.amps.get(c1 + c2, 0)
                els.setdefault(c1 + c2, np.zeros(4, dtype=complex))[col] = amp
    return els                                                              # each: a row vector ⟨φ| scaled


els = type2_povm()
bell = {"Φ+": np.array([1, 0, 0, 1]) / math.sqrt(2), "Φ−": np.array([1, 0, 0, -1]) / math.sqrt(2),
        "Ψ+": np.array([0, 1, 1, 0]) / math.sqrt(2), "Ψ−": np.array([0, 1, -1, 0]) / math.sqrt(2)}
t2 = {}
for k, v in els.items():
    n = np.linalg.norm(v)
    best = max(bell, key=lambda b: abs(np.vdot(bell[b], v / n)))
    t2["".join(map(str, k))] = {"bell_state": best, "overlap": float(abs(np.vdot(bell[best], v / n)) ** 2), "weight": float(n ** 2)}
res["type2"] = {"success_outcomes": t2, "success_probability_on_any_input_avg": float(sum(x["weight"] for x in t2.values()) / 4)}
print("type II:", json.dumps(res["type2"], ensure_ascii=False))
json.dump(res, open(OUT / "fusion-checks.json", "w"), indent=1, ensure_ascii=False)
