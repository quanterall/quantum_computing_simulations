"""Phases 1–2: dual-rail gates against the qubit engine, the HOM dip, the KLM NS gate and CZ, with power checks.
    ~/quanterall/.venv-tools/bin/python tools/optics/klm_checks.py site/data"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim.fock import FockState, beamsplitter, ns_matrix  # noqa: E402
from qsim.engine import ry, rz  # noqa: E402

OUT = Path(sys.argv[1])
rng = np.random.default_rng(20261029)
res = {}

# ── dual-rail single-qubit gates against the qubit engine: B(θ) = Ry(−2θ), a phase on the |1⟩ rail = Rz up to phase
worst = 0.0
for _ in range(50):
    th, ph = rng.uniform(-math.pi, math.pi, 2)
    a, b = rng.normal(size=2) + 1j * rng.normal(size=2); nrm = math.hypot(abs(a), abs(b)); a, b = a / nrm, b / nrm
    st = (FockState.basis((1, 0)).scale(a) + FockState.basis((0, 1)).scale(b)).bs(0, 1, th).phase(1, ph)
    out = np.array([st.amps.get((1, 0), 0), st.amps.get((0, 1), 0)])
    ref = (np.exp(1j * ph / 2) * rz(ph)) @ ry(-2 * th) @ np.array([a, b])
    worst = max(worst, float(np.abs(out - ref).max()))
res["dual_rail_vs_qubit_engine"] = {"cases": 50, "max_diff": worst}
print("dual-rail gates vs the qubit engine, 50 random cases: max difference", f"{worst:.1e}")

# ── HOM: photons in ports a, b, each with an internal mode (0, 1); photon b overlaps photon a by o = e^{−γ|τ|/2}
#    (exponential pulses, Eq. 43); modes (a0, a1, b0, b1); a balanced beamsplitter on a↔b for each internal mode
rows = []
for gt in np.linspace(-4, 4, 81):
    o = math.exp(-abs(gt) / 2)
    st = FockState.basis((1, 0, 1, 0)).scale(o) + FockState.basis((1, 0, 0, 1)).scale(math.sqrt(1 - o * o))
    st = st.bs(0, 2, math.pi / 4).bs(1, 3, math.pi / 4)
    coinc = sum(p for k, p in st.probabilities().items() if k[0] + k[1] == 1 and k[2] + k[3] == 1)
    rows.append({"gamma_tau": round(float(gt), 3), "coincidence": coinc, "formula": 0.5 * (1 - math.exp(-abs(gt)))})
res["hom"] = {"points": len(rows), "max_diff_from_formula": max(abs(r["coincidence"] - r["formula"]) for r in rows),
              "at_zero_delay": rows[40]["coincidence"]}
print("HOM dip, 81 delays: max difference from ½(1 − e^{−γ|τ|})", f"{res['hom']['max_diff_from_formula']:.1e}", "· at τ = 0:", f"{rows[40]['coincidence']:.1e}")
with open(OUT / "hom-dip.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ── NS: mode 0 the signal (α0|0⟩ + α1|1⟩ + α2|2⟩), modes 1, 2 the ancilla |1, 0⟩; herald (1, 0)
def ns_run(alpha, **kw):
    st = sum((FockState.basis((n, 1, 0)).scale(alpha[n]) for n in range(1, 3)), FockState.basis((0, 1, 0)).scale(alpha[0]))
    h = st.apply(ns_matrix(3, 0, (1, 2), **kw)).project([1, 2], [1, 0])
    return np.array([h.amps.get((n,), 0) for n in range(3)])

worst_f, probs = 1.0, []
for _ in range(50):
    al = rng.normal(size=3) + 1j * rng.normal(size=3); al /= np.linalg.norm(al)
    for kw in ({"use_klm": True}, {}):
        out = ns_run(al, **kw)
        p = float(np.vdot(out, out).real); probs.append(p)
        target = al * np.array([1, 1, -1])
        worst_f = min(worst_f, abs(np.vdot(target, out)) ** 2 / p)
res["ns"] = {"cases": 50, "both_forms": True, "success_min": min(probs), "success_max": max(probs), "worst_fidelity": worst_f}
print(f"NS gate, 50 random inputs, KLM's unitary and the beamsplitter form: success {min(probs):.12f}–{max(probs):.12f}, worst fidelity 1 − {1 - worst_f:.1e}")

# ── CZ: qubit 1 on modes (0, 1), qubit 2 on (2, 3) (|1⟩_L = the photon in the second mode), NS ancillas (4, 5), (6, 7);
#    a balanced beamsplitter on the |1⟩_L modes 1 and 3, NS on each, the inverse beamsplitter; herald (1, 0, 1, 0)
def cz_map(**kw):
    K = np.zeros((4, 4), dtype=complex)
    leak = 0.0
    for col, (q1, q2) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
        occ = [1 - q1, q1, 1 - q2, q2, 1, 0, 1, 0]
        st = FockState.basis(tuple(occ)).bs(1, 3, math.pi / 4)
        st = st.apply(ns_matrix(8, 1, (4, 5), **kw)).apply(ns_matrix(8, 3, (6, 7), **kw)).bs(1, 3, -math.pi / 4)
        h = st.project([4, 5, 6, 7], [1, 0, 1, 0])
        for row, (r1, r2) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
            K[row, col] = h.amps.get((1 - r1, r1, 1 - r2, r2), 0)
        leak += h.norm2() - float(np.sum(np.abs(K[:, col]) ** 2))
    return K, leak

def proc_fid(K, T):
    return abs(np.trace(T.conj().T @ K)) ** 2 / (4 * np.trace(K.conj().T @ K).real)

CZ = np.diag([1, 1, 1, -1])
K, leak = cz_map()
res["cz"] = {"process_fidelity": proc_fid(K, CZ), "success": float(np.trace(K.conj().T @ K).real / 4),
             "leakage_out_of_dual_rail": leak, "heralded_map_over_quarter": np.round(K / 0.25, 9).real.tolist()}
print(f"CZ from two NS gates: process fidelity 1 − {1 - res['cz']['process_fidelity']:.1e}, success {res['cz']['success']:.12f} (1/16 = {1/16})")

# ── power checks
d = 22.5
K1, _ = cz_map(angles=(22.5, 65.5302 + 0.5, -22.5))
K2, _ = cz_map(phase_pi=False)
res["power"] = {"theta2_plus_half_degree": {"process_fidelity_to_CZ": proc_fid(K1, CZ)},
                "without_pi_phase": {"process_fidelity_to_CZ": proc_fid(K2, CZ), "process_fidelity_to_identity": proc_fid(K2, np.eye(4)),
                                     "success": float(np.trace(K2.conj().T @ K2).real / 4)}}
print("power: θ2 + 0.5° → fidelity to CZ", f"{proc_fid(K1, CZ):.4f}", "· no π phase → fidelity to CZ", f"{proc_fid(K2, CZ):.4f}",
      "to identity", f"{proc_fid(K2, np.eye(4)):.4f}", "success", f"{res['power']['without_pi_phase']['success']:.4f}")
json.dump(res, open(OUT / "klm-checks.json", "w"), indent=1)
