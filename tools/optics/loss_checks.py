"""Phase 5: detectors and loss (Kok & Lovett 2010; the review's Eqs. 80–81).
A detector of efficiency η is a beamsplitter of transmissivity η onto a fresh vacuum (environment) mode, which is
then traced out, followed by an ideal detector. Checks:
  1. the POVMs: for |n⟩, k counts with probability C(n,k) ηᵏ (1−η)ⁿ⁻ᵏ (resolving); no click with (1−η)ⁿ (bucket);
  2. the heralded NS and CZ gates against η, with resolving and with bucket detectors on the heralds;
  3. loss on a dual-rail qubit: ρ → ηρ + (1−η)|vac⟩⟨vac| exactly — an erasure, heralded by the missing photon.
    ~/quanterall/.venv-tools/bin/python tools/optics/loss_checks.py site/data"""
import csv
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim.fock import FockState, ns_matrix  # noqa: E402

OUT = Path(sys.argv[1])
res = {}


def lossy(st: FockState, modes, eta) -> FockState:
    """each mode in `modes` sends a fraction 1 − η of its photons into a new environment mode (appended)"""
    m0 = st.m
    big = FockState(m0 + len(modes), {k + (0,) * len(modes): a for k, a in st.amps.items()})
    t = math.acos(math.sqrt(eta))                          # cos²t = η transmitted
    for i, k in enumerate(modes):
        big = big.bs(k, m0 + i, t)
    return big


def reduced(st: FockState, keep):
    """density matrix over the occupations of `keep`, tracing the other modes"""
    rest = [k for k in range(st.m) if k not in keep]
    groups = {}
    for o, a in st.amps.items():
        groups.setdefault(tuple(o[k] for k in rest), []).append((tuple(o[k] for k in keep), a))
    rho = {}
    for _, items in groups.items():
        for (k1, a1), (k2, a2) in itertools.product(items, items):
            rho[(k1, k2)] = rho.get((k1, k2), 0) + a1 * np.conj(a2)
    return rho


# ── 1. POVMs ───────────────────────────────────────────────────────────────────────────────────────────────
worst = 0.0
for eta in (0.3, 0.7, 0.95):
    for n in range(0, 5):
        st = lossy(FockState.basis((n,)), [0], eta)
        for k in range(0, n + 1):
            p = st.project([0], [k]).norm2()
            worst = max(worst, abs(p - math.comb(n, k) * eta ** k * (1 - eta) ** (n - k)))
        worst = max(worst, abs(st.project([0], [0]).norm2() - (1 - eta) ** n))
res["povm"] = {"cases": "η = 0.3, 0.7, 0.95; n = 0..4", "max_diff_from_eqs_80_81": worst}
print("detector POVMs vs Eqs. 80–81:", f"{worst:.1e}")

# ── 2. heralded gates against detector efficiency ──────────────────────────────────────────────────────────
def herald_patterns(n_modes, want, bucket):
    """detected patterns that count as the herald `want`, for resolving or bucket detectors"""
    if not bucket:
        return [tuple(want)]
    # a bucket detector only says click / no click: any count ≥ 1 where want has 1
    return [p for p in itertools.product(range(5), repeat=n_modes) if all((x >= 1) == (w >= 1) for x, w in zip(p, want))]   # up to 4 photons a mode


def ns_vs_eta(eta, bucket):
    """input α(|0⟩+|1⟩+|2⟩)/√3 on mode 0, ancilla |1,0⟩; heralds on modes 1, 2 through efficiency η.
    Returns the herald probability and the fidelity of the heralded mode-0 state with the ideal NS output."""
    al = np.ones(3) / math.sqrt(3)
    st = FockState.basis((0, 1, 0)).scale(al[0]) + FockState.basis((1, 1, 0)).scale(al[1]) + FockState.basis((2, 1, 0)).scale(al[2])
    st = lossy(st.apply(ns_matrix(3, 0, (1, 2))), [1, 2], eta)              # modes 3, 4: the environment
    target = {(0,): al[0], (1,): al[1], (2,): -al[2]}
    tn = math.sqrt(sum(abs(v) ** 2 for v in target.values()))
    p_h, f_num = 0.0, 0.0
    for pat in herald_patterns(2, (1, 0), bucket):
        h = st.project([1, 2], list(pat))                                   # modes 0, env3, env4
        rho = reduced(h, [0])
        p = sum(rho.get(((k,), (k,)), 0).real for k in range(4))
        p_h += p
        f_num += sum((np.conj(target.get(k1, 0)) * rho.get((k1, k2), 0) * target.get(k2, 0)).real for (k1, k2) in rho) / tn ** 2
    return p_h, (f_num / p_h if p_h else 0.0)


def cz_vs_eta(eta, bucket):
    """input |+⟩|+⟩ (dual rail), heralds on the four ancilla modes through efficiency η; the heralded output's
    fidelity with CZ|++⟩ and the herald probability"""
    st = FockState(8, {})
    for q1, q2 in itertools.product(range(2), repeat=2):
        st = st + FockState.basis((1 - q1, q1, 1 - q2, q2, 1, 0, 1, 0)).scale(0.5)
    st = st.bs(1, 3, math.pi / 4).apply(ns_matrix(8, 1, (4, 5))).apply(ns_matrix(8, 3, (6, 7))).bs(1, 3, -math.pi / 4)
    st = lossy(st, [4, 5, 6, 7], eta)                                       # env 8..11
    target = {}
    for q1, q2 in itertools.product(range(2), repeat=2):
        target[(1 - q1, q1, 1 - q2, q2)] = 0.5 * (-1 if q1 and q2 else 1)
    p_h, f_num = 0.0, 0.0
    for pat in herald_patterns(4, (1, 0, 1, 0), bucket):
        h = st.project([4, 5, 6, 7], list(pat))                             # modes 0..3, env
        rho = reduced(h, [0, 1, 2, 3])
        p = sum(v.real for (k1, k2), v in rho.items() if k1 == k2)
        p_h += p
        f_num += sum((np.conj(target.get(k1, 0)) * v * target.get(k2, 0)).real for (k1, k2), v in rho.items())
    return p_h, (f_num / p_h if p_h else 0.0)


rows = []
for eta in (1.0, 0.99, 0.97, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5):
    row = {"eta": eta}
    for bucket in (False, True):
        tag = "bucket" if bucket else "resolving"
        row[f"ns_herald_{tag}"], row[f"ns_fidelity_{tag}"] = ns_vs_eta(eta, bucket)
        row[f"cz_herald_{tag}"], row[f"cz_fidelity_{tag}"] = cz_vs_eta(eta, bucket)
    rows.append(row)
    print(f"η {eta:.2f}: NS fidelity {row['ns_fidelity_resolving']:.4f} (bucket {row['ns_fidelity_bucket']:.4f}), CZ fidelity {row['cz_fidelity_resolving']:.4f} (bucket {row['cz_fidelity_bucket']:.4f}), CZ herald {row['cz_herald_resolving']:.4f}")
with open(OUT / "heralded-vs-efficiency.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
res["heralded_vs_eta"] = rows

# ── 3. loss on a dual-rail qubit is an erasure ─────────────────────────────────────────────────────────────
rng = np.random.default_rng(3)
worst = 0.0
for _ in range(30):
    a, b = rng.normal(size=2) + 1j * rng.normal(size=2); nrm = math.hypot(abs(a), abs(b)); a, b = a / nrm, b / nrm
    eta = rng.uniform(0.2, 0.99)
    st = lossy(FockState.basis((1, 0)).scale(a) + FockState.basis((0, 1)).scale(b), [0, 1], eta)
    rho = reduced(st, [0, 1])
    ideal = {((1, 0), (1, 0)): eta * abs(a) ** 2, ((1, 0), (0, 1)): eta * a * np.conj(b),
             ((0, 1), (1, 0)): eta * b * np.conj(a), ((0, 1), (0, 1)): eta * abs(b) ** 2, ((0, 0), (0, 0)): 1 - eta}
    keys = set(rho) | set(ideal)
    worst = max(worst, max(abs(rho.get(k, 0) - ideal.get(k, 0)) for k in keys))
res["erasure_channel"] = {"cases": 30, "max_diff_from_eta_rho_plus_vacuum": worst}
print("dual-rail loss = ηρ + (1−η)|vac⟩⟨vac|:", f"{worst:.1e}")
json.dump(res, open(OUT / "loss-checks.json", "w"), indent=1)
