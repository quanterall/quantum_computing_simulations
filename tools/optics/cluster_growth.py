"""Phase 3b: the type-I fusion graph rule on state vectors, and the cost of growing a linear cluster.

Graph rule (Browne & Rudolph 2005): fusing the end qubits of linear clusters of n and m qubits gives, on success, a
linear cluster of n + m − 1 (the "−" outcome needs a Z on the fused qubit); failure measures both in Z, removing
them. Checked exactly for n, m ≤ 5 with the fusion Kraus operators found at the photon level (fusion_checks.py).

Growth: a long chain grows by fusing pieces — k-qubit linear clusters — onto its end. A piece is itself built from
Bell pairs (2-qubit clusters) by type-I fusion. On a failure at the chain, the chain loses its end qubit and the
piece its first; the rest of the piece (k − 1 qubits) can be reused or thrown away. Monte Carlo of the Bell pairs
used per qubit the chain gains, for k = 2 … 8, with and without reuse.
    ~/quanterall/.venv-tools/bin/python tools/optics/cluster_growth.py site/data"""
import json
import math
import sys
from functools import reduce
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1])
res = {}
I2, Z = np.eye(2), np.diag([1.0, -1.0])
X = np.array([[0, 1], [1, 0.0]])


def kron(*ops):
    return reduce(np.kron, ops)


def cluster(n):
    """the linear cluster state on n qubits, qubit 0 the leftmost tensor factor"""
    psi = np.ones(2 ** n) / 2 ** (n / 2)
    for i in range(n - 1):
        for idx in range(2 ** n):
            if (idx >> (n - 1 - i)) & 1 and (idx >> (n - 2 - i)) & 1:
                psi[idx] *= -1
    return psi


def apply_two_to_one(psi, n, i, K):
    """apply a 2-qubit → 1-qubit operator K (2 × 4) to adjacent qubits i, i+1 of an n-qubit state"""
    t = psi.reshape([2] * n)
    t = np.moveaxis(t, [i, i + 1], [0, 1]).reshape(4, -1)
    r = (K @ t).reshape([2] + [2] * (n - 2))
    return np.moveaxis(r, 0, i).reshape(-1)


Kp = np.array([[1, 0, 0, 0], [0, 0, 0, 1.0]]) / math.sqrt(2)
Km = np.array([[1, 0, 0, 0], [0, 0, 0, -1.0]]) / math.sqrt(2)
worst = 1.0
for n in range(2, 6):
    for m in range(2, 6):
        psi = np.kron(cluster(n), cluster(m))
        for K, corr in ((Kp, I2), (Km, Z)):
            out = apply_two_to_one(psi, n + m, n - 1, K)
            p = np.vdot(out, out).real
            out = kron(*([I2] * (n - 1) + [corr] + [I2] * (m - 1))) @ (out / math.sqrt(p))
            worst = min(worst, abs(np.vdot(cluster(n + m - 1), out)) ** 2)
res["graph_rule"] = {"cases": "n, m = 2..5, both outcomes", "worst_fidelity_to_cluster_n_plus_m_minus_1": worst,
                     "success_probability_each": 0.5}
print(f"graph rule: fusing clusters of n and m (2–5) gives the cluster of n + m − 1 with fidelity ≥ 1 − {1 - worst:.1e}")


# ── growth Monte Carlo ──────────────────────────────────────────────────────────────────────────────────────
def build_piece(k, rng, start=2, fresh=True):
    """grow a piece to k qubits from `start`, fusing Bell pairs onto its end; returns Bell pairs used (a fresh start
    is one pair; a reused leftover costs nothing)"""
    used, a = (1 if fresh else 0), start
    while a < k:
        used += 1                                            # one more Bell pair
        if rng.random() < 0.5:
            a += 1                                           # success: a + 2 − 1
        else:
            a -= 1                                           # failure: both ends measured; the pair is lost
            if a < 2:
                used += 1; a = 2                             # start again from a fresh pair
    return used


def grow(k, reuse, attempts, rng):
    pairs, gain, spare = 0, 0, None
    for _ in range(attempts):
        if spare is not None and spare >= 2:
            pairs += build_piece(k, rng, start=spare, fresh=False)
        else:
            pairs += build_piece(k, rng)
        spare = None
        if rng.random() < 0.5:
            gain += k - 1                                    # success: L + k − 1
        else:
            gain -= 1                                        # failure: the chain loses its end qubit
            if reuse:
                spare = k - 1                                # the piece loses its first qubit; the rest is kept
    return pairs, gain


rng = np.random.default_rng(20261029)
rows = []
for k in range(2, 9):
    for reuse in (False, True):
        per = []
        for _ in range(20):
            p, g = grow(k, reuse, 20000, rng)
            per.append(p / g if g > 0 else float("inf"))
        finite = [x for x in per if math.isfinite(x)]
        mean = float(np.mean(finite)) if len(finite) == len(per) else float("inf")
        se = float(np.std(finite, ddof=1) / math.sqrt(len(finite))) if len(finite) == len(per) else float("nan")
        rows.append({"piece_qubits": k, "reuse_leftover": reuse, "pairs_per_added_qubit": mean, "std_error": se,
                     "expected_gain_per_attempt": (k - 1) * 0.5 - 0.5})
        print(f"k = {k}, reuse {reuse!s:5}: {mean:.3f} ± {se:.3f} Bell pairs per added qubit (gain per attempt {(k - 2) / 2:+.1f})")
res["growth"] = {"method": "chain grown by type-I fusion of k-qubit pieces onto its end; each piece built from Bell pairs by type-I fusion; on a failure at the chain the chain loses one qubit and the piece its first; the leftover reused or not; 20 runs × 20,000 attempts",
                 "rows": rows}
json.dump(res, open(OUT / "cluster-growth.json", "w"), indent=1)
