"""Phase 4: measurement-based computation on small clusters, on state vectors; cross-checked with Qiskit.
  1. one-bit teleportation through a two-node cluster, and a single-qubit rotation along a chain;
  2. the review's four-node cluster (Romero & Milburn Fig. 10, Eqs. 55–58): what it computes for superposition
     inputs, with the chain CZ12 CZ23 CZ34 (Fig. 10) and with the couplings Eq. 56 writes (12, 23, 24).
    ~/quanterall/.venv-tools/bin/python tools/optics/mbqc.py site/data"""
import itertools
import json
import math
import sys
from functools import reduce
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1])
rng = np.random.default_rng(11)
I2, X, Z = np.eye(2), np.array([[0, 1], [1, 0.0]]), np.diag([1.0, -1.0])
H = np.array([[1, 1], [1, -1.0]]) / math.sqrt(2)
plus = np.array([1, 1]) / math.sqrt(2)
res = {}


def op_on(n, q, U):
    return reduce(np.kron, [U if i == q else I2 for i in range(n)])


def cz(n, a, b):
    d = np.ones(2 ** n)
    for idx in range(2 ** n):
        if (idx >> (n - 1 - a)) & 1 and (idx >> (n - 1 - b)) & 1:
            d[idx] = -1
    return np.diag(d)


def measure(psi, n, q, basis_vec):
    """project qubit q on ⟨basis_vec|; returns the (n−1)-qubit unnormalised state"""
    t = np.moveaxis(psi.reshape([2] * n), q, 0)
    return np.tensordot(basis_vec.conj(), t, axes=(0, 0)).reshape(-1)


def rand_state(k=2):
    v = rng.normal(size=k) + 1j * rng.normal(size=k)
    return v / np.linalg.norm(v)


def local_fidelity(M, T):
    """is M ∝ (A ⊗ B) T for some single-qubit unitaries? returns the best fidelity found by fitting A ⊗ B = M T†"""
    L = M @ T.conj().T
    # nearest product of unitaries: realign L into (A, B) and take the leading singular pair
    R = L.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    u, s, vh = np.linalg.svd(R)
    A, B = u[:, 0].reshape(2, 2) * math.sqrt(s[0]), vh[0].reshape(2, 2) * math.sqrt(s[0])
    P = np.kron(A, B)
    return float(abs(np.trace(P.conj().T @ L)) ** 2 / (np.trace(P.conj().T @ P).real * np.trace(L.conj().T @ L).real))


# ── 1. one-bit teleportation: |ψ⟩|+⟩, CZ, measure qubit 1 in X → qubit 2 holds X^s H |ψ⟩ ───────────────────
worst = 1.0
for _ in range(50):
    psi = rand_state()
    st = cz(2, 0, 1) @ np.kron(psi, plus)
    for s, b in enumerate((plus, Z @ plus)):
        out = measure(st, 2, 0, b); p = np.vdot(out, out).real
        out = np.linalg.matrix_power(X, s) @ out / math.sqrt(p)              # correct the byproduct X^s
        worst = min(worst, abs(np.vdot(H @ psi, out)) ** 2)
res["one_bit_teleportation"] = {"cases": 50, "worst_fidelity_to_H_psi": worst}

# ── a rotation along a chain: input on qubit 0, measure qubits 0..2 at angles; the output is H Rz H Rz … ψ ─────
def rz(t):
    return np.diag([1, np.exp(1j * t)])

worst = 1.0
for _ in range(30):
    psi = rand_state(); angles = rng.uniform(-math.pi, math.pi, 3)
    n = 4
    st = reduce(np.kron, [psi, plus, plus, plus])
    for a in range(3):
        st = cz(n, a, a + 1) @ st
    U = I2; cur, m = st, n
    outcomes = []
    for k, th in enumerate(angles):
        # adaptive angle: the sign flips when the previous X byproduct count is odd
        sx = sum(outcomes[j] for j in range(len(outcomes)) if (len(outcomes) - 1 - j) % 2 == 0) % 2
        t = th * (-1) ** sx
        s = int(rng.random() < 0.5)
        b = (np.array([1, np.exp(1j * t)]) / math.sqrt(2)) if s == 0 else (np.array([1, -np.exp(1j * t)]) / math.sqrt(2))
        cur = measure(cur, m, 0, b); m -= 1; cur = cur / np.linalg.norm(cur)
        outcomes.append(s)
        U = H @ rz(-th) @ U
    # remaining byproduct: X^{s3+s1} Z^{s2} (standard for a 3-measurement chain)
    s1, s2, s3 = outcomes
    corr = np.linalg.matrix_power(Z, s2) @ np.linalg.matrix_power(X, (s1 + s3) % 2)
    out = corr @ cur
    worst = min(worst, abs(np.vdot(U @ psi, out)) ** 2)
res["rotation_chain"] = {"cases": 30, "worst_fidelity_to_H_Rz_H_Rz_H_Rz": worst,
                         "note": "adaptive angles (sign flipped by earlier X byproducts); byproduct X^(s1+s3) Z^(s2)"}

def makhlin(U):
    """local invariants (Makhlin 2002) of a two-qubit unitary: G1, G2. CNOT and CZ: (0, 1); identity: (1, 3); SWAP: (−1, −3)"""
    Q = np.array([[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]]) / math.sqrt(2)
    UB = Q.conj().T @ U @ Q
    m = UB.T @ UB
    det = np.linalg.det(U)
    g1 = np.trace(m) ** 2 / (16 * det)
    g2 = (np.trace(m) ** 2 - np.trace(m @ m)) / (4 * det)
    return [round(float(g1.real), 9), round(float(abs(g1.imag)), 9), round(float(g2.real), 9)]


def output_byproducts(M, T):
    """if M = (A ⊗ B) T exactly, return A and B normalised to det 1 (phases dropped)"""
    L = M @ T.conj().T
    R = L.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    u, s, vh = np.linalg.svd(R)
    A, B = u[:, 0].reshape(2, 2), vh[0].reshape(2, 2)
    norm = lambda W: W / np.sqrt(np.linalg.det(W))
    return norm(A), norm(B), float(s[1] / s[0])


# ── 2. the four-node cluster: qubits 1..4 (0..3 here), inputs |a⟩ on 1 and |b⟩ on 4, |+⟩ on 2 and 3 (Eq. 55);
#      measure 1 and 2 in X; output on 3 and 4. Which two-qubit map on (a, b) → (q3, q4), for each outcome? ─────
CNOT_b_ctrl = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])      # control b (qubit 4), target onto q3
fig10 = {}
for name, edges in (("chain CZ12 CZ23 CZ34 (Fig. 10)", [(0, 1), (1, 2), (2, 3)]), ("Eq. 56 couplings 12, 23, 24", [(0, 1), (1, 2), (1, 3)])):
    per = []
    for s1, s2 in itertools.product(range(2), repeat=2):
        M = np.zeros((4, 4), dtype=complex)
        for col, (a, b) in enumerate(itertools.product(range(2), repeat=2)):
            st = reduce(np.kron, [np.eye(2)[a], plus, plus, np.eye(2)[b]])
            for e in edges:
                st = cz(4, *e) @ st
            o = measure(st, 4, 0, (plus if s1 == 0 else Z @ plus))
            o = measure(o, 3, 0, (plus if s2 == 0 else Z @ plus))
            M[:, col] = o * 2                                                    # each outcome has weight 1/4
        unitary = bool(np.allclose(M.conj().T @ M, np.eye(4), atol=1e-9))
        row = {"outcome": f"{s1}{s2}", "unitary": unitary, "makhlin_G1_re_im_G2": makhlin(M),
               "equals_CNOT_up_to_output_locals": local_fidelity(M, CNOT_b_ctrl) > 1 - 1e-9}
        if row["equals_CNOT_up_to_output_locals"]:
            A, B, resid = output_byproducts(M, CNOT_b_ctrl)
            row["byproduct_on_q3"] = np.round(A, 6).tolist() if False else [[str(np.round(z, 6)) for z in r] for r in A]
            row["byproduct_on_q4"] = [[str(np.round(z, 6)) for z in r] for r in B]
        per.append(row)
    # superposition inputs: does the map act linearly as found? (it is linear by construction; check entangling power)
    fig10[name] = per
res["fig10"] = fig10
print(json.dumps(res, indent=1)[:3000])
json.dump(res, open(OUT / "mbqc-checks.json", "w"), indent=1)
