"""Linear optics with photons, in Fock space (the quanterall tenant's photonic module; the same module serves
quanterall.com's checks from ~/quanterall/qsim). Conventions: docs/optics_conventions.md in ~/quanterall.

A state is a sparse map from occupation tuples to complex amplitudes. An interferometer is a unitary m × m matrix
M acting on creation operators, a_k† → Σ_j M[j, k] a_j†. Two independent ways to apply it, which must agree:

  * apply(M)       — expand the creation operators photon by photon (a polynomial in the a_j†);
  * apply_perm(M)  — every output amplitude as a permanent, Perm(M[s, t]) / √(Π s! Π t!), by Ryser's formula.
"""
from __future__ import annotations

import math
from itertools import combinations

import numpy as np


class FockState:
    def __init__(self, modes: int, amps: dict | None = None):
        self.m = modes
        self.amps: dict[tuple, complex] = amps if amps is not None else {(0,) * modes: 1.0 + 0j}

    @classmethod
    def basis(cls, occ) -> "FockState":
        return cls(len(occ), {tuple(occ): 1.0 + 0j})

    def copy(self) -> "FockState":
        return FockState(self.m, dict(self.amps))

    def norm2(self) -> float:
        return float(sum(abs(a) ** 2 for a in self.amps.values()))

    def __add__(self, other: "FockState") -> "FockState":
        out = dict(self.amps)
        for k, a in other.amps.items():
            out[k] = out.get(k, 0) + a
        return FockState(self.m, out)

    def scale(self, c: complex) -> "FockState":
        return FockState(self.m, {k: c * a for k, a in self.amps.items()})

    # ── linear optics, by expanding creation operators ──────────────────────────────────────────────────
    def apply(self, M: np.ndarray, tol: float = 1e-14) -> "FockState":
        out: dict[tuple, complex] = {}
        for t, a in self.amps.items():
            # Π_k (a_k†)^t_k / √(t_k!) |vac⟩ with a_k† → Σ_j M[j,k] a_j†: grow the polynomial photon by photon,
            # keeping monomial coefficients; a monomial Π (a_j†)^s_j |vac⟩ is √(Π s_j!) |s⟩
            poly = {(0,) * self.m: complex(a) / math.sqrt(math.prod(math.factorial(x) for x in t))}
            for k, tk in enumerate(t):
                for _ in range(tk):
                    nxt: dict[tuple, complex] = {}
                    for mono, c in poly.items():
                        for j in range(self.m):
                            if M[j, k] != 0:
                                s = list(mono); s[j] += 1; s = tuple(s)
                                nxt[s] = nxt.get(s, 0) + c * M[j, k]
                    poly = nxt
            for s, c in poly.items():
                out[s] = out.get(s, 0) + c * math.sqrt(math.prod(math.factorial(x) for x in s))
        return FockState(self.m, {k: v for k, v in out.items() if abs(v) > tol})

    # ── linear optics, by permanents (Ryser) ─────────────────────────────────────────────────────────────
    def apply_perm(self, M: np.ndarray, tol: float = 1e-14) -> "FockState":
        out: dict[tuple, complex] = {}
        for t, a in self.amps.items():
            n = sum(t)
            cols = [k for k, c in enumerate(t) for _ in range(c)]
            for s in occupations(self.m, n):
                rows = [j for j, c in enumerate(s) for _ in range(c)]
                v = permanent(M[np.ix_(rows, cols)]) / math.sqrt(math.prod(math.factorial(x) for x in s + t))
                if abs(v) > tol:
                    out[s] = out.get(s, 0) + a * v
        return FockState(self.m, {k: v for k, v in out.items() if abs(v) > tol})

    # ── elements ─────────────────────────────────────────────────────────────────────────────────────────
    def bs(self, i: int, j: int, theta: float) -> "FockState":
        return self.apply(beamsplitter(self.m, i, j, theta))

    def phase(self, k: int, phi: float) -> "FockState":
        return FockState(self.m, {o: a * np.exp(1j * phi * o[k]) for o, a in self.amps.items()})

    # ── measurement ──────────────────────────────────────────────────────────────────────────────────────
    def project(self, modes: list[int], counts: list[int]) -> "FockState":
        """the unnormalised conditional state of the other modes, given `counts` on `modes` (resolving detectors);
        its norm² is the outcome's probability"""
        keep = [k for k in range(self.m) if k not in modes]
        out: dict[tuple, complex] = {}
        for o, a in self.amps.items():
            if all(o[k] == c for k, c in zip(modes, counts)):
                r = tuple(o[k] for k in keep)
                out[r] = out.get(r, 0) + a
        return FockState(len(keep), out)

    def probabilities(self) -> dict[tuple, float]:
        return {k: abs(a) ** 2 for k, a in self.amps.items()}


def occupations(m: int, n: int):
    """every way to put n photons into m modes"""
    if m == 1:
        yield (n,)
        return
    for k in range(n, -1, -1):
        for rest in occupations(m - 1, n - k):
            yield (k,) + rest


def permanent(A: np.ndarray) -> complex:
    """Ryser's formula with Gray-code-free subsets; fine for n ≤ 8"""
    n = A.shape[0]
    if n == 0:
        return 1.0 + 0j
    total = 0j
    for r in range(1, n + 1):
        for S in combinations(range(n), r):
            total += (-1) ** r * np.prod(A[:, S].sum(axis=1))
    return (-1) ** n * total


def beamsplitter(m: int, i: int, j: int, theta: float) -> np.ndarray:
    """a_i† → cos θ a_i† − sin θ a_j†, a_j† → sin θ a_i† + cos θ a_j† (columns: where an input photon goes)"""
    M = np.eye(m, dtype=complex)
    c, s = math.cos(theta), math.sin(theta)
    M[i, i], M[j, i], M[i, j], M[j, j] = c, -s, s, c
    return M


def phase_matrix(m: int, k: int, phi: float) -> np.ndarray:
    M = np.eye(m, dtype=complex)
    M[k, k] = np.exp(1j * phi)
    return M


def random_unitary(m: int, rng: np.random.Generator) -> np.ndarray:
    """Haar-random, by QR of a complex Gaussian matrix"""
    z = (rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))) / math.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


# ── the KLM gates (Knill, Laflamme, Milburn, Nature 409, 46 (2001), Eq. 9) ─────────────────────────────────
S2 = math.sqrt(2)
NS_U = np.array([[1 - S2, 2 ** -0.25, math.sqrt(3 / S2 - 2)],
                 [2 ** -0.25, 0.5, 0.5 - 1 / S2],
                 [math.sqrt(3 / S2 - 2), 0.5 - 1 / S2, S2 - 0.5]])


def ns_matrix(m: int, signal: int, anc: tuple[int, int], angles=(22.5, 65.5302, -22.5), phase_pi=True, use_klm=False) -> np.ndarray:
    """the NS gate on mode `signal` with ancilla modes `anc` (to hold |1, 0⟩): KLM's unitary, or the beamsplitter
    form B(a1,a2, θ1) B(sig,a1, θ2) B(a1,a2, θ3) with a π phase on the signal (the review's Fig. 8, θ3 = −22.5°)"""
    idx = [signal, anc[0], anc[1]]
    if use_klm:
        U3 = NS_U
    else:
        d = math.pi / 180
        U3 = beamsplitter(3, 1, 2, angles[2] * d) @ beamsplitter(3, 0, 1, angles[1] * d) @ beamsplitter(3, 1, 2, angles[0] * d)
        if phase_pi:
            U3 = U3 @ np.diag([-1, 1, 1])
    M = np.eye(m, dtype=complex)
    M[np.ix_(idx, idx)] = U3
    return M
