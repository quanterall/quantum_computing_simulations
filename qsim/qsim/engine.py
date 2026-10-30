"""A state-vector quantum simulator on NumPy.

The state of n qubits is a complex array of 2**n amplitudes. Basis index i holds qubit q in bit q of i
(qubit 0 is the least significant bit), and a measured outcome is written with qubit 0 as the RIGHTMOST
character -- the same convention as the BEAM simulator it replaces, so every demo's counts read alike.

A gate touches the whole array at once: the state is viewed as an n-dimensional 2x2x...x2 tensor (qubit q
is axis n-1-q) and a single-qubit gate is one tensordot on its axis; CX and CZ act on slices. That keeps a
gate O(2**n) instead of the O(4**n) of copying a tuple per amplitude, which is what limited the BEAM
version to about nine qubits.
"""
from __future__ import annotations

import time

import numpy as np

MAX_QUBITS = 24          # 2**24 complex128 = 256 MB: the most one request may ask for

INV_SQRT2 = 1.0 / np.sqrt(2.0)
H = np.array([[INV_SQRT2, INV_SQRT2], [INV_SQRT2, -INV_SQRT2]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def rx(theta: float) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)


def ry(theta: float) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rz(theta: float) -> np.ndarray:
    return np.array([[np.exp(-0.5j * theta), 0], [0, np.exp(0.5j * theta)]], dtype=complex)


class Circuit:
    """A state that gates are applied to in order; each method returns self, so a circuit reads as a chain."""

    def __init__(self, n: int):
        if not 1 <= n <= MAX_QUBITS:
            raise ValueError(f"{n} qubits: between 1 and {MAX_QUBITS}")
        self.n = n
        self.psi = np.zeros((2,) * n, dtype=complex)
        self.psi[(0,) * n] = 1.0

    def _ax(self, q: int) -> int:
        if not 0 <= q < self.n:
            raise ValueError(f"qubit {q} outside 0..{self.n - 1}")
        return self.n - 1 - q

    def gate(self, u: np.ndarray, q: int) -> "Circuit":
        ax = self._ax(q)
        self.psi = np.moveaxis(np.tensordot(u, self.psi, axes=([1], [ax])), 0, ax)
        return self

    def h(self, q): return self.gate(H, q)
    def x(self, q): return self.gate(X, q)
    def y(self, q): return self.gate(Y, q)
    def z(self, q): return self.gate(Z, q)
    def rx(self, q, theta): return self.gate(rx(theta), q)
    def ry(self, q, theta): return self.gate(ry(theta), q)
    def rz(self, q, theta): return self.gate(rz(theta), q)

    def cx(self, ctrl: int, tgt: int) -> "Circuit":
        """flip the target where the control is 1"""
        if ctrl == tgt:
            raise ValueError("control and target must differ")
        c, t = self._ax(ctrl), self._ax(tgt)
        idx = [slice(None)] * self.n
        idx[c] = 1
        sub = self.psi[tuple(idx)]                       # the control-1 half; the target axis shifts if after c
        self.psi[tuple(idx)] = np.flip(sub, axis=t - (1 if t > c else 0))
        return self

    def cz(self, ctrl: int, tgt: int) -> "Circuit":
        """a phase of -1 where both qubits are 1"""
        if ctrl == tgt:
            raise ValueError("control and target must differ")
        idx = [slice(None)] * self.n
        idx[self._ax(ctrl)] = 1
        idx[self._ax(tgt)] = 1
        self.psi[tuple(idx)] *= -1
        return self

    def swap(self, a: int, b: int) -> "Circuit":
        """exchange qubits a and b"""
        return self.cx(a, b).cx(b, a).cx(a, b)

    def cswap(self, c: int, a: int, b: int) -> "Circuit":
        """Fredkin: exchange a and b where the control is 1"""
        idx = [slice(None)] * self.n
        idx[self._ax(c)] = 1
        sub = self.psi[tuple(idx)]                        # the control-1 half; the other axes shift if after c
        ax = lambda q: self._ax(q) - (1 if self._ax(q) > self._ax(c) else 0)
        self.psi[tuple(idx)] = np.swapaxes(sub, ax(a), ax(b)).copy()
        return self

    def cp(self, c: int, t: int, theta: float) -> "Circuit":
        """controlled phase: e^{iθ} where both qubits are 1"""
        idx = [slice(None)] * self.n
        idx[self._ax(c)] = 1
        idx[self._ax(t)] = 1
        self.psi[tuple(idx)] *= np.exp(1j * theta)
        return self

    def probabilities(self) -> np.ndarray:
        p = np.abs(self.psi.reshape(-1)) ** 2
        return p / p.sum()

    def sample(self, shots: int, rng: np.random.Generator | None = None) -> dict[str, int]:
        """measure every qubit `shots` times: {bitstring (qubit 0 rightmost): count}"""
        rng = rng or np.random.default_rng()
        hits = np.bincount(rng.choice(2 ** self.n, size=shots, p=self.probabilities()), minlength=2 ** self.n)
        return {format(i, f"0{self.n}b"): int(c) for i, c in enumerate(hits) if c}


    def measure_reset(self, q: int, rng: np.random.Generator | None = None) -> int:
        """measure qubit q, collapse the state onto the outcome, and put the qubit back to |0⟩ (an ancilla reused)"""
        rng = rng or np.random.default_rng()
        ax = self._ax(q)
        one = [slice(None)] * self.n
        one[ax] = 1
        p1 = float(np.sum(np.abs(self.psi[tuple(one)]) ** 2) / np.sum(np.abs(self.psi) ** 2))
        bit = int(rng.random() < p1)
        keep = [slice(None)] * self.n
        keep[ax] = bit
        kept = self.psi[tuple(keep)] / np.sqrt(p1 if bit else 1.0 - p1)
        self.psi = np.zeros_like(self.psi)
        zero = [slice(None)] * self.n
        zero[ax] = 0
        self.psi[tuple(zero)] = kept                       # the outcome's half, moved to the qubit's |0⟩
        return bit

    def sample_low(self, k: int, shots: int, rng: np.random.Generator | None = None) -> dict[str, int]:
        """measure the first k qubits `shots` times, ignoring the others: {bitstring (qubit 0 rightmost): count}"""
        rng = rng or np.random.default_rng()
        hits = np.bincount(rng.choice(2 ** self.n, size=shots, p=self.probabilities()) % (2 ** k), minlength=2 ** k)
        return {format(i, f"0{k}b"): int(c) for i, c in enumerate(hits) if c}


def run(n: int, shots: int, build, rng: np.random.Generator | None = None) -> tuple[dict[str, int], float]:
    """build(Circuit) applies the gates; returns (counts, milliseconds) for building and sampling"""
    t0 = time.perf_counter()
    c = Circuit(n)
    build(c)
    counts = c.sample(shots, rng)
    return counts, (time.perf_counter() - t0) * 1000.0
