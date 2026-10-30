"""The engine against a dense reference: every gate as a full 2**n x 2**n matrix built from Kronecker
products, so the tensor axes, the bit order and the CX/CZ slicing are checked independently."""
import numpy as np
import pytest

from qsim.engine import Circuit, H, X, Y, Z, rx, rz

I2 = np.eye(2, dtype=complex)


def full_1q(u, q, n):
    """u on qubit q of n; qubit 0 is the least significant bit, i.e. the LAST Kronecker factor"""
    m = np.array([[1]], dtype=complex)
    for k in reversed(range(n)):
        m = np.kron(m, u if k == q else I2)
    return m


def full_cx(c, t, n):
    m = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(2 ** n):
        j = i ^ (1 << t) if (i >> c) & 1 else i
        m[j, i] = 1
    return m


def full_cz(c, t, n):
    return np.diag([-1 if (i >> c) & 1 and (i >> t) & 1 else 1 for i in range(2 ** n)]).astype(complex)


@pytest.mark.parametrize("seed", range(20))
def test_random_circuits_match_the_dense_reference(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(2, 7))
    c = Circuit(n)
    ref = np.zeros(2 ** n, dtype=complex); ref[0] = 1
    for _ in range(25):
        kind = rng.choice(["h", "x", "y", "z", "rx", "rz", "cx", "cz"])
        if kind in ("cx", "cz"):
            a, b = (int(v) for v in rng.choice(n, 2, replace=False))
            getattr(c, kind)(a, b)
            ref = (full_cx if kind == "cx" else full_cz)(a, b, n) @ ref
        else:
            q = int(rng.integers(n))
            if kind in ("rx", "rz"):
                th = float(rng.uniform(-np.pi, np.pi))
                getattr(c, kind)(q, th)
                ref = full_1q((rx if kind == "rx" else rz)(th), q, n) @ ref
            else:
                getattr(c, kind)(q)
                ref = full_1q({"h": H, "x": X, "y": Y, "z": Z}[kind], q, n) @ ref
    assert np.allclose(c.psi.reshape(-1), ref, atol=1e-12)
    assert np.isclose(np.sum(np.abs(ref) ** 2), 1.0)


def test_bit_order_qubit0_is_rightmost():
    c = Circuit(3).x(0)
    assert c.sample(10) == {"001": 10}
    assert Circuit(3).x(2).sample(5) == {"100": 5}


def test_bell_is_half_00_half_11():
    p = Circuit(2).h(0).cx(0, 1).probabilities()
    assert np.allclose(p, [0.5, 0, 0, 0.5])


def test_cx_both_directions():
    assert Circuit(2).x(1).cx(1, 0).sample(4) == {"11": 4}      # control above the target
    assert Circuit(2).x(0).cx(0, 1).sample(4) == {"11": 4}      # control below the target
    assert Circuit(2).cx(1, 0).sample(4) == {"00": 4}           # control 0: nothing


def test_sampling_is_unbiased():
    counts = Circuit(1).h(0).sample(100_000, np.random.default_rng(1))
    assert abs(counts["0"] / 100_000 - 0.5) < 0.01


@pytest.mark.parametrize("bad", [0, 25])
def test_qubit_count_is_bounded(bad):
    with pytest.raises(ValueError):
        Circuit(bad)


def test_gates_reject_bad_qubits():
    with pytest.raises(ValueError):
        Circuit(2).h(2)
    with pytest.raises(ValueError):
        Circuit(2).cx(1, 1)
