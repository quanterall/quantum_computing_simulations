"""One optical mode, continuous variables: squeezed vacuum, GKP states (the review's Eq. 38), cat states (Eqs. 40, 75)
and the homodyne-heralded comb (Eqs. 66–73) — Phase 7 of ~/quanterall/docs/photonic-plan.md.

Conventions (Eq. 4 of the review, hbar = 1): q = (a + a†)/√2, p = −i(a − a†)/√2, [q, p] = i, vacuum Var(q) = 1/2.
S(ζ) = exp[(ζ* a² − ζ a†²)/2] with ζ = r e^{−2iϕ} (Eq. 22); with ϕ = 0 it squeezes q: Var(q) = e^{−2r}/2.
GKP qubit (n = 2): |0̃⟩ has q-peaks at 2sα, |1̃⟩ at (2s + 1)α, α = √π; Z̄ = e^{iπq/α}, X̄ = e^{−ipα}.

Two independent routes, compared in tools/optics/cv_checks.py: the truncated Fock basis (matrices, Hermite functions)
and wavefunctions on a q-grid (closed forms, FFT to p). numpy and math only."""
from __future__ import annotations

import math

import numpy as np

SQRT_PI = math.sqrt(math.pi)


# --- grid -------------------------------------------------------------------------------------------------------

def grid(half_width: float = 48.0, n: int = 1 << 14) -> np.ndarray:
    return np.linspace(-half_width, half_width, n, endpoint=False)


def normalise(psi: np.ndarray, q: np.ndarray) -> np.ndarray:
    return psi / math.sqrt(float(np.sum(np.abs(psi) ** 2)) * (q[1] - q[0]))


def inner(a: np.ndarray, b: np.ndarray, q: np.ndarray) -> complex:
    return complex(np.vdot(a, b) * (q[1] - q[0]))


def to_p(psi: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ψ(p) = (2π)^{−1/2} ∫ e^{−ipq} ψ(q) dq on the FFT grid, normalised on that grid"""
    dq, n = q[1] - q[0], len(q)
    p = np.fft.fftshift(np.fft.fftfreq(n, d=dq)) * 2 * math.pi
    # ψ(p_k) = dq/√(2π) Σ_j ψ(q_j) e^{−i p_k q_j}, with q_j = q_0 + j dq
    ph =np.fft.fftshift(np.fft.fft(psi)) * np.exp(-1j * p * q[0]) * dq / math.sqrt(2 * math.pi)
    return p, ph


def translate(psi: np.ndarray, q: np.ndarray, shift: float) -> np.ndarray:
    """e^{−ip·shift}ψ: ψ(q − shift), spectrally (exact for band-limited states on the periodic grid)"""
    k = np.fft.fftfreq(len(q), d=q[1] - q[0]) * 2 * math.pi
    return np.fft.ifft(np.fft.fft(psi) * np.exp(-1j * k * shift))


# --- Fock basis -------------------------------------------------------------------------------------------------

def annihilation(cut: int) -> np.ndarray:
    return np.diag(np.sqrt(np.arange(1, cut)), 1).astype(complex)


def expm_antihermitian(g: np.ndarray) -> np.ndarray:
    """exp(G) for anti-Hermitian G, through the Hermitian iG (no scipy)"""
    w, v = np.linalg.eigh(1j * g)
    return (v * np.exp(-1j * w)) @ v.conj().T


def squeeze_op(r: float, phi: float, cut: int) -> np.ndarray:
    a = annihilation(cut)
    zeta = r * np.exp(-2j * phi)
    return expm_antihermitian((np.conj(zeta) * a @ a - zeta * a.conj().T @ a.conj().T) / 2)


def displace_op(beta: complex, cut: int) -> np.ndarray:
    a = annihilation(cut)
    return expm_antihermitian(beta * a.conj().T - np.conj(beta) * a)


def coherent_fock(beta: complex, cut: int) -> np.ndarray:
    n = np.arange(cut)
    logf = np.array([math.lgamma(k + 1) for k in n]) / 2
    with np.errstate(divide="ignore"):
        mag = np.exp(-abs(beta) ** 2 / 2 + n * (math.log(abs(beta)) if beta != 0 else 0) - logf)
    if beta == 0:
        mag = (n == 0).astype(float)
    return mag * np.exp(1j * n * np.angle(beta))


def cat_fock(beta: complex, parity: int, cut: int) -> np.ndarray:
    """N±^{−1/2}(|β⟩ ± |−β⟩), N± = 2 ± 2e^{−2|β|²}: Eq. 75's normalisation (Eq. 40's N^{−1} would not normalise)"""
    s = 1 if parity == 0 else -1
    v = coherent_fock(beta, cut) + s * coherent_fock(-beta, cut)
    return v / math.sqrt(2 + s * 2 * math.exp(-2 * abs(beta) ** 2))


def hermite_functions(cut: int, q: np.ndarray) -> np.ndarray:
    """φ_n(q) = ⟨q|n⟩, n < cut, by the stable recurrence"""
    out = np.zeros((cut, len(q)))
    out[0] = math.pi ** -0.25 * np.exp(-q * q / 2)
    if cut > 1:
        out[1] = math.sqrt(2) * q * out[0]
    for k in range(1, cut - 1):
        out[k + 1] = math.sqrt(2 / (k + 1)) * q * out[k] - math.sqrt(k / (k + 1)) * out[k - 1]
    return out


def fock_to_q(c: np.ndarray, q: np.ndarray) -> np.ndarray:
    return c @ hermite_functions(len(c), q)


def quadratures(cut: int) -> tuple[np.ndarray, np.ndarray]:
    """X1 = a + a†, X2 = −i(a − a†) — the review's Fig. 3 and Eq. 25 (vacuum variance 1)"""
    a = annihilation(cut)
    return a + a.conj().T, -1j * (a - a.conj().T)


# --- closed forms on the q-grid -----------------------------------------------------------------------------------

def squeezed_vacuum_q(q: np.ndarray, r: float, phi: float = 0.0) -> np.ndarray:
    """⟨q|S(ζ)|0⟩ ∝ exp[−(q²/2)(1 + e^{−2iϕ} tanh r)/(1 − e^{−2iϕ} tanh r)]"""
    t = np.exp(-2j * phi) * math.tanh(r)
    return normalise(np.exp(-(q * q / 2) * (1 + t) / (1 - t)), q)


def coherent_q(q: np.ndarray, beta: complex) -> np.ndarray:
    """⟨q|β⟩ = π^{−1/4} exp[−(q − q0)²/2 + i p0 q − i q0 p0/2], q0 = √2 Re β, p0 = √2 Im β"""
    q0, p0 = math.sqrt(2) * beta.real, math.sqrt(2) * beta.imag
    return math.pi ** -0.25 * np.exp(-(q - q0) ** 2 / 2 + 1j * p0 * q - 1j * q0 * p0 / 2)


def cat_q(q: np.ndarray, beta: complex, parity: int) -> np.ndarray:
    s = 1 if parity == 0 else -1
    return (coherent_q(q, beta) + s * coherent_q(q, -beta)) / math.sqrt(2 + s * 2 * math.exp(-2 * abs(beta) ** 2))


def gkp_q(q: np.ndarray, mu: int, delta2: float, peak_var: float, alpha: float = SQRT_PI,
          envelope_index=None) -> np.ndarray:
    """Eq. 38: Σ_s e^{−½Δ²[(2s + μ)α]²} T[(2s + μ)α]|ψ0⟩, ψ0 a Gaussian whose density has variance peak_var
    (a q-squeezed vacuum: peak_var = e^{−2r}/2); normalised on the grid. envelope_index lets a test use another
    weighting (Fig. 5's left panel matches e^{−Δ̃(sα)²})."""
    out = np.zeros_like(q, dtype=complex)
    smax = int(q[-1] / (2 * alpha)) + 2
    for s in range(-smax, smax + 1):
        c = (2 * s + mu) * alpha
        w = math.exp(-0.5 * delta2 * c * c) if envelope_index is None else envelope_index(s, c)
        out += w * np.exp(-(q - c) ** 2 / (4 * peak_var))
    return normalise(out, q)


def logical_q(q: np.ndarray, a0: complex, a1: complex, delta2: float, peak_var: float) -> np.ndarray:
    return normalise(a0 * gkp_q(q, 0, delta2, peak_var) + a1 * gkp_q(q, 1, delta2, peak_var), q)


# --- GKP readout under Gaussian shift noise -------------------------------------------------------------------------

def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def wrong_parity_prob(x: float, sigma: float, spacing: float, want_odd: bool) -> float:
    """P(round((x + u)/spacing) has the wrong parity), u ~ N(0, σ²); σ = 0 is the plain decision"""
    if sigma == 0:
        return float((round(x / spacing) % 2 == 1) != want_odd)
    k0 = math.floor(x / spacing)
    tot = 0.0
    span = int(8 * sigma / spacing) + 2
    for k in range(k0 - span, k0 + span + 1):
        if (k % 2 == 1) == want_odd:
            continue
        lo, hi = (k - 0.5) * spacing - x, (k + 0.5) * spacing - x
        tot += _phi(hi / sigma) - _phi(lo / sigma)
    return tot


def readout_error_exact(density: np.ndarray, x: np.ndarray, sigma: float, spacing: float, want_odd: bool) -> float:
    """∫ ρ(x) P(wrong parity | x) dx — exact up to the grid's quadrature; with σ = 0 the decision is a step, integrated
    interval by interval from the cumulative (trapezoid) distribution, so the steps fall between grid points correctly"""
    dx = x[1] - x[0]
    if sigma == 0:
        cdf = np.concatenate([[0.0], np.cumsum((density[1:] + density[:-1]) / 2) * dx])
        F = lambda v: float(np.interp(v, x, cdf))
        k_lo, k_hi = math.floor(x[0] / spacing) - 1, math.ceil(x[-1] / spacing) + 1
        return sum(F((k + 0.5) * spacing) - F((k - 0.5) * spacing) for k in range(k_lo, k_hi + 1) if (k % 2 == 1) != want_odd)
    keep = density > 1e-16 * float(density.max())          # the rest cannot move the sum
    g = np.array([wrong_parity_prob(float(v), sigma, spacing, want_odd) for v in x[keep]])
    return float(np.sum(density[keep] * g) * dx)


def readout_error_closed(sigma_eff: float, spacing: float) -> float:
    """an ideal comb (infinitely narrow peaks, flat envelope) with Gaussian shifts of σ_eff:
    Σ_k [Φ((2k + 3/2)a/σ) − Φ((2k + 1/2)a/σ)] over all k, a = spacing"""
    if sigma_eff == 0:
        return 0.0
    tot, k = 0.0, 0
    while True:
        t = (_phi((2 * k + 1.5) * spacing / sigma_eff) - _phi((2 * k + 0.5) * spacing / sigma_eff))
        t += (_phi((-2 * k - 0.5) * spacing / sigma_eff) - _phi((-2 * k - 1.5) * spacing / sigma_eff))
        tot += t
        k += 1
        if t < 1e-17 or k > 10000:
            return tot


def sample_density(density: np.ndarray, x: np.ndarray, shots: int, rng: np.random.Generator) -> np.ndarray:
    dx = x[1] - x[0]
    cdf = np.cumsum(density) * dx
    cdf /= cdf[-1]
    u = rng.random(shots)
    i = np.searchsorted(cdf, u)
    i = np.clip(i, 1, len(x) - 1)
    c0, c1 = cdf[i - 1], cdf[i]
    return x[i - 1] + (u - c0) / np.maximum(c1 - c0, 1e-300) * dx


def readout_mc(density, x, sigma, spacing, want_odd, shots, rng) -> int:
    v = sample_density(density, x, shots, rng) + rng.normal(0, sigma, shots) if sigma > 0 else sample_density(density, x, shots, rng)
    return int(np.sum((np.round(v / spacing).astype(np.int64) % 2 == 1) != want_odd))


# --- the homodyne-heralded comb (Eqs. 66–73) ------------------------------------------------------------------------

def upsilon(q: np.ndarray, y: float, beta: float, theta: float, phi: float) -> np.ndarray:
    """Eq. 69 with its overloaded r renamed: β the probe amplitude, θ the coupling; the Kraus amplitude"""
    c = theta * q - phi
    return (2 * math.pi) ** -0.25 * np.exp(1j * y * beta * np.sin(c) - (y - 2 * beta * np.cos(c)) ** 2 / 4)


def upsilon_sq_printed(q: np.ndarray, y: float, beta: float, theta: float, phi: float) -> np.ndarray:
    """Eq. 70 as printed: (2π)^{−1/2} exp[−(y − 2β cos(θq − ϕ))²] — not the square of Eq. 69"""
    return (2 * math.pi) ** -0.5 * np.exp(-(y - 2 * beta * np.cos(theta * q - phi)) ** 2)


def heralded_comb(q, psi_in, beta, theta, phi, y_grid):
    """P(y) (Eq. 73) on y_grid and a function giving the normalised output for one y (Eq. 72)"""
    dens = np.abs(psi_in) ** 2
    dq = q[1] - q[0]
    py = np.array([float(np.sum(np.abs(upsilon(q, y, beta, theta, phi)) ** 2 * dens) * dq) for y in y_grid])

    def out(y):
        return normalise(upsilon(q, y, beta, theta, phi) * psi_in, q)
    return py, out


# --- Wigner function from a wavefunction -----------------------------------------------------------------------------

def wigner(psi: np.ndarray, q: np.ndarray, qs: np.ndarray, ps: np.ndarray) -> np.ndarray:
    """W(q, p) = (1/π) ∫ ψ*(q + y) ψ(q − y) e^{2ipy} dy, by interpolation on the grid"""
    dq = q[1] - q[0]
    ys = np.arange(-len(q) // 4, len(q) // 4) * dq
    re = lambda x: np.interp(x, q, psi.real, left=0, right=0) + 1j * np.interp(x, q, psi.imag, left=0, right=0)
    W = np.zeros((len(ps), len(qs)))
    ph = np.exp(2j * np.outer(ps, ys))
    for j, qq in enumerate(qs):
        f = np.conj(re(qq + ys)) * re(qq - ys)
        W[:, j] = (ph @ f).real * dq / math.pi
    return W
