"""Phase 7: one optical mode — squeezing (Eqs. 22–26), GKP states (Eqs. 38–39, Fig. 5), cat states (Eqs. 40, 75) and the
homodyne-heralded comb (Eqs. 66–73) of the review (Romero and Milburn, arXiv:2404.03367), on qsim/qsim/cv.py
(the same file as TechUPGRADE's sim-manager/quantum/cv.py). Every closed form is checked against an independent
route (the Fock basis against the q-grid), and every check has a power check that must fail.
    ~/quanterall/.venv-tools/bin/python tools/optics/cv_checks.py site/data [path/to/2404.03367.pdf]"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim import cv  # noqa: E402

OUT = Path(sys.argv[1])
PDF = Path(sys.argv[2]) if len(sys.argv) > 2 else None
rng = np.random.default_rng(20261030)
res = {}
A = cv.SQRT_PI


def check(name, ok, **kw):
    res[name] = {"pass": bool(ok), **kw}
    print(("PASS " if ok else "FAIL ") + name, {k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in kw.items()})
    assert ok, name


q = cv.grid(24.0, 1 << 13)
dq = q[1] - q[0]

# ── squeezing: the Fock route (expm of the generator, Hermite functions) against the closed-form wavefunction ──
CUT = 320                     # tanh(1.2)^(n/2) < 1e-12 well before the cut
worst, worst_flip = 0.0, 0.0
for _ in range(10):
    r, ph = rng.uniform(0.0, 1.2), rng.uniform(-math.pi, math.pi)
    vac = np.zeros(CUT, complex); vac[0] = 1
    f = cv.fock_to_q(cv.squeeze_op(r, ph, CUT) @ vac, q)
    for tgt, key in ((cv.squeezed_vacuum_q(q, r, ph), 0), (cv.squeezed_vacuum_q(q, r, -ph), 1)):
        ov_ = cv.inner(tgt, f, q)
        d = float(np.abs(f - tgt * ov_ / abs(ov_)).max())           # up to the global phase
        worst, worst_flip = (max(worst, d), worst_flip) if key == 0 else (worst, max(worst_flip, d))
check("squeezed_vacuum_fock_vs_closed_form", worst < 1e-9 and worst_flip > 1e-2, cases=10, max_diff=worst, power_phi_flipped=worst_flip)

# Eq. 24: S†aS = a cosh r − a† e^{−2iϕ} sinh r (the low-photon block of a large cut-off)
CUT2, BLK = 400, 10            # squeezed number states spread to ~n cosh 2r: keep the block well inside the cut
a = cv.annihilation(CUT2)
worst, worst_flip = 0.0, 0.0
for _ in range(5):
    r, ph = rng.uniform(0.1, 1.2), rng.uniform(-math.pi, math.pi)
    S = cv.squeeze_op(r, ph, CUT2)
    lhs = (S.conj().T @ a @ S)[:BLK, :BLK]
    rhs = (a * math.cosh(r) - a.conj().T * np.exp(-2j * ph) * math.sinh(r))[:BLK, :BLK]
    bad = (a * math.cosh(r) - a.conj().T * np.exp(2j * ph) * math.sinh(r))[:BLK, :BLK]
    worst, worst_flip = max(worst, float(np.abs(lhs - rhs).max())), max(worst_flip, float(np.abs(lhs - bad).max()))
check("eq24_canonical_transformation", worst < 1e-9 and worst_flip > 1e-2, cases=5, max_diff=worst, power_conjugate_phase=worst_flip)

# Eq. 26: in S(ζ)|0⟩ the rotated quadratures have variances e^{−2r} and e^{2r} — but only when the rotation of Eq. 25
# and the phase of ζ in Eq. 22 match: with ζ = r e^{−2iϕ} (Eq. 22 as printed) it takes Y1 + iY2 = (X1 + iX2)e^{+iϕ};
# Eq. 25's e^{−iϕ} needs ζ = r e^{+2iϕ} (as in Walls and Milburn). As printed together they hold only at ϕ = 0.
worst, printed_dev = 0.0, 0.0
for _ in range(5):
    r, ph = rng.uniform(0.1, 1.2), rng.uniform(0.3, 1.2) * rng.choice([-1, 1])
    vac = np.zeros(CUT2, complex); vac[0] = 1
    v = cv.squeeze_op(r, ph, CUT2) @ vac                       # ζ = r e^{−2iϕ}, Eq. 22
    var = lambda O: float(np.vdot(v, O @ O @ v).real - np.vdot(v, O @ v).real ** 2)
    for rot, consistent in ((1, True), (-1, False)):
        Y1 = a * np.exp(rot * 1j * ph) + a.conj().T * np.exp(-rot * 1j * ph)
        Y2 = -1j * (a * np.exp(rot * 1j * ph) - a.conj().T * np.exp(-rot * 1j * ph))
        dev = max(abs(var(Y1) - math.exp(-2 * r)), abs(var(Y2) - math.exp(2 * r)))
        worst, printed_dev = (max(worst, dev), printed_dev) if consistent else (worst, max(printed_dev, dev))
check("eq26_quadrature_variances", worst < 1e-9 and printed_dev > 0.1, cases=5, max_diff_consistent_signs=worst,
      max_diff_eq22_with_eq25_as_printed=printed_dev)

# ── cat states (Eqs. 40, 75) ──
b = 2.0
c_even, c_odd = cv.cat_fock(b, 0, CUT), cv.cat_fock(b, 1, CUT)
par = np.diag((-1.0) ** np.arange(CUT))
d_q = max(float(np.abs(cv.fock_to_q(c, q) - cv.cat_q(q, b, k)).max()) for k, c in ((0, c_even), (1, c_odd)))
check("cat_fock_vs_q_grid", d_q < 1e-10, beta=b, max_diff=d_q)
check("cat_normalised_and_parity", abs(np.vdot(c_even, c_even) - 1) < 1e-12 and abs(np.vdot(c_odd, c_odd) - 1) < 1e-12
      and abs(np.vdot(c_even, par @ c_even) - 1) < 1e-12 and abs(np.vdot(c_odd, par @ c_odd) + 1) < 1e-12
      and float(np.abs(c_even[1::2]).max()) < 1e-15 and float(np.abs(c_odd[0::2]).max()) < 1e-15
      and abs(np.vdot(c_even, c_odd)) < 1e-15,
      even_parity=float(np.vdot(c_even, par @ c_even).real), odd_parity=float(np.vdot(c_odd, par @ c_odd).real))
Np = 2 + 2 * math.exp(-2 * b * b)
eq40 = (cv.coherent_fock(b, CUT) + cv.coherent_fock(-b, CUT)) / Np
check("eq40_normalisation_as_printed", abs(float(np.vdot(eq40, eq40).real) - 1 / Np) < 1e-12,
      norm_squared_with_N_inverse=float(np.vdot(eq40, eq40).real), expected_one_over_N=1 / Np)
ov = abs(np.vdot(cv.coherent_fock(b, CUT), cv.coherent_fock(-b, CUT))) ** 2
check("coherent_overlap_e_minus_4_alpha_sq", abs(ov - math.exp(-4 * b * b)) < 1e-18, alpha=b, overlap=float(ov))
# Wigner at the origin is ⟨Π⟩/π: +1/π for the even cat, −1/π for the odd — through the q-grid Wigner integral
W0 = [float(cv.wigner(cv.cat_q(q, b, k), q, np.array([0.0]), np.array([0.0]))[0, 0]) for k in (0, 1)]
check("cat_wigner_origin_is_parity_over_pi", abs(W0[0] - 1 / math.pi) < 1e-9 and abs(W0[1] + 1 / math.pi) < 1e-9, even=W0[0], odd=W0[1])

# ── GKP (Eq. 38): 10 dB of squeezing in the peaks, Δ² = 0.01 in the envelope ──
qg = cv.grid(48.0, 1 << 14)
DB = 10.0
pv = 0.5 * 10 ** (-DB / 10)                      # density variance of each peak: Var(q) = e^{−2r}/2, 10 dB below vacuum
D2 = 0.01
g0, g1 = cv.gkp_q(qg, 0, D2, pv), cv.gkp_q(qg, 1, D2, pv)
ov01 = abs(cv.inner(g0, g1, qg))
zbar = lambda g: cv.inner(g, np.exp(1j * math.pi * qg / A) * g, qg).real


def gauss_sum(mu_a, mu_b, kq):
    """⟨g_a|e^{ik q}|g_b⟩ for two Eq. 38 states as an analytic double sum over their Gaussian peaks (no grid):
    ψ_c1 ψ_c2 = (2πv)^{−1/2} e^{−(c1−c2)²/8v} e^{−(q−m)²/2v}, m = (c1 + c2)/2, so ∫ e^{ikq} … = e^{ikm − k²v/2}"""
    S_ = range(-40, 41)
    ca = np.array([(2 * s + mu_a) * A for s in S_]); cb = np.array([(2 * s + mu_b) * A for s in S_])
    wa, wb = np.exp(-0.5 * D2 * ca ** 2), np.exp(-0.5 * D2 * cb ** 2)
    C1, C2 = np.meshgrid(ca, cb, indexing="ij")
    M = np.outer(wa, wb) * np.exp(-(C1 - C2) ** 2 / (8 * pv)) * np.exp(1j * kq * (C1 + C2) / 2 - kq * kq * pv / 2)
    return complex(M.sum())


n0, n1 = gauss_sum(0, 0, 0).real, gauss_sum(1, 1, 0).real
z_pred0, z_pred1 = gauss_sum(0, 0, math.pi / A).real / n0, gauss_sum(1, 1, math.pi / A).real / n1
ov_pred = abs(gauss_sum(0, 1, 0)) / math.sqrt(n0 * n1)
check("gkp_logical_z", abs(zbar(g0) - z_pred0) < 1e-10 and abs(zbar(g1) - z_pred1) < 1e-10 and abs(ov01 - ov_pred) < 1e-10,
      Z_on_0=zbar(g0), Z_on_1=zbar(g1), predicted=[z_pred0, z_pred1], single_peak_value=math.exp(-math.pi ** 2 * pv / (2 * A * A)),
      overlap_0_1=ov01, overlap_predicted=ov_pred)

# X̄ = e^{−ipα} translates q by α: X̄|0̃⟩ has |1̃⟩'s peaks but |0̃⟩'s envelope (centred α off), so the fidelity is below 1;
# the grid (a spectral shift) against the analytic double sum
xg0 = cv.translate(g0, qg, A)
fx = abs(cv.inner(g1, xg0, qg)) ** 2
cs_ = np.arange(-40, 41)
c1, w1 = (2 * cs_ + 1) * A, np.exp(-0.5 * D2 * ((2 * cs_ + 1) * A) ** 2)
cx, wx = 2 * cs_ * A + A, np.exp(-0.5 * D2 * (2 * cs_ * A) ** 2)
ovx = lambda ca, wa, cb, wb: float((np.outer(wa, wb) * np.exp(-np.subtract.outer(ca, cb) ** 2 / (8 * pv))).sum())
fx_pred = ovx(c1, w1, cx, wx) ** 2 / (ovx(c1, w1, c1, w1) * ovx(cx, wx, cx, wx))
check("gkp_logical_x_shift", abs(fx - fx_pred) < 1e-10, fidelity_X0_with_1=fx, predicted=fx_pred)
# the p-side: |0̃⟩ is a comb of spacing π/α in p, |±̃⟩ of spacing 2π/α offset by 0 / π/α
p, g0p = cv.to_p(g0, qg)
plus = cv.normalise(g0 + g1, qg); minus = cv.normalise(g0 - g1, qg)
_, plp = cv.to_p(plus, qg); _, mip = cv.to_p(minus, qg)
dp = p[1] - p[0]
peaks = lambda ph: sorted(p[i] for i in range(1, len(p) - 1) if abs(ph[i]) > abs(ph[i - 1]) and abs(ph[i]) >= abs(ph[i + 1]) and abs(ph[i]) > 0.3 * np.abs(ph).max())
pk0, pkp, pkm = peaks(g0p), peaks(plp), peaks(mip)
sp0 = float(np.median(np.diff(pk0))); spp = float(np.median(np.diff(pkp)))
off = min(abs(x) for x in pkm)
check("gkp_p_combs", abs(sp0 - math.pi / A) < 2 * dp and abs(spp - 2 * math.pi / A) < 2 * dp and abs(off - math.pi / A) < 2 * dp,
      spacing_0=sp0, spacing_plus=spp, minus_offset=off, expected=[math.pi / A, 2 * math.pi / A, math.pi / A], grid_dp=float(dp))
pn = float(np.sum(np.abs(g0p) ** 2) * dp)
check("fft_unitary", abs(pn - 1) < 1e-9, norm_in_p=pn)

# ── readout under Gaussian shift noise: exact integral, Monte Carlo, closed form (σ_eff² = σ² + peak variance) ──
rows = []
SHOTS = 200_000
d0 = np.abs(g0) ** 2
dplus_p = np.abs(plp) ** 2
var_p_peak = D2 / 2                                  # the q-envelope's Fourier width: each p-peak's density variance
for sigma in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8):
    ex_q = cv.readout_error_exact(d0, qg, sigma, A, False)
    mc_q = cv.readout_mc(d0, qg, sigma, A, False, SHOTS, rng) / SHOTS
    cf_q = cv.readout_error_closed(math.sqrt(sigma ** 2 + pv), A)
    ex_p = cv.readout_error_exact(dplus_p, p, sigma, math.pi / A, False)
    mc_p = cv.readout_mc(dplus_p, p, sigma, math.pi / A, False, SHOTS, rng) / SHOTS
    cf_p = cv.readout_error_closed(math.sqrt(sigma ** 2 + var_p_peak), math.pi / A)
    rows.append({"sigma": sigma, "q_exact": ex_q, "q_mc": mc_q, "q_closed": cf_q, "p_exact": ex_p, "p_mc": mc_p, "p_closed": cf_p})
    print(f"  σ {sigma:.1f}: Z readout of |0̃⟩ exact {ex_q:.5f} MC {mc_q:.5f} closed {cf_q:.5f};  X readout of |+̃⟩ exact {ex_p:.5f} MC {mc_p:.5f} closed {cf_p:.5f}")
mc_ok = all(abs(r[k + "_mc"] - r[k + "_exact"]) < 5 * math.sqrt(max(r[k + "_exact"], 1e-6) / SHOTS) + 2e-5 for r in rows for k in ("q", "p"))
cf_dev = max(abs(r[k + "_closed"] - r[k + "_exact"]) for r in rows for k in ("q", "p"))
check("gkp_readout_mc_vs_exact", mc_ok, shots=SHOTS, points=len(rows) * 2)
check("gkp_readout_closed_form_vs_exact", cf_dev < 2e-3, max_abs_diff=cf_dev,
      note="closed form neglects the envelope; differences are the envelope's")
# power: decoding on a spacing 5 % off, and the closed form without the peak width, must both disagree
bad_sp = cv.readout_error_exact(d0, qg, 0.3, 1.05 * A, False)
bad_cf = cv.readout_error_closed(0.3, A)
check("gkp_readout_power", abs(bad_sp - rows[3]["q_exact"]) > 0.01 and abs(bad_cf - rows[3]["q_exact"]) > 1e-3,
      exact_at_0_3=rows[3]["q_exact"], wrong_spacing=bad_sp, closed_without_peak_width=bad_cf)
with open(OUT / "gkp-readout.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# Eq. 39: a fixed shift u is corrected while |u| < α/2 (q) and |u| < π/(2α) (p); beyond, the logical value flips
# a hard decision on narrow p-peaks needs a fine p-grid: |+̃⟩ on a q-grid four times wider (dp four times finer)
q4 = cv.grid(192.0, 1 << 16)
plus4 = cv.normalise(cv.gkp_q(q4, 0, D2, pv) + cv.gkp_q(q4, 1, D2, pv), q4)
p4 = cv.to_p(plus4, q4)[0]
shift_rows = []
for u in np.linspace(0, 1.6, 81):
    e_q = cv.readout_error_exact(np.abs(cv.translate(g0, qg, u * A)) ** 2, qg, 0.0, A, False)          # exact shifts:
    e_p = cv.readout_error_exact(np.abs(cv.to_p(plus4 * np.exp(1j * u * math.pi / A * q4), q4)[1]) ** 2,   # in q spectrally,
                                 p4, 0.0, math.pi / A, False)                                              # in p by e^{iuq}
    shift_rows.append({"shift_over_bound": float(u) / 0.5, "q_error": e_q, "p_error": e_p})
with open(OUT / "gkp-shift.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(shift_rows[0])); w.writeheader(); w.writerows(shift_rows)
# with finite peaks the step at the bound is smeared by the peak width: against the same decision for one Gaussian
# peak (variance pv in q, D2/2 in p; the rest is the envelope and the neighbouring peaks), and one half at the bound
dev = max(max(abs(r["q_error"] - cv.wrong_parity_prob(r["shift_over_bound"] * A / 2, math.sqrt(pv), A, False)),
              abs(r["p_error"] - cv.wrong_parity_prob(r["shift_over_bound"] * math.pi / A / 2, math.sqrt(var_p_peak), math.pi / A, False)))
          for r in shift_rows)
at_bound = next(r for r in shift_rows if abs(r["shift_over_bound"] - 1) < 1e-9)
check("eq39_correctable_shifts", dev < 5e-3 and abs(at_bound["q_error"] - 0.5) < 5e-3 and abs(at_bound["p_error"] - 0.5) < 5e-3,
      max_dev_from_single_peak=dev, q_error_at_bound=at_bound["q_error"], p_error_at_bound=at_bound["p_error"])


# ── Fig. 5: the vector curve read from the paper's PDF, against Eq. 38 ──
fig5 = {}
if PDF and PDF.exists():
    import pymupdf
    page = pymupdf.open(PDF)[10]
    words = [w for w in page.get_text("words") if w[1] < 330]
    curves = [d for d in page.get_drawings() if sum(1 for it in d["items"] if it[0] in ("l", "c")) > 50]
    for name, side in (("zero", lambda x0: x0 < 300), ("one", lambda x0: x0 >= 300)):
        ticks = {w[4]: (w[0] + w[2]) / 2 for w in words if side(w[0]) and w[4] in ("-40", "-20", "20", "40") and w[1] > 190}
        ylab = sorted((w[1], float(w[4])) for w in words if side(w[0]) and w[4] in ("0.2", "0.4", "0.6", "0.8", "1.0"))
        xs = (ticks["40"] - ticks["-40"]) / 80; x0 = (ticks["40"] + ticks["-40"]) / 2
        ys = (ylab[-1][0] - ylab[0][0]) / (ylab[0][1] - ylab[-1][1])  # points per unit (labels run top-down)
        pts = []
        for d in curves:
            if not side(d["rect"].x0):
                continue
            for it in d["items"]:
                if it[0] == "l":
                    pts += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
                elif it[0] == "c":
                    pts += [(it[1].x, it[1].y), (it[4].x, it[4].y)]
        arr = np.array(sorted(set(pts)))
        base = arr[:, 1].max()
        x, y = (arr[:, 0] - x0) / xs, (base - arr[:, 1]) / ys
        pk = [(float(x[i]), float(y[i])) for i in range(1, len(y) - 1) if y[i] >= y[i - 1] and y[i] >= y[i + 1] and y[i] > 0.02]
        px, ph = np.array([t[0] for t in pk]), np.array([t[1] for t in pk])
        cfit = np.linalg.lstsq(np.vstack([np.ones_like(px), px ** 2]).T, np.log(ph), rcond=None)[0]
        fig5[name] = {"peaks": [[round(a_, 3), round(b_, 3)] for a_, b_ in pk], "spacing": float(np.median(np.diff(px))),
                      "envelope_c": float(-cfit[1]), "curve": [[round(float(a_), 3), round(float(b_), 4)] for a_, b_ in zip(x, y)]}
    # unnormalised Gaussian peaks, fitted variance of each peak's density (grid search: no scipy needed)
    def fit_var(curve, c, odd):
        cx, cy = np.array(curve)[:, 0], np.array(curve)[:, 1]
        cs = (np.arange(-14, 15) * 2 + (1 if odd else 0)) * A
        best = None
        for V in np.linspace(0.02, 0.3, 281):
            m = sum(np.exp(-c * x0 ** 2) * np.exp(-(cx - x0) ** 2 / (2 * V)) / math.sqrt(2 * math.pi * V) for x0 in cs)
            e = float(np.sqrt(np.mean((m - cy) ** 2)))
            best = min(best or (e, V), (e, V))
        return best
    for name, odd in (("zero", False), ("one", True)):
        e, V = fit_var(fig5[name]["curve"], fig5[name]["envelope_c"], odd)
        fig5[name].update({"peak_variance_fit": V, "fit_rms": e})
    ratio = fig5["one"]["envelope_c"] / fig5["zero"]["envelope_c"]
    check("fig5_digitised", abs(fig5["zero"]["spacing"] - 2 * A) < 0.02 and abs(fig5["one"]["spacing"] - 2 * A) < 0.02,
          spacing_zero=fig5["zero"]["spacing"], spacing_one=fig5["one"]["spacing"], two_sqrt_pi=2 * A,
          envelope_c_zero=fig5["zero"]["envelope_c"], envelope_c_one=fig5["one"]["envelope_c"], ratio=ratio,
          peak_var_zero=fig5["zero"]["peak_variance_fit"], peak_var_one=fig5["one"]["peak_variance_fit"],
          caption_e_plus_2r=math.exp(2.4), e_minus_2r=math.exp(-2.4))
    res["fig5_reading"] = (
        f"peaks at multiples of √π on the plotted axis (α = √π on q); each peak's variance fits {fig5['zero']['peak_variance_fit']:.2f} "
        "(the plotted polyline is coarse at the tips) — near e^{−2r} = 0.09 for r = 1.2, nowhere near the caption's e^{2r} = 11; the |1̃⟩ envelope falls off four times faster than the |0̃⟩ one "
        f"(ratio {ratio:.2f}), which one Δ in Eq. 38 cannot give. With the caption's Δ̃ = 0.01 both are matched by "
        "amplitudes e^{−Δ̃(sα)²} for |0̃⟩ and e^{−Δ̃[(2s+1)α]²} for |1̃⟩: the |0̃⟩ panel uses s where Eq. 38 has 2s.")
    json.dump(fig5, open(OUT / "gkp-fig5-digitised.json", "w"))

# the states drawn per Eq. 38 for the site's chart (Fig. 5's peak width and its |1̃⟩ envelope, both panels the same Δ)
xs_plot = np.linspace(-40, 40, 3201)
V5 = math.exp(-2.4)
D5 = 2 * 0.00499          # the |0̃⟩ panel's envelope read as Eq. 38's Δ² (amplitude e^{−½Δ²x²} ⇒ intensity e^{−Δ²x²})
unnorm = lambda mu, D: sum(math.exp(-0.5 * D * ((2 * s + mu) * A) ** 2) * np.exp(-(xs_plot - (2 * s + mu) * A) ** 2 / (4 * V5)) for s in range(-15, 16))
with open(OUT / "gkp-states.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["q", "zero_eq38", "one_eq38"])
    z, o = np.abs(unnorm(0, D5 / 2)) ** 2 / math.sqrt(2 * math.pi * V5), np.abs(unnorm(1, D5 / 2)) ** 2 / math.sqrt(2 * math.pi * V5)
    for i in range(0, len(xs_plot)):
        w.writerow([round(float(xs_plot[i]), 3), round(float(z[i]), 5), round(float(o[i]), 5)])

# ── the homodyne-heralded comb (Eqs. 66–73) ──
qh = cv.grid(30.0, 1 << 13)
BETA, THETA, PHI = 3.0, A, 0.0                  # probe amplitude (Eq. 69's r), coupling chosen so the y = 0 comb spacing is √π
psi_in = cv.squeezed_vacuum_q(qh, 1.5, math.pi / 2)   # phase-squeezed: q anti-squeezed, Var(q) = e^{3}/2
ys = np.linspace(-12, 12, 4801)
dy = ys[1] - ys[0]
compl = [float(np.sum([abs(cv.upsilon(np.array([qq]), y, BETA, THETA, PHI)[0]) ** 2 for y in ys]) * dy) for qq in (-1.3, 0.0, 0.7, 2.2)]
printed = [float(np.sum([cv.upsilon_sq_printed(np.array([qq]), y, BETA, THETA, PHI)[0] for y in ys]) * dy) for qq in (-1.3, 0.0, 0.7, 2.2)]
check("eq69_povm_complete_eq70_not", max(abs(c - 1) for c in compl) < 1e-9 and max(abs(c - 1 / math.sqrt(2)) for c in printed) < 1e-9,
      eq69_integral=compl, eq70_printed_integral=printed)
py, out = cv.heralded_comb(qh, psi_in, BETA, THETA, PHI, ys)
pyn = float(np.sum(py) * (ys[1] - ys[0]))
check("eq73_probability_normalised", abs(pyn - 1) < 1e-8, integral=pyn)
cdf = np.cumsum(py); cdf /= cdf[-1]
y_samples = [float(np.interp(u, cdf, ys)) for u in rng.random(5)]
comb_rows = []
for y in [0.0] + y_samples:
    o = out(y); dens = np.abs(o) ** 2
    pk = [float(qh[i]) for i in range(1, len(qh) - 1) if dens[i] > dens[i - 1] and dens[i] >= dens[i + 1] and dens[i] > 0.05 * dens.max()]
    gaps = sorted({round(float(g), 2) for g in np.diff(pk)})
    comb_rows.append({"y": y, "peaks": len(pk), "gaps": gaps[:4]})
    print(f"  heralded y = {y:+.3f}: {len(pk)} peaks, gaps {gaps[:4]}")
# for any y the peaks sit at θq − ϕ = ±arccos(y/2β) + 2πk: two interleaved gaps, 2 arccos(y/2β)/θ and its complement
gap_dev = 0.0
for r_ in comb_rows:
    if abs(r_["y"]) < 2 * BETA:
        ac = math.acos(r_["y"] / (2 * BETA))
        want = sorted({2 * ac / THETA, (2 * math.pi - 2 * ac) / THETA})
        gap_dev = max(gap_dev, max(min(abs(g - w_) for w_ in want) for g in r_["gaps"]))
check("heralded_comb_gaps", gap_dev < 0.03,   # maxima pulled slightly inwards by the input envelope
      max_dev=gap_dev, samples=len(comb_rows))
y0 = comb_rows[0]
check("heralded_comb_y0_uniform", all(abs(g - math.pi / THETA) < 0.02 for g in y0["gaps"]), gaps=y0["gaps"], expected=math.pi / THETA)
res["heralded_comb"] = {"beta": BETA, "theta": THETA, "phi": PHI, "input": "S(1.5) phase-squeezed (q anti-squeezed)", "samples": comb_rows}
with open(OUT / "heralded-comb.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["q", "input"] + [f"y={r['y']:+.2f}" for r in comb_rows[:3]])
    dens_in = np.abs(psi_in) ** 2
    outs = [np.abs(out(r["y"])) ** 2 for r in comb_rows[:3]]
    for i in range(0, len(qh), 4):
        if abs(qh[i]) <= 10:
            w.writerow([round(float(qh[i]), 4), round(float(dens_in[i]), 6)] + [round(float(o[i]), 6) for o in outs])
with open(OUT / "heralded-comb-py.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["y", "P"]); [w.writerow([round(float(ys[i]), 3), round(float(py[i]), 7)]) for i in range(0, len(ys), 8)]

# Wigner functions for the chart: the odd cat (β = 2) and |0̃⟩ (Δ² = 0.05, 6 dB)
wq, wp = np.linspace(-5, 5, 161), np.linspace(-5, 5, 161)
Wc = cv.wigner(cv.cat_q(q, 2.0, 1), q, wq, wp)
Wg = cv.wigner(cv.gkp_q(q, 0, 0.05, 0.5 * 10 ** -0.6), q, wq, wp)
json.dump({"q": wq.round(4).tolist(), "p": wp.round(4).tolist(), "cat_odd_beta2": Wc.round(5).tolist(), "gkp0_6dB": Wg.round(5).tolist(),
           "cat_min": float(Wc.min()), "gkp_min": float(Wg.min())}, open(OUT / "cv-wigner.json", "w"))
print(f"  Wigner minima: odd cat {Wc.min():.4f} (−1/π = {-1 / math.pi:.4f}), GKP |0̃⟩ {Wg.min():.4f}")

json.dump(res, open(OUT / "cv-checks.json", "w"), indent=1)
print("all", len(res), "entries written")
