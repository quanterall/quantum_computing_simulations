"""Phase 7 cross-check: QxCV (Elixir, the BEAM service) against qsim/qsim/cv.py (Python, TechUPGRADE's quantum/cv.py),
written separately — GKP Z readout on the same q-grid, X readout by an analytic Fourier sum (BEAM) against an FFT
(Python), cat photon statistics under loss by the generating function's distribution (BEAM) against the Fock basis
thinned binomially (Python), the heralded comb on the same grid.
    ~/quanterall/.venv-tools/bin/python tools/optics/xcheck_cv.py site/data [qx_cv.ex copy to compile in the container]"""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "qsim"))
from qsim import cv  # noqa: E402

OUT = Path(sys.argv[1])
A = cv.SQRT_PI
subprocess.run(["docker", "exec", "quanterall_quantum", "mkdir", "-p", "/tmp/cv"], check=True)
subprocess.run(["docker", "cp", str(HERE / "cv_beam.exs"), "quanterall_quantum:/tmp/cv/cv_beam.exs"], check=True)
args = []
if len(sys.argv) > 2:
    subprocess.run(["docker", "cp", sys.argv[2], "quanterall_quantum:/tmp/cv/qx_cv.ex"], check=True)
    args = ["/tmp/cv/qx_cv.ex"]
out = subprocess.run(["docker", "exec", "-w", "/app", "-e", "MIX_ENV=prod", "quanterall_quantum", "mix", "run", "--no-start",
                      "/tmp/cv/cv_beam.exs", *args], capture_output=True, text=True).stdout
B = json.loads([l for l in out.splitlines() if l.startswith("{")][-1])
res = {}


def check(name, ok, **kw):
    res[name] = {"pass": bool(ok), **kw}
    print(("PASS " if ok else "FAIL ") + name, {k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in kw.items()})
    assert ok, name


q = cv.grid(48.0, 8192)
d0 = np.abs(cv.gkp_q(q, 0, 0.01, 0.05)) ** 2
dd = float(np.abs(d0[::97] - np.array(B["d0_sample"])).max())
check("gkp_density_same_grid", dd < 1e-12, max_diff=dd)
z_py = [cv.readout_error_exact(d0, q, s, A, False) for s in (0.0, 0.2, 0.4, 0.6)]
dz = max(abs(a - b) for a, b in zip(z_py, B["z"]))
check("gkp_z_readout", dz < 1e-12, python=z_py, beam=B["z"], max_diff=dz)

# X: BEAM sums the Fourier transform over the peaks on its own p-grid; Python takes the FFT of the q-wavefunction
q4 = cv.grid(192.0, 1 << 16)
plus = cv.normalise(cv.gkp_q(q4, 0, 0.01, 0.05) + cv.gkp_q(q4, 1, 0.01, 0.05), q4)
p, ph = cv.to_p(plus, q4)
x_py = [cv.readout_error_exact(np.abs(ph) ** 2, p, s, math.pi / A, False) for s in (0.0, 0.2, 0.4, 0.6)]
dx = max(abs(a - b) for a, b in zip(x_py, B["x"]))
check("gkp_x_readout_fourier_sum_vs_fft", dx < 1e-6, python=x_py, beam=B["x"], max_diff=dx)
dc = max(abs(cv.readout_error_closed(s, A) - b) for s, b in zip((0.2, 0.4), B["closed"]))
check("closed_form", dc < 1e-14, max_diff=dc)
mc = B["mc"]
z = (mc["wrong"] / mc["shots"] - mc["exact"]) / math.sqrt(mc["exact"] * (1 - mc["exact"]) / mc["shots"])
check("beam_sampled_readout_vs_exact", abs(z) < 4, wrong=mc["wrong"], shots=mc["shots"], exact=mc["exact"], z_score=z, time_us=mc["us"])

# cats under loss: Fock amplitudes, then binomial thinning
CUT = 120
worst = 0.0
for c, (b, par, eta) in zip(B["cats"], [(2.0, 0, 1.0), (2.0, 1, 0.9), (1.3, 0, 0.7)]):
    pn = np.abs(cv.cat_fock(b, par, CUT)) ** 2
    outp = np.zeros(CUT)
    for n in range(CUT):
        for m in range(n + 1):
            outp[m] += pn[n] * math.comb(n, m) * eta ** m * (1 - eta) ** (n - m)
    worst = max(worst, float(np.abs(outp[:20] - np.array(c["dist"])).max()), abs(c["parity"] - c["closed"]))
    pi_py = float(np.sum(outp * (-1.0) ** np.arange(CUT)))
    worst = max(worst, abs(pi_py - c["closed"]))
check("cat_loss_fock_vs_generating_function", worst < 1e-12, cases=3, max_diff=worst)

# the comb at y = 0, the same grid
qc = cv.grid(16.0, 3072)
din = np.abs(cv.squeezed_vacuum_q(qc, 1.5, math.pi / 2)) ** 2
dens = np.abs(cv.upsilon(qc, 0.0, 3.0, A, 0.0)) ** 2 * din
dens /= np.sum(dens) * (qc[1] - qc[0])
bd = np.array(B["comb_y0"]["dens"])
dcomb = float(np.abs(dens[::4] - bd[:, 1]).max())
check("heralded_comb_y0", dcomb < 1e-5 and abs(B["comb_y0"]["p_y"] - 1) < 1e-5, max_diff_rounded_output=dcomb,
      p_y_integral=B["comb_y0"]["p_y"], gaps=B["comb_y0"]["gaps"])

json.dump(res, open(OUT / "cv-xcheck.json", "w"), indent=1)
print("cross-check written")
