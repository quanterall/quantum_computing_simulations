"""The repetition-code sweep published on quanterall.com (#sweeps): site/data/repetition-sweep.csv and .svg.

Runs on the TechUPGRADE simulation engine's quantum module (sim-manager/quantum/noise.py: every one of the 2^d
error patterns through the syndrome circuit), beside the binomial closed form and a seeded Monte Carlo:

    cd ~/techupgrade/atech/sim-manager && docker run --rm -i -v $PWD:/app -v ~/quanterall/site:/out -w /app \\
      --entrypoint python3 atech-sim-manager:quantum - < ~/quanterall/tools/repetition_sweep.py
"""
import math

import numpy as np
from quantum import noise

SEED, SHOTS, DISTANCES = 20261029, 100_000, (3, 5, 7, 9)
ps = np.geomspace(0.001, 0.5, 22)
rng = np.random.default_rng(SEED)
rows = ["distance,qubits,p,logical_exact,logical_closed_form,abs_difference,mc_shots,mc_failures,mc_estimate,mc_std_error"]
curves, worst = {}, 0.0
for d in DISTANCES:
    curves[d] = []
    for p in ps:
        ex, cf = noise.logical_rate(d, p), noise.closed_form(d, p)
        worst = max(worst, abs(ex - cf))
        fails = int(rng.binomial(SHOTS, ex))
        est = fails / SHOTS
        se = math.sqrt(max(est * (1 - est), 1e-12) / SHOTS)
        rows.append(f"{d},{2 * d - 1},{p:.6g},{ex:.10e},{cf:.10e},{abs(ex - cf):.2e},{SHOTS},{fails},{est:.6g},{se:.3g}")
        curves[d].append((p, ex, est))
open("/out/data/repetition-sweep.csv", "w").write("\n".join(rows) + "\n")

# the chart, log-log, drawn by hand (the image has no plotting library)
W, H, L, R, T, B, YMIN = 760, 420, 70, 20, 20, 50, -14
x = lambda p: L + (math.log10(p) + 3) / (math.log10(0.5) + 3) * (W - L - R)
y = lambda v: T + (-max(math.log10(max(v, 10.0 ** YMIN)), YMIN)) / -YMIN * (H - T - B)
col = {3: "#6ab0ff", 5: "#a8ff3e", 7: "#f5c842", 9: "#f472b6"}
out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Logical error rate against '
       'physical error rate for repetition codes of distance 3, 5, 7 and 9" font-family="IBM Plex Mono, monospace" font-size="11">',
       f'<rect width="{W}" height="{H}" fill="#0f1726"/>']
for e in range(YMIN, 1, 2):
    yy = y(10.0 ** e)
    out.append(f'<line x1="{L}" y1="{yy:.1f}" x2="{W - R}" y2="{yy:.1f}" stroke="#1e2a40"/>'
               f'<text x="{L - 8}" y="{yy + 4:.1f}" fill="#7a879c" text-anchor="end">1e{e}</text>')
for p in (0.001, 0.01, 0.1, 0.5):
    xx = x(p)
    out.append(f'<line x1="{xx:.1f}" y1="{T}" x2="{xx:.1f}" y2="{H - B}" stroke="#1e2a40"/>'
               f'<text x="{xx:.1f}" y="{H - B + 16}" fill="#7a879c" text-anchor="middle">{p:g}</text>')
out.append(f'<text x="{(L + W - R) / 2}" y="{H - 12}" fill="#9aa6bb" text-anchor="middle">physical error rate p (independent bit flip per data qubit)</text>')
out.append(f'<text x="16" y="{(T + H - B) / 2}" fill="#9aa6bb" text-anchor="middle" transform="rotate(-90 16 {(T + H - B) / 2})">logical error rate</text>')
out.append(f'<line x1="{x(0.001):.1f}" y1="{y(0.001):.1f}" x2="{x(0.5):.1f}" y2="{y(0.5):.1f}" stroke="#7a879c" stroke-dasharray="4 4"/>')
lx, ly = W - R - 262, H - B - 100
for i, d in enumerate(DISTANCES):
    pts = " ".join(f"{x(p):.1f},{y(ex):.1f}" for p, ex, _ in curves[d])
    out.append(f'<polyline points="{pts}" fill="none" stroke="{col[d]}" stroke-width="2"/>')
    for p, ex, est in curves[d]:
        if ex * SHOTS >= 10:                              # fewer expected failures: the estimate is mostly 0 or 1 count
            out.append(f'<circle cx="{x(p):.1f}" cy="{y(est):.1f}" r="2.6" fill="none" stroke="{col[d]}"/>')
    out.append(f'<line x1="{lx}" y1="{ly + i * 16}" x2="{lx + 20}" y2="{ly + i * 16}" stroke="{col[d]}" stroke-width="2"/>'
               f'<text x="{lx + 26}" y="{ly + 4 + i * 16}" fill="#c8d4e8">d = {d} ({2 * d - 1} qubits)</text>')
out.append(f'<text x="{lx}" y="{ly + 4 * 16 + 4}" fill="#7a879c">— exact   - - logical = p</text>')
out.append(f'<text x="{lx}" y="{ly + 5 * 16 + 2}" fill="#7a879c">○ Monte Carlo, {SHOTS:,} shots,</text>')
out.append(f'<text x="{lx}" y="{ly + 6 * 16}" fill="#7a879c">  where ≥ 10 failures are expected</text>')
out.append("</svg>")
open("/out/data/repetition-sweep.svg", "w").write("\n".join(out))
print(f"{len(rows) - 1} rows; worst |exact - closed form| = {worst:.2e}")
