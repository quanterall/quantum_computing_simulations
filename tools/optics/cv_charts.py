"""Phase 7 charts for the site, from the data tools/optics/cv_checks.py writes:
gkp-fig5.svg (the paper's Fig. 5 read from its PDF, against Eq. 38), gkp-readout.svg (logical error under shift noise),
heralded-comb.svg (Eqs. 66–73), cv-wigner.png (the odd cat and a GKP |0̃⟩).
    ~/quanterall/.venv-tools/bin/python tools/optics/cv_charts.py site/data"""
import csv
import json
import math
import sys
from pathlib import Path

D = Path(sys.argv[1])
BG, GRID, TICK, LAB = "#0f1726", "#1e2a40", "#7a879c", "#9aa6bb"
BLUE, GREEN, PINK, YELLOW, GREY = "#6ab0ff", "#a8ff3e", "#f472b6", "#f5c842", "#8b97ab"


class Panel:
    def __init__(self, s, x0, y0, w, h, xr, yr, xt, yt, xlabel, ylabel=None, title=None):
        self.s, self.x0, self.y0, self.w, self.h, self.xr, self.yr = s, x0, y0, w, h, xr, yr
        for v in yt:
            s.append(f'<line x1="{x0}" y1="{self.y(v):.1f}" x2="{x0 + w}" y2="{self.y(v):.1f}" stroke="{GRID}"/>'
                     f'<text x="{x0 - 6}" y="{self.y(v) + 4:.1f}" fill="{TICK}" text-anchor="end">{v:g}</text>')
        for v in xt:
            s.append(f'<line x1="{self.x(v):.1f}" y1="{y0}" x2="{self.x(v):.1f}" y2="{y0 + h}" stroke="{GRID}"/>'
                     f'<text x="{self.x(v):.1f}" y="{y0 + h + 15}" fill="{TICK}" text-anchor="middle">{v:g}</text>')
        s.append(f'<text x="{x0 + w / 2}" y="{y0 + h + 33}" fill="{LAB}" text-anchor="middle">{xlabel}</text>')
        if ylabel:
            cy = y0 + h / 2
            s.append(f'<text x="{x0 - 44}" y="{cy}" fill="{LAB}" text-anchor="middle" transform="rotate(-90 {x0 - 44} {cy})">{ylabel}</text>')
        if title:
            s.append(f'<text x="{x0 + 6}" y="{y0 + 14}" fill="{LAB}">{title}</text>')

    def x(self, v):
        return self.x0 + (v - self.xr[0]) / (self.xr[1] - self.xr[0]) * self.w

    def y(self, v):
        v = min(max(v, self.yr[0]), self.yr[1])
        return self.y0 + (self.yr[1] - v) / (self.yr[1] - self.yr[0]) * self.h

    def line(self, pts, color, width=1.6, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        pts = [(a, b) for a, b in pts if self.xr[0] <= a <= self.xr[1]]
        self.s.append(f'<polyline points="{" ".join(f"{self.x(a):.1f},{self.y(b):.1f}" for a, b in pts)}" fill="none" stroke="{color}" stroke-width="{width}"{d}/>')

    def dots(self, pts, color, r=3):
        for a, b in pts:
            self.s.append(f'<circle cx="{self.x(a):.1f}" cy="{self.y(b):.1f}" r="{r}" fill="none" stroke="{color}"/>')

    def text(self, x, y, t, color, anchor="start"):
        self.s.append(f'<text x="{x:.1f}" y="{y:.1f}" fill="{color}" text-anchor="{anchor}">{t}</text>')


def svg(W, H, aria):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="{aria}" font-family="IBM Plex Mono, monospace" font-size="11">',
            f'<rect width="{W}" height="{H}" fill="{BG}"/>']


# ── Fig. 5 against Eq. 38 ──
fig5 = json.load(open(D / "gkp-fig5-digitised.json"))
rows = list(csv.DictReader(open(D / "gkp-states.csv")))
s = svg(760, 350, "The review's Figure 5 read from its PDF against Eq. 38 drawn with one envelope for both states: the zero state matches, the one state in the figure is twice as narrow")
for k, (key, col, name) in enumerate((("zero", "zero_eq38", "|0̃⟩"), ("one", "one_eq38", "|1̃⟩"))):
    P = Panel(s, 56 + k * 360, 26, 320, 240, (-40, 40), (0, 1.25), (-40, -20, 0, 20, 40), (0, 0.4, 0.8, 1.2), "q", "probability density" if k == 0 else None, title=name)
    P.line([tuple(p) for p in fig5[key]["curve"]], BLUE, 3.2)
    P.line([(float(r["q"]), float(r[col])) for r in rows], YELLOW, 1.0)
P.text(56, 318, "— Eq. 38 with one Δ for both states (Δ fitted to the |0̃⟩ panel), peak variance e^(−2r)", YELLOW)
P.text(56, 336, "▬ the paper's Fig. 5, read from the vector curve in its PDF", BLUE)
s.append("</svg>")
(D / "gkp-fig5.svg").write_text("\n".join(s))

# ── readout error under shift noise ──
rows = list(csv.DictReader(open(D / "gkp-readout.csv")))
s = svg(760, 380, "GKP logical readout error against Gaussian shift noise: the exact integral, Monte Carlo and the closed form agree, for Z (q) and X (p)")
P = Panel(s, 64, 20, 676, 300, (0, 0.8), (0, 0.3), (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
          "shift-noise standard deviation σ (vacuum: σ = 0.71)", "logical error")
for k, c in (("q", BLUE), ("p", PINK)):
    P.line([(float(r["sigma"]), float(r[f"{k}_closed"])) for r in rows], c, 1.2, "5 4")
    P.line([(float(r["sigma"]), float(r[f"{k}_exact"])) for r in rows], c, 2)
    P.dots([(float(r["sigma"]), float(r[f"{k}_mc"])) for r in rows], c)
P.text(80, 50, "Z readout of |0̃⟩ (q, spacing √π)", BLUE)
P.text(80, 68, "X readout of |+̃⟩ (p, spacing √π)", PINK)
P.text(80, 92, "— exact   ○ Monte Carlo, 200,000 shots   - - ideal comb, σ² + peak variance", LAB)
s.append("</svg>")
(D / "gkp-readout.svg").write_text("\n".join(s))

# ── the heralded comb ──
rows = list(csv.DictReader(open(D / "heralded-comb.csv")))
keys = [k for k in rows[0] if k.startswith("y=")]
s = svg(760, 390, "Homodyne-heralded comb: a phase-squeezed input, and the output for the outcome y = 0 (evenly spaced peaks) and for another outcome (pairs)")
top = max(max(float(r[k]) for r in rows) for k in keys[:2])
P = Panel(s, 64, 20, 676, 270, (-10, 10), (0, math.ceil(top * 10) / 10), (-10, -5, 0, 5, 10), [0, 0.2, 0.4, 0.6, 0.8],
          "q", "probability density")
P.line([(float(r["q"]), float(r["input"]) * 6) for r in rows], GREY, 1.2, "4 3")
P.line([(float(r["q"]), float(r[keys[0]])) for r in rows], YELLOW, 1.6)
P.line([(float(r["q"]), float(r[keys[1]])) for r in rows], GREEN, 1.2)
P.text(64, 344, "- - input S(1.5)|0⟩, q anti-squeezed (× 6)", GREY)
P.text(64, 362, f"— output, {keys[0].replace('y=', 'y = ')}: peaks every π/θ = √π", YELLOW)
P.text(64, 380, f"— output, {keys[1].replace('y=', 'y = ')} (a sampled outcome): pairs, gaps 2 arccos(y/2β)/θ and the rest", GREEN)
s.append("</svg>")
(D / "heralded-comb.svg").write_text("\n".join(s))

# ── Wigner functions ──
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
W = json.load(open(D / "cv-wigner.json"))
cmap = LinearSegmentedColormap.from_list("w", [PINK, BG, BLUE])
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.6), dpi=150, facecolor=BG)
for ax, key, title in ((axs[0], "cat_odd_beta2", "odd cat, β = 2"), (axs[1], "gkp0_6dB", "GKP |0̃⟩, 6 dB, Δ² = 0.05")):
    m = max(abs(min(min(r) for r in W[key])), max(max(r) for r in W[key]))
    im = ax.imshow(W[key], origin="lower", extent=(W["q"][0], W["q"][-1], W["p"][0], W["p"][-1]), cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-m, vcenter=0, vmax=m))
    ax.set_facecolor(BG); ax.set_title(title, color=LAB, fontsize=9, family="monospace")
    ax.set_xlabel("q", color=LAB, fontsize=8); ax.set_ylabel("p", color=LAB, fontsize=8)
    ax.tick_params(colors=TICK, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(colors=TICK, labelsize=6); cb.outline.set_edgecolor(GRID)
fig.tight_layout()
fig.savefig(D / "cv-wigner.png", facecolor=BG)
print("charts written")
