"""The HOM-dip chart for the site, from site/data/hom-dip.csv (tools/optics/klm_checks.py)."""
import csv, sys
from pathlib import Path
rows = list(csv.DictReader(open(Path(sys.argv[1]) / "hom-dip.csv")))
W, H, L, R, T, B = 760, 360, 64, 20, 20, 46
x = lambda g: L + (float(g) + 4) / 8 * (W - L - R)
y = lambda p: T + (0.55 - float(p)) / 0.55 * (H - T - B)
s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Hong–Ou–Mandel dip: coincidence probability against delay, simulated and the formula one half times one minus e to the minus gamma tau" font-family="IBM Plex Mono, monospace" font-size="11">',
     f'<rect width="{W}" height="{H}" fill="#0f1726"/>']
for p in (0, 0.1, 0.2, 0.3, 0.4, 0.5):
    s.append(f'<line x1="{L}" y1="{y(p):.1f}" x2="{W-R}" y2="{y(p):.1f}" stroke="#1e2a40"/><text x="{L-8}" y="{y(p)+4:.1f}" fill="#7a879c" text-anchor="end">{p:g}</text>')
for g in (-4, -2, 0, 2, 4):
    s.append(f'<line x1="{x(g):.1f}" y1="{T}" x2="{x(g):.1f}" y2="{H-B}" stroke="#1e2a40"/><text x="{x(g):.1f}" y="{H-B+16}" fill="#7a879c" text-anchor="middle">{g}</text>')
s.append(f'<text x="{(L+W-R)/2}" y="{H-10}" fill="#9aa6bb" text-anchor="middle">delay γτ</text>')
s.append(f'<text x="14" y="{(T+H-B)/2}" fill="#9aa6bb" text-anchor="middle" transform="rotate(-90 14 {(T+H-B)/2})">coincidence probability</text>')
s.append('<polyline points="' + " ".join(f"{x(r['gamma_tau']):.1f},{y(r['formula']):.1f}" for r in rows) + '" fill="none" stroke="#f5c842" stroke-width="2"/>')
for r in rows[::2]:
    s.append(f'<circle cx="{x(r["gamma_tau"]):.1f}" cy="{y(r["coincidence"]):.1f}" r="3" fill="none" stroke="#6ab0ff"/>')
s.append(f'<text x="{x(1.2):.1f}" y="{y(0.14):.1f}" fill="#f5c842">— ½(1 − e^(−γ|τ|))</text><text x="{x(1.2):.1f}" y="{y(0.14)+18:.1f}" fill="#6ab0ff">○ simulated (QxFock, fock.py)</text>')
s.append("</svg>")
(Path(sys.argv[1]) / "hom-dip.svg").write_text("\n".join(s))
