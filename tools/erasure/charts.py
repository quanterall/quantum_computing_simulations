"""Charts for the erasure sweeps: site/data/erasure-pure.svg and erasure-mixed.svg (from the CSVs)."""
import csv
import sys
from pathlib import Path

D = Path(sys.argv[1])
col = {3: "#6ab0ff", 5: "#a8ff3e", 7: "#f472b6", 9: "#f5c842"}


def chart(path, rows, xkey, xlabel, xr, aria, legend):
    W, H, L, R, T, B = 760, 380, 64, 20, 20, 46
    x = lambda v: L + (float(v) - xr[0]) / (xr[1] - xr[0]) * (W - L - R)
    y = lambda v: T + (0.45 - float(v)) / 0.45 * (H - T - B)
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="{aria}" font-family="IBM Plex Mono, monospace" font-size="11">',
         f'<rect width="{W}" height="{H}" fill="#0f1726"/>']
    for v in (0, 0.1, 0.2, 0.3, 0.4):
        s.append(f'<line x1="{L}" y1="{y(v):.1f}" x2="{W - R}" y2="{y(v):.1f}" stroke="#1e2a40"/><text x="{L - 8}" y="{y(v) + 4:.1f}" fill="#7a879c" text-anchor="end">{v:g}</text>')
    for v in sorted({float(r[xkey]) for r in rows}):
        s.append(f'<line x1="{x(v):.1f}" y1="{T}" x2="{x(v):.1f}" y2="{H - B}" stroke="#1e2a40"/><text x="{x(v):.1f}" y="{H - B + 16}" fill="#7a879c" text-anchor="middle">{v:g}</text>')
    s.append(f'<text x="{(L + W - R) / 2}" y="{H - 10}" fill="#9aa6bb" text-anchor="middle">{xlabel}</text>')
    s.append(f'<text x="14" y="{(T + H - B) / 2}" fill="#9aa6bb" text-anchor="middle" transform="rotate(-90 14 {(T + H - B) / 2})">logical error</text>')
    groups = {}
    for r in rows:
        groups.setdefault(r["group"], []).append(r)
    for rs in groups.values():
        rs.sort(key=lambda r: float(r[xkey]))
        c = col[int(rs[0]["d"])]
        dash = ' stroke-dasharray="5 4"' if rs[0].get("dash") else ""
        s.append(f'<polyline points="{" ".join(f"{x(r[xkey]):.1f},{y(r["qxdecode_logical"]):.1f}" for r in rs)}" fill="none" stroke="{c}" stroke-width="2"{dash}/>')
        for r in rs:
            s.append(f'<circle cx="{x(r[xkey]):.1f}" cy="{y(r["pymatching_logical"]):.1f}" r="3.5" fill="none" stroke="{c}"/>')
    for i, (lab, c, dash) in enumerate(legend):
        s.append(f'<line x1="{L + 16}" y1="{T + 14 + i * 16}" x2="{L + 36}" y2="{T + 14 + i * 16}" stroke="{c}" stroke-width="2"{dash}/><text x="{L + 42}" y="{T + 18 + i * 16}" fill="#c8d4e8">{lab}</text>')
    s.append(f'<text x="{L + 16}" y="{T + 18 + len(legend) * 16 + 2}" fill="#7a879c">— QxDecode   ○ PyMatching, erased edges at weight 0</text>')
    s.append("</svg>")
    (D / path).write_text("\n".join(s))


pure = list(csv.DictReader(open(D / "erasure-pure.csv")))
for r in pure:
    r["group"] = r["d"]
chart("erasure-pure.svg", pure, "erasure", "erasure probability per qubit (pure erasure, perfect syndrome)", (0.3, 0.6),
      "Surface code under pure erasure: logical error against erasure probability, d = 3, 5, 7, 9; the curves cross at one half",
      [(f"d = {d}", col[d], "") for d in (3, 5, 7, 9)])
mixed = list(csv.DictReader(open(D / "erasure-mixed.csv")))
sel = [r for r in mixed if r["erasure"] in ("0.0", "0.2")]
for r in sel:
    r["group"] = r["d"] + "_" + r["erasure"]
    r["dash"] = r["erasure"] == "0.2"
chart("erasure-mixed.svg", sel, "flip", "bit-flip probability on the other qubits", (0.02, 0.12),
      "Surface code with erasure and bit flips: logical error against the flip probability, without erasure and with 20 percent erasure",
      [(f"d = {d}, no erasure", col[d], "") for d in (3, 5, 7)] + [(f"d = {d}, 20 % erasure", col[d], ' stroke-dasharray="5 4"') for d in (3, 5, 7)])
print("charts written")
