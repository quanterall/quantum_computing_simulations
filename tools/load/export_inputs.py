"""Shared inputs for the decoder load test (tools/load/): every decoder reads the same graph and syndromes.

  graph.txt      line 1: detectors edges; then per edge: u v length flips_observable weight (v = detectors
                 for the boundary). length = max(4, floor(4·w/w_min + ½)), w = ln((1−p)/p) — as QxDecode
  syndromes.bin  per shot: uint16 count, then count uint16 detector ids; after all shots, one byte per shot:
                 the observable's actual flip

    ~/quanterall/.venv-tools/bin/python tools/load/export_inputs.py GRAPH.json SAMPLE.bin SHOTS OUT_DIR
"""
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np

g = json.load(open(sys.argv[1]))
shots, out = int(sys.argv[3]), Path(sys.argv[4])
out.mkdir(parents=True, exist_ok=True)
n = g["detectors"]
ws = [math.log((1 - p) / p) for _, _, p, _ in g["edges"]]
wmin = min(ws)
with open(out / "graph.txt", "w") as f:
    f.write(f"{n} {len(g['edges'])}\n")
    for (u, v, p, o), w in zip(g["edges"], ws):
        f.write(f"{u} {n if v < 0 else v} {max(4, math.floor(4 * w / wmin + 0.5))} {o} {w!r}\n")
raw = np.fromfile(sys.argv[2], dtype=np.uint8).reshape(n + 1, -1)
bits = np.unpackbits(raw, axis=1, bitorder="little")[:, :shots].T.astype(bool)
with open(out / "syndromes.bin", "wb") as f:
    for k in range(shots):
        ids = np.where(bits[k, :n])[0]
        f.write(struct.pack(f"<H{len(ids)}H", len(ids), *ids))
    f.write(bits[:, n].astype(np.uint8).tobytes())
print(f"{n} detectors, {len(g['edges'])} edges, {shots} syndromes, mean {bits[:, :n].sum() / shots:.2f} fired")
