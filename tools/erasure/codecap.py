"""Phase 6a/6b: erasure decoding at code capacity — the rotated surface code, perfect syndrome measurement.

Each data qubit is erased with probability e (it then holds a random Pauli: an X component with probability ½,
the herald says which qubit) and otherwise flips (X) with probability p. Memory in Z: the Z checks see the X
components; the logical Z is a row of data qubits. The decoding graph has one edge per data qubit (two Z checks,
or one and the boundary) — so an erased qubit is an erased edge.

    codecap.py graph D OUT.json               the graph, in tools/dem_graph.py's format (one edge per qubit)
    codecap.py sample D E P SHOTS SEED DIR    samples: syndromes, erased edges, the logical flip (binary)
    codecap.py pymatching D E P DIR           PyMatching with the erased edges' weights set to 0, shot by shot
"""
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np


def layout(d):
    """data qubits (row-major), Z checks (each a list of data qubits), the logical Z row — as QxStab"""
    data = [(2 * x + 1, 2 * y + 1) for y in range(d) for x in range(d)]
    di = {q: i for i, q in enumerate(data)}
    z = []
    for y in range(d + 1):
        for x in range(d + 1):
            mx, my, kind = 2 * x, 2 * y, "x" if (x + y) % 2 == 0 else "z"
            nb = [di[q] for q in ((mx + 1, my + 1), (mx + 1, my - 1), (mx - 1, my + 1), (mx - 1, my - 1)) if q in di]
            if kind == "z" and (len(nb) == 4 or (len(nb) == 2 and mx in (0, 2 * d))):
                z.append(nb)
    xs = []
    for y in range(d + 1):
        for x in range(d + 1):
            mx, my, kind = 2 * x, 2 * y, "x" if (x + y) % 2 == 0 else "z"
            nb = [di[q] for q in ((mx + 1, my + 1), (mx + 1, my - 1), (mx - 1, my + 1), (mx - 1, my - 1)) if q in di]
            if kind == "x" and (len(nb) == 4 or (len(nb) == 2 and my in (0, 2 * d))):
                xs.append(nb)
    for y in range(d):
        row = [di[(2 * x + 1, 2 * y + 1)] for x in range(d)]
        if all(sum(q in row for q in c) % 2 == 0 for c in xs):
            return data, z, set(row)
    raise RuntimeError("no logical row")


def graph(d, p=0.01):
    data, z, obs = layout(d)
    n = len(z)
    edges = []
    for q in range(len(data)):
        ch = [i for i, c in enumerate(z) if q in c]
        u, v = (ch[0], -1) if len(ch) == 1 else (ch[0], ch[1])
        edges.append([u, v, p, int(q in obs)])
    return {"detectors": n, "edges": edges}


def sample(d, e, p, shots, seed):
    g = graph(d)
    n, E = g["detectors"], len(g["edges"])
    rng = np.random.default_rng(seed)
    erased = rng.random((shots, E)) < e
    xerr = np.where(erased, rng.random((shots, E)) < 0.5, rng.random((shots, E)) < p)
    syn = np.zeros((shots, n), dtype=bool)
    for q, (u, v, _, _) in enumerate(g["edges"]):
        syn[:, u] ^= xerr[:, q]
        if v >= 0:
            syn[:, v] ^= xerr[:, q]
    obs = np.zeros(shots, dtype=bool)
    for q, (_, _, _, o) in enumerate(g["edges"]):
        if o:
            obs ^= xerr[:, q]
    return g, syn, erased, obs


def write(path, syn, erased, obs):
    with open(path, "wb") as f:
        for k in range(len(obs)):
            a, b = np.where(syn[k])[0], np.where(erased[k])[0]
            f.write(struct.pack(f"<H{len(a)}H", len(a), *a))
            f.write(struct.pack(f"<H{len(b)}H", len(b), *b))
        f.write(obs.astype(np.uint8).tobytes())


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "graph":
        g = graph(int(sys.argv[2]))
        json.dump(g, open(sys.argv[3], "w"))
        print(g["detectors"], "checks,", len(g["edges"]), "edges")
    elif cmd == "sample":
        d, e, p, shots, seed, out = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), Path(sys.argv[7])
        out.mkdir(parents=True, exist_ok=True)
        g, syn, erased, obs = sample(d, e, p, shots, seed)
        json.dump(graph(d, max(p, 1e-6)), open(out / f"g_d{d}.json", "w"))
        write(out / f"s_d{d}_e{e}_p{p}.bin", syn, erased, obs)
    elif cmd == "pymatching":
        import pymatching
        import scipy.sparse as sp
        d, e, p, shots, seed = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        g, syn, erased, obs = sample(d, e, p, shots, seed)
        n, E = g["detectors"], len(g["edges"])
        rows, cols = [], []
        for q, (u, v, _, _) in enumerate(g["edges"]):
            rows.append(u); cols.append(q)
            if v >= 0:
                rows.append(v); cols.append(q)
        H = sp.csc_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(n, E))
        L = np.array([o for *_, o in g["edges"]], dtype=np.uint8)
        w_p = math.log((1 - p) / p) if p > 0 else 1e4
        fails = 0
        for k in range(shots):
            w = np.where(erased[k], 0.0, w_p)
            m = pymatching.Matching(H, weights=w)
            corr = m.decode(syn[k].astype(np.uint8))
            fails += int((corr @ L) % 2 != obs[k])
        print(json.dumps({"d": d, "e": e, "p": p, "shots": shots, "pymatching_failures": fails}))
