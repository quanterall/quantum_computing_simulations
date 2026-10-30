"""Phase 6c: decoding circuit-level heralded erasure — QxDecode (BEAM) with the heralds, PyMatching with the heralds
(the erased edges reweighted to 0, shot by shot), and PyMatching ignoring them (erasure treated as Pauli noise of
its probability), on the same samples from QxStab (tools/erasure/circuit_runs.exs).

The decoding graph comes from Stim's detector error model (with approximate_disjoint_errors, which HERALDED_ERASE
needs): the mechanisms without a herald are the ordinary edges; each herald's mechanisms give the edges its
erasure can flip — erased when the herald fires.
    ~/quanterall/.venv-tools/bin/python tools/erasure/circuit.py CL_DIR site/data"""
import csv
import json
import math
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import pymatching
import scipy.sparse as sp
import stim

CL, OUT = Path(sys.argv[1]), Path(sys.argv[2])


def graph_from_dem(c, heralds):
    dem = c.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    hs = set(heralds)
    keep = [i for i in range(c.num_detectors) if i not in hs]
    ren = {d: i for i, d in enumerate(keep)}
    base, er_p, herald_edges = {}, {}, {}

    def key(dets):
        dets = [ren[d] for d in dets]
        return (dets[0], -1) if len(dets) == 1 else tuple(sorted(dets))

    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        comps, cur, obs = [], [], False
        for t in inst.targets_copy():
            if t.is_separator():
                comps.append((cur, obs)); cur, obs = [], False
            elif t.is_relative_detector_id():
                cur.append(t.val)
            elif t.is_logical_observable_id():
                obs = not obs
        comps.append((cur, obs))
        hset = [d for cdets, _ in comps for d in cdets if d in hs]
        for cdets, o in comps:
            dets = [d for d in cdets if d not in hs]
            if not dets or len(dets) > 2:
                continue
            k = key(dets)
            if hset:
                for h in hset:
                    herald_edges.setdefault(h, set()).add(k)
                q = er_p.get(k, (0.0, o))[0]
                er_p[k] = (q * (1 - p) + p * (1 - q), er_p.get(k, (0.0, o))[1])
            else:
                q, oo = base.get(k, (0.0, o))
                base[k] = (q * (1 - p) + p * (1 - q), oo if q >= p else o)
    keys = sorted(set(base) | set(er_p))
    idx = {k: i for i, k in enumerate(keys)}
    edges = []
    for k in keys:
        pb, ob = base.get(k, (0.0, er_p.get(k, (0, False))[1]))
        edges.append([k[0], k[1], max(pb, 1e-9), int(ob)])
    return len(keep), keep, edges, idx, {h: sorted(idx[k] for k in ks) for h, ks in herald_edges.items()}, er_p


rows = []
for stim_path in sorted(CL.glob("cl_d*_pe*.stim")):
    name = stim_path.stem
    d, pe = int(name.split("_d")[1].split("_")[0]), float(name.split("_pe")[1])
    c = stim.Circuit(stim_path.read_text())
    heralds = json.loads((CL / f"{name}.heralds.json").read_text())
    S = 20000
    raw = np.fromfile(CL / f"{name}.bin", dtype=np.uint8).reshape(c.num_detectors + 1, -1)
    bits = np.unpackbits(raw, axis=1, bitorder="little")[:, :S].T.astype(bool)
    det, obs = bits[:, :c.num_detectors], bits[:, c.num_detectors]
    n, keep, edges, idx, hedges, er_p = graph_from_dem(c, heralds)
    g = {"detectors": n, "edges": edges}
    (CL / f"{name}.graph.json").write_text(json.dumps(g))
    # samples in decode_erasure.exs's format
    with open(CL / f"{name}.samples.bin", "wb") as f:
        for k in range(S):
            fired = [i for i, dd in enumerate(keep) if det[k, dd]]
            er = sorted({e for h in heralds if det[k, h] for e in hedges.get(h, [])})
            f.write(struct.pack(f"<H{len(fired)}H", len(fired), *fired))
            f.write(struct.pack(f"<H{len(er)}H", len(er), *er))
        f.write(obs.astype(np.uint8).tobytes())
    subprocess.run(["docker", "cp", str(CL / f"{name}.graph.json"), f"quanterall_quantum:/tmp/erasure/cl/{name}.graph.json"], check=True)
    subprocess.run(["docker", "cp", str(CL / f"{name}.samples.bin"), f"quanterall_quantum:/tmp/erasure/cl/{name}.samples.bin"], check=True)
    o = subprocess.run(["docker", "exec", "-w", "/app", "-e", "MIX_ENV=prod", "quanterall_quantum", "mix", "run", "--no-start",
                        "/tmp/erasure/decode_erasure.exs", f"/tmp/erasure/cl/{name}.graph.json", f"/tmp/erasure/cl/{name}.samples.bin", str(S)],
                       capture_output=True, text=True).stdout
    uf = json.loads([l for l in o.splitlines() if l.startswith("{")][-1])["failures"]
    # PyMatching: the check matrix of the graph; with heralds, per-shot weights; without, erasure as Pauli noise
    E = len(edges)
    r_, c_ = [], []
    for j, (u, v, _, _) in enumerate(edges):
        r_.append(u); c_.append(j)
        if v >= 0:
            r_.append(v); c_.append(j)
    H = sp.csc_matrix((np.ones(len(r_), dtype=np.uint8), (r_, c_)), shape=(n, E))
    Lv = np.array([o for *_, o in edges], dtype=np.uint8)
    w0 = np.array([math.log((1 - p) / p) for *_, p, _ in [(u, v, p, o) for u, v, p, o in edges]])
    syn = det[:, keep].astype(np.uint8)
    fail_h = 0
    for k in range(S):
        w = w0.copy()
        for h in heralds:
            if det[k, h]:
                w[hedges.get(h, [])] = 0.0
        corr = pymatching.Matching(H, weights=w).decode(syn[k])
        fail_h += int((corr @ Lv) % 2 != obs[k])
    p_tot = np.array([1 - (1 - p) * (1 - er_p.get(k, (0.0, 0))[0]) for (k, (u, v, p, o)) in zip(sorted(idx, key=idx.get), edges)])
    m = pymatching.Matching(H, weights=np.log((1 - p_tot) / p_tot))
    fail_n = int(np.sum((m.decode_batch(syn) @ Lv) % 2 != obs))
    row = {"d": d, "erasure": pe, "pauli": 0.001, "rounds": d, "shots": S, "heralds_fired_per_shot": float(det[:, heralds].sum() / S) if heralds else 0.0,
           "qxdecode_with_heralds": uf / S, "pymatching_with_heralds": fail_h / S, "pymatching_ignoring_heralds": fail_n / S}
    rows.append(row)
    print(f"d={d} pe={pe}: QxDecode+heralds {uf / S:.4f}, PyMatching+heralds {fail_h / S:.4f}, PyMatching ignoring {fail_n / S:.4f}", flush=True)
with open(OUT / "erasure-circuit.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
