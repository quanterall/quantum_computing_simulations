"""The Python side of the decoder load test: PyMatching, and union-find in plain Python (the same algorithm as
QxDecode and tools/load/uf.cpp), under the same load as the others.

A dispatcher (this process) releases syndromes at their arrival times into a queue; worker processes, each with
its own decoder, take them; latency = completion − scheduled arrival (open loop) or − submission (closed).
Python's queues cost tens of microseconds per item: a no-op decoder measures that overhead.

    python load_py.py DECODER DIR capacity WORKERS
    python load_py.py DECODER DIR open RATE COUNT WORKERS SEED
    python load_py.py DECODER DIR closed INFLIGHT COUNT WORKERS
    python load_py.py pyuf DIR predict          (decode all in one process; logical error and µs per syndrome)
DECODER: pymatching | pyuf | noop
"""
import heapq
import json
import math
import multiprocessing as mp
import struct
import sys
import time
from pathlib import Path


def load_syndromes(path):
    b = Path(path).read_bytes()
    recs, i = [], 0
    while len(b) - i != len(recs):
        (c,) = struct.unpack_from("<H", b, i)
        recs.append(list(struct.unpack_from(f"<{c}H", b, i + 2)))
        i += 2 + 2 * c
    return recs, list(b[i:])


class PyUF:
    """union-find, weighted, peeled over a minimum-weight spanning tree per cluster — as QxDecode"""

    def __init__(self, path):
        lines = Path(path).read_text().split("\n")
        self.n, m = map(int, lines[0].split())
        self.edges = []
        self.adj = [[] for _ in range(self.n + 1)]
        for e in range(m):
            u, v, ln, o, w = lines[1 + e].split()
            u, v, ln, o, w = int(u), int(v), int(ln), int(o), float(w)
            self.edges.append((u, v, ln, o, w))
            self.adj[u].append(e)
            if v != self.n:
                self.adj[v].append(e)

    def decode(self, fired):
        if not fired:
            return 0
        B, E = self.n, self.edges
        parent = {v: v for v in fired}
        cl = {v: [True, False, [v]] for v in fired}                      # odd, boundary, vertices
        support = {}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v

        while True:
            active = [r for r, c in cl.items() if c[0] and not c[1]]
            if not active:
                break
            fused = []
            for r in active:
                for v in cl[r][2]:
                    for e in self.adj[v]:
                        s = support.get(e, 0)
                        if s < E[e][2]:
                            support[e] = s + 1
                            if s + 1 >= E[e][2]:
                                fused.append(e)
            for e in fused:
                u, v = E[e][0], E[e][1]
                if v == B:
                    cl[find(u)][1] = True
                    continue
                for x in (u, v):
                    if x not in parent:
                        parent[x] = x
                        cl[x] = [False, False, [x]]
                ru, rv = find(u), find(v)
                if ru == rv:
                    continue
                if len(cl[ru][2]) < len(cl[rv][2]):
                    ru, rv = rv, ru
                a, b = cl[ru], cl.pop(rv)
                parent[rv] = ru
                cl[ru] = [a[0] != b[0], a[1] or b[1], a[2] + b[2]]
        flips = 0
        defects = set(fired)
        for r, (odd, bnd, verts) in cl.items():
            inside = set(verts)
            cand = sorted({e for v in verts for e in self.adj[v]
                           if support.get(e, 0) >= E[e][2] and E[e][0] in inside and (E[e][1] == B or E[e][1] in inside)},
                          key=lambda e: E[e][4])
            dsu, tree = {}, {}

            def df(x):
                while dsu.get(x, x) != x:
                    x = dsu[x]
                return x

            for e in cand:
                a, b = df(E[e][0]), df(E[e][1])
                if a != b:
                    dsu[a] = b
                    tree.setdefault(E[e][0], []).append(e)
                    tree.setdefault(E[e][1], []).append(e)
            root = B if bnd else verts[0]
            order, q, seen, h = [], [root], {root}, 0
            while h < len(q):
                v = q[h]; h += 1
                for e in tree.get(v, []):
                    w = E[e][1] if E[e][0] == v else E[e][0]
                    if w not in seen:
                        seen.add(w); q.append(w); order.append((w, e, v))
            marks = {v for v in verts if v in defects}
            for v, e, par in reversed(order):
                if v in marks:
                    marks.discard(v)
                    flips ^= E[e][3]
                    if par != B:
                        marks.symmetric_difference_update({par})
        return flips


def make_decoder(kind, d):
    if kind == "pyuf":
        return PyUF(Path(d) / "graph.txt").decode
    if kind == "pymatching":
        import pymatching
        import stim
        c = stim.Circuit(open(Path(d) / "circuit.stim").read())
        m = pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True))
        nd = c.num_detectors
        import numpy as np

        def dec(fired):
            z = np.zeros(nd, dtype=np.uint8); z[fired] = 1
            return int(m.decode(z)[0])
        return dec
    return lambda fired: 0


def worker(kind, d, inq, outq):
    dec = make_decoder(kind, d)
    syn, _ = load_syndromes(Path(d) / "syndromes.bin")
    n = len(syn)
    outq.put(("ready", 0))
    while True:
        item = inq.get()
        if item is None:
            return
        k, _t = item
        dec(syn[k % n])
        outq.put((k, time.perf_counter_ns()))


def main():
    kind, d, mode = sys.argv[1], sys.argv[2], sys.argv[3]
    if mode == "predict":
        dec = make_decoder(kind, d)
        syn, obs = load_syndromes(Path(d) / "syndromes.bin")
        t0 = time.perf_counter()
        pred = [dec(s) for s in syn]
        secs = time.perf_counter() - t0
        print(json.dumps({"mode": "predict", "decoder": kind, "count": len(syn), "us_per_syndrome": secs / len(syn) * 1e6,
                          "logical": sum(p != o for p, o in zip(pred, obs)) / len(syn)}))
        return
    args = sys.argv[4:]
    W = int(args[-1]) if mode == "capacity" else int(args[2])
    inq, outq = mp.Queue(), mp.Queue()
    ps = [mp.Process(target=worker, args=(kind, d, inq, outq)) for _ in range(W)]
    for p in ps:
        p.start()
    for _ in ps:
        outq.get()                                                       # every worker loaded
    sent_at = {}
    if mode == "capacity":
        count = 20000
        t0 = time.perf_counter_ns()
        for k in range(count):
            inq.put((k, 0))
        for _ in range(count):
            outq.get()
        secs = (time.perf_counter_ns() - t0) / 1e9
        out = {"mode": "capacity", "decoder": kind, "workers": W, "count": count, "seconds": secs, "throughput": count / secs}
    else:
        count = int(args[1])
        lat = []
        t0 = time.perf_counter_ns()
        if mode == "open":
            import random
            rate, rng = float(args[0]), random.Random(int(args[3]))
            at = 0.0
            due = []
            for k in range(count):
                at += rng.expovariate(rate)
                due.append(t0 + int(at * 1e9))
            done = 0
            for k in range(count):
                while time.perf_counter_ns() < due[k]:
                    pass
                inq.put((k, due[k]))
                while not outq.empty():                                  # collect as we go
                    kk, tdone = outq.get()
                    lat.append((tdone - due[kk]) / 1e3); done += 1
            while done < count:
                kk, tdone = outq.get(); lat.append((tdone - due[kk]) / 1e3); done += 1
            extra = {"rate": rate}
        else:
            inflight = int(args[0])
            sent = 0
            for k in range(min(inflight, count)):
                sent_at[k] = time.perf_counter_ns(); inq.put((k, 0)); sent += 1
            for _ in range(count):
                kk, tdone = outq.get()
                lat.append((tdone - sent_at.pop(kk)) / 1e3)
                if sent < count:
                    sent_at[sent] = time.perf_counter_ns(); inq.put((sent, 0)); sent += 1
            extra = {"inflight": inflight}
        secs = (time.perf_counter_ns() - t0) / 1e9
        lat.sort()
        pct = lambda p: lat[min(len(lat) - 1, int(p * len(lat)))]
        out = {"mode": mode, "decoder": kind, "workers": W, "count": count, "seconds": secs, "throughput": count / secs,
               "p50_us": pct(0.5), "p99_us": pct(0.99), "p999_us": pct(0.999), "max_us": lat[-1],
               "over_1ms": sum(x > 1000 for x in lat) / len(lat), **extra}
    for _ in ps:
        inq.put(None)
    for p in ps:
        p.join()
    print(json.dumps(out))


if __name__ == "__main__":
    main()
