"""The decoding graph for our decoder, from Stim's detector error model (the page says so: there we depend on Stim).

Each error mechanism of the circuit, decomposed into graph-like parts (one or two detectors, and whether it flips
the logical observable), becomes an edge; a part with one detector is an edge to the boundary. Parallel edges are
merged (their probabilities combined as independent flips). Written as JSON: the number of detectors and
[u, v, p, flips_observable] per edge, v = -1 for the boundary.

    ~/quanterall/.venv-tools/bin/python tools/dem_graph.py CIRCUIT.stim OUT.json
"""
import json
import sys

import stim

circuit = stim.Circuit(open(sys.argv[1]).read())
dem = circuit.detector_error_model(decompose_errors=True)
edges = {}


def add(dets, obs, p):
    if not dets:
        return                                                  # touches no detector: invisible to any decoder
    u, v = (dets[0], -1) if len(dets) == 1 else tuple(sorted(dets))
    key = (u, v)
    if key in edges:
        q, o = edges[key]
        # two independent flips on the same edge: it flips if exactly one does; keep the likelier observable effect
        edges[key] = (q * (1 - p) + p * (1 - q), o if q >= p else obs)
    else:
        edges[key] = (p, obs)


for inst in dem.flattened():
    if inst.type != "error":
        continue
    p = inst.args_copy()[0]
    dets, obs = [], False
    for t in inst.targets_copy():
        if t.is_separator():
            add(dets, obs, p)
            dets, obs = [], False
        elif t.is_relative_detector_id():
            dets.append(t.val)
        elif t.is_logical_observable_id():
            obs = not obs
    add(dets, obs, p)

out = {"detectors": dem.num_detectors, "edges": [[u, v, p, int(o)] for (u, v), (p, o) in sorted(edges.items())]}
json.dump(out, open(sys.argv[2], "w"))
print(f"{len(out['edges'])} edges over {out['detectors']} detectors ({sum(1 for e in out['edges'] if e[1] == -1)} to the boundary)")
