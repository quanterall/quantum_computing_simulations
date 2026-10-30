"""The decoder load test, end to end (quanterall.com, The research question).

Every decoder on the same 8 fast cores (5–8, 15–18 of a 20-core Arm machine; the rest left to the machine's
other services): QxDecode on the BEAM (8 schedulers, a process per syndrome), union-find in C++ (a dispatcher and
7 worker threads), union-find in Python and PyMatching (a dispatcher and 7 worker processes). Same graph, same
syndromes (d = 3, 5, 7 at p = 0.4 %).

  capacity     the decode-only throughput (reported), and the capacity the loads are set from: the throughput
               through the decoder's own harness with 10,000 syndromes outstanding — what its dispatcher and
               workers can actually sustain (a decode-only figure would ask the dispatcher for more arrivals
               than it can release, and measure the dispatcher)
  open         Poisson arrivals at 50, 80 and 95 % of that capacity; 10 seeds at d = 5, 3 at d = 3 and 7;
               p50, p99, p99.9 of completion − arrival; the share over the 1 ms deadline
  noop         the same harness with a decoder that does nothing, at the same rates: the harness's own cost
  closed       1, 100 and 10,000 syndromes always outstanding (the page's first design; by Little's law its
               mean latency is in-flight ÷ throughput)

    ~/quanterall/.venv-tools/bin/python tools/load/run_all.py LOAD_DIR site/data
"""
import json
import random
import statistics
import subprocess
import sys
from pathlib import Path

LOAD, OUT = Path(sys.argv[1]), Path(sys.argv[2])
CORES = "5-8,15-18"
PY = str(Path(sys.executable))
HERE = Path(__file__).resolve().parent


def run(cmd):
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=1800).stdout
    lines = [l for l in out.splitlines() if l.startswith("{")]
    if not lines:
        raise RuntimeError(f"no result from {' '.join(cmd)}:\n{out[-2000:]}")
    return json.loads(lines[-1])


def cmd(decoder, d, *args):
    dd = LOAD / f"d{d}"
    if decoder == "cpp":
        return ["taskset", "-c", CORES, str(LOAD / "uf"), str(dd), *map(str, args)]
    if decoder == "beam":
        return ["docker", "exec", "-w", "/app", "-e", "MIX_ENV=prod", "-e", "ELIXIR_ERL_OPTIONS=+S 8:8", "quanterall_quantum",
                "taskset", "-c", CORES, "mix", "run", "--no-start", "/tmp/load/beam_load.exs", f"/tmp/load/d{d}", *map(str, args)]
    return ["taskset", "-c", CORES, PY, str(HERE / "load_py.py"), decoder, str(dd), *map(str, args)]


def mode_args(decoder, mode, *a, noop=False):
    """the harnesses take slightly different argument lists; workers: 7 (+ the dispatcher) where they take them"""
    if mode == "capacity":
        base = ["capacity"] if decoder == "beam" else ["capacity", 7]
    elif mode == "open":
        base = ["open", float(a[0]), a[1], a[2]] if decoder == "beam" else ["open", a[0], a[1], 7, a[2]]
    else:
        base = ["closed", a[0], a[1]] if decoder == "beam" else ["closed", a[0], a[1], 7]
    return base + (["noop"] if noop and decoder in ("cpp", "beam") else [])


results = []
for d in (3, 5, 7):
    seeds = range(1, 11) if d == 5 else range(1, 4)
    for decoder in ("cpp", "beam", "pymatching", "pyuf"):
        raw = run(cmd(decoder, d, *mode_args(decoder, "capacity")))["throughput"]
        closed = {}
        for inflight in (1, 100, 10_000):
            count = 50_000 if decoder in ("cpp", "beam") else 20_000
            r = run(cmd(decoder, d, *mode_args(decoder, "closed", inflight, count)))
            r.update({"d": d, "decoder": decoder, "inflight": inflight}); results.append(r); closed[inflight] = r
            print(f"   closed {inflight}: p99 {r['p99_us']:,.0f} µs, throughput {r['throughput']:,.0f}/s", flush=True)
        cap = closed[10_000]["throughput"]
        print(f"d={d} {decoder}: decode-only {raw:,.0f}/s; through the harness {cap:,.0f}/s", flush=True)
        results.append({"d": d, "decoder": decoder, "mode": "capacity", "throughput": cap, "decode_only": raw})
        for load in (0.5, 0.8, 0.95):
            rate = round(cap * load, 1)
            count = int(min(200_000, max(20_000, rate * 1.0)))
            if decoder in ("pymatching", "pyuf"):
                count = int(min(count, 40_000))
            p99s = []
            for seed in seeds:
                r = run(cmd(decoder, d, *mode_args(decoder, "open", rate, count, seed)))
                r.update({"d": d, "decoder": decoder, "load": load, "seed": seed})
                results.append(r); p99s.append(r["p99_us"])
            if decoder in ("cpp", "beam"):
                r = run(cmd(decoder, d, *mode_args(decoder, "open", rate, count, 1, noop=True)))
            else:
                r = run(["taskset", "-c", CORES, PY, str(HERE / "load_py.py"), "noop", str(LOAD / f"d{d}"), "open", str(rate), str(count), "7", "1"])
            r.update({"d": d, "decoder": decoder, "load": load, "noop": True}); results.append(r)
            print(f"   load {load:.0%} ({rate:,.0f}/s): p99 median {statistics.median(p99s):,.0f} µs over {len(p99s)} seeds; no-op p99 {r['p99_us']:,.0f} µs", flush=True)
        json.dump(results, open(OUT / "decoder-load-raw.json", "w"))

# the summary: per decoder, d and load, the median p99 over seeds with a bootstrap 95 % interval
rng = random.Random(7)
summary = []
for d in (3, 5, 7):
    for decoder in ("cpp", "beam", "pymatching", "pyuf"):
        cap = next(r["throughput"] for r in results if r["d"] == d and r["decoder"] == decoder and r["mode"] == "capacity")
        dec_only = next(r["decode_only"] for r in results if r["d"] == d and r["decoder"] == decoder and r["mode"] == "capacity")
        row = {"d": d, "decoder": decoder, "capacity": round(cap), "decode_only": round(dec_only)}
        for load in (0.5, 0.8, 0.95):
            runs = [r for r in results if r["d"] == d and r["decoder"] == decoder and r.get("load") == load and not r.get("noop")]
            p99 = [r["p99_us"] for r in runs]
            boots = sorted(statistics.median(rng.choices(p99, k=len(p99))) for _ in range(2000))
            noop = next(r for r in results if r["d"] == d and r["decoder"] == decoder and r.get("load") == load and r.get("noop"))
            row[f"p99_us_{int(load*100)}"] = round(statistics.median(p99), 1)
            row[f"p99_ci_{int(load*100)}"] = [round(boots[50], 1), round(boots[1949], 1)]
            row[f"p50_us_{int(load*100)}"] = round(statistics.median(r["p50_us"] for r in runs), 1)
            row[f"over_1ms_{int(load*100)}"] = round(statistics.median(r["over_1ms"] for r in runs), 5)
            row[f"noop_p99_us_{int(load*100)}"] = round(noop["p99_us"], 1)
        for inflight in (1, 100, 10_000):
            r = next(r for r in results if r["d"] == d and r["decoder"] == decoder and r.get("inflight") == inflight)
            row[f"closed_{inflight}_p99_us"] = round(r["p99_us"], 1)
        summary.append(row)
json.dump({"summary": summary, "cores": CORES, "note": __doc__}, open(OUT / "decoder-load.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
