"""Phase 6: the code-capacity erasure sweeps — QxDecode (on the BEAM) against PyMatching (erased edges at weight 0,
shot by shot), on the same samples. (a) pure erasure, d = 3..9, e = 0.30..0.60: the curves should cross near 0.5;
(b) erasure plus bit flips, e ∈ {0, 0.1, 0.2}, p swept: the Pauli threshold as the erased share grows.
    ~/quanterall/.venv-tools/bin/python tools/erasure/sweep.py WORK_DIR site/data"""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
work, out = Path(sys.argv[1]), Path(sys.argv[2])
work.mkdir(parents=True, exist_ok=True)
SHOTS = 20000


def run(d, e, p, seed):
    subprocess.run([PY, str(HERE / "codecap.py"), "sample", str(d), str(e), str(p), str(SHOTS), str(seed), str(work)], check=True)
    subprocess.run(["docker", "cp", f"{work}/.", "quanterall_quantum:/tmp/erasure/"], check=True)
    o = subprocess.run(["docker", "exec", "-w", "/app", "-e", "MIX_ENV=prod", "quanterall_quantum", "mix", "run", "--no-start",
                        "/tmp/erasure/decode_erasure.exs", f"/tmp/erasure/g_d{d}.json", f"/tmp/erasure/s_d{d}_e{e}_p{p}.bin", str(SHOTS)],
                       capture_output=True, text=True).stdout
    uf = json.loads([l for l in o.splitlines() if l.startswith("{")][-1])
    pm = json.loads(subprocess.run([PY, str(HERE / "codecap.py"), "pymatching", str(d), str(e), str(p), str(SHOTS), str(seed)],
                                   capture_output=True, text=True).stdout.strip().splitlines()[-1])
    f1, f2 = uf["failures"] / SHOTS, pm["pymatching_failures"] / SHOTS
    row = {"d": d, "erasure": e, "flip": p, "shots": SHOTS, "seed": seed, "qxdecode_failures": uf["failures"],
           "pymatching_failures": pm["pymatching_failures"], "qxdecode_logical": f1, "pymatching_logical": f2,
           "qxdecode_err": math.sqrt(max(f1 * (1 - f1), 1e-12) / SHOTS), "pymatching_err": math.sqrt(max(f2 * (1 - f2), 1e-12) / SHOTS)}
    print(f"d={d} e={e} p={p}: QxDecode {f1:.4f}, PyMatching {f2:.4f}", flush=True)
    return row


rows = []
for d in (3, 5, 7, 9):
    for i, e in enumerate((0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)):
        rows.append(run(d, e, 0.0, 1000 * d + i))
with open(out / "erasure-pure.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
mixed = []
for e in (0.0, 0.1, 0.2):
    for d in (3, 5, 7):
        for i, p in enumerate((0.02, 0.04, 0.06, 0.08, 0.10, 0.12)):
            mixed.append(run(d, e, p, 50000 + 1000 * d + 100 * int(e * 10) + i))
with open(out / "erasure-mixed.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(mixed[0])); w.writeheader(); w.writerows(mixed)
