"""QxFock (Elixir, the BEAM service) against fock.py (Python) on random interferometers, every output amplitude,
both of QxFock's routes. Writes cases, runs them in the container, compares.
    ~/quanterall/.venv-tools/bin/python tools/optics/xcheck_beam.py"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim.fock import FockState, random_unitary  # noqa: E402

rng = np.random.default_rng(7)
cases = []
for m in range(4, 9):
    for n in range(1, 5):
        U = random_unitary(m, rng)
        t = [0] * m
        for k in rng.choice(m, size=n):
            t[k] += 1
        amps = FockState.basis(tuple(t)).apply(U).amps
        cases.append({"m": m, "t": t, "U": [[[float(z.real), float(z.imag)] for z in row] for row in U],
                      "py": {",".join(map(str, k)): [v.real, v.imag] for k, v in amps.items()}})
Path("/tmp/fockx.json").write_text(json.dumps(cases))
subprocess.run(["docker", "cp", "/tmp/fockx.json", "quanterall_quantum:/tmp/fockx.json"], check=True)
code = r'''
cases = File.read!("/tmp/fockx.json") |> JSON.decode!()
worst = Enum.reduce(cases, 0.0, fn c, w ->
  u = Enum.map(c["U"], fn row -> Enum.map(row, fn [a, b] -> {a, b} end) end)
  st = QxFock.basis(List.to_tuple(c["t"]))
  {a, p} = {QxFock.linear(st, u), QxFock.apply_perm(st, u)}
  py = Map.new(c["py"], fn {k, [x, y]} -> {k |> String.split(",") |> Enum.map(&String.to_integer/1) |> List.to_tuple(), {x, y}} end)
  keys = Enum.uniq(Map.keys(py) ++ Map.keys(a.amps))
  d = keys |> Enum.map(fn k ->
    {x, y} = Map.get(a.amps, k, {0.0, 0.0}); {x2, y2} = Map.get(py, k, {0.0, 0.0}); {x3, y3} = Map.get(p.amps, k, {0.0, 0.0})
    Enum.max([abs(x - x2), abs(y - y2), abs(x - x3), abs(y - y3)]) end) |> Enum.max()
  max(w, d) end)
IO.puts(JSON.encode!(%{cases: length(cases), max_diff: worst}))
'''
out = subprocess.run(["docker", "exec", "-w", "/app", "-e", "MIX_ENV=prod", "quanterall_quantum", "mix", "run", "--no-start", "-e", code],
                     capture_output=True, text=True).stdout
print([l for l in out.splitlines() if l.startswith("{")][-1])
