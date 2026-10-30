"""Phase 1 cross-check: our Fock module's two routes (creation operators; Ryser permanents) against each other and
against Perceval (Quandela), on Haar-random interferometers of 4–8 modes with 1–4 photons; every output amplitude.
Power check: one beamsplitter sign flipped in our matrix must be caught.
    ~/quanterall/.venv-tools/bin/python tools/optics/xcheck_linear.py site/data"""
import json
import sys
from pathlib import Path

import numpy as np
import perceval as pcvl

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "qsim"))
from qsim.fock import FockState, beamsplitter, random_unitary  # noqa: E402

backend = pcvl.BackendFactory.get_backend("SLOS")


def perceval_amps(M, t):
    """Perceval's amplitudes for input t through the unitary M (its convention established below)"""
    backend.set_circuit(pcvl.Unitary(pcvl.Matrix(M)))
    backend.set_input_state(pcvl.BasicState(list(t)))
    return backend


# 1. the convention: a non-symmetric 3-mode unitary, one photon in mode 0 — where does it go?
rng = np.random.default_rng(1)
U = random_unitary(3, rng)
ours = FockState.basis((1, 0, 0)).apply(U).amps
b = perceval_amps(U, (1, 0, 0))
pc = {k: complex(b.prob_amplitude(pcvl.BasicState(list(k)))) for k in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]}
same = all(abs(ours.get(k, 0) - pc[k]) < 1e-9 for k in pc)
print("Perceval's Unitary uses our convention (column k: where a photon entering k goes):", same)
conv = (lambda M: M) if same else (lambda M: M.T)

# 2. random interferometers
rows, worst_routes, worst_pcvl = [], 0.0, 0.0
for m in range(4, 9):
    for n in range(1, 5):
        for c in range(3):
            U = random_unitary(m, rng)
            t = [0] * m
            for k in rng.choice(m, size=n, replace=True):
                t[k] += 1
            st = FockState.basis(tuple(t))
            a, p = st.apply(U), st.apply_perm(U)
            keys = set(a.amps) | set(p.amps)
            d_routes = max(abs(a.amps.get(k, 0) - p.amps.get(k, 0)) for k in keys)
            b = perceval_amps(conv(U), tuple(t))
            d_pcvl = max(abs(a.amps.get(k, 0) - complex(b.prob_amplitude(pcvl.BasicState(list(k))))) for k in keys)
            worst_routes, worst_pcvl = max(worst_routes, d_routes), max(worst_pcvl, d_pcvl)
            rows.append({"modes": m, "photons": n, "outputs": len(keys), "max_diff_routes": d_routes, "max_diff_perceval": d_pcvl})
print(f"{len(rows)} interferometers (4–8 modes, 1–4 photons): largest amplitude difference, our two routes {worst_routes:.1e}, ours vs Perceval {worst_pcvl:.1e}")

# 3. power check: one beamsplitter's sign flipped
M = beamsplitter(4, 0, 1, 0.3) @ beamsplitter(4, 1, 2, 0.7)
bad = beamsplitter(4, 0, 1, -0.3) @ beamsplitter(4, 1, 2, 0.7)
st = FockState.basis((1, 1, 0, 0))
b = perceval_amps(conv(M), (1, 1, 0, 0))
caught = max(abs(st.apply(bad).amps.get(k, 0) - complex(b.prob_amplitude(pcvl.BasicState(list(k))))) for k in st.apply(M).amps)
print(f"power check, one sign flipped: largest difference {caught:.2f} (must be far above {worst_pcvl:.0e})")
json.dump({"interferometers": len(rows), "max_diff_routes": worst_routes, "max_diff_perceval": worst_pcvl,
           "perceval_same_convention": same, "power_check_diff": caught, "perceval": pcvl.__version__, "rows": rows},
          open(Path(sys.argv[1]) / "fock-vs-perceval.json", "w"), indent=1)
