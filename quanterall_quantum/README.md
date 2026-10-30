# Quanterall Quantum Demo Suite

**10 live quantum circuit demonstrations running natively on the BEAM**  
Built with Erlang, Elixir and Nova 0.13. The simulator is `QxSim` in `lib/qx_demo.ex`: a pure-Elixir state-vector simulator on Erlang tuples, with no simulator library underneath.

Also here: `QxStab` (`lib/qx_stab.ex`), a stabilizer simulator — a tableau for the noiseless reference and a
bit-packed Pauli-frame sampler with gate-level noise, batches of shots in parallel — with repetition- and
surface-code memory circuits; `to_stim/1` writes any circuit in Stim's format. Endpoints: `/quantum/shor` (N = 15 or 21, a), `/quantum/grover/search`
(n, marked, iterations), `/quantum/stabilizer/ghz` (n ≤ 2000), `/quantum/stabilizer/memory` (code, d ≤ 7, rounds, p).
The published study and benchmarks: `tools/beam_runs.exs` (BEAM side) and `tools/stab_study.py` (Stim, PyMatching,
Qiskit Aer); `tools/dynamic_xval.py` (QxSim's measured circuits against Aer's dynamic circuits), `tools/stim_exact.py`
(exact tests against Stim), `tools/qx_compare.py` (the comparison with qx).

Origin and credit: the demos began on [qx](https://github.com/richarc/qx) by Craig Richards (Apache-2.0), version 0.5.1, which earlier versions of this README named as "qx_sim 0.5.1". `QxSim` is a separate implementation: compared on 29 October 2026 with every qx release from 0.2.0 to 0.11.0 (all built on Nx tensors), it shares no code with them. See CREDITS.md.

---

## What's Inside

| # | Endpoint | Demo | Type |
|---|---|---|---|
| 1 | `/quantum/bell` | Bell State Entanglement | Quantum-AI |
| 2 | `/quantum/teleport` | Quantum Teleportation | Quantum-AI |
| 3 | `/quantum/grover` | Grover's Search | Quantum-AI |
| 4 | `/quantum/qml` | Quantum Classifier | Quantum-AI |
| 5 | `/quantum/shor15` | Mini Shor's (N=15) | Quantum-AI |
| 6 | `/quantum/qec/bitflip` | 3-Qubit Bit-Flip Code | QEC |
| 7 | `/quantum/qec/phaseflip` | 3-Qubit Phase-Flip Code | QEC |
| 8 | `/quantum/qec/shor9` | Shor's 9-Qubit Code | QEC |
| 9 | `/quantum/qec/steane7` | Steane [[7,1,3]] Code | QEC |
| 10 | `/quantum/qec/surface3` | Distance-3 Surface Code | QEC |

---

## Project Structure

```
quanterall_quantum/
├── Dockerfile
├── docker-compose.yml
├── mix.exs
├── config/
│   ├── config.exs
│   └── runtime.exs
├── lib/
│   ├── application.ex
│   └── qx_demo.ex          ← all 10 quantum demos (Elixir)
├── src/
│   ├── quanterall_quantum_router.erl
│   └── controllers/
│       ├── quantum_controller.erl   ← 10 Nova route handlers
│       └── static_controller.erl
├── priv/static/
│   └── index.html          ← interactive demo UI
└── test_all_quantum_demos.sh
```

---

## Quick Start

### Docker (recommended)

```bash
docker compose up --build
```

Open `http://localhost:8080` in your browser.

### Manual

```bash
mix deps.get
mix compile
mix run --no-halt
```

---

## API

All endpoints accept `?shots=N` (default 1024, max 10000).  
QEC endpoints also accept `?error_qubit=N` (0-indexed, or `none`).

### Example requests

```bash
# Bell state
curl "http://localhost:8080/quantum/bell?shots=512"

# Bit-flip correction with error on qubit 1
curl "http://localhost:8080/quantum/qec/bitflip?shots=1024&error_qubit=1"

# Health check
curl "http://localhost:8080/health"
```

### Response schema (AI demos)

```json
{
  "status": "ok",
  "name": "bell",
  "qubits": 2,
  "shots": 1024,
  "execution_time_ms": 12,
  "counts": { "00": 512, "11": 512 },
  "visualisation_hint": "histogram",
  "plots": { "counts": { /* Vega-Lite spec */ } },
  "circuit_qasm": "OPENQASM 3.0; ...",
  "metadata": { "description": "...", "version": "2.5" }
}
```

### Response schema (QEC demos)

QEC responses additionally include:  
`counts_before`, `corrected_counts`, `syndrome_counts`, `logical_fidelity`, `error_type`

---

## Run Tests

```bash
chmod +x test_all_quantum_demos.sh
./test_all_quantum_demos.sh http://localhost:8080
```

---

## Notes

- **Shor-15**: QPE initialisation skeleton only. Full modular-exponentiation requires hardware beyond current simulators.
- **Steane-7 and Surface-3**: Simplified demonstration circuits. Full generator-matrix encoders omitted for clarity.
- All circuits run on the BEAM process scheduler — no external quantum hardware required.
- CORS headers are set to `*` for easy frontend integration.

---

## Architecture

```
Browser / curl
     │
     ▼
Nova (cowboy) HTTP server
     │
     ▼
quantum_controller.erl   (Erlang — route dispatch)
     │
     ▼
QxDemo.run_*/2           (Elixir — circuit construction + simulation)
     │
     ▼
QxSim                    (state-vector simulator in lib/qx_demo.ex, pure Elixir)
     │
     ▼
JSON response
```

The Erlang controller delegates to the Elixir module via the BEAM's transparent cross-language function calls (`'Elixir.QxDemo':run_bell/1`). No FFI, no processes crossing language boundaries — it's all one VM.
