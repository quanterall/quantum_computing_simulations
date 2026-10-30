# Photonic quantum computing on our simulators — implementation plan

*Source read: J. Romero and G. Milburn, "Photonic Quantum Computing", arXiv:2404.03367v2 (9 Oct 2025), a 29-page
review (encyclopedia article). Plan written 29 October 2026, after reading the paper and two external comments.*

## 0. What this is, and what it is not

The paper is a **review**: it reports results from primary papers and adds no new results of its own. So we do not
"implement the paper". We **reproduce, by simulation, the schemes it reviews**, from their primary sources:

| Scheme | Primary source |
|---|---|
| Linear-optics gates, NS gate, CZ from two NS gates | Knill, Laflamme, Milburn, *Nature* 409, 46 (2001) |
| Type-I and type-II fusion, cluster growth | Browne, Rudolph, *PRL* 95, 010501 (2005) |
| Fusion-based architecture, the 2.7 % loss figure | Bartolucci et al., *Nat. Commun.* 14, 912 (2023) |
| GKP encoding | Gottesman, Kitaev, Preskill, *PRA* 64, 012310 (2001) |
| Cat-state code | Ralph et al., *PRA* 68, 042319 (2003) |
| HOM interference | Hong, Ou, Mandel, *PRL* 59, 2044 (1987) |
| Detector POVMs | Kok, Lovett, *Introduction to Optical Quantum Information Processing* (2010) |

Wording on the site: **"simulation of schemes from the literature"**, under the existing *Reproduction* heading —
never photonic hardware results. No text or figures copied from the review; equations re-derived; primary papers cited.

The one part that **extends our own research** rather than reproducing: photon loss is a *heralded erasure*, which
union-find decoders handle natively (Phase 6). It connects directly to QxDecode and the load-test question.

## 1. The review's equations: our own check

Checked while reading the paper's text (not yet against the primary papers — Phase 0 confirms each there):

| Where | What the review prints | What it should be | Status |
|---|---|---|---|
| Eq. 13 | P₂(t) = 1 − η\|ξ(t)\|² | (1 − η)\|ξ(t)\|² (a rate proportional to the pulse) | confirmed by context |
| Eq. 22 | S(ζ)\|0⟩ = e^{(ζ*a² − ζa†²)/2} | … e^{…}\|0⟩ — the vacuum is missing | confirmed |
| Eq. 40 | \|0⟩ = N₊^{−1}(\|α⟩+\|−α⟩), N± = 2 ± 2e^{−2\|α\|²} | N₊^{−1/2}; Eq. 75 writes it correctly | confirmed (Eq. 75 agrees with the fix) |
| p. 11 | overlap e^{−4α²} "for α ≤ 2 … ∼ 10⁻⁷" | α ≥ 2 (e^{−16} ≈ 1.1 × 10⁻⁷) | confirmed |
| Eq. 56 | σz⁽¹⁾σz⁽²⁾ + σz⁽²⁾σz⁽³⁾ + σz⁽²⁾σz⁽⁴⁾ | the chain CZ₁₂CZ₂₃CZ₃₄ (Fig. 10) needs σz⁽³⁾σz⁽⁴⁾; and e^{−iπ/4·ZZ} is CZ only up to local Z rotations | confirmed inconsistency |
| Eqs. 57–58 | the byproduct U_con | determine numerically against CNOT, for superposition inputs too (Eq. 55 uses basis inputs) | to derive |
| Fig. 8 vs text | text θ₃ = −22.5°, caption θ₃ = 22.5°; caption θ₂ = 65.5° vs text 65.5302° | settle the sign numerically (success 1/4 and the −α₂ sign) and from KLM 2001 | **new, not in the comments** |
| Eq. 70 | \|Υ(q\|y)\|² = (2π)^{−1/2} exp[−(y − 2r cos(θq − ϕ))²] | exp[−(…)²/2], since Eq. 69 has exp[−(…)²/4] in the amplitude | **new, not in the comments** |
| p. 19, fusion | "does not increase the size on average" and then "+½ qubit", "2 × 4 − 1 = 7 … thus 6.5 Bell pairs per added qubit" | contradictory; settle by simulating the growth strategy (Phase 3) | confirmed contradiction |
| p. 27 | "can tolerate a 2.7 % probability of suffering a single-photon loss" | encoded 6-ring network, Bell-pair-boosted fusion, Eqs. B1–B2 of Bartolucci et al. | **reproduced** (tools/optics/fbqc_loss.py) |

A `NOTES.md` in the repository keeps this table, with the outcome of each check.

## 2. Where each piece runs

| Component | Our simulator (quanterall.com, Elixir) | TechUPGRADE engine (sim-manager, Python) | Independent reference |
|---|---|---|---|
| Linear optics in Fock space | **QxFock** (`lib/qx_fock.ex`, new) | **quantum/fock.py** (new, NumPy) — a second implementation | **Perceval** (Quandela; installable here) and closed forms |
| Dual-rail single-qubit gates | QxFock vs QxSim | fock.py vs engine.py | Qiskit |
| Measurement-based computation | QxSim (cluster primitives) | engine.py (same) | Qiskit dynamic circuits |
| Cluster growth by fusion | QxStab (graph-state rules) | — | brute-force QxSim on small clusters |
| Loss as erasure, circuit level | QxStab (`:heralded_erase`) | — | **Stim 1.16 `HERALDED_ERASE`** (present) |
| Erasure decoding | QxDecode (erased edges pre-grown) | — | PyMatching with per-shot reweighting |

Two implementations of the Fock core — one in each engine, written from the same conventions file — plus Perceval
give three-way agreement, the practice we already use (QxSim vs NumPy vs Qiskit; QxStab vs Stim).

## 3. Phases

Sizes for one engineer: S ≈ 1–2 days, M ≈ 3–5 days, L ≈ 1–2 weeks. Every phase ends with: an independent
cross-check, a **power check** (a deliberately broken version must fail), seeds, published data and scripts.

### Phase 0 — conventions and sources (S)
- `docs/optics_conventions.md`: dual-rail \|0⟩_L = \|1,0⟩, \|1⟩_L = \|0,1⟩; beamsplitter â₁ → cos θ â₁ + sin θ â₂,
  â₂ → cos θ â₂ − sin θ â₁ (Eqs. 47–48); phase gate (Eq. 51); mode order; outcome labels; polarisation H/V ↔ modes.
- Read the primary papers for every number we will check; record each in the conventions file with its source.
- Install Perceval in `.venv-tools`; one test reproducing a beamsplitter matrix in both conventions.

*Done when* a reviewer can reproduce any sign from the file alone.

### Phase 1 — Fock-space core in both engines (M)
- State: a sparse map from occupation tuples to complex amplitudes; ≤ 12 modes, ≤ 6 photons (8 modes and
  4 photons are 330 basis states — the CZ gate).
- Gates: beamsplitter, phase shifter, polarising beamsplitter, waveplate; each by transforming creation operators,
  and a second path through permanents (Ryser) for whole interferometers; the two must agree.
- Measurement: number-resolving and bucket detectors on a subset of modes, returning outcome, probability and the
  heralded (unnormalised) state of the rest.
- Checks: norm preserved; dual-rail single-photon gates equal QxSim's / engine.py's single-qubit gates; \|1,1⟩ on a
  50/50 beamsplitter gives no coincidences; random 4–8-mode interferometers with 1–4 photons, every output
  probability equal in QxFock, fock.py and Perceval. Power check: one beamsplitter sign flipped must be caught.
- **HOM dip** with partial distinguishability (an orthogonal internal mode per photon): coincidence probability
  ½(1 − e^{−γ\|τ\|}) for the exponential pulse of Eq. 43 (Eq. 45) — published as a curve, data and endpoint.

### Phase 2 — KLM gates (M)
- **NS gate**: modes 1–3, ancilla \|1⟩₂\|0⟩₃, angles from Phase 0 (θ₃'s sign settled), phase π on mode 1. Checks: for
  random α, heralding (1, 0) gives α₀\|0⟩ + α₁\|1⟩ − α₂\|2⟩ with probability 1/4, independent of α.
- **CZ from two NS gates** (8 modes, 4 photons): the heralded map on the dual-rail two-qubit space equals CZ with
  process fidelity 1; success 1/16.
- Power checks: one angle off by 0.5° (fidelity must drop; report how far); no π phase (CZ must become the identity).

### Phase 3 — fusion and cluster growth (M–L)
- **Type-I fusion** in the Fock core (PBS, 45° rotation, one number-resolving detector): success on exactly one
  click, probability 1/2; the Kraus operators \|0⟩⟨00\| ± \|1⟩⟨11\|; failure equal to Z measurements of both qubits.
- **Type-II fusion** (two detectors): XX and ZZ parity projection on success, probability 1/2; the failure effect.
- **Graph rules** (success merges two end vertices: n + m − 1; failure removes both) checked against QxStab on
  small clusters for every outcome sequence.
- **Growth Monte Carlo** of the Browne–Rudolph strategy for linear and 2D clusters: Bell pairs per added qubit
  with error bars, and the strategy stated exactly — this settles the review's arithmetic on p. 19.

### Phase 4 — measurement-based computation (S–M)
- QxSim / engine.py primitives: cluster states, rotated-basis measurement, byproduct tracking.
- Checks with random inputs: one-bit teleportation; an arbitrary single-qubit rotation along a chain; the
  four-node cluster of Fig. 10 — find numerically what it implements for superposition inputs, and the byproduct
  (Eqs. 57–58) — reported as found, including if a CNOT holds only for basis inputs. Cross-check against Qiskit
  dynamic circuits. Power check: a wrong byproduct fails.

### Phase 5 — loss and detectors (M)
- Detector POVMs (Eqs. 80–81) checked against a loss beamsplitter with vacuum environment modes; E₀ + E₁ = 1.
- Heralded-gate fidelity against detector efficiency η for NS and CZ: false heralds from loss; a published curve.
- Loss on a dual-rail qubit is exactly ρ → ηρ + (1 − η)\|vac⟩⟨vac\| — checked; this justifies the erasure model.

### Phase 6 — erasure decoding (L) — the research-relevant phase
- **QxStab**: `:heralded_erase` (the erased qubit fully depolarised, the herald recorded), matching Stim's
  `HERALDED_ERASE`; the sampler checked against Stim on the same circuits (exact tests, as for Pauli noise).
- **QxDecode**: erased edges start fully grown; peeling on the erased subgraph is the linear-time maximum-likelihood
  decoder for pure erasure (Delfosse–Zémor); mixed erasure + Pauli noise through the existing growth.
- Experiments: (a) code capacity, pure erasure — the surface-code curves must cross near 50 % (bond percolation on
  the square lattice): the validation target; (b) mixed erasure and Pauli noise, QxDecode vs PyMatching with
  per-shot weights (erased edges at weight 0); (c) circuit level, gate noise plus erasure, the threshold as the
  erased fraction grows; (d) the 2.7 % figure, once Phase 0 has established what it refers to.
- Tie to the research question: erasure syndromes in the load test (`tools/load/`), with batches per process.

### Phase 7 — optional: GKP and cat states (M)
**Done 29 Oct 2026** (user: "go with phase 7"). qsim/qsim/cv.py (= TechUPGRADE quantum/cv.py: Fock basis + q-grid) and
lib/qx_cv.ex (q-grid, p by analytic Fourier sum); tools/optics/cv_checks.py (23 checks), xcheck_cv.py (BEAM vs Python),
cv_charts.py; /quantum/optics/{gkp,cat,comb} on the BEAM, /quantum/photonic/{gkp,cat,comb} on TechUPGRADE. Added beyond the
plan: cat parity under loss (closed form), the heralded comb (Eqs. 66–73). New errata: Eqs. 22/25 sign, Fig. 5 caption and
envelopes. GBS (Eq. 74) still left out.
- One mode on a truncated Fock basis or a q-grid: squeezed-state variances e^{±2r}; the Eq. 38 comb; the GKP
  logical error rate under Gaussian shift noise against the closed form; cat-state normalisation and overlap (the
  corrected forms). Gaussian boson sampling is left out (exponential cost, weak case next to the rest).

### Phase 8 — publication, TechUPGRADE, review (M, overlapping)
- **quanterall.com**: a "Reproduction · photonic schemes" section — HOM curve, NS/CZ numbers, fusion statistics,
  the erasure threshold; endpoints `/quantum/optics/{hom,ns,cz,fusion,loss}` with seeds and caps; the
  "which code runs what" table extended; data and charts; the errata table as `NOTES.md` in the repository.
- **TechUPGRADE (quanterall tenant)**: sim-manager `quantum/fock.py` with endpoints `/quantum/photonic/*`; a
  *Photonics* tab in the Quantum panel (HOM dip with a delay slider, NS/CZ heralding statistics, fusion outcomes,
  detector-efficiency curve, erasure threshold chart); tests in `tests/test_quantum.py`; a smoke step and contract
  cases (quanterall only, 403 for other tenants); released by `plants/quanterall/go-live-beam.sh`.
- **Review**: an optics physicist (IAPS may know one) for Phases 1–3 — outside our team's own physics.
- **Release**: tag, Zenodo DOI; a short preprint on Phase 6 (erasure decoding with union-find on the BEAM), the
  part that is ours.

## 4. Order, and what comes first

1. Phase 0 + Phase 1 in both engines (with the HOM curve) — exact, quick, gives the cross-check machinery.
2. Phase 2 (NS, CZ).
3. Phase 6 in parallel with Phase 3 — Phase 6 does not depend on the Fock core, and it is the research part.
4. Phases 3, 4, 5; then Phase 8's TechUPGRADE tab once 1–3 are stable; Phase 7 only if time allows.

Rough total: Phases 0–6 and 8 ≈ 6–9 weeks for one engineer; Phase 7 adds about a week.

## 5. Risks

- **Conventions** (signs, port labels) are the commonest failure — the conventions file and three-way agreement
  exist for this.
- **Fock-space size**: stay within 12 modes and 6 photons in QxFock; qubit-level work goes to QxSim / QxStab.
- **Claims drift**: resource counts such as "300 CZ₉/₁₆ gates for 95 %" or "52 Bell states" belong to the primary
  sources; we reproduce a number only by simulating it, otherwise we cite it.
- **The review's slips**: every number is taken from the primary paper or derived, never copied from the review.

## 6. Decisions for Quanterall

1. Order: start with Phases 0–2, or with Phase 6 (erasure, the research part)?
2. Whether to ask IAPS for an optics reviewer now (Phases 1–3).
3. Phase 7 (GKP, cats): in or out.
4. Whether the TechUPGRADE *Photonics* tab ships with Phase 2, or waits for Phase 6.
