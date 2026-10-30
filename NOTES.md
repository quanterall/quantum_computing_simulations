# Notes on Romero & Milburn, "Photonic Quantum Computing", arXiv:2404.03367v2

Equation issues found while implementing (docs/photonic-plan.md); each number we use comes from the primary
paper or our own derivation, never from the review alone.

| Where | Printed | Correct | How we know |
|---|---|---|---|
| Eq. 13 | P₂(t) = 1 − η\|ξ(t)\|² | (1 − η)\|ξ(t)\|² | the rate must be proportional to the pulse; P₁ + P₂ = \|ξ\|² |
| Eq. 22 | S(ζ)\|0⟩ = e^{(ζ*a² − ζa†²)/2} | … e^{…}\|0⟩ | the operator needs the vacuum to act on |
| Eq. 40 | N₊^{−1}, N± = 2 ± 2e^{−2\|α\|²} | N₊^{−1/2} (as Eq. 75 has) | normalisation of \|α⟩ ± \|−α⟩ |
| p. 11 | overlap ∼ 10⁻⁷ "for α ≤ 2" | α ≥ 2: e^{−4α²} = e^{−16} ≈ 1.1 × 10⁻⁷ at α = 2 | arithmetic |
| Eq. 56 vs Fig. 10 | Eq. 56 couples 12, 23, 24 (a star around qubit 2); Fig. 10's caption a chain CZ₁₂CZ₂₃CZ₃₄ | **the star is right for the CNOT claim**: with CZ gates it gives CNOT (qubit 4 control, qubit 3 = a ⊕ b) for every outcome and superposition input, with Pauli byproducts X^{s₂}Z^{s₁} on 3 and Z^{s₁} on 4; the chain gives a CZ between the teleported input and qubit 4 — CNOT-class (Makhlin invariants 0, 1) but needing Hadamards before as well as after. (An earlier note here, from an external comment, guessed σz⁽³⁾σz⁽⁴⁾: wrong.) | tools/optics/mbqc.py |
| Fig. 8 | caption θ₃ = 22.5°, θ₂ = 65.5° | θ₃ = −22.5° (the text), θ₂ = 65.5302° | only −22.5° reproduces KLM's Eq. 9 unitary (numerically, tools/optics/ns_check.py) |
| Eq. 70 | \|Υ\|² = (2π)^{−1/2} exp[−(y − 2r cos(θq − ϕ))²] | exp[−(…)²/2] | the square of Eq. 69's exp[−(…)²/4] |
| p. 19 | "does not increase … on average", then "+½ qubit", "2 × 4 − 1 = 7 … 6.5 per added qubit" | Bell pairs alone: zero drift; 3-qubit pieces: 8 pairs per qubit, 7 with the leftover reused (the "2 × 4 − 1"); 6.5 not found in any strategy simulated (best 7.0) | tools/optics/cluster_growth.py, Monte Carlo |
| p. 27 | 2.7 % single-photon loss tolerated | Bartolucci et al., Nat. Commun. 14, 912 (2023), Sec. VI.2 + App. B: 6-ring network, (2,2)-Shor-encoded resource states, fusion boosted with a Bell pair (p_fail 1/4, 4 photons), loss on every photon; their simulated marginal erasure threshold 11.98 % through Eqs. B1–B2 gives 2.708 % per photon = 10.4 % per fusion (and 43.2 % failure threshold at zero loss) — **confirmed**, the threshold itself taken from the paper | tools/optics/fbqc_loss.py |
| Eqs. 57–58 | byproduct U_con, with σ_y rotations | with CZ gates the byproducts are Paulis (above); the review's non-Pauli U_con comes from writing the entangler as exp(−iπ/4·ZZ), which is CZ only up to local Z phases — not checked beyond that | tools/optics/mbqc.py |
| Eqs. 22, 25 | ζ = r e^{−2iϕ} with Y1 + iY2 = (X1 + iX2)e^{−iϕ} | one sign flipped: with Eq. 22 as printed the rotation is e^{+iϕ} (or ζ = r e^{+2iϕ} with Eq. 25) | as printed, Eq. 26's variances fail for ϕ ≠ 0 (off by up to 5.6); either flip gives e^{∓2r} to 10⁻¹⁴ — tools/optics/cv_checks.py; Eq. 24 is right for Eq. 22 as printed |
| Fig. 5 caption | variance of ψ0 "e^{2r}", r = 1.2 | e^{−2r} | read from the PDF's vector curve: peak variance fits 0.11 (e^{−2r} = 0.09; e^{2r} = 11) |
| Fig. 5 panels | one Δ in Eq. 38 for both states | the \|1̃⟩ panel's envelope falls off 4.00× faster than \|0̃⟩'s; with Δ̃ = 0.01, amplitudes e^{−Δ̃(sα)²} (\|0̃⟩) and e^{−Δ̃[(2s+1)α]²} (\|1̃⟩) fit: s where Eq. 38 has 2s | peaks at multiples of √π (α = √π); tools/optics/cv_checks.py, site/data/gkp-fig5-digitised.json |
| Eq. 40 | N±^{−1} | confirmed numerically: norm² 1/N± (0.49983 at α = 2); N±^{−1/2} as Eq. 75 | tools/optics/cv_checks.py |
| Eq. 70 | as printed | confirmed numerically: ∫\|Υ\|² dy = 1/√2 as printed, 1 for the square of Eq. 69 | tools/optics/cv_checks.py |
| Eqs. 80–81 | detector POVMs | confirmed: a lossy beamsplitter onto vacuum plus an ideal detector gives them to 4 × 10⁻¹⁶ | tools/optics/loss_checks.py |
