# Optics conventions (QxFock, quantum/fock.py)

Every sign in the photonic simulations follows this file. Both engines implement it; Perceval (Quandela) checks it.

## States
- A state is a sparse map from occupation tuples (n₀, n₁, …, n_{m−1}) to complex amplitudes; modes numbered from 0.
- The occupation basis vector is |n⟩ = Π_k (a_k†)^{n_k} / √(n_k!) |vac⟩.
- **Dual rail.** A qubit on modes (2j, 2j+1): |0⟩_L = |1,0⟩, |1⟩_L = |0,1⟩ (Romero–Milburn Eqs. 28–29; KLM's
  "bosonic qubit"). Polarisation: mode 2j = H, 2j+1 = V.

## Passive linear optics
- An interferometer is a unitary m × m matrix M acting on creation operators: a_k† → Σ_j M_{jk} a_j†.
  (Column k is where a photon entering mode k goes.) Output amplitudes for input t and output s:
  ⟨s|U(M)|t⟩ = Perm(M[s,t]) / √(Π s_j! Π t_k!), rows of M repeated s_j times, columns t_k times.
- **Beamsplitter** B(i, j, θ): a_i† → cos θ a_i† − sin θ a_j†, a_j† → sin θ a_i† + cos θ a_j† — the matrix
  [[cos θ, sin θ], [−sin θ, cos θ]] on rows/columns (i, j) in the column convention above. On a dual-rail qubit
  (i, j) = (2q, 2q+1): |0⟩_L → cos θ |0⟩_L − sin θ |1⟩_L, |1⟩_L → sin θ |0⟩_L + cos θ |1⟩_L — the qubit
  matrix [[cos θ, sin θ], [−sin θ, cos θ]] = e^{iθY}.
  (KLM write B_θ = [[cos θ, −sin θ], [sin θ, cos θ]] for the mode operators themselves; the review's Eqs. 47–48
  are the same up to which side is called input. The NS gate below comes out right in both readings — checked.)
- **Phase shifter** P(k, φ): a_k† → e^{iφ} a_k†.
- **Balanced beamsplitter**: B(i, j, π/4). Hong–Ou–Mandel: |1,1⟩ → (|2,0⟩ − |0,2⟩)/√2, no coincidences.
- **Polarising beamsplitter** (PBS) on two spatial ports with polarisation modes (1H, 1V, 2H, 2V): H is transmitted,
  V reflected — the V modes swap ports: (1H, 1V, 2H, 2V) → (1H, 2V, 2H, 1V).
- **Half-wave plate at angle α** on (H, V): the rotation B(H, V, 2α) followed by P(V, π) (a reflection);
  a 22.5° plate maps H → D = (H + V)/√2. (Used for the fusion gates' 45° polarisation rotation.)

## Measurement
- **Number-resolving detector** on mode k: outcome n with the Kraus projector |n⟩⟨n|_k; the remaining modes keep
  the unnormalised conditional state; its squared norm is the outcome's probability.
- **Bucket detector**: outcome 0 (n = 0) or 1 (n ≥ 1).
- **Detector efficiency η**: a beamsplitter of transmissivity η onto a fresh vacuum mode, which is traced out;
  POVMs E₀ = Σ (1 − η)ⁿ |n⟩⟨n|, E₁ = 1 − E₀ (bucket); E_n = Σ_{k≥n} C(k, n) ηⁿ (1 − η)^{k−n} |k⟩⟨k| (resolving).

## The KLM gates (from Knill, Laflamme, Milburn, Nature 409, 46 (2001), Eq. 9)
- **NS gate** on mode 1 with ancilla modes 2, 3 prepared |1, 0⟩: the 3 × 3 symmetric unitary
  U = [[1 − √2, 2^{−1/4}, (3/√2 − 2)^{1/2}], [2^{−1/4}, 1/2, 1/2 − 1/√2], [(3/√2 − 2)^{1/2}, 1/2 − 1/√2, √2 − 1/2]];
  accept when modes 2, 3 read |1, 0⟩. Heralded amplitudes for 0, 1, 2 photons in mode 1: ½, ½, −½ — NS with
  probability 1/4. The review's beamsplitter form reproduces it: B(2,3, 22.5°) B(1,2, 65.5302°) B(2,3, −22.5°) with
  a π phase on mode 1 (θ₃ = −22.5° as in the review's text; its Fig. 8 caption's +22.5° does not — checked).
- **CZ** on dual-rail qubits (modes 1,2 and 3,4, the |1⟩_L modes 1 and 3 meeting): balanced beamsplitter on the two
  |1⟩_L modes, NS on each, the inverse beamsplitter; success 1/16.
