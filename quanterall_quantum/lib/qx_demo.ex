defmodule QxSim do
  @moduledoc """
  Quantum state-vector simulator. Pure Elixir, zero external dependencies.

  Fix 1 (performance): State vector stored as an Erlang tuple — O(1) random
    access via elem/2 and put_elem/3, replacing O(n) list traversal.

  Fix 3 (correctness): Every gate reads ALL needed amplitudes from the
    ORIGINAL tuple first, then writes all results back atomically.
    No gate ever reads from a partially-modified state.

  Sampling uses binary search O(shots·log·2ⁿ) over a cumulative
    probability tuple, replacing O(shots·2ⁿ) linear scan.
  """
  import Bitwise

  # ── Types ────────────────────────────────────────────────────────────────────
  # State: an Erlang tuple of {2^n} elements, each element = {re::float, im::float}

  # ── Initialisation ───────────────────────────────────────────────────────────

  def init(n) do
    dim = 1 <<< n
    # Build as list then convert once — faster than repeated put_elem on large tuples
    List.to_tuple(
      [{1.0, 0.0} | List.duplicate({0.0, 0.0}, dim - 1)]
    )
  end

  # ── Single-qubit gate (Fix 3: read first, write after) ───────────────────────
  # mat = {{a_re,a_im}, {b_re,b_im}, {c_re,c_im}, {d_re,d_im}}
  # corresponds to 2×2 complex matrix [[a,b],[c,d]]

  # Every gate builds the new state in one pass over the old one (O(2ⁿ)): put_elem copies the whole
  # tuple, so writing amplitudes one by one cost O(4ⁿ) per gate.
  defp sg(state, q, {ar, ai}, {br, bi}, {cr, ci}, {dr, di}) do
    dim  = tuple_size(state)
    mask = 1 <<< q

    for i <- 0..(dim - 1) do
      if (i &&& mask) == 0 do
        {re_i, im_i} = elem(state, i)
        {re_j, im_j} = elem(state, i ||| mask)
        {ar * re_i - ai * im_i + br * re_j - bi * im_j, ar * im_i + ai * re_i + br * im_j + bi * re_j}
      else
        {re_i, im_i} = elem(state, Bitwise.bxor(i, mask))
        {re_j, im_j} = elem(state, i)
        {cr * re_i - ci * im_i + dr * re_j - di * im_j, cr * im_i + ci * re_i + dr * im_j + di * re_j}
      end
    end
    |> List.to_tuple()
  end


  defp cx_gate(state, ctrl, tgt) do
    dim    = tuple_size(state)
    m_ctrl = 1 <<< ctrl
    m_tgt  = 1 <<< tgt
    for i <- 0..(dim - 1) do
      if (i &&& m_ctrl) != 0, do: elem(state, Bitwise.bxor(i, m_tgt)), else: elem(state, i)
    end
    |> List.to_tuple()
  end


  defp cz_gate(state, ctrl, tgt) do
    dim  = tuple_size(state)
    both = (1 <<< ctrl) ||| (1 <<< tgt)
    for i <- 0..(dim - 1) do
      {re, im} = elem(state, i)
      if (i &&& both) == both, do: {-re, -im}, else: {re, im}
    end
    |> List.to_tuple()
  end


  @doc """
  Measures qubit q in the computational basis: returns {bit, state}, the state collapsed onto the outcome
  and renormalised.
  """
  def measure(state, _n, q) do
    dim  = tuple_size(state)
    mask = 1 <<< q
    p1 = Enum.reduce(0..(dim - 1), 0.0, fn i, acc ->
      if (i &&& mask) != 0 do
        {re, im} = elem(state, i)
        acc + re * re + im * im
      else
        acc
      end
    end)
    bit = if :rand.uniform() < p1, do: 1, else: 0
    norm = :math.sqrt(if bit == 1, do: p1, else: 1.0 - p1)
    collapsed =
      for i <- 0..(dim - 1) do
        if ((i &&& mask) != 0) == (bit == 1) do
          {re, im} = elem(state, i)
          {re / norm, im / norm}
        else
          {0.0, 0.0}
        end
      end
      |> List.to_tuple()
    {bit, collapsed}
  end

  @doc "Measures qubit q and puts it back to |0⟩, so it can be used again (an ancilla)."
  def measure_reset(state, n, q) do
    {bit, s} = measure(state, n, q)
    {bit, if(bit == 1, do: x(s, n, q), else: s)}
  end

  @doc "Samples shots of the first k qubits (qubit 0 rightmost), ignoring the others."
  def sample_low(state, k, shots) do
    cdf   = build_cdf(state)
    size  = tuple_size(cdf)
    total = elem(cdf, size - 1)
    mask  = (1 <<< k) - 1
    Enum.reduce(1..shots, %{}, fn _, acc ->
      idx = bisect(cdf, :rand.uniform() * total, 0, size - 1) &&& mask
      Map.update(acc, Integer.to_string(idx, 2) |> String.pad_leading(k, "0"), 1, &(&1 + 1))
    end)
  end

  defp inv2, do: 1.0 / :math.sqrt(2.0)

  def h(s, _n, q) do
    f = inv2()
    sg(s, q, {f, 0.0}, {f, 0.0}, {f, 0.0}, {-f, 0.0})
  end

  def x(s, _n, q), do: sg(s, q, {0.0,0.0}, {1.0,0.0}, {1.0,0.0}, {0.0,0.0})
  def y(s, _n, q), do: sg(s, q, {0.0,0.0}, {0.0,-1.0}, {0.0,1.0}, {0.0,0.0})
  def z(s, _n, q), do: sg(s, q, {1.0,0.0}, {0.0,0.0}, {0.0,0.0}, {-1.0,0.0})

  def rx(s, _n, q, t) do
    c = :math.cos(t / 2.0)
    sn = :math.sin(t / 2.0)
    sg(s, q, {c, 0.0}, {0.0, -sn}, {0.0, -sn}, {c, 0.0})
  end

  def ry(s, _n, q, t) do
    c = :math.cos(t / 2.0)
    sn = :math.sin(t / 2.0)
    sg(s, q, {c, 0.0}, {-sn, 0.0}, {sn, 0.0}, {c, 0.0})
  end

  def rz(s, _n, q, t) do
    c = :math.cos(t / 2.0)
    sn = :math.sin(t / 2.0)
    sg(s, q, {c, -sn}, {0.0, 0.0}, {0.0, 0.0}, {c, sn})
  end

  def cx(s, _n, ctrl, tgt), do: cx_gate(s, ctrl, tgt)

  @doc "Exchanges qubits a and b."
  def swap(s, _n, a, b), do: involution(s, fn i -> if bit(i, a) != bit(i, b), do: Bitwise.bxor(i, (1 <<< a) ||| (1 <<< b)), else: i end)

  @doc "Fredkin: exchanges a and b where the control is 1."
  def cswap(s, _n, c, a, b),
    do: involution(s, fn i -> if bit(i, c) == 1 and bit(i, a) != bit(i, b), do: Bitwise.bxor(i, (1 <<< a) ||| (1 <<< b)), else: i end)

  @doc "Controlled phase: e^{iθ} where both qubits are 1."
  def cp(s, _n, c, t, theta) do
    both = (1 <<< c) ||| (1 <<< t)
    {pc, ps} = {:math.cos(theta), :math.sin(theta)}
    for i <- 0..(tuple_size(s) - 1) do
      {re, im} = elem(s, i)
      if (i &&& both) == both, do: {re * pc - im * ps, re * ps + im * pc}, else: {re, im}
    end
    |> List.to_tuple()
  end

  @doc """
  Controlled multiplication by c mod nn on the work register `work` (its first qubit the lowest bit): where
  the control is 1, a work value y < nn becomes c·y mod nn; values ≥ nn are left alone. c must be coprime to
  nn, so this is a permutation of basis states — the modular-exponentiation oracle of Shor's algorithm,
  applied as one unitary rather than decomposed into elementary gates.
  """
  def cmulmod(s, _n, ctrl, work, c, nn) do
    cinv = Enum.find(1..(nn - 1), &(rem(&1 * c, nn) == 1))
    shift = hd(work)
    w = length(work)
    wmask = ((1 <<< w) - 1) <<< shift
    for i <- 0..(tuple_size(s) - 1) do
      y = (i &&& wmask) >>> shift
      if bit(i, ctrl) == 1 and y < nn do
        src = (i &&& bnot(wmask)) ||| (rem(cinv * y, nn) <<< shift)          # the amplitude that lands here
        elem(s, src)
      else
        elem(s, i)
      end
    end
    |> List.to_tuple()
  end

  @doc "Multi-controlled Z on the qubits in `qs`: a phase of −1 where all of them are 1."
  def mcz(s, _n, qs) do
    mask = Enum.reduce(qs, 0, fn q, m -> m ||| (1 <<< q) end)
    for i <- 0..(tuple_size(s) - 1) do
      {re, im} = elem(s, i)
      if (i &&& mask) == mask, do: {-re, -im}, else: {re, im}
    end
    |> List.to_tuple()
  end

  @doc "The probability of basis state k (qubit 0 the lowest bit)."
  def prob(s, k), do: (fn {re, im} -> re * re + im * im end).(elem(s, k))

  defp bit(i, q), do: (i >>> q) &&& 1

  # a gate that permutes basis states and is its own inverse: the new amplitude at i is the old one at f(i)
  defp involution(s, f), do: (for i <- 0..(tuple_size(s) - 1), do: elem(s, f.(i))) |> List.to_tuple()
  def cz(s, _n, ctrl, tgt), do: cz_gate(s, ctrl, tgt)

  # ── Sampling (Fix 1: O(shots·log·2ⁿ) binary search) ─────────────────────────

  defp build_cdf(state) do
    dim  = tuple_size(state)
    probs = for i <- 0..(dim - 1) do
      {re, im} = elem(state, i)
      re * re + im * im
    end
    # Build cumulative distribution as a tuple for O(1) index access
    {cum, _} = Enum.map_reduce(probs, 0.0, fn p, acc -> {acc + p, acc + p} end)
    List.to_tuple(cum)
  end

  defp bisect(cdf, r, lo, hi) when lo >= hi, do: lo
  defp bisect(cdf, r, lo, hi) do
    mid = div(lo + hi, 2)
    if elem(cdf, mid) < r,
      do:   bisect(cdf, r, mid + 1, hi),
      else: bisect(cdf, r, lo, mid)
  end

  defp sample(state, n, shots) do
    cdf  = build_cdf(state)
    size = tuple_size(cdf)
    total = elem(cdf, size - 1)

    Enum.reduce(1..shots, %{}, fn _, acc ->
      r   = :rand.uniform() * total
      idx = bisect(cdf, r, 0, size - 1)
      k   = Integer.to_string(idx, 2) |> String.pad_leading(n, "0")
      Map.update(acc, k, 1, &(&1 + 1))
    end)
  end

  # ── Public API ────────────────────────────────────────────────────────────────

  def run(n, shots, circuit_fn) do
    t1 = System.monotonic_time(:microsecond)
    counts = n |> init() |> circuit_fn.() |> sample(n, shots)
    {counts, Float.round((System.monotonic_time(:microsecond) - t1) / 1000, 3)}
  end
end

defmodule QxDemo do
  @moduledoc """
  Quanterall's quantum demos: algorithms and protocols (Bell, teleportation, Grover on 2 to 10 qubits, a
  variational classifier, Shor for N = 15) and five error-correcting codes with full correction cycles.
  Uses QxSim (tuple-based state-vector simulator); no external quantum dependencies.
  """
  import Bitwise
  defp bit(i, q), do: (i >>> q) &&& 1

  @default_shots 1024

  # ============================================================
  # QUANTUM-AI DEMOS
  # ============================================================

  def run_bell(shots \\ @default_shots) do
    n = 2
    {counts, ms} = QxSim.run(n, shots, fn s ->
      s |> QxSim.h(n, 0) |> QxSim.cx(n, 0, 1)
    end)
    resp("bell", n, shots, ms, counts,
      "Bell state |Φ⁺⟩ = (|00⟩+|11⟩)/√2 — maximally entangled pair",
      "OPENQASM 3.0;\nqubit[2] q; bit[2] c;\nh q[0];\ncx q[0],q[1];\nc[0]=measure q[0]; c[1]=measure q[1];\n")
  end

  def run_teleport(shots \\ @default_shots) do
    # A random state ψ = Rz(φ)·Ry(θ)|0⟩ on Alice's q0 (θ, φ from the seed); a Bell pair q1–q2; Alice's Bell
    # measurement; Bob's corrections X (from q1) and Z (from q0) as controlled gates (deferred measurement).
    # Bob then undoes the preparation, Ry(−θ)·Rz(−φ), on q2: if ψ arrived intact, q2 reads 0 in every shot. A
    # bit copy would not pass this — ψ has a phase and a superposition, not just a value.
    n = 3
    theta = :math.acos(1 - 2 * :rand.uniform())
    phi = 2 * :math.pi() * :rand.uniform()
    {counts, ms} = QxSim.run(n, shots, fn s ->
      s
      |> QxSim.ry(n, 0, theta) |> QxSim.rz(n, 0, phi)
      |> QxSim.h(n, 1) |> QxSim.cx(n, 1, 2)
      |> QxSim.cx(n, 0, 1) |> QxSim.h(n, 0)
      |> QxSim.cx(n, 1, 2) |> QxSim.cz(n, 0, 2)
      |> QxSim.rz(n, 2, -phi) |> QxSim.ry(n, 2, -theta)
    end)
    {td, pd} = {Float.round(theta * 180 / :math.pi(), 2), Float.round(phi * 180 / :math.pi(), 2)}
    ok = counts |> Enum.filter(fn {k, _} -> String.first(k) == "0" end) |> Enum.map(&elem(&1, 1)) |> Enum.sum()
    resp("teleport", n, shots, ms, counts,
      "Quantum teleportation of a random state (θ = #{td}°, φ = #{pd}°): Bob undoes the preparation, so q2 reads 0 if ψ arrived",
      "OPENQASM 3.0;\nqubit[3] q; bit[3] c;\nry(#{Float.round(theta, 6)}) q[0]; rz(#{Float.round(phi, 6)}) q[0];   // ψ, random\nh q[1]; cx q[1],q[2];        // the Bell pair\ncx q[0],q[1]; h q[0];        // Alice's Bell measurement\ncx q[1],q[2]; cz q[0],q[2];  // Bob's corrections X and Z (deferred measurement)\nrz(#{Float.round(-phi, 6)}) q[2]; ry(#{Float.round(-theta, 6)}) q[2];   // Bob undoes the preparation\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2];   // q[2] = 0 every time\n")
    |> Map.merge(%{state: %{theta_deg: td, phi_deg: pd}, bob_zero_fraction: Float.round(ok / shots, 4)})
  end

  def run_grover(shots \\ @default_shots) do
    n = 2
    {counts, ms} = QxSim.run(n, shots, fn s ->
      s
      |> QxSim.h(n, 0) |> QxSim.h(n, 1)
      |> QxSim.cz(n, 0, 1)
      |> QxSim.h(n, 0) |> QxSim.h(n, 1)
      |> QxSim.z(n, 0) |> QxSim.z(n, 1)
      |> QxSim.cz(n, 0, 1)
      |> QxSim.h(n, 0) |> QxSim.h(n, 1)
    end)
    resp("grover", n, shots, ms, counts,
      "Grover's search — oracle marks |11⟩, amplitude amplification boosts it",
      "OPENQASM 3.0;\nqubit[2] q; bit[2] c;\nh q[0]; h q[1];\ncz q[0],q[1];\nh q[0]; h q[1]; z q[0]; z q[1];\ncz q[0],q[1]; h q[0]; h q[1];\nc[0]=measure q[0]; c[1]=measure q[1];\n")
  end

  def run_qml(shots \\ @default_shots) do
    n = 4
    {counts, ms} = QxSim.run(n, shots, fn s ->
      s
      |> QxSim.rx(n, 0, 0.5) |> QxSim.rx(n, 1, 0.3)
      |> QxSim.cx(n, 0, 2)   |> QxSim.cx(n, 1, 3)
      |> QxSim.rx(n, 2, 0.7) |> QxSim.rx(n, 3, 0.7)
      |> QxSim.cx(n, 2, 3)
    end)
    resp("qml", n, shots, ms, counts,
      "Variational quantum classifier — Rx rotations + CNOT entanglement layer",
      "OPENQASM 3.0;\nqubit[4] q; bit[4] c;\nrx(0.5) q[0]; rx(0.3) q[1];\ncx q[0],q[2]; cx q[1],q[3];\nrx(0.7) q[2]; rx(0.7) q[3]; cx q[2],q[3];\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2]; c[3]=measure q[3];\n")
  end

  def run_shor15(shots \\ @default_shots) do
    # Shor's algorithm for N = 15 with a = 7. Counting register q0..q2 (q0 lowest), work register q3..q6 = |1⟩.
    # Controlled-U^(2^j), U|y⟩ = |7y mod 15⟩: U from q0 (three controlled swaps and four CNOTs), U² = ×4 from q1
    # (two controlled swaps), U⁴ = ×1, the identity, from q2. Then the inverse QFT on the counting register.
    # A measured k gives the phase k/8 = s/r; the period r of 7^x mod 15 is 4, so k ∈ {0, 2, 4, 6}; with r even
    # and 7^(r/2) ≠ −1 mod 15, gcd(7^(r/2) ± 1, 15) gives the factors 3 and 5.
    n = 7
    {counts, ms} = QxSim.run(n, shots, fn s ->
      s
      |> QxSim.h(n, 0) |> QxSim.h(n, 1) |> QxSim.h(n, 2)
      |> QxSim.x(n, 3)
      |> QxSim.cswap(n, 0, 3, 4) |> QxSim.cswap(n, 0, 4, 5) |> QxSim.cswap(n, 0, 5, 6)
      |> QxSim.cx(n, 0, 3) |> QxSim.cx(n, 0, 4) |> QxSim.cx(n, 0, 5) |> QxSim.cx(n, 0, 6)
      |> QxSim.cswap(n, 1, 4, 6) |> QxSim.cswap(n, 1, 3, 5)
      |> QxSim.swap(n, 0, 2)
      |> QxSim.h(n, 0)
      |> QxSim.cp(n, 0, 1, -:math.pi() / 2) |> QxSim.h(n, 1)
      |> QxSim.cp(n, 0, 2, -:math.pi() / 4) |> QxSim.cp(n, 1, 2, -:math.pi() / 2) |> QxSim.h(n, 2)
    end)
    # the counting register only (q0..q2); the work register is not read
    counting = Enum.reduce(counts, %{}, fn {bits, c}, acc -> Map.update(acc, String.slice(bits, -3, 3), c, &(&1 + c)) end)
    resp("shor15", 3, shots, ms, counting,
      "Shor's algorithm, N = 15, a = 7: the complete circuit compiled for a = 7 (controlled ×7 and ×4 mod 15), the inverse QFT and the classical step — 15 = 3 × 5",
      "OPENQASM 3.0;\n// Shor N=15, a=7: counting q[0..2], work q[3..6]\nqubit[7] q; bit[3] c;\nh q[0]; h q[1]; h q[2]; x q[3];\n// controlled ×7 mod 15 (U) from q[0]\ncswap q[0],q[3],q[4]; cswap q[0],q[4],q[5]; cswap q[0],q[5],q[6];\ncx q[0],q[3]; cx q[0],q[4]; cx q[0],q[5]; cx q[0],q[6];\n// controlled ×4 mod 15 (U²) from q[1]; U⁴ = identity\ncswap q[1],q[4],q[6]; cswap q[1],q[3],q[5];\n// inverse QFT on q[0..2]\nswap q[0],q[2];\nh q[0];\ncp(-pi/2) q[0],q[1]; h q[1];\ncp(-pi/4) q[0],q[2]; cp(-pi/2) q[1],q[2]; h q[2];\nc[0]=measure q[0]; c[1]=measure q[1]; c[2]=measure q[2];\n")
    |> Map.put(:classical, shor_classical(counting))
    |> Map.put(:qubits_total, n)
    |> then(fn r ->
      # how many shots give the factors: directly (k/8 in lowest terms has the period as denominator), and with
      # the classical check of the period's multiples (which rescues k = 4)
      direct = Enum.sum(for c <- r.classical, c.k in [2, 6], do: c.count)
      checked = Enum.sum(for c <- r.classical, c.factors != [], do: c.count)
      Map.merge(r, %{success_direct: Float.round(direct / shots, 4), success_with_check: Float.round(checked / shots, 4)})
    end)
  end

  # each measured k: the phase k/8, the period it gives, and the factors (or why it gives none)
  defp shor_classical(counting) do
    for {bits, c} <- Enum.sort(counting) do
      k = String.to_integer(bits, 2)
      g = Integer.gcd(k, 8)
      {r, note} =
        cond do
          k == 0 -> {nil, "phase 0: no information, run again"}
          true ->
            r0 = div(8, g)                                              # k/8 = s/r in lowest terms
            r = Enum.find([r0, 2 * r0, 4 * r0], &(rem(Integer.pow(7, &1), 15) == 1))
            {r, if(r == r0, do: "period #{r}", else: "#{k}/8 = 1/#{r0}: 7^#{r0} ≠ 1 mod 15, so the period is a multiple: #{r}")}
        end
      factors =
        if r != nil and rem(r, 2) == 0 and rem(Integer.pow(7, div(r, 2)), 15) != 14,
          do: [Integer.gcd(Integer.pow(7, div(r, 2)) - 1, 15), Integer.gcd(Integer.pow(7, div(r, 2)) + 1, 15)],
          else: []
      %{measured: bits, k: k, count: c, phase: "#{k}/8", period: r || :null, factors: factors, note: note}   # :null is JSON null
    end
  end

  @doc """
  Grover's search on n qubits (2..10) for one marked state. Each iteration is the oracle (X on the marked
  state's 0 bits, a multi-controlled Z, X back) and the diffuser (H, X, multi-controlled Z, X, H on every
  qubit). Returns the success probability after every iteration up to past the optimum, computed exactly from
  the state vector, beside the theory sin²((2k+1)θ), sin θ = 1/√N, and shots sampled after `iterations`
  (by default the optimum, round(π/(4θ) − 1/2)).
  """
  def run_grover_search(n, marked, iterations, shots) do
    t0 = System.monotonic_time(:microsecond)
    nn = 1 <<< n
    all = Enum.to_list(0..(n - 1))
    theta = :math.asin(1 / :math.sqrt(nn))
    best = max(round(:math.pi() / (4 * theta) - 0.5), 1)
    k_run = if iterations == :auto, do: best, else: iterations
    k_max = max(k_run, 2 * best + 1)
    flip0 = fn st -> Enum.reduce(all, st, fn q, a -> if bit(marked, q) == 0, do: QxSim.x(a, n, q), else: a end) end
    hs = fn st -> Enum.reduce(all, st, &QxSim.h(&2, n, &1)) end
    xs = fn st -> Enum.reduce(all, st, &QxSim.x(&2, n, &1)) end
    step = fn st ->
      st |> flip0.() |> QxSim.mcz(n, all) |> flip0.()                                   # oracle
         |> hs.() |> xs.() |> QxSim.mcz(n, all) |> xs.() |> hs.()                       # diffuser
    end
    s0 = n |> QxSim.init() |> hs.()
    {curve, states} =
      Enum.map_reduce(0..k_max, s0, fn k, st ->
        nxt = if k < k_max, do: step.(st), else: st
        {%{k: k, p: Float.round(QxSim.prob(st, marked), 6), theory: Float.round(:math.pow(:math.sin((2 * k + 1) * theta), 2), 6), state: st}, nxt}
      end)
    _ = states
    run_state = Enum.find(curve, &(&1.k == k_run)).state
    counts = QxSim.sample_low(run_state, n, shots)
    ms = Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3)
    target = Integer.to_string(marked, 2) |> String.pad_leading(n, "0")
    %{
      status: "ok", name: "grover_search", qubits: n, shots: shots, execution_time_ms: ms,
      marked: target, iterations: k_run, optimal_iterations: best,
      success: Float.round(Map.get(counts, target, 0) / shots, 4),
      curve: Enum.map(curve, &Map.delete(&1, :state)),
      counts: counts, visualisation_hint: "histogram", plots: %{counts: vega(counts)},
      circuit_qasm: "OPENQASM 3.0;\n// Grover, #{n} qubits, marked |#{target}⟩, #{k_run} iterations\nqubit[#{n}] q; bit[#{n}] c;\nh q;   // all qubits\n// repeated #{k_run} times:\n//   oracle: x on the marked state's 0 bits; a multi-controlled z on all qubits; x back\n//   diffuser: h q; x q; multi-controlled z; x q; h q;\nc = measure q;\n",
      metadata: %{description: "Grover's search on #{n} qubits: success probability per iteration, exact and in theory", version: "4.0", date: now()}
    }
  end

  @doc """
  Shor's algorithm for N = 15 or 21 and any a coprime to N, with nothing about the period put in: t = 2·⌈log₂N⌉
  counting qubits, a work register |1⟩, controlled multiplication by a^(2^j) mod N from counting qubit j
  (the classical powers by repeated squaring), the inverse QFT, then the classical step: continued fractions
  of k/2^t give candidate periods, each checked by a^r ≡ 1 (mod N), and gcd(a^(r/2) ± 1, N) the factors.
  """
  def run_shor(nn, a, shots) do
    t0 = System.monotonic_time(:microsecond)
    g = Integer.gcd(a, nn)
    if g > 1 do
      # no circuit needed: a shares a factor with N
      %{status: "ok", name: "shor", N: nn, a: a, shots: 0, execution_time_ms: 0.0, circuit_run: false,
        factors: Enum.sort([g, div(nn, g)]), success: 1.0, period_true: :null,
        note: "gcd(#{a}, #{nn}) = #{g}: a shares a factor with N, which gives #{g} × #{div(nn, g)} at once, without the quantum circuit",
        classical: [], counts: %{}, failures: %{}, metadata: %{date: now()}}
    else
      run_shor_circuit(nn, a, shots, t0)
    end
  end

  defp run_shor_circuit(nn, a, shots, t0) do
    w = if nn > 16, do: 5, else: 4
    t = 2 * w
    n = t + w
    work = Enum.to_list(t..(n - 1))
    powers = Enum.map_reduce(0..(t - 1), rem(a, nn), fn _j, c -> {c, rem(c * c, nn)} end) |> elem(0)
    circuit = fn s ->
      s = Enum.reduce(0..(t - 1), s, &QxSim.h(&2, n, &1)) |> QxSim.x(n, t)
      s = powers |> Enum.with_index() |> Enum.reduce(s, fn {c, j}, s -> QxSim.cmulmod(s, n, j, work, c, nn) end)
      s = Enum.reduce(0..(div(t, 2) - 1), s, fn j, s -> QxSim.swap(s, n, j, t - 1 - j) end)
      Enum.reduce(0..(t - 1), s, fn j, s ->
        s = Enum.reduce(0..(j - 1)//1, s, fn m, s -> QxSim.cp(s, n, m, j, -:math.pi() / :math.pow(2, j - m)) end)
        QxSim.h(s, n, j)
      end)
    end
    state = n |> QxSim.init() |> circuit.()
    # the exact distribution of the counting register, from the state vector
    mask = (1 <<< t) - 1
    exact =
      Enum.reduce(0..(tuple_size(state) - 1), %{}, fn i, acc ->
        pr = QxSim.prob(state, i)
        if pr > 1.0e-12, do: Map.update(acc, i &&& mask, pr, &(&1 + pr)), else: acc
      end)
    counting = QxSim.sample_low(state, t, shots)
    true_r = Enum.find(1..nn, &(rem(Integer.pow(a, &1), nn) == 1))
    classify = fn k ->
      r = shor_period(k, t, a, nn)
      cond do
        r == nil -> {nil, [], "no period found"}
        rem(r, 2) == 1 -> {r, [], "period odd"}
        rem(Integer.pow(a, div(r, 2)), nn) == nn - 1 -> {r, [], "a^(r/2) ≡ −1 (mod N)"}
        true ->
          f = Enum.sort([Integer.gcd(Integer.pow(a, div(r, 2)) - 1, nn), Integer.gcd(Integer.pow(a, div(r, 2)) + 1, nn)])
          if Enum.all?(f, &(&1 not in [1, nn])), do: {r, f, "factors"}, else: {r, [], "trivial factors"}
      end
    end
    rows =
      for {bits, c} <- Enum.sort_by(counting, fn {_, c} -> -c end) do
        k = String.to_integer(bits, 2)
        {r, f, outcome} = classify.(k)
        %{measured: bits, k: k, count: c, p_exact: Float.round(Map.get(exact, k, 0.0), 6), phase: "#{k}/#{1 <<< t}",
          period: r || :null, factors: f, outcome: outcome}
      end
    failures = Enum.reduce(rows, %{}, fn r, acc -> Map.update(acc, r.outcome, r.count, &(&1 + r.count)) end)
    # the exact success probability: the same classification over the exact distribution
    exact_success = Enum.reduce(exact, 0.0, fn {k, pr}, acc -> if elem(classify.(k), 2) == "factors", do: acc + pr, else: acc end)
    found = Map.get(failures, "factors", 0)
    ms = Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3)
    %{status: "ok", name: "shor", N: nn, a: a, qubits: n, counting_qubits: t, shots: shots, execution_time_ms: ms, circuit_run: true,
      period_true: true_r, success: Float.round(found / shots, 4), success_exact: Float.round(exact_success, 4),
      factors: (Enum.find(rows, &(&1.factors != [])) || %{factors: []}).factors,
      failures: failures, outcomes: length(rows), classical: Enum.take(rows, 32), counts: counting,
      metadata: %{description: "Shor's algorithm, N = #{nn}, a = #{a}: #{t} counting + #{w} work qubits. The controlled multiplication by a^(2^j) mod N is a permutation of basis states computed classically (values ≥ N left alone) — an oracle, not a gate-level arithmetic circuit, so no gate counts or resource estimates follow from it. Continued fractions on k/2^#{t}.",
                  date: now()}}
  end

  # the smallest denominator r ≤ N among the continued-fraction convergents of k/2^t with a^r ≡ 1 (mod N),
  # trying small multiples of each (a convergent can give a divisor of the period)
  defp shor_period(0, _t, _a, _nn), do: nil
  defp shor_period(k, t, a, nn) do
    convergents(k, 1 <<< t)
    |> Enum.map(&elem(&1, 1))
    |> Enum.filter(&(&1 > 0 and &1 <= nn))
    |> Enum.flat_map(fn q -> for m <- 1..4, q * m <= nn, do: q * m end)
    |> Enum.find(&(rem(Integer.pow(a, &1), nn) == 1))
  end

  # h₋₁ = 1, h₋₂ = 0 and k₋₁ = 0, k₋₂ = 1: the convergents hₙ/kₙ of p/q
  defp convergents(p, q), do: convergents(p, q, {1, 0}, {0, 1}, [])
  defp convergents(_p, 0, _h, _k, acc), do: Enum.reverse(acc)
  defp convergents(p, q, {h1, h2}, {k1, k2}, acc) do
    a = div(p, q)
    {h, k} = {a * h1 + h2, a * k1 + k2}
    convergents(q, p - a * q, {h, h1}, {k, k1}, [{h, k} | acc])
  end

  # ============================================================
  # QEC DEMOS
  # ============================================================

  def run_bitflip_repetition(shots \\ @default_shots, eq \\ "none") do
    n = 3; ns = 5
    {before_c, ms1} = QxSim.run(n, shots, fn s -> inject(s, n, eq, :x) end)
    {syn_c, ms2}    = QxSim.run(ns, shots, fn s ->
      s
      |> inject(ns, eq, :x)
      |> QxSim.cx(ns, 0, 3) |> QxSim.cx(ns, 1, 3)
      |> QxSim.cx(ns, 1, 4) |> QxSim.cx(ns, 2, 4)
    end)
    qec_resp("bitflip", n, shots, ms1 + ms2, before_c, syn_c, "bit_flip", eq,
      "OPENQASM 3.0;\n// 3-qubit bit-flip repetition code\nqubit[5] q; bit[2] s;\ncx q[0],q[3]; cx q[1],q[3];\ncx q[1],q[4]; cx q[2],q[4];\ns[0]=measure q[3]; s[1]=measure q[4];\n")
  end

  # ── Full correction cycles: phase-flip, Shor-9, Steane-7 ──────────────────────────────────────────────
  # encode logical |0⟩ → inject the error → measure every stabiliser check with one ancilla, reset after each
  # check → look the syndrome up in the decoder → apply the correction → read the logical qubit out in the
  # code's basis. The same circuit without the correction gives "before". With noiseless gates a single Pauli
  # error's syndrome is certain, so it is measured once on the state vector and holds for every shot.

  def run_phaseflip_repetition(shots \\ @default_shots, eq \\ "none") do
    # |0_L⟩ = |+++⟩: a repetition code in the X basis. Checks X0X1 and X1X2 find a Z error; Z corrects it.
    code = %{
      name: "phaseflip", data: 3, error: :z,
      encode: fn s, n -> s |> QxSim.cx(n, 0, 1) |> QxSim.cx(n, 0, 2) |> hs(n, [0, 1, 2]) end,
      checks: [{:x, [0, 1]}, {:x, [1, 2]}],
      decode: fn [a, b] -> flip(:z, %{{1, 0} => 0, {1, 1} => 1, {0, 1} => 2}[{a, b}]) end,
      readout: fn s, n -> hs(s, n, [0, 1, 2]) end,                      # X basis: |+⟩ reads 0
      logical: &majority/1, keep_data: true,
      qasm: "OPENQASM 3.0;\n// 3-qubit phase-flip code: encode |0_L⟩ = |+++⟩, checks X0X1, X1X2 on a reused ancilla q[3]\nqubit[4] q; bit[2] s; bit[3] c;\ncx q[0],q[1]; cx q[0],q[2]; h q[0]; h q[1]; h q[2];\n// error: z q[e]\nh q[3]; cx q[3],q[0]; cx q[3],q[1]; h q[3]; s[0]=measure q[3]; reset q[3];\nh q[3]; cx q[3],q[1]; cx q[3],q[2]; h q[3]; s[1]=measure q[3]; reset q[3];\n// correction: z on the qubit the syndrome names (10→0, 11→1, 01→2)\nh q[0]; h q[1]; h q[2]; c = measure q[0:2];\n"
    }
    qec_cycle(code, shots, eq, "phase_flip")
  end

  def run_shor9(shots \\ @default_shots, eq \\ "none") do
    # Shor's code: a phase-flip code over three bit-flip codes. Six Z checks (two per block) find the bit
    # flip, two X checks across blocks find the phase flip. The injected error is Y: both at once.
    zsyn = %{{1, 0} => 0, {1, 1} => 1, {0, 1} => 2}
    code = %{
      name: "shor9", data: 9, error: :y,
      encode: fn s, n ->
        s |> QxSim.cx(n, 0, 3) |> QxSim.cx(n, 0, 6) |> hs(n, [0, 3, 6])
          |> QxSim.cx(n, 0, 1) |> QxSim.cx(n, 0, 2) |> QxSim.cx(n, 3, 4)
          |> QxSim.cx(n, 3, 5) |> QxSim.cx(n, 6, 7) |> QxSim.cx(n, 6, 8)
      end,
      checks: [{:z, [0, 1]}, {:z, [1, 2]}, {:z, [3, 4]}, {:z, [4, 5]}, {:z, [6, 7]}, {:z, [7, 8]},
               {:x, [0, 1, 2, 3, 4, 5]}, {:x, [3, 4, 5, 6, 7, 8]}],
      decode: fn [a0, b0, a1, b1, a2, b2, t1, t2] ->
        bit = [{a0, b0}, {a1, b1}, {a2, b2}] |> Enum.with_index()
              |> Enum.flat_map(fn {ab, blk} -> case zsyn[ab] do nil -> []; q -> [{:x, 3 * blk + q}] end end)
        phase = case zsyn[{t1, t2}] do nil -> []; blk -> [{:z, 3 * blk}] end
        bit ++ phase
      end,
      readout: fn s, n -> hs(s, n, Enum.to_list(0..8)) end,             # X basis; logical = parity of all 9
      logical: &parity/1, keep_data: false,
      qasm: "OPENQASM 3.0;\n// Shor 9-qubit code: encode |0_L⟩, six Z checks and two X checks on a reused ancilla q[9]\nqubit[10] q; bit[8] s; bit[9] c;\ncx q[0],q[3]; cx q[0],q[6]; h q[0]; h q[3]; h q[6];\ncx q[0],q[1]; cx q[0],q[2]; cx q[3],q[4]; cx q[3],q[5]; cx q[6],q[7]; cx q[6],q[8];\n// error: y q[e]\n// Z checks (bit flips), e.g. Z0Z1:\ncx q[0],q[9]; cx q[1],q[9]; s[0]=measure q[9]; reset q[9];\n// ... Z1Z2, Z3Z4, Z4Z5, Z6Z7, Z7Z8 alike\n// X checks (phase flips), X0..X5 and X3..X8:\nh q[9]; cx q[9],q[0]; cx q[9],q[1]; cx q[9],q[2]; cx q[9],q[3]; cx q[9],q[4]; cx q[9],q[5]; h q[9]; s[6]=measure q[9]; reset q[9];\nh q[9]; cx q[9],q[3]; cx q[9],q[4]; cx q[9],q[5]; cx q[9],q[6]; cx q[9],q[7]; cx q[9],q[8]; h q[9]; s[7]=measure q[9]; reset q[9];\n// correction: x on the flipped qubit of a block, z on the first qubit of the flipped block\nh q[0:8]; c = measure q[0:8];   // logical bit = parity of c\n"
    }
    qec_cycle(code, shots, eq, "bit_phase")
  end

  def run_steane7(shots \\ @default_shots, eq \\ "none") do
    # Steane's [[7,1,3]] code from the Hamming code: three Z checks and three X checks on the same supports.
    # Each syndrome, read as a binary number, is the position (1..7) of the error. The error is Y.
    g = [[0, 2, 4, 6], [1, 2, 5, 6], [3, 4, 5, 6]]                     # positions 1,3,5,7 · 2,3,6,7 · 4,5,6,7
    pos = fn [b0, b1, b2] -> b0 + 2 * b1 + 4 * b2 end
    code = %{
      name: "steane7", data: 7, error: :y,
      encode: fn s, n ->
        s |> hs(n, [0, 1, 3])
          |> QxSim.cx(n, 3, 4) |> QxSim.cx(n, 3, 5) |> QxSim.cx(n, 3, 6)
          |> QxSim.cx(n, 1, 2) |> QxSim.cx(n, 1, 5) |> QxSim.cx(n, 1, 6)
          |> QxSim.cx(n, 0, 2) |> QxSim.cx(n, 0, 4) |> QxSim.cx(n, 0, 6)
      end,
      checks: Enum.map(g, &{:z, &1}) ++ Enum.map(g, &{:x, &1}),
      decode: fn bits ->
        {zb, xb} = Enum.split(bits, 3)
        (if pos.(zb) > 0, do: [{:x, pos.(zb) - 1}], else: []) ++
        (if pos.(xb) > 0, do: [{:z, pos.(xb) - 1}], else: [])
      end,
      readout: fn s, _n -> s end,                                        # Z basis; logical = parity of all 7
      logical: &parity/1, keep_data: false,
      qasm: "OPENQASM 3.0;\n// Steane [[7,1,3]] code: encode |0_L⟩, three Z and three X checks on a reused ancilla q[7]\nqubit[8] q; bit[6] s; bit[7] c;\nh q[0]; h q[1]; h q[3];\ncx q[3],q[4]; cx q[3],q[5]; cx q[3],q[6]; cx q[1],q[2]; cx q[1],q[5]; cx q[1],q[6]; cx q[0],q[2]; cx q[0],q[4]; cx q[0],q[6];\n// error: y q[e]\n// Z checks on {0,2,4,6}, {1,2,5,6}, {3,4,5,6}, e.g.:\ncx q[0],q[7]; cx q[2],q[7]; cx q[4],q[7]; cx q[6],q[7]; s[0]=measure q[7]; reset q[7];\n// X checks on the same supports, e.g.:\nh q[7]; cx q[7],q[0]; cx q[7],q[2]; cx q[7],q[4]; cx q[7],q[6]; h q[7]; s[3]=measure q[7]; reset q[7];\n// correction: x at position s[0]+2s[1]+4s[2], z at position s[3]+2s[4]+4s[5]\nc = measure q[0:6];   // logical bit = parity of c\n"
    }
    qec_cycle(code, shots, eq, "bit_phase")
  end

  defp qec_cycle(code, shots, eq, error_type) do
    t0 = System.monotonic_time(:microsecond)
    k = code.data; n = k + 1; anc = k
    err = parse_eq(eq, k)
    s0 = n |> QxSim.init() |> code.encode.(n)
    # a code made by measurement: one round of checks first, its outcomes the reference
    {ref, s0} =
      if Map.get(code, :reference, false),
        do: Enum.map_reduce(code.checks, s0, fn chk, s -> check(s, n, anc, chk) end),
        else: {List.duplicate(0, length(code.checks)), s0}
    s1 = if err, do: pauli(s0, n, code.error, err), else: s0
    {raw, s2} = Enum.map_reduce(code.checks, s1, fn chk, s -> check(s, n, anc, chk) end)
    syn = Enum.zip_with(raw, ref, &Bitwise.bxor/2)
    fixes = code.decode.(syn)
    s3 = Enum.reduce(fixes, s2, fn {gate, q}, s -> pauli(s, n, gate, q) end)
    before = s2 |> code.readout.(n) |> QxSim.sample_low(k, shots)
    fixed  = s3 |> code.readout.(n) |> QxSim.sample_low(k, shots)
    logical = fn counts -> Enum.reduce(counts, %{}, fn {b, c}, a -> Map.update(a, code.logical.(b), c, &(&1 + c)) end) end
    {lb, la} = {logical.(before), logical.(fixed)}
    {cb, ca} = if code.keep_data, do: {before, fixed}, else: {lb, la}
    syn_s = Enum.join(syn)
    ms = Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3)
    %{
      status: "ok", name: code.name, qubits: n, shots: shots, execution_time_ms: ms,
      error_type: error_type, error_qubit: eq,
      counts_before: cb, syndrome_counts: %{syn_s => shots}, corrected_counts: ca,
      logical_before: lb, logical_after: la, data_before: before, data_after: fixed,
      logical_fidelity: Float.round(Map.get(la, "0", 0) / shots, 4),
      syndrome: syn_s,
      correction: Enum.map(fixes, fn {g, q} -> "#{g} q[#{q}]" end),
      visualisation_hint: "qec",
      plots: %{counts_before: vega(cb), corrected_counts: vega(ca), syndrome: vega(%{syn_s => shots})},
      circuit_qasm: code.qasm,
      metadata: %{
        description: "#{code.name}: encode |0_L⟩, inject #{code.error}, measure the stabiliser checks with a reused ancilla, decode the syndrome, correct, read the logical qubit",
        code_distance: 3, error_correction: "syndrome decoding (lookup table)",
        readout: if(code.keep_data, do: "data qubits in the X basis", else: "logical qubit: 0 or 1"),
        version: "4.0", date: now()
      }
    }
  end

  # one stabiliser check on the ancilla: Z-type copies the parity in; X-type runs it in the Hadamard basis
  defp check(s, n, anc, {:z, qs}) do
    s |> then(&Enum.reduce(qs, &1, fn q, acc -> QxSim.cx(acc, n, q, anc) end)) |> QxSim.measure_reset(n, anc)
  end
  defp check(s, n, anc, {:x, qs}) do
    s |> QxSim.h(n, anc)
      |> then(&Enum.reduce(qs, &1, fn q, acc -> QxSim.cx(acc, n, anc, q) end))
      |> QxSim.h(n, anc) |> QxSim.measure_reset(n, anc)
  end

  defp hs(s, n, qs), do: Enum.reduce(qs, s, fn q, acc -> QxSim.h(acc, n, q) end)
  defp pauli(s, n, :x, q), do: QxSim.x(s, n, q)
  defp pauli(s, n, :y, q), do: QxSim.y(s, n, q)
  defp pauli(s, n, :z, q), do: QxSim.z(s, n, q)
  defp flip(_gate, nil), do: []
  defp flip(gate, q), do: [{gate, q}]
  defp parity(bits), do: if(rem(Enum.count(String.graphemes(bits), &(&1 == "1")), 2) == 0, do: "0", else: "1")
  defp majority(bits), do: if(Enum.count(String.graphemes(bits), &(&1 == "1")) * 2 > String.length(bits), do: "1", else: "0")

  def run_surface3(shots \\ @default_shots, eq \\ "none") do
    # The rotated surface code of distance 3: data q0..q8 on a 3×3 grid (row by row), eight stabilizers —
    # X on {0,1,3,4} {4,5,7,8} {1,2} {6,7}, Z on {1,2,4,5} {3,4,6,7} {0,3} {5,8} — measured one after another on
    # one ancilla, q9, reset after each. Logical |0⟩ is made by measurement: from |0…0⟩ (every Z check +1, the
    # logical Z = Z0Z1Z2 = +1) one round of checks projects onto the code; its X-check outcomes are random and
    # kept as the reference, and a syndrome is a change from it. The injected error is Y. Each single-qubit
    # error has its own check signature (or one equal up to a stabilizer), so a lookup table decodes it.
    xcs = [[0, 1, 3, 4], [4, 5, 7, 8], [1, 2], [6, 7]]
    zcs = [[1, 2, 4, 5], [3, 4, 6, 7], [0, 3], [5, 8]]
    sig = fn cs, q -> Enum.map(cs, &if(q in &1, do: 1, else: 0)) end
    xfix = Enum.reduce(8..0//-1, %{}, fn q, m -> Map.put(m, sig.(zcs, q), q) end)     # Z checks see X errors
    zfix = Enum.reduce(8..0//-1, %{}, fn q, m -> Map.put(m, sig.(xcs, q), q) end)     # X checks see Z errors
    code = %{
      name: "surface3", data: 9, error: :y, reference: true,
      encode: fn s, _n -> s end,
      checks: Enum.map(xcs, &{:x, &1}) ++ Enum.map(zcs, &{:z, &1}),
      decode: fn bits ->
        {xb, zb} = Enum.split(bits, 4)
        (if Enum.any?(zb, &(&1 == 1)), do: [{:x, xfix[zb]}], else: []) ++
        (if Enum.any?(xb, &(&1 == 1)), do: [{:z, zfix[xb]}], else: [])
      end,
      readout: fn s, _n -> s end,                                         # Z basis; logical = parity of q0, q1, q2
      logical: fn bits -> parity(String.slice(bits, -3, 3)) end, keep_data: false,
      qasm: "OPENQASM 3.0;\n// Rotated surface code, d = 3: data q[0..8] (3×3, row by row), one reused ancilla q[9]\nqubit[10] q; bit[16] s; bit[9] c;\n// X checks {0,1,3,4} {4,5,7,8} {1,2} {6,7}; Z checks {1,2,4,5} {3,4,6,7} {0,3} {5,8}\n// round 1 (reference, projects |0…0⟩ onto logical |0⟩), e.g. the X check {0,1,3,4}:\nh q[9]; cx q[9],q[0]; cx q[9],q[1]; cx q[9],q[3]; cx q[9],q[4]; h q[9]; s[0]=measure q[9]; reset q[9];\n// and the Z check {1,2,4,5}:\ncx q[1],q[9]; cx q[2],q[9]; cx q[4],q[9]; cx q[5],q[9]; s[4]=measure q[9]; reset q[9];\n// error: y q[e]; round 2 alike, s[8..15]; syndrome = round 2 xor round 1\n// correction: x at the qubit whose Z-check signature matches, z at the one whose X-check signature matches\nc = measure q[0:8];   // logical bit = c[0] xor c[1] xor c[2]\n"
    }
    qec_cycle(code, shots, eq, "bit_phase")
  end

  # ============================================================
  # PRIVATE HELPERS
  # ============================================================

  defp inject(state, n, eq, gate) do
    case parse_eq(eq, n) do
      nil -> state
      q   ->
        case gate do
          :x -> QxSim.x(state, n, q)
          :z -> QxSim.z(state, n, q)
        end
    end
  end

  defp parse_eq("none", _n), do: nil
  defp parse_eq(q, n) when is_binary(q) do
    case Integer.parse(q) do
      {i, _} when i >= 0 and i < n -> i
      _ -> nil
    end
  end
  defp parse_eq(_, _), do: nil

  defp majority_decode(counts) do
    Enum.reduce(counts, %{}, fn {bits, count}, acc ->
      width = String.length(bits)
      ones  = Enum.count(String.graphemes(bits), &(&1 == "1"))
      key   = if ones > div(width, 2),
        do:   String.duplicate("1", width),
        else: String.duplicate("0", width)
      Map.update(acc, key, count, &(&1 + count))
    end)
  end

  defp fidelity(corrected) do
    total = corrected |> Map.values() |> Enum.sum()
    if total == 0 do
      0.0
    else
      width = corrected |> Map.keys() |> List.first("") |> String.length()
      zero  = String.duplicate("0", width)
      Float.round(Map.get(corrected, zero, 0) / total, 4)
    end
  end

  defp vega(counts) do
    %{
      mark: "bar",
      data: %{values: Enum.map(counts, fn {s, c} -> %{state: s, count: c} end)},
      encoding: %{
        x: %{field: "state", type: "nominal",      title: "State"},
        y: %{field: "count", type: "quantitative", title: "Counts"}
      }
    }
  end

  defp now, do: DateTime.utc_now() |> DateTime.to_iso8601()

  defp resp(name, qubits, shots, ms, counts, description, qasm) do
    %{
      status: "ok",
      name: name,
      qubits: qubits,
      shots: shots,
      execution_time_ms: ms,
      counts: counts,
      visualisation_hint: "histogram",
      plots: %{counts: vega(counts)},
      circuit_qasm: qasm,
      metadata: %{
        description: description,
        version: "3.0",
        improvements: ["tuple_state_vector", "atomic_gate_update", "binary_search_sampling"],
        date: now()
      }
    }
  end

  defp qec_resp(name, qubits, shots, ms, before_c, syn_c, error_type, error_qubit, qasm) do
    corrected = majority_decode(before_c)
    fid       = fidelity(corrected)
    %{
      status: "ok",
      name: name,
      qubits: qubits,
      shots: shots,
      execution_time_ms: ms,
      error_type: error_type,
      error_qubit: error_qubit,
      counts_before: before_c,
      syndrome_counts: syn_c,
      corrected_counts: corrected,
      logical_fidelity: fid,
      visualisation_hint: "qec",
      plots: %{
        counts_before:    vega(before_c),
        corrected_counts: vega(corrected),
        syndrome:         vega(syn_c)
      },
      circuit_qasm: qasm,
      metadata: %{
        description: "#{name} QEC demo",
        code_distance: 3,
        error_correction: "majority_vote",
        version: "3.0",
        date: now()
      }
    }
  end
end
