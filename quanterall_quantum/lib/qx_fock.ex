defmodule QxFock do
  @moduledoc """
  Linear optics with photons, in Fock space — for simulating photonic schemes from the literature (KLM gates,
  fusion, detectors). Conventions: docs/optics_conventions.md. The Python twin is quantum/fock.py in the
  TechUPGRADE engine (and ~/quanterall/qsim); Perceval checks both (tools/optics/).

  A state is `%{modes: m, amps: %{occupation_tuple => {re, im}}}`. An interferometer is an m × m unitary M, given
  as a list of rows of `{re, im}`, acting on creation operators: a_k† → Σ_j M[j][k] a_j†. `linear/2` expands the
  creation operators photon by photon; `apply_perm/2` computes each output amplitude as a permanent (Ryser) —
  two routes that must agree.
  """

  # ── complex arithmetic on {re, im} ─────────────────────────────────────────────────────────────────
  defp cmul({a, b}, {c, d}), do: {a * c - b * d, a * d + b * c}
  defp cadd({a, b}, {c, d}), do: {a + c, b + d}
  defp cscale({a, b}, s), do: {a * s, b * s}
  defp abs2({a, b}), do: a * a + b * b

  def basis(occ), do: %{modes: tuple_size(occ), amps: %{occ => {1.0, 0.0}}}

  def add(%{modes: m, amps: a}, %{amps: b}), do: %{modes: m, amps: Map.merge(a, b, fn _k, x, y -> cadd(x, y) end)}
  def scale(%{amps: a} = s, c), do: %{s | amps: Map.new(a, fn {k, v} -> {k, cmul(v, c)} end)}
  def norm2(%{amps: a}), do: a |> Map.values() |> Enum.map(&abs2/1) |> Enum.sum()

  defp fact(0), do: 1
  defp fact(n), do: n * fact(n - 1)
  defp sqfact(occ), do: occ |> Tuple.to_list() |> Enum.map(&fact/1) |> Enum.product() |> :math.sqrt()

  defp prune(amps, tol \\ 1.0e-14), do: amps |> Enum.reject(fn {_k, v} -> abs2(v) < tol * tol end) |> Map.new()

  @doc "Applies the interferometer M (linear optics) by expanding the creation operators photon by photon."
  def linear(%{modes: m, amps: amps} = st, mat) do
    mt = mat |> Enum.map(&List.to_tuple/1) |> List.to_tuple()
    out =
      Enum.reduce(amps, %{}, fn {t, a}, acc ->
        start = %{List.to_tuple(List.duplicate(0, m)) => cscale(a, 1 / sqfact(t))}
        poly =
          t |> Tuple.to_list() |> Enum.with_index()
          |> Enum.reduce(start, fn {tk, k}, poly ->
            Enum.reduce(1..tk//1, poly, fn _, poly ->
              # one more photon from input mode k: it can leave by any output mode j with amplitude M[j][k]
              Enum.reduce(poly, %{}, fn {mono, c}, nxt ->
                Enum.reduce(0..(m - 1), nxt, fn j, nxt ->
                  mjk = elem(elem(mt, j), k)
                  if abs2(mjk) > 0,
                    do: Map.update(nxt, put_elem(mono, j, elem(mono, j) + 1), cmul(c, mjk), &cadd(&1, cmul(c, mjk))),
                    else: nxt
                end)
              end)
            end)
          end)
        Enum.reduce(poly, acc, fn {s, c}, acc ->
          v = cscale(c, sqfact(s))
          Map.update(acc, s, v, &cadd(&1, v))
        end)
      end)
    %{st | amps: prune(out)}
  end

  @doc "Applies M by computing every output amplitude as a permanent (Ryser's formula)."
  def apply_perm(%{modes: m, amps: amps} = st, mat) do
    mt = mat |> Enum.map(&List.to_tuple/1) |> List.to_tuple()
    out =
      Enum.reduce(amps, %{}, fn {t, a}, acc ->
        n = t |> Tuple.to_list() |> Enum.sum()
        cols = expand(t)
        Enum.reduce(occupations(m, n), acc, fn s, acc ->
          rows = expand(s)
          sub = for r <- rows, do: for(c <- cols, do: elem(elem(mt, r), c))
          v = cscale(permanent(sub), 1 / (sqfact(s) * sqfact(t)))
          if abs2(v) < 1.0e-28, do: acc, else: Map.update(acc, s, cmul(a, v), &cadd(&1, cmul(a, v)))
        end)
      end)
    %{st | amps: prune(out)}
  end

  defp expand(occ), do: occ |> Tuple.to_list() |> Enum.with_index() |> Enum.flat_map(fn {c, i} -> List.duplicate(i, c) end)

  @doc "Every way to put n photons into m modes."
  def occupations(1, n), do: [{n}]
  def occupations(m, n), do: for(k <- n..0//-1, rest <- occupations(m - 1, n - k), do: Tuple.insert_at(rest, 0, k))

  @doc "Ryser's formula: Perm(A) = (−1)ⁿ Σ_S (−1)^|S| Π_i Σ_{j∈S} A[i][j]."
  def permanent([]), do: {1.0, 0.0}
  def permanent(a) do
    n = length(a)
    total =
      Enum.reduce(1..(Bitwise.bsl(1, n) - 1), {0.0, 0.0}, fn mask, acc ->
        cols = for j <- 0..(n - 1), Bitwise.band(mask, Bitwise.bsl(1, j)) != 0, do: j
        prod = Enum.reduce(a, {1.0, 0.0}, fn row, p -> cmul(p, Enum.reduce(cols, {0.0, 0.0}, &cadd(&2, Enum.at(row, &1)))) end)
        sign = if rem(length(cols), 2) == 0, do: 1.0, else: -1.0
        cadd(acc, cscale(prod, sign))
      end)
    cscale(total, if(rem(n, 2) == 0, do: 1.0, else: -1.0))
  end

  # ── elements ─────────────────────────────────────────────────────────────────────────────────────────
  @doc "Identity m × m, with a 2 × 2 or 3 × 3 block placed on the given modes."
  def embed(m, idx, block) do
    pos = idx |> Enum.with_index() |> Map.new()
    for j <- 0..(m - 1) do
      for k <- 0..(m - 1) do
        case {Map.get(pos, j), Map.get(pos, k)} do
          {nil, _} -> if j == k, do: {1.0, 0.0}, else: {0.0, 0.0}
          {_, nil} -> if j == k, do: {1.0, 0.0}, else: {0.0, 0.0}
          {r, c} -> Enum.at(Enum.at(block, r), c)
        end
      end
    end
  end

  @doc "Beamsplitter: a_i† → cos θ a_i† − sin θ a_j†, a_j† → sin θ a_i† + cos θ a_j†."
  def beamsplitter(m, i, j, theta) do
    {c, s} = {:math.cos(theta), :math.sin(theta)}
    embed(m, [i, j], [[{c, 0.0}, {s, 0.0}], [{-s, 0.0}, {c, 0.0}]])
  end

  def bs(%{modes: m} = st, i, j, theta), do: linear(st, beamsplitter(m, i, j, theta))

  def phase(%{amps: a} = st, k, phi) do
    %{st | amps: Map.new(a, fn {o, v} -> {o, cmul(v, {:math.cos(phi * elem(o, k)), :math.sin(phi * elem(o, k))})} end)}
  end

  @doc "The unnormalised state of the other modes given `counts` on `modes`; its norm² is the outcome's probability."
  def project(%{modes: m, amps: a}, modes, counts) do
    keep = Enum.reject(0..(m - 1), &(&1 in modes))
    out =
      for {o, v} <- a, Enum.zip(modes, counts) |> Enum.all?(fn {k, c} -> elem(o, k) == c end), reduce: %{} do
        acc -> Map.update(acc, List.to_tuple(for(k <- keep, do: elem(o, k))), v, &cadd(&1, v))
      end
    %{modes: length(keep), amps: out}
  end

  # ── the KLM gates (Nature 409, 46 (2001), Eq. 9; conventions file) ───────────────────────────────────
  @deg :math.pi() / 180

  @doc "The NS gate on `signal` with ancilla modes `anc` = {a1, a2} (holding |1, 0⟩), in the beamsplitter form."
  def ns_matrix(m, signal, {a1, a2}, opts \\ []) do
    {t1, t2, t3} = Keyword.get(opts, :angles, {22.5, 65.5302, -22.5})
    b = fn i, j, t -> beamsplitter(3, i, j, t * @deg) end
    u3 = b.(1, 2, t3) |> matmul(b.(0, 1, t2)) |> matmul(b.(1, 2, t1))
    u3 = if Keyword.get(opts, :phase_pi, true), do: matmul(u3, embed(3, [0], [[{-1.0, 0.0}]])), else: u3
    embed(m, [signal, a1, a2], u3)
  end

  def matmul(a, b) do
    bt = b |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
    for row <- a, do: for(col <- bt, do: Enum.zip(row, col) |> Enum.reduce({0.0, 0.0}, fn {x, y}, acc -> cadd(acc, cmul(x, y)) end))
  end

  # ── demos behind the endpoints ───────────────────────────────────────────────────────────────────────

  @doc "HOM: coincidence probability for photons overlapping by e^{−γ|τ|/2} (exponential pulses), at delay γτ."
  def hom(gamma_tau) do
    o = :math.exp(-abs(gamma_tau) / 2)
    st = add(scale(basis({1, 0, 1, 0}), {o, 0.0}), scale(basis({1, 0, 0, 1}), {:math.sqrt(1 - o * o), 0.0}))
    st = st |> bs(0, 2, :math.pi() / 4) |> bs(1, 3, :math.pi() / 4)
    c = for {k, v} <- st.amps, elem(k, 0) + elem(k, 1) == 1, reduce: 0.0, do: (acc -> acc + abs2(v))
    %{gamma_tau: gamma_tau, coincidence: c, formula: 0.5 * (1 - :math.exp(-abs(gamma_tau)))}
  end

  @doc "Endpoint: HOM coincidence at one delay, and the curve around it."
  def run_hom(gt) do
    t0 = System.monotonic_time(:microsecond)
    curve = for x <- -40..40, do: hom(x / 10)
    %{status: "ok", name: "optics_hom", point: hom(gt), curve: curve,
      execution_time_ms: Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3),
      metadata: %{description: "Hong–Ou–Mandel: two photons on a balanced beamsplitter, overlapping by e^(−γ|τ|/2) (exponential pulses); coincidence probability against ½(1 − e^(−γ|τ|)). A simulation of the scheme, not a measurement.",
                  modes: "(a0, a1, b0, b1): two ports × two internal modes"}}
  end

  @doc "Endpoint: the NS gate on a random input α0|0⟩ + α1|1⟩ + α2|2⟩ (from the seed): output, success, fidelity."
  def run_ns(opts \\ []) do
    al = for _ <- 0..2, do: {:rand.normal(), :rand.normal()}
    nrm = :math.sqrt(Enum.sum(Enum.map(al, &abs2/1)))
    al = Enum.map(al, &cscale(&1, 1 / nrm))
    st = Enum.with_index(al) |> Enum.map(fn {a, n} -> scale(basis({n, 1, 0}), a) end) |> Enum.reduce(&add/2)
    h = st |> linear(ns_matrix(3, 0, {1, 2}, opts)) |> project([1, 2], [1, 0])
    out = for n <- 0..2, do: Map.get(h.amps, {n}, {0.0, 0.0})
    p = Enum.sum(Enum.map(out, &abs2/1))
    target = Enum.zip_with(al, [1.0, 1.0, -1.0], &cscale/2)
    ov = Enum.zip(target, out) |> Enum.reduce({0.0, 0.0}, fn {{a, b}, y}, acc -> cadd(acc, cmul({a, -b}, y)) end)
    pair = fn {re, im} -> [Float.round(re, 9), Float.round(im, 9)] end
    %{status: "ok", name: "optics_ns", input: Enum.map(al, pair), output_heralded: Enum.map(out, pair),
      success: p, fidelity: abs2(ov) / p,
      metadata: %{description: "KLM's nonlinear sign gate: ancilla |1, 0⟩, herald (1, 0); the output should be α0|0⟩ + α1|1⟩ − α2|2⟩ with probability 1/4"}}
  end

  @doc """
  Endpoint: the CZ circuit sampled like an experiment — for a logical input (00, 01, 10, 11), `shots` runs, each
  reading the four ancilla detectors; heralded runs give the logical output. Beside it the exact process fidelity
  and success. `opts`: angles (a power check), phase_pi.
  """
  def run_cz(input, shots, opts \\ []) do
    t0 = System.monotonic_time(:microsecond)
    exact = cz(opts)
    {q1, q2} = {String.to_integer(String.at(input, 0)), String.to_integer(String.at(input, 1))}
    st = basis({1 - q1, q1, 1 - q2, q2, 1, 0, 1, 0}) |> bs(1, 3, :math.pi() / 4)
    st = st |> linear(ns_matrix(8, 1, {4, 5}, opts)) |> linear(ns_matrix(8, 3, {6, 7}, opts)) |> bs(1, 3, -:math.pi() / 4)
    outcomes = st.amps |> Enum.map(fn {k, v} -> {k, abs2(v)} end) |> Enum.sort()
    {cdf, _} = Enum.map_reduce(outcomes, 0.0, fn {k, p}, acc -> {{k, acc + p}, acc + p} end)
    samples = for _ <- 1..shots, do: (r = :rand.uniform_real(); Enum.find(cdf, fn {_k, c} -> r <= c end) || List.last(cdf)) |> elem(0)
    heralded = Enum.filter(samples, fn k -> {elem(k, 4), elem(k, 5), elem(k, 6), elem(k, 7)} == {1, 0, 1, 0} end)
    logical = fn k ->
      case {elem(k, 0), elem(k, 1), elem(k, 2), elem(k, 3)} do
        {1, 0, 1, 0} -> "00"; {1, 0, 0, 1} -> "01"; {0, 1, 1, 0} -> "10"; {0, 1, 0, 1} -> "11"; _ -> "outside dual rail"
      end
    end
    %{status: "ok", name: "optics_cz", input: input, shots: shots, heralded: length(heralded),
      herald_rate: length(heralded) / shots, heralded_outputs: Enum.frequencies_by(heralded, logical),
      exact: exact, execution_time_ms: Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3),
      metadata: %{description: "CZ from two KLM NS gates on dual-rail qubits (8 modes, 4 photons): sampled detector readings; heralded when the ancillas read 1, 0, 1, 0. Expected: herald rate 1/16; the logical value unchanged (CZ only adds a sign)"}}
  end

  @doc """
  Type-I fusion (Browne & Rudolph 2005) on two polarisation photons in a random state (from the seed), sampled:
  a polarising beamsplitter (V modes swap ports), a 45° rotation on port 2, a polarisation-resolving detector there.
  One photon on port 2 is success — the fused qubit stays on port 1, acted on by (|0⟩⟨00| ± |1⟩⟨11|)/√2; zero or two
  is failure (a Z measurement of both). Returns the outcome counts and, exactly, the success probability against the
  input's weight on |HH⟩ and |VV⟩.
  """
  def run_fusion(shots) do
    psi = for _ <- 0..3, do: {:rand.normal(), :rand.normal()}
    nrm = :math.sqrt(Enum.sum(Enum.map(psi, &abs2/1)))
    psi = Enum.map(psi, &cscale(&1, 1 / nrm))
    occ = [{1, 0, 1, 0}, {1, 0, 0, 1}, {0, 1, 1, 0}, {0, 1, 0, 1}]                 # HH, HV, VH, VV
    st = Enum.zip(occ, psi) |> Enum.map(fn {o, a} -> scale(basis(o), a) end) |> Enum.reduce(&add/2)
    pbs = [[{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}], [{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}],
           [{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}], [{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}]]
    st = st |> linear(pbs) |> bs(2, 3, :math.pi() / 4)
    probs = st.amps |> Enum.map(fn {k, v} -> {k, abs2(v)} end) |> Enum.sort()
    {cdf, _} = Enum.map_reduce(probs, 0.0, fn {k, p}, acc -> {{k, acc + p}, acc + p} end)
    label = fn k ->
      case {elem(k, 2), elem(k, 3)} do
        {1, 0} -> "success +"
        {0, 1} -> "success −"
        {a, b} when a + b == 0 -> "failure: no photon"
        _ -> "failure: two photons"
      end
    end
    counts = (for _ <- 1..shots, do: (r = :rand.uniform_real(); Enum.find(cdf, fn {_k, c} -> r <= c end) || List.last(cdf)) |> elem(0) |> label.())
             |> Enum.frequencies()
    exact = probs |> Enum.filter(fn {k, _} -> elem(k, 2) + elem(k, 3) == 1 end) |> Enum.map(&elem(&1, 1)) |> Enum.sum()
    %{status: "ok", name: "optics_fusion", shots: shots, counts: counts, success_exact: exact,
      hh_plus_vv_weight: abs2(Enum.at(psi, 0)) + abs2(Enum.at(psi, 3)),
      metadata: %{description: "Type-I fusion of two polarisation photons in a random two-qubit state: the success probability equals the state's weight on HH and VV (1/2 for cluster qubits); failure measures both in Z"}}
  end

  @doc """
  CZ from two NS gates on dual-rail qubits (modes 0,1 and 2,3; ancillas 4,5 and 6,7): for each logical input, the
  heralded output (ancillas reading 1,0,1,0) as a 4 × 4 map; its process fidelity to CZ and its success probability.
  """
  def cz(opts \\ []) do
    inputs = [{0, 0}, {0, 1}, {1, 0}, {1, 1}]
    cols =
      for {q1, q2} <- inputs do
        st = basis({1 - q1, q1, 1 - q2, q2, 1, 0, 1, 0}) |> bs(1, 3, :math.pi() / 4)
        st = st |> linear(ns_matrix(8, 1, {4, 5}, opts)) |> linear(ns_matrix(8, 3, {6, 7}, opts)) |> bs(1, 3, -:math.pi() / 4)
        h = project(st, [4, 5, 6, 7], [1, 0, 1, 0])
        for {r1, r2} <- inputs, do: Map.get(h.amps, {1 - r1, r1, 1 - r2, r2}, {0.0, 0.0})
      end
    # K[row][col]: cols are inputs; CZ = diag(1, 1, 1, −1)
    tr = Enum.with_index(cols) |> Enum.reduce({0.0, 0.0}, fn {col, i}, acc -> cadd(acc, cscale(Enum.at(col, i), if(i == 3, do: -1.0, else: 1.0))) end)
    tot = cols |> List.flatten() |> Enum.map(&abs2/1) |> Enum.sum()
    %{map: Enum.zip_with(cols, & &1) |> Enum.map(fn row -> Enum.map(row, fn {re, _im} -> Float.round(re, 9) end) end),
      process_fidelity: abs2(tr) / (4 * tot), success: tot / 4}
  end
end
