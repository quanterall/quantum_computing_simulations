defmodule QxStab do
  @moduledoc """
  A stabilizer (Clifford) simulator with gate-level noise, on the BEAM.

  Two parts, the method Stim uses:

    * a **tableau** (Aaronson–Gottesman, "CHP") runs the circuit once without noise and gives a reference
      sample of every measurement — exact, for circuits of hundreds of qubits;
    * a **Pauli-frame sampler** then runs the shots: each shot carries only the difference from the
      reference (a Pauli frame), pushed through the gates. Frames are bit-packed across shots: qubit q's
      X part is one integer whose bit k belongs to shot k, so a CNOT on 1,024 shots is two XORs of 1,024-bit
      integers. Batches of shots run in parallel, one BEAM process per batch.

  A circuit is a list of instructions, the same as Stim's text format (`to_stim/1` writes it out, so the
  exact circuit can be checked with Stim):

      {:r, qs}  {:h, qs}  {:s, qs}  {:cx, [{c, t}, ...]}  {:m, qs}            gates and measurements
      {:x_error, qs, p}  {:depolarize1, qs, p}  {:depolarize2, pairs, p}       noise channels
      {:heralded_erase, qs, p}   with probability p a uniformly random Pauli (I, X, Y, Z), and a herald bit
                                 (1 when it happened) appended to the measurement record — Stim's HERALDED_ERASE
      {:detector, [rec, ...]}  {:observable, [rec, ...]}  :tick                 annotations

  `rec` is an absolute measurement index (0 = the first measurement of the circuit).
  """
  import Bitwise

  # ── the tableau ────────────────────────────────────────────────────────────────────────────────────────
  # rows 0..n-1 are destabilizers, n..2n-1 stabilizers; a row is {x, z, r}: X and Z parts as bitmasks
  # over the qubits, r the sign bit. Starts as |0…0⟩: destabilizer i = X_i, stabilizer i = Z_i.

  defp tableau(n) do
    rows = (for i <- 0..(n - 1), do: {1 <<< i, 0, 0}) ++ (for i <- 0..(n - 1), do: {0, 1 <<< i, 0})
    %{n: n, rows: List.to_tuple(rows)}
  end

  defp map_rows(t, f), do: %{t | rows: t.rows |> Tuple.to_list() |> Enum.map(f) |> List.to_tuple()}

  defp b(v, q), do: (v >>> q) &&& 1

  defp t_h(t, q) do
    map_rows(t, fn {x, z, r} ->
      xq = b(x, q); zq = b(z, q)
      {x |> bxor((xq ^^^ zq) <<< q), z |> bxor((xq ^^^ zq) <<< q), r ^^^ (xq &&& zq)}
    end)
  end

  defp t_s(t, q) do
    map_rows(t, fn {x, z, r} ->
      xq = b(x, q); zq = b(z, q)
      {x, bxor(z, xq <<< q), r ^^^ (xq &&& zq)}
    end)
  end

  defp t_cx(t, c, tq) do
    map_rows(t, fn {x, z, r} ->
      xc = b(x, c); zc = b(z, c); xt = b(x, tq); zt = b(z, tq)
      {bxor(x, xc <<< tq), bxor(z, zt <<< c), r ^^^ (xc &&& zt &&& (xt ^^^ zc ^^^ 1))}
    end)
  end

  # the sign of the product of two Pauli rows (Aaronson–Gottesman's g, summed with bit counts)
  defp rowsum({x1, z1, r1}, {x2, z2, r2}) do
    y1 = x1 &&& z1; xo = x1 &&& bnot(z1); zo = bnot(x1) &&& z1
    plus = pop(y1 &&& z2 &&& bnot(x2)) + pop(xo &&& z2 &&& x2) + pop(zo &&& x2 &&& bnot(z2))
    minus = pop(y1 &&& x2 &&& bnot(z2)) + pop(xo &&& z2 &&& bnot(x2)) + pop(zo &&& x2 &&& z2)
    phase = Integer.mod(2 * r1 + 2 * r2 + plus - minus, 4)
    {bxor(x1, x2), bxor(z1, z2), if(phase == 2, do: 1, else: 0)}
  end

  defp pop(v) when v <= 0, do: 0
  defp pop(v), do: for(<<bit::1 <- :binary.encode_unsigned(v)>>, reduce: 0, do: (acc -> acc + bit))

  # measure qubit q in Z; a random outcome is taken as `choice` (the reference sample uses 0)
  defp t_measure(t, q, choice) do
    n = t.n
    p = Enum.find(n..(2 * n - 1), fn i -> b(elem(elem(t.rows, i), 0), q) == 1 end)
    if p do
      rp = elem(t.rows, p)
      rows =
        for i <- 0..(2 * n - 1) do
          row = elem(t.rows, i)
          cond do
            i == p -> {0, 1 <<< q, choice}
            i == p - n -> rp
            b(elem(row, 0), q) == 1 -> rowsum(row, rp)
            true -> row
          end
        end
      {choice, true, %{t | rows: List.to_tuple(rows)}}
    else
      scratch =
        Enum.reduce(0..(n - 1), {0, 0, 0}, fn i, acc ->
          if b(elem(elem(t.rows, i), 0), q) == 1, do: rowsum(acc, elem(t.rows, i + n)), else: acc
        end)
      {elem(scratch, 2), false, t}
    end
  end

  defp t_reset(t, q) do
    {bit, _random, t} = t_measure(t, q, 0)
    if bit == 1, do: t |> t_h(q) |> t_s(q) |> t_s(q) |> t_h(q), else: t                    # X = H S S H
  end

  @doc """
  The noiseless reference: every measurement's outcome, random ones taken as 0. Returns a tuple of bits.
  """
  def reference(circuit, n) do
    {_t, rec} =
      Enum.reduce(circuit, {tableau(n), []}, fn
        {:r, qs}, {t, rec} -> {Enum.reduce(qs, t, &t_reset(&2, &1)), rec}
        {:h, qs}, {t, rec} -> {Enum.reduce(qs, t, &t_h(&2, &1)), rec}
        {:s, qs}, {t, rec} -> {Enum.reduce(qs, t, &t_s(&2, &1)), rec}
        {:cx, pairs}, {t, rec} -> {Enum.reduce(pairs, t, fn {c, x}, acc -> t_cx(acc, c, x) end), rec}
        {:m, qs}, {t, rec} ->
          Enum.reduce(qs, {t, rec}, fn q, {t, rec} -> {bit, _, t} = t_measure(t, q, 0); {t, [bit | rec]} end)
        {:heralded_erase, qs, _p}, {t, rec} -> {t, List.duplicate(0, length(qs)) ++ rec}
        _, acc -> acc
      end)
    rec |> Enum.reverse() |> List.to_tuple()
  end

  # ── the Pauli-frame sampler ────────────────────────────────────────────────────────────────────────────

  @doc """
  Samples `shots` runs of the noisy circuit. Returns `%{detectors: [integer], observables: [integer], shots}`:
  detector j's integer has bit k set when shot k's detection event j fired (likewise the observables).
  Runs in batches of `batch` shots, in parallel; `seed` makes it repeatable.
  """
  def sample(circuit, n, shots, opts \\ []) do
    batch = Keyword.get(opts, :batch, 1024)
    seed = Keyword.get(opts, :seed, :rand.uniform(1 <<< 32) - 1)
    ref = reference(circuit, n)
    nb = div(shots + batch - 1, batch)
    results =
      0..(nb - 1)
      |> Task.async_stream(
        fn i ->
          size = min(batch, shots - i * batch)
          :rand.seed(:exsss, {seed, i, 7})
          frame_batch(circuit, n, size, ref)
        end,
        max_concurrency: System.schedulers_online(), timeout: :infinity, ordered: true)
      |> Enum.map(fn {:ok, r} -> r end)
    # stitch the batches: shot k of batch i becomes bit (i*batch + k)
    stitch = fn key ->
      cols = results |> Enum.map(&List.to_tuple(Map.fetch!(&1, key))) |> Enum.with_index()
      count = tuple_size(elem(hd(cols), 0))
      for j <- 0..(count - 1)//1 do
        Enum.reduce(cols, 0, fn {c, i}, acc -> acc ||| (elem(c, j) <<< (i * batch)) end)
      end
    end
    %{detectors: (if results == [], do: [], else: stitch.(:detectors)),
      observables: (if results == [], do: [], else: stitch.(:observables)), shots: shots, seed: seed}
  end

  defp frame_batch(circuit, n, size, ref) do
    full = (1 <<< size) - 1
    rnd = fn -> rand_bits(size) end
    x = :array.new(n, default: 0)
    z = Enum.reduce(0..(n - 1), :array.new(n, default: 0), fn q, a -> :array.set(q, rnd.(), a) end)
    st = %{x: x, z: z, meas: [], size: size, full: full}
    st =
      Enum.reduce(circuit, st, fn
        {:r, qs}, st -> Enum.reduce(qs, st, fn q, st -> %{st | x: :array.set(q, 0, st.x), z: :array.set(q, rnd.(), st.z)} end)
        {:h, qs}, st ->
          Enum.reduce(qs, st, fn q, st -> xq = :array.get(q, st.x); zq = :array.get(q, st.z)
            %{st | x: :array.set(q, zq, st.x), z: :array.set(q, xq, st.z)} end)
        {:s, qs}, st ->
          Enum.reduce(qs, st, fn q, st -> %{st | z: :array.set(q, bxor(:array.get(q, st.z), :array.get(q, st.x)), st.z)} end)
        {:cx, pairs}, st ->
          Enum.reduce(pairs, st, fn {c, t}, st ->
            %{st | x: :array.set(t, bxor(:array.get(t, st.x), :array.get(c, st.x)), st.x),
                   z: :array.set(c, bxor(:array.get(c, st.z), :array.get(t, st.z)), st.z)} end)
        {:m, qs}, st ->
          # a shot's outcome is the reference outcome flipped where its frame has an X here; keep the flips
          Enum.reduce(qs, st, fn q, st ->
            %{st | meas: [:array.get(q, st.x) | st.meas], z: :array.set(q, rnd.(), st.z)}   # the measurement collapses Z
          end)
        {:x_error, qs, p}, st -> Enum.reduce(qs, st, fn q, st -> %{st | x: :array.set(q, bxor(:array.get(q, st.x), bern(size, p)), st.x)} end)
        {:depolarize1, qs, p}, st ->
          Enum.reduce(qs, st, fn q, st ->
            {ex, ez} = pick(size, p, 3, fn r -> {r &&& 1, r >>> 1} end)          # 1 = X, 2 = Z, 3 = Y
            %{st | x: :array.set(q, bxor(:array.get(q, st.x), ex), st.x), z: :array.set(q, bxor(:array.get(q, st.z), ez), st.z)}
          end)
        {:depolarize2, pairs, p}, st ->
          Enum.reduce(pairs, st, fn {a, c}, st ->
            {ax, az, cx, cz} = pick4(size, p)
            %{st | x: st.x |> then(&:array.set(a, bxor(:array.get(a, &1), ax), &1)) |> then(&:array.set(c, bxor(:array.get(c, &1), cx), &1)),
                   z: st.z |> then(&:array.set(a, bxor(:array.get(a, &1), az), &1)) |> then(&:array.set(c, bxor(:array.get(c, &1), cz), &1))}
          end)
        {:heralded_erase, qs, p}, st ->
          Enum.reduce(qs, st, fn q, st ->
            hit = bern(size, p)
            {rx, rz} = {rand_bits(size) &&& hit, rand_bits(size) &&& hit}             # I, X, Z, Y each ¼
            %{st | x: :array.set(q, bxor(:array.get(q, st.x), rx), st.x), z: :array.set(q, bxor(:array.get(q, st.z), rz), st.z),
                   meas: [hit | st.meas]}                                              # the herald fires where it happened
          end)
        _, st -> st
      end)
    flips = st.meas |> Enum.reverse() |> List.to_tuple()
    # a detector (or the observable) fires where the parity of its measurements differs from the reference
    # parity, that is, where an odd number of them flipped
    fired = fn recs -> Enum.reduce(recs, 0, fn k, acc -> bxor(acc, elem(flips, k)) end) end
    _ = {ref, full}
    %{detectors: for({:detector, recs} <- circuit, do: fired.(recs)),
      observables: for({:observable, recs} <- circuit, do: fired.(recs))}
  end

  defp rand_bits(0), do: 0
  defp rand_bits(size) do
    bytes = div(size + 7, 8)
    :binary.decode_unsigned(:rand.bytes(bytes)) &&& ((1 <<< size) - 1)
  end

  # positions of a Bernoulli(p) mask over `size` shots, by geometric gaps
  defp positions(_size, p) when p <= 0, do: []
  defp positions(size, p) when p >= 1, do: Enum.to_list(0..(size - 1))
  defp positions(size, p) do
    lq = :math.log(1 - p)
    Stream.unfold(-1, fn i ->
      j = i + 1 + floor(:math.log(1.0 - :rand.uniform_real()) / lq)
      if j < size, do: {j, j}, else: nil
    end)
    |> Enum.to_list()
  end

  defp bern(size, p), do: Enum.reduce(positions(size, p), 0, fn j, acc -> acc ||| (1 <<< j) end)

  # each hit picks one of k non-identity Paulis; f maps 1..k to {x bit, z bit}
  defp pick(size, p, k, f) do
    Enum.reduce(positions(size, p), {0, 0}, fn j, {ex, ez} ->
      {xb, zb} = f.(:rand.uniform(k))
      {ex ||| (xb <<< j), ez ||| (zb <<< j)}
    end)
  end

  # two-qubit depolarizing: one of the 15 non-identity Pauli pairs
  defp pick4(size, p) do
    Enum.reduce(positions(size, p), {0, 0, 0, 0}, fn j, {ax, az, cx, cz} ->
      r = :rand.uniform(15)
      {ax ||| ((r &&& 1) <<< j), az ||| (((r >>> 1) &&& 1) <<< j), cx ||| (((r >>> 2) &&& 1) <<< j), cz ||| (((r >>> 3) &&& 1) <<< j)}
    end)
  end

  # ── circuits ───────────────────────────────────────────────────────────────────────────────────────────

  @doc "GHZ on n qubits: H, a CNOT chain, measure all. The outcomes are all 0 or all 1."
  def ghz(n), do: [{:r, Enum.to_list(0..(n - 1))}, {:h, [0]}, {:cx, for(i <- 0..(n - 2)//1, do: {i, i + 1})}, {:m, Enum.to_list(0..(n - 1))}]

  @doc """
  A repetition-code memory: d data qubits (even indices), d−1 ancillas between them, `rounds` rounds of ZZ
  checks, then the data measured. Circuit-level noise p: after resets (flip), after each gate
  (depolarizing), before measurements (flip), on the data at the start of each round (depolarizing).
  """
  def repetition_memory(d, rounds, p) do
    data = for i <- 0..(d - 1), do: 2 * i
    anc = for i <- 0..(d - 2), do: 2 * i + 1
    n = 2 * d - 1
    round = fn ->
      [{:depolarize1, data, p}, {:r, anc}, {:x_error, anc, p}, :tick,
       {:cx, for(a <- anc, do: {a - 1, a})}, {:depolarize2, for(a <- anc, do: {a - 1, a}), p}, :tick,
       {:cx, for(a <- anc, do: {a + 1, a})}, {:depolarize2, for(a <- anc, do: {a + 1, a}), p}, :tick,
       {:x_error, anc, p}, {:m, anc}]
    end
    memory(n, data, Enum.map(anc, &{:z, &1, [&1 - 1, &1 + 1]}), rounds, p, round, [0])
  end

  @doc """
  The rotated surface code of distance d, memory in Z: d² data qubits, d² − 1 measure qubits (X and Z
  checks), `rounds` rounds, circuit-level noise p as in `repetition_memory/3`. Data at odd grid points,
  measure qubits at even ones; the CNOT order avoids hook errors along the logical operator.
  """
  def surface_memory(d, rounds, p, pe \\ 0.0, herald_dets \\ false) do
    dpos = for y <- 0..(d - 1), x <- 0..(d - 1), do: {2 * x + 1, 2 * y + 1}
    # a measure qubit at every even grid point with 4 data neighbours, and on the boundary with 2: X checks on
    # the top and bottom edges, Z checks on the left and right
    mpos =
      (for y <- 0..d, x <- 0..d, do: {2 * x, 2 * y, if(rem(x + y, 2) == 0, do: :x, else: :z)})
      |> Enum.filter(fn {mx, my, kind} ->
        nb = Enum.count([{mx + 1, my + 1}, {mx + 1, my - 1}, {mx - 1, my + 1}, {mx - 1, my - 1}], &(&1 in dpos))
        nb == 4 or (nb == 2 and ((kind == :x and my in [0, 2 * d]) or (kind == :z and mx in [0, 2 * d])))
      end)
      |> Enum.map(fn {mx, my, kind} -> {{mx, my}, kind} end)
    idx = (dpos ++ Enum.map(mpos, &elem(&1, 0))) |> Enum.with_index() |> Map.new()
    data = Enum.map(dpos, &idx[&1])
    # neighbour order per layer. A fault halfway through a check spreads to its last two data qubits (a hook);
    # X checks end on a horizontal pair and Z checks on a vertical one, across the logical errors' direction,
    # so a hook never counts as two steps of a logical error (Stim's distance check confirms d)
    order = %{x: [{1, 1}, {-1, 1}, {1, -1}, {-1, -1}], z: [{1, 1}, {1, -1}, {-1, 1}, {-1, -1}]}
    checks = for {{mx, my} = m, kind} <- mpos do
      {kind, idx[m], order[kind] |> Enum.map(fn {ox, oy} -> {mx + ox, my + oy} end) |> Enum.filter(&(&1 in dpos)) |> Enum.map(&idx[&1])}
    end
    # CNOT layer k: every check with a data neighbour at its k-th offset (X checks: ancilla controls)
    layer = fn k ->
      Enum.flat_map(mpos, fn {{mx, my} = m, kind} ->
        {ox, oy} = Enum.at(order[kind], k)
        q = {mx + ox, my + oy}
        cond do
          q not in dpos -> []
          kind == :x -> [{idx[m], idx[q]}]
          true -> [{idx[q], idx[m]}]
        end
      end)
    end
    xs = for {:x, a, _} <- checks, do: a
    anc = Enum.map(checks, &elem(&1, 1))
    round = fn ->
      (if pe > 0, do: [{:heralded_erase, data, pe}], else: []) ++
      [{:depolarize1, data, p}, {:r, anc}, {:x_error, anc, p}, :tick, {:h, xs}, {:depolarize1, xs, p}, :tick] ++
        Enum.flat_map(0..3, fn k -> [{:cx, layer.(k)}, {:depolarize2, layer.(k), p}, :tick] end) ++
        [{:h, xs}, {:depolarize1, xs, p}, :tick, {:x_error, anc, p}, {:m, anc}]
    end
    # the logical Z: a row of data qubits that meets every X check an even number of times
    obs =
      Enum.find_value(0..(d - 1), fn y ->
        row = for x <- 0..(d - 1), do: idx[{2 * x + 1, 2 * y + 1}]
        if Enum.all?(for({:x, _, nb} <- checks, do: rem(Enum.count(nb, &(&1 in row)), 2) == 0)), do: row
      end) ||
      Enum.find_value(0..(d - 1), fn x ->
        col = for y <- 0..(d - 1), do: idx[{2 * x + 1, 2 * y + 1}]
        if Enum.all?(for({:x, _, nb} <- checks, do: rem(Enum.count(nb, &(&1 in col)), 2) == 0)), do: col
      end)
    memory(map_size(idx), data, checks, rounds, p, round, obs, (if pe > 0, do: length(data), else: 0), herald_dets)
  end

  # rounds of checks, detectors comparing each check with its previous round (Z checks also with the
  # initial |0⟩ state and with the final data measurement), the logical observable from the final data
  # nh: herald bits written at the start of each round (heralded erasure on the data), before the ancillas'
  # measurements; herald_dets: declare each herald a detector (for comparing samplers, not for decoding)
  defp memory(n, data, checks, rounds, p, round, obs_qubits, nh \\ 0, herald_dets \\ false) do
    anc = Enum.map(checks, &elem(&1, 1))
    na = length(anc)
    per = na + nh
    init = [{:r, Enum.to_list(0..(n - 1))}, {:x_error, Enum.to_list(0..(n - 1)), p}, :tick]
    {body, _} =
      Enum.map_reduce(1..rounds, 0, fn r, base ->
        heralds = if herald_dets, do: for(j <- 0..(nh - 1)//1, do: {:detector, [base + j]}), else: []
        dets = for {{kind, _, _}, i} <- Enum.with_index(checks), kind == :z or r > 1 do
          now = base + nh + i
          if r == 1, do: {:detector, [now]}, else: {:detector, [now, now - per]}
        end
        {round.() ++ heralds ++ dets, base + per}
      end)
    last = (rounds - 1) * per + nh
    mdata = rounds * per
    pos = data |> Enum.with_index() |> Map.new()
    final =
      [{:x_error, data, p}, {:m, data}] ++
        (for {{:z, _, nb}, i} <- Enum.with_index(checks), do: {:detector, [last + i | Enum.map(nb, &(mdata + pos[&1]))]}) ++
        [{:observable, Enum.map(obs_qubits, &(mdata + pos[&1]))}]
    {init ++ List.flatten(body) ++ final, n}
  end

  @doc """
  Writes a sample as bytes: each detector, then each observable, as `ceil(shots/8)` little-endian bytes,
  bit k = shot k (a numpy `unpackbits(bitorder="little")` reads it back).
  """
  def write_packed(path, %{detectors: dets, observables: obs, shots: shots}) do
    nbytes = div(shots + 7, 8)
    File.write!(path, for(v <- dets ++ obs, into: <<>>, do: <<v::little-size(nbytes * 8)>>))
  end

  # ── the demos behind the public endpoints ───────────────────────────────────────────────────────────

  @doc "GHZ on n qubits, sampled: every shot should be all 0 or all 1 (detectors compare neighbours)."
  def run_ghz(n, shots) do
    t0 = System.monotonic_time(:microsecond)
    circuit = ghz(n) ++ (for i <- 0..(n - 2)//1, do: {:detector, [i, i + 1]}) ++ [{:observable, [0]}]
    r = sample(circuit, n, shots)
    bad = Enum.reduce(r.detectors, 0, &bor/2)
    ones = hd(r.observables) &&& bnot(bad)
    nbad = pop(bad)
    n1 = pop(ones)
    %{status: "ok", name: "stabilizer_ghz", qubits: n, shots: shots,
      execution_time_ms: Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3),
      counts: %{"all 0" => shots - nbad - n1, "all 1" => n1, "other" => nbad}, seed: r.seed,
      metadata: %{description: "GHZ state on #{n} qubits on the stabilizer simulator (a state vector would need 2^#{n} amplitudes)",
                  engine: "QxStab: tableau + bit-packed Pauli frames (Elixir)"}}
  end

  @doc """
  A memory experiment under circuit-level noise p: the repetition code or the rotated surface code, distance
  d, `rounds` rounds. Returns how often the detectors fire (per round and overall), how many shots saw any,
  and how often the logical observable flipped before any decoding; the circuit in Stim's format for d ≤ 5.
  """
  def run_memory(code, d, rounds, p, shots) do
    t0 = System.monotonic_time(:microsecond)
    {circuit, n} = if code == :surface, do: surface_memory(d, rounds, p), else: repetition_memory(d, rounds, p)
    r = sample(circuit, n, shots)
    ndet = length(r.detectors)
    any = Enum.reduce(r.detectors, 0, &bor/2)
    rate = if ndet == 0, do: 0.0, else: Enum.sum(Enum.map(r.detectors, &pop/1)) / (ndet * shots)
    %{status: "ok", name: "stabilizer_memory", code: Atom.to_string(code), distance: d, rounds: rounds, p: p,
      qubits: n, detectors: ndet, shots: shots, seed: r.seed,
      execution_time_ms: Float.round((System.monotonic_time(:microsecond) - t0) / 1000, 3),
      detection_rate: Float.round(rate, 6),
      shots_with_detections: pop(any),
      observable_flip_rate_undecoded: Float.round(pop(hd(r.observables)) / shots, 6),
      stim_circuit: if(d <= 5, do: to_stim(circuit), else: nil),
      metadata: %{description: "#{code} code memory, d = #{d}, #{rounds} rounds, circuit-level noise p = #{p}; decoding is not done here: the published results decode with PyMatching",
                  noise: "after resets and before measurements: X flip p; after each gate: depolarizing p (1- and 2-qubit); data at each round's start: depolarizing p",
                  engine: "QxStab: tableau + bit-packed Pauli frames (Elixir), batches in parallel"}}
  end

  # ── Stim's text format ─────────────────────────────────────────────────────────────────────────────────

  @doc "The circuit in Stim's text format (measurement references made relative), to check it with Stim."
  def to_stim(circuit) do
    {lines, _} =
      Enum.map_reduce(circuit, 0, fn
        {:m, qs}, k -> {"M " <> js(qs), k + length(qs)}
        {:r, qs}, k -> {"R " <> js(qs), k}
        {:h, qs}, k -> {"H " <> js(qs), k}
        {:s, qs}, k -> {"S " <> js(qs), k}
        {:cx, ps}, k -> {"CX " <> js(Enum.flat_map(ps, &Tuple.to_list/1)), k}
        {:x_error, qs, p}, k -> {"X_ERROR(#{p}) " <> js(qs), k}
        {:depolarize1, qs, p}, k -> {"DEPOLARIZE1(#{p}) " <> js(qs), k}
        {:depolarize2, ps, p}, k -> {"DEPOLARIZE2(#{p}) " <> js(Enum.flat_map(ps, &Tuple.to_list/1)), k}
        {:heralded_erase, qs, p}, k -> {"HERALDED_ERASE(#{p}) " <> js(qs), k + length(qs)}
        {:detector, recs}, k -> {"DETECTOR " <> Enum.map_join(recs, " ", &"rec[#{&1 - k}]"), k}
        {:observable, recs}, k -> {"OBSERVABLE_INCLUDE(0) " <> Enum.map_join(recs, " ", &"rec[#{&1 - k}]"), k}
        :tick, k -> {"TICK", k}
      end)
    lines |> Enum.reject(&String.ends_with?(&1, " ")) |> Enum.join("\n")
  end

  defp js(qs), do: Enum.join(qs, " ")
end
