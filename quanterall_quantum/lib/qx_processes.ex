defmodule QxProcesses do
  @moduledoc """
  Quantum concepts modelled as BEAM processes -- for teaching and for demonstrating the concepts.

  ## Entanglement as linked processes

  Each round spawns three processes: a *pair* process that holds the shared two-qubit state |Φ⁺⟩, and
  Alice's and Bob's processes. Alice and Bob never message each other. Each asks the pair process to
  measure *its own* qubit along an angle in the X-Z plane; whoever asks first collapses the shared state,
  and the other then reads a result correlated with it. Over many rounds the four angle settings of the
  CHSH test give S = E(a0,b0) − E(a0,b1) + E(a1,b0) + E(a1,b1) ≈ 2√2 ≈ 2.83.

  The same processes with only a shared classical secret (a hidden angle λ, local hidden variables)
  cannot exceed S = 2: they are run beside it for comparison.

  ## Teleportation as process migration between two BEAM nodes

  Node A is this VM; node B is a second VM started with OTP's :peer (QxProcesses.NodeB). Alice's process on
  node A holds an unknown qubit |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩. The three-qubit state (ψ, Alice's
  half of a Bell pair, Bob's half) is held by a pair process on node A: a classical simulation must keep
  the entangled state somewhere. Alice entangles ψ with her half and measures: two classical bits; her
  process ends and ψ is gone from node A. What reaches node B is Bob's half of the pair (what physics
  would already have put there) and the two bits; a new process on node B applies X^m2 then Z^m1 and holds
  ψ again. The fidelity |⟨ψ|ψ_B⟩|² is computed on node B against the input sent for checking only.
  """

  @angles_a [0.0, :math.pi() / 2]
  @angles_b [:math.pi() / 4, 3 * :math.pi() / 4]

  # ── entanglement ─────────────────────────────────────────────────────────────────────────────────────

  def run_entangle(rounds \\ 400) do
    t0 = System.monotonic_time(:microsecond)
    trace_rounds = 3

    results =
      for r <- 1..rounds do
        ia = :rand.uniform(2) - 1
        ib = :rand.uniform(2) - 1
        round_quantum(Enum.at(@angles_a, ia), Enum.at(@angles_b, ib), r <= trace_rounds)
        |> Map.merge(%{a: ia, b: ib, round: r})
      end

    classical =
      for _ <- 1..rounds do
        ia = :rand.uniform(2) - 1
        ib = :rand.uniform(2) - 1
        Map.merge(round_classical(Enum.at(@angles_a, ia), Enum.at(@angles_b, ib)), %{a: ia, b: ib})
      end

    ms = (System.monotonic_time(:microsecond) - t0) / 1000

    %{
      status: "ok", name: "entangle", rounds: rounds,
      execution_time_ms: Float.round(ms, 1),
      settings: %{alice: Enum.map(@angles_a, &deg/1), bob: Enum.map(@angles_b, &deg/1)},
      quantum: chsh(results),
      classical: chsh(classical),
      trace: results |> Enum.filter(&Map.has_key?(&1, :events)) |> Enum.map(&Map.take(&1, [:round, :events, :alice, :bob])),
      # the first rounds as they came out, for an animation: a, b the setting index of each side, outcomes ±1
      pairs: results |> Enum.take(24) |> Enum.map(&Map.take(&1, [:a, :b, :alice, :bob])),
      bound: %{classical: 2.0, quantum: Float.round(2 * :math.sqrt(2), 4)},
      metadata: %{description: "Entanglement as linked processes: three BEAM processes per round (pair, Alice, Bob); CHSH test",
                  version: "1.0", engine: "BEAM processes (Elixir)", date: DateTime.utc_now() |> DateTime.to_iso8601()}
    }
  end

  defp round_quantum(theta_a, theta_b, trace?) do
    parent = self()
    t0 = System.monotonic_time(:microsecond)
    s = 1 / :math.sqrt(2)
    pair = spawn(fn -> pair_loop({s, 0.0, 0.0, s}) end)             # |Φ⁺⟩ = (|00⟩ + |11⟩)/√2, real amplitudes
    party = fn name, q, theta ->
      spawn(fn ->
        send(pair, {:measure, q, theta, self()})
        receive do
          {:outcome, o, first?} -> send(parent, {name, o, first?, self(), System.monotonic_time(:microsecond) - t0})
        end
      end)
    end
    # who asks first is left to the scheduler: either order gives the same statistics
    pa = party.(:alice, 0, theta_a)
    pb = party.(:bob, 1, theta_b)
    {oa, fa, ta} = receive do {:alice, o, f, ^pa, t} -> {o, f, t} end
    {ob, _fb, tb} = receive do {:bob, o, f, ^pb, t} -> {o, f, t} end
    send(pair, :stop)

    base = %{alice: oa, bob: ob}
    if trace? do
      first = if fa, do: "alice", else: "bob"
      Map.put(base, :events, [
        %{t_us: 0, who: "pair", what: "spawned #{inspect(pair)} holding |Φ⁺⟩"},
        %{t_us: 0, who: "alice", what: "spawned #{inspect(pa)}: measure qubit 0 at #{deg(theta_a)}°"},
        %{t_us: 0, who: "bob", what: "spawned #{inspect(pb)}: measure qubit 1 at #{deg(theta_b)}°"},
        %{t_us: min(ta, tb), who: first, what: "measured first: #{if first == "alice", do: oa, else: ob} — the shared state collapses"},
        %{t_us: max(ta, tb), who: if(first == "alice", do: "bob", else: "alice"),
          what: "measured: #{if first == "alice", do: ob, else: oa} — read from the collapsed state, no message from #{first}"}
      ])
    else
      base
    end
  end

  # the shared state: measures one qubit along an angle, collapses, answers ±1
  defp pair_loop(state, measured \\ 0) do
    receive do
      {:measure, q, theta, from} ->
        {outcome, new_state} = measure(state, q, theta)
        send(from, {:outcome, outcome, measured == 0})
        pair_loop(new_state, measured + 1)

      :stop -> :ok
    end
  end

  # rotate qubit q by Ry(−θ) so that the direction θ becomes |0⟩, then measure Z; amplitudes {a00,a01,a10,a11}
  # with index bit q (qubit 0 = least significant bit)
  defp measure(state, q, theta) do
    c = :math.cos(theta / 2)
    s = :math.sin(theta / 2)
    amps = Tuple.to_list(state)
    rot =
      for i <- 0..3 do
        bit = Bitwise.band(Bitwise.bsr(i, q), 1)
        partner = Bitwise.bxor(i, Bitwise.bsl(1, q))
        a0 = Enum.at(amps, if(bit == 0, do: i, else: partner))
        a1 = Enum.at(amps, if(bit == 0, do: partner, else: i))
        if bit == 0, do: c * a0 + s * a1, else: -s * a0 + c * a1
      end
    p0 = rot |> Enum.with_index() |> Enum.filter(fn {_, i} -> Bitwise.band(Bitwise.bsr(i, q), 1) == 0 end) |> Enum.map(fn {a, _} -> a * a end) |> Enum.sum()
    bit = if :rand.uniform() < p0, do: 0, else: 1
    norm = :math.sqrt(if bit == 0, do: p0, else: 1 - p0)
    collapsed =
      for {a, i} <- Enum.with_index(rot) do
        if Bitwise.band(Bitwise.bsr(i, q), 1) == bit, do: a / norm, else: 0.0
      end
    # back to the computational basis: Ry(+θ)
    back =
      for i <- 0..3 do
        b = Bitwise.band(Bitwise.bsr(i, q), 1)
        partner = Bitwise.bxor(i, Bitwise.bsl(1, q))
        a0 = Enum.at(collapsed, if(b == 0, do: i, else: partner))
        a1 = Enum.at(collapsed, if(b == 0, do: partner, else: i))
        if b == 0, do: c * a0 - s * a1, else: s * a0 + c * a1
      end
    {if(bit == 0, do: 1, else: -1), List.to_tuple(back)}
  end

  # local hidden variables: both share a random angle λ, each answers from it alone
  defp round_classical(theta_a, theta_b) do
    lambda = :rand.uniform() * 2 * :math.pi()
    o = fn theta -> if :math.cos(theta - lambda) >= 0, do: 1, else: -1 end
    %{alice: o.(theta_a), bob: o.(theta_b)}
  end

  defp chsh(results) do
    e =
      for ia <- 0..1, ib <- 0..1, into: %{} do
        rs = Enum.filter(results, &(&1.a == ia and &1.b == ib))
        n = length(rs)
        corr = if n == 0, do: 0.0, else: Enum.sum(Enum.map(rs, &(&1.alice * &1.bob))) / n
        {"a#{ia}b#{ib}", %{n: n, e: Float.round(corr, 4), same: Enum.count(rs, &(&1.alice == &1.bob))}}
      end
    s = e["a0b0"].e - e["a0b1"].e + e["a1b0"].e + e["a1b1"].e
    %{correlations: e, s: Float.round(s, 4)}
  end

  defp deg(rad), do: round(rad * 180 / :math.pi())

  # ── teleportation ────────────────────────────────────────────────────────────────────────────────────

  @doc "one teleportation with a random input state, plus statistics over `runs`"
  def run_teleport(runs \\ 100) do
    t0 = System.monotonic_time(:microsecond)
    node_b = QxProcesses.NodeB.info()
    first = teleport_once(true)
    rest = for _ <- 2..max(runs, 2), do: teleport_once(false)
    all = [first | rest]
    ms = (System.monotonic_time(:microsecond) - t0) / 1000
    bits = Enum.frequencies_by(all, fn r -> "#{r.m1}#{r.m2}" end)
    fids = Enum.map(all, & &1.fidelity)

    %{
      status: "ok", name: "teleport", runs: length(all),
      execution_time_ms: Float.round(ms, 1),
      node_a: %{os_pid: to_string(:os.getpid()), vm: "this node"},
      node_b: node_b,
      sample: Map.drop(first, [:fidelity_raw]),
      # the first runs without their events, for an animation
      samples: all |> Enum.take(8) |> Enum.map(&Map.take(&1, [:m1, :m2, :input, :output, :fidelity])),
      bits_counts: bits,
      fidelity_min: Float.round(Enum.min(fids), 6),
      fidelity_mean: Float.round(Enum.sum(fids) / length(fids), 6),
      metadata: %{description: "Teleportation as process migration between two BEAM nodes (node B: a second VM via :peer)",
                  version: "1.0", engine: "BEAM processes on two nodes (Elixir)", date: DateTime.utc_now() |> DateTime.to_iso8601()}
    }
  end

  defp teleport_once(trace?) do
    t0 = System.monotonic_time(:microsecond)
    theta = :math.acos(1 - 2 * :rand.uniform())                     # uniform on the Bloch sphere
    phi = :rand.uniform() * 2 * :math.pi()
    psi = {{:math.cos(theta / 2), 0.0}, {:math.sin(theta / 2) * :math.cos(phi), :math.sin(theta / 2) * :math.sin(phi)}}

    # node A: Alice's process with ψ, and the pair process with the three-qubit state
    parent = self()
    alice =
      spawn(fn ->
        {m1, m2, bob_half} = alice_measures(psi)
        send(parent, {:alice_done, self(), m1, m2, bob_half})          # ψ no longer exists on node A
      end)

    {m1, m2, bob_half} = receive do {:alice_done, ^alice, a, b, h} -> {a, b, h} end
    t_alice = System.monotonic_time(:microsecond) - t0

    # node B: a new process, given only Bob's half of the pair and the two classical bits
    {:ok, result} = QxProcesses.NodeB.call(QxProcesses, :bob_receives, [bob_half, m1, m2, psi])
    t_bob = System.monotonic_time(:microsecond) - t0

    base = %{m1: m1, m2: m2, fidelity: Float.round(result.fidelity, 9),
             input: %{theta_deg: Float.round(theta * 180 / :math.pi(), 2), phi_deg: Float.round(phi * 180 / :math.pi(), 2)},
             output: result.bloch}

    if trace? do
      Map.put(base, :events, [
        %{t_us: 0, where: "A", what: "Alice's process #{inspect(alice)} holds |ψ⟩ (θ = #{base.input.theta_deg}°, φ = #{base.input.phi_deg}°)"},
        %{t_us: 0, where: "A", what: "Bell pair made; Bob's half belongs to node B"},
        %{t_us: t_alice, where: "A", what: "Alice entangles ψ with her half and measures: bits #{m1}#{m2}; her process ends, ψ is gone from node A"},
        %{t_us: t_alice, where: "A→B", what: "classical channel: 2 bits (#{m1}, #{m2})"},
        %{t_us: t_bob, where: "B", what: "process #{result.pid} on node B (OS pid #{result.os_pid}) applies #{corr(m1, m2)} and holds ψ again: fidelity #{Float.round(result.fidelity, 6)}"}
      ])
    else
      base
    end
  end

  defp corr(m1, m2) do
    case {m2, m1} do
      {0, 0} -> "no correction"
      {1, 0} -> "X"
      {0, 1} -> "Z"
      {1, 1} -> "X then Z"
    end
  end

  # the three-qubit state: q0 = ψ, q1 = Alice's half, q2 = Bob's half; complex amplitudes {re, im}
  defp alice_measures({a, b}) do
    z = {0.0, 0.0}
    st = List.to_tuple([a, b, z, z, z, z, z, z])                     # ψ ⊗ |0⟩ ⊗ |0⟩ (q0 the low bit)
    st = st |> h(1) |> cx(1, 2) |> cx(0, 1) |> h(0)
    {m1, st} = measure_z(st, 0)
    {m2, st} = measure_z(st, 1)
    base = m1 + 2 * m2
    bob_half = {elem(st, base), elem(st, base + 4)}                  # q2 = 0, 1 with q0 = m1, q1 = m2
    {m1, m2, normalise(bob_half)}
  end

  @doc false
  # runs on node B
  def bob_receives({b0, b1}, m1, m2, {p0, p1}) do
    {b0, b1} = if m2 == 1, do: {b1, b0}, else: {b0, b1}             # X
    {b0, b1} = if m1 == 1, do: {b0, cmul({-1.0, 0.0}, b1)}, else: {b0, b1}   # Z
    # |⟨ψ|ψ_B⟩|²: the input is passed only to check the result, not used to build it
    ov = cadd(cmul(conj(p0), b0), cmul(conj(p1), b1))
    fid = cabs2(ov)
    %{pid: inspect(self()), os_pid: to_string(:os.getpid()), fidelity: fid, bloch: bloch({b0, b1})}
  end

  defp bloch({a, b}) do
    # global phase removed: a real and non-negative
    pa = :math.atan2(elem(a, 1), elem(a, 0))
    rot = {:math.cos(-pa), :math.sin(-pa)}
    a2 = cmul(a, rot)
    b2 = cmul(b, rot)
    theta = 2 * :math.acos(min(1.0, max(-1.0, elem(a2, 0))))
    phi = :math.atan2(elem(b2, 1), elem(b2, 0))
    phi = if phi < 0, do: phi + 2 * :math.pi(), else: phi
    %{theta_deg: Float.round(theta * 180 / :math.pi(), 2), phi_deg: Float.round(phi * 180 / :math.pi(), 2)}
  end

  # ── small complex state-vector helpers (index bit q = qubit q) ──
  defp cadd({a, b}, {c, d}), do: {a + c, b + d}
  defp cmul({a, b}, {c, d}), do: {a * c - b * d, a * d + b * c}
  defp conj({a, b}), do: {a, -b}
  defp cabs2({a, b}), do: a * a + b * b
  defp cscale({a, b}, k), do: {a * k, b * k}

  defp normalise({x, y}) do
    n = :math.sqrt(cabs2(x) + cabs2(y))
    {cscale(x, 1 / n), cscale(y, 1 / n)}
  end

  defp h(st, q) do
    k = 1 / :math.sqrt(2)
    n = tuple_size(st)
    for i <- 0..(n - 1) do
      j = Bitwise.bxor(i, Bitwise.bsl(1, q))
      if Bitwise.band(Bitwise.bsr(i, q), 1) == 0,
        do: cscale(cadd(elem(st, i), elem(st, j)), k),
        else: cscale(cadd(elem(st, j), cscale(elem(st, i), -1.0)), k)
    end
    |> List.to_tuple()
  end

  defp cx(st, c, t) do
    n = tuple_size(st)
    for i <- 0..(n - 1) do
      if Bitwise.band(Bitwise.bsr(i, c), 1) == 1, do: elem(st, Bitwise.bxor(i, Bitwise.bsl(1, t))), else: elem(st, i)
    end
    |> List.to_tuple()
  end

  defp measure_z(st, q) do
    n = tuple_size(st)
    p1 = Enum.sum(for i <- 0..(n - 1), Bitwise.band(Bitwise.bsr(i, q), 1) == 1, do: cabs2(elem(st, i)))
    bit = if :rand.uniform() < p1, do: 1, else: 0
    norm = :math.sqrt(if bit == 1, do: p1, else: 1 - p1)
    st =
      for i <- 0..(n - 1) do
        if Bitwise.band(Bitwise.bsr(i, q), 1) == bit, do: cscale(elem(st, i), 1 / norm), else: {0.0, 0.0}
      end
      |> List.to_tuple()
    {bit, st}
  end
end

defmodule QxProcesses.NodeB do
  @moduledoc """
  Node B for the teleportation model: a second BEAM node (a separate OS process) started with OTP's :peer
  over standard I/O -- no Erlang distribution, so no network port and no cookie. Restarted when it dies.
  """
  use GenServer

  def start_link(_), do: GenServer.start_link(__MODULE__, nil, name: __MODULE__)

  def call(m, f, a), do: GenServer.call(__MODULE__, {:call, m, f, a}, 15_000)
  def info, do: GenServer.call(__MODULE__, :info, 15_000)

  @impl true
  def init(_), do: {:ok, %{peer: nil, os_pid: nil}}

  @impl true
  def handle_call({:call, m, f, a}, _from, st) do
    st = ensure(st)
    {:reply, {:ok, :peer.call(st.peer, m, f, a, 10_000)}, st}
  end

  def handle_call(:info, _from, st) do
    st = ensure(st)
    {:reply, %{os_pid: st.os_pid, vm: "second BEAM node (OTP :peer)"}, st}
  end

  @impl true
  def handle_info({:EXIT, _, _}, st), do: {:noreply, %{st | peer: nil}}
  def handle_info(_, st), do: {:noreply, st}

  defp ensure(%{peer: p} = st) when is_pid(p) do
    if Process.alive?(p), do: st, else: ensure(%{st | peer: nil})
  end

  defp ensure(st) do
    Process.flag(:trap_exit, true)
    paths = Enum.flat_map(:code.get_path(), fn p -> [~c"-pa", p] end)
    peer =
      case :peer.start_link(%{connection: :standard_io, args: paths}) do
        {:ok, p} -> p
        {:ok, p, _node} -> p
      end
    os_pid = to_string(:peer.call(peer, :os, :getpid, []))
    %{st | peer: peer, os_pid: os_pid}
  end
end
