# The BEAM side of the benchmarks and the noise study (tools/stab_study.py, tools/sv_bench.py run it):
#   docker cp tools/beam_runs.exs quanterall_quantum:/tmp/ && docker exec -w /app quanterall_quantum \
#     sh -c 'MIX_ENV=prod mix run --no-start /tmp/beam_runs.exs <job> <out dir>'
# job "study": surface-code memories under circuit-level noise, sampled; job "sv": state-vector timings.
[job, out] = System.argv()
File.mkdir_p!(out)

case job do
  "study" ->
    plan = for d <- [3, 5, 7], p <- [0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012], do: {d, p}
    manifest =
      for {d, p} <- plan do
        shots = 100_000
        {c, n} = QxStab.surface_memory(d, d, p)
        name = "surface_d#{d}_p#{p}"
        File.write!(Path.join(out, name <> ".stim"), QxStab.to_stim(c))
        {us, r} = :timer.tc(fn -> QxStab.sample(c, n, shots, seed: 20261029 + d * 1000 + round(p * 1.0e5)) end)
        QxStab.write_packed(Path.join(out, name <> ".bin"), r)
        IO.puts("#{name}: #{shots} shots in #{div(us, 1000)} ms")
        %{name: name, d: d, p: p, rounds: d, shots: shots, qubits: n, beam_ms: us / 1000, seed: r.seed}
      end
    File.write!(Path.join(out, "study.json"), JSON.encode!(manifest))

  "sv" ->
    # random layers: Rx on every qubit (angles fixed), then CNOTs on even pairs, then odd; n layers; 1,024
    # shots. Timed: building the state and sampling, inside the VM (no HTTP); one warm-up, then 5 runs (3 from
    # 17 qubits on), median and minimum. Stops after a run over 60 s; 20 qubits exceed the 1 GB container.
    rows =
      for n <- 4..19, reduce: [] do
        acc ->
          last = List.first(acc)
          if last && last.median_ms > 60_000 do
            acc
          else
            circ = fn s ->
              Enum.reduce(0..(n - 1), s, fn l, s ->
                s = Enum.reduce(0..(n - 1), s, fn q, s -> QxSim.rx(s, n, q, 0.1 * (q + 1) + 0.01 * l) end)
                s = Enum.reduce(0..(n - 2)//2, s, fn q, s -> QxSim.cx(s, n, q, q + 1) end)
                Enum.reduce(1..(n - 2)//2, s, fn q, s -> QxSim.cx(s, n, q, q + 1) end)
              end)
            end
            reps = if n >= 17, do: 3, else: 5
            QxSim.run(n, 1024, circ)
            ts = for _ <- 1..reps, do: elem(:timer.tc(fn -> QxSim.run(n, 1024, circ) end), 0) / 1000
            med = Enum.at(Enum.sort(ts), div(reps, 2))
            IO.puts("sv n=#{n}: median #{round(med)} ms")
            [%{n: n, median_ms: med, min_ms: Enum.min(ts), runs: reps} | acc]
          end
      end
    File.write!(Path.join(out, "sv.json"), JSON.encode!(Enum.reverse(rows)))

  "xval" ->
    # random circuits for the cross-check against Qiskit Aer (tools/stab_study.py): for n = 2..12, 5 circuits
    # of 6n gates drawn from h, x, rx, ry, rz, cx, cz, cp, swap with angles; the final state vector, written out
    :rand.seed(:exsss, {2026, 10, 30})
    cases =
      for n <- 2..12, c <- 1..5 do
        gates =
          for _ <- 1..(6 * n) do
            g = Enum.random(~w(h x rx ry rz cx cz cp swap))
            [a, b] = Enum.take_random(0..(n - 1), 2)
            t = Float.round(:rand.uniform() * 2 * :math.pi(), 6)
            case g do
              g when g in ~w(h x) -> [g, a]
              g when g in ~w(rx ry rz) -> [g, a, t]
              "cp" -> [g, a, b, t]
              _ -> [g, a, b]
            end
          end
        state =
          Enum.reduce(gates, QxSim.init(n), fn
            ["h", a], s -> QxSim.h(s, n, a)
            ["x", a], s -> QxSim.x(s, n, a)
            ["rx", a, t], s -> QxSim.rx(s, n, a, t)
            ["ry", a, t], s -> QxSim.ry(s, n, a, t)
            ["rz", a, t], s -> QxSim.rz(s, n, a, t)
            ["cx", a, b], s -> QxSim.cx(s, n, a, b)
            ["cz", a, b], s -> QxSim.cz(s, n, a, b)
            ["cp", a, b, t], s -> QxSim.cp(s, n, a, b, t)
            ["swap", a, b], s -> QxSim.swap(s, n, a, b)
          end)
        %{n: n, case: c, gates: gates, amps: state |> Tuple.to_list() |> Enum.map(fn {re, im} -> [re, im] end)}
      end
    File.write!(Path.join(out, "xval.json"), JSON.encode!(cases))
    IO.puts("xval: #{length(cases)} circuits")

  "stab" ->
    # 100,000 shots of each memory circuit at p = 0.001. Timed: the reference tableau and the sampling (batches
    # of 1,024 in parallel processes on all cores), inside the VM; one warm-up, then 5 runs, median and minimum.
    rows =
      for {code, d} <- [{:repetition, 9}, {:surface, 3}, {:surface, 5}, {:surface, 7}, {:surface, 9}] do
        {c, n} = if code == :surface, do: QxStab.surface_memory(d, d, 0.001), else: QxStab.repetition_memory(d, d, 0.001)
        name = "#{code}_d#{d}"
        File.write!(Path.join(out, name <> ".stim"), QxStab.to_stim(c))
        QxStab.sample(c, n, 100_000, seed: 1)
        ts = for k <- 1..5, do: elem(:timer.tc(fn -> QxStab.sample(c, n, 100_000, seed: k) end), 0) / 1000
        med = Enum.at(Enum.sort(ts), 2)
        IO.puts("#{name}: median #{round(med)} ms")
        %{name: name, code: code, d: d, qubits: n, shots: 100_000, median_ms: med, min_ms: Enum.min(ts), runs: 5,
          schedulers: System.schedulers_online()}
      end
    File.write!(Path.join(out, "stab.json"), JSON.encode!(rows))
end
