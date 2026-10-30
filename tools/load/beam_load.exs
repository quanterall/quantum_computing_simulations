# The BEAM side of the decoder load test: QxDecode, every syndrome its own process.
#   docker exec -w /app quanterall_quantum sh -c 'MIX_ENV=prod mix run --no-start /tmp/load/beam_load.exs DIR MODE ARGS...'
#   DIR: graph.json (tools/dem_graph.py) and syndromes.bin (tools/load/export_inputs.py)
#   capacity                          throughput, syndromes decoded in chunks on all schedulers
#   open RATE COUNT SEED [noop]       Poisson arrivals at RATE per second; a dispatcher spins on the clock and
#                                     spawns one process per arrival; latency = completion − scheduled arrival
#   closed INFLIGHT COUNT [noop]      INFLIGHT always outstanding; latency = completion − submission
# Garbage collections (count, words reclaimed) are reported for the run.
import Bitwise
[dir, mode | args] = System.argv()
g = QxDecode.load(File.read!(Path.join(dir, "graph.json")))

# records until the remaining bytes equal the number read (the observable bytes follow the records)
parse = fn
  parse, b, acc, cnt when byte_size(b) == cnt -> {acc |> Enum.reverse() |> List.to_tuple(), b}
  parse, <<c::little-16, ids::binary-size(c * 2), more::binary>>, acc, cnt ->
    parse.(parse, more, [for(<<i::little-16 <- ids>>, do: i) | acc], cnt + 1)
end
read_syn = fn -> parse.(parse, File.read!(Path.join(dir, "syndromes.bin")), [], 0) end
{syn, _obs} = read_syn.()
n = tuple_size(syn)
noop = List.last(args) == "noop"
# the graph and the syndromes in persistent_term: a spawned process reads them without copying (a closure
# capturing them would copy them into every process)
:persistent_term.put(:load_g, g)
:persistent_term.put(:load_syn, syn)
work = fn k ->
  if noop, do: 0, else: QxDecode.decode(:persistent_term.get(:load_g), elem(:persistent_term.get(:load_syn), rem(k, n)))
end
now = fn -> :erlang.monotonic_time(:nanosecond) end
gc = fn -> {c, w, _} = :erlang.statistics(:garbage_collection); {c, w} end

report = fn mode, lats, secs, count, extra ->
  s = lats |> Enum.sort() |> List.to_tuple()
  pct = fn p -> elem(s, min(tuple_size(s) - 1, floor(p * tuple_size(s)))) end
  over = Enum.count(lats, &(&1 > 1000.0))
  Map.merge(%{mode: mode, count: count, seconds: secs, throughput: count / secs, p50_us: pct.(0.5), p99_us: pct.(0.99),
              p999_us: pct.(0.999), max_us: elem(s, tuple_size(s) - 1), over_1ms: over / length(lats), noop: noop}, extra)
  |> JSON.encode!() |> IO.puts()
end

case mode do
  "capacity" ->
    syns = Tuple.to_list(syn)
    QxDecode.decode_all(g, Enum.take(syns, 5000))
    t0 = now.()
    QxDecode.decode_all(g, syns ++ syns)
    secs = (now.() - t0) / 1.0e9
    IO.puts(JSON.encode!(%{mode: "capacity", schedulers: System.schedulers_online(), count: 2 * n, seconds: secs, throughput: 2 * n / secs}))

  "open" ->
    [rate, count, seed | _] = args
    {rate, count, seed} = {String.to_float(rate), String.to_integer(count), String.to_integer(seed)}
    :rand.seed(:exsss, {seed, 1, 2})
    {arrivals, _} = Enum.map_reduce(1..count, 0.0, fn _, at -> at = at - :math.log(1.0 - :rand.uniform_real()) / rate; {at, at} end)
    me = self()
    collector = spawn(fn ->
      lats = for _ <- 1..count, do: (receive do {:done, l} -> l end)
      send(me, {:lats, lats})
    end)
    {gc0c, gc0w} = gc.()
    t0 = now.()
    spin = fn spin, target -> if now.() < target, do: spin.(spin, target), else: :ok end
    arrivals
    |> Enum.with_index()
    |> Enum.each(fn {at, k} ->
      due = t0 + round(at * 1.0e9)
      spin.(spin, due)
      spawn(fn -> work.(k); send(collector, {:done, (now.() - due) / 1000}) end)
    end)
    lats = receive do {:lats, l} -> l end
    secs = (now.() - t0) / 1.0e9
    {gc1c, gc1w} = gc.()
    report.("open", lats, secs, count, %{rate: rate, gc_count: gc1c - gc0c, gc_words: gc1w - gc0w, schedulers: System.schedulers_online()})

  "closed" ->
    [inflight, count | _] = args
    {inflight, count} = {String.to_integer(inflight), String.to_integer(count)}
    me = self()
    go = fn k -> t = now.(); spawn(fn -> work.(k); send(me, {:done, (now.() - t) / 1000}) end) end
    {gc0c, gc0w} = gc.()
    t0 = now.()
    Enum.each(0..(min(inflight, count) - 1), go)
    lats =
      Enum.map(0..(count - 1), fn i ->
        l = receive do {:done, l} -> l end
        if i + inflight < count, do: go.(i + inflight)
        l
      end)
    secs = (now.() - t0) / 1.0e9
    {gc1c, gc1w} = gc.()
    report.("closed", lats, secs, count, %{inflight: inflight, gc_count: gc1c - gc0c, gc_words: gc1w - gc0w, schedulers: System.schedulers_online()})
end
