# QxDecode on code-capacity erasure samples (tools/erasure/codecap.py's format): per shot uint16 count + fired
# checks, uint16 count + erased edges; then one byte per shot, the actual logical flip. Prints failures as JSON.
#   docker exec -w /app -e MIX_ENV=prod quanterall_quantum mix run --no-start /tmp/erasure/decode_erasure.exs GRAPH SAMPLES SHOTS
[gpath, spath, shots] = System.argv()
shots = String.to_integer(shots)
g = QxDecode.load(File.read!(gpath))
parse = fn
  _parse, rest, acc, 0 -> {Enum.reverse(acc), rest}
  parse, <<a::little-16, f::binary-size(a * 2), b::little-16, e::binary-size(b * 2), more::binary>>, acc, k ->
    parse.(parse, more, [{for(<<i::little-16 <- f>>, do: i), for(<<i::little-16 <- e>>, do: i)} | acc], k - 1)
end
{recs, obs} = parse.(parse, File.read!(spath), [], shots)
:persistent_term.put(:er_g, g)
t0 = System.monotonic_time(:microsecond)
preds =
  recs
  |> Enum.chunk_every(1000)
  |> Task.async_stream(fn ch -> Enum.map(ch, fn {f, e} -> QxDecode.decode(:persistent_term.get(:er_g), f, e) end) end,
       max_concurrency: System.schedulers_online(), timeout: :infinity)
  |> Enum.flat_map(fn {:ok, r} -> r end)
us = System.monotonic_time(:microsecond) - t0
fails = Enum.zip(preds, :binary.bin_to_list(obs)) |> Enum.count(fn {p, o} -> p != o end)
IO.puts(JSON.encode!(%{shots: shots, failures: fails, seconds: us / 1.0e6}))
