# Circuit-level erasure: surface-code memories (d rounds, Pauli noise p) with heralded erasure pe on the data each
# round, sampled by QxStab. Writes NAME.stim, NAME.bin (QxStab.write_packed: detectors incl. heralds, observable)
# and NAME.heralds.json (which detector indices are heralds).
#   mix run --no-start /tmp/erasure/circuit_runs.exs OUT SHOTS
[out, shots] = System.argv()
shots = String.to_integer(shots)
File.mkdir_p!(out)
for d <- [3, 5, 7], pe <- [0.0, 0.01, 0.02, 0.04, 0.06, 0.08] do
  {c, n} = QxStab.surface_memory(d, d, 0.001, pe, true)
  name = "cl_d#{d}_pe#{pe}"
  File.write!(Path.join(out, name <> ".stim"), QxStab.to_stim(c))
  # the detectors that are heralds: those whose single record is a herald bit
  {herald_recs, _} =
    Enum.reduce(c, {MapSet.new(), 0}, fn
      {:m, qs}, {h, k} -> {h, k + length(qs)}
      {:heralded_erase, qs, _}, {h, k} -> {Enum.reduce(0..(length(qs) - 1), h, &MapSet.put(&2, k + &1)), k + length(qs)}
      _, acc -> acc
    end)
  heralds = for {{:detector, [r]}, i} <- Enum.with_index(Enum.filter(c, &match?({:detector, _}, &1))), MapSet.member?(herald_recs, r), do: i
  File.write!(Path.join(out, name <> ".heralds.json"), JSON.encode!(heralds))
  r = QxStab.sample(c, n, shots, seed: 777 + d * 100 + round(pe * 1000))
  QxStab.write_packed(Path.join(out, name <> ".bin"), r)
  IO.puts("#{name}: #{length(heralds)} heralds, #{length(r.detectors)} detectors")
end
