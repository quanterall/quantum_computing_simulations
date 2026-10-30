# Our union-find decoder (lib/qx_decode.ex) on the threshold study's samples: for each graph in GRAPHS, the
# syndromes of the matching sample in BENCH, decoded on all cores; predictions written for tools/decoder_compare.py.
#   docker exec -w /app quanterall_quantum sh -c 'MIX_ENV=prod mix run --no-start /tmp/decode_runs.exs BENCH GRAPHS OUT'
[bench, graphs, out] = System.argv()
File.mkdir_p!(out)
study = bench |> Path.join("study.json") |> File.read!() |> JSON.decode!()
rows =
  for run <- study, File.exists?(Path.join(graphs, run["name"] <> ".json")) do
    g = QxDecode.load(File.read!(Path.join(graphs, run["name"] <> ".json")))
    {syn, obs} = QxDecode.read_packed(Path.join(bench, run["name"] <> ".bin"), g.detectors, run["shots"])
    QxDecode.decode_all(g, Enum.take(syn, 2000))                                   # warm-up
    {us, pred} = :timer.tc(fn -> QxDecode.decode_all(g, syn) end)
    {us1, _} = :timer.tc(fn -> Enum.each(Enum.take(syn, 10_000), &QxDecode.decode(g, &1)) end)   # one process
    fails = pred |> Enum.with_index() |> Enum.count(fn {p, k} -> p != (if MapSet.member?(obs, k), do: 1, else: 0) end)
    File.write!(Path.join(out, run["name"] <> ".pred"), :erlang.list_to_binary(pred))
    IO.puts("#{run["name"]}: logical #{fails / run["shots"]}, #{round(run["shots"] / (us / 1.0e6))} syndromes/s on #{System.schedulers_online()} schedulers, #{round(10_000 / (us1 / 1.0e6))}/s in one process")
    %{name: run["name"], d: run["d"], p: run["p"], shots: run["shots"], uf_failures: fails, uf_all_cores_s: us / 1.0e6,
      uf_one_process_per_syndrome_us: us1 / 10_000, schedulers: System.schedulers_online()}
  end
File.write!(Path.join(out, "uf.json"), JSON.encode!(rows))
