# QxCV's numbers for tools/optics/xcheck_cv.py, as one JSON line. Run in the container:
#   mix run --no-start /tmp/cv/cv_beam.exs [path/to/qx_cv.ex — to test a copy before the image is rebuilt]
case System.argv() do
  [f | _] -> Code.compile_file(f)
  _ -> :ok
end
:rand.seed(:exsss, 7)
q = QxCV.grid(48.0, 8192)
d0 = QxCV.gkp_density(0, 0.01, 0.05, q)
p = QxCV.grid(16.0, 4096)
dp = QxCV.gkp_plus_p_density(0.01, 0.05, p, 48.0)
z = for s <- [0.0, 0.2, 0.4, 0.6], do: QxCV.readout_exact(d0, q, s, QxCV.alpha(), false)
x = for s <- [0.0, 0.2, 0.4, 0.6], do: QxCV.readout_exact(dp, p, s, :math.pi() / QxCV.alpha(), false)
closed = for s <- [0.2, 0.4], do: QxCV.readout_closed(s, QxCV.alpha())
cats = for {b, par, eta} <- [{2.0, 0, 1.0}, {2.0, 1, 0.9}, {1.3, 0, 0.7}] do
  dist = QxCV.cat_distribution(b, par, eta)
  pi_ = dist |> Enum.with_index() |> Enum.reduce(0.0, fn {pp, n}, acc -> acc + if(rem(n, 2) == 0, do: pp, else: -pp) end)
  %{dist: Enum.take(dist, 20), parity: pi_, closed: QxCV.parity_closed(b, par, eta), total: Enum.sum(dist)}
end
comb = QxCV.run_comb(1.5, 3.0, :math.sqrt(:math.pi()), 0.0, 0.0)
mc = QxCV.run_gkp(10.0, 0.01, 0.4, :z, 20000)
IO.puts(:thoas.encode(%{z: z, x: x, closed: closed, cats: cats, d0_sample: Enum.take_every(d0, 97), dp_sample: Enum.take_every(dp, 37),
  comb_y0: %{gaps: comb.gaps, pred: comb.gaps_predicted, p_y: comb.p_y_integral, dens: comb.density},
  mc: %{wrong: mc.wrong, shots: mc.shots, exact: mc.error_exact, closed: mc.error_closed_form, us: mc.time_us}}))
