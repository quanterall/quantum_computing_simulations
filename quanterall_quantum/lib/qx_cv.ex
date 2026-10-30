defmodule QxCV do
  @moduledoc """
  One optical mode, continuous variables, on a q-grid — Phase 7 of docs/photonic-plan.md: GKP states (the review's
  Eq. 38) read out under Gaussian shift noise, cat states (Eq. 75) and their parity under photon loss, and the
  homodyne-heralded comb (Eqs. 66–73). The Python twin is qsim/qsim/cv.py (TechUPGRADE's quantum/cv.py); this one is
  written separately: the p-side of a GKP state is its Fourier transform summed analytically over the Gaussian peaks
  (no FFT), the cat's photon statistics come from the generating function. Conventions: hbar = 1,
  q = (a + a†)/√2, vacuum Var(q) = 1/2, α = √π.
  """

  @alpha :math.sqrt(:math.pi())

  def alpha, do: @alpha

  def grid(half, n), do: for(j <- 0..(n - 1), do: -half + 2 * half * j / n)

  # ── GKP (Eq. 38) ──
  defp centres(mu, half) do
    smax = trunc(half / (2 * @alpha)) + 2
    for s <- -smax..smax, do: (2 * s + mu) * @alpha
  end

  @doc "the density |⟨q|μ̃⟩|² on the grid; μ = :plus gives (|0̃⟩ + |1̃⟩)/√2 — peaks of variance pv, envelope Δ²"
  def gkp_density(mu, delta2, pv, q) do
    half = -hd(q)
    cs = case mu do
      :plus -> centres(0, half) ++ centres(1, half)
      m -> centres(m, half)
    end
    ws = Enum.map(cs, fn c -> {c, :math.exp(-0.5 * delta2 * c * c)} end)
    amp = Enum.map(q, fn x -> Enum.reduce(ws, 0.0, fn {c, w}, acc -> d = x - c; if abs(d) > 12 * :math.sqrt(pv) + 1, do: acc, else: acc + w * :math.exp(-d * d / (4 * pv)) end) end)
    normalise(Enum.map(amp, &(&1 * &1)), q)
  end

  @doc """
  the p-density of (|0̃⟩ + |1̃⟩)/√2: each q-peak e^{−(q−c)²/4v} has the transform ∝ e^{−ipc} e^{−v p²}, so
  ψ(p) ∝ e^{−v p²} Σ_c w_c e^{−ipc}, summed exactly over the peaks
  """
  def gkp_plus_p_density(delta2, pv, p, q_half) do
    ws = Enum.map(centres(0, q_half) ++ centres(1, q_half), fn c -> {c, :math.exp(-0.5 * delta2 * c * c)} end)
    dens = Enum.map(p, fn k ->
      {re, im} = Enum.reduce(ws, {0.0, 0.0}, fn {c, w}, {r, i} -> {r + w * :math.cos(k * c), i - w * :math.sin(k * c)} end)
      g = :math.exp(-pv * k * k)
      g * g * (re * re + im * im)
    end)
    normalise(dens, p)
  end

  defp normalise(dens, x) do
    dx = Enum.at(x, 1) - hd(x)
    tot = Enum.sum(dens) * dx
    Enum.map(dens, &(&1 / tot))
  end

  defp phi(x), do: 0.5 * (1 + :math.erf(x / :math.sqrt(2)))

  @doc "P(round((x + u)/spacing) has the wrong parity), u ~ N(0, σ²), σ > 0"
  def wrong_parity_prob(x, sigma, spacing, want_odd) do
    k0 = floor(x / spacing)
    span = trunc(8 * sigma / spacing) + 2
    Enum.reduce((k0 - span)..(k0 + span), 0.0, fn k, acc ->
      if (Integer.mod(k, 2) == 1) == want_odd, do: acc,
      else: acc + phi(((k + 0.5) * spacing - x) / sigma) - phi(((k - 0.5) * spacing - x) / sigma)
    end)
  end

  @doc "∫ ρ(x) P(wrong parity | x) dx on the grid; σ = 0 integrates the step from the cumulative distribution"
  def readout_exact(dens, x, sigma, spacing, want_odd) do
    dx = Enum.at(x, 1) - hd(x)
    if sigma == 0 do
      cdf = dens |> Enum.chunk_every(2, 1, :discard) |> Enum.scan(0.0, fn [a, b], acc -> acc + (a + b) / 2 * dx end)
      xs = List.to_tuple(x); cs = List.to_tuple([0.0 | cdf])
      f = fn v -> interp(xs, cs, v) end
      {klo, khi} = {floor(hd(x) / spacing) - 1, ceil(List.last(x) / spacing) + 1}
      for k <- klo..khi, (Integer.mod(k, 2) == 1) != want_odd, reduce: 0.0 do
        acc -> acc + f.((k + 0.5) * spacing) - f.((k - 0.5) * spacing)
      end
    else
      Enum.zip(dens, x) |> Enum.reduce(0.0, fn {d, v}, acc -> if d < 1.0e-16, do: acc, else: acc + d * wrong_parity_prob(v, sigma, spacing, want_odd) end) |> Kernel.*(dx)
    end
  end

  defp interp(xs, cs, v) do
    n = tuple_size(xs)
    cond do
      v <= elem(xs, 0) -> elem(cs, 0)
      v >= elem(xs, n - 1) -> elem(cs, n - 1)
      true ->
        dx = elem(xs, 1) - elem(xs, 0)
        i = min(trunc((v - elem(xs, 0)) / dx), n - 2)
        t = (v - elem(xs, i)) / dx
        elem(cs, i) + t * (elem(cs, i + 1) - elem(cs, i))
    end
  end

  @doc "an ideal comb with Gaussian shifts of σ_eff: Σ_k [Φ((2k + 3/2)a/σ) − Φ((2k + 1/2)a/σ)]"
  def readout_closed(sigma_eff, spacing) do
    Enum.reduce_while(Stream.iterate(0, &(&1 + 1)), 0.0, fn k, acc ->
      t = phi((2 * k + 1.5) * spacing / sigma_eff) - phi((2 * k + 0.5) * spacing / sigma_eff) +
          phi((-2 * k - 0.5) * spacing / sigma_eff) - phi((-2 * k - 1.5) * spacing / sigma_eff)
      if t < 1.0e-17 or k > 10_000, do: {:halt, acc + t}, else: {:cont, acc + t}
    end)
  end

  defp sampler(dens, x) do
    dx = Enum.at(x, 1) - hd(x)
    {cdf, tot} = Enum.map_reduce(dens, 0.0, fn d, acc -> {acc + d * dx, acc + d * dx} end)
    {List.to_tuple(Enum.map(cdf, &(&1 / tot))), List.to_tuple(x), dx}
  end

  defp draw({cdf, xs, dx}) do
    u = :rand.uniform_real()
    i = bsearch(cdf, u, 0, tuple_size(cdf) - 1)
    lo = if i == 0, do: 0.0, else: elem(cdf, i - 1)
    elem(xs, i) - dx + (u - lo) / max(elem(cdf, i) - lo, 1.0e-300) * dx
  end

  defp bsearch(_c, _u, lo, hi) when lo >= hi, do: lo
  defp bsearch(c, u, lo, hi) do
    mid = div(lo + hi, 2)
    if elem(c, mid) < u, do: bsearch(c, u, mid + 1, hi), else: bsearch(c, u, lo, mid)
  end

  @doc """
  GKP readout under Gaussian shift noise: Z (q, |0̃⟩, spacing α) or X (p, |+̃⟩, spacing π/α). The exact integral,
  `shots` sampled readouts, and the closed form for an ideal comb with σ_eff² = σ² + the peak variance.
  """
  def run_gkp(db, delta2, sigma, basis, shots) do
    pv = 0.5 * :math.pow(10, -db / 10)
    {dens, x, spacing, peak_var, half} =
      case basis do
        :x ->
          p = grid(16.0, 4096)
          {gkp_plus_p_density(delta2, pv, p, 48.0), p, :math.pi() / @alpha, delta2 / 2, 16.0}
        _ ->
          q = grid(48.0, 8192)
          {gkp_density(0, delta2, pv, q), q, @alpha, pv, 48.0}
      end
    t0 = System.monotonic_time(:microsecond)
    exact = readout_exact(dens, x, sigma, spacing, false)
    s = sampler(dens, x)
    wrong = Enum.count(1..shots, fn _ ->
      v = draw(s) + if(sigma > 0, do: sigma * :rand.normal(), else: 0.0)
      Integer.mod(round(v / spacing), 2) == 1
    end)
    closed = readout_closed(:math.sqrt(sigma * sigma + peak_var), spacing)
    step = div(length(x), 400)
    plot = x |> Enum.zip(dens) |> Enum.take_every(step) |> Enum.filter(fn {v, _} -> abs(v) <= min(half, 24.0) end)
           |> Enum.map(fn {v, d} -> [Float.round(v, 3), Float.round(d, 5)] end)
    %{status: "ok", name: "optics_gkp", basis: Atom.to_string(basis), squeezing_db: db, delta2: delta2, sigma: sigma,
      shots: shots, wrong: wrong, error_sampled: wrong / shots, error_exact: exact, error_closed_form: closed,
      peak_variance: peak_var, spacing: spacing, density: plot, time_us: System.monotonic_time(:microsecond) - t0,
      metadata: %{description: "A GKP qubit (Eq. 38: peaks at multiples of √π, " <> "#{db} dB" <> ") read out after a Gaussian shift of σ: the exact error, sampled readouts and the closed form for an ideal comb with σ_eff² = σ² + peak variance"}}
  end

  # ── cat states (Eq. 75) and photon loss ──
  @doc """
  photon-number distribution of the even (parity 0) or odd cat of real amplitude β after loss (transmissivity η):
  P(n) = e^{−β²} β^{2n}/n! · 2/N± on the right parity; loss thins it binomially. ⟨Π⟩ = cosh(β²(1−2η))/cosh β² (even),
  sinh(β²(1−2η))/sinh β² (odd), from the generating function — checked against the thinned distribution.
  """
  def cat_distribution(beta, parity, eta, nmax \\ 80) do
    b2 = beta * beta
    lf = log_factorials(nmax)
    nrm = 2 + (if parity == 0, do: 2, else: -2) * :math.exp(-2 * b2)
    pn = for n <- 0..nmax, do: (if rem(n, 2) == parity, do: 2 * :math.exp(-b2 + n * :math.log(b2) - elem(lf, n)) * 2 / nrm, else: 0.0)
    for m <- 0..nmax do
      Enum.with_index(pn) |> Enum.reduce(0.0, fn {p, n}, acc ->
        if n < m or p == 0.0, do: acc,
        else: acc + p * :math.exp(elem(lf, n) - elem(lf, m) - elem(lf, n - m) + m * safe_log(eta) + (n - m) * safe_log(1 - eta))
      end)
    end
  end

  defp safe_log(x) when x <= 0, do: -1.0e300
  defp safe_log(x), do: :math.log(x)

  defp log_factorials(nmax) do
    Enum.scan(0..nmax, 0.0, fn n, acc -> if n < 2, do: 0.0, else: acc + :math.log(n) end) |> List.to_tuple()
  end

  def parity_closed(beta, parity, eta) do
    b2 = beta * beta
    if parity == 0, do: :math.cosh(b2 * (1 - 2 * eta)) / :math.cosh(b2), else: :math.sinh(b2 * (1 - 2 * eta)) / :math.sinh(b2)
  end

  def run_cat(beta, parity, eta, shots) do
    t0 = System.monotonic_time(:microsecond)
    dist = cat_distribution(beta, parity, eta)
    par = dist |> Enum.with_index() |> Enum.reduce(0.0, fn {p, n}, acc -> acc + (if rem(n, 2) == 0, do: p, else: -p) end)
    cdf = dist |> Enum.scan(&(&1 + &2)) |> List.to_tuple()
    counts = for(_ <- 1..shots, do: bsearch(cdf, :rand.uniform_real() * elem(cdf, tuple_size(cdf) - 1), 0, tuple_size(cdf) - 1)) |> Enum.frequencies()
    flipped = counts |> Enum.filter(fn {n, _} -> rem(n, 2) != parity end) |> Enum.map(&elem(&1, 1)) |> Enum.sum()
    mean = dist |> Enum.with_index() |> Enum.reduce(0.0, fn {p, n}, acc -> acc + n * p end)
    %{status: "ok", name: "optics_cat", beta: beta, parity: if(parity == 0, do: "even", else: "odd"), eta: eta, shots: shots,
      distribution: dist |> Enum.take(31) |> Enum.map(&Float.round(&1, 6)), counts: Map.new(counts, fn {k, v} -> {Integer.to_string(k), v} end),
      parity_exact: par, parity_closed_form: parity_closed(beta, parity, eta), parity_flip_sampled: flipped / shots,
      parity_flip_exact: (1 - (if parity == 0, do: par, else: -par)) / 2, mean_photons: mean,
      overlap_alpha_minus_alpha: :math.exp(-4 * beta * beta), time_us: System.monotonic_time(:microsecond) - t0,
      metadata: %{description: "A cat state (Eq. 75) counted photon by photon after loss η: loss of a single photon flips the parity, the cat code's logical error"}}
  end

  # ── the homodyne-heralded comb (Eqs. 66–73) ──
  @doc """
  a phase-squeezed input S(r)|0⟩ (q anti-squeezed, Var(q) = e^{2r}/2), a probe of amplitude β coupled with strength θ,
  a homodyne outcome y sampled from P(y) (Eq. 73); the output density |Υ(q|y)|²|ψ_in(q)|² normalised (Eq. 72), with
  |Υ|² = (2π)^{−1/2} exp[−(y − 2β cos(θq − ϕ))²/2] — the square of Eq. 69 (Eq. 70 as printed lacks the /2)
  """
  def run_comb(r, beta, theta, phi, y_fixed \\ nil) do
    t0 = System.monotonic_time(:microsecond)
    q = grid(16.0, 3072)
    dq = Enum.at(q, 1) - hd(q)
    vin = :math.exp(2 * r) / 2
    din = Enum.map(q, fn x -> :math.exp(-x * x / (2 * vin)) / :math.sqrt(2 * :math.pi() * vin) end)
    u2 = fn x, y -> c = 2 * beta * :math.cos(theta * x - phi); :math.exp(-(y - c) * (y - c) / 2) / :math.sqrt(2 * :math.pi()) end
    ys = for j <- 0..600, do: -2 * beta - 6 + j * (4 * beta + 12) / 600
    dy = Enum.at(ys, 1) - hd(ys)
    py = Enum.map(ys, fn y -> Enum.zip(q, din) |> Enum.reduce(0.0, fn {x, d}, acc -> acc + u2.(x, y) * d end) |> Kernel.*(dq) end)
    y = y_fixed || (fn ->
      {cdf, _} = Enum.map_reduce(py, 0.0, fn p, acc -> {acc + p * dy, acc + p * dy} end)
      u = :rand.uniform_real() * List.last(cdf)
      i = Enum.find_index(cdf, &(&1 >= u)) || 0
      Enum.at(ys, i)
    end).()
    out = normalise(Enum.zip(q, din) |> Enum.map(fn {x, d} -> u2.(x, y) * d end), q)
    peaks = Enum.zip([q, out, tl(out) ++ [0.0], [0.0 | out]]) |> Enum.filter(fn {_x, o, nx, pr} -> o > pr and o >= nx end)
    top = Enum.max(out)
    pk = peaks |> Enum.filter(fn {_, o, _, _} -> o > 0.05 * top end) |> Enum.map(&elem(&1, 0))
    gaps = pk |> Enum.chunk_every(2, 1, :discard) |> Enum.map(fn [a, b] -> Float.round(b - a, 2) end) |> Enum.uniq() |> Enum.sort()
    pred = if abs(y) < 2 * beta, do: (ac = :math.acos(y / (2 * beta)); Enum.sort([2 * ac / theta, (2 * :math.pi() - 2 * ac) / theta])), else: []
    step = 4
    %{status: "ok", name: "optics_comb", r: r, beta: beta, theta: theta, phi: phi, y: y, p_y_integral: Enum.sum(py) * dy,
      peaks: length(pk), gaps: gaps, gaps_predicted: pred,
      density: q |> Enum.zip(out) |> Enum.take_every(step) |> Enum.map(fn {x, d} -> [Float.round(x, 3), Float.round(d, 5)] end),
      input: q |> Enum.zip(din) |> Enum.take_every(step) |> Enum.map(fn {x, d} -> [Float.round(x, 3), Float.round(d, 5)] end),
      time_us: System.monotonic_time(:microsecond) - t0,
      metadata: %{description: "Eqs. 66–73: a homodyne outcome y on a probe heralds a comb in the signal; y near 0 gives evenly spaced peaks (π/θ apart), other y pairs of peaks"}}
  end
end
