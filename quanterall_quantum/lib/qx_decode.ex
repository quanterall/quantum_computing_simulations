defmodule QxDecode do
  @moduledoc """
  A union-find decoder (Delfosse and Nickerson, "Almost-linear time decoding algorithm for topological codes",
  Quantum 5, 595, 2021), weighted, on a decoding graph from Stim's detector error model (tools/dem_graph.py).

  For one syndrome — the list of detectors that fired:

    1. **Growth.** Every fired detector starts a cluster. In rounds, every cluster with an odd number of fired
       detectors that does not yet touch the boundary grows by one unit along each of its edges. An edge is
       fully grown when its growth reaches its length, and the clusters at its ends merge; an edge to the
       boundary makes its cluster neutral. Lengths come from the error probabilities: an edge of probability p
       has weight ln((1−p)/p), and its length is that weight relative to the lightest edge, times 4 (so growth
       follows the weights in steps of a quarter of the lightest edge).
    2. **Peeling.** In each cluster a minimum-weight spanning tree of its grown edges (Kruskal, on the true
       weights; the boundary one vertex of it), rooted at the boundary if the cluster touches it; leaves are
       removed one by one, and an edge goes into the correction when its leaf carries an odd number of fired
       detectors (which passes on to its parent). Ties in growth are thus settled by weight, not by order.
    3. **Prediction.** The correction flips the logical observable when an odd number of its edges do.

  `decode/2` returns 0 or 1: whether the logical observable is predicted to have flipped.

  **Erasures** (`decode/3`): edges known to have been erased (a lost photon, a heralded leak) carry an error with
  probability ½ — no information. They start fully grown, their clusters already formed, and weigh nothing in the
  peeling tree. For pure erasure this is the peeling decoder of Delfosse and Zémor, which is maximum-likelihood;
  with Pauli noise as well, ordinary growth continues from the erased clusters.
  """
  import Bitwise

  @doc "The graph from tools/dem_graph.py's JSON: %{detectors: n, edges: tuple, adj: tuple (vertex n = the boundary)}"
  def load(json) do
    %{"detectors" => n, "edges" => es} = JSON.decode!(json)
    ws = Enum.map(es, fn [_u, _v, p, _o] -> :math.log((1 - p) / p) end)
    wmin = Enum.min(ws)
    edges =
      Enum.zip(es, ws)
      |> Enum.map(fn {[u, v, _p, o], w} -> {u, if(v < 0, do: n, else: v), max(4, round(4 * w / wmin)), o, w} end)
      |> List.to_tuple()
    adj =
      edges
      |> Tuple.to_list()
      |> Enum.with_index()
      |> Enum.reduce(:array.new(n + 1, default: []), fn {{u, v, _l, _o, _w}, e}, a ->
        a = :array.set(u, [e | :array.get(u, a)], a)
        if v == n, do: a, else: :array.set(v, [e | :array.get(v, a)], a)
      end)
    adj = for v <- 0..n, do: :array.get(v, adj)
    %{detectors: n, edges: edges, adj: List.to_tuple(adj)}
  end

  @doc "Decodes one syndrome (the fired detectors); 1 if the logical observable is predicted to have flipped."
  def decode(g, fired), do: decode(g, fired, [])

  def decode(_g, [], []), do: 0
  def decode(g, fired, erased) do
    b = g.detectors
    st = %{
      parent: Map.new(fired, &{&1, &1}),
      cl: Map.new(fired, &{&1, %{odd: true, bnd: false, verts: [&1]}}),
      support: %{},
      defect: MapSet.new(fired),
      erased: MapSet.new(erased)
    }
    # erased edges: fully grown from the start, their ends joined into clusters
    st = Enum.reduce(erased, st, fn e, st ->
      {u, _v, len, _o, _w} = elem(g.edges, e)
      st = st |> join(u) |> Map.update!(:support, &Map.put(&1, e, len))
      fuse(g, b, e, st)
    end)
    st = grow(g, b, st)
    peel(g, b, st)
  end

  # ── growth ────────────────────────────────────────────────────────────────────────────────────────────
  defp grow(g, b, st) do
    active = for {r, c} <- st.cl, c.odd and not c.bnd, do: r
    if active == [] do
      st
    else
      {support, fused} =
        Enum.reduce(active, {st.support, []}, fn r, acc ->
          Enum.reduce(st.cl[r].verts, acc, fn v, acc ->
            Enum.reduce(elem(g.adj, v), acc, fn e, {sup, fu} ->
              len = elem(elem(g.edges, e), 2)
              s = Map.get(sup, e, 0)
              cond do
                s >= len -> {sup, fu}
                s + 1 >= len -> {Map.put(sup, e, s + 1), [e | fu]}
                true -> {Map.put(sup, e, s + 1), fu}
              end
            end)
          end)
        end)
      st = Enum.reduce(fused, %{st | support: support}, &fuse(g, b, &1, &2))
      grow(g, b, st)
    end
  end

  defp fuse(g, b, e, st) do
    {u, v, _l, _o, _w} = elem(g.edges, e)
    cond do
      v == b ->
        {r, st} = find(st, u)
        put_in(st, [:cl, r, :bnd], true)
      true ->
        st = st |> join(u) |> join(v)
        {ru, st} = find(st, u)
        {rv, st} = find(st, v)
        if ru == rv do
          st
        else
          cu = st.cl[ru]
          cv = st.cl[rv]
          {big, small, cb, cs} = if length(cu.verts) >= length(cv.verts), do: {ru, rv, cu, cv}, else: {rv, ru, cv, cu}
          merged = %{odd: cb.odd != cs.odd, bnd: cb.bnd or cs.bnd, verts: cs.verts ++ cb.verts}
          %{st | parent: Map.put(st.parent, small, big), cl: st.cl |> Map.delete(small) |> Map.put(big, merged)}
        end
    end
  end

  # a vertex reached by a fully grown edge joins as a neutral one-vertex cluster, if it is in none yet
  defp join(st, v) do
    if Map.has_key?(st.parent, v),
      do: st,
      else: %{st | parent: Map.put(st.parent, v, v), cl: Map.put(st.cl, v, %{odd: false, bnd: false, verts: [v]})}
  end

  defp find(st, v) do
    p = Map.fetch!(st.parent, v)
    if p == v do
      {v, st}
    else
      {r, st} = find(st, p)
      {r, %{st | parent: Map.put(st.parent, v, r)}}                      # path compression
    end
  end

  # ── peeling ───────────────────────────────────────────────────────────────────────────────────────────
  defp peel(g, b, st) do
    grown = for {e, s} <- st.support, s >= elem(elem(g.edges, e), 2), into: MapSet.new(), do: e
    Enum.reduce(st.cl, 0, fn {_r, c}, parity ->
      bxor(parity, peel_cluster(g, b, grown, c, st.defect, Map.get(st, :erased, MapSet.new())))
    end)
  end

  defp peel_cluster(g, b, grown, c, defect, erased) do
    inside = MapSet.new(c.verts)
    # the cluster's grown edges (to the boundary too), lightest first
    cand =
      c.verts
      |> Enum.flat_map(fn v -> elem(g.adj, v) end)
      |> Enum.uniq()
      |> Enum.filter(fn e ->
        {x, y, _l, _o, _w} = elem(g.edges, e)
        MapSet.member?(grown, e) and MapSet.member?(inside, x) and (y == b or MapSet.member?(inside, y))
      end)
      |> Enum.sort_by(fn e -> if MapSet.member?(erased, e), do: 0.0, else: elem(elem(g.edges, e), 4) end)
    # Kruskal: a minimum-weight spanning forest, the boundary one vertex
    {tree, _uf} =
      Enum.reduce(cand, {%{}, %{}}, fn e, {tree, uf} ->
        {x, y, _l, _o, _w} = elem(g.edges, e)
        {rx, uf} = uf_find(uf, x)
        {ry, uf} = uf_find(uf, y)
        if rx == ry do
          {tree, uf}
        else
          tree = tree |> Map.update(x, [{y, e}], &[{y, e} | &1]) |> Map.update(y, [{x, e}], &[{x, e} | &1])
          {tree, Map.put(uf, rx, ry)}
        end
      end)
    root = if c.bnd, do: b, else: hd(c.verts)
    order = tree_order(tree, [root], MapSet.new([root]), [])
    # peel from the leaves: a vertex with an odd mark puts its parent edge in the correction
    {_marks, flips} =
      Enum.reduce(order, {MapSet.new(Enum.filter(c.verts, &MapSet.member?(defect, &1))), 0}, fn {v, e, par}, {marks, flips} ->
        if MapSet.member?(marks, v) do
          o = elem(elem(g.edges, e), 3)
          marks = MapSet.delete(marks, v)
          marks = cond do
            par == b -> marks
            MapSet.member?(marks, par) -> MapSet.delete(marks, par)
            true -> MapSet.put(marks, par)
          end
          {marks, bxor(flips, o)}
        else
          {marks, flips}
        end
      end)
    flips
  end

  defp uf_find(uf, x) do
    case Map.get(uf, x) do
      nil -> {x, uf}
      p -> {r, uf} = uf_find(uf, p); {r, Map.put(uf, x, r)}
    end
  end

  # breadth-first over the tree from the root; returned leaves first: {vertex, edge to parent, parent}
  defp tree_order(_tree, [], _seen, acc), do: acc
  defp tree_order(tree, frontier, seen, acc) do
    {next, seen, acc} =
      Enum.reduce(frontier, {[], seen, acc}, fn v, {nx, seen, acc} ->
        Enum.reduce(Map.get(tree, v, []), {nx, seen, acc}, fn {w, e}, {nx, seen, acc} ->
          if MapSet.member?(seen, w), do: {nx, seen, acc}, else: {[w | nx], MapSet.put(seen, w), [{w, e, v} | acc]}
        end)
      end)
    tree_order(tree, next, seen, acc)
  end

  # ── batches ───────────────────────────────────────────────────────────────────────────────────────────

  @doc """
  The syndromes of a packed sample (QxStab.write_packed's format): a list, one per shot, of fired detectors;
  and the observable's bit per shot.
  """
  def read_packed(path, detectors, shots) do
    nbytes = div(shots + 7, 8)
    data = File.read!(path)
    rows = for j <- 0..detectors, do: binary_part(data, j * nbytes, nbytes)
    {det_rows, [obs_row]} = Enum.split(rows, detectors)
    lists =
      det_rows
      |> Enum.with_index()
      |> Enum.reduce(:array.new(shots, default: []), fn {row, j}, a ->
        Enum.reduce(set_bits(row), a, fn k, a -> if k < shots, do: :array.set(k, [j | :array.get(k, a)], a), else: a end)
      end)
    {for(k <- 0..(shots - 1), do: Enum.reverse(:array.get(k, lists))), set_bits(obs_row) |> Enum.filter(&(&1 < shots)) |> MapSet.new()}
  end

  defp set_bits(row) do
    for <<byte <- row>>, reduce: {0, []} do
      {i, acc} ->
        acc = if byte == 0, do: acc, else: Enum.reduce(0..7, acc, fn k, acc -> if (byte >>> k &&& 1) == 1, do: [i * 8 + k | acc], else: acc end)
        {i + 1, acc}
    end
    |> elem(1)
    |> Enum.reverse()
  end

  @doc "Decodes many syndromes on all cores, in chunks; returns the predictions in order."
  def decode_all(g, syndromes, chunk \\ 1000) do
    syndromes
    |> Enum.chunk_every(chunk)
    |> Task.async_stream(fn ch -> Enum.map(ch, &decode(g, &1)) end, max_concurrency: System.schedulers_online(), timeout: :infinity)
    |> Enum.flat_map(fn {:ok, r} -> r end)
  end
end
