// The compiled baseline for the decoder load test: the same union-find as QxDecode (lib/qx_decode.ex) —
// weighted growth in integer lengths, peeling over a minimum-weight spanning tree per cluster — in C++, with
// the load harness built in. Reads tools/load/export_inputs.py's graph.txt and syndromes.bin.
//
//   g++ -O2 -std=c++17 -pthread -o uf tools/load/uf.cpp
//   ./uf DIR predict OUT              decode all, one thread; write predictions (one byte each)
//   ./uf DIR capacity THREADS         throughput with THREADS workers taking syndromes as fast as they can
//   ./uf DIR open RATE COUNT WORKERS SEED [noop]    Poisson arrivals at RATE per second, a dispatcher thread
//   ./uf DIR closed INFLIGHT COUNT WORKERS [noop]   INFLIGHT syndromes always outstanding
// Latency = completion − arrival (open) or completion − submission (closed); printed as JSON.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

using Clock = std::chrono::steady_clock;
struct Edge { int u, v, len, obs; double w; };

struct Graph {
  int n = 0;                                  // detectors; vertex n is the boundary
  std::vector<Edge> edges;
  std::vector<std::vector<int>> adj;
};

static Graph load_graph(const std::string& path) {
  Graph g; std::ifstream f(path); int m; f >> g.n >> m;
  g.edges.resize(m); g.adj.assign(g.n + 1, {});
  for (int e = 0; e < m; e++) {
    Edge& x = g.edges[e]; f >> x.u >> x.v >> x.len >> x.obs >> x.w;
    g.adj[x.u].push_back(e); if (x.v != g.n) g.adj[x.v].push_back(e);
  }
  return g;
}

struct Syndromes { std::vector<std::vector<uint16_t>> fired; std::vector<uint8_t> obs; };

static Syndromes load_syndromes(const std::string& path, int shots_hint = -1) {
  std::ifstream f(path, std::ios::binary); std::vector<char> buf((std::istreambuf_iterator<char>(f)), {});
  Syndromes s; size_t i = 0;
  // shots are not stored: read records until the remaining bytes equal the number of records read (the obs bytes)
  while (i + 2 <= buf.size()) {
    if (buf.size() - i == s.fired.size()) break;
    uint16_t c; std::memcpy(&c, &buf[i], 2); i += 2;
    std::vector<uint16_t> ids(c); std::memcpy(ids.data(), &buf[i], 2 * c); i += 2 * c;
    s.fired.push_back(std::move(ids));
  }
  s.obs.assign(buf.begin() + i, buf.end());
  (void)shots_hint;
  return s;
}

// one decoder per thread: buffers sized to the graph, reset through "touched" lists
struct Decoder {
  const Graph& g;
  std::vector<int> parent, support, odd, bnd, mark, dsu, seen_e;
  std::vector<char> vseen;
  std::vector<std::vector<int>> verts, tree;
  std::vector<int> touched_v, touched_e, roots;
  explicit Decoder(const Graph& gr) : g(gr) {
    int V = g.n + 1, E = (int)g.edges.size();
    parent.assign(V, -1); support.assign(E, 0); odd.assign(V, 0); bnd.assign(V, 0); mark.assign(V, 0);
    dsu.assign(V, -1); seen_e.assign(E, 0); verts.assign(V, {}); tree.assign(V, {}); vseen.assign(V, 0);
  }
  int find(int v) { while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; } return v; }
  void join(int v) { if (parent[v] < 0) { parent[v] = v; verts[v] = {v}; odd[v] = 0; bnd[v] = 0; touched_v.push_back(v); roots.push_back(v); } }
  int dfind(int v) { while (dsu[v] >= 0 && dsu[v] != v) v = dsu[v]; return v; }

  int decode(const std::vector<uint16_t>& fired) {
    if (fired.empty()) return 0;
    const int B = g.n;
    for (int v : fired) { parent[v] = v; verts[v] = {v}; odd[v] = 1; bnd[v] = 0; mark[v] = 1; touched_v.push_back(v); roots.push_back(v); }
    std::vector<int> fused;
    for (;;) {                                                  // growth
      std::vector<int> active;
      for (int r : roots) if (parent[r] == r && odd[r] && !bnd[r]) active.push_back(r);
      if (active.empty()) break;
      fused.clear();
      for (int r : active)
        for (int v : verts[r])
          for (int e : g.adj[v]) {
            int& s = support[e];
            if (s >= g.edges[e].len) continue;
            if (s == 0) touched_e.push_back(e);
            if (++s >= g.edges[e].len) fused.push_back(e);
          }
      for (int e : fused) {
        const Edge& x = g.edges[e];
        if (x.v == B) { bnd[find(x.u)] = 1; continue; }
        join(x.u); join(x.v);
        int ru = find(x.u), rv = find(x.v);
        if (ru == rv) continue;
        if (verts[ru].size() < verts[rv].size()) std::swap(ru, rv);
        parent[rv] = ru; odd[ru] ^= odd[rv]; bnd[ru] |= bnd[rv];
        verts[ru].insert(verts[ru].end(), verts[rv].begin(), verts[rv].end()); verts[rv].clear();
      }
    }
    int flips = 0;                                              // peeling, per cluster
    for (int r : roots) {
      if (parent[r] != r || verts[r].empty()) continue;
      std::vector<int> cand;
      for (int v : verts[r])
        for (int e : g.adj[v]) {
          const Edge& x = g.edges[e];
          if (seen_e[e] || support[e] < x.len) continue;
          if (x.v != B && find(x.v) != r) continue;
          if (find(x.u) != r) continue;
          seen_e[e] = 1; cand.push_back(e);
        }
      std::stable_sort(cand.begin(), cand.end(), [&](int a, int b) { return g.edges[a].w < g.edges[b].w; });
      std::vector<int> dt;                                       // Kruskal
      for (int e : cand) {
        const Edge& x = g.edges[e];
        for (int v : {x.u, x.v}) if (dsu[v] < 0) { dsu[v] = v; dt.push_back(v); }
        int a = dfind(x.u), b = dfind(x.v);
        if (a == b) continue;
        dsu[a] = b; tree[x.u].push_back(e); tree[x.v].push_back(e);
      }
      int root = bnd[r] ? B : verts[r][0];
      // breadth-first order from the root, then peel leaves first
      std::vector<std::array<int, 3>> order;                    // vertex, edge to parent, parent
      std::vector<int> q = {root}; mark[B] = 0; vseen[root] = 1;
      size_t h = 0;
      while (h < q.size()) {
        int v = q[h++];
        for (int e : tree[v]) {
          const Edge& x = g.edges[e]; int w = x.u == v ? x.v : x.u;
          if (vseen[w]) continue;
          vseen[w] = 1; q.push_back(w); order.push_back({w, e, v});
        }
      }
      for (int v : q) vseen[v] = 0;
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int v = (*it)[0], e = (*it)[1], par = (*it)[2];
        if (!mark[v]) continue;
        mark[v] = 0; flips ^= g.edges[e].obs;
        if (par != B) mark[par] ^= 1;
      }
      for (int v : dt) { dsu[v] = -1; tree[v].clear(); }
      for (int e : cand) seen_e[e] = 0;
    }
    for (int v : touched_v) { parent[v] = -1; verts[v].clear(); odd[v] = bnd[v] = mark[v] = 0; }
    for (int e : touched_e) support[e] = 0;
    touched_v.clear(); touched_e.clear(); roots.clear(); mark[B] = 0;
    return flips;
  }
};

// ── a simple blocking queue ─────────────────────────────────────────────────────────────────────────────
struct Item { int idx; Clock::time_point t0; };
struct Queue {
  std::mutex m; std::condition_variable cv; std::deque<Item> q; bool done = false;
  void push(Item it) { { std::lock_guard<std::mutex> l(m); q.push_back(it); } cv.notify_one(); }
  bool pop(Item& it) {
    std::unique_lock<std::mutex> l(m); cv.wait(l, [&] { return !q.empty() || done; });
    if (q.empty()) return false; it = q.front(); q.pop_front(); return true;
  }
  void finish() { { std::lock_guard<std::mutex> l(m); done = true; } cv.notify_all(); }
};

static void report(const char* mode, std::vector<double> lat_us, double secs, long count, const char* extra) {
  std::sort(lat_us.begin(), lat_us.end());
  auto pct = [&](double p) { return lat_us[std::min(lat_us.size() - 1, (size_t)(p * lat_us.size()))]; };
  long over = std::count_if(lat_us.begin(), lat_us.end(), [](double x) { return x > 1000.0; });
  std::printf("{\"mode\":\"%s\",%s\"count\":%ld,\"seconds\":%.6f,\"throughput\":%.1f,\"p50_us\":%.2f,\"p99_us\":%.2f,\"p999_us\":%.2f,\"max_us\":%.2f,\"over_1ms\":%.6f}\n",
              mode, extra, count, secs, count / secs, pct(0.5), pct(0.99), pct(0.999), lat_us.back(), (double)over / lat_us.size());
}

int main(int argc, char** argv) {
  std::string dir = argv[1], mode = argv[2];
  Graph g = load_graph(dir + "/graph.txt");
  Syndromes s = load_syndromes(dir + "/syndromes.bin");
  const long N = s.fired.size();
  bool noop = argc > 0 && std::string(argv[argc - 1]) == "noop";

  if (mode == "predict") {
    Decoder d(g); std::vector<uint8_t> pred(N);
    auto t0 = Clock::now();
    for (long k = 0; k < N; k++) pred[k] = d.decode(s.fired[k]);
    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    std::ofstream(argv[3], std::ios::binary).write((char*)pred.data(), N);
    long fails = 0; for (long k = 0; k < N; k++) fails += pred[k] != s.obs[k];
    std::printf("{\"mode\":\"predict\",\"count\":%ld,\"seconds\":%.6f,\"us_per_syndrome\":%.3f,\"logical\":%.6f}\n", N, secs, secs / N * 1e6, (double)fails / N);
    return 0;
  }
  if (mode == "capacity") {
    int T = std::stoi(argv[3]); std::atomic<long> next{0}; long total = 5 * N;
    std::vector<std::thread> th; auto t0 = Clock::now();
    for (int t = 0; t < T; t++) th.emplace_back([&] { Decoder d(g); long k; while ((k = next++) < total) d.decode(s.fired[k % N]); });
    for (auto& x : th) x.join();
    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    std::printf("{\"mode\":\"capacity\",\"threads\":%d,\"count\":%ld,\"seconds\":%.6f,\"throughput\":%.1f}\n", T, total, secs, total / secs);
    return 0;
  }
  if (mode == "open" || mode == "closed") {
    // workers claim items by an atomic counter and spin until the dispatcher has published them (the BEAM's
    // schedulers also spin before sleeping): no locks, no sleeping threads, so the harness adds little
    bool open = mode == "open";
    double rate = open ? std::stod(argv[3]) : 0; long inflight = open ? 0 : std::stol(argv[3]);
    long count = std::stol(argv[4]); int W = std::stoi(argv[5]);
    std::vector<Clock::time_point> t0s(count);
    std::vector<double> lat(count);
    std::atomic<long> published{0}, claimed{0}, done{0};
    std::vector<std::thread> workers;
    for (int w = 0; w < W; w++)
      workers.emplace_back([&] {
        Decoder d(g);
        for (;;) {
          long k = claimed.fetch_add(1);
          if (k >= count) return;
          while (published.load(std::memory_order_acquire) <= k) {}
          if (!noop) d.decode(s.fired[k % N]);
          lat[k] = std::chrono::duration<double, std::micro>(Clock::now() - t0s[k]).count();
          done.fetch_add(1, std::memory_order_release);
        }
      });
    auto t0 = Clock::now();
    if (open) {
      std::mt19937_64 rng(std::stoul(argv[6])); std::exponential_distribution<double> gap(rate);
      double at = 0;
      for (long k = 0; k < count; k++) {
        at += gap(rng);
        auto when = t0 + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(at));
        while (Clock::now() < when) {}                           // spin: arrivals on time to the microsecond
        t0s[k] = when; published.store(k + 1, std::memory_order_release);
      }
    } else {
      for (long k = 0; k < count; k++) {
        while (k - done.load(std::memory_order_acquire) >= inflight) {}   // at most INFLIGHT outstanding
        t0s[k] = Clock::now(); published.store(k + 1, std::memory_order_release);
      }
    }
    while (done.load() < count) {}
    double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    for (auto& x : workers) x.join();
    char extra[160];
    std::snprintf(extra, sizeof extra, "\"workers\":%d,\"%s\":%.1f,\"noop\":%s,", W, open ? "rate" : "inflight", open ? rate : (double)inflight, noop ? "true" : "false");
    report(mode.c_str(), lat, secs, count, extra);
    return 0;
  }
  return 1;
}
