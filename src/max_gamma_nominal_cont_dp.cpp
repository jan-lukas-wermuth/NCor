// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>     // hardware_concurrency

// Try to get TBB global_control if available (RcppParallel >= ~5.1.6 w/oneTBB)
#if __has_include(<RcppParallel/TBB.h>)
#include <RcppParallel/TBB.h>
#define HAS_TBB_GLOBAL_CONTROL 1
#else
#define HAS_TBB_GLOBAL_CONTROL 0
#endif

using namespace Rcpp;
using namespace RcppParallel;

static inline int ctz_u32(unsigned int x) { return __builtin_ctz(x); }
static inline int popcnt_u32(unsigned int x) { return __builtin_popcount(x); }

// -------------------- Build H from (x,y) with Y-tie tolerance ----------------
struct HRes {
  std::vector<std::vector<long long>> H; // K_used x K_used
  long long T{0};
  std::vector<int> old_from_new;         // map new idx -> original (0-based)
};

static HRes build_H_nominal_cont(const IntegerVector& x,
                                 const NumericVector& y,
                                 int K,
                                 double tol) {
  const int n = x.size();
  std::vector<int> good; good.reserve(n);
  for (int i = 0; i < n; ++i) {
    int xi = x[i];
    double yi = y[i];
    if (xi != NA_INTEGER && xi >= 1 && xi <= K && R_finite(yi)) good.push_back(i);
  }

  HRes res;
  const int m = (int)good.size();
  if (m == 0) return res;

  std::vector<int> present(K, 0);
  for (int idx : good) present[x[idx]-1] = 1;

  std::vector<int> old2new(K, -1);
  for (int j = 0; j < K; ++j) if (present[j]) {
    old2new[j] = (int)res.old_from_new.size();
    res.old_from_new.push_back(j);
  }
  const int K_used = (int)res.old_from_new.size();
  if (K_used == 0) return res;

  res.H.assign(K_used, std::vector<long long>(K_used, 0LL));

  // stable sort by y
  std::vector<int> ord(m);
  for (int i = 0; i < m; ++i) ord[i] = i;
  std::stable_sort(ord.begin(), ord.end(),
                   [&](int a, int b){ return y[good[a]] < y[good[b]]; });

  std::vector<long long> seen(K_used, 0LL);

  int i = 0;
  while (i < m) {
    int s = i;
    double vy = y[ good[ord[i]] ];
    int e = i;
    while (e + 1 < m && std::fabs(y[ good[ord[e+1]] ] - vy) <= tol) ++e;

    // frequency in this tie group
    std::vector<long long> grp(K_used, 0LL);
    for (int t = s; t <= e; ++t) {
      int old0 = x[ good[ ord[t] ] ] - 1;
      int ni   = old2new[old0];
      if (ni >= 0) ++grp[ni];
    }

    // cross-group pairs: seen -> grp
    for (int p = 0; p < K_used; ++p) {
      long long sp = seen[p];
      if (!sp) continue;
      for (int q = 0; q < K_used; ++q) {
        long long gq = grp[q];
        if (!gq) continue;
        res.H[p][q] += sp * gq;
      }
    }
    for (int k2 = 0; k2 < K_used; ++k2) seen[k2] += grp[k2];
    i = e + 1;
  }

  // zero diag, compute T
  long long T = 0;
  for (int p = 0; p < K_used; ++p) {
    res.H[p][p] = 0;
    for (int q = 0; q < K_used; ++q) if (p != q) T += res.H[p][q];
  }
  res.T = T;
  return res;
}

// -------------------- Parallel precompute add[j,S] ---------------------------
struct AddPrecompute : public Worker {
  const std::vector<std::vector<long long>>& H;
  int K_used;
  unsigned int N;
  long long* add; // flat buffer: size K_used * N

  AddPrecompute(const std::vector<std::vector<long long>>& H_,
                int K_used_, unsigned int N_, std::vector<long long>& addVec)
    : H(H_), K_used(K_used_), N(N_), add(addVec.data()) {}

  void operator()(std::size_t j_begin, std::size_t j_end) {
    for (std::size_t j = j_begin; j < j_end; ++j) {
      size_t base = j * N;
      add[base + 0] = 0;
      for (unsigned int S = 1u; S < N; ++S) {
        unsigned int lsb = S & -S;
        int i = ctz_u32(lsb);
        add[base + S] = add[base + (S ^ lsb)] + H[i][(int)j];
      }
    }
  }
};

// -------------------- Parallel layer DP over target masks --------------------
struct LayerDPWorker : public Worker {
  const std::vector<unsigned int>& layer_T;              // masks of size ℓ+1
  const std::vector<long long>& dp_prev;                 // dp on layer ℓ
  std::vector<long long>& dp_next;                       // write: layer ℓ+1
  std::vector<int>& parent;                              // write parent[T]
  std::vector<int>& last;                                // write last[T]
  const std::vector<std::vector<long long>>& H;          // for on-the-fly add
  const std::vector<long long>* add;                     // null if on-the-fly
  unsigned int N;

  LayerDPWorker(const std::vector<unsigned int>& layer_T_,
                const std::vector<long long>& dp_prev_,
                std::vector<long long>& dp_next_,
                std::vector<int>& parent_,
                std::vector<int>& last_,
                const std::vector<std::vector<long long>>& H_,
                const std::vector<long long>* add_,
                unsigned int N_)
    : layer_T(layer_T_), dp_prev(dp_prev_), dp_next(dp_next_),
      parent(parent_), last(last_), H(H_), add(add_), N(N_) {}

  inline long long add_on_the_fly(int j, unsigned int S) const {
    long long s = 0;
    unsigned int R = S;
    while (R) {
      unsigned int bit = R & -R;
      int i = ctz_u32(bit);
      s += H[i][j];
      R ^= bit;
    }
    return s;
  }

  void operator()(std::size_t begin, std::size_t end) {
    const long long NEG = std::numeric_limits<long long>::min()/4;
    for (std::size_t idx = begin; idx < end; ++idx) {
      unsigned int T = layer_T[idx];
      long long best = NEG;
      int bestS = -1, bestj = -1;

      unsigned int R = T;
      while (R) {
        unsigned int bit = R & -R;
        int j = ctz_u32(bit);
        unsigned int S = T ^ bit;              // predecessor on layer ℓ
        long long base = dp_prev[S];
        if (base != NEG) {
          long long inc = (add ? (*add)[(size_t)j * N + S] : add_on_the_fly(j, S));
          long long cand = base + inc;
          if (cand > best) { best = cand; bestS = (int)S; bestj = j; }
        }
        R ^= bit;
      }
      dp_next[T] = best;
      parent[T]  = bestS;
      last[T]    = bestj;
    }
  }
};

// -------------------- Exported function --------------------------------------
// [[Rcpp::export]]
List max_gamma_nominal_cont_dp_parallel(const IntegerVector& x,
                                        const NumericVector& y,
                                        int K,
                                        double tol = 0.0,
                                        int num_threads = 0) {
  // best-effort thread control (optional)
#if HAS_TBB_GLOBAL_CONTROL
  std::unique_ptr<tbb::global_control> gc;
  if (num_threads > 0)
    gc.reset(new tbb::global_control(tbb::global_control::max_allowed_parallelism,
                                     (size_t)num_threads));
#endif

  HRes hr = build_H_nominal_cont(x, y, K, tol);
  const int K_used = (int)hr.old_from_new.size();

  if (K_used == 0 || hr.T == 0) {
    int threads_used = (num_threads>0) ? num_threads
    : (int)std::max(1u, std::thread::hardware_concurrency());
    return List::create(_["value"]=NumericVector::create(NA_REAL),
                        _["order"]=IntegerVector(0),
                        _["C"]=NA_REAL, _["T"]=(double)hr.T,
                        _["K_used"]=K_used, _["method"]="subset_dp_end_agnostic_parallel",
                          _["precomputed_add"]=false, _["threads"]=threads_used);
  }

  const unsigned int N = 1u << K_used;

  // group masks by popcount
  std::vector<std::vector<unsigned int>> layers(K_used + 1);
  layers[0].push_back(0u);
  for (unsigned int S = 1u; S < N; ++S) {
    int pc = popcnt_u32(S);
    layers[pc].push_back(S);
  }

  // decide add-table precompute
  const double bytes_needed = (double)K_used * (double)N * sizeof(long long);
  const double BYTES_BUDGET = 256.0 * 1024.0 * 1024.0; // ~256 MB
  bool precompute_add = (bytes_needed <= BYTES_BUDGET);

  std::vector<long long> add;
  if (precompute_add) {
    add.assign((size_t)K_used * (size_t)N, 0LL);
    AddPrecompute worker(hr.H, K_used, N, add);
    parallelFor((size_t)0, (size_t)K_used, worker);
  }

  // DP arrays
  const long long NEG = std::numeric_limits<long long>::min()/4;
  std::vector<long long> dp_prev(N, NEG), dp_next(N, NEG);
  std::vector<int> parent(N, -1), last(N, -1);

  dp_prev[0] = 0;

  for (int s = 0; s < K_used; ++s) {
    LayerDPWorker w(layers[s+1], dp_prev, dp_next, parent, last, hr.H,
                    precompute_add ? &add : nullptr, N);
    parallelFor((size_t)0, (size_t)layers[s+1].size(), w);
    dp_prev.swap(dp_next);
  }

  unsigned int full = N - 1u;
  long long Cbest = dp_prev[full];

  // reconstruct
  std::vector<int> order_new; order_new.reserve(K_used);
  for (unsigned int S = full; S; S = (unsigned int)parent[S]) order_new.push_back(last[S]);
  std::reverse(order_new.begin(), order_new.end());

  IntegerVector order(order_new.size());
  for (int t = 0; t < (int)order_new.size(); ++t) {
    int old0 = hr.old_from_new[ order_new[t] ];
    order[t] = old0 + 1;
  }

  double gamma = (hr.T == 0) ? NA_REAL : ( (2.0 * (double)Cbest - (double)hr.T) / (double)hr.T );
  int threads_used = (num_threads>0) ? num_threads
  : (int)std::max(1u, std::thread::hardware_concurrency());

  return List::create(
    _["value"]          = gamma,
    _["order"]          = order,
    _["method"]         = "subset_dp_end_agnostic_parallel",
    _["precomputed_add"]= precompute_add,
    _["threads"]        = threads_used
  );
}
