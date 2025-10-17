// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdint>

// Use TBB's global_control (available in oneTBB) to cap threads per call.
#ifdef RCPP_PARALLEL_USE_TBB
  #include <tbb/global_control.h>
  #include <memory>
#endif

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// Helper: compute factorials up to n using uint64_t, error if overflow
static inline std::vector<uint64_t> factorials_up_to(int n) {
  std::vector<uint64_t> fact((size_t)n + 1);
  fact[0] = 1ULL;
  for (int i = 1; i <= n; ++i) {
    if (fact[i - 1] > std::numeric_limits<uint64_t>::max() / (uint64_t)i) {
      Rcpp::stop("Factorial overflow at %d!; reduce problem size.", i);
    }
    fact[i] = fact[i - 1] * (uint64_t)i;
  }
  return fact;
}

// -----------------------------------------------------------------------------
// Helper: k-th permutation via Lehmer code
static inline void kth_perm_lehmer(int n, uint64_t k,
                                   const std::vector<uint64_t>& fact,
                                   std::vector<int>& out) {
  out.resize((size_t)n);
  std::vector<int> pool((size_t)n);
  for (int i = 0; i < n; ++i) pool[i] = i;
  for (int pos = 0; pos < n; ++pos) {
    const int remaining = n - pos - 1;
    const uint64_t base = fact[remaining];
    const size_t idx = (size_t)(remaining >= 0 ? (k / base) : 0ULL);
    if (idx >= pool.size())
      Rcpp::stop("k-th permutation index out of range (idx=%zu, pool=%zu).", idx, pool.size());
    out[pos] = pool[idx];
    pool.erase(pool.begin() + idx);
    k %= base;
  }
}

// -----------------------------------------------------------------------------
// gamma_from_perm (same core as before)
static inline double gamma_from_perm(const arma::mat& M,
                                     const int R, const int C,
                                     const int* rperm0,
                                     const int* cperm0,
                                     const double eps)
{
  std::vector<double> below((size_t)R * C);
  std::vector<double> row_below_tot(R, 0.0);

  for (int jpos = 0; jpos < C; ++jpos) {
    const int jc = cperm0[jpos];
    double suffix = 0.0;
    for (int ipos = R - 1; ipos >= 0; --ipos) {
      const int ir = rperm0[ipos];
      const size_t idx = (size_t)ipos * C + jpos;
      below[idx] = suffix;
      suffix += M(ir, jc);
      row_below_tot[ipos] += below[idx];
    }
  }

  double Csum = 0.0, Dsum = 0.0;
  for (int ipos = 0; ipos < R; ++ipos) {
    const int ir = rperm0[ipos];
    const size_t base = (size_t)ipos * C;
    double rights_suffix = 0.0;
    for (int jpos = C - 1; jpos >= 0; --jpos) {
      const int jc = cperm0[jpos];
      const double aij = M(ir, jc);
      Csum += aij * rights_suffix;
      const double lefts = row_below_tot[ipos] - rights_suffix - below[base + jpos];
      Dsum += aij * lefts;
      rights_suffix += below[base + jpos];
    }
  }
  const double den = Csum + Dsum;
  if (std::abs(den) <= eps) return std::numeric_limits<double>::quiet_NaN();
  return (Csum - Dsum) / den;
}

// -----------------------------------------------------------------------------
// Worker
struct MaxWorkerLehmer : public Worker {
  const arma::mat& Cont;
  const int R, C;
  const double eps;

  const std::vector<uint64_t> factR, factC;
  const uint64_t n_row_perms, n_col_perms;

  double best;
  uint64_t best_i, best_j;

  MaxWorkerLehmer(const arma::mat& Cont_,
                  double eps_,
                  const std::vector<uint64_t>& factR_,
                  const std::vector<uint64_t>& factC_,
                  uint64_t n_row_perms_,
                  uint64_t n_col_perms_)
    : Cont(Cont_), R(Cont_.n_rows), C(Cont_.n_cols), eps(eps_),
      factR(factR_), factC(factC_),
      n_row_perms(n_row_perms_), n_col_perms(n_col_perms_),
      best(-std::numeric_limits<double>::infinity()),
      best_i(0), best_j(0) {}

  MaxWorkerLehmer(MaxWorkerLehmer& w, Split)
    : Cont(w.Cont), R(w.R), C(w.C), eps(w.eps),
      factR(w.factR), factC(w.factC),
      n_row_perms(w.n_row_perms), n_col_perms(w.n_col_perms),
      best(-std::numeric_limits<double>::infinity()),
      best_i(0), best_j(0) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> rperm0, cperm0;
    rperm0.reserve((size_t)R);
    cperm0.reserve((size_t)C);

    for (uint64_t i = (uint64_t)begin; i < (uint64_t)end; ++i) {
      kth_perm_lehmer(R, i, factR, rperm0);
      for (uint64_t j = 0; j < n_col_perms; ++j) {
        kth_perm_lehmer(C, j, factC, cperm0);
        const double val = gamma_from_perm(Cont, R, C, rperm0.data(), cperm0.data(), eps);
        if (std::isfinite(val) && val > best) {
          best   = val;
          best_i = i;
          best_j = j;
        }
      }
    }
  }

  void join(const MaxWorkerLehmer& rhs) {
    if (rhs.best > best) {
      best   = rhs.best;
      best_i = rhs.best_i;
      best_j = rhs.best_j;
    }
  }
};

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List fill_matrix_max_tbb_lehmer_cpp(const Rcpp::NumericMatrix& ContScaled_,
                                          double eps = 1e-15,
                                          Rcpp::Nullable<Rcpp::NumericVector> max_perms_ = R_NilValue,
                                          int num_threads = -1)
{
  const arma::mat ContScaled(REAL(ContScaled_), ContScaled_.nrow(), ContScaled_.ncol(), false, true);
  const int R = ContScaled_.nrow();
  const int C = ContScaled_.ncol();
  if (R <= 0 || C <= 0) Rcpp::stop("ContScaled must have positive dimensions.");

  // Optional per-call thread cap: oneTBB only.
  #ifdef RCPP_PARALLEL_USE_TBB
    std::unique_ptr<tbb::global_control> tbb_ctrl;
    if (num_threads > 0) {
      tbb_ctrl.reset(new tbb::global_control(
        tbb::global_control::max_allowed_parallelism,
        (size_t)num_threads
      ));
    }
  #else
    if (num_threads > 0) {
      Rcpp::Rcout << "Warning: per-call thread limiting is only available with TBB; "
                     "current backend does not support it. Use RcppParallel::setThreadOptions() in R instead.\n";
    }
  #endif

  // Factorials and totals
  const std::vector<uint64_t> factR = factorials_up_to(R);
  const std::vector<uint64_t> factC = factorials_up_to(C);
  uint64_t total_row_perms = factR[(size_t)R];
  uint64_t total_col_perms = factC[(size_t)C];

  uint64_t n_row_perms = total_row_perms;
  uint64_t n_col_perms = total_col_perms;
  if (max_perms_.isNotNull()) {
    Rcpp::NumericVector mp(max_perms_.get());
    if (mp.size() >= 1 && !NumericVector::is_na(mp[0])) {
      double mr = mp[0];
      if (mr < 1) Rcpp::stop("max_perms[1] must be >= 1");
      n_row_perms = std::min<uint64_t>((uint64_t)mr, total_row_perms);
    }
    if (mp.size() >= 2 && !NumericVector::is_na(mp[1])) {
      double mc = mp[1];
      if (mc < 1) Rcpp::stop("max_perms[2] must be >= 1");
      n_col_perms = std::min<uint64_t>((uint64_t)mc, total_col_perms);
    }
  }

  MaxWorkerLehmer worker(ContScaled, eps, factR, factC, n_row_perms, n_col_perms);
  parallelReduce((std::size_t)0, (std::size_t)n_row_perms, worker);

  // Recover maximizing perms (1-based for return)
  IntegerVector row_perm(R), col_perm(C);
  std::vector<int> rperm0, cperm0;
  kth_perm_lehmer(R, worker.best_i, factR, rperm0);
  kth_perm_lehmer(C, worker.best_j, factC, cperm0);
  for (int r = 0; r < R; ++r) row_perm[r] = rperm0[r] + 1;
  for (int c = 0; c < C; ++c) col_perm[c] = cperm0[c] + 1;

  arma::mat CT_best(R, C);
  for (int i = 0; i < R; ++i) {
    const int ir = row_perm[i] - 1;
    for (int j = 0; j < C; ++j) {
      const int jc = col_perm[j] - 1;
      CT_best(i, j) = ContScaled(ir, jc);
    }
  }

  return List::create(_["value"]    = worker.best,
                      _["table"]    = wrap(CT_best),
                      _["row_perm"] = row_perm,
                      _["col_perm"] = col_perm,
                      _["n_row_perms_searched"] = Rcpp::wrap((double)worker.best_i + 1.0),
                      _["n_col_perms_searched_each"] = Rcpp::wrap((double)n_col_perms));
}
