// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>

namespace ncor_u1 {

using namespace Rcpp;
using namespace RcppParallel;

// ------------------ Lehmer decoding (factoradic -> permutation) ---------------
static inline void lehmer_decode_order(std::size_t idx,
                                       int R,
                                       const std::vector<std::size_t>& fact,
                                       std::vector<int>& out_perm) {
  // out_perm[new_row] = original_row  (0..R-1)
  out_perm.resize(R);
  std::vector<int> pool(R);
  std::iota(pool.begin(), pool.end(), 0);
  for (int pos = 0; pos < R; ++pos) {
    std::size_t f = fact[R - 1 - pos];         // factorial for remaining digits
    std::size_t q = (f ? idx / f : 0);
    idx = (f ? idx % f : 0);
    out_perm[pos] = pool[(std::size_t)q];
    pool.erase(pool.begin() + (std::size_t)q);
  }
}

// ----------------------------- Inputs bundle ----------------------------------
struct Inputs {
  int R, C;
  std::size_t P;                    // number of permutations = R!
  // counts (original, row-major r*C + c)
  std::vector<long long> cnt;
  long long n; double nd;
  // column sums & prefix (Y-order is fixed)
  std::vector<long long> colSum, colPref; // size C, C+1
  // row sums (original rows)
  std::vector<long long> rowSum;
  // weights p_{r,c} = cnt/n
  std::vector<double> w;            // length R*C
  // tie prob (invariant)
  double tie_prob;
  // factorial table for Lehmer
  std::vector<std::size_t> fact;    // length R (fact[k] = k!)
};

// Build Inputs from a contingency table
static Inputs makeInputs_lehmer(const IntegerMatrix& ContTable) {
  Inputs I;
  I.R = ContTable.nrow();
  I.C = ContTable.ncol();

  // factorials and P
  I.fact.resize(I.R);
  I.fact[0] = 1;
  for (int k = 1; k < I.R; ++k) {
    if (I.fact[k-1] > std::numeric_limits<std::size_t>::max() / (std::size_t)k)
      stop("R! overflows size_t.");
    I.fact[k] = I.fact[k-1] * (std::size_t)k;
  }
  // R! = fact[R-1] * R
  if (I.fact[I.R-1] > std::numeric_limits<std::size_t>::max() / (std::size_t)I.R)
    stop("R! overflows size_t.");
  I.P = I.fact[I.R-1] * (std::size_t)I.R;

  // counts and totals
  I.cnt.assign((size_t)I.R * (size_t)I.C, 0);
  long long n = 0;
  for (int r=0; r<I.R; ++r)
    for (int c=0; c<I.C; ++c) {
      int v = ContTable(r,c);
      if (v < 0) stop("ContTable must be nonnegative counts.");
      I.cnt[(size_t)r*I.C + c] = (long long)v;
      n += v;
    }
    if (n < 2) stop("Total sample size n must be >= 2.");
    I.n = n; I.nd = (double)n;

    // row/col sums + col prefix
    I.rowSum.assign(I.R, 0);
    I.colSum.assign(I.C, 0);
    for (int r=0; r<I.R; ++r)
      for (int c=0; c<I.C; ++c) {
        long long v = I.cnt[(size_t)r*I.C + c];
        I.rowSum[r] += v; I.colSum[c] += v;
      }
      I.colPref.assign(I.C+1, 0);
    for (int c=1; c<=I.C; ++c) I.colPref[c] = I.colPref[c-1] + I.colSum[c-1];

    // tie prob (invariant)
    double X_Tie=0.0, Y_Tie=0.0, XY_Tie=0.0;
    for (int r=0; r<I.R; ++r) { double pr=(double)I.rowSum[r]/I.nd; X_Tie += pr*pr; }
    for (int c=0; c<I.C; ++c) { double pc=(double)I.colSum[c]/I.nd; Y_Tie += pc*pc; }
    for (int r=0; r<I.R; ++r)
      for (int c=0; c<I.C; ++c) {
        double p=(double)I.cnt[(size_t)r*I.C + c]/I.nd; XY_Tie += p*p;
      }
      I.tie_prob = X_Tie + Y_Tie - XY_Tie;

    // weights
    I.w.assign((size_t)I.R*I.C, 0.0);
    for (size_t i=0; i<I.w.size(); ++i) I.w[i] = (double)I.cnt[i]/I.nd;

    return I;
}

// ------------------------- A1 helper (two entries) ----------------------------
static inline void setA1(std::vector<double>& A1, std::size_t P, std::size_t i,
                         double inv1mTie, double tau, double inv1mTie2) {
  size_t rowoff = i * (size_t)(2*P);
  A1[rowoff + (2*i + 0)] = inv1mTie;
  A1[rowoff + (2*i + 1)] = tau * inv1mTie2;
}

// ----------------------- Per-permutation worker -------------------------------
struct PermWorker : public Worker {
  const Inputs& I;
  std::vector<double>& A1vec;     // P x (2P)
  std::vector<double>& Sigvec;    // (2P) x (2P), we fill only diag block here
  std::vector<double>& Abig;      // P x (R*C)
  std::vector<double>& Bbig;      // P x (R*C)

  PermWorker(const Inputs& I_,
             std::vector<double>& A1vec_,
             std::vector<double>& Sigvec_,
             std::vector<double>& Abig_,
             std::vector<double>& Bbig_)
    : I(I_), A1vec(A1vec_), Sigvec(Sigvec_), Abig(Abig_), Bbig(Bbig_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int R = I.R, C = I.C;
    const double invn = 1.0 / I.nd;
    const double inv1mTie  = 1.0 / (1.0 - I.tie_prob);
    const double inv1mTie2 = inv1mTie * inv1mTie;

    // scratch
    std::vector<int> perm; perm.reserve(R);
    std::vector<long long> rowSumPerm(R,0), rowPref(R+1,0);
    std::vector<long long> cntPerm((size_t)R*C, 0);
    std::vector<long long> cum((size_t)(R+1)*(C+1), 0);
    std::vector<double> GX(R,0.0), GY(C,0.0);

    // GY (fixed by Y)
    for (int j=1; j<=C; ++j) {
      double pjm1=(double)I.colPref[j-1]*invn, pj=(double)I.colPref[j]*invn;
      GY[j-1] = 0.5*(pjm1 + pj);
    }

    for (std::size_t i = begin; i < end; ++i) {
      // --- decode permutation i (0..P-1) into physical row order
      lehmer_decode_order(i, R, I.fact, perm); // perm[new_row] = original_row

      // row sums & prefix in perm order
      for (int r=0; r<R; ++r) rowSumPerm[r] = I.rowSum[ perm[r] ];
      rowPref[0] = 0;
      for (int r=1; r<=R; ++r) rowPref[r] = rowPref[r-1] + rowSumPerm[r-1];
      for (int r=0; r<R; ++r) {
        double pm1=(double)rowPref[r]*invn, pr=(double)rowPref[r+1]*invn;
        GX[r] = 0.5*(pm1 + pr);
      }

      // permuted counts and 2D prefix
      for (int r=0; r<R; ++r) {
        int orow = perm[r];
        for (int c=0; c<C; ++c)
          cntPerm[(size_t)r*C + c] = I.cnt[(size_t)orow*C + c];
      }
      std::fill(cum.begin(), cum.end(), 0);
      for (int r=1; r<=R; ++r) {
        long long rs = 0;
        for (int c=1; c<=C; ++c) {
          rs += cntPerm[(size_t)(r-1)*C + (c-1)];
          cum[(size_t)r*(C+1) + c] = cum[(size_t)(r-1)*(C+1) + c] + rs;
        }
      }

      // Kendall's tau (C-D)/choose(n,2) via prefix
      long long Cnum=0, Dnum=0, total=I.n;
      for (int r=1; r<=R; ++r) {
        for (int c=1; c<=C; ++c) {
          long long v = cntPerm[(size_t)(r-1)*C + (c-1)];
          if (!v) continue;
          long long sum_row_to_r = cum[(size_t)r*(C+1) + C];
          long long sum_col_to_c = cum[(size_t)R*(C+1) + c];
          long long upto_rc      = cum[(size_t)r*(C+1) + c];

          long long NE = total - sum_row_to_r - sum_col_to_c + upto_rc;
          long long col_lt_c_all  = cum[(size_t)R*(C+1) + (c-1)];
          long long col_lt_c_to_r = cum[(size_t)r*(C+1) + (c-1)];
          long long SW = col_lt_c_all - col_lt_c_to_r;

          Cnum += v * NE;
          Dnum += v * SW;
        }
      }
      double choose2 = 0.5 * (double)I.n * ((double)I.n - 1.0);
      double tau = (choose2>0) ? ((double)Cnum - (double)Dnum) / choose2 : NA_REAL;

      // A1 row i
      setA1(A1vec, I.P, i, inv1mTie, tau, inv1mTie2);

      // Build A,B; accumulate diag moments; store A,B in ORIGINAL coords
      double meanA2=0.0, meanB2=0.0, meanAB=0.0;
      size_t offAB = i * (size_t)(I.R * I.C);
      double* Ai = Abig.data() + offAB;
      double* Bi = Bbig.data() + offAB;

      for (int r=1; r<=R; ++r) {
        double GXr = GX[r-1];
        double px_eq = (double)rowSumPerm[r-1] * invn;
        for (int c=1; c<=C; ++c) {
          double A = (double)cum[(size_t)r*(C+1)+c]*invn;
          double B = (double)cum[(size_t)r*(C+1)+(c-1)]*invn;
          double Cc= (double)cum[(size_t)(r-1)*(C+1)+c]*invn;
          double D = (double)cum[(size_t)(r-1)*(C+1)+(c-1)]*invn;
          double GXY = 0.25*(A+B+Cc+D);
          double GYr = GY[c-1];

          long long v = cntPerm[(size_t)(r-1)*C + (c-1)];
          double Ay = 4.0*GXY - 2.0*(GXr + GYr) + 1.0 - tau;
          double By = px_eq + (double)I.colSum[c-1]*invn - (double)v*invn - I.tie_prob;

          double p_rc = (double)v * invn;
          meanA2 += p_rc * (Ay*Ay);
          meanB2 += p_rc * (By*By);
          meanAB += p_rc * (Ay*By);

          int orow = perm[r-1];
          size_t oidx = (size_t)orow*I.C + (c-1);
          Ai[oidx] = Ay;
          Bi[oidx] = By;
        }
      }

      // place diagonal 2x2 block
      size_t N2 = (size_t)(2*I.P);
      size_t rAA = (size_t)(2*i),   cAA = (size_t)(2*i);
      size_t rBB = (size_t)(2*i+1), cBB = (size_t)(2*i+1);
      size_t rBA = (size_t)(2*i+1), cBA = (size_t)(2*i);
      Sigvec[rAA*N2 + cAA] = 4.0*meanA2;
      Sigvec[rBB*N2 + cBB] = 4.0*meanB2;
      Sigvec[rBA*N2 + cBA] = 4.0*meanAB; // lower-left; mirror later
    }
  }
};

// ----------------------- Pairwise lower-triangle worker -----------------------
struct PairWorker : public Worker {
  const Inputs& I;
  const std::vector<double>& Abig;
  const std::vector<double>& Bbig;
  const std::vector<double>& w;
  std::vector<double>& Sigvec;

  PairWorker(const Inputs& I_,
             const std::vector<double>& Abig_,
             const std::vector<double>& Bbig_,
             const std::vector<double>& w_,
             std::vector<double>& Sigvec_)
    : I(I_), Abig(Abig_), Bbig(Bbig_), w(w_), Sigvec(Sigvec_) {}

  void operator()(std::size_t j0_begin, std::size_t j0_end) {
    const size_t RC = (size_t)I.R * (size_t)I.C;
    const size_t N2 = (size_t)(2*I.P);

    for (std::size_t j0 = j0_begin; j0 < j0_end; ++j0) {
      const double* Aj = Abig.data() + j0 * RC;
      const double* Bj = Bbig.data() + j0 * RC;

      for (std::size_t i0 = 0; i0 < j0; ++i0) {
        const double* Ai = Abig.data() + i0 * RC;
        const double* Bi = Bbig.data() + i0 * RC;

        double sAA=0.0, sAB=0.0, sBB=0.0, sBA=0.0;
        for (size_t t=0; t<RC; ++t) {
          double wt = w[t];
          sAA += wt * (Ai[t] * Aj[t]);
          sAB += wt * (Ai[t] * Bj[t]);
          sBB += wt * (Bi[t] * Bj[t]);
          sBA += wt * (Bi[t] * Aj[t]);
        }
        sAA *= 4.0; sAB *= 4.0; sBB *= 4.0; sBA *= 4.0;

        size_t rA = (size_t)(2*j0),     rB = (size_t)(2*j0 + 1);
        size_t cA = (size_t)(2*i0),     cB = (size_t)(2*i0 + 1);

        Sigvec[rA*N2 + cA] = sAA; // (2j-1, 2i-1) 1-based
        Sigvec[rB*N2 + cA] = sAB; // (2j,   2i-1)
        Sigvec[rA*N2 + cB] = sBA; // (2j-1, 2i  )
        Sigvec[rB*N2 + cB] = sBB; // (2j,   2i  )
      }
    }
  }
};
} // namespace ncor_u1

// ---------------------------- Exported function -------------------------------
// Enumerate *all* permutations of the R rows in lexicographic order via
// Lehmer decoding; do not store any permutation matrix.
// [[Rcpp::export]]
Rcpp::List build_A1_SigmaU1_lehmer(const Rcpp::IntegerMatrix& ContTable) {

  ncor_u1::Inputs I = ncor_u1::makeInputs_lehmer(ContTable);
  const std::size_t P = I.P;
  const int R = I.R, C = I.C;

  // memory budget (still need A1, Sigma, and the A/B caches)
  const size_t A1_rows = P, A1_cols = 2*P, Sig_dim = 2*P;
  const double bytes_A1  = (double)(A1_rows*A1_cols) * sizeof(double);
  const double bytes_Sig = (double)(Sig_dim*Sig_dim) * sizeof(double);
  const double bytes_AB  = (double)P * (double)(R*C) * sizeof(double) * 2.0;
  const double BUDGET    = 10.0 * 1024.0 * 1024.0 * 1024.0; // ~10 GiB
  if (bytes_A1 + bytes_Sig + bytes_AB > BUDGET)
    Rcpp::stop("Requested dimensions too large for memory (A1 + Sigma_U1 + AB cache > ~6 GiB).");

  // flat buffers
  std::vector<double> A1vec(A1_rows*A1_cols, 0.0);
  std::vector<double> Sigvec(Sig_dim*Sig_dim, NA_REAL);
  std::vector<double> Abig((size_t)P * (size_t)(R*C), 0.0);
  std::vector<double> Bbig((size_t)P * (size_t)(R*C), 0.0);

  // per-permutation stage
  {
    ncor_u1::PermWorker w(I, A1vec, Sigvec, Abig, Bbig);
    parallelFor((size_t)0, (size_t)P, w);
  }
  // pairwise lower/upper triangle
  {
    ncor_u1::PairWorker w(I, Abig, Bbig, I.w, Sigvec);
    parallelFor((size_t)1, (size_t)P, w); // j0 = 1..P-1
  }

  // symmetrize lower -> upper
  for (size_t r = 0; r < Sig_dim; ++r) {
    if (R_IsNA(Sigvec[r*Sig_dim + r])) Sigvec[r*Sig_dim + r] = 0.0;
    for (size_t c = 0; c < r; ++c)
      Sigvec[c*Sig_dim + r] = Sigvec[r*Sig_dim + c];
  }

  // wrap into R matrices (column-major safe)
  Rcpp::NumericMatrix A1((int)A1_rows, (int)A1_cols);
  for (size_t r=0; r<A1_rows; ++r)
    for (size_t c=0; c<A1_cols; ++c)
      A1(r,c) = A1vec[r*A1_cols + c];

  Rcpp::NumericMatrix Sigma_U1((int)Sig_dim, (int)Sig_dim);
  for (size_t r=0; r<Sig_dim; ++r)
    for (size_t c=0; c<Sig_dim; ++c)
      Sigma_U1(r,c) = Sigvec[r*Sig_dim + c];

  return Rcpp::List::create(
    Rcpp::_["A1"]       = A1,
    Rcpp::_["Sigma_U1"] = Sigma_U1,
    Rcpp::_["n"]        = (double)I.n,
    Rcpp::_["tie_prob"] = I.tie_prob,
    Rcpp::_["dims"]     = Rcpp::IntegerVector::create(R, C),
    Rcpp::_["P"]        = (double)P,            // may exceed 32-bit int
    Rcpp::_["method"]   = "Lehmer on-the-fly, counts-based, TBB-parallel"
  );
}
