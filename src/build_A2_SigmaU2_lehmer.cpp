// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <thread>

namespace ncor_u2 {

using namespace Rcpp;
using namespace RcppParallel;

// ---------- Lehmer decode: index -> permutation (as "order of originals") ----
static inline void lehmer_decode(std::size_t idx, int K,
                                 const std::vector<std::size_t>& fact,
                                 std::vector<int>& out) {
  out.resize(K);
  std::vector<int> pool(K);
  std::iota(pool.begin(), pool.end(), 0);         // 0..K-1 originals
  for (int pos = 0; pos < K; ++pos) {
    std::size_t f = fact[K-1-pos];                // (K-1-pos)!
    std::size_t q = (f ? idx / f : 0);
    idx = (f ? idx % f : 0);
    out[pos] = pool[(std::size_t)q];                   // new[pos] = original
    pool.erase(pool.begin() + (std::size_t)q);
  }
}

// ------------------------------ Inputs bundle --------------------------------
struct Inputs {
  int R, C;                   // #row cats, #col cats
  std::size_t P, Q, S;        // P = R!, Q = C!, S = P*Q
  std::vector<long long> cnt; // R*C, row-major in ORIGINAL coords
  long long n; double nd;
  std::vector<long long> rowSum, colSum;
  std::vector<double> w;      // weights p_rc = cnt/n in ORIGINAL coords
  double tie_prob;            // invariant
  std::vector<std::size_t> factR, factC; // factorial tables
};

// Build Inputs and invariants from contingency table (counts)
static Inputs makeInputs(const IntegerMatrix& CT) {
  Inputs I;
  I.R = CT.nrow();
  I.C = CT.ncol();

  // factorials & number of permutations
  I.factR.resize(I.R); I.factC.resize(I.C);
  I.factR[0] = 1; for (int k=1;k<I.R;++k){
    if (I.factR[k-1] > std::numeric_limits<std::size_t>::max() / (std::size_t)k)
      stop("R! overflows size_t");
    I.factR[k] = I.factR[k-1]*(size_t)k;
  }
  I.factC[0] = 1; for (int k=1;k<I.C;++k){
    if (I.factC[k-1] > std::numeric_limits<std::size_t>::max() / (std::size_t)k)
      stop("C! overflows size_t");
    I.factC[k] = I.factC[k-1]*(size_t)k;
  }
  if (I.factR[I.R-1] > std::numeric_limits<std::size_t>::max()/(std::size_t)I.R)
    stop("R! overflows size_t");
  if (I.factC[I.C-1] > std::numeric_limits<std::size_t>::max()/(std::size_t)I.C)
    stop("C! overflows size_t");
  I.P = I.factR[I.R-1]*(size_t)I.R;
  I.Q = I.factC[I.C-1]*(size_t)I.C;
  if (I.P > std::numeric_limits<std::size_t>::max()/I.Q)
    stop("R!*C! overflows size_t");
  I.S = I.P * I.Q;

  // counts, totals, row & col sums
  I.cnt.assign((size_t)I.R*(size_t)I.C, 0);
  long long n=0;
  I.rowSum.assign(I.R,0); I.colSum.assign(I.C,0);
  for (int r=0;r<I.R;++r){
    for (int c=0;c<I.C;++c){
      int v = CT(r,c);
      if (v<0) stop("ContTable must be nonnegative counts.");
      I.cnt[(size_t)r*I.C + c] = (long long)v;
      I.rowSum[r]+=v; I.colSum[c]+=v; n+=v;
    }
  }
  if (n<2) stop("Total sample size n must be >= 2.");
  I.n = n; I.nd = (double)n;

  // tie probability (invariant)
  double X2=0.0, Y2=0.0, XY2=0.0;
  for (int r=0;r<I.R;++r){ double pr=(double)I.rowSum[r]/I.nd; X2+=pr*pr; }
  for (int c=0;c<I.C;++c){ double pc=(double)I.colSum[c]/I.nd; Y2+=pc*pc; }
  for (int r=0;r<I.R;++r) for (int c=0;c<I.C;++c){
    double p=(double)I.cnt[(size_t)r*I.C+c]/I.nd; XY2+=p*p;
  }
  I.tie_prob = X2 + Y2 - XY2;

  // weights in ORIGINAL coords
  I.w.assign((size_t)I.R*(size_t)I.C, 0.0);
  for (size_t t=0;t<I.w.size();++t) I.w[t] = (double)I.cnt[t]/I.nd;

  return I;
}

// ---------------------------- A2 row helper ----------------------------------
static inline void setA2(std::vector<double>& A2buf, std::size_t S, std::size_t s,
                         double inv1mTie, double tau, double inv1mTie2) {
  size_t rowoff = s * (size_t)(2*S);
  A2buf[rowoff + (2*s + 0)] = inv1mTie;
  A2buf[rowoff + (2*s + 1)] = tau * inv1mTie2;
}

// ----------------------- Per-(rowPerm,colPerm) worker ------------------------
struct ComboWorker : public Worker {
  const Inputs& I;

  // outputs
  std::vector<double>& A2buf;     // S x (2S), row-major
  std::vector<double>& Sigbuf;    // (2S) x (2S), row-major; here we write only diagonal 2x2 blocks
  std::vector<double>& Abig;      // S x (R*C), row-major; A per ORIGINAL (r,c)
  std::vector<double>& Bbig;      // S x (R*C), row-major; B per ORIGINAL (r,c)
  const std::size_t N2;
  const std::size_t RC;

  ComboWorker(const Inputs& I_,
              std::vector<double>& A2buf_,
              std::vector<double>& Sigbuf_,
              std::vector<double>& Abig_,
              std::vector<double>& Bbig_,
              std::size_t N2_, std::size_t RC_)
    : I(I_), A2buf(A2buf_), Sigbuf(Sigbuf_), Abig(Abig_), Bbig(Bbig_), N2(N2_), RC(RC_) {}

  void operator()(std::size_t s_begin, std::size_t s_end) {
    const int R=I.R, C=I.C;
    const double invn = 1.0/I.nd;
    const double inv1m = 1.0/(1.0 - I.tie_prob);
    const double inv1m2 = inv1m*inv1m;

    // scratch
    std::vector<int> permR, permC; permR.reserve(R); permC.reserve(C);
    std::vector<int> posR(R), posC(C);
    std::vector<long long> rowSumPerm(R), colSumPerm(C);
    std::vector<long long> rowPref(R+1), colPref(C+1);
    std::vector<double> GX(R), GY(C);

    std::vector<long long> cntPerm((size_t)R*(size_t)C, 0);
    std::vector<long long> cum((size_t)(R+1)*(size_t)(C+1), 0);

    for (std::size_t s = s_begin; s < s_end; ++s) {
      // decode s -> (i,j)
      std::size_t i = s / I.Q;      // row permutation index
      std::size_t j = s % I.Q;      // col permutation index
      lehmer_decode(i, R, I.factR, permR);  // permR[new] = original row
      lehmer_decode(j, C, I.factC, permC);  // permC[new] = original col

      for (int r=0;r<R;++r) posR[ permR[r] ] = r;
      for (int c=0;c<C;++c) posC[ permC[c] ] = c;

      // row/col sums and mid-CDFs in permuted order
      for (int r=0;r<R;++r) rowSumPerm[r] = I.rowSum[ permR[r] ];
      for (int c=0;c<C;++c) colSumPerm[c] = I.colSum[ permC[c] ];
      rowPref[0]=0; for (int r=1;r<=R;++r) rowPref[r] = rowPref[r-1] + rowSumPerm[r-1];
      colPref[0]=0; for (int c=1;c<=C;++c) colPref[c] = colPref[c-1] + colSumPerm[c-1];
      for (int r=0;r<R;++r) { GX[r] = 0.5*((double)rowPref[r]*invn + (double)rowPref[r+1]*invn); }
      for (int c=0;c<C;++c) { GY[c] = 0.5*((double)colPref[c]*invn + (double)colPref[c+1]*invn); }

      // counts in permuted grid (r',c') and its 2D prefix
      for (int r=0;r<R;++r){
        int orow = permR[r];
        for (int c=0;c<C;++c){
          int ocol = permC[c];
          cntPerm[(size_t)r*C + c] = I.cnt[(size_t)orow*I.C + ocol];
        }
      }
      std::fill(cum.begin(), cum.end(), 0);
      for (int r=1;r<=R;++r){
        long long rs=0;
        for (int c=1;c<=C;++c){
          rs += cntPerm[(size_t)(r-1)*C + (c-1)];
          cum[(size_t)r*(C+1) + c] = cum[(size_t)(r-1)*(C+1) + c] + rs;
        }
      }

      // Kendall's tau via C-D count on permuted grid
      long long Cnum=0, Dnum=0, total=I.n;
      for (int r=1;r<=R;++r){
        for (int c=1;c<=C;++c){
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

      // A2 row (two entries)
      setA2(A2buf, I.S, s, inv1m, tau, inv1m2);

      // Build A,B at ORIGINAL (ro,co), and diagonal moments
      double meanA2=0.0, meanB2=0.0, meanAB=0.0;
      size_t off = s * RC;
      double* As = Abig.data() + off;
      double* Bs = Bbig.data() + off;

      for (int ro=0; ro<R; ++ro){
        int r = posR[ro]; double GXr = GX[r];
        double px_eq = (double)I.rowSum[ro]*invn;
        for (int co=0; co<C; ++co){
          int c = posC[co];
          // mid-CDF GXY at perm coords (r,c)
          double A = (double)cum[(size_t)(r+1)*(C+1) + (c+1)] * invn;
          double B = (double)cum[(size_t)(r+1)*(C+1) + (c  )] * invn;
          double Cc= (double)cum[(size_t)(r  )*(C+1) + (c+1)] * invn;
          double D = (double)cum[(size_t)(r  )*(C+1) + (c  )] * invn;
          double GXY = 0.25 * (A + B + Cc + D);

          double GYr = GY[c];
          double y_eq = (double)I.colSum[co]*invn;

          long long v = cntPerm[(size_t)r*C + c];
          double p_rc = (double)v * invn;           // = cnt/n
          double Ay = 4.0*GXY - 2.0*(GXr + GYr) + 1.0 - tau;
          double By = px_eq + y_eq - p_rc - I.tie_prob;

          meanA2 += p_rc * (Ay*Ay);
          meanB2 += p_rc * (By*By);
          meanAB += p_rc * (Ay*By);

          size_t oidx = (size_t)ro*I.C + co;
          As[oidx] = Ay;
          Bs[oidx] = By;
        }
      }

      // place diagonal 2x2 block into Sigbuf (row-major)
      size_t rAA = (size_t)(2*s),     cAA = (size_t)(2*s);
      size_t rBB = (size_t)(2*s + 1), cBB = (size_t)(2*s + 1);
      size_t rBA = (size_t)(2*s + 1), cBA = (size_t)(2*s);
      Sigbuf[rAA*N2 + cAA] = 4.0*meanA2;
      Sigbuf[rBB*N2 + cBB] = 4.0*meanB2;
      Sigbuf[rBA*N2 + cBA] = 4.0*meanAB;  // lower-left; mirror later
    }
  }
};

// --------------------- Pairwise lower-triangle worker ------------------------
struct PairWorker : public Worker {
  const Inputs& I;
  const std::vector<double>& Abig;
  const std::vector<double>& Bbig;
  const std::vector<double>& w;  // ORIGINAL weights p_rc
  std::vector<double>& Sigbuf;
  const std::size_t N2;
  const std::size_t RC;

  PairWorker(const Inputs& I_,
             const std::vector<double>& Abig_,
             const std::vector<double>& Bbig_,
             const std::vector<double>& w_,
             std::vector<double>& Sigbuf_,
             std::size_t N2_, std::size_t RC_)
    : I(I_), Abig(Abig_), Bbig(Bbig_), w(w_), Sigbuf(Sigbuf_), N2(N2_), RC(RC_) {}

  void operator()(std::size_t j0_begin, std::size_t j0_end) {

    for (std::size_t j0 = j0_begin; j0 < j0_end; ++j0) {   // second index (row block), 0..S-1
      const double* Aj = Abig.data() + j0 * RC;
      const double* Bj = Bbig.data() + j0 * RC;

      for (std::size_t i0 = 0; i0 < j0; ++i0) {            // first index (col block), 0..j0-1
        const double* Ai = Abig.data() + i0 * RC;
        const double* Bi = Bbig.data() + i0 * RC;

        double sAA=0.0, sAB=0.0, sBB=0.0, sBA=0.0;
        for (size_t t=0; t<RC; ++t) {
          double wt = w[t];
          sAA += wt * (Ai[t] * Aj[t]);   // E[A_i A_j]
          sAB += wt * (Ai[t] * Bj[t]);   // E[A_i B_j]
          sBB += wt * (Bi[t] * Bj[t]);   // E[B_i B_j]
          sBA += wt * (Bi[t] * Aj[t]);   // E[B_i A_j]
        }
        sAA *= 4.0; sAB *= 4.0; sBB *= 4.0; sBA *= 4.0;

        // Write the 2x2 lower-triangle block for (j0,i0)
        size_t rA = (size_t)(2*j0),     rB = (size_t)(2*j0 + 1); // rows of A_j, B_j
        size_t cA = (size_t)(2*i0),     cB = (size_t)(2*i0 + 1); // cols of A_i, B_i

        Sigbuf[rA*N2 + cA] = sAA;  // (2j-1, 2i-1)
        Sigbuf[rB*N2 + cA] = sAB;  // (2j,   2i-1)
        Sigbuf[rA*N2 + cB] = sBA;  // (2j-1, 2i  )
        Sigbuf[rB*N2 + cB] = sBB;  // (2j,   2i  )
      }
    }
  }
};
} // namespace ncor_u2

// --------------------------------- Export ------------------------------------
// Build A2 and Sigma_U2 for two nominal variables (rows & columns permuted),
// enumerating all permutations on the fly via Lehmer.
// [[Rcpp::export]]
Rcpp::List build_A2_SigmaU2_lehmer(const Rcpp::IntegerMatrix& ContTable) {

  ncor_u2::Inputs I = ncor_u2::makeInputs(ContTable);
  const size_t S = I.S;      // total (rowPerm, colPerm) pairs
  const int R = I.R, C = I.C;
  const std::size_t N2 = 2 * S;
  const size_t RC = (size_t)R*(size_t)C;

  // Memory guard (Σ is 2S x 2S; A/B caches are S x (R*C))
  const size_t A2_rows = S, A2_cols = 2*S, Sig_dim = 2*S;
  const double bytes_A2  = (double)(A2_rows*A2_cols) * sizeof(double);
  const double bytes_Sig = (double)(Sig_dim*Sig_dim) * sizeof(double);
  const double bytes_AB  = (double)S * (double)(R*C) * sizeof(double) * 2.0;
  const double BUDGET    = 10.0 * 1024.0 * 1024.0 * 1024.0; // ~10 GiB
  if (bytes_A2 + bytes_Sig + bytes_AB > BUDGET) {
    Rcpp::stop("Requested dimensions too large for memory (A2 + Sigma_U2 + AB cache > ~10 GiB).\n"
           "Reduce the number of categories or use a sub-sampling of permutations.");
  }

  // Flat buffers (row-major for convenience)
  std::vector<double> A2buf(A2_rows*A2_cols, 0.0);
  std::vector<double> Sigbuf(Sig_dim*Sig_dim, NA_REAL);
  std::vector<double> Abig((size_t)S * (size_t)(R*C), 0.0);
  std::vector<double> Bbig((size_t)S * (size_t)(R*C), 0.0);

  // Stage 1: per-(rowPerm, colPerm) work (tau, A2 row, diagonal Σ block, and A/B cache)
  {
    ncor_u2::ComboWorker w(I, A2buf, Sigbuf, Abig, Bbig, N2, RC);
    parallelFor((size_t)0, (size_t)S, w);
  }

  // Stage 2: pairwise lower/upper triangle of Σ_U2
  {
    ncor_u2::PairWorker w(I, Abig, Bbig, I.w, Sigbuf, N2, RC);
    parallelFor((size_t)1, (size_t)S, w);  // j0 = 1..S-1 (0-based); j0=0 would be a no-op
  }

  // Symmetrize Σ: copy lower -> upper, set diag to 0 if still NA (shouldn't be)
  for (size_t r=0; r<Sig_dim; ++r) {
    if (R_IsNA(Sigbuf[r*Sig_dim + r])) Sigbuf[r*Sig_dim + r] = 0.0;
    for (size_t c=0; c<r; ++c)
      Sigbuf[c*Sig_dim + r] = Sigbuf[r*Sig_dim + c];
  }

  // Wrap into R matrices (column-major safe)
  Rcpp::NumericMatrix A2((int)A2_rows, (int)A2_cols);
  for (size_t r=0; r<A2_rows; ++r)
    for (size_t c=0; c<A2_cols; ++c)
      A2(r,c) = A2buf[r*A2_cols + c];

  Rcpp::NumericMatrix Sigma_U2((int)Sig_dim, (int)Sig_dim);
  for (size_t r=0; r<Sig_dim; ++r)
    for (size_t c=0; c<Sig_dim; ++c)
      Sigma_U2(r,c) = Sigbuf[r*Sig_dim + c];

  return Rcpp::List::create(
    Rcpp::_["A2"]        = A2,
    Rcpp::_["Sigma_U2"]  = Sigma_U2,
    Rcpp::_["n"]         = (double)I.n,
    Rcpp::_["tie_prob"]  = I.tie_prob,
    Rcpp::_["dims"]      = Rcpp::IntegerVector::create(R, C),
    Rcpp::_["P_rows"]    = (double)I.P,
    Rcpp::_["Q_cols"]    = (double)I.Q,
    Rcpp::_["S_pairs"]   = (double)I.S,
    Rcpp::_["method"]    = "Lehmer on-the-fly (rows & cols), counts-based, TBB-parallel"
  );
}

