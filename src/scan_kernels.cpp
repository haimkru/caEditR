// Compiled kernels for caNRD_editQTL(engine = "scan") (R/caNRD_editQTL_scan.R). They reproduce the R code they replace
// (same formulas, same order of operations where it matters, same singularity rule) and only remove R-level loops,
// temporaries and repeated passes over the genotype matrix.
// Some R builds (e.g. this cluster's R 4.3.3) compile packages without any -O flag; request optimisation for these
// small numeric kernels only, on GCC (other compilers ignore this block).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("O3")
#endif
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Per-variant statistics in ONE pass over G[rows, cols] (1-based indices): number of missing values, sum of rounded
// dosages, carriers (round >= 1), non-hom-alt (round <= 1), sum and sum of squares (missing values skipped).
// round() is R's round(x, 0) (IEEE round half to even, as nearbyint).
// [[Rcpp::export(.sc_variant_stats_cpp)]]
NumericMatrix sc_variant_stats_cpp(NumericMatrix G, IntegerVector rows, IntegerVector cols) {
  const int nr = rows.size(), nc = cols.size(), ldg = G.nrow();
  std::vector<double> acc(6 * (size_t) nr, 0.0);                         // row-major accumulators (contiguous writes)
  const double* g0 = G.begin(); const int* rw = rows.begin();
  for (int j = 0; j < nc; ++j) {
    const double* col = g0 + (size_t) (cols[j] - 1) * ldg;
    double* a = acc.data();
    for (int r = 0; r < nr; ++r, a += 6) {
      const double g = col[rw[r] - 1];
      if (g != g) { a[0] += 1; continue; }                                // NaN / NA
      double gr = std::floor(g + 0.5);                                    // round half to even (as R's round(x))
      if (gr - g == 0.5 && std::fmod(gr, 2.0) != 0.0) gr -= 1.0;
      a[1] += gr; a[2] += (gr >= 1); a[3] += (gr <= 1); a[4] += g; a[5] += g * g;
    }
  }
  NumericMatrix out(nr, 6);
  for (int r = 0; r < nr; ++r) for (int k = 0; k < 6; ++k) out(r, k) = acc[6 * (size_t) r + k];
  return out;
}

// Per-variant algebra of .sc_scan_block after the BLAS cross-products.
// Mg0: nv x (K*q0), column (k*q0 + m) = sum_i w_i phi_ik g_i X0_im; Mgg: nv x nut (ut order); u: nv x K;
// M00i: q0 x q0; a0: q0; kidx: K x K (1-based index into the nut upper-triangle columns); ut: nut x 2 (1-based).
// Returns inv (nv x nut, R^-1 entries), ok, beta, se, stat, vif, mu (nv x K) and, if want_cov, cmb (nv x K x K array,
// [v, j, l] = sum_k T_k[v, j] Sinv(k, l)) and cmm (nv x nut).
// [[Rcpp::export(.sc_block_post_cpp)]]
List sc_block_post_cpp(NumericMatrix Mg0, NumericMatrix Mgg, NumericMatrix u, NumericMatrix M00i, NumericVector a0,
                       IntegerMatrix kidx, IntegerMatrix ut, bool want_cov, double piv_tol) {
  const int nv = Mg0.nrow(), K = u.ncol(), q0 = M00i.nrow(), nut = ut.nrow();
  NumericMatrix inv(nv, nut), beta(nv, K), se(nv, K), vif(nv, K), mu(nv, K), cmm(want_cov ? nv : 0, want_cov ? nut : 0);
  NumericVector stat(nv), cmb(want_cov ? nv * K * K : 0);
  LogicalVector ok(nv);
  std::vector<double> T(K * q0), Sf(nut), D(K), L(K * K), Li(K * K), Si(K * K), TS(K * K);
  for (int v = 0; v < nv; ++v) {
    // T_k[j] = sum_m Mg0[v, k*q0+m] * M00i[m, j]
    for (int k = 0; k < K; ++k) for (int j = 0; j < q0; ++j) {
      double s = 0; for (int m = 0; m < q0; ++m) s += Mg0(v, k * q0 + m) * M00i(m, j);
      T[k * q0 + j] = s;
    }
    for (int x = 0; x < nut; ++x) {
      const int k = ut(x, 0) - 1, l = ut(x, 1) - 1;
      double s = 0; for (int m = 0; m < q0; ++m) s += T[k * q0 + m] * Mg0(v, l * q0 + m);
      Sf[x] = Mgg(v, x) - s;
    }
    // correlation-scaled Cholesky with pivot tolerance (as .sc_binv)
    bool okv = true;
    for (int a = 0; a < K; ++a) { double d = Sf[kidx(a, a) - 1]; d = d > 0 ? std::sqrt(d) : 0;
      if (!(std::isfinite(d) && d > 0)) { okv = false; d = 1; } D[a] = d; }
    auto R = [&](int i, int j) { return Sf[kidx(i, j) - 1] / (D[i] * D[j]); };
    for (int j = 0; j < K; ++j) {
      double s = R(j, j); for (int k = 0; k < j; ++k) s -= L[k * K + j] * L[k * K + j];
      if (!(s > piv_tol)) { okv = false; s = 1; }
      const double ljj = std::sqrt(s); L[j * K + j] = ljj;
      for (int i = j + 1; i < K; ++i) {
        double t = R(j, i); for (int k = 0; k < j; ++k) t -= L[k * K + i] * L[k * K + j];
        L[j * K + i] = t / ljj;
      }
    }
    for (int c = 0; c < K; ++c) for (int i = c; i < K; ++i) {
      double s = (i == c) ? 1 : 0; for (int k = c; k < i; ++k) s -= L[k * K + i] * Li[c * K + k];
      Li[c * K + i] = s / L[i * K + i];
    }
    for (int a = 0; a < K; ++a) for (int b = a; b < K; ++b) {
      double s = 0; for (int i = b; i < K; ++i) s += Li[a * K + i] * Li[b * K + i];
      const double val = s / (D[a] * D[b]); Si[a * K + b] = val; Si[b * K + a] = val;
    }
    ok[v] = okv;
    for (int x = 0; x < nut; ++x) { const int a = ut(x, 0) - 1, b = ut(x, 1) - 1; inv(v, x) = okv ? Si[a * K + b] : NA_REAL; }
    if (!okv) {
      for (int k = 0; k < K; ++k) { beta(v, k) = NA_REAL; se(v, k) = NA_REAL; vif(v, k) = R_PosInf; mu(v, k) = NA_REAL; }
      stat[v] = NA_REAL;
      if (want_cov) { for (int x = 0; x < nut; ++x) cmm(v, x) = NA_REAL;
        for (int j = 0; j < K; ++j) for (int l = 0; l < K; ++l) cmb[v + nv * (j + K * l)] = NA_REAL; }
      continue;
    }
    double st = 0;
    for (int k = 0; k < K; ++k) {
      double b = 0; for (int l = 0; l < K; ++l) b += Si[k * K + l] * u(v, l);
      beta(v, k) = b; st += b * u(v, k);
      const double sd = Si[k * K + k]; se(v, k) = std::sqrt(sd > 0 ? sd : 0);
      vif(v, k) = Mgg(v, kidx(k, k) - 1) * sd;
    }
    stat[v] = st;
    for (int j = 0; j < K; ++j) {                                         // mu_j = a0_j - sum_k beta_k T_k[j]
      double a = a0[j]; for (int k = 0; k < K; ++k) a -= beta(v, k) * T[k * q0 + j];
      mu(v, j) = a;
    }
    if (want_cov) {
      for (int l = 0; l < K; ++l) for (int j = 0; j < K; ++j) {
        double s = 0; for (int k = 0; k < K; ++k) s += T[k * q0 + j] * Si[k * K + l];
        TS[l * K + j] = s; cmb[v + nv * (j + K * l)] = s;
      }
      for (int x = 0; x < nut; ++x) {
        const int j = ut(x, 0) - 1, jp = ut(x, 1) - 1;
        double s = M00i(j, jp); for (int l = 0; l < K; ++l) s += TS[l * K + j] * T[l * q0 + jp];
        cmm(v, x) = s;
      }
    }
  }
  if (want_cov) cmb.attr("dim") = IntegerVector::create(nv, K, K);
  return List::create(_["inv"] = inv, _["ok"] = ok, _["beta"] = beta, _["se"] = se, _["stat"] = stat, _["vif"] = vif,
                      _["mu"] = mu, _["cmb"] = cmb, _["cmm"] = cmm);
}
