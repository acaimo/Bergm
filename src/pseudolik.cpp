#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline double bergm_log1pexp(double x) {
  return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

static inline double bergm_expit(double x) {
  if (x >= 0.0) {
    double e = std::exp(-x);
    return 1.0 / (1.0 + e);
  }
  double e = std::exp(x);
  return e / (1.0 + e);
}

static inline void adjust(const double *theta, const double *W, const double *tmle,
                          const double *tpl, int p, double *out) {
  for (int a = 0; a < p; ++a) {
    double s = tpl[a];
    for (int b = 0; b < p; ++b) s += W[a + p * b] * (theta[b] - tmle[b]);
    out[a] = s;
  }
}

static inline double logpl_at(const double *X, const double *y, const double *w,
                              int n, int p, const double *t) {
  double ll = 0.0;
  for (int i = 0; i < n; ++i) {
    double eta = 0.0;
    for (int b = 0; b < p; ++b) eta += X[i + n * b] * t[b];
    ll += w[i] * (y[i] * eta - bergm_log1pexp(eta));
  }
  return ll;
}

static inline double logprior_at(const double *theta, const double *mu,
                                 const double *P, int p, std::vector<double> &tmp) {
  for (int a = 0; a < p; ++a) tmp[a] = theta[a] - mu[a];
  double q = 0.0;
  for (int a = 0; a < p; ++a) {
    double s = 0.0;
    for (int b = 0; b < p; ++b) s += P[a + p * b] * tmp[b];
    q += tmp[a] * s;
  }
  return -0.5 * q;
}

// [[Rcpp::export]]
double cpp_logpl(NumericMatrix X, NumericVector y, NumericVector w, NumericVector theta) {
  return logpl_at(X.begin(), y.begin(), w.begin(), X.nrow(), X.ncol(), theta.begin());
}

// [[Rcpp::export]]
List cpp_adjpl_batch(NumericMatrix X, NumericVector y, NumericVector w,
                     NumericMatrix Theta, NumericMatrix W,
                     NumericVector theta_mle, NumericVector theta_pl,
                     bool grad) {
  const int n = X.nrow(), p = X.ncol(), N = Theta.nrow();
  if (Theta.ncol() != p || W.nrow() != p || W.ncol() != p)
    stop("dimension mismatch in cpp_adjpl_batch");
  NumericVector ll(N);
  NumericMatrix G(grad ? N : 0, grad ? p : 0);
  std::vector<double> th(p), t(p), gt(p);
  const double *Xp = X.begin(), *yp = y.begin(), *wp = w.begin(), *Wp = W.begin();
  for (int r = 0; r < N; ++r) {
    for (int b = 0; b < p; ++b) th[b] = Theta(r, b);
    adjust(th.data(), Wp, theta_mle.begin(), theta_pl.begin(), p, t.data());
    double s = 0.0;
    if (grad) std::fill(gt.begin(), gt.end(), 0.0);
    for (int i = 0; i < n; ++i) {
      double eta = 0.0;
      for (int b = 0; b < p; ++b) eta += Xp[i + n * b] * t[b];
      s += wp[i] * (yp[i] * eta - bergm_log1pexp(eta));
      if (grad) {
        double res = wp[i] * (yp[i] - bergm_expit(eta));
        for (int b = 0; b < p; ++b) gt[b] += Xp[i + n * b] * res;
      }
    }
    ll[r] = s;
    if (grad) {
      for (int a = 0; a < p; ++a) {
        double g = 0.0;
        for (int b = 0; b < p; ++b) g += Wp[b + p * a] * gt[b];
        G(r, a) = g;
      }
    }
  }
  return List::create(_["logpl"] = ll, _["grad"] = G);
}

// [[Rcpp::export]]
List cpp_rwm_adjpl(NumericMatrix X, NumericVector y, NumericVector w,
                   NumericMatrix W, NumericVector theta_mle, NumericVector theta_pl,
                   double temperature, NumericVector prior_mean, NumericMatrix prior_prec,
                   NumericVector init, NumericMatrix L,
                   int burnin, int iters, int thin) {
  const int n = X.nrow(), p = X.ncol();
  if (thin < 1) thin = 1;
  const int nkeep = iters / thin;
  NumericMatrix out(nkeep, p);
  RNGScope scope;

  const double *Xp = X.begin(), *yp = y.begin(), *wp = w.begin();
  const double *Wp = W.begin(), *tm = theta_mle.begin(), *tp = theta_pl.begin();
  const double *mu = prior_mean.begin(), *P = prior_prec.begin(), *Lp = L.begin();

  std::vector<double> cur(init.begin(), init.end()), prop(p), t(p), z(p), tmp(p);
  adjust(cur.data(), Wp, tm, tp, p, t.data());
  double lcur = temperature * logpl_at(Xp, yp, wp, n, p, t.data()) +
                logprior_at(cur.data(), mu, P, p, tmp);
  if (!std::isfinite(lcur)) stop("initial value has non-finite log density");

  long accepted = 0;
  int kept = 0;
  const int total = burnin + iters;
  for (int it = 0; it < total; ++it) {
    if ((it & 1023) == 0) checkUserInterrupt();
    for (int b = 0; b < p; ++b) z[b] = norm_rand();
    for (int a = 0; a < p; ++a) {
      double s = cur[a];
      for (int b = 0; b <= a; ++b) s += Lp[a + p * b] * z[b];
      prop[a] = s;
    }
    adjust(prop.data(), Wp, tm, tp, p, t.data());
    double lprop = temperature * logpl_at(Xp, yp, wp, n, p, t.data()) +
                   logprior_at(prop.data(), mu, P, p, tmp);
    if (std::log(unif_rand()) < lprop - lcur) {
      cur.swap(prop);
      lcur = lprop;
      if (it >= burnin) ++accepted;
    }
    if (it >= burnin && ((it - burnin + 1) % thin) == 0 && kept < nkeep) {
      for (int b = 0; b < p; ++b) out(kept, b) = cur[b];
      ++kept;
    }
  }
  return List::create(_["draws"] = out,
                      _["accept"] = (iters > 0) ? (double)accepted / iters : NA_REAL,
                      _["last"] = NumericVector(cur.begin(), cur.end()));
}
