// C++ (RcppArmadillo) implementation of the EM algorithm for the
// conditional logistic normal (CLN) model and CLN regression.
// Used by cln.mle2() and cln.reg2().
//
// Rows of y are split into the complete ones (no zeros) and the ones with
// at least one zero. The latter are grouped by their zero pattern, because
// Q, SA = Q S Q', K = S Q' SA^{-1} and Var(z | obs) depend only on the
// pattern, so they are computed once per pattern and iteration.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <map>
#include <string>
#include <vector>

using namespace arma;

namespace {

const double LOG2PI = std::log(2.0 * M_PI);

struct Pattern {
  uvec obs;        // observed (non-zero) components, 0-based, in increasing order
  int c;           // number of observed log-ratios = |obs| - 1
  mat Q;           // c x d: maps the alr (w.r.t. part D) to the alr of the subcomposition
  uvec rows;       // positions of the rows in the zero-containing block (0..n2-1)
  mat B;           // n_p x c observed log-ratios
};

struct Data {
  int n, D, d, n1, n2;
  mat full;                  // n1 x d alr of the complete rows
  uvec idx1, idx2;           // row indices (0-based) of complete / zero-containing rows
  double sly;                // sum of log(y) over all positive entries
  std::vector<Pattern> pats;
};

Data prepare(const mat& y) {
  Data dat;
  dat.n = y.n_rows;  dat.D = y.n_cols;  dat.d = dat.D - 1;
  const int D = dat.D, d = dat.d;

  std::vector<uword> i1, i2;
  dat.sly = 0.0;
  for ( int i = 0; i < dat.n; ++i ) {
    int nz = 0;
    for ( int j = 0; j < D; ++j ) {
      if ( y(i, j) > 0 ) { ++nz;  dat.sly += std::log( y(i, j) ); }
    }
    if ( nz == D ) i1.push_back(i); else i2.push_back(i);
  }
  dat.idx1 = conv_to<uvec>::from(i1);
  dat.idx2 = conv_to<uvec>::from(i2);
  dat.n1 = dat.idx1.n_elem;  dat.n2 = dat.idx2.n_elem;

  dat.full.set_size(dat.n1, d);
  for ( int i = 0; i < dat.n1; ++i ) {
    const uword r = dat.idx1[i];
    const double lD = std::log( y(r, D - 1) );
    for ( int j = 0; j < d; ++j )  dat.full(i, j) = std::log( y(r, j) ) - lD;
  }

  // group the zero-containing rows by their pattern
  std::map<std::string, int> key2id;
  std::vector< std::vector<uword> > members;
  for ( int i = 0; i < dat.n2; ++i ) {
    const uword r = dat.idx2[i];
    std::string key(D, '0');
    for ( int j = 0; j < D; ++j )  if ( y(r, j) > 0 )  key[j] = '1';
    std::map<std::string, int>::iterator it = key2id.find(key);
    int id;
    if ( it == key2id.end() ) {
      id = dat.pats.size();
      key2id[key] = id;
      Pattern p;
      std::vector<uword> ob;
      for ( int j = 0; j < D; ++j )  if ( key[j] == '1' )  ob.push_back(j);
      p.obs = conv_to<uvec>::from(ob);
      p.c = (int)p.obs.n_elem - 1;
      if ( p.c < 0 )  Rcpp::stop("A row of y has all its values equal to zero.");
      // Q = F(c) Sm F(d)' H(d)^{-1}. With the alr w.r.t. part D (z_D = 0), this
      // is simply z_{obs_k} - z_{obs_last}, k = 1, ..., c.
      p.Q.zeros(p.c, d);
      const uword last = p.obs[p.c];
      for ( int k = 0; k < p.c; ++k ) {
        if ( p.obs[k] < (uword)d )  p.Q(k, p.obs[k]) += 1.0;
        if ( last < (uword)d )      p.Q(k, last) -= 1.0;
      }
      dat.pats.push_back(p);
      members.push_back( std::vector<uword>() );
    } else id = it->second;
    members[id].push_back(i);
  }

  for ( size_t k = 0; k < dat.pats.size(); ++k ) {
    Pattern& p = dat.pats[k];
    p.rows = conv_to<uvec>::from(members[k]);
    p.B.set_size(p.rows.n_elem, p.c);
    const uword last = p.obs[p.c];
    for ( uword i = 0; i < p.rows.n_elem; ++i ) {
      const uword r = dat.idx2[ p.rows[i] ];
      const double lc = std::log( y(r, last) );
      for ( int j = 0; j < p.c; ++j )  p.B(i, j) = std::log( y(r, p.obs[j]) ) - lc;
    }
  }
  return dat;
}

// Upper Cholesky factor; stops with an informative error if S is not p.d.
mat chol_or_stop(const mat& S, const char* what) {
  mat R;
  if ( !chol(R, symmatu(S)) )  Rcpp::stop("The %s is not positive definite.", what);
  return R;
}

// log-likelihood of the complete rows with residuals res (n1 x d)
double ll_complete(const mat& res, const mat& S) {
  const int n1 = res.n_rows, d = res.n_cols;
  if ( n1 == 0 )  return 0.0;
  const mat R = chol_or_stop(S, "covariance matrix");
  const double logdet = 2.0 * sum( log( R.diag() ) );
  const mat Z = solve( trimatl( R.t() ), res.t() );   // d x n1
  return -0.5 * n1 * (d * LOG2PI + logdet) - 0.5 * accu( square(Z) );
}

// E-step for every pattern, given the means Mu2 (n2 x d) of the zero-containing rows.
// Fills Ez (n2 x d), V[k] = Var(z | obs) for pattern k, and returns the
// log-likelihood contribution of the zero-containing rows (without -sum log y).
double estep(const Data& dat, const mat& Mu2, const mat& S, mat& Ez, std::vector<mat>& V) {
  double ll = 0.0;
  for ( size_t k = 0; k < dat.pats.size(); ++k ) {
    const Pattern& p = dat.pats[k];
    const mat Mu = Mu2.rows(p.rows);                   // n_p x d
    if ( p.c == 0 ) {                                  // a single non-zero part: no information
      Ez.rows(p.rows) = Mu;
      V[k] = S;
      continue;
    }
    const mat SQt = S * p.Q.t();                       // d x c
    const mat SA = p.Q * SQt;                          // c x c
    const mat R = chol_or_stop(SA, "covariance matrix of an observed subcomposition");
    const mat Rinv = inv( trimatu(R) );
    const mat SAinv = Rinv * Rinv.t();
    const mat K = SQt * SAinv;                         // d x c
    V[k] = S - K * SQt.t();
    V[k] = 0.5 * ( V[k] + V[k].t() );
    const mat E = p.B - Mu * p.Q.t();                  // n_p x c residuals of the observed part
    Ez.rows(p.rows) = Mu + E * K.t();
    const mat Z = E * Rinv;                            // rows: R^{-T} e_i
    ll += -0.5 * p.rows.n_elem * ( p.c * LOG2PI + 2.0 * sum( log( R.diag() ) ) ) - 0.5 * accu( square(Z) );
  }
  return ll;
}

}  // namespace


// [[Rcpp::export(".cln_mle_cpp", rng = false)]]
Rcpp::List cln_mle_cpp(const arma::mat& y, const double tol, const int maxit) {
  const Data dat = prepare(y);
  const int n = dat.n, d = dat.d, n1 = dat.n1, n2 = dat.n2;
  if ( n1 < 2 )  Rcpp::stop("At least two compositions without zeros are required.");

  const rowvec mfull = sum(dat.full, 0);
  const mat sfull = dat.full.t() * dat.full;
  rowvec m = mfull / n1;
  mat S = ( sfull - n1 * m.t() * m ) / n1;           // ((n1 - 1)/n1) * var(full)

  mat Ez(n2, d, fill::zeros);
  std::vector<mat> V( dat.pats.size() );

  // loglik at the current parameters; the E-step quantities come for free
  double ll_old = ll_complete( dat.full.each_row() - m, S ) +
                  estep(dat, repmat(m, n2, 1), S, Ez, V) - dat.sly;
  double ll_new = ll_old;
  int it;
  for ( it = 1; it <= maxit; ++it ) {
    // M-step
    mat EzzSum = Ez.t() * Ez;
    for ( size_t k = 0; k < dat.pats.size(); ++k )  EzzSum += dat.pats[k].rows.n_elem * V[k];
    m = ( mfull + sum(Ez, 0) ) / n;
    S = ( sfull + EzzSum ) / n - m.t() * m;
    S = 0.5 * ( S + S.t() );
    // log-likelihood at the new parameters + E-step for the next iteration
    ll_new = ll_complete( dat.full.each_row() - m, S ) +
             estep(dat, repmat(m, n2, 1), S, Ez, V) - dat.sly;
    if ( std::abs(ll_new - ll_old) < tol )  break;
    ll_old = ll_new;
  }
  if ( it > maxit )  it = maxit;

  return Rcpp::List::create( Rcpp::Named("m") = m, Rcpp::Named("S") = S,
                             Rcpp::Named("loglik") = ll_new, Rcpp::Named("iters") = it );
}


// Louis (1982) observed information for (B, S); returns the standard errors of B.
static mat louis_se(const Data& dat, const mat& B, const mat& S, const mat& X, const mat& x2,
                    const std::vector<mat>& V, const mat& res2) {
  const int p = B.n_rows, d = B.n_cols, n = X.n_rows;
  const int nb = p * d, ns = d * (d + 1) / 2;
  const mat Sinv = inv_sympd( symmatu(S) );

  // basis of the symmetric matrices, in the order: for j, for i >= j
  std::vector<uword> bi(ns), bj(ns);
  int idx = 0;
  for ( int j = 0; j < d; ++j )  for ( int i = j; i < d; ++i ) { bi[idx] = i;  bj[idx] = j;  ++idx; }

  // complete-data information
  const mat IcBB = kron( Sinv, X.t() * X );
  mat IcSS(ns, ns);
  {
    // tr(Sinv Da Sinv Db), with Da = e_i e_j' + e_j e_i' (or e_i e_i')
    std::vector<mat> SDa(ns);
    for ( int a = 0; a < ns; ++a ) {
      mat Da(d, d, fill::zeros);
      Da(bi[a], bj[a]) = 1.0;  Da(bj[a], bi[a]) = 1.0;
      SDa[a] = Sinv * Da;
    }
    for ( int a = 0; a < ns; ++a )  for ( int b = a; b < ns; ++b )
      IcSS(a, b) = IcSS(b, a) = 0.5 * n * accu( SDa[a] % SDa[b].t() );
  }

  // variance of the complete-data score, given the observed data
  mat VarSB(nb, nb, fill::zeros), VarSS(ns, ns, fill::zeros), CovBS(nb, ns, fill::zeros);
  std::vector<mat> M(ns);
  for ( int a = 0; a < ns; ++a ) {
    mat Da(d, d, fill::zeros);
    Da(bi[a], bj[a]) = 1.0;  Da(bj[a], bi[a]) = 1.0;
    M[a] = Sinv * Da * Sinv;
  }
  for ( size_t k = 0; k < dat.pats.size(); ++k ) {
    const Pattern& pt = dat.pats[k];
    const mat& Vk = V[k];
    const mat SVS = Sinv * Vk * Sinv;
    const mat Xp = x2.rows(pt.rows);
    VarSB += kron( SVS, Xp.t() * Xp );
    // 0.5 * tr(Ma V Mb V), identical for all rows of the pattern
    std::vector<mat> MV(ns);
    for ( int a = 0; a < ns; ++a )  MV[a] = M[a] * Vk;
    const double np = pt.rows.n_elem;
    for ( int a = 0; a < ns; ++a )  for ( int b = a; b < ns; ++b ) {
      const double v = 0.5 * np * accu( MV[a] % MV[b].t() );
      VarSS(a, b) += v;  if ( b != a )  VarSS(b, a) += v;
    }
    const mat SV = Sinv * Vk;
    mat W(d, ns);
    for ( uword i = 0; i < pt.rows.n_elem; ++i ) {
      const uword r = pt.rows[i];
      const vec ri = res2.row(r).t();
      for ( int a = 0; a < ns; ++a )  W.col(a) = M[a] * ri;   // w_a = Ma r_i
      VarSS += W.t() * Vk * W;                                 // r' Ma V Mb r
      CovBS += kron( SV * W, x2.row(r).t() );                 // x_i' Ea Sinv V Mb r_i
    }
  }

  mat Iobs(nb + ns, nb + ns, fill::zeros);
  Iobs.submat(0, 0, nb - 1, nb - 1) = IcBB - VarSB;
  Iobs.submat(nb, nb, nb + ns - 1, nb + ns - 1) = IcSS - VarSS;
  Iobs.submat(0, nb, nb - 1, nb + ns - 1) = -CovBS;
  Iobs.submat(nb, 0, nb + ns - 1, nb - 1) = -CovBS.t();

  mat Vfull;
  if ( !inv(Vfull, Iobs) )  Rcpp::stop("The observed information matrix is singular.");
  mat se = reshape( sqrt( Vfull.submat(0, 0, nb - 1, nb - 1).diag() ), p, d );
  return se;
}


// [[Rcpp::export(".cln_reg_cpp", rng = false)]]
Rcpp::List cln_reg_cpp(const arma::mat& y, const arma::mat& x, const double tol,
                       const int maxit, const bool seb) {
  const Data dat = prepare(y);
  const int n = dat.n, d = dat.d, n1 = dat.n1, n2 = dat.n2;
  const int p = x.n_cols;
  if ( n1 <= p )
    Rcpp::stop("The compositions without zeros must be more than the columns of the design matrix.");

  const mat x1 = x.rows(dat.idx1), x2 = x.rows(dat.idx2);
  const mat X = join_cols(x1, x2);
  const mat& full = dat.full;

  // starting values: OLS on the complete cases
  mat B = solve( x1.t() * x1, x1.t() * full, solve_opts::likely_sympd );
  mat res1 = full - x1 * B;
  mat S = res1.t() * res1 / n1;
  S = 0.5 * ( S + S.t() );

  const mat XtX = X.t() * X;
  const mat Xtfull = x1.t() * full;

  mat Ez(n2, d, fill::zeros), Ez_prev;
  std::vector<mat> V( dat.pats.size() ), V_prev;

  double ll_old = ll_complete(res1, S) + estep(dat, x2 * B, S, Ez, V) - dat.sly;
  double ll_new = ll_old;
  int it;
  for ( it = 1; it <= maxit; ++it ) {
    // M-step
    mat EzzVarSum(d, d, fill::zeros);
    for ( size_t k = 0; k < dat.pats.size(); ++k )  EzzVarSum += dat.pats[k].rows.n_elem * V[k];
    B = solve( XtX, Xtfull + x2.t() * Ez, solve_opts::likely_sympd );
    res1 = full - x1 * B;
    const mat res2 = Ez - x2 * B;
    S = ( res1.t() * res1 + res2.t() * res2 + EzzVarSum ) / n;
    S = 0.5 * ( S + S.t() );
    // keep the E-step quantities used for this update (needed for the standard errors)
    Ez_prev = Ez;  V_prev = V;
    // log-likelihood at the new parameters + E-step for the next iteration
    ll_new = ll_complete(res1, S) + estep(dat, x2 * B, S, Ez, V) - dat.sly;
    if ( std::abs(ll_new - ll_old) < tol )  break;
    ll_old = ll_new;
  }
  if ( it > maxit )  it = maxit;

  Rcpp::List out = Rcpp::List::create( Rcpp::Named("beta") = B, Rcpp::Named("S") = S,
                                       Rcpp::Named("loglik") = ll_new, Rcpp::Named("iters") = it );
  if ( seb ) {
    if ( maxit < 1 ) { Ez_prev = Ez;  V_prev = V; }
    const mat res2 = Ez_prev - x2 * B;
    out["se.beta"] = louis_se(dat, B, S, X, x2, V_prev, res2);
  }
  return out;
}
