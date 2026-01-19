#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Logaritmo del determinante vía Cholesky
// [[Rcpp::export]]
double log_determinante(const arma::mat& A) {
  arma::mat L = arma::chol(A, "lower");
  return 2.0 * arma::sum(arma::log(L.diag()));
}

// Simulación: inversa (suponiendo definida positiva)
// [[Rcpp::export]]
arma::mat inversa_sim_def_positiva(const arma::mat& A) {
  return arma::inv_sympd(A);
}

// [[Rcpp::export]]
double dt_mvstudent_rpp_inv_compare_centered_opt(
  const arma::vec& x,
  double log_k,
  const arma::mat& Sigma_inv_candidate,
  const arma::mat& Sigma_inv_old,
  double nu,
  double a_0_b_0,
  bool logd = true) {
  
  int p = x.n_elem;
  arma::vec LcTx = arma::chol(Sigma_inv_candidate).t() * x;
  arma::vec LoTx = arma::chol(Sigma_inv_old).t() * x;
  double quad_candidate = arma::dot(LcTx, LcTx) * a_0_b_0;
  double quad_old = arma::dot(LoTx, LoTx) * a_0_b_0;
  
  double log_density = log_k - 0.5 * (nu + p) *
    (std::log1p(quad_candidate / nu) - std::log1p(quad_old / nu));
  
  return logd ? log_density : std::exp(log_density);
}
