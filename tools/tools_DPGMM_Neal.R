library(Rcpp)

# t-student localización escala ----
dt_ls_log <- function(x, mu, sigma, nu){
    dt((x-mu)/sigma, df = nu, log = log) - log(sigma)
}

dt_ls <- function(x, mu, sigma, nu){
  dt((x-mu)/sigma, nu)/sigma
}



rt_ls <- function(n, mu, sigma, nu){
  rt(n, df = nu) * sigma + mu
}

pt_ls <- function(x, mu, sigma, nu){
  pt((x-mu)/sigma, nu)
}

# Segunda función: implementación manual
dt_loc_esc <- function(x, mu = 0, sigma = 1, nu = 1, log = FALSE ) {
  coef <- gamma((nu + 1) / 2) / (gamma(nu / 2) * sqrt(nu * pi) * sigma)
  z <- (x - mu) / sigma
  dens <- coef * (1 + (z^2) / nu)^(-(nu + 1) / 2)
  if(log){return (log(dens))}
  return(dens)
}


cppFunction('
double dt_ls_rcpp(double x, double mu, double sigma, double nu) {
  double z = (x - mu)/sigma;
  return R::dt(z, nu, 0)/sigma;
}')

cppFunction('
NumericVector dt_ls_rcpp_fast(NumericVector x, double mu, double sigma, double nu, bool log = false ) {
  int n = x.size();
  NumericVector dens(n);
  
  if (sigma <= 0 || nu <= 0) {
    stop("sigma y nu deben ser positivos.");
  }

  // Logaritmo del coeficiente
  double log_coef = std::lgamma((nu + 1.0) / 2.0) - 
                    std::lgamma(nu / 2.0) - 
                    0.5 * std::log(nu * M_PI) - 
                    std::log(sigma);
  
  for (int i = 0; i < n; ++i) {
    double z = (x[i] - mu) / sigma;
    double log_kernel = -0.5 * (nu + 1.0) * std::log(1.0 + (z * z) / nu);
    double log_dens = log_coef + log_kernel;
    dens[i] = log ? log_dens : std::exp(log_dens);
  }
  
  return dens;
}
')

cppFunction('
double rcpp_mean(NumericVector Y) {
  
  double sum = 0.0;
  int n = Y.size();
  
  for (int i = 0; i < n; ++i) {
    sum += Y[i];
  }
  
  return sum / n;
}
')

cppFunction('
double S_rcpp_mean(NumericVector Y, double mean_Y) {
  int n = Y.size();
  double sum_sq = 0.0;

  for (int i = 0; i < n; ++i) {
    double diff = Y[i] - mean_Y;
    sum_sq += diff * diff;
  }

  return sum_sq;
}
')

# Funciones de Sample Rcpp -----
cppFunction('
int sample_rcpp(IntegerVector x, NumericVector prob) {
  prob = prob/sum(prob);
  NumericVector prob_cumsum = cumsum(prob);
  NumericVector u = runif(2);
  int n = prob.size();

  for (int i = 0; i < n; ++i) {
    if (u[1] <= prob_cumsum[i]) {
      return x[i];
    }
  }
  return x[n - 1]; // por seguridad numérica
}
')

cppFunction(code = '
#include <random>

int sample_rcpp_stl(IntegerVector x, NumericVector prob) {
  static std::mt19937 gen(std::random_device{}());
  std::discrete_distribution<> dist(prob.begin(), prob.end());
  return x[dist(gen)];
}
', plugins = c("cpp11")
)
