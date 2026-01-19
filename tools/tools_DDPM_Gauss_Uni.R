
library(progress)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

## stick_breaking_cpp -----


cppFunction('

NumericVector stick_breaking_cpp(NumericVector v) {
  int n = v.size();
  NumericVector proportions(n);
  double remaining_stick = 1.0;

  for(int i = 0; i < n; ++i) {
    proportions[i] = remaining_stick * v[i];
    remaining_stick *= (1 - v[i]);
  }

  return proportions;
}
')

## z_sample_rnorm_cpp ----
cppFunction('
  IntegerVector z_sample_rnorm_XB_cpp(NumericVector Y ,NumericVector w, NumericMatrix X, NumericMatrix Matrix_beta, NumericVector sigma2_vec, int H) {
    int n = Y.size();
    int p = X.ncol();
    NumericMatrix prob_zij(n, H);  // Matriz de probabilidades
    NumericMatrix cum_probs(n, H); // Matriz para las CDF acumuladas
    double productoPunto = 0.0;
    // Llenar la matriz de probabilidades prob_zij
    for (int j = 0; j < H; j++) {
      for (int i = 0; i < n; i++) {
      productoPunto = 0.0;
      for (int m = 0; m < p; m++) {
          productoPunto += X(i,m) *Matrix_beta(m,j);
      }
      
        prob_zij(i, j) = w[j] * R::dnorm(Y[i], productoPunto, std::sqrt(sigma2_vec[j]), false) + 1e-100;
      }
    }
  
    // Evitar problemas numéricos
    for (int i = 0; i < n; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < H; ++j) {
        row_sum += prob_zij(i, j);
      }
      for (int j = 0; j < H; ++j) {
        prob_zij(i, j) /= row_sum;
      }
    }
  NumericVector random_vals = runif(n);
    IntegerVector z_sampled(n);
    
    for (int i = 0; i < n; ++i) {
      double cumsum = 0.0;
      for (int j = 0; j < H; ++j) {
        cumsum += prob_zij(i, j);
        if (random_vals[i] <= cumsum) {
          z_sampled[i] = j + 1;
          break;
        }
      }
    }
    return z_sampled;
  }
')

cppFunction(depends = "RcppArmadillo", plugins = "cpp11", code = '
arma::vec mat_vec_mult(const arma::mat& A, const arma::vec& x) {
  return A * x;
}
')

# Función 'inversa_rapida Matrix' con RcppArmadillo ----
# Define la función 'inversa_rapida' con RcppArmadillo
cppFunction(depends = "RcppArmadillo",  plugins = c("cpp11"),code = '
  arma::mat inversa_rapida(const arma::mat& A) {
    return inv(A);  // Inversa general (para cualquier matriz cuadrada invertible)
  }
')
cppFunction(depends = "RcppArmadillo",plugins = c("cpp11"), code = '
  arma::mat inversa_sim_def_positiva(const arma::mat& A) {
    return inv_sympd(A);  // Solo para matrices simétricas definidas positivas
  }
')

cppFunction(depends = "RcppArmadillo",plugins = c("cpp11"), code = '
  double log_determinante(const arma::mat& A) {
  if (A.n_rows != A.n_cols) {
    stop("La matriz debe ser cuadrada.");
  }

  // Descomposición de Cholesky: A = L * Lᵗ
  arma::mat L = arma::chol(A, "lower");

  // log(|A|) = 2 * sum(log(diag(L)))
  return 2.0 * arma::sum(arma::log(L.diag()));
}
')

dt_mvstudent_manual <- function(x, mu, Sigma, nu) {
  p <- length(x)
  Sigma_inv <- solve(Sigma)
  det_Sigma <- det(Sigma)
  
  # Constante de normalización
  log_k <- lgamma((nu + p)/2) - lgamma(nu/2) - (p/2)*log(nu*pi) - 0.5*log(det_Sigma)
  k <- exp(log_k)
  
  # Cálculo de la densidad
  cuadratica <- as.numeric(t(x - mu) %*% Sigma_inv %*% (x- mu))
  densidad <- k * (1 + (1/nu) * cuadratica)^(-(nu + p)/2)
  
  return(densidad)
}

# Densidad T-Multivariada con RcppArmadillo ----
cppFunction(depends = "RcppArmadillo",  plugins = c("cpp11"),code ='
    double dt_mvstudent_rcpp(const arma::vec& x, const arma::vec& mu, 
    const arma::mat& Sigma, double nu, bool logd = false) {
    
    int p = x.n_elem;
    const arma::vec x_centered = x - mu;
    // Calcula la inversa y determinante de Sigma (forma estable)
    arma::mat Sigma_inv;
    double det_Sigma;
    bool success = arma::inv_sympd(Sigma_inv, Sigma);
    
    if (!success) {
        Rcpp::stop("Sigma no es definida positiva");
    }
    det_Sigma = arma::det(Sigma);
    
    // Constante de normalización
    double log_k = R::lgammafn((nu + p)* 0.5) - R::lgammafn(nu* 0.5) - 
                  (p* 0.5) * std::log(nu * M_PI) - 0.5 * std::log(det_Sigma);

    
    // Forma cuadrática eficiente
    double quadratic = arma::dot(x_centered, Sigma_inv * x_centered);
    
   
    // Densidad (log o escala original)
    if (logd) {
        return log_k - ((nu + p)* 0.5) * std::log1p( quadratic/nu);
    } else {
        return std::exp(log_k) * std::pow(1.0 + quadratic/nu, -(nu + p)* 0.5);
    }
}
')



# Densidad T-Multivariada_inv centrada en 0 con RcppArmadillo ----
cppFunction(depends = "RcppArmadillo",  plugins = c("cpp11"),code ='
    double dt_mvstudent_rcpp_inv_centered(const arma::vec& x,
    const arma::mat& Sigma_inv, double nu, bool logd = false) {
    
    int p = x.n_elem;
    
    // Calcula la inversa y determinante de Sigma (forma estable)
    double det_Sigma_inv;
    
    arma::mat L_Sigma_inv = arma::chol(Sigma_inv, "lower");
    
    // Constante de normalización
    double log_k = R::lgammafn((nu + p)*0.5) - R::lgammafn(nu*0.5) - 
                  (p*0.5) * std::log(nu * M_PI) +  arma::sum(arma::log(L_Sigma_inv.diag()));

    
    // Forma cuadrática eficiente
    arma::vec LcTx = L_Sigma_inv.t() * x;
    double quadratic = arma::dot(LcTx, LcTx);
    
    // Densidad (log o escala original)
   double log_density = log_k - (nu + p)*0.5 * std::log1p(quadratic/nu);

    return logd ? log_density : std::exp(log_density);
    
}
')

# Densidad T-Multivariada_inv comparacion log(t_candidate) - log(t_old) con RcppArmadillo ----
cppFunction(
  depends = "RcppArmadillo", plugins = "cpp11", code = '
double dt_mvstudent_rpp_inv_compare_centered_opt(
    const arma::vec& x_candidate,
    const arma::vec& x_old,
    const arma::mat& Sigma_inv_candidate,
    const arma::mat& Sigma_inv_old,
    double log_k,
    double a_0_b_0,
    double nu,
    bool logd = false) {

    int p = x_candidate.n_elem;

    // Cuadráticas
    double quad_candidate = arma::as_scalar(x_candidate.t() * (Sigma_inv_candidate * x_candidate)) * a_0_b_0;
    double quad_old = arma::as_scalar(x_old.t() * (Sigma_inv_old * x_old)) * a_0_b_0;

    double log_density = log_k - 0.5 * (nu + p) * 
        (std::log1p(quad_candidate / nu) - std::log1p(quad_old / nu)); // log1p(x) = log(1 + x)

    return logd ? log_density : std::exp(log_density);
}
')


cppFunction(depends = "RcppArmadillo", plugins = "cpp11", code = '
double dt_mvstudent_rpp_inv_compare_centered_chol(
    const arma::vec& x,
    const arma::mat& Sigma_inv_candidate,
    const arma::mat& Sigma_inv_old,
    double a_0_b_0,
    double nu,
    bool logd = false) {

    int p = x.n_elem;

    // Descomposición de Cholesky
    arma::mat L_candidate = arma::chol(Sigma_inv_candidate, "lower");
    arma::mat L_old = arma::chol(Sigma_inv_old, "lower");

    // log-determinante: log(|A|) = 2 * sum(log(diag(L)))
    // log_k = 0.5*log(det(Sigma_inv_candidate)) - log(det(Sigma_inv_old))
    double log_k =  arma::sum(arma::log(L_candidate.diag())) - arma::sum(arma::log(L_old.diag()));

    // Forma cuadrática: || Lᵗ x ||²
    arma::vec LcTx = L_candidate.t() * x;
    arma::vec LoTx = L_old.t() * x;

    double quad_candidate = arma::dot(LcTx, LcTx) * a_0_b_0;
    double quad_old = arma::dot(LoTx, LoTx) * a_0_b_0;

    double log_density = log_k - 0.5 * (nu + p) *
        (std::log1p(quad_candidate / nu) - std::log1p(quad_old / nu));

    return logd ? log_density : std::exp(log_density);
}
')


# Forma Cuadratica Xi*Sigma*Xi^t
cppFunction(depends = "RcppArmadillo",plugins = c("cpp11"), code = '
  Rcpp::NumericVector quad_xi_sigma_xi(const arma::vec& x, const arma::cube& Sigma) {
    int C = Sigma.n_slices;
    Rcpp::NumericVector result(C);

    for (int c = 0; c < C; ++c) {
      result(c) = arma::as_scalar(x.t() * Sigma.slice(c) * x);
    }

    return result;
  }
')





# Densidad de la mezcla de normales univariado -------
densidad_mezcla_normales_uni_X <- function(x, X_new , betas, sigmas2, proporciones) {
  # Validar entradas
  k <- dim(betas)[2]
  if (length(sigmas2) != k | length(proporciones) != k) {
    stop("Las listas de medias, sigmas y proporciones deben tener la misma longitud")
  }
  if (abs(sum(proporciones) - 1) > 1e-5) {
    stop("Las proporciones deben sumar 1")
  }
  
  dens <- 0
  
  # Sumar las densidades ponderadas de cada componente
  for (i in 1:k) {
    dens <- dens + proporciones[i] * dnorm(x, mean = X_new %*%betas[,i], sd = sqrt(sigmas2[i]))
  }
  
  return(dens)
}


