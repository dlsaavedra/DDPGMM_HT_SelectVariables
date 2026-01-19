# Función para realizar el stick breaking
library(MASS)
library(tidyverse)
library(progress)
library(Rcpp)


# Densidad de la mezcla de normales univariado -------
densidad_mezcla_normales_uni <- function(x, medias, sigmas2, proporciones) {
  # Validar entradas
  k <- length(medias)
  if (length(sigmas2) != k | length(proporciones) != k) {
    stop("Las listas de medias, sigmas y proporciones deben tener la misma longitud")
  }
  if (abs(sum(proporciones) - 1) > 1e-5) {
    stop("Las proporciones deben sumar 1")
  }
  
  dens <- 0
  
  # Sumar las densidades ponderadas de cada componente
  for (i in 1:k) {
    dens <- dens + proporciones[i] * dnorm(x, mean = medias[[i]], sd = sqrt(sigmas2[[i]]))
  }
  
  return(dens)
}

generar_mezcla_normales_uni <- function(n, medias, sigmas2, proporciones) {
  # Validar entradas
  k <- ncol(medias)
  if (length(sigmas2[1,1,]) != k | length(proporciones) != k) {
    stop("Las listas de medias, sigmas y proporciones deben tener la misma longitud")
  }
  if (abs(sum(proporciones) -1) > 1e-10) {
    stop("Las proporciones deben sumar 1")
  }
  
  # Generar muestras para cada componente
  muestras <- list()
  indicadora <- c()
  for (i in 1:k) {
    ni <- round(n * proporciones[i])
    muestras[[i]] <- rnorm(ni, medias[i], sqrt(sigmas2[i]))
    indicadora <- c(indicadora, rep(i,ni))
  }
  # Combinar todas las muestras
  todas_muestras <- do.call(rbind, muestras)
  return(list(muestra = todas_muestras, ind = indicadora))
}


# Tools DPM ----

suma_matrices = function(lista_matrices){
  suma <- matrix(0, nrow = nrow(lista_matrices[[1]]), ncol = ncol(lista_matrices[[1]]))
  # Sumar cada matriz en la lista
  for (mat in lista_matrices) {
    suma <- suma + mat
  }
  return(suma)
}

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
IntegerVector z_sample_rnorm_XB_cpp(NumericVector Y ,NumericVector w, NumericVector Xt_beta_vec, NumericVector sigma2_vec, int H) {
  int n = Y.size();
  NumericMatrix prob_zij(n, H);  // Matriz de probabilidades
  NumericMatrix cum_probs(n, H); // Matriz para las CDF acumuladas

  // Llenar la matriz de probabilidades prob_zij
  for (int j = 0; j < H; ++j) {
    for (int i = 0; i < n; ++i) {
      prob_zij(i, j) = w[j] * R::dnorm(Y[i], Xt_beta_vec[j], std::sqrt(sigma2_vec[j]), false) + 1e-100;
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

  // Calcular las CDF acumuladas
  for (int i = 0; i < n; ++i) {
    cum_probs(i, 0) = prob_zij(i, 0);
    for (int j = 1; j < H; ++j) {
      cum_probs(i, j) = cum_probs(i, j - 1) + prob_zij(i, j);
    }
  }

  // Generar valores aleatorios
  NumericVector random_vals = runif(n);

  // Encontrar la primera categoría cuyo valor acumulado exceda al valor aleatorio
  IntegerVector z_sampled(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < H; ++j) {
      if (random_vals[i] <= cum_probs(i, j)) {
        z_sampled[i] = j + 1;
        break;
      }
    }
  }

  return z_sampled;
}')



# DPM Gaussian Univariado con covariables Y = X^tB + error,  error = Mezcla gaussiana ------
DPM_Gaussian_Uni <- function(Y, X, param_chain, param_init ){
  
  # Priori NIG(mu, Sigma| mu_0, lambda,Psi, eta) = N(mu|mu_0, 1/lambda Sigma) * IG(Sigma|  eta/2, Psi/2)
  ## Parametros MCMC ----
  burnit = param_chain$burnit
  n_chain = param_chain$n_chain
  thining = param_chain$thining
  Psi_0 = param_init$Psi_0
  Psi_inv_0 = solve(Psi_0)
  mu_0 = param_init$mu_0
  H = param_init$H
  sigma2_0 = param_init$sigma2_0
  eta_0 = param_init$eta_0
  alpha = param_init$alpha
  p = ncol(X)
  n = nrow(Y)
  Xt = t(X)
  
  # Crear las matrices de los parámetros a guardar ----
  Matrix = list()
  Matrix$beta = array(dim = c(H, n_chain))
  Matrix$beta = matrix(rnorm(p*n_chain), nrow = p, ncol = n_chain)
  Matrix$Sigma2 = array(dim = c(H, n_chain))
  Matrix$v = array(dim = c(H, n_chain))
  Matrix$v[H, ] = 1
  Matrix$w = array(dim = c(H, n_chain))
  Matrix$y_pred = array(dim = c(1, n_chain))
  Matrix$y_pred = matrix(rnorm(p*n_chain), nrow = 1, ncol = n_chain)
  Matrix$z_pred = array(dim = c(n_chain))
  frame_categorias = data.frame(categoria = 1:H)
  D = diag(n)
  
  ## Puntos iniciales ----
  Sigma2_array_inicial <- replicate(H, 1/rgamma(1, shape = eta_0, rate  = Psi_0))
  #array(rep(Sigma_inicial, ncol = 2), H), dim = c(2, 2, 3))
  v_inicial = 
  w_inicial = 
  
  Matrix$beta[,1] = Matrix$beta[,1] * chol(Psi_0) +  mu_0
  Matrix$Sigma2[,1] = Sigma2_array_inicial
  Matrix$v[, 1] = c(rbeta(H - 1,1,alpha),1)
  Matrix$w[, 1] = stick_breaking_cpp(v_inicial)
  
  ## Variables inicializadas al comienzo ----
  sum_yj <- matrix(0, nrow = length(1:H), ncol = ncol(Y))
  rownames(sum_yj) <- 1:H
  z = rep(0,n)
  nj = rep(0,H)
  cum_nj = rep(0,H)
  matrix_Y_muj = 0
  matrix_muj_mo = 0
  
  # Inicializar la barra de progreso
  #pb <- txtProgressBar(min = 0, max = n_chain, style = 3)
  # Crear la barra de progreso
  pb <- progress_bar$new(
    format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
    total = n_chain/100,
    clear = FALSE,
    width = 60
  )
  for (iter in 1:(n_chain-1)){
    
    # Vector indicadora de la categoria de la observación i ----
    z = z_sample_rnorm_XB_cpp(Y, w = Matrix$w[, iter],
                           Xt_beta_vec = X%*%Matrix$beta[, iter],
                           sigma2_vec = Matrix$Sigma2[, iter], H)
    # Suma del número de componentes
    nj = as.vector(table(factor(z,levels = 1:H)))
    # Suma acumulada
    cum_nj = sum(nj) - cumsum(nj)
    
    # Actualización v y w ----
    #Matrix$v[, iter + 1] = c(mapply(function(x, y) rbeta(1, 1 + x, alpha + y), nj[-H], cum_nj[-H]), 1)
    for (i in 1:(H-1)) {
      Matrix$v[i, iter + 1] <- rbeta(1, 1 + nj[i], alpha + cum_nj[i])
    }
    Matrix$w[, iter + 1] = stick_breaking_cpp(Matrix$v[, iter + 1])
    
    # Actualizamos los betas
    
    Matrix$beta[, iter + 1] = Matrix$beta[, iter + 1] * chol(solve(Xt%*%D%*%X + Psi_inv_0))
    
    #Sigma (como cada componente es independiente acualizaremos mu_j y sigma_j) -----
    
    
    matrix_muj_mo = (mu_0 - Matrix$mu[, iter + 1])^2
    for (j in 1:H){
      # Calculos para actualizar el Sigma
      if (nj[j] == 0){
        matrix_Y_muj = 0
      }else{
        matrix_Y_muj = sum((Y[z ==j] - Matrix$mu[j,iter + 1])^2)
        
      }
      #matrix_muj_mo = (lambda * nj[j])/(lambda + nj[j])* (mean(y[z == j]) - mu_0)^2 # Chat GPT
      
      Matrix$Sigma[j, iter + 1] = 1/rgamma(1, shape = (nj[j] + eta + 1)/2, rate = (matrix_Y_muj + (lambda*nj[j])/(nj[j] + lambda)*matrix_muj_mo[j] + Psi)/2)
    }
    ## Predicción -----
    
    Matrix$z_pred[iter + 1] = sample(1:H, 1, prob =  Matrix$w[, iter + 1])
    Matrix$y_pred[iter + 1] = Matrix$y_pred[iter + 1]* sqrt(Matrix$Sigma[Matrix$z_pred[iter + 1], iter + 1]) + Matrix$mu[Matrix$z_pred[iter + 1]]
    
    
    # Actualizar la barra de progreso
    #setTxtProgressBar(pb, iter)
    if (iter %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
      pb$tick()
    }
  }
  # Cerrar la barra de progreso al completar
  #close(pb)
  
  ## Burnint y Thining ----
  secuencia = seq(burnit + 1, n_chain, by = thining)
  Matrix$mu = Matrix$mu[, secuencia]
  Matrix$Sigma =  Matrix$Sigma[, secuencia]
  Matrix$v = Matrix$v[, secuencia]
  Matrix$w =  Matrix$w[, secuencia]
  Matrix$y_pred =  Matrix$y_pred[secuencia]
  Matrix$z_pred =  Matrix$z_pred[secuencia]
  
  return(Matrix)
}


# DPM Gaussian Univariado priori alpha ------
DPM_Gaussian_Uni_EW <- function(Y, param_chain, param_init ){
  
  # Priori NIG(mu, Sigma| mu_0, lambda,Psi, eta) = N(mu|mu_0, 1/lambda Sigma) * IG(Sigma|  eta/2, Psi/2)
  ## Parametros MCMC ----
  burnit = param_chain$burnit
  n_chain = param_chain$n_chain
  thining = param_chain$thining
  Psi_0 = param_init$Psi_0
  Psi_inv_0 = solve(Psi_0)
  mu_0 = param_init$mu_0
  H = param_init$H
  sigma2_0 = param_init$sigma2_0
  eta_0 = param_init$eta_0
  alpha = param_init$alpha
  p = ncol(X)
  
  # Crear las matrices de los parámetros a guardar ----
  Matrix = list()
  Matrix$beta = array(dim = c(H, n_chain))
  Matrix$beta = matrix(rnorm(p*n_chain), nrow = p, ncol = n_chain)
  Matrix$Sigma2 = array(dim = c(H, n_chain))
  Matrix$v = array(dim = c(H, n_chain))
  Matrix$v[H, ] = 1
  Matrix$w = array(dim = c(H, n_chain))
  Matrix$y_pred = array(dim = c(1, n_chain))
  Matrix$y_pred = matrix(rnorm(p*n_chain), nrow = 1, ncol = n_chain)
  Matrix$z_pred = array(dim = c(n_chain))
  frame_categorias = data.frame(categoria = 1:H)
  
  ## Puntos iniciales ----
  Sigma2_array_inicial <- replicate(H, 1/rgamma(1, shape = eta_0, rate  = Psi_0))
  #array(rep(Sigma_inicial, ncol = 2), H), dim = c(2, 2, 3))
  v_inicial = 
    w_inicial = 
    
    Matrix$beta[,1] = Matrix$beta[,1] * chol(Psi_0) +  mu_0
  Matrix$Sigma2[,1] = Sigma2_array_inicial
  Matrix$v[, 1] = c(rbeta(H - 1,1,alpha),1)
  Matrix$w[, 1] = stick_breaking_cpp(v_inicial)
  
  ## Variables inicializadas al comienzo ----
  sum_yj <- matrix(0, nrow = length(1:H), ncol = ncol(Y))
  rownames(sum_yj) <- 1:H
  z = rep(0,n)
  nj = rep(0,H)
  cum_nj = rep(0,H)
  matrix_Y_muj = 0
  matrix_muj_mo = 0
  
  # Inicializar la barra de progreso
  #pb <- txtProgressBar(min = 0, max = n_chain, style = 3)
  # Crear la barra de progreso
  pb <- progress_bar$new(
    format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
    total = n_chain/100,
    clear = FALSE,
    width = 60
  )
  for (iter in 1:(n_chain-1)){
    
    # Vector indicadora de la categoria de la observación i ----
    z = z_sample_rnorm_XB_cpp(Y, w = Matrix$w[, iter],
                           Xt_beta_vec = X%*%Matrix$beta[, iter],
                           sigma2_vec = Matrix$Sigma2[, iter], H)
    # Suma del número de componentes
    nj = as.vector(table(factor(z,levels = 1:H)))
    # Suma acumulada
    cum_nj = sum(nj) - cumsum(nj)
    
    # Actualización v y w ----
    #Matrix$v[, iter + 1] = c(mapply(function(x, y) rbeta(1, 1 + x, Matrix$alpha[iter] + y), nj[-H], cum_nj[-H]), 1)
    for (i in 1:(H-1)) {
      Matrix$v[i, iter + 1] <- rbeta(1, 1 + nj[i], Matrix$alpha[iter] + cum_nj[i])
    }
    Matrix$w[, iter + 1] = stick_breaking_cpp(Matrix$v[, iter + 1])
    
    # Actualizamos de nu y alpha ----
    # Actualizar alpha (Escobar y West 1995)
    Matrix$nu[iter + 1] <- rbeta(1, Matrix$alpha[iter] + 1, n)
    
    Matrix$alpha[iter + 1] <- rgamma(1,shape = a + H - rbinom(1,1,(n * (b - log(Matrix$nu[iter + 1]))) / (n * (b - log(Matrix$nu[iter + 1])) + a + H - 1)), 
                                     rate = b - log(Matrix$nu[iter + 1]))
    
    
    # Actualizamos mu y Sigma (como cada componente es independiente acualizaremos mu_j y sigma_j)-----
    for (j in 1:H){
      # Actualizamos el mu
      
      #Matrix$mu[j, iter + 1] = rnorm(1,mean = (lambda * mu_0 + sum(Y[z == j]))/(nj[j] + lambda),sd = sqrt(1/(nj[j] + lambda)* Matrix$Sigma[j, iter]))
      Matrix$mu[j, iter + 1]  = Matrix$mu[j, iter + 1] * sqrt(1/(nj[j] + lambda)* Matrix$Sigma[j, iter]) + (lambda * mu_0 + sum(Y[z == j]))/(nj[j] + lambda)
    }
    
    matrix_muj_mo = (mu_0 - Matrix$mu[, iter + 1])^2
    for (j in 1:H){
      # Calculos para actualizar el Sigma
      if (nj[j] == 0){
        matrix_Y_muj = 0
      }else{
        matrix_Y_muj = sum((Y[z ==j] - Matrix$mu[j,iter + 1])^2)
        
      }
      #matrix_muj_mo = (lambda * nj[j])/(lambda + nj[j])* (mean(y[z == j]) - mu_0)^2 # Chat GPT
      
      Matrix$Sigma[j, iter + 1] = 1/rgamma(1, shape = (nj[j] + eta + 1)/2, rate = (matrix_Y_muj + (lambda*nj[j])/(nj[j] + lambda)*matrix_muj_mo[j] + Psi)/2)
    }
    
    
    ## Predicción -----
    
    Matrix$z_pred[iter + 1] = sample(1:H, 1, prob =  Matrix$w[, iter + 1])
    Matrix$y_pred[iter + 1] = Matrix$y_pred[iter + 1]* sqrt(Matrix$Sigma[Matrix$z_pred[iter + 1], iter + 1]) + Matrix$mu[Matrix$z_pred[iter + 1]]
    
    
    # Actualizar la barra de progreso
    #setTxtProgressBar(pb, iter)
    if (iter %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
      pb$tick()
    }
  }
  # Cerrar la barra de progreso al completar
  #close(pb)
  
  ## Burnint y Thining ----
  secuencia = seq(burnit + 1, n_chain, by = thining)
  Matrix$mu = Matrix$mu[, secuencia]
  Matrix$Sigma =  Matrix$Sigma[, secuencia]
  Matrix$v = Matrix$v[, secuencia]
  Matrix$w =  Matrix$w[, secuencia]
  Matrix$alpha =  Matrix$alpha[secuencia]
  Matrix$nu =  Matrix$nu[secuencia]
  Matrix$y_pred =  Matrix$y_pred[secuencia]
  Matrix$z_pred =  Matrix$z_pred[secuencia]
  
  return(Matrix)
}

