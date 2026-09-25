source("tools/WomackPrior.R")
source("tools/tools_DPGMM_Neal.R")
source("tools/tools_DDPM_Gauss_Uni.R")
library(progress)

DDPM_Gaussian_Uni_Neal_Selection <- function(Y, X, param_chain, param_init, 
                                   n_muestra = FALSE, Y_scale = T, stream = F){
  
  ## Parametros MCMC ----
  burnit = param_chain$burnit
  n_chain = param_chain$n_chain
  thining = param_chain$thining
  
  ## Hiperparámetros ----
  p = dim(X)[2]
  mu_0 = param_init$mu_0 # 0
  tau = param_init$tau # 10
  a_sigma2 = param_init$a_sigma2 #2 #2
  b_sigma2 = param_init$b_sigma2 #.5 # E(sigma^2)= b/(a-1)
  a_alpha = param_init$a_alpha #.1
  b_alpha = param_init$b_alpha #.1
  
  n = length(Y)
  
  mu_0_tau = mu_0/tau
  nu_0 = 2*a_sigma2 
  zeta = param_init$zeta
  gamma_inicio = param_init$gamma_inicio
  tau_inv = 1/tau
  tau_inv_matrix_tot = diag(tau_inv, p)
  tau_matrix_tot = diag(tau,p)
  ones_matrix = diag(1, n)
  
  womack_prior = womack(p -1, zeta)
  
  H_max = param_init$H_max #50
  n_cluster_ini = param_init$n_cluster_ini #5
  
  
  # Crear las matrices de los parámetros a guardar ----
  Matrix = list()
  Matrix$param_init = param_init
  Matrix$alpha = array(dim = c(n_chain))
  Matrix$H = array(dim = c(n_chain))
  Matrix$ind_cluster = array(dim = c(n, n_chain))
  Matrix_prob_cluster <- rep(0,H_max)
  Matrix$gamma = array(dim = c(p, n_chain))
  
  # Puntos iniciales
  Matrix$gamma[,1] = gamma_inicio
  Matrix$alpha[1] = rgamma(1, shape = a_alpha, rate = b_alpha)
  nu <- rbeta(1, Matrix$alpha[1] + 1, n)
  
  # Asignar clusters basados en cuantiles
  ind_cluster <- cut(Y, breaks = quantile(Y,
                                          probs = seq(0, 1,length.out = n_cluster_ini + 1)), 
                     labels = FALSE, include.lowest = TRUE)
  Matrix$ind_cluster[,1] = ind_cluster
  Matrix$H[1] <- max(ind_cluster)
  
  # Escalar los datos
  if (Y_scale){
    Matrix$mean_Y_or = mean(Y)
    Matrix$sd_Y_or = sd(Y)
    Y = (Y - Matrix$mean_Y_or)/Matrix$sd_Y_or
  }else{
    Matrix$mean_Y_or = 0
    Matrix$sd_Y_or = 1
  }
  Matrix$Y = Y
  Matrix$X = X
  Xt = t(X)
  
  nj = tabulate(factor(ind_cluster, levels = 1:Matrix$H[1]) , 
                nbins = Matrix$H[1])
  
  Y_t_c <- rep(0, H_max)
  aux_Sigma_t_c  <- array(dim = c(p, p, H_max))
  Sigma_t_c  <- array(dim = c(p, p, H_max))
  XY_c <- array(dim = c(p, H_max))
  mu_t_c  <- array(dim = c(p, H_max))
  a_t_c = a_sigma2 + nj/2
  b_t_c = rep(0, H_max)
  change_gamma = T
  dic_dt_ls_Y_0 <- new.env(hash = TRUE, parent = emptyenv())
  # Inicializar la barra de progreso
  # Crear la barra de progreso
  if (stream){
    pb <- progress_bar$new(
      format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
      total = n_chain/100,
      clear = FALSE,
      width = 60
    )
  }
  
  for (iter in 1:(n_chain-1)){
    
    H_iter = Matrix$H[iter]
    n_denom = 1/(n-1+Matrix$alpha[iter])
    alpha_iter_n_denom = Matrix$alpha[iter] * n_denom
    
    if( change_gamma){
      index_gamma_1 = Matrix$gamma[, iter] == 1
      index_gamma_0 = Matrix$gamma[, iter] == 0
      p_activos = sum(index_gamma_1)
      X_gamma = X[,index_gamma_1, drop = F]
      mu_0_2_tau = crossprod(mu_0[index_gamma_1]) * tau_inv
      XY_c[index_gamma_1,1:Matrix$H[iter]] = t(rowsum(X_gamma*Y, group = ind_cluster)) + mu_0_tau[index_gamma_1] # X^tY + mu_0/tau
      #prob cada dato de ser de la medida base
      if(toString(index_gamma_1) %in% ls(dic_dt_ls_Y_0)){
        dt_ls_Y_0 = dic_dt_ls_Y_0[[toString(index_gamma_1)]]
      }
      else{
        dt_ls_Y_0 = dt_ls(x = Y, mu = X_gamma%*%mu_0[index_gamma_1], sigma = sqrt((b_sigma2/a_sigma2)*(1 + tau*diag(tcrossprod(X_gamma)))), nu = nu_0) 
        dic_dt_ls_Y_0[[toString(index_gamma_1)]] = dt_ls_Y_0
      }
      for (c in 1:Matrix$H[iter]){
        idx = ind_cluster == c
        aux_Sigma_t_c[index_gamma_1,index_gamma_1,c] = crossprod(X_gamma[idx,, drop = F]) + tau_inv_matrix_tot[1:p_activos, 1:p_activos]
        Sigma_t_c[index_gamma_1,index_gamma_1,c] = inversa_sim_def_positiva(aux_Sigma_t_c[,,c][index_gamma_1,index_gamma_1, drop = F])
        mu_t_c[index_gamma_1,c] = Sigma_t_c[index_gamma_1,index_gamma_1,c] %*% XY_c[index_gamma_1, c]
        Y_t_c[c] =crossprod(Y[idx])
        b_t_c[c] = b_sigma2 + 0.5 * (Y_t_c[c] + mu_0_2_tau - crossprod(mu_t_c[index_gamma_1,c], XY_c[index_gamma_1,c]))
      
        }
      
      }
    
    ## Actualizamos el cluster para cada dato ----
    for (i in 1:n){
      
      X_i = X_gamma[i,]
      Y_i = Y[i]
      Matrix_prob_cluster[1:(H_iter + 1)] = 0
      idx_i = ind_cluster[i]
      c = (1:H_iter)[-idx_i]
      Matrix_prob_cluster[c] = nj[c]*n_denom * dt_ls(x = Y_i, mu = X_i %*%mu_t_c[index_gamma_1,c], 
                                                     #sigma = sqrt(b_t_c[c]/a_t_c[c] * (1 + apply(Sigma_t_c[index_gamma_1,index_gamma_1,c, drop = F],3, function (x) X_i%*% x%*%X_i))) , 
                                                     sigma = sqrt(b_t_c[c]/a_t_c[c] * (1 + quad_xi_sigma_xi(X_i, Sigma_t_c[index_gamma_1,index_gamma_1,c, drop = F]))) ,
                                                     nu = 2*a_t_c[c])
      
        
      if (nj[idx_i] > 1){
        
        n_ic = nj[idx_i] - 1
        aux_Sigma_ic = aux_Sigma_t_c[index_gamma_1,index_gamma_1,idx_i] - tcrossprod(X_i)
        Sigma_ic = inversa_sim_def_positiva(aux_Sigma_ic)
        XY_ic = XY_c[index_gamma_1, idx_i] - X_i*Y_i
        mu_ic = Sigma_ic %*% XY_ic
        a_ic= a_t_c[idx_i] - 0.5 
        Y_ic = Y_t_c[idx_i] - Y_i**2
        b_ic = b_sigma2 + 0.5 *  (Y_ic + mu_0_2_tau - crossprod(mu_ic, XY_ic))
        
        Matrix_prob_cluster[idx_i] = n_ic * n_denom * dt_ls_rcpp(x = Y_i, mu = X_i%*%mu_ic, sigma = sqrt(b_ic/a_ic * (1 + X_i%*% Sigma_ic%*%X_i)), nu = 2*a_ic)
        Matrix_prob_cluster[H_iter + 1] =  alpha_iter_n_denom * dt_ls_Y_0[i] 
        aux_cluster =  sample_rcpp_stl(1:(H_iter + 1) ,prob = Matrix_prob_cluster[1:(H_iter + 1)]) 
      }else{
        
        # un cluster se quedo solo con un elemento
        Matrix_prob_cluster[idx_i] = alpha_iter_n_denom * dt_ls_Y_0[i] 
        
        n_ic = 0
        aux_Sigma_ic = tau_inv_matrix_tot[1:p_activos, 1:p_activos]
        Sigma_ic = tau_matrix_tot[1:p_activos, 1:p_activos]
        XY_ic = mu_0_tau[index_gamma_1]
        mu_ic = mu_0[index_gamma_1]
        a_ic= a_sigma2 
        Y_ic = 0
        b_ic = b_sigma2
        aux_cluster =  sample_rcpp_stl(1:(H_iter ) ,prob = Matrix_prob_cluster[1:(H_iter)]) 
      }
      
      
      
      # Si cambia de cluster calculamos los nuevos estadisticos para ambos cluster que cambiaron (el nuevo y antiguo)
      if (aux_cluster != idx_i){
        
        # cluster antiguo
        nj[idx_i] = n_ic
        aux_Sigma_t_c[index_gamma_1,index_gamma_1, idx_i] = aux_Sigma_ic
        Sigma_t_c[index_gamma_1,index_gamma_1, idx_i] = Sigma_ic
        XY_c[index_gamma_1, idx_i] = XY_ic 
        mu_t_c[index_gamma_1, idx_i] = mu_ic
        a_t_c[idx_i] = a_ic
        Y_t_c[idx_i] = Y_ic
        b_t_c[idx_i] = b_ic
        
        
        if(aux_cluster == H_iter + 1){
          
          if(H_iter + 1 > H_max){

            aux_Sigma_t_c_aux  <- array(dim = c(p, p, 2*H_max))
            aux_Sigma_t_c_aux[,,1:H_max] <- aux_Sigma_t_c
            Sigma_t_c_aux  <- array(dim = c(p, p, 2*H_max))
            Sigma_t_c_aux[,,1:H_max] <- Sigma_t_c
            XY_c_aux <- array(dim = c(p, 2*H_max))
            XY_c_aux[,1:H_max] <- XY_c
            mu_t_c_aux  <- array(dim = c(p, 2*H_max))
            mu_t_c_aux[,1:H_max] <- mu_t_c
            
            aux_Sigma_t_c = aux_Sigma_t_c_aux
            Sigma_t_c = Sigma_t_c_aux
            XY_c = XY_c_aux 
            mu_t_c = mu_t_c_aux
            
            H_max = 2*H_max
            
            
          }
          
          nj[aux_cluster] = 1
          aux_Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster] = tcrossprod(X_i) + tau_inv_matrix_tot[1:p_activos, 1:p_activos]
          Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster]  = inversa_sim_def_positiva(aux_Sigma_t_c[,,aux_cluster][index_gamma_1,index_gamma_1, drop = F])
          XY_c[index_gamma_1, aux_cluster] = X_i * Y_i + mu_0_tau[index_gamma_1]
          mu_t_c[index_gamma_1, aux_cluster] =  Sigma_t_c[index_gamma_1,index_gamma_1,aux_cluster] %*% XY_c[index_gamma_1, aux_cluster]
          a_t_c[aux_cluster] = a_sigma2 + 0.5
          Y_t_c[aux_cluster] = Y_i**2 
          #b_t_c[aux_cluster] = b_sigma2 + 0.5* max(0, Y_i**2 + mu_0_2_tau - crossprod(mu_t_c[, aux_cluster], XY_c[, aux_cluster]))
          b_t_c[aux_cluster] = b_sigma2 + 0.5* (Y_t_c[aux_cluster] + mu_0_2_tau - crossprod(mu_t_c[index_gamma_1, aux_cluster], XY_c[index_gamma_1, aux_cluster]))
 
          H_iter =  H_iter + 1
          
        }else{
          
          nj[aux_cluster] = nj[aux_cluster] + 1
          
          aux_Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster] = aux_Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster] + tcrossprod(X_i) 
          Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster] = inversa_sim_def_positiva(aux_Sigma_t_c[,,aux_cluster][index_gamma_1,index_gamma_1, drop = F])
          XY_c[index_gamma_1, aux_cluster] = XY_c[index_gamma_1, aux_cluster] + X_i * Y_i
          mu_t_c[index_gamma_1, aux_cluster] = Sigma_t_c[index_gamma_1,index_gamma_1, aux_cluster] %*% XY_c[index_gamma_1, aux_cluster]
          a_t_c[aux_cluster]= a_t_c[aux_cluster] + 0.5
          Y_t_c[aux_cluster] = Y_t_c[aux_cluster] + Y_i**2 
          b_t_c[aux_cluster] = b_sigma2 + 0.5 * (Y_t_c[aux_cluster] + mu_0_2_tau - crossprod(mu_t_c[index_gamma_1, aux_cluster], XY_c[index_gamma_1, aux_cluster]))
          
        }
        
        ind_cluster[i] = aux_cluster
        
        
        # Reasignar nuevas etiquetas consecutivas
        if (nj[idx_i] == 0){
          
          # Encontrar etiquetas únicas ordenadas
          etiquetas_unicas <- sort(unique(ind_cluster))
          ind_cluster <- match(ind_cluster, etiquetas_unicas)
          
          # Reetiqueto todos los demás cluster
          nj = nj[etiquetas_unicas]
          
          aux_Sigma_t_c[,, 1:length(etiquetas_unicas)] <-  aux_Sigma_t_c[,, etiquetas_unicas]
          Sigma_t_c[,,1:length(etiquetas_unicas)] <- Sigma_t_c[,,etiquetas_unicas]
          XY_c[, 1:length(etiquetas_unicas)] <- XY_c[, etiquetas_unicas]
          mu_t_c[, 1:length(etiquetas_unicas)] <- mu_t_c[,etiquetas_unicas]
          a_t_c[1:length(etiquetas_unicas)] <- a_t_c[etiquetas_unicas]
          b_t_c[1:length(etiquetas_unicas)] <- b_t_c[etiquetas_unicas]
          Y_t_c[1:length(etiquetas_unicas)] <- Y_t_c[etiquetas_unicas]
          
          
          H_iter =  H_iter - 1
        }
        
        #Matrix$H[iter + 1] <- max(ind_cluster)
      }
      
    }
    Matrix$H[iter + 1] = H_iter
    
    # Guardar los nuevos cluster
    Matrix$ind_cluster[, iter + 1] = ind_cluster
    
    
    # Actualizamos de nu y alpha ----
    # Actualizar alpha (Escobar y West 1995)
    nu <- rbeta(1, Matrix$alpha[iter] + 1, n)
    
    Matrix$alpha[iter + 1] <- rgamma(1,shape = a_alpha + Matrix$H[iter + 1] - rbinom(1,1,(n * (b_alpha - log(nu))) / (n * (b_alpha - log(nu)) + a_alpha + Matrix$H[iter + 1] - 1)), 
                                     rate = b_alpha - log(nu))
    
    
    # Actualizar gamma ----
    
    candidate_gamma = Matrix$gamma[, iter]
    index_change_gamma = sample.int(p-1, size = 1) + 1
    candidate_gamma[index_change_gamma] = 1 - candidate_gamma[index_change_gamma]
    index_gamma_candidate_1 = candidate_gamma == 1
    p_activos_candidate = sum(index_gamma_candidate_1)
    #womack_prior[p_activos - 1] # p(gamma)
    #womack_prior[sum(candidate_gamma) - 1] # p(gamma_new)
    
    A = log(womack_prior[p_activos_candidate]) - log(womack_prior[p_activos])
    
    aux_A = (p_activos_candidate - p_activos)*log(tau)
    
    for (h in 1:H_iter){
      
      idx = ind_cluster == h
 
      #Sigma_Y_old_inv = inversa_sim_def_positiva(diag_nj + tau * XtX_old)
      Sigma_Y_old_inv = (ones_matrix[1:nj[h], 1:nj[h]] -  X[idx, index_gamma_1 , drop  = FALSE] %*% Sigma_t_c[index_gamma_1,index_gamma_1, h] %*% Xt[index_gamma_1, idx, drop  = FALSE])
      #Sigma_Y_candidate_inv = inversa_sim_def_positiva(diag_nj + tau * XtX_candidate)
      aux_Sigma_candidate = crossprod(X[idx, index_gamma_candidate_1 , drop  = FALSE]) + tau_inv_matrix_tot[1:p_activos_candidate, 1:p_activos_candidate]
      Sigma_Y_candidate_inv =  (ones_matrix[1:nj[h], 1:nj[h]]  -  X[idx, index_gamma_candidate_1 , drop  = FALSE] %*%inversa_sim_def_positiva(aux_Sigma_candidate) %*% Xt[index_gamma_candidate_1, idx, drop  = FALSE])
      
      log_k = 0.5*(-log_determinante(aux_Sigma_candidate) + log_determinante(aux_Sigma_t_c[,,h][index_gamma_1,index_gamma_1, drop = F]) - aux_A)
      
      A = A + dt_mvstudent_rpp_inv_compare_centered_opt(x_old = Y[idx] - X[idx, index_gamma_1 , drop  = FALSE]%*% mu_0[index_gamma_1], 
                                                        x_candidate = Y[idx] - X[idx, index_gamma_candidate_1 , drop  = FALSE]%*% mu_0[index_gamma_candidate_1],
                                                        log_k = log_k,
                                                        Sigma_inv_candidate = Sigma_Y_candidate_inv, 
                                                        Sigma_inv_old = Sigma_Y_old_inv, nu= 2*a_sigma2,
                                                        a_0_b_0 = (a_sigma2/b_sigma2),logd = T)
      
      #A = A + dt_mvstudent_rpp_inv_compare_centered_chol(Y[idx] - X[idx, index_gamma_1 , drop  = FALSE]%*% mu_0[index_gamma_1], 
       #                                             Sigma_inv_candidate = Sigma_Y_candidate_inv, 
        #                                            Sigma_inv_old = Sigma_Y_old_inv, nu= 2*a_sigma2,
         #                                           a_0_b_0 = (a_sigma2/b_sigma2),logd = T)
      
    }
    
    A = min(1,exp(A))
    if (runif(1) <= A){
      Matrix$gamma[, iter + 1] = candidate_gamma
      change_gamma = T
    }else{
      Matrix$gamma[, iter + 1] = Matrix$gamma[, iter ]
      change_gamma = F
    }
    
    
    
    if (stream){
      # Actualizar la barra de progreso
      if (iter %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
        pb$tick()
      }
    }
  }
  
  # Burnint y Thining
  secuencia = seq(burnit + 1, n_chain, by = thining)
  Matrix$alpha =  Matrix$alpha[secuencia]
  Matrix$H =  Matrix$H[secuencia]
  Matrix$ind_cluster =  Matrix$ind_cluster[,secuencia]
  Matrix$gamma =  Matrix$gamma[,secuencia]
  
  if(n_muestra){
    
    Matrix$y_pred <- DPM_Gaussian_Uni_Neal_muestra(Y, Matrix, n_muestra)
    return(Matrix)
  }
  
  return(Matrix)
}


DDPM_Gaussian_Uni_Neal_muestra <- function(X_new, Matrix, n_muestra, stream = F){
  
  Y = Matrix$Y
  X = Matrix$X
  p = dim(X)[2]
  n = dim(X)[1]
  mu_0 = Matrix$param_init$mu_0 # 0
  tau = Matrix$param_init$tau # 10
  a_sigma2 = Matrix$param_init$a_sigma2 #2 #2
  b_sigma2 = Matrix$param_init$b_sigma2 #.5 # E(sigma^2)= b/(a-1)
  a_alpha = Matrix$param_init$a_alpha #.1
  b_alpha = Matrix$param_init$b_alpha #.1
  tau_inv = diag(1/tau, p)
  mu_0_tau = mu_0/tau
  nu_0 = 2*a_sigma2 
  #sigma_0 = sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new)))
  
  y_pred = array(0 ,dim=c(n_muestra))
  l_chain = dim(Matrix$alpha)
  index = rep(1:l_chain, ceiling(n_muestra/l_chain))[1:n_muestra]
  
  i = 1
  if (stream){
    pb <- progress_bar$new(
      format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
      total = length(index)/100,
      clear = FALSE,
      width = 60
    )
  }
  # Calcular la predictiva
  for (iter in index){
    
    index_gamma_1 = Matrix$gamma[, iter] == 1
    
    nj = tabulate(factor(Matrix$ind_cluster[, iter], levels = 1:Matrix$H[iter]) , 
                  nbins = Matrix$H[iter])
    
    ind_pred = sample_rcpp_stl(1:(Matrix$H[iter]+1), 
                               prob = c(sapply(1:Matrix$H[iter], 
                                               function(x)  nj[x]/(n + Matrix$alpha[iter])),
                                        Matrix$alpha[iter]/(n + Matrix$alpha[iter])))
    
    
    if(ind_pred <= Matrix$H[iter]){ 
      
      idx <- Matrix$ind_cluster[, iter] == ind_pred
      
      a_t = a_sigma2 + nj[ind_pred]/2
      aux_Sigma = crossprod(X[idx,index_gamma_1, drop = FALSE]) + tau_inv[index_gamma_1, index_gamma_1]
      Sigma_t = inversa_sim_def_positiva(aux_Sigma)
      mu_t = Sigma_t%*%(crossprod(X[idx,index_gamma_1, drop = FALSE], Y[idx, drop = FALSE]) + mu_0_tau[index_gamma_1])
      b_t = b_sigma2 + 0.5 * (crossprod(Y[idx]) + crossprod(mu_0[index_gamma_1])/tau - t(mu_t) %*% aux_Sigma %*% mu_t)
      
      
      y_pred[i] =  rt_ls(1, mu = X_new[index_gamma_1] %*%mu_t, 
                         sigma = sqrt(b_t/a_t * (1 +  X_new[index_gamma_1]%*% Sigma_t %*%X_new[index_gamma_1])) , 
                         nu = 2*a_t) 
    } else y_pred[i] = rt_ls(1, mu = X_new[index_gamma_1]%*%mu_0[index_gamma_1], sigma = sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new[index_gamma_1]))), nu = nu_0)
    i = i + 1
    
    if (stream){
      # Actualizar la barra de progreso
      if (i %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
        pb$tick()
      }
    }
  }
  
  
  return(y_pred * Matrix$sd_Y_or + Matrix$mean_Y_or)
}



DDPM_Gaussian_Uni_Neal_dens = function(X_new, Matrix, grilla, n_iter = 1e3, stream = F){
  
  
  Y = Matrix$Y
  X = Matrix$X
  p = dim(X)[2]
  n = dim(X)[1]
  mu_0 = Matrix$param_init$mu_0 # 0
  tau = Matrix$param_init$tau # 10
  a_sigma2 = Matrix$param_init$a_sigma2 #2 #2
  b_sigma2 = Matrix$param_init$b_sigma2 #.5 # E(sigma^2)= b/(a-1)
  a_alpha = Matrix$param_init$a_alpha #.1
  b_alpha = Matrix$param_init$b_alpha #.1
  tau_inv = diag(1/tau, p)
  mu_0_tau = mu_0/tau
  nu_0 = 2*a_sigma2 
  #sigma_0 = sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new)))
  
  l_secuencia = dim(Matrix$alpha)
  n_utilizado = n_iter
  if (l_secuencia < n_utilizado){
    warning("n_iter debe ser menor a la cantidad de muestras en Matrix")
    return()
  }
  dens = array(0 ,dim=c(length(grilla), n_utilizado))
  
  if (stream){
    pb <- progress_bar$new(
      format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
      total = n_utilizado/100,
      clear = FALSE,
      width = 60
    )
  }
  
  for (iter in (l_secuencia - n_utilizado + 1):l_secuencia){
    
    index_gamma_1 = Matrix$gamma[, iter] == 1
    
    nj = tabulate(factor(Matrix$ind_cluster[, iter], levels = 1:Matrix$H[iter]) , 
                  nbins = Matrix$H[iter])
    dens[, iter - (l_secuencia - n_utilizado)] = Matrix$alpha[iter]/(n + Matrix$alpha[iter]) * dt_ls(x = (grilla - Matrix$mean_Y_or)/Matrix$sd_Y_or , 
                                                                                                     mu = c(X_new[index_gamma_1]%*%mu_0[index_gamma_1]), sigma = as.vector(sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new[index_gamma_1])))), 
                                                                                                     nu = nu_0) / Matrix$sd_Y_or
    for (c in 1:Matrix$H[iter]){
      
      n_c = nj[c]
      idx <- Matrix$ind_cluster[, iter] == c
      
      a_t = a_sigma2 + n_c/2
      aux_Sigma = crossprod(X[idx,index_gamma_1, drop = FALSE]) + tau_inv[index_gamma_1, index_gamma_1]
      Sigma_t = inversa_sim_def_positiva(aux_Sigma)
      mu_t = Sigma_t%*%(crossprod(X[idx,index_gamma_1, drop = FALSE], Y[idx, drop = FALSE]) + mu_0_tau[index_gamma_1])
      b_t = b_sigma2 + 0.5 * (crossprod(Y[idx]) + crossprod(mu_0[index_gamma_1])/tau - t(mu_t) %*% aux_Sigma %*% mu_t)
      
      dens[, iter - (l_secuencia - n_utilizado)] =  dens[, iter - (l_secuencia - n_utilizado)] + n_c/(n +Matrix$alpha[iter]) * dt_ls(x = (grilla - Matrix$mean_Y_or)/Matrix$sd_Y_or, mu = c(X_new[index_gamma_1] %*%mu_t), 
                                                                                                                                     sigma = as.vector(sqrt(b_t/a_t * (1 +  X_new[index_gamma_1]%*% Sigma_t %*%X_new[index_gamma_1]))) , 
                                                                                                                                     nu = 2*a_t) / Matrix$sd_Y_or
    }
    
    if (stream){
      if (iter %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
        pb$tick()
      }
    }
  }
  return(dens)  
}




DDPM_Gaussian_Uni_Neal_dist = function(X_new, Matrix, grilla, n_iter = 1e3, stream = F){
  
  
  Y = Matrix$Y
  X = Matrix$X
  p = dim(X)[2]
  n = dim(X)[1]
  mu_0 = Matrix$param_init$mu_0 # 0
  tau = Matrix$param_init$tau # 10
  a_sigma2 = Matrix$param_init$a_sigma2 #2 #2
  b_sigma2 = Matrix$param_init$b_sigma2 #.5 # E(sigma^2)= b/(a-1)
  a_alpha = Matrix$param_init$a_alpha #.1
  b_alpha = Matrix$param_init$b_alpha #.1
  tau_inv = diag(1/tau, p)
  mu_0_tau = mu_0/tau
  nu_0 = 2*a_sigma2 
  #sigma_0 = sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new)))
  
  l_secuencia = dim(Matrix$alpha)
  n_utilizado = n_iter
  if (l_secuencia < n_utilizado){
    warning("n_iter debe ser menor a la cantidad de muestras en Matrix")
    return()
  }
  dist = array(0 ,dim=c(length(grilla), n_utilizado))
  
  if (stream){
    pb <- progress_bar$new(
      format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
      total = n_utilizado/100,
      clear = FALSE,
      width = 60
    )
  }
  
  for (iter in (l_secuencia - n_utilizado + 1):l_secuencia){
    
    index_gamma_1 = Matrix$gamma[, iter] == 1
    
    nj = tabulate(factor(Matrix$ind_cluster[, iter], levels = 1:Matrix$H[iter]) , 
                  nbins = Matrix$H[iter])
    dist[, iter - (l_secuencia - n_utilizado)] = Matrix$alpha[iter]/(n + Matrix$alpha[iter]) * pt_ls(x = (grilla - Matrix$mean_Y_or)/Matrix$sd_Y_or , 
                                                                                                     mu = c(X_new[index_gamma_1]%*%mu_0[index_gamma_1]), sigma = as.vector(sqrt((b_sigma2/a_sigma2)*(1 + tau*crossprod(X_new[index_gamma_1])))), 
                                                                                                     nu = nu_0)
    for (c in 1:Matrix$H[iter]){
      
      n_c = nj[c]
      idx <- Matrix$ind_cluster[, iter] == c
      
      a_t = a_sigma2 + n_c/2
      aux_Sigma = crossprod(X[idx,index_gamma_1, drop = FALSE]) + tau_inv[index_gamma_1, index_gamma_1]
      Sigma_t = inversa_sim_def_positiva(aux_Sigma)
      mu_t = Sigma_t%*%(crossprod(X[idx,index_gamma_1, drop = FALSE], Y[idx, drop = FALSE]) + mu_0_tau[index_gamma_1])
      b_t = b_sigma2 + 0.5 * (crossprod(Y[idx]) + crossprod(mu_0[index_gamma_1])/tau - t(mu_t) %*% aux_Sigma %*% mu_t)
      
      dist[, iter - (l_secuencia - n_utilizado)] =  dist[, iter - (l_secuencia - n_utilizado)] + n_c/(n +Matrix$alpha[iter]) * pt_ls(x = (grilla - Matrix$mean_Y_or)/Matrix$sd_Y_or, mu = c(X_new[index_gamma_1] %*%mu_t), 
                                                                                                                                     sigma = as.vector(sqrt(b_t/a_t * (1 +  X_new[index_gamma_1]%*% Sigma_t %*%X_new[index_gamma_1]))) , 
                                                                                                                                     nu = 2*a_t)
    }
    if (stream){
      if (iter %% 100 == 0) {  # Actualiza solo cada 100 iteraciones
        pb$tick()
      }
    }
  }
  return(dist)  
}
