rm(list = ls())
library(VGAM)

# LOAD FUNCTIONS ----------------------------------------------------------
source("../DDPGMM_Uni_Neal_Selection.R")



#####  Generalized pareto distribution ##
# Tamaño de muestra
n <- 800  #200, 400, 800
for (n in c(200,400,800)){
  set.seed(1000)
  # Parámetros de la muestra ----
  p = 6  #5 predictores - 3 relevantes, 100 - 20 relevantes.
  rademacher = sample(c(-1, 1), size = p, replace = TRUE)
  betas_or <- rademacher #matrix(rademacher*runif(p*n_cluster,.5,2), ncol = n_cluster) 
  X = cbind(1,matrix(rnorm(n*(p-1), mean = 0, sd = 1), nrow = n))
  list_p_activos = as.matrix(cbind(1,expand.grid(rep(list(0:1), p-1))), ncol = p)
  list_betas = list_p_activos *  rep(betas_or, each = nrow(list_p_activos))
  sigma_true <- 1     # alpha = 1/xi
  mu_true <- 0
  for (iter_beta in 1:nrow(list_betas)){
    betas = t(list_betas[iter_beta,, drop = F])
    log_xi<- X %*% betas
    xi_true <- sqrt(exp(log_xi))
    print(c(min(1/xi_true), max(1/xi_true)))
  }


# Hiperparametros ----
mu_0 = rep(0, p)
tau = 10
a_sigma2 = 2 #2
b_sigma2 = .01 #.5 # E(sigma^2)= b/(a-1)
a_alpha = .1
b_alpha = .1
H_max = 4
n_cluster_ini = 2
zeta = 1
gamma_inicio = c(1, rep(1,p-1))

param_init = list(mu_0 = mu_0, tau = tau,
                  a_sigma2 = a_sigma2, b_sigma2 = b_sigma2, 
                  a_alpha = a_alpha, b_alpha = b_alpha, H_max = H_max,
                  n_cluster_ini = n_cluster_ini,
                  zeta = zeta, gamma_inicio = gamma_inicio)

# Parametros de la cadena ----
burnit = 0
n_chain = 10000
thining = 1
param_chain = list(burnit= burnit, n_chain = n_chain, thining = thining )
output_dir = paste0("../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n",n,"/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# Parallel ----
library(parallel)

# --- Configuración de clusters ---
n_cores <- 10  # Usa todos los núcleos menos uno
n_mc = 100
#1:nrow(list_betas)
for (iter_beta in 1:nrow(list_betas)){ 
  
  betas = t(list_betas[iter_beta,, drop = F])
  print(betas)
  log_xi<- X %*% betas
  xi_true <- sqrt(exp(log_xi))
  
  #cl <- makeCluster(detectCores() - 1)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    source("../DDPGMM_Uni_Neal_Selection.R")  # Se ejecuta 1 vez por worker
    library(VGAM)
  })
  
  clusterExport(cl, c("n", "mu_true", "sigma_true", "xi_true", "X", "param_chain", 
                      "param_init", "list_p_activos", "iter_beta", "output_dir"))
  
  # --- Ejecución paralela ---
  resultados <- parLapply(cl, 1:n_mc, function(iter_mc) {
    set.seed(iter_mc)
    
    y = rgpd(n = n, location = mu_true,  scale = sigma_true, shape = xi_true)
    Matrix <- DDPM_Gaussian_Uni_Neal_Selection(log(y), X, param_chain, param_init, Y_scale = TRUE, stream = FALSE)
    save(Matrix, file = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata"))
  })
  # --- Cerrar clusters ---
  stopCluster(cl)
}
}
# Si paramos el codigo ejecutar: 
#stopCluster(cl)

# Creación de Matrix ----

output_dir = "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n800/"
burnit = 2000
thining = 3
secuencia = seq(burnit + 1, n_chain, by = thining)
n_seq = length(secuencia)
beta_paste = sapply(1:nrow(list_betas),function(x) paste(list_p_activos[x,], collapse = ","))
matrix_result = array(0, dim  = c(nrow(list_betas), nrow(list_betas), n_mc))

for (iter_beta in 1:nrow(list_betas)){
  print(iter_beta)
  result = c()
  for (iter_mc in 1:n_mc){
    
    file_path = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata")
    load(file = file_path)
    
    vec_str <- apply(Matrix$gamma[,secuencia], 2, paste, collapse = ",")
    matrix_result[iter_beta, , iter_mc] = table(factor(vec_str, levels = beta_paste))/n_seq
  }
}


matrix_result_n800 = matrix_result
save(matrix_result_n800, file = paste0(output_dir,"matrix_result_n800_prob_mc", ".Rdata"))
matrix_result_n800_mean  = apply(matrix_result, c(1, 2), mean)
save(matrix_result_n800_mean, file = paste0(output_dir,"matrix_result_n800_prob_mean", ".Rdata"))

## Guardar la matrix Máxima ----
matrix_result_n800_max <- matrix(0, nrow = nrow(list_betas), ncol = nrow(list_betas))

for (k in 1:dim(matrix_result)[3]) {
  pos_max <- apply(matrix_result[,,k], 1, which.max)
  B <- matrix(0, nrow(list_betas), nrow(list_betas))
  B[cbind(1:nrow(list_betas), pos_max)] <- 1
  matrix_result_n800_max <- matrix_result_n800_max + B
}
save(matrix_result_n800_max, file = paste0(output_dir,"matrix_result_n800_max", ".Rdata"))


## Reordenar matrices ----

output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n200/"
load(file = paste0(output_dir,"matrix_result_n200_prob_mean", ".Rdata"))
output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n400/"
load(file = paste0(output_dir,"matrix_result_n400_prob_mean", ".Rdata"))
output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n800/"
load(file = paste0(output_dir,"matrix_result_n800_prob_mean", ".Rdata"))

matrix_result_n200_mean
matrix_result_n400_mean
matrix_result_n800_mean
orden = order(rowSums(list_p_activos), decreasing = T)
matrix_result_n200_mean = matrix_result_n200_mean[orden, rev(orden)]
matrix_result_n400_mean = matrix_result_n400_mean[orden, rev(orden)]
matrix_result_n800_mean = matrix_result_n800_mean[orden, rev(orden)]
M_ordenada <- list_p_activos[orden, ]
nombres <- apply(M_ordenada, 1, function(x) paste(x, collapse = ","))
# Asignar nombres a filas y columnas
rownames(matrix_result_n200_mean) <- nombres
colnames(matrix_result_n200_mean) <- rev(nombres)
# Asignar nombres a filas y columnas
rownames(matrix_result_n400_mean) <- nombres
colnames(matrix_result_n400_mean) <- rev(nombres)
# Asignar nombres a filas y columnas
rownames(matrix_result_n800_mean) <- nombres
colnames(matrix_result_n800_mean) <- rev(nombres)

output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n200/"
save(matrix_result_n200_mean, file = paste0(output_dir,"matrix_result_n200_prob_mean_2", ".Rdata"))

output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n400/"
save(matrix_result_n400_mean, file = paste0(output_dir,"matrix_result_n400_prob_mean_2", ".Rdata"))

output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n800/"
save(matrix_result_n800_mean, file = paste0(output_dir,"matrix_result_n800_prob_mean_2", ".Rdata"))



# Pixel Plot ----

output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n800/"
load(file = paste0(output_dir,"matrix_result_n800_prob_mean_2", ".Rdata"))
output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n400/"
load(file = paste0(output_dir,"matrix_result_n400_prob_mean_2", ".Rdata"))
output_dir =  "../Example_Pareto_Generalized_shape/MC_study_p5/MC_study_n200/"
load(file = paste0(output_dir,"matrix_result_n200_prob_mean_2", ".Rdata"))

matrix_result_n200_mean
matrix_result_n400_mean
matrix_result_n800_mean


library(ggplot2)
library(reshape2)
library(patchwork)  # Para disposición horizontal elegante

# Suponiendo que matrix_result_n200 y matrix_result_n400 son matrices ya cargadas
# Convertir ambas matrices a formato largo
df1 <- melt(matrix_result_n200_mean)
df2 <- melt(matrix_result_n400_mean)
df3 <- melt(matrix_result_n800_mean)

# Nombres para los ejes, como factores para mantener el orden
df1$Var1 <- factor(df1$Var1, levels = rev(rownames(matrix_result_n200_mean)))
df1$Var2 <- factor(df1$Var2, levels = colnames(matrix_result_n200_mean))

df2$Var1 <- factor(df2$Var1, levels = rev(rownames(matrix_result_n400_mean)))
df2$Var2 <- factor(df2$Var2, levels = colnames(matrix_result_n400_mean))

df3$Var1 <- factor(df3$Var1, levels = rev(rownames(matrix_result_n800_mean)))
df3$Var2 <- factor(df3$Var2, levels = colnames(matrix_result_n800_mean))


# Crear ambos plots
p1 <- ggplot(df1, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1)) +
  labs(title = "n = 200", x = "Model Estimate", y = "GP Distribution \n True Model") +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 20)
  )


p2 <- ggplot(df2, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1)) +
  labs(title = "n = 400", x = "Model Estimate", y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20)
  )

p3 <- ggplot(df3, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1)) +
  labs(title = "n = 800", x = "Model Estimate", y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20)
  )


# Combinar horizontalmente
final_plot = p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "right")
path_save = "../Example_Pareto_Generalized_shape//MC_study_p5/Img_Study/"
ggsave(filename = paste(path_save,"Pixel_plot", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)


#Combinacion GPD Pixel Plot -----




# 
# # Densidad predictiva estimada --
# l_betas_pred = 1:16#c(1,7,16)
# X_pred = c(1,-1,1,-1,1)
# 
# 
# list_grid = list()
# #list_grid[[1]] = seq(-5, 10, length.out = 1000)
# #list_grid[[7]] = seq(-10, 25, length.out = 1000)
# #list_grid[[16]] = seq(-25, 55,length.out = 1000)
# for (i in l_betas_pred){
#   list_grid[[i]] = seq(-25, 55,length.out = 1000)
# }
# 
# burnit = 2000
# thining = 3
# n_mc = 100
# n_utilizado = 1000
# 
# ## n200 --
# output_dir = "../Example_Pareto_Generalized_shape/MC_study_n200/"
# 
# result_dens_betas_n200 = list()
# for (iter_beta in l_betas_pred){
#   
#   result_dens = matrix(0, nrow = n_mc, ncol = length(list_grid[[iter_beta]]))
#   for (iter_mc in 1:n_mc){
#     
#     file_path = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata")
#     load(file = file_path)
#     secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
#     
#     dens = DDPM_Gaussian_Uni_Neal_dens(X_new = X_pred,
#                                        Matrix, grilla = list_grid[[iter_beta]],
#                                        n_iter = n_utilizado)
#     result_dens[iter_mc, ] = rowMeans(dens)
#   }
#   
#   mean_dens = colMeans(result_dens)
#   lower_dens <- apply(result_dens, 2, quantile, probs = 0.025)
#   upper_dens <- apply(result_dens, 2, quantile, probs = 0.975)
#   result_dens_betas_n200[[iter_beta]] = data.frame(mean_dens = mean_dens, 
#                                                    lower_dens = lower_dens,
#                                                    upper_dens = upper_dens)
#   
# }
# 
# save(result_dens_betas_n200, file = paste0(output_dir,"result_dens_n200", ".Rdata"))
# 
# ## n400 --
# 
# output_dir = "../Example_Pareto_Generalized_shape/MC_study_n400/"
# 
# result_dens_betas_n400 = list()
# for (iter_beta in l_betas_pred){
#   
#   result_dens = matrix(0, nrow = n_mc, ncol = length(list_grid[[iter_beta]]))
#   for (iter_mc in 1:n_mc){
#     
#     file_path = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata")
#     load(file = file_path)
#     secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
#     
#     dens = DDPM_Gaussian_Uni_Neal_dens(X_new = X_pred,
#                                        Matrix, grilla = list_grid[[iter_beta]],
#                                        n_iter = n_utilizado)
#     result_dens[iter_mc, ] = rowMeans(dens)
#   }
#   
#   mean_dens = colMeans(result_dens)
#   lower_dens <- apply(result_dens, 2, quantile, probs = 0.025)
#   upper_dens <- apply(result_dens, 2, quantile, probs = 0.975)
#   result_dens_betas_n400[[iter_beta]] = data.frame(mean_dens = mean_dens, 
#                                                    lower_dens = lower_dens,
#                                                    upper_dens = upper_dens)
#   
# }
# 
# 
# save(result_dens_betas_n400, file = paste0(output_dir,"result_dens_n400", ".Rdata"))
# 
# ## n800 --
# 
# output_dir = "../Example_Pareto_Generalized_shape/MC_study_n800/"
# 
# result_dens_betas_n800 = list()
# for (iter_beta in l_betas_pred){
#   
#   result_dens = matrix(0, nrow = n_mc, ncol = length(list_grid[[iter_beta]]))
#   for (iter_mc in 1:n_mc){
#     
#     file_path = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata")
#     load(file = file_path)
#     secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
#     
#     dens = DDPM_Gaussian_Uni_Neal_dens(X_new = X_pred,
#                                        Matrix, grilla = list_grid[[iter_beta]],
#                                        n_iter = n_utilizado)
#     result_dens[iter_mc, ] = rowMeans(dens)
#   }
#   
#   mean_dens = colMeans(result_dens)
#   lower_dens <- apply(result_dens, 2, quantile, probs = 0.025)
#   upper_dens <- apply(result_dens, 2, quantile, probs = 0.975)
#   result_dens_betas_n800[[iter_beta]] = data.frame(mean_dens = mean_dens, 
#                                                    lower_dens = lower_dens,
#                                                    upper_dens = upper_dens)
#   
# }
# 
# save(result_dens_betas_n800, file = paste0(output_dir,"result_dens_n800", ".Rdata"))
# 
# 
# 
# # Grafico densidad predictiva ---
# 
# output_dir =  "../Example_Pareto_Generalized_shape/MC_study_n800/"
# load(file = paste0(output_dir,"result_dens_n800", ".Rdata"))
# output_dir =  "../Example_Pareto_Generalized_shape/MC_study_n400/"
# load(file = paste0(output_dir,"result_dens_n400", ".Rdata"))
# output_dir =  "../Example_Pareto_Generalized_shape/MC_study_n200/"
# load(file = paste0(output_dir,"result_dens_n200", ".Rdata"))
# 
# 
# 
# for (iter_beta in 1:dim(list_betas)[1]){
#   path_save = "../Example_Pareto_Generalized_shape/Img_MC_study/Density/"
#   beta_pred = list_betas[iter_beta,]
#   alpha_pred <- sqrt(exp(X_pred %*% beta_pred))
#   label_y <- paste(list_p_activos[iter_beta,], collapse = "")
#   grid_vals <- list_grid[[iter_beta]]
#   
#   # Densidades teóricas y estimadas para n = 200
#   dens_200 <- result_dens_betas_n200[[iter_beta]]
#   dens_df_200 <- data.frame(
#     x = grid_vals,
#     mean = dens_200$mean_dens,
#     lower = dens_200$lower_dens,
#     upper = dens_200$upper_dens,
#     type = "n = 200"
#   )
#   
#   # Densidades teóricas y estimadas para n = 400
#   dens_400 <- result_dens_betas_n400[[iter_beta]]
#   dens_df_400 <- data.frame(
#     x = grid_vals,
#     mean = dens_400$mean_dens,
#     lower = dens_400$lower_dens,
#     upper = dens_400$upper_dens,
#     type = "n = 400"
#   )
#   
#   # Densidades teóricas y estimadas para n = 800
#   dens_800 <- result_dens_betas_n800[[iter_beta]]
#   dens_df_800 <- data.frame(
#     x = grid_vals,
#     mean = dens_800$mean_dens,
#     lower = dens_800$lower_dens,
#     upper = dens_800$upper_dens,
#     type = "n = 800"
#   )
#   
#   # Densidad teórica (misma en ambos)
#   dens_teorica <- data.frame(
#     x = grid_vals,
#     y = exp(grid_vals) * dfrechet(exp(grid_vals), location = mu_true, scale = sigma_true, shape = as.numeric(alpha_pred))
#   )
#   
#   # Función para crear un gráfico
#   make_plot <- function(dens_df, title, y_null = F) {
#     ggplot(dens_df, aes(x = x)) +
#       geom_line(data = dens_teorica, aes(y = y), color = "red", size = 1) +
#       geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
#       geom_line(aes(y = mean), color = "blue", linetype = "dashed", size = 1.2) +
#       coord_cartesian(ylim = c(0, 0.6)) +
#       labs(title = title, y = label_y, x = NULL) +
#       theme_minimal(base_size = 14) +
#       theme(
#         plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_text(size = 10)
#       )
#   }
#   
#   
#   
#   # Crear los dos gráficos
#   p1 <- make_plot(dens_df_200, "n = 200")
#   p2 <- make_plot(dens_df_400, "n = 400") +
#     labs(y = NULL) +
#     theme(axis.title.x = element_blank())
#   p3 <- make_plot(dens_df_800, "n = 800") +
#     labs(y = NULL) +
#     theme(axis.title.x = element_blank())
#   # Combinarlos horizontalmente
#   final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
#   
#   # Mostrar
#   print(final_plot)
#   
#   # Guardar en ruta personalizada
#   ggsave(filename = paste(path_save,"densities_n200_n400_n800_",iter_beta, ".png", sep= ""), plot = final_plot,
#          width = 10, height = 4.5, dpi = 300)
# }