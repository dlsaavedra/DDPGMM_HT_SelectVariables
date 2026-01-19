rm(list = ls())
library(VGAM)

# LOAD FUNCTIONS ----------------------------------------------------------
source("../DDPGMM_Uni_Neal_Selection.R")


# Betas y Matriz de diseño ----

#####  Generalized pareto distribution ##
# Tamaño de muestra
n <- 200  #200, 400, 800
for (n in c(200,400,800)){
set.seed(1000)
# Parámetros de la muestra ----
p = 6  #5 predictores - 3 relevantes, 100 - 20 relevantes.
rademacher = sample(c(-1, 1), size = p, replace = TRUE)
betas_or <- rademacher #matrix(rademacher*runif(p*n_cluster,.5,2), ncol = n_cluster) 
X = cbind(1,matrix(rnorm(n*(p-1), mean = 0, sd = 1), nrow = n))
list_p_activos = as.matrix(cbind(1,expand.grid(rep(list(0:1), p-1))), ncol = p)
list_betas = list_p_activos *  rep(betas_or, each = nrow(list_p_activos))
alpha_true <- 0.9     # alpha = 1/xi
mu_true <- 0
for (iter_beta in 1:nrow(list_betas)){
  betas = t(list_betas[iter_beta,, drop = F])
  log_sigma <- X %*% betas
  sigma_true <- sqrt(exp(log_sigma))
  print(c(min(sigma_true), max(sigma_true)))
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
output_dir = paste0("MC_study_p5/MC_study_n",n,"/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# Parallel ----
library(parallel)

# --- Configuración de clusters ---
n_cores <- 10  # Usa todos los núcleos menos uno
n_mc = 100
#1:nrow(list_betas)
for (iter_beta in 1:nrow(list_betas)){
  
  
  betas = t(list_betas[iter_beta,, drop = F])
  log_sigma <- X %*% betas
  sigma_true <- sqrt(exp(log_sigma))
  print(betas)
  #cl <- makeCluster(detectCores() - 1)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    source("../DDPGMM_Uni_Neal_Selection.R")  # Se ejecuta 1 vez por worker
    library(VGAM)
  })
  
  clusterExport(cl, c("n", "mu_true", "sigma_true", "alpha_true", "X", "param_chain", 
                      "param_init", "list_p_activos", "iter_beta", "output_dir"))
  
  # --- Ejecución paralela ---
  resultados <- parLapply(cl, 1:n_mc, function(iter_mc) {
    set.seed(iter_mc)
    y <- as.vector(rfrechet(n = n, location = mu_true, scale = sigma_true, shape = alpha_true))
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


burnit = 2000
thining = 3
secuencia = seq(burnit + 1, n_chain, by = thining)
n_seq = length(secuencia)
beta_paste = sapply(1:nrow(list_betas),function(x) paste(list_p_activos[x,], collapse = ","))

# n = 200 ----
output_dir = "MC_study_p5/MC_study_n200/"
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

matrix_result_n200 = matrix_result
save(matrix_result_n200, file = paste0(output_dir,"matrix_result_n200_prob_mc", ".Rdata"))
matrix_result_n200_mean  = apply(matrix_result, c(1, 2), mean)
save(matrix_result_n200_mean, file = paste0(output_dir,"matrix_result_n200_prob_mean", ".Rdata"))

# n = 400 ----
output_dir = "MC_study_p5/MC_study_n400/"
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

matrix_result_n400 = matrix_result
save(matrix_result_n400, file = paste0(output_dir,"matrix_result_n400_prob_mc", ".Rdata"))
matrix_result_n400_mean  = apply(matrix_result, c(1, 2), mean)
save(matrix_result_n400_mean, file = paste0(output_dir,"matrix_result_n400_prob_mean", ".Rdata"))


# n = 800 ----
output_dir = "MC_study_p5/MC_study_n800/"
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
#matrix_result_max <- matrix(0, nrow = dim(matrix_result)[1], ncol = dim(matrix_result)[2])

#for (k in 1:dim(matrix_result_n800)[3]) {
#  pos_max <- apply(matrix_result[,,k], 1, which.max)
#  B <- matrix(0, dim(matrix_result)[1], dim(matrix_result)[2])
#  B[cbind(1:dim(matrix_result)[1], pos_max)] <- 1
#  matrix_result_max <- matrix_result_max + B
#}
#save(matrix_result_max, file = paste0(output_dir,"matrix_result_n800_max", ".Rdata"))


## Reordenar matrices
output_dir =  "MC_study_p5/MC_study_n200/"
load(file = paste0(output_dir,"matrix_result_n200_prob_mean", ".Rdata"))
output_dir =  "MC_study_p5/MC_study_n400/"
load(file = paste0(output_dir,"matrix_result_n400_prob_mean", ".Rdata"))
output_dir =  "MC_study_p5/MC_study_n800/"
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

output_dir =  "MC_study_p5/MC_study_n200/"
save(matrix_result_n200_mean, file = paste0(output_dir,"matrix_result_n200_prob_mean_2", ".Rdata"))

output_dir =  "MC_study_p5/MC_study_n400/"
save(matrix_result_n400_mean, file = paste0(output_dir,"matrix_result_n400_prob_mean_2", ".Rdata"))

output_dir =  "MC_study_p5/MC_study_n800/"
save(matrix_result_n800_mean, file = paste0(output_dir,"matrix_result_n800_prob_mean_2", ".Rdata"))



# Pixel Plot ----

output_dir =  "MC_study_p5/MC_study_n800/"
load(file = paste0(output_dir,"matrix_result_n800_prob_mean_2", ".Rdata"))
output_dir =  "MC_study_p5/MC_study_n400/"
load(file = paste0(output_dir,"matrix_result_n400_prob_mean_2", ".Rdata"))
output_dir =  "MC_study_p5/MC_study_n200/"
load(file = paste0(output_dir,"matrix_result_n200_prob_mean_2", ".Rdata"))

#matrix_result_n200_mean
#matrix_result_n400_mean
#matrix_result_n800_mean


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
  labs(title = "n = 200", x = "Model Estimate", y = "Frechet Distribution \n True Model") +
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
final_plot
path_save = "MC_study_p5/Img_Study/"
ggsave(filename = paste(path_save,"Pixel_plot", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)






