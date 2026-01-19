rm(list = ls())

source("../DDPGMM_Uni_Neal_Selection.R")

library(VGAM)

# Mixture Frechet ----
beta1 <- c(1, -3) 
beta2 <- c(-1.5, 3)
sigma_true1 <- 1    
sigma_true2 <- 1  
mu_true1 <- 0
mu_true2 <- 0
w_1 = 0.6
w_2 = 1 - w_1


# Funciones de la mixtura Frechet ----
indicator <- function(x, m) {
  ifelse(x > m, 1, 0)
}


rmixfrechet = function(n, seed = seed){
  
  set.seed(seed)
  # Paso 1: Generar un vector que indica de qué componente viene cada observación
  component <- rbinom(n, 1, w_1)  # 1 para el primer componente, 0 para el segundo
  # Paso 2: Generar observaciones según el componente seleccionado
  
  alpha_true = alpha_true1*(component == 1) + alpha_true2*(component == 0)
  mu_true = mu_true1*(component == 1) + mu_true2*(component == 0)
  sigma_true = sigma_true1*(component == 1) + sigma_true2*(component == 0)
  
  data <- rfrechet(n = n, location = mu_true, scale = sigma_true, shape = alpha_true)
  
  return(data)
}

dmixfrechet = function(x, X_new){
  
  alpha1 = sqrt(exp( X_new %*% beta1))
  alpha2 = sqrt(exp( X_new %*% beta2))
  return(
    w_1 * dfrechet(x, location = mu_true1, scale = sigma_true1, shape = alpha1) +
      w_2 * dfrechet(x, location = mu_true2, scale = sigma_true2, shape = alpha2)
  )
}

pmixfrechet = function(x, X_new){
  
  alpha1 = sqrt(exp( X_new %*% beta1))
  alpha2 = sqrt(exp( X_new %*% beta2))
  return(
    w_1 * pfrechet(x, location = mu_true1, scale = sigma_true1, shape = alpha1) +
      w_2 * pfrechet(x, location = mu_true2, scale = sigma_true2, shape = alpha2)
  )
}


#V1 = Vectorize(function(x){dmixfrechet(x, X_new = c(1,0))})

#hist(rmixfrechet(n=n, seed = 4)[X[,2]==0], freq = F, breaks = c(seq(0,20,by = .1), Inf), right = F, xlim = c(0,10))
#curve(V1(x), col = "red", 
#      xlim = c(0, 10), n = 1e3,
#      ylab = "y", xlab = "t", add = T)

#hist(log(rmixfrechet(n=n, seed = 4)[X[,2]==0]), freq = F, breaks = 50,  xlim = c(-10, 10))
#curve(V1(exp(x))*exp(x), col = "red", n = 1e4,
#      ylab = "y", xlab = "t", add = T)


#V2 = Vectorize(function(x){dmixfrechet(x, X_new = c(1,1))})

#curve(V2(x), col = "red", 
#      xlim = c(0, 5), n = 1e5,
#      ylab = "y", xlab = "t")

#hist(log(rmixfrechet(n=n, seed = 2)[X[,2]==1]), freq = F, breaks = 50,  xlim = c(-10, 10))
#curve(V2(exp(x))*exp(x), col = "red",
#      n = 1e3,
#      ylab = "y", xlab = "t", add = T)

for (n in c(200)){
#for (n in c(200,400,800, 1600)){
  set.seed(1000)
  # Covariable continua
  s = sample(c(0,1),n, replace = T)
  X <- model.matrix(~ as.factor(s))
  # Parámetros dependientes de covariables
  log_alpha1 <- X %*% beta1
  log_alpha2 <- X %*% beta2
  alpha_true1 <- sqrt(exp(log_alpha1))
  alpha_true2 <- sqrt(exp(log_alpha2))
  
  print(c(min(alpha_true1), max(alpha_true1)))
  print(c(min(alpha_true2), max(alpha_true2)))
  
# Hiperparametros ----
p = dim(X)[2]
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
output_dir = paste0("../Example_MixFrechet_Binario/MC_study_n",n,"/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# Parallel ----
library(parallel)

# --- Configuración de clusters ---
n_cores <- 10  # Usa todos los núcleos menos uno
n_mc = 100

#cl <- makeCluster(detectCores() - 1)
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  source("../DDPGMM_Uni_Neal_Selection.R")  # Se ejecuta 1 vez por worker
  library(VGAM)
  
})

clusterExport(cl, c("n", "mu_true1", "mu_true2", "sigma_true1", "sigma_true2", "alpha_true1", "alpha_true2",
                    "w_1","w_2", "rmixfrechet", "X", "param_chain", 
                    "param_init", "output_dir"))

# --- Ejecución paralela ---
resultados <- parLapply(cl, 1:n_mc, function(iter_mc) {
 
  y <- as.vector(rmixfrechet(n=n, seed = iter_mc))
  Matrix <- DDPM_Gaussian_Uni_Neal_Selection(log(y), X, param_chain, param_init, Y_scale = TRUE, stream = FALSE)
  save(Matrix, file = paste0(output_dir, "result_", iter_mc, ".Rdata"))
})
# --- Cerrar clusters ---
stopCluster(cl)
}

# Creación de Matrix ----
library(progress)

for (n in c(200)){
#for (n in c(200,400,800, 1600)){
output_dir = paste0("../Example_MixFrechet_Binario/MC_study_n", n,"/")

n_chain = 10000
burnit = 4000
thining = 3
secuencia = seq(burnit + 1, n_chain, by = thining)
n_seq = length(secuencia)
matrix_result = array(0, dim  = c(2, n_mc))
beta_paste = c("1,0", "1,1")

pb <- progress_bar$new(
  format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
  total = n_mc,
  clear = FALSE,
  width = 60
)
for (iter_mc in 1:n_mc){
  
  file_path = paste0(output_dir, "result_", iter_mc, ".Rdata")
  load(file = file_path)
  
  vec_str <- apply(Matrix$gamma[,secuencia], 2, paste, collapse = ",")
  matrix_result[ , iter_mc] = table(factor(vec_str, levels = beta_paste))/n_seq
  
  pb$tick()
}


matrix_result_mean  = apply(matrix_result, c(1), mean);print(matrix_result_mean)
matrix_result_max <-  table(factor(apply(matrix_result, 2, which.max), levels = 1:2))/n_mc;print(matrix_result_max)
}

# Densidad espacio logaritmo -----

#

for (n in c(200)){ 
#for (n in c(200,400,800, 1600)){
output_dir = paste0("../Example_MixFrechet_Binario/MC_study_n", n,"/")

grid_vals = seq(-10, 15,length.out = 1000)
list_X_pred = list(c(1,0),c(1,1))

burnit = 2000
thining = 3
n_mc = 100
n_utilizado = 1000


result_dens_list= list()
for (iter_pred in 1:2){
  
  pb <- progress_bar$new(
    format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
    total = n_mc,
    clear = FALSE,
    width = 60
  )
  result_dens = matrix(0, nrow = n_mc, ncol = length(grid_vals))
  for (iter_mc in 1:n_mc){
    
    file_path = paste0(output_dir, "result_", iter_mc, ".Rdata")
    load(file = file_path)
    secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
    Matrix$alpha = Matrix$alpha[secuencia]
    Matrix$H = Matrix$H[secuencia]
    Matrix$ind_cluster = Matrix$ind_cluster[,secuencia]
    Matrix$gamma = Matrix$gamma[,secuencia]
    
    dens = DDPM_Gaussian_Uni_Neal_dens(X_new = list_X_pred[[iter_pred]],
                                       Matrix, grilla = grid_vals,
                                       n_iter = n_utilizado, stream = F)
    result_dens[iter_mc, ] = rowMeans(dens)
    
    pb$tick()
  }
  
  mean_dens = colMeans(result_dens)
  lower_dens <- apply(result_dens, 2, quantile, probs = 0.025)
  upper_dens <- apply(result_dens, 2, quantile, probs = 0.975)
  result_dens_list[[iter_pred]] = data.frame(x = grid_vals,
                                             mean_dens = mean_dens,
                                             lower_dens = lower_dens,
                                             upper_dens = upper_dens)
  
}


save(result_dens_list, file = paste0(output_dir,"result_dens_n", n, ".Rdata"))
}

#Distribución Espacio Original cola -----


grid_vals_or_tail = seq(qfrechet(.8, location = mu_true1, scale = sigma_true1, shape = max(c(sqrt(exp(c(beta1[1], sum(beta1)))), 
                                                                                             sqrt(exp(c(beta2[1], sum(beta2))))))), 
                   qfrechet(.99, location = mu_true1, scale = sigma_true1, shape = min(c(sqrt(exp(c(beta1[1], sum(beta1)))), 
                                                                                         sqrt(exp(c(beta2[1], sum(beta2))))))),
                   length.out = 5000)

for (n in c(200)){
  output_dir = paste0("../Example_MixFrechet_Binario/MC_study_n", n,"/")
result_dist_list= list()
for (iter_pred in 1:2){
  
  pb <- progress_bar$new(
    format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
    total = n_mc,
    clear = FALSE,
    width = 60
  )
  result_dist = matrix(0, nrow = n_mc, ncol = length(grid_vals_or_tail))
  for (iter_mc in 1:n_mc){
    
    file_path = paste0(output_dir, "result_", iter_mc, ".Rdata")
    load(file = file_path)
    secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
    Matrix$alpha = Matrix$alpha[secuencia]
    Matrix$H = Matrix$H[secuencia]
    Matrix$ind_cluster = Matrix$ind_cluster[,secuencia]
    Matrix$gamma = Matrix$gamma[,secuencia]
    
    dist = DDPM_Gaussian_Uni_Neal_dist(X_new = list_X_pred[[iter_pred]],
                                         Matrix, grilla = log(grid_vals_or_tail),
                                         n_iter = n_utilizado, stream = F)
    result_dist[iter_mc, ] = rowMeans(dist)
    
    pb$tick()
  }
  
  mean_dist = colMeans(result_dist)
  lower_dist <- apply(result_dist, 2, quantile, probs = 0.025)
  upper_dist <- apply(result_dist, 2, quantile, probs = 0.975)
  result_dist_list[[iter_pred]] = data.frame(x = grid_vals_or_tail,
                                             mean_dens = mean_dist,
                                             lower_dens = lower_dist,
                                             upper_dens = upper_dist)
  
}


save(result_dist_list, file = paste0(output_dir,"result_dist_or_tail_n", n, ".Rdata"))
}

#Densidad Espacio Original Cuerpo -----


grid_vals_or_body = seq(1e-6, 
                        qfrechet(.8, location = mu_true1, scale = sigma_true1, shape = min(c(sqrt(exp(c(beta1[1], sum(beta1)))), 
                                                                                             sqrt(exp(c(beta2[1], sum(beta2))))))),
                        length.out = 5000)

for (n in c(200)){
#for (n in c(200,400,800, 1600)){
output_dir = paste0("../Example_MixFrechet_Binario/MC_study_n", n,"/")

result_dens_list= list()
for (iter_pred in 1:2){
  
  pb <- progress_bar$new(
    format = " Procesando [:bar] :percent en :elapsedfull con ETA :eta",
    total = n_mc,
    clear = FALSE,
    width = 60
  )
  result_dens = matrix(0, nrow = n_mc, ncol = length(grid_vals_or_body))
  for (iter_mc in 1:n_mc){
    
    file_path = paste0(output_dir, "result_", iter_mc, ".Rdata")
    load(file = file_path)
    secuencia = seq(burnit + 1, length(Matrix$H), by = thining)
    Matrix$alpha = Matrix$alpha[secuencia]
    Matrix$H = Matrix$H[secuencia]
    Matrix$ind_cluster = Matrix$ind_cluster[,secuencia]
    Matrix$gamma = Matrix$gamma[,secuencia]
    
    dens = DDPM_Gaussian_Uni_Neal_dens(X_new = list_X_pred[[iter_pred]],
                                       Matrix, grilla = log(grid_vals_or_body),
                                       n_iter = n_utilizado, stream = F)/grid_vals_or_body
    result_dens[iter_mc, ] = rowMeans(dens)
    
    pb$tick()
  }
  
  mean_dens = colMeans(result_dens)
  lower_dens <- apply(result_dens, 2, quantile, probs = 0.025)
  upper_dens <- apply(result_dens, 2, quantile, probs = 0.975)
  result_dens_list[[iter_pred]] = data.frame(x = grid_vals_or_body,
                                             mean_dens = mean_dens,
                                             lower_dens = lower_dens,
                                             upper_dens = upper_dens)
  
}

save(result_dens_list, file = paste0(output_dir,"result_dens_or_body_n", n, ".Rdata"))
}



# Plot Densidad espacio logaritmo ----
library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales)

path_save = "../Example_MixFrechet_Binario/Img_MC_study/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)


dens_teorica1 <- data.frame(
  x = grid_vals,
  y = dmixfrechet(exp(grid_vals), X_new = c(1,0))*exp(grid_vals)
)
dens_teorica2 <- data.frame(
  x = grid_vals,
  y = dmixfrechet(exp(grid_vals), X_new = c(1,1))*exp(grid_vals)
)

# Función para crear un gráfico
make_plot <- function(dist_df1, dist_df2, title) {
  
  dist_df1$grupo <- "control"
  
  dist_df2$grupo <- "treatment"
  
  df_total <- bind_rows(dist_df1, dist_df2)
  
  p <- ggplot(df_total, aes(x = x, fill = grupo, color = grupo)) +
    # Bandas de confianza
    geom_ribbon(aes(ymin = lower_dens, ymax = upper_dens), alpha = 0.2, color = NA) +
    # Línea segmentada
    geom_line(aes(y = mean_dens), linetype = "dotted", size = 1) +
    # Línea teórica 1 (roja)
    geom_line(data = dens_teorica1, aes(x = x, y = y), inherit.aes = FALSE, color = "red", size = .5) +
    # Línea teórica 2 (azul)
    geom_line(data = dens_teorica2, aes(x = x, y = y), inherit.aes = FALSE, color = "blue", size = .5) +
    xlim(-5,10) +
    ylim(0,0.6) +
    scale_fill_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    scale_color_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    labs(title = title, y = NULL, x = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}


output_dir = paste0("../Example_MixFrechet_Binario/")
n=200
load(paste0(output_dir,"MC_study_n",n,"/result_dens_n", n, ".Rdata"))
p1 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 200") +
  theme(axis.title.y = element_text(size = 20),
        legend.position = "none")
print(p1)

n=400
load(paste0(output_dir,"MC_study_n",n,"/result_dens_n", n, ".Rdata"))
p2 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 400") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
print(p2)


n=800
load(paste0(output_dir,"MC_study_n",n,"/result_dens_n", n, ".Rdata"))
p3 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 800") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
print(p3)

n=1600
load(paste0(output_dir,"MC_study_n",n,"/result_dens_n", n, ".Rdata"))
p4 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 1600") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )
print(p4)

final_plot <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 4) +
  theme(plot.title = element_text(hjust = 0.5))
print(final_plot)

ggsave(filename = paste(path_save,"density_log_space_n200_n400_n800_n1600", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)


# Plot Survival espacio original -----

dist_or_teorica1 <- data.frame(
  x = grid_vals_or_tail,
  y = Vectorize(function(x){pmixfrechet(x, X_new = c(1,0))})(grid_vals_or_tail)
)
dist_or_teorica2 <- data.frame(
  x = grid_vals_or_tail,
  y = Vectorize(function(x){pmixfrechet(x, X_new = c(1,1))})(grid_vals_or_tail)
)

x_lims_global = c(qfrechet(.9, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))),
     qfrechet(.99, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))))
# Función para crear un gráfico
make_plot_surv <- function(dist_df1, dist_df2, title, xlims, ylims) {
  
  dist_df1 <- dist_df1 %>%
    mutate(
      mean_surv  = 1 - mean_dens,
      lower_surv = 1 - upper_dens,
      upper_surv = 1 - lower_dens,
      grupo = "control"  # Grupo rojo
    )
  
  dist_df2 <- dist_df2 %>%
    mutate(
      mean_surv  = 1 - mean_dens,
      lower_surv = 1 - upper_dens,
      upper_surv = 1 - lower_dens,
      grupo = "treatment"  # Grupo azul
    )
  
  df_total <- bind_rows(dist_df1, dist_df2)
  
  p <- ggplot(df_total, aes(x = x, fill = grupo, color = grupo)) +
    # Bandas de confianza
    geom_ribbon(aes(ymin = log(lower_surv), ymax = log(upper_surv)), alpha = 0.2, color = NA) +
    # Línea segmentada
    geom_line(aes(y = log(mean_surv)), linetype = "dotted", size = 1) +
    # Línea teórica 1 (roja)
    geom_line(data = dist_or_teorica1, aes(x = x, y = log(1-y)), inherit.aes = FALSE, color = "red", size = .5) +
    # Línea teórica 2 (azul)
    geom_line(data = dist_or_teorica2, aes(x = x, y = log(1-y)), inherit.aes = FALSE, color = "blue", size = .5) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    scale_fill_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    scale_color_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    labs(title = title, y = NULL, x = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}

n=200
load(paste0(output_dir,"MC_study_n", n,"/result_dist_or_tail_n", n, ".Rdata"))
p1 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 200", xlims = x_lims_global, ylims = c(-15,-1)) +
  labs(y = "Log(Survival)")+
  theme(axis.title.y = element_text(size = 20),
        legend.position = "none")
print(p1)

n=400
load(paste0(output_dir,"MC_study_n", n,"/result_dist_or_tail_n", n, ".Rdata"))
p2 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 400", xlims = x_lims_global, ylims = c(-15,-1)) +
  labs(y = NULL) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
  ) 
print(p2)

n=800
load(paste0(output_dir,"MC_study_n", n,"/result_dist_or_tail_n", n, ".Rdata"))
p3 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 800", xlims = x_lims_global, ylims = c(-15,-1)) +
  labs(y = NULL) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
print(p3)

n=1600
load(paste0(output_dir,"MC_study_n", n,"/result_dist_or_tail_n", n, ".Rdata"))
p4 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 1600", xlims = x_lims_global, ylims = c(-15,-1)) +
  labs(y = NULL) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )
print(p4)


final_plot <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 4)
print(final_plot)

ggsave(filename = paste(path_save,"Log_Survival_n200_n400_n800_n1600", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)


# Plot dens espacio original ----


dens_teorica1 <- data.frame(
  x = grid_vals_or_body,
  y = dmixfrechet(grid_vals_or_body, X_new = c(1,0))
)
dens_teorica2 <- data.frame(
  x = grid_vals_or_body,
  y = dmixfrechet(grid_vals_or_body, X_new = c(1,1))
)

# Función para crear un gráfico
make_plot <- function(dist_df1, dist_df2, title) {
  
  dist_df1$grupo <- "control"
  
  dist_df2$grupo <- "treatment"
  
  df_total <- bind_rows(dist_df1, dist_df2)
  
  p <- ggplot(df_total, aes(x = x, fill = grupo, color = grupo)) +
    # Bandas de confianza
    geom_ribbon(aes(ymin = lower_dens, ymax = upper_dens), alpha = 0.2, color = NA) +
    # Línea segmentada
    geom_line(aes(y = mean_dens), linetype = "dotted", size = 1) +
    # Línea teórica 1 (roja)
    geom_line(data = dens_teorica1, aes(x = x, y = y), inherit.aes = FALSE, color = "red", size = .5) +
    # Línea teórica 2 (azul)
    geom_line(data = dens_teorica2, aes(x = x, y = y), inherit.aes = FALSE, color = "blue", size = .5) +
    xlim(0,4) +
    ylim(0,1) +
    scale_fill_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    scale_color_manual(name = "Cases", values = c("control" = "red", "treatment" = "blue")) +
    labs(title = title, y = NULL, x = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}


output_dir = paste0("../Example_MixFrechet_Binario/")

n = 200
load(paste0(output_dir,"MC_study_n", n, "/result_dens_or_body_n", n, ".Rdata"))
p1 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 200") +
  theme(axis.title.y = element_text(size = 20),
        legend.position = "none")
print(p1)


n = 400
load(paste0(output_dir,"MC_study_n", n, "/result_dens_or_body_n", n, ".Rdata"))
p2 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 400") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
print(p2)


n = 800
load(paste0(output_dir,"MC_study_n", n, "/result_dens_or_body_n", n, ".Rdata"))
p3 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 800") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
print(p3)


n = 1600
load(paste0(output_dir,"MC_study_n", n, "/result_dens_or_body_n", n, ".Rdata"))
p4 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 1600") +
  labs(y = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )
print(p4)


final_plot <- (p1 + p2 + p3 + p4) +
  plot_layout(ncol = 4) +
  theme(plot.title = element_text(hjust = 0.5))
print(final_plot)

ggsave(filename = paste(path_save,"density_or_space_n200_n400_n800_n1600", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)

