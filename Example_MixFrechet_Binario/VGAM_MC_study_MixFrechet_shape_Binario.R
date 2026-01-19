rm(list = ls())
library(VGAM)


## Generalized pareto distribution

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


for (n in c(200,400,800,1600)){
set.seed(1000)
# Covariable continua
s = sample(c(0,1),n, replace = T)
X <- model.matrix(~ as.factor(s))
colnames(X) = c("Intercept", "x1")

log_alpha1 <- X %*% beta1
log_alpha2 <- X %*% beta2
alpha_true1 <- sqrt(exp(log_alpha1))
alpha_true2 <- sqrt(exp(log_alpha2))

print(c(min(alpha_true1), max(alpha_true1)))
print(c(min(alpha_true2), max(alpha_true2)))


n_mc = 100
result_scale = array(dim = c(2, n_mc))
result_shape = array(dim = c(2, n_mc))

for (iter_mc in 1:n_mc){
  
  
  y <- as.vector(rmixfrechet(n=n, seed = iter_mc))
  
  #y[y<.01] = .01
  dat = data.frame(
    y = y
  )
  dat = data.frame(cbind(y,X)) 

  success <- FALSE
  while (!success) {
    fit_try <- tryCatch(
      {
  modelo_completo <- vglm(y ~ x1 + 1,
                          frechet(location = 0, lscale = "loglink", lshape = "loglink"),
                          data = dat,
                          trace = F, coefstart = rep(0,4))
      },
  error = function(e) NULL
    )
    
    if (!is.null(fit_try)) {
      modelo_completo <- fit_try
      success <- TRUE
    }
  }
  success <- FALSE
  while (!success) {
    fit_try <- tryCatch(
      {
        modelo_incompleto <-vglm(y ~1,
                                 frechet(location = 0, lscale = "loglink",
                                         lshape = "loglink"),
                                 trace = F, data = dat,  coefstart = rep(0,2))
      },
      error = function(e) NULL
    )
    
    if (!is.null(fit_try)) {
      modelo_incompleto <- fit_try
      success <- TRUE
    }
  }
  
  
  if(BIC(modelo_completo) < BIC(modelo_incompleto)){
    
    result_scale[, iter_mc] = modelo_completo@coefficients[c(1,3)]
    result_shape[, iter_mc] = modelo_completo@coefficients[c(2,4)]
  }else{
    result_scale[, iter_mc] = c(modelo_incompleto@coefficients[c(1)],0)
    result_shape[, iter_mc] = c(modelo_incompleto@coefficients[c(2)],0)
  }
  
}

path_save = "../Example_MixFrechet_Binario/VGAM/"
save(result_scale, file = paste0(path_save,"result_scale_n",n, ".Rdata"))
save(result_shape, file = paste0(path_save, "result_shape_n",n, ".Rdata"))

}


#S = summary(modelo_completo);S
#exp(S@coefficients[c(1,3)]%*%c(1,1)) # Sigma
#exp(S@coefficients[c(2,4)]%*%c(1,0))
#exp(S@coefficients[c(2,4)]%*%c(1,1))

## Density estimation ----
grid_vals_or_body = seq(1e-6, 
                        qfrechet(.8, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))),
                        length.out = 5000)


for (n in c(200,400,800,1600)){

  
load(paste0(path_save,"result_scale_n",n, ".Rdata"))
load(paste0(path_save,"result_shape_n",n, ".Rdata"))
result_dens1 = array(dim = c(length(grid_vals_or_body), n_mc))
result_dens2 = array(dim = c(length(grid_vals_or_body), n_mc))

for (iter_mc in 1:n_mc){
  
  result_dens1[, iter_mc] = dfrechet(grid_vals_or_body, location = 0, 
                                     scale = exp(result_scale[, iter_mc] %*% c(0,1)),
                                     shape = exp(result_shape[, iter_mc] %*% c(0,1)))
  
  result_dens2[, iter_mc] = dfrechet(grid_vals_or_body, location = 0, 
                                     scale = exp(result_scale[, iter_mc] %*% c(1,1)),
                                     shape = exp(result_shape[, iter_mc] %*% c(1,1)))
}

result_dens_list = list()
mean_dens1 = rowMeans(result_dens1)
lower_dens1 <- apply(result_dens1, 1, quantile, probs = 0.025)
upper_dens1 <- apply(result_dens1, 1, quantile, probs = 0.975)
result_dens_list[[1]] = data.frame(x = grid_vals_or_body,
                                   mean_dens = mean_dens1,
                                   lower_dens = lower_dens1,
                                   upper_dens = upper_dens1)
mean_dens2 = rowMeans(result_dens2)
lower_dens2 <- apply(result_dens2, 1, quantile, probs = 0.025)
upper_dens2 <- apply(result_dens2, 1, quantile, probs = 0.975)
result_dens_list[[2]] = data.frame(x = grid_vals_or_body,
                                   mean_dens = mean_dens2,
                                   lower_dens = lower_dens2,
                                   upper_dens = upper_dens2)

save(result_dens_list, file = paste0("../Example_MixFrechet_Binario/VGAM/dens_n",n, ".Rdata"))
}

## Distribution estimation ----

grid_vals_or_tail = seq(qfrechet(.8, location = mu_true1, scale = sigma_true1, shape = max(c(alpha_true1, alpha_true2))), 
                        qfrechet(.99, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))),
                        length.out = 5000)


for (n in c(200,400,800,1600)){
  

  load(paste0(path_save,"result_scale_n",n, ".Rdata"))
  load(paste0(path_save,"result_shape_n",n, ".Rdata"))
  result_dist1 = array(dim = c(length(grid_vals_or_tail), n_mc))
  result_dist2 = array(dim = c(length(grid_vals_or_tail), n_mc))
  
  for (iter_mc in 1:n_mc){
    
    result_dist1[, iter_mc] = pfrechet(grid_vals_or_tail, location = 0, 
                                       scale = exp(result_scale[, iter_mc] %*% c(0,1)),
                                       shape = exp(result_shape[, iter_mc] %*% c(0,1)))
    
    result_dist2[, iter_mc] = pfrechet(grid_vals_or_tail, location = 0, 
                                       scale = exp(result_scale[, iter_mc] %*% c(1,1)),
                                       shape = exp(result_shape[, iter_mc] %*% c(1,1)))
  }
  
  result_dist_list = list()
  mean_dist1 = rowMeans(result_dist1)
  lower_dist1 <- apply(result_dist1, 1, quantile, probs = 0.025)
  upper_dist1 <- apply(result_dist1, 1, quantile, probs = 0.975)
  result_dist_list[[1]] = data.frame(x = grid_vals_or_tail,
                                     mean_dist = mean_dist1,
                                     lower_dist = lower_dist1,
                                     upper_dist = upper_dist1)
  mean_dist2 = rowMeans(result_dist2)
  lower_dist2 <- apply(result_dist2, 1, quantile, probs = 0.025)
  upper_dist2 <- apply(result_dist2, 1, quantile, probs = 0.975)
  result_dist_list[[2]] = data.frame(x = grid_vals_or_tail,
                                     mean_dist = mean_dist2,
                                     lower_dist = lower_dist2,
                                     upper_dist = upper_dist2)
  
  save(result_dist_list, file = paste0("../Example_MixFrechet_Binario/VGAM/dist_n",n, ".Rdata"))
}



# Plot ----
library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales)

## Densidad body ------
# Función para crear un gráfico
dens_teorica1 <- data.frame(
  x = grid_vals_or_body,
  y = dmixfrechet(grid_vals_or_body, X_new = c(1,0))
)
dens_teorica2 <- data.frame(
  x = grid_vals_or_body,
  y = dmixfrechet(grid_vals_or_body, X_new = c(1,1))
)

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
    xlim(0,5) +
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


#output_dir = paste0("../Example_MixFrechet_Binario/")




n = 200
load(paste0("../Example_MixFrechet_Binario/VGAM/dens_n",n, ".Rdata"))
p1 <- make_plot(result_dens_list[[1]], result_dens_list[[2]], "n = 200") +
  theme(axis.title.y = element_text(size = 20),
        legend.position = "none")
print(p1)


n = 400
load(paste0("../Example_MixFrechet_Binario/VGAM/dens_n",n, ".Rdata"))
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
load(paste0("../Example_MixFrechet_Binario/VGAM/dens_n",n, ".Rdata"))
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
load(paste0("../Example_MixFrechet_Binario/VGAM/dens_n",n, ".Rdata"))
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

path_save = "../Example_MixFrechet_Binario/VGAM/"
ggsave(filename = paste(path_save,"density_or_space_n200_n400_n800_n1600_VGAM", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)



## Survival plot ------
# Función para crear un gráfico

dist_or_teorica1 <- data.frame(
  x = grid_vals_or_tail,
  y = Vectorize(function(x){pmixfrechet(x, X_new = c(1,0))})(grid_vals_or_tail)
)
dist_or_teorica2 <- data.frame(
  x = grid_vals_or_tail,
  y = Vectorize(function(x){pmixfrechet(x, X_new = c(1,1))})(grid_vals_or_tail)
)
make_plot_surv <- function(dist_df1, dist_df2, title, xlims, ylims) {
  
  dist_df1 <- dist_df1 %>%
    mutate(
      mean_surv  = 1 - mean_dist,
      lower_surv = 1 - upper_dist,
      upper_surv = 1 - lower_dist,
      grupo = "control"  # Grupo rojo
    )
  
  dist_df2 <- dist_df2 %>%
    mutate(
      mean_surv  = 1 - mean_dist,
      lower_surv = 1 - upper_dist,
      upper_surv = 1 - lower_dist,
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

x_lims_global = c(qfrechet(.9, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))),
                  qfrechet(.99, location = mu_true1, scale = sigma_true1, shape = min(c(alpha_true1, alpha_true2))))

y_lims_global = c(-15,-1)
n=200
load(paste0("../Example_MixFrechet_Binario/VGAM/dist_n",n, ".Rdata"))
p1 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 200", xlims = x_lims_global, ylims =y_lims_global ) +
  labs(y = "Log(Survival)")+
  theme(axis.title.y = element_text(size = 20),
        legend.position = "none")
print(p1)

n=400
load(paste0("../Example_MixFrechet_Binario/VGAM/dist_n",n, ".Rdata"))
p2 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 400", xlims = x_lims_global, ylims = y_lims_global) +
  labs(y = NULL) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
  ) 
print(p2)

n=800
load(paste0("../Example_MixFrechet_Binario/VGAM/dist_n",n, ".Rdata"))
p3 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 800", xlims = x_lims_global, ylims = y_lims_global) +
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
load(paste0("../Example_MixFrechet_Binario/VGAM/dist_n",n, ".Rdata"))
p4 <- make_plot_surv(result_dist_list[[1]], result_dist_list[[2]], "n = 1600", xlims = x_lims_global, ylims = y_lims_global) +
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

ggsave(filename = paste(path_save,"Log_Survival_n200_n400_n800_n1600_VGAM", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)


