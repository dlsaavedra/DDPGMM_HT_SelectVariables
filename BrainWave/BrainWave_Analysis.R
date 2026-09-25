rm(list = ls())

# LOAD FUNCTIONS ----------------------------------------------------------
source("../DDPGMM_Uni_Neal_Selection.R")


load(file = "brainwave.RData")


# Exploratorio

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # Para combinar gráficos

# Nombres de condiciones
cond_names <- c("Mathematics", "Relaxation", "Music", "Color", "Video", "Relax and think")

# Transformar de wide a long para identificar condiciones
long_df <- brainwave %>%
  pivot_longer(cols = starts_with("x"), names_to = "condition", values_to = "active") %>%
  filter(active == 1) %>%
  mutate(
    condition = factor(condition,
                       levels = paste0("x", 1:6),
                       labels = cond_names)
  )

# Crear data frame largo con bandas alpha y beta
plot_df <- long_df %>%
  pivot_longer(cols = c(Y1, Y2), names_to = "Band", values_to = "Power") %>%
  mutate(Band = factor(Band, levels = c("Y1", "Y2"), labels = c("Alpha", "Beta")))

# Crear un gráfico por cada condición y banda
make_plot <- function(cond, band) {
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  ggplot(df_sub, aes(x = time, y = Power)) +
    geom_line(color = "black") +
    ylim(0, 200)+
    labs(title = ifelse(band == "Alpha", cond, ""),
         x = "Time (s)",
         y = expression(paste(Power, " (", mu, "V"^2, ")"))) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

# Crear los 12 gráficos (6 condiciones × 2 bandas)
all_plots <- lapply(cond_names, function(cond) {
  list(
    make_plot(cond, "Alpha"),
    make_plot(cond, "Beta")
  )
}) %>% unlist(recursive = FALSE)

# Organizar en una grilla: 3 filas × 4 columnas
layout <- wrap_plots(all_plots, ncol = 4, axes = "collect_y")

# Mostrar
print(layout)

# Guardar la figura
ggsave("brainwave_power_plot.png", layout, width = 12, height = 9, dpi = 300)

#Histogramas -----

# Estimar la máxima densidad (frecuencia relativa) para fijar un ylim común
# usando la misma binwidth que se usará en el gráfico

# Función para graficar histogramas con frecuencia relativa
make_histogram <- function(cond, band) {
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  ggplot(df_sub, aes(x = Power)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 2,
                   fill = "gray30", color = "white") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.08)) +
    labs(
      title = ifelse(band == "Alpha", cond, ""),
      x = expression(paste(Power, " (", mu, "V"^2, ")")),
      y = "Relative Frequency"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
  
  
}

# Generar los 12 histogramas (6 condiciones × 2 bandas)
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Alpha")
    #make_histogram(cond, "Beta")
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar
ggsave("brainwave_histograms.png",
       layout_hist, width = 10, height = 4.5, dpi = 300)

# Violin Plot ------


# Crear violin plots agrupados

# Asegurar niveles
plot_df <- plot_df %>% 
  mutate(
    condition = factor(condition, levels = c("Mathematics", "Relaxation", "Music",
                                             "Color", "Video", "Relax and think")),
    Band      = factor(Band, levels = c("Alpha", "Beta"))
  )

band_colors <- c(Alpha = "#4da6ff", Beta = "#ff6666")

# --- Función para un solo panel ---
make_violin_dot <- function(band, cond) {
  df <- plot_df %>% filter(Band == band, condition == cond)
  ggplot(df, aes(x = "", y = Power)) +
    geom_violin(aes(fill = Band), alpha = .70, color = NA) +
    scale_fill_manual(values = band_colors) +
    geom_dotplot(binaxis = "y",
                 stackdir = "centerwhole",
                 binwidth = 2,
                 dotsize = .2,
                 alpha = 0.2,
                 fill = band_colors[band]) +
    coord_cartesian(ylim = c(0, 200)) +
    labs(title = ifelse(band == "Alpha", cond, ""),  # título solo en la fila Alpha
         x = NULL, y = expression(paste(Power, " (", mu, "V"^2, ")"))) +
    theme_minimal() +               # ← sin theme_void() para poder dibujar el eje
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = .5),
          plot.margin = margin(4, 4, 4, 4))
}

# --- Crear lista de gráficos ---
plots <- lapply(levels(plot_df$condition), function(cond) {
  list(
    make_violin_dot("Alpha", cond),
    make_violin_dot("Beta",  cond)
  )
}) |> unlist(recursive = FALSE)

# --- Ensamblar en cuadrícula 6×2 ---
grid <- wrap_plots(plots, ncol = 4, byrow = TRUE, axes = "collect_y")

# --- Etiquetas globales ---
grid_final <- grid +
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    theme = theme(plot.margin = margin(10, 10, 10, 10),
                  plot.title = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(guides = "collect")  # compartir leyendas si las hubiera
  

# Mostrar
grid_final
# Guardar
ggsave("brainwave_violinplot.png", grid_final, width = 15, height = 10, dpi = 500)



#Histogramas Espacio Logaritmo -----

# Estimar la máxima densidad (frecuencia relativa) para fijar un ylim común
# usando la misma binwidth que se usará en el gráfico

# Función para graficar histogramas con frecuencia relativa
make_histogram <- function(cond, band) {
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  ggplot(df_sub, aes(x = log(Power))) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = .2,
                   fill = "gray30", color = "white") +
    coord_cartesian(xlim = c(-3, log(200)), ylim = c(0, 0.6)) +
    labs(
      title = ifelse(band == "Alpha", cond, ""),
      x = expression(paste(Power, " (", mu, "V"^2, ")")),
      y = "Relative Frequency"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

# Generar los 12 histogramas (6 condiciones × 2 bandas)
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Alpha"),
    make_histogram(cond, "Beta")
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 4, axes = "collect_y")

# Mostrar
print(layout_hist)





#Ajustar el modelo -----
## Hiperparametros ----
p = 6
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
n_chain = 50000
thining = 1
param_chain = list(burnit= burnit, n_chain = n_chain, thining = thining )


X = as.matrix(cbind(1,brainwave[,4:8]))


Matrix_Y1 <- DDPM_Gaussian_Uni_Neal_Selection(log(brainwave$Y1), X, param_chain, param_init, Y_scale = TRUE, stream = T)
save(Matrix_Y1, file = "Matrix_Y1.Rdata")

Matrix_Y2 <- DDPM_Gaussian_Uni_Neal_Selection(log(brainwave$Y2), X, param_chain, param_init, Y_scale = TRUE, stream = T)
save(Matrix_Y2, file = "Matrix_Y2.Rdata")

# Analisis ----
load("Matrix_Y1.Rdata")
load("Matrix_Y2.Rdata")

n_chain = length(Matrix_Y1$H)
secuencia = seq(10001, n_chain, by = 4) 
table(apply(Matrix_Y1$gamma[,secuencia], 2, paste, collapse = ","))
table(apply(Matrix_Y2$gamma[,secuencia], 2, paste, collapse = ","))


Matrix_Y1_burn = Matrix_Y1
Matrix_Y1_burn$alpha =  Matrix_Y1_burn$alpha[secuencia]
Matrix_Y1_burn$H =  Matrix_Y1_burn$H[secuencia]
Matrix_Y1_burn$ind_cluster =  Matrix_Y1_burn$ind_cluster[,secuencia]
Matrix_Y1_burn$gamma =  Matrix_Y1_burn$gamma[,secuencia]

table(apply(Matrix_Y1_burn$gamma, 2, paste, collapse = ","))
ts.plot(colSums(Matrix_Y1_burn$gamma))
ts.plot(Matrix_Y1_burn$gamma[4,])

apply(Matrix_Y1_burn$gamma, 1, mean)

Matrix_Y2_burn = Matrix_Y2
Matrix_Y2_burn$alpha =  Matrix_Y2_burn$alpha[secuencia]
Matrix_Y2_burn$H =  Matrix_Y2_burn$H[secuencia]
Matrix_Y2_burn$ind_cluster =  Matrix_Y2_burn$ind_cluster[,secuencia]
Matrix_Y2_burn$gamma =  Matrix_Y2_burn$gamma[,secuencia]

table(apply(Matrix_Y2_burn$gamma, 2, paste, collapse = ","))
ts.plot(colSums(Matrix_Y2_burn$gamma))
ts.plot(Matrix_Y2_burn$gamma[4,])

apply(Matrix_Y2_burn$gamma, 1, mean)

#Densidad Espacio Original-----
n_utilizado = 1000
grid_vals = seq(1e-2,200,length.out = 1000)
band_names <- c("Alpha", "Beta")

x_new_conditions <- data.frame(
  condition = cond_names,
  X_new = I(list(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,0,1,0,0,0),
    c(1,0,0,1,0,0),
    c(1,0,0,0,1,0),
    c(1,0,0,0,0,1)
  ))
)

# 3. y 4. Iterar y crear el data.frame final
density_df <- list()

for (cond in cond_names) {
  for (band in band_names) {
    
    # Selecciona la matriz y X_new correspondientes
    matrix_to_use <- if (band == "Alpha") {
      Matrix_Y1_burn
    } else {
      Matrix_Y2_burn
    }
    
    x_new_to_use <- x_new_conditions %>%
      filter(condition == cond) %>%
      pull(X_new) %>%
      `[[`(1)
    
    # Llama a la función
    result <- DDPM_Gaussian_Uni_Neal_dens(
      X_new = x_new_to_use,
      Matrix = matrix_to_use,
      grilla = log(grid_vals),
      n_iter = n_utilizado,
      stream = TRUE
    )
    # Crea un data.frame temporal con los resultados
    temp_df <- data.frame(
      x = grid_vals,
      mean_dens = rowMeans(result)/grid_vals,
      lower_dens = apply(result, 1, function(x) quantile(x, 0.025))/grid_vals,
      upper_dens = apply(result, 1, function(x) quantile(x, 0.975))/grid_vals,
      condition = cond,
      Band = band
    )
    
    density_df[[paste(cond, band, sep = "_")]] <- temp_df
  }
}

# Combina los resultados
density_df <- bind_rows(density_df)


# Función para graficar histogramas con la densidad y sus intervalos
make_histogram <- function(cond, band, xlims, ylims ) {
  
  # Filtra los datos del histograma
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  
  # Filtra los datos de la densidad y los intervalos
  density_sub <- density_df %>% filter(condition == cond, Band == band)
  
  ggplot(df_sub, aes(x = Power)) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = 2,
                   fill = "gray30", color = "white") +
    
    # Añade la banda de credibilidad (Intervalo de Confianza)
    geom_ribbon(data = density_sub, aes(x = x, ymin = lower_dens, ymax = upper_dens), 
                fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
    
    # Añade la línea de la densidad media
    geom_line(data = density_sub, aes(x = x, y = mean_dens), 
              color = "blue", size = 0.5, inherit.aes = FALSE) +
    
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title = cond,
      x = expression(paste(Power, " (", mu, "V"^2, ")")),
      y = "Relative Frequency"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

## Banda Alpha ----
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Alpha", xlims = c(0,100), ylims = c(0, 0.08))
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"density_estimation_alpha", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(histograms[i], ncol = 1)
  ggsave(filename = paste(path_save,"density_estimation_alpha_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}

## Banda Beta ----
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Beta", xlims = c(0,100), ylims = c(0, 0.1))
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"density_estimation_beta", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(histograms[i], ncol = 1)
  ggsave(filename = paste(path_save,"density_estimation_beta_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}




# Densidad Espacio Logaritmo ----

n_utilizado = 1000
log_grid_vals = seq(-3, 8,length.out = 1000)
band_names <- c("Alpha", "Beta")

x_new_conditions <- data.frame(
  condition = cond_names,
  X_new = I(list(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,0,1,0,0,0),
    c(1,0,0,1,0,0),
    c(1,0,0,0,1,0),
    c(1,0,0,0,0,1)
  ))
)

# 3. y 4. Iterar y crear el data.frame final
density_df <- list()

for (cond in cond_names) {
  for (band in band_names) {
    
    # Selecciona la matriz y X_new correspondientes
    matrix_to_use <- if (band == "Alpha") {
      Matrix_Y1_burn
    } else {
      Matrix_Y2_burn
    }
    
    x_new_to_use <- x_new_conditions %>%
      filter(condition == cond) %>%
      pull(X_new) %>%
      `[[`(1)
    
    # Llama a la función
    result <- DDPM_Gaussian_Uni_Neal_dens(
      X_new = x_new_to_use,
      Matrix = matrix_to_use,
      grilla = log_grid_vals,
      n_iter = n_utilizado,
      stream = TRUE
    )
    # Crea un data.frame temporal con los resultados
    temp_df <- data.frame(
      x = log_grid_vals,
      mean_dens = rowMeans(result),
      lower_dens = apply(result, 1, function(x) quantile(x, 0.025)),
      upper_dens = apply(result, 1, function(x) quantile(x, 0.975)),
      condition = cond,
      Band = band
    )
    
    density_df[[paste(cond, band, sep = "_")]] <- temp_df
  }
}

# Combina los resultados
density_df <- bind_rows(density_df)


# Función para graficar histogramas con la densidad y sus intervalos
make_histogram <- function(cond, band, xlims, ylims) {
  
  # Filtra los datos del histograma
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  
  # Filtra los datos de la densidad y los intervalos
  density_sub <- density_df %>% filter(condition == cond, Band == band)
  
  ggplot(df_sub, aes(x = log(Power))) +
    geom_histogram(aes(y = after_stat(density)),
                   binwidth = .2,
                   fill = "gray30", color = "white") +
    
    # Añade la banda de credibilidad (Intervalo de Confianza)
    geom_ribbon(data = density_sub, aes(x = x, ymin = lower_dens, ymax = upper_dens), 
                fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
    
    # Añade la línea de la densidad media
    geom_line(data = density_sub, aes(x = x, y = mean_dens), 
              color = "blue", size = 0.5, inherit.aes = FALSE) +
    
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title = cond,
      x = expression(paste(Power, " (", mu, "V"^2, ")")),
      y = "Relative Frequency"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

## Banda Alpha ----
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Alpha", xlims = c(-3, 7), ylims = c(0,0.5))
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"density_estimation_log_alpha", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(histograms[i], ncol = 1)
  ggsave(filename = paste(path_save,"density_estimation_log_alpha_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}


## Banda Beta ----
histograms <- lapply(cond_names, function(cond) {
  list(
    make_histogram(cond, "Beta", xlims = c(-3, 7), ylims = c(0,0.7))
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(histograms, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"density_estimation_log_beta", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(histograms[i], ncol = 1)
  ggsave(filename = paste(path_save,"density_estimation_log_beta_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}





# Sobrevivencia Espacio Original ----
n_utilizado = 1000
grid_vals_surv = seq(10,500,length.out = 2000)

x_new_conditions <- data.frame(
  condition = cond_names,
  X_new = I(list(
    c(1,0,0,0,0,0),
    c(1,1,0,0,0,0),
    c(1,0,1,0,0,0),
    c(1,0,0,1,0,0),
    c(1,0,0,0,1,0),
    c(1,0,0,0,0,1)
  ))
)

# 3. y 4. Iterar y crear el data.frame final
surv_df <- list()

for (cond in cond_names) {
  for (band in band_names) {
    
    # Selecciona la matriz y X_new correspondientes
    matrix_to_use <- if (band == "Alpha") {
      Matrix_Y1_burn
    } else {
      Matrix_Y2_burn
    }
    
    x_new_to_use <- x_new_conditions %>%
      filter(condition == cond) %>%
      pull(X_new) %>%
      `[[`(1)
    
    # Llama a la función
    result <- DDPM_Gaussian_Uni_Neal_dist(
      X_new = x_new_to_use,
      Matrix = matrix_to_use,
      grilla = log(grid_vals_surv),
      n_iter = n_utilizado,
      stream = TRUE
    )
    # Crea un data.frame temporal con los resultados
    temp_df <- data.frame(
      x = grid_vals_surv,
      mean_dist = 1 - rowMeans(result),
      lower_dist = 1 - apply(result, 1, function(x) quantile(x, 0.025)),
      upper_dist = 1 - apply(result, 1, function(x) quantile(x, 0.975)),
      condition = cond,
      Band = band
    )
    
    surv_df[[paste(cond, band, sep = "_")]] <- temp_df
  }
}

# Combina los resultados
surv_df <- bind_rows(surv_df)


# Función para graficar histogramas con la densidad y sus intervalos
make_plot_surv <- function(cond, band, xlims, ylims) {
  
  # Filtra los datos del histograma
  df_sub <- plot_df %>% filter(condition == cond, Band == band)
  
  # Filtra los datos de la densidad y los intervalos
  surv_sub <- surv_df %>% filter(condition == cond, Band == band)
  
  ggplot(df_sub, aes(x = Power)) +
    # Agrega la capa geom_rug para mostrar la distribución de los datos df_sub
  geom_rug(sides = "b", alpha = 0.5) +
  # Añade la banda de credibilidad (Intervalo de Confianza)
  geom_ribbon(data = surv_sub, aes(x = x, ymin = log(lower_dist), ymax = log(upper_dist)), 
              fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  
  # Añade la línea de la densidad media
  geom_line(data = surv_sub, aes(x = x, y = log(mean_dist)), 
            color = "blue", size = 0.5, inherit.aes = FALSE) +
  
  coord_cartesian(xlim = xlims, ylim = ylims) +
  labs(
    title = cond,
    x = expression(paste(Power, " (", mu, "V"^2, ")")),
    y = "Log-Survival"
  ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
}

## Banda Alpha ----
pp <- lapply(cond_names, function(cond) {
  list(
    make_plot_surv(cond, "Alpha",
                   xlims = c(quantile(brainwave$Y1,probs = .9), quantile(brainwave$Y1,probs = .999)),
                   ylims = c(-8, -1)
                   )
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(pp, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"Log_Survival_estimation_alpha", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(pp[i], ncol = 1)
  ggsave(filename = paste(path_save,"Log_Survival_estimation_alpha_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}


## Banda Beta ----
pp <- lapply(cond_names, function(cond) {
  list(
    make_plot_surv(cond, "Beta",
                   xlims = c(quantile(brainwave$Y2,probs = .9), quantile(brainwave$Y2,probs = .999)),
                   ylims = c(-8, -1)
    )
  )
}) %>% unlist(recursive = FALSE)

# Organizar en 3 filas x 4 columnas
layout_hist <- wrap_plots(pp, ncol = 3, axes = "collect_y")

# Mostrar
print(layout_hist)

# Guardar 
path_save = "Img_Result/"
if (!dir.exists(path_save)) dir.create(path_save, recursive = TRUE)

ggsave(filename = paste(path_save,"Log_Survival_estimation_beta", ".png", sep= ""), plot = layout_hist,
       width = 10, height = 4.5, dpi = 300)


for(i in 1:6){
  plot = wrap_plots(pp[i], ncol = 1)
  ggsave(filename = paste(path_save,"Log_Survival_estimation_beta_", cond_names[i], ".png", sep= ""), plot = plot,
         width = 10, height = 4.5, dpi = 300)
}


# VGAM ----


back_BIC <- function(dat){
  # Ajustar el modelo completo
  modelo_completo <- vglm(y ~ x2 + x3 + x4 + x5 + x6 + 1,
                          frechet(location = 0, lscale = "loglink",
                                  lshape = "loglink"),
                          trace = F, data = dat)
  
  # Calcular el BIC del modelo completo
  bic_actual <- BIC(modelo_completo);bic_actual
  # Obtener los nombres de las variables predictoras
  variables <- colnames(X)[-1]
  
  # Bucle para realizar la selección hacia atrás
  repetir_proceso <- TRUE
  while (repetir_proceso) {
    bic_mejor_paso <- bic_actual
    variable_a_eliminar <- NULL
    
    # Probar la eliminación de cada variable
    for (variable in variables) {
      # Crear la fórmula sin la variable actual
      variables_restantes <- setdiff(variables, variable)
      formula_nueva <- paste("y ~ ", paste(variables_restantes, collapse = " + "), "+1")
      
      # Ajustar el nuevo modelo
      modelo_temporal <- vglm(as.formula(formula_nueva),
                              frechet(location = 0, lscale = "loglink",
                                      lshape = "loglink"),
                              trace = FALSE, data = dat)
      
      # Calcular el BIC
      bic_temporal <- BIC(modelo_temporal)
      #print(bic_temporal)
      
      # Comprobar si este modelo es mejor (menor BIC)
      if (bic_temporal < bic_mejor_paso) {
        bic_mejor_paso <- bic_temporal
        variable_a_eliminar <- variable
      }
    }
    
    # Si se encontró una variable cuya eliminación mejora el BIC...
    if (!is.null(variable_a_eliminar)) {
      # Actualizar el modelo y el BIC
      bic_actual <- bic_mejor_paso
      variables <- setdiff(variables, variable_a_eliminar)
      #print(paste("Se eliminó:", variable_a_eliminar, "Nuevo BIC:", bic_actual))
    } else {
      # Si ninguna eliminación mejora el BIC, detener el proceso
      repetir_proceso <- FALSE
      #print("El proceso de selección hacia atrás ha finalizado.")
    }
  }
  
  # El modelo final se puede ajustar con las 'variables' restantes
  formula_final <- paste("y ~", paste(variables, collapse = " + "), "+ 1")
  modelo_final <- vglm(as.formula(formula_final),
                       frechet(location = 0, lscale = "loglink",
                               lshape = "loglink"),
                       trace = F, data = dat)
  return(modelo_final)
  
}

dat = brainwave[c(1,3:8)]
colnames(dat)[1] = "y"
back_BIC(dat)

modelo_completo <- vglm(y ~ x2 + x3 + x4 + x5 + x6 + 1,
                        frechet(location = 0, lscale = "loglink",
                                lshape = "loglink"),
                        trace = T, data = dat)


