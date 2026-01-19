rm(list = ls())
library(VGAM)

# Betas y Matriz de diseño ----
set.seed(1000)
#####  Generalized pareto distribution ##
# Tamaño de muestra



step_AIC <- function(dat, variables_iniciales = NULL){
  # Todas las posibles variables predictoras
  todas <- colnames(dat)[colnames(dat) != "y"]
  
  # Variables iniciales: si no se especifican, comienza con todas
  if (is.null(variables_iniciales)) {
    variables <- todas
  } else {
    variables <- variables_iniciales
  }
  
  # Ajustar el modelo inicial
  formula_ini <- paste("y ~", paste(variables, collapse = " + "), "+ 1")
  modelo_actual <- vglm(as.formula(formula_ini),
                        frechet(location = 0, lscale = "loglink",
                                lshape = "loglink"),
                        trace = FALSE, data = dat, coefstart = rep(0, (length(variables) + 1)*2))
  bic_actual <- AIC(modelo_actual)
  
  repetir_proceso <- TRUE
  while (repetir_proceso) {
    bic_mejor_paso <- bic_actual
    mejor_accion <- NULL
    mejor_variable <- NULL
    
    # ----- 1. Intentar ELIMINAR cada variable -----
    if (length(variables) > 0) {
      for (var in variables) {
        vars_restantes <- setdiff(variables, var)
        if (length(vars_restantes) == 0) {
          formula_nueva <- "y ~ 1"
        } else {
          formula_nueva <- paste("y ~", paste(vars_restantes, collapse = " + "), "+ 1")
        }
        modelo_temp <- vglm(as.formula(formula_nueva),
                            frechet(location = 0, lscale = "loglink",
                                    lshape = "loglink"),
                            trace = FALSE, data = dat, coefstart = rep(0, (length(vars_restantes) + 1)*2))
        bic_temp <- AIC(modelo_temp)
        
        if (bic_temp < bic_mejor_paso) {
          bic_mejor_paso <- bic_temp
          mejor_accion <- "eliminar"
          mejor_variable <- var
        }
      }
    }
    
    # ----- 2. Intentar AGREGAR variables faltantes -----
    vars_fuera <- setdiff(todas, variables)
    if (length(vars_fuera) > 0) {
      for (var in vars_fuera) {
        vars_nuevas <- c(variables, var)
        formula_nueva <- paste("y ~", paste(vars_nuevas, collapse = " + "), "+ 1")
        modelo_temp <- vglm(as.formula(formula_nueva),
                            frechet(location = 0, lscale = "loglink",
                                    lshape = "loglink"),
                            trace = FALSE, data = dat, coefstart = rep(0, (length(vars_nuevas) + 1)*2))
        bic_temp <- AIC(modelo_temp)
        
        if (bic_temp < bic_mejor_paso) {
          bic_mejor_paso <- bic_temp
          mejor_accion <- "agregar"
          mejor_variable <- var
        }
      }
    }
    
    # ----- 3. Actualizar según el mejor movimiento -----
    if (!is.null(mejor_accion)) {
      if (mejor_accion == "eliminar") {
        variables <- setdiff(variables, mejor_variable)
      } else if (mejor_accion == "agregar") {
        variables <- c(variables, mejor_variable)
      }
      bic_actual <- bic_mejor_paso
      # print(paste("Se", mejor_accion, "la variable:", mejor_variable, "Nuevo BIC:", bic_actual))
    } else {
      repetir_proceso <- FALSE
      # print("El proceso stepwise BIC ha finalizado.")
    }
  }
  
  # Modelo final
  if (length(variables) == 0) {
    formula_final <- "y ~ 1"
  } else {
    formula_final <- paste("y ~", paste(variables, collapse = " + "), "+ 1")
  }
  
  modelo_final <- vglm(as.formula(formula_final),
                       frechet(location = 0, lscale = "loglink",
                               lshape = "loglink"),
                       trace = FALSE, data = dat,  coefstart = rep(0, (length(variables) + 1)*2))
  
  return(modelo_final)
}
step_BIC <- function(dat, variables_iniciales = NULL){
  # Todas las posibles variables predictoras
  todas <- colnames(dat)[colnames(dat) != "y"]
  
  # Variables iniciales: si no se especifican, comienza con todas
  if (is.null(variables_iniciales)) {
    variables <- todas
  } else {
    variables <- variables_iniciales
  }
  
  # Ajustar el modelo inicial
  formula_ini <- paste("y ~", paste(variables, collapse = " + "), "+ 1")
  modelo_actual <- vglm(as.formula(formula_ini),
                        frechet(location = 0, lscale = "loglink",
                                lshape = "loglink"),
                        trace = FALSE, data = dat, coefstart = rep(0, (length(variables) + 1)*2))
  bic_actual <- BIC(modelo_actual)
  
  repetir_proceso <- TRUE
  while (repetir_proceso) {
    bic_mejor_paso <- bic_actual
    mejor_accion <- NULL
    mejor_variable <- NULL
    
    # ----- 1. Intentar ELIMINAR cada variable -----
    if (length(variables) > 0) {
      for (var in variables) {
        vars_restantes <- setdiff(variables, var)
        if (length(vars_restantes) == 0) {
          formula_nueva <- "y ~ 1"
        } else {
          formula_nueva <- paste("y ~", paste(vars_restantes, collapse = " + "), "+ 1")
        }
        modelo_temp <- vglm(as.formula(formula_nueva),
                            frechet(location = 0, lscale = "loglink",
                                    lshape = "loglink"),
                            trace = FALSE, data = dat, coefstart = rep(0, (length(vars_restantes) + 1)*2))
        bic_temp <- BIC(modelo_temp)
        
        if (bic_temp < bic_mejor_paso) {
          bic_mejor_paso <- bic_temp
          mejor_accion <- "eliminar"
          mejor_variable <- var
        }
      }
    }
    
    # ----- 2. Intentar AGREGAR variables faltantes -----
    vars_fuera <- setdiff(todas, variables)
    if (length(vars_fuera) > 0) {
      for (var in vars_fuera) {
        vars_nuevas <- c(variables, var)
        formula_nueva <- paste("y ~", paste(vars_nuevas, collapse = " + "), "+ 1")
        modelo_temp <- vglm(as.formula(formula_nueva),
                            frechet(location = 0, lscale = "loglink",
                                    lshape = "loglink"),
                            trace = FALSE, data = dat, coefstart = rep(0, (length(vars_nuevas) + 1)*2))
        bic_temp <- BIC(modelo_temp)
        
        if (bic_temp < bic_mejor_paso) {
          bic_mejor_paso <- bic_temp
          mejor_accion <- "agregar"
          mejor_variable <- var
        }
      }
    }
    
    # ----- 3. Actualizar según el mejor movimiento -----
    if (!is.null(mejor_accion)) {
      if (mejor_accion == "eliminar") {
        variables <- setdiff(variables, mejor_variable)
      } else if (mejor_accion == "agregar") {
        variables <- c(variables, mejor_variable)
      }
      bic_actual <- bic_mejor_paso
      # print(paste("Se", mejor_accion, "la variable:", mejor_variable, "Nuevo BIC:", bic_actual))
    } else {
      repetir_proceso <- FALSE
      # print("El proceso stepwise BIC ha finalizado.")
    }
  }
  
  # Modelo final
  if (length(variables) == 0) {
    formula_final <- "y ~ 1"
  } else {
    formula_final <- paste("y ~", paste(variables, collapse = " + "), "+ 1")
  }
  
  modelo_final <- vglm(as.formula(formula_final),
                       frechet(location = 0, lscale = "loglink",
                               lshape = "loglink"),
                       trace = FALSE, data = dat,  coefstart = rep(0, (length(variables) + 1)*2))
  
  return(modelo_final)
}

p = 6  #5 predictores - 3 relevantes, 100 - 20 relevantes.
n_mc = 100
n_cluster = 1

for (n in c(200, 400, 800)){
  set.seed(1000)
  # Parámetros de la muestra ----

rademacher = sample(c(-1, 1), size = p*n_cluster, replace = TRUE)
betas_or <- matrix(rademacher*runif(p*n_cluster,.5,2), ncol = n_cluster) 
X = cbind(1,matrix(rnorm(n*(p-1), mean = 0, sd = 1), nrow = n))
colnames(X) = c("Intercept", "x1", "x2", "x3", "x4", "x5")
list_p_activos = as.matrix(cbind(1,expand.grid(rep(list(0:1), p-1))), ncol = p)
list_betas = list_p_activos *  rep(betas_or, each = nrow(list_p_activos))
sigma_true <- 1     # alpha = 1/xi
mu_true <- 0
for (iter_beta in 1:nrow(list_betas)){
  betas = t(list_betas[iter_beta,, drop = F])
  log_alpha <- X %*% betas
  alpha_true <- sqrt(exp(log_alpha))
  print(c(min(alpha_true), max(alpha_true)))
}




output_dir_AIC = paste0("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_AIC_n", n, "/")
if (!dir.exists(output_dir_AIC)) dir.create(output_dir_AIC, recursive = TRUE)
output_dir_BIC = paste0("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_BIC_n", n, "/")
if (!dir.exists(output_dir_BIC)) dir.create(output_dir_BIC, recursive = TRUE)


#1:nrow(list_betas)
for (iter_beta in 1: nrow(list_betas)){ 
  
  betas = t(list_betas[iter_beta,, drop = F])
  print(betas)
  log_alpha <- X %*% betas
  alpha_true <- sqrt(exp(log_alpha))
  
  # --- Ejecución paralela ---
  for(iter_mc in 1:n_mc){
    
    set.seed(iter_mc)
    y <- as.vector(rfrechet(n = n, location = mu_true, scale = sigma_true, shape = alpha_true))
    y[y<0.1] = 0.1
    dat = data.frame(
      y = y
    )
    # Código para ajustar el modelo
    modelo_final_AIC <- tryCatch(
      {
        step_AIC(dat = cbind(dat, X[,-1]))
      },
      # Bloque 'error': Se ejecuta si vglm falla
      error = function(e) {
        message("Se ha producido un error al ajustar el modelo vglm:")
        message(e)
        # Devolver NA o NULL para indicar que el modelo falló
        return(NULL)
      }
    )
    modelo_final_BIC <- tryCatch(
      {
        step_BIC(dat = cbind(dat, X[,-1]))
      },
      # Bloque 'error': Se ejecuta si vglm falla
      error = function(e) {
        message("Se ha producido un error al ajustar el modelo vglm:")
        message(e)
        # Devolver NA o NULL para indicar que el modelo falló
        return(NULL)
      }
    )
    save(modelo_final_AIC, file = paste0(output_dir_AIC, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata"))
    save(modelo_final_BIC, file = paste0(output_dir_BIC, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata"))
  }
}
}




# Creaeción Matrix ----
get_variable <- function(modelo){
  if (is.null(modelo)){return(NULL)}
  vars = names(coef(modelo))
  nombres <- unique(gsub(":.*", "", vars))
  nombres <- nombres[nombres != "(Intercept)"]
  ind = as.integer(sub("x", "", nombres)) + 1
  vec = integer(p)
  vec[c(1,ind)] = 1
  
  return(vec)
}

orden = order(rowSums(list_p_activos), decreasing = T)
M_ordenada <- list_p_activos[orden, ]
nombres <- apply(M_ordenada, 1, function(x) paste(x, collapse = ","))

for (n in c(200,400,800)){
  
  output_dir = paste0("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_BIC_n", n , "/")
  matrix_result = array(0, dim  = c(nrow(list_betas), nrow(list_betas) + 1))
  
  
  for (iter_beta in 1:nrow(list_betas)){
    print(iter_beta)
    result = c()
    for (iter_mc in 1:n_mc){
      
      file_path = paste0(output_dir, "result_", paste(list_p_activos[iter_beta,], collapse = ""), "_", iter_mc, ".Rdata")
      load(file = file_path)
      vec = get_variable(modelo_final_BIC)
      if (is.null(vec)){
        
        matrix_result[iter_beta, nrow(list_betas) + 1] = matrix_result[iter_beta,nrow(list_betas) + 1] + 1
      }
      else{
      matrix_result[iter_beta, which(apply(list_p_activos, 1, function(x) all(x == vec)))] = matrix_result[iter_beta, which(apply(list_p_activos, 1, function(x) all(x == vec)))] + 1
      }
    }
    matrix_result[iter_beta, ] = matrix_result[iter_beta, ]/sum(matrix_result[iter_beta, ])
  }
  matrix_result = matrix_result[orden, c(rev(orden), 33)]
  rownames(matrix_result) <- nombres
  colnames(matrix_result) <- c(rev(nombres), "NULL")
  save(matrix_result, file = paste0(output_dir, "Matrix_n", n, ".Rdata"))
  
}

# Extraer nombres sin ":n"

load("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_BIC_n200/Matrix_n200.Rdata")
matrix_result_n200 = matrix_result
load("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_BIC_n400/Matrix_n400.Rdata")
matrix_result_n400 = matrix_result
load("../Example_Frechet_shape/MC_study_p5/MC_study_VGAM_BIC_n800/Matrix_n800.Rdata")
matrix_result_n800 = matrix_result

library(ggplot2)
library(reshape2)
library(patchwork)  # Para disposición horizontal elegante

# Suponiendo que matrix_result_n200 y matrix_result_n400 son matrices ya cargadas
# Convertir ambas matrices a formato largo
df1 <- melt(matrix_result_n200)
df2 <- melt(matrix_result_n400)
df3 <- melt(matrix_result_n800)

# Nombres para los ejes, como factores para mantener el orden
df1$Var1 <- factor(df1$Var1, levels = rev(rownames(matrix_result_n200)))
df1$Var2 <- factor(df1$Var2, levels = colnames(matrix_result_n200))

df2$Var1 <- factor(df2$Var1, levels = rev(rownames(matrix_result_n400)))
df2$Var2 <- factor(df2$Var2, levels = colnames(matrix_result_n400))

df3$Var1 <- factor(df3$Var1, levels = rev(rownames(matrix_result_n800)))
df3$Var2 <- factor(df3$Var2, levels = colnames(matrix_result_n800))


# Crear ambos plots
p1 <- ggplot(df1, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.001)) +
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
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.001)) +
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
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.001)) +
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
path_save = "../Example_Frechet_shape/MC_study_p5/Img_Study/"
ggsave(filename = paste(path_save,"Pixel_plot_VGAM_BIC", ".png", sep= ""), plot = final_plot,
       width = 10, height = 4.5, dpi = 300)

