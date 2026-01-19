rm(list = ls())
library(VGAM)


load(file = "brainwave.RData")



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

# Datos

X = as.matrix(brainwave[,4:8])
Y1 = brainwave$Y1
Y2 = brainwave$Y2


dat1 = data.frame(
  y = Y1
)
dat2 = data.frame(
  y = Y2
)

# Código para ajustar el modelo
modelo1_final_AIC <- tryCatch(
  {
    step_AIC(dat = cbind(dat1, X))
  },
  # Bloque 'error': Se ejecuta si vglm falla
  error = function(e) {
    message("Se ha producido un error al ajustar el modelo vglm:")
    message(e)
    # Devolver NA o NULL para indicar que el modelo falló
    return(NULL)
  }
)
modelo1_final_BIC <- tryCatch(
  {
    step_BIC(dat = cbind(dat1, X))
  },
  # Bloque 'error': Se ejecuta si vglm falla
  error = function(e) {
    message("Se ha producido un error al ajustar el modelo vglm:")
    message(e)
    # Devolver NA o NULL para indicar que el modelo falló
    return(NULL)
  }
)

# Código para ajustar el modelo
modelo2_final_AIC <- tryCatch(
  {
    step_AIC(dat = cbind(dat2, X))
  },
  # Bloque 'error': Se ejecuta si vglm falla
  error = function(e) {
    message("Se ha producido un error al ajustar el modelo vglm:")
    message(e)
    # Devolver NA o NULL para indicar que el modelo falló
    return(NULL)
  }
)
modelo2_final_BIC <- tryCatch(
  {
    step_BIC(dat = cbind(dat2, X))
  },
  # Bloque 'error': Se ejecuta si vglm falla
  error = function(e) {
    message("Se ha producido un error al ajustar el modelo vglm:")
    message(e)
    # Devolver NA o NULL para indicar que el modelo falló
    return(NULL)
  }
)


