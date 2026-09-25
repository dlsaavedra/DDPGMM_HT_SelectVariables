rm(list = ls())

source("../DDPGMM_Uni_Neal_Selection.R")
set.seed(1234)

# # Crear Muestra ------
# Ejemplo de uso Normal Univariada
# Definir medias
#betas <- matrix(c(1,1,0,5), ncol = 1, byrow = F)
betas <- matrix(c(c(1,10), c(-3,14)), ncol = 2, byrow = F)
#betas <- matrix(c(c(1,10,0), c(-3,14,0), c(-8, 18, 0)), ncol = 3, byrow = F)
#betas <- matrix(c(c(5,10,7),c(15, 20, 5),c(-20,-15, 3),c(-10, -5, 1)), ncol = 4 , byrow = F)
# Definir matrices de covarianza
#sigmas2 <- c(1)
sigmas2 <- c(.5, 1)
#sigmas2 <- c(.5,.8, .4)
#sigmas2 <- c(1,2,3,4)
# Definir proporciones
#proporciones <- c(1)
proporciones <- c(0.6, 0.4)
#proporciones <- c(0.4, 0.3, 0.3)
#proporciones <- c(0.4, 0.3,0.2,0.1)
# Número total de muestras a generar
n <- 8e2
p = dim(betas)[1]
s = sample(c(0,1),n, replace = T)
X <- model.matrix(~ as.factor(s) - 1)#matrix(rnorm(n * (p-1)), nrow = n))
#X = cbind(X, rnorm(n), rnorm(n))
#X = cbind(1,matrix(rnorm(n * (p-1)), nrow = n))

# Calcular los tamaños de los subgrupos
ni <- round(proporciones * n)
# Ajustar el último tamaño para asegurar que sume n
ni[length(ni)] <- n - sum(ni[-length(ni)])
aux_ni = c(0,cumsum(ni))
# Generar muestras para cada componente
muestras = list()
indicadora <- c()
for (i in 1:dim(betas)[2]) {
  muestras[[i]] <- rnorm(ni[i], mean = X[(aux_ni[i] + 1) :aux_ni[i +1],]%*%betas[,i], sd = sqrt(sigmas2[i]))
  indicadora <- c(indicadora, rep(i,ni[i]))
}
Y = unlist(muestras)

hist(Y[X[,1]==1], freq = F, breaks = 20, ylim = c(0, 0.5), xlim = c(-15, 25))
X_new1 = c(1,0) #rep(-.1,9))
mu_new = 5
dens_teo1 = Vectorize(function (x) {densidad_mezcla_normales_uni_X(x, X_new1, betas, sigmas2, proporciones)})
grid = seq(mu_new - 20, mu_new + 20, length.out = 2000)
lines(grid,dens_teo1(grid), col = 3, type = "l", lwd = 3)

hist(Y[X[,1]==0], freq = F, breaks = 50, ylim = c(0, 0.5), xlim = c(-15, 25), add = T)
X_new2 = c(0,1) #rep(-.1,9))
mu_new = 5
dens_teo2 = Vectorize(function (x) {densidad_mezcla_normales_uni_X(x, X_new2, betas, sigmas2, proporciones)})
grid = seq(mu_new - 20, mu_new + 20, length.out = 2000)
lines(grid,dens_teo2(grid), col = 4, type = "l", lwd = 3)




# Hiperparametros ----
mu_0 = rep(0, p)
tau = 100
a_sigma2 = 2 #2
b_sigma2 = .01 #.5 # E(sigma^2)= b/(a-1)
a_alpha = .1
b_alpha = .1
H_max = 10
n_cluster_ini = 10

zeta = 1
gamma_inicio = c(1, rep(1,p-1))
#gamma_inicio = c(1, rep(0,p-1))


param_init = list(mu_0 = mu_0, tau = tau,
                  a_sigma2 = a_sigma2, b_sigma2 = b_sigma2, 
                  a_alpha = a_alpha, b_alpha = b_alpha, H_max = H_max,
                  n_cluster_ini = n_cluster_ini,
                  zeta = zeta, gamma_inicio = gamma_inicio)

# Parametros de la cadena ----
burnit = 2000
n_chain = 20000
thining = 2
param_chain = list(burnit= burnit, n_chain = n_chain, thining = thining )



#library(profvis)
#profvis({
Matrix = DDPM_Gaussian_Uni_Neal_Selection(Y, X, 
                                             param_chain = param_chain, 
                                             param_init = param_init, Y_scale = T, 
                                          stream = T)
#})


ts.plot(Matrix$H)
ts.plot(colSums(Matrix$gamma))
rowMeans(Matrix$gamma)
betas

vec_str <- apply(Matrix$gamma, 2, paste, collapse = ",")

# Contamos la frecuencia de cada cadena
tab <- table(vec_str);tab

# Extraemos la cadena más frecuente
moda_str <- names(which.max(tab))

# Convertimos la cadena de nuevo a vector numérico
moda_vec <- as.numeric(strsplit(moda_str, ",")[[1]])

# Resultado:
print(rbind(Posición = 1:length(moda_vec), Gamma = moda_vec))
which(moda_vec == 1)
#sort(p_activos)


# Densidad dado X -----



Y = Matrix$Y * Matrix$sd_Y_or + Matrix$mean_Y_or
X_new1 = c(1,0,0)#rep(-.1,9))
mu_new = 5
dens_teo1 = Vectorize(function (x) {densidad_mezcla_normales_uni_X(x, X_new1, betas, sigmas2, proporciones)})
grid = seq(mu_new - 20, mu_new + 20, length.out = 2000)

y_pred = DDPM_Gaussian_Uni_Neal_muestra(X_new = X_new1, 
                                        Matrix = Matrix, n_muestra = 1e4)
hist(y_pred[X[,1] == 1], breaks = 100, freq = F, xlim = c(min(Y)-5, max(Y)+5))
lines(grid,dens_teo1(grid), col = 3, type = "l", lwd = 3)





l_secuencia = dim(Matrix$alpha)
n_utilizado = 1000
dens1 = DDPM_Gaussian_Uni_Neal_dens(X_new = X_new1,
                                    Matrix, grilla = grid,
                                    n_iter = n_utilizado)

X_new2 = c(0,1,0)
dens_teo2 = Vectorize(function (x) {densidad_mezcla_normales_uni_X(x, X_new2, betas, sigmas2, proporciones)})

dens2 = DDPM_Gaussian_Uni_Neal_dens(X_new = X_new2,
                                    Matrix, grilla = grid,
                                    n_iter = n_utilizado)

mean_dens1 = rowMeans(dens1)
lower_dens1 <- apply(dens1, 1, quantile, probs = 0.025)
upper_dens1 <- apply(dens1, 1, quantile, probs = 0.975)


mean_dens2 = rowMeans(dens2)
lower_dens2 <- apply(dens2, 1, quantile, probs = 0.025)
upper_dens2 <- apply(dens2, 1, quantile, probs = 0.975)


#hist(y_pred, freq = F, breaks = 1000,xlim = c(min(Y_or)-10, max(Y_or)+10))
hist(Y[X[,1] == 1], freq = F, breaks = 20,xlim = c(-15,25), ylim = c(0,.5), 
     main = "X = c(1,0,10)", xlab = "y")
lines(grid, dens_teo1(grid), col  = "red", lty = 1, lwd= 3)
lines(grid,mean_dens1, col = "blue", lty = 2, lwd = 3, xlim = c(min(Y)-5, max(Y)+5), ylim = c(0,1))

# Intervalo de credibilidad (sombreado)
polygon(c(grid, rev(grid)),
        c(lower_dens1, rev(upper_dens1)),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)

legend("topright", legend=c("Muestra Empírica", "Densidad Teórica", "Densidad Estimada") ,
       col=c(1,2,"blue"), lwd=2, lty= c(1,1), cex = 0.6)

hist(Y[X[,1] == 0], freq = F, breaks = 20,xlim = c(-10, 25), 
     ylim = c(0,.5), main = "X = c(0,1,10)", xlab = "y")
lines(grid, dens_teo2(grid), col  = "red", lty = 1, lwd= 3)
lines(grid,mean_dens2, col = "blue", lty = 2, lwd = 3, xlim = c(min(Y)-5, max(Y)+5), ylim = c(0,1))

# Intervalo de credibilidad (sombreado)
polygon(c(grid, rev(grid)),
        c(lower_dens2, rev(upper_dens2)),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)
legend("topright", legend=c("Muestra Empírica", "Densidad Teórica", "Densidad Estimada") ,
       col=c(1,2,"blue"), lwd=2, lty= c(1,1), cex = 0.6)


#save(Matrix, file = "../Example_Gauss/n100_2modas_Selection.RData")
#load("../Example_Gauss/n1000_2modas_Selection.RData")
