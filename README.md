# Example: Univariate DDP-GMM with Variable Selection (Neal Algorithm)

This example illustrates the implementation of a **Dirichlet Dependent
Process Gaussian Mixture Model (DDP-GMM)** for univariate responses with
**Bayesian variable selection**, using Neal's sampling algorithm. The
script demonstrates data simulation, model fitting, posterior inference,
and predictive density estimation.

------------------------------------------------------------------------

## Requirements

Before running the example, install and load the following R libraries:

``` r
library(progress)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
```

------------------------------------------------------------------------

## Project Structure

    ├── DDPGMM_Uni_Neal_Selection.R
    ├── Example_Gauss
    │   └── Ex_DDPGMM_Uni_Neal_Selection.R

-   **DDPGMM_Uni_Neal_Selection.R**\
    Contains the main implementation of the Bayesian nonparametric
    model.

-   **Example_Gauss/Ex_DDPGMM_Uni_Neal_Selection.R**\
    Provides a complete working example including data generation, model
    fitting, and visualization.

------------------------------------------------------------------------

## Running the Example

From the `Example_Gauss` directory, run:

``` r
source("../DDPGMM_Uni_Neal_Selection.R")
source("Ex_DDPGMM_Uni_Neal_Selection.R")
```

------------------------------------------------------------------------

## Example Description

The script performs the following steps:

### 1. Simulated Data Generation

Observations are generated from a mixture of univariate Gaussian
regression models. The data are constructed using:

-   Regression parameters (`betas`)
-   Component variances (`sigmas2`)
-   Mixture proportions (`proporciones`)
-   Covariate matrix (`X`)
-   Response variable (`Y`)

The total sample size is controlled by:

``` r
betas <- matrix(c(c(1,10), c(-3,14)), ncol = 2, byrow = F)
sigmas2 <- c(.5, 1)
proporciones <- c(0.6, 0.4)
n <- 8e2
p = dim(betas)[1]
s = sample(c(0,1),n, replace = T)
X <- model.matrix(~ as.factor(s) - 1)
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
```
------------------------------------------------------------------------

### 3. Visualization Data


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


------------------------------------------------------------------------

### 3. Hyperparameter Specification

The Bayesian model uses prior distributions defined through:

``` r
param_init = list(
  mu_0 = mu_0,
  tau = tau,
  a_sigma2 = a_sigma2,
  b_sigma2 = b_sigma2,
  a_alpha = a_alpha,
  b_alpha = b_alpha,
  H_max = H_max,
  n_cluster_ini = n_cluster_ini,
  zeta = zeta,
  gamma_inicio = gamma_inicio
)
```

------------------------------------------------------------------------

### 3. MCMC Configuration

The Markov Chain Monte Carlo (MCMC) settings are defined as:

``` r
param_chain = list(
  burnit = burnit,
  n_chain = n_chain,
  thining = thining
)
```

------------------------------------------------------------------------

### 4. Model Fitting

The model is estimated using:

``` r
Matrix = DDPM_Gaussian_Uni_Neal_Selection(
  Y, X,
  param_chain = param_chain,
  param_init = param_init,
  Y_scale = TRUE,
  stream = TRUE
)
```

The object `Matrix` contains posterior samples and latent model
structures.

------------------------------------------------------------------------

### 5. Variable Selection

Variable selection is performed through the binary vector `gamma`. The
example includes:

-   Frequency analysis of selected models
-   Identification of the modal model
-   Detection of relevant covariates

------------------------------------------------------------------------

### 6. Predictive Inference

#### Predictive Sampling

``` r
y_pred = DDPM_Gaussian_Uni_Neal_muestra(
  X_new = X_new1,
  Matrix = Matrix,
  n_muestra = 1e4
)
```

#### Predictive Density Estimation

``` r
dens1 = DDPM_Gaussian_Uni_Neal_dens(
  X_new = X_new1,
  Matrix,
  grilla = grid,
  n_iter = n_utilizado
)
```

The example compares:

-   Theoretical density
-   Estimated posterior density
-   Credible intervals

------------------------------------------------------------------------

## Visualization

The script produces graphical diagnostics including:

-   Cluster number evolution
-   Variable selection summaries
-   Comparison between empirical, theoretical, and estimated densities
-   Credible interval bands

------------------------------------------------------------------------

## Customization

Users can easily modify:

-   Number of mixture components
-   Regression dimension
-   Bayesian hyperparameters
-   MCMC configuration
-   Prediction covariates

------------------------------------------------------------------------

## Saving and Loading Results

Results can be saved using:

``` r
save(Matrix, file = "results.RData")
```

And loaded later with:

``` r
load("results.RData")
```

------------------------------------------------------------------------

## Methodological Background

This implementation is based on Bayesian nonparametric mixture models
using dependent Dirichlet processes and Neal-type MCMC sampling
strategies.
