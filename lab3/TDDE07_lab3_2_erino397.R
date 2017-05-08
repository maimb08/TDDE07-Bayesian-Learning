require(mvtnorm)
require(msm)

grid_w = 5
grid_h = 4

# Lab 3 - Assignment 2

# Read data
data <- read.table("./data/WomenWork.dat", header = TRUE)

y <- as.vector(data$Work)
X <- as.matrix(data[, 2:ncol(data)])

# Data spec.
n_features = ncol(X)
n_samples = nrow(X)

tau = 10

# Beta prior parameters
mu0 = rep(0, n_features)
covar0 = tau^2 * diag(n_features)

draw_beta = function(y) {
  
  X_X = t(X) %*% X
  
  # Least squares approximate of beta
  beta_hat = ginv(X_X) %*% t(X) %*% y
  
  # Posterior parameters
  mu_n = ginv(X_X + covar0)%*%(X_X%*%beta_hat+covar0%*%mu0)
  covar_n = X_X + covar0
  
  # Assuming sigma_sq = 1
  b_draw = rmvnorm(1, mean=mu_n, sigma=ginv(covar_n))
  
  return(b_draw)
}

draw_u = function(beta) {
  
  # Mean of predictive distr.
  regr_mean = X %*% t(beta)
  
  u = rep(0, n_samples)
  for(i in 1:n_samples) {
    y_i = y[i]
    
    if(y_i == 0){
      # Truncate [-inf, 0]
      u[i] = rtnorm(n=1, mean=regr_mean[i], sd=1, upper=0)
    }else{
      # Truncate [0, inf]
      u[i] = rtnorm(n=1, mean=regr_mean[i], sd=1, lower=0)
    }
  }
  
  return(u)
}

# Initial prediction
u = rnorm(n_samples, covar0)

n_draws = 100
beta_draws = matrix(0, n_draws, n_features)
u_draws = matrix(0, n_draws, n_samples)
for(i in 1:n_draws) {
  beta = draw_beta(u)
  u = draw_u(beta)
  beta_draws[i,] = beta
  u_draws[i,] = u
}

# Predict as y_i = 1 if u_i > 0, else y_i = 0
y_pred = rep(0, n_samples)
y_pred[colMeans(u_draws) > 0] = 1


