require(mvtnorm)
require(msm)
require(MASS)
library(LaplacesDemon)

grid_w = 6
grid_h = 5

# -------
#  Lab 3
# -------


# Read data
data <- read.table("./data/WomenWork.dat", header = TRUE)

feature_labels = colnames(data)

y = as.vector(data$Work)
X = as.matrix(data[, 2:ncol(data)])

# Data spec.
n_features = ncol(X)
n_samples = nrow(X)

tau = 10


# -------------
#  (a) and (b)
# -------------


# Beta prior parameters
mu0 = rep(0, n_features)
covar0 = diag(tau^2, n_features)

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
      # Truncate [-inf, 0)
      u[i] = rtnorm(n=1, mean=regr_mean[i], sd=1, upper=0)
    }else{
      # Truncate (0, inf]
      u[i] = rtnorm(n=1, mean=regr_mean[i], sd=1, lower=0)
    }
  }
  
  return(u)
}

# Initial prediction
u = rnorm(n_samples, covar0)

n_draws = 1500
beta_draws = matrix(0, n_draws, n_features)
u_draws = matrix(0, n_draws, n_samples)
for(i in 1:n_draws) {
  beta = draw_beta(u)
  u = draw_u(beta)
  beta_draws[i,] = beta
  u_draws[i,] = u
}

# Avoid first 10% of the draws
burn_in = floor(n_draws / 10)
beta_draws = beta_draws[burn_in:nrow(beta_draws),]

# -----
#  (c)
# -----


# Calculate the log posterior
LogPosteriorProbit <- function(betas, y, X, mu, Sigma){
  
  # Multiply data by parameters to get predictions
  predictions <- X%*%betas;
  
  # Log likelihood (for probit)                             
  log_likelihood = sum(y*pnorm(predictions, log.p = TRUE) + 
                         (1-y)*pnorm(predictions, log.p = TRUE, lower.tail = FALSE))
  
  # Log prior
  log_prior <- dmvnorm(betas, mu0, covar0, log=TRUE);
  
  # Sum of log likelihood and log prior is log posterior
  return(log_likelihood + log_prior)
}

log_posterior = LogPosteriorProbit

# Initialize as zeros
init_betas = rep(0, n_features)

opt_results = optim(init_betas,
                    log_posterior,
                    gr=NULL,
                    y,
                    X,
                    mu0,
                    covar0,
                    method=c("BFGS"),
                    control=list(fnscale=-1),
                    hessian=TRUE)

# Posterior mode (beta hat)
post_mode = opt_results$par
# Posterior covariance (J^-1(beta hat))
post_cov = -solve(opt_results$hessian)

pdf("plots/3_2_3_norm_gibbs_comp.pdf")

par(mfrow=c(4,2))

beta_grid = seq(-1.5, 1.5, 0.001)
for (i in 1:n_features) {
  
  # Build histogram of Gibbs draws of beta_i
  h = hist(beta_draws[,i], breaks=30, plot=FALSE)
  
  # Get the normal approximation for beta_i via optim.
  mean = post_mode[i]
  std_dev = sqrt(post_cov[i,i]/n_features)
  norm_approx = dnorm(x=beta_grid, mean=mean, sd=std_dev)
  
  # Find x- and y-limits for plot
  min_x = min(c(min(beta_draws[,i]), mean-4*std_dev))
  max_x = max(c(max(beta_draws[,i]), mean+4*std_dev))
  max_y = max(c(max(norm_approx), max(h$density)))
  
  # Plot the histogram
  plot(h,
       freq=FALSE,
       xlim=c(min_x,max_x), 
       ylim=c(0,max_y),
       border="gray",
       xlab=expression(theta),
       main=feature_labels[i+1])
  
  # Plot the normal approximation
  lines(beta_grid, norm_approx)
}

dev.off()
