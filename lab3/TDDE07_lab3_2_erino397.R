require(mvtnorm)
require(msm)
library(LaplacesDemon)

grid_w = 6
grid_h = 5

# -------
#  Lab 3
# -------


# Read data
data <- read.table("./data/WomenWork.dat", header = TRUE)

y <- as.vector(data$Work)
X <- as.matrix(data[, 2:ncol(data)])

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

n_draws = 500
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

beta_mean = colMeans(beta_draws)
beta_cov = cov(beta_draws)

# -----
#  (c)
# -----


# Calculate the log posterior
LogPosteriorProbit <- function(betas, y, X, mu, Sigma){
  
  # Multiply data by parameters to get predictions
  predictions <- X%*%betas;
  
  # Log likelihood (for probit)                             
  log_likelihood = sum(y*pnorm(predictions, log.p = TRUE) + (1-y)*pnorm(predictions, log.p = TRUE, lower.tail = FALSE))
  
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

# Sample to predict
sample = c(constant=1,
           husband=10,
           edu_years=8,
           exp_years1=10,
           exp_years2=(10/10^2),
           age=40,
           n_small_child=1,
           n_big_child=1)

get_pred = function(beta){
  e = exp(sample%*%beta)
  # Calculate the probability (bernoulli parameter)
  p = e / (1 + e)
  # Draw a y prediction
  y_draw = rbern(n=1, prob=p)
}

n_draws = 100
y_draws_2 = c() # As in lab 2
y_draws_3 = c() # As in lab 3
for (i in 1:n_draws){
  # Get prediction according to lab 2
  beta1 = as.vector(rmvnorm(n=1, mean=post_mode, sigma=post_cov))
  y_draw = get_pred(beta1)
  y_draws_2 = c(y_draws_2, y_draw)

  # Get prediction according to lab 3
  beta2 = as.vector(rmvnorm(n=1, mean=beta_mean, sigma=beta_cov))
  y_draw = get_pred(beta1)
  y_draws_3 = c(y_draws_3, y_draw)
}

pdf("plots/3_2_3_pred_distr_comp.pdf", width=grid_w, height=grid_h)

prob_density = density(y_draws_2)
plot(prob_density,
     type="l",
     lwd=2,
     xlim=c(0,1),
     ylim=c(0,2.5),
     ylab="Density",
     xlab="Prediction (0 = not working, 1 = working)",
     main="Predictive distribution for sample",
     col="black",
     cex.main=.9,
     cex.lab=.9,
     cex.axis=.8)

prob_density = density(y_draws_3)
lines(prob_density, col="gray", lwd=2)

legend("topright", 
       legend = c("Normal Approximation","Gibbs Probit"),
       fill = c("black", "gray"),
       inset = 0.02)

dev.off()