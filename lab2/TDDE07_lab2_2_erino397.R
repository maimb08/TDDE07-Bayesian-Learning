require(MASS)
require(geoR)
library(mvtnorm)

# Lab 2 - Assignment 2

df = read.table("data/WomenWork.dat.txt", header=TRUE)

# (a)

glmModel <- glm(Work ~ 0 + ., data = df,family = binomial)
summary(glmModel)

# (b)

y = df[["Work"]]
X = as.matrix(df[,2:9])
headers = colnames(df)
n_params = dim(X)[2]

tau_sq = 100 # tau = 10

# Initialize prior hyper parameters
mu_prior <- as.vector(rep(0,n_params))
sigma_prior = diag(tau_sq, n_params, n_params)

# Calculate the log posterior
LogPosteriorLogistic <- function(betas, y, X, mu, Sigma){
  
  # Multiply data by parameters to get predictions
  predictions <- X%*%betas;
  
  # Log likelihood                                   
  log_likelihood <- sum(predictions*y - log(1 + exp(predictions)));
  if (abs(log_likelihood) == Inf) log_likelihood = -20000;
  
  # Log prior
  log_prior <- dmvnorm(betas, mu_prior, sigma_prior, log=TRUE);
  
  # Sum of log likelihood and log prior is log posterior
  return(log_likelihood + log_prior)
}

log_posterior = LogPosteriorLogistic

# Initialize as zeros
init_betas = rep(0, n_params)

opt_results = optim(init_betas,
                    log_posterior,
                    gr=NULL,
                    y,
                    X,
                    mu_prior,
                    sigma_prior,
                    method=c("BFGS"),
                    control=list(fnscale=-1),
                    hessian=TRUE)

# Posterior mode (beta hat)
post_mode = opt_results$par
# Posterior covariance (J^-1(beta hat))
post_cov = -solve(opt_results$hessian)
approx_post_std_dev = sqrt(diag(post_cov))

# Plots of the marginal distributions of the parameters
# par(mfrow = c(2,2))
# for (k in 1:7){
#   beta_grid <- seq(0, post_mode[k] + 4*approx_post_std_dev[k], length = 1000)
#   plot(beta_grid, dnorm(x = beta_grid, mean = post_mode[k], sd = approx_post_std_dev[k]), type = "l", lwd = 2, main = names(post_mode)[k], ylab = '', xlab = headers[k+1])
# }

# Plot NSmallChild parameter
pmode = post_mode[7]
pstd = approx_post_std_dev[7]
beta_grid = seq(pmode - 4*pstd, pmode + 4 * pstd, length=1000)
eti = qnorm(c(0.025, 0.975), pmode, pstd)
dn = dnorm(x=beta_grid, mean=pmode, sd=pstd)
plot(beta_grid, dn, type = "l", lwd = 2, main="ETI for NSmallChild parameter", ylab = '', xlab=headers[8])
lines(eti, rep(0.04, 2), col="black", lwd=2)

# (c)

# Sample to predict
husband = 1
edu_years = 8
exp_years1 = 10
exp_years2 = (exp_years1/10)^2
age = 40
n_small_child = 1
n_big_child = 1

sample = c(1,
           husband,
           edu_years,
           exp_years1,
           exp_years2,
           age,
           n_small_child,
           n_big_child)

# Normal prior and normal posterior => Normal predictive distribution
sample_cov = diag(1, n_params) # ?
mu_n = post_mode
covar_n = sample_cov + post_cov

distr = dmvnorm(sample, mean=mu_n, sigma=covar_n)
