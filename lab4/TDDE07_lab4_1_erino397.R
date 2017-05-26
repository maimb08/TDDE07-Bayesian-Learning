require(MASS)
require(geoR)
require(mvtnorm)
require(LaplacesDemon)
library(xtable)
require(coda)

# -------
#  Lab 4
# -------

data = read.table("data/eBayNumberOfBidderData.dat", header=TRUE)


n = length(data)
n_features = ncol(data) - 1 # Except y and const

feature_labels = colnames(data[,2:ncol(data)])
  
y = data$nBids
X = as.matrix(data[,2:ncol(data)])

X_X = t(X)%*%X

# -----
#  (a)
# -----

glm_model = glm(nBids ~ 0 + ., data = data, family = poisson)

pdf("./plots/4_1_1_mle_beta.pdf", width=7, height=7)

par(oma = c(0, 0, 3, 0))
layout(matrix(c(0,1,1,0,2,3,4,5,6,7,8,9), 3, 4, byrow = TRUE))
for (i in 1:ncol(X)){
  mean = glm_model$coefficients[i]
  std_dev = summary(glm_model)[["coefficients"]][,2][i]
  x_grid = seq(mean-4*std_dev, mean+4*std_dev, 0.001)
  plot(x_grid,
       dnorm(x_grid, mean=mean, sd=std_dev),
       type="l",
       ylab="Density",
       xlab=expression(beta),
       main=feature_labels[i])
}
title("Normal approximation of MLE of beta", outer=TRUE, cex=1.5)

dev.off()

# -----
#  (b)
# -----

# Beta prior (Zellner’s g-prior)
mu0 = rep(0, n_features)
covar0 = 100 * ginv(X_X)
init_beta = mvrnorm(n=1, mu0, covar0)

# This is the log of the Poisson model
logPostPoiNorm <- function(betas, X, y){
  
  log_prior = dmvnorm(betas, mu0, covar0, log=TRUE)
  
  lambda = exp(X%*%betas)
  
  # Assume independence among samples and take the sum of
  # log(p(y_i|lambda)), where lambda is exp(X.dot(beta)) and p ~ Poisson
  log_lik = sum(dpois(y, lambda, log=TRUE))
  
  return (log_lik + log_prior)
}

log_post = logPostPoiNorm
opt_results = optim(init_beta,
                    log_post,
                    gr=NULL,
                    X,
                    y,
                    method=c("BFGS"),
                    control=list(fnscale=-1),
                    hessian=TRUE)

# MLE beta
post_mode = opt_results$par
# Covariance (J^-1(beta hat))
post_cov = -solve(opt_results$hessian)


# -----
#  (c)
# -----

Sigma = post_cov
c = .6

n_draws = 5000

metropolisHastings = function(logPostFunc, theta, c, ...){
  theta_draws = matrix(0, n_draws, length(theta))
  # Set initial 
  theta_c = mvrnorm(n=1, theta, c*Sigma) 
  prob_sum = 0
  accepted = 0
  for(i in 1:n_draws){
    # 1: Draw new proposal theta
    theta_p = mvrnorm(n=1, theta_c, c*Sigma)
    # 2: Determine the acceptance probability
    p_prev = logPostFunc(theta_c, ...)
    p_new = logPostFunc(theta_p, ...)
    acc_prob = min(c(1, exp(p_new - p_prev)))
    prob_sum = prob_sum + acc_prob
    # 3: Set new value with prob = acc_prob
    if(rbern(n=1, p=acc_prob)==1){
      theta_c = theta_p
      accepted = accepted + 1
    }
    theta_draws[i,] = theta_c
  }
  prob_acc = accepted / n_draws
  print(paste('Avg. acc. prob1. = ', round(prob_sum/n_draws, 2)))
  print(paste('Avg. acc. prob2. = ', round(prob_acc, 2)))
  
  return (theta_draws)
}

init_beta = mvrnorm(n=1, mu0, covar0)
beta_draws = metropolisHastings(logPostPoiNorm, init_beta, c, X, y)


# Calculate mean of batches of 2 draws to visualize the
# auto correlation between sequential draws
mean_draws = matrix(0, n_draws/2, n_features)
for (i in seq(2,n_draws,2)){
  mean_draws[i/2,] = colMeans(beta_draws[c(i-1,i),])
}

# Avoid first 10% of the draws
burn_in = floor(n_draws / 10)
beta_draws = beta_draws[burn_in:nrow(beta_draws),]

beta_means = colMeans(beta_draws)

pdf("./plots/4_1_2_beta_conv.pdf", width=7, height=7)

par(oma = c(0, 0, 3, 0))
layout(matrix(c(0,1,1,0,2,3,4,5,6,7,8,9), 3, 4, byrow = TRUE))
x_grid = 1:nrow(mean_draws)
for (i in 1:ncol(X)){
  # traceplot(mcmc(beta_draws[,i]),
  #      ylab="",
  #      xlab="Iteration",
  #      col="lightgray",
  #      main=feature_labels[i])
  plot(x_grid,
       mean_draws[,i],
       type="l",
       ylab="",
       xlab="Iteration",
       col="lightgray",
       main=feature_labels[i])
}
title("Convergence of beta during Metropolis Hastings", outer=TRUE, cex=1.5)

dev.off()

# -----
#  (d)
# -----

sample = c(
  Constant = 1,
  PowerSeller = 1,
  VerifyID = 1,
  Sealed = 1,
  MinBlem = 0,
  MajBlem = 0,
  LargNeg = 0,
  LogBook = 1,
  MinBidShare = 0.5
)

# Calculate lambda of pred. dens.
lambda = exp(beta_means%*%sample)

# Determine the predictive density of the sample
beta_grid = 0:max(y)
pred_dens = dpois(beta_grid, lambda)
names(pred_dens) = beta_grid

# Remove dependent variables that have prob. < .1%
pred_dens = pred_dens[pred_dens > .001]

# Probability that the sample has 0 bidders
prob = pred_dens[1]

pdf("./plots/4_1_3_pred_distr.pdf", width=5, height=4)

# Plot the predictive distribution
pred_plot = barplot(pred_dens,
                    col="white",
                    xaxt="n",
                    xlab="Number of bidders",
                    ylab="Probability",
                    main="Predictive distribution of sample")

axis(1, at=pred_plot, labels=names(pred_dens))

dev.off()