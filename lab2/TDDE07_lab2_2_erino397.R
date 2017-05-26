require(MASS)
require(geoR)
library(mvtnorm)
library(LaplacesDemon)


grid_w = 5
grid_h = 4


# Lab 2 - Assignment 2

df = read.table("data/WomenWork.dat", header=TRUE)

# (a)

glmModel <- glm(Work ~ 0 + ., data = df,family = binomial)
summary(glmModel)

# (b)

y = df[["Work"]]
X = as.matrix(df[,2:9])
headers = colnames(df)
n = dim(X)[1]
n_params = dim(X)[2]

tau = 10

# Initialize prior hyper parameters
mu0 = as.vector(rep(0,n_params))
sigma_sq0 = diag(tau^2, n_params, n_params)

# Calculate the log posterior
LogPosteriorLogistic <- function(betas, y, X, mu, Sigma){
  
  # Multiply data by parameters to get predictions
  predictions <- X%*%betas;
  
  # Log likelihood                                   
  log_likelihood <- sum(predictions*y - log(1 + exp(predictions)));
  if (abs(log_likelihood) == Inf) log_likelihood = -20000;
  
  # Log prior
  log_prior <- dmvnorm(betas, mu0, sigma_sq0, log=TRUE);
  
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
                    mu0,
                    sigma_sq0,
                    method=c("BFGS"),
                    control=list(fnscale=-1),
                    hessian=TRUE)

# Posterior mode (beta hat)
post_mode = opt_results$par
# Posterior covariance (J^-1(beta hat))
post_cov = -solve(opt_results$hessian)
approx_post_std_dev = sqrt(diag(post_cov))

pdf("plots/2_2_nsmallchild_cred_interval.pdf", width=grid_w, height=grid_h)

# Plot NSmallChild parameter
pmode = post_mode[7]
pstd = approx_post_std_dev[7]
beta_grid = seq(pmode - 4*pstd, pmode + 4 * pstd, length=1000)
eti = qnorm(c(0.025, 0.975), pmode, pstd)
dn = dnorm(x=beta_grid, mean=pmode, sd=pstd)
plot(beta_grid, 
     dn, 
     type = "l", 
     lwd = 2, 
     main="ETI for NSmallChild parameter", 
     ylab = 'Density', xlab=headers[8],
     cex.main=.9, 
     cex.lab=.9, 
     cex.axis=.8)
lines(eti, rep(0.04, 2), col="black", lwd=2)

dev.off()

# (c)

# Sample to predict
sample = c(constant=1,
           husband=10,
           edu_years=8,
           exp_years1=10,
           exp_years2=(10/10^2),
           age=40,
           n_small_child=1,
           n_big_child=1)

y_draws = c()
n_draws = 3000
for (i in 1:n_draws){
  # Draw a beta
  beta_draw = as.vector(rmvnorm(n=1, mean=post_mode, sigma=post_cov))
  e = exp(sample%*%beta_draw)
  # Calculate the probability (bernoulli parameter)
  p = e / (1 + e)
  # Draw a y prediction
  y_draw = rbern(n=1, prob=p)
  y_draws = c(y_draws, y_draw)
}

outcomes = table(y_draws)
n_working = outcomes[names(outcomes)==1]
p = n_working / length(y_draws)

probs = c(1-p, p) # Probabilities of wether working
names(probs) = 100 * c((1-p), p) # Prob -> percentage

pdf("plots/2_3_pred_distr.pdf", width=grid_w, height=grid_h)

prob = hist(y_draws, breaks=2, plot=FALSE)
barplot(probs, 
        col=c("lightgray", "white"),
        legend=c("Not working", "Working"),
       ylab="Probability", 
       xlab="Working", 
       main="Predictive distribution for sample")

dev.off()