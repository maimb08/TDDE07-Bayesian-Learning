require(MASS)
require(geoR)


grid_w = 5
grid_h = 4


# Lab 2 - Assignment 1

df = read.table("data/TempLinkoping2016.txt", header=TRUE)

n = 366
n_params = 3

# (a)

# Linear model (test)
fit = lm(temp ~ time + I(time^2), data=df)
summary(fit)


# Set low v0 as we are unsure about the variance of the beta parameters
mu0 = c(-10, 90, -80)
covar0 = diag(c(1, 1, 1))
v0 = 10
sigma_sq0 = 100

# (b)

t = df[["time"]]
temp = df[["temp"]]

pdf("plots/1_2_prior_draws.pdf", width=grid_w, height=grid_h)

plot(t, 
     temp,
     type="p", 
     col="lightgray", 
     ylim=c(-30, 30), 
     xlim=c(0,1), 
     xlab="Fraction of year", 
     ylab="Temperature",
     main="Prior prediction",
     cex.main=.9, 
     cex.lab=.9, 
     cex.axis=.8)

n_draws = 10
for (iter in 1:n_draws) {
  sigma_sq = rinvchisq(n=1, df=v0, scale=sigma_sq0)
  beta = mvrnorm(n=1, mu=mu0, Sigma=sigma_sq*ginv(covar0))
  error = rnorm(n=1, 0, 1)
  preds = beta[1] + beta[2]*t + beta[3]*I(t)^2 + error
  lines(t, preds)
}

dev.off()

# (c)

# 366x3 (1 t t^2)
X = t(rbind(rep(1,n), t, I(t)^2))
y = as.vector(temp)

X_X = t(X) %*% X

# Least squares approximate of beta
beta_hat = ginv(X_X) %*% t(X) %*% y

# Posterior parameters
mu_n = ginv(X_X + covar0)%*%(X_X%*%beta_hat+covar0%*%mu0)
covar_n = X_X + covar0
v_n = v0 + n
sigma_sq_n = as.double((1/v_n)*(v0%*%sigma_sq0 + (t(y)%*%y + t(mu0)%*%covar0%*%mu0 - t(mu_n)%*%covar_n%*%mu_n)))

# pdf("plots/1_3_posterior_draws.pdf", width=grid_w, height=grid_h)

plot(t, 
     temp, 
     type="p", 
     col="lightgray", 
     ylim=c(-30, 30), 
     xlim=c(0,1), 
     xlab="Fraction of year", 
     ylab="Temperature",
     main="Posterior prediction with credibility interval",
     cex.main=.9, 
     cex.lab=.9, 
     cex.axis=.8)

beta_1s = c()
beta_2s = c()
beta_3s = c()

n_draws = 5000
for (iter in 1:n_draws) {
  sigma_sq = rinvchisq(n=1, df=v_n, scale=sigma_sq_n)
  beta = mvrnorm(n=1, mu=mu_n, Sigma=sigma_sq*ginv(covar_n))
  error = rnorm(n=1, 0, sigma_sq)
  # Add beta draws
  beta_1s = c(beta_1s, beta[1])
  beta_2s = c(beta_2s, beta[2])
  beta_3s = c(beta_3s, beta[3])
}

error = 0

# Line using mean values from simulation
beta_1 = mean(beta_1s)
beta_2 = mean(beta_2s)
beta_3 = mean(beta_3s)
preds_mean = beta_1 + beta_2*t + beta_3*I(t)^2 + error
lines(t, preds_mean, lwd=2)

# Equal tail interval
beta_1_eti = quantile(beta_1s, probs=c(0.025, 0.975))
beta_2_eti = quantile(beta_2s, probs=c(0.025, 0.975))
beta_3_eti = quantile(beta_3s, probs=c(0.025, 0.975))

# Lower
lower_beta_1 = beta_1_eti[1]
lower_beta_2 = beta_2_eti[1]
lower_beta_3 = beta_3_eti[1]
preds_low = lower_beta_1 + lower_beta_2*t + lower_beta_3*I(t)^2 + error
lines(t, preds_low, col="blue", lwd=2)

# Higher
higher_beta_1 = beta_1_eti[2]
higher_beta_2 = beta_2_eti[2]
higher_beta_3 = beta_3_eti[2]
preds_high = higher_beta_1 + higher_beta_2*t + higher_beta_3*I(t)^2 + error
lines(t, preds_high, col="red", lwd=2)

legend("bottomright", 
       legend = c("Mean","Lower", "Upper"),
       fill = c("black", "blue", "red"),
       inset = 0.02)

# dev.off()

# (d)

# pdf("plots/1_4_warm_days.pdf", width=grid_w, height=grid_h)

# Could also solve for time by derivation set to zero
i = which.max(preds_mean)
max_time = t[i]
max_day = round(366 * max_time)

# Std dev chosen by prediction plot
xGrid = seq(0, 366, 1)
n = dnorm(xGrid, mean=max_time*366, sd=0.1*366)
plot(xGrid, 
     n, 
     type="l", 
     lwd=2, 
     xlab="Days", 
     ylab="Probability", 
     main="Prob. distribution of day with hottest temperature",
     cex.main=.9, 
     cex.lab=.9, 
     cex.axis=.8)

# dev.off()

# (e)

# Chooses new mu0 and covar0 as the previous posterior hyper parameters and the additional
# betas set to zero (to combat overfitting).
mu0 = c(mu_n, 0, 0, 0, 0)
covar0 = matrix(0, nrow=8, ncol=8)
covar0[1:3, 1:3] = covar_n

# Might need to set variance of other parameters to smaall value to 
# allow for a posterior value for those parameters other than zero
