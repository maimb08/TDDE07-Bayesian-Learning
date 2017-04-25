require(MASS)
require(geoR)
# Lab 2 - Assignment 1

df = read.table("data/TempLinkoping2016.txt", header=TRUE)

n = 366
n_params = 3

# (a)

# Linear model (test)
fit = lm(temp ~ time + I(time^2), data=df)
summary(fit)


mu0 = c(-10, 90, -80)
covar0 = diag(c(.5, .5, .5))
v0 = n - n_params
sigma_sq0 = 1

# (b)

t = df[["time"]]
temp = df[["temp"]]

plot(t, 
     temp, 
     type="p", 
     col="lightgray", 
     ylim=c(-30, 30), 
     xlim=c(0,1), 
     xlab="Fraction of year", 
     ylab="Temperature",
     main="Prior prediction")

n_draws = 10
for (iter in 1:n_draws) {
  sigma_sq = rinvchisq(n=1, df=v0, scale=sigma_sq0)
  beta = mvrnorm(n=1, mu=mu0, Sigma=sigma_sq*ginv(covar0))
  error = rnorm(n=1, 0, 1)
  preds = beta[1] + beta[2]*t + beta[3]*I(t)^2 + error
  lines(t, preds)
}

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

plot(t, 
     temp, 
     type="p", 
     col="lightgray", 
     ylim=c(-30, 30), 
     xlim=c(0,1), 
     xlab="Fraction of year", 
     ylab="Temperature",
     main="Posterior prediction")

beta_1s = c()
beta_2s = c()
beta_3s = c()

n_draws = 100
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
preds = beta_1 + beta_2*t + beta_3*I(t)^2 + error
lines(t, preds, lwd=2)

# Equal tail interval
beta_1_eti = quantile(beta_1s, probs=c(0.025, 0.975))
beta_2_eti = quantile(beta_2s, probs=c(0.025, 0.975))
beta_3_eti = quantile(beta_3s, probs=c(0.025, 0.975))

# Lower
lower_beta_1 = beta_1_eti[1]
lower_beta_2 = beta_2_eti[1]
lower_beta_3 = beta_3_eti[1]
preds = lower_beta_1 + lower_beta_2*t + lower_beta_3*I(t)^2 + error
lines(t, preds, col="blue", lwd=2)

# Higher
higher_beta_1 = beta_1_eti[2]
higher_beta_2 = beta_2_eti[2]
higher_beta_3 = beta_3_eti[2]
preds = higher_beta_1 + higher_beta_2*t + higher_beta_3*I(t)^2 + error
lines(t, preds, col="red", lwd=2)

legend("bottomright", 
       legend = c("Mean","5%", "95%"),
       fill = c("black", "blue", "red"),
       inset = 0.02)

# (d)



# (e)

# Chooses new mu0 and covar0 as the previous posterior hyper parameters
mu0 = mu_n
covar0 = covar_n