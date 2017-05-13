require(MASS)
require(geoR)

grid_w = 5
grid_h = 5

# -------
#  Lab 4
# -------

data = read.table("data/eBayNumberOfBidderData.dat", header=TRUE)


n = length(data)
n_features = ncol(data) - 2 # Except y and const

feature_labels = colnames(data[,2:ncol(data)])
  
y = data$nBids
X = as.matrix(data[,2:ncol(data)])

X_X = t(X)%*%X

# -----
#  (a)
# -----

glm_model = glm(nBids ~ 0 + ., data = data, family = poisson)

pdf("./plots/4_1_1_mle_beta.pdf", height=grid_h)

par(mfrow=c(2,4), oma = c(0, 0, 3, 0))
for (i in 2:ncol(X)){
  mean = glm_model$coefficients[i]
  std_dev = summary(glm_model)[["coefficients"]][,2][i]
  x_grid = seq(mean-4*std_dev, mean+4*std_dev, 0.001)
  plot(x_grid,
       dnorm(x_grid, mean=mean, sd=std_dev),
       type="l",
       lwd=2,
       ylab="Density",
       xlab=expression(beta),
       main=feature_labels[i])
}
title("Normal approximation of beta", outer=TRUE, cex=1.5)

dev.off()

# -----
#  (b)
# -----

# Beta prior parameters (Zellner’s g-prior)
mu0 = 0
covar0 = 100 * ginv(X_X)

# This is the log posterior density of the beta(s+a,f+b) density
# LogPostBernBeta <- function(theta, s, f, a, b){
#   logPost <- (s+a-1)*log(theta) + (f+b-1)*log(1-theta)
#   return(logPost)
# }
