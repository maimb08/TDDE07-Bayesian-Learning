require(MASS)
require(geoR)

grid_w = 5
grid_h = 4

# -------
#  Lab 4
# -------

data = read.table("data/eBayNumberOfBidderData.dat", header=TRUE)


n = length(data)
n_features = ncol(data) - 2 # Except y and const

feature_labels = colnames(data[,2:n_features])
  
y = data$nBids
X = as.matrix(data[,3:ncol(data)])

X_X = t(X)%*%X

# -----
#  (a)
# -----

glm_model = glm(nBids ~ 0 + ., data = data, family = poisson)

# -----
#  (b)
# -----

# Beta prior parameters (Zellner’s g-prior)
mu0 = 0
covar0 = 100 * ginv(X_X)

# This is the log posterior density of the beta(s+a,f+b) density
LogPostPoisson = function(theta, s, f, a, b){
  log_prior = log(dnorm(mu0, covar0))
  log_likelihood = 
  return(logPost)
}

# This is the log posterior density of the beta(s+a,f+b) density
# LogPostBernBeta <- function(theta, s, f, a, b){
#   logPost <- (s+a-1)*log(theta) + (f+b-1)*log(1-theta)
#   return(logPost)
# }
