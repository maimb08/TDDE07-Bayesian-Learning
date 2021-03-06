
# -----------------------------
#  Lab 3 - Assignment 2 (b)
# -----------------------------


# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linköping University. http://mattiasvillani.com

##########    BEGIN USER INPUT #################
# Data options
data(faithful)
rawData <- faithful
x <- as.matrix(rawData['eruptions'])

# Lab 3 - 
data = read.table("data/rainfall.txt", header=FALSE)[,1]
x = as.matrix(data)

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(0,nComp) # Prior mean of theta
tau2Prior <- rep(10,nComp) # Prior std theta
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 100 # Number of Gibbs sampling draws

# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  thetaDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    thetaDraws[j] <- rgamma(1,param[j],1)
  }
  # Diving every column of ThetaDraws by the sum of the elements in that column.
  thetaDraws = thetaDraws/sum(thetaDraws) 
  return(thetaDraws)
}

# Simple function that converts between two different 
# representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
# nObs-by-nComp matrix with component allocations.
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) 
theta <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))

gibbs_thetas = matrix(0,nIter,2)
gibbs_sigmas = matrix(0,nIter,2)
for (k in 1:nIter){
  message(paste('Iteration number:',k))
  # Just a function that converts between different representations 
  # of the group allocations
  alloc <- S2alloc(S) 
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  w <- rDirichlet(alpha + nAlloc)
  
  # Update theta's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    theta[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  gibbs_thetas[k, ] = theta
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], 
       scale = (nu0[j]*sigma2_0[j] + 
                  sum((x[alloc == j] - theta[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  gibbs_sigmas[k,] = sigma2
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- w[j]*dnorm(x[i], mean = theta[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), 
         main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,theta[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + w[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topleft", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}


# Calculate mean of batches of 2 draws to visualize the
# auto correlation between sequential draws
t1 = c()
t2 = c()
s1 = c()
s2 = c()
for (i in 1:nIter){
  if(i%%2 == 0){
    t1 = c(t1, mean(gibbs_thetas[,1][i-1:i]))
    t2 = c(t2, mean(gibbs_thetas[,2][i-1:i]))
    s1 = c(s1, mean(gibbs_sigmas[,1][i-1:i]))
    s2 = c(s2, mean(gibbs_sigmas[,2][i-1:i]))
  }
}


# Plots displaying convergence of the Normal hyper 
# parameters during sampling

pdf("plots/3_1_2_gibbs_mixt.pdf")

par(mfrow=c(3,1))

# Plot comparison between kernel density, mixture of normals and a normal
# approximation
hist(x, breaks = 20, 
     cex=.1, 
     border="lightgray", 
     freq = FALSE, 
     xlim = c(xGridMin,xGridMax), 
     xlab="Precipitation", 
     ylab="Density", 
     main = "Rainfall: Mixture of Normals")

lines(xGrid, 
      mixDensMean, 
      type = "l", 
      lwd = 2, 
      lty = 4, 
      col = "black")

lines(xGrid, 
      dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), 
      type = "l", 
      lwd = 2, 
      col = "gray")

legend("topright", 
       box.lty = 1, 
       legend = c("Data histogram","Mixture density","Normal density"), 
       col=c("lightgray","black","gray"), 
       lwd = 2)

# Traceplot (Convergence)

# Theta conv.
y_min = min(gibbs_thetas)
y_max = max(gibbs_thetas)
traceplot(mcmc(gibbs_thetas[,1]), 
          ylim=c(y_min,y_max),
          main="Gibbs sampling: Theta convergence")
lines(mcmc(gibbs_thetas[,2]), col="gray")

legend("topright",
       box.lty = 1,
       legend = c(expression(paste(theta, " (1)", sep=" ")),
                  expression(paste(theta, " (2)", sep=" "))),
       col=c("black","gray"),
       lwd = 2)

# Sigma conv.
y_min = min(gibbs_sigmas)
y_max = max(gibbs_sigmas)
traceplot(mcmc(gibbs_sigmas[,1]), 
          ylim=c(y_min,y_max),
          main="Gibbs sampling: Sigma convergence")
lines(mcmc(gibbs_sigmas[,2]), col="gray")

legend("topright",
       box.lty = 1,
       legend = c(expression(paste(sigma, " (1)", sep=" ")),
                  expression(paste(sigma, " (2)", sep=" "))),
       col=c("black","gray"),
       lwd = 2)

# # Plot the auto correlation (convergence) between draws of mu
# min_t = min(c(min(t1), min(t2)))
# max_t = max(c(max(t1), max(t2)))
# plot(t1, 
#      type="l", 
#      ylim=c(min_t, max_t), 
#      cex=.1,
#      lwd=2,
#      main=expression(
#        paste("Convergence of Gibbs Sampling ", "(", theta, ")", sep=" ")
#        ),
#      xlab="Batches of sequential draws",
#      ylab=expression(paste("Mean of seq. draws of ", theta, sep=" ")))
# 
# lines(t2, lwd=2, col="gray")
# 
# legend("topright", 
#        box.lty = 1, 
#        legend = c(expression(paste(theta, " (1)", sep=" ")),
#                   expression(paste(theta, " (2)", sep=" "))), 
#        col=c("black","gray"), 
#        lwd = 2)

# Plot the auto correlation (convergence) between draws of sigma
# min_s = min(c(min(s1), min(s2)))
# max_s = max(c(max(s1), max(s2)))
# plot(s1, 
#      type="l", 
#      ylim=c(min_s, max_s), 
#      cex=.1,
#      lwd=2,
#      main=expression(
#        paste("Convergence of Gibbs Sampling ", "(", sigma^2, ")", sep=" ")
#        ),
#      xlab="Batches of sequential draws",
#      ylab=expression(paste("Mean of seq. draws of ", sigma^2, sep=" ")))
# 
# lines(s2, lwd=2, col="gray")
# 
# legend("topright", 
#        box.lty = 1, 
#        legend = c(expression(paste(sigma^2, " (1)", sep=" ")),
#                   expression(paste(sigma^2, " (2)", sep=" "))), 
#        col=c("black","gray"), 
#        lwd = 2)


dev.off()

#########################    Helper functions    #######################