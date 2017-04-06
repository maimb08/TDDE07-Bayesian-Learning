# 3.1

y = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mean = 2.39
k = seq(0.001,10,0.01)

calc_prob = function(k_val, y_val){
  prob = exp(k_val * cos(y_val - mean)) / (2*pi*besselI(k_val, 0))
  
  return (prob)
}

calc_post = function(k_val){
  probabilities = sapply(y, calc_prob, k_val=k_val)
  prior = dexp(k_val)
  posterior = prod(probabilities) * prior
  
  return (posterior)
}

posterior = sapply(k, calc_post)

plot(k, posterior, type="l", col="black")

# 3.2

# Find the k value which maximizes posterior
mode = k[which.max(posterior)]

# Plot vertical line where k maximizes posterior
abline(v=mode, col="gray")

# Legend for posterior and mode
legend('topright', c('Posterior', 'Mode'), fill=c("black", "gray"), inset=0.02)
