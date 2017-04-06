# 3.1

y = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mean = 2.39
k = seq(0.001,10,0.01)

calc_prob = function(k_val, y_val){
  prob = exp(k_val * cos(y_val - mean)) / (2*pi*besselI(k_val, 0))
  
  return (prob)
}

calc_post = function(k_val){
  prob = sapply(y, calc_prob, k_val=k_val)
  post = prod(prob) * dexp(k_val)
  
  return (post)
}

posterior = sapply(k, calc_post)

plot(k, posterior, type="l")

# 3.2

#get_mode = function(v) {
#  uniqv = unique(v)
#  uniqv[which.max(tabulate(match(v, uniqv)))]
#}
#mode = get_mode(posterior)

mode = max(posterior)
