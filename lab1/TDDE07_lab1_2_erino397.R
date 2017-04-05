require("geoR")

# 2.1

y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
mean = 3.5
n = length(y)
n_draws = 10000

# Sample variance
tausq = sum((log(y) - mean)^2) / n
# X ~ chi(n)
X_draws = rchisq(n_draws, n)
# This is a draw from inv_chi(n, tausq)
var = n * tausq / X_draws
interval = seq(min(var), max(var), 0.001)
invchisq = dinvchisq(interval, n, tau_sq)

# Plot draws against density
hist(var, 500, border="gray", prob=TRUE)
lines(interval, invchisq)


# 2.2

z = sqrt(var / 2)
G = 2*pnorm(z)-1
hist(G, 100)

# 2.3

get_mode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
mode = get_mode(G)

