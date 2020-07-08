# Test Code

family = poisson(link = "log")

n = 100
p = 10
p1 = 2

x = mvtnorm::rmvnorm(n = n,
                     mean = base::rep(x = 0,
                                      times = p),
                     sigma = base::diag(p))

beta = base::c(base::rep(x = 1,
                         times = p1),
               base::rep(x = 0,
                         times = p - p1))
eta = x %*% beta
mu = base::exp(eta)

y = stats::rpois(n = n, lambda = mu)


#####

# Example

data = base::data.frame(x,
                        y)
B = 10
m = n/2

RobStabExamp = test = model_space(data = data,
                                  B = B,
                                  m = m,
                                  family = family)
