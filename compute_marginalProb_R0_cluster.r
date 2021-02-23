get_marginalised_prob <- function(r, w, p, log = FALSE, verbose = FALSE){
  ## Get Pr(S = 0 | r, w, p)
  a <- r/w
  b <- r/(w * (1 + r/w)^(w + 1)) * (1-p)
  solve_for_u <- function(b, w){
    ## We know b < 1 because a and omega are positive
    ## it suffices to find the positive root smaller than 1/omega, which
    ## is the positive root of the objective function obj_fun
    b_of_u <- function(u) u/(1 + u)^{w + 1}
    obj_fun <- function(cand) (b-b_of_u(cand))^2
    u <- optimise(obj_fun, c(0, 1/w), tol = 1e-25)$minimum
    return(u)
  }
  ustar <- solve_for_u(b, w = w)
  if(verbose) cat("u* is:", ustar, "\n")
  ans <- log1p(a)-log(a) + log(ustar) - log1p(ustar)
  if(!log) ans <- exp(ans)
  return(ans)
}
dR0 <- function(y, r, w, log = FALSE){
  if(y == 0){
    dens <- -Inf
  }else{
    c1 <- lgamma(w*y + y-1)
    c2 <- lgamma(w*y)
    c3 <- lgamma( y+1)
    c4 <- (y-1) * (log(r) - log(w))
    c5 <- (w*y + y -1) * log(1 + r/w)
    dens <- c1 - (c2+c3) + (c4-c5)
  }
  if(!log){
    dens <- exp(dens)
  }
  return(dens)
}
dR0 <- Vectorize(dR0)
################
source("truncation_package.r")
source("truncation_aux.r")

R0_cluster_marginal_detection <- function(k, theta){
  R0 <- theta[1]
  omega <- theta[2]
  x <- theta[3]
  p <- theta[4]
  ans <- ifelse(k < x, -Inf,
                dR0(y = k, r = R0, w = omega, log = TRUE) + dbinom(x = x, size = k, prob = p, log = TRUE))
  return(ans)
}

R0 <- .82
Omega <- .05
detProb <- .1
obsX <- 0
Theta <- c(R0, Omega, obsX, detProb)
Eps <- 1E-20
lgL <- log(R0) + log1p(-detProb) + (1 + Omega) * ( log1p(Omega) - log(R0 + Omega) )

N <- 1000
lps <- R0_cluster_marginal_detection(obsX:N, theta = Theta)
plot(diff(lps),
      main= "Ratios *increase* to L",
     ylab = expression(log(a[n+1]/a[n])), xlab = expression(n))
abline(h = lgL, lwd = 2, lty = 2)

(TrueValue <- get_marginalised_prob(r = R0, w = Omega, p = detProb, log = TRUE) )
log_sum_exp(lps)
robust_difference(TrueValue, log_sum_exp(lps))

result <- compare_approximations(R0_cluster_marginal_detection, theta = Theta,
                       exact = TrueValue, epsilon = Eps, max_iter = 1E5, logL = lgL)
result

library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(R0_cluster_marginal_detection, theta = Theta, epsilon = Eps, n0 = obsX, max_iter = mit),
  naive_threshold = approx_naive_under(R0_cluster_marginal_detection, theta = Theta, epsilon = Eps, n0 = obsX, max_iter = mit),
  doubling = approx_doubling(R0_cluster_marginal_detection, theta = Theta, epsilon = Eps, n0 = obsX, max_iter = mit),
  adaptive = approx_adaptive(R0_cluster_marginal_detection, theta = Theta, epsilon = Eps, n0 = obsX, max_iter = mit)
)

result
timing

