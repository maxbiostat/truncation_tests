library(sumR)
source("../aux/aux.r")
#######
negativeBinomial_marginalised_lpmf_R <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  x <- theta[4]
  ans <- ifelse(k < x,  -Inf, 
                dnbinom(k, mu = mu, size = phi, log = TRUE) +
                  dbinom(x = x, size = k, prob = eta, log = TRUE))
  return(ans)
}

Rcpp::cppFunction(code='
  NumericVector negativeBinomial_marginalised_lpmf_C(IntegerVector n, NumericVector p)
  {
  NumericVector output(n.size());
  for (int i = 0; i < n.size(); i++) {
    if (p[4] > n[i]) output[i] = -INFINITY;
    else
      output[i] = R::dnbinom_mu(n[i], p[1], p[0], true) +
        R::dbinom(p[3], n[i], p[2], true);
  }
  return output;
  }')

#######
Mu <- 100
Phi <- .05
Eta <- .025
obsX <- 20
Theta <- c(Mu, Phi, Eta, obsX)
lgL <- log(Mu) - matrixStats::logSumExp(c(log(Mu), log(Phi))) + log1p(-Eta)
Eps <- .Machine$double.eps
M <- 1E5
TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE)
#######
negativeBinomial_marginalised_lpmf_R(k = obsX + 2, theta = Theta)
negativeBinomial_marginalised_lpmf_C(n = obsX + 2, p = Theta)
#######

result.R <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_R,
                                 theta = Theta,
                                 exact = TrueValue,
                                 eps = Eps,
                                 logL = lgL)

result.C <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)

result.R
result.C

( error.bound <- Eps * exp(lgL - log1m_exp(lgL)) ) ## This should be theoretical bound

exp(lgL)

####
##$ Table
spit_approx <- function(loc, disp, detp, obx, epsilon){
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = negativeBinomial_marginalised_lpmf_C,
      theta =  c(loc, disp, detp, obx),
      logL = log(loc) - matrixStats::logSumExp(c(log(loc), log(disp))) + log1p(-detp),
      exact = dnbinom(x = obx, mu = detp*loc, size = disp, log = TRUE),
      max_iter = 3E5,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    mu = loc,
    phi = disp,
    eta = detp,
    obs_x = obx,
    L = exp(log(loc) - matrixStats::logSumExp(c(log(loc), log(disp))) + log1p(-detp)),
    target_error = ans$target_error,
    error = ans$error,
    success = abs(ans$error) <= ans$target_error,
    n_iter = ans$n_evaluations,
    method = ans$Method
  )
  return(out)
}

approximation.grid <- expand.grid(
  Mu = c(1, 10, 100),
  Phi = c(.1, .5, 1, 10),
  Eta = c(.01, .1, .5, .75),
  X = c(0, 5, 10),
  eps = 10^(0:6)*.Machine$double.eps
)

approximation.list <- lapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(loc = linha$Mu,
                disp = linha$Phi,
                detp = linha$Eta,
                obx = linha$X,
                epsilon = linha$eps)      
  })

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt$L_half <- ifelse(approximation.dt$L < 1/2, 1, 0)
approximation.dt$theo_bound <- approximation.dt$target_error * approximation.dt$L /(1-approximation.dt$L)
approximation.dt$theo_success <- abs(approximation.dt$error) <= approximation.dt$theo_bound 

aggregate(success~method, approximation.dt, mean)
aggregate(success~method+target_error, approximation.dt, mean)
aggregate(success~method+target_error+L_half, approximation.dt, mean)
aggregate(success~method+L_half, approximation.dt, mean)

aggregate(theo_success~method, approximation.dt, mean)
aggregate(theo_success~method+L_half, approximation.dt, mean)

#### Investigating the instances in which adaptive/naive fails
#### but the sum with a huge (3E5) number of terms does not.

approximation.dt$row_id <- paste(approximation.dt$mu,
                                 approximation.dt$phi,
                                 approximation.dt$eta,
                                 approximation.dt$obs_x,
                                 approximation.dt$target_error,
                                 sep = "_")
( fail.adapt <- subset(approximation.dt, method == "Adaptive" & !success) )
( fail.naive <- subset(approximation.dt, method == "Naive" & !success) )
( fail.huge <- subset(approximation.dt, method == "Fixed_3e+05" & !success) )

interesting.naive <- subset(approximation.dt,
       row_id %in% setdiff(fail.naive$row_id, fail.huge$row_id))

interesting.adaptive <- subset(approximation.dt,
       row_id %in% setdiff(fail.adapt$row_id, fail.huge$row_id))

interesting.adaptive
tail(interesting.adaptive)
