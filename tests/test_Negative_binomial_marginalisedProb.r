library(sumR)
library(tidyverse)
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
    if (n[i] < p[3]) output[i] = -INFINITY;
    else
      output[i] = R::dnbinom_mu(n[i], p[1], p[0], true) +
        R::dbinom(p[3], n[i], p[2], true);
  }
  return output;
  }')

#######
Mu <- 1
Phi <- .1
Eta <- .01
obsX <- 0
Theta <- c(Mu, Phi, Eta, obsX)
lgL <- log(Mu) - matrixStats::logSumExp(c(log(Mu), log(Phi))) + log1p(-Eta)
Eps <- .Machine$double.eps
M <- 1E5
TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE)
#######
negativeBinomial_marginalised_lpmf_R(k = obsX, theta = Theta)
negativeBinomial_marginalised_lpmf_C(n = obsX, p = Theta)
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

approximation.list <- parallel::mclapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(loc = linha$Mu,
                disp = linha$Phi,
                detp = linha$Eta,
                obx = linha$X,
                epsilon = linha$eps)      
  }, mc.cores = 8)

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt$L_half <- ifelse(approximation.dt$L < 1/2, 1, 0)
approximation.dt$theo_bound <- approximation.dt$target_error * approximation.dt$L /(1-approximation.dt$L)
approximation.dt$theo_success <- abs(approximation.dt$error) <= approximation.dt$theo_bound 

approximation.dt$row_id <- paste(approximation.dt$mu,
                                 approximation.dt$phi,
                                 approximation.dt$eta,
                                 approximation.dt$obs_x,
                                 approximation.dt$target_error,
                                 sep = "_")

approximation.dt <- approximation.dt %>%
  group_by(row_id) %>%
  mutate(comparative_error = abs(error)/abs(error[method == 'Fixed_3e+05']))

approximation.dt <- approximation.dt %>%
  mutate(comparative_success = comparative_error <= 1)

approximation.dt <- approximation.dt %>%
  mutate(actual_success = as.logical(success + comparative_success))

aggregate(success~method, approximation.dt, mean)
aggregate(comparative_success~method, approximation.dt, mean)
aggregate(actual_success~method, approximation.dt, mean)

aggregate(success~method+target_error, approximation.dt, mean)
aggregate(comparative_success~method+target_error, approximation.dt, mean)
aggregate(actual_success~method+target_error, approximation.dt, mean)

aggregate(success~method+target_error+L_half, approximation.dt, mean)
aggregate(comparative_success~method+target_error+L_half, approximation.dt, mean)
aggregate(actual_success~method+target_error+L_half, approximation.dt, mean)

aggregate(success~method+L_half, approximation.dt, mean)
aggregate(comparative_success~method+L_half, approximation.dt, mean)
aggregate(actual_success~method+L_half, approximation.dt, mean)

aggregate(theo_success~method, approximation.dt, mean)

aggregate(theo_success~method+L_half, approximation.dt, mean)

subset(approximation.dt,
       !actual_success
       & method == "Threshold")%>% print(n = 20)

subset(approximation.dt,
       !actual_success 
       & method == "BoundingPair") %>% print(n = 10)
