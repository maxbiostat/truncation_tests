library(sumR)
library(bbmle)
source("Erlang_aux.r")
source("aux/aux.r")
###########

Mu <- 1500
B <- .56
J <- 5
data <- simulate_obsdata(n = J, mu = Mu, b = B,
                         seed = 1234)
pars <- c(Mu, B)

ll.bessel.precomp <- marginal_loglikelihood(data = data, pars = c(Mu, B),
                                            eps = .Machine$double.eps,
                                            verbose = TRUE)

ll.bessel.fast <- marginal_loglikelihood_fast(data = data, pars = c(Mu, B),
                                              eps = .Machine$double.eps,
                                              verbose = TRUE)

ll.full <- marginal_loglikelihood_full(data = data, pars = c(Mu, B),
                                       eps = .Machine$double.eps,
                                       verbose = TRUE)

ll.full.fast <- marginal_loglikelihood_full_fast(data = data, pars = c(Mu, B),
                                                 eps = .Machine$double.eps,
                                                 verbose = TRUE)

robust_difference(ll.bessel.fast, ll.bessel.precomp)
robust_difference(ll.bessel.precomp, ll.full)
robust_difference(ll.full, ll.full.fast)

#####################
library(numDeriv)

bessel_mu <- function(x){
  marginal_loglikelihood(data = data, pars = c(x, B),
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
bessel_mu <- Vectorize(bessel_mu)

bessel_beta <- function(x){
  marginal_loglikelihood(data = data, pars = c(Mu, x),
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
bessel_beta <- Vectorize(bessel_beta)

full_mu <- function(x){
  marginal_loglikelihood_full_fast(data = data, pars = c(x, B),
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
full_mu <- Vectorize(full_mu)

full_beta <- function(x){
  marginal_loglikelihood_full_fast(data = data, pars = c(Mu, x),
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
full_beta <- Vectorize(full_beta)

Mus <- seq(Mu-Mu/2, Mu+Mu/2, length.out = 5)

plot(Mus, bessel_mu(Mus))
lines(Mus, full_mu(Mus))
abline(v = Mu, lwd = 2)

Betas <- seq(B-B/2, B+B/2, length.out = 5)

plot(Betas, bessel_beta(Betas))
lines(Betas, full_beta(Betas))
abline(v = B, lwd = 2)

sapply(1:5, function(i){
  robust_difference(bessel_mu(Mus[i]), full_mu(Mus[i]))
})
sapply(1:5, function(i){
  robust_difference(bessel_beta(Betas[i]), full_beta(Betas[i]))
})

bessel_pars <- function(x){
  marginal_loglikelihood(data = data, pars = x,
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
bessel_pars <- Vectorize(bessel_pars)

full_pars <- function(x){
  marginal_loglikelihood_full_fast(data = data, pars = x,
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}
full_pars <- Vectorize(full_pars)

sapply(Mus, function(m){
  numDeriv::hessian(
    func = bessel_mu,
    x = m
  )  
})

sapply(Mus, function(m){
  numDeriv::hessian(
    func = full_mu,
    x = m
  )  
})

num.details <- list(eps = 1E-10, d = 1E-5, r = 6,
                    zero.tol = .Machine$double.eps)

sapply(Betas, function(b){
  numDeriv::grad(
    func = bessel_beta,
    x = b,
    method.args = num.details
  )  
})

sapply(Betas, function(b){
  numDeriv::grad(
    func = full_beta,
    x = b,
    method.args = num.details
  )  
})

sapply(Betas, function(b){
  numDeriv::hessian(
    func = bessel_beta,
    x = b,
    method.args = num.details
  )  
})

sapply(Betas, function(b){
  numDeriv::hessian(
    func = full_beta,
    x = b,
    method.args = num.details
  )  
})

####################

nll_adaptive_full_fast <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full_fast(data,
                                    pars = theta,
                                    eps =.Machine$double.eps,
                                    verbose = FALSE)
}

nll_adaptive <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood(data,
                          pars = theta,
                          eps =.Machine$double.eps,
                          verbose = FALSE)
}
get_estimates <- function(data, minusloglik_fun, prof_ci = FALSE){
  
  inits <- log(
    c(runif(1, 1, 100),
      runif(1, .1, 50))
  )
  running.time <- system.time(
    MLE <- bbmle::mle2(minuslogl = minusloglik_fun,
                       start = list(log_mu = inits[1],
                                    log_beta = inits[2]),
                       method = "L-BFGS-B"
    )
  )
  
  summy <- summary(MLE)
  
  point.ests <- as.numeric(summy@coef[, 1])
  MLE.se <- as.numeric(summy@coef[, 2])
  conf.level <- 0.95
  Z <- qnorm(p = (1 + conf.level)/2)
  
  MLE.lower <- point.ests - Z*MLE.se
  MLE.upper <- point.ests + Z*MLE.se
  
  MLE.ests <- tibble::tibble(
    parameter = c("mu", "beta"),
    true = c(Mu, B),
    point =  exp(point.ests),
    lwr = exp(MLE.lower),
    upr = exp(MLE.upper),
    time = running.time[3]
  )
  out <- list(
    obj = MLE,
    results = MLE.ests
  )
  if(prof_ci){
    out$CI <- confint(MLE)
  }
  return(out)
}

est.adaptive.full.fast <- get_estimates(data = data,
                                        nll_adaptive_full_fast, 
                                        prof_ci = TRUE)

est.adaptive.full.fast$results
exp(est.adaptive.full.fast$CI)


est.adaptive.bessel.fast <- get_estimates(data = data,
                                          nll_adaptive,
                                          prof_ci = TRUE)

est.adaptive.bessel.fast$results
exp(est.adaptive.bessel.fast$CI)

source("Erlang_bootstrap.r")

Nboot <- 50
full.boot.nonpar <- erlang_bootstrap(dt = data,
                                     Nrep = Nboot,
                                     ncores = 8,
                                     mll = nll_adaptive_full_fast,
                                     mle = est.adaptive.full.fast$results$point,
                                     type = "nonparametric")
full.boot.par <- erlang_bootstrap(dt = data,
                                     Nrep = Nboot,
                                     ncores = 8,
                                     mll = nll_adaptive_full_fast,
                                     mle = est.adaptive.full.fast$results$point,
                                     type = "parametric")

