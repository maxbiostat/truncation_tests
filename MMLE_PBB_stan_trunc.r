library(sumR)
library(extraDistr)
library(rstan)
source("Poisson_Betabinomial_trunc_aux.r")
model <- rstan::stan_model(file = "stan/Poisson_BetaBinomial_MLE.stan",
                           allow_undefined = TRUE)
get_p <- function(pars){
  pars[2]/(pars[2] + pars[3])
}
####################
Mu <- 1000
A <- 1
B <- 99
J <- 5000
simu <- simulate_obsdata(n = J, mu = Mu, a = A, b = B,
                         seed = 9561362)
data <- compress_counts(simu$obs_xprime)

mean(simu$obs_y)
mean(simu$obs_xprime)
##############
target.eps <- .Machine$double.eps
M <- 1E3

PBB.data <- list(
  K = data$K,
  x = data$count,
  n = data$freq,
  epsilon = target.eps,
  max_iter = M
)

Opt.MLE <- rstan::optimizing(obj = model,
                             data = PBB.data,
                             hessian = TRUE,
                             draws = 500,
                             verbose = TRUE)

point.ests <- as.numeric(Opt.MLE$par[1:3])
MLE.se <- as.numeric(sqrt(diag(solve(-Opt.MLE$hessian))))
conf.level <- 0.95
Z <- qnorm(p = (1 + conf.level)/2)

MLE.lower <- point.ests - Z*MLE.se
MLE.upper <- point.ests + Z*MLE.se

draws.mean <- colMeans(Opt.MLE$theta_tilde)[1:3]
draws.lwr <-  apply(Opt.MLE$theta_tilde, 2,
                    quantile, probs = (1 - conf.level)/2)[1:3]
draws.upr <- apply(Opt.MLE$theta_tilde, 2,
                   quantile, probs = (1 + conf.level)/2)[1:3]

draws.p <- apply(Opt.MLE$theta_tilde[, 1:3], 1,  get_p)

pl <- quantile(draws.p, probs = (1 - conf.level)/2)
pu <- quantile(draws.p, probs = (1 + conf.level)/2)  
  
MLE.ests <- tibble::tibble(
  parameter = c("mu", "alpha", "beta", "a/(a+b)"),
  true = c(Mu, A, B, A/(A+B)),
  point =  c(point.ests, get_p(point.ests)),
  # lwr = c(MLE.lower, get_p(MLE.lower)),
  # upr = c(MLE.upper, get_p(MLE.upper)),
  # point = c(draws.mean, get_p(draws.mean)),
  lwr = c(draws.lwr, pl),
  upr = c(draws.upr, pu)
)
MLE.ests

