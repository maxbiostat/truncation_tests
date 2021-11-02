library(rstan)
library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
source("R0_Betabinomial_truc_aux.r")
####################
R0 <- .67
Disp <- .5
A <- 1
B <- 10
J <- 1000
simu <- simulate_obsdata(n = J,
                         r = R0, w = Disp,
                         a = A, b = B)
data <- compress_counts(simu$obs_xprime)
data$K
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(R0, Disp, A, B),
                       verbose = TRUE)

cluster.data <- list(D = data$K,
                     s = data$count,
                     n = data$freq,
                     epsilon = 1e6 * .Machine$double.eps,
                     max_iter = 1E5
)

R0_model <- cmdstan_model("stan/R0_clusters_BetaBinomial_MLE.stan",
                          include_paths = "./stan/")

system.time(
  opt.MAP <-  R0_model$optimize(data = cluster.data,
                                algorithm = "lbfgs",
                                tol_obj = 1E-20,
                                tol_rel_obj = 50000,
                                tol_grad = 1e-12,
                                tol_param = 1e-13
                                )
)
opt.MAP$mle()[1:4]
c(R0, Disp, A, B)

use_rstan <- FALSE

if(use_rstan){
  
  model <- rstan::stan_model(file = "stan/R0_clusters_BetaBinomial_MLE.stan",
                             allow_undefined = TRUE)
  
  
  Opt.MLE <- rstan::optimizing(obj = model,
                               data = cluster.data,
                               hessian = TRUE,
                               verbose = TRUE)
  
  point.ests <- as.numeric(Opt.MLE$par[1:4])
  
  FisherInfo <- solve(-Opt.MLE$hessian)
  MLE.sds <- sqrt(diag(FisherInfo))
  
  conf.level <- 0.95
  Z <- qnorm(p = (1 + conf.level)/2)
  MLE.lower <- point.ests - Z*MLE.sds/sqrt(compressed$J)
  MLE.upper <- point.ests + Z*MLE.sds/sqrt(compressed$J)
  
  MLE.ests <- tibble::tibble(
    parameter = c("R0", "omega", "psi"),
    true = c(true.R0, true.omega, true.psi),
    point = c(point.ests[1], guessDispersion, point.ests[2]),
    lwr = c(MLE.lower[1], NA, MLE.lower[2]),
    upr = c(MLE.upper[1], NA, MLE.upper[2])
  )
  MLE.ests
  
  Opt.MLE$value
}

