library(gamlss.dist)
library(rstan)
source("double_poisson_binomial_aux.r")

Mu <- 40
Phi <- 1
PP <- .5
J <- 1000
simu <- simulate_obsdata(n = J, mu = Mu, phi = Phi, p = PP)
data <- compress_counts(simu$obs_x)


stanData <- list(
  K = data$K,
  n = data$freq,
  x = data$count,
  x_int = data$count,
  s_mu = .01,
  r_mu = .01,
  nu_sd = 1,
  a_p = 1,
  b_p = 1,
  eps = .Machine$double.eps,
  M = 1e5,
  bayesian = 0
)

##########
use_Rstan <- FALSE

if(use_Rstan){
  model <- rstan::stan_model(file = "stan/double_Poisson_adaptive.stan",
                             allow_undefined = TRUE)
  
  
  Opt.MLE <- rstan::optimizing(obj = model,
                               data = stanData,
                               hessian = TRUE,
                               verbose = TRUE)
  
  
  FisherInfo <- solve(Opt.MLE$hessian)
  MLE.sds <- sqrt(diag(FisherInfo))
  
  conf.level <- 0.95
  Z <- qnorm(p = (1 + conf.level)/2)
  MLE.lower <- Opt.MLE$par[1:3] - Z*MLE.sds
  MLE.upper <- Opt.MLE$par[1:3] + Z*MLE.sds
  
  MLE.ests <- tibble::tibble(
    parameter = names(Opt.MLE$par)[1:3],
    mean =  Opt.MLE$par[1:3],
    lwr = MLE.lower,
    upr = MLE.upper
  )
  MLE.ests
}

########

library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

adaptive_impl <- cmdstanr::cmdstan_model("stan/double_Poisson_adaptive.stan",
                                         include_paths = "./stan/")


opt_adaptive <- adaptive_impl$optimize(data = stanData)
opt_adaptive$mle()

mcmc_iter <- 500
stanData$bayesian <- 1
adaptive.raw <-
  adaptive_impl$sample(data = stanData,
                   refresh = floor(mcmc_iter/5),
                   chains = 4,
                   parallel_chains = 4,
                   iter_warmup = mcmc_iter,
                   adapt_delta = .90,
                   max_treedepth = 12,
                   iter_sampling = mcmc_iter,
                   show_messages = FALSE)
adaptive.mcmc <- stanfit(adaptive.raw)

adaptive.mcmc
check_hmc_diagnostics(adaptive.mcmc)
adaptive.raw$time()
