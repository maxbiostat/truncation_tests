map_Rn <- function(x){
  return(
    c(
      arm::logit(x[1]),
      log(x[2]),
      arm::logit(x[3])
    )
  )
}
map_back <- function(y){
  return(
    c(
      arm::invlogit(y[1]),
      exp(y[2]),
      arm::invlogit(y[3])
    )
  )
}
is_in <- function(x, l, u){
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}
organise_output <- function(fit, truePars, name, conf.level = 0.95){
  summy <- summary(fit)
  point.ests <- as.numeric(coef(fit))
  MLE.se <- as.numeric(summy@coef[, 2])
  if(length(MLE.se) == 2) MLE.se <- c(MLE.se[1], 0, MLE.se[2])
  Z <- qnorm(p = (1 + conf.level)/2)
  MLE.lower <- point.ests - Z*MLE.se
  MLE.upper <- point.ests + Z*MLE.se
  
  CENTER <- map_back(point.ests)
  LWR <- map_back(MLE.lower)
  UPR <- map_back(MLE.upper)
  
  MLE.ests <- tibble::tibble(
    parameter = c("R_0", "omega", "nu"),
    point = CENTER,
    lwr = LWR,
    upr = UPR,
    true = c(truePars[1], truePars[2], truePars[3]),
    covers = is_in(truePars, l = LWR, u = UPR),
    bias = truePars-CENTER,
    method = name
  )
  return(MLE.ests)
}

##
library(sumR)
library(extraDistr)
library(stats4)
source("R0_sentinel_trunc_aux.r")
####################
R0 <- .95
Disp <- .5
obsPr <- .8
J <- 20000
simu <- simulate_sentinel_data(n = J,
                               r = R0, w = Disp,
                               nu = obsPr,
                               seed = 666)
data <- compress_counts(simu$obs_xprime)
dd <- compress_counts(simu$obs_x)

mean(simu$obs_y)
##############
## Testing
pars <- c(R0, Disp, obsPr)
## Ratio, L
logLL <- log(pars[1]) +
  (1 + pars[2]) * (log1p(pars[2]) - log(pars[1] + pars[2]) ) + log1p(-pars[3])

lps <- sentinel_summand(k = 0:1E5, theta = pars)
tail(diff(lps))
logLL

## Pr(X = 0)
lp0 <- sumR::infiniteSum(logFunction = sentinel_summand,
                         parameters = c(pars, 0),
                         logL = logLL,
                         epsilon = .Machine$double.eps,
                         n0 = 1,
                         maxIter = 1E5)
lp0$n
exp(lp0$sum)
get_marginalised_prob(r = R0, w = Disp, p = obsPr)
empirical.probs.x <- dd$freq/sum(dd$freq)
empirical.probs.x[1]

## Pr(X'= k)

lpp0 <- log1p(-exp(lp0$sum))
raw.lpmfs <- sapply(data$count, function(k) sentinel_lpmf(k , theta = pars))
theo.probs <- exp(raw.lpmfs - lpp0)

probs.xprime <- tibble::tibble(
  x_prime = data$count,
  empirical_pr = data$freq/sum(data$freq),
  theo_pr = theo.probs,
  unn_theo = exp(raw.lpmfs)
) 
probs.xprime

##############
nll_adaptive <- function(logit_R0, log_omega, logit_nu){
  theta <- map_back(c(logit_R0, log_omega, logit_nu))
  -marginal_loglikelihood_adaptive(data, pars = theta)
}

nll_exact <- function(logit_R0, log_omega, logit_nu){
  theta <- map_back(c(logit_R0, log_omega, logit_nu))
  -marginal_loglikelihood_exact(data, pars = theta)
}


lower_l <- map_Rn(c(.001, .01, .0001))
upper_l <- map_Rn(c(.9999, 20, .9999))
inits <- map_Rn(c(.5, .5, .3))

names(lower_l) <- names(upper_l) <- names(inits) <- c("logit_R0", "log_omega",
                                                      "logit_nu")

MLE.adaptive <- bbmle::mle2(minuslogl = nll_adaptive,
                            start = list(logit_R0 = inits[1],
                                         log_omega = inits[2],
                                         logit_nu = inits[3]),
                            lower = lower_l,
                            upper = upper_l,
                            method = "L-BFGS-B"
)

MLE.exact <- bbmle::mle2(minuslogl = nll_exact,
                         start = list(logit_R0 = inits[1],
                                      log_omega = inits[2],
                                      logit_nu = inits[3]),
                         lower = lower_l,
                         upper = upper_l,
                         method = "L-BFGS-B"
)


MLE.adaptive.fixed <- bbmle::mle2(minuslogl = nll_adaptive,
                            start = list(logit_R0 = inits[1],
                                         logit_nu = inits[3]),
                            fixed = list(log_omega = log(Disp)),
                            lower = lower_l,
                            upper = upper_l,
                            method = "L-BFGS-B"
)

res.adaptive <- organise_output(fit = MLE.adaptive,
                                truePars = c(R0, Disp, obsPr),
                                name = "adaptive")

res.exact <- organise_output(fit = MLE.exact,
                             truePars = c(R0, Disp, obsPr),
                             name = "exact")

res.adaptive.fixed <- organise_output(fit = MLE.adaptive.fixed,
                             truePars = c(R0, Disp, obsPr),
                             name = "adaptive_fixed")

rbind(res.adaptive, res.exact, res.adaptive.fixed)



MLE.adaptive@min
MLE.exact@min
MLE.adaptive.fixed@min
## Likelihood at 'true' values:
marginal_loglikelihood_adaptive(data = data, pars = c(R0, Disp, obsPr),
                                verbose = TRUE)
marginal_loglikelihood_exact(data = data, pars = c(R0, Disp, obsPr))
############################
### Now, how good is the fitting on the uncorrupted data?

perfect.data <- compress_counts(simu$obs_y)

marginal_loglikelihood_direct <- function(data, pars){
  logliks <- unlist(
    lapply(1:data$K, function(i){
      data$freq[i] * dR0(y = data$count[i], r = pars[1], w = pars[2], log = TRUE)
    })
  )
  logLikelihood <- sum(logliks)
  return(logLikelihood)
}

nll_perfect <- function(logit_R0, log_omega){
  theta <- map_back(c(logit_R0, log_omega, 0))
  -marginal_loglikelihood_direct(data = perfect.data, pars = theta)
}

MLE.perfect <- bbmle::mle2(minuslogl = nll_perfect,
                           start = list(logit_R0 = inits[1],
                                        log_omega = inits[2]),
                           method = "L-BFGS-B"
)

map_back(c(coef(MLE.perfect), 0))
