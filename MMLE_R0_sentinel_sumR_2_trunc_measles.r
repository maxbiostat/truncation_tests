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
organise_output <- function(fit, name, conf.level = 0.95){
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
    parameter = c("R[0]", "omega", "nu"),
    point = CENTER,
    lwr = LWR,
    upr = UPR,
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
Disp <- .5

USA <- list(
  count = c(1, 2, 3, 4, 5, 6, 8, 9, 11, 13, 15, 33),
  freq = c(122, 13, 10, 6, 5, 2, 2, 1, 1, 1, 1, 1)
)
USA$K <- length(USA$count)

Canada <- list(
  count = c(1, 2, 3, 4, 6, 8, 17, 30, 155),
  freq = c(35, 5, 3, 1, 1, 1, 1, 1, 1)
)
Canada$K <- length(Canada$count)

country <- "USA"
if(country == "Canada"){
  data <- Canada
}else{
  data <- USA
}

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

res.adaptive <- organise_output(fit = MLE.adaptive,
                                name = "adaptive")

res.exact <- organise_output(fit = MLE.exact,
                             name = "exact")



rbind(res.adaptive, res.exact)



MLE.adaptive@min
MLE.exact@min

Dispersions <- c(seq(0.01, .3, length.out = 10), .4, .5,
                 1, 1.5, 2)
fits <- lapply(Dispersions, function(w){
  MLE.adaptive.fixed <- bbmle::mle2(minuslogl = nll_adaptive,
                                    start = list(logit_R0 = inits[1],
                                                 logit_nu = inits[3]),
                                    fixed = list(log_omega = log(w)),
                                    lower = lower_l,
                                    upper = upper_l,
                                    method = "L-BFGS-B"
  )
  res.adaptive.fixed <- organise_output(fit = MLE.adaptive.fixed,
                                        name = "adaptive_fixed")
  res.adaptive.fixed$method <- NULL
  res.adaptive.fixed$omega <- w
  return(res.adaptive.fixed)
})

sensitivity.results <- do.call(rbind, fits)
subset(sensitivity.results, parameter == "R_0")
subset(sensitivity.results, parameter == "nu")

library(ggplot2)

ggplot() + 
  geom_pointrange(data = subset(sensitivity.results, parameter != "omega"),
                  mapping = aes(x = omega, y = point,
                                ymin = lwr, ymax = upr)) +
  geom_vline(xintercept = as.numeric(res.adaptive[2, 2:4]),
             linetype = "dotted") + 
  scale_x_continuous(expression(omega)) +
  scale_y_continuous("") +
  ggtitle(country) +
  facet_grid(parameter~., scales = "free_y", labeller = label_parsed) +
  theme_bw(base_size = 20)

