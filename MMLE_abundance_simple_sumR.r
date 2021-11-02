library(sumR)
library(gamlss.dist)
library(stats4)
source("double_poisson_binomial_aux.r")

### 
map_Rn <- function(x){
  return(
    c(
      log(x[1]),
      log(x[2]),
      arm::logit(x[3])
    )
  )
}
map_back <- function(y){
  return(
    c(
      exp(y[1]),
      exp(y[2]),
      arm::invlogit(y[3])
    )
  )
}
####################
Mu <- 40
Phi <- 1
PP <- .22
J <- 1000
simu <- simulate_obsdata(n = J, mu = Mu, phi = Phi, p = PP)
data <- compress_counts(simu$obs_x)

##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(Mu, Phi, PP))

## Profile likelihood: Mu
f_mu <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(x, Phi, PP))
}
f_mu <- Vectorize(f_mu)
mus <- seq(from = 0.1*Mu, to = 10*Mu, length.out = 20)
plot(mus, f_mu(mus))
abline(v = Mu, lwd = 2, lty = 2)

## Profile likelihood: Phi
f_phi <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(Mu, x, PP))
}
f_phi <- Vectorize(f_phi)
phis <- seq(from = 0.1*Phi, to = 5*Phi, length.out = 20)
plot(phis, f_phi(phis))
abline(v = Phi, lwd = 2, lty = 2)

## Profile likelihood: p
f_p <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(Mu, Phi, x))
}
f_p <- Vectorize(f_p)
ps <- seq(from = 0.1, to = 1, length.out = 20)
plot(ps, f_p(ps))
abline(v = PP, lwd = 2, lty = 2)

####

lik_adaptive <- function(log_mu, log_phi, logit_p){
  theta <- map_back(c(log_mu, log_phi, logit_p))
  -marginal_loglikelihood(data, pars = theta)
}
  
inits <- map_Rn(c(runif(1, 1, 100),
                  runif(1, .1, 10),
                  runif(1))
                )
MLE <- stats4::mle(minuslog = lik_adaptive,
            start = list(log_mu = inits[1],
                         log_phi = inits[2],
                         logit_p = inits[3]),
            lower = map_Rn(c(.1, .1, .01)),
            upper = map_Rn(c(1E4, 10, .99)))

MLE.alt <- bbmle::mle2(minuslogl = lik_adaptive,
                       start = list(log_mu = inits[1],
                                    log_phi = inits[2],
                                    logit_p = inits[3]),
                       method = "L-BFGS-B"
                       )

summary(MLE)

map_back(coef(MLE))
map_back(coef(MLE.alt))
c(Mu, Phi, PP)
MLE@details

# CIs <- confint(MLE)
