library(sumR)
library(bbmle)
source("Erlang_aux.r")
source("aux/aux.r")
###########
Mu <- 1500
B <- 1
J <- 5
data <- simulate_obsdata(n = J, mu = Mu, b = B,
                         seed = 1234)

ll.bessel.precomp <- marginal_loglikelihood(data = data,
                                            pars = c(Mu, B+B/2),
                                            eps = .Machine$double.eps,
                                            verbose = TRUE)

ll.bessel.fast <- marginal_loglikelihood_fast(data = data,
                                              pars = c(Mu, B+B/2),
                                              eps = .Machine$double.eps,
                                              verbose = TRUE)

ll.full <- marginal_loglikelihood_full(data = data,
                                       pars = c(Mu, B+B/2),
                                       eps = .Machine$double.eps,
                                       verbose = TRUE)

ll.full.fast <- marginal_loglikelihood_full_fast(data = data,
                                                 pars = c(Mu, B+B/2),
                                                 eps = .Machine$double.eps,
                                                 verbose = TRUE)
ll.bessel.fixed <- marginal_loglikelihood(data = data,
                                          pars = c(Mu, B+B/2),
                                          adaptive = FALSE,
                                          N_fix = 4E5,
                                          eps = .Machine$double.eps,
                                          verbose = TRUE)

ll.bessel.precomp
ll.bessel.fast
ll.full
ll.full.fast
ll.bessel.fixed

robust_difference(ll.bessel.fast, ll.bessel.precomp)
robust_difference(ll.bessel.precomp, ll.full)
robust_difference(ll.full, ll.full.fast)