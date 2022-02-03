library(Rmpfr)

source("../aux/aux_rmpfr.r")
source("../aux/Poisson_aux.r")

Mu <- 11.8
R <- 2
Mu.prec <- Rmpfr::mpfr(Mu, 1024) ## try 11.8 and 11.9 and 12.8 and be amazed
R.prec <- Rmpfr::mpfr(R, 1024)
Theta <- c(Mu, R)
TrueValue <- R.prec*log(Mu.prec)
Eps <- .Machine$double.eps
lgL <- log(0)

result.R <- compare_approximations(compute_lterm = poisson_lfactmom_R,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)

result.C <- compare_approximations(compute_lterm = poisson_lfactmom_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)

result.preComp <- compare_approximations(
  compute_lterm = "poisson_fact_moment",
  theta = Theta,
  exact = TrueValue,
  eps = Eps)

result.R
result.C
result.preComp

res <- result.preComp$error
names(res) <- result.preComp$Method
res
