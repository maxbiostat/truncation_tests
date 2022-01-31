source("../aux/aux.r")
source("../aux/Poisson_aux.r")

Mu <- 11.8
R <- 2
Theta <- c(Mu, R)
TrueValue <- R*log(Mu)
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