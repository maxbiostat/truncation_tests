library(sumR)
source("../aux/aux.r")
source("../aux/powerLaw_aux.r")
delta <- .Machine$double.eps
Phi0 <- .1
Phi <- c(Phi0, 0)
Theta <- c(2, 1,  Phi)
lgL <- log(1-delta)
Eps <- delta
TrueValue <- logDiffExp(0, -Phi[1]) + 2*log(pi) - log(6)

result.R <- compare_approximations(compute_lterm = PL_diff_function_direct,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)


result.R