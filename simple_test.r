# remotes::install_github("GuidoAMoreira/sumR")
library(sumR)
source("compare_approximations.r")
######################################
log_geom_series <- function(k, p) return(k * log1p(-p[1]))

p <- .01
Theta <- p
Eps <- 1E-6
M <- 2E5
x0 <- 0
lgL <- log1p(-p)
Ns <- 10

TrueValue <- -log(p)

result <- compare_approximations(compute_lterm = log_geom_series,
                                 theta = Theta,
                                 exact = TrueValue,
                                 n0 = x0,
                                 logL = lgL,
                                 Nstart = Ns,
                                 eps = Eps,
                                 max_iter = M/2)
result