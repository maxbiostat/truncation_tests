# remotes::install_github("GuidoAMoreira/sumR")
library(sumR)
source("compare_approximations.r")
######################################
plaw_diff_lpmf <- function(k, theta){
  a <- theta[1]
  m <- theta[2]
  p0 <- theta[3]
  p1 <- theta[4]
  ans <- sapply(k, function(k) ifelse(k < m,  -Inf,
                                  -a*log(k) + log_diff_exp(0, -p0 - p1* (k-m)))
  )
  return(ans)
}

Alpha <- 2
Xm <- 1
Phi0 <- 0.05
Phi1 <- 0.0 # This can't really be 0.0 when Xm > 1,
#             because the Lerch doesn't converge

Theta <- c(Alpha, Xm, Phi0, Phi1)
Eps <- 1E-16
M <- 2E5
x0 <- Xm
lgL <- log(.99)
Ns <- 10

is_special_case <- all(Alpha == 2, Xm == 1, Phi1 == 0)

if(is_special_case){
  cat("Doing special case \n")
  TrueValue <- log_diff_exp(0, -Phi0) + 2*log(pi) - log(6)
}else{
  ## "True" value using specialised functions from VGAM
  if(Phi1 == 0) stop("Phi1 0 can't be zero when x_min > 1")
  logHZeta <- log(VGAM::zeta(Alpha, shift = Xm))
  logLerch <- log(VGAM::lerch(x = exp(-Phi1),
                              s = Alpha, v = Xm, tolerance = Eps))- Xm*Phi1  
  TrueValue <-  log_diff_exp(logHZeta, -Phi0 + Phi1*Xm + logLerch)  
}

lps <- plaw_diff_lpmf(k = 0:M, theta = Theta)
matrixStats::logSumExp(lps)


TrueValue
result <- compare_approximations(compute_lterm = plaw_diff_lpmf,
                                 theta = Theta,
                                 exact = TrueValue,
                                 n0 = x0, logL = lgL,
                                 Nstart = Ns,
                                 eps = Eps,
                                 max_iter = M/2)
result
