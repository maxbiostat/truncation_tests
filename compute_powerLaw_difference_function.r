source("truncation_package.r")
source("truncation_aux.r")

PL_diff_function_direct<- function(k, theta){
  alpha <- theta[1]
  xmin <- theta[2]
  phi <- theta[3:4]
  ans <- sapply(k, function(k) ifelse(k < xmin,  -Inf,
                                      -alpha*log(k) + log_diff_exp(0, -phi[1] - phi[2] * (k-xmin)))
  )
  return(ans)
}

Lerch_term <- function(k, theta){
  alpha <- theta[1]
  xmin <- theta[2]
  phi <- theta[3:4]
  ans <- sapply(k, function(k) -alpha*log(k) - phi[2] * k)
  return(ans)
}

is_special_case <- TRUE
phi0 <- .25
if(is_special_case){
  Alpha <- 2
  Phi <- c(phi0, 0)
  minX <- 1
}else{
  Alpha <- 2.2
  Phi <- c(phi0, .02)
  minX <- 2
}

Theta <- c(Alpha, minX, Phi)
Eps <- 1E-10


if(is_special_case){
  TrueValue <- log_diff_exp(0, -Phi[1]) + 2*log(pi) - log(6)
}else{
  ## "True" value using specialised functions from VGAM
  logHZeta <- log(VGAM::zeta(Alpha, shift = minX))
  logLerch <- log(VGAM::lerch(x = exp(-Phi[2]), s = Alpha, v = minX, tolerance = Eps))- minX*Phi[2]  
  TrueValue <-  log_diff_exp(logHZeta, -Phi[1] + Phi[2]*minX + logLerch)  
}
##
TrueValue

result.direct <- compare_approximations(PL_diff_function_direct, theta = Theta,
                       exact = TrueValue, epsilon = Eps, n0 = minX, max_iter = 2E5, logL = log(.9999))

result.direct
