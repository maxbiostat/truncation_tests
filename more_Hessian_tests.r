library(sumR)
library(numDeriv)
source("aux/aux.r")
source("Erlang_aux.r")
source("erlang_hessian.R")

logf_x_mu <- function(m){
  sum_erlang(c(m, bb, xx))  
}
logf_x_mu <- Vectorize(logf_x_mu)

exactLogDeriv_f10 <- function(m){
  logSi <- sum_erlang(parameters = c(m, bb, xx),
                      epsilon = 1E-16)
  logS1i <- S1(c(m, bb, xx))
  
  return(logDiffExp(logS1i, logSi)-logSi)
}
exactLogDeriv_f10 <- Vectorize(exactLogDeriv_f10)

exactLogDeriv_f20 <- function(m){
  logSi <- sum_erlang(parameters = c(m, bb, xx),
                      epsilon = 1E-16)
  logS1i <- S1(c(m, bb, xx))
  logS2i <- S2(c(m, bb, xx))
  
  return(logDiffExp(logS2i, logS1i-log(m))-logSi )
}
exactLogDeriv_f20 <- Vectorize(exactLogDeriv_f20)

logf_x_beta <- function(b){
  sum_erlang(c(mm, b, xx))  
}
logf_x_beta <- Vectorize(logf_x_beta)

exactLogDeriv_f01 <- function(b){
  logSi <- sum_erlang(parameters = c(mm, b, xx),
                      epsilon = 1E-16)
  logS3i <- S3(c(mm, b, xx))
  return(logDiffExp(logS3i-log(b), log(xx) + logSi) - logSi )
}
exactLogDeriv_f10 <- Vectorize(exactLogDeriv_f10)

exactLogDeriv_f02 <- function(b){
  logSi <- sum_erlang(parameters = c(mm, b, xx),
                      epsilon = 1E-16)
  logS3i <- S3(c(mm, b, xx))
  logS4i <- S4(c(mm, b, xx))
  res <- logDiffExp(logS4i, logS3i-2*log(b)) - logSi 
  return(res)
}
exactLogDeriv_f02 <- Vectorize(exactLogDeriv_f02)

mm <- 1000
xx <- 1000
bb <- 2.1

nn10 <- numDeriv::grad(logf_x_mu, x = mm)
log(abs(nn10))
exactLogDeriv_f10(mm)

nn20 <- numDeriv::hessian(logf_x_mu, x = mm)
log(abs(nn20))
exactLogDeriv_f20(mm)

nn01 <- numDeriv::grad(logf_x_beta, x = bb)
log(abs(nn01))
exactLogDeriv_f01(bb)

nn02 <- numDeriv::hessian(logf_x_beta, x = bb)
log(abs(nn02))
exactLogDeriv_f02(bb)
