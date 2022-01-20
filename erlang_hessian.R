#### Hessian ####
Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double erlang(long n, double *p)
{
  return n < 1 ? -INFINITY :
    R::dpois(n, p[0], 1) + R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double sum_erlang(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(erlang, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double f10_term(long n, double *p)
{
  return n < 1 ? -INFINITY :
    R::dpois(n-1, p[0], 1) + R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double S1(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(f10_term, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double f20_term(long n, double *p)
{
  return n < 1 ? -INFINITY :
    -2*log(p[0]) + R::dpois(n, p[0], 1) + log((p[0]-n)*(p[0]-n)) +
    R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double S2(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(f20_term, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double f01_term(long n, double *p)
{
  return n < 1 ? -INFINITY :
     log(n) + R::dpois(n, p[0], 1) + R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double S3(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(f01_term, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double f02_term(long n, double *p)
{
  return n < 1 ? -INFINITY :
    log((p[2]*p[1]-n)*(p[2]*p[1]-n))-2*log(p[1]) + R::dpois(n, p[0], 1) +
    R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double S4(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(f02_term, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

erlangHSingle <- function(xi, mu, beta){
  
  logF <-  sum_erlang(parameters = c(mu, beta, xi),
                    epsilon = 1E-16)
  
  logS1i <- S1(c(mu, beta, xi))
  logS2i <- S2(c(mu, beta, xi))
  logS3i <- S3(c(mu, beta, xi))
  logS4i <- S4(c(mu, beta, xi))
  #
  logFF10 <- logDiffExp(logF, logS1i)-logF
  logFF20 <- logDiffExp(logS2i, logS1i-log(mu))-logF
  logFF01 <- logDiffExp(logS3i-log(beta), log(xi) + logF) - logF 
  logFF02 <- logDiffExp(logS4i, logS3i-2*log(beta)) - logF 
  logFF11 <- logDiffExp(logFF10, logFF01)
  #
  Ai <- exp(
    logDiffExp(logF + logFF20, 2*logFF10)-2*logF
  ) 
  Bi <- exp( logFF11 - logF )
  Ci <- exp(
    logDiffExp(logF + logFF02, 2*logFF01)-2*logF
  ) 
  result <- matrix(c(Ai, Bi, Bi, Ci), 2, 2)
  return(result)
}
erlangHessian <- function(xs, muMLE, betaMLE) {
  matlist <- lapply(xs, function(x){
    erlangHSingle(xi = x, mu = muMLE, beta = betaMLE)
  })
  output <- Reduce('+', matlist)
  return(output)
  
}
