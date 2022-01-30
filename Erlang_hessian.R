#### Hessian ####
Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double sss(long n, double *p)
{
  return n * (logl(p[0]) + logl(p[1]) + logl(p[2])) - logl(n) - 2 * lgammal(n);
}

// [[Rcpp::export]]
double series(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];
  
  r = infiniteSum(sss, parameter, -INFINITY, 0, epsilon, maxIter, 1, &n);
  
  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double s1(long n, double *p)
{
  return n * (logl(p[0]) + logl(p[1]) + logl(p[2])) - 2 * lgammal(n);
}

// [[Rcpp::export]]
double series1(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(s1, parameter, -INFINITY, 0, epsilon, maxIter, 1, &n);
  
  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double s2(long n, double *p)
{
  return n * (logl(p[0]) + logl(p[1]) + logl(p[2])) + logl(n - 1) - 2 * lgammal(n);
}

// [[Rcpp::export]]
double series2(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(s2, parameter, -INFINITY, 0, epsilon, maxIter, 1, &n);
  
  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double s3(long n, double *p)
{
  return n * (logl(p[0]) + logl(p[1]) + logl(p[2])) + logl(n) - 2 * lgammal(n);
}

// [[Rcpp::export]]
double series3(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(s3, parameter, -INFINITY, 0, epsilon, maxIter, 1, &n);
  
  return (double)r;
}
')

erlangHessian <- function(x, mu, beta) {
  set.seed(123)
  mu2 <- mu * mu
  beta2 <- beta * beta
  expmu <- exp(-mu)
  
  serie <- sapply(x, function(xx) series(c(xx, mu, beta)))
  
  addSer <- function(s) sapply(1:length(x), function(i)
    exp(s(c(x[i], mu, beta)) - serie[i]))
  output <- matrix(0, 2, 2)
  
  s12 <- sum(addSer(series1) ^ 2)
  s2 <- sum(addSer(series2))
  output[1, 1] <- - length(x) * expmu / (1 - expmu) ^ 2 + (s2 - s12) / mu2
  output[2, 2] <- (s2 - s12) / beta2
  output[1, 2] <- (sum(addSer(series3)) - s12) / (mu * beta)
  output[2, 1] <- output[1, 2]
  
  output
}
