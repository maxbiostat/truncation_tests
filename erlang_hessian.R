#### Hessian ####
Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumRAPI.h>

long double sss(long n, double *p)
{
  return n * (logl(p[0]) + logl(p[1]) + logl(p[2])) - logl(n) -
    2 * lgammal(n);
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

#include <sumRAPI.h>

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

#include <sumRAPI.h>

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

#include <sumRAPI.h>

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

erlangHessian <- function(x, muMLE, betaMLE) {
  muMLE2 <- muMLE * muMLE
  betaMLE2 <- betaMLE * betaMLE
  expmu <- exp(-muMLE)
  onemexpmu <- 1 - expmu
  k <- expmu / onemexpmu + (expmu / onemexpmu) ^ 2
  output <- matrix(0, 2, 2)
  output[1, 1] <- length(x) * k
  for (i in 1:length(x)) {
    sss <- series(c(x[i], muMLE, betaMLE))
    s1Oversss2 <- (series1(c(x[i], muMLE, betaMLE)) / sss) ^ 2
    s2Oversss <- series2(c(x[i], muMLE, betaMLE)) / sss
    s3Oversss <- series3(c(x[i], muMLE, betaMLE)) / (sss * muMLE * betaMLE)
    output[1, 1] <- output[1, 1] - s1Oversss2 / muMLE2 + s2Oversss / muMLE2
    output[2, 2] <- output[2, 2] - s1Oversss2 / betaMLE2 + s2Oversss / betaMLE2
    output[1, 2] <- output[1, 2] - s1Oversss2 / (muMLE * betaMLE) + s3Oversss
  }
  output[2, 1] <-  output[1, 2]
  output
}
