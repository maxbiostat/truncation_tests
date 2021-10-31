// Sum-To-Threshold inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
real[] infiniteSumToCap_2(real[] p, int Nmax, int maxIter, int n0) {
  vector[maxIter + 1] storeVal;
  int n = 2;
  int n0_ = n0;

  // Setting up first iterations
  storeVal[1] = logFunction_2(n0_, p);
  n0_ += 1;
  storeVal[2] = logFunction_2(n0_, p);
  
  // Find the maximum
  while (storeVal[n] > storeVal[n - 1]) {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction_2(n0_, p);
    if (n >= maxIter) return({log_sum_exp(storeVal[:n]), 1. * n}); // Return if maxIter is reached
  }
  // Start testing convergence after the maximum
  while ( n < Nmax )  {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction_2(n0_, p);
    if (n >= maxIter) break;
  }
  return {log_sum_exp(sort_asc(storeVal[:n])), 1. * n};
}
