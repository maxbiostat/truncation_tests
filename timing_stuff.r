timing <- bench::mark(
  naive = sumR::infiniteSum(logFunction = COMP_lpmf,
                                 parameters = Theta,
                                 eps = Eps, maxIter = M/2,
                                 n0 = x0, forceAlgorithm = 1),
  doubling = sumR::infiniteSum_cFolding_C(logFunction = COMP_lpmf,
                               parameters = Theta,
                               eps = Eps,
                               maxIter = M/2,
                               n0 = x0,
                               N_start = Ns),
  adaptive = sumR::infiniteSum(logFunction = COMP_lpmf,
                               parameters = Theta,
                               eps = Eps,
                               logL = lgL,
                               maxIter = M/2,
                               n0 = x0,
                               forceAlgorithm = 2),
  automatic = sumR::infiniteSum(logFunction = COMP_lpmf,
                               parameters = Theta,
                               eps = Eps,
                               logL = lgL,
                               maxIter = M/2,
                               n0 = x0),
  check = FALSE
)

timing
