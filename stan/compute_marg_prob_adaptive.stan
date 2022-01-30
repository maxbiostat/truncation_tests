real[] margProb_DOP_adaptive(real mu, real phi, real p,
                             real x, int x_i, real eps, int M) {
  return(infiniteSumToThreshold_2({mu, phi, p, x}, eps, M, x_i));
}
