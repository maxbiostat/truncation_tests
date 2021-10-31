real[] log_Z_DOP_adaptive(real mu, real phi, real eps, int M) {
  return( infiniteSumToThreshold({mu, phi}, eps, M, 0) );
}
