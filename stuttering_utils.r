dR0 <- function(y, r, w, log = FALSE){
  c1 = lgamma(w*y + y-1)
  c2 = lgamma(w*y)
  c3 = lgamma( y+1)
  c4 = (y-1) * (log(r) - log(w))
  c5 = (w*y + y -1) * log(1 + r/w)
  dens = c1 - (c2+c3) + (c4-c5)
  if(!log){
    dens <- exp(dens)
  }
  return(dens)
}
#
rR0 <- function(n, r, w, D, UpperBound = 1e4){
  ys <- 1:UpperBound
  Ps <- dR0(y = ys, r = r, w = w)
  as.vector(
    matrix(
      sample(ys, n * D, prob = Ps, replace = TRUE),
      ncol = D, nrow = n
    )  
  )
}
get_marginalised_prob <- function(r, w, p, verbose = FALSE){
  ## Get Pr(S = 0 | r, w, p)
  a <- r/w
  b <- r/(w * (1 + r/w)^(w + 1)) * (1-p)
  solve_for_u <- function(b, w){
    ## We know b < 1 because a and omega are positive
    ## it suffices to find the positive root smaller than 1/omega, which
    ## is the positive root of the objective function obj_fun
    b_of_u <- function(u) u/(1 + u)^{w + 1}
    obj_fun <- function(cand) (b-b_of_u(cand))^2
    u <- optimise(obj_fun, c(0, 1/w), tol = 1e-25)$minimum
    return(u)
  }
  ustar <- solve_for_u(b, w = w)
  if(verbose) cat("u* is:", ustar, "\n")
  return( (a + 1)/a * ustar/(1 + ustar))
}
##
getMLE <- function(C) 1 - (1/mean(C)) ## maximum likelihood estimator using the mean of C
##
# Get (theoretical) mean and variance of the observed sequence clusters
get_mean_s_binomial <- function(r, w, p){
  p0 <- get_marginalised_prob(r, w, p)
  p/((1-p0)*(1-r))  
}
#
get_var_s_binomial <- function(r, w, p){
  p0 <- get_marginalised_prob(r, w, p)
  f <- 1/(1-r)
  g <- r*(1 + r/w)/(1-r)^3
  ans <-((1-p0)*(f*p*(1-p) + p^2*g + p^2*f^2) -p^2*f^2)/(1-p0)^2
  return(ans)
}
#
get_bias_R0_binomial <- function(r, w, p){
  ## if r, w and p are the true parameters, what is the bias if I estimate R0 by MLE from S (excluding zeroes)
  (1- 1/get_mean_s_binomial(r, w, p)) - r
}
#
get_mean_s_sentinel <- function(r, w, nu){
  get_expect_term_sentinel <- function(r, w, v){
    a <- r/w
    b <- r/(w * (1 + r/w)^(w + 1)) * (1-v)
    solve_for_eta <- function(b, w){
      ## We know b < 1 because a and omega are positive
      ## it suffices to find the positive root smaller than 1/omega, which
      ## is the positive root of the objective function obj_fun
      b_of_eta <- function(eta) eta/(1 + eta)^{w + 1}
      obj_fun <- function(cand) (b-b_of_eta(cand))^2
      eta <- optimise(obj_fun, c(0, 1/w), tol = 1e-25)$minimum
      return(eta)
    }
    eta_star <- solve_for_eta(b, w = w)
    # if(verbose) cat("v* is:", eta_star, "\n")
    return( (w + r)/r * eta_star/(1 + eta_star) * 1/(1 - w *eta_star) )
  }
  ##
  rho0 <- get_marginalised_prob(r, w, nu)
  expect <- get_expect_term_sentinel(r, w, nu)
  return( 1/((1-rho0)*(1-r)) - expect/(1-rho0) )
}
#
get_bias_R0_sentinel <- function(r, w, nu){
  ## if r, w and p are the true parameters, what is the bias if I estimate R0 by MLE from S (excluding zeroes)
  (1- 1/get_mean_s_sentinel(r, w, nu)) - r
}
##
compress_data <- function(simu){
  raw <- table(simu$s)
  ss <- as.numeric(names(raw))
  ns <- as.numeric(raw)
  return(list(s = ss, n = ns, D = length(ss), J = length(simu$s)))
}
##
generate_stuttering_data_one_level <- function(r0, omega, psi, K, seed = NULL){
  ## one-stage sampling scheme: sequence at rate psi
  if(!is.null(seed)) set.seed(seed)
  
  cs <- rR0(n = 1, r = r0, w = omega, D = K)
  true.N <- sum(cs) ## total number of cases
  
  all.x <- sapply(cs, function(x) rbinom(1, size = x, prob = psi))
  sequenced.pos <- which(all.x > 0)
  obs.x <- all.x[sequenced.pos]
  obs.S <- sum(obs.x)  ## number of cases sequenced
  obs.J <- length(obs.x)  
  return(
    list(
      trueN = true.N, # how many cases
      true.K = K,# how many clusters
      y = cs, ## all case clusters 
      sequenced = sequenced.pos,  ## which clusters were sequenced
      s = obs.x,  ## the observed sequence clusters
      Nseq = obs.S,  ## how many sequences
      Nclu = obs.J ## how many sequence clusters
    )
  )  
}
####
get_confidences <- function(map_est){
  ## For MAP estimates with a Hessian matrix, return confidence intervals
  ## assuming estimates are asymptotically normal
  p <- ncol(map_est$hessian)
  mus <- map_est$par[1:p]
  sds <- sqrt(diag(solve(-map_est$hessian)))
  return(data.frame(lwr = qnorm(c(.025), mean = mus, sd = sds),
                    mean = mus,
                    upr = qnorm(c(.975), mean = mus, sd = sds),
                    sd = sds))
}
################
R0_cluster_marginal_detection <- function(k, theta){
  R0 <- theta[1]
  omega <- theta[2]
  p <- theta[3]
  x <- theta[4]
  ans <- ifelse(k < x, -Inf,
                dR0(y = k, r = R0, w = omega, log = TRUE) +
                  dbinom(x = x, size = k, prob = p, log = TRUE))
  return(ans)
}

## Pr(X = x | R0, w, psi)
### pars = c(R0, w, psi)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  
  P0 <- get_marginalised_prob(
    r = pars[1], w = pars[2], p = pars[3]
  )
  log1mP0 <- log1p(-P0)
  
  logLL <- log(pars[1]) + log1p(-pars[3]) +
    (1 + pars[2]) * (log1p(pars[2]) - log(pars[1] + pars[2]) )
  if(adaptive){
    logMargProbKernel <- sumR::infiniteSum(logFunction =
                                             R0_cluster_marginal_detection,
                                           parameters = c(pars, x),
                                           logL = logLL,
                                           epsilon = eps,
                                           maxIter = 1E5)
    if(verbose) cat("Marginalised probability took", logMargProbKernel$n,
                    "iterations \n")
    
  }else{
    logMargProbKernel <- sumR::finiteSum(logFunction =
                                           R0_cluster_marginal_detection,
                                         parameters = c(pars, x),
                                         n = N_fix)
  }
  ans <- logMargProbKernel$sum - log1mP0
  return(ans)
}

marginal_loglikelihood <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE, N_fix = NULL,
                                   verbose = FALSE){
  
  logliks <- unlist(
    lapply(1:data$D, function(i){
      data$n[i] * marg_lik(x = data$s[i], pars = pars,
                              eps = eps, adaptive = adaptive,
                              N_fix = N_fix, verbose = verbose)
    })
  )
  return(sum(logliks))
}