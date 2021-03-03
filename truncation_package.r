log_sum_exp <- function(x){
  # log(sum(exp(x - max(x)))) + max(x) ## unstable
  ans <- matrixStats::logSumExp(lx = x)
  return(ans)
}
##
log_diff_exp <- function(x, y) {
  # if(x <= y) stop("computing the log of a negative number")
  if(y == -Inf){
    return (x)
  }else{    
    return(x + log (-expm1(y-x)))
  }
}
##
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- log_diff_exp(max(x, y), min(x, y))
  }else{
    ans <- sgn * exp(log_diff_exp(max(x, y), min(x, y)))
  }
  return(ans)
}
##
relative_difference <- function(x, y){
  obj <- all.equal(x, y, tolerance = 0)
  ans <- suppressWarnings(as.numeric(gsub("Mean relative difference: ", "", obj)))
  if(is.na(ans)) ans <- 0
  return(ans)
}
approx_naive <- function(compute_lterm, theta, epsilon, n0 = 0, max_iter = 1E5){
  lterms <- rep(NA, max_iter + 1) ## avoid creating objects of variable size
  n <- n0
  lterm <- compute_lterm(n, theta) 
  lterms[n-n0 + 1] <- lterm
  n <- n + 1
  lterm <- compute_lterm(n, theta) 
  lterms[n-n0 + 1] <- lterm
  while(lterm > log(epsilon)){
    n <- n + 1
    lterm <- compute_lterm(n, theta) 
    lterms[n-n0 + 1] <- lterm
    if(n > max_iter){
      warning("Exceeded maximum number of iterations")
      ans <- log_sum_exp(na.omit(lterms))
      return(list(iter = n-n0, sum = ans)) 
    }
  }
  ans <- log_sum_exp(na.omit(lterms))
  return(list(iter = n-n0, sum = ans))
}
#
approx_naive_under <- function(compute_lterm, theta, epsilon, n0 = 0, max_iter = 1E5){
  lterms <- rep(NA, max_iter + 1) ## avoid creating objects of variable size
  n <- n0
  lterm <- compute_lterm(n, theta) 
  lterms[n-n0 + 1] <- lterm
  n <- n + 1
  lterm <- compute_lterm(n, theta)
  lterms[n-n0 + 1] <- lterm
  while(lterm > -700){
    n <- n + 1
    lterm <- compute_lterm(n, theta)
    lterms[n-n0 + 1] <- lterm
    if(n > max_iter){
      warning("Exceeded maximum number of iterations")
      ans <- log_sum_exp(na.omit(lterms))
      return(list(iter = n-n0, sum = ans)) 
    }
  }
  ans <- log_sum_exp(na.omit(lterms))
  return(list(iter = n-n0, sum = ans))
}
##
approx_doubling <- function(compute_lterm, theta, epsilon,
                            N_start = 20, c = 2, n0 = 0, max_iter = 1E5){
  ## TODO: make multicore with parallel.
  if(c <=1) stop("c needs to be greater than 1")
  if(round(c*N_start) > max_iter) stop("c*N_start already exceeds max_iter")
  leps <- log(epsilon)
## Phase 1: sum up the first batch.
  N <- N_start
  S1 <- log_sum_exp(compute_lterm(k = n0:N, theta))
  N <- round(c*N)
  if(N > max_iter/2){
    warning("Exceeded maximum number of iterations")
    return(list(iter = N-n0, sum = S2)) 
  }
  S2 <- log_sum_exp(compute_lterm(k = n0:N, theta))
  diff <- log_diff_exp(S2, S1)
  ## Phase 2: multiply ("double") until difference is below target.
  while(diff > leps){
    S1 <- S2
    N <- round(c*N)
    if(N > max_iter/2){
      warning("Exceeded maximum number of iterations")
      return(list(iter = N-n0, sum = S2)) 
    }
    S2 <- log_sum_exp(compute_lterm(k = n0:N, theta))
    diff <- log_diff_exp(S2, S1)
  }
  return(list(iter = N, sum = S2))
}
##
get_logz <- function(la, lap1){
  lnum <- la + lap1
  ldenom <- log_diff_exp(la, lap1)
  return(lnum-ldenom)  
}
##
get_delta <- function(lZ, lan, lR){
  lS <- lan + lR
  ## if lZ > lS, the sequence ratios decreases to L.
  ## if lS > lZ, the sequence rations increase to L.
  ans <- log_diff_exp(max(lZ, lS), min(lZ, lS))
  return(ans)
}
##
get_lratio <- function(la, lap1){
  ans <- lap1 - la
  return(ans)
}
##
approx_adaptive <- function(compute_lterm, theta, epsilon, n0 = 0, max_iter = 1E5, logL = -Inf){
  lterms.bulk <- rep(NA, max_iter + 1) ## avoid creating objects of variable size
  leps <- log(epsilon)
  target.leps  <- leps + log(2)
  logR <- logL - log_diff_exp(0, logL)
  #######
  n <- n0
  lterm <- compute_lterm(n, theta)
  lterms.bulk[n-n0+1] <- lterm
  old.lterm <- lterm
  n <- n + 1 ## how many calls to compute_lterm there have been.
  lterm <- compute_lterm(n, theta)
  lterms.bulk[n-n0+1] <- lterm
  lratio <- get_lratio(la = old.lterm, lap1 = lterm)
  ## Phase 1: sum up to the mode.
  while(lratio >= 0){
    old.lterm <- lterm
    n <- n + 1
    lterm <- compute_lterm(n, theta)
    lterms.bulk[n-n0+1] <- lterm
    lratio <- get_lratio(la = old.lterm, lap1 = lterm)
    if(n > max_iter){
      warning("Exceeded maximum number of iterations in the bulk phase")
      S.bulk <- log_sum_exp(na.omit(lterms.bulk))
      return(list(iter = n, sum = )) 
    }
  }
  nmode <- length(na.omit(lterms.bulk))
  nt <- n-nmode+1
  S.bulk <- log_sum_exp(na.omit(lterms.bulk))
  ## Phase 2: sum down the tail.
  lterms.tail <- rep(NA, max_iter-nmode + 1) ## avoid creating objects of variable size
  old.lterm <- lterm
  n <- n+1
  lterm <- compute_lterm(k = n, theta)
  lterms.tail[n-nmode-n0+1] <- lterm
  lz <- get_logz(la = old.lterm, lap1 = lterm)
  ldelta <- get_delta(lZ = lz, lan = old.lterm, lR = logR)
  while(ldelta > target.leps){
    old.lterm <- lterm
    n <- n+1
    nt <- n-nmode
    lterm <- compute_lterm(k = n, theta)
    lterms.tail[nt-n0+1] <- lterm
    lz <- get_logz(la = old.lterm, lap1 = lterm)
    ldelta <- get_delta(lZ = lz, lan = old.lterm, lR = logR)
    if(n > max_iter){
      warning("Exceeded maximum number of iterations in the tail phase")
      S.tail <- log_sum_exp(na.omit(lterms.tail))
      lVn <- log_sum_exp(c(lz, old.lterm + logR))
      S.tail.hat <- log_sum_exp(c(S.tail, lVn-log(2)))   
      ans <-  log_sum_exp(c(S.bulk, S.tail.hat)) 
      return(list(iter = n, sum = ans, ld = ldelta)) 
    }
  }
  S.tail <- log_sum_exp(na.omit(lterms.tail))
  lVn <- log_sum_exp(c(lz, old.lterm + logR))
  S.tail.hat <- log_sum_exp(c(S.tail, lVn-log(2)))
  # ans <- log_sum_exp(c(S.bulk, S.tail.hat))
  ans <- log_sum_exp(c(na.omit(lterms.bulk), na.omit(lterms.tail), lz - log(2), old.lterm + logR - log(2)))
  return(list(mode = nmode, Sb = S.bulk, St = S.tail, St.est = S.tail.hat,
              lz = lz, lh = old.lterm + logR, ld = ldelta,
              iter = n-n0, sum = ans))
}
