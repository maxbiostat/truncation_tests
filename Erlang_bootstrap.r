erlang_bootstrap_once <- function(dt, 
                                  mll,
                                  mle,
                                  type = "nonparametric"){
  if(type == "nonparametric"){
    indices <- sample(seq_len(dt$K),
                      size = dt$K,
                      replace = TRUE)
    dd <- list(
      obs_x = dt$obs_x[indices],
      K = dt$K
    )
  }else{
    dd <- simulate_obsdata(n = dt$K,
                           mu = mle[1], b = mle[2])
  }
  obj <- get_estimates(data = dd, minusloglik_fun = mll)
  return(
    data.frame(
      mu_est = obj$results$point[1],
      beta_est = obj$results$point[2]
    )
  )
}

erlang_bootstrap <- function(dt, Nrep = 200, ncores,
                             mll,
                             mle,
                             type = "nonparametric"){
  result <- parallel::mclapply(
    1:Nrep,
    function(i){
      data.frame(
        erlang_bootstrap_once(dt = dt, mll = mll,
                              mle = mle, type = type),
        replicate = i
      )
    },
    mc.cores = ncores
  )
  return(do.call(rbind, result))
}