library(tidyverse)
library(sumR)
source("../aux/aux.r")
source("../aux/powerLaw_aux.r")
####
##$ Table
run <- TRUE

Miter <- 5E5
name.huge <- paste0('Fixed_', Miter)

spit_approx <- function(p0, L, epsilon){
  B <- max(L/(1-L), 1)
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = PL_diff_function_direct,
      theta =  c(2, 1, p0, 0),
      logL = log(L),
      exact = logDiffExp(0, -p0) + 2*log(pi) - log(6),
      batch_size = B + 1,
      max_iter = Miter,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    phi_0 = p0,
    L = L,
    target_error = ans$target_error,
    error = ans$error,
    error_max = ans$error_max,
    success = abs(ans$error) <= ans$target_error,
    success_max = abs(ans$error_max) <= ans$target_error,
    n_iter = ans$n_evaluations,
    method = ans$Method
  )
  return(out)
}

approximation.grid <- expand.grid(
  phi0 = c(.1, 1, 10, 100),
  L = c(0,.5, .9, .99, 1-.Machine$double.eps),
  eps = 10^(c(0, 1, 4))*.Machine$double.eps
)

if(run){
  approximation.list <- parallel::mclapply(
    1:nrow(approximation.grid),
    function(i){
      linha <- approximation.grid[i, ] ## linha == Portuguese for row
      spit_approx(p0 = linha$phi0,
                  L = linha$L,
                  epsilon = linha$eps)      
    }, mc.cores = 10)
  
  save(approximation.list,
            file = "../saved_data/errorAnalysis_powerLaw.RData")
}else{
  load("../saved_data/errorAnalysis_powerLaw.RData")
}

approximation.dt <- do.call(rbind, approximation.list)

#### Investigating the instances in which adaptive/naive fails
#### but the sum with a huge (3E5) number of terms does not.

approximation.dt$row_id <- paste(approximation.dt$phi_0,
                                 approximation.dt$L,
                                 approximation.dt$target_error,
                                 sep = "_")


approximation.dt <- approximation.dt %>%
  group_by(row_id) %>%
  mutate(comparative_error = abs(error)/
           max(abs(error[method == name.huge]), target_error))

approximation.dt <- approximation.dt %>%
  mutate(comparative_success = comparative_error <= 1)

approximation.dt <- approximation.dt %>%
  mutate(actual_success = as.logical(success + comparative_success))

aggregate(success~method, approximation.dt, mean)
aggregate(comparative_success~method, approximation.dt, mean)

aggregate(success~method+target_error, approximation.dt, mean)
aggregate(comparative_success~method+target_error, approximation.dt, mean)

aggregate(comparative_success~method+L,
          approximation.dt,
          mean)

subset(approximation.dt,
       !actual_success 
       & method == "Threshold") %>% print(n = 45)

subset(approximation.dt,
       !actual_success
       & method == "BoundingPair") %>% print(n = 28)
