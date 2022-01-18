library(tidyverse)
library(sumR)
source("../aux/aux.r")
#######
PL_diff_function_direct <- function(k, theta){
  alpha <- theta[1]
  xmin <- theta[2]
  phi <- theta[3:4]
  ans <- sapply(k, function(k) ifelse(k < xmin,  -Inf,
                                      -alpha*log(k) + logDiffExp(0, -phi[1] - phi[2] * (k-xmin)))
  )
  return(ans)
}

#######
delta <- .Machine$double.eps
Phi0 <- .1
Phi <- c(Phi0, 0)
Theta <- c(2, 1,  Phi)
lgL <- log(1-delta)
Eps <- delta
TrueValue <- logDiffExp(0, -Phi[1]) + 2*log(pi) - log(6)

result.R <- compare_approximations(compute_lterm = PL_diff_function_direct,
                                 theta = Theta,
                                 exact = TrueValue,
                                 eps = Eps,
                                 logL = lgL)


result.R

####
##$ Table
spit_approx <- function(p0, L, epsilon){
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = PL_diff_function_direct,
      theta =  c(2, 1, p0, 0),
      logL = log(L),
      exact = logDiffExp(0, -p0) + 2*log(pi) - log(6),
      max_iter = 3E5,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    phi_0 = p0,
    L = L,
    target_error = ans$target_error,
    error = ans$error,
    success = abs(ans$error) <= ans$target_error,
    n_iter = ans$n_evaluations,
    method = ans$Method
  )
  return(out)
}

approximation.grid <- expand.grid(
  phi0 = c(.1, 1, 10, 100),
  L = c(0,.5, .9, .99, 1-delta),
  eps = 10^(0:6)*delta
)

run <- FALSE 

if(run){
  approximation.list <- parallel::mclapply(
    1:nrow(approximation.grid),
    function(i){
      linha <- approximation.grid[i, ] ## linha == Portuguese for row
      spit_approx(p0 = linha$phi0,
                  L = linha$L,
                  epsilon = linha$eps)      
    }, mc.cores = 10)
  
  approximation.dt <- do.call(rbind, approximation.list)
  
  write.csv(approximation.dt,
            file = "../saved_data/errorAnalysis_powerLaw.csv",
            row.names = FALSE)
}else{
  approximation.dt <- read.csv("../saved_data/errorAnalysis_powerLaw.csv")
}

#### Investigating the instances in which adaptive/naive fails
#### but the sum with a huge (3E5) number of terms does not.

approximation.dt$row_id <- paste(approximation.dt$phi_0,
                                 approximation.dt$L,
                                 approximation.dt$target_error,
                                 sep = "_")


approximation.dt <- approximation.dt %>%
  group_by(row_id) %>%
  mutate(comparative_error = abs(error)/
           max(abs(error[method == 'Fixed_3e+05']), target_error))

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
