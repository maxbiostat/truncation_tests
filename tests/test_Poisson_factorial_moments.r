library(sumR)
library(tidyverse)
source("../aux/aux.r")
source("../aux/Poisson_aux.r")

#######
##$ Table

Miter <- 5E5
name.huge <- paste0('Fixed_', Miter)

spit_approx <- function(lambda, order, epsilon){
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = "poisson_fact_moment",
      theta = c(lambda, order),
      exact = order*log(lambda),
      max_iter = Miter,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    mu = lambda,
    r = order,
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
  Mu = c(0.5, 1, 10, 100),
  R = c(2, 5, 10),
  eps = 10^(c(0, 1, 4))*.Machine$double.eps
)

approximation.list <- parallel::mclapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(lambda = linha$Mu,
                order = linha$R,
                epsilon = linha$eps)      
  }, mc.cores = 8)

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt$row_id <- paste(approximation.dt$mu,
                                  approximation.dt$r,
                                  approximation.dt$target_error,
                                  sep = "_")

approximation.dt <- approximation.dt %>%
  group_by(row_id) %>%
  mutate(comparative_error = abs(error)/abs(error[method == name.huge]))

approximation.dt <- approximation.dt %>%
  mutate(comparative_success = comparative_error <= 1)

approximation.dt <- approximation.dt %>%
  mutate(actual_success = as.logical(success + comparative_success))

success <- aggregate(success~method, approximation.dt, mean)
success.max <- aggregate(success_max~method, approximation.dt, mean)
comparative <- aggregate(comparative_success~method, approximation.dt, mean)
actual <- aggregate(actual_success~method, approximation.dt, mean)

success.results <- plyr::join_all(list(success, success.max, comparative, actual),
                          by = "method",
                          type = "left",
                          match = "all")

success.results

to.present <- success.results
to.present$success_max <- NULL
print(xtable::xtable(to.present), include.rownames = FALSE)

aggregate(success~method+target_error, approximation.dt, mean)
aggregate(success_max~method+target_error, approximation.dt, mean)
aggregate(comparative_success~method+target_error, approximation.dt, mean)
aggregate(actual_success~method+target_error, approximation.dt, mean)

failThresh <- subset(approximation.dt,
       !actual_success
       & method == "Threshold")%>% print(n = 58)

failThresh

failBPair <- subset(approximation.dt,
       !actual_success 
       & method == "BoundingPair") %>% print(n = 28)

failBPair

suspicious <- subset(approximation.dt,
               !success 
               & comparative_error > 1.01)

subset(approximation.dt, row_id %in% suspicious$row_id)
