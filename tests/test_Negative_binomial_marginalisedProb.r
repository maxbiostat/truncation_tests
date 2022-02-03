library(sumR)
library(tidyverse)
source("../aux/aux.r")
source("../aux/NegBinomial_aux.r")

####
##$ Table
Miter <- 5E5
name.huge <- paste0('Fixed_', Miter)

spit_approx <- function(loc, disp, detp, obx, epsilon){
  
  logL <- log(loc) - matrixStats::logSumExp(c(log(loc), log(disp))) + log1p(-detp)
  L <- exp(logL)
  B <- max(L/(1-L), 1)
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = negativeBinomial_marginalised_lpmf_C,
      theta =  c(loc, disp, detp, obx),
      logL = logL,
      exact = dnbinom(x = obx, mu = detp*loc, size = disp, log = TRUE),
      max_iter = Miter,
      batch_size = B + 1,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    mu = loc,
    phi = disp,
    eta = detp,
    obs_x = obx,
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
  Mu = c(1, 10, 100),
  Phi = c(.1, .5, 1, 10),
  Eta = c(.01, .1, .5, .75),
  X = c(0, 5, 10),
  eps = 10^(c(0, 1, 4))*.Machine$double.eps
)

approximation.list <- parallel::mclapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(loc = linha$Mu,
                disp = linha$Phi,
                detp = linha$Eta,
                obx = linha$X,
                epsilon = linha$eps)      
  }, mc.cores = 8)

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt$row_id <- paste(approximation.dt$mu,
                                 approximation.dt$phi,
                                 approximation.dt$eta,
                                 approximation.dt$obs_x,
                                 approximation.dt$target_error,
                                 sep = "_")

approximation.dt <- approximation.dt %>%
  group_by(row_id) %>%
  mutate(comparative_error = abs(error)/abs(error[method == name.huge]))

approximation.dt <- approximation.dt %>% 
  mutate(L_half = ifelse(L > 1/2, 1, 0))

approximation.dt <- approximation.dt %>% 
  mutate(theo_bound = target_error * L /(1-L))

approximation.dt <- approximation.dt %>% 
  mutate(theo_success = abs(error) <= theo_bound )

approximation.dt <- approximation.dt %>%
  mutate(comparative_success = comparative_error <= 1)

approximation.dt <- approximation.dt %>%
  mutate(actual_success = as.logical(success + comparative_success))

success <- aggregate(success~method+L_half, approximation.dt, mean)
success.max <- aggregate(success_max~method+L_half, approximation.dt, mean)
comparative <- aggregate(comparative_success~method+L_half, approximation.dt, mean)
actual <- aggregate(actual_success~method+L_half, approximation.dt, mean)

success.results <- plyr::join_all(list(success, success.max,
                                       comparative, actual),
                                  by = c("method", "L_half"),
                                  type = "left",
                                  match = "all")

success.results

to.present <- success.results
to.present$success_max <- NULL
print(xtable::xtable(to.present), include.rownames = FALSE)

aggregate(success~method+target_error, approximation.dt, mean)
aggregate(comparative_success~method+target_error, approximation.dt, mean)
aggregate(actual_success~method+target_error, approximation.dt, mean)

aggregate(success~method+target_error+L_half, approximation.dt, mean)
aggregate(comparative_success~method+target_error+L_half, approximation.dt, mean)
aggregate(actual_success~method+target_error+L_half, approximation.dt, mean)

aggregate(success~method+L_half, approximation.dt, mean)
aggregate(comparative_success~method+L_half, approximation.dt, mean)
aggregate(actual_success~method+L_half, approximation.dt, mean)

aggregate(theo_success~method, approximation.dt, mean)

aggregate(theo_success~method+L_half, approximation.dt, mean)

subset(approximation.dt,
       !actual_success
       & method == "Threshold")%>% print(n = 20)

subset(approximation.dt,
       !actual_success 
       & method == "BoundingPair") %>% print(n = 10)

suspicious <- subset(approximation.dt, success & !success_max)
suspicious

subset(approximation.dt, row_id == suspicious$row_id[1])
