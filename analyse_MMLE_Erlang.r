library(tidyverse)

is_in <- function(x, l, u){
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}

MM <- 1500
BB <- .1
J <- 50 

simus <- read_csv(paste0("Mu=", MM, "_Beta=", BB, "_K=", J, ".csv"))

# load(
#   paste0("Mu=", MM, "_Beta=", BB, "_K=", J, ".RData")
# )

# setdiff(1:500, unique(simus$replicate))

simus <- simus %>% mutate(squared_deviation = (point-true)^2,
                              .before = "time")

simus <- simus %>% mutate(covers = is_in(x = true,
                                         l = lwr, u = upr),
                          .after = "squared_deviation")

simus <- simus %>% mutate(covers_exact = is_in(x = true,
                                               l = lwr_exact, u = upr_exact),
                          .after = "covers")

simus <- simus %>% mutate(CI_width = upr-lwr,
                          .after = "covers_exact")

simus <- simus %>% mutate(CI_width_exact = upr_exact-lwr_exact,
                          .after = "CI_width")
##############

aggregate(time~implementation+method, simus, mean)

mus <- subset(simus, parameter == "mu")

betas <- subset(simus, parameter == "beta")

# in what follows, 'parameter' is supefluous, but kept for clarity

## Mu
aggregate(sqrt(squared_deviation)/true~implementation+method+parameter, mus, mean)
aggregate(squared_deviation~implementation+method+parameter, mus,
          function(x) sqrt(mean(x)))
aggregate(covers~implementation+method+parameter, mus, mean)
aggregate(covers_exact~implementation+method+parameter, mus, mean)

aggregate(CI_width~implementation+method+parameter, mus, mean)
aggregate(CI_width_exact~implementation+method+parameter, mus, mean)

## Beta
aggregate(sqrt(squared_deviation)/true~implementation+method+parameter,
          betas, mean)
aggregate(squared_deviation~implementation+method+parameter, betas,
          function(x) sqrt(mean(x)))
aggregate(covers~implementation+method+parameter, betas, mean)
aggregate(covers_exact~implementation+method+parameter, betas, mean)

aggregate(CI_width~implementation+method+parameter, betas, mean)
aggregate(CI_width_exact~implementation+method+parameter, betas, mean)
