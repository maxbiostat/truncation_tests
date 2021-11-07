library(sumR)
library(extraDistr)
source("poisson_Betabinomial_trunc_aux.r")

spit_approx <- function(Mu, Alpha, Beta, X, Eps, M){
  Theta <- c(Mu, Alpha, Beta, X)
  
  naive <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                             parameters = Theta,
                             epsilon = Eps,
                             maxIter = M,
                             n0 = X,
                             forceAlgorithm = 1)
  
  adaptive <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                                parameters = Theta,
                                logL = log(0),
                                epsilon = Eps,
                                maxIter = M,
                                n0 = X,
                                forceAlgorithm = 2)
  results <- tibble::tibble(
    rbind(
      data.frame(
        mu = Mu, 
        alpha = Alpha,
        beta = Beta,
        x = X,
        epsilon = Eps,
        logPr = naive$sum,
        n_iter = naive$n,
        algorithm = "threshold"
      ),
      data.frame(
        mu = Mu, 
        alpha = Alpha,
        beta = Beta,
        x = X,
        epsilon = Eps,
        logPr = adaptive$sum,
        n_iter = adaptive$n,
        algorithm = "bounding_pair"
      )
    )
  )
  return(results)
}


#############

target.eps <- .Machine$double.eps
maxIter <- 3E5

spit_approx(Mu = 100, Alpha = .1, Beta = 10, X = 10,
            Eps = target.eps, M = maxIter)

approximation.grid <- expand.grid(
  Mu = c(10, 100, 1000, 10000),
  Alpha = c(0.1, 1, 10, 20),
  Beta = c(0.1, 1, 10, 20),
  X = c(0, 1, 5, 10, 100)
)

approximation.list <- lapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(Mu = linha$Mu,
                Alpha = linha$Alpha,
                Beta = linha$Beta,
                X =  linha$X,
                Eps = target.eps,
                M = maxIter)      
  })

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt

## Sum-to-threshold = STT

STT.only <- subset(approximation.dt, algorithm == "threshold")

ns <- STT.only$n_iter
STT.sorted <- STT.only[order(ns), ]
tail(STT.sorted, 20)

k1 <- 2
k2 <- 2

tab <- rbind(
  head(STT.sorted, k1),
  tail(subset(STT.sorted, mu == 1000), k1),
  tail(STT.sorted, k2)
)

tab$epsilon <- NULL
tab$x <- NULL
tab$algorithm <- NULL

print(xtable::xtable(tab), include.rownames = FALSE)
