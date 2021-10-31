# remove.packages("sumR")
# remotes::install_github("GuidoAMoreira/sumR")
library(sumR)
### Here we will make a table to match Figure 5 in Benson & Friel (2021)
######################################
COMP_lpmf_2 <- function(k, theta){
  mu <- theta[1]
  nu <- theta[2]
  return(
    nu* (k * log(mu) - lfactorial(k))
  )
}

spit_approx <- function(Mu, Nu, Eps, M){
  Theta <- c(Mu, Nu)
  
  naive <- sumR::infiniteSum(logFunction = COMP_lpmf_2,
                             parameters = Theta,
                             epsilon = Eps,
                             maxIter = M,
                             n0 = 0,
                             forceAlgorithm = 1)
  
  adaptive <- sumR::infiniteSum(logFunction = COMP_lpmf_2,
                                parameters = Theta,
                                logL = log(0),
                                epsilon = Eps,
                                maxIter = M,
                                n0 = 0,
                                forceAlgorithm = 2)
  results <- tibble::tibble(
    rbind(
      data.frame(
        mu = Mu, 
        nu = Nu,
        epsilon = Eps,
        logZ = naive$sum,
        n_iter = naive$n,
        algorithm = "threshold"
      ),
      data.frame(
        mu = Mu, 
        nu = Nu,
        epsilon = Eps,
        logZ = adaptive$sum,
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

approximation.grid <- data.frame(
  Mu = c(10, 100, 1000, 10000),
  Nu = c(0.1, 0.01, 0.001, 0.0001)
)

approximation.list <- lapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(Mu = linha$Mu,
                Nu = linha$Nu,
                Eps = target.eps,
                M = maxIter)      
  })

approximation.dt <- do.call(rbind, approximation.list)

approximation.dt
