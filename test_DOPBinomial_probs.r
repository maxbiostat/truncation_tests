source("double_poisson_binomial_aux.r")
Mu <- 4
Phi <- 1/2
PP <- .5
J <- 5000
simu <- simulate_obsdata(n = J, mu = Mu, phi = Phi, p = PP)
data <- compress_counts(simu$obs_x)
data

lps <- sapply(data$count,
              function(xx) marg_lik(x = xx, pars = c(Mu, Phi, PP)))
  
pr_pois <- sapply(data$count,
                     function(xx) dpois(x = xx, lambda = Mu*PP))

prob.tab <- tibble::tibble(
  x = data$count,
  pr_x = data$freq/sum(data$freq),
  theo_pr_x = exp(lps),
  theo_pr_x_pois = pr_pois
)
prob.tab
tail(prob.tab)
