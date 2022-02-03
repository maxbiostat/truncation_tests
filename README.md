# Test scripts for infinite series truncation algorithms

Accompanying code to the paper "Adaptive truncation of infinite sums: applications to Statistics" by [Luiz Max Carvalho](https://github.com/maxbiostat) and [Guido A. Moreira](https://github.com/GuidoAMoreira).

Note that the code provided here makes heavy use of our package, **sumR**, which can be found [here](https://github.com/GuidoAMoreira/sumR).

- The importance of selecting the right algorithm: as [this](https://github.com/maxbiostat/truncation_tests/blob/main/tests/adaptive_versus_threshold.r) script demonstrates, when L > 1/2 one really ought to use the error-bounding ("adaptive") algorithm or Batches with a suitably chosen `batch_size`. 

- Ever wondered how many iterations you would need to correctly approximate the normalising constant of the [Conway-Maxwell Poisson](https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution) distribution? [This](https://github.com/maxbiostat/truncation_tests/blob/main/COMP_normalisingConstant_table.r) script computes it for a few values, extending the results in Figure 5 of [Benson & Friel (2021)](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-3/Bayesian-Inference-Model-Selection-and-Likelihood-Estimation-using-Fast-Rejection/10.1214/20-BA1230.full) .

- [This](https://github.com/maxbiostat/truncation_tests/blob/main/MMLE_Erlang_sumR.r) script implements the marginal maximum likelihood estimation (MMLE) example with a toy queuing model. The basic idea is that you have a Poisson number of calls with exponential duration, but you only observe the total call time (i.e. their sum). A simulation study is implemented [here](https://github.com/maxbiostat/truncation_tests/blob/main/MMLE_Erlang_sumR_simuStudy.r) and has a companion [script](https://github.com/maxbiostat/truncation_tests/blob/main/analyse_MMLE_Erlang.r) to analyse it.
 
- Tests where the true answer is known in closed form can be found in [this](https://github.com/maxbiostat/truncation_tests/tree/main/tests) folder. In particular, we provide tests for the [factorial moments](https://en.wikipedia.org/wiki/Factorial_moment) of a Poisson random variable ([script](https://github.com/maxbiostat/truncation_tests/blob/main/tests/Poisson_factorial_moments.r)) and the marginal probability in a size-independent (binomial) observation error model with a negative binomial count-generating distribution ([script](https://github.com/maxbiostat/truncation_tests/blob/main/tests/Negative_binomial_obsError.r)).
