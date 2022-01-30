# Test scripts for infinite series truncation algorithms
Preliminary code to test truncation algorithms for infinite series and their applications in Statistics.

- The importance of selecting the right algorithm: as [this](https://github.com/maxbiostat/truncation_tests/blob/main/tests/adaptive_versus_threshold.r) script demonstrates, when L > 1/2 one really ought to use the error-bounding ("adaptive") algorithm.

- Ever wondered how many iterations you would need to correctly approximate the normalising constant of the [Conway-Maxwell Poisson](https://en.wikipedia.org/wiki/Conway%E2%80%93Maxwell%E2%80%93Poisson_distribution) distribution? [This](https://github.com/maxbiostat/truncation_tests/blob/main/COMP_normalisingConstant_table.r) script computes it for a few values, extending the results in Figure 5 of [Benson & Friel (2021)](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-3/Bayesian-Inference-Model-Selection-and-Likelihood-Estimation-using-Fast-Rejection/10.1214/20-BA1230.full) .

- [This](https://github.com/maxbiostat/truncation_tests/blob/main/MMLE_Erlang_sumR.r) script implements the marginal maximum likelihood estimation (MMLE) example with a toy queuing model. The basic idea is that you have a Poisson number of calls with exponential duration, but you only observe the total call time (i.e. their sum). A simulation study is implemented [here](https://github.com/maxbiostat/truncation_tests/blob/main/MMLE_Erlang_sumR_simuStudy.r) and has a companion [script](https://github.com/maxbiostat/truncation_tests/blob/main/analyse_MMLE_Erlang.r) to analyse it.
 
