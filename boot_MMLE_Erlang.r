load("bootstrap_mu=1500_beta=1_K=50_B=5.RData")
library(boot)

boot.ci(boot.fixed.bessel)
plot(boot.fixed.bessel, index = 1)
# plot(boot.fixed.bessel, index = 2)

boot.ci(boot.adaptive.bessel)
plot(boot.adaptive.bessel, index = 1)
# plot(boot.adaptive.bessel, index = 2)

boot.ci(boot.adaptive.full)
plot(boot.adaptive.full, index = 1)
# plot(boot.adaptive.full, index = 2)
