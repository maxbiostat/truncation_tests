library(tidyverse)
library(ggplot2)

results <- readr::read_csv("results/PBB_MMLE_A.csv")


results
true.theta <- results$true[1:4] 

ggplot() + 
  geom_pointrange(data = subset(results, parameter != "n_iter"),
                  mapping = aes(x = replicate, y = point,
                                ymin = lwr, ymax = upr)) +
  # geom_hline(yintercept = I.quad$value,
  #            linetype = "dotted") + 
  scale_x_continuous("") +
  scale_y_continuous("") +
  facet_grid(parameter~method, labeller = label_parsed, scales = "free_y") +
  theme_bw(base_size = 20)

results

aggregate((point-true)/true~method+parameter, mean, data = results)
aggregate(point~method+parameter, var, data = results)
aggregate(covers~method+parameter, mean, data = results)

hist(subset(results, parameter == "n_iter" & method == "adaptive")$point)
quantile(subset(results, parameter == "n_iter"& method == "adaptive")$point,
         probs = c(.9,.99))

res.conv <- subset(results, point < 2E4)
aggregate((point-true)/true~method+parameter, mean, data = res.conv)
aggregate(point~method+parameter, var, data = res.conv)
aggregate(covers~method+parameter, mean, data = res.conv)

mean(subset(results, parameter == "n_iter"& method == "adaptive")$point < 1000)
