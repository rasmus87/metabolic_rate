library(tidyverse)
library(gridExtra)
library(coda)
library(MCMCglmm)
library(ggpmisc)

i = 1
chain.1 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain1.rds"))
chain.2 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain2.rds"))
chain.3 <- readRDS(paste0("builds/mcmcglmms/tree", i, ".chain3.rds"))

### Checking only the 3 chains
sol <- bind_rows(as.data.frame(chain.1$Sol[, 1:4]), 
                 as.data.frame(chain.2$Sol[, 1:4]), 
                 as.data.frame(chain.3$Sol[, 1:4]))
sol["chain"] <- gl(3, 333)
sol["sample"] <- rep(1:333, 3)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) +
  geom_line() + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, lty = "dotted", col = "black") +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme_bw() +
  theme(legend.position="none") + 
  ylab("") + 
  stat_poly_eq(aes(label = ..adj.rr.label..), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3, vstep = 0, hstep = 0.20)
right <- ggplot(sol, aes(x = value, col = chain)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x = "", y = "")
p.main <- grid.arrange(left, right, nrow = 1)
ggsave("output/appendix2_figX1.png", p.main, width = 25.6, height = 28.8, units = "cm")

VCV <- bind_rows(as.data.frame(chain.1$VCV), 
                 as.data.frame(chain.2$VCV), 
                 as.data.frame(chain.3$VCV))
VCV["chain"] <- gl(3, 333)
VCV["sample"] <- rep(1:333, 3)
VCV <- gather(VCV, key = "variable", value = "value", -chain, -sample)

left <- ggplot(VCV, aes(x = sample, y = value, col = chain)) +
  geom_line() + 
  geom_smooth(formula = y ~ x, method = "lm", se = TRUE, lty = "dotted", col = "black") +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  ylab("") + 
  stat_poly_eq(aes(label = ..adj.rr.label..), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3, vstep = 0, hstep = 0.20)
right <- ggplot(VCV, aes(x = value, col = chain)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  theme_bw() +
  theme(legend.position="none") + 
  labs(x = "", y = "")
p.random <- grid.arrange(left, right, nrow = 1)
ggsave("output/appendix2_figX2.png", p.random, width = 25.6, height = 14.4, units = "cm")


# Checking convergence for our fixed factors
gelman.diag(mcmc.list(chain.1$Sol[, 1:4], chain.2$Sol[, 1:4], chain.3$Sol[, 1:4]), autoburnin = FALSE)

# Checking convergence for our random terms
gelman.diag(mcmc.list(chain.1$VCV, chain.2$VCV, chain.3$VCV), autoburnin = FALSE)


# Gelman plot:
gelman.plot(mcmc.list(chain.1$VCV, chain.2$VCV, chain.3$VCV), autoburnin = FALSE)
gelman.plot(mcmc.list(chain.1$Sol[, 1:2], chain.2$Sol[, 1:2], chain.3$Sol[, 1:2]), autoburnin = FALSE)


### Checking effective sample size
chain.1.2.3.Sol <- list(chain.1$Sol[, 1:2], chain.2$Sol[, 1:2], chain.3$Sol[, 1:2])
chain.1.2.3.VCV <- list(chain.1$VCV, chain.2$VCV, chain.3$VCV)
effectiveSize(chain.1.2.3.Sol)/3
# G/R structure
effectiveSize(chain.1.2.3.VCV)/3

# acf plot for the fixed estimates
acf(chain.1$Sol[, 1], lag.max = 20)
acf(chain.1$Sol[, 2], lag.max = 20)
acf(chain.1$Sol[, 3], lag.max = 20)
acf(chain.1$Sol[, 4], lag.max = 20)
# acf plot for the first random term in our model (the phyl term)
acf(chain.1$VCV[, 1], lag.max = 20)
acf(chain.1$VCV[, 2], lag.max = 20)

lambda <- chain.1$VCV[,'Binomial.1.2'] / (chain.1$VCV[,'Binomial.1.2'] + chain.1$VCV[,'units'])
mean(lambda)

posterior.mode(lambda)
HPDinterval(lambda)

summary(chain.1)

n <- nrow(mr)
sigma1 <- sqrt(sum((predict(chain.1) - mr$log10MR)^2)/(n-2))
sigma1
sigma2 <- sqrt(sum((predict(chain.2) - mr$log10MR)^2)/(n-2))
sigma2
sigma3 <- sqrt(sum((predict(chain.3) - mr$log10MR)^2)/(n-2))
sigma3
mean(c(sigma1,sigma2,sigma3))
