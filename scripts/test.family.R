# Compile test imputation set-up on one tree

# Load libraries
library(tidyverse)
library(MCMCglmm)
library(ape)
library(readr)

# Load data ---------------------------------------------------------------

# Load MR data
mr <- read_csv("builds/metabolic_rate_data.csv", col_types = cols())

mr <- mr %>% 
  select(-Source) %>% 
  mutate(dataset = "mr") %>% 
  as.data.frame()

# Load forest
forest <- readRDS("builds/forest.rds")

# Set options -------------------------------------------------------------
# Number of mcmc samples per (1000 trees)
# Run 333 for good chains for testing convergence
mcmc.samples <- 333

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))
thin <- 75
burnin <- 1000
nitt <- mcmc.samples * thin + burnin
tree <- forest[[1]]
inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)

set.seed(42)

mr

chain.0 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE)


m0 <- glm(log10MR ~ log10BM * MR.type, data = mr)
mr0 <- mr
mr0[c("pred", "se")] <- bind_cols(predict(m0, se.fit = TRUE)[1:2])
alpha <- 0.95
sc <- qnorm(0.975)  ## Normal approx. to likelihood
ggplot(mr0, aes(log10BM, log10MR, col = MR.type)) +
  geom_point(alpha = .3) +
  geom_line(aes(y = pred)) +
  geom_ribbon(aes(ymin = pred-se*qnorm(0.975), ymax = pred+se*qnorm(0.975), x = log10BM), alpha = .3) +
  theme_bw()

m1 <- glm(MR ~ log10BM * MR.type, data = mr, family = Gamma(link = log))
mr1 <- mr
mr1[c("pred0", "se")] <- bind_cols(predict(m1, se.fit = TRUE, type = "link")[1:2])
linkinv <- family(m1)$linkinv
mr1$pred <- linkinv(mr1$pred0)

alpha <- 0.95
sc <- qnorm(0.975)  ## Normal approx. to likelihood
ggplot(mr1, aes(log10BM, log10MR, col = MR.type)) +
  geom_point(alpha = .3) +
  geom_line(aes(y = log10(pred))) +
  geom_ribbon(aes(ymin = log10(linkinv(pred0-se*qnorm(0.975))), ymax = log10(linkinv(pred0+se*qnorm(0.975))), x = log10BM), alpha = .3) +
  theme_bw()


chain.1 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                    family = Gamma(link = log), ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE)
