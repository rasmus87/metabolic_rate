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
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  ylab("")
right <- ggplot(sol, aes(x = value, col = chain)) +
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "")
grid.arrange(left, right, nrow = 1)

# Checking convergence for our fixed factors
gelman.diag(mcmc.list(chain.1$Sol[, 1:4], chain.2$Sol[, 1:4], chain.3$Sol[, 1:4]), autoburnin = FALSE)

# Checking convergence for our random terms
gelman.diag(mcmc.list(chain.1$VCV, chain.2$VCV, chain.3$VCV), autoburnin = FALSE)

### Checking only the first chain
effectiveSize(chain.1$Sol[, 1:4])
effectiveSize(chain.1$VCV)

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
