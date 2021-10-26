library(tidyverse)
library(MCMCglmm)
library(ape)
library(gridExtra)
library(tictoc)

## Set options:
# Number of mcmc samples per (1000 trees)
# Run 333 for good chains for testing convergence
# Run 3 samples for actual data is enough
mcmc.samples <- 333

mr <- read_csv("builds/mr.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

mr <- mr %>% 
  select(-source) %>% 
  mutate(dataset = "mr")
mr$Binomial.1.2 %>% unique %>% length
mr$Family.1.2 %>% unique %>% length
mr$Order.1.2 %>% unique %>% length
mr %>% count(MR)
mr %>% filter(MR == "BMR") %>% pull(Binomial.1.2) %>% unique() %>% length()
mr %>% filter(MR == "FMR") %>% pull(Binomial.1.2) %>% unique() %>% length()
cut(10^mr$log10BM/1000, breaks = c(0,1,10,100,1000,10000)) %>% table(useNA = "a")
mr <- as.data.frame(mr) # MR dataset for imputation

# Select all species we want prediction for
mam <- mam %>% 
  mutate(log10BM = log10(Mass.g), MR = "BMR", log10MR = NA)
mam <- mam %>% bind_rows(mutate(mam, MR = "FMR"))
mam <- mam %>% mutate(dataset = "mam")
mam <- mam %>% select(names(mr))

# Combine with imputation dataset for prediction
n.mam <- nrow(mam)
df <- rbind(mam, mr)
df <- as.data.frame(df)

forest <- readRDS("builds/forest.rds")

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))
thin <- 75
burnin <- 1000
nitt <- mcmc.samples * thin + burnin
tree <- forest[[1]]
inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)

set.seed(42)

chain.1 <- MCMCglmm(log10MR ~ log10BM * MR, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE)
chain.2 <- MCMCglmm(log10MR ~ log10BM * MR, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE)
chain.3 <- MCMCglmm(log10MR ~ log10BM * MR, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE)

saveRDS(chain.1, paste0("builds/mcmcglmms/tree", i, ".chain1.rds"), compress = FALSE)
saveRDS(chain.2, paste0("builds/mcmcglmms/tree", i, ".chain2.rds"), compress = FALSE)
saveRDS(chain.3, paste0("builds/mcmcglmms/tree", i, ".chain3.rds"), compress = FALSE)
