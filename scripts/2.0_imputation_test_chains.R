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
# Run 333 for good chains for testing convergence
mcmc.samples <- 1000

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))
thin <- 75
burnin <- 1000 * 2
nitt <- mcmc.samples * thin + burnin
tree <- forest[[1]]
inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)

set.seed(42)

chain.1 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)
chain.2 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)
chain.3 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                    family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                    prior = prior,
                    data = mr, nitt = nitt, burnin = burnin, thin = thin,
                    pr = TRUE,
                    verbose = FALSE)

write_rds(chain.1, paste0("builds/mcmcglmms/tree1.chain1.rds"))
write_rds(chain.2, paste0("builds/mcmcglmms/tree1.chain2.rds"))
write_rds(chain.3, paste0("builds/mcmcglmms/tree1.chain3.rds"))
