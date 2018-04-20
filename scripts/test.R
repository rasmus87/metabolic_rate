library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

## Set options:
# Set parralell cluster size
cluster.size <- 6
#cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
n.trees <- 6
#n.trees <- 1000

mr <- read_csv("builds/mr.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

terrestrial <- mam %>% filter(Terrestrial == 1) %>% pull(Binomial.1.2)

# Linear model:
mr <- mr %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  select(-source) %>% 
  mutate(dataset = "mr")
mr <- as.data.frame(mr)

mam <- mam %>% 
  filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(log10BM = log10(Mass.g), MR = "BMR", log10MR = NA)
mam <- mam %>% bind_rows(mutate(mam, MR = "FMR"))
mam <- mam %>% mutate(dataset = "mam")
mam <- mam %>% select(names(mr))

n.mam <- nrow(mam)

df <- rbind(mam, mr)
df <- as.data.frame(df)

forest <- readRDS("builds/forest.rds")
species <- forest[[1]]$tip.label
drop.species <- species[!species %in% terrestrial]
forest <- lapply(forest, drop.tip, tip = drop.species)

prior <- list(G = list(G1 = list(V = 1, nu = 0.02)), 
              R = list(V = 1, nu = 0.02))
thin <- 75
burnin <- thin * 10
nitt <- 333 * thin + burnin
i = 1
mcmc.regression <- function(i) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
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
  
  gc()
  pred1 <- predict(chain.1, df, marginal = NULL, interval = "confidence")[1:n.mam, ]
  pred2 <- predict(chain.2, df, marginal = NULL, interval = "confidence")[1:n.mam, ]
  pred3 <- predict(chain.3, df, marginal = NULL, interval = "confidence")[1:n.mam, ]

  pred1 <- as.data.frame(pred1)
  pred2 <- as.data.frame(pred2)
  pred3 <- as.data.frame(pred3)
  
  pred1["Binomial.1.2"] <- mam$Binomial.1.2
  pred2["Binomial.1.2"] <- mam$Binomial.1.2
  pred3["Binomial.1.2"] <- mam$Binomial.1.2
  
  pred1["Binomial.1.2"] <- mam$MR
  pred2["Binomial.1.2"] <- mam$MR
  pred3["Binomial.1.2"] <- mam$MR
  
  pred <- rbind(pred1, pred2, pred3)
  pred["tree"] <- i
  return(pred)
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
timestamp()
tic()
pb <- txtProgressBar(max = 3, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
imputed <- foreach(i = 1:n.trees,
                   .packages = c('MCMCglmm'),
                   .inorder = FALSE,
                   .options.snow = opts,
                   .combine = rbind) %dopar% mcmc.regression(i)
toc()
stopCluster(cl)
gc()
