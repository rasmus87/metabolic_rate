library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)

## Set options:
# Set parralell cluster size
# cluster.size <- 2
cluster.size <- 6
#cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
# n.trees <- 1
# n.trees <- 6
# n.trees <- 18*10
n.trees <- 1000

# Number of mcmc samples per (1000 trees)
# Run 333 for good chains for testing convergence
# Run 3 samples for actual data is enough
mcmc.samples <- 3

mr <- read_csv("builds/mr.csv", col_types = cols()) # deprecated
mr <- read_csv("builds/metabolic_rate_data.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

mr <- mr %>% 
  select(-Source) %>% 
  mutate(dataset = "mr")
mr$Binomial.1.2 %>% unique %>% length
mr$Family.1.2 %>% unique %>% length
mr$Order.1.2 %>% unique %>% length
mr %>% count(MR)
mr %>% filter(MR.type == "BMR") %>% pull(Binomial.1.2) %>% unique() %>% length()
mr %>% filter(MR.type == "FMR") %>% pull(Binomial.1.2) %>% unique() %>% length()
cut(10^mr$log10BM/1000, breaks = c(0,1,10,100,1000,10000)) %>% table(useNA = "a")
mr <- as.data.frame(mr) # MR dataset for imputation

# Select all species we want prediction for
mam <- mam %>% 
  mutate(log10BM = log10(Mass.g), MR.type = "BMR", log10MR = NA)
mam <- mam %>% bind_rows(mutate(mam, MR.type = "FMR"))
mam <- mam %>% mutate(dataset = "mam")
mam <- mam %>% select(c("Binomial.1.2",
                        "Order.1.2",
                        "Family.1.2",
                        "MR.type",
                        "log10BM",
                        "log10MR",
                        "MR.type",
                        "dataset"))

# Combine with imputation dataset for prediction
n.mam <- nrow(mam)
df <- bind_rows(mam, mr)
df <- as.data.frame(df)

forest <- readRDS("builds/forest.rds")

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))
thin <- 75
burnin <- 1000
nitt <- mcmc.samples * thin + burnin
mcmc.regression <- function(i) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = mr, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  chain.2 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = mr, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  chain.3 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = mr, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE)
  gc()

  pred1 <- MCMC.predict(chain.1, df)
  pred2 <- MCMC.predict(chain.2, df)
  pred3 <- MCMC.predict(chain.3, df)
  
  post.pred <- rbind(pred1, pred2, pred3)
  
  solution <- rbind(chain.1$Sol[, 1:4],
                    chain.2$Sol[, 1:4],
                    chain.3$Sol[, 1:4])
  solution <- as.data.frame(solution)
  random.effects <- 5:ncol(chain.1$Sol)
  solution["random.effect"] <- c(rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.1$Sol[, random.effects]),
                                 rowMeans(chain.1$Sol[, random.effects]))
  solution["tree"] <- i
  solution["chain"] <- rep(1:3, each = mcmc.samples)
  
  return(list(solution, post.pred))
}

MCMC.predict <- function(object, newdata) {
  object2 <- MCMCglmm(fixed=object$Fixed$formula, 
                      random=object$Random$formula, 
                      rcov=object$Residual$formula, 
                      family=object$Residual$original.family,
                      data=newdata, 
                      nitt=1, 
                      thin=1,
                      burnin=0, 
                      ginverse=object$ginverse, 
                      verbose=FALSE, 
                      pr=TRUE)
  
  W <- cbind(object2$X, object2$Z)
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:n.mam]
  
  colnames(post.pred) <- mam$Binomial.1.2
  
  return(post.pred)
}

comb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}

cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
timestamp()
tic()
pb <- txtProgressBar(max = n.trees, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
imputed <- foreach(i = 1:n.trees,
                   .packages = c('MCMCglmm'),
                   .inorder = FALSE,
                   .options.snow = opts,
                   .combine = comb,
                   .multicombine = TRUE) %dopar% mcmc.regression(i)
toc()
stopCluster(cl)
gc()

write_csv(as_tibble(imputed[[1]]), paste0("builds/", mcmc.samples ,"_mr_fit.solution.csv"))
write_csv(as_tibble(imputed[[2]]), paste0("builds/", mcmc.samples ,"_mr_post.pred.csv"))