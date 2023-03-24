# Run imputation for all trees

# Load libraries
library(tidyverse)
library(MCMCglmm)
library(ape)
library(doSNOW)
library(gridExtra)
library(tictoc)


# Load data ---------------------------------------------------------------

# Load MR data
mr <- read_csv("builds/metabolic_rate_data.csv", col_types = cols())

mr <- mr %>% 
  select("Binomial.1.2", "Order.1.2", "Family.1.2", "MR.type", "log10BM", "log10MR") %>% 
  mutate(dataset = "mr") %>% 
  as.data.frame(mr)

# Load PHYLACINE
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Add needed columns for prediction
mam2 <- mam %>% 
  mutate(log10BM = log10(Mass.g), 
         MR.type = "BMR", 
         log10MR = NA,
         dataset = "mam")
# Duplicate dataset for predictions on FMR too
mam2 <- mam2 %>% 
  bind_rows(mutate(mam2, MR.type = "FMR"))

# Select needed columns
mam2 <- mam2 %>% 
  select(c("Binomial.1.2",
           "Order.1.2",
           "Family.1.2",
           "MR.type",
           "log10BM",
           "log10MR",
           "MR.type",
           "dataset"))

# Combine with imputation dataset for prediction
n.mam2 <- nrow(mam2)
df <- bind_rows(mam2, mr)
df <- as.data.frame(df)

# Cross validation: 5-fold
# Set seed and make a random ordered data vector
set.seed(42)
n <- nrow(mr)
mr <- mr[sample(n), ]
mr <- mr %>% distinct(Binomial.1.2, MR.type, .keep_all = TRUE)
n <- nrow(mr)

mr %>% count(MR.type)
mr %>% distinct(Binomial.1.2) %>% nrow

# Create 10 equally size folds
folds <- cut(1:n, breaks = 5, labels = FALSE)

# Load forest
forest <- readRDS("builds/forest.rds")


# Set options -------------------------------------------------------------

# Set parallel cluster size
cluster.size <- 2
# cluster.size <- 6
#cluster.size <- 20
# How many trees do you want to run this for? 2-1000?
# n.trees <- 1
n.trees <- 2
# n.trees <- 18*10
# n.trees <- 1000

# Number of mcmc samples per (1000 trees)
mcmc.samples <- 1

prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))
thin <- 75
burnin <- 1000 * 2
nitt <- mcmc.samples * thin + burnin


mcmc.regression <- function(i, train.data, test.data) {
  tree <- forest[[i]]
  inv.phylo <- inverseA(tree, nodes = "ALL", scale = TRUE)
  chain.1 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.2 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  chain.3 <- MCMCglmm(log10MR ~ log10BM * MR.type, random = ~Binomial.1.2,
                      family = "gaussian", ginverse = list(Binomial.1.2 = inv.phylo$Ainv), 
                      prior = prior,
                      data = train.data, nitt = nitt, burnin = burnin, thin = thin,
                      pr = TRUE,
                      verbose = FALSE)
  gc()
  
  pred1 <- MCMC.predict(chain.1, test.data, train.data)
  pred2 <- MCMC.predict(chain.2, test.data, train.data)
  pred3 <- MCMC.predict(chain.3, test.data, train.data)
  
  post.pred <- rbind(pred1, pred2, pred3)
  
  return(post.pred)
}

MCMC.predict <- function(object, test.data, train.data) {
  newdata <- rbind(test.data, train.data)
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
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:nrow(test.data)]
  
  names(post.pred) <- test.data$Binomial.1.2
  
  return(post.pred)
}


cl <- makeCluster(cluster.size)
registerDoSNOW(cl)
timestamp()
tic()

# Perform 5 fold cross validation
res <- tibble()
for(fold in 1:5) {
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds == fold)
  
  train.data <- mr[-testIndexes, ]
  
  # Prediction dataset
  test.data <- mr[testIndexes, ] %>% 
    mutate(log10MR = NA, dataset = "mam")
  
  # Impute
  pb <- txtProgressBar(max = n.trees, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  imputed <- foreach(i = 1:n.trees,
                     .packages = c('MCMCglmm', 'tibble'),
                     .inorder = FALSE,
                     .options.snow = opts,
                     .combine = rbind,
                     .multicombine = TRUE) %dopar% mcmc.regression(i, train.data, test.data)
  imputed.mean <- colMeans(imputed)
  imputed.long <- tibble(Binomial.1.2 = names(imputed.mean),
                         log10.mr.mean = imputed.mean,
                         fold = fold,
                         MR.type = test.data$MR.type,
                         log10.mr.empirical = mr$log10MR[testIndexes])
  res <- rbind(res, imputed.long)
}

toc()
stopCluster(cl)
gc()


write_csv(res, "builds/mr_5xcross_val.csv")
res <- read_csv("builds/mr_5xcross_val.csv")


# Persons R-squared and RMSE
res %>% 
  group_by(MR.type, fold) %>% 
  summarise(persons.r2 = cor(log10.mr.mean, log10.mr.empirical)^2,
            rmse = sqrt(mean((log10.mr.empirical - log10.mr.mean)^2))) %>% 
  summarise_at(c("persons.r2", "rmse"), mean)

res %>% 
  group_by(MR.type, fold) %>% 
  summarise(persons.r2 = cor(10^log10.mr.mean, 10^log10.mr.empirical)^2,
            rmse = sqrt(mean((10^log10.mr.empirical - 10^log10.mr.mean)^2))) %>% 
  summarise_at(c("persons.r2", "rmse"), mean)

# total.res <- res %>%
#   mutate(error = abs(log10.mr.empirical - log10.mr.mean)^2,
#          error10 = abs(10^log10.mr.empirical - 10^log10.mr.mean)^2)
# ggplot(total.res, aes(log10.mr.mean, log10.mr.empirical, col = error)) +
#   geom_point() +
#   geom_smooth() +
#   geom_abline(slope = 1) +
#   scale_color_viridis_c()

