# field.metabolic.rate.imputation
# Based on the package Rphylopars

## Load required libraries
library(tidyverse)
library(Rphylopars)
library(tictoc)
library(doSNOW)

## Set options:
# Set parralell cluster size
cluster.size <- 2
# How many trees do you want to run this for? 2-1000?
n.trees <- 2


## Load data
fmr <- read_csv("builds/fmr.csv", col_types = cols())
bmr <- read_csv("builds/bmr.csv", col_types = cols())

mr <- bind_rows(fmr, bmr)
mr <- mr %>% transmute(species = Binomial.1.2,
                       log10BM,
                       log10BMR,
                       log10FMR)
# Load 1000 trees
forest <- read.nexus("../PHYLACINE_1.1/Data/Phylogenies/Complete_phylogeny.nex")

# Load trait dataset
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
terrestrial <- mam %>% filter(Terrestrial == 1) %>% pull(Binomial.1.2)
mam <- mam %>% filter(Binomial.1.2 %in% terrestrial)

# Make data ready
df <- mam %>% transmute(species = Binomial.1.2,
                        log10BM = log10(Mass.g))
df <- bind_rows(df, mr)
df <- df %>% filter(species %in% terrestrial)
df <- as.data.frame(df)

tree <- forest[[1]]
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% terrestrial])
res <- phylopars(df, tree)

anc.rec <- mam %>% left_join(as_data_frame(res$anc_recon[1:nrow(df), ],
                                           rownames = "Binomial.1.2"),
                             by = "Binomial.1.2")
test <- anc.rec %>% transmute(Binomial.1.2, m0=Mass.g, m1 = 10^log10BM, diff = m0-m1,
                              diff2 = diff/m0*100)
test %>% ggplot(aes(diff)) + geom_density() + geom_rug() + scale_x_log10()
test %>% ggplot(aes(diff2)) + geom_density() + geom_rug()
test %>% filter(abs(diff2) > 100)

res.df <- as.data.frame(res$anc_recon[1:length(tree$tip.label), ])
res.df$Binomial.1.2 <- tree$tip.label
res.df <- left_join(res.df, mam %>% select(Binomial.1.2, Order.1.2, Family.1.2))
res.df.var <- as.data.frame(res$anc_var[1:length(tree$tip.label), ])
names(res.df.var) <- paste0(names(res.df.var), ".var")
res.df <- bind_cols(res.df, res.df.var)


ggplot(res.df, aes(log10BM, log10FMR)) +
  geom_point(data = fmr, col = "grey") + 
  geom_errorbar(aes(ymin = log10FMR - log10FMR.var,
                    ymax = log10FMR + log10FMR.var,
                    col = Order.1.2)) +
  geom_errorbarh(aes(xmin = log10BM - log10BM.var,
                     xmax = log10BM + log10BM.var,
                     col = Order.1.2)) +
  geom_point(aes(col = Order.1.2)) +
  theme_bw()


hc <- readxl::read_xlsx("data/halls cave clean.xlsx", sheet = "hall cave clean")
hc <- hc %>% mutate(Binomial.1.2 = paste(genus.new, 
                                         ifelse(is.na(likely.sp), species.new, likely.sp),
                                         sep = "_"))
hc <- left_join(hc, mam %>% dplyr::select(Binomial.1.2:Family.1.2))


## Impute mising data
n.rows <- nrow(df)
cl <- makeCluster(cluster.size)
registerDoSNOW(cl)

timestamp()
tic()
pb <- txtProgressBar(max = n.trees, style = 3)

progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
imputed <- foreach(i = 1:n.trees, 
                   .packages=c('Rphylopars'), 
                   .combine = rbind,
                   .inorder = FALSE,
                   .options.snow = opts) %dopar% 
                   {
                     p_BM <- phylopars(trait_data = mam.traits, tree = forest[[i]])
                     res <- p_BM$anc_recon[1:n.rows, ]
                     res <- as.data.frame(res)
                     res$species <- rownames(res)
                     rownames(res) <- NULL
                     return(res)
                   }

toc()

stopCluster(cl)
gc()