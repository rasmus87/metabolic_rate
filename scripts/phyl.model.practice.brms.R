# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html

library(tidyverse)
library(ape)
library(brms)

# Load phylo
forest <- readRDS("builds/forest.rds")
phylo <- forest[[444]]

# Load dataset
mr <- read_csv("builds/mr.csv")
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
mam <- mam %>% filter(Terrestrial == 1, Marine == 0, Freshwater == 0, Aerial == 0) 
terrestrial <- mam %>% pull(Binomial.1.2)
mr <- mr %>% filter(Binomial.1.2 %in% terrestrial)

non.terrestrial <- phylo$tip.label[!phylo$tip.label %in% terrestrial]

phylo <- drop.tip(phylo, non.terrestrial)

# Prepare phylo
library(tictoc)
tic()
inv.phylo <- MCMCglmm::inverseA(phylo, nodes = "TIPS", scale = TRUE)
toc()

tic()
A <- Matrix::solve(inv.phylo$Ainv)
toc()
rownames(A) <- rownames(inv.phylo$Ainv)

tic()
model_simple <- brm(log10MR ~ log10BM * MR * (1|Binomial.1.2), data = mr, 
                    family = gaussian(), cov_ranef = list(Binomial.1.2 = A),
                    chains = 3, iter = 1000, warmup = floor(1000/6), cores = 3)
toc()
model_simple
plot(model_simple)

predict(model_simple) %>% head()

plot(marginal_effects(model_simple), ask = FALSE)

WAIC(model_simple)
pp_check(model_simple)

bmr <- mam %>% transmute(Binomial.1.2, log10BM = log10(Mass.g), MR = "BMR")
fmr <- mam %>% transmute(Binomial.1.2, log10BM = log10(Mass.g), MR = "FMR")
new.data <- bind_rows(bmr, fmr)

predicted <- bind_cols(new.data, predict(model_simple, newdata = new.data))

#library(phytools)

set.seed(1)
tree <- pbtree(n=1000, scale = 100)
plot(tree)
inv.phylo <- MCMCglmm::inverseA(tree, nodes = "ALL", scale = TRUE)
tic()
A1 <- base::solve(inv.phylo$Ainv)
toc()
tic()
A2 <- Matrix::solve(inv.phylo$Ainv)
toc()

all.equal(as.numeric(A1), as.numeric(A2))
