library(ape)
library(caper)
library(MCMCglmm)
# library(devtools)
# install_github("TGuillerme/mulTree", ref = "release")
library(mulTree)

# Data have been log transformed, mean centered and expressed in units of
# standard deviation.
data(lifespan)
head(lifespan_volant)

# The 10Ktrees from Khun et al (2011) gives a set of trees to represent
# different polytomies. For now let just take one.
mammal_tree <- trees_mammalia[[1]]
plot(mammal_tree, cex = 0.3)

# The number of species
Ntip(mammal_tree)

# We can also check that its ultrametric
is.ultrametric(mammal_tree)

# Subset for mammals
lifespan_mammals <- lifespan_volant[lifespan_volant$class == "Mammalia",]

# lets define our fixed factors
formula_a <- longevity ~ mass + volant

#  and run a simple glm
glm_mod <- glm(formula = formula_a, family = "gaussian",
               data = lifespan_mammals)
summary(glm_mod)

ggplot(lifespan_mammals, aes(mass, longevity, col = volant)) +
  geom_point() +
  geom_smooth(method = "lm")

# Create the comparative data
### DOES NOT HANDLE MULTIPLE MEASURES PER SPECIES ###
comp_data <- comparative.data(phy = mammal_tree, data =lifespan_volant,
                              names.col = species, vcv = TRUE)
head(comp_data$data)

# Note that in the comp_data$data that there are now no birds;
# these have been dropped
head(comp_data$dropped)

# We have the formula and the comparative.data object comp_data which contains
# but the phylogeny and the data. Lets set the lambda in this case to 1. 
pgls_l1 <- pgls(formula = formula_a, data = comp_data, lambda = c(1))
pgls_l0 <- pgls(formula = formula_a, data = comp_data, lambda = c(0.01))
summary(pgls_l1)
summary(pgls_l0)

# Finally we also need to set the lambda in this case to ML. This means the we
# will using Maximum Likelihood to calculate the lambda.
pgls_mod <- pgls(formula = formula_a, data = comp_data, lambda = "ML")
summary(pgls_mod)

# How phyl autocorrelated is our data? 0 means not at all 1 means fully
pgls_mod$param["lambda"]

mod_profile <- pgls.profile(pgls_mod)
plot(mod_profile)

ggplot(lifespan_mammals, aes(mass, longevity, col = volant)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  geom_abline(slope = 0.416752, intercept = 0.076286 + c(0, 0.988252))

## MCMCglmm

# For the random effect variance prior we will use an inverse-Gamma
# distribution. In MCMCglmm this is described by two parameters: nu and V. These
# terms are related to the shape (alpha) and scale (beta) parameters on an
# inverse-Gamma with alpha = nu/2, and Beta = (nu*V)/2. As we don’t want our
# estimates to be heavily influenced by our prior we will use weakly informative
# prior values such as descripted as V = 1 and nu = 0.002

prior <- list(R = list(V=1, nu=0.002),
              G = list(G1 = list(V=1, nu=0.002)))

# To save time we will only run this model over 12000 iterations (however, much
# larger nitt is often required).

# Number of interations
nitt <- 12000

# Length of burnin
burnin <- 2000

# Amount of thinning
thin <- 5

#Matched data
mcmc_data <- comp_data$data

#As MCMCglmm requires a colume named animal for it to identify it as a phylo
# model we include an extra colume with the species names in it.
mcmc_data <- cbind(animal = rownames(mcmc_data), mcmc_data)
mcmc_tree <- comp_data$phy

mod_mcmc <- MCMCglmm(fixed = formula_a, 
                     random = ~ animal, 
                     family = "gaussian",
                     pedigree = mcmc_tree, 
                     data = mcmc_data,
                     nitt = nitt,
                     burnin = burnin,
                     thin = thin,
                     prior = prior)
# Looks good:
plot(mod_mcmc$Sol)
plot(mod_mcmc$VCV)
summary(mod_mcmc)

inv.phylo <- inverseA(mammal_tree, nodes = "ALL", scale = TRUE)
lifespan_volant <- lifespan_volant[lifespan_volant$species %in% mammal_tree$tip.label, ]
mod_mcmc_v2 <- MCMCglmm(fixed = formula_a, 
                        random = ~species,
                        family = "gaussian",
                        ginverse = list(species = inv.phylo$Ainv),
                        data = lifespan_volant,
                        nitt = nitt,
                        burnin = burnin,
                        thin = thin,
                        prior = prior)
summary(mod_mcmc_v2)

# Looks bad:
mod_mcmc_short_run <- MCMCglmm(fixed = formula_a, 
                               random= ~ animal, 
                               family="gaussian",
                               pedigree = mcmc_tree, 
                               data = mcmc_data,
                               nitt = c(1000),
                               burnin = c(1),
                               thin = c(1),
                               prior = prior,
                               verbose=FALSE)

traceplot(mod_mcmc_short_run$VCV[,2])

autocorr.diag(mod_mcmc$Sol)
autocorr.diag(mod_mcmc$VCV)

# acf plot for the first fixed estimate in our model (the intercept)
acf(mod_mcmc$Sol[,1], lag.max = 20)

# acf plot for the first random term in our model (the animal term)
acf(mod_mcmc$VCV[,1], lag.max = 20)


# Fix autocorrelation:
nitt2 <- 240000
burnin2 = 40000
thin2 = 100

mod_mcmc_long <- MCMCglmm(fixed = formula_a, 
                          random= ~ animal, 
                          family="gaussian",
                          pedigree = mcmc_tree, 
                          data = mcmc_data,
                          nitt = nitt2,
                          burnin = burnin2,
                          thin = thin2,
                          prior = prior,
                          verbose = FALSE)

acf(mod_mcmc_long$VCV[,1], lag.max = 20)
effectiveSize(mod_mcmc_long$Sol)
effectiveSize(mod_mcmc_long$VCV)

mod_mcmc_2 <- MCMCglmm(fixed = formula_a, 
                       random = ~ animal, 
                       family ="gaussian",
                       pedigree = mcmc_tree, 
                       data = mcmc_data,
                       nitt = nitt2,
                       burnin = burnin2,
                       thin = thin2,
                       prior = prior,
                       verbose = FALSE)

# We can now check the convergence of the two chains using the Gelman and Rubin
# Multiple Sequence Diagnostic. This calculates the within-chain and
# between-chain variance of the chains and then gives a scale reduced factor,
# (see here for more). When this number is close to one (say below 1.1) the
# chains are indistinguishable and hence can be considered to be converged.

# Checking convergence for our fixed factors
gelman.diag(mcmc.list(mod_mcmc_long$Sol, mod_mcmc_2$Sol))

# Checking convergence for our random terms
gelman.diag(mcmc.list(mod_mcmc_long$VCV, mod_mcmc_2$VCV))

summary(mod_mcmc_long)
median(mod_mcmc_long$Sol[,1])

H <- ( var(mod_mcmc_long$VCV[,"animal"]))/
  (var(mod_mcmc_long$VCV[,"animal"]) + var(mod_mcmc_long$VCV[,"units"]))
H


data(lifespan)
# Data have been log transformed, mean centered and expressed in units of
# standard deviation.
head(lifespan_volant)

# Subset for mammals
lifespan_mammals <- lifespan_volant[lifespan_volant$class == "Mammalia",]



##lets package all the data up into one mulTree object 
mulTree_data <- as.mulTree(data = lifespan_mammals, tree = trees_mammalia,
                           taxa = "species")

is.ultrametric(trees_mammalia)

# The formula
mul_formula <- longevity ~ mass + volant
# The MCMC parameters (iterations, thining, burnin)
mul_parameters <- c(nitt2, thin2, burnin2)
# The MCMCglmm priors
mul_priors <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002)))

getwd()



mulTree(mulTree.data = mulTree_data, formula = mul_formula, priors = mul_priors,
        parameters = mul_parameters, output = "longevity_example", ESS = 1000,
        chains = 2)


# Reading only one specific model
one_model <- read.mulTree("longevity_example-tree1_chain1", model = TRUE)

# This model is a normal MCMCglmm object that has been ran on one single tree
class(one_model) ; names(one_model)

# Reading the convergence diagnosis test to see if the two chains converged for
# each tree
read.mulTree("longevity_example", convergence = TRUE)
# As indicated here, the chains converged for both chains!

# Reading all the models to perform the MCMCglmm analysis on multiple trees
all_models <- read.mulTree("longevity_example")
str(all_models)
# This object contains 39600 estimations of the Intercept and the terms!

# If you want to remove the chains from the current directory run the following:
# file.remove(list.files(pattern="longevity_example"))
# However when doing your actual analysis you should keep all your models stored
# somewhere!

summarised_results <- summary(all_models, use.hdr = FALSE, cent.tend = mean,
                              prob = c(75, 25))

plot(summarised_results, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-2,2), cex.terms = 0.5,
     terms = c("Intercept", "Body Mass", "Volancy", "Phylogeny", "Residuals"),
     col = "grey", cex.main = 0.8)
