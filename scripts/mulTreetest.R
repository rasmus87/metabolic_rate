library(mulTree)
data(lifespan)

combined_trees <- tree.bind(trees_mammalia, trees_aves, sample = 2,
                            root.age = 250)
is.ultrametric(combined_trees)

sessionInfo()

data(lifespan)
# Data have been log transformed, mean centered and expressed in units of
# standard deviation.
head(lifespan_volant)

##lets package all the data up into one mulTree object 
mulTree_data <- as.mulTree(data = lifespan_volant, tree = combined_trees,
                           taxa = "species")


# The formula
mul_formula <- longevity ~ mass + volant
# The MCMC parameters (iterations, thining, burnin)
mul_parameters <- c(nitt, thin, burnin)
# The MCMCglmm priors
mul_priors <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002)))


mulTree(mulTree.data = mulTree_data, formula = mul_formula, priors = mul_priors,
        parameters = mul_parameters, output = "longevity_example", ESS = 1000,
        chains = 2)
