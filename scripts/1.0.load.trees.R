# Loads the big forest and stores it as an R-object for faster handling
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ape)

# Load PHYLACINE 1.2.1 posterior phylogeny distribution
forest <- read.nexus("../PHYLACINE_1.2/Data/Phylogenies/Complete_phylogeny.nex")

# Store for fast loading
write_rds(forest, "builds/forest.rds")