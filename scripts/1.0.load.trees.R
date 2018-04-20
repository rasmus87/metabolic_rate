# Loads the big forest and stores it as an R-object for faster handling

require(ape)

forest <- read.nexus("../PHYLACINE_1.1/Data/Phylogenies/Complete_phylogeny.nex")

saveRDS(forest, "builds/forest.rds")