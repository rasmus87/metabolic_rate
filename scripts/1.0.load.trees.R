# Loads the big forest and stores it as an R-object for faster handling

require(ape)

forest <- read.nexus("../PHYLACINE_1.2/Data/Phylogenies/Complete_phylogeny.nex")

saveRDS(forest, "builds/forest.rds")
