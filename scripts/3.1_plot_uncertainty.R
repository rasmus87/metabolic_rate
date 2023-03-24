### Show family uncertainty ###


# Setup -------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ape)
library(ggtree)
# BiocManager::install("ggtreeExtra")
library(ggtreeExtra)

# Load forest
forest <- read_rds("builds/forest.rds")

# Subset to tree one
tree <- forest[[1]]

# Load PHYLACINE
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Load metabolic rate
mr <- read_csv("builds/Imputed metabolic rate.csv")


# Create family tree ------------------------------------------------------

# Replace all species names in tree with family names
index <- match(tree$tip.label, mam$Binomial.1.2)
tree$tip.label[] <- mam$Family.1.2[index]

# Subset to one tip per family
family.tree <- drop.tip(tree, which(duplicated(tree$tip.label)))


# Plot tree ---------------------------------------------------------------
p.1 <- ggtree(family.tree, layout = "fan", ladderize = TRUE) +
  geom_tiplab(offset = 5)

# Highlight orders --------------------------------------------------------
node <- sapply(unique(mam$Order.1.2), function(Order) getMRCA(family.tree, unique(mam$Family.1.2[which(mam$Order.1.2 == Order)]))) %>% unlist
node <- as_tibble(node, rownames = "Order.1.2")

order.of.orders <- unique(mam$Order.1.2[match(get_taxa_name(p.1), mam$Family.1.2)])
node$Order.1.2 <- as_factor(node$Order.1.2)
node$Order.1.2 <- fct_reorder(node$Order.1.2, match(node$Order.1.2, order.of.orders))

p.2 <- p.1 +
  geom_hilight(data=node, aes(node=value, fill=Order.1.2))

### INFO FOR PLOTTING DENSITY
family.order <- rev(get_taxa_name(p.2))
fill <- levels(node$Order.1.2)
plot.info <- list(family.order = family.order, fill = fill)
write_rds(plot.info, "builds/uncertainty.plot.info.rds")


# Add point with color for mean uncertainty of FMR in family -------------
fam.mr.sd <- mr %>% group_by(Family.1.2) %>% summarise(fam.sd = mean(sd.fmr))

p.3 <- p.2  %<+% fam.mr.sd +
  geom_tippoint(aes(color=fam.sd), size=3, position = position_nudge(x = 3)) +
  scale_color_continuous("Mean logFMR standard deviation", low = "black", high = "red")
p.3


# Save plot --------------------------------------------------------------
ggsave("output/appendix2_FMR_uncertainty.png", p.3, width = 35, height = 30, units = "cm")
