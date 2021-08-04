# Compile FMR dataset from various sources
# 30/07-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ggrepel)
library(ggpmisc)

# Load PHYLACINE 1.2.1
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())



# Load Hudson.2013 --------------------------------------------------------

# Load data
hudson <- read_csv("data/jane12086-sup-0003-AppendixS5.csv")
hudson <- hudson %>%
  filter(Class == "Mammalia") %>% 
  mutate(Binomial.Source = paste(Genus, Species, sep = "_"),
         BM = `M (kg)` * 1000, 
         FMR = `FMR (kJ / day)`, 
         log10BM = log10(BM), 
         log10FMR = log10(FMR),
         Binomial.1.2 = str_replace(Binomial.Source, "Lama_glama", "Lama_guanicoe"),
         Source = "Hudson.2013")
# Check that all are in there
stopifnot(all(hudson$Binomial.1.2 %in% mam$Binomial.1.2))

# Add Family and Order from PHYLACINE
hudson <- left_join(hudson, mam, by = "Binomial.1.2")
# Make tibble names uniform
hudson <- hudson %>% transmute(Binomial.1.2,
                               Order.1.2,
                               Family.1.2, 
                               Binomial.Source, 
                               BM, 
                               FMR, 
                               log10BM, 
                               log10FMR, 
                               Source)


# Load Nagy.1999 --------------------------------------------------------

# Load data
nagy <- read_csv("data/Nagy1999.csv")

# Load Nagy 1999 tax solver
tax.nagy <- read_csv("data/Nagy1999_tax.solver.csv")
nagy <- left_join(nagy, tax.nagy, by = c("Species", "Common.name"))
nagy <- nagy %>% mutate(Binomial.Source = Species, 
                        Species = coalesce(Binomial.1.2, Species),
                        Binomial.1.2 = str_replace(Species, " ", "_"))
# Check that all are in there
stopifnot(all(nagy$Binomial.1.2 %in% mam$Binomial.1.2))

# Add Family and Order from PHYLACINE
nagy <- left_join(nagy, mam, by = "Binomial.1.2")
# Make tibble names uniform
nagy <- nagy %>% transmute(Binomial.1.2,
                           Order.1.2, 
                           Family.1.2,
                           Binomial.Source,
                           BM = `Mass.(g)`, 
                           FMR = `FMR.(kJ/day)`, 
                           log10BM = log10(BM),
                           log10FMR = log10(FMR),
                           Source = "Nagy.1999")

# Load Capellini.2010 --------------------------------------------------------

# Load data
cap <- read_tsv("data/Capellini2010.txt", col_types = cols(), na = "-9999")
# Use names in paranthesis where available and replace spaces with underscore as PHYLACINE
cap <- cap %>% mutate(Species = str_replace(Species, "(.* \\()", ""),
                      Species = str_replace(Species, "\\)", ""),
                      Species = str_replace(Species, " ", "_"))

# Load Capellini 2010 taxonomy solver
cap.tax <- read_csv("data/Capellini2010_tax.solver.csv", col_types = cols())
cap <- left_join(cap, cap.tax %>% select(Species, Binomial.1.2), by = "Species")
# Update names where appropriate
cap <- cap %>% mutate(Binomial.1.2 = coalesce(Binomial.1.2, Species))
# Check that we got them all
stopifnot(all(cap$Binomial.1.2 %in% mam$Binomial.1.2))

# Add Family and Order from PHYLACINE
cap <- left_join(cap, mam, by = "Binomial.1.2")
# Make tibble names uniform and remove data with missing FMR
cap <- cap %>% 
  transmute(Binomial.1.2, 
            Order.1.2, 
            Family.1.2,
            Binomial.Source = Species,
            BM = `Body mass for FMR (gr)`, 
            FMR = `FMR (kJ/day)`, 
            log10BM = log10(BM), 
            log10FMR = log10(FMR),
            Source = "Capellini.2010") %>% 
  filter(!is.na(log10FMR))

# Add extra data from smaller sources -------------------------------------

# Load additions
add <- read_csv("data/additions.fmr.csv", col_types = cols())
# Transform additions so they align
add <- add %>% mutate(BM = mass.g,
                      FMR = fmr.kJ.day,
                      log10BM = log10(BM),
                      log10FMR = log10(FMR)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  transmute(Binomial.1.2, 
            Order.1.2, 
            Family.1.2,
            Binomial.Source = Binomial.1.2, 
            BM,
            FMR, 
            log10BM, 
            log10FMR, 
            Source = FMR.source, 
            Comment = FMR.comment)

# Combine datasets --------------------------------------------------------

# Combine
fmr <- bind_rows(hudson, cap, nagy, add)

# Non terrestrials and humans
bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")
humans <- "Homo"

# Terrestrial species list:
terrestrial <- mam %>% 
  filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
         !Family.1.2 %in% c(whale.families, seal.families),
         !Binomial.1.2 %in% marine.carnivores,
         !Genus.1.2 %in% humans) %>% 
  pull(Binomial.1.2)

# Filter BMR dataset to the terrestrial list
fmr <- fmr %>% filter(Binomial.1.2 %in% terrestrial)

# Check for exact duplicated values and remove
fmr <- fmr %>% 
  distinct(Binomial.1.2, log10BM, log10FMR, .keep_all = TRUE)

# Visually inspect data
ggplot(fmr, aes(log10BM, log10FMR)) +
  geom_point() + 
  geom_smooth(method = lm) +
  stat_dens2d_filter(aes(label = paste(Source, Binomial.1.2, sep = ": ")),
                     geom = "text_repel", keep.fraction = 0.005)


# Check for pseudo duplicated values
# I.e. values that we think within data certainty might be duplicated
# Removes datapoins that are within 1 % of the mean datapoints for the species
d <- fmr$Binomial.1.2[duplicated(fmr$Binomial.1.2)] %>% unique()
for(species in d) {
  select <- which(fmr$Binomial.1.2 == species)
  sub <- fmr[select,]
  mean.sp.sq.val <- mean(sqrt(sub[, "log10BM"]^2 + sub[, "log10FMR"]^2)[[1]])
  difference <- dist(sub[, c("log10BM", "log10FMR")])
  difference <- difference / mean.sp.sq.val * 100
  difference <- as.matrix(difference)
  difference[upper.tri(difference, diag = TRUE)] <- NA
  difference <- (difference < 1) * 1
  duplicates <- select[which(rowSums(difference, na.rm = T) > 0)]
  if(length(duplicates) > 0) fmr <- fmr[-duplicates, ]
}

# Write data out and check visually that it seems OK ----------------------

# Write data
write_csv(fmr, "builds/fmr_data.csv")

# Display data regarding data source:
ggplot(fmr, aes(log10BM, log10FMR)) + 
  geom_point(aes(col = Source)) +
  geom_smooth(method = "loess", col = "black") +
  geom_smooth(method = "lm", col = "black", lty = "dashed")

# Display data regading Order:
ggplot(fmr, aes(log10BM, log10FMR)) + 
  geom_point(aes(col = Order.1.2)) +
  geom_smooth(method = "loess", col = "black", lty = "dotted") +
  geom_smooth(method = "lm", col = "red")
