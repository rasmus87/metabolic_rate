# Compile BMR dataset from various sources
# 04/08-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ggrepel)
library(ggpmisc)

# Load PHYLACINE 1.2.1
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# # Load taxonomy solver
# syn <- read_csv("../PHYLACINE_1.2/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv")
# syn <- syn %>% 
#   transmute(Binomial.1.2,
#             Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
#             Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
#             Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
#   filter(!str_detect(Binomial.1.2, "000"))
# syn <- syn %>% 
#   gather(System, Binomial, -Binomial.1.2) %>% 
#   filter(Binomial.1.2 != Binomial) %>% 
#   filter(!duplicated(Binomial)) %>% 
#   bind_rows(mam %>% transmute(Binomial.1.2,
#                               System = "Binomial.1.2",
#                               Binomial = Binomial.1.2))


# Load Genoud 2018 -----------------------------------------------------

# Load data
gen <- readxl::read_xlsx("data/Genoud2018.xlsx", na = c("", "NA"))

# Clean
gen <- gen %>% 
  transmute(Binomial.Source = str_replace_all(`Species W&R`, " ", "_"),
            BM = `body mass`,
            BMR = BMR,
            log10BM = log10(BM),
            log10BMR = log10(BMR),
            Source = "Genoud.2018",
            Accepted = Accepted)

# Load Genoud 2018 taxonomy solver
gen.tax <- read_csv("data/Genoud2018_tax.solver.csv", col_types = cols())
gen <- left_join(gen, gen.tax, by = c("Binomial.Source" = "Binomial"))
# Update names where appropriate
gen <- gen %>% mutate(Binomial.1.2 = coalesce(Binomial.1.2, Binomial.Source))
# Check which are missing
gen[which(!gen$Binomial.1.2 %in% mam$Binomial.1.2),]
# Remove non matched
gen <- gen %>% filter(Binomial.1.2 %in% mam$Binomial.1.2)

# Add Data from PHYLACINE
gen <- gen %>% left_join(mam, by = "Binomial.1.2")

# Remove all unaccepted values from genera which already have a good genus match
accepted.genera <- gen %>% 
  filter(Accepted == 1) %>% 
  pull(Genus.1.2) %>% 
  unique()
# Remove all not accepted species that are already represented in Genoud
gen <- gen %>%
  filter(Accepted == 1 | !Genus.1.2 %in% accepted.genera)


# Make tibble names uniform and remove data with missing BMR
gen <- gen %>% 
  transmute(Binomial.1.2, 
            Genus.1.2,
            Order.1.2, 
            Family.1.2,
            Binomial.Source,
            BM,
            BMR, 
            log10BM, 
            log10BMR,
            Source,
            Accepted) %>% 
  filter(!is.na(log10BMR), !is.na(log10BM))



# Combine datasets --------------------------------------------------------

bmr <- gen

# Check for exact duplicated values and remove
bmr <- bmr %>% 
  distinct(Binomial.1.2, log10BM, log10BMR, .keep_all = TRUE)

# Visually inspect data
ggplot(bmr, aes(log10BM, log10BMR, col = Accepted == 1)) +
  geom_point() + 
  geom_smooth(method = lm) +
  stat_dens2d_filter(aes(label = paste(Source, Binomial.1.2, sep = ": ")),
                     geom = "text_repel", keep.fraction = 0.005)

# Subtract 12 % from all unaccepted species
bmr <- bmr %>% 
  mutate(BMR = ifelse(Accepted == 1, BMR, BMR * 0.88),
         log10BMR = log10(BMR))
# Visually inspect data
ggplot(bmr, aes(log10BM, log10BMR, col = Accepted == 1)) +
  geom_point() + 
  geom_smooth(method = lm) +
  stat_dens2d_filter(aes(label = paste(Source, Binomial.1.2, sep = ": ")),
                     geom = "text_repel", keep.fraction = 0.005)

# Show no seperation
ggplot(bmr, aes(log10BM, log10BMR)) +
  geom_point() + 
  geom_smooth(method = lm) +
  stat_dens2d_filter(aes(label = paste(Source, Binomial.1.2, sep = ": ")),
                     geom = "text_repel", keep.fraction = 0.005)


# Convert from [mLO2 / hour] to [kJ / day] --------------------------------

# Notes on estimating calories from O2
# [0.001 L / 60*60 s = mL/hr]
# Calories from O2
# https://www.csun.edu/~bby44411/346pdf/Ch4notes.pdf
# 5.05 kcal / L O2 # Assuming pure carb burn
# 4.69 kcal / L O2 # Assuming pure fat burn
# 4.87 kcal / L O2 # Assuming mixed burn
# 4184 J / kcal
# 4184 * c(5.05, 4.87, 4.69)

# Transform our data from mlO2 / hour to kJ/day
# kJ from O2
# http://www.jbc.org/content/59/1/41.full.pdf
# [kcal] = (3.815 + 1.2321 * RQ) * [L O2]
bmr.mlO2.hour <- bmr$BMR
bmr.lO2.hour <- bmr.mlO2.hour / 1000
RQ = .8 # Resting burn around .8 [.7-1]
bmr.kcal.hour <- bmr.lO2.hour * (3.815 + 1.2321 * RQ)
bmr.kcal.day <- bmr.kcal.hour * 24
bmr.kj.day <- bmr.kcal.day * 4.184

bmr$BMR <- bmr.kj.day
bmr$log10BMR <- log10(bmr.kj.day)

# Add extra data from smaller sources -------------------------------------

# Load additions
add <- read_csv("data/additions.bmr.csv", col_types = cols())
# Transform additions so they align
add <- add %>% 
  mutate(BM = mass.g,
         BMR = bmr.kJ.day,
         log10BM = log10(BM),
         log10BMR = log10(BMR)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  transmute(Binomial.1.2,
            Order.1.2, 
            Family.1.2, 
            Binomial.Source = Binomial.1.2,
            BM, 
            BMR, 
            log10BM, 
            log10BMR, 
            Source = BMR.source, 
            Comment = BMR.comment)

# # Tursiops_truncatus
# mlO2_kg_min <- 6.53
# mass.kg <- 148.6
# mlO2_kg_min / 1000 * 60 * mass.kg * (3.815 + 1.2321 * RQ) * 24 * 4.184
# 
# # Leptonychotes_weddelli
# mlO2_kg_min <- 3.58
# mass.kg <- 388.5
# mlO2_kg_min / 1000 * 60 * mass.kg * (3.815 + 1.2321 * RQ) * 24 * 4.184
# 
# # Orcinus_orca
# LO2_min <- 13.5
# LO2_min * 60 * (3.815 + 1.2321 * RQ) * 24 * 4.184
# 
# # Delphinapterus_leucas
# LO2_min <- 1.9
# LO2_min * 60 * (3.815 + 1.2321 * RQ) * 24 * 4.184
# 
# # Dromiciops_gliroides
# mlO2_g_hour <- 0.79
# mass.g <- 40.2
# mlO2_g_hour / 1000 * mass.g * (3.815 + 1.2321 * RQ) * 24 * 4.184

# Any of add species already in the dataset
add$Binomial.1.2[which(add$Binomial.1.2 %in% bmr$Binomial.1.2)]
bmr[which(bmr$Binomial.1.2 %in% add$Binomial.1.2), ]

# Add additions to the dataset
bmr <- bind_rows(bmr, add)


# # Check for pseudo duplicated values
# # I.e. values that we think within data certainty might be duplicated
# # Removes datapoins that are within 1 % of the mean datapoints for the species
# duplicated.species <- bmr %>% 
#   filter(duplicated(Binomial.1.2)) %>%
#   pull(Binomial.1.2) %>% 
#   unique()
# for(species in duplicated.species) {
#   select <- which(bmr$Binomial.1.2 == species)
#   sub <- bmr[select, ]
#   mean.sp.sq.val <- mean(sqrt(sub[, "log10BM"]^2 + sub[, "log10BMR"]^2)[[1]])
#   difference <- dist(sub[, c("log10BM", "log10BMR")])
#   difference <- difference / mean.sp.sq.val * 100
#   difference <- as.matrix(difference)
#   difference[upper.tri(difference, diag = TRUE)] <- NA
#   difference <- (difference < 1) * 1
#   duplicates <- select[which(rowSums(difference, na.rm = T) > 0)]
#   if(length(duplicates) > 0) bmr <- bmr[-duplicates, ]
# }


# Which families are not represented
mam %>% filter(!IUCN.Status.1.2 %in% c("EP", "DD", "EX", "EW"),
               !Family.1.2 %in% bmr$Family.1.2) %>% 
  select(1:3) %>% count(Family.1.2) %>% arrange(-n)

ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = Source %in% c("Genoud.2018"))) +
  geom_smooth(method = "loess", col = "black", lty = "dotted") +
  geom_smooth(method = "lm", col = "red")

ggplot(bmr, aes(log10BM, log10BMR)) +
  geom_point(aes(col = Binomial.1.2 == "Lagenorhynchus_obliquidens")) +
  geom_smooth(method = "loess", col = "black", lty = "dotted") +
  geom_smooth(method = "lm", col = "red")


# Write data out and check visually that it seems OK ----------------------

# Write data
write_csv(bmr, "builds/bmr_data.csv")

# Display data regading Order:
ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = Order.1.2, shape = Source)) +
  geom_smooth(method = "loess", col = "black", lty = "dotted") +
  geom_smooth(method = "lm", col = "red")

# Display data regarding data source:
ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = Source)) +
  geom_smooth(method = "loess", col = "black", lty = "dotted") +
  geom_smooth(method = "lm", col = "black")
