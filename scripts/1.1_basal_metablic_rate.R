# Compile BMR dataset from various sources
# 04/08-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
library(ggrepel)
library(ggpmisc)

# Load PHYLACINE 1.2.1
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Load taxonomy solver
syn <- read_csv("../PHYLACINE_1.2/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv")
syn <- syn %>% 
  transmute(Binomial.1.2,
            Binomial.EltonTraits.1.0 = paste(EltonTraits.1.0.Genus, EltonTraits.1.0.Species, sep = "_"),
            Binomial.1.0 = paste(Genus.1.0, Species.1.0, sep = "_"),
            Binomial.1.1 = paste(Genus.1.1, Species.1.0, sep = "_")) %>% 
  filter(!str_detect(Binomial.1.2, "000"))
syn <- syn %>% 
  gather(System, Binomial, -Binomial.1.2) %>% 
  filter(Binomial.1.2 != Binomial) %>% 
  filter(!duplicated(Binomial)) %>% 
  bind_rows(mam %>% transmute(Binomial.1.2,
                              System = "Binomial.1.2",
                              Binomial = Binomial.1.2))

# Load Capellini 2010 -----------------------------------------------------

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
cap <- cap %>% left_join(mam, by = "Binomial.1.2")
# Make tibble names uniform and remove data with missing BMR
cap <- cap %>% 
  transmute(Binomial.1.2, 
            Order.1.2, 
            Family.1.2,
            Binomial.Source = Species,
            BM = `Body mass for BMR (gr)`, 
            BMR = `BMR (mlO2/hour)`, 
            log10BM = log10(BM), 
            log10BMR = log10(BMR),
            Source = "Capellini.2010") %>% 
   filter(!is.na(log10BMR))


# Load PanTHERIA -----------------------------------------------------

# Load data
pantheria <- read_tsv("data/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
pantheria <- pantheria %>% 
  filter(!is.na(`18-1_BasalMetRate_mLO2hr`)) %>% 
  transmute(Binomial.Source = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            BM = `5-2_BasalMetRateMass_g`, 
            BMR = `18-1_BasalMetRate_mLO2hr`, 
            log10BM = log10(BM),
            log10BMR = log10(BMR),
            Source = "PanTHERIA.2008")

# Check that all PANTHERIA can be matched in syn-table
stopifnot(all(pantheria$Binomial.Source %in% syn$Binomial))
pantheria <- left_join(pantheria, syn %>% select(-System), by = c("Binomial.Source" = "Binomial"))
pantheria <- pantheria %>% left_join(mam, by = "Binomial.1.2")

# Make tibble names uniform and remove data with missing BMR
pantheria <- pantheria %>% 
  transmute(Binomial.1.2, 
            Order.1.2, 
            Family.1.2,
            Binomial.Source, 
            BM, 
            BMR, 
            log10BM, 
            log10BMR,
            Source)


# Load White & Seymour 2003 -----------------------------------------------------

# Load data
ws <- read_csv("data/White_Seymour2003.csv")
ws <- ws %>% 
  transmute(Binomial.Source = str_replace_all(binomial, " ", "_"),
            BM = mass.g,
            BMR = bmr.ml.O2.h,
            log10BM = log10(BM),
            log10BMR = log10(BMR),
            Source = "WhiteSeymore.2003")

# Load White & Seymour 2003 taxonomy solver
ws.tax <- read_csv("data/White_Seymour2003_tax.solver.csv")
ws <- ws %>% left_join(ws.tax, by = c("Binomial.Source" = "Species"))
# Update names where appropriate
ws <- ws %>% mutate(Binomial.1.2 = coalesce(Binomial.1.2, Binomial.Source))

# Add Family and Order from PHYLACINE
ws <- ws %>% left_join(mam, by = "Binomial.1.2")
# Make tibble names uniform and remove data with missing BMR
ws <- ws %>% 
  transmute(Binomial.1.2,
            Order.1.2, 
            Family.1.2,
            Binomial.Source,
            BM,
            BMR, 
            log10BM,
            log10BMR, 
            Source)



# Combine datasets --------------------------------------------------------

# Combine
bmr <- bind_rows(ws, pantheria, cap)

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
bmr <- bmr %>% filter(Binomial.1.2 %in% terrestrial)

# Visually inspect data
ggplot(bmr, aes(log10BM, log10BMR)) +
  geom_point() + 
  geom_smooth(method = lm) +
  stat_dens2d_filter(aes(label = paste(Source, Binomial.1.2, sep = ": ")),
                     geom = "text_repel", keep.fraction = 0.005)

# PanTHERIA 2008 has a crazy outlier for a single species which we remove
bmr <- bmr %>% filter(!(Source == "PanTHERIA.2008" & Binomial.1.2 == "Acrobates_pygmaeus"))

# Check for exact duplicated values and remove
bmr <- bmr %>% 
  distinct(Binomial.1.2, log10BM, log10BMR, .keep_all = TRUE)

# Check for pseudo duplicated values
# I.e. values that we think within data certainty might be duplicated
# Removes datapoins that are within 1 % of the mean datapoints for the species
duplicated.species <- bmr %>% 
  filter(duplicated(Binomial.1.2)) %>%
  pull(Binomial.1.2) %>% 
  unique()
for(species in duplicated.species) {
  select <- which(bmr$Binomial.1.2 == species)
  sub <- bmr[select, ]
  mean.sp.sq.val <- mean(sqrt(sub[, "log10BM"]^2 + sub[, "log10BMR"]^2)[[1]])
  difference <- dist(sub[, c("log10BM", "log10BMR")])
  difference <- difference / mean.sp.sq.val * 100
  difference <- as.matrix(difference)
  difference[upper.tri(difference, diag = TRUE)] <- NA
  difference <- (difference < 1) * 1
  duplicates <- select[which(rowSums(difference, na.rm = T) > 0)]
  if(length(duplicates) > 0) bmr <- bmr[-duplicates, ]
}



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

# Add additions to the dataset
bmr <- bind_rows(bmr, add)
bmr <- bmr %>% filter(Binomial.1.2 %in% terrestrial)


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