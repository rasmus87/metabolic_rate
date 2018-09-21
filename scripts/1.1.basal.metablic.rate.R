# Compile BMR

library(tidyverse)

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

cap <- read_tsv("data/Capellini2010.txt", col_types = cols(), na = "-9999")
names(cap) <- make.names(names(cap))
cap <- cap %>% mutate(Species = str_replace(Species, "(.* \\()", ""))
cap <- cap %>% mutate(Species = str_replace(Species, "\\)", ""))
cap <- cap %>% mutate(Species = str_replace(Species, " ", "_"))
# Load taxonomy solver
cap.tax <- read_csv("data/Capellini2010_tax.solver.csv", col_types = cols())
cap <- left_join(cap, cap.tax %>% select(Species, Binomial.1.2), by = "Species")
cap <- cap %>% mutate(Binomial.1.2 = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
stopifnot(nrow(anti_join(cap, mam, by = "Binomial.1.2")) == 0)

cap <- left_join(cap, mam %>% dplyr::select(Binomial.1.2:Family.1.2))
cap <- cap %>% mutate(log10BM = log10(Body.mass.for.BMR..gr.), log10BMR = log10(BMR..mlO2.hour.))
cap <- cap %>% filter(!is.na(log10BMR))
cap <- cap %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM) %>% 
  mutate(BMR.source = "Capellini.2010")


pantheria <- read_tsv("data/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
names(pantheria) <- make.names(names(pantheria))
pantheria <- pantheria %>% 
  filter(!is.na(X18.1_BasalMetRate_mLO2hr)) %>% 
  transmute(Binomial.Pantheria = paste(MSW05_Genus, MSW05_Species, sep = "_"),
            log10BM = log10(X5.2_BasalMetRateMass_g),
            log10BMR = log10(X18.1_BasalMetRate_mLO2hr),
            BMR.source = "PanTHERIA.2008")
# Load taxonomy solver
syn <- read_csv("../PHYLACINE_1.1/Data/Taxonomy/Synonymy_table_with_unaccepted_species.csv", col_types = cols(), guess_max = 5000)
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
stopifnot(!anti_join(pantheria, syn, by = c("Binomial.Pantheria" = "Binomial")) %>% nrow)

pantheria <- left_join(pantheria, syn %>% select(-System), by = c("Binomial.Pantheria" = "Binomial"))
pantheria <- pantheria %>% left_join(mam)

pantheria <- pantheria %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM, BMR.source)


ws <- read_csv("data/White_Seymour2003.csv")
ws <- ws %>% 
  transmute(Binomial.White.Seymour = str_replace_all(binomial, " ", "_"),
            log10BM = log10(mass.g),
            log10BMR = log10(bmr.ml.O2.h),
            BMR.source = "WhiteSeymore.2003")

# Check taxonomy
ws.tax <- read_csv("data/White_Seymour2003_tax.solver.csv")
ws <- ws %>% left_join(ws.tax, by = c("Binomial.White.Seymour" = "Species"))
ws$Binomial.1.2 <- ifelse(is.na(ws$Binomial.1.2), ws$Binomial.White.Seymour, ws$Binomial.1.2)

ws <- ws %>% left_join(mam)
ws <- ws %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM, BMR.source)
######################


bmr <- bind_rows(ws, pantheria, cap)
bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Physeteridae", "Phocoenidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

terrestrial <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                              !Family.1.2 %in% c(whale.families, seal.families),
                              !Binomial.1.2 %in% marine.carnivores) %>% pull(Binomial.1.2)
bmr <- bmr %>% filter(Binomial.1.2 %in% terrestrial)

## PanTHERIA 2008 has a crazy outlier for a single species which we remove:
bmr <- bmr %>% filter(!(BMR.source == "PanTHERIA.2008" & Binomial.1.2 == "Acrobates_pygmaeus"))

d <- which(duplicated(bmr[c("Binomial.1.2", "log10BM", "log10BMR")]))
bmr <- bmr[-d, ]

d <- bmr$Binomial.1.2[duplicated(bmr$Binomial.1.2)] %>% unique()
for(species in d) {
  select <- which(bmr$Binomial.1.2 == species)
  sub <- bmr[select,]
  difference <- dist(sub[, 4:5]) / sqrt(sub[1,4]^2 + sub[1,5]^2) * 100
  difference <- as.matrix(difference)
  difference[upper.tri(difference, diag = TRUE)] <- NA
  difference <- (difference < 1)*1
  duplicates <- select[which(rowSums(difference, na.rm = T) > 0)]
  if(length(duplicates) > 0) bmr <- bmr[-duplicates, ]
}

# More notes:
# 0.001 L / 60*60 s = mL/hr
# Calories from O2
# https://www.csun.edu/~bby44411/346pdf/Ch4notes.pdf
# 5.05 kcal / L O2 # Assuming pure carb burn
# 4.69 kcal / L O2 # Assuming pure fat burn
# 4.87 kcal / L O2 # Assuming mixed burn
# 4184 J / kcal
# 4184 * c(5.05, 4.87, 4.69)

#http://www.jbc.org/content/59/1/41.full.pdf
# [kcal] = (3.815 + 1.2321 * RQ) * [L O2]
bmr.mlO2.hour <- 10^bmr$log10BMR
bmr.lO2.hour <- bmr.mlO2.hour / 1000
RQ = .8 # Resting burn around .8 [.7-1]
bmr.kcal.hour <- bmr.lO2.hour * (3.815 + 1.2321 * RQ)
bmr.kcal.day <- bmr.kcal.hour * 24
bmr.kj.day <- bmr.kcal.day * 4.184
bmr$log10BMR <- log10(bmr.kj.day)

add <- read_csv("data/additions.bmr.csv", col_types = cols())
add <- add %>% mutate(log10BM = log10(mass.g), log10BMR = log10(bmr.kJ.day)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM, BMR.source)

bmr <- bind_rows(bmr, add)
bmr <- bmr %>% filter(Binomial.1.2 %in% terrestrial)

write_csv(bmr, "builds/bmr.csv")

# Display data:
ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = Order.1.2, shape = BMR.source)) +
  geom_smooth(method = "loess", col = "blue") +
  geom_smooth(method = "lm", col = "red")

# Display data:
ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = BMR.source)) +
  geom_smooth(method = "loess", col = "blue") +
  geom_smooth(method = "lm", col = "red")
