# Compile FMR: kJ/day

library(tidyverse)
library(stringr)

mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

## Load Hudson.2013
hudson <- read_csv("data/jane12086-sup-0003-AppendixS5.csv", col_types = cols())
names(hudson) <- make.names(names(hudson))
hudson <- hudson %>%
  filter(Class == "Mammalia") %>% 
  mutate(Binomial.Source = paste(Genus, Species, sep = "_")) %>% 
  mutate(BM = M..kg.*1000, log10BM = log10(M..kg.*1000), 
         FMR = FMR..kJ...day., log10FMR = log10(FMR..kJ...day.)) %>% 
  mutate(Binomial.1.2 = str_replace(Binomial.Source, "Lama_glama", "Lama_guanicoe")) %>% 
  mutate(FMR.source = "Hudson.2013")
stopifnot(anti_join(hudson, mam, by = "Binomial.1.2") %>% nrow == 0)

hudson <- left_join(hudson, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
hudson <- hudson %>% transmute(Binomial.1.2, Order.1.2, Family.1.2, Binomial.Source, BM, FMR, log10BM, log10FMR, Source = FMR.source)

## Load Hudson.1999
nagy <- read_csv("data/Nagy1999.csv", col_types = cols())
names(nagy) <- make.names(names(nagy))
tax.nagy <- read_csv("data/Nagy1999_tax.solver.csv", col_types = cols())
nagy <- left_join(nagy, tax.nagy, by = c("Species", "Common.name"))
nagy <- nagy %>% mutate(Binomial.Source = Species, Species = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
nagy <- nagy %>% mutate(Binomial.1.2 = str_replace(Species, " ", "_"))
stopifnot(anti_join(nagy, mam, by = "Binomial.1.2") %>% nrow == 0)

nagy <- left_join(nagy, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
nagy <- nagy %>% mutate(BM = Mass..g., log10BM = log10(Mass..g.), FMR = FMR..kJ.day., log10FMR = log10(FMR..kJ.day.))
nagy <- nagy %>% transmute(Binomial.1.2, Order.1.2, Family.1.2, Binomial.Source, BM, FMR, log10BM, log10FMR, Source = "Nagy.1999")

## Load Capellini.2010
cap <- read_tsv("data/Capellini2010.txt", col_types = cols(), na = "-9999")
names(cap) <- make.names(names(cap))
cap <- cap %>% mutate(Species = str_replace(Species, "(.* \\()", ""))
cap <- cap %>% mutate(Species = str_replace(Species, "\\)", ""))
cap <- cap %>% mutate(Species = str_replace(Species, " ", "_"))
# Load taxonomy solver
cap.tax <- read_csv("data/Capellini2010_tax.solver.csv", col_types = cols())
cap <- left_join(cap, cap.tax %>% select(Species, Binomial.1.2), by = "Species")
cap <- cap %>% mutate(Binomial.Source = Species, Binomial.1.2 = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
stopifnot(nrow(anti_join(cap, mam, by = "Binomial.1.2")) == 0)

cap <- left_join(cap, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
cap <- cap %>% mutate(BM = Body.mass.for.FMR..gr., log10BM = log10(Body.mass.for.FMR..gr.), FMR = FMR..kJ.day., log10FMR = log10(FMR..kJ.day.))
cap <- cap %>% filter(!is.na(log10FMR))
cap <- cap %>% transmute(Binomial.1.2, Order.1.2, Family.1.2, Binomial.Source, BM, FMR, log10BM, log10FMR, Source = "Capellini.2010")

## Load additions
add <- read_csv("data/additions.fmr.csv", col_types = cols())
add <- add %>% mutate(BM = mass.g,
                      log10BM = log10(mass.g),
                      FMR = fmr.kJ.day,
                      log10FMR = log10(fmr.kJ.day)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  transmute(Binomial.1.2, Order.1.2, Family.1.2, Binomial.Source = Binomial.1.2, BM, FMR, log10BM, log10FMR, Source = FMR.source, Comment = FMR.comment)

fmr <- bind_rows(hudson, cap, nagy, add)

bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

terrestrial <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                              !Family.1.2 %in% c(whale.families, seal.families),
                              !Binomial.1.2 %in% marine.carnivores) %>% pull(Binomial.1.2)
fmr <- fmr %>% filter(Binomial.1.2 %in% terrestrial)

d <- which(duplicated(fmr[c("Binomial.1.2", "log10BM", "log10FMR")]))
fmr <- fmr[-d, ]

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

write_csv(fmr, "builds/fmr2.csv")


################################################################################

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
