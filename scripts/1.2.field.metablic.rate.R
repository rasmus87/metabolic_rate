# Compile FMR: kJ/day

library(tidyverse)
library(stringr)

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

## Load Hudson.2013
hudson <- read_csv("data/jane12086-sup-0003-AppendixS5.csv", col_types = cols())
names(hudson) <- make.names(names(hudson))
hudson <- hudson %>%
  filter(Class == "Mammalia") %>% 
  mutate(Binomial = paste(Genus, Species, sep = "_")) %>% 
  mutate(log10BM = log10(M..kg.*1000), log10FMR = log10(FMR..kJ...day.)) %>% 
  mutate(Binomial.1.2 = str_replace(Binomial, "Lama_glama", "Lama_guanicoe")) %>% 
  mutate(FMR.source = "Hudson.2013")
stopifnot(anti_join(hudson, mam, by = "Binomial.1.2") %>% nrow == 0)

hudson <- left_join(hudson, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
hudson <- hudson %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM, FMR.source)

## Load Hudson.1999
nagy <- read_csv("data/Nagy1999.csv", col_types = cols())
names(nagy) <- make.names(names(nagy))
tax.nagy <- read_csv("data/Nagy1999_tax.solver.csv", col_types = cols())
nagy <- left_join(nagy, tax.nagy, by = c("Species", "Common.name"))
nagy <- nagy %>% mutate(Species = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
nagy <- nagy %>% mutate(Binomial.1.2 = str_replace(Species, " ", "_"))
stopifnot(anti_join(nagy, mam, by = "Binomial.1.2") %>% nrow == 0)

nagy <- left_join(nagy, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
nagy <- nagy %>% mutate(log10BM = log10(Mass..g.), log10FMR = log10(FMR..kJ.day.))
nagy <- nagy %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM) %>% 
  mutate(FMR.source = "Nagy.1999")

## Load Capellini.2010
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

cap <- left_join(cap, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
cap <- cap %>% mutate(log10BM = log10(Body.mass.for.FMR..gr.), log10FMR = log10(FMR..kJ.day.))
cap <- cap %>% filter(!is.na(log10FMR))
cap <- cap %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM) %>% 
  mutate(FMR.source = "Capellini.2010")

## Load additions
add <- read_csv("data/additions.fmr.csv", col_types = cols())
add <- add %>% mutate(log10BM = log10(mass.g), log10FMR = log10(fmr.kJ.day)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM, FMR.source)

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

d <- fmr$Binomial.1.2[duplicated(fmr$Binomial.1.2)] %>% unique()
for(species in d) {
  select <- which(fmr$Binomial.1.2 == species)
  sub <- fmr[select,]
  difference <- dist(sub[, 4:5]) / sqrt(sub[1,4]^2 + sub[1,5]^2) * 100
  difference <- as.matrix(difference)
  difference[upper.tri(difference, diag = TRUE)] <- NA
  difference <- (difference < 1)*1
  duplicates <- select[which(rowSums(difference, na.rm = T) > 0)]
  if(length(duplicates) > 0) fmr <- fmr[-duplicates, ]
}

write_csv(fmr, "builds/fmr.csv")


################################################################################

# Look at the data:
ggplot(fmr, aes(log10BM, log10FMR)) + 
  geom_point(aes(col = FMR.source)) +
  geom_smooth(method = "loess", col = "black") +
  geom_smooth(method = "lm", col = "black", lty = "dashed")
