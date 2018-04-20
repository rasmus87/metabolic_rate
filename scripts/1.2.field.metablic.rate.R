# Compile FMR: kJ/day

library(tidyverse)
library(stringr)

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

fmr <- read_csv("data/MR/jane12086-sup-0003-AppendixS5.csv", col_types = cols())
names(fmr) <- make.names(names(fmr))
fmr <- fmr %>%
  filter(Class == "Mammalia") %>% 
  mutate(Binomial = paste(Genus, Species, sep = "_")) %>% 
  mutate(log10BM = log10(M..kg.*1000), log10FMR = log10(FMR..kJ...day.)) %>% 
  mutate(Binomial.1.2 = str_replace(Binomial, "Lama_glama", "Lama_guanicoe")) %>% 
  mutate(FMR.source = "Hudson.2013")

stopifnot(anti_join(fmr, mam, by = "Binomial.1.2") %>% nrow == 0)

fmr <- left_join(fmr, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
fmr <- fmr %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM, FMR.source)

nagy <- read_csv("data/MR/Nagy1999.csv", col_types = cols())
names(nagy) <- make.names(names(nagy))
tax.nagy <- read_csv("data/MR/Nagy1999_tax.solver.csv", col_types = cols())
nagy <- left_join(nagy, tax.nagy, by = c("Species", "Common.name"))
nagy <- nagy %>% mutate(Species = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
nagy <- nagy %>% mutate(Binomial.1.2 = str_replace(Species, " ", "_"))
stopifnot(anti_join(nagy, mam, by = "Binomial.1.2") %>% nrow == 0)

nagy <- nagy %>% filter(!Binomial.1.2 %in% fmr$Binomial.1.2)
nagy <- left_join(nagy, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
nagy <- nagy %>% mutate(log10BM = log10(Mass..g.), log10FMR = log10(FMR..kJ.day.))
nagy <- nagy %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM) %>% 
  mutate(FMR.source = "Nagy.1999")

fmr <- rbind(fmr, nagy)


cap <- read_tsv("data/MR/Capellini2010.txt", col_types = cols(), na = "-9999")
names(cap) <- make.names(names(cap))
cap <- cap %>% mutate(Species = str_replace(Species, "(.* \\()", ""))
cap <- cap %>% mutate(Species = str_replace(Species, "\\)", ""))
cap <- cap %>% mutate(Species = str_replace(Species, " ", "_"))
# Load taxonomy solver
cap.tax <- read_csv("data/MR/Capellini2010_tax.solver.csv", col_types = cols())
cap <- left_join(cap, cap.tax %>% select(Species, Binomial.1.2), by = "Species")
cap <- cap %>% mutate(Binomial.1.2 = ifelse(is.na(Binomial.1.2), Species, Binomial.1.2))
stopifnot(nrow(anti_join(cap, mam, by = "Binomial.1.2")) == 0)

cap <- cap %>% filter(!Binomial.1.2 %in% fmr$Binomial.1.2)
cap <- left_join(cap, mam %>% dplyr::select(Binomial.1.2:Family.1.2), by = "Binomial.1.2")
cap <- cap %>% mutate(log10BM = log10(Body.mass.for.FMR..gr.), log10FMR = log10(FMR..kJ.day.))
cap <- cap %>% filter(!is.na(log10FMR))
cap <- cap %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM) %>% 
  mutate(FMR.source = "Capellini.2010")

fmr <- rbind(fmr, cap)

add <- read_csv("data/MR/additions.csv", col_types = cols())
add <- add %>% mutate(log10BM = log10(mass.g), log10FMR = log10(fmr.kJ.day)) %>% 
  left_join(mam, by = "Binomial.1.2") %>% 
  select(Binomial.1.2, Order.1.2, Family.1.2, log10FMR, log10BM, FMR.source)

fmr <- rbind(fmr, add)

write_csv(fmr, "builds/fmr.csv")


################################################################################

# Look at the data:
ggplot(fmr, aes(log10BM, log10FMR)) + 
  geom_point(aes(shape = FMR.source == "Pagano.2018",
                 col = Family.1.2)) +
  geom_smooth(method = "loess", col = "black") +
  geom_smooth(method = "lm", col = "black", lty = "dashed")

# Model and calculate prediction interval
m <- lm(log10FMR ~ log10BM, data = fmr)
hc.fmr <- hc %>% left_join(mam) %>% transmute(log10BM = log10(Mass.g),
                                              log10FMR = predict(m, data.frame(log10BM)),
                                              lwr = predict(m, data.frame(log10BM), interval = "predict")[, 2],
                                              upr = predict(m, data.frame(log10BM), interval = "predict")[, 3])

# Merge with Halls Cave data to look at how well they are represented
hc <- readxl::read_xlsx("data/halls cave clean.xlsx", sheet = "hall cave clean")
hc <- hc %>% mutate(Binomial.1.2 = paste(genus.new, 
                                         ifelse(is.na(likely.sp), species.new, likely.sp),
                                         sep = "_"))
hc <- left_join(hc, mam %>% dplyr::select(Binomial.1.2:Family.1.2), col_types = cols())
anti_join(hc, fmr, by = "Order.1.2") %>% select(Order.1.2, Family.1.2, Binomial.1.2)

ggplot(fmr, aes(log10BM, log10FMR, col = FMR.source)) +
  geom_ribbon(data = hc.fmr, aes(col = NULL, shape = NULL, ymin = lwr, ymax = upr, group = NULL), alpha = 0.2) +
  geom_point(data = hc.fmr, aes(col = NULL, shape = NULL), col = "red", pch = 3, stroke = 1.3) +
  geom_line(data = hc.fmr, aes(col = NULL, shape = NULL), col = "red") +
  geom_point()



# Not run: Exclude some things
# fmr <-  fmr %>%
#   # Order Cetacea (Whales s.l.)
#   filter(!Family.1.2 %in% c("Balaenidae", "Balaenopteridae", "Delphinidae",
#                             "Eschrichtiidae", "Iniidae", "Monodontidae",
#                             "Neobalaenidae", "Phocoenidae", "Physeteridae",
#                             "Ziphiidae")) %>%
#   # Families in the clade Pinnipedia (Seal s.l.):
#   # Odobenidae (walruses)
#   # Otariidae (fur seals and sea lions)
#   # Phocidae (true seals)
#   filter(!Family.1.2 %in% c("Odobenidae", "Otariidae", "Phocidae")) %>%
#   # Order Sirenia (Sea cows s.l.):
#   filter(Order.1.2 != "Sirenia") %>%
#   # Order Chiroptera (Bats)
#   filter(Order.1.2 != "Chiroptera")