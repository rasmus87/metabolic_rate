# Compile BMR

library(tidyverse)

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

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

cap <- left_join(cap, mam %>% dplyr::select(Binomial.1.2:Family.1.2))
cap <- cap %>% mutate(log10BM = log10(Body.mass.for.BMR..gr.), log10BMR = log10(BMR..mlO2.hour.))
cap <- cap %>% filter(!is.na(log10BMR))
cap <- cap %>% select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM) %>% 
  mutate(BMR.source = "Capellini.2010")


pantheria <- read_tsv("data/MR/PanTHERIA_1-0_WR05_Aug2008.txt", col_types = cols())
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
  select(Binomial.1.2, Order.1.2, Family.1.2, log10BMR, log10BM, BMR.source) %>% 
  filter(!Binomial.1.2 %in% cap$Binomial.1.2)

# More notes:
# 0.001 L / 60*60 s = mL/hr
# Calories from O2
# https://www.csun.edu/~bby44411/346pdf/Ch4notes.pdf
# 5.05 kcal / L O2 # Assuming pure carb burn
# 4.69 kcal / L O2 # Assuming pure fat burn
# 4.87 kcal / L O2 # Assuming mixed burn
# 4184 J / kcal
# 4184*c(5.05, 4.87, 4.69)

bmr <- rbind(cap, pantheria)
#http://www.jbc.org/content/59/1/41.full.pdf
# [kcal] = (3.815 + 1.2321 * RQ) * [L O2]
bmr.mlO2.hour <- 10^bmr$log10BMR
bmr.lO2.hour <- bmr.mlO2.hour / 1000
RQ = .8 # Resting burn around .8 [.7-1]
bmr.kcal.hour <- bmr.lO2.hour * (3.815 + 1.2321 * RQ)
bmr.kcal.day <- bmr.kcal.hour * 24
bmr.kj.day <- bmr.kcal.day * 4.184
bmr$log10BMR <- log10(bmr.kj.day)

write_csv(bmr, "builds/bmr.csv")


################################################################################

# Look at out dataset

hc <- readxl::read_xlsx("data/halls cave clean.xlsx", sheet = "hall cave clean")
hc <- hc %>% mutate(Binomial.1.2 = paste(genus.new, 
                                         ifelse(is.na(likely.sp), species.new, likely.sp),
                                         sep = "_"))
hc <- left_join(hc, mam %>% dplyr::select(Binomial.1.2:Family.1.2))
anti_join(hc, bmr, by = "Order.1.2") %>% select(Order.1.2, Family.1.2, Binomial.1.2)

ggplot(bmr, aes(log10BM, log10BMR)) + 
  geom_point(aes(col = Order.1.2, shape = BMR.source)) +
  geom_smooth(method = "loess", col = "blue") +
  geom_smooth(method = "lm", col = "red")

m <- lm(log10BMR ~ log10BM, data = bmr)
hc.fmr <- hc %>% left_join(mam) %>% transmute(log10BM = log10(Mass.g),
                                              log10BMR = predict(m, data.frame(log10BM)),
                                              lwr = predict(m, data.frame(log10BM), interval = "predict")[, 2],
                                              upr = predict(m, data.frame(log10BM), interval = "predict")[, 3])


ggplot(bmr, aes(log10BM, log10BMR)) +
  geom_ribbon(data = hc.fmr, aes(ymin = lwr, ymax = upr, group = NULL), alpha = 0.2) +
  geom_point(data = hc.fmr, col = "red", pch = 3, stroke = 1.3) +
  geom_line(data = hc.fmr, col = "red") +
  geom_point(aes(col = BMR.source))
