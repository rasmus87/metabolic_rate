# All mammal nutrient uptake

library(tidyverse)
library(stringr)
library(readxl)

n.samples = 1000

# Load Phylacine
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
mam <- mam %>% filter(Terrestrial == 1)
mam <- mam %>% filter(Diet.Plant >= 80)
#mam <- mam %>% filter(Order.1.2 != "Proboscidea")

# Load FMR and sample distribution
fmr <- read_csv("../metabolic_rate/builds/imputed.metabolic.rate.csv", col_types = cols())
fmr <- fmr[match(mam$Binomial.1.2, fmr$Binomial.1.2), ]

df.fmr <- 100
fmr <- fmr %>% mutate(se.fmr = (log10FMR_uprCI - log10FMR_lwrCI)/2/qt(0.975, df.fmr))

log10FMR_samples <- sapply(1:nrow(fmr), 
                           FUN = function(i) {
                             rnorm(n = n.samples, 
                                   mean = fmr$log10FMR_est[i], 
                                   sd = fmr$se.fmr[i])
                           }
)
log10FMR_samples <- t(log10FMR_samples)

# Load animal density and sample distribution
dens <- read_csv("../predict_density/output/animal.density.km2.csv", col_types = cols())
dens <- dens[match(mam$Binomial.1.2, dens$Binomial.1.2), ]

log10dens_samples <- sapply(1:nrow(dens), 
                            FUN = function(i) {
                              rnorm(n = n.samples, 
                                    mean = dens$log10dens.est[i], 
                                    sd = dens$se.fit[i])
                            }
)
log10dens_samples <- t(log10dens_samples)


# Calculate metabolic requirements:
# MR = Assimilation.fraction * Food.intake * energy.content <=>
# Food.intake = MR / (Assimilation.fraction * energy.content)
# Energy content (Dry matter):
# Protein and carbohydrates: 4 kcal/g
# Fat: 9 kcal/g
# For assimilation fraction refs see:
# Nagy, K. A. (1987). Field Metabolic Rate and Food Requirement Scaling in
# Mammals and Birds. Source: Ecological Monographs Ecological Monographs,
# 57(572), 111–128.
# Howler monkeys are around 35 %:
# Nagy, K. A., Milton, K., Nagy, K. A., & Milton, K. (2016). Energy Metabolism
# and Food Consumption by Wild Howler Monkeys, 60(3), 475–480.
Assimilation.fraction = 35/100

#Degen, A. A., Benjamin, R. W., Abdraimov, S. A., & Sarbasov, T. I. (2002).
#Browse selection by Karakul sheep in relation to plant composition and
#estimated metabolizable energy content. Journal of Agricultural Science,
#139(3), 353–358. https://doi.org/10.1017/S0021859602002551
#Sheep ME; 8.5 # kJ / g , sd = 1.4

# Food component	Energy density
# http://www.fao.org/docrep/006/Y5022E/y5022e04.htm
#               kJ/g	kcal/g
# Fat	          37 	  9
# Ethanol	      29	  7
# Proteins	    17	  4
# Carbohydrates	17	  4
# Organic acids	13	  3
# Polyols     	10	  2.4
# Fiber       	 8	  2

# Assuming almost pure carb/protein diet.
Energy.content = 17 # kJ/g
Energy.content = Energy.content * 1000 # kJ / kg

ME = (Energy.content * Assimilation.fraction)
# From sheeps diet:
ME = 8500 # kJ / kg
ME.samples <- matrix(rep(rnorm(n.samples, 8.5, 1.4) * 1000, nrow(mam)), byrow=T, nrow=nrow(mam))

biomass.consumption.kgC.yr.samples <- 10^log10FMR_samples * 365.25 / ME.samples

density.samples <- 10^log10dens_samples
Q.samples = density.samples * biomass.consumption.kgC.yr.samples
colnames(Q.samples) <- paste0("Q.sample.", 1:n.samples)
consumption.samples <- bind_cols(Binomial.1.2 = mam$Binomial.1.2, as_data_frame(Q.samples))

biomass.consumption.kgC.yr <- exp((log(10) * fmr$se.fmr)^2/2) * 10^fmr$log10FMR_est * 365.25 / ME # New-New corrected version
density <- exp((log(10) * dens$se.fit)^2/2) * 10^dens$log10dens.est # New-New corrected version
Q = density * biomass.consumption.kgC.yr
consumption <- bind_cols(Binomial.1.2 = mam$Binomial.1.2, Q = Q)

mam <- left_join(mam, consumption, by = "Binomial.1.2")

ggplot(mam, aes(Mass.g, Q, col = Order.1.2)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(aes(col = NA))



which(consumption.samples$Binomial.1.2 == "Capreolus_capreolus")
which(consumption.samples$Binomial.1.2 == "Mammuthus_primigenius")

X.roe <- consumption.samples[1248, -1]/consumption.samples[318, -1]
summary(t(X.roe))
X.0 <- 10^log10dens_samples[318,]
X.1 <- 10^log10dens_samples[318,]*X.roe %>% as.numeric()
X.3 <- 10^log10dens_samples[1248,]

ggplot(data.frame(), aes()) +
  geom_boxplot(aes(x = "0", y = X.0)) +
  geom_boxplot(aes(x = "1", y = X.1)) +
  geom_boxplot(aes(x = "3", y = X.3)) +
  scale_y_log10()
