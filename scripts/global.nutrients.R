# All mammal nutrient uptake

library(tidyverse)
library(stringr)
library(readxl)

n.samples = 100

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





mam.sample <- left_join(mam, consumption.samples, by = "Binomial.1.2")
mam.sample <- gather(mam.sample, Q.sample, Q.est, Q.sample.1:paste0("Q.sample.", n.samples))
consumption.tot.sample <- mam.sample %>%
  group_by(Q.sample) %>%
  summarise(consumption = sum(Q.est))

consumption.tot.sample %>% group_by(time) %>% summarise(mean(consumption))

ci95 <- consumption.tot.sample %>%
  group_by(time) %>%
  summarise(q025 = quantile(consumption, 0.025),
            q975 = quantile(consumption, 0.975),
            q25 = quantile(consumption, 0.25),
            q75 = quantile(consumption, 0.75),
            median = median(consumption),
            mean = mean(consumption))

#y.max = consumption.tot.sample %>% filter(time == 13.8045) %>% pull(consumption) %>% quantile(0.975)


ggplot(consumption.tot, aes(time, consumption)) +
  #geom_line(data = consumption.tot.sample, aes(time, consumption, group = Q.sample), col = alpha("blue", .1)) +
  geom_point() +
  geom_line() +
  scale_x_reverse("Time (kaBP)") +
  #scale_y_continuous(bquote("Q: Consumption"~(MgCarbon / km^2 / year)), labels = function(x)x/1000) +
  scale_y_log10(bquote("Q: Consumption"~(MgCarbon / km^2 / year)), labels = function(x)x/1000) +
  theme_bw() +
  geom_hline(yintercept = npp.hc.kgC_km2_yr, lty = 2) +
  geom_text(y = npp.hc.kgC_km2_yr, x = 0, label = "Current NPP", hjust = 1, vjust = -1) +
  geom_line(data = ci95, aes(y = median), lwd = 1, col = "red") +
  geom_line(data = ci95, aes(y = q025), col = "red", lty = 2) +
  geom_line(data = ci95, aes(y = q975), col = "red", lty = 2) +
  geom_line(data = ci95, aes(y = q25), col = "red", lty = 3) +
  geom_line(data = ci95, aes(y = q75), col = "red", lty = 3) +
  geom_line(data = ci95, aes(y = mean), lwd = 1, col = "blue")

