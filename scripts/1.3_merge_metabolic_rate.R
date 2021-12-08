# Compile BMR and FMR datasets into one
# 04/08-2021 Rasmus Ø Pedersen

# Load libraries
library(tidyverse)
 

# Load datasets -----------------------------------------------------------

# Load BMR
bmr <- read_csv("builds/bmr_data.csv")
bmr <- bmr %>% 
  mutate(MR = BMR,
         log10MR = log10BMR,
         MR.type = "BMR") %>% 
  select(-log10BMR, -BMR)
paste(nrow(bmr), "BMR datapoints for", length(unique(bmr$Binomial.1.2)), "unique species")

# Load FMR
fmr <- read_csv("builds/fmr_data.csv")
fmr <- fmr %>% 
  mutate(MR = FMR,
         log10MR = log10FMR,
         MR.type = "FMR") %>% 
  select(-log10FMR, -FMR)
paste(nrow(fmr), "FMR datapoints for", length(unique(fmr$Binomial.1.2)), "unique species")


# Combine datasets --------------------------------------------------------

mr <- bind_rows(bmr, fmr)
paste(nrow(mr), "MR datapoints for", length(unique(mr$Binomial.1.2)), "unique species")

mr <- mr %>% 
  transmute(Binomial.1.2,
            Order.1.2, 
            Family.1.2, 
            Binomial.Source, 
            BM, 
            MR, 
            log10BM,
            log10MR,
            MR.type, 
            Source, 
            Comment)


# Write data -----------------------------------------------------

write_csv(mr, "builds/metabolic_rate_data.csv")



# Plot test data and make simple not phylo models------------------

# Plot data:
ggplot(mr, aes(log10BM, log10MR, col = MR.type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# Load PHYLACINE 1.2.1
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Make summary per species
mr.species.summary <- mr %>% 
  group_by(Binomial.1.2, MR.type) %>% 
  summarise_at(c("log10BM", "log10MR"), mean)

# Fit model
m <- lm(log10MR ~ log10BM + MR.type, data = mr.species.summary)
summary(m)

# Model and calculate prediction interval
mam <- mam %>% 
  bind_rows(mam) %>% 
  mutate(MR.type = c(rep("BMR", n()/2), rep("FMR", n()/2)),
         log10BM = log10(Mass.g)) %>% 
  mutate(log10MR = predict(m, .),
         lwr = predict(m, ., interval = "predict")[, 2],
         upr = predict(m, ., interval = "predict")[, 3])

# Plot and extent to full mammalian range
segment <- c(7.5, 7.5)
segment <- c(segment, predict(m, list(log10BM = segment[1:2], MR.type = c("BMR", "FMR"))))
ggplot(mam, aes(log10BM, log10MR, group = MR.type)) +
  geom_point(data = mr, aes(log10BM, log10MR), col = "grey", pch = 21) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = MR.type), alpha = 0.2, show.legend = FALSE) +
  geom_line(aes(y = log10MR, col = MR.type), lwd = 1) +
  geom_point(data = mr.species.summary, aes(col = MR.type), pch = 21) +
  theme_bw() +
  geom_segment(aes(x = segment[1],
                   xend = segment[2], 
                   y = segment[3],
                   yend = segment[4]),
               arrow = arrow(angle = 90, length = unit(0.2, "cm"), ends = "both")) +
  geom_text(aes(x = segment[1], y = mean(segment[3:4]), 
                label = paste0("×", round(10^coef(m)[3], 2))),
            hjust = -.1) +
  scale_x_continuous("log10 Body mass [kg]") +
  scale_y_continuous("log10 Metabolic rate [kJ / day]") +
  scale_color_discrete(name = NULL, 
                     breaks = c("FMR", "BMR"),
                     labels = c("Field metabolic rate", "Basal metabolic rate")) +
  theme(legend.position=c(.75, .25),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6))

# Calc confidence interval on the factor offset between FMR and BMR
res <- coef(summary(m))
diff <- 10^(res[3, 1])
diff <- round(c(diff, diff - res[3, 2], diff + res[3, 2]), 2)
paste0("FMR is BMR x ", diff[1], " (", diff[2], "-", diff[3], ")")


# m1 simple model with offset between fmr and bmr
m1 <- lm(log10MR ~ log10BM + MR.type, data = mr.species.summary)
# m2 interaction i.e. change in slope between fmr and bmr?
m2 <- lm(log10MR ~ log10BM * MR.type, data = mr.species.summary)
# m3 3/4 scaling law?
m3 <- lm(log10MR ~ log10BM + MR.type + offset(3/4 * log10BM), data = mr.species.summary)
# m4 2/3 scaling law?
m4 <- lm(log10MR ~ log10BM + MR.type + offset(2/3 * log10BM), data = mr.species.summary)
# m5 mean between the two scaling laws?
m5 <- lm(log10MR ~ log10BM + MR.type + offset(8.5/12 * log10BM), data = mr.species.summary)

anova(m1, m2)
anova(m1, m2, m3, m4, m5)

library(AICcmodavg)
aictab(list(m1, m2, m3, m4, m5), paste0("m", 1:5))
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

