library(tidyverse)

fmr <- read_csv("builds/fmr2.csv", col_types = cols())
fmr <- fmr %>% mutate(MR = FMR,
                      log10MR = log10FMR,
                      MR.type = "FMR") %>% 
  select(-log10FMR, -FMR)
nrow(fmr)
length(unique(fmr$Binomial.1.2))

bmr <- read_csv("builds/bmr2.csv", col_types = cols())
bmr <- bmr %>% mutate(MR = BMR,
                      log10MR = log10BMR,
                      MR.type = "BMR") %>% 
  select(-log10BMR, -BMR)
nrow(bmr)
length(unique(bmr$Binomial.1.2))

mr <- bind_rows(fmr, bmr)
nrow(mr)
length(unique(mr$Binomial.1.2))

mr <- mr %>% transmute(Binomial.1.2, Order.1.2, Family.1.2, Binomial.Source, BM, MR, log10BM, log10MR, MR.type, Source, Comment)

write_csv(mr, "builds/mr2.csv")

# Plot data:
ggplot(mr, aes(log10BM, log10MR, col = MR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# Load PHYLACINE
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

mr.all <- mr  
mr <- mr %>% group_by(Binomial.1.2, MR) %>% summarise_at(c("log10BM", "log10MR"), mean)

m <- lm(log10MR ~ log10BM + MR, data = mr)
summary(m)

# Model and calculate prediction interval
mam <- mam %>% bind_rows(mam) %>% filter(Binomial.1.2 %in% terrestrial) %>% 
  mutate(MR = c(rep("BMR", n()/2), rep("FMR", n()/2)),
         log10BM = log10(Mass.g)) %>% 
  mutate(log10MR = predict(m, .),
         lwr = predict(m, ., interval = "predict")[, 2],
         upr = predict(m, ., interval = "predict")[, 3])



# Plot and extent to full mammalian range
segment <- predict(m, list(log10BM = c(6.5, 6.5), MR = c("BMR", "FMR")))
ggplot(mam, aes(log10BM, log10MR, group = MR)) +
  geom_point(data = mr.all, aes(log10BM, log10MR), col = "grey") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(y = log10MR, col = MR), lwd = 1) +
  geom_point(data = mr, aes(col = MR)) +
  theme_bw() +
  geom_segment(aes(x = 6.5,
                   xend = 6.5, 
                   y = segment[1],
                   yend = segment[2]),
               arrow = arrow(angle = 90, length = unit(0.2, "cm"), ends = "both")) +
  geom_text(aes(x = 6.5, y = mean(segment), 
                label = paste0("×", round(10^coef(m)[3], 2))),
            hjust = -.1) +
  scale_x_continuous("log10 Body mass [kg]") +
  scale_y_continuous("log10 Metabolic rate [kJ / day]") +
  scale_color_brewer(name = NULL, palette = "Dark2", 
                     breaks = c("FMR", "BMR"),
                     labels = c("Field metabolic rate", "Basal metabolic rate")) +
  theme(legend.position=c(.75, .25),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6))

# FMR is ~ 2.65 times higher than BMR
# Calc confidence interval
res <- coef(summary(m))
round(10^(res[3, 1]), 2)
round(10^(res[3, 1]-res[3, 2]), 2)
round(10^(res[3, 1]+res[3, 2]), 2)


# m1 simple model with offset between fmr and bmr
m1 <- lm(log10MR ~ log10BM + MR, data = mr)
# m2 interaction i.e. change in slope between fmr and bmr?
m2 <- lm(log10MR ~ log10BM * MR, data = mr)
# m3 3/4 scaling law?
m3 <- lm(log10MR ~ log10BM + MR + offset(3/4 * log10BM), data = mr)
# m4 2/3 scaling law?
m4 <- lm(log10MR ~ log10BM + MR + offset(2/3 * log10BM), data = mr)
# m5 mean between the two scaling laws?
m5 <- lm(log10MR ~ log10BM + MR + offset(8.5/12 * log10BM), data = mr)

anova(m1, m2)
anova(m1, m2, m3, m4, m5)

library(AICcmodavg)
aictab(list(m1, m2, m3, m4, m5), paste0("m", 1:5))
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

# Model 2 is worse in AIC values
# Both m3 and m4 have sigificant slope differences to our model