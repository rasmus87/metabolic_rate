library(tidyverse)

# imputed <- read_csv("builds/imputed.metabolic.rate_all.samples.csv")
imputed <- read_csv("builds/test_imputed.metabolic.rate_all.samples.csv")

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())
mam <- mam %>% filter(Terrestrial == 1)

BMR <- imputed %>% 
  as.tbl %>% 
  filter(MR == "BMR") %>%
  transmute(Binomial.1.2, log10BMR_est = fit, 
            log10BMR_lwrCI = lwrCI, log10BMR_uprCI = uprCI,
            log10BMR_lwrPI = lwrPI, log10BMR_uprPI = uprPI) %>% 
  group_by(Binomial.1.2) %>%
  summarise_if(is.numeric, mean)
FMR <- imputed %>% 
  as.tbl %>% 
  filter(MR == "FMR") %>%
  transmute(Binomial.1.2, log10FMR_est = fit,
            log10FMR_lwrCI = lwrCI, log10FMR_uprCI = uprCI,
            log10FMR_lwrPI = lwrPI, log10FMR_uprPI = uprPI) %>% 
  group_by(Binomial.1.2) %>%
  summarise_if(is.numeric, mean)

mam.mr <- mam %>% 
  transmute(Binomial.1.2, Order.1.2, Family.1.2, log10BM = log10(Mass.g)) %>% 
  left_join(BMR) %>% 
  left_join(FMR)

ggplot(mam.mr, aes(x = log10BM, col = Order.1.2)) +
  geom_linerange(aes(ymin = log10FMR_lwrCI, ymax = log10FMR_uprCI)) +
  geom_linerange(aes(ymin = log10FMR_lwrPI, ymax = log10FMR_uprPI), lty = 3) +
  geom_point(aes(y = log10FMR_est)) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Field Metabolic rate (kJ / day)")

ggplot(mam.mr, aes(x = log10BM)) +
  geom_linerange(aes(ymin = log10FMR_lwrCI, ymax = log10FMR_uprCI, col = "FMR")) +
  geom_linerange(aes(ymin = log10FMR_lwrPI, ymax = log10FMR_uprPI, col = "FMR"), lty = 3) +
  geom_point(aes(y = log10FMR_est, col = "FMR"), pch = 1, col = "blue") +
  geom_linerange(aes(ymin = log10BMR_lwrCI, ymax = log10BMR_uprCI, col = "BMR")) +
  geom_linerange(aes(ymin = log10BMR_lwrPI, ymax = log10BMR_uprPI, col = "BMR"), lty = 3) +
  geom_point(aes(y = log10BMR_est, col = "BMR"), pch = 1, col = "red") +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Metabolic rate (kJ / day)") +
  scale_color_discrete(c("FMR", "BMR"))

write_csv(mam.mr, "builds/imputed.metabolic.rate.csv")

ggplot(mam.mr, aes(x = log10BM)) +
  geom_point(aes(y = log10FMR_est, col = Order.1.2), pch = 1) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Field metabolic rate (kJ / day)")
