library(tidyverse)

dataset <- read_csv("builds/mr.csv", col_types = cols())

imputed <- read_csv("builds/333_mr_post.pred.csv")

mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

bat.order <- "Chiroptera"
sea.cow.order <- "Sirenia"
whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
                    "Neobalaenidae", "Delphinidae", "Monodontidae", 
                    "Eschrichtiidae", "Iniidae", "Physeteridae", 
                    "Phocoenidae", "Platanistidae")
seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")

mam <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
                      !Family.1.2 %in% c(whale.families, seal.families),
                      !Binomial.1.2 %in% marine.carnivores)

imputed.bmr <- imputed[, 1:nrow(mam)]
imputed.fmr <- imputed[, (1+nrow(mam)):(2*nrow(mam))]

bmr <- imputed.bmr %>% gather("Binomial.1.2", "log10bmr")
fmr <- imputed.fmr %>% gather("Binomial.1.2", "log10fmr")

imputed.bmr <- as.matrix(imputed.bmr)
attr(imputed.bmr, "class") <- "mcmc"
imputed.bmr.ci <- coda::HPDinterval(imputed.bmr)
imputed.bmr.ci <- imputed.bmr.ci %>% as_tibble(rownames = "Binomial.1.2")
imputed.bmr.ci <- imputed.bmr.ci %>% transmute(Binomial.1.2,
                                       log10.bmr.lower.95hpd = lower,
                                       log10.bmr.upper.95hpd = upper)

imputed.fmr <- as.matrix(imputed.fmr)
attr(imputed.fmr, "class") <- "mcmc"
imputed.fmr.ci <- coda::HPDinterval(imputed.fmr)
imputed.fmr.ci <- imputed.fmr.ci %>% as_tibble(rownames = "Binomial.1.2")
imputed.fmr.ci <- imputed.fmr.ci %>% transmute(Binomial.1.2,
                                       log10.fmr.lower.95hpd = lower,
                                       log10.fmr.upper.95hpd = upper)

bmr.summary <- bmr %>%
  group_by(Binomial.1.2) %>% 
  summarise(log10.bmr.median = median(log10bmr),
            log10.bmr.mean = mean(log10bmr),
            sd = sd(log10bmr)) %>% 
  left_join(imputed.bmr.ci, by = "Binomial.1.2")

fmr.summary <- fmr %>%
  group_by(Binomial.1.2) %>% 
  summarise(log10.fmr.median = median(log10fmr),
            log10.fmr.mean = mean(log10fmr),
            sd = sd(log10fmr)) %>% 
  left_join(imputed.fmr.ci, by = "Binomial.1.2")

mr.summary <- bind_cols(bmr.summary, fmr.summary[, -1])

mam.mr <- mam %>% 
  transmute(Binomial.1.2, Order.1.2, Family.1.2, log10BM = log10(Mass.g)) %>% 
  left_join(mr.summary, by = "Binomial.1.2")
mam.mr <- mam.mr %>% mutate(bmr.median= 10^log10.bmr.median,
                            bmr.mean = 10^log10.bmr.mean,
                            bmr.lower.95hpd = 10^log10.bmr.lower.95hpd,
                            bmr.upper.95hpd = 10^log10.bmr.upper.95hpd,
                            fmr.median= 10^log10.fmr.median,
                            fmr.mean = 10^log10.fmr.mean,
                            fmr.lower.95hpd = 10^log10.fmr.lower.95hpd,
                            fmr.upper.95hpd = 10^log10.fmr.upper.95hpd)

write_csv(mam.mr, "builds/imputed.metabolic.rate.csv")


ggplot(mam.mr, aes(x = log10BM, col = Order.1.2)) +
  geom_linerange(aes(ymin = log10.bmr.lower.95hpd, ymax = log10.bmr.upper.95hpd), lty = 3) +
  geom_point(aes(y = log10.bmr.median)) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Field Metabolic rate (kJ / day)")

ggplot(mam.mr, aes(x = log10BM)) +
  geom_linerange(aes(ymin = log10.bmr.lower.95hpd, ymax = log10.bmr.upper.95hpd, col = "BMR"), lty = 3) +
  geom_point(aes(y = log10.bmr.median, col = "BMR")) +
  geom_linerange(aes(ymin = log10.fmr.lower.95hpd, ymax = log10.fmr.upper.95hpd, col = "FMR"), lty = 3) +
  geom_point(aes(y = log10.fmr.median, col = "FMR")) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Metabolic rate (kJ / day)") +
  scale_color_discrete(c("FMR", "BMR"))

ggplot(mam.mr, aes(x = log10BM)) +
  geom_point(aes(y = log10.fmr.median, col = Order.1.2), pch = 1) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Field metabolic rate (kJ / day)")

col9 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6")
col27 <- rep(col9, each = 3)
pch3 <- c(1, 2, 3)
pch27 <- rep(pch3, times = 9)

ggplot(mam.mr, aes(x = log10BM, col = Order.1.2, shape = Order.1.2)) +
  geom_point(aes(y = log10.fmr.median)) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^2))) +
  scale_color_manual(values = col27) +
  scale_shape_manual(values = pch27) +
  theme(legend.position = c(.2,.8), legend.background = element_rect(linetype="solid", colour = "black")) +
  guides(col=guide_legend(nrow=9))


ggplot(mam.mr, aes(x = log10BM, col = Binomial.1.2 %in% dataset$Binomial.1.2)) +
  geom_point(aes(y = log10.fmr.median), pch = 19) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Density~(km^-2))) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "none")


full.data.bmr <- dataset %>%
  filter(MR == "BMR") %>% 
  transmute(Binomial.1.2, log10.bmr.data = log10MR) %>% 
  right_join(mam.mr)
full.data.fmr <- dataset %>%
  filter(MR == "FMR") %>% 
  transmute(Binomial.1.2, log10.fmr.data = log10MR) %>% 
  right_join(mam.mr)

ggplot(full.data.bmr, aes(x = log10.bmr.data, y = log10.bmr.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_bw() +
  xlab(expression(log[10]~Empirical~BMR~(kJ/day))) +
  ylab(expression(log[10]~Imputed~median~BMR~(kJ/day))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal()

ggplot(full.data.fmr, aes(x = log10.fmr.data, y = log10.fmr.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_bw() +
  xlab(expression(log[10]~Empirical~FMR~(kJ/day))) +
  ylab(expression(log[10]~Imputed~median~FMR~(kJ/day))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal()
