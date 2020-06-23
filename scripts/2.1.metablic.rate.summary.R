library(tidyverse)
library(ggtext)

dataset <- read_csv("builds/mr.csv", col_types = cols())

imputed <- read_csv("builds/3_mr_post.pred.csv")

mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

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
write_csv(imputed.bmr, "builds/3_bmr_post.pred.csv")
imputed.fmr <- imputed[, (1+nrow(mam)):(2*nrow(mam))]
names(imputed.fmr) <- names(imputed.bmr)
write_csv(imputed.fmr, "builds/3_fmr_post.pred.csv")

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

# write_csv(mam.mr, "builds/Table S5 Imputed metabolic rate.csv")
# mam.mr <- read_csv("builds/Table S5 Imputed metabolic rate.csv")


theme_R <- function() {
  theme_bw() %+replace% 
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
}

ggplot(mam.mr, aes(x = log10BM, col = Order.1.2)) +
  geom_linerange(aes(ymin = log10.bmr.lower.95hpd, ymax = log10.bmr.upper.95hpd), lty = 3) +
  geom_point(aes(y = log10.bmr.median)) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Field Metabolic rate (kJ / day)")


col9 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6")
col27 <- rep(col9, each = 3)
pch3 <- c(1, 2, 3)
pch27 <- rep(pch3, times = 9)
orders <- sort(unique(mam.mr$Order.1.2))
labels <- paste0("<span style='color:", col27, "'>", orders, "</span>")

ggplot(mam.mr, aes(x = log10BM, col = Order.1.2, shape = Order.1.2)) +
  geom_point(aes(y = log10.fmr.median)) +
  theme_R() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Imputed~median~FMR~(kJ/day))) +
  scale_color_manual(values = col27, breaks = orders, labels = labels) +
  scale_shape_manual(values = pch27, breaks = orders, labels = labels) +
  theme(legend.text = element_markdown()) +
  theme(legend.position = c(0, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(0,1)) +
  guides(col=guide_legend(nrow=9))
ggsave("output/appendix2_fig1_FMR_plot.png", width = 25.6, height = 14.4, units = "cm")


mr.types <- dataset %>% 
  group_by(Binomial.1.2, MR) %>% 
  summarise()
mr.types$MR[mr.types$MR == "BMR"] <- 1
mr.types$MR[mr.types$MR == "FMR"] <- 2
mr.types <- mr.types %>% 
  group_by(Binomial.1.2) %>% 
  summarise(mr.type = sum(as.numeric(MR)))
  
mam.mr <- mam.mr %>% 
  left_join(mr.types)
mam.mr$mr.type[is.na(mam.mr$mr.type)] <- 0
mam.mr$mr.type <- as_factor(mam.mr$mr.type)
ggplot(mam.mr, aes(x = log10BM, col = mr.type)) +
  geom_point(aes(y = log10.fmr.median), pch = 19, cex = 1) +
  theme_bw() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Imputed~median~FMR~(kJ/day))) +
  scale_color_brewer(palette = "Paired", name = "Imputed metabolic rate for species",
                     labels = c("without any known metabolic rate",
                                "with empirical basal metabolic rate",
                                "with empirical field metabolic rate",
                                "with empirical basal and field metabolic rate")) +
  theme_R() +
  theme(legend.position = c(0, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(-0.5, 1.5)) +
ggsave("output/appendix2_fig2_FMR_plot_known_data.png", width = 25.6, height = 14.4, units = "cm")



full.data.bmr <- dataset %>%
  filter(MR == "BMR") %>% 
  transmute(Binomial.1.2, log10.bmr.data = log10MR) %>% 
  right_join(mam.mr)
full.data.fmr <- dataset %>%
  filter(MR == "FMR") %>% 
  transmute(Binomial.1.2, log10.fmr.data = log10MR) %>% 
  right_join(mam.mr)

p1 <- ggplot(full.data.bmr %>% filter(!is.na(log10.bmr.data)),
             aes(x = log10.bmr.data, y = log10.bmr.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(formula = y ~ x, method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_R() +
  xlab(expression(log[10]~Empirical~BMR~(kJ/day))) +
  ylab(expression(log[10]~Imputed~median~BMR~(kJ/day))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal(xlim = c(.5, 5.5), ylim = c(.5, 5.5)) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*\' ,  \'*")), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3)

p2 <- ggplot(full.data.fmr %>% filter(!is.na(log10.fmr.data)), 
             aes(x = log10.fmr.data, y = log10.fmr.median, col = log10BM)) +
  geom_point(pch = 19) +
  geom_abline(slope = 1, lty = 2, lwd = .7) +
  geom_smooth(formula = y ~ x, method = lm, se = F, lty = 1, col = "black", lwd = .9) +
  theme_R() +
  xlab(expression(log[10]~Empirical~FMR~(kJ/day))) +
  ylab(expression(log[10]~Imputed~median~FMR~(kJ/day))) +
  scale_color_continuous(expression(log[10]~Body~mass~(g))) +
  coord_equal(xlim = c(.5, 5.5), ylim = c(.5, 5.5)) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*\' ,  \'*")), 
               label.x.npc = "left", label.y.npc = "top",
               formula = y ~ x, parse = TRUE, size = 3)

library(ggpubr)
main.plot <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("output/appendix2_fig3_data_vs_modelled.png", main.plot, width = 25.6, height = 14.4, units = "cm")


ggplot(mam.mr, aes(x = log10BM)) +
  geom_linerange(aes(ymin = log10.bmr.lower.95hpd, ymax = log10.bmr.upper.95hpd, col = "BMR"), lty = 3) +
  geom_point(aes(y = log10.bmr.median, col = "BMR"), cex = 1) +
  geom_linerange(aes(ymin = log10.fmr.lower.95hpd, ymax = log10.fmr.upper.95hpd, col = "FMR"), lty = 3) +
  geom_point(aes(y = log10.fmr.median, col = "FMR"), cex = 1) +
  theme_R() +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Metabolic~rate~(kJ/day))) +
  scale_color_discrete(name = NULL,
                     labels = c("Basal metabolic rate", "Field metabolic rate"),
                     breaks = c("BMR", "FMR")) +
  theme_R() +
  theme(legend.position = c(0, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(-0.5, 1.5))
ggsave("output/appendix2_fig5_FMR_BMR_plot.png", width = 25.6, height = 14.4, units = "cm")
