library(tidyverse)

# solution <- read_csv("builds/metabolic.rate_fit.solution.csv")
solution <- read_csv("builds/3_mr_fit.solution.csv")

# Borrowed from the coda::HPDinterval() function
HPDinterval.mcmc <- function(obj, prob = 0.95, ...) {
  obj <- as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)],
               vals[cbind(inds + gap, 1:npar)])
  dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
  attr(ans, "Probability") <- gap/nsamp
  ans
}

# Result summary
# Posterier mean, and highest posterier density interval (95 %)
solution.summary <- cbind(mean = colMeans(solution[, 1:5]), HPDinterval.mcmc(solution[, 1:5]))
solution.summary <- cbind(median = apply(solution[, 1:5], 2, median), HPDinterval.mcmc(solution[, 1:5]))
signif(solution.summary, 2)

# In real values how much larger is FMR than BMR
round(10^solution.summary[3,], 2)

# Test effect of slope:
intercept = solution.summary[1, 1] 
slope = solution.summary[2, 1]
fmr = solution.summary[3, 1]
slope.change = solution.summary[4, 1]
random.effect = solution.summary[5, 1]

# log.bm <- 0:7
# y0 = intercept + random.effect + fmr + slope * log.bm
# y1 = intercept + random.effect + fmr + (slope + slope.change) * log.bm
# plot(log.bm, y0, type = "l")
# lines(log.bm, y1)
# cbind(log10BM.g = log.bm, difference = (1 - y1/y0) * 100) %>% as_data_frame()


# Plot chains
# Since we are now sampling a lot less gather all results per 10 trees
solution.gathered <- solution %>% 
  select(-random.effect) %>% 
  mutate(id = tree %/% 10) %>% 
  select(-tree, -chain)
solution.gathered <- gather(solution.gathered, key = "variable", value = "value", -id)
p <- ggplot(solution.gathered, aes(x = value, group = id)) +
  geom_line(stat = "density", col = alpha("blue", .10)) + # For fewer trees
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "") +
  geom_line(stat = "density", aes(group = NULL), lwd = 1) +
  theme_bw()
p
# ggsave("figures/metabolic.fit2.png", plot = p, height = 18, width = 12, dpi = 600, units = "cm")
ggsave("output/appendix2_fig6_avg_chains.png", width = 25.6, height = 21.6, units = "cm")

# Plot data and regression line:
# Load data
mr <- read_csv("builds/mr.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

# # Filter for terrestial
# bat.order <- "Chiroptera"
# sea.cow.order <- "Sirenia"
# whale.families <- c("Balaenidae", "Balaenopteridae", "Ziphiidae", 
#                     "Neobalaenidae", "Delphinidae", "Monodontidae", 
#                     "Eschrichtiidae", "Iniidae", "Physeteridae", 
#                     "Phocoenidae", "Platanistidae")
# seal.families <- c("Otariidae", "Phocidae", "Odobenidae")
# marine.carnivores <- c("Enhydra_lutris", "Lontra_felina", "Ursus_maritimus")
# 
# terrestrial <- mam %>% filter(!Order.1.2 %in% c(bat.order, sea.cow.order),
#                               !Family.1.2 %in% c(whale.families, seal.families),
#                               !Binomial.1.2 %in% marine.carnivores) %>% pull(Binomial.1.2)
# mr <- mr %>% 
#   filter(Binomial.1.2 %in% terrestrial)

mr.all <- mr  
mr <- mr %>% group_by(Binomial.1.2, MR) %>% summarise_at(c("log10BM", "log10MR"), mean)

ss <- as_tibble(t(solution.summary))

mam.mr <- read_csv("builds/Table S5 Imputed metabolic rate.csv")
bm.range <- range(mam.mr$log10BM)

ggplot(mr.all, aes(x = log10BM, y = log10MR)) +
  geom_point(col = "grey", alpha = .5, cex = 1) +
  geom_point(data = mr, aes(col = MR), alpha = .75, cex = 1, pch = 19) + 
  geom_line(data = data.frame(),
            aes(x = bm.range,
                y = ss$`(Intercept)`[1] + ss$random.effect[1] + ss$log10BM[1] * bm.range, col = "BMR")) +
  geom_ribbon(data = data.frame(),
              aes(x = bm.range,
                  y = NULL,
                  ymin = ss$`(Intercept)`[2] + ss$random.effect[1] + ss$log10BM[2] * bm.range,
                  ymax = ss$`(Intercept)`[3] + ss$random.effect[1] + ss$log10BM[3] * bm.range,
                  fill = "BMR"), alpha = .1) +
  geom_line(data = data.frame(),
            aes(x = bm.range,
                y = ss$`(Intercept)`[1] + ss$random.effect[1] + ss$MRFMR[2] + (ss$log10BM[1] + ss$`log10BM:MRFMR`[1]) * bm.range, col = "FMR")) +
  geom_ribbon(data = data.frame(),
              aes(x = bm.range,
                  y = NULL,
                  ymin = ss$`(Intercept)`[2] + ss$MRFMR[2] + ss$random.effect[1] + (ss$log10BM[2] + ss$`log10BM:MRFMR`[2]) * bm.range,
                  ymax = ss$`(Intercept)`[3] + ss$MRFMR[3] + ss$random.effect[1] + (ss$log10BM[3] + ss$`log10BM:MRFMR`[3]) * bm.range,
                  fill = "FMR"), alpha = .1) +
  xlab(expression(log[10]~Body~mass~(g))) +
  ylab(expression(log[10]~Metabolic~rate~(kJ/day))) +
  scale_color_discrete(name = NULL,
                       labels = c("Basal metabolic rate", "Field metabolic rate"),
                       breaks = c("BMR", "FMR")) +
  scale_fill_discrete(name = NULL,
                       labels = c("Basal metabolic rate", "Field metabolic rate"),
                       breaks = c("BMR", "FMR")) +
  theme_R() +
  theme(legend.position = c(0, 1), 
        legend.background = element_rect(linetype = "solid", colour = "black"),
        legend.justification = c(-0.5, 1.5)) +
  scale_y_continuous(limits = c(0, 6.5))
ggsave("output/appendix2_fig4_FMR_BMR_plot.png", width = 25.6, height = 14.4, units = "cm")
