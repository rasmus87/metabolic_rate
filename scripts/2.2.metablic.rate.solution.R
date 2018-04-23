library(tidyverse)

# solution <- read_csv("builds/metabolic.rate_fit.solution.csv")
solution <- read_csv("builds/test_metabolic.rate_fit.solution.csv")

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
solution.summary

# In real values how much larger is FMR than BMR
round(10^solution.summary[3,], 2)

# Plot chains
solution.gathered <- solution %>% 
  select(-random.effect) %>% 
  mutate(id = paste0(tree, ".", chain)) %>% 
  select(-tree, -chain)
solution.gathered <- gather(solution.gathered, key = "variable", value = "value", -id)
ggplot(solution.gathered, aes(x = value, group = id)) +
  geom_line(stat = "density", col = alpha("blue", .1)) +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "") +
  geom_line(stat = "density", aes(group = NULL), lwd = 1)

# Plot data and regression line:
# Load data
mr <- read_csv("builds/mr.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

# Filter for terrestial
terrestrial <- mam %>% filter(Terrestrial == 1) %>% pull(Binomial.1.2)
mr <- mr %>% 
  filter(Binomial.1.2 %in% terrestrial)

mr.all <- mr  
mr <- mr %>% group_by(Binomial.1.2, MR) %>% summarise_at(c("log10BM", "log10MR"), mean)

ss <- as_data_frame(t(solution.summary))

mam.mr <- read_csv("builds/imputed.metabolic.rate.csv")
bm.range <- range(mam.mr$log10BM)

ggplot(mr.all, aes(x = log10BM, y = log10MR)) +
  geom_point(col = "grey", alpha = .5) +
  geom_point(data = mr, aes(col = MR, fill = MR), alpha = .75) + 
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
                  fill = "FMR"), alpha = .1)

ggplot(mam.mr, aes(x = log10BM)) +
  theme_bw() +
  xlab("log10 Body mass (g)") +
  ylab("log10 Metabolic rate (kJ / day)") +
  scale_color_discrete(c("FMR", "BMR")) +
  # geom_point(data = mr.all, aes(x = log10BM, y = log10MR), col = "grey") +
  # geom_point(data = mr, aes(x = log10BM, y = log10MR, col = MR)) + 
  geom_point(aes(y = log10FMR_est, col = "FMR"), pch = 1) +
  geom_point(aes(y = log10BMR_est, col = "BMR"), pch = 1) +
  geom_line(aes(y = ss$`(Intercept)`[1] + ss$random.effect[1] + ss$log10BM[1] * log10BM, col = "BMR")) +
  geom_ribbon(aes(ymin = ss$`(Intercept)`[2] + ss$random.effect[1] + ss$log10BM[2] * log10BM,
                  ymax = ss$`(Intercept)`[3] + ss$random.effect[1] + ss$log10BM[3] * log10BM,
                  fill = "BMR"), alpha = .1) +
  geom_line(aes(y = ss$`(Intercept)`[1] + ss$MRFMR[2] + ss$random.effect[1] + (ss$log10BM[1] + ss$`log10BM:MRFMR`[1]) * log10BM, col = "FMR")) +
  geom_ribbon(aes(ymin = ss$`(Intercept)`[2] + ss$MRFMR[2] + ss$random.effect[1] + (ss$log10BM[2] + ss$`log10BM:MRFMR`[2]) * log10BM,
                  ymax = ss$`(Intercept)`[3] + ss$MRFMR[3] + ss$random.effect[1] + (ss$log10BM[3] + ss$`log10BM:MRFMR`[3]) * log10BM,
                  fill = "FMR"), alpha = .1)


