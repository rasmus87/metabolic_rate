library(tidyverse)

solution <- read_csv("builds/metabolic.rate_fit.solution.csv")

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
solution.summary <- cbind(mean = colMeans(solution[, 1:4]), HPDinterval.mcmc(solution[, 1:4]))
solution.summary

# In real values how much larger is FMR than BMR
round(10^solution.summary[3,], 2)

solution.gathered <- solution %>% 
  mutate(id = paste0(tree, ".", chain)) %>% 
  select(-tree, -chain)
solution.gathered <- gather(solution.gathered, key = "variable", value = "value", -id)
ggplot(solution.gathered, aes(x = value, group = id)) +
  geom_line(stat = "density", col = alpha("blue", .1)) +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "") +
  geom_line(stat = "density", aes(group = NULL), lwd = 1)
