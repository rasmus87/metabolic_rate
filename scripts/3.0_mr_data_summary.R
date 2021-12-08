# Run imputation for all trees

# Load libraries
library(tidyverse)

# Load data
mr <- read_csv("builds/metabolic_rate_data.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.2/Data/Traits/Trait_data.csv", col_types = cols())

# Summarise sources
mr %>% count(MR.type, Source)

# Summarise species orders families
mr %>% count(Binomial.1.2)
mr %>% count(Family.1.2)
mr %>% count(Order.1.2)

mr %>% 
  group_by(MR.type) %>% 
  summarise(N.species = length(unique(Binomial.1.2)), 
            N.datapoints = n())

# Summarise species info
mr %>% group_by(MR.type) %>% 
  summarise(n.species = length(unique(Binomial.1.2)),
            n.min = min(table(Binomial.1.2)),
            n.max = max(table(Binomial.1.2)),
            n.median = median(table(Binomial.1.2)) %>% as.integer,
            log10BM.mean = mean(log10BM),
            log10BM.min = min(log10BM),
            log10BM.max = max(log10BM))

mr %>% 
  group_by(MR.type) %>% 
  count(Binomial.1.2) %>%
  count(n) %>% 
  arrange(MR.type, n) %>% 
  rename(samples.pr.species = n,
         n.species = nn)

