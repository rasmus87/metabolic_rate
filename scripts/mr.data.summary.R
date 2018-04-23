library(tidyverse)

# Load data
mr <- read_csv("builds/mr.csv", col_types = cols())
mam <- read_csv("../PHYLACINE_1.1/Data/Traits/Trait_data.csv", col_types = cols())

# Filter for terrestial
terrestrial <- mam %>% filter(Terrestrial == 1) %>% pull(Binomial.1.2)
mr <- mr %>% 
  filter(Binomial.1.2 %in% terrestrial)

# Summarise sources
mr %>% count(MR, source)

# Summarise species info
mr %>% group_by(MR) %>% 
  summarise(n.species = length(unique(Binomial.1.2)),
            n.min = min(table(Binomial.1.2)),
            n.max = max(table(Binomial.1.2)),
            n.median = median(table(Binomial.1.2)) %>% as.integer,
            log10BM.mean = mean(log10BM),
            log10BM.min = min(log10BM),
            log10BM.max = max(log10BM))

mr %>% 
  group_by(MR) %>% 
  count(Binomial.1.2) %>%
  count(n) %>% 
  arrange(MR, n) %>% 
  rename(samples.pr.species = n,
         n.species = nn)
