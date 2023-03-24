# Comple tax resolution list

library(tidyverse)

gen.tax <- read_csv("data/Genoud2018_tax.solver.csv", col_types = cols()) %>% 
  transmute(Binomial.Source = Binomial,
            Binomial.1.2,
            Source = "Genoud.2018")

tax.nagy <- read_csv("data/Nagy1999_tax.solver.csv", col_types = cols()) %>% 
  transmute(Binomial.Source = str_replace(Species, " ", "_"),
            Binomial.1.2 = str_replace(Binomial.1.2, " ", "_"),
            Source = "Nagy.1999")

cap.tax <- read_csv("data/Capellini2010_tax.solver.csv", col_types = cols()) %>% 
  transmute(Binomial.Source = Species,
            Binomial.1.2,
            Source = "Capellini.2010")

hud.tax <- tibble(Binomial.Source = "Lama_glama",
                  Binomial.1.2 = "Lama_guanicoe",
                  Source = "Hudson.2013")

final.tax <- bind_rows(gen.tax, tax.nagy, cap.tax, hud.tax)

write_csv(final.tax, "output/taxonmy.resolution.csv")

