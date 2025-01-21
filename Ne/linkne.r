library(tidyverse)

linkne <- read_table("separate_pops/linkne.results")

linkne <- linkne %>% pivot_longer(-Generations, names_to = "Population", values_to = "NE") %>% 
  mutate(NE = raster::clamp(NE, lower = 0, upper = Inf))

linkne %>% 
  ggplot(aes(x = Generations, y = NE)) +
  geom_line(color = "grey70") +
  geom_point(aes(color = Population))

linkne %>% filter(Generations == 1) %>% 
  ggplot(aes(x = Population, y = NE)) +
  geom_point(color = "dodgerblue", size = 4)

