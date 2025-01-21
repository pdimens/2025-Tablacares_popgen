## Figure out which SNPs are shared between the loci used to
## construct the linkage map, and those used in the popgen study

library(tidyverse)
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/Ne")

# make the linkage map all nice
linkmap <- read_delim(
  "../../../../Linkage Map/LepAnchor/miss90_both/miss90_both.map",
  col_names = FALSE
) %>%
  select(1,2,5,6) %>% 
  rename("lg" = X1, "position" = X2, "male" = X5, "female" = X6) %>% 
  arrange(lg, position) %>% 
  rowwise() %>% mutate(cm = round(mean(c(male, female)), digits = 4)) %>% 
  mutate(male = round(male, digits = 4), female = round(female, digits = 4))

# save it for forever
#write.table(linkmap, file = "sexaverage.map", row.names = FALSE, quote = FALSE)

linkmap <- read_delim("sexaverage.map")
snp_pos <- read_delim("../data/scaffold_positions.snp")
# get the index of the snp positions for later filtering
snp_pos$idx <- 1:nrow(snp_pos)
# subset the positions to extract only snps on scaffolds
scaffold_pos <- snp_pos[grepl("LG", snp_pos$scaffold),]
names(scaffold_pos) <- c("lg", "position", "idx")

recomb_rate <- linkmap %>%  
  mutate(position = position / 1e6) %>%
  group_by(lg) %>%
  mutate(recomb = max(cm) / (max(position) - min(position))) %>% 
  select(lg, recomb) %>% unique

scaff_conversion <- merge(scaffold_pos, recomb_rate, by = "lg")

scaff_conversion <- scaff_conversion %>% 
  mutate(cm = position / 1e6 * recomb)

# used to filter the genepop via PopGen.jl
write.table(scaff_conversion, file = "conversion_to_cm.snp", row.names = FALSE, quote = FALSE)

# add snp names
snps <- read.table("names.snp", header = TRUE)
linkne_file <- scaff_conversion %>% 
  mutate(locus = snps$loci, chromosome = as.numeric(gsub("LG", "",lg))) %>% 
  select(locus, chromosome, cm) %>% 
  mutate(cm = round(cm, digits = 4))

write.table(linkne_file, file = "linkne.map", row.names = FALSE, quote= FALSE)
