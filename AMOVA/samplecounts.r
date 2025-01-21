library(tidyverse)

meta <- read_table("metadata.txt")

counts <- meta %>% group_by(population, year) %>% summarize(n = n())

meta2 <- meta
meta2$newpop <-  sapply(meta2$name, function(x) {strsplit(x, "_")[[1]][1]})
meta2$newpop <- gsub("TXL", "TX",meta2$newpop)

counts <- meta2 %>% group_by(newpop, year) %>% summarize(n = n())

widecnt <- counts %>% pivot_wider(id_cols = newpop, names_from = year, values_from = n)

write.csv(widecnt, file = "passing_samplecounts.txt", row.names = FALSE, quote= FALSE, na = "0")
