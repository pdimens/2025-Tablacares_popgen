library(tidyverse)

setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/all_samples/kinship/king_pcrelate/simulated_ci")

kinship <- read.table("PCRelate_sims.txt") %>% filter(grepl("off1", ID1) & grepl("off2", ID2))

simsplit <- function(x){
  .tmp <- sapply(x, function(y){strsplit(y, "_")[[1]][1]})
  as.character(unlist(.tmp))
}

fs <- kinship %>% 
  filter(grepl("_fs", ID1) & grepl("_fs", ID2)) %>% 
  filter(simsplit(ID1) == simsplit(ID2)) %>% 
  mutate(relationship = "fullsib")

hs <- kinship %>% 
  filter(grepl("_hs", ID1) & grepl("_hs", ID2)) %>% 
  filter(simsplit(ID1) == simsplit(ID2)) %>% 
  mutate(relationship = "halfsib")

un <- kinship %>% 
  filter(grepl("_un", ID1) & grepl("_un", ID2)) %>% 
  filter(simsplit(ID1) == simsplit(ID2)) %>% 
  mutate(relationship = "unrelated")

kinship_merged <- rbind(fs, hs, un)
write.table(kinship_merged, file = "simulated_pcrelate.txt", quote = FALSE)

kinship_merged <- read.table("simulated_pcrelate.txt", header = T)
mycolors <- c("#8a556e", "#f4cf30", "#bbbbbb")
mycolors <- c("#bbbbbb", "dodgerblue", "#8a556e")
kinship_merged %>% 
  ggplot(aes(x=k0, y = kin)) +
  geom_point(aes(color = relationship), alpha = 0.5) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  scale_color_manual(values = mycolors) +
  labs(title = "PCRelate kinship estimates for simulated sibship pairs") +
  ylab("kinship estimate") +
  xlab("Probability of sharing 0 alleles") +
  theme_bw()

ggsave('PCRelate_sims.png', height = 8, width = 10, units = "in")
