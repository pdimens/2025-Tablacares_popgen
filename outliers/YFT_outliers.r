## Load packages
library("ggplot2")
library("dplyr")
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/outliers")

### outflank dataframe ####
all_loci <- read.csv("outFLANK/outFLANK_kinless_FST.csv")%>% 
  select(LocusName, He, FST, FSTNoCorr, OutlierFlag) %>% 
  mutate(LocusName = gsub("\\.[0-9]+", "", LocusName))

all_outflank <- filter(all_loci, OutlierFlag == TRUE)
all_outflank_putatives <- all_outflank$LocusName

##### Outflank and Bayescan results ####
# snps
out_bayescan <- readLines("bayescan/bayescan_outliers.txt")
bayescan_idx <- c()
# snps
out_outflank <- read.csv("outFLANK/outFLANK_kinless_putatives.csv")$LocusName
out_outflank <- gsub("\\.[0-9]+", "", out_outflank)
outflank_idx <- c()

# all shared putative outliers regardless of He <- used in graphic
all_shared <- intersect(out_bayescan, all_outflank_putatives)
all_onlybayes <- setdiff(out_bayescan, all_shared)
all_onlyoutflnk <- setdiff(all_outflank_putatives, all_shared)
all_nonout <- setdiff(all_loci$LocusName, c(all_onlybayes, all_onlyoutflnk, all_shared))

# only shared putative outliers (He threshold) <- used in publication
out_shared <- intersect(out_bayescan, out_outflank)
out_only_baye <- setdiff(out_bayescan, out_shared)
out_only_outf <- setdiff(out_outflank, out_shared)
non_out <- setdiff(all_loci$LocusName, c(out_bayescan, out_outflank, out_shared))

## fancy ggplot ##
all_loci <-all_loci %>% 
  mutate(
    outlier = case_when(
      (LocusName %in% all_onlybayes) ~ "bayescan",
      (LocusName %in% all_onlyoutflnk) ~ "outflank",
      (LocusName %in% all_shared) ~ "bayescan + outflank",
      (LocusName %in% all_nonout ) ~ "neutral"
    )
  )


mycolors <- c("#8a556e", "#bbbbbb", "#f4cf30")

all_loci %>% 
  ggplot(x = He, y = FST) +
  geom_point(aes(x = He, y = FST, col = outlier), alpha = 0.8,  size = 2.1) + 
  geom_vline(xintercept = 0.1, alpha = 0.8, linetype = "dashed", size = 0.3) +
  labs(title = "Putative Outlier Loci", x = "Heterozygosity", y = "FST", color = "Outlier Detection",
       subtitle = "All Bayescan-detected loci were detected by outFLANK") +
  scale_color_manual(values = mycolors)
ggsave("outlierplots.png", height = 6, width = 12, units = "in")
save.image(file = "outliers.rdata")
write.table(all_loci, file = "bscan_outflank_outliers.txt", row.names = FALSE, quote = FALSE)
