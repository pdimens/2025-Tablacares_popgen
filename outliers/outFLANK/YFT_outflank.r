## Load packages
library("dartR")
library("ggplot2")
library("dplyr")
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin_haplo/outFLANK")
# ----------------- #
#
# OutFLANK
#
# ----------------- #

#source("outflank.rdata")

##### Outflank and Bayescan results ####
# xx haplotypes
bayescan_idx <- c()
# xx haplotypes
outflank_idx <- c()

#### Running outflank #### 
#infile <- "../data/YFT.pcrelate.kinrm.haplo.gen"
infile <- "../data/YFT.snp.kinrm.pcrelate.gen"


## Import data
full_dataset <- import2genind(infile, ncode = 3L)

pop_names <- c("ATL", "GOM", "IVC", "GOM", "SEN", "TX", "VZ")


## Change group labels
popNames(full_dataset) <- pop_names
  
## Run OutFLANK using dartR wrapper script
full_outflnk <- gl.outflank(
                  full_dataset,
                  Hmin= 0.1,
                  qthreshold = 0.01,
                  LeftTrimFraction = 0.05,
                  RightTrimFraction = 0.1
                )



#### Processing results ####
## Outliers
full_outflnk <- full_outflnk$outflank$results

## Remove duplicated rows for each SNP
toRemove <- seq(1, nrow(full_outflnk), by=2)
full_outflnk <- full_outflnk[-toRemove, ]
#summary(full_outflnk)

# check for linear conformity, remove loci that deviate substantially
plot(full_outflnk$FST, full_outflnk$FSTNoCorr)

#full_outflnk$putative <- FALSE
full_outflnk <- full_outflnk %>% 
  mutate(putative = 
    case_when(
      (He >= 0.1 & OutlierFlag == TRUE) ~ TRUE,
      (He >= 0.1 & OutlierFlag == FALSE) ~ FALSE,
      (He < 0.1) ~ FALSE
    )
  )

full_outflnk %>% 
ggplot(aes(x = He, y = FST)) + 
  geom_point(aes(color = putative)) +
  scale_color_manual(values = c("grey49", "dodgerblue")) +
  labs(title = "outFlank outliers" , x = "Heterozygosity", y = "FST", color = "Detection")

ggsave(filename = "outflank_outliers.kinless.snps.png")

# pull out outliers w/ He > 0.1
putatives <- full_outflnk %>% filter(He >= 0.1 & OutlierFlag == TRUE)
write.csv(putatives, file = "outFLANK_kinless_putatives.csv", row.names = FALSE, quote = FALSE)


## Get indexes for outliers
out_index <- which(full_outflnk$OutlierFlag==TRUE)
outflank_names <- locNames(full_dataset)[out_index]
outlier_sub <- full_outflnk[out_index,]
bayes_sub <- full_outflnk[bayescan_idx,]

#### Write output files ####
write.csv(full_outflnk, file = "outFLANK_kinless_FST.csv", row.names = FALSE, quote = FALSE)

full_out.outflnk <- locNames(full_dataset)[out_index]
full_out.outflnk

# output only outlier loci
write.csv(full_out.outflnk, file="outlier_kinless_OutFLANK.csv", row.names = FALSE, quote = FALSE, col.names = FALSE)

#### Plots ####
## qvalue hist ##
hist(full_outflnk$qvalues, breaks = 100, 
     xlab = "Q values", main = "Q Value distribution of SNPs")
plot(x = full_outflnk$qvalues)


## fancy ggplot ##
full_outflnk$outlier <- "Neither"
outlier_both <- intersect(bayescan_idx, out_index)
full_outflnk$outlier[out_index] <- "outFLANK"
bayescan_safe <- bayescan_idx[bayescan_idx <= length(full_outflnk$LocusName)]
full_outflnk$outlier[bayescan_safe] <- "BayeScan"
full_outflnk$outlier[outlier_both] <- "Both"
full_outflnk$outlier <- factor(full_outflnk$outlier, levels = c("Neither", "outFLANK", "BayeScan", "Both"), ordered = TRUE)

mycolors <- c("#bbbbbb", "#4095b5", "#8a556e", "#f4cf30")
ggplot(data = full_outflnk, x = He, y = FST) + 
  theme_classic() +
  geom_point(aes(x = He, y = FST, col = outlier), alpha = 0.8,  size = 2.1) + 
  geom_point(data = subset(full_outflnk, outlier %in% "BayeScan"), aes(x = He, y = FST), color = "#8a556e") +
  geom_vline(xintercept = 0.1, alpha = 0.8, linetype = "dashed", size = 0.3) +
  labs(x = "Heterozygosity", y = "FST", color = "Detection") +
  scale_color_manual(values = mycolors) +
  ggtitle("Outlier Loci") +
  theme(plot.title = element_text(hjust = 0.5))


bayescan_out <- locNames(full_dataset)[bayescan_idx]
outflank_out <- locNames(full_dataset)[out_index]

write.table(outflank_out, file = "outflank_outliers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("outflank.rdata")
