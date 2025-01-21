library(tidyverse)

outliers <- as.vector(read.table("outliers.txt")$V1)
outlier_idx <- gsub("snp_", "", outliers) |> as.numeric()

markers <- read.table("out.hwe", header = TRUE) %>%  select(1, 2)
indexed_outliers <- markers[outlier_idx,]

write.table(indexed_outliers, file= "outliers_vcf.txt", row.names = F, col.names = F, quote = F)
