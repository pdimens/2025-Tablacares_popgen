##
#  find the sites/contigs associated with the outlier loci
##
setwd("C:/Users/pdime/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/outlier_analyses/outFLANK")
sitematrix <- read.table("../../../VCF filtering/sitenames.txt", header = TRUE)
outliers <- read.table("strict_outliers.txt")$V1
# split out the snp indices from PGDSpider's naming scheme (i.e. the number after the underscore)
indices <- as.numeric(sapply(outliers, function(x) {strsplit(x, "_")[[1]][2]}))

# index the VCF snps based on the indices
outlier_positions <- sitematrix[indices, ]

write.table(outlier_positions, file = "outlier_contig_IDs.txt", row.names = FALSE, quote = FALSE)

save.image("outlier_contig_ID.rdata")
