library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(ggplot2)
library(dplyr)
library(SeqVarTools)

setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/all_samples/kinship/king_pcrelate")
#setwd(getwd())

#radiator::genomic_converter(
#    "~/simulated_sibship.gen",
#    output = "gds"
#)

#SNPRelate::snpgdsVCF2GDS(vcf.fn = vcf_file, 
#                         out.fn = gds_file,
#                         method = 'biallelic.only',
#                         verbose = TRUE)
#

gdsfile <- "../../data/YFT_maf01.snps.gds"
gdsobj <- snpgdsOpen(gdsfile)

# pruning
snpset <- snpgdsLDpruning(gdsobj, method="corr", autosome.only = FALSE, slide.max.bp=10e6, ld.threshold=sqrt(0.1))
pruned <- unlist(snpset)
gdsfile <- "kinship.ldpruned.gds"
gdsobj <- snpgdsOpen(gdsfile)
# KING PC-based kinship
king <- snpgdsIBDKING(gdsobj, autosome.only = FALSE)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
kinship <- snpgdsIBDSelection(king)

ggplot(kinship, aes(IBS0, kinship)) +
    geom_point(alpha=0.5, color = "dodgerblue") +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    ylab("kinship estimate")

#ggsave("king_kinship.png")



# GENESIS PC-AIR
# based on KING values
# thresholds set for "unrelated is less than first cousins"
sampset <- pcairPartition(
    kinobj=kingMat, kin.thresh=2^(-9/2),
    divobj=kingMat, div.thresh=-2^(-9/2)
)

# pca on unrelated individuals
pca.unrel <- snpgdsPCA(gdsobj, sample.id=sampset$unrels, autosome.only = FALSE)
# project values for relatives
snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj=gdsobj)
samp.load <- snpgdsPCASampLoading(snp.load, gdsobj=gdsobj, sample.id=sampset$rels)

# combine unrelated and related PCs and order as in GDS file
pcs <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(pcs) <- c(pca.unrel$sample.id, samp.load$sample.id)
sample.id <- rownames(kingMat)
samp.ord <- match(sample.id, rownames(pcs))
pcs <- pcs[samp.ord,]

# add population information
annot <- read.table("../../data/popmap", header = FALSE)
pc.df <- as.data.frame(pcs)
names(pc.df) <- 1:ncol(pcs)
pc.df$sample.id <- row.names(pcs)
names(annot) <- c("sample.id", "population")
pc.df <- left_join(pc.df, annot, by="sample.id")

library(GGally)
library(RColorBrewer)
popn <- length(unique(annot$population))
pop.cols <- setNames(brewer.pal(popn, "Paired"), unique(annot$population))
ggparcoord(pc.df, columns=1:12, groupColumn="population", scale="uniminmax") +
    scale_color_manual(values=pop.cols) +
    xlab("PC") + ylab("")

# 3 seems to separate things well
#ggsave("parallel_coords.png")


# Finally, pcrelate
showfile.gds(closeall=TRUE)
geno <- GdsGenotypeReader(filename = gdsfile)
genodata <- GenotypeData(geno)
genodata <- GenotypeBlockIterator(genodata)

pcrel <- pcrelate(genodata, pcs=pcs[,1:3], training.set=sampset$unrels, 
                  sample.include=sample.id)

pcrelMat <- pcrelateToMatrix(pcrel, scaleKin = 1)

pca <- pcair(genodata,
             kinobj=pcrelMat, kin.thresh=2^(-9/2),
             divobj=kingMat, div.thresh=-2^(-9/2),
             sample.include=sample.id, 
             autosome.only = FALSE)


pcs <- pca$vectors
pc.df <- as.data.frame(pcs)
names(pc.df) <- paste0("PC", 1:ncol(pcs))
pc.df$sample.id <- row.names(pcs)
pc.df <- left_join(pc.df, annot, by="sample.id")

ggplot(pc.df, aes(PC1, PC2, color=population)) + geom_point() +
    scale_color_manual(values=pop.cols)

#ggsave("post_pcrelate.png")

# perform pcrelate again
pcrel <- pcrelate(genodata, pcs=pcs[,1:3], training.set=pca$unrels, 
                  sample.include=sample.id)

kinship <- pcrel$kinBtwn
pcrelate_filter <- kinship %>% 
    mutate(
        relationship = case_when(
            (kin <= 0.1767) ~ "unrelated",
            (kin >=0.1767 & kin < 0.3535) ~ "halfsib",
            (kin >= 0.3535) ~ "fullsib"
        )
    )
write.table(pcrelate_filter, file = "allkinship_pcrelate.txt", row.names = FALSE, quote = FALSE)
pcrelate_filter <- read.table("allkinship_pcrelate.txt", header = T)
mycolors <- c("#bbbbbb", "dodgerblue", "#8a556e")

ggplot(pcrelate_filter, aes(k0, kin)) +
    geom_point(alpha=0.5, aes(color = relationship)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    scale_color_manual(values = mycolors) +
    labs(title = "PCRelate pairwise kinship estimates") +
    ylab("kinship estimate") +
    xlab("probability of sharing 0 alleles") +
    theme_bw()

ggsave("kinship_pcarelate.png", height = 8, width = 10, units = "in")



kintable <- pcrelate_filter %>% 
                filter(relationship != "unrelated") %>% 
                arrange(desc(kin))

write.table(kintable, file = "kin_pcrelate.txt", row.names = FALSE, quote = FALSE)
