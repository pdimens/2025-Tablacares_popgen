library(adegenet)
library(poppr)
library(hierfstat)

setwd("/mnt/Win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/data/")

data <- read.genepop("YFT.snp.kinrm.pcrelate.gen", ncode = 3L)

popNames(data) <- c("wATL","GOM","IVC","GOM","SEN","TX", "VZ")

WC_lxl <- data.frame(pop1 = character(), pop2 = character(), locus = character(), FST = numeric())

pnames <- popNames(data)
for(p1 in 1:length(pnames)-1) {
    pn1 <- pnames[p1]
    for(p2 in 2:length(pnames)) {
        pn2 <- pnames[p2]
        tmpgen <- genind2hierfstat(popsub(data, sublist = c(pn1, pn2)))
        tmpwc <- wc(tmpgen)$per.loc$FST
        WC_lxl <- rbind(WC_lxl, data.frame(pop1 = p1, pop2= p2, locus = locNames(data), FST = tmpwc))
    }
}

write.table(WC_lxl, file = "~/WC_lxl.txt", row.names = F, quote = F)

WC_lxl_global <- wc(genind2hierfstat(data))$per.loc$FST
write.table(data.frame(locus = locNames(data), FST = WC_lxl_global), file = "~/WC_lxl.global.txt", row.names = F, quote = F)
