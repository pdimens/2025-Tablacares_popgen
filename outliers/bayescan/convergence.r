library(coda)
library(adegenet)
source("plot_R.r")
setwd("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/bayescan")
infile <- "yft_bscan_fst.txt"
yft_bscan_fst <- read.csv(infile, sep="")
genfile <- "../../data_maf_01/BFT_final80_maf01.gen"
bft <- read.genepop(genfile, ncode = 3L)

# convergence checks
chain <- read.table("yft_bscan.sel", header = TRUE)
chain <- chain[-c(1)]
mc_chain <- mcmc(chain, thin = 10)
plot(mc_chain)
summary(mc_chain)
autocorr.diag(mc_chain)
effectiveSize(mc_chain)
geweke.diag(mc_chain, frac1 = 0.1, frac2 = 0.5)
heidel.diag(mc_chain, eps=0.1, pvalue = 0.05)


plt <- plot_bayescan(infile,0,FDR=0.05)
out_ <- yft_bscan_fst[yft_bscan_fst$qval<=0.05,]

bayescan_idx <- plt$outliers

locinames <- readLines("../outFLANK/YFT.snp.kinrm.pcrelate.gen", (dim(yft_bscan_fst)[1] + 1))
locinames <- locinames[2:length(locinames)]

bayescan_out <- locinames[bayescan_idx]
write.csv(bayescan_out, "bayescan_outliers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("bayescan.convergence.rdata")
