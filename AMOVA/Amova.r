### AMOVA ####

library(adegenet)
library(poppr)
library(pegas)
library(ade4)


setwd("~/Omega/USM PhD//Projects//Active//Yellowfin Tuna//Popgen//Analyses//nonkin/AMOVA/")
#setwd(getwd())
infile_outlier <- "../DAPC/YFT.snp.kinrm.pcrelate.outlier.gen"
infile_neutral <- "../DAPC/YFT.snp.kinrm.pcrelate.neutral.gen"

pops <- c("ATL","GOM","IVC", "SEN", "GOM", "VZ")

neut <- read.genepop(infile_neutral, ncode = 3L)
outl <- read.genepop(infile_outlier, ncode = 3L)
popNames(neut) <- pops
popNames(outl) <- pops

#### Add Strata and perform AMOVA ####

.strata <- read.delim("../DAPC/metadata.txt", stringsAsFactors = TRUE)
# rename TX pops to GOM
.strata$population <- gsub("TX", "GOM", .strata$population)
strata(neut) <- .strata
strata(outl) <- .strata

amova_results_neut <- poppr.amova(
  neut,
  hier = ~population/year,
  clonecorrect = FALSE,
  within = TRUE,
  squared = TRUE,
  correction = "quasieuclid",
  algorithm = "farthest_neighbor",
  threads = 4,
  missing = "loci",
  cutoff = 0.10,
  quiet = FALSE,
  method = "pegas",
  nperm = 50000
)

amova_results_out <- poppr.amova(
  outl,
  hier = ~population/year,
  clonecorrect = FALSE,
  within = TRUE,
  squared = TRUE,
  correction = "quasieuclid",
  algorithm = "farthest_neighbor",
  threads = 4,
  missing = "loci",
  cutoff = 0.10,
  quiet = FALSE,
  method = "pegas",
  nperm = 50000
)
