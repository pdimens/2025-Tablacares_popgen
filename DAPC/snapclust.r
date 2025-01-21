#! /usr/bin/env Rscript

library("adegenet")
library(tidyverse)
#setwd(getwd())
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/DAPC")

##### Load in data and create backup object
i <- "neutral"
infile <- paste0("YFT.snp.kinrm.pcrelate.", i, ".gen")
data <- read.genepop(infile, ncode = 3L)
#backup <- data
PopNames <- c("ATL","GOM","IVC", "SEN", "TX", "VZ")
popNames(data) <- PopNames
metadata <- read_tsv("metadata.txt", col_names = TRUE)


##### Cross validation
set.seed(6969)
loc_year <- paste(metadata$population, metadata$year, sep = "_")
par(mfrow = c(1,3))
# by year
by_year <- dapc(data, pop = all_meta$year, n.pca = 300, n.da = 6)
scatter(by_year)
by_year$grp <- as.factor(all_meta$population)
scatter(by_year)
by_year$grp <- as.factor(loc_year)
scatter(by_year)

# by location
by_location <- dapc(data, pop = all_meta$population, n.pca = 300, n.da = 5)
scatter(by_location)
by_location$grp <- as.factor(all_meta$year)
scatter(by_location)
by_location$grp <- as.factor(loc_year)
scatter(by_location)

# by location-year
by_locyear <- dapc(data, pop = loc_year, n.pca = 300, n.da = 5)
scatter(by_locyear)
by_locyear$grp <- as.factor(metadata$year)
scatter(by_locyear)
by_locyear$grp <- as.factor(metadata$population)
scatter(by_locyear)

par(mfrow = c(1,1))
location_clust <- snapclust(data, k = 4, dim.ini = 298)
compoplot(location_clust)
