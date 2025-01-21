#! /usr/bin/env Rscript

library(stringr)
library(randomForest)
library(dplyr)
library(BEDMatrix)
setwd(getwd())

#load in ped data
ped <- read.table("../data/YFT.randomforest.neutral.ped", sep="\t", header=FALSE)
#.ped <- BEDMatrix("../data/YFT.neutral.kinless.bed")
# 401 x 11270
# isolate sample names


# rename populations
ped$V1 <- as.factor(ped$V1)
levels(ped$V1) <- c("sATL", "nATL", "sATL", "GOM", "sATL",  "nATL", "GOM",  "sATL" )
ped$V1 <- as.character(ped$V1)

# load in snp info
map <- as.data.frame(read.table("../data/YFT.neutral.kinless.bim", sep="\t", header=FALSE))
snps <- read.delim("../data/neutral_markers.txt", header = FALSE)$V1
map$V2 <- snps
POPID <- factor(ped[,1])
# 5632 x 4

# convert genotype information
f <- ped[,7:ncol(ped)]
names(f) <- unlist(lapply(map$V2, function(x) rep(x, 2)))

f[] <- lapply(f, function(x) str_replace_all(x, c("A"="1", "C"="2", "G"="3", "T"="4")))

features <- as.matrix(apply(f, MARGIN=2, FUN=function(x) as.numeric(x)))

#rewriting 1-4 as 0 and 1 within each column (2)
change <- function(x) {                    
        occur<-as.data.frame(sort(table(x),decreasing=TRUE))
        x[x==occur[1,1]]<-1
        x[x==occur[2,1]]<-0
        x
}

rewrite <- apply(features, 2, FUN=change)
paired <- rewrite[,order(colnames(rewrite))]

#combining columns into proper single feature so 0/0=0 0/1=0.5 and 1/1=1
oddindex <- c(((1:(ncol(paired)/2))*2-1))
evenindex <- c(((1:(ncol(paired)/2))*2))

even<-paired[,evenindex]
odd<-paired[,oddindex]
new_features<-(odd+even)/2

dim(new_features)
#  401 5632

BFT.rf <- randomForest(new_features, POPID, ntree=3500, nodesize=8, importance=TRUE, proximity=TRUE, do.trace=250)
bft_import <- importance(BFT.rf)
write.table(bft_import, file = "importance.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

suppressMessages(pdf("importanceplot.pdf", height = 11, width = 8.5))
varImpPlot(BFT.rf)
graphics.off()


predicted_df <- data.frame(sample = sample_names, population = as.vector(unlist(ped[1])), predicted = BFT.rf$predicted)
predicted_df$correct <- predicted_df$population == predicted_df$predicted
sum(predicted_df$correct == TRUE) / length(predicted_df$correct)
# 31.9% accurately predicted
write.table(predicted_df, file = "predicted.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

save.image(file = "randomforest.RData")