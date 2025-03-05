library(tidyverse) |> suppressPackageStartupMessages()
library(circlize) |> suppressPackageStartupMessages()

setwd("/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/all_samples/kinship/king_pcrelate/")
#### circos pairwise ####
infile <- "/mnt/win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/all_samples/kinship/king_pcrelate/kin_pcrelate.txt"
kinship = read.table(infile, header = T) %>% 
  select("ID1","ID2", "relationship")

set.seed(6969)
samplenames <- read.table("samplenames.txt") %>% rename( "name" = V1) %>% 
  mutate(population = gsub("_[0-9]+", "" ,name), sample = sample(seq_along(name))) %>% 
  mutate(population = gsub("TXL", "TX", population))


## rearrange samplenames ##
smps <- samplenames[samplenames$population == "ATL",]
poporder <- c("GOM", "MSAL","LA","TX","VZ","IVC","SEN")
for (i in poporder){
  smps <- rbind(smps, samplenames[samplenames$population == i,])
}

kintable <- kinship
kintable$pop1 <- "1"
kintable$pop2 <- "1"
kintable$n1 <- 1
kintable$n2 <- 2

for (i in 1:dim(kintable)[1]){
  .row <- smps[smps$name == kintable$ID1[i],]
  kintable$pop1[i] <- .row[2]
  kintable$n1[i] <- .row[3]
  .row <- smps[smps$name == kintable$ID2[i],]
  kintable$pop2[i] <- .row[2]
  kintable$n2[i] <- .row[3]
}

samplenames$population <- factor(
    samplenames$population, 
    levels = c("ATL","GOM", "MSAL","LA","TX","VZ","IVC","SEN"),
    ordered = TRUE)

smps <- samplenames

fs <- kintable[kintable$relationship == "fullsib",]
fs$pop1 <- unlist(fs$pop1)
fs$pop2 <- unlist(fs$pop2)
fs$n1 <- unlist(fs$n1)
fs$n2 <- unlist(fs$n2)
hs <- kintable[kintable$relationship == "halfsib",]
hs$pop1 <- unlist(hs$pop1)
hs$pop2 <- unlist(hs$pop2)
hs$n1 <- unlist(hs$n1)
hs$n2 <- unlist(hs$n2)

#png(file="kinship.circos.png", width=1000, height=1000)
pal <- RColorBrewer::brewer.pal(n=3, "Dark2")
pal <- c(rep(pal[1], 6), rep(pal[2], 2))
circos.clear()
circos.par("track.height" = 0.1, gap.after = c(rep(3,8)))
circos.initialize(smps$population, x = smps$sample)
circos.track(smps$population, y = smps$sample, bg.col = pal, 
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(3), 
                           CELL_META$sector.index,
                           cex = 1.1, niceFacing = TRUE)
               circos.axis(labels.cex = 0.6, major.tick = F, minor.ticks = 0, labels = FALSE)
             }, bg.border = NA)

for(i in 1:51){
  circos.link(fs$pop1[i], fs$n1[i], fs$pop2[i], fs$n2[i], col = "dodgerblue", h = 0.7, lwd = 1)
}

for(i in 1:29){
  circos.link(hs$pop1[i], hs$n1[i], hs$pop2[i], hs$n2[i], col = "grey", lwd = 2)
}

#### kinship results scatterplot ####
kinship <- read.table("allkinship_pcrelate.txt", header = T)
mycolors <- c("#8a556e", "dodgerblue", "#bbbbbb")

kinplot <- ggplot(kinship, aes(k0, kin)) +
  geom_point(alpha=0.5, aes(color = relationship)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  scale_color_manual(values = mycolors) +
  labs(title = "PCRelate pairwise kinship estimates") +
  ylab("kinship estimate") +
  xlab("probability of sharing 0 alleles")


simkin <- read.table("simulated_ci/simulated_pcrelate.txt", header = T)

mycolors <- c("#8a556e", "#f4cf30", "#bbbbbb")

simkinplot <- simkin %>% 
  ggplot(aes(x=k0, y = kin)) +
  geom_point(alpha = 0.5, aes(color = relationship)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  scale_color_manual(values = mycolors) +
  labs(title = "PCRelate kinship estimates for simulated sibship pairs") +
  ylab("kinship estimate") +
  xlab("Probability of sharing 0 alleles")


#### kinship network graph ####
library(network)
library(sna)
library(igraph)
library(intergraph)
library(ggnet)
library(ggpubr)


### number of sibships per sample ###
infile <- "file:///mnt/Win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/all_samples/kinship/king_pcrelate/kin_pcrelate.txt"
tmp = read.table(infile, header = T) %>% 
  select("ID1","ID2", "relationship")

# get counts
.tmp <- tmp %>% select(2,3) %>% rename("ID1" = ID2)
counts <- rbind(select(tmp, 1,3), .tmp) %>% table %>% as.data.frame.matrix
counts$ID <- rownames(counts)
counts <- pivot_longer(counts, -ID, names_to = "relationship", values_to = "count")
# stacked barplot
pal <- c("grey70", "dodgerblue")
brplt <- ggplot(counts, aes(fill=relationship, y=count, x=ID)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = pal) +
  xlab("Sample ID") +
  ylab("Count") +
  labs(title = "Occurrence in sibling pairs") +
  coord_cartesian(ylim = c(0,12), expand = FALSE, clip = "off") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank())
ggsave(filename = "pcrelate.occurrence.png", height = 8, width = 10, units = "in")


kingraph <- graph.edgelist(as.matrix(tmp[,c(1,2)], directed=FALSE))
relatn <- tmp$relationship
relatn[relatn == "fullsib"] <- 1
relatn[relatn == "halfsib"] <- 2
relatn <- as.numeric(relatn)
E(kingraph)$weight <- relatn
kingraph <- as.undirected(kingraph)
mtx <- as_adjacency_matrix(kingraph, attr = "weight") %>% as.matrix
net <- as.network.matrix(mtx, matrix.type = "adjacency", directed = FALSE)
nodes <- as.edgelist(net)
net %v% "relationship" <- tmp$relationship
pops <- c()
for(i in 1:23){
  pops <- c(pops, net$val[i][[1]]$vertex.names)  
}
pops <- gsub("_\\d+", "", pops)
pops <- gsub("TXL", "TX", pops)
net %v% "population" <- pops

nwrk <- ggnet2(net, 
       color = "population", 
       palette = "Set2", 
       edge.color = pal[relatn],
       legend.size = 12, 
       size = 6, 
       edge.size = 0.5, 
       legend.position = "right",
       mode = "kamadakawai",
       color.legend = "location") +
  labs(title = "     Network of Sibling Pairs") +
  geom_point(aes(color = color), size = 6, color = "white") +
  geom_point(aes(color = color), size = 6, alpha = 0.5) +
  geom_point(aes(color = color), size = 4)

tops <- ggarrange(
    kinplot,
    simkinplot,
    ncol = 2, labels = c("A", "B")
)
  
bottoms <- ggarrange(  
  brplt, nwrk,
  ncol = 2,
  labels = c("C", "D")
)

ggarrange(tops, bottoms, nrow = 2)

ggsave("kinship.multiplot.png", width = 2048, height = 1080, units = "px")
