#! /usr/bin/env Rscript

library(adegenet)
library(tidyverse) |> suppressPackageStartupMessages()
library(ggplot2)
library(ggpubr)
setwd(getwd())
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/DAPC")


##### Load in data and create backup object
i <- "neutral"
infile <- paste0("YFT.snp.kinrm.pcrelate.", i, ".gen")
data <- read.genepop(infile, ncode = 3L)
#backup <- data
PopNames <- c("ATL","GOM","IVC", "SEN", "TX", "VZ")
popNames(data) <- PopNames
#metadata <- read_delim("metadata.txt", delim = " ", col_names = TRUE)
metadata <- read_table("metadata.txt", col_names = TRUE)

#meta_arranged <- merge(data.frame(name = indNames(data)), metadata , by = "name", sort = FALSE)
#write.table(meta_arranged, file = "sample_metadata.txt")

### kmeans clustering ###

maxK <- 7
myMat <- matrix(nrow=50, ncol=maxK)
colnames(myMat) <- 1:maxK
for(i in 1:nrow(myMat)){
  grp <- find.clusters(data, n.pca = 298, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
  #print(grp$Kstat)
}

k_df <- reshape2::melt(myMat)
colnames(k_df)[1:3] <- c("Group", "K", "BIC")
k_df$K <- as.factor(k_df$K)

k_plot <- k_df %>% 
  ggplot(aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "K-means Clustering") +
  ylab("Bayesian Information Criterion (BIC)") +
  xlab("Number of groups (K)")

# prefered k=4

##### Cross validation
set.seed(6969)
#pdf(paste0("outlier_cross_validation_", i, ".pdf"))
#cat(paste0("Performing cross validation on ", i))
pramx <- xvalDapc(tab(data, NA.method="mean"), pop(data), training.set = 0.6, n.da = 4, n.pca=200:300, parallel = "multicore", ncpus = 4)
#scatter(pramx$DAPC, posi.da="bottomright", bg="white", pch=17:22)
#save.image(paste0("dapc_kinless.", i, ".rdata"))

#loc_year <- paste(metadata$population, metadata$year, sep = "_")

#par(mfrow = c(1, 3))
#####

##### by year ####
#pop(data) <- metadata$year
#by_year <- dapc(data, pop = metadata$year, perc.pca = 85, n.da = 6)
#scatter(by_year)
#by_year$grp <- as.factor(metadata$population)
#scatter(by_year)
#by_year$grp <- as.factor(loc_year)
#scatter(by_year)

##### by location ####
by_location <- dapc(data, pop = metadata$population, perc.pca = 85, n.da = 4)
scatter(by_location)
#by_location$grp <- as.factor(metadata$year)
#scatter(by_location)
#by_location$grp <- as.factor(loc_year)
#scatter(by_location)

##### by location-year ####
#by_locyear <- dapc(data, pop = loc_year, perc.pca = 85, n.da = 4)
#scatter(by_locyear)
#by_location$grp <- as.factor(metadata$year)
#scatter(by_locyear)
#by_location$grp <- as.factor(loc_year)
#scatter(by_locyear)
#scatter(by_locyear, 1, 3)

##### 2D ggplots ######
my_pal <- RColorBrewer::brewer.pal(n=6, name = "Dark2")

dapc_l <- by_location
my_df <- as.data.frame(dapc_l$ind.coord)
my_df$Group <- dapc_l$grp

df <- as.data.frame(by_location$ind.coord)
df <- cbind(metadata$name, metadata$population, df)
names(df) <- c("name", "population", "LD1", "LD2", "LD3", "LD4")
df$population <- df$population <- factor(df$population, levels = c("GOM", "TX", "VZ", "ATL", "SEN", "IVC"), ordered = TRUE)


dapc_plot <- df %>% 
  ggplot(aes(x = LD1, y = LD2, color = population, fill = population)) +
  geom_point(size = 4, shape = 21) +
  stat_ellipse(show.legend = FALSE) +
  labs(title = "DAPC (neutral markers)") +
  xlab("Linear Discriminant 1") +
  ylab("Linear Discriminant 2") +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) +
  scale_fill_manual(values=c(paste(my_pal, "96", sep = "")))


## structure plots ###

tmp <- as.data.frame(dapc_l$posterior)
tmp <- cbind(metadata$name, metadata$population, tmp)[,1:8] #drop redundant col
names(tmp) <- c("name", "origin", "ATL","GOM","IVC","SEN","TX","VZ")
tmp <- pivot_longer(tmp, c(-name, -origin), names_to = "population", values_to = "posterior")
tmp$population <- factor(tmp$population, levels = c("GOM", "TX", "VZ", "ATL", "SEN", "IVC"), ordered = TRUE)
tmp$origin <- factor(tmp$origin, levels = c("GOM", "TX", "VZ", "ATL", "SEN", "IVC"), ordered = TRUE)


posterior_plot <- tmp %>% 
  ggplot(aes(x = name, y = posterior, fill = population)) +
  geom_bar(stat = "identity", width = 1.0) +
  scale_fill_manual(values=c(my_pal)) +
  labs(title = "Posterior Membership Probability (K=4)") +
  ylab("Posterior membership probability") +
  xlab("Samples") +
  facet_grid(~origin, scales = "free_x", space = "free" ) + 
  theme_bw() +
  guides(fill=guide_legend(title="Membership")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off")

#### combining the plots ####
ggarrange(
  ggarrange(
    k_plot,
    dapc_plot,
    ncol = 2, labels = c("A", "B")
  ),
  posterior_plot,
  nrow = 2,
  labels = c("", "C"),
  heights = c(2, 1.5)
)

ggsave("DAPC_neutral.png", height = 6, width = 11, units = "in")

#### 3d Interactive plots #####
library(plotly)
colors <- c("black", "#E69F00", "#56B4E9", "dodgerblue", "grey50", "red")

locdf <- as.data.frame(by_location$ind.coord)
locdf$population <- as.factor(metadata$population)                       
locdf$year <- as.factor(metadata$year)

fig <- plot_ly(locdf, x = ~LD1, y = ~LD2, z = ~LD3, color = ~population, colors = colors, size = 15)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'LD1'),
                                   yaxis = list(title = 'LD2'),
                                   zaxis = list(title = 'LD3')))

fig
htmlwidgets::saveWidget(as_widget(fig), paste0("dapc_", i, ".html"))


fig2 <- plot_ly(locdf, x = ~LD1, y = ~LD2, z = ~LD3, color = ~year, colors = colors, size = 15)
fig2 <- fig2 %>% add_markers()
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = 'LD1'),
                                     yaxis = list(title = 'LD2'),
                                     zaxis = list(title = 'LD3')))
fig2

htmlwidgets::saveWidget(as_widget(fig2), paste0("dapc_year_", i, ".html"))