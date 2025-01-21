library(tidyverse) |> suppressPackageStartupMessages()
library(qqman)
setwd("/mnt/Win10/Users/Pavel/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/FST")
## pairwise locus-by-locus Hudson fst ##
fst <- read.csv("locbyloc.pairwise.hudson.fst")
fst$pop1 <- factor(fst$pop1, levels = c("GOM", "VZ", "IVC", "SEN"), ordered = TRUE)
fst$pop2 <- factor(fst$pop2, levels = c("wATL", "GOM", "IVC", "SEN"), ordered = TRUE)
fst$locus <- as.factor(fst$locus)

## snp info ##
positions <- read.table("snp_pos.bim", header = F)
names(positions) <- c("chrom", "locus", "ignore", "bp", "allele1", "allele2")
positions <- select(positions, 1,2,4)
positions$locus <- tolower(positions$locus)

## global loc-by-loc fst ##
globalfst <- read.csv("locbyloc.global.nei.fst")


## merge fst with marker info ##
# pairwise
tb <- merge(fst, positions, by = "locus")
tb <- filter(tb, chrom != 0) %>% select(1,5,6,2,3,4)
outloci <- c("snp_1571","snp_1962","snp_2565","snp_4665")
outliers <- tb[tb$locus]

# global
gtb <- merge(globalfst, positions, by = "locus")
gtb <- filter(gtb, chrom != 0) %>% select(1,12,13, 8,9)
names(gtb) <- c("locus", "chrom", "bp", "fst", "fstp")

#### ggplot global ####
## get cumulative snp position rather than absolute
don <- gtb %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gtb, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, bp) %>%
  mutate( bp_cum=bp+tot)

# setup the X axis
axisdf = don %>% group_by(chrom) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )

ggplot(don, aes(x=bp_cum, y=fstp)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "grey50"), length(unique(don$chrom)) )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Linkage Group") +
  ylab(expression("F"[ST]))+
  labs(title = expression(paste("Global F"[ST], " of ddRAD markers across the genome"))) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("manhattan.globalFST.png", height = 8, width = 12, units = "in")
#### pairwise FST ####

pon <- gtb %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(tb, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, bp) %>%
  mutate( bp_cum=bp+tot, popval=paste(pop1,pop2, sep = "-"))
  pon$popval <- as.factor(pon$popval)


axispdf = pon %>% group_by(chrom) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )
pon$pop1 <- factor(pon$pop1, levels = c("GOM", "IVC", "SEN", "VZ"), ordered = TRUE)
pon$pop2 <- factor(pon$pop2, levels = c("wATL", "GOM", "IVC", "SEN"), ordered = TRUE)


ggplot(pon, aes(x=bp_cum, y=fst)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chrom)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "grey50"), length(unique(don$chrom)) )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Linkage Group") +
  ylab(expression("F"[ST]))+
  labs(title = expression(paste("Pairwise F"[ST], " of ddRAD markers across the genome"))) +
  # Custom the theme:
  facet_grid(vars(pop1), vars(pop2)) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("manhattan.pairwiseFST.png", height = 8, width = 19, units = "in")
