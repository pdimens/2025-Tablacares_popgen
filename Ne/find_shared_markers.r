## Figure out which SNPs are shared between the loci used to
## construct the linkage map, and those used in the popgen study

library(tidyverse)
setwd("~/Omega/USM PhD/Projects/Active/Yellowfin Tuna/Popgen/Analyses/nonkin/Ne")

# make the linkage map all nice
linkmap <- read_delim(
  "../../../../Linkage Map/LepAnchor/miss90_both/miss90_both.map",
  col_names = FALSE
) %>%
  select(1,2,5,6) %>% 
  rename("lg" = X1, "position" = X2, "male" = X5, "female" = X6) %>% 
  arrange(lg, position) %>% 
  rowwise() %>% mutate(cm = round(mean(c(male, female)), digits = 4)) %>% 
  mutate(male = round(male, digits = 4), female = round(female, digits = 4))

# save it for forever
#write.table(linkmap, file = "sexaverage.map", row.names = FALSE, quote = FALSE)

linkmap <- read_delim("sexaverage.map")
snp_pos <- read_delim("../data/scaffold_positions.snp")
# get the index of the snp positions for later filtering
snp_pos$idx <- 1:nrow(snp_pos)
# subset the positions to extract only snps on scaffolds
scaffold_pos <- snp_pos[grepl("LG", snp_pos$scaffold),]
names(scaffold_pos) <- c("lg", "position", "idx")

# nothing :()
match_pos <- merge(scaffold_pos, linkmap, by = c("lg", "position"))

approx_markers <- function(snps, lnkmap, thresh){
  lgs <- unique(snps$lg)
  # instantiate output vectors
  out.lg <- c()
  out.pos <- c()
  out.linkpos <- c()
  out.cm <- c()
  out.idx <- c()
  for(lg in lgs){
    # subset by linkage group
    submap <- lnkmap[lnkmap$lg == lg,]
    subsnp <- snps[snps$lg == lg,]
    # get length - 1 to prevent overflow of while loop
    l <- nrow(submap) - 1
    j <- 1  # get the snp index
    for(snp in subsnp$position){
      # setup while conditions
      i <- 0  # while loop iterator
      cond <- FALSE
      while (cond == FALSE){
        i <- i + 1
        # check if the snp position is between the pair of adjacent markers
        cond <- between(snp, submap$position[i],submap$position[i+1])
        # break if reaching the end
        if(i == l){break()}
      }
      if(cond == TRUE){
        # check to see if the marker falls within the threshold
        diff_lower <- abs(snp - submap$position[i])
        diff_upper <- abs(submap$position[i+1] - snp)
        # if one of them does, add the marker to the output
        if(diff_lower <= thresh){
          out.lg <- c(out.lg, lg)
          out.pos <- c(out.pos, snp)
          out.linkpos <- c(out.linkpos, submap$position[i])
          out.cm <- c(out.cm, submap$cm[i])
          out.idx <- c(out.idx, subsnp$idx[j])
        } else if(diff_upper <= thresh){
          out.lg <- c(out.lg, lg)
          out.pos <- c(out.pos, snp)
          out.linkpos <- c(out.linkpos, submap$position[i+1])
          out.cm <- c(out.cm, submap$cm[i+1])
          out.idx <- c(out.idx, subsnp$idx[j])
        }
      }
    j <- j+1
    }
  }
  data.frame(lg = out.lg, snp.pos = out.pos, link.pos = out.linkpos, cm = out.cm, idx = out.idx)
}

matchy <- approx_markers(scaffold_pos, linkmap, 5000)



recomb_rate <- linkmap %>%  
  mutate(position = position / 1e6) %>%
  group_by(lg) %>%
  mutate(recomb = max(cm) / (max(position) - min(position))) %>% 
  select(lg, recomb) %>% unique

scaff_conversion <- merge(scaffold_pos, recomb_rate, by = "lg")
