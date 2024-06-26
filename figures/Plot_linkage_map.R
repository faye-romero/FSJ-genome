## Plotting a linkage map with LinkageMapView ##
## Faye Romero. Feb 20, 2024 ##

#### Set-up ####
setwd('/Users/Faye/Documents/FSJ_genome/')

# Load packages
library(tidyverse)
#install.packages("LinkageMapView")
library(LinkageMapView)

# Read in linkage map
load("./LOD5anonmap.rdata")

#### cleaning ####
# Re-arrange to correct format: linakge group, cM position, locus
map.simple <- select(anonmap, LinkageGroup, sexAvg, SNP)
colnames(map.simple) <- c("Chromosome", "cMPosition", "Marker") # Re-name

# I just want marker ticks (no locus/marker labels). Make a dummy color scale so that I can plot a density map, but have no color.
sectcoldf <- lmvdencolor(map.simple,wsize = 50, colorin =
                           colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(25))
sectcoldf$col <- c("#FFFFFF")

#### Plot and save ####
lmv.linkage.plot(map.simple, "./FIGURES/linkage_map_v3.pdf",
                 denmap = TRUE,
                 cex.axis = 2,
                 cex.lgtitle = 2.5,
                 #lgperrow = 17,
                 sectcoldf = sectcoldf,
                 mapthese = c("1", "1A", "2","3", "4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31","34", "Z"))
