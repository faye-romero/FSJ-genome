## Plot a pie chart of RepeatMasker .tbl summary ##
## Faye Romero. 31 March 2024 ##

# Set wd and load libraries
setwd('/Users/Faye/Documents/FSJ_genome/FIGURES/')
library(tidyverse)
library(RColorBrewer)

# Load in data... this is just the .tbl output from RepeatMasker, but transformed into a .csv
data <- read.table("./repeat_summary_AphCoe_V3_internal_Mar2024.csv", header = T, stringsAsFactors = F, sep = ",")

# Manually make a data frame with desired info for ease of plotting
percents <- c(0.1, 4.17, 10.52, 0.1, 3.82, 81.29)
TE_fams <- c("DNA transposons", "LINE", "LTR", "Non-LTR", "Unspecified", "Non-repetitive content")
colors <- c("#ffd700", "#fa8675", "#e95f94", "#cd34b5", "#bebebe", "#878787")

# Plot
pie(percents, col = colors)