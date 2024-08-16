## Looking at windowed coverage across the chimeric scaffold, W_unloc_scaf_3 ##
## Faye Romero. University of Rochester. ##
## 15 July 2024. ##

# Set wd, load packages
setwd('/Users/Faye/Documents/FSJ_genome/Revisions/W_unloc_scaf_3/')
library(tidyverse)
library(ggview)

# Load in data
load("./AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.Rdata")

#### Just FYI: Generate windowed coverage data frame (AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.Rdata) ####
# Define a function to read each file and add the SampleID column
read_and_label <- function(file_path) {
  sample_id <- str_extract(basename(file_path), "H\\d+")
  data <- read_tsv(file_path)
  data <- data %>%
    mutate(SampleID = sample_id)
  data$Window <- 1:nrow(data)
  return(data)
}

# Get a list of all relevant files in the directory
files <- list.files(pattern = "H\\d+_AvgCov_Wunloc3_AllPos_1kbwin\\.txt")

# Read and combine all files into one data frame
combined_data <- map_df(files, read_and_label)
final <- as.data.frame(combined_data) # convert into data frame
rm(combined_data)

# Merge in sex data
sex <- read.table("/Users/Faye/Documents/FSJ_genome/Revisions/weirdZ/AllHighCov_08March2024_Sex.tsv", sep = "\t", stringsAsFactors = F)
colnames(sex) <- c("SampleID", "Sex")

# Merge
Wcov <- merge(final, sex, by = "SampleID")
write.table(Wcov, '/Users/Faye/Documents/FSJ_genome/Revisions/W_unloc_scaf_3/AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.csv', sep = ",", quote = F, row.names = F, col.names = T)
save(Wcov, file = "AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.Rdata")

#### Plot ####
p1 <- ggplot(data = Wcov, aes(x = Window, y = AverageCoverage, color = Sex)) +
  geom_line(alpha = 0.8) +
  labs(x = "Position (1 Kb windows)", y = "Average Coverage") +
  scale_color_manual(values = c("#F2B25F", "#639CF6")) +
  theme_classic() +
  annotate("text", x = 0, y = 38, label = "N = 50", size = 4, color = "black") +
  theme(text = element_text(family = "sans"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

# View and plot
ggview(plot = p1, units = "in", width = 10, height = 4)
ggsave("/Users/Faye/Documents/FSJ_genome/Revisions_v3/Wunloc3_coverage.pdf", p1, units = "in", width = 10, height = 4)
