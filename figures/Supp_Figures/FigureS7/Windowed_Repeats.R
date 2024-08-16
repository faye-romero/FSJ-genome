### Plotting windowed repetitive content per chromosome ##
## Faye Romero. University of Rochester. 02 Aug 2024. ##

# Set wd and load packages
setwd('/Users/Faye/Documents/FSJ_genome/Revisions/')
library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(ggview)

# Load in RepeatMasker .out
repeats_raw <- fread("/Users/Faye/Documents/FSJ_genome/Revisions_v3/data/FSJgenome_July2024_FINAL.fasta.out",sep=" ",fill = TRUE,skip=3)
colnames(repeats_raw) <- c("SW_score", "Perc_div", "Perc_del", "Perc_ins", "Query", "Pos_in_Query_Start", "Pos_in_Query_End", "Pos_in_Query_left", "Strand", "Matching_repeat", "Repeat_class_family", "Pos_in_Repeat_Start", "Pos_in_Repeat_End", "Pos_in_Repeat_left", "ID", "Star")

# Retain relevant columns
repeats_df <- select(repeats_raw, Query, Pos_in_Query_Start, Pos_in_Query_End, Matching_repeat, Repeat_class_family, ID)

# Group names
repeats_df <- repeats_df[repeats_df$Repeat_class_family != "ARTEFACT",]
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("Unknown","Unspecified","Other")] <- "Unspecified"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("LTR?","LTR", "LTR/ERV1", "LTR/ERV3", "LTR/ERVK", "LTR/ERVL", "LTR/ERVL?", "LTR/Unknown", "LTR/ERV")] <- "LTR"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("SINE","SINE2","SINE?", "Retroposon", "SINE2/tRNA", "SINE/MIR", "SINE/Deu")] <- "Non-LTR"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("LINE","LINE?", "LINE/CR1", "LINE/CR1?", "LINE/L2", "LINE/R2")] <- "LINE"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("DNA","Tc1-Mariner","RC","RC?","DNA?")] <- "DNA transposon"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("tRNA","rRNA","snRNA","scRNA")] <- "housekeeping RNAs"
repeats_df$Repeat_class_family[repeats_df$Repeat_class_family %in% c("Low_complexity", "Simple_repeat")] <- "Simple repeat"

# Load in chromosome lengths (a tab-separated file with 2 columns: scaffold name, length)
lengths_df <- read.table("/Users/Faye/Documents/FSJ_genome/Revisions_v3/data/scaffold.lengths_FSJgenome_July2024_FINAL.txt", sep = "\t", stringsAsFactors = F)
colnames(lengths_df) <- c("Query", "Total_Length")

# Function that calculates repeat content in windows for one chromosome, in 500 Kb windows. Output will be in MB.
compute_total_base_pairs <- function(query, repeats_df, lengths_df, window_size = 500000) {
  total_length <- lengths_df %>% filter(Query == query) %>% pull(Total_Length)
  windows <- seq(1, total_length, by = window_size)
  
  result <- data.frame()
  
  for (start in windows) {
    end <- start + window_size - 1
    window_repeats <- repeats_df %>%
      filter(Query == query, 
             Pos_in_Query_Start <= end, 
             Pos_in_Query_End >= start)
    
    if (nrow(window_repeats) > 0) {
      repeat_families <- window_repeats %>%
        mutate(Overlap_Start = pmax(Pos_in_Query_Start, start),
               Overlap_End = pmin(Pos_in_Query_End, end),
               Overlap_Length = Overlap_End - Overlap_Start + 1) %>%
        group_by(Repeat_class_family) %>%
        summarise(Total_Base_Pairs = sum(Overlap_Length)) %>%
        mutate(Total_Base_Pairs_in_MB = Total_Base_Pairs / 1000000,
               Query = query,
               Window_Start_MB = start / 1000000,
               Window_End_MB = end / 1000000) %>%
        select(Query, Window_Start_MB, Window_End_MB, Repeat_class_family, Total_Base_Pairs_in_MB)
      
      result <- bind_rows(result, repeat_families)
    }
  }
  
  return(result)
}

# Apply the function to each chromosome, then combine into one data frame
all_results <- data.frame()
queries <- unique(repeats_df$Query)

for (query in queries) {
  query_results <- compute_total_base_pairs(query, repeats_df, lengths_df)
  all_results <- bind_rows(all_results, query_results)
}

# Optional: Use convert file (a 2-column file with chromosome name, shortened chromosome name) to change chromosome names to readable
convert <- read.table("//Users/Faye/Documents/FSJ_genome/Revisions_v3/data/FSJgenome_July2024_FINAL.convert.txt", header = FALSE,  sep = "\t")
colnames(convert) <- c("NewNames", "Query")
final <- merge(all_results, convert, by = "Query")
final$Query <- NULL

# Order chromosomes via factor
final$NewNames = factor(final$NewNames, levels=c("1","1A","2", "3", "4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28","29", "30", "31", "34", "W", "Z"))

# Plot
filter(final, Repeat_class_family %in% c("LTR", "Non-LTR", "LINE", "DNA transposon", "housekeeping RNAs", "Simple repeat", "Unspecified")) %>%
  ggplot(., aes(x = Window_Start_MB, y = Total_Base_Pairs_in_MB, fill = Repeat_class_family)) +
  geom_area(position = "stack") +
  facet_wrap(~NewNames, strip.position = "bottom", scales = "free", ncol = 5) +
  labs(x = "Position (Mb)", y = "Mbp") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9),
        legend.text = element_text(size = 8.5),
        strip.background = element_rect(fill = "white"),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "white", fill = NA, linewidth = 0)) +
        scale_fill_manual("", values = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "grey"))

# View plot and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 9, height = 9)
ggsave("/Users/Faye/Documents/FSJ_genome/Revisions_v3/Windowed_Repeats.pdf", ggplot2::last_plot(), units = "in", width = 9, height = 9)
