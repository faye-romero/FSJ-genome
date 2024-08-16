## Plotting proportion of each chromosome made of repetitive elements ##
## Faye Romero. University of Rochester. 01 Aug 2024 ##

setwd('/Users/Faye/Documents/FSJ_genome/Revisions_v3/')
library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(ggview)

# Load in RepeatMasker .out
repeats_raw <- fread("/Users/Faye/Documents/FSJ_genome/Revisions_v3/data/FSJgenome_July2024_FINAL.fasta.out", sep=" ", fill = TRUE, skip=3)
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

# Function that computes proportion of each chromosome made up of each repeat family
compute_proportion <- function(query, repeats_df, lengths_df) {
  total_length <- lengths_df %>% filter(Query == query) %>% pull(Total_Length)
  
  repeat_families <- repeats_df %>%
    filter(Query == query) %>%
    mutate(Overlap_Length = Pos_in_Query_End - Pos_in_Query_Start + 1) %>%
    group_by(Repeat_class_family) %>%
    summarise(Total_Base_Pairs = sum(Overlap_Length)) %>%
    mutate(Proportion = Total_Base_Pairs / total_length,
           Query = query) %>%
    select(Query, Repeat_class_family, Proportion)
  
  return(repeat_families)
}

# Apply the function to each chromosome, then combine into one data frame
all_results <- data.frame()
queries <- unique(repeats_df$Query)

for (query in queries) {
  query_results <- compute_proportion(query, repeats_df, lengths_df)
  all_results <- bind_rows(all_results, query_results)
}

# Retain only the Repeat_class_family elements we want
all_results_filt <- filter(all_results, Repeat_class_family %in% c("Unspecified", "LTR", "Non-LTR", "LINE", "DNA transposon", "housekeeping RNAs", "Simple repeat"))

# Optional: Use convert file (a 2-column file with chromosome name, shortened chromosome name) to change chromosome names to readable
convert <- read.table("//Users/Faye/Documents/FSJ_genome/Revisions_v3/data/FSJgenome_July2024_FINAL.convert.txt", header = FALSE,  sep = "\t")
colnames(convert) <- c("NewNames", "Query")
final <- merge(all_results_filt, convert, by = "Query")
final$Query <- NULL

# Order chromosomes via factor
final$NewNames = factor(final$NewNames, levels=c("1","1A","2", "3", "4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28","29", "30", "31", "34", "W", "Z"))

# Plot 
ggplot(data = final, aes(x = NewNames, y = Proportion, fill = Repeat_class_family)) +
  geom_col(position = "stack") +
  scale_fill_manual("", values = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "grey")) +
  theme_classic() +
  labs(x = "Chromosome", y = "Proportion of Chromosome", fill = "") +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14))

# View and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 12, height = 5)
ggsave("/Users/Faye/Documents/FSJ_genome/Revisions_v3/Repeats_Prop.pdf", ggplot2::last_plot(), units = "in", width = 12, height = 5)
