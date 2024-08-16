## Plotting average read depth per chromosome across 25 male and 25 female Florida Scrub-Jays ##
## Faye Romero. University of Rochester. 01 Aug 2024 ##

# Load packages
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggstatsplot)
library(ggview)

# Load in data frame with average read depth per scaffold per sample
load("/Users/Faye/Documents/FSJ_genome/Revisions_v3/Indiv_Cov_Per_Scaf.Rdata")

## Statistical tests ##
## Bonferroni correction ##
# Use a t-test for each scaffold
scaffold_test_results <- by(Indiv.Cov.Per.Scaf, Indiv.Cov.Per.Scaf$Scaffold, function(subset) {
  # Check if there are enough observations in both groups
  if (length(unique(subset$Sex)) == 2 && all(table(subset$Sex) >= 2)) {
    # Check for variability in Average_Coverage
    if (length(unique(subset$Average_Coverage)) > 1) {
      t_test_result <- t.test(Average_Coverage ~ Sex, data = subset)
      return(data.frame(Scaffold = unique(subset$Scaffold), p.value = t_test_result$p.value))
    } else {
      # If no variability in Average_Coverage, return NA
      return(data.frame(Scaffold = unique(subset$Scaffold), p.value = NA))
    }
  } else {
    # If not enough observations, return NA
    return(data.frame(Scaffold = unique(subset$Scaffold), p.value = NA))
  }
})

# Convert into a single data frame
scaffold_test_results_df <- do.call(rbind, scaffold_test_results)

# Apply Bonferroni correction
scaffold_test_results_df$adjusted_p.value <- p.adjust(scaffold_test_results_df$p.value, method = "bonferroni")

# Save results that have an adjusted p-value of < 0.05 into a new data frame, apply significance codes
sig_results <- filter(scaffold_test_results_df, adjusted_p.value <= 0.05)
sig_results$sig[sig_results$adjusted_p.value <= 0.01] <- "*"
sig_results$sig[sig_results$adjusted_p.value <= 0.001] <- "**"
sig_results$sig[sig_results$adjusted_p.value <= 0.0001] <- "***"

# How many comparisons did you make?
sig_results[1,3]/sig_results[1,2] #329

# Merge scaffold_test_results_df with final_data based on Scaffold
master_df <- merge(scaffold_test_results_df, Indiv.Cov.Per.Scaf, by = "Scaffold")

# Load in key with shortened scaffold names for plotting purposes... this is just a 2-column file with scaffold name, shortened scaffold name
key <- read.csv("/Users/Faye/Documents/FSJ_genome/key_FSJgenome_July2024_FINAL_shortnames.csv", header = T, stringsAsFactors = F)
final <- merge(master_df, key, by = "Scaffold", all = T)

# Retain chosen chromosomes/scaffolds for plotting
plot <- filter(final, Plot_Scaffold_Name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "Z", "Z_unloc_scaf_1", "Z_unloc_scaf_2", "Z_unloc_scaf_3", "W_unloc_scaf_1", "W_unloc_scaf_2", "W_unloc_scaf_3", "W_unloc_scaf_4", "W_unloc_scaf_5"))

# Specify order
order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "Z", "Z_unloc_scaf_1", "Z_unloc_scaf_2", "Z_unloc_scaf_3", "W_unloc_scaf_1", "W_unloc_scaf_2", "W_unloc_scaf_3", "W_unloc_scaf_4", "W_unloc_scaf_5")
plot$Plot_Scaffold_Name <- factor(plot$Plot_Scaffold_Name, levels = order)

# Create box plot
ggplot(plot, aes(x = Plot_Scaffold_Name, y = Average_Coverage, fill = Sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F2B25F", "#639CF6")) +
  labs(x = "Scaffold", y = "Average Read Depth", fill = "Sex") +
  theme_bw() +
  ylim(0,38) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  #Add significance codes based on t-test sig_results table
  annotate("text", x = 11, y = 28, label = "***", size = 6, color = "black") + 
  annotate("text", x = 12, y = 32, label = "***", size = 6, color = "black") +
  annotate("text", x = 13, y = 23, label = "***", size = 6, color = "black") +
  annotate("text", x = 14, y = 38, label = "***", size = 6, color = "black") +
  annotate("text", x = 15, y = 15, label = "***", size = 6, color = "black") +
  annotate("text", x = 16, y = 15, label = "**", size = 6, color = "black") +
  annotate("text", x = 17, y = 15, label = "***", size = 6, color = "black") +
  annotate("text", x = 18, y = 15, label = "***", size = 6, color = "black") +
  annotate("text", x = 19, y = 15, label = "***", size = 6, color = "black")

# View and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 12, height = 6)
ggsave("/Users/Faye/Documents/FSJ_genome/Revisions_v3/FINAL-FILES-v3/Figure4.pdf", plot = last_plot(), units = "in", width = 12, height = 6)
