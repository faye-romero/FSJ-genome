## Plot sex-linked scaffolds ##
## Faye Romero. 31 March 2024 ##

#### Set-up ####
# Load libraries
library(dplyr)
library(tidyr)
library(tidyverse)

# Create combined data frame
# Set the directory path. This directory has files describing the average coverage of each scaffold for each individual.
directory_path <- "/Users/Faye/Documents/FSJ_genome/avg_cov_out_two"

# Make list of all files in the directory
file_list <- list.files(directory_path, full.names = TRUE)

# Initialize empty list to store the data tables
data_tables <- list()

# Loop through each file and read it into a data table
for (file_path in file_list) {
  # Extract the sample ID from the file name
  sample_id <- gsub("_.*", "", basename(file_path))
  
  # Read the file into a data table (assuming it's a space-separated text file)
  data <- read.table(file_path, header = FALSE, col.names = c("Scaffold", "Average_Coverage"))
  
  # Add the SampleID column
  data$SampleID <- sample_id
  
  # Store the data table in the list
  data_tables[[sample_id]] <- data
}

# Combine all data frames into one large data frame
combined_data <- do.call(rbind, data_tables)

#### Incoporate sex data ####
sex <- read.table("/Users/Faye/Documents/FSJ_genome/random_samples_13Feb2024.txt", header = TRUE, sep = "\t")
temp <- merge(combined_data, sex, by = "SampleID")
merged_data <- filter(temp, !SampleID == "H1422") #Remove random female, so that we have an even 25 males and 25 females

# Calculate genome-wide average coverage for each individual
merged_data_2 <- merged_data %>%
  group_by(SampleID) %>%
  reframe(SampleID = SampleID, Scaffold = Scaffold, Average_Coverage = Average_Coverage, Genome_Mean_Average_Coverage = mean(Average_Coverage), Sex = Sex)

#### Statistical tests #####
# Bonferroni correction
# Use a t-test for each scaffold
final_data <- merged_data_2
scaffold_test_results <- by(final_data, final_data$Scaffold, function(subset) {
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

#sig_results[1,3]/sig_results[1,2]... I made 442 comparisons:
#0.01/442... old sig p-value of 0.01 = new sig p-value of 0.00002262443439 (2.26 * 10-5)
#0.001/442...old sig p-value of 0.001 = new sig p-value of 0.000002262443439 (2.26 * 10-6)
#0.0001/442... old sig p-value of 0.0001 = new sig p-value of 0.000000226244344 (2.26 * 10-7)

# Save results that have an adjusted p-value of < 0.05
sig_results <- filter(scaffold_test_results_df, adjusted_p.value <= 0.05)

#### More data wrangling ####
# Merge scaffold_test_results_df with final_data based on Scaffold
master_df <- merge(scaffold_test_results_df, final_data, by = "Scaffold")

# Load in key with new names for plotting
key <- read.csv("/Users/Faye/Documents/FSJ_genome/key_AphCoe_V3_internal_Mar2024-VERSUS-FSJv3_internal_Feb2024.csv", header = T, stringsAsFactors = F)
colnames(key) <- c("New_Scaffold", "Scaffold")
final <- merge(master_df, key, by = "Scaffold", all = T)

# Retain chosen chromosomes/scaffolds for plotting (those that are large enough to double-check with separate synteny analyses)
plot <- filter(final, New_Scaffold %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "ChrZ", "ChrZ_unlocalized_scaffold_1", "ChrZ_unlocalized_scaffold_2", "ChrW_unlocalized_scaffold_1", "ChrW_unlocalized_scaffold_2", "ChrW_unlocalized_scaffold_3", "ChrW_unlocalized_scaffold_4"))

# Specify order
order <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "ChrZ", "ChrZ_unlocalized_scaffold_1", "ChrZ_unlocalized_scaffold_2", "ChrW_unlocalized_scaffold_1", "ChrW_unlocalized_scaffold_2", "ChrW_unlocalized_scaffold_3", "ChrW_unlocalized_scaffold_4")
plot$New_Scaffold <- factor(plot$New_Scaffold, levels = order)

#### Plot ####
# Create box plot...
ggplot(plot, aes(x = New_Scaffold, y = Average_Coverage, fill = Sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F2B25F", "#639CF6")) +
  labs(x = "Scaffold", y = "Average Coverage", fill = "Sex") +
  theme_bw() +
  ylim(0,32) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  annotate("text", x = 11, y = 28, label = "***", size = 5, color = "black") +
  annotate("text", x = 12, y = 32, label = "***", size = 5, color = "black") +
  annotate("text", x = 13, y = 25, label = "***", size = 5, color = "black") +
  annotate("text", x = 14, y = 15, label = "***", size = 5, color = "black") +
  annotate("text", x = 15, y = 15, label = "***", size = 5, color = "black") +
  annotate("text", x = 16, y = 15, label = "***", size = 5, color = "black") +
  annotate("text", x = 17, y = 15, label = "***", size = 5, color = "black")
