## TE lengths

setwd('/Users/Faye/Documents/FSJ_genome/Revisions_v3/')
library(tidyverse)
library(ggview)

# Load TE data
load("TE_lengths_v123_July2024.Rdata")
load("Plot_TE_lengths_v123_July2024.Rdata")
load("TE_counts_v23_July2024.Rdata")

#### TE median lengths ####

# Remove v1, as we're not considering it for this project
plot2 <- filter(plot, genome != "v1")
final2 <- filter(final, genome != "v1")

# Perform Wilcoxon tests (clunky; there's a better way!)
final2 %>%
  filter(repeat_class.family == "DNA transposon") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

final2 %>%
  filter(repeat_class.family == "housekeeping RNAs") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value = 0.0001246 ***

final2 %>%
  filter(repeat_class.family == "LINE") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

final2 %>%
  filter(repeat_class.family == "LTR") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

final2 %>%
  filter(repeat_class.family == "Non-LTR") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value = 0.3262 NS

final2 %>%
  filter(repeat_class.family == "Satellite") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

final2 %>%
  filter(repeat_class.family == "Simple repeat") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

final2 %>%
  filter(repeat_class.family == "Unspecified") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16 ***

# Create a boxplot given the above significances
order <- c("DNA transposon", "housekeeping RNAs", "LINE", "LTR", "Non-LTR", "Satellite", "Simple repeat", "Unspecified")
plot$repeat_class.family <- factor(plot$repeat_class.family, levels = order)
ggplot(plot2, aes(x = repeat_class.family, y = TE_length, fill = genome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"), labels = c("Short-read", "Long-read")) +
  labs(x = "TE superfamily", y = "Median TE Length", fill = "Genome") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  coord_flip() +
  ylim(0,1300) +
  geom_segment(aes(x = 0.6, y = 170, xend = 1.4, yend = 170)) +
  annotate("text", x = 1, y = 220, label = "***", size = 6, color = "black") + #DNA transposon
  geom_segment(aes(x = 1.6, y = 120, xend = 2.4, yend = 120)) +
  annotate("text", x = 2, y = 170, label = "***", size = 6, color = "black") + #housekeeping RNAs
  geom_segment(aes(x = 2.6, y = 260, xend = 3.4, yend = 260)) +
  annotate("text", x = 3, y = 310, label = "***", size = 6, color = "black") + #LINE
  geom_segment(aes(x = 3.6, y = 570, xend = 4.4, yend = 570)) +
  annotate("text", x = 4, y = 620, label = "***", size = 6, color = "black") + #LTR
  geom_segment(aes(x = 5.6, y = 1250, xend = 6.4, yend = 1250)) +
  annotate("text", x = 6, y = 1300, label = "***", size = 6, color = "black") + #Satellite
  geom_segment(aes(x = 6.6, y = 80, xend = 7.4, yend = 80)) +
  annotate("text", x = 7, y = 130, label = "***", size = 6, color = "black") + #Simple repeat
  geom_segment(aes(x = 7.6, y = 200, xend = 8.4, yend = 200)) +
  annotate("text", x = 8, y = 250, label = "***", size = 6, color = "black") + #Unspecified
  theme(text = element_text(family = "sans"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

# View and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 11, height = 5)
ggsave("Median_TElength.pdf", plot = last_plot(), units = "in", width = 11, height = 5)

#### TE counts ####
#Use a proportion test for each category to see whether counts differ significantly
test_data <- table(final2$genome[final2$repeat_class.family == "Unspecified"])
prop.test(test_data)
#All TE families returned p < 0.001, so all ***

#Plotting counts
ggplot(counts, aes(x = repeat_class.family, y = count, fill = genome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"), labels = c("Short-read", "Long-read")) +
  labs(x = "TE superfamily", y = "Count", fill = "Genome") +
  theme_classic() +
  geom_segment(aes(x = 0.6, y = 26000, xend = 1.4, yend = 26000)) +
  annotate("text", x = 1, y = 28000, label = "***", size = 6, color = "black") + #DNA transposon
  geom_segment(aes(x = 1.6, y = 10000, xend = 2.4, yend = 10000)) +
  annotate("text", x = 2, y = 12000, label = "***", size = 6, color = "black") + #housekeeping RNAs
  geom_segment(aes(x = 2.6, y = 210000, xend = 3.4, yend = 210000)) +
  annotate("text", x = 3, y = 211000, label = "***", size = 6, color = "black") + #LINE
  geom_segment(aes(x = 3.6, y = 190000, xend = 4.4, yend = 190000)) +
  annotate("text", x = 4, y = 193000, label = "***", size = 6, color = "black") + #LTR
  geom_segment(aes(x = 4.6, y = 16000, xend = 5.4, yend = 16000)) +
  annotate("text", x = 5, y = 17000, label = "***", size = 6, color = "black") + #Non-LTR
  geom_segment(aes(x = 5.6, y = 60000, xend = 6.4, yend = 60000)) +
  annotate("text", x = 6, y = 61000, label = "***", size = 6, color = "black") + #Satellite
  geom_segment(aes(x = 6.6, y = 290000, xend = 7.4, yend = 290000)) +
  annotate("text", x = 7, y = 291000, label = "***", size = 6, color = "black") + #Simple repeat
  geom_segment(aes(x = 7.6, y = 75000, xend = 8.4, yend = 75000)) +
  annotate("text", x = 8, y = 76000, label = "***", size = 6, color = "black") + #Unspecified
  theme(text = element_text(family = "sans"),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))

# View and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 11, height = 6)
ggsave("TE_counts.pdf", plot = last_plot(), units = "in", width = 11, height = 6)

#### Just FYI: calculate TE counts/lengths for FSJ genome V1 (Feng et al. 2020) ####
# Read in repeat tables (.out from RepeatMasker)
v1 <- read.csv("./final.assembly.fasta.out.csv", header = T, stringsAsFactors = F, sep = ",")

v1 <- v1[v1$repeat_class.family != "ARTEFACT",]
v1$repeat_class.family[v1$repeat_class.family %in% c("Unspecified", "Unknown")] <- "Unspecified"
v1$repeat_class.family[v1$repeat_class.family %in% c("LTR/Copia", "LTR/ERV", "LTR/ERV1", "LTR/ERV4", "LTR/ERVK", "LTR/ERVL", "LTR/Gypsy", "LTR/Unknown", "LTR/?", "LTR?", "LTR/ERVL?", "LTR/ERV?")] <- "LTR"
v1$repeat_class.family[v1$repeat_class.family %in% c("SINE/5S", "SINE/5S-Deu-L2", "SINE/MIR", "SINE/tRNA", "SINE/tRNA-CR1")] <- "Non-LTR"
v1$repeat_class.family[v1$repeat_class.family %in% c("LINE/CR1", "LINE/I-Jockey", "LINE/L1", "LINE/L2", "LINE/Penelope", "LINE/R2", "LINE/RTE-BovB", "LINE/CR1?", "LINE?", "LINEE/L1?", "LINE?/Penelope?")] <- "LINE"
v1$repeat_class.family[v1$repeat_class.family %in% c("DNA","DNA transposon", "DNA/CMC-EnSpm", "DNA/CMC-Transib", "DNA/Crypton", "DNA/hAT", "DNA/hAT-Charlie", "DNA/Merlin", "DNA/P", "DNA/PIF-Harbing", "DNA/TcMar-Tc1", "RC/Helitron", "DNA?", "DNA/DNA?", "DNA/DNA", "DNA/TcMar?", "DNA/hAT?", "RC/Helitron?", "RC?/Helitron?")] <- "DNA transposon"
v1$repeat_class.family[v1$repeat_class.family %in% c("house RNAs", "rRNA", "tRNA", "rRNA", "snRNA", "scRNA")] <- "housekeeping RNAs"
v1$repeat_class.family[v1$repeat_class.family %in% c("Simple repeat", "Simple_repeat", "Low_complexity")] <- "Simple repeat"
v1$repeat_class.family[v1$repeat_class.family %in% c("Satellite", "Satellite/acro","Satellite/macro")] <- "Satellite"
keep <- c("DNA transposon", "housekeeping RNAs", "LINE", "LTR", "Non-LTR", "Satellite", "Simple repeat", "Unspecified")
v1_filt <- subset(v1, repeat_class.family %in% keep)
v1_filt$pos_in_query_end <- as.numeric(v1_filt$pos_in_query_end)
v1_filt$pos_in_query_begin <- as.numeric(v1_filt$pos_in_query_begin)

v1_length <- mutate(v1_filt, TE_length = pos_in_query_end - pos_in_query_begin)
v1_length$genome <- "v1"

#V1 tallies
v1_counts <- v1_length %>%
  group_by(repeat_class.family) %>%
  summarize(count = n())

#Median lengths
v1_median <- aggregate(TE_length ~ repeat_class.family, data = v1_length, FUN = median)
v1_median$genome <- "v1"

# Create a boxplot with ggplot2
ggplot(v1_length, aes(x = repeat_class.family, y = TE_length)) +
  geom_boxplot(fill = "lightblue") +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", position = position_dodge(width = 0.75)) +
  labs(title = "Boxplot of Length by Family", x = "Family", y = "Length") +
  theme_minimal() +
  ylim(0,300)

#### Just FYI: calculate TE counts/lengths for FSJ genome V2 (Driscoll and Beaudry et al. 2024) ####

v2 <- read.csv("./FSJ.v2_1.fasta.out.July2024.csv", header = T, stringsAsFactors = F, sep = ",")

v2 <- v2[v2$repeat_class.family != "ARTEFACT",]
v2$repeat_class.family[v2$repeat_class.family %in% c("Unspecified", "Unknown")] <- "Unspecified"
v2$repeat_class.family[v2$repeat_class.family %in% c("LTR/Copia", "LTR/ERV", "LTR/ERV1", "LTR/ERV4", "LTR/ERVK", "LTR/ERVL", "LTR/Gypsy", "LTR/Unknown", "LTR/?", "LTR?", "LTR/ERVL?", "LTR/ERV?")] <- "LTR"
v2$repeat_class.family[v2$repeat_class.family %in% c("SINE/5S", "SINE/5S-Deu-L2", "SINE/MIR", "SINE/tRNA", "SINE/tRNA-CR1")] <- "Non-LTR"
v2$repeat_class.family[v2$repeat_class.family %in% c("LINE/CR1", "LINE/I-Jockey", "LINE/L1", "LINE/L2", "LINE/Penelope", "LINE/R2", "LINE/RTE-BovB", "LINE/CR1?", "LINE?", "LINEE/L1?", "LINE?/Penelope?")] <- "LINE"
v2$repeat_class.family[v2$repeat_class.family %in% c("DNA","DNA transposon", "DNA/CMC-EnSpm", "DNA/CMC-Transib", "DNA/Crypton", "DNA/hAT", "DNA/hAT-Charlie", "DNA/Merlin", "DNA/P", "DNA/PIF-Harbing", "DNA/TcMar-Tc1", "RC/Helitron", "DNA?", "DNA/DNA?", "DNA/DNA", "DNA/TcMar?", "DNA/hAT?", "RC/Helitron?", "RC?/Helitron?")] <- "DNA transposon"
v2$repeat_class.family[v2$repeat_class.family %in% c("house RNAs", "rRNA", "tRNA", "rRNA", "snRNA", "scRNA")] <- "housekeeping RNAs"
v2$repeat_class.family[v2$repeat_class.family %in% c("Simple repeat", "Simple_repeat", "Low_complexity")] <- "Simple repeat"
v2$repeat_class.family[v2$repeat_class.family %in% c("Satellite", "Satellite/acro","Satellite/macro")] <- "Satellite"
keep <- c("DNA transposon", "housekeeping RNAs", "LINE", "LTR", "Non-LTR", "Satellite", "Simple repeat", "Unspecified")
v2_filt <- subset(v2, repeat_class.family %in% keep)
v2_filt$pos_in_query_end <- as.numeric(v2_filt$pos_in_query_end)
v2_filt$pos_in_query_begin <- as.numeric(v2_filt$pos_in_query_begin)

v2_length <- mutate(v2_filt, TE_length = pos_in_query_end - pos_in_query_begin)
v2_length$genome <- "v2"

#v2 tallies
v2_counts <- v2_length %>%
  group_by(repeat_class.family) %>%
  summarize(count = n())
v2_counts$genome <- "v2"

test <- v2_length %>%
  filter(str_detect(query_seq, "ch")) %>%
  group_by(query_seq, repeat_class.family) %>%
  summarise(total_te_length = sum(TE_length))

#Median lengths
v2_median <- aggregate(TE_length ~ repeat_class.family, data = v2_length, FUN = median)
v2_median$genome <- "v2"

# Create a boxplot with ggplot2
ggplot(v2_length, aes(x = repeat_class.family, y = TE_length)) +
  geom_boxplot(fill = "lightblue") +
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "red", position = position_dodge(width = 0.75)) +
  labs(title = "Boxplot of Length by Family", x = "Family", y = "Length") +
  theme_minimal() +
  ylim(0,300)

#### Just FYI: calculate TE counts/lengths for FSJ genome V3 (this study) ####
v3 <- read.csv("./FSJgenome_July2024_FINAL.fasta.out.csv", header = T, stringsAsFactors = F, sep = ",")

v3 <- v3[v3$repeat_class.family != "ARTEFACT",]
v3$repeat_class.family[v3$repeat_class.family %in% c("Unspecified", "Unknown")] <- "Unspecified"
v3$repeat_class.family[v3$repeat_class.family %in% c("LTR/Copia", "LTR/ERV", "LTR/ERV1", "LTR/ERV4", "LTR/ERVK", "LTR/ERVL", "LTR/Gypsy", "LTR/Unknown", "LTR/?", "LTR?", "LTR/ERVL?", "LTR/ERV?")] <- "LTR"
v3$repeat_class.family[v3$repeat_class.family %in% c("SINE/5S", "SINE/5S-Deu-L2", "SINE/MIR", "SINE/tRNA", "SINE/tRNA-CR1")] <- "Non-LTR"
v3$repeat_class.family[v3$repeat_class.family %in% c("LINE/CR1", "LINE/I-Jockey", "LINE/L1", "LINE/L2", "LINE/Penelope", "LINE/R2", "LINE/RTE-BovB", "LINE/CR1?", "LINE?", "LINEE/L1?", "LINE?/Penelope?")] <- "LINE"
v3$repeat_class.family[v3$repeat_class.family %in% c("DNA","DNA transposon", "DNA/CMC-EnSpm", "DNA/CMC-Transib", "DNA/Crypton", "DNA/hAT", "DNA/hAT-Charlie", "DNA/Merlin", "DNA/P", "DNA/PIF-Harbing", "DNA/TcMar-Tc1", "RC/Helitron", "DNA?", "DNA/DNA?", "DNA/DNA", "DNA/TcMar?", "DNA/hAT?", "RC/Helitron?", "RC?/Helitron?")] <- "DNA transposon"
v3$repeat_class.family[v3$repeat_class.family %in% c("house RNAs", "rRNA", "tRNA", "rRNA", "snRNA", "scRNA")] <- "housekeeping RNAs"
v3$repeat_class.family[v3$repeat_class.family %in% c("Simple repeat", "Simple_repeat", "Low_complexity")] <- "Simple repeat"
v3$repeat_class.family[v3$repeat_class.family %in% c("Satellite", "Satellite/acro","Satellite/macro")] <- "Satellite"
keep <- c("DNA transposon", "housekeeping RNAs", "LINE", "LTR", "Non-LTR", "Satellite", "Simple repeat", "Unspecified")
v3_filt <- subset(v3, repeat_class.family %in% keep)
v3_filt$pos_in_query_end <- as.numeric(v3_filt$pos_in_query_end)
v3_filt$pos_in_query_begin <- as.numeric(v3_filt$pos_in_query_begin)

v3_length <- mutate(v3_filt, TE_length = pos_in_query_end - pos_in_query_begin)
v3_length_2 <- filter(v3_length, TE_length != -99994598)
v3_length_2$genome <- "v3"

#v3 tallies
v3_counts <- v3_length_2 %>%
  group_by(repeat_class.family) %>%
  summarize(count = n())
v3_counts$genome <- "v3"

#Median lengths
v3_median <- aggregate(TE_length ~ repeat_class.family, data = v3_length_2, FUN = median)
v3_median$genome <- "v3"

final <- rbind(v1_length, v2_length, v3_length_2)
#save(final, file = "TE_lengths_v123_July2024.Rdata")
plot <- rbind(v1_median, v2_median, v3_median)
#save(plot, file = "Plot_TE_lengths_v123_July2024.Rdata")
counts <- rbind(v2_counts, v3_counts)
#save(counts, file = "TE_counts_v23_July2024.Rdata")
