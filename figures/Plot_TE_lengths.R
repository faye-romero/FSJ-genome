## Calculate TE counts and median TE lengths from RepeatMasker outputs ##
## Faye Romero. 31 March 2024 ##

#### Set-up ####
# Set wd and load libraries
setwd('/Users/Faye/Documents/FSJ_genome/FIGURES/')
library(tidyverse)

# Load in TE data
load("TE_lengths.Rdata")
load("TE_lengths_plottable.Rdata")
load("TE_counts.Rdata")

# Remove v1, because I'm not using for these analyses
plot2 <- filter(plot, genome != "v1")
final2 <- filter(final, genome != "v1")

#### Median TE lengths ####
# Create a barplot
order <- c("DNA transposon", "house RNAs", "LINE", "LTR", "Non-LTR", "Satellite", "Simple repeat", "Unspecified")
plot$repeat_class.family <- factor(plot$repeat_class.family, levels = order)
ggplot(plot2, aes(x = repeat_class.family, y = TE_length, fill = genome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A")) +
  labs(x = "TE superfamily", y = "Median TE Length", fill = "Genome") +
  ylim(0,845) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  geom_segment(aes(x = 0.6, y = 150, xend = 1.4, yend = 150)) +
  annotate("text", x = 1, y = 170, label = "***", size = 6, color = "black") + #DNA transposon
  geom_segment(aes(x = 1.6, y = 535, xend = 2.4, yend = 535)) +
  annotate("text", x = 2, y = 555, label = "***", size = 6, color = "black") + #house RNAs
  geom_segment(aes(x = 2.6, y = 230, xend = 3.4, yend = 230)) +
  annotate("text", x = 3, y = 250, label = "***", size = 6, color = "black") + #LINE
  geom_segment(aes(x = 3.6, y = 385, xend = 4.4, yend = 385)) +
  annotate("text", x = 4, y = 405, label = "***", size = 6, color = "black") + #LTR
  geom_segment(aes(x = 4.6, y = 120, xend = 5.4, yend = 120)) +
  annotate("text", x = 5, y = 135, label = "**", size = 6, color = "black") + #Non-LTR
  geom_segment(aes(x = 5.6, y = 820, xend = 6.4, yend = 820)) +
  annotate("text", x = 6, y = 840, label = "***", size = 6, color = "black") + #Satellite
  geom_segment(aes(x = 6.6, y = 55, xend = 7.4, yend = 55)) +
  annotate("text", x = 7, y = 75, label = "***", size = 6, color = "black") + #Simple repeat
  geom_segment(aes(x = 7.6, y = 295, xend = 8.4, yend = 295)) +
  annotate("text", x = 8, y = 315, label = "***", size = 6, color = "black")#Unspecified

# Perform Wilcoxon test on each of the TE superfamilies (there's probably a more efficient way to do this)
final2 %>%
  filter(repeat_class.family == "DNA transposon") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "house RNAs") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "LINE") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "LTR") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value = 5.969e-07

final2 %>%
  filter(repeat_class.family == "Non-LTR") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value = 0.007715 **

final2 %>%
  filter(repeat_class.family == "Satellite") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "Simple repeat") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "Unspecified") %>%
  wilcox.test(TE_length ~ genome, data = .)
#p-value < 2.2e-16

final2 %>%
  filter(repeat_class.family == "house RNAs") %>%
  group_by(genome) %>%
  summarise(median = median(TE_length))

#### TE counts ####
#Plotting counts
ggplot(counts, aes(x = repeat_class.family, y = count, fill = genome)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A")) +
  labs(x = "TE superfamily", y = "Count", fill = "Genome") +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15)) +
  geom_segment(aes(x = 0.6, y = 20000, xend = 1.4, yend = 20000)) +
  annotate("text", x = 1, y = 21000, label = "***", size = 6, color = "black") + #DNA transposon
  geom_segment(aes(x = 1.6, y = 8000, xend = 2.4, yend = 8000)) +
  annotate("text", x = 2, y = 9000, label = "***", size = 6, color = "black") + #house RNAs
  geom_segment(aes(x = 2.6, y = 210000, xend = 3.4, yend = 210000)) +
  annotate("text", x = 3, y = 211000, label = "***", size = 6, color = "black") + #LINE
  geom_segment(aes(x = 3.6, y = 280000, xend = 4.4, yend = 280000)) +
  annotate("text", x = 4, y = 281000, label = "***", size = 6, color = "black") + #LTR
  geom_segment(aes(x = 4.6, y = 13000, xend = 5.4, yend = 13000)) +
  annotate("text", x = 5, y = 14000, label = "***", size = 6, color = "black") + #Non-LTR
  geom_segment(aes(x = 5.6, y = 90000, xend = 6.4, yend = 90000)) +
  annotate("text", x = 6, y = 91000, label = "***", size = 6, color = "black") + #Satellite
  geom_segment(aes(x = 6.6, y = 295000, xend = 7.4, yend = 295000)) +
  annotate("text", x = 7, y = 296000, label = "***", size = 6, color = "black") + #Simple repeat
  geom_segment(aes(x = 7.6, y = 125000, xend = 8.4, yend = 125000)) +
  annotate("text", x = 8, y = 126000, label = "***", size = 6, color = "black")#Unspecified

# Conduct a proportion test for each TE superfamily. E.g., for DNA transposon:
test_data <- table(final2$genome[final2$repeat_class.family == "DNA transposon"])
prop.test(test_data)

#DNA transposon: p-value < 2.2e-16
#house RNAs: p-value = 1.379e-05
#LINE: p-value < 2.2e-16
#LTR: p-value < 2.2e-16
#Non-LTR: p-value < 2.2e-16
#Satellite: p-value < 2.2e-16
#Simple repeat: p-value < 2.2e-16
#Unspecified: p-value < 2.2e-16