## Create a Marey Map ##
## Faye Romero. 31 March 2024 ##

#### Set-up ####
setwd('/Users/Faye/Documents/FSJ_genome/Revisions_v3') #Set wd

# Load libraries
library(tidyverse)
library(gridExtra)
library(scales)
library(ggview)

map <- read.table("../LOD3LOD5_anon_map.csv", sep = ",", header = T, stringsAsFactors = F) #Read in linkage map. Must have physical and genetic position of markers

#### Data cleaning ####
#Remove the "Chr" prefix for cleanliness
map$Scaffold <- gsub("Chr", "", map$Scaffold)

# Order by factor
final <- map %>% mutate(across(Scaffold, ~factor(., levels=c("1", "1A", "2", "3", "4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28","29", "30", "31", "33","34", "Z")))) %>%
  filter(!SNP %in% c("10620","5852", "6953", "5122", "7796", "8133", "10414", "365", "10269", "10945", "9043")) #Remove outlier markers

# Pivot longer for plotting
final_long <- final %>%
  pivot_longer(cols = c(SexAveraged_cM, Female_cM, Male_cM),
               names_to = "Sex",
               values_to = "cM.Position") %>%
  mutate(Sex = gsub("_cM", "", Sex))

#### Plot ####
final_long %>%
  mutate(SNPpos_scaled = Position / 1000000) %>%  # Create a new column with scaled values
  arrange(cM.Position) %>%
  ggplot(., aes(x = SNPpos_scaled, y = cM.Position, color = Sex)) +
  scale_color_manual(values = c("#F2B25F", "#639CF6", "#101010")) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~ Scaffold, scales = 'free') +
  theme_classic() +
  labs(x = "Physical position (Mb)", y = "Genetic position (cM)") +
  theme(axis.title = element_text(size = 16), strip.text.x = element_text(size = 13),
        axis.text.x = element_text(size = 9)) +
  guides(color = guide_legend(title = NULL, label.theme = element_text(size = 12)))

# View and save
ggview(plot = ggplot2::last_plot(), units = "in", width = 11.5, height = 9)
ggsave("./Figure2.pdf", plot = last_plot(), units = "in", width = 11.5, height = 9)
