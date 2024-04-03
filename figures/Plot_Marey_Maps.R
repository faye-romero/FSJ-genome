## Create a Marey Map ##
## Faye Romero. 31 March 2024 ##

#### Set-up ####
setwd('/Users/Faye/Documents/FSJ_genome/') #Set wd

# Load libraries
library(tidyverse)
library(gridExtra)
library(scales)

load("LOD5anonmap.rdata") #Read in linkage map. Must have physical and genetic position of markers

#### Data cleaning ####
# Re-name chromosomes to just numbers for cleanliness
anonmap$NewScaff <- gsub("\\bChr2\\b", "2", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr3\\b", "3", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr1\\b", "1", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr29\\b", "29", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr4\\b", "4", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChrZ\\b", "Z", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr5\\b", "5", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr7\\b", "7", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr6\\b", "6", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr8\\b", "8", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr11\\b", "11", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr9\\b", "9", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr12\\b", "12", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr31\\b", "31", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr10\\b", "10", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr13\\b", "13", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr14\\b", "14", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr20\\b", "20", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr15\\b", "15", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr18\\b", "18", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr19\\b", "19", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr17\\b", "17", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr21\\b", "21", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr26\\b", "26", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr23\\b", "23", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr24\\b", "24", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr28\\b", "28", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr27\\b", "27", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr22\\b", "22", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr25\\b", "25", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr30\\b", "30", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr34\\b", "34", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr31\\b", "31", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr1A\\b", "1A", anonmap$NewScaff)
anonmap$NewScaff <- gsub("\\bChr4A\\b", "4A", anonmap$NewScaff)

#Remove small unlocalized scaffolds
plot <- filter(anonmap, !NewScaff %in% c("Chr2_unlocalized_s1", "Chr2_unlocalized_s2", "ChrZ_unlocalized_s3"))

#Order by factor
final <- plot %>% mutate(across(NewScaff, ~factor(., levels=c("1", "1A", "2", "3", "4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28","29", "30", "31", "33","34", "Z")))) %>%
  filter(!SNP.Name %in% c("s2722p51223","s3416p45084", "s54p1216454", "s8544p300", "s638p23447", "s581p328910", "s1256p156067", "s3531p28650", "s4014p34966", "s7669p254", "s373p97829", "s1351p15391")) #Remove outlier markers

# Make the data frame plottable by sex
femalemap <- select(final, SNP.Name, chr, SNP, cMPosition.Female, NewScaff, Position, LinkageGroup)
femalemap$Sex <- "Female"
colnames(femalemap)[4] <- "cMPosition"
malemap <- select(final, SNP.Name, chr, SNP, cMPosition.Male, NewScaff, Position, LinkageGroup)
malemap$Sex <- "Male"
colnames(malemap)[4] <- "cMPosition"
totalmap <- select(final, SNP.Name, chr, SNP, sexAvg, NewScaff, Position, LinkageGroup)
totalmap$Sex <- "Sex-averaged"
colnames(totalmap)[4] <- "cMPosition"

# Bind into one data frame
Marey <- rbind(femalemap, malemap, totalmap)

#### Plot ####
Marey %>%
  mutate(Position_scaled = Position / 1000000) %>%  # Create a new column with scaled values
  arrange(cMPosition) %>%
  ggplot(., aes(x = Position_scaled, y = cMPosition, color = Sex)) +
  scale_color_manual(values = c("#F2B25F", "#639CF6", "#101010")) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~ NewScaff, scales = 'free') +
  theme_classic() +
  labs(x = "Physical position (Mb)", y = "Genetic position (cM)") +
  theme(axis.title = element_text(size = 18), strip.text.x = element_text(size = 15),
        axis.text.x = element_text(size = 9)) +
  guides(color = guide_legend(title = NULL, label.theme = element_text(size = 12)))