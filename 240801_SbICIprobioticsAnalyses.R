# AUTHOR INFORMATION ###########################################################
#                                                                              #
# Statistical Analyses for Gorushi et al., 'A yeast-based oral therapeutic     #  
# platform for delivery of immune checkpoint inhibitors reduces intestinal     #
# tumor burden. (2024)                                                         #
#                                                                              #
# Bioconductor version 3.14 (BiocManager 1.30.18), R 4.1.3 (2022-03-10)        #
#                                                                              #
# Contact: Jerome Prusa (prusa@wustl.edu)                                      #
#                                                                              #   
#                                                                              #
# LOAD PACKAGES ################################################################

install.packages("BiocManager")
library(dada2)
library(phyloseq)
library(reshape2)
library(dplyr)
library(ggplot2)
library(NatParksPalettes)
library(pals)
library(vegan)
library(labdsv)
library(ggpubr)
library(rstatix)
library(Maaslin2)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(labdsv)
library(Hmisc)

# LOAD DATA ####################################################################


#  metadata_d36_d50_2           =  Metadata for the 49 and 48 mice that stool was 
#                               collected from days 36 and 50, respectively.
#  genus_abund_d36_d50_2        =  Genus relative abundance table, unfiltered.
#  fam_abund_d36_d50_2          =  Family relative abundance table, unfiltered
#  phyl_abund_d36_d50_2         =  Phylum relative abundance table, unfiltered

metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
fam_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/fam_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
phyl_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/phyl_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)

# Filtering genera at a relative abundance cutoff of 0.1% and ordering them by relative abundance

genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]


# Generating a plot showing the top 15 genera for day 50 samples
top15genera <- colnames(genus_abund_d36_d50)[1:15]

genus_abund_d36_d50_top15 <- data.frame(genus_abund_d36_d50[,colnames(genus_abund_d36_d50) %in% top15genera],Others=rowSums(genus_abund_d36_d50[,!colnames(genus_abund_d36_d50) %in% top15genera]))
genus_abund_d36_d50_top15 <- melt(as.matrix(genus_abund_d36_d50_top15))
colnames(genus_abund_d36_d50_top15)[1] <- "Sample"
colnames(genus_abund_d36_d50_top15)[2] <- "Genus"
colnames(genus_abund_d36_d50_top15)[3] <- "Value"
genus_abund_d36_d50_top15 <- merge(genus_abund_d36_d50_top15, metadata_d36_d50, by = "Sample")
genus_abund_d36_d50_top15 <- genus_abund_d36_d50_top15[,c(1,2,3,4,6,7,8,9,10)]
genus_abund_d36_d50_top15 <- genus_abund_d36_d50_top15 %>% group_by(Genus, Day, Treatment, Sex, ID, Cage) %>% summarise(Value = mean(Value))

genus_abund_d36_d50_top15 <- genus_abund_d36_d50_top15 %>%
  filter(!(Day %in% c("36")))
genus_abund_d36_d50_top15$Day <- as.character(genus_abund_d36_d50_top15$Day)
genus_abund_d36_d50_top15$Day <- factor(genus_abund_d36_d50_top15$Day, levels = c("50"))

genus_abund_d36_d50_top15$Treatment <- as.character(genus_abund_d36_d50_top15$Treatment)
genus_abund_d36_d50_top15$Treatment <- factor(genus_abund_d36_d50_top15$Treatment, levels = c("PBS","Sb_haPD1","Sb_OG539","aPDL-1"))

NatParksPalettes <- list(
  Acadia = list(c("#212E52", "#444E7E", "#8087AA", "#B7ABBC", "#F9ECE8", "#FCC893", "#FEB424", "#FD8700", "#D8511D"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=TRUE),
  Arches = list(c("#1A3D82", "#0C62AF", "#4499F5", "#8FCAFD", "#F2F2F2", "#F0AC7D", "#CD622E", "#B14311", "#832B0F"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=TRUE),
  Arches2 = list(c("#3A1F46", "#7F4B89", "#B46DB3", "#E3A5D6", "#F3DAE4"), c(1, 2, 3, 4, 5), colorblind=TRUE),
  Banff = list(c("#006475", "#00A1B7", "#55CFD8", "#586028", "#898928", "#616571", "#9DA7BF"), c(2, 5, 1, 6, 3, 7, 4), colorblind=FALSE),
  BryceCanyon = list(c("#882314", "#C0532B", "#CF932C", "#674D53", "#8C86A0", "#724438", "#D5AB85"), c(1, 5, 2, 7, 4, 3, 6), colorblind=FALSE),
  CapitolReef = list(c("#291919", "#532A34", "#7C5467", "#878195", "#AEB2B7", "#D4D9DD"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
  Chamonix = list(c("#008FF8", "#B6AA0D", "#E2C2A2", "#E23B0E", "#F2C621", "#196689"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  CraterLake = list(c("#1D4A79", "#794C23", "#6B7444", "#6089B5", "#BF9785", "#275E4D", "#807B7F"), c(1, 2, 3, 4, 5, 6, 7), colorblind=FALSE),
  Cuyahoga = list(c("#E07529", "#FAAE32", "#7F7991", "#A84A00", "#5D4F36", "#B39085"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
  DeathValley = list(c("#8C2B0E", "#C5692D", "#FEB359", "#132F5B", "#435F90", "#68434E", "#B47E83"), c(1, 5, 7, 2, 6, 3, 4), colorblind=TRUE),
  Denali = list(c("#20223E", "#3F3F7B", "#278192", "#00B089", "#2EEA8C", "#8FF7BD"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE), 
  Everglades = list(c("#345023", "#596C0B", "#83A102", "#003B68", "#426F86", "#7A712F"), c(3, 4, 1, 6, 5, 2), colorblind=FALSE),
  Glacier = list(c("#01353D", "#088096", "#58B3C7", "#7AD4E4", "#B8FCFC"), c(1, 2, 3, 4, 5), colorblind=TRUE),
  GrandCanyon = list(c("#521E0F", "#9C593E", "#DDA569", "#3F4330", "#8E7E3C", "#2A4866", "#6592B0"), c(2, 6, 3, 4, 7, 1, 5), colorblind=FALSE),
  Halekala = list(c("#722710", "#A3844D", "#675243", "#A85017", "#838BAA"), c(1, 2, 3, 4, 5), colorblind=TRUE),
  IguazuFalls = list(c("#415521", "#97AD3D", "#4C3425", "#7F6552", "#5A8093", "#9FBAD3"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  KingsCanyon = list(c("#613921", "#A77652", "#F2C27B", "#AAC9ED", "#44637D", "#8E949F"), c(1, 5, 6, 3, 2, 4), colorblind=TRUE),
  LakeNakuru = list(c("#D76E9A", "#A1ACC8", "#AD3C36", "#332627", "#EACACF", "#AA6B77"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  Olympic = list(c("#3A4330", "#426737", "#75871B", "#BAB97D", "#FAF3CE", "#FDE16A", "#F9B40E", "#E88C23", "#A25933"), c(1, 2, 3, 4, 5, 6, 7, 8, 9), colorblind=FALSE),
  Redwood = list(c("#5E3B49", "#9B5F6B", "#BA817D", "#325731", "#6A9741", "#5F4E2F"), c(2, 5, 6, 3, 4, 1), colorblind=FALSE),
  RockyMtn = list(c("#274C31", "#A3AEB5", "#2F4B6A", "#8F8081", "#3F7156", "#6F89A7", "#5B5443"), c(1, 2, 3, 4, 5, 6, 7), colorblind=FALSE),
  Saguaro = list(c("#127088", "#C85729", "#92874B", "#CD8A39", "#AC3414", "#57643C"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  SmokyMtns = list(c("#42511A", "#889D35", "#D3D175", "#B50200", "#DA6C41", "#7C6E66", "#BCAFA6"), c(1, 4, 2, 6, 3, 5, 7), colorblind=FALSE),
  SouthDowns = list(c("#948D2A", "#D5B44D", "#89A4BF", "#F1D6B6", "#9B8358", "#577291"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  Torres = list(c("#2F397A", "#7391BD", "#894846", "#E9988C", "#535260", "#B7A7A6", "#785838", "#C68D61", "#4F6008", "#93995C"), c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), colorblind=FALSE),
  Triglav = list(c("#386EC2", "#B5B5B2", "#990006", "#625D0A", "#B9741F", "#213958"), c(1, 2, 3, 4, 5, 6), colorblind=TRUE),
  WindCave = list(c("#2F100E", "#6C3322", "#B07159", "#C9A197", "#E0CDCD"), c(1, 2, 3, 4, 5), colorblind=TRUE),
  Volcanoes = list(c("#082544", "#1E547D", "#79668C", "#DE3C37", "#F2DC7E"), c(1, 2, 3, 4, 5), colorblind=TRUE),
  Yellowstone = list(c("#0067A2", "#DFCB91", "#CB7223", "#289A84", "#7FA4C2", "#AF7E56"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE),
  Yosemite = list(c("#293633", "#3D5051", "#6B7F7F", "#87A1C7", "#516B95", "#304F7D"), c(1, 2, 3, 4, 5, 6), colorblind=FALSE))

Genus_order <- c("Lachnospiraceae.NK4A136.group", "Lactobacillus", "Clostridium.sensu.stricto.1", "Faecalibaculum", "Akkermansia", "Bacteroides", "Turicibacter", 
                     "Roseburia", "Lachnoclostridium", "Muribaculum", "X.Eubacterium..xylanophilum.group", "Oscillibacter", "Lachnospiraceae.FCS020.group", "Alistipes", "Marvinbryantia", "Others")
  
genus_abund_d36_d50_top15$Genus <- factor(genus_abund_d36_d50_top15$Genus, levels = Genus_order)
  
p <- ggplot(data = genus_abund_d36_d50_top15, aes(x = Day, y = Value, fill = Genus)) + geom_bar(stat = "identity", position = "fill") + facet_wrap(~ Treatment, ncol = 4) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = natparks.pals("DeathValley", 16)) + 
  ylab("Relative abundance")
p


# Filtering family at a relative abundance cutoff of 0.1% and ordering them by relative abundance

fam_abund_d36_d50 <- fam_abund_d36_d50[rowSums(fam_abund_d36_d50) > 0,]
fam_abund_d36_d50 <- fam_abund_d36_d50/rowSums(fam_abund_d36_d50)
fam_abund_d36_d50[fam_abund_d36_d50 < 0.001] <- 0
fam_abund_d36_d50 <- fam_abund_d36_d50/rowSums(fam_abund_d36_d50)
fam_abund_d36_d50 <- fam_abund_d36_d50[,order(colSums(fam_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]


# Generating a plot showing the top 10 family for day 50 samples
top10family <- colnames(fam_abund_d36_d50)[1:10]

fam_abund_d36_d50_top10 <- data.frame(fam_abund_d36_d50[,colnames(fam_abund_d36_d50) %in% top10family],Others=rowSums(fam_abund_d36_d50[,!colnames(fam_abund_d36_d50) %in% top10family]))
fam_abund_d36_d50_top10 <- melt(as.matrix(fam_abund_d36_d50_top10))
colnames(fam_abund_d36_d50_top10)[1] <- "Sample"
colnames(fam_abund_d36_d50_top10)[2] <- "Family"
colnames(fam_abund_d36_d50_top10)[3] <- "Value"
fam_abund_d36_d50_top10 <- merge(fam_abund_d36_d50_top10, metadata_d36_d50, by = "Sample")
fam_abund_d36_d50_top10 <- fam_abund_d36_d50_top10[,c(1,2,3,4,6,7,8,9,10)]
fam_abund_d36_d50_top10 <- fam_abund_d36_d50_top10 %>% group_by(Family, Day, Treatment, Sex, ID, Cage) %>% summarise(Value = mean(Value))

fam_abund_d36_d50_top10 <- fam_abund_d36_d50_top10 %>%
  filter(!(Day %in% c("36")))
fam_abund_d36_d50_top10$Day <- as.character(fam_abund_d36_d50_top10$Day)
fam_abund_d36_d50_top10$Day <- factor(fam_abund_d36_d50_top10$Day, levels = c("50"))

fam_abund_d36_d50_top10$Treatment <- as.character(fam_abund_d36_d50_top10$Treatment)
fam_abund_d36_d50_top10$Treatment <- factor(fam_abund_d36_d50_top10$Treatment, levels = c("PBS","Sb_haPD1","Sb_OG539","aPDL-1"))

Family_order <- c("Muribaculaceae", "Lachnospiraceae", "Erysipelotrichaceae", "Lactobacillaceae", "Clostridiaceae", "Oscillospiraceae", "Akkermansiaceae", 
                 "Bacteroidaceae", "Ruminococcaceae", "Rikenellaceae", "Others")

fam_abund_d36_d50_top10$Family <- factor(fam_abund_d36_d50_top10$Family, levels = Family_order)

p <- ggplot(data = fam_abund_d36_d50_top10, aes(x = Day, y = Value, fill = Family)) + geom_bar(stat = "identity", position = "fill") + facet_wrap(~ Treatment, ncol = 4) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = natparks.pals("DeathValley", 11)) + 
  ylab("Relative abundance")
p



# Filtering phyla at a relative abundance cutoff of 0.1% and ordering them by relative abundance

phyl_abund_d36_d50 <- phyl_abund_d36_d50[rowSums(phyl_abund_d36_d50) > 0,]
phyl_abund_d36_d50 <- phyl_abund_d36_d50/rowSums(phyl_abund_d36_d50)
phyl_abund_d36_d50[phyl_abund_d36_d50 < 0.001] <- 0
phyl_abund_d36_d50 <- phyl_abund_d36_d50/rowSums(phyl_abund_d36_d50)
phyl_abund_d36_d50 <- phyl_abund_d36_d50[,order(colSums(phyl_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]


# Generating a plot showing the top 4 phyla for day 50 samples
top5phyla <- colnames(phyl_abund_d36_d50)[1:5]

phyl_abund_d36_d50_top5 <- data.frame(phyl_abund_d36_d50[, colnames(phyl_abund_d36_d50) %in% top5phyla])
phyl_abund_d36_d50_top5 <- melt(as.matrix(phyl_abund_d36_d50_top5))
colnames(phyl_abund_d36_d50_top5)[1] <- "Sample"
colnames(phyl_abund_d36_d50_top5)[2] <- "Phylum"
colnames(phyl_abund_d36_d50_top5)[3] <- "Value"
phyl_abund_d36_d50_top5 <- merge(phyl_abund_d36_d50_top5, metadata_d36_d50, by = "Sample")
phyl_abund_d36_d50_top5 <- phyl_abund_d36_d50_top5[,c(1,2,3,4,6,7,8,9,10)]
phyl_abund_d36_d50_top5 <- phyl_abund_d36_d50_top5 %>% group_by(Phylum, Day, Treatment, Sex, ID, Cage) %>% summarise(Value = mean(Value))

phyl_abund_d36_d50_top5 <- phyl_abund_d36_d50_top5 %>%
  filter(!(Day %in% c("36")))
phyl_abund_d36_d50_top5$Day <- as.character(phyl_abund_d36_d50_top5$Day)
phyl_abund_d36_d50_top5$Day <- factor(phyl_abund_d36_d50_top5$Day, levels = c("50"))

phyl_abund_d36_d50_top5$Treatment <- as.character(phyl_abund_d36_d50_top5$Treatment)
phyl_abund_d36_d50_top5$Treatment <- factor(phyl_abund_d36_d50_top5$Treatment, levels = c("PBS","Sb_haPD1","Sb_OG539","aPDL-1"))

p <- ggplot(data = phyl_abund_d36_d50_top5, aes(x = Day, y = Value, fill = Phylum)) + geom_bar(stat = "identity", position = "fill") + facet_wrap(~ Treatment, ncol = 4) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = natparks.pals("DeathValley", 6)) + 
  ylab("Relative abundance")
p


###Calculating alpha diversity metrics (richness and Shannon diversity index)
richness_d36_d50 <- as.data.frame(specnumber(genus_abund_d36_d50))
richness_d36_d50$Sample <- rownames(richness_d36_d50)
colnames(richness_d36_d50)[1] <-"Richness"

shannon_div_d36_d50 <- as.data.frame(diversity(genus_abund_d36_d50, index = "shannon"))
shannon_div_d36_d50$Sample <- row.names(shannon_div_d36_d50)
colnames(shannon_div_d36_d50)[1] <- "ShannonDiv"

alpha_div_d36_d50 <- merge(richness_d36_d50, shannon_div_d36_d50, by = "Sample")
alpha_div_d36_d50 <- merge(alpha_div_d36_d50, metadata_d36_d50[, c(1,2,3,4,5,6,7,8)], by = "Sample")
alpha_div_d36_d50$Day <- factor(alpha_div_d36_d50$Day, levels = c("36", "50"))
alpha_div_d36_d50$Treatment <- factor(alpha_div_d36_d50$Treatment, levels = c("PBS","Sb_haPD1","Sb_OG539","aPDL-1"))


###Calculating pairwise Bray-Curtis dissimilarity values between all samples
all_Genus_BrayCurtis <- vegdist(genus_abund_d36_d50, method = "bray")
BC_pairwise <- as.matrix(all_Genus_BrayCurtis)
BC_pairwise <- melt(BC_pairwise)
BC_pairwise <- BC_pairwise[!duplicated(apply(BC_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
BC_pairwise <- BC_pairwise[which(BC_pairwise$Var1 != BC_pairwise$Var2),]
BC_pairwise <- merge(BC_pairwise, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")
colnames(BC_pairwise)[4] <- "Day1"
colnames(BC_pairwise)[5] <- "Treatment1"
BC_pairwise <- merge(BC_pairwise, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var2", by.y = "Sample")
colnames(BC_pairwise)[6] <- "Day2"
colnames(BC_pairwise)[7] <- "Treatment2"
BC_pairwise$DayPairs <- ""
BC_pairwise$TreatmentPairs <- ""

for (i in seq(1,dim(BC_pairwise)[1])) {
  BC_pairwise$DayPairs[i] <- paste(sort(c(BC_pairwise$Day1[i], BC_pairwise$Day2[i]))[1], sort(c(BC_pairwise$Day1[i], BC_pairwise$Day2[i]))[2])
}
for (i in seq(1,dim(BC_pairwise)[1])) {
  BC_pairwise$TreatmentPairs[i] <- paste(sort(c(BC_pairwise$Treatment1[i], BC_pairwise$Treatment2[i]))[1], sort(c(BC_pairwise$Treatment1[i], BC_pairwise$Treatment2[i]))[2])
}

####Subsetting the Bray-Curtis dissimilarities to only include dissimilarity values between day 36 and day 50 samples within each treatment group####
DayPairsToKeep <- c("36 50")
TreatmentPairsToKeep <- c("PBS PBS","Sb_haPD1 Sb_haPD1","Sb_OG539 Sb_OG539","aPDL-1 aPDL-1")
BC_pairwise_subset  <- BC_pairwise[which(BC_pairwise$DayPairs %in% DayPairsToKeep & BC_pairwise$TreatmentPairs %in% TreatmentPairsToKeep),]
BC_pairwise_subset$DayPairs <- factor(BC_pairwise_subset$DayPairs, levels = c("36 50"))
colnames(BC_pairwise_subset)[3] <- "Bray-Curtis Dissimilarity"


###Calculating PCoA vectors with Bray-Curtis dissimilarities for all day 36 and day 50 samples
pcoa_genus_d36_d50_BC <- pco(all_Genus_BrayCurtis, k = 4)
pcoa_genus_d36_d50_BC_eigen <- as.data.frame(pcoa_genus_d36_d50_BC$eig*100/sum(pcoa_genus_d36_d50_BC$eig))
pcoa_genus_d36_d50_BC <- as.data.frame(pcoa_genus_d36_d50_BC$points)
pcoa_genus_d36_d50_BC <- rownames_to_column(pcoa_genus_d36_d50_BC, var = "Var1")
metadata_d36_d50 <- rownames_to_column(metadata_d36_d50, var = "Sample")
pcoa_genus_d36_d50_BC <- merge(pcoa_genus_d36_d50_BC, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting day 36 and day 50 time points, and adding visualization aides including color and ellipses
ggplot() + 
  geom_point(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 36),], aes(x = V1, y = V2), color = "#949494", shape = 17) +
  geom_point(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 50),], aes(x = V1, y = V2, color = Treatment)) +
  stat_ellipse(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 50),], aes(x = V1, y = V2, color = Treatment), type = "norm", level = 0.95) +
  scale_color_manual(
    values = c("PBS" = "#000000", "Sb_OG539" = "#F78F26", "Sb_haPD1" = "#409CDC", "aPDL-1" = "#AC55EA"),
    breaks = c("PBS", "Sb_OG539", "Sb_haPD1", "aPDL-1")) +  
  xlab("PCoA 1 (40.6%)") + 
  ylab("PCoA 2 (14.1%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA))

###Merging the metadata with genus abundance data to subset for treatment groups for PERMANOVA comparisons between day 36 and day 50 and then removing unneeded columns from the merged data frame
genus_abund_d36_d50 <- rownames_to_column(genus_abund_d36_d50, var = "Sample")
genus_abund_d36_d50 <- merge(genus_abund_d36_d50, metadata_d36_d50, by = "Sample")
genus_abund_d36_d50 <- genus_abund_d36_d50[, -c(82:85, 87:88)]
rownames(genus_abund_d36_d50) <- genus_abund_d36_d50[, 1]
genus_abund_d36_d50 <- genus_abund_d36_d50[, -1]


###Calculating Bray-Curtis distances for d36 and d50 PBS samples and plotting PCoA with day 36 and day 50 samples differentiated by color, and ellipses drawn around all day 36 samples and another ellipse drawn around all day 50 samples###
PBS_d36_d50 <- c("PBS")
genus_abund_PBS_d36_d50 <- genus_abund_d36_d50[genus_abund_d36_d50$Treatment %in% PBS_d36_d50, ]
genus_abund_PBS_d36_d50_BrayCurtis <- vegdist(genus_abund_PBS_d36_d50[, c(-81)], method = "bray")
pcoa_PBS_d36_d50 <- pco(genus_abund_PBS_d36_d50_BrayCurtis, k = 4)
pcoa_PBS_d36_d50_eigen <- as.data.frame(pcoa_PBS_d36_d50$eig*100/sum(pcoa_PBS_d36_d50$eig))
pcoa_PBS_d36_d50 <- as.data.frame(pcoa_PBS_d36_d50$points)
pcoa_PBS_d36_d50 <- rownames_to_column(pcoa_PBS_d36_d50, var = "Var1")
pcoa_PBS_d36_d50 <- merge(pcoa_PBS_d36_d50, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting PCoA Bray-Curtis points for day 36 and day 50 PBS samples###
p <- ggplot() + 
  geom_point(data = pcoa_PBS_d36_d50[which(pcoa_PBS_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36")) +
  stat_ellipse(data = pcoa_PBS_d36_d50[which(pcoa_PBS_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36"), type = "norm", level = 0.95) +
  geom_point(data = pcoa_PBS_d36_d50[which(pcoa_PBS_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50")) +
  stat_ellipse(data = pcoa_PBS_d36_d50[which(pcoa_PBS_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50"), type = "norm", level = 0.95) +
  scale_color_manual(values = c("Day 36" = "#949494", "Day 50" = "#000000")) +
  xlab("PCoA 1 (38.1%)") + 
  ylab("PCoA 2 (26.7%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "PBS")
p



###Calculating Bray-Curtis distances for d36 and d50 Sb_OG539 samples and plotting PCoA with day 36 and day 50 samples differentiated by color, and ellipses drawn around all day 36 samples and another ellipse drawn around all day 50 samples###
OG539_d36_d50 <- c("Sb_OG539")
genus_abund_OG539_d36_d50 <- genus_abund_d36_d50[genus_abund_d36_d50$Treatment %in% OG539_d36_d50, ]
genus_abund_OG539_d36_d50_BrayCurtis <- vegdist(genus_abund_OG539_d36_d50[, c(-81)], method = "bray")
pcoa_OG539_d36_d50 <- pco(genus_abund_OG539_d36_d50_BrayCurtis, k = 4)
pcoa_OG539_d36_d50_eigen <- as.data.frame(pcoa_OG539_d36_d50$eig*100/sum(pcoa_OG539_d36_d50$eig))
pcoa_OG539_d36_d50 <- as.data.frame(pcoa_OG539_d36_d50$points)
pcoa_OG539_d36_d50 <- rownames_to_column(pcoa_OG539_d36_d50, var = "Var1")
pcoa_OG539_d36_d50 <- merge(pcoa_OG539_d36_d50, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting PCoA Bray-Curtis points for day 36 and day 50 Sb_OG539 samples###
p <- ggplot() + 
  geom_point(data = pcoa_OG539_d36_d50[which(pcoa_OG539_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36")) +
  stat_ellipse(data = pcoa_OG539_d36_d50[which(pcoa_OG539_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36"), type = "norm", level = 0.95) +
  geom_point(data = pcoa_OG539_d36_d50[which(pcoa_OG539_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50")) +
  stat_ellipse(data = pcoa_OG539_d36_d50[which(pcoa_OG539_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50"), type = "norm", level = 0.95) +
  scale_color_manual(values = c("Day 36" = "#949494", "Day 50" = "#F78F26")) +
  xlab("PCoA 1 (55.9%)") + 
  ylab("PCoA 2 (25.2%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "Sb_OG539")
p



###Calculating Bray-Curtis distances for d36 and d50 Sb_haPD1 samples and plotting PCoA with day 36 and day 50 samples differentiated by color, and ellipses drawn around all day 36 samples and another ellipse drawn around all day 50 samples###
haPD1_d36_d50 <- c("Sb_haPD1")
genus_abund_haPD1_d36_d50 <- genus_abund_d36_d50[genus_abund_d36_d50$Treatment %in% haPD1_d36_d50, ]
genus_abund_haPD1_d36_d50_BrayCurtis <- vegdist(genus_abund_haPD1_d36_d50[, c(-81)], method = "bray")
pcoa_haPD1_d36_d50 <- pco(genus_abund_haPD1_d36_d50_BrayCurtis, k = 4)
pcoa_haPD1_d36_d50_eigen <- as.data.frame(pcoa_haPD1_d36_d50$eig*100/sum(pcoa_haPD1_d36_d50$eig))
pcoa_haPD1_d36_d50 <- as.data.frame(pcoa_haPD1_d36_d50$points)
pcoa_haPD1_d36_d50 <- rownames_to_column(pcoa_haPD1_d36_d50, var = "Var1")
pcoa_haPD1_d36_d50 <- merge(pcoa_haPD1_d36_d50, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting PCoA Bray-Curtis points for day 36 and day 50 Sb_haPD1 samples###
p <- ggplot() + 
  geom_point(data = pcoa_haPD1_d36_d50[which(pcoa_haPD1_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36")) +
  stat_ellipse(data = pcoa_haPD1_d36_d50[which(pcoa_haPD1_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36"), type = "norm", level = 0.95) +
  geom_point(data = pcoa_haPD1_d36_d50[which(pcoa_haPD1_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50")) +
  stat_ellipse(data = pcoa_haPD1_d36_d50[which(pcoa_haPD1_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50"), type = "norm", level = 0.95) +
  scale_color_manual(values = c("Day 36" = "#949494", "Day 50" = "#409CDC")) +
  xlab("PCoA 1 (44.7%)") + 
  ylab("PCoA 2 (16.3%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "Sb_haPD1")
p



###Calculating Bray-Curtis distances for d36 and d50 aPDL-1 samples and plotting PCoA with day 36 and day 50 samples differentiated by color, and ellipses drawn around all day 36 samples and another ellipse drawn around all day 50 samples###
aPDL1_d36_d50 <- c("aPDL-1")
genus_abund_aPDL1_d36_d50 <- genus_abund_d36_d50[genus_abund_d36_d50$Treatment %in% aPDL1_d36_d50, ]
genus_abund_aPDL1_d36_d50_BrayCurtis <- vegdist(genus_abund_aPDL1_d36_d50[, c(-81)], method = "bray")
pcoa_aPDL1_d36_d50 <- pco(genus_abund_aPDL1_d36_d50_BrayCurtis, k = 4)
pcoa_aPDL1_d36_d50_eigen <- as.data.frame(pcoa_aPDL1_d36_d50$eig*100/sum(pcoa_aPDL1_d36_d50$eig))
pcoa_aPDL1_d36_d50 <- as.data.frame(pcoa_aPDL1_d36_d50$points)
pcoa_aPDL1_d36_d50 <- rownames_to_column(pcoa_aPDL1_d36_d50, var = "Var1")
pcoa_aPDL1_d36_d50 <- merge(pcoa_aPDL1_d36_d50, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting PCoA Bray-Curtis points for day 36 and day 50 aPDL-1 samples###
p <- ggplot() + 
  geom_point(data = pcoa_aPDL1_d36_d50[which(pcoa_aPDL1_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36")) +
  stat_ellipse(data = pcoa_aPDL1_d36_d50[which(pcoa_aPDL1_d36_d50$Day == 36),], aes(x = V1, y = V2, color = "Day 36"), type = "norm", level = 0.95) +
  geom_point(data = pcoa_aPDL1_d36_d50[which(pcoa_aPDL1_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50")) +
  stat_ellipse(data = pcoa_aPDL1_d36_d50[which(pcoa_aPDL1_d36_d50$Day == 50),], aes(x = V1, y = V2, color = "Day 50"), type = "norm", level = 0.95) +
  scale_color_manual(values = c("Day 36" = "#949494", "Day 50" = "#AC55EA")) +
  xlab("PCoA 1 (35.4%)") + 
  ylab("PCoA 2 (27.0%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "aPDL-1")
p



###PERMANOVA analysis comparing treatment groups at day 50###
#Comparing PBS and Sb_haPD-1 at d50#
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
PBSvSb_haPD1 <- c("PBS", "Sb_haPD1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% PBSvSb_haPD1, ]
PBSvSb_haPD1_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
                               permutations = 9999, method = "bray")


#Comparing Sb_haPD-1 and aPDL-1 at d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
aPDL1vSb_haPD1 <- c("aPDL-1", "Sb_haPD1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% aPDL1vSb_haPD1, ]
aPDL1vSb_haPD1_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
                               permutations = 9999, method = "bray")


#Comparing PBS and Sb_OG539 at d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
PBSvOG539 <- c("PBS", "Sb_OG539")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% PBSvOG539, ]
PBSvOG539_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
                               permutations = 9999, method = "bray")


#Comparing PBS and aPDL-1 at d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
PBSvPDL1 <- c("PBS", "aPDL-1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% PBSvPDL1, ]
PBSvPDL1_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
                            permutations = 9999, method = "bray")


#Comparing Sb_haPD1 and Sb_OG539 at d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
OG539vhaPD1 <- c("Sb_OG539", "Sb_haPD1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% OG539vhaPD1, ]
haPD1vOG539_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
                           permutations = 9999, method = "bray")



###PERMANOVA analysis comparing d36 and d50 for each treatment group###
#PBS d36 versus d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

TreatmentToKeep <- c("PBS")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% TreatmentToKeep, ]
PBS_36v50_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Day,
                               permutations = 9999, method = "bray")


#Sb_haPD1 d36 versus d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

TreatmentToKeep <- c("Sb_haPD1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% TreatmentToKeep, ]
Sb_haPD1_36v50_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Day,
                            permutations = 9999, method = "bray")


#Sb_OG539 d36 versus d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

TreatmentToKeep <- c("Sb_OG539")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% TreatmentToKeep, ]
Sb_OG539_36v50_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Day,
                                 permutations = 9999, method = "bray")


#aPDL1 d36 versus d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

TreatmentToKeep <- c("aPDL-1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% TreatmentToKeep, ]
aPDL1_36v50_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Day,
                                 permutations = 9999, method = "bray")


####Running Maaslin2 on the data####
#Re-inputting metadata and sample
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)

#Turning raw read count into relative abundance
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50[, -1]) > 0, ]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50, na.rm = TRUE)

#Turning .csv metadata to .txt
output_file <- "input_metadata"
write.table(metadata_d36_d50, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
metadata_d36_d50 <- read.table(output_file, header = TRUE)
rownames(metadata_d36_d50) <- metadata_d36_d50$Sample
metadata_d36_d50 <- metadata_d36_d50[, -1]

#Specifying day(s)/timepoint(s) for intput_data (ONLY for single experiments)
#timepoint <- c("36", "50") #Change this
input_metadata <- metadata_d36_d50

#choosing PBS as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="PBS"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0.1,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "240801_PBS_d36vd50_genus", 
                    fixed_effects  = c("Day"),
                    random_effects = c("Cage", "Sex"),
                    reference      = c("Day,36"))

#choosing Sb_haPD1 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="Sb_haPD1"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0.1,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "240801_haPD1_d36vd50_genus", 
                    fixed_effects  = c("Day"),
                    random_effects = c("Cage", "Sex"),
                    reference      = c("Day,36"))


#choosing Sb_OG539 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="Sb_OG539"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0.1,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "240801_OG539_d36vd50_genus", 
                    fixed_effects  = c("Day"),
                    random_effects = c("Cage", "Sex"),
                    reference      = c("Day,36"))
                    

#choosing aPDL-1 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="aPDL-1"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0.1,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "240801_aPDL1_d36vd50_genus", 
                    fixed_effects  = c("Day"),
                    random_effects = c("Cage", "Sex"),
                    reference      = c("Day,36"))
          


####Spearman correlation of just Sb_haPD1 d50 samples####
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
metadata_haPD1 <- subset(input_metadata, (Treatment =="Sb_haPD1"))
metadata_haPD1_d50 <- subset(metadata_haPD1, (Day =="50"))
metadata_haPD1_d50 <- rownames_to_column(metadata_haPD1_d50, var = "Sample")
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]
genus_abund_d36_d50 <- rownames_to_column(genus_abund_d36_d50, var = "Sample")
haPD1_d50_data <- merge(genus_abund_d36_d50, metadata_haPD1_d50, by = "Sample")
rownames(haPD1_d50_data) <- haPD1_d50_data$Sample
haPD1_d50_data <- haPD1_d50_data[, -1]
####Spearman with d50 haDP1 samples
burden_column <- haPD1_d50_data[, 84]
genus_columns <- haPD1_d50_data[, 1:80]
burden_genus_haPD1_d50 <- cbind(burden_column, genus_columns)
columns_to_keep <- apply(burden_genus_haPD1_d50, 2, function(col) !all(col == 0))
burden_genus_haPD1_d50_filtered <- burden_genus_haPD1_d50[, columns_to_keep]
genus_columns <- burden_genus_haPD1_d50_filtered[, 2:48]

cor_results_burden_genus_haPD1_d50_filtered <- rcorr(as.matrix(burden_genus_haPD1_d50_filtered), type = "spearman")
spearman_correlations_burden_genus_haPD1_d50_filtered <- cor_results_burden_genus_haPD1_d50_filtered$r[1, -1]
spearman_pvalues_burden_genus_haPD1_d50_filtered <- cor_results_burden_genus_haPD1_d50_filtered$P[1, -1]
spearman_burden_genus_haPD1_d50_filtered <- data.frame(
  Variable = colnames(genus_columns),
  Correlation = spearman_correlations_burden_genus_haPD1_d50_filtered,
  PValue = spearman_pvalues_burden_genus_haPD1_d50_filtered
)
adjusted_pvalues <- p.adjust(spearman_burden_genus_haPD1_d50_filtered$PValue, method = "BH")
spearman_burden_genus_haPD1_d50_filtered$AdjustedPValue <- adjusted_pvalues



####Spearman correlation of just Sb_OG539 d50 samples####
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
metadata_OG539 <- subset(input_metadata, (Treatment =="Sb_OG539"))
metadata_OG539_d50 <- subset(metadata_OG539, (Day =="50"))
metadata_OG539_d50 <- rownames_to_column(metadata_OG539_d50, var = "Sample")
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]
genus_abund_d36_d50 <- rownames_to_column(genus_abund_d36_d50, var = "Sample")
OG539_d50_data <- merge(genus_abund_d36_d50, metadata_OG539_d50, by = "Sample")
rownames(OG539_d50_data) <- OG539_d50_data$Sample
OG539_d50_data <- OG539_d50_data[, -1]
####Spearman with d50 OG539 samples
burden_column <- OG539_d50_data[, 84]
genus_columns <- OG539_d50_data[, 1:80]
burden_genus_OG539_d50 <- cbind(burden_column, genus_columns)
columns_to_keep <- apply(burden_genus_OG539_d50, 2, function(col) !all(col == 0))
burden_genus_OG539_d50_filtered <- burden_genus_OG539_d50[, columns_to_keep]
genus_columns <- burden_genus_OG539_d50_filtered[, 2:46]

cor_results_burden_genus_OG539_d50_filtered <- rcorr(as.matrix(burden_genus_OG539_d50_filtered), type = "spearman")
spearman_correlations_burden_genus_OG539_d50_filtered <- cor_results_burden_genus_OG539_d50_filtered$r[1, -1]
spearman_pvalues_burden_genus_OG539_d50_filtered <- cor_results_burden_genus_OG539_d50_filtered$P[1, -1]
spearman_burden_genus_OG539_d50_filtered <- data.frame(
  Variable = colnames(genus_columns),
  Correlation = spearman_correlations_burden_genus_OG539_d50_filtered,
  PValue = spearman_pvalues_burden_genus_OG539_d50_filtered
)
adjusted_pvalues <- p.adjust(spearman_burden_genus_OG539_d50_filtered$PValue, method = "BH")
spearman_burden_genus_OG539_d50_filtered$AdjustedPValue <- adjusted_pvalues



####Spearman correlation of just PBS d50 samples####
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
metadata_PBS <- subset(input_metadata, (Treatment =="PBS"))
metadata_PBS_d50 <- subset(metadata_PBS, (Day =="50"))
metadata_PBS_d50 <- rownames_to_column(metadata_PBS_d50, var = "Sample")
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]
genus_abund_d36_d50 <- rownames_to_column(genus_abund_d36_d50, var = "Sample")
PBS_d50_data <- merge(genus_abund_d36_d50, metadata_PBS_d50, by = "Sample")
rownames(PBS_d50_data) <- PBS_d50_data$Sample
PBS_d50_data <- PBS_d50_data[, -1]
####Spearman with d50 PBS samples
burden_column <- PBS_d50_data[, 84]
genus_columns <- PBS_d50_data[, 1:80]
burden_genus_PBS_d50 <- cbind(burden_column, genus_columns)
columns_to_keep <- apply(burden_genus_PBS_d50, 2, function(col) !all(col == 0))
burden_genus_PBS_d50_filtered <- burden_genus_PBS_d50[, columns_to_keep]
genus_columns <- burden_genus_PBS_d50_filtered[, 2:46]

cor_results_burden_genus_PBS_d50_filtered <- rcorr(as.matrix(burden_genus_PBS_d50_filtered), type = "spearman")
spearman_correlations_burden_genus_PBS_d50_filtered <- cor_results_burden_genus_PBS_d50_filtered$r[1, -1]
spearman_pvalues_burden_genus_PBS_d50_filtered <- cor_results_burden_genus_PBS_d50_filtered$P[1, -1]
spearman_burden_genus_PBS_d50_filtered <- data.frame(
  Variable = colnames(genus_columns),
  Correlation = spearman_correlations_burden_genus_PBS_d50_filtered,
  PValue = spearman_pvalues_burden_genus_PBS_d50_filtered
)
adjusted_pvalues <- p.adjust(spearman_burden_genus_PBS_d50_filtered$PValue, method = "BH")
spearman_burden_genus_PBS_d50_filtered$AdjustedPValue <- adjusted_pvalues



####Spearman correlation of just aPDL-1 d50 samples####
metadata_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/metadata_d36_d50_2.csv", header = TRUE)
metadata_aPDL1 <- subset(input_metadata, (Treatment =="aPDL-1"))
metadata_aPDL1_d50 <- subset(metadata_aPDL1, (Day =="50"))
metadata_aPDL1_d50 <- rownames_to_column(metadata_aPDL1_d50, var = "Sample")
genus_abund_d36_d50 <- read.csv("~/Desktop/Cancer Probiotics Manuscript 2/genus_abund_d36_d50_2.csv", row.names=1, check.names=FALSE)
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50) > 0,]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50[genus_abund_d36_d50 < 0.001] <- 0
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50)
genus_abund_d36_d50 <- genus_abund_d36_d50[,order(colSums(genus_abund_d36_d50, na.rm = TRUE),decreasing=TRUE)]
genus_abund_d36_d50 <- rownames_to_column(genus_abund_d36_d50, var = "Sample")
aPDL1_d50_data <- merge(genus_abund_d36_d50, metadata_aPDL1_d50, by = "Sample")
rownames(aPDL1_d50_data) <- aPDL1_d50_data$Sample
aPDL1_d50_data <- aPDL1_d50_data[, -1]
####Spearman with d50 aPDL-1 samples
burden_column <- aPDL1_d50_data[, 84]
genus_columns <- aPDL1_d50_data[, 1:80]
burden_genus_aPDL1_d50 <- cbind(burden_column, genus_columns)
columns_to_keep <- apply(burden_genus_aPDL1_d50, 2, function(col) !all(col == 0))
burden_genus_aPDL1_d50_filtered <- burden_genus_aPDL1_d50[, columns_to_keep]
genus_columns <- burden_genus_aPDL1_d50_filtered[, 2:48]

cor_results_burden_genus_aPDL1_d50_filtered <- rcorr(as.matrix(burden_genus_aPDL1_d50_filtered), type = "spearman")
spearman_correlations_burden_genus_aPDL1_d50_filtered <- cor_results_burden_genus_aPDL1_d50_filtered$r[1, -1]
spearman_pvalues_burden_genus_aPDL1_d50_filtered <- cor_results_burden_genus_aPDL1_d50_filtered$P[1, -1]
spearman_burden_genus_aPDL1_d50_filtered <- data.frame(
  Variable = colnames(genus_columns),
  Correlation = spearman_correlations_burden_genus_aPDL1_d50_filtered,
  PValue = spearman_pvalues_burden_genus_aPDL1_d50_filtered
)
adjusted_pvalues <- p.adjust(spearman_burden_genus_aPDL1_d50_filtered$PValue, method = "BH")
spearman_burden_genus_aPDL1_d50_filtered$AdjustedPValue <- adjusted_pvalues


