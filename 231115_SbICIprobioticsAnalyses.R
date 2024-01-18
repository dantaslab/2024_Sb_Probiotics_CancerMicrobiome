# AUTHOR INFORMATION ###########################################################
#                                                                              #
# Statistical Analyses for Gorushi et al., 'A yeast-based oral therapeutic     #  
# platform for delivery of immune checkpoint inhibitors reduces intestinal     #
# tumor burden. (2023)                                                         #
#                                                                              #
# Bioconductor version 3.14 (BiocManager 1.30.18), R 4.1.3 (2022-03-10)        #
#                                                                              #
# Contact: Jerome Prusa (prusa@wustl.edu)                                      #
#                                                                              #   
#                                                                              #
# LOAD PACKAGES ################################################################

library(dada2)
library(phyloseq)
library(reshape2)
library(dplyr)
library(ggplot2)
library(NatParksPalettes)
library(vegan)
library(labdsv)
library(ggpubr)
library(rstatix)
library(Maaslin2)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)

# LOAD DATA ####################################################################


#  metadata_d36_d50           =  Metadata for the 49 and 48 mice that stool was 
#                           collected from days 36 and 50, respectively.
#  genus_abund_d36_d50        =  Genus relative abundance table, unfiltered.
#  fam_abund_d36_d50          =  Family relative abundance table, unfiltered
#  phyl_abund_d36_d50         =  Phylum relative abundance table, unfiltered

metadata_d36_d50 <- read.csv("~/Desktop/metadata_d36_d50.csv", header = TRUE)
genus_abund_d36_d50 <- read.csv("~/Desktop/genus_abund_d36_d50.csv", row.names=1, check.names=FALSE)
fam_abund_d36_d50 <- read.csv("~/Desktop/fam_abund_d36_d50.csv", row.names=1, check.names=FALSE)
phyl_abund_d36_d50 <- read.csv("~/Desktop/phyl_abund_d36_d50.csv", row.names=1, check.names=FALSE)

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
top4phyla <- colnames(phyl_abund_d36_d50)[1:4]

phyl_abund_d36_d50_top4 <- data.frame(phyl_abund_d36_d50[,colnames(phyl_abund_d36_d50) %in% top4phyla],Others=rowSums(phyl_abund_d36_d50[,!colnames(phyl_abund_d36_d50) %in% top4phyla]))
phyl_abund_d36_d50_top4 <- melt(as.matrix(phyl_abund_d36_d50_top4))
colnames(phyl_abund_d36_d50_top4)[1] <- "Sample"
colnames(phyl_abund_d36_d50_top4)[2] <- "Phylum"
colnames(phyl_abund_d36_d50_top4)[3] <- "Value"
phyl_abund_d36_d50_top4 <- merge(phyl_abund_d36_d50_top4, metadata_d36_d50, by = "Sample")
phyl_abund_d36_d50_top4 <- phyl_abund_d36_d50_top4[,c(1,2,3,4,6,7,8,9,10)]
phyl_abund_d36_d50_top4 <- phyl_abund_d36_d50_top4 %>% group_by(Phylum, Day, Treatment, Sex, ID, Cage) %>% summarise(Value = mean(Value))

phyl_abund_d36_d50_top4 <- phyl_abund_d36_d50_top4 %>%
  filter(!(Day %in% c("36")))
phyl_abund_d36_d50_top4$Day <- as.character(phyl_abund_d36_d50_top4$Day)
phyl_abund_d36_d50_top4$Day <- factor(phyl_abund_d36_d50_top4$Day, levels = c("50"))

phyl_abund_d36_d50_top4$Treatment <- as.character(phyl_abund_d36_d50_top4$Treatment)
phyl_abund_d36_d50_top4$Treatment <- factor(phyl_abund_d36_d50_top4$Treatment, levels = c("PBS","Sb_haPD1","Sb_OG539","aPDL-1"))

p <- ggplot(data = phyl_abund_d36_d50_top4, aes(x = Day, y = Value, fill = Phylum)) + geom_bar(stat = "identity", position = "fill") + facet_wrap(~ Treatment, ncol = 4) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = natparks.pals("DeathValley", 5)) + 
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

###Plotting Shannon Diversity###
p <- ggplot(data = alpha_div_d36_d50, aes(x = Treatment, y = ShannonDiv, color = Day)) + 
  geom_violin(position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  theme_test() + ylab("Shannon Diversity (H)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=natparks.pals("Cuyahoga", 2), name = "Day") +
  labs(color = "Day")
p

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

####Plotting Bray-Curtis dissimilarities between day 36 and day 50 samples, subset by treatment group####
p <- ggplot(data = BC_pairwise_subset, aes(x = DayPairs, y = `Bray-Curtis Dissimilarity`, color = Treatment1)) + geom_boxplot() + scale_color_manual(values=natparks.pals("DeathValley", 4))
p

###Calculating PCoA vectors with Bray-Curtis dissimilarities for all day 36 and day 50 samples
pcoa_genus_d36_d50_BC <- pco(all_Genus_BrayCurtis, k = 4)
pcoa_genus_d36_d50_BC_eigen <- as.data.frame(pcoa_genus_d36_d50_BC$eig*100/sum(pcoa_genus_d36_d50_BC$eig))
pcoa_genus_d36_d50_BC <- as.data.frame(pcoa_genus_d36_d50_BC$points)
pcoa_genus_d36_d50_BC <- rownames_to_column(pcoa_genus_d36_d50_BC, var = "Var1")
pcoa_genus_d36_d50_BC <- merge(pcoa_genus_d36_d50_BC, metadata_d36_d50[,c("Sample", "Day", "Treatment")], by.x = "Var1", by.y = "Sample")

###Plotting day 36 and day 50 time points, and adding visualization aides including color and ellipses
ggplot() + 
  geom_point(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 36),], aes(x = V1, y = V2, color = "#949494")) +
  stat_ellipse(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 36),], aes(x = V1, y = V2, color = "#949494"), type = "norm", level = 0.95) +
  geom_point(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 50),], aes(x = V1, y = V2, color = Treatment)) +
  stat_ellipse(data = pcoa_genus_d36_d50_BC[which(pcoa_genus_d36_d50_BC$Day == 50),], aes(x = V1, y = V2, color = Treatment), type = "norm", level = 0.95) +
  scale_color_manual(
    values = c("PBS" = "#000000", "Sb_OG539" = "#F78F26", "Sb_haPD1" = "#409CDC", "aPDL-1" = "#AC55EA"),
    breaks = c("PBS", "Sb_OG539", "Sb_haPD1", "aPDL-1")) +  
  xlab("PCoA 1 (40.3%)") + 
  ylab("PCoA 2 (14.2%)") +
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
  xlab("PCoA 1 (39.7%)") + 
  ylab("PCoA 2 (23.7%)") +
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
  xlab("PCoA 1 (32.6.7%)") + 
  ylab("PCoA 2 (25.0%)") +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "aPDL-1")
p



###PERMANOVA analysis comparing treatment groups at day 50###
#Comparing PBS and Sb_haPD-1 at d50#
adonis_genus_abund <- merge(genus_abund_d36_d50, metadata_d36_d50[,c(1,4,6)], by.x = "row.names", by.y = "Sample")
row.names(adonis_genus_abund) <- adonis_genus_abund$Row.names
adonis_genus_abund <- adonis_genus_abund[,c(-1)]

DaysToKeep <- c("50")
PBSvSb_haPD1 <- c("PBS", "Sb_haPD1")
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Day %in% DaysToKeep, ]
adonis_genus_abund <- adonis_genus_abund[adonis_genus_abund$Treatment %in% PBSvSb_haPD1, ]
PBSvSb_haPD1_result <- adonis2(adonis_genus_abund[, c(-81, -82)] ~ adonis_genus_abund$Treatment,
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
metadata_d36_d50 <- read.csv("~/Desktop/metadata_d36_d50.csv", header = TRUE)
genus_abund_d36_d50 <- read.csv("~/Desktop/genus_abund_d36_d50.csv", row.names=1, check.names=FALSE)

#Turning raw read count into relative abundance
genus_abund_d36_d50 <- genus_abund_d36_d50[rowSums(genus_abund_d36_d50[, -1]) > 0, ]
genus_abund_d36_d50 <- genus_abund_d36_d50/rowSums(genus_abund_d36_d50, na.rm = TRUE)

#Turning .csv metadata to .txt
output_file <- "input_metadata"
write.table(metadata_d36_d50, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
metadata_d36_d50 <- read.table(output_file, header = TRUE)
rownames(metadata_d36_d50) <- metadata_d36_d50$Samples
metadata_d36_d50 <- metadata_d36_d50[, -1]

#Specifying day(s)/timepoint(s) for intput_data (ONLY for single experiments)
#timepoint <- c("36", "50") #Change this
input_metadata <- metadata_d36_d50

#choosing PBS as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="PBS"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "PBS_d36vd50", 
                    fixed_effects  = c("Day"),
                    reference      = c("Day,36"),
                    random_effects = NULL)

#choosing Sb_haPD1 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="Sb_haPD1"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "SbhaPD1_d36vd50", 
                    fixed_effects  = c("Day"),
                    reference      = c("Day,36"),
                    random_effects = NULL)

#choosing Sb_OG539 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="Sb_OG539"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "SbOG539_d36vd50", 
                    fixed_effects  = c("Day"),
                    reference      = c("Day,36"),
                    random_effects = NULL)

#choosing aPDL-1 as the treatment group to compare across day 36 and day 50
input_metadata_treatment <- subset(input_metadata, (Treatment =="aPDL-1"))

#Running MaAsLin2
fit_data = Maaslin2(input_data     = genus_abund_d36_d50, 
                    input_metadata = input_metadata_treatment, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "aPDL1_d36vd50", 
                    fixed_effects  = c("Day"),
                    reference      = c("Day,36"),
                    random_effects = NULL)
