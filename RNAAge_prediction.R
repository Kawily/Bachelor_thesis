library(RColorBrewer)
require(dplyr)
library(RNAAgeCalc)
library(readr)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(ggsci)
library(ggthemes)
library(stringr)
library(janitor)
library(tidyverse)
library(dplyr)
library(tidyr)

gene_lengh <- (x <- matrix(data = 1000L, nrow = 58736, ncol = 1))
row.names(gene_lengh) <- rawnames
genl <-as.numeric(gene_lengh)


rawcounts <- read_tsv("featurecounts.count.gene.tsv")
rawcounts <- as.data.frame(rawcounts)

raw <- rawcounts[,-1]
rawnames <- rawcounts[,1]
as.vector(rawnames)
row.names(raw) <- rawnames
colnames(raw) <- c("Untreated60a","Untreated60b","Untreated60c","DMSO60a","DMSO60b","DMSO60c","Dox60a","Dox60b","Dox60c","Untreated90a","Untreated90b","Untreated90c","DMSO90a","DMSO90b","DMSO90c","Dox90a","Dox90b","Dox90c")


cpm <- read_tsv("featurecounts.cpm.gene.tsv")
cpm <- as.data.frame(cpm)

cpm_as_FPKM <- cpm[,-1]
cpmnames <- cpm[,1]
as.vector(cpmnames)
row.names(cpm_as_FPKM) <- cpmnames
colnames(cpm_as_FPKM) <- c("Untreated60a","Untreated60b","Untreated60c","DMSO60a","DMSO60b","DMSO60c","Dox60a","Dox60b","Dox60c","Untreated90a","Untreated90b","Untreated90c","DMSO90a","DMSO90b","DMSO90c","Dox90a","Dox90b","Dox90c")

##Predict Age
cpm_predict_age <- predict_age(exprdata = cpm_as_FPKM, exprtype = "FPKM")
raw_predict_age <- predict_age(exprdata = raw, genelength =  genl, exprtype = "counts")


cpm_predict_age$Group <- rep(1:(nrow(cpm_predict_age)/3), each = 3)

cpm_predict_age$Group <- factor(cpm_predict_age$Group)

num_groups <- length(unique(cpm_predict_age$Group))

group_colors <- brewer.pal(num_groups, "Set3")

ggplot(cpm_predict_age, aes(x = rownames(cpm_predict_age), y = RNAAge, group = Group, color = Group)) +
  geom_point(size = 4, shape = 16, color = "black") +
  geom_point(size = 3) +
  labs(title = "Across Tissue Age", x = "Sample", y = "RNAAge") +
  scale_color_manual(values = group_colors, 
                     breaks = levels(cpm_predict_age$Group),
                     labels = c("Untreated60","DMSO60","DOX60","Untreated90", "DMSO90","DOX90")) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis labels
#pdf 7 x 4


raw_brain_age <- predict_age(
  raw,
  tissue = c("brain"),
  exprtype = c("counts"),
  idtype = c("ENSEMBL"),
  stype = c("caucasian"),
  signature = NULL,
  genl,
  chronage = NULL,
  maxp = NULL
)


Annotation <- read.csv("parameter.orig.csv")
Sample_treatments <- Annotation[,1]
raw_brain_age$treatment <- Sample_treatments


ggplot(raw_brain_age, aes(x = treatment, y = RNAAge, fill = treatment, shape = treatment)) + 
  geom_point(size = 4) +
  labs(x = "treatment", y = "RNAAge", fill = "Treatment", shape = "Treatment") +
  theme_minimal() +
  scale_fill_manual(values = c("darkred", "darkorange", "darkgreen", "navy", "purple", "darkcyan")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 19)) +
  theme(axis.text.y = element_text(size = 12))+
  ggtitle("RNAAge Brain Tissue Prediction Caucasian")

#12 x 6

raw_brain_age_all <- predict_age(
  raw,
  tissue = c("brain"),
  exprtype = c("counts"),
  idtype = c("ENSEMBL"),
  stype = c("all"),
  signature = NULL,
  genl,
  chronage = NULL,
  maxp = NULL
)


Annotation <- read.csv("parameter.orig.csv")
Sample_treatments <- Annotation[,1]
raw_brain_age_all$treatment <- Sample_treatments


ggplot(raw_brain_age_all, aes(x = treatment, y = RNAAge, fill = treatment, shape = treatment)) + 
  geom_point(size = 4) +
  labs(x = "treatment", y = "RNAAge", fill = "Treatment", shape = "Treatment") +
  theme_minimal() +
  scale_fill_manual(values = c("darkred", "darkorange", "darkgreen", "navy", "purple", "darkcyan")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 19)) +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("RNAAge Brain Tissue Prediction All")


library(readxl)
Deseqfun <- read_excel("Deseqfun.xlsx")


str(Deseqdata)
Deseqdata <- as.data.frame(Deseqfun)

Deseqdata_1 <- Deseqdata[,-1]
Deseqnames <- Deseqdata[,1]
as.vector(Deseqnames)
row.names(Deseqdata_1) <- Deseqnames

Deseqcharacter <- as.character(Deseqdata_1)

Deseqtry <- predict_age(
  raw,
  tissue = c("brain"),
  exprtype = c("counts"),
  idtype = c("ENSEMBL"),
  stype = c("caucasian"),
  signature = NULL,
  genl,
  chronage = NULL,
  maxp = NULL
)

Annotation <- read.csv("parameter.orig.csv")
Sample_treatments <- Annotation[,1]
Deseqtry$treatment <- Sample_treatments


ggplot(Deseqtry, aes(x = treatment, y = RNAAge, fill = treatment, shape = treatment)) + 
  geom_point(size = 4) +
  labs(x = "treatment", y = "RNAAge", fill = "Treatment", shape = "Treatment") +
  theme_minimal() +
  scale_fill_manual(values = c("darkred", "darkorange", "darkgreen", "navy", "purple", "darkcyan")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 19)) +
  theme(axis.text.y = element_text(size = 12)) +
  ggtitle("RNAAge Brain Tissue Prediction with own DESeq2")

Deseqtry %in% raw_brain_age

Lisas_Annotation <- read.csv("parameter.orig.csv")
Sample_treatments <- Lisas_Annotation[,1]
Deseqtry$treatment <- Sample_treatments
