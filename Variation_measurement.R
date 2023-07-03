#Variation measurement clean
library(readr)
library(openxlsx)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

setwd("C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit/Lisas transcriptom/CO.UMIcor.iiRNAseq_ii.GRCh38.Progerin_20220407/data/EXPR_summarized")
cpm <- read_tsv("featurecounts.cpm.gene.tsv")
cpm <- as.data.frame(cpm)
raw1 <- cpm[,-1]
rawnames1 <- cpm[,1]
as.vector(rawnames1)
row.names(raw1) <- rawnames1

row_sd2 <- apply(raw1, 1, sd)


day60untreated <- select(raw1, c('019','020','021'))
day60DMSO <- select(raw1, c('022','023','024'))
day60DOX <- select(raw1, c('025','026','027'))

day90untreated <- select(raw1, c('028','029','030'))
day90DMSO <- select(raw1, c('031','032','033'))
day90DOX <- select(raw1, c('034','035','036'))


sd_day60untreated <- as.data.frame(apply(day60untreated, 1, sd))
sd_day60DMSO <- as.data.frame(apply(day60DMSO, 1, sd))
sd_day60DOX <- as.data.frame(apply(day60DOX, 1, sd))

sd_day90untreated <- as.data.frame(apply(day90untreated, 1, sd))
sd_day90DMSO <- as.data.frame(apply(day90DMSO, 1, sd))
sd_day90DOX  <- as.data.frame(apply(day90DOX, 1, sd))


#colSums(sd_day60untreated)
#colSums(sd_day60DMSO)
#colSums(sd_day60DOX)

#colSums(sd_day90untreated)
#colSums(sd_day90DMSO)
#colSums(sd_day90DOX)

sdlog <- data.frame(sd_day60DMSO,sd_day90DOX)
sdlog1 <- sd_day90DOX/sd_day60DMSO
sdlog2 <- log2(sdlog1)
str(sdlog2)
as.numeric(unlist (sdlog2))
sdlog3 <- unlist(sdlog2)
sdlog3[!is.finite(sdlog3)] <- NA

day_90_DOX_vs_day_60_DMSO  <- sdlog3
colSums(sdlog3)


#table(sd_day60DMSO)
#table(sd_day90DOX)


plot(day_90_DOX_vs_day_60_DMSO,)
abline(h=0, col="blue")
title(" Day 90 Dox / Day 60 DMSO", size = 28)
xlab("Genes")
ylab("log2")


#as these are the same processes as before it is copied with the same names!
#Make sure to clean the scripts data before or just know that it will be overwritten

# Delete Zero values so the variation becomes more comparable
setwd("C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit/Lisas transcriptom/CO.UMIcor.iiRNAseq_ii.GRCh38.Progerin_20220407/data/EXPR_summarized")

cpm_on_zero<- read_tsv("featurecounts.cpm.gene.tsv")
cpm_on_zero[cpm_on_zero == 0] <- NA
cpm_no_zero <- na.omit(cpm_on_zero)
setwd("C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit")
write.csv(cpm_no_zero,"C:/Users/kawil/OneDrive - uibk.ac.at/WS22_23/Bachelorarbeit/cpm_clean.csv")

#Standardisation and norming CPM again
Standard_cpm <- cpm_no_zero[-1]
Colsumsresults <- colSums(Standard_cpm)
Colsumsresults <- as.data.frame(Colsumsresults)

Colstand <- 1000000 / Colsumsresults
colstand_data <- as.data.frame(t(Colstand))

cpm_norm<- mapply(`*`, colstand_data, Standard_cpm)
cpm_norm <-  as.data.frame(cpm_norm)
colSums(cpm_norm)


row_sd2 <- apply(cpm_norm, 1, sd)

day60untreated <- select(cpm_norm, c('019','020','021'))
day60DMSO <- select(cpm_norm, c('022','023','024'))
day60DOX <- select(cpm_norm, c('025','026','027'))

day90untreated <- select(cpm_norm, c('028','029','030'))
day90DMSO <- select(cpm_norm, c('031','032','033'))
day90DOX <- select(cpm_norm, c('034','035','036'))




sd_day60untreated <- as.data.frame(apply(day60untreated, 1, sd))
sd_day60DMSO <- as.data.frame(apply(day60DMSO, 1, sd))
sd_day60DOX <- as.data.frame(apply(day60DOX, 1, sd))

sd_day90untreated <- as.data.frame(apply(day90untreated, 1, sd))
sd_day90DMSO <- as.data.frame(apply(day90DMSO, 1, sd))
sd_day90DOX  <- as.data.frame(apply(day90DOX, 1, sd))

#colSums(sd_day60untreated)
#colSums(sd_day60DMSO)
#colSums(sd_day60DOX)

#colSums(sd_day90untreated)
#colSums(sd_day90DMSO)
#colSums(sd_day90DOX)

sd_all <- data.frame(sd_day60DMSO,sd_day60DOX, sd_day60untreated, sd_day90DMSO, sd_day90DOX,sd_day90untreated)

#mean(sd_all$apply.day60DMSO..1..sd.)
#mean(sd_all$apply.day90DOX..1..sd.)

#median(sd_all$apply.day60DMSO..1..sd.)
#median(sd_all$apply.day90DOX..1..sd.)


#Example code

#Dox60 DMSO60
DOX60_DMSO60_div <- sd_day60DOX/sd_day60DMSO
DOX60_DMSO60_log2 <- log2(DOX60_DMSO60_div)
as.numeric(unlist (DOX60_DMSO60_log2))
DOX60_DMSO60 <- unlist(DOX60_DMSO60_log2)

DOX60_DMSO60_sorted <- as.data.frame(DOX60_DMSO60) %>% arrange(desc(DOX60_DMSO60))
DOX60_DMSO60_subset <- as.data.frame(DOX60_DMSO60_sorted[seq(10, nrow(DOX60_DMSO60_sorted), by = 10),])


DOX60_DMSO60_subset_n <-  as.numeric(unlist(DOX60_DMSO60_subset))

DOX60_DMSO60_randomized <- sample(DOX60_DMSO60_subset_n)

DOX60_untreated60_div <- sd_day60DOX / sd_day60untreated
DOX60_untreated60_log2 <- log2(DOX60_untreated60_div)
DOX60_untreated60 <- as.numeric(unlist(DOX60_untreated60_log2))
DOX60_untreated60_sorted <- as.data.frame(DOX60_untreated60) %>% arrange(desc(DOX60_untreated60))
DOX60_untreated60_subset <- as.data.frame(DOX60_untreated60_sorted[seq(10, nrow(DOX60_untreated60_sorted), by = 10),])
DOX60_untreated60_subset_n <- as.numeric(unlist(DOX60_untreated60_subset))
DOX60_untreated60_randomized <- sample(DOX60_untreated60_subset_n)

DOX90_DMSO90_div <- sd_day90DOX / sd_day90DMSO
DOX90_DMSO90_log2 <- log2(DOX90_DMSO90_div)
DOX90_DMSO90 <- as.numeric(unlist(DOX90_DMSO90_log2))
DOX90_DMSO90_sorted <- as.data.frame(DOX90_DMSO90) %>% arrange(desc(DOX90_DMSO90))
DOX90_DMSO90_subset <- as.data.frame(DOX90_DMSO90_sorted[seq(10, nrow(DOX90_DMSO90_sorted), by = 10),])
DOX90_DMSO90_subset_n <- as.numeric(unlist(DOX90_DMSO90_subset))
DOX90_DMSO90_randomized <- sample(DOX90_DMSO90_subset_n)

# Example code for DOX90 untreated90
DOX90_untreated90_div <- sd_day90DOX / sd_day90untreated
DOX90_untreated90_log2 <- log2(DOX90_untreated90_div)
DOX90_untreated90 <- as.numeric(unlist(DOX90_untreated90_log2))
DOX90_untreated90_sorted <- as.data.frame(DOX90_untreated90) %>% arrange(desc(DOX90_untreated90))
DOX90_untreated90_subset <- as.data.frame(DOX90_untreated90_sorted[seq(10, nrow(DOX90_untreated90_sorted), by = 10),])
DOX90_untreated90_subset_n <- as.numeric(unlist(DOX90_untreated90_subset))
DOX90_untreated90_randomized <- sample(DOX90_untreated90_subset_n)

# Example code for DMSO90 untreated90
DMSO90_untreated90_div <- sd_day90DMSO / sd_day90untreated
DMSO90_untreated90_log2 <- log2(DMSO90_untreated90_div)
DMSO90_untreated90 <- as.numeric(unlist(DMSO90_untreated90_log2))
DMSO90_untreated90_sorted <- as.data.frame(DMSO90_untreated90) %>% arrange(desc(DMSO90_untreated90))
DMSO90_untreated90_subset <- as.data.frame(DMSO90_untreated90_sorted[seq(10, nrow(DMSO90_untreated90_sorted), by = 10),])
DMSO90_untreated90_subset_n <- as.numeric(unlist(DMSO90_untreated90_subset))
DMSO90_untreated90_randomized <- sample(DMSO90_untreated90_subset_n)

DOX90_DMSO60_div <- sd_day90DOX / sd_day60DMSO
DOX90_DMSO60_log2 <- log2(DOX90_DMSO60_div)
DOX90_DMSO60 <- as.numeric(unlist(DOX90_DMSO60_log2))
DOX90_DMSO60_sorted <- as.data.frame(DOX90_DMSO60) %>% arrange(desc(DOX90_DMSO60))
DOX90_DMSO60_subset <- as.data.frame(DOX90_DMSO60_sorted[seq(10, nrow(DOX90_DMSO60_sorted), by = 10),])
DOX90_DMSO60_subset_n <- as.numeric(unlist(DOX90_DMSO60_subset))
DOX90_DMSO60_randomized <- sample(DOX90_DMSO60_subset_n)

# Example code for DOX60 DOX90
DOX60_DOX90_div <- sd_day60DOX / sd_day90DOX
DOX60_DOX90_log2 <- log2(DOX60_DOX90_div)
DOX60_DOX90 <- as.numeric(unlist(DOX60_DOX90_log2))
DOX60_DOX90_sorted <- as.data.frame(DOX60_DOX90) %>% arrange(desc(DOX60_DOX90))
DOX60_DOX90_subset <- as.data.frame(DOX60_DOX90_sorted[seq(10, nrow(DOX60_DOX90_sorted), by = 10),])
DOX60_DOX90_subset_n <- as.numeric(unlist(DOX60_DOX90_subset))
DOX60_DOX90_randomized <- sample(DOX60_DOX90_subset_n)

# Example code for DMSO60 DMSO90
DMSO60_DMSO90_div <- sd_day60DMSO / sd_day90DMSO
DMSO60_DMSO90_log2 <- log2(DMSO60_DMSO90_div)
DMSO60_DMSO90 <- as.numeric(unlist(DMSO60_DMSO90_log2))
DMSO60_DMSO90_sorted <- as.data.frame(DMSO60_DMSO90) %>% arrange(desc(DMSO60_DMSO90))
DMSO60_DMSO90_subset <- as.data.frame(DMSO60_DMSO90_sorted[seq(10, nrow(DMSO60_DMSO90_sorted), by = 10),])
DMSO60_DMSO90_subset_n <- as.numeric(unlist(DMSO60_DMSO90_subset))
DMSO60_DMSO90_randomized <- sample(DMSO60_DMSO90_subset_n)

# Example code for untreated60 untreated90
untreated60_untreated90_div <- sd_day60untreated / sd_day90untreated
untreated60_untreated90_log2 <- log2(untreated60_untreated90_div)
untreated60_untreated90 <- as.numeric(unlist(untreated60_untreated90_log2))
untreated60_untreated90_sorted <- as.data.frame(untreated60_untreated90) %>% arrange(desc(untreated60_untreated90))
untreated60_untreated90_subset <- as.data.frame(untreated60_untreated90_sorted[seq(10, nrow(untreated60_untreated90_sorted), by = 10),])
untreated60_untreated90_subset_n <- as.numeric(unlist(untreated60_untreated90_subset))
untreated60_untreated90_randomized <- sample(untreated60_untreated90_subset_n)

# Example code for DMSO60 Untreated60
DMSO60_Untreated60_div <- sd_day60DMSO / sd_day60untreated
DMSO60_Untreated60_log2 <- log2(DMSO60_Untreated60_div)
DMSO60_Untreated60 <- as.numeric(unlist(DMSO60_Untreated60_log2))
DMSO60_Untreated60_sorted <- as.data.frame(DMSO60_Untreated60) %>% arrange(desc(DMSO60_Untreated60))
DMSO60_Untreated60_subset <- as.data.frame(DMSO60_Untreated60_sorted[seq(10, nrow(DMSO60_Untreated60_sorted), by = 10),])
DMSO60_Untreated60_subset_n <- as.numeric(unlist(DMSO60_Untreated60_subset))
DMSO60_Untreated60_randomized <- sample(DMSO60_Untreated60_subset_n)


# Set the number of plots to create
num_plots1 <- 8

#num_plots <- 3

# Set the different data for each plot
graph_list1 <- list(
  data.frame(value = DOX60_DMSO60_randomized),
  data.frame(value = DOX90_DMSO90_randomized),
  data.frame(value = DOX90_DMSO60_randomized),
  data.frame(value = DOX60_untreated60_randomized),
  data.frame(value = DOX90_untreated90_randomized),
  data.frame(value = DMSO60_Untreated60_randomized),
  data.frame(value = DMSO90_untreated90_randomized),
  data.frame(value = DOX90_DMSO60_randomized)
)
title_list1 <- c(
  "DOX60/DMSO60",
  "DOX90/DMSO90",
  "DOX90/DMSO60",
  "DOX60/Untreated60",
  "DOX90/Untreated90",
  "DMSO60/Untreated60",
  "DMSO90/Untreated90",
  "Dox90/DMSO60"
)
hline_list <- c(
  mean(DOX60_DMSO60),
  mean(DOX90_DMSO90),
  mean(DOX90_DMSO90),
  mean(DOX60_untreated60),
  mean(DOX90_untreated90),
  mean(DMSO60_Untreated60),
  mean(DMSO90_untreated90),
  mean(DOX90_DMSO60)
  
)

# Create a list to store the plots
plot_list1 <- list()

num_plots1 <- 8
# Loop through and create the plots
for (i in 1:num_plots1) {
  
  # Create the plot
  my_plot <- ggplot(data = graph_list1[[i]], aes(x = 1:length(value), y = value)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 0.8) +
    geom_hline(yintercept = hline_list[i], linetype = "dashed", color = "blue", size = 1.2) +
    lims(y = c(-6, 6)) +
    labs(title = title_list1[i],
         x = "every 10th dot",
         y = "log2") +
    theme(plot.title = element_text(size = 18))
  
  
  
  # Add the plot to the list
  plot_list1[[i]] <- my_plot
}

# Plot the list of plots using grid.arrange
gridExtra::grid.arrange(grobs = plot_list1, nrow = 4, ncol = 2)

print(plot_list1)
#Controls

graph_list2 <- list(
  data.frame(value = DOX60_DOX90_randomized),
  data.frame(value = DMSO60_DMSO90_randomized),
  data.frame(value = untreated60_untreated90_randomized)
)
title_list2 <- c(
  "DOX60/DOX90",
  "DMSO60/DMSO90",
  "Untreated60/Untreated90"
)
hline_list1 <- c(
  mean(DOX60_DOX90),
  mean(DMSO60_DMSO90),
  mean(untreated60_untreated90)
)

# Create a list to store the plots
plot_list2 <- list()

num_plots2 <- 3
# Loop through and create the plots
for (i in 1:num_plots2) {
  
  # Create the plot
  my_plot1 <- ggplot(data = graph_list2[[i]], aes(x = 1:length(value), y = value)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 0.8) +
    geom_hline(yintercept = hline_list1[i], linetype = "dashed", color = "blue", size = 1.2) +
    lims(y = c(-6, 6)) +
    labs(title = title_list2[i],
         x = "every 10th dot",
         y = "log2") +
    theme(plot.title = element_text(size = 18))
  
  
  
  # Add the plot to the list
  plot_list2[[i]] <- my_plot1
}

# Plot the list of plots using grid.arrange
gridExtra::grid.arrange(grobs = plot_list2, nrow = 3, ncol = 1)

print(plot_list2)

#Vioplots

Log2_Vioplots <- data.frame(DOX60_DMSO60,DOX60_untreated60,DMSO60_Untreated60,DOX90_DMSO90, DOX90_untreated90, DMSO90_untreated90, DOX90_DMSO60)
Log2_Vioplots2 <-data.frame(id = 1:13853, DOX60_DMSO60,DOX60_untreated60,DMSO60_Untreated60,DOX90_DMSO90, DOX90_untreated90, DMSO90_untreated90, DOX90_DMSO60)
Log2_Vioplots3 <-data.frame(id = 1:13853, DOX60_DOX90, DMSO60_DMSO90, untreated60_untreated90)

data_long <- melt(Log2_Vioplots2, id.vars = "id")


vio1 <- ggplot(data_long, aes(x = variable, y = value, fill = variable)) +
  geom_violin() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

vio2 <- vio1 + labs(x = "Sample", y = "log2", fill = "variation1/variation2") +
  labs(title = "Variation compared") +
  theme(plot.title = element_text(size = 25))

print(vio2)

data_long1 <- melt(Log2_Vioplots3, id.vars = "id")


vio3 <- ggplot(data_long1, aes(x = variable, y = value, fill = variable)) +
  geom_violin() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

vio4 <- vio3 +labs(x = "Sample", y = "log2", fill= "Day60/Day90")

print(vio4)
#Scatterplot

Day90DOX_log <-log2(sd_day90DOX)
day60DMSO_log <- log2(sd_day60DMSO)


df <- data.frame(x = Day90DOX_log, y = day60DMSO_log)

names(df) <-c("Day90DOX_log2", "Day60DMSO_log2")

ggplot(df, aes(x = Day90DOX_log2, y = Day60DMSO_log2)) + 
  geom_point(color = "#0072B2", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "#D55E00", size = 1.5) +
  geom_density_2d(alpha = 0.5, color = "black") +
  xlab("Day 90 DOX log2") + ylab("Day 60 DMSO log2") +
  ggtitle("Scatterplot of DOX90 and DMSO60") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  coord_fixed(ratio = 1) +
  xlim(-5, 10) + ylim(-5, 10)


