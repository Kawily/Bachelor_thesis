library(readr)
library(openxlsx)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)

age_groups <- read_xlsx("Altersgruppen_Sortiert.xlsx")

age_groups <- as.data.frame(t(age_groups))

Counts_brain <- read.csv("FPKM_brain.csv")

age_minus05 <- select(Counts_brain , c('SRR1554537' , 'SRR1554538', 'SRR1554541'
                                 ,'SRR1554566'
                                 ,'SRR1554567'
                                 ,'SRR1554568'
                                 ,'SRR2071348'
                                 ,'SRR2071349'
                                 ,'SRR2071352'
                                 ,'SRR2071377','SRR2071378' , 'SRR2071379'))

age_60 <- select(Counts_brain, c('SRR1554533',
                           'SRR1554555',
                           'SRR1554557',
                           'SRR1554560',
                           'SRR1554563',
                           'SRR2071341',
                           'SRR2071366',
                           'SRR2071368',
                           'SRR2071371',
                           'SRR2071374'))

age_40 <- select(Counts_brain, c('SRR1554534',
                           'SRR1554535',
                           'SRR1554536',
                           'SRR1554561',
                           'SRR2071345',
                           'SRR2071346',
                           'SRR2071347',
                           'SRR2071372'))

age_4 <- select(Counts_brain, c('SRR1554542',
                          'SRR1554543',
                          'SRR1554544',
                          'SRR1554545',
                          'SRR1554546',
                          'SRR1554549',
                          'SRR1554551',
                          'SRR1554552',
                          'SRR1554553',
                          'SRR1554554',
                          'SRR1554564',
                          'SRR1554565',
                          'SRR2071353',
                          'SRR2071354',
                          'SRR2071355',
                          'SRR2071356',
                          'SRR2071357',
                          'SRR2071360',
                          'SRR2071362',
                          'SRR2071363',
                          'SRR2071364',
                          'SRR2071365',
                          'SRR2071375',
                          'SRR2071376'))

age_16 <- select(Counts_brain, c('SRR1554540',
                           'SRR1554547',
                           'SRR1554548',
                           'SRR1554550',
                           'SRR1554558',
                           'SRR1554562',
                           'SRR2071351',
                           'SRR2071358',
                           'SRR2071359',
                           'SRR2071361',
                           'SRR2071373'))


sd_05 <- as.data.frame(apply(age_minus05, 1, sd))
sd_4 <- as.data.frame(apply(age_4, 1, sd))
sd_16 <- as.data.frame(apply(age_16, 1, sd))
sd_40 <- as.data.frame(apply(age_40, 1, sd))
sd_60 <- as.data.frame(apply(age_60, 1, sd))


colSums(na.omit(sd_05))
colSums(na.omit(sd_4))
colSums(na.omit(sd_16))
colSums(na.omit(sd_40))
colSums(na.omit(sd_60))

### Normalization!!!

Cut_off_point <- 0.001
Test_1 <-  Counts_brain[rowSums(Counts_brain < Cut_off_point) == 0, ]
Test_2 <- Counts_brain[Counts_brain >=Cut_off_point,]
Na_free_data <- na.omit(Test_2)


#write.csv(NA_free_data,"your path")

Colsumsresults <- as.data.frame(Test_1[-1])
Ready_for_formatting <- Test_1[-1]


Ready_for_formatting <- na.omit(Ready_for_formatting)


sums <- apply(Ready_for_formatting, 2, sum) # calculate the sum of each column
normalized_data <- sweep(Ready_for_formatting, 2, 1000000/sums, "*") # normalize each column to sum up to 1000000


colSums(normalized_data)


age_minus05_n <- select(normalized_data , c('SRR1554537' , 'SRR1554538', 'SRR1554541'
                                            ,'SRR1554566'
                                            ,'SRR1554567'
                                            ,'SRR1554568'
                                            ,'SRR2071348'
                                            ,'SRR2071349'
                                            ,'SRR2071352'
                                            ,'SRR2071377','SRR2071378' , 'SRR2071379'))

age_60_n <- select(normalized_data, c('SRR1554533',
                                      'SRR1554555',
                                      'SRR1554557',
                                      'SRR1554560',
                                      'SRR1554563',
                                      'SRR2071341',
                                      'SRR2071366',
                                      'SRR2071368',
                                      'SRR2071371',
                                      'SRR2071374'))

age_40_n <- select(normalized_data, c('SRR1554534',
                                      'SRR1554535',
                                      'SRR1554536',
                                      'SRR1554561',
                                      'SRR2071345',
                                      'SRR2071346',
                                      'SRR2071347',
                                      'SRR2071372'))

age_4_n <- select(normalized_data, c('SRR1554542',
                                     'SRR1554543',
                                     'SRR1554544',
                                     'SRR1554545',
                                     'SRR1554546',
                                     'SRR1554549',
                                     'SRR1554551',
                                     'SRR1554552',
                                     'SRR1554553',
                                     'SRR1554554',
                                     'SRR1554564',
                                     'SRR1554565',
                                     'SRR2071353',
                                     'SRR2071354',
                                     'SRR2071355',
                                     'SRR2071356',
                                     'SRR2071357',
                                     'SRR2071360',
                                     'SRR2071362',
                                     'SRR2071363',
                                     'SRR2071364',
                                     'SRR2071365',
                                     'SRR2071375',
                                     'SRR2071376'))

age_16_n <- select(normalized_data, c('SRR1554540',
                                      'SRR1554547',
                                      'SRR1554548',
                                      'SRR1554550',
                                      'SRR1554558',
                                      'SRR1554562',
                                      'SRR2071351',
                                      'SRR2071358',
                                      'SRR2071359',
                                      'SRR2071361',
                                      'SRR2071373'))


sd_05_n <- as.data.frame(apply(age_minus05_n, 1, sd))
sd_4_n <- as.data.frame(apply(age_4_n, 1, sd))
sd_16_n <- as.data.frame(apply(age_16_n, 1, sd))
sd_40_n <- as.data.frame(apply(age_40_n, 1, sd))
sd_60_n <- as.data.frame(apply(age_60_n, 1, sd))


colSums(na.omit( sd_05_n))
colSums(na.omit(sd_4_n))
colSums(na.omit(sd_16_n))
colSums(na.omit(sd_40_n))
colSums(na.omit(sd_60_n))

nullo <- log2(sd_05_n)
vier <- log2(sd_4_n)
sechzehn <- log2(sd_16_n)
vierzig <- log2(sd_40_n)
sechzig <- log2(sd_60_n)


age_groups_merged <- cbind(id = 1:23144, nullo, vier, sechzehn, vierzig, sechzig )

colnames(age_groups_merged) <- c("id" ,"before birth", "0,25-4 years", "14-18 years","40-45 years", "66-74 years")
data_long <- melt(age_groups_merged, id.vars = "id")

ggplot(data_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Groups", y = "log2", title = "Variation compared", fill = "Age groups")
theme(axis.title.y = element_text(size = 28))


