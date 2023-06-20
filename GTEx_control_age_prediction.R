library(RNAAgeCalc)
library(readr)
library(ggplot2)
library(tidyverse)
library(tibble)
library(dplyr)
library(ggpubr)
library(extrafont)
library(ggsci)
library(ggthemes)
library(data.table)
library('recount')


url <- download_study('SRP012682')
load(file.path('SRP012682', 'rse_gene.Rdata'))

## Scale counts
rse <- scale_counts(rse_gene)


counts_brain <- rse@assays@data@listData[["counts"]]
write.csv(counts_brain,"(GTEX_data1.csv" )
rcounts <- read.csv("(GTEX_data1.csv")
new_data <- rcounts[ c( "SRR1069188", "SRR1071289",
                        "SRR1071880",
                        "SRR1072178",
                        "SRR1072367",
                        "SRR1072504",
                        "SRR1081741",
                        "SRR1084842",
                        "SRR1085495",
                        "SRR1310136",
                        "SRR1313642",
                        "SRR1316254",
                        "SRR1320071",
                        "SRR1335400")]

namesdata <- rcounts [1]

rawcounts <- read_xlsx("Book1.xlsx")

rnamesversionfree <-pivot_longer(namesdata, cols = everything()) %>%
  separate(value, into = c('gene_id', 'beta'), sep = "\\.") %>%
  select(-name)

rnames <- rnamesversionfree[,1]
rnamesdata <- as.data.frame(rnames)
rawnamesdataunique <-make.names(rnamesdata[,1], unique = TRUE)
row.names(new_data) <- rawnamesdataunique


cp <- predict_age(
  new_data,
  tissue = c("brain"),
  exprtype = c("counts"),
  idtype = c("ENSEMBL"),
  stype = c("all"),
  signature = NULL,
  genelength = NULL,
  chronage = NULL,
  maxp = NULL
)

write.csv(cp,"GTEX_RNAAGE.csv" )

Agerange <- t(data.frame(35,55,45,65,65,65,55,65,65,45,45,75,75,25 ))
df <- rownames_to_column(Agerange, var = "X")

Age_r <- merge(Agerange,Sample_names)
Sample_names <- t(data.frame("SRR1335400",
                             "SRR1081741",
                             "SRR1084842",
                             "SRR1085495",
                             "SRR1069188",
                             "SRR1071289",
                             "SRR1071880",
                             "SRR1072178",
                             "SRR1072367",
                             "SRR1072504",
                             "SRR1310136",
                             "SRR1313642",
                             "SRR1316254",
                             "SRR1320071"
))

df3 <- (data.frame(names =c("SRR1069188"
                            ,"SRR1071289"
                            ,"SRR1071880"
                            ,"SRR1072178"
                            ,"SRR1072367"
                            ,"SRR1072504"
                            ,"SRR1081741"
                            ,"SRR1084842"
                            ,"SRR1085495",
                            "SRR1310136",
                            "SRR1313642"
                            ,"SRR1316254"
                            ,"SRR1320071"
                            ,"SRR1335400"
                            
                      ), age =c(65,65,55,65,65,45,55,45,65,45,75,75,25,35 ),
                   RNAAGE = c(57.01823,
                   
                   52.93410
                   ,
                   41.68924
                   ,
                   56.01721
                   ,
                   55.07727
                   ,
                   43.92025
                   ,
                   47.60838
                   ,
                   42.83530
                   ,
                   59.59275
                   ,
                   43.84962
                   ,
                   54.96432
                   ,
                   61.57459
                   ,
                   22.29486
                   ,
                   41.56827)))



df2 <- df1[-1]
rnames1 <- df1[,1]
rnamesdata1 <- as.data.frame(rnames1)

row.names(df2) <- rnamesdata1[1]


Happy <- data.frame(cp,c(
  "SRR1081741",
  "SRR1084842",
  "SRR1085495",
  "SRR1069188",
  "SRR1071289",
  "SRR1071880",
  "SRR1072178",
  "SRR1072367",
  "SRR1072504",
  "SRR1310136",
  "SRR1313642",
  "SRR1316254",
  "SRR1320071",
  "SRR1335400"))

setDT(Agerange)
#Plotting
