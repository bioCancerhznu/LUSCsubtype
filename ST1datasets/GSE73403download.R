#=======================================================

#=======================================================

#load packages
rm(list = ls())

setwd("D:\\St1datasets\\GSE73403")

library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(ggplot2)
library(AnnoProbe)

#=======================================================

#=======================================================

#setting dataset name
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 16)
gsename <- "GSE73403"
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================

#annotation
ids <- idmap('GPL6480', type = 'soft')
ids <- na.omit(ids)

#=======================================================

#=======================================================

#expression matrix
exprSet <- as.data.frame(exprs(gse$GSE73403_series_matrix.txt.gz))

exprSet$ID <- rownames(exprSet)
express <- merge(x = ids, y = exprSet, by = "ID")
exprSet <- NULL  
express$ID <- NULL

express[which(is.na(express), arr.ind = TRUE)] <- 0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$symbol),
                     FUN = max)


exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]

#===================================================

#===================================================

#data normalization
exprSet <- log2(exprSet+1)
min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

#processing clinical data
pd <- pData(gse$GSE73403_series_matrix.txt.gz)
pd <- subset(pd, select=c("geo_accession", 
                          "events (death = 1, alive = 0):ch1",
                          "survival (year):ch1", 
                          "age:ch1", "gender:ch1",
                           "Stage:ch1","tnm:ch1"))
head(pd)
names(pd) <- c("Sample", "OS", "OS_Time", "Age", "Gender", 
               "Stage", "TNM")
table(pd$OS)
str(pd)

pd$OS <- as.numeric(as.character(pd$OS))
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- round(pd$OS_Time, 1)
str(pd)
pd <- subset(pd, pd$OS_Time>0)


pd$Age <- as.numeric(as.character(pd$Age))

table(pd$Gender)
pd$Gender[pd$Gender=="Female"] <- "F"
pd$Gender[pd$Gender=="Male"] <- "M"
table(pd$Gender)


table(pd$Stage)
pd$Stage <- NULL


head(pd)

pd$Tstage <- gsub(".*(T[0-9]+).*", "\\1", pd$TNM)
pd$Nstage <- gsub(".*(N[0-9]+).*", "\\1", pd$TNM)
pd$Mstage<- gsub(".*(M[0-9]+).*", "\\1", pd$TNM)

table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)

pd$TNM <- NULL

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]

#=======================================================

#=======================================================

#common samples
samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)

pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]


#=======================================================

#=======================================================


str(pd)

mean(pd$Age, na.rm = TRUE)
table(pd$Gender)
table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)
table(pd$OS)
mean(pd$OS_Time)


#=======================================================

#=======================================================


write.csv(pd, file = "clinical.csv")
pd <- subset(pd, select=c("OS", "OS_Time"))
head(pd)

#saving data files
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")
