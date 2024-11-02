#=======================================================

#=======================================================

#load packages
rm(list = ls())
setwd("D:\\St1datasets\\GSE30219")

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
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
gsename <- "GSE30219"
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================

#annotation
ids <- idmap('GPL570', type = 'soft')
ids <- ids[nchar(ids[,2]) > 1, ]
ids1 <- ids[grepl('///', ids[,2]), ]
ids2 <- ids[!grepl('///', ids[,2]), ]

tmp <- do.call(rbind, apply(ids1, 1, function(x) {
  x[1]
  x[2]
  data.frame(ID = x[1], symbol = strsplit(x[2], ' /// ')[[1]])
}))

ids <- rbind(ids2, tmp)
anno <- annoGene(ids$symbol, "SYMBOL")
ids <- merge(ids, anno, by.x = 'symbol', by.y = 'SYMBOL', all.x = TRUE)
ids <- subset(ids, select = c("ID", "symbol"))

#=======================================================

#=======================================================

#expression matrix
exprSet <- as.data.frame(exprs(gse$GSE30219_series_matrix.txt.gz))
str(exprSet)


exprSet$ID <- rownames(exprSet)
express <- merge(x = ids, y = exprSet, by = "ID")
exprSet <- NULL  
express$ID <- NULL


express[which(is.na(express), arr.ind = TRUE)] <- 0 
express[1:4, 1:4]

exprSet <- aggregate(x = express[, 2:ncol(express)],
                     by = list(express$symbol),
                     FUN = max)

exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4, 1:4]

#===================================================

#===================================================

#data normalization
exprSet <- log2(exprSet+1)
min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

#processing clinical data
pd <- pData(gse$GSE30219_series_matrix.txt.gz)
head(pd)
table(pd$"histology:ch1")
pd <- subset(pd, pd$"histology:ch1" == "SQC")
head(pd)

pd <- subset(pd, select=c("geo_accession", 
                          "status:ch1", "follow-up time (months):ch1",
                          "age at surgery:ch1", "gender:ch1", 
                          "pt stage:ch1" , "pn stage:ch1",  "pm stage:ch1"))

head(pd)
names(pd) <- c("Sample", "OS", "OS_Time",  "Age", "Gender", 
               "Tstage", "Nstage", "Mstage")
head(pd)

table(pd$Gender)

table(pd$OS)
pd$OS[pd$OS == "ALIVE"] <- 0
pd$OS[pd$OS == "DEAD"] <- 1
table(pd$OS)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
pd$OS <- as.numeric(as.character(pd$OS))
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))

pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 1)
pd <- subset(pd, pd$OS_Time>0)
str(pd)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]


#=======================================================

#=======================================================

#common samples
samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)

pd <- pd[which(pd$Sample %in% samplesname), ]
exprSet <- exprSet[, which(colnames(exprSet) %in% samplesname)]

head(pd)


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
