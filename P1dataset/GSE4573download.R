#=======================================================

#=======================================================
# 清除所有变量
rm(list = ls())

# 设置工作目录
setwd("D:\\S3LUSC\\P1dataset\\GSE4573")

# 加载包
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

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 16)

gsename <- "GSE4573"

# 下载基因芯片数据
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================

# 获取基因 ID
ids <- idmap('GPL96', type = 'soft')
ids <- ids[nchar(ids[,2]) > 1, ]
ids1 <- ids[grepl('///', ids[,2]), ]
ids2 <- ids[!grepl('///', ids[,2]), ]

# 处理多个基因名的情况
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

# 构建表达矩阵
exprSet <- as.data.frame(exprs(gse$GSE4573_series_matrix.txt.gz))

# 将 ID 转换为基因名
exprSet$ID <- rownames(exprSet)
express <- merge(x = ids, y = exprSet, by = "ID")
exprSet <- NULL  # 释放内存
express$ID <- NULL

# 处理缺失值
express[which(is.na(express), arr.ind = TRUE)] <- 0 
express[1:4,1:4]

# 去除重复基因并获取最大值
exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$symbol),
                     FUN = max)

# 清洗数据
exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]

# ===================================================
# 
# ===================================================


# 
# # 计算每行的变异系数
# cv <- apply(exprSet, 1, function(x) sd(x) / mean(x))
# cv <- as.data.frame(cv)
# names(cv)[1] <- 'value'
# dim(cv)
# cv <- subset(cv, cv$value > 0.001)
# dim(cv)
# exprSet <-  exprSet[which(rownames(exprSet) %in% rownames(cv)),]
# 
# 
# exprSet[1:4,1:4]
# exprSet <- as.data.frame(t(exprSet))
# 
# data_binary <- apply(exprSet, 2, function(x) {
#   median_x <- median(x)  # 计算中位数
#   return(ifelse(x > median_x, 1, 0))  # 大于中位数的赋值为1，否则为0
# })  
# 
# rownames(data_binary) <- rownames(exprSet)
# exprSet <- data_binary
# exprSet[1:4,1:4]
# exprSet <- as.data.frame(t(exprSet))
# exprSet[1:4,1:4]



#=======================================================

#=======================================================

pd <- pData(gse$GSE4573_series_matrix.txt.gz)
pd <- subset(pd, select=c("title", "geo_accession"))
dt <- fread("clinical.txt")
dt <- subset(dt, select=c("RNA array SCC ID", "Date of death", "Survival time (mo)",
                          "AGE", "SEX", "STAGE", "T", "N", "M"))
head(dt)
names(dt) <-  c("title", "OS", "OS_Time", "Age", "Gender", "Stage", "Tstage", "Nstage", "Mstage")
dt <- as.data.frame(dt)
table(dt$OS)
str(dt)



dt$OS[dt$OS==""] <- 0
table(dt$OS)
dt$OS[dt$OS!=0] <- 1
table(dt$OS)

pd <- merge(pd, dt, by="title")
pd$title <- NULL

str(pd)

pd$OS <- as.numeric(as.character(pd$OS))




pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 1)


head(pd)


names(pd)[1] <- "Sample"

table(pd$Stage)
pd$Stage[pd$Stage=="Ia"] <- "I"
pd$Stage[pd$Stage=="Ib"] <- "I"
pd$Stage[pd$Stage=="IIa"] <- "II"
pd$Stage[pd$Stage=="IIb"] <- "II"
pd$Stage[pd$Stage=="IIIa"] <- "III"
pd$Stage[pd$Stage=="IIIb"] <- "III"
table(pd$Stage)


table(pd$Tstage)
pd$Tstage <- paste0("T", pd$Tstage)
table(pd$Tstage)


table(pd$Nstage)
pd$Nstage <- paste0("N", pd$Nstage)
table(pd$Nstage)

table(pd$Mstage)
pd$Mstage <- paste0("M", pd$Mstage)
table(pd$Mstage)

pd$Mstage[pd$Mstage=="MNA"] <- "Not available"


head(pd)

rownames(pd) <- pd$Sample


str(pd)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]

str(pd)

table(pd$OS)

mean(pd$OS_Time)

#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)

pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]


write.csv(pd, file = "clinical.csv")
pd <- subset(pd, select=c("OS", "OS_Time"))
head(pd)

#=======================================================

#=======================================================

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")
