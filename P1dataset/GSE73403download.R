#=======================================================

#=======================================================
# 清除所有变量
rm(list = ls())

# 设置工作目录
setwd("D:\\S3LUSC\\P1dataset\\GSE73403")

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

gsename <- "GSE73403"

# 下载基因芯片数据
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================

# 获取基因 ID
ids <- idmap('GPL6480', type = 'soft')
ids <- na.omit(ids)

#=======================================================

#=======================================================

# 构建表达矩阵
exprSet <- as.data.frame(exprs(gse$GSE73403_series_matrix.txt.gz))

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


table(pd$Gender)
pd$Gender[pd$Gender=="Female"] <- "F"
pd$Gender[pd$Gender=="Male"] <- "M"
table(pd$Gender)


table(pd$Stage)
pd$Stage[pd$Stage=="1"] <- "I"
pd$Stage[pd$Stage=="2"] <- "II"
pd$Stage[pd$Stage=="3"] <- "III"
table(pd$Stage)

head(pd)


pd$Tstage <- gsub(".*(T[0-9]+).*", "\\1", pd$TNM)
pd$Nstage <- gsub(".*(N[0-9]+).*", "\\1", pd$TNM)
pd$Mstage<- gsub(".*(M[0-9]+).*", "\\1", pd$TNM)

table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)

pd$TNM <- NULL

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

