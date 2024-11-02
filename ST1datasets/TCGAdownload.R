# ===================================================

# ===================================================

#load packages
rm(list = ls()) 
library(dplyr)
library(tidyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(stringr)
library(openxlsx)

setwd('D:\\St1datasets\\TCGA')

#===================================================

#===================================================

cancer <- TCGAbiolinks:::getGDCprojects()$project_id
cancer <- str_subset(cancer, "TCGA")
cancer <- sort(cancer)
cancer_select <- "TCGA-LUSC"


query <- GDCquery(
  project = cancer_select,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 50)

# ===================================================

# ===================================================

expr <- GDCprepare(query = query)
TPM <- as.data.frame(assay(expr, i = "tpm_unstrand"))

#Gene annotation
anno <- as.data.frame(expr@rowRanges@elementMetadata)
anno <- subset(anno, select = c('gene_name', 'gene_id'))
TPM$gene_id <- rownames(TPM)
TPM <- merge(anno, TPM, by = 'gene_id')
TPM$gene_id <- NULL

#Missing values
TPM[which(is.na(TPM), arr.ind = TRUE)] <- 0
TPM[1:5,1:5]

#genes with multiple probes
exprSet <- aggregate(x = TPM[, 2:ncol(TPM)], 
                     by = list(TPM$gene_name), 
                     FUN = max)
exprSet[1:5,1:5]

names(exprSet)[1] <- "ID"
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]
max(exprSet)
min(exprSet)

exprSet <- na.omit(exprSet)
colnames(exprSet) <- chartr(old = '-', new = "_", colnames(exprSet))
exprSet[1:4,1:4]

#===================================================

#===================================================

#data normalization
exprSet <- log2(exprSet+1)
min(exprSet)
max(exprSet)

#===================================================

#===================================================

#remove normal samples
metad <- data.frame(names(exprSet))
colnames(metad)[1] <- 'sample'
metad$id <- substr(metad$sample, start = 1, stop = 12)
metad$tape <- substr(metad$sample, start = 14, stop = 16)
metad <- metad[order(metad$id), ]
metad <- subset(metad, metad$tape == "01A")
metad <- metad %>% distinct(id, .keep_all = TRUE)
exprSet <- exprSet[, which(colnames(exprSet) %in% metad$sample)]
colnames(exprSet) <- substr(colnames(exprSet), start = 1, stop = 12)
exprSet[1:5,1:5]

#===================================================

#===================================================

#downloading survival information
my_data <- read.xlsx("mmc1.xlsx")
my_data <- subset(my_data, my_data$type=="LUSC")

my_data <- my_data %>%
  dplyr::select("bcr_patient_barcode",
                "OS", "OS.time",
                "age_at_initial_pathologic_diagnosis",
                "gender","ajcc_pathologic_tumor_stage" )
names(my_data) <- c("Sample",  "OS", "OS_Time",
                    "Age", "Gender",
                    "Stage")

str(my_data)
my_data$Gender <- ifelse(my_data$Gender == 'FEMALE', "F", "M")
my_data$OS_Time <- my_data$OS_Time / 365
my_data$OS_Time <- round(my_data$OS_Time, 1)
str(my_data)

clin <- GDCquery_clinic("TCGA-LUSC", 
                        type = "clinical", 
                        save.csv = FALSE)

clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_t",
                                "ajcc_pathologic_n",  "ajcc_pathologic_m" ))

names(clin) <- c("Sample",  "Tstage", "Nstage", "Mstage")

pd <- merge(my_data, clin, by = "Sample")

head(pd)

#===================================================

#===================================================

#preprocessing clinical information
table(pd$Stage)
# pd$Stage[pd$Stage=="Stage I"] <- "I"
# pd$Stage[pd$Stage=="Stage IA"] <- "I"
# pd$Stage[pd$Stage=="Stage IB"] <- "I"
# pd$Stage[pd$Stage=="Stage II"] <- "II"
# pd$Stage[pd$Stage=="Stage IIA"] <- "II"
# pd$Stage[pd$Stage=="Stage IIB"] <- "II"
# pd$Stage[pd$Stage=="Stage III"] <- "III"
# pd$Stage[pd$Stage=="Stage IIIA"] <- "III"
# pd$Stage[pd$Stage=="Stage IIIB"] <- "III"
# pd$Stage[pd$Stage=="Stage IV"] <- "IV"
# pd$Stage[pd$Stage=="[Discrepancy]"] <- NA
pd$Stage <- NULL

table(pd$Tstage)
pd$Tstage[pd$Tstage=="T1"] <- "T1"
pd$Tstage[pd$Tstage=="T1a"] <- "T1"
pd$Tstage[pd$Tstage=="T1b"] <- "T1"
pd$Tstage[pd$Tstage=="T2"] <- "T2"
pd$Tstage[pd$Tstage=="T2a"] <- "T2"
pd$Tstage[pd$Tstage=="T2b"] <- "T2"
pd$Tstage[pd$Tstage=="T3"] <- "T3"
pd$Tstage[pd$Tstage=="T4"] <- "T4"
table(pd$Tstage)


table(pd$Nstage)
pd$Nstage[pd$Nstage=="NX"] <- NA
table(pd$Nstage)

table(pd$Mstage)
pd$Mstage[pd$Mstage=="M1a"] <- "M1"
pd$Mstage[pd$Mstage=="M1b"] <- "M1"
pd$Mstage[pd$Mstage=="MX"] <- NA
table(pd$Mstage)


str(pd)

dim(pd)
pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
dim(pd)

pd <- subset(pd, pd$OS_Time > 0)
dim(pd)

head(pd)

pd$Sample <- chartr(old = '-', new = "_", pd$Sample)
rownames(pd) <- pd$Sample

#===================================================

#===================================================

#common samples
samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]


#===================================================

#===================================================

str(pd)

mean(pd$Age, na.rm = TRUE)
table(pd$Gender)
table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)
table(pd$OS)
mean(pd$OS_Time)


#===================================================

#===================================================

#save data 
write.csv(pd, file = "clinical.csv")
pd <- subset(pd, select=c("OS", "OS_Time"))
head(pd)

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

