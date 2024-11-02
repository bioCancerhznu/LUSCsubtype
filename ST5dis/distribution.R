#=======================================================
# 
#=======================================================

rm(list = ls())
library(dplyr)
library(TCGAbiolinks)
library(officer)
library(flextable)
library(compareGroups)

#=======================================================
# 
#=======================================================

setwd("D:\\St1datasets\\TCGA")

query <- GDCquery(project = "TCGA-LUSC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")

GDCdownload(query)

clinical_drug_data <- GDCprepare(query)
head(clinical_drug_data)
names(clinical_drug_data)

#=======================================================
# 
#=======================================================

dt1 <- clinical_drug_data$clinical_drug_lusc
head(dt1)
dt1 <- dt1[-c(1:2),]
dt1 <- as.data.frame(dt1)
dt1 <- subset(dt1, select = c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name",
                              "pharmaceutical_therapy_type"))


names(dt1) <- c("Sample", "Name", "Chemotherapy")
dt1$Sample <- chartr(old = '-', new = '_', x=dt1$Sample)

dt1 <- dt1 %>%
  distinct(Sample, .keep_all = TRUE)

dt1$Name <- NULL

#=======================================================
# 
#=======================================================


setwd('D:\\St1datasets\\TCGA')
sur_data <- read.csv("clinical.csv", header = T, row.names = 1)
sur_data$Sample <- rownames(sur_data)

setwd("D:\\St2subtype")
subtype <-  read.csv("Train_subtype.csv", header = T, row.names = 1)

data <- merge(subtype, sur_data, by="Sample")

#=======================================================
# 
#=======================================================

result <- dplyr::left_join(data, dt1, by = "Sample")

#=======================================================
# 
#=======================================================

dt5 <- clinical_drug_data$clinical_radiation_lusc
head(dt5)
dt5 <- dt5[-c(1:2),]
dt5 <- as.data.frame(dt5)
names(dt5)
dt5 <- subset(dt5, select = c("bcr_patient_barcode", "radiation_therapy_type"))


names(dt5) <- c("Sample",  "Radiation")
dt5$Sample <- chartr(old = '-', new = '_', x=dt5$Sample)
dt5$Radiation <- "Yes"

dt5 <- dt5 %>%
  distinct(Sample, .keep_all = TRUE)

#=======================================================
# 
#=======================================================


result <- dplyr::left_join(result, dt5, by = "Sample")

head(result)

#=======================================================
# 
#=======================================================

table(result$Chemotherapy, useNA = "ifany")

result <- result %>%
  mutate(Chemotherapy = ifelse(Chemotherapy == "Chemotherapy", "Yes", "No")) %>%
  mutate(Chemotherapy = ifelse(is.na(Chemotherapy), "No", Chemotherapy))

table(result$Chemotherapy, useNA = "ifany")


#=======================================================
# 
#=======================================================

head(result)

table(result$Radiation, useNA = "ifany")

result <- result %>%
  mutate(Radiation = ifelse(Radiation == "Yes", "Yes", "No")) %>%
  mutate(Radiation = ifelse(is.na(Radiation), "No", Radiation))

table(result$Radiation, useNA = "ifany")

#=======================================================
# 
#=======================================================

table(result$OS, useNA = "ifany")
result <- result %>%
  mutate(OS = ifelse(OS == 1, "Dead", "Alive"))

table(result$OS, useNA = "ifany")

#=======================================================
# 
#=======================================================

head(result)


