#=======================================================

#=======================================================

library(survival)
library(dplyr)
library(randomForestSRC)
library(survivalsvm)
library(CoxBoost)
library(gbm)
library(survcomp)
library(pec)
library(timeROC)
library(tidyr)
library(pheatmap)
library(tibble)
library(ggplot2)
library(survivalROC)
rm(list = ls())

#=======================================================

#=======================================================


setwd("D:St1datasets\\GSE73403")
expr_or <- read.csv("exprSet.csv", header = T, row.names = 1)
expr_or <- as.data.frame(t(expr_or))

#=======================================================

#=======================================================

repNum <- 10
num_signatures <- 15

signature_names <- paste0("GeneSg", 1:num_signatures)


CindexValue <- data.frame(matrix(0, nrow = repNum, ncol = num_signatures))
colnames(CindexValue) <- signature_names
rownames(CindexValue) <- paste0("Traindata", 1:repNum)


aucValue1 <- data.frame(matrix(0, nrow = repNum, ncol = num_signatures))
colnames(aucValue1) <- signature_names
rownames(aucValue1) <- paste0("Traindata", 1:repNum)


aucValue3 <- data.frame(matrix(0, nrow = repNum, ncol = num_signatures))
colnames(aucValue3) <- signature_names
rownames(aucValue3) <- paste0("Traindata", 1:repNum)

aucValue5 <- data.frame(matrix(0, nrow = repNum, ncol = num_signatures))
colnames(aucValue5) <- signature_names
rownames(aucValue5) <- paste0("Traindata", 1:repNum)

#=======================================================

#=======================================================

setwd("D:\\st3valSubtype")
signature15 <- read.csv("top9_importance.csv", header = T, row.names = 1)
signature15 <-  signature15$Variable

signature1 <- c("MMP12", "PLAU", "EDNRB", "APLN", "AGTR2", "NR4A1")
signature2 <- c("MGST3", "TMED3", "PPIB", "GEMIN6")
signature3 <- c('POLD4', 'MRPL40', 'ITPA', 'ERCC3', 'TK2', 
                'POLR3GL', 'VPS28', 'CANT1', 'SDCBP', 'CCNO')
signature4 <- c('IL6', 'NOD1', 'CASP4')
signature5 <- c("ANGPTL5", "CD1B", "CD1E", "CNTFR", "CTSG", "EDN3", "IL12B", "IL2")
signature6 <- c("CCL15", "CXCL7", "VAV2")
signature7 <- c("CD1B", "CHRNA6", "CLEC12B", "CLEC17A", "CLNK", "INHA", "SLC14A2")
signature8 <- c("FGG", "C3", "FGA", "JUN", "CST3", "CPSF4", "HIST1H2BH")
signature9 <- c("CALCB", "GCGR", "HTR3A", "AMH", "VGF", "SEMA3B", "NRTN", "ENG", "ACVRL1", "NR4A1")
signature10 <- c("CXCL2", "SMAD7", "HELLS", "IL1B")
signature11 <- c("SREBF2", "GP2", "BMX", "NR1H4", "DDX41", "GOPC")
signature12 <- c("S100P", "PLAU", "NOD1", "TRAV39")
signature13 <- c("ABCC5", "CLDN1", "CSTA")
signature14 <- c("ALOX5", "DPP4", "PHKG2", "FADS2", "NOX1")

geneslist <- list(gsg1  = signature1,  gsg2  = signature2,  gsg3  = signature3,  
                  gsg4 =  signature4,  gsg5 =  signature5,
                  gsg6  = signature6,  gsg7  = signature7,  gsg8  = signature8,  
                  gsg9 =  signature9,  gsg10 = signature10,
                  gsg11 = signature11, gsg12 = signature12, gsg13 = signature13, 
                  gsg14 = signature14, gsg15 = signature15)

#=======================================================

#=======================================================

for (g in 1:num_signatures) {
  
  genesnam <- geneslist[[g]]
  print(g)
  
  expr <- expr_or[, which(colnames(expr_or) %in% genesnam)]
  expr$Sample <- rownames(expr)
  
  clincadata <- read.csv("D:St1datasets\\GSE73403\\pd.csv", row.names = 1)
  clincadata$Sample <- rownames(clincadata)
  head(clincadata)
  
  selected_data <- merge(clincadata, expr, by = "Sample")
  rownames(selected_data) <- selected_data$Sample
  selected_data$Sample <- NULL
  
  head(selected_data)
  
  results_list <- list()
  
  iterations <- repNum
  
  for (i in 1:iterations) {
    
    set.seed(1234 + i)
    selected_samples <- selected_data[sample(1:nrow(selected_data), nrow(selected_data)/2), ]
    results_list[[i]] <- selected_samples
    
  }
  
  
  cindex_list <- numeric(length(results_list))
  
  
  for (i in 1:length(results_list)) {
    
    
    data <- results_list[[i]]
    signature_vars <- data[, !(colnames(data) %in% c("OS", "OS_Time"))]
    surv_obj <- Surv(time = data$OS_Time, event = data$OS)
    cox_model <- coxph(surv_obj ~ ., data = signature_vars)
    risk_score <- predict(cox_model, type = "lp")
    cindex <- concordance.index(risk_score, surv.time = data$OS_Time, 
                                surv.event = data$OS)$c.index
    
    cindex_list[i] <- cindex
  }
  
  CindexValue[,g] <- cindex_list
  
  # =========================================================
  
  timeauc_list_1year <- numeric(length(results_list))
  
  
  one_year <- 1
  
  for (i in 1:length(results_list)) {
    
    data <- results_list[[i]]
    signature_vars <- data[, !(colnames(data) %in% c("OS", "OS_Time"))]
    surv_obj <- Surv(time = data$OS_Time, event = data$OS)
    cox_model <- coxph(surv_obj ~ ., data = signature_vars)
    risk_score <- predict(cox_model, type = "lp")
    survival_roc <- survivalROC(Stime = data$OS_Time, status = data$OS,
                                marker = risk_score, predict.time = one_year, method = "KM")
    
    timeauc_list_1year[i] <- survival_roc$AUC
    
  }
  
  aucValue1[,g] <- timeauc_list_1year
  
  # =========================================================
  
  timeauc_list_3years <- numeric(length(results_list))
  
  three_years <- 3
  
  for (i in 1:length(results_list)) {
    
    data <- results_list[[i]]
    signature_vars <- data[, !(colnames(data) %in% c("OS", "OS_Time"))]
    surv_obj <- Surv(time = data$OS_Time, event = data$OS)
    cox_model <- coxph(surv_obj ~ ., data = signature_vars)
    risk_score <- predict(cox_model, type = "lp")
    survival_roc <- survivalROC(Stime = data$OS_Time, status = data$OS,
                                marker = risk_score, predict.time = three_years, method = "KM")
    
    timeauc_list_3years[i] <- survival_roc$AUC
    
  }
  
  aucValue3[,g] <- timeauc_list_3years
  
  # =========================================================
  
  timeauc_list_5years <- numeric(length(results_list))
  
  five_years <- 5
  
  for (i in 1:length(results_list)) {
    
    data <- results_list[[i]]
    signature_vars <- data[, !(colnames(data) %in% c("OS", "OS_Time"))]
    surv_obj <- Surv(time = data$OS_Time, event = data$OS)
    cox_model <- coxph(surv_obj ~ ., data = signature_vars)
    risk_score <- predict(cox_model, type = "lp")
    
    survival_roc <- survivalROC(Stime = data$OS_Time, status = data$OS,
                                marker = risk_score, predict.time = five_years, method = "KM")
    
    timeauc_list_5years[i] <- survival_roc$AUC
    
  }
  
  aucValue5[,g] <- timeauc_list_5years
  
}

#=======================================================

#=======================================================

plot_boxplot <- function(data, title, y_label) {
  data_long <- pivot_longer(data, cols = -dataset, 
                            names_to = "Signature", values_to = "Value")
  data_long <- as.data.frame(data_long)
  data_long <- data_long %>%
    dplyr::arrange(Signature)
  
  data_long$Signature <- factor(data_long$Signature, 
                                levels = paste0("GeneSg", 1:15))
  
  
  p <- ggplot(data_long, aes(x = Signature, y = Value)) +
    geom_boxplot(outlier.shape = NA, fill = "white", color = "#686D76") + 
    stat_summary(fun = mean, geom = "point", shape = 20, size = 4, color = "#FDAF7B") +
    stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 3)), 
                 vjust = -2, color = "#627254", size = 4) +
    labs(x = "", y = y_label, title = title) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "#1E201E", 
                                     size = 12),
          axis.text.y = element_text(color = "#1E201E", size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, hjust = 0.5),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.border = element_rect(color = "#1E201E", fill = NA, size = 0.5)) + 
    coord_cartesian(ylim = c(min(data_long$Value) - 0.1, max(data_long$Value) + 0.1))
  
  print(p)
}


#=======================================================

#=======================================================

CindexValue$dataset <- rownames(CindexValue)
aucValue1$dataset <- rownames(aucValue1)
aucValue3$dataset <- rownames(aucValue3)
aucValue5$dataset <- rownames(aucValue5)


setwd("D:\\st6compare\\GSE73403")


p1 <- plot_boxplot(CindexValue, title = "C-index Values by Signature", y_label = "C-index Value")
ggsave("CindexValue_Boxplot.pdf", plot = p1, width = 6, height = 5)


p2 <- plot_boxplot(aucValue1, title = "1 Year AUC Values by Signature", y_label = "AUC Value")
ggsave("AUC1Year_Boxplot.pdf", plot = p2, width = 6, height = 5)


p3 <- plot_boxplot(aucValue5, title = "5 Year AUC Values by Signature", y_label = "AUC Value")
ggsave("AUC5Year_Boxplot.pdf", plot = p3, width = 6, height = 5)


getwd()

