#=======================================================

#=======================================================

rm(list = ls())
library(survival)
library(dplyr)
library(randomForestSRC)
library(CoxBoost)
library(gbm)
library(survcomp)
library(timeROC)
library(tidyr)
library(tibble)
library(survivalROC)
library(Hmisc)
library(reshape2)
library(data.table)


#=======================================================

#=======================================================


process_data <- function(clinical_path, expr_path, sig_path) {
  clinical <- read.csv(clinical_path, header = TRUE, row.names = 1)
  clinical <- subset(clinical, clinical$OS_Time > 0)
  clinical <- na.omit(clinical)
  clinical$Stage <- NULL
  clinical$Gender <- ifelse(clinical$Gender == "M", 1, 0)
  clinical$Tstage <- as.numeric(factor(clinical$Tstage, levels = c("T1", "T2", "T3", "T4")))
  clinical$Nstage <- ifelse(clinical$Nstage == "N1", 1, 0)
  clinical$Mstage <- ifelse(clinical$Mstage == "M0", 0, 1)
  
  
  expr_or <- read.csv(expr_path, header = TRUE, row.names = 1)
  expr_or <- as.data.frame(t(expr_or))
  
  sig <- read.csv(sig_path, header = TRUE, row.names = 1)$Variable
  expr <- expr_or[, which(colnames(expr_or) %in% sig), drop = FALSE]
  
  
  data_binary <- scale(expr)
  data_binary <- as.data.frame(data_binary)
  colnames(data_binary) <- colnames(expr) 
  rownames(data_binary) <- rownames(expr) 
  
  
  
  expr <- data_binary
  expr$Sample <- rownames(expr)
  
  data <- merge(clinical, expr, by = "Sample")
  rownames(data) <- data$Sample
  data$Sample <- NULL
  return(data)
}


setwd("D:\\S3LUSC")

# TCGA 
tcga_data <- process_data(
  clinical_path = "St1datasets\\TCGA\\clinical.csv",
  expr_path = "St1datasets\\TCGA\\exprSet.csv",
  sig_path = "st3valSubtype\\top9_importance.csv"
)

# GSE73403
gse73403_data <- process_data(
  clinical_path = "St1datasets\\GSE73403\\clinical.csv",
  expr_path = "St1datasets\\GSE73403\\exprSet.csv",
  sig_path = "st3valSubtype\\top9_importance.csv"
)

# GSE30219 
gse30219_data <- process_data(
  clinical_path = "St1datasets\\GSE30219\\clinical.csv",
  expr_path = "St1datasets\\GSE30219\\exprSet.csv",
  sig_path = "st3valSubtype\\top9_importance.csv"
)

str(tcga_data)
str(gse73403_data)
str(gse30219_data)


train_data <- tcga_data
test_data <- gse30219_data

str(train_data)


#=======================================================

#=======================================================

set.seed(123) 


train_index <- sample(seq_len(nrow(train_data)), size = 0.8 * nrow(train_data))

train_split <- train_data[train_index, ]     
validation_split <- train_data[-train_index, ]

#=======================================================

#=======================================================


fit_RF1 <- rfsrc(Surv(OS_Time, OS) ~ ., 
                 data = train_split, keep.data = FALSE, importance = TRUE,
                 ntree = 100, mtry = sqrt(ncol(train_split)))
fit_RF2 <- rfsrc(Surv(OS_Time, OS) ~ ., data = train_split, keep.data = FALSE, ntree = 200, mtry = log(ncol(train_split)))
fit_RF3 <- rfsrc(Surv(OS_Time, OS) ~ ., data = train_split, keep.data = FALSE, ntree = 500, mtry = ncol(train_split) / 2)


fit_CoxBoost1 <- CoxBoost(time = train_split$OS_Time, 
                          status = train_split$OS,
                          x = as.matrix(train_split[, -c(1, 2)]), penalty = 50)
fit_CoxBoost2 <- CoxBoost(time = train_split$OS_Time, 
                          status = train_split$OS, 
                          x = as.matrix(train_split[, -c(1, 2)]), penalty = 100)
fit_CoxBoost3 <- CoxBoost(time = train_split$OS_Time,
                          status = train_split$OS, 
                          x = as.matrix(train_split[, -c(1, 2)]), penalty = 200)


fit_gbm1 <- gbm(Surv(OS_Time, OS) ~ ., data = train_split, 
                distribution = 'coxph', n.trees = 100)

fit_gbm2 <- gbm(Surv(OS_Time, OS) ~ ., data = train_split, 
                distribution = 'coxph', n.trees = 200)

fit_gbm3 <- gbm(Surv(OS_Time, OS) ~ ., data = train_split,
                distribution = 'coxph', n.trees = 500)


#=======================================================

#=======================================================


list_of_dfs <- list()


loopNum <- 10

for (h in 1:loopNum) {
  
  print(h)

  set.seed(1234 + h)
  
  resampled_test_data <- validation_split[sample(1:nrow(validation_split), nrow(validation_split)*0.5), ]
  
  
  resampled_x_test <- model.matrix(Surv(OS_Time, OS) ~ ., data = resampled_test_data)[, -1]
  
  
  risk_scores_rfSRC1 <- as.numeric(predict(fit_RF1, resampled_test_data)$predicted)
  risk_scores_rfSRC2 <- as.numeric(predict(fit_RF2, resampled_test_data)$predicted)
  risk_scores_rfSRC3 <- as.numeric(predict(fit_RF3, resampled_test_data)$predicted)
  
  risk_scores_CoxBoost1 <- as.numeric(predict(fit_CoxBoost1, newdata = as.matrix(resampled_test_data[, -c(1, 2)]), type = "lp"))
  risk_scores_CoxBoost2 <- as.numeric(predict(fit_CoxBoost2, newdata = as.matrix(resampled_test_data[, -c(1, 2)]), type = "lp"))
  risk_scores_CoxBoost3 <- as.numeric(predict(fit_CoxBoost3, newdata = as.matrix(resampled_test_data[, -c(1, 2)]), type = "lp"))
  
  risk_scores_GBM1  <- predict(fit_gbm1, newdata = resampled_test_data, type = "link")
  risk_scores_GBM2  <- predict(fit_gbm2, newdata = resampled_test_data, type = "link")
  risk_scores_GBM3  <- predict(fit_gbm3, newdata = resampled_test_data, type = "link")
  
  scores_df <- data.frame(
    rfSRC1 = risk_scores_rfSRC1,
    rfSRC2 = risk_scores_rfSRC2,
    rfSRC3 = risk_scores_rfSRC3,
    
    CoxB1 = risk_scores_CoxBoost1,
    CoxB2 = risk_scores_CoxBoost2,
    CoxB3 = risk_scores_CoxBoost3,
    
    GBM1 = risk_scores_GBM1,
    GBM2 = risk_scores_GBM2,
    GBM3 = risk_scores_GBM3
  )
  
  head(scores_df)
  
  
  rownames(scores_df) <- rownames(resampled_test_data)
  scores_df$OS <- resampled_test_data$OS
  scores_df$OS_Time <- resampled_test_data$OS_Time
  
  scores_df <- scores_df %>%
    dplyr::select(OS, OS_Time, everything())
  
  #C-index
  cindex_list <- list()
  
  for (col in colnames(scores_df)[3:ncol(scores_df)]) {
    
    surv_obj <- Surv(scores_df$OS_Time, scores_df$OS)
    risk_scores <- scores_df[[col]]
    c_index <- rcorr.cens(risk_scores, surv_obj)
    c_index_value <- as.numeric(c_index["C Index"])
    c_index_value <- ifelse(c_index_value < 0.5, 1 - c_index_value, c_index_value)
    cindex_list[[col]]  <- c_index_value
    
  }
  
  cindex_dt1 <- as.data.frame(cindex_list)
  
  list_of_dfs[[h]] <- cindex_dt1
  
}


list_of_dfs

cindex_dt <- do.call(rbind, list_of_dfs)


mean_cindex <- as.data.frame(colMeans(cindex_dt))

mean_cindex <- t(mean_cindex)

mean_cindex <-  as.data.frame(mean_cindex)

mean_cindex

#=======================================================

#=======================================================


rf_df <- data.frame(
  ntree = c(100, 200, 500),
  cindex = c(mean_cindex[1, "rfSRC1"], mean_cindex[1, "rfSRC2"], mean_cindex[1, "rfSRC3"])
)


coxboost_df <- data.frame(
  penalty = c(50, 100, 200),
  cindex = c(mean_cindex[1, "CoxB1"], mean_cindex[1, "CoxB2"], mean_cindex[1, "CoxB3"])
)


gbm_df <- data.frame(
  n_trees = c(100, 500, 1000),
  cindex = c(mean_cindex[1, "GBM1"], mean_cindex[1, "GBM2"], mean_cindex[1, "GBM3"])
)

best_rf <- rf_df[which.max(rf_df$cindex), ]
best_coxboost <- coxboost_df[which.max(coxboost_df$cindex), ]
best_gbm <- gbm_df[which.max(gbm_df$cindex), ]

cat("Best parameters based on cindex:\n")
cat("Random Forest SRC: ntree =", best_rf$ntree, "with cindex =", best_rf$cindex, "\n")
cat("CoxBoost: penalty =", best_coxboost$penalty, "with cindex =", best_coxboost$cindex, "\n")
cat("GBM: n_trees =", best_gbm$n_trees, "with cindex =", best_gbm$cindex, "\n")
