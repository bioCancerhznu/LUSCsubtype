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
library(ggplot2)
library(survivalROC)
library(Hmisc)


#=======================================================

#=======================================================


process_data <- function(clinical_path, expr_path, sig_path) {
  clinical <- read.csv(clinical_path, header = TRUE, row.names = 1)
  clinical <- subset(clinical, clinical$OS_Time > 0)
  clinical$Stage <- NULL
  clinical$Gender <- ifelse(clinical$Gender == "M", 1, 0)
  clinical$Tstage <- as.numeric(factor(clinical$Tstage, levels = c("T1", "T2", "T3", "T4")))
  clinical$Nstage <- ifelse(clinical$Nstage == "N1", 1, 0)
  clinical$Mstage <- ifelse(clinical$Mstage == "M0", 0, 1)
  
  
  get_mode <- function(x) {
    unique_x <- unique(na.omit(x))
    unique_x[which.max(tabulate(match(x, unique_x)))]
  }
  
  
  clinical <- clinical %>% 
    mutate(across(everything(), ~ ifelse(is.na(.), get_mode(.), .)))
  
  
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


fit_RF <- rfsrc(Surv(OS_Time, OS) ~ ., data = train_data, keep.data = FALSE, ntree = 500)

fit_CoxBoost <- CoxBoost(time = train_data$OS_Time,
                         status = train_data$OS, 
                         penalty = 50,
                         x = as.matrix(train_data[, -c(1, 2)]))

fit_gbm <- gbm(Surv(OS_Time, OS) ~ ., data = train_data, distribution = 'coxph',
               n.trees = 100)


#=======================================================

#=======================================================


list_of_dfs <- list()

loopNum <- 10

for (h in 1:loopNum) {
  
  print(h)
  
  set.seed(123 + h)

  resampled_test_data <- test_data[sample(1:nrow(test_data), nrow(test_data)*0.5), ]
  
  
  resampled_x_test <- model.matrix(Surv(OS_Time, OS) ~ ., data = resampled_test_data)[, -1]
  
  
  risk_scores_rfSRC <- as.numeric(predict(fit_RF, resampled_test_data)$predicted)
  risk_scores_CoxBoost <- as.numeric(predict(fit_CoxBoost, newdata = as.matrix(resampled_test_data[, -c(1, 2)]), type = "lp"))
  risk_scores_GBM  <- predict(fit_gbm, newdata = resampled_test_data, type = "link")
  
  scores_df <- data.frame(
    rfSRC = risk_scores_rfSRC,
    CoxB = risk_scores_CoxBoost,
    GBM = risk_scores_GBM
  )
  
  scores_df$ensemble <- scores_df$rfSRC +  scores_df$CoxB + scores_df$GBM

  head(scores_df)
  
  
  rownames(scores_df) <- rownames(resampled_test_data)
  scores_df$OS <- resampled_test_data$OS
  scores_df$OS_Time <- resampled_test_data$OS_Time
  
  scores_df <- scores_df %>%
    dplyr::select(OS, OS_Time, everything())
  
  auc_list <- list()
  
  for (col in colnames(scores_df)[3:ncol(scores_df)]) {
    
    risk_scores <- scores_df[[col]]
    
    time_dependent_auc <- timeROC(
      T = scores_df$OS_Time,     
      delta = scores_df$OS,      
      marker = risk_scores,  
      cause = 1,           
      times = seq(1, max(scores_df$OS_Time), by = 1) 
    )
    
    
    auc_val <- mean(time_dependent_auc$AUC, na.rm = TRUE)
    
    auc_list[[col]] <- auc_val
  }
  
  auc_dt1 <- as.data.frame(auc_list)
  
  list_of_dfs[[h]] <- auc_dt1
  
}


#=======================================================

#=======================================================


list_of_dfs

auc_dt <- do.call(rbind, list_of_dfs)

auc_dt



auc_dt <- pmin(auc_dt, 1)
auc_dt <- pmax(auc_dt, 0.5)
auc_dt <- auc_dt %>% 
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

auc_dt$iter <- 1:loopNum

#=======================================================

#=======================================================

auc_long <- melt(auc_dt, id.vars = "iter", 
                 variable.name = "Algorithm", value.name = "AUC")

head(auc_long)

auc_long$iter <- NULL


algorithm_order <- auc_long %>%
  group_by(Algorithm) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE)) %>%
  arrange(mean_auc) %>%
  pull(Algorithm)


auc_long$Algorithm <- factor(auc_long$Algorithm, levels = algorithm_order)


#=======================================================

#=======================================================

p <- ggplot(auc_long, aes(x = Algorithm, y = AUC)) +
  geom_boxplot(width=0.5, color='#C4D7FF', fill="#C4D7FF", outlier.shape = NA) +
  stat_boxplot(geom="errorbar", width=0.5, color='#C4D7FF') + 
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 3)), 
               vjust = 0.1, color = "#CB6040", size = 3.5) +
  labs(title = "AUC Boxplot by Algorithm", 
       x = "Algorithm", 
       y = "AUC Value") +
  coord_flip() +
  theme_bw()

print(p)

getwd()

setwd("D:\\S3LUSC\\st8survML\\GSE30219")
pdf(file = "AUC.pdf", height = 5, width = 5)
print(p)
dev.off()



