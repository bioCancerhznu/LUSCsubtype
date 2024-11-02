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
library(pheatmap)
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


setwd("D:\\DuanWork\\S3LUSC")

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
  
  set.seed(1234 + h)
  
  resampled_test_data <- test_data[sample(1:nrow(test_data), nrow(test_data)*0.8), ]
  
  
  resampled_x_test <- model.matrix(Surv(OS_Time, OS) ~ ., data = resampled_test_data)[, -1]
  
  
  risk_scores_rfSRC <- as.numeric(predict(fit_RF, resampled_test_data)$predicted)
  risk_scores_CoxBoost <- as.numeric(predict(fit_CoxBoost, newdata = as.matrix(resampled_test_data[, -c(1, 2)]), type = "lp"))
  risk_scores_GBM  <- predict(fit_gbm, newdata = resampled_test_data, type = "link")
  
  scores_df <- data.frame(
    rfSRC = risk_scores_rfSRC,
    CoxB = risk_scores_CoxBoost,
    GBM = risk_scores_GBM
  )
  
  scores_df$ensemble <- scores_df$rfSRC + scores_df$CoxB + scores_df$GBM
  
  head(scores_df)
  
  
  rownames(scores_df) <- rownames(resampled_test_data)
  scores_df$OS <- resampled_test_data$OS
  scores_df$OS_Time <- resampled_test_data$OS_Time
  
  scores_df <- scores_df %>%
    dplyr::select(OS, OS_Time, everything())
  
  cindex_list <- list()
  
  for (col in colnames(scores_df)[3:ncol(scores_df)]) {
    
    risk_scores <- scores_df[[col]]
    cindex <- concordance.index(risk_scores, surv.time = scores_df$OS_Time, 
                                surv.event = scores_df$OS)$c.index
    cindex_list[[col]] <- cindex
  }
  
  cindex_dt1 <- as.data.frame(cindex_list)
  
  list_of_dfs[[h]] <- cindex_dt1
  
}


list_of_dfs

cindex_dt <- do.call(rbind, list_of_dfs)

cindex_dt$iter <- 1:loopNum

#=======================================================

#=======================================================

cindex_long <- melt(cindex_dt, id.vars = "iter", 
                    variable.name = "Algorithm", value.name = "cindex")

head(cindex_long)

cindex_long$iter <- NULL

algorithm_order <- cindex_long %>%
  group_by(Algorithm) %>%
  summarise(mean_cindex = mean(cindex, na.rm = TRUE)) %>%
  arrange(mean_cindex) %>%
  pull(Algorithm)


cindex_long$Algorithm <- factor(cindex_long$Algorithm, levels = algorithm_order)


p <- ggplot(cindex_long, aes(x = Algorithm, y = cindex)) +
  geom_boxplot(width=0.5, color='#C4D7FF', fill="#C4D7FF", outlier.shape = NA) +
  stat_boxplot(geom="errorbar", width=0.5, color='#C4D7FF') + 
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 3)), 
               vjust = 0.1, color = "#CB6040", size = 3.5) +
  labs(title = "Cindex Boxplot by Algorithm", 
       x = "Algorithm", 
       y = "Cindex Value") +
  coord_flip() +
  theme_bw()

print(p)

getwd()

setwd("D:\\DuanWork\\S3LUSC\\st8survML\\GSE30219")
pdf(file = "cindex.pdf", height = 5, width = 4)
print(p)
dev.off()



