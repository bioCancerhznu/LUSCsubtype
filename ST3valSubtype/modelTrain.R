#=======================================================

#=======================================================

rm(list = ls())
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(survival)
library(survminer)
library(caret)
library(randomForest)

#=======================================================

#=======================================================

setwd("D:St2subtype")
expr <- read.csv("data_train.csv", header = T, row.names = 1)
expr$Sample <- rownames(expr)
groupdata <- read.csv("Train_subtype.csv", row.names = 1)

data <- merge(groupdata, expr, by = "Sample")

rownames(data) <- data$Sample
data$Sample <- NULL

setwd("D:\\st3valSubtype")
write.csv(data, file = "dataTain.csv")

#=======================================================

#=======================================================

set.seed(123)

data$group <- as.factor(data$group)

train_control <- trainControl(method = "cv", number = 5) 

rf_model <- train(
  group ~ .,       
  data = data,   
  method = "rf",    
  trControl = train_control 
)


print(rf_model)

save(rf_model, file = "rf_model.Rdata")

#=======================================================

#=======================================================

# Calculate Variable Importance and Save Top 9
importance <- varImp(rf_model, scale = FALSE)
print(importance)

importance_df <- as.data.frame(importance$importance)
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]

top9_importance <- head(importance_df[, c("Variable", "Overall")], 9)
print(top9_importance)
write.csv(top9_importance, file = "top9_importance.csv")

top9_genes <- top9_importance$Variable  

top9_genes


#=======================================================

#=======================================================

top9_importance

plot <- ggplot(top9_importance, aes(x = Overall, y = reorder(Variable, Overall))) +
  geom_bar(stat = "identity", fill = "#FFDBB5") +
  xlim(0, 4) + 
  geom_text(aes(label = round(Overall, 3)), hjust = -0.2) +
  labs(x = "Importance Values", y = "Variables") +
  theme_bw()


ggsave("top9_importance_plot.pdf",
       plot = plot, width = 5, height = 5)

#=======================================================

#=======================================================

top9_genes

rf_model <- randomForest(group ~ ., data = data, ntree = 500) 
head(data)[1:6,1:6]


top9_genes[1]

pdf_width <- 5
pdf_height <- 5

setwd("D:st3valSubtype\\PDP")

for (i in 1:9) {

  pdf_filename <- paste0("pdp", top9_genes[i], ".pdf")
  pdf(pdf_filename, width = pdf_width, height = pdf_height)
  partialPlot(rf_model, pred.data = data, 
              x.var = top9_genes[i], 
              which.class = "C1",
              main = paste(top9_genes[i], "- Class C1"))
  dev.off()
  
}
