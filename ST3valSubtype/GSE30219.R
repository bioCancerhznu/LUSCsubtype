#=======================================================

#=======================================================

rm(list = ls())
library(tibble)
library(survival)
library(survminer)

#=======================================================

#=======================================================

setwd("D:\\St2subtype")
seleGnes <- read.csv(file = "data_train.csv", row.names = 1)

setwd("D:St1datasets\\GSE30219")
testdata <- read.csv("exprSet.csv", header = T, row.names = 1)
testdata <- testdata[which(rownames(testdata) %in% colnames(seleGnes)),]
testdata <- as.data.frame(t(testdata))


preprocess_data <- apply(testdata, 2, function(x) {
  q20 <- quantile(x, 0.2)
  q80 <- quantile(x, 0.8)
  x <- ifelse(x <= q20, 0, x)     
  x <- ifelse(x >= q80, 1, x)       
  x <- ifelse(x > q20 & x < q80, (x - q20) / (q80 - q20), x) 
  
  return(x)
})

data_binary <- as.data.frame(preprocess_data)
colnames(data_binary) <- colnames(testdata)
rownames(data_binary) <- rownames(testdata)
data_binary[1:5,1:5]
dim(data_binary)

testdata <- data_binary


setwd("D:\\st3valSubtype")
load("rf_model.Rdata")

#=======================================================

#=======================================================

group  <- predict(rf_model, newdata = testdata)

subtype <- as.data.frame(group)

subtype$Sample <- rownames(data_binary)

head(subtype)

setwd('D:St1datasets\\GSE30219')

sur_data <- read.csv("pd.csv", header = T, row.names = 1)
sur_data <- subset(sur_data, sur_data$OS_Time > 0)
sur_data$Sample <- rownames(sur_data)
data <- merge(subtype, sur_data, by="Sample")
head(data)
data <- na.omit(data)


str(data)

head(data)

table(data$group)

res.cox <- coxph(Surv(OS_Time, OS) ~ group, data = data)



#=======================================================

#=======================================================


weight_C1 <- nrow(data[data$group == "C2", ]) / nrow(data[data$group == "C1", ])
weight_C2 <- 1 


data$weight <- ifelse(data$group == "C1", weight_C1, weight_C2)
cox_model <- coxph(Surv(OS_Time, OS) ~ group, data = data, weights = weight)
summary(cox_model)
p_value <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
p_value

#=======================================================

#=======================================================



fit <- survfit(Surv(OS_Time, OS) ~ group, data = data)

ggsurv <- ggsurvplot(fit, data = data, pval = TRUE, conf.int = FALSE, 
                     legend.title = "Subtype", xlim = c(0, 5), 
                     xlab = "Time in years", break.time.by = 1, 
                     pval.size = 8, risk.table = "absolute", 
                     risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, 
                     risk.table.fontsize = 5, 
                     font.y = c(16, "bold"), font.x = c(16, "bold"), 
                     legend = "top", font.legend = c(16, "bold"),
                     palette = c("#E8751A", "#5356FF"))

ggsurv

setwd("D:st3valSubtype")


pdf(file = "survGSE30219.pdf", height = 5, width = 5)
print(ggsurv)
dev.off()

save(p_value, file = "p_value.Rdata")

