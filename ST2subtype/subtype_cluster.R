#=======================================================

#=======================================================

rm(list = ls())

library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(factoextra)
library(cluster)
library(survival)
library(survminer)
library(NbClust)

#=======================================================

#=======================================================

setwd('D:St1datasets\\TCGA')
sur_data <- read.csv("pd.csv", header = T, row.names = 1)
sur_data$Sample <- rownames(sur_data)

setwd("D:St2subtype")
data_sub <- read.csv("data_sub.csv", header = T, row.names = 1)

min(data_sub)
max(data_sub)

write.csv(data_sub, file = "data_train.csv")

#=======================================================

#=======================================================

data_sub <- as.matrix(data_sub)

#NbClust
nb_results <- NbClust(data_sub, min.nc=2, max.nc=10,
                      method="kmeans")


barplot(table(nb_results$Best.n[1, ]), 
        xlab="Number of Clusters", 
        ylab="Number of Criteria",
        main="Optimal Number of Clusters by NbClust")


#
best_nc <- nb_results$Best.nc

# n = 2 
n2_indices <- which(best_nc[1, ] == 2)
n2_algorithms_names <- colnames(best_nc)[n2_indices]
cat("The following algorithms determined n=2 as the optimal number of clusters:\n")
print(n2_algorithms_names)


#CH, Silhouette, McClain
ch_values <- nb_results$All.index[,"CH"]
silhouette_values <- nb_results$All.index[,"Silhouette"]
mcclain_values <- nb_results$All.index[,"McClain"]


df <- data.frame(
  Clusters = 2:10,
  CH = ch_values,
  Silhouette = silhouette_values,
  McClain = mcclain_values
)

#ggplot2
df_long <- reshape2::melt(df, id.vars = "Clusters", 
                          variable.name = "Metric",
                          value.name = "Value")


# CH
p2 <- ggplot(df, aes(x = Clusters, y = CH)) +
  geom_line(color = "#424769", size = 1.2) +
  geom_point(color = "#B19470", size = 3) +
  ggtitle("CH Index for Different Clusters") +
  theme_light() +  
  theme(plot.title = element_text(hjust = 0.5))

# Silhouette
p3 <- ggplot(df, aes(x = Clusters, y = Silhouette)) +
  geom_line(color = "#424769", size = 1.2) +
  geom_point(color = "#B19470", size = 3) +
  ggtitle("Silhouette Index for Different Clusters") +
  theme_light() +  
  theme(plot.title = element_text(hjust = 0.5))

# McClain
p4 <- ggplot(df, aes(x = Clusters, y = McClain)) +
  geom_line(color = "#424769", size = 1.2) +
  geom_point(color = "#B19470", size = 3) +
  ggtitle("McClain Index for Different Clusters") +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
print(p3)
print(p4)

pdf(file = "CH.pdf", height = 5, width = 5)
print(p2)
dev.off()


pdf(file = "Silhouette.pdf", height = 5, width = 5)
print(p3)
dev.off()


pdf(file = "McClain.pdf", height = 5, width = 5)
print(p4)
dev.off()


#=======================================================

#=======================================================

# Initialize a dataframe to store k and the corresponding average silhouette score
silhouette_data <- data.frame(k = integer(), avg_silhouette = numeric())


cluster_list <- list()

k_range = 2:10


for (k in k_range) {
  
  kmeans_result <- kmeans(data_sub, centers = k, nstart = 10)
  cluster_list[[as.character(k)]] <- kmeans_result$cluster
  sil_widths <- silhouette(kmeans_result$cluster, dist(data_sub))
  avg_sil <- mean(sil_widths[, "sil_width"])
  silhouette_data <- rbind(silhouette_data, data.frame(k = k, avg_silhouette = avg_sil))
}




cluster_assignments <- do.call(cbind, cluster_list)
row.names(cluster_assignments) <- row.names(data_sub)
colnames(cluster_assignments) <- paste0("sub", colnames(cluster_assignments))


cluster_assignments_df <- as.data.frame(cluster_assignments)
cluster_assignments_df$Sample <- row.names(cluster_assignments_df)

subtype <- cluster_assignments_df
subtype <- subset(subtype, select = c("sub2"))
names(subtype)[1] <- 'group'
subtype$Sample <- rownames(subtype)
table(subtype$group)
subtype$group <- paste0("C", subtype$group)
table(subtype$group)

write.csv(subtype, file = "Train_subtype.csv")

#=======================================================

#=======================================================

data <- merge(subtype, sur_data, by="Sample")
head(data)
data <- na.omit(data)

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

save(p_value, file = "p_value.Rdata")


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


pdf(file = "OSsubtype_Train.pdf", height = 6, width = 6)
ggsurv 
dev.off()

