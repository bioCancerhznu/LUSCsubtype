#=======================================================

#=======================================================
rm(list = ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

#=======================================================

#=======================================================

setwd("D:\\st4TME")

df <- read.csv("epic.csv", header = TRUE, row.names = 1)
df <- as.data.frame(t(df)) 
  
print(df[1:5, 1:5])
  
df <- as.data.frame(t(df))
df$Sample <- rownames(df)
head(df)[1:5, 1:5]

df$otherCells <- NULL

#=======================================================

#=======================================================

setwd("D:\\St2subtype")

train <- read.csv("Train_subtype.csv", row.names = 1)

head(train)

data <- merge(train, df, by = "Sample")

rownames(data) <- data$Sample
data$Sample <- NULL

head(data)[1:5, 1:5]

data_long <- melt(data, id.vars = "group", 
                  variable.name = "CellType", value.name = "Value")

head(data_long)

#=======================================================

#=======================================================

p <- ggboxplot(data_long, 
               x = "CellType", 
               y = "Value", 
               color = "group", 
               palette = c("#E8751A", "#5356FF"),
               xlab = "Cell Type", 
               ylab = "Value",
               width = 0.7,           
               outlier.shape = NA) +  
  stat_compare_means(aes(group = group), 
                     label = "p.signif", 
                     method = "t.test",
                     label.y = 0.5) +   
  theme_light() + ylim(-1,4) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), 
    plot.title = element_text(hjust = 0.5)                            
  ) +
  labs(title = "Boxplot of Cell Types by Group")         


print(p)


setwd("D:\\st4TME")
pdf(file = "cells.pdf", height = 5, width = 5)
print(p)
dev.off()


