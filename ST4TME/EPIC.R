#=======================================================

#=======================================================

rm(list = ls())

library(EPIC)

#=======================================================

#=======================================================

setwd("D:\\St1datasets\\TCGA")

expression_data <- read.csv("exprSet.csv", header = TRUE, row.names = 1)

epic_result <- EPIC(bulk = expression_data)

cell_proportions <- epic_result$cellFractions

cell_proportions <- as.data.frame(cell_proportions)

normalized_result <- as.data.frame(scale(cell_proportions))

rownames(normalized_result) <- rownames(cell_proportions)

#=======================================================

#=======================================================

setwd("D:S3LUSC\\st4TME")

write.csv(normalized_result, file = "epic.csv")
