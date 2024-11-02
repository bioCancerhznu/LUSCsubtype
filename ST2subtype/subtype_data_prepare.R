#=======================================================

#=======================================================

rm(list = ls())

library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(cluster)
library(survival)
library(survminer)
library(biomaRt)

#=======================================================

#=======================================================

setwd('D:\\St1datasets\\TCGA')

sur_data <- read.csv("pd.csv", header = T, row.names = 1)
sur_data$Sample <- rownames(sur_data)


#finding common genes
TCGA <- read.csv("exprSet.csv", header = T, row.names = 1)
TCGA[1:4,1:4]

setwd("D:\\St1datasets\\GSE30219")
GSE30219 <- read.csv("exprSet.csv", header = T, row.names = 1)


setwd("D:\\St1datasets\\GSE73403")
GSE73403 <- read.csv("exprSet.csv", header = T, row.names = 1)


SeleGenes <- intersect(rownames(TCGA), rownames(GSE30219))
SeleGenes <- intersect(SeleGenes, rownames(GSE73403))

data <- TCGA[which(rownames(TCGA) %in% SeleGenes),]
data[1:4,1:4]
dim(data)

#=======================================================

#=======================================================

#finding protein coding genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "transcript_biotype", "transcript_length"),
  filters = "biotype",
  values = "protein_coding",
  mart = mart
)

mRNA_genes <- subset(genes, transcript_biotype == "protein_coding")
head(mRNA_genes)
mRNA_genes$external_gene_name <- gsub("[^[:alnum:]]", "", mRNA_genes$external_gene_name)
save(mRNA_genes, file = "mRNA_genes.Rdata")

#=======================================================

#=======================================================

#selecting protein coding genes
load("mRNA_genes.Rdata")

protenGenes <- unique(mRNA_genes$external_gene_name)
data_sig <- data[which(rownames(data) %in% protenGenes),]
data_sig <- as.data.frame(t(data_sig))
data_sig[1:5,1:5]
dim(data_sig)

#=======================================================

#=======================================================


setwd("D:\\St2subtype")
siggenes <- read.csv("GeneCardsSearchResults.csv", header = T)
siggenes <- subset(siggenes, siggenes$Relevance.score > 1)
siggenes <- siggenes$Gene.Symbol
data_sig <- data_sig[,which(colnames(data_sig) %in% siggenes)]
dim(data_sig)

#=======================================================

#=======================================================

#remove constant expression genes
remove_constant_columns <- function(df, threshold = 0.8) {
  df[, sapply(df, function(col) {
    max(table(col) / length(col)) < threshold
  })]
}

data_sig_filtered <- remove_constant_columns(data_sig)
data_sig <- data_sig_filtered
data_sig[1:5,1:5]

min(data_sig)

#=======================================================

#=======================================================


preprocess_data <- apply(data_sig, 2, function(x) {
  q20 <- quantile(x, 0.2)
  q80 <- quantile(x, 0.8)
  x <- ifelse(x <= q20, 0, x)     
  x <- ifelse(x >= q80, 1, x)       
  x <- ifelse(x > q20 & x < q80, (x - q20) / (q80 - q20), x) 
  
  return(x)
})

data_binary <- as.data.frame(preprocess_data)

colnames(data_binary) <- colnames(data_sig)
rownames(data_binary) <- rownames(data_sig)

#=======================================================

#=======================================================

#univariate Cox analysis
data_binary$Sample <- rownames(data_binary)
data_loop <- merge(sur_data, data_binary, by="Sample")
rownames(data_loop) <- data_loop$Sample
data_loop$Sample <- NULL
data_loop[1:5,1:5]


surv_obj <- with(data_loop, Surv(OS_Time, OS))


results <- data.frame(
  Gene = character(), 
  HR = numeric(), 
  p_value = numeric(), 
  lower_95 = numeric(), 
  upper_95 = numeric(), 
  stringsAsFactors = FALSE
)


for (j in 3:ncol(data_loop)) {
  
  print(j)
  current_formula <- as.formula(paste("surv_obj ~", names(data_loop)[j]))
  fit <- coxph(current_formula, data = data_loop)
  summary_fit <- summary(fit)
  
  coef_info <- summary_fit$coefficients[1, ]  
  conf_int <- summary_fit$conf.int[1, ]
  
  results <- rbind(results, data.frame(
    Gene = names(data_loop)[j],
    HR = coef_info["exp(coef)"],
    p_value = coef_info["Pr(>|z|)"],
    lower_95 = conf_int["lower .95"],
    upper_95 = conf_int["upper .95"],
    stringsAsFactors = FALSE
  ))
}


head(results)


initial_threshold <- 0.05
results1 <- results[results$p_value < initial_threshold, ]

dim(results)
dim(results1)

#=======================================================

#=======================================================

setwd("D:\\St2subtype")

data_sub <- data_loop[,which(colnames(data_loop) %in% results1$Gene)]

write.csv(data_sub, file = "data_sub.csv")
