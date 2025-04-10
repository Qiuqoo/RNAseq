
library(Hmisc)
library(dplyr)


generate_interactions <- function(fullExpr, targetExpr, threshold_r = 0.8, threshold_p = 0.01, topN = NULL) {

  common_samples <- intersect(colnames(fullExpr), colnames(targetExpr))
  if(length(common_samples) == 0) {
    stop("fullExpr and targetExpr do not have common samples, please check the data!")
  }
  

  fullExpr <- fullExpr[, common_samples, drop = FALSE]
  targetExpr <- targetExpr[, common_samples, drop = FALSE]
  

  fullExpr <- fullExpr[apply(fullExpr, 1, sd) > 0, ]
  

  cor_res <- rcorr(as.matrix(t(fullExpr)), as.matrix(t(targetExpr)), type = "pearson")
  cor_matrix <- cor_res$r
  p_matrix   <- cor_res$P
  

  target_genes <- rownames(targetExpr)
  

  coexpressed_genes <- list()
  

  for (gene in target_genes) {

    gene_cor <- cor_matrix[, gene]
    gene_p   <- p_matrix[, gene]
    

    selected <- names(which(abs(gene_cor) > threshold_r & gene_p < threshold_p))
    selected <- setdiff(selected, target_genes)
    
    if(length(selected) > 0) {
      gene_interactions <- data.frame(
        TargetGene      = gene,
        CoexpressedGene = selected,
        Correlation     = gene_cor[selected],
        PValue          = gene_p[selected],
        stringsAsFactors = FALSE
      )
      coexpressed_genes[[gene]] <- gene_interactions
    }
  }
  
  all_interactions <- do.call(rbind, coexpressed_genes)
  all_interactions <- unique(all_interactions)
  
  if(!is.null(topN)) {
    sorted_interactions <- all_interactions[order(-abs(all_interactions$Correlation)), ]
    all_interactions <- head(sorted_interactions, topN)
  }
  
  return(all_interactions)
}





WOXExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/WOX.csv", row.names = 1)
fullExpr <- read.csv("/data3/users/Qiuqoo/CM_ZJU/article/GRN/RNAseq_counts_mean.csv", row.names = 1)


result_interactions <- generate_interactions(fullExpr, WOXExpr, threshold_r = 0.8, threshold_p = 0.01, topN = 300)


head(result_interactions)
dim(result_interactions)

write.csv(result_interactions, file = "/data3/users/Qiuqoo/CM_ZJU/article/GRN/WOX_top300_interactions.csv", row.names = FALSE)

