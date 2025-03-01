source("convertToNumeric.R")
library(dplyr)
library(ggrepel)
library(pheatmap)

bulk_expression_df_micro <- bulk_expression_df[1:100, 58685:58785]
bulk_counts <- bulk_expression_df_micro[, 1:96]
bulk_counts_t <- t(bulk_counts)
View(convertToNumeric(bulk_counts_t))

bulk_expression_df_micro_t <- bulk_expression_df_micro[, 1:96] %>% t %>% convertToNumeric %>% filter(rowSums(across(everything(), ~ .x != 0)) > 0)
View(bulk_expression_df_micro_t)

pca_res <- prcomp(t(bulk_expression_df_micro_t), center=TRUE, scale.=TRUE)

pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], 
                     Patient = bulk_expression_df_micro[,"Patient"], 
                     Visit = bulk_expression_df_micro[, "Visit"], 
                     Diagnosis = bulk_expression_df_micro[, "Diagnosis"])

pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=as.factor(Visit), shape=as.factor(Diagnosis))) +
              geom_point(size=3, alpha=0.7) +
              labs(title = "PCA Gene Expression across Visits",
              x = "PC1",
              y = "PC2",
              color = "Visit",
              shape = "Diagnosis") + theme_minimal()

pca_plot

pca_plot + geom_polygon(data = pca_df %>% group_by(Patient) %>% slice(chull(PC1, PC2)), 
                        aes(group=Patient), fill=NA, color="gray50", linetype="dashed", size=0.1) + 
  labs(title = "PCA - sample points from same participants")

pca_plot + stat_ellipse(aes(fill=as.factor(Visit)), alpha=0.1, geom="polygon") + 
  labs(title = "PCA - spread per visit category")

pca_res <- prcomp(t(bulk_expression_df_micro_t), center=TRUE, scale.=TRUE)
gene_loadings <- data.frame(Gene = rownames(pca_res$rotation),
                            PC1 = abs(pca_res$rotation[,1]),
                            PC2 = abs(pca_res$rotation[,2]))

top_genes <- gene_loadings %>% rowwise() %>% mutate(Max_Contribution = max(PC1, PC2)) %>%
  arrange(desc(Max_Contribution))

top_20 <- top_genes[1:20,]

ggplot(top_20, aes(x = reorder(Gene, -Max_Contribution), y = Max_Contribution)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + 
  labs("Top Variable Genes", x = "Gene", y = "Max Contribution") + theme_minimal()

subset <- bulk_expression_df_micro[,top_20$Gene]
subset <- cbind(subset, Diagnosis = bulk_expression_df_micro[,"Diagnosis"])
subset <- cbind(subset, Visit = bulk_expression_df_micro[,"Visit"])
subset <- cbind(subset, Patient = bulk_expression_df_micro[,"Patient"])
subset <- cbind(subset, Gender = bulk_expression_df_micro[,"Gender"])
View(subset)

subset_counts <- subset[,1:20]
subset_scaled <- scale(t(convertToNumeric(subset_counts)))
pheatmap(subset_scaled, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         scale="none", 
         show_rownames=TRUE, 
         show_colnames=FALSE, 
         color=colorRampPalette(c("blue", "white", "red"))(50),
         main="Heatmap of most variable genes in PCA")

# bulk_expression_df_micro <- bulk_expression_df[1:100, 58685:58785]
# 
# bulk_expression_df_micro <- convertToNumeric(bulk_expression_df_micro)
# 
# cohort_counts <- bulk_expression_df_micro %>% group_by(Diagnosis) %>% summarise(Sample_Count = n(), .groups="drop")
# participant_counts <- bulk_expression_df_micro %>% group_by(Patient, Diagnosis) %>% summarise(Sample_Count = n(), .groups="drop")
# cohort_counts
# participant_counts
# 
# sample_summary <- participant_counts %>% group_by(Diagnosis, Sample_Count) %>% summarise(Num_Participants = n(), .groups="drop")
# sample_summary
# 
# ggplot(sample_summary, aes(x = as.factor(Sample_Count), y = Num_Participants, fill = Diagnosis)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   labs(title = "Number of samples per participant by cohort", x = "Samples per participant", y = "Count") +
#   theme_minimal()
# 
# expression_df <- as.data.frame(bulk_expression_df_micro)
# multi_sampled <- expression_df %>% group_by(Patient) %>% filter(n() > 1)
# multi_sampled %>% arrange(desc(Patient)) -> multi_sampled
# 
# counts_cols <- !names(multi_sampled) %in% (c("Diagnosis", "Patient", "Sample", "Visit", "Gender"))
# expression_matrix_df <- multi_sampled[, counts_cols]
# 
# non_zero_var_genes <- multi_sampled %>% filter(rowSums(across(everything(), ~ .x != 0)) > 0)
# multi_sampled_t <- convertToNumeric(expression_matrix_df) %>% t
# pca_res <- prcomp(multi_sampled_t, scale.=TRUE, center=TRUE)
# pca_df <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], Gene = colnames(multi_sampled[,1:96]))
# 
# View(pca_df)
