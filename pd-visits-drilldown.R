library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(data.table)

source("tidyDataframe.R")

ppmi_counts <- read.csv("ppmi_counts_matrix.csv")
bulk_expression_df <- tidyDataframe(ppmi_counts)

metadata<-read.csv("metaDataIR3.csv")

# Keep only passed samples - filter out the failed ones
metadata %>% filter(QCflagIR3 == 'pass') -> metadata
sani_samples <- str_replace_all(metadata$Sample, "PPMI.Phase\\d{1}.IR\\d{1}.", "")
sani_samples <- str_replace_all(sani_samples, "-", ".")
metadata$Sample <- sani_samples

bulk_expression_df <- cbind(bulk_expression_df, Diagnosis = metadata$DIAGNOSIS)
bulk_expression_df <- cbind(bulk_expression_df, Sample = metadata$Sample)
bulk_expression_df <- cbind(bulk_expression_df, Visit = metadata$CLINICAL_EVENT)
bulk_expression_df <- cbind(bulk_expression_df, Gender = metadata$GENDER)
bulk_expression_df <- cbind(bulk_expression_df, Patient = metadata$PATNO)

write.csv(bulk_expression_df, file="./ppmi_counts_meta_dataset.csv")

cohort_counts <- bulk_expression_df %>% group_by(Diagnosis) %>% summarise(Sample_Count = n(), .groups="drop")
participant_counts <- bulk_expression_df %>% group_by(Patient, Diagnosis) %>% summarise(Sample_Count = n(), .groups="drop")
cohort_counts
participant_counts

sample_summary <- participant_counts %>% group_by(Diagnosis, Sample_Count) %>% summarise(Num_Participants = n(), .groups="drop")
sample_summary

ggplot(sample_summary, aes(x = as.factor(Sample_Count), y = Num_Participants, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Number of samples per participant by cohort", x = "Samples per participant", y = "Count") +
  theme_minimal()

expression_df <- bulk_expression_df
expression_df$Sample <- rownames(bulk_expression_df)
multi_sampled <- expression_df %>% group_by(Patient) %>% filter(n() > 1)
multi_sampled %>% arrange(desc(Patient)) -> multi_sampled

# TODO
# maybe later, once genes with high variation emerge from PCA
ggplot(multi_sampled, aes(x = Visit, y = "LRRK2", fill = Visit)) +
  geom_boxplot() +
  labs(title = "Expression Variation Across Visits", x = "Visit", y = "LRRK2 Expression") +
  theme_minimal()

counts_cols <- !names(multi_sampled) %in% (c("Diagnosis", "Patient", "Sample", "Visit", "Gender"))
expression_matrix_df <- multi_sampled[, counts_cols]
expression_matrix_dt <- as.data.table(expression_matrix_df)
expression_matrix <- data.frame(transpose(expression_matrix_dt))
expression_matrix[] <- lapply(expression_matrix, type.convert, as.is=TRUE)
expression_matrix <- as.data.frame(expression_matrix)
rownames(expression_matrix) <- colnames(expression_matrix_df)
colnames(expression_matrix) <- rownames(expression_matrix_df)

non_zero_var_genes <- expression_matrix %>% filter(rowSums(across(everything(), ~ .x != 0)) > 0)

pca_res <- prcomp(non_zero_var_genes, scale.=FALSE, center=TRUE)

pca_df <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], Patient=multi_sampled$Patient)
