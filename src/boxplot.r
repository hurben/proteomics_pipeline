library(ggplot2)
library(dplyr)
library(tidyr)

library(ggtext)
library(tibble)


output_dir <- "/Users/m221138/proteomics_core/analysis/P24-188"
if (!dir.exists(file.path(output_dir, "boxplots"))) {
    dir.create(file.path(output_dir, "boxplots"))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript boxplot.r <class_column>")
}

class_column <- args[1]

# Load data
regression_file <- file.path(output_dir, "linear_model", paste0("linear_regression.", class_column, ".tsv"))
regression_df <- read.csv(regression_file, sep = "\t", row.names = 1, check.names = FALSE)

concat_file <- file.path(output_dir, "concat_data_matrix.tsv")
concat_df <- read.csv(concat_file, sep = "\t", row.names = 1, check.names = FALSE)

# Filter significant proteins
sig_df <- regression_df %>% filter(p_value_of_class_adj_covariates < 0.05)

# Top 5 up- and down-regulated proteins
top_pos <- sig_df %>% filter(log2_fold_change > 0) %>% arrange(p_value_of_class_adj_covariates) %>% head(5)
top_neg <- sig_df %>% filter(log2_fold_change < 0) %>% arrange(p_value_of_class_adj_covariates) %>% head(5)
top_genes <- unique(c(rownames(top_pos), rownames(top_neg)))

cat("Selected", nrow(top_pos), "positive and", nrow(top_neg), "negative DA proteins.\n")

# Extract protein expression
protein_cols <- intersect(colnames(concat_df), top_genes)
expr_df <- concat_df[, protein_cols, drop = FALSE]
expr_df[[class_column]] <- as.factor(concat_df[[class_column]])
expr_df$SampleID <- rownames(concat_df)

# Reshape to long format
long_df <- expr_df %>%
  pivot_longer(cols = all_of(top_genes), names_to = "Protein", values_to = "Abundance")

# Prepare colored labels with P-values
label_df <- regression_df[rownames(regression_df) %in% top_genes, , drop = FALSE]
label_df$direction <- ifelse(label_df$log2_fold_change > 0, "up", "down")
label_df$formatted_p <- sprintf("%.3g", label_df$p_value_of_class_adj_covariates)

label_df$colored_label <- ifelse(
  label_df$log2_fold_change > 0,
  paste0("<span style='color:red'>", rownames(label_df), " (P = ", label_df$formatted_p, ")</span>"),
  paste0("<span style='color:blue'>", rownames(label_df), " (P = ", label_df$formatted_p, ")</span>")
)

# Merge colored labels into long_df
label_df <- label_df %>% rownames_to_column("Protein")
long_df <- long_df %>%
  left_join(label_df %>% select(Protein, colored_label), by = "Protein")

# Plot
pdf(file = file.path(output_dir, "boxplots", paste0("boxplots.", class_column, ".pdf")), width = 12, height = 6)

ggplot(long_df, aes(x = .data[[class_column]], y = Abundance, fill = .data[[class_column]])) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(aes(color = .data[[class_column]]), width = 0.2, size = 0.8, alpha = 0.6) +
  facet_wrap(~ colored_label, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("#B57623", "#78AF3F")) +
  scale_color_manual(values = c("#B57623", "#78AF3F")) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_markdown(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = paste0("Top Differentially Abundant Proteins (", length(top_genes), " total)"),
       x = class_column, y = "Protein Abundance")

dev.off()