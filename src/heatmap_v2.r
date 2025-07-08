library(pheatmap)

output_dir <- "/Users/m221138/proteomics_core/analysis/P24-188"
if (!dir.exists(file.path(output_dir, "heatmap"))) {
    dir.create(file.path(output_dir, "heatmap"))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript heatmap.r <column of interest> <group_a_str> <group_b_str>")
}

class_column <- args[1]

# Load data
regression_file <- paste0("/Users/m221138/proteomics_core/analysis/P24-188/linear_model/linear_regression.", class_column, ".tsv")
regression_df <- read.csv(regression_file, sep="\t", row.names=1, check.names=FALSE)

concat_file <- "/Users/m221138/proteomics_core/analysis/P24-188/concat_data_matrix.tsv"
concat_df <- read.csv(concat_file, sep="\t", row.names=1, check.names=FALSE)

# Filter significant proteins
sig_df <- regression_df[regression_df$p_value_of_class_adj_covariates < 0.05, ]

# Top 5 positive fold-change proteins
pos_df <- sig_df[sig_df$log2_fold_change > 0, ]
top_pos_genes <- rownames(head(pos_df[order(pos_df$p_value_of_class_adj_covariates), ], 5))

# Top 5 negative fold-change proteins
neg_df <- sig_df[sig_df$log2_fold_change < 0, ]
top_neg_genes <- rownames(head(neg_df[order(neg_df$p_value_of_class_adj_covariates), ], 5))

# Combine selected genes
top_genes <- unique(c(top_pos_genes, top_neg_genes))
cat("Selected", length(top_pos_genes), "positive and", length(top_neg_genes), "negative DA proteins.\n")

# Subset expression data and transpose
protein_cols <- intersect(colnames(concat_df), top_genes)
heatmap_data <- t(concat_df[, protein_cols])

# Z-score normalization per protein
heatmap_data_zscore <- t(scale(t(heatmap_data)))

# Clip Z-scores to ±3 for visualization
heatmap_data_zscore[heatmap_data_zscore > 3] <- 3
heatmap_data_zscore[heatmap_data_zscore < -3] <- -3

# Create annotation dataframe
annotation_col <- data.frame(as.factor(concat_df[[class_column]]))
colnames(annotation_col) <- class_column
rownames(annotation_col) <- rownames(concat_df)

# Reorder columns by class
grouping <- annotation_col[[class_column]]
ordered_grouping <- factor(grouping, levels = unique(grouping))
col_order <- order(ordered_grouping)
heatmap_data_zscore <- heatmap_data_zscore[, col_order]
annotation_col <- annotation_col[col_order, , drop=FALSE]

# Column gap (between group A and B)
gaps_col <- which(diff(as.numeric(ordered_grouping[col_order])) != 0)

# Annotation colors
group_levels <- levels(annotation_col[[class_column]])
group_colors <- setNames(c("#B57623", "#78AF3F")[seq_along(group_levels)], group_levels)
annotation_colors <- setNames(list(group_colors), class_column)

# Color palette and breaks
my_palette <- colorRampPalette(c("steelblue4", "white", "firebrick4"))(100)
my_breaks <- seq(-3, 3, length.out = 101)

# Output file
figure_title <- paste0("Top DA Proteins (", length(top_pos_genes), " UP-regulated / ", length(top_neg_genes), " DOWN-regulated)")
output_pdf <- paste0(output_dir, "/heatmap/heatmap.", class_column, ".pdf")

# Save heatmap to PDF
pdf(output_pdf, width = 12, height = 6)

pheatmap(heatmap_data_zscore,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         fontsize_row = 6,
         fontsize_col = 1,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         gaps_col = gaps_col,
         scale = "none",
         color = my_palette,
         breaks = my_breaks,
         main = figure_title)

dev.off()