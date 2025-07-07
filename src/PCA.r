args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript pca_plot_script.R <input_file>")
}
file_dir <- args[1]

library(dplyr)
library(ggplot2)
library(factoextra)

data_matrix_df <- read.table(file_dir, header = TRUE, sep = "\t", row.names = 2, check.names = FALSE)
protein_data_df <- data_matrix_df[, 23:ncol(data_matrix_df)]

pca_result <- prcomp(protein_data_df, center = TRUE)

# PCA outlier detection
pca_coords <- as.data.frame(pca_result$x)
pca_coords$distance <- sqrt(pca_coords$PC1^2 + pca_coords$PC2^2)
pca_coords$sample <- rownames(pca_coords)
outliers <- pca_coords %>% arrange(desc(distance)) %>% head(5)

# Plotting function
plot_pca_by_group <- function(label, remission_col, response_col) {
  grouping <- data_matrix_df %>%
    transmute(Group = paste(.data[[remission_col]], .data[[response_col]], sep = "_")) %>%
    pull(Group)

  p <- fviz_pca_ind(pca_result,
                    geom.ind = "point",
                    col.ind = grouping,
                    palette = "jco",
                    addEllipses = TRUE,
                    legend.title = "Group") +
    geom_text(data = outliers,
              aes(x = PC1, y = PC2, label = sample),
              color = "black",
              size = 3,
              vjust = -1,
              aspect.ratio = 1) +
    theme_bw() +
    labs(
      title = paste0("PCA: ", label)
    )

  ggsave(file.path(dirname(file_dir), paste0("PCA/", label, ".png")),
         plot = p, width = 10, height = 7, dpi = 300)
}

# Run all 4 plots
plot_pca_by_group("QIDS_week4", "RemissionQIDS_W4", "QIDS_Response2_W4")
plot_pca_by_group("QIDS_week8", "RemissionQIDS_W8", "QIDS_Response2_W8")
plot_pca_by_group("HRSD_week4", "HRSDRemission_W4", "HRSD_Response2_W4")
plot_pca_by_group("HRSD_week8", "HRSDRemission_W8", "HRSD_Response2_W8")
