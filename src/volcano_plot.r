# MAKE_volcano_plot_Rscript
# This script generates volcano plots from linear regression results

library(ggplot2)    # For plotting
library(ggrepel)    # For non-overlapping text labels


output_dir <- "/Users/m221138/proteomics_core/analysis/P24-188"
if (!dir.exists(file.path(output_dir, "volcano_plot"))) {
    dir.create(file.path(output_dir, "volcano_plot"))
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript volcano_plot.r <column of interest> <group_a_str> <group_b_str>")
}

data_file <- paste0("/Users/m221138/proteomics_core/analysis/P24-188/linear_model/linear_regression.", args[1], ".tsv")
condition_a <- args[2]  # e.g., "Response"
condition_b <- args[3]  # e.g., "Non-Response"

input_df = read.csv(data_file, sep="\t", header=TRUE, row.names=1)

# Set x and y axes for the volcano plot
x_axis <- input_df$log2_fold_change
y_axis <- -log10(input_df$p_value_of_class_adj_covariates)  # Convert p-values to –log10 scale

# Set axis labels and colors based on condition
xaxis_label <- paste0("log2(geometric mean ", condition_a, " / geometric mean ", condition_b, ")")
numerator_color = "firebrick4"        
denominator_color = "steelblue4"      

gene_list <- rownames(input_df)  # Feature names
df <- data.frame(
    log2_fold_change = x_axis,
    log10_p_value_of_class_adj_covariates = y_axis,
    row.names = gene_list
)
print (dim(df))
df$genes <- gene_list  # add gene names as a column too (for ggplot)

# Set y-axis limits and p-value threshold line
log2pval_threshold = 1.30103  # = -log10(0.05)

# Filter significant features (above p-value threshold)
sig_subset <- subset(df, log10_p_value_of_class_adj_covariates > log2pval_threshold)

# Subset for up/down regulated features
sig_red_subset <- subset(sig_subset, log2_fold_change > 0)
sig_blue_subset <- subset(sig_subset, log2_fold_change < -0)

# Same subset for labeling
sig_red_text_subset <- sig_red_subset
sig_blue_text_subset <- sig_blue_subset

# Debug printout
print ('---------------')
print (paste("up: ", nrow(sig_red_subset), sep=""))
print (paste("down: ", nrow(sig_blue_subset), sep=""))
print ('---------------')
        
# Set figure title and output path
figure_title = paste(args[1], sep="")
output_pdf <- paste0("/Users/m221138/proteomics_core/analysis/P24-188/volcano_plot/volcano_plot.", args[1], ".pdf")

# # Save volcano plot to PDF
pdf(output_pdf)
# Draw volcano plot

xlim_range <- max(abs(df$log2_fold_change)) 

plot_pdf <- ggplot(df, aes(x=log2_fold_change, y=log10_p_value_of_class_adj_covariates)) + 
    xlim(-xlim_range, xlim_range) +
    geom_point(colour="#DCDCDC", size = 2.5, stroke = 0, alpha = 0.5) +  # Background dots
    geom_hline(yintercept = log2pval_threshold, colour="#BEBEBE", linetype="dashed") +
    geom_point(data = sig_red_subset, colour=numerator_color, size = 2.5, stroke = 0, alpha = 0.7) +
    geom_point(data = sig_blue_subset, colour=denominator_color, size = 2.5, stroke = 0, alpha = 0.5) +
    geom_text_repel(data = sig_red_text_subset, aes(log2_fold_change, log10_p_value_of_class_adj_covariates, label=genes),
                    colour=numerator_color, size=2) +
    geom_text_repel(data = sig_blue_text_subset, aes(log2_fold_change, log10_p_value_of_class_adj_covariates, label=genes),
                    colour=denominator_color, size=2) +
    ylab("-log10 (P-value)") + xlab(xaxis_label) +  
    ggtitle(figure_title) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

print (plot_pdf)
dev.off()



