# Load necessary libraries
library(lme4)          # For linear mixed-effects models (not used in this script)
library(lmerTest)      # Adds p-values to lmer models (not used in this script)
library(stringr)       # String handling (not used in this script)
library(effects)       # For visualizing model effects (not used here)
library(dplyr)         # For data wrangling (used for filtering)

main <- function (tmp_df, column_of_interest, group_a_str, group_b_str){
    output_dir <- "/Users/m221138/proteomics_core/analysis/P24-188"
    if (!dir.exists(file.path(output_dir, "linear_model"))) {
        dir.create(file.path(output_dir, "linear_model"))
    }

    results <- data.frame()

    # Convert relevant covariates to factors once
    tmp_df$GENDER_CODE <- as.factor(tmp_df$GENDER_CODE)
    tmp_df$Med_Start_Med <- as.factor(tmp_df$Med_Start_Med)

    tmp_df[[column_of_interest]] <- as.factor(tmp_df[[column_of_interest]])

    for (i in 23:ncol(tmp_df)) {
        feature_name <- colnames(tmp_df)[i]
        y <- tmp_df[[feature_name]]

        # Define comparison groups dynamically
        group_a <- tmp_df[[column_of_interest]] == group_a_str
        group_b <- tmp_df[[column_of_interest]] == group_b_str
        condition_a_list <- tmp_df[group_a, i]
        condition_b_list <- tmp_df[group_b, i]

        log2fc <- mean(condition_a_list, na.rm = TRUE) - mean(condition_b_list, na.rm = TRUE)

        # Marginal model
        formula_marginal <- as.formula(paste("y ~", column_of_interest))
        marginal_model <- lm(formula_marginal, data = tmp_df)
        raw_pval <- coef(summary(marginal_model))[,4][2]
        raw_coef <- coef(summary(marginal_model))[,1][2]

        # Full model
        full_formula <- paste(
        "y ~", column_of_interest,
        "+ AGE + GENDER_CODE + Med_Start_Med + `B_HRS-D17/QIDS16_HRS-D17 Score` + `B_HRS-D17/QIDS16_QIDS16 Score`"
        )
        final_model <- lm(as.formula(full_formula), data = tmp_df)
        adj_pval <- coef(summary(final_model))[,4][2]
        adj_coef <- coef(summary(final_model))[,1][2]

        results <- rbind(results, data.frame(
        feature = feature_name,
        log2_fold_change = log2fc,
        coefficients_of_class= raw_coef,
        p_value_of_class = raw_pval,
        coefficients_of_class_adj_covariates = adj_coef,
        p_value_of_class_adj_covariates = adj_pval,
        stringsAsFactors = FALSE
        ))
    }
    # Adjust p-values
    results$FDR_value_of_class_adj_covariates <- p.adjust(results$p_value_of_class_adj_covariates, method = "fdr")
    output_file <- paste0(output_dir, "/linear_model/linear_regression.", column_of_interest, ".tsv")
    write.table(results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

input_file <- "/Users/m221138/proteomics_core/analysis/P24-188/concat_data_matrix.tsv"
input_df <- read.csv(input_file, header = TRUE, sep = "\t", check.names = FALSE)
tmp_df <- input_df 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript linear_regression.r <column of interest> <group_a_str> <group_b_str>")
}
column_of_interest <- args[1]
group_a_str <- args[2]
group_b_str <- args[3]

main(tmp_df, column_of_interest, group_a_str, group_b_str)

