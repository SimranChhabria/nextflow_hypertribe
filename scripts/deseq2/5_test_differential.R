# Load libraries ---------------------------
# Import libraries
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--counts_file", nargs = "+", help = "Path to the counts file")
parser$add_argument("--cond", help = "Condition for deseq2")
parser$add_argument("--ctrl", help = "Condition for deseq2")

# Parse the arguments
args <- parser$parse_args()
counts_path_list  <- args$counts_file
control           <- args$ctrl
condition         <- args$cond


# Read formatted dataframe
counts_df <- read.csv((counts_path_list), check.names = FALSE)

# Perform statistical tests ---------------------------
#' Run DESeq2 on data for given comparison
#'
#' This function creates a DESeqDataSet object for the given dataframe,
#' runs DESeq on the object and retrive the results for the comparison. The
#' results are formatted and saved.
#'
#' @param ctrl_name The name of the control sample.
#' @param test_name The name of the test sample.
#' @param counts_df The counts dataframe.


run_deseq <- function(ctrl_name, test_name, counts_df) {

  #Filter count_df by condition and sample
  counts_ctrl_df <- counts_df %>% dplyr::select(starts_with(paste0(ctrl_name, "-")))
  counts_test_df <- counts_df %>% dplyr::select(starts_with(paste0(test_name, "-")))
  counts_comp_df <- cbind(counts_ctrl_df, counts_test_df)
  
  # Generate column metadata table
  col_comp_df <- DataFrame(sample=c(names(counts_ctrl_df),
                                    names(counts_test_df)),
                           condition=as.factor(c(rep("ctrl", ncol(counts_ctrl_df)),
                                                 rep("test", ncol(counts_test_df)))))
  
  # Generate DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_comp_df,
    colData = col_comp_df,
    design = ~condition
  )

  # Run DESeq
  dds <- DESeq(dds, quiet = TRUE)

  # Get results
  res <- results(dds, contrast = c("condition", "test", "ctrl"))

  # Get normalized counts
  norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE))
  ctrl_mean <- rowMeans(norm_counts_df %>% select(starts_with(ctrl_name)))
  test_mean <- rowMeans(norm_counts_df %>% select(starts_with(test_name)))

  # Format dataframe
  df <- counts_df[, 1:10]
  df$diff_mean <- res$log2FoldChange
  df$ctrl_mean <- ctrl_mean
  df$test_mean <- test_mean
  df$pval <- res$pvalue
  df$padj <- p.adjust(res$pvalue, "BH")
  df <- cbind(df, norm_counts_df)
  return (df)
  # Write results to file
  
}

# Run DESeq for different comparisons
df <- run_deseq(control, condition, counts_df)

file_name <- paste0(
    "5_diff_", condition, "_vs_", control, ".csv"
  )
write.csv(df, file_name, row.names = FALSE)