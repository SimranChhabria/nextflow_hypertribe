if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}

if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("plyr")
}

# Load libraries ---------------------------
# Import libraries
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(plyr)
library(stringr)

# Read parameters ---------------------------
args = commandArgs(trailingOnly=TRUE)
# lab <- args[1]
project <- "RiboSTAMP MOLM"
analysis <- args[1]
padj_threshold <- as.double(args[2])
diff_threshold <- as.double(args[3])
selection_name <- args[4]

# analysis <- "deseq2"
# padj_threshold <- as.double("0.05")
# diff_threshold <- as.double("1")
# selection_name <- NULL

# Load data ---------------------------
# Data directories
output_data_folder <- paste0("../../output_data/", analysis, "/5_test_differential/")
selection_df <- NULL
if (selection_name != "NULL") {
  selection_df <- read.csv(paste0("selection/", selection_name, ".csv"),
                           header = FALSE)
}

# Get list of result files
diff_file_list <- grep("_diff_", list.files(output_data_folder), value = TRUE)

# Read formatted output
formatted_output_file_list <- grep("_formatted_output", list.files(output_data_folder), value = TRUE)
formatted_output_df <- read.csv(paste0(output_data_folder, formatted_output_file_list[1]), check.names = FALSE)


# Plot PCA plot ---------------------------
if (analysis == "deseq2" | analysis == "erv") {
  suffix = "counts"

  raw_df <-t(as.matrix(formatted_output_df %>% dplyr::select(ends_with(suffix))))
  raw_df <- raw_df[ , which(apply(raw_df, 2, var) != 0)]
  pca <- prcomp(raw_df, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  var_list <- as.vector(summary(pca)$importance[2,1:2]) * 100
  
  pca_df$condition <- as.vector(sapply(row.names(pca_df), function(x) {tmp <- str_split(x, "_")[[1]]; paste(tmp[1:(length(tmp) - 2)], collapse = "_")}))
  pca_df$replicate <- as.vector(sapply(row.names(pca_df), function(x) {tmp <- str_split(x, "_")[[1]]; tmp[length(tmp) - 1]}))
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = condition, shape = replicate)) +
    ggtitle(paste0(project, " | ", analysis, " | ", suffix)) +
    xlab(paste0("PC1 (", var_list[1], "%)")) +
    ylab(paste0("PC2 (", var_list[2], "%)")) +
    theme(axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(paste0(output_data_folder, "pca.jpeg"),
         pca_plot, height = 7, width = 7)
}

# Plot volcano plot ---------------------------
#' Plot volcano plot for the differential analysis
#'
#' @param diff_df The differential testing results.
#' @param ctrl Name of the control sample.
#' @param test Name of the test sample.
#' @param analysis Type of analysis.
#' @param padj_threshold Adjusted p-value threshold.
#' @param diff_threshold Differential threshold.
#' @param selection Plot subtitle.
#'
plot_volcano <- function(diff_df,
                         ctrl,
                         test,
                         analysis,
                         padj_threshold,
                         diff_threshold,
                         subtitle) {
  
  # Filter out genes without gene symbols
  diff_df <- diff_df %>% filter(gene_symbol != "")
  
  # Filter for selected genes and select colors
  colCustom <- NULL
  col_list <- c("blue", "red", "orange", "green", "purple", "brown")
  if (!is.null(selection_df)) {
    diff_df <- merge(x = diff_df, y = selection_df,
                     by.x = "gene_symbol", by.y = "V1",
                     all.x = FALSE, all.y = TRUE)
    if (ncol(selection_df) > 1) {
      colCustom <- col_list[as.factor(diff_df$V2)]
      names(colCustom) <- diff_df$V2
      nonsig_list <- !(abs(diff_df$diff_mean) >= diff_threshold & diff_df$padj <= padj_threshold)
      colCustom[nonsig_list] <- "gray"
      names(colCustom)[nonsig_list] <- "Non-DE"
    }
  }
  
  # Select colors for ERVs
  if (analysis == "erv" & is.null(selection_df)) {
    sig_df <- diff_df %>% filter(abs(diff_df$diff_mean) >= diff_threshold & padj <= padj_threshold)
    feature_list <- unique(sig_df$feature)
    label_list <- diff_df$feature
    label_list[!label_list %in% feature_list] <- "Non-DE"
    nonsig_list <- !(abs(diff_df$diff_mean) >= diff_threshold & diff_df$padj <= padj_threshold)
    label_list[nonsig_list] <- "Non-DE"
    colCustom <- col_list[as.factor(label_list)]
    colCustom[nonsig_list] <- "gray"
    names(colCustom) <- label_list
  }
  
  # Determine p-value threshold
  pval_threshold <- tail(diff_df %>% arrange(pval) %>%
                         filter(padj <= padj_threshold) %>%
                         dplyr::select(pval), n=1)$pval
    
  # Determine plot configurations
  if (analysis == "hypertribe") {
    xlim = 1
    xlab = bquote(Freq[.(test)] - Freq[.(ctrl)])
    xtick = seq(-1, 1, 0.2)
  } else if (analysis == "deseq2" | analysis == "erv") {
    xlim = ceiling(max(abs(diff_df$diff_mean[!is.na(diff_df$diff_mean)])))
    xlab = bquote(log[2]*"FC ("*.(test)*"/"*.(ctrl)*")")
    if (xlim >= 10) {
      xlim = round_any(xlim, 2, f = ceiling)
      xtick = seq(-xlim, xlim, 2)
    } else {
      xtick = seq(-xlim, xlim, 1)
    }
  } else if (analysis == "labrat") {
    xlim = 1
    xlab = bquote(Psi[.(test)] - Psi[.(ctrl)])
    xtick = seq(-1, 1, 0.2)
  }
  
  if (analysis == "hypertribe") {
    analysis_name = "Diff. Editing"
  } else if (analysis == "deseq2") {
    analysis_name = "DGE"
  } else if (analysis == "erv") {
    analysis_name = "Diff. ERV Expression"
  } else if (analysis == "labrat") {
    analysis_name = "Diff. APA"
  }

  # Generate volcano plot
  #max_nlog_pval <- max(-log10(diff_df %>% filter(!is.na(pval)) %>% filter(abs(diff_mean) >= diff_threshold) %>% pull(pval)))
  max_nlog_pval <- max(-log10(diff_df %>% filter(!is.na(pval)) %>% pull(pval)))
  if (max_nlog_pval <= 20) {
    ylim <- round_any(max_nlog_pval, 2, f = ceiling)
    ytick <- seq(0, ylim, 2)
  } else if (max_nlog_pval <= 100) {
    ylim <- round_any(max_nlog_pval, 10, f = ceiling)
    ytick <- seq(0, ylim, 10)
  } else {
    ylim <- round_any(max_nlog_pval, 100, f = ceiling)
    ytick <- seq(0, ylim, 100)
  }
  p <- EnhancedVolcano(
    diff_df,
    lab = diff_df$gene_symbol,
    x = "diff_mean",
    y = "pval",
    xlab = xlab,
    ylab = expression(paste(-log[10], "(", italic(P), " value)")),
    colAlpha = 1,
    axisLabSize = 20,
    title = paste0(project, " | ", ctrl, " VS ", test, " | ", analysis_name),
    subtitle = subtitle,
    caption = NULL,
    titleLabSize = 20,
    subtitleLabSize = 16,
    pCutoff = pval_threshold,
    FCcutoff = diff_threshold,
    colCustom = colCustom,
    labSize = 6.5,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    border = "full",
    legendLabSize=20
  ) + 
    coord_cartesian(xlim=c(-xlim, xlim), ylim=c(0, ylim)) + 
    scale_x_continuous(breaks=xtick) +
    scale_y_continuous(breaks=ytick) +
    theme(axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          text = element_text(family = "Arial"))
  
  if (is.null(colCustom)) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom")
  }
  p$layers[[4]]$aes_params$fontface <- 1
  p$layers[[4]]$aes_params$family <- "Arial"
  p
}

# Plot volcano plots
for (i in c(1:length(diff_file_list))) {
  # Read table
  diff_file <- diff_file_list[i]
  diff_df <- read.csv(paste0(output_data_folder, diff_file), check.names = FALSE)
  min_pval <- min(diff_df %>% filter(!is.na(padj)) %>% filter(padj > 0) %>% pull(pval))
  min_padj <- min(diff_df %>% filter(!is.na(padj)) %>% filter(padj > 0) %>% pull(padj))
  if (nrow(diff_df %>% filter(padj == 0)) > 0) {
    diff_df <- diff_df %>% filter(padj == 0) %>% mutate(pval = min_pval) %>% mutate(padj = min_padj) %>% bind_rows(diff_df %>% filter(padj > 0))
  }
  
  # Set parameters
  tmp <- str_split(diff_file, "_vs_")[[1]]
  ctrl <- substr(tmp[1], start = 8, stop = str_length(tmp[1]))
  test <- substr(tmp[2], start = 1, stop = str_length(tmp[2]) - 4)
  padj_threshold <- padj_threshold
  diff_threshold <- diff_threshold
  plot_path <- paste0(output_data_folder, ctrl, "_vs_", test, "_volcano.jpeg")
  
  # Plot volcano plot
  if (is.null(selection_df)) {
    subtitle <- paste0("diff_mean \u2265 ", diff_threshold, " AND padj \u2264 ", padj_threshold)
  } else {
    subtitle <- paste0("diff_mean \u2265 ", diff_threshold, " AND padj \u2264 ", padj_threshold, " AND SELECT ", selection_name)
  }
  volcano_plot <- plot_volcano(diff_df,
                               ctrl,
                               test,
                               analysis,
                               padj_threshold,
                               diff_threshold,
                               subtitle)
  ggsave(paste0(output_data_folder, ctrl, "_vs_", test, "_volcano.jpeg"),
         volcano_plot, height = 10, width = 10)
}





