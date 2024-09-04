# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# if (!requireNamespace("DESeq2", quietly = TRUE)) {
#   BiocManager::install("DESeq2")
# }

# if (!requireNamespace("dplyr", quietly = TRUE)) {
#   install.packages("dplyr")
# }

# if (!requireNamespace("argparse", quietly = TRUE)) {
#   install.packages("argparse")
# }

# if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
#   BiocManager::install("EnhancedVolcano")
# }

# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#   install.packages("ggplot2")
# }

# if (!requireNamespace("stringr", quietly = TRUE)) {
#   install.packages("stringr")
# }

# if (!requireNamespace("stringr", quietly = TRUE)) {
#   install.packages("plyr")
# }

# Load libraries ---------------------------
# Import libraries
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(plyr))

parser <- ArgumentParser()
parser$add_argument("--deseq_file", nargs = "+", help = "Path to the differential analysis results file")
parser$add_argument("--formatted_file", nargs = "+", help = "Path to the formatted results file")
parser$add_argument("--cond", help = "Condition for deseq2")
parser$add_argument("--ctrl", help = "Condition for deseq2")
parser$add_argument("--analysis", help = "Type of Analysis (deseq2 or hypertribe)")
parser$add_argument("--padj_threshold", help = "Threshold for padj")
parser$add_argument("--diff_threshold", help = "Threshold for log2FC")
parser$add_argument("--selection_name", help = "Subsetting based on genes")
parser$add_argument("--project_name", help= "Name of the project")

# Parse the arguments
args <- parser$parse_args()
deseq_file        <- args$deseq_file
formatted_file    <- args$formatted_file
ctrl              <- args$ctrl
test              <- args$cond
analysis          <- args$analysis
padj_threshold    <- as.double(args$padj_threshold)
diff_threshold    <- as.double(args$diff_threshold)
selection_name    <- args$selection_name
project           <- args$project_name


#analysis <- "deseq2"
#padj_threshold <- as.double("0.05")
#diff_threshold <- as.double("1")
#selection_name <- "NULL"
#project <- "project_15619"


# Read formatted dataframe
selection_df <- NULL
if (selection_name != "NULL") {
  selection_df <- read.csv(paste0(selection_name, ".csv"),
                           header = FALSE)
}

formatted_output_df <- read.csv((formatted_file),check.names = FALSE)
diff_df <- read.csv((deseq_file),check.names = FALSE)

# Plot PCA plot ---------------------------
if (analysis == "deseq2" | analysis == "erv") {
  suffix = "counts"
  
  raw_df <-t(as.matrix(formatted_output_df %>% dplyr::select(ends_with(suffix))))
  raw_df <- raw_df[ , which(apply(raw_df, 2, var) != 0)]
  pca <- prcomp(raw_df, center = TRUE, scale. = TRUE)
  
  raw_pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df <- raw_pca_df[grepl(paste0("^", ctrl), row.names(raw_pca_df)) | grepl(paste0("^", test), row.names(raw_pca_df)), ]
  var_list <- as.vector(summary(pca)$importance[2,1:2]) * 100
  
  #pca_df$condition <- as.vector(sapply(row.names(pca_df), function(x) {tmp <- str_split(x, "_")[[1]]; paste(tmp[1:(length(tmp) - 2)], collapse = "_")}))
  pca_df$condition <- as.vector(sapply(row.names(pca_df), function(x) {
    tmp <- str_split(x, "_")[[1]]
    base_name <- paste(tmp[1:(length(tmp) - 2)], collapse = "_")
    str_remove(base_name, "-\\d+$")}))
  
  pca_df$replicate <- as.vector(sapply(row.names(pca_df), function(x) {tmp <- str_split(x, "_")[[1]]; tmp[length(tmp) - 1]}))
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) + 
    geom_point(aes(color = condition, shape = replicate)) +
    ggtitle(paste0(analysis, " | ", suffix)) +
    xlab(paste0("PC1 (", var_list[1], "%)")) +
    ylab(paste0("PC2 (", var_list[2], "%)")) +
    theme(axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(paste0(project, "-", test, "-", ctrl, "-", analysis, "-pca.jpeg"),
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
    title = paste0(project, " | ", test, " VS ", ctrl, " | ", analysis_name),
    #subtitle = subtitle,
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

#--- Volcano Plots ---#s
min_pval <- min(diff_df %>% filter(!is.na(padj)) %>% filter(padj > 0) %>% pull(pval))
min_padj <- min(diff_df %>% filter(!is.na(padj)) %>% filter(padj > 0) %>% pull(padj))
if (nrow(diff_df %>% filter(padj == 0)) > 0) {
  diff_df <- diff_df %>% filter(padj == 0) %>% mutate(pval = min_pval) %>% mutate(padj = min_padj) %>% bind_rows(diff_df %>% filter(padj > 0))
}
padj_threshold <- padj_threshold
diff_threshold <- diff_threshold

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
ggsave(paste0(project, "_", test, "_vs_", ctrl,  "-", analysis, "_volcano.jpeg"),
       volcano_plot, height = 10, width = 10)