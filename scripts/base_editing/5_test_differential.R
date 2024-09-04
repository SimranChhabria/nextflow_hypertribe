# Load libraries ---------------------------
# Import libraries
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))
library(bbmle)
library(parallel)
library(VGAM)

parser <- ArgumentParser()
parser$add_argument("--formatted_file", nargs = "+", help = "Path to the counts file")
parser$add_argument("--cond", help = "Condition for deseq2")
parser$add_argument("--ctrl", help = "Condition for deseq2")


# Parse the arguments
args <- parser$parse_args()
formatted_path_list  <- args$formatted_file
ctrl_name            <- args$ctrl
test_name            <- args$cond


formatted_df <- read.csv((formatted_path_list), check.names = FALSE)


# Likelihood models ---------------------------
#' Likelihood model for the null hypothesis
#'
#' This function computes the likelihood of the data under the null hypothesis.
#'
#' @param ref The number of edits in the control samples.
#' @param alt The number of edits in the treatment samples.
#'
#' @return The likelihood.
mle_custom_h0 <- function(ref, alt, debug = FALSE) {
  x <- alt
  size <- ref + alt

  bbll <- function(prob, rho) {
    if ((prob > 0) & (rho > 0) & (prob < 1) & (rho < 1)) {
      -sum(dbetabinom(x, size, prob, rho, log = TRUE))
    } else {
      NA
    }
  }

  fit <- mle2(bbll,
    start = list(prob = .1, rho = .5),
    method = "Nelder-Mead",
    control = list(maxit = 1e5, trace = as.integer(debug))
  )
}

#' Likelihood model for the alternative hypothesis
#'
#' This function computes the likelihood of the data under the alternative
#' hypothesis.
#'
#' @param ref The number of edits in the control samples.
#' @param alt The number of edits in the treatment samples.
#'
#' @return The likelihood.
mle_custom_h1 <- function(ref1, alt1, ref2, alt2, debug = FALSE) {
  x1 <- alt1
  size1 <- ref1 + alt1
  prob1_init <- mean(x1 / size1) + 1e-3
  if (is.na(prob1_init) | prob1_init >= 1 | prob1_init <= 0) {
    prob1_init <- .05
  }
  x2 <- alt2
  size2 <- ref2 + alt2
  prob2_init <- mean(x2 / size2) + 1e-3

  if (is.na(prob2_init) | prob2_init >= 1 | prob2_init <= 0) {
    prob2_init <- .05
  }

  bbll <- function(prob1, prob2, rho) {
    if ((prob1 > 0) & (prob2 > 0) & (rho > 0) &
      (prob1 < 1) & (prob2 < 1) & (rho < 1)) {
      -(sum(dbetabinom(x1, size1, prob1, rho, log = TRUE)) +
        sum(dbetabinom(x2, size2, prob2, rho, log = TRUE)))
    } else {
      NA
    }
  }

  fit <- mle2(bbll,
    start = list(prob1 = prob1_init, prob2 = prob2_init, rho = .1),
    method = "Nelder-Mead",
    control = list(maxit = 1e5, trace = as.integer(debug))
  )
}

# Computing the likelihood of data ---------------------------
#' Compute the differential likelihood of the data to be under the alternative
#' hypothesis
#'
#' This function computes the differential likelihood of the data to be under
#' the alternative hypothesis.
#'
#' @param i Row number of the formatted dataframe.
#' @param ctrl_ref_df Number of reference alleles in the control samples
#' @param ctrl_alt_df Number of alternative alleles in the control samples
#' @param test_ref_df Number of reference alleles in the treatment samples
#' @param test_alt_df Number of alternative alleles in the treatment samples
#'
#' @return The p value.
fit_mle <- function(i, ctrl_ref_df, ctrl_alt_df, test_ref_df, test_alt_df) {

  # Get the counts of row i
  ctrl_ref_list <- as.numeric(as.vector(ctrl_ref_df[i,]))
  ctrl_alt_list <- as.numeric(as.vector(ctrl_alt_df[i,]))
  test_ref_list <- as.numeric(as.vector(test_ref_df[i,]))
  test_alt_list <- as.numeric(as.vector(test_alt_df[i,]))

  # Null hypothesis: ctrl = test
  fit0 <- mle_custom_h0(c(ctrl_ref_list, test_ref_list),
                        c(ctrl_alt_list, test_alt_list))

  # Tested hypothesis: ctrl < test
  fit1 <- mle_custom_h1(ctrl_ref_list, ctrl_alt_list,
                        test_ref_list, test_alt_list)

  # Compute p value, fold change and convergence
  pval <- pchisq(2 * (fit0@min - fit1@min), 1, lower.tail = FALSE)
  value <- fit1@min
  c(pval = pval, value = value)
}

col_list <- names(formatted_df)
ctrl_ref_col_list <- col_list[startsWith(col_list, ctrl_name) & endsWith(col_list, "_ref_counts")]
ctrl_alt_col_list <- col_list[startsWith(col_list, ctrl_name) & endsWith(col_list, "_alt_counts")]
test_ref_col_list <- col_list[startsWith(col_list, test_name) & endsWith(col_list, "_ref_counts")]
test_alt_col_list <- col_list[startsWith(col_list, test_name) & endsWith(col_list, "_alt_counts")]
  
# Select sites that have edits
ref_df <- cbind(formatted_df[ctrl_ref_col_list], formatted_df[test_ref_col_list])
alt_df <- cbind(formatted_df[ctrl_alt_col_list], formatted_df[test_alt_col_list])
mask_list <- rowSums(ref_df) != 0 & rowSums(alt_df != 0) > 1 & (rowMeans(ref_df) + rowMeans(alt_df) >= 5)
ctrl_ref_df <- formatted_df[mask_list, ctrl_ref_col_list, drop = FALSE]
ctrl_alt_df <- formatted_df[mask_list, ctrl_alt_col_list, drop = FALSE]
test_ref_df <- formatted_df[mask_list, test_ref_col_list, drop = FALSE]
test_alt_df <- formatted_df[mask_list, test_alt_col_list, drop = FALSE]
  
# Run statistical test
res_list <- mclapply(seq_len(nrow(ctrl_ref_df)),
                       fit_mle,
                       ctrl_ref_df, ctrl_alt_df,
                       test_ref_df, test_alt_df,
                       mc.cores = 32)
res_df <- data.frame(do.call("rbind", res_list))
pos_mask_list <- res_df$value > 0
  
# Compute editing frequency
ctrl_freq_df <- (ctrl_alt_df / (ctrl_ref_df + ctrl_alt_df))[pos_mask_list,, drop = FALSE]
names(ctrl_freq_df) <- paste0(ctrl_name, "_", c(1:length(ctrl_ref_col_list)), "_freq")
test_freq_df <- (test_alt_df / (test_ref_df + test_alt_df))[pos_mask_list,, drop = FALSE]
names(test_freq_df) <- paste0(test_name, "_", c(1:length(test_ref_col_list)), "_freq")
ctrl_mean_list <- rowMeans(ctrl_freq_df, na.rm = TRUE)
test_mean_list <- rowMeans(test_freq_df, na.rm = TRUE)
stats_df <- data.frame(
    diff_mean = test_mean_list - ctrl_mean_list,
    ctrl_mean = ctrl_mean_list,
    test_mean = test_mean_list,
    pval = res_df$pval[pos_mask_list],
    padj = p.adjust(res_df$pval[pos_mask_list], "BH")
  )
stats_df <- cbind(stats_df, ctrl_freq_df, test_freq_df)
  
# Combine annotations and results
anno_df <- formatted_df[mask_list, 1:10][pos_mask_list,]
output_df <- cbind(formatted_df[mask_list, 1:10][pos_mask_list,],
                     stats_df,
                     ctrl_ref_df[pos_mask_list,],
                     ctrl_alt_df[pos_mask_list,],
                     test_ref_df[pos_mask_list,],
                     test_alt_df[pos_mask_list,])
                     #formatted_df[mask_list, paste0(ctrl_name, "_fpkm"), drop = FALSE][pos_mask_list,, drop = FALSE],
                     #formatted_df[mask_list, paste0(test_name, "_fpkm"), drop = FALSE][pos_mask_list,, drop = FALSE])
output_df <- merge(output_df, as.data.frame(table(output_df$ensg_id)), by.x="ensg_id", by.y="Var1")
output_df$gene_num_events <- output_df$Freq
output_df <- output_df[, -which(names(output_df) == "Freq")]
  
# Save dataframe
write.csv(output_df,
            paste0("5_diff_", test_name, "_vs_", ctrl_name, ".csv"),
            row.names = FALSE
)