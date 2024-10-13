# Load libraries ---------------------------
# Import libraries
# Load libraries ---------------------------
# Import libraries

# if (!require("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")}

# if (!requireNamespace("biomaRt", quietly = TRUE)) {
#   BiocManager::install("biomaRt")
# }

# if (!requireNamespace("BiocParallel", quietly = TRUE)) {
#   BiocManager::install("BiocParallel")
# }
# if (!requireNamespace("DESeq2", quietly = TRUE)) {
#   BiocManager::install("DESeq2")
# }
# if (!requireNamespace("dplyr", quietly = TRUE)) {
#   install.packages("dplyr")
# }
# if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
#   BiocManager::install("GenomicAlignments")
# }
# if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
#   BiocManager::install("GenomicFeatures")
# }
# if (!requireNamespace("parallel", quietly = TRUE)) {
#   install.packages("parallel")
# }
# if (!requireNamespace("reshape2", quietly = TRUE)) {
#   install.packages("reshape2")
# }
# if (!requireNamespace("Rsamtools", quietly = TRUE)) {
#   BiocManager::install("Rsamtools")
# }
# if (!requireNamespace("stringr", quietly = TRUE)) {
#   install.packages("stringr")
# }
# if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
#   BiocManager::install("VariantAnnotation")
# }
# if (!requireNamespace("argparse", quietly = TRUE)) {
#   install.packages("argparse")
# }
suppressMessages(library(biomaRt))
suppressMessages(library(BiocParallel))
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(parallel))
suppressMessages(library(reshape2))
suppressMessages(library(Rsamtools))
suppressMessages(library(stringr))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(argparse))
# Load data ---------------------------
# Data directories


# Create argument parser
parser <- ArgumentParser()
parser$add_argument("--samples", help = "Path to the sample list CSV file")
parser$add_argument("--gtf_path", help = "Path to the GTF file")
parser$add_argument("--bam_path_list", nargs = "+", help = "List of BAM file paths")
parser$add_argument("--vcf_path_list", nargs = "+", help = "List of VCF file paths")
parser$add_argument("--genome_build", default = "hg19", help = "Genome build")
parser$add_argument("--editor", default = "ADAR", help = "Editor (ADAR/APOBEC)")
parser$add_argument("--working_dir", default = getwd(), help = "Current working directory")
#add resume flag
parser$add_argument("--resume", default = FALSE, action = "store_true", help = "Resume from existing files")

# Parse the arguments
args <- parser$parse_args()

working_dir <- args$working_dir
# Set the working directory
setwd(working_dir)

sample_list <- strsplit(args$samples, ",")[[1]]
gtf_path <- args$gtf_path
bam_path_list <- strsplit(args$bam_path_list, ",")[[1]]
cat('bam_path_list: ', bam_path_list)
vcf_path_list <- strsplit(args$vcf_path_list, ",")[[1]]
cat('vcf_path_list: ', bam_path_list)
genome_build <- args$genome_build
editor  <- args$editor

# genome_build <- "hg38"
# genome_dir <- paste0("/data/kharas/baalii/nextflow_pipeline/genome/", genome_build)
# gtf_path <- paste0(genome_dir,'/', genome_build, '.gtf')
# 
# editor="ADAR"
# 
# sample_list <- "MSI2-ADAR-3,MSI2-ADAR-2,MSI2-ADAR-1,MIG-2,MSI1-ADAR-3,MSI1-ADAR-2,MSI1-ADAR-1,MIG-3,MIG-1"
# sample_list <- strsplit(sample_list, ",")[[1]][c(3, 7, 9)]
# 
# vcf_path_list <- "/lila/data/kharas/baalii/nextflow_pipeline/work/b5/2e36b708a3d0a6ca2cef083404b7a9/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/4f/ff887ab100dff906c5a765fbfb5d8c/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/1f/b01a472dd524dff28552f03ae9c43f/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/03/7830ed538e1d37f179725edc4c8312/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/19/263c9a2776f0c9615fb50f468d620e/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/43/6d74e40404bf4df8f5813613e2163e/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/ad/8cc380e5f36aae085fb3c64fb1b26e/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/79/b5c5a1aed5539b43d9a0ff2505bc0e/final.vcf, /lila/data/kharas/baalii/nextflow_pipeline/work/56/e1069fe517acf7b58b0d3668f27265/final.vcf"
# vcf_path_list <- gsub(" ", "", vcf_path_list)
# vcf_path_list <- gsub("/lila", "", vcf_path_list)
# vcf_path_list <- strsplit(vcf_path_list, ",")[[1]][c(3, 7, 9)]
# 
# bam_path_list <- "/lila/data/kharas/baalii/nextflow_pipeline/work/7e/b29aa9bf482d83cde0154a66640adc/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/93/f6c5d3f0ec240049a7dc13197a56e1/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/5c/599a4b1cb3b3d49a62ea83450faefc/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/4d/1fabd0932f8160ae532357ee3e0944/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/37/e4171544692e1c966da9d186a98cb5/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/42/7b6c4647b97f5d7f9c3600f2cca528/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/7c/d7f93a5a48f6a28f0ed1f7db47e726/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/63/89478d230af76b8070120ad243ee51/recalibrated.bam, /lila/data/kharas/baalii/nextflow_pipeline/work/87/2c3985d01eb7a32a14c460ab98b04c/recalibrated.bam"
# bam_path_list <- gsub(" ", "", bam_path_list)
# bam_path_list <- gsub("/lila", "", bam_path_list)
# bam_path_list <- strsplit(bam_path_list, ",")[[1]][c(3, 7, 9)]


names(bam_path_list) <- sample_list
names(vcf_path_list) <- sample_list


# Read transcripts from GTF file
TxDb <- makeTxDbFromGFF(gtf_path, format = "gtf")
txbygene_grl <- transcriptsBy(TxDb, "gene")




# Get read counts per transcript ---------------------------
read_counts_gr <- summarizeOverlaps(
  features = txbygene_grl,
  reads = BamFileList(bam_path_list,
    yieldSize = 1e6
  ),
  mode = "Union",
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE,
  BPPARAM = SnowParam(workers = 3)
)
# saveRDS(read_counts_gr, "4a_read_counts_gr.rds")

#read_counts_gr <- readRDS("/data/morrisq/simranch/nextflow/scripts/base_editing/4a_read_counts_gr.rds")


# Filter VCF files ---------------------------
#' Read a VCF file and filter entries
#'
#' This function reads a VCF file, and determines the strand that the edit is
#' on. It also filters out edits that are potential SNPs and those that are
#' not present within annotated transcripts.
#'
#' @param vcf_path The path to the VCF file.
#' @param genome The genome symbol.
#' @param editor The editor (ADAR/APOBEC).
#'
#' @return The filtered VCF table as a GenomicRanges object.
filter_vcf <- function(vcf_path, genome, editor) {

  cat ("Reading VCF file: ", vcf_path, "\n")
  # Read VCF file
  vcf <- readVcf(vcf_path, genome)

  cat("Number of edits in VCF file:", length(vcf), "\n")

  # cat (info(vcf)$DB, "\n")
  # Filter out edits that are potential SNPs
  filtered_vcf <- vcf
  # [!info(vcf)$DB]
  
  # Determine reference and alternative allele
  if (editor == "ADAR") {
    ref_pos <- "A"
    alt_pos <- "G"
    ref_neg <- "T"
    alt_neg <- "C"
  }
  if (editor == "APOBEC") {
    ref_pos <- "C"
    alt_pos <- "T"
    ref_neg <- "G"
    alt_neg <- "A"
  }

  # Determine the strand that the edit is on
  pos_mask <- rowRanges(filtered_vcf)$REF == ref_pos &
    sapply(rowRanges(filtered_vcf)$ALT, function(seq) {
      return(alt_pos %in% seq)
    })
  pos_gr <- rowRanges(filtered_vcf)[pos_mask]
  strand(pos_gr) <- "+"
  neg_mask <- rowRanges(filtered_vcf)$REF == ref_neg &
    sapply(rowRanges(filtered_vcf)$ALT, function(seq) {
      return(alt_neg %in% seq)
    })
  neg_gr <- rowRanges(filtered_vcf)[neg_mask]
  strand(neg_gr) <- "-"
  filtered_vcf_gr <- c(pos_gr, neg_gr)


  # Only keep edits found in transcripts
  filtered_vcf_gr <- filtered_vcf_gr[countOverlaps(
    filtered_vcf_gr,
    txbygene_grl
  ) > 0]

  cat("Number of edits after filtering for transcripts:", length(filtered_vcf_gr), "\n")
  
  return(filtered_vcf_gr)  # Add this line to return the filtered_vcf_gr object
}

names(vcf_path_list) <- sample_list
# Read and filter VCF files
filtered_vcf_gr_list <- lapply(vcf_path_list, function(vcf_path) {
  tryCatch(
    filter_vcf(vcf_path, genome_build, editor),
    error = function(e) {
      message(paste("Error filtering VCF:", vcf_path, e$message))
      NULL
    }
  )
})

filtered_vcf_gr_list <- filtered_vcf_gr_list[!sapply(filtered_vcf_gr_list, is.null)]
names(filtered_vcf_gr_list) <- sample_list
filtered_vcf_grl <- GRangesList(filtered_vcf_gr_list)

# # Save VCF dataframe
# saveRDS(filtered_vcf_grl, "4b_filtered_vcf_grl.rds")

#filtered_vcf_grl <- readRDS("/data/morrisq/simranch/nextflow/scripts/base_editing/4b_filtered_vcf_grl.rds")

# Get all edit sites that only match to one transcript
unique_vcf_gr <- sort(unique(unlist(filtered_vcf_grl, use.names = FALSE)))
subset_vcf_gr <- unique_vcf_gr[which(countOverlaps(unique_vcf_gr, txbygene_grl) == 1)]

# Generate pileup dataframes ---------------------------
#' Generate pileup dataframe for all edit sites for a sample
#'
#' This function reads the BAM file for a sample and generate a pileup
#' dataframe for the union of edit sites of all samples.
#'
#' @param bam_path The path to the BAM file.
#' @param unique_vcf_gr The GenomicRanges object describing the location
#'   of the union of edit sites of all samples.
#'
#' @return The pileup dataframe.
generate_pileup_df <- function(bam_path, unique_vcf_gr) {
  
  # Read BAM file
  bai_path <- paste0(bam_path, ".bai")
  bam <- BamFile(bam_path, bai_path)

  # Generate pileup dataframe
  bam_param <- ScanBamParam(
    flag = scanBamFlag(
      hasUnmappedMate = FALSE,
      isProperPair = TRUE,
      isDuplicate = FALSE
    ),
    which = unique_vcf_gr
  )
  pileup_param <- PileupParam(
    distinguish_strands = FALSE,
    min_base_quality = 10,
    max_depth = 1e4
  )
  pileup_df <- pileup(bam, unique_vcf_gr,
    scanBamParam = bam_param, pileupParam = pileup_param
  )
}
 cat ("Generating pileup dataframes for all samples\n")
# Generate pileup dataframes for all samples
pileup_df_list <- mclapply(bam_path_list, generate_pileup_df, unique_vcf_gr,
                           mc.cores = length(bam_path_list))
names(pileup_df_list) <- sample_list

# Save list of pileup dataframes
# saveRDS(pileup_df_list, "4c_pileup_df_list.rds")

#pileup_df_list <- readRDS("/data/morrisq/simranch/nextflow/scripts/base_editing/4c_pileup_df_list.rds")


cat ("Pileup dataframes saved\n")

# Process pileup dataframes ---------------------------
#' Process pileup dataframe for a sample
#'
#' This function adds strand information to and sort a pileup dataframe
#'
#' @param pileup_df The pileup dataframe.
#' @param unique_vcf_gr The GenomicRanges object describing the location
#'   of the union of edit sites of all samples.
#' @param editor The editor (ADAR/APOBEC).
#'
#' @return The processed pileup dataframe.
process_pileup_df <- function(pileup_df, unique_vcf_gr, editor) {

  # Melt table by nucleotide coordinates
  pileup_df <- dcast(pileup_df, which_label ~ nucleotide,
    value.var = "count",
    fill = 0, drop = FALSE
  )

  # Add strand information
  pileup_df$strand <- as.character(strand(unique_vcf_gr))

  # Check if there are missing rows
  snp_id_list <- sprintf(
    "%s:%d-%d",
    seqnames(unique_vcf_gr),
    start(unique_vcf_gr),
    start(unique_vcf_gr)
  )
  stopifnot(all(snp_id_list == pileup_df$which_label))
  row.names(pileup_df) <- pileup_df$which_label
  
  # Determine reference and alternative allele
  if (editor == "ADAR") {
    ref_pos <- "A"
    alt_pos <- "G"
    ref_neg <- "T"
    alt_neg <- "C"
  }
  if (editor == "APOBEC") {
    ref_pos <- "C"
    alt_pos <- "T"
    ref_neg <- "G"
    alt_neg <- "A"
  }

  # Split and reorganize dataframe based on strand
  pileup_df <- split(pileup_df, pileup_df$strand)
  pileup_df$`+` <- data.frame(
    ref_count = pileup_df$`+`[ref_pos],
    alt_count = pileup_df$`+`[alt_pos],
    row.names = row.names(pileup_df$`+`)
  )
  names(pileup_df$`+`) <- c("ref", "alt")
  pileup_df$`-` <- data.frame(
    ref_count = pileup_df$`-`[ref_neg],
    alt_count = pileup_df$`-`[alt_neg],
    row.names = row.names(pileup_df$`-`)
  )
  names(pileup_df$`-`) <- c("ref", "alt")

  # Combine and sort data frames
  rbind(pileup_df$`+`, pileup_df$`-`)[snp_id_list, ]
}

cat ("Processing pileup dataframes for all samples\n")

# Process pileup dataframe for all samples
processed_pileup_df_list <- lapply(
  pileup_df_list, process_pileup_df,
  unique_vcf_gr,
  editor
)

# # Save list of processed pileup dataframes
# saveRDS(
#   processed_pileup_df_list,
#   "4d_processed_pileup_df_list.rds"
# )

#processed_pileup_df_list <- readRDS("/data/morrisq/simranch/nextflow/scripts/base_editing/4d_processed_pileup_df_list.rds")

cat ("Processed pileup dataframes saved\n")

# Get genomic data ---------------------------
#' Obtain genomic information for all edit sites
#'
#' This function queries the Ensembl database and returns the gene symbol,
#' for the genes of all edit sites. It also queries for the genomic feature
#' where the edit site is at.
#'
#' @param subset_vcf_gr All edit sites that are present in only one transcript.
#' @param txbygene_grl The list of all genes and their transcripts.
#' @param gtf_path The path to the GTF file.
#'
#' @return A dataframe consisting of all genomic information of the edits.
get_anno_df <- function(subset_vcf_gr, txbygene_grl, gtf_path, genome_build) {

  # Get matches between edits and transcripts
  hits <- findOverlaps(subset_vcf_gr, txbygene_grl)

  # Get names of transcripts
  ensg_id_list <- names(txbygene_grl)[subjectHits(hits)]
  
  # Determine genome database
  if (genome_build == "hg38") {
    biomart_dataset = "hsapiens_gene_ensembl"
    gene_symbol = "hgnc_symbol"
  }
  if (genome_build == "mm10") {
    biomart_dataset = "mmusculus_gene_ensembl"
    gene_symbol = "mgi_symbol"
  }

  # Load Ensembl database
  ensembl <- useEnsembl(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = biomart_dataset,
    mirror = "uswest"
  )

  # Get gene symbols for all ENSG IDs
  genomic_df <- getBM(
    attributes = c("ensembl_gene_id", gene_symbol),
    filters = "ensembl_gene_id",
    values = ensg_id_list,
    mart = ensembl,
    useCache = FALSE
  )

  # Remove entries with duplicated ENSG IDs
  genomic_df <- genomic_df %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

  # Generate annotation dataframe
  anno_df <- DataFrame(ensg_id = ensg_id_list)
  anno_df <- merge(anno_df, genomic_df,
    by.x = "ensg_id", by.y = "ensembl_gene_id",
    all.x = TRUE, sort = FALSE
  )
  anno_df <- anno_df[match(ensg_id_list, anno_df$ensg_id), ]
  names(anno_df) <- c("ensg_id", "gene_symbol")

  # Determine whether edit is in intro, exon, 5'UTR or 3'UTR
  gtf <- read.table(gtf_path, header = FALSE, sep = "\t")

  exon_gtf <- gtf[gtf$V3 == "exon", c(1, 4, 5, 7)]
  names(exon_gtf) <- c("seqname", "start", "end", "strand")
  exon_gr <- makeGRangesFromDataFrame(exon_gtf)

  utr5_gtf <- gtf[gtf$V3 == "5UTR", c(1, 4, 5, 7)]
  names(utr5_gtf) <- c("seqname", "start", "end", "strand")
  utr5_gr <- makeGRangesFromDataFrame(utr5_gtf)

  utr3_gtf <- gtf[gtf$V3 == "3UTR", c(1, 4, 5, 7)]
  names(utr3_gtf) <- c("seqname", "start", "end", "strand")
  utr3_gr <- makeGRangesFromDataFrame(utr3_gtf)

  anno_df$feature <- "intron"
  anno_df$feature[countOverlaps(subset_vcf_gr, exon_gr) > 0] <- "exon"
  anno_df$feature[countOverlaps(subset_vcf_gr, utr5_gr) > 0] <- "utr5"
  anno_df$feature[countOverlaps(subset_vcf_gr, utr3_gr) > 0] <- "utr3"
  anno_df
}

anno_df <- get_anno_df(subset_vcf_gr, txbygene_grl, gtf_path, genome_build)
cat ("Genomic information obtained\n")

# Combine data ---------------------------
#' Adds counts to all edit sites
#'
#' This function adds counts to the locations of all edit sites.
#'
#' @param pileup_df The pileup dataframe.
#' @param subset_vcf_gr All edit sites that are present in only one transcript.
#'
#' @return A dataframe consisting of both the counts and location of edit sites.
add_counts <- function(pileup_df, subset_vcf_gr) {
  snp_id_list <- sprintf(
    "%s:%d-%d",
    seqnames(subset_vcf_gr),
    start(subset_vcf_gr),
    start(subset_vcf_gr)
  )
  pileup_df <- pileup_df[snp_id_list, ]
  formatted_df <- cbind(as.data.frame(subset_vcf_gr), pileup_df)
  formatted_df
}

# Add genomic information
mcols(subset_vcf_gr) <- anno_df

# Add counts to all pileup dataframes
formatted_df <- lapply(processed_pileup_df_list, add_counts, subset_vcf_gr)

# Get genomic information
formatted_anno_df <- formatted_df[[1]][(1:(ncol(formatted_df[[1]]) - 2))]

# Get counts
formatted_counts_df <- do.call(
  "cbind",
  lapply(
    formatted_df,
    function(df) {
      df[-(1:(ncol(formatted_df[[1]]) - 2))]
    }
  )
)


# Combine annotations and count dataframes
formatted_df <- cbind(formatted_anno_df, formatted_counts_df)

# Format data ---------------------------
# Add number of edits per gene
event_id_list <- row.names(formatted_df)
gene_num_edits_df <- as.data.frame(table(formatted_df$ensg_id))
names(gene_num_edits_df) <- c("ensg_id", "gene_num_events")
formatted_df <- merge(formatted_df, gene_num_edits_df,
  by = "ensg_id",
  keep.x = TRUE, keep.y = FALSE, sort = FALSE
)
formatted_df$event_id <- event_id_list


# Add FPKM
colData(read_counts_gr)$condition <- factor(sample_list)
fpkm_df <- fpkm(DESeqDataSet(read_counts_gr, ~condition))
prefix_list <- unique(sub("\\_[0-9]", "", colnames(fpkm_df)))
fpkm_df <- as.data.frame(sapply(prefix_list, function(x) rowMeans(fpkm_df[, startsWith(colnames(fpkm_df), x), drop = FALSE])))
colnames(fpkm_df) <- paste0(colnames(fpkm_df), "_fpkm")
fpkm_df$ensg_id <- row.names(fpkm_df)
formatted_df <- merge(formatted_df, fpkm_df, by = "ensg_id")

# Rename columns
renamed_col_list <- names(formatted_df)
renamed_col_list <- sub("seqnames", "chr", renamed_col_list)
renamed_col_list <- sub("width", "length", renamed_col_list)
renamed_col_list <- sub(".ref", "_ref_counts", renamed_col_list)
renamed_col_list <- sub(".alt", "_alt_counts", renamed_col_list)
names(formatted_df) <- renamed_col_list

# Order columns and rows
ordered_col_list <- c(
  "ensg_id", "gene_symbol", "chr", "strand", "start",
  "end", "length", "feature", "event_id", "gene_num_events",
  renamed_col_list[endsWith(renamed_col_list, "_counts")],
  renamed_col_list[endsWith(renamed_col_list, "_fpkm")])
formatted_df <- formatted_df %>%
  dplyr::select(all_of(ordered_col_list)) %>%
  arrange(ensg_id)

# Save dataframe
write.csv(formatted_df,
  "4_formatted_output.csv",
  row.names = FALSE
)

cat ("Formatted output saved\n")
