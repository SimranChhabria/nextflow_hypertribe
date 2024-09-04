# Load libraries ---------------------------
# Import libraries


suppressMessages(library(dplyr))
suppressMessages(library(biomaRt))
suppressMessages(library(stringr))

suppressMessages(library(argparse))


# Create argument parser
parser <- ArgumentParser()
parser$add_argument("--samples", help = "Path to the sample list CSV file")
parser$add_argument("--gtf_path", help = "Path to the GTF file")
parser$add_argument("--genome_build", default = "hg19", help = "Genome build")
parser$add_argument("--counts_path_list", nargs = "+", help = "Path to the counts file")

# Parse the arguments
args <- parser$parse_args()


genome <- args$genome_build
sample_list <- strsplit(args$samples, ",")[[1]]
print(sample_list)
gtf_path <- args$gtf_path
counts_path_list <- strsplit(args$counts_path_list, ",")[[1]]

# Load data ---------------------------
# Data directories
#counts_folder <- paste0("../../projects/", lab, "/", project, "/output_data/deseq2/3_counts/")
#formatted_output_path <- paste0("../../projects/", lab, "/", project, "/output_data/deseq2/4_formatted_output.csv")


# Filter data ---------------------------
#' Merge count files into dataframe
#'
#' This function reads count files and merge them into a dataframe.
#'
#' @param counts_folder The folder that stores all count files.
#' @param sample_list The sample names.
#'
#' @return The dataframe containing counts from all files.
merge_counts <- function(counts_path_list, sample_list) {

  # Get list of paths of count files
  #counts_path_list <- paste(counts_folder, sample_list, ".counts", sep = "")

  # Read list of count files
  counts_df_list <- lapply(
    counts_path_list,
    function(x) {
      read.table(file = x, sep = "\t")
    }
  )

  # Merge list of count tables
  counts_df <- suppressWarnings(Reduce(function(x, y) {
    merge(x, y, by = "V1")
  }, counts_df_list))
  counts_list <- paste0(sample_list, "_counts")
  names(counts_df) <- c("ensg_id", counts_list)
  
  # Ignore rows that do not represent genes
  counts_df <- counts_df %>%
    filter(!startsWith(ensg_id, "__")) %>%
    arrange(ensg_id)
  counts_df
}

# Merge count files
counts_df <- merge_counts(counts_path_list, sample_list)

# Get genomic data ---------------------------
#' Obtain genomic information for ENSG IDs
#'
#' This function queries the Ensembl database and returns the gene symbol,
#' chromosome name, strand, starting and ending positions of given ENSG IDs.
#'
#' @param ensg_id_list The list of ENSG IDs to be queried.
#' @param genome The genome symbol.
#'
#' @return A dataframe consisting of all genomic information of the ENSG IDs.
get_genomic_info <- function(ensg_id_list, genome) {

  # Determine genome database
  if (genome == "hg38") {
    biomart_dataset = "hsapiens_gene_ensembl"
    gene_symbol = "hgnc_symbol"
  }
  if (genome == "mm10") {
    biomart_dataset = "mmusculus_gene_ensembl"
    gene_symbol = "mgi_symbol"
  }
  
  # Load Ensembl database
  ensembl <- useEnsembl(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = biomart_dataset
  )

  # Get gene symbols for all ENSG IDs
  attribute_list <- c(
    "ensembl_gene_id", gene_symbol, "chromosome_name",
    "strand", "start_position", "end_position"
  )
  genomic_df <- getBM(
    attributes = attribute_list,
    filters = "ensembl_gene_id",
    values = ensg_id_list,
    mart = ensembl,
    useCache = FALSE
  )

  # Rename columns
  names(genomic_df) <- c(
    "ensg_id", "gene_symbol", "chr", "strand",
    "start", "end"
  )

  # Function to convert strand to UCSC format
  convert_strand <- function(strand) {
    if (strand == 1) {
      "+"
    } else {
      "-"
    }
  }

  # Remove entries with special chromosomes
  # Rename chromosomes in UCSC format
  # Rename strands in UCSC format
  # Compute gene length
  # Add feature
  # Add event_id
  # Add gene_num_events
  # Remove entries with duplicated ENSG IDs
  # Sort rows by ENSG IDs
  genomic_df <- genomic_df %>%
    filter(chr %in% c(1:22, "X", "Y")) %>%
    mutate(chr = paste0("chr", chr)) %>%
    mutate(strand = sapply(strand, convert_strand)) %>%
    mutate(length = end - start + 1) %>%
    mutate(feature = "gene") %>%
    mutate(event_id = paste0(ensg_id, "_1")) %>%
    mutate(gene_num_events = 1) %>%
    distinct(ensg_id, .keep_all = TRUE) %>%
    arrange(ensg_id)
}

genomic_df <- get_genomic_info(counts_df$ensg_id, genome)

# Format data ---------------------------
# Ordered list of samples
counts_list <- paste0(sample_list, "_counts")

# List of output column names
formatted_column_list <- c(
  "ensg_id", "gene_symbol", "chr", "strand",
  "start", "end", "length", "feature",
  "event_id", "gene_num_events", counts_list
)

# Merge counts with genomic information
formatted_df <- merge(counts_df, genomic_df, by = "ensg_id", all.x = TRUE) %>%
  dplyr::select(all_of(formatted_column_list))
formatted_df[!startsWith(formatted_df$ensg_id, "ENS"), "gene_symbol"] <- formatted_df[!startsWith(formatted_df$ensg_id, "ENS"), "ensg_id"]


# Save dataframe
write.csv(formatted_df,
  "4_formatted_output_counts.csv",
  row.names = FALSE
)

