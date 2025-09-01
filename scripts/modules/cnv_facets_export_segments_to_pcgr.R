#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-06-12
## Updated: 2025-06-12
## Description: Converts a segment TSV file for PCGR input.
## to a custom TSV format with Chromosome, Start, End, nMajor, nMinor.
##############################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
        library(fs)
        library(tidyverse)
    })
)

# Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character", 
        default = NULL,
        help = "Input segment TSV file path (e.g., from FACETS VCF-to-TSV)", 
        metavar = "FILE"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", 
        default = NULL,
        help = "Output custom TSV file path (columns: Chromosome, Start, End, nMajor, nMinor)", 
        metavar = "FILE"
    )
)

# Create an option parser
opt_parser <- OptionParser(
    option_list = option_list, 
    description = "Converts a segment TSV file to a custom format with nMajor/nMinor."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input TSV file path must be specified with -i or --input", call. = FALSE)
}

if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output TSV file path must be specified with -o or --output", call. = FALSE)
}

if (!file.exists(opt$input)) {
    stop(paste("Input file not found:", opt$input), call. = FALSE)
}

# Read input, treating "." as NA for relevant columns
# https://github.com/dariober/cnv_facets
required_input_cols <- c("CHROM", "POS", "END", "TCN_EM", "LCN_EM")

input_tsv <- readr::read_tsv(
    opt$input, 
    show_col_types = FALSE
)

missing_cols <- required_input_cols[!required_input_cols %in% names(input_tsv)]
    
if (length(missing_cols) > 0) {
    stop(paste("Input TSV file is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

output_tsv <- input_tsv |> 
            dplyr::select(all_of(required_input_cols)) |> 
            ## The not estimateable Lesser (minor) copy number segments were filtered out
            dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |>
            dplyr::rename(
                Chromosome = CHROM,
                Start = POS,
                End = END
            ) |> 
            mutate(
                Chromosome=str_replace(Chromosome, "chr", ""),
                Start = as.integer(Start),
                End = as.integer(End),
                TCN_EM = as.integer(TCN_EM),
                LCN_EM = as.integer(LCN_EM)
            ) |> 
            mutate(
                nMajor = as.integer(round(TCN_EM - LCN_EM)),
                nMinor = as.integer(round(LCN_EM))
            ) |>
            dplyr::select(
                Chromosome, Start, End, nMajor, nMinor
            ) |> 
            ## Make sure the the segments to be start > end
            mutate(
                Start_fixed = if_else(Start > End, End, Start),
                End_fixed = if_else(Start > End, Start, End)
            ) |> 
            dplyr::select(
                Chromosome, Start = Start_fixed, End = End_fixed, nMajor, nMinor
            )

# Write output to TSV file    
readr::write_tsv(
    output_tsv, 
    file = opt$output, 
    na = "NA"
)
