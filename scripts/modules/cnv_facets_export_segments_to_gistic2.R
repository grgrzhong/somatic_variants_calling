#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-06-12
## Updated: 2025-06-12
## Description: Converts a segment TSV file (e.g., from FACETS VCF-to-TSV conversion)
##              to a custom TSV format with 
##              Sample Chromosome Start End Num_Probes Segment_Mean
##############################################################################

# Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(optparse)
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
# Note: Seg.CN = CNLR.MEDIAN - dipLogR.
# https://github.com/dariober/cnv_facets/issues/52
required_input_cols <- c("CHROM", "POS", "END", "NUM_MARK", "CNLR_MEDIAN")

input_tsv <- readr::read_tsv(
    opt$input, 
    show_col_types = FALSE
)

missing_cols <- required_input_cols[!required_input_cols %in% names(input_tsv)]
    
if (length(missing_cols) > 0) {
    stop(paste("Input TSV file is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

## Convert to custom format, filtering out neutral segments
output_tsv <- input_tsv |> 
            filter(SVTYPE != "NEUTR") |> 
            ## The not estimateable Lesser (minor) copy number segments were filtered out
            dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |>
            dplyr::select(all_of(required_input_cols)) |> 
            dplyr::rename(
                Chromosome = CHROM,
                Start = POS,
                End = END,
                Num_Probes = NUM_MARK,
                Segment_Mean = CNLR_MEDIAN
            ) |> 
            mutate(
                Chromosome=str_replace(Chromosome, "chr", ""),
                Start = as.integer(Start),
                End = as.integer(End)
            ) |> 
            ## Make sure the the segments to be start < end
            mutate(
                Start_fixed = if_else(Start > End, End, Start),
                End_fixed = if_else(Start > End, Start, End)
            ) |> 
            dplyr::select(
                Chromosome, Start = Start_fixed, End = End_fixed, 
                Num_Probes, Segment_Mean
            ) |> 
            ## Keep data types consistent
            mutate(
                Chromosome = as.character(Chromosome),
                Start = as.integer(Start),
                End = as.integer(End),
                Num_Probes = as.integer(Num_Probes),
                Segment_Mean = as.numeric(Segment_Mean)
            )

# Write output to TSV file    
readr::write_tsv(
    output_tsv, 
    file = opt$output, 
    na = "NA",
    quote = "none"
)
