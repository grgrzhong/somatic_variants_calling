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
                chrom = CHROM,
                segment_start = POS,
                segment_end = END,
                num_mark = NUM_MARK,
                log2ratio = CNLR_MEDIAN
            ) |> 
            mutate(
                chrom=str_replace(chrom, "chr", ""),
                segment_start = as.integer(segment_start),
                segment_end = as.integer(segment_end)
            ) |> 
            ## Make sure the the segments to be start < end
            mutate(
                segment_start_fixed = if_else(
                    segment_start > segment_end, segment_end, segment_start
                ),
                segment_end_fixed = if_else(
                    segment_start > segment_end, segment_start, segment_end
                )
            ) |> 
            dplyr::select(
                chrom, segment_start = segment_start_fixed, segment_end = segment_end_fixed, 
                num_mark, log2ratio
            ) |> 
            ## Keep data types consistent
            mutate(
                chrom = as.character(chrom),
                segment_start = as.integer(segment_start),
                segment_end = as.integer(segment_end),
                num_mark = as.integer(num_mark),
                log2ratio = as.numeric(log2ratio)
            )

# Write output to TSV file    
readr::write_tsv(
    output_tsv, 
    file = opt$output, 
    na = "NA",
    quote = "none"
)

output_tsv <- output_tsv |> 
    mutate(sample = "DFSP-010-T-P1", .before = chrom)
