#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-06-12
## Description: Converts a segment TSV file (e.g., from FACETS VCF-to-TSV conversion)
##              to a custom TSV format
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
        help = "Output custom .seg file path for converting to bedtools format", 
        metavar = "FILE"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    stop("Input file or directory must be specified with -i or --input")
}

if (is.null(opt$output)) {
    stop("Output file path must be specified with -o or --output")
}

if (!file.exists(opt$input)) {
    stop(paste("Input file not found:", opt$input), call. = FALSE)
}

input_tsv <- readr::read_tsv(
    opt$input, 
    show_col_types = FALSE
)

## The not estimateable Lesser (minor) copy number segments were filtered out
output_tsv <- input_tsv |> 
            dplyr::filter(!(LCN_EM %in% c(".", NA, "NA"))) |> 
            dplyr::filter(SVTYPE != "NEUTR") |> 
            ## 0-based for bed format
            mutate(POS = POS - 1)

## Write output to .seg file    
readr::write_tsv(
    output_tsv, 
    file = opt$output, 
    na = "NA"
)