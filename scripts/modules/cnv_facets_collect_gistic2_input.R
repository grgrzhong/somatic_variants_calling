#!/usr/bin/env Rscript
##############################################################################
## Author:  Zhong Guorui
## Created: 2025-07-17
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
        help = "Input directory containing *.gistic2.tsv files", 
        metavar = "character"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", 
        default = NULL,
        help = "Output TSV file path for combined data", 
        metavar = "character"
    ),
    make_option(
        c("-p", "--pattern"),
        type = "character", 
        default = "\\.gistic2$",
        help = "Pattern to match GISTIC2 files (default: \".gistic2\")", 
        metavar = "character"
    )
)

## Create an option parser
opt_parser <- OptionParser(
    option_list = option_list,
    description = "Collects and combines GISTIC2 segment files from a directory."
)

opt <- parse_args(opt_parser)

## Check if input directory is provided
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input directory must be specified with -i or --input")
}

## Check if output file is provided
if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output file path must be specified with -o or --output")
}

## Check if input directory exists
if (!dir.exists(opt$input)) {
    stop(paste("Input directory not found:", opt$input), call. = FALSE)
}

## Fina all GISTIC2 segment files in the input directory
files <- dir_ls(
    opt$input,
    glob = "*.gistic2.tsv",
    recurse = TRUE
)

## Load and combine all GISTIC2 segment files
gistic_list <- list()
for (file in files) {
    
    ## Get the sample id from the file name
    sample_id <- path_file(file) |> 
        str_remove("\\.gistic2\\.tsv$")
    
    message(paste0(" - Processing ", sample_id, " = ", file))

    gistic_data <- read_tsv(
        file,
        show_col_types = FALSE,
        col_types = cols(
            Chromosome = col_character(),
            Start = col_integer(),
            End = col_integer(),
            Num_Probes = col_integer(),
            Segment_Mean = col_double()
        ) 
    ) |> 
    mutate(Sample = sample_id, .before = Chromosome)

    gistic_list[[sample_id]] <- gistic_data

}

## Combine all data frames into one
gistic_tbl <- list_rbind(gistic_list)

## Save cohort level segment data for GISTIC2 input
write.table(
    gistic_tbl,
    file = opt$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

## example
# opt$input <- "/mnt/f/projects/sarcoma_multiomics/data/wes/CNV/cnv_facets"
# opt$output <- "/mnt/f/projects/sarcoma_multiomics/data/wes/Processed/cnv_facets_DFSP_cohort_gistic2.tsv"
