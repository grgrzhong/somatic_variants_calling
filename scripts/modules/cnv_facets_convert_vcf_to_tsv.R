#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-05-20
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
    make_option(c("-i", "--input"),
        type = "character", 
        default = NULL,
        help = "Input VCF file path", 
        metavar = "character"
    ),
    make_option(
        c("-o", "--output"),
        type = "character", 
        default = NULL,
        help = "Output TSV file path (required)", 
        metavar = "character"
    )
)

# Create an option parser
opt_parser <- OptionParser(
    option_list = option_list,
    description = "Converts a VCF file to a custom TSV format with required columns."
)

opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    stop("Input VCF file must be specified with -i or --input")
}

if (is.null(opt$output)) {
    stop("Output TSV file path must be specified with -o or --output")
}

if (!file.exists(opt$input)) {
    stop(paste("Input VCF file not found:", opt$input), call. = FALSE)
}

# Read VCF file
con <- if (grepl("\\.gz$", opt$input)) {
    gzfile(opt$input, "r")
} else {
    file(opt$input, "r")
}

lines <- readLines(con)
close(con)

# Extract sample-level metadata from headers (purity, ploidy, etc.)
metadata <- list()

header_lines <- lines[grepl("^##", lines)]

suppressWarnings(

    for (line in header_lines) {
        if (grepl("^##purity=", line)) {

            metadata[["purity"]] <- as.numeric(sub("^##purity=", "", line))

        } else if (grepl("^##ploidy=", line)) {

            metadata[["ploidy"]] <- as.numeric(sub("^##ploidy=", "", line))

        } else if (grepl("^##dipLogR=", line)) {
            
            metadata[["dipLogR"]] <- as.numeric(sub("^##dipLogR=", "", line))

        } else if (grepl("^##est_insert_size=", line)) {
            
            metadata[["est_insert_size"]] <- as.numeric(
                sub("^##est_insert_size=", "", line)
            )

        } else if (grepl("^##emflags=", line)) {
            
            metadata[["emflags"]] <- sub("^##emflags=", "", line)
        }
    }
)

# Extract INFO field descriptions to understand data types
info_fields <- list()
for (line in header_lines) {
    if (grepl("^##INFO=<ID=", line)) {
        # Extract field name
        id <- sub("^##INFO=<ID=([^,]+),.+", "\\1", line)
        info_fields[[id]] <- line
    }
}

# Filter out header lines and parse remaining lines
data_lines <- lines[!grepl("^#", lines)]

# Check if we have data
if (length(data_lines) == 0) {
    message("No data found in ", opt$input)
    return(NULL)
}

# Parse VCF data lines
records <- list()

for (line in data_lines) {

    fields <- strsplit(line, "\t")[[1]]

    # Extract basic fields (CHROM, POS)
    CHROM <- fields[1]
    POS <- as.integer(fields[2])

    # Process INFO field - careful parsing to handle all field types correctly
    info_dict <- list()
    info_fields <- strsplit(fields[8], ";")[[1]]

    for (info in info_fields) {
    
        if (grepl("=", info)) {
    
            key_value <- strsplit(info, "=", fixed = TRUE)[[1]]
    
            if (length(key_value) == 2) {
    
                # Try to convert to numeric if appropriate
                value <- key_value[2]
    
                if (grepl("^-?[0-9.]+$", value)) {
                    value <- as.numeric(value)
                }
    
                info_dict[[key_value[1]]] <- value
            }

        } else {
            
            info_dict[[info]] <- TRUE
        }
    }

    # Required fields for SEG format
    if (!all(c("END", "SVTYPE") %in% names(info_dict))) {
        next
    }

    END <- as.integer(info_dict[["END"]])

    SVTYPE <- info_dict[["SVTYPE"]]

    # Extract NUM_MARK (Num_Markers)
    NUM_MARK <- if ("NUM_MARK" %in% names(info_dict)) {
        as.integer(info_dict[["NUM_MARK"]])
    } else {
        0
    }

    # Extract CNLR_MEDIAN (Seg.CN)
    CNLR_MEDIAN <- if ("CNLR_MEDIAN" %in% names(info_dict)) {
        as.numeric(info_dict[["CNLR_MEDIAN"]])
    } else {
        0
    }

    # Create record with the required fields in the specified order
    record <- list(
        "CHROM" = CHROM,
        "POS" = POS,
        "END" = END,
        "NUM_MARK" = NUM_MARK,
        "CNLR_MEDIAN" = CNLR_MEDIAN,
        "SVTYPE" = SVTYPE
    )

    # Add all other INFO fields
    for (key in names(info_dict)) {
        if (!key %in% c("END", "SVTYPE", "NUM_MARK", "CNLR_MEDIAN")) {
            record[[key]] <- info_dict[[key]]
        }
    }

    # Add sample-level metadata
    for (key in names(metadata)) {
        if (!key %in% names(record)) {
            record[[key]] <- metadata[[key]]
        }
    }

    records[[length(records) + 1]] <- record
}

# Check if we have records
if (length(records) == 0) {
    message("No valid records found in ", opt$input)
    return(NULL)
}

# Get all keys from the records
all_keys <- unique(unlist(lapply(records, names)))

# Define fixed columns - these will always be present in the output in this exact order
fixed_cols <- c(
    # Required SEG format columns
    "CHROM", "POS", "END", "NUM_MARK", "CNLR_MEDIAN", "SVTYPE",

    # Important sample metadata
    "purity", "ploidy", "dipLogR", "est_insert_size", "emflags",

    # Other CNV facet info tags
    "SVLEN", "NHET", "CNLR_MEDIAN", "CNLR_MEDIAN_CLUST", "MAF_R",
    "MAF_R_CLUST", "SEGCLUST", "CF_EM", "TCN_EM", "LCN_EM", "CNV_ANN"
)

# Get any remaining columns not in the fixed list
remaining_cols <- setdiff(all_keys, fixed_cols)

# Final column order: fixed columns first, then any remaining columns
col_order <- c(fixed_cols, remaining_cols)

# Create a data frame with all needed columns
df_cols <- union(fixed_cols, all_keys)
df <- data.frame(matrix(NA, nrow = length(records), ncol = length(df_cols)))
names(df) <- df_cols

# Fill data frame with records
for (i in seq_along(records)) {
    for (key in names(records[[i]])) {
        df[i, key] <- records[[i]][[key]]
    }
}

# Get the final columns that actually exist in our data
final_cols <- intersect(col_order, names(df))

# Apply column ordering
df <- df[, final_cols]

# Write to output TSV file
write.table(
    df,
    file = opt$output,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

