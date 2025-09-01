#!/usr/bin/env Rscript

##############################################################################
## Author:  Zhong Guorui
## Created: 2025-05-20
##############################################################################

# Load required libraries
library(optparse)

# Parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Input VCF file or directory containing VCF files", metavar = "character"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file path (required)", metavar = "character"
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

# Function to parse a VCF file and convert to SEG format
process_vcf_to_seg <- function(vcf_file, sample_name = NULL) {
    # Extract sample name from filename if not provided
    if (is.null(sample_name)) {
        sample_name <- basename(vcf_file)
        sample_name <- gsub("\\.vcf\\.gz$", "", sample_name)
    }

    # Read VCF file
    con <- if (grepl("\\.gz$", vcf_file)) {
        gzfile(vcf_file, "r")
    } else {
        file(vcf_file, "r")
    }

    lines <- readLines(con)
    close(con)

    # Extract sample-level metadata from headers (purity, ploidy, etc.)
    metadata <- list()
    header_lines <- lines[grepl("^##", lines)]

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
        message("No data found in ", vcf_file)
        return(NULL)
    }

    # Parse VCF data lines
    records <- list()

    for (line in data_lines) {
        fields <- strsplit(line, "\t")[[1]]

        # Extract basic fields (CHROM, POS)
        chrom <- fields[1]
        pos <- as.integer(fields[2])

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

        end <- as.integer(info_dict[["END"]])
        svtype <- info_dict[["SVTYPE"]]

        # Extract NUM_MARK (Num_Markers)
        num_mark <- if ("NUM_MARK" %in% names(info_dict)) {
            as.integer(info_dict[["NUM_MARK"]])
        } else {
            0
        }

        # Extract CNLR_MEDIAN (Seg.CN)
        cnlr_median <- if ("CNLR_MEDIAN" %in% names(info_dict)) {
            as.numeric(info_dict[["CNLR_MEDIAN"]])
        } else {
            0
        }

        # Create record with the required fields in the specified order
        record <- list(
            "Sample.ID" = sample_name,
            "Chromosome" = chrom,
            "Start_Position" = pos,
            "End_Position" = end,
            "Num_Markers" = num_mark,
            "Seg.CN" = cnlr_median,
            "SVTYPE" = svtype
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
        message("No valid records found in ", vcf_file)
        return(NULL)
    }

    # Get all keys from the records
    all_keys <- unique(unlist(lapply(records, names)))

    # Define fixed columns - these will always be present in the output in this exact order
    fixed_cols <- c(
        # Required SEG format columns
        "Chromosome", "Start_Position", "End_Position", "Num_Markers", "Seg.CN", "SVTYPE", "Sample.ID", 

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

    df$Start_Position <- as.integer(df$Start_Position - 1)  # Convert to 0-based for bed file
    df$End_Position <- as.integer(df$End_Position)

    # Clean it to bed format
    # df$Chromosome <- gsub("chr", "", df$Chromosome)
    # df$Chromosome <- ifelse(
    #     grepl("^[0-9]+$", df$Chromosome),
    #     as.numeric(df$Chromosome),
    #     df$Chromosome
    # )

    # Return the processed data frame instead of writing to file
    return(df[df$SVTYPE != "NEUTR", ])
}

# Process input (file or directory) and combine results
process_files <- function(input_vcf, output_file) {
    # Check if input is a valid file
    if (!file.exists(input_vcf)) {
        stop("Input file does not exist: ", input_vcf)
    }
    
    # Process single VCF file
    df <- suppressWarnings(process_vcf_to_seg(input_vcf))
    
    # Write to output file if data exists
    if (!is.null(df) && nrow(df) > 0) {
        write.table(
            df,
            file = output_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
        )
        
        message("Saving segment file: ", output_file)
    } else {
        stop("No valid data found to process in ", input_vcf)
    }
}

# If script is being run directly (not sourced)
if (!interactive()) {
    input_vcf <- opt$input
    output_file <- opt$output
    process_files(input_vcf, output_file)
}

# test
# input_dir <- "/home/zhonggr/projects/250224_sarcoma_multiomics/data/wes/variant_calling/cnv/facets"
# input_dir <- "/home/zhonggr/projects/250224_sarcoma_multiomics/data/wes/variant_calling/cnv/facets_test_normal"
# process_files(input_dir)

## Delete all .txt files in the directory
# library(fs)
# dir_ls(input_dir, glob = "*.txt", recurse = TRUE) |> 
# file_delete()