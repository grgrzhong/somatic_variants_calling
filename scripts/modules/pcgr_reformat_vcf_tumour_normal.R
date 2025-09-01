#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(vcfR)
  library(argparse)
})

#' Reformat VCF file to add tumor/normal depth and allele frequency information
#' 
#' @param vcf_file Path to input VCF file
#' @param output_file Path to output VCF file
#' @return None (writes output file)
reformat_vcf <- function(vcf_file, output_file) {
  # Read VCF file
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  
  # Extract genotype data
  gt_data <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)
  af_data <- extract.gt(vcf, element = 'AF', as.numeric = TRUE)
  
  # Check if we have at least 2 samples (normal and tumor)
  if (ncol(gt_data) < 2) {
    stop("VCF file must contain at least 2 samples (normal and tumor)")
  }
  
  # Assuming sample order: [0] = normal, [1] = tumor
  normal_dp <- gt_data[, 1]
  tumor_dp <- gt_data[, 2]
  normal_af <- af_data[, 1]
  tumor_af <- af_data[, 2]
  
  # Add INFO fields to header
  info_header <- vcf@meta[grepl("^##INFO", vcf@meta)]
  new_info <- c(
    '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor sample depth">',
    '##INFO=<ID=NDP,Number=1,Type=Integer,Description="Normal sample depth">',
    '##INFO=<ID=TAF,Number=1,Type=Float,Description="Tumor sample AF">',
    '##INFO=<ID=NAF,Number=1,Type=Float,Description="Normal sample AF">'
  )
  
  # Insert new INFO fields after existing INFO fields
  info_end <- max(which(grepl("^##INFO", vcf@meta)))
  vcf@meta <- c(vcf@meta[1:info_end], new_info, vcf@meta[(info_end + 1):length(vcf@meta)])
  
  # Add new INFO fields to the INFO column
  for (i in 1:nrow(vcf@fix)) {
    info_string <- vcf@fix[i, "INFO"]
    
    # Handle missing values
    tdp_val <- if (is.na(tumor_dp[i])) "." else as.character(tumor_dp[i])
    ndp_val <- if (is.na(normal_dp[i])) "." else as.character(normal_dp[i])
    taf_val <- if (is.na(tumor_af[i])) "." else sprintf("%.6f", tumor_af[i])
    naf_val <- if (is.na(normal_af[i])) "." else sprintf("%.6f", normal_af[i])
    
    # Append new INFO fields
    new_info_fields <- paste0(";TDP=", tdp_val, ";NDP=", ndp_val, ";TAF=", taf_val, ";NAF=", naf_val)
    vcf@fix[i, "INFO"] <- paste0(info_string, new_info_fields)
  }
  
  # Write modified VCF
  write.vcf(vcf, file = output_file)
  cat("Successfully reformatted VCF file:", vcf_file, "-> output:", output_file, "\n")
}

# Create argument parser
parser <- ArgumentParser(description = 'Process input and output VCF files.')

# Add arguments
parser$add_argument('-I', '--input', type = 'character', required = TRUE,
                    help = 'Path to the input VCF file')
parser$add_argument('-O', '--output', type = 'character', required = TRUE,
                    help = 'Path to the output VCF file')

# Parse command-line arguments
args <- parser$parse_args()

# Access the values of the arguments
input_vcf <- args$input
output_vcf <- args$output

# Call the main function
tryCatch({
  reformat_vcf(input_vcf, output_vcf)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  quit(status = 1)
})

