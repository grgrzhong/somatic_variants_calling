import pysam
import argparse

#input_bam = "/media/Data/WES/BAM/GLI-001-N/GLI-001-N.bam"  # Change this to your aligned BAM file
#output_bam = "/media/Data/WES/BAM/GLI-001-N/GLI-001-N_umi.bam"  # Output file with UMI tag


# Create an argument parser object
parser = argparse.ArgumentParser(description='Process input and output files.')

# Add arguments to the parser
parser.add_argument('-I', '--input', type=str, required=True, help='Path to the input file')
parser.add_argument('-O', '--output', type=str, required=True, help='Path to the output file')

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the arguments
input_bam = args.input
output_bam = args.output

with pysam.AlignmentFile(input_bam, "rb") as infile, pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
    for read in infile:
        read_name = read.query_name
        umi = read_name.split(":")[-1]  # Extract UMI from read name
        umi = umi.replace("_", "-") # Replace the UMI delimiter
        read.query_name = ":".join(read.query_name.split(":")[:-1]) # Remove UMI from read names
        read.set_tag("RX", umi)  # Set UMI as a custom tag "RX"
        outfile.write(read)

