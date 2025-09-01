#!/usr/bin/env python3

import argparse
from pysam import VariantFile


def reformat_vcf(vcf_file, out):
    with VariantFile(vcf_file) as fr:
        header = fr.header
        header.info.add('TDP', number=1, type='Integer', description='Tumor sample depth')
        header.info.add('NDP', number=1, type='Integer', description='Normal sample depth')
        header.info.add('TAF', number=1, type='Float', description='Tumor sample AF')
        header.info.add('NAF', number=1, type='Float', description='Normal sample AF')

        with VariantFile(out, 'w', header=header) as fw:
            for record in fr:
                record.info['TDP'] = record.samples[1]['DP']
                record.info['NDP'] = record.samples[0]['DP']
                record.info['TAF'] = record.samples[1]['AF'][0]
                record.info['NAF'] = record.samples[0]['AF'][0]
                fw.write(record)
                
                
# Create an argument parser object
parser = argparse.ArgumentParser(description='Process input and output files.')

# Add arguments to the parser
parser.add_argument('-I', '--input', type=str, required=True, help='Path to the input file')
parser.add_argument('-O', '--output', type=str, required=True, help='Path to the output file')

# Parse the command-line arguments
args = parser.parse_args()


# Access the values of the arguments
input_vcf = args.input
output_vcf = args.output

reformat_vcf(input_vcf, output_vcf)

