import sys
import pysam

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]

reader = pysam.VariantFile(input_vcf)
reader.header.formats.add('GT', '1', 'String', 'Genotype')

with pysam.VariantFile(output_vcf, 'w', header=reader.header) as writer:
    for record in reader:
        for sample in record.samples:
            if 'GT' not in record.samples[sample].keys():
                gt = './.'
                if 'DP' in record.samples[sample].keys() and record.samples[sample]['DP'] > 0:
                    gt = '0/0'
                record.samples[sample]['GT'] = gt
        writer.write(record)

reader.close()
