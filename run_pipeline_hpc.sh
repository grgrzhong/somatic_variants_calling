#!/bin/bash

sbatch ./somatic_calling_pipeline.sh \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data/Raw \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data \
    2 \
    /lustre1/g/path_my/Reference

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/MFS/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/MFS \
#     16 \
#     /lustre1/g/path_my/Reference