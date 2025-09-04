#!/bin/bash

sbatch ./somatic_calling_pipeline.sh \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data/Raw \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/results \
    2 \
    2 \
    /lustre1/g/path_my/Reference