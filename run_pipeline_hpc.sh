#!/bin/bash

sbatch ./somatic_calling_pipeline.sh \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/GLI1/Raw \
    /lustre1/g/path_my/pipeline/somatic_variants_calling/data/GLI1 \
    6 \
    /lustre1/g/path_my/Reference \
    /lustre1/g/path_my/Software

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/BLMS/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/BLMS \
#     2 \
#     /lustre1/g/path_my/Reference \
#     /lustre1/g/path_my/Software

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/AFH/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/AFH \
#     10 \
#     /lustre1/g/path_my/Reference \
#     /lustre1/g/path_my/Software

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/SARC/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/SARC \
#     20 \
#     /lustre1/g/path_my/Reference \
#     /lustre1/g/path_my/Software

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/IM/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/IM \
#     20 \
#     /lustre1/g/path_my/Reference \
#     /lustre1/g/path_my/Software

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/DFSP \
#     8 \
#     /lustre1/g/path_my/Reference

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/MFS/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/MFS \
#     16 \
#     /lustre1/g/path_my/Reference

# sbatch ./somatic_calling_pipeline.sh \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data/Raw \
#     /lustre1/g/path_my/pipeline/somatic_variants_calling/data/test_data \
#     2 \
#     /lustre1/g/path_my/Reference \
#     /lustre1/g/path_my/Software