
## Load required libraries and functions
library(maftools)
source(here::here("scripts/lib/study_lib.R"))

## Merge all the annovar annoated variants -------------
data_dir <- "data/processed"
figure_dir <- "figures/wes"

## "=========================================================================="
## Collect the annovar annotated variants
## "=========================================================================="
annovar_dir <- "data/wes/Mutect2"
annovar_tbl <- CollectAnnovarData(dir = annovar_dir)
filename <- "wes_annovar_DFSP_cohort_merged_tbl"
SaveData(annovar_tbl, dir = data_dir, filename = filename)

## Filter variants
annovar_tbl <- LoadData(dir = data_dir, filename = filename)

filtered_variants <- FilterAnnovarData(data = annovar_tbl)

## Find the enriched variants in FST ----
maf_obj <- read.maf(maf = filtered_variants)

## Sample groups
sample_groups <- list(
    `U-DFSP` = c("Classic", "Myxoid", "Pigmented"),
    `Pre_FST` = c(
        "Pretransformed classic",
        "Pretransformed myxoid",
        "Paired classic",
        "Paired myxoid"
    ),
    `Post_FST` = c(
        "Posttransformed FST",
        "Paired FST",
        "Paired Pleomorphic"
    ),
    `FS_DFSP` = c("Unpaired FST")
)

sample_info <- LoadSampleInfo() |> 
    filter(Specimen.Class == "Tumour") |> 
    select(
        Sample.ID, Diagnosis, Specimen.Class, Specimen.Nature, Histology.Nature,
        Somatic.Status,purity,ploidy
    ) |> 
    mutate(
        sample_group = case_when(
            Histology.Nature %in% sample_groups$`U-DFSP` ~ "U-DFSP",
            Histology.Nature %in% sample_groups$`Pre_FST` ~ "Pre-FST",
            Histology.Nature %in% sample_groups$`Post_FST` ~ "Post-FST",
            Histology.Nature %in% sample_groups$`FS_DFSP` ~ "FS-DFSP",
            TRUE ~ "Other"
        )
    ) |> 
    rename(Tumor_Sample_Barcode = Sample.ID)

## Two samples were not appeared in the maf data: "DFSP-139-T", "DFSP-294-T-M1"
# setdiff(
#     sort(sample_info |> pull(Tumor_Sample_Barcode)),
#     sort(maf_obj@clinical.data$Tumor_Sample_Barcode)
# )

## Add the sample group to the maf data
maf_obj@clinical.data <- maf_obj@clinical.data |> 
    select(Tumor_Sample_Barcode) |>
    left_join(sample_info, by = "Tumor_Sample_Barcode")

annotation_colors <- list(
    sample_group = c(
        "U-DFSP"    = "#3498db", # Blue - for untransformed DFSP
        "Pre-FST"   = "#2ecc71", # Green - for pre-transformation samples
        "Post-FST"  = "#e74c3c", # Red - for post-transformation samples
        "FS-DFSP"   = "#9b59b6"  # Purple - for unpaired FST
    ),
    Specimen.Nature = c(
        "Primary"       = "#ff9800", # Orange
        "Recurrence"    = "#009688", # Teal
        "Metastasis"    = "#795548", # Brown
        "Residual"      = "#607d8b"  # Blue-gray
    )
)

n_samples <- nrow(getSampleSummary(maf_obj))
top_n_genes <- 30
clinical_features <- c("sample_group", "Specimen.Nature")

oncoplot(
    maf = maf_obj,
    clinicalFeatures = clinical_features,
    annotationColor = annotation_colors,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)