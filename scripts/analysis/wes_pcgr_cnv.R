#!/usr/bin/env Rscript
##############################################################################
## Description: This script processes PCGR data and performs cohort analysis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required libraries and functions
suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(circlize)
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(GenomicRanges)
        library(ggsci)
        library(ComplexHeatmap)
        library(reshape2)
    })
)
source(here::here("scripts/lib/study_lib.R"))
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Preprocess PCGR CNV data --------------
## "==========================================================================="


## "========================================================================="
## Different and shared cytobands---- 
## "========================================================================="
pcgr_cnv_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_cnv_tbl"
)

## Analyze the somatic matched samples
clinical_info <- LoadClinicalInfo()
cohort_size <- clinical_info |> nrow()

pcgr_cytoband_tbl <- pcgr_cnv_tbl |> 
    group_by(sample_id, cytoband) |> 
    summarise(
        variant_class = paste(unique(cn_status), collapse = ","),
        .groups = "drop"
    ) |> 
    pivot_wider(
        names_from = sample_id,
        values_from = variant_class,
        values_fill = ""
    )

pcgr_cytoband_stat_res <- list()
for (gc in names(group_comparisons)) {

    group1 <- group_comparisons[[gc]][["group1"]]
    group2 <- group_comparisons[[gc]][["group2"]]

    if (gc == "Pre-FST_vs_Post-FST") {

        ## Paired  samples
        group1_samples <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(FST.Group %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(FST.Group %in% group2) |> 
            pull(Sample.ID)

    } else if (gc == "Primary_vs_Metastasis") {
        
        ## Get the case that have metastatsis
        case_ids <- clinical_info |> 
            filter(Number.of.samples.per.patient >= 2) |> 
            filter(grepl("Metastasis", Specimen.Nature)) |> 
            pull(Case.ID)

        group1_samples <- clinical_info |> 
            filter(Case.ID %in% case_ids) |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(Specimen.Nature %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(Case.ID %in% case_ids) |> 
            filter(Number.of.samples.per.patient >= 2) |>
            filter(Specimen.Nature %in% group2) |> 
            pull(Sample.ID)

    } else {

        group1_samples <- clinical_info |> 
            filter(FST.Group %in% group1) |> 
            pull(Sample.ID)
        
        group2_samples <- clinical_info |> 
            filter(FST.Group %in% group2) |> 
            pull(Sample.ID)
    }

    pcgr_cytoband_stat_res[[gc]] <- GetGroupStatDiff(
        data = pcgr_cytoband_tbl,
        group1_samples = group1_samples,
        group2_samples = group2_samples,
        group_comparison = gc,
        alternative = "two.sided",
        p_adjust_method = "BH"
    )
}

data <- cytoband_tbl
group1 <- "U-DFSP"
group2 <- "FS-DFSP"

group1_samples <- clinical_info |> 
    filter(FST.Group %in% group1) |> 
    pull(Sample.ID)
        
group2_samples <- clinical_info |> 
    filter(FST.Group %in% group2) |> 
    pull(Sample.ID)

group_comparison <- "U-DFSP_vs_FS-DFSP"

    filter(sample_id %in% clinical_info$Sample.ID) |> 
    left_join(clinical_info, by = c("sample_id" = "Sample.ID")) |> 
    relocate(FST.Group, .after = sample_id)

## Get the significantly different cytobands across groups
all_diff_cytobands <- map_dfr(
    FST_group_comparisons,
    function(comp) {
        message(paste0(" - Processing ", comp$name))

        GetPCGRCNVGroupDiffCytoband(
            data = cna_data,
            var = "FST.Group",
            cn_column = "cn_status",
            event_type = NULL,
            somatic_matched = TRUE,
            group1 = comp$group1,
            group2 = comp$group2,
            p_adjust_method = "fdr",
            min_altered_samples = 3,
            group_comparison = comp$name
        )
    }
) |>
    mutate(
        cytoband_alteration = paste0(cytoband, "_", cn_status),
        .after = "cn_status"
    )

all_diff_cytobands <- all_diff_cytobands |> 
    mutate(
        chrom = sapply(str_split(cytoband, ":"), function(x) x[1]),
        .before = "cytoband"
    ) |> 
    mutate(chrom = str_replace(chrom, "chr", "")) |> 
    filter(!(chrom %in% c("X", "Y"))) |>
    mutate(chrom = as.numeric(chrom))

## Fisher exact test p value threshold
p_val_thres <- 0.05

## At least 15% samples are mutated
freq_thres <- 0.15

sig_diff_cytobands <- all_diff_cytobands |> 
    filter(fisher_p < p_val_thres) |>
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |>
    group_by(cytoband_alteration, cytoband, cn_status) |> 
    summarise(
        n_comparisons_present = n(),
        comparisons_present = paste(group_comparison, collapse = ", "),
        .groups = "drop"
    ) |> 
    arrange(desc(n_comparisons_present)) |> 
        mutate(
        chrom = sapply(str_split(cytoband, ":"), function(x) x[1]),
        .before = "cytoband"
    ) |> 
    mutate(chrom = str_replace(chrom, "chr", "")) |> 
    filter(!(chrom %in% c("X", "Y"))) |>
    mutate(chrom = as.numeric(chrom))

## Test some peaks
selected_peaks <- c(
    ## U-DFSP -> Pre-FST
    "chr7:q36.1_DEL",
    "chr17:q21.33-q22_AMP",
    "chr22:q11.21-q11.23_GAIN",

    ## Pre-FST -> Post-FST
    "chr17:q23.1-q22_GAIN",
    "chr19:q13.41-q13.33_DEL",

    ## U-DFSP -> Post-FST
    "chr1:p36.33_DEL",
    "chr14:q32.33_DEL",
    "chr17:q21.2_DEL",
    "chr20:q13.33_DEL",

    ## Pre-FST -> FS-DFSP
    "chr16:p13.3_GAIN",
    "chr17:p12-p13.1_GAIN",
    "chr2:q37.3_DEL",
    "chr6:p22.1_DEL",

    ## Post-FST -> FS-DFSP
    "chr17:q25.3_GAIN",
    "chr19:p13.11_DEL",
    "chr19:q13.12-q13.11_DEL",
    "chr7:p22.3_GAIN",

    ## U-DFSP -> FS-DFSP
    "chr14:q11.2-q12_DEL",
    "chr22:q11.21-q11.23_AMP",
    "chr2:p25.1_HOMDEL",
    "chr5:q31.3_GAIN",
    "chr9:p21.3_HOMDEL",

    ## U-DFSP -> Pre-FST
    "chr1:p36.32_DEL",
    "chr1:p36.33_DEL",
    "chr2:p25.1_HOMDEL",
    "chr2:q37.3_DEL",
    "chr3:q27.1_DEL",
    "chr5:p15.33_GAIN",
    "chr5:q31.3_AMP",
    "chr5:q31.3_GAIN",
    "chr6:p21.33-p21.32_HOMDEL",
    "chr6:p22.1-p21.33_DEL",
    "chr6:q21_DEL",
    "chr6:q22.2-q23.2_DEL",
    "chr6:p22.1_DEL",
    "chr7:q22.1_AMP",
    "chr7:p22.3_GAIN",
    "chr7:q36.1_DEL"
)

all_diff_cytobands |> 
    filter(fisher_p < p_val_thres) |>
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |> 
    filter(
        cytoband_alteration %in% c(
            "chr7:q36.1_DEL",
            "chr17:q21.33-q22_AMP",
            "chr22:q11.21-q11.23_GAIN"
        )
    )

## Get the shared cytobands across groups
shared_cytobands <- GetPCGRCNVGroupShareCytoband(
    data = cna_data,
    var = "FST.Group",
    cn_column = "cn_status",
    somatic_matched = TRUE,
    min_freq_per_group = 0.3
)

shared_cytobands$shared_cytobands <- shared_cytobands$shared_cytobands |> 
    mutate(
        chrom = sapply(
            str_split(cytoband, ":"), 
            function(x) x[1]
        ),
        .after = "cytoband"
    ) |> 
    mutate(chrom = str_replace(chrom, "chr", "")) |> 
    filter(!(chrom %in% c("X", "Y"))) |>
    mutate(chrom = as.numeric(chrom))

## Save the cytoband stat data
SaveData(
    list(
        all_diff_cytobands = all_diff_cytobands,
        sig_diff_cytobands = sig_diff_cytobands,
        shared_cytobands = shared_cytobands
    ),
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_somatic_matched_cytoband_stats"
)

## Save the different and shared regions across groups
write_xlsx(
    list(
        all_diff_cytobands = all_diff_cytobands,
        sig_diff_cytobands = sig_diff_cytobands,
        shared_cytobands = shared_cytobands$shared_cytobands
    ),
    here(
        "outputs", "wes_pcgr_DFSP_cohort_somatic_matched_diff_shared_cytoband.xlsx"
    )
)

## "========================================================================"
## Cytobands oncoplot ----
## "========================================================================"
## Cytoband freq data
cytoband_freq <- cna_data |> 
    select(sample_id, cytoband) |>
    distinct() |>
    group_by(cytoband) |> 
    summarise(
        n_samples_altered = n_distinct(sample_id),
        .groups = "drop"
    ) |> 
    mutate(
        cohort_freq = n_samples_altered / cohort_total_samples
    ) |> 
    arrange(desc(cohort_freq))

all_diff_cytobands |> 
    filter(fisher_p < p_val_thres) |>
    filter(
        cytoband  %in% c(
            # "chr5:p15.33"
            "chr9:p21.3"
            # "chr22:q13.1-q12.3"
        )
    ) |> 
    filter(group_comparison == "U-DFSP_vs_FS-DFSP")

## Prepare the oncoplot data
cytoband_oncoplot_data <- cna_data  |> 
    select(sample_id, cytoband, cn_status) |>
    distinct() |>
    group_by(sample_id, cytoband) |>
    summarise(
        alterations = paste(unique(cn_status), collapse = ","),
        .groups = "drop"
    ) |> 
    arrange(desc(alterations)) |> 
    mutate(
        final_alteration = case_when(
            ## Count number of different alteration types (split by comma)
            str_count(alterations, ",") >= 1 ~ "MULTI",  # Multiple different alterations
            ## Single alterations
            alterations == "HOMDEL" ~ "HOMDEL",
            alterations == "DEL" ~ "DEL", 
            alterations == "GAIN" ~ "GAIN",
            alterations == "AMP" ~ "AMP",
            TRUE ~ alterations
        )
    ) |> 
    distinct()

## Select the cytobands to plot
selected_cytobands <- all_diff_cytobands |> 
    filter(fisher_p < p_val_thres) |> 
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |> 
    mutate(cytoband_alteration = paste0(cytoband, "_", cn_status)) |> 
    pull(cytoband) |> 
    unique()

## Prepare the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    cytoband = selected_cytobands
) |> 
    left_join(
        cytoband_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "cytoband")
    )

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, cytoband, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("cytoband") |>
    as.matrix()

dim(onco_mat)

## Plot the oncoplot
plot_para <- list(
    entire_cohort  = list(
        sample_sorted_by = "FST.Group",
        sample_sorted_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"),
        sample_annotation = c(
            "FST.Group", "Specimen.Nature", "Metastasis",
            "Meth.Subtype.Main.2Class", "HRD.cat", "MSI.cat"
        )
    ),
    paired_samples = list(
        sample_sorted_by = NULL,
        sample_sorted_level = NULL,
        sample_annotation = c(
            "FST.Group", "Specimen.Nature", "Metastasis", 
            "Meth.Subtype.Main.2Class", "HRD.cat", "MSI.cat"
        )
    ),
    EPIC_meth_group = list(
        sample_sorted_by = "Meth.Subtype.Main.2Class",
        sample_sorted_level = c("Meth1", "Meth2"),
        sample_annotation = c(
            "Meth.Subtype.Main.2Class", 
            "FST.Group", "Specimen.Nature", "Metastasis", 
            "HRD.cat", "MSI.cat"
        )
    ),
    HRD_group = list(
        sample_sorted_by = "HRD.cat",
        sample_sorted_level = c("Low", "High"),
        sample_annotation = c(
            "HRD.cat", "FST.Group", "Specimen.Nature", "Metastasis", 
            "Meth.Subtype.Main.2Class", "MSI.cat"
        )
    ),
    MSI_group = list(
        sample_sorted_by = "MSI.cat",
        sample_sorted_level = c("Low", "High"),
        sample_annotation = c(
            "MSI.cat", "FST.Group", "Specimen.Nature", "Metastasis",
            "Meth.Subtype.Main.2Class", "HRD.cat"
        )
    )
)

sample_info <- LoadClinicalInfo() |> 
    filter(Somatic.Status == "Matched")

for (i in names(plot_para)) {

    ## all samples
    all_sample_ids <- colnames(onco_mat)

    if (i == "paired_samples") {

        sample_ids <- sample_info |> 
            filter(grepl("Paired", Histology.Nature)) |> 
            pull(Sample.ID)
        
        samples <- all_sample_ids[all_sample_ids %in% sample_ids]

        n_samples <- length(samples)

        column_title <- paste0(
                "Somatic matched samples, N = ", n_samples, ", Paired samples only",
                "\n", "(Sorted by Case.ID)"
        )

        filename <- "wes_cnvfacets_pcgr_oncoplot_somatic_matched_cytobands_paired_only"

    } else if (i == "EPIC_meth_group") {

        sample_ids <- sample_info |> 
            filter(!is.na(Meth.Subtype.Main.2Class)) |> 
            pull(Sample.ID)
        
        samples <- all_sample_ids[all_sample_ids %in% sample_ids]

        n_samples <- length(samples)

        column_title <- paste0(
                "Somatic matched samples, N = ", n_samples,
                "\n", "(Sorted by ", plot_para[[i]]$sample_sorted_by, ")"
        )

        filename <- paste0(
            "wes_cnvfacets_pcgr_oncoplot_somatic_matched_cytobands_", 
            plot_para[[i]]$sample_sorted_by
        )

    } else {
    
        samples <- all_sample_ids
        
        n_samples <- length(samples)

        column_title <- paste0(
                "Somatic matched samples, N = ", n_samples,
                "\n", "(Sorted by ", plot_para[[i]]$sample_sorted_by, ")"
        )

        filename <- paste0(
            "wes_cnvfacets_pcgr_oncoplot_somatic_matched_cytobands_", 
            plot_para[[i]]$sample_sorted_by
        )
    }
    
    ## Subset the matrix according to groups
    mat <- onco_mat[, samples, drop = FALSE]
    
    ## Plot the oncoplot
    PCGRComplexOncoplot(
        mat = mat,
        sample_annotation = plot_para[[i]]$sample_annotation,
        sample_sorted_by = plot_para[[i]]$sample_sorted_by,
        sample_sorted_level = plot_para[[i]]$sample_sorted_level,
        column_title = column_title,
        width = 20,
        height = 10,
        dir = "figures/oncoplot",
        filename = filename
    )
}

## "==========================================================================="
## Cytoband chromosomal plot --------------
## "==========================================================================="
# peak_label_list <- list(
    
#     "U-DFSP" = all_diff_cytobands |> 
#         filter(
#             group_comparison %in% c("U-DFSP_vs_Pre-FST")
#         ) |> 
#         pull(cytoband_alteration) |> 
#         unique(),
        
#     "Pre-FST" = all_diff_cytobands |> 
#         filter(
#             group_comparison %in% c("U-DFSP_vs_Pre-FST")
#         ) |> 
#         pull(cytoband_alteration) |> 
#         unique(),

#     "Post-FST" = all_diff_cytobands |> 
#         filter(
#             group_comparison %in% c(
#                 "U-DFSP_vs_Post-FST",
#                 "U-DFSP_vs_FS-DFSP",
#                 "Pre-FST_vs_Post-FST",
#                 "Pre-FST_vs_FS-DFSP",
#                 "Post-FST_vs_FS-DFSP",
#                 "U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP"
#             )
#         ) |> 
#         pull(cytoband_alteration) |> 
#         unique(),
    
#     "FS-DFSP" = all_diff_cytobands |> 
#         filter(
#             group_comparison %in% c(
#                 "U-DFSP_vs_Post-FST",
#                 "U-DFSP_vs_FS-DFSP",
#                 "Pre-FST_vs_Post-FST",
#                 "Pre-FST_vs_FS-DFSP",
#                 "Post-FST_vs_FS-DFSP",
#                 "U-DFSP+Pre-FST_vs_Post-FST+FS-DFSP"
#             )
#         ) |> 
#         pull(cytoband_alteration) |> 
#         unique()
# )

gistic_dir <- "data/wes/GISTIC2/somatic_matched"

plots <- list()
for (group in FST.Group) {
    
    ## Load the gistic score
    gistic_score <- LoadGisticScore(
        gistic_dir = gistic_dir,
        group = group
    )

    ## the sample info in entire cohort
    n_samples <- LoadClinicalInfo() |> 
        filter(Somatic.Status == "Matched") |> 
        filter(FST.Group == group) |> 
        nrow()

    # peak_label <- peak_label_list[[group]]
    
    ## Get the peak coordinates
    peaks_data <- GetPCGRGisticOverlap2(
        peak_label = selected_peaks,
        pcgr_cna_data = cna_data, 
        gistic_score_data = gistic_score,
        min_overlap = 0.1
    )

    peaks_data <- AddPeaksPCGRCountData(
        peaks_data,
        pcgr_cna_data = cna_data,
        var = "FST.Group",
        group = group,
        is_somatic_matched = TRUE
    )

    ## Plot the enriched regions
    plots[[group]] <- GisticChromPlot(
        data = gistic_score,
        peaks_data = peaks_data,
        y_lim = c(-2, 3),
        x_expand = c(0.02, 0),
        title = paste0(group, " (", n_samples, ")"),
        show_horizontal_grid = TRUE,
        show_vertical_grid = TRUE,
        grid_line_type = "dotted",
        grid_line_color = "grey80",
        grid_line_alpha = 0.5,
        grid_line_width = 0.2,
        show_chromosome_ideogram = TRUE,
        ideogram_label_size = 4,
        rel_heights = c(4, 0.3)
    )
}

SavePlot(
    plot = wrap_plots(plots, ncol = 2) +
        plot_layout(guides = "collect"),
    width = 12,
    height = 5,
    dir = "figures/chromplot",
    filename = "wes_facet_pcgr_DFSP_cohort_gistic_plots_FST.Group"
)

## "========================================================================="
## Different and shared CNV genes -------
## "========================================================================="
pcgr_cnv_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_cnv_filtered"
)

## Analyze the somatic matched samples
clinical_info <- LoadClinicalInfo() |>
    filter(Somatic.Status %in% c("Matched"))

cohort_total_samples <- clinical_info |> nrow()

cna_data <- pcgr_cnv_tbl |> 
    filter(sample_id %in% clinical_info$Sample.ID) |> 
    left_join(clinical_info, by = c("sample_id" = "Sample.ID")) |> 
    relocate(FST.Group, .after = sample_id) |> 
    mutate(
        gene_alteration = paste0(symbol, "_", cn_status),
        .after = symbol
    )

## Get the significantly different genes across groups
all_diff_genes <- map_dfr(
    FST_group_comparisons,
    function(comp) {
        GetPCGRCNVGroupDiffGene(
            data = cna_data,
            var = "FST.Group",
            cn_column = "cn_status",
            event_type = NULL,
            somatic_matched = TRUE,
            group1 = comp$group1,
            group2 = comp$group2,
            p_adjust_method = "fdr",
            min_altered_samples = 3,
            group_comparison = comp$name
        )
    }
)

all_diff_genes |> 
    filter(symbol %in% c("RB1", "TP53")) |> 
    filter(fisher_p < p_val_thres)

sig_diff_genes <- all_diff_genes |> 
    filter(fisher_p < p_val_thres) |>
    filter(group1_freq > freq_thres | group2_freq > freq_thres) |>
    mutate(gene_alteration = paste0(symbol, "_", cn_status)) |>
    group_by(gene_alteration, symbol, cn_status) |>
    summarise(
        n_comparisons_present = n(),
        comparisons_present = paste(group_comparison, collapse = ", "),
        .groups = "drop"
    ) |> 
    arrange(desc(n_comparisons_present))

## Get the shared genes across groups
shared_genes <- GetPCGRCNVGroupShareGene(
    data = cna_data,
    var = "FST.Group",
    cn_column = "cn_status",
    somatic_matched = TRUE,
    min_freq_per_group = 0.3
)

## Save the gene stat data
SaveData(
    list(
        all_diff_genes = all_diff_genes,
        sig_diff_genes = sig_diff_genes,
        shared_genes = shared_genes
    ),
    dir = "data/processed",
    filename = "wes_pcgr_DFSP_cohort_somatic_matched_gene_stats"
)

## Save the different and shared regions across groups
write_xlsx(
    list(
        all_diff_genes = all_diff_genes,
        sig_diff_genes = sig_diff_genes,
        shared_genes = shared_genes$shared_genes
    ),
    here(
        "outputs", "wes_pcgr_DFSP_cohort_somatic_matched_diff_shared_gene.xlsx"
    )
)

## "========================================================================"
## CNV genes oncoplot public cohort ----
## "========================================================================"
## Calculate the frequency for all CNV genes
gene_freq <- cna_data |>
    select(sample_id, symbol) |>
    distinct() |> 
    group_by(symbol) |> 
    summarise(
        n_samples_altered = n_distinct(sample_id),
        .groups = "drop"
    ) |> 
    mutate(cohort_freq = n_samples_altered / cohort_total_samples) |>
    arrange(desc(cohort_freq))

## Prepare the oncoplot data
gene_oncoplot_data <- cna_data  |> 
    select(sample_id, symbol, cn_status) |>
    distinct() |>
    group_by(sample_id, symbol) |>
    summarise(
        alterations = paste(unique(cn_status), collapse = ","),
        .groups = "drop"
    ) |> 
    arrange(desc(alterations)) |> 
    mutate(
        final_alteration = case_when(
            ## Count number of different alteration types (split by comma)
            str_count(alterations, ",") >= 1 ~ "MULTI",  # Multiple different alterations
            ## Single alterations
            alterations == "HOMDEL" ~ "HOMDEL",
            alterations == "DEL" ~ "DEL", 
            alterations == "GAIN" ~ "GAIN",
            alterations == "AMP" ~ "AMP",
            TRUE ~ alterations
        )
    ) |> 
    distinct()

## cBioportal sarcoma CNV genes
sarcoma_cna_file <- here("data/public/cBioPortal_sarcoma_CNA_Genes.txt")
col_names <- c(
    "symbol", 
    "gistic_q_value",
    "cytoband",
    "cn_status",
    "profiled_samples",
    "#",
    "freq",
    "is_oncoKB"
)
sarcoma_cna_data <- read_tsv(sarcoma_cna_file, show_col_types = FALSE)
colnames(sarcoma_cna_data) <- col_names

sarcoma_cna_data <- sarcoma_cna_data |> 
    mutate(freq = parse_number(freq)) |> 
    arrange(desc(freq)) |>
    filter(freq >= 10) |>
    distinct(symbol, .keep_all = TRUE)

## Define public cohorts
dfsp_public_cohorts <- list(
    cbioportal_sarcoma = list(
        title = "Top mutated genes (cBioPortal Sarcoma, Freq > 1%)",
        cnv = sarcoma_cna_data$symbol
    ),
    ## Modern pathology
    smith_2025 = list(
        title = "Top mutated genes (Smith et al. 2025 Reported, n=53)",
        cnv = c(
            "PRKAR1A", "H3F3B", "MSI2", 
            "SRSF2", "BRIP1", "DDX5", "CLTC", 
            "GNA13", "CRKL", "MCL1", "MYH9",
            "SEPT9", "SOX10",
            "ATP1A1", "CANT1", "CLTCL1", "DAXX", "ERG", "EZR",
            "FGFR1OP",
            "LIFR",
            "MAF",
            "NOTCH2",
            "PDGFRB",
            "PRRX1",
            "RAC1",
            "RICTOR",
            "RNF43"
        )
    ),
    ## BJD
    peng_2022 = list(
        title = "Top mutated genes (Peng et al. 2022 Reported, n=59)",
        cnv = c(
            "TERT", "CDKN2A", "CDKN2B", 
            "AKT1", "NFKBIA", "BTBD7", "SPHK1", 
            "ITGB4", "COL1A1", "PDGFB", "PDGFD"
        )
    )
)

## Genes to plot
selected_genes <- c(
    dfsp_public_cohorts$smith_2025$cnv,
    dfsp_public_cohorts$peng_2022$cnv,
    ## well known oncodriver gene
    "TP53", "PDGFRB", "AKT"
) |> 
    unique()

## Construct the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    symbol = selected_genes
) |> 
    left_join(
        gene_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "symbol")
    ) |> 
    arrange(final_alteration)

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, symbol, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("symbol") |>
    as.matrix()

## Plot the oncoplot
sample_sorted_by <- "FST.Group"  # Define which variable to sort by
sample_sorted_level <- c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
column_title <- paste0(
    "Somatic matched samples, N = ", cohort_total_samples, 
    "\n", "Public cohort (Peng et al. 2022; Smith et al. 2025)"
)

filename <- paste0(
    "wes_cnvfacets_pcgr_oncoplot_public_DFSP_cohort_CNV_genes_somatic_matched_", 
    sample_sorted_by
)

PCGRComplexOncoplot(
    mat = onco_mat,
    sample_annotation = c(
        "FST.Group", "Metastasis", "Specimen.Nature",
        "Meth.Subtype.Main.2Class", "HRD.cat", "MSI.cat"
    ),
    sample_sorted_by = sample_sorted_by,
    sample_sorted_level = sample_sorted_level,
    column_title = column_title,
    width = 20,
    height = 10,
    dir = "figures/oncoplot",
    filename = filename
)

## "========================================================================"
## CNV HRD genes oncoplot ----
## "========================================================================"
## Genes to plot
selected_genes <- HRD_genes

## Construct the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    symbol = selected_genes
) |> 
    left_join(
        gene_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "symbol")
    ) |> 
    arrange(final_alteration)

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, symbol, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("symbol") |>
    as.matrix()

## Plot the oncoplot
sample_sorted_by <- "FST.Group"
sample_sorted_level <- c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
# sample_sorted_by <- "HRD.cat"
# sample_sorted_level <- c("Low", "High")
sample_annotation = c(
    "FST.Group", "Specimen.Nature", "Metastasis", 
    "Meth.Subtype.Main.2Class", "MSI.cat", "HRD.cat"
)

column_title <- paste0(
    "Somatic matched samples, N = ", cohort_total_samples, 
    "\n", "HRD genes"
)

filename <- paste0(
    "wes_cnvfacets_pcgr_oncoplot_CNV_HRD_genes_somatic_matched_", 
    sample_sorted_by
)

PCGRComplexOncoplot(
    mat = onco_mat,
    sample_annotation = sample_annotation,
    sample_sorted_by = sample_sorted_by,
    sample_sorted_level = sample_sorted_level,
    column_title = column_title,
    width = 20,
    height = 10,
    dir = "figures/oncoplot",
    filename = filename
)

## "========================================================================"
## CNV MMR genes oncoplot ----
## "========================================================================"
## Genes to plot
selected_genes <- MMR_genes

## Construct the oncoplot data
onco_data <- expand_grid(
    sample_id = clinical_info$Sample.ID,
    symbol = selected_genes
) |> 
    left_join(
        gene_oncoplot_data |> select(-alterations), 
        by = c("sample_id", "symbol")
    ) |> 
    arrange(final_alteration)

onco_mat <- onco_data |> 
    mutate(
        final_alteration = if_else(
            is.na(final_alteration), "", final_alteration
        )
    ) |> 
    ## Handle duplicates by keeping distinct combinations
    distinct(sample_id, symbol, .keep_all = TRUE) |>
    pivot_wider(
        names_from = sample_id,
        values_from = final_alteration,
        values_fill = ""
    ) |>
    column_to_rownames("symbol") |>
    as.matrix()

## Plot the oncoplot
sample_sorted_by <- "FST.Group"
sample_sorted_level <- c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
# sample_sorted_by <- "MSI.cat"
# sample_sorted_level <- c("Low", "High")
sample_annotation = c(
    "MSI.cat", "FST.Group", "Specimen.Nature", "Metastasis",
    "Meth.Subtype.Main.2Class", "HRD.cat"
)

column_title <- paste0(
    "Somatic matched samples, N = ", cohort_total_samples, 
    "\n", "MMR genes"
)

filename <- paste0(
    "wes_cnvfacets_pcgr_oncoplot_CNV_MMR_genes_somatic_matched_", 
    sample_sorted_by
)

PCGRComplexOncoplot(
    mat = onco_mat,
    sample_annotation = sample_annotation,
    sample_sorted_by = sample_sorted_by,
    sample_sorted_level = sample_sorted_level,
    column_title = column_title,
    width = 20,
    height = 10,
    dir = "figures/oncoplot",
    filename = filename
)
