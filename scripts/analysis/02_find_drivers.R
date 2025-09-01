#!/usr/bin/env Rscript
##############################################################################
## Description: Find genomic drivers for FST and Metastasis
## Author:  Zhong Guorui
## Created: 2025-06-15
##############################################################################
## Load required functions and configurations
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Find diff and shared events  -----
## "==========================================================================="
stat_res_list <- LoadData(
    dir = "data/processed",
    filename = "wes_snv_indels_cnv_stat_res"
)

## No significant genes found if use the p_adj_thres
p_val_thres <- 0.05
p_adj_thres <- 0.1
min_altered_samples <- 4

stat_res_sig <- list()
for (type in names(stat_res_list$diff)) {

    stat_res_sig[[type]] <- map(
        stat_res_list$diff[[type]],
        ~ .x |>
            dplyr::filter(fisher_p_val < p_val_thres) |>
            # filter(fisher_p_adj < p_adj_thres) |>
            dplyr::filter(
                group1_altered >= min_altered_samples |
                    group2_altered >= min_altered_samples
            )
    )

    # output_dir <- "outputs/group_comparison"
    output_dir <- "outputs/group_comparison"
    dir_create(here(output_dir))

    if (type == "cnv_gene_bedtools") {

        for (gc in names(stat_res_sig[[type]])) {
            
            data <- stat_res_sig[[type]][[gc]]
            stat_res_sig[[type]][[gc]] <- data |> 
                left_join(
                    cnv_facet_gene_tbl |> 
                        dplyr::select(chr, start, end, symbol, cn_status) |> 
                        # unite("symbol", c("symbol", "cn_status"), sep = "_") |> 
                        distinct(),
                    by = c("symbol", "alteration_type" = "cn_status")
                ) |> 
                relocate(
                    chr, start, end, .after = symbol
                )
        }
        
        for (gc in names(stat_res_list$diff[[type]])) {

            data <- stat_res_list$diff[[type]][[gc]]
            stat_res_list$diff[[type]][[gc]] <- data |> 
                left_join(
                    cnv_facet_gene_tbl |> 
                        dplyr::select(chr, start, end, symbol, cn_status) |> 
                        # unite("symbol", c("symbol", "cn_status"), sep = "_") |> 
                        distinct(),
                    by = c("symbol", "alteration_type" = "cn_status")
                ) |> 
                dplyr::relocate(
                    chr, start, end, .after = symbol
                )
        }

    }

    ## Save signficant results
    write_xlsx(
        stat_res_sig[[type]],
        path = here(
            output_dir,
            paste0("wes_group_comparison_stat_res_sig_", type, ".xlsx")
        )
    )

    ## Save all results
    write_xlsx(
        stat_res_list$diff[[type]],
        path = here(
            output_dir,
            paste0("wes_group_comparison_stat_res_all_", type, ".xlsx")
        )
    )
}

## "========================================================================="
## Figure1: shared events across groups ----
## "========================================================================="
### Load data ----
stat_res_list <- LoadData(
    dir = "data/processed",
    filename = "wes_snv_indels_cnv_stat_res"
)

pcgr_snv_indel_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_snv_indels_tbl"
)

pcgr_cnv_tbl <- LoadData(
    dir = "data/processed",
    filename = "wes_pcgr_cnv_tbl"
)

cnv_facet_gene_tbl <- LoadData(
    dir = "data/processed",
    filename = "cnv_facet_gene_tbl"
)

snv_indel_freq_thres <- 0.05
cnv_gene_freq_thres <- 0.3
cnv_cytoband_freq_thres <- 0.2

## Load cBioportal sarcoma data and DFSP studies
sarcoma_data <- LoadCBioPortalSarcomaData(
    mutate_gene_freq_thres = snv_indel_freq_thres,
    cna_gene_freq_thres = 0.05,
    structural_variants_freq_thres = 0.01
)

sarcoma_cnv_gene <- sarcoma_data$cnv_gene_data |> 
    filter(n_altered >= 3 & freq > cnv_gene_freq_thres) |> 
    pull(symbol)

## Other public DFSP studies
DFSP_study_data <- LoadDFSPStudyData()

public_cnv_gene_DFSP <- DFSP_study_data$smith_2025

## "========================================================================="
### SNV/Indel genes -----
## "========================================================================="
pcgr_snv_indel_mat <- GetVariantMat(
    data = pcgr_snv_indel_tbl,
    variant_type = "symbol", 
    variant_class = "variant_class", 
    is_somatic_matched = TRUE
)

snv_indel_freq <- GetVariantFreq(
    mat = pcgr_snv_indel_mat,
    variant_type = "symbol"
)

snv_indel_gene_cohort <- stat_res_list$share$mutated_gene |>
    filter(n_altered >= 3 & freq > snv_indel_freq_thres) |>
    pull(symbol)

snv_indel_gene_sarcoma <- sarcoma_data$mutated_gene_data |> 
    pull(symbol)

snv_indel_gene_DFSP <- intersect(
    DFSP_study_data$smith_2025$mutated_gene,
    DFSP_study_data$peng_2022$mutated_gene
)

snv_indel_gene_list <- c(
    snv_indel_gene_cohort,
    snv_indel_gene_sarcoma
    # snv_indel_gene_DFSP
) |> unique()

snv_indel_gene_list <- snv_indel_gene_list[
    snv_indel_gene_list %in% rownames(pcgr_snv_indel_mat)
]

snv_indel_gene_list <- snv_indel_freq |> 
    filter(symbol %in% snv_indel_gene_list) |> 
    filter(freq > 0.02) |> 
    pull(symbol)

## Prepare the matrix
pcgr_snv_indel_mat <- pcgr_snv_indel_mat[snv_indel_gene_list, ]

## "========================================================================="
### CNV cytobands -----
## "========================================================================="
cnv_cytoband_mat <- GetVariantTypeMat(
    data = pcgr_cnv_tbl,
    variant_type = "cytoband",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cnv_cytoband_freq <- stat_res_list$share$cnv_cytoband |> 
    unite("cytoband", c(cytoband, alteration_type), sep = "_")

top_cytobands <- cnv_cytoband_freq |> 
    filter(n_altered >=3 & freq > cnv_cytoband_freq_thres)|>
    pull(cytoband) |> 
    unique()

cnv_cytoband_mat <- cnv_cytoband_mat[top_cytobands, ]

## "========================================================================="
### CNV genes -----
## "========================================================================="
## For CNV genes, we fouced on deep deletion and high-level amplification
alteration_types <- c("HOMDEL", "AMP")
cnv_facet_gene_mat <- GetVariantMat(
    data = cnv_facet_gene_tbl |> 
        filter(cn_status %in% alteration_types),
    variant_type = "symbol",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cnv_facet_gene_freq <- stat_res_list$share$cnv_gene_bedtools |> 
    filter(alteration_type %in% alteration_types) |> 
    filter(n_altered >= 3 & freq > cnv_gene_freq_thres)
    # unite("symbol", c(symbol, alteration_type), sep = "_")

top_cnv_genes <- cnv_facet_gene_freq |> 
    slice_max(freq, n = 10, with_ties = TRUE) |> 
    pull(symbol)

## top cytoband genes in the oncokb list
cytobands_gene_data <- pcgr_cnv_tbl |> 
    mutate(cytoband = paste0(cytoband, "_", cn_status)) |>
    filter(cytoband %in% top_cytobands) |> 
    dplyr::select(cytoband, symbol) |> 
    distinct()

top_cnv_genes <- cytobands_gene_data |> 
    filter(symbol %in% top_cnv_genes) |> 
    pull(symbol) |> 
    unique()

cnv_gene_list <- cytobands_gene_data |> 
    filter(symbol %in% cnv_facet_gene_freq$symbol) |> 
    pull(symbol)

## Known cancer driver genes
oncokb_genes <- LoadOncoKBGeneList()
cancer_drivers <- intersect(cnv_gene_list, oncokb_genes)

## cell cycle and DNA repair
cell_cycle_dna_repair <- intersect(cnv_gene_list, c(
    "CDK3", "BIRC5", "TK1", "RAD51C", "BRIP1", "PPM1D", "RECQL5",
    "EME1", "FAAP100", "TBX2", "TBX4", "PRKAR1A", "SKA2", "SIRT7"
))

## Tier 3: Transcription factors and chromatin regulators
transcription_chromatin = intersect(cnv_gene_list, c(
    "SOX9", "FOXK2", "FOXJ1", "CBX2", "CBX4", "CBX8", "HLF", 
    "TBX2", "TBX4", "SMARCD2", "BPTF", "JMJD6", "MED13"
))

## Tier 4: Signal transduction and kinases
signaling_kinases = intersect(cnv_gene_list, c(
    "MAP2K6", "MAP3K3", "GRB2", "RAC3", "RPTOR", "RPS6KB1", 
    "CSNK1D", "SPHK1", "TLK2", "AATK"
))

cnv_genes <- c(
    cell_cycle_dna_repair,
    transcription_chromatin,
    signaling_kinases,
    top_cnv_genes
) |> unique()
    
cnv_gene_mat <- cnv_facet_gene_mat[cnv_genes, ]

## "========================================================================="
### Plot the oncoplot ---------
## "========================================================================="
VerticalCombinedOncoplot(
    pcgr_snv_indel_mat = pcgr_snv_indel_mat,
    cnv_gene_mat = cnv_gene_mat,
    cnv_cytoband_mat = cnv_cytoband_mat,
    is_somatic_matched = TRUE,
    main_type = "cytoband",
    sample_annotation = c("FST.Group", "Metastasis", "Specimen.Nature"),
    column_title = NULL,
    text_size = 6,
    title_size = 7,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "snv_indel_cnv_cytoband_gene_combined_freq_oncoplot"
)

## "========================================================================="
## Figure2:  Pre-FST vs Post-FST ---- 
## "========================================================================="
## No significant mutated gene
## CNV arm: 5q, 20p enriched in Post-FST
## Total 37 cases
## Total samples = 80, Pre-FST=41, Post-FST=39
### CNV cytoband matrix --------------
clinical_info <- LoadClinicalInfo() |> 
    dplyr::filter(FST.Group %in% c("Pre-FST", "Post-FST")) |>
    dplyr::filter(Number.of.samples.per.patient >= 2)

type <- "cnv_gene_bedtools"

n_case <- length(unique(clinical_info$Case.ID))
clinical_info |> 
    group_by(FST.Group) |> 
    summarise(n_cases = n())

cnv_cytoband_mat <- GetVariantTypeMat(
    data = pcgr_cnv_tbl,
    variant_type = "cytoband",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

share_cytobands <- stat_res_list$share$cnv_cytoband[["Pre-FST_vs_Post-FST"]] |> 
    arrange(desc(group2_freq)) |> 
    filter(group1_freq > 0.3 & group2_freq > 0.3) |> 
    mutate(cytoband = paste0(cytoband, "_", alteration_type)) |>
    pull(cytoband)

diff_cytobands <- stat_res_sig$cnv_cytoband[["Pre-FST_vs_Post-FST"]] |> 
    arrange(desc(group2_freq)) |> 
    mutate(cytoband = paste0(cytoband, "_", alteration_type)) |>
    pull(cytoband)

case_event_diff <- cnv_cytoband_mat[diff_cytobands, ] |> 
    as_tibble(rownames = "cytoband") |> 
    pivot_longer(
        -cytoband,
        names_to = "Sample.ID",
        values_to = "cnv_status"
    ) |> 
    left_join(clinical_info, by = "Sample.ID") |> 
    # filter(cytoband %in% diff_cytobands) |> 
    filter(FST.Group  %in% c("Pre-FST", "Post-FST")) |> 
    dplyr::filter(Number.of.samples.per.patient >= 2) |> 
    group_by(Case.ID, cytoband) |> 
    summarise(
        pre_fst_status = cnv_status[FST.Group == "Pre-FST"][1],
        post_fst_status = cnv_status[FST.Group == "Post-FST"][1],
        .groups = "drop"
    )  |> 
    filter(!is.na(pre_fst_status) & !is.na(post_fst_status)) |>
    mutate(
        event_type = case_when(
            pre_fst_status == "" & post_fst_status != "" ~ "Acquired",
            pre_fst_status != "" & post_fst_status == "" ~ "Lost", 
            pre_fst_status != "" & post_fst_status != "" ~ "Maintained",
            pre_fst_status == "" & post_fst_status == "" ~ "Absent"
        )
    ) |> 
    group_by(cytoband, event_type) |>
    summarise(
        n_cases = n(),
        .groups = "drop"
    )  |> 
    pivot_wider(
        names_from = event_type, 
        values_from = n_cases, 
        values_fill = 0
    ) |>
    arrange(desc(Acquired))

case_event_share <- cnv_cytoband_mat[share_cytobands, ] |> 
    as_tibble(rownames = "cytoband") |> 
    pivot_longer(
        -cytoband,
        names_to = "Sample.ID",
        values_to = "cnv_status"
    ) |> 
    left_join(clinical_info, by = "Sample.ID") |> 
    # filter(cytoband %in% diff_cytobands) |> 
    filter(FST.Group  %in% c("Pre-FST", "Post-FST")) |> 
    dplyr::filter(Number.of.samples.per.patient >= 2) |> 
    group_by(Case.ID, cytoband) |> 
    summarise(
        pre_fst_status = cnv_status[FST.Group == "Pre-FST"][1],
        post_fst_status = cnv_status[FST.Group == "Post-FST"][1],
        .groups = "drop"
    )  |> 
    filter(!is.na(pre_fst_status) & !is.na(post_fst_status)) |>
    mutate(
        event_type = case_when(
            pre_fst_status == "" & post_fst_status != "" ~ "Acquired",
            pre_fst_status != "" & post_fst_status == "" ~ "Lost", 
            pre_fst_status != "" & post_fst_status != "" ~ "Maintained",
            pre_fst_status == "" & post_fst_status == "" ~ "Absent"
        )
    ) |> 
    group_by(cytoband, event_type) |>
    summarise(
        n_cases = n(),
        .groups = "drop"
    )  |> 
    pivot_wider(
        names_from = event_type, 
        values_from = n_cases, 
        values_fill = 0
    ) |>
    arrange(desc(Maintained))

cnv_cytobands  <- unique(
    c(
        share_cytobands, 
        diff_cytobands[!diff_cytobands %in% c("chr17:p11.2_DEL", "chr17:q24.3-q25.1_GAIN")]
    )
)

cnv_cytoband_mat <- cnv_cytoband_mat[cnv_cytobands, ]

### CNV gene matrix -------
chroms <- sapply(
    cnv_cytobands,
    function(x) strsplit(x, ":")[[1]][1]
)

cnv_gene_share <- stat_res_list[["share"]][[type]][["Pre-FST_vs_Post-FST"]] |> 
    filter(group1_freq > 0.3 & group2_freq > 0.3) |>
    arrange(desc(group2_freq))

cnv_gene_diff <- stat_res_sig[[type]][["Pre-FST_vs_Post-FST"]] |> 
    # filter(grepl("chr6", chr)) |>
    filter(fisher_p_val < 0.01) |> 
    mutate(symbol = paste0(symbol, "_", alteration_type)) |>
    filter(chr  %in% chroms) |>
    group_by(chr) |> 
    slice_head(prop = 0.5) |> 
    pull(symbol)

cnv_gene_mat <- GetVariantTypeMat(
    data = cnv_facet_gene_tbl,
    variant_type = "symbol",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cnv_gene_mat <- cnv_gene_mat[cnv_gene_diff, ]

## Plot the oncoplot
PairedPatientOncoplot(
    cnv_gene_mat = cnv_gene_mat,
    cnv_cytoband_mat = cnv_cytoband_mat,
    is_somatic_matched = TRUE,
    paired_group = c("Pre-FST", "Post-FST"),
    column_sort = c("Case.ID", "FST.Group"),
    sample_annotation = c("FST.Group", "Metastasis", "Specimen.Nature"),
    column_title = NULL,
    title_size = 7,
    text_size = 6,
    width = 14,
    height = 6,
    dir = "figures/oncoplot",
    filename = "paired_samples_Pre-FST_vs_Post-FST_oncoplot"
)

## "========================================================================="
## Figure3: signficant events across groups ---- 
## "========================================================================="
## All mutated genes tends to not signficant if use padj < 0.05
cnv_cytoband_mat <- GetVariantTypeMat(
    data = pcgr_cnv_tbl,
    variant_type = "cytoband",
    variant_class = "cn_status",
    is_somatic_matched = TRUE
)

cnv_cytobands <- list_rbind(stat_res_sig$cnv_cytoband) |> 
    dplyr::select(cytoband, alteration_type) |>
    unite("cytoband", c(cytoband, alteration_type), sep = "_") |>
    distinct() |> 
    pull(cytoband)

cnv_cytoband_mat <- cnv_cytoband_mat[cnv_cytobands, ]

sort_group_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")

for (i in 2:length(sort_group_level)) {

    ComplexGroupOncoplot(
        mat = cnv_cytoband_mat,
        is_somatic_matched = TRUE,
        # split_by_group = "FST.Type",
        split_by_group = "FST.Group",
        sort_group_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"),
        # sort_group_level = c("Classic", "FST"),
        main_group = sort_group_level[[i]],
        sample_annotation = c("FST.Group", "Metastasis", "Specimen.Nature"),
        column_title = NULL,
        width = 18,
        height = 12,
        dir = "figures/oncoplot",
        filename = paste0("freq_sort_by_", sort_group_level[i])
    )
}

# gene_entrez <- bitr(
#     geneID = unique(top_cytobands_genes$symbol),
#     fromType = "SYMBOL",
#     toType = "ENTREZID",
#     OrgDb = "org.Hs.eg.db"
# ) |> 
#     filter(!is.na(ENTREZID))

# go_bp <- enrichGO(
#     gene = gene_entrez$ENTREZID,
#     OrgDb = org.Hs.eg.db,
#     ont = "BP",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.1,
#     readable = TRUE
# )
# head(go_bp@result, 10)




