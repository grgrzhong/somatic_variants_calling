#!/usr/bin/env Rscript

suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(ComplexHeatmap)
        library(circlize)
    })
)

source(here::here("scripts/lib/study_lib.R"))
source(here::here("conf/study_conf.R"))

## "==========================================================================="
## Oncoplot genomic driver event for FST and metastasis ----
## "==========================================================================="
cytoband_stat_res <- LoadData(
    cytoband_stat_res,
    dir = "data/processed",
    filename = "wes_cnvfacets_gistic2_group_comparisons_stat_res_all"
)

## At least 15% samples are mutated
p_val_thres <- 0.05
p_adj_thres <- 0.1
min_altered_samples <- 3

# Combine all group_comparisons into one dataframe
sig_stat_res <- map(
    cytoband_stat_res,
    function(gc) {
        gc |>
            filter(fisher_p_value < p_val_thres) |> 
            # filter(fisher_p_adj < p_adj_thres) |>
            filter(
                group1_altered >= min_altered_samples |
                group2_altered >= min_altered_samples
            )
    }
)

write_xlsx(
    sig_stat_res,
    path = here(
        "outputs",
        "wes_cnvfacets_gistic2_group_comparisons_stat_res_sig.xlsx"
    )
)

## Run the fisher's exact test for each FST group comparison
# stat_res <- list()

# for (gc in names(FST_group_comparisons)) {

#     group1 <- FST_group_comparisons[[gc]][["group1"]]
#     group2 <- FST_group_comparisons[[gc]][["group2"]]

#     message(sprintf(" - Processing %s", gc))

#     stat_res[[gc]] <- GetGisticCytobandGroupStatRes(
#         gistic_dir = gistic_dir,
#         var = "FST.Group",
#         group1 = group1,
#         group2 = group2,
#         group_comparison = gc,
#         is_somatic_matched = TRUE
#     )
# }

# ## At least 15% samples are mutated
# freq_thres <- 0.15
# p_val_thres <- 0.05
# p_adj_thres <- 0.1
# min_altered_samples <- 3

# # Combine all group_comparisons into one dataframe
# sig_cytoband <- map_dfr(
#     names(stat_res),
#     function(gc) {
#         stat_res[[gc]] |>
#             filter(fisher_p_value < p_val_thres) |> 
#             filter(fisher_p_adj < p_adj_thres) |>
#             filter(
#                 group1_altered >= min_altered_samples |
#                 group2_altered >= min_altered_samples
#             )
#     }
# )

# sig_cytoband_summary <- sig_cytoband |> 
#     group_by(Unique_Name, Cytoband, Variant_Classification) |> 
#     summarise(
#         n_comparisons_present = n(),
#         comparisons_present = paste(
#             group_comparison, collapse = ","
#         ),
#         .groups = "drop"
#     ) |> 
#     arrange(desc(n_comparisons_present)) |> 
#         mutate(
#         Chromosome = sapply(str_split(Cytoband, "q|p"), function(x) x[1]),
#         .before = "Cytoband"
#     ) |> 
#     filter(!(Chromosome %in% c("X", "Y"))) |>
#     mutate(Chromosome = as.numeric(Chromosome)) |> 
#     arrange(Chromosome)

# cytoband_stat_res <- list(
#     sig_cytoband = sig_cytoband,
#     sig_cytoband_summary = sig_cytoband_summary
# )

# ## Save the results
# SaveData(
#     cytoband_stat_res,
#     dir = "data/processed",
#     filename = "wes_cnvfacets_gistic2_group_stat_res_list"
# )

# write_xlsx(
#     cytoband_stat_res,
#     path = here(
#         "outputs",
#         "wes_cnvfacets_gistic2_cytoband_group_stat_res_somatic_matched.xlsx"
#     )
# )

sig_cytoband_data <- list_rbind(sig_stat_res)
sig_cytobands <- unique(sig_cytoband_data$Unique_Name)

## Get gistic copy number data
gistic_data <- GetGisticCNData(
    gistic_dir = gistic_dir,
    group = "all_tumors",
    gistic_qval_thres = 0.5
)

cn_mat <- gistic_data$cn_data |> 
    filter(Unique_Name %in% sig_cytobands) |> 
    column_to_rownames("Unique_Name") |> 
    as.matrix()

## Plot the oncoplot
GisticComplexGroupOncoplot(
    mat = cn_mat,
    is_somatic_matched = TRUE,
    split_by_group = "FST.Group",
    sort_group_level = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP"),
    main_group = "FS-DFSP",
    sample_annotation = c("FST.Group", "Metastasis", "Specimen.Nature"),
    column_title = NULL,
    width = 18,
    height = 12,
    dir = "figures/oncoplot",
    filename = "wes_cnvfacets_oncoplot_entire_cohort_cytoband"
)

## "==========================================================================="
## Gistic oncoplot all tumors ----
## "==========================================================================="
gistic_dirs <- list(
    somatic_matched = here("data/wes/GISTIC2/somatic_matched"),
    somatic_unmatched = here("data/wes/GISTIC2/somatic_unmatched")
)

gistic_dir <- here("data/wes/GISTIC2/somatic_matched")

gistic_obj <- readGistic(
    gisticAllLesionsFile = here(gistic_dir, "all_tumors", "all_lesions.conf_99.txt"),
    gisticAmpGenesFile = here(gistic_dir, "all_tumors", "amp_genes.conf_99.txt"),
    gisticDelGenesFile = here(gistic_dir, "all_tumors", "del_genes.conf_99.txt"),
    gisticScoresFile = here(gistic_dir, "all_tumors", "scores.gistic"),
    verbose = FALSE
)

sig_cytobands <- sig_cytoband_data |> 
    pull(cytoband_alteration) |> 
    unique()

all_tumors_cytobands_data <- gistic_obj@data |> 
    as_tibble() |> 
    mutate(
        cytoband = sapply(
            Cytoband,
            function(x) str_split(x, ":")[[1]][2]
        )
    ) |> 
    mutate(
        cytoband_status = paste(Variant_Classification, cytoband, sep = "_")
    ) |> 
    filter(cytoband %in% all_cytobands$Cytoband) |> 
    pull(Cytoband) |> 
    unique() |> 
    sort()

table(sig_cytobands %in% all_tumors_cytobands_data$cytoband_status)

gistic_obj@gis.scores
gistic_obj@numericMatrix

clinical_data <- LoadClinicalInfo() |> 
    mutate(
        Tumor_Sample_Barcodes = Sample.ID
    ) |> 
    select(
        Tumor_Sample_Barcodes, 
        Specimen.Nature, 
        FST.Group, 
        Metastasis
    ) |> 
    as.data.frame()

gisticOncoPlot(
    gistic = gistic_obj,
    # bands = show_cytobands,
    clinicalData = clinical_data,
    clinicalFeatures = "FST.Group",
    showTumorSampleBarcodes = TRUE,
    sortByAnnotation = TRUE,
    # colors = c(
    #     Amp = "#D95F02",
    #     Del = "#1B9E77"
    # ),
    SampleNamefontSize = 0.6,
    fontSize = 0.8,
    legendFontSize = 0.7,
    annotationFontSize = 1.2,
    borderCol = "white",
    bgCol = "#CCCCCC"
)

## "==========================================================================="
## Gistic plot group comparsions ---------
## "==========================================================================="
plot_para <- list(
    plot_dir = "figures/gistic2",
    gistic_data = list(
        epic_cnv = list(
            gistic_dir = here("data/epic/GISTIC2"),
            plot_filename = "epic_cnv_gistic2"
        ),
        somatic_matched_samples = list(
            gistic_dir = here("data/wes/GISTIC2/somatic_matched"),
            plot_filename = "wes_cnvfacets_gistic2_somatic_matched_samples"
        ),
        all_samples = list(
            gistic_dir = here("data/wes/GISTIC2/somatic_unmatched"),
            plot_filename = "wes_cnvfacets_gistic2_all_samples"   
        )
    ),
    group_list = list(
        Metastatsis = c("Non-Meta", "Meta"),
        FST.subtype = c("U-DFSP", "Pre-FST", "Post-FST", "FS-DFSP")
    )
)

for (i in names(plot_para$gistic_data)) {

    if (i == "somatic_matched") {

        is_somatic_matched <- TRUE

    } else {
        
        is_somatic_matched <- FALSE
    }

    for (group in names(plot_para$group_list)) {
        
        message(
            paste0("\nPlotting GISTIC2 results for: ", group)
        )

        gistic_dir <- plot_para$gistic_data[[i]]$gistic_dir

        groups <- plot_para$group_list[[group]]

        filename <- paste0(
            plot_para$gistic_data[[i]]$plot_filename, "_", group
        )

        plot_dir <- plot_para$plot_dir

        PlotGisticGroupComparsion(
            gistic_dir = gistic_dir,
            groups = groups,
            is_somatic_matched = is_somatic_matched,
            fdrCutOff = 0.25,
            markBands = show_cytobands,
            color = c("#c82118", "#2c5496"),
            ref.build = "hg38",
            cytobandOffset = 0.05,
            txtSize = 0.8,
            cytobandTxtSize = 0.7,
            y_lims = NULL,
            width = 8,
            height = 8,
            plot_dir = plot_dir,
            filename = filename
        )
    }
}

## "==========================================================================="
## Explore the GISTIC2 results ----
## "==========================================================================="
clinical_info <- LoadClinicalInfo()

## Output directories
gistic_dir <- "data/wes/GISTIC2/somatic_matched"
# group_name <- "all_tumors"

sample_groups <- LoadSampleGroupInfo()
# groups <- dir_ls(gistic_dir) |> path_file()
groups <- names(sample_groups)[1:11]

for (group in groups) {

    ## Read the GISTIC2 data
    gistic_obj <- readGistic(
        gisticAllLesionsFile = here(
            gistic_dir, group, "all_lesions.conf_99.txt"
        ),
        gisticAmpGenesFile = here(gistic_dir, group, "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, group, "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, group, "scores.gistic"),
        verbose = TRUE
    )

    n_sample <- length(sample_groups[[group]])

    ## Plot parameters
    width <- 8
    height <- 4
    out_dir <- "figures/wes/gistic2/cnv_facets"

    dir_create(out_dir)

    ## Plot the chromosomal plot
    for (img in c("png", "pdf")) {
        file <- here(
            out_dir,
            paste0("wes_cnvfacets_gistic2_", group, "_chromplot.", img)
        )

        if (img == "png") {
            CairoPNG(
                file = file, width = width, height = height, res = 300,
                fonts = "Arial", units = "in"
            )
        } else {
            CairoPDF(
                file = file, width = width, height = height, fonts = "Arial"
            )
        }

        gisticChromPlot(
            gistic = gistic_obj,
            fdrCutOff = 0.25,
            txtSize = 0.6,
            cytobandTxtSize = 0.5,
            color = c("#c82118", "#2c5496"),
            markBands = "all",
            ref.build = "hg38",
            y_lims = c(-2, 5)
        )
        title(
            main = paste0(group, " (", n_sample, ")"),
            family = "Arial"
        )
        
        message(paste0("Saving plot: ", file))
        dev.off()
    }
}
    # ## Bubble plot
    # for (img in c("png", "pdf")) {
    #     file <- here(
    #         out_dir,
    #         paste0("wes_cnvfacets_gistic2_", group, "_bubbleplot.", img)
    #     )

    #     if (img == "png") {
    #         CairoPNG(
    #             file = file, width = width, height = height, res = 300,
    #             fonts = "Arial", units = "in"
    #         )
    #     } else {
    #         CairoPDF(
    #             file = file, width = width, height = height, fonts = "Arial"
    #         )
    #     }

    #     gisticBubblePlot(
    #         gistic = gistic_obj,
    #         color = c("#D95F02", "#1B9E77"),
    #         markBands = NULL,
    #         log_y = TRUE,
    #         fdrCutOff = 0.25,
    #         txtSize = 0.6
    #     )

    #     message(paste0("Saving plot: ", file))
    #     dev.off()
    # }

    # ## Oncoplot
    # for (img in c("png", "pdf")) {
    #     file <- here(
    #         out_dir,
    #         paste0("wes_cnvfacets_gistic2_", group, "_oncoplot.", img)
    #     )

    #     if (img == "png") {
    #         CairoPNG(
    #             file = file, width = width, height = height, res = 300,
    #             fonts = "Arial", units = "in"
    #         )
    #     } else {
    #         CairoPDF(
    #             file = file, width = width, height = height, fonts = "Arial"
    #         )
    #     }

    #     gisticOncoPlot(
    #         gistic = gistic_obj,
    #         sortByAnnotation = FALSE,
    #         top = 10,
    #         gene_mar = 10,
    #         barcode_mar = 10,
    #         sepwd_genes = 0.5,
    #         bandsToIgnore = NULL,
    #         removeNonAltered = TRUE,
    #         colors = c(
    #             Amp = "#D95F02",
    #             Del = "#1B9E77"
    #         ),
    #         SampleNamefontSize = 0.6,
    #         fontSize = 0.8,
    #         legendFontSize = 0.7,
    #         annotationFontSize = 1.2,
    #         borderCol = "white",
    #         bgCol = "#CCCCCC"
    #     )

    #     message(paste0("Saving plot: ", file))
    #     dev.off()
    # }


## "==========================================================================="
## Non-Meta vs Meta chromosomal plot ---------
## "==========================================================================="
## Non-Meta vs Meta - Combined Plot with Aligned Axes
for (img in c("png", "pdf")) {
    file <- here(
        out_dir,
        paste0("wes_cnvfacets_gistic2_NonMeta_vs_Meta_combined.", img)
    )

    if (img == "png") {
        CairoPNG(
            file = file, width = 8, height = 8, res = 300,
            fonts = "Arial", units = "in"
        )
    } else {
        CairoPDF(
            file = file, width = 8, height = 8, fonts = "Arial"
        )
    }

    # Set up 2 rows, 1 column layout
    par(mfrow = c(2, 1), mar = c(3, 4, 3, 2))

    # Read GISTIC objects for both groups
    gistic_obj_nonmeta <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "Non-Meta", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "Non-Meta", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "Non-Meta", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "Non-Meta", "scores.gistic"),
        verbose = FALSE
    )
    
    gistic_obj_meta <- readGistic(
        gisticAllLesionsFile = here(gistic_dir, "Meta", "all_lesions.conf_99.txt"),
        gisticAmpGenesFile = here(gistic_dir, "Meta", "amp_genes.conf_99.txt"),
        gisticDelGenesFile = here(gistic_dir, "Meta", "del_genes.conf_99.txt"),
        gisticScoresFile = here(gistic_dir, "Meta", "scores.gistic"),
        verbose = FALSE
    )

    n_sample_nonmeta <- length(sample_groups[["Non-Meta"]])

    n_sample_meta <- length(sample_groups[["Meta"]])

    # Plot 1: Non-Meta (top panel, no x-axis labels)
    gisticChromPlot(
        gistic = gistic_obj_nonmeta,
        fdrCutOff = 0.25,
        txtSize = 0.6,
        cytobandTxtSize = 0.5,
        color = c("#c82118", "#2c5496"),
        ref.build = "hg38",
        y_lims = c(-2, 5)  # Same y-limits for both plots
    )
    title(main = paste0("Non-Meta (n=", n_sample_nonmeta, ")"), family = "Arial", cex.main = 1.2)

    # Plot 2: Meta (bottom panel, with x-axis labels)
    par(mar = c(5, 4, 1, 2))  # More bottom margin for x-axis labels
    gisticChromPlot(
        gistic = gistic_obj_meta,
        fdrCutOff = 0.25,
        txtSize = 0.6,
        cytobandTxtSize = 0.5,
        color = c("#c82118", "#2c5496"),
        ref.build = "hg38",
        y_lims = c(-2, 5)  # Same y-limits for both plots
    )
    title(
        main = paste0("Meta (n=", n_sample_meta, ")"), 
        family = "Arial", cex.main = 1.2
    )

    # Reset par settings
    par(mfrow = c(1, 1))
    
    message(paste0("Saving combined plot: ", file))
    dev.off()
}


# PlotGistic2Chrom <- 

# scores <- "/mnt/f/projects/sarcoma_multiomics/data/wes/GISTIC2/cnv_facets/FS-DFSP/scores.gistic"

# library(BSgenome.Hsapiens.UCSC.hg38)
# chrom_info <- tibble(
#     chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38),
#     chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
# )
# chrom_info$chromNum <- 1:length(chrom_info$chromName)
# chrom_info <- chrom_info[1:22, ]

# ## Calculate cumulative chromosome lengths
# chrom_info$chromlengthCumsum <- cumsum(as.numeric(chrom_info$chromlength))

# ## Calculate chromosome start positions from 0
# chrom_info$chromStartPosFrom0 <- c(0, chrom_info$chromlengthCumsum[-nrow(chrom_info)])

# ## Calculate chromosome middle positions from 0
# tmp_middle <- diff(c(0, chrom_info$chromlengthCumsum)) / 2
# chrom_info$chromMiddlePosFrom0 <- chrom_info$chromStartPosFrom0 + tmp_middle

# scores <- read.table(
#     file = scores, 
#     header = TRUE, 
#     sep = "\t", 
#     stringsAsFactors = FALSE
# )
# ## Add chromosome ID and start/end positions
# chromID <- scores$Chromosome
# scores$StartPos <- scores$Start + chrom_info$chromStartPosFrom0[chromID]
# scores$EndPos <- scores$End + chrom_info$chromStartPosFrom0[chromID]

# range(scores$G.score)
# ## Adjust G-scores for deletions for visualization
# # scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
# scores <- scores |> 
#     mutate(
#         G.score = case_when(
#             Type == "Del" ~ G.score * -1,
#             Type == "Amp" ~ G.score
#         )
#     )

# library(ggplot2)
# library(ggsci)

# ggplot(scores, aes(StartPos, G.score)) +
#     geom_area(aes(group = Type, fill = factor(Type, levels = c("Del", "Amp")))) +
#     scale_fill_lancet(guide = guide_legend(reverse = T)) +
#     geom_vline(data = chrom_info, mapping = aes(xintercept = chromlengthCumsum), linetype = 2) +
#     geom_text(data = chrom_info, aes(x = chromMiddlePosFrom0, y = 0.2, label = chromName)) +
#     scale_x_continuous(expand = c(0, -1000), limits = c(0, 2.9e9)) +
#     ylim(-0.3, 0.3) +
#     theme_minimal()

# chrom_info$ypos <- rep(c(0.2, 0.25), 11)

# ggplot(scores, aes(StartPos, G.score)) +
#     geom_area(aes(group = Type, fill = factor(Type, levels = c("Del", "Amp")))) +
#     scale_fill_lancet(guide = guide_legend(reverse = T), name = "Type") +
#     geom_vline(data = chrom_info, mapping = aes(xintercept = chromlengthCumsum), linetype = 2) +
#     geom_text(data = chrom_info, aes(x = chromMiddlePosFrom0, y = ypos, label = chromName)) +
#     scale_x_continuous(expand = c(0, -1000), limits = c(0, 2.9e9), name = NULL, labels = NULL) +
#     ylim(-0.3, 0.3) +
#     theme_minimal() +
#     theme(
#         legend.position = "top",
#         axis.text.y = element_text(color = "black", size = 14),
#         axis.title.y = element_text(color = "black", size = 16)
#     )
