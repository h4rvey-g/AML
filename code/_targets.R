library(targets)
library(crew)
tar_source(c(
    "code/R/101.load_data.R",
    "code/R/102.cluster_annotate.R",
    "code/R/103.DEG.R",
    "code/R/104.fine_analysis.R",
    "code/R/105.myeloid.R",
    "code/R/106.Tcell.R"
))
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "rliger", "Seurat", "SeuratExtend", "tidyseurat", "scCustomize", "patchwork", "tidyplots",
        "liana", "SingleCellExperiment", "scDblFinder", "foreach", "doParallel", "sccomp"
    ),
    controller = crew_controller_local(workers = 20, seconds_timeout = 6000),
    format = "qs",
    storage = "worker", retrieval = "worker",
    seed = 42
)
tar_config_set(
    script = "code/_targets.R",
    store = "data/_targets"
)
options(max.print = 12, spe = "human")
list(
    tar_target(sc_raw, load_sc()),
    tar_target(sc_raw_filtered, filter_low_quality_cells(sc_raw)),
    tar_target(sc_int_preview, process_sc_preview(sc_raw_filtered)),
    tar_target(sc_doublet, find_doublet(sc_int_preview)), # 添加新的 target
    tar_target(sc_int, process_sc(sc_raw_filtered, sc_doublet)),
    tar_target(sc_cluster, cluster_data(sc_int)),
    tar_target(cell_type_path, "results/102.cluster_annotate/cell_type.tsv", format = "file"),
    tar_target(sc_annotate, annotate_data(sc_cluster, cell_type_path)),
    tar_target(sc_sub_stero, sub_cluster_steroidogenic(sc_annotate)),
    tar_target(sc_final, final_annotation(sc_annotate, sc_sub_stero)),
    tar_target(anno_DEG_path, find_DEGs_pseudo_bulk(sc_final), format = "file"),
    # tar_target(anndata, save_annotate(sc_final)),
    tar_target(cell_types, unique(sc_final$cell_type_dtl)),
    tar_target(DEG_path, DEG_annotation(sc_final, cell_types), format = "file", pattern = map(cell_types)),
    tar_target(GSEA_path, run_GSEA(sc_final), format = "file"),
    tar_target(composition_test, test_distribution(sc_final)),
    tar_target(distribution_path, plot_distribution(sc_final, composition_test), format = "file"),
    tar_target(parents, {
        options(max.print = 12, spe = "human")
        GeneSetAnalysisGO() %>% names()
    }),
    tar_target(
        sc_filtered,
        sc_final %>% filter(cell_type == "Plasticity")
    ),
    tar_target(GSEA_Plasticity_path, run_GSEA_Plasticity(sc_filtered, parents),
        format = "file",
        pattern = map(parents)
    ),
    tar_target(GSEA_Plasticity_hallmark_path, compare_cell_types_hallmark50(sc_filtered),
        format = "file"
    ),
    tar_target(calculate_DEG_mCAF_vs_Others_path, calculate_DEG_mCAF_vs_Others(sc_filtered),
        format = "file"
    ),
    tar_target(cell_communication_path, cell_communication(sc_final)),
    tar_target(sc_mye_clust, sub_cluster_myeloid(sc_final)),
    tar_target(sc_mye, myeloid_annotate(sc_mye_clust)),
    tar_target(mye_distribution_path, plot_myeloid_distribution(sc_mye), format = "file"),
    tar_target(mye_GSEA_path, run_myeloid_GSEA(sc_mye), format = "file"),
    tar_target(mye_DEG_path, myeloid_DEG_tumor_vs_normal(sc_mye), format = "file"),
    tar_target(myeloid_plot, plot_myeloid(sc_mye)),
    tar_target(sc_tcell_clust, sub_cluster_tcell(sc_final)),
    tar_target(sc_tcell_clust_impro, Tcell_annotation_analysis(sc_tcell_clust)),
    tar_target(sc_tcell, Tcell_annotate(sc_tcell_clust_impro)),
    tar_target(t_composition_test, test_tcell_composition(sc_tcell)),
    tar_target(tcell_distribution_path, plot_tcell_distribution(sc_tcell, t_composition_test), format = "file"),
    tar_target(tcell_GSEA_path, run_tcell_GSEA(sc_tcell), format = "file"),
    tar_target(tcell_DEG_path, tcell_DEG_tumor_vs_normal(sc_tcell), format = "file")
)
