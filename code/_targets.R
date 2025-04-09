library(targets)
library(crew)
tar_source(c(
    "code/R/101.load_data.R",
    "code/R/102.cluster_annotate.R",
    "code/R/103.DEG.R",
    "code/R/104.fine_analysis.R",
    "code/R/105.myeloid.R",
    "code/R/106.Tcell.R",
    "code/R/107.paper_plot.R",
    "code/R/108.adipo.R",
    "code/R/102.cluster_annotate copy.R"
))
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "rliger", "Seurat", "reticulate", "SeuratExtend", "tidyseurat", "scCustomize", "patchwork", "tidyplots",
        "liana", "SingleCellExperiment", "scDblFinder", "foreach", "doParallel", "scLink", "sccomp", "ggalign", "escape",
        "RColorBrewer"
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
    tar_target(sc_int_2, integrate_data_2(sc_int)),
    tar_target(sc_cluster_2, cluster_data_2(sc_int_2)),
    tar_target(cell_type_path, "results/102.cluster_annotate/cell_type.tsv", format = "file"),
    tar_target(sc_annotate_2, annotate_data_2(sc_cluster_2)),
    tar_target(sc_primary_2, primary_annotation(sc_annotate_2)),
    tar_target(sc_primary_opt, optimize_annotation(sc_primary_2)),
    tar_target(sc_final_2, final_annotation_2(sc_primary_opt)),
    tar_target(sc_final, sc_final_2),
    tar_target(anno_DEG_path, find_DEGs_pseudo_bulk(sc_final), format = "file"),
    # tar_target(anndata, save_annotate(sc_final)),
    # tar_target(DEG_path, DEG_annotation(sc_final, cell_types), format = "file", pattern = map(cell_types)),
    # tar_target(GSEA_path, run_GSEA(sc_final), format = "file"),
    tar_target(composition_test, test_distribution(sc_final)),
    tar_target(distribution_path, plot_distribution(sc_final, composition_test), format = "file"),
    tar_target(parents, {
        options(max.print = 12, spe = "human")
        GeneSetAnalysisGO() %>% names()
    }),
    # tar_target(
    #     sc_filtered,
    #     sc_final %>% filter(cell_type == "Plasticity")
    # ),
    # tar_target(GSEA_Plasticity_path, run_GSEA_Plasticity(sc_filtered, parents),
    #     format = "file",
    #     pattern = map(parents)
    # ),
    # tar_target(GSEA_Plasticity_hallmark_path, compare_cell_types_hallmark50(sc_filtered),
    #     format = "file"
    # ),
    # tar_target(calculate_DEG_Fib_vs_Others_path, calculate_DEG_mCAF_vs_Others(sc_filtered),
    #     format = "file"
    # ),
    tar_target(cell_communication_path, cell_communication(sc_final)),
    tar_target(sc_mye_clust, sub_cluster_myeloid(sc_final)),
    tar_target(sc_mye_annotation, sub_annotation_myeloid(sc_mye_clust, sc_final)),
    tar_target(sc_mye, myeloid_annotate(sc_mye_annotation)),
    tar_target(mye_cluster, unique(sc_mye$cell_type_dtl)),
    tar_target(mye_DEG_cluster_path, myeloid_DEG_vs_cluster(sc_mye, mye_cluster), format = "file", pattern = map(mye_cluster)),
    tar_target(mye_composition_test, test_myeloid_composition(sc_mye)),
    tar_target(mye_distribution_path, plot_myeloid_distribution(sc_mye, mye_composition_test), format = "file"),
    tar_target(mye_GSEA_path, run_myeloid_GSEA(sc_mye), format = "file"),
    tar_target(mye_DEG_path, myeloid_DEG_tumor_vs_normal(sc_mye, mye_cluster), pattern = map(mye_cluster)),
    tar_target(mye_trajectory_path, run_myeloid_trajectory(sc_mye), format = "file"),
    tar_target(sc_tcell_clust_list, sub_cluster_tcell(sc_final)),
    tar_target(sc_tcell_clust_impro_list, Tcell_annotation_analysis(sc_tcell_clust_list, sc_final)),
    tar_target(sc_tcell, Tcell_annotate(sc_tcell_clust_impro_list, sc_final)),
    tar_target(t_composition_test, test_tcell_composition(sc_tcell)),
    tar_target(tcell_distribution_path, plot_tcell_distribution(sc_tcell, t_composition_test), format = "file"),
    tar_target(tcell_GSEA_path, run_tcell_GSEA(sc_tcell), format = "file"),
    tar_target(t_cell_type, unique(sc_tcell$cell_type_dtl)),
    tar_target(tcell_DEG_path, tcell_DEG_tumor_vs_normal(sc_tcell, t_cell_type), format = "file", pattern = map(t_cell_type)),
    tar_target(tcell_DEG_all_path, tcell_DEG_tumor_vs_normal_all(sc_tcell), format = "file"),
    tar_target(Tex_TRM_correlations, analyze_CD8_Tex_TRM_correlations(sc_tcell)),
    tar_target(GSEA_Tex_TRM_path, run_GSEA_Tex_TRM(sc_tcell), format = "file"),
    tar_target(sc_adipo, cluster_adipo(sc_final)),
    tar_target(adipo_DEG_plot_path, adipo_DEG_plot(sc_adipo), format = "file"),

    # paper plot
    tar_target(paper_final_annotation_path, paper_final_annotation(sc_final), format = "file"),
    tar_target(paper_clc_expression_by_sample_path, paper_clc_expression_by_sample(sc_final), format = "file"),
    tar_target(paper_myeloid_annotate_path, paper_myeloid_annotate(sc_mye), format = "file"),
    tar_target(paper_myeloid_lipid_DEG_path, paper_myeloid_lipid_DEG(sc_mye), format = "file"),
    tar_target(paper_myeloid_GSEA_path, paper_myeloid_GSEA(sc_mye), format = "file"),
    tar_target(paper_tcell_exhaustion_path, paper_tcell_exhaustion(sc_tcell), format = "file"),
    tar_target(paper_tcell_fate_DEG_path, paper_tcell_fate_DEG(sc_tcell), format = "file"),
    
    # re cluster
    tar_target(cell_types, unique(sc_final_2$cell_type_dtl)),
    tar_target(DEG_cluster_path, DEG_cluster(sc_final_2, cell_types), format = "file", pattern = map(cell_types)),
    tar_target(compare_neural_clc_path, compare_neural_clc(sc_final_2), format = "file")
)
