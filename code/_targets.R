library(targets)
library(crew)
source("code/R/101.load_data.R")
source("code/R/102.cluster_annotate.R")
source("code/R/103.DEG.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "rliger", "Seurat", "SeuratExtend", "tidyseurat", "scCustomize", "patchwork", "tidyplots"
    ),
    controller = crew_controller_local(workers = 20, seconds_timeout = 6000),
    format = "qs",
    storage = "worker", retrieval = "worker"
)
tar_config_set(
    script = "code/_targets.R",
    store = "data/_targets"
)
options(max.print = 12, spe = "human")
list(
    tar_target(sc_raw, load_sc()),
    tar_target(sc_int, process_sc(sc_raw)),
    tar_target(sc_cluster, cluster_data(sc_int)),
    tar_target(cell_type_path, "results/102.cluster_annotate/cell_type.tsv", format = "file"),
    tar_target(sc_annotate, annotate_data(sc_cluster, cell_type_path)),
    tar_target(sc_sub_stero, sub_cluster_steroidogenic(sc_annotate)),
    tar_target(sc_final, final_annotation(sc_annotate, sc_sub_stero)),
    tar_target(anno_DEG_path, find_DEGs_pseudo_bulk(sc_final), format = "file"),
    # tar_target(anndata, save_annotate(sc_final)),
    tar_target(DEG_path, DEG_annotation(sc_final), format = "file"),
    tar_target(GSEA_path, run_GSEA(sc_final), format = "file"),
    tar_target(distribution_path, plot_distribution(sc_final), format = "file"),
    tar_target(parents, {
        options(max.print = 12, spe = "human")
        GeneSetAnalysisGO() %>% names()
    }),
    tar_target(
        sc_filtered,
        sc_final %>% filter(cell_type == "Fib-Cap-Adipo")
    ),
    tar_target(GSEA_Fib_Cap_Adipo_path, run_GSEA_Fib_Cap_Adipo(sc_filtered, parents),
        format = "file",
        pattern = map(parents)
    ),
    tar_target(GSEA_Fib_Cap_Adipo_hallmark_path, compare_cell_types_hallmark50(sc_filtered),
        format = "file"
    ),
    tar_target(calculate_DEG_mCAF_vs_Others_path, calculate_DEG_mCAF_vs_Others(sc_filtered),
        format = "file"
    )
)
