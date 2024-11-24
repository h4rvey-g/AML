library(targets)
library(crew)
source("code/R/101.load_data.R")
source("code/R/102.cluster_annotate.R")
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "rliger", "Seurat", "SeuratExtend", "tidyseurat", "scCustomize"
    ),
    controller = crew_controller_local(workers = 20, seconds_timeout = 6000),
    format = "qs",
    storage = "worker", retrieval = "worker"
)
tar_config_set(
    script = "code/_targets.R",
    store = "data/_targets"
)
list(
    tar_target(sc_raw, load_sc()),
    tar_target(sc_int, process_sc(sc_raw)),
    tar_target(sc_cluster, cluster_data(sc_int)),
    tar_target(cell_type_path, "results/102.cluster_annotate/cell_type.tsv", format = "file"),
    tar_target(sc_annotate, annotate_data(sc_cluster, cell_type_path))
    # tar_target(anndata, save_annotate(sc_annotate))
)
