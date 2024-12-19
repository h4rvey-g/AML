run_palantir <- function(sc_final) {
    library(SeuratExtend)
    library(Seurat)
    library(tidyseurat)
    tar_load(sc_final)
    sc_final <- sc_final %>%
        filter(cell_type %in% c("Fib-Cap-Adipo"))
    sc_final <- Palantir.RunDM(sc_final, conda_env = "base")
    # create a var, combine cell_type_dtl and group
    sc_final <- sc_final %>%
        mutate(cell_type_group = paste(cell_type_dtl, group, sep = "_"))
    p <- DimPlot2(sc_final, reduction = "ms", group.by = "cell_type_group", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/palantir.png", p, width = 8, height = 6)
}
