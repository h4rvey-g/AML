cluster_data <- function(sc_int) {
    sc_int <- FindClusters(sc_int, resolution = 0.3)
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = c("seurat_clusters", "group")) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/umap_integrated.png", p, width = 14, height = 7)
    sc_int
}