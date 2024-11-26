cluster_data <- function(sc_int) {
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.2, algorithm = 4, method = "igraph")
    sc_int[["umap_integrated_adj"]] <- CollapseEmbeddingOutliers(
        sc_int,
        reduction = "umap_integrated",
        dims = 1:2,
        group.by = "seurat_clusters",
        outlier.sd = 2
    )
    p <- DimPlot2(sc_int, reduction = "umap_integrated_adj", group.by = "seurat_clusters", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/cluster.png", p, width = 7, height = 7)
    sc_int
}

annotate_data <- function(sc_int, cell_type_path) {
    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2", "HSD3B2", "CYP17A1", "CYP21A2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1"),
        "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
        "Med" = c("TH", "CHGA", "CHGB"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP"),
        "Fib" = c("LUM", "PDGFRA", "ACTA2", "TAGLN", "MGP"),
        "VSM" = c("ACTA2", "MYH11"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1", "THY1"),
        "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R")
    )
    #     marker_genes <- c(
    #     "CD19", "CD20", "MS4A1", "CD79A", "CD79B", # B cell
    #     "CD3D", "CD3E", "CD3G", "CD3", "TRBC2", # T cell
    #     "NKG7", "KLRD1", "KLRF1", "GNLY", "NCR1", # NK cell
    #     "PECAM1", "CD34", "CD31", "KDR", # endothelial cell
    #     "MZB1", "IGHA1", "IGHG1", # plasma cell
    #     "CD14", "CD16", "LYZ", # monocyte
    #     "CD68", "C1QA", # macrophage
    #     "CD1C", "CLEC4C", "IRF8", # dendritic cell
    #     "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R", # Granulocyte
    #     "TPSB2", "TPSAB1", "MS4A2", # Mast cell
    #     "ACTA2", "COL1A1", "COL1A2", "TAGLN", # Fibroblast
    #     "PDZK1IP1", "LRP2", "ALDOB", # proximal tubule cell
    #     "SLC4A4", "SLC5A12", "SLC22A19", # kidney cell
    #     "CA9", "NDUFA4L2", "VEGFA" # RCC cell
    # )
    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "zscore", order = "value"
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p1 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_zscore.png", p1, width = 7, height = 7)

    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "logFC", order = "value", exp.transform = TRUE
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p2 <- Heatmap(toplot, lab_fill = "logFC", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_logFC.png", p2, width = 7, height = 7)

    p3 <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "seurat_clusters", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))

    p <- (p1 + p2) / p3
    ggsave("results/102.cluster_annotate/heatmap.png", p, width = 14, height = 19)

    cell_type <- read_tsv(cell_type_path) %>%
        mutate(cluster = factor(cluster))
    sc_int <- sc_int %>%
        left_join(cell_type, by = c("seurat_clusters" = "cluster"))
    sc_int
}

sub_cluster_steroidogenic <- function(sc_int) {
    sc_int <- sc_int %>%
        filter(cell_type == "Steroidogenic")
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.15, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_steroidogenic")
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "seurat_clusters_steroidogenic", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/sub_cluster_steroidogenic.png", p, width = 7, height = 7)

    sc_int
    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2", "HSD3B2", "CYP17A1", "CYP21A2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1")
    )
    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "zscore", order = "value"
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p1 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/sub_cluster_steroidogenic_heatmap_zscore.png", p1, width = 7, height = 7)

    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "logFC", order = "value", exp.transform = TRUE
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p2 <- Heatmap(toplot, lab_fill = "logFC", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/sub_cluster_steroidogenic_heatmap_logFC.png", p2, width = 7, height = 7)

    p3 <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "seurat_clusters_steroidogenic", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    p <- (p1 + p2) / p3
    ggsave("results/102.cluster_annotate/sub_cluster_steroidogenic_heatmap.png", p, width = 14, height = 19)
    sc_anno_df <- sc_int %>%
        select(cell, seurat_clusters_steroidogenic) %>%
        as_tibble()
    sc_anno_df
}

final_annotation <- function(sc_int, sc_sub_stero) {
    sc_int <- sc_int %>%
        left_join(sc_sub_stero, by = c("cell" = "cell")) %>%
        # change cell_type_dtl based on seurat_clusters_steroidogenic, 1 is ZR, 2 is ZG, 3 is ZF, keep the rest unchanged
        mutate(cell_type_dtl = case_when(
            seurat_clusters_steroidogenic == 1 ~ "ZR",
            seurat_clusters_steroidogenic == 2 ~ "ZG",
            seurat_clusters_steroidogenic == 3 ~ "ZF",
            TRUE ~ cell_type_dtl
        )) %>%
        select(-seurat_clusters_steroidogenic)
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "cell_type_dtl", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/final_annotation.png", p, width = 7, height = 7)
    sc_int
}

save_annotate <- function(sc_final) {
    scCustomize::as.anndata(
        x = sc_final, file_path = "data/104.RNA_velocity", file_name = "anndata.h5ad",
        main_layer = "counts", other_layers = c("data")
    )
    "data/104.RNA_velocity/anndata.h5ad"
}
