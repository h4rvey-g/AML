cluster_data <- function(sc_int) {
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.3, algorithm = 4, method = "igraph")
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
        # "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
        "Med" = c("TH", "CHGA", "CHGB", "KIT"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
        "Fib" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "mCAFs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "ENG", "MCAM"),
        "ASC_preadipocyte" = c("PPARG", "LPL", "CD36", "PDGFRA", "PDGFRB", "LEP"),
        # "Fib_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
        # "Fib_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
        "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
        "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
    )
    # markers <- list(
    #     "Fib_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
    #     "Ret_Fib" = c("THY1", "FMO1", "MYOC", "LSP1",  "ACTA2", "PPARG", "CD36"),
    #     "mCAFs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC"),
    #     "iCAFs" = c("MMP1", "MMP3", "IL6", "CXCL8", "IDO1"),
    #     "RGS5_positive" = c("RGS5", "KCNJ8", "MCAM", "TAGLN"),
    #     "vSMCs" = c("RERGL", "MYH11", "CNN1", "MCAM")
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
        method = "tscore", order = "value", exp.transform = TRUE
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p2 <- Heatmap(toplot, lab_fill = "tscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_tscore.png", p2, width = 7, height = 7)

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
        method = "mean", order = "value", exp.transform = TRUE
    )
    # see which is in markers but not in rownames(toplot)
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    # remove missing genes in each element of markers
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    # order gene_groups by rownames(toplot)
    gene_groups <- gene_groups[rownames(toplot)]
    p2 <- Heatmap(toplot, lab_fill = "mean", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/sub_cluster_steroidogenic_heatmap_mean.png", p2, width = 7, height = 7)

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
            seurat_clusters_steroidogenic == 3 ~ "ZR",
            seurat_clusters_steroidogenic %in% c(2, 4) ~ "ZG",
            seurat_clusters_steroidogenic %in% c(1, 5) ~ "ZF",
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

find_DEGs_pseudo_bulk <- function(sc_final) {
    # 聚合表达数据
    # sc_pseudo <- AggregateExpression(sc_final, assays = "RNA", return.seurat = TRUE, group.by = c("group", "dataset", "cell_type_dtl"))
    # Idents(sc_pseudo) <- "cell_type_dtl"
    cell_types <- c("11", "13")

    # 创建目录
    dir.create("results/102.cluster_annotate/cell_types_DEG/normal", showWarnings = FALSE, recursive = TRUE)
    dir.create("results/102.cluster_annotate/cell_types_DEG/tumor", showWarnings = FALSE, recursive = TRUE)

    for (group_1 in c("normal", "tumor")) {
        sc_final_sub <- sc_final %>%
            filter(group == group_1)
        sc_pseudo <- AggregateExpression(sc_final_sub, assays = "RNA", return.seurat = TRUE, group.by = c("dataset", "seurat_clusters"))
        Idents(sc_pseudo) <- "seurat_clusters"
        cell_type <- ifelse(group_1 == "normal", "11", "13")
        deg <- FindMarkers(
            sc_pseudo,
            ident.1 = cell_type,
            test.use = "DESeq2",
            min.cells.group = 2
        ) %>%
            # Add_Pct_Diff() %>%
            rownames_to_column("gene") %>%
            as_tibble() %>%
            filter(p_val_adj < 0.05) %>%
            arrange(desc(abs(avg_log2FC)))
        write_tsv(deg, paste0("results/102.cluster_annotate/cell_types_DEG/", group_1, "/DEG_", gsub("/", "_", cell_type), ".tsv"))
    }
    sc_pseudo <- AggregateExpression(sc_final, assays = "RNA", return.seurat = TRUE, group.by = c("dataset", "cell_type_dtl"))
    Idents(sc_pseudo) <- "cell_type_dtl"
    deg <- FindMarkers(
        sc_pseudo,
        ident.1 = "mCAF",
        ident.2 = "MSC",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        # Add_Pct_Diff() %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        arrange(desc(abs(avg_log2FC)))
    write_tsv(deg, paste0("results/102.cluster_annotate/cell_types_DEG/Fib_Cap_DEG.tsv"))
}

# 示例调用
# find_DEGs_pseudo_bulk(sc_final)
