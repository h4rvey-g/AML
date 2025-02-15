sub_cluster_tcell <- function(sc_final) {
    # 提取T细胞
    sc_int <- sc_final %>%
        filter(cell_type_dtl == "Tcell")

    # 重新聚类
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.8, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster")
    sc_int <- RunUMAP(sc_int,
        dims = 1:30,
        reduction = "harmony",
        reduction.name = "umap_integrated",
        min.dist = 0.05, # Decreased further for tighter clusters
        n.neighbors = 50, # Increased for more connections
        spread = 0.5
    ) # Decreased for more compactness

    # 可视化聚类结果
    p1 <- DimPlot2(sc_int,
        reduction = "umap_integrated",
        group.by = "immune_subcluster", label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    # dir.create("results/108.Tcell", showWarnings = FALSE, recursive = TRUE)
    ggsave("results/108.Tcell/sub_cluster_immune.png", p1, width = 7, height = 7)
    sc_int
}
Tcell_annotation_analysis <- function(sc_tcell_clust) {
    # Extract cluster 5 cells
    cluster5_cells <- subset(sc_tcell_clust, immune_subcluster == "5")

    # Re-cluster cluster 5
    cluster5_cells <- cluster5_cells %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "subcluster")

    # Extract cluster 9 cells
    cluster9_cells <- subset(sc_tcell_clust, immune_subcluster == "9")

    # Re-cluster cluster 9
    cluster9_cells <- cluster9_cells %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "subcluster")

    # Add sub-cluster labels back to main object
    sc_tcell_clust$immune_subcluster <- as.character(sc_tcell_clust$immune_subcluster)

    # Update cluster 5 labels
    new_labels_5 <- paste0("5-", cluster5_cells$subcluster)
    sc_tcell_clust$immune_subcluster[colnames(cluster5_cells)] <- new_labels_5

    # Update cluster 9 labels
    new_labels_9 <- paste0("9-", cluster9_cells$subcluster)
    sc_tcell_clust$immune_subcluster[colnames(cluster9_cells)] <- new_labels_9

    sc_tcell_clust$immune_subcluster <- as.factor(sc_tcell_clust$immune_subcluster)

    # Merge clusters 2 and 3
    sc_tcell_clust$immune_subcluster <- as.character(sc_tcell_clust$immune_subcluster)
    sc_tcell_clust$immune_subcluster[sc_tcell_clust$immune_subcluster == "3"] <- "2"
    sc_tcell_clust$immune_subcluster <- as.factor(sc_tcell_clust$immune_subcluster)
    # Create a mapping from old to new cluster numbers
    current_levels <- sort(unique(as.character(sc_tcell_clust$immune_subcluster)))
    new_levels <- as.character(seq_along(current_levels))
    names(new_levels) <- current_levels

    # Apply the mapping to reorder clusters
    sc_tcell_clust$immune_subcluster <- plyr::mapvalues(sc_tcell_clust$immune_subcluster,
        from = names(new_levels),
        to = new_levels
    )
    sc_tcell_clust$immune_subcluster <- factor(sc_tcell_clust$immune_subcluster, levels = new_levels)
    markers <- list(
        # Core T Cell Markers
        "T_Cells" = c("CD3D", "CD3E", "CD3G"),
        "CD4_T_Cells" = c("CD4"),
        "CD8_T_Cells" = c("CD8A", "CD8B"),

        # Specialized T Cell Subsets
        "Tregs" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
        "gdT_Cells" = c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2"),
        "NK_Cells" = c("NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A", "NCR1", "NCR3", "KIR2DL1", "KIR2DL3", "KIR3DL1"),

        # Naive and Memory Markers
        "Naive_T" = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "SKAP1", "THEMIS", "PTPRC"),
        "Memory_T" = c("CD44", "IL7R", "CD45RO", "PRKCQ", "STAT4"),

        # Tissue-Resident Memory
        "Trm" = c("ITGAE", "CD69", "CXCR6"),

        # Effector Functions
        "Effector" = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "IFNG", "TNF", "FASLG"),

        # Exhaustion Markers
        "Exhausted" = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "ENTPD1"),

        # T Helper Subtypes and Key Transcription Factors
        "Th1" = c("IFNG", "CXCR3", "TBX21"), # TBX21 (T-bet)
        "Th2" = c("IL4", "IL5", "IL13", "CCR4", "GATA3"),
        "Th17" = c("IL17A", "IL17F", "IL22", "CCR6", "RORC")
    )


    # 计算并绘制热图
    toplot <- CalcStats(sc_tcell_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap.png", p2, width = 8, height = 15)
    write_csv(
        toplot %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/108.Tcell/immune_subcluster_expression.csv"
    )

    sc_tcell_clust
    # # 按照数据集和亚群聚合进行DEG分析
    # sc_pseudo <- AggregateExpression(sc_tcell_clust,
    #     assays = "RNA",
    #     return.seurat = TRUE,
    #     group.by = c("dataset", "immune_subcluster")
    # )

    # # 设置比较分组
    # Idents(sc_pseudo) <- "immune_subcluster"
    # clusters <- unique(sc_tcell_clust$immune_subcluster)
    # dir.create("results/108.Tcell/cluster_DEG", showWarnings = FALSE, recursive = TRUE)

    # # # Setup parallel processing
    # # library(foreach)
    # # library(doParallel)
    # num_cores <- parallel::detectCores() - 1
    # registerDoParallel(cores = num_cores)

    # # Compare each cluster with all other clusters in parallel
    # foreach(cluster = clusters, .packages = c("Seurat", "tidyverse", "MAST")) %dopar% {
    #     # DESeq2 analysis
    #     deg <- FindMarkers(
    #         sc_pseudo,
    #         ident.1 = cluster,
    #         test.use = "DESeq2",
    #         min.cells.group = 2
    #     ) %>%
    #         rownames_to_column("gene") %>%
    #         as_tibble() %>%
    #         filter(p_val_adj < 0.05) %>%
    #         filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
    #         arrange(desc(abs(avg_log2FC)))
    #     write_tsv(deg, sprintf(
    #         "results/108.Tcell/cluster_DEG/%s_vs_rest_DEG_psuedo.tsv",
    #         cluster
    #     ))

    #     # MAST analysis
    #     deg <- FindMarkers(
    #         sc_tcell_clust,
    #         ident.1 = cluster,
    #         test.use = "MAST",
    #         min.cells.group = 2
    #     ) %>%
    #         rownames_to_column("gene") %>%
    #         as_tibble() %>%
    #         filter(p_val_adj < 0.05) %>%
    #         filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
    #         arrange(desc(abs(avg_log2FC)))

    #     write_tsv(deg, sprintf(
    #         "results/108.Tcell/cluster_DEG/%s_vs_rest_DEG_MAST.tsv",
    #         cluster
    #     ))
    # }

    # # Stop parallel processing
    # stopImplicitCluster()
}

Tcell_annotate <- function(sc_tcell_clust) {
    cluster_map <- c(
        "1" = "CD4+ Naive/CM T", # Naive or Central Memory T cells
        "2" = "CD69+ Treg", # CD69+ Tregs
        "3" = "CD4+ Naive/CM T", # Naive or Central Memory T cells
        "4" = "CD4+ Th17/Th2", # Th17 or Th2 cells
        "5" = "CD8+ NK-like CTL", # Cytotoxic T cells
        "6" = "CD8+ NK-like gdT", # gamma delta T cells
        "7" = "Unknown", # Unknown T cell type
        "8" = "CD8+ Mem T", # Memory T cells
        "9" = "Th17 ex", # Exhausted Th17 cells
        "10" = "CD4+ Th17/Th2", # Th17 or Th2 cells
        "11" = "CD8+ FOXP+ gdT" # Tissue-resident CD8+ gamma delta T cells
    )


    # Add new annotation column
    sc_tcell_clust <- sc_tcell_clust %>%
        mutate(cell_type_dtl = case_when(
            immune_subcluster %in% names(cluster_map) ~ cluster_map[immune_subcluster],
            TRUE ~ as.character(immune_subcluster)
        ))

    # Run UMAP
    sc_tcell_clust <- RunUMAP(sc_tcell_clust,
        dims = 1:30, reduction = "harmony", reduction.name = "umap_tcell",
        min.dist = 0.05,
        spread = 0.3
    )

    # Visualize the new annotations
    p <- DimPlot2(sc_tcell_clust,
        reduction = "umap_tcell",
        group.by = "cell_type_dtl",
        label = TRUE,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))

    # Save the plot
    ggsave("results/108.Tcell/tcell_annotated.png", p, width = 8, height = 7)
    sc_tcell_clust
}
