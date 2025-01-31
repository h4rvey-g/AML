sub_cluster_myeloid <- function(sc_int) {
    # 提取巨噬细胞和粒细胞
    sc_int <- sc_int %>%
        filter(cell_type_dtl == "Myeloid")

    # 重新聚类
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.4, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster")

    # 可视化聚类结果
    p1 <- DimPlot2(sc_int,
        reduction = "umap_integrated",
        group.by = "immune_subcluster", label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    dir.create("results/107.myeloid", showWarnings = FALSE, recursive = TRUE)
    ggsave("results/107.myeloid/sub_cluster_immune.png", p1, width = 7, height = 7)

    markers <- list(
        # Monocytes
        "Classical_Mono" = c("CD14", "VCAN", "CCR2", "LYZ"),
        "NonClassical_Mono" = c("FCGR3A", "CX3CR1", "LST1", "MS4A7"),

        # Macrophages
        "Mac_M1" = c("CD80", "CD86", "STAT1", "IL1B", "CXCL10", "HLA-DRA", "NOS2", "IL6", "IL12B"),
        "Mac_M2" = c("CD163", "MRC1", "FOLR2", "CCL18", "MAFB"),
        "Mac_LAM" = c("TREM2", "APOE", "LPL", "FABP4"),

        # Granulocytes
        "Neutrophil" = c("CSF3R", "CEACAM8", "CXCR2", "LCN2"),
        "Eosinophil" = c("CLC", "PRG2", "IL5RA", "EPX", "ALOX15"),
        "Mast" = c("KIT", "TPSAB1", "CMA1", "CPA3"),

        # Activation/Proliferation
        "Activation" = c("CD83", "CD40", "LAMP1"),
        "Proliferation" = c("MKI67", "TOP2A", "CENPF"),

        # Tissue-Resident Macrophages
        "Microglia" = c("P2RY12", "TMEM119"),
        "Kupffer" = c("CLEC4F", "ID1"),
        "Alveolar_Mac" = c("MARCO", "PPARG")
    )

    # 计算并绘制热图
    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "zscore", order = "value"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/sub_cluster_immune_heatmap.png", p2, width = 8, height = 15)
    write_tsv(
        toplot %>% as.data.frame() %>% rownames_to_column("gene"),
        "results/107.myeloid/immune_subcluster_expression.tsv"
    )
    # 创建DEG结果目录
    dir.create("results/107.myeloid/immune_DEG", showWarnings = FALSE, recursive = TRUE)

    # 按照数据集和亚群聚合进行DEG分析
    sc_pseudo <- AggregateExpression(sc_int,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c("dataset", "immune_subcluster")
    )

    # 设置比较分组
    Idents(sc_pseudo) <- "immune_subcluster"
    clusters <- unique(sc_int$immune_subcluster)
    dir.create("results/107.myeloid/immune_DEG", showWarnings = FALSE, recursive = TRUE)

    # Compare each cluster with all other clusters
    for (cluster in clusters) {
        deg <- FindMarkers(
            sc_pseudo,
            ident.1 = cluster,
            test.use = "DESeq2",
            min.cells.group = 2
        ) %>%
            rownames_to_column("gene") %>%
            as_tibble() %>%
            filter(p_val_adj < 0.05) %>%
            arrange(desc(abs(avg_log2FC)))

        write_tsv(deg, sprintf(
            "results/107.myeloid/immune_DEG/%s_vs_rest_DEG.tsv",
            cluster
        ))
    }

    # Compare clusters pairwise
    for (i in 1:(length(clusters) - 1)) {
        for (j in (i + 1):length(clusters)) {
            cluster1 <- clusters[i]
            cluster2 <- clusters[j]

            deg <- FindMarkers(
                sc_pseudo,
                ident.1 = cluster1,
                ident.2 = cluster2,
                test.use = "DESeq2",
                min.cells.group = 2
            ) %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                filter(p_val_adj < 0.05) %>%
                arrange(desc(abs(avg_log2FC)))

            write_tsv(deg, sprintf(
                "results/107.myeloid/immune_DEG/%s_vs_%s_DEG.tsv",
                cluster1, cluster2
            ))
        }
    }
    sc_int
}

myeloid_annotate <- function(sc_mye_clust) {
    # Create mapping of cluster numbers to cell type annotations
    cluster_map <- c(
        "1" = "TREM2_Mac",
        "2" = "FOLR2_Mac",
        "3" = "M2_Mac",
        "4" = "Mast",
        "5" = "Eosinophil",
        "6" = "Foam_Cell"
    )

    # Add new annotation column
    sc_mye_clust <- sc_mye_clust %>%
        # mutate cell_type_dtl = related cell type based on immune_subcluster
        mutate(cell_type_dtl = case_when(
            immune_subcluster %in% names(cluster_map) ~ cluster_map[immune_subcluster],
            TRUE ~ as.character(immune_subcluster)
        ))
    # run umap
    sc_mye_clust <- RunUMAP(sc_mye_clust, dims = 1:30, reduction = "harmony", reduction.name = "umap_mye")
    # Visualize the new annotations
    p <- DimPlot2(sc_mye_clust,
        reduction = "umap_mye",
        group.by = "cell_type_dtl",
        label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))

    # Save the plot
    ggsave("results/107.myeloid/myeloid_annotated.png", p, width = 8, height = 7)

    return(sc_mye_clust)
}

plot_myeloid_distribution <- function(sc_mye) {
    dir.create("results/107.myeloid/distribution", showWarnings = FALSE, recursive = TRUE)
    
    df1 <- ClusterDistrBar(origin = sc_mye$dataset, cluster = sc_mye$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")
    
    p <- tidyplot(
        df1 %>%
            mutate(group = ifelse(str_starts(origin, "N"), "N", "T")),
        x = cluster, y = percentage, color = group
    ) %>%
        add_mean_bar() %>%
        add_sem_errorbar() %>%
        adjust_colors(c(colors_discrete_metro, colors_discrete_seaside)) %>%
        adjust_size(width = 170, height = 100)
    
    save_plot(p, "results/107.myeloid/distribution/myeloid_distribution.png")
    writexl::write_xlsx(list(df1 = df1), "results/107.myeloid/distribution/myeloid_distribution.xlsx")
    
    return("results/107.myeloid/distribution")
}
