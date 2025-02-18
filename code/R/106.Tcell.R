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
        "5" = "CD8+ NKL CTL", # Cytotoxic T cells
        "6" = "CD8+ NKL gdT", # gamma delta T cells
        "7" = "Unknown", # Unknown T cell type
        "8" = "CD8+ Mem T", # Memory T cells
        "9" = "Th17 ex", # Exhausted Th17 cells
        "10" = "CD4+ Th17/Th2", # Th17 or Th2 cells
        "11" = "CD8+ FOXP+ gdT" # CD8+ FOXP+ gamma delta T cells
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
    ggsave("results/108.Tcell/tcell_annotated.png", p, width = 9, height = 6)
    markers <- list(
        "CD4/8+" = c("CD4", "CD8A", "CD8B"),
        "CD4+ Naive/CM T" = c("CD4", "CCR7", "TCF7", "LEF1", "SELL"),
        "CD69+ Treg" = c("CD69", "NCAM1", "IFNG"),
        "CD4+ Th17/Th2" = c("RORC", "GATA3", "IL17F", "TCF7", "IL7R"),
        "CD8+ NKL CTL" = c("GZMB", "PRF1", "NKG7"),
        "CD8+ NKL gdT" = c("TRDV1", "NCR1"),
        "Unknown" = c(NA, NA, NA, NA, NA), # Placeholder - needs further analysis
        "CD8+ Mem T" = c("SKAP1", "THEMIS", "CD44"),
        "Th17 ex" = c("IL17A", "HAVCR2", "ENTPD1", "ITGAE"),
        "CD8+ FOXP+ gdT" = c("FOXP3", "TRDC", "TRDV2", "CD8A")
    )
    # 计算并绘制热图
    toplot <- CalcStats(sc_tcell_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_heatmap.png", p2, width = 15, height = 8)

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
        group.by = "cell_type_dtl"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]
    write_csv(
        toplot %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/108.Tcell/tcell_annotated_expression.csv"
    )
    sc_tcell_clust
}

test_tcell_composition <- function(sc_tcell) {
    # Create directory
    dir.create("results/108.Tcell/distribution", showWarnings = FALSE, recursive = TRUE)

    # Create count data frame for sccomp with proper dataset column
    # tcell_counts <- sc_tcell@meta.data %>%
    #     select(cell_type_dtl, dataset, group) %>%
    #     mutate(
    #         count = 1
    #     ) %>%
    #     group_by(cell_type_dtl, dataset, group) %>%
    #     summarize(count = n(), .groups = "drop")

    # Run composition test with updated parameters
    sc_tcell <- sc_tcell %>%
        AddMetaData(
            object = .,
            metadata = .$dataset,
            col.name = "sample"
        ) %>%
        AddMetaData(
            object = .,
            metadata = .$cell_type_dtl,
            col.name = "cell_group"
        )
    composition_test <- sc_tcell %>%
        sccomp_estimate(
            formula_composition = ~group,
            # formula_variability = ~1,
            .sample = sample,
            .cell_group = cell_group,
            bimodal_mean_variability_association = TRUE,
            # .abundance = count, # Updated from .count to .abundance
            cores = 20,
            verbose = TRUE
        ) %>%
        sccomp_remove_outliers(cores = 31, verbose = FALSE) %>%
        sccomp_test()


    # Generate plots
    p1 <- composition_test %>%
        sccomp_boxplot(factor = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save fold changes
    fold_changes <- composition_test %>%
        sccomp_proportional_fold_change(
            formula_composition = ~group,
            from = "normal",
            to = "tumor"
        ) %>%
        select(cell_group, statement)

    write_tsv(fold_changes, "results/108.Tcell/distribution/fold_changes.tsv")

    # Save test results
    write_tsv(
        composition_test %>%
            select(-count_data),
        "results/108.Tcell/distribution/test_results.tsv"
    )

    composition_test
}

plot_tcell_distribution <- function(sc_tcell, t_composition_test) {
    dir.create("results/108.Tcell/distribution", showWarnings = FALSE, recursive = TRUE)

    # Calculate distribution percentages
    df1 <- ClusterDistrBar(origin = sc_tcell$dataset, cluster = sc_tcell$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")

    # Calculate cell counts
    cell_counts <- sc_tcell@meta.data %>%
        group_by(cell_type_dtl, dataset) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(cell_type_dtl) %>%
        summarise(
            N_count = sum(count[str_starts(dataset, "N")]),
            T_count = sum(count[str_starts(dataset, "T")]),
            .groups = "drop"
        )

    # Add cell counts to cluster labels
    df1 <- df1 %>%
        left_join(cell_counts, by = c("cluster" = "cell_type_dtl")) %>%
        mutate(cluster = sprintf("%s\n(N=%d,T=%d)", cluster, N_count, T_count))

    t_composition_test_res <- t_composition_test %>%
        filter(factor == "group", c_FDR < 0.05) %>%
        # convert c_FDR to stars
        mutate(stars = case_when(
            c_FDR < 0.001 ~ "***",
            c_FDR < 0.01 ~ "**",
            c_FDR < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        select(cell_group, stars)

    # Add stars to cluster labels based on significance
    df1 <- df1 %>%
        mutate(cluster = map_chr(str_extract(cluster, ".*(?=\\n)"), ~ {
            stars <- t_composition_test_res$stars[t_composition_test_res$cell_group == .x]
            count_info <- sprintf("(N=%d,T=%d)", cell_counts$N_count[cell_counts$cell_type_dtl == .x], cell_counts$T_count[cell_counts$cell_type_dtl == .x])
            if (length(stars) > 0 && stars != "") {
                paste0(.x, "\n", count_info, "\n", stars)
            } else {
                paste0(.x, "\n", count_info, "\n")
            }
        }))

    p <- tidyplot(
        df1 %>%
            mutate(group = ifelse(str_starts(origin, "N"), "N", "T")),
        x = cluster, y = percentage, color = group
    ) %>%
        add_mean_bar() %>%
        add_sem_errorbar() %>%
        adjust_colors(c(colors_discrete_metro, colors_discrete_seaside)) %>%
        adjust_size(width = 170, height = 100) %>%
        adjust_x_axis(rotate_labels = 45)

    tidyplots::save_plot(p, "results/108.Tcell/distribution/tcell_distribution.png")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/108.Tcell/distribution/tcell_distribution.tsv"
    )

    return("results/108.Tcell/distribution")
}

run_tcell_GSEA <- function(sc_tcell) {
    options(max.print = 12, spe = "human")
    dir.create("results/108.Tcell/GSEA", showWarnings = FALSE, recursive = TRUE)

    # Run GO analysis
    parents <- GeneSetAnalysisGO() %>% names()
    for (parent in parents) {
        dir.create(paste0("results/108.Tcell/GSEA/", parent), showWarnings = FALSE, recursive = TRUE)
        cat("Running GSEA for", parent, "\n")
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 50, parent = parent, n.min = 3)
        matr <- sc_tcell@misc$AUCell$GO[[parent]]
        matr <- RenameGO(matr, add_id = FALSE)

        # Overall heatmap for all subclusters
        toplot <- CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10)
        p <- Heatmap(
            toplot,
            lab_fill = "zscore", text.size = 5
        )
        ggsave(paste0("results/108.Tcell/GSEA/", parent, "/Heatmap_GO_", parent, ".png"), p, width = 14, height = nrow(toplot) * 0.4)
        write_tsv(
            toplot %>% as.data.frame() %>% rownames_to_column("pathway"),
            paste0("results/108.Tcell/GSEA/", parent, "/Heatmap_GO_", parent, ".tsv")
        )

        # Compare tumor vs normal within each subcluster
        subclusters <- unique(sc_tcell$cell_type_dtl)
        for (subcluster in subclusters) {
            sc_temp <- sc_tcell %>%
                filter(cell_type_dtl == subcluster)
            Idents(sc_temp) <- "group"

            if (length(unique(sc_temp$group)) == 2) {
                p <- WaterfallPlot(
                    matr,
                    f = sc_temp$group,
                    ident.1 = "tumor",
                    ident.2 = "normal",
                    top.n = 20,
                    color = "p",
                    length = "logFC",
                    title = paste0(parent, ": ", subcluster)
                )
                ggsave(paste0(
                    "results/108.Tcell/GSEA/", parent, "/Waterfall_GO_",
                    parent, "_", gsub("/", "_", subcluster), "_tumor_vs_normal.png"
                ), p, width = 15, height = 15)
            }
        }
    }

    # Run hallmark analysis
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = hall50$human, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for hallmark pathways
    dir.create("results/108.Tcell/GSEA/hallmark50", showWarnings = FALSE, recursive = TRUE)
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10),
        lab_fill = "zscore", text.size = 5
    )
    ggsave("results/108.Tcell/GSEA/hallmark50/Heatmap_hallmark50.png", p, width = 14, height = 14)
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA/hallmark50/Heatmap_hallmark50.tsv"
    )

    # Run for each GO root category
    for (root in c("BP", "MF", "CC")) {
        cat("Processing", root, "...\n")
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 50, root = root, n.min = 3)
        matr <- sc_tcell@misc$AUCell$GO[[root]]

        top30_pathways <- data.frame()
        for (cell_type in unique(sc_tcell$cell_type_dtl)) {
            stats <- CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore") %>%
                as.data.frame() %>%
                dplyr::select(all_of(cell_type)) %>%
                arrange(desc(.[[1]])) %>%
                head(30) %>%
                rownames_to_column("pathway")

            top30_pathways <- bind_rows(top30_pathways, stats)
        }

        top30_pathways <- top30_pathways %>%
            pull(pathway) %>%
            unique()

        matr_filt <- matr[top30_pathways, ]
        matr_filt <- RenameGO(matr_filt, add_id = FALSE)
        rownames(matr_filt) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr_filt)))
        toplot <- CalcStats(matr_filt, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 30)

        dir.create(file.path("results/108.Tcell/GSEA/GO", root), showWarnings = FALSE, recursive = TRUE)

        p <- Heatmap(
            CalcStats(matr_filt, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10),
            lab_fill = "zscore", text.size = 5
        )
        ggsave(file.path("results/108.Tcell/GSEA/GO", root, "Heatmap_top10.png"), p, width = 14, height = 14)
        write_tsv(
            toplot %>%
                as.data.frame() %>% rownames_to_column("pathway"),
            file.path("results/108.Tcell/GSEA/GO", root, "Heatmap_top10.tsv")
        )
    }

    return("results/108.Tcell/GSEA")
}

tcell_DEG_tumor_vs_normal <- function(sc_tcell) {
    sc_pseudo <- AggregateExpression(sc_tcell, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    cell_types <- unique(sc_pseudo$cell_type_dtl)

    dir.create("results/108.Tcell/DEG_tumor_vs_normal", showWarnings = FALSE, recursive = TRUE)

    # DESeq2 analysis
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_pseudo) && ident2 %in% Idents(sc_pseudo))) {
            next
        }
        deg <- tryCatch(
            {
                FindMarkers(
                    sc_pseudo,
                    ident.1 = ident1,
                    ident.2 = ident2,
                    test.use = "DESeq2",
                    min.cells.group = 2
                ) %>%
                    Add_Pct_Diff() %>%
                    rownames_to_column("gene") %>%
                    as_tibble() %>%
                    filter(p_val_adj < 0.05) %>%
                    filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                    arrange(desc(abs(avg_log2FC)))
            },
            error = function(e) {
                message("Error in FindMarkers for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )

        if (!is.null(deg)) {
            write_tsv(
                deg %>%
                    mutate(across(where(is.numeric), ~ ifelse(abs(.) < 0.001, signif(., 3), round(., 3)))),
                paste0("results/108.Tcell/DEG_tumor_vs_normal/", gsub("/", "_", cell_type), "_DEG.tsv")
            )
        }
    }

    # MAST analysis with parallel processing
    sc_tcell$cell_type_group <- paste(sc_tcell$cell_type_dtl, sc_tcell$group, sep = "_")
    Idents(sc_tcell) <- "cell_type_group"

    cores <- min(parallel::detectCores() - 1, length(cell_types))
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    foreach(cell_type = cell_types, .packages = c("Seurat", "tidyverse", "dplyr")) %dopar% {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if ((ident1 %in% Idents(sc_tcell) && ident2 %in% Idents(sc_tcell))) {
            deg <- FindMarkers(
                sc_tcell,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                min.cells.group = 2
            ) %>%
                scCustomize::Add_Pct_Diff() %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                filter(p_val_adj < 0.05) %>%
                filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                arrange(desc(abs(avg_log2FC)))

            write_tsv(
                deg %>%
                    mutate(across(where(is.numeric), ~ ifelse(abs(.) < 0.001, signif(., 3), round(., 3)))),
                paste0("results/108.Tcell/DEG_tumor_vs_normal/", gsub("/", "_", cell_type), "_DEG_MAST.tsv")
            )
        }
    }

    parallel::stopCluster(cl)
    return("results/108.Tcell/DEG_tumor_vs_normal")
}
