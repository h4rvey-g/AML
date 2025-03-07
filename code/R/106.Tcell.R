sub_cluster_tcell <- function(sc_final) {
    # 提取T细胞
    sc_int <- sc_final %>%
        filter(cell_type_dtl == "Tcell")

    # 重新聚类
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.9, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster")
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
Tcell_annotation_analysis <- function(sc_tcell_clust, sc_final) {
    # # Extract cluster 5 cells
    # cluster5_cells <- subset(sc_tcell_clust, immune_subcluster == "6")

    # # Re-cluster cluster 5
    # cluster5_cells <- cluster5_cells %>%
    #     FindNeighbors(dims = 1:30, reduction = "harmony") %>%
    #     FindClusters(resolution = 0.9, algorithm = 4, method = "igraph", cluster.name = "subcluster")

    # # Extract cluster 9 cells
    # cluster9_cells <- subset(sc_tcell_clust, immune_subcluster == "9")

    # # Re-cluster cluster 9
    # cluster9_cells <- cluster9_cells %>%
    #     FindNeighbors(dims = 1:30, reduction = "harmony") %>%
    #     FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "subcluster")

    # # Add sub-cluster labels back to main object
    # sc_tcell_clust$immune_subcluster <- as.character(sc_tcell_clust$immune_subcluster)

    # # Update cluster 5 labels
    # new_labels_5 <- paste0("6-", cluster5_cells$subcluster)
    # sc_tcell_clust$immune_subcluster[colnames(cluster5_cells)] <- new_labels_5

    # # # Update cluster 9 labels
    # # new_labels_9 <- paste0("9-", cluster9_cells$subcluster)
    # # sc_tcell_clust$immune_subcluster[colnames(cluster9_cells)] <- new_labels_9

    # sc_tcell_clust$immune_subcluster <- as.factor(sc_tcell_clust$immune_subcluster)

    # Merge clusters 2 and 3
    sc_tcell_clust$immune_subcluster <- as.character(sc_tcell_clust$immune_subcluster)
    sc_tcell_clust$immune_subcluster[sc_tcell_clust$immune_subcluster == "2"] <- "1"
    sc_tcell_clust$immune_subcluster <- as.factor(sc_tcell_clust$immune_subcluster)
    # # Create a mapping from old to new cluster numbers
    # current_levels <- sort(unique(as.character(sc_tcell_clust$immune_subcluster)))
    # new_levels <- as.character(seq_along(current_levels))
    # names(new_levels) <- current_levels

    # # Apply the mapping to reorder clusters
    # sc_tcell_clust$immune_subcluster <- plyr::mapvalues(sc_tcell_clust$immune_subcluster,
    #     from = names(new_levels),
    #     to = new_levels
    # )
    # sc_tcell_clust$immune_subcluster <- factor(sc_tcell_clust$immune_subcluster, levels = new_levels)
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
    # 计算并绘制热图
    toplot <- CalcStats(sc_tcell_clust %>% filter(group == "normal"),
        features = markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_normal.png", p2, width = 8, height = 15)
    write_csv(
        toplot %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/108.Tcell/immune_subcluster_expression_normal.csv"
    )

    sc_tcell_clust
}

Tcell_annotate <- function(sc_tcell_clust) {
    cluster_map <- c(
        "1" = "CD4_Th17",
        # "2" = "CD4_Treg",
        "3" = "NKL_resting_T",
        "4" = "CD8_Trm",
        "5" = "CD4_TN", # Naive CD4+ T cells
        "6" = "CD8_IFNG_gdT", # γδ T cells producing IFNγ
        "7" = "CTLA4_Act_T (early)",
        "8" = "CD8_Th2L_gdT" # γδ T cells with Th2-like features
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
        # Core T cell subtype markers
        "CD4/8" = c("CD4", "CD8A", "CD8B"),

        # Cluster 1: CD4+ Th17 cells
        "CD4_Th17" = c("RORC", "CCR6", "IL17F"),

        # Cluster 5: Naive CD4+ T cells
        "CD4_TN" = c("CCR7", "LEF1", "SELL"),
        "Act_T (early)" = c("CTLA4"),
        "NKL/Innate_T" = c("NCR3", "FCGR3A"),
        # Cluster
        "Trm" = c("CD69", "CXCR3", "ITGAE"),

        # Cluster 6: γδ T cells producing IFNγ
        "CD8_IFNG_gdT" = c("TRDV1", "IFNG"),

        # Cluster 8: γδ T cells with Th2-like features
        "CD8_Th2L_gdT" = c("TRDV2", "IL13", "IL4")
        # Cluster 3: Ambiguous T cells

        # # Cluster 4: CD8+ Regulatory T cells
        # "CD8_Treg" = c("FOXP3", "CD69", "THEMIS", "CXCR3"),
    )
    # 计算并绘制热图
    toplot <- CalcStats(sc_tcell_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )
    toplot <- toplot[, c("CD4_Th17", "CD4_TN", "CTLA4_Act_T (early)", "NKL_resting_T", "CD8_Trm", "CD8_IFNG_gdT", "CD8_Th2L_gdT")]

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p2 <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_heatmap.png", p2, width = 15, height = 8)

    p <- DotPlot2(sc_tcell_clust,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_dotplot.png", p, width = 15, height = 8)
    # 计算并绘制热图
    toplot <- CalcStats(sc_tcell_clust %>% filter(group == "normal"),
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )
    toplot <- toplot[, c("CD4_Th17", "CD4_TN", "CTLA4_Act_T (early)", "NKL_resting_T", "CD8_Trm", "CD8_IFNG_gdT", "CD8_Th2L_gdT")]

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p2 <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_normal_heatmap.png", p2, width = 15, height = 8)
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
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 20, parent = parent, n.min = 3)
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
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 20, root = root, n.min = 3)
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

analyze_CD8_Tex_TRM_correlations <- function(sc_tcell) {
    # 创建结果目录
    dir.create("results/108.Tcell/correlation", showWarnings = FALSE, recursive = TRUE)

    # 提取CD8_Tex_TRM细胞
    sc_tex <- subset(sc_tcell, cell_type_dtl == "CD8_Tex_TRM")

    # 定义要分析的基因列表
    gene_sets <- list(
        "Adenosine_Receptors" = c("ADORA1", "ADORA2A", "ADORA2B", "ADORA3"),
        "Cytotoxicity" = c("IFNG", "GZMB", "PRF1", "TNF", "GNLY"),
        "Treg_Markers" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3"),
        "Exhaustion" = c("PDCD1", "HAVCR2", "TIGIT"),
        "Tissue_Residency" = c("ITGAE", "CD69", "CXCR6")
    )

    # 初始化结果列表
    correlation_results <- list()

    # Only analyze tumor group
    sc_group <- sc_tex %>%
        filter(group == "tumor")

    # 提取感兴趣基因的表达矩阵
    genes_to_analyze <- unique(c("ENTPD1", unlist(gene_sets)))
    expr_matrix <- GetAssayData(sc_group, slot = "counts", assay = "RNA") %>%
        as.matrix() %>%
        t() %>%
        sclink_norm(scale.factor = 1e6, filter.genes = FALSE, gene.names = genes_to_analyze)

    # 使用scLink计算相关性
    corr <- sclink_net(expr = expr_matrix, ncores = 30)
    correlation_results <- corr$cor

    # Process correlation matrix
    process_correlation_matrix <- function(corr_matrix) {
        # Only keep correlations with ENTPD1
        entpd1_corr <- corr_matrix["ENTPD1", ]

        # Convert to data frame
        df <- data.frame(
            Gene = names(entpd1_corr),
            Correlation = entpd1_corr
        ) %>%
            # Remove self-correlation
            filter(Gene != "ENTPD1")

        # Add gene set information
        df <- df %>%
            mutate(Gene_Set = case_when(
                Gene %in% gene_sets$Adenosine_Receptors ~ "Adenosine_Receptors",
                Gene %in% gene_sets$Cytotoxicity ~ "Cytotoxicity",
                Gene %in% gene_sets$Treg_Markers ~ "Treg_Markers",
                Gene %in% gene_sets$Exhaustion ~ "Exhaustion",
                Gene %in% gene_sets$Tissue_Residency ~ "Tissue_Residency"
            ))

        return(df)
    }

    all_correlations <- process_correlation_matrix(correlation_results)

    # Save detailed results
    write_csv(
        all_correlations,
        "results/108.Tcell/correlation/CD8_Tex_TRM_ENTPD1_correlations_tumor.csv"
    )

    # Create visualization
    p2 <- ggplot(
        all_correlations,
        aes(
            x = reorder(Gene, Correlation),
            y = Correlation
        )
    ) +
        geom_bar(stat = "identity", aes(fill = Correlation > 0)) +
        scale_fill_manual(values = c("#4169E1", "#FF6B6B")) +
        geom_text(aes(label = sprintf("%.2f", Correlation)),
            vjust = ifelse(all_correlations$Correlation > 0, -0.5, 1.5)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        ) +
        labs(
            x = "Gene",
            y = "scLink Correlation with ENTPD1",
            title = "Correlation of ENTPD1 with Selected Genes in CD8_Tex_TRM (Tumor)"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_grid(~Gene_Set, scales = "free_x", space = "free")
    ggsave("results/108.Tcell/correlation/CD8_Tex_TRM_ENTPD1_correlations_tumor.png",
        p2,
        width = 15, height = 7
    )
}
run_GSEA_Tex_TRM <- function(sc_tcell) {
    options(max.print = 12, spe = "human")
    dir.create("results/108.Tcell/GSEA_Tex_TRM", showWarnings = FALSE, recursive = TRUE)

    # Create custom gene set of adenosine
    adenosine_genes <- SearchDatabase(
        item = "adenosine",
        type = "SetName",
        return = "genelist",
        database = "GO",
        n.min = 3
    )

    # Run GSEA for adenosine gene set
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = adenosine_genes, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for adenosine pathways
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20),
        lab_fill = "zscore"
    ) +
        theme(text = element_text(size = 14))
    ggsave("results/108.Tcell/GSEA_Tex_TRM/Heatmap_adenosine.png", p, width = 14, height = 14)
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA_Tex_TRM/Heatmap_adenosine.tsv"
    )
    # Create custom gene set of hypoxia
    hypoxia_genes <- SearchDatabase(
        item = "hypoxia",
        type = "SetName",
        return = "genelist",
        database = "GO",
        n.min = 3
    )

    # Run GSEA for hypoxia gene set
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = hypoxia_genes, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for hypoxia pathways
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20),
        lab_fill = "zscore"
    ) +
        theme(text = element_text(size = 14))
    ggsave("results/108.Tcell/GSEA_Tex_TRM/Heatmap_hypoxia.png", p, width = 14, height = 14)
    # write tsv files
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA_Tex_TRM/Heatmap_hypoxia.tsv"
    )
    "results/108.Tcell/GSEA_Tex_TRM"
}

run_tcell_trajectory <- function(sc_tcell) {
    dir.create("results/108.Tcell/trajectory", showWarnings = FALSE, recursive = TRUE)
    py_run_string("")

    # 1. scVelo Analysis
    # Convert Seurat object to AnnData for scVelo
    dir.create("data/108.Tcell/trajectory", showWarnings = FALSE, recursive = TRUE)
    Seu2Adata(sc_tcell, save.adata = "data/108.Tcell/trajectory/tcell.h5ad", conda_env = "base")
    adata_path <- "data/108.Tcell/trajectory/tcell.h5ad"
    seu <- sc_tcell
    # seu <- scVelo.SeuratToAnndata(
    #     sc_tcell,
    #     filename = adata_path,
    #     velocyto.loompath = "data/103.self_workflow/velocyto_combined.loom",
    #     prefix = "",
    #     postfix = "",
    #     remove_duplicates = TRUE,
    #     conda_env = "base"
    # )

    # # Generate velocity plots
    # scVelo.Plot(
    #     color = "cell_type_dtl",
    #     basis = "umap_tcell_cell_embeddings",
    #     save = "results/108.Tcell/trajectory/velocity_umap.png",
    #     figsize = c(7, 6),
    #     conda_env = "base"
    # )

    # 2. Palantir Analysis
    # Run diffusion map
    set.seed(42)
    seu <- Palantir.RunDM(seu,
        reduction = "harmony",
        conda_env = "base"
    )
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    cells_to_remove_1 <- CellSelector(p)
    # Remove selected cells
    seu <- subset(seu, cells = setdiff(colnames(seu), cells_to_remove_1))
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    # Create plot again after removing cells
    cells_to_remove_2 <- CellSelector(p)
    seu <- subset(seu, cells = setdiff(colnames(seu), cells_to_remove_2))
    p <- DimPlot2(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/trajectory/plantir_dm.png", p, width = 7, height = 6)
    # Get cells from seu and update adata to only include those cells
    cells <- colnames(seu)
    py_run_string(sprintf("adata = adata[adata.obs.index.isin(%s)]", paste0("[", paste0("'", cells, "'", collapse = ","), "]")))
    # Add ms dimension reduction to AnnData object
    adata.AddDR(seu, dr = "ms", scv.graph = FALSE, conda_env = "base")

    # Plot velocity using ms coordinates
    # scVelo.Plot(
    #     color = "cell_type_dtl",
    #     basis = "ms",
    #     save = "results/108.Tcell/trajectory/velocity_ms.png",
    #     figsize = c(7, 6),
    #     conda_env = "base"
    # )

    # Select start cell (CD4_TN cluster)
    p <- DimPlot(seu %>% subset(cell_type_dtl == "CD4_TN"), reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    start_cell <- CellSelector(p)
    start_cell <- colnames(seu)[which(seu$cell_type_dtl == "CD4_TN")[1]]

    # Calculate pseudotime
    seu <- Palantir.Pseudotime(seu, start_cell = start_cell, conda_env = "base", n_jobs = 10)

    # Get pseudotime data
    ps <- seu@misc$Palantir$Pseudotime
    colnames(ps)[3:4] <- c("fate1", "fate2")
    seu@meta.data[, colnames(ps)] <- ps

    # Plot pseudotime and entropy
    p1 <- DimPlot2(seu,
        features = colnames(ps),
        reduction = "ms",
        cols = list(Entropy = "D"),
        theme = NoAxes()
    ) + theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/trajectory/palantir_pseudotime.png", p1, width = 12, height = 12)

    # # 3. MAGIC for gene expression smoothing
    # sc_tcell <- Palantir.Magic(sc_tcell)
    # sc_tcell <- NormalizeData(sc_tcell)

    # # Plot gene expression before and after MAGIC
    # key_genes <- c("CD4", "CD8A", "FOXP3", "PDCD1")
    # p2 <- DimPlot2(sc_tcell,
    #     features = c(key_genes, paste0("magic_", key_genes)),
    #     cols = "A",
    #     theme = NoAxes()
    # )
    # ggsave("results/108.Tcell/trajectory/magic_expression.png", p2, width = 12, height = 12)

    # 4. CellRank Analysis
    # Add pseudotime to adata
    # Add pseudotime to adata
    adata.AddMetadata(seu, col = colnames(ps), conda_env = "base")

    # Run CellRank
    Cellrank.Compute(time_key = "Pseudotime", conda_env = "base")

    # Generate CellRank plots
    Cellrank.Plot(
        color = "cell_type_dtl",
        basis = "ms",
        save = "results/108.Tcell/trajectory/cellrank_ms.png",
        conda_env = "base"
    )

    # 5. Gene Expression Dynamics
    # Generate trend curves for key genes
    key_genes <- c("CD4", "CD8A", "FOXP3", "PDCD1", "CTLA4", "LAG3")
    p3 <- GeneTrendCurve.Palantir(
        seu,
        pseudotime.data = ps,
        features = key_genes,
        point = FALSE,
        se = TRUE,
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/gene_trends.png", p3, width = 10, height = 8)

    fate1_genes <- c("ENTPD1", "HAVCR2", "ITGAE", "PDCD1", "LAG3", "TIGIT", "CXCR6", "CD69", "PRDM1", "IL2RA", "FASLG", "NCR3") %>%
        unique()
    fate2_genes <- c(
        "FOXP3",
        "IFNG", # Key cytotoxic cytokine
        "GZMB", # Granzyme B effector molecule
        "PRF1", # Perforin
        "GNLY", # Granulysin
        "NKG7", # Cytolytic granule component
        "TRDC", # γδ-TCR chain
        "TRGC1", # γδ-TCR constant region
        "TRDV1", # αβ-TCR chain
        "TRBC1", # αβ-TCR constant region
        "KLRD1", # CD94 (NK cell activation)
        "FCGR3A", # CD16 (ADCC potential)
        "TNF", # Pro-inflammatory cytokine
        "TBX21", # T-bet (CTL master regulator)
        "CX3CR1", # Tissue surveillance receptor
        "GZMH", # Early activation granzyme
        "CCL5", # Chemokine for immune recruitment
        "STAT4", # IL-12 signaling
        "CD8A" # Core CTL marker
    )
    # Generate heatmap for fate1 genes
    p4 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate1_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/gene_heatmap_fate1.png", p4, width = 12, height = 8)


    p5 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate2_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/gene_heatmap_fate2.png", p5, width = 12, height = 8)

    # 6. Slingshot Analysis
    sc_tcell <- RunSlingshot(sc_tcell,
        group.by = "cell_type_dtl",
        start.clus = "CD4_TN"
    )

    # Add Slingshot pseudotime to metadata
    sling <- sc_tcell@misc$slingshot$PCA$SlingPseudotime
    sc_tcell@meta.data[, colnames(sling)] <- as.data.frame(sling)

    # Plot Slingshot results
    p5 <- DimPlot2(sc_tcell,
        features = colnames(sling),
        cols = "C",
        theme = NoAxes()
    )
    ggsave("results/108.Tcell/trajectory/slingshot_pseudotime.png", p5, width = 12, height = 5)

    return(sc_tcell)
}
