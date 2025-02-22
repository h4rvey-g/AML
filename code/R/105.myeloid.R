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
            filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
            arrange(desc(abs(avg_log2FC)))
        if (nrow(deg) == 0) {
            deg <- FindMarkers(
                sc_int,
                ident.1 = cluster,
                test.use = "MAST",
                min.cells.group = 2
            ) %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                filter(p_val_adj < 0.05) %>%
                filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                arrange(desc(abs(avg_log2FC)))
        }

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
                filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
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
        "1" = "M2_Mac",
        "2" = "FOLR2_Mac",
        "3" = "Mast",
        "4" = "TREM2_Mac",
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
    markers <- list(
        "TREM2_Mac" = c("CD40", "STAT1", "TREM2"), # TREM2+ Macrophages
        "FOLR2_Mac" = c("FOLR2", "BAG5", "TRIM47"), # FOLR2+ Macrophages
        "M2_Mac" = c("CCL18", "CD163", "MRC1"), # M2 Macrophages
        "Mast" = c("TPSAB1", "LYZ", "CCR2"), # Mast Cells
        "Eosino" = c("IL5RA", "EPO", "CXCR2"), # Eosinophils
        "Foam" = c("APOE", "ID1", "IL12B") # Foam Cells
    )

    # 计算并绘制热图
    toplot <- CalcStats(sc_mye_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/myeloid_annotate_heatmap.png", p2, width = 10, height = 6)

    return(sc_mye_clust)
}

test_myeloid_composition <- function(sc_mye) {
    # Create directory
    dir.create("results/107.myeloid/distribution", showWarnings = FALSE, recursive = TRUE)

    # Prepare data for sccomp
    sc_mye <- sc_mye %>%
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

    # Run composition test
    composition_test <- sc_mye %>%
        sccomp_estimate(
            formula_composition = ~group,
            .sample = sample,
            .cell_group = cell_group,
            bimodal_mean_variability_association = TRUE,
            cores = 20,
            verbose = TRUE
        ) %>%
        sccomp_remove_outliers(cores = 31, verbose = FALSE) %>%
        sccomp_test()

    # Generate plots
    p1 <- composition_test %>%
        sccomp_boxplot(factor = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save fold changes and test results
    fold_changes <- composition_test %>%
        sccomp_proportional_fold_change(
            formula_composition = ~group,
            from = "normal",
            to = "tumor"
        ) %>%
        select(cell_group, statement)

    write_tsv(fold_changes, "results/107.myeloid/distribution/fold_changes.tsv")

    write_tsv(
        composition_test %>%
            select(-count_data),
        "results/107.myeloid/distribution/test_results.tsv"
    )

    composition_test
}

plot_myeloid_distribution <- function(sc_mye, m_composition_test) {
    dir.create("results/107.myeloid/distribution", showWarnings = FALSE, recursive = TRUE)

    # Calculate distribution percentages
    df1 <- ClusterDistrBar(origin = sc_mye$dataset, cluster = sc_mye$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")

    # Calculate cell counts
    cell_counts <- sc_mye@meta.data %>%
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
    # Get significance stars from composition test
    m_composition_test_res <- m_composition_test %>%
        filter(factor == "group", c_FDR < 0.05) %>%
        mutate(stars = case_when(
            c_FDR < 0.001 ~ "***",
            c_FDR < 0.01 ~ "**",
            c_FDR < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        select(cell_group, stars)

    # Add significance stars and cell counts to labels
    df1 <- df1 %>%
        mutate(cluster = map_chr(str_extract(cluster, ".*(?=\\n)"), ~ {
            stars <- m_composition_test_res$stars[m_composition_test_res$cell_group == .x]
            count_info <- sprintf("(N=%d,T=%d)", cell_counts$N_count[cell_counts$cell_type_dtl == .x], cell_counts$T_count[cell_counts$cell_type_dtl == .x])
            if (length(stars) > 0 && stars != "") {
                paste0(.x, "\n", count_info, "\n", stars)
            } else {
                paste0(.x, "\n", count_info, "\n")
            }
        }))

    # Create and save plot
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

    tidyplots::save_plot(p, "results/107.myeloid/distribution/myeloid_distribution.png")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/107.myeloid/distribution/myeloid_distribution.tsv"
    )

    return("results/107.myeloid/distribution")
}

run_myeloid_GSEA <- function(sc_mye) {
    options(max.print = 12, spe = "human")
    dir.create("results/107.myeloid/GSEA", showWarnings = FALSE, recursive = TRUE)

    # Run GO analysis
    parents <- GeneSetAnalysisGO() %>% names()
    for (parent in parents) {
        dir.create(paste0("results/107.myeloid/GSEA/", parent), showWarnings = FALSE, recursive = TRUE)
        cat("Running GSEA for", parent, "\n")
        sc_mye <- GeneSetAnalysisGO(sc_mye, nCores = 50, parent = parent)
        matr <- sc_mye@misc$AUCell$GO[[parent]]
        matr <- RenameGO(matr, add_id = FALSE)

        # Overall heatmap for all subclusters
        toplot <- CalcStats(matr, f = sc_mye$cell_type_dtl, method = "zscore", order = "p", n = 10)
        p <- Heatmap(
            toplot,
            lab_fill = "zscore", text.size = 5
        )
        ggsave(paste0("results/107.myeloid/GSEA/", parent, "/Heatmap_GO_", parent, ".png"), p, width = 14, height = nrow(toplot) * 0.4)
        write_tsv(
            toplot %>% as.data.frame() %>% rownames_to_column("pathway"),
            paste0("results/107.myeloid/GSEA/", parent, "/Heatmap_GO_", parent, ".tsv")
        )

        # Compare tumor vs normal within each subcluster
        subclusters <- unique(sc_mye$cell_type_dtl)
        for (subcluster in subclusters) {
            # Create temporary subset for this subcluster
            sc_temp <- sc_mye %>%
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
                    "results/107.myeloid/GSEA/", parent, "/Waterfall_GO_",
                    parent, "_", gsub("/", "_", subcluster), "_tumor_vs_normal.png"
                ), p, width = 15, height = 15)
            }
        }
    }

    # Run hallmark analysis
    sc_mye <- GeneSetAnalysis(sc_mye, genesets = hall50$human, nCores = 30)
    matr <- sc_mye@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for hallmark pathways
    dir.create("results/107.myeloid/GSEA/hallmark50", showWarnings = FALSE, recursive = TRUE)
    p <- Heatmap(
        CalcStats(matr, f = sc_mye$cell_type_dtl, method = "zscore", order = "p", n = 10),
        lab_fill = "zscore", text.size = 5
    )
    ggsave("results/107.myeloid/GSEA/hallmark50/Heatmap_hallmark50.png", p, width = 14, height = 14)
    write_tsv(
        CalcStats(matr, f = sc_mye$cell_type_dtl, method = "zscore", order = "p", n = 10) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/107.myeloid/GSEA/hallmark50/Heatmap_hallmark50.tsv"
    )

    # Compare tumor vs normal within each subcluster for hallmark pathways
    subclusters <- unique(sc_mye$cell_type_dtl)
    for (subcluster in subclusters) {
        # Create temporary subset for this subcluster
        sc_temp <- sc_mye %>%
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
                title = paste0("Hallmark50: ", subcluster)
            )
            ggsave(paste0(
                "results/107.myeloid/GSEA/hallmark50/Waterfall_hallmark50_",
                gsub("/", "_", subcluster), "_tumor_vs_normal.png"
            ), p, width = 15, height = 15)
        }
    }

    # Run for each GO root category
    for (root in c("BP", "MF", "CC")) {
        cat("Processing", root, "...\n")

        sc_mye <- GeneSetAnalysisGO(sc_mye, nCores = 50, root = root)
        matr <- sc_mye@misc$AUCell$GO[[root]]

        top30_pathways <- data.frame()
        for (cell_type in unique(sc_mye$cell_type_dtl)) {
            stats <- CalcStats(matr, f = sc_mye$cell_type_dtl, method = "zscore") %>%
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
        toplot <- CalcStats(matr_filt, f = sc_mye$cell_type_dtl, method = "zscore", order = "p", n = 30)

        # Create directory for each root category
        dir.create(file.path("results/107.myeloid/GSEA/GO", root), showWarnings = FALSE, recursive = TRUE)

        # Save heatmap and data
        p <- Heatmap(
            CalcStats(matr_filt, f = sc_mye$cell_type_dtl, method = "zscore", order = "p", n = 10),
            lab_fill = "zscore", text.size = 5
        )
        ggsave(file.path("results/107.myeloid/GSEA/GO", root, "Heatmap_top10.png"), p, width = 14, height = 14)
        write_tsv(
            toplot %>%
                as.data.frame() %>% rownames_to_column("pathway"),
            file.path("results/107.myeloid/GSEA/GO", root, "Heatmap_top10.tsv")
        )
    }


    return("results/107.myeloid/GSEA")
}
run_myeloid_customGSEA <- function(sc_mye) {
    dir.create("results/107.myeloid/customGSEA", showWarnings = FALSE, recursive = TRUE)

    # Define pathways of interest for each cell type
    pathways <- list(
        "Foam" = c(
            "Fatty Acid Beta-Oxidation",
            "Brown Fat Cell Differentiation",
            "Contractile Actin Filament Bundle",
            "Stress Fiber"
        ),
        "TREM2_Mac" = c(
            "Complement-Mediated Synapse Pruning",
            "Clathrin-Coated Endocytic Vesicle Membrane",
            "Membrane Hyperpolarization",
            "Smooth Muscle Cell-Matrix Adhesion",
            "Mono-ADP-D-Ribose Binding"
        ),
        "FOLR2_Mac" = c(
            "Positive Regulation of Cardiac Muscle Relaxation",
            "Proteasome Regulatory Particle",
            "Response to Carbon Dioxide",
            "Cellular Response to Staurosporine",
            "V1b Vasopressin Receptor Binding"
        ),
        "M2_Mac" = c(
            "Canonical Wnt Signaling Pathway Involved In Positive Regulation Of Epithelial To Mesenchymal Transition",
            "Regulation of Response to Oxidative Stress",
            "Lipid Droplet",
            "Monocyte Chemotaxis",
            "Phagolysosome Assembly"
        ),
        "Mast" = c(
            "Cell Surface Pattern Recognition Receptor Signaling Pathway",
            "Myeloid Progenitor Cell Differentiation",
            "Fc-Gamma Receptor I Complex Binding",
            "Regulation of Fibroblast Growth Factor Receptor Signaling Pathway",
            "Interleukin-18 Receptor Complex"
        ),
        "Eosino" = c(
            "Ion Channel Complex",
            "Transmitter-Gated Ion Channel Activity Involved In Regulation Of Postsynaptic Membrane Potential",
            "Presynapse Assembly",
            "Gephyrin Clustering Involved In Postsynaptic Density Assembly"
        )
    )

    options(max.print = 12, spe = "human")
    # Get custom gene sets for each pathway
    custom_sets <- lapply(unlist(pathways), function(x) {
        SearchDatabase(item = x, type = "SetName", return = "genelist", database = "GO")
    })
    # Filter each sublist to keep only the longest gene set
    custom_sets <- lapply(custom_sets, function(x) {
        # Find the longest gene set
        lengths <- sapply(x, length)
        x[which.max(lengths)]
    })
    names(custom_sets) <- rep(names(pathways), lengths(pathways))
    # Get lengths and sum up pathway lengths for each cell type using pipes
    # Calculate pathway lengths using base R operations
    pathway_lengths <- tapply(sapply(custom_sets, length), 
                            factor(names(custom_sets), 
                                  levels = unique(names(custom_sets))), 
                            sum)
    # Create repeated names based on lengths
    pathway_lengths <- rep(names(pathway_lengths), pathway_lengths)

    # Combine all gene sets into a single list
    custom_sets <- Reduce(c, custom_sets)
    names(pathway_lengths) <- names(custom_sets) %>%
        RenameGO(add_id = FALSE)

    # Run AUCell analysis with custom gene sets
    sc_mye <- GeneSetAnalysis(sc_mye, genesets = custom_sets, nCores = 30, n.min = 3)
    matr <- sc_mye@misc$AUCell$genesets
    toplot <- CalcStats(matr %>% RenameGO(add_id = FALSE),
        f = sc_mye$cell_type_dtl,
        method = "zscore",
        order = "value"
    )

    p2 <- Heatmap(toplot,
        lab_fill = "zscore",
        facet_row = pathway_lengths[rownames(toplot)]
    ) +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/107.myeloid/customGSEA/grouped_heatmap.png", p2, width = 14, height = 8)

    # Compare tumor vs normal for each cell type
    for (cell_type in names(pathways)) {
        sc_subset <- subset(sc_mye, cell_type_dtl == cell_type)
        if (length(unique(sc_subset$group)) == 2) {
            relevant_pathways <- pathways[[cell_type]]
            matr_subset <- matr[relevant_pathways, ]

            p2 <- WaterfallPlot(
                matr_subset,
                f = sc_subset$group,
                ident.1 = "tumor",
                ident.2 = "normal",
                color = "p",
                length = "logFC",
                title = paste0("Custom Pathways: ", cell_type)
            )
            ggsave(paste0(
                "results/107.myeloid/customGSEA/",
                gsub("/", "_", cell_type),
                "_tumor_vs_normal.png"
            ), p2, width = 10, height = 8)
        }
    }

    return("results/107.myeloid/customGSEA")
}
myeloid_DEG_tumor_vs_normal <- function(sc_mye) {
    sc_pseudo <- AggregateExpression(sc_mye, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    cell_types <- unique(sc_pseudo$cell_type_dtl)
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

        if (is.null(deg)) next
        dir.create("results/107.myeloid/DEG_tumor_vs_normal", showWarnings = FALSE, recursive = TRUE)
        write_tsv(deg, paste0("results/107.myeloid/DEG_tumor_vs_normal/", cell_type, "_DEG.tsv"))
    }
    sc_mye$cell_type_group <- paste(sc_mye$cell_type_dtl, sc_mye$group, sep = "_")
    Idents(sc_mye) <- "cell_type_group"
    cell_types <- unique(sc_mye$cell_type_dtl)
    # Set up parallel processing
    cores <- min(parallel::detectCores() - 1, length(cell_types))
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    # Parallel foreach loop
    results <- foreach(cell_type = cell_types, .packages = c("Seurat", "tidyverse", "dplyr")) %dopar% {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if ((ident1 %in% Idents(sc_mye) && ident2 %in% Idents(sc_mye))) {
            deg <- FindMarkers(
                sc_mye,
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
            write_tsv(deg, paste0("results/107.myeloid/DEG_tumor_vs_normal/", cell_type, "_DEG_MAST.tsv"))
        }
    }

    # Stop cluster
    parallel::stopCluster(cl)
    return("results/107.myeloid/DEG_tumor_vs_normal")
}

plot_myeloid <- function(sc_mye) {
    # Create directory for plots
    dir.create("results/107.myeloid/plots", showWarnings = FALSE, recursive = TRUE)

    # Define B cell related genes to plot by category
    b_cell_genes_by_category <- list(
        "Steroidogenesis" = c("STAR", "CYP11B1", "HSD3B2"),
        "Lipid Metabolism" = c("SCD5", "CLU", "PEBP1"),
        "Immunoglobulins" = c("IGHG1", "IGHG4", "IGLC2")
    )

    # Subset TREM2+ macrophages and calculate expression statistics
    trem2_cells <- sc_mye[, sc_mye$cell_type_dtl == "TREM2_Mac"]

    # Extract expression matrix for B cell genes
    all_genes <- unlist(b_cell_genes_by_category)
    expr_mat <- GetAssayData(trem2_cells, slot = "data")[all_genes, ]

    # Scale the expression matrix
    expr_mat_scaled <- t(scale(t(expr_mat)))

    # Get tumor/normal status for each cell
    group_order <- trem2_cells$group %>% factor(levels = c("normal", "tumor"))
    column_order <- order(group_order)

    # Create annotation for tumor/normal
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        Group = trem2_cells$group[column_order],
        col = list(Group = c("normal" = "#00BFC4", "tumor" = "#F8766D"))
    )

    # Create gene category split
    gene_categories <- rep(names(b_cell_genes_by_category), sapply(b_cell_genes_by_category, length))
    names(gene_categories) <- unlist(b_cell_genes_by_category)

    # Create heatmap with split columns by group and split rows by category
    p <- ComplexHeatmap::Heatmap(as.matrix(expr_mat[, column_order]),
        name = "Expression",
        col = viridis::viridis(100),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_title = "TREM2+ Macrophages",
        row_title = "Genes related to immunoglobulins, lipid metabolism",
        top_annotation = anno_col,
        column_split = trem2_cells$group[column_order],
        row_split = gene_categories[rownames(expr_mat)],
        row_gap = unit(5, "mm")
    )

    # Save the heatmap
    png("results/107.myeloid/plots/trem2_mac_heatmap.png", width = 1000, height = 800, res = 100)
    ComplexHeatmap::draw(p)
    dev.off()

    # 定义基因组
    mast_genes <- list(
        "Steroidogenesis" = c("STAR", "CYP11B1", "HSD3B2"),
        "Immune" = c("IL7R", "ANGPT1", "C3", "RARRES2"),
        "Detoxification" = c("GSTA1", "GSTA4", "AOX1", "SULT2A1"),
        "ECM_Adhesion" = c("NOV", "SPON2", "DCN", "SLIT2", "SORBS2", "CDH2", "ITGA1", "LAMA2")
    )

    # 提取Mast细胞
    mast_cells <- sc_mye[, sc_mye$cell_type_dtl == "Mast"]

    # 提取基因表达矩阵
    all_genes <- unlist(mast_genes)
    expr_mat <- GetAssayData(mast_cells, slot = "data")[all_genes, ]

    # 标准化表达矩阵
    expr_mat_scaled <- t(scale(t(expr_mat)))

    # 获取tumor/normal分组信息
    group_order <- mast_cells$group %>% factor(levels = c("normal", "tumor"))
    column_order <- order(group_order)

    # 创建分组注释
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        Group = mast_cells$group[column_order],
        col = list(Group = c("tumor" = "#F8766D", "normal" = "#00BFC4"))
    )

    # 创建基因类别分组
    gene_categories <- rep(names(mast_genes), sapply(mast_genes, length))
    names(gene_categories) <- unlist(mast_genes)

    # 创建热图
    p <- ComplexHeatmap::Heatmap(as.matrix(expr_mat[, column_order]),
        name = "Expression",
        col = viridis::viridis(100),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_title = "Mast Cells",
        row_title = "Steroidogenesis and Immune Related Genes",
        top_annotation = anno_col,
        column_split = mast_cells$group[column_order],
        row_split = gene_categories[rownames(expr_mat)],
        row_gap = unit(5, "mm")
    )

    # 保存热图
    png("results/107.myeloid/plots/mast_cells_heatmap.png",
        width = 1000, height = 800, res = 100
    )
    ComplexHeatmap::draw(p)
    dev.off()
    # 定义基因组
    m2_genes <- list(
        "Steroidogenesis" = c("STAR", "CYP11B1", "HSD3B2"),
        "Metabolic" = c("GSTA1", "MEG3", "MEG8", "DCN"),
        "Immune" = c("TIMD4", "ITGAD", "VCAM1")
    )

    # 提取M2巨噬细胞
    m2_cells <- sc_mye[, sc_mye$cell_type_dtl == "M2_Mac"]

    # 获取基因表达矩阵
    all_genes <- unlist(m2_genes)
    expr_mat <- GetAssayData(m2_cells, slot = "data")[all_genes, ]

    # 标准化表达矩阵
    expr_mat_scaled <- t(scale(t(expr_mat)))

    # 获取tumor/normal分组信息
    group_order <- m2_cells$group %>% factor(levels = c("normal", "tumor"))
    column_order <- order(group_order)

    # 创建分组注释
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        Group = m2_cells$group[column_order],
        col = list(Group = c("tumor" = "#F8766D", "normal" = "#00BFC4"))
    )

    # 创建基因类别分组
    gene_categories <- rep(names(m2_genes), sapply(m2_genes, length))
    names(gene_categories) <- unlist(m2_genes)

    # 创建热图
    p1 <- ComplexHeatmap::Heatmap(as.matrix(expr_mat_scaled[, column_order]),
        name = "Expression",
        col = viridis::viridis(100),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_title = "M2 Macrophages",
        row_title = "Key M2 Macrophage Genes",
        top_annotation = anno_col,
        column_split = m2_cells$group[column_order],
        row_split = gene_categories[rownames(expr_mat)],
        row_gap = unit(5, "mm")
    )

    # 保存热图
    png("results/107.myeloid/plots/m2_macrophage_heatmap.png",
        width = 1000, height = 800, res = 100
    )
    ComplexHeatmap::draw(p1)
    dev.off()
}
