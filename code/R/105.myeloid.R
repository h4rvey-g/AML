sub_cluster_myeloid <- function(sc_int) {
    # 提取巨噬细胞和粒细胞
    sc_int <- sc_int %>%
        filter(cell_type_dtl == "Myeloid")

    # 重新聚类
    sc_int <- sc_int %>%
        FindNeighbors(dims = 1:30, reduction = "harmony") %>%
        FindClusters(resolution = 0.4, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster")

    sc_int
}
sub_annotation_myeloid <- function(sc_int, sc_final) {
    # # remove cluster 7 cells first
    # sc_int <- sc_int %>% filter(immune_subcluster != "7")
    # 可视化聚类结果
    p1 <- DimPlot2(sc_int,
        reduction = "umap_integrated",
        group.by = "immune_subcluster", label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    dir.create("results/107.myeloid", showWarnings = FALSE, recursive = TRUE)
    ggsave("results/107.myeloid/sub_cluster_immune.png", p1, width = 7, height = 7)

    markers <- list(
        # Tissue-resident & specialized macrophages
        "LYVE1_Mac" = c("LYVE1", "CD163", "MRC1", "FOLR2"), # Removed MRC1L1 (non-human ortholog)
        "Kupffer_like" = c("CLEC4F", "ID1"),
        "Microglia_like" = c("P2RY12", "TMEM119"),
        "Alveolar_Mac" = c("MARCO", "PPARG"),

        # Monocytes
        "Classical_Mono" = c("CD14", "VCAN", "CCR2", "LYZ", "S100A8", "S100A9"),
        "NonClassical_Mono" = c("FCGR3A", "CX3CR1", "LST1", "MS4A7", "NR4A1"),

        # Inflammatory/activated macrophages
        "Mac_M1" = c("CD80", "CD86", "STAT1", "IL1B", "CXCL10", "HLA-DRA", "NOS2", "IL6", "IL12B"),
        "Mac_M2" = c("CD163", "MRC1", "FOLR2", "CCL18", "MAFB"),

        # Lipid-associated macrophages / foam cells
        "Mac_LAM" = c("TREM2", "APOE", "LPL", "FABP4"),
        "Foam_Cell" = c("CD36", "MSR1", "PLIN2", "ADRP"),

        # Erythroid lineage (exclude these if only immune focus)
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2"),

        # Granulocytes
        "Neutrophil" = c("CSF3R", "CEACAM8", "CXCR2", "LCN2", "MPO"),
        "Eosinophil" = c("CLC", "PRG2", "IL5RA", "EPX", "ALOX15"),
        "Mast" = c("KIT", "TPSAB1", "CMA1", "CPA3"),

        # Dendritic cells
        "cDC1" = c("CLEC9A", "XCR1", "BATF3", "IRF8"),
        "cDC2" = c("CD1C", "FCER1A", "IRF4", "ITGAX"),
        "mo_DC" = c("CD1C", "FCGR3A", "CCR7", "CD80", "CD86", "LAMP3"),

        # Activation / proliferation
        "Activation" = c("CD83", "CD40", "LAMP1"),
        "Proliferation" = c("MKI67", "TOP2A", "CENPF")
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
    p <- DotPlot2(sc_int,
        features = markers,
        group.by = "immune_subcluster",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/sub_cluster_immune_dotplot.png", p, width = 10, height = 20)
    # Get cells from immune_subcluster 1 and 7
    cells_1 <- rownames(sc_int@meta.data[sc_int$immune_subcluster == "1", ])

    # Update cell_type_dtl in sc_final for these cells
    sc_final$cell_type_dtl[cells_1] <- "unknown_1"
    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2", "HSD3B2", "CYP17A1", "CYP21A2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1"),
        # "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
        "CLC" = c("TH", "CHGA", "CHGB", "KIT", "SYT1"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
        "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "THY1", "ENG", "MCAM"),
        "MSC" = c("NT5E", "CD105", "NES", "CD106", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3"),
        # "PSC_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
        # "PSC_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
        # "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Myeloid" = c("CD33", "CSF1R", "CD68", "S100A8", "CSF3R", "C1QA", "MPO"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        # "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
        "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21", "NR2F2", "ID1"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
    )
    toplot <- CalcStats(sc_final,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
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
    ggsave("results/107.myeloid/unknown_heatmap.png", p1, width = 8, height = 15)
    write_tsv(
        toplot %>% as.data.frame() %>% rownames_to_column("gene"),
        "results/107.myeloid/unknown_expression.tsv"
    )
    # 创建DEG结果目录
    dir.create("results/107.myeloid/immune_DEG", showWarnings = FALSE, recursive = TRUE)
    sc_int
}

myeloid_annotate <- function(sc_mye_clust) {
    # Create mapping of cluster numbers to cell type annotations
    cluster_map <- c(
        "1" = "Tissue-resident_Mac",
        "2" = "TREM2_LAM",
        "3" = "DC",
        "4" = "Mast",
        "5" = "Eosinophil",
        "6" = "APOE_LAM",
        "7" = "Eryth"
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
        "TREM2_LAM" = c("CD40", "STAT1", "TREM2"), # TREM2+ Macrophages
        "DC" = c("PPARG", "CD83", "CD86"), # Foam Macrophages
        "APOE_LAM" = c("APOE", "ID1"),
        "Mast" = c("TPSAB1", "LYZ", "CCR2"), # Mast Cells
        "Eosinophil" = c("IL5RA", "EPO", "CXCR2"), # Eosinophils
        "Tissue-resident" = c("LYVE1", "MRC1", "CD163"), # Perivascular Macrophages
        "Erythrocyte" = c("GYPA", "HBB", "SLC4A1") # Erythrocytes
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

    markers <- list(
        # Tissue-resident & specialized macrophages
        "LYVE1_Mac" = c("LYVE1", "CD163", "MRC1", "FOLR2"), # Removed MRC1L1 (non-human ortholog)
        "Kupffer_like" = c("CLEC4F", "ID1"),
        "Microglia_like" = c("P2RY12", "TMEM119"),
        "Alveolar_Mac" = c("MARCO", "PPARG"),

        # Monocytes
        "Classical_Mono" = c("CD14", "VCAN", "CCR2", "LYZ", "S100A8", "S100A9"),
        "NonClassical_Mono" = c("FCGR3A", "CX3CR1", "LST1", "MS4A7", "NR4A1"),

        # Inflammatory/activated macrophages
        "Mac_M1" = c("CD80", "CD86", "STAT1", "IL1B", "CXCL10", "HLA-DRA", "NOS2", "IL6", "IL12B"),
        "Mac_M2" = c("CD163", "MRC1", "FOLR2", "CCL18", "MAFB"),

        # Lipid-associated macrophages / foam cells
        "Mac_LAM" = c("TREM2", "APOE", "LPL", "FABP4"),
        "Foam_Cell" = c("CD36", "MSR1", "PLIN2", "ADRP"),

        # Erythroid lineage (exclude these if only immune focus)
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2"),

        # Granulocytes
        "Neutrophil" = c("CSF3R", "CEACAM8", "CXCR2", "LCN2", "MPO"),
        "Eosinophil" = c("CLC", "PRG2", "IL5RA", "EPX", "ALOX15"),
        "Mast" = c("KIT", "TPSAB1", "CMA1", "CPA3"),

        # Dendritic cells
        "cDC1" = c("CLEC9A", "XCR1", "BATF3", "IRF8"),
        "cDC2" = c("CD1C", "FCER1A", "IRF4", "ITGAX"),
        "mo_DC" = c("CD1C", "FCGR3A", "CCR7", "CD80", "CD86", "LAMP3"),

        # Activation / proliferation
        "Activation" = c("CD83", "CD40", "LAMP1"),
        "Proliferation" = c("MKI67", "TOP2A", "CENPF")
    )
    toplot <- CalcStats(sc_mye_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    write_tsv(
        toplot %>% as.data.frame() %>% rownames_to_column("gene"),
        "results/107.myeloid/myeloid_annotate_zscore.tsv"
    )
    p <- DotPlot2(sc_mye_clust,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE,
        flip = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/myeloid_annotate_dotplot.png", p, width = 8, height = 10)

    return(sc_mye_clust)
}

myeloid_DEG_vs_cluster <- function(sc_mye, myeloid_cluster) {
    # 按照数据集和亚群聚合进行DEG分析
    sc_pseudo <- AggregateExpression(sc_mye,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c("dataset", "cell_type_dtl")
    )

    # 设置比较分组
    Idents(sc_pseudo) <- "cell_type_dtl"
    dir.create("results/107.myeloid/myeloid_DEG", showWarnings = FALSE, recursive = TRUE)
    # if myeloid_cluster contains _, replace with -
    myeloid_cluster <- gsub("_", "-", myeloid_cluster)

    deg <- FindMarkers(
        sc_pseudo,
        ident.1 = myeloid_cluster,
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
        arrange(desc(abs(avg_log2FC)))
    if (nrow(deg) == 0) {
        myeloid_cluster <- sub("-([^-]*)$", "_\\1", myeloid_cluster)
        Idents(sc_mye) <- "cell_type_dtl"
        deg <- FindMarkers(
            sc_mye,
            ident.1 = myeloid_cluster,
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
        "results/107.myeloid/myeloid_DEG/%s_vs_rest_DEG.tsv",
        myeloid_cluster
    ))
    "results/107.myeloid/myeloid_DEG"
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
        sc_mye <- GeneSetAnalysisGO(sc_mye, nCores = 20, parent = parent)
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

        sc_mye <- GeneSetAnalysisGO(sc_mye, nCores = 20, root = root)
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
myeloid_DEG_tumor_vs_normal <- function(sc_mye, cell_type) {
    sc_pseudo <- AggregateExpression(sc_mye, assays = "RNA", return.seurat = TRUE, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"

    ident1 <- paste0(cell_type, "_tumor")
    ident2 <- paste0(cell_type, "_normal")

    if (!(ident1 %in% Idents(sc_pseudo) && ident2 %in% Idents(sc_pseudo))) {
        return(NULL)
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
                arrange(desc(abs(avg_log2FC))) %>%
                filter(!str_starts(gene, "CYP"))
        },
        error = function(e) {
            message("Error in FindMarkers for ", cell_type, ": ", e$message)
            return(NULL)
        }
    )

    if (is.null(deg)) {
        sc_mye$cell_type_group <- paste(sc_mye$cell_type_dtl, sc_mye$group, sep = "_")
        Idents(sc_mye) <- "cell_type_group"

        if (!(ident1 %in% Idents(sc_mye) && ident2 %in% Idents(sc_mye))) {
            return(NULL)
        }

        deg <- tryCatch(
            {
                FindMarkers(
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
                    arrange(desc(abs(avg_log2FC))) %>%
                    filter(!str_starts(gene, "CYP"))
            },
            error = function(e) {
                message("Error in FindMarkers using MAST for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )
    }

    dir.create("results/107.myeloid/DEG_tumor_vs_normal", showWarnings = FALSE, recursive = TRUE)
    if (!is.null(deg)) {
        write_tsv(deg, paste0("results/107.myeloid/DEG_tumor_vs_normal/", cell_type, "_DEG.tsv"))
    }
    return(paste0("results/107.myeloid/DEG_tumor_vs_normal/", cell_type, "_DEG.tsv"))
}

run_myeloid_trajectory <- function(sc_mye) {
    dir.create("results/107.myeloid/trajectory", showWarnings = FALSE, recursive = TRUE)
    tar_load(sc_mye)
    library(reticulate)
    py_run_string("")
    library(SeuratExtend)
    library(tidyseurat)
    # remove eosinophils and mast
    sc_mye <- sc_mye %>%
        filter(!cell_type_dtl %in% c("Eosinophil", "Mast", "Eryth", "DC"))

    # 1. scVelo Analysis
    # Convert Seurat object to AnnData for scVelo
    dir.create("data/107.myeloid/trajectory", showWarnings = FALSE, recursive = TRUE)
    Seu2Adata(sc_mye, save.adata = "data/107.myeloid/trajectory/myeloid.h5ad", conda_env = "base")
    adata_path <- "data/107.myeloid/trajectory/myeloid.h5ad"
    seu <- scVelo.SeuratToAnndata(
        sc_mye,
        filename = adata_path,
        velocyto.loompath = "data/103.self_workflow/velocyto_combined.loom",
        prefix = "",
        postfix = "",
        remove_duplicates = TRUE,
        conda_env = "base"
    )

    # Generate velocity plots
    scVelo.Plot(
        color = "cell_type_dtl",
        basis = "umap_mye_cell_embeddings",
        save = "results/107.myeloid/trajectory/velocity_umap.png",
        figsize = c(7, 6),
        conda_env = "base"
    )

    # 2. Palantir Analysis
    # Run diffusion map
    seu <- Palantir.RunDM(seu,
        reduction = "harmony",
        conda_env = "base"
    )
    # p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
    #     theme(plot.background = element_rect(fill = "white"))
    # cells_to_remove_1 <- CellSelector(p)
    # Remove selected cells
    # seu <- subset(seu, cells = setdiff(colnames(seu), cells_to_remove_1))

    p <- DimPlot2(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/trajectory/plantir_dm.png", p, width = 7, height = 6)

    # Get cells from seu and update adata
    cells <- colnames(seu)
    py_run_string(sprintf("adata = adata[adata.obs.index.isin(%s)]", paste0("[", paste0("'", cells, "'", collapse = ","), "]")))

    # Add ms dimension reduction to AnnData object
    adata.AddDR(seu, dr = "ms", scv.graph = TRUE, conda_env = "base")

    # Plot velocity using ms coordinates
    scVelo.Plot(
        color = "cell_type_dtl",
        basis = "ms",
        save = "results/107.myeloid/trajectory/velocity_ms.png",
        figsize = c(7, 6),
        conda_env = "base"
    )

    # Select start cell (Eosinophil cluster as starting point)
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    start_cell <- CellSelector(p)
    fate1 <- "N1_CELL71521_N1"
    fate2 <- "N1_CELL78136_N1"
    # fate3 <- "T4_CELL47246_N1"
    # start_cell <- colnames(seu)[which(seu$cell_type_dtl == "Tissue-resident_Mac")[1]]

    # Calculate pseudotime
    seu <- Palantir.Pseudotime(seu,
        start_cell = start_cell, conda_env = "base",
        terminal_states = c("fate1" = fate1, "fate2" = fate2)
    )

    # Get pseudotime data
    ps <- seu@misc$Palantir$Pseudotime
    # colnames(ps)[3:4] <- c("fate1", "fate2")
    seu@meta.data[, colnames(ps)] <- ps

    # Plot pseudotime and entropy
    p1 <- DimPlot2(seu,
        features = colnames(ps),
        reduction = "ms",
        cols = list(Entropy = "D"),
        theme = NoAxes()
    ) + theme(plot.background = element_rect(fill = "white"))
    ggsave("results/107.myeloid/trajectory/palantir_pseudotime.png", p1, width = 12, height = 12)

    # 4. CellRank Analysis
    # Add pseudotime to adata
    adata.AddMetadata(seu, col = colnames(ps), conda_env = "base")

    # Run CellRank
    Cellrank.Compute(time_key = "Pseudotime", conda_env = "base")

    # Generate CellRank plots
    Cellrank.Plot(
        color = "cell_type_dtl",
        basis = "ms",
        save = "results/107.myeloid/trajectory/cellrank_ms.png",
        conda_env = "base"
    )

    # # 5. Gene Expression Dynamics

fate1_genes <- c(
    "TREM2", "TYROBP",     # Lipid & damage sensing
    "SPP1",                # Osteopontin: matrix remodeling
    "CD9", "LGALS3", "LGALS1",   # Galectins & signatures, immunoreg.
    "CLEC7A",              # Damage sensing
    "AXL", "MERTK",        # Apoptotic cell clearance
    "CTSB", "CTSL",        # Lysosomal function
    "GPNMB",               # Fibrosis/immunomod
    "LPL",                 # Lipoprotein uptake
    "FABP4", "FABP5"       # Fatty acid chaperones
)


fate2_genes <- c(
    "APOE", "ABCA1", "ABCG1",
    "PPARG", "NR1H3", "PLTP",
    "MRC1",                # Homeostatic/anti-inflammatory
    "IRAK1", "TRAF6",      # miR-146a NF-κB regulatory axis
    "MERTK",               # Efferocytosis, resolution via ApoE
    "SELENOP"      # Antioxidant defense
)


    # Generate heatmaps for each fate
    p4 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate1_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/107.myeloid/trajectory/gene_heatmap_fate1.png", p4, width = 12, height = 8)

    p5 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate2_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/107.myeloid/trajectory/gene_heatmap_fate2.png", p5, width = 12, height = 8)

    # # 6. Slingshot Analysis
    # sc_mye <- RunSlingshot(sc_mye,
    #     group.by = "cell_type_dtl",
    #     start.clus = "Eosinophil"
    # )

    # # Add Slingshot pseudotime to metadata
    # sling <- sc_mye@misc$slingshot$PCA$SlingPseudotime
    # sc_mye@meta.data[, colnames(sling)] <- as.data.frame(sling)

    # # Plot Slingshot results
    # p6 <- DimPlot2(sc_mye,
    #     features = colnames(sling),
    #     cols = "C",
    #     theme = NoAxes()
    # )
    # ggsave("results/107.myeloid/trajectory/slingshot_pseudotime.png", p6, width = 12, height = 5)

    return(sc_mye)
}