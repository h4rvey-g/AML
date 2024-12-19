DEG_annotation <- function(sc_annotate) {
    sc_pseudo <- AggregateExpression(sc_annotate, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    cell_types <- unique(sc_pseudo$cell_type_dtl)
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_pseudo) && ident2 %in% Idents(sc_pseudo))) {
            next
        }
        deg <- FindMarkers(
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
            arrange(desc(abs(avg_log2FC)))
        write_tsv(deg, paste0("results/103.DEG/DEG_", gsub("/", "_", cell_type), ".tsv"))
    }
    "results/103.DEG"
}

run_GSEA <- function(sc_final) {
    options(max.print = 12, spe = "human")
    parents <- GeneSetAnalysisGO() %>% names()
    # run analysis on every parent category
    for (parent in parents) {
        cat("Running GSEA for", parent, "\n")
        sc_final <- GeneSetAnalysisGO(sc_final, nCores = 50, parent = parent)
        matr <- sc_final@misc$AUCell$GO[[parent]]
        matr <- RenameGO(matr, add_id = FALSE)
        p <- Heatmap(
            CalcStats(matr, f = sc_final$cell_type_dtl, method = "zscore", order = "p", n = 10),
            lab_fill = "zscore"
        )
        ggsave(paste0("results/104.GSEA/key/Heatmap_GO_", parent, ".png"), p, width = 30, height = 30)
        # create a new column in sc_final, the combination of cell_type_dtl and group
        sc_final$cell_type_group <- paste(sc_final$cell_type_dtl, sc_final$group, sep = "_")
        # run WaterfallPlot for each cell type between tumor and normal
        cell_types <- unique(sc_final$cell_type_dtl)
        Idents(sc_final) <- "cell_type_group"
        for (cell_type in cell_types) {
            ident1 <- paste0(cell_type, "_tumor")
            ident2 <- paste0(cell_type, "_normal")
            if (!(ident1 %in% Idents(sc_final) && ident2 %in% Idents(sc_final))) {
                next
            }
            p <- WaterfallPlot(
                matr,
                f = sc_final$cell_type_group,
                ident.1 = ident1,
                ident.2 = ident2,
                top.n = 20
            )
            ggsave(paste0("results/104.GSEA/key/Waterfall_GO_", parent, "_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
        }
    }
    sc_final <- GeneSetAnalysis(sc_final, genesets = hall50$human, nCores = 30)
    matr <- sc_final@misc$AUCell$genesets
    p <- Heatmap(CalcStats(matr, f = sc_final$cell_type_dtl, method = "zscore", order = "p", n = 10), lab_fill = "zscore")
    ggsave("results/104.GSEA/key/Heatmap_hallmark50.png", p, width = 14, height = 14)
    sc_final$cell_type_group <- paste(sc_final$cell_type_dtl, sc_final$group, sep = "_")
    cell_types <- unique(sc_final$cell_type_dtl)
    Idents(sc_final) <- "cell_type_group"
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_final) && ident2 %in% Idents(sc_final))) {
            next
        }
        p <- WaterfallPlot(
            matr,
            f = sc_final$cell_type_group,
            ident.1 = ident1,
            ident.2 = ident2,
            top.n = 20
        )
        ggsave(paste0("results/104.GSEA/key/Waterfall_hallmark50_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
    }
    "results/104.GSEA"
}

run_GSEA_Fib_Cap_Adipo <- function(sc_filtered, parent) {
    # Create directories
    # dirs <- c(
    #     "results/104.GSEA/key/overall",
    #     "results/104.GSEA/key/t_vs_n",
    #     "results/104.GSEA/key/one_to_one"
    # )
    # for (d in dirs) {
    #     dir.create(d, recursive = TRUE, showWarnings = FALSE)
    # }
    options(max.print = 12, spe = "human")
    cat("Running GSEA for", parent, "\n")
    sc_filtered_local <- sc_filtered
    sc_filtered_local <- GeneSetAnalysisGO(sc_filtered_local, nCores = 1, parent = parent)
    matr <- sc_filtered_local@misc$AUCell$GO[[parent]]
    matr <- RenameGO(matr, add_id = FALSE)
    toplot <- CalcStats(matr, f = sc_filtered_local$cell_type_dtl, method = "zscore", order = "p", n = 10)
    p <- Heatmap(
        toplot,
        lab_fill = "zscore"
    )
    ggsave(paste0("results/104.GSEA/key/overall/Heatmap_GO_", parent, ".png"), p, width = 10, height = 10)
    # save to tsv of toplot
    toplot %>%
        as.data.frame() %>%
        rownames_to_column("pathway") %>%
        write_tsv(paste0("results/104.GSEA/key/overall/Heatmap_GO_", parent, ".tsv"))
    sc_filtered_local$cell_type_group <- paste(sc_filtered_local$cell_type_dtl, sc_filtered_local$group, sep = "_")
    cell_types <- unique(sc_filtered_local$cell_type_dtl)
    Idents(sc_filtered_local) <- "cell_type_group"
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_filtered_local) && ident2 %in% Idents(sc_filtered_local))) {
            next
        }
        p <- WaterfallPlot(
            matr,
            f = sc_filtered_local$cell_type_group,
            ident.1 = ident1,
            ident.2 = ident2,
            top.n = 20
        )
        ggsave(paste0("results/104.GSEA/key/t_vs_n/Waterfall_GO_", parent, "_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
    }
    cell_type_combinations <- combn(cell_types, 2, simplify = FALSE)
    Idents(sc_filtered_local) <- "cell_type_dtl"
    for (combo in cell_type_combinations) {
        ident1 <- combo[1]
        ident2 <- combo[2]
        p <- WaterfallPlot(
            matr,
            f = sc_filtered_local$cell_type_dtl,
            ident.1 = ident1,
            ident.2 = ident2,
            top.n = 20
        )
        ggsave(paste0(
            "results/104.GSEA/key/one_to_one/Waterfall_GO_",
            parent, "_", gsub("/", "_", combo[1]), "_vs_", gsub("/", "_", combo[2]), ".png"
        ), p, width = 15, height = 15)
    }
    "results/104.GSEA"
}
compare_cell_types_hallmark50 <- function(sc_filtered) {
    # Run analysis using hallmark50 dataset
    options(max.print = 12, spe = "human")
    sc_filtered <- GeneSetAnalysis(sc_filtered, genesets = hall50$human, nCores = 30)
    matr <- sc_filtered@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    # rownames to lower case, but keep first letter uppercase
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))
    toplot <- CalcStats(matr, f = sc_filtered$cell_type_dtl, method = "zscore", order = "p", n = 10)
    p <- Heatmap(toplot, lab_fill = "zscore")
    ggsave("results/104.GSEA/key/overall/Heatmap_hallmark50.png", p, width = 14, height = 14)
    toplot %>%
        as.data.frame() %>%
        rownames_to_column("pathway") %>%
        write_tsv("results/104.GSEA/key/overall/Heatmap_hallmark50.tsv")

    # Compare tumor vs normal for each cell type
    sc_filtered$cell_type_group <- paste(sc_filtered$cell_type_dtl, sc_filtered$group, sep = "_")
    cell_types <- unique(sc_filtered$cell_type_dtl)
    Idents(sc_filtered) <- "cell_type_group"
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_filtered) && ident2 %in% Idents(sc_filtered))) {
            next
        }
        p <- WaterfallPlot(
            matr,
            f = sc_filtered$cell_type_group,
            ident.1 = ident1,
            ident.2 = ident2,
            top.n = 20
        )
        ggsave(paste0("results/104.GSEA/key/t_vs_n/Waterfall_hallmark50_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
    }
    "results/104.GSEA"
}

plot_distribution <- function(sc_final) {
    p1 <- ClusterDistrBar(origin = sc_final$orig.ident, cluster = sc_final$cell_type_dtl)
    p2 <- ClusterDistrBar(origin = sc_final$group, cluster = sc_final$cell_type_dtl)
    p <- wrap_plots(list(p1, p2), nrow = 2, heights = c(4, 1))
    ggsave("results/105.Distribution/ClusterDistrBar.png", p, width = 10, height = 13)
    p1 <- ClusterDistrBar(origin = sc_final$orig.ident, cluster = sc_final$cell_type_dtl, rev = TRUE, normalize = TRUE)
    p2 <- ClusterDistrBar(origin = sc_final$group, cluster = sc_final$cell_type_dtl, rev = TRUE, normalize = TRUE)
    p <- wrap_plots(list(p1, p2), nrow = 2, heights = c(4, 1))
    ggsave("results/105.Distribution/ClusterDistrBar_rev.png", p, width = 10, height = 13)

    df1 <- ClusterDistrBar(origin = sc_final$orig.ident, cluster = sc_final$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")
    p1 <- tidyplot(df1, x = origin, y = percentage, color = cluster) %>%
        add_areastack_absolute() %>%
        adjust_colors(c(colors_discrete_metro, colors_discrete_seaside))
    df2 <- ClusterDistrBar(origin = sc_final$group, cluster = sc_final$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")
    p2 <- tidyplot(df2, x = origin, y = percentage, color = cluster) %>%
        add_areastack_absolute() %>%
        adjust_colors(c(colors_discrete_metro, colors_discrete_seaside))
    p <- wrap_plots(list(p1, p2), nrow = 1)
    save_plot(p, "results/105.Distribution/ClusterDistrBar_area.png")
    # combine df1 and df2 in one excel, two sheets
    writexl::write_xlsx(list(df1 = df1, df2 = df2), "results/105.Distribution/ClusterDistr.xlsx")
    "results/105.Distribution"
}

# sc_final_filt <- sc_final %>%
#     filter(seurat_clusters %in% c(11, 13))
# df1 <- ClusterDistrBar(origin = sc_final_filt$orig.ident, cluster = sc_final_filt$seurat_clusters, plot = FALSE) %>%
#     as.data.frame() %>%
#     rownames_to_column("cluster") %>%
#     pivot_longer(-cluster, names_to = "origin", values_to = "percentage")
# p1 <- tidyplot(df1, x = origin, y = percentage, color = cluster) %>%
#     add_areastack_absolute() %>%
#     adjust_colors(c(colors_discrete_metro, colors_discrete_seaside))
# df2 <- ClusterDistrBar(origin = sc_final_filt$group, cluster = sc_final_filt$seurat_clusters, plot = FALSE) %>%
#     as.data.frame() %>%
#     rownames_to_column("cluster") %>%
#     pivot_longer(-cluster, names_to = "origin", values_to = "percentage")
# p2 <- tidyplot(df2, x = origin, y = percentage, color = cluster) %>%
#     add_areastack_absolute() %>%
#     adjust_colors(c(colors_discrete_metro, colors_discrete_seaside))
# p <- patchwork::wrap_plots(list(p1, p2), nrow = 1)
# save_plot(p, "results/105.Distribution/ClusterDistrBar_area_cluster.png")
