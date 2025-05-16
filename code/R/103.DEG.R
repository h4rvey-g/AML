DEG_annotation <- function(sc_annotate, cell_types) {
    # Compare between tumor and normal for each cell type
    sc_pseudo <- AggregateExpression(sc_annotate, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    for (cell_type in cell_types) {
        tryCatch(
            {
                ident1 <- paste0(cell_type, "_tumor")
                ident2 <- paste0(cell_type, "_normal")
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
                    filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                    arrange(desc(abs(avg_log2FC)))
                write_tsv(deg, paste0("results/103.DEG/tumor_vs_normal/DEG_", gsub("/", "_", cell_type), "_DESeq2.tsv"))
            },
            error = function(e) {
                message("Error in DESeq2 analysis for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )
        tryCatch(
            {
                ident1 <- "tumor"
                ident2 <- "normal"
                Idents(sc_annotate) <- "group"
                deg <- FindMarkers(
                    sc_annotate %>%
                        filter(cell_type_dtl == cell_type),
                    ident.1 = ident1,
                    ident.2 = ident2,
                    test.use = "MAST",
                    min.cells.group = 2
                ) %>%
                    Add_Pct_Diff() %>%
                    rownames_to_column("gene") %>%
                    as_tibble() %>%
                    filter(p_val_adj < 0.05) %>%
                    filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                    arrange(desc(abs(avg_log2FC)))
                write_tsv(deg, paste0("results/103.DEG/tumor_vs_normal/DEG_", gsub("/", "_", cell_type), "_MAST.tsv"))
            },
            error = function(e) {
                message("Error in MAST analysis for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )
    }

    # Compare each cell type with all others
    sc_pseudo2 <- AggregateExpression(sc_annotate, assays = "RNA", return.seurat = TRUE, group.by = c("dataset", "cell_type_dtl"))
    Idents(sc_pseudo2) <- "cell_type_dtl"
    dir.create("results/103.DEG/one_vs_rest", showWarnings = FALSE, recursive = TRUE)

    for (cell_type in cell_types) {
        tryCatch(
            {
                deg <- FindMarkers(
                    sc_pseudo2,
                    ident.1 = cell_type,
                    test.use = "DESeq2",
                    min.cells.group = 2
                ) %>%
                    Add_Pct_Diff() %>%
                    rownames_to_column("gene") %>%
                    as_tibble() %>%
                    filter(p_val_adj < 0.05) %>%
                    filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                    arrange(desc(abs(avg_log2FC)))
                write_tsv(deg, sprintf(
                    "results/103.DEG/one_vs_rest/%s_vs_rest_DEG.tsv",
                    gsub("/", "_", cell_type)
                ))
            },
            error = function(e) {
                message("Error in DESeq2 analysis for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )
        tryCatch(
            {
                Idents(sc_annotate) <- "cell_type_dtl"
                deg <- FindMarkers(
                    sc_annotate,
                    ident.1 = cell_type,
                    test.use = "MAST",
                    min.cells.group = 2
                ) %>%
                    Add_Pct_Diff() %>%
                    rownames_to_column("gene") %>%
                    as_tibble() %>%
                    filter(p_val_adj < 0.05) %>%
                    filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                    arrange(desc(abs(avg_log2FC)))

                write_tsv(deg, sprintf(
                    "results/103.DEG/one_vs_rest/%s_vs_rest_MAST_DEG.tsv",
                    gsub("/", "_", cell_type)
                ))
            },
            error = function(e) {
                message("Error in MAST analysis for ", cell_type, ": ", e$message)
                return(NULL)
            }
        )
    }

    "results/103.DEG"
}

run_GSEA <- function(sc_final) {
    options(max.print = 12, spe = "human")
    parents <- GeneSetAnalysisGO() %>% names()
    # run analysis on every parent category
    for (parent in parents) {
        cat("Running GSEA for", parent, "\n")
        sc_final <- GeneSetAnalysisGO(sc_final, nCores = 20, parent = parent)
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
                top.n = 20,
                color = "p"
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
            top.n = 20,
            color = "p"
        )
        ggsave(paste0("results/104.GSEA/key/Waterfall_hallmark50_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
    }
    "results/104.GSEA"
}

run_GSEA_Plasticity <- function(sc_filtered, parent) {
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
            top.n = 20,
            color = "p"
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
            top.n = 20,
            color = "p"
        )
        ggsave(paste0(
            "results/104.GSEA/key/one_to_one/Waterfall_GO_",
            parent, "_", gsub("/", "_", combo[1]), "_vs_", gsub("/", "_", combo[2]), ".png"
        ), p, width = 15, height = 15)
    }
    "results/104.GSEA"
}
compare_cell_types_hallmark50 <- function(sc_filtered) {
    options(max.print = 12, spe = "human")
    sc_filtered <- GeneSetAnalysis(sc_filtered, genesets = hall50$human, nCores = 30)
    matr <- sc_filtered@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
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

    # Original tumor vs normal comparison
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
            top.n = 20,
            color = "p"
        )
        ggsave(paste0("results/104.GSEA/key/t_vs_n/Waterfall_hallmark50_", gsub("/", "_", cell_type), ".png"), p, width = 15, height = 15)
    }

    # New: Compare between each pair of cell types
    cell_type_combinations <- combn(cell_types, 2, simplify = FALSE)
    Idents(sc_filtered) <- "cell_type_dtl"
    for (combo in cell_type_combinations) {
        ident1 <- combo[1]
        ident2 <- combo[2]
        p <- WaterfallPlot(
            matr,
            f = sc_filtered$cell_type_dtl,
            ident.1 = ident1,
            ident.2 = ident2,
            top.n = 20,
            color = "p"
        )
        ggsave(paste0(
            "results/104.GSEA/key/one_to_one/Waterfall_hallmark50_",
            gsub("/", "_", combo[1]), "_vs_", gsub("/", "_", combo[2]), ".png"
        ), p, width = 15, height = 15)
    }

    "results/104.GSEA"
}

test_distribution <- function(sc_final) {
    # Create directory
    dir.create("results/105.Distribution", showWarnings = FALSE, recursive = TRUE)

    # Add required metadata
    sc_final <- sc_final %>%
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
    composition_test <- sc_final %>%
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
    ggsave("results/105.Distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/105.Distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/105.Distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save test results
    write_tsv(
        composition_test %>% select(-count_data),
        "results/105.Distribution/test_results.tsv"
    )

    return(composition_test)
}

plot_distribution <- function(sc_final, composition_test) {
    # Calculate distribution percentages and cell counts
    df1 <- ClusterDistrBar(origin = sc_final$dataset, cluster = sc_final$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")

    cell_counts <- sc_final@meta.data %>%
        group_by(cell_type_dtl, dataset) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(cell_type_dtl) %>%
        summarise(
            N_count = sum(count[str_starts(dataset, "N")]),
            T_count = sum(count[str_starts(dataset, "T")]),
            .groups = "drop"
        )

    # Add significance stars
    composition_test_res <- composition_test %>%
        filter(factor == "group", c_FDR < 0.05) %>%
        mutate(stars = case_when(
            c_FDR < 0.001 ~ "***",
            c_FDR < 0.01 ~ "**",
            c_FDR < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        select(cell_group, stars)

    # Add cell counts and stars to cluster labels
    df1 <- df1 %>%
        mutate(cluster = map_chr(str_extract(cluster, ".*"), ~ {
            stars <- composition_test_res$stars[composition_test_res$cell_group == .x]
            count_info <- sprintf(
                "(N=%d,T=%d)",
                cell_counts$N_count[cell_counts$cell_type_dtl == .x],
                cell_counts$T_count[cell_counts$cell_type_dtl == .x]
            )
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

    save_plot(p, "results/105.Distribution/cell_distribution.png")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/105.Distribution/cell_distribution.tsv"
    )

    return("results/105.Distribution")
}

calculate_DEG_Fib_vs_Others <- function(sc_filtered) {
    sc_pseudo <- AggregateExpression(sc_filtered, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    comparisons <- list(c("Fib_tumor", "PSC_normal"), c("Fib_tumor", "MSC_normal"))
    for (comp in comparisons) {
        ident1 <- comp[1]
        ident2 <- comp[2]
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
        write_tsv(
            deg,
            paste0("results/103.DEG/DEG_", ident1, "_vs_", ident2, ".tsv")
        )
    }
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
