paper_sub_annotation_adipo <- function(sc_adipo) {
    # Create UMAP visualization with new annotations
    p_annotated_clusters <- DimPlot2(sc_adipo,
        reduction = "adipo_umap",
        group.by = "cell_type_dtl",
        label = TRUE
    ) +
        ggtitle("Annotated Adipocyte and Stromal cells") +
        theme(plot.background = element_rect(fill = "white"))

    # Save plot
    ggsave("results/109.paper/Fig2/dimplot_adipo.tiff", p_annotated_clusters, width = 8, height = 7)

    # Create dotplot of markers by annotated cell types
    marker_genes <- list(
        "MSC/Progenitor" = c("MGP", "PDGFRA", "ENG", "CD14", "THY1", "NT5E"),
        "Pericyte" = c("PDGFRB", "CSPG4", "RGS5", "ACTA2", "MCAM"),
        "CAR-like" = c("CXCL12", "LEPR", "FOXC1"),
        # "Negative Markers" = c("PECAM1", "CD34", "PTPRC", "EPCAM"),
        "Adipocyte" = c("ADIPOQ", "FABP4", "PPARG")
    )

    # Combine all markers into a single vector
    all_markers <- unlist(marker_genes)
    # Generate heatmap of marker expression by subtype
    toplot <- CalcStats(sc_adipo,
        features = all_markers,
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    gene_groups <- rep(names(marker_genes), lengths(marker_genes)) %>%
        setNames(marker_genes %>% unlist())
    gene_groups <- gene_groups[rownames(toplot)]

    p_heatmap <- Heatmap(t(toplot),
        lab_fill = "zscore",
        facet_col = gene_groups
    ) +
        ggtitle("Stromal cell type marker expression") +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/109.paper/Fig2/heatmap_adipo.tiff",
        p_heatmap,
        width = 14, height = 7
    )

    "results/109.paper/Fig2"
}
paper_plot_adipo_distribution <- function(sc_adipo, m_composition_test) {
    dir.create("results/109.paper/Fig2", showWarnings = FALSE, recursive = TRUE)

    # Calculate distribution percentages
    df1 <- ClusterDistrBar(origin = sc_adipo$dataset, cluster = sc_adipo$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")

    # Calculate cell counts
    cell_counts <- sc_adipo@meta.data %>%
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

    tidyplots::save_plot(p, "results/109.paper/Fig2/adipo_distribution.tiff")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/109.paper/Fig2/adipo_distribution.tsv"
    )

    return("results/109.paper/Fig2")
}
