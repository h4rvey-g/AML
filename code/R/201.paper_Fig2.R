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
paper_plot_DEG_Car_vs_MSC <- function(DEG_res) {
    # Create directory for results
    dir.create("results/109.paper/Fig2", showWarnings = FALSE, recursive = TRUE)
    
    # Get DEG results from pseudobulk analysis
    DEG_pseudobulk <- DEG_res$DEG_pseudobulk
    
    # Define genes to highlight with their corresponding categories
    highlight_genes <- list(
        # CAR-like cell markers (Myelolipoma)
        "CAR-like markers" = c("LEPR", "CXCL12", "THY1", "VCAM1"),
        # Adipogenesis markers
        "Adipogenesis" = c("PPARG", "LPL", "CFD"),
        # ECM and stromal support
        "Stromal/ECM" = c("COL1A1", "FBN1", "PTX3", "GAS6", "ANGPT1"),
        # Hematopoietic factors
        "Hematopoiesis" = c("RUNX1", "IKZF1"),
        # Normal adrenal MSC markers
        "Adrenal MSC" = c("NR2F1", "NR2F1-AS1", "HHIP", "PTCH1", "RSPO3", "MEG3", "PDGFRA")
    )
    
    # Combine all genes into a single vector
    all_genes_to_highlight <- unlist(highlight_genes)
    
    # Filter DEG results to include only our specified genes
    highlight_degs <- DEG_pseudobulk %>%
        filter(gene %in% all_genes_to_highlight)
    
    # Create a category column for the genes
    highlight_degs <- highlight_degs %>%
        mutate(category = sapply(gene, function(g) {
            for(cat_name in names(highlight_genes)) {
                if(g %in% highlight_genes[[cat_name]]) return(cat_name)
            }
            return("Other")
        }))
    
    # Assign direction (up in CAR-like or up in MSC)
    highlight_degs <- highlight_degs %>%
        mutate(direction = ifelse(avg_log2FC > 0, "Up in CAR-like cells", "Up in MSC"))
    
    # Create enhanced volcano plot with category-based legend
    # Create a volcano plot using ggplot2
    p_volcano <- ggplot(highlight_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = category), size = 3, alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c(
            "CAR-like markers" = "#E41A1C", 
            "Adipogenesis" = "#377EB8", 
            "Stromal/ECM" = "#4DAF4A", 
            "Hematopoiesis" = "#984EA3", 
            "Adrenal MSC" = "#FF7F00"
        )) +
        geom_text_repel(
            aes(label = gene, color = category),
            box.padding = 0.5,
            point.padding = 0.3,
            segment.color = "grey50",
            segment.size = 0.5,
            force = 2,
            max.overlaps = 30
        ) +
        labs(
            title = "DEG: CAR-like cells (Tumor) vs MSC (Normal)",
            subtitle = "Pseudobulk DESeq2 analysis",
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-value",
            color = "Gene Category"
        ) +
        theme_bw() +
        theme(
            plot.background = element_rect(fill = "white"),
            legend.position = "right",
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()
        )
    
    # Save the volcano plot
    ggsave("results/109.paper/Fig2/CAR_vs_MSC_volcano_plot.tiff", 
           p_volcano, 
           width = 10, height = 8)
    
    # Create a lollipop plot
    # Prepare the data
    lollipop_data <- highlight_degs %>%
        arrange(category, avg_log2FC)
    
    # Set factor levels for proper ordering
    lollipop_data$gene <- factor(lollipop_data$gene, levels = lollipop_data$gene)
    
    # Create lollipop plot with faceting by gene category
    p_lollipop <- ggplot(lollipop_data, aes(x = avg_log2FC, y = gene, color = category)) +
        geom_segment(aes(x = 0, xend = avg_log2FC, y = gene, yend = gene, 
            color = category), size = 1.5) +
        geom_point(aes(size = -log10(p_val_adj))) +
        facet_grid(category ~ ., scales = "free_y", space = "free") +
        scale_color_manual(values = c(
            "CAR-like markers" = "#E41A1C", 
            "Adipogenesis" = "#377EB8", 
            "Stromal/ECM" = "#4DAF4A", 
            "Hematopoiesis" = "#984EA3", 
            "Adrenal MSC" = "#FF7F00"
        )) +
        scale_size_continuous(range = c(3, 8), name = "-log10(p-adj)") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(
            title = "Differential Expression: CAR-like cells (Tumor) vs MSC (Normal)",
            x = "Log2 Fold Change",
            y = "",
            color = "Gene Category"
        ) +
        theme_bw() +
        theme(
            plot.background = element_rect(fill = "white"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "right",
            strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(face = "bold")
        )
    
    # Save the lollipop plot
    ggsave("results/109.paper/Fig2/CAR_vs_MSC_lollipop_plot.tiff", 
           p_lollipop, 
           width = 10, height = 12)
    # Return the file paths to the generated figures
    "results/109.paper/Fig2"
}