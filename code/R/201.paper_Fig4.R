paper_tcell_annotate <- function(sc_tcell) {
    # Create directory if it doesn't exist
    dir.create("results/109.paper/Fig4", showWarnings = FALSE, recursive = TRUE)
    p <- DimPlot2(sc_tcell,
        reduction = "umap_tcell",
        group.by = "cell_type_dtl",
        label = TRUE,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig4/dimplot_tcell.tiff", p, width = 7, height = 5)

    markers <- list(
        "proliferation" = c("MKI67"),
        "Tex" = c("SRGN", "DUSP4"),
        "Th1" = c("CCL5", "CXCR3"),
        "Th17" = c("IL7R", "CCR6"),
        "Th22" = c("IL22", "AHR"),
        "Tn" = c("SELL", "TCF7", "CCR7"),
        "Treg" = c("FOXP3", "IL2RA")
    )
    # calculate for CD4, CD8, and DN separately
    toplot_CD4 <- CalcStats(
        sc_tcell %>%
            # filter cell_type_dtl starts with CD4
            filter(grepl("^CD4", cell_type_dtl)),
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    gene_groups_CD4 <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_CD4)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p1 <- Heatmap(toplot_CD4 %>% t(), lab_fill = "zscore", facet_col = gene_groups_CD4, text.size = 5) +
        labs(title = "CD4 T Cell Subtypes") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig4/heatmap_CD4_tcell.tiff", p1, width = 15, height = 8)

    markers <- list(
        "proliferation" = c("MKI67", "HMGB2", "HMGB1"),
        "quiescence" = c("RUNX1"),
        "Tem" = c("CST7", "EOMES", "GZMK", "GZMA"),
        "Temra" = c("GNLY", "FGFBP2"),
        "Tex" = c("LAG3", "TIGIT", "HAVCR2", "PDCD1"),
        "Tn" = c("LEF1", "SELL", "TCF7", "CCR7")
    )
    toplot_CD8 <- CalcStats(
        sc_tcell %>%
            # filter cell_type_dtl starts with CD4
            filter(grepl("^CD8", cell_type_dtl)),
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    gene_groups_CD8 <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_CD8)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p2 <- Heatmap(toplot_CD8 %>% t(), lab_fill = "zscore", facet_col = gene_groups_CD8, text.size = 5) +
        labs(title = "CD8 T Cell Subtypes") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig4/heatmap_CD8_tcell.tiff", p2, width = 15, height = 8)

    markers <- list(
        # "CD3" = c("CD3D", "CD3E"),
        "CD4/8" = c("CD4", "CD8A", "CD8B"),
        "gdT_1" = c("TRDV1"),
        "gdT_2" = c("TRDV2"),
        "gdT_constant" = c("TRGC1", "TRGC2", "TRDC")
    )
    sc_tcell <- sc_tcell %>%
        mutate(cell_type_dtl_DN = case_when(
            grepl("CD4", cell_type_dtl) ~ "CD4_Tcell",
            grepl("CD8", cell_type_dtl) ~ "CD8_Tcell",
            TRUE ~ cell_type_dtl
        ))
    toplot_DN <- CalcStats(
        sc_tcell,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl_DN"
    )
    gene_groups_DN <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_DN)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p3 <- Heatmap(toplot_DN %>% t(), lab_fill = "zscore", facet_col = gene_groups_DN) +
        labs(title = "CD4+, CD8+, DN T Cell Subtypes (z-score)") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig4/heatmap_broad_tcell.tiff", p3, width = 15, height = 8)

    markers <- list(
        "CD3" = c("CD3D", "CD3E"),
        "CD4/8" = c("CD4", "CD8A", "CD8B"),
        "gdT_1" = c("TRDV1"),
        "gdT_2" = c("TRDV2"),
        "gdT_constant" = c("TRGC1", "TRGC2", "TRDC")
    )
    # Creating the same plot but with mean values
    toplot_DN_mean <- CalcStats(
        sc_tcell,
        features = markers %>% unlist(),
        method = "mean", group.by = "cell_type_dtl_DN"
    )
    gene_groups_DN <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_DN_mean)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p3_mean <- Heatmap(toplot_DN_mean %>% t(), lab_fill = "mean expression", facet_col = gene_groups_DN) +
        labs(title = "CD4+, CD8+, DN T Cell Subtypes (mean expression)") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig4/heatmap_broad_tcell_mean.tiff", p3_mean, width = 15, height = 8)

    p <- p1 +
        p2 +
        plot_layout(ncol = 1)
    ggsave("results/109.paper/Fig4/heatmap_combined_tcell.tiff", p, width = 15, height = 15)
    "results/109.paper/Fig4"
}
compare_tcell_bcell_counts <- function(sc_final, sc_tcell) {
    # Create directory if it doesn't exist
    dir.create("results/109.paper/Fig4", showWarnings = FALSE, recursive = TRUE)

    # Filter for tumor samples (orig.ident starts with T)
    sc_tumor <- sc_final[, grepl("^T", sc_final$orig.ident)]

    # Check if we have tumor samples
    if (ncol(sc_tumor) == 0) {
        message("No tumor samples found with orig.ident starting with 'T'")
        return(NULL)
    }

    # First part: Analyze Tcell and Bcell distribution in sc_final
    cell_counts <- data.frame(
        sample = sc_tumor$orig.ident,
        cell_type_dtl = sc_tumor$cell_type_dtl,
        stringsAsFactors = FALSE
    )

    # Filter to only include T cells and B cells
    cell_counts <- cell_counts %>%
        filter(grepl("^Tcell|^Bcell", cell_type_dtl))

    # Summarize T and B cell counts
    broad_cell_summary <- cell_counts %>%
        group_by(sample, cell_type_dtl) %>%
        summarize(count = n(), .groups = "drop") %>%
        mutate(percentage = count / sum(count) * 100, .by = sample) %>%
        arrange(sample, desc(count))

    # Write detailed T and B cell type information to TSV
    write.table(
        broad_cell_summary,
        file = "results/109.paper/Fig4/tcell_bcell_types.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Create visualization for broad cell types
    p_broad <- ggplot(broad_cell_summary) +
        geom_bar(aes(x = sample, y = percentage, fill = cell_type_dtl),
            position = "stack",
            stat = "identity"
        ) +
        labs(
            title = "Broad cell type distribution by tumor sample",
            x = "Sample",
            y = "Percentage (%)",
            fill = "Cell type"
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )

    # Save broad cell type plot
    ggsave("results/109.paper/Fig4/broad_cell_types_distribution.tiff",
        p_broad,
        width = 12,
        height = 8
    )

    # Second part: Focus specifically on T cells and B cells
    # Define cell type patterns for T cells and B cells
    t_cell_pattern <- "^Tcell"
    b_cell_pattern <- "^Bcell"

    # Calculate T and B cell counts by sample
    cell_tb_summary <- cell_counts %>%
        group_by(sample) %>%
        summarize(
            total_cells = n(),
            t_cells = sum(grepl(t_cell_pattern, cell_type_dtl)),
            b_cells = sum(grepl(b_cell_pattern, cell_type_dtl)),
            t_cell_percent = t_cells / total_cells * 100,
            b_cell_percent = b_cells / total_cells * 100
        ) %>%
        arrange(sample)

    # Write to TSV
    write.table(
        cell_tb_summary,
        file = "results/109.paper/Fig4/tcell_bcell_counts.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Create visualization - count plot for T and B cells
    p1 <- ggplot(cell_tb_summary %>%
        pivot_longer(
            cols = c("t_cells", "b_cells"),
            names_to = "cell_type_dtl",
            values_to = "count"
        )) +
        geom_bar(aes(x = sample, y = count, fill = cell_type_dtl),
            position = "dodge",
            stat = "identity"
        ) +
        labs(
            title = "T cell and B cell counts by tumor sample",
            x = "Sample",
            y = "Cell count",
            fill = "Cell Type"
        ) +
        scale_fill_manual(
            values = c("t_cells" = "#FF7F0E", "b_cells" = "#1F77B4"),
            labels = c("t_cells" = "T cells", "b_cells" = "B cells")
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    # Third part: T cell subtypes from sc_tcell
    # Get T cell barcodes from sc_final to match with sc_tcell
    t_cell_barcodes <- colnames(sc_tumor)[grepl(t_cell_pattern, sc_tumor$cell_type)]
    
    # Filter sc_tcell to include only T cells from tumor samples
    sc_tcell_tumor <- sc_tcell[, colnames(sc_tcell) %in% t_cell_barcodes]
    
    # Create a data frame with sample and T cell subtypes
    t_cell_subtypes_df <- data.frame(
        sample = sc_tcell_tumor$orig.ident,
        cell_type_dtl = sc_tcell_tumor$cell_type_dtl,
        stringsAsFactors = FALSE
    )

    # Summarize T cell subtypes
    t_cell_subtypes <- t_cell_subtypes_df %>%
        group_by(sample, cell_type_dtl) %>%
        summarize(count = n(), .groups = "drop") %>%
        mutate(percentage = count / sum(count) * 100, .by = sample) %>%
        arrange(sample, desc(count))

    # Write detailed T cell subtype information to TSV
    write.table(
        t_cell_subtypes,
        file = "results/109.paper/Fig4/tcell_subtypes.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Create visualization for T cell subtypes
    p_tcell_subtypes <- ggplot(t_cell_subtypes) +
        geom_bar(aes(x = sample, y = percentage, fill = cell_type_dtl),
            position = "stack",
            stat = "identity"
        ) +
        labs(
            title = "T cell subtype distribution by tumor sample",
            x = "Sample",
            y = "Percentage (%)",
            fill = "T cell subtype"
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
        )

    # Save T cell subtype plot
    ggsave("results/109.paper/Fig4/tcell_subtypes_distribution.tiff",
        p_tcell_subtypes,
        width = 12,
        height = 8
    )

    # Combine all plots
    combined_plot <- p_broad / p1 / p_tcell_subtypes + plot_layout(ncol = 1, heights = c(1, 0.8, 1))
    
    # Save combined plot
    ggsave("results/109.paper/Fig4/cell_type_distribution_combined.tiff",
        combined_plot,
        width = 12,
        height = 18
    )

    # Return the directory path
    return("results/109.paper/Fig4")
}
