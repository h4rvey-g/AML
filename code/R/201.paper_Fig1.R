paper_final_annotation <- function(sc_final, sc_adipo) {
    dir.create("results/109.paper/Fig1", showWarnings = FALSE)
    # Update the cell_type_dtl in sc_final from sc_adipo
    cell_mapping <- setNames(sc_adipo$cell_type_dtl, colnames(sc_adipo))
    matching_cells <- intersect(Cells(sc_final), Cells(sc_adipo))
    sc_final$cell_type_dtl[matching_cells] <- cell_mapping[matching_cells]
    # Define a custom color palette for all cell types that will be used consistently
    cell_types <- unique(sc_final$cell_type_dtl) %>% sort()
    set.seed(42) # For reproducible color generation
    custom_colors <- c(
        "#dc467fcc", "#5cb248cc", "#b25ac9cc", "#b4b540cc", "#6b82d9cc", "#ce4f32cc",
        "#49afcfcc", "#cd8f44cc", "#6f56a2cc", "#54a778cc", "#a34673cc", "#767933cc",
        "#d489c4cc", "#c2685fcc", "#4b9fe3cc"
    )
    names(custom_colors) <- cell_types

    # Get the subset of cell types present in tumor and normal samples
    tumor_cell_types <- unique((sc_final %>% filter(group == "tumor"))$cell_type_dtl)
    normal_cell_types <- unique((sc_final %>% filter(group == "normal"))$cell_type_dtl)

    # Create subset color palettes (maintaining the same colors for each cell type)
    tumor_colors <- custom_colors[names(custom_colors) %in% tumor_cell_types]
    normal_colors <- custom_colors[names(custom_colors) %in% normal_cell_types]

    # Plot all cells with custom colors
    p <- DimPlot2(sc_final,
        reduction = "umap_integrated",
        group.by = "cell_type_dtl",
        label = TRUE,
        cols = custom_colors,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/dimplot_all.tiff", p, width = 7, height = 7, device = "tiff", compression = "lzw")

    # Plot tumor cells with the same custom colors (matching subset)
    p_tumor <- DimPlot2(sc_final %>% filter(group == "tumor"),
        reduction = "umap_integrated",
        group.by = "cell_type_dtl",
        label = TRUE,
        cols = tumor_colors,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white")) +
        ggtitle("Tumor")
    ggsave("results/109.paper/Fig1/dimplot_tumor.tiff", p_tumor, width = 7, height = 7, device = "tiff", compression = "lzw")

    # Plot normal cells with the same custom colors (matching subset)
    p_normal <- DimPlot2(sc_final %>% filter(group == "normal"),
        reduction = "umap_integrated",
        group.by = "cell_type_dtl",
        label = TRUE,
        cols = normal_colors,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white")) +
        ggtitle("Normal")
    ggsave("results/109.paper/Fig1/dimplot_normal.tiff", p_normal, width = 7, height = 7, device = "tiff", compression = "lzw")
    markers <- list(
        "ZG" = c("NR5A1", "DACH1", "CYP11B2"),
        "ZG/ZF" = c("NR5A1", "NOV", "NCAM1"),
        "ZR" = c("NR5A1", "CYB5A", "SULT2A1"),
        "CLC" = c("TH", "CHGA"),
        "Endo" = c("PECAM1", "EMCN"),
        "MSC/Progenitor" = c("PDGFRA", "ENG", "NT5E", "MGP"),
        "PSC" = c("PDGFRB", "CSPG4", "RGS5", "ACTA2", "MCAM"),
        "CAR-like" = c("CXCL12", "LEPR", "FOXC1"),
        "Adipo" = c("ADIPOQ", "FABP4", "PPARG"),
        "Tcell" = c("CD3D", "CD3E", "TRBC1"),
        "Bcell" = c("CD19", "CD79A", "MS4A1"),
        "Myeloid" = c("ITGAM", "CD33"),
        "Neural" = c("PTPRD", "RBFOX1", "HCN1", "SYT16"),
        "Plasma" = c("CD38", "SDC1", "IGHG1"),
        "LEC" = c("PDPN", "PROX1", "NR2F2")
        # "Eryth" = c("HBB", "GYPA", "SLC4A1")
    )

    # Combined plot for all samples
    toplot <- CalcStats(sc_final,
        features = markers %>% unlist(),
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist())
    gene_groups <- gene_groups[rownames(toplot)]
    p_all <- Heatmap(t(toplot),
        lab_fill = "zscore",
        facet_col = gene_groups
    ) +
        ggtitle("All Samples") +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/109.paper/Fig1/final_annotation_heatmap.tiff", p_all, width = 14, height = 7)

    # get Number of cells of each cell type in normal and tumor samples, save to tsv
    cell_type_counts <- sc_final %>%
        group_by(group, cell_type_dtl) %>%
        summarise(cell_count = n(), .groups = "drop")

    write_tsv(cell_type_counts, "results/109.paper/Fig1/cell_type_counts.tsv")

    sc_final
}

paper_cell_distribution <- function(sc_final) {
    # Create a pie chart showing distribution of cells in tumor vs normal
    cell_distribution <- sc_final %>%
        group_by(group) %>%
        summarise(cell_count = n()) %>%
        mutate(percentage = cell_count / sum(cell_count) * 100)

    # Create a cell count summary by sample
    sample_distribution <- sc_final %>%
        group_by(orig.ident) %>%
        summarise(cell_count = n()) %>%
        mutate(
            percent = round(cell_count / sum(cell_count) * 100, 2)
        ) %>%
        arrange(orig.ident)

    # Save the data as TSV
    write_tsv(sample_distribution, "results/109.paper/Fig1/sample_cell_counts.tsv")

    # Create a connected stacked bar chart for cell type distribution between tumor and normal
    # Calculate cell type percentages for each group
    cell_type_dist <- sc_final %>%
        group_by(group, cell_type_dtl) %>%
        summarise(cell_count = n(), .groups = "drop") %>%
        group_by(group) %>%
        mutate(percentage = cell_count / sum(cell_count) * 100) %>%
        arrange(group, desc(percentage)) %>%
        ungroup()

    # Define a custom color palette for all cell types that will be used consistently
    cell_types <- unique(sc_final$cell_type_dtl) %>% sort()
    set.seed(42) # For reproducible color generation
    custom_colors <- c(
        "#dc467fcc", "#5cb248cc", "#b25ac9cc", "#b4b540cc", "#6b82d9cc", "#ce4f32cc",
        "#49afcfcc", "#cd8f44cc", "#6f56a2cc", "#54a778cc", "#a34673cc", "#767933cc",
        "#d489c4cc", "#c2685fcc"
    )
    names(custom_colors) <- cell_types

    # Get the subset of cell types present in tumor and normal samples
    tumor_cell_types <- unique((sc_final %>% filter(group == "tumor"))$cell_type_dtl)
    normal_cell_types <- unique((sc_final %>% filter(group == "normal"))$cell_type_dtl)

    # Create subset color palettes (maintaining the same colors for each cell type)
    tumor_colors <- custom_colors[names(custom_colors) %in% tumor_cell_types]
    normal_colors <- custom_colors[names(custom_colors) %in% normal_cell_types]

    # Calculate the percentage change for each cell type (tumor - normal)
    cell_type_changes <- cell_type_dist %>%
        select(group, cell_type_dtl, percentage) %>%
        pivot_wider(names_from = group, values_from = percentage, values_fill = 0) %>%
        mutate(
            percentage_change = tumor - normal,
            abs_change = abs(percentage_change)
        ) %>%
        # Order by percentage change (ascending, so largest negative change is first)
        arrange(percentage_change)

    # Create ordered factor based on percentage change
    ordered_cell_types <- cell_type_changes$cell_type_dtl

    # Process data for the connected bar chart using the new ordering
    connector_data <- data.frame()

    for (g in c("normal", "tumor")) {
        group_data <- cell_type_dist %>%
            filter(group == g) %>%
            # Use the ordering based on percentage change
            mutate(cell_type_dtl = factor(cell_type_dtl, levels = ordered_cell_types)) %>%
            arrange(cell_type_dtl)

        # Add any missing cell types with 0 percentage
        missing_types <- setdiff(ordered_cell_types, group_data$cell_type_dtl)
        if (length(missing_types) > 0) {
            missing_data <- data.frame(
                group = g,
                cell_type_dtl = factor(missing_types, levels = ordered_cell_types),
                cell_count = 0,
                percentage = 0
            )
            group_data <- bind_rows(group_data, missing_data)
        }

        # Calculate y positions (cumulative percentages)
        group_data <- group_data %>%
            arrange(cell_type_dtl) %>%
            mutate(
                y_bottom = lag(cumsum(percentage), default = 0),
                y_top = cumsum(percentage),
                y_mid = (y_bottom + y_top) / 2
            )

        connector_data <- bind_rows(connector_data, group_data)
    }

    # Find cell types that exist in both normal and tumor
    cell_types_in_both <- intersect(
        unique((connector_data %>% filter(group == "normal", percentage > 0))$cell_type_dtl),
        unique((connector_data %>% filter(group == "tumor", percentage > 0))$cell_type_dtl)
    )

    # Create the connected stacked bar plot
    p_connected <- ggplot() +
        # Add bars (width reduced from 0.8 to 0.4)
        geom_rect(
            data = connector_data,
            aes(
                xmin = as.numeric(factor(group)) - 0.2, # Changed from -0.4 to -0.2
                xmax = as.numeric(factor(group)) + 0.2, # Changed from +0.4 to +0.2
                ymin = y_bottom,
                ymax = y_top,
                fill = cell_type_dtl
            )
        ) +
        # Add connecting lines ONLY for cell types that exist in both normal and tumor
        # Also update x coordinates for the connecting lines
        geom_segment(
            data = connector_data %>%
                select(group, cell_type_dtl, y_bottom, y_top) %>%
                pivot_wider(names_from = group, values_from = c(y_bottom, y_top)) %>%
                filter(cell_type_dtl %in% cell_types_in_both) %>% # Only include cell types in both
                filter(!is.na(y_bottom_normal) & !is.na(y_bottom_tumor)),
            aes(
                x = 1 + 0.2, xend = 2 - 0.2, # Changed from +0.4/-0.4 to +0.2/-0.2
                y = y_bottom_normal, yend = y_bottom_tumor
            ),
            color = "grey50", linetype = "dashed", size = 0.5
        ) +
        geom_segment(
            data = connector_data %>%
                select(group, cell_type_dtl, y_bottom, y_top) %>%
                pivot_wider(names_from = group, values_from = c(y_bottom, y_top)) %>%
                filter(cell_type_dtl %in% cell_types_in_both) %>% # Only include cell types in both
                filter(!is.na(y_top_normal) & !is.na(y_top_tumor)),
            aes(
                x = 1 + 0.2, xend = 2 - 0.2, # Changed from +0.4/-0.4 to +0.2/-0.2
                y = y_top_normal, yend = y_top_tumor
            ),
            color = "grey50", linetype = "dashed", size = 0.5
        ) +
        # Add text labels in the center of each bar segment
        geom_text(
            data = connector_data %>%
                filter(percentage >= 2), # Only label segments that are at least 2% (adjust as needed)
            aes(
                x = as.numeric(factor(group)),
                y = (y_bottom + y_top) / 2, # Center of the bar segment
                label = cell_type_dtl,
                color = "white" # All text white
            ),
            fontface = "bold", # Make all text bold
            size = 3, # Adjust text size as needed
            check_overlap = TRUE # Helps prevent overlapping labels
        ) +
        # Labels and theme
        scale_x_continuous(
            breaks = 1:2, labels = c("normal", "tumor"),
            limits = c(0.5, 2.5)
        ) +
        scale_fill_manual(values = custom_colors) +
        scale_color_identity() + # Use the colors directly specified in geom_text
        labs(
            title = "Cell Type Distribution Changes Between Normal and Tumor",
            x = "Sample Type",
            y = "Percentage (%)"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none", # Remove legend
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white")
        )
    ggsave("results/109.paper/Fig1/cell_type_distribution_connected.tiff", p_connected, width = 10, height = 8)

    # Add a percentage change plot to visualize the shifts
    p_change <- ggplot(cell_type_changes, aes(x = percentage_change, y = cell_type_dtl, fill = percentage_change > 0)) +
        geom_bar(stat = "identity") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_manual(
            values = c("TRUE" = "#E64B35FF", "FALSE" = "#4DBBD5FF"),
            labels = c("TRUE" = "Increase", "FALSE" = "Decrease"),
            name = "Change Direction"
        ) +
        labs(
            title = "Change in Cell Type Proportions (Tumor vs Normal)",
            x = "Percentage Point Change",
            y = "Cell Type"
        ) +
        theme_minimal() +
        theme(
            legend.position = "top",
            plot.background = element_rect(fill = "white")
        )
    ggsave("results/109.paper/Fig1/cell_type_percentage_change.tiff", p_change, width = 8, height = 6)

    # Original stacked bar plot with new ordering
    p_bar <- ggplot(connector_data, aes(x = group, y = percentage, fill = cell_type_dtl)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = custom_colors) +
        labs(
            title = "Cell Type Distribution by Sample",
            x = "Sample Type",
            y = "Percentage",
            fill = "Cell Type"
        ) +
        theme_bw() +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/109.paper/Fig1/cell_type_distribution_bar.tiff", p_bar, width = 8, height = 6)
    "results/109.paper/Fig1"
}