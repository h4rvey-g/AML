paper_final_annotation <- function(sc_final, sc_mye_clust) {
    # get cells in sc_mye_clust where immune_subcluster == 7
    steroid_cells <- sc_mye_clust %>%
        filter(immune_subcluster == 7) %>%
        pull(cell)
    # assign steroidogenic to these cells
    sc_final <- sc_final %>% mutate(
        cell_type = if_else(cell %in% steroid_cells, "Steroidogenic", cell_type),
        cell_type_dtl = if_else(cell %in% steroid_cells, "ZF", cell_type_dtl)
    )
    dir.create("results/109.paper/Fig1", showWarnings = FALSE)
    p <- DimPlot2(sc_final, reduction = "umap_integrated", group.by = "cell_type_dtl", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/final_annotation.png", p, width = 7, height = 7)
    markers <- list(
        "ZG" = c("NR5A1", "DACH1", "CYP11B2"),
        "ZG/ZF" = c("NR5A1", "NOV", "NCAM1"),
        "ZR" = c("NR5A1", "CYB5A", "SULT2A1"),
        "CLC" = c("TH", "CHGA"),
        "Endo" = c("PECAM1", "EMCN"),
        "Fib" = c("COL1A1", "COL3A1", "THY1"),
        "PSC" = c("RGS5", "PDGFRB"),
        "Adipo" = c("ADIPOQ", "FABP4", "PPARG"),
        "Tcell" = c("CD3D", "CD3E", "TRBC1"),
        "Bcell" = c("CD19", "CD79A", "MS4A1"),
        "Myeloid" = c("ITGAM", "CD33"),
        "Plasma" = c("CD38", "SDC1", "IGHG1"),
        "LEC" = c("PDPN", "PROX1", "NR2F2"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1")
    )

    p <- DotPlot2(sc_final,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/final_annotation_dotplot.png", p, width = 7, height = 7)
    # Get expression stats for tumor samples
    toplot_tumor <- CalcStats(sc_final %>% filter(group == "tumor"),
        features = markers %>% unlist(),
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    # Get expression stats for normal samples
    toplot_normal <- CalcStats(sc_final %>% filter(group == "normal"),
        features = markers %>% unlist(),
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    # Create gene group annotations
    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist())

    # For tumor plot
    gene_groups_tumor <- gene_groups[rownames(toplot_tumor)]
    p_tumor <- Heatmap(t(toplot_tumor),
        lab_fill = "zscore",
        facet_col = gene_groups_tumor
    ) +
        ggtitle("Tumor") +
        theme(plot.background = element_rect(fill = "white"))

    # For normal plot
    gene_groups_normal <- gene_groups[rownames(toplot_normal)]
    p_normal <- Heatmap(t(toplot_normal),
        lab_fill = "zscore",
        facet_col = gene_groups_normal
    ) +
        ggtitle("Normal") +
        theme(plot.background = element_rect(fill = "white"))

    # Combined plot for all samples
    toplot <- CalcStats(sc_final,
        features = markers %>% unlist(),
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    gene_groups_all <- gene_groups[rownames(toplot)]
    p_all <- Heatmap(t(toplot),
        lab_fill = "zscore",
        facet_col = gene_groups_all
    ) +
        ggtitle("All Samples") +
        theme(plot.background = element_rect(fill = "white"))
    # Create a combined plot (vertically stacked)
    p_combined <- p_tumor / p_normal

    # Save all plots
    ggsave("results/109.paper/Fig1/final_annotation_heatmap_tumor.png", p_tumor, width = 14, height = 5)
    ggsave("results/109.paper/Fig1/final_annotation_heatmap_normal.png", p_normal, width = 14, height = 5)
    ggsave("results/109.paper/Fig1/final_annotation_heatmap_combined.png", p_combined, width = 14, height = 10)
    ggsave("results/109.paper/Fig1/final_annotation_heatmap.png", p_all, width = 14, height = 7)
    "results/109.paper/Fig1"
}
paper_clc_expression_by_sample <- function(sc_final) {
    # Create directory for results
    dir.create("results/109.paper/Fig1/CLC", recursive = TRUE, showWarnings = FALSE)

    # Define CLC markers
    clc_markers <- c("TH", "CHGA", "CHGB")
    p <- VlnPlot2(sc_final %>% filter(cell_type_dtl == "CLC"),
        features = clc_markers,
        group.by = "dataset",
        violin = FALSE, pt.style = "quasirandom", ncol = 1
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/CLC/clc_expression_by_sample_violinplot.png", p, width = 6, height = 15)
    p <- VlnPlot2(sc_final,
        features = clc_markers,
        group.by = "dataset",
        violin = FALSE, pt.style = "quasirandom", ncol = 1
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/CLC/clc_expression_by_sample_all_violinplot.png", p, width = 6, height = 15)
    p <- VlnPlot2(sc_final %>% filter(group == "tumor"),
        features = clc_markers,
        group.by = "cell_type_dtl",
        violin = FALSE, pt.style = "quasirandom", ncol = 1
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/CLC/clc_expression_by_sample_tumor_violinplot.png", p, width = 6, height = 15)
    p <- VlnPlot2(sc_final %>% filter(group == "normal"),
        features = clc_markers,
        group.by = "cell_type_dtl",
        violin = FALSE, pt.style = "quasirandom", ncol = 1
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig1/CLC/clc_expression_by_sample_normal_violinplot.png", p, width = 6, height = 15)
    # Find DEGs for CLC vs rest in normal and tumor samples separately
    # Create directory for DEG results
    deg_dir <- file.path("results/109.paper/Fig1/CLC", "DEGs")
    dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

    # Set identity to cell_type_dtl for comparison
    Idents(sc_final) <- "cell_type_dtl"

    # Find DEGs for CLC vs rest in normal samples
    message("Finding DEGs for CLC vs rest in normal samples...")
    normal_degs <- FindMarkers(
        sc_final %>% filter(group == "normal"),
        ident.1 = "CLC",
        test.use = "MAST",
        min.cells.group = 2
    ) %>%
        Add_Pct_Diff() %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
        arrange(desc(abs(avg_log2FC)))
    # Save normal DEGs
    write.csv(normal_degs, file.path(deg_dir, "CLC_vs_rest_normal_DEGs.csv"), row.names = FALSE)

    # Find DEGs for CLC vs rest in tumor samples
    message("Finding DEGs for CLC vs rest in tumor samples...")
    tumor_degs <- FindMarkers(
        sc_final %>% filter(group == "tumor"),
        ident.1 = "CLC",
        test.use = "MAST",
        min.cells.group = 2
    ) %>%
        Add_Pct_Diff() %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
        arrange(desc(abs(avg_log2FC)))

    # Save tumor DEGs
    write.csv(tumor_degs, file.path(deg_dir, "CLC_vs_rest_tumor_DEGs.csv"), row.names = FALSE)

    # Add pseudobulk analysis
    message("Performing pseudobulk analysis for more robust DEG identification...")

    # Create pseudobulk for normal samples
    normal_pseudo <- sc_final %>%
        filter(group == "normal") %>%
        AggregateExpression(
            assays = "RNA",
            return.seurat = TRUE,
            group.by = c("dataset", "cell_type_dtl")
        )

    # Set identity for comparison
    Idents(normal_pseudo) <- "cell_type_dtl"

    # Find DEGs using DESeq2 on pseudobulk data for normal samples
    normal_pseudo_degs <- FindMarkers(
        normal_pseudo,
        ident.1 = "CLC",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        arrange(desc(abs(avg_log2FC)))

    # Save normal pseudobulk DEGs
    write.csv(normal_pseudo_degs, file.path(deg_dir, "CLC_vs_rest_normal_pseudobulk_DEGs.csv"), row.names = FALSE)

    # Create pseudobulk for tumor samples
    tumor_pseudo <- sc_final %>%
        filter(group == "tumor") %>%
        AggregateExpression(
            assays = "RNA",
            return.seurat = TRUE,
            group.by = c("dataset", "cell_type_dtl")
        )

    # Set identity for comparison
    Idents(tumor_pseudo) <- "cell_type_dtl"

    # Find DEGs using DESeq2 on pseudobulk data for tumor samples
    tumor_pseudo_degs <- FindMarkers(
        tumor_pseudo,
        ident.1 = "CLC",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        arrange(desc(abs(avg_log2FC)))

    # Save tumor pseudobulk DEGs
    write.csv(tumor_pseudo_degs, file.path(deg_dir, "CLC_vs_rest_tumor_pseudobulk_DEGs.csv"), row.names = FALSE)

    # Save tumor DEGs
    # write.csv(tumor_degs, file.path(deg_dir, "CLC_vs_rest_tumor_DEGs.csv"), row.names = FALSE)

    # Return results directory
    return("results/109.paper/Fig1/CLC")
}
paper_myeloid_annotate <- function(sc_mye) {
    markers <- list(
        "PMP" = c("AXL", "LYVE1", "CD34", "HLA-DQB1"), # Tissue-Resident Macrophages
        "TREM2_LAM" = c("CD40", "STAT1", "TREM2"), # TREM2+ Macrophages
        "Foam" = c("CCL18", "PPARG", "PLIN2"), # Foam Macrophages
        "Mast" = c("TPSAB1", "LYZ", "CCR2"), # Mast Cells
        "Eosinophil" = c("IL5RA", "EPO", "CXCR2"), # Eosinophils
        "Dendritic" = c("CD83", "CD86") # Dendritic Cells
    )

    toplot <- CalcStats(sc_mye,
        features = markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "cell_type_dtl"
    )

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)]

    p2 <- Heatmap(toplot %>% t(), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig2/myeloid_annotation_heatmap.png", p2, width = 7, height = 5)
    p <- DotPlot2(sc_mye,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig2/myeloid_annotation_dotplot.png", p, width = 7, height = 7)
}

paper_myeloid_lipid_DEG <- function(sc_mye) {
    markers <- list(
        "Lipid metabolism" = c("GML", "PPARG", "MEG3", "SORBS2", "DHCR24", "SCARB1", "PNPLA3", "FADS2", "FADS1", "SCD5", "CD36"),
        "Immune/Inflammatory" = c("CCL18", "HMOX1", "VCAM1", "RHOH", "CD247"),
        "ECM Remodeling" = c("MMP19", "SDC2", "DCN")
    )
    p <- DotPlot2(sc_mye,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig2/myeloid_lipid_dotplot.png", p, width = 7, height = 7)

    p <- DotPlot2(sc_mye_select,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig2/myeloid_lipid_select_dotplot.png", p, width = 7, height = 7)

    # --- Step 2: Prepare the gene list ---
    # Flatten the list of markers into a single vector
    all_markers <- unlist(markers)

    # Create a data frame to store gene categories
    gene_categories <- data.frame(
        gene = all_markers,
        category = rep(names(markers), times = sapply(markers, length))
    )

    # --- Step 3 & 4: Calculate Average Expression, Percentage, and LogFC ---

    # Create pseudobulk and find DEGs
    # Set identities for comparison
    sc_mye$cell_type_group <- paste(sc_mye$cell_type_dtl, sc_mye$group, sep = "_")
    Idents(sc_mye) <- "cell_type_group"

    # Function to extract DEGs for a specific cell type
    get_cell_type_degs <- function(cell_type) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")

        if ((ident1 %in% Idents(sc_mye) && ident2 %in% Idents(sc_mye))) {
            FindMarkers(
                sc_mye,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                min.cells.group = 2
            ) %>%
                Add_Pct_Diff() %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                mutate(cell_type = cell_type)
        } else {
            warning(paste("Either", ident1, "or", ident2, "not found in the dataset"))
            # Create an empty tibble with required structure
            tibble(
                gene = character(), avg_log2FC = numeric(),
                pct.1 = numeric(), pct.2 = numeric(), cell_type = cell_type
            )
        }
    }

    # Get DEGs for both cell types
    deg_trem2 <- get_cell_type_degs("TREM2_LAM")
    deg_foam <- get_cell_type_degs("Foam")

    # Combine DEGs from both cell types
    all_degs <- bind_rows(deg_trem2, deg_foam)

    # Filter to our markers of interest and prepare data for plotting
    plot_data <- all_degs %>%
        filter(gene %in% all_markers) %>%
        mutate(
            logFC = avg_log2FC,
            pct_expr_tumor = pct.1 * 100,
            pct_expr_normal = pct.2 * 100
        )

    # Add category information
    plot_data <- merge(plot_data, gene_categories, by = "gene")

    # Convert to long format for plotting
    plot_data_long <- pivot_longer(
        plot_data,
        cols = c(pct_expr_tumor, pct_expr_normal),
        names_to = "condition",
        values_to = "pct_expr"
    )

    # Clean up condition names
    plot_data_long$condition <- gsub("pct_expr_", "", plot_data_long$condition)

    # Order genes within each category by logFC
    plot_data_long <- plot_data_long %>%
        group_by(cell_type, category) %>%
        mutate(gene = factor(gene, levels = unique(gene[order(logFC)]))) %>%
        ungroup()

    # Set category order
    plot_data_long$category <- factor(plot_data_long$category, levels = names(markers))

    # Sort the data frame so that dots with lower pct_expr appear later (on top)
    plot_data_long <- plot_data_long %>%
        arrange(gene, category, desc(pct_expr))

    # Create plots for each cell type
    for (ct in c("TREM2_LAM", "Foam")) {
        cell_data <- plot_data_long %>% filter(cell_type == ct)

        # Skip if no data for this cell type
        if (nrow(cell_data) == 0) {
            warning(paste("No DEG data available for", ct))
            next
        }

        p <- ggplot(cell_data, aes(x = logFC, y = gene, size = pct_expr, color = condition)) +
            geom_point(alpha = 0.8, position = position_identity()) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
            scale_size_continuous(name = "% Expressing", range = c(1, 10)) +
            scale_color_manual(
                name = "Condition",
                values = c("tumor" = "#fc8d59", "normal" = "#91bfdb"),
                labels = c("tumor" = "Tumor", "normal" = "Normal")
            ) +
            facet_grid(category ~ ., scales = "free_y", space = "free_y") +
            labs(x = "Log2 Fold Change (Tumor/Normal)", y = "", title = paste(ct, "Macrophages")) +
            theme_bw() +
            theme(
                strip.background = element_rect(fill = "lightgrey"),
                strip.text = element_text(face = "bold"),
                panel.grid.minor = element_blank(),
                legend.position = "right",
                plot.background = element_rect(fill = "white")
            )

        filename <- paste0("results/109.paper/Fig2/myeloid_lipid_", tolower(ct), "_dotplot.png")
        ggsave(filename, p, width = 7, height = 7)
    }

    # Create a combined plot showing both cell types
    p_combined <- ggplot(plot_data_long, aes(x = logFC, y = gene, size = pct_expr, color = condition, shape = cell_type)) +
        geom_point(alpha = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        scale_size_continuous(name = "% Expressing", range = c(1, 10)) +
        scale_color_manual(
            name = "Condition",
            values = c("tumor" = "#fc8d59", "normal" = "#91bfdb")
        ) +
        scale_shape_manual(
            name = "Cell Type",
            values = c("TREM2_LAM" = 16, "Foam" = 17)
        ) +
        facet_grid(category ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Log2 Fold Change (Tumor/Normal)", y = "") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = "lightgrey"),
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "right",
            plot.background = element_rect(fill = "white")
        )
    ggsave("results/109.paper/Fig2/myeloid_lipid_combined_dotplot.png", p_combined, width = 8, height = 9)

    "results/109.paper/Fig2"
}

paper_myeloid_GSEA <- function(sc_mye) {
    # Create directory for results
    dir.create("results/109.paper/Fig2", showWarnings = FALSE)
    # Get gene sets related to lipid metabolism
    # First try to get all collections and filter for lipid-related terms
    all_gs <- getGeneSets(library = c("H", "C2", "C5"))

    # Filter gene sets related to lipid metabolism using keywords
    lipid_keywords <- c(
        "lipid", "triglyceride", "fatty acid", "cholesterol",
        "lipoprotein", "sterol", "steroid", "adipocyte", "fat"
    )

    # Find gene sets containing these keywords
    lipid_gs <- list()
    for (gs_name in names(all_gs)) {
        if (any(sapply(lipid_keywords, function(x) grepl(x, tolower(gs_name))))) {
            lipid_gs[[gs_name]] <- all_gs[[gs_name]]
        }
    }

    # Print how many gene sets we're using
    cat("Running GSEA with", length(lipid_gs), "lipid-related gene sets\n")

    # Run GSEA on myeloid cells using UCell (fastest method and handles dropouts well)
    sc_mye <- runEscape(sc_mye,
        gene.sets = lipid_gs,
        method = "UCell",
        min.size = 5,
        new.assay.name = "GSEA_lipid",
        BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores() - 1)
    )

    # Perform normalization to account for dropout effects
    sc_mye <- performNormalization(sc_mye,
        assay = "GSEA_lipid",
        gene.sets = lipid_gs,
        make.positive = TRUE
    )

    # Identify top differentially enriched pathways between tumor and normal
    # Need to set identities for proper comparison
    sc_mye$cell_type_group <- paste(sc_mye$cell_type_dtl, sc_mye$group, sep = "_")
    Idents(sc_mye) <- "cell_type_group"

    # Function to perform differential enrichment for a cell type
    diff_enrichment <- function(cell_type) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")

        if ((ident1 %in% Idents(sc_mye)) && (ident2 %in% Idents(sc_mye))) {
            results <- FindMarkers(
                sc_mye,
                assay = "GSEA_lipid_normalized",
                ident.1 = ident1,
                ident.2 = ident2,
                min.pct = 0,
                logfc.threshold = 0
            ) %>%
                rownames_to_column("pathway") %>%
                as_tibble() %>%
                mutate(cell_type = cell_type) %>%
                filter(p_val_adj < 0.05) %>%
                arrange(desc(abs(avg_log2FC)))

            return(results)
        } else {
            return(NULL)
        }
    }

    # Get differential enrichment for TREM2_LAM and Foam cells
    diff_trem2 <- diff_enrichment("TREM2_LAM")
    diff_foam <- diff_enrichment("Foam")

    # Combine results
    diff_combined <- bind_rows(diff_trem2, diff_foam)

    # Save differential enrichment results
    write.csv(diff_combined, "results/109.paper/Fig2/myeloid_lipid_differential_enrichment.csv", row.names = FALSE)
    # Create lists of pathways to filter
    {
        downregulated_pathways <- c(
            "GOCC-SPHERICAL-HIGH-DENSITY-LIPOPROTEIN-PARTICLE",
            "GOBP-CHOLESTEROL-CATABOLIC-PROCESS",
            "GOBP-HIGH-DENSITY-LIPOPROTEIN-PARTICLE-CLEARANCE",
            "GOBP-VERY-LOW-DENSITY-LIPOPROTEIN-PARTICLE-CLEARANCE",
            "GOBP-TRIGLYCERIDE-RICH-LIPOPROTEIN-PARTICLE-CLEARANCE",
            # Additional lipoprotein processing pathways
            "GOBP-REGULATION-OF-VERY-LOW-DENSITY-LIPOPROTEIN-PARTICLE-REMODELING",
            # Additional cholesterol homeostasis pathways
            "GOBP-CELLULAR-RESPONSE-TO-CHOLESTEROL",
            "GOBP-NEGATIVE-REGULATION-OF-CHOLESTEROL-EFFLUX"
            # Additional fatty acid metabolism pathways
        )

        uptake_storage_pathways <- c(
            "GOMF-LIPID-BINDING",
            "GOMF-TRIGLYCERIDE-BINDING",
            "GOCC-LIPID-DROPLET",
            # Additional pathways related to cholesterol metabolism
            "GOBP-POSITIVE-REGULATION-OF-CHOLESTEROL-BIOSYNTHETIC-PROCESS",
            # Additional pathways for fatty acid regulation
            "GOBP-NEGATIVE-REGULATION-OF-FATTY-ACID-OXIDATION"
        )

        fatty_acid_pathways <- c(
            "REACTOME-CHOLESTEROL-BIOSYNTHESIS",
            "GOBP-POSITIVE-REGULATION-OF-CHOLESTEROL-METABOLIC-PROCESS",
            "GOBP-CELLULAR-LIPID-METABOLIC-PROCESS",
            "GOBP-NEGATIVE-REGULATION-OF-FATTY-ACID-BETA-OXIDATION",
            "GOBP-LIPID-METABOLIC-PROCESS",
            "REACTOME-FATTY-ACID-METABOLISM",
            "GOBP-FATTY-ACID-METABOLIC-PROCESS",
            "GOMF-FATTY-ACID-TRANSMEMBRANE-TRANSPORTER-ACTIVITY",
            "GOMF-LONG-CHAIN-FATTY-ACID-TRANSPORTER-ACTIVITY",
            # Additional fatty acid pathways
            "GOMF-LONG-CHAIN-FATTY-ACID-BINDING"
        )

        # All pathways to filter - removing duplicates
        all_filter_pathways <- unique(c(downregulated_pathways, uptake_storage_pathways, fatty_acid_pathways))

        # Filter the differential enrichment results for specific pathways
        filtered_pathways <- diff_combined %>%
            filter(pathway %in% all_filter_pathways) %>%
            mutate(pathway_category = case_when(
                pathway %in% downregulated_pathways ~ "Downregulated cholesterol/lipoprotein pathways",
                pathway %in% uptake_storage_pathways ~ "Enhanced lipid uptake and storage",
                pathway %in% fatty_acid_pathways ~ "Altered fatty acid metabolism",
                TRUE ~ "Other"
            ))

        # Save the filtered results
        write.csv(filtered_pathways, "results/109.paper/Fig2/myeloid_lipid_filtered_pathways.csv", row.names = FALSE)
        # Create lollipop plots for each cell type
        # For each cell type, create a lollipop plot
        cell_types <- unique(filtered_pathways$cell_type)

        for (ct in cell_types) {
            # Filter data for this cell type
            ct_data <- filtered_pathways %>%
                filter(cell_type == ct) %>%
                # Shorten pathway names for better display
                mutate(
                    pathway_short = gsub("^(GOBP|GOMF|GOCC|HALLMARK|REACTOME|WP)-", "", pathway),
                    pathway_short = gsub("-", " ", pathway_short),
                    # Capitalize first letter of each word
                    pathway_short = tools::toTitleCase(tolower(pathway_short))
                ) %>%
                # Order by pathway category first, then by logFC
                arrange(pathway_category, avg_log2FC) %>%
                mutate(
                    # Create factor with proper ordering
                    pathway_short = factor(pathway_short, levels = unique(pathway_short))
                )

            # Create lollipop plot
            p <- ggplot(ct_data, aes(x = avg_log2FC, y = pathway_short, color = pathway_category)) +
                geom_segment(aes(x = 0, xend = avg_log2FC, y = pathway_short, yend = pathway_short),
                    color = "grey50"
                ) +
                geom_point(size = 4, alpha = 0.8) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
                # Use size to represent significance
                scale_size_continuous(range = c(2, 6)) +
                # Custom color palette for pathway categories
                scale_color_brewer(palette = "Set1") +
                labs(
                    x = "Log2 Fold Change (Tumor vs Normal)",
                    y = "",
                    title = paste0(ct, " Lipid Pathway Enrichment"),
                    color = "Pathway Category"
                ) +
                theme_bw() +
                theme(
                    panel.grid.minor = element_blank(),
                    panel.grid.major.y = element_blank(),
                    legend.position = "right",
                    axis.text.y = element_text(size = 9),
                    plot.background = element_rect(fill = "white")
                )

            # Save the plot
            ggsave(
                paste0(
                    "results/109.paper/Fig2/myeloid_lipid_lollipop_",
                    gsub(" ", "_", tolower(ct)), ".png"
                ),
                p,
                width = 10, height = 0.3 * nrow(ct_data) + 2
            )
        }
    }


    # Create heatmap showing these specific pathways
    # if(nrow(filtered_pathways) > 0) {
    #     # First organize data for heatmap
    #     filtered_pathways <- filtered_pathways %>%
    #         arrange(pathway_category, desc(avg_log2FC))

    #     # Make a plot focusing on these pathways
    #     p_filtered <- heatmapEnrichment(sc_mye,
    #                          group.by = "cell_type_dtl",
    #                          facet.by = "group",
    #                          assay = "GSEA_lipid",
    #                          scale = TRUE,
    #                          gene.sets = filtered_pathways$pathway,
    #                          cluster.rows = FALSE,
    #                          cluster.columns = TRUE)

    #     ggsave("results/109.paper/Fig2/myeloid_lipid_GSEA_filtered_heatmap.png",
    #            p_filtered, width = 10, height = 8)

    #     # Create ridgeplots for these specific pathways
    #     for (pathway in filtered_pathways$pathway) {
    #         p_ridge <- ridgeEnrichment(sc_mye,
    #                            assay = "GSEA_lipid",
    #                            gene.set = pathway,
    #                            group.by = "cell_type_dtl",
    #                            facet.by = "group",
    #                            add.rug = TRUE)

    #         ggsave(paste0("results/109.paper/Fig2/ridge_filtered_",
    #                       gsub("[^[:alnum:]]", "_", pathway), ".png"),
    #                p_ridge, width = 8, height = 6)
    #     }
    # }

    # Return the path to the results directory
    return("results/109.paper/Fig2")
}

paper_tcell_exhaustion <- function(sc_tcell) {
    ext_marker <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "ENTPD1")
    p <- DotPlot2(sc_tcell,
        features = ext_marker,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig3/tcell_exhaustion_dotplot.png", p, width = 10, height = 8)
    "results/109.paper/Fig3"
}

paper_tcell_fate_DEG <- function(sc_tcell) {
    # Set identity to cell_type_dtl for the comparison
    sc_tcell <- sc_tcell %>% filter(group == "tumor")
    Idents(sc_tcell) <- "cell_type_dtl"

    # Perform differential expression analysis between CD8_IFNG_gdT and CD4_Th17
    deg_results <- FindMarkers(
        sc_tcell,
        ident.1 = "CD8_Cyto_gdT",
        ident.2 = "CD4_Treg/Th17",
        test.use = "MAST",
        min.pct = 0.1,
        logfc.threshold = 0.25
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05) %>%
        filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        Add_Pct_Diff()

    # Save full results to CSV
    write.csv(deg_results, "results/109.paper/Fig3/CD8_Cyto_gdT_vs_CD4_TregTh17_MAST_DEG.csv", row.names = FALSE)
    # Create a pseudobulk object for more robust differential expression analysis
    sc_pseudo <- AggregateExpression(sc_tcell,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c("dataset", "cell_type_dtl")
    )

    Idents(sc_pseudo) <- "cell_type_dtl"

    # Define identities for comparison
    ident1 <- "CD8-Cyto-gdT"
    ident2 <- "CD4-Treg/Th17"

    # Run DESeq2 analysis on the pseudobulk data
    pseudo_deg <-
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
        # filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        filter(!str_starts(gene, "CYP")) %>%
        filter(!str_starts(gene, "MT-"))
    write_tsv(
        pseudo_deg %>%
            mutate(across(where(is.numeric), ~ ifelse(abs(.) < 0.001, signif(., 3), round(., 3)))),
        "results/109.paper/Fig3/CD8_Cyto_gdT_vs_CD4_TregTh17_pseudobulk_DEG.tsv"
    )

    # Save pseudobulk DEG results if successful

    # # Select top DEGs for visualization (top 10 up and down)
    # top_degs <- bind_rows(
    #     deg_results %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% head(10),
    #     deg_results %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% head(10)
    # )

    # # Create a volcano plot
    # p_volcano <- EnhancedVolcano(deg_results,
    #     lab = deg_results$gene,
    #     x = "avg_log2FC",
    #     y = "p_val_adj",
    #     title = "CD8_IFNG_gdT vs CD4_Th17",
    #     subtitle = "Differential Expression",
    #     pCutoff = 0.05,
    #     FCcutoff = 0.5,
    #     pointSize = 3.0,
    #     labSize = 4.0,
    #     colAlpha = 0.5,
    #     legendPosition = "right",
    #     drawConnectors = TRUE,
    #     widthConnectors = 0.5,
    #     colConnectors = "grey50"
    # ) +
    #     theme(plot.background = element_rect(fill = "white"))

    # ggsave("results/109.paper/Fig3/tcell_fate_volcano.png", p_volcano, width = 10, height = 8)

    # # Create a heatmap of top DEGs
    # top_genes <- c(
    #     top_degs %>% filter(avg_log2FC > 0) %>% pull(gene),
    #     top_degs %>% filter(avg_log2FC < 0) %>% pull(gene)
    # )

    # # Calculate scaled expression of top DEGs
    # sc_subset <- subset(sc_tcell, cell_type_dtl %in% c("CD8_IFNG_gdT", "CD4_Th17"))
    # exp_data <- AverageExpression(sc_subset,
    #     features = top_genes,
    #     group.by = "cell_type_dtl",
    #     assays = "RNA"
    # )$RNA

    # exp_data <- scale(exp_data)

    # # Set up annotation indicating which genes are upregulated in which cell type
    # gene_group <- ifelse(top_genes %in% (top_degs %>% filter(avg_log2FC > 0) %>% pull(gene)),
    #     "Higher in CD8_IFNG_gdT", "Higher in CD4_Th17"
    # )

    # # Create a heatmap
    # ha <- rowAnnotation(
    #     Group = gene_group,
    #     col = list(Group = c(
    #         "Higher in CD8_IFNG_gdT" = "#E41A1C",
    #         "Higher in CD4_Th17" = "#377EB8"
    #     ))
    # )

    # hm <- ComplexHeatmap::Heatmap(
    #     exp_data,
    #     name = "Scaled Expression",
    #     cluster_rows = TRUE,
    #     cluster_columns = FALSE,
    #     show_row_names = TRUE,
    #     show_column_names = TRUE,
    #     row_split = gene_group,
    #     right_annotation = ha
    # )

    # # Save the heatmap
    # png("results/109.paper/Fig3/tcell_fate_heatmap.png", width = 8, height = 10, units = "in", res = 300)
    # draw(hm)
    # dev.off()

    # # Generate a dot plot for selected genes
    # # Select biologically relevant genes for T cell fate and function
    # selected_genes <- list(
    #     "Effector function" = c("IFNG", "GZMB", "GZMA", "PRF1", "NKG7"),
    #     "T cell activation" = c("CD69", "IL2RA", "IL7R", "IL2RB"),
    #     "T cell lineage" = c("CD8A", "CD8B", "CD4", "TRDV1", "TRDV2"),
    #     "T helper" = c("IL17A", "IL22", "RORC", "TBX21", "GATA3")
    # )

    # p_dotplot <- DotPlot2(sc_tcell,
    #     features = selected_genes,
    #     group.by = "cell_type_dtl",
    #     idents = c("CD8_IFNG_gdT", "CD4_Th17"),
    #     show_grid = FALSE
    # ) +
    #     theme(
    #         plot.background = element_rect(fill = "white"),
    #         axis.text.x = element_text(angle = 45, hjust = 1)
    #     )

    # ggsave("results/109.paper/Fig3/tcell_fate_dotplot.png", p_dotplot, width = 10, height = 6)

    # Return the path to the results directory
    return("results/109.paper/Fig3")
}

paper_tcell_full_dotplot <- function(sc_tcell) {
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
        "Naive_T" = c("CCR7", "SELL", "TCF7", "LEF1", "SKAP1", "THEMIS", "PTPRC"),
        "Memory_T" = c("CD44", "IL7R", "PRKCQ", "STAT4"),

        # Tissue-Resident Memory
        "Trm" = c("ITGAE", "CD69", "CXCR6"),

        # Effector Functions
        "Effector" = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "TNF", "FASLG"),

        # Exhaustion Markers
        "Exhausted" = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "ENTPD1"),

        # T Helper Subtypes and Key Transcription Factors
        "Th1" = c("IFNG", "CXCR3", "TBX21"), # TBX21 (T-bet)
        "Th2" = c("IL4", "IL5", "IL13", "CCR4", "GATA3"),
        "Th17" = c("IL17F", "IL22", "CCR6", "RORC")
    )
    p <- DotPlot2(sc_tcell,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        # split.by.method = "color",
        # split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig3/tcell_full_dotplot.png", p, width = 8, height = 15)
}
paper_hormone_receptor_expression <- function(sc_final) {
    # Create directory for results
    dir.create("results/109.paper/HormoneReceptors", recursive = TRUE, showWarnings = FALSE)

    # Define the receptor genes
    receptors <- c("MC2R", "EPOR") # MC2R is the ACTH receptor, EPOR is the EPO receptor

    # Create violin plots to visualize expression across cell types
    p_violin <- VlnPlot2(sc_final,
        features = receptors,
        group.by = "cell_type_dtl",
        pt.size = 0,
        ncol = 1
    ) +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_violin.png", p_violin, width = 12, height = 10)

    # Create dot plots for better visualization of expression pattern
    p_dot <- DotPlot2(sc_final,
        features = receptors,
        group.by = "cell_type_dtl",
        split.by = "group"
    ) +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_dotplot.png", p_dot, width = 10, height = 8)

    # Calculate average expression by cell type
    avg_exp <- AverageExpression(sc_final,
        features = receptors,
        assays = "RNA",
        group.by = "cell_type_dtl"
    )$RNA %>%
        as.data.frame() %>%
        rownames_to_column("receptor") %>%
        pivot_longer(-receptor, names_to = "cell_type", values_to = "avg_expression")

    # Sort and identify top expressing cell types
    top_expressing <- avg_exp %>%
        group_by(receptor) %>%
        arrange(desc(avg_expression)) %>%
        slice_head(n = 5) %>%
        ungroup()

    # Save the results
    write.csv(avg_exp, "results/109.paper/HormoneReceptors/hormone_receptors_avg_expression.csv", row.names = FALSE)
    write.csv(top_expressing, "results/109.paper/HormoneReceptors/hormone_receptors_top_expressing.csv", row.names = FALSE)

    # Create bar plot of top expressing cell types
    p_bar <- ggplot(top_expressing, aes(x = reorder(cell_type, avg_expression), y = avg_expression, fill = receptor)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~receptor, scales = "free_y") +
        labs(
            title = "Top 5 Cell Types by Hormone Receptor Expression",
            x = "Cell Type",
            y = "Average Expression"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        )

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_top_expressing_barplot.png", p_bar, width = 12, height = 8)

    # Create a detailed feature plot to visualize expression on UMAP
    p_feature <- DimPlot2(sc_final,
        features = receptors,
        split.by = "group",
        reduction = "umap_integrated",
        ncol = 1,
        pt.size = 0.5
    ) +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_featureplot.png", p_feature, width = 10, height = 12)

    # Calculate the percentage of cells expressing each receptor by cell type using magic imputation
    percent_expr <- data.frame()
    
    # Create a temporary copy of the Seurat object for imputation
    message("Performing MAGIC imputation for hormone receptors...")
    sc_final_imp <- sc_final
    
    # Run MAGIC imputation for the receptor genes
    # Using the Rmagic package which interfaces with the Python MAGIC algorithm
    library(Rmagic)
    
    # Extract counts for imputation
    expression_matrix <- GetAssayData(sc_final, assay = "RNA", slot = "data")
    
    # Perform MAGIC imputation on the expression matrix (focusing on receptors to save computation)
    # This will help recover dropout events in sparse scRNA-seq data
    genes_to_impute <- c(receptors, sample(rownames(expression_matrix), 100)) # Include some random genes to improve imputation
    imputed_matrix <- magic(expression_matrix[genes_to_impute,], genes=receptors, 
                           knn=5, t="auto", npca=20, seed=42)
    
    # For each receptor, analyze both original and imputed expression
    for (receptor in receptors) {
        # Extract expression data from the normalized RNA assay (original)
        orig_expr <- GetAssayData(sc_final, assay = "RNA", slot = "data")[receptor, ]
        
        # Get imputed expression
        imputed_expr <- imputed_matrix[receptor, ]
        
        # Combine with cell type information
        expr_df <- data.frame(
            cell = names(orig_expr),
            orig_expression = orig_expr,
            imputed_expression = imputed_expr[names(orig_expr)],
            cell_type = sc_final$cell_type_dtl[names(orig_expr)],
            group = sc_final$group[names(orig_expr)]
        )
        
        # Calculate percentage using both original and imputed values
        pct_result <- expr_df %>%
            group_by(cell_type, group) %>%
            summarize(
                total_cells = n(),
                # Original expression metrics
                orig_expressing_cells = sum(orig_expression > 0),
                orig_percent_expressing = round(100 * sum(orig_expression > 0) / n(), 2),
                # Imputed expression metrics (with a slightly higher threshold to account for imputation bias)
                imputed_expressing_cells = sum(imputed_expression > 0.1),
                imputed_percent_expressing = round(100 * sum(imputed_expression > 0.1) / n(), 2),
                # Average expression
                avg_orig_expression = mean(orig_expression),
                avg_imputed_expression = mean(imputed_expression)
            ) %>%
            ungroup() %>%
            mutate(receptor = receptor)
        
        percent_expr <- bind_rows(percent_expr, pct_result) %>%
            arrange(receptor, desc(imputed_percent_expressing))
    }

    # Save the percentage results
    write.csv(percent_expr, "results/109.paper/HormoneReceptors/hormone_receptors_percent_expressing.csv", row.names = FALSE)

    # Create heatmap showing percentage of cells expressing each receptor by cell type
    p_heatmap <- ggplot(percent_expr, aes(x = cell_type, y = receptor, fill = percent_expressing)) +
        geom_tile() +
        facet_grid(. ~ group) +
        scale_fill_gradient(low = "white", high = "#E41A1C", name = "% Expressing") +
        labs(
            title = "Percentage of Cells Expressing Hormone Receptors",
            x = "Cell Type",
            y = "Receptor"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        )

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_percent_heatmap.png", p_heatmap, width = 14, height = 5)

    # Create a bar plot showing percentage of cells expressing each receptor
    p_bar <- ggplot(percent_expr, aes(x = reorder(cell_type, percent_expressing), y = percent_expressing, fill = group)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~receptor, scales = "free_y", ncol = 1) +
        labs(
            title = "Percentage of Cells Expressing Hormone Receptors by Cell Type",
            x = "Cell Type",
            y = "% Cells Expressing"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        )

    ggsave("results/109.paper/HormoneReceptors/hormone_receptors_percent_barplot.png", p_bar, width = 12, height = 10)

    # Get the top 5 cell types with highest expression percentage for each receptor and group
    top5_percent <- percent_expr %>%
        group_by(receptor, group) %>%
        arrange(desc(percent_expressing)) %>%
        slice_head(n = 5) %>%
        ungroup()

    # Save the top expressing results
    write.csv(top5_percent, "results/109.paper/HormoneReceptors/hormone_receptors_top5_percent.csv", row.names = FALSE)
    # Return path to the results directory
    return("results/109.paper/HormoneReceptors")
}
