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
        "CLC" = c("TH", "CHGA", "CHGB"),
        "Endo" = c("PECAM1", "EMCN"),
        "Fib" = c("COL1A1", "COL3A1", "THY1"),
        "PSC" = c("RGS5", "PDGFRB"),
        "Adipo" = c("ADIPOQ", "FABP4", "PPARG"),
        "Tcell" = c("CD3D", "CD3E", "TRBC1"),
        "Bcell" = c("CD19", "CD79A", "MS4A1"),
        "Myeloid" = c("ITGAM", "CD33", "ANPEP"),
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
    "results/109.paper/Fig1"
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
