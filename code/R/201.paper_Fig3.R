paper_myeloid_annotate <- function(sc_mye) {
    markers <- list(
        "TREM2_LAM" = c("CD40", "STAT1", "TREM2"), # TREM2+ Macrophages
        "DC" = c("PPARG", "CD83", "CD86"), # Foam Macrophages
        "APOE_LAM" = c("APOE", "ID1"),
        "Mast" = c("TPSAB1", "LYZ", "CCR2"), # Mast Cells
        "Eosinophil" = c("IL5RA", "EPO", "CXCR2"), # Eosinophils
        "Tissue-resident_Mac" = c("LYVE1", "MRC1", "CD163"), # Perivascular Macrophages
        "Erythrocyte" = c("GYPA", "HBB", "SLC4A1") # Erythrocytes
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
    ggsave("results/109.paper/Fig3/heatmap_myeloid.tiff", p2, width = 9, height = 5)
    p <- DimPlot2(sc_mye,
        reduction = "umap_mye",
        group.by = "cell_type_dtl",
        label = TRUE,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/109.paper/Fig3/dimplot_myeloid.tiff", p, width = 7, height = 5)
    # p <- DotPlot2(sc_mye,
    #     features = markers,
    #     group.by = "cell_type_dtl",
    #     split.by = "group",
    #     split.by.method = "color",
    #     split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
    #     show_grid = FALSE
    # ) +
    #     theme(plot.background = element_rect(fill = "white"))
    # ggsave("results/109.paper/Fig3/myeloid_annotation_dotplot.tiff", p, width = 7, height = 7)
    "results/109.paper/Fig3"
}
paper_myeloid_cell_counts <- function(sc_mye) {
    # Ensure directory exists
    dir.create("results/109.paper/Fig3", showWarnings = FALSE, recursive = TRUE)
    
    # Extract cell type and dataset information
    cell_data <- data.frame(
        cell_type = sc_mye$cell_type_dtl,
        dataset = sc_mye$dataset,
        group = sc_mye$group,
        stringsAsFactors = FALSE
    )
    
    # Count cells per cell type per dataset
    cell_counts <- cell_data %>%
        group_by(dataset, cell_type, group) %>%
        summarize(count = n(), .groups = "drop") %>%
        # Ensure all cell types are present for all datasets (fill with 0)
        complete(dataset, cell_type, fill = list(count = 0)) %>%
        # Join back with group information
        left_join(unique(cell_data[,c("dataset", "group")]), by = "dataset") %>%
        # Fix the group column issue by using group.y as the definitive group
        mutate(group = coalesce(group.y, group.x)) %>%
        select(-group.x, -group.y)  # Remove the redundant columns
    
    # Order cell types by their overall abundance
    cell_type_order <- cell_data %>%
        dplyr::count(cell_type) %>%
        arrange(desc(n)) %>%
        pull(cell_type)
    
    cell_counts$cell_type <- factor(cell_counts$cell_type, levels = cell_type_order)
    
    # Create the bar plot
    p <- ggplot(cell_counts, aes(x = dataset, y = count, fill = cell_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_grid(. ~ group, scales = "free_x", space = "free") +
        scale_fill_brewer(palette = "Set2", name = "Cell Type") +
        labs(
            title = "Number of Cells per Myeloid Subpopulation by dataset",
            x = "dataset",
            y = "Cell Count"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            plot.background = element_rect(fill = "white"),
            strip.background = element_rect(fill = "lightgrey"),
            strip.text = element_text(face = "bold")
        )
    
    # Save the plot
    ggsave("results/109.paper/Fig3/myeloid_cell_counts_per_sample.tiff", p, width = 12, height = 6, dpi = 300)
    # Save cell counts to TSV file - grouped by group and cell_type
    cell_counts_summary <- cell_counts %>%
        group_by(group, cell_type) %>%
        summarize(total_count = sum(count), .groups = "drop") %>%
        arrange(group, desc(total_count))
    
    write.table(cell_counts_summary, 
                file = "results/109.paper/Fig3/myeloid_cell_counts_summary.tsv", 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)

    # Create a stacked percentage plot to show relative composition
    cell_pct <- cell_counts %>%
        group_by(dataset) %>%
        mutate(percentage = count / sum(count) * 100) %>%
        ungroup()
    
    p2 <- ggplot(cell_pct, aes(x = dataset, y = percentage, fill = cell_type)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_grid(. ~ group, scales = "free_x", space = "free") +
        scale_fill_brewer(palette = "Set2", name = "Cell Type") +
        labs(
            title = "Percentage of Cells per Myeloid Subpopulation by Sample",
            x = "Sample",
            y = "Percentage (%)"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            plot.background = element_rect(fill = "white"),
            strip.background = element_rect(fill = "lightgrey"),
            strip.text = element_text(face = "bold")
        )
    
    # Save the percentage plot
    ggsave("results/109.paper/Fig3/myeloid_cell_percentage_per_sample.tiff", p2, width = 12, height = 6, dpi = 300)
    
    # Return path to the results
    return("results/109.paper/Fig3")
}
paper_myeloid_lipid_DEG <- function(sc_mye) {
    # markers 按功能分组
    apoe_markers <- list(
        "Established lipid & immune regulators" = c("ABCA8", "ABCA1", "ABCG1", "GULP1", "RORA", "FST"),
        "Additional novel/underappreciated" = c("GGT5", "PRKG1", "KALRN", "AC069544.1", "MEIS1-AS2", "LINC01798")
    )
    trem2_markers <- list(
        "Complement activation" = c("C3"),
        "Cell stress & mTOR" = c("DDIT4"),
        "Immune checkpoints/regulatory" = c("ENTPD1", "SIGLEC8", "PALD1"),
        "Downregulated migratory/exhausted" = c("NR4A1", "NR4A2", "NR4A3", "RGS1", "CD36", "PPARG"),
        "Downregulated pro-inflammatory" = c("NLRP3", "FOSB"),
        "Downregulated kinases" = c("FYN", "BACH2", "MMP19")
    )
    # 不再绘制原始 dotplot

    # --- Step 2: Prepare the gene list ---
    # Flatten the list of markers into a single vector
    # 合并所有 marker
    all_markers <- unique(unlist(c(apoe_markers, trem2_markers)))

    # gene 分组信息
    gene_categories <- data.frame(
        gene = unlist(apoe_markers),
        category = rep(names(apoe_markers), times = sapply(apoe_markers, length))
    )
    gene_categories <- rbind(
        gene_categories,
        data.frame(
            gene = unlist(trem2_markers),
            category = rep(names(trem2_markers), times = sapply(trem2_markers, length))
        )
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
    deg_foam <- get_cell_type_degs("APOE_LAM")

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

    # 添加 category 信息
    plot_data <- merge(plot_data, gene_categories, by = "gene")

    # 转为长表
    plot_data_long <- pivot_longer(
        plot_data,
        cols = c(pct_expr_tumor, pct_expr_normal),
        names_to = "condition",
        values_to = "pct_expr"
    )
    plot_data_long$condition <- gsub("pct_expr_", "", plot_data_long$condition)
    plot_data_long <- plot_data_long %>%
        group_by(cell_type, category) %>%
        mutate(gene = factor(gene, levels = unique(gene[order(logFC)]))) %>%
        ungroup() %>%
        arrange(gene, category, desc(pct_expr))

    # Create plots for each cell type
    for (ct in c("TREM2_LAM", "APOE_LAM")) {
        if (ct == "APOE_LAM") {
            marker_set <- unlist(apoe_markers)
            marker_levels <- apoe_markers
        } else if (ct == "TREM2_LAM") {
            marker_set <- unlist(trem2_markers)
            marker_levels <- trem2_markers
        }
        cell_data <- plot_data_long %>% filter(cell_type == ct, gene %in% marker_set)
        # 仅对 trem2_LAM，按每个 category 内 logFC 降序排序 gene
        if (ct == "TREM2_LAM") {
            cell_data <- cell_data %>%
                group_by(category) %>%
                mutate(gene = factor(gene, levels = unique(gene[order(logFC)]))) %>%
                ungroup()
        }
        # 只保留 APOE_LAM 的肿瘤样本
        if (ct == "APOE_LAM") {
            cell_data <- cell_data %>% filter(condition == "tumor")
        }

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

        filename <- paste0("results/109.paper/Fig3/myeloid_lipid_", tolower(ct), "_dotplot.tiff")
        ggsave(filename, p, width = 7, height = 7)
    }

    # 可选：如需合并图，可仿照上述方式，仅保留新 marker

    "results/109.paper/Fig3"
}

paper_TREM2_LAM_violin <- function(sc_mye) {
    # 筛选TREM2_LAM细胞
    cells <- WhichCells(sc_mye, expression = cell_type_dtl == "TREM2_LAM")
    sc_sub <- subset(sc_mye, cells = cells)
    # 绘制小提琴图
    p <- VlnPlot2(
        sc_sub,
        features = c("CD36", "PPARG"),
        split.by = "group",
        nrow = 1
    ) +
        theme_minimal(base_size = 14) +
        labs(
            title = "CD36, PPARG expression in TREM2_LAM",
            x = "Gene",
            y = "Expression"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "top",
            plot.background = element_rect(fill = "white")
        )
    ggsave("results/109.paper/Fig3/TREM2_LAM_violin.tiff", p, width = 10, height = 10)

    # Add comparison between APOE_LAM and other myeloid cells
    # Create a new column to identify APOE_LAM vs other myeloid cells
    sc_mye$comparison_group <- ifelse(sc_mye$cell_type_dtl == "APOE_LAM",
        "APOE_LAM",
        "Other Myeloid"
    )

    # Plot violin plots for ABCA8 and RORA
    p2 <- VlnPlot2(
        sc_mye,
        features = c("ABCA8", "RORA"),
        group.by = "comparison_group",
        ncol = 1
    ) +
        theme_minimal(base_size = 14) +
        labs(
            title = "ABCA8, RORA expression in APOE_LAM vs other myeloid cells",
            x = "Cell Type",
            y = "Expression"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "top",
            plot.background = element_rect(fill = "white")
        )
    ggsave("results/109.paper/Fig3/APOE_LAM_vs_others_violin.tiff", p2, width = 10, height = 10)
    "results/109.paper/Fig3/TREM2_LAM_violin.tiff"
}
paper_myeloid_GSEA <- function(sc_mye) {
    # Create directory for results
    dir.create("results/109.paper/Fig3", showWarnings = FALSE)
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
    write.csv(diff_combined, "results/109.paper/Fig3/myeloid_lipid_differential_enrichment.csv", row.names = FALSE)
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
        write.csv(filtered_pathways, "results/109.paper/Fig3/myeloid_lipid_filtered_pathways.csv", row.names = FALSE)
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
                    "results/109.paper/Fig3/myeloid_lipid_lollipop_",
                    gsub(" ", "_", tolower(ct)), ".tiff"
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

    #     ggsave("results/109.paper/Fig3/myeloid_lipid_GSEA_filtered_heatmap.tiff",
    #            p_filtered, width = 10, height = 8)

    #     # Create ridgeplots for these specific pathways
    #     for (pathway in filtered_pathways$pathway) {
    #         p_ridge <- ridgeEnrichment(sc_mye,
    #                            assay = "GSEA_lipid",
    #                            gene.set = pathway,
    #                            group.by = "cell_type_dtl",
    #                            facet.by = "group",
    #                            add.rug = TRUE)

    #         ggsave(paste0("results/109.paper/Fig3/ridge_filtered_",
    #                       gsub("[^[:alnum:]]", "_", pathway), ".tiff"),
    #                p_ridge, width = 8, height = 6)
    #     }
    # }

    # Return the path to the results directory
    return("results/109.paper/Fig3")
}

paper_macrophage_proportion_change <- function(sc_mye) {
    # 确保目录存在
    dir.create("results/109.paper/Fig3", showWarnings = FALSE, recursive = TRUE)

    # 筛选所有巨噬细胞类型（APOE_LAM、TREM2_LAM和Tissue-resident_Mac）
    # 按照注释，这些细胞类型都是巨噬细胞亚群
    mac_cells <- WhichCells(sc_mye, expression = cell_type_dtl %in% c("APOE_LAM", "TREM2_LAM", "Tissue-resident_Mac"))
    sc_mac <- subset(sc_mye, cells = mac_cells)

    # 计算每个组(正常/肿瘤)中各巨噬细胞亚群的数量和比例
    mac_counts <- table(sc_mac$cell_type_dtl, sc_mac$group)
    mac_props <- prop.table(mac_counts, margin = 2) * 100

    # 将数据转换为长格式，便于绘图
    mac_data <- as.data.frame(mac_props)
    colnames(mac_data) <- c("Macrophage_Type", "Condition", "Percentage")

    # 确保Macrophage_Type列是因子，设置显示顺序
    mac_data$Macrophage_Type <- factor(mac_data$Macrophage_Type,
        levels = c("TREM2_LAM", "APOE_LAM", "Tissue-resident_Mac")
    )

    # 确保正常和肿瘤条件的顺序正确
    mac_data$Condition <- factor(mac_data$Condition, levels = c("normal", "tumor"))

    # 绘制堆叠条形图
    p1 <- ggplot(mac_data, aes(x = Condition, y = Percentage, fill = Macrophage_Type)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_brewer(palette = "Set2") +
        labs(
            title = "Macrophage Subtype Distribution",
            x = "Condition",
            y = "Percentage (%)",
            fill = "Macrophage Type"
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            legend.position = "right",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 16)
        )

    # 保存堆叠条形图
    ggsave("results/109.paper/Fig3/macrophage_proportion_stacked.tiff", p1, width = 7, height = 5)

    # 绘制分组条形图，每种巨噬细胞类型单独显示
    p2 <- ggplot(mac_data, aes(x = Macrophage_Type, y = Percentage, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("normal" = "#4DBBD5FF", "tumor" = "#E64B35FF")) +
        labs(
            title = "Macrophage Subtype Changes from Normal to Tumor",
            x = "Macrophage Type",
            y = "Percentage (%)",
            fill = "Condition"
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            legend.position = "right",
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 16)
        )

    # 保存分组条形图
    ggsave("results/109.paper/Fig3/macrophage_proportion_grouped.tiff", p2, width = 7, height = 5)

    # 计算每种巨噬细胞在肿瘤与正常组之间的变化倍数
    fold_change <- mac_props[, "tumor"] / mac_props[, "normal"]
    fc_data <- data.frame(
        Macrophage_Type = names(fold_change),
        FoldChange = as.numeric(fold_change)
    )
    fc_data$Macrophage_Type <- factor(fc_data$Macrophage_Type,
        levels = c("TREM2_LAM", "APOE_LAM", "Tissue-resident_Mac")
    )

    # 绘制倍数变化条形图
    p3 <- ggplot(fc_data, aes(x = Macrophage_Type, y = FoldChange, fill = Macrophage_Type)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +
        scale_fill_brewer(palette = "Set2") +
        labs(
            title = "Fold Change in Macrophage Proportions (Tumor/Normal)",
            x = "Macrophage Type",
            y = "Fold Change",
            fill = "Macrophage Type"
        ) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            legend.position = "right",
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 16)
        )

    # 保存倍数变化条形图
    ggsave("results/109.paper/Fig3/macrophage_proportion_foldchange.tiff", p3, width = 7, height = 5)

    # 统计检验：比较正常和肿瘤样本中巨噬细胞亚型比例的差异
    # 使用卡方检验比较巨噬细胞亚型在正常和肿瘤之间的分布差异
    chisq_result <- chisq.test(mac_counts)

    # 保存检验结果
    sink("results/109.paper/Fig3/macrophage_proportion_stats.txt")
    cat("Chi-square test for macrophage subtype distribution between normal and tumor:\n")
    print(chisq_result)

    # 添加每种巨噬细胞类型的百分比信息
    cat("\nPercentage of macrophage subtypes:\n")
    print(mac_props)

    # 添加倍数变化信息
    cat("\nFold change (tumor/normal):\n")
    print(fold_change)
    sink()

    "results/109.paper/Fig3"
}
