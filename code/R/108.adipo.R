subcluster_stromal_cells <- function(sc_final) {
    # Create directory for results
    dir.create("results/110.adipo/subclusters", showWarnings = FALSE, recursive = TRUE)

    # Define marker genes for different stromal cell types
    marker_genes <- list(
        "Fibroblast" = c("COL1A1", "COL1A2", "MGP", "PDGFRA"),
        "MSC/Progenitor" = c("ID1", "NR2F2", "GLI1", "NGFR"),
        "MSC_2" = c("ENG", "CD14", "VCAM1", "MCAM"),
        "Pericyte" = c("PDGFRB", "CSPG4", "RGS5", "ACTA2", "MCAM"),
        "CAR-like" = c("CXCL12", "LEPR", "FOXC1", "ADIPOQ"),
        "Negative Markers" = c("PECAM1", "CD34", "PTPRC", "EPCAM")
    )

    # Combine all markers into a single vector
    all_markers <- unlist(marker_genes)

    # Subset to only include Fibroblasts and PSCs
    combined_stromal_cells <- sc_final %>%
        filter(cell_type_dtl %in% c("Fib", "PSC"))

    # Run clustering on combined data
    combined_stromal_cells <- combined_stromal_cells %>%
        RunUMAP(dims = 1:30, reduction = "scvi", reduction.name = "combined_stromal_umap") %>%
        FindNeighbors(dims = 1:30, reduction = "scvi", graph.name = "combined_snn") %>%
        FindClusters(
            resolution = 0.5, algorithm = 4, method = "igraph",
            graph.name = "combined_snn", cluster.name = "combined_stromal_clusters"
        )

    # Create UMAP visualization colored by cluster
    p_combined_clusters <- DimPlot2(combined_stromal_cells,
        reduction = "combined_stromal_umap",
        group.by = "combined_stromal_clusters",
        label = TRUE
    ) +
        ggtitle("Combined tumor and normal stromal subclusters") +
        theme(plot.background = element_rect(fill = "white"))

    # Create UMAP visualization colored by group (tumor vs normal)
    p_combined_group <- DimPlot2(combined_stromal_cells,
        reduction = "combined_stromal_umap",
        group.by = "group",
        cols = c("normal" = "#4DBBD5FF", "tumor" = "#E64B35FF")
    ) +
        ggtitle("Combined stromal subclusters by sample type") +
        theme(plot.background = element_rect(fill = "white"))

    # Save plots
    ggsave("results/110.adipo/subclusters/combined_umap_clusters.png", p_combined_clusters, width = 8, height = 7)
    ggsave("results/110.adipo/subclusters/combined_umap_group.png", p_combined_group, width = 8, height = 7)

    # Create dotplot of all markers across clusters
    p_combined_dotplot <- DotPlot2(combined_stromal_cells,
        features = all_markers,
        group.by = "combined_stromal_clusters"
    ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/110.adipo/subclusters/combined_dotplot.png", p_combined_dotplot, width = 12, height = 6)

    # Generate heatmap of marker expression by subtype
    toplot_combined <- CalcStats(combined_stromal_cells,
        features = all_markers,
        group.by = "combined_stromal_clusters",
        method = "zscore",
        order = "value"
    )

    gene_groups_combined <- rep(names(marker_genes), lengths(marker_genes)) %>%
        setNames(marker_genes %>% unlist())
    gene_groups_combined <- gene_groups_combined[rownames(toplot_combined)]

    p_combined_heatmap <- Heatmap(t(toplot_combined),
        lab_fill = "zscore",
        facet_col = gene_groups_combined
    ) +
        ggtitle("Combined stromal clusters marker expression") +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/110.adipo/subclusters/combined_marker_heatmap_zscores.png",
        p_combined_heatmap,
        width = 14, height = 7
    )

    write_tsv(
        toplot_combined %>% as.data.frame() %>% rownames_to_column("gene"),
        "results/110.adipo/subclusters/combined_marker_heatmap_zscores.tsv"
    )

    return(combined_stromal_cells)
}
sub_annotation_adipo <- function(sc_adipo_clust, sc_final) {
    # Create directory for results
    dir.create("results/110.adipo/subclusters", showWarnings = FALSE, recursive = TRUE)

    # Define new annotations for the combined stromal clusters
    # Define annotations for the combined stromal clusters
    annotations <- c(
        "1" = "MSC",
        "2" = "CAR-like cells",
        "3" = "MSC",
        "4" = "MSC",
        "5" = "MSC",
        "6" = "PSC"
    )

    # Add the annotations to the subclustered object
    sc_adipo_clust$cell_type_dtl_new <- unname(annotations[as.character(sc_adipo_clust$combined_stromal_clusters)])

    # Create a mapping from cell IDs to new annotations
    cell_annotations <- data.frame(
        cell_id = Cells(sc_adipo_clust),
        cell_type_dtl = sc_adipo_clust$cell_type_dtl_new,
        stringsAsFactors = FALSE
    )

    # Update the cell_type_dtl in sc_final for cells in the subcluster
    sc_final$cell_type_dtl <- as.character(sc_final$cell_type_dtl)
    matching_cells <- intersect(Cells(sc_final), cell_annotations$cell_id)
    cell_mapping <- setNames(cell_annotations$cell_type_dtl, cell_annotations$cell_id)
    sc_final$cell_type_dtl[matching_cells] <- cell_mapping[matching_cells]

    sc_adipo <- sc_final %>%
        filter(cell_type_dtl %in% c("Fibroblasts", "MSC", "CAR-like cells", "PSC", "Adipo"))
    # do UMAP
    sc_adipo <- RunUMAP(sc_adipo, dims = 1:30, reduction = "scvi", reduction.name = "adipo_umap")
    # Create UMAP visualization with new annotations
    p_annotated_clusters <- DimPlot2(sc_adipo,
        reduction = "adipo_umap",
        group.by = "cell_type_dtl",
        label = TRUE
    ) +
        ggtitle("Annotated Adipocyte and Stromal cells") +
        theme(plot.background = element_rect(fill = "white"))

    # Save plot
    ggsave("results/110.adipo/annotated_umap_clusters.png", p_annotated_clusters, width = 8, height = 7)

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

    p_annotated_dotplot <- DotPlot2(sc_adipo,
        features = all_markers,
        group.by = "cell_type_dtl"
    ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/110.adipo/annotated_dotplot.png", p_annotated_dotplot, width = 12, height = 6)
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

    ggsave("results/110.adipo/heatmap_zscores.png",
        p_heatmap,
        width = 14, height = 7
    )

    write_tsv(
        toplot %>% as.data.frame() %>% rownames_to_column("gene"),
        "results/110.adipo/heatmap_zscores.tsv"
    )

    sc_adipo
}

test_adipo_composition <- function(sc_adipo) {
    # Create directory
    dir.create("results/110.adipo/distribution", showWarnings = FALSE, recursive = TRUE)

    # Prepare data for sccomp
    sc_adipo <- sc_adipo %>%
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
    composition_test <- sc_adipo %>%
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
    ggsave("results/110.adipo/distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save fold changes and test results
    fold_changes <- composition_test %>%
        sccomp_proportional_fold_change(
            formula_composition = ~group,
            from = "normal",
            to = "tumor"
        ) %>%
        select(cell_group, statement)

    write_tsv(fold_changes, "results/110.adipo/distribution/fold_changes.tsv")

    write_tsv(
        composition_test %>%
            select(-count_data),
        "results/110.adipo/distribution/test_results.tsv"
    )

    composition_test
}

plot_adipo_distribution <- function(sc_adipo, m_composition_test) {
    dir.create("results/110.adipo/distribution", showWarnings = FALSE, recursive = TRUE)

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

    tidyplots::save_plot(p, "results/110.adipo/distribution/adipo_distribution.png")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/110.adipo/distribution/adipo_distribution.tsv"
    )

    return("results/110.adipo/distribution")
}
plot_stromal_markers <- function(sc_adipo, sc_final) {
    # Create directory for results
    dir.create("results/110.adipo/stromal_markers", showWarnings = FALSE, recursive = TRUE)

    # Update the cell_type_dtl in sc_final from sc_adipo
    cell_mapping <- setNames(sc_adipo$cell_type_dtl, colnames(sc_adipo))
    matching_cells <- intersect(Cells(sc_final), Cells(sc_adipo))
    sc_final$cell_type_dtl[matching_cells] <- cell_mapping[matching_cells]

    # Define marker genes based on the categories provided
    marker_genes <- list(
        "Fibroblast Marker" = c("MCAM", "ENG", "POSTN", "COL1A1", "TNC", "COL1A2", "THY1", "COL3A1"),
        "MSC Marker" = c("STRIP1", "KLF5", "MTAP", "MFAP5", "EFNA5", "NT5E", "NES"),
        "Fibro + MSC Marker" = c(
            "ENG", "THY1", "ITGAM", "CD14", "CD19", "CD34",
            "PTPRC", "CD79A", "HLA-DRA", "ACTA2"
        )
    )

    # Combine all markers into a single vector
    all_markers <- unlist(marker_genes)

    # Create gene groups mapping for heatmap faceting
    gene_groups <- rep(names(marker_genes), lengths(marker_genes)) %>%
        setNames(marker_genes %>% unlist())

    # Define datasets to analyze
    datasets <- list(
        stromal = list(
            data = sc_final %>% filter(cell_type_dtl %in% c("Fib", "PSC", "MSC", "CAR-like cells", "Adipo")),
            output_prefix = "stromal",
            title = "Stromal cell"
        ),
        all = list(
            data = sc_final,
            output_prefix = "all_cell_types",
            title = "All cell types"
        )
    )

    # Create plots for each dataset
    for (dataset_name in names(datasets)) {
        dataset <- datasets[[dataset_name]]
        sc_data <- dataset$data
        prefix <- dataset$output_prefix
        title_prefix <- dataset$title
        
        # Basic dotplot
        p_dotplot <- DotPlot2(sc_data,
            features = all_markers,
            group.by = "cell_type_dtl"
        ) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            theme(plot.background = element_rect(fill = "white"))

        ggsave(paste0("results/110.adipo/stromal_markers/", prefix, "_markers_dotplot.png"), 
               p_dotplot, width = ifelse(dataset_name == "all", 14, 12), height = ifelse(dataset_name == "all", 8, 6))

        # Heatmap
        # Z-score heatmap
        toplot <- CalcStats(sc_data,
            features = all_markers,
            group.by = "cell_type_dtl",
            method = "zscore",
            order = "value"
        )

        # Filter gene_groups to only include genes present in the toplot data
        filtered_gene_groups <- gene_groups[rownames(toplot)]

        p_heatmap <- Heatmap(t(toplot),
            lab_fill = "zscore",
            facet_col = filtered_gene_groups
        ) +
            ggtitle(paste0(title_prefix, " marker expression (z-score)")) +
            theme(plot.background = element_rect(fill = "white"))

        ggsave(paste0("results/110.adipo/stromal_markers/", prefix, "_marker_heatmap_zscores.png"),
            p_heatmap,
            width = ifelse(dataset_name == "all", 16, 14), height = ifelse(dataset_name == "all", 10, 7)
        )
        
        # Raw mean expression heatmap
        toplot_mean <- CalcStats(sc_data,
            features = all_markers,
            group.by = "cell_type_dtl",
            method = "mean",
            order = "value"
        )

        p_heatmap_mean <- Heatmap(t(toplot_mean),
            lab_fill = "expression",
            facet_col = filtered_gene_groups
        ) +
            ggtitle(paste0(title_prefix, " marker expression (mean)")) +
            theme(plot.background = element_rect(fill = "white"))

        ggsave(paste0("results/110.adipo/stromal_markers/", prefix, "_marker_heatmap_mean.png"),
            p_heatmap_mean,
            width = ifelse(dataset_name == "all", 16, 14), height = ifelse(dataset_name == "all", 10, 7)
        )

        # Split by group dotplot
        p_split_dotplot <- DotPlot2(sc_data,
            features = all_markers,
            group.by = "cell_type_dtl",
            split.by = "group",
            split.by.method = "color",
            split.by.colors = c("#4DBBD5FF", "#E64B35FF")
        ) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            theme(plot.background = element_rect(fill = "white"))

        ggsave(paste0("results/110.adipo/stromal_markers/", prefix, "_markers_dotplot_by_group.png"), 
               p_split_dotplot, width = ifelse(dataset_name == "all", 16, 14), height = ifelse(dataset_name == "all", 10, 7))
    }

    # Calculate and save expression data for the stromal markers
    expr_data <- AverageExpression(datasets$stromal$data,
        features = all_markers,
        group.by = "cell_type_dtl",
        slot = "data"
    )

    write.csv(expr_data$RNA, "results/110.adipo/stromal_markers/stromal_marker_expression.csv")

    # Return the path to results
    return("results/110.adipo/stromal_markers")
}
DEG_Car_vs_MSC <- function(sc_adipo) {
    # Create directory for results
    dir.create("results/110.adipo/DEG_Car_vs_MSC", showWarnings = FALSE, recursive = TRUE)
    
    # Create a new column for cell type and condition combination
    sc_adipo$cell_type_group <- paste0(sc_adipo$cell_type_dtl, "_", sc_adipo$group)
    
    # Filter for only the cell types we're interested in
    sc_filtered <- sc_adipo[, sc_adipo$cell_type_group %in% c("CAR-like cells_tumor", "MSC_normal")]
    
    # Set the identity to the combined cell type and condition
    Idents(sc_filtered) <- "cell_type_group"
    
    # Method 1: MAST
    message("Running MAST for differential expression between CAR-like cells (tumor) vs MSC (normal)...")
    deg_mast <- FindMarkers(
        sc_filtered,
        ident.1 = "CAR-like cells_tumor",
        ident.2 = "MSC_normal",
        test.use = "MAST",
        min.pct = 0.1,
        logfc.threshold = 0
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
        arrange(desc(abs(avg_log2FC)))
    
    # Save MAST results
    write_tsv(
        deg_mast,
        "results/110.adipo/DEG_Car_vs_MSC/DEG_CAR_tumor_vs_MSC_normal_MAST.tsv"
    )
    
    # Method 2: Pseudobulk with DESeq2
    message("Running pseudobulk analysis for differential expression between CAR-like cells (tumor) vs MSC (normal)...")
    
    # Create pseudobulk data using AggregateExpression
    sc_pseudo <- AggregateExpression(
        sc_filtered,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c("dataset", "cell_type_group")
    )
    
    # Set identity for comparison
    Idents(sc_pseudo) <- "cell_type_group"
    
    # Run DESeq2 analysis
    deg_pseudobulk <- FindMarkers(
        sc_pseudo,
        ident.1 = "CAR-like cells-tumor",
        ident.2 = "MSC-normal",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
        arrange(desc(abs(avg_log2FC)))
    
    # Save pseudobulk results
    write_tsv(
        deg_pseudobulk,
        "results/110.adipo/DEG_Car_vs_MSC/DEG_CAR_tumor_vs_MSC_normal_pseudobulk_DESeq2.tsv"
    )
    
    
    # Return a list with both results
    return(list(
        DEG_MAST = deg_mast,
        DEG_pseudobulk = deg_pseudobulk
    ))
}
save_adipo_to_h5ad <- function(sc_adipo) {
    # Create directory if it doesn't exist
    dir.create("data/110.adipo", showWarnings = FALSE, recursive = TRUE)
    
    # Path for saving the h5ad file
    file_path <- "data/110.adipo"
    file_name <- "sc_adipo.h5ad"
    
    # Save the Seurat object as AnnData h5ad file
    scCustomize::as.anndata(
        x = sc_adipo, 
        file_path = file_path, 
        file_name = file_name
    )
    
    return(file.path(file_path, file_name))
}
plot_fib_volcano <- function(sc_adipo) {
    # Create directory if it doesn't exist
    dir.create("results/110.adipo/DEG", showWarnings = FALSE, recursive = TRUE)

    # Define genes to label
    genes_to_label <- c(
        # Upregulated
        "PPARG", "LPL", "CXCL12", "VEGFC", "VCAN", "ADAMTS17", "ADAMTS2",
        # Downregulated
        "LSAMP", "IGFBP6"
    )

    # Path to the pre-calculated DEG results
    deg_file_path <- "/workspaces/AML/results/103.DEG/tumor_vs_normal/DEG_Fib_DESeq2.tsv"

    # Check if the DEG file exists
    if (!file.exists(deg_file_path)) {
        warning(paste("DEG file not found:", deg_file_path))
        message("Skipping Volcano plot generation for Fib.")
        return(NULL)
    }

    # Read the DEG results
    deg_results <- read_tsv(deg_file_path)

    # Ensure required columns are present
    if (!all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(deg_results))) {
        warning(paste("Required columns (gene, avg_log2FC, p_val_adj) not found in", deg_file_path))
        message("Skipping Volcano plot generation for Fib.")
        return(NULL)
    }

    p_volcano_fib <- EnhancedVolcano(
        deg_results,
        lab = deg_results$gene,
        x = "avg_log2FC",
        y = "p_val_adj",
        selectLab = genes_to_label, # Label only the specified genes
        title = "Volcano Plot: Fib (Tumor vs Normal)",
        subtitle = "Differential Expression (DESeq2)",
        pCutoff = 0.05, # Adjust p-value cutoff as needed
        FCcutoff = 1, # Adjust logFC cutoff as needed
        pointSize = ifelse(deg_results$gene %in% genes_to_label, 4.0, 2.0), # Make labeled points bigger
        labSize = 4.0,
        colAlpha = 0.5,
        legendPosition = "none", # Remove legend
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        colConnectors = "grey50",
        boxedLabels = TRUE, # Make labels boxed
        # Define custom colors: non-significant, significant FC, significant p, significant both
        # Make non-labeled points grey, labeled points orange
        colCustom = data.frame(
            key = deg_results$gene,
            color = ifelse(deg_results$gene %in% genes_to_label, "orange", "grey30")
        ) %>% tibble::deframe()
    ) + theme(plot.background = element_rect(fill = "white")) # Ensure white background

    # Save the plot
    ggsave("results/110.adipo/DEG/fib_volcano_plot.png", p_volcano_fib, width = 8, height = 7)

    message("Fib volcano plot saved to results/110.adipo/DEG/fib_volcano_plot.png")

    # Return the path to the saved plot
    return("results/110.adipo/DEG/fib_volcano_plot.png")
}


cluster_adipo <- function(sc_final) {
    dir.create("results/110.adipo", showWarnings = FALSE, recursive = TRUE)
    # Filter for plasticity cells only (Fib, PSC, adipo)
    sc_adipo <- sc_final %>%
        filter(cell_type %in% c("Adipo", "Fib", "PSC"))

    # Run UMAP
    sc_adipo <- RunUMAP(sc_adipo, dims = 1:30, reduction = "scvi", verbose = FALSE)

    # Plot UMAP
    p <- DimPlot2(sc_adipo, group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/adipo_umap.png", p, width = 7, height = 6)

    return(sc_adipo)
}
run_adipo_trajectory <- function(sc_adipo) {
    dir.create("results/110.adipo/trajectory", showWarnings = FALSE, recursive = TRUE)
    tar_load(sc_adipo)
    reticulate::use_virtualenv("./scveloenv")
    library(reticulate)
    py_config()
    library(SeuratExtend)
    library(tidyseurat)
    # filter for plasticity cells only (Fib, PSC, adipo)

    # 1. scVelo Analysis
    # Convert Seurat object to AnnData for scVelo
    dir.create("data/110.adipo/trajectory", showWarnings = FALSE, recursive = TRUE)
    Seu2Adata(sc_adipo, save.adata = "data/110.adipo/trajectory/adipo.h5ad", conda_env = "base")
    adata_path <- "data/110.adipo/trajectory/adipo.h5ad"
    seu <- scVelo.SeuratToAnndata(
        sc_adipo,
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
        basis = "adipo_umap_cell_embeddings",
        save = "results/110.adipo/trajectory/velocity_umap.png",
        figsize = c(7, 6),
        conda_env = "base"
    )

    # Create a combined cell_type_group column in Seurat object
    seu$cell_type_group <- paste(seu$cell_type_dtl, seu$group, sep = "_")

    # Add the new column to adata
    adata.AddMetadata(seu, col = "cell_type_group", conda_env = "base")

    # You can also use python directly to verify/modify if needed
    py_run_string("print(f'Added cell_type_group column to adata with {len(adata.obs.cell_type_group.unique())} unique values')")
    # 2. Palantir Analysis
    # Run diffusion map
    seu <- Palantir.RunDM(seu,
        reduction = "harmony",
        conda_env = "base"
    )

    p <- DimPlot2(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/trajectory/plantir_dm.png", p, width = 7, height = 6)

    # Get cells from seu and update adata
    cells <- colnames(seu)
    py_run_string(sprintf("adata = adata[adata.obs.index.isin(%s)]", paste0("[", paste0("'", cells, "'", collapse = ","), "]")))

    # Add ms dimension reduction to AnnData object
    adata.AddDR(seu, dr = "ms", scv.graph = TRUE, conda_env = "base")

    # Plot velocity using ms coordinates
    scVelo.Plot(
        color = "cell_type_group",
        basis = "ms",
        save = "results/110.adipo/trajectory/velocity_ms.png",
        figsize = c(7, 6),
        conda_env = "base"
    )

    # Select start cell (Fib cluster as starting point)
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_group", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    start_cell <- CellSelector(p)
    # fate1_cell <- CellSelector(p)
    # fate2_cell <- CellSelector(p)
    # fate1_cell <- fate1_cell[1]
    # fate2_cell <- fate2_cell[2]
    # start_cell <- colnames(seu)[which(seu$cell_type_dtl == "Fib")[1]]

    # Calculate pseudotime
    seu <- Palantir.Pseudotime(seu,
        start_cell = start_cell, conda_env = "base"
        # terminal_states = c("fate1" = fate1_cell, "fate2" = fate2_cell)
    )

    # Get pseudotime data
    ps <- seu@misc$Palantir$Pseudotime
    colnames(ps)[3:4] <- c("Adipo_fate", "Fib_fate")
    seu@meta.data[, colnames(ps)] <- ps

    # Plot pseudotime and entropy
    p1 <- DimPlot2(seu,
        features = colnames(ps),
        reduction = "ms",
        cols = list(Entropy = "D"),
        theme = NoAxes()
    ) + theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/trajectory/palantir_pseudotime.png", p1, width = 12, height = 12)

    # 4. CellRank Analysis
    # Add pseudotime to adata
    adata.AddMetadata(seu, col = colnames(ps), conda_env = "base")

    # Run CellRank
    Cellrank.Compute(time_key = "Pseudotime", conda_env = "base")
    # SeuratExtend::adata.Save("data/110.adipo/trajectory/adipo.h5ad", conda_env = "base")
    # Generate CellRank plots
    Cellrank.Plot(
        color = "cell_type_group",
        basis = "ms",
        save = "results/110.adipo/trajectory/cellrank_ms.png",
        conda_env = "base"
    )

    return(seu)
}

classify_adipocyte_type <- function(sc_adipo) {
    # Create directory if it doesn't exist
    dir.create("results/110.adipo/adipocyte_types", showWarnings = FALSE, recursive = TRUE)

    sc_adipo <- sc_adipo %>%
        filter(cell_type_dtl == "Adipo")
    # Define marker genes for each adipocyte type based on the provided image
    adipocyte_markers <- list(
        "Brown adipocyte markers" = c("UCP1", "UCP", "PPARGC1A", "PPARGC1", "CIDEA", "EBF2", "PRDM16", "DIO2", "ZIC1", "LHX8"),
        "White adipocyte markers" = c("LEP", "HOXC9", "HOX3", "HOX3B", "TCF21", "LPL", "NRIP1", "RB1", "FABP4"),
        "Beige/Brite adipocyte markers" = c("UCP1", "UCP", "PPARGC1A", "PPARGC1", "TNFRSF9", "CITED1", "HOXC9", "HOX3", "HOX3B", "TMEM26")
    )

    # Create a dotplot showing expression of adipocyte markers across cell types
    p_markers <- DotPlot2(sc_adipo,
        features = adipocyte_markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(title = "Adipocyte Type Marker Expression")

    ggsave("results/110.adipo/adipocyte_types/adipocyte_markers_dotplot.png", p_markers, width = 10, height = 6)

    # Calculate adipocyte type scores for each cell
    sc_adipo <- AddModuleScore(
        sc_adipo,
        features = list(
            Brown = c("UCP1", "UCP", "PPARGC1A", "PPARGC1", "CIDEA", "EBF2", "PRDM16", "DIO2", "ZIC1", "LHX8"),
            White = c("LEP", "HOXC9", "HOX3", "HOX3B", "TCF21", "LPL", "NRIP1", "RB1", "FABP4"),
            Beige = c("UCP1", "UCP", "PPARGC1A", "PPARGC1", "TNFRSF9", "CITED1", "HOXC9", "HOX3", "HOX3B", "TMEM26")
        ),
        name = "AdipocyteType_",
        seed = 123
    )

    # Normalize scores to 0-1 range for better comparison
    score_cols <- c("AdipocyteType_1", "AdipocyteType_2", "AdipocyteType_3")

    # Min-max scaling for each score
    for (col in score_cols) {
        min_val <- min(sc_adipo@meta.data[[col]])
        max_val <- max(sc_adipo@meta.data[[col]])
        sc_adipo@meta.data[[col]] <- (sc_adipo@meta.data[[col]] - min_val) / (max_val - min_val)
    }

    # Plot adipocyte type scores on UMAP
    p_scores <- FeaturePlot(
        sc_adipo,
        features = score_cols,
        cols = c("lightgrey", "brown", "navy", "gold"),
        order = TRUE,
        ncol = 3
    ) & theme(plot.title = element_text(size = 10))

    # Rename the plots for clarity
    for (i in 1:length(p_scores)) {
        if (i == 1) p_scores[[i]] <- p_scores[[i]] + labs(title = "Brown Adipocyte Score")
        if (i == 2) p_scores[[i]] <- p_scores[[i]] + labs(title = "White Adipocyte Score")
        if (i == 3) p_scores[[i]] <- p_scores[[i]] + labs(title = "Beige/Brite Adipocyte Score")
    }

    ggsave("results/110.adipo/adipocyte_types/adipocyte_type_scores.png", p_scores, width = 12, height = 4)

    # Classify cells based on highest adipocyte type score
    sc_adipo@meta.data$adipocyte_type <- apply(
        sc_adipo@meta.data[, score_cols], 1,
        function(x) {
            scores <- as.numeric(x)
            types <- c("Brown", "White", "Beige")

            # Check if scores exceed minimum threshold
            min_threshold <- 0.2
            max_score <- max(scores)

            if (max_score < min_threshold) {
                return("Non-adipocyte")
            } else {
                # Sort scores in descending order and get indices
                sorted_indices <- order(scores, decreasing = TRUE)
                top_score <- scores[sorted_indices[1]]
                second_score <- scores[sorted_indices[2]]

                # Define threshold for considering mixed type
                # If difference between top scores is less than 20% of the top score
                score_diff_threshold <- 0.2

                if ((top_score - second_score) < (score_diff_threshold * top_score)) {
                    # It's a mixed type, combine the two types
                    return(paste(types[sorted_indices[1]], types[sorted_indices[2]], sep = "/"))
                } else {
                    # Clear dominant score
                    return(types[sorted_indices[1]])
                }
            }
        }
    )

    # Visualize the classification results on UMAP
    p_classification <- DimPlot(
        sc_adipo,
        group.by = "adipocyte_type",
        cols = c("Brown" = "#8C3310", "White" = "#2166AC", "Beige/Brite" = "#D9A441", "Non-adipocyte" = "grey"),
        label = TRUE,
        repel = TRUE
    ) +
        labs(title = "Adipocyte Type Classification") +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/110.adipo/adipocyte_types/adipocyte_classification.png", p_classification, width = 8, height = 6)

    # Generate a Venn diagram showing overlap between adipocyte types
    # First, determine which cells express each type (using a threshold)
    # Create cell lists based on adipocyte type classifications
    brown_cells <- rownames(sc_adipo@meta.data)[sc_adipo$adipocyte_type %in% c("Brown", "Brown/Beige", "Brown/White", "Beige/Brown", "White/Brown")]
    white_cells <- rownames(sc_adipo@meta.data)[sc_adipo$adipocyte_type %in% c("White", "White/Beige", "White/Brown", "Beige/White", "Brown/White")]
    beige_cells <- rownames(sc_adipo@meta.data)[sc_adipo$adipocyte_type %in% c("Beige", "Beige/Brown", "Beige/White", "Brown/Beige", "White/Beige")]

    # Create list for ggVennDiagram package
    library(ggVennDiagram)

    # Create the cell lists for the Venn diagram
    cell_lists <- list(
        "Brown" = brown_cells,
        "White" = white_cells,
        "Beige/Brite" = beige_cells
    )

    # Create and save the Venn diagram using ggVennDiagram
    p_venn <- ggVennDiagram(
        cell_lists,
        category.names = c("Brown", "White", "Beige/Brite")
        # set_color = c("#8C3310", "#2166AC", "#D9A441"),
    ) +
        scale_fill_gradient(low = "white", high = "#2166ac") +
        labs(title = "Overlap of Adipocyte Type Markers") +
        theme(plot.background = element_rect(fill = "white"))

    # Save the plot
    ggsave("results/110.adipo/adipocyte_types/adipocyte_venn_diagram.png",
        p_venn,
        width = 8, height = 8, dpi = 300
    )

    # Generate a violin plot showing adipocyte type scores by cell type
    score_data <- sc_adipo@meta.data %>%
        select(cell_type_dtl, group, AdipocyteType_1, AdipocyteType_2, AdipocyteType_3, adipocyte_type) %>%
        pivot_longer(
            cols = c(AdipocyteType_1, AdipocyteType_2, AdipocyteType_3),
            names_to = "score_type",
            values_to = "score"
        ) %>%
        mutate(score_type = case_when(
            score_type == "AdipocyteType_1" ~ "Brown Score",
            score_type == "AdipocyteType_2" ~ "White Score",
            score_type == "AdipocyteType_3" ~ "Beige Score"
        ))

    p_violin <- ggplot(score_data, aes(x = cell_type_dtl, y = score, fill = score_type)) +
        geom_violin(scale = "width", adjust = 1.5, trim = FALSE) +
        facet_wrap(~group, ncol = 1) +
        scale_fill_manual(values = c("Brown Score" = "#8C3310", "White Score" = "#2166AC", "Beige Score" = "#D9A441")) +
        theme_bw() +
        labs(x = "Cell Type", y = "Adipocyte Type Score", title = "Adipocyte Type Scores by Cell Type") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        )

    ggsave("results/110.adipo/adipocyte_types/adipocyte_score_violin.png", p_violin, width = 7, height = 8)

    p <- DotPlot2(sc_adipo %>% filter(adipocyte_type != "Non-adipocyte"),
        features = adipocyte_markers,
        group.by = "adipocyte_type"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/adipocyte_types/adipocyte_type_markers_dot.png", p, width = 7, height = 6)
    # Generate summary statistics
    classification_summary <- sc_adipo@meta.data %>%
        group_by(cell_type_dtl, group, adipocyte_type) %>%
        summarise(
            count = n(),
            mean_brown_score = mean(AdipocyteType_1),
            mean_white_score = mean(AdipocyteType_2),
            mean_beige_score = mean(AdipocyteType_3),
            mean_mito_content = mean(percent.mt)
        ) %>%
        ungroup()

    write.csv(classification_summary, "results/110.adipo/adipocyte_types/adipocyte_classification_summary.csv", row.names = FALSE)

    # Return the updated Seurat object with adipocyte type classification
    return(sc_adipo)
}
