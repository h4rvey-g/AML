adipo_DEG_plot <- function(sc_adipo) {
    # Create directory if it doesn't exist
    dir.create("results/110.adipo/DEG", showWarnings = FALSE, recursive = TRUE)

    markers <- list(
        "Adipogenic Program" = c("PPARG", "EBF3", "LPL", "LEPR", "CD36"),
        "Angiogenesis" = c("ANGPT1", "VEGFC", "ADAMTS2", "TEK"),
        "ECM Remodeling" = c("VCAN", "MFAP5", "ADAMTS2", "ADAMTS17", "ITIH5"),
        "Downregulated Developmental" = c("LSAMP", "ZFPM2", "KIAA1217", "NR2F1", "MEIS1", "ADAMTSL3"),
        "Upregulated Developmental" = c("SOX5", "EYA1", "RUNX1", "NCKAP5", "CHSY3")
    )

    p <- DotPlot2(sc_adipo,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/DEG/adipo_dotplot.png", p, width = 7, height = 7)

    # Filter for PSC and Fib cells only
    sc_adipo_select <- sc_adipo %>%
        filter(cell_type_dtl %in% c("PSC", "Fib"))

    p <- DotPlot2(sc_adipo_select,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        split.by.method = "color",
        split.by.colors = c("#4DBBD5FF", "#E64B35FF"),
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/DEG/adipo_select_dotplot.png", p, width = 7, height = 7)

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
    sc_adipo$cell_type_group <- paste(sc_adipo$cell_type_dtl, sc_adipo$group, sep = "_")
    Idents(sc_adipo) <- "cell_type_group"

    # Function to extract DEGs for a specific cell type
    get_cell_type_degs <- function(cell_type) {
        output_file <- paste0("results/110.adipo/DEG/", cell_type, "_DEGs.tsv")

        # Check if the file already exists
        if (file.exists(output_file)) {
            message(paste0("Loading existing DEGs from ", output_file))
            return(read_tsv(output_file))
        }

        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")

        if ((ident1 %in% Idents(sc_adipo) && ident2 %in% Idents(sc_adipo))) {
            results <- FindMarkers(
                sc_adipo,
                ident.1 = ident1,
                ident.2 = ident2,
                test.use = "MAST",
                min.cells.group = 2
            ) %>%
                Add_Pct_Diff() %>%
                filter(p_val_adj < 0.05) %>%
                rownames_to_column("gene") %>%
                as_tibble() %>%
                mutate(cell_type = cell_type)

            # Save the results to a TSV file
            write_tsv(results, output_file)

            return(results)
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
    deg_psc <- get_cell_type_degs("PSC")
    deg_fib <- get_cell_type_degs("Fib")

    # Combine DEGs from both cell types
    all_degs <- bind_rows(deg_psc, deg_fib)

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
    for (ct in c("PSC", "Fib")) {
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
            labs(x = "Log2 Fold Change (Tumor/Normal)", y = "", title = paste(ct)) +
            theme_bw() +
            theme(
                strip.background = element_rect(fill = "lightgrey"),
                strip.text = element_text(face = "bold"),
                panel.grid.minor = element_blank(),
                legend.position = "right",
                plot.background = element_rect(fill = "white")
            )

        filename <- paste0("results/110.adipo/DEG/adipo_", tolower(ct), "_dotplot.png")
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
            values = c("PSC" = 16, "Fib" = 17)
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
    ggsave("results/110.adipo/DEG/adipo_combined_dotplot.png", p_combined, width = 8, height = 9)

    return("results/110.adipo/DEG")
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
        FCcutoff = 0.5, # Adjust logFC cutoff as needed
        pointSize = 2.0,
        labSize = 4.0,
        colAlpha = 0.5,
        legendPosition = "none", # Remove legend
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        colConnectors = "grey50",
        boxedLabels = TRUE, # Make labels boxed
        # Define custom colors for labeled points (e.g., highlight in orange)
        # This overrides the default coloring for the selected labels
        colCustom = data.frame(
            key = deg_results$gene,
            color = ifelse(deg_results$gene %in% genes_to_label, 'orange', 'black')
        ) %>% tibble::deframe()
        # Removed the redundant 'col' argument here
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
    py_run_string("")
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
        basis = "umap_cell_embeddings",
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
    fate1_cell <- CellSelector(p)
    fate2_cell <- CellSelector(p)
    fate1_cell <- fate1_cell[1]
    fate2_cell <- fate2_cell[2]
    # start_cell <- colnames(seu)[which(seu$cell_type_dtl == "Fib")[1]]

    # Calculate pseudotime
    seu <- Palantir.Pseudotime(seu,
        start_cell = start_cell, conda_env = "base",
        terminal_states = c("fate1" = fate1_cell, "fate2" = fate2_cell)
    )

    # Get pseudotime data
    ps <- seu@misc$Palantir$Pseudotime
    # colnames(ps)[3, 4] <- c("fate1", "fate2")
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
                    return(paste(types[sorted_indices[1]], types[sorted_indices[2]], sep="/"))
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
           p_venn, width = 8, height = 8, dpi = 300)

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
