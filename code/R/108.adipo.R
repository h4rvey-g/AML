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

cluster_adipo <- function(sc_final) {
    dir.create("results/110.adipo", showWarnings = FALSE, recursive = TRUE)
    # Filter for plasticity cells only (Fib, PSC, adipo)
    sc_adipo <- sc_final %>%
        filter(cell_type %in% c("Plasticity"))

    # Run UMAP
    sc_adipo <- RunUMAP(sc_adipo, dims = 1:30, reduction = "harmony", verbose = FALSE)

    # Plot UMAP
    p <- DimPlot2(sc_adipo, group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/110.adipo/adipo_umap.png", p, width = 7, height = 6)

    return(sc_adipo)
}
run_adipo_trajectory <- function(sc_adipo) {
    dir.create("results/110.adipo/trajectory", showWarnings = FALSE, recursive = TRUE)
    tar_load(sc_adipo)
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