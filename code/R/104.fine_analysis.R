run_palantir <- function(sc_final) {
    library(SeuratExtend)
    library(Seurat)
    library(tidyseurat)
    tar_load(sc_final)
    sc_final <- sc_final %>%
        filter(cell_type %in% c("Plasticity"))
    sc_final <- Palantir.RunDM(sc_final, conda_env = "base", reduction = "harmony")
    # create a var, combine cell_type_dtl and group
    sc_final <- sc_final %>%
        mutate(cell_type_group = paste(cell_type_dtl, group, sep = "_"))
    p <- DimPlot2(sc_final, reduction = "ms", group.by = c("cell_type_group", "cell_type_dtl"), label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/palantir.png", p, width = 12, height = 6)
    p <- DimPlot(sc_final, reduction = "ms", group.by = c("cell_type_dtl"))
    cells <- CellSelector(p)
    sc_final <- Palantir.Pseudotime(sc_final, start_cell = cells, conda_env = "base")
    ps <- sc_final@misc$Palantir$Pseudotime
    head(ps)
    colnames(ps)[3:4] <- c("fate1", "fate2")
    sc_final@meta.data[, colnames(ps)] <- ps
    p <- DimPlot2(sc_final,
        features = colnames(ps), reduction = "ms",
        cols = list(Entropy = "D"), theme = NoAxes()
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/palantir_pseudotime.png", p, width = 12, height = 12)
    adata.AddDR(sc_final,
        dr = "ms", scv.graph = TRUE, load.adata = "data/104.RNA_velocity/anndata_filtered.h5ad",
        save.adata = "data/106.fine_analysis/adata_palantir.h5ad",
        conda_env = "base"
    )
    adata.AddMetadata(sc_final,
        col = "cell_type_group", load.adata = "data/106.fine_analysis/adata_palantir.h5ad",
        save.adata = "data/106.fine_analysis/adata_palantir.h5ad",
        conda_env = "base"
    )
    scVelo.Plot(
        load.adata = "data/106.fine_analysis/adata_palantir.h5ad", basis = "ms",
        color = "cell_type_dtl", save = "results/106.fine_analysis/velocity.png", figsize = c(5, 4),
        conda_env = "base"
    )
    scVelo.Plot(
        load.adata = "data/106.fine_analysis/adata_palantir.h5ad", basis = "ms",
        color = "cell_type_group", save = "results/106.fine_analysis/velocity_group.png", figsize = c(5, 4),
        conda_env = "base"
    )
    genes_up <- c(
        "LEPR", # Leptin receptor, matrix remodeling
        "CNTNAP2", # Cell adhesion and signaling
        "EYA1", # Transcriptional regulator
        "PDZRN4", # Protein-protein interactions
        "COL11A1", # ECM organization
        "DCDC1", # Cell differentiation
        "PPARG", # Key transcription factor
        "RUNX1", # Master regulator of differentiation
        "VCAM1", # Vascular cell adhesion
        "PTK2B", # Focal adhesion signaling
        "VEGFC", # Angiogenesis
        "THY1", # CAF marker
        "COL1A1", # ECM production
        "FBN1", # ECM organization
        "MARCKS"
    ) # Cell adhesion and motility

    # Top 15 downregulated genes contributing to Fib differentiation
    genes_down <- c(
        "LINC02388", # Regulatory RNA
        "CYP11B1", # Metabolic regulation
        "ADAMTSL2", # ECM modulation
        "LSAMP", # Cell adhesion
        "PTH1R", # Hormone signaling
        "CYP17A1", # Steroid biosynthesis
        "RALYL", # RNA binding
        "GUCY1A2", # Signaling pathway
        "HHIP", # Hedgehog pathway inhibitor
        "GATA6", # Transcription factor
        "WFDC1", # Tumor suppressor
        "NR2F1", # Nuclear receptor
        "ADAMTS12", # ECM organization
        "SEMA3B", # Tumor suppressor
        "MEIS1"
    ) # Transcription factor
    p1 <- GeneTrendHeatmap.Palantir(
        sc_final,
        features = genes_up,
        pseudotime.data = ps,
        magic = FALSE,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/106.fine_analysis/gene_trend_up.png", p1, width = 12, height = 6)
    p2 <- GeneTrendHeatmap.Palantir(
        sc_final,
        features = genes_down,
        pseudotime.data = ps,
        magic = FALSE,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/106.fine_analysis/gene_trend_down.png", p2, width = 12, height = 6)
    ggsave("results/106.fine_analysis/gene_trend.png", p1 / p2, width = 12, height = 12)

    # add col ps
    adata.AddMetadata(sc_final,
        col = colnames(ps), load.adata = "data/106.fine_analysis/adata_palantir.h5ad",
        save.adata = "data/106.fine_analysis/adata_palantir.h5ad", conda_env = "base"
    )
    Cellrank.Compute(
        load.adata = "data/106.fine_analysis/adata_palantir.h5ad", time_key = "Pseudotime", conda_env = "base"
    )
    Cellrank.Plot(
        load.adata = "data/106.fine_analysis/adata_palantir.h5ad", basis = "ms",
        color = "cell_type_dtl", save = "results/106.fine_analysis/cellrank.png", figsize = c(5, 4),
        conda_env = "base"
    )
    Cellrank.Plot(
        load.adata = "data/106.fine_analysis/adata_palantir.h5ad", basis = "ms",
        color = "cell_type_group", save = "results/106.fine_analysis/cellrank_group.png", figsize = c(5, 4),
        conda_env = "base"
    )
    qs::qsave(sc_final, "data/106.fine_analysis/sc_final_palantir.qs")
}

cell_communication <- function(sc_final, sc_mye, sc_tcell) {
    sc_final <- sc_final %>%
        filter(group == "tumor")
    library(liana)
    library(magrittr)
    library(nichenetr)

    # Load myeloid and T cell data if needed
    tar_load(sc_mye)
    tar_load(sc_tcell)

    # Extract cell_type_dtl from myeloid and T cell objects
    mye_meta <- sc_mye@meta.data %>%
        dplyr::select(cell_type_dtl) %>%
        tibble::rownames_to_column("cell_id")

    tcell_meta <- sc_tcell@meta.data %>%
        dplyr::select(cell_type_dtl) %>%
        tibble::rownames_to_column("cell_id")

    # Get current metadata
    current_meta <- sc_final@meta.data %>%
        tibble::rownames_to_column("cell_id")

    # Join the detailed cell types
    updated_meta <- current_meta %>%
        dplyr::left_join(mye_meta, by = "cell_id", suffix = c("", "_mye")) %>%
        dplyr::left_join(tcell_meta, by = "cell_id", suffix = c("", "_tcell"))

    # Update cell_type_dtl where needed
    updated_meta <- updated_meta %>%
        dplyr::mutate(cell_type_dtl = case_when(
            !is.na(cell_type_dtl_tcell) ~ cell_type_dtl_tcell,
            !is.na(cell_type_dtl_mye) ~ cell_type_dtl_mye,
            TRUE ~ cell_type_dtl
        )) %>%
        dplyr::select(-cell_type_dtl_mye, -cell_type_dtl_tcell) %>%
        tibble::column_to_rownames("cell_id")

    # Update sc_final metadata
    sc_final@meta.data <- updated_meta
    liana_test <- liana_wrap(sc_final, idents_col = "cell_type_dtl")
    liana_test <- liana_test %>%
        liana_aggregate()
    library(furrr)
    plan(multisession, workers = 20)

    unique_cell_types <- unique(sc_final$cell_type_dtl)

    liana_dotplot <- function(liana_res,
                              source_groups = NULL,
                              target_groups = NULL,
                              ntop = NULL,
                              specificity = "natmi.edge_specificity",
                              magnitude = "sca.LRscore",
                              y.label = "Interactions (Ligand -> Receptor)",
                              size.label = "Interaction\nSpecificity",
                              colour.label = "Expression\nMagnitude",
                              show_complex = TRUE,
                              size_range = c(2, 10),
                              invert_specificity = FALSE,
                              invert_magnitude = FALSE,
                              facet_by = "source",
                              invert_function = function(x) -log10(x + 1e-10)) {
        if (show_complex) {
            entities <- c("ligand.complex", "receptor.complex")
        } else {
            entities <- c("ligand", "receptor")
        }

        # Modify for the plot
        liana_mod <- liana_res %>%
            # Filter to only the cells of interest
            `if`(
                !is.null(source_groups),
                filter(., source %in% source_groups),
                .
            ) %>%
            `if`(
                !is.null(target_groups),
                filter(., target %in% target_groups),
                .
            )


        if (!is.null(ntop)) {
            # Subset to the X top interactions
            top_int <- liana_mod %>%
                distinct_at(entities) %>%
                head(ntop)
            liana_mod %<>% inner_join(top_int, by = entities)
        }

        if (invert_magnitude) {
            liana_mod %<>% mutate(!!magnitude := invert_function(.data[[magnitude]]))
        }
        if (invert_specificity) {
            liana_mod %<>% mutate(!!specificity := invert_function(.data[[specificity]]))
        }

        liana_mod %<>%
            dplyr::rename(magnitude = !!magnitude) %>%
            dplyr::rename(specificity = !!specificity) %>%
            unite(entities, col = "interaction", sep = " -> ") %>%
            unite(c("source", "target"), col = "source_target", remove = FALSE)



        # ensure levels & order is kept the plot
        interactions_order <- liana_mod %>%
            pull("interaction") %>%
            unique()
        liana_mod %<>%
            mutate(interaction = factor(interaction, levels = rev(interactions_order))) %>%
            mutate(across(where(is.character), as.factor))

        # colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
        cbPalette <- c(
            "#E69F00", "#56B4E9",
            "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#DF69A7"
        )

        # plot
        suppressWarnings(
            ggplot(
                liana_mod,
                aes(
                    x = if (facet_by == "source") target else source,
                    y = interaction,
                    colour = magnitude,
                    size = specificity,
                    group = if (facet_by == "source") target else source
                )
            ) +
                geom_point() +
                scale_color_gradientn(colours = viridis::viridis(20)) +
                scale_size_continuous(range = size_range) +
                facet_grid(. ~ get(facet_by),
                    space = "free",
                    scales = "free",
                    switch = "y"
                ) +
                labs(
                    y = y.label,
                    colour = colour.label,
                    size = size.label,
                    x = if (facet_by == "source") "Target" else "Source",
                    title = stringr::str_to_title(facet_by)
                ) +
                theme_bw(base_size = 20) +
                theme(
                    legend.text = element_text(size = 16),
                    axis.text.x = element_text(
                        colour =
                            cbPalette[1:length(
                                unique(if (facet_by == "source") liana_mod$source else liana_mod$target)
                            )],
                        face = "bold",
                        size = 23
                    ),
                    axis.title.x = element_text(colour = "gray6"),
                    axis.text.y = element_text(
                        size = 18,
                        vjust = 0.5
                    ),
                    legend.title = element_text(size = 18),
                    panel.spacing = unit(0.1, "lines"),
                    strip.background = element_rect(fill = NA),
                    plot.title = element_text(vjust = 0, hjust = 0.5, colour = "gray6"),
                    strip.text = element_text(size = 24, colour = "gray6") # ,
                    # strip.text.y.left = element_text(angle = 0)
                )
        )
    }

    options(future.globals.maxSize = 2000 * 1024^2)
    future_map(unique_cell_types, function(cell_type) {
        # Replace any slashes with underscores in the cell type name for the filename
        cell_type_safe <- gsub("/", "_", cell_type)
        tryCatch(
            {
                p <- liana_test %>%
                    liana_dotplot(
                        source_groups = c(cell_type),
                        ntop = 20
                    ) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(paste0("results/106.fine_analysis/liana_dotplot_send_", cell_type_safe, ".png"), p, width = 12, height = 10)
            },
            error = function(e) {
                message(paste("Error processing cell type", cell_type, ":", e$message))
            }
        )
        # do receiver
        tryCatch(
            {
                p <- liana_test %>%
                    liana_dotplot(
                        target_groups = c(cell_type),
                        ntop = 20,
                        facet_by = "target"
                    ) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                ggsave(paste0("results/106.fine_analysis/liana_dotplot_receive_", cell_type_safe, ".png"), p, width = 12, height = 10)
            },
            error = function(e) {
                message(paste("Error processing cell type", cell_type, ":", e$message))
            }
        )
    })
    key_types <- c("Adipo", "Fib", "PSC", "MSC", "Endo")
    p <- liana_test %>%
        liana_dotplot(
            source_groups = c(key_types),
            target_groups = c(key_types),
            ntop = 20
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("results/106.fine_analysis/liana_dotplot_key.png", p, width = 24, height = 10)
    liana_trunc <- liana_test %>%
        # only keep interactions concordant between methods
        filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

    png("results/106.fine_analysis/liana_heatmap.png", width = 12, height = 10, units = "in", res = 300)
    heat_freq(liana_trunc)
    dev.off()
    liana_trunc_filtered <- liana_test %>%
        filter(aggregate_rank <= 0.01, source %in% key_types, target %in% key_types)
    png("results/106.fine_analysis/liana_heatmap_filtered.png", width = 7, height = 6, units = "in", res = 300)
    heat_freq(liana_trunc_filtered)
    dev.off()
    p <- chord_freq(
        liana_trunc
        # source_groups = c("CD8 T", "NK", "B"),
        # target_groups = c("CD8 T", "NK", "B")
    )
    png("results/106.fine_analysis/liana_chord.png", width = 12, height = 10, units = "in", res = 300)
    p
    dev.off()
    p <- chord_freq(
        liana_trunc_filtered
        # source_groups = c("CD8 T", "NK", "B"),
        # target_groups = c("CD8 T", "NK", "B")
    )
    png("results/106.fine_analysis/liana_chord_filtered.png", width = 7, height = 6, units = "in", res = 300)
    p
    dev.off()
    "data/106.fine_analysis"
}

run_nichenet <- function(sc_final) {
    # 1. Define receiver and sender cell populations
    receiver <- "Endo"
    sender <- "Fib"
    Idents(sc_final) <- "cell_type_dtl"
    expressed_genes_receiver <- get_expressed_genes(receiver, sc_final, pct = 0.05)

    lr_network <- readRDS("data/reference_data/lr_network_human_21122021.rds")
    ligand_target_matrix <- readRDS("data/reference_data/ligand_target_matrix_nsga2r_final.rds")
    weighted_networks <- readRDS("data/reference_data/weighted_networks_nsga2r_final.rds")
    # Get potential ligands for sender-agnostic approach
    all_receptors <- unique(lr_network$to)
    expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
    potential_ligands <- lr_network %>%
        filter(to %in% expressed_receptors) %>%
        pull(from) %>%
        unique()

    # Get expressed genes in sender cells for sender-focused approach
    expressed_genes_sender <- get_expressed_genes(sender, sc_final, pct = 0.05)
    potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

    # 2. Define gene set of interest (DE genes in receiver)
    seurat_obj_receiver <- sc_final %>%
        filter(cell_type_dtl == receiver)
    DE_table_receiver <- FindMarkers(
        object = seurat_obj_receiver,
        ident.1 = "tumor",
        ident.2 = "normal",
        group.by = "group",
        min.pct = 0.05
    ) %>%
        rownames_to_column("gene")

    geneset_oi <- DE_table_receiver %>%
        filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
        pull(gene)
    geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

    # 3. Define background genes
    background_expressed_genes <- expressed_genes_receiver %>%
        .[. %in% rownames(ligand_target_matrix)]

    # 4. Perform NicheNet ligand activity analysis
    ligand_activities <- predict_ligand_activities(
        geneset = geneset_oi,
        background_expressed_genes = background_expressed_genes,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
    )

    ligand_activities <- ligand_activities %>%
        arrange(-aupr_corrected) %>%
        mutate(rank = rank(desc(aupr_corrected)))

    # Get top ligands
    best_upstream_ligands <- ligand_activities %>%
        top_n(30, aupr_corrected) %>%
        arrange(-aupr_corrected) %>%
        pull(test_ligand)

    # 5. Visualization of ligand activity scores
    vis_ligand_aupr <- ligand_activities %>%
        filter(test_ligand %in% best_upstream_ligands) %>%
        column_to_rownames("test_ligand") %>%
        select(aupr_corrected) %>%
        arrange(aupr_corrected) %>%
        as.matrix()

    p1 <- make_heatmap_ggplot(vis_ligand_aupr,
        "Prioritized ligands", "Ligand activity",
        legend_title = "AUPR",
        color = "darkorange"
    ) +
        theme(axis.text.x.top = element_blank(), plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/ligand_activity.png", p1, width = 9, height = 6)
    # 6. Get active target genes
    active_ligand_target_links_df <- best_upstream_ligands %>%
        lapply(get_weighted_ligand_target_links,
            geneset = geneset_oi,
            ligand_target_matrix = ligand_target_matrix,
            n = 100
        ) %>%
        bind_rows() %>%
        drop_na()

    active_ligand_target_links <- prepare_ligand_target_visualization(
        ligand_target_df = active_ligand_target_links_df,
        ligand_target_matrix = ligand_target_matrix,
        cutoff = 0.33
    )

    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
    order_targets <- active_ligand_target_links_df$target %>%
        unique() %>%
        intersect(rownames(active_ligand_target_links))

    vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

    p2 <- make_heatmap_ggplot(vis_ligand_target,
        "Prioritized ligands", "Predicted target genes",
        color = "purple",
        legend_title = "Regulatory potential"
    ) +
        scale_fill_gradient2(low = "whitesmoke", high = "purple") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/ligand_target.png", p2, width = 24, height = 6)

    # 7. Get receptor information
    ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
        best_upstream_ligands,
        expressed_receptors,
        lr_network,
        weighted_networks$lr_sig
    )

    vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
        ligand_receptor_links_df,
        best_upstream_ligands,
        order_hclust = "both"
    )

    p3 <- make_heatmap_ggplot(t(vis_ligand_receptor_network),
        y_name = "Ligands",
        x_name = "Receptors",
        color = "mediumvioletred",
        legend_title = "Prior interaction potential"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/ligand_receptor.png", p3, width = 15, height = 6)

    # 8. Expression of ligands in sender cells
    DE_table_sender <- FindMarkers(
        object = sc_final %>% filter(cell_type_dtl == sender),
        ident.1 = "tumor",
        ident.2 = "normal",
        group.by = "group",
        min.pct = 0,
        logfc.threshold = 0,
        features = best_upstream_ligands
    ) %>%
        rownames_to_column("gene")

    vis_ligand_lfc <- DE_table_sender %>%
        column_to_rownames("gene") %>%
        select(avg_log2FC) %>%
        as.matrix()

    p4 <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
        "Prioritized ligands",
        "Log2 Fold Change in Sender",
        low_color = "midnightblue",
        mid_color = "white",
        mid = 0,
        high_color = "red",
        legend_title = "Log2 FC"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/ligand_expression.png", p4, width = 4, height = 6)
    # Print results
    print("Top 10 prioritized ligands:")
    print(head(ligand_activities, 10))
    "data/106.fine_analysis"
}
