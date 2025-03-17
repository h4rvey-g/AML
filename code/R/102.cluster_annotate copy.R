integrate_data_2 <- function(sc_int) {
    # reticulate::virtualenv_create("scvienv", packages = c("scvi-tools", "anndata"))
    # Ensure required libraries are loaded
    reticulate::use_virtualenv("./scvienv")

    # Use Seurat for variable gene selection
    sc_int <- NormalizeData(sc_int, normalization.method = "LogNormalize", scale.factor = 10000)
    sc_int <- FindVariableFeatures(sc_int, selection.method = "vst", nfeatures = 2000)
    top2000 <- head(VariableFeatures(sc_int), 2000)

    # Filter genes for scVI
    sc_int_filtered <- sc_int[top2000]

    # Convert to AnnData format
    adata <- scCustomize::as.anndata(sc_int_filtered,
        file_path = "data/111.annotate",
        file_name = "anndata.h5ad",
        main_layer = "counts",
        other_layers = c("data")
    )

    # Import scvi from Python via reticulate
    scvi <- reticulate::import("scvi")

    # Convert sparse matrix to CSR format for faster training
    scipy <- reticulate::import("scipy.sparse")
    if (scipy$issparse(adata$X)) {
        adata$X <- scipy$csr_matrix(adata$X)
    }

    # Set up anndata with batch key as dataset
    scvi$model$SCVI$setup_anndata(adata, batch_key = "dataset")
    scvi$settings$dl_num_workers <- 4L

    # Create and train the model
    model <- scvi$model$SCVI(adata, n_layers = as.integer(2), n_latent = as.integer(30), gene_likelihood = "nb")

    # Set training parameters to improve performance with more workers
    model$train(
        batch_size = as.integer(256)
    )

    # Get the latent representation
    latent <- model$get_latent_representation()

    # Put it back in Seurat object
    latent <- as.matrix(latent)
    rownames(latent) <- colnames(sc_int)
    sc_int[["scvi"]] <- CreateDimReducObject(
        embeddings = latent, key = "scvi_",
        assay = DefaultAssay(sc_int)
    )

    # Run UMAP on scVI embeddings
    sc_int <- RunUMAP(sc_int,
        dims = 1:30, reduction = "scvi",
        reduction.name = "umap_integrated",
        reduction.key = "umapin_"
    )

    # Create UMAP plot by dataset to check integration
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "dataset") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/integration_by_dataset.png", p, width = 7, height = 7)

    # Create UMAP plot by group
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/integration_by_group.png", p, width = 7, height = 7)

    # Return integrated Seurat object
    sc_int
}
cluster_data_2 <- function(sc_int) {
    # Separate clustering for normal samples
    sc_int_normal <- sc_int %>% filter(group == "normal")
    sc_int_normal <- sc_int_normal %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.3, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_normal")

    # Separate clustering for tumor samples
    sc_int_tumor <- sc_int %>% filter(group == "tumor")
    sc_int_tumor <- sc_int_tumor %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.3, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_tumor")

    # Add the tumor and normal specific clusters back to the main object
    sc_int$seurat_clusters_normal <- "NA"
    sc_int$seurat_clusters_tumor <- "NA"
    sc_int <- AddMetaData(sc_int,
        metadata = sc_int_normal$seurat_clusters_normal,
        col.name = "seurat_clusters_normal"
    )
    sc_int <- AddMetaData(sc_int,
        metadata = sc_int_tumor$seurat_clusters_tumor,
        col.name = "seurat_clusters_tumor"
    )

    # Plot normal clusters
    p_normal <- DimPlot2(sc_int_normal, reduction = "umap_integrated", group.by = "seurat_clusters_normal", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/cluster_normal.png", p_normal, width = 7, height = 7)

    # Plot tumor clusters
    p_tumor <- DimPlot2(sc_int_tumor, reduction = "umap_integrated", group.by = "seurat_clusters_tumor", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/cluster_tumor.png", p_tumor, width = 7, height = 7)

    sc_int
}

annotate_data_2 <- function(sc_int, cell_type_path) {
    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2", "HSD3B2", "CYP17A1", "CYP21A2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1"),
        # "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
        "CLC" = c("TH", "CHGA", "CHGB", "KIT"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
        "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
        "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM"),
        "MSC" = c("ENG", "NES", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3"),
        # "PSC_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
        # "PSC_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
        "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
        "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
    )

    # Process all data together
    toplot <- CalcStats(sc_int,
        features = markers %>% unlist(),
        method = "zscore", order = "value"
    )
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot))
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    gene_groups <- gene_groups[rownames(toplot)]
    p1 <- Heatmap(toplot, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_zscore.png", p1, width = 7, height = 14)

    # Process normal samples
    sc_int_normal <- sc_int %>% filter(group == "normal")
    toplot_normal <- CalcStats(sc_int_normal,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "seurat_clusters_normal"
    )
    # make toplog_normal colnames order as number
    toplot_normal <- toplot_normal[, order(as.numeric(colnames(toplot_normal)))]
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot_normal))
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot_normal)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    gene_groups <- gene_groups[rownames(toplot_normal)]
    p_normal <- Heatmap(toplot_normal, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_zscore_normal.png", p_normal, width = 7, height = 14)
    write_csv(
        toplot_normal %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/102.cluster_annotate/heatmap_zscore_normal.csv"
    )


    # Process tumor samples
    sc_int_tumor <- sc_int %>% filter(group == "tumor")
    toplot_tumor <- CalcStats(sc_int_tumor,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "seurat_clusters_tumor"
    )
    # make toplog_tumor colnames order as number
    toplot_tumor <- toplot_tumor[, order(as.numeric(colnames(toplot_tumor)))]
    missing_genes <- setdiff(markers %>% unlist(), rownames(toplot_tumor))
    markers_after <- markers %>% map(~ .x[.x %in% rownames(toplot_tumor)])
    gene_groups <- rep(names(markers_after), lengths(markers_after)) %>%
        setNames(markers_after %>% unlist())
    gene_groups <- gene_groups[rownames(toplot_tumor)]
    p_tumor <- Heatmap(toplot_tumor, lab_fill = "zscore", facet_row = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/heatmap_zscore_tumor.png", p_tumor, width = 7, height = 14)
    write_csv(
        toplot_tumor %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/102.cluster_annotate/heatmap_zscore_tumor.csv"
    )

    # Create comparison plot
    p_combined <- (p_normal + p_tumor) / p1
    ggsave("results/102.cluster_annotate/heatmap_zscore_comparison.png", p_combined, width = 14, height = 30)

    sc_int
}

primary_annotation <- function(sc_int) {
    # Add primary annotation based on scvi clusters
    cell_type_normal <- list(
        "Steroidogenic" = c("1", "2", "3", "8", "9", "17"),
        "Endo" = c("4", "13"),
        "PSC" = c("5"),
        "Myeloid" = c("6"),
        "Fib" = c("7"),
        "Tcell" = c("10"),
        "LEC" = c("11"),
        "CLC" = c("12"),
        "Neural/glial" = c("14"),
        "Plasma" = c("15"),
        "Adipo" = c("16")
    )
    cell_type_tumor <- list(
        "Tcell" = c("1"),
        "Bcell" = c("2"),
        "Stero/medulla" = c("4", "11"),
        "Fib/PSC" = c("5"),
        "Endo" = c("7"),
        "Myeloid" = c("3", "8"),
        "Eryth" = c("10"),
        "Plasma" = c("9"),
        "Adipo" = c("6")
    )
    # convert to mapping table
    # Create mapping tables for normal and tumor cell types
    cell_type_normal_map <- tibble(
        cluster = unlist(cell_type_normal),
        cell_type = rep(names(cell_type_normal), lengths(cell_type_normal))
    )

    cell_type_tumor_map <- tibble(
        cluster = unlist(cell_type_tumor),
        cell_type = rep(names(cell_type_tumor), lengths(cell_type_tumor))
    )

    # Add cell type annotations to Seurat object
    sc_int$cell_type_normal <- "Unknown"
    sc_int$cell_type_tumor <- "Unknown"

    sc_int$cell_type_normal <- cell_type_normal_map$cell_type[
        match(sc_int$seurat_clusters_normal, cell_type_normal_map$cluster)
    ]
    sc_int$cell_type_tumor <- cell_type_tumor_map$cell_type[
        match(sc_int$seurat_clusters_tumor, cell_type_tumor_map$cluster)
    ]

    # Create unified cell type column
    sc_int$cell_type <- ifelse(sc_int$group == "normal",
        sc_int$cell_type_normal,
        sc_int$cell_type_tumor
    )
    # Create UMAP plot of cell types
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "cell_type", label = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/primary_annotation.png", p, width = 7, height = 7)

    # Create separate UMAP plots for normal and tumor samples
    p_normal <- DimPlot2(sc_int %>% filter(group == "normal"),
        reduction = "umap_integrated",
        group.by = "cell_type_normal",
        label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/primary_annotation_normal.png", p_normal, width = 7, height = 7)

    p_tumor <- DimPlot2(sc_int %>% filter(group == "tumor"),
        reduction = "umap_integrated",
        group.by = "cell_type_tumor",
        label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/primary_annotation_tumor.png", p_tumor, width = 7, height = 7)

    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2", "HSD3B2", "CYP17A1", "CYP21A2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1"),
        # "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
        "CLC" = c("TH", "CHGA", "CHGB", "KIT"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
        "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
        "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM"),
        "MSC" = c("ENG", "NES", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3"),
        # "PSC_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
        # "PSC_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
        "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
        "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
    )
    p_normal <- DotPlot2(sc_int %>% filter(group == "normal"),
        features = markers,
        group.by = "cell_type",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/primary_annotation_normal_dotplot.png", p_normal, width = 14, height = 14)
    p_tumor <- DotPlot2(sc_int %>% filter(group == "tumor"),
        features = markers,
        group.by = "cell_type",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/primary_annotation_tumor_dotplot.png", p_tumor, width = 14, height = 14)

    # Return the annotated Seurat object
    return(sc_int)
}
optimize_annotation <- function(sc_int) {
    # Sub-cluster CLC cells in normal samples
    sc_int_CLC <- sc_int %>%
        filter(group == "normal" & cell_type == "CLC")

    if (nrow(sc_int_CLC@meta.data) > 0) {
        sc_int_CLC <- sc_int_CLC %>%
            FindNeighbors(dims = 1:30, reduction = "scvi") %>%
            FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_CLC")

        # Plot CLC sub-clusters
        p_CLC <- DimPlot2(sc_int_CLC,
            reduction = "umap_integrated",
            group.by = "seurat_clusters_CLC", label = TRUE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_CLC.png", p_CLC, width = 7, height = 7)

        # Generate marker heatmap for CLC sub-clusters
        markers_CLC <- list(
            "CLC" = c("TH", "CHGA", "CHGB", "KIT", "SYT1"),
            "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
        )

        toplot_CLC <- CalcStats(sc_int_CLC,
            features = markers_CLC %>% unlist(),
            group.by = "seurat_clusters_CLC",
            method = "zscore", order = "value"
        )

        gene_groups <- rep(names(markers_CLC), lengths(markers_CLC)) %>%
            setNames(markers_CLC %>% unlist())
        gene_groups <- gene_groups[rownames(toplot_CLC)]

        p_CLC_heat <- Heatmap(t(toplot_CLC), lab_fill = "zscore", facet_col = gene_groups) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_CLC_heatmap.png", p_CLC_heat, width = 10, height = 7)
        # Add the CLC sub-clusters back to the main object
        sc_int$seurat_clusters_CLC <- "NA"
        sc_int <- AddMetaData(sc_int,
            metadata = sc_int_CLC$seurat_clusters_CLC,
            col.name = "seurat_clusters_CLC"
        )
    }

    # Sub-cluster Stero/medulla cells in tumor samples
    sc_int_Stero <- sc_int %>%
        filter(group == "tumor" & cell_type == "Stero/medulla")

    if (nrow(sc_int_Stero@meta.data) > 0) {
        sc_int_Stero <- sc_int_Stero %>%
            FindNeighbors(dims = 1:30, reduction = "scvi") %>%
            FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_stero_medulla")

        # Plot steroidogenic sub-clusters
        p_stero <- DimPlot2(sc_int_Stero,
            reduction = "umap_integrated",
            group.by = "seurat_clusters_stero_medulla", label = TRUE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_stero_medulla.png", p_stero, width = 7, height = 7)

        # Generate marker heatmap for steroidogenic sub-clusters
        markers_stero <- list(
            "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2"),
            "ZG" = c("DACH1", "VSNL1"),
            "ZG/ZF" = c("CCN3", "NCAM1"),
            "ZF" = c("CYP11B1", "ABCB1"),
            "ZR" = c("CYB5A", "SULT2A1"),
            # "Cap" = c("NR2F2", "RSPO3", "SPARCL1", "COL1A1", "COL1A2"),
            "CLC" = c("TH", "CHGA", "CHGB", "KIT"),
            "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
            "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
            "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM"),
            "MSC" = c("ENG", "NES", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3"),
            # "PSC_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
            # "PSC_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
            "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
            "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
            "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
            "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
            "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
            "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
            "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21"),
            "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2"),
            "Mast" = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT")
        )

        toplot_stero <- CalcStats(sc_int_Stero,
            features = markers_stero %>% unlist(),
            group.by = "seurat_clusters_stero_medulla",
            method = "zscore", order = "value"
        )

        gene_groups <- rep(names(markers_stero), lengths(markers_stero)) %>%
            setNames(markers_stero %>% unlist())
        gene_groups <- gene_groups[rownames(toplot_stero)]

        p_stero_heat <- Heatmap(t(toplot_stero), lab_fill = "zscore", facet_col = gene_groups) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_stero_medulla_heatmap.png", p_stero_heat, width = 20, height = 7)

        p_stero_dot <- DotPlot2(sc_int_Stero,
            features = markers_stero,
            group.by = "seurat_clusters_stero_medulla",
            show_grid = FALSE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_stero_medulla_dotplot.png", p_stero_dot, width = 14, height = 14)

        # Find markers for each steroidogenic-medulla subcluster
        stero_medulla_clusters <- unique(sc_int_Stero$`seurat_clusters_stero_medulla`)
        message("Finding markers for ", length(stero_medulla_clusters), " steroidogenic-medulla subclusters")

        # Create directory for marker gene results
        dir.create("results/102.cluster_annotate/stero_medulla_markers", recursive = TRUE, showWarnings = FALSE)

        # Set identity for finding markers
        Idents(sc_int_Stero) <- "seurat_clusters_stero_medulla"

        for (cluster in stero_medulla_clusters) {
            message("Finding markers for cluster ", cluster)

            # Find markers for this cluster vs all others
            # Use all tumor cells as the reference for comparison
            Idents(sc_int_Stero) <- "seurat_clusters_stero_medulla"

            # Get all tumor cells for comparison
            all_tumor_cells <- sc_int %>%
                filter(group == "tumor" &
                    !colnames(sc_int) %in% colnames(sc_int_Stero[, sc_int_Stero$seurat_clusters_stero_medulla == cluster]))

            # Create a combined object for the comparison
            combined_obj <- merge(
                sc_int_Stero[, sc_int_Stero$seurat_clusters_stero_medulla == cluster],
                all_tumor_cells
            )

            # Join layers after merging
            combined_obj <- JoinLayers(combined_obj)

            # Set identities for the comparison
            combined_obj$comparison_group <- "other_tumor"
            combined_obj$comparison_group[1:length(which(sc_int_Stero$seurat_clusters_stero_medulla == cluster))] <- "target_cluster"
            Idents(combined_obj) <- "comparison_group"

            # Find markers comparing this cluster against all other tumor cells
            cluster_markers <- FindMarkers(
                combined_obj,
                ident.1 = "target_cluster",
                ident.2 = "other_tumor",
                min.pct = 0.2,
                test.use = "MAST",
                only.pos = TRUE
            ) %>%
                rownames_to_column("gene") %>%
                filter(p_val_adj < 0.05) %>%
                Add_Pct_Diff() %>%
                arrange(desc(avg_log2FC))

            # Save to TSV
            write_tsv(
                cluster_markers,
                file = paste0("results/102.cluster_annotate/stero_medulla_markers/cluster_", cluster, "_markers.tsv")
            )
        }

        # Also generate a combined visualization of top markers
        top_markers <- list()
        for (cluster in stero_medulla_clusters) {
            markers_file <- paste0("results/102.cluster_annotate/stero_medulla_markers/cluster_", cluster, "_markers.tsv")
            if (file.exists(markers_file)) {
                markers_df <- read_tsv(markers_file)
                top_markers[[paste0("Cluster_", cluster)]] <- head(markers_df$gene, 10)
            }
        }

        # Create dotplot of top markers for each cluster
        if (length(top_markers) > 0) {
            p_markers <- DotPlot2(sc_int_Stero,
                features = top_markers,
                group.by = "seurat_clusters_stero_medulla",
                show_grid = FALSE
            ) +
                theme(
                    plot.background = element_rect(fill = "white"),
                    axis.text.x = element_text(angle = 45, hjust = 1)
                )

            ggsave("results/102.cluster_annotate/stero_medulla_markers/top_markers_dotplot.png",
                p_markers,
                width = 16, height = 8
            )
        }
        # Store steroidogenic sub-cluster labels
        sc_int$seurat_clusters_stero_medulla <- "NA"
        sc_int$seurat_clusters_stero_medulla[match(colnames(sc_int_Stero), colnames(sc_int))] <-
            sc_int_Stero$seurat_clusters_stero_medulla
    }
    # Sub-cluster Fib/PSC cells in tumor samples
    sc_int_FibPSC <- sc_int %>%
        filter(group == "tumor" & cell_type == "Fib/PSC")

    if (nrow(sc_int_FibPSC@meta.data) > 0) {
        sc_int_FibPSC <- sc_int_FibPSC %>%
            FindNeighbors(dims = 1:30, reduction = "scvi") %>%
            FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_FibPSC")

        # Plot Fib/PSC sub-clusters
        p_fibpsc <- DimPlot2(sc_int_FibPSC,
            reduction = "umap_integrated",
            group.by = "seurat_clusters_FibPSC", label = TRUE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_FibPSC.png", p_fibpsc, width = 7, height = 7)

        # Generate marker heatmap for Fib/PSC sub-clusters
        markers_fibpsc <- list(
            "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
            "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM")
        )

        toplot_fibpsc <- CalcStats(sc_int_FibPSC,
            features = markers_fibpsc %>% unlist(),
            group.by = "seurat_clusters_FibPSC",
            method = "zscore", order = "value"
        )

        gene_groups <- rep(names(markers_fibpsc), lengths(markers_fibpsc)) %>%
            setNames(markers_fibpsc %>% unlist())
        gene_groups <- gene_groups[rownames(toplot_fibpsc)]

        p_fibpsc_heat <- Heatmap(t(toplot_fibpsc), lab_fill = "zscore", facet_col = gene_groups) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_FibPSC_heatmap.png", p_fibpsc_heat, width = 10, height = 7)

        # Store Fib/PSC sub-cluster labels
        sc_int$seurat_clusters_FibPSC <- "NA"
        sc_int$seurat_clusters_FibPSC[match(colnames(sc_int_FibPSC), colnames(sc_int))] <-
            sc_int_FibPSC$seurat_clusters_FibPSC
    }
    # Select the cells from normal PSC, normal Fib, and tumor Fib/PSC
    sc_int_mesen <- sc_int %>%
        filter((group == "normal" & cell_type %in% c("PSC", "Fib")) |
            (group == "tumor" & cell_type == "Fib/PSC"))

    # Special combined group for visualization
    sc_int_mesen$mesen_group <- paste0(sc_int_mesen$cell_type, "_", sc_int_mesen$group)

    # Define markers for mesenchymal cell comparison
    markers_mesen <- list(
        "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
        "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM"),
        "MSC" = c("ENG", "NES", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3")
    )

    # Calculate expression statistics by mesenchymal group
    toplot_mesen <- CalcStats(sc_int_mesen,
        features = markers_mesen %>% unlist(),
        group.by = "mesen_group",
        method = "zscore", order = "value"
    )

    # Create gene group annotations
    gene_groups <- rep(names(markers_mesen), lengths(markers_mesen)) %>%
        setNames(markers_mesen %>% unlist())
    gene_groups <- gene_groups[rownames(toplot_mesen)]

    # Generate heatmap
    p_mesen_heat <- Heatmap(t(toplot_mesen), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/mesenchymal_comparison_heatmap.png", p_mesen_heat, width = 10, height = 7)
    write_tsv(
        toplot_mesen %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/102.cluster_annotate/mesenchymal_comparison_heatmap.tsv"
    )

    # Generate dotplot
    p_mesen_dot <- DotPlot2(sc_int_mesen,
        features = markers_mesen,
        group.by = "mesen_group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/mesenchymal_comparison_dotplot.png", p_mesen_dot, width = 12, height = 8)

    # UMAP visualization of mesenchymal cells
    p_mesen_umap <- DimPlot2(sc_int_mesen,
        reduction = "umap_integrated",
        group.by = "mesen_group",
        label = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/mesenchymal_comparison_umap.png", p_mesen_umap, width = 7, height = 7)

    # Sub-cluster Steroidogenic cells in normal samples
    sc_int_Steroidogenic <- sc_int %>%
        filter(cell_type == "Steroidogenic")

    if (nrow(sc_int_Steroidogenic@meta.data) > 0) {
        sc_int_Steroidogenic <- sc_int_Steroidogenic %>%
            FindNeighbors(dims = 1:30, reduction = "scvi") %>%
            FindClusters(resolution = 0.5, algorithm = 4, method = "igraph", cluster.name = "seurat_clusters_steroidogenic_normal")

        # Plot Steroidogenic sub-clusters
        p_steroidogenic <- DimPlot2(sc_int_Steroidogenic,
            reduction = "umap_integrated",
            group.by = "seurat_clusters_steroidogenic_normal", label = TRUE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_steroidogenic_normal.png", p_steroidogenic, width = 7, height = 7)

        # Generate marker heatmap for Steroidogenic sub-clusters
        markers_steroidogenic <- list(
            "ZG" = c("NR5A1", "DACH1", "CYP11B2"),
            "ZG/ZF" = c("NR5A1", "NOV", "NCAM1"),
            "ZR" = c("NR5A1", "CYB5A", "SULT2A1")
        )

        toplot_steroidogenic <- CalcStats(sc_int_Steroidogenic,
            features = markers_steroidogenic %>% unlist(),
            group.by = "seurat_clusters_steroidogenic_normal",
            method = "zscore", order = "value"
        )

        gene_groups <- rep(names(markers_steroidogenic), lengths(markers_steroidogenic)) %>%
            setNames(markers_steroidogenic %>% unlist())
        gene_groups <- gene_groups[rownames(toplot_steroidogenic)]

        p_steroidogenic_heat <- Heatmap(t(toplot_steroidogenic), lab_fill = "zscore", facet_col = gene_groups) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_steroidogenic_heatmap.png", p_steroidogenic_heat, width = 15, height = 7)

        # Also create a dotplot for clearer visualization
        p_steroidogenic_dot <- DotPlot2(sc_int_Steroidogenic,
            features = markers_steroidogenic,
            group.by = "seurat_clusters_steroidogenic_normal",
            show_grid = FALSE
        ) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/sub_clusters_steroidogenic_dotplot.png", p_steroidogenic_dot, width = 12, height = 8)

        # Store Steroidogenic sub-cluster labels
        sc_int$seurat_clusters_steroidogenic_normal <- "NA"
        sc_int$seurat_clusters_steroidogenic_normal[match(colnames(sc_int_Steroidogenic), colnames(sc_int))] <-
            sc_int_Steroidogenic$seurat_clusters_steroidogenic_normal
    }
    # Compare Stero/medulla cells in tumor samples vs CLC cells in normal samples
    # Extract the relevant cells
    stero_medulla_cells <- sc_int %>% filter(group == "tumor" & cell_type == "Stero/medulla")
    clc_cells <- sc_int %>% filter(group == "normal" & cell_type == "CLC")

    if (nrow(stero_medulla_cells@meta.data) > 0 && nrow(clc_cells@meta.data) > 0) {
        # Merge the two cell populations
        comparison_obj <- merge(stero_medulla_cells, clc_cells)

        # Join layers if needed
        comparison_obj <- JoinLayers(comparison_obj)

        # Create a new grouping variable for DEG analysis
        comparison_obj$comparison_group <- ifelse(
            comparison_obj$group == "tumor" & comparison_obj$cell_type == "Stero/medulla",
            "Stero_medulla_tumor",
            "CLC_normal"
        )

        # Set identity for comparison
        Idents(comparison_obj) <- "comparison_group"

        # Find DEGs
        neural_vs_clc_degs <- FindMarkers(
            comparison_obj,
            ident.1 = "Stero_medulla_tumor",
            ident.2 = "CLC_normal",
            min.pct = 0.2,
            test.use = "MAST"
        ) %>%
            rownames_to_column("gene") %>%
            filter(p_val_adj < 0.05) %>%
            arrange(desc(avg_log2FC))

        # Save results
        write_tsv(
            neural_vs_clc_degs,
            file = "results/102.cluster_annotate/neural_vs_CLC.tsv"
        )

        message("Saved DEGs between Stero/medulla in tumor vs CLC in normal to neural_vs_CLC.tsv")
    } else {
        message("Not enough cells for Stero/medulla in tumor vs CLC in normal comparison")
    }
    return(sc_int)
}
final_annotation_2 <- function(sc_int) {
    # # Process CLC subclusters - change subclusters 1 and 2 to Eryth
    # clc_cells <- which(sc_int$cell_type == "CLC" & !is.na(sc_int$seurat_clusters_CLC))
    # if (length(clc_cells) > 0) {
    #     # Identify which cells are in subclusters 1 and 2
    #     eryth_cells <- which(sc_int$cell_type == "CLC" &
    #         sc_int$seurat_clusters_CLC %in% c("3"))

    #     # Change these cells to Eryth
    #     if (length(eryth_cells) > 0) {
    #         sc_int$cell_type[eryth_cells] <- "Eryth"
    #     }
    # }
    # Process steroidogenic subclusters in tumor samples
    if (any(grepl("Fib/PSC", sc_int$cell_type))) {
        # Change all Fib/PSC cells in tumor samples to Fib
        fib_psc_cells <- which(sc_int$cell_type == "Fib/PSC" & sc_int$group == "tumor")
        if (length(fib_psc_cells) > 0) {
            sc_int$cell_type[fib_psc_cells] <- "Fib"
        }
    }
    # Remove Neural/glial cells completely
    neural_glial_cells <- which(sc_int$cell_type == "Neural/glial")
    if (length(neural_glial_cells) > 0) {
        # Log how many cells will be removed
        message("Removing ", length(neural_glial_cells), " Neural/glial cells from the dataset")

        # Remove cells from Seurat object
        sc_int <- sc_int[, -neural_glial_cells]
    }
    # Change Eryth cells to Myeloid
    eryth_cells <- which(sc_int$cell_type == "Eryth")
    if (length(eryth_cells) > 0) {
        message("Changing ", length(eryth_cells), " cells from Eryth to Myeloid")
        sc_int$cell_type[eryth_cells] <- "Myeloid"
    }
    # Process Stero/medulla subclusters in tumor samples
    if ("seurat_clusters_stero_medulla" %in% colnames(sc_int@meta.data)) {
        # Update the cell types based on the subclusters

        # Change subclusters 1, 2 to Neural-like
        neural_like_cells <- which(sc_int$group == "tumor" &
            !is.na(sc_int$seurat_clusters_stero_medulla) &
            sc_int$seurat_clusters_stero_medulla %in% c("1", "2"))

        if (length(neural_like_cells) > 0) {
            message("Changing ", length(neural_like_cells), " cells from Stero/medulla (subclusters 1, 2) to Neural-like")
            sc_int$cell_type[neural_like_cells] <- "Neural-like"
        }

        # Change subclusters 3, 4 to Myeloid
        myeloid_cells <- which(sc_int$group == "tumor" &
            !is.na(sc_int$seurat_clusters_stero_medulla) &
            sc_int$seurat_clusters_stero_medulla %in% c("3", "4"))

        if (length(myeloid_cells) > 0) {
            message("Changing ", length(myeloid_cells), " cells from Stero/medulla (subclusters 3, 4) to Myeloid")
            sc_int$cell_type[myeloid_cells] <- "Myeloid"
        }
    }
    # Create a detailed cell type column based on cell_type and apply specific refinements
    sc_int$cell_type_dtl <- sc_int$cell_type

    # Handle Steroidogenic cell subtypes based on clusters
    if ("seurat_clusters_steroidogenic_normal" %in% colnames(sc_int@meta.data)) {
        # For normal samples with Steroidogenic cells, assign specific zones
        # ZG: clusters 1, 5, 7, 9
        # ZF: clusters 2, 4
        # ZR: clusters 3, 6, 8
        zg_clusters <- c("5", "7", "9")
        zf_clusters <- c("1", "2", "4")
        zr_clusters <- c("3", "6", "8")

        # Assign ZG cells
        zg_cells <- which(sc_int$cell_type == "Steroidogenic" &
            sc_int$seurat_clusters_steroidogenic_normal %in% zg_clusters)
        if (length(zg_cells) > 0) {
            sc_int$cell_type_dtl[zg_cells] <- "ZG"
        }

        # Assign ZF cells
        zf_cells <- which(sc_int$cell_type == "Steroidogenic" &
            sc_int$seurat_clusters_steroidogenic_normal %in% zf_clusters)
        if (length(zf_cells) > 0) {
            sc_int$cell_type_dtl[zf_cells] <- "ZF"
        }

        # Assign ZR cells
        zr_cells <- which(sc_int$cell_type == "Steroidogenic" &
            sc_int$seurat_clusters_steroidogenic_normal %in% zr_clusters)
        if (length(zr_cells) > 0) {
            sc_int$cell_type_dtl[zr_cells] <- "ZR"
        }
    }

    # Create heatmap to verify detailed steroidogenic cell types
    steroidogenic_markers <- list(
        "ZG" = c("NR5A1", "DACH1", "CYP11B2", "VSNL1"),
        "ZF" = c("NR5A1", "CYP11B1", "ABCB1"),
        "ZR" = c("NR5A1", "CYB5A", "SULT2A1")
    )

    # Plot only for steroidogenic subtypes
    sc_int_stero <- sc_int %>% filter(cell_type_dtl %in% c("ZG", "ZF", "ZR"))

    if (nrow(sc_int_stero@meta.data) > 0) {
        toplot_stero <- CalcStats(sc_int_stero,
            features = steroidogenic_markers %>% unlist(),
            group.by = "cell_type_dtl",
            method = "zscore",
            order = "value"
        )

        gene_groups <- rep(names(steroidogenic_markers), lengths(steroidogenic_markers)) %>%
            setNames(steroidogenic_markers %>% unlist())
        gene_groups <- gene_groups[rownames(toplot_stero)]

        p_stero <- Heatmap(t(toplot_stero), lab_fill = "zscore", facet_col = gene_groups) +
            theme(plot.background = element_rect(fill = "white"))
        ggsave("results/102.cluster_annotate/detailed_steroidogenic_verification.png", p_stero, width = 10, height = 7)
    }

    # Update plotting data to reflect the new cell type assignments

    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/final_annotation_2.png", p, width = 7, height = 7)
    p1 <- DimPlot2(sc_int %>% filter(group == "normal"),
        reduction = "umap_integrated", group.by = "cell_type_dtl",
        label = TRUE, repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white")) +
        ggtitle("Normal Cells - Detailed Cell Types")
    ggsave("results/102.cluster_annotate/final_annotation_normal_2.png", p1, width = 7, height = 7)
    p2 <- DimPlot2(sc_int %>% filter(group == "tumor"),
        reduction = "umap_integrated", group.by = "cell_type_dtl",
        label = TRUE, repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white")) +
        ggtitle("Tumor Cells - Detailed Cell Types")
    ggsave("results/102.cluster_annotate/final_annotation_tumor_2.png", p2, width = 7, height = 7)
    markers <- list(
        "Steroidogenic" = c("NR5A1", "CYP11A1", "CYP11B2"),
        "ZG" = c("DACH1", "VSNL1"),
        "ZG/ZF" = c("CCN3", "NCAM1"),
        "ZF" = c("CYP11B1", "ABCB1"),
        "ZR" = c("CYB5A", "SULT2A1"),
        "Mast" = c("CPA3", "TPSB2"),
        "CLC" = c("TH", "CHGA", "CHGB", "KIT"),
        "Endo" = c("PECAM1", "EMCN", "PLVAP", "CDH5", "KDR"),
        "PSC" = c("LUM", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP", "PDGFRB", "RGS5"),
        "Fibs" = c("COL1A1", "COL1A2", "COL3A1", "POSTN", "TNC", "NT5E", "THY1", "MCAM"),
        "MSC" = c("ENG", "NES", "STRIP1", "MFAP5", "KLF5", "EFNA5", "EMILIN3"),
        # "PSC_2_bd" = c("CCL19", "APOE", "CXCL2", "CXCL3", "EFEMP1"),
        # "PSC_3_gd" = c("LUM", "PDGFRA", "TAGLN", "MGP", "AIFM2", "S100A4", "FAP"),
        "Adipo" = c("ADIPOQ", "CD34", "FABP4", "ICAM1"),
        "Mac" = c("CD68", "CSF1R", "C1QA", "C1QC"),
        "Tcell" = c("CD4", "TRBC1", "CD3D", "CD3E", "CD3G", "TRBC2", "CD8A"),
        "Bcell" = c("CD19", "MS4A1", "CD79A", "CD79B"),
        "Plasma" = c("CD38", "SDC1", "MZB1", "IGHA1", "IGHG1"),
        "Granul" = c("S100A8", "CCR3", "CD33", "IL5RA", "S100A9", "CSF3R"),
        "LEC" = c("PDPN", "PROX1", "LYVE1", "CCL21"),
        "Eryth" = c("HBB", "GYPA", "SLC4A1", "ALAS2")
    )
    # Select top classical markers for each cell type
    markers_specific <- list(
        "ZG" = c("NR5A1", "DACH1", "CYP11B2"),
        "ZG/ZF" = c("NR5A1", "NOV", "NCAM1"),
        "ZR" = c("NR5A1", "CYB5A", "SULT2A1"),
        "CLC" = c("TH", "CHGA"),
        "Endo" = c("PECAM1", "EMCN"),
        "Fib" = c("COL1A1", "COL3A1", "THY1"),
        "PSC" = c("RGS5", "PDGFRB", "ACTA2"),
        "Adipo" = c("ADIPOQ", "FABP4", "PPARG"),
        "Tcell" = c("CD3D", "CD3E", "TRBC1"),
        "Bcell" = c("CD19", "CD79A", "MS4A1"),
        "Myeloid" = c("ITGAM", "CD33"),
        "Neural" = c("PTPRD", "RBFOX1", "HCN1", "SYT16"),
        "Plasma" = c("CD38", "SDC1", "IGHG1"),
        "LEC" = c("PDPN", "PROX1", "NR2F2")
        # "Eryth" = c("HBB", "GYPA", "SLC4A1")
    )

    toplot <- CalcStats(sc_int,
        features = markers_specific %>% unlist(),
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    gene_groups <- rep(names(markers_specific), lengths(markers_specific)) %>%
        setNames(markers_specific %>% unlist())
    gene_groups <- gene_groups[rownames(toplot)]

    p <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/final_annotation_heatmap_zscore_classical_2.png", p, width = 14, height = 7)

    p_normal <- DotPlot2(sc_int %>% filter(group == "normal"),
        features = markers,
        group.by = "cell_type_dtl",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/final_annotation_normal_dotplot.png", p_normal, width = 14, height = 14)
    p_tumor <- DotPlot2(sc_int %>% filter(group == "tumor"),
        features = markers,
        group.by = "cell_type_dtl",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/final_annotation_tumor_dotplot.png", p_tumor, width = 14, height = 14)
    sc_int
}

save_annotate <- function(sc_final) {
    scCustomize::as.anndata(
        x = sc_final, file_path = "data/104.RNA_velocity", file_name = "anndata.h5ad",
        main_layer = "counts", other_layers = c("data")
    )
    "data/104.RNA_velocity/anndata.h5ad"
}

DEG_cluster <- function(sc_final, cell_type_name) {
    # 聚合表达数据
    # sc_pseudo <- AggregateExpression(sc_final, assays = "RNA", return.seurat = TRUE, group.by = c("group", "dataset", "cell_type_dtl"))
    # Idents(sc_pseudo) <- "cell_type_dtl"


    # 创建目录
    dir.create("results/102.cluster_annotate/cell_types_DEG/normal", showWarnings = FALSE, recursive = TRUE)
    dir.create("results/102.cluster_annotate/cell_types_DEG/tumor", showWarnings = FALSE, recursive = TRUE)

    for (sub_group in c("normal", "tumor")) {
        sc_final_sub <- sc_final %>% filter(group == sub_group)

        # Check if the cell type exists in this group
        if (!(cell_type_name %in% unique(sc_final_sub$cell_type_dtl))) {
            message(paste("Cell type", cell_type_name, "not found in", sub_group, "group. Skipping."))
            next
        }

        sc_pseudo <- AggregateExpression(sc_final_sub, assays = "RNA", return.seurat = TRUE, group.by = c("dataset", "cell_type_dtl"))
        Idents(sc_pseudo) <- "cell_type_dtl"

        deg <- FindMarkers(
            sc_pseudo,
            ident.1 = cell_type_name,
            test.use = "DESeq2",
            min.cells.group = 2
        ) %>%
            rownames_to_column("gene") %>%
            as_tibble() %>%
            filter(p_val_adj < 0.05) %>%
            arrange(desc(abs(avg_log2FC)))

        write_tsv(deg, paste0(
            "results/102.cluster_annotate/cell_types_DEG/", sub_group, "/",
            gsub("/", "_", cell_type_name), "_DEG.tsv"
        ))
    }
    # sc_pseudo <- AggregateExpression(sc_final, assays = "RNA", return.seurat = TRUE, group.by = c("dataset", "cell_type_dtl"))
    # Idents(sc_pseudo) <- "cell_type_dtl"
    # deg <- FindMarkers(
    #     sc_pseudo,
    #     ident.1 = "Fib",
    #     ident.2 = "MSC",
    #     test.use = "DESeq2",
    #     min.cells.group = 2
    # ) %>%
    #     # Add_Pct_Diff() %>%
    #     rownames_to_column("gene") %>%
    #     as_tibble() %>%
    #     filter(p_val_adj < 0.05) %>%
    #     arrange(desc(abs(avg_log2FC)))
    # write_tsv(deg, paste0("results/102.cluster_annotate/cell_types_DEG/PSC_Cap_DEG.tsv"))
    "results/102.cluster_annotate/cell_types_DEG"
}

# 示例调用
compare_neural_clc <- function(sc_int) {
    # Create directory for results
    dir.create("results/102.cluster_annotate/neural_vs_clc", showWarnings = FALSE, recursive = TRUE)

    # Set identity for DEG analysis
    Idents(sc_int) <- "cell_type_dtl"
    # Find differentially expressed genes at single-cell level
    neural_vs_clc_degs <- FindMarkers(
        sc_int,
        ident.1 = "Neural-like",
        ident.2 = "CLC",
        min.pct = 0.2,
        test.use = "MAST"
    ) %>%
        rownames_to_column("gene") %>%
        filter(p_val_adj < 0.05) %>%
        Add_Pct_Diff() %>%
        arrange(desc(abs(avg_log2FC)))

    # Save results to TSV
    write_tsv(
        neural_vs_clc_degs,
        file = "results/102.cluster_annotate/neural_vs_clc/neural_vs_clc_degs_MAST.tsv"
    )

    # Pseudobulk analysis
    # Aggregate expression data by cell type within each sample
    sc_pseudo <- AggregateExpression(
        sc_int,
        assays = "RNA",
        return.seurat = TRUE,
        group.by = c("dataset", "cell_type_dtl")
    )

    # Set identity for pseudobulk DEG analysis
    Idents(sc_pseudo) <- "cell_type_dtl"

    # Find markers using DESeq2 on pseudobulk data
    neural_vs_clc_degs_pseudo <- FindMarkers(
        sc_pseudo,
        ident.1 = "Neural-like",
        ident.2 = "CLC",
        test.use = "DESeq2",
        min.cells.group = 2
    ) %>%
        rownames_to_column("gene") %>%
        filter(p_val_adj < 0.05) %>%
        Add_Pct_Diff() %>%
        arrange(desc(abs(avg_log2FC)))

    # Save pseudobulk results
    write_tsv(
        neural_vs_clc_degs_pseudo,
        file = "results/102.cluster_annotate/neural_vs_clc/neural_vs_clc_degs_pseudobulk.tsv"
    )

    # Define gene markers for visualization
    markers <- list(
        "Chrom ID" = c("DLK1", "GRAMD1B", "AOX1", "CHGA", "CHGB", "TH"),
        "Imm/Lymph" = c("IGHM", "IGHG1", "IGLC2", "CD247", "RHOH", "IL7R", "IKZF1", "IKZF3", "PTPRC", "CELF2"),
        "Cell Cyc" = c("CCND3"),
        "Memb Trans" = c("TMC6", "TMC8", "PPARG")
    )

    # Step 1: Filter dataset to only include Neural-like and CLC cells
    sc_neural_clc <- sc_int %>%
        filter(cell_type_dtl %in% c("Neural-like", "CLC"))

    # Step 2: Prepare gene list and categories
    all_markers <- unlist(markers)
    gene_categories <- data.frame(
        gene = all_markers,
        category = rep(names(markers), times = sapply(markers, length))
    )

    # Step 3: Set up comparison groups
    sc_neural_clc$cell_type_group <- paste(sc_neural_clc$cell_type_dtl, sc_neural_clc$group, sep = "_")
    Idents(sc_neural_clc) <- "cell_type_dtl"

    # Step 4: Get DEGs
    neural_vs_clc <- FindMarkers(
        sc_neural_clc,
        ident.1 = "Neural-like",
        ident.2 = "CLC",
        test.use = "MAST",
        min.cells.group = 2
    ) %>%
        Add_Pct_Diff() %>%
        rownames_to_column("gene") %>%
        as_tibble()

    # Step 5: Filter to markers of interest and prepare data for plotting
    plot_data <- neural_vs_clc %>%
        filter(gene %in% all_markers) %>%
        mutate(
            logFC = avg_log2FC,
            pct_expr_neural = pct.1 * 100,
            pct_expr_clc = pct.2 * 100
        )

    # Add category information
    plot_data <- merge(plot_data, gene_categories, by = "gene")

    # Convert to long format for plotting
    plot_data_long <- pivot_longer(
        plot_data,
        cols = c(pct_expr_neural, pct_expr_clc),
        names_to = "condition",
        values_to = "pct_expr"
    )

    # Clean up condition names
    plot_data_long$condition <- gsub("pct_expr_", "", plot_data_long$condition)

    # Order genes within each category by logFC
    plot_data_long <- plot_data_long %>%
        group_by(category) %>%
        mutate(gene = factor(gene, levels = unique(gene[order(logFC)]))) %>%
        ungroup()

    # Set category order
    plot_data_long$category <- factor(plot_data_long$category, levels = names(markers))

    # Sort data so dots with lower pct_expr appear later (on top)
    plot_data_long <- plot_data_long %>%
        arrange(gene, category, desc(pct_expr))

    # Create the plot
    p <- ggplot(plot_data_long, aes(x = logFC, y = gene, size = pct_expr, color = condition)) +
        geom_point(alpha = 0.8, position = position_identity()) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
        scale_size_continuous(name = "% Expressing", range = c(1, 10)) +
        scale_color_manual(
            name = "Cell Type",
            values = c("neural" = "#fc8d59", "clc" = "#91bfdb"),
            labels = c("neural" = "Neural-like (Tumor)", "clc" = "CLC (Normal)")
        ) +
        facet_grid(category ~ ., scales = "free_y", space = "free_y") +
        labs(x = "Log2 Fold Change (Neural-like/CLC)", y = "", title = "Neural-like vs CLC Cell Comparison") +
        theme_bw() +
        theme(
            strip.background = element_rect(fill = "lightgrey"),
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "right",
            plot.background = element_rect(fill = "white")
        )

    # Save the plot
    dir.create("results/102.cluster_annotate/neural_vs_clc", showWarnings = FALSE, recursive = TRUE)
    ggsave("results/102.cluster_annotate/neural_vs_clc/neural_vs_clc_dotplot.png", p, width = 9, height = 10)

    # Create a heatmap of expression patterns
    neural_clc_genes <- intersect(all_markers, rownames(sc_neural_clc))
    toplot <- CalcStats(sc_neural_clc,
        features = neural_clc_genes,
        group.by = "cell_type_dtl",
        method = "zscore",
        order = "value"
    )

    # Add gene categories as row annotations
    gene_groups <- setNames(gene_categories$category, gene_categories$gene)
    gene_groups <- gene_groups[rownames(toplot)]

    p_heat <- Heatmap(t(toplot), lab_fill = "zscore", facet_col = gene_groups) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/102.cluster_annotate/neural_vs_clc/neural_vs_clc_heatmap.png", p_heat, width = 10, height = 7)

    "results/102.cluster_annotate/neural_vs_clc"
}
