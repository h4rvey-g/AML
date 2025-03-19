sub_cluster_tcell <- function(sc_final) {
    # Extract T cells
    sc_int <- sc_final %>%
        filter(cell_type_dtl == "Tcell")

    # Create directory for results
    dir.create("results/108.Tcell", showWarnings = FALSE, recursive = TRUE)

    # Split by group
    sc_int_normal <- sc_int %>% filter(group == "normal")
    sc_int_tumor <- sc_int %>% filter(group == "tumor")

    # Process normal group
    sc_int_normal <- sc_int_normal %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.9, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster") %>%
        RunUMAP(
            dims = 1:30,
            reduction = "scvi",
            reduction.name = "umap_integrated",
            min.dist = 0.05,
            n.neighbors = 50,
            spread = 0.5
        )

    # Process tumor group
    sc_int_tumor <- sc_int_tumor %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.9, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster") %>%
        RunUMAP(
            dims = 1:30,
            reduction = "scvi",
            reduction.name = "umap_integrated",
            min.dist = 0.05,
            n.neighbors = 50,
            spread = 0.5
        )

    # Visualize normal group clustering
    p1 <- DimPlot2(sc_int_normal,
        reduction = "umap_integrated",
        group.by = "immune_subcluster",
        label = TRUE
    ) +
        ggtitle("Normal") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_normal.png", p1, width = 7, height = 7)

    # Visualize tumor group clustering
    p2 <- DimPlot2(sc_int_tumor,
        reduction = "umap_integrated",
        group.by = "immune_subcluster",
        label = TRUE
    ) +
        ggtitle("Tumor") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_tumor.png", p2, width = 7, height = 7)

    # Return all three objects as a list
    list(
        normal = sc_int_normal,
        tumor = sc_int_tumor,
        combined = sc_int
    )
}
TCellAnnotation_CellTypist <- function(sc_tcell_clust_list) {
    # Create directory for results
    dir.create("results/108.Tcell/celltypist", showWarnings = FALSE, recursive = TRUE)

    # Extract normal and tumor objects
    sc_tcell_clust_normal <- sc_tcell_clust_list$normal
    sc_tcell_clust_tumor <- sc_tcell_clust_list$tumor

    # Define function to run celltypist
    run_celltypist <- function(seurat_obj, output_prefix) {
        # Convert to AnnData
        adata_path <- paste0("results/108.Tcell/celltypist/", output_prefix, "_data.h5ad")
        as.anndata(seurat_obj,
            file_path = dirname(adata_path),
            file_name = basename(adata_path),
            main_layer = "counts",
            other_layers = c("data")
        )

        # Run celltypist using reticulate
        reticulate::use_condaenv("r-reticulate", required = TRUE)
        celltypist <- reticulate::import("celltypist")

        # Load the model
        model <- celltypist$models$Model$load(model = "Immune_All_Low.pkl")

        # Run prediction
        adata <- reticulate::import("anndata")$read_h5ad(adata_path)
        predictions <- celltypist$annotate(adata, model = model, majority_voting = TRUE)

        # Extract predictions
        pred_df <- reticulate::py_to_r(predictions$predicted_labels)
        pred_df <- as.data.frame(pred_df)

        # Save predictions
        write_csv(pred_df, paste0("results/108.Tcell/celltypist/", output_prefix, "_predictions.csv"))

        # Add annotations to Seurat object
        seurat_obj$celltypist_celltype <- pred_df$majority_voting[match(colnames(seurat_obj), rownames(pred_df))]

        # Create UMAP plot with celltypist annotations
        p <- DimPlot2(seurat_obj,
            reduction = "umap_integrated",
            group.by = "celltypist_celltype",
            label = TRUE,
            label.size = 3
        ) +
            ggtitle(paste0(output_prefix, " - CellTypist Annotations")) +
            theme(plot.background = element_rect(fill = "white"))

        ggsave(paste0("results/108.Tcell/celltypist/", output_prefix, "_umap.png"), p, width = 12, height = 10)

        return(seurat_obj)
    }

    # Run celltypist on normal and tumor data
    sc_tcell_clust_normal <- run_celltypist(sc_tcell_clust_normal, "normal")
    sc_tcell_clust_tumor <- run_celltypist(sc_tcell_clust_tumor, "tumor")

    # Compare celltypist annotations with immune_subcluster
    compare_annotations <- function(seurat_obj, output_prefix) {
        comparison <- table(
            CellTypist = seurat_obj$celltypist_celltype,
            Cluster = seurat_obj$immune_subcluster
        )

        # Save comparison table
        write.csv(
            comparison,
            paste0("results/108.Tcell/celltypist/", output_prefix, "_comparison.csv")
        )

        # Create heatmap
        pheatmap::pheatmap(
            log1p(comparison),
            display_numbers = TRUE,
            filename = paste0("results/108.Tcell/celltypist/", output_prefix, "_comparison_heatmap.png"),
            width = 12,
            height = 10
        )
    }

    compare_annotations(sc_tcell_clust_normal, "normal")
    compare_annotations(sc_tcell_clust_tumor, "tumor")

    # Return updated list
    list(
        normal = sc_tcell_clust_normal,
        tumor = sc_tcell_clust_tumor,
        combined = sc_tcell_clust_list$combined
    )
}
Tcell_annotation_analysis <- function(sc_tcell_clust_list, sc_final) {
    sc_tcell_clust_normal <- sc_tcell_clust_list$normal
    sc_tcell_clust_tumor <- sc_tcell_clust_list$tumor
    markers <- list(
        # Core T Cell Markers
        "T_Cells" = c("CD3D", "CD3E", "CD3G"),
        "CD4_T_Cells" = c("CD4"),
        "CD8_T_Cells" = c("CD8A", "CD8B"),

        # Specialized T Cell Subsets
        "Tregs" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "TNFRSF9", "TNFRSF4", "IL10", "ENTPD1"),
        "gdT_Cells" = c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2"),
        "NK_Cells" = c("NKG7", "GNLY", "KLRD1", "NCAM1", "FCGR3A", "NCR1", "NCR3", "KIR2DL1", "KIR2DL3", "KIR3DL1"),
        "MAIT_NKT" = c("ZBTB16", "SLC4A10", "KLRB1", "TRAV1-2", "DPP4"),

        # Naive and Memory Markers
        "Naive_T" = c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "SKAP1", "THEMIS", "PTPRC", "S1PR1", "KLF2"),
        "Memory_T" = c("CD44", "IL7R", "CD45RO", "PRKCQ", "STAT4", "CD27", "CD28"),

        # Tissue-Resident Memory
        "Trm" = c("ITGAE", "CD69", "CXCR6", "ITGA1"),

        # Migration and Tissue Homing
        "Homing" = c("CCR9", "CXCR4", "CCR2", "CCR10", "SELL", "ITGB7"),

        # Effector Functions and Cytotoxicity
        "Effector" = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "IFNG", "TNF", "FASLG", "LAMP1", "KLRG1", "CX3CR1", "IL32", "TNFSF10", "EOMES"),

        # Exhaustion Markers
        "Exhausted" = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "ENTPD1", "TOX", "PRDM1", "BTLA"),

        # Activation Markers
        "Activation" = c("ICOS", "HLA-DRA", "HLA-DRB1", "CD69", "MKI67", "PCNA", "TOP2A", "CDK1"),

        # T Helper Subtypes and Key Transcription Factors
        "Th1" = c("IFNG", "CXCR3", "TBX21", "CCR5", "GZMK", "STAT1"),
        "Th2" = c("IL4", "IL5", "IL13", "CCR4", "GATA3", "CCR3", "STAT6"),
        "Th17" = c("IL17A", "IL17F", "IL22", "CCR6", "RORC", "IL23R", "IL1R1", "CCL20", "BATF"),
        "Tfh" = c("BCL6", "CXCR5", "ICOS", "MAF"),
        "Th22" = c("IL22", "AHR"),

        # T Cell Fate Decision Markers
        "Fate_Regulators" = c("ID2", "ID3", "TOX", "BATF", "MAF")
    )


    # 计算并绘制热图
    # Calculate and plot heatmap for normal group
    toplot_normal <- CalcStats(sc_tcell_clust_normal,
        features = markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster"
    )

    gene_groups_normal <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_normal)]

    p2_normal <- Heatmap(toplot_normal, lab_fill = "zscore", facet_row = gene_groups_normal) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_normal.png", p2_normal, width = 8, height = 15)
    p2_normal_dot <- DotPlot2(sc_tcell_clust_normal,
        features = markers,
        group.by = "immune_subcluster"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_normal.png", p2_normal_dot, width = 8, height = 15)
    write_csv(
        toplot_normal %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/108.Tcell/immune_subcluster_expression_normal.csv"
    )

    # Calculate and plot heatmap for tumor group
    toplot_tumor <- CalcStats(sc_tcell_clust_tumor,
        features = markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster"
    )

    gene_groups_tumor <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_tumor)]

    p2_tumor <- Heatmap(toplot_tumor, lab_fill = "zscore", facet_row = gene_groups_tumor) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_tumor.png", p2_tumor, width = 8, height = 15)
    p2_tumor_dot <- DotPlot2(sc_tcell_clust_tumor,
        features = markers,
        group.by = "immune_subcluster"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_tumor.png", p2_tumor_dot, width = 8, height = 15)
    write_csv(
        toplot_tumor %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 2)),
        "results/108.Tcell/immune_subcluster_expression_tumor.csv"
    )

    CD4_clusters_tumor <- c(1, 2, 4)
    CD8_clusters_tumor <- c(3, 5, 6)
    DN_clusters_tumor <- c(7)
    sc_tcell_CD4_tumor <- sc_tcell_clust_tumor %>%
        filter(immune_subcluster %in% CD4_clusters_tumor) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 1.5, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_CD4")
    sc_tcell_CD8_tumor <- sc_tcell_clust_tumor %>%
        filter(immune_subcluster %in% CD8_clusters_tumor) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 1.8, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_CD8")
    sc_tcell_DN_tumor <- sc_tcell_clust_tumor %>%
        filter(immune_subcluster %in% DN_clusters_tumor) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 1, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_DN")
    CD4_clusters_normal <- c(1)
    CD8_clusters_normal <- c(2, 3)
    DN_clusters_normal <- c(4, 5)
    sc_tcell_CD4_normal <- sc_tcell_clust_normal %>%
        filter(immune_subcluster %in% CD4_clusters_normal) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 0.9, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_CD4")
    sc_tcell_CD8_normal <- sc_tcell_clust_normal %>%
        filter(immune_subcluster %in% CD8_clusters_normal) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 1, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_CD8")
    sc_tcell_DN_normal <- sc_tcell_clust_normal %>%
        filter(immune_subcluster %in% DN_clusters_normal) %>%
        FindNeighbors(dims = 1:30, reduction = "scvi") %>%
        FindClusters(resolution = 1, algorithm = 4, method = "igraph", cluster.name = "immune_subcluster_DN")

    # remove the first 3 vectors in markers
    cd4_markers <- list(
        # T cell states
        "activated" = c("TNFRSF18", "TNFRSF4"),
        "Tc" = c("GZMA", "GZMK"),
        "proliferation" = c("MKI67"),
        "Tcm" = c("S100A10", "LTB"),
        "Tem" = c("CXCR4", "KLRB1"),
        "Temra" = c("NKG7", "GNLY"),
        "Tex" = c("SRGN", "DUSP4"),
        "Tfh" = c("PDCD1", "ICOS"),
        "Th1" = c("CCL5", "CXCR3"),
        "Th17" = c("IL7R", "CCR6"),
        "Th22" = c("IL22", "AHR"),
        "Ttsg" = c("ISG15", "STAT1"),
        "Tn" = c("SELL", "TCF7", "CCR7"),
        "Treg" = c("FOXP3", "IL2RA"),
        "Trm" = c("ITGA4", "ITGB2"),
        "Tstr" = c("HSP90AA1", "HSPA1B")
    )


    toplot_CD4_tumor <- CalcStats(sc_tcell_CD4_tumor,
        features = cd4_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_CD4"
    )
    gene_groups_CD4_tumor <- rep(names(cd4_markers), lengths(cd4_markers)) %>%
        setNames(cd4_markers %>% unlist()) %>%
        .[rownames(toplot_CD4_tumor)]
    p2_CD4_tumor <- Heatmap(toplot_CD4_tumor, lab_fill = "zscore", facet_row = gene_groups_CD4_tumor) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_CD4_tumor.png", p2_CD4_tumor, width = 8, height = 15)
    p2_CD4_tumor_dot <- DotPlot2(sc_tcell_CD4_tumor,
        features = cd4_markers,
        group.by = "immune_subcluster_CD4"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_CD4_tumor.png", p2_CD4_tumor_dot, width = 8, height = 15)
    write_csv(
        toplot_CD4_tumor %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_CD4_tumor.csv"
    )

    toplot_CD4_normal <- CalcStats(sc_tcell_CD4_normal,
        features = cd4_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_CD4"
    )
    gene_groups_CD4_normal <- rep(names(cd4_markers), lengths(cd4_markers)) %>%
        setNames(cd4_markers %>% unlist()) %>%
        .[rownames(toplot_CD4_normal)]
    p2_CD4_normal <- Heatmap(toplot_CD4_normal, lab_fill = "zscore", facet_row = gene_groups_CD4_normal) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_CD4_normal.png", p2_CD4_normal, width = 8, height = 15)
    p2_CD4_normal_dot <- DotPlot2(sc_tcell_CD4_normal,
        features = cd4_markers,
        group.by = "immune_subcluster_CD4"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_CD4_normal.png", p2_CD4_normal_dot, width = 8, height = 15)
    write_csv(
        toplot_CD4_normal %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_CD4_normal.csv"
    )

    cd8_markers <- list(
        "NKT-like" = c("NKG7", "KLRD1", "FCGR3A"),
        "proliferation" = c("MKI67", "HMGB2", "HMGB1"),
        "quiescence" = c("RUNX1"),
        "senescence" = c("ATM", "TSC2"),
        "Tc" = c("CX3CR1", "GZMH"),
        "Tcm" = c("CREM", "IL7R", "CCR7"),
        "Tem" = c("CST7", "EOMES", "GZMK", "GZMA"),
        "Temra" = c("GNLY", "FGFBP2"),
        "Tex" = c("LAG3", "TIGIT", "HAVCR2", "PDCD1"),
        "Tn" = c("LEF1", "SELL", "TCF7", "CCR7"),
        # "Tn_quiescence" = c("BTG1", "KLF2", "BTG2"),
        # "Tn_IFN_response" = c("STAT1", "IFITM1", "IRF1", "IFITM3"),
        "Trm" = c("ITGA1", "ITGB1"),
        "Tstr" = c("HSPA5", "HSP90AA1", "HSP90AB1")
    )

    toplot_CD8_tumor <- CalcStats(sc_tcell_CD8_tumor,
        features = cd8_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_CD8"
    )
    gene_groups_CD8_tumor <- rep(names(cd8_markers), lengths(cd8_markers)) %>%
        setNames(cd8_markers %>% unlist()) %>%
        .[rownames(toplot_CD8_tumor)]
    p2_CD8_tumor <- Heatmap(toplot_CD8_tumor, lab_fill = "zscore", facet_row = gene_groups_CD8_tumor) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_CD8_tumor.png", p2_CD8_tumor, width = 8, height = 15)
    p2_CD8_tumor_dot <- DotPlot2(sc_tcell_CD8_tumor,
        features = cd8_markers,
        group.by = "immune_subcluster_CD8"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_CD8_tumor.png", p2_CD8_tumor_dot, width = 8, height = 15)
    write_csv(
        toplot_CD8_tumor %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_CD8_tumor.csv"
    )
    toplot_CD8_normal <- CalcStats(sc_tcell_CD8_normal,
        features = cd8_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_CD8"
    )
    gene_groups_CD8_normal <- rep(names(cd8_markers), lengths(cd8_markers)) %>%
        setNames(cd8_markers %>% unlist()) %>%
        .[rownames(toplot_CD8_normal)]
    p2_CD8_normal <- Heatmap(toplot_CD8_normal, lab_fill = "zscore", facet_row = gene_groups_CD8_normal) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_CD8_normal.png", p2_CD8_normal, width = 8, height = 15)
    p2_CD8_normal_dot <- DotPlot2(sc_tcell_CD8_normal,
        features = cd8_markers,
        group.by = "immune_subcluster_CD8"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_CD8_normal.png", p2_CD8_normal_dot, width = 8, height = 15)
    write_csv(
        toplot_CD8_normal %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_CD8_normal.csv"
    )

    dn_markers <- list(
        # TCR Types
        "alphabeta_TCR" = c("TRAV", "TRBV", "TRAJ", "TRBJ"),
        "gammadelta_TCR" = c("TRGC1", "TRGC2", "TRDC", "TRDV1", "TRDV2", "TRGV9", "TRDV3"),

        # DN T Cell Developmental Stages
        "DN1" = c("CD44", "KIT", "IL7R", "FLT3", "CD24", "BCL11A", "GATA3"),
        "DN2" = c("CD44", "KIT", "IL7R", "CD24", "BCL11B", "GATA3", "NOTCH1", "HES1"),
        "DN3" = c("CD24", "PTCRA", "RAG1", "RAG2", "NOTCH1", "HES1", "DTX1", "BCL11B"),
        "DN4" = c("CD24", "ITM2A", "LEF1", "ID3", "THEMIS"),

        # Unconventional T Cell Subtypes
        "NKT_Cells" = c("ZBTB16", "CD160", "KLRB1", "CD244", "KLRC1", "GZMA", "VA24-JA18", "Vβ11"),
        "MAIT_Cells" = c("SLC4A10", "TRAV1-2", "RORC", "ZBTB16", "KLRB1", "DPP4", "IL18R1", "CXCR6", "CCR6"),

        # iNKT Developmental Stages
        "iNKT_Stage0" = c("LEF1", "KLF2", "GATA3", "IL4", "IL7R", "CCR7"),
        "iNKT_Stage1" = c("GATA3", "IL4", "IFNG", "ZBTB16"),
        "iNKT_Stage2" = c("TBXT", "IFNG", "IL4", "RORC"),

        # Functional DN Subsets
        "DN_Regulatory" = c("FOXP3", "IKZF2", "EGR2", "IL10", "CTLA4", "TNFRSF18", "ENTPD1", "TGFB1", "LAG3", "NRP1"),
        "DN_Cytotoxic" = c("GZMB", "PRF1", "FASLG", "TNFSF10", "IFNG", "LAMP1", "NKG7"),
        "DN_Helper" = c("IL2", "IL4", "IL17A", "IL21", "CD40LG"),

        # Tissue-Specific DN T Cells
        "Intestinal_DN" = c("ITGAE", "IL22", "IL17A", "RORC", "AHR", "CCR9", "ITGA4", "ITGB7"),
        "Skin_DN" = c("CLA", "CCR4", "CCR8", "CCR10", "IL22", "IL17A"),
        "Liver_DN" = c("CXCR3", "CXCR6", "CCR5", "S1PR1", "KLRB1", "CD69"),

        # Pathological DN T Cells
        "Autoimmune_DN" = c("IL17A", "IL17F", "IL23R", "CCR6", "RORC", "STAT3", "IL1R1"),
        "Senescent_DN" = c("CDKN1A", "CDKN2A", "GLB1", "TP53BP1", "MAPK14", "CD57"),

        # DN T Cell Activation/Differentiation
        "Naive_DN" = c("SELL", "CCR7", "CD27", "IL7R", "BACH2", "SATB1", "LEF1"),
        "Memory_DN" = c("S100A4", "ANXA1", "CD101", "IL7R", "CD44"),
        "Activated_DN" = c("CD69", "HLA-DRA", "HLA-DRB1", "MKI67", "TNFRSF9", "IL2RA"),

        # DN-Specific Transcription Factors
        "DN_TF" = c("ID3", "SOX4", "SOX13", "ETV5", "NFIL3", "PLZF", "RUNX1", "RUNX3", "TOX2", "MYB")
    )

    toplot_DN_tumor <- CalcStats(sc_tcell_DN_tumor,
        features = dn_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_DN"
    )
    gene_groups_DN_tumor <- rep(names(dn_markers), lengths(dn_markers)) %>%
        setNames(dn_markers %>% unlist()) %>%
        .[rownames(toplot_DN_tumor)]
    p2_DN_tumor <- Heatmap(toplot_DN_tumor, lab_fill = "zscore", facet_row = gene_groups_DN_tumor) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_DN_tumor.png", p2_DN_tumor, width = 8, height = 15)
    p2_DN_tumor_dot <- DotPlot2(sc_tcell_DN_tumor,
        features = dn_markers,
        group.by = "immune_subcluster_DN"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_DN_tumor.png", p2_DN_tumor_dot, width = 8, height = 15)
    write_csv(
        toplot_DN_tumor %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_DN_tumor.csv"
    )
    toplot_DN_normal <- CalcStats(sc_tcell_DN_normal,
        features = dn_markers %>% unlist(),
        method = "zscore", order = "value",
        group.by = "immune_subcluster_DN"
    )
    gene_groups_DN_normal <- rep(names(dn_markers), lengths(dn_markers)) %>%
        setNames(dn_markers %>% unlist()) %>%
        .[rownames(toplot_DN_normal)]
    p2_DN_normal <- Heatmap(toplot_DN_normal, lab_fill = "zscore", facet_row = gene_groups_DN_normal) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_heatmap_DN_normal.png", p2_DN_normal, width = 8, height = 15)
    p2_DN_normal_dot <- DotPlot2(sc_tcell_DN_normal,
        features = dn_markers,
        group.by = "immune_subcluster_DN"
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/sub_cluster_immune_dotplot_DN_normal.png", p2_DN_normal_dot, width = 8, height = 15)
    write_csv(
        toplot_DN_normal %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            mutate(across(where(is.numeric), round, 4)),
        "results/108.Tcell/immune_subcluster_expression_DN_normal.csv"
    )

    # Add immune_subcluster_CD4 back to the tumor Seurat object
    sc_tcell_clust_tumor$immune_subcluster_CD4 <- NA
    sc_tcell_clust_tumor$immune_subcluster_CD4[colnames(sc_tcell_CD4_tumor)] <- sc_tcell_CD4_tumor$immune_subcluster_CD4

    # Add immune_subcluster_CD8 back to the tumor Seurat object
    sc_tcell_clust_tumor$immune_subcluster_CD8 <- NA
    sc_tcell_clust_tumor$immune_subcluster_CD8[colnames(sc_tcell_CD8_tumor)] <- sc_tcell_CD8_tumor$immune_subcluster_CD8

    # Add immune_subcluster_DN back to the tumor Seurat object
    sc_tcell_clust_tumor$immune_subcluster_DN <- NA
    sc_tcell_clust_tumor$immune_subcluster_DN[colnames(sc_tcell_DN_tumor)] <- sc_tcell_DN_tumor$immune_subcluster_DN

    # Add immune_subcluster_CD4 back to the normal Seurat object
    sc_tcell_clust_normal$immune_subcluster_CD4 <- NA
    sc_tcell_clust_normal$immune_subcluster_CD4[colnames(sc_tcell_CD4_normal)] <- sc_tcell_CD4_normal$immune_subcluster_CD4

    # Add immune_subcluster_CD8 back to the normal Seurat object
    sc_tcell_clust_normal$immune_subcluster_CD8 <- NA
    sc_tcell_clust_normal$immune_subcluster_CD8[colnames(sc_tcell_CD8_normal)] <- sc_tcell_CD8_normal$immune_subcluster_CD8

    # Add immune_subcluster_DN back to the normal Seurat object
    sc_tcell_clust_normal$immune_subcluster_DN <- NA
    sc_tcell_clust_normal$immune_subcluster_DN[colnames(sc_tcell_DN_normal)] <- sc_tcell_DN_normal$immune_subcluster_DN

    list(
        "normal" = sc_tcell_clust_normal,
        "tumor" = sc_tcell_clust_tumor
    )
}


Tcell_annotate <- function(sc_tcell_clust_list, sc_final) {
    sc_tcell_clust_normal <- sc_tcell_clust_list$normal
    sc_tcell_clust_tumor <- sc_tcell_clust_list$tumor

    # Define cluster maps for tumor and normal groups
    # Define cluster maps for tumor and normal groups
    cluster_map_tumor_CD4 <- c(
        "1" = "CD4_Th22",
        "2" = "CD4_Tn",
        "3" = "CD4_Th1/Th17",
        "4" = "CD4_Tn",
        "5" = "CD4_Tn",
        "6" = "CD4_Th1/Th17",
        "7" = "CD4_Th1/Th17",
        "8" = "CD4_Treg",
        "9" = "CD4_Tex"
    )
    cluster_map_tumor_CD8 <- c(
        "1" = "CD8_Temra",
        "2" = "CD8_Tex",
        "3" = "CD8_Tex",
        "4" = "CD8_Tn",
        "5" = "CD8_Tn",
        "6" = "CD8_Temra",
        "7" = "CD8_Tn",
        "8" = "CD8_Tn",
        "9" = "CD8_Tn"
    )
    cluster_map_tumor_DN <- c(
        "1" = "DN_gdT_2",
        "2" = "DN_gdT_1",
        "3" = "DN_gdT_2"
    )

    cluster_map_normal_CD4 <- c(
        "1" = "CD4_Tex",
        "2" = "CD4_Th1/Th17",
        "3" = "CD4_proliferation"
    )
    cluster_map_normal_CD8 <- c(
        "1" = "CD8_quiescence/gdT_1",
        "2" = "CD8_proliferation",
        "3" = "CD8_Tem/gdT_2",
        "4" = "CD8_Temra"
    )
    cluster_map_normal_DN <- c(
        "1" = "DN_gdT_2",
        "2" = "DN_gdT_1"
    )

    # Add cell types to sc_tcell_clust_tumor based on the subcluster assignments
    sc_tcell_clust_tumor$T_cell_type_CD4 <- NA
    sc_tcell_clust_tumor$T_cell_type_CD8 <- NA
    sc_tcell_clust_tumor$T_cell_type_DN <- NA

    # Map CD4 cells in tumor
    cd4_cells_tumor <- !is.na(sc_tcell_clust_tumor$immune_subcluster_CD4)
    sc_tcell_clust_tumor$T_cell_type_CD4[cd4_cells_tumor] <- cluster_map_tumor_CD4[as.character(sc_tcell_clust_tumor$immune_subcluster_CD4[cd4_cells_tumor])]

    # Map CD8 cells in tumor
    cd8_cells_tumor <- !is.na(sc_tcell_clust_tumor$immune_subcluster_CD8)
    sc_tcell_clust_tumor$T_cell_type_CD8[cd8_cells_tumor] <- cluster_map_tumor_CD8[as.character(sc_tcell_clust_tumor$immune_subcluster_CD8[cd8_cells_tumor])]

    # Map DN cells in tumor
    dn_cells_tumor <- !is.na(sc_tcell_clust_tumor$immune_subcluster_DN)
    sc_tcell_clust_tumor$T_cell_type_DN[dn_cells_tumor] <- cluster_map_tumor_DN[as.character(sc_tcell_clust_tumor$immune_subcluster_DN[dn_cells_tumor])]

    # Add cell types to sc_tcell_clust_normal based on the subcluster assignments
    sc_tcell_clust_normal$T_cell_type_CD4 <- NA
    sc_tcell_clust_normal$T_cell_type_CD8 <- NA
    sc_tcell_clust_normal$T_cell_type_DN <- NA

    # Map CD4 cells in normal
    cd4_cells_normal <- !is.na(sc_tcell_clust_normal$immune_subcluster_CD4)
    sc_tcell_clust_normal$T_cell_type_CD4[cd4_cells_normal] <- cluster_map_normal_CD4[as.character(sc_tcell_clust_normal$immune_subcluster_CD4[cd4_cells_normal])]

    # Map CD8 cells in normal
    cd8_cells_normal <- !is.na(sc_tcell_clust_normal$immune_subcluster_CD8)
    sc_tcell_clust_normal$T_cell_type_CD8[cd8_cells_normal] <- cluster_map_normal_CD8[as.character(sc_tcell_clust_normal$immune_subcluster_CD8[cd8_cells_normal])]

    # Map DN cells in normal
    dn_cells_normal <- !is.na(sc_tcell_clust_normal$immune_subcluster_DN)
    sc_tcell_clust_normal$T_cell_type_DN[dn_cells_normal] <- cluster_map_normal_DN[as.character(sc_tcell_clust_normal$immune_subcluster_DN[dn_cells_normal])]

    # Add a combined cell_type_dtl field by taking the first non-NA value from the T_cell_type fields
    sc_tcell_clust_tumor <- sc_tcell_clust_tumor %>%
        mutate(cell_type_dtl = case_when(
            !is.na(T_cell_type_CD4) ~ T_cell_type_CD4,
            !is.na(T_cell_type_CD8) ~ T_cell_type_CD8,
            !is.na(T_cell_type_DN) ~ T_cell_type_DN,
            TRUE ~ "Undefined"
        ))

    sc_tcell_clust_normal <- sc_tcell_clust_normal %>%
        mutate(cell_type_dtl = case_when(
            !is.na(T_cell_type_CD4) ~ T_cell_type_CD4,
            !is.na(T_cell_type_CD8) ~ T_cell_type_CD8,
            !is.na(T_cell_type_DN) ~ T_cell_type_DN,
            TRUE ~ "Undefined"
        ))
    # Merge the tumor and normal objects for the combined object
    sc_tcell_clust <- sc_final %>% filter(cell_type == "Tcell")
    sc_tcell_clust@meta.data$cell_type_dtl <- NA
    sc_tcell_clust@meta.data[colnames(sc_tcell_clust_normal), "cell_type_dtl"] <- sc_tcell_clust_normal$cell_type_dtl
    sc_tcell_clust@meta.data[colnames(sc_tcell_clust_tumor), "cell_type_dtl"] <- sc_tcell_clust_tumor$cell_type_dtl


    # Run UMAP
    sc_tcell_clust <- RunUMAP(sc_tcell_clust,
        dims = 1:30, reduction = "scvi", reduction.name = "umap_tcell",
        min.dist = 0.05,
        spread = 0.3
    )

    # Visualize the new annotations
    p <- DimPlot2(sc_tcell_clust,
        reduction = "umap_tcell",
        group.by = "cell_type_dtl",
        label = TRUE,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))

    # Save the plot
    ggsave("results/108.Tcell/tcell_annotated.png", p, width = 9, height = 6)

    # Plot by group
    p_split <- DimPlot2(sc_tcell_clust,
        reduction = "umap_tcell",
        group.by = "cell_type_dtl",
        split.by = "group",
        label = TRUE,
        repel = TRUE
    ) +
        theme(plot.background = element_rect(fill = "white"))

    ggsave("results/108.Tcell/tcell_annotated_by_group.png", p_split, width = 15, height = 6)

    markers <- list(
        # Core T cell subtype markers
        "CD4/8" = c("CD4", "CD8A", "CD8B"),
        "CD4_Tc" = c("GZMA", "GZMK"),
        "CD4_Tcm" = c("S100A10", "LTB"),
        # "Tem" = c("CXCR4", "KLRB1"),
        "CD4_Temra" = c("NKG7", "GNLY"),
        "CD4_Tex" = c("SRGN", "DUSP4"),
        "Tn" = c("LEF1", "SELL", "TCF7", "CCR7"),
        "Tfh" = c("PDCD1", "ICOS"),
        "CD4_Th1" = c("CCL5", "CXCR3"),
        "CD4_Th17" = c("IL7R", "CCR6"),
        "CD4_Th22" = c("IL22", "AHR"),
        # "Treg" = c("FOXP3", "IL2RA"),
        # "Trm" = c("ITGA4", "ITGB2"),
        # "NKT-like" = c("NKG7", "KLRD1", "FCGR3A"),
        "CD8_proliferation" = c("MKI67", "HMGB2", "HMGB1"),
        # "quiescence" = c("RUNX1"),
        # "senescence" = c("ATM", "TSC2"),
        "CD8_Tc" = c("CX3CR1", "GZMH"),
        # "Tcm" = c("CREM", "IL7R", "CCR7"),
        "CD8_Tem" = c("CST7", "EOMES", "GZMK", "GZMA"),
        "CD8_Temra" = c("GNLY", "FGFBP2"),
        "CD8_Tex" = c("LAG3", "TIGIT", "HAVCR2", "PDCD1"),
        "CD8_Tn_quiescence" = c("BTG1", "KLF2", "BTG2"),
        "CD8_Tn_IFN_response" = c("STAT1", "IFITM1", "IRF1", "IFITM3"),
        "DN_gdT" = c("TRGC1", "TRGC2", "TRDC", "TRDV1", "TRDV2")
    )

    # Calculate and plot heatmap
    toplot <- CalcStats(sc_tcell_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )

    # toplot <- toplot[, c("CD4_Treg/Th17", "CD4_TN/Th17", "CD4_TN", "CD8_FOXP3_Treg", "CD8_Cyto_gdT", "CD8_Cyto_Ext_gdT", "Resting_T", "CTLA4_T")]
    # # Get unique cell types from both groups
    # cell_types <- unique(sc_tcell_clust$cell_type_dtl)

    gene_groups <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))

    p2 <- Heatmap(toplot %>% t(), lab_fill = "zscore", facet_col = gene_groups, text.size = 5) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_heatmap.png", p2, width = 18, height = 8)

    p <- DotPlot2(sc_tcell_clust,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_dotplot.png", p, width = 15, height = 8)

    markers <- list(
        # "activated" = c("TNFRSF18", "TNFRSF4"),
        # "Tc" = c("GZMA", "GZMK"),
        "proliferation" = c("MKI67"),
        # "Tcm" = c("S100A10", "LTB"),
        # "Tem" = c("CXCR4", "KLRB1"),
        # "Temra" = c("NKG7", "GNLY"),
        "Tex" = c("SRGN", "DUSP4"),
        # "Tfh" = c("PDCD1", "ICOS"),
        "Th1" = c("CCL5", "CXCR3"),
        "Th17" = c("IL7R", "CCR6"),
        "Th22" = c("IL22", "AHR"),
        # "Ttsg" = c("ISG15", "STAT1"),
        "Tn" = c("SELL", "TCF7", "CCR7"),
        "Treg" = c("FOXP3", "IL2RA")
        # "Trm" = c("ITGA4", "ITGB2"),
        # "Tstr" = c("HSP90AA1", "HSPA1B")
    )
    # calculate for CD4, CD8, and DN separately
    toplot_CD4 <- CalcStats(
        sc_tcell_clust %>%
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
    ggsave("results/108.Tcell/tcell_annotated_CD4_heatmap.png", p1, width = 15, height = 8)
    p <- DotPlot2(
        sc_tcell_clust %>%
            filter(grepl("^CD4", cell_type_dtl)),
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_CD4_dotplot.png", p, width = 15, height = 8)


    markers <- list(
        "proliferation" = c("MKI67", "HMGB2", "HMGB1"),
        "quiescence" = c("RUNX1"),
        # "senescence" = c("ATM", "TSC2"),
        # "Tc" = c("CX3CR1", "GZMH"),
        # "Tcm" = c("CREM", "IL7R", "CCR7"),
        "Tem" = c("CST7", "EOMES", "GZMK", "GZMA"),
        "Temra" = c("GNLY", "FGFBP2"),
        "Tex" = c("LAG3", "TIGIT", "HAVCR2", "PDCD1"),
        "Tn" = c("LEF1", "SELL", "TCF7", "CCR7")
        # "Tn_quiescence" = c("BTG1", "KLF2", "BTG2"),
        # "Tn_IFN_response" = c("STAT1", "IFITM1", "IRF1", "IFITM3")
    )
    toplot_CD8 <- CalcStats(
        sc_tcell_clust %>%
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
    ggsave("results/108.Tcell/tcell_annotated_CD8_heatmap.png", p2, width = 15, height = 8)
    p <- DotPlot2(
        sc_tcell_clust %>%
            filter(grepl("^CD8", cell_type_dtl)),
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_CD8_dotplot.png", p, width = 15, height = 8)

    markers <- list(
        "CD4/8" = c("CD4", "CD8A", "CD8B"),
        "gdT_1" = c("TRDV1"),
        "gdT_2" = c("TRDV2"),
        "gdT_constant" = c("TRGC1", "TRGC2", "TRDC")
    )
    toplot_DN <- CalcStats(
        sc_tcell_clust,
        features = markers %>% unlist(),
        method = "zscore", order = "value", group.by = "cell_type_dtl"
    )
    gene_groups_DN <- rep(names(markers), lengths(markers)) %>%
        setNames(markers %>% unlist()) %>%
        .[rownames(toplot_DN)] %>%
        factor(levels = unique(rep(names(markers), lengths(markers))))
    p3 <- Heatmap(toplot_DN %>% t(), lab_fill = "zscore", facet_col = gene_groups_DN) +
        labs(title = "CD4+, CD8+, DN T Cell Subtypes") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_DN_heatmap.png", p3, width = 15, height = 8)
    p <- DotPlot2(sc_tcell_clust,
        features = markers,
        group.by = "cell_type_dtl",
        split.by = "group",
        show_grid = FALSE
    ) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/tcell_annotated_DN_dotplot.png", p, width = 15, height = 8)

    p <- p1 +
        p2 +
        p3 +
        plot_layout(ncol = 1)
    ggsave("results/108.Tcell/tcell_annotated_heatmap_combined.png", p, width = 15, height = 25)
    # Return the updated list with annotated objects
    sc_tcell_clust
}

test_tcell_composition <- function(sc_tcell) {
    # Create directory
    dir.create("results/108.Tcell/distribution", showWarnings = FALSE, recursive = TRUE)

    # Create count data frame for sccomp with proper dataset column
    # tcell_counts <- sc_tcell@meta.data %>%
    #     select(cell_type_dtl, dataset, group) %>%
    #     mutate(
    #         count = 1
    #     ) %>%
    #     group_by(cell_type_dtl, dataset, group) %>%
    #     summarize(count = n(), .groups = "drop")

    # Run composition test with updated parameters
    sc_tcell <- sc_tcell %>%
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
    composition_test <- sc_tcell %>%
        sccomp_estimate(
            formula_composition = ~group,
            # formula_variability = ~1,
            .sample = sample,
            .cell_group = cell_group,
            bimodal_mean_variability_association = TRUE,
            # .abundance = count, # Updated from .count to .abundance
            cores = 20,
            verbose = TRUE
        ) %>%
        sccomp_remove_outliers(cores = 31, verbose = FALSE) %>%
        sccomp_test()


    # Generate plots
    p1 <- composition_test %>%
        sccomp_boxplot(factor = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/boxplot.png", p1, width = 10, height = 8)

    p2 <- composition_test %>%
        plot_1D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/effect_size.png", p2, width = 8, height = 6)

    p3 <- composition_test %>%
        plot_2D_intervals() +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/distribution/abundance_variability.png", p3, width = 8, height = 6)

    # Save fold changes
    fold_changes <- composition_test %>%
        sccomp_proportional_fold_change(
            formula_composition = ~group,
            from = "normal",
            to = "tumor"
        ) %>%
        select(cell_group, statement)

    write_tsv(fold_changes, "results/108.Tcell/distribution/fold_changes.tsv")

    # Save test results
    write_tsv(
        composition_test %>%
            select(-count_data),
        "results/108.Tcell/distribution/test_results.tsv"
    )

    composition_test
}

plot_tcell_distribution <- function(sc_tcell, t_composition_test) {
    dir.create("results/108.Tcell/distribution", showWarnings = FALSE, recursive = TRUE)

    # Calculate distribution percentages
    df1 <- ClusterDistrBar(origin = sc_tcell$dataset, cluster = sc_tcell$cell_type_dtl, plot = FALSE) %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        pivot_longer(-cluster, names_to = "origin", values_to = "percentage")

    # Calculate cell counts
    cell_counts <- sc_tcell@meta.data %>%
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

    t_composition_test_res <- t_composition_test %>%
        filter(factor == "group", c_FDR < 0.05) %>%
        # convert c_FDR to stars
        mutate(stars = case_when(
            c_FDR < 0.001 ~ "***",
            c_FDR < 0.01 ~ "**",
            c_FDR < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        select(cell_group, stars)

    # Add stars to cluster labels based on significance
    df1 <- df1 %>%
        mutate(cluster = map_chr(str_extract(cluster, ".*(?=\\n)"), ~ {
            stars <- t_composition_test_res$stars[t_composition_test_res$cell_group == .x]
            count_info <- sprintf("(N=%d,T=%d)", cell_counts$N_count[cell_counts$cell_type_dtl == .x], cell_counts$T_count[cell_counts$cell_type_dtl == .x])
            if (length(stars) > 0 && stars != "") {
                paste0(.x, "\n", count_info, "\n", stars)
            } else {
                paste0(.x, "\n", count_info, "\n")
            }
        }))

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

    tidyplots::save_plot(p, "results/108.Tcell/distribution/tcell_distribution.png")
    write_tsv(
        df1 %>%
            mutate(cluster = str_extract(cluster, ".*(?=\\n)")),
        "results/108.Tcell/distribution/tcell_distribution.tsv"
    )

    return("results/108.Tcell/distribution")
}

run_tcell_GSEA <- function(sc_tcell) {
    options(max.print = 12, spe = "human")
    dir.create("results/108.Tcell/GSEA", showWarnings = FALSE, recursive = TRUE)

    # Run GO analysis
    parents <- GeneSetAnalysisGO() %>% names()
    for (parent in parents) {
        dir.create(paste0("results/108.Tcell/GSEA/", parent), showWarnings = FALSE, recursive = TRUE)
        cat("Running GSEA for", parent, "\n")
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 20, parent = parent, n.min = 3)
        matr <- sc_tcell@misc$AUCell$GO[[parent]]
        matr <- RenameGO(matr, add_id = FALSE)

        # Overall heatmap for all subclusters
        toplot <- CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10)
        p <- Heatmap(
            toplot,
            lab_fill = "zscore", text.size = 5
        )
        ggsave(paste0("results/108.Tcell/GSEA/", parent, "/Heatmap_GO_", parent, ".png"), p, width = 14, height = nrow(toplot) * 0.4)
        write_tsv(
            toplot %>% as.data.frame() %>% rownames_to_column("pathway"),
            paste0("results/108.Tcell/GSEA/", parent, "/Heatmap_GO_", parent, ".tsv")
        )

        # Compare tumor vs normal within each subcluster
        subclusters <- unique(sc_tcell$cell_type_dtl)
        for (subcluster in subclusters) {
            sc_temp <- sc_tcell %>%
                filter(cell_type_dtl == subcluster)
            Idents(sc_temp) <- "group"

            if (length(unique(sc_temp$group)) == 2) {
                p <- WaterfallPlot(
                    matr,
                    f = sc_temp$group,
                    ident.1 = "tumor",
                    ident.2 = "normal",
                    top.n = 20,
                    color = "p",
                    length = "logFC",
                    title = paste0(parent, ": ", subcluster)
                )
                ggsave(paste0(
                    "results/108.Tcell/GSEA/", parent, "/Waterfall_GO_",
                    parent, "_", gsub("/", "_", subcluster), "_tumor_vs_normal.png"
                ), p, width = 15, height = 15)
            }
        }
    }

    # Run hallmark analysis
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = hall50$human, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for hallmark pathways
    dir.create("results/108.Tcell/GSEA/hallmark50", showWarnings = FALSE, recursive = TRUE)
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10),
        lab_fill = "zscore", text.size = 5
    )
    ggsave("results/108.Tcell/GSEA/hallmark50/Heatmap_hallmark50.png", p, width = 14, height = 14)
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA/hallmark50/Heatmap_hallmark50.tsv"
    )

    # Run for each GO root category
    for (root in c("BP", "MF", "CC")) {
        cat("Processing", root, "...\n")
        sc_tcell <- GeneSetAnalysisGO(sc_tcell, nCores = 20, root = root, n.min = 3)
        matr <- sc_tcell@misc$AUCell$GO[[root]]

        top30_pathways <- data.frame()
        for (cell_type in unique(sc_tcell$cell_type_dtl)) {
            stats <- CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore") %>%
                as.data.frame() %>%
                dplyr::select(all_of(cell_type)) %>%
                arrange(desc(.[[1]])) %>%
                head(30) %>%
                rownames_to_column("pathway")

            top30_pathways <- bind_rows(top30_pathways, stats)
        }

        top30_pathways <- top30_pathways %>%
            pull(pathway) %>%
            unique()

        matr_filt <- matr[top30_pathways, ]
        matr_filt <- RenameGO(matr_filt, add_id = FALSE)
        rownames(matr_filt) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr_filt)))
        toplot <- CalcStats(matr_filt, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 30)

        dir.create(file.path("results/108.Tcell/GSEA/GO", root), showWarnings = FALSE, recursive = TRUE)

        p <- Heatmap(
            CalcStats(matr_filt, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 10),
            lab_fill = "zscore", text.size = 5
        )
        ggsave(file.path("results/108.Tcell/GSEA/GO", root, "Heatmap_top10.png"), p, width = 14, height = 14)
        write_tsv(
            toplot %>%
                as.data.frame() %>% rownames_to_column("pathway"),
            file.path("results/108.Tcell/GSEA/GO", root, "Heatmap_top10.tsv")
        )
    }

    return("results/108.Tcell/GSEA")
}

tcell_DEG_tumor_vs_normal <- function(sc_tcell, t_cell_type) {
    sc_pseudo <- AggregateExpression(sc_tcell, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"

    dir.create("results/108.Tcell/DEG_tumor_vs_normal", showWarnings = FALSE, recursive = TRUE)

    # DESeq2 analysis
    ident1 <- paste0(t_cell_type, "_tumor")
    ident2 <- paste0(t_cell_type, "_normal")
    # if (!(ident1 %in% Idents(sc_pseudo) && ident2 %in% Idents(sc_pseudo))) {
    #     next
    # }
    deg <- tryCatch(
        {
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
                filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
                arrange(desc(abs(avg_log2FC))) %>%
                filter(!str_starts(gene, "CYP")) %>%
                filter(!str_starts(gene, "MT-"))
        },
        error = function(e) {
            message("Error in FindMarkers for ", t_cell_type, ": ", e$message)
            return(NULL)
        }
    )

    if (!is.null(deg)) {
        write_tsv(
            deg %>%
                mutate(across(where(is.numeric), ~ ifelse(abs(.) < 0.001, signif(., 3), round(., 3)))),
            paste0("results/108.Tcell/DEG_tumor_vs_normal/", gsub("/", "_", t_cell_type), "_DEG.tsv")
        )
    }

    # MAST analysis with parallel processing
    sc_tcell$cell_type_group <- paste(sc_tcell$cell_type_dtl, sc_tcell$group, sep = "_")
    Idents(sc_tcell) <- "cell_type_group"

    ident1 <- paste0(t_cell_type, "_tumor")
    ident2 <- paste0(t_cell_type, "_normal")
    if ((ident1 %in% Idents(sc_tcell) && ident2 %in% Idents(sc_tcell))) {
        deg <- FindMarkers(
            sc_tcell,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            min.cells.group = 2
        ) %>%
            scCustomize::Add_Pct_Diff() %>%
            rownames_to_column("gene") %>%
            as_tibble() %>%
            filter(p_val_adj < 0.05) %>%
            filter((avg_log2FC > 0 & pct.1 > 0.2) | (avg_log2FC < 0 & pct.2 > 0.2)) %>%
            arrange(desc(abs(avg_log2FC))) %>%
            filter(!str_starts(gene, "CYP")) %>%
            filter(!str_starts(gene, "MT-"))

        write_tsv(
            deg %>%
                mutate(across(where(is.numeric), ~ ifelse(abs(.) < 0.001, signif(., 3), round(., 3)))),
            paste0("results/108.Tcell/DEG_tumor_vs_normal/", gsub("/", "_", t_cell_type), "_DEG_MAST.tsv")
        )
    }

    return("results/108.Tcell/DEG_tumor_vs_normal")
}

analyze_CD8_Tex_TRM_correlations <- function(sc_tcell) {
    # 创建结果目录
    dir.create("results/108.Tcell/correlation", showWarnings = FALSE, recursive = TRUE)

    # 提取CD8_Tex_TRM细胞
    sc_tex <- subset(sc_tcell, cell_type_dtl == "CD8_Tex_TRM")

    # 定义要分析的基因列表
    gene_sets <- list(
        "Adenosine_Receptors" = c("ADORA1", "ADORA2A", "ADORA2B", "ADORA3"),
        "Cytotoxicity" = c("IFNG", "GZMB", "PRF1", "TNF", "GNLY"),
        "Treg_Markers" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3"),
        "Exhaustion" = c("PDCD1", "HAVCR2", "TIGIT"),
        "Tissue_Residency" = c("ITGAE", "CD69", "CXCR6")
    )

    # 初始化结果列表
    correlation_results <- list()

    # Only analyze tumor group
    sc_group <- sc_tex %>%
        filter(group == "tumor")

    # 提取感兴趣基因的表达矩阵
    genes_to_analyze <- unique(c("ENTPD1", unlist(gene_sets)))
    expr_matrix <- GetAssayData(sc_group, slot = "counts", assay = "RNA") %>%
        as.matrix() %>%
        t() %>%
        sclink_norm(scale.factor = 1e6, filter.genes = FALSE, gene.names = genes_to_analyze)

    # 使用scLink计算相关性
    corr <- sclink_net(expr = expr_matrix, ncores = 30)
    correlation_results <- corr$cor

    # Process correlation matrix
    process_correlation_matrix <- function(corr_matrix) {
        # Only keep correlations with ENTPD1
        entpd1_corr <- corr_matrix["ENTPD1", ]

        # Convert to data frame
        df <- data.frame(
            Gene = names(entpd1_corr),
            Correlation = entpd1_corr
        ) %>%
            # Remove self-correlation
            filter(Gene != "ENTPD1")

        # Add gene set information
        df <- df %>%
            mutate(Gene_Set = case_when(
                Gene %in% gene_sets$Adenosine_Receptors ~ "Adenosine_Receptors",
                Gene %in% gene_sets$Cytotoxicity ~ "Cytotoxicity",
                Gene %in% gene_sets$Treg_Markers ~ "Treg_Markers",
                Gene %in% gene_sets$Exhaustion ~ "Exhaustion",
                Gene %in% gene_sets$Tissue_Residency ~ "Tissue_Residency"
            ))

        return(df)
    }

    all_correlations <- process_correlation_matrix(correlation_results)

    # Save detailed results
    write_csv(
        all_correlations,
        "results/108.Tcell/correlation/CD8_Tex_TRM_ENTPD1_correlations_tumor.csv"
    )

    # Create visualization
    p2 <- ggplot(
        all_correlations,
        aes(
            x = reorder(Gene, Correlation),
            y = Correlation
        )
    ) +
        geom_bar(stat = "identity", aes(fill = Correlation > 0)) +
        scale_fill_manual(values = c("#4169E1", "#FF6B6B")) +
        geom_text(aes(label = sprintf("%.2f", Correlation)),
            vjust = ifelse(all_correlations$Correlation > 0, -0.5, 1.5)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white")
        ) +
        labs(
            x = "Gene",
            y = "scLink Correlation with ENTPD1",
            title = "Correlation of ENTPD1 with Selected Genes in CD8_Tex_TRM (Tumor)"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_grid(~Gene_Set, scales = "free_x", space = "free")
    ggsave("results/108.Tcell/correlation/CD8_Tex_TRM_ENTPD1_correlations_tumor.png",
        p2,
        width = 15, height = 7
    )
}
run_GSEA_Tex_TRM <- function(sc_tcell) {
    options(max.print = 12, spe = "human")
    dir.create("results/108.Tcell/GSEA_Tex_TRM", showWarnings = FALSE, recursive = TRUE)

    # Create custom gene set of adenosine
    adenosine_genes <- SearchDatabase(
        item = "adenosine",
        type = "SetName",
        return = "genelist",
        database = "GO",
        n.min = 3
    )

    # Run GSEA for adenosine gene set
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = adenosine_genes, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for adenosine pathways
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20),
        lab_fill = "zscore"
    ) +
        theme(text = element_text(size = 14))
    ggsave("results/108.Tcell/GSEA_Tex_TRM/Heatmap_adenosine.png", p, width = 14, height = 14)
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA_Tex_TRM/Heatmap_adenosine.tsv"
    )
    # Create custom gene set of hypoxia
    hypoxia_genes <- SearchDatabase(
        item = "hypoxia",
        type = "SetName",
        return = "genelist",
        database = "GO",
        n.min = 3
    )

    # Run GSEA for hypoxia gene set
    sc_tcell <- GeneSetAnalysis(sc_tcell, genesets = hypoxia_genes, nCores = 30)
    matr <- sc_tcell@misc$AUCell$genesets
    matr <- RenameGO(matr, add_id = FALSE)
    rownames(matr) <- stringr::str_to_title(stringr::str_to_lower(rownames(matr)))

    # Overall heatmap for hypoxia pathways
    p <- Heatmap(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20),
        lab_fill = "zscore"
    ) +
        theme(text = element_text(size = 14))
    ggsave("results/108.Tcell/GSEA_Tex_TRM/Heatmap_hypoxia.png", p, width = 14, height = 14)
    # write tsv files
    write_tsv(
        CalcStats(matr, f = sc_tcell$cell_type_dtl, method = "zscore", order = "p", n = 20) %>%
            as.data.frame() %>% rownames_to_column("pathway"),
        "results/108.Tcell/GSEA_Tex_TRM/Heatmap_hypoxia.tsv"
    )
    "results/108.Tcell/GSEA_Tex_TRM"
}

run_tcell_trajectory <- function(sc_tcell) {
    dir.create("results/108.Tcell/trajectory", showWarnings = FALSE, recursive = TRUE)
    tar_load(sc_tcell)
    library(reticulate)
    py_run_string("")
    library(SeuratExtend)
    library(tidyseurat)
    sc_tcell <- sc_tcell %>% filter(group == "tumor")

    # 1. scVelo Analysis
    # Convert Seurat object to AnnData for scVelo
    dir.create("data/108.Tcell/trajectory", showWarnings = FALSE, recursive = TRUE)
    Seu2Adata(sc_tcell, save.adata = "data/108.Tcell/trajectory/tcell.h5ad", conda_env = "base")
    adata_path <- "data/108.Tcell/trajectory/tcell.h5ad"
    seu <- sc_tcell
    seu <- scVelo.SeuratToAnndata(
        sc_tcell,
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
        basis = "umap_tcell_cell_embeddings",
        save = "results/108.Tcell/trajectory/velocity_umap.png",
        figsize = c(7, 6),
        conda_env = "base"
    )

    # 2. Palantir Analysis
    # Run diffusion map
    set.seed(42)
    seu <- Palantir.RunDM(seu,
        reduction = "scvi",
        conda_env = "base",
        n_components = 30
    )
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    cells_to_remove_1 <- CellSelector(p)
    # Remove selected cells
    seu <- subset(seu, cells = setdiff(colnames(seu), cells_to_remove_1))
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    # Create plot again after removing cells
    cells_to_remove_2 <- CellSelector(p)
    seu <- subset(seu, cells = setdiff(colnames(seu), cells_to_remove_2))
    p <- DimPlot2(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/trajectory/plantir_dm.png", p, width = 7, height = 6)
    # Get cells from seu and update adata to only include those cells
    cells <- colnames(seu)
    py_run_string(sprintf("adata = adata[adata.obs.index.isin(%s)]", paste0("[", paste0("'", cells, "'", collapse = ","), "]")))
    # Add ms dimension reduction to AnnData object
    adata.AddDR(seu, dr = "ms", scv.graph = FALSE, conda_env = "base")

    # Plot velocity using ms coordinates
    # scVelo.Plot(
    #     color = "cell_type_dtl",
    #     basis = "ms",
    #     save = "results/108.Tcell/trajectory/velocity_ms.png",
    #     figsize = c(7, 6),
    #     conda_env = "base"
    # )

    # Select start cell (CD4_TN cluster)
    p <- DimPlot(seu, reduction = "ms", group.by = "cell_type_dtl", label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    start_cell <- CellSelector(p)
    start_cell <- colnames(seu)[which(seu$cell_type_dtl == "CD4_TN")[1]]
    # fate1_cell <- CellSelector(p)
    # fate2_cell <- CellSelector(p)
    # fate3_cell <- CellSelector(p)

    # Calculate pseudotime
    seu <- Palantir.Pseudotime(seu,
        start_cell = start_cell, conda_env = "base", n_jobs = 10
        # terminal_states = c("fate1" = fate1_cell, "fate2" = fate2_cell, "fate3" = fate3_cell)
    )

    # Get pseudotime data
    ps <- seu@misc$Palantir$Pseudotime
    colnames(ps)[3:4] <- c("fate1", "fate2")
    seu@meta.data[, colnames(ps)] <- ps

    # Plot pseudotime and entropy
    p1 <- DimPlot2(seu,
        features = colnames(ps),
        reduction = "ms",
        cols = list(Entropy = "D"),
        theme = NoAxes()
    ) + theme(plot.background = element_rect(fill = "white"))
    ggsave("results/108.Tcell/trajectory/palantir_pseudotime.png", p1, width = 15, height = 12)

    # # 3. MAGIC for gene expression smoothing
    # sc_tcell <- Palantir.Magic(sc_tcell)
    # sc_tcell <- NormalizeData(sc_tcell)

    # # Plot gene expression before and after MAGIC
    # key_genes <- c("CD4", "CD8A", "FOXP3", "PDCD1")
    # p2 <- DimPlot2(sc_tcell,
    #     features = c(key_genes, paste0("magic_", key_genes)),
    #     cols = "A",
    #     theme = NoAxes()
    # )
    # ggsave("results/108.Tcell/trajectory/magic_expression.png", p2, width = 12, height = 12)

    # 4. CellRank Analysis
    # Add pseudotime to adata
    # Add pseudotime to adata
    adata.AddMetadata(seu, col = colnames(ps), conda_env = "base")

    # Run CellRank
    Cellrank.Compute(time_key = "Pseudotime", conda_env = "base")

    # Generate CellRank plots
    Cellrank.Plot(
        color = "cell_type_dtl",
        basis = "ms",
        save = "results/108.Tcell/trajectory/cellrank_ms.png",
        conda_env = "base"
    )

    # 5. Gene Expression Dynamics
    # Generate trend curves for key genes
    # key_genes <- c("CD4", "CD8A", "FOXP3", "PDCD1", "CTLA4", "LAG3")
    # p3 <- GeneTrendCurve.Palantir(
    #     seu,
    #     pseudotime.data = ps,
    #     features = key_genes,
    #     point = FALSE,
    #     se = TRUE,
    #     conda_env = "base"
    # )
    # ggsave("results/108.Tcell/trajectory/gene_trends.png", p3, width = 10, height = 8)

    # fate1_genes <- c(
    #     "PRF1",        # Key cytolytic protein for target cell killing
    #     "GZMK",        # Granzyme K, serine protease involved in cytotoxicity
    #     "GNLY",        # Granulysin, antimicrobial protein in cytotoxic cells
    #     "KLRC3",       # NK cell receptor family member
    #     "KLRC4-KLRK1", # NK cell activating receptor
    #     "CD244",       # NK cell receptor that regulates cytotoxicity
    #     "RUNX3",       # Master regulator of CD8+ T cell differentiation
    #     "TOX",         # Critical for effector and exhausted CD8+ T cell development
    #     "CCL5"         # Chemokine expressed by activated T cells
    # )
    # fate2_genes <- c(
    #     "TSHZ2",      # Associated with regulatory T cell function
    #     "SERINC5",    # Membrane protein expressed in CD4+ T cells
    #     "PRKCA",      # Protein kinase C alpha, involved in T cell activation pathways
    #     "RASGRF2",    # Involved in T cell receptor signaling
    #     "TIAM1"      # Regulator of T cell migration
    # )
    # Define gene sets that highlight fate differences
    tf_genes <- c(
        "RORC", # Th17 transcription factor
        "FOXP3", # Regulatory T cell master regulator
        "TBX21", # T-bet, Th1/cytotoxic transcription factor
        "GATA3" # Th2 cells transcription factor
    )

    effector_genes <- c(
        "IL17F", # Th17 cytokine
        "IFNG", # Th1/cytotoxic cytokine
        "PRF1", # Perforin, cytotoxic effector
        "GZMA", # Granzyme A, cytotoxic effector
        "GZMB" # Granzyme B, cytotoxic effector
    )

    surface_genes <- c(
        "CD4", # Helper T cells marker
        "CD8A", # Cytotoxic T cells marker
        "TRDV1", # γδ T cell receptor marker
        "TRDV2", # γδ T cell receptor marker
        "ENTPD1", # CD39, Treg suppressive function
        "IL2RA" # CD25, Treg marker
    )

    immunoregulatory_genes <- c(
        "CTLA4", # Checkpoint inhibitor, high in Tregs
        "PDCD1", # PD-1, exhaustion marker
        "TIGIT", # Inhibitory receptor
        "NKG7", # Natural killer cell granule protein 7
        "GNLY", # Granulysin, cytotoxic molecule
        "CCR6" # Chemokine receptor for Th17 cell migration
    )
    # Define core exhaustion gene set
    core_exhaustion_genes <- c(
        "PDCD1", # PD-1
        "CTLA4", # CTLA-4
        "HAVCR2", # TIM-3
        "LAG3", # Lymphocyte activation gene 3
        "TIGIT", # T cell immunoreceptor with Ig and ITIM domains
        "ENTPD1" # CD39
    )

    # Check if genes exist in the dataset
    existing_genes <- intersect(core_exhaustion_genes, rownames(seu))

    # For fate 1
    p_exh_fate1 <- GeneTrendHeatmap.Palantir(
        seu,
        features = existing_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/exhaustion_heatmap_fate1.png",
        p_exh_fate1,
        width = 10, height = 6
    )

    # For fate 2
    p_exh_fate2 <- GeneTrendHeatmap.Palantir(
        seu,
        features = existing_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/exhaustion_heatmap_fate2.png",
        p_exh_fate2,
        width = 10, height = 6
    )

    # Generate trend curves
    p_exh_trends <- GeneTrendCurve.Palantir(
        seu,
        pseudotime.data = ps,
        features = existing_genes,
        point = FALSE,
        se = TRUE,
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/exhaustion_genes_trends.png",
        p_exh_trends,
        width = 12, height = 8
    )

    # Generate trend heatmaps for each gene set on both fates
    # Transcription factors
    p_tf_fate1 <- GeneTrendHeatmap.Palantir(
        seu,
        features = tf_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/tf_genes_heatmap_fate1.png", p_tf_fate1, width = 10, height = 6)

    p_tf_fate2 <- GeneTrendHeatmap.Palantir(
        seu,
        features = tf_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/tf_genes_heatmap_fate2.png", p_tf_fate2, width = 10, height = 6)

    # Effector molecules
    p_eff_fate1 <- GeneTrendHeatmap.Palantir(
        seu,
        features = effector_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/effector_genes_heatmap_fate1.png", p_eff_fate1, width = 10, height = 6)

    p_eff_fate2 <- GeneTrendHeatmap.Palantir(
        seu,
        features = effector_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/effector_genes_heatmap_fate2.png", p_eff_fate2, width = 10, height = 6)

    # Surface markers
    p_surface_fate1 <- GeneTrendHeatmap.Palantir(
        seu,
        features = surface_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/surface_genes_heatmap_fate1.png", p_surface_fate1, width = 10, height = 7)

    p_surface_fate2 <- GeneTrendHeatmap.Palantir(
        seu,
        features = surface_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/surface_genes_heatmap_fate2.png", p_surface_fate2, width = 10, height = 7)

    # Immunoregulatory genes
    p_imreg_fate1 <- GeneTrendHeatmap.Palantir(
        seu,
        features = immunoregulatory_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/immunoregulatory_genes_heatmap_fate1.png", p_imreg_fate1, width = 10, height = 7)

    p_imreg_fate2 <- GeneTrendHeatmap.Palantir(
        seu,
        features = immunoregulatory_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/immunoregulatory_genes_heatmap_fate2.png", p_imreg_fate2, width = 10, height = 7)

    # Generate trend curves for key discriminating genes
    key_discriminating_genes <- c("RORC", "TBX21", "ENTPD1", "TRDV1", "IL17F", "GZMB", "CD4", "CD8A")
    p_trends <- GeneTrendCurve.Palantir(
        seu,
        pseudotime.data = ps,
        features = key_discriminating_genes,
        point = FALSE,
        se = TRUE,
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/key_discriminating_genes_trends.png", p_trends, width = 12, height = 8)
    # Generate heatmap for fate1 genes
    p4 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate1_genes,
        pseudotime.data = ps,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/gene_heatmap_fate1.png", p4, width = 12, height = 8)


    p5 <- GeneTrendHeatmap.Palantir(
        seu,
        features = fate2_genes,
        pseudotime.data = ps,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/108.Tcell/trajectory/gene_heatmap_fate2.png", p5, width = 12, height = 8)

    # 6. Slingshot Analysis
    sc_tcell <- RunSlingshot(sc_tcell,
        group.by = "cell_type_dtl",
        start.clus = "CD4_TN"
    )

    # Add Slingshot pseudotime to metadata
    sling <- sc_tcell@misc$slingshot$PCA$SlingPseudotime
    sc_tcell@meta.data[, colnames(sling)] <- as.data.frame(sling)

    # Plot Slingshot results
    p5 <- DimPlot2(sc_tcell,
        features = colnames(sling),
        cols = "C",
        theme = NoAxes()
    )
    ggsave("results/108.Tcell/trajectory/slingshot_pseudotime.png", p5, width = 12, height = 5)

    return(sc_tcell)
}
