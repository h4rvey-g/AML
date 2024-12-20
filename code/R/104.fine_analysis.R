run_palantir <- function(sc_final) {
    library(SeuratExtend)
    library(Seurat)
    library(tidyseurat)
    tar_load(sc_final)
    sc_final <- sc_final %>%
        filter(cell_type %in% c("Fib-Cap-Adipo"))
    sc_final <- Palantir.RunDM(sc_final, conda_env = "base", reduction = "harmony")
    # create a var, combine cell_type_dtl and group
    sc_final <- sc_final %>%
        mutate(cell_type_group = paste(cell_type_dtl, group, sep = "_"))
    p <- DimPlot2(sc_final, reduction = "ms", group.by = c("cell_type_group", "cell_type_dtl"), label = TRUE, repel = TRUE) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/106.fine_analysis/palantir.png", p, width = 12, height = 6)
    p <- DimPlot(sc_final, reduction = "ms", group.by = c("cell_type_dtl"))
    cells <- CellSelector(p)
    sc_final <- Palantir.Pseudotime(sc_final, start_cell = cells)
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

    # Top 15 downregulated genes contributing to mCAF differentiation
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
    p <- GeneTrendHeatmap.Palantir(
        sc_final,
        features = genes_up,
        pseudotime.data = ps,
        magic = FALSE,
        lineage = "fate2",
        conda_env = "base"
    )
    ggsave("results/106.fine_analysis/gene_trend_up.png", p, width = 12, height = 6)
    p <- GeneTrendHeatmap.Palantir(
        sc_final,
        features = genes_down,
        pseudotime.data = ps,
        magic = FALSE,
        lineage = "fate1",
        conda_env = "base"
    )
    ggsave("results/106.fine_analysis/gene_trend_down.png", p, width = 12, height = 6)
    ggsave("results/106.fine_analysis/gene_trend.png", p, width = 12, height = 6)

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
