DEG_annotation <- function(sc_annotate) {
    sc_pseudo <- AggregateExpression(sc_annotate, assays = "RNA", return.seurat = T, group.by = c("group", "dataset", "cell_type_dtl"))
    sc_pseudo$cell_type_group <- paste(sc_pseudo$cell_type_dtl, sc_pseudo$group, sep = "_")
    Idents(sc_pseudo) <- "cell_type_group"
    # DEG_adipo <- FindMarkers(sc_pseudo,
    #     ident.1 = "Adipo_tumor",
    #     ident.2 = "Adipo_normal",
    #     test.use = "DESeq2",
    #     min.cells.group = 2
    # ) %>%
    #     Add_Pct_Diff() %>%
    #     rownames_to_column("gene") %>%
    #     as_tibble() %>%
    #     filter(p_val_adj < 0.05) %>%
    #     arrange(desc(abs(avg_log2FC)))
    # DEG_Fibro <- FindMarkers(sc_pseudo,
    #     ident.1 = "Fib_tumor",
    #     ident.2 = "Fib_normal",
    #     test.use = "DESeq2",
    #     min.cells.group = 2
    # ) %>%
    #     Add_Pct_Diff() %>%
    #     rownames_to_column("gene") %>%
    #     as_tibble() %>%
    #     filter(p_val_adj < 0.05) %>%
    #     arrange(desc(abs(avg_log2FC)))
    # write_tsv(DEG_adipo, "results/103.DEG/DEG_adipo.tsv")
    # write_tsv(DEG_Fibro, "results/103.DEG/DEG_Fibro.tsv")
    cell_types <- unique(sc_pseudo$cell_type_dtl)
    for (cell_type in cell_types) {
        ident1 <- paste0(cell_type, "_tumor")
        ident2 <- paste0(cell_type, "_normal")
        if (!(ident1 %in% Idents(sc_pseudo) && ident2 %in% Idents(sc_pseudo))) {
            next
        }
        deg <- FindMarkers(
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
            arrange(desc(abs(avg_log2FC)))
        write_tsv(deg, paste0("results/103.DEG/DEG_", gsub("/", "_", cell_type), ".tsv"))
    }
    "results/103.DEG"
}
