load_sc <- function() {
    samples <- c("N1", "N2", "N4", "T1", "T2", "T4")
    sc_raw_list <- imap(samples, ~ Read10X(paste0("data/103.self_workflow/", .x, "/output/filter_matrix"), gene.column = 1))
    names(sc_raw_list) <- samples
    sc_raw <- createLiger(sc_raw_list)
    sc_raw
}

process_sc <- function(sc_raw) {
    sc <- sc_raw %>% ligerToSeurat()
    sc <- sc %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap_unintegrated")
    sc <- sc %>%
        mutate(group = ifelse(grepl("^N", orig.ident), "normal", "tumor"))

    p <- DimPlot2(sc, reduction = "umap_unintegrated", group.by = c("orig.ident"), split.by = "group") +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/101.load_data/umap_unintegrated.png", p, width = 14, height = 7)

    sc_int <- IntegrateLayers(sc, dims = 1:30, orig.reduction = "pca", method = HarmonyIntegration, new.reduction = "harmony")
    sc_int[["RNA"]] <- JoinLayers(sc_int[["RNA"]])
    sc_int <- sc_int %>% RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap_integrated")
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = c("orig.ident", "group")) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/101.load_data/umap_integrated.png", p, width = 14, height = 7)
    sc_int
}
