load_sc <- function() {
    samples <- c("N1", "N2", "N4", "T1", "T2", "T4")
    sc_raw_list <- imap(samples, ~ Read10X(paste0("data/102.expr_data/", .x), gene.column = 1))
    names(sc_raw_list) <- samples
    sc_raw_list <- list(
        sc_raw_n = createLiger(sc_raw_list[c("N1", "N2", "N4")]),
        sc_raw_t = createLiger(sc_raw_list[c("T1", "T2", "T4")])
    )
    sc_raw_list
}

process_sc <- function(sc_raw_list) {
    sc_list <- sc_raw_list %>% map(~ ligerToSeurat(.x))
    sc_list <- map(sc_list, ~ .x %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap_unintegrated"))

    p_list <- map(sc_list, ~ DimPlot2(.x, reduction = "umap_unintegrated", group.by = c("orig.ident")) +
        theme(plot.background = element_rect(fill = "white")) +
        labs(title = paste("UMAP un-integrated", .x$orig.ident[1])))
    # use patchwork to combine multiple plots
    p <- patchwork::wrap_plots(p_list, ncol = 2)
    ggsave("results/101.load_data/umap_unintegrated.png", p, width = 14, height = 7)

    sc_int_list <- map(sc_list, ~ IntegrateLayers(.x, dims = 1:30, orig.reduction = "pca", method = HarmonyIntegration, new.reduction = "harmony"))
    sc_int_list[[1]][["RNA"]] <- JoinLayers(sc_int_list[[1]][["RNA"]])
    sc_int_list[[2]][["RNA"]] <- JoinLayers(sc_int_list[[2]][["RNA"]])
    sc_int_list <- map(sc_int_list, ~ RunUMAP(.x, dims = 1:30, reduction = "harmony", reduction.name = "umap_integrated"))
    p <- map(sc_int_list, ~ DimPlot2(.x, reduction = "umap_integrated", group.by = c("orig.ident")) +
        theme(plot.background = element_rect(fill = "white")) +
        labs(title = paste("UMAP integrated", .x$orig.ident[1]))) %>% patchwork::wrap_plots(ncol = 2)
    ggsave("results/101.load_data/umap_integrated.png", p, width = 14, height = 7)
    sc_int <- merge(sc_int_list[[1]], sc_int_list[[2]], merge.dr = TRUE, merge.data = TRUE)
    sc_int <- sc_int %>%
        RunUMAP(dims = 1:30, reduction = "harmony", reduction.name = "umap_integrated")
    sc_int <- sc_int %>%
        # if orig.ident is starting with "N", then it is normal, otherwise tumor
        mutate(group = ifelse(grepl("^N", orig.ident), "normal", "tumor"))
    p <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = c("orig.ident", "group")) +
        theme(plot.background = element_rect(fill = "white"))
    ggsave("results/101.load_data/umap_integrated_merge.png", p, width = 14, height = 7)
    sc_int
}
