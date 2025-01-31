load_sc <- function() {
  samples <- c("N1", "N2", "N4", "N7", "T1", "T2", "T4", "T7")
  sc_raw_list <- imap(samples, ~ Read10X(paste0("data/103.self_workflow/", .x, "/output/filter_matrix"), gene.column = 1))
  names(sc_raw_list) <- samples
  sc_raw <- createLiger(sc_raw_list)
  sc_raw
}

process_sc_preview <- function(sc_raw) {
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

find_doublet <- function(sc_int) {
  # Convert Seurat to SingleCellExperiment
  sce <- as.SingleCellExperiment(sc_int)

  # Run scDblFinder
  set.seed(42)
  sce <- scDblFinder(sce, samples = "orig.ident", dbr = 0.03)

  # Transfer results back to Seurat object
  sc_int$scDblFinder.score <- sce$scDblFinder.score
  sc_int$scDblFinder.class <- sce$scDblFinder.class

  # Plot UMAP with doublet scores
  p1 <- DimPlot2(sc_int, reduction = "umap_integrated", group.by = "scDblFinder.class") +
    ggtitle("Doublets Classification") +
    theme(plot.background = element_rect(fill = "white"))

  p2 <- DimPlot2(sc_int, features = "scDblFinder.score", reduction = "umap_integrated") +
    ggtitle("Doublet Scores") +
    theme(plot.background = element_rect(fill = "white"))

  p <- p1 + p2
  ggsave("results/101.load_data/doublets.png", p, width = 14, height = 7)

  sc_int <- sc_int %>%
    filter(scDblFinder.class == "singlet")
  sc_int
}

process_sc <- function(sc_raw, sc_doublet) {
  sc <- sc_raw %>% ligerToSeurat()
  cells_to_keep <- sc_doublet %>%
    filter(scDblFinder.class == "singlet") %>%
    pull(cell)
  sc <- sc %>%
    filter(cell %in% cells_to_keep)
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
