load_sc <- function() {
  samples <- c("N1", "N2", "N4", "N7", "T1", "T2", "T4", "T7")
  sc_raw_list <- imap(samples, ~ Read10X(paste0("data/103.self_workflow/", .x, "/output/filter_matrix"), gene.column = 1))
  names(sc_raw_list) <- samples
  sc_raw <- createLiger(sc_raw_list)
  sc_raw
}

filter_low_quality_cells <- function(sc) {
  sc <- sc %>% ligerToSeurat()
  # 计算MAD并设置阈值
  calculate_mad_threshold <- function(x, nmads = 5) {
    median_x <- median(x)
    mad_x <- mad(x)
    lower_bound <- median_x - nmads * mad_x
    upper_bound <- median_x + nmads * mad_x
    return(c(lower_bound, upper_bound))
  }

  # 计算各指标的MAD阈值
  nUMI_bounds <- calculate_mad_threshold(log1p(sc$nUMI))
  nGene_bounds <- calculate_mad_threshold(log1p(sc$nGene))
  mito_bounds <- calculate_mad_threshold(sc$mito, nmads = 3)

  # 过滤细胞
  cells_to_keep <- sc@meta.data %>%
    mutate(
      pass_nUMI = log1p(nUMI) >= nUMI_bounds[1] & log1p(nUMI) <= nUMI_bounds[2],
      pass_nGene = log1p(nGene) >= nGene_bounds[1] & log1p(nGene) <= nGene_bounds[2],
      pass_mito = mito <= min(mito_bounds[2], 8)
    ) %>%
    mutate(pass_qc = pass_nUMI & pass_nGene & pass_mito)

  # 输出过滤统计
  message("原始细胞数: ", ncol(sc))
  message("过滤后细胞数: ", sum(cells_to_keep$pass_qc))

  # 绘制QC指标图
  p1 <- FeatureScatter(sc,
    feature1 = "nUMI",
    feature2 = "nGene",
    cols = ifelse(cells_to_keep$pass_qc, "gray", "red")
  ) +
    ggtitle(paste0("QC Metrics (", round(100*sum(cells_to_keep$pass_qc)/nrow(cells_to_keep), 1), "% cells passed)")) +
    theme(plot.background = element_rect(fill = "white"))

  ggsave("results/101.load_data/qc_metrics.png", p1, width = 10, height = 8)

  # 过滤数据
  sc <- subset(sc, cells = rownames(cells_to_keep)[cells_to_keep$pass_qc])

  return(sc)
}

process_sc_preview <- function(sc) {
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

process_sc <- function(sc, sc_doublet) {
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
