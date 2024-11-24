library(pagoda2)
library(Seurat)
library(velocyto.R)

ReadveloMatrix <- function(mex_dir = NULL,
                           barcode.path = NULL,
                           feature.path = NULL,
                           matrix.path = NULL,
                           use_10X = FALSE) {
    if (is.null(mex_dir) && is.null(barcode.path) && is.null(feature.path) &&
        is.null(matrix.path)) {
        stop("No matrix set.")
    }
    if (!is.null(mex_dir) && !file.exists(mex_dir)) {
        stop(paste0(mex_dir, " does not exist."))
    }
    if (is.null(barcode.path) && is.null(feature.path) && is.null(matrix.path)) {
        barcode.path <- paste0(mex_dir, "/barcodes.tsv.gz")
        feature.path <- paste0(mex_dir, "/features.tsv.gz")
        matrix.path <- paste0(mex_dir, "/matrix.mtx.gz")
    }
    spliced.path <- paste0(mex_dir, "/spliced.mtx.gz")
    unspliced.path <- paste0(mex_dir, "/unspliced.mtx.gz")
    spanning.path <- paste0(mex_dir, "/spanning.mtx.gz")

    if (!file.exists(barcode.path) || !file.exists(feature.path)) {
        stop(paste0("No expression file found at ", mex_dir))
    }

    .ReadPISA0 <- function(barcode.path, feature.path, matrix.path, use_10X) {
        mat <- Matrix::readMM(file = matrix.path)
        feature.names <- read.delim(feature.path,
            header = FALSE,
            stringsAsFactors = FALSE
        )
        barcode.names <- read.delim(barcode.path,
            header = FALSE,
            stringsAsFactors = FALSE
        )
        colnames(mat) <- barcode.names$V1
        if (use_10X == TRUE) {
            rownames(mat) <- make.unique(feature.names$V2)
        } else {
            rownames(mat) <- make.unique(feature.names$V1)
        }
        mat
    }

    if (!file.exists(spliced.path) && file.exists(matrix.path)) {
        return(.ReadPISA0(barcode.path, feature.path, matrix.path, use_10X))
    }
    mat <- list()
    cat("Load spliced matrix ...\n")
    mat$spliced <- Matrix::readMM(file = spliced.path)
    mat$spliced <- as(mat$spliced, "dgCMatrix")
    cat("Load unspliced matrix ...\n")
    mat$unspliced <- Matrix::readMM(file = unspliced.path)
    mat$unspliced <- as(mat$unspliced, "dgCMatrix")
    cat("Load spanning matrix ...\n")
    mat$spanning <- Matrix::readMM(file = spanning.path)
    mat$spanning <- as(mat$spanning, "dgCMatrix")

    feature.names <- read.delim(feature.path,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    barcode.names <- read.delim(barcode.path,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    colnames(mat$spliced) <- barcode.names$V1
    rownames(mat$spliced) <- make.unique(feature.names$V1)
    colnames(mat$unspliced) <- barcode.names$V1
    rownames(mat$unspliced) <- make.unique(feature.names$V1)
    colnames(mat$spanning) <- barcode.names$V1
    rownames(mat$spanning) <- make.unique(feature.names$V1)
    mat
}

EC <- ReadveloMatrix("data/103.self_workflow/N1/output/attachment/RNAvelocity_matrix")

emat <- EC$spliced
nmat <- EC$unspliced
smat <- EC$spanning
hist(log10(colSums(emat)), col = "wheat", xlab = "cell size")

emat <- emat[, colSums(emat) >= 1e3]

r <- Pagoda2$new(emat, modelType = "plain", trim = 10, log.scale = T)
r$adjustVariance(plot = T, do.par = T, gam.k = 10)

r$calculatePcaReduction(nPcs = 100, n.odgenes = 3e3, maxit = 300)
r$makeKnnGraph(k = 30, type = "PCA", center = T, distance = "cosine")
r$getKnnClusters(method = multilevel.community, type = "PCA", name = "multilevel")
r$getEmbedding(type = "PCA", embeddingType = "tSNE", perplexity = 50, verbose = T)

par(mfrow = c(1, 2))
r$plotEmbedding(type = "PCA", embeddingType = "tSNE", show.legend = F, mark.clusters = T, min.group.size = 10, shuffle.colors = F, mark.cluster.cex = 1, alpha = 0.3, main = "cell clusters")
r$plotEmbedding(type = "PCA", embeddingType = "tSNE", colors = r$depth, main = "depth")

emat <- EC$spliced
nmat <- EC$unspliced
emat <- emat[, rownames(r$counts)]
nmat <- nmat[, rownames(r$counts)]
cluster.label <- r$clusters$PCA$multilevel
library(sccore)
cell.colors <- fac2col(cluster.label)
emb <- r$embeddings$PCA$tSNE

cell.dist <- as.dist(1 - armaCor(t(r$reductions$PCA)))
emat <- filter.genes.by.cluster.expression(emat, cluster.label, min.max.cluster.average = 0.2)
nmat <- filter.genes.by.cluster.expression(nmat, cluster.label, min.max.cluster.average = 0.05)
length(intersect(rownames(emat), rownames(nmat)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 25, cell.dist = cell.dist, fit.quantile = fit.quantile)


library(ggplot2)
pdf("cell_velocity.pdf", height = 6, width = 8)
show.velocity.on.embedding.cor(emb, rvel.cd, n = 200, scale = "sqrt", cell.colors = ac(cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = F, cell.border.alpha = 0.1)
dev.off()
