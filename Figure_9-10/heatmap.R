library(tgcc)
library(pheatmap)

# generate FS model ground truth
set.seed(2024)
FSmodel <- tgcc:::make.fourspherical(n=400)
order <- order(FSmodel$label)
data <- FSmodel$data[order, ]
label <- FSmodel$label[order]
data0 <- FSmodel$groundtruth[order, ]

# tuning parameter of lambda and gamma
threshold <- 1e-05
maxIter <- 100
gammaSeq <- c(1, 1.5, 1.9)
lambdaSeq <- c(30, 60, 120)

# fit the spTGCC model
tgccFit <- tgcc:::spTGCC(
  data = data,
  lambdaSeq = lambdaSeq,
  gammaSeq = gammaSeq,
  threshold = threshold,
  maxIter = maxIter)

predLabel <- tgcc:::clusterLabel(tgccFit, numClusters=4)
1-mclust::classError(label, predLabel)$errorRate

# set color
range <- quantile(data[!is.na(data)], c(0.0, 0.5, 1))
colfun <- circlize::colorRamp2(range, c("blue", "white", "red"))
breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 256)
colvec <- colfun(breaks)

# show results
showHeatmap <- function(data, colvec) {
  pheatmap::pheatmap(
    mat = data,
    color = colvec,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    na_col = "grey65"
  )
}

# 500 * 500 sized
hplot0 <- showHeatmap(data0, colvec)
hplot1 <- showHeatmap(data, colvec)
hplot2 <- showHeatmap(tgccFit$theta[[1]], colvec)
hplot3 <- showHeatmap(tgccFit$theta[[2]], colvec)
hplot3 <- showHeatmap(tgccFit$theta[[3]], colvec)

# show colorbar
# 150 * 350 sized
nlevel <- 256
col_level <- seq(from = range[1], to = range[3], length.out = nlevel)
image(1, col_level, matrix(col_level, nrow = 1),
  col = gplots::bluered(nlevel),
  xlab = NA, ylab = NA, xaxt = "n")

################################################################################
# generate check board model ground truth
set.seed(2024)
CBmodel <- tgcc:::make.checkerboard(n = 400)
data <- CBmodel$data
label <- CBmodel$label
data0 <- CBmodel$groundTruth

# tuning parameter of lambda and gamma
lambdaSeq <- gammaSeq <- c(20, 50, 150)

# fit the biTGCC model
tgccFit <-
  tgcc:::biTGCC(
    data,
    lambdaSeq,
    gammaSeq,
    threshold = 1e-05 * nrow(data) * ncol(data),
    maxIter = 100)

# calculate the error rate
predLabel <- tgcc:::clusterLabel(tgccFit, numClusters=4)
1-mclust::classError(label, predLabel)$errorRate

# set color
range <- quantile(data[!is.na(data)], c(0, 0.5, 1))
colfun <- circlize::colorRamp2(range, c("blue", "white", "red"))
breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 256)
colvec <- colfun(breaks)

# show results
# 500 * 500 sized
hplot0 <- showHeatmap(data0, colvec)
hplot1 <- showHeatmap(data, colvec)
hplot2 <- showHeatmap(tgccFit$theta[[1]], colvec)
hplot3 <- showHeatmap(tgccFit$theta[[2]], colvec)
hplot4 <- showHeatmap(tgccFit$theta[[3]], colvec)

# show colorbar (150 * 450 sized)
nlevel <- 256
col_level <- seq(from = range[1], to = range[3], length.out = nlevel)
image(1, col_level, matrix(col_level, nrow = 1),
  col = gplots::bluered(nlevel),
  xlab = NA, ylab = NA, xaxt = "n")

