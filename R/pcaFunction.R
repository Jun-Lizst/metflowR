#' @title PCAanalysis
#' @description PCA analysis for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param QC Use qc data for PCA analyis or not? Default is FALSE.
#' @param log.scale Which transformation methd you want to use? 2, 10  or "e",
#' defaulit is FALSE mean don't transformation.
#' @param scale.method Which scale methd you want to use? "auto" or "pareto",
#' defaulit is "auto".
#' @param path Work Directory
#' @param color Color for different class.
#' @param shape Shape for different class. please see the help of par(pch).
#' @param cex.lab cex of lab.
#' @param cex.axis cex of axis.
#' @param cex cex.
#' @param cex.text cex of text.
#' @param width width of plot.
#' @param height height of plot.
#' @param text Add text in PCA score plot or not? Deefault is FALSE.
#' @param ellipse Add ellipse in PCA score plot or not? Deefault is TRUE.
#' @param xlim1 x axis limitation. Default is NULL.
#' @param ylim1 y axis limitation. Default is NULL.
#' @return PCA score plot.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "metflowR")
#' data(sample.information, package = "metflowR")
#'
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'
#' # metflowR process
#' metflowR(#ImportData para
#' data = "data.csv",
#' sample.information = "sample.information.csv",
#' polarity = "positive",
#' #DataNormalization
#' method = "svr",
#' threads = 2)
#'## run
#'PCAanalysis(met.data.after.pre, QC = TRUE, path = "Demo for PCAanalysis")
#' }

PCAanalysis <- function(MetFlowData,
                        QC = FALSE,
                        log.scale = FALSE,
                        scale.method = "auto",
                        path = ".",
                        color = c(
                          "palegreen",
                          "royalblue",
                          "firebrick1",
                          "tan1",
                          "deepskyblue",
                          "cyan",
                          "gray48",
                          "chocolate4",
                          "darkmagenta",
                          "indianred1"
                        ),
                        shape = c(17, 19, 15, 18, 2, 8, 11, 13, 12, 14),
                        cex.lab = 1.3,
                        cex.axis = 1.3,
                        cex = 1.3,
                        cex.text = 1,
                        width = 7,
                        height = 7,
                        text = FALSE,
                        ellipse = TRUE,
                        xlim1 = NULL,
                        ylim1 = NULL) {

  if (path != ".") {
    dir.create(path)
  }
  subject <- MetFlowData@subject
  subject.info <- MetFlowData@subject.info
  group <- subject.info[, "group"]
  group.unique <- sort(unique(group))
  subject.name <- subject.info[, 1]
  qc <- MetFlowData@qc
  if (is.null(qc)) {
    QC <- FALSE
  }


  info <- list()
  for (i in 1:length(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }

  names(info) <- group.unique

  #select the subject in info and need QC or not
  index <- NULL
  for (i in seq_along(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index, index1)
  }

  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }

  index <- index[!is.na(index)]
  index <- match(index, colnames(subject))
  index <- index[!is.na(index)]
  subject <- subject[, index]


  ##discard the subject's name who is not in the subject data
  for (i in seq_along(info)) {
    idx <- as.character(info[[i]])
    idx <- match(idx, colnames(subject))
    idx <- idx[!is.na(idx)]
    info[[i]] <- colnames(subject)[idx]
  }

  if (QC) {
    int <- cbind(subject, qc)
  }
  else {
    int <- subject
  }

  ifelse(QC, int <- cbind(subject, qc) , int <- subject)
  name <- colnames(int)

  q <- grep("QC", name)

  ##log transformation
  if (log.scale == FALSE) {
    int <- int
  }

  if (log.scale == "e") {
    int <- log(int + 1)
  }

  if (log.scale != FALSE & log.scale != "e") {
    int <- log(int + 1, as.numeric(log.scale))
  }


  if (scale.method == "auto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x)) / sd(x)
    })
    int <- t(int)
  }
  if (scale.method == "pareto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    })
    int <- t(int)
  }

  if (scale.method == "no") {
    int <- int
  }
  if (scale.method == "center") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE))
      int <- t(int)
    })
  }

  ## PCA analysis
  int.pca <-
    prcomp(t(data.frame(int)),
           retx = TRUE,
           center = FALSE,
           scale = FALSE)
  sample.pca <- int.pca

  loading <- summary(int.pca)$rotation
  pov <- summary(int.pca)$importance[2,]
  sd <- summary(int.pca)$importance[1,]
  cp <- summary(int.pca)$importance[3,]
  pc <- int.pca$x

  pc1 <- pov[1]
  pc2 <- pov[2]
  pc3 <- pov[3]

  x <- pc[, 1]
  y <- pc[, 2]
  z <- pc[, 3]

  xmin <- 1.2 * min(x)
  xmax <- 1.2 * max(x)
  ymin <- 1.2 * min(y)
  ymax <- 1.2 * max(y)
  zmin <- 1.2 * min(z)
  zmax <- 1.2 * max(z)

  name <- colnames(subject)
  label <- list()
  for (i in seq_along(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
  }

  legend <- NULL
  for (i in seq_along(label)) {
    legend[label[[i]]] <- names(info)[i]
  }

  if (QC) {
    legend[(ncol(subject) + 1):(ncol(subject) + ncol(qc))] <- "QC"
  }

  colour <- NULL
  if (length(color) < length(info))
    stop("Color list is not enough")

  colourlist <- color
  for (i in seq_along(label)) {
    colour[label[[i]]] <- colourlist[i]
  }
  if (QC) {
    colour[(ncol(subject) + 1):(ncol(subject) + ncol(qc))] <-
      colourlist[length(info) + 1]
  }

  pcha <- NULL
  if (length(shape) < length(info))
    stop("Shape list is not enough")
  pchalist <- shape
  for (i in seq_along(label)) {
    pcha[label[[i]]] <- pchalist[i]
  }
  if (QC) {
    pcha[(ncol(subject) + 1):(ncol(subject) + ncol(qc))] <-
      pchalist[length(info) + 1]
  }

  #PCA 2D
  pdf(file.path(path, "PCA score plot.pdf"),
      width = width,
      height = height)
  #t1 vs t2 plot
  if (is.null(xlim1)) {
    xlim = c(xmin, xmax)
  } else {
      xlim = xlim1
  }
  if (is.null(ylim1)) {
    ylim = c(ymin, ymax)
  } else {
    ylim = ylim1
  }

  par(mar = c(5, 5, 4, 2))
  plot(
    x,
    y,
    xlim = xlim,
    ylim = ylim,
    col = colour,
    pch = pcha,
    xlab = paste("PC1:", pc1, sep = ""),
    ylab = paste("PC2:", pc2, sep = ""),
    cex = cex,
    cex.axis = cex.axis,
    cex.lab = cex.lab
  )
  if (text)
  {
    text(x, y, name, pos = 4, cex = cex.text)
  }

  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)

  if (ellipse)
  {
    # require(ellipse)
    lines(ellipse::ellipse(
      0,
      scale = c(sd(x), sd(y)),
      centre = c(mean(x), mean(y))
    ), lty = 2)
  }

  if (QC) {
    legend(
      "topleft",
      c(names(info), "QC"),
      pch = pchalist[1:(length(info) + 1)],
      col = colourlist[1:(length(info) + 1)],
      bty = "n",
      cex = 1.1
    )
  }
  else{
    legend(
      "topleft",
      names(info),
      pch = pchalist[1:length(info)],
      col = colourlist[1:length(info)],
      bty = "n",
      cex = 1.1
    )
  }

  #t2 vs t3 plot

  par(mar = c(5, 5, 4, 2))
  plot(
    y,
    z,
    xlim = xlim,
    ylim = ylim,
    col = colour,
    pch = pcha,
    xlab = paste("PC2:", pc2, sep = ""),
    ylab = paste("PC3:", pc3, sep = ""),
    cex = cex,
    cex.axis = cex.axis,
    cex.lab = cex.lab
  )
  if (text)
  {
    text(y, z, name, pos = 4, cex = cex.text)
  }
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  if (ellipse)
  {
    lines(ellipse::ellipse(
      0,
      scale = c(sd(y), sd(z)),
      centre = c(mean(y), mean(z))
    ), lty = 2)
  }


  if (QC) {
    legend(
      "topleft",
      c(names(info), "QC"),
      pch = pchalist[1:(length(info) + 1)],
      col = colourlist[1:(length(info) + 1)],
      bty = "n",
      cex = 1.1
    )
  }
  else{
    legend(
      "topleft",
      names(info),
      pch = pchalist[1:length(info)],
      col = colourlist[1:length(info)],
      bty = "n",
      cex = 1.1
    )
  }

  #t1 vs t3 plot
  plot(
    x,
    z,
    xlim = xlim,
    ylim = ylim,
    col = colour,
    pch = pcha,
    xlab = paste("PC1:", pc1, sep = ""),
    ylab = paste("PC3:", pc3, sep = ""),
    cex = 1.3,
    cex.axis = 1.3,
    cex.lab = 1.3
  )
  if (text)
  {
    text(x, z, name, pos = 4, cex = cex.text)
  }
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  if (ellipse)
  {
    lines(ellipse::ellipse(
      0,
      scale = c(sd(x), sd(z)),
      centre = c(mean(x), mean(z))
    ), lty = 2)
  }


  if (QC) {
    legend(
      "topleft",
      c(names(info), "QC"),
      pch = pchalist[1:(length(info) + 1)],
      col = colourlist[1:(length(info) + 1)],
      bty = "n",
      cex = 1.1
    )
  }
  else{
    legend(
      "topleft",
      names(info),
      pch = pchalist[1:length(info)],
      col = colourlist[1:length(info)],
      bty = "n",
      cex = 1.1
    )
  }

  #PCA 3D
  scatterplot3d::scatterplot3d(
    x,
    y,
    z,
    color = colour,
    xlab = paste("PC1:", pc1, sep = ""),
    ylab = paste("PC2:", pc2, sep = ""),
    zlab = paste("PC3:", pc3, sep = ""),
    angle = 40,
    pch = pcha,
    box = FALSE,
    cex.symbol = 1,
    cex.lab = 1.3,
    cex.axis = 1.3,
    xlim = xlim,
    ylim = ylim,
    zlim = c(zmin, zmax)
  )

  if (QC) {
    legend(
      "topleft",
      c(names(info), "QC"),
      pch = pchalist[1:(length(info) + 1)],
      col = colourlist[1:(length(info) + 1)],
      bty = "n",
      cex = 1.5
    )
  }
  else{
    legend(
      "topleft",
      names(info),
      pch = pchalist[1:length(info)],
      col = colourlist[1:length(info)],
      bty = "n",
      cex = 1.5
    )
  }
  dev.off()

}






SXTpca <- function(subject = NULL,
                   qc = NULL,
                   info = NULL,
                   #used data
                   QC = FALSE,
                   scale.method = "auto",
                   path = ".") {
  ##remove the peaks with zero ratio > 50
  temp.idx1 <- apply(subject, 1, function(x) sum(x == 0)/ncol(subject)) > 0.5
  temp.idx1 <- which(temp.idx1)
  if(length(temp.idx1) > 0){
    subject <- subject[-temp.idx1, ]
    if(!is.null(qc)) {
      qc <- qc[-temp.idx1,]
    }
  }

  if (path != ".") {
    dir.create(path)
  }
  if (any(is.na(subject)) |
      any(is.na(qc)))
    stop("Please impute MV in subject or QC samples.")
  if (is.null(subject))
    stop("Subject sample is NULL")
  if (!is.null(qc)) {
    if (nrow(subject) != nrow(qc))
      stop("ThSe row number of Subject and QC must same")
  }
  if (is.null(qc) & QC)
    stop("QC shoud be FALSE because qc is NULL")
  if (is.null(info))
    stop("Info must not be NULL")

  #select the subject in info and need QC or not
  index <- NULL
  for (i in seq_along(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index, index1)
  }

  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }

  index <- index[!is.na(index)]
  index <- match(index, colnames(subject))
  index <- index[!is.na(index)]
  subject <- subject[, index]


  ##discard the subject's name who is not in the subject data
  for (i in seq_along(info)) {
    idx <- as.character(info[[i]])
    idx <- match(idx, colnames(subject))
    idx <- idx[!is.na(idx)]
    info[[i]] <- colnames(subject)[idx]
  }

  if (QC) {
    int <- cbind(subject, qc)
  } else {
    int <- subject
  }

  ifelse(QC, int <- cbind(subject, qc) , int <- subject)
  name <- colnames(int)
  #
  q <- grep("QC", name)

  if (scale.method == "auto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x)) / sd(x)
    })
    int <- t(int)
  }
  if (scale.method == "pareto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    })
    int <- t(int)
  }

  if (scale.method == "no") {
    int <- int
  }
  if (scale.method == "center") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE))
      int <- t(int)
    })
  }



  ## PCA analysis
  int.pca <-
    prcomp(t(data.frame(int)),
           retx = TRUE,
           center = FALSE,
           scale = FALSE)
  sample.pca <- int.pca

  SXTpcaData <- list(
    sample.pca = sample.pca,
    subject = subject,
    qc = qc,
    info = info,
    QC = QC,
    scale.method = scale.method
  )
  class(SXTpcaData) <- "SXTpcaData"
  return(SXTpcaData)
}






##plot functioin for SXTpcaData
SXTpcaPlot <- function(SXTpcaData,
                       plot = "score",
                       loading.plot.name = "PCA loading",
                       score.plot.name = "PCA score",
                       color = c("palegreen","royalblue","firebrick1",
                                 "tan1","deepskyblue",
                                 "cyan","gray48", "chocolate4",
                                 "darkmagenta","indianred1"),
                       shape = c(17,19,15,18,2,8,11,13,12,14),
                       cexlab = 1.3,
                       cexaxis = 1.3,
                       cexa = 1.3,
                       cextext = 1,
                       width = 7,
                       height = 7,
                       text = FALSE,
                       ellipse = FALSE,
                       path = NULL) {

  #
  if (is.null(path)) path <- getwd()
  ## get data
  subject <- SXTpcaData[["subject"]]
  qc <- SXTpcaData[["qc"]]
  sample.pca = SXTpcaData[["sample.pca"]]
  info = SXTpcaData[["info"]]
  QC = SXTpcaData[["QC"]]
  int.pca <- sample.pca


  loading <- summary(int.pca)$rotation
  pov <- summary(int.pca)$importance[2, ]
  sd <- summary(int.pca)$importance[1, ]
  cp <- summary(int.pca)$importance[3, ]
  pc <- int.pca$x

  pc1 <- pov[1]
  pc2 <- pov[2]
  pc3 <- pov[3]

  x <- pc[, 1]
  y <- pc[, 2]
  z <- pc[, 3]

  xmin <- 1.2 * min(x)
  xmax <- 1.2 * max(x)
  ymin <- 1.2 * min(y)
  ymax <- 1.2 * max(y)
  zmin <- 1.2 * min(z)
  zmax <- 1.2 * max(z)

  name <- colnames(subject)
  label <- list()
  for (i in seq_along(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
  }

  legend <- NULL
  for (i in seq_along(label)) {
    legend[label[[i]]] <- names(info)[i]
  }

  if (QC) {
    legend[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <- "QC"
  }

  colour <- NULL
  if (length(color) < length(info))
    stop("Color list is not enough")

  colourlist <- color
  for (i in seq_along(label)) {
    colour[label[[i]]] <- colourlist[i]
  }
  if (QC) {
    colour[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <-
      colourlist[length(info) + 1]
  }

  pcha <- NULL
  if (length(shape) < length(info))
    stop("shape list is not enough")
  pchalist <- shape
  for (i in seq_along(label)) {
    pcha[label[[i]]] <- pchalist[i]
  }
  if (QC) {
    pcha[(ncol(subject)+1):(ncol(subject)+ncol(qc))] <-
      pchalist[length(info) + 1]
  }

  if (plot == "loading") {
    #laoding plot
    pdf(file.path(path, paste(loading.plot.name,".pdf", sep = "")),
        width = width,
        height = height)
    plot(
      loading[, 1],
      loading[, 2],
      pch = 20,
      xlab = "Component 1",
      ylab = "Component 2",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    plot(
      loading[, 2],
      loading[, 3],
      pch = 20,
      xlab = "Component 2",
      ylab = "Component 3",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    plot(
      loading[, 1],
      loading[, 3],
      pch = 20,
      xlab = "Component 1",
      ylab = "Component 3",
      cex.lab = cexlab,
      cex.axis = cexaxis
    )
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    #loading plot 3d
    scatterplot3d::scatterplot3d(
      loading[, 1],
      loading[, 2],
      loading[, 3],
      xlab = "Component 1",
      ylab = "Component 2",
      zlab = "Component 3",
      angle = 40,
      pch = 20,
      box = FALSE,
      cex.symbol = 1,
      cex.lab = 1.3,
      cex.axis = 0.8
    )
    dev.off()
  }

  if (plot == "score") {
    #PCA 2D
    pdf(file.path(path, paste(score.plot.name,".pdf", sep = "")),
        width = width,
        height = height)
    #t1plot
    par(mar = c(5, 5, 4, 2))
    plot(
      y,
      ylim = c(ymin, ymax),
      col = colour,
      pch = pcha,
      xlab = "QC injection order",
      ylab = paste("PC1:", pc1, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )

    abline(h = 0, lty = 2)


    #t1 vs t2 plot
    par(mar = c(5, 5, 4, 2))
    plot(
      x,
      y,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      col = colour,
      pch = pcha,
      xlab = paste("PC1:", pc1, sep = ""),
      ylab = paste("PC2:", pc2, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )
    if (text)
    {
      text(x, y, name, pos = 4, cex = cextext)
    }

    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)

    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(x), sd(y)),
        centre = c(mean(x), mean(y))
      ), lty = 2)
    }

    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #t2 vs t3 plot

    par(mar = c(5, 5, 4, 2))
    plot(
      y,
      z,
      xlim = c(ymin, ymax),
      ylim = c(zmin, zmax),
      col = colour,
      pch = pcha,
      xlab = paste("PC2:", pc2, sep = ""),
      ylab = paste("PC3:", pc3, sep = ""),
      cex = cexa,
      cex.axis = cexaxis,
      cex.lab = cexlab
    )
    if (text)
    {
      text(y, z, name, pos = 4, cex = cextext)
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(y), sd(z)),
        centre = c(mean(y), mean(z))
      ), lty = 2)
    }


    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #t1 vs t3 plot
    plot(x, z,
         xlim = c(xmin, xmax),
         ylim = c(zmin, zmax),
         col = colour,
         pch = pcha,
         xlab = paste("PC1:", pc1, sep = ""),
         ylab = paste("PC3:", pc3, sep = ""),
         cex = 1.3,
         cex.axis = 1.3,
         cex.lab = 1.3
    )
    if (text)
    {
      text(x, z, name, pos = 4, cex = cextext)
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    if (ellipse)
    {
      lines(ellipse::ellipse(
        0,
        scale = c(sd(x), sd(z)),
        centre = c(mean(x), mean(z))
      ), lty = 2)
    }


    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.1
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.1
      )
    }

    #PCA 3D
    scatterplot3d::scatterplot3d(x, y, z,
                                 color = colour,
                                 xlab = paste("PC1:", pc1, sep = ""),
                                 ylab = paste("PC2:", pc2, sep = ""),
                                 zlab = paste("PC3:", pc3, sep = ""),
                                 angle = 40,
                                 pch = pcha,
                                 box = FALSE,
                                 cex.symbol = 1,
                                 cex.lab = 1.3,
                                 cex.axis = 1.3,
                                 xlim = c(xmin, xmax),
                                 ylim = c(ymin, ymax),
                                 zlim = c(zmin, zmax)
    )

    if (QC) {
      legend(
        "topleft",
        c(names(info), "QC"),
        pch = pchalist[1:(length(info) + 1)],
        col = colourlist[1:(length(info) + 1)],
        bty = "n",
        cex = 1.5
      )
    }
    else{
      legend(
        "topleft",
        names(info),
        pch = pchalist[1:length(info)],
        col = colourlist[1:length(info)],
        bty = "n",
        cex = 1.5
      )
    }
    dev.off()
  }
}