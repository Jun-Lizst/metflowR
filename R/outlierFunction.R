#' @title QCOutlierFilter
#' @description Using PCA to find QC outliers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param CI confidence interval.
#' @param qc.outlier.filter Remove outlier or not.
#' @param path Work directory.
#' @return MetFlowData whose qc outliers have been removed.
#' @seealso \code{\link{SubjectOutlierFilter}}
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
#'
#'
#'## run
#'new.met.data <- QCOutlierFilter(met.data.after.pre)
#' }

QCOutlierFilter <- function(MetFlowData,
                            CI = 0.95,
                            qc.outlier.filter = TRUE,
                            path = ".") {

  QCOutlierFiderData <- QCOutlierFinder(MetFlowData = MetFlowData,
                                        CI = CI,
                                        path = path)

  metData <- QCOutlierFiderData[[1]]
  obs.remove <- QCOutlierFiderData[[2]]

  data <- SplitBatch(MetFlowData = metData)
  qc1 <- data[[2]]

  if(qc.outlier.filter){
    for (i in seq_along(qc1)) {
      cat(paste("Batch", i))
      cat("\n")
      cat("-----------------------\n")
      temp.qc <- qc1[[i]]
      temp.idx <- obs.remove[[i]]
      if (length(temp.idx) != 0) {
        cat("QC shoulde be removed are:\n")
        cat(temp.idx)
        cat("\n")
        temp.idx <- temp.idx
          # readline(
          #   "Which QC you want to remove(please type the index of QC sample,
          #   and separate them using comma,
          #   if you don't want to remove any QC, please type n):"
          # )
# if (temp.idx == "n") {
#   temp.qc <- temp.qc
# } else {
  # temp.idx <- strsplit(temp.idx, split = ",")
  # temp.idx <- as.numeric(temp.idx[[1]])
  # temp.idx <- as.numeric(temp.idx)
  temp.qc <- temp.qc[, -temp.idx]
# }
      } else {
        temp.qc <- temp.qc
      }
      qc1[[i]] <- temp.qc
    }
  }


  qc2 <- qc1[[1]]
  if (length(qc1) > 1) {
    for (i in 2:length(qc1)) {
      qc2 <- cbind(qc2, qc1[[i]])
    }
  }

  ##remove QC information who have been removed from data
  qc.name <- colnames(qc2)
  qc.info <- metData@qc.info
  qc.order <- metData@qc.order
  qc.index <- which(is.na(match(qc.info[, 1], qc.name)))

  if (length(qc.index) != 0) {
    qc.info <- qc.info[-qc.index, ]
    qc.order <- qc.order[-qc.index]
  }

  metData@qc.info <- as.matrix(qc.info)
  metData@qc <- as.matrix(qc2)
  metData@qc.order <- as.numeric(qc.order)
  metData@qc.outlier.filter <- "yes"
  return(metData)

  }





# temp <- QCOutlierFinder(MetFlowData = met.data)
#QC outlier filtering according to zero ratio and PCA
QCOutlierFinder <- function(MetFlowData,
                            CI = 0.95,
                            path = ".") {
  options(warn = -1)
  if (path != ".") {
    dir.create(path)
  }
  qc <- MetFlowData@qc
  qc.info <- MetFlowData@qc.info
  tags <- MetFlowData@tags

  data <- SplitBatch(MetFlowData = MetFlowData)
  qc1 <- data[[2]]

  obs.remove <- list()
  ## PCA analysis
  for (i in seq_along(qc1)) {
    info <- list("QC" = colnames(qc1[[i]]))
    SXTpcaData <-
      SXTpca(
        subject = qc1[[i]],
        info = info,
        QC = FALSE,
        scale.method = "auto"
      )
    index2 <- SXTpcaFindOutlier(
      SXTpcaData = SXTpcaData,
      CI = CI,
      plot.name = paste("Batch", i, "outliers"),
      output.plot = TRUE,
      path = path
    )
    #
    ## The first and last QC shouldn't be removed
    if (any(c(1, ncol(qc1[[i]])) %in% index2)) {
      index2 <-
        index2[-match(c(1, ncol(qc1[[i]])), index2)[!is.na(match(c(1, ncol(qc1[[i]])), index2))]]
    }

    obs.remove[[i]] <- index2
  }
  QCOutlierFinderData <- list(MetFlowData = MetFlowData,
                              obs.remove = obs.remove)
  class(QCOutlierFinderData) <- "QCOutlierFinderData"
  options(warn = 1)
  return(QCOutlierFinderData)
}


#' @title SubjectOutlierFilter
#' @description Using PCA to filter subject outliers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param CI Confidence interva.
#' @param subject.outlier.filter Remove outliers or not.
#' @param path Work directory.
#' @return MetFlowData whose subject outliers have been removed.
#' @seealso \code{\link{QCOutlierFilter}}
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
#'new.met.data <- SubjectOutlierFilter(met.data.after.pre,
#'                                     path = "Demo for SubjectOutlierFilter")
#' }

SubjectOutlierFilter <- function(MetFlowData,
                                 CI = 0.95,
                                 subject.outlier.filter = TRUE,
                                 path = ".") {
  SubjectOutlierFiderData <- SubjectOutlierFinder(MetFlowData = MetFlowData,
                                                  CI = CI,
                                                  path = path)

  metData <- SubjectOutlierFiderData[[1]]
  obs.remove <- SubjectOutlierFiderData[[2]]

  data <- SplitBatch(MetFlowData = metData)
  subject1 <- data[[1]]


  if(subject.outlier.filter){
    for (i in seq_along(subject1)) {
      cat(paste("Batch",i))
      cat("\n")
      cat("-------------------------------------------\n")
      temp.subject <- subject1[[i]]
      temp.idx <- obs.remove[[i]]
      if (length(temp.idx) != 0) {
        #
        cat("Subject shoulde be removed are:")
        cat(temp.idx)
        cat("\n")
        temp.idx <- temp.idx
        # readline(
        #   "Which subject you want to remove(please type the index of subject sample,
        #   and separate them using comma,
        #   if you don't want to remove any subject, please type n):"
        # )
        # if (temp.idx == "n") {
        #   temp.subject <- temp.subject
        # } else {
          # temp.idx <- strsplit(temp.idx, split = ",")
          # temp.idx <- as.numeric(temp.idx[[1]])
          # temp.idx <- as.numeric(temp.idx)
          temp.subject <- temp.subject[, -temp.idx]
        # }
      } else {
        temp.subject <- temp.subject
      }
      subject1[[i]] <- temp.subject
    }
  }

  subject2 <- subject1[[1]]
  if (length(subject1) > 1) {
    for (i in 2:length(subject1)) {
      subject2 <- cbind(subject2, subject1[[i]])
    }
  }

  ##remove subject information who have been removed from data
  subject.name <- colnames(subject2)
  subject.info <- metData@subject.info
  subject.order <- metData@subject.order
  subject.index <- which(is.na(match(subject.info[, 1], subject.name)))

  if (length(subject.index) != 0) {
    subject.info <- subject.info[-subject.index, ]
    subject.order <- subject.order[-subject.index]
  }

  metData@subject.info <- as.matrix(subject.info)
  metData@subject <- as.matrix(subject2)
  metData@subject.order <- as.numeric(subject.order)
  metData@subject.outlier.filter <- "yes"
  return(metData)

  }











#subject outlier filtering according to zero ratio and PCA
SubjectOutlierFinder <- function(MetFlowData,
                                 CI = 0.95,
                                 path = "."){
  options(warn = -1)
  #
  if (path != ".") {
    dir.create(path)
  }
  subject <- MetFlowData@subject
  subject.info <- MetFlowData@subject.info
  tags <- MetFlowData@tags

  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]

  obs.remove <- list()
  ## PCA analysis
  for (i in seq_along(subject1)) {
    info <- list("Subject" = colnames(subject1[[i]]))
    SXTpcaData <- SXTpca(subject = subject1[[i]],
                         info = info, QC = FALSE, scale.method = "auto")
    index2 <- SXTpcaFindOutlier(SXTpcaData = SXTpcaData,
                                CI = CI,
                                plot.name = paste("Batch",i,"outliers"),
                                output.plot = TRUE,
                                path = path)
    obs.remove[[i]] <- index2
  }
  SubjectOutlierFinderData <- list(MetFlowData = MetFlowData,
                                   obs.remove = obs.remove)
  class(SubjectOutlierFinderData) <- "SubjectOutlierFinderData"
  options(warn = 0)
  return(SubjectOutlierFinderData)
}







##find outlier functions using SXTpcaData
SXTpcaFindOutlier <- function(SXTpcaData,
                              CI = 0.95,
                              output.plot = TRUE,
                              plot.name = "PCA score plot for outliers",
                              path = NULL) {


  if (is.null(path))path <- getwd()
  sample.pca <- SXTpcaData[["sample.pca"]]
  data <- SXTpcaData[["subject"]]

  loading <- summary(sample.pca)$rotation
  pov <- summary(sample.pca)$importance[2,]
  sd <- summary(sample.pca)$importance[1,]
  cp <- summary(sample.pca)$importance[3,]
  pc <- sample.pca$x

  pc1 <- round(pov[1], 2)
  pc2 <- round(pov[2], 2)
  pc3 <- round(pov[3], 2)

  x <- pc[, 1]
  y <- pc[, 2]
  z <- pc[, 3]

  xmin <- 1.2 * min(x)
  xmax <- 1.2 * max(x)
  ymin <- 1.2 * min(y)
  ymax <- 1.2 * max(y)
  zmin <- 1.2 * min(z)
  zmax <- 1.2 * max(z)

  ellipse.data <-
    SXTellipse(as.numeric(cor(x, y)),
               scale = c(sd(x), sd(y)),
               centre = c(mean(x), mean(y)))
  data.for.plot <- ellipse.data[["ellipse.data"]]
  t = ellipse.data[["t"]]
  a = ellipse.data[["a"]]
  d = ellipse.data[["d"]]
  scale = ellipse.data[["scale"]]
  centre = ellipse.data[["centre"]]

  ## right point
  point1 <- c(t * scale[1] * cos(0) + centre[1], t * scale[2] *
                cos(-d) + centre[2])
  ## top point
  point2 <- c(t * scale[1] * cos(pi / 2) + centre[1], t * scale[2] *
                cos(pi / 2 - d) + centre[2])
  ## left point
  point3 <- c(t * scale[1] * cos(pi) + centre[1], t * scale[2] *
                cos(pi - d) + centre[2])

  ellipse.a <- sqrt(sum((point1 - centre) ^ 2))
  ellipse.b <- sqrt(sum((point2 - centre) ^ 2))
  ellipse.c <- sqrt(ellipse.a ^ 2 - ellipse.b ^ 2)

  ## get the focus points

  lm.reg <- lm(c(point1[2], centre[2]) ~ c(point1[1], centre[1]))
  a <- lm.reg$coefficients[[2]]
  b <- lm.reg$coefficients[[1]]

  foo.f.a <- a ^ 2 + 1
  foo.f.b <-  2 * a * (b - centre[2]) - 2 * centre[1]
  foo.f.c <- centre[1] ^ 2 + (b - centre[2]) ^ 2 - ellipse.c ^ 2

  foo.f <- function(x,
                    a = foo.f.a,
                    b = foo.f.b,
                    c = foo.f.c)
  {
    a * x ^ 2 + b * x + c
  }
  result1 <-
    uniroot(
      foo.f,
      c(0, 10000),
      a = foo.f.a,
      b = foo.f.b,
      c = foo.f.c,
      tol = 0.0001
    )
  result2 <-
    uniroot(
      foo.f,
      c(-10000, 0),
      a = foo.f.a,
      b = foo.f.b,
      c = foo.f.c,
      tol = 0.0001
    )

  p1 <- c(result1$root, foo.f(x = result1$root))
  p2 <- c(result2$root, foo.f(x = result2$root))

  x1 <- data.for.plot[, 1]
  y1 <- data.for.plot[, 2]
  ellipse.standard <-
    mean(sqrt((x1 - p1[1]) ^ 2 + (y1 - p1[2]) ^ 2) +
           sqrt((x1 - p2[1]) ^ 2 + (y1 - p2[2]) ^ 2))

  distance <-
    sqrt((x - p1[1]) ^ 2 + (y - p1[2]) ^ 2) +
    sqrt((x - p2[1]) ^ 2 + (y - p2[2]) ^2)

  outlier.index <- which(distance > ellipse.standard)
  if (length(outlier.index) > 0) {
    cat(names(x)[outlier.index], " are outliers.\n")
  }

  if (output.plot) {
    colour <- rep(NA, length(x))
    colour[outlier.index] <- "firebrick1"
    colour[is.na(colour)] <- "black"
    pdf(file.path(path, paste(plot.name, ".pdf", sep = "")),
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    plot(
      x,
      y,
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      xlab = paste("PC1:", pc1),
      ylab = paste("PC2:", pc2),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colour,
      main = "PCA score plot for outliers"
    )
    if (length(outlier.index) > 0) {
      text(x = x[outlier.index],
           y = y[outlier.index],
           names(x)[outlier.index],
           pos = 4)
    } else
      (text(x = x,
            y = y,
            names(x),
            pos = 4))
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    lines(data.for.plot, lty = 2)
    dev.off()
  }

  return(outlier.index)

}
