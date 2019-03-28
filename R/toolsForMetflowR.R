#' @title BatchEffectOverview
#' @description Using PCA score plot to view the batch effect.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData.before MetFlowData before normalization or integration.
#' @param MetFlowData.after MetFlowData after normalization or integration.
#' @param path work directory
#' @return Give PCA score plot for QC and subject.


### Batch effect for multiple batch datasets
BatchEffectOverview <- function(MetFlowData.before,
                                MetFlowData.after,
                                path = ".") {
  options(warn = -1)
  #
  if (path != ".") {
    dir.create(path)
  }

  hasQC <- MetFlowData.before@"hasQC"
  subject.info <- MetFlowData.before@"subject.info"
  qc.info <- MetFlowData.before@"qc.info"

  subject.bef <- MetFlowData.before@"subject"
  qc.bef <- MetFlowData.before@"qc"

  subject.aft <- MetFlowData.after@"subject"
  qc.aft <- MetFlowData.after@"qc"

  if ((sum(is.na(subject.bef)) + sum(is.na(qc.bef))) != 0) {
    stop("The data has MV, please do MV imputation first.")
  }

  data.bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject.bef1 <- data.bef[[1]]
  qc.bef1 <- data.bef[[2]]
  subject.info.bef1 <- data.bef[[3]]
  qc.info.bef1 <- data.bef[[4]]

  data.aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject.aft1 <- data.aft[[1]]
  qc.aft1 <- data.aft[[2]]
  subject.info.aft1 <- data.aft[[3]]
  qc.info.aft1 <- data.aft[[4]]

  if (hasQC != "no") {
    ##PCA analysis  for QC samples
    qc.name <- qc.info[, 1]
    qc.batch <- as.numeric(qc.info[, 4])

    qc.info <- list()
    for (i in seq_along(qc.bef1)) {
      qc.info[[i]] <- colnames(qc.bef1[[i]])
    }


    names(qc.info) <- paste("Batch", c(1:length(qc.info)), sep = "")
    ##before
    qc.SXTpcaData <- SXTpca(
      subject = qc.bef,
      info = qc.info,
      QC = FALSE,
      scale.method = "auto"
    )

    SXTpcaPlot(
      SXTpcaData = qc.SXTpcaData,
      score.plot.name = "Before Batch effect in QC PCA",
      ellipse = TRUE,
      path = path
    )
    ##after
    qc.SXTpcaData <- SXTpca(
      subject = qc.aft,
      info = qc.info,
      QC = FALSE,
      scale.method = "auto"
    )

    SXTpcaPlot(
      SXTpcaData = qc.SXTpcaData,
      score.plot.name = "After Batch effect in QC PCA",
      ellipse = TRUE,
      path = path
    )
  }
  ##PCA analysis for subject samples
  subject.info <- list()
  for (i in seq_along(subject.bef1)) {
    subject.info[[i]] <- colnames(subject.bef1[[i]])
  }
  names(subject.info) <-
    paste("Batch", c(1:length(subject.info)), sep = "")

  ##before
  subject.SXTpcaData <- SXTpca(
    subject = subject.bef,
    info = subject.info,
    QC = FALSE,
    scale.method = "auto"
  )

  SXTpcaPlot(
    SXTpcaData = subject.SXTpcaData,
    score.plot.name = "Before Batch effect in Subject PCA",
    ellipse = TRUE,
    path = path
  )

  ##after
  subject.SXTpcaData <- SXTpca(
    subject = subject.aft,
    info = subject.info,
    QC = FALSE,
    scale.method = "auto"
  )

  SXTpcaPlot(
    SXTpcaData = subject.SXTpcaData,
    score.plot.name = "After Batch effect in Subject PCA",
    ellipse = TRUE,
    path = path
  )

  if (hasQC != "no") {
    ## QC sum intensity distribution
    colours.list <-
      c(
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
      )
    colours <- NULL
    for (i in seq_along(qc.bef1)) {
      colours[match(colnames(qc.bef1[[i]]), colnames(qc.bef))] <-
        colours.list[i]
    }

    pdf(
      file.path(path, "QC total intensity distribution.pdf"),
      width = 14,
      height = 7
    )
    layout(matrix(c(1:2), ncol = 2))
    par(mar = c(5, 5, 4, 2))
    ##before
    plot(
      colSums(qc.bef),
      col = colours,
      pch = 19,
      xlab = "QC injection order",
      ylab = "Total intensity",
      cex.lab = 1.3,
      cex.axis = 1.3,
      cex = 1.3,
      ylim = c(0.5 * min(colSums(qc.bef)), 1.5 * max(colSums(qc.bef))),
      main = "Before"
    )

    qc.scale <- t(apply(qc.bef, 1, function(x) {
      (x - mean(x)) / sd(x)
    }))
    boxplot(
      qc.scale,
      col = colours,
      xlab = "QC index",
      ylab = "Intensity (auto scaled)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      notch = FALSE,
      outline = FALSE,
      main = "Before"
    )

    ##after
    plot(
      colSums(qc.aft),
      col = colours,
      pch = 19,
      xlab = "QC injection order",
      ylab = "Total intensity",
      cex.lab = 1.3,
      cex.axis = 1.3,
      cex = 1.3,
      ylim = c(0.5 * min(colSums(qc.bef)), 1.5 * max(colSums(qc.bef))),
      main = "After"
    )

    qc.scale <- t(apply(qc.aft, 1, function(x) {
      (x - mean(x)) / sd(x)
    }))
    boxplot(
      qc.scale,
      col = colours,
      xlab = "QC index",
      ylab = "Intensity (auto scaled)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      notch = FALSE,
      outline = FALSE,
      main = "After"
    )

    dev.off()
  }
  options(warn = -1)
}


#' @title DataOverview
#' @description Give a overview of MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData
#' @param feature.distribution Draw a rt vs mz vs intensity plot or not.
#' Default is TRUE.
#' @param path Work directory.
#' @return Data overview_RT vs mz vs intensity.pdf:
#' A RT vs mz vs intensity plot.
#' @return Data overview.txt: A overview information for MetFlowData.
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
#'## run DataOverview
#'DataOverview(MetFlowData = met.data.after.pre, path = "Demo for DataOverview")
#'}

### DataOverview for MeeFlowData
DataOverview <- function(MetFlowData,
                         feature.distribution = TRUE,
                         path = ".") {
  if (path != ".") {
    dir.create(path)
  }
  #
  hasQC <- MetFlowData@hasQC
  subject <- MetFlowData@subject
  qc <- MetFlowData@qc
  tags <- MetFlowData@tags
  subject.info <- MetFlowData@subject.info
  qc.info <- MetFlowData@qc.info
  subject.order <- as.numeric(MetFlowData@subject.order)
  qc.order <- as.numeric(MetFlowData@qc.order)

  #### a 3D figure: mz vs RT vs intensity
  if (feature.distribution) {
    mz <- as.numeric(tags[, "mz"])
    rt <- as.numeric(tags[, "rt"])
    if (hasQC != "no") {
      sample <- data.frame(qc, subject)
    }
    else {
      sample <- subject
    }

    int.log <-
      log(apply(sample, 1, function(x) {
        mean(x, na.rm = TRUE)
      }) + 10, 10)


    rt.mz.int <- data.frame(rt, mz, int.log)

    # library(ggplot2)
    par(mar = c(5, 5, 4, 2))
    rt.mz.int <-
      ggplot2::ggplot(data = rt.mz.int,
                      ggplot2::aes(x = rt, y = mz, colour = int.log)) +
      ggplot2::geom_point(alpha = 0.3) +
      ggplot2::scale_color_gradient(low = "green", high = "red") +
      ggplot2::labs(x = "Retention time (RT)",
           y = "Mass to charge ratio (m/z)",
           colour = "log10(intensity)") +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14),
            axis.title.y = ggplot2::element_text(size = 14)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_text(size = 10)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 14)) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 10)) +
      ggplot2::ggtitle("RT vs mz vs Intensity") +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      filename = file.path(path, "Data overview_RT vs mz vs intensity.pdf"),
      plot = rt.mz.int,
      width = 8,
      height = 6
    )

  }


  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  file.name.for.txt <- file.path(path, "Data overview.txt")

  ## begin output data overview
  cat("The overview of data\n", file = file.name.for.txt, append = FALSE)
  cat("---------------\n", file = file.name.for.txt, append = TRUE)
  cat('\n', file = file.name.for.txt, append = TRUE)
  cat('\n', file = file.name.for.txt, append = TRUE)

  ## batch information
  cat("There are",
      length(subject1),
      "batches",
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  ##
  subject.number <- unlist(lapply(subject1, length))
  names(subject.number) <- paste("Batch", seq_along(subject1))
  if (hasQC == "no") {
    qc.number <- NULL
  }
  else {
    qc.number <- unlist(lapply(qc1, length))
    names(qc.number) <- paste("Batch", seq_along(subject1))
  }

  cat("Subject number in each Batch:\n",
      file = file.name.for.txt,
      append = TRUE)
  options(warn = -1)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat(subject.number, file = file.name.for.txt, append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat(qc.number, file = file.name.for.txt, append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  options(warn = 0)

  ## peak information
  cat("---------------\n", file = file.name.for.txt, append = TRUE)
  cat(paste("Peak number:", nrow(subject)),
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)

  cat("The tags information contains:\n",
      file = file.name.for.txt,
      append = TRUE)
  cat(colnames(tags), file = file.name.for.txt, append = TRUE)
  cat('\n', file = file.name.for.txt, append = TRUE)
  cat('\n', file = file.name.for.txt, append = TRUE)

  ## subject sample information
  cat("---------------\n", file = file.name.for.txt, append = TRUE)
  cat("Subject sample info:\n", file = file.name.for.txt, append = TRUE)
  cat(paste("Subject sample number:", ncol(subject)),
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)

  cat("Subject sample name:\n", file = file.name.for.txt, append = TRUE)
  cat(colnames(subject), file = file.name.for.txt, append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)

  if (hasQC != "no") {
    ## QC sample information
    cat("\n", file = file.name.for.txt, append = TRUE)
    cat("---------------\n", file = file.name.for.txt, append = TRUE)
    cat("QC sample info:\n", file = file.name.for.txt, append = TRUE)
    cat(paste("QC sample number:", ncol(qc)),
        file = file.name.for.txt,
        append = TRUE)
    cat("\n", file = file.name.for.txt, append = TRUE)

    cat("QC sample name:\n", file = file.name.for.txt, append = TRUE)
    cat(colnames(qc), file = file.name.for.txt, append = TRUE)
    cat("\n", file = file.name.for.txt, append = TRUE)
  }
  else {
    cat("QC sample info:\n", file = file.name.for.txt, append = TRUE)
    cat("No QC in data", file = file.name.for.txt, append = TRUE)
    cat("\n", file = file.name.for.txt, append = TRUE)
  }

  ## other processing information
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("---------------\n", file = file.name.for.txt, append = TRUE)
  cat("Some processing information\n",
      file = file.name.for.txt,
      append = TRUE)
  cat("MV imputation:",
      MetFlowData@mv.imputation,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Imputation method:",
      MetFlowData@imputation.method,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Zero filter:",
      MetFlowData@zero.filter,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Zero filter criteria:",
      MetFlowData@zero.filter.criteria,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("QC outlier filter:",
      MetFlowData@qc.outlier.filter,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Normalization:",
      MetFlowData@normalization,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Normalization method:",
      MetFlowData@normalization.method,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Data integration:",
      MetFlowData@data.integration,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
  cat("Data integration method:",
      MetFlowData@data.integration.method,
      file = file.name.for.txt,
      append = TRUE)
  cat("\n", file = file.name.for.txt, append = TRUE)
}












#' @title ExportData
#' @description Export MetFlowData as csv.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param data.name The name of the data you want to output. Default is
#' "data_new".
#' @param subject.info.name The name for subject information you want to
#' output. Default is "subject.info".
#' @param qc.info.name The name for QC information you want to output.
#' Default is "qc.info".
#' @param path Work directory.
#' @return Write csv data.
#' @export
#' @seealso \code{\link{ImportData}}
#' @examples
#' #load the demo data
#' data(met.data, package = "metflowR")
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'## run
#'ExportData(MetFlowData = met.data, path = "Demo for ExportData")

ExportData <- function(MetFlowData,
                       data.name = "data_new",
                       subject.info.name = "subject.info",
                       qc.info.name = "qc.info",
                       path = ".") {

  if (path != ".") {
    dir.create(path)
  }

  subject <- MetFlowData@subject
  qc <- MetFlowData@qc
  tags <- MetFlowData@tags
  subject.info <- MetFlowData@subject.info
  qc.info <- MetFlowData@qc.info

  has.qc <- MetFlowData@hasQC

  if (has.qc == "yes") {
    write.csv(cbind(tags, subject, qc),
              file.path(path, paste(data.name, ".csv", sep = "")),
              row.names = FALSE)
    write.csv(subject.info,
              file.path(path, paste(subject.info.name, ".csv", sep = "")),
              row.names = FALSE)
    write.csv(qc.info,
              file.path(path, paste(qc.info.name, ".csv", sep = "")),
              row.names = FALSE)
  } else{
    write.csv(cbind(tags, subject),
              file.path(path, paste(data.name, ".csv", sep = "")),
              row.names = FALSE)
    write.csv(subject.info,
              file.path(path, paste(subject.info.name, ".csv", sep = "")),
              row.names = FALSE)
  }
}


#------------------------------------------------------------------------------
#' @title ImportData
#' @description Import data and retrun a standard MetProcessor data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data Data name for analysis. Default is "data.csv". Data are csv
#' format from XCMS, MZmine or other software.
#'  Please see the demo data in example.
#' @param sample.information Sample information name for analysis. Default is
#' "sample.information.csv". Column 1 is sample.name, column 2 is
#' injection.order, column 3 is class ("Subject" or "QC"), column 4 is
#' batch information, column 5 is group ("control" or "case"), other columns
#' are information for sample. Please see demo data in example.
#' @param polarity The polarity of data, "positive", "negative"" or "none"",
#' default is positive.
#' @param hasQC The data has QC samples or not? Default is "yes".
#' @param hasIS The data has IS or not? Default is "no".
#' @param path Work directory.
#' @param peak.identification The data has identification result or not?
#' Default is "no".
#' @return  Return a standard MetProcesser dataset.
#' @export
#' @examples
#' #load the demo data
#' data(data, package = "metflowR")
#' data(sample.information, package = "metflowR")
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'#Import data
#'met.data <- ImportData()

ImportData <- function(data = "data.csv",
                       sample.information = "sample.information.csv",
                       polarity = c("positive", "negative", "none"),
                       # posfix = NULL,
                       # ordered.qc = FALSE,
                       # worklist.from = "manual",
                       hasIS = c("yes", "no"),
                       hasQC = c("yes", "no"),
                       peak.identification = c("yes", "no"),
                       path = ".") {
  polarity <- match.arg(polarity)
  hasIS <- match.arg(hasIS)
  hasQC <- match.arg(hasQC)
  peak.identification <- match.arg(peak.identification)
  #
  if (path != ".") {
    dir.create(path)
  }

  # if (worklist.from != "GetWorklist")
  # {
  #   data <- ChangeSampleName(
  #     data = data,
  #     sample.information = sample.information,
  #     polarity = polarity,
  #     posfix = posfix,
  #     ordered.qc = ordered.qc,
  #     output = FALSE,
  #     path = path
  #   )
  #   sample.information <-
  #     read.csv(
  #       "sample.information1.csv",
  #       stringsAsFactors = FALSE,
  #       check.names = FALSE
  #     )
  # } else {
  data <- readr::read_csv(file.path(path, data),
                          col_types = readr::cols(),
                          progress = FALSE)
  data <- as.data.frame(data)

  # data <-
  #   read.csv(file.path(path, data),
  #            stringsAsFactors = FALSE,
  #            check.names = FALSE)
  if (sum(duplicated(colnames(data))) > 0) {
    stop("There are duplicated samples (names) in you data.")
  }
  sample.information <-
    read.csv(
      file.path(path, sample.information),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  # }

  ## read sample information
  sample.information <-
    sample.information[!is.na(sample.information[, 1]), ]

  ## sort sample information according to sample order
  sample.information <-
    sample.information[order(as.numeric(sample.information[, 2])), ]


  ## sort sample in data according to sample order
  ##
  sample <- data[,match(sample.information$sample.name, colnames(data))]

  # sample.index <- grep("Sample", colnames(data))
  # sample <- data[, sample.index]
  # tags <- data[, -sample.index]
  tags <- data[,-match(sample.information$sample.name, colnames(data))]
  sample.name <- colnames(sample)
  # sample.order <-
  #   as.numeric(substr(
  #     x = sample.name,
  #     start = 7,
  #     stop = unlist(lapply(gregexpr("_", sample.name), function(x) {
  #       x[1][1]
  #     })) - 1
  #   ))
  # sample <- sample[, order(sample.order)]
  data <- cbind(tags, sample)

  ## get subject, qc and tags
  # data.name <- colnames(data)
  # sample.index <- grep("Sample", data.name)
  # sample <- data[, sample.index]
  # tags <- data[, -sample.index]
  # sample.name <- colnames(sample)

  qc.index <- grep("QC", sample.information$class)
  ## has QC or not?
  if (length(qc.index) == 0) {
    hasQC = "no"
  }else{
    hasQC = "yes"
  }

  if (hasQC == "no") {
    qc <- NULL
    subject <- sample
    qc.name <- NULL
  }else {
    qc <- sample[, qc.index]
    subject <- sample[, -qc.index]
    qc.name <- colnames(qc)
  }

  subject.name <- colnames(subject)

  if (hasQC == "no") {
    qc.order <- NULL
  }else{
    qc.order <- sample.information$injection.order[grep("QC", sample.information$class)]
  }

  subject.order <- sample.information$injection.order[grep("Subject", sample.information$class)]

  subject.name <- colnames(subject)

  if (hasQC == "no") {
    qc.name <- NULL
  } else {
    qc.name <-
      qc.name <- colnames(qc)
  }
  colnames(subject) <- subject.name
  names(subject.order) <- subject.name
  if (hasQC != "no") {
    colnames(qc) <- qc.name
    names(qc.order) <- qc.name
  }

  subject.info <- sample.information[sample.information[, 3] == "Subject", ]
  if (hasQC == "no") {
    qc.info <- NULL
  } else {
    qc.info <- sample.information[sample.information[, 3] == "QC", ]
  }

  if ((sum(is.na(subject)) + sum(is.na(qc))) == 0) {
    mv.imputation <- "yes"
    imputation.method <- "default"
  } else {
    mv.imputation <- "no"
    imputation.method <- "no"
  }

  name <- tags[, "name"]
  if (polarity == "positive"){
    tags[, "name"] <- paste(name, "POS", sep = "_")
    pol <- rep("POS", nrow(tags))
    tags <- data.frame(pol, tags)
    colnames(tags)[1] <- "polarity"
  }
  if (polarity == "negative") {
    tags[, "name"] <- paste(name, "NEG", sep = "_")
    pol <- rep("NEG", nrow(tags))
    tags <- data.frame(pol, tags)
    colnames(tags)[1] <- "polarity"
  }
  #

  MetFlowData <- new("MetFlowData",
                     subject = as.matrix(subject),
                     qc = as.matrix(qc),
                     tags = as.matrix(tags),
                     tags.old = as.matrix(tags),
                     subject.info = as.matrix(subject.info),
                     qc.info = as.matrix(qc.info),
                     subject.order = as.numeric(subject.order),
                     qc.order = as.numeric(qc.order),
                     ## some preprocessing information
                     mv.imputation = mv.imputation,
                     imputation.method = imputation.method,
                     zero.filter = "no",
                     zero.filter.criteria = "no",
                     normalization = "no",
                     normalization.method = "no",
                     data.integration = "no",
                     data.integration.method = "no",
                     hasIS = hasIS,
                     hasQC = hasQC,
                     peak.identification = peak.identification,
                     foldchange = "temp",
                     marker.selection.condition = "no",
                     mv.filter = "no",
                     mv.filter.criteria = "no",
                     univariate.test = "no",
                     qc.outlier.filter = "no",
                     subject.outlier.filter = "no")

  MetFlowData <- MetFlowData
}





#' #------------------------------------------------------------------------------
#' #' @title ImportData
#' #' @description Import data and retrun a standard MetProcessor data.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@sioc.ac.cn}
#' #' @param data Data name for analysis. Default is "data.csv". Data are csv
#' #' format from XCMS, MZmine or other software.
#' #'  Please see the demo data in example.
#' #' @param sample.information Sample information name for analysis. Default is
#' #' "sample.information.csv". Column 1 is sample.name, column 2 is
#' #' injection.order, column 3 is class ("Subject" or "QC"), column 4 is
#' #' batch information, column 5 is group ("control" or "case"), other columns
#' #' are information for sample. Please see demo data in example.
#' #' @param polarity The polarity of data, "positive", "negative"" or "none"",
#' #' default is positive.
#' #' @param hasQC The data has QC samples or not? Default is "yes".
#' #' @param hasIS The data has IS or not? Default is "no".
#' #' @param posfix Default is NULL.
#' #' @param ordered.qc Default is FALSE.
#' #' @param worklist.from Default is "manual".
#' #' @param path Work directory.
#' #' @param peak.identification The data has identification result or not?
#' #' Default is "no".
#' #' @return  Return a standard MetProcesser dataset.
#' #' @export
#' #' @examples
#' #' #load the demo data
#' #' data(data, package = "metflowR")
#' #' data(sample.information, package = "metflowR")
#' #' # export the demo data as csv
#' #' write.csv(data, "data.csv", row.names = FALSE)
#' #' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#' #'#Import data
#' #'met.data <- ImportData()
#'
#' ImportData <- function(data = "data.csv",
#'                        sample.information = "sample.information.csv",
#'                        polarity = c("positive", "negative", "none"),
#'                        posfix = NULL,
#'                        ordered.qc = FALSE,
#'                        worklist.from = "manual",
#'                        hasIS = "no",
#'                        hasQC = "yes",
#'                        peak.identification = "no",
#'                        path = ".") {
#'   polarity <- match.arg(polarity)
#'   #
#'   if (path != ".") {
#'     dir.create(path)
#'   }
#'   if (worklist.from != "GetWorklist")
#'   {
#'     data <- ChangeSampleName(
#'       data = data,
#'       sample.information = sample.information,
#'       polarity = polarity,
#'       posfix = posfix,
#'       ordered.qc = ordered.qc,
#'       output = FALSE,
#'       path = path
#'     )
#'     sample.information <-
#'       read.csv(
#'         "sample.information1.csv",
#'         stringsAsFactors = FALSE,
#'         check.names = FALSE
#'       )
#'   }
#'
#'   else {
#'     data <-
#'       readr::read_csv(file.path(path, data),
#'                       col_types = readr::cols(),
#'                       progress = FALSE)
#'     data <- as.data.frame(data)
#'
#'     # data <-
#'     #   read.csv(file.path(path, data),
#'     #            stringsAsFactors = FALSE,
#'     #            check.names = FALSE)
#'     if (sum(duplicated(colnames(data))) > 0) {
#'       stop("There are duplicated samples (names) in you data.")
#'     }
#'     sample.information <-
#'       read.csv(
#'         file.path(path, sample.information),
#'         stringsAsFactors = FALSE,
#'         check.names = FALSE
#'       )
#'   }
#'
#'   ## read sample information
#'   sample.information <-
#'     sample.information[!is.na(sample.information[, 1]), ]
#'
#'   ## sort sample information according to sample order
#'   sample.information <-
#'     sample.information[order(as.numeric(sample.information[, 2])), ]
#'   write.csv(sample.information,
#'             file.path(path, "sample.information1.csv"),
#'             row.names = FALSE)
#'
#'   ## sort sample in data according to sample order
#'   sample.index <- grep("Sample", colnames(data))
#'   sample <- data[, sample.index]
#'   tags <- data[, -sample.index]
#'   sample.name <- colnames(sample)
#'   sample.order <-
#'     as.numeric(substr(
#'       x = sample.name,
#'       start = 7,
#'       stop = unlist(lapply(gregexpr("_", sample.name), function(x) {
#'         x[1][1]
#'       })) - 1
#'     ))
#'   sample <- sample[, order(sample.order)]
#'   data <- cbind(tags, sample)
#'
#'   ## get subject, qc and tags
#'   data.name <- colnames(data)
#'   sample.index <- grep("Sample", data.name)
#'   sample <- data[, sample.index]
#'   tags <- data[, -sample.index]
#'   sample.name <- colnames(sample)
#'
#'   qc.index <- grep("QC", sample.name)
#'   ## has QC or not?
#'   if (length(qc.index) == 0) {
#'     hasQC = "no"
#'   }
#'
#'   if (hasQC == "no") {
#'     qc <- NULL
#'     subject <- sample
#'     qc.name <- NULL
#'   }else {
#'     qc <-
#'       sample[, qc.index]
#'     subject <- sample[, -qc.index]
#'     qc.name <- colnames(qc)
#'   }
#'
#'   subject.name <- colnames(subject)
#'
#'   ## the sample name from sample information vs the sample name from data
#'   a <- setdiff(sample.information[, 1], sample.name)
#'   b <- setdiff(sample.name, sample.information[, 1])
#'   if (length(a) != 0)
#'   {
#'     stop(paste(
#'       paste(a, collapse = " "),
#'       "are in sample information but not in data."
#'     ))
#'   }
#'
#'   if (length(b) != 0)
#'   {
#'     stop(paste(
#'       paste(b, collapse = " "),
#'       "are in data but not in sample information."
#'     ))
#'   }
#'
#'   ## the sample order form sample information
#'   ## vs the sample order form sample name
#'   sample.name.from.sample.information <- sample.information[, 1]
#'   sample.order.from.sample.information <-
#'     substr(x = sample.name.from.sample.information,
#'            start = 7,
#'            stop = unlist(lapply(gregexpr("_", sample.name.from.sample.information),
#'                                 function(x) {
#'                                   x[1][1]
#'                                 })) - 1)
#'
#'   sample.different.index <-
#'     which(sample.order.from.sample.information != sample.information[, 2])
#'   if (length(sample.different.index) != 0)
#'   {
#'     print(sample.information[sample.different.index, c(1, 2)])
#'     stop("The sample order from data and from
#'          sample information are different.")
#'   }
#'
#'   if (hasQC == "no") {
#'     qc.order <- NULL
#'   }else {
#'     qc.order <-
#'       substr(x = qc.name,
#'              start = 7,
#'              stop = unlist(lapply(gregexpr("_", qc.name), function(x) {
#'                x[1][1]
#'              })) - 1)
#'   }
#'
#'   subject.order <-
#'     substr(x = subject.name,
#'            start = 7,
#'            stop = unlist(lapply(gregexpr("_", subject.name), function(x) {
#'              x[1][1]
#'            })) - 1)
#'
#'   subject.name <-
#'     substr(subject.name,
#'            start = unlist(lapply(gregexpr("_", subject.name), function(x) {
#'              x[1][1]
#'            })) + 1,
#'            stop = unlist(lapply(gregexpr("_", subject.name), function(x) {
#'              x[2][1]
#'            })) - 1)
#'
#'   if (hasQC == "no") {
#'     qc.name <- NULL
#'   } else {
#'     qc.name <-
#'       substr(qc.name,
#'              start = unlist(lapply(gregexpr("_", qc.name), function(x) {
#'                x[1][1]
#'              })) + 1,
#'              stop = unlist(lapply(gregexpr("_", qc.name), function(x) {
#'                x[2][1]
#'              })) - 1)
#'   }
#'   colnames(subject) <- subject.name
#'   names(subject.order) <- subject.name
#'   if (hasQC != "no") {
#'     colnames(qc) <- qc.name
#'     names(qc.order) <- qc.name
#'   }
#'   ## remove sample order form sample.name.from.data and
#'   ## sample.name.from.sample.information
#'   sample.name.from.data <- colnames(sample)
#'   sample.name.from.sample.information <- sample.information[, 1]
#'   sample.name.from.data <-
#'     substr(sample.name.from.data,
#'            start = unlist(lapply(gregexpr("_", sample.name.from.data),
#'                                  function(x) {
#'                                    x[1][1]
#'                                  })) + 1,
#'            stop = unlist(lapply(gregexpr("_", sample.name.from.data),
#'                                 function(x) {
#'                                   x[2][1]
#'                                 })) - 1)
#'
#'   sample.name.from.sample.information <-
#'     substr(
#'       sample.name.from.sample.information,
#'       start = unlist(lapply(gregexpr("_", sample.name.from.sample.information),
#'                             function(x) {
#'                               x[1][1]
#'                             })) + 1,
#'       stop = unlist(lapply(gregexpr("_", sample.name.from.sample.information),
#'                            function(x) {
#'                              x[2][1]
#'                            })) - 1
#'     )
#'
#'   colnames(sample) <- sample.name.from.data
#'   sample.information[, 1] <- sample.name.from.sample.information
#'
#'   subject.info <-
#'     sample.information[sample.information[, 3] == "Subject", ]
#'   if (hasQC == "no") {
#'     qc.info <- NULL
#'   } else {
#'     qc.info <- sample.information[sample.information[, 3] == "QC", ]
#'   }
#'
#'   if ((sum(is.na(subject)) + sum(is.na(qc))) == 0) {
#'     mv.imputation <- "yes"
#'     imputation.method <- "default"
#'   } else {
#'     mv.imputation <- "no"
#'     imputation.method <- "no"
#'   }
#'
#'   name <- tags[, "name"]
#'   if (polarity == "positive"){
#'     tags[, "name"] <- paste(name, "POS", sep = "_")
#'     pol <- rep("POS", nrow(tags))
#'     tags <- data.frame(pol, tags)
#'     colnames(tags)[1] <- "polarity"
#'   }
#'   if (polarity == "negative") {
#'     tags[, "name"] <- paste(name, "NEG", sep = "_")
#'     pol <- rep("NEG", nrow(tags))
#'     tags <- data.frame(pol, tags)
#'     colnames(tags)[1] <- "polarity"
#'   }
#'   #
#'
#'   MetFlowData <- new("MetFlowData",
#'                      subject = as.matrix(subject),
#'                      qc = as.matrix(qc),
#'                      tags = as.matrix(tags),
#'                      tags.old = as.matrix(tags),
#'                      subject.info = as.matrix(subject.info),
#'                      qc.info = as.matrix(qc.info),
#'                      subject.order = as.numeric(subject.order),
#'                      qc.order = as.numeric(qc.order),
#'                      ## some preprocessing information
#'                      mv.imputation = mv.imputation,
#'                      imputation.method = imputation.method,
#'                      zero.filter = "no",
#'                      zero.filter.criteria = "no",
#'                      normalization = "no",
#'                      normalization.method = "no",
#'                      data.integration = "no",
#'                      data.integration.method = "no",
#'                      hasIS = hasIS,
#'                      hasQC = hasQC,
#'                      peak.identification = peak.identification,
#'                      foldchange = "temp",
#'                      marker.selection.condition = "no",
#'                      mv.filter = "no",
#'                      mv.filter.criteria = "no",
#'                      univariate.test = "no",
#'                      qc.outlier.filter = "no",
#'                      subject.outlier.filter = "no")
#'
#'   MetFlowData <- MetFlowData
#'   }






#' @title MetabolitePlot
#' @description Give scatter plot for each feature.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData.before MetFlowData before normalization or integration.
#' @param MetFlowData.after MetFlowData after normalization or integration.
#' @param path Work directory.
#' @param figure Figure type you want to draw. jpeg ot pdf, default is jpeg.
#' @return Return the metabolite plot before and after processing.
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
#'MetabolitePlot(met.data.after.pre,
#'               met.data.after.pre,
#'               path = "Demo for MetabolitePlot")
#'               }

MetabolitePlot <- function(MetFlowData.before,
                           MetFlowData.after,
                           path = ".",
                           figure = "jpeg") {
  if (path != ".") {
    dir.create(path)
  }

  #
  qc_bef <- MetFlowData.before@qc
  subject_bef <- MetFlowData.before@subject
  tags_bef <- MetFlowData.before@tags
  data_bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject_bef1 <- data_bef[[1]]
  qc_bef1 <- data_bef[[2]]
  subject.info_bef1 <- data_bef[[3]]
  qc.info_bef1 <- data_bef[[4]]


  qc_aft <- MetFlowData.after@qc
  subject_aft <- MetFlowData.after@subject
  tags_aft <- MetFlowData.after@tags
  data_aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject_aft1 <- data_aft[[1]]
  qc_aft1 <- data_aft[[2]]
  subject.info_aft1 <- data_aft[[3]]
  qc.info_aft1 <- data_aft[[4]]

  subject.order1 <- as.numeric(MetFlowData.before@subject.order)
  qc.order1 <- as.numeric(MetFlowData.before@qc.order)

  subject.order2 <- as.numeric(MetFlowData.after@subject.order)
  qc.order2 <- as.numeric(MetFlowData.after@qc.order)
  feature.name <- as.character(tags_bef[, "name"])

  cat("Draw metabolite plot (%)\n")
  for (i in 1:nrow(subject_bef)) {
    if (figure == "jpeg")
    {
      jpeg(file.path(path, paste(feature.name[i], ".jpeg", sep = "")),
           width = 960,
           height = 480)
    }
    else
    {
      pdf(file.path(path, paste(feature.name[i], ".pdf", sep = "")),
          width = 12,
          height = 6)
    }
    layout(matrix(c(1:2), ncol = 2))
    par(mar = c(5, 5, 4, 2))
    ## before
    plot(
      subject.order1,
      as.numeric(subject_bef[i, ]),
      xlab = "Injection order",
      ylab = "Intensity",
      pch = 19,
      col = "royalblue",
      cex.lab = 1.5,
      cex.axis = 1.3,
      main = "Before"
    )
    points(qc.order1,
           as.numeric(qc_bef[i, ]),
           pch = 19,
           col = "firebrick1")
    ## add lines
    for (j in seq_along(qc.info_bef1)) {
      abline(v = max(qc.info_bef1[[j]][, 2]), lty = 2)
    }
    ## after
    plot(
      subject.order2,
      as.numeric(subject_aft[i, ]),
      xlab = "Injection order",
      ylab = "Intensity",
      pch = 19,
      col = "royalblue",
      cex.lab = 1.5,
      cex.axis = 1.3,
      main = "After"
    )
    points(qc.order2,
           as.numeric(qc_aft[i, ]),
           pch = 19,
           col = "firebrick1")
    ## add lines
    for (j in seq_along(qc.info_aft1)) {
      abline(v = max(qc.info_aft1[[j]][, 2]), lty = 2)
    }

    dev.off()
    ##
    count <- floor(nrow(subject_bef) * c(seq(0, 1, 0.01)))
    if (any(i == count)) {
      cat(ceiling(i * 100 / nrow(subject_bef)))
      cat(" ")
    }

  }

  layout(1)
}







#' @title ReChangeGroup
#' @description Change the group information in MetFlowData for statistical
#' analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param new.group New Group information name. Default is new.group.csv.
#' It can be use the sample.information which is only changed the group column.
#' @return Return a standard MetProcesser data which is changed the
#' group informatio.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "metflowR")
#' data(sample.information, package = "metflowR")
#' data(new.group, package = "metflowR")
#'
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#' write.csv(new.group, "new.group.csv", row.names = FALSE)
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
#' #run
#' new.met.data <- ReChangeGroup(met.data.after.pre)
#' }

ReChangeGroup <- function(MetFlowData,
                          new.group = "new.group.csv") {
  #
  new.group <-
    read.csv(new.group, stringsAsFactors = FALSE, check.names = FALSE)

  subject.info <- MetFlowData@subject.info
  subject <- MetFlowData@subject
  subject.order <- MetFlowData@subject.order

  ##which sample you want to remove from the dataset
  remove.name <- as.character(new.group[,1][which(is.na(new.group[,"group"]))])
  if(length(remove.name) != 0) {
    cat("The samples you want to remove from dataset are:\n")
    cat(remove.name)
    right <- readline("Right(y) or wrong(n)?")
    if (right == "n") {
      cat("Please change your new group information again.\n")
      return(MetFlowData)}
  }

  #

  ##remove the NA from new.group inforamtion
  new.group <- new.group[which(!is.na(new.group[,"group"])), ]
  new.subject.info <- new.group[new.group[, "class"] == "Subject", ]
  new.subject.name <- new.subject.info[, 1]

  ##remove samples from MetFlowData
  if (length(remove.name) != 0) {
    remove.idx <- match(remove.name, subject.info[,1])
    remove.idx <- remove.idx[!is.na(remove.idx)]
    subject <- subject[,-remove.idx]
    subject.info <- subject.info[-remove.idx,]
    subject.order <- subject.order[-remove.idx]
  }

  subject.name <- subject.info[, 1]

  ##change group information
  index <- match(subject.name, new.subject.name)
  new.subject.info <- new.subject.info[index, ]

  MetFlowData@subject.info <- as.matrix(new.subject.info)
  MetFlowData@subject <- as.matrix(subject)
  MetFlowData@subject.order <- as.numeric(subject.order)
  return(MetFlowData)
}











#' @title RSDfilter
#' @description Filter features according to RSD in QC samples.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param rsd.cutoff The cutoff value of RSD. Default is 30. It means that for
#' a feature, if its RSD larger than 30\%, it will be removed form the dataset.
#' @return Return a MetFlowData whose features have been removed.
#' @details \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:abe41368-d08d-4806-871f-3aa035d21743}{Dunn}
#' recommen the cutoff value of RSD should be set as 20\% in LC-MS
#' data and 30\% in GC-MS data.
#' @export
#' @examples
#' data(met.data, package = "metflowR")
#' new.met.data <- RSDfilter(met.data)

RSDfilter <- function(MetFlowData,
                      rsd.cutoff = 30) {
  ##RSD filtering
  qc <- MetFlowData@qc
  qc.rsd <- apply(qc, 1, function(x) {sd(x)*100/mean(x)})

  var.index <- which(qc.rsd <= rsd.cutoff)
  MetFlowData@qc <- MetFlowData@qc[var.index,]
  MetFlowData@tags <- MetFlowData@tags[var.index,]
  MetFlowData@subject <- MetFlowData@subject[var.index,]
  return(MetFlowData)
}


#S4 class
setClass("MetFlowData",
         representation(
           subject = "matrix",
           qc = "matrix",
           tags = "matrix",
           tags.old = "matrix",
           subject.info = "matrix",
           qc.info = "matrix",
           subject.order = "numeric",
           qc.order = "numeric",
           ## some preprocessing information
           mv.imputation = "character",
           imputation.method = "character",
           zero.filter = "character",
           zero.filter.criteria = "character",
           normalization = "character",
           normalization.method = "character",
           data.integration = "character",
           data.integration.method = "character",
           hasIS = "character",
           hasQC = "character",
           peak.identification = "character",
           foldchange = "character",
           marker.selection.condition = "character",
           mv.filter = "character",
           mv.filter.criteria = "character",
           univariate.test = "character",
           qc.outlier.filter = "character",
           subject.outlier.filter = "character"
         )
)

















#' @title RSDoverview
#' @description Evaluate the RSD information before and after processing.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData.before MetFlowData.before.
#' @param MetFlowData.after MetFlowData.after.
#' @param path Work directory.
#' @return RSD comparation.
#' @export
#' @examples
#' #load the demo data
#' data(met.data, package = "metflowR")
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'RSDoverview(met.data, met.data, path = "Demo for RSDoverview")

# RSDoverview(
#   MetFlowData.before = met.data.zero.filter,
#   MetFlowData.after = met.data,
#   path = file.path(path, "10 RSD overview")
# )

RSDoverview <- function(MetFlowData.before,
                        MetFlowData.after,
                        path = ".") {
  if (path != ".") {
    dir.create(path)
  }

  hasQC <- MetFlowData.before@hasQC
  if (hasQC == 'no') {
    stop("Data has no QC.")
  }
  ##before
  subject.bef <- MetFlowData.before@subject
  qc.bef <- MetFlowData.before@qc
  subject.info.bef <- MetFlowData.before@subject.info
  qc.info.bef <- MetFlowData.before@qc.info

  ##split data
  data.bef <- SplitBatch(MetFlowData = MetFlowData.before)
  subject.bef1 <- data.bef[[1]]
  qc.bef1 <- data.bef[[2]]

  ##after
  subject.aft <- MetFlowData.after@subject
  qc.aft <- MetFlowData.after@qc
  subject.info.aft <- MetFlowData.after@subject.info
  qc.info.aft <- MetFlowData.after@qc.info

  ##split data
  data.aft <- SplitBatch(MetFlowData = MetFlowData.after)
  subject.aft1 <- data.aft[[1]]
  qc.aft1 <- data.aft[[2]]


  pdf(
    file.path(path, "RSD distribution in different batch.pdf"),
    width = 14,
    height = 7
  )
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5, 5, 4, 2))
  for (i in seq_along(qc.bef1)) {
    #before
    QC <- qc.bef1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {
      sd(x) * 100 / mean(x)
    })
    par(mar = c(5, 5, 4, 2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd > 30] <- "tomato"
    colours1[is.na(colours1)] <- "grey"

    plot(
      QC.rsd,
      xlab = "Feature index",
      ylab = "Relative Standard Deviation (RSD, %)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colours1,
      main = paste("Batch", i, "QC (Before)"),
      cex.main = 1.3
    )
    abline(h = 30, lty = 2)
    legend("topleft",
           paste("RSD<30%: ", round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
           bty = "n",
           cex = 1.3)

    #after
    QC <- qc.aft1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {
      sd(x) * 100 / mean(x)
    })
    par(mar = c(5, 5, 4, 2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd > 30] <- "tomato"
    colours1[is.na(colours1)] <- "grey"

    plot(
      QC.rsd,
      xlab = "Feature index",
      ylab = "Relative Standard Deviation (RSD, %)",
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      col = colours1,
      main = paste("Batch", i, "QC (After)"),
      cex.main = 1.3
    )
    abline(h = 30, lty = 2)
    legend("topleft",
           paste("RSD<30%: ", round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
           bty = "n",
           cex = 1.3)

  }
  dev.off()

  pdf(
    file.path(path, "RSD distribution in all batches.pdf"),
    width = 14,
    height = 7
  )
  layout(matrix(c(1:2), ncol = 2))
  par(mar = c(5, 5, 4, 2))
  # before
  qc.rsd <- apply(qc.bef, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd > 30] <- "tomato"
  colours1[is.na(colours1)] <- "grey"
  plot(
    qc.rsd,
    xlab = "Feature index",
    ylab = "Relative Standard Deviation (RSD, %)",
    cex.lab = 1.3,
    cex.axis = 1.3,
    pch = 19,
    col = colours1,
    cex.main = 1.3,
    main = "Before"
  )
  abline(h = 30, lty = 2)
  legend("topleft",
         paste("RSD<30%: ", round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  # after
  qc.rsd <- apply(qc.aft, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd > 30] <- "tomato"
  colours1[is.na(colours1)] <- "grey"
  plot(
    qc.rsd,
    xlab = "Feature index",
    ylab = "Relative Standard Deviation (RSD, %)",
    cex.lab = 1.3,
    cex.axis = 1.3,
    pch = 19,
    col = colours1,
    cex.main = 1.3,
    main = "After"
  )
  abline(h = 30, lty = 2)
  legend("topleft",
         paste("RSD<30%: ", round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  dev.off()
  layout(1)
}
















#' @title SplitBatch
#' @description Split MetFlowData accoding to different batch.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @return Return a data (list), subject, qc, subject.info and qc.info.
#' @export
#' @examples
#' #load the demo data
#' data(met.data, package = "metflowR")
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'new.data <- SplitBatch(met.data)


SplitBatch <- function(MetFlowData) {
  ## split batch
  #
  hasQC <- MetFlowData@"hasQC"
  qc <- MetFlowData@"qc"
  subject <- MetFlowData@"subject"
  qc.info <- MetFlowData@"qc.info"
  subject.info <- MetFlowData@"subject.info"
  tags <- MetFlowData@"tags"

  subject.name <- as.character(subject.info[, 1])
  if (hasQC != "no") {
    qc.name <- as.character(qc.info[, 1])
  }
  else {
    qc.name <- NULL
  }

  subject.batch <- as.numeric(subject.info[, 4])
  if (hasQC != "no") {
    qc.batch <- as.numeric(qc.info[, 4])
  }
  else {
    qc.batch <- NULL
  }

  subject1 <- list()
  qc1 <- list()
  subject.info1 <- list()
  qc.info1 <- list()
  #
  for (i in seq_along(unique(subject.batch))) {
    subject.name.for.this.batch <- subject.name[subject.batch == i]
    if (hasQC == "no") {
      qc.name.for.this.batch <- NULL
    }
    else {
      qc.name.for.this.batch <- qc.name[qc.batch == i]
    }
    subject.index.for.this.batch <-
      match(subject.name.for.this.batch, colnames(subject))
    subject.index.for.this.batch <-
      subject.index.for.this.batch[!is.na(subject.index.for.this.batch)]
    if (hasQC == "no") {
      qc.index.for.this.batch <- NULL
    }
    else {
      qc.index.for.this.batch <-
        match(qc.name.for.this.batch, colnames(qc))
      qc.index.for.this.batch <-
        qc.index.for.this.batch[!is.na(qc.index.for.this.batch)]
    }

    subject1[[i]] <- subject[, subject.index.for.this.batch]
    subject.info1[[i]] <- subject.info[subject.batch == i, ]
    if (hasQC == "no") {
      qc1[[i]] <- NULL
      qc.info1[[i]] <- NULL
    }
    else {
      qc1[[i]] <- qc[, qc.index.for.this.batch]
      qc.info1[[i]] <- qc.info[qc.batch == i, ]
    }
  }

  data <- list(subject1, qc1, subject.info1, qc.info1)
  return(data)
}
