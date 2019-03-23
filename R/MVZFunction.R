#' @title MVimputation
#' @description Impute MV in data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param imputation.method Which imputation method you want to use? It
#' contains "knn", "rf" (missForest), "mean", "median", "zero", "minium",
#' "bpca" (BPCA), "svd" (SVD) and "ppca" (PPCA). Default is "knn".
#' The detial of
#' this method can be find in detail and reference paperes.
#' @param k See ?impute.knn
#' @param rowmax See ?impute.knn
#' @param colmax See ?impute.knn
#' @param maxp See ?impute.knn
#' @param rng.seed See ?impute.knn
#' @param maxiter See ?missForest
#' @param ntree See ?missForest
#' @param decreasing See ?missForest
#' @param replace See ?missForest
#' @param classwt See ?missForest
#' @param cutoff See ?missForest
#' @param strata See ?missForest
#' @param sampsize See ?missForest
#' @param nodesize See ?missForest
#' @param maxnodes See ?missForest
#' @param xtrue See ?missForest
#' @param parallelize See ?missForest
#' @param nPcs See ?bpca
#' @param maxSteps See ?bpca
#' @param threshold See ?bpca
#' @return Return a MetFlowData whose MVs have been imputated.
#' @seealso The MV imputation methods can see in
#' \code{\link[impute]{impute.knn}}, \code{\link[missForest]{missForest}},
#' \code{\link[pcaMethods]{bpca}}, \code{\link[pcaMethods]{ppca}} and
#' \code{\link[pcaMethods]{svdImpute}}.
#' @references The MV imputation in metabolomics data can see in
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:c9d05d0f-e945-43d0-bb4a-50ea0f90338e}{Guida's} paper.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "MetCleaning")
#' data(sample.information, package = "MetCleaning")
#'
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'#Import data
#'met.data <- ImportData(data = "data.csv",
#'                       sample.information = "sample.information.csv",
#'                       polarity = "positive")
#'#MV filtering
#'met.data <- MZfilter(MetFlowData = met.data,
#'                     obs.per.cutoff = 0.5,
#'                     var.per.cutoff = 0.5,
#'                     what = "mv",
#'                     path = "Demo for MV filter")
#'run
#'new.met.data <- MVimputation(met.data, rowmax = 0.9, colmax = 0.9)
#'}


MVimputation <- function(MetFlowData,
                         ##MV imputation method
                         imputation.method = "knn",
                         # knn parameters
                         k = 10,
                         rowmax = 0.5,
                         colmax = 0.8,
                         maxp = 1500,
                         rng.seed = 362436069,
                         # missForest parameters
                         maxiter = 10,
                         ntree = 100,
                         decreasing = FALSE,
                         replace = TRUE,
                         classwt = NULL,
                         cutoff = NULL,
                         strata = NULL,
                         sampsize = NULL,
                         nodesize = NULL,
                         maxnodes = NULL,
                         xtrue = NA,
                         parallelize = 'no',
                         #BPCA PPCA, and SVD parameters
                         nPcs = 2,
                         maxSteps = 100,
                         threshold = 1e-04) {
  options(warn = -1)
  #
  #### MV imputation
  if ((sum(is.na(MetFlowData@subject)) + sum(is.na(MetFlowData@qc))) == 0)
  {
    warning("MVs have been imputed.")
    return(MetFlowData)
  }
  qc <- MetFlowData@qc
  subject <- MetFlowData@subject
  qc.info <- MetFlowData@qc.info
  subject.info <- MetFlowData@subject.info
  tags <- MetFlowData@tags

  subject.name <- subject.info[, 1]
  qc.name <- qc.info[, 1]

  subject.batch <- subject.info[, 4]
  qc.batch <- qc.info[, 4]

  data <- SplitBatch(MetFlowData =  MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]

  # var.index <- list()
  for (i in seq_along(subject1)) {
    temp <- cbind(qc1[[i]], subject1[[i]])
    temp <- SXTMVimputation(
      data = temp,
      method = imputation.method,
      # knn parameters
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      maxp = maxp,
      rng.seed = rng.seed,
      # missForest parameters
      maxiter = maxiter,
      ntree = ntree,
      decreasing = decreasing,
      replace = replace,
      classwt = classwt,
      cutoff = cutoff,
      strata = strata,
      sampsize = sampsize,
      nodesize = nodesize,
      maxnodes = maxnodes,
      xtrue = xtrue,
      parallelize = parallelize,
      #BPCA PPCA, and SVD parameters
      nPcs = nPcs,
      maxSteps = maxSteps,
      threshold = threshold
    )
    qc1[[i]] <- temp[,(1:ncol(qc1[[i]]))]
    subject1[[i]] <- temp[,-c(1:ncol(qc1[[i]]))]
  }

  subject2 <- subject1[[1]]
  qc2 <- qc1[[1]]

  if (length(subject1) > 1) {
    for (i in 2:length(subject1)) {
      subject2 <- cbind(subject2, subject1[[i]])
      qc2 <- cbind(qc2, qc1[[i]])
    }
  }

  MetFlowData@subject <- as.matrix(subject2)
  MetFlowData@qc <- as.matrix(qc2)
  MetFlowData@mv.imputation <- "yes"
  MetFlowData@imputation.method <- imputation.method
  options(warn = 0)
  return(MetFlowData)
}



#' @title MZfilter
#' @description Filter feature and samples according to MV/zero ratio.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param var.per.cutoff The MV/zero ratio cutoff value of features,
#' default is 0.5. It means that for a feature, non-MV/zero ratio must be
#' larger than obs.per.cutoff.
#' @param obs.per.cutoff The MV/zero ratio cutoff value of samples,
#' default is 0.5. It means that for a sample, non-MV ratio/zero must be larger
#' than obs.per.cutoff.
#' @param what Filter missing values ("mv") or zero values ("zero")?
#' @param path Work directory.
#' @return Return a MetFlowData which has been filtered according MV/zero ratio.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "MetCleaning")
#' data(sample.information, package = "MetCleaning")
#'
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'#Import data
#'met.data <- ImportData(data = "data.csv",
#'                       sample.information = "sample.information.csv",
#'                       polarity = "positive")
#'#MV filtering
#'new.met.data <- MZfilter(MetFlowData = met.data,
#'                         obs.per.cutoff = 0.5,
#'                         var.per.cutoff = 0.5,
#'                         what = "mv",
#'                         path = "Demo for MV filter")
#'                         }


MZfilter <- function(MetFlowData,
                     obs.per.cutoff = 0.5,
                     var.per.cutoff = 0.5,
                     what = "mv",
                     path = ".") {
  if (path != ".") {
    dir.create(path)
  }
  options(warn = -1)

  MZfinderData <- MZfinder(
    MetFlowData = MetFlowData,
    obs.per.cutoff = obs.per.cutoff,
    var.per.cutoff = var.per.cutoff,
    what = what,
    path = path
  )

  MetFlowData <- MZfinderData[[1]]
  feature.remove <- MZfinderData[[2]]
  qc.remove <- MZfinderData[[3]]
  subject.remove <- MZfinderData[[4]]

  qc <- MetFlowData@qc
  subject <- MetFlowData@subject
  qc.info <- MetFlowData@qc.info
  subject.info <- MetFlowData@subject.info
  tags <- MetFlowData@tags
  hasQC <- MetFlowData@hasQC

  if (length(feature.remove) != 0) {
    if (hasQC == "yes")  {
      qc <- qc[-feature.remove, ]
    }
    subject <- subject[-feature.remove, ]
    tags <- tags[-feature.remove, ]
  }

  if (!is.null(qc.remove) & length(qc.remove) != 0) {
    cat("QC shoulde be removed are: \n")
    cat(qc.remove)
    cat("\n")
    cat("\n")
    qc.remove <-
      readline(
        "Which QC you want to remove(please type the index of QC sample,
        and separate them using comma,
        if you don't want to remove any QC, please type n):"
      )

    if (qc.remove == "n") {
      qc <- qc
    } else {
      qc.remove <- strsplit(qc.remove, split = ",")
      qc.remove <- as.numeric(qc.remove[[1]])
      qc.remove <- as.numeric(qc.remove)
      qc.remove.name <- colnames(qc)[qc.remove]
      qc <- qc[, -qc.remove]
    }
  }

  if (!is.null(subject.remove) & length(subject.remove) != 0) {
    cat("Subject shoulde be removed are: \n")
    cat(subject.remove)
    cat("\n")
    subject.remove <-
      readline(
        "Which subject you want to remove(
        please type the index of subject sample,
        and separate them using comma,
        if you don't want to remove any subject, please type n):"
      )
    if (subject.remove == "n") {
      subject <- subject
    }
    else {
      subject.remove <- strsplit(subject.remove, split = ",")
      subject.remove <- as.numeric(subject.remove[[1]])
      subject.remove <-
        as.numeric(subject.remove)
      subject <- subject[, -subject.remove]
    }
  }


  subject.name <- colnames(subject)
  if (hasQC == 'yes') {
    qc.name <- colnames(qc)
  }
  else {
    qc.name <- NULL
  }

  subject.index <-
    which(is.na(match(subject.info[, 1], subject.name)))
  if (length(subject.index) != 0) {
    subject.info <- subject.info[-subject.index,]
    MetFlowData@subject.order <-
      MetFlowData@subject.order[-subject.index]
  }

  if (hasQC == 'yes') {
    qc.index <- which(is.na(match(qc.info[, 1], qc.name)))
    if (length(qc.index) != 0) {
      qc.info <- qc.info[-qc.index,]
      MetFlowData@qc.order <-
        MetFlowData@qc.order[-qc.index]
    }
  }

  MetFlowData@subject.info <- as.matrix(subject.info)
  MetFlowData@qc.info <- as.matrix(qc.info)
  MetFlowData@subject <- as.matrix(subject)
  MetFlowData@qc <- as.matrix(qc)
  MetFlowData@tags <- as.matrix(tags)
  MetFlowData@mv.filter <- "yes"
  MetFlowData@mv.filter.criteria <-
    paste("varibale:",
          var.per.cutoff,
          "observation:",
          obs.per.cutoff)
  options(warn = 0)
  return(MetFlowData)
  }






##remove peaks whose MV ratio > threshold in QC or subject.
MZfinder <- function(MetFlowData,
                     obs.per.cutoff = 0.5,
                     var.per.cutoff = 0.5,
                     what = "mv",
                     path = ".") {
  #
  options(warn = -1)
  if (path != ".") {
    dir.create(path)
  }
  hasQC <- MetFlowData@hasQC
  if (what == "mv") {
    path1 <- file.path(path, "mv finder")
  } else {
    path1 <- file.path(path, "zero finder")
  }
  dir.create(path1)

  qc <- MetFlowData@qc
  subject <- MetFlowData@subject
  qc.info <- MetFlowData@qc.info
  subject.info <- MetFlowData@subject.info
  tags <- MetFlowData@tags
  subject.order <- as.numeric(MetFlowData@subject.order)
  qc.order <- as.numeric(MetFlowData@qc.order)

  subject.name <- subject.info[, 1]
  if (hasQC != "no") {
    qc.name <- qc.info[, 1]
  } else {
    qc.name <- NULL
  }

  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  var.index <- as.list(rep(NA, length(subject1)))
  if (hasQC == "yes") {
    ##remove peak whose MV/zero ratio more than 50%
    for (i in seq_along(subject1)) {
      if (what == "mv") {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          qc1[[i]],
          filter.item = "mv",
          filter.rule = "intersect",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      } else {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          qc1[[i]],
          filter.item = "zero",
          filter.rule = "intersect",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      }

      var.index[[i]] <- SXTMinifracData[["var.index"]]
    }
  } else {
    for (i in seq_along(subject1)) {
      if (what == "mv") {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          filter.item = "mv",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      } else {
        SXTMinifracData <- SXTMinifrac(
          subject1[[i]],
          filter.item = "zero",
          filter.rule = "intersect",
          minifrac.variable = var.per.cutoff,
          minifrac.observation = 0
        )
      }
      var.index[[i]] <- SXTMinifracData[["var.index"]]
    }
  }

  var.index1 <- var.index[[1]]
  if (length(var.index) > 1) {
    for (i in 2:length(var.index)) {
      var.index1 <- intersect(var.index1, var.index[[i]])
    }
  }

  subject1 <- lapply(subject1, function(x) {
    x[var.index1,]
  })

  if (hasQC == "yes") {
    qc1 <- lapply(qc1, function(x) {
      x[var.index1,]
    })
  } else {
    qc1 <- NULL
  }

  if (hasQC == "yes") {
    subject2 <- subject1[[1]]
    qc2 <- qc1[[1]]

    if (length(var.index) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
        qc2 <- cbind(qc2, qc1[[i]])
      }
    }
  } else {
    subject2 <- subject1[[1]]
    if (length(var.index) > 1) {
      for (i in 2:length(subject1)) {
        subject2 <- cbind(subject2, subject1[[i]])
      }
    }
    qc <- NULL
  }

  ##remove sample whose MV/zero more than 50%
  if (hasQC == "yes") {
    q.sample.per <- apply(qc2, 2, function(x) {
      ifelse(what == "mv", sum(!is.na(x)), sum(x != 0)) / nrow(qc2)
    })
    qc.index <- which(q.sample.per > obs.per.cutoff)
  } else {
    qc.index <- NULL
  }
  s.sample.per <-
    apply(subject2, 2, function(x) {
      ifelse(what == "mv", sum(!is.na(x)), sum(x != 0)) / nrow(subject2)
    })
  subject.index <- which(s.sample.per > obs.per.cutoff)

  ## output some information about feature and samples
  feature.remove <- setdiff(c(1:nrow(subject)), var.index1)
  if (hasQC == "yes") {
    qc.remove <- setdiff(c(1:ncol(qc)), qc.index)
  } else {
    qc.remove <- NULL
  }
  subject.remove <- setdiff(c(1:ncol(subject)), subject.index)

  if (length(feature.remove) == 0) {
    cat("No features should be removed.\n")
  } else {
    write.csv(tags[feature.remove,],
              file.path(path1, "Feature to remove information.csv"))
  }

  if (hasQC != "no") {
    if (length(qc.remove) == 0) {
      cat("No QC should be removed.\n")
    } else {
      cat(colnames(qc2)[qc.remove], "sholud be removed.\n")
    }
  }

  if (length(subject.remove) == 0) {
    cat("No subject should be removed.\n")
  } else {
    cat(colnames(subject2)[subject.remove], "sholud be removed.\n")
  }
  options(warn = 0)
  MZfinderData <- list(
    MetFlowData = MetFlowData,
    feature.remove = feature.remove,
    qc.remove = qc.remove,
    subject.remove = subject.remove
  )
  class(MZfinderData) <- "MZfinderData"
  return(MZfinderData)
}





#' @title MZOverview
#' @description Evaluate the missing values or zero value information
#' in MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path Work Directory.
#' @param what Missing value ("mv") or zero values ("zero").
#' @return Batch i feature MV/zero ratio.csv: Batch i feature MV/zero ratio
#' distribution.
#' @return Batch i QC feature MV/zero ratio.pdf: Batch i QC feature MV.zero
#' ratio distribution.
#' @return Batch i Subject feature MV/zero ratio.pdf: Batch i Subject feature
#'  MV/zero ratio distribution.
#' @return Sample MV/zero distribution.csv: Sample MV/zero ratio distribution.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "MetCleaning")
#' data(sample.information, package = "MetCleaning")
#'
#' ##create a folder for demo
#' dir.create("demo")
#' setwd("demo")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'#Import data
#'met.data <- ImportData(data = "data.csv",
#'                       sample.information = "sample.information.csv",
#'                       polarity = "positive")
#'#run
#'MZoverview(MetFlowData = met.data,
#'                     what = "mv",
#'                     path = "Demo for MZoverview")
#'                     }

MZoverview <- function(MetFlowData,
                       path = ".",
                       what = "mv") {

  has.qc <- MetFlowData@hasQC

  if (path != ".") {
    dir.create(path)
  }


  subject <- MetFlowData@subject
  qc <- MetFlowData@qc

  if(is.null(qc)) {has.qc <- "no"}

  tags <- MetFlowData@tags
  subject.order <- as.numeric(MetFlowData@subject.order)

  if(has.qc == "no") {
    qc.order <- NULL
  }else{
    qc.order <- as.numeric(MetFlowData@qc.order)
  }

  if (what == "zero" & (sum(is.na(qc)) + sum(is.na(subject))) != 0){
    stop("Please impute MV first.")
  }

  ##variable for MV record
  data <- SplitBatch(MetFlowData = MetFlowData)
  subject1 <- data[[1]]
  qc1 <- data[[2]]
  subject.info1 <- data[[3]]
  qc.info1 <- data[[4]]

  if(has.qc == "yes") {
    s.number.all <- NULL
    s.per.all <- NULL
    q.number.all <- NULL
    q.per.all <- NULL

    b <- list()
    for (i in seq_along(subject1)) {
      temp.subject <- subject1[[i]]
      temp.qc <- qc1[[i]]

      if (what == "mv") {
        s.number <- sum(is.na(temp.subject))
        q.number <- sum(is.na(temp.qc))
      } else {
        s.number <- sum(temp.subject == 0)
        q.number <- sum(temp.qc == 0)
      }

      s.per <-
        round(s.number * 100 / (nrow(temp.subject) * ncol(temp.subject)), 2)
      q.per <-
        round(q.number * 100 / (nrow(temp.qc) * ncol(temp.qc)), 2)

      ## record
      s.number.all[i] <- s.number
      s.per.all[i] <- s.per

      q.number.all[i] <- q.number
      q.per.all[i] <- q.per

      #ratio in each feature and sample
      if (what == "mv") {
        s.feature.per <-
          apply(subject1[[i]], 1, function(x) {
            sum(is.na(x)) * 100 / ncol(subject1[[i]])
          })
        s.sample.per <-
          apply(subject1[[i]], 2, function(x) {
            sum(is.na(x)) * 100 / nrow(subject1[[i]])
          })

        q.feature.per <-
          apply(qc1[[i]], 1, function(x) {
            sum(is.na(x)) * 100 / ncol(qc1[[i]])
          })
        q.sample.per <-
          apply(qc1[[i]], 2, function(x) {
            sum(is.na(x)) * 100 / nrow(qc1[[i]])
          })
      } else {
        s.feature.per <-
          apply(subject1[[i]], 1, function(x) {
            sum(x == 0) * 100 / ncol(subject1[[i]])
          })
        s.sample.per <-
          apply(subject1[[i]], 2, function(x) {
            sum(x == 0) * 100 / nrow(subject1[[i]])
          })

        q.feature.per <-
          apply(qc1[[i]], 1, function(x) {
            sum(x == 0) * 100 / ncol(qc1[[i]])
          })
        q.sample.per <-
          apply(qc1[[i]], 2, function(x) {
            sum(x == 0) * 100 / nrow(qc1[[i]])
          })
      }

      ## plot
      pdf(file.path(
        path,
        paste(
          "Batch",
          i,
          "feature",
          ifelse(what == "mv", "MV", "zero") ,
          "distribution.pdf"
        )
      ), width = 14)
      layout(matrix(c(1:2), ncol = 2))
      par(mar = c(5, 5, 4, 2))
      #subject
      plot(
        s.feature.per,
        xlab = "Feature index",
        ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        pch = 19,
        main = paste(
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (Subject)"
        ),
        cex.main = 1.3
      )
      abline(
        h = 30,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )
      #QC
      plot(
        q.feature.per,
        xlab = "Feature index",
        ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        pch = 19,
        main = paste(
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (QC)"
        ),
        cex.main = 1.3
      )
      abline(
        h = 30,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )

      #subject
      hist(
        s.feature.per,
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        main = paste(
          "Histogram of",
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (Subject)"
        ),
        cex.main = 1.3
      )
      #QC
      hist(
        q.feature.per,
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        main = paste(
          "Histogram of",
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (QC)"
        ),
        cex.main = 1.3
      )

      #subject
      plot(
        x = sort(s.feature.per),
        y = c(1:length(s.feature.per)) * 100 / length(s.feature.per),
        type = "l",
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        ylab = "Cumulative feature percentage (%)",
        cex.lab = 1.3,
        cex.axis = 1.3,
        lwd = 2,
        main = paste(
          "Cumulative",
          ifelse(what == "mv", "MV", "zero"),
          "ratio (Subject)"
        )
      )
      a <-
        round(sum(s.feature.per < 50) * 100 / length(s.feature.per), 2)
      abline(
        v = 50,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )
      legend(
        "topleft",
        legend = paste(a, "%"),
        title = paste(ifelse(what == "mv", "MV", "zero"), "ratio < 30%"),
        bty = "n"
      )
      #QC
      plot(
        x = sort(q.feature.per),
        y = c(1:length(q.feature.per)) * 100 / length(q.feature.per),
        type = "l",
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        ylab = "Cumulative feature percentage (%)",
        cex.lab = 1.3,
        cex.axis = 1.3,
        lwd = 2,
        main = paste(
          "Cumulative",
          ifelse(what == "mv", "MV", "zero"),
          "ratio (QC)"
        )
      )
      a <-
        round(sum(q.feature.per < 50) * 100 / length(q.feature.per), 2)
      abline(
        v = 50,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )
      legend(
        "topleft",
        legend = paste(a, "%"),
        title = paste(ifelse(what == "mv", "MV", "zero"), "ratio < 30%"),
        bty = "n"
      )
      dev.off()

      b[[i]] <- data.frame(s.feature.per, q.feature.per)
      colnames(b[[i]]) <-
        c(
          paste(
            "Batch",
            i,
            "Subject",
            ifelse(what == "mv", "MV", "zero"),
            "MV ratio"
          ),
          paste(
            "Batch",
            i,
            "QC",
            ifelse(what == "mv", "MV", "zero"),
            "MV ratio"
          )
        )
    }

    if (length(b) > 1) {
      b1 <- b[[1]]
      for (i in 2:length(b)) {
        b1 <- cbind(b1, b[[i]])
      }
    } else{
      b1 <- b[[1]]
    }
    b1 <- cbind(tags[, "name"], b1)
    colnames(b1)[1] <- "Feature name"
    write.csv(b1,
              file.path(path,
                        paste(
                          "feature",
                          ifelse(what == "mv", "MV", "zero") , "ratio.csv"
                        )))

    ## whole data
    if (what == "mv") {
      q.sample.ratio <-
        apply(qc, 2, function(x) {
          sum(is.na(x)) * 100 / nrow(qc)
        })
      s.sample.ratio <-
        apply(subject, 2, function(x) {
          sum(is.na(x)) * 100 / nrow(subject)
        })
    } else {
      q.sample.ratio <-
        apply(qc, 2, function(x) {
          sum(x == 0) * 100 / nrow(qc)
        })
      s.sample.ratio <-
        apply(subject, 2, function(x) {
          sum(x == 0) * 100 / nrow(subject)
        })
    }
    pdf(file.path(path, paste("Sample", ifelse(what == "mv", "MV", "zero") ,
                              "distribution.pdf")))
    par(mar = c(5, 5, 4, 2))
    plot(
      subject.order,
      s.sample.ratio,
      xlab = "Injection order",
      ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      main = paste(ifelse(what == "mv", "MV", "zero"), "ratio distribution"),
      cex.main = 1.3,
      col = "royalblue",
      ylim = c(0, max(c(
        s.sample.ratio, q.sample.ratio
      ))),
      xlim = c(1, max(qc.order))
    )

    abline(
      h = 50,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )

    #add text
    idx <- which(s.sample.ratio >= 50)
    if (length(idx) >= 1) {
      text(
        x = subject.order[idx],
        y = s.sample.ratio[idx],
        labels = colnames(subject)[idx],
        pos = 1
      )
    }
    points(qc.order, q.sample.ratio, pch = 19, col = "firebrick1")

    #add text
    idx <- which(q.sample.ratio >= 50)
    if (length(idx) >= 1) {
      text(
        x = qc.order[idx],
        y = q.sample.ratio[idx],
        labels = colnames(qc)[idx],
        pos = 1
      )
    }
    v.index <- qc.order[which(diff(qc.order) == 1)]
    for (i in v.index) {
      abline(v = i, lty = 2, lwd = 2)
    }
    dev.off()
    c <- data.frame(c(q.sample.ratio, s.sample.ratio))
    colnames(c) <- paste(ifelse(what == "mv", "MV", "zero"), "ratio")

    write.csv(c, file.path(path,
                           paste(
                             "Sample",
                             ifelse(what == "mv", "MV", "zero"),
                             "distribution.csv"
                           )))

    ## MV information in each batch
    batch.info <-
      data.frame(s.number.all,
                 s.per.all,
                 q.number.all,
                 q.per.all)
    colnames(batch.info) <-
      c(
        paste("Subject", ifelse(what == "mv", "MV", "zero"), "number"),
        paste("Subject", ifelse(what == "mv", "MV", "zero"), "percentage(%)"),
        paste("QC", ifelse(what == "mv", "MV", "zero"), "number"),
        paste("QC", ifelse(what == "mv", "MV", "zero"), "percentage(%)")
      )
    write.csv(batch.info,
              file.path(path,
                        paste(
                          ifelse(what == "mv", "MV", "zero"),
                          "information in each batch.csv"
                        )))
  }else{##if has no QC
    s.number.all <- NULL
    s.per.all <- NULL

    b <- list()
    for (i in seq_along(subject1)) {
      temp.subject <- subject1[[i]]

      if(what == "mv") {
        s.number <- sum(is.na(temp.subject))
      } else {
        s.number <- sum(temp.subject == 0)
      }

      s.per <-
        round(s.number * 100 / (nrow(temp.subject) * ncol(temp.subject)), 2)

      ## record
      s.number.all[i] <- s.number
      s.per.all[i] <- s.per

      #ratio in each feature and sample
      if (what == "mv") {
        s.feature.per <-
          apply(subject1[[i]], 1, function(x) {
            sum(is.na(x)) * 100 / ncol(subject1[[i]])
          })
        s.sample.per <-
          apply(subject1[[i]], 2, function(x) {
            sum(is.na(x)) * 100 / nrow(subject1[[i]])
          })

      } else {
        s.feature.per <-
          apply(subject1[[i]], 1, function(x) {
            sum(x == 0) * 100 / ncol(subject1[[i]])
          })
        s.sample.per <-
          apply(subject1[[i]], 2, function(x) {
            sum(x == 0) * 100 / nrow(subject1[[i]])
          })
      }

      ## plot
      pdf(file.path(
        path,
        paste(
          "Batch",
          i,
          "feature",
          ifelse(what == "mv", "MV", "zero") ,
          "distribution.pdf"
        )), width = 7)
      par(mar = c(5, 5, 4, 2))
      #subject
      plot(
        s.feature.per,
        xlab = "Feature index",
        ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        pch = 19,
        main = paste(
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (Subject)"
        ),
        cex.main = 1.3
      )
      abline(
        h = 30,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )

      #subject
      hist(
        s.feature.per,
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        cex.lab = 1.3,
        cex.axis = 1.3,
        main = paste(
          "Histogram of",
          ifelse(what == "mv", "MV", "zero"),
          "ratio distribution (Subject)"
        ),
        cex.main = 1.3
      )

      #subject
      plot(
        x = sort(s.feature.per),
        y = c(seq_along(s.feature.per)) * 100 / length(s.feature.per),
        type = "l",
        xlab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
        ylab = "Cumulative feature percentage (%)",
        cex.lab = 1.3,
        cex.axis = 1.3,
        lwd = 2,
        main = paste(
          "Cumulative",
          ifelse(what == "mv", "MV", "zero"),
          "ratio (Subject)"
        )
      )
      a <-
        round(sum(s.feature.per < 50) * 100 / length(s.feature.per), 2)
      abline(
        v = 50,
        lty = 2,
        col = "firebrick1",
        lwd = 2
      )
      legend(
        "topleft",
        legend = paste(a, "%"),
        title = paste(ifelse(what == "mv", "MV", "zero"), "ratio < 30%"),
        bty = "n"
      )
      dev.off()

      # b[[i]] <- data.frame(s.feature.per, q.feature.per)
      # colnames(b[[i]]) <-
      #   c(
      #     paste(
      #       "Batch",
      #       i,
      #       "Subject",
      #       ifelse(what == "mv", "MV", "zero"),
      #       "MV ratio"
      #     ),
      #     paste(
      #       "Batch",
      #       i,
      #       "QC",
      #       ifelse(what == "mv", "MV", "zero"),
      #       "MV ratio"
      #     )
      #   )
    }

    # if (length(b) > 1) {
    #   b1 <- b[[1]]
    #   for (i in 2:length(b)) {
    #     b1 <- cbind(b1, b[[i]])
    #   }
    # } else{
    #   b1 <- b[[1]]
    # }
    # b1 <- cbind(tags[, "name"], b1)
    # colnames(b1)[1] <- "Feature name"
    # write.csv(b1,
    #           file.path(path,
    #                     paste(
    #                       "feature",
    #                       ifelse(what == "mv", "MV", "zero") , "ratio.csv"
    #                     )))

    ## whole data
    if (what == "mv") {
      s.sample.ratio <-
        apply(subject, 2, function(x) {
          sum(is.na(x)) * 100 / nrow(subject)
        })
    } else {
      s.sample.ratio <-
        apply(subject, 2, function(x) {
          sum(x == 0) * 100 / nrow(subject)
        })
    }
    pdf(file.path(path, paste("Sample", ifelse(what == "mv", "MV", "zero"),
                              "distribution.pdf")))
    par(mar = c(5, 5, 4, 2))
    plot(
      subject.order,
      s.sample.ratio,
      xlab = "Injection order",
      ylab = paste(ifelse(what == "mv", "MV", "zero"), "ratio (%)"),
      cex.lab = 1.3,
      cex.axis = 1.3,
      pch = 19,
      main = paste(ifelse(what == "mv", "MV", "zero"), "ratio distribution"),
      cex.main = 1.3,
      col = "royalblue",
      ylim = c(0, max(c(s.sample.ratio))))

    abline(
      h = 50,
      lty = 2,
      col = "firebrick1",
      lwd = 2
    )

    #add text
    idx <- which(s.sample.ratio >= 50)
    if (length(idx) >= 1) {
      text(
        x = subject.order[idx],
        y = s.sample.ratio[idx],
        labels = colnames(subject)[idx],
        pos = 1
      )
    }
    dev.off()
    #
    # c <- data.frame(c(q.sample.ratio, s.sample.ratio))
    # colnames(c) <- paste(ifelse(what == "mv", "MV", "zero"), "ratio")

    # write.csv(c, file.path(path,
    #                        paste(
    #                          "Sample",
    #                          ifelse(what == "mv", "MV", "zero"),
    #                          "distribution.csv"
    #                        )))

    ## MV information in each batch
    # batch.info <-
    #   data.frame(s.number.all,
    #              s.per.all,
    #              q.number.all,
    #              q.per.all)
    # colnames(batch.info) <-
    #   c(
    #     paste("Subject", ifelse(what == "mv", "MV", "zero"), "number"),
    #     paste("Subject", ifelse(what == "mv", "MV", "zero"), "percentage(%)"),
    #     paste("QC", ifelse(what == "mv", "MV", "zero"), "number"),
    #     paste("QC", ifelse(what == "mv", "MV", "zero"), "percentage(%)")
    #   )
    # write.csv(batch.info,
    #           file.path(path,
    #                     paste(
    #                       ifelse(what == "mv", "MV", "zero"),
    #                       "information in each batch.csv"
    #                     )))



  }

}








SXTMVimputation <- function(data,
                            method = "knn",
                            ## knn parameters
                            k = 10,
                            rowmax = 0.5,
                            colmax = 0.8,
                            maxp = 1500,
                            rng.seed = 362436069,
                            ## missForest parameters
                            maxiter = 10,
                            ntree = 100,
                            decreasing =FALSE,
                            mtry = floor(sqrt(nrow(data))),
                            replace = TRUE,
                            classwt = NULL,
                            cutoff = NULL,
                            strata = NULL,
                            sampsize = NULL,
                            nodesize = NULL,
                            maxnodes = NULL,
                            xtrue = NA,
                            parallelize = 'no',
                            ##BPCA PPCA, and SVD parameters
                            nPcs = 2,
                            maxSteps = 100,
                            threshold = 1e-04
) {
  #

  ## KNN method
  if (method == "knn") {
    # library(impute)
    if(exists(".Random.seed")) rm(.Random.seed)
    data.knn <- impute::impute.knn(as.matrix(data),
                                   k = k,
                                   rowmax = rowmax,
                                   colmax = colmax,
                                   maxp = maxp,
                                   rng.seed = rng.seed)
    data.knn <- data.knn[["data"]]
    return(data.knn)
  }

  #missForest
  if (method=="rf") {
    # library(missForest)
    data.rf <- missForest::missForest(t(data),
                                      maxiter = maxiter,
                                      ntree = ntree,
                                      decreasing = decreasing,
                                      mtry = mtry,
                                      replace = replace,
                                      classwt = classwt,
                                      cutoff = cutoff,
                                      strata = strata,
                                      sampsize = sampsize,
                                      nodesize = nodesize,
                                      maxnodes = maxnodes,
                                      xtrue = xtrue,
                                      parallelize = 'no')
    data.rf <- t(data.rf$ximp)
    return(data.rf)
  }

  ## mean imputation
  if (method == "mean") {
    data.mean <-
      apply(data,1,function(x) {x <- ifelse (is.na(x), mean(x, na.rm = TRUE), x)})
    return(t(data.mean))
  }

  ## median imputation
  if (method == "median") {
    data.median <-
      apply(data,1,function(x) {x <- ifelse (is.na(x), median(x, na.rm = TRUE), x)})
    return(t(data.median))
  }

  ## zero imputation
  if (method == "zero") {
    data.zero <- data
    data.zero[is.na(data.zero)] <- 0
    return(data.zero)
  }

  ## minimum imputation
  if (method == "minimum") {
    data.minimum <-
      apply(data,1,function(x) {x <- ifelse (is.na(x), min(x, na.rm = TRUE), x)})
    return(t(data.minimum))
  }

  ##BPCA
  if (method == "bpca") {
    data.bpca <- pcaMethods::pca(t(data),
                                 method="bpca",
                                 nPcs = nPcs,
                                 maxSteps = maxSteps,
                                 threshold = threshold)
    data.bpca <- t(pcaMethods::completeObs(data.bpca))
    return(data.bpca)
  }

  ##SVD imputation
  if (method == "svd") {
    data.svd <- pcaMethods::pca(t(data),
                                method="svdImpute",
                                nPcs=nPcs)
    data.svd <- t(pcaMethods::completeObs(data.svd))
    return(data.svd)
  }

  ##PPCA imputation
  if (method == "ppca") {
    data.ppca <- pcaMethods::pca(t(data),
                                 method = "ppca",
                                 nPcs = nPcs)
    data.ppca <- t(pcaMethods::completeObs(data.ppca))
    return(data.ppca)
  }

}