#' @title metflowR
#' @description A whole work flow for high throughput
#' MS based metabolomics data cleaning.
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
#' @param combine.mz.tol Batch alignment mz tolerance.
#' @param combine.rt.tol Batch alignment RT tolerance.
#' @param mv.filter Filter missing values or not? Default is TRUE.
#' @param zero.filter Filter zero values or not? Default is TRUE.
#' @param normalization Data normalization or not? Default is TRUE.
#' @param integration Data integration or not? Default is TRUE.
#' @param qc.outlier.filter Filter QC outliers or not? Default is TRUE.
#' @param subject.outlier.filter Filter subject outliers or not?
#' Default is TRUE.
#' @param obs.mv.cutoff The observation MV ratio cutoff.
#' @param var.mv.cutoff The variable MV ratio cutoff.
#' @param obs.zero.cutoff The observation zero ratio cutoff.
#' @param var.zero.cutoff The variable zero ratio cutoff.
#' @param imputation.method Which imputation method you want to use? It
#' contains "knn", "rf" (missForest), "mean", "median", "zero", "minium",
#' "bpca" (BPCA), "svd" (SVD) and "ppca" (PPCA). Default is "knn".
#' The detial of this method can be find in detail and reference paperes.
#' @param k See ?impute.knn
#' @param rowmax See ?impute.knn
#' @param colmax See ?impute.knn
#' @param maxp See ?impute.knn
#' @param method Normalization method, mean, median, total svr or loess,
#' default is svr. Please see the details.
#' @param multiple see ?SXTsvrNor.
#' @param threads Thread number.
#' @param hmdb.matching Default is FALSE.
#' @param mass.tolerance mz tolerance.
#' @param show Default is 5.
#' @param mz.tolerance mz tolerance for ms1 and ms2 data matching.
#' @param rt.tolerance RT tolerance for ms1 and ms2 data matching.
#' @param met.plot Scatter of peak.
#' @param path Work directory.
#' @return  All the results can be got form other functions and instruction.
#' @export
#' @details The manual of metflowR can be found in \href{https://github.com/jaspershen/metflowR/blob/master/vignettes/metflowR.pdf}{github}
#' or in \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:d8b5fcff-c725-4689-a97d-7ff106322fb6}{my library}.
#' @examples
#' \donttest{
#' #load the demo data
#' data(data, package = "metflowR")
#' data(sample.information, package = "metflowR")
#'
#' ##create a folder for metflowR demo
#' dir.create("Demo for metflowR")
#' setwd("Demo for metflowR")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'
#' #run metflowR
#' metflowR(#ImportData para
#' data = "data.csv",
#' sample.information = "sample.information.csv",
#' polarity = "positive",
#' #DataNormalization
#' method = "svr",
#' threads = 2)
#' }
#
# setwd("D:/study/R/my git/metCube/version1.1/user_data/test/test")
#
# one.user.path <- file.path(getwd(), "one_data_cleaning")
# load("one_data_cleaning/temp.file")
# load("one_data_cleaning/oneDCparameter")
#
#
# metflowR.info <- try(expr = {metflowR(
#   data = grep("batch", x = temp.file, value = TRUE),
#   sample.information = "sample.info.csv",
#   polarity = oneDCparameter[[1]],
#   hasIS = "no",
#   hasQC = "yes",
#   #batch alignment
#   combine.mz.tol = oneDCparameter[[2]],
#   combine.rt.tol = oneDCparameter[[3]],
#   #MVFilter para
#   mv.filter = TRUE,
#   obs.mv.cutoff = 0,
#   var.mv.cutoff = 1 - oneDCparameter[[3]]/100,
#   #MVimputation
#   imputation.method = oneDCparameter[[5]],
#   k = oneDCparameter[[6]],
#   rowmax = oneDCparameter[[7]],
#   colmax = oneDCparameter[[8]],
#   maxp = 1500,
#   #ZeroFilter para
#   zero.filter = TRUE,
#   obs.zero.cutoff = 0,
#   var.zero.cutoff = 1 - oneDCparameter[[12]]/100,
#   #DataNormalization
#   normalization = TRUE,
#   method = oneDCparameter[[14]],
#   multiple = oneDCparameter[[19]],
#   threads = 2,
#   #PeakIdentification
#   hmdb.matching = FALSE,
#   show = 5,
#   mass.tolerance = 30,
#   mz.tolerance = oneDCparameter[[29]],
#   rt.tolerance = oneDCparameter[[30]],
#   #DataOverview para
#   met.plot = FALSE,
#   path = one.user.path,
#   # worklist.from = "manual",
#   #other slection
#   qc.outlier.filter = TRUE,
#   subject.outlier.filter = TRUE,
#   integration = TRUE)},silent = FALSE)




setGeneric(name = "metflowR",
           def = function(#ImportData para
             data = c("batch1.demo.csv","batch2.demo.csv"),
             sample.information = "sample.info.demo.csv",
             polarity = c("positive", "negative", "none"),
             hasIS = c("no", "yes"),
             hasQC = c("yes", "no"),
             #batch alignment
             combine.mz.tol = 25,
             combine.rt.tol = 30,
             #MVFilter para
             mv.filter = TRUE,
             obs.mv.cutoff = 0.5,
             var.mv.cutoff = 0.5,
             #MVimputation
             imputation.method = "knn",
             k = 10,
             rowmax = 0.5,
             colmax = 0.5,
             maxp = 1500,
             #ZeroFilter para
             zero.filter = TRUE,
             obs.zero.cutoff = 0,
             var.zero.cutoff = 0.5,
             #DataNormalization
             normalization = TRUE,
             method = c("svr", "loess"),
             multiple = 5,
             threads = 2,
             #PeakIdentification
             hmdb.matching = FALSE,
             show = 5,
             mass.tolerance = 30,
             mz.tolerance = 30,
             rt.tolerance = 180,
             #DataOverview para
             met.plot = FALSE,
             path = ".",
             # worklist.from = "manual",
             #other slection
             qc.outlier.filter = TRUE,
             subject.outlier.filter = TRUE,
             integration = TRUE){

             polarity <- match.arg(polarity)
             hasIS <- match.arg(hasIS)
             hasQC <- match.arg(hasQC)
             method <- match.arg(method)


             if (path != ".") {
               dir.create(path)
             }

             path.inter <- file.path(path, "intermediate")
             dir.create(path.inter)

             options(warn = -1)
             #--------------------------------------------------------------------------
             #check data
             cat("\n")
             cat("Check data...\n")
             stat <- checkData(data = data,
                               sample.info = sample.information, path = path)
             if(any(as.numeric(stat[,4]) > 0)){
               stop("Error is in you data.")
             }


             #-------------------------------------------------------------------------
             #batch alignment
             cat("\n")
             cat("Batch alignment...\n")
             if (all(dir(path.inter) != "met.data.raw")) {
              ##read data first
              peak.table <- lapply(data, function(x){
                as.data.frame(readr::read_csv(file.path(path, x),
                                              col_types = readr::cols()))
              })

              roughMatchResult <- roughAlign(peak.table = peak.table,
                                             combine.mz.tol = combine.mz.tol,
                                             combine.rt.tol = combine.rt.tol)

              accurateMatchResult <- accurateAlign(peak.table = peak.table,
                                             simple.data = roughMatchResult)

              ##output data
              readr::write_csv(accurateMatchResult, file.path(path, "peak.table.after.batch.alignment.csv"))
             }else{
               cat("Using previous result.\n")
             }

             #--------------------------------------------------------------------------
             #read data
             cat("\n")
             cat("Importing data...\n")

             if (all(dir(path.inter) != "met.data.raw")) {
               met.data <- ImportData(
                 data = 'peak.table.after.batch.alignment.csv',
                 sample.information = sample.information,
                 polarity = polarity,
                 hasIS = hasIS,
                 hasQC = hasQC, path = path)
               #save data
               met.data.raw <- met.data
               save(met.data.raw, file = file.path(path.inter, "met.data.raw"))
             } else {
               load(file.path(path.inter, "met.data.raw"))
               met.data <- met.data.raw
             }

             batch <- unique(met.data@subject.info[, 4])
             subject <- met.data@subject

             #-------------------------------------------------------------------------
             if (sum(is.na(subject)) != 0) {
               MZoverview(
                 MetFlowData = met.data,
                 path = file.path(path, "1 MV overview"),
                 what = "mv"
               )

               if (mv.filter) {
                 cat(
                   "------------------------------------------------------------------\n"
                 )
                 cat("Missing values filter...\n")
                 #filter mv
                 met.data <- MZfilter(
                   MetFlowData = met.data,
                   obs.per.cutoff = obs.mv.cutoff,
                   var.per.cutoff = var.mv.cutoff,
                   what = "mv",
                   path = file.path(path, "2 MV filter")
                 )
                 #save data
                 met.data.mv.filter <- met.data
                 save(met.data.mv.filter,
                      file = file.path(path.inter, "met.data.mv.filter"))
               }

               cat("---------------------------------------------------------------\n")
               cat("Missing values imputation...\n")
               #mv imputation
               met.data <- MVimputation(
                 MetFlowData = met.data,
                 ##MV imputation method
                 imputation.method = imputation.method,
                 # knn parameters
                 k = k,
                 rowmax = rowmax,
                 colmax = colmax,
                 maxp = maxp
               )
               #save data
               met.data.mv.imputation <- met.data
               save(met.data.mv.imputation,
                    file = file.path(path.inter, "met.data.mv.imputation"))
             }

             subject <- met.data@subject

             #zero distribution
             MZoverview(
               MetFlowData = met.data,
               what = "zero",
               path = file.path(path, "3 Zero overview")
             )

             if (zero.filter) {
               cat("------------------------------------------------------------------\n")
               cat("Zero filter...\n")
               #filter zero
               met.data <- MZfilter(
                 MetFlowData = met.data,
                 obs.per.cutoff = obs.zero.cutoff,
                 var.per.cutoff = var.zero.cutoff,
                 what = "zero",
                 path = file.path(path, "4 Zero filter")
               )
               #save data
               met.data.zero.filter <- met.data
               save(met.data.zero.filter,
                    file = file.path(path.inter, "met.data.zero.filter"))
             }

             cat("------------------------------------------------------------------\n")
             cat("Peak identification...\n")
             #peak identification
             if (any(dir(path) == "peak identification")) {
               if (all(dir(path.inter) != "met.data.peak.iden")) {
                 met.data <- PeakIdentification(
                   MetFlowData = met.data,
                   ##parameters for matching
                   mz.tolerance = mz.tolerance,
                   rt.tolerance = rt.tolerance
                 )
                 #save data
                 met.data.peak.iden <- met.data
                 save(met.data.peak.iden,
                      file = file.path(path.inter, "met.data.peak.iden"))
               } else{
                 load(file.path(path.inter, "met.data.peak.iden"))
                 met.data <- met.data.peak.iden
               }
             }

             if (hmdb.matching) {
               #mass identification

               if (all(dir(path.inter) != "met.data.mass.iden")) {
                 met.data <- MassIdentification(
                   MetFlowData = met.data,
                   mass.tolerance = mass.tolerance,
                   polarity = "positive",
                   show = 5
                 )
                 #save data
                 met.data.mass.iden <- met.data
                 save(met.data.mass.iden,
                      file = file.path(path.inter, "met.data.mass.iden"))
               } else{
                 load(file.path(path.inter, "met.data.mass.iden"))
                 met.data <- met.data.mass.iden
               }
             }


               cat("------------------------------------------------------------------\n")
               cat("QC outlier filtering...\n")
               #QC outlier
               if (all(dir(path.inter) != "met.data.qc.outlier.filter")) {
                 met.data <- QCOutlierFilter(MetFlowData = met.data,
                                             CI = 0.95,
                                             qc.outlier.filter = qc.outlier.filter,
                                             path = file.path(path, "5 QC outlier filter"))
                 met.data.qc.outlier.filter <- met.data
                 save(
                   met.data.qc.outlier.filter,
                   file = file.path(path.inter, "met.data.qc.outlier.filter")
                 )
               } else {
                 load(file.path(path.inter, "met.data.qc.outlier.filter"))
                 met.data <- met.data.qc.outlier.filter
               }


             if (normalization) {
               cat("------------------------------------------------------------------\n")
               cat("Data normalization...\n")
               ##Data Normalization
               if (length(batch) > 1) {
                 peak.plot = FALSE
               } else{
                 peak.plot = FALSE
               }
               if (all(dir(path.inter) != "met.data.nor")) {
                 met.data <- DataNormalization(
                   MetFlowData = met.data,
                   method = method,
                   multiple = multiple,
                   threads = threads,
                   path = path,
                   peakplot = peak.plot
                 )
                 #save data
                 met.data.nor <- met.data
                 save(met.data.nor, file = file.path(path.inter, "met.data.nor"))
               } else{
                 load(file.path(path.inter, "met.data.nor"))
                 met.data <- met.data.nor
               }
             }

             # if (subject.outlier.filter) {
               cat("------------------------------------------------------------------\n")
               cat("Subject outlier filtering...\n")
               #subject outlier
               if (all(dir(path.inter) != "met.data.subject.outlier.filter")) {
                 met.data <- SubjectOutlierFilter(
                   MetFlowData = met.data,
                   CI = 0.95,
                   subject.outlier.filter = subject.outlier.filter,
                   path = file.path(path, "6 Subject outlier finder")
                 )
                 met.data.subject.outlier.filter <- met.data
                 save(
                   met.data.subject.outlier.filter,
                   file = file.path(path.inter, "met.data.subject.outlier.filter")
                 )
               } else {
                 load(file.path(path.inter, "met.data.subject.outlier.filter"))
                 met.data <- met.data.subject.outlier.filter
               }
             # }

             if (integration) {
               cat("------------------------------------------------------------------\n")
               cat("Data integration...\n")

               if (length(batch) > 1) {
                 ## data integration
                 met.data <- DataIntegration(MetFlowData = met.data)
                 #save data
                 met.data.integration <- met.data
                 save(met.data.integration,
                      file = file.path(path.inter, "met.data.integration"))
               }
             }

             if (zero.filter) {
               if (length(batch) > 1) {
                 #batch effect
                 BatchEffectOverview(
                   MetFlowData.before = met.data.zero.filter,
                   MetFlowData.after = met.data,
                   path = file.path(path, "8 Batch effect")
                 )
               }

               if (met.plot) {
                 #
                 cat(
                   "------------------------------------------------------------------\n"
                 )
                 cat("Metabolite plot...\n")
                 #metabolite plot
                 MetabolitePlot(
                   MetFlowData.before = met.data.zero.filter,
                   MetFlowData.after = met.data,
                   path = file.path(path, "9 metabolite plot")
                 )
               }
             }

               #ouput data
               ExportData(
                 MetFlowData = met.data,
                 data.name = "peak.table.after.data.cleaning",
                 subject.info.name = "subject.info",
                 qc.info.name = "qc.info", path = path
               )

             cat("\n")
             if (zero.filter) {
               cat("------------------------------------------------------------------\n")
               cat("RSD overview...\n")
               #RSD distribution
               RSDoverview(
                 MetFlowData.before = met.data.zero.filter,
                 MetFlowData.after = met.data,
                 path = file.path(path, "10 RSD overview")
               )
             }

             #data overview
             cat("------------------------------------------------------------------\n")
             cat("Data overview...\n")
             DataOverview(
               MetFlowData = met.data,
               feature.distribution = TRUE,
               path = file.path(path, "11 Data overview")
             )

             #save data
             met.data.after.pre <- met.data
             save(met.data.after.pre, file = file.path(path, "met.data.after.pre"))
             cat("------------------------------------------------------------------\n")
             cat("metflowR is done.\n")
             options(warn = 0)

           })

##RSDoverview function
RSDoverview <- function(MetFlowData.before = MetFlowData1,
                        MetFlowData.after = MetFlowData2,
                        path = NULL) {
  if (is.null(path)) {
    path <- getwd()
  } else{
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
    colours1[QC.rsd > 30] <- "firebrick1"
    colours1[is.na(colours1)] <- "black"

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
           paste("RSD<30%: ",
                 round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
           bty = "n",
           cex = 1.3)

    #after
    QC <- qc.aft1[[i]]
    QC.rsd <- apply(QC, 1, function(x) {
      sd(x) * 100 / mean(x)
    })
    par(mar = c(5, 5, 4, 2))
    colours1 <- rep(NA, length(QC.rsd))
    colours1[QC.rsd > 30] <- "firebrick1"
    colours1[is.na(colours1)] <- "black"

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
           paste("RSD<30%: ",
                 round(sum(QC.rsd < 30) / length(QC.rsd), 4) * 100, "%"),
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
  colours1[qc.rsd > 30] <- "firebrick1"
  colours1[is.na(colours1)] <- "black"
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
         paste("RSD<30%: ",
               round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  # after
  qc.rsd <- apply(qc.aft, 1, function(x) {
    sd(x) * 100 / mean(x)
  })
  colours1 <- rep(NA, length(qc.rsd))
  colours1[qc.rsd > 30] <- "firebrick1"
  colours1[is.na(colours1)] <- "black"
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
         paste("RSD<30%: ",
               round(sum(qc.rsd < 30) / length(qc.rsd), 4) * 100, "%"),
         bty = "n",
         cex = 1.3)

  dev.off()
  layout(1)
}



.onAttach <- function(libname, pkgname){
packageStartupMessage("metflowR
Authors: Xiaotao Shen and Dr. Zhengjiang Zhu
Maintainer: Xiaotao Shen.\n2019-03-23
Version 0.99.00
--------------
o First release version.")
}

packageStartupMessage("metflowR
Authors: Xiaotao Shen and Dr. Zhengjiang Zhu
Maintainer: Xiaotao Shen.\n2019-03-23
Version 0.99.00
--------------
o First release version.")



