#' Select markers using fold change and p values.
#'
#' @title MarkerSelection
#' @description Select markers using fold change and p values.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path The directory you want to write results.
#' @param foldchange which foldchange information your want to use?
#' @param p which p value information your want to use?
#' @param foldchange.cutoff The cutoff value of fold change.
#' @param p.cutoff The cutoff of p value.
#' @param vip.cutoff The cutoff of VIP value.
#' @param vip VIP.
#' @return MetFlowData: MetFlowData which has been added is.
#' marker information into tags.
#' @return maker.information.csv: The marker information.
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
#'# FoldChange
#'met.data <- FoldChange(MetFlowData = met.data.after.pre,
#'                       to = c("1", "0"))
#'##PLSanalysis
#'PLSanalysis(met.data,
#'            plsmethod = "plsr",
#'            path = "Demo for PLSanalysis")
#'
#'load("Demo for PLSanalysis/vip")
#'vip <- apply(vip, 2, mean)
#'##VIP
#'tags <- met.data[["tags"]]
#'tags <- data.frame(tags, vip)
#'met.data[["tags"]] <- tags
#'##UnivariateTest
#'met.data <- UnivariateTest(met.data)
#'## run
#'new.met.data <- MarkerSelection(met.data)
#'}

MarkerSelection <- function(MetFlowData,
                            foldchange = "foldchange",
                            p = "p.correct",
                            vip = "vip",
                            foldchange.cutoff = c(4 / 3, 3 / 4),
                            p.cutoff = 0.05,
                            vip.cutoff = 1,
                            path = ".") {
  #
  if (path != ".") {
    dir.create(path)
  }

  tags <- MetFlowData@tags
  if (any(colnames(tags) == "is.marker")) {
    warning("Markers have been selected.")
  }

  f.cutoff1 <- as.numeric(foldchange.cutoff[1])
  f.cutoff2 <- as.numeric(foldchange.cutoff[2])
  p.cutoff <- as.numeric(p.cutoff)
  vip.cutoff <- as.numeric(vip.cutoff)

  ##foldchange selection
  if (foldchange != FALSE) {
    if (all(colnames(tags) != "foldchange")) {
      stop("Please calculate fold change first.")
    }
    foldchange <- as.numeric(tags[, foldchange])
    marker.index1 <-
      which(foldchange > f.cutoff1 | foldchange < f.cutoff2)
  }
  else {
    marker.index1 <- c(1:nrow(tags))
  }

  ##vip selection
  if (vip != FALSE) {
    if (all(colnames(tags) != "vip")) {
      stop("Please calculate VIP first.")
    }
    vip <- as.numeric(tags[, "vip"])
    marker.index2 <- which(vip > vip.cutoff)
  } else {
    marker.index2 <- c(1:nrow(tags))
  }

  ##p selection
  if (p != FALSE) {
    if (all(colnames(tags) != "p")) {
      stop("Please calculate p value first.")
    }
    p <- as.numeric(tags[, p])
    marker.index3 <- which(p < p.cutoff)
  } else {
    marker.index3 <- c(1:nrow(tags))
  }

  marker.index <-
    intersect(intersect(marker.index1, marker.index2), marker.index3)
  if (length(marker.index) == 0) {
    stop("No marker are selected, please change canditios and try again.")
  }
  cat(paste(
    "There are",
    length(marker.index),
    "variables are selected as marker.\n"
  ))
  marker.info <- tags[marker.index, ]
  write.csv(marker.info, file.path(path, "marker.info.csv"), row.names = FALSE)
  is.marker <- rep(NA, nrow(tags))
  is.marker[marker.index] <- "yes"
  is.marker[is.na(is.marker)] <- "no"

  if (any(colnames(tags) == "is.marker")) {
    tags[, "is.marker"] <- is.marker
  } else{
    tags <- cbind(tags, is.marker)
  }

  MetFlowData@tags <- as.matrix(tags)
  MetFlowData@marker.selection.condition <-
    paste(
      "foldchange:",
      foldchange,
      "fc:",
      paste(foldchange.cutoff, collapse = ","),
      "p:",
      p,
      "pc:",
      p.cutoff,
      "vip:",
      vip,
      "vc:",
      vip.cutoff
    )
  return(MetFlowData)
}





#' @title MarkerShow
#' @description Draw boxplot for markers.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path The directory you want to write results.
#' @param beeswarm Do you want draw beeswarm on the boxplot? Deafult is TRUE.
#' @return Box plot of markers.
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
#'# FoldChange
#'met.data <- FoldChange(MetFlowData = met.data.after.pre,
#'                       to = c("1", "0"))
#'##PLSanalysis
#'PLSanalysis(met.data,
#'            plsmethod = "plsr",
#'            path = "Demo for PLSanalysis")
#'
#'load("Demo for PLSanalysis/vip")
#'vip <- apply(vip, 2, mean)
#'##VIP
#'tags <- met.data[["tags"]]
#'tags <- data.frame(tags, vip)
#'met.data[["tags"]] <- tags
#'##UnivariateTest
#'met.data <- UnivariateTest(met.data)
#'## run
#'new.met.data <- MarkerSelection(met.data)
#'run
#'MarkerShow(new.met.data, path = "Demo for MarkerShow")
#'}

MarkerShow <- function(MetFlowData,
                       beeswarm = TRUE,
                       path = ".") {
  #
  options(warn = -1)

  if (path != ".") {
    dir.create(path)
  }

  # library(beeswarm)
  tags <- MetFlowData@tags
  subject <- MetFlowData@subject
  subject.info <- MetFlowData@subject.info
  if (all(colnames(tags) != "is.marker")) {
    stop("Please select marker first(use MarkerSelection function).")
  }

  is.marker <- tags[, "is.marker"]
  marker.index <- which(is.marker == "yes")
  feature.name <- tags[, "name"]
  group <- subject.info[, "group"]
  subject.name <- subject.info[, 1]
  group.unique <- sort(unique(group))

  info <- list()
  for (i in seq_along(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }
  names(info) <- group.unique

  for (i in marker.index) {
    temp.data <- list()
    for (j in seq_along(info)) {
      temp.data[[j]] <-
        as.numeric(subject[i,])[match(info[[j]], subject.name)]
    }
    names(temp.data) <- names(info)
    pdf(file.path(path, paste(feature.name[i], "boxplot.pdf")))
    par(mar = c(5, 5, 4, 2))
    boxplot(
      temp.data,
      xlab = "Group",
      ylab = "Intensity",
      cex.lab = 1.5,
      cex.axis = 1.3
    )
    if (beeswarm) {
      beeswarm::beeswarm(temp.data,
                         pch = 19,
                         add = TRUE,
                         col = "grey")
    }
    dev.off()
  }
  options(warn = 0)
}
