#' @title MetStat
#' @description A whole work flow for high throughput MS based metabolomics
#' data statistical analysis.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param new.group New group information or not?
#' @param rsd.cutoff RSD cutoff.
#' @param log.scale Data transformation method. Default is FALSE.
#' @param QC Use qc data for PCA analyis or not? Default is FALSE.
#' @param scale.method Which scale methd you want to use? "auto" or "pareto",
#' defaulit is "auto".
#' @param xlim1 x axis limitation. Default is NULL.
#' @param ylim1 y axis limitation. Default is NULL.
#' @param plsmethod Default is "plsreg".
#' @param fc Default is TRUE.
#' @param to Which the ratio you want? default is case/control.
#' @param ratio Which ratio you want to use to calculate fold
#'  change. median ot mean.
#' @param test.method Which test you want to use? "t" means
#' stutent t test and "wilcox" mean wilcoxon test.
#' @param adjust.method p value correction method. See p.adjust function.
#' @param foldchange which foldchange information your want to use?
#' @param p which p value information your want to use?
#' @param foldchange.cutoff The cutoff value of fold change.
#' @param p.cutoff The cutoff of p value.
#' @param vip.cutoff The cutoff of VIP value.
#' @param vip VIP.
#' @param x x axis factor.
#' @param y y axis factor.
#' @param z z axis factor.
#' @param col Colour for markers and non-markers.
#' @param variable "all" or "marker" for heatmap.
#' @param Group group for heatmap.
#' @param beeswarm Do you want draw beeswarm on the boxplot? Deafult is TRUE.
#' @param pls.analysis Do PLS or not.
#' @param pca.analysis Do PCA or not.
#' @param uni.test Do univariate test or not.
#' @param path The directory you want to write results.
#' @return All the results can be got form other functions and instruction.
#' @export
#' @details The manual of MetCleaning can be found in \href{https://github.com/jaspershen/MetCleaning/blob/master/vignettes/MetCleaning.pdf}{github}
#' or in \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:d8b5fcff-c725-4689-a97d-7ff106322fb6}{my library}.
#' @examples
#' \donttest{
#' #' #load the demo data
#' data(data, package = "MetCleaning")
#' data(sample.information, package = "MetCleaning")
#'
#' ##create a folder for MetCleaning demo
#' dir.create("Demo for MetCleaning")
#' setwd("Demo for MetCleaning")
#'
#' # export the demo data as csv
#' write.csv(data, "data.csv", row.names = FALSE)
#' write.csv(sample.information, "sample.information.csv", row.names = FALSE)
#'
#' #run MetCleaning
#' MetCleaning(#ImportData para
#' data = "data.csv",
#' sample.information = "sample.information.csv",
#' polarity = "positive",
#' #DataNormalization
#' method = "svr",
#' threads = 2)
#'
#' ## load the demo data
#'data(new.group, package = "MetCleaning")
#'load("met.data.after.pre")
#'##create a folder for MetStat demo
#'dir.create("Demo for MetStat")
#'setwd("Demo for MetStat")
#'## export the demo data as csv
#'write.csv(new.group, "new.group.csv", row.names = FALSE)
#'
#'## run MetStat
#'MetStat(MetFlowData = met.data.after.pre,
#'new.group = TRUE,
#'Group = c("0", "1"),
#'to = c("1", "0"))
#' }

MetStat <- function(MetFlowData,
                    new.group = TRUE,
                    rsd.cutoff = 30,
                    #transformation para
                    log.scale = FALSE,
                    #PCA analysis para
                    QC = TRUE,
                    scale.method = "auto",
                    xlim1 = NULL,
                    ylim1 = NULL,
                    #PLS analysis para
                    plsmethod = "plsr",
                    #FoldChange para
                    fc = TRUE,
                    to = c("case", "control"),
                    ratio = "median",
                    #UnivariateTest para
                    test.method = "t",
                    adjust.method = "fdr",
                    #MarkerSelection para
                    foldchange = "foldchange",
                    p = "p",
                    vip = "vip",
                    foldchange.cutoff = c(4 / 3, 3 / 4),
                    p.cutoff = 0.05,
                    vip.cutoff = 0,
                    #VolcanoPlot para
                    x = "foldchange",
                    y = "p",
                    z = "vip",
                    col = c("black", "firebrick1"),
                    #heatmap para
                    variable = "all",
                    Group = c("control", "case"),
                    #MarkerShow para
                    beeswarm = TRUE,
                    pls.analysis = TRUE,
                    pca.analysis = TRUE,
                    uni.test = TRUE,
                    path = ".") {
  options(warn = -1)
  #
  if (path != ".") {
    dir.create(path)
  }

  options(warn = -1)
  path.inter <- file.path(path, "intermediate")
  dir.create(path.inter)

  ##RSD filtering
  met.data <-
    RSDfilter(MetFlowData = MetFlowData, rsd.cutoff = rsd.cutoff)
  #save data
  met.data.rsd.filter <- met.data
  save(met.data.rsd.filter,
       file = file.path(path.inter, "met.data.rsd.filter"))

  ##ReChangeGroup
  if (new.group) {
    cat("-----------------------------------------------------------------\n")
    cat("Change group information\n")
    if (all(dir() != "new.group.csv"))
      stop("No new.group information.")
    met.data <- ReChangeGroup(MetFlowData = MetFlowData)
    #save data
    met.data.new.group <- met.data
    save(met.data.new.group,
         file = file.path(path.inter, "met.data.new.group"))
  }

  if (pca.analysis) {
    cat("------------------------------------------------------------------\n")
    cat("PCA analysis\n")
    ##PCA analysis
    PCAanalysis(
      MetFlowData = met.data,
      QC = QC,
      log.scale = log.scale,
      scale.method = scale.method,
      path = file.path(path, "12PCA analysis"),
      xlim1 = xlim1,
      ylim1 = ylim1
    )
  }

  if (pls.analysis) {
    cat("------------------------------------------------------------------\n")
    cat("PLS analysis\n")
    ##PLS analysis
    PLSanalysis(
      MetFlowData = met.data,
      log.scale = log.scale,
      scalemethod = scale.method,
      path = file.path(path, "13PLS analysis"),
      plsmethod = plsmethod,
      xlim1 = xlim1,
      ylim1 = ylim1
    )

    load(file.path(path, "13PLS analysis", "vip"))
    vip <- apply(vip, 2, mean)
  }
#
  ##VIP
  load(file.path(path, "13PLS analysis", "vip"))
  vip <- apply(vip, 2, mean)
  tags <- met.data@tags
  tags <- data.frame(tags, vip)
  met.data@tags <- as.matrix(tags)

  #save data
  met.data.vip <- met.data
  save(met.data.vip, file = file.path(path.inter, "met.data.vip"))

  cat("------------------------------------------------------------------\n")
  cat("Heat map...\n")
  #
  HeatMap(
    MetFlowData = met.data,
    log.scale = log.scale,
    variable = variable,
    Group = Group,
    scale.method = scale.method,
    path = file.path(path, "14heat map")
  )

  if (fc) {
    ##FoldChange
    cat("------------------------------------------------------------------\n")
    met.data <- FoldChange(MetFlowData = met.data,
                           to = to,
                           ratio = ratio)
    #save data
    met.data.fc <- met.data
    save(met.data.fc, file = file.path(path.inter, "met.data.fc"))
  }

  if (uni.test) {
    ##UnivariateTest
    cat("------------------------------------------------------------------\n")
    cat("Univariate test\n")
    met.data <- UnivariateTest(
      MetFlowData = met.data,
      test.method = test.method,
      adjust.method = adjust.method,
      log.scale = log.scale
    )
    #save data
    met.data.uni.test <- met.data
    save(met.data.uni.test,
         file = file.path(path.inter, "met.data.uni.test"))
  }

  ##MarkerSelection
  cat("------------------------------------------------------------------\n")
  cat("Marker selection\n")
  met.data <- MarkerSelection(
    MetFlowData = met.data,
    foldchange = foldchange,
    p = p,
    vip = vip,
    foldchange.cutoff = foldchange.cutoff,
    p.cutoff = p.cutoff,
    vip.cutoff = vip.cutoff,
    path = file.path(path, "15marker selection")
  )
  #save data
  met.data.after.stat <- met.data
  save(met.data.after.stat,
       file = file.path(path.inter, "met.data.after.stat"))

  ##Export data
  ExportData(
    MetFlowData = met.data,
    data.name = "data.after.stat",
    subject.info.name = "subject.info",
    qc.info.name = "qc.info"
  )

  #save data
  met.data.after.stat <- met.data
  save(met.data.after.stat,
       file = file.path(path.inter, "met.data.after.stat"))

  ##VolcanoPlot
  cat("------------------------------------------------------------------\n")
  cat("Draw volcano plot\n")
  VolcanoPlot(
    MetFlowData = met.data,
    x = x,
    y = y,
    z = z,
    col = c("black", "firebrick1"),
    foldchange.cutoff = foldchange.cutoff,
    p.cutoff = p.cutoff,
    vip.cutoff = vip.cutoff,
    path = file.path(path, "15marker selection")
  )

  ##MarkerShow
  cat("------------------------------------------------------------------\n")
  cat("Marker show\n")
  MarkerShow(
    MetFlowData = met.data,
    beeswarm = beeswarm,
    path = file.path(path, "15marker selection")
  )
  options(warn = 0)

  cat("------------------------------------------------------------------\n")
  cat("MetStat is done.\n")
  options(warn = 0)
}