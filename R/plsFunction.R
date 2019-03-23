#' @title PLSanalysis
#' @description PLS analysis for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param log.scale Which transformation methd you want to use? 2, 10  or "e",
#' defaulit is FALSE mean don't transformation.
#' @param scalemethod Whihch scale methd you want to use? "auto" or "pareto",
#' defaulit is "auto".
#' @param plsmethod Default is "plsreg".
#' @param path Work directory.
#' @param width width of plot.
#' @param height height of plot.
#' @param text Add text in PCA score plot or not? Deefault is FALSE.
#' @param ellipse Add ellipse in PCA score plot or not? Deefault is TRUE.
#' @param color Color for different class.
#' @param shape Shape for different class.
#' @param cexa cex.
#' @param xlim1 x axis limitation. Default is NULL.
#' @param ylim1 y axis limitation. Default is NULL.
#' @return PLS score plot: PLS score plot.
#' @return permutation test plot: Permutation test plot.
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
#'met.data <- FoldChange(MetFlowData = met.data.after.pre)
#'##PLSanalysis
#'PLSanalysis(met.data,
#'            plsmethod = "plsr",
#'            path = "Demo for PLSanalysis")
#' }

PLSanalysis <- function(MetFlowData,
                        #used data
                        log.scale = 10,
                        scalemethod = "auto",
                        plsmethod = "plsreg",
                        path = ".",
                        width = 7,
                        height = 7,
                        text = FALSE,
                        ellipse = TRUE,
                        color = c("palegreen",
                                  "royalblue",
                                  "firebrick1",
                                  "yellow",
                                  "black",
                                  "cyan",
                                  "gray48"),
                        shape = c(17, 19, 15, 18, 2, 8, 11),
                        cexa = 1,
                        xlim1 = NULL,
                        ylim1 = NULL)
#parameter setting
{
  #
  # requireNamespace("pls")
  if (path != ".") {
    dir.create(path)
  }
  options(warn = -1)

  subject <- MetFlowData@subject
  qc <- MetFlowData@qc
  subject.info <- MetFlowData@subject.info
  group <- subject.info[, "group"]
  group.unique <- sort(unique(group))
  subject.name <- subject.info[, 1]

  if (is.null(qc)) {
    QC <- FALSE
  }

  info <- list()
  for (i in seq_along(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }

  names(info) <- group.unique

  int <- t(subject)
  index <- NULL
  for (i in seq_along(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index, index1)
  }
  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }


  index <- index[!is.na(index)]
  index <- match(index, rownames(int))
  index <- index[!is.na(index)]
  int <- int[index, ]

  #######
  name <- rownames(int)
  #
  Y <- NULL
  label <- list()
  for (i in seq_along(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
    Y[label[[i]]] <- i - 1
  }

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

  int.scale <- SXTscale(int, method = scalemethod)
  # int.Y<-SXTscale(Y,method=scalemethod)
  int.Y <- Y
  #
  ncompa <- min(nrow(int), ncol(int)) - 1

  if (plsmethod == "plsr") {
    pls1 <-
      pls::plsr(
        int.Y ~ int.scale,
        scale = FALSE,
        validation = "CV",
        ncomp = ncompa,
        method = "oscorespls"
      )
    save(pls1, file = "pls1")

    #########select the number of compents#################
    msep <- pls::MSEP(pls1)
    save(msep, file = file.path(path, "msep"))
    msep <- msep$val[, ,]

    yesnot <- "y"
    while (yesnot == "y" | yesnot == "") {
      comps.number <- readline("How many comps do you want to see?")
      while (!exists("comps.number") | comps.number == "")
      {
        cat("You must give a comps number to continute.\n")
        comps.number <-
          readline("How many comps do you want to see?")
      }
      comps.number <- as.numeric(comps.number)
      par(mar = c(5,5,4,2))
      plot(
        x = c(1:comps.number),
        y = msep[1, 2:(comps.number + 1)],
        type = "b",
        col = "tomato",
        pch = 19,
        xlab = "ncomp",
        ylab = "MSEP",
        cex.lab = 1.5,
        cex.axis = 1.3,
        ylim= c(0.8 * min(msep[c(1, 2), 2:(comps.number + 1)]),
                1.2 * max(msep[c(1, 2), 2:(comps.number + 1)]))
      )
      points(
        x = c(1:comps.number),
        y = msep[2, 2:(comps.number + 1)],
        type = "b",
        pch = 17,
        col = "grey"
      )
      legend(
        "top",
        legend = c("CV", "adjCV"),
        col = c("tomato", "grey"),
        pch = c(19, 19),
        bty = "n",
        cex = 1.3,
        pt.cex = 1.3
      )
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    pdf(file.path(path, "MSEP plot.pdf"))
    par(mar = c(5,5,4,2))
    plot(
      x = c(1:comps.number),
      y = msep[1, 2:(comps.number + 1)],
      type = "b",
      col = "tomato",
      pch = 19,
      xlab = "ncomp",
      ylab = "MSEP",
      cex.lab = 1.5,
      cex.axis = 1.3,
      ylim= c(0.99 * min(msep[c(1, 2), 2:(comps.number + 1)]),
              1.1 * max(msep[c(1, 2), 2:(comps.number + 1)]))
    )
    points(
      x = c(1:comps.number),
      y = msep[2, 2:(comps.number + 1)],
      type = "b",
      pch = 17,
      col = "grey"
    )
    legend(
      "top",
      legend = c("CV", "adjCV"),
      col = c("tomato", "grey"),
      pch = c(19, 19),
      bty = "n",
      cex = 1.3,
      pt.cex = 1.3
    )
    dev.off()

    number <-
      readline("Please type number and press Enter to continute:")
    while (!exists("number") | number == "") {
      cat("You must give a number to continute.\n")
      number <-
    readline("Please type comps number value and press Enter  to continute:")
    }
    number <- as.numeric(number)

    ##################construct final pls model###################
    pls2 <-
      pls::plsr(
        int.Y ~ int.scale,
        scale = FALSE,
        validation = "CV",
        ncomp = number,
        method = "oscorespls"
      )
    save(pls2, file = file.path(path, "pls2"))
    vip <- SXTvip(pls2)
    save(vip, file = file.path(path, "vip"))

    scores <- scores(pls2)
    x <- scores[, 1]
    y <- scores[, 2]
    if (number > 2) {
      z <- scores[, 3]
      zmin <- 1.2 * min(z)
      zmax <- 1.2 * max(z)
    }

    xmin <- 1.2 * min(x)
    xmax <- 1.2 * max(x)
    ymin <- 1.2 * min(y)
    ymax <- 1.2 * max(y)
  }

  else {
    #
    dummy <- SXTdummy(Y)
    # int.dummy<-SXTscale(dummy,method=scalemethod)
    int.dummy <- dummy
    # ncompa = nrow(int.scale) - 1
    ncompa <- min(nrow(int), ncol(int)) - 1
    pls1 <- plsdepot::plsreg1(int.scale, Y, comps = ncompa)
    save(pls1, file = file.path(path, "pls1"))
    #########select the number of compents#################
    Q2cum <- pls1$Q2[, 5]
    Q2cum[is.nan(Q2cum)] <- 1
    yesnot <- "y"
    while (yesnot == "y" | yesnot == "") {
      comps.number <- readline("How many comps do you want to see?")
      while (!exists("comps.number") |
             comps.number == "") {
        cat("You must give a comps number to continute.\n")
        comps.number <-
          readline("How many comps do you want to see?")
      }
      comps.number <- as.numeric(comps.number)
      par(mar = c(5,5,4,2))
      a <-
        barplot(
          Q2cum[1:comps.number],
          xlab = "ncomp",
          ylab = "Q2cum",
          cex.lab = 1.5,
          cex.axis = 1.3
        )
      abline(h = 0)
      points(
        a,
        Q2cum[1:comps.number],
        type = "b",
        col = "tomato",
        pch = 20,
        cex = 2
      )
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }
    pdf(file.path(path, "Q2cum plot.pdf"),
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    a <- barplot(
      Q2cum[1:comps.number],
      xlab = "ncomp",
      ylab = "Q2cum",
      cex.lab = 1.5,
      cex.axis = 1.3
    )
    abline(h = 0)
    points(
      a,
      Q2cum[1:comps.number],
      type = "b",
      col = "tomato",
      pch = 20,
      cex = 2
    )
    dev.off()

    number <-
      readline("Please type number and press Enter to continute:")
    while (!exists("number") |
           number == "") {
      cat("You must give a number to continute.\n")
      number <-
    readline("Please type comps number value and press Enter to continute:")
    }
    number <- as.numeric(number)

    ##################construct final pls model###################
    cat(paste(
      "Construct PLS model with all peaks using",
      number,
      "comps ...",
      "\n"
    ))
    pls2 <- plsdepot::plsreg1(int.scale, Y, comps = number)
    save(pls2, file = file.path(path, "pls2"))
    pls.temp <- plsdepot::plsreg2(int.scale, int.dummy, comps = number)
    vip <- pls.temp$VIP
    Q2cum <- pls2$Q2[, 5]
    R2cum <- cumsum(pls2$R2)

    write.csv(cbind(R2cum, Q2cum),
              file.path(path, "R2Q2.csv"),
              row.names = FALSE)

    ##draw barplot of Q2cum, R2Xcum and R2Ycum
    Q2R2 <- cbind(R2cum, Q2cum)
    colnames(Q2R2) <- c("R2cum", "Q2cum")
    pdf(file.path(path, "Q2R2cum.pdf"),
        width = 8,
        height = 6)
    par(mar = c(5,5,4,2))
    barplot(
      t(Q2R2),
      beside = TRUE,
      col = c("royalblue", "tomato"),
      cex.lab = 1.5,
      cex.axis = 1.3,
      cex.names = 1.5
    )
    abline(h = 0)
    legend(
      "topleft",
      legend = c("R2Ycum", "Q2cum"),
      pch = 15,
      col = c("royalblue", "tomato"),
      cex = 1.3,
      pt.cex = 1.3,
      bty = "n"
    )
    dev.off()

    save(vip, file = file.path(path, "vip"))
    save(Q2cum, file = file.path(path, "Q2cum"))
    save(R2cum, file = file.path(path, "R2cum"))

    x <- pls2$x.scores[, 1]
    y <- pls2$x.scores[, 2]
    if (number > 2) {
      z <- pls2$x.scores[, 3]
      zmin <- 1.2 * min(z)
      zmax <- 1.2 * max(z)
    }

    xmin <- 1.2 * min(x)
    xmax <- 1.2 * max(x)
    ymin <- 1.2 * min(y)
    ymax <- 1.2 * max(y)
  }

  legend <- NULL
  for (i in seq_along(label)) {
    legend[label[[i]]] <- names(info)[i]
  }

  colour <- NULL
  colourlist <- color
  for (i in seq_along(label)) {
    colour[label[[i]]] <- colourlist[i]
  }


  pcha <- NULL
  pchalist <- shape
  for (i in seq_along(label)) {
    pcha[label[[i]]] <- pchalist[i]
  }

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

  #PLS 2D
  pdf(file.path(path, "pls score plot.pdf"),
      width = width,
      height = height)
  par(mar = c(5,5,4,2))
  plot(
    x,
    y,
    xlim = xlim,
    ylim = ylim,
    col = colour,
    pch = pcha,
    xlab = "t[1]",
    ylab = "t[2]",
    cex = cexa,
    cex.axis = 1.3,
    cex.lab = 1.35
  )
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  if (text) {
    text(x, y, rownames(int), pos = 4)
  }
  if (ellipse) {
    lines(ellipse::ellipse(
      0,
      scale = c(sd(x), sd(y)),
      centre = c(mean(x), mean(y))
    ), lty = 2)
  }

  legend(
    "topleft",
    names(info),
    pch = pchalist[seq_along(info)],
    col = colourlist[seq_along(info)],
    bty = "n",
    cex = 1.5
  )
  dev.off()

  #PLS 3D
  # if (number > 2) {
  #   pdf(file.path(path, "pls score plot 3d.pdf"),
  #       width = width,
  #       height = height)
  #   scatterplot3d::scatterplot3d(
  #     x,
  #     y,
  #     z,
  #     color = colour,
  #     xlab = "t[1]",
  #     ylab = "t[2]",
  #     zlab = "t[3]",
  #     angle = 50,
  #     pch = pcha,
  #     box = FALSE,
  #     cex.symbol = cexa,
  #     cex.lab = 1.5,
  #     cex.axis = 1.3,
  #     xlim = xlim,
  #     ylim = ylim,
  #     zlim = c(zmin, zmax)
  #   )
  #   legend(
  #     "topleft",
  #     names(info),
  #     pch = pchalist[1:length(info)],
  #     col = colourlist[1:length(info)],
  #     bty = "n",
  #     cex = 1.5
  #   )
  #   dev.off()
  # }
  #
  PLSpermutation(
    data = t(subject),
    log.scale = log.scale,
    path = path,
    info = info,
    repeats = 200,
    ncomp = number,
    scalemethod = scalemethod
  )
  cat("\n")
  cat("PLS analysis is done\n")
}








PLSpermutation <- function(data = NULL,
                           log.scale = FALSE,
                           info = NULL,
                           repeats = 200,
                           ncomp = 3,
                           scalemethod = "auto",
                           path = NULL) {
  options(warn = -1)
  #
  if (is.null(path)) {
    path <- getwd()
  }
  else{
    dir.create(path)
  }

  if (is.null(data))
    stop("sample is NULL")
  if (is.null(info))
    stop("info must not be NULL")

  int <- data
  index <- unlist(info)

  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }

  index <- index[!is.na(index)]
  index <- match(index, rownames(int))
  index <- index[!is.na(index)]
  int <- int[index,]

  #######
  name <- rownames(int)
  #
  Y <- NULL
  label <- as.list(rep(NA, length(info)))
  for (i in seq_along(info)) {
    label[[i]] <- match(as.character(info[[i]]), name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
    Y[label[[i]]] <- i - 1
  }


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

  int.scale <- SXTscale(int, method = scalemethod)
  # int.Y<-SXTscale(Y,method=scalemethod)

  # int.dummy<-SXTscale(dummy,method=scalemethod)
  pls <- plsdepot::plsreg1(int.scale, Y, comps = ncomp)
  Q2 <- pls$Q2[ncomp, 5]
  R2 <- sum(pls$R2)

  save(pls, file = file.path(path, "pls"))
  save(Q2, file = file.path(path, "Q2"))
  save(R2, file = file.path(path, "R2"))

  ##begin repeat
  q2 <- rep(NA, repeats)
  r2 <- rep(NA, repeats)
  cor <- rep(NA, repeats)
  cat("Permutation test...\n")
  for (i in 1:repeats) {
    temp.Y <- Y[order(sample(1:length(Y), length(Y)))]
    temp.pls <- plsdepot::plsreg1(int.scale, temp.Y, comps = ncomp)
    q2[i] <- temp.pls$Q2[ncomp, 5]
    r2[i] <- sum(temp.pls$R2)
    cor[i] <- abs(cor(Y, temp.Y))
    cat(i)
    cat(" ")
  }

  save(q2, file = file.path(path, "q2_200"))
  save(r2, file = file.path(path, "r2_200"))
  save(cor, file = file.path(path, "cor"))

  ##draw perumtation test
  pdf(file.path(path, "Permutation test.pdf"))
  par(xpd = FALSE)
  par(mar = c(5, 5, 4, 2))
  plot(
    x = 0,
    y = 0,
    xlim = c(0, 1),
    ylim = c(min(c(q2, r2)), 1),
    col = "white",
    xlab = "Correlation",
    ylab = "Values",
    cex.axis = 1.3,
    cex.lab = 1.5
  )
  abline(h = 0, lty = 2)

  points(x = cor,
         y = q2,
         col = "tomato",
         pch = 19)
  points(x = cor,
         y = r2,
         col = "royalblue",
         pch = 19)

  points(x = 1,
         y = Q2,
         col = "tomato",
         pch = 19)
  points(x = 1,
         y = R2,
         col = "royalblue",
         pch = 19)

  lm.r2 <- lm(c(R2, r2) ~ c(1, cor))
  lm.q2 <- lm(c(Q2, q2) ~ c(1, cor))

  intercept.q2 <- lm.q2$coefficients[1]
  intercept.r2 <- lm.r2$coefficients[1]

  segments(
    x0 = 0,
    y0 = intercept.q2,
    x1 = 1,
    y1 = Q2,
    lty = 2,
    lwd = 2
  )
  segments(
    x0 = 0,
    y0 = intercept.r2,
    x1 = 1,
    y1 = R2,
    lty = 2,
    lwd = 2
  )


  legend(
    "bottomright",
    title = "Intercepts",
    legend = c(paste("Q2", round(intercept.q2, 2), sep = ": "),
               paste("R2", round(intercept.r2, 2), sep = ": ")),
    col = c("tomato", "royalblue"),
    pch = 19,
    pt.cex = 1.3,
    cex = 1.5,
    bty = "n"
  )
  par(xpd = TRUE)
  dev.off()
  options(warn = 1)
}

