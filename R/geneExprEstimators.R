#' @title Estimation of the M and A Values
#'
#' @description Estimates the intensity values M and A, considering p spots/genes and n microarray slides,
#' by using the conventional method \insertCite{yang2002comparison}{omicsMA} or the one that takes into account
#' pixel-level variability \insertCite{ribeiro2019omicsMA}{omicsMA}.
#'
#' @param R.mean A p x n matrix of mean values for the red (test) foreground intensities.
#' @param G.mean A p x n matrix of mean values for the green (reference) foreground intensities.
#' @param R.var A p x n matrix of variances of the red foreground intensities.
#' @param G.var A p x n matrix of variances of the green foreground intensities.
#' @param R.bckg A p x n matrix of mean values of the red background intensities.
#' @param G.bckg A p x n matrix of mean values of the green background intensities.
#' @param RG.cov A p x n matrix of covariances between the green and red foreground intensities.
#' @param estimator A character string indicating the estimation method. The options are 'conventional' and 'second-order'.
#' If estimator is set as 'second-order', then an estimator that takes into account pixel-level variability is used.
#' @param bgcorr A character string indicating the background correction method. The options are 'none' and 'normexp'.
#' @param offset A numeric parameter for the normexp background correction.
#' @param putaway.perc A scalar between 0 and 100 indicating percentage of genes expressions.
#' with high variance to be discarded. Default is 0.05.
#' @param array.ids A character vector with identifiers of the arrays.
#' @param gene.ids A character vector with identifiers of the genes.
#'
#' @return A limma's MAlist object with the M and A values
#'
#' @examples
#' \donttest{
#' data(metaplasia)
#' MA <- estimateMAValues(metaplasia$R.mean, metaplasia$G.mean,
#'                        metaplasia$R.var, metaplasia$G.var, metaplasia$RG.cov,
#'                        metaplasia$R.bckg, metaplasia$G.bckg,
#'                        array.ids=metaplasia$array.ids, gene.ids=metaplasia$gene.ids)
#'}
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @importFrom limma backgroundCorrect
#' @importFrom stats quantile
#' @importFrom methods new
#'
#' @export estimateMAValues
estimateMAValues <- function(R.mean, G.mean, R.var, G.var, RG.cov,
      R.bckg, G.bckg, estimator='second-order',
      bgcorr='normexp', offset=50, putaway.perc = 0.05,
      array.ids=NULL, gene.ids=NULL)
{
  if (bgcorr == 'normexp') {
    RG <-  new("RGList")
    RG$R <- R.mean
    RG$G <- G.mean
    RG$Rb <- R.bckg
    RG$Gb <- G.bckg

    # applying the background correction
    RGc <- backgroundCorrect(RG, method="normexp", offset=offset)
    R.mean <- RGc$R
    G.mean <- RGc$G
  }

  # conventional estimations
  M.mean <- log2(R.mean/G.mean)
  A.mean <- (log2(G.mean) + log2(R.mean))/2

  if (estimator == 'second-order') {
    # improving estimations with second-order terms of the Taylor expansion
    M.mean <- M.mean + 1/(2 * log(2)) * (-R.var/(R.mean^2) + G.var/(G.mean^2))
    A.mean <- A.mean - 1/(4*log(2)) * (G.var/(G.mean^2) + R.var/(R.mean^2))
  }

  M.var <- 1/(log(2)^2) * (R.var/(R.mean^2) + G.var/(G.mean^2) - 2*RG.cov/(R.mean*G.mean))
  if (putaway.perc > 0) {
    #A.var <- 1/(4*(log(2)^2)) * (R.var/(R.mean)^2 + G.var/(G.mean)^2 + 2*RG.cov/(R.mean*G.mean))
    for (arrayId in 1:ncol(M.var)) {
      badSpots <- which(M.var[,arrayId] > quantile(M.var[,arrayId], prob=((100-putaway.perc)/100)))
      M.mean[badSpots,arrayId] <- NA
      A.mean[badSpots,arrayId] <- NA
    }
  }

  rownames(M.mean) <- rownames(A.mean) <- gene.ids
  colnames(M.mean) <- colnames(A.mean) <- array.ids

  return(list(M.mean, A.mean))
}

#' @title Optimal LOWESS Within-Slide Normalization of a MAList
#'
#' @description Applies LOWESS for within-slide normalization of the MA values.
#' The parameters are set by a data-driven parameter selection method, proposed by
#' \insertCite{ribeiro2019omicsMA}{omicsMA},
#' which is parsimonious and considers intrinsic characteristics of microarray data,
#' such as heteroskedasticity. Particularly, the best bandwidth is selected according
#' to the HRCp criterion \insertCite{liu2013heteroscedasticity}{omicsMA}.
#'
#' @references{
#'   \insertAllCited
#' }
#'
#' @param MA A MAList object, as in the limma R package, with the non-normalized MA values.
#' @param eva.values A vector with values between 0 and 1 corresponding to the bandwith values to be considered by the LOWESS parameter selection method.
#' @param debug A logical value indication if you want to view logs of the execution.
#' @param save.objs A logical value indicating if you want to save objects with partial results.
#' @param save.plots A logical value indicating if you want to generate the M plots
#' illustrating the bandwidth parameter selection process.
#' @param dir.to.save Path to the folder you want to save the output objects.
#'
#' @return A MAList object, as in the limma R package,
#' with the MA values after applying within-slide normalization by LOWESS with optimal parameter settings.
#'
#' @examples
#' \donttest{
#' data(metaplasia)
#' MA <- estimateMAValues(metaplasia$R.mean, metaplasia$G.mean,
#'                        metaplasia$R.var, metaplasia$G.var, metaplasia$RG.cov,
#'                        metaplasia$R.bckg, metaplasia$G.bckg,
#'                        array.ids=metaplasia$array.ids, gene.ids=metaplasia$gene.ids)
#' normMA <- normalizeWithinArraysByOptimalLowess(MA)
#' }
#'
#' @importFrom Rdpack reprompt
#'
#' @export normalizeWithinArraysByOptimalLowess
normalizeWithinArraysByOptimalLowess <- function(MA, eva.values=seq(0.2,1.0,by=0.05),
                                                 debug=TRUE, save.objs=TRUE, save.plots=TRUE,
                                                 dir.to.save='./')
{
  normM <- MA[[1]]
  normA <- MA[[2]]

  for (arrayId in 1:ncol(MA[[1]])) {
    M.mean <- MA[[1]][,arrayId]
    A.mean <- MA[[2]][,arrayId]
    array.id =colnames(MA[[1]])[arrayId]

    badMASpots <- badSpots(M.mean, A.mean)

    if (debug) {
      cat(paste0("Number of bad spots in slide ", array.id, ": ", length(badMASpots), "\n"))
    }

    if (length(badMASpots) > 0) {
      M.mean <- M.mean[-badMASpots]
      A.mean <- A.mean[-badMASpots]
    }
    out <- getOptimalLowessFvalue(M.mean, A.mean,
                                  eva.values=eva.values,
                                  debug=debug, save.objs=save.objs,
                                  save.plots=save.plots,
                                  dir.to.save=dir.to.save, array.id=array.id)
    optimal_value <- out[[1]]
    optimal_lowess_fit <- out[[4]][[optimal_value]]

    if (length(badMASpots) > 0) {
      normM[-badMASpots,arrayId] <- residuals(optimal_lowess_fit,
        type = "raw")
      normM[badMASpots,arrayId] <- NA
    } else {
      normM[,arrayId] <- residuals(optimal_lowess_fit, type = "raw")
    }
  }

  return(list(normM, normA))
}


#' @title Optimal Selection of the LOWESS Bandwidth Parameter Using the HRCp Criterion
#'
#' @description LOWESS parameter selection method, proposed by \insertCite{ribeiro2019omicsMA}{omicsMA},
#' that is parsimonious and considers intrinsic characteristics of microarray data,
#' such as heteroskedasticity. Particularly, the best bandwidth is selected according to
#' the HRCp criterion \insertCite{liu2013heteroscedasticity}{omicsMA}.
#'
#' @references{
#'    \insertAllCited
#' }
#'
#' @param M.mean A numeric vector with the M values for the p genes.
#' @param A.mean A numeric vector with the A values for the p genes.
#' @param eva.values A numeric vector with bandwith values to be considered by the LOWESS parameter selection method.
#' @param debug A logical value indication if you want to view logs of the execution.
#' @param save.objs A logical value indicating if you want to save objects with partial results.
#' @param save.plots A logical value indicating if you want to generate the M plots
#' illustrating the bandwidth parameter selection process.
#' @param dir.to.save Path to the folder you want to save the output objects.
#' @param array.id A character string identifying the array.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{opt.value}{The optimal bandwith.}
#' \item{criterion.val}{A list with the mean square errors for the estimates obtained for each evaluated bandwith.}
#' \item{criterion.df}{A list with the effective degrees of freedom for the estimates obtained for each evaluated bandwith.}
#' \item{fits}{A list with the estimates obtained for each evaluated bandwith.}
#' }
#'
#' @importFrom Rdpack reprompt
#' @importFrom locfit locfit
#' @importFrom grDevices dev.off
#' @importFrom stats fitted predict residuals update median
#'
#' @export getOptimalLowessFvalue
getOptimalLowessFvalue <- function(M.mean, A.mean,
                                   eva.values=seq(0.2,1.0,by=0.05),
                                   debug=TRUE, save.objs=TRUE, save.plots=TRUE,
                                   dir.to.save='./', array.id ='array')
{
  locfit.degree=1
  locfit.iter=3

  criterion.df <- rep(NA, length(eva.values))
  names(criterion.df) <- as.character(eva.values)
  criterion.val <- criterion.df
  fits <- list()

  hrcp.name <- 'HRCp'

  # variance fit
  if (debug)
    cat(paste(Sys.time(), " -- ", array.id, " - Processing variance fit... ", "\n", sep=""))

  fit.unbiased <- locfit(M.mean ~ A.mean, deg=locfit.degree,
                         alpha=0.1, maxk=500, kern="tricube")
  y <- -2 * predict(fit.unbiased, what="lik", where="data")
  w <- predict(fit.unbiased, what="rdf", where="data")
  rm(fit.unbiased)

  f_vars <- seq(1, 0.05, by=-0.05)
  max.values <- rep(NA, length(f_vars))
  names(max.values) <- f_vars
  for (f_var in f_vars) {
    max.values[as.character(f_var)] <- tryCatch({
      fit.variance <- locfit(y ~ A.mean, weights=w,
                             family="gamma", alpha=f_var)
      max.fitted.var <- max(fitted(fit.variance))
      max.fitted.var
    },
    error = function(e) {
      NA
    })
  }
  f_var <- names(which(max.values == min(max.values, na.rm=T)))

  # variance fit
  if (debug) {
    cat(paste(Sys.time(), " -- ", array.id, " - f_var = ", f_var, " ; max.value = ", max.values[as.character(f_var)], "\n", sep=""))
    print(max.values)
    cat("\n")
  }

  fit.variance <- locfit(y ~ A.mean, weights=w,
                         family="gamma", alpha=f_var)
  fitted.var <- fitted(fit.variance)

  if (save.objs) {
    save(fit.variance, file=paste(dir.to.save, "fit.variance_" ,
                                  array.id, ".RData", sep=""))
    save(fitted.var, file=paste(dir.to.save, "fitted.var_" ,
                                array.id, ".RData", sep=""))
  }
  rm(fit.variance)

  criterionRObj <- paste(dir.to.save,
                         hrcp.name, "_", array.id, ".RData", sep="")
  if (!file.exists(criterionRObj)) {
    for(f.value in eva.values) {
      hrcp.name = hrcp.name
      ret <- HRCPCriterion(M.mean, A.mean, fitted.var,
                           locfit.degree, f.value, locfit.iter, debug=TRUE)

      criterion.val[as.character(f.value)] <- ret[[1]]
      criterion.df[as.character(f.value)] <- ret[[2]]
      fits[[as.character(f.value)]] <- ret[[3]]

      if (debug)
        cat(paste(Sys.time(), " -- ", array.id, " - ",
                  hrcp.name, " - f.value=", f.value, ": val=",
                  criterion.val[as.character(f.value)],
                  " df=", criterion.df[as.character(f.value)], "\n", sep=""))
    }
  } else {
    load(criterionRObj)
    criterion.df.back <- criterion.df
    criterion.df <- criterion.val
    criterion.val <- criterion.df.back
  }

  rm(fitted.var)

  if (save.objs)
    save(criterion.df, criterion.val, file=criterionRObj)

  # selecting optimized bandwidth value
  range = max(criterion.val) - min(criterion.val)
  c = range*0.05

  opt_idx <- 1
  pnt1 = pnt2 = pnt3 = 0
  for(idx in length(eva.values):1) {
    pnt1 = criterion.val[idx]

    if (abs(pnt3-min(criterion.val)) < c &&
        abs(pnt3 - pnt2) < c  &&
        abs(pnt3 - pnt1) < c) {
      opt_idx <- idx
      break;
    }
    pnt3 = pnt2
    pnt2 = pnt1
  }

  opt.value <- eva.values[opt_idx]

  if (debug)
    print(paste("optimal value: ", opt.value))

  if (save.plots) {
    if (debug)
      cat(paste(Sys.time(), " -- ", "Plotting: ", array.id,
                "_degr_", locfit.degree, "_criterion_",
                hrcp.name, ".png", "\n", sep=""))

    my.plot <- lattice::xyplot(criterion.val ~ criterion.df,
                               main  =  paste("Selection of LOWESS f paramater for slide: ", array.id, sep=""),
                               xlab  = "Degrees of Freedom",
                               ylab  = hrcp.name,
                               xlim  = range(1 : (max(criterion.df)+1)),
                               pch = ifelse(criterion.df == criterion.df[opt_idx], 19, 19),
                               col = ifelse(criterion.df == criterion.df[opt_idx], "black",
                                            ifelse(criterion.val <= (min(criterion.val)+c),
                                                   "turquoise4", "red4")),
                               cex = ifelse(criterion.df == criterion.df[opt_idx], 2.5, 1.5)
    )

    my.plot <- update(my.plot, panel = function(...) {
      lattice::panel.text(criterion.df[opt_idx]+1,
                          criterion.val[opt_idx]+1.5*c,
                          labels=c(as.expression(bquote(alpha == .(opt.value)))),
                          col="black", cex=2.0)
      lattice::panel.grid(h = 19, #subscripts=subscripts,
                          v = ceiling(max(criterion.df))-2)
      lattice::panel.xyplot(...) })


    lattice::trellis.device(device="png", filename=paste(dir.to.save,
                                                         "criterion_", hrcp.name, "_degr_", locfit.degree, "_",
                                                         array.id, ".png", sep=""))
    print(my.plot)
    dev.off()
  }

  return(list(opt.value, criterion.val, criterion.df, fits))
}


badSpots <- function(M, A) {
  badMSpots <- c(which(M==-Inf),which(M==Inf),which(is.na(M)))
  badASpots <- c(which(A==-Inf),which(A==Inf),which(is.na(A)))
  badMASpots <- unique(badMSpots, badASpots)
  return(badMASpots)
}


HRCPCriterion <- function(M.mean, A.mean, fitted.var, locfit.degree,
                          locfit.alpha, locfit.iter, debug=TRUE) {
  if (debug)
    print(paste(Sys.time(), " -- ", "Processing hatmatrix... ", "\n", sep=""))

  lfr.wt <- rep(1, length(M.mean))
  for(i in 0:(locfit.iter-1)) {
    fit <- locfit(M.mean ~ A.mean, weights=lfr.wt, deg=locfit.degree,
                  alpha=0.1, maxk=500,  kern="tricube")
    res <- residuals(fit, type = "raw")
    s <- median(abs(res))
    lfr.wt <- pmax(1 - (res/(6 * s))^2, 0)^2
  }
  rm(fit)

  fit_robust_alpha <- locfit(M.mean ~ A.mean, weights=lfr.wt,
                             deg=locfit.degree, alpha=locfit.alpha, maxk=500,  kern="tricube")

  hatmatrix <- locfit(M.mean ~ A.mean, weights=lfr.wt,
                      deg=locfit.degree, ev=locfit::dat(), alpha=locfit.alpha,
                      maxk=500,  kern="tricube", geth=1)

  traceomegahm <- 0
  for (obs in 1: length(M.mean)) {
    traceomegahm <- traceomegahm + fitted.var[obs] * hatmatrix[obs,obs]
  }
  rm(hatmatrix)

  if (debug)
    cat(paste(Sys.time(), " -- ", "Processing HRCp... ", "\n", sep=""))

  HRCP.val <- (sum(residuals(fit_robust_alpha)^2)) + 2 * traceomegahm
  HRCP.df <- fit_robust_alpha$dp["df2"]
  rm(traceomegahm)

  return(list(HRCP.val, HRCP.df, fit_robust_alpha))
}


