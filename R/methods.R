
.comboGroups <- function(truths)
  # Function that returns a list of vectors of indices,
  # where each vector refers to the rows with the same
  # combination of TRUE/FALSE values in 'truths'.
  #
  # written by Aaron Lun
  # Created 24 October 2014
{
  #	Integer packing will only work for 31 libraries at a time.
  assembly <- list()
  collected <- 0L
  step <- 31L
  bits <- as.integer(2^(1:step-1L))

  while (collected < ncol(truths)) {
    upper <- pmin(ncol(truths) - collected, step)
    keys <- t(truths[,collected+1:upper,drop=FALSE]) * bits[1:upper]
    assembly[[length(assembly)+1L]] <- as.integer(colSums(keys))
    collected <- collected + step
  }

  #	Figuring out the unique components.
  o <- do.call(order, assembly)
  nr <- nrow(truths)
  is.different <- logical(nr)
  for (i in 1:length(assembly)) {
    is.different <- is.different | c(TRUE, diff(assembly[[i]][o])!=0L)
  }
  first.of.each <- which(is.different)
  last.of.each <- c(first.of.each[-1]-1L, nr)

  #	Returning the groups.
  output <- list()
  for (u in 1:length(first.of.each)) {
    output[[u]] <- o[first.of.each[u]:last.of.each[u]]
  }
  return(output)
}



.residDF <- function(zero, design)
  #	Effective residual degrees of freedom after adjusting for exact zeros
  #	Gordon Smyth and Aaron Lun
  #	Created 6 Jan 2014.  Last modified 2 Sep 2014
{
  nlib <- ncol(zero)
  ncoef <- ncol(design)
  nzero <- as.integer(rowSums(zero))

  #	Default is no zero
  DF <- rep(nlib-ncoef,length(nzero))

  #	All zero case
  DF[nzero==nlib] <- 0L

  #	Anything in between?
  somezero <- nzero>0L & nzero<nlib
  if(any(somezero)) {
    zero2 <- zero[somezero,,drop=FALSE]
    groupings <- .comboGroups(zero2)

    #		Identifying the true residual d.f. for each of these rows.
    DF2 <- nlib-nzero[somezero]
    for (u in 1:length(groupings)) {
      i <- groupings[[u]]
      zeroi <- zero2[i[1],]
      DF2[i] <- DF2[i]-qr(design[!zeroi,,drop=FALSE])$rank
    }
    DF2 <- pmax(DF2, 0L)
    DF[somezero] <- DF2
  }

  DF
}



.estimateDispWeighted = function (y, design = NULL, prior.df = NULL, trend.method = "locfit", tagwise = TRUE, span = NULL, min.row.sum = 5, grid.length = 21,
                                 grid.range = c(-10, 10), robust = FALSE, winsor.tail.p = c(0.05,
                                                                                            0.1), tol = 1e-06, weights=NULL)
{
  #adjusted by Koen VdB on 04 March 2016
  if(!is.null(y$weights)) weights=y$weights
  if (!is(y, "DGEList"))
    stop("y must be a DGEList")
  trend.method <- match.arg(trend.method, c("none", "loess",
                                            "locfit", "movingave"))
  ntags <- nrow(y$counts)
  nlibs <- ncol(y$counts)
  offset <- getOffset(y)
  AveLogCPM <- aveLogCPM(y)
  offset <- expandAsMatrix(offset, dim(y))
  sel <- rowSums(y$counts) >= min.row.sum
  spline.pts <- seq(from = grid.range[1], to = grid.range[2],
                    length = grid.length)
  spline.disp <- 0.1 * 2^spline.pts
  grid.vals <- spline.disp/(1 + spline.disp)
  l0 <- matrix(0, sum(sel), grid.length)
  if (is.null(design)) {
    cat("Design matrix not provided. Switch to the classic mode.\n")
    group <- y$samples$group <- as.factor(y$samples$group)
    if (length(levels(group)) == 1)
      design <- matrix(1, nlibs, 1)
    else design <- model.matrix(~group)
    if (all(tabulate(group) <= 1)) {
      warning("There is no replication, setting dispersion to NA.")
      y$common.dispersion <- NA
      return(y)
    }
    pseudo.obj <- y[sel, ]
    q2q.out <- equalizeLibSizes(y[sel, ], dispersion = 0.01)
    pseudo.obj$counts <- q2q.out$pseudo
    ysplit <- splitIntoGroups(pseudo.obj)
    delta <- optimize(commonCondLogLikDerDelta, interval = c(1e-04,
                                                             100/(100 + 1)), tol = tol, maximum = TRUE, y = ysplit,
                      der = 0)
    delta <- delta$maximum
    disp <- delta/(1 - delta)
    q2q.out <- equalizeLibSizes(y[sel, ], dispersion = disp)
    pseudo.obj$counts <- q2q.out$pseudo
    ysplit <- splitIntoGroups(pseudo.obj)
    for (j in 1:grid.length) for (i in 1:length(ysplit)) l0[,
                                                            j] <- condLogLikDerDelta(ysplit[[i]], grid.vals[j],
                                                                                     der = 0) + l0[, j]
  }
  else {
    design <- as.matrix(design)
    if (ncol(design) >= ncol(y$counts)) {
      warning("No residual df: setting dispersion to NA")
      y$common.dispersion <- NA
      return(y)
    }
    glmfit <- glmFit(y$counts[sel, ], design, offset = offset[sel,
                                                              ], dispersion = 0.05, prior.count = 0, weights=weights[sel,]) ###www
    zerofit <- (glmfit$fitted.values < 1e-04) & (glmfit$counts <
                                                   1e-04)
    by.group <- .comboGroups(zerofit)
    for (subg in by.group) {
      cur.nzero <- !zerofit[subg[1], ]
      if (!any(cur.nzero)) {
        next
      }
      if (all(cur.nzero)) {
        redesign <- design
      }
      else {
        redesign <- design[cur.nzero, , drop = FALSE]
        QR <- qr(redesign)
        redesign <- redesign[, QR$pivot[1:QR$rank], drop = FALSE]
        if (nrow(redesign) == ncol(redesign)) {
          next
        }
      }
      last.beta <- NULL
      for (i in 1:grid.length) {
        out <- adjustedProfileLik(spline.disp[i], y = y$counts[sel,
                                                               ][subg, cur.nzero, drop = FALSE], design = redesign,
                                  offset = offset[sel, ][subg, cur.nzero, drop = FALSE],
                                  start = last.beta, get.coef = TRUE, weights=weights[sel,][subg, cur.nzero, drop = FALSE]) ###www
        l0[subg, i] <- out$apl
        last.beta <- out$beta
      }
    }
  }
  out.1 <- WLEB(theta = spline.pts, loglik = l0, covariate = AveLogCPM[sel],
                trend.method = trend.method, span = span, individual = FALSE,
                m0.out = TRUE)
  y$common.dispersion <- 0.1 * 2^out.1$overall
  disp.trend <- 0.1 * 2^out.1$trend
  y$trended.dispersion <- rep(disp.trend[which.min(AveLogCPM[sel])],
                              ntags)
  y$trended.dispersion[sel] <- disp.trend
  y$trend.method <- trend.method
  y$AveLogCPM <- AveLogCPM
  y$span <- out.1$span
  if (!tagwise)
    return(y)
  if (is.null(prior.df)) {
    glmfit <- glmFit(y$counts[sel, ], design, offset = offset[sel,
                                                              ], dispersion = disp.trend, prior.count = 0, weights=weights[sel,]) ###www
    df.residual <- glmfit$df.residual
    zerofit <- (glmfit$fitted.values < 1e-04) & (glmfit$counts <
                                                   1e-04)
    df.residual <- .residDF(zerofit, design)
    s2 <- glmfit$deviance/df.residual
    s2[df.residual == 0] <- 0
    s2 <- pmax(s2, 0)
    s2.fit <- squeezeVar(s2, df = df.residual, covariate = AveLogCPM[sel],
                         robust = robust, winsor.tail.p = winsor.tail.p)
    prior.df <- s2.fit$df.prior
  }
  ncoefs <- ncol(design)
  prior.n <- prior.df/(nlibs - ncoefs)
  if (trend.method != "none") {
    y$tagwise.dispersion <- y$trended.dispersion
  }
  else {
    y$tagwise.dispersion <- rep(y$common.dispersion, ntags)
  }
  too.large <- prior.n > 1e+06
  if (!all(too.large)) {
    temp.n <- prior.n
    if (any(too.large)) {
      temp.n[too.large] <- 1e+06
    }
    out.2 <- WLEB(theta = spline.pts, loglik = l0, prior.n = temp.n,
                  covariate = AveLogCPM[sel], trend.method = trend.method,
                  span = span, overall = FALSE, trend = FALSE, m0 = out.1$shared.loglik)
    if (!robust) {
      y$tagwise.dispersion[sel] <- 0.1 * 2^out.2$individual
    }
    else {
      y$tagwise.dispersion[sel][!too.large] <- 0.1 * 2^out.2$individual[!too.large]
    }
  }
  if (!robust) {
    y$prior.df <- prior.df
    y$prior.n <- prior.n
  }
  else {
    y$prior.df <- y$prior.n <- rep(Inf, ntags)
    y$prior.df[sel] <- prior.df
    y$prior.n[sel] <- prior.n
  }
  y
}



.pvalueAdjustment_kvdb <- function(baseMean, filter, pValue,
                                   theta, alpha=0.05, pAdjustMethod="BH") {
  ## this function has been adapted from the pValueAdjustment function from the DESeq2 Bioconductor package and was originally written by Mike Love.
  require(genefilter)
  # perform independent filtering
  if (missing(filter)) {
    filter <- baseMean
  }
  if (missing(theta)) {
    lowerQuantile <- mean(filter == 0)
    if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
    theta <- seq(lowerQuantile, upperQuantile, length=50)
  }

  # do filtering using genefilter
  stopifnot(length(theta) > 1)
  filtPadj <- filtered_p(filter=filter, test=pValue,
                         theta=theta, method=pAdjustMethod)
  numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
  # prevent over-aggressive filtering when all genes are null,
  # by requiring the max number of rejections is above a fitted curve.
  # If the max number of rejection is not greater than 10, then don't
  # perform independent filtering at all.
  lo.fit <- lowess(numRej ~ theta, f=1/5)
  if (max(numRej) <= 10) {
    j <- 1
  } else {
    residual <- if (all(numRej==0)) {
      0
    } else {
      numRej[numRej > 0] - lo.fit$y[numRej > 0]
    }
    thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
    j <- if (any(numRej > thresh)) {
      which(numRej > thresh)[1]
    } else {
      1
    }
  }
  padj <- filtPadj[, j, drop=TRUE]
  cutoffs <- quantile(filter, theta)
  filterThreshold <- cutoffs[j]
  filterNumRej <- data.frame(theta=theta, numRej=numRej)
  filterTheta <- theta[j]

  return(list(padj=padj, filterThreshold=filterThreshold, filterTheta=filterTheta, filterNumRej = filterNumRej, lo.fit=lo.fit, alpha=alpha))

}

#' Estimate ZINB count component posterior probabilities
#'
#' Estimate posterior probabilities to belong to the count component according to a zero-inflated negative binomial (ZINB) model. Internally, edgeR is used for the estimation of the NB component.
#'
#' @param counts A count matrix with feature-wise expression values. Values in this matrix must be integers.
#' @param design Design matrix specifying the experimental design.
#' @param designZI The design for the zero-excess model. If \code{NULL}, the effective library size (defined as the sequencing depth multiplied by the normalization factors) is used by default.
#' @param normalization The normalization method to use. Can be one of \code{"TMM"}, \code{"DESeq2"} or \code{"phyloseq"}. If none of the methods are of interest, global normalization factors can also be given as input in the \code{normFactors} argument.
#'      If \code{"TMM"}, the trimmed mean of M-values (Robinson & Oshlack, 2010) normalization is used as implemented in \code{edgeR}.
#'      If \code{"DESeq2"}, the default median-of-ratios method from the \code{DESeq2} package (Love et al., 2014) is used for normalization.
#'      If \code{"DESeq2_poscounts"}, an adapted median-of-ratios method now implemented in \code{DESeq2}. This method was originally first implemented in the \code{phyloseq} package (McMurdie & Holmes, 2013). The adaptation ensures that genes with zero counts can be used for the purpose of normalization.
#' @param normFactors A vector of user-supplied global normalization factors for every sample. The normalization factors should be sorted according to the samples in the count matrix.
#' @param colData Only applicable if \code{normalization="DESeq2"} or \code{normalization="phyloseq"}. The \code{colData} with pheno data for constructing a \code{\link[DESeq2]{DESeqDataSet-class}} object.
#' @param designFormula Only applicable if \code{normalization="DESeq2"} or \code{normalization="phyloseq"}. The design formula required for constructing a \code{\link[DESeq2]{DESeqDataSet-class}} object.
#' @param maxit The number of iterations for the EM-algorithm. 200 by default, but larger may be useful for large datasets (many samples). Convergence of the posterior probabilities can be checked by following the distribution of posterior probabilities over iterations with \code{plotW}. The EM-algorithm will automatically stop if convergence is achieved before the maximum number of iterations.
#' @param plot Logical. Should the BCV plot be plotted in every iteration?
#' @param plotW Logical. Should the distribution of posterior probabilities for all zeros in the count matrix be plotted in every iteration?
#' @name zeroWeightsLS
#' @rdname zeroWeightsLS
#' @references
#'
#' Robinson MD and Oshlack A (2010). "A scaling normalization method for differential expression analysis of RNA-seq data." Genome Biology, 11, pp. 25.
#'
#' Love MI, Huber W and Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, pp. 550.
#'
#' McMurdie PJ and Holmes S (2013). "phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data." PLoS ONE, 8(4), pp. e61217.
#'
#' @examples
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' zeroWeights=zeroWeightsLS(counts=islam, design=design, maxit=200)
#' @export
zeroWeightsLS <- function(counts, design, maxit=200, normalization="TMM", colData, designFormula, normFactors=NULL, plot=FALSE, plotW=FALSE, verbose=TRUE, designZI=NULL, llTol=1e-4){

    if (!is.integer(counts)) {
      if (any(round(counts) != counts)) {
        stop("some values in assay are not integers")
      }
      message("converting counts to integer mode")
      mode(counts) <- "integer"
    }

  #### normalization
    require(edgeR)
    if(plot | plotW) par(mfrow=c(1,plot+plotW))
    counts <- DGEList(counts)
    if(normalization=="TMM"){
      counts <- suppressWarnings(edgeR::calcNormFactors(counts))
    } else if(normalization=="DESeq2"){
      dse = DESeqDataSetFromMatrix(counts$counts, colData=colData, design=designFormula)
      dse = DESeq2::estimateSizeFactors(dse)
      counts$samples$norm.factors = 1/dse$sizeFactor
    } else if(normalization=="DESeq2_poscounts"){
      dse = DESeqDataSetFromMatrix(counts$counts, colData=colData, design=designFormula)
      dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
      counts$samples$norm.factors = 1/dse$sizeFactor
    }
    # user input normalization factors
    if(!is.null(normFactors)) counts$samples$norm.factors = normFactors
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)

    ## initialize EM
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts), dimnames=list(c(1:nrow(counts)), NULL))
    ## starting values based on P(zero) in the library
    for(k in 1:ncol(w)) w[counts$counts[,k]==0,k] <- 1-mean(counts$counts[,k]==0)
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))
    converged=FALSE
    j=0

    for(i in 1:maxit){
      j=j+1
      zeroId <- counts$counts==0
      counts$weights <- w

      ### M-step counts
      #only estimate dispersions after weight convergence
      if(i==1 | converged){
        counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
        counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
      }
      if(plot) plotBCV(counts)
      fit <- glmFit(counts, design)
      likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)

      ### M-step mixture parameter: model zero probability
      successes <- colSums(1-w) #P(zero)
      failures <- colSums(w) #1-P(zero)
      if(is.null(designZI)){
        zeroFit <- suppressWarnings(glm(cbind(successes,failures) ~ logEffLibSize, family="binomial"))} else{
          zeroFit <- suppressWarnings(glm(cbind(successes,failures) ~-1+designZI, family="binomial"))}
      pi0Hat <- predict(zeroFit,type="response")

      ## E-step: Given estimated parameters, calculate expected value of weights
      pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
      w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)

      ## data log-likelihood
      if(i>1) llOld=ll
      ll <- log(pi0HatMat*zeroId + (1-pi0HatMat)*likC)

      delta <- (rowSums(ll)-rowSums(llOld))/(rowSums(llOld)+llTol)
      if(mean(abs(delta) < llTol)>.999){ #if 99.9% has converged
        if(j==1 & mean(abs(delta) < llTol)>.999){ #final convergence?
          cat(paste0("converged. \n")) ; return(w)}
        j=0
        converged=TRUE} else {converged=FALSE}
      if(verbose) cat(paste0("iteration: ",i,". mean conv.: ",round(mean(abs(delta) < llTol),5),"\n"))
      if(plotW) hist(w[zeroId],main=paste0("iteration: ",i,". mean conv.: ",mean(abs(delta) < llTol)))
    }
    return(w)
  }



#' Estimate weights and weighted dispersions in edgeR analysis
#'
#' This function estimates posterior probabilities to belong to the count component in a zero-inflated negative binomial (ZINB) model and integrates these in dispersion estimation with \code{\link[edgeR]{estimateDisp}}.
#'
#' @param d DGEList object containing the gene expression counts.
#' @param design The design matrix specifying the experimental design.
#' @param maxit The number of iterations for the EM-algorithm. 30 by default but larger may be useful for large datasets (many samples). Convergence of the posterior probabilities can be checked by following the distribution of posterior probabilities over iterations with \code{plotW}.
#' @param plot Logical. Should the BCV plot be plotted in every iteration?
#' @param plotW Logical. Should the distribution of posterior probabilities for all zeros in the count matrix be plotted in every iteration?
#' @param designZI. The design for the zero-excess model. If \code{NULL}, the effective library size is used by default.
#' @param ... Other arguments passed on to a weighted version of \code{\link[edgeR]{estimateDisp}}
#' @seealso
#' \code{\link[edgeR]{calcNormFactors}},  \code{\link[edgeR]{DGEList}}, \code{\link[edgeR]{estimateDisp}}
#' @examples
#' #' library(edgeR)
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' d=DGEList(islam)
#' d=calcNormFactors(d)
#' d=estimateWeightedDispersions(d,design, maxit=200)
#' fit=glmFit(d,design)
#' lrt=glmWeightedF(fit,coef=2)
#' @name estimateWeightedDispersions
#' @rdname estimateWeightedDispersions
#' @export
estimateWeightedDispersions <- function(d, design, maxit=30, plot=FALSE, plotW=FALSE, verbose=TRUE, designZI=NULL, llTol=1e-4, normalization="TMM", normFactors=NULL, ...){
  #estimate weights
  #use in .estimateDispWeighted
  #return DGEList for edgeR that can be used as input in weighted glmLRT.
  require(edgeR)
  weights <- zeroWeightsLS(counts=d$counts, design=design, maxit=maxit, plot=plot, plotW=plotW, designZI=designZI, verbose=verbose, llTol=llTol, normalization=normalization, normFactors=normFactors)
  d$weights <- weights
  d <- .estimateDispWeighted(d,design, ...)
  return(d)
}



#' Zero-inflation adjusted statistical tests for assessing differential expression.
#'
#' This function recycles an old version of the \code{\link[edgeR]{glmLRT}} method that allows an F-test with adjusted denominator degrees of freedom to account for the downweighting in the zero-inflation model.
#'
#' @param fit a \code{\link[edgeR]{DGEGLM-class}} object, usually output from \code{\link[edgeR]{glmFit}}.
#' @param coef integer or character vector indicating which coefficients of the linear model are to be tested equal to zero. Values must be columns or column names of design. Defaults to the last coefficient. Ignored if \code{contrast} is specified.
#' @param contrast numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero. Number of rows must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.
#' @param ZI Logical, specifying whether the degrees of freedom in the statistical test should be adjusted according to the weights in the \code{fit} object to account for the downweighting. Defaults to TRUE and this option is highly recommended.
#' @references
#' McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40, 4288-4297.
#' @seealso \code{\link[edgeR]{glmLRT}}
#' @examples
#' library(edgeR)
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' d=DGEList(islam)
#' d=calcNormFactors(d)
#' d=estimateWeightedDispersions(d,design, maxit=200)
#' fit=glmFit(d,design)
#' lrt=glmWeightedF(fit,coef=2)
#' @name glmWeightedF
#' @rdname glmWeightedF
#' @export
glmWeightedF <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, test="F", ZI=TRUE, independentFiltering=TRUE, filter=NULL)
  ## original function obtained from https://github.com/Bioconductor-mirror/edgeR/blob/release-3.0/R/glmfit.R
  #	Tagwise likelihood ratio tests for DGEGLM
  #	Gordon Smyth, Davis McCarthy and Yunshun Chen.
  #	Created 1 July 2010.  Last modified 22 Nov 2013.
{
  #	Check glmfit
  if(!is(glmfit,"DGEGLM")) {
    if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
      stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
    }
    stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
  }
  if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
  nlibs <- ncol(glmfit)

  #	Check test
  test <- match.arg(test,c("F","f","chisq"))
  if(test=="f") test <- "F"

  #	Check design matrix
  design <- as.matrix(glmfit$design)
  nbeta <- ncol(design)
  if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
  coef.names <- colnames(design)

  #	Evaluate logFC for coef to be tested
  #	Note that contrast takes precedence over coef: if contrast is given
  #	then reform design matrix so that contrast of interest is last column.
  if(is.null(contrast)) {
    if(length(coef) > 1) coef <- unique(coef)
    if(is.character(coef)) {
      check.coef <- coef %in% colnames(design)
      if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
      coef.name <- coef
      coef <- match(coef, colnames(design))
    }
    else
      coef.name <- coef.names[coef]
    logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
  } else {
    contrast <- as.matrix(contrast)
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    if(ncontrasts==0) stop("contrasts are all zero")
    coef <- 1:ncontrasts
    if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[coef]]
    logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
    if(ncontrasts>1) {
      coef.name <- paste("LR test of",ncontrasts,"contrasts")
    } else {
      contrast <- drop(contrast)
      i <- contrast!=0
      coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
    }
    Dvec <- rep.int(1,nlibs)
    Dvec[coef] <- diag(qrc$qr)[coef]
    Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
    design <- design %*% Q
  }
  if(length(coef)==1) logFC <- as.vector(logFC)

  #	Null design matrix
  design0 <- design[,-coef,drop=FALSE]

  #	Null fit
  fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,prior.count=0)

  #	Likelihood ratio statistic
  LR <- fit.null$deviance - glmfit$deviance
  ### ADDED
  if(ZI) fit.null$df.residual <- rowSums(fit.null$weights)-ncol(design0)
  if(ZI) glmfit$df.residual <- rowSums(glmfit$weights)-ncol(design)
  ## END ADDED
  df.test <- fit.null$df.residual - glmfit$df.residual ## okay

  #	Chisquare or F-test
  # LRT.pvalue <- switch(test,
  #                      "F" = {
  #                        phi <- quantile(glmfit$dispersion,p=0.5)
  #                        mu <- quantile(glmfit$fitted.values,p=0.5)
  #                        gamma.prop <- (phi*mu/(1 + phi*mu))^2
  #                        prior.df <- glmfit$prior.df
  #                        if(is.null(prior.df)) prior.df <- 20
  #                        glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
  #                        pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
  #                      },
  #                      "chisq" = pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)
  # )
  LRT.pvalue <- {
                         phi <- quantile(glmfit$dispersion,p=0.5)
                         mu <- quantile(glmfit$fitted.values,p=0.5)
                         gamma.prop <- (phi*mu/(1 + phi*mu))^2
                         prior.df <- glmfit$prior.df
                         if(is.null(prior.df)) prior.df <- 20
                         glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
                         pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
                       }

  rn <- rownames(glmfit)
  if(is.null(rn))
    rn <- 1:nrow(glmfit)
  else
    rn <- make.unique(rn)
  tab <- data.frame(
    logFC=logFC,
    logCPM=glmfit$AveLogCPM,
    LR=LR,
    PValue=LRT.pvalue,
    row.names=rn
  )
  glmfit$counts <- NULL
  glmfit$table <- tab
  glmfit$comparison <- coef.name
  glmfit$df.test <- df.test
  res <- new("DGELRT",unclass(glmfit))
  if(independentFiltering){
    if(is.null(filter)) filter=rowMeans(glmfit$fitted.values) #aprrox. linear w\ basemean
    res <- independentFiltering(res,filter=filter, objectType="edgeR")
  } else return(res)
}

#' Perform independent filtering in differential expression analysis.
#'
#' This function uses the \code{DESeq2} independent filtering method to increase detection power in high throughput gene expression studies.
#'
#' @param object Either a \code{\link[edgeR]{DGELRT-class}} object or a \code{\link{data.frame}} with differential expression results.
#' @param filter The characteristic to use for filtering, usually a measure of normalized mean expression for the features.
#' @param objectType Either \code{"edgeR"} or \code{"limma"}. If \code{"edgeR"}, it is assumed that \code{object} is of class \code{\link[edgeR]{DGELRT-class}}, the output of \code{\link[edgeR]{glmLRT}}. If \code{"limma"}, it is assumed that \code{object} is a \code{\link{data.frame}} and the output of a limma-voom analysis.
#' @seealso \code{\link[DESeq2]{results}}
#' @references
#' Michael I Love, Wolfgang Huber, and Simon Anders. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12):550, dec 2014.
#' @name independentFiltering
#' @rdname independentFiltering
#' @examples
#' library(limma)
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' d=DGEList(islam)
#' nf=calcNormFactors(islam)
#' y=zeroWeightedVoom(d,design,nf=nf,maxit=200)
#' fit=lmWeightedFit(y,design)
#' fit=eBayes(fit)
#' tt=topTable(fit,coef=2,sort.by="none",number=nrow(fit))
#' baseMean=unname(rowMeans(sweep(d$counts,2,nf,FUN="*")))
#' ttFiltered=independentFiltering(tt,filter=baseMean, objectType="limma")
#' @export
independentFiltering <- function(object, filter, objectType=c("edgeR","limma")){
  if(objectType=="edgeR"){
    hlp <- .pvalueAdjustment_kvdb(filter=filter, pValue=object$table$PValue)
    object$table$padjFilter <- hlp$padj
    return(object)
  } else if(objectType=="limma"){
    hlp <- .pvalueAdjustment_kvdb(filter=filter, pValue=object$P.Value)
    object$padjFilter <- hlp$padj
    return(object)
  } else stop("objectType must be either one of 'edgeR' or 'limma'.")
}
