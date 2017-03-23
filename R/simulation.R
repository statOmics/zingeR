.mytrun <- function (par = c(0), family = "NO",
                     name = "tr", type = c("left", "right", "both"),
                     varying = FALSE, ...)
{
  ## Function adapted from gen.trun because of issue with namespaces.
  ## Function written by Joris Meys on 2017/02/24.

  type <- match.arg(type)
  fam <- as.gamlss.family(family)
  fname <- fam$family[[1]]
  dfun <- paste(paste("d", fname, sep = ""), name, sep = "")
  pfun <- paste(paste("p", fname, sep = ""), name, sep = "")
  qfun <- paste(paste("q", fname, sep = ""), name, sep = "")
  rfun <- paste(paste("r", fname, sep = ""), name, sep = "")
  fun <- paste(fname, name, sep = "")
  alldislist <- c(dfun, pfun, qfun, rfun, fun)

  d.f <- trun.d(par, family = fname, type = type, varying = varying,
                ...)
  p.f <- trun.p(par, family = fname, type = type, varying = varying,
                ...)
  q.f <- trun.q(par, family = fname, type = type, varying = varying,
                ...)
  r.f <- trun.r(par, family = fname, type = type, varying = varying,
                ...)
  f <- trun(par, family = fname, type = type,
            name = name, local = FALSE, varying = varying, ...)

  out <- list(d.f, p.f, q.f, r.f, f)
  names(out) <- alldislist
  return(out)
}



.getParamsZTNB <- function(counts, offset, design=NULL, ztnbFun=NULL) {
  libSize=offset
  #fit a ZTNB model only on positive counts part
  countsModel = counts[counts>0]
  if(length(countsModel)<2) stop("Need at least two positive counts")
  libSizeModel = libSize[counts>0]
  if(is.null(design)){designFit=matrix(1,nrow=length(countsModel),ncol=1)} else {designFit=design[counts>0,]}
  fit=suppressWarnings(try(gamlss(formula=countsModel~-1+designFit+offset(log(libSizeModel)), family=ztnbFun$NBIZeroTruncated, control=gamlss.control(trace=FALSE, n.cyc=300)),silent=TRUE))

  if(class(fit)[1]=="try-error") return(c(dispersion=NA, lambda=NA))
  #lambda=exp(fit$mu.coefficients)
  lambda=mean(countsModel/libSizeModel)
  #lambda=exp(mean(fit$mu.coefficients)) #geometric mean
  dispersion=exp(fit$sigma.coefficients)
  return(c(dispersion=dispersion,lambda=lambda))
}

#' Estimate feature-wise parameters according to a ZTNB distribution.
#'
#' This function fits a zero-truncated negative binomial distribution to the positive counts for each feature (row). Before running this function, a zero-truncated negative binomial distribution family must be created with the \code{gamlss} package as: \code{gamlss.tr::gen.trun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)}
#'
#' @param counts A numeric matrix containing gene expression counts.
#' @param design The design of the experiments with rows corresponding to samples and columns corresponding to coefficients.
#' @param drop.extreme.dispersion Either a numeric value between $0$ and $1$, stating the proportion of genes with extreme (high) dispersions to remove for simulation, or FALSE (default), if no dispersions should be removed for the analysis.
#' @param offset The offset to use (typically the sequencing depth) when estimating gene-wise means and dispersions in the zero-truncated negative binomial model. These parameters will be used as a basis for the simulation.
#' @seealso \code{\link[zingeR]{NBsimSingleCell}}
#' @examples
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' params = getDatasetZTNB(counts=islam, design=design)
#' @name getDatasetZTNB
#' @rdname getDatasetZTNB
#' @export
getDatasetZTNB = function(counts, design, drop.extreme.dispersion = FALSE, offset=NULL){

  #create ZTNB distribution
  ztnbFun = .mytrun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)

  #### estimate lambda and overdispersion based on ZTNB.
  d <- edgeR::DGEList(counts)
  cp <- edgeR::cpm(d,normalized.lib.sizes=TRUE)
  dFiltered=d
  dFiltered <- edgeR::calcNormFactors(dFiltered)
  dFiltered$AveLogCPM <- edgeR::aveLogCPM(dFiltered)
  if(is.null(offset)) offset=dFiltered$samples$lib.size
  library(gamlss) ; library(gamlss.tr)
  params=t(apply(dFiltered$counts,1,function(x) .getParamsZTNB(counts=x, offset=offset, design=design, ztnbFun=ztnbFun)))
  rmRows = which(params[,2]>1) #impossibly high lambda
  rmRows2 = which(params[,2]==0) #zero lambda
  naRows = which(apply(params,1, function(row) any(is.na(row)))) #not fitted
  nonZeroDispRows = which(params[,1]<0 | params[,1]==0) #negative dispersion
  throwRows = c(rmRows,rmRows2,naRows,nonZeroDispRows)
  params = params[-throwRows,]

  ### estimate logistic GAM P(zero) ~ s(aveLogCPM) + logLibSize
  ### use unfiltered data for this model.
  propZero = colMeans(counts==0)
  propZeroGene = rowMeans(counts==0)
  d <- edgeR::DGEList(counts)
  d <- edgeR::calcNormFactors(d)
  avCpm <- edgeR::aveLogCPM(d)
  cpmHist = hist(avCpm, breaks=150, plot=FALSE)
  breaks = cpmHist$breaks
  mids = cpmHist$mids
  midsHlp=rep(mids,ncol(d$counts))
  logLibSize = log(colSums(counts))
  logLibHlp=rep(logLibSize,each=length(mids))
  binHlp=sapply(breaks[-length(breaks)],function(x) avCpm>x)
  binId=apply(binHlp,1,function(x) max(which(x)))
  nonNullCounts = t(sapply(1:length(mids), function(bin){
    binRows <- binId==bin
    if(sum(binRows)==0) rep(0,ncol(counts)) else
      if(sum(binRows)==1) (counts[which(binRows),]!=0)*1 else
        colSums(counts[which(binRows),]!=0)
  }))
  nullCounts = t(sapply(1:length(mids), function(bin){
    binRows <- binId==bin
    if(sum(binRows)==0) rep(0,ncol(counts)) else
      if(sum(binRows)==1) (counts[which(binRows),]==0)*1 else
        colSums(counts[which(binRows),]==0)
  }))
  expectCounts=cbind(c(nullCounts),c(nonNullCounts))
  #zeroFit=mgcv::gam(expectCounts~s(midsHlp)+logLibHlp,family=binomial)
  zeroFit=mgcv::gam(expectCounts~s(midsHlp,by=logLibHlp),family=binomial)

  ### drop extreme dispersions
  dFiltered$AveLogCPM <- edgeR::aveLogCPM(dFiltered)
  dFiltered$AveLogCPM <- dFiltered$AveLogCPM[-throwRows]
  propZeroGene = propZeroGene[-throwRows]
  params=data.frame(dispersion=params[,1], lambda=params[,2], aveLogCpm=dFiltered$AveLogCPM, propZeroGene=propZeroGene)
  dispersion <- params$dispersion
  AveLogCPM <- params$aveLogCpm
  lambda <- params$lambda
  propZeroGene <- params$propZeroGene

  if(is.numeric(drop.extreme.dispersion))
  {
    bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE, na.rm=TRUE)
    ids <- dispersion <= bad
    AveLogCPM <- AveLogCPM[ids]
    dispersion <- dispersion[ids]
    lambda <- lambda[ids]
    propZeroGene <- propZeroGene[ids]
    params <- params[ids,]
    dFiltered <- dFiltered[ids,]
  }
  #lambda=lambda/sum(lambda) #make sure they sum to 1
  dataset.AveLogCPM <- AveLogCPM
  dataset.dispersion <- dispersion
  dataset.lambda <- lambda
  dataset.propZeroGene <- propZeroGene
  dataset.lib.size <- d$samples$lib.size
  dataset.nTags <- nrow(d)
  list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags, dataset.propZeroFit=zeroFit, dataset.lambda=lambda, dataset.propZeroGene=propZeroGene, dataset.breaks = breaks, dataset.cpm=cpm)
}

#' Simulate a scRNA-seq experiment
#'
#' This function simulates a count matrix for a scRNA-seq experiment based on parameters estimated from a real scRNA-seq count matrix. Global patterns of zero abundance as well as feature-specific mean-variance relationships are retained in the simulation.
#'
#' @param dataset An expression matrix representing the dataset on which the simulation is based.
#' @param group Group indicator specifying the attribution of the samples to the different conditions of interest that are being simulated.
#' @param nTags The number of features (genes) to simulate. $1000$ by default
#' @param nlibs The number of samples to simulate. Defaults to \code{length(group)}.
#' @param lib.size The library sizes for the simulated samples. If \code{NULL} (default), library sizes are resampled from the original datset.
#' @param pUp Numeric value between $0$ and $1$ ($0.5$ by default) specifying the proportion of differentially expressed genes that show an upregulation in the second condition.
#' @param foldDiff The fold changes used in simulating the differentially expressed genes. Either one numeric value for specifying the same fold change for all DE genes, or a vector of the same length as \code{ind} to specify fold changes for all differentially expressed genes. Note that fold changes above $1$ should be used as input of which a fraction will be inversed (i.e. simulation downregulation) according to `pUp`. Defaults to $3$.
#' @param verbose Logical, stating whether progress be printed.
#' @param ind Integer vector specifying the rows of the count matrix that represent differential features.
#' @param params An object containing feature-wise parameters used for simulation as created by \code{\link[zingeR]{getDatasetZTNB}}. If \code{NULL}, parameters are estimated from the dataset provided.
#' @param  randomZero A numeric value between $0$ and $1$ specifying the random fraction of cells that are set to zero after simulating the expression count matrix. Defaults to $0$.
#' @param min.dispersion The minimum dispersion value to use for simulation. $0.1$ by default.
#' @param max.dipserion The maximum dispersion value to use for simulation. $400$ by default.
#' @param drop.extreme.dispersion Only applicable if \code{params=NULL} and used as input to \code{\link[zingeR]{getDatasetZTNB}}. Numeric value between $0$ and $1$ specifying the fraction of highest dispersion values to remove after estimating feature-wise parameters.
#' @seealso \code{\link[zingeR]{getDatasetZTNB}}
#' @examples
#' data(islamEset,package="zingeR")
#' islam=exprs(islamEset)[1:2000,]
#' design=model.matrix(~pData(islamEset)[,1])
#' gamlss.tr::gen.trun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)
#' params = getDatasetZTNB(counts=islam, design=design)
#' nSamples=80
#' grp=as.factor(rep(0:1, each = nSamples/2)) #two-group comparison
#' nTags=2000 #nr of features
#' set.seed(436)
#' DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #10% differentially expressed
#' fcSim=(2 + rexp(length(DEind), rate = 1/2)) #fold changes
#' libSizes=sample(colSums(islam),nSamples,replace=TRUE) #library sizes
#' simDataIslam <- NBsimSingleCell(foldDiff=fcSim, ind=DEind, dataset=islam, nTags=nTags, group=grp, verbose=TRUE, params=params, lib.size=libSizes)
#' @export
NBsimSingleCell <- function(dataset, group, nTags = 10000, nlibs = length(group), lib.size = NULL, drop.extreme.dispersion = FALSE, pUp=.5, foldDiff=3, verbose=TRUE, ind=NULL, params=NULL, randomZero=0, max.dispersion=400, min.dispersion=0.1)
{
  require(edgeR)
  group = as.factor(group)
  expit=function(x) exp(x)/(1+exp(x))
  logit=function(x) log(x/(1-x))

  sample.fun <- function(object)
  {
    nlibs <- object$nlibs
    nTags <- object$nTags
    AveLogCPM <- object$dataset$dataset.AveLogCPM
    dispersion <- object$dataset$dataset.dispersion
    lambda <- object$dataset$dataset.lambda
    #lambda <- (2^AveLogCPM)/1e6
    propZeroGene <- dat$dataset$dataset.propZeroGene
    id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
    object$AveLogCPM <- AveLogCPM[id_r]
    Lambda <- lambda[id_r]
    #Lambda <- Lambda/sum(Lambda) #normalize so they all sum to 1
    Dispersion <- dispersion[id_r]
    Dispersion[Dispersion>max.dispersion] = max.dispersion
    Dispersion[Dispersion<min.dispersion] = min.dispersion
    propZeroGene <- propZeroGene[id_r]
    Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
    object$Lambda <- Lambda
    Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
    object$Dispersion <- Dispersion
    object$propZeroGene <- propZeroGene
    object
  }
  diff.fun <- function(object)
  {
    group <- object$group
    pUp <-  object$pUp
    foldDiff <- object$foldDiff
    Lambda <- object$Lambda
    nTags <- object$nTags
    g <- group == levels(group)[1]
    if(length(ind)>0 & !all(foldDiff==1)) {
      fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
      Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
      Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir))
      object$Lambda <- Lambda
      object$indDE <- ind
      object$indNonDE <- (1:nTags)[-ind]
      foldDiff[fcDir==1] <- 1/foldDiff[fcDir==1]
      object$foldDiff <- foldDiff #group2 / group1
    }
    if(all(foldDiff==1)) object$indDE <- NA
    object
  }
  sim.fun <- function(object)
  {
    Lambda <- object$Lambda
    Dispersion <- object$Dispersion
    nTags <- object$nTags
    nlibs <- object$nlibs
    lib.size <- object$lib.size
    zeroFit <- dat$dataset$dataset.propZeroFit
    propZeroGene <- dat$propZeroGene
    propZeroGene[propZeroGene==1] <- 1-1e-4
    propZeroGene[propZeroGene==0] <- 1e-4
    design <- object$design
    avLogCpm <- object$AveLogCPM
    mids <- object$dataset$dataset.mids
    breaks <- object$dataset$dataset.breaks
    ## get matrix of zero probabilities
    libPredict=rep(log(lib.size),each=length(avLogCpm))
    cpmPredict=rep(avLogCpm,length(lib.size))
    ## no noise
      zeroProbMatLink = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="link"), byrow=FALSE, ncol=nlibs)
      meanDiff = rowMeans(zeroProbMatLink)-logit(propZeroGene)
      zeroProbMat = expit(sweep(zeroProbMatLink,1,meanDiff,"-"))

    ## introduce random zeroes for lower count genes
    zeroProbMat[sample(1:length(zeroProbMat), floor(randomZero*length(zeroProbMat)))]=1-1e-5

    ## adjust mu for adding zeroes
    mu=sweep(Lambda,2,lib.size,"*")
    mu[mu<1e-16] = 1
    adjustment = zeroProbMat*mu
    mu=mu+adjustment

    ## simulate counts acc to a zero-adjusted NB model
    counts = rZANBI(n=nTags*nlibs, mu=mu, sigma=Dispersion, nu=zeroProbMat)

    ## the rZANBI function rarely simulates Inf values for very low mu estimates. Resimulate for these genes using same params, if present
    ## also, resample features with all zero counts
    zeroCountsId <- which(rowSums(counts)==0)
    infId <- which(apply(counts,1,function(row) any(is.infinite(row))))
    while(length(zeroCountsId)>0 | length(infId)>0){
      if(length(zeroCountsId)>0){ #resimulate
        counts[zeroCountsId,] = rZANBI(n=length(zeroCountsId)*nlibs, mu=mu[zeroCountsId,], sigma=Dispersion[zeroCountsId,], nu=zeroProbMat[zeroCountsId,])
      }
      if(length(infId)>0){ #replace Inf values by resampling
        counts[infId,] <- rZANBI(n=length(infId)*nlibs, mu=mu[infId,], sigma=Dispersion[infId,], nu=zeroProbMat[infId,])
      }
      zeroCountsId <- which(rowSums(counts)==0)
      infId <- which(apply(counts,1,function(row) any(is.infinite(row))))
    }

    rownames(counts) <- paste("ids", 1:nTags, sep = "")
    object$counts <- counts
    object
  }

  if(verbose) message("Preparing dataset.\n")
  if(is.null(params)){
    dataset <- getDatasetZTNB(counts = dataset, drop.extreme.dispersion = drop.extreme.dispersion)
  } else {
    dataset <- params
  }
  dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pUp = pUp, foldDiff = foldDiff))


  if(is.null(dat$lib.size)){
    dat$lib.size <- sample(dataset$dataset.lib.size, nlibs, replace=TRUE)}
  if(is.null(nTags)) dat$nTags <- dat$dataset$dataset.nTags
  if(verbose) message("Sampling.\n")
  dat <- sample.fun(dat)
  if(verbose) message("Calculating differential expression.\n")
  dat <- diff.fun(dat)
  if(verbose) message("Simulating data.\n")
  dat <- sim.fun(dat)
  dat
}

