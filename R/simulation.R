
.getParamsZTNB <- function(counts, offset, design=NULL) {

  libSize=offset
  #fit a ZTNB model only on positive counts part
  countsModel = counts[counts>0]
  if(length(countsModel)<2) stop("Need at least two positive counts")
  libSizeModel = libSize[counts>0]
  if(is.null(design)){designFit=matrix(1,nrow=length(countsModel),ncol=1)} else {designFit=design[counts>0,]}
  fit=suppressWarnings(try(gamlss(formula=countsModel~-1+designFit+offset(log(libSizeModel)), family="NBIZeroTruncated", control=gamlss.control(trace=FALSE, n.cyc=300)),silent=TRUE))


  if(class(fit)[1]=="try-error") return(c(dispersion=NA, lambda=NA))
  #lambda=exp(fit$mu.coefficients)
  lambda=mean(countsModel/libSizeModel)
  dispersion=exp(fit$sigma.coefficients)
  return(c(dispersion=dispersion,lambda=lambda))
}

#' Estimate feature-wise parameters according to a ZTNB distribution.
#'
#' @param counts the expression matrix
#' @export
getDatasetZTNB = function(counts, design, drop.extreme.dispersion = FALSE){

  #loadNamespace("gamlss.tr")
  #loadNamespace("gamlss")
  #loadNamespace("mgcv")
  #loadNamespace("MASS")
  #loadNamespace("edgeR")

  #create ZTNB distribution
  suppressPackageStartupMessages({library(gamlss)}) #for loading NBI family.
  #gamlss.tr::gen.trun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)

  #### estimate lambda and overdispersion based on ZTNB.
  d <- edgeR::DGEList(counts)
  cp <- edgeR::cpm(d,normalized.lib.sizes=TRUE)
  dFiltered=d
  dFiltered <- edgeR::calcNormFactors(dFiltered)
  dFiltered$AveLogCPM <- edgeR::aveLogCPM(dFiltered)
  params=t(apply(dFiltered$counts,1,function(x) .getParamsZTNB(counts=x,offset=dFiltered$samples$lib.size, design=design)))
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
  zeroFit=mgcv::gam(expectCounts~s(midsHlp)+logLibHlp,family=binomial)

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



