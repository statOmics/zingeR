% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{NBsimSingleCell}
\alias{NBsimSingleCell}
\title{Simulate a scRNA-seq experiment}
\usage{
NBsimSingleCell(dataset, group, nTags = 10000, nlibs = length(group),
  lib.size = NULL, drop.extreme.dispersion = FALSE, pUp = 0.5,
  foldDiff = 3, verbose = TRUE, ind = NULL, params = NULL,
  randomZero = 0, max.dispersion = 400, min.dispersion = 0.01)
}
\arguments{
\item{dataset}{An expression matrix representing the dataset on which the simulation is based.}

\item{group}{Group indicator specifying the attribution of the samples to the different conditions of interest that are being simulated.}

\item{nTags}{The number of features (genes) to simulate. $1000$ by default}

\item{nlibs}{The number of samples to simulate. Defaults to \code{length(group)}.}

\item{lib.size}{The library sizes for the simulated samples. If \code{NULL} (default), library sizes are resampled from the original datset.}

\item{drop.extreme.dispersion}{Only applicable if \code{params=NULL} and used as input to \code{\link[zingeR]{getDatasetZTNB}}. Numeric value between $0$ and $1$ specifying the fraction of highest dispersion values to remove after estimating feature-wise parameters.}

\item{pUp}{Numeric value between $0$ and $1$ ($0.5$ by default) specifying the proportion of differentially expressed genes that show an upregulation in the second condition.}

\item{foldDiff}{The fold changes used in simulating the differentially expressed genes. Either one numeric value for specifying the same fold change for all DE genes, or a vector of the same length as \code{ind} to specify fold changes for all differentially expressed genes. Note that fold changes above $1$ should be used as input of which a fraction will be inversed (i.e. simulation downregulation) according to `pUp`. Defaults to $3$.}

\item{verbose}{Logical, stating whether progress be printed.}

\item{ind}{Integer vector specifying the rows of the count matrix that represent differential features.}

\item{params}{An object containing feature-wise parameters used for simulation as created by \code{\link[zingeR]{getDatasetZTNB}}. If \code{NULL}, parameters are estimated from the dataset provided.}

\item{randomZero}{A numeric value between $0$ and $1$ specifying the random fraction of cells that are set to zero after simulating the expression count matrix. Defaults to $0$.}

\item{min.dispersion}{The minimum dispersion value to use for simulation. $0.1$ by default.}

\item{max.dipserion}{The maximum dispersion value to use for simulation. $400$ by default.}
}
\description{
This function simulates a count matrix for a scRNA-seq experiment based on parameters estimated from a real scRNA-seq count matrix. Global patterns of zero abundance as well as feature-specific mean-variance relationships are retained in the simulation.
}
\examples{
data(islamEset,package="zingeR")
islam=exprs(islamEset)[1:2000,]
design=model.matrix(~pData(islamEset)[,1])
gamlss.tr::gen.trun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)
params = getDatasetZTNB(counts=islam, design=design)
nSamples=80
grp=as.factor(rep(0:1, each = nSamples/2)) #two-group comparison
nTags=2000 #nr of features
set.seed(436)
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #10\% differentially expressed
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #fold changes
libSizes=sample(colSums(islam),nSamples,replace=TRUE) #library sizes
simDataIslam <- NBsimSingleCell(foldDiff=fcSim, ind=DEind, dataset=islam, nTags=nTags, group=grp, verbose=TRUE, params=params, lib.size=libSizes)
}
\seealso{
\code{\link[zingeR]{getDatasetZTNB}}
}

