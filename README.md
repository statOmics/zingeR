zingeR: Zero-Inflated Negative binomial Gene Expression in R.

Note, that the simulation framework in zingeR has been updated, and we recommend interested users to use the updated simulation framework, which can be found on the GitHub repository of our updated manuscript: https://github.com/statOmics/zinbwaveZinger. We will therefore no longer support the simulation framework present in this repository.

For installation, please make sure you have recent versions of edgeR (v >= 3.19) and DESeq2 (v >= 1.17.6) installed along with their required dependencies. You should be able to install zingeR and its dependencies through:

```R
source("https://bioconductor.org/biocLite.R")
install.packages("Formula")
biocLite("DESeq2") #install DESeq2 dependencies
library(devtools)
devtools::install_git("https://git.bioconductor.org/packages/DESeq2”) #install most recent DESeq2 version
biocLite("edgeR") #install edgeR dependencies
devtools::install_git("https://git.bioconductor.org/packages/edgeR") #install most recent edgeR version
install_github("statOmics/zingeR")


```

For bug reports, issues or requested extensions, please raise an issue on the GitHub or send an e-mail to koen.vandenberge@ugent.be.
