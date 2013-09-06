% --- Source file: ASSET.Rd ---
\name{ASSET}
\alias{ASSET}
\docType{package}
\title{
Association analysis for SubSETs
}
\description{
This package is for subset-based association analysis of heterogeneous but possibly related traits. 
}
\details{
The package consists of two main functions: (1) \code{\link{h.traits}} and (2) \code{\link{h.types}}. The function \code{\link{h.traits}}
is suitable for conducting meta-analysis of studies of possibly different traits when summary level data are available from individual
studies. The function allows correlation among different studies/traits, which, for example, may arise due to shared subjects across studies.
The function can also be used to conduct "meta-analysis" across multiple correlated traits on the same individuals by appropriately specifying the correlation
matrix for the multivariate trait. The method, however, is not optimized yet (from a power perspective) for analyzing multivariate traits measured on the same 
individuals. The function \code{\link{h.types}} is suitable for analysis of case-control studies when cases consist of distinct disease subtypes. This
function assumes individual level data are available. 
The functions \code{\link{h.summary}} and \code{\link{h.forestPlot}} are useful for summarizing results and displaying forest plots. 
The helper functions
\code{\link{z.max}} and \code{\link{p.dlm}} are generic functions called internally for obtaining the maximized subset-based test-statistics and the corresponding
p-values approximated by the Discrete Local Maximization (DLM) method. These functions can be further customized for specific applications. 
For example, the default options of these functions currently assume all possible subsets to search. For analysis of case-control studies with ordered 
diseased subtypes (e.g. stages of a cancer), however, it may be more meaningful to restrict the subset search to incorporate ordering 
constraints among the disease subtypes. In such a situation, one can pass a function argument \code{sub.def} to \code{z.max} and \code{p.dlm} for performing restricted
subset searches.
}
\author{Samsiddhi Bhattacharjee, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}
\references{ Bhattacharjee S, Chatterjee N and others. A subset-based approach improves power and interpretation for combined-analysis of genetic association studies of heterogeneous traits. Am J Hum Genet, 2012, 90(5):821-35 \cr
}

