% --- Source file: ASSET.Rd ---
\name{ASSET}
\alias{ASSET}
\docType{package}
\title{
Association analysis for SubSETs
}
\description{
ASSET is a suite of statistical tools specifically designed to be powerful for pooling association signals across multiple studies when true effects may exist only in a subset of the studies and possibly in opposite directions across studies. The method explores all possible subsets (or a restricted set if the user specifies so) of studies and evaluates fixed-effect meta-analysis-type test-statistics for each subset. The final test-statistic is obtained by maximizing the subset-specific test-statistics over all possible subsets and then evaluating its significance after efficient adjustment for multiple-testing, taking into account the correlation between test-statistics across different subsets due to overlapping subjects. The method not only returns a p-value for significance for the overall evidence of association of a SNP across studies, but also outputs the "best subset" containing the studies that contributed to the overall association signal. For detection of association signals with effects in opposite directions, ASSET allows subset search separately for positively- and negatively- associated studies and then combines association signals from two directions using a chi-square test-statistic. The method can take into account correlation due to overlapping subjects across studies (e.g. shared controls). Although the method is originally developed for conducting genetic association scans, it can also be applied for analysis of non-genetic risk factors. 
Since version 2.0.0, a new option has been added to enable
pre-screening of studies based on marginal association tests followed by subset-based 
meta-analysis to achieve an increase computational speed while analyzing a large number of studies.
}
\details{
\bold{The new version includes option to pre-screen phenotypes based on significance of marginal association tests. 
It increases computational speed for ASSET when analyzing larger number of phenotypes.}

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
For example, the default options of these functions currently assume all possible subsets are to be searched. 
  For analysis of case-control studies with ordered 
diseased subtypes (e.g. stages of a cancer), however, it may be more meaningful to restrict the subset search to incorporate ordering 
constraints among the disease subtypes. In such a situation, one can pass a function argument \code{sub.def} to \code{z.max} and \code{p.dlm} for performing restricted
subset searches.
}
\author{Samsiddhi Bhattacharjee, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}
\references{Samsiddhi Bhattacharjee, Preetha Rajaraman, Kevin B. Jacobs, William A. Wheeler, Beatrice S. Melin, Patricia Hartge, GliomaScan Consortium, Meredith Yeager, Charles C. Chung, Stephen J. Chanock, Nilanjan Chatterjee. A subset-based approach improves power and interpretation for combined-analysis of genetic association studies of heterogeneous traits. Am J Hum Genet, 2012, 90(5):821-35 \cr

Bhattacharjee et al. Pre-screening and Meta-analysis based on SubSETs: A Fast and Powerful Approach
 to Pleiotropic Analysis Across a Large Number of Traits (In Preparation)
}

