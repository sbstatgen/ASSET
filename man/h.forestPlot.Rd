% --- Source file: man/h.forestPlot.Rd ---
\name{h.forestPlot}
\alias{h.forestPlot}
\title{Forest plot for meta-analysis of heterogenerous traits or types.}
\description{
Forest Plot for meta-analysis of heterogenerous traits or types.
}
\usage{
h.forestPlot(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2, calc.miss=FALSE)
}
\arguments{
  \item{rlist}{The list of results returned by \code{\link{h.traits}} or \code{\link{h.types}}. SNPs other than \code{snp.var} are ignored.} 
  \item{snp.var}{A character string giving the name of the SNP variable to be plotted. No default.}
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} confidence intervals.}
  \item{p.adj}{Logical. Whether to report Bonferroni adjusted p-values for each individual subtype. Default is TRUE.}  
  \item{digits}{Number of significant digits to display the odds ratios in the plot.}
  \item{calc.miss}{Logical. If TRUE, first 3 components of rlist are calculated and filled (if NULL) so that all 3 summarizations can be 
        shown in the forest plot.}
}
\value{
Forest plot for a SNP showing regression coefficients (e.g. log-odds-ratio for case-control studies) for individual studies/traits and 
confidence intervals, estimate of an overall regression coefficient and confidence interval based on standard fixed-effect meta-analysis
and estimate of regression coefficient(s) and confidence intervals associated with the identified best subset(s).
}
\seealso{\code{\link{h.summary}}, \code{\link{h.traits}}, \code{\link{h.types}}}

\examples{
 # Use the example data
 data(ex_trait, package="ASSET")
 data

 # Define the input arguments to h.traits
 snps       <- as.vector(data[, "SNP"])
 traits.lab <- paste("Trait_", 1:6, sep="")
 beta.hat   <- as.matrix(data[, paste(traits.lab, ".Beta", sep="")])
 sigma.hat  <- as.matrix(data[, paste(traits.lab, ".SE", sep="")])
 cor        <- list(N11=N11, N00=N00, N10=N10)
 ncase      <- diag(N11)
 ncntl      <- diag(N00)

 # Now let us call h.traits on these summary data. 
 res <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=cor, 
                 cor.numr=FALSE, search=NULL, side=2, meta=TRUE, zmax.args=NULL, 
                 meth.pval="DLM")

 h.forestPlot(res, "SNP_1", digits=3)

}
