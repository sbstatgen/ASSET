% --- Source file: man/h.summary.Rd ---
\name{h.summary}
\alias{h.summary}
\title{Summary results from subset-search.}
\description{
This function produces summary results from subset-based association analysis.
}
\usage{
h.summary(rlist, level = 0.05, digits = 3)
}
\arguments{
  \item{rlist}{List returned by \code{\link{h.traits}} or \code{\link{h.types}} }
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} confidence intervals. The confidence intervals are obtained by inverting the 
corresponding multiple-testing adjusted p-values.}
  \item{digits}{Number of significant digits to retain in odds ratios and confidence intervals in the summary table}
}
\value{
A list of data frames, one for each of the methods specified the original call of the functions \code{\link{h.traits}} or \code{\link{h.types}}. Each row of a data frame corresponds to a SNP and 
the values include p-values for overall association (including component-wise p-values for two-sided search), names of phenotypes or disease subtypes included in
the best-subset, summary regression coefficients (e.g. log-odds-ratio for case-control studies) representing strength of association of a SNP with the
identified subset of traits/subtype  and corresponding confidence intervals.
}
\details{ 
Returns a list of data frames containing p-values, odds-ratios, confidence intervals and the traits/types for each analysis. 
The number of data frames in the list will depend on which function (\code{\link{h.traits}} or \code{\link{h.types}})
was called and on the function options specified. 
}
\seealso{\code{\link{h.forestPlot}}, \code{\link{h.traits}}, \code{\link{h.types}} }

\examples{
 # Use the example data
 data(ex_trait, package="ASSET")

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
                 cor.numr=FALSE, search=NULL, side=2, meta=TRUE, 
                 zmax.args=NULL, meth.pval="DLM")
 
 h.summary(res)

}
