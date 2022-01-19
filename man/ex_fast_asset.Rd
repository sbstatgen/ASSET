% --- Source file: ex_trait.Rd ---
\name{ex_fast_asset}
\alias{ex_fast_asset}
\alias{betahat}
\alias{SE}
\alias{traits}
\alias{SNP}
\alias{Neff}
\alias{ldscintmat}


%\docType{package}
\title{
  Data for the fast_asset example
}
\description{
 Data for \code{\link{fast_asset}}
}
\details{ The object data contains estimated log odds-ratios and their standard errors
  for 1 snp and 116 traits in \code{betahat} and \code{SE}. The trait names are in \code{traits}
  and the snp name in \code{SNP}. The \code{Neff} vector contains effective sample sizes of the 
  individual traits. The \code{ldscintmat} matrix contains correlation matrix between pairs of 
  traits as estimated using genome-wide LD-score regression.
}
\examples{

 data(ex_fast_asset, package="ASSET")

 # Display the parts of the data objects
 head(betahat)
 head(SE)
 head(traits)
 SNP
 head(Neff)
 ldscintmat[1:6, 1:6]

}
\keyword{fast_asset}
\seealso{h.traits, h.types}

