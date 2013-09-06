% --- Source file: man/p.dlm.Rd ---
\name{p.dlm}
\alias{p.dlm}
\title{Discrete Local Maxima approximate p-value.}
\description{
Function to obtain Discrete Local Maxima based estimates of p-values for z-scores maximized over subsets (of traits or subtypes), 
with possible restrictions and weights. Should not be called directly. See details.
}
\usage{
p.dlm(t.vec, k, search, side, cor.def=NULL, cor.args=NULL, sizes=rep(1, k), 
      sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL)
}
\arguments{
  \item{t.vec}{Numeric vector of (positive) points for which to calculate p-values, i.e. general observed Z-max values. No default.}
  \item{k}{Integer (currently less than 30). The number of studies (traits) or subtypes being analyzed. No default.}  
  \item{search}{0, 1 or 2. Search option, with 0 indicating subtype analysis, 1 and 2 denote one-sided and two-sided subset-search. No default.}
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximized), \code{side} should be 2. For one-tailed tests, 
  \code{side} should be 1 (positive tail assumed). No default. Ignored when search is 2.}
  \item{cor.def}{A function with at least 3 arguments which calculates correlation between its first argument (a subset) and its second argument (subsets such
  as its neighbors). The third argument is the number of traits/subtypes and the function should return a vector of correlations with the neighbors.
 If NULL or a non-function value is specified, internal default 
  functions for the corresponding search option are used.}
  \item{cor.args}{Other arguments to be passed to \code{cor.def}. These can include sample sizes and overlaps of different studies or subtypes and analysis option
  such as case-control or case-complement that affect the correlation structure.
  If \code{cor.def} is NULL, then \code{ncase} and \code{ncntl} must be specified in this list.
  }
  \item{sizes}{Sizes of equivalence classes of traits. By default, no two traits or studies are equivalent. This argument is for internal use.}
  \item{sub.def}{A function to restrict subsets, e.g., order restrictions in subtype analysis. Should accept a subset (a logical vector of size k) as
  its first argument and should return TRUE if the subset satisfies restrictions and FALSE otherwise. Default is NULL implying all (2^k - 1) subsets
  are considered in the maximum.}
  \item{sub.args}{Other arguments to be passed to \code{sub.def} as list. Default is NULL (i.e. none).}
  \item{wt.def}{A function that gives the weight of one subset with respect to another. Should accept two subsets as the first two arguments and return
  a single positive weight. Default NULL. Currently this option is not implemented and the argument is ignored.}
  \item{wt.args}{Other arguments to be passed to \code{wt.def} as a list . Default NULL. Currently ignored.}
  }
\value{
A numeric vector of estimated p-values.
}
\details{ 
The function is vectorized to handle blocks of SNPs at a time. Currently weight options are ignored.
This is a helper function that is called internally by \code{\link{h.traits}} and \code{\link{h.types}} 
and should not be called directly. The arguments of this function that have defaults, can be customized using
the argument \code{pval.args} in \code{\link{h.traits}} and \code{\link{h.types}}.
}

\examples{

  # A function to define the correlations between a subset and its neighbors
  # Returned values should not exceed the value of 1
  cor.def <- function(subset, neighbors, k, ncase, ncntl) {
    n <- ncol(neighbors)
    mat <- matrix(subset, nrow=k, ncol=n, byrow=FALSE)
    cor <- (mat + neighbors)*(1:k)/(k^2)
    cor <- colSums(cor)
    cor <- cor/max(cor)
    dim(cor) <- c(n, 1)

    cor  
  }

  # Subset definition
  sub.def <- function(logicalVec) {
    # Only allow the cummulative subsets:
    # TRUE FALSE FALSE FALSE ...
    # TRUE TRUE FALSE FALSE ...
    # TRUE TRUE TRUE FALSE ...
    # etc
    sum <- sum(logicalVec)  
    ret <- all(logicalVec[1:sum])

    ret
  }

  k     <- 5
  t.vec <- 1:k

  p.dlm(t.vec, k, 1, 2, cor.def=cor.def, sub.def=sub.def,
         cor.args=list(ncase=rep(1000, k), ncntl=rep(1000,k)))

}