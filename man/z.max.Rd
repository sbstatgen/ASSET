% --- Source file: man/z.max.Rd ---
\name{z.max}
\alias{z.max}
\title{Z-score Maximization}
\description{
Function to maximize z-scores over subsets of traits or subtypes, with possible restrictions and weights. 
Should not be called directly. See details. 
}

\usage{
z.max(k, snp.vars, side, meta.def, meta.args, th=NULL, z.sub=rep(1, k), 
      sub.def=NULL, sub.args=NULL)
}
\arguments{
  \item{k}{Single integer. The total number of traits or studies or subtypes being analyzed. No default}
  \item{snp.vars}{Vector of integers or string labels for the SNPs being analyzed. No default.}  
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximzed), \code{side} should be 2. For one-tailed tests, 
  \code{side} should be 1 (positive tail is assumed). Default is 2, ignored when search is 2.}
  \item{meta.def}{Function that calculates z-scores for a given subset, for all the SNPs. Should accept a subset (logical vector of length k) as its
  first argument, followed by a list of SNPs (subset of snp.vars) as its second argument. 
  Should return a named list with at least the name "z", which is the vector of z-scores. The length of the vector should be the same
    length as \code{snp.vars}. 
  Missing z-scores if any are treated as zero, in the maximization. No default.}
  \item{meta.args}{Other arguments to be passed to \code{meta.def} as a named list. For example, this could include an entire data frame containing individual level data as in case of     
  subtype analysis, or sample sizes and correlation matrix in case of meta-analysis of heterogeneous traits. No default.}
  \item{th}{A vector of thresholds for each SNP, beyond which to stop maximization for that SNP. Default is a threshold of -1 for each SNP
  , implying no threshold. This argument is for internal use.}
  \item{sub.def}{A function to restrict subsets, e.g., order restrictions in subtype analysis. Should accept a subset (a logical vector of size k) as
  its first argument and should return TRUE if the subset satisfies restrictions and FALSE otherwise. Default is NULL implying all (2^k - 1) subsets
  are considered in the maximum.}
  \item{sub.args}{Other arguments to be passed to \code{sub.def} as list. Default is NULL (i.e. none).}
  %\item{wt.def}{A function that gives weight of one subset with respect to another. Should accept two subsets as first two argumets and return
  %a single positive weight. Default NULL. Currently this option is not implemented and the argument is ignored.}
  %\item{wt.args}{Other arguments to be passed to \code{wt.def} as a list. Default NULL. Currently ignored.}
  \item{z.sub}{Subset of traits/subtypes over whose subsets the maximization should be restricted. Default is all traits/subtypes. (i.e. none).}
}
\value{
   A list with two components. A vector of optimized z-scores (opt.z) and a logical matrix (opt.s) of dimension length(snp.vars) by \code{k}.
   Each row of (opt.s) has indicators of each trait/subtype being included in the best (optimal) subset.
}
\details{
This function loops through all possible (2^k - 1) subsets of (k) studies (or traits or subtypes), skips subsets that are not valid (e.g. 
that do not satisfy order restrictions), and maximizes the z-scores or re-weighted z-scores if weights are specified.  
The function is vectorized to handle blocks of SNPs at a time. 
This is a helper function that is called internally by \code{\link{h.traits}} and \code{\link{h.types}} 
and should not be called directly. The arguments of this function that have defaults, can be customized using
the argument \code{zmax.args} in \code{\link{h.traits}} and \code{\link{h.types}}. Specifying a a subset of traits/subtypes \code{z.sub} to be considered for maximization, helps
in speeding up the code when the number of traits or subtypes is relatively large. For example if \code{p.bound=0.25} is chosen in \code{h.traits}, on an average (under the null)
only a quarter of the traits will be maximized, allowing more traits to be analyzed in a computationally feasible manner. Note that the studies being 
maximized over will vary from SNP to SNP, and appropriate multiple-testing adjustment is done internally to account for this pre-selection.
}

\examples{
  set.seed(123)

  # Define the function to calculate the z-scores
  meta.def <- function(logicalVec, SNP.list, arg.beta, arg.sigma) {

    # Get the snps and subset to use
    beta <- as.matrix(arg.beta[SNP.list, logicalVec])
    se   <- as.matrix(arg.sigma[SNP.list, logicalVec])
    test <- (beta/se)^2
    ret  <- apply(test, 1, max)
    list(z=ret) 
  }

  # Define the function to determine which subsets to consider
  sub.def <- function(logicalVec, args) {
    # Only allow the cummulative subsets:
    # TRUE FALSE FALSE FALSE ...
    # TRUE TRUE FALSE FALSE ...
    # TRUE TRUE TRUE FALSE ...
    # etc
    sum <- sum(logicalVec)  
    ret <- all(logicalVec[1:sum])
    ret
  }

  # Assume there are 10 subtypes and 3 SNPs
  k        <- 10
  snp.vars <- 1:3
  
  # Generate some data 
  nsnp     <- length(snp.vars)
  beta     <- matrix(-0.5 + runif(k*nsnp), nrow=nsnp)
  sigma    <- matrix(runif(k*nsnp)^2, nrow=nsnp) 

  meta.args <- list(arg.beta=beta, arg.sigma=sigma)
  
  z.max(k, snp.vars, 2, meta.def, meta.args, sub.def=sub.def)
}