% --- Source file: man/h.traits.Rd ---
\name{h.traits}
\alias{h.traits}
\title{Heterogeneous traits or studies}
\description{
Performs one-sided or two-sided subset-based meta-analysis of  heteregeneous studies/traits.
}
\usage{
h.traits(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl, 
  cor=NULL, cor.numr=FALSE, search=NULL, side=2, meta=FALSE, 
  zmax.args=NULL, meth.pval=c("DLM","IS","B"), pval.args=NULL)
}
\arguments{
  \item{snp.vars}{A character vector giving the SNP names to be analyzed. No default.}
  \item{traits.lab}{A character vector giving the names/identifiers of the \code{k} studies/traits being analyzed.
        The order of this vector must match the columns of \code{beta.hat} and \code{sigma.hat} No default.}
  \item{beta.hat}{A matrix of dimension length(\code{snp.vars}) by (\code{k}) (or a vector of length \code{k} when 1 SNP is passed). Each row gives the coefficients obtained
  from the analysis of that SNP across the \code{k} studies/traits. No default.}
  \item{sigma.hat}{A vector or matrix of same dimension as beta.hat, giving the corresponding standard errors. No default.}
  \item{ncase}{The number of cases in each of the \code{k} studies. This can be same for each SNP, in which case \code{ncase} is a vector of length \code{k}.
  Alternatively if the number of non-missing cases analyzed for each SNP is known, \code{ncase} can be a length(\code{snp.vars}) by (\code{k}) matrix. No default.}
  \item{ncntl}{Same as \code{ncase} (above) for controls. No default.}
  \item{cor}{Either a \code{k} by \code{k} matrix of inter-study correlations or a list containing three case/control overlap matrices (for
case-control studies)
  named \code{N11}, \code{N00} and \code{N10}. See details. Default is NULL, so that studies are assumed to be independent.
The rows and columns of all matrices  needs to be in the same order as \code{traits.lab} and columns 
of \code{beta.hat} and \code{sigma.hat} }
  \item{cor.numr}{Logical. When specified as TRUE the correlation information is used for optimal weighting of the studies in the definition of the meta-analysis test-statistics. 
The default is FALSE. In either case, the correlation is accounted for variance calculation of the meta-analysis test-statistic.}
  \item{search}{1, 2 or NULL. Search option 1 and 2 indicate one-sided and two-sided subset-searches respectively. The default option is NULL that
automatically returns both one-sided and two-sided subset searches. 
  Default is NULL.}
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximized), \code{side} should be 2. For one-tailed tests, 
  \code{side} should be 1 (positive tail is assumed). Default is 2. The option is ignored when search=2 since the two-sided subset search is automatically a two-sided
test.}
  \item{meta}{Logical. When specified as TRUE, standard fixed effect meta-analysis results are returned together with results from subset-based meta-analysis. 
The Default is FALSE.}
  \item{zmax.args}{Optional arguments to be passed to \code{\link{z.max}} as a named list. This option can be useful if the user wants to
restrict subset searches in some structured way, for example, incorporating some ordering constraints.}
 \item{pval.args}{Optional arguments  to be
  passed to \code{\link{p.dlm}} or \code{\link{p.tube}} as a named list. This option can be useful if the user wants to
restrict subset searches in some structured way, for example, incorporating some ordering constraints.}
 \item{meth.pval}{The method of p-value computation. Currently the options are DLM (Discrete Local
  Maximum), IS(Exact Importance Sampling) and B (Bonferroni) with the default option being DLM. The IS method is currently computationally
feasible for analysis of at most k=10 studies/traits}
}
\value{
   A list containing 3 main component lists named:
    
   (1) "Meta" (Results from standard fixed effect meta-analysis of all  studies/traits). 
	   This list is non-null when \code{meta} is TRUE and contains 3 vectors named (pval, beta, sd) of length same as snp.vars.
   
   (2) "Subset.1sided" (one-sided subset search):
	   This list is non-null when \code{search} is NULL or 1 and contains, 4 vectors named (pval, beta, sd, sd.meta) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each phenotype. For a particular SNP and phenotype, the entry
	   has "TRUE" if this phenotype was in the selected subset for that SNP. In the output, the p-value is automatically adjusted for multiple testing 
           due to subset search. The beta  and sd correspond to the standard fixed-effect meta-analysis estimate and corresponding standard error estimate 
for the regression coefficient of a SNP based only on those studies/traits that are included in the identified subset.
The vector sd.meta gives the meta-analysis standard errors for estimates of beta based on studies in the identified subset ignoring the randomness of the subset.


	   
   (3) Subset.2sided (two-sided subset search)
       This list is non-null when \code{search} is NULL or 2 and contains 9 vectors named 
  (pval, pval.1, pval.2, beta.1, sd.1, beta.2, sd.2, sd.1.meta, sd.2.meta) of length
	   same as snp.vars and two matrices named "pheno.1" and "pheno.2" giving logical indicators of a phenotype being among the positively or negatively
	   associated subsets (respectively) as identified by 2-sided subset search. In the output, while pval provides the significance of the overall test-statistics that
combined association signals from two directions, pval.1 and and pval.2 return the corresponding level of significance for each of the component one-sided
test-statistics in the positive and negative directions.  The values (beta.1, sd.1, beta.2, sd.2) denote the corresponding meta-analysis estimate of regression coefficients
and standard errors for the identifed subsets of traits/studies that show association in positive and negative directions, respectively.
The vector sd.1.meta and sd.2.meta give the meta-analysis standard errors for estimates of beta based on studies in the identified subset ignoring the randomness of the subset in positive and negative directions, respectively.

The other objects in the list are the input arguments passed into \code{h.traits}.
}
\details{
The one-sided subset search maximizes the standard fixed-effect meta-analysis test-statistics over all possible subsets (or over a restricted set of 
subsets if such an option is specified) of studies to detect the best possible association signals. The p-value returned for the maximum test-statistics automatically
accounts for multiple testing penalty due to subset search and can be taken as an evidence of an overall association for the SNP accross the \code{k} studies/traits. 
The one-sided method automatically guarantees identification of studies/traits that have associations in the same direction and thus is useful in applications where
it is desirable to identify SNPs that shows effects in the same direction across multiple traits/studies. The two-sided subset search, applies one-side subset search separately for positively and negatively associated traits for a given SNP and then combines the association
signals from two directions into a single combined chi-square type statistic. The method is sensitive in detecting SNPs that may be associated with different
traits in different directions.

The methods allow for accounting for correlation among studies/subject that might arise due to shared subjects across distinct studies or due to correlation
among related traits in the same study. For application of the method for meta-analysis of case-control studies, the matrices \code{N11}, \code{N10}
 and \code{N00} denote the number subjects that are shared between studies by case-control status. By defintion, the
diagonals of the matrices \code{N11} and \code{N00} contain the number of cases and controls, respectively, in the \code{k} studies. Also, by definition, 
the diagonal of \code{N10} is zero since cases cannot serve as controls and vice versa in the same study. The most common situation
may involve shared controls accross studies, ie non-zero off-diagonal elements of the matrix N00. 

The output standard errors are approximate (based on inverting p-values) and are used for constructing confidence
intervals in \code{\link{h.summary}} and \code{\link{h.forestPlot}}.
}

\references{Samsiddhi Bhattacharjee, Preetha Rajaraman, Kevin B. Jacobs, William A. Wheeler, Beatrice S. Melin, Patricia  Hartge, GliomaScan Consortium, Meredith Yeager, Charles C. Chung, Stephen J. Chanock, Nilanjan Chatterjee. A subset-based  approach improves power and interpretation for combined-analysis of genetic association studies of heterogeneous traits. Am  J Hum Genet, 2012, 90(5):821-35 \cr
}

\seealso{ \code{\link{h.summary}}, \code{\link{h.forestPlot}} }

\examples{
 # Use the example data
 data(ex_trait, package="ASSET")

 # Display the data, and case/control overlap matrices
 data
 N00
 N11
 N10
 
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
 

\keyword{ traits } 
