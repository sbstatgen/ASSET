% --- Source file: man/h.types.Rd ---
\name{h.types}
\alias{h.types}
\title{Heterogeneous Subtype analysis}
\description{
Subset-based analysis of case-control studies with heterogeneous disease subtypes.
}
\usage{
h.types(dat, response.var, snp.vars, adj.vars, types.lab, cntl.lab, 
        subset=NULL, method=NULL, side=2, logit=FALSE, test.type="Score", 
        zmax.args=NULL, pval.args=NULL, p.bound = 1, 
        NSAMP=5000, NSAMP0=50000)
}
\arguments{
  \item{dat}{A data frame containing individual level data for phenotype (disease status/subtype information), covariate data and
   SNPs. No default.}
  \item{response.var}{Variable name or position of the response variable column in the data frame.  
        This variable needs to contain disease status/subtype information in the data frame. No default.}
\item{snp.vars}{A character or numeric vector giving the variable names or positions of the SNP variables. Missing values for SNP genotypes are indicated by NA. No default.}
 \item{adj.vars}{A character or numeric vector containing the variable names or positions of the columns in the data frame that would be used as adjusting covariates in the analysis. 
   Use NULL if no covariates are used for adjustment.}
\item{types.lab}{NULL or a character vector giving the names/identifiers of the disease subtypes in \code{response.var} to be
                 included in the analysis. If NULL, then all subtypes will be included. No default.}
\item{cntl.lab}{A single character string giving the name/identifier of controls (disease-free subjects) in \code{response.var}. No default.}
  \item{subset}{A logical vector with length=nrow(\code{dat}) indicating the subset of rows of the data frame to be included in the analysis. Default is NULL, all rows are used.}
  \item{method}{A single character string indicating the choice of method as "case-control" or "case-complement".  
  The Default option is NULL which will carry out both types of analysis. For the case-complement analysis of disease subtype \code{i}, 
  the set of control subjects is formed by taking the complement of disease subtype \code{i}, ie the original controls and the cases
  not defined by disease subtype \code{i}.}
  \item{side}{A numeric value of either 1 or 2 indicating whether one or two-tailed p-values should be computed, respectively. The default is 2.}
  \item{logit}{If TRUE, results are returned from an overall case-control analysis using standard logistic regression. Default is FALSE.}
  \item{test.type}{A character string indicating the type of tests to be performed. The current options are "Score" and "Wald". The default is "Score."}
 \item{zmax.args}{Optional arguments to be passed to \code{\link{z.max}} as a named list. This option can be useful if the user wants to
   restrict subset searches in some structured way, for example, incorporating ordering constraints.}
%\item{meth.pval}{A character string indicating the method of evaluating the p-value. Currently the only option is "DLM" (Discrete Local
%  Maximum)
%"IS" (Exact Importance Sampling) and "B" (Bonferroni) with the default option being DLM. 
%The IS method is currently computationally feasible for analysis of at most k=10 studies/traits}
 \item{pval.args}{Optional arguments  to be
  passed to \code{\link{p.dlm}} as a named list. This option can be useful if the user wants to
restrict subset searches in some structured way, for example, incorporating ordering constraints.}
  \item{p.bound}{P-value threshold for screening studies based on marginal association before
  performing subset search. Default is 1, that is all studies are included in the subset search.
  The p-value for the overall procedure accounts for this pre-screening step. See details.}
 % \item{add.enrich}{TRUE or FALSE to retrieve some of the lost information using binomial enrichment. Default is FALSE.}
 % \item{cond.all}{TRUE or FALSE to indicate whether p-values are computed conditional on all z-scores being above or
 % below the screening cutoff during pre-screening.
 %   That is condition on Z1>c, Z2>c, Z3<c, Z4<c etc, or to only condition on those retained for subset search. 
 %   i.e condition on Z1>c and Z2>c and search among subsets of those defined by cond.all=FALSE. The default is TRUE.}
  \item{NSAMP}{Number of samples from a truncated multivariate normal distribution used to compute the DLM p-value.
              The default is 5000.}
  \item{NSAMP0}{Number of samples from truncated multivariate normal distribution used to calculate
  the probability of the truncation region in DLM p-value calculation. For 1-sided subset search
  this is ignored unless \code{p.bound} < 1. The default is 50000. See details.}
}
\value{
   A list containing 3 component lists named:
    
   (1) "Overall.Logistic" (output for overall case-control analysis using standard logistic regression): 
	   This list is non-null when \code{logit} is TRUE and contains 3 vectors named (pval, beta, sd) of length same as snp.vars.
   
   (2) "Subset.Case.Control" (output for subset-based case-control analysis):
	   This list is non-null when \code{method} is NULL or "case-control". The output contains, 3 vectors named (pval, beta, sd) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each disease subtype. For a particular SNP and disease-subtype, the
	   corresponding entry is "TRUE" if that disease subtype is included the best subset of disease subtypes that is identified to be associated with
the SNP in the subset-based case-control analysis. In the output, the p-value is automatically adjusted for multiple testing 
           due to subset search. The beta and sd corresponds to estimate of log-odds-ratio and standard error for a SNP from a logistic regression analysis involving
the cases of the identified disease subtypes and the controls.  

   (3) "Subset.Case.Complement" (output for subset-based case-complement analysis):
       This list is non-null when \code{method} is NULL or "case-complement". The output contains, 3 vectors named (pval, beta, sd) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each disease subtype. For a particular SNP and disease-subtype,  the
	   corresponding entry is "TRUE" if that disease subtype is included the best subset of disease subtypes that is identified to be associated with
the SNP in the subset-based case-complement analysis. In the output, the p-value is automatically adjusted for multiple testing 
           due to subset search. The beta and sd corresponds to estimate of log-odds-ratio and standard error for the SNP from a logistic regression analysis involving
the cases of the selected disease subtypes and the whole complement set of subjects that includes original controls and the cases of unselected disease subtypes. 
}
\details{
The output standard errors are approximate (based on inverting DLM pvalues) and are used for constructing confidence intervals 
in \code{\link{h.summary}} and \code{\link{h.forestPlot}}. For a particular SNP, if any of the genotypes are missing, then those
subjects will be removed from the analysis for that SNP.

Currently ASSET calculates p-values by a stochastic approximation to the DLM formula as described in
Bhattacharjee et al. (In Preparation). The method works by simulating truncated multivariate normal
variates by importance sampling to estimate the probability term appearing in the DLM formula.
Since version 2.0.0, the previous \code{meth.pval}="DLM" option to calculate upper bound 
p-values (as in Bhattacharjee et al. 2012) has been dropped as the current stochastic 
approximation is expected to be more accurate in all cases although slightly slower.
The new p-value method also enables pre-screening of traits by the \code{p.bound} argument.

Specifying a p-value upper bound through \code{p.bound}, helps in speeding up the code when the number of traits or subtypes is relatively large. For example if \code{p.bound=0.25} is 
chosen, on an average (under the null) only a quarter of the traits will be used for subset search, 
 allowing more traits to be analyzed in a computationally feasible manner. 
%Note that the studies being maximized over will vary from SNP to SNP, and appropriate multiple-testing adjustment is done internally
% to account for this pre-selection. 
%A caveat is that for \code{p.dlm} when thisoption is set to a non-default value, i.e. \code{p.bound<1} 
%and studies are not independent, final p-values will be stochastic. 
%This is because function \code{pmvnorm} in package \code{mvtnorm} based on a stochastic algorithm is used in estimating the
% probability of the p-value based pre-selection.
The arguments \code{NSAMP} and \code{NSAMP0} give the number of importance sampling replicates to be generated. 
Either of these can be increased to achieve more accuracy at the cost of computational speed or vice versa.
}

\seealso{ \code{\link{h.summary}}, \code{\link{h.forestPlot}} }

\examples{
 # Use the example data
 data(ex_types, package="ASSET")

 # Display the first 10 rows of the data and a table of the subtypes
 data[1:10, ]
 table(data[, "TYPE"])
 
 # Define the input arguments to h.types. 
 snps     <- paste("SNP_", 1:3, sep="")
 adj.vars <- c("CENTER_1", "CENTER_2", "CENTER_3")
 types <- paste("SUBTYPE_", 1:5, sep="")

 # SUBTYPE_0 will denote the controls
 res <- h.types(data, "TYPE", snps, adj.vars, types, "SUBTYPE_0", subset=NULL, 
        method="case-control", side=2, logit=FALSE, test.type="Score", 
        zmax.args=NULL, pval.args=NULL)

 
 h.summary(res)

}

\references{Samsiddhi Bhattacharjee, Preetha Rajaraman, Kevin B. Jacobs, William A. Wheeler, Beatrice S. Melin, Patricia  Hartge, GliomaScan Consortium, Meredith Yeager, Charles C. Chung, Stephen J. Chanock, Nilanjan Chatterjee. A subset-based  approach improves power and interpretation for combined-analysis of genetic association studies of heterogeneous traits. Am  J Hum Genet, 2012, 90(5):821-35 \cr
}
