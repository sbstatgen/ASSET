# Nov 01 2011  Add function h.forest


myInv <- function(A)
{
	Ainv <- try(solve(A))
	if (inherits(Ainv, "try-error")) Ainv <- ginv(A)
	Ainv
}

mySolve <- function(A, B)
{
	AinvB <- try(solve(A, B))
	if (inherits(AinvB, "try-error")) AinvB <- ginv(A) %*% B 
	AinvB
}

corr.mat.logit <- function(nmat11, nmat00, nmat10=NULL, nmat01=NULL)
{
	k <- nrow(nmat11)
	if(is.null(nmat10)) nmat10 <- 0
	if(is.null(nmat01)) nmat01 <- 0
	
	nvec0 <- diag(nmat00)
	nvec1 <- diag(nmat11)
	nvec <- nvec0 + nvec1
	
	neff <- (nvec0 * nvec1)/nvec
	rneff <- sqrt(neff)
	
	rmat <- (nmat11 * outer(1/nvec1, 1/nvec1) + nmat00 * outer(1/nvec0, 1/nvec0)
			 -nmat10 * outer(1/nvec1, 1/nvec0) - nmat01 * outer(1/nvec0, 1/nvec1)) * sqrt(outer(neff, neff))
	rmat
}

# Function to identify rows with missing values of certain variables. 
findMiss.vars <- function(x, vars=NULL, miss=NA) 
{	
	if(is.null(dim(x))) dim(x) <- c(length(x), 1)
	if (is.null(vars)) vars <- seq_len(ncol(x))
	
  ##SMB: changed this a bit to speed it up
	temp <- rep(TRUE, times=nrow(x)) 
  if(!is.na(miss)) x[x == miss] <- NA
	for (var in vars) temp <- temp & !(is.na(x[, var]))
	##for (var in vars) temp <- temp & !(x[, var] %in% miss)
	
	temp
	
} # END: findMiss.vars


z.max <- function(k, snp.vars, side, meta.def, meta.args, th = rep(-1, length(snp.vars))
			, sub.def = NULL, sub.args = NULL, wt.def=NULL, wt.args=NULL)
{
	nsnp <- length(snp.vars)
	snp.sub <- seq_len(nsnp)
	
	if(length(th) == 1) th = rep(th, nsnp)
	if(length(th) != nsnp) stop("Expected one threshold or a threshold for each SNP")
	
	x <- rep(0, k)
	opt.s <- matrix(FALSE, nsnp, k)
	opt.s[, k] <- 1
	
	opt.z <- rep(0, nsnp)
	nopt <- rep(0, nsnp)
	ztmp <- opt.z
	
	i <- 1
	ipos <- k
    subsetCount <- 0
	while(i <= ((2^k) - 1))
	{	
		ipos <- max(which(x == 0))
		x[ipos:k] <- 0
		x[ipos] <- 1
		set <- as.logical(x)

		if(!is.null(sub.def)) if(! do.call(sub.def, c(list(set), sub.args))){ i <- i + 1 ; next }
	    subsetCount <- subsetCount + 1
#		sub <- which(set)
		nsub <- sum(set)

		ztmp <- do.call(meta.def, c(list(as.logical(set), snp.vars[snp.sub]), meta.args))$z
#		print(c(i, set, ztmp, opt.z))
	
		ztmp[is.na(ztmp) | is.nan(ztmp)] <- 0
		
		lrr <- 0
		
		#if(side == 2) pos <- (pnorm(abs(ztmp), 0, 1, lower=FALSE, log=TRUE) - pnorm(abs(opt.z[snp.sub]), 0, 1, lower=FALSE, log=TRUE) < lrr)
		#else pos <- (pnorm(ztmp, 0, 1, lower=FALSE, log=TRUE) - pnorm(opt.z[snp.sub], 0, 1, lower=FALSE, log=TRUE) < lrr)
        if(side == 2) {
          pos <- pnorm(abs(ztmp), 0, 1, lower.tail=FALSE, log.p=TRUE) <
                 pnorm(abs(opt.z[snp.sub]), 0, 1, lower.tail=FALSE, log.p=TRUE)
        } else {
          pos <- pnorm(ztmp, 0, 1, lower.tail=FALSE, log.p=TRUE) <
                 pnorm(opt.z[snp.sub], 0, 1, lower.tail=FALSE, log.p=TRUE)
        }

		npos <- sum(pos)

	
		if(npos > 0)
		{
            jj          <- snp.sub[pos]   
			opt.z[jj]   <- ztmp[pos]
			opt.s[jj, ] <- matrix(set, npos, k, byrow=TRUE)
			nopt[jj]    <- nsub

		}
		
		if(side == 2) ltmp <- (th[snp.sub] < 0 | (abs(opt.z[snp.sub]) < th[snp.sub]))
		else ltmp <- (th[snp.sub] < 0 | (opt.z[snp.sub] < th[snp.sub]))
	
		snp.sub <- snp.sub[ltmp]
		if(length(snp.sub) < 1) break
		
		i <- i + 1
	}
#	print(c(i, set, ztmp, opt.z))
#	stop()

    # Change here
	#opt.z[opt.z == 0] <- NA
    
    opt.z[!is.finite(opt.z)] <- 0
    if (all(opt.z == 0)) opt.s[] <- 0

	ret <- list(opt.z = opt.z, opt.s = opt.s, subsetCount=subsetCount)

	ret
	
}

getCI <- function(beta, sd, level=0.05)
{
	half.width <- qnorm(1 - level/2) * sd
	low <- beta - half.width
	high <- beta + half.width
	return(data.frame(mid = exp(beta), low=exp(low), high=exp(high)))
}

h.summary <- function(rlist, level = 0.05, digits = 3)
{
	
	frlist <- list()
	nmlist <- NULL
	
	rnames <- names(rlist)
	if(is.null(rnames)) rnames <- paste("Results-", seq_along(rlist), sep="")

	for(i in seq_along(rlist))
	{
		if(!is.null(rlist[[i]]))
		{
			res <- rlist[[i]]
			#if(!is.list(res)) stop("Expected a list")
            if (!is.list(res)) next

			len <- length(res$pval)
			#if(len == 0) { warning("Empty p-value vector") ; next }
            if(len == 0) next 

			svars <- names(res$pval)
			if(is.null(svars)) svars <- paste("SNP-", seq_len(len), sep="")
			spheno <- NULL
			
			if(!is.null(res$beta))
			{
				if(length(res$beta) != len || length(res$sd) != len) stop("beta and sd expected")
				CI <- getCI(res$beta, res$sd, level = level)
				#mat <- signif(cbind(res$pval, CI$mid, CI$low, CI$high), digits=digits)
                mid11 <- round(CI$mid, digits=digits)
                low11 <- round(CI$low, digits=digits)
                high11 <- round(CI$high, digits=digits)
                #mat <- cbind(res$pval, CI$mid, CI$low, CI$high)
                mat <- cbind(res$pval, mid11, low11, high11)
				colnames(mat) <- c("Pvalue", "OR", "CI.low", "CI.high")
			}

			if(!is.null(res$beta.1) && !is.null(res$beta.2))
			{
				if(length(res$pval.1) != len || length(res$beta.1) != len || length(res$sd.1) != len) stop("pval.1, beta.1 and sd.1 expected")
				CI.1 <- getCI(res$beta.1, res$sd.1, level = level)
				if(length(res$pval.2) != len || length(res$beta.2) != len || length(res$sd.2) != len) stop("pval.2, beta.2 and sd.2 expected")
				CI.2 <- getCI(res$beta.2, res$sd.2, level = level)

				#mat <- signif(cbind(res$pval, res$pval.1, res$pval.2, CI.1$mid, CI.1$low, CI.1$high
				#					, CI.2$mid, CI.2$low, CI.2$high), digits=digits)
                mat <- cbind(res$pval, res$pval.1, res$pval.2, CI.1$mid, CI.1$low, CI.1$high
									, CI.2$mid, CI.2$low, CI.2$high)
				colnames(mat) <- c("Pvalue", "Pvalue.1", "Pvalue.2", "OR.1", "CI.low.1", "CI.high.1", "OR.2", "CI.low.2", "CI.high.2")
                for (ii in 4:9) mat[, ii] <- round(mat[, ii], digits=digits)
			}

			if(!is.null(res$pheno))
			{
				if(nrow(res$pheno) != len) stop("pheno expected for all SNPs")
				types.lab <- colnames(res$pheno)
				if(is.null(types.lab)) types.lab <- seq_len(ncol(res$pheno))
				spheno <- matrix(apply(res$pheno, 1, function(x) paste(types.lab[x], collapse=",")), ncol = 1)
				colnames(spheno) <- "Pheno"
			}

			if(!is.null(res$pheno.1) && !is.null(res$pheno.2))
			{
				if(nrow(res$pheno.1) != len || nrow(res$pheno.2) != len) stop("pheno.1 and pheno.2 expected for all SNPs")
				types.lab <- colnames(res$pheno.1)
				if(is.null(types.lab)) types.lab <- seq_len(ncol(res$pheno.1))
				spheno.1 <- apply(res$pheno.1, 1, function(x) paste(types.lab[x], collapse=","))
				spheno.2 <- apply(res$pheno.2, 1, function(x) paste(types.lab[x], collapse=","))
				spheno <- cbind(spheno.1, spheno.2)
				colnames(spheno) <- c("Pheno.1", "Pheno.2")
			}
			
			if(is.null(spheno)) fres <- cbind(data.frame(SNP = svars), data.frame(mat))
			else fres <- cbind(data.frame(SNP = svars), data.frame(mat), data.frame(spheno))
			
			rownames(fres) <- NULL
			
			frlist <- c(frlist, list(fres))
			nmlist <- c(nmlist, rnames[i])
		}
	}
	names(frlist) <- nmlist
	frlist
}


h.forest <- function(k, snp.var, t.lab, rlist, res, side, level, p.adj, digits)
{
	k <- length(t.lab)
	
	res0 <- rlist[[1]]
	res1 <- rlist[[2]]
	res2 <- rlist[[3]]
	
	i0 <- pmatch(snp.var, names(res0$pval))
	i1 <- pmatch(snp.var, names(res1$pval))
	i2 <- pmatch(snp.var, names(res2$pval))
	
	CI0 <- getCI(res0$beta[snp.var], res0$sd[snp.var], level = level)
	CI1 <- getCI(res1$beta[snp.var], res1$sd[snp.var], level = level)	
	if(side == 1)
	{
		CI2 <- getCI(res2$beta[snp.var], res2$sd[snp.var], level = level)
	} else
	{
		CI2.1 <- getCI(res2$beta.1[snp.var], res2$sd.1[snp.var], level = level)
		CI2.2 <- getCI(res2$beta.2[snp.var], res2$sd.2[snp.var], level = level)
	}

	CI0.str <- paste("(", round(CI0[1, 2], digits=digits), ", ", 
                          round(CI0[1, 3], digits=digits), ")", sep="")
	CI1.str <- paste("(", round(CI1[1, 2], digits=digits), ", ", 
                          round(CI1[1, 3], digits=digits), ")", sep="")
	if(side == 1)
	{
		CI2.str <- paste("(", round(CI2[1, 2], digits=digits), ", ", 
                              round(CI2[1, 3], digits=digits), ")", sep="")
	} else
	{
		CI2.1.str <- paste("(", round(CI2.1[1, 2], digits=digits), ", ", 
                                round(CI2.1[1, 3], digits=digits), ")", sep="")
		CI2.2.str <- paste("(", round(CI2.2[1, 2], digits=digits), ", ", 
                                round(CI2.2[1, 3], digits=digits), ")", sep="")
	}
	
	CI <- getCI(res$beta, res$sd, level = level)
	CI.str <- apply(data.matrix(CI), 1, function(x)
					{
					paste("(", round(x[2], digits=digits), ", ", 
                               round(x[3], digits=digits), ")", sep="")
					})
	if(side == 1)
	{
		pos <- (res$beta > 0 & res2$pheno[i2, ])
		neg <- (res$beta < 0 & res2$pheno[i2, ])
		null <- (!res2$pheno[i2, ])
	}
	else
	{
		pos <- (res$beta > 0 & res2$pheno.1[i2, ])
		neg <- (res$beta < 0 & res2$pheno.2[i2, ])
		null <- (!res2$pheno.1[i2, ] & !res2$pheno.2[i2, ])
	}

	mid <- c(NA, NA, CI$mid[neg], NA, NA, CI$mid[null]
			 , NA, NA, CI$mid[pos], NA, CI0$mid[1], CI1$mid[1])

	if(side == 1) mid <- c(mid, CI2$mid[1])
	else mid <- c(mid, NA, CI2.1$mid[1], CI2.2$mid[1])

	low <- c(NA, NA, CI$low[neg], NA, NA, CI$low[null]
			 , NA, NA, CI$low[pos], NA, CI0$low[1], CI1$low[1])
	
	if(side == 1) low <- c(low, CI2$low[1])
	else low <- c(low, NA, CI2.1$low[1], CI2.2$low[1])
	
	high <- c(NA, NA, CI$high[neg], NA, NA, CI$high[null]
			 , NA, NA, CI$high[pos], NA, CI0$high[1], CI1$high[1])
	
	if(side == 1) high <- c(high, CI2$high[1])
	else high <- c(high, NA, CI2.1$high[1], CI2.2$high[1])

	
	if(p.adj) res$pval <- pmin(res$pval * k, 1)
	
	pvalues <- c(NA, NA, res$pval[neg], NA, NA, res$pval[null], NA, NA, res$pval[pos], NA, res0$pval[i0]
				 , res1$pval[i1], res2$pval[i2])
	
	if(side == 2) pvalues <- c(pvalues, res2$pval.1[i2], res2$pval.2[i2])
	
	pval.str <- format(pvalues, digits=2, scientific=TRUE)
	pval.str[is.na(pvalues)] <- NA_character_

	pheno.col <- c("Phenotype", NA, "Negative", t.lab[neg], NA, "Null", t.lab[null], NA, "Positive"
				   , t.lab[pos], NA, names(rlist))
	if(side == 2) pheno.col <- c(pheno.col, "    Positive", "    Negative")
	
	CI.vec <- c(NA, NA, CI.str[neg], NA
				, NA, CI.str[null], NA, NA, CI.str[pos], NA, CI0.str, CI1.str)

	if(side == 1) CI.vec <- c(CI.vec, CI2.str)
	else CI.vec <- c(CI.vec, c(NA, CI2.1.str, CI2.1.str))

	tabletext <- cbind(pheno.col, c("OR", round(mid, digits=digits))
					   , rep("", length(mid) + 1)
					   , c(paste(round(1 - level, digits=3), "% CI", sep=""), CI.vec)
					   , rep("", length(mid) + 1)
					   , c(paste("P-value", if(p.adj) "(Adj)" else "(Unadj)"), pval.str)
					   , rep("", length(mid) + 1)					   
					   )

#	par(omi = c(1, 0.05, 0.5, 0.05))
	is.summary <-  c(rep(TRUE, 3), rep(FALSE, sum(neg)), FALSE, TRUE, rep(FALSE, sum(null))
					 , FALSE, TRUE, rep(FALSE, sum(pos)), FALSE, rep(TRUE, 2))
	if(side == 1) is.summary <- c(is.summary, TRUE)
	if(side == 2) is.summary <- c(is.summary, c(TRUE, TRUE, TRUE))
	rmeta::forestplot(tabletext, c(NA, log(mid)), c(NA, log(low)), c(NA, log(high)), zero=0, is.summary=is.summary
				, clip=c(log(1/5),log(5)), xlog=TRUE, col=meta.colors(box="royalblue", lines="darkblue", summary="royalblue"))

	title(main=snp.var)
}

# Function to compute the ordinary meta-analysis standard error for estimate of beta based 
# on studies selected as part of s ignoring the randomness of the set s
meta.se <- function(rmat, sigma, subset=NULL) {
  # rmat   Correlation matrix for beta
  # sigma  Vector of standard errors
  
  dim(sigma) <- NULL
  if (!is.null(subset)) {
    sigma <- sigma[subset]
    rmat  <- as.matrix(rmat[subset, subset])
  }
  n          <- length(sigma)
  if (!n) return(NA)
  if (n == 1) {
    temp <- sigma
    dim(temp) <- c(1, 1)
  } else {
    temp      <- diag(sigma)
  }
  COV        <- temp %*% rmat %*% temp
  s2inv      <- 1/(sigma*sigma)
  sum        <- sum(s2inv)
  denom      <- sum*sum
  dim(s2inv) <- c(1, n)
  se <- sqrt((s2inv %*% COV %*% t(s2inv))/denom)
  dim(se) <- NULL
  se

} # END: meta.se

# Top-level function to make a forest plot
h.forestPlot <- function(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2) {

  which <- rlist[["which", exact=TRUE]]
  if (which == "h.types") {
    types.forest(rlist, snp.var, level=level, p.adj=p.adj, digits=digits)
  } else if (which == "h.traits") {
    traits.forest(rlist, snp.var, level=level, p.adj=p.adj, digits=digits)
  } else {
    stop("rlist is not a valid return list from h.types() or h.traits()")
  }

  NULL

} # END: h.forest

# Function to check that variables exist in the input data set
checkData.vars <- function(data, vars) {

  if (is.null(vars)) return(NULL)

  cnames <- colnames(data)
  cflag  <- is.character(vars)
  if (cflag) {
    temp <- !(vars %in% cnames)
    if (any(temp)) {
      miss <- vars[temp]
      print(miss)
      emsg <- paste("ERROR: The above variables were not found in the input data")
      stop(emsg) 
    } 
  } else {
    nc <- ncol(data)
    temp <- !(vars %in% seq_len(nc))
    if (any(temp)) {
      miss <- vars[temp]
      print(miss)
      emsg <- paste("ERROR: The above variable positions are not valid")
      stop(emsg) 
    } 
  } 

  NULL

} # END: checkData.vars

# Function to check that variables exist in the input data set
checkData.rownames <- function(data, vars) {

  if (is.null(vars)) return(NULL)
  if (is.null(dim(data))) return(NULL)

  cnames <- rownames(data)
  cflag  <- is.character(vars)
  c2flag <- !is.null(cnames)
  if (cflag && c2flag) {
    temp <- !(vars %in% cnames)
    if (any(temp)) {
      miss <- vars[temp]
      print(miss)
      emsg <- paste("ERROR: The above row names were not found in the input data")
      stop(emsg) 
    } 
  } else {
    nc <- nrow(data)
    temp <- !(vars %in% seq_len(nc))
    if (any(temp)) {
      miss <- vars[temp]
      print(miss)
      emsg <- paste("ERROR: The above row positions are not valid")
      stop(emsg) 
    } 
  } 

  NULL

} # END: checkData.vars

# Function to convert positions to variable names
convertPosToVars <- function(data, pos) {

  if (is.null(pos)) return(NULL)
  if (is.character(pos)) return(pos)

  ret <- colnames(data)[pos]
  ret

} # END: convertPosToVars 

