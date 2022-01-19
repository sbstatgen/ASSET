# Nov 01 2011  Add function h.forest


myInv <- function(A){

	Ainv <- try(solve(A))
	if (inherits(Ainv, "try-error")) Ainv <- ginv(A)
	Ainv
}

mySolve <- function(A, B){

	AinvB <- try(solve(A, B))
	if (inherits(AinvB, "try-error")) AinvB <- ginv(A) %*% B 
	AinvB
}

corr.mat.logit <- function(nmat11, nmat00, nmat10=NULL, nmat01=NULL){

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
findMiss.vars <- function(x, vars=NULL, miss=NA) {

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
			, sub.def = NULL, sub.args = NULL, wt.def=NULL, wt.args=NULL) {
			
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

getCI <- function(beta, sd, level=0.05) {

	half.width <- qnorm(1 - level/2) * sd
	low <- beta - half.width
	high <- beta + half.width
	return(data.frame(mid = exp(beta), low=exp(low), high=exp(high)))
}

h.summary <- function(rlist, level = 0.05, digits = 3){

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


h.forest <- function(k, snp.var, t.lab, rlist, res, side, level, p.adj, digits){

	k <- length(t.lab)
	ind <- c(NA, NA, NA)
	
	res0 <- rlist[[1]]
	res1 <- rlist[[2]]
	res2 <- rlist[[3]]

	pos <- (res$beta > 0 )
	neg <- (res$beta < 0 )
	null <- rep(FALSE, k)
	for(i in seq_len(3))
	{
		if(!is.null(rlist[[i]])) ind[i] <- pmatch(snp.var, names((rlist[[i]])$pval))
	}

	if(!is.null(res2))
	{
		if(side == 1)
		{
			pos <- (pos & res2$pheno[ind[3], ])
			neg <- (neg & res2$pheno[ind[3], ])
			null <- (null | !res2$pheno[ind[3], ])
		} else
		{
			pos <- (pos & res2$pheno.1[ind[3], ])
			neg <- (neg & res2$pheno.2[ind[3], ])
			null <- (null | (!res2$pheno.1[ind[3], ] & !res2$pheno.2[ind[3], ]))		
		}
	}
	
	CI <- getCI(res$beta, res$sd, level = level)

	CI.str <- apply(data.matrix(CI), 1, function(x)
					{
					paste("(", round(x[2], digits=digits), ", ", 
	                               round(x[3], digits=digits), ")", sep="")
					})
	
	if(p.adj) res$pval <- pmin(res$pval * k, 1)
	pval.str <- format(res$pval, digits=2, scientific=TRUE)
	pval.str[is.na(res$pval)] <- NA_character_
	
	
					
	nrow <- 1 
	nrow <- nrow + (if(any(neg)) (sum(neg) + 2) else 0) 
	nrow <- nrow + (if(any(null)) (sum(null) + 2) else 0)
	nrow <- nrow + (if(any(pos)) (sum(pos) + 2) else 0)
	nrow <- nrow + 1
	if(!is.null(res0)) nrow <- nrow + 1
	if(!is.null(res1)) nrow <- nrow + 1
	if(!is.null(res2) && side==1) nrow <- nrow + 1
	if(!is.null(res2) && side==2) nrow <- nrow + 3 
	
	
	mdat <- matrix(NA, nrow, 4)
	colnames(mdat) <- c("Mid", "Low", "High", "Pvalue")
	tab <- matrix(as.character(NA), nrow, 7)
	is.summary <- rep(FALSE, nrow)
	is.summary[1] <- TRUE
	colnames(tab) <- c("Phenotype", "OR", "GAP1", "CIstr", "GAP2", "PVALstr", "GAP3")
	tab[, "GAP1"] <- "    "
	tab[, "GAP2"] <- "    "
	tab[, "GAP3"] <- "  "
	
	dir <- list(neg, null, pos)
	names(dir) <- c("Negative", "Null", "Positive")
	cur <- 1
	for(i in seq_len(3))
	{
		clen <- 0
		if(any(dir[[i]]))
		{
			len <- sum(dir[[i]])
			v <- seq_len(len)
			mdat[cur + clen + 2 + v,"Mid"] <- log(CI$mid[dir[[i]]])
			mdat[cur + clen + 2 + v,"Low"] <- log(CI$low[dir[[i]]])
			mdat[cur + clen + 2 + v,"High"] <- log(CI$high[dir[[i]]])
			mdat[cur + clen + 2 + v,"Pvalue"] <- res$pval[dir[[i]]]

			tab[cur + clen + 2,"Phenotype"] <- names(dir)[i]
			tab[cur + clen + 2 + v,"Phenotype"] <- t.lab[dir[[i]]]
			tab[cur + clen + 2 + v,"OR"] <- round(CI$mid[dir[[i]]], digits=digits)
			tab[cur + clen + 2 + v,"CIstr"] <- CI.str[dir[[i]]]
			tab[cur + clen + 2 + v,"PVALstr"] <- pval.str[dir[[i]]]

			is.summary[cur + clen + 2] <- TRUE

			clen <- clen + len + 2
		}
		cur <- cur + clen
	}
	cur <- cur + 1
	
	for(i in seq_len(3))
	{
		res.i <- rlist[[i]]
		if(!is.null(res.i))
		{
			if(i != 3 || side !=2)
			{
				CI.i <- getCI(res.i$beta[snp.var], res.i$sd[snp.var], level = level)
				CI.i.str <- paste("(", round(CI.i[1, 2], digits=digits), ", "
					, round(CI.i[1, 3], digits=digits), ")", sep="")
				pval.i <- res.i$pval[ind[i]]
				pval.i.str <- if(is.na(pval.i)) NA_character_ else format(pval.i, digits=2, scientific=TRUE)

				mdat[cur + 1, "Mid"] <- log(CI.i$mid[1])
				mdat[cur + 1, "Low"] <- log(CI.i$low[1])
				mdat[cur + 1, "High"] <- log(CI.i$high[1])
				mdat[cur + 1, "Pvalue"] <- pval.i

				tab[cur + 1,"Phenotype"] <- names(rlist)[i]
				tab[cur + 1, "OR"] <- round(CI.i$mid[1], digits=digits)
				tab[cur + 1,"CIstr"] <- CI.i.str
				tab[cur + 1,"PVALstr"] <- pval.i.str

				is.summary[cur + 1] <- TRUE

				cur <- cur + 1
			} else
			{
				CI2.1 <- getCI(res.i$beta.1[snp.var], res.i$sd.1[snp.var], level = level)
				CI2.2 <- getCI(res.i$beta.2[snp.var], res.i$sd.2[snp.var], level = level)
				CI2.1.str <- paste("(", round(CI2.1[1, 2], digits=digits), ", ", 
			                        round(CI2.1[1, 3], digits=digits), ")", sep="")
				CI2.2.str <- paste("(", round(CI2.2[1, 2], digits=digits), ", ", 
		                        round(CI2.2[1, 3], digits=digits), ")", sep="")
				pval.i <- res.i$pval[ind[i]]
				pval.i.str <- if(is.na(pval.i)) NA_character_ else format(pval.i, digits=2, scientific=TRUE)
				pval.i.1 <- res.i$pval.1[ind[i]]
				pval.i.2 <- res.i$pval.2[ind[i]]
				pval.i.1.str <- if(is.na(pval.i.1)) NA_character_ else format(pval.i.1, digits=2, scientific=TRUE)
				pval.i.2.str <- if(is.na(pval.i.2)) NA_character_ else format(pval.i.2, digits=2, scientific=TRUE)

				v <- seq_len(3)
				mdat[cur + v, "Mid"] <- log(c(NA, CI2.1$mid[1], CI2.2$mid[1]))
				mdat[cur + v, "Low"] <- log(c(NA, CI2.1$low[1], CI2.2$low[1]))
				mdat[cur + v, "High"] <- log(c(NA, CI2.1$high[1], CI2.2$high[1]))
				mdat[cur + v, "Pvalue"] <- c(res.i$pval[ind[i]], res.i$pval.1[ind[i]], res.i$pval.2[ind[i]])

				tab[cur + v,"Phenotype"] <- c(names(rlist)[3], "    Positive", "    Negative")
				tab[cur + v, "OR"] <- round(c(NA, CI2.1$mid[1], CI2.2$mid[1]), digits=digits)
				tab[cur + v,"CIstr"] <- c(NA, CI2.1.str, CI2.1.str)
				tab[cur + v,"PVALstr"] <- c(pval.i.str, pval.i.1.str, pval.i.2.str)
				
				is.summary[cur + v] <- c(TRUE, TRUE, TRUE)
				
				cur <- cur + 3
			}
		}
	}
	hdr <- c("Phenotype", "OR", "", paste(round(100*(1 - level), digits=1), "% CI", sep="")
			, "", paste("P-value", if(p.adj) "(Adj)" else "(Unadj)"), "")
	tab[1,] <- hdr


#	par(omi = c(1, 0.05, 0.5, 0.05))

	rmeta::forestplot(tab, mdat[, "Mid"], mdat[, "Low"], mdat[, "High"], zero=0, is.summary=is.summary
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
h.forestPlot <- function(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2, calc.miss=FALSE) {

  which <- rlist[["which", exact=TRUE]]
  if (which == "h.types") {
    types.forest(rlist, snp.var, level=level, p.adj=p.adj, digits=digits, calc.miss=calc.miss)
  } else if (which == "h.traits") {
    traits.forest(rlist, snp.var, level=level, p.adj=p.adj, digits=digits, calc.miss=calc.miss)
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
 

