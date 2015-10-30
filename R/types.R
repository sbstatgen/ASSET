# Nov 01 2011  Pass the subset option in h.types to types.wald
# Nov 03 2011  Allow types.lab to be NULL
# Apr 24 2014  Change sub.def1 function definition

h.types <- function(dat, response.var, snp.vars, adj.vars, types.lab, cntl.lab, subset=NULL, method=NULL, side=2
, logit=FALSE, test.type = "Score", zmax.args = NULL, pval.args = NULL,
  p.bound = 1, NSAMP=5000, NSAMP0=50000)

{
    meth.pval  <- "DLM"
    add.enrich <- FALSE
    cond.all   <- TRUE

    if (!is.data.frame(dat)) stop("ERROR: dat must be a data frame")
    if ((is.null(side)) || (!(side %in% c(-1, 1, 2)))) stop("ERROR: side must be -1, 1 or 2")
    if ((is.null(meth.pval)) || (!(meth.pval %in% c("DLM", "IS", "B")))) stop("ERROR: meth.pval must be DLM, IS, or B")
    if ((is.null(logit)) || (!(logit %in% 0:1))) stop("ERROR: logit must be TRUE or FALSE")
    if ((is.null(test.type)) || (!(test.type %in% c("Score", "Wald")))) stop("ERROR: test.type must be Score or Wald")
    if (!is.null(pval.args)) {
      if (!is.list(pval.args)) stop("ERROR: pval.args must be a named list")
      temp <- names(pval.args)
      if (is.null(temp)) stop("ERROR: pval.args must be a named list")
    }

	if(p.bound>1 || p.bound <0) stop("Expected a p-value bound between 0 and 1")

    # Check that the variables exist
    checkData.vars(dat, response.var)
    checkData.vars(dat, snp.vars)
    checkData.vars(dat, adj.vars)

    # Convert positions to vars
    response.var <- convertPosToVars(dat, response.var)
    snp.vars     <- convertPosToVars(dat, snp.vars)
    adj.vars     <- convertPosToVars(dat, adj.vars)

    if (!is.null(method)) {
      if (length(method) != 1) stop("ERROR: The method option must be NULL, case.control, or case.complement")
      method <- tolower(method)
      if (method == -1) {
        pool <= -1
      } else {
        temp <- c("case-control", "case-complement") %in% method
        if (!any(temp)) stop("ERROR: The method option must be NULL, case-control, or case-complement")
        pool <- temp[2]
      }
    } else {
      pool <- NULL
    }
    meth <- test.type
	if(!is.null(subset) && sum(subset) == 0) stop("Empty subset of rows")
	meth.pval <- meth.pval[1]

    # Check cntl.lab, types.lab. 
    if (is.null(subset)) {
      labs <- unique(dat[, response.var])
    } else {
      labs <- unique(dat[subset, response.var])
    }

    if (!(cntl.lab %in% labs)) stop(paste("ERROR: cntl.lab = ", cntl.lab, " is not a valid label", sep=""))
    if (is.null(types.lab)) types.lab <- labs[!(labs %in% cntl.lab)]
    temp      <- !(types.lab %in% cntl.lab)
    types.lab <- types.lab[temp]
    k         <- length(types.lab)
    if (!k) stop("ERROR: with types.lab. No disease subtypes are being included.")
    temp <- types.lab %in% labs
    if (!all(temp)) {
      temp <- types.lab[!temp]
      print(temp)
      stop("ERROR: with types.lab. The above disease subtypes are not valid.")
    }
    if (is.null(subset)) {
      subset <- dat[, response.var] %in% c(cntl.lab, types.lab)
    } else {
      subset <- subset & (dat[, response.var] %in% c(cntl.lab, types.lab))
    }
    
	pval.l <- pval1 <- pval2 <- NULL
	par.l <- par1 <- par2 <- NULL
	sigma.l <- sigma1 <- sigma2 <- NULL
	spheno1 <- spheno2 <- NULL
	
	SC <- (is.null(pool) || pool==FALSE)
	CC <- (is.null(pool) || pool==TRUE)

     logit.all <- CC.all <- SC.all <- list()
     NULL.pheno <- rep(FALSE, k)

     for (ii in 1:length(snp.vars)) { 
	  logit.res <- CC.res <- SC.res <- NULL
       SNP <- snp.vars[ii]
	  if(logit)
	  { 
		res <- try(types.wald(sub=rep(TRUE, k), snp.vars=SNP, dat=dat, response.var=response.var, adj.vars=adj.vars
					, types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, pool=FALSE, side=side))
		if(inherits(res, "try-error")) warning(paste("Error in Overall Logistic:", res))
		else logit.res <- res
	  }

	  if(SC)
	  { 
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 0))
		else pval.args$search <- 0

		res <- try(h.types0(k, dat, response.var=response.var, snp.vars=SNP, adj.vars=adj.vars, types.lab=types.lab
						, cntl.lab=cntl.lab, subset=subset, pool=FALSE, meth = meth
							, side=side, meta.C = TRUE, zmax.args=zmax.args
							, meth.pval=meth.pval, pval.args=pval.args, p.bound=p.bound, add.enrich=add.enrich,
                                       cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0))
		
		if(inherits(res, "try-error")) warning(paste("Error in Subset-search (Case-Control):", res))
		else SC.res <- res
		
	  }
	  if(CC)
	  { 
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 0))
		else pval.args$search <- 0
		
		res <- h.types0(k, dat, response.var=response.var, snp.vars=SNP, adj.vars=adj.vars, types.lab=types.lab
						  , cntl.lab=cntl.lab, subset=subset, pool=TRUE, meth = meth
						, side=side, meta.C = TRUE, zmax.args=zmax.args
						, meth.pval=meth.pval, pval.args=pval.args, p.bound=p.bound, add.enrich=add.enrich,
                                  cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)

		if(inherits(res, "try-error")) warning(paste("Error in Subset-search (Case-Control):", res))
		else CC.res <- res
	  }
      
       # Combine results
       logit.all <- htypes.combine(logit.all, logit.res, 1, SNP, k) 
       SC.all    <- htypes.combine(SC.all, SC.res, 2, SNP, k)
       CC.all    <- htypes.combine(CC.all, CC.res, 2, SNP, k)

     }    

	ret <- list(Overall.Logistic=logit.all, Subset.Case.Control=SC.all, Subset.Case.Complement=CC.all,
                data=dat, response.var=response.var, adj.vars=adj.vars, types.lab=types.lab,
                cntl.lab=cntl.lab, subset=subset, method=method, side=side, test.type=test.type,
                zmax.args=zmax.args, meth.pval=meth.pval, pval.args=pval.args, logit=logit,
                snp.vars=snp.vars, which="h.types", cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)

	ret
}

# Function to combine results
htypes.combine <- function(all, new, which, SNP, k) {

  if (is.null(new)) {
    x        <- NA
    names(x) <- SNP
    z <- zopt <- beta <- sd <- pval <- x
    pheno <- matrix(data=FALSE, nrow=1, ncol=k)
    rownames(pheno) <- SNP 
  } else {
    z     <- new[["z", exact=TRUE]]
    zopt  <- new[["zopt", exact=TRUE]]
    beta  <- new[["beta", exact=TRUE]]
    sd    <- new[["sd", exact=TRUE]]
    pval  <- new[["pval", exact=TRUE]]
    pheno <- new[["pheno", exact=TRUE]]
  }
  
  if (which == 1) {
    all$z    <- c(all$z, z)
    all$beta <- c(all$beta, beta)
    all$sd   <- c(all$sd, sd)
    all$pval <- c(all$pval, pval)
  } else {
    all$pval  <- c(all$pval, pval)
    all$beta  <- c(all$beta, beta)
    all$sd    <- c(all$sd, sd)
    all$zopt  <- c(all$zopt, zopt)
    all$pheno <- rbind(all$pheno, pheno)
  } 
 
  all

} # END: htypes.combine

h.types0 <- function(k, dat, response.var, snp.vars, adj.vars, types.lab, cntl.lab, pool, meth = "score"
, subset=NULL, side=2, meta.C = FALSE, zmax.args=NULL, meth.pval="DLM", pval.args=NULL, p.bound=1, 
   add.enrich=FALSE, cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{
	
	if(k == 0) return(list(pval = 1, pheno = NULL, beta = NA, sd = NA, zopt=NA, subsetCount=0))

	k <- length(types.lab)
	nsnp <- length(snp.vars)
	N <- nrow(dat)
		
	pheno <- matrix(data=FALSE, nrow=1, ncol=k)
     rownames(pheno) <- snp.vars
	colnames(pheno) <- types.lab

	pval <- rep(NA, nsnp) ; names(pval) <- snp.vars
	par <- pval ; sigma <- pval
     zopt <- pval

	z0 <- pvals0 <- rep(NA, k)
	for(t in 1:k)
	{
		meta.def <- switch(meth, Score=types.score, Wald=types.wald)
		res <- meta.def(sub=((1:k)==t), snp.vars=snp.vars[1], dat=dat, response.var=response.var, 
              adj.vars=adj.vars, types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, pool=pool)
		z0[t] <- res$z
		pvals0[t] <- if(abs(side) == 1) pnorm(abs(z0[t]), lower.tail=FALSE) else 2 * pnorm(abs(z0[t]), lower.tail=FALSE)
	}
	
	z.sub <- (pvals0 <= p.bound) * sign(z0)
	k1 <- sum(z.sub != 0)

	if(k1 == 0) return(list(pval=pval, pheno=pheno, beta=par, sd=sigma, zopt=zopt, subsetCount=0))

     # Determine if the user passed in subset functions
     ###############################################
     sub.def  <- zmax.args[["sub.def", exact=TRUE]] 
     sub.args <- zmax.args[["sub.args", exact=TRUE]]
     psub.def  <- pval.args[["sub.def", exact=TRUE]] 
     psub.args <- pval.args[["sub.args", exact=TRUE]]
     ###############################################

     # For zmax
     if (!is.null(sub.def)) {
	  sub.def1 <- function(x, op)
	  {
           z.sub <- op$z.sub
		x0    <- (z.sub != 0) & x
		ret   <- sub.def(x0, op)
		ret
	  }
	  if(k1!=k)
  	  {
           zmax.args$sub.def  <- sub.def1
           zmax.args$sub.args <- c(sub.args, list(z.sub=z.sub))
	  }
     }

     # For p.dlm or p.tube
     if (!is.null(psub.def)) {
	  psub.def1 <- function(x, op)
	  {
           x0          <- (op$z.sub != 0)
           pos         <- which(x0)
		x0[pos[!x]] <- FALSE
		ret         <- psub.def(x0, op)
		ret
	  }
	  if(k1!=k)
	  {
           pval.args$sub.def  <- psub.def1
           pval.args$sub.args <- c(psub.args, list(z.sub=z.sub))
	  }
     }


	if(k1 != k && !is.null(pval.args$sizes))
	{
		sizes <- pval.args$sizes
		labs <- rep(1:length(sizes), sizes)
		sizes1 <- NULL
		for(i in 1:length(sizes))
		{
			num.skip <- sum((z.sub == 0) & labs == i)
			sizes1 <- c(sizes1, c(sizes[i] - num.skip))
		}
	pval.args$sizes <- sizes1
	}

	if(is.null(zmax.args) || is.null(zmax.args$z.sub)) zmax.args <- c(zmax.args, list(z.sub = z.sub))
	if(is.null(pval.args) || is.null(pval.args$p.bound)) pval.args <- c(pval.args, list(p.bound = p.bound))
	
	meta.def <- switch(meth, Score=types.score, Wald=types.wald)

	res <- do.call(z.max, c(list(k=k, snp.vars=snp.vars, side=side, meta.def = meta.def
				, meta.args = list(dat=dat, response.var=response.var, adj.vars = adj.vars
				, types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, pool=pool)), zmax.args))

	zopt <- as.double(res$opt.z)
	pheno <- matrix(as.logical(res$opt.s), ncol=k, byrow=FALSE)
	subsetCount <- res$subsetCount
	
	names(zopt) <- snp.vars
	dimnames(pheno)[[1]] <- as.list(snp.vars)
	dimnames(pheno)[[2]] <- as.list(types.lab)

	
	ncase <- matrix(0, k, nsnp)
	ncntl <- rep(0, nsnp)

	if(is.null(subset)) subset <- rep(TRUE, nrow(dat))
	subset1 <- findMiss.vars(dat, vars=c(response.var, adj.vars), miss = NA)
	subset0 <- if(is.null(subset)) subset1 else (subset & subset1)
	
	for(j in 1:nsnp)	
	{
		subset0g <- subset0 & findMiss.vars(dat, vars=snp.vars[j], miss = NA)
		
		tab <- table(dat[subset0g, response.var], useNA = "no")
		pos <- pmatch(c(cntl.lab, types.lab), names(tab), nomatch = 0)

		ncase[, j] <- ifelse(pos[-1] > 0, tab[pos[-1]], 0)
		ncntl[j] <- if(pos[1] > 0) tab[pos[1]] else 0
	}

	if(meth.pval == "DLM")
	{
		#if(is.null(pval.args) || !("cor.def" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.def=NULL))
	     #if(!("cor.args" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.args = list(ncase=matrix(ncase[(z.sub != 0), ], ncol=nsnp), ncntl=ncntl, pool=pool)))
		#pval <- do.call(p.dlm, c(list(t.vec=abs(zopt), z.sub=z.sub, side=abs(side), cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0), pval.args))
          cor.args <- pval.args[["cor.args", exact=TRUE]]
          if (is.null(cor.args))  cor.args <- list(ncase=matrix(ncase[(z.sub != 0), ], ncol=nsnp), ncntl=ncntl, pool=pool)

          pval <- p.dlm(abs(zopt), z.sub, 0, abs(side), cor.def=pval.args[["cor.def", exact=TRUE]], 
                  cor.args=cor.args, sizes=pval.args[["sizes", exact=TRUE]], p.bound=p.bound, 
                  sub.def=pval.args[["sub.def", exact=TRUE]], sub.args=pval.args[["sub.args", exact=TRUE]],
                  NSAMP=NSAMP, NSAMP0=NSAMP0)

	}

	if(meth.pval == "Boot")
	{
		p.boot(t0=abs(zopt), z.sub=z.sub, search=0, side = abs(side), ncase0=ncase, ncntl0=ncntl, pool=pool, p.bound=p.bound,
                  cond.all=cond.all, NSAMP0=NSAMP0)
	}


	if(meth.pval == "IS") pval <- do.call(p.tube, c(list(t.vec=abs(zopt), k=k1
									, side = abs(side), ncase=matrix(ncase[(z.sub != 0), ], ncol=nsnp), ncntl=ncntl
									, pool=pool), pval.args))
	
	if(meth.pval == "B") pval <- p.bon(abs(zopt), subsetCount, search = 0, side = abs(side))
	
	names(pval) <- snp.vars
	
	par <- sigma <- NULL

#	Estimate Parameters
	if(any(pheno))
	{
		res <- types.wald(sub=pheno, snp.vars=snp.vars, dat=dat, response.var=response.var
							, adj.vars=adj.vars, types.lab=types.lab, cntl.lab=cntl.lab, pool=pool, side = side)
		par <- res$beta
		names(par) <- snp.vars
		if(side == 2) sigma <- abs(par)/qnorm(pval/2, lower.tail=FALSE, log.p=FALSE)
		else sigma <- abs(par)/qnorm(pval, lower.tail=(side == -1), log.p=FALSE)
	}
	if(add.enrich && p.bound < 1)
	{
		pval1 <- calcP1(p.bound, k1, k, search=0, side=side, rmat=diag(1, k), sizes=sizes)
		pval.c  <- pchisq(-2 * (log(pval)+ log(pval1)), df=4, lower.tail=FALSE)
		pval <- pval.c
	}
	list(pval=pval, beta=par, sd=sigma, zopt = zopt, pheno=pheno)
}

types.score <- function(sub, snp.vars, dat, response.var, adj.vars, types.lab, cntl.lab, subset=NULL, pool=FALSE)
{
	
	sub.case <- types.lab[sub]
	sub.cntl <- if(pool == FALSE) cntl.lab else c(types.lab[!sub], cntl.lab)
	
	nsnp <- length(snp.vars)
	
#	valid row numbers
	subset1 <- dat[, response.var] %in% c(sub.case, sub.cntl)
	subset01 <- if(is.null(subset)) subset1 else subset & subset1
	subset01x <- subset01 & findMiss.vars(dat, vars = c(response.var, adj.vars), miss = NA)
	
	if(is.null(adj.vars)) p <- 1
	else p <- length(adj.vars) + 1
	geno.flag <- (p == 1 && all(dat[, snp.vars] %in% c(0, 1, 2, NA)))

	if(is.null(adj.vars)) adj.vars <- ""
	if(length(adj.vars) > 1) adj.vars <- paste(adj.vars, collapse="+")	
	fmla <- "~1"
	if(adj.vars != "") fmla <- paste(fmla, "+" , adj.vars, sep="")
	xmat <- model.matrix(as.formula(fmla), dat[subset01x, ], na.action="na.pass")
	
	N <- sum(subset01x)
	
	d.vec <- rep(NA, N)
	d.vec[dat[subset01x, response.var] %in% sub.case] <- 1
	d.vec[dat[subset01x, response.var] %in% sub.cntl] <- 0
	
	
	p.vec <- rep(-1, N)
	if(geno.flag) { 
		mat <- matrix(table(d.vec), ncol=2)
		res <- try(glm(mat ~ 1, family=binomial(link="logit")))
	} else res <- try(glm(d.vec ~ xmat, family=binomial(link="logit")))
	
	if(inherits(res, "try-error") || (!res$converged))
	{
		warning("error in glm")
		p.vec <- rep(mean(d.vec, na.rm = TRUE), N)
	} else
	{
		if(geno.flag) p.vec <- rep(res$fitted, N)
		else p.vec <- res$fitted
	}
	
	g.vec <- as.vector(data.matrix(dat[subset01x, snp.vars]))
	nmiss.ind <- !(is.na(g.vec))
	 
#	nmiss.ind <- matrix(!tmp, ncol=nsnp)
#	g.vec <- as.vector(g.mat)
	
	xinfo <- t(xmat) %*% ( p.vec * (1 - p.vec) * xmat)
	xcov <- myInv(xinfo)
	
	res <- score.compute(d.vec, p.vec, g.vec, xmat, nmiss.ind, N, nsnp, p, xcov=xcov)
	numr <- res$numr
	denr <- res$denr
	
#	
#	If this flag is set, all SNPs with negative denominator are refit (using 
#	individuals with non-missing genotype).
#	
	REFIT <- FALSE
	if(REFIT && any(denr < 0))
	{
		for(j in which(denr < 0))
		{
			warning(paste("Negative denominator before refitting:", snp.vars[j]))
			g.vec.j <- g.vec[, (j - 1) * N + (1:N)]
			nmiss.j <- nmiss.ind[(j - 1) * N + (1:N)]
			NN <- sum(nmiss.j)
			p.vec.j <- rep(-1, N)
			if(geno.flag) { 
				mat <- matrix(table(d.vec[nmiss.j]), ncol=2)
				res <- try(glm(mat ~ 1, family=binomial(link="logit")))
			} else res <- try(glm(d.vec ~ xmat, subset = nmiss.j, family=binomial(link="logit")))
			
			if(inherits(res, "try-error") || (!res$converged)) warning("error in glm")
			else
			{
				if(geno.flag) p.vec[nmiss.j] <- rep(res$fitted, NN)
				else p.vec[nmiss.j] <- res$fitted
			}
			
			res <- score.compute(d.vec, p.vec.j, g.vec.j, xmat, nmiss.j, N, 1, p, xcov=NULL)
			numr[j] <- res$numr
			denr[j] <- res$denr
			if(denr[j] < 0) warning(paste("Negative denominator after refitting:", snp.vars[j]))
		}
	}
#	print(c(sub))
#	print(rbind(numr, denr))
	
    # Changed July 31, 2013
	#ztmp <- ifelse(is.na(denr) | is.nan(denr) | denr <= 0, 0, numr/sqrt(denr))	
    ztmp <- ifelse(is.na(denr) | is.nan(denr), 0, numr/sqrt(denr))
	list(z = ztmp)
}

score.compute <- function(d.vec, p.vec, g.vec, xmat, nmiss, N, nsnp, p, xcov=NULL)
{
	q.vec <- (1 - p.vec)
	tmp.mat <- matrix(g.vec * nmiss * (d.vec - p.vec), N, nsnp, byrow=FALSE)
	numr0 <- apply(tmp.mat, 2, sum, na.rm=TRUE)
	tmp.mat <- matrix(nmiss * g.vec * g.vec * p.vec * q.vec, N, nsnp, byrow=FALSE)
	denr0 <- apply(tmp.mat, 2, sum, na.rm=TRUE)
	
	XtDg <- matrix(sapply(1:p, function(j)
						  {
						  tmp.mat <- matrix(g.vec * nmiss * p.vec * q.vec * xmat[, j], N, nsnp, byrow=FALSE)
						  apply(tmp.mat, 2, sum, na.rm=TRUE)				   
						  }), ncol=p, byrow=FALSE)
	XtD <- matrix(sapply(1:p, function(j)
						 {
						 tmp.mat <- matrix(nmiss * (d.vec - p.vec) * xmat[, j], N, nsnp, byrow=FALSE)
						 apply(tmp.mat, 2, sum, na.rm=TRUE)
						 }), ncol=p, byrow=FALSE)
	rm(tmp.mat)
	gc()
	
	if(! is.null(xcov))
	{
		numr <- numr0 - sapply(1:nsnp, function(j){ 
							   crossprod(matrix(XtDg[j, ], ncol = 1), xcov %*% matrix(XtD[j, ], ncol=1))
							   })
		
		denr <- denr0 - apply(XtDg, 1, function(x){ 
							  mx <- matrix(x, ncol=1)
							  crossprod(mx, xcov %*% mx)
							  })
	}
	else
	{
		xinfo <- t(xmat) %*% (nmiss * p.vec * q.vec * xmat)
		numr <- numr0 - sapply(1:nsnp, function(j){ crossprod(matrix(XtDg[j, ], ncol = 1), mySolve(xinfo, matrix(XtD[j, ], ncol=1)))})
		denr <- denr0 - apply(XtDg, 1, function(x){ mx <- matrix(x, ncol=1) ; crossprod(mx, mySolve(xinfo, mx))})
	}
	return(list(numr=numr, denr=denr))
}

types.wald <- function(sub, snp.vars, dat, response.var, adj.vars, types.lab, cntl.lab, subset=NULL, pool=FALSE, side=2)
{
	nsnp <- length(snp.vars)
	k <- length(types.lab)
	N <- nrow(dat)

	if(length(sub) == (k * nsnp))
	{
		if(k == 1 && is.null(dim(sub))) dim(sub) <- c(nsnp, 1)
		if(nsnp == 1 && is.null(dim(sub))) dim(sub) <- c(1, k)
		if(nrow(sub) != nsnp || ncol(sub) != k) stop("sub should be of dimension nsnp X k")
	}
	if(is.null(dim(sub)) && length(sub) == k) sub <- matrix(sub, nsnp, k, byrow=TRUE)		
	
	sub.case <- types.lab[sub[1, ]]
	sub.cntl <- if(pool == FALSE) cntl.lab else c(types.lab[!sub[1, ]], cntl.lab)
	
	if(is.null(adj.vars)) p <- 1
	else p <- length(adj.vars) + 1
	geno.flag <- (p == 1 && all(dat[, snp.vars] %in% c(0, 1, 2, NA)))
	
	if(is.null(adj.vars)) adj.form <- ""
	if(is.null(snp.vars)) snp.form <- ""
	
  ##SMB: changed this from > 1 to >= 1, since I was getting an error when there
  ##was only one variable that needed adjusting for
	if(length(adj.vars) >= 1) adj.form <- paste(adj.vars, collapse="+")
	##if(length(adj.vars) > 1) adj.form <- paste(adj.vars, collapse="+")
	
	fmla <- paste(response.var, "~ 1", sep="")					 
	if(snp.vars[1] != "") fmla <- paste(fmla, "+" , snp.vars[1], sep="")
	if(adj.form != "") fmla <- paste(fmla, "+" , adj.form, sep="")
	ndat <- model.frame(as.formula(fmla), dat, na.action="na.pass")
	
	dx <- dat[,response.var]
	dx1 <- rep(0, N)
	
	dx1[dx %in% sub.cntl] <- 0
	dx1[dx %in% sub.case] <- 1

	subset1 <- (dx %in% c(sub.cntl, sub.case))
	subset01 <- if(is.null(subset)) subset1 else subset & subset1
	
	dx1[!subset01] <- NA
	
	ndat[, response.var] <- dx1
	ndat[, response.var] <- as.integer(ndat[, response.var])
	
	stat <- sapply(1:nsnp, function(j)
					 {

					if(nrow(sub) == nsnp)
					{
						sub.case <- types.lab[sub[j, ]]
						sub.cntl <- if(pool == FALSE) cntl.lab else c(types.lab[!sub[j, ]], cntl.lab)

						dx1 <- rep(0, N)
						dx1[dx %in% sub.cntl] <- 0
						dx1[dx %in% sub.case] <- 1
						subset1 <- (dx %in% c(sub.cntl, sub.case))
						subset01 <- if(is.null(subset)) subset1 else subset & subset1
						dx1[!subset01] <- NA
						ndat[, response.var] <- dx1
					}
					snp <- snp.vars[j]
					g.vec <- dat[, snp]
					ndat[,2] <- g.vec
					colnames(ndat)[2] <- snp
					 
					fmla <- paste(response.var, "~ 1", sep="")
					if(snp != "") fmla <- paste(fmla, "+" , snp, sep="")
					if(adj.form != "") fmla <- paste(fmla, "+" , adj.vars, sep="")

					
					ret <- c(NA, NA)
					if(geno.flag)
					{
						mat <- table(g.vec, 1 - dx1, useNA="no") ; geno <- (0:2)
						res <- try(glm(mat ~ geno, family=binomial(link="logit")))
					} else {
                                          res <- try(glm(as.formula(fmla), data=ndat, subset = subset01, family=binomial(link="logit")))
                                   }

					if(inherits(res, "try-error") || (!res$converged)) 
                                   {
                                          warning("error in glm")
                                   } else {
						coef <- summary(res)$coef
						pos <- pmatch(if(geno.flag) "geno" else snp, rownames(coef))
						if(!is.na(pos)) ret <- coef[pos, 1:2]
					}
					ret
					})

	colnames(stat) <- snp.vars

	beta <- stat[1, ]
	names(beta) <- snp.vars
	
	sd <- stat[2, ]
	names(sd) <- snp.vars
	
	z <- ifelse(is.na(sd) | is.nan(sd) | sd <= 0, 0, beta/sd)

	if(side == 2) pval <- 2 * pnorm(abs(z), lower.tail=FALSE)
	if(side == 1) pval <- pnorm(z, lower.tail=FALSE)
	if(side == -1) pval <- pnorm(z, lower.tail=TRUE)
	
	list(z=z, beta=beta, sd=sd, pval=pval)
}

types.forest0 <- function(dat, snp.var, response.var, adj.vars, types.lab, cntl.lab, subset=NULL, rlist = NULL
, level = 0.05, p.adj = TRUE, digits = 2)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")
	if(!(snp.var %in% colnames(dat))) stop("colnames of dat does not have snp.var")
	k <- length(types.lab)
	nsub <- matrix(FALSE, k, k)
	diag(nsub) <- TRUE
	
	if(!is.null(rlist))
	{
		if(is.null(rlist$Overall.Logistic)) stop("Missing Component in rlist: Overall.Logistic")
		if(is.null(rlist$Subset.Case.Control)) stop("Missing Component in rlist: Subset.Case.Control")
		if(is.null(rlist$Subset.Case.Complement)) stop("Missing Component in rlist: Subset.Case.Complement")
	} else
	{
		rlist <- h.types(dat=dat, snp.vars=snp.var, response.var=response.var, adj.vars=adj.vars
						 , types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, method = NULL,
                         logit = TRUE, test.type = "Wald")
	}
 	
	res <- 	types.wald(sub = nsub, snp.vars=rep(snp.var, k), dat=dat, response.var=response.var, adj.vars=adj.vars
					   , types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, pool=FALSE)

	h.forest(k, snp.var, types.lab, rlist, res, side = 1, level=level, p.adj = p.adj, digits=digits)
}

types.forest <- function(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")
	if(!(snp.var %in% colnames(rlist$data))) stop("colnames of data does not have snp.var")
    if (!(snp.var %in% rlist$snp.vars)) stop("snp.var was not a SNP analyzed")

    types.lab    <- rlist$types.lab
    cntl.lab     <- rlist$cntl.lab
    response.var <- rlist$response.var
    adj.vars     <- rlist$adj.vars
    subset       <- rlist$subset
    method       <- rlist$method
    side         <- rlist$side
    test.type    <- rlist$test.type
    zmax.args    <- rlist$zmax.args
    meth.pval    <- rlist$meth.pval
    pval.args    <- rlist$pval.args 
    cond.all     <- rlist$cond.all
    NSAMP        <- rlist$NSAMP
    NSAMP0       <- rlist$NSAMP0

	k <- length(types.lab)
	nsub <- matrix(FALSE, k, k)
	diag(nsub) <- TRUE
	
    ov      <- rlist[["Overall.Logistic", exact=TRUE]]
    cc      <- rlist[["Subset.Case.Control", exact=TRUE]]
    cp      <- rlist[["Subset.Case.Complement", exact=TRUE]]
    ov.flag <- !is.null(ov)
    cc.flag <- !is.null(cc)
    cp.flag <- !is.null(cp)
    if ((!ov.flag) || (!cc.flag) || (!cp.flag)) {
      if (ov.flag) {
        logit <- FALSE
      } else {
        logit <- TRUE
      }

      method <- "-1"
      if ((!cc.flag) && (!cp.flag)) {
        method <- NULL
      } else if (!cc.flag) {
        method <- "case-control"
      } else if (!cp.flag) {
        method <- "case-complement"
      }

      ret <- h.types(rlist$data, response.var, snp.var, adj.vars, types.lab, cntl.lab, subset=subset,
              method=method, side=side, logit=logit, test.type=test.type, zmax.args=zmax.args, 
              pval.args=pval.args, NSAMP=NSAMP, NSAMP0=NSAMP0)
      if (!ov.flag) ov <- ret[["Overall.Logistic", exact=TRUE]]
      if (!cc.flag) cc <- ret[["Subset.Case.Control", exact=TRUE]]
      if (!cp.flag) cp <- ret[["Subset.Case.Complement", exact=TRUE]]
    }
    
    newlist <- list(Overall.Logistic=ov, Subset.Case.Control=cc, Subset.Case.Complement=cp)
    
	res <- 	types.wald(sub = nsub, snp.vars=rep(snp.var, k), dat=rlist$data, response.var=response.var, 
               adj.vars=adj.vars, types.lab=types.lab, cntl.lab=cntl.lab, subset=subset, pool=FALSE)

	h.forest(k, snp.var, types.lab, newlist, res, side = 1, level=level, p.adj = p.adj, digits=digits)
}


