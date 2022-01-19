
h.traits <- function(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=NULL
	, cor.numr=FALSE, search=NULL, side=2, meta=FALSE, zmax.args=NULL
	, meth.pval="DLM", pval.args=NULL)
{
    if (!is.null(search)) {
      if (!(search %in% c(1, 2))) stop("search must be NULL, 1 or 2")
    }
    if ((is.null(side)) || (!(side %in% c(1, 2)))) stop("side must be 1 or 2")
    if ((is.null(meth.pval)) || (!(meth.pval %in% c("DLM", "IS", "B")))) stop("meth.pval must be DLM, IS, or B")
    if ((is.null(meta)) || (!(meta %in% 0:1))) stop("meta must be TRUE or FALSE")
    if ((is.null(cor.numr)) || (!(cor.numr %in% 0:1))) stop("cor.numr must be TRUE or FALSE")
    if (!is.null(pval.args)) {
      if (!is.list(pval.args)) stop("pval.args must be a named list")
      temp <- names(pval.args)
      if (is.null(temp)) stop("pval.args must be a named list")
    }
    if (!is.null(zmax.args)) {
      if (!is.list(zmax.args)) stop("zmax.args must be a named list")
      temp <- names(zmax.args)
      if (is.null(temp)) stop("zmax.args must be a named list")
    }

	k         <- length(traits.lab)
	meth.pval <- meth.pval[1]
	nsnp      <- length(snp.vars)
	beta.hat  <- matrix(beta.hat, ncol=k, byrow=FALSE)
	sigma.hat <- matrix(sigma.hat, ncol=k, byrow=FALSE)
	
       #SB: traits.forest calling for 1 snp with whole matrix. So this was generating error.
       #if(nrow(beta.hat) > nsnp) stop("nrow(beta.hat) must equal length(snp.vars)")	
       #if(nrow(sigma.hat) > nsnp) stop("nrow(sigma.hat) must equal length(snp.vars)")

	if(nrow(beta.hat) > nsnp) beta.hat <- matrix(beta.hat[snp.vars, drop=FALSE], nrow = nsnp, ncol=k, byrow=FALSE)
	if(nrow(sigma.hat) > nsnp) sigma.hat <- matrix(sigma.hat[snp.vars, drop=FALSE], nrow = nsnp, ncol=k, byrow=FALSE)

	rownames(beta.hat) <- snp.vars
	rownames(sigma.hat) <- snp.vars
	colnames(beta.hat) <- traits.lab
	colnames(sigma.hat) <- traits.lab
	
	if(!((length(ncase) == k) && is.null(dim(ncase))) && !(is.matrix(ncase) && (ncol(ncase) == k) && (nrow(ncase) == nsnp)))
		stop("ncase should be a vector of size k or a matrix of dimension nsnp X k")
	if(is.null(dim(ncase))) dim(ncase) <- c(1, k)
	if(!((length(ncntl) == k) && is.null(dim(ncntl))) && !(is.matrix(ncntl) && (ncol(ncntl) == k) && (nrow(ncntl) == nsnp)))
		stop("ncntl should be a vector of size k or a matrix of dimension nsnp X k")
	if(is.null(dim(ncntl))) dim(ncntl) <- c(1, k)

	if(is.null(cor)) cor <- diag(k)
	if(!is.matrix(cor) && !is.list(cor)) stop("cor should be a matrix or a list")
	if(is.list(cor))
	{
		#if(!all(c("N11", "N00", "N10", "N01") %in% names(cor)))
	    #		stop("Expected four matrices in list cor with names N11, N00, N10 and N01")
		#if(!all(vapply(cor[c("N11", "N00", "N10", "N01")]
			#, function(m) { is.matrix(m) && nrow(m) == k  && ncol(m) == k }, TRUE)))
	    # 		stop("N11, N00, N10 and N01 should be square matrices of dimension k")

       # Changed Aug 08 2013
       if(!all(c("N11", "N00", "N10") %in% names(cor))) {
	     stop("Expected three matrices in list cor with names N11, N00 and N10")
       }
	   if(!all(vapply(cor[c("N11", "N00", "N10")]
		, function(m) { is.matrix(m) && nrow(m) == k  && ncol(m) == k }, TRUE)))
	     		stop("N11, N00 and N10 should be square matrices of dimension k")
       cor$N01 <- t(cor$N10)

	   cor <- corr.mat.logit(cor$N11, cor$N00, cor$N10, cor$N01)
	}

	if(!is.null(colnames(cor)) && any(colnames(cor) != traits.lab)) stop("cor should have column names in same order as traits.lab")
	if(!is.null(rownames(cor)) && any(rownames(cor) != traits.lab)) stop("cor should have row names in same order as traits.lab")
	
	
	pval <- pval1 <- pval2 <- rep(NA, nsnp)
	names(pval) <- names(pval1) <- names(pval2) <- snp.vars

	beta <- beta1 <- beta2 <- rep(NA, nsnp)
	names(beta) <- names(beta1) <- names(beta2) <- snp.vars

	sd <- sd1 <- sd2 <- rep(NA, nsnp)
	names(sd) <- names(sd1) <- names(sd2) <- snp.vars
	
	pheno <- pheno1 <- pheno2 <- matrix(FALSE, nsnp, k)
	colnames(pheno) <- colnames(pheno1) <- colnames(pheno2) <- traits.lab
	rownames(pheno) <- rownames(pheno1) <- rownames(pheno2) <- snp.vars

    # Vectors for new meta-anlaysis standard errors
	sd.meta <- sd1.meta <- sd2.meta <- rep(NA, nsnp)
	names(sd.meta) <- names(sd1.meta) <- names(sd2.meta) <- snp.vars


	OS <- (is.null(search) || search==1)
	TS <- (is.null(search) || search==2)
	meta.res <- sub1.res <- sub2.res <- NULL
	if(meta)
	{ 
		for(j in seq_len(nsnp))
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
				nsub <- sum(sub)
				rmat <- matrix(cor[sub, sub], nsub, nsub)
				res <- traits.meta(sub, snp.vars[j], beta.hat[j, ], sigma.hat[j,]
							, ncase[jj, ], ncntl[jj, ], rmat = cor, side = side
							, cor.numr = cor.numr, wt.sigma = TRUE)
				pval[j] <- res$pval
				beta[j] <- res$beta
				sd[j] <- res$sd

            			# Compute new meta-analysis standard error
            			#sd.meta[j] <- meta.se(cor, sigma.hat[j, ])
			}
		}
		meta.res <- list(pval=pval, beta=beta, sd=sd)
	}

	subset0 <- data.frame(pval=rep(NA, nsnp), pheno=rep("NA", nsnp), stringsAsFactors=FALSE)
	subset00 <- data.frame(pval=rep(NA, nsnp), pheno1=rep("NA", nsnp), pheno2=rep("NA", nsnp), 
                           stringsAsFactors=FALSE)

	if(OS)
	{ 	
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 1))
		else pval.args$search <- 1

		for(j in seq_len(nsnp))
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
				nsub <- sum(sub)
				rmat <- matrix(cor[sub, sub], nsub, nsub)
				res <- h.traits1(nsub, beta.hat[j, sub], sigma.hat[j, sub], ncase[jj, sub], ncntl[jj, sub]
							, rmat = rmat, cor.numr=cor.numr, side = side, zmax.args=zmax.args
						 	, meth.pval=meth.pval, pval.args=pval.args)

				pval[j]       <- res$pval
				pheno[j, sub] <- res$pheno
				beta[j]       <- res$beta
				sd[j]         <- res$sd

				# Compute new meta-analysis standard error
				sd.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res$pheno)

			}
			
		}
		sub1.res <- list(pval=pval, beta=beta, sd=sd, pheno=pheno, sd.meta=sd.meta)
	}

	if(TS)
	{ 
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 2))
		else pval.args$search <- 2

		for(j in seq_len(nsnp))
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			res <- NULL
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
				nsub <- sum(sub)
				rmat <- matrix(cor[sub, sub], nsub, nsub)

				res <- h.traits2(nsub, beta.hat[j, sub], sigma.hat[j, sub], ncase[jj, sub], ncntl[jj, sub], 
                                 rmat=rmat, cor.numr=cor.numr, side = side, zmax.args=zmax.args, 
                                 meth.pval=meth.pval, pval.args=pval.args)

                pval[j]        <- res$pval
			    pval1[j]       <- res$pval.1
		     	pval2[j]       <- res$pval.2
			    pheno1[j, sub] <- res$pheno.1
			    pheno2[j, sub] <- res$pheno.2
			    beta1[j]       <- res$beta.1
			    beta2[j]       <- res$beta.2
			    sd1[j]         <- res$sd.1
			    sd2[j]         <- res$sd.2

                # Compute new meta-analysis standard error
                sd1.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res$pheno.1)
                sd2.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res$pheno.2)

			}
			
			
		}

		sub2.res <- list(pval=pval, pval.1 = pval1, beta.1=beta1, sd.1=sd1, pheno.1=pheno1,
						 pval.2=pval2, beta.2=beta2, sd.2=sd2, pheno.2=pheno2,
                         sd.1.meta=sd1.meta, sd.2.meta=sd2.meta)
	}

	#ret <- list(Meta=meta.res, Subset.1sided=sub1.res, Subset.2sided=sub2.res)
    ret <- list(Meta=meta.res, Subset.1sided=sub1.res, Subset.2sided=sub2.res,
           snp.vars=snp.vars, traits.lab=traits.lab, beta.hat=beta.hat, sigma.hat=sigma.hat,
           ncase=ncase, ncntl=ncntl, cor=cor, cor.numr=cor.numr, search=search,
           side=side, meta=meta, zmax.args=zmax.args, meth.pval=meth.pval,
           pval.args=pval.args, which="h.traits")
 
	ret
}


h.traits1 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, side=2, zmax.args=NULL, meth.pval="DLM", pval.args=NULL)
{
	if(k == 0) return(list(pval = 1, pheno = NULL, beta = NA, sd = NA, subsetCount=0))
	
##### NOTE: meta.def is set to traits.meta ######
	res <- do.call(z.max, c(list(k, 1, side=side, meta.def=traits.meta, meta.args=list(beta.hat=beta.hat, sigma.hat=sigma.hat
				 , ncase=ncase, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr, wt.sigma=FALSE)), zmax.args))

	zopt <- as.double(res$opt.z)
	pheno <- as.logical(res$opt.s)
	subsetCount <- res$subsetCount
	rm(res)

	if(meth.pval == "DLM")
	{
		if(is.null(pval.args) || !("cor.def" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.def=NULL))
		if(!("cor.args" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.args = list(ncase=ncase, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr)))

		pval <- do.call(p.dlm, c(list(t.vec=abs(zopt), k=k, side = side), pval.args))
	}
	if(meth.pval == "IS")
	{
		pval <- do.call(p.tube, c(list(t.vec=abs(zopt), k=k, side = side, ncase=ncase
								, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr), pval.args))
	}
	if(meth.pval=="B") pval <- p.bon(abs(zopt), subsetCount, search = 1, side = side)
    

	beta <- sd <- NA
		
	res <- traits.meta(pheno, 1, beta.hat, sigma.hat, ncase, ncntl, rmat=rmat, wt.sigma=TRUE)
	
	beta <- res$beta
	if(side == 2) sd <- abs(beta)/qnorm(pval/2, lower.tail=FALSE, log.p = FALSE)
	else sd <- beta/qnorm(pval, lower.tail=FALSE, log.p = FALSE)
	ret <- list(pval = pval, pheno = pheno, beta = beta, sd = sd, subsetCount=subsetCount)
}

h.traits2 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, side=2, zmax.args=NULL, meth.pval="DLM", pval.args=NULL)
{
	sub1 <- which(beta.hat >= 0)
	sub2 <- which(beta.hat < 0)

	k1 <- length(sub1)
	k2 <- length(sub2)

	pval.args1 <- pval.args
	pval.args2 <- pval.args

	if(!is.null(pval.args) && !is.null(pval.args$sizes))
	{
		sizes <- pval.args$sizes
		csizes <- cumsum(sizes)
		csizes0 <- c(0, cumsum(sizes[-length(sizes)]))
		sizes1 <- vapply(seq_along(sizes), function(i) { sum(sub1 > csizes0[i] & sub1 <= csizes[i]) }, 0L)
		sizes2 <- vapply(seq_along(sizes), function(i) { sum(sub2 > csizes0[i] & sub2 <= csizes[i]) }, 0L)
		sizes1 <- sizes1[sizes1 > 0]
		sizes2 <- sizes2[sizes2 > 0]
		pval.args1$sizes <- sizes1
		pval.args2$sizes <- sizes2
	}

	res1 <- h.traits1(k1, beta.hat=beta.hat[sub1], sigma.hat=sigma.hat[sub1], ncase=ncase[sub1]
					  , ncntl=ncntl[sub1], rmat=rmat[sub1, sub1]
					  , cor.numr=cor.numr, side=side, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args1)

	res2 <- h.traits1(k2, beta.hat=beta.hat[sub2], sigma.hat=sigma.hat[sub2], ncase=ncase[sub2]
					  , ncntl=ncntl[sub2], rmat[sub2, sub2]
					  , cor.numr=cor.numr, side=side, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args2)

	#zopt <- if ( (!is.finite(res1$pval)) || (!is.finite(res2$pval)) || 
    #             (res1$pval <= 0) || (res2$pval <= 0) ) NA else -2 * (log(res1$pval) + log(res2$pval))		
    zopt <- -2*(log(res1$pval) + log(res2$pval))		

    subsetCount <- res1$subsetCount + res2$subsetCount


	if(k1 > 0 && k2 > 0) 
	{ 
		if(meth.pval != "B") pval <- pchisq(zopt, df = 4, lower.tail = FALSE)
		#if(meth.pval == "B") pval <- pchisq(zopt, df = 4, lower.tail = FALSE) * ((2^k1 - 1) * (2^k2 - 1))
        if(meth.pval == "B") pval <- pchisq(zopt, df = 4, lower.tail = FALSE) * subsetCount
	}
	if(k1 == 0 || k2 == 0 )
	{
		if(meth.pval != "B") pval <- pchisq(zopt, df = 2, lower.tail = FALSE)
		#if(meth.pval == "B") pval <- pchisq(zopt, df = 2, lower.tail = FALSE) * ((2^k - 1))
        if(meth.pval == "B") pval <- pchisq(zopt, df = 2, lower.tail = FALSE) * subsetCount
	}

	pheno1 <- rep(FALSE, k) ; pheno1[sub1] <- res1$pheno
	pheno2 <- rep(FALSE, k) ; pheno2[sub2] <- res2$pheno
	
	ret <- list(pval = pval, pval.1=res1$pval, pval.2=res2$pval, pheno.1 = pheno1, pheno.2 = pheno2
				, beta.1 = res1$beta, sd.1 = res1$sd, beta.2 = res2$beta, sd.2 = res2$sd)
	ret
}	

traits.meta <- function(sub, snp.vars, beta.hat, sigma.hat, ncase, ncntl, rmat, side=2, cor.numr=FALSE, wt.sigma=FALSE)
{	
	k <- length(sub)
	if(is.null(dim(rmat))) dim(rmat) <- c(k, k)
	nsnp <- length(snp.vars)
	if(nsnp > 1) stop("Expected nsnp=1")

	z <- beta <- sd <- pval <- rep(NA, nsnp)
	nsub <- sum(sub)
	
	if(!wt.sigma)
	{
		zvec <- beta.hat[sub]/sigma.hat[sub]
		rneff <- sqrt((ncase[sub] * ncntl[sub])/(ncase[sub] + ncntl[sub]))
		if(cor.numr) rneff <- mySolve(rmat[sub, sub], rneff)
		numr <- sum(zvec * rneff)
		denr <- sum(outer(rneff, rneff) * matrix(rmat[sub, sub], nsub, nsub))

        # Changed July 31, 2013
		#z <- ifelse(is.na(denr) | is.nan(denr) | denr == 0, NA, numr/sqrt(denr))
        z <- ifelse(is.na(denr) | is.nan(denr), NA, numr/sqrt(denr))
	} else
	{
		Sigma <- diag(sigma.hat[sub], nsub) %*% matrix(rmat[sub, sub], nsub, nsub) %*% diag(sigma.hat[sub], nsub)
		if(cor.numr) wt <- mySolve(Sigma, rep(1, nsub))
		else wt <- (1/sigma.hat[sub]^2)
		
		numr  <- sum(beta.hat[sub] * wt)
		denr  <- sum(outer(wt, wt) * matrix(Sigma, nsub, nsub))

#		sumwt <- sum(wt)
#		beta  <- numr/sumwt
		beta <- numr/denr		
#		sd    <- sqrt(denr)/sumwt
		sd <- (1/denr)
        # Changed July 31, 2013
		#z <- ifelse(is.na(denr) | is.nan(denr) | denr == 0, NA, numr/sqrt(denr))
        z <- ifelse(is.na(denr) | is.nan(denr), NA, numr/sqrt(denr))
	}
	if(side == 2) pval <- 2 * pnorm(abs(numr/sqrt(denr)), lower.tail=FALSE)
	else pval <- pnorm(numr/sqrt(denr), lower.tail=FALSE)

	list(z = z, beta = beta, sd = sd, pval = pval)
}

traits.forest0 <- function(snp.var, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor = NULL, rlist = NULL
, level = 0.05, p.adj = TRUE, digits = 2)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")	
	if(!(snp.var %in% rownames(beta.hat))) stop("rownames of beta.hat does not have snp.var")
	k <- length(traits.lab)
	
	if(!is.null(rlist))
	{
		if(is.null(rlist$Meta)) stop("Missing Component in rlist: Meta")
		if(is.null(rlist$Subset.1sided)) stop("Missing Component in rlist: Subset.1sided")
		if(is.null(rlist$Subset.2sided)) stop("Missing Component in rlist: Subset.2sided")
	} else
	{
		rlist <- h.traits(snp.vars = snp.var, traits.lab = traits.lab, beta.hat = beta.hat, sigma.hat = sigma.hat
						  , ncase = ncase, ncntl = ncntl, cor = cor, search = NULL, meta = TRUE)
	}
	
	beta <- beta.hat[snp.var, ]
	sd <- sigma.hat[snp.var, ]
	res <- list(beta = beta, sd = sd, pval = 2 * pnorm(abs(beta/sd), lower.tail=FALSE))
	
	h.forest(k, snp.var, traits.lab, rlist, res, side = 2, level=level, p.adj = p.adj, digits=digits)
}

traits.forest <- function(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2, calc.miss=FALSE)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")
	if(!(snp.var %in% rlist$snp.vars)) stop("snp.var was not a SNP analyzed")
    traits.lab <- rlist$traits.lab
    beta.hat   <- rlist$beta.hat[snp.var, ]
    sigma.hat  <- rlist$sigma.hat[snp.var, ]
    ncase      <- rlist$ncase[1, ]
    ncntl      <- rlist$ncntl[1, ]
    cor        <- rlist$cor
    cor.numr   <- rlist$cor.numr
    search     <- rlist$search
    side       <- rlist$side
    zmax.args  <- rlist$zmax.args
    meth.pval  <- rlist$meth.pval
    pval.args  <- rlist$pval.args 

	k <- length(traits.lab)
	nsub <- matrix(FALSE, k, k)
	diag(nsub) <- TRUE
	
    ov      <- rlist[["Meta", exact=TRUE]]
    cc      <- rlist[["Subset.1sided", exact=TRUE]]
    cp      <- rlist[["Subset.2sided", exact=TRUE]]
    ov.flag <- !is.null(ov)
    cc.flag <- !is.null(cc)
    cp.flag <- !is.null(cp)
    if (calc.miss && ((!ov.flag) || (!cc.flag) || (!cp.flag))) 
    {
	if (ov.flag)
	{
		meta <- FALSE
      	} else 
      	{
      		meta <- TRUE
      	}
      	search <- NULL
      	if ((!cc.flag) && (!cp.flag)) 
      	{
      		search <- NULL
      	} else if (!cc.flag) 
      	{
		search <- 1
	} else if (!cp.flag) 
	{
		search <- 2
	}
    	ret <- h.traits(snp.var, traits.lab, beta.hat, sigma.hat, ncase, ncntl, 
                     cor=cor, cor.numr=cor.numr, search=search, side=side, 
                     meta=meta, zmax.args=zmax.args, meth.pval=meth.pval, 
                     pval.args=pval.args)
	if (!ov.flag) ov <- ret[["Meta", exact=TRUE]]
	if (!cc.flag) cc <- ret[["Subset.1sided", exact=TRUE]]
	if (!cp.flag) cp <- ret[["Subset.2sided", exact=TRUE]]
    }
    
    newlist <- list(Meta=ov, Subset.1sided=cc, Subset.2sided=cp)
    
    beta <- beta.hat
    sd   <- sigma.hat
    res  <- list(beta=beta, sd=sd, pval=2*pnorm(abs(beta/sd), lower.tail=FALSE))
	
    h.forest(k, snp.var, traits.lab, newlist, res, side=2, level=level, 
             p.adj=p.adj, digits=digits)

}


