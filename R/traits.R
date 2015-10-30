
h.traits <- function(snp.vars, traits.lab, beta.hat, sigma.hat, ncase=NULL, ncntl=NULL, cor=NULL
	, cor.numr=FALSE, search=NULL, side=2, meta=FALSE, zmax.args=NULL, pval.args=NULL, p.bound=1, 
      NSAMP=5000, NSAMP0=50000)
{

    meth.pval  <- "DLM"
    cond.all   <- TRUE
    add.enrich <- FALSE

    if (is.null(ncase)) ncase <- 1/(sigma.hat^2)
    if (is.null(ncntl)) ncntl <- 1/(sigma.hat^2)

    if (!is.null(search)) {
      if (!(search %in% 1:2)) stop("ERROR: search must be NULL, 1 or 2")
    }
    if ((is.null(side)) || (!(side %in% c(-1, 1, 2)))) stop("ERROR: side must be -1, 1 or 2")
    if ((is.null(meth.pval)) || (!(meth.pval %in% c("DLM", "IS", "B", "Boot")))) stop("ERROR: meth.pval must be DLM, IS, or B")
    if ((is.null(meta)) || (!(meta %in% 0:1))) stop("ERROR: meta must be TRUE or FALSE")
    if ((is.null(cor.numr)) || (!(cor.numr %in% 0:1))) stop("ERROR: cor.numr must be TRUE or FALSE")
    if (!is.null(pval.args)) {
      if (!is.list(pval.args)) stop("ERROR: pval.args must be a named list")
      temp <- names(pval.args)
      if (is.null(temp)) stop("ERROR: pval.args must be a named list")
    }
    if (!is.null(zmax.args)) {
      if (!is.list(zmax.args)) stop("ERROR: zmax.args must be a named list")
      temp <- names(zmax.args)
      if (is.null(temp)) stop("ERROR: zmax.args must be a named list")
    }

	k         <- length(traits.lab)
	meth.pval <- meth.pval[1]
	nsnp      <- length(snp.vars)

	if(p.bound>1 || p.bound <0) stop("Expected a p-value bound between 0 and 1")	

	beta.hat  <- matrix(beta.hat, ncol=k, byrow=FALSE)
	sigma.hat <- matrix(sigma.hat, ncol=k, byrow=FALSE)
	
     if(nrow(beta.hat) > nsnp) stop("ERROR: nrow(beta.hat) must equal length(snp.vars)")	
     if(nrow(sigma.hat) > nsnp) stop("ERROR: nrow(sigma.hat) must equal length(snp.vars)")

	#if(nrow(beta.hat) > nsnp) beta.hat <- matrix(beta.hat[snp.vars, ], nrow = nsnp, ncol=k, byrow=FALSE)
	#if(nrow(sigma.hat) > nsnp) sigma.hat <- matrix(sigma.hat[snp.vars, ], nrow = nsnp, ncol=k, byrow=FALSE)

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
		#if(!all(sapply(cor[c("N11", "N00", "N10", "N01")], function(m) { is.matrix(m) && nrow(m) == k  && ncol(m) == k })))
	    # 		stop("N11, N00, N10 and N01 should be square matrices of dimension k")

       # Changed Aug 08 2013
       if(!all(c("N11", "N00", "N10") %in% names(cor))) {
	     stop("Expected three matrices in list cor with names N11, N00 and N10")
       }
	   if(!all(sapply(cor[c("N11", "N00", "N10")], function(m) { is.matrix(m) && nrow(m) == k  && ncol(m) == k })))
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
		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
				nsub <- sum(sub)
				rmat <- matrix(cor[sub, sub], nsub, nsub)
				res <- try(traits.meta(sub, snp.vars[j], beta.hat[j, ], sigma.hat[j,]
							, ncase[jj, ], ncntl[jj, ], rmat = cor, side = side
							, cor.numr = cor.numr, wt.sigma = TRUE), silent=TRUE)
                     if ("try-error" %in% class(res)) next
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
          pval <- pval1 <- pval2 <- rep(NA, nsnp)
	     names(pval) <- names(pval1) <- names(pval2) <- snp.vars
	     beta <- beta1 <- beta2 <- rep(NA, nsnp)
	     names(beta) <- names(beta1) <- names(beta2) <- snp.vars
	     sd <- sd1 <- sd2 <- rep(NA, nsnp)
	     names(sd) <- names(sd1) <- names(sd2) <- snp.vars
	     pheno <- pheno1 <- pheno2 <- matrix(FALSE, nsnp, k)
	     colnames(pheno) <- colnames(pheno1) <- colnames(pheno2) <- traits.lab
	     rownames(pheno) <- rownames(pheno1) <- rownames(pheno2) <- snp.vars

		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 1))
		else pval.args$search <- 1

		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
				nsub <- sum(sub)
				rmat <- matrix(cor[sub, sub], nsub, nsub)
				res <- try(h.traits1(nsub, beta.hat[j, sub], sigma.hat[j, sub], ncase[jj, sub], ncntl[jj, sub]
							, rmat = rmat, cor.numr=cor.numr, search=1, side=side, zmax.args=zmax.args
						 	, meth.pval=meth.pval, pval.args=pval.args, p.bound=p.bound, add.enrich=add.enrich,
                                       cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0), silent=TRUE)
                     if ("try-error" %in% class(res)) next
                     res.pheno <- res[["pheno", exact=TRUE]]
                     PFLAG     <- !is.null(res.pheno)

				pval[j]       <- res$pval
				if (PFLAG) pheno[j, sub] <- res.pheno
				beta[j]       <- res$beta
				sd[j]         <- res$sd

				# Compute new meta-analysis standard error
				if (PFLAG) sd.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res.pheno)

			}
			
		}
		sub1.res <- list(pval=pval, beta=beta, sd=sd, pheno=pheno, sd.meta=sd.meta)
	}

	t00 <- proc.time()
	if(TS)
	{ 
           pval <- pval1 <- pval2 <- rep(NA, nsnp)
	     names(pval) <- names(pval1) <- names(pval2) <- snp.vars
	     beta <- beta1 <- beta2 <- rep(NA, nsnp)
	     names(beta) <- names(beta1) <- names(beta2) <- snp.vars
	     sd <- sd1 <- sd2 <- rep(NA, nsnp)
	     names(sd) <- names(sd1) <- names(sd2) <- snp.vars
	     pheno <- pheno1 <- pheno2 <- matrix(FALSE, nsnp, k)
	     colnames(pheno) <- colnames(pheno1) <- colnames(pheno2) <- traits.lab
	     rownames(pheno) <- rownames(pheno1) <- rownames(pheno2) <- snp.vars

		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 2))
		else pval.args$search <- 2

		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			res <- NULL
			sub <- !(is.na(beta.hat[j, ]) | is.na(sigma.hat[j, ]) | is.na(ncase[jj, ]) | is.na(ncntl[jj, ]))
			if(any(sub))
			{
			    nsub <- sum(sub)
			    rmat <- matrix(cor[sub, sub], nsub, nsub)

			    res <- try(h.traits2(nsub, beta.hat[j, sub], sigma.hat[j, sub], ncase[jj, sub], ncntl[jj, sub], 
                                 rmat=rmat, cor.numr=cor.numr, side = side, zmax.args=zmax.args, 
                                 meth.pval=meth.pval, pval.args=pval.args, p.bound=p.bound, add.enrich=add.enrich,
                                 cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0), silent=TRUE)
                    if ("try-error" %in% class(res)) next 
                    res.pheno.1 <- res[["pheno.1", exact=TRUE]]
                    PFLAG1      <- !is.null(res.pheno.1)
                    res.pheno.2 <- res[["pheno.2", exact=TRUE]]
                    PFLAG2      <- !is.null(res.pheno.2)

                    pval[j]        <- res$pval
			    pval1[j]       <- res$pval.1
		     	    pval2[j]       <- res$pval.2
			    if (PFLAG1) pheno1[j, sub] <- res.pheno.1
			    if (PFLAG2) pheno2[j, sub] <- res.pheno.2
			    beta1[j]       <- res$beta.1
			    beta2[j]       <- res$beta.2
			    sd1[j]         <- res$sd.1
			    sd2[j]         <- res$sd.2

                  # Compute new meta-analysis standard error
                  if (PFLAG1) sd1.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res.pheno.1)
                  if (PFLAG2) sd2.meta[j] <- meta.se(cor, sigma.hat[j, ], subset=res.pheno.2)

			}
			
			
		}

		sub2.res <- list(pval=pval, pval.1 = pval1, beta.1=beta1, sd.1=sd1, pheno.1=pheno1,
						 pval.2=pval2, beta.2=beta2, sd.2=sd2, pheno.2=pheno2,
                         sd.1.meta=sd1.meta, sd.2.meta=sd2.meta)
	}

    ret <- list(Meta=meta.res, Subset.1sided=sub1.res, Subset.2sided=sub2.res,
           snp.vars=snp.vars, traits.lab=traits.lab, beta.hat=beta.hat, sigma.hat=sigma.hat,
           ncase=ncase, ncntl=ncntl, cor=cor, cor.numr=cor.numr, search=search,
           side=side, meta=meta, zmax.args=zmax.args, meth.pval=meth.pval,
           pval.args=pval.args, which="h.traits", cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)
 
	ret
}


h.traits1 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, search=1, side=2, zmax.args=NULL, 
                meth.pval="DLM", pval.args=NULL, p.bound=1, add.enrich=FALSE,
                cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{
#print("BEGIN h.traits1")
	if(k == 0) return(list(pval = 1, pheno = NULL, beta = NA, sd = NA, subsetCount=0))
	if (length(rmat) == 1) dim(rmat) <- c(1, 1)

	z0 <- beta.hat/sigma.hat
	pvals0 <- if(search < 2 && abs(side) == 1) pnorm(abs(z0), lower.tail=FALSE) else (2 * pnorm(abs(z0), lower.tail=FALSE))

	z.sub <- (pvals0 <= p.bound) * sign(z0)
	k1 <- sum(z.sub != 0)
	if(k1 == 0) return(list(pval = 1, pheno = NULL, beta = NA, sd = NA, subsetCount=0))
 
     # Determine if the user passed in subset functions
     ###############################################
     sub.def   <- zmax.args[["sub.def", exact=TRUE]] 
     sub.args  <- zmax.args[["sub.args", exact=TRUE]]
     psub.def  <- pval.args[["sub.def", exact=TRUE]] 
     psub.args <- pval.args[["sub.args", exact=TRUE]]
     ###############################################

     # For zmax
     if (!is.null(sub.def)) {
	  sub.def1 <- function(x, op)
	  {
           z.sub <- op$z.sub
		x0    <- (z.sub!=0) & x
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
	  if(k1 != k)
	  {
           pval.args$sub.def  <- psub.def1
           pval.args$sub.args <- c(psub.args, list(z.sub=z.sub))
	  }
     }


	if((k1 != k) && !is.null(pval.args$sizes))
	{
		sizes <- pval.args$sizes
		labs <- rep(1:length(sizes), sizes)
		sizes1 <- NULL
		for(i in 1:length(sizes))
		{
			num.skip <- sum((z.sub==0) & labs == i)
			sizes1 <- c(sizes1, c(sizes[i] - num.skip))
		}
		pval.args$sizes <- sizes1
	}

	if(is.null(zmax.args) || is.null(zmax.args$z.sub)) zmax.args <- c(zmax.args, list(z.sub = z.sub))
	if(is.null(pval.args) || is.null(pval.args$p.bound)) pval.args <- c(pval.args, list(p.bound = p.bound))

##### NOTE: meta.def is set to traits.meta ######
#print("BEGIN z.max")
	res <- do.call(z.max, c(list(k, 1, side=side, meta.def=traits.meta, meta.args=list(beta.hat=beta.hat, sigma.hat=sigma.hat
				 , ncase=ncase, ncntl=ncntl, rmat=rmat, side=side, cor.numr=cor.numr, wt.sigma=FALSE)), zmax.args))
#print("END z.max")
	zopt <- as.double(res$opt.z)
	pheno <- as.logical(res$opt.s)
	subsetCount <- res$subsetCount
	rm(res)

	if(meth.pval == "DLM")
	{
		if(is.null(pval.args) || !("cor.def" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.def=NULL))
		if(!("cor.args" %in% names(pval.args))) {
             #pval.args <- c(pval.args, 
             #  list(cor.args=list(ncase=ncase[(z.sub!=0)], ncntl=ncntl[(z.sub != 0)], rmat=rmat[(z.sub != 0), (z.sub != 0)], cor.numr=cor.numr)))
             pval.args <- c(pval.args, list(cor.args=list(ncase=ncase, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr)))
           }
		#pval <- do.call(p.dlm, c(list(t.vec=abs(zopt), z.sub=z.sub, search=search, side = abs(side)), pval.args))      
#print("BEGIN p.dlm")    
           pval <- p.dlm(abs(zopt), z.sub, search, abs(side), cor.def=pval.args[["cor.def", exact=TRUE]], 
                  cor.args=pval.args[["cor.args", exact=TRUE]], sizes=pval.args[["sizes", exact=TRUE]], p.bound=p.bound, 
                  sub.def=pval.args[["sub.def", exact=TRUE]], sub.args=pval.args[["sub.args", exact=TRUE]],
                  NSAMP=NSAMP, NSAMP0=NSAMP0)
#print("END p.dlm")
	}
	if(meth.pval == "IS")
	{
		pval <- do.call(p.tube, c(list(t.vec=abs(zopt), k=k1, side = abs(side), ncase=ncase
					, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr, cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0), pval.args))
	}
	if(meth.pval == "Boot")
	{
		pval <- p.boot(t0=abs(zopt), z.sub=z.sub, search=search, side=abs(side), ncase0=ncase, ncntl0=ncntl, rmat0=rmat, 
                   cor.numr=cor.numr, p.bound=p.bound, cond.all=cond.all, NSAMP0=NSAMP0)
	}
	if(meth.pval=="B") pval <- p.bon(abs(zopt), subsetCount, search = 1, side = abs(side))


	beta <- sd <- NA
		
       if (any(pheno)) 
       {

	  res <- traits.meta(pheno, 1, beta.hat, sigma.hat, ncase, ncntl, rmat=rmat, side=side, cor.numr=cor.numr, wt.sigma=TRUE)

	  beta <- res$beta
	  if(side == 2) sd <- abs(beta)/qnorm(pval/2, lower.tail=FALSE, log.p = FALSE)
	  else sd <- abs(beta/qnorm(pval, lower.tail= (side == -1), log.p = FALSE))
       }
	if(add.enrich && p.bound < 1)
	{
		pval1 <- calcP1(p.bound, k1, k, search=search, side=side, rmat=rmat, sizes=sizes)
		pval.c  <- pchisq(-2 * (log(pval)+ log(pval1)), df=4, lower.tail=FALSE)
		pval <- pval.c
	}
#print("END h.traits1")
	list(pval = pval, pheno = pheno, beta = beta, sd = sd, subsetCount=subsetCount)
}

h.traits2 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, side=2, zmax.args=NULL, meth.pval="DLM", 
                      pval.args=NULL, p.bound=1, add.enrich=FALSE, cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{
#print("BEGIN h.traits2")
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
		sizes1 <- sapply(1:length(sizes), function(i) { sum(sub1 > csizes0[i] & sub1 <= csizes[i]) })
		sizes2 <- sapply(1:length(sizes), function(i) { sum(sub2 > csizes0[i] & sub2 <= csizes[i]) })
		sizes1 <- sizes1[sizes1 > 0]
		sizes2 <- sizes2[sizes2 > 0]
		pval.args1$sizes <- sizes1
		pval.args2$sizes <- sizes2
	}

	res1 <- h.traits1(k1, beta.hat=beta.hat[sub1], sigma.hat=sigma.hat[sub1], ncase=ncase[sub1]
					  , ncntl=ncntl[sub1], rmat=rmat[sub1, sub1]
					  , cor.numr=cor.numr, search=2, side=1, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args1, p.bound=p.bound, add.enrich=add.enrich,
                              cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)

	res2 <- h.traits1(k2, beta.hat=beta.hat[sub2], sigma.hat=sigma.hat[sub2], ncase=ncase[sub2]
					  , ncntl=ncntl[sub2], rmat[sub2, sub2]
					  , cor.numr=cor.numr, search=2, side=-1, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args2, p.bound=p.bound, add.enrich=add.enrich,
                               cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)

	#zopt <- if ( (!is.finite(res1$pval)) || (!is.finite(res2$pval)) || 
    #             (res1$pval <= 0) || (res2$pval <= 0) ) NA else -2 * (log(res1$pval) + log(res2$pval))		
    zopt <- -2*(log(res1$pval) + log(res2$pval))		

    subsetCount <- res1$subsetCount + res2$subsetCount


       if(k1 > 0 && k2 > 0) 
	{ 
		pval <- pchisq(zopt, df = 4, lower.tail = FALSE)
	}
	if(k1 == 0 || k2 == 0 )
	{
		pval <- pchisq(zopt, df = 2, lower.tail = FALSE)
	}

	pheno1 <- rep(FALSE, k)
      if (!is.null(res1$pheno)) pheno1[sub1] <- res1$pheno
	pheno2 <- rep(FALSE, k)
      if (!is.null(res2$pheno)) pheno2[sub2] <- res2$pheno
	
	ret <- list(pval = pval, pval.1=res1$pval, pval.2=res2$pval, pheno.1 = pheno1, pheno.2 = pheno2
				, beta.1 = res1$beta, sd.1 = res1$sd, beta.2 = res2$beta, sd.2 = res2$sd)

#print("END h.traits2")

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
	if(side == 2) pval <- 2 * pnorm(abs(z), lower.tail=FALSE)
	if(side == 1) pval <- pnorm(z, lower.tail=FALSE)
	if(side == -1) pval <- pnorm(z, lower.tail=TRUE)
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

traits.forest <- function(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")
	if(!(snp.var %in% rlist$snp.vars)) stop("snp.var was not a SNP analyzed")
    traits.lab <- rlist$traits.lab
    beta.hat   <- rlist$beta.hat
    sigma.hat  <- rlist$sigma.hat
    ncase      <- rlist$ncase
    ncntl      <- rlist$ncntl
    cor        <- rlist$cor
    cor.numr   <- rlist$cor.numr
    search     <- rlist$search
    side       <- rlist$side
    zmax.args  <- rlist$zmax.args
    meth.pval  <- rlist$meth.pval
    pval.args  <- rlist$pval.args 
    cond.all   <- rlist$cond.all
    NSAMP      <- rlist$NSAMP
    NSAMP0     <- rlist$NSAMP0

	k <- length(traits.lab)
	nsub <- matrix(FALSE, k, k)
	diag(nsub) <- TRUE
	
    ov      <- rlist[["Meta", exact=TRUE]]
    cc      <- rlist[["Subset.1sided", exact=TRUE]]
    cp      <- rlist[["Subset.2sided", exact=TRUE]]
    ov.flag <- !is.null(ov)
    cc.flag <- !is.null(cc)
    cp.flag <- !is.null(cp)
    if ((!ov.flag) || (!cc.flag) || (!cp.flag)) {
      if (ov.flag) {
        meta <- FALSE
      } else {
        meta <- TRUE
      }

      search <- NULL
      if ((!cc.flag) && (!cp.flag)) {
        search <- NULL
      } else if (!cc.flag) {
        search <- 1
      } else if (!cp.flag) {
        search <- 2
      }
      ret <- h.traits(snp.var, traits.lab, beta.hat, sigma.hat, ncase, ncntl, 
                      cor=cor, cor.numr=cor.numr, search=search, side=side, 
                      meta=meta, zmax.args=zmax.args, 
                      pval.args=pval.args, NSAMP=NSAMP, NSAMP0=NSAMP0)

      if (!ov.flag) ov <- ret[["Meta", exact=TRUE]]
      if (!cc.flag) cc <- ret[["Subset.1sided", exact=TRUE]]
      if (!cp.flag) cp <- ret[["Subset.2sided", exact=TRUE]]
    }
    
    newlist <- list(Meta=ov, Subset.1sided=cc, Subset.2sided=cp)
    
	beta <- beta.hat[snp.var, ]
	sd   <- sigma.hat[snp.var, ]
	res  <- list(beta=beta, sd=sd, pval=2*pnorm(abs(beta/sd), lower.tail=FALSE))
	
	h.forest(k, snp.var, traits.lab, newlist, res, side=2, level=level, 
             p.adj=p.adj, digits=digits)

}


