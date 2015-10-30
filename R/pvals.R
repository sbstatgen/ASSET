
cargs.proc <- function(search, side, ncase0, ncntl0, rmat0=NULL, cor.numr=FALSE, pool=FALSE)
{
     if (is.null(cor.numr)) cor.numr <- FALSE
	k0 <- length(ncase0)
	srch <- search
	cnumr <- cor.numr

	if(search == 1 || search == 2)
	{
		if(is.null(rmat0)) rmat0 <- diag(1, k0)
		if((length(ncase0) != k0) || (length(ncntl0) != k0)) stop("Length of ncase and ncntl should be k")
		neff0 <- (ncase0 * ncntl0)/(ncase0 + ncntl0)
	}
	if(search == 0)
	{
		if((length(ncase0) != k0) || (length(ncntl0) != 1)) stop("Length of ncase should be k and Length of ncntl should be 1")
		N0 <- sum(ncase0) + ncntl0
		ncntl1 <- rep(ncntl0, k0)
		if(pool == FALSE)
		{
			neff0       <- (ncase0 * ncntl1)/(ncase0 + ncntl1)
                vec         <- sqrt(ncase0/(ncase0 + ncntl0))
                dim(vec)    <- NULL
			rmat0       <- outer(vec, vec)
			diag(rmat0) <- 1
		}
		if(pool == TRUE)
		{
			N1          <- ncase0 + ncntl1
			neff0       <- (ncase0 * (N1 - ncase0))/N1
                vec         <- ncase0/sqrt(neff0)
                dim(vec)    <- NULL
			rmat0       <- (-1)*outer(vec, vec)/N0
			diag(rmat0) <- 1
		}
	
#		Note: hard coded
		cnumr <- FALSE
#		srch <- 1
	}

	list(neff0=neff0, rmat0=rmat0, cor.numr=cnumr, search=srch)
}


p.boot <- function(t0, z.sub, search, side, ncase0, ncntl0, rmat0=NULL, cor.numr=FALSE, pool=FALSE, 
             p.bound=1, cond.all=TRUE, NSAMP0=5e4)
{
	k <- sum(z.sub != 0)
	k0 <- length(z.sub)
	Z0 <- if(search < 2 && side == 1) qnorm(p.bound, lower.tail=FALSE) else qnorm(p.bound/2, lower.tail=FALSE)

	ret <- cargs.proc(search, side, ncase0, ncntl0, rmat0=rmat0, cor.numr=cor.numr, pool=pool)
	rmat0 <- ret$rmat0 ; neff0 <- ret$neff0 ; cor.numr <- ret$cor.numr ; search <- ret$search
	rmat <- matrix(rmat0[(z.sub != 0), (z.sub != 0)], k, k)
	ncase <- ncase0[z.sub != 0]
	ncntl <- if(search == 0) ncntl0 else ncntl0[z.sub != 0]
	neff <- neff0[z.sub != 0]

	low <- rep(-Inf, k0)
	high <- rep(Inf, k0)
	means <- rep(0, k0)
	imod <- rep(FALSE, k0)
	if(side==2)
	{
		if(cond.all) imod[1:k0] <- TRUE
		else imod[which(z.sub != 0)] <- TRUE
	}		
	if(k > 0) low[which(z.sub != 0)] <- Z0
	if(cond.all)
	{
		if(k < k0) high[which(z.sub == 0)] <- Z0
	}

	res <- tmvnsim(nsamp=NSAMP0, k=k0, lower=low, upper=high, imod=imod, means=means, sigma=rmat0)
	zz <- res$samp
	wts <- res$wts
	num <- 0
	for(i in 1:NSAMP0)
	{
		res <- z.max(k0, 1, side=side, meta.def=traits.meta, meta.args=list(beta.hat=zz[i, ], sigma.hat = rep(1, k0)
				 , ncase=neff0, ncntl=neff0, rmat=rmat0, side=side, cor.numr=cor.numr, wt.sigma=FALSE)
				, th = NULL, z.sub = z.sub)
		if(side == 1 && res$opt.z > t0) num <- num + wts[i]
		if(side == 2 && abs(res$opt.z) > t0) num <- num + wts[i]
	}
	ans <- num/sum(wts)
	ans
}

p.bon <- function(t.vec, M, search, side)
{
	# M    number of subsets considered

	if(any(t.vec < 0)) warning("Negative input thresholds! Using absolute value.")
	t.vec <- abs(t.vec)
	if(search < 2 & side == 2) pval.b <- pmin(2 * pnorm(t.vec, lower.tail=FALSE) * M, 1)
	if(search < 2 & side == 1) pval.b <- pmin(pnorm(t.vec, lower.tail=FALSE) * M, 1)
	if(search == 2) stop("Seach option 2 not yet implemented")

	pval.b
}

p.tube <- function(t.vec, z.sub, search, side, ncase, ncntl, rmat=NULL, cor.numr=FALSE, pool=FALSE, sizes=rep(1, k)
	, p.bound=1, nsamp=50, sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL,
       cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{
	k <- sum(z.sub != 0)
	k0 <- length(z.sub)
	if(!is.null(sub.def) && !is.function(sub.def)) stop("'sub.def' not recognized")
	
	if(k <= 0 || (k%%1 != 0)) stop("Number of studies k, not a positive integer!")
	if(!search %in% c(0, 1, 2)) stop("Invalid search option")
	if(search < 2 && !(side %in% c(1, 2))) stop("side should be 1 or 2")
	if(any(t.vec < 0)) warning("Negative input thresholds! Using absolute value.")

	nsnp <- length(t.vec)
	if(!is.null(dim(ncase)) && ncol(ncase) != 1 && ncol(ncase) != nsnp) stop("Unexpected dim of ncase")
	if(!is.null(dim(ncntl)) && ncol(ncntl) != 1 && ncol(ncntl) != nsnp) stop("Unexpected dim of ncntl")
	
	tube.pval(t.vec, z.sub=z.sub, search=search, side=side, ncase0=ncase, ncntl0=ncntl
			  , rmat0=rmat, cor.numr=cor.numr, pool=pool, sizes=sizes, p.bound=p.bound, nsamp=nsamp
			  , sub.def=sub.def, sub.args=sub.args, cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)
}


p.dlm <- function(t.vec, z.sub, search, side, cor.def = NULL, cor.args=NULL, sizes=NULL, p.bound=1, sub.def=NULL, sub.args=NULL
		, NSAMP=5000, NSAMP0=5e4)
{
     wt.def   <- NULL
     wt.args  <- NULL
     cond.all <- TRUE     

	k <- sum(z.sub != 0)
	k0 <- length(z.sub)
	if(!is.null(sub.def) && !is.function(sub.def)) stop("'sub.def' not recognized")
     if (is.null(sizes)) sizes <- rep(1, k)

	
	if(k <= 0 || (k%%1 != 0)) stop("Number of studies k, not a positive integer!")
	if(!search %in% c(0, 1, 2)) stop("Invalid search option")
	if(search < 2 && !(side %in% c(1, 2))) stop("side should be 1 or 2")
	if(search == 2) side <- 1

	if(any(t.vec < 0)) { warning("Negative input thresholds! Using absolute value.") ; t.vec <- abs(t.vec) ; }
	dlm.pval(t.vec=t.vec, z.sub=z.sub, search=search, side=side, cor.def=cor.def, cor.args=cor.args, sizes=sizes, 
              p.bound=p.bound, sub.def=sub.def, sub.args=sub.args,
              cond.all=cond.all, NSAMP=NSAMP, NSAMP0=NSAMP0)
}

calcP1 <- function(p.bound, k1, k, search=1, side=2, rmat=diag(1, k), sizes=rep(1, k))
{
	pbinom(k1, k, prob=p.bound, lower.tail=FALSE)	
}

prob.sim <- function(x10, X20, z.sub, t0, nsamp, rmat0, rneff0, search=1, side=2, cor.numr=FALSE, 
            p.bound=1, Prob0=1, cond.all=TRUE)
{
#print("BEGIN prob.sim")
	k0 <- length(x10)
	k <- sum(z.sub != 0)
	nn <- ncol(X20)

	Z0 <- if(search < 2 && side == 1) qnorm(p.bound, lower.tail=FALSE) else qnorm(p.bound/2, lower.tail=FALSE)

	x10 <- as.logical(x10)
	if(cor.numr)
	{
		n1 <- sum(x10)
		rneff1 <- rep(0, k0)
		rneff1[x10] <- drop(mySolve(matrix(rmat0[x10, x10], n1, n1) , rneff0[x10]))
	} else rneff1 <- (rneff0 * x10)

	denr1 <- drop(matrix(rneff1, nrow=1) %*% rmat0 %*% matrix(rneff1, ncol=1))
	w <- rneff1/ sqrt(denr1) 

	dim(w) <- c(k0,1)
	sigma <- matrix(0, k0 + 1, k0 + 1)
	rw <- drop(rmat0 %*% w)
	sigma[(1:k0)+1, (1:k0)+1] <- rmat0
	sigma[(1:k0)+1, 1] <- rw
	sigma[1, (1:k0)+1] <- rw
	sigma[1, 1] <- sum(w * rw)

	means <- rep(0, k0 + 1)
	low <- rep(-Inf, k0+1)
	high <- rep(Inf, k0+1)
	imod <- rep(FALSE, k0 + 1)
	low[1] <- t0
	if(side==2)
	{
		if(cond.all) imod[c(1, 1 + (1:k0))] <- TRUE
		else imod[c(1, 1 + which(z.sub != 0))] <- TRUE
	}		
	if(k > 0) low[1 + which(z.sub != 0)] <- Z0
	if(cond.all)
	{
		if(k < k0) high[1 + which(z.sub == 0)] <- Z0
	}
	if(any(is.na(low)) || any(is.na(high)) || any(is.nan(low)) || any(is.nan(high))) stop("NA/NaN !!")

#print(z.sub)
#save(nsamp, k0, low, high, imod, means, sigma, file="/spin1/users/wheelerwi/nilanjan/wga/ASSET/simulation/oct14_2014/data.rda")

	res <- tmvnsim(nsamp=nsamp, k=(1+k0), lower=low, upper=high, imod=imod, means=means, sigma=sigma)
	Z <- res$samp
	wts <- res$wts
	Prob1 <- mean(wts)

	z.func <- function(x20)
	{	
		x20 <- as.logical(x20)
#		p20 <- which(x20)
		k2 <- sum(x20)
		rneff2 <- rep(0, k0)
		if(cor.numr) rneff2[x20] <- mySolve(matrix(rmat0[x20, x20], k2, k2), rneff0[x20])
		else rneff2 <- (rneff0 * x20)
		denr2 <- drop(matrix(rneff2, nrow=1) %*% rmat0 %*% matrix(rneff2, ncol=1))
		rr2 <- rneff2/ sqrt(denr2) 
#		z <- drop(matrix(Z[, 1+p20], nsamp, k2) %*% matrix(rr2[x20], ncol = 1))
		z <- drop(matrix(Z[, -1], nsamp, k0) %*% matrix(rr2, ncol = 1))
		if(any(is.na(z))) stop("NA")

		z
	}
     if (nn) {
 	  zvals  <- matrix(apply(X20, 2, z.func), nsamp, nn)
	  zval0  <- Z[, 1]
	  check1 <- sapply(1:nsamp, 
                 function(ii){ if(side == 2) all(abs(zvals[ii, ]) <= abs(zval0[ii])) else all(zvals[ii, ] <= zval0[ii])})
     } else{
       check1 <- rep(1, nsamp)
     }
	prob  <- sum(check1 * wts)/sum(wts)
	tprob <- min(Prob1/Prob0, 1)
	ret   <- (prob * tprob)
	if(is.na(ret) || is.nan(ret)) stop("NA/NaN !!") 

#print("END prob.sim")

	ret
}

calcP0 <- function(Z0, z.sub, search, side, p.bound, rmat=diag(1, k), cond.all=TRUE, NSAMP0=5e4)
{
	k0 <- length(z.sub)
	k <- sum(z.sub != 0)

	means <- rep(0, k0)
	low <- rep(-Inf, k0)
	high <- rep(Inf, k0)
	imod <- rep(FALSE, k0)
	if(side==2)
	{
		if(cond.all) imod[1:k0] <- TRUE
		else imod[which(z.sub != 0)] <- TRUE
	}		
	if(k > 0) low[which(z.sub != 0)] <- Z0
	if(cond.all)
	{
		if(k < k0) high[which(z.sub == 0)] <- Z0
	}

	res <- tmvnsim(nsamp=NSAMP0, k=k0, lower=low, upper=high, imod=imod, means=means, sigma=rmat)
	Prob0 <- mean(res$wts)

	Prob0
}


dlm.pval <- function(t.vec, z.sub, search, side, cor.def, cor.args, sizes = NULL, p.bound=1
				, sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL, 
                       cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{

	NT <- length(t.vec)
	k0 <- length(z.sub)
	k <- sum(z.sub != 0)
     if (is.null(sizes)) sizes <- rep(1, k)
	
	NS <- (prod(sizes + 1) - 1)	
	if(any(t.vec < 0)) t.vec <- abs(t.vec)

	low <- t.vec
	high <- rep(Inf, NT)
		
	r.vec2 <- r.vec3 <- NULL

	Z0 <- if(search < 2 && side == 1) qnorm(p.bound, lower.tail=FALSE) else qnorm(p.bound/2, lower.tail=FALSE)
	rmat0 <- cor.args$rmat
	if(is.null(rmat0)) rmat0 <- diag(k0)
	ncase0 <- cor.args$ncase
	ncntl0 <- cor.args$ncntl
	cor.numr <- cor.args$cor.numr
	pool <- cor.args$pool


	ret <- cargs.proc(search, side, ncase0, ncntl0, rmat0=rmat0, cor.numr=cor.numr, pool=pool)
	rmat0 <- ret$rmat0 ; neff0 <- ret$neff0 ; cor.numr <- ret$cor.numr ; search <- ret$search

	rmat <- matrix(rmat0[(z.sub != 0), (z.sub != 0)], k, k)
	ncase <- ncase0[z.sub != 0]
	ncntl <- if(search == 0) ncntl0 else ncntl0[z.sub != 0]
	neff <- neff0[z.sub != 0]

	Prob0 <- calcP0(Z0=Z0, z.sub=z.sub, search=search, side=side, p.bound=p.bound, rmat=rmat0, 
                     cond.all=cond.all, NSAMP0=NSAMP0)

	x <- rep(0, k)
	kk <- length(sizes)
	xx <- rep(0, kk)
	cc <- c(0, cumsum(sizes[-kk]))
	
	ss <- rep(0, NT)
	i <- 1
#print("BEGIN while loop")
	while (i <= NS)
	{

		pos <- max(which(xx < sizes))
		if(pos < kk)
		{
			xx[(pos + 1):kk] <- 0
			x[(cc[pos + 1] + 1):k] <- 0
		}
		xx[pos] <- xx[pos] + 1
#print(cc[pos])
#print(xx[pos])
#print(kk)
		x[cc[pos] + xx[pos]] <- 1
		set <- as.logical(x)
		x0 <- rep(0, k0)
#print(c(k0, k, length(x0), length(x), sum(z.sub != 0)))
#print(c(pos, kk))

		x0[z.sub != 0] <- x
		NXX <- prod(choose(sizes, xx))
#		DEBUG
#		print(c(i, pos, cc[pos], xx, NXX))
		

#		Skip this subset
		#if(!is.null(sub.def) && !do.call(sub.def, c(list(set), sub.args))) { i <- i + 1 ; next }
		if (!is.null(sub.def)) {
			if (!sub.def(set, sub.args)) { i <- i + 1 ; next }
		}
#		Matrix with neighboring subsets of x as columns
		nn.x <- (matrix(x, k, k) + diag(1, k)) %% 2

		nn.sub <- rep(TRUE, k)
		
#		Note: All zeros is not a valid neighbor
 		if(sum(x) == 1)
		{
			nn.x <- nn.x[, (x == 0)]
			if(is.null(dim(nn.x))) nn.x <- matrix(nn.x, nrow = k)
			nn.sub <- rep(TRUE, k - 1)
		}
		nn.x0 <- matrix(0, k0, ncol(nn.x))

		nn.x0[z.sub != 0, ] <- nn.x
#		Valid neighbors
		if(!is.null(sub.def)) 
		{
		  nn.set0 <- as.logical(nn.x0)
		  dim(nn.set0) <- dim(nn.x0)
  		  # nn.sub <- nn.sub & apply(nn.set, 2, function(set1) do.call(sub.def, c(list(set1), sub.args)))
		  nn.sub <- nn.sub & apply(nn.set0, 2, function(set1) sub.def(set1, sub.args) )
		}
#		Add a term to points with integral still below 1
		pos <- which(ss < 1)
		if(length(pos) == 0) break
		t.sub <- (1:NT)[pos]
#print("BEGIN sapply")

		int <- sapply(t.sub, function(j)
					{
						prob.sim(x0, matrix(nn.x0[, nn.sub], nrow=k0, byrow=FALSE), z.sub, 
                                        low[j], NSAMP, rmat0, sqrt(neff0), search=search, 
                                         side=side, cor.numr=cor.numr, p.bound=p.bound, Prob0=Prob0, cond.all=cond.all)
					})
#print("END sapply")
		ss[t.sub] <- ss[t.sub] + NXX * int
		i <- i + 1
	}
#print("END while loop")

	ret <- pmin(1, ss)

#print("END dlm.pval")

	ret
}


tubeSim <- function(w, a, b, inner, nsamp, rmat)
{
	k <- length(w)

	Z <- matrix(0, c(k, nsamp))
	r1 <- (which(w != 0))[1]
	r2 <- (1:k)[-r1]

	A <- diag(k)
	A[r1, ] <- w
	A <- A[c(r1, r2), ]

	Sigma <- A %*% rmat %*% t(A)
	Sig <- Sigma[(2:k), (2:k)] -  ((1/Sigma[1, 1]) * outer(Sigma[(2:k), 1], Sigma[(2:k), 1]))	
	ev <- eigen(Sig, symmetric = TRUE)

	if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1])))
	{
		warning("Sigma is numerically not positive definite")
	}
	Sigroot <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
	if(inner) {
		W1 <- rtnorm(nsamp, mean = 0, sd = sqrt(Sigma[1, 1]), lower = a, upper=b)
	} else {
		p1 <- pnorm(b, lower.tail=FALSE)
		p2 <- pnorm(a)
		del <- rbinom(nsamp, 1, prob=(p1/(p1+p2)))
		#del <- rep(nsamp, 1)
		#del[1:floor(nsamp * (p1/(p1 + p2))] <- 1
		W11 <- rtnorm(nsamp, mean = 0, sd = sqrt(Sigma[1, 1]), lower = b)
		W12 <- rtnorm(nsamp, mean = 0, sd = sqrt(Sigma[1, 1]), upper = a)
		W1 <- ifelse(del == 1, W11, W12)
	}
	mu <- outer(Sigma[(2:k), 1]/Sigma[1, 1], W1)
		
	W2 <- t(matrix(rnorm(nsamp * (k - 1)), nrow = nsamp) %*% Sigroot) + mu				
	
	W <- rbind(W1, W2)
	Z <- solve(A, W)
	Z
}


checkSim <- function(Z, w.mat, a.vec, b.vec, inner.vec, rmat)
{
	k <- nrow(Z)

	nsamp <- ncol(Z)
	con <- length(b.vec)
	a.mat <- matrix(a.vec, con, nsamp, byrow=FALSE)
	b.mat <- matrix(b.vec, con, nsamp, byrow=FALSE)
	inner.mat <- matrix(inner.vec, con, nsamp, byrow=FALSE)

	wZ <-  w.mat %*% Z
	f1 <- (wZ > a.mat)
	f2 <- (wZ < b.mat)
	fm <- colSums(ifelse(inner.mat, f1 & f2, ((!f1) | (!f2))))

	gm <- mean(ifelse(fm > 0, 1/fm, 0))

	gm
}



tube.pval <- function(t0, z.sub, search, side, ncase0, ncntl0, rmat0 = NULL
			, cor.numr = FALSE, pool=FALSE, sizes = rep(1, k), p.bound = 1, nsamp = 50, sub.def=NULL, 
                  sub.args=NULL, cond.all=TRUE, NSAMP=5000, NSAMP0=5e4)
{
	if(length(t0) != 1) stop("t0 must have length 1")
	t0 <- abs(t0)

	k <- sum(z.sub != 0)
	k0 <- length(z.sub)

#	if(k == 1) return(side * pnorm(t.vec, lower.tail=FALSE))

	Z0 <- if(search < 2 && side == 1) qnorm(p.bound, lower.tail=FALSE) else qnorm(p.bound/2, lower.tail=FALSE)

	if(is.null(rmat0)) rmat0 <- diag(1, k0) 
	ret <- cargs.proc(search, side, ncase0=ncase0, ncntl0=ncntl0, rmat0=rmat0, cor.numr=cor.numr, pool=pool)
	rmat0 <- ret$rmat0 ; rneff0 <- sqrt(ret$neff0) ; cor.numr <- ret$cor.numr ; search <- ret$search
	rmat <- matrix(rmat0[(z.sub != 0), (z.sub != 0)], k, k)
	ncase <- ncase0[z.sub != 0]
	ncntl <- if(search == 0) ncntl0 else ncntl0[z.sub != 0]
	rneff <- rneff0[z.sub != 0]

	x <- rep(0, k)
	kk <- length(sizes)
	xx <- rep(0, kk)
	cc <- c(0, cumsum(sizes[-kk]))
	
	nn <- (prod(sizes + 1) - 1)
	
	jj <- 1
	w.mat <- NULL
	a.vec <- NULL
	b.vec <- NULL
	inner.vec <- NULL
	prob.vec <- NULL
	while (jj <= nn)
	{
		pos <- max(which(xx < sizes))
		if(pos < kk)
		{
			xx[(pos + 1):kk] <- 0
			x[(cc[pos + 1] + 1):k] <- 0
		}

		xx[pos] <- xx[pos] + 1
		x[cc[pos] + xx[pos]] <- 1
		x0 <- rep(0, k0)
		x0[z.sub != 0] <- x
		set <- as.logical(x0)

#		Skip this subset
		if (!is.null(sub.def)) {
			if (!sub.def(set, sub.args)) { jj <- jj + 1 ; next }
           	}
		
		mm <- prod(choose(sizes, xx))
		rr <- if(!cor.numr ) (rneff0[set]) else mySolve(rmat0[set, set], rneff0[set])
		w <- rep(0, k0)
		w[set] <- rr/sqrt(matrix(rr, nrow=1) %*% as.matrix(rmat0[set, set]) %*% matrix(rr, ncol=1))
		w.mat <- rbind(w.mat, w)
		if(side == 1) { a <- t0 ; b <-  Inf ; inner <- TRUE }
		if(side == 2) { a <- -t0 ; b <-  t0 ; inner <- FALSE }

		prob <- min(pnorm(b, lower.tail=FALSE) + pnorm(a), 1)
		if(inner) prob <- max(0, 1 - prob)
		a.vec <- c(a.vec, a)
		b.vec <- c(b.vec, b)
		inner.vec <- c(inner.vec, inner)
		prob.vec <- c(prob.vec, prob * mm)

		jj <- jj + 1
	}
	con <- length(a.vec)
	if(search == 2 || p.bound < 1)
	{
		w.mat1 <- NULL
		a.vec1 <- NULL
		b.vec1 <- NULL
		inner.vec1 <- NULL
		prob.vec1 <- NULL
		for(j in 1:k0)
		{
			if(z.sub[j] == 0) next
			w <- rep(0, k0)
			w[j] <- 1
			w.mat1 <- rbind(w.mat1, w)
#			if(z.sub[j] != 0)
			{
				if(search < 2 && side == 1) { a <- (-Inf) ; b <- Z0 ; inner <- TRUE }
				if(search < 2 && side == 2) { a <- (-Z0) ; b <- Z0 ; inner <- TRUE }
				if(search == 2) { a <- (-Inf) ; b <- max(Z0, 0) ; inner <- TRUE }
			}
#			if(FALSE && z.sub[j] == 0)
#			{
#				if(search < 2 && side == 1) { a <- Z0 ; b <- Inf ; inner <- TRUE }
#				if(search < 2 && side == 2) { a <- (-Z0) ; b <- Z0 ; inner <- FALSE }
#				if(search == 2) { a <- 0 ; b <- Z0 ; inner <- TRUE }
#			}
			prob <- min(pnorm(b, lower.tail=FALSE) + pnorm(a), 1)
			if(inner) prob <- max(0, 1 - prob)
			a.vec1 <- c(a.vec1, a)
			b.vec1 <- c(b.vec1, b)
			inner.vec1 <- c(inner.vec1, inner)
			prob.vec1 <- c(prob.vec1, prob)
		}
		w.mat <- rbind(w.mat, w.mat1)
		a.vec <- c(a.vec, a.vec1)
		b.vec <- c(b.vec, b.vec1)
		inner.vec <- c(inner.vec, inner.vec1)
		prob.vec <- c(prob.vec, prob.vec1)
	}
	pv <- 0

	for(j in 1:con)
	{
		Z <- tubeSim(w.mat[j, ], a.vec[j], b.vec[j], inner=inner.vec[j], nsamp=nsamp, rmat=rmat0)
		pv <- pv + checkSim(Z, w.mat, a.vec, b.vec, inner.vec=inner.vec, rmat=rmat0) * prob.vec[j]
	}
	if(search == 2 || p.bound < 1)
	{
		pv1 <- 0
		con1 <- length(a.vec1)
		for(j in 1:con1)
		{
			Z1 <- tubeSim(w.mat1[j, ], a.vec1[j], b.vec1[j], inner=inner.vec1[j], nsamp=nsamp, rmat=rmat0)
			pv <- pv + checkSim(Z1, w.mat, a.vec, b.vec, inner.vec=inner.vec, rmat=rmat0) * prob.vec1[j]
			pv1 <- pv1 + checkSim(Z1, w.mat1, a.vec1, b.vec1, inner.vec=inner.vec1, rmat=rmat0) * prob.vec1[j]
		}

#		Prob0 <- (1 - pv1)
		Prob0 <- calcP0(Z0=Z0, z.sub=z.sub, search=search, side=side, p.bound=p.bound, rmat=rmat0,
                    cond.all=cond.all, NSAMP0=NSAMP0)

		if(pv1 > pv) {
			pv1 <- min((1 - Prob0), pv)
		}
		pv <- ((pv - pv1)/(1 - pv1))
	}
	pv
}

