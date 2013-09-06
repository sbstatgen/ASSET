

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

p.tube <- function(t.vec, k, search, side, ncase, ncntl, pool, rmat, cor.numr, sizes=rep(1, k)
	, nsamp=50, sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL)
{
	if(!is.null(sub.def) && !is.function(sub.def)) stop("'sub.def' not recognized")
	
	if(k <= 0 || (k%%1 != 0)) stop("Number of studies k, not a positive integer!")
	if(!search %in% c(0, 1, 2)) stop("Invalid search option")
	if(search < 2 && !(side %in% c(1, 2))) stop("side should be 1 or 2")
	if(any(t.vec < 0)) warning("Negative input thresholds! Using absolute value.")
	
	nsnp <- length(t.vec)
	if(!is.null(dim(ncase)) && ncol(ncase) != 1 && ncol(ncase) != nsnp) stop("Unexpected dim of ncase")
	if(!is.null(dim(ncntl)) && ncol(ncntl) != 1 && ncol(ncntl) != nsnp) stop("Unexpected dim of ncntl")
	
	tube.pval(t.vec, k, search, side, ncase = ncase, ncntl = ncntl
			  , rmat = rmat, pool = pool, cor.numr = cor.numr, sizes = rep(1, k), nsamp = nsamp
			  , sub.def = sub.def, sub.args = sub.args)
}


p.dlm <- function(t.vec, k, search, side, cor.def = NULL, cor.args=NULL, sizes=rep(1, k), sub.def=NULL, sub.args=NULL
		, wt.def=NULL, wt.args=NULL)
{
	if(is.null(cor.def) || !is.function(cor.def))
	{
		if(search == 0) cor.def <- cor.types
		if(search == 1) cor.def <- cor.meta
		if(search == 2) cor.def <- coef.meta
	}
	if(!is.null(sub.def) && !is.function(sub.def)) stop("'sub.def' not recognized")
	
	if(k <= 0 || (k%%1 != 0)) stop("Number of studies k, not a positive integer!")
	if(!search %in% c(0, 1, 2)) stop("Invalid search option")
	if(search < 2 && !(side %in% c(1, 2))) stop("side should be 1 or 2")
	if(search == 2) side <- 1

	if(any(t.vec < 0)) warning("Negative input thresholds! Using absolute value.")

	dlm.pval(t.vec, k, search, side, cor.def, cor.args, sizes=sizes, sub.def=sub.def, sub.args=sub.args)
}

dlm.pval <- function(t.vec, k, search, side, cor.def, cor.args, sizes = rep(1, k)
				, sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL)
{
	if(k == 1) return(2 * pnorm(t.vec, lower.tail=FALSE))

	NT <- length(t.vec)
	
	NS <- (prod(sizes + 1) - 1)	
	if(any(t.vec < 0)) t.vec <- abs(t.vec)

	low <- t.vec
	high <- rep(Inf, NT)
		
	r.vec2 <- r.vec3 <- NULL
	
	x <- rep(0, k)	
	kk <- length(sizes)
	xx <- rep(0, kk)
	cc <- c(0, cumsum(sizes[-kk]))
	
	ss <- rep(0, NT)
	i <- 1
	while (i <= NS)
	{

		pos <- max(which(xx < sizes))
		if(pos < kk)
		{
			xx[(pos + 1):kk] <- 0
			x[(cc[pos + 1] + 1):k] <- 0			
		}
		xx[pos] <- xx[pos] + 1
		x[cc[pos] + xx[pos]] <- 1
		set <- as.logical(x)
		NXX <- prod(choose(sizes, xx))
		
#		DEBUG
#		print(c(i, pos, cc[pos], xx, NXX))

#		Skip this subset
		if(!is.null(sub.def) && !do.call(sub.def, c(list(set), sub.args))) { i <- i + 1 ; next }

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
#		Valid neighbors
		if(!is.null(sub.def)) 
		{
		nn.set <- as.logical(nn.x)
		dim(nn.set) <- dim(nn.x)
		nn.sub <- nn.sub & apply(nn.set, 2, function(set1) do.call(sub.def, c(list(set1), sub.args)))
		}
#		Add a term to points with integral still below 1
		t.sub <- (1:NT)[ss < 1]

		if(length(t.sub) > 0 && search < 2)
		{
			cor.args$ncase <- matrix(cor.args$ncase, nrow=k, ncol=NT, byrow=FALSE)
			int <- sapply(t.sub, function(j)
						{
							ncor.args <- cor.args
							ncor.args$ncase <- cor.args$ncase[, j]
							ncor.args$ncntl <- ncor.args$ncntl[j]
							qxd.prod.int(low[j], high[j], search, qxd.prod.cor, side=side, x, matrix(nn.x[, nn.sub]
							, nrow=k, byrow=FALSE), cor.def, ncor.args)
						})
			ss[t.sub] <- ss[t.sub] + NXX * int
		}
		if(length(t.sub) > 0 && search == 2)
		{
			
#			Older code 2-sided search with 3 correlations instead of coefficients
#			nn.x2 <- diag(1, k)
#			sgn <- 1 - 2 * x
#			if(sum(x) == 1)
#			{
#				nn.x2 <- nn.x2[, (x == 0)]
#				sgn <- sgn[x == 0]
#				if(is.null(dim(nn.x2))) nn.x2 <- matrix(nn.x2, nrow = k)
#			}	
#
#			r.vec2 <- sapply(1:length(nn.sub), function(j) do.call(cor.def, c(list(x, nn.x2[, j], k, side), cor.args)))
#			r.vec3 <- r.vec * r.vec2 + sgn[nn.sub] * sqrt((1 - r.vec * r.vec) * (1 - r.vec2 * r.vec2))
#			rho.mat <- cbind(r.vec, r.vec2, r.vec3)

			cor.args$ncase <- matrix(cor.args$ncase, nrow=k, ncol=NT, byrow=FALSE)
			int <- sapply(t.sub, function(j)
						{
							ncor.args <- cor.args
							ncor.args$ncase <- cor.args$ncase[, j]
							ncor.args$ncntl <- ncor.args$ncntl[j]
							qxd.prod.int(low[j], high[j], search, qxd.prod.coef, side=side, x=x, nn.x=matrix(nn.x[, nn.sub]
								, nrow=k, byrow=FALSE), cor.def, cor.args)
						})
						
			ss[t.sub] <- ss[t.sub] + NXX * int
		}
		if(length(t.sub) == 0) break
		i <- i + 1
	}

	ret <- pmin(1, ss)

	ret
}
	
qxd.prod.int <- function(low1, high1, search, int.func, side, x, nn.x, cor.def, cor.args)
{
	val <- 0

	if(low1 < Inf)
	{
		int <- integrate(function(zv)
						{
					 
						 prob1 <- side * int.func(zv=zv, search, side, x, nn.x, cor.def, cor.args)

						 (prob1 * dnorm(zv))
						}, lower=low1, upper=high1)
		if(int$message != "OK") { warning(int$message) ; val <- NA }
		else val <- int$value
	}
	val
}

qxd.prod.coef <- function(zv, search, side, x, nn.x, cor.def, cor.args)
{
	k <- length(x)
	npts <- length(zv)
	if(is.null(nn.x) || ncol(nn.x) == 0) return(rep(1, npts))	
	beta.mat <- matrix(drop(do.call(cor.def, c(list(x, nn.x, k), cor.args))), ncol=2, byrow=FALSE)


	if(any(is.na(beta.mat))) stop("NA in beta.mat")
	if(is.null(dim(beta.mat)) || ncol(beta.mat) != 2) stop("Expected two columns in beta.mat")
	
	kk <- nrow(beta.mat)
	b.s <- beta.mat[, 1]
	b.k <- beta.mat[, 2]

	lrr <- 0
	zv1 <- zv
	
	T1 <- (1/b.k) %o% zv1
	T2 <- (b.s/b.k) %o% zv1
	
	b.k2 <- b.k %o% rep(1, npts)
	
	r.vec <- (1 - b.s^2 - b.k^2)/(2 * b.s * b.k)
#	rho <- ((1 - b.k^2) + b.s^2)/(2 * b.s)
	
	mu.z <- r.vec %o% zv
	sd <- sqrt(1 - r.vec^2)
	
	if(search < 2)
	{
		if(side == 2)
		{
			low <- ifelse(b.k2 > 0, -T1 - T2, T1 - T2)
			high <- ifelse(b.k2 > 0, T1 - T2, - T1 - T2)
		}
		else
		{
			low <- ifelse(b.k2 > 0, -Inf, T1 - T2)
			high <- ifelse(b.k2 > 0, T1 - T2, Inf)
		}
		p1 <- ifelse(T1 == Inf, 0, pnorm(high, mean = mu.z, sd = sd))
		p0 <- ifelse(T1 == Inf, 0, pnorm(low, mean = mu.z, sd = sd))

	}
	if(search == 2)
	{
		low <- ifelse(b.k2 > 0, 0, pmax(0, T1 - T2))
		high <- ifelse(b.k2 > 0, T1 - T2, Inf)
		p1 <- ifelse(T1 == Inf | (b.k2 < 0 & T1 < T2), 0, pnorm(high, mean = mu.z, sd = sd))
		p0 <- ifelse(T1 == Inf | (b.k2 < 0 & T1 < T2), 0, pnorm(low, mean = mu.z, sd = sd))
	}

	pp <- ifelse(T1 == Inf, 1, (p1 - p0))		
	ss <- apply(pp, 2, function(x) { if(search == 2) (prod(2 * x)) else prod(x) })
	
	
#	DEBUG
#	print(rbind(zv,ss))
	
	if(any(is.na(ss))) stop("NA in integrand")
	ss
}

qxd.prod.cor <- function(zv, search, side, x, nn.x, cor.def, cor.args)
{
	k <- length(x)
	npts <- length(zv)
	if(is.null(nn.x) || ncol(nn.x) == 0) return(rep(1, npts))

	rho.mat <- do.call(cor.def, c(list(x, nn.x, k), cor.args))
	kk <- nrow(rho.mat)

	if(is.null(dim(rho.mat))) stop("rho.mat should be a matrix")
	if(search == 2 && ncol(rho.mat) != 3) stop("Expected 3 columns in rho.mat")
	if(search < 2 && ncol(rho.mat) != 1) stop("Expected 1 column in rho.mat")
	
	r.vec <- rho.mat[, 1]

	lrr <- 0
	zv1 <- zv
	
	if(search == 2)
	{
		zv11 <- rep(1, kk) %o% zv1
		r.vec2 <- rho.mat[, 2]
		r.vec3 <- rho.mat[, 3]
		beta <- (r.vec3 - r.vec * r.vec2)/(1 - r.vec * r.vec)
		beta2 <- beta %o% rep(1, kk)
	}
		
	mu.z <- r.vec %o% zv
	if(search == 2) mu.z1 <- (r.vec + beta * r.vec2) %o% zv

	sd <- sqrt(1 - r.vec^2)

	if(search < 2)
	{
		low <- rep(1, kk) %o% -abs(zv1)
		high <- rep(1, kk) %o% abs(zv1)
		p1 <- pnorm(high, mean = mu.z, sd = sd)
		if(side == 2) p0 <- pnorm(low, mean = mu.z, sd = sd)
		else p0 <- rep(0, kk, npts)
	}
	if(search == 2)
	{
		low <- ifelse(beta2 < 0, mu.z1, -Inf)
		high <- ifelse(beta2 < 0, abs(zv11), pmin(abs(zv11), mu.z1))
		p1 <- pnorm(high, mean = mu.z, sd = sd)
		p0 <- pnorm(low, mean = mu.z, sd = sd)
	}
	
	pp <- matrix((p1 - p0), nrow=kk, ncol=npts, byrow=FALSE)
	
	ss <- apply(pp, 2, function(x) { if(search == 2) (prod(2 * x)) else prod(x) })
	ss <- ifelse(zv1 == Inf, 1, ss)

#	DEBUG
#	print(rbind(zv,ss))

	if(any(is.na(ss))) stop("NA in integrand")
	ss
}


coef.meta <- function(x1, X2, k, ncase, ncntl, rmat=NULL, cor.numr=FALSE)
{
	if(is.null(rmat)) rmat <- diag(1, k)
	rneff <- sqrt(ncase * ncntl/(ncase + ncntl))

	
	nsnp <- length(ncase) %/% k
	if(is.null(dim(rneff))) dim(rneff) <- c(k, nsnp)
	
	x1 <- as.logical(x1)
	X2 <- as.logical(X2)

	nn <- length(X2)/k
	if(is.null(dim(X2))) dim(X2) <- c(k, nn)

	beta <- array(NA, c(nn, 2, nsnp))
	
	if(is.null(rmat)) rmat <- diag(1, k)
	
	xcom <- which(x1)
	xdel <- apply(X2, 2, function(x2) (1:k)[xor(x1, x2)])
	sgn <- (1 - 2 * x1[xdel]) %o% rep(1, nsnp)
	ncom <- length(xcom)
	ndel <- length(xdel)

	if(cor.numr)
	{
		rneff.com <- mySolve(matrix(rmat[xcom, xcom], ncom, ncom), matrix(rneff[xcom, ], ncom, nsnp, byrow=FALSE))
		rneff.del <- mySolve(matrix(rmat[xdel, xdel], ndel, ndel), matrix(rneff[xdel, ], ndel, nsnp, byrow=FALSE))
	} else 
	{
		rneff.com <- matrix(rneff[xcom, ], ncom, nsnp, byrow=FALSE)
		rneff.del <- matrix(rneff[xdel, ], ndel, nsnp, byrow=FALSE)
	}
	
	if(!is.vector(xdel)) stop("Coefficients defined for immediate neighbors only")
		
	sxx <- diag(t(rneff.com) %*% matrix(rmat[xcom, xcom], ncom, ncom) %*% rneff.com)
	sxd <- (matrix(rmat[xdel, xcom], ndel, ncom) %*% rneff.com) * (sgn * rneff.del)
	sdd <- rneff.del * rneff.del

	sall <- (rep(1, nn) %o% sxx) + 2 * sxd + sdd

	beta[, 1, ] <- ifelse(sall == 0, NA, sqrt(sxx/sall))
	beta[, 2, ] <- ifelse(sall == 0, NA, sgn * sqrt(sdd/sall))
	if(is.null(dim(beta))) dim(beta) <- c(nn, 2, nsnp)

#	rho <- rep(NA, ndel)
#	if(sxx != 0) rho <- ifelse(sall == 0, NA, (sxx + sxd)/sqrt(sxx * sall))
#	rho

	beta
}

cor.meta <- function(x1, X2, k, ncase, ncntl, rmat=NULL, cor.numr=FALSE)
{
	if(is.null(rmat)) rmat <- diag(1, k)
	rneff <- sqrt(ncase * ncntl/(ncase + ncntl))
	nsnp <- length(ncase) %/% k
	if(is.null(dim(rneff))) dim(rneff) <- c(k, nsnp)

	x1 <- as.logical(x1)
	n1 <- sum(x1)
	if(cor.numr) rneff1 <- mySolve(matrix(rmat[x1, x1], n1, n1), matrix(rneff[x1, ], n1, nsnp, byrow=FALSE))
	else rneff1 <- matrix(rneff[x1, ], n1, nsnp, byrow=FALSE)

	cor.func <- function(x2)
	{	
		x2 <- as.logical(x2)
		n2 <- sum(x2)
		if(cor.numr) rneff2 <- mySolve(matrix(rmat[x2, x2], n2, n2), matrix(rneff[x2, ], n2, nsnp, byrow=FALSE))
		else rneff2 <- rneff[x2, ]
		sx12 <- t(matrix(rneff1, ncol=nsnp)) %*% matrix(rmat[x1, x2], n1, n2) %*% matrix(rneff2, ncol = nsnp)
		sx11 <- t(matrix(rneff1, ncol=nsnp)) %*% matrix(rmat[x1, x1], n1, n1) %*% matrix(rneff1, ncol = nsnp)
		sx22 <- t(matrix(rneff2, ncol=nsnp)) %*% matrix(rmat[x2, x2], n2, n2) %*% matrix(rneff2, ncol = nsnp)

		rho <- ifelse(sx11 == 0 | sx22 == 0, NA, sx12/sqrt(sx11 * sx22))
		rho
	}

	rho.mat <- matrix(apply(X2, 2, cor.func), ncol = nsnp, byrow = TRUE)

	rho.mat
}

cor.types <- function(x1, X2, k, ncase, ncntl, pool)
{
	nsnp <- length(ncntl)
	if(is.null(dim(ncase))) dim(ncase) <- c(k, nsnp)

	x1 <- as.logical(x1)
	x2 <- as.logical(X2)
	
	if(is.null(dim(X2))) dim(X2) <- c(k, length(X2)/k)
	nn <- ncol(X2)
	
#	if(any(x1 & !x2) && any(!x1 & x2)) stop("Correlation defined for subset/superset only")
	
	ncase.com <- apply(X2, 2, function(x2) colSums(matrix(ncase[x1 & x2,], ncol=nsnp)))
	ncase.all <- apply(X2, 2, function(x2) colSums(matrix(ncase[x1 | x2,], ncol=nsnp)))
	N <- (colSums(ncase) + ncntl) %o% rep(1, nn)
	
	if(pool == FALSE) rho <- sqrt((ncase.com * (ncase.all + ncntl))/(ncase.all * (ncase.com + ncntl)))
	else rho <- sqrt((ncase.com * (N - ncase.all))/((ncase.all) * (N - ncase.com)))

	rho <- matrix(rho, nn, nsnp, byrow=TRUE)
	rho
}

tubeSim <- function(x, t0, nsamp, rmat0, neff0, cor.numr=FALSE)
{
	k <- length(x)
	nsnp <- length(t0)

	Z <- array(0, c(k, nsamp, nsnp))
	x <- as.logical(x)
	r1 <- (which(x))[1]
	r2 <- (1:k)[-r1]

	for(j in 1:nsnp)
	{
		if(j == 1 || ncol(rmat0) > 1 || ncol(neff0) > 1)
		{
			rmat <- matrix(rmat0[, (j - 1) %% ncol(rmat0) + 1], k, k)
			rneff <- sqrt(neff0[, (j - 1) %% ncol(neff0) + 1])
			if(cor.numr)
			{
				rr <- rep(0, k)
				rr[x] <- mySolve(rmat[x, x], rneff[x])
			}
			else rr <- (rneff * x)
		
			A <- diag(k)
			A[r1, ] <- rr/sqrt(sum(matrix(rr, nrow=1) %*% rmat %*% matrix(rr, ncol=1)))
			A <- A[c(r1, r2), ]

			Sigma <- A %*% rmat %*% t(A)
			Sig <- Sigma[(2:k), (2:k)] -  ((1/Sigma[1, 1]) * outer(Sigma[(2:k), 1], Sigma[(2:k), 1]))	

			ev <- eigen(Sig, symmetric = TRUE)

			if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1])))
			{
				warning("Sigma is numerically not positive definite")
			}
			Sigroot <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
		}
	
		W1 <- rtnorm(nsamp, mean = 0, sd = sqrt(Sigma[1, 1]), lower = t0[j])
		mu <- outer(Sigma[(2:k), 1]/Sigma[1, 1], W1)
		
		W2 <- t(matrix(rnorm(nsamp * (k - 1)), nrow = nsamp) %*% Sigroot) + mu				
		
		W <- rbind(W1, W2)
	
		Z[, , j] <- solve(A) %*% W
	}
	Z
}


checkSim1 <- function(Z, t1, t2, rmat0, neff0, cor.numr=FALSE, sub.def=NULL, sub.args=NULL)
{
	k <- dim(Z)[1]
	nsamp <- dim(Z)[2]
	nsnp <- dim(Z)[3]
					

	x <- rep(0, k)
	fm <- matrix(0, nsamp, nsnp)
	i <- 1
	while(i <= ((2^k)-1))
	{		
		pos <- max(which(x == 0))
		x[pos:k] <- 0
		x[pos] <- 1
		set <- as.logical(x)

#		Skip this subset
		if(!is.null(sub.def) && !do.call(sub.def, c(list(set), sub.args))) { i <- i + 1 ; next }
		
		for(j in 1:nsnp)
		{
			if(cor.numr)
			{
				if(j == 1 || ncol(rmat0) > 1 || ncol(neff0) > 1)
				{
					rmat <- matrix(rmat0[, (j - 1) %% ncol(rmat0) + 1], k, k)
					rneff <- sqrt(neff0[, (j - 1) %% ncol(neff0) + 1])
					denr <- sqrt(matrix(rneff[set], nrow=1) %*% mySolve(rmat[set, set], matrix(rneff[set], ncol=1)))
				}
				numr <- drop(matrix(rneff[set], nrow = 1) %*% mySolve(rmat[set, set], Z[set, , j]))
			}
			else
			{
				if(j == 1 || ncol(rmat0) > 1 || ncol(neff0) > 1)
				{
					rmat <- matrix(rmat0[, (j - 1) %% ncol(rmat0) + 1], k, k)
					rneff <- sqrt(neff0[, (j - 1) %% ncol(neff0) + 1])

					denr <- sqrt(matrix(rneff[set], nrow=1) %*% rmat[set, set] %*% matrix(rneff[set], ncol=1))

				}

				numr <- drop(matrix(rneff[set], nrow = 1) %*% Z[set, , j])
			}
	
			f1 <- ((numr/denr) > t1[j])
			f2 <- ((-numr/denr) > t2[j])

			fm[, j] <- fm[, j] + (f1 + f2)
		}

		i <- i + 1		
	}

	gm <- apply(fm, 2, function(x) mean(ifelse(x > 0, 1/x, 0)))

	gm
	
}

checkSim2 <- function(Z1, Z2, t1, t2, rmat0, neff0, cor.numr=FALSE, sub.def=NULL, sub.args=NULL)
{
	k <- dim(Z1)[1]
	nsamp <- dim(Z1)[2]
	nsnp <- dim(Z1)[3]
						
	x <- rep(0, k)
	
	fm1 <- matrix(0, nsamp, nsnp)
	fm2 <- matrix(0, nsamp, nsnp)
	fm3 <- matrix(0, nsamp, nsnp)
	
	i <- 1
	while(i <= ((2^k)-1))
	{		
		pos <- max(which(x == 0))
		x[pos:k] <- 0
		x[pos] <- 1
		set <- as.logical(x)
					
#		Skip this subset
		if(!is.null(sub.def) && !do.call(sub.def, c(list(set), sub.args))) { i <- i + 1 ; next }

		for(j in 1:nsnp)
		{
			if(cor.numr)
			{
				if(j == 1 || ncol(rmat0) > 1 || ncol(neff0) > 1)
				{
					rmat <- matrix(rmat0[, (j - 1) %% ncol(rmat0) + 1], k, k)
					rneff <- sqrt(neff0[, (j - 1) %% ncol(neff0) + 1])
					denr <- sqrt(matrix(rneff[set], nrow=1) %*% mySolve(rmat[set, set], matrix(rneff[set], ncol=1)))
				}
				numr1 <- drop(matrix(rneff[set], nrow = 1) %*% mySolve(rmat[set, set], Z1[set, , j]))
				numr2 <- drop(matrix(rneff[set], nrow = 1) %*% mySolve(rmat[set, set], Z2[set, , j]))
			}
			else
			{
				if(j == 1 || ncol(rmat0) > 1 || ncol(neff0) > 1)
				{
					rmat <- matrix(rmat0[, (j - 1) %% ncol(rmat0) + 1], k, k)
					rneff <- sqrt(neff0[, (j - 1) %% ncol(neff0) + 1])
					denr <- sqrt(matrix(rneff[set], nrow=1) %*% rmat[set, set] %*% matrix(rneff[set], ncol=1))
				}
				numr1 <- drop(matrix(rneff[set], nrow = 1) %*% Z1[set, , j])
				numr2 <- drop(matrix(rneff[set], nrow = 1) %*% Z2[set, , j])
			}
			
			f11 <- ((numr1/denr) > t1[j, 1])
			f12 <- ((-numr1/denr) > t2[j, 1])
			fm1[, j] <- fm1[, j] + (f11 + f12)
			
			f21 <- ((numr2/denr) > t1[j, 2])
			f22 <- ((-numr2/denr) > t2[j, 2])
			fm2[, j] <- fm2[, j] + (f21 + f22)
			
			f31 <- ((numr2/denr) > t1[j, 3])
			f32 <- ((-numr2/denr) > t2[j, 3])
			fm3[, j] <- fm3[, j] + (f31 + f32)			
		}
		
		i <- i + 1		
	}
	gm1 <- apply(fm1, 2, function(x) mean(ifelse(x > 0, 1/x, 0)))
	gm2 <- apply(fm2, 2, function(x) mean(ifelse(x > 0, 1/x, 0)))
	gm3 <- apply(fm3, 2, function(x) mean(ifelse(x > 0, 1/x, 0)))
		
	cbind(gm1, gm2, gm3)
}


tube.pval <- function(t.vec, k, search, side, ncase, ncntl, pool, rmat = NULL
					  , cor.numr = FALSE, sizes = rep(1, k), nsamp = 50, sub.def=NULL, sub.args=NULL)
{	
	if(k == 1) return(side * pnorm(t.vec, lower.tail=FALSE))
	
	nsnp <- length(t.vec)

	if(is.null(dim(ncase))) ncase <- matrix(ncase, ncol = nsnp, byrow=FALSE)
	if(is.null(dim(ncntl))) ncntl <- matrix(ncntl, ncol = nsnp, byrow=FALSE)

	if(search == 1 || search == 2)
	{
		neff <- (ncase * ncntl)/(ncase + ncntl)
		if(is.null(rmat)) rmat <- diag(k)
		rmat <- matrix(as.vector(rmat), ncol=1)
	}

	if(search == 0)
	{
		N <- apply(ncase, 2, sum) + ncntl
		ncntl1 <- matrix(ncntl, k, nsnp, byrow=TRUE)
		if(pool == FALSE)
		{
			neff <- (ncase * ncntl1)/(ncase + ncntl1)
			rmat <- sapply(1:nsnp
						, function(j)
						   {
						   rtmp <- outer(sqrt(ncase[,j]/(ncase[,j] + ncntl[1, j])), sqrt(ncase[,j]/(ncase[,j]+ncntl[1, j])))
						   diag(rtmp) <- 1
						   as.vector(rtmp)
						   })
		}
		if(pool == TRUE) 
		{
			N1 <- ncase + ncntl1
			neff <- (ncase * (N1 - ncase))/N1
			rmat <- sapply(1:nsnp
						, function(j)
						   {
								rtmp <- (-1) * outer(ncase[, j]/sqrt(neff[, j]), ncase[, j]/sqrt(neff[, j]))/N[j]
								diag(rtmp) <- 1
								as.vector(rtmp)
							})
		}
#		Note: hard coded
		cor.numr <- TRUE
		search <- 1
	}
	
	kk <- length(sizes)
	cc <- c(0, cumsum(sizes[-kk]))
	
	nn <- (prod(sizes + 1) - 1)
	
	x <- rep(0, k)
	xx <- rep(0, kk)
	
	prob1 <- pnorm(t.vec, lower.tail=FALSE)
	zeros <- rep(0, nsnp)
	infs <- rep(Inf, nsnp)				
					
	jj <- 1
	pv <- rep(0, nsnp)
	pv0 <- rep(0, nsnp)

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
		set <- as.logical(x)

#		Skip this subset
		if(!is.null(sub.def) && !do.call(sub.def, c(list(set), sub.args))) { jj <- jj + 1 ; next }
				
		mm <- prod(choose(sizes, xx))
	
		if(search == 1)
		{
			Z1 <- tubeSim(x = x, t0 = t.vec, nsamp = nsamp, rmat0 = rmat, neff0 = neff)

			p1 <- checkSim1(Z1, t1 = t.vec, t2 = t.vec, rmat0 = rmat, neff0 = neff
							, sub.def = sub.def, sub.args = sub.args) * (prob1 * mm)

			p2 <- if(side == 2) p1 else rep(0, nsnp)
			pv <- pv + (p1 + p2)
		}
		if(search == 2)
		{
			Z1 <- tubeSim(x = (x), t0 = t.vec, nsamp = nsamp, rmat0 = rmat, neff0 = neff)
			
			Z2 <- tubeSim(x = (x), t0 = zeros, nsamp = nsamp, rmat0 = rmat, neff0 = neff)

			res <- checkSim2(Z1, Z2, t1 = cbind(t.vec, zeros, zeros), t2 = cbind(zeros, t.vec, infs)
							 , rmat0 = rmat, neff0 = neff, sub.def = sub.def, sub.args = sub.args)
			p1 <- res[, 1] * (prob1 * mm)
			p2 <- res[, 2]  * (0.5 * mm)
			p0 <- res[, 3] * (0.5 * mm)
			
			pv <- pv + (p1 + (p2 - p0))
			pv0 <- pv0 + p0
		}
		
		jj <- jj + 1
	}
#	pv0 <- 1 - (1/2^k)	
#   print(c(pv, pv0, 1 - (1/(2^k))))
	
	if(search == 2) pv <- pmax(pv, 0)/pmax((1/2^k), 1 - pv0)
	
	pv
}

