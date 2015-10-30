test_z.max <- function() {

  meta.def <- function(logicalVec, SNP.list, arg.beta, arg.sigma) {
    beta <- as.matrix(arg.beta[SNP.list, logicalVec])
    se   <- as.matrix(arg.sigma[SNP.list, logicalVec])
    test <- (beta/se)^2
    ret  <- apply(test, 1, max)
    list(z=ret) 
  }
  sub.def <- function(logicalVec, sub.args) {
    sum <- sum(logicalVec)  
    ret <- all(logicalVec[1:sum])
    ret
  }

  k        <- 5
  snp.vars <- 1:3
  nsnp     <- length(snp.vars)
  beta     <- matrix(1:(k*nsnp), nrow=nsnp)
  sigma    <- matrix(rep(1, k*nsnp), nrow=nsnp) 

  meta.args <- list(arg.beta=beta, arg.sigma=sigma)
  
  res1 <- z.max(k, snp.vars, 2, meta.def, meta.args, sub.def=sub.def)

  vec1 <- c(169, 196, 225)

  #checkEquals(sum(abs(res1$opt.z - vec1) < 1e-6), nsnp)
 
}
