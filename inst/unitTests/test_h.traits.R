test_h.traits <- function() {

  data(ex_trait, package="ASSET")
  data       <- data[1, , drop=FALSE]
  snps       <- as.vector(data[, "SNP"])
  traits.lab <- paste("Trait_", 1:6, sep="")
  beta.hat   <- as.matrix(data[, paste(traits.lab, ".Beta", sep="")])
  sigma.hat  <- as.matrix(data[, paste(traits.lab, ".SE", sep="")])
  cor        <- list(N11=N11, N00=N00, N10=N10)
  ncase      <- diag(N11)
  ncntl      <- diag(N00)
  res        <- h.traits(snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=cor, 
                 cor.numr=FALSE, search=NULL, side=2, meta=TRUE, 
                 zmax.args=NULL, meth.pval="DLM")
 
  lvec1 <- c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  lvec2 <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
  lvec3 <- c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)


  checkEqualsNumeric(res$Meta$pval, 0.26756223, tolerance=1.0e-6)
  checkEqualsNumeric(res$Subset.1sided$pval, 1.877814e-02, tolerance=1.0e-6)
  checkEquals(sum(res$Subset.1sided$pheno == lvec1), 6)
  checkEqualsNumeric(res$Subset.2sided$pval, 1.248431e-03, tolerance=1.0e-6)
  checkEqualsNumeric(res$Subset.2sided$pval.1, 0.0539668, tolerance=1.0e-6)
  checkEqualsNumeric(res$Subset.2sided$pval.2, 0.002316304, tolerance=1.0e-6)
  checkEquals(sum(res$Subset.2sided$pheno.1 == lvec2), 6)
  checkEquals(sum(res$Subset.2sided$pheno.2 == lvec3), 6)
}
