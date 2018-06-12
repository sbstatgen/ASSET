test_h.types <- function() {
 
  data(ex_types, package="ASSET")
  snps     <- "SNP_1"
  adj.vars <- c("CENTER_1", "CENTER_2", "CENTER_3")
  types    <- paste("SUBTYPE_", 1:5, sep="")
  res      <- h.types(data, "TYPE", snps, adj.vars, types, "SUBTYPE_0", subset=NULL, 
        method="case-control", side=2, logit=TRUE, test.type="Score", 
        zmax.args=NULL, meth.pval="DLM", pval.args=NULL)

  lvec1 <- c(TRUE, FALSE, FALSE, FALSE, FALSE)

  checkEqualsNumeric(res$Overall.Logistic$pval, 0.003239651, tolerance=1.0e-6)
  checkEqualsNumeric(res$Subset.Case.Control$pval, 1.223364e-05, tolerance=1.0e-6)
  checkEquals(sum(res$Subset.Case.Control$pheno == lvec1), 5)

}
