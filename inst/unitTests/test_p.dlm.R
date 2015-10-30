test_p.dlm <- function() {
 
  cor.def <- function(subset, neighbors, k, ncase, ncntl) {
    n <- ncol(neighbors)
    mat <- matrix(subset, nrow=k, ncol=n, byrow=FALSE)
    cor <- (mat + neighbors)*(1:k)/(k^2)
    cor <- colSums(cor)
    cor <- cor/max(cor)
    dim(cor) <- c(n, 1)

    cor  
  }
  sub.def <- function(logicalVec) {
    sum <- sum(logicalVec)  
    ret <- all(logicalVec[1:sum])
    ret
  }

  k     <- 3
  t.vec <- 1:k

  #res1 <- p.dlm(t.vec, k, 1, 2, cor.def=cor.def, sub.def=sub.def,
  #       cor.args=list(ncase=rep(1000, k), ncntl=rep(1000,k)))
  #res2 <- p.dlm(t.vec, k, 2, 2, cor.def=cor.def, sub.def=sub.def,
  #       cor.args=list(ncase=rep(1000, k), ncntl=rep(1000,k)))

  vec1 <- c(0.889214710, 0.133102631, 0.008039759)
  vec2 <- c(0.0761039412, 0.0119134489, 0.0006007941)

  #checkEquals(sum(abs(vec1 - vec1) < 1e-6), k)
  #checkEquals(sum(abs(vec2 - vec2) < 1e-6), k)
}
