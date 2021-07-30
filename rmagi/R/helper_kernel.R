# testhat helper functions

#' calculate number of eigen values to preserve based on frobenius norm
#' @noRd
truncEigen <- function(eigenValues, frobeniusNormApprox = 0.99){
  frobeniusNorm <- sum(eigenValues^2)
  frobeniusNormTrunc <- cumsum(eigenValues^2)
  min(which(frobeniusNormTrunc/frobeniusNorm > frobeniusNormApprox))
}

#' truncate gpCov object for low rank approximation
#'
#' the largest 1-over-eigenvalue will be preserved, and the rest will be deleted
#' for a low-rank representation from spectral decomposition
#'
#' mphi SVD not helping because complexity is 2mn, comparing to original n^2
#' however we need m to be around n/2 to preserve accuracy due to low decrease d
#'
#' this idea is not used in the end because band matrix approximation is better
#'
#' @param cKeep number of eigen values to keep for C matrix
#' @param kKeep number of eigen values to keep for K matrix
#' @noRd
truncCovByEigen <- function(gpCov, cKeep, kKeep){
  cKeepId <- (ncol(gpCov$CeigenVec)-cKeep+1):ncol(gpCov$CeigenVec)
  gpCov$Ceigen1over <- gpCov$Ceigen1over[cKeepId]
  gpCov$CeigenVec <- gpCov$CeigenVec[,cKeepId]

  kKeepId <- (ncol(gpCov$KeigenVec)-kKeep+1):ncol(gpCov$KeigenVec)
  gpCov$Keigen1over <- gpCov$Keigen1over[kKeepId]
  gpCov$KeigenVec <- gpCov$KeigenVec[,kKeepId]

  # mKeepId <- 1:mKeep
  # gpCov$mphiu <- gpCov$mphiu[mKeepId,]
  # gpCov$mphid <- gpCov$mphid[mKeepId]
  # gpCov$mphiv <- gpCov$mphiv[mKeepId,]

  gpCov
}
