load("data/largeExperimentSummary.rda")

organizeOutput <- function(maternDf, noise, nobs){
  indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                       length(noise.candidates), 
                                       length(nobs.candidates),
                                       length(filllevel.candidates) ))
  indicatorArray[kernel.candidates == ifelse(maternDf == 2.5, "matern", "generalMatern"), 
                 noise.candidates == noise, 
                 nobs.candidates == nobs,
                 ] <- TRUE
  indicatorArray[apply(indicatorArray, 1, any),
                 apply(indicatorArray, 2, any),
                 apply(indicatorArray, 3, any),
                 (nobs.candidates[apply(indicatorArray, 3, any)]-1)*2^filllevel.candidates+1 > 801
                 ] <- FALSE
  
  ndiscretize <- (nobs.candidates[apply(indicatorArray, 3, any)]-1)*2^filllevel.candidates+1
  ndiscretize <- ndiscretize[ndiscretize <= 801]
  
  indicatorVec <- which(indicatorArray)
  outSize <- do.call(rbind, sizeSummary[indicatorVec])
  outCoverage <- resultSummary[indicatorVec]
  rownames(outSize) <- names(outCoverage) <- paste0("discretize-", ndiscretize)
  outCoverage <- sapply(outCoverage, identity, simplify = "array")
  dimnames(outCoverage)[[2]] <- c("a", "b", "c")  
  list(repetitionSize = outSize, coverage = outCoverage)
}
