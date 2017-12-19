load("data/hes1-largeExperimentSummary.rda", envir = .GlobalEnv)

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
  dimnames(outCoverage)[[2]] <- letters[1:dim(outCoverage)[[2]]]
  list(repetitionSize = outSize, coverage = outCoverage)
}

physicalSystem.choices <- c(
  "FitzHugh-Nagumo (FN) Model",
  "Oscillatory expression of the Hes1",
  "HIV model"
)

# for debug purpose ---------------------------------------------------------
input <- list(
  maternDf = 2.01,
  noise = 0.1,
  nobs = 51,
  pdfBaseDir = "/Users/shihaoyang/GoogleDrive/results"
)
