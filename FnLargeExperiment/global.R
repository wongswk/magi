load("data/FN-largeExperimentSummary.rda", envir = .GlobalEnv)


organizeOutput <- function(temperPrior, noise, nobs){
  indicatorArray <- array(FALSE, dim=c(length(temperPrior.candidates), 
                                       length(noise.candidates), 
                                       length(nobs.candidates),
                                       length(filllevel.candidates) ))
  indicatorArray[temperPrior.candidates == temperPrior, 
                 noise.candidates == noise, 
                 nobs.candidates == nobs,
                 ] <- TRUE
  indicatorArray[apply(indicatorArray, 1, any),
                 apply(indicatorArray, 2, any),
                 apply(indicatorArray, 3, any),
                 (nobs.candidates[apply(indicatorArray, 3, any)]-1)*2^filllevel.candidates+1 > 401
                 ] <- FALSE
  
  ndiscretize <- (nobs.candidates[apply(indicatorArray, 3, any)]-1)*2^filllevel.candidates+1
  ndiscretize <- ndiscretize[ndiscretize <= 401]
  
  indicatorVec <- which(indicatorArray)
  outSize <- do.call(rbind, sizeSummary[indicatorVec])
  outCoverage <- resultSummary[indicatorVec]
  rownames(outSize) <- names(outCoverage) <- paste0("discretize-", ndiscretize)
  outCoverage <- sapply(outCoverage, identity, simplify = "array")
  dimnames(outCoverage)[[2]] <- c("a", "b", "c")  
  list(repetitionSize = outSize, coverage = outCoverage)
}


physicalSystem.choices <- c(
  "FitzHugh-Nagumo (FN) Model"
  # "Oscillatory expression of the Hes1",
  # "HIV model"
)

# for debug purpose ---------------------------------------------------------
input <- list(
  temperPrior = TRUE,
  noise = 0.1,
  nobs = 51,
  pdfBaseDir = "/Users/shihaoyang/Workspace/DynamicSys/results/n/regal/kou_lab/shihaoyang/DynamicSys/results",
  phaseType = "phase1",
  variablePrint = c("avg 2.5% quantile" = "2.5%",
                    "avg median" = "50%",
                    "avg 97.5% quantile" = "97.5%",
                    "avg mean" = "mean",
                    "avg coverage" = "coverage")
)
