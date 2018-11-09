library(gpds)

nobs.candidates <- c(41, 81, 121)
filllevel.candidates <- 0:2
temperPrior.candidates <- c(TRUE, FALSE)

indicatorArray <- array(FALSE, dim=c(length(temperPrior.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))

resultSummary <- list()
resultSdSummary <- list()
sizeSummary <- list()
files2zip <- c()
ciAll <- list()

for(arg in 1:length(indicatorArray)){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 2, any)],
    noise = c(0.15, 0.07),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 1e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = filllevel.candidates[apply(indicatorArray, 3, any)],
    modelName = "FN",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = FALSE,
    useGPmean = TRUE,
    forseTrueMean = TRUE,
    phase2 = FALSE,
    phase3 = FALSE,
    temperPrior = temperPrior.candidates[apply(indicatorArray, 1, any)],
    max.epoch = 10,
    epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  if(config$ndis > 201) next
  
  if(grepl("/n/",getwd())){
    baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
    config$seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9 # random seed on cluster
  }else{
    baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
  }
  outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                                "-ndis",ndis,"-",ifelse(config$temperPrior, "temperPrior", "unitHeatPrior"),"/"))
  
  allf <- list.files(outDir)
  allf <- allf[grep("rds", allf)]
  ci <- lapply(allf, function(f) try(readRDS(file.path(outDir, f))))
  names(ci) <- allf
  ci <- ci[sapply(ci, class) != "try-error"]
  ci <- sapply(ci, identity, simplify = "array")
  label <- dimnames(ci)[[3]]
  label <- gsub(".*-(phase[0-9]+)-.*", "\\1", label)
  
  resultSummary[[arg]] <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                                 identity, simplify = "array")
  resultSdSummary[[arg]] <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, sd)),
                                   identity, simplify = "array")
  sizeSummary[[arg]] <- table(label)
  ciAll[[arg]] <- ci
}
save.image("FN-largeExperimentSummary-forpaper.rda")


for(arg in 1:length(indicatorArray)){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 2, any)],
    noise = c(0.15, 0.07),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 1e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = filllevel.candidates[apply(indicatorArray, 3, any)],
    modelName = "FN",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = FALSE,
    useGPmean = TRUE,
    forseTrueMean = TRUE,
    phase2 = FALSE,
    phase3 = FALSE,
    temperPrior = temperPrior.candidates[apply(indicatorArray, 1, any)],
    max.epoch = 10,
    epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  if(config$ndis > 201) next
  
  flag <- with(config, paste0("nobs",nobs,"-ndis",ndis,"-",ifelse(config$temperPrior, "temperPrior", "unitHeatPrior")))
  
  names(resultSummary)[arg] <- flag
  names(resultSdSummary)[arg] <- flag
  names(sizeSummary)[arg] <- flag
  names(ciAll)[arg] <- flag
}

resultSummary[["nobs41-ndis41-temperPrior"]][,,"phase1"]
resultSummary[["nobs41-ndis41-temperPrior"]][,,"phase4"]

resultSummary[["nobs41-ndis41-unitHeatPrior"]][,,"phase1"]
resultSummary[["nobs41-ndis41-unitHeatPrior"]][,,"phase4"]

resultSummary[["nobs41-ndis81-temperPrior"]][,,c("phase1", "phase4")]
resultSummary[["nobs41-ndis81-unitHeatPrior"]][,,c("phase1", "phase4")]

resultSummary[["nobs41-ndis161-temperPrior"]][,,c("phase1", "phase4")]

printr <- function(x) format(round(x, 4), nsmall=4)
tablizeEstErr <- function(est, err){
  paste(format(round(est, 4), nsmall=4), "\\pm", format(round(err, 4), nsmall=4))
}

tab <- rbind(
  c("41", tablizeEstErr(
    resultSummary[["nobs41-ndis161-temperPrior"]][,,c("phase1")]["mean",],
    resultSdSummary[["nobs41-ndis161-temperPrior"]][,,c("phase1")]["mean",]
  )),
  c("81", tablizeEstErr(
    resultSummary[["nobs81-ndis161-temperPrior"]][,,c("phase1")]["mean",],
    resultSdSummary[["nobs81-ndis161-temperPrior"]][,,c("phase1")]["mean",]
  )),
  c("121", tablizeEstErr(
    resultSummary[["nobs121-ndis121-temperPrior"]][,,c("phase1")]["mean",],
    resultSdSummary[["nobs121-ndis121-temperPrior"]][,,c("phase1")]["mean",]
  ))
)
tab <- data.frame(tab)
colnames(tab) <- c("Samples", "a", "b", "c")
rownames(tab) <- NULL
print(xtable(tab), include.rownames=FALSE)


tab <- rbind(
  c("41", printr(resultSummary[["nobs41-ndis161-temperPrior"]][,,c("phase1")]["coverage",])),
  c("81", printr(resultSummary[["nobs81-ndis161-temperPrior"]][,,c("phase1")]["coverage",])),
  c("121", printr(resultSummary[["nobs121-ndis121-temperPrior"]][,,c("phase1")]["coverage",]))
)
tab <- data.frame(tab)
colnames(tab) <- c("Samples", "a", "b", "c")
rownames(tab) <- NULL
print(xtable(tab), include.rownames=FALSE)

