library(gpds)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:6
temperPrior.candidates <- c(TRUE, FALSE)

indicatorArray <- array(FALSE, dim=c(length(temperPrior.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))

resultSummary <- list()
sizeSummary <- list()
files2zip <- c()

for(arg in 1:length(indicatorArray)){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 3, any)],
    noise = rep(noise.candidates[apply(indicatorArray, 2, any)], 2),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 2e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
    modelName = "FN",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = FALSE,
    useGPmean = TRUE,
    forseTrueMean = FALSE,
    phase2 = TRUE,
    temperPrior = temperPrior.candidates[apply(indicatorArray, 1, any)]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  if(config$ndis > 401) next
  
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster
  outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                                "-ndis",ndis,"-",ifelse(config$temperPrior, "temperPrior", "unitHeatPrior"),"/"))
  
  allf <- list.files(outDir)
  allf <- allf[grep("rds", allf)]
  ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), 
               simplify = "array")
  label <- dimnames(ci)[[3]]
  label <- gsub(".*-(phase[0-9])-.*", "\\1", label)
  
  resultSummary[[arg]] <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                          identity, simplify = "array")
  sizeSummary[[arg]] <- table(label)
}
save.image("FN-largeExperimentSummary.rda")


length(sizeSummary)
length(which(sapply(sizeSummary, is.null)))
length(indicatorArray)
summary(do.call(rbind, sizeSummary))
dim(do.call(rbind, sizeSummary))

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

organizeOutput(temperPrior = TRUE, noise = 0.1, nobs = 51)

organizeOutput(temperPrior = FALSE, noise = 0.1, nobs = 51)

organizeOutput(temperPrior = TRUE, noise = 0.1, nobs = 201)

organizeOutput(temperPrior = FALSE, noise = 0.1, nobs = 401) 

organizeOutput(temperPrior = TRUE, noise = 0.5, nobs = 201)

organizeOutput(temperPrior = FALSE, noise = 1.0, nobs = 201)

# summarize refilled priorTempered first phase --------------------------------
library(gpds)
library(parallel)
nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:4
kernel.candidates <- c("generalMatern", "matern") # basically df in matern

indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))

resultSummary <- list()
sizeSummary <- list()
files2zip <- c()

allSummaryList <- mclapply(1:length(indicatorArray), function(arg){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 3, any)],
    noise = noise.candidates[apply(indicatorArray, 2, any)],
    kernel = kernel.candidates[apply(indicatorArray, 1, any)],
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 5,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 200,
    n.iter = 10000,
    burninRatio = 0.1,
    stepSizeFactor = 0.1,
    filllevel = filllevel.candidates[apply(indicatorArray, 4, any)]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  if(config$ndis > 801) return(NULL)
  
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster
  outDir <- with(config, paste0(baseDir, loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))
  
  allf <- list.files(outDir)
  allf <- allf[grep("priorTempered", allf)]
  allf.priorTemperedPhase2 <- allf[grep("priorTemperedPhase2\\.pdf", allf)]
  seed.priorTemperedPhase2 <- gsub(".*-([0-9]+)-.*", "\\1", allf.priorTemperedPhase2)
  seed.priorTemperedPhase2 <- head(sort(as.numeric(seed.priorTemperedPhase2)), 3)
  seed.priorTemperedPhase2 <- as.character(seed.priorTemperedPhase2)
  allpdff <- allf[grep("pdf", allf)]
  allpdff <- allpdff[gsub(".*-([0-9]+)-.*", "\\1", allpdff) %in% seed.priorTemperedPhase2]
  
  allf <- allf[grep("rds", allf)]
  ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), 
               simplify = "array")
  label <- dimnames(ci)[[3]]
  label <- strsplit(label, "-")
  label <- sapply(label, function(x) x[2])
  
  resultSummaryEach <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                                 identity, simplify = "array")
  sizeSummaryEach <- table(label)
  
  return(list(
    resultSummaryEach,
    sizeSummaryEach,
    file.path(outDir, allpdff)
  ))
}, mc.cores = 12)

resultSummary <- lapply(allSummaryList, function(x) x[[1]])
sizeSummary <- lapply(allSummaryList, function(x) x[[2]])
files2zip <- unlist(lapply(allSummaryList, function(x) x[[3]]))

save.image("largeExperimentSummary-priorTempered.rda")
cat(files2zip, file="file_list-priorTempered.txt", sep = "\n")
system("tar -czv -T file_list-priorTempered.txt -f pdfVisuals-priorTempered.tar.gz")
# load("../results/2017-11-28/largeExperimentSummary.rda")

# consolidate to previous rda file --------------------------------------------
library(abind)
load("../results/2017-12-07/largeExperimentSummary-priorTempered.rda")
resultSummaryPriorTempered <- resultSummary
sizeSummaryPriorTempered <- sizeSummary
load("../results/2017-11-28/largeExperimentSummary.rda")
for(arg in 1:length(indicatorArray)){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 3, any)],
    noise = noise.candidates[apply(indicatorArray, 2, any)],
    kernel = kernel.candidates[apply(indicatorArray, 1, any)],
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 5,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 200,
    n.iter = 10000,
    burninRatio = 0.1,
    stepSizeFactor = 0.1,
    filllevel = filllevel.candidates[apply(indicatorArray, 4, any)]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  if(config$ndis > 801) next
  
  resultSummary[[arg]] <- abind(resultSummary[[arg]], resultSummaryPriorTempered[[arg]])
  sizeSummary[[arg]] <- c(sizeSummary[[arg]], sizeSummaryPriorTempered[[arg]])
}
save.image("../results/2017-12-07/largeExperimentSummary.rda")

# see if coverage of df 2.5 is universally better ---------------------------
load("../results/2017-12-07/largeExperimentSummary.rda")

# in terms of coverage
# exclude trueMu as it is unfair
# exclude priorTemperedPhase2 as I only ran three repetitions

# use sample counting table
comparisonMaternDf <- sapply(noise.candidates, function(noise){
  sapply(nobs.candidates, function(nobs){
    out <- organizeOutput(maternDf = 2.01, noise = noise, nobs = nobs)$coverage["coverage",,c(-5),] -
      organizeOutput(maternDf = 2.5, noise = noise, nobs = nobs)$coverage["coverage",,c(-5),]
    apply((out > 0.05) - (out < -0.05), 2, mean)
    apply(out, 2, mean)
    # apply(out[,,(dim(out)[3]-1):dim(out)[3]], 2:3, mean)
  }, simplify = "array")
}, simplify = "array")

dim(comparisonMaternDf)
dimnames(comparisonMaternDf)[[2]] <- nobs.candidates
dimnames(comparisonMaternDf)[[3]] <- noise.candidates
mean(comparisonMaternDf)
apply(comparisonMaternDf, 1, mean)
apply(comparisonMaternDf, 2, mean)
apply(comparisonMaternDf, 3, mean)
# df = 2.5 is only better in the very low observation case
# df = 2.01 is particularly good for low noise high observation case high discretization


# in terms of bias
comparisonMaternDf <- sapply(noise.candidates, function(noise){
  sapply(nobs.candidates, function(nobs){
    out <- 
      abs(organizeOutput(maternDf = 2.01, noise = noise, nobs = nobs)$coverage["mean",,c(-5),]-c(0.2,0.2,3)) -
      abs(organizeOutput(maternDf = 2.5, noise = noise, nobs = nobs)$coverage["mean",,c(-5),]-c(0.2,0.2,3))
    pmin(1,pmax(-1, apply(out, 2, mean)))
    # apply(out[,,(dim(out)[3]-1):dim(out)[3]], 2:3, mean)
  }, simplify = "array")
}, simplify = "array")

dim(comparisonMaternDf)
dimnames(comparisonMaternDf)[[2]] <- nobs.candidates
dimnames(comparisonMaternDf)[[3]] <- noise.candidates
mean(comparisonMaternDf)
apply(comparisonMaternDf, 1, mean)
apply(comparisonMaternDf, 2, mean)
apply(comparisonMaternDf, 3, mean)
# df = 2.01 is particularly good for low noise high observation case high discretization

# see if coverage of prior tempered phase 2 is universally better --------------


# in terms of coverage
comparisonMaternDf <- sapply(noise.candidates, function(noise){
  sapply(nobs.candidates, function(nobs){
    out <- organizeOutput(maternDf = 2.01, noise = noise, nobs = nobs)$coverage["coverage",,"priorTemperedPhase2",] -
      organizeOutput(maternDf = 2.01, noise = noise, nobs = nobs)$coverage["coverage",,"phase2",]
    mean(out)
  }, simplify = "array")
}, simplify = "array")
dimnames(comparisonMaternDf)[[1]] <- nobs.candidates
dimnames(comparisonMaternDf)[[2]] <- noise.candidates

comparisonMaternDf
# coverage is almost universally better

# in terms of bias
comparisonMaternDf <- sapply(noise.candidates, function(noise){
  sapply(nobs.candidates, function(nobs){
    coverageInfo <- organizeOutput(maternDf = 2.01, noise = noise, nobs = nobs)$coverage
    out <- abs(coverageInfo["mean",,"priorTemperedPhase2",] - c(0.2, 0.2, 3)) -
      abs(coverageInfo["mean",,"phase2",] - c(0.2, 0.2, 3))
    mean(out)
  }, simplify = "array")
}, simplify = "array")
dimnames(comparisonMaternDf)[[1]] <- nobs.candidates
dimnames(comparisonMaternDf)[[2]] <- noise.candidates

comparisonMaternDf
# prior Tempered phase 2 bias is smaller at reasonable range of noise / nobs

# summarize priorTempered hes1 model --------------------------------
library(colorout)
library(gpds)
library(parallel)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:4

indicatorArray <- array(FALSE, dim=c(length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))

resultSummary <- list()
sizeSummary <- list()
files2zip <- c()

allSummaryList <- mclapply(1:length(indicatorArray), function(arg){
  indicatorArray[] <- FALSE
  indicatorArray[arg] <- TRUE
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 2, any)],
    noise = noise.candidates[apply(indicatorArray, 1, any)] * c(4, 1, 8),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 200,
    n.iter = 10000,
    burninRatio = 0.2,
    stepSizeFactor = 0.1,
    filllevel = filllevel.candidates[apply(indicatorArray, 3, any)],
    modelName = "hes1"
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  
  if(config$ndis > 801) return(NULL)
  
  if(grepl("/n/",getwd())){
    baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
  }else{
    baseDir <- "~/Workspace/DynamicSys/results/"  
  }
  outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                                "-ndis",ndis,"/"))
  
  allf <- list.files(outDir)
  
  seedOutputPdf <- allf[grep("trueMu\\.pdf", allf)]
  seedOutputPdf <- gsub(".*-([0-9]+)-.*", "\\1", seedOutputPdf)
  seedOutputPdf <- head(sort(as.numeric(seedOutputPdf)), 3)
  seedOutputPdf <- as.character(seedOutputPdf)
  allpdff <- allf[grep("pdf", allf)]
  allpdff <- allpdff[gsub(".*-([0-9]+)-.*", "\\1", allpdff) %in% seedOutputPdf]
  
  allf <- allf[grep("rds", allf)]
  ci <- sapply(allf, function(f) {
    tryCatch(readRDS(file.path(outDir, f)),
             error = function(e) NULL)
  }, simplify = "array")
  if(class(ci) == "list"){
    ci <- sapply(ci[-which(sapply(ci, is.null))], identity, simplify = "array")
  }
  label <- dimnames(ci)[[3]]
  label <- strsplit(label, "-")
  label <- sapply(label, function(x) x[2])
  
  resultSummaryEach <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                              identity, simplify = "array")
  sizeSummaryEach <- table(label)
  
  return(list(
    resultSummaryEach,
    sizeSummaryEach,
    file.path(outDir, allpdff)
  ))
}, mc.cores = 12)

resultSummary <- lapply(allSummaryList, function(x) x[[1]])
sizeSummary <- lapply(allSummaryList, function(x) x[[2]])
files2zip <- unlist(lapply(allSummaryList, function(x) x[[3]]))

save.image("hes1-largeExperimentSummary.rda")
cat(files2zip, file="hes1-file_list.txt", sep = "\n")
system("tar -czv -T hes1-file_list.txt -f pdfVisuals-hes1.tar.gz")
# load("../results/2017-11-28/largeExperimentSummary.rda")
