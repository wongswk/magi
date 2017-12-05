library(gpds)

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
  
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster
  outDir <- with(config, paste0(baseDir, loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))
  
  allf <- list.files(outDir)
  allpdff <- head(sort(allf[grep("pdf", allf)]), 9)
  allf <- allf[grep("rds", allf)]
  ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), 
               simplify = "array")
  label <- dimnames(ci)[[3]]
  label <- strsplit(label, "-")
  label <- sapply(label, function(x) x[2])
  
  resultSummary[[arg]] <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                          identity, simplify = "array")
  sizeSummary[[arg]] <- table(label)
  files2zip <- c(files2zip, file.path(outDir, allpdff))
}
save.image("largeExperimentSummary.rda")
cat(files2zip, file="file_list.txt", sep = "\n")
system("tar -czv -T file_list.txt -f pdfVisuals.tar.gz")
# load("../results/2017-11-28/largeExperimentSummary.rda")

length(sizeSummary)
length(which(sapply(sizeSummary, is.null)))
length(indicatorArray)
summary(do.call(rbind, sizeSummary))
dim(do.call(rbind, sizeSummary))

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

organizeOutput(maternDf = 2.01, noise = 0.1, nobs = 51)

organizeOutput(maternDf = 2.5, noise = 0.1, nobs = 51)

organizeOutput(maternDf = 2.01, noise = 0.1, nobs = 201)

organizeOutput(maternDf = 2.01, noise = 0.1, nobs = 401) 

organizeOutput(maternDf = 2.01, noise = 0.5, nobs = 201)

organizeOutput(maternDf = 2.01, noise = 1.0, nobs = 201)

# summarize refilled priorTempered first phase --------------------------------
library(gpds)

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

for(arg in 270:length(indicatorArray)){
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
  
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster
  outDir <- with(config, paste0(baseDir, loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))
  
  allf <- list.files(outDir)
  allf <- allf[grep("priorTempered", allf)]
  allpdff <- head(sort(allf[grep("pdf", allf)]), 3)
  allf <- allf[grep("rds", allf)]
  ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), 
               simplify = "array")
  label <- dimnames(ci)[[3]]
  label <- strsplit(label, "-")
  label <- sapply(label, function(x) x[2])
  
  resultSummary[[arg]] <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                                 identity, simplify = "array")
  sizeSummary[[arg]] <- table(label)
  files2zip <- c(files2zip, file.path(outDir, allpdff))
  print(arg)
}
save.image("largeExperimentSummary-priorTempered.rda")
cat(files2zip, file="file_list-priorTempered.txt", sep = "\n")
system("tar -czv -T file_list-priorTempered.txt -f pdfVisuals-priorTempered.tar.gz")
# load("../results/2017-11-28/largeExperimentSummary.rda")

# consolidate to previous rda file --------------------------------------------
library(abind)
load("../results/2017-12-05/largeExperimentSummary-priorTempered.rda")
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
save.image("../results/2017-12-05/largeExperimentSummary-added-priorTempered.rda")
