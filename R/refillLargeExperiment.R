library(gpds)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:4
kernel.candidates <- c("generalMatern", "matern") # basically df in matern

indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1

indicatorArray[arg] <- TRUE

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = noise.candidates[apply(indicatorArray, 2, any)],
  kernel = kernel.candidates[apply(indicatorArray, 1, any)],
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
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/"  
}
outDir <- with(config, paste0(baseDir, loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))

if(config$ndis <= 801){
  fileLists <- list.files(outDir)
  seedLists <- gsub(".*-([0-9]+)[-\\.].*", "\\1", fileLists)
  seedLists <- as.numeric(unique(seedLists))
  if(config$kernel == "matern") seedLists <- head(sort(seedLists), 3)
  for(oldSeed in seedLists){
    config$seed <- oldSeed
    
    rdsFile <- paste0(
      outDir,
      config$loglikflag,"-priorTemperedPhase2-",config$kernel,"-",config$seed,".rds")
    
    pdfFile <- paste0(
      outDir,
      config$kernel,"-",config$seed,"-priorTemperedPhase2.pdf")
    
    if(file.exists(rdsFile) && file.exists(pdfFile)){
      next
    }
    source("R/priorTempered-repeated-sample.R")
  }
}
