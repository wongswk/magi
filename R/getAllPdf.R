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
  allf <- head(sort(allf[grep("pdf", allf)]), 9)
  
  files2zip <- c(files2zip, file.path(outDir, allf))
  print(arg)
}
cat(files2zip, file="file_list.txt", sep = "\n")
system(paste("tar -czvf resultPdf.tar.gz ", files2zip))
