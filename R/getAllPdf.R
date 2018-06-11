baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster
allDirs <- list.dirs(baseDir, recursive = FALSE)
allDirs <- allDirs[grep("FN-withmeanBand-generalMatern", allDirs)]

files2zip <- c()

for(outDir in allDirs){
  allf <- list.files(outDir)
  allf <- head(sort(allf[grep("pdf", allf)]), 6)
  
  files2zip <- c(files2zip, file.path(outDir, allf))
}
cat(files2zip, file="file_list.txt", sep = "\n")
system("tar -czvf resultPdf.tar.gz -T file_list.txt")
