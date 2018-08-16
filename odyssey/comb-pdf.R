file_lists <- list.files(".")
file_lists <- file_lists[grep("_kernelmat.pdf", file_lists)]
filename_lists <- gsub("_kernelmat.pdf", "", file_lists)
for (filename in filename_lists){
  # system(paste0("pdftk ",
  #               filename, "_m-likelihoodmove-stan.pdf ",
  #               " cat 8 output ",
  #               filename, "_m-likelihoodtrue-stan.pdf "
  # ))
  
  system(paste0("pdftk ", filename, "_gpds.pdf ", 
                filename, "_kernelmat.pdf ",
                " cat output ",
                filename, "_comb.pdf"
  ))
  
}
