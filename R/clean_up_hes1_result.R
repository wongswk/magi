# Summarize the results
library(gpds)

# remove results that don't have common seed
subdirs <- c(
  "../results/for_paper/7param//variablephi-notemper",
  "../results/for_paper/7param//variablephi-temper-warmstart",
  "../results/for_paper/7param//variablephi-temper-warmstart-updatephi",
  "../results/for_paper/7param//ramsay",
  "../results/for_paper/fullobs//variablephi-notemper", 
  "../results/for_paper/fullobs//variablephi-temper-coldstart", 
  "../results/for_paper/fullobs//variablephi-temper-warmstart",
  "../results/for_paper/fullobs//ramsay"
)
all_files <- lapply(subdirs, list.files)
all_files <- lapply(all_files, function(x) x[grep(".rda", x)])
sapply(all_files, length)
all_seeds <- lapply(all_files, function(x) gsub(".*log-([0-9]+)-.*rda", "\\1", x))
common_seeds <- all_seeds[[1]]
for (i in 2:length(all_seeds)){
  common_seeds <- intersect(common_seeds, all_seeds[[i]])  
}

for(i in 1:length(subdirs)){
  to_delete <- all_files[[i]][rowSums(sapply(common_seeds, function(s) grepl(s, all_files[[i]]))) == 0]
  to_delete <- paste0(subdirs[i], "/", to_delete)
  for(each_d in to_delete){
    system(paste0("rm ", each_d))
  }
}

