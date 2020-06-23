all_files = list.files('./')
all_files = all_files[grep("txt", all_files)]

run_time <- lapply(all_files, function(f) read.table(f, sep=" ", header = FALSE))

run_time <- do.call(rbind, run_time)
summary(run_time[,3] / run_time[,1])
summary(run_time[,2] / run_time[,1])
