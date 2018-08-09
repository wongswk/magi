phase8 = list(
  list(mu = curCov[[1]]$mu, dotmu = curCov[[1]]$dotmu),
  list(mu = curCov[[2]]$mu, dotmu = curCov[[2]]$dotmu)
)
saveRDS(phase8, paste0(outDir, "phase8-mu.rds"))


curCov[[1]]$mu = phase8[[1]]$mu
curCov[[1]]$dotmu = phase8[[1]]$dotmu
curCov[[2]]$mu = phase8[[2]]$mu
curCov[[2]]$dotmu = phase8[[2]]$dotmu
