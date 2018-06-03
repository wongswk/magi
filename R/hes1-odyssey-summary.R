library(gpds)

outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-log-withmeanBand-generalMatern-nobs33-noise0.15_0.15_0.15-ndis33/"
outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-log-withmeanBand-generalMatern-nobs33-noise0.15_0.15_0.1-ndis33/"
outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-withmeanBand-generalMatern-nobs33-noise0.8_0.4_1-ndis33/"

allf <- list.files(outDir)
allf <- allf[grep("\\.rds", allf)]
# allf <- allf[grep("fully-observed", allf)]
allf <- allf[grep("partial-observed", allf)]

ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), simplify = "array")
label <- dimnames(ci)[[3]]
label <- strsplit(label, "-generalMatern-")
label <- sapply(label, function(x) x[1])
label <- unique(label)

output <- round(t(apply(ci, 1:2, mean)), 4)
rownames(output) <- letters[1:7]

print(label)
print(output)
xtable::xtable(output, caption = label,  digits = 4)
