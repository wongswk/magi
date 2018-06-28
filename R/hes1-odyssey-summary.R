library(gpds)

outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-log-withmeanBand-generalMatern-nobs33-noise0.15_0.15_0.15-ndis33/"
outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-log-withmeanBand-generalMatern-nobs33-noise0.15_0.15_0.1-ndis33/"
outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-withmeanBand-generalMatern-nobs33-noise0.8_0.4_1-ndis33/"
outDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/Hes1-withmeanBand-generalMatern-nobs33-noise0.8_0.4_1.5-ndis33/"

allf <- list.files(outDir)
allf <- allf[grep("\\.rds", allf)]
# allf <- allf[grep("fully-observed", allf)]
# allf <- allf[grep("partial-observed", allf)]

ci <- sapply(allf, function(f) readRDS(file.path(outDir, f)), simplify = "array")
label <- dimnames(ci)[[3]]
label <- strsplit(label, "-generalMatern-")
label <- sapply(label, function(x) x[1])

resultSummary <- sapply(tapply(1:length(label), label, function(id) apply(ci[,,id], 1:2, mean)),
                               identity, simplify = "array")
sizeSummary <- table(label)

label_numeric <- gsub(".*-phase([0-9]+).*", "\\1", names(sizeSummary))
label_numeric <- as.numeric(label_numeric)
label_numeric[is.na(label_numeric)] <- 1

sizeSummary[order(label_numeric)]
resultSummary[,,order(label_numeric)]
outDir <- substr(outDir, 1, nchar(outDir)-1)
system(paste0("mv ", outDir, " ", outDir, format(Sys.time(), "%Y-%m-%d_%H-%M-%S")))


output <- round(t(apply(ci, 1:2, mean)), 4)
rownames(output) <- letters[1:7]

print(label)
print(output)
xtable::xtable(output, caption = label,  digits = 4)
