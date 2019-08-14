#' calculate frequency of data series and give recommanded bandwidth parameter phi2
#' 
#' phi2 rule of thumb should be half of the periodicity
#' 
#' 
#' @export
getFrequencyBasedPrior <- function(x, showplot=FALSE){
  z <- fft(x)
  zmod <- Mod(z)
  zmodEffective <- zmod[-1]
  zmodEffective <- zmodEffective[1:(length(zmodEffective)/2)]
  if(showplot) {
    names(zmodEffective) <- 1:length(zmodEffective)
    barplot(zmodEffective)
    boxplot(zmodEffective)
  }
  outliers <- boxplot.stats(zmodEffective)$out 
  if(length(outliers)==0){
    freq <- which.max(zmodEffective)
  }else{
    outliers <- outliers[outliers > median(zmodEffective)]
    whichOutliers <- which(zmodEffective %in% outliers)
    freq <- max(whichOutliers)
  }
  meanFactor <- 0.5 / freq
  sdFactor <- (1 - meanFactor) / 3
  c(meanFactor=meanFactor, sdFactor=sdFactor)
}
