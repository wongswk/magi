# testhat helper functions

#' calculate frequency of data series and give recommanded bandwidth parameter phi2
#' 
#' phi2 rule of thumb should be half of the periodicity
#' 
#' 
#' @noRd
getFrequencyBasedPrior <- function(x, showplot=FALSE){
  z <- fft(x)
  zmod <- Mod(z)
  zmodEffective <- zmod[-1]
  zmodEffective <- zmodEffective[1:(length(zmodEffective)/2)]
  if(showplot) {
    plot(as.numeric(x), xlab="index", ylab="x")
    names(zmod) <- NULL
    plot(zmod, col=c(1,rep(2, (length(zmod)-1)/2),rep(1, (length(zmod)-1)/2)))
    title("modulus of fft", line=7)
    
    names(zmodEffective) <- 1:length(zmodEffective)
    
    upperQuarter = sort(zmodEffective)[ceiling(length(zmodEffective) * 0.75)]
    lowerQuarter = sort(zmodEffective)[floor(length(zmodEffective) * 0.25)]
    iqr = upperQuarter - lowerQuarter
    outliers = zmodEffective[zmodEffective > upperQuarter + 1.5 * iqr]
    if (length(outliers) == 0) {
      freq <- which.max(zmodEffective)
      abline(v=1+freq, col=4)
    } else {
      outliers <- outliers[outliers > median(zmodEffective)]
      whichOutliers <- which(zmodEffective %in% outliers)
      abline(v=1+whichOutliers, col=4)
      freq <- max(whichOutliers)
    }
    meanFactor <- 0.5/freq
    legend("top", c("range of frequencies with large loadings", 
                    "modulus of effective frequency loading"), 
           lty=c(1, NA), pch=c(NA, 1), col=c(4, 2))
    msg <- paste0("WANT: frequency with largest loading = ", which.max(zmodEffective), 
                  ", corresponding prior factor = ", round(0.5/which.max(zmodEffective), digits = 4))
    freq_wtd <- weighted.mean(whichOutliers, outliers^2)
    msg <- paste0(msg, "\nWANT: mod^2 weighted average frequency with large loadings = ", freq_wtd, ", corresponding prior factor = ", round(0.5/freq_wtd, digits = 4))
    freq_wtd <- weighted.mean(1:length(zmodEffective), zmodEffective^2)
    msg <- paste0(msg, "\nWANT: mod^2 weighted average frequency among all = ", freq_wtd, ", corresponding prior factor = ", round(0.5/freq_wtd, digits = 4))
    period_wtd <- weighted.mean(0.5/whichOutliers, outliers^2)
    msg <- paste0(msg, "\nWANT: mod^2 weighted average half periodicity (i.e. prior factor) with large loadings = ", period_wtd)
    period_wtd <- weighted.mean(0.5/(1:length(zmodEffective)), zmodEffective^2)
    msg <- paste0(msg, "\nWANT: mod^2 weighted average half periodicity (i.e. prior factor) among all = ", period_wtd)
    msg <- paste0(msg, "\nCURRENT: highest frequency with large loadings = ", freq, ", corresponding prior factor = ", round(0.5/freq, digits = 4))
    mtext(msg)
  }
  
  freq <- weighted.mean(1:length(zmodEffective), zmodEffective^2)
  meanFactor <- 0.5 / freq
  sdFactor <- (1 - meanFactor) / 3
  c(meanFactor=meanFactor, sdFactor=sdFactor)
}
