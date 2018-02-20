getFourierFrequency <- function(x, showplot=FALSE){
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
    return(which.max(zmodEffective))
  }else{
    outliers <- outliers[outliers > median(zmodEffective)]
    whichOutliers <- which(zmodEffective %in% outliers)
    return(max(whichOutliers))
  }
}

N <- 100
Tall <- 2*pi
x1 <- sin(Tall*(1:N)/N)
x2 <- sin(2*Tall*(1:N)/N)
x3 <- sin(5*Tall*(1:N)/N)
plot(x1)
plot(x2)
plot(x3)
x <- (x1+x2+x3+rnorm(N)*0.5)
plot(x)
z <- fft(x)
barplot(Mod(z))

