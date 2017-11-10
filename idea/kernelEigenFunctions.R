# convergence of eigen value and eigen vector ----------------------------------
maxT <- 10
tAll <- seq(0, maxT, length=201)
eigenList <- lapply(10:201, function(ndis){
  tDis <- seq(0, maxT, length=ndis)
  signedDist <- outer(tDis, tDis, '-')
  gpcov <- calCov(c(1,1), abs(signedDist), -sign(signedDist))
  outX <- apply(gpcov$CeigenVec[,1:10], 2, function(y) {
    y <- y * sqrt(ndis / maxT)
    outy <- approx(tDis, y, tAll)$y
    outy <- outy*sign(outy[1])
    outy
  })
  list(head(1/gpcov$Ceigen1over / ndis, 10), outX)
})
matplot(t(sapply(eigenList, function(x) x[[1]])), type="l")
eigenFun <- sapply(eigenList, function(x) x[[2]], simplify = "array")

mycolor <- rev(gray.colors(dim(eigenFun)[3]))
matplot(eigenFun[,1,], type="l", lty=1, col=mycolor)
matplot(eigenFun[,2,], type="l", lty=1, col=mycolor)
matplot(eigenFun[,3,], type="l", lty=1, col=mycolor)

Kfunc <- function(r) calCovMatern(c(1,1), abs(r), NULL, complexity = 0)$C
plot.function(Kfunc, from = -5, to = 5, n=1e4)

