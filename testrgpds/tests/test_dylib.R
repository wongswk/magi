library(testrgpds)
# segfault
testrgpds:::phisigllikRcpp(c(1, 0.5, 2),
                           as.matrix(1:10),
                           as.matrix(dist(1:10)),
                           "matern")
