delta <- 0.01
df <- 2.5
x <- seq(1, 3, delta)
y <- besselK(x, df)
dy <- -(besselK(x, df-1) + besselK(x, df+1))/2
# plot(x, y, type="l")
# plot(x[-1], diff(y)/delta, type="l")
# lines(x, dy, col=2)
# lines(x, -y*df/x, col=3)
testthat::expect_equal((besselK(x, df-1) - besselK(x, df+1))/2, -y*df/x)
testthat::expect_equal((dy[-1]+dy[-length(dy)])/2, diff(y)/delta, tolerance=1e-3)
