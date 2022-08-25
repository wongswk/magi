#' @useDynLib magi
#' @import Rcpp 
NULL

.onAttach <- function( ... )
{
  # magiLib <- dirname(system.file(package = "magi"))
  # pkgdes <- as.list(packageDescription("magi", lib.loc = magiLib))
  # version <- pkgdes$Version
  # BuildDate <- pkgdes$Date

  foo <- paste(#"## \n##  magi (Version ", version, ", Build Date: ", BuildDate, ")\n",
               "##  See vignette https://arxiv.org/abs/2203.06066 for detailed usage instructions, \n",
               "##  and https://github.com/wongswk/magi for additional documentation.\n",
               "##  \n",
               "##  Please cite MAGI method as:\n",
               "##  Yang, S., Wong, S.W. and Kou, S.C., 2021.\n",
               "##  Inference of dynamic systems from noisy and sparse data via manifold-constrained Gaussian processes.\n",
               "##  Proceedings of the National Academy of Sciences, 118(15), e2020397118.\n",
               sep = "")
  packageStartupMessage(foo)
}
