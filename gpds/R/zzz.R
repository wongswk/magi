#' @useDynLib gpds
#' @import Rcpp 
#' @import microbenchmark
#' @import RcppArmadillo
NULL

.onAttach <- function( ... )
{
  # gpdsLib <- dirname(system.file(package = "gpds"))
  # pkgdes <- as.list(packageDescription("gpds", lib.loc = gpdsLib))
  # version <- pkgdes$Version
  # BuildDate <- pkgdes$Date

  foo <- paste(#"## \n##  gpds (Version ", version, ", Build Date: ", BuildDate, ")\n", 
               "##  See https://www.overleaf.com/9732586bkctqjywbxny for additional documentation.\n",
               "##  Please cite software as:\n",
               "##   our citation\n",
               "##   goes here\n",
               "##   Journal of Statistical Software, XX(X): X-XX. \n##\n",
               sep = "")
  packageStartupMessage(foo)
}
