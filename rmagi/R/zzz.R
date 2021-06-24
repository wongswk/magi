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
               "##  See https://www.overleaf.com/9732586bkctqjywbxny for additional documentation.\n",
               "##  COMPILING_INFORMATION_HERE\n",
               "##  Please cite software as:\n",
               "##   our citation\n",
               "##   goes here\n",
               "##   Journal of Statistical Software, XX(X): X-XX. \n##\n",
               sep = "")
  packageStartupMessage(foo)
}
