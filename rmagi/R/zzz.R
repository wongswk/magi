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
               "##  See https://doi.org/10.18637/jss.v109.i04 for a user guide, \n",
               "##  https://doi.org/10.1073/pnas.2020397118 for details of the MAGI method, \n",
               "##  and https://github.com/wongswk/magi for additional documentation.\n",
               "##  \n",
               '##  To cite the magi software, please see citation("magi")\n',
               sep = "")
  packageStartupMessage(foo)
}
