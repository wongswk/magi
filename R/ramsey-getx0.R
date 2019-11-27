### Use this to extract fitted x0 from B-splines

getfdx0 <- function (x, y, Lfdobj = 0, href = TRUE, titles = NULL, xlim = NULL, 
          ylim = NULL, xlab = NULL, ylab = NULL, ask = FALSE, nx = NULL, 
          axes = NULL, ...) 
{
  fdobj <- x
  if (!(inherits(fdobj, "fd"))) 
    stop("First argument is not a functional data object.")
  {
    if (is.null(axes)) {
      if (is.null(fdobj$basis$axes)) {
        Axes <- TRUE
        axFun <- FALSE
      }
      else {
        if (!inherits(fdobj$basis$axes, "list")) 
          stop("fdobj$basis$axes must be a list;  ", 
               "class(fdobj$basis$axes) = ", class(fdobj$basis$axes))
        if (!(inherits(fdobj$basis$axes[[1]], "character") || 
              inherits(fdobj$basis$axes[[1]], "function"))) 
          stop("fdobj$basis$axes[[1]] must be either a function or the ", 
               "name of a function;  class(fdobj$basis$axes[[1]]) = ", 
               class(fdobj$basis$axes[[1]]))
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(fdobj$basis$axes, ...)
      }
    }
    else {
      if (is.logical(axes)) {
        Axes <- axes
        axFun <- FALSE
      }
      else {
        if (!inherits(axes, "list")) 
          stop("axes must be a logical or a list;  class(axes) = ", 
               class(axes))
        if (!(inherits(axes[[1]], "character") || inherits(axes[[1]], 
                                                           "function"))) 
          stop("axes[[1]] must be either a function or the ", 
               "name of a function;  class(axes[[1]]) = ", 
               class(axes[[1]]))
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(axes, ...)
      }
    }
    }
  Lfdobj <- int2Lfd(Lfdobj)
  if (!inherits(Lfdobj, "Lfd")) 
    stop("Second argument is not a linear differential operator.")
  coef <- fdobj$coefs
  coefd <- dim(coef)
  ndim <- length(coefd)
  #show(ndim)
  nbasis <- coefd[1]
  if (is.null(nx)) 
    nx <- max(c(501, 10 * nbasis + 1))
  nrep <- coefd[2]
  if (ndim > 2) 
    nvar <- coefd[3]
  else nvar <- 1
  basisobj <- fdobj$basis
  rangex <- basisobj$rangeval
  if (missing(y)) {
    y <- nx
  }
  else {
    if (is.numeric(y)) 
      y <- as.vector(y)
  }
  Y <- y
  if (length(y) == 1) {
    if (y >= 1) {
      y <- seq(rangex[1], rangex[2], len = round(y))
    }
    else {
      stop("'y' a single number less than one.")
    }
  }
  if (min(y) < rangex[1] || max(y) > rangex[2]) 
    stop("Values in Y are outside the basis range.")
  if (is.null(xlim)) {
    xlim <- rangex
  }
  else {
    rangex[1] <- max(rangex[1], xlim[1])
    rangex[2] <- min(rangex[2], xlim[2])
    if (length(Y) == 1) 
      y <- seq(rangex[1], rangex[2], len = round(Y))
  }
  fdmat <- eval.fd(y, fdobj, Lfdobj)
  rangey <- range(fdmat, na.rm = TRUE)
  if (is.null(ylim)) 
    ylim <- rangey
  fdnames = fdobj$fdnames
  fdlabelslist = fdlabels(fdnames, nrep, nvar)
  xlabel = fdlabelslist$xlabel
  ylabel = fdlabelslist$ylabel
  casenames = fdlabelslist$casenames
  varnames = fdlabelslist$varnames
  if (is.null(xlab)) 
    xlab <- xlabel
  if (is.null(ylab)) 
    ylab <- ylabel
  if (ndim < 2) {
    plot(y, fdmat, type = "l", xlim = xlim, ylim = ylim, 
         xlab = xlab, ylab = ylab, axes = Axes, ...)
    if (axFun) 
      do.call(axList[[1]], axList[-1])
    if (zerofind(fdmat) && href) 
      abline(h = 0, lty = 2)
  }
  if (ndim == 2) {
    if (!ask) {
      
      return(fdmat[1,]) #### these are the intial values!
      
      matplot(y, fdmat, type = "l", xlim = xlim, ylim = ylim, 
              xlab = xlab, ylab = ylab, axes = Axes, ...)
      if (axFun) 
        do.call(axList[[1]], axList[-1])
      if (zerofind(fdmat) && href) 
        abline(h = 0, lty = 2)
    }
    else {
      op <- par(ask = FALSE)
      on.exit(par(op))
      cat("Multiple plots:  Click in the plot to advance to the next")
      for (irep in 1:nrep) {
        plot(y, fdmat[, irep], type = "l", xlim = xlim, 
             ylim = ylim, xlab = xlab, ylab = ylab, axes = Axes, 
             ...)
        if (axFun) 
          do.call(axList[[1]], axList[-1])
        if (irep < 2) 
          par(ask = ask)
        if (!is.null(casenames)) 
          title(casenames[irep])
        else title(paste("Case", irep))
        if (zerofind(ylim) && href) 
          abline(h = 0, lty = 2)
      }
    }
  }
  if (ndim == 3) {
    if (!ask) {
      for (ivar in 1:nvar) {
        matplot(y, fdmat[, , ivar], type = "l", xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, ask = FALSE, 
                axes = Axes, ...)
        if (axFun) 
          do.call(axList[[1]], axList[-1])
        if (!is.null(varnames)) 
          title(varnames[ivar])
        else title(paste("Variable", ivar))
        if (zerofind(ylim) && href) 
          abline(h = 0, lty = 2)
      }
    }
    else {
      op <- par(ask = FALSE)
      on.exit(par(op))
      cat("Multiple plots:  Click in the plot to advance to the next")
      for (irep in 1:nrep) {
        for (ivar in 1:nvar) {
          plot(y, fdmat[, irep, ivar], type = "l", xlim = xlim, 
               ylim = ylim, xlab = xlab, ylab = ylab, axes = Axes, 
               ...)
          if (axFun) 
            do.call(axList[[1]], axList[-1])
          if (!is.null(casenames)) 
            titlestr = casenames[irep]
          else titlestr = paste("Case", irep)
          if (!is.null(varnames)) {
            titlestr = paste(titlestr, "  ", varnames[ivar])
          }
          else {
            titlestr = paste(titlestr, "  ", "Variable", 
                             ivar)
          }
          title(titlestr)
          if (zerofind(ylim) && href) 
            abline(h = 0, lty = 2)
        }
      }
    }
  }
  "done"
}