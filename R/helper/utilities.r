# UTILITY ROUTINES FOR WRITING UPDATE AND LPR FUNCTIONS.
#
# Radford M. Neal, 2012.


# CHECK WHETHER OR NOT A VALUE IS OUT OF BOUNDS.  The value and the lower
# and upper bounds are either vectors or lists.  The result if TRUE if the
# value is within bounds, FALSE if not.

value_out_of_bounds <- function (value, lower, upper)
{ 
  if (is.list(value))
  { out.of.bounds <- FALSE
    if (!is.null(lower))
    { for (e in names(lower))
      { if (!all(value[[e]]>=lower[[e]]))
        { out.of.bounds <- TRUE
          break
        }
      }
    }
    if (!out.of.bounds && !is.null(upper))
    { for (e in names(upper))
      { if (!all(value[[e]]<=upper[[e]])) 
        { out.of.bounds <- TRUE
          break
        }
      }
    }
    out.of.bounds
  }
  else
    !is.null(lower) && !all(value>=lower) ||
    !is.null(upper) && !all(value<=upper)
}


# PROCESS THE REP ARGUMENT OF AN UPDATE FUNCTION.  Checks that it is a
# numeric scalar, rounds it to an integer, and checks that it is at
# least one.  Returns the argument, which may have changed from rounding.

process_rep_argument <- function (rep)
{
  if (!is.numeric(rep) || length(rep)!=1 || !(rep <- round(rep)) >= 1)
  { stop("Invalid rep argument")
  }

  rep
}


# PROCESS THE STEP AND RAND.STEP ARGUMENTS OF AN UPDATE FUNCTION.  Arguments
# are the length of the (vector) state, the step argument, and the rand.step
# argument.

process_step_arguments <- function (n, step, rand.step)
{
  if (!is.numeric(step) || length(step)!=1 && length(step)!=n 
       || any(!is.finite(step)) || any(step<0))
  { stop("Invalid step argument")
  }

  if (!is.numeric(rand.step) 
       || length(rand.step)!=1 && length(rand.step)!=length(step) 
       || any(!is.finite(rand.step)) || any(rand.step<0))
  { stop("Invalid rand.step argument")
  }

  if (any(rand.step!=0))
  { 
    step <- step * exp (runif (length(rand.step), -rand.step, rand.step))
  }

  step
}


# PROCESS THE NSTEPS ARGUMENT OF AN UPDATE FUNCTION.  The argument should
# be the number of steps (eg, in a trajectory), which must be at least one.  
# Stops with an error message if it's invalid.

process_nsteps_argument <- function (nsteps)
{
  if (!is.numeric(nsteps) || length(nsteps)!=1 
   || nsteps!=floor(nsteps) || nsteps<1)
  { stop("Invalid nsteps argument")
  }
}


# PROCESS LOWER OR UPPER BOUNDS.  The arguments are the dimensionality and 
# a numeric scalar or vector giving bounds, which if not scalar must match the 
# dimensionality.  Stops with an error message if it's invalid, except that 
# NULL is also allowed.  Returns its argument, repeated to the dimensionality 
# if it's a scalar.

process_bounds <- function (n, bounds)
{
  if (is.null(bounds)) return (bounds)

  if (!is.numeric(bounds) || length(bounds)!=1 && length(bounds)!=n)
  { stop("Invalid bounds")
  }

  if (length(bounds)==n) bounds else rep(bounds,length=n)
}


# CHECK GRADIENT COMPUTATIONS BY FINITE DIFFERENCES.

check_grad <- function (lpr, value, eps)
{
  g <- attr (lpr(value,grad=TRUE), "grad")

  if (is.list(g))
  { f <- list()
    for (e in names(g))
    { f[[e]] <- numeric(length(g[[e]]))
      for (i in 1:length(g[[e]]))
      { vp <- vm <- value
        vp[[e]][i] <- vp[[e]][i] + eps
        vm[[e]][i] <- vm[[e]][i] - eps
        f[[e]][i] <- (lpr(vp,grad=TRUE) - lpr(vm,grad=TRUE)) / (2*eps)
      }
    }
    mapply(function(x,y)rbind((x),(y),x-y),g,f,SIMPLIFY=FALSE)
  }
  else
  { f <- numeric(length(g))
    for (i in 1:length(g))
    { vp <- vm <- value
      vp[i] <- vp[i] + eps
      vm[i] <- vm[i] - eps
      f[i] <- (lpr(vp,grad=TRUE) - lpr(vm,grad=TRUE)) / (2*eps)
    }
    rbind(g,f,g-f)
  }
}


# ADD NUMBERS REPRESENTED BY THEIR LOGS.  Computes the log of the sum of
# the exponentials of its arguments.  The arguments may be vectors or matrices,
# in which case this operation is carried out separately on each element.
# The arguments must all have the same dimensions, or be scalar.  The 
# computation is done in a manner that guarantees that the result is valid 
# even when directly computing the exponentials would result in overflow or 
# underflow.

log_add_exp <- function (...)
{
  m <- pmax(...)
  x <- 0

  for (a in list(...))
  { x <- x + exp(a-m)
  }

  m + log(x)
}


# SUM VECTORS OF NUMBERS REPRESENTED BY THEIR LOGS.  Computes the log of the 
# sum of the exponentials of all the elements in all its arguments.  The 
# computation is done in a manner that guarantees that the result is valid 
# even when directly computing the exponentials would result in overflow or 
# underflow.

log_sum_exp <- function (...)
{
  m <- max(...)

  m + log(sum(exp(c(...)-m)))
}


# COMPUTE A WEIGHTED AVERAGE OF VALUES REPRESENTED BY THEIR LOGS.  Takes a
# vector of weights and a vector of values of the same length, and returns
# the log of the weighted average of the exponentials of the values.  This
# is done in a manner that works even if exponentiating the values would
# result in overflow or underflow.  The weights must be non-negative and not
# all zero, but need not sum to one (they are normalized to sum to one here).

log_average_exp <- function (log.values, weights)
{ m <- max(log.values)
  m + log (sum (weights*exp(log.values-m)) / sum(weights))
}
