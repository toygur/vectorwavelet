#' Check the format of multivariate time series
#'
#' @param y time series y in matrix format (\code{m} rows x 2 columns). The
#' first column should contain the time steps and the second column should
#' contain the values.
#' @param x multivariate time series x in matrix format (\code{m} rows x (1 + (n-1)) columns). The
#' first column should contain the time steps and the other columns should
#' contain the values.
#'
#' @return Returns a named list containing:
#' \item{t}{time steps}
#' \item{dt}{size of a time step}
#' \item{n.obs}{number of observations}
#'
#' @author Tunc Oygur (info@tuncoygur.com.tr)
#'
#' Code based on biwavelet package written by Tarik C. Gouhier.
#'
#' @examples
#' #Example 1:
#' t1 <- cbind(1:100, rnorm(100))
#' n.check.data(y = t1)
#'
#' #Example 2:
#' t1 <- cbind(1:100, rnorm(100))
#' t2 <- cbind(1:100, rnorm(100), rnorm(100), rnorm(100))
#' n.check.data(y = t1, x = t2)
#'
#'
#' @export
n.check.data <- function (y, x = NULL) {

  if (is.null(y)) {
    stop("The time series cannot be NULL")
  }

  y.check <- n.check.datum(y)

  if(!is.null(x)) {
    for(i in 2:ncol(x)) {

      xi <- x[,c(1,i)]
      xi.check <- n.check.datum(xi)

      if (any(diff(y[, 1]) != diff(xi[, 1]))) {
        stop("The time series must have the same step size")
      }
      if (y.check$n.obs != xi.check$n.obs) {
        stop("The time series must have the same length (see merge command)")
      }
    }
  }
  return(list(y = y.check))
}

#' Helper function
#' @param x matrix
#' @return list(t, dt, n.obs)
#' @note This function is not exported
n.check.datum <- function (x) {
  if (NCOL(x) > 1) {
    t <- x[, 1]
    diffs <- diff(t)
    dt <- as.numeric(diffs[1])
    epsilon <- 0.1 * dt
    if (any(abs(diff(t) - dt) > epsilon)) {
      stop("The step size must be constant ",
           "(see approx function to interpolate)")
    }
    else {
      if (class(t) == "Date") {
        t <- seq_len(NROW(t))
        dt <- diff(t)[1]
      }
    }
  }
  else {
    stop("Error: time steps have to be in column 1")
  }
  return(list(t = t, dt = dt, n.obs = NROW(x)))
}
