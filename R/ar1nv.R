#' AR1NV - Estimate the parameters for an AR(1) model
#'
#' @param x One dimensional time series vector
#'
#' @return Return a list containing:
#' \item{g}{estimate of the lag-one autocorrelation.}
#' \item{a}{estimate of the noise variance.}
#'
#' @author Tunc Oygur (info@tuncoygur.com.tr)
#'
#' Code based on a cross wavelet and wavelet coherence toolbox MATLAB package written by Eric Breitenberger
#'
#' @references
#' SGrinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross
#' wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics} 11:561-566.
#'
ar1nv <- function(x){

  N <- length(x)
  m <- mean(x)
  x <- x-m

  # Lag zero and one covariance estimates:
  c0 <- t(x) %*% (x/N)
  c1 <- x[1:N-1] %*% (x[2:N]/(N-1))
  g <- c1/c0
  a <- sqrt((1-g^2)*c0)

  return(list(g=as.numeric(g), a=as.numeric(a)))
}
