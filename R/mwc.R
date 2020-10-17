#' Compute multiple wavelet coherence
#'
#' @param y time series 1 in matrix format (\code{m} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param x1 time series 2 in matrix format (\code{m} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param x2 time series 3 in matrix format (\code{m} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param pad pad the values will with zeros to increase the speed of the
#'   transform. Default is TRUE.
#' @param dj spacing between successive scales. Default is 1/12.
#' @param s0 smallest scale of the wavelet. Default is \code{2*dt}.
#' @param J1 number of scales - 1.
#' @param max.scale maximum scale. Computed automatically if left unspecified.
#' @param mother type of mother wavelet function to use. Can be set to
#'   \code{morlet}, \code{dog}, or \code{paul}. Default is \code{morlet}.
#'   Significance testing is only available for \code{morlet} wavelet.
#' @param param nondimensional parameter specific to the wavelet function.
#' @param lag1 vector containing the AR(1) coefficient of each time series.
#' @param sig.level significance level. Default is \code{0.95}.
#' @param sig.test type of significance test. If set to 0, use a regular
#'   \eqn{\chi^2} test. If set to 1, then perform a time-average test. If set to
#'   2, then do a scale-average test.
#' @param nrands number of Monte Carlo randomizations. Default is 300.
#' @param quiet Do not display progress bar. Default is \code{FALSE}
#'
#' @return Return a \code{vectorwavelet} object containing:
#' \item{coi}{matrix containg cone of influence}
#' \item{rsq}{matrix of wavelet coherence}
#' \item{phase}{matrix of phases}
#' \item{period}{vector of periods}
#' \item{scale}{vector of scales}
#' \item{dt}{length of a time step}
#' \item{t}{vector of times}
#' \item{xaxis}{vector of values used to plot xaxis}
#' \item{s0}{smallest scale of the wavelet }
#' \item{dj}{spacing between successive scales}
#' \item{mother}{mother wavelet used}
#' \item{type}{type of \code{vectorwavelet} object created (\code{mwc})}
#' \item{signif}{matrix containg \code{sig.level} percentiles of wavelet coherence
#'               based on the Monte Carlo AR(1) time series}
#'
#' @author Tunc Oygur (info@tuncoygur.com.tr)
#'
#' Code based on MWC MATLAB package written by Eric K. W. Ng and Johnny C. L. Chan.
#'
#' @references
#' T. Oygur, G. Unal.. Vector wavelet coherence for multiple time series.
#' \emph{Int. J. Dynam. Control (2020).}
#'
#' T. Oygur, G. Unal. 2017. The large fluctuations of the stock
#' return and financial crises evidence from Turkey: using wavelet
#' coherency and VARMA modeling to forecast stock return.
#' \emph{Fluctuation and Noise Letters}
#'
#' Ng, Eric KW and Chan, Johnny CL. 2012. Geophysical applications of partial
#' wavelet coherence and multiple wavelet coherence. \emph{Journal of Atmospheric
#' and Oceanic Technology} 29-12:1845--1853.
#'
#' @examples
#' old.par <- par(no.readonly=TRUE)
#'
#' t <- (-100:100)
#'
#' y <- sin(t*2*pi)+sin(t*2*pi/4)+sin(t*2*pi/8)+sin(t*2*pi/16)+sin(t*2*pi/32)+sin(t*2*pi/64)
#' x1 <- sin(t*2*pi/8)
#' x2 <- sin(t*2*pi/32)
#'
#' y <- cbind(t,y)
#' x1 <- cbind(t,x1)
#' x2 <- cbind(t,x2)
#'
#' ## Multiple wavelet coherence
#' result <- mwc(y, x1, x2, nrands = 10)
#' \donttest{
#' result <- mwc(y, x1, x2)
#' }
#'
#' ## Plot wavelet coherence and make room to the right for the color bar
#' ## Note: plot function can be used instead of plot.vectorwavelet
#' par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1,  pin = c(3,3))
#' plot.vectorwavelet(result, plot.cb = TRUE, main = "Plot multiple wavelet coherence")
#'
#' par(old.par)
#'
#' @keywords wavelet
#' @keywords coherence
#' @keywords continuous wavelet transform
#' @keywords multiple wavelet coherence
#'
#' @importFrom stats sd
#' @importFrom biwavelet check.data wt smooth.wavelet wtc.sig
#' @export
mwc <- function (y, x1, x2, pad = TRUE, dj = 1/12, s0 = 2 * dt, J1 = NULL,
                  max.scale = NULL, mother = "morlet", param = -1, lag1 = NULL,
                  sig.level = 0.95, sig.test = 0, nrands = 300, quiet = FALSE) {

  mother <- match.arg(tolower(mother), c("morlet", "paul", "dog"))

  # Check data format
  checked <- check.data(y = y, x1 = x1, x2 = x2)
  xaxis <- y[, 1]
  dt <- checked$y$dt

  t <- checked$y$t
  n <- checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt # automatic maxscale
    }
    J1 <- round(log2(max.scale/s0)/dj)
  }

  # Get AR(1) coefficients for each time series
  if (is.null(lag1)) {

    y.ar1 <- ar1nv(y[,2])$g
    x1.ar1 <- ar1nv(x1[,2])$g
    x2.ar1 <- ar1nv(x2[,2])$g

    lag1 <- c(y.ar1, x1.ar1, x2.ar1)
  }

  # Get CWT of each time series
  wt.y <- wt(d = y, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale,
             mother = mother, param = param, sig.level = sig.level,
             sig.test = sig.test, lag1 = lag1[1])
  wt.x1 <- wt(d = x1, pad = pad, dj = dj, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, param = param,
              sig.level = sig.level, sig.test = sig.test, lag1 = lag1[2])
  wt.x2 <- wt(d = x2, pad = pad, dj = dj, s0 = s0, J1 = J1,
              max.scale = max.scale, mother = mother, param = param,
              sig.level = sig.level, sig.test = sig.test, lag1 = lag1[3])

  # Standard deviation for each time series
  y.sigma <- sd(y[, 2], na.rm = TRUE)
  x1.sigma <- sd(x1[, 2], na.rm = TRUE)
  x2.sigma <- sd(x2[, 2], na.rm = TRUE)

  s.inv <- 1/t(wt.y$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt.y$wave))

  smooth.wt_y <- smooth.wavelet(s.inv*(abs(wt.y$wave)^2), dt, dj, wt.y$scale)
  smooth.wt_x1 <- smooth.wavelet(s.inv*(abs(wt.x1$wave)^2), dt, dj, wt.x1$scale)
  smooth.wt_x2 <- smooth.wavelet(s.inv*(abs(wt.x2$wave)^2), dt, dj, wt.x2$scale)

  coi <- pmin(wt.y$coi, wt.x1$coi, wt.x2$coi, na.rm = T)

  # Cross-wavelet computation
  cw.yx1 <- wt.y$wave * Conj(wt.x1$wave)
  cw.yx2 <- wt.y$wave * Conj(wt.x2$wave)
  cw.x2x1 <- wt.x2$wave * Conj(wt.x1$wave)

  smooth.cw_yx1 <- smooth.wavelet(s.inv*(cw.yx1), dt, dj, wt.y$scale)
  smooth.cw_yx2 <- smooth.wavelet(s.inv*(cw.yx2), dt, dj, wt.y$scale)
  smooth.cw_x2x1 <- smooth.wavelet(s.inv*(cw.x2x1), dt, dj, wt.y$scale)

  rsq.yx1 <- abs(smooth.cw_yx1)^2/(smooth.wt_y * smooth.wt_x1)
  rsq.yx2 <- abs(smooth.cw_yx2)^2/(smooth.wt_y * smooth.wt_x2)
  rsq.x2x1 <- abs(smooth.cw_x2x1)^2/(smooth.wt_x2 * smooth.wt_x1)

  r.yx1 <- sqrt(rsq.yx1)
  r.yx2 <- sqrt(rsq.yx2)
  r.x2x1 <- sqrt(rsq.x2x1)

  # Wavelet coherence
  norm.rsq <- (1 - r.x2x1^2)
  rsq <- (r.yx1^2+r.yx2^2-2*Re(r.yx1*Conj(r.yx2)*Conj(r.x2x1)))/norm.rsq

  # Phase difference
  phase <- atan2(Im(cw.yx1), Re(cw.yx1))

  if (nrands > 0) {
    signif <- wtc.sig(nrands = nrands, lag1 = lag1,
                      dt = dt, n, pad = pad, dj = dj, J1 = J1, s0 = s0,
                      max.scale = max.scale, mother = mother, sig.level = sig.level,
                      quiet = quiet)
  }
  else {
    signif <- NA
  }

  results <- list(coi = coi,
                  rsq = rsq,
                  phase = phase,
                  period = wt.y$period,
                  scale = wt.y$scale,
                  dt = dt,
                  t = t,
                  xaxis = xaxis,
                  s0 = s0,
                  dj = dj,
                  mother = mother,
                  type = "mwc",
                  signif = signif)

  class(results) <- "vectorwavelet"
  return(results)

}
