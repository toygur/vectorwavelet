#' Compute n-dimensional vector wavelet coherence
#'
#' @param y time series y in matrix format (\code{m} rows x 2 columns). The
#'   first column should contain the time steps and the second column should
#'   contain the values.
#' @param x multivariate time series x in matrix format (\code{m} rows x n columns).
#'   The first column should contain the time steps and the other columns should
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
#' \item{type}{type of \code{vectorwavelet} object created (\code{vwc})}
#' \item{signif}{matrix containg \code{sig.level} percentiles of wavelet coherence
#'               based on the Monte Carlo AR(1) time series}
#'
#' @author Tunc Oygur (info@tuncoygur.com.tr)
#'
#' @references
#' T. Oygur, G. Unal.. Vector wavelet coherence for multiple time series.
#' \emph{Int. J. Dynam. Control (2020).}
#'
#' @examples
#' old.par <- par(no.readonly=TRUE)
#'
#' t <- (-100:100)
#'
#' y <- sin(t*2*pi)+sin(t*2*pi/4)+sin(t*2*pi/8)+sin(t*2*pi/16)+sin(t*2*pi/32)+sin(t*2*pi/64)
#' x1 <- sin(t*2*pi/8)
#' x2 <- sin(t*2*pi/16)
#' x3 <- sin(t*2*pi/32)
#' x4 <- sin(t*2*pi/64)
#'
#' y <- cbind(t,y)
#' x <- cbind(t,x1,x2,x3,x4)
#'
#' ## n-dimensional multiple wavelet coherence
#' result <- vwc(y, x, nrands = 10)
#' \donttest{
#' result <- vwc(y, x)
#' }
#'
#' ## Plot wavelet coherence and make room to the right for the color bar
#' ## Note: plot function can be used instead of plot.vectorwavelet
#' par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1,  pin = c(3,3))
#' plot.vectorwavelet(result, plot.cb = TRUE, main = "Plot n-dimensional vwc (n=5)")
#'
#' par(old.par)
#'
#' @keywords wavelet
#' @keywords coherence
#' @keywords continuous wavelet transform
#' @keywords n-dimensional wavelet coherence
#'
#' @importFrom stats sd
#' @importFrom biwavelet check.data wt smooth.wavelet wtc.sig
#' @export
vwc <- function (y, x, pad = TRUE, dj = 1/12, s0 = 2 * dt, J1 = NULL,
                  max.scale = NULL, mother = "morlet", param = -1, lag1 = NULL,
                  sig.level = 0.95, sig.test = 0, nrands = 300, quiet = FALSE) {

  mother <- match.arg(tolower(mother), c("morlet", "paul", "dog"))

  # Check data format
  checked <- n.check.data(y = y, x = x)
  xaxis <- y[, 1]
  dt <- checked$y$dt

  t <- checked$y$t
  n <- checked$y$n.obs

  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale <- (n * 0.17) * 2 * dt
    }
    J1 <- round(log2(max.scale/s0)/dj)
  }

  # Get AR(1) coefficients
  if (is.null(lag1)) {
    y.ar1 <- ar1nv(y[, 2])$g
    x.ar1 <- as.numeric(apply(x[,-1], 2, function(x) {ar1nv(x)$g }))

    lag1 <- c(y.ar1, x.ar1)
  }

  # Standard deviation
  y.sigma <- sd(y[, 2], na.rm = TRUE)
  x.sigma <- sd(x[, 2], na.rm = TRUE)

  df <- cbind(y,x[,-1])
  colnames(df) <- c("t","y",paste0("x",1:(ncol(x)-1)))

  n_dim <- ncol(df)-1

  # Get CWT of each time series
  wt.res <- list()
  for(i in 2:ncol(df)) {
    xi <- df[,c(1,i)]
    temp.col <- colnames(df)[i]
    wt.res[[temp.col]] <- wt(d = xi, pad = pad, dj = dj, s0 = s0, J1 = J1, max.scale = max.scale, mother = mother,
                             param = param, sig.level = sig.level,sig.test = sig.test, lag1 = lag1[i-1])
    rm(xi);rm(temp.col)
  }
  rm(i)

  s.inv <- 1/t(wt.res[["y"]]$scale)
  s.inv <- matrix(rep(s.inv, n), nrow = NROW(wt.res[["y"]]$wave))

  s.wt.res <- list()
  for(i in 2:ncol(df)) {
    temp.col <- colnames(df)[i]
    s.wt.res[[temp.col]] <- smooth.wavelet(s.inv*(abs(wt.res[[temp.col]]$wave)^2),  dt, dj, wt.res[[temp.col]]$scale)
    rm(temp.col)
  }
  rm(i)

  coi.res <- matrix(NA, nrow = n, ncol = (ncol(df)-1))
  for(i in 2:ncol(df)) {
    temp.col <- colnames(df)[i]
    coi.res[,(i-1)] <- wt.res[[temp.col]]$coi
    rm(temp.col)
  }
  rm(i)
  coi <- apply(coi.res, 1, function(x) min(x, na.rm = T))
  rm(coi.res)

  # Cross-wavelet computation
  cw <- list()
  smooth.cw <- list()
  rsq <- list()
  r <- list()
  for(i in 2:(ncol(df)-1)){
    for(j in (i+1):ncol(df)){
      temp.col.i <- colnames(df)[i]
      temp.col.j <- colnames(df)[j]
      temp.cw <- wt.res[[temp.col.i]]$wave * Conj(wt.res[[temp.col.j]]$wave)
      temp.smooth.cw <- smooth.wavelet(s.inv*(temp.cw), dt, dj, wt.res[["y"]]$scale)
      temp.rsq <- abs(temp.smooth.cw)^2/(s.wt.res[[temp.col.i]] * s.wt.res[[temp.col.j]])
      temp.r <- sqrt(temp.rsq)

      cw[[paste0(temp.col.i,"-",temp.col.j)]] <- temp.cw
      smooth.cw[[paste0(temp.col.i,"-",temp.col.j)]] <- temp.smooth.cw
      rsq[[paste0(temp.col.i,"-",temp.col.j)]] <- temp.rsq
      r[[paste0(temp.col.i,"-",temp.col.j)]] <- temp.r

      rm(temp.col.i);rm(temp.col.j);rm(temp.cw);rm(temp.smooth.cw);rm(temp.rsq);rm(temp.r)
    }
    rm(j)
  }
  rm(i)

  ######################################################################################################################
  ######################################################################################################################
  # Wavelet coherence
  r_tilde <- list()
  for(i in 1:(ncol(df)-1)) {
    for(j in 1:(ncol(df)-1)) {
      temp.col.i <- colnames(df)[i+1]
      temp.col.j <- colnames(df)[j+1]
      if(i==j) {
        r_tilde[[paste0(temp.col.i,"-",temp.col.j)]] <- 1
      } else if (j > i) {
        r_tilde[[paste0(temp.col.i,"-",temp.col.j)]] <- r[[paste0(temp.col.i,"-",temp.col.j)]]
      } else {
        r_tilde[[paste0(temp.col.i,"-",temp.col.j)]] <- Conj(r[[paste0(temp.col.j,"-",temp.col.i)]])
      }
    }
    rm(j)
  }
  rm(i)

  ##Cofactor function ##################################################################################################
  cofactor.wavelogy <- function(ii, jj, del_ii, del_jj, level, n_dim,order) {

    order.sign <- if(order %% 2 == 1) 1 else -1
    temp.col.i <- colnames(df)[ii+1]
    temp.col.j <- colnames(df)[jj+1]

    res <- order.sign * r_tilde[[paste0(temp.col.i,"-",temp.col.j)]]

    ii_vec <- setdiff(c(1:n_dim),del_ii)
    jj_vec <- setdiff(c(1:n_dim),del_jj)

    if(n_dim-level != 2) {

      temp.res <- 0
      for(kk in 1:(n_dim-level)) {
        temp.res <- temp.res + cofactor.wavelogy(ii=min(ii_vec), jj=jj_vec[kk], del_ii=c(del_ii,min(ii_vec)), del_jj=c(del_jj,jj_vec[kk]), level=level+1, n_dim, order=kk)
      }
      res <- res * temp.res
    } else {

      temp.col.i_1 <-  colnames(df)[ii_vec[1]+1]
      temp.col.j_1 <-  colnames(df)[jj_vec[1]+1]

      temp.col.i_2 <-  colnames(df)[ii_vec[2]+1]
      temp.col.j_2 <-  colnames(df)[jj_vec[2]+1]

      res <- res * (r_tilde[[paste0(temp.col.i_1,"-",temp.col.j_1)]]*r_tilde[[paste0(temp.col.i_2,"-",temp.col.j_2)]] -
                      r_tilde[[paste0(temp.col.i_2,"-",temp.col.j_1)]]*r_tilde[[paste0(temp.col.i_1,"-",temp.col.j_2)]])
    }
    return(res)
  }

  ##Cxx ####################################################################################################
  Cxx <- 0
  for(k in 1:n_dim) {
    Cxx <- Cxx + cofactor.wavelogy(ii=1, jj=k, del_ii=1, del_jj=k, level=1, n_dim, order=k)
  }

  C11 <- cofactor.wavelogy(ii=1, jj=1, del_ii=1, del_jj=1, level=1, n_dim, order=1)

  ##Rsq ####################################################################################################
  rsq <- 1 - (Cxx / C11)

  ######################################################################################################################
  ######################################################################################################################
  # Phase difference
  phase <- atan2(Im(cw[["y-x1"]]), Re(cw[["y-x1"]]))

  if (nrands > 0) {
    signif <- wtc.sig(nrands = nrands, lag1 = c(y.ar1, x.ar1),
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
                  period = wt.res[["y"]]$period,
                  scale = wt.res[["y"]]$scale,
                  dt = dt,
                  t = t,
                  xaxis = xaxis,
                  s0 = s0,
                  dj = dj,
                  mother = mother,
                  type = "vwc",
                  signif = signif)

  class(results) <- "vectorwavelet"
  return(results)
}
