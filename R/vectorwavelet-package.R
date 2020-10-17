#' @docType package
#' @name vectorwavelet-package
#' @aliases vectorwavelet
#' @exportPattern ^[[:alpha:]]+
#' @importFrom stats as.dist filter qchisq quantile rnorm sd ts var weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @title
#' Vector wavelet coherence for multiple time series
#'
#' @description
#' Description: This package can be used to perform multiple wavelet coherence (mwc),
#'  quadruple wavelet coherence (qmwc), and n-dimensional vector wavelet coherence (vwc) analyses.
#'
#' @author Tunc Oygur, Gazanfer Unal
#'
#' Maintainer: Tunc Oygur <info@tuncoygur.com.tr>
#'
#' Code based on biwavelet package written by Tarik C. Gouhier, Aslak Grinsted, Viliam Simko.
#'
#' @references
#' T. Oygur, G. Unal.. Vector wavelet coherence for multiple time series.
#' \emph{Int. J. Dynam. Control (2020).}
#'
#' T. Oygur, G. Unal.. The large fluctuations of the stock
#' return and financial crises evidence from Turkey: using wavelet
#' coherency and VARMA modeling to forecast stock return.
#' \emph{Fluctuation and Noise Letters, 2017}
#'
#' T.C. Gouhier, A. Grinstead and V. Simko. 2016. \emph{biwavelet:
#' Conduct univariate and bivariate wavelet analyses (Version 0.20.15).}
#' Available from http://github.com/tgouhier/biwavelet
#'
#' Ng, Eric KW and Chan, Johnny CL. 2012. Geophysical applications of partial
#' wavelet coherence and multiple wavelet coherence. \emph{Journal of Atmospheric
#' and Oceanic Technology} 29-12:1845--1853.
#'
#' Grinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross
#' wavelet transform and wavelet coherence to geophysical time series.
#' \emph{Nonlinear Processes in Geophysics} 11:561-566.
#'
#' Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis.
#' \emph{Bulletin of the American Meteorological Society} 79:61-78.
#'
#' @keywords wavelet
#' @keywords coherence
#' @keywords continuous wavelet transform
#' @keywords multiple wavelet coherence
#' @keywords quadruple wavelet coherence
#' @keywords vector wavelet coherence
NULL

.onAttach <- function(libname, pkgname) {

  # just to show a startup message
  message <- paste("vectorwavelet", utils::packageVersion("vectorwavelet"), "loaded.")
  packageStartupMessage(message, appendLF = TRUE)
}
