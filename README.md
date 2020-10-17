<!-- README.md is generated from README.Rmd. Please edit that file -->

vectorwavelet R package
=======================

[![Support
badge](https://img.shields.io/badge/support-vectorwavelet-green.svg)](http://stackoverflow.com/questions/tagged/vectorwavelet)

Download and Install
--------------------

To download the development version of the package, type the following
at the R command line:

``` r
install.packages("devtools")
devtools::install_github("toygur/vectorwavelet")
```

About vectorwavelet
-------------------

The vectorwavelet R package was created for n-dimension wavelet
coherence analysis for multivariate time series. Biwavelet R package was
used for univariate and bivariate analysis and components.

Besides, multiple wavelet coherence (mwc) in the package is the R
transformation of Matlab codes created by Ng, Eric KW, and Chan.

Quadruple (qmwc) and n-dimensional vector wavelet coherence (vwc) are
used for multidimensional wavelet analysis with a vectorwavelet package.

How to cite
-----------

``` r
citation("vectorwavelet")
#> 
#> To cite package 'vectorwavelet' in publications use:
#> 
#>   Tunc Oygur and Gazanfer Unal (2020). vectorwavelet: Vector wavelet
#>   coherence for multiple time series. R package version 0.1.0.
#>   https://github.com/toygur/vectorwavelet
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {vectorwavelet: Vector wavelet coherence for multiple time series},
#>     author = {Tunc Oygur and Gazanfer Unal},
#>     year = {2020},
#>     note = {R package version 0.1.0},
#>     url = {https://github.com/toygur/vectorwavelet},
#>   }
#> 
#> ATTENTION: This citation information has been auto-generated from the
#> package DESCRIPTION file and may need manual editing, see
#> 'help("citation")'.
```

References
----------

<a id="Oygur2020"/> T. Oygur, G. Unal. 2020. **Vector wavelet coherence
for multiple time series.**. Int. J. Dynam. Control.

<a id="Oygur2017"/> T. Oygur, G. Unal. 2017. **The large fluctuations of
the stock return and financial crises evidence from Turkey: using
wavelet coherency and VARMA modeling to forecast stock return.**.
Fluctuation and Noise Letters 16/02

<a id="biwavelet"/> T.C. Gouhier, A. Grinstead and V. Simko. 2016.
**biwavelet: Conduct univariate and bivariate wavelet analyses**.
Available from github.com/tgouhier/biwavelet.

<a id="Eric2012"/> Ng, Eric KW and Chan, Johnny CL. 2012. **Geophysical
applications of partial wavelet coherence and multiple wavelet
coherence.**. Journal of Atmospheric and Oceanic Technology
29-12:1845–1853.

<a id="Grinsted2004"/> Grinsted, A., J. C. Moore, and S. Jevrejeva.
2004. **Application of the cross wavelet transform and wavelet coherence
to geophysical time series**. Nonlinear Processes in Geophysics
11:561–566.

<a id="TorrenceCompo1998"/> Torrence, C., and G. P. Compo. 1998. **A
practical guide to wavelet analysis**. Bulletin of the American
Meteorological Society 79:61–78.

<a id="TorrenceWebster1998"/> Torrence, C., and P. J. Webster. 1998.
**The annual cycle of persistence in the El Niño/Southern Oscillation**.
Quarterly Journal of the Royal Meteorological Society 124:1985–2004.
