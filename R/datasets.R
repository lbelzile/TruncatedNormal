#' Latent membranous Lupus Nephritis dataset
#'
#' The data represents two clinical measurements (covariates), which are used to predict the occurrence of latent membranous lupus nephritis. The dataset consists of measurements on 55 patients of which 18 have been diagnosed with latent membranous lupus.
#'  
#' @name lupus
#' @docType data
#' @source The data were transcribed from Table 1, page 22, of Dyk and Meng (2001).
#' @references D. A. van Dyk and X.-L. Meng (2001) \emph{The art of data augmentation (with discussion)}. Journal of Computational and  Graphical Statistics, volume 10, pages 1-50.
#' @format a data frame with columns "response", "const", "x1" and "x2"
#' @seealso The dataset is used in the examples of \code{\link{mvrandn}}
NULL

#' Women wage dataset from Mroz (1987)
#'
#' The data are from the Panel Study of Income Dynamics (PSID) longitudinal study, 1976 wave.
#' They give the number of work hours of married women along with socio-economic variables and the number of children.
#'  
#' @name mroz
#' @docType data
#' @source W. Greene's website, accessed 17.12.2019 at <http://www.stern.nyu.edu/~wgreene/Text/Edition7/TableF5-1.csv>.
#' @references T. A. Mroz, 1987. \emph{The Sensitivity of an Empirical Model of Married Women's Hours of Work to Economic and Statistical Assumptions}, Econometrica, \bold{55}(4), pp. 765-799
#' @seealso \code{\link[carData]{Mroz}}
#' @format a data frame containing the following variables:
#' \itemize{
#' \item \code{whrs}: hours of work
#' \item \code{kidslt6}: number of children aged 5 and below years old in household
#' \item \code{kidsge6}: number of children between age of 6 and 18 in household
#' \item \code{age}: age (in years)
#' \item \code{educ}: number of years in school
#' \item \code{hearn}: hourly earnings
#' \item \code{exp}: years of previous labor market experience
#' }
NULL
