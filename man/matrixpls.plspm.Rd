\name{matrixpls.plspm}
\alias{matrixpls.plspm}
\title{A plspm compatibility wrapper for matrixpls}
\usage{
  matrixpls.plspm(Data, inner_matrix, outer_list,
    modes = NULL, scheme = "centroid", scaled = TRUE,
    tol = 1e-05, iter = 100, boot.val = FALSE, br = NULL,
    plsr = FALSE, dataset = TRUE)
}
\arguments{
  \item{Data}{A numeric matrix or data frame containing the
  manifest variables.}

  \item{inner_matrix}{A square (lower triangular) boolean
  matrix representing the inner model (i.e. the path
  relationships betwenn latent variables).}

  \item{outer_list}{List of vectors with column indices
  from \code{x} indicating the sets of manifest variables
  asociated to the latent variables (i.e. which manifest
  variables correspond to the latent variables). Length of
  \code{outer_list} must be equal to the number of rows of
  \code{inner_matrix}.}

  \item{modes}{A character vector indicating the type of
  measurement for each latent variable. \code{"A"} for
  reflective measurement or \code{"B"} for formative
  measurement (\code{NULL} by default). The length of
  \code{modes} must be equal to the length of
  \code{outer_list}).}

  \item{scheme}{A string of characters indicating the type
  of inner weighting scheme. Possible values are
  \code{"centroid"}, \code{"factor"}, or \code{"path"}.}

  \item{scaled}{A logical value indicating whether scaling
  data is performed When (\code{TRUE} data is scaled to
  standardized values (mean=0 and variance=1) The variance
  is calculated dividing by \code{N} instead of
  \code{N-1}).}

  \item{tol}{Decimal value indicating the tolerance
  criterion for the iterations (\code{tol=0.00001}). Can be
  specified between 0 and 0.001.}

  \item{iter}{An integer indicating the maximum number of
  iterations (\code{iter=100} by default). The minimum
  value of \code{iter} is 100.}

  \item{boot.val}{A logical value indicating whether
  bootstrap validation is performed (\code{FALSE} by
  default).}

  \item{br}{An integer indicating the number bootstrap
  resamples. Used only when \code{boot.val=TRUE}. When
  \code{boot.val=TRUE}, the default number of re-samples is
  100, but it can be specified in a range from 100 to
  1000.}

  \item{plsr}{A logical value indicating whether pls
  regression is applied to calculate path coefficients
  (\code{FALSE} by default).}

  \item{dataset}{A logical value indicating whether the
  data matrix should be retrieved (\code{TRUE} by
  default).}
}
\value{
  An object of class \code{\link[plspm]{plspm}}.
}
\description{
  \code{matrixpls.plspm} mimics \code{\link[plspm]{plspm}}
  function of the \code{plspm} package. The arguments and
  their default values and the output of the function are
  identical with \code{\link[plspm]{plspm}} function, but
  internally the function uses matrixpls estimation.
}
\details{
  The function \code{matrixpls.plspm} calculates indicator
  weights and estimates a model identically to the
  \code{\link[plspm]{plspm}} function. In contrast to the
  \code{\link{matrixpls}} function that provides only
  weights and parameter estimates, this function also
  reports multiple post-estimation statistics. Because of
  this \code{matrixpls.plspm} is substantially less
  efficient than the \code{\link{matrixpls}} function.

  The argument \code{inner_matrix} is a matrix of zeros and
  ones that indicates the structural relationships between
  composites. This must be a lower triangular matrix.
  \code{inner_matrix} will contain a 1 when column \code{j}
  affects row \code{i}, 0 otherwise. \cr
}
\examples{
library(plspm)

# Run the example from plspm package

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
# outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply plspm
plspm.res <- plspm(satisfaction, sat_inner, sat_outer, sat_mod,
                   scaled=FALSE, boot.val=FALSE)

# apply matrixpls
matrixpls.res <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod,
                                 scaled=FALSE, boot.val=FALSE)

# print the resuls

summary(plspm.res)
summary(matrixpls.res)

# If RUnit is installed check that the results are identical

if(is.element("RUnit", installed.packages()[,1])){
	library(RUnit)
	checkEquals(plspm.res, matrixpls.res, tol = 0.001)
}
}
\references{
  Sanchez, G. (2013). \emph{PLS Path Modeling with R.}
  Retrieved from http://www.gastonsanchez.com/PLS Path
  Modeling with R.pdf
}
\seealso{
  \code{\link[plspm]{plspm}}
}

