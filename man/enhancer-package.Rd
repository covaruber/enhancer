\name{enhancer-package}
\alias{enhancer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

Mixed Effect Model Enhancing Functions

}
\description{

enhancer is a collection of functions useful in the fitting of linear mixed effect models
and general data manipulation. As more of my packages use the same special functions I
decided to have the special functions as its own package to facilitate the maintenance.

}

\section{Functions for genetic analysis}{

The package provides kernels to the estimate additive (\code{\link{A.matr}})
relationship matrix for diploid and polyploid organisms. A good converter from 
letter code to numeric format is implemented in the function \code{\link{atcg1234}} 
and \code{\link{atcg1234BackTransform}}, 
which supports higher ploidy levels than diploid. Additional functions for genetic 
analysis have been included such as build a genotypic hybrid marker matrix 
(\code{\link{build.HMM}}). If you need to use pedigree you need to convert your 
pedigree into a relationship matrix (use the `getA` function from the pedigreemm package).

Also the \code{\link{LD.decay}} and the \code{\link{neMarker}} functions can be used
to calculate linkage disequilibrium and effective population size.

}

\section{Functions for trial analysis}{

Recently, spatial modeling has been added added to enhancer using the two-dimensional
spline \code{\link{tps}} function. Based on \code{\link{bbasis}}

}

\section{Functions for special mixed models}{

Available are \code{\link{add.diallel.vars}}, \code{\link{leg}},
\code{\link{overlay}}, \code{\link{rrm}}, \code{\link{redmm}},
\code{\link{simGECorMat}} function.

}

\section{Functions for data manipulation and others}{

Available are \code{\link{adiag1}}, \code{\link{imputev}}, \code{\link{propMissing}}
\code{\link{replaceValues}}, \code{\link{stan}}, \code{\link{logspace}}
function.

Simple coloring functions available are \code{\link{jet.colors}}, \code{\link{bathy.colors}},
\code{\link{transp}}

}

\section{Keeping enhancer updated}{

The enhancer package is updated on CRAN every 4-months due to CRAN policies 
but you can find the latest source at https://github.com/covaruber/enhancer. 
This can be easily installed typing the following in the R console:

library(devtools)

install_github("covaruber/enhancer")

This is recommended if you reported a bug, was fixed and was immediately pushed 
to GitHub but not in CRAN until the next update.
}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com 
or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable amount of
time developing this package. Please cite this package in your publication. 
Type 'citation("enhancer")' to know how to cite it.


}
\keyword{ R package }