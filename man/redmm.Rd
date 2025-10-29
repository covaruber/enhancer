\name{redmm}
\alias{redmm}
\title{Reduced Model Matrix}
\description{
  `redmm` reduces a model matrix by performing a singular value decomposition or Cholesky on an incidence matrix.
}
\usage{
  redmm(x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE)
}

\arguments{
  \item{x}{as vector with values to form a model matrix or the complete incidence matrix itself for an effect of interest.}
  \item{M}{an optional matrix of features explaining the levels of x. If not provided is assumed that the entire incidence matrix has been provided in x. But if provided, the decomposition occurs in the matrix M. }
  \item{Lam}{a matrix of loadings in case is already available to avoid recomputing it.}
  \item{nPC}{number of principal components to keep from the matrix of loadings to form the model matrix.}
  \item{cholD}{should a Cholesky or a Singular value decomposition should be used. The default is the SVD.}
  \item{returnLam}{should the function return the loading matrix in addition to the incidence matrix. Default is FALSE.}
}
\value{
  \describe{
    
    \item{$S3}{ A list with 3 elements:
    
    1) The model matrix to be used in the mixed modeling.
    
    2) The reduced matrix of loadings (nPC columns).
    
    3) The full matrix of loadings.
    
    }
  }
}
\references{
  
Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
# data set including the covariate to be reduced
DT <- data.frame(Env=sort(rep(paste("e",1:5,sep="_"), 2)))
DT
# GxE matrix
H0 <- matrix(rnorm(30), nrow=5, ncol=6) # 5 envs, 6 features
rownames(H0) <- paste("e",1:5,sep="_")
colnames(H0) <- paste("x",1:6,sep="_")
H0
# get compression
Z <- with(DT,  redmm(Env, M = H0, nPC = 3))
Z


}
