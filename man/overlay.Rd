\name{overlay}
\alias{overlay}
\title{Overlay Matrix}
\description{
  `overlay` adds r times the design matrix for model term t to the existing design matrix. Specifically, if the model up to this point has p effects and t has a effects, the a columns of the design matrix for t are multiplied by the scalar r (default value 1.0). This can be used to force a correlation of 1 between two terms as in a diallel analysis.
}
\usage{
  overlay(..., rlist=NULL, prefix=NULL, sparse=FALSE)
}

\arguments{
  \item{...}{as many vectors as desired to overlay.}
  \item{rlist}{a list of scalar values indicating the times that each incidence matrix overlayed should be multiplied by. By default r=1.}
  \item{prefix}{a character name to be added before the column names of the final overlay matrix. This may be useful if you have entries with names starting with numbers which programs such as asreml doesn't like, or for posterior extraction of parameters, that way 'grep'ing is easier.}
  \item{sparse}{a TRUE/FALSE statement specifying if the matrices should be built as sparse or regular matrices.}
}
\value{
  \describe{
    
    \item{$S3}{ an incidence matrix with as many columns levels in the vectors provided to build the incidence matrix.}
  }
}
\references{
  Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.
  
Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
femalef <- as.factor(c("a","b","c"))
malef <- as.factor(c("d","e","a"))
overlay(femalef,malef, sparse = TRUE)
overlay(femalef,malef, sparse = FALSE)

}
