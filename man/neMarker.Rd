\name{neMarker}
\alias{neMarker}
\title{Effective population size based on marker matrix}
\description{
  `neMarker` uses a marker matrix to approximate the effective population size (Ne) by discovering how many individuals are needed to sample all possible alleles in a population.
}
\usage{
  neMarker(M, neExplore=NULL, maxMarker=1000, nSamples=5)
}

\arguments{
  \item{M}{marker matrix coded in a numerical faashion (any allele dosage is fine).}
  \item{neExplore}{a vector of numbers with the effective population sizes to be explored.}
  \item{maxMarker}{maximum number of markers to use for the analysis.}
  \item{nSamples}{number of individuals to sample for the Ne calculation.}
}
\value{
  \describe{
    
    \item{$S3}{ A vector with allele coverage based on different number of individuals  }
  }
}
\references{
  
  Not based on any theory published yet but in a solid intuition on what is really important for a breeding program when we ask what is the effective population size
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
nInds=300
nMarks=1000

#random population of nInds lines with nMarks markers
M <- matrix(rep(0,nInds*nMarks),nInds,nMarks)
for (i in 1:nInds) {
  M[i,] <- ifelse(runif(nMarks)<0.5,-1,1)
}

# run the function
ne <- neMarker(M, neExplore = seq(2,300,50), nSamples = 10)
plot(ne$Ne)

}