\name{replaceValues}
\alias{replaceValues}

\title{
Replace the values of a vector
}
\description{
Simple function to map a vector of values to another.
}
\usage{
replaceValues(Source, Search, Replace)
}

\arguments{
  \item{Source}{
A source of values to change.
}

  \item{Search}{
The values to search.
}

  \item{Replace}{
The equivalent values to replace.
}

}
\details{

Simple function to map a vector of values another.

}
\value{
\describe{
\item{$res}{new values from a second vector.}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}

\examples{


x <- letters[1:5]
y <- 1:5
z <- rep(letters[1:5],3)

replaceValues(Source=z, Search=x, Replace=y)


}
