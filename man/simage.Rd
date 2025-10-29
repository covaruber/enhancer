\name{simage}
\alias{simage}
\title{Image of sparsity between two variables}
\description{\code{image} for sparsity.
}
\usage{
simage(data, Var1=NULL, Var2=NULL, \dots)
}
\arguments{
 \item{data}{model data of class \code{"data.frame"}}
 \item{Var1}{variable to set in the x axis}
 \item{Var2}{variable to set in the y axis}
 \item{\dots}{Further arguments to be passed to the image function.}
}
\value{
  vector of image
 }
\author{Giovanny Covarrubias}
\examples{

# row x column combinations existing

DT <- data.frame(Row=sort(rep(1:10,5)), Col=rep(1:5,10))
DT$Rowf <- as.factor(DT$Row)
DT$Colf <- as.factor(DT$Col)
DT <- DT[-1,] # remove one intentionally
simage(data=DT, Var1 = "Rowf", Var2 = "Colf")

# you can also inpect which cross combinations exist

}
\seealso{\code{\link{image}} }
\keyword{models}

