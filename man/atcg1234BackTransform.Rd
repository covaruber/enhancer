\name{atcg1234BackTransform}
\alias{atcg1234BackTransform}
\title{
Letter to number converter
}
\description{
This function was designed to help users back transform the numeric marker matrices from the function atcg1234 into letters. 

}
\usage{
atcg1234BackTransform(marks, refs)
}
\arguments{
  \item{marks}{a centered marker matrix coming from atcg1234.}
  \item{refs}{a 2 x m matrix for m markers (columns) and 2 rows where the reference and alternate alleles for each marker are indicated.}
  
}
\value{
\describe{

\item{markers}{ a new marker matrix leter coded according to the reference allele matrix.}

}
}

\author{
Giovanny Covarrubias-Pazaran
}

\examples{

genotypes <- rbind(
  c("AAAG"  ,   "AAGG"  ,   "TTCC"  ,   "TTTT"  ,   "CCCC"),
  c("AAGG"  ,  "AGGG"   ,  "TTCC"   ,  "TTTT"   ,  "TCCC"),
  c("AAGG"  ,   "GGGG"  ,   "TTCC"  ,   "TTTT"  ,   NA),
  c("AAAA"  ,   "GGGG"  ,   "CCCC"  ,   "TTTT"  ,   NA),
  c("AAAA"  ,   "GGGG"  ,   "CCCC"  ,   "TTCC"  ,   NA)
)
rownames(genotypes) <- letters[1:5]
colnames(genotypes) <- paste0("m",1:5)

## convert markers to numeric format polyploid potatoes
numo <- atcg1234(data=genotypes, ploidy=4)
numo

## convert markers back to letters
Xb <- atcg1234BackTransform(marks= numo$M, refs= numo$ref.alleles)
Xb


}