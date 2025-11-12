\name{DT_legendre}
\alias{DT_legendre}
\docType{data}
\title{
Simulated data for random regression
}
\description{
A data frame with 4 columns; SUBJECT, X, Xf and Y to show how to use the Legendre polynomials in the lmeb function using a numeric variable X and a response variable Y. 

}
\usage{data("DT_legendre")}
\format{
  The format is:
 chr "DT_legendre"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_legendre)
DT <- DT_legendre
head(DT)

\donttest{

########################################
############### sommer #################
########################################
if(requireNamespace("sommer")){
library(sommer)
library(orthopolynom)
mRR2<-mmes(Y~ 1 + Xf
           , random=~ vsm(usm(leg(X,1)),ism(SUBJECT))
           , rcov=~units
           , data=DT)
summary(mRR2)$varcomp

}

##############################################
############### lme4breeding #################
##############################################
if(requireNamespace("lme4breeding")){
library(lme4breeding)
library(orthopolynom)
Z <- with(DT, smm(leg(X,1)) )
for(i in 1:ncol(Z)){DT[,colnames(Z)[i]] <- Z[,i]}
## diagonal random regression Y ~ Xf + (0+leg0+leg1|| SUBJECT)
## unstructured random regression Y ~ Xf + (0+leg0+leg1| SUBJECT)
mRR2b<-lmeb(Y ~ Xf + (0+leg0+leg1| SUBJECT),
                , data=DT)
vc <- VarCorr(mRR2b); print(vc,comp=c("Variance"))
sigma(mRR2b)^2 # error variance

BLUP <- ranef(mRR2b, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

}

}

}
\keyword{datasets}
