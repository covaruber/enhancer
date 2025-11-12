\name{DT_ige}
\alias{DT_ige}
\alias{A_ige}
\docType{data}
\title{
Data to fit indirect genetic effects.
}
\description{

This dataset contains phenotpic data for 98 individuals where they are measured with the purpose of identifying the effect of the neighbour in a focal individual. 

}
\usage{data("DT_ige")}
\format{
  The format is:
 chr "DT_ige"
}
\source{
This data was masked from a shared study.
}
\references{

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_ige)
DT <- DT_ige

\donttest{
########################################
############### sommer #################
########################################
if(requireNamespace("sommer")){
library(sommer)
# Indirect genetic effects model without covariance between DGE and IGE
modIGE <- mmes(trait ~ block, dateWarning = FALSE,
               random = ~ focal + neighbour,
               rcov = ~ units, nIters=100,
              data = DT)
summary(modIGE)$varcomp
pmonitor(modIGE)

# Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
modIGE <- mmes(trait ~ block, dateWarning = FALSE,
               random = ~ covm( vsm(ism(focal)), vsm(ism(neighbour)) ),
               rcov = ~ units, nIters=100,
              data = DT)
summary(modIGE)$varcomp
pmonitor(modIGE)

# form relationship matrix
Ai <- solve(A_ige + diag(1e-5, nrow(A_ige),nrow(A_ige) ))
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, 'inverse')=TRUE
# Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
modIGE <- mmes(trait ~ block, dateWarning = FALSE,
               random = ~ covm( vsm(ism(focal), Gu=Ai), vsm(ism(neighbour), Gu=Ai) ),
               rcov = ~ units, nIters=100,
              data = DT)
summary(modIGE)$varcomp
pmonitor(modIGE)

}
##############################################
############### lme4breeding #################
##############################################
if(requireNamespace("lme4breeding")){
library(lme4breeding)
# Indirect genetic effects model without covariance between DGE and IGE
modIGE <- lmeb(trait ~ block + (1|focal) + (1|neighbour),
              data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))

## Add relationship matrices
A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
modIGE <- lmeb(trait ~ block + (1|focal) + (1|neighbour),
                   relmat = list(focal=A_ige,
                                 neighbour=A_ige),
                   data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))

## Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
## Relationship matrix
A_ige <- A_ige + diag(1e-4, ncol(A_ige), ncol(A_ige) )
## Define 2 dummy variables to make a fake covariance
## for two different random effects
DT$fn <- DT$nn <- 1
## Create the incidence matrix for the first random effect
Zf <- Matrix::sparse.model.matrix( ~ focal-1, data=DT )
colnames(Zf) <- gsub("focal","", colnames(Zf))
## Create the incidence matrix for the second random effect
Zn <- Matrix::sparse.model.matrix( ~ neighbour-1, data=DT )
colnames(Zn) <- gsub("neighbour","", colnames(Zn))
## Fit the model
modIGE <- lmeb(trait ~ block + (0+fn+nn|both),
                   addmat = list(both=list(Zf,Zn)),
                   relmat = list(both=A_ige),
                   data = DT)
vc <- VarCorr(modIGE); print(vc,comp=c("Variance"))
sigma(modIGE)^2 # error variance

blups <- ranef(modIGE)
condVAR <- lapply(blups, function(x){attr(x, which="postVar")}) # take sqrt() for SEs
pairs(blups$both)
cov2cor(vc$both)

}

########## end ###########

}


}
\keyword{datasets}

