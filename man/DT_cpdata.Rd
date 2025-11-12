\name{DT_cpdata}
\alias{DT_cpdata}
\alias{GT_cpdata}
\alias{MP_cpdata}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population
}
\description{
A CP population or F1 cross is the designation for a cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection.

This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("DT_cpdata")}
\format{
  The format is:
 chr "DT_cpdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.


}
\examples{

data(DT_cpdata)
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
## create the variance-covariance matrix
A <- A.matr(GT) # additive relationship matrix
A <- A + diag(1e-4, ncol(A), ncol(A))
## look at the data and fit the model
head(DT)

\donttest{

######## sommer #########

if(requireNamespace("sommer")){
library(sommer)
mix1 <- mmes(Yield~1, henderson=FALSE,
              random=~vsm(ism(id),Gu=A)
                      + Rowf + Colf,
              rcov=~units,
              data=DT)
summary(mix1)$varcomp

## mmec uses the inverse of the relationship matrix
Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, 'inverse')=TRUE
mix2 <- mmes(Yield~1, henderson=TRUE,
              random=~vsm(ism(id),Gu=Ai)
                      + Rowf + Colf,
              rcov=~units,
              data=DT)
summary(mix2)$varcomp

vg <- summary(mix2)$varcomp[1,1] # genetic variance
G <- A*vg # genetic variance-covariance
Ci <- mix2$Ci # coefficient matrix
ind <- as.vector(mix2$partitions$`vsm(ism(id), Gu = Ai)`)
ind <- seq(ind[1],ind[2])
Ctt <- Ci[ind,ind] # portion of Ci for genotypes
R2 <- (G - Ctt)/G # reliability matrix
mean(diag(R2)) # average reliability of the trial
####====================####
#### multivariate model ####
####     2 traits       ####
####====================####
traits <- c("color","Yield")
DT[,traits] <- apply(DT[,traits],2,scale)
DTL <- reshape(DT[,c("id", traits)],
               idvar = c("id"),
               varying = traits,
               v.names = "value", direction = "long",
               timevar = "trait", times = traits )
DTL <- DTL[with(DTL, order(trait)), ]
head(DTL)

# if using mmes=TRUE you need to provide the inverse
Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, 'inverse')=TRUE
#### be patient take some time
ansm <- mmes( value ~ trait, # henderson=TRUE,
               random=~ vsm(usm(trait), ism(id), Gu=A), # Ai if henderson
               rcov=~ vsm(dsm(trait), ism(units)),
               data=DTL)
cov2cor(ansm$theta[[1]])

}

######## lme4breeding #########

if(requireNamespace("lme4breeding")){
library(lme4breeding)
mix1 <- lmeb(Yield~ (1|id) + (1|Rowf) + (1|Colf),
                 relmat=list(id=A),
                 data=DT)
vc <- VarCorr(mix1); print(vc,comp=c("Variance"))
sigma(mix1)^2 # error variance

# run one last iteration with imputed data
# to make sure you get predictions for every level
DT2 <- DT
DT2$Yield <- imputev(DT2$Yield)
mix2 <- update(mix1, suppressOpt = TRUE,
               start=getME(mix1, "theta"), 
               data=DT2)
predsMix2 <- ranef(mix2)
condVAR <- lapply(predsMix2, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

# if you don't want the imputed vector to have an effect in
# the predictions you can use the getMME function to use
# the extended model and get predictions without including the 
# imputed data (I know is a bit messy)
preds <- getMME(object=mix2, # extended model
                vc=VarCorr(mix1), # variance components
                recordsToUse = which(!is.na(DT$Yield)) # records to use for MME
                )
# now you could compare between both types of predictions, the last ones are in 
# theory the correct ones.
plot(preds$bu[2:364,], predsMix2$id[,1])

}

######## evola #########

if(requireNamespace("evola")){
library(evola)
DT$occ <- 1
DT$Yield <- imputev(DT$Yield)
A <- A[DT$id,DT$id]
# get best 20 individuals weighting variance by ~0.5=(30*pi)/180
res<-evolafit(formula=cbind(Yield, occ)~id, dt= DT, 
              # constraints: if sum is greater than this ignore 
              constraintsUB = c(Inf,20), 
              # constraints: if sum is smaller than this ignore
              constraintsLB= c(-Inf,-Inf), 
              # weight the traits for the selection
              b = c(1,0), 
              # population parameters
              nCrosses = 100, nProgeny = 10, 
              recombGens=1, nChr=1, mutRateAllele=0,
              # coancestry parameters
              D=A, lambda= (30*pi)/180 , nQtlStart = 20, 
              # selection parameters
              propSelBetween = 0.5, propSelWithin =0.5, 
              nGenerations = 40) 
              
Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
best = bestSol(res$pop)[,"fitness"];best
qa = (Q \%*\% DT$Yield)[best,]; qa 
qDq = Q[best,] \%*\% A \%*\% Q[best,]; qDq 
sum(Q[best,]) # total # of inds selected

evolmonitor(res)

plot(DT$Yield, col=as.factor(Q[best,]), 
     pch=(Q[best,]*19)+1)

pareto(res)

}

############ end ############

}

}
\keyword{datasets}

