\name{DT_technow}
\alias{DT_technow}
\alias{Ad_technow}
\alias{Af_technow}
\alias{Md_technow}
\alias{Mf_technow}
\docType{data}
\title{
Genotypic and Phenotypic data from single cross hybrids (Technow et al.,2014)
}
\description{
This dataset contains phenotpic data for 2 traits measured in 1254 single cross hybrids coming from the cross of Flint x Dent heterotic groups. In addition contains the genotipic data (35,478 markers) for each of the 123 Dent lines and 86 Flint lines. The purpose of this data is to demosntrate the prediction of unrealized crosses (9324 unrealized crosses, 1254 evaluated, total 10578 single crosses). We have added the additive relationship matrix (A) but can be easily obtained using the A.matr function on the marker data. Please if using this data for your own research cite Technow et al. (2014) publication (see References).

}
\usage{data("DT_technow")}
\format{
  The format is:
 chr "DT_technow"
}
\source{
This data was extracted from Technow et al. (2014).
}
\references{
If using this data for your own research please cite:

Technow et al. 2014. Genome properties and prospects of genomic predictions of hybrid performance in a Breeding program of maize. Genetics 197:1343-1355.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

Giovanny Covarrubias-Pazaran (2024).  evola: a simple evolutionary algorithm for complex problems. To be submitted to Bioinformatics.

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Chen GK, Marjoram P, Wall JD (2009). Fast and Flexible Simulation of DNA Sequence Data. Genome Research, 19, 136-142. http://genome.cshlp.org/content/19/1/136.

}
\examples{

data(DT_technow)
DT <- DT_technow

Md <- apply(Md_technow,2,as.numeric)
rownames(Md) <- rownames(Md_technow)
Mf <- apply(Mf_technow,2,as.numeric)
rownames(Mf) <- rownames(Mf_technow)

Md <- (Md*2) - 1
Mf <- (Mf*2) - 1

Ad <- A.matr(Md)
Af <- A.matr(Mf)
Ad <- Ad + diag(1e-4, ncol(Ad), ncol(Ad))
Af <- Af + diag(1e-4, ncol(Af), ncol(Af))

\donttest{

################# sommer ####################
if(requireNamespace("sommer")){
library(sommer)
ans2 <- mmes(GY~1,
             random=~vsm(ism(dent),Gu=Ad) + vsm(ism(flint),Gu=Af),
             rcov=~units,
             data=DT)
summary(ans2)$varcomp

Adi <- solve(Ad + diag(1e-4,ncol(Ad),ncol(Ad)))
Adi <- as(as(as( Adi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Adi, 'inverse')=TRUE
Afi <- solve(Af + diag(1e-4,ncol(Af),ncol(Af)))
Afi <- as(as(as( Afi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Afi, 'inverse')=TRUE
####====================####
#### multivariate model ####
####     2 traits       ####
####====================####
head(DT)

traits <- c("GY","GM")
DT[,traits] <- apply(DT[,traits],2,scale)
DTL <- reshape(DT[,c("hybrid","dent","flint", traits)],
               idvar = c("hybrid","dent","flint"),
               varying = traits,
               v.names = "value", direction = "long",
               timevar = "trait", times = traits )
DTL <- DTL[with(DTL, order(trait,hybrid)), ]
head(DTL)

M <- rbind(Md,Mf)
A <- A.mat(M)
Ai <- solve(A + diag(1e-4,ncol(A),ncol(A)))
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, 'inverse')=TRUE

ans3 <- mmes(value~trait, henderson=TRUE,
             random=~vsm(usm(trait),ism(overlay(dent,flint)),Gu=Ai),
             rcov=~ vsm(dsm(trait), ism(units)),
             data=DTL)
summary(ans3)
cov2cor(ans3$theta[[1]])

}

################# lme4breeding ####################

if(requireNamespace("lme4breeding")){
library(lme4breeding)
## simple model
ans2 <- lmeb(GY ~ (1|dent) + (1|flint),
                 data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
BLUP <- ranef(ans2, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

### with relationship matrices
ans2 <- lmeb(GY ~ (1|dent) + (1|flint),
             relmat = list(dent=Ad,
                           flint=Af),
             data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))

### overlayed model
M <- rbind(Md,Mf)
A <- A.matr(M)
A <- A + diag(1e-4,ncol(A), ncol(A))
Z <- with(DT, overlay(dent,flint) )
Z = Z[which(!is.na(DT$GY)),]
#### model using overlay without relationship matrix
ans2 <- lmeb(GY ~ (1|fema),
                 addmat = list(fema=Z),
                 relmat = list(fema=A),
                 data=DT)
vc <- VarCorr(ans2); print(vc,comp=c("Variance"))
sigma(ans2)^2 # error variance
BLUP <- ranef(ans2, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

### rotated model for hybrids
H <- kronecker(Ad,Af, make.dimnames = TRUE)
H <- H[which(colnames(H)\%in\% DT$hy),which(colnames(H)\%in\% DT$hy)]

ans3 <- lmeb(GY ~ (1|hy),
             relmat = list(hy=H), 
             rotation=TRUE,
             data=DT)
vc <- VarCorr(ans3); print(vc,comp=c("Variance"))
BLUP <- ranef(ans3, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

}
################# evola ####################

if(requireNamespace("evola")){
library(evola)
DT <- DT_technow
  DT$occ <- 1; DT$occ[1]=0
  combos <- build.HMM(Md,Mf, return.combos.only = TRUE)
  combos <- combos$data.used[which(combos$data.used$hybrid \%in\%DT$hy),]
  M <- build.HMM(Md,Mf, custom.hyb = combos)
  A <- A.matr(M$HMM.add)
  A <- A[DT$hy,DT$hy]
  # run the genetic algorithm
  # we assig a weight to x'Dx of (20*pi)/180=0.34
  res<-evolafit(formula = c(GY, occ)~hy,
                dt= DT, 
                # constraints: if sum is greater than this ignore
                constraintsUB = c(Inf,100), 
                # constraints: if sum is smaller than this ignore
                constraintsLB= c(-Inf,-Inf),
                # weight the traits for the selection
                b = c(1,0), 
                # population parameters
                nCrosses = 100, nProgeny = 10, 
                recombGens=1, nChr=1, mutRateAllele=0,
                # coancestry parameters
                D=A, lambda= (20*pi)/180 , nQtlStart = 90, 
                # selection parameters
                propSelBetween = 0.5, propSelWithin =0.5, 
                nGenerations = 20) 
  
  Q <- pullQtlGeno(res$pop, simParam = res$simParam, trait=1); Q <- Q/2
  best = bestSol(res$pop)[,"fitness"]
  qa = (Q \%*\% DT$GY)[best,]; qa 
  qAq = Q[best,] \%*\% A \%*\% Q[best,]; qAq 
  sum(Q[best,]) # total # of inds selected
  
  evolmonitor(res)
  plot(DT$GY, col=as.factor(Q[best,]), 
       pch=(Q[best,]*19)+1)
       
  pareto(res)
  
}
################## end ###################
}

}
\keyword{datasets}
