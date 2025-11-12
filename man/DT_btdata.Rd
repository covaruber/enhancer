\name{DT_btdata}
\alias{DT_btdata}
\docType{data}
\title{
Blue Tit Data for a Quantitative Genetic Experiment
}
\description{

a data frame with 828 rows and 7 columns, with variables tarsus length (tarsus) and colour (back) measured on 828 individuals (animal). The mother of each is also recorded (dam) together with the foster nest (fosternest) in which the chicks were reared. The date on which the first egg in each nest hatched (hatchdate) is recorded together with the sex (sex) of the individuals.

}
\usage{data("DT_btdata")}
\format{
  The format is:
 chr "DT_btdata"
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

}
\examples{

data(DT_btdata)
DT <- DT_btdata
head(DT)

\donttest{

############# sommer ################

if(requireNamespace("sommer")){
library(sommer)
mix4 <- mmes(tarsus ~ sex,
             random = ~ dam + fosternest,
             rcov=~units,
             data = DT)
summary(mix4)$varcomp

# MULTI-TRAIT EXAMPLE

traits <- c("tarsus","back","hatchdate")
DT[,traits] <- apply(DT[,traits],2,scale)
DTL <- reshape(DT[,c("animal","sex","dam","fosternest", traits)],
               idvar = c("animal","sex","dam","fosternest"),
               varying = traits,
               v.names = "value", direction = "long",
               timevar = "trait", times = traits )
DTL <- DTL[with(DTL, order(trait,animal)), ]
head(DTL)

mix3 <- mmes(value ~ trait:sex - 1, henderson=TRUE,
                random = ~ vsm(usm(trait),ism(dam)) +
                           vsm(usm(trait), ism(fosternest)),
                rcov= ~ vsm(dsm(trait),ism(units)),
                data = DTL)


summary(mix3)$varcomp
#### calculate the genetic correlation
cov2cor(mix3$theta[[1]])
cov2cor(mix3$theta[[2]])

}

############# lme4breeding ################

if(requireNamespace("lme4breeding")){
library(lme4breeding)
mix4 <- lmeb(tarsus ~ sex + (1|dam) + (1|fosternest),
             data = DT)
vc <- VarCorr(mix4); print(vc,comp=c("Variance"))
sigma(mix4)^2 # error variance
BLUP <- ranef(mix4, condVar=TRUE)
PEV <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

### multi-trait model
traits <- c("tarsus" ,"back", "hatchdate")
for(iTrait in traits){DT[,iTrait] <- scale(DT[,iTrait])}
DTL <- reshape(DT[,c("animal", traits)], idvar = "animal", varying = traits,
               v.names = "value", direction = "long",
               timevar = "trait", times = traits )
DTL <- merge(DTL, unique(DT[,c("animal","dam","fosternest","sex")]), 
             by="animal", all.x = TRUE)
DTL <- DTL[with(DTL, order(trait)), ]
head(DTL)

system.time(
  mix <- lmeb(value ~ (0+trait|dam),
                  # relmat = list(geno=A),
                  # rotation = TRUE,
                  data=DTL)
)
vc <- VarCorr(mix); print(vc,comp=c("Variance"))
cov2cor(vc$dam)
sigma(mix)^2 # error variance

}

############# end ################

}

}
\keyword{datasets}
