\name{DT_sleepstudy}
\alias{DT_sleepstudy}
\docType{data}
\title{
Reaction times in a sleep deprivation study
}
\description{
The average reaction time per day for subjects in a sleep deprivation study. On day 0 the subjects had their normal amount of sleep. Starting that night they were restricted to 3 hours of sleep per night. The observations represent the average reaction time on a series of tests given each day to each subject. Data from sleepstudy to see how lme4 models can be translated in sommer.

}
\usage{data("DT_sleepstudy")}
\format{
  The format is:
 chr "DT_sleepstudy"
}
\source{
These data are from the study described in Belenky et al. (2003), for the sleep deprived group and for the first 10 days of the study, up to the recovery period.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Giovanny Covarrubias-Pazaran (2024).  lme4breeding: enabling genetic evaluation in the age of genomic data. To be submitted to Bioinformatics.

Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.

Gregory Belenky et al. (2003) Patterns of performance degradation and restoration during sleep restrictions and subsequent recovery: a sleep dose-response study. Journal of Sleep Research 12, 1-12.


}
\examples{


data(DT_sleepstudy)
DT <- DT_sleepstudy
head(DT)

\donttest{

##############################################
###############  sommer.     #################
##############################################
if(requireNamespace("sommer")){
library(sommer)
##################################
fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject, 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject + vsm(ism(Days), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject + vsm(ism(Days), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

##################################
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(ism(Days), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
summary(fm2)$varcomp

}

##############################################
############### lme4breeding #################
##############################################
if(requireNamespace("lme4breeding")){
library(lme4breeding)
##################################
fm1 <- lmeb(Reaction ~ Days + (1 | Subject),
                data=DT)
vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
sigma(fm1)^2 # error variance

##################################
fm1 <- lmeb(Reaction ~ Days + (Days || Subject), data=DT)
vc <- VarCorr(fm1); print(vc,comp=c("Variance"))

##################################
fm1 <- lmeb(Reaction ~ Days + (Days | Subject), data=DT)
vc <- VarCorr(fm1); print(vc,comp=c("Variance"))

##################################
fm1 <- lmeb(Reaction ~ Days + (0 + Days | Subject), data=DT)
vc <- VarCorr(fm1); print(vc,comp=c("Variance"))

BLUP <- ranef(fm1, condVar=TRUE)
condVAR <- lapply(BLUP, function(x){attr(x, which="postVar")}) # take sqrt() for SEs

}
}

}
\keyword{datasets}
