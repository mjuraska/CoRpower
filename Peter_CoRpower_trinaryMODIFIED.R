# Main programs computepower() and computepower.n()
# for correlate of risk (CoR) power calculations and
# sample size calculations as described in Gilbert, Janes, and Huang (2015)
#
# September, 2014
# updated to new notation October, 2015
#
# The notation in this program matches that in the article
#####################################
# Trinary immune response S = 0,1,2
#####################################
# P2 = P(S=2) = P(S=hi)
# P1 = P(S=1) = P(S=med) (==0 for binary S)
# P0 = P(S=0) = P(S=lo)

# Suppose 3 latent baseline subgroups X = higher protection, medium protection, lower protection (2, 1, 0)

# Let Plat2, Plat1, Plat0 be the prevalence of the latent categories

# risk1(s,x) = infection risk for the subgroup S(1)=s,X=x if assigned vaccine
# risk0(s,x) = infection risk for the subgroup S(1)=s,X=x if assigned placebo

# VE(s,x) = 1 - risk1(s,x)/risk0(s,x)
# RR(s,x) = 1 - VE(s,x)

# Assume for simplicity risk0(s,x) = risk0 for all s,x
# Also assume P(Y=1|latent stratum x,s,Z=1) = P(Y=1|latent stratum x,Z=1)
# i.e., risk1(s,x) = risk1(x)

# Let the CoR effect size in vaccinees RR_t = risk1(2)/risk1(0),
# This is the x-axis of the power curves

# The program considers CoR analysis in the vaccine group only, where the goal is to test
# H0: risk1(2)=risk1(1)=risk1(0) vs. H1: risk1(2) < risk1(0), i.e.
# H0: RR_t = 1 vs. H1: RR_t < 1
# based on the case-control sample from vaccine recipients

# Sensitivity, Specificity, false negatives and false positives
# Sens \equiv P(S(1)=2|X=2) (==1 for S with no measurement error)
# Spec \equiv P(S(1)=0|X=0) (==1 for S with no measurement error)
# FP^0 \equiv P(S(1)=2|X=0) (==0 for S with no measurement error)
# FN^2 \equiv P(S(1)=0|X=2)  (==0 for S with no measurement error)
# FP^1 \equiv P(S(1)=2|X=1) (==0 for S with no measurement error)
# FN^1 \equiv P(S(1)=0|X=1) (==0 for S with no measurement error)

########################################################################
# End of unique notations for trinary immune response S = 0,1,2
# The remaining notations below are the same for trinary (Approach 2) and continous S
# except that PlatVElowest and VElowestvect only apply to continous S
#######################################################################

# The biomarker is measured at tau months and subjects are followed through taumax months

# RRoverall = 1 - overall vaccine efficacy
# annincinfectionplac = annual HIV infection incidence in placebo group
# annincdropout = annual dropout rate assumed the same in both groups

# numAtRiskTauCases     All cases in the vaccine group at-risk at tau and a case by taumax
#                             (regardless of whether the biomarker is measured)
# numAtRiskTauControls  All controls in the vaccine group at-risk at tau and not diseased at the end of follow-up taumax
#                             (regardless of whether the biomarker is measured)
# numAtRiskTauCasesPhase2  As above and also have the biomarker measured (i.e., in Phase 2)
#
# For computepower.n() vectors are inputted to repeat power calculations at multiple sample sizes
#     numAtRiskTauCasesVectALL
#     numAtRiskTauControlsVectALL
#     numAtRiskTauCasesVectPhase2

# sigma2obs  observed variance of the continuous marker S*
# rhos       vector of rho, the proportion of between vaccine recipient variability of S* that is
#            potentially protection relevant
# risk1      estimated probability that a vaccine recipient at-risk at tau experiences
#            the clinical endpoint by taumax.  It may be estimated differently for different studies.
# risk0      Same as risk1 for the placebo group
# This risk0 is used for both computepower and computepower.n

# PlatVElowest The percentage of vaccine recipients with the lowest value of VE (VElowestvect).

# VElowestvect A vector of the lowest possible value of vaccine efficacy.
#              Typical applications will range VElowestvect from 0 to VEoverall

# M                number of iterations of simulated clinical trials
# alpha            2-sided alpha level of CoR tests
# controlCaseRatio ratio of number of controls sampled per case in the vaccine arm
#######################################################################

###################################################################

# Inputs:
# Fixed sample size
# VE, risk0, P^{lat}_{0}, P^{lat}_{2} (single values)
# Ranges of values for VE^{lat}_{0} and VE^{lat}_{1}
# One value each of P_{0} and P_{2} values, typically set to
# P_{2} = P^{lat}_{0} and P_{2} = P^{lat}_{2}
# Other analyses could specify alternative values of P_{0} and P_{2} to study how power varies with
# the selected thresholds for a high and low response.

# Determined from the above fixed values: P^{lat}_{1}, VE^{lat}_{2}

# This program considers the scenario RRlat1 = RRoverall, which yields the formula
# VEoverall = (Plat0*VElat0 + Plat2*VElat2)/(Plat0 + Plat2)
# For fixed values of Plat0 and Plat2, this links VElat0 to VElat2 by the
# formula VElat2 = (VEoverall*(Plat0+Plat2) - Plat0*VElat0)/Plat2, which
# is used in the computepower() function
#
# Following Step 7 in the article, the program for a trichotomous marker accounts for assay noise in one of two ways.
# The first approach specifies Spec, Sens, FP0, FN2 which determine FP1 and FN1 from equations (8) and (9). Spec, Sens, 
# FP0, and FN2 must all be vectors of the same length.
#
# The second approach specifies sigma2obs and rho; this approach assumes the normal measurement error model (4) in the article.
# Specifying Spec=NULL, FP0=NULL, Sens=NULL, FN2=NULL defaults to approach 2, which is used in the manuscript.

# Output: Power

#' Sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement
#' error and treatment efficacy
#'
#' Performs sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement error and treatment efficacy [Gilbert, Janes, and Huang (2015).
#' ``Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials.'']

#'
#' @param numAtRiskTauCases  Number of subjects in the vaccine group at-risk at tau and with the clinical event (cases) by taumax (regardless of whether the biomarker is measured).
#' @param numAtRiskTauCasesPhase2 Number of subjects in the vaccine group at-risk at tau and with the clinical event (cases) by taumax and with the biomarker measured (i.e., in Phase 2).
#' @param numAtRiskTauControls Number of subjects in the vaccine group at-risk at tau and without the clinical event (controls) by taumax (regardless of whether the biomarker is measured).
#' @param risk0 Estimated probability that a placebo recipient at-risk at tau experiences the clinical event by taumax.
#' @param RRoverall One minus the overall vaccine efficacy (VE).
#' @param Plat0 For a trichotomous biomarker (or binary as a special case), probability that the (latent) true biomarker takes the lowest (low) value.
#' @param Plat2 For a trichotomous biomarker (or binary as a special case), probability that the (latent) true biomarker takes the highest (high) value.
#' @param P0 For a trichotomous biomarker (or binary as a special case), probability that the measured/observed biomarker takes the lowest (low) value.  If unspecified, this parameter is set to \code{Plat0}.
#' @param P2 For a trichotomous biomarker (or binary as a special case), probability that the measured/observed biomarker takes the highest (high) value.  If unspecified, this parameter is set to \code{Plat2}.
#' @param RRlat0 For a trichotomous biomarker (or binary as a special case), a vector of values corresponding to one minus VE for the lowest (low) true biomarker subgroup.  Each value of \code{RRlat0} corresponds to one unique effect size (RR_t).
#' @param RRlat1 For a trichotomous biomarker, a vector of values corresponding to one minus VE for the middle true biomarker subgroup.  For a binary biomarker, \code{RRlat1} can be left unspecified.
#' @param PlatVElowest For a continuous biomarker, the percentage of vaccine recipients with the lowest value of VE.
#' @param VElowestvect For a continuous bioarker, a vector of values corresponding to the lowest possible value of VE.  Typical applications will range \code{VElowestvect} from 0 to 1 - \code{RRoverall}.
#' @param controlCaseRatio Number of controls sampled per case in the vaccine arm (i.e. sampled in to Phase 2).
#' @param M Number of simulated clinical trials.
#' @param alpha Two-sided type-I error rate for CoR hypothesis tests.
#' @param sigma2obs For a continuous biomarker, or for a trichotomous or binary biomarker simulated using `approach 2', the variance of the continuous measured/observed biomarker.
#' @param rhos For a continuous biomarker, or for a trichotomous or binary biomarker simulated using 'approach 2', a vector of length 4 with values for the fraction of protection-relevant variability in the measured/observed continuous biomarker.  The first element of this vector should be 1, corresponding to the case of no measurement error.
#' @param Spec For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the sensitivity of the measured biomarker. Specifying \code{Spec=NULL}, \code{FP0=NULL}, \code{Sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param FP0 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false positive rate (FP^1) of the measured/observed biomarker. Specifying \code{Spec=NULL}, \code{FP0=NULL}, \code{Sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param Sens For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the specificity of the measured/observed biomarker. Specifying \code{Spec=NULL}, \code{FP0=NULL}, \code{Sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param FN2 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false negative rate (FN^2) of the measured/observed biomarker. Specifying \code{Spec=NULL}, \code{FP0=NULL}, \code{Sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param tpsMethod Character denoting method for fitting the logistic regression model. Choose from "PL" for pseudo-likelihood (default), "ML" for maximum likelihood, and "WL" for weighted likelihood. 
#' @param biomType Type of biomarker that is used. The default is "continuous"; other choices are "trichotomous" and "binary".
#' @param cohort Sampling design to be used. Default is \code{FALSE}, specifying case-control sampling design. If \code{TRUE}, case-cohort sampling is used. 
#' @param p For case-cohort sampling design, probability that a subject will be in the cohort. 
#' 
#' @details
#' This function performs simulations to calculate the power for testing whether a trichotomous or continuous biomarker
#' is a correlate of risk (CoR).
#'
#' This program considers the scenario \code{RRlat1 = RRoverall}, which yields the formula
#' \code{1-RRoverall = (Plat0*VElat0 + Plat2*VElat2)/(Plat0 + Plat2)}.
#' For fixed values of \code{Plat0} and \code{Plat2}, this links \code{1-RRlat0} to \code{1-RRlat2} by the
#' formula \code{1-RRlat2 = (1-RRoverall*(Plat0+Plat2) - Plat0*(1-RRlat0))/Plat2}, which
#' is what the function uses to obtain \code{RRlat2}.
#'
#' Following Step 7 in the manuscript, the measurement error in a trichotomous (or binary) biomarker is accounted for in one of two ways.
#' Approach 1 specifies \code{Spec}, \code{Sens}, \code{FP0}, and \code{FN2} which determine
#' \code{FP1} and \code{FN1} from equations (7) and (8).  Four values are required for each input parameter, to allow the evaluation of biomarkers with different levels of measurement error.
#'
#' Approach 2 for a trichotomous (or binary) biomarker specifies \code{sigma2obs} and \code{rho}; this approach assumes the normal measurement error model (4) in the manuscript.  
#' Specifying \code{Spec=NULL}, \code{FP0=NULL}, \code{Sens=NULL}, and \code{FN2=NULL}
#' defaults to approach 2, which is what is used in illustrations in the manuscript.
#'
#' For a continuous biomarker, \code{VElowestvect}, \code{sigma2obs} and \code{rho} must be specified.  Setting \code{VElowestvect = NULL}
#' indicates that the biomarker is not continuous (it is trichotomous).
#'
#' This program implements a scenario with without-replacement-sampling (e.g., typically used in case-control and 2-phase sampling)
#'
#' @return Power- the fraction of simulated trials in which the null hypothesis H_0 (expression (14) of the manuscript for a trichotomous (or binary) biomarker and expression (16) for a continuous biomarker) is rejected.
#'
#' @examples
#' ## 'Global' parameters (independent of marker)
#' VEoverall <- 0.26   # VE in the at-risk-month-tau cohort
#' RRoverall <- 1 - VEoverall
#' numAtRiskTauCases <- 41
#' numAtRiskTauControls <- 7662
#' numAtRiskTauCasesPhase2 <- 41
#' risk1 <- numAtRiskTauCases/(numAtRiskTauCases + numAtRiskTauControls)
#' risk0 <- risk1/RRoverall # risk in placebo
#'
#' ## Parameters used for the trichotomous or binary biomarker calculations, Approach 1
#' ## For no measurement error scenario, set
#' ## Spec[1]=1, FP0[1]=0, Sens[1]=1, FN2[1]=0
#' ## Note the other elements of at least one of these four parameter vectors need
#' ## to be set to a different value to tell the program to use Approach 1
#' Spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' Sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' ## Parameters used for continuous biomarker calculations, or
#' ## for trichotomous/binary biomarkers under Approach 2
#' ## Note these values are needed but are immaterial trichotomous/binary Approach 1
#' sigma2obs <- 1
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowestvect <- seq (0, 1-RRoverall,len=100)
#'
#' #################################################
#' ## use the function to perform power calculations
#' #################################################
#'
#' ## Binary biomarker, Approach 1 ##
#' Spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' Sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' ## i.e. Plat1=P1=0 and Plat0 = 1- Plat2
#' Plat0 <- 0.1
#' Plat2 <- 1-Plat0
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rhos, Spec, FP0, Sens, FN2)
#'
#'
#'## Trichotomous biomarker, Approach 1 ##
#' Spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' Sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,
#' VElowestvect=NULL, controlCaseRatio, M, alpha=0.05, sigma2obs, rhos, Spec, FP0, Sens, FN2)
#'
#'
#' ## Binary biomarker, Approach 2 ##
#' Spec <- Sens <- FP0 <- FN2 <- NULL
#'
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' Plat0 <- 0.1
#' Plat2 <- 1-Plat0
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls, risk0,
#' RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL, controlCaseRatio,
#' M, alpha=0.05, sigma2obs, rhos, Spec, FP0, Sens, FN2)
#'
#'
#' ## Trichotomous biomarker, Approach 2 ##
#' Spec <- Sens <- FP0 <- FN2 <- NULL
#'
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve

#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rhos, Spec, FP0, Sens, FN2)
#'
#'
#' ## Continuous biomarker ##
#' sigma2obs <- 1
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowestvect <- seq (0, 1-RRoverall,len=100)
#'
#' ##Note these values are needed but are immaterial for continous markers
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest,VElowestvect,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rhos)
#'
#' #######################
#' ## plotting the results
#' #######################
#'
#' ## trichotomous biomarker Approach 1
#' Spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' Sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#' sigma2obs <- 1
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case


#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,
#' VElowestvect=NULL, controlCaseRatio, M, alpha=0.05, sigma2obs=NULL, rhos=NULL, Spec, FP0, Sens, FN2)

#' ## plot power vs. CoR risk ratio in vaccine group (hi vs. lo)  (Figure 4)
#' ## file name = paste("powerstrinary",P2,P1,controlCaseRatio,".dat")
#' powstrin40105 <- matrix(scan('powerstrinary0.40.15.dat'),ncol=4,byrow=T)

#' VEoverall <- scan('VEoverallCoRpower.dat')
#' alpha <- scan('alpha.dat')
#' rhos <- scan('rhos.dat')
#' RRlat2 <- scan('RRlat2.dat')
#' RRlat0 <- scan('RRlat0.dat')
#' N <- scan('sampsizeALL.dat')
#' nPhase2 <- scan('numbeventsPhase2.dat')

#' corrr <- RRlat2/RRlat0
#'
#' postscript("powertrichotomous_CoR.eps",horizontal=T)
#' par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,mar=c(5.1,4.6,3.1,1.1),oma=c(5,3,5,5),las=1)
#' xlims <- c(min(corrr),1)
#' plot(corrr,powstrin40105[,1],xlim=xlims,ylim=c(0,1),type='n',xlab="CoR Risk Ratio RR_t(rho=1) in Vaccine Group (Hi vs. Lo)",ylab="Power")
#' title("P0 = Plat0 = 0.1; P2 = Plat2 = 0.4")
#'
#' axis(1)
#' axis(2)
#' lines(corrr,powstrin40105[,1],lty=1,col="blue",lwd=4)
#' lines(corrr,powstrin40105[,2],lty=2,col="orange",lwd=4)
#' lines(corrr,powstrin40105[,3],lty=3,col="green",lwd=4)
#' lines(corrr,powstrin40105[,4],lty=4,col="black",lwd=4)
#'
#' legend(x="topright",legend=c(paste("Power rho=",rhos[1],sep=""),paste("Power rho=",rhos[2],sep=""),paste("Power rho=",rhos[3],sep=""),paste("Power rho=",rhos[4],sep="")),
#' lty=c(1,2,3,4),col=c("blue","orange","green","black"),lwd=2)
#' mtext(paste("Power to Detect a Trichotomous CoR in Vaccine Recipients [2-sided alpha = ",alpha,"]"),outer=T,cex=1.3)
#' mtext(paste("Overall VE = ",VEoverall,"; Number controls  = ",round(controlCaseRatio*nPhase2),"; Number cases = ",round(nPhase2),"; Controls:cases = ",
#' controlCaseRatio,":1"),side=1,line=0.7,outer=T,cex=1.3)
#' mtext(paste("VElat_0 varies from ",VEoverall," to 0 as VElat_2 varies from ",VEoverall," to ",round(2*VEoverall,2)),side=1,line=3,outer=T,cex=1.3)
#' ## Note: the upper limit is VEoverall*(PlatloVE+PlathiVE)/PlathiVE, in the special case of this plot with
#' ## PlatloVE = PlathiVE, this simplifies to 2*VE_0.
#' dev.off()
#'
#' ## continuous biomarker
#' sigma2obs <- 1
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowestvect <- seq (0, 1-RRoverall,len=100)
#' ##Note these values are needed but are immaterial for continous markers
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest,VElowestvect,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rhos)


#' ## Plot Power vs. CoR Relative Risk per SD Increase in X*  (Figure 2)
#' ## Power to Detect a Normally Distributed CoR in Vaccine Recipients [ 2-sided alpha  0.05]
#' VEoverall <- scan('VEoverallCoRpower.dat')
#' alpha <- scan('alpha.dat')
#' rhos <- scan('rhos.dat')
#' RRlat2 <- scan('RRlat2.dat')
#' RRlat0 <- scan('RRlat0.dat')
#' N <- scan('sampsizeALL.dat')
#' nPhase2 <- scan('numbeventsPhase2.dat')
#' PlatVElowest <- scan('PlatVElowest.dat')
#' RRs <- scan('RRs.dat')
#' reverseRRs <- RRs[length(RRs):1]
#' ## The 5 refers to a controlCaseRatio of 5:1
#' powscont5 <- matrix(scan('powerscont5.dat'),ncol=4,byrow=T)
#' reversepowscont <- powscont5[nrow(powscont5):1,]
#'
#' sigma2tr <- 1
#' sigma2erho1 <- ((1-rhos[1])/rhos[1])*sigma2tr
#' sigma2erho2 <- ((1-rhos[2])/rhos[2])*sigma2tr
#' sigma2erho3 <- ((1-rhos[3])/rhos[3])*sigma2tr
#' sigma2erho4 <- ((1-rhos[4])/rhos[4])*sigma2tr
#' benchmarklab="V2"
#' benchmarkestRRrho1 <- 0.57 # From Haynes et al. 2012
#' benchmarkestRRrho2 <- benchmarkestRRrho1^(1/sqrt(rhos[2]))
#' benchmarkestRRrho3 <- benchmarkestRRrho1^(1/sqrt(rhos[3]))
#' benchmarkestRRrho4 <- benchmarkestRRrho1^(1/sqrt(rhos[4]))
#'
#' postscript("powercontinuous_CoR.eps",horizontal=T)
#' par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,oma=c(3,3,4,4))
#' plot(RRs,powscont5[,1],ylim=c(0,1),type='n',xlab="Vaccine Group CoR Relative Risk RR_c per 1 SD Increase in the True Biomarker X*",ylab="Power",axes=FALSE)
#' box()
#'
#' m <-length(RRs)
#' RRsgrid <- seq(RRs[1],RRs[m],len=10)
#' axis(1,at=RRsgrid,labels=round(RRsgrid,2))
#' axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0),labels=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0))
#'
#' lines(RRs,powscont5[,1],lty=1,col="blue",lwd=4)
#' lines(RRs,powscont5[,2],lty=2,col="orange",lwd=4)
#' lines(RRs,powscont5[,3],lty=3,col="green",lwd=4)
#' lines(RRs,powscont5[,4],lty=4,col="black",lwd=4)
#'
#' ind <- c(1:m)[(RRs - benchmarkestRRrho1) <= 0 & (RRs - benchmarkestRRrho1)==max((RRs - benchmarkestRRrho1)[(RRs - benchmarkestRRrho1)<=0])]
#' if (length(ind) > 0) {
#'  text(RRs[ind],powscont5[ind,1],benchmarklab,cex=1.3,col="blue") }
#'
#' ind <- c(1:m)[(RRs - benchmarkestRRrho2) <= 0 & (RRs - benchmarkestRRrho2)==max((RRs - benchmarkestRRrho2)[(RRs - benchmarkestRRrho2)<=0])]
#' if (length(ind) > 0) {
#'  text(RRs[ind],powscont5[ind,2],benchmarklab,cex=1.3,col="orange") }
#'
#' ind <- c(1:m)[(RRs - benchmarkestRRrho3) <= 0 & (RRs - benchmarkestRRrho3)==max((RRs - benchmarkestRRrho3)[(RRs - benchmarkestRRrho3)<=0])]
#' if (length(ind) > 0) {
#'  text(RRs[ind],powscont5[ind,3],benchmarklab,cex=1.3,col="green") }
#'
#' ind <- c(1:m)[(RRs - benchmarkestRRrho4) <= 0 & (RRs - benchmarkestRRrho4)==max((RRs - benchmarkestRRrho4)[(RRs - benchmarkestRRrho4)<=0])]
#' if (length(ind) > 0) {
#'  text(RRs[ind],powscont5[ind,4],benchmarklab,cex=1.3,col="black") }
#'
#' abline(h=alpha/2,lty=3)
#'
#' text(RRs[round(m/3)],.95,paste("Fraction of subjects with VE=0: ",PlatVElowest),cex=1.45)
#'
#' legend(x="topright",legend=c(paste("Power rho=",rhos[1],sep=""),paste("Power rho=",rhos[2],sep=""),paste("Power rho=",rhos[3],sep=""),paste("Power rho=",rhos[4],sep="")),
#'       lty=c(1,2,3,4),col=c("blue","orange","green","black"),lwd=2,cex=1.32)
#'
#' title(paste("Power to Detect a Normally Distributed CoR in Vaccine Recipients [2-sided alpha = ",alpha,"]"))
#' mtext(paste("Overall VE = ",VEoverall,"; Number controls  = ",round(nPhase2*controlCaseRatio),"; Number cases = ",round(nPhase2),"; Controls:cases = ",
#' controlCaseRatio,":1"),side=1,line=2,outer=T,cex=1.3)
#' dev.off()
#'

#' @import survival
#' @import osDesign
#' @export
computepower <- function(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
                         risk0, RRoverall,
                         Plat0, Plat2,
                         P0=Plat0, P2=Plat2,
                         RRlat0=seq(1,RRoverall,len=20), RRlat1=rep(RRoverall,20),
                         PlatVElowest, VElowestvect,
                         controlCaseRatio=5,
                         M=100,
                         alpha=0.05,
                         sigma2obs=1,
                         rhos=c(1,0.9,0.7,0.5),
                         Spec=NULL, FP0=NULL, Sens=NULL, FN2=NULL,
                         tpsMethod=c("PL", "ML","WL"),
                         biomType=c("continuous", "trichotomous", "binary"),
                         cohort=FALSE, p=NULL) {
  
  # sigma2tr is the variance of the true biomarker X
  # rhos must be a vector with 4 values and is for the continuous and binary biomarker correlates correlations
  # RRlat0 is the span of true relative risks (vaccine vs. placebo) in the latent lower protected subgroup
  # RRlat1 is the span of true relative risks (vaccine vs. placebo) in the latent medium protected subgroup
  # The power calculations should always include rho=1 in the first element of the rhos vector,
  # as the best case scenario, and the plotting functions assume this.
  
  # VElowestvect is used for a continuous biomarker- a vector of fixed value of VE(x) for
  # the subgroup of subjects with lowest X^* values, where this subgroup has prevalence PlatVElowest

  
  tpsMethod <- match.arg(tpsMethod, choices = c("PL","ML","WL"))
  biomType <- match.arg(biomType, choices = c("continuous", "trichotomous", "binary"))
  
  # check sampling design input parameters are specified and valid
  if(cohort==TRUE) {  #case-cohort 
    
    if (is.null(p)==TRUE) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is unspecified")
    } else if (p < 0 | p > 1) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is not a valid probability")
    }
    
  } else if (is.null(controlCaseRatio)==TRUE) {  #case-control
    stop("Case-control sampling was chosen and controlCaseRatio is unspecified")
  }

  # check biomarker type and input parameters match
  if(biomType=="binary" & (P0+P2 != 1)){
    stop("Binary biomarker was specified but P0 and P2 do not add up to 1")
  }
  if(biomType=="continuous" & is.null(VElowestvect)){
    stop("Continuous biomarker was specified but VElowestvect is NULL")
  }
  
  
  nCases <- numAtRiskTauCases
  nPhase2 <- numAtRiskTauCasesPhase2
  # nPhase2 = n_{cases} in manuscript
  # Overall denominator: number observed to be at risk when the immune response is measured (N = N in manuscript):
  N <- nCases + numAtRiskTauControls
  
  # Compute VElat2:
  VEoverall <- 1 - RRoverall
  VElat0 <- 1 - RRlat0
  VElat1 <- 1 - RRlat1
  VElat2 <- (VEoverall*(Plat0+Plat2) - Plat0*VElat0)/Plat2  #adjust function for this??
  # This formula for VElat0 assumes VElat1 = VEoverall
  RRlat2 <- 1-VElat2
  Plat1 <- 1 - Plat2 - Plat0
  P1 <- 1 - P0 - P2
  
  ### Stop function to 
  # Checks out that all values of RRlat2 are between 0 and 1 and that PlatVElowest meets bounds.
  # If there are incompatible values of RRlat2, consider making Plat0 smaller,
  # and/or making VElat0 larger, and then check again if all values of 
  # RRlat2 are properly between 0 and 1.
  checkprobabilityviolation <- function(VEoverall,VElat0,Plat0,Plat2,PlatVElowest,VElowestvect) {
    VElat2 <- (VEoverall*(Plat0+Plat2) - Plat0*VElat0)/Plat2
    RRlat2 <- 1 - VElat2
    cat(paste("RRlat2="),"\n")
    cat(paste(round(RRlat2,3)))
    cat("\n")
    if (length(RRlat2[RRlat2 < 0 ])>0) {stop("Input parameters violate probability constraints for trichotomous marker calculations")}
    if (min(VElowestvect)==0 & PlatVElowest > 1 - VEoverall) {stop("Input parameters PlatVElowest and VElowestvect violate probability constraints for normal marker calculations") }
  }
  
  # checkprobabilityviolation(VEoverall,VElat0,Plat0,Plat2,PlatVElowest,VElowestvect)
  
  sigma2e <- (1-rhos)*sigma2obs
  
  sigma2tr <- rhos*sigma2obs
  
  #################################################
  # Computations for a trinary biomarker
  if(biomType=="trichotomous" | biomType=="binary") {
  
    Approach2 <- (all(is.null(Spec), is.null(Sens), is.null(FP0), is.null(FN2)))
    
    if (Approach2) {  # Default choice
      
      # Compute Sens, Spec, FP0, FP1, FN2, FN1
      
      ans <- computeSensSpecFPFN(sigma2obs,rhos,Plat0,Plat2,P0,P2)
      Sens <- unlist(lapply(ans, function(x) x[[1,10]])) 
      Spec <- unlist(lapply(ans, function(x) x[[1,11]]))
      FP0 <- unlist(lapply(ans, function(x) x[[1,12]])) 
      FP1 <- unlist(lapply(ans, function(x) x[[1,13]])) 
      FN2 <- unlist(lapply(ans, function(x) x[[1,14]]))
      FN1 <- unlist(lapply(ans, function(x) x[[1,15]])) 
      
      
      # Write out a vector that can be used to make a table mapping (rho,sigma2obs) to the Sens etc. parameters
      
      x1 <- unlist(c(sigma2obs[1],rhos[1],Plat0,P0,Plat2,P2,Sens[1],Spec[1],FP0[1],FN2[1],FP1[1],FN1[1]))
      #cat("\n")
      #cat(paste(round(x1,3)),"\n")
      
      x2 <- unlist(c(sigma2obs[2],rhos[2],Plat0,P0,Plat2,P2,Sens[2],Spec[2],FP0[2],FN2[2],FP1[2],FN1[2]))
      #cat("\n")
      #cat(paste(round(x2,3)),"\n")
      
      x3 <- unlist(c(sigma2obs[3],rhos[3],Plat0,P0,Plat2,P2,Sens[3],Spec[3],FP0[3],FN2[3],FP1[3],FN1[3]))
      #cat("\n")
      #cat(paste(round(x3,3)),"\n")
      
      x4 <- unlist(c(sigma2obs[4],rhos[4],Plat0,P0,Plat2,P2,Sens[4],Spec[4],FP0[4],FN2[4],FP1[4],FN1[4]))
      #cat("\n")
      #cat(paste(round(x4,3)),"\n")
  
    }
    
    # Approach 1 in the manuscript:
    if (!Approach2) { #*# use given Sens, Spec, FP0, and FN2 params
      
      checkParamLengthsMatch <- function(Sens, Spec, FP0, FN2){
        lengths <- sapply(list(Sens,Spec,FP0,FN2), length)
        if(max(lengths) != min(lengths)){
          stop("Vector lengths differ for Sens, Spec, FP0, FN2")
        }
      }
      checkParamLengthsMatch(Sens,Spec,FP0,FN2)
      
      # Apply formula (8) in the manuscript
      FN1 <- (P0 - Spec*Plat0 - FN2*Plat2)/Plat1   #*#P0, Plat0, Plat2 given params
      
      # Apply formula (9) in the manuscript
      FP1 <- (P2 - Sens*Plat2 - FP0*Plat0)/Plat1
      
      # Check if an error in the ranges of values due to an out of
      # bounds input parameter
      
      if (any(FN1 < 0 | FN1 > 1 | FP1 < 0 | FP1 > 1)){
        stop("Approach 1 was used and one of the parameters Sens, Spec, FP0, FN2 is out of range")
      }
          # if (FN11 < 0 | FN11 > 1 | FN12 < 0 | FN12 > 1 | FN13 < 0 | FN13 > 1 | FN14 < 0 | FN14 > 1 | FP11 < 0 | FP11 > 1 | FP12 < 0 | FP12 > 1 | FP13 < 0 | FP13 > 1 | FP14 < 0 | FP14 > 1) {
          #   cat(paste("Approach 1 was used and one of the parameters Sens, Spec, FP0, FN2 is out of range"),"\n")
          # }
      
    }
    
    # Binary biomarker special case (to remove small values of P1x)
    if (biomType=="binary") {
      P1 <- 0
      P2 <- 1 - P0
    }
    
    # Compute the marginal risks:
    # Made it to the end of follow-up HIV negative
    risk1 <- RRoverall*risk0
    
    # Observed risks P(Y(1)=1|S(1)=0, 1, or 2)  #*# for diff values of rho; using Bayes' rule
    
    probX0_cond_S2 <- FP0*Plat0/P2
    probX1_cond_S2 <- FP1*Plat1/P2
    probX2_cond_S2 <- Sens*Plat2/P2
    risk1_2 <- (probX0_cond_S2 %o% RRlat0 + probX1_cond_S2 %o% RRlat1 + probX2_cond_S2 %o% RRlat2 )*risk0 # use outer product to get matrix with nrow=length(rho), ncol=length(RRlat0)
    probX0_cond_S0 <- Spec*Plat0/P0
    probX1_cond_S0 <- FN1*Plat1/P0
    probX2_cond_S0 <- FN2*Plat2/P0
    risk1_0 <- (probX0_cond_S0 %o% RRlat0 + probX1_cond_S0 %o% RRlat1 + probX2_cond_S0 %o% RRlat2)*risk0
    risk1_1 <- (risk1 - risk1_0*P0 - risk1_2*P2)/P1
    
    # Note: For the binary biomarker special case, the risk1medx are NA
    
    esvect <- risk1_2/risk1_0 # matrix with nrow=length(rho) and ncol=length(RRlat0)
    
    # Vaccine risks within the latent subgroups (independent of rho of course)
    
    risk1lat_2 <- RRlat2*risk0
    risk1lat_1 <- RRlat1*risk0
    risk1lat_0 <- RRlat0*risk0
    
  } else if (biomType=="continuous") {  
  
    #################################################
    # Computations for a continuous biomarker
    # Define the truebetas indexed by the user-specified vector VElowestvect
  
    o <- length(VElowestvect)
    nus <- sqrt(rhos*sigma2obs)*qnorm(PlatVElowest)
    truebetas <- rep(NA,o)
    alphalatvect <- rep(NA,o)
    
    for (l in 1:o) {
      
      risk1latnu <- (1-VElowestvect[l])*risk0
      
      f <- function(alpha) {
        
        g <- function(x) {
          rho <- 1
          piece1 <- exp(alpha*(1 - x/nus[1]))*(risk1latnu^(x/nus[1]))
          piece2 <- (1-risk1latnu)^(x/nus[1]) + piece1
          piece3 <- dnorm(x/(sqrt(rho*sigma2obs)))
          kernel <- (piece1/piece2)*piece3
          return(kernel) 
        }
        
        # nus[1] corresponds to rho=1
        logitterm <- integrate(g,lower=nus[1],upper=6)$value  #*# logit term in eqn (13)
        
        ans <- 1-VEoverall - (PlatVElowest*risk1latnu + logitterm)/risk0
        return(ans) 
      }
      
      alphalatvect[l] <- uniroot(f,lower=-10,upper=10)$root  #*# solve for alphalat using eqn (12)
      
      # Second solve for betalat:
      D <- risk1latnu
      truebetas[l] <- (log(D/(1-D)) - alphalatvect[l])/nus[1]
    }
  }  

  # Function for computing the infection probabilities of vaccinees
  # risk1cont = risk_1^{lat}(x*) in the manuscript
  risk1cont <- function(x,alphalat,betalat) {
    linpart <- alphalat + betalat*x
    ans <- exp(linpart)/(1+exp(linpart))
    return(ans) 
  }
  
  ###################################################
  # Power calculations repeated for M simulations
  
  powerstrinary <- matrix(0, nrow=nrow(esvect), ncol=ncol(esvect))
  rownames(powerstrinary) <- paste0(rep("rho"), seq(1,nrow(powerstrinary)))
  powerscont <- matrix(0, nrow=length(rhos), ncol=o)
  rownames(powerscont) <- paste0(rep("rho"), seq(1,nrow(powerscont)))
  
  
  for (i in 1:M) {
    
    if(biomType=="trichotomous" | biomType=="binary") {
    # Trinary biomarker:
      for (j in 1:ncol(esvect)) { # for each value of RRlat0
        
        # Use Bayes rule to compute
        # P(X=x|Y(1)=1) = rrlatx for x=lo, med, hi
        # These formulas require the simplifying assumption
        # risk_0(x,s_1) = risk_0(x) = risk_0 for all x, s_1
        rrlat0 <- risk1lat_0[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        rrlat1 <- risk1lat_1[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        rrlat2 <- risk1lat_2[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        denominat <- Plat0*rrlat0 + Plat1*rrlat1 + Plat2*rrlat2
        P0case <- (Plat0*rrlat0)/denominat  #*# success probabilities for trinomial random variable
        P1case <- (Plat1*rrlat1)/denominat
        P2case <- 1 - P0case - P1case
        
        # Deal with rare crashes of rmultinom due to numerical problems where the
        # program treats probability 0 as a small negative number:
        
        adjustprob <- function(prob) {
          # First break ties:
          if (prob[1]==prob[2] & prob[1]==prob[3] & prob[2]==prob[3]) { prob <- prob+ c(-0.000005,0.000005,0) }
          if (prob[1]==prob[2]) { prob <- prob + c(-0.000005,0.000005,0) }
          if (prob[1]==prob[3]) { prob <- prob + c(-0.000005,0,0.000005) }
          if (prob[2]==prob[3]) { prob <- prob + c(0,-0.000005,0.000005) }
          
          pmin <- min(prob)
          pmax <- max(prob)
          pmiddle <- 1-pmin-pmax
          if (prob[1]==pmin) { prob[1] <- prob[1] + 0.00001 }
          if (prob[2]==pmin) { prob[2] <- prob[2] + 0.00001 }
          if (prob[3]==pmin) { prob[3] <- prob[3] + 0.00001 }
          
          if (prob[1]==pmax) { prob[1] <- prob[1] - 0.00001 }
          if (prob[2]==pmax) { prob[2] <- prob[2] - 0.00001 }
          if (prob[3]==pmax) { prob[3] <- prob[3] - 0.00001 }
          
          prob[prob < 0] <- 0
          return(prob) 
        }
        
        inds <- rmultinom(nCases,1,adjustprob(c(P0case,P1case,P2case)))
        # Computes the numbers of cases in the lo, med, hi latent groups
        nCases0 <- length(inds[1,][inds[1,]==1])
        nCases2 <- length(inds[3,][inds[3,]==1])
        nCases1 <- nCases - nCases0 - nCases2
        N0 <- round(Plat0*N)
        N2 <- round(Plat2*N)
        N1 <- N - N0 - N2
        # Address rounding that could make N1 negative in the dichotomous marker case
        # Keep N fixed at a constant
        if (N1==-1) {
          N0 <- N0 + 1
          N1 <- 0 
        }
        # Also keep  nCases fixed at a constant
        if (nCases1==-1) {
          nCases0 <- nCases0 - 1 
          nCases1 <- 0 
        }
        if (nCases1==1 & N1==0) { 
          nCases0 <- nCases0 + 1 
          nCases1 <- 0 
        }
        #*# total number at risk at tau, separated into cases and controls
        Y <- c(rep(1,nCases0),rep(0,N0-nCases0),rep(1,nCases1),rep(0,N1-nCases1),rep(1,nCases2),rep(0,N - N0 - N1 - nCases2))
        
        
        # Simulate the trinary surrogate with 0,1,2 = lo,med,hi
        # Formulas (9) and (10) in the manuscript:
        
        # First simulate not conditioning on Y=1 or 0 to determine S=0,1,2:
        
        # Given specifications for Spec, FP0, Sens, and FN2 and a logical value indicating if the 
        # biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2),
        # where each subject is assigned a specific level
        AssignBiomarkerLevels <- function(SpecSens, binary){
          Spec <- SpecSens[1]
          Sens <- SpecSens[2]
          FP0 <- SpecSens[3]
          FN2 <- SpecSens[4]
          FP1 <- SpecSens[5]
          FN1 <- SpecSens[6]
          if(binary==TRUE){
            Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP0-Spec,FP0))),
                             rmultinom(N2,1,adjustprob(c(FN2,1-FN2-Sens,Sens))))
          } else{
            Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP0-Spec,FP0))),
                             rmultinom(N1,1,adjustprob(c(FN1,1-FP1-FN1,FP1))),
                             rmultinom(N2,1,adjustprob(c(FN2,1-FN2-Sens,Sens))))
          }
          rownames(Svalues) <- c("S=0","S=1","S=2")
          Svalues <- ifelse(Svalues[1,]==1,0,ifelse(Svalues[2,]==1,1,2))
          return(Svalues)
        }
        
        SpecSens <- cbind(Spec,Sens,FP0,FN2,FP1,FN2)
        
        if (biomType=="binary") { # binary case 
  
          S <- t(apply(SpecSens, 1, function(x) AssignBiomarkerLevels(x, binary=TRUE))) # each row is a set of Sens, Spec, etc. parameters
          
        } else { # trichotomous
          
          S <- t(apply(SpecSens, 1, function(x) AssignBiomarkerLevels(x, binary=FALSE))) # each row is a set of Sens, Spec, etc. parameters
              
        }
        
        # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
        BiomSubset <- function(Y, N, nPhase2, controlCaseRatio, p, cohort){
          
          if (cohort==TRUE) {  # case-cohort sampling design
            
            # subset of subjects with biomarker measured is obtained by drawing a Bernoulli random sample from all at-risk observations 
            # to form the cohort, then augmenting the cohort with all cases
            R <- numeric(length(Y))
            R <- ifelse(rbinom(N, 1, p)==1, 1, R)  # from N, draw Bernoulli sample with sampling probability, p
            R <- ifelse(Y==1, 1, R)  # augment all cases
            keepinds <- c(1:N)[R==1]
            
          } else {  # case-control sampling design
            
            # Keep the S's in nPhase2 of the cases (deleting the rest) and in controlCaseRatio*nPhase2 controls
            casesinds <- c(1:N)[Y==1]
            keepcasesinds <- sample(casesinds,nPhase2,replace=FALSE)
            controlinds <- c(1:N)[Y==0]
            keepcontrolinds <- sample(controlinds,controlCaseRatio*nPhase2,replace=FALSE)
            keepinds <- sort(c(keepcasesinds,keepcontrolinds))
          }
          return(keepinds)
        }
        
        # Those with biomarker data:
        keepinds <- BiomSubset(Y, N, nPhase2, controlCaseRatio, p, cohort)
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(Sens)
  
        
        ##############################################################
        # Now analyze with osDesign
        # (first check if there are 'zeros', in which case Fisher's exact test for the lo vs. hi categories is used. Otherwise,
        # osDesign logistic regression with the pseudo-likelihood method is used as an ordered score test
        
        for(k in 1:nrow(Scc)){
          lodim <- dim(table(Ycc,Scc[k,]))[2]<2 #*# check there are at least two biomarker categories (columns)
          zerosflag <-  lodim
          if (dim(table(Ycc,Scc[k,]))[2]==3) { #*# check if any categories have zero entries
            zerosflag <- table(Ycc,Scc[k,])[1,1]==0 | table(Ycc,Scc[k,])[1,2]==0 | table(Ycc,Scc[k,])[1,3]==0 | table(Ycc,Scc[k,])[2,1]==0 | table(Ycc,Scc[k,])[2,2]==0 | table(Ycc,Scc[k,])[2,3]==0 
          }
          
          if (zerosflag) {
            if (lodim) { pval <- 1}
            if (!lodim) { #*# there are zeros, so Fisher's exact test is used
              pval <- fisher.test(table(Ycc,Scc[k,])[,c(1,dim(table(Ycc,Scc[k,]))[2])])$p.value 
            }
            if (pval <= alpha & length(Ycc[Scc[k,]==2&Ycc==1])/length(Scc[k,][Scc[k,]==2]) < length(Ycc[Scc[k,]==0&Ycc==1])/length(Scc[k,][Scc[k,]==0])) {
              powerstrinary[k,j] <- powerstrinary[k,j] + 1
            }
          }
          
          if (!zerosflag) {
            fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
            pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[k,j] <- powerstrinary[k,j] + 1}
          }
        }
      }
    } else if (biomType=="continuous") {  
    
      # Continuous biomarker
      
      # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
      # using the function risk1cont() above
      
      for (j in 1:o) {
        beta <- truebetas[j]
        alphalat <- alphalatvect[j]
        
        # These simulations condition on n (i.e., number of infections in vaccine arm) and
        # also on the number of controls fixed at controlCaseRatio*n
        # This matches what was done for a binary correlate
        
        # Arbitrarily put the cases first and controls second
        # The numbers of cases and controls are fixed, e.g., a typical retrospective design
        Y <- c(rep(1,nCases),rep(0,N-nCases))
        
        # Compute the denominator of the density of X|Y=1 when Y|X follows a logistic regression model with the truncated part
        # associated with VElowestvect and X is normal
        # with mean zero and standard deviation sqrt(sigma2tr)
        #
        # rhos[1]:
        for(k in 1:length(nus)){
          f <- function(x) {
            ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
            return(ans)
          }
          denomdensityXcases <- integrate(f,lower=nus[k],upper=5)$value
          denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowestvect[j])*risk0
          
          numerdensXcases <- function(x) {
            num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
            num[x <= nus[k]] <- PlatVElowest*(1-VElowestvect[j])*risk0
            return(num)
          }
          
          numerdensXcontrols <- function(x) {
            num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr[k]))
            num[x <= nus[k]] <- PlatVElowest*(1-(1-VElowestvect[j])*risk0)
            return(num)
          }
          
          Xpoints <- seq(-3.5,3.5,len=25000)
          probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
          probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)
          
          Xcases <-    sample(Xpoints,size=nCases,prob=probscases,replace=TRUE)
          Xcontrols <- sample(Xpoints,size=N-nCases,prob=probscontrols,replace=TRUE)
          X <- c(Xcases,Xcontrols)
        }  
        
        # Create the 4 immune response variables for the 4 degrees of measurement error
        error <- t(sapply(sigma2e, function(x) rnorm(N,mean=0,sd=sqrt(x))))
        S <- X + error
   
        # Those with data:
        keepinds <- BiomSubset(Y, N, nPhase2, controlCaseRatio, p, cohort)
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rhos)

        for(k in 1:nrow(Scc)){
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
          pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}
        }
      }
    }    
  }
  
  powerstrinary <- powerstrinary/M
  powerscont <- powerscont/M

  if (o > 0)
    ans <- list(RRlat2,t(powerstrinary),exp(truebetas),t(powerscont))
  else
    ans <- list(RRlat2,t(powerstrinary))
  write(N,file="sampsizeALL.dat")
  write(nCases,file="numbeventsALL.dat")
  write(nPhase2,file="numbeventsPhase2.dat")
  write(1-RRoverall,file="VEoverallCoRpower.dat")
  write(alpha,file="alpha.dat")
  write(rhos,file="rhos.dat",ncolumns=1,append=FALSE)
  write(RRlat2,file="RRlat2.dat",ncolumns=1,append=FALSE)
  write(RRlat0,file="RRlat0.dat",ncolumns=1,append=FALSE)
  write(powerstrinary,file=paste("powerstrinary",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=4,append=FALSE)
  # RRs the relative risks that are the effect sizes RR_c that
  # need to be on the x-axis of powerplots
  if(o>1){
    write(exp(truebetas),file="RRs.dat",ncolumns=1,append=FALSE)
    write(powerscont,file=paste("powerscont",controlCaseRatio,".dat",sep=""),ncolumns=4,append=FALSE)
    write(PlatVElowest,file="PlatVElowest.dat")
    write(VElowestvect,file="VElowestvect.dat",ncolumns=1,append=FALSE)
    write(truebetas,file="truebetas.dat")
  }
  write(controlCaseRatio,file="controlCaseRatio.dat")
  write(P2,file="P2.dat")
  
  # Print out the CoR effect sizes
  write(risk1_0,file=paste("vaccineriskslo",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=4,append=FALSE)
  write(risk1_2,file=paste("vaccineriskshi",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=4,append=FALSE)
  
  # write out alpha intercept as logit(Y=1|s=0) for trinary/binary case
  write(c(t(logit(risk1_0))), file="trinaryalpha.dat",ncolumns=1,append=FALSE)
  # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/binary case
  write(c(t(logit(risk1_2)-logit(risk1_0))), file="trinarybeta.dat",ncolumns=1,append=FALSE)
  
  return(ans)
  
}


# Repeat the compute power function with a variation in the variable input

# Inputs:
# Range of fixed sample sizes
# VE, risk0, P^{lat}_{0} (a single value)
# Single values for VE^{lat}_{0}, VE^{lat}_{1}, VE^{lat}_{2}
# One value each of P_{0} and P_{2}
# Other analyses could specify alternative values of P_{0} and P_{2} to study how power varies with
# the selected thresholds for a high and low response.

# Determined/computed from formulas in the article:
# P^{lat}_{1}, P^{lat}_{2}

# As for computepower(),
# following Step 7 in the article, the program for a trichotomous marker accounts for assay noise in one of two ways.
# The first approach specifies Spec, Sens, FP0, FN2 which determine FP1 and FN2 from equations (8) and (9).  4 settings
# must be specified.
#
# The second approach specifies sigma2obs and rho, again requiring 4 settings for rho; this approach assumes the normal
# measurement error model (4) in the article.  Specifying
# the vectors Spec=1, Sens=1, FP0=1, FN2=1 defaults to approach 2, which is used in the manuscript.

# Output: Power

#' @describeIn computepower Vectors are inputted to repeat power calculations at multiple sample sizes.
#' @param numAtRiskTauCasesVect  Vector of the number of cases in the vaccine group at-risk at tau and a case by taumax (regardless of whether the biomarker is measured).
#' @param numAtRiskTauCasesPhase2Vect Vector of the number of cases in the vaccine group at-risk at tau and a case by taumax and with the biomarker measured (i.e., in Phase 2).
#' @param numAtRiskTauControlsVect Vector of the number of controls in the vaccine group at-risk at tau and not diseased at the end of follow-up taumax (regardless of whether the biomarker is measured).
#' @param RRlat0point For a trichotomous biomarker, one minus VE for the latent low biomarker subgroup.
#' @param RRlat1point For a trichotomous biomarker, one minus VE for the latent middle biomarker subgroup.
#' @param RRlat2point For a trichotomous biomarker, one minus VE for the latent high biomarker subgroup.
#' @param VElowest The lowest possible value of vaccine efficacy
#'
#' @export
computepower.n <- function(numAtRiskTauCasesVect, numAtRiskTauCasesPhase2Vect, numAtRiskTauControlsVect,
                           risk0, RRoverall,
                           Plat0, P2, P0,
                           RRlat2point, RRlat1point, RRlat0point,
                           PlatVElowest, VElowest,
                           controlCaseRatio=5,
                           M=100,
                           alpha=0.05,
                           sigma2obs=1,
                           rhos=c(1,.9,.7,.5),
                           Spec=rep(1,4), FP0=rep(0,4), Sens=rep(1,4), FN2=rep(0,4),
                           tpsMethod=c("PL", "ML","WL"),
                           biomType=c("continuous", "trichotomous", "binary"),
                           cohort=FALSE, p=NULL) {
  
  # Like computepower except RRlat2 (now RRlat2point) and VElowest are singulars, and
  # numAtRiskTauCasesVect and numAtRiskTauControlsVect are vectors
  # reflecting different sample sizes
  
  # check sampling design input parameters are specified and valid
  if(cohort==TRUE) {  #case-cohort 
    
    if (is.null(p)==TRUE) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is unspecified")
    } else if (p < 0 | p > 1) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is not a valid probability")
    }
    
  } else if (is.null(controlCaseRatio)==TRUE) {  #case-control
    stop("Case-control sampling was chosen and controlCaseRatio is unspecified")
  }
  
  # check biomarker type and input parameters match
  if(biomType=="binary" & (P0+P2 != 1)){
    stop("Binary biomarker was specified but P0 and P2 do not add up to 1")
  }
  if(biomType=="continuous" & is.null(VElowestvect)){
    stop("Continuous biomarker was specified but VElowestvect is NULL")
  }
  tpsMethod <- match.arg(tpsMethod, choices = c("PL","ML","WL"))
  biomType <- match.arg(biomType, choices = c("continuous", "trichotomous", "binary"))
  
  
  # Compute Plat2:
  VEoverall <- 1 - RRoverall
  VElat2point <- 1 - RRlat2point
  VElat1point <- 1 - RRlat1point
  VElat0point <- 1 - RRlat0point
  num <- VEoverall - VElat1point - Plat0*(VElat0point - VElat1point)
  den <- VElat2point - VElat1point
  Plat2 <- num/den
  Plat1 <- 1 - Plat2 - Plat0
  P1 <- 1 - P0 - P2
  
  nCasesVect <- numAtRiskTauCasesVect
  nPhase2Vect <- numAtRiskTauCasesPhase2Vect
  # Overall denominator: number observed to be at risk when the immune response is measured:
  NVect <- nCasesVect + numAtRiskTauControlsVect
  
  sigma2e <- (1-rhos)*sigma2obs

  sigma2tr <- rhos*sigma2obs
  
  if(biomType=="trichotomous" | biomType=="binary") {
    ###############################################################
    # Computations for a trinary biomarker
    
    Approach2 <- (all(is.null(Spec), is.null(Sens), is.null(FP0), is.null(FN2)))
    
    if (Approach2) {
      # Default choice
      
      # Compute Sens, Spec, FP0, FP1, FN2, FN2
      ans <- computeSensSpecFPFN(sigma2obs,rhos,Plat0,Plat2,P0,P2)
      Sens <- unlist(lapply(ans, function(x) x[[1,10]])) 
      Spec <- unlist(lapply(ans, function(x) x[[1,11]]))
      FP0 <- unlist(lapply(ans, function(x) x[[1,12]])) 
      FP1 <- unlist(lapply(ans, function(x) x[[1,13]])) 
      FN2 <- unlist(lapply(ans, function(x) x[[1,14]]))
      FN1 <- unlist(lapply(ans, function(x) x[[1,15]]))
    }
    
    # Approach 1 in the manuscript:
    if (!Approach2) {
      
      # Apply formula (8) in the manuscript
      FN1 <- (P0 - Spec*Plat0 - FN2*Plat2)/Plat1   #*#P0, Plat0, Plat2 given params
      
      # Apply formula (9) in the manuscript
      FP1 <- (P2 - Sens*Plat2 - FP0*Plat0)/Plat1
      
      #Check if an error in the ranges of values due to an out of
      # bounds input parameter
      if (any(FN1 < 0 | FN1 > 1 | FP1 < 0 | FP1 > 1)){
        stop("Approach 1 was used and one of the parameters Sens, Spec, FP0, FN2 is out of range")
      }
    }
    
    # Binary biomarker special case (to remove small values of P1x)
    if (biomType=="binary") {
      P1 <- 0
      P2 <- 1 - P0
    }
    
    ################## OLD on the chopping block
    ## Compute the marginal risks:
    ## Made it to the end of follow-up HIV negative
    #risk1 <- RRoverall*risk0
    #
    ## Observed risks P(Y(1)=1|S(1)=lo, med, or hi)
    #risk1hi1 <- (RRlat0point*FP01 + RRlat1point*FP11 + RRlat2point*Sens1)*risk0
    #risk1lo1 <- (RRlat0point*Spec1 + RRlat1point*FN11 + RRlat2point*FN21)*risk0
    #risk1hi2 <- (RRlat0point*FP02 + RRlat1point*FP12 + RRlat2point*Sens2)*risk0
    #risk1lo2 <- (RRlat0point*Spec2 + RRlat1point*FN12 + RRlat2point*FN22)*risk0
    #risk1hi3 <- (RRlat0point*FP03 + RRlat1point*FP13 + RRlat2point*Sens3)*risk0
    #risk1lo3 <- (RRlat0point*Spec3 + RRlat1point*FN13 + RRlat2point*FN23)*risk0
    #risk1hi4 <- (RRlat0point*FP04 + RRlat1point*FP14 + RRlat2point*Sens4)*risk0
    #risk1lo4 <- (RRlat0point*Spec4 + RRlat1point*FN14 + RRlat2point*FN24)*risk0
    
    # Compute the marginal risks:
    # Made it to the end of follow-up HIV negative
    risk1 <- RRoverall*risk0
    
    # Observed risks P(Y(1)=1|S(1)=0, 1, or 2)  #*# for diff values of rho; using Bayes' rule
    
    probX0_cond_S2 <- FP0*Plat0/P2
    probX1_cond_S2 <- FP1*Plat1/P2
    probX2_cond_S2 <- Sens*Plat2/P2
    risk1_2 <- (probX0_cond_S2 * RRlat0point + probX1_cond_S2 * RRlat1point + probX2_cond_S2 * RRlat2point)*risk0 #vector with length=length(rhos)
    probX0_cond_S0 <- Spec*Plat0/P0
    probX1_cond_S0 <- FN1*Plat1/P0
    probX2_cond_S0 <- FN2*Plat2/P0
    risk1_0 <- (probX0_cond_S0 * RRlat0point + probX1_cond_S0 * RRlat1point + probX2_cond_S0 * RRlat2point)*risk0
    risk1_1 <- (risk1 - risk1_0*P0 - risk1_2*P2)/P1
    
    # Note: For the binary biomarker special case, the risk1medx are NA
    #       They are never used so it is irrelevant
    
    es <- risk1_2/risk1_0
    
    # Vaccine risks within the latent subgroups (independent of rho of course)
    risk1lat_2 <- RRlat2point*risk0
    risk1lat_1 <- RRlat1point*risk0
    risk1lat_0 <- RRlat0point*risk0
    #################################################
  }  else if (biomType=="continuous") {  
  
    #################################################
    # Computations for a continuous biomarker
    # Define the truebeta indexed by the user-specified vector VElowest
      
    nus <- sqrt(rhos*sigma2obs)*qnorm(PlatVElowest)
    
    risk1latnu <- (1-VElowest)*risk0
    
    f <- function(alpha) {
      
      g <- function(x) {
        rho <- 1
        piece1 <- exp(alpha*(1 - x/nus[1]))*(risk1latnu^(x/nus[1]))
        piece2 <- (1-risk1latnu)^(x/nus[1]) + piece1
        piece3 <- dnorm(x/(sqrt(rho*sigma2obs)))
        kernel <- (piece1/piece2)*piece3
        return(kernel) 
      }
      
      # nus[1] corresponds to rho=1
      logitterm <- integrate(g,lower=nus[1],upper=6)$value
      
      ans <- 1-VEoverall - (PlatVElowest*risk1latnu + logitterm)/risk0
      return(ans) 
    }
    
    alphalat <- uniroot(f,lower=-10,upper=10)$root
    
    # Second solve for betalat:
    D <- risk1latnu
    truebeta <- (log(D/(1-D)) - alphalat)/nus[1]
    
  }
  
  # Function for computing the infection probabilities of vaccinees
  risk1cont <- function(x,alphalat,betalat) {
    linpart <- alphalat + betalat*x
    ans <- exp(linpart)/(1+exp(linpart))
    return(ans) 
  }
  ###################################################
  
  powerstrinary <- matrix(0, nrow=length(rhos), ncol=length(NVect))
  powerscont <- matrix(0, nrow=length(rhos), ncol=length(NVect))
  
  for (i in 1:M) {
    
    if(biomType=="trichotomous" | biomType=="binary"){
    
      # Trinary biomarker:
      for (j in 1:length(NVect)) {
        # Fix the number of cases and controls, putting the cases first and controls
        # second for each subgroup:
        N <- NVect[j]
        nCases <- nCasesVect[j]
        nPhase2 <- nPhase2Vect[j]
        
        rrlat0 <- risk1lat_0/(risk1lat_0+risk1lat_1+risk1lat_2)
        rrlat1 <- risk1lat_1/(risk1lat_0+risk1lat_1+risk1lat_2)
        rrlat2 <- risk1lat_2/(risk1lat_0+risk1lat_1+risk1lat_2)
        denominat <- Plat0*rrlat0 + Plat1*rrlat1 + Plat2*rrlat2
        P0case <- (Plat0*rrlat0)/denominat
        P1case <- (Plat1*rrlat1)/denominat
        P2case <- 1 - P0case - P1case
        
        # Deal with rare crashes of rmultinom due to numerical problems where the
        # program treats probability 0 as a small negative number:
        
        adjustprob <- function(prob) {
          
          # First break ties:
          if (prob[1]==prob[2] & prob[1]==prob[3] & prob[2]==prob[3]) { prob <- prob + c(-0.000005,0.000005,0) }
          if (prob[1]==prob[2]) { prob <- prob + c(-0.000005,0.000005,0) }
          if (prob[1]==prob[3]) { prob <- prob + c(-0.000005,0,0.000005) }
          if (prob[2]==prob[3]) { prob <- prob + c(0,-0.000005,0.000005) }
          
          pmin <- min(prob)
          pmax <- max(prob)
          pmiddle <- 1-pmin-pmax
          if (prob[1]==pmin) { prob[1] <- prob[1] + 0.00001 }
          if (prob[2]==pmin) { prob[2] <- prob[2] + 0.00001 }
          if (prob[3]==pmin) { prob[3] <- prob[3] + 0.00001 }
          
          if (prob[1]==pmax) { prob[1] <- prob[1] - 0.00001 }
          if (prob[2]==pmax) { prob[2] <- prob[2] - 0.00001 }
          if (prob[3]==pmax) { prob[3] <- prob[3] - 0.00001 }
          
          prob[prob < 0] <- 0
          return(prob) 
        }
        
        inds <- rmultinom(nCases,1,adjustprob(c(P0case,P1case,P2case)))
        #inds <- rmultinom(nCases,1,c(Plat0*risk1lat_0/(Plat0*risk1lat_0+Plat1*risk1lat_1+Plat2*risk1lat_2),Plat1*risk1lat_1/(Plat0*risk1lat_0+Plat1*risk1lat_1+Plat2*risk1lat_2),Plat2*risk1lat_2/(Plat0*risk1lat_0+Plat1*risk1lat_1+Plat2*risk1lat_2)))
        
        # Number of cases in the lo, med, hi latent groups
        nCases0 <- length(inds[1,][inds[1,]==1])
        nCases2 <- length(inds[3,][inds[3,]==1])
        nCases1 <- nCases - nCases0 - nCases2
        N0 <- round(Plat0*N)
        N2 <- round(Plat2*N)
        N1 <- N - N0 - N2
        cat(paste("1 N1 = ",N1),"\n")
        cat(paste("1 nCases1 = ",nCases1),"\n")
        # Address rounding that could make N1 negative in the dichotomous marker case
        # Keep N fixed at a constant
        if (N1==-1) {
          N0 <- N0 + 1
          N1 <- 0 
        }
        # Also keep nCases fixed at a constant
        if (nCases1==-1) {
          nCases0 <- nCases0 - 1 
          nCases1 <- 0 
        }
        if (nCases1==1 & N1==0) { 
          nCases0 <- nCases0 + 1 
          nCases1 <- 0 
        }
        cat(paste("2 N1 = ",N1),"\n")
        cat(paste("2 nCases1 = ",nCases1),"\n")
        Y <- c(rep(1,nCases0),rep(0,N0-nCases0),rep(1,nCases1),rep(0,N1-nCases1),rep(1,nCases2),rep(0,N - N0 - N1 - nCases2))
        
        # Simulate the trinary surrogate with 0,1,2 = lo,med,hi
        # Formulas (12) and (13) in the manuscript:
        
        # Given specifications for Spec, FP0, Sens, and FN2 and a logical value indicating if the 
        # biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2),
        # where each subject is assigned a specific level
        AssignBiomarkerLevels <- function(SpecSens, binary){
          Spec <- SpecSens[1]
          Sens <- SpecSens[2]
          FP0 <- SpecSens[3]
          FN2 <- SpecSens[4]
          FP1 <- SpecSens[5]
          FN1 <- SpecSens[6]
          if(binary==TRUE){
            Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP0-Spec,FP0))),
                             rmultinom(N2,1,adjustprob(c(FN2,1-FN2-Sens,Sens))))
          } else{
            Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP0-Spec,FP0))),
                             rmultinom(N1,1,adjustprob(c(FN1,1-FP1-FN1,FP1))),
                             rmultinom(N2,1,adjustprob(c(FN2,1-FN2-Sens,Sens))))
          }
          rownames(Svalues) <- c("S=0","S=1","S=2")
          Svalues <- ifelse(Svalues[1,]==1,0,ifelse(Svalues[2,]==1,1,2))
          return(Svalues)
        }
        
        SpecSens <- cbind(Spec,Sens,FP0,FN2,FP1,FN1)
        
        if (biomType=="binary") { # binary case only
          
          S <- t(apply(SpecSens, 1, function(x) AssignBiomarkerLevels(x, binary=TRUE))) # each row is a set of Sens, Spec, etc. parameters
        
        } else { # trichotomous
          
          S <- t(apply(SpecSens, 1, function(x) AssignBiomarkerLevels(x, binary=FALSE))) # each row is a set of Sens, Spec, etc. parameters
        
        }
        
        # Those with data:
        keepinds <- BiomSubset(Y, N, nPhase2, controlCaseRatio, p, cohort)
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1,function(x) x[keepinds])) #nrow=length(rhos)
        
        ##############################################################
        # Now analyze with osDesign
        # (first check if there are 'zeros', in which case Fisher's exact test for the lo vs. hi categories is used. Otherwise,
        # osDesign logistic regression is used as an ordered score test
        
        for(k in 1:nrow(Scc)){
          lodim <- dim(table(Ycc,Scc[k,]))[2]<2 #*# check there are at least two biomarker categories (columns)
          zerosflag <-  lodim
          if (dim(table(Ycc,Scc[k,]))[2]==3) { #*# check if any categories have zero entries
            zerosflag <- table(Ycc,Scc[k,])[1,1]==0 | table(Ycc,Scc[k,])[1,2]==0 | table(Ycc,Scc[k,])[1,3]==0 | table(Ycc,Scc[k,])[2,1]==0 | table(Ycc,Scc[k,])[2,2]==0 | table(Ycc,Scc[k,])[2,3]==0 
          }
          
          if (zerosflag) {
            if (lodim) { pval <- 1}
            if (!lodim) { #*# there are zeros, so Fisher's exact test is used
              pval <- fisher.test(table(Ycc,Scc[k,])[,c(1,dim(table(Ycc,Scc[k,]))[2])])$p.value 
            }
            if (pval <= alpha & length(Ycc[Scc[k,]==2&Ycc==1])/length(Scc[k,][Scc[k,]==2]) < length(Ycc[Scc[k,]==0&Ycc==1])/length(Scc[k,][Scc[k,]==0])) {
              powerstrinary[k,j] <- powerstrinary[k,j] + 1
            }
          }
          
          if (!zerosflag) {
            fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
            pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[k,j] <- powerstrinary[k,j] + 1}
          }
        }
      }
    } else if (biomType=="continuous") {  
    
      # Continuous biomarker

      # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
      # using the function risk1cont() above
      for (j in 1:length(NVect)) {
        N <- NVect[j]
        nCases <- nCasesVect[j]
        nPhase2 <- nPhase2Vect[j]
        beta <- truebeta
        
        # These simulations condition on n (i.e., number of infections in vaccine arm) and
        # also on the number of controls fixed at controlCaseRatio*n
        # This matches what was done for a binary correlate
        
        # Arbitrarily put the cases first and controls second
        # The numbers of cases and controls are fixed, e.g., a typical retrospective design
        Y <- c(rep(1,nCases),rep(0,N-nCases))
        
        # Compute the denominator of the density of X|Y=1 when Y|X follows a log. regr model with the truncated part
        # associated with VElowestvect and X is normal
        # with mean zero and standard deviation sqrt(sigma2tr)
        #
        # rhos[1]:
        for(k in 1:length(nus)){
          f <- function(x) {
            ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
            return(ans)
          }
          denomdensityXcases <- integrate(f,lower=nus[k],upper=5)$value
          denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowestvect[j])*risk0
          
          numerdensXcases <- function(x) {
            num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
            num[x <= nus[k]] <- PlatVElowest*(1-VElowestvect[j])*risk0
            return(num)
          }
          numerdensXcontrols <- function(x) {
            num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr[k]))
            num[x <= nus[k]] <- PlatVElowest*(1-(1-VElowestvect[j])*risk0)
            return(num)
          }
          
          Xpoints <- seq(-3.5,3.5,len=25000)
          probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
          probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)
          
          Xcases <-    sample(Xpoints,size=nCases,prob=probscases,replace=TRUE)
          Xcontrols <- sample(Xpoints,size=N-nCases,prob=probscontrols,replace=TRUE)
          X <- c(Xcases,Xcontrols)
        }
        
        # Create the 4 immune response variables for the 4 degrees of measurement error
        error <- t(sapply(sigma2e, function(x) rnorm(N,mean=0,sd=sqrt(x))))
        S <- X + error
        
        # Those with data:
        keepinds <- BiomSubset(Y, N, nPhase2, controlCaseRatio, p, cohort)
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rhos)
        
        for(k in 1:nrow(Scc)){
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
          pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}
        }
      }
    }  
  }
  
  powerstrinary <- powerstrinary/M
  powerscont <- powerscont/M
  
  ans <- list(RRlat0point,RRlat2point,t(powerstrinary),t(powerscont))
  write(rhos,file="rhos.dat",ncolumns=1,append=FALSE)
  write(1-RRoverall,file="VEoverallpower2.dat")
  write(alpha,file="alpha2.dat")
  write(controlCaseRatio,file="controlCaseRatio.dat")
  write(1-RRlat2point,file=paste("RRlat2point",P2,".dat",sep=""))
  write(1-RRlat0point,file=paste("RRlat0point",P2,".dat",sep=""))
  write(P2,file=paste("P2",P2,".dat",sep=""))
  write(numAtRiskTauCasesVect,"samplesizescasesALL.dat",ncolumns=1,append=FALSE)
  write(numAtRiskTauCasesPhase2Vect,"samplesizescases2Phase.dat",ncolumns=1,append=FALSE)
  write(Plat2,file=paste("Plat2point",P2,".dat",sep=""))
  write(PlatVElowest,file="PlatVElowest.dat")
  write(VElowest,file="VElowest.dat")
  write(truebeta,file="truebeta.dat")
  write(risk1lat_2,file="risk1lat_2.dat")
  write(risk1lat_1,file="risk1lat_1.dat")
  write(risk1lat_0,file="risk1lat_0.dat")
  write(powerstrinary,file=paste("powerstrinary",1-RRlat2point,1-RRlat0point,P2,P0,controlCaseRatio,"2ss.dat",sep=""),ncolumns=4,append=FALSE)
  write(powerscont,file="powerscont2ss.dat",ncolumns=4,append=FALSE)
  # write out alpha intercept as logit(Y=1|s=0) for trinary/binary case
  write(c(t(logit(risk1_0))), file="trinaryalpha.dat",ncolumns=1,append=FALSE)
  # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/binary case
  write(c(t(logit(risk1_2)-logit(risk1_0))), file="trinarbeta.dat",ncolumns=1,append=FALSE)
  
  return(ans)
  
}


#######################################################
### computeSensSpecFPFN() function from Peter       ###
### extracted from runCoRpower_trinary_manuscript.R ###
#######################################################


#############
# computeSensSpecFPFN is a function for mapping input parameters to Sensitivity and Specificity,
# FP0, FP1, FN2, FN1 (defined above)
# S a trichotomous biomarker S = 2 if S* > tauhi for a fixed tauhi and S=0 if S* <= taulo for a fixed taulo,
# and S = 1 if S* is in between taulo and tauhi.
# X* is the 'true' biomarker with X* = 2 if X* > thetahiVE and X* = 0 if X* <= thetaloVE
# for fixed thetaloVE and thetahiVE that are solved for.
#
# Variance of S* = sigma2obs
# from a classical measurement error model S* = X* + e  where e~N(0,sigma2e), X* ~ N(0,sigma2tr)
#
# rho is the protection relevant fraction of the variability of S*
# rho = 1 - sigma2e/sigma2obs = sigma2tr/sigma2ob
# sigma2tr = Var(X*) rho*sigma2obs
#
# This function also applies for a binary biomarker in which case only Sensitivity and Specificity
# are relevant (FP0, FP1, FN2, FN1 are not used in the calculations)
##################################################

# Original function that requires Plat2 to be a scalar
computeSensSpecFPFN <- function(sigma2obs,rhos,Plat0,Plat2,P0,P2) {
  # sigma2tr = Var(X) = Var(Str) = rho*sigma2obs
  # If Plat0 + Plat2 = 1 then the method collapses to a binary biomarker,
  # and FP0, FP1, FN2, FN1 are irrelevant; this function simply returns 0's in that scenario
  
  # P2 may be a scalar or a vector, which should include one value equal to
  #      Plat2 and values straddling either side
  # P0 may be a scalar or a vector, which should include one value equal to
  #      Plat0 and values straddling either side
  
  sigma2e <- (1-rhos)*sigma2obs
  sigma2tr <- rhos*sigma2obs
  thetahiVE <- qnorm(1-Plat2)*sqrt(sigma2tr)
  thetaloVE <- qnorm(Plat0)*sqrt(sigma2tr)
  
  Plat1 <- 1 - Plat0 - Plat2
  m <- length(P2)
  ans <- list()
  
  for(i in 1:length(rhos)){
    Sens <- rep(1,m)
    Spec <- rep(1,m)
    FP0 <- rep(0,m)
    FP1 <- rep(0,m)
    FN2 <- rep(0,m)
    FN1 <- rep(0,m)
    tauhisolution <- rep(0,m)
    taulosolution <- rep(0,m)
    if (rhos[i] < 1) {  #*# if rho=1, then Sens=1, Spec=1, FP0=0, FP1=0, FN2=0, FN1=0
      # Stochastic integration
      set.seed(1)
      X <- rnorm(20000,0,sqrt(sigma2tr[i]))
      set.seed(2)
      S <- X + rnorm(20000,0,sqrt(sigma2e[i]))
      
      Phi <- sum(X>thetahiVE[i])/length(X)
      Plo <- sum(X<=thetaloVE[i])/length(X)
      Pmed <- 1 - Phi - Plo
      
          # # Compute a kernel over a grid of tau values:
          # n <- 10000
          # tauhi <- seq(-2.5,2.5,len=n)
          # taulo <- seq(-2.5,2.5,len=n)
          # Sensvec <- rep(1,n)
          # Specvec <- rep(1,n)
          # FP0vec <- rep(0,n)
          # FP1vec <- rep(0,n)
          # FN2vec <- rep(0,n)
          # FN1vec <- rep(0,n)
          # 
          # for (j in 1:length(tauhi)) {
          #   Sensvec[j] <- (sum(S>tauhi[j] & X > thetahiVE[i])/length(S))/Phi
          #   Specvec[j] <- (sum(S<=taulo[j] & X <= thetaloVE[i])/length(S))/Plo
          #   if (Pmed==0) { #*# if binary biomarker, 0's for FP1, FP0, FN2, FN1
          #     FP1vec[j] <- 0
          #     FP0vec[j] <- 0
          #     FN2vec[j] <- 0
          #     FN1vec[j] <- 0 }
          #   if (Pmed > 0) {
          #     FP1vec[j] <- (sum(S>tauhi[j] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
          #     FP0vec[j] <- (sum(S>tauhi[j] & X <= thetaloVE[i])/length(S))/Plo
          #     FN2vec[j] <- (sum(S<=taulo[j] & X > thetahiVE[i])/length(S))/Phi
          #     FN1vec[j] <- (sum(S<=taulo[j] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
          #   }
          # }
      
      for (l in 1:m) {
        
        # Find the cut points tauhi and taulo by solving the following equations:
        #   0 = Sensvec*Plat2 + FP1vec*Plat1 + FP0vec*Plat0 - P2  (eqn 8 in manuscript)
        #   0 = Specvec*Plat0 + FN1vec*Plat1 + FN2vec*Plat2 - P0  (eqn 7 in manuscript)
        # where 
        #   Sensvec <- (sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi
        #   Specvec <- (sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo
        # if binary biomarker, 
        #   FP1vec <- 0
        #   FP0vec <- 0
        #   FN2vec <- 0
        #   FN1vec <- 0
        # if trichotomous biomarker,
        #   FP1vec <- (sum(S>tauhi & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        #   FP0vec <- (sum(S>tauhi & X <= thetaloVE[i])/length(S))/Plo
        #   FN2vec <- (sum(S<=taulo & X > thetahiVE[i])/length(S))/Phi
        #   FN1vec <- (sum(S<=taulo & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        
        if (Pmed==0){  # binary 
          f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi)*Plat2 - P2[l]
          f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo)*Plat0 - P0[l]
        } else {  # trichotomous
          f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi)*Plat2 +
                                ((sum(S>tauhi & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed)*Plat1 +
                                ((sum(S>tauhi & X <= thetaloVE[i])/length(S))/Plo)*Plat0 - P2[l]
          f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo)*Plat0 +
                                ((sum(S<=taulo & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed)*Plat1 +
                                ((sum(S<=taulo & X > thetahiVE[i])/length(S))/Phi)*Plat2 - P0[l]
        }
        tauhisol <- uniroot(f2, interval=c(-2.5,2.5))$root
        taulosol <- uniroot(f0, interval=c(-2.5,2.5))$root
        
            # kernelhi <- Sensvec*Plat2 + FP1vec*Plat1 + FP0vec*Plat0 - P2[l]
            # kernelhi <- kernelhi^2
            # 
            # kernello <- Specvec*Plat0 + FN1vec*Plat1 + FN2vec*Plat2 - P0[l]
            # kernello <- kernello^2
            # 
            # indhi <- c(1:length(kernelhi))[kernelhi==min(kernelhi)] #*# equiv to which(kernelhi==min(kernelhi))
            # tauhisol <- tauhi[indhi]
            # indlo <- c(1:length(kernello))[kernello==min(kernello)]
            # taulosol <- taulo[indlo]
            # 
            # if (length(tauhisol)>1) {tauhisol <- tauhisol[1]} #*# if more than one tau given, choose first one
            # if (length(taulosol)>1) {taulosol <- taulosol[1]}
        
        tauhisolution[l] <- tauhisol
        taulosolution[l] <- taulosol
        
        Sens[l] <- sum(S>tauhisolution[l] & X > thetahiVE[i])/sum(X>thetahiVE[i])
        Spec[l] <- sum(S<=taulosolution[l] & X <= thetaloVE[i])/sum(X<=thetaloVE[i])
        if (Pmed==0) {  #*# if binary biomarker, 0's for FP1, FP0, FN2, FN1
          FP1[l] <- 0
          FP0[l] <- 0
          FN2[l] <- 0
          FN1[l] <- 0 }
        if (Pmed > 0) {
          FP1[l] <- (sum(S>tauhisolution[l] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
          FP0[l] <- (sum(S>tauhisolution[l] & X <= thetaloVE[i])/length(S))/Plo
          FN2[l] <- (sum(S<=taulosolution[l] & X > thetahiVE[i])/length(S))/Phi
          FN1[l] <- (sum(S<=taulosolution[l] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        }
      }
    }
    ans[[i]] <- cbind(rep(thetaloVE[i],m),rep(thetahiVE[i],m),rep(Plat0,m),rep(Plat1,m),rep(Plat2,m),P0,P2,
                      taulosolution,tauhisolution,Sens,Spec,FP0,FP1,FN2,FN1)
  }  
  return(ans)
}
