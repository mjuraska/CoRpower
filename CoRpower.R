# Main program computepower() for correlate of risk (CoR) power calculations and
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

# sensitivity, specificity, false negatives and false positives
# sens \equiv P(S(1)=2|X=2) (==1 for S with no measurement error)
# spec \equiv P(S(1)=0|X=0) (==1 for S with no measurement error)
# FP^0 \equiv P(S(1)=2|X=0) (==0 for S with no measurement error)
# FN^2 \equiv P(S(1)=0|X=2)  (==0 for S with no measurement error)
# FP^1 \equiv P(S(1)=2|X=1) (==0 for S with no measurement error)
# FN^1 \equiv P(S(1)=0|X=1) (==0 for S with no measurement error)

########################################################################
# End of unique notations for trinary immune response S = 0,1,2
# The remaining notations below are the same for trinary (Approach 2) and continous S
# except that PlatVElowest and VElowest only apply to continous S
#######################################################################

# The biomarker is measured at tau months and subjects are followed through taumax months

# RRoverall = 1 - overall vaccine efficacy
# annincinfectionplac = annual HIV infection incidence in placebo group
# annincdropout = annual dropout rate assumed the same in both groups

# nAtRiskTauCases     All cases in the vaccine group at-risk at tau and a case by taumax
#                             (regardless of whether the biomarker is measured)
# nAtRiskTauControls  All controls in the vaccine group at-risk at tau and not diseased at the end of follow-up taumax
#                             (regardless of whether the biomarker is measured)
# nAtRiskTauCasesPhase2  As above and also have the biomarker measured (i.e., in Phase 2)

# sigma2obs  observed variance of the continuous marker S*
# rho       vector of rho, the proportion of between vaccine recipient variability of S* that is
#            potentially protection relevant
# risk1      estimated probability that a vaccine recipient at-risk at tau experiences
#            the clinical endpoint by taumax.  It may be estimated differently for different studies.
# risk0      Same as risk1 for the placebo group
# This risk0 is used for both computepower and computepower.n

# PlatVElowest The percentage of vaccine recipients with the lowest value of VE (VElowest).

# VElowest A vector of the lowest possible value of vaccine efficacy.
#              Typical applications will range VElowest from 0 to VEoverall

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
# The first approach specifies spec, sens, FP0, FN2 which determine FP1 and FN1 from equations (8) and (9). spec, sens, 
# FP0, and FN2 must all be vectors of the same length.
#
# The second approach specifies sigma2obs and rho; this approach assumes the normal measurement error model (4) in the article.
# specifying spec=NULL, FP0=NULL, sens=NULL, FN2=NULL defaults to approach 2, which is used in the manuscript.

# Output: Power

# checks that sampling design input parameters are specified and valid, throwing an error if case-cohort
# sampling is chosen but p is unspecified or p is not a valid probability, or if case-control sampling 
# is chosen but controlCaseRatio is unspecified
checkSamplingDesign <- function(cohort, p, controlCaseRatio) {
  if(cohort==TRUE) {  #case-cohort 
    if (is.null(p)==TRUE) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is unspecified")
    } else if (p < 0 | p > 1) {
      stop("Case-cohort sampling was chosen and sampling probability, p, is not a valid probability")
    }
  } else if (is.null(controlCaseRatio)==TRUE) {  #case-control because cohort==FALSE
    stop("Case-control sampling was chosen and controlCaseRatio is unspecified")
  }
}

# checks that biomarker type and input parameters match, throwing an error if the biomarker type is binary
# but P0 + P2 != 1, or if the biomarker type is continuous but VElowest is NULL
checkBiomarkerType <- function(biomType, P0, P2, VElowest, PlatVElowest) {
  if((biomType=="binary") & (P0+P2 != 1)){
    stop("Binary biomarker was specified but P0 and P2 do not add up to 1")
  }
  if((biomType=="continuous") & (is.null(VElowest) | is.null(PlatVElowest))) {
    stop("Continuous biomarker was specified but VElowest and PlatVElowest are not specified")
  } 
}

# checks that sample size input parameters are valid, throwing an error if sample size inputs are vectors
# but rho is not scalar, or if sample size input vectors are of different lengths
checkSampleSizeParams <- function(sampleLengths, rho) {
  if(max(sampleLengths) > 1) {
    if(length(rho)>1) {
      stop("If multiple sample sizes are specified, input parameter rho must be scalar")
    } else if (max(sampleLengths) != min(sampleLengths)) {
      stop("Vector lengths differ for nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls")
    }
  }  
}

### Stop function
# Checks out that all values of RRlat2 are between 0 and 1 and that PlatVElowest meets bounds.
# If there are incompatible values of RRlat2, consider making Plat0 smaller,
# and/or making VElat0 larger, and then check again if all values of 
# RRlat2 are properly between 0 and 1.
checkProbabilityViolation <- function(VEoverall,RRlat2,PlatVElowest,VElowest, biomType) {
  if(biomType=="continuous") {
    if (min(VElowest)==0 & PlatVElowest > 1 - VEoverall) {
      stop("Input parameters PlatVElowest and VElowest violate probability constraints for normal marker calculations") 
    }
  } else if (any(RRlat2 < 0)) {
    cat(paste("RRlat2="),"\n")
    cat(paste(round(RRlat2,3)))
    cat("\n")
    stop("Input parameters violate probability constraints for trichotomous marker calculations. Consider making Plat0 smaller and/or RRlat0 smaller.")
  }
}

#############
# computeSensSpecFPFN is a function for mapping input parameters to sensitivity and specificity,
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
# This function also applies for a binary biomarker in which case only sensitivity and specificity
# are relevant (FP0, FP1, FN2, FN1 are not used in the calculations)
##################################################

# Original function that requires Plat2 to be a scalar
computeSensSpecFPFN <- function(sigma2obs,rho,Plat0,Plat2,P0,P2) {
  # sigma2tr = Var(X) = Var(Str) = rho*sigma2obs
  # If Plat0 + Plat2 = 1 then the method collapses to a binary biomarker,
  # and FP0, FP1, FN2, FN1 are irrelevant; this function simply returns 0's in that scenario
  
  # P2 may be a scalar or a vector, which should include one value equal to
  #      Plat2 and values straddling either side
  # P0 may be a scalar or a vector, which should include one value equal to
  #      Plat0 and values straddling either side
  
  sigma2e <- (1-rho)*sigma2obs
  sigma2tr <- rho*sigma2obs
  thetahiVE <- qnorm(1-Plat2)*sqrt(sigma2tr)
  thetaloVE <- qnorm(Plat0)*sqrt(sigma2tr)
  
  Plat1 <- 1 - Plat0 - Plat2
  m <- length(P2)
  ans <- list()
  
  for(i in 1:length(rho)){
    sens <- rep(1,m)
    spec <- rep(1,m)
    FP0 <- rep(0,m)
    FP1 <- rep(0,m)
    FN2 <- rep(0,m)
    FN1 <- rep(0,m)
    tauhisolution <- rep(0,m)
    taulosolution <- rep(0,m)
    if (rho[i] < 1) {  #*# if rho=1, then sens=1, spec=1, FP0=0, FP1=0, FN2=0, FN1=0
      # Stochastic integration
      X <- rnorm(20000,0,sqrt(sigma2tr[i]))
      S <- X + rnorm(20000,0,sqrt(sigma2e[i]))
      
      Phi <- sum(X>thetahiVE[i])/length(X)
      Plo <- sum(X<=thetaloVE[i])/length(X)
      Pmed <- 1 - Phi - Plo
      
      for (l in 1:m) {
        
        # Find the cut points tauhi and taulo by solving the following equations:
        #   0 = sensvec*Plat2 + FP1vec*Plat1 + FP0vec*Plat0 - P2  (f2 below; eqn 8 in manuscript)
        #   0 = specvec*Plat0 + FN1vec*Plat1 + FN2vec*Plat2 - P0  (f0 below; eqn 7 in manuscript)
        # where 
        #   sensvec <- (sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi
        #   specvec <- (sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo
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
        tauhisolution[l] <- tauhisol
        taulosolution[l] <- taulosol
        
        sens[l] <- sum(S>tauhisolution[l] & X > thetahiVE[i])/sum(X>thetahiVE[i])
        spec[l] <- sum(S<=taulosolution[l] & X <= thetaloVE[i])/sum(X<=thetaloVE[i])
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
                      taulosolution,tauhisolution,sens,spec,FP0,FP1,FN2,FN1)
  }  
  return(ans)
}

# check lengths of sens, spec, FP0, and FN2 vectors are equal
checkParamLengthsMatch <- function(sens, spec, FP0, FN2){
  lengths <- sapply(list(sens,spec,FP0,FN2), length)
  if(max(lengths) != min(lengths)){
    stop("Vector lengths differ for sens, spec, FP0, FN2")
  }
}

# Computes kernel of D(x, alphalat) (more information in Appendix B), which is the logit term in the zero-equation involving alphalat
computeKernel <- function(x, alpha, nus, risk1latnu, sigma2obs){
  rho <- 1
  piece1 <- exp(alpha*(1 - x/nus[1]))*(risk1latnu^(x/nus[1]))
  piece2 <- (1-risk1latnu)^(x/nus[1]) + piece1
  piece3 <- dnorm(x/(sqrt(rho*sigma2obs)))
  kernel <- (piece1/piece2)*piece3
  return(kernel)
}

# Equation whose root gives alphalat. Labeled U(alphalat) in Appendix B and based on eqn 13 in the manuscript
alphaLatEqn <- function(alpha, nus, risk1latnu, sigma2obs, VEoverall, PlatVElowest, risk0){
  logitterm <- integrate(computeKernel, lower=nus[1], upper=6, alpha=alpha, nus=nus, risk1latnu=risk1latnu, sigma2obs=sigma2obs)$value
  ans <- 1-VEoverall - (PlatVElowest*risk1latnu + logitterm)/risk0
  return(ans)
}

# Function for computing the infection probabilities of vaccinees
# risk1cont = risk_1^{lat}(x*) in the manuscript
risk1cont <- function(x,alphalat,betalat) {
  linpart <- alphalat + betalat*x
  ans <- exp(linpart)/(1+exp(linpart))
  return(ans) 
}

# Deal with rare crashes of rmultinom due to numerical problems where the
# program treats probability 0 as a small negative number:
adjustProb <- function(prob) {
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

# Given specifications for spec, sens, FP0, FN2, FP1, FN1, and a logical value indicating if the 
# biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2),
# where each subject is assigned a specific level
assignBiomarkerLevels <- function(specSens, binary, N0, N1, N2){
  spec <- specSens[1]
  sens <- specSens[2]
  FP0 <- specSens[3]
  FN2 <- specSens[4]
  FP1 <- specSens[5]
  FN1 <- specSens[6]
  if(binary==TRUE){
    Svalues <- cbind(rmultinom(N0,1,adjustProb(c(spec,1-FP0-spec,FP0))),
                     rmultinom(N2,1,adjustProb(c(FN2,1-FN2-sens,sens))))
  } else{
    Svalues <- cbind(rmultinom(N0,1,adjustProb(c(spec,1-FP0-spec,FP0))),
                     rmultinom(N1,1,adjustProb(c(FN1,1-FP1-FN1,FP1))),
                     rmultinom(N2,1,adjustProb(c(FN2,1-FN2-sens,sens))))
  }
  rownames(Svalues) <- c("S=0","S=1","S=2")
  Svalues <- ifelse(Svalues[1,]==1,0,ifelse(Svalues[2,]==1,1,2))
  return(Svalues)
}

# Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
biomSubset <- function(Y, N, nCasesPhase2, controlCaseRatio, p, cohort, pDropout){
  
  if (cohort==TRUE) {  # case-cohort sampling design
    
    # subset of subjects with biomarker measured is obtained by drawing a Bernoulli random sample from all at-risk observations 
    # to form the cohort, then augmenting the cohort with all cases
    
    R <- numeric(length(Y))
    if (!is.null(pDropout)){
      R <- ifelse(rbinom(N, 1, p)==1 & rbinom(N, 1, pDropout)==0, 1, R) # from N, draw Bernoulli sample with sampling probability p, and measure biomarker in those that do not drop out
    } else {
      R <- ifelse(rbinom(N, 1, p)==1, 1, R)
    }    
    R <- ifelse(Y==1, 1, R)  # augment all cases
    keepinds <- which(R==1)
    
  } else {  # case-control sampling design
    
    # Keep the S's in nCasesPhase2 of the cases (deleting the rest) and in controlCaseRatio*nCasesPhase2 controls
    casesinds <- which(Y==1)
    keepcasesinds <- sample(casesinds,nCasesPhase2,replace=FALSE)
    controlinds <- which(Y==0)
    keepcontrolinds <- sample(controlinds,controlCaseRatio*nCasesPhase2,replace=FALSE)
    keepinds <- sort(c(keepcasesinds,keepcontrolinds))
  }
  return(keepinds)
}

#' Sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement
#' error and treatment efficacy
#'
#' Performs sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement error and treatment efficacy [Gilbert, Janes, and Huang (2015).
#' ``Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials.'']
#'
#' @param nAtRiskTauCases  Number of subjects in the vaccine group at-risk at tau and with the clinical event (cases) by taumax (regardless of whether the biomarker is measured).
#' @param nAtRiskTauControls Number of subjects in the vaccine group at-risk at tau and without the clinical event (controls) by taumax (regardless of whether the biomarker is measured).
#' @param nAtRiskTauCasesPhase2 Number of subjects in the vaccine group at-risk at tau and with the clinical event (cases) by taumax and with the biomarker measured (i.e., in Phase 2).
#' @param controlCaseRatio Number of controls sampled per case in the vaccine arm (i.e. sampled in to Phase 2).
#' @param VEoverall Overall vaccine efficacy.
#' @param risk0 Estimated probability that a placebo recipient at-risk at tau experiences the clinical event by taumax.
#' @param VElat0 For a trichotomous biomarker (or binary as a special case), a vector of values corresponding to VE for the lowest (low) true biomarker subgroup.  Each value of \code{VElat0} corresponds to one unique effect size (RR_t).
#' @param VElat1 For a trichotomous biomarker, a vector of values corresponding to VE for the middle true biomarker subgroup.  For a binary biomarker, \code{VElat1} can be left unspecified.
#' @param VElowest For a continuous bioarker, a vector of values corresponding to the lowest possible value of VE.  Typical applications will range \code{VElowest} from 0 to 1 - \code{RRoverall}.
#' @param Plat0 For a trichotomous biomarker (or binary as a special case), probability that the (latent) true biomarker takes the lowest (low) value.
#' @param Plat2 For a trichotomous biomarker (or binary as a special case), probability that the (latent) true biomarker takes the highest (high) value.
#' @param P0 For a trichotomous biomarker (or binary as a special case), probability that the measured/observed biomarker takes the lowest (low) value.  If unspecified, this parameter is set to \code{Plat0}.
#' @param P2 For a trichotomous biomarker (or binary as a special case), probability that the measured/observed biomarker takes the highest (high) value.  If unspecified, this parameter is set to \code{Plat2}.
#' @param PlatVElowest For a continuous biomarker, the percentage of vaccine recipients with the lowest value of VE.
#' @param sens For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the specificity of the measured/observed biomarker. specifying \code{spec=NULL}, \code{FP0=NULL}, \code{sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param spec For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the sensitivity of the measured biomarker. specifying \code{spec=NULL}, \code{FP0=NULL}, \code{sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param FP0 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false positive rate (FP^1) of the measured/observed biomarker. specifying \code{spec=NULL}, \code{FP0=NULL}, \code{sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param FN2 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false negative rate (FN^2) of the measured/observed biomarker. specifying \code{spec=NULL}, \code{FP0=NULL}, \code{sens=NULL}, and \code{FN2=NULL} indicates that approach 2 is used.
#' @param M Number of simulated clinical trials.
#' @param alpha Two-sided type-I error rate for CoR hypothesis tests.
#' @param sigma2obs For a continuous biomarker, or for a trichotomous or binary biomarker simulated using `approach 2', the variance of the continuous measured/observed biomarker.
#' @param rho For a continuous biomarker, or for a trichotomous or binary biomarker simulated using 'approach 2', a vector of length 4 with values for the fraction of protection-relevant variability in the measured/observed continuous biomarker.  The first element of this vector should be 1, corresponding to the case of no measurement error.
#' @param biomType Type of biomarker that is used. The default is "continuous"; other choices are "trichotomous" and "binary".
#' @param cohort Sampling design to be used. Default is \code{FALSE}, specifying case-control sampling design. If \code{TRUE}, case-cohort sampling is used. 
#' @param p For case-cohort sampling design, probability that a subject will be in the cohort. 
#' @param pDropout For case-cohort sampling design, probability that a subject will drop out before time taumax. Defined as $P(\Delta = 0)$, where $\Delta$ is the indicator that $Y$ is observed. Default is \code{NULL}.
#' @param tpsMethod Character denoting method for fitting the logistic regression model. Choose from "PL" for pseudo-likelihood (default), "ML" for maximum likelihood, and "WL" for weighted likelihood. 
#' @param saveDir Character denoting the directory that the function output is to be saved in. Default is \code{NULL}. 
#' @param saveFile Character denoting the name of the file the function output will be saved in. Output will be saved as an .RData file. Default is \code{NULL}.
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
#' Approach 1 specifies \code{spec}, \code{sens}, \code{FP0}, and \code{FN2} which determine
#' \code{FP1} and \code{FN1} from equations (7) and (8).  Four values are required for each input parameter, to allow the evaluation of biomarkers with different levels of measurement error.
#'
#' Approach 2 for a trichotomous (or binary) biomarker specifies \code{sigma2obs} and \code{rho}; this approach assumes the normal measurement error model (4) in the manuscript.  
#' specifying \code{spec=NULL}, \code{FP0=NULL}, \code{sens=NULL}, and \code{FN2=NULL}
#' defaults to approach 2, which is what is used in illustrations in the manuscript.
#'
#' For a continuous biomarker, \code{VElowest}, \code{sigma2obs} and \code{rho} must be specified.  Setting \code{VElowest = NULL}
#' indicates that the biomarker is not continuous (it is trichotomous).
#'
#' This program implements a scenario with without-replacement-sampling (e.g., typically used in case-control and 2-phase sampling) and a scenario with case-cohort sampling
#' 
#' If \code{nAtRiskTauCases}, \code{nAtRiskTauControls}, and \code{nAtRiskTauCasesPhase2} are vectors, then \code{rho} must be scalar.
#'
#' @return Power- the fraction of simulated trials in which the null hypothesis H_0 (expression (14) of the manuscript for a trichotomous (or binary) biomarker and expression (16) for a continuous biomarker) is rejected.
#'
#' @examples
#' ## 'Global' parameters (independent of marker)
#' VEoverall <- 0.26   # VE in the at-risk-month-tau cohort
#' RRoverall <- 1 - VEoverall
#' nAtRiskTauCases <- 41
#' nAtRiskTauControls <- 7662
#' nAtRiskTauCasesPhase2 <- 41
#' risk1 <- nAtRiskTauCases/(nAtRiskTauCases + nAtRiskTauControls)
#' risk0 <- risk1/RRoverall # risk in placebo
#'
#' ## Parameters used for the trichotomous or binary biomarker calculations, Approach 1
#' ## For no measurement error scenario, set
#' ## spec[1]=1, FP0[1]=0, sens[1]=1, FN2[1]=0
#' ## Note the other elements of at least one of these four parameter vectors need
#' ## to be set to a different value to tell the program to use Approach 1
#' spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' ## Parameters used for continuous biomarker calculations, or
#' ## for trichotomous/binary biomarkers under Approach 2
#' ## Note these values are needed but are immaterial trichotomous/binary Approach 1
#' sigma2obs <- 1
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowest <- seq (0, 1-RRoverall,len=100)
#'
#' #################################################
#' ## use the function to perform power calculations
#' #################################################
#'
#' ## Binary biomarker, Approach 1 ##
#' spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' sens <- c(1, 0.9, 0.8, 0.7)
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
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowest=NULL,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rho, spec, FP0, sens, FN2)
#'
#'
#'## Trichotomous biomarker, Approach 1 ##
#' spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' sens <- c(1, 0.9, 0.8, 0.7)
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
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,
#' VElowest=NULL, controlCaseRatio, M, alpha=0.05, sigma2obs, rho, spec, FP0, sens, FN2)
#'
#'
#' ## Binary biomarker, Approach 2 ##
#' spec <- sens <- FP0 <- FN2 <- NULL
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
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls, risk0,
#' RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowest=NULL, controlCaseRatio,
#' M, alpha=0.05, sigma2obs, rho, spec, FP0, sens, FN2)
#'
#'
#' ## Trichotomous biomarker, Approach 2 ##
#' spec <- sens <- FP0 <- FN2 <- NULL
#'
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve

#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' ## Note these values are needed but are immaterial for trichotomous/binary markers
#' sigma2obs <- 1 #
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#'
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowest=NULL,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rho, spec, FP0, sens, FN2)
#'
#'
#' ## Continuous biomarker ##
#' sigma2obs <- 1
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowest <- seq (0, 1-RRoverall,len=100)
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
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest,VElowest,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rho)
#'
#' #######################
#' ## plotting the results
#' #######################
#'
#' ## trichotomous biomarker Approach 1
#' spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- rep(0,4)
#' sens <- c(1, 0.9, 0.8, 0.7)
#' FN2 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#' sigma2obs <- 1
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case


#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' M <- 1000
#' controlCaseRatio <- 5
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,
#' VElowest=NULL, controlCaseRatio, M, alpha=0.05, sigma2obs=NULL, rho=NULL, spec, FP0, sens, FN2)

#' ## plot power vs. CoR risk ratio in vaccine group (hi vs. lo)  (Figure 4)
#' ## file name = paste("powerstrinary",P2,P1,controlCaseRatio,".dat")
#' powstrin40105 <- matrix(scan('powerstrinary0.40.15.dat'),ncol=4,byrow=T)

#' VEoverall <- scan('VEoverallCoRpower.dat')
#' alpha <- scan('alpha.dat')
#' rho <- scan('rho.dat')
#' RRlat2 <- scan('RRlat2.dat')
#' RRlat0 <- scan('RRlat0.dat')
#' N <- scan('sampsizeALL.dat')
#' nCasesPhase2 <- scan('numbeventsPhase2.dat')

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
#' legend(x="topright",legend=c(paste("Power rho=",rho[1],sep=""),paste("Power rho=",rho[2],sep=""),paste("Power rho=",rho[3],sep=""),paste("Power rho=",rho[4],sep="")),
#' lty=c(1,2,3,4),col=c("blue","orange","green","black"),lwd=2)
#' mtext(paste("Power to Detect a Trichotomous CoR in Vaccine Recipients [2-sided alpha = ",alpha,"]"),outer=T,cex=1.3)
#' mtext(paste("Overall VE = ",VEoverall,"; Number controls  = ",round(controlCaseRatio*nCasesPhase2),"; Number cases = ",round(nCasesPhase2),"; Controls:cases = ",
#' controlCaseRatio,":1"),side=1,line=0.7,outer=T,cex=1.3)
#' mtext(paste("VElat_0 varies from ",VEoverall," to 0 as VElat_2 varies from ",VEoverall," to ",round(2*VEoverall,2)),side=1,line=3,outer=T,cex=1.3)
#' ## Note: the upper limit is VEoverall*(PlatloVE+PlathiVE)/PlathiVE, in the special case of this plot with
#' ## PlatloVE = PlathiVE, this simplifies to 2*VE_0.
#' dev.off()
#'
#' ## continuous biomarker
#' sigma2obs <- 1
#' rho <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowest <- seq (0, 1-RRoverall,len=100)
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
#' ans <- computepower(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest,VElowest,
#' controlCaseRatio, M, alpha=0.05, sigma2obs, rho)


#' ## Plot Power vs. CoR Relative Risk per SD Increase in X*  (Figure 2)
#' ## Power to Detect a Normally Distributed CoR in Vaccine Recipients [ 2-sided alpha  0.05]
#' VEoverall <- scan('VEoverallCoRpower.dat')
#' alpha <- scan('alpha.dat')
#' rho <- scan('rho.dat')
#' RRlat2 <- scan('RRlat2.dat')
#' RRlat0 <- scan('RRlat0.dat')
#' N <- scan('sampsizeALL.dat')
#' nCasesPhase2 <- scan('numbeventsPhase2.dat')
#' PlatVElowest <- scan('PlatVElowest.dat')
#' RRs <- scan('RRs.dat')
#' reverseRRs <- RRs[length(RRs):1]
#' ## The 5 refers to a controlCaseRatio of 5:1
#' powscont5 <- matrix(scan('powerscont5.dat'),ncol=4,byrow=T)
#' reversepowscont <- powscont5[nrow(powscont5):1,]
#'
#' sigma2tr <- 1
#' sigma2erho1 <- ((1-rho[1])/rho[1])*sigma2tr
#' sigma2erho2 <- ((1-rho[2])/rho[2])*sigma2tr
#' sigma2erho3 <- ((1-rho[3])/rho[3])*sigma2tr
#' sigma2erho4 <- ((1-rho[4])/rho[4])*sigma2tr
#' benchmarklab="V2"
#' benchmarkestRRrho1 <- 0.57 # From Haynes et al. 2012
#' benchmarkestRRrho2 <- benchmarkestRRrho1^(1/sqrt(rho[2]))
#' benchmarkestRRrho3 <- benchmarkestRRrho1^(1/sqrt(rho[3]))
#' benchmarkestRRrho4 <- benchmarkestRRrho1^(1/sqrt(rho[4]))
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
#' legend(x="topright",legend=c(paste("Power rho=",rho[1],sep=""),paste("Power rho=",rho[2],sep=""),paste("Power rho=",rho[3],sep=""),paste("Power rho=",rho[4],sep="")),
#'       lty=c(1,2,3,4),col=c("blue","orange","green","black"),lwd=2,cex=1.32)
#'
#' title(paste("Power to Detect a Normally Distributed CoR in Vaccine Recipients [2-sided alpha = ",alpha,"]"))
#' mtext(paste("Overall VE = ",VEoverall,"; Number controls  = ",round(nCasesPhase2*controlCaseRatio),"; Number cases = ",round(nCasesPhase2),"; Controls:cases = ",
#' controlCaseRatio,":1"),side=1,line=2,outer=T,cex=1.3)
#' dev.off()
#'
#'
#' @import survival
#' @import osDesign
#' @export
computePower <- function(nAtRiskTauCases, nAtRiskTauControls, nAtRiskTauCasesPhase2,
                         controlCaseRatio=NULL,
                         VEoverall, risk0, 
                         VElat0=seq(0, VEoverall, len=20), VElat1=rep(VEoverall, 20),
                         VElowest=NULL,
                         Plat0=0.2, Plat2=0.6,
                         P0=Plat0, P2=Plat2,
                         PlatVElowest=NULL, 
                         spec=NULL, FP0=NULL, sens=NULL, FN2=NULL,
                         M=100,
                         alpha=0.05,
                         sigma2obs=1, rho=1,
                         biomType=c("continuous", "trichotomous", "binary"),
                         cohort=FALSE, p=NULL, pDropout=NULL,
                         tpsMethod=c("PL", "ML","WL"),
                         saveDir=NULL, saveFile=NULL) {
  
  # sigma2tr is the variance of the true biomarker X
  # rho must be a vector with 4 values and is for the continuous and binary biomarker correlates correlations
  # RRlat0 is the span of true relative risks (vaccine vs. placebo) in the latent lower protected subgroup
  # RRlat1 is the span of true relative risks (vaccine vs. placebo) in the latent medium protected subgroup
  # The power calculations should always include rho=1 in the first element of the rho vector,
  # as the best case scenario, and the plotting functions assume this.
  
  # VElowest is used for a continuous biomarker- a vector of fixed value of VE(x) for
  # the subgroup of subjects with lowest X^* values, where this subgroup has prevalence PlatVElowest
  
  
  tpsMethod <- match.arg(tpsMethod, choices = c("PL","ML","WL"))
  biomType <- match.arg(biomType, choices = c("continuous", "trichotomous", "binary"))
  
  # check sampling design input parameters are specified and valid
  checkSamplingDesign(cohort, p, controlCaseRatio)
  # check biomarker type and input parameters match
  checkBiomarkerType(biomType, P0, P2, VElowest, PlatVElowest)
  
  # check sample size parameters are valid
  sampleLengths <- sapply(list(nAtRiskTauCases, nAtRiskTauCasesPhase2, nAtRiskTauControls), length)
  checkSampleSizeParams(sampleLengths, rho)
  
  
  nCases <- nAtRiskTauCases
  nCasesPhase2 <- nAtRiskTauCasesPhase2
  nControls <- nAtRiskTauControls
  # Overall denominator: number observed to be at risk when the immune response is measured (N in manuscript):
  N <- nCases + nControls
  
  # Compute VElat2:
  RRoverall <- 1 - VEoverall
  RRlat0 <- 1 - VElat0
  RRlat1 <- 1 - VElat1
  VElat2 <- (VEoverall*(Plat0+Plat2) - Plat0*VElat0)/Plat2  # This formula assumes VElat1 = VEoverall
  RRlat2 <-round(1-VElat2, 10)   # rounded to avoid problems when 0 is treated as a small negative number
  Plat1 <- 1 - Plat2 - Plat0
  P1 <- 1 - P0 - P2
  
  # check all values of RRlat2 are between 0 and 1 and that PlatVElowest meets bounds
  checkProbabilityViolation(VEoverall,RRlat2,PlatVElowest,VElowest, biomType)
  
  sigma2e <- (1-rho)*sigma2obs
  sigma2tr <- rho*sigma2obs
  
  #################################################
  # Computations for a trinary biomarker
  if(biomType=="trichotomous" | biomType=="binary") {
    
    # Compute VElat2:
    
    Approach2 <- (all(is.null(spec), is.null(sens), is.null(FP0), is.null(FN2)))
    
    if (Approach2) {  # Default choice
      
      # Compute sens, spec, FP0, FP1, FN2, FN1
      
      ans <- computeSensSpecFPFN(sigma2obs,rho,Plat0,Plat2,P0,P2)
      sens <- unlist(lapply(ans, function(x) x[[1,10]])) 
      spec <- unlist(lapply(ans, function(x) x[[1,11]]))
      FP0 <- unlist(lapply(ans, function(x) x[[1,12]])) 
      FP1 <- unlist(lapply(ans, function(x) x[[1,13]])) 
      FN2 <- unlist(lapply(ans, function(x) x[[1,14]]))
      FN1 <- unlist(lapply(ans, function(x) x[[1,15]])) 
      
      
      # dataframe of rho, sens, spec, etc.
      # used to create Table 1: mapping of sigma2obs and rho to the sens, spec, etc. parameters
      table1 <- as.data.frame(round(cbind(rho, Plat0, P0, Plat2, P2, sens, spec, FP0, FN2, FP1, FN1),3))
      
    }
    
    # Approach 1 in the manuscript:
    if (!Approach2) { #*# use given sens, spec, FP0, and FN2 params
      
      # check lengths of sens, spec, FP0, and FN2 vectors are equal
      checkParamLengthsMatch(sens,spec,FP0,FN2)
      
      # Apply formula (7) in the manuscript
      FN1 <- (P0 - spec*Plat0 - FN2*Plat2)/Plat1   #*#P0, Plat0, Plat2 given params
      # Apply formula (8) in the manuscript
      FP1 <- (P2 - sens*Plat2 - FP0*Plat0)/Plat1
      
      # Check if an error in the ranges of values due to an out of
      # bounds input parameter
      
      if (any(FN1 < 0 | FN1 > 1 | FP1 < 0 | FP1 > 1)){
        stop("Approach 1 was used and one of the parameters sens, spec, FP0, FN2 is out of range")
      }
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
    probX2_cond_S2 <- sens*Plat2/P2
    # use outer product to get matrix with nrow=length(rho), ncol=length(RRlat0)
    risk1_2 <- (probX0_cond_S2 %o% RRlat0 + probX1_cond_S2 %o% RRlat1 + probX2_cond_S2 %o% RRlat2 )*risk0  
    probX0_cond_S0 <- spec*Plat0/P0
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
    
    # initialize power calculation matrix
    if (max(sampleLengths) > 1) {
      powerstrinary <- matrix(0, nrow=length(N), ncol=ncol(esvect))
      rownames(powerstrinary) <- paste0(rep("N"), seq(1,nrow(powerstrinary)))
    } else {
      powerstrinary <- matrix(0, nrow=nrow(esvect), ncol=ncol(esvect))
      rownames(powerstrinary) <- paste0(rep("rho"), seq(1,nrow(powerstrinary)))
    }
    
    
    ###################################################
    # Power calculations repeated for M simulations
    
    for (i in 1:M) {  # M simulations
      
      for (j in 1:ncol(esvect)) {  # for each value of RRlat0
        
        # Fix the number of cases and controls, putting the cases first and controls
        # second for each subgroup:
        
        rrlat0 <- risk1lat_0[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        rrlat1 <- risk1lat_1[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        rrlat2 <- risk1lat_2[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        denominat <- Plat0*rrlat0 + Plat1*rrlat1 + Plat2*rrlat2
        P0case <- (Plat0*rrlat0)/denominat  #*# success probabilities for trinomial random variable
        P1case <- (Plat1*rrlat1)/denominat
        P2case <- 1 - P0case - P1case
        
        # loops through the different sample sizes
        for(k in 1:length(N)) {
          # Deal with rare crashes of rmultinom due to numerical problems where the
          # program treats probability 0 as a small negative number:
          inds <- rmultinom(nCases[k],1,adjustProb(c(P0case,P1case,P2case)))
          
          # Number of cases in the lo, med, hi latent groups
          nCases0 <- length(inds[1,][inds[1,]==1])
          nCases2 <- length(inds[3,][inds[3,]==1])
          nCases1 <- nCases[k] - nCases0 - nCases2
          N0 <- round(Plat0*N[k])
          N2 <- round(Plat2*N[k])
          N1 <- N[k] - N0 - N2
          
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
          
          Y <- c(rep(1,nCases0),rep(0,N0-nCases0),rep(1,nCases1),rep(0,N1-nCases1),rep(1,nCases2),rep(0, N[k] - N0 - N1 - nCases2))
          
          # Simulate the trinary surrogate with 0,1,2 = lo,med,hi
          # Formulas (12) and (13) in the manuscript:
          
          # Given specifications for spec, FP0, sens, and FN2 and a logical value indicating if the 
          # biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2),
          # where each subject is assigned a specific level
          specSens <- cbind(spec,sens,FP0,FN2,FP1,FN1)
          
          if (biomType=="binary") { # binary case only
            S <- t(apply(specSens, 1, function(x) assignBiomarkerLevels(x, binary=TRUE, N0, N1, N2))) # each row is a set of sens, spec, etc. parameters
          } else { # trichotomous
            S <- t(apply(specSens, 1, function(x) assignBiomarkerLevels(x, binary=FALSE, N0, N1, N2))) # each row is a set of sens, spec, etc. parameters
          }
          
          # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
          keepinds <- biomSubset(Y, N[k], nCasesPhase2[k], controlCaseRatio, p, cohort, pDropout)
          
          # Those with biomarker data:
          Ycc <- Y[keepinds]
          Scc <- t(apply(S,1,function(x) x[keepinds])) #nrow=length(rho)
          
          ##############################################################
          # Now analyze with osDesign
          # (first check if there are 'zeros', in which case Fisher's exact test for the lo vs. hi categories is used. Otherwise,
          # osDesign logistic regression is used as an ordered score test
          
          if(max(sampleLengths)>1) {
            lodim <- dim(table(Ycc,Scc))[2]<2 #*# check there are at least two biomarker categories (columns)
            zerosflag <-  lodim
            if (dim(table(Ycc,Scc))[2]==3) { #*# check if any categories have zero entries
              zerosflag <- table(Ycc,Scc)[1,1]==0 | table(Ycc,Scc)[1,2]==0 | table(Ycc,Scc)[1,3]==0 | table(Ycc,Scc)[2,1]==0 | table(Ycc,Scc)[2,2]==0 | table(Ycc,Scc)[2,3]==0 
            }
            
            if (zerosflag) {
              if (lodim) { pval <- 1}
              if (!lodim) { #*# there are zeros, so Fisher's exact test is used
                pval <- fisher.test(table(Ycc,Scc)[,c(1,dim(table(Ycc,Scc))[2])])$p.value 
              }
              if (pval <= alpha & length(Ycc[Scc==2&Ycc==1])/length(Scc[Scc==2]) < length(Ycc[Scc==0&Ycc==1])/length(Scc[Scc==0])) {
                powerstrinary[k,j] <- powerstrinary[k,j] + 1
              }
            }
            
            if (!zerosflag) {
              fit <- tps(Ycc~Scc[1,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
              pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
              if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[k,j] <- powerstrinary[k,j] + 1}
            }
          } else {
            for(l in 1:nrow(Scc)){
              lodim <- dim(table(Ycc,Scc[l,]))[2]<2 #*# check there are at least two biomarker categories (columns)
              zerosflag <-  lodim
              if (dim(table(Ycc,Scc[l,]))[2]==3) { #*# check if any categories have zero entries
                zerosflag <- table(Ycc,Scc[l,])[1,1]==0 | table(Ycc,Scc[l,])[1,2]==0 | table(Ycc,Scc[l,])[1,3]==0 | table(Ycc,Scc[l,])[2,1]==0 | table(Ycc,Scc[l,])[2,2]==0 | table(Ycc,Scc[l,])[2,3]==0 
              }
              
              if (zerosflag) {
                if (lodim) { pval <- 1}
                if (!lodim) { #*# there are zeros, so Fisher's exact test is used
                  pval <- fisher.test(table(Ycc,Scc[l,])[,c(1,dim(table(Ycc,Scc[l,]))[2])])$p.value 
                }
                if (pval <= alpha & length(Ycc[Scc[l,]==2&Ycc==1])/length(Scc[l,][Scc[l,]==2]) < length(Ycc[Scc[l,]==0&Ycc==1])/length(Scc[l,][Scc[l,]==0])) {
                  powerstrinary[l,j] <- powerstrinary[l,j] + 1
                }
              }
              
              if (!zerosflag) {
                fit <- tps(Ycc~Scc[l,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
                pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
                if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[l,j] <- powerstrinary[l,j] + 1}
              }
            }
          }
          
        }  
      }  
    }  
    # power calculations
    power <- powerstrinary/M
    
    # write out alpha intercept as logit(Y=1|s=0) for trinary/binary case
    alphaLat <- c(logit(risk1_0))
    # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/binary case
    betaLat <- c(logit(risk1_2)-logit(risk1_0))
    # CoR effect sizes
    RRt <- risk1_2/risk1_0
    
    pwr <- list("power"=power, "RRt"=RRt, "risk1_2"=risk1_2, "risk1_0"=risk1_0, "VElat2"=VElat2, "VElat0"=VElat0, "Plat2"=Plat2, "Plat0"=Plat0, 
                "P2"=P2, "P0"=P0, "alphaLat"=alphaLat, "betaLat"=betaLat, "sens"=sens, "spec"=spec, "FP0"=FP0, "FN2"=FN2)
    
    write(RRlat2,file="RRlat2.dat",ncolumns=1,append=FALSE)
    write(RRlat0,file="RRlat0.dat",ncolumns=1,append=FALSE)
    write(power,file=paste("powerstrinary",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=nrow(powerstrinary),append=FALSE)
    write(P2,file="P2.dat")
    # Print out the CoR effect sizes
    write(risk1_0,file=paste("vaccineriskslo",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=length(rho),append=FALSE)
    write(risk1_2,file=paste("vaccineriskshi",P2,P0,controlCaseRatio,".dat",sep=""),ncolumns=length(rho),append=FALSE)
    # write out alpha intercept as logit(Y=1|s=0) for trinary/binary case
    write(c(t(logit(risk1_0))), file="trinaryalpha.dat",ncolumns=1,append=FALSE)
    # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/binary case
    write(c(t(logit(risk1_2)-logit(risk1_0))), file="trinarybeta.dat",ncolumns=1,append=FALSE)
    
  } else if (biomType=="continuous") {  
    
    #################################################
    # Computations for a continuous biomarker
    # Define the truebetas indexed by the user-specified vector VElowest
    
    o <- length(VElowest)
    nus <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
    truebetas <- rep(NA,o)
    alphalatvect <- rep(NA,o)
    
    for (l in 1:o) {
      
      # find solutions alphalat and betalat by solving eqn (4) in Appendix B
      risk1latnu <- (1-VElowest[l])*risk0
      
      alphalatvect[l] <- uniroot(alphaLatEqn, lower=-10, upper=10, nus=nus, risk1latnu=risk1latnu, sigma2obs=sigma2obs, VEoverall=VEoverall, PlatVElowest=PlatVElowest, risk0=risk0)$root
      
      # Second solve for betalat:
      D <- risk1latnu
      truebetas[l] <- (log(D/(1-D)) - alphalatvect[l])/nus[1]
    }
    
    # initialize power calculation matrix
    if(max(sampleLengths)>1) {
      powerscont <- matrix(0, nrow=length(N), ncol=length(VElowest))
      rownames(powerscont) <- paste0(rep("N"), seq(1,nrow(powerscont)))
    } else {
      powerscont <- matrix(0, nrow=length(rho), ncol=length(VElowest))
      rownames(powerscont) <- paste0(rep("rho"), seq(1,nrow(powerscont)))
    }
    
    
    ###################################################
    # Power calculations repeated for M simulations
    
    for (i in 1:M) {
      
      # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
      # using the function risk1cont() above
      
      # If multiple sample size specified, loops through the different sample sizes;
      # else, loops through the different RRlat0's
      
      for (j in 1:o) {
        
        beta <- truebetas[j]
        alphalat <- alphalatvect[j]
        
        # These simulations condition on n (i.e., number of infections in vaccine arm) and
        # also on the number of controls fixed at controlCaseRatio*n
        # This matches what was done for a binary correlate
        
        if(max(sampleLengths)>1){
          
          for(k in 1:length(N)){
            
            # Arbitrarily put the cases first and controls second
            # The numbers of cases and controls are fixed, e.g., a typical retrospective design
            Y <- c(rep(1,nCases[j]),rep(0,N[j]-nCases[j]))
            
            # Compute the denominator of the density of X|Y=1 when Y|X follows a log. regr model with the truncated part
            # associated with VElowest and X is normal
            # with mean zero and standard deviation sqrt(sigma2tr)
            #
            # rho[1]:
            
            f <- function(x) {
              ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
              return(ans)
            }
            denomdensityXcases <- integrate(f,lower=nus,upper=5)$value
            denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowest[j])*risk0
            
            numerdensXcases <- function(x) {
              num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
              num[x <= nus] <- PlatVElowest*(1-VElowest[j])*risk0
              return(num)
            }
            numerdensXcontrols <- function(x) {
              num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr))
              num[x <= nus] <- PlatVElowest*(1-(1-VElowest[j])*risk0)
              return(num)
            }
            
            Xpoints <- seq(-3.5,3.5,len=25000)
            probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
            probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)
            
            Xcases <-    sample(Xpoints,size=nCases[k],prob=probscases,replace=TRUE)
            Xcontrols <- sample(Xpoints,size=N[k]-nCases[k],prob=probscontrols,replace=TRUE)
            X <- c(Xcases,Xcontrols)
            
            # Create the immune response variables for the different degrees of measurement error
            error <- rnorm(N[k],mean=0,sd=sqrt(sigma2e))
            S <- X + error
            
            # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
            keepinds <- biomSubset(Y, N[k], nCasesPhase2, controlCaseRatio, p, cohort, pDropout)
            
            # Those with biomarker data:
            Ycc <- Y[keepinds]
            Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rho)
            
            fit <- tps(Ycc~Scc[1,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
            pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}
            
          }  
          
        } else {
          # Arbitrarily put the cases first and controls second
          # The numbers of cases and controls are fixed, e.g., a typical retrospective design
          Y <- c(rep(1,nCases),rep(0,N-nCases))
          
          # Compute the denominator of the density of X|Y=1 when Y|X follows a log. regr model with the truncated part
          # associated with VElowest and X is normal
          # with mean zero and standard deviation sqrt(sigma2tr)
          #
          # rho[1]:
          
          for(k in 1:length(nus)){
            f <- function(x) {
              ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
              return(ans)
            }
            denomdensityXcases <- integrate(f,lower=nus[k],upper=5)$value
            denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowest[j])*risk0
            
            numerdensXcases <- function(x) {
              num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
              num[x <= nus[k]] <- PlatVElowest*(1-VElowest[j])*risk0
              return(num)
            }
            numerdensXcontrols <- function(x) {
              num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr[k]))
              num[x <= nus[k]] <- PlatVElowest*(1-(1-VElowest[j])*risk0)
              return(num)
            }
            
            Xpoints <- seq(-3.5,3.5,len=25000)
            probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
            probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)
            
            Xcases <-    sample(Xpoints,size=nCases,prob=probscases,replace=TRUE)
            Xcontrols <- sample(Xpoints,size=N-nCases,prob=probscontrols,replace=TRUE)
            X <- c(Xcases,Xcontrols)
            
            # Create the immune response variables for the different degrees of measurement error
            error <- rnorm(N,mean=0,sd=sqrt(sigma2e[k]))
            S <- X + error
            
            # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
            keepinds <- biomSubset(Y, N, nCasesPhase2, controlCaseRatio, p, cohort, pDropout)
            
            # Those with biomarker data:
            Ycc <- Y[keepinds]
            Scc <- S[keepinds] 
            
            fit <- tps(Ycc~Scc,nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
            pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}
          }
        }
      }
    } 
    # power calculations
    power <- powerscont/M
    
    # RRc the relative risks that are the effect sizes RR_c that need to be on the x-axis of powerplots
    RRc <- exp(truebetas)
    
    pwr <- list("power"=power, "RRc"=RRc, "betaLat"=truebetas, "PlatVElowest"=PlatVElowest, "VElowest"=VElowest, "sigma2obs"=sigma2obs)
    
    # RRs the relative risks that are the effect sizes RR_c that
    # need to be on the x-axis of powerplots
    write(exp(truebetas),file="RRs.dat",ncolumns=1,append=FALSE)
    write(power,file=paste("powerscont",controlCaseRatio,".dat",sep=""),ncolumns=nrow(powerscont),append=FALSE)
    write(PlatVElowest,file="PlatVElowest.dat")
    write(VElowest,file="VElowest.dat",ncolumns=1,append=FALSE)
    write(truebetas,file="truebetas.dat")
  }
  
  # VEoverall <- 1-RRoverall
  # ans <- c(ans, list(N), list(nCases), list(nCasesPhase2), VEoverall, alpha, list(rho), controlCaseRatio)
  pwr$N <- N
  pwr$nCases <- nCases
  pwr$nCasesPhase2 <- nCasesPhase2
  pwr$VEoverall <- 1-RRoverall
  pwr$alpha <- alpha
  pwr$rho <- rho
  pwr$controlCaseRatio <- controlCaseRatio
  write(N,file="sampsizeALL.dat")
  write(nCases,file="numbeventsALL.dat")
  write(nCasesPhase2,file="numbeventsPhase2.dat")
  write(1-RRoverall,file="VEoverallCoRpower.dat")
  write(alpha,file="alpha.dat")
  write(rho,file="rho.dat",ncolumns=1,append=FALSE)
  write(controlCaseRatio,file="controlCaseRatio.dat")
  if(!is.null(saveDir) & !is.null(saveFile)) {
    save(pwr, file=paste0(file.path(saveDir, saveFile),".RData"))
  }
  
  return(pwr)
  
}
