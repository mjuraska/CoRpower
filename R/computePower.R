# Main program computePower() for correlate of risk (CoR) power calculations and
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
# P1 = P(S=1) = P(S=med) (==0 for dichotomous S)
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

# Let the CoR effect size in vaccinees be RR_t = risk1(2)/risk1(0),
# This is the x-axis of the power curves

# The program considers CoR analysis in the vaccine group only, where the goal is to test
# H0: risk1(2)=risk1(1)=risk1(0) vs. H1: risk1(2) < risk1(0), i.e.
# H0: RR_t = 1 vs. H1: RR_t < 1
# based on the case-control sample from vaccine recipients

# sensitivity, specificity, false negatives and false positives
# sens \equiv P(S(1)=2|X=2) (==1 for S with no measurement error)
# spec \equiv P(S(1)=0|X=0) (==1 for S with no measurement error)
# FP^0 \equiv P(S(1)=2|X=0) (==0 for S with no measurement error)
# FN^2 \equiv P(S(1)=0|X=2) (==0 for S with no measurement error)
# FP^1 \equiv P(S(1)=2|X=1) (==0 for S with no measurement error)
# FN^1 \equiv P(S(1)=0|X=1) (==0 for S with no measurement error)

########################################################################
# End of unique notations for trinary immune response S = 0,1,2
# The remaining notations below are the same for trinary (Approach 2) and continous S
# except that PlatVElowest and VElowest only apply to continuous S
#######################################################################

# The biomarker is measured at time tau and subjects are followed through time taumax

# RRoverall = 1 - overall vaccine efficacy
# annincinfectionplac = annual HIV infection incidence in placebo group
# annincdropout = annual dropout rate assumed the same in both groups

# nCases     All cases in the vaccine group at-risk at tau and a case by taumax
#                             (regardless of whether the biomarker is measured)
# nControls  All controls in the vaccine group at-risk at tau and not diseased at the end of follow-up taumax
#                             (regardless of whether the biomarker is measured)
# nCasesWithS  As above and also have the biomarker measured (i.e., in Phase 2)

# sigma2obs  observed variance of the continuous marker S*
# rho        vector of rho, the proportion of between vaccine recipient variability of S* that is
#            potentially protection relevant
# risk1      estimated probability that a vaccine recipient at-risk at tau experiences
#            the clinical endpoint by taumax.  It may be estimated differently for different studies.
# risk0      Same as risk1 for the placebo group

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

# Output: Power


# Internal Functions

checkSamplingDesign <- function(cohort, p, controlCaseRatio) {
  # Checks that sampling design input parameters are specified and valid.
  #
  # Args:
  #   cohort: If TRUE, case-cohort sampling is used; if not, case-control sampling is used. Default is FALSE.
  #   p: Probability a subject will be in the cohort
  #   controlCaseRatio: Ratio of controls to cases
  #
  # Returns:
  #   Error if case-cohort sampling is chosen but p is unspecified or p is not a valid probability,
  #   or if case-control sampling is chosen but controlCaseRatio is unspecified

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

checkBiomarkerType <- function(biomType, P0, P2, VElowest, PlatVElowest) {
  # Checks that biomarker type and input parameters match.
  #
  # Args:
  #   biomType: Character string specifying type of biomarker used. Default is "continuous";
  #             other choices are "trichotomous" and "dichotomous"
  #   P0: Probability of low biomarker response
  #   P2: Probability of high biomarker response
  #   VElowest: Lowest VE level for true biomarker X*
  #   PlatVElowest: Prevalence of VElowest
  #
  # Returns:
  #   Error if the biomarker type is dichotomous but P0 + P2 != 1,
  #   or if the biomarker type is continuous but VElowest is NULL

  if((biomType=="dichotomous") & (P0+P2 != 1)){
    stop("dichotomous biomarker was specified but P0 and P2 do not add up to 1")
  }
  if((biomType=="continuous") & (is.null(VElowest) | is.null(PlatVElowest))) {
    stop("Continuous biomarker was specified but VElowest and PlatVElowest are not specified")
  }
}

checkSampleSizeParams <- function(sampleLengths, rho) {
  # Checks that sample size input parameters are valid.
  #
  # Args:
  #   sampleLengths: Numeric vector denoting lengths of the sample size input parameters
  #   rho: protecion-relevant fraction of the variance of S or S*
  #
  # Returns:
  #   Error if sample size inputs are vectors but rho is not scalar,
  #   or if sample size input vectors are of different lengths

  if(max(sampleLengths) > 1) {
    if(length(rho)>1) {
      stop("If multiple sample sizes are specified, input parameter rho must be scalar")
    } else if (max(sampleLengths) != min(sampleLengths)) {
      stop("Vector lengths differ for nCases, nCasesWithS, nControls")
    }
  }
}

checkVElat1violation <- function(VElat0, VElat1, biomType) {
  # Checks that 'VElat0' and 'VElat1' are valid and that specifications for 'biomType' and 'VElat1' match.
  #
  # Args:
  #   VElat0: a numeric vector specifying a grid of treatment efficacy levels in the latent lower protected
  #           subgroup for a dichotomous or trichotomous biomarker
  #   VElat1: a numeric vector specifying a grid of treatment efficacy levels in the latent medium protected
  #           subgroup for a trichotomous biomarker (NULL by default for a dichotomous biomarker)
  #   biomType: a character string specifying the biomarker type
  #
  # Returns:
  #   Error if there is a discrepancy between 'VElat1' and 'biomType' or if an element in 'VElat1' is less
  #   than its corresponding element in 'VElat0'
  if(biomType == "dichotomous" & !(is.null(VElat1))) {
    stop("There is a discrepancy between 'VElat1' and 'biomType'")
  } else if(any(VElat1 < VElat0) == "TRUE") {
    stop("An element in 'VElat1' is less than its corresponding element in 'VElat0'")
  }
}

checkProbabilityViolation <- function(VEoverall,RRlat2,PlatVElowest,VElowest, biomType) {
  # Checks that all values of RRlat2 are between 0 and 1 and that PlatVElowest and VElowest meet bounds.
  #
  # Args:
  #   VEoverall: overall vaccine efficacy
  #   RRlat2: grid of relative rsisk (vaccine/placebo) among higher protected latent subgroup
  #   PlatVElowest: Prevalence of VElowest
  #   VElowest: Lowest VE level for true biomarker X*
  #   biomType: Character string specifying type of biomarker used. Default is "continuous";
  #             other choices are "trichotomous" and "dichotomous"
  #
  # Returns:
  #   Error if there are incompatible values of RRlat2, or if values of PlatVElowest and
  #   VElowest violate probabiliy constraints for normal marker caclualtions

  if(biomType=="continuous") {
    if (min(VElowest)==0 & PlatVElowest > 1 - VEoverall) {
      stop("Input parameters PlatVElowest and VElowest violate probability constraints for normal biomarker calculations")
    }
  } else if (any(RRlat2 < 0)) {
    cat(paste("RRlat2="),"\n")
    cat(paste(round(RRlat2,3)))
    cat("\n")
    stop("Input parameters violate probability constraints for trichotomous biomarker calculations.
         Consider making Plat0 smaller and/or VElat0 larger.")
  }
}

computeSensSpecFPFN <- function(sigma2obs,rho,Plat0,Plat2,P0,P2) {
  # For trichotomous biomarker specified using Approach 2, maps input parameters rho and sigma2obs
  # to sensitivity, specificity, FP0, FP1, FN2, and FN1 values (defined above).
  # Can also be used for a dichotomous biomarker (Plat0 + Plat2 = 1), in which case only sensitivity and specificity
  # are relevant (FP0, FP1, FN2, FN1 are not used in the calculations).
  #
  # Args:
  #   sigma2obs: Variance of observed biomarker S*
  #   rho: Protection-relevant fraction of the variance of S*
  #   Plat0: Prevalence of lower protected latent subgroup.
  #   Plat2: Prevalence of higher protected latent subgroup.
  #   P0: Probability of low biomarker response.
  #   P2: Probability of high biomarker response.
  #
  # Returns:
  #   Matrix with each row corresponding to a value of rho and with the following columns:
  #   thetaloVE*, thetahiVE*, Plat0, Plat1, Plat2, P0, P2, taulosolution**, tauhisolution**,
  #   sens, spec, FP0, FP1, FN2, and FN1.
  #   If dichotomous biomarker, returns 0's for FP0, FP1, FN2, and FN1, which are irrelavant.
  #
  #   *Let X* denote the true latent subgroups. X* = 2 if X* > thetahiVE, X* = 0 if X <= thetaloVE, and
  #   X* = 1 if X* is in between thetaloVE and thetahiVE, for fixed thetaloVE and thetahiVE that are solved for.
  #   **Let S* denote a trichotomous biomarker. S* = 2 if S* > tauhisolution and S* = 0 if S* <= taulosolution,
  #   and S* = 1 if S* is in between taulosolution and tauhisolution, for fixed taulosolution and tauhisolution
  #   that are solved for.

  # Based upon classical measurement error model S* = X* + e  where e ~ N(0,sigma2e), X* ~ N(0,sigma2tr)
  # sigma2obs = Var(S*) = sigma2tr + sigma2e
  # rho = 1 - sigma2e/sigma2obs = sigma2tr/sigma2ob
  # sigma2tr = Var(X*) = rho*sigma2obs
  sigma2e <- (1-rho)*sigma2obs
  sigma2tr <- rho*sigma2obs
  thetahiVE <- qnorm(1-Plat2)*sqrt(sigma2tr)
  thetaloVE <- qnorm(Plat0)*sqrt(sigma2tr)

  Plat1 <- 1 - Plat0 - Plat2
  ans <- list()

  for(i in 1:length(rho)){
    sens <- 1
    spec <- 1
    FP0 <- 0
    FP1 <- 0
    FN2 <- 0
    FN1 <- 0
    tauhisolution <- 0
    taulosolution <- 0
    if (rho[i] < 1) {  # if rho=1, then sens=1, spec=1, FP0=0, FP1=0, FN2=0, FN1=0
      # Stochastic integration
      X <- rnorm(20000,0,sqrt(sigma2tr[i]))
      S <- X + rnorm(20000,0,sqrt(sigma2e[i]))

      Phi <- sum(X>thetahiVE[i])/length(X)
      Plo <- sum(X<=thetaloVE[i])/length(X)
      Pmed <- 1 - Phi - Plo

      # Find the cut points tauhi and taulo by solving the following equations:
      #   0 = sensvec*Plat2 + FP1vec*Plat1 + FP0vec*Plat0 - P2  (f2 below; eqn 8 in manuscript)
      #   0 = specvec*Plat0 + FN1vec*Plat1 + FN2vec*Plat2 - P0  (f0 below; eqn 7 in manuscript)
      # where
      #   sensvec <- (sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi
      #   specvec <- (sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo
      # if dichotomous biomarker,
      #   FP1vec <- 0
      #   FP0vec <- 0
      #   FN2vec <- 0
      #   FN1vec <- 0
      # if trichotomous biomarker,
      #   FP1vec <- (sum(S>tauhi & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
      #   FP0vec <- (sum(S>tauhi & X <= thetaloVE[i])/length(S))/Plo
      #   FN2vec <- (sum(S<=taulo & X > thetahiVE[i])/length(S))/Phi
      #   FN1vec <- (sum(S<=taulo & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed

      if (Pmed==0){  # dichotomous
        f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi)*Plat2 - P2
        f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo)*Plat0 - P0
      } else {  # trichotomous
        f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE[i])/length(S))/Phi)*Plat2 +
          ((sum(S>tauhi & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed)*Plat1 +
          ((sum(S>tauhi & X <= thetaloVE[i])/length(S))/Plo)*Plat0 - P2
        f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE[i])/length(S))/Plo)*Plat0 +
          ((sum(S<=taulo & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed)*Plat1 +
          ((sum(S<=taulo & X > thetahiVE[i])/length(S))/Phi)*Plat2 - P0
      }
      tauhisolution <- uniroot(f2, interval=c(-2.5,2.5))$root
      taulosolution <- uniroot(f0, interval=c(-2.5,2.5))$root

      sens <- sum(S>tauhisolution & X > thetahiVE[i])/sum(X>thetahiVE[i])
      spec <- sum(S<=taulosolution & X <= thetaloVE[i])/sum(X<=thetaloVE[i])
      if (Pmed==0) {  # if dichotomous biomarker, 0's for FP1, FP0, FN2, FN1
        FP1 <- 0
        FP0 <- 0
        FN2 <- 0
        FN1 <- 0 }
      if (Pmed > 0) {
        FP1 <- (sum(S>tauhisolution & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        FP0 <- (sum(S>tauhisolution & X <= thetaloVE[i])/length(S))/Plo
        FN2 <- (sum(S<=taulosolution & X > thetahiVE[i])/length(S))/Phi
        FN1 <- (sum(S<=taulosolution & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
      }
    }
    ans[[i]] <- cbind(thetaloVE[i],thetahiVE[i],Plat0,Plat1,Plat2,P0,P2,
                      taulosolution,tauhisolution,sens,spec,FP0,FP1,FN2,FN1)
  }
  return(ans)
}


checkParamLengthsMatch <- function(sens, spec, FP0, FN2){
  # Checks that the lengths of sens, spec, FP0, and FN2 are equal.
  #
  # Args:
  #   sens: Numeric scalar or vector corresponding to P(S=2|X=2)
  #   spec: Numeric scalar or vector corresponding to P(S=0|X=0)
  #   FP0: Numeric scalar or vector corresponding to P(S=2|X=0)
  #   FN2: Numeric scalar or vector corresponding to P(S=0|X=2)
  #
  # Returns:
  #   Error if the vector lengths differ for sens, spec, FP0, and FN2

  lengths <- sapply(list(sens,spec,FP0,FN2), length)
  if(max(lengths) != min(lengths)){
    stop("Vector lengths differ for sens, spec, FP0, FN2")
  }
}

computeKernel <- function(x, alphaLat, nu, risk1latnu, sigma2obs){
  # Computes the kernel of D(x, alphalat), which is the logit term in the
  # zero-equation involving alphaLat (more information in Appendix B).
  #
  # Args:
  #   x: Variable to be integrated over
  #   alphaLat: variable in logistic regression model: logit(risk1lat(x*)) = alphaLat + betaLat x*
  #   nu: Value that denotes the cut-off for lowest X* values. Fraction PlatVElowest of subjects have x* value below this
  #       and therefore have VE = VElowest. nu = sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  #   risk1latnu: Vaccine-group endpoint risk for true biomarker x* <= nu
  #   sigma2obs: Variance of observed biomarker S*
  #
  # Returns:
  #   Kernel of logit term in zero-equation involving alphaLat

  rho <- 1
  piece1 <- exp(alphaLat*(1 - x/nu[1]))*(risk1latnu^(x/nu[1]))
  piece2 <- (1-risk1latnu)^(x/nu[1]) + piece1
  piece3 <- dnorm(x/(sqrt(rho*sigma2obs)))
  kernel <- (piece1/piece2)*piece3
  return(kernel)
}

alphaLatEqn <- function(alphaLat, nu, risk1latnu, sigma2obs, VEoverall, PlatVElowest, risk0){
  # Computes alphalat by solving for the root of the zero-equation labeled U(alphalat) in Appendix B and
  # based on eqn 13 in the manuscript.
  #
  # Args:
  #   alphaLat: Variable to be used in uniroot() equation
  #   nu: Value that denotes the cut-off for lowest X* values. Fraction PlatVElowest of subjects have x* value below this
  #       and therefore have VE = VElowest. nu = sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  #   risk1latnu: Vaccine-group endpoint risk for true biomarker x* <= nu
  #   sigma2obs: Variance of observed biomarker S*
  #   VEoverall: Overal vaccine efficacy
  #   PlatVElowest: Prevalence of VElowest
  #   risk0: Placebo-group endpoint risk
  #
  # Returns:
  #   alphaLat as the solution to zero-equation involving alphaLat

  logitterm <- integrate(computeKernel, lower=nu[1], upper=6, alphaLat=alphaLat, nu=nu,
                         risk1latnu=risk1latnu, sigma2obs=sigma2obs)$value
  ans <- 1-VEoverall - (PlatVElowest*risk1latnu + logitterm)/risk0
  return(ans)
}

risk1cont <- function(x,alphalat,betalat) {
  # Computes vaccine-group endpoint risk for true biomarker x* > nu as logit^{-1}(alphalat+betalat*x) = risk1cont,
  # where risk1cont = risk_1^{lat}(x*) for x* > nu in the manuscript.
  #
  # Args:
  #   x: Variable
  #   alphalat: Numeric value that is solved for in other functions. Parameter in logistic regression model.
  #   betalat: Numeric value that is solved for in other functions. Parameter in logistic regression model.
  #
  # Returns:
  #   risk1cont: vaccine-group endpoint risk for true biomarker x* > nu

  linpart <- alphalat + betalat*x
  ans <- exp(linpart)/(1+exp(linpart))
  return(ans)
}

adjustProb <- function(prob) {
  # Adjusts probabilities to prevent rare crashes of rmultinom due to numerical problems where the program treats
  # probability 0 as a small negative number.
  #
  # Args:
  #   prob: Success probabilities for trinomial/binomial random variable during simulation process.
  #         Consists of linear combination of sens, spec, FP0, FP1, FN2, and FN1.
  #
  # Returns:
  #   Adjusted probability where ties are broken and problematic values are corrected

  # Break ties:
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

assignBiomarkerLevels <- function(specSens, dichotomous, N0, N1, N2){
  # Assigns an observed biomarker level (S=0, 1, or 2) to each subject.
  #
  # Args:
  #   specSens: Numeric matrix where each row is a set of spec, sens, FP0, FN2, FP1, and FN1 values
  #   dichotomous: If TRUE, indicates biomarker is dichotomous; if FALSE, indicates biomarker is trichotomous
  #   N0: Number of subjects at risk at tau in the lower protected latent subgroup, excluding dropouts
  #   N1: Number of subjects at risk at tau in the medium protected latent subgroup, excluding dropouts
  #   N2: Number of subjects at risk at tau in the higher protected latent subgroup, excluding dropouts
  #
  # Returns:
  #   Vector composed of biomarker levels (0,1,2), where each subject is assigned a specific level

  spec <- specSens[1]
  sens <- specSens[2]
  FP0 <- specSens[3]
  FN2 <- specSens[4]
  FP1 <- specSens[5]
  FN1 <- specSens[6]
  if(dichotomous==TRUE){
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

biomSubset <- function(Y, Ncomplete, nCasesWithS, controlCaseRatio, p, cohort){
  # Selects subset of subjects that have biomarker S or S* measured (R_i=1) according to
  # a case-cohort or case-control sampling design.
  #
  # Args:
  #   Y: Numeric vector indicating cases and controls (1 vs. 0), with length = Ncomplete
  #   Ncomplete: Total number of subjects at risk at tau, excluding dropouts
  #   nCasesWithS: Number of observed cases between tau and taumax with measured S or S*
  #   controlCaseRatio: Ratio of controls to cases in case-control sampling design
  #   p: Probability a subject will be in the cohort
  #   cohort: If TRUE, indicates case-cohort sampling design; if FALSE, indicates case-control sampling design
  #
  # Returns:
  #   Indices of the subjects selected to have biomarker measured

  if (cohort==TRUE) {  # case-cohort sampling design

    # Subset of subjects with biomarker measured is obtained by drawing a Bernoulli random sample from all at-risk observations
    # to form the cohort, then augmenting the cohort with all cases
    R <- numeric(length(Y))
    R <- ifelse(rbinom(Ncomplete, 1, p)==1, 1, R) # from (Ncomplete=nCases+nControls), draw Bernoulli sample with sampling probability p
    R <- ifelse(Y==1, 1, R)  # augment all cases
    keepinds <- which(R==1)

  } else {  # case-control sampling design

    # Keep the S's in nCasesWithS of the cases (deleting the rest) and in controlCaseRatio*nCasesWithS controls
    casesinds <- which(Y==1)
    keepcasesinds <- sample(casesinds,nCasesWithS,replace=FALSE)
    controlinds <- which(Y==0)
    keepcontrolinds <- sample(controlinds,controlCaseRatio*nCasesWithS,replace=FALSE)
    keepinds <- sort(c(keepcasesinds,keepcontrolinds))
  }
  return(keepinds)
}


#' Power Calculations for Assessing Intermediate Biomarkers as Correlates of Risk in the Active Treatment Group in Clinical Efficacy Trials, Accounting for Biomarker's Measurement Error and Treatment Efficacy
#'
#' Performs a power calculation for assessing a univariate dichotomous, trichotomous, or continuous intermediate biomarker response as a correlate of risk
#' in the active treatment group in a clinical efficacy trial, accounting for the biomarker's measurement error and treatment efficacy. The statistical methods are described in [Gilbert, Janes, and Huang (2016).
#' ``Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials.'']
#'
#' @param nCases an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group (a numeric vector of multiple counts/scenarios is allowed)
#' @param nControls an integer value specifying the number of controls observed (or projected) to complete follow-up through \eqn{\tau_{max}} endpoint-free in the active treatment group (a numeric vector of multiple counts/scenarios is allowed)
#' @param nCasesWithS an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response (a numeric vector of multiple counts/scenarios is allowed)
#' @param controlCaseRatio an integer value specifying the number of controls sampled per case for biomarker measurement in the without replacement case-control sampling design
#' @param VEoverall a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#' @param risk0 a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#' @param VElat0 a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lower protected subgroup for a dichotomous or trichotomous biomarker. Each value of \code{VElat0} corresponds to one unique effect size (\eqn{RR_t}). It typically ranges from \code{VEoverall} (\eqn{H_0}) to 0 (maximal \eqn{H_1} not allowing harm by treatment).
#' @param VElat1 a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent medium protected subgroup for a trichotomous biomarker (\code{NULL} by default for a dichotomous biomarker)
#' @param VElowest a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup for a continuous biomarker. It typically ranges from \code{VEoverall} (\eqn{H_0}) to 0 (maximal \eqn{H_1} not allowing harm by treatment).
#' @param Plat0 a numeric value specifying the prevalence of the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#' @param Plat2 a numeric value specifying the prevalence of the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#' @param P0 a numeric value specifying the probability of low biomarker response for a dichotomous or trichotomous biomarker. If unspecified, it is set to \code{Plat0}.
#' @param P2 a numeric value specifying the probability of high biomarker response for a dichotomous or trichotomous biomarker. If unspecified, it is set to \code{Plat2}.
#' @param PlatVElowest a numeric value specifying the prevalence of the latent lowest-efficacy subgroup for a continuous biomarker
#' @param sens a numeric vector specifying the sensitivity, i.e., the probability of high biomarker response conditional on membership in the higher protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'.
#' @param spec a numeric vector specifying the specificity, i.e., the probability of low biomarker response conditional on membership in the lower protected subgroup, of a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'.
#' @param FP0 a numeric vector specifying the false positive rate, i.e., the probability of high biomarker response conditional on membership in the lower protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'.
#' @param FN2 a numeric vector specifying the false negative rate, i.e., the probability of low biomarker response conditional on membership in the higher protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'.
#' @param M an integer value specifying the number of simulated clinical trials
#' @param alpha a numeric value specifying the two-sided Wald test type-I error rate
#' @param sigma2obs a numeric value specifying the variance of the observed continuous biomarker or of the dichotomous or trichotomous biomarker simulated using 'approach 2'
#' @param rho a numeric vector specifying distinct protection-relevant fractions of \code{sigma2obs}
#' @param biomType a character string specifying the biomarker type. Default is \code{continuous}; other choices are \code{dichotomous} and \code{trichotomous}.
#' @param cohort a logical value for whether a case-cohort Bernoulli sampling design is to be used. If \code{FALSE} (default), the case-control without replacement sampling is used.
#' @param p a numeric value specifying the probability of sampling into the subcohort in the case-cohort design
#' @param tpsMethod a character string specifying the estimation method in the inverse probability weighted logistic regression model fit by the \code{tps} function in the \code{osDesign} package. The options are \code{PL} for pseudo-likelihood (default), \code{ML} for maximum likelihood, and \code{WL} for weighted likelihood.
#' @param saveDir a character string specifying the path for a directory in which the output is to be saved. If \code{NULL} (default), the output is returned only.
#' @param saveFile a character string specifying the name of the \code{.RData} file storing the output. If \code{NULL} (default), the output is returned only.
#'
#' @details
#' If \code{nCases}, \code{nControls}, and \code{nCasesWithS} are vectors (of the same length), then \code{rho} must be a scalar.
#'
#' To save output in an \code{.RData} file, both \code{saveDir} and \code{saveFile} must be specified.
#'
#' Parameters independent of biomarker type and sampling design: \code{nCases}, \code{nControls}, \code{nCasesWithS}, \code{VEoverall}, \code{risk0},
#' \code{M}, \code{alpha}, \code{tpsMethod}, \code{saveDir}, \code{saveFile}.
#'
#' Parameters for trichotomous (or dichotomous) biomarker: \code{VElat0}, \code{VElat1}, \code{Plat0}, \code{Plat2}, \code{P0},
#' \code{P2}, \code{biomType="trichotomous"} (or \code{"dichotomous"})
#'   \itemize{
#'     \item Parameters for Approach 1: \code{sens}, \code{spec}, \code{FP0}, \code{FN2}
#'     \item Parameters for Approach 2: \code{sigma2obs}, \code{rho}
#'   }
#'
#' Parameters for continuous biomarker: \code{VElowest}, \code{PlatVElowest}, \code{sigma2obs}, \code{rho}, \code{biomType = "continuous"}
#'
#' Parameters for a case-control without replacement sampling design: \code{controlCaseRatio}
#'
#' Parameters for a case-cohort Bernoulli sampling design: \code{cohort=TRUE}, \code{p}
#'
#' @return If \code{saveFile} and \code{saveDir} are both specified, the output list (named \code{pwr}) is saved as an \code{.RData} file; otherwise it is returned only.
#' For a dichotomous or trichotomous biomarker, the output list has the following components:
#' \itemize{
#'   \item \code{power}: a matrix of fractions of simulated trials in which the null hypothesis \eqn{H_0} is rejected. Rows represent calculations for different values of \code{rho}, \code{sens}, or \code{nCases}, depending on which is a vector. Columns represent calculations for the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{RRt}: a matrix of correlate-of-risk relative-risk effect sizes. Rows represent different values of \code{rho}, \code{sens}, or \code{nCases}, depending on which is a vector. Columns represent the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{risk1_2}: a matrix of conditional endpoint risks given a high biomarker response in the active treatment group. Rows represent different values of \code{rho}, \code{sens}, or \code{nCases}, depending on which is a vector. Columns represent the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{risk1_0}: a matrix of conditional endpoint risks given a low biomarker response in the active treatment group. Rows represent different values of \code{rho}, \code{sens}, or \code{nCases}, depending on which is a vector. Columns represent the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{VElat2}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{VElat0}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{Plat2}: a numeric value specifying the prevalence of the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{Plat0}: a numeric value specifying the prevalence of the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{P2}: a numeric value specifying the probability of high biomarker response for a dichotomous or trichotomous biomarker
#'   \item \code{P0}: a numeric value specifying the probability of low biomarker response for a dichotomous or trichotomous biomarker
#'   \item \code{alphaLat}: a numeric vector of the log odds of the clinical endpoint in the subgroup of active treatment recipients with the latent \eqn{x^{\ast}=0} (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{betaLat}: a numeric vector of the log odds ratio of the clinical endpoint comparing two subgroups of active treatment recipients differing in the latent \eqn{x*} by 1 (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{sens}: a numeric vector of sensitivities (i.e., the probability of high biomarker response conditional on membership in the higher protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{spec}: a numeric vector of specificities (i.e., the probability of low biomarker response conditional on membership in the lower protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{FP0}: a numeric vector of false positive rates (i.e., the probability of high biomarker response conditional on membership in the lower protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{FN2}: a numeric vector of false negative rates (i.e., the probability of low biomarker response conditional on membership in the higher protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{Ncomplete}: an integer value specifying \code{nCases} + \code{nControls}, i.e., the number, observed or projected, of active treatment recipients at risk at \eqn{\tau} with an observed endpoint or a completed follow-up through \eqn{\tau_{max}}
#'   \item \code{nCases}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group
#'   \item \code{nCasesWithS}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response
#'   \item \code{controlCaseRatio}: an integer specifying the number of controls sampled per case for biomarker measurement in the without replacement case-control sampling design
#'   \item \code{VEoverall}: a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{risk0}: a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{alpha}: a numeric value specifying the two-sided Wald test type-I error rate
#'   \item \code{rho}: a numeric vector specifying distinct protection-relevant fractions of the variance of the observed biomarker
#' }
#'
#' For a continuous biomarker, a list with the following components:
#' \itemize{
#'   \item \code{power}: a matrix of fractions of simulated trials in which the null hypothesis \eqn{H_0} is rejected. Rows represent calculations for different values of \code{rho} or \code{nCases}, depending on which is a vector. Columns represent calculations for the grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup, specified by \code{VElowest}.
#'   \item \code{RRc}: a numeric vector of correlate-or-risk relative-risk effect sizes as a function of the grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup, specified by \code{VElowest}
#'   \item \code{betaLat}: a numeric vector specifying the log odds ratio of the clinical endpoint comparing two subgroups of active treatment recipients differing in the latent \eqn{x*} by 1 (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{alphaLat}: a numeric vector specifying the the log odds of the clinical endpoint in the subgroup of active treatment recipients with the latent \eqn{x^{\ast}=0} (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{PlatVElowest}: a numeric value specifying the prevalence of the latent lowest-efficacy subgroup for a continuous biomarker
#'   \item \code{VElowest}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup for a continuous biomarker
#'   \item \code{sigma2obs}: a numeric value specifying the variance of the observed continuous biomarker or of the dichotomous or trichotomous biomarker simulated using 'approach 2'
#'   \item \code{Ncomplete}: an integer value specifying \code{nCases} + \code{nControls}, i.e., the number, observed or projected, of active treatment recipients at risk at \eqn{\tau} with an observed endpoint or a completed follow-up through \eqn{\tau_{max}}
#'   \item \code{nCases}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group
#'   \item \code{nCasesWithS}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response
#'   \item \code{VEoverall}: a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{alpha}: a numeric value specifying the two-sided Wald test type-I error rate
#'   \item \code{rho}: a numeric vector specifying distinct protection-relevant fractions of the variance of the observed biomarker
#'   \item \code{controlCaseRatio}: an integer value specifying the number of controls sampled per case for biomarker measurement in the without replacement case-control sampling design
#'   \item \code{risk0}: a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#' }
#'
#' @examples
#'
#'## Trichotomous biomarker, Approach 1, varying sens and spec ##
#'## Specify sens, spec, FP0, FN2
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' VElat0 <- seq(0, VEoverall, len=20)  # 20 data points for the power curve
#' VElat1 <- rep(VEoverall, 20)
#' Plat0 <- 0.2
#' Plat2 <- 0.6
#' P0 <- Plat0  # different values of P0 can be set
#' P2 <- Plat2  # different values of P2 can be set
#' sens <- spec <- c(1, 0.9, 0.8, 0.7)
#' FP0 <- FN2 <- rep(0, 4)
#' M <- 5
#' alpha <- 0.05
#' biomType <- "trichotomous"
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              VElat0=VElat0, VElat1=VElat1, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'              M=M, alpha=alpha, spec=spec, FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' \dontrun{
#' ## Trichotomous biomarker, Approach 2, varying rho ##
#' ## Specify rho and sigma2obs
#'
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' VElat0 <- seq(0, VEoverall, len=20)
#' VElat1 <- rep(VEoverall, 20)
#' Plat0 <- 0.2
#' Plat2 <- 0.6
#' P0 <- Plat0
#' P2 <- Plat2
#' M <- 5
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.9, 0.7, 0.5)
#' biomType <- "trichotomous"
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              VElat0=VElat0, VElat1=VElat1, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#'
#' ## dichotomous biomarker, Approach 2, varying rho ##
#' ## Plat0 + Plat2 = 1
#'
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' VElat0 <- seq(0, VEoverall, len=20)  # 20 data points for the power curve
#' VElat1 <- rep(0, 20)  # will not be used by function
#' Plat0 <- 0.2
#' Plat2 <- 1 - Plat0
#' P0 <- Plat0
#' P2 <- Plat2
#' M <- 5
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.9, 0.7, 0.5)
#' biomType <- "dichotomous"
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              VElat0=VElat0, VElat1=VElat1, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#'
#' ## Continuous biomarker, varying rho ##
#'
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' PlatVElowest <- 0.2
#' VElowest <- seq(0, VEoverall, len=20)
#' M <- 5
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.9, 0.7, 0.5)
#' biomType <- "continuous"
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              PlatVElowest=PlatVElowest, VElowest=VElowest, M=M, alpha=alpha,
#'              sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#'
#' ## Continuous biomarker, case-cohort sampling design, varying p ##
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' PlatVElowest <- 0.2
#' VElowest <- seq(0, VEoverall, len=20)
#' M <- 5
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- 0.9
#' biomType <- "continuous"
#' cohort <- TRUE
#' p <- 0.01
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              VEoverall=VEoverall, risk0=risk0, PlatVElowest=PlatVElowest, VElowest=VElowest,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType, cohort=cohort,
#'              p=p)
#' p <- 0.02
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              VEoverall=VEoverall, risk0=risk0, PlatVElowest=PlatVElowest, VElowest=VElowest,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType, cohort=cohort, p=p)
#' p <- 0.03
#' computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'              VEoverall=VEoverall, risk0=risk0, PlatVElowest=PlatVElowest, VElowest=VElowest,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType, cohort=cohort, p=p)
#'
#'
#' ## Continuous biomarker, saving output, varying sample sizes ##
#'
#' nCases <- 32
#' nControls <- 1000
#' nCasesWithS <- 32
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' PlatVElowest <- 0.2
#' VElowest <- seq(0, VEoverall, len=20)
#' M <- 5
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.9, 0.7, 0.5)
#' biomType <- "continuous"
#' saveDir <- "~/myDir"
#' saveFile <- "MyFile"
#' computePower(nCases=nCases, nCasesWithS=nCasesWithS, nControls=nControls,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              PlatVElowest=PlatVElowest, VElowest=VElowest, M=M, alpha=alpha,
#'              sigma2obs=sigma2obs, rho=rho, biomType=biomType, saveDir=saveDir, saveFile=saveFile)
#' }
#'
#' @seealso \code{\link{computeN}}
#'
#' @import survival
#' @import osDesign
#' @importFrom stats dnorm fisher.test integrate pexp pnorm qnorm rbinom rmultinom rnorm uniroot
#'
#' @export
computePower <- function(nCases, nControls, nCasesWithS,
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
                         biomType=c("continuous", "trichotomous", "dichotomous"),
                         cohort=FALSE, p=NULL,
                         tpsMethod=c("PL", "ML","WL"),
                         saveDir=NULL, saveFile=NULL) {


  tpsMethod <- match.arg(tpsMethod, choices = c("PL","ML","WL"))
  biomType <- match.arg(biomType, choices = c("continuous", "trichotomous", "dichotomous"))

  # check sampling design input parameters are specified and valid
  checkSamplingDesign(cohort, p, controlCaseRatio)
  # check biomarker type and input parameters match
  checkBiomarkerType(biomType, P0, P2, VElowest, PlatVElowest)

  # check sample size parameters are valid
  sampleLengths <- sapply(list(nCases, nCasesWithS, nControls), length)
  checkSampleSizeParams(sampleLengths, rho)

  # Overall number observed to be at risk when the immune response is measured and that do not drop out (smaller than N):
  Ncomplete <- nCases + nControls

  # Compute VElat2:
  RRoverall <- 1 - VEoverall
  RRlat0 <- 1 - VElat0
  RRlat1 <- 1 - VElat1
  Plat1 <- 1 - Plat0 - Plat2
  P1 <- 1 - P0 - P2
  VElat2 <- (VEoverall - (Plat0*VElat0 + Plat1*VElat1))/Plat2  # This formula assumes VElat1 = VEoverall
  RRlat2 <-round(1-VElat2, 10)   # rounded to avoid problems when 0 is treated as a small negative number

  # check VElat0 and VElat1 are valid and specifications for biomType and VElat1 match
  checkVElat1violation(VElat0, VElat1, biomType)

  # check all values of RRlat2 are between 0 and 1 and that PlatVElowest meets bounds
  checkProbabilityViolation(VEoverall,RRlat2,PlatVElowest,VElowest, biomType)


  sigma2e <- (1-rho)*sigma2obs
  sigma2tr <- rho*sigma2obs  # variance of true biomarker X


  #################################################
  # Computations for a trinary biomarker
  if(biomType=="trichotomous" | biomType=="dichotomous") {

    Approach2 <- (all(is.null(spec), is.null(sens), is.null(FP0), is.null(FN2)))

    # Approach 2 in the manuscript (default choice):
    if (Approach2) {

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
    if (!Approach2) {  # use given sens, spec, FP0, and FN2 params

      # check lengths of sens, spec, FP0, and FN2 vectors are equal
      checkParamLengthsMatch(sens,spec,FP0,FN2)

      if (biomType=="dichotomous") {  # if dichotomous biomarker, FN1 and FP1 are irrelevant
        FN1 <- 0
        FP1 <- 0
      } else {
        # Apply formula (7) in the manuscript
        FN1 <- (P0 - spec*Plat0 - FN2*Plat2)/Plat1   #P0, Plat0, Plat2 given params
        # Apply formula (8) in the manuscript
        FP1 <- (P2 - sens*Plat2 - FP0*Plat0)/Plat1
      }

      # Check if there is an error in the ranges of values due to an out of bounds input parameter
      if (any(FN1 < 0 | FN1 > 1 | FP1 < 0 | FP1 > 1)){
        stop("Approach 1 was used and one of the parameters sens, spec, FP0, FN2 is out of range")
      }
    }

    # dichotomous biomarker special case (to remove small values of P1)
    if (biomType=="dichotomous") {
      P1 <- 0
      P2 <- 1 - P0
    }

    # Compute the marginal risks:
    # Made it to the end of follow-up HIV negative
    risk1 <- RRoverall*risk0

    # Observed risks P(Y(1)=1|S(1)=0, 1, or 2)
    # for diff values of rho; using Bayes' rule
    probX0_cond_S2 <- FP0*Plat0/P2
    probX1_cond_S2 <- FP1*Plat1/P2
    probX2_cond_S2 <- sens*Plat2/P2
    # use outer product to get matrix with nrow=length(rho), ncol=length(RRlat0)
    risk1_2 <- (probX0_cond_S2 %o% RRlat0 + probX1_cond_S2 %o% RRlat1 + probX2_cond_S2 %o% RRlat2 )*risk0
    probX0_cond_S0 <- spec*Plat0/P0
    probX1_cond_S0 <- FN1*Plat1/P0
    probX2_cond_S0 <- FN2*Plat2/P0
    risk1_0 <- (probX0_cond_S0 %o% RRlat0 + probX1_cond_S0 %o% RRlat1 + probX2_cond_S0 %o% RRlat2)*risk0
    risk1_1 <- (risk1 - risk1_0*P0 - risk1_2*P2)/P1  # Note: For the dichotomous biomarker special case, the risk1medx are NA

    esvect <- risk1_2/risk1_0  # matrix with nrow=length(rho) and ncol=length(RRlat0)

    # Vaccine risks within the latent subgroups (independent of rho of course)
    risk1lat_2 <- RRlat2*risk0
    risk1lat_1 <- RRlat1*risk0
    risk1lat_0 <- RRlat0*risk0

    # initialize power calculation matrix
    if (max(sampleLengths) > 1) {
      powerstrinary <- matrix(0, nrow=length(Ncomplete), ncol=ncol(esvect))
      rownames(powerstrinary) <- paste0(rep("N"), seq(1,nrow(powerstrinary)))
    } else {
      powerstrinary <- matrix(0, nrow=nrow(esvect), ncol=ncol(esvect))
      if(length(rho)>1) {
        rownames(powerstrinary) <- paste0(rep("rho"), seq(1,nrow(powerstrinary)))
      } else if(length(sens)>1) {
        rownames(powerstrinary) <- paste0(rep("sens/spec"), seq(1,nrow(powerstrinary)))
      }
    }


    ###################################################
    # Power calculations repeated for M simulations
    for (i in 1:M) {  # M simulations; Step 8 in the manuscript

      for (j in 1:ncol(esvect)) {  # for each value of RRlat0

        # Determine success probabilities for trinomial random variable:
        # P(X=0|Y=1, Y^tau=0, Z=1), P(X=1|Y=1, Y^tau=0, Z=1), P(X=2|Y=1, Y^tau=0, Z=1),
        # using Bayes rule to express them in terms of Platx and risk1(x), and risk1
        rrlat0 <- risk1lat_0[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])  # risk1(0)
        rrlat1 <- risk1lat_1[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        rrlat2 <- risk1lat_2[j]/(risk1lat_0[j]+risk1lat_1[j]+risk1lat_2[j])
        denominat <- Plat0*rrlat0 + Plat1*rrlat1 + Plat2*rrlat2  # risk1
        P0case <- (Plat0*rrlat0)/denominat  # success probabilities for trinomial random variable
        P1case <- (Plat1*rrlat1)/denominat
        P2case <- 1 - P0case - P1case

        # loops through the different sample sizes
        for(k in 1:length(Ncomplete)) {

          # Draw from trinomial random variable with success probabilities defined above.
          # adjustProb() function deals with rare crashes of rmultinom due to numerical problems
          # where the program treats probability 0 as a small negative number
          inds <- rmultinom(nCases[k],1,adjustProb(c(P0case,P1case,P2case)))

          # Number of cases in the 0, 1, 2 latent groups
          nCases0 <- length(inds[1,][inds[1,]==1])
          nCases2 <- length(inds[3,][inds[3,]==1])
          nCases1 <- nCases[k] - nCases0 - nCases2
          N0 <- round(Plat0*Ncomplete[k])
          N2 <- round(Plat2*Ncomplete[k])
          N1 <- Ncomplete[k] - N0 - N2

          # Address rounding that could make N1 negative in the dichotomous marker case
          # Keep Ncomplete fixed at a constant
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

          # Fix the number of cases and controls, putting the cases first and controls second for each subgroup:
          Y <- c(rep(1,nCases0),rep(0,N0-nCases0),rep(1,nCases1),rep(0,N1-nCases1),rep(1,nCases2),rep(0, Ncomplete[k] - N0 - N1 - nCases2))

          # Simulate the trinary surrogate
          # Formulas (12) and (13) in the manuscript:

          # Given specifications for spec, FP0, sens, and FN2 and a logical value indicating if the biomarker is dichotomous
          # or not, the function assignBiomarkerLevels() returns a vector composed of biomarker levels (S=0,1,2),
          # where each subject is assigned a specific level
          specSens <- cbind(spec,sens,FP0,FN2,FP1,FN1)
          if (biomType=="dichotomous") { # dichotomous
            S <- t(apply(specSens, 1, function(x) assignBiomarkerLevels(x, dichotomous=TRUE, N0, N1, N2))) # each row is set of sens, etc. parameters
          } else { # trichotomous
            S <- t(apply(specSens, 1, function(x) assignBiomarkerLevels(x, dichotomous=FALSE, N0, N1, N2))) # each row is set of sens, etc. parameters
          }

          # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
          keepinds <- biomSubset(Y, Ncomplete[k], nCasesWithS[k], controlCaseRatio, p, cohort)

          # Those with biomarker data:
          Ycc <- Y[keepinds]
          Scc <- t(apply(S,1,function(x) x[keepinds])) #nrow=length(rho)

          ##############################################################
          # Now analyze with osDesign
          # (first check if there are 'zeros', in which case Fisher's exact test for the lo vs. hi categories is used.
          # Otherwise, osDesign logistic regression is used as an ordered score test

          if(max(sampleLengths)>1) {
            lodim <- dim(table(Ycc,Scc))[2]<2 # check there are at least two biomarker categories (columns)
            zerosflag <-  lodim
            if (dim(table(Ycc,Scc))[2]==3) { # check if any categories have zero entries
              zerosflag <- table(Ycc,Scc)[1,1]==0 | table(Ycc,Scc)[1,2]==0 | table(Ycc,Scc)[1,3]==0 | table(Ycc,Scc)[2,1]==0 |
                           table(Ycc,Scc)[2,2]==0 | table(Ycc,Scc)[2,3]==0
            }

            if (zerosflag) {
              if (lodim) { pval <- 1}
              if (!lodim) { # there are zeros, so Fisher's exact test is used
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
              lodim <- dim(table(Ycc,Scc[l,]))[2]<2 # check there are at least two biomarker categories (columns)
              zerosflag <-  lodim
              if (dim(table(Ycc,Scc[l,]))[2]==3) { # check if any categories have zero entries
                zerosflag <- table(Ycc,Scc[l,])[1,1]==0 | table(Ycc,Scc[l,])[1,2]==0 | table(Ycc,Scc[l,])[1,3]==0 |
                             table(Ycc,Scc[l,])[2,1]==0 | table(Ycc,Scc[l,])[2,2]==0 | table(Ycc,Scc[l,])[2,3]==0
              }

              if (zerosflag) {
                if (lodim) { pval <- 1}
                if (!lodim) { # there are zeros, so Fisher's exact test is used
                  pval <- fisher.test(table(Ycc,Scc[l,])[,c(1,dim(table(Ycc,Scc[l,]))[2])])$p.value
                }
                if (pval <= alpha & length(Ycc[Scc[l,]==2&Ycc==1])/length(Scc[l,][Scc[l,]==2]) <
                    length(Ycc[Scc[l,]==0&Ycc==1])/length(Scc[l,][Scc[l,]==0])) {
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

    # write out alpha intercept as logit(Y=1|s=0) for trinary/dichotomous case
    alphaLat <- c(logit(risk1_0))
    # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/dichotomous case
    betaLat <- c(logit(risk1_2)-logit(risk1_0))
    # CoR effect sizes
    RRt <- risk1_2/risk1_0

    # output list for trichotomous/dichotomous biomarker
    pwr <- list("power"=power, "RRt"=RRt, "risk1_2"=risk1_2, "risk1_0"=risk1_0, "VElat2"=VElat2, "VElat0"=VElat0,
                "Plat2"=Plat2, "Plat0"=Plat0, "P2"=P2, "P0"=P0, "alphaLat"=alphaLat, "betaLat"=betaLat,
                "sens"=sens, "spec"=spec, "FP0"=FP0, "FN2"=FN2)

  } else if (biomType=="continuous") {

    #################################################
    # Computations for a continuous biomarker

    # Define the truebetas (betaLat) indexed by the user-specified vector VElowest.
    # VElowest: a vector of fixed values of VE(x) for the subgroup of subjects with lowest X^* values,
    # where this subgroup has prevalence PlatVElowest

    o <- length(VElowest)
    nu <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
    truebetas <- rep(NA,o)
    alphalatvect <- rep(NA,o)

    for (l in 1:o) {

      # find solutions alphalat and betalat by solving eqn (4) in Appendix B using functions kernel() and alphaLatEqn()
      risk1latnu <- (1-VElowest[l])*risk0

      alphalatvect[l] <- uniroot(alphaLatEqn, lower=-10, upper=10, nu=nu, risk1latnu=risk1latnu, sigma2obs=sigma2obs,
                                 VEoverall=VEoverall, PlatVElowest=PlatVElowest, risk0=risk0)$root

      # Second solve for betalat:
      D <- risk1latnu
      truebetas[l] <- (log(D/(1-D)) - alphalatvect[l])/nu[1]
    }

    # initialize power calculation matrix
    if(max(sampleLengths)>1) {
      powerscont <- matrix(0, nrow=length(Ncomplete), ncol=length(VElowest))
      rownames(powerscont) <- paste0(rep("N"), seq(1,nrow(powerscont)))
    } else {
      powerscont <- matrix(0, nrow=length(rho), ncol=length(VElowest))
      rownames(powerscont) <- paste0(rep("rho"), seq(1,nrow(powerscont)))
    }


    ###################################################
    # Power calculations repeated for M simulations

    for (i in 1:M) {

      # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
      # using the function risk1cont()
      # Step 4 for continuous biomarker in manuscript

      for (j in 1:o) {  # loop through each value of VElowest

        beta <- truebetas[j]
        alphalat <- alphalatvect[j]

        # These simulations condition on nCases (i.e., number of infections in vaccine arm) and
        # also on the number of controls fixed at controlCaseRatio*nCases

        # If multiple sample size specified, loop through the different sample sizes;
        # else, loop through the different RRlat0's
        if(max(sampleLengths)>1){
          for(k in 1:length(Ncomplete)){

            # Arbitrarily put the cases first and controls second
            Y <- c(rep(1,nCases[j]),rep(0,Ncomplete[j]-nCases[j]))

            # Compute the denominator of the density of X|Y=1 when Y|X follows a logistic regression model
            # with the truncated part associated with VElowest and X is normal with mean zero and
            # standard deviation sqrt(sigma2tr)
            f <- function(x) {
              ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
              return(ans)
            }
            denomdensityXcases <- integrate(f,lower=nu,upper=5)$value
            denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowest[j])*risk0

            numerdensXcases <- function(x) {
              num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
              num[x <= nu] <- PlatVElowest*(1-VElowest[j])*risk0
              return(num)
            }
            numerdensXcontrols <- function(x) {
              num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr))
              num[x <= nu] <- PlatVElowest*(1-(1-VElowest[j])*risk0)
              return(num)
            }

            # From a sequence of x* ranging from -3.5 to 3.5, sample with replacement nCases
            # with probability probscases determined by the pdf. Do the same for controls.
            Xpoints <- seq(-3.5,3.5,len=25000)
            probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
            probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)

            Xcases <-    sample(Xpoints,size=nCases[k],prob=probscases,replace=TRUE)
            Xcontrols <- sample(Xpoints,size=Ncomplete[k]-nCases[k],prob=probscontrols,replace=TRUE)
            X <- c(Xcases,Xcontrols)

            # Create the immune response variables for the different degrees of measurement error
            error <- rnorm(Ncomplete[k],mean=0,sd=sqrt(sigma2e))
            S <- X + error

            # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
            keepinds <- biomSubset(Y, Ncomplete[k], nCasesWithS, controlCaseRatio, p, cohort)

            # Those with biomarker data:
            Ycc <- Y[keepinds]
            Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rho)

            # osDesign logistic regression
            fit <- tps(Ycc~Scc[1,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
            pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}

          }

        } else { # loop through the different RRlat0's

          # Arbitrarily put the cases first and controls second
          # The numbers of cases and controls are fixed, e.g., a typical retrospective design
          Y <- c(rep(1,nCases),rep(0,Ncomplete-nCases))

          # Compute the denominator of the density of X|Y=1 when Y|X follows a logistic regression model
          # with the truncated part associated with VElowest and X is normal with mean zero and
          # standard deviation sqrt(sigma2tr)
          for(k in 1:length(nu)){
            f <- function(x) {
              ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
              return(ans)
            }
            denomdensityXcases <- integrate(f,lower=nu[k],upper=5)$value
            denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowest[j])*risk0

            numerdensXcases <- function(x) {
              num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[k]))
              num[x <= nu[k]] <- PlatVElowest*(1-VElowest[j])*risk0
              return(num)
            }
            numerdensXcontrols <- function(x) {
              num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr[k]))
              num[x <= nu[k]] <- PlatVElowest*(1-(1-VElowest[j])*risk0)
              return(num)
            }

            # From a sequence of x* ranging from -3.5 to 3.5, sample with replacement nCases
            # with probability probscases determined by the pdf. Do the same for controls.
            Xpoints <- seq(-3.5,3.5,len=25000)
            probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
            probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)

            Xcases <-    sample(Xpoints,size=nCases,prob=probscases,replace=TRUE)
            Xcontrols <- sample(Xpoints,size=Ncomplete-nCases,prob=probscontrols,replace=TRUE)
            X <- c(Xcases,Xcontrols)

            # Create the immune response variables for the different degrees of measurement error
            error <- rnorm(Ncomplete,mean=0,sd=sqrt(sigma2e[k]))
            S <- X + error

            # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
            keepinds <- biomSubset(Y, Ncomplete, nCasesWithS, controlCaseRatio, p, cohort)

            # Those with biomarker data:
            Ycc <- Y[keepinds]
            Scc <- S[keepinds]

            # osDesign logistic regression
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

    # output list for continuous biomarker
    pwr <- list("power"=power, "RRc"=RRc, "betaLat"=truebetas, "alphaLat"=alphalatvect, "PlatVElowest"=PlatVElowest,
                "VElowest"=VElowest, "sigma2obs"=sigma2obs)
  }

  # global outputs
  pwr$Ncomplete <- Ncomplete
  pwr$nCases <- nCases
  pwr$nCasesWithS <- nCasesWithS
  pwr$VEoverall <- 1-RRoverall
  pwr$alpha <- alpha
  pwr$rho <- rho
  pwr$controlCaseRatio <- controlCaseRatio
  pwr$risk0 <- risk0

  # If saveDir and saveFile specified, save output list to .RData file with given file name and location
  if(!is.null(saveDir) & !is.null(saveFile)) {
    save(pwr, file=paste0(file.path(saveDir, saveFile),".RData"))
  }

  return(pwr)

}
