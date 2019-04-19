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

# nCasesTx     All cases in the vaccine group at-risk at tau and a case by taumax
#                             (regardless of whether the biomarker is measured)
# nControlsTx  All controls in the vaccine group at-risk at tau and not diseased at the end of follow-up taumax
#                             (regardless of whether the biomarker is measured)
# nCasesTxWithS  As above and also have the biomarker measured (i.e., in Phase 2)

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

checkFileSavingParams <- function(saveFile, saveDataFile, varyingParam, varyingParamName) {
  # Checks that the input parameters for file names include ".RData"
  #
  # Args:
  #   saveFile: Character string specifying the name of the file storing the power calculation output
  #   saveDataFile: Character string specifying the name of the file storing the simulated data
  #
  # Returns:
  #   Error if saveFile or saveDataFile is a vector whose length does not match the length of the varying argument
  #   or if a component of saveFile or saveDataFile does not include ".RData"

  if (length(saveDataFile) != 1) {
    if(varyingParamName == "controlCaseRatio" | varyingParamName == "p") {
      stop("If the varying argument is controlCaseRatio or p, saveDataFile must be a character string")
    } else if (length(saveDataFile) != length(varyingParam[[1]])) {
      stop("saveDataFile is a vector and its length is not equal to the length of the varying argument")
    }
  }
  if (length(saveFile) != 1 & length(saveFile) != length(varyingParam[[1]])) {
    stop("saveFile is a vector and its length is not equal to the length of the varying argument")
  }
  for (i in length(saveFile)) {
    fileName <- saveFile[i]
    if (substr(fileName, start = nchar(fileName) - 5, stop = nchar(fileName)) != ".RData") {
      stop("All components of vector for saveFile must include .RData")
    }
  }
  for (i in length(saveDataFile)) {
    fileName <- saveDataFile[i]
    if (substr(fileName, start = nchar(fileName) - 5, stop = nchar(fileName)) != ".RData") {
      stop("All components of saveDataFile must include .RData")
    }
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

checkVectorParamsMatch <- function(vectorParamLength, varyingParamName) {
  # Checks that vector parameter lengths are consistent
  #
  # Args:
  #   vectorParamLength: numeric vector denoting the lengths of the varying input parameters
  #   varyingParamName: character vector of the names of the varying parameters
  #
  # Returns:
  #   Error if input varying parameter components are vectors of different lengths

  if (max(vectorParamLength) > 1) {
    if (max(vectorParamLength) != min(vectorParamLength)) {
      stop(paste0("Vector lengths differ for ", paste0(varyingParamName, collapse=", ")))
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
      stop("Input arguments PlatVElowest and VElowest violate probability constraints for normal biomarker calculations")
    }
  } else if (any(RRlat2 < 0)) {
    cat(paste("RRlat2="),"\n")
    cat(paste(round(RRlat2,3)))
    cat("\n")
    stop("Input arguments violate probability constraints for trichotomous biomarker calculations.
         Consider making Plat0 smaller and/or VElat0 larger.")
  }
  }

checkSaveDataParams <- function(Approach2, biomType, corr, nCasesPla, nControlsPla, sensBIP, specBIP, FP0BIP, FN2BIP, P0BIP, P2BIP) {
  # Checks that sample size input parameters are valid.
  #
  # Args:
  #   Approach2: logical value denoting whether Approach 2 for trichotomous/dichotomous biomarker is to be used
  #   biomType: Character spring specifying type of biomarker used. Default is "continuous";
  #             other choices are "trichotomous" and "dichotomous"
  #   corr: correlation between BIP and S
  #   nCasesPla: all cases in the placebo group at-risk at tau and a case by taumax
  #   nControlsPla: all controls in the placebo group at-risk at tau and not diseased at the end of follow-up taumax
  #   sensBIP: Numeric scalar or vector corresponding to P(BIP=2|S=2)
  #   specBIP: Numeric scalar or vector corresponding to P(BIP=0|S=0)
  #   FP0BIP: Numeric scalar or vector corresponding to P(BIP=2|S=0)
  #   FN2BIP: Numeric scalar or vector corresponding to P(BIP=0|S=2)
  #   P0BIP: Probability of a low value of the BIP
  #   P2BIP: Probability of a high value of the BIP
  #
  # Returns:
  #   Error if nCasesPla, or nControlsPla is not specified,
  #   or if Approach 1 is to be used and sensBIP, specBIP, etc. are not specified,
  #   or if more than one of the input parameters (excluding VElat0, VElat1, and VElowest) is vectorized,
  #   or if Approach 2 or a continuous biomarker is to be used and corr is not specified or is not a valid correlation

  if(is.null(nCasesPla) | is.null(nControlsPla)) {
    stop("If full data is to be saved, input arguments nCasesPla, and nControlsPla must be specified.")
  }

  if(!Approach2) {
    if(is.null(sensBIP) | is.null(specBIP) | is.null(FP0BIP) | is.null(FN2BIP) | is.null(P0BIP) | is.null(P2BIP)) {
      stop("If full data is to be saved for a trichotomous/dichotomous biomarker and Approach 1 is to be used, input arguments sensBIP, specBIP, FP0BIP, FN2BIP, P0BIP, and P2BIP must be specified.")
    }
    sensSpecBIPlengths <- sapply(list(sensBIP, specBIP, FP0BIP, FN2BIP), length)
    PxBIPlengths <- sapply(list(P0BIP, P2BIP), length)
    checkVectorParamsMatch(sensSpecBIPlengths, c("sensBIP, specBIP, FP0BIP, FN2BIP"))
    checkVectorParamsMatch(PxBIPlengths, c("P0BIP, P2BIP"))
    if(sum(c(max(sensSpecBIPlengths), max(PxBIPlengths), length(corr)) > 1) > 1) {
      stop("Only one out of the three groups: 1) sensBIP, specBIP, FP0BIP, FN2BIP, 2) P0BIP, P2BIP, and 3) corr may be vectorized at a time.")
    }
  } else {
    if (is.null(corr)) {
      if (biomType == "continuous") {
        stop("If full data is to be saved for a continuous biomarker, input argument corr must be specified.")
      } else {
        stop("If full data is to be saved for a trichotomous/dichotomous biomarker and Approach 2 is to be used, input argument corr must be specified.")
      }
    }
    for(i in 1:length(corr)) {
      if(corr[i] < -1 | corr[i] > 1) {
        stop("Input argument corr is not a valid correlation.")
      }
    }
  }
}

getVaryingParam <- function(sampleSizes, Platx, sensSpec, controlCaseRatio, p, rho, PlatVElowest) {
  # Gets the input parameter whose values represent distinct scenarios for power assessment
  #
  # Args:
  #   sampleSizes: list containing nCasesTx, nControlsTx, and nCasesTxWith S
  #     (also includes nCasesPla and nControlsPla if saveDataFile is specified)
  #   Platx: list containing Plat0, Plat2, P0, and P2
  #   sensSpec: list containing sens, spec, FP0, and FN2
  #   controlCaseRatio: Ratio of controls to cases
  #   p: Probability a subject will be in the cohort
  #   rho: Protection-relevant fraction of the variance of S*
  #   PlatVElowest: Prevalence of VElowest
  #
  # Returns:
  #   A list of the name(s) of the input parameter(s) whose vector components represent distinct scenarios,
  #   and a list of the values of the input parameter(s) whose vector components represent distinct scenarios.
  #   If no such input parameter exists, returns '3' as a placeholder element.

  sampleSizesLengths <- sapply(sampleSizes, length)
  PlatxLengths <- sapply(Platx, length)
  sensSpecLengths <- sapply(sensSpec, length)

  allLengths <- c(max(sampleSizesLengths), max(PlatxLengths), max(sensSpecLengths), length(controlCaseRatio), length(p), length(rho),
                  length(PlatVElowest))
  if(sum(allLengths > 1) > 1) {
    stop("Excluding VElat0, VElat1, and VElowest, only one group of input arguments may be varied at a time")
  }

  varyingParamName <- ""
  varyingParam <- list()
  if (max(sampleSizesLengths) > 1) {
    varyingParamName <- names(sampleSizes)
    varyingParam <- sampleSizes
    checkVectorParamsMatch(sampleSizesLengths, varyingParamName)
  } else if (max(PlatxLengths) > 1) {
    varyingParamName <- names(Platx)
    varyingParam <- Platx
    checkVectorParamsMatch(PlatxLengths, varyingParamName)
  } else if (max(sensSpecLengths) > 1) {
    varyingParamName <- names(sensSpec)
    varyingParam <- sensSpec
    checkVectorParamsMatch(sensSpecLengths, varyingParamName)
  } else if (length(controlCaseRatio) > 1) {
    varyingParamName <- "controlCaseRatio"
    varyingParam$controlCaseRatio <- controlCaseRatio
  } else if (length(p) > 1) {
    varyingParamName <- "p"
    varyingParam$p <- p
  } else if (length(rho) > 1) {
    varyingParamName <- "rho"
    varyingParam$rho <- rho
  } else if (length(PlatVElowest) > 1) {
    varyingParamName <- "PlatVElowest"
    varyingParam$PlatVElowest <- PlatVElowest
  } else {
    varyingParamName <- 3
    varyingParam <- 3
  }

  return(list("varyingParamName" = varyingParamName, "varyingParam" = varyingParam))
}


computeRisks <- function(biomType, RRoverall, risk0, sens, spec, FP0, FN2, FP1, FN1,
                         Plat0, Plat1, Plat2, P0, P1, P2, RRlat0, RRlat1, RRlat2) {
  # Computes the endpoint risks conditional on each latent subgroup and each observed subgroup in the active treatment arm
  #
  # Args:
  #   FP1: Numeric scalar or vector corresponding to P(S=2|X=1)
  #   FN1: Numeric scalar or vector corresponding to P(S=0|X=1)
  #   P1: P(S=1)
  #   RRoverall: 1 - VEoverall
  #   RRlat0: 1 - VElat0
  #   RRlat1: 1 - VElat1
  #   RRlat2: 1 - VElat2
  #
  # Returns:
  #   A named list containing the endpoint risks conditional on each observed subgroup in the active treatment arm
  #   and the endpoint risks conditional on each latent protected subgroup in the active treatment arm

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
  # get vectors with length=length(RRlat0)
  risk1_2 <- (probX0_cond_S2 * RRlat0 + probX1_cond_S2 * RRlat1 + probX2_cond_S2 * RRlat2 )*risk0
  probX0_cond_S0 <- spec*Plat0/P0
  probX1_cond_S0 <- FN1*Plat1/P0
  probX2_cond_S0 <- FN2*Plat2/P0
  risk1_0 <- (probX0_cond_S0 * RRlat0 + probX1_cond_S0 * RRlat1 + probX2_cond_S0 * RRlat2)*risk0
  risk1_1 <- (risk1 - risk1_0*P0 - risk1_2*P2)/P1  # Note: For the dichotomous biomarker special case, the risk1medx are NA

  esvect <- risk1_2/risk1_0  # vector with length=length(RRlat0)

  # Vaccine risks within the latent subgroups (independent of rho of course)
  risk1lat_2 <- RRlat2*risk0
  risk1lat_1 <- RRlat1*risk0
  risk1lat_0 <- RRlat0*risk0

  return(list("risk1_0" = risk1_0, "risk1_1" = risk1_1, "risk1_2" = risk1_2,
              "risk1lat_0" = risk1lat_0, "risk1lat_1" = risk1lat_1, "risk1lat_2" = risk1lat_2))
}




simTrich <- function(risk1lat_0, risk1lat_1, risk1lat_2, Plat0, Plat1, Plat2, sens, spec, FP0, FN2, FP1, FN1,
                     nCasesTx, nCasesTxWithS, NcompleteTx, nCasesPla, NcompletePla,
                     biomType, controlCaseRatio, p, cohort, tpsMethod, alpha, saveDataDir,
                     sensBIP, specBIP, FP0BIP, FN2BIP, FP1BIP, FN1BIP, varyFullData) {
  # Simulates data and conducts a hypothesis test for a trichotomous or dichotomous biomarker
  #
  # Args:
  #   risk1lat_0: endpoint risk conditional on the latent lower protected response in the active treatment arm
  #   risk1lat_1: endpoint risk conditional on the latent lower protected response in the active treatment arm
  #   risk1lat_2: endpoint risk conditional on the latent lower protected response in the active treatment arm
  #   varyFullData: logical value indicating if full data is to be simulated multiple times (FALSE if varyingParam
  #     is controlCaseRatio or p)
  #
  # Returns:
  #   A list containing a binary integer denoting the outcome of the hypothesis test.
  #   If simulated data is to be saved, the returned list also includes a data frame containing the simulated
  #   data, including placebo group and BIP data

  # Determine success probabilities for trinomial random variable:
  # P(X=0|Y=1, Y^tau=0, Z=1), P(X=1|Y=1, Y^tau=0, Z=1), P(X=2|Y=1, Y^tau=0, Z=1),
  # using Bayes rule to express them in terms of Platx and risk1(x), and risk1
  rrlat0 <- risk1lat_0/(risk1lat_0+risk1lat_1+risk1lat_2)  # risk1(0)
  rrlat1 <- risk1lat_1/(risk1lat_0+risk1lat_1+risk1lat_2)
  rrlat2 <- risk1lat_2/(risk1lat_0+risk1lat_1+risk1lat_2)
  denominat <- Plat0*rrlat0 + Plat1*rrlat1 + Plat2*rrlat2  # risk1
  P0case <- (Plat0*rrlat0)/denominat  # success probabilities for trinomial random variable
  P1case <- (Plat1*rrlat1)/denominat
  P2case <- 1 - P0case - P1case

  # Draw from trinomial random variable with success probabilities defined above.
  # adjustProb() function deals with rare crashes of rmultinom due to numerical problems
  # where the program treats probability 0 as a small negative number
  indsTx <- rmultinom(nCasesTx,1,adjustProb(c(P0case,P1case,P2case)))

  # Number of cases in the 0, 1, 2 latent groups
  nCasesTx0 <- length(indsTx[1,][indsTx[1,]==1])
  nCasesTx2 <- length(indsTx[3,][indsTx[3,]==1])
  nCasesTx1 <- nCasesTx - nCasesTx0 - nCasesTx2
  Ntx0 <- round(Plat0*NcompleteTx)
  Ntx2 <- round(Plat2*NcompleteTx)
  Ntx1 <- NcompleteTx - Ntx0 - Ntx2

  # Address rounding that could make Ntx1 negative in the dichotomous marker case
  # Keep NcompleteTx fixed at a constant
  if (Ntx1==-1) {
    Ntx0 <- Ntx0 + 1
    Ntx1 <- 0
  }
  # Also keep nCasesTx fixed at a constant
  if (nCasesTx1==-1) {
    nCasesTx0 <- nCasesTx0 - 1
    nCasesTx1 <- 0
  }
  if (nCasesTx1==1 & Ntx1==0) {
    nCasesTx0 <- nCasesTx0 + 1
    nCasesTx1 <- 0
  }

  # Latent subgroup assignments in active treatment arm
  Xtx <- c(rep(0, Ntx0), rep(1, Ntx1), rep(2, Ntx2))

  # Fix the number of cases and controls, putting the cases first and controls second for each subgroup:
  Ytx <- c(rep(1,nCasesTx0),rep(0,Ntx0-nCasesTx0),rep(1,nCasesTx1),rep(0,Ntx1-nCasesTx1),rep(1,nCasesTx2),rep(0, NcompleteTx - Ntx0 - Ntx1 - nCasesTx2))

  # Simulate the trinary surrogate
  # Formulas (12) and (13) in the manuscript:

  # Given specifications for spec, FP0, sens, and FN2 and a logical value indicating if the biomarker is dichotomous
  # or not, the function assignBiomarkerLevels() returns a vector composed of biomarker levels (S=0,1,2),
  # where each subject is assigned a specific level
  if (biomType=="dichotomous") { # dichotomous
    Stx <- assignBiomarkerLevels(sens, spec, FP0, FN2, FP1, FN1, dichotomous=TRUE, Ntx0, Ntx1, Ntx2)
  } else { # trichotomous
    Stx <- assignBiomarkerLevels(sens, spec, FP0, FN2, FP1, FN1, dichotomous=FALSE, Ntx0, Ntx1, Ntx2)
  }

  # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
  keepinds <- biomSubset(Ytx, NcompleteTx, nCasesTxWithS, controlCaseRatio, p, cohort)

  # Those with biomarker data:
  Ycc <- Ytx[keepinds]
  Scc <- Stx[keepinds]

  ##############################################################
  # Now analyze with osDesign
  # (first check if there are 'zeros', in which case Fisher's exact test for the lo vs. hi categories is used.
  # Otherwise, osDesign logistic regression is used as an ordered score test

  addPower <- 0
  lodim <- dim(table(Ycc,Scc))[2]<2 # check there are at least two biomarker categories (columns)
  zerosflag <-  lodim
  if (dim(table(Ycc,Scc))[2]==3) { # check if any categories have zero entries
    zerosflag <- table(Ycc,Scc)[1,1]==0 | table(Ycc,Scc)[1,2]==0 | table(Ycc,Scc)[1,3]==0 |
      table(Ycc,Scc)[2,1]==0 | table(Ycc,Scc)[2,2]==0 | table(Ycc,Scc)[2,3]==0
  }

  if (zerosflag) {
    if (lodim) { pval <- 1}
    if (!lodim) { # there are zeros, so Fisher's exact test is used
      pval <- fisher.test(table(Ycc,Scc)[,c(1,dim(table(Ycc,Scc))[2])])$p.value
    }
    if (pval <= alpha & length(Ycc[Scc==2&Ycc==1])/length(Scc[Scc==2]) <
        length(Ycc[Scc==0&Ycc==1])/length(Scc[Scc==0])) {
      addPower <- 1
      # powerstrinary <- powerstrinary + 1
    }
  }

  if (!zerosflag) {
    fit <- tps(Ycc~Scc,nn0=length(Ytx[Ytx==0]),nn1=length(Ytx[Ytx==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
    pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
    if (pval <= alpha & fit$coef[2] < 0) {
      addPower <- 1
      # powerstrinary <- powerstrinary + 1
    }
  }

  output <- list("addPower" = addPower)

  if(!is.null(saveDataDir) & varyFullData) {
    ################################################
    # Generate simulated X, Y, and S for placebo group.

    # Draw from trinomial random variable with success probabilities Plat0, Plat1, and Plat2
    indsPla <- rmultinom(nCasesPla,1,adjustProb(c(Plat0,Plat1,Plat2)))

    # Number of cases in the 0, 1, 2 latent groups
    nCasesPla0 <- length(indsPla[1,][indsPla[1,]==1])
    nCasesPla2 <- length(indsPla[3,][indsPla[3,]==1])
    nCasesPla1 <- nCasesPla - nCasesPla0 - nCasesPla2

    # Number of subjects in the 0, 1, 2 latent groups
    Npla0 <- round(Plat0*NcompletePla)
    Npla2 <- round(Plat2*NcompletePla)
    Npla1 <- NcompletePla - Npla0 - Npla2

    # Address rounding that could make Npla1 negative in the dichotomous marker case
    # Keep NcompletePla fixed at a constant
    if (Npla1==-1) {
      Npla0 <- Npla0 + 1
      Npla1 <- 0
    }
    # Also keep nCasesPla fixed at a constant
    if (nCasesPla1==-1) {
      nCasesPla0 <- nCasesPla0 - 1
      nCasesPla1 <- 0
    }
    if (nCasesPla1==1 & Npla1==0) {
      nCasesPla0 <- nCasesPla0 + 1
      nCasesPla1 <- 0
    }

    # Latent subgroup assignments in placebo arm
    Xpla <- c(rep(0, Npla0), rep(1, Npla1), rep(2, Npla2))

    # Endpoint indicator variable for the placebo group
    # Fix the number of cases and controls, putting the cases first and controls second for each subgroup:
    Ypla <- c(rep(1,nCasesPla0),rep(0,Npla0-nCasesPla0),rep(1,nCasesPla1),rep(0,Npla1-nCasesPla1),rep(1,nCasesPla2),rep(0, NcompletePla - Npla0 - Npla1 - nCasesPla2))

    # Simulate the trinary surrogate
    # Formulas (12) and (13) in the manuscript:

    # Given specifications for spec, FP0, sens, and FN2 and a logical value indicating if the biomarker is dichotomous
    # or not, the function assignBiomarkerLevels() returns a vector composed of biomarker levels (S=0,1,2),
    # where each subject is assigned a specific level
    if (biomType=="dichotomous") { # dichotomous
      Spla <- assignBiomarkerLevels(sens, spec, FP0, FN2, FP1, FN1, dichotomous=TRUE, Npla0, Npla1, Npla2)
    } else { # trichotomous
      Spla <- assignBiomarkerLevels(sens, spec, FP0, FN2, FP1, FN1, dichotomous=FALSE, Npla0, Npla1, Npla2)
    }

    #### Simulate a trichotomous BIP ####

    # S1 denotes biomarker observed under assignment to treatment (either at randomization or after crossover)
    S1 <- c(Stx, Spla)

    # Number of subjects in the 0, 1, 2 observed biomarker groups
    Ns0 <- sum(S1==0)
    Ns1 <- sum(S1==1)
    Ns2 <- sum(S1==2)

    # If sensBIP, specBIP, etc. are vectors, simulate multiple BIPs,
    # with each BIP corresponding to one round of sensitivity/specificity parameters
    BIP <- matrix(0, nrow=NcompleteTx+NcompletePla, ncol=length(sensBIP))
    colnames(BIP) <- paste0("BIP", 1:length(sensBIP))
    for (m in 1:length(sensBIP)) {

      # Given specifications for sensBIP, specBIP, FP0BIP, FN2BIP, FP1BIP, and FN1BIP and a logical value indicating
      # if the biomarker is dichotomous or not. The function assignBiomarkerLevels() returns a vector composed of
      # the levels of the BIP (BIP=0,1,2), where each subject is assigned a specific level
      if (biomType=="dichotomous") { # dichotomous
        BIPdata <- assignBiomarkerLevels(sensBIP[m], specBIP[m], FP0BIP[m], FN2BIP[m], FP1BIP[m], FN1BIP[m], dichotomous=TRUE, Ns0, Ns1, Ns2)
      } else { # trichotomous
        BIPdata <- assignBiomarkerLevels(sensBIP[m], specBIP[m], FP0BIP[m], FN2BIP[m], FP1BIP[m], FN1BIP[m], dichotomous=FALSE, Ns0, Ns1, Ns2)
      }
      BIP[S1==0, m] <- BIPdata[1:Ns0]
      BIP[S1==1, m] <- BIPdata[(Ns0+1):(Ns0+Ns1)]
      BIP[S1==2, m] <- BIPdata[-(1:(Ns0+Ns1))]
    }

    #### Gather all data needed for full data output ####
    X <- c(Xtx, Xpla)
    Y <- c(Ytx, Ypla)
    Z <- c(rep(1, NcompleteTx), rep(0, NcompletePla))
    simData <- data.frame(X, Y, Z, S1, BIP)
    output$simData <- simData
  }

  return(output)
}

simCont <- function(nCasesTx, NcompleteTx, nCasesTxWithS, nCasesPla, NcompletePla,
                    alphalat, beta, sigma2obs, sigma2tr, sigma2e, nu, PlatVElowest, VElowest,
                    risk0, controlCaseRatio, p, cohort, tpsMethod, alpha, corr, saveDataDir, varyFullData) {
  # Simulates data and conducts a hypothesis test for a continuous biomarker
  #
  # Args:
  #   alphalat: parameter in logistic regression model: logit(risk1lat(x*)) = alphalat + beta x*
  #   beta: parameter in logistic regression model: logit(risk1lat(x*)) = alphalat + beta x*
  #   sigma2obs: Variance of observed biomarker S*
  #   sigma2tr: Variance of latent response variable X*
  #   sigma2e: Variance of the error term in the classical measurement error model: S* = X* + e
  #   nu: Value that denotes the cut-off for lowest X* values. Fraction PlatVElowest of subjects have x* value below this
  #       and therefore have VE = VElowest. nu = sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  #   varyFullData: logical value indicating if full data is to be simulated multiple times (FALSE if varyingParam
  #     is controlCaseRatio or p)
  #
  # Returns:
  #   A list containing a binary integer denoting the outcome of the hypothesis test.
  #   If simulated data is to be saved, the returned list also includes a data frame containing the simulated
  #   data, including placebo group and BIP data

  # Arbitrarily put the cases first and controls second
  # The numbers of cases and controls are fixed, e.g., a typical retrospective design
  Ytx <- c(rep(1,nCasesTx),rep(0,NcompleteTx-nCasesTx))

  # Compute the denominator of the density of X|Y=1 when Y|X follows a logistic regression model
  # with the truncated part associated with VElowest and X is normal with mean zero and
  # standard deviation sqrt(sigma2tr)
  f <- function(x) {
    ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
    return(ans)
  }
  denomdensityXcases <- integrate(f,lower=nu,upper=5)$value
  denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowest)*risk0

  numerdensXcases <- function(x) {
    num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr))
    num[x <= nu] <- PlatVElowest*(1-VElowest)*risk0
    return(num)
  }
  numerdensXcontrols <- function(x) {
    num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr))
    num[x <= nu] <- PlatVElowest*(1-(1-VElowest)*risk0)
    return(num)
  }

  # From a sequence of x* ranging from -3.5 to 3.5, sample with replacement nCasesTx
  # with probability probscases determined by the pdf. Do the same for controls.
  Xpoints <- seq(-3.5,3.5,len=25000)
  probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
  probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)

  Xcases <-    sample(Xpoints,size=nCasesTx,prob=probscases,replace=TRUE)
  Xcontrols <- sample(Xpoints,size=NcompleteTx-nCasesTx,prob=probscontrols,replace=TRUE)
  Xtx <- c(Xcases,Xcontrols)

  # Create the immune response variables for the different degrees of measurement error
  error <- rnorm(NcompleteTx,mean=0,sd=sqrt(sigma2e))
  Stx <- Xtx + error

  # Select subset of subjects with biomarker measured (R_i=1) according to case-cohort or case-control sampling design
  keepinds <- biomSubset(Ytx, NcompleteTx, nCasesTxWithS, controlCaseRatio, p, cohort)

  # Those with biomarker data:
  Ycc <- Ytx[keepinds]
  Scc <- Stx[keepinds]

  addPower <- 0
  # osDesign logistic regression
  fit <- tps(Ycc~Scc,nn0=length(Ytx[Ytx==0]),nn1=length(Ytx[Ytx==1]),group=rep(1,length(Ycc)), method=tpsMethod, cohort=cohort)
  pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
  if (pval <= alpha & fit$coef[2] < 0) {
    addPower <- 1
    # powerscont[k,j] <- powerscont[k,j] + 1
  }

  output <- list("addPower" = addPower)

  if(!is.null(saveDataDir) & varyFullData) {
    ##############################
    # Simulations for placebo group

    # Arbitrarily put the cases first and controls second
    Ypla <- c(rep(1,nCasesPla),rep(0,NcompletePla-nCasesPla))

    # From a sequence of x* ranging from -3.5 to 3.5, sample with replacement NcompletePla
    # with probability determined by P(X* = x*).
    Xpla <- sample(Xpoints, size=NcompletePla, prob=dnorm(Xpoints/sqrt(sigma2tr)), replace=TRUE)

    # Create the biomarker response variable
    errorPla <- rnorm(NcompletePla,mean=0,sd=sqrt(sigma2e))
    Spla <- Xpla + errorPla

    ### Simulate the baseline immunogenicity predictor (BIP)

    # S1 denotes biomarker observed under assignment to treatment (either at randomization or after crossover)
    S1 <- c(Stx, Spla)

    BIP <- matrix(0, nrow=NcompleteTx+NcompletePla, ncol=length(corr))
    colnames(BIP) <- paste0("BIP", 1:length(corr))
    for (m in 1:length(corr)) {
      sigma2d <- (sigma2obs / corr[m]^2) - sigma2tr - sigma2e
      BIP[,m] <- S1 + rnorm(NcompleteTx + NcompletePla, mean = 0, sd = sqrt(sigma2d))
    }

    ### Full data output
    X <- c(Xtx, Xpla)
    Y <- c(Ytx, Ypla)
    Z <- c(rep(1, NcompleteTx), rep(0, NcompletePla))
    simData <- data.frame(X, Y, Z, S1, BIP)
    output$simData <- simData
  }
  return(output)
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

assignBiomarkerLevels <- function(sens, spec, FP0, FN2, FP1, FN1, dichotomous, N0, N1, N2){
  # Assigns an observed biomarker level (S=0, 1, or 2) to each subject.
  #
  # Args:
  #   dichotomous: If TRUE, indicates biomarker is dichotomous; if FALSE, indicates biomarker is trichotomous
  #   N0: Number of subjects at risk at tau in the lower protected latent subgroup, excluding dropouts
  #   N1: Number of subjects at risk at tau in the medium protected latent subgroup, excluding dropouts
  #   N2: Number of subjects at risk at tau in the higher protected latent subgroup, excluding dropouts
  #
  # Returns:
  #   Vector composed of biomarker levels (0,1,2), where each subject is assigned a specific level

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

biomSubset <- function(Y, NcompleteTx, nCasesTxWithS, controlCaseRatio, p, cohort){
  # Selects subset of subjects that have biomarker S or S* measured (R_i=1) according to
  # a case-cohort or case-control sampling design.
  #
  # Args:
  #   Y: Numeric vector indicating cases and controls (1 vs. 0), with length = NcompleteTx
  #   NcompleteTx: Total number of subjects at risk at tau, excluding dropouts
  #   nCasesTxWithS: Number of observed cases between tau and taumax with measured S or S*
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
    R <- ifelse(rbinom(NcompleteTx, 1, p)==1, 1, R) # from (NcompleteTx=nCasesTx+nControlsTx), draw Bernoulli sample with sampling probability p
    R <- ifelse(Y==1, 1, R)  # augment all cases
    keepinds <- which(R==1)

  } else {  # case-control sampling design

    # Keep the S's in nCasesTxWithS of the cases (deleting the rest) and in controlCaseRatio*nCasesTxWithS controls
    casesinds <- which(Y==1)
    keepcasesinds <- sample(casesinds,nCasesTxWithS,replace=FALSE)
    controlinds <- which(Y==0)
    keepcontrolinds <- sample(controlinds,controlCaseRatio*nCasesTxWithS,replace=FALSE)
    keepinds <- sort(c(keepcasesinds,keepcontrolinds))
  }
  return(keepinds)
}


#' Power Calculations for Assessing Intermediate Biomarkers as Correlates of Risk in the Active Treatment Group in Clinical Efficacy Trials, Accounting for Biomarker's Measurement Error and Treatment Efficacy
#'
#' Performs a power calculation for assessing a univariate dichotomous, trichotomous, or continuous intermediate biomarker response as a correlate of risk
#' in the active treatment group in a clinical efficacy trial, accounting for the biomarker's measurement error and treatment efficacy. The statistical methods are described in [Gilbert, Janes, and Huang (2016).
#' "Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials."] Simulated data sets, extended to include placebo group and baseline immunogenicity predictor data, can be exported for
#' harmonized assessment of biomarker-specific treatment efficacy.
#'
#' @param nCasesTx an integer vector specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group. Each value represents a distinct scenario for power assessment.
#' @param nControlsTx an integer vector specifying the number of controls observed (or projected) to complete follow-up through \eqn{\tau_{max}} endpoint-free in the active treatment group. Each value represents a distinct scenario for power assessment. The ordering in \code{nCasesTx} and \code{nControlsTx} must match.
#' @param nCasesTxWithS an integer vector specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response. Each value represents a distinct scenario for power assessment. The ordering must match \code{nCasesTx} and \code{nControlsTx}.
#' @param controlCaseRatio an integer vector specifying the number of controls sampled per case for biomarker measurement in the without replacement case-control sampling design. Each value represents a distinct scenario for power assessment.
#' @param VEoverall a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#' @param risk0 a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#' @param VElat0 a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lower protected subgroup for a dichotomous or trichotomous biomarker. Each value of \code{VElat0} corresponds to one unique effect size (\eqn{RR_t}). Default ranges from \code{VEoverall} (\eqn{H_0}) to 0 (maximal \eqn{H_1} not allowing harm by treatment).
#' @param VElat1 a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent medium protected subgroup for a trichotomous biomarker. Each value corresponds to one unique effect size (\eqn{RR_t}). The ordering must match \code{VElat0}. Set to \code{VEoverall} by default, but must be set to \code{NULL} for a dichotomous biomarker.
#' @param VElowest a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup for a continuous biomarker. Default ranges from \code{VEoverall} (\eqn{H_0}) to 0 (maximal \eqn{H_1} not allowing harm by treatment).
#' @param Plat0 a numeric vector specifying the prevalence of the latent lower protected subgroup for a dichotomous or trichotomous biomarker. Each value represents a distinct scenario for power assessment. The ordering in \code{Plat0}, \code{Plat2}, \code{P0}, and \code{P2} must match.
#' @param Plat2 a numeric vector specifying the prevalence of the latent higher protected subgroup for a dichotomous or trichotomous biomarker. Each value represents a distinct scenario for power assessment. The ordering in \code{Plat0}, \code{Plat2}, \code{P0}, and \code{P2} must match.
#' @param P0 a numeric vector specifying the probability of low biomarker response for a dichotomous or trichotomous biomarker (set to \code{Plat0} by default). Each value represents a distinct scenario for power assessment. The ordering in \code{Plat0}, \code{Plat2}, \code{P0}, and \code{P2} must match.
#' @param P2 a numeric vector specifying the probability of high biomarker response for a dichotomous or trichotomous biomarker (set to \code{Plat2} by default). Each value represents a distinct scenario for power assessment. The ordering in \code{Plat0}, \code{Plat2}, \code{P0}, and \code{P2} must match.
#' @param PlatVElowest a numeric vector specifying the prevalence of the latent lowest-efficacy subgroup for a continuous biomarker. Each value represents a distinct scenario for power assessment.
#' @param sens a numeric vector specifying the sensitivity, i.e., the probability of high biomarker response conditional on membership in the higher protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'. Each value represents a distinct scenario for power assessment. The ordering in \code{sens}, \code{spec}, \code{FP0}, and \code{FN2} must match.
#' @param spec a numeric vector specifying the specificity, i.e., the probability of low biomarker response conditional on membership in the lower protected subgroup, of a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'. Each value represents a distinct scenario for power assessment. The ordering in \code{sens}, \code{spec}, \code{FP0}, and \code{FN2} must match.
#' @param FP0 a numeric vector specifying the false positive rate, i.e., the probability of high biomarker response conditional on membership in the lower protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'. Each value represents a distinct scenario for power assessment. The ordering in \code{sens}, \code{spec}, \code{FP0}, and \code{FN2} must match.
#' @param FN2 a numeric vector specifying the false negative rate, i.e., the probability of low biomarker response conditional on membership in the higher protected subgroup, for a dichotomous or trichotomous biomarker. Default is \code{NULL}, which indicates the use of 'approach 2'. Each value represents a distinct scenario for power assessment. The ordering in \code{sens}, \code{spec}, \code{FP0}, and \code{FN2} must match.
#' @param M an integer value specifying the number of simulated clinical trials. Default is \code{100}.
#' @param alpha a numeric value specifying the two-sided Wald test type-I error rate. Default is \code{0.05}.
#' @param sigma2obs a numeric value specifying the variance of the observed continuous biomarker or of the dichotomous or trichotomous biomarker simulated using 'approach 2' (set to \code{1} by default).
#' @param rho a numeric vector specifying distinct protection-relevant fractions of \code{sigma2obs}. Each value represents a distinct scenario for power assessment.
#' @param biomType a character string specifying the biomarker type. Default is \code{continuous}; other choices are \code{dichotomous} and \code{trichotomous}.
#' @param cohort a logical value for whether a case-cohort Bernoulli sampling design is to be used. If \code{FALSE} (default), the case-control without replacement sampling is used.
#' @param p a numeric vector specifying the probability of sampling into the subcohort in the case-cohort design (\code{NULL} by default). Each value represents a distinct scenario for power assessment.
#' @param tpsMethod a character string specifying the estimation method in the inverse probability weighted logistic regression model fit by the \code{tps} function in the \code{osDesign} package. The options are \code{PL} for pseudo-likelihood (default), \code{ML} for maximum likelihood, and \code{WL} for weighted likelihood.
#' @param saveDir a character string specifying the path for a directory in which the output of the power calculation is to be saved. If \code{NULL} (default), the output is returned only.
#' @param saveFile a character vector specifying the name(s) of the \code{.RData} file(s) storing the output of the power calculation, used only if \code{saveDir} is not \code{NULL}. All file names must include ".RData" at the end. Default is \code{CoRpower.RData}.
#' @param saveDataDir a character string specifying the path for a directory in which the simulated data, including placebo group and baseline immunogenicity predictor (BIP) data, are to be saved. If \code{NULL} (default), the simulated data are not saved.
#' @param saveDataFile a character vector specifying the name(s) of the \code{.RData} file(s) in which the simulated data, including placebo group and BIP data, are to be saved; used only if \code{saveDataDir} is not \code{NULL}. All file names must include ".RData" at the end. Default is \code{fullData.RData}.
#' @param corr a numeric vector in \eqn{[-1,1]} specifying the correlation between a continuous baseline immunogenicity predictor (BIP) and the (underlying) continuous intermediate biomarker response (\code{NULL} by default). Each value represents a distinct scenario for power assessment. A useful BIP is highly correlated with the biomarker response at \eqn{\tau}. It must be provided if \code{saveDataDir} is specified and a trichotomous biomarker under 'approach 2' or a continuous biomarker is considered.
#' @param nCasesPla an integer vector specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the placebo group. Each value represents a distinct scenario matching \code{nCasesTx}. Default is \code{NULL}. It must be provided if \code{saveDataDir} is specified.
#' @param nControlsPla an integer vector specifying the number of controls observed (or projected) to complete follow-up through \eqn{\tau_{max}} endpoint-free in the placebo group. Each value represents a distinct scenario matching \code{nControlsTx}. Default is \code{NULL}. It must be provided if \code{saveDataDir} is specified.
#' @param sensBIP a numeric vector specifying "the sensitivity" of a dichotomous or trichotomous BIP, i.e., the probability of a high value of the BIP conditional on high biomarker response. Default is \code{NULL}, which indicates the use of 'approach 2'. It must be provided if \code{saveDataDir} is specified and 'approach 1' is to be used. Each value results in generating a separate BIP variable in the output data.
#' @param specBIP a numeric vector specifying "the specificity" of a dichotomous or trichotomous BIP, i.e., the probability of a low value of the BIP conditional on low biomarker response. Default is \code{NULL}, which indicates the use of 'approach 2'. It must be provided if \code{saveDataDir} is specified and 'approach 1' is to be used. Each value results in generating a separate BIP variable in the output data.
#' @param FP0BIP a numeric vector specifying "the false positive rate" of a dichotomous or trichotomous BIP, i.e., the probability of a high value of the BIP conditional on low biomarker response. Default is \code{NULL}, which indicates the use of 'approach 2'. It must be provided if \code{saveDataDir} is specified and 'approach 1' is to be used. Each value results in generating a separate BIP variable in the output data.
#' @param FN2BIP a numeric vector specifying "the false negative rate" of a dichotomous or trichotomous BIP, i.e., the probability of a low value of the BIP conditional on high biomarker response. Default is \code{NULL}, which indicates the use of 'approach 2'. It must be provided if \code{saveDataDir} is specified and 'approach 1' is to be used. Each value results in generating a separate BIP variable in the output data.
#' @param P0BIP a numeric vector specifying the probability of a low value of a dichotomous or trichotomous BIP. If unspecified, it is set to \code{P0}. Each value results in generating a separate BIP variable in the output data.
#' @param P2BIP a numeric vector specifying the probability of a high value of a dichotomous or trichotomous BIP. If unspecified, it is set to \code{P2}. Each value results in generating a separate BIP variable in the output data.
#'
#' @details
#' A number of calling arguments can be specified as vectors with each component specifying a distinct scenario for power assessment (saved in a separate \code{.RData} file).
#' These are referred to as "varying arguments."
#' Some varying arguments occur in group, where the length and order of all specified vectors in the group must match; others are the only varying argument in their group.
#' Only arguments belonging to a single group may be varied at a time; if two or more groups contain vector inputs, the function will treat such inputs as an error.
#' The following are the groups of varying arguments that can be vectorized:
#'   \itemize{
#'     \item \code{nCasesTx}, \code{nControlsTx}, and \code{nCasesTxWithS} (together with \code{nCasesPla} and \code{nControlsPla} if simulated data sets are to be saved)
#'     \item \code{Plat0}, \code{Plat2}, \code{P0}, and \code{P2}
#'     \item \code{sens}, \code{spec}, \code{FP0}, and \code{FN2}
#'     \item \code{controlCaseRatio}
#'     \item \code{rho}
#'     \item \code{p}
#'   }
#'
#' Arguments independent of biomarker type and sampling design: \code{nCasesTx}, \code{nControlsTx}, \code{nCasesTxWithS}, \code{VEoverall}, \code{risk0},
#' \code{M}, \code{alpha}, \code{tpsMethod}, \code{saveDir}, \code{saveFile}.
#'
#' Arguments specific to a trichotomous (or dichotomous) biomarker response: \code{VElat0}, \code{VElat1}, \code{Plat0}, \code{Plat2}, \code{P0},
#' \code{P2}, \code{biomType = "trichotomous"} (or \code{"dichotomous"})
#'   \itemize{
#'     \item Arguments for Approach 1: \code{sens}, \code{spec}, \code{FP0}, \code{FN2}
#'     \item Arguments for Approach 2: \code{sigma2obs}, \code{rho}
#'   }
#'
#' Arguments specific to a continuous biomarker response: \code{VElowest}, \code{PlatVElowest}, \code{sigma2obs}, \code{rho}, \code{biomType = "continuous"}
#'
#' Arguments for a case-control without replacement sampling design: \code{controlCaseRatio}
#'
#' Arguments for a case-cohort Bernoulli sampling design: \code{cohort = TRUE}, \code{p}
#'
#' To save output from the power calculations in an \code{.RData} file, \code{saveDir} must be specified. The default file name is \code{CoRpower.RData};
#' a different file name may be specified by \code{saveFile} as a single character string, to which the value of the varying argument(s) will be appended for descriptive file naming purposes,
#' or, alternatively, a character vector may be specified with full file names (a single file will be produced for each value of the varying argument(s)).
#'
#' To link power calculations for detecting a correlate of risk and a correlate of treatment efficacy, simulated data sets used in the power calculations can be exported
#' with extensions including placebo-group and BIP data for harmonized use by methods assessing biomarker-specific treatment efficacy.
#' The vignette "Algorithms for Simulating Placebo Group and Baseline Immunogenicity Predictor Data" provides more information on the algorithms and underlying assumptions for
#' simulating placebo-group and BIP data.
#' The exported data sets include full rectangular data to allow the user to consider various biomarker sub-sampling designs. To generate and export such data,
#' \code{saveDataDir}, \code{nCasesPla}, and \code{nControlsPla} must be specified. In addition, if the biomarker is trichotomous and Approach 1 is used,
#' \code{sensBIP}, \code{specBIP}, \code{FP0BIP}, \code{FN2BIP}, \code{P0BIP}, and \code{P2BIP} must be specified;
#' if the biomarker is trichotomous and Approach 2 is used, or if it is continuous, \code{corr} must be specified.
#' Arguments pertaining to the simulated data export may also be vectorized in the following groups:
#'    \itemize{
#'     \item \code{nCasesPla}, \code{nControlsPla}
#'     \item \code{sensBIP}, \code{specBIP}, \code{FP0BIP}, \code{FN2BIP}
#'     \item \code{P0BIP}, \code{P2BIP}
#'     \item \code{corr}
#'   }
#' \code{nCasesPla} and \code{nControlsPla} must be vectorized together with \code{nCasesTx}, \code{nControlsTx}, and \code{nCasesTxWithS},
#' but the other three groups can be vectorized independently of the above varying arguments defining the power calculation scenarios in the active treatment group.
#' Only one of the four groups may be vectors at a time, otherwise an error will be generated. Each component of these vectors will result in the generation of
#' a separate BIP variable, in the same order, in the output data.
#'
#' The default file name for the outputted data sets is \code{fullData.RData}. A different file name may be specified by \code{saveDataFile}
#' as a single character string, to which the value of the "varying argument" for the power calculations will be appended for descriptive file naming purposes,
#' or, alternatively, a character vector may be specified with full file names (a single file will be produced for each value of the varying argument(s)).
#' Note: if the "varying argument" is \code{controlCaseRatio} or \code{p}, only one file will be generated because these arguments do not affect
#' the simulation of the full data; therefore, \code{saveDataFile} must be a character string in these cases.
#'
#' @return If \code{saveDir} is specified, an output list (named \code{pwr}) for each power scenario is saved as an \code{.RData} file. Otherwise, the function returns a list of lists,
#' where the outer list ranges over specified values of the varying argument(s) whose components denote distinct scenarios, and the inner list is the output list for each power scenario.
#' For a dichotomous or trichotomous biomarker, each output list has the following components:
#' \itemize{
#'   \item \code{power}: a numeric vector of fractions of simulated trials in which the null hypothesis \eqn{H_0} is rejected. Each value of the vector corresponds to a value in the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{RRt}: a numeric vector of correlate-of-risk relative-risk effect sizes. Each value of the vector corresponds to a value in the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{risk1_2}: a numeric vector of conditional endpoint risks given a high biomarker response in the active treatment group. Each value of the vector corresponds to a value in the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{risk1_0}: a numeric vector of conditional endpoint risks given a low biomarker response in the active treatment group. Each value of the vector corresponds to a value in the grid of treatment (vaccine) efficacies specified by \code{VElat0} and \code{VElat1}.
#'   \item \code{VElat2}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{VElat0}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{Plat2}: a numeric value specifying the prevalence of the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{Plat0}: a numeric value specifying the prevalence of the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#'   \item \code{P2}: a numeric value specifying the probability of high biomarker response for a dichotomous or trichotomous biomarker
#'   \item \code{P0}: a numeric value specifying the probability of low biomarker response for a dichotomous or trichotomous biomarker
#'   \item \code{alphaLat}: a numeric vector of the log odds of the clinical endpoint in the subgroup of active treatment recipients with the latent \eqn{x^{\ast}=0} (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{betaLat}: a numeric vector of the log odds ratio of the clinical endpoint comparing two subgroups of active treatment recipients differing in the latent \eqn{x^{\ast}} by 1 (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{sens}: a numeric vector of sensitivities (i.e., the probability of high biomarker response conditional on membership in the higher protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{spec}: a numeric vector of specificities (i.e., the probability of low biomarker response conditional on membership in the lower protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{FP0}: a numeric vector of false positive rates (i.e., the probability of high biomarker response conditional on membership in the lower protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{FN2}: a numeric vector of false negative rates (i.e., the probability of low biomarker response conditional on membership in the higher protected subgroup) of the observed dichotomous or trichotomous biomarker as a function of \code{rho}
#'   \item \code{NcompleteTx}: an integer value specifying \code{nCasesTx} + \code{nControlsTx}, i.e., the number, observed or projected, of active treatment recipients at risk at \eqn{\tau} with an observed endpoint or a completed follow-up through \eqn{\tau_{max}}
#'   \item \code{nCasesTx}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group
#'   \item \code{nCasesTxWithS}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response
#'   \item \code{controlCaseRatio}: an integer specifying the number of controls sampled per case for
#'   biomarker measurement in the without replacement case-control sampling design
#'   \item \code{VEoverall}: a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{risk0}: a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{alpha}: a numeric value specifying the two-sided Wald test type-I error rate
#'   \item \code{rho}: a numeric vector specifying distinct protection-relevant fractions of the variance of the observed biomarker
#' }
#'
#' For a continuous biomarker, each output list has the following components:
#' \itemize{
#'   \item \code{power}: a numeric vector of fractions of simulated trials in which the null hypothesis \eqn{H_0} is rejected. Rows represent calculations for different values of \code{rho} or \code{nCasesTx}, depending on which is a vector. Columns represent calculations for the grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup, specified by \code{VElowest}.
#'   \item \code{RRc}: a numeric vector of correlate-or-risk relative-risk effect sizes as a function of the grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup, specified by \code{VElowest}
#'   \item \code{betaLat}: a numeric vector specifying the log odds ratio of the clinical endpoint comparing two subgroups of active treatment recipients differing in the latent \eqn{x^{\ast}} by 1 (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{alphaLat}: a numeric vector specifying the the log odds of the clinical endpoint in the subgroup of active treatment recipients with the latent \eqn{x^{\ast}=0} (this coefficient estimate applies to a continuous biomarker)
#'   \item \code{PlatVElowest}: a numeric value specifying the prevalence of the latent lowest-efficacy subgroup for a continuous biomarker
#'   \item \code{VElowest}: a numeric vector specifying a grid of treatment (vaccine) efficacy levels in the latent lowest-efficacy subgroup for a continuous biomarker
#'   \item \code{sigma2obs}: a numeric value specifying the variance of the observed continuous biomarker or of the dichotomous or trichotomous biomarker simulated using 'approach 2'
#'   \item \code{NcompleteTx}: an integer value specifying \code{nCasesTx} + \code{nControlsTx}, i.e., the number, observed or projected, of active treatment recipients at risk at \eqn{\tau} with an observed endpoint or a completed follow-up through \eqn{\tau_{max}}
#'   \item \code{nCasesTx}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group
#'   \item \code{nCasesTxWithS}: an integer value specifying the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response
#'   \item \code{controlCaseRatio}: an integer value specifying the number of controls sampled per case for biomarker measurement in the without replacement case-control sampling design
#'   \item \code{VEoverall}: a numeric value specifying the overall treatment (vaccine) efficacy between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{risk0}: a numeric value specifying the overall placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#'   \item \code{alpha}: a numeric value specifying the two-sided Wald test type-I error rate
#'   \item \code{rho}: a numeric vector specifying distinct protection-relevant fractions of the variance of the observed biomarker
#' }
#'
#' If \code{saveDataDir} is specified, the simulated data, including placebo group and BIP data, are saved in one or more \code{.RData} file(s)
#' containing a list of lists of data frames.
#' The components of the outer list consist each of one Monte-Carlo iteration of simulated data for all values of \code{VElat0} or \code{VElat1} if
#' the biomarker is trichotomous, or of \code{VElowest} if the biomarker is continuous. Each data frame corresponds to one simulated trial.
#'
#' @examples
#'
#'## Trichotomous biomarker, Approach 1, varying sens and spec ##
#'## Specify sens, spec, FP0, FN2
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall,
#'              risk0=risk0, VElat0=VElat0, VElat1=VElat1, Plat0=Plat0,
#'              Plat2=Plat2, P0=P0, P2=P2, M=M, alpha=alpha, spec=spec,
#'              FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' \dontrun{
#' ## Trichotomous biomarker, Approach 2, varying rho ##
#' ## Saving simulated data (including placebo and BIP data)
#' ## Specify rho, sigma2obs, saveDataDir, saveDataFile, corr
#'
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' saveDataDir <- "~/myDir"
#' saveDataFile <- "myDataFile.RData"
#' corr <- 0.7
#' computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              VElat0=VElat0, VElat1=VElat1, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType,
#'              saveDataDir=saveDataDir, saveDataFile=saveDataFile, corr=corr)
#'
#'
#' ## dichotomous biomarker, Approach 2, varying rho ##
#' ## Plat0 + Plat2 = 1
#'
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              VElat0=VElat0, VElat1=VElat1, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#'
#' ## Continuous biomarker, varying rho ##
#'
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall, risk0=risk0,
#'              PlatVElowest=PlatVElowest, VElowest=VElowest, M=M, alpha=alpha,
#'              sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#'
#' ## Continuous biomarker, case-cohort sampling design, varying p ##
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' p <- c(0.01, 0.02, 0.03)
#' computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'              VEoverall=VEoverall, risk0=risk0, PlatVElowest=PlatVElowest,
#'              VElowest=VElowest, M=M, alpha=alpha, sigma2obs=sigma2obs,
#'              rho=rho, biomType=biomType, cohort=cohort, p=p)
#'
#' ## Continuous biomarker, saving output, varying sample sizes ##
#'
#' nCasesTx <- 32
#' nControlsTx <- 1000
#' nCasesTxWithS <- 32
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
#' saveFile <- "MyFile.RData"
#' computePower(nCasesTx=nCasesTx, nCasesTxWithS=nCasesTxWithS, nControlsTx=nControlsTx,
#'              controlCaseRatio=controlCaseRatio, VEoverall=VEoverall,
#'              risk0=risk0, PlatVElowest=PlatVElowest, VElowest=VElowest,
#'              M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho,
#'              biomType=biomType, saveDir=saveDir, saveFile=saveFile)
#' }
#'
#' @seealso \code{\link{computeN}}, \code{\link{plotPowerTri}}, \code{\link{plotPowerCont}}
#'
#' @import survival
#' @import osDesign
#' @importFrom stats dnorm fisher.test integrate pexp pnorm qnorm rbinom rmultinom rnorm uniroot
#'
#' @export
computePower <- function(nCasesTx, nControlsTx, nCasesTxWithS,
                         controlCaseRatio,
                         VEoverall, risk0,
                         VElat0=seq(0, VEoverall, len=20), VElat1=rep(VEoverall, 20),
                         VElowest=NULL,
                         Plat0, Plat2,
                         P0=Plat0, P2=Plat2,
                         PlatVElowest=NULL,
                         sens=NULL, spec=NULL, FP0=NULL, FN2=NULL,
                         M=100,
                         alpha=0.05,
                         sigma2obs=1, rho=1,
                         biomType=c("continuous", "trichotomous", "dichotomous"),
                         cohort=FALSE, p=NULL,
                         tpsMethod=c("PL", "ML","WL"),
                         saveDir=NULL, saveFile="CoRpower.RData",
                         saveDataDir=NULL, saveDataFile="fullData.RData",
                         corr=NULL, nCasesPla=NULL, nControlsPla=NULL,
                         sensBIP=NULL, specBIP=NULL, FP0BIP=NULL, FN2BIP=NULL,
                         P0BIP=P0, P2BIP=P2) {


  tpsMethod <- match.arg(tpsMethod, choices = c("PL","ML","WL"))
  biomType <- match.arg(biomType, choices = c("continuous", "trichotomous", "dichotomous"))

  # check if Approach 1 or Approach 2 is being used for trichotomous or dichotomous biomarker
  Approach2 <- (all(is.null(spec), is.null(sens), is.null(FP0), is.null(FN2)))

  # check sample size parameters are valid
  sampleSizes <- list("nCasesTx" = nCasesTx, "nCasesTxWithS" = nCasesTxWithS, "nControlsTx" = nControlsTx)
  if(!is.null(saveDataDir)){
    sampleSizes$nCasesPla <- nCasesPla
    sampleSizes$nControlsPla <- nControlsPla
  }

  # Find varying parameter
  Platx <- list("Plat0" = Plat0, "Plat2" = Plat2, "P0" = P0, "P2" = P2)
  sensSpec <- list("sens" = sens, "spec" = spec, "FP0" = FP0, "FN2" = FN2)
  vary <- getVaryingParam(sampleSizes, Platx, sensSpec, controlCaseRatio, p, rho, PlatVElowest)
  varyingParam <- vary$varyingParam
  varyingParamName <- vary$varyingParamName

  # logical value checking if multiple simulations of full data are needed
  # will switch to false later if controlCaseRatio or p is the varying parameter
  varyFullData <- TRUE

  # If full data (X, Y, S1, Z, and a BIP for treatment and placebo) is to be outputted, check for errors and input violations
  if(!is.null(saveDataDir)) {
    # check parameters are valid for saving full data as an output
    checkSaveDataParams(Approach2, biomType, corr, nCasesPla, nControlsPla, sensBIP, specBIP, FP0BIP, FN2BIP, P0BIP, P2BIP)
  }

  # Check file name parameters
  if (saveFile[[1]] != "CoRpower" | saveDataFile[[1]] != "fullData") {
    checkFileSavingParams(saveFile, saveDataFile, varyingParam, varyingParamName)
  }

  pwrAll <- list()

  for(i in 1:length(varyingParam[[1]])) {

    # If full data (X, Y, S1, Z, and a BIP for treatment and placebo) is to be outputted, initialize output list
    if(!is.null(saveDataDir) & varyFullData) {
      fullData <- list()
    }

    if ("nCasesTx" %in% varyingParamName) {
      nCasesTx <-  varyingParam$nCasesTx[i]
      nControlsTx <-varyingParam$nControlsTx[i]
      nCasesTxWithS <- varyingParam$nCasesTxWithS[i]
      if(!is.null(saveDataDir)) {
        nCasesPla <- varyingParam$nCasesPla[i]
        nControlsPla <- varyingParam$nControlsPla[i]
      }
    } else if ("Plat0" %in% varyingParamName) {
      Plat0 <- varyingParam$Plat0[i]
      Plat2 <- varyingParam$Plat2[i]
      P0 <- varyingParam$P0[i]
      P2 <- varyingParam$P2[i]
    } else if ("sens" %in% varyingParamName) {
      sens <- varyingParam$sens[i]
      spec <- varyingParam$spec[i]
      FP0 <- varyingParam$FP0[i]
      FN2 <- varyingParam$FN2[i]
    } else if (varyingParamName == "controlCaseRatio") {
      controlCaseRatio <- varyingParam$controlCaseRatio[i]
    } else if (varyingParamName == "p") {
      p <- varyingParam$p[i]
    } else if (varyingParamName == "rho") {
      rho <- varyingParam$rho[i]
    } else if (varyingParamName == "PlatVElowest") {
      PlatVElowest <- varyingParam$PlatVElowest[i]
    }

    # check sampling design input parameters are specified and valid
    checkSamplingDesign(cohort, p, controlCaseRatio)
    # check biomarker type and input parameters match
    checkBiomarkerType(biomType, P0, P2, VElowest, PlatVElowest)

    # Overall number in the treatment group observed to be at risk when the immune response is measured and that do not drop out (smaller than N):
    NcompleteTx <- nCasesTx + nControlsTx

    # Calculate NcompletePla: number in the placebo group observed to be at risk when the immune response is measured and that do not drop out
    if(!is.null(saveDataDir)) {
      NcompletePla <- nCasesPla + nControlsPla
    }

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

      # initialize power calculation vector
      powerstrinary <- numeric(length(VElat0))

      # Approach 2 in the manuscript (default choice):
      if (Approach2) {

        # Compute sens, spec, FP0, FP1, FN2, FN1
        ans <- computeSensSpecFPFN(sigma2obs, rho, Plat0, Plat2, P0, P2)

        # extract values for sens, spec, FP0, FP1, FN2, FN1
        sens <- sapply(ans, function(x) x[[1,10]])
        spec <- sapply(ans, function(x) x[[1,11]])
        FP0 <- sapply(ans, function(x) x[[1,12]])
        FP1 <- sapply(ans, function(x) x[[1,13]])
        FN2 <- sapply(ans, function(x) x[[1,14]])
        FN1 <- sapply(ans, function(x) x[[1,15]])

        # dataframe of rho, sens, spec, etc.
        # used to create Table 1: mapping of sigma2obs and rho to the sens, spec, etc. parameters
        table1 <- as.data.frame(round(cbind(rho, Plat0, P0, Plat2, P2, sens, spec, FP0, FN2, FP1, FN1),3))

        if (!is.null(saveDataDir) & varyFullData) {
          sigma2d <- (sigma2obs / corr^2) - sigma2obs

          sigma2BIP <- sigma2obs + sigma2d
          rhoBIP <- sigma2obs / (sigma2obs + sigma2d)

          # Compute sens, spec, FP0, FP1, FN2, FN1 for the BIP
          ansBIP <- computeSensSpecFPFN(sigma2BIP, rhoBIP, P0, P2, P0BIP, P2BIP)

          # extract values for sens, spec, FP0, FP1, FN2, FN1
          sensBIP <- sapply(ansBIP, function(x) x[[1,10]])
          specBIP <- sapply(ansBIP, function(x) x[[1,11]])
          FP0BIP <- sapply(ansBIP, function(x) x[[1,12]])
          FP1BIP <- sapply(ansBIP, function(x) x[[1,13]])
          FN2BIP <- sapply(ansBIP, function(x) x[[1,14]])
          FN1BIP <- sapply(ansBIP, function(x) x[[1,15]])
        }

      } else { # Approach 1 in the manuscript: use given sens, spec, FP0, and FN2 params

        # check lengths of sens, spec, FP0, and FN2 vectors are equal
        checkParamLengthsMatch(sens,spec,FP0,FN2)

        if (biomType=="dichotomous") {  # if dichotomous biomarker, FN1 and FP1 are irrelevant
          FN1 <- 0
          FP1 <- 0
          if ((P0 != (spec*Plat0 + FN2*Plat2)) | (P2 != (sens*Plat2 + FP0*Plat0))) {
            stop("Approach 1 was used for a dichotomous biomarker and the parameters Plat0, Plat2, P0, P2, sens, spec, FP0, and FN2 are inconsistent with each other")
          }
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

        # If treatment and placebo data is to be saved and Approach 1 is used, calculate sigma2tr and sigma2e
        # to be used to generate a BIP
        if(!is.null(saveDataDir) & varyFullData) {

          # check lengths of sensBIP, specBIP, FP0BIP, and FN2BIP vectors are equal
          checkParamLengthsMatch(sensBIP, specBIP, FP0BIP, FN2BIP)

          if (biomType=="dichotomous") {  # if dichotomous biomarker, FN1 and FP1 are irrelevant
            FN1BIP <- 0
            FP1BIP <- 0
            if ((P0BIP != (specBIP*P0 + FN2BIP*P2)) | (P2BIP != (sensBIP*P2 + FP0BIP*P0))) {
              stop("Approach 1 was used for a dichotomous biomarker and the BIP-associated parameters P0, P2, P0BIP, P2BIP, sensBIP, specBIP, FP0BIP, and FN2BIP are inconsistent with each other")
            }
          } else {
            # Apply formula (7) in the manuscript
            FN1BIP <- (P0BIP - specBIP*P0 - FN2BIP*P2)/P1   #P0, Plat0, Plat2 given params
            # Apply formula (8) in the manuscript
            FP1BIP <- (P2BIP - sensBIP*P2 - FP0BIP*P0)/P1
          }

          # Check if there is an error in the ranges of values due to an out of bounds input parameter
          if (any(FN1BIP < 0 | FN1BIP > 1 | FP1BIP < 0 | FP1BIP > 1)){
            stop("Approach 1 was used and one of the parameters sensBIP, specBIP, FP0BIP, FN2BIP is out of range")
          }
        }
      }

      # compute risks for biomarker subgroups and for latent subgroups
      risks <- computeRisks(biomType, RRoverall, risk0, sens, spec, FP0, FN2, FP1, FN1,
                            Plat0, Plat1, Plat2, P0, P1, P2, RRlat0, RRlat1, RRlat2)
      risk1_0 <- risks$risk1_0
      risk1_1 <- risks$risk1_1
      risk1_2 <- risks$risk1_2
      risk1lat_0 <- risks$risk1lat_0
      risk1lat_1 <- risks$risk1lat_1
      risk1lat_2 <- risks$risk1lat_2

      for (j in 1:M) {

        # If full data is to be outputted, initialize list of data frames for a single trial iteration
        if(!is.null(saveDataDir) & varyFullData) {
          fullDataIter <- list()
        }

        for (k in 1:length(VElat0)) {

          data <- simTrich(risk1lat_0[k], risk1lat_1[k], risk1lat_2[k], Plat0, Plat1, Plat2,
                           sens, spec, FP0, FN2, FP1, FN1,
                           nCasesTx, nCasesTxWithS, NcompleteTx, nCasesPla, NcompletePla,
                           biomType, controlCaseRatio, p, cohort, tpsMethod, alpha, saveDataDir,
                           sensBIP, specBIP, FP0BIP, FN2BIP, FP1BIP, FN1BIP, varyFullData)
          powerstrinary[k] <- powerstrinary[k] + data$addPower

          if(!is.null(saveDataDir) & varyFullData) {
            fullDataIter[[k]] <- data$simData
          }
        }

        if(!is.null(saveDataDir) & varyFullData) {
          fullData[[j]] <- fullDataIter
        }
      }

      # write out alpha intercept as logit(Y=1|s=0) for trinary/dichotomous case
      alphaLat <- c(logit(risk1_0))
      # write out beta coefficient as the log odds ratio: logit(Y=1|S=2)-logit(Y=1|s=0) for trinary/dichotomous case
      betaLat <- c(logit(risk1_2)-logit(risk1_0))
      # CoR effect sizes
      RRt <- risk1_2/risk1_0

      power <- powerstrinary / M

      # output list for trichotomous/dichotomous biomarker
      pwr <- list("power" = power, "RRt"=RRt, "risk1_2"=risk1_2, "risk1_0"=risk1_0, "VElat2"=VElat2, "VElat0"=VElat0,
                  "Plat2"=Plat2, "Plat0"=Plat0, "P2"=P2, "P0"=P0, "alphaLat"=alphaLat, "betaLat"=betaLat,
                  "sens"=sens, "spec"=spec, "FP0"=FP0, "FN2"=FN2, "NcompleteTx"=NcompleteTx, "nCasesTx"=nCasesTx,
                  "nCasesTxWithS"=nCasesTxWithS, "controlCaseRatio"=controlCaseRatio, "VEoverall"=VEoverall,
                  "risk0"=risk0, "alpha"=alpha, "rho"=rho)

      pwrAll[[i]] <- pwr

    } else if (biomType=="continuous") {


      # initialize power calculation vector
      powerscont <- numeric(length(VElowest))

      #################################################
      # Computations for a continuous biomarker

      # Define the truebetas (betaLat) indexed by the user-specified vector VElowest.
      # VElowest: a vector of fixed values of VE(x) for the subgroup of subjects with lowest X^* values,
      # where this subgroup has prevalence PlatVElowest

      o <- length(VElowest)
      nu <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
      truebetas <- rep(NA,o)
      alphalatvect <- rep(NA,o)

      ###################################################
      # Power calculations repeated for M simulations

      for (j in 1:M) {

        # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
        # using the function risk1cont()
        # Step 4 for continuous biomarker in manuscript

        # If full data is to be outputted, initialize list of data frames for a single trial iteration
        if(!is.null(saveDataDir) & varyFullData) {
          fullDataIter <- list()
        }

        for (k in 1:o) {  # loop through each value of VElowest

          # find solutions alphalat and betalat by solving eqn (4) in Appendix B using functions kernel() and alphaLatEqn()
          risk1latnu <- (1-VElowest[k])*risk0

          alphalat <- uniroot(alphaLatEqn, lower=-10, upper=10, nu=nu, risk1latnu=risk1latnu, sigma2obs=sigma2obs,
                              VEoverall=VEoverall, PlatVElowest=PlatVElowest, risk0=risk0)$root
          alphalatvect[k] <- alphalat

          # Second solve for betalat:
          D <- risk1latnu
          beta <- (log(D/(1-D)) - alphalat)/nu
          truebetas[k] <- beta

          # These simulations condition on nCasesTx (i.e., number of infections in vaccine arm) and
          # also on the number of controls fixed at controlCaseRatio*nCasesTx

          data <- simCont(nCasesTx, NcompleteTx, nCasesTxWithS, nCasesPla, NcompletePla,
                          alphalat, beta, sigma2obs, sigma2tr, sigma2e, nu, PlatVElowest, VElowest[k],
                          risk0, controlCaseRatio, p, cohort, tpsMethod, alpha, corr, saveDataDir, varyFullData)
          powerscont[k] <- powerscont[k] + data$addPower

          if(!is.null(saveDataDir) & varyFullData) {
            fullDataIter[[k]] <- data$simData
          }
        }

        if(!is.null(saveDataDir) & varyFullData) {
          fullData[[j]] <- fullDataIter
        }
      }

      # RRc the relative risks that are the effect sizes RR_c that need to be on the x-axis of powerplots
      RRc <- exp(truebetas)

      power <- powerscont / M

      # output list for continuous biomarker
      pwr <- list("power"=power, "RRc"=RRc, "betaLat"=truebetas, "alphaLat"=alphalatvect, "PlatVElowest"=PlatVElowest,
                  "VElowest"=VElowest, "sigma2obs"=sigma2obs, "NcompleteTx"=NcompleteTx, "nCasesTx"=nCasesTx,
                  "nCasesTxWithS"=nCasesTxWithS, "controlCaseRatio"=controlCaseRatio, "VEoverall"=VEoverall,
                  "risk0"=risk0, "alpha"=alpha, "rho"=rho)

      pwrAll[[i]] <- pwr

    }



    fileName <- ""
    if ("nCasesTx" %in% varyingParamName) {
      name <- varyingParamName[-2]
      param <- varyingParam[-2]
      paramValues <- numeric()
      for(j in 1:length(param)) {
        paramValues <- c(paramValues, param[[j]][i])
      }
      fileName <- paste0("_", paste0(name,"_",paramValues, collapse="_"))
    } else if (is.list(varyingParam)) {
      paramValues <- numeric()
      for(j in 1:length(varyingParam)) {
        paramValues <- c(paramValues, varyingParam[[j]][i])
      }
      fileName <- paste0("_", paste0(varyingParamName,"_",paramValues, collapse="_"))
    }

    # If saveDir is specified, save output list to .RData file with given location.
    # If saveFile is specified, use given file name; otherwise, file name will be "CoRpower"
    if(!is.null(saveDir)) {
      if (saveFile[[1]] != "CoRpower") {
        if(length(saveFile)==1) {
          save(pwr, file=paste0(file.path(saveDir, paste0(substr(saveFile, start = 1, stop = nchar(saveFile) - 6),fileName)), ".RData"))
        } else {
          save(pwr, file=paste0(file.path(saveDir, saveFile[i])))
        }
      } else {
        save(pwr, file=paste0(file.path(saveDir, paste0("CoRpower", fileName)), ".RData"))
      }
    }

    # If saveDataDir is specified, save output list to .RData file with given file location
    # If saveDataFile is specified, use given file name; otherwise, file name will be "fullData"
    if(!is.null(saveDataDir) & varyFullData) {
      if (varyingParamName == "controlCaseRatio" | varyingParamName == "p") {
        fileName <- ""
      }
      if (saveDataFile[[1]] != "fullData.RData") {
        if(length(saveDataFile)==1) {
          save(fullData, file=paste0(file.path(saveDataDir, paste0(substr(saveDataFile, 1, nchar(saveDataFile) - 6),fileName)), ".RData"))
        } else {
          save(fullData, file=paste0(file.path(saveDataDir, saveDataFile[i])))
        }
      } else {
        save(fullData, file=paste0(file.path(saveDataDir, paste0("fullData",fileName)),".RData"))
      }
      rm(fullData)
    }

    if (varyingParamName == "controlCaseRatio" | varyingParamName == "p") {
      varyFullData <- FALSE
    }
  }

  return(pwrAll)
}
