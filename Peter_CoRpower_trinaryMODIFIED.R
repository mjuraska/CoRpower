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
#' @param Spec For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the sensitivity of the measured biomarker. Specifying \code{Spec = rep(1,4)}, \code{FP1 = rep(0,4)}, \code{Sens = rep(1,4)}, and \code{FN1 = rep(0,4)} indicates that approach 2 is used.
#' @param FP1 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false positive rate (FP^1) of the measured/observed biomarker. Specifying \code{Spec = rep(1,4)}, \code{FP1 = rep(0,4)}, \code{Sens = rep(1,4)}, and \code{FN1 = rep(0,4)}  indicates that approach 2 is used.
#' @param Sens For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the specificity of the measured/observed biomarker. Specifying \code{Spec = rep(1,4)}, \code{FP1 = rep(0,4)}, \code{Sens = rep(1,4)}, and \code{FN1 = rep(0,4)}  indicates that approach 2 is used.
#' @param FN1 For a trichotomous biomarker (or binary as a special case), simulated using 'approach 1', a vector of length 4 with values for the first false negative rate (FN^1) of the measured/observed biomarker. Specifying \code{Spec = rep(1,4)}, \code{FP1 = rep(0,4)}, \code{Sens = rep(1,4)}, and \code{FN1 = rep(0,4)}  indicates that approach 2 is used.
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
#' Approach 1 specifies \code{Spec}, \code{Sens}, \code{FP1}, and \code{FN1} which determine
#' \code{FP2} and \code{FN2} from equations (7) and (8).  Four values are required for each input parameter, to allow the evaluation of biomarkers with different levels of measurement error.
#'
#' Approach 2 for a trichotomous (or binary) biomarker specifies \code{sigma2obs} and \code{rho}, again requiring four settings
#' for \code{rho}; this approach assumes the normal measurement error model (4) in the manuscript.  Specifying \code{Spec = rep(1,4)}, \code{FP1 = rep(0,4)}, \code{Sens = rep(1,4)}, and \code{FN1 = rep(0,4)}
#' defaults to approach 2, which is what is used in illustrations in the manuscript.
#'
#' For a continuous biomarker, \code{VElowestvect}, \code{sigma2obs} and \code{rho} must be specified.  Setting \code{VElowestvect = NULL}
#' indicates that the biomarker is not continuous (it is trichotomous).
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
#' ## Spec[1]=1, FP1[1]=0, Sens[1]=1, FN1[1]=0
#' ## Note the other elements of at least one of these four parameter vectors need
#' ## to be set to a different value to tell the program to use Approach 1
#' Spec <- c(1, 0.9, 0.8, 0.7)
#' FP1 <- rep(0,4)
#' Sens <- c(1, 0.9, 0.8, 0.7)
#' FN1 <- rep(0,4)
#' RRlat1 <- rep(0,100) # will be turned into NA in computepower() for binary case
#' RRlat0 <- seq(1,RRoverall,len=100) # 100 data points for the power curve
#'
#' ## Parameters used for for continuous biomarker calculations, or
#' ## for trichotomous/binary biomarkers under Approach 2
#' ## Note these are immaterial trichotomous/binary Approach 1
#' sigma2obs <- 1
#' rhos <- c(1,0.9,0.7,0.5) # rho = 1 corresponds to no measurement error case
#' PlatVElowest <- 0.40
#' VElowestvect <- seq (0, 1-RRoverall,len=100)
#'
#' ## Binary biomarker, Approach 1
#' ## i.e. Plat1=P1=0 and Plat0 = 1- Plat2
#' Plat0 <- 0.1
#' Plat2 <- 1-Plat0
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL,
#' controlCaseRatio=5, M=1000, alpha=0.05, sigma2obs, rhos, Spec, FP1, Sens, FN1)
#'
#'## Trichotomous biomarker, Approach 1
#' Plat0 <- 0.1
#' Plat2 <- 0.4
#' P2 <- Plat2 # different values of P2 can be set
#' P0 <- Plat0 # different values of P0 can be set
#' computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,
#' VElowestvect=NULL, controlCaseRatio=5, M=1000, alpha=0.05, sigma2obs, rhos, Spec, FP1, Sens, FN1)
#'
#' ## Binary biomarker, Approach 2
#' computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls, risk0,
#' RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL, controlCaseRatio=5,
#' M=1000, alpha=0.05, sigma2obs, rhos, Spec=rep(1,4), FP1=rep(0,4), Sens=rep(1,4), FN1=rep(0,4))
#'
#' ## Trichotomous biomarker, Approach 2
#' computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest=0,VElowestvect=NULL,
#' controlCaseRatio=5, M=1000, alpha=0.05, sigma2obs, rhos, Spec=rep(1,4), FP1=rep(0,4), Sens=rep(1,4), FN1=rep(0,4))
#'
#' ## Continuous biomarker

#' computepower(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
#' risk0, RRoverall, Plat0,Plat2, P0,P2, RRlat0,RRlat1, PlatVElowest,VElowestvect,
#' controlCaseRatio=5, M=1000, alpha=0.05, sigma2obs, rhos)
#'
#' @import survival
#' @import osDesign
#' @export
computepower <- function(numAtRiskTauCases, numAtRiskTauCasesPhase2, numAtRiskTauControls,
                         risk0, RRoverall,
                         Plat0, Plat2,
                         P0=0.3, P2=0.2,
                         RRlat0=seq(1,RRoverall,len=20), RRlat1=rep(RRoverall,20),
                         PlatVElowest, VElowestvect,
                         controlCaseRatio=5,
                         M=100,
                         alpha=0.05,
                         sigma2obs=1,
                         rhos=c(1,.9,.7,.5),
                         Spec=rep(1,4), FP1=rep(0,4), Sens=rep(1,4), FN1=rep(0,4)) {
  
  # sigma2tr is the variance of the true biomarker X
  # rhos must be a vector with 4 values and is for the continuous and binary biomarker correlates correlations
  # RRlat0 is the span of true relative risks (vaccine vs. placebo) in the latent lower protected subgroup
  # RRlat1 is the span of true relative risks (vaccine vs. placebo) in the latent medium protected subgroup
  # The power calculations should always include rho=1 in the first element of the rhos vector,
  # as the best case scenario, and the plotting functions assume this.
  
  # VElowestvect is used for a continuous biomarker- a vector of fixed value of VE(x) for
  # the subgroup of subjects with lowest X^* values, where this subgroup has prevalence PlatVElowest
  
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
    if (min(VElowestvect)==0 & PlatVElowest > 1 - VEoverall) {stop("Error: Input parameters PlatVElowest and VElowestvect violate probability constraints for normal marker calculations") }
  }
  
  checkprobabilityviolation(VEoverall,VElat0,Plat0,Plat2,PlatVElowest,VElowestvect)
  
  sigma2e <- (1-rhos)*sigma2obs
  
  sigma2tr <- rhos*sigma2obs
  
  #################################################
  # Computations for a trinary biomarker
  
  Approach2 <- Spec==rep(1,4) & FP1==rep(0,4) & Sens==rep(1,4) & FN1==rep(0,4)   #*# if App2, = TRUE TRUE TRUE TRUE
  Approach2 <- length(Approach2[Approach2])==4
  if (Approach2) {
    # Default choice
    
    # Compute Sens, Spec, FP1, FP2, FN1, FN2
    ans <- computeSensSpecFPFN(sigma2obs,rhos,Plat0,Plat2,P0,P2)
    Sens <- unlist(lapply(ans, function(x) x[[1,10]])) 
    Spec <- unlist(lapply(ans, function(x) x[[1,11]]))
    FP1 <- unlist(lapply(ans, function(x) x[[1,12]])) 
    FP2 <- unlist(lapply(ans, function(x) x[[1,13]])) 
    FN1 <- unlist(lapply(ans, function(x) x[[1,14]]))
    FN2 <- unlist(lapply(ans, function(x) x[[1,15]])) 
    
    # Write out a vector that can be used to make a table mapping (rho,sigma2obs) to the Sens etc. parameters
    
    x1 <- unlist(c(sigma2obs[1],rhos[1],Plat0,P0,Plat2,P2,Sens[1],Spec[1],FP1[1],FN1[1],FP2[1],FN2[1]))
    #cat("sigma2obs,rho,Plat0,P0,Plat2,P2,Sens,Spec,FP1,FN1,FP2,FN2 in order of the vector rhos")
    #cat("\n")
    #cat(paste(round(x1,3)),"\n")
    
    x2 <- unlist(c(sigma2obs[2],rhos[2],Plat0,P0,Plat2,P2,Sens[2],Spec[2],FP1[2],FN1[2],FP2[2],FN2[2]))
    #cat("\n")
    #cat(paste(round(x2,3)),"\n")
    
    x3 <- unlist(c(sigma2obs[3],rhos[3],Plat0,P0,Plat2,P2,Sens[3],Spec[3],FP1[3],FN1[3],FP2[3],FN2[3]))
    #cat("\n")
    #cat(paste(round(x3,3)),"\n")
    
    x4 <- unlist(c(sigma2obs[4],rhos[4],Plat0,P0,Plat2,P2,Sens[4],Spec[4],FP1[4],FN1[4],FP2[4],FN2[4]))
    #cat("\n")
    #cat(paste(round(x4,3)),"\n")

  }
  
  # Approach 1 in the manuscript:
  if (!Approach2) { #*# use given Sens, Spec, FP1, and FN1 params
    
    # Apply formula (8) in the manuscript
    FN2 <- (P0 - Spec*Plat0 - FN1*Plat2)/Plat1   #*#P0, Plat0, Plat2 given params
    
    # Apply formula (9) in the manuscript
    FP2 <- (P2 - Sens*Plat2 - FP1*Plat0)/Plat1
    
        # Sens1 <- Sens[1] #*# don't need these because we can treat Sens as a vector and it is given
        # Sens2 <- Sens[2]
        # Sens3 <- Sens[3]
        # Sens4 <- Sens[4]
        # Spec1 <- Spec[1]
        # Spec2 <- Spec[2]
        # Spec3 <- Spec[3]
        # Spec4 <- Spec[4]
        # FP11 <- FP1[1]
        # FP12 <- FP1[2]
        # FP13 <- FP1[3]
        # FP14 <- FP1[4]
        # FN11 <- FN1[1]
        # FN12 <- FN1[2]
        # FN13 <- FN1[3]
        # FN14 <- FN1[4]
        
        
        # # Apply formula (8) in the manuscript
        # FN21 <- (P0 - Spec1*Plat0 - FN11*Plat2)/Plat1 #*# don't need b/c FN2 is a vector
        # FN22 <- (P0 - Spec2*Plat0 - FN12*Plat2)/Plat1
        # FN23 <- (P0 - Spec3*Plat0 - FN13*Plat2)/Plat1
        # FN24 <- (P0 - Spec4*Plat0 - FN14*Plat2)/Plat1
        # 
        # # Apply formula (9) in the manuscript
        # FP21 <- (P2 - Sens1*Plat2 - FP11*Plat0)/Plat1
        # FP22 <- (P2 - Sens2*Plat2 - FP12*Plat0)/Plat1
        # FP23 <- (P2 - Sens3*Plat2 - FP13*Plat0)/Plat1
        # FP24 <- (P2 - Sens4*Plat2 - FP14*Plat0)/Plat1
    
    
    # Check if an error in the ranges of values due to an out of
    # bounds input parameter
    
    if (any(FN2 < 0 | FN2 > 1 | FP2 < 0 | FP2 > 1)){
      stop("Approach 1 was used and one of the parameters Sens, Spec, FP1, FN1 is out of range")
    }
        # if (FN21 < 0 | FN21 > 1 | FN22 < 0 | FN22 > 1 | FN23 < 0 | FN23 > 1 | FN24 < 0 | FN24 > 1 | FN21 < 0 | FP21 > 1 | FP22 < 0 | FP22 > 1 | FP23 < 0 | FP23 > 1 | FP24 < 0 | FP24 > 1) {
        #   cat(paste("Approach 1 was used and one of the parameters Sens, Spec, FP1, FN1 is out of range"),"\n")
        # }
    
  }
  
  # Binary biomarker special case (to remove small values of P1x)
  if ((P0+P2)==1) {
    P1 <- 0
    P2 <- 1 - P0
  }
  
  # Compute the marginal risks:
  # Made it to the end of follow-up HIV negative
  risk1 <- RRoverall*risk0
  
  # Observed risks P(Y(1)=1|S(1)=0, 1, or 2)  #*# for diff values of rho; using Bayes' rule
  
  probX0_cond_S2 <- FP1*Plat0/P2
  probX1_cond_S2 <- FP2*Plat1/P2
  probX2_cond_S2 <- Sens*Plat2/P2
  risk1_2 <- (probX0_cond_S2 %o% RRlat0 + probX1_cond_S2 %o% RRlat1 + probX2_cond_S2 %o% RRlat2 )*risk0 # use outer product to get matrix with nrow=length(rho), ncol=length(RRlat0)
  probX0_cond_S0 <- Spec*Plat0/P0
  probX1_cond_S0 <- FN2*Plat1/P0
  probX2_cond_S0 <- FN1*Plat2/P0
  risk1_0 <- (probX0_cond_S0 %o% RRlat0 + probX1_cond_S0 %o% RRlat1 + probX2_cond_S0 %o% RRlat2)*risk0
  risk1_1 <- (risk1 - risk1_0*P0 - risk1_2*P2)/P1
  
      # probX0_cond_S2 <- FP12*Plat0/P2 #*# don't need these anymore because the variables are vectors now
      # probX1_cond_S2 <- FP22*Plat1/P2
      # probX2_cond_S2 <- Sens2*Plat2/P2
      # risk1hi2 <- (RRlat0*probX0_cond_S2 + RRlat1*probX1_cond_S2 + RRlat2*probX2_cond_S2)*risk0
      # probX0_cond_S0 <- Spec2*Plat0/P0
      # probX1_cond_S0 <- FN22*Plat1/P0
      # probX2_cond_S0 <- FN12*Plat2/P0
      # risk1lo2 <- (RRlat0*probX0_cond_S0 + RRlat1*probX1_cond_S0 + RRlat2*probX2_cond_S0)*risk0
      # risk1med2 <- (risk1 - risk1lo2*P0 - risk1hi2*P2)/P1
      # 
  
  # Note: For the binary biomarker special case, the risk1medx are NA
  
  esvect <- risk1_2/risk1_0 # matrix with nrow=length(rho) and ncol=length(RRlat0)
  
      # esvect1 <- risk1hi1/risk1lo1  # the index effect size because rho = 1
      # esvect2 <- risk1hi2/risk1lo2
      # esvect3 <- risk1hi3/risk1lo3
      # esvect4 <- risk1hi4/risk1lo4
  
  # Vaccine risks within the latent subgroups (independent of rho of course)
  
  risk1lat_2 <- RRlat2*risk0
  risk1lat_1 <- RRlat1*risk0
  risk1lat_0 <- RRlat0*risk0
  
  
  #################################################
  # Computations for a continuous biomarker
  # Define the truebetas indexed by the user-specified vector VElowestvect

  o <- length(VElowestvect)
  if (o>0){
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
  
      # powerstrinaryrho1 <- rep(0,length(esvect1))
      # powerstrinaryrho2 <- rep(0,length(esvect1))
      # powerstrinaryrho3 <- rep(0,length(esvect1))
      # powerstrinaryrho4 <- rep(0,length(esvect1))
      # powerscontrho1  <- rep(0,o)
      # powerscontrho2  <- rep(0,o)
      # powerscontrho3  <- rep(0,o)
      # powerscontrho4  <- rep(0,o)
  
  for (i in 1:M) {
    
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
      
      # Given specifications for Spec, FP1, Sens, and FN1 and a logical value indicating if the 
      # biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2)
      biomarkerGroups <- function(SpecSens, binary){
        Spec <- SpecSens[1]
        Sens <- SpecSens[2]
        FP1 <- SpecSens[3]
        FN1 <- SpecSens[4]
        FP2 <- SpecSens[5]
        FN2 <- SpecSens[6]
        if(binary==TRUE){
          Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP1-Spec,FP1))),
                           rmultinom(N2,1,adjustprob(c(FN1,1-FN1-Sens,Sens))))
        } else{
          Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP1-Spec,FP1))),
                           rmultinom(N1,1,adjustprob(c(FN2,1-FP2-FN2,FP2))),
                           rmultinom(N2,1,adjustprob(c(FN1,1-FN1-Sens,Sens))))
        }
        rownames(Svalues) <- c("S=0","S=1","S=2")
        Svalues <- ifelse(Svalues[1,]==1,0,ifelse(Svalues[2,]==1,1,2))
        return(Svalues)
      }
      
      SpecSens <- cbind(Spec,Sens,FP1,FN1,FP2,FN2)
      
      if ((P0+P2)==1) { # binary case only

        S <- t(apply(SpecSens, 1, function(x) biomarkerGroups(x, binary=TRUE))) # each row is a set of Sens, Spec, etc. parameters
        
            # Srho1 <- cbind(rmultinom(N0,1,adjustprob(c(Spec1,1-FP11-Spec1,FP11))),
            #                rmultinom(N2,1,adjustprob(c(FN11,1-FN11-Sens1,Sens1))))
            # Srho1 <- ifelse(Srho1[1,]==1,0,ifelse(Srho1[2,]==1,1,2))
            # 
        
      } else { # trichotomous
        
        S <- t(apply(SpecSens, 1, function(x) biomarkerGroups(x, binary=FALSE))) # each row is a set of Sens, Spec, etc. parameters
        
            # Srho1 <- cbind(rmultinom(N0,1,adjustprob(c(Spec1,1-FP11-Spec1,FP11))),
            #                rmultinom(N1,1,adjustprob(c(FN21,1-FP21-FN21,FP21))),
            #                rmultinom(N2,1,adjustprob(c(FN11,1-FN11-Sens1,Sens1))))
            # Srho1 <- ifelse(Srho1[1,]==1,0,ifelse(Srho1[2,]==1,1,2))
            
      }
      
      # Now keep the S's in nPhase2 of the cases (deleting the rest) and in controlCaseRatio*nPhase2 controls
      casesinds <- c(1:N)[Y==1]
      keepcasesinds <- sample(casesinds,nPhase2,replace=FALSE)
      controlinds <- c(1:N)[Y==0]
      keepcontrolinds <- sample(controlinds,controlCaseRatio*nPhase2,replace=FALSE)
      keepinds <- sort(c(keepcasesinds,keepcontrolinds))
      
      # Those with data:
      Ycc <- Y[keepinds]
      Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(Sens)
          # Sccrho1 <- Srho1[keepinds]
          # Sccrho2 <- Srho2[keepinds]
          # Sccrho3 <- Srho3[keepinds]
          # Sccrho4 <- Srho4[keepinds]
      
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
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
          pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[k,j] <- powerstrinary[k,j] + 1}
        }
      }
      
      
          # lodim <- dim(table(Ycc,Sccrho1))[2]<2 #*# check there are at least two biomarker categories (columns)
          # zerosflag <-  lodim
          # if (dim(table(Ycc,Sccrho1))[2]==3) { #*# check if any categories have zero entries
          #   zerosflag <- table(Ycc,Sccrho1)[1,1]==0 | table(Ycc,Sccrho1)[1,2]==0 | table(Ycc,Sccrho1)[1,3]==0 | table(Ycc,Sccrho1)[2,1]==0 | table(Ycc,Sccrho1)[2,2]==0 | table(Ycc,Sccrho1)[2,3]==0 }
          # 
          # if (zerosflag) {
          #   if (lodim) { pval <- 1}
          #   if (!lodim) { #*# there are zeros, so Fisher's exact test is used
          #     pval <- fisher.test(table(Ycc,Sccrho1)[,c(1,dim(table(Ycc,Sccrho1))[2])])$p.value }
          #   if (pval <= alpha & length(Ycc[Sccrho1==2&Ycc==1])/length(Sccrho1[Sccrho1==2]) < length(Ycc[Sccrho1==0&Ycc==1])/length(Sccrho1[Sccrho1==0])) {
          #     powerstrinaryrho1[j] <- powerstrinaryrho1[j] + 1}}
          # 
          # if (!zerosflag) {
          #   fit <- tps(Ycc~Sccrho1,nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
          #   pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          #   if (pval <= alpha & fit$coef[2] < 0) { powerstrinaryrho1[j] <- powerstrinaryrho1[j] + 1}}
          # 
          # 
      
    }
    
    # Continuous biomarker
    
    # Simulate the infection indicators of all vaccine recipients, from a logistic regression model
    # using the function risk1cont() above
    if (o>0){ #*# at least one lowestVE specified (TRUE if continuous approach is to be used)
      
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

            # f <- function(x) {
            #   ans <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[1]))
            #   return(ans)
            # }
            # denomdensityXcases <- integrate(f,lower=nus[1],upper=5)$value
            # denomdensityXcases <- denomdensityXcases + PlatVElowest*(1-VElowestvect[j])*risk0
            # 
            # numerdensXcases <- function(x) {
            #   num <- risk1cont(x,alphalat,beta)*dnorm(x/sqrt(sigma2tr[1]))
            #   num[x <= nus[1]] <- PlatVElowest*(1-VElowestvect[j])*risk0
            #   return(num)
            # }
            # 
            # numerdensXcontrols <- function(x) {
            #   num <- (1-risk1cont(x,alphalat,beta))*dnorm(x/sqrt(sigma2tr[1]))
            #   num[x <= nus[1]] <- PlatVElowest*(1-(1-VElowestvect[j])*risk0)
            #   return(num)
            # }
            # 
            # Xpoints <- seq(-3.5,3.5,len=25000)
            # probscases <-    numerdensXcases(Xpoints)/denomdensityXcases
            # probscontrols <- numerdensXcontrols(Xpoints)/(1-denomdensityXcases)
            # 
            # Xcasesrho1 <-    sample(Xpoints,size=nCases,prob=probscases,replace=TRUE)
            # Xcontrolsrho1 <- sample(Xpoints,size=N-nCases,prob=probscontrols,replace=TRUE)
            # Xrho1 <- c(Xcasesrho1,Xcontrolsrho1)
        
        
        # Create the 4 immune response variables for the 4 degrees of measurement error
        error <- t(sapply(sigma2e, function(x) rnorm(N,mean=0,sd=sqrt(x))))
        S <- X + error
            # Srho1 <- Xrho1 + rnorm(N,mean=0,sd=sqrt(sigma2erho1))
            # Srho2 <- Xrho2 + rnorm(N,mean=0,sd=sqrt(sigma2erho2))
            # Srho3 <- Xrho3 + rnorm(N,mean=0,sd=sqrt(sigma2erho3))
            # Srho4 <- Xrho4 + rnorm(N,mean=0,sd=sqrt(sigma2erho4))
        
        # Now keep the S's in nPhase2 of the cases (deleting the rest) and in controlCaseRatio*nPhase2 of the controls
        casesinds <- c(1:N)[Y==1]
        keepcasesinds <- sample(casesinds,nPhase2,replace=FALSE)
        controlinds <- c(1:N)[Y==0]
        keepcontrolinds <- sample(controlinds,controlCaseRatio*nPhase2,replace=FALSE)
        keepinds <- sort(c(keepcasesinds,keepcontrolinds))
        
        # Those with data:
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rhos)
            # Sccrho1 <- Srho1[keepinds]
            # Sccrho2 <- Srho2[keepinds]
            # Sccrho3 <- Srho3[keepinds]
            # Sccrho4 <- Srho4[keepinds]

        for(k in 1:nrow(Scc)){
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
          pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          if (pval <= alpha & fit$coef[2] < 0) { powerscont[k,j] <- powerscont[k,j] + 1}
        }

            # fit <- tps(Ycc~Sccrho1,nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
            # pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
            # if (pval <= alpha & fit$coef[2] < 0) { powerscontrho1[j] <- powerscontrho1[j] + 1}
      }
    }
  }
  
  powerstrinary <- powerstrinary/M
      # powerstrinaryrho2 <- powerstrinaryrho2/M
      # powerstrinaryrho3 <- powerstrinaryrho3/M
      # powerstrinaryrho4 <- powerstrinaryrho4/M
  powerscont <- powerscont/M
      # powerscontrho2 <- powerscontrho2/M
      # powerscontrho3 <- powerscontrho3/M
      # powerscontrho4 <- powerscontrho4/M
  
      # anstrin <- cbind(powerstrinaryrho1,powerstrinaryrho2,powerstrinaryrho3,powerstrinaryrho4)
      # anscont <- cbind(powerscontrho1,powerscontrho2,powerscontrho3,powerscontrho4)
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
# The first approach specifies Spec, Sens, FP1, FN1 which determine FP2 and FN2 from equations (8) and (9).  4 settings
# must be specified.
#
# The second approach specifies sigma2obs and rho, again requiring 4 settings for rho; this approach assumes the normal
# measurement error model (4) in the article.  Specifying
# the vectors Spec=1, Sens=1, FP1=1, FN1=1 defaults to approach 2, which is used in the manuscript.

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
                           Spec=rep(1,4), FP1=rep(0,4), Sens=rep(1,4), FN1=rep(0,4)) {
  
  # Like computepower except RRlat2 (now RRlat2point) and VElowest are singulars, and
  # noobsatriskmotaucasesvect and noobsatriskmotaucontrolsvect are vectors
  # reflecting different sample sizes
  
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
  
  ###############################################################
  # Computations for a trinary biomarker
  
  Approach2 <- Spec==rep(1,4) & FP1==rep(0,4) & Sens==rep(1,4) & FN1==rep(0,4)
  Approach2 <- length(Approach2[Approach2])==4
  
  if (Approach2) {
    # Default choice
    
    # Compute Sens, Spec, FP1, FP2, FN1, FN2
    ans <- computeSensSpecFPFN(sigma2obs,rhos,Plat0,Plat2,P0,P2)
    Sens <- unlist(lapply(ans, function(x) x[[1,10]])) 
    Spec <- unlist(lapply(ans, function(x) x[[1,11]]))
    FP1 <- unlist(lapply(ans, function(x) x[[1,12]])) 
    FP2 <- unlist(lapply(ans, function(x) x[[1,13]])) 
    FN1 <- unlist(lapply(ans, function(x) x[[1,14]]))
    FN2 <- unlist(lapply(ans, function(x) x[[1,15]]))
  }
  
  # Approach 1 in the manuscript:
  if (!Approach2) {
    
    # Apply formula (8) in the manuscript
    FN2 <- (P0 - Spec*Plat0 - FN1*Plat2)/Plat1   #*#P0, Plat0, Plat2 given params
    
    # Apply formula (9) in the manuscript
    FP2 <- (P2 - Sens*Plat2 - FP1*Plat0)/Plat1
    
    #Check if an error in the ranges of values due to an out of
    # bounds input parameter
    if (any(FN2 < 0 | FN2 > 1 | FP2 < 0 | FP2 > 1)){
      stop("Approach 1 was used and one of the parameters Sens, Spec, FP1, FN1 is out of range")
    }
  }
  
  # Binary biomarker special case (to remove small values of P1x)
  if ((P0+P2)==1) {
    P1 <- 0
    P2 <- 1 - P0
  }
  
  ################## OLD on the chopping block
  ## Compute the marginal risks:
  ## Made it to the end of follow-up HIV negative
  #risk1 <- RRoverall*risk0
  #
  ## Observed risks P(Y(1)=1|S(1)=lo, med, or hi)
  #risk1hi1 <- (RRlat0point*FP11 + RRlat1point*FP21 + RRlat2point*Sens1)*risk0
  #risk1lo1 <- (RRlat0point*Spec1 + RRlat1point*FN21 + RRlat2point*FN11)*risk0
  #risk1hi2 <- (RRlat0point*FP12 + RRlat1point*FP22 + RRlat2point*Sens2)*risk0
  #risk1lo2 <- (RRlat0point*Spec2 + RRlat1point*FN22 + RRlat2point*FN12)*risk0
  #risk1hi3 <- (RRlat0point*FP13 + RRlat1point*FP23 + RRlat2point*Sens3)*risk0
  #risk1lo3 <- (RRlat0point*Spec3 + RRlat1point*FN23 + RRlat2point*FN13)*risk0
  #risk1hi4 <- (RRlat0point*FP14 + RRlat1point*FP24 + RRlat2point*Sens4)*risk0
  #risk1lo4 <- (RRlat0point*Spec4 + RRlat1point*FN24 + RRlat2point*FN14)*risk0
  
  # Compute the marginal risks:
  # Made it to the end of follow-up HIV negative
  risk1 <- RRoverall*risk0
  
  # Observed risks P(Y(1)=1|S(1)=0, 1, or 2)  #*# for diff values of rho; using Bayes' rule
  
  probX0_cond_S2 <- FP1*Plat0/P2
  probX1_cond_S2 <- FP2*Plat1/P2
  probX2_cond_S2 <- Sens*Plat2/P2
  risk1_2 <- (probX0_cond_S2 * RRlat0point + probX1_cond_S2 * RRlat1point + probX2_cond_S2 * RRlat2point)*risk0 #vector with length=length(rhos)
  probX0_cond_S0 <- Spec*Plat0/P0
  probX1_cond_S0 <- FN2*Plat1/P0
  probX2_cond_S0 <- FN1*Plat2/P0
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
  
  #################################################
  # Computations for a continuous biomarker
  # Define the truebeta indexed by the user-specified vector VElowest
  
  if (!is.null(VElowest)){
    
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
      
      # Given specifications for Spec, FP1, Sens, and FN1 and a logical value indicating if the 
      # biomarker is binary or not, the function returns a vector composed of biomarker levels (S=0,1,2)
      biomarkerGroups <- function(SpecSens, binary){
        Spec <- SpecSens[1]
        Sens <- SpecSens[2]
        FP1 <- SpecSens[3]
        FN1 <- SpecSens[4]
        FP2 <- SpecSens[5]
        FN2 <- SpecSens[6]
        if(binary==TRUE){
          Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP1-Spec,FP1))),
                           rmultinom(N2,1,adjustprob(c(FN1,1-FN1-Sens,Sens))))
        } else{
          Svalues <- cbind(rmultinom(N0,1,adjustprob(c(Spec,1-FP1-Spec,FP1))),
                           rmultinom(N1,1,adjustprob(c(FN2,1-FP2-FN2,FP2))),
                           rmultinom(N2,1,adjustprob(c(FN1,1-FN1-Sens,Sens))))
        }
        rownames(Svalues) <- c("S=0","S=1","S=2")
        Svalues <- ifelse(Svalues[1,]==1,0,ifelse(Svalues[2,]==1,1,2))
        return(Svalues)
      }
      
      SpecSens <- cbind(Spec,Sens,FP1,FN1,FP2,FN2)
      
      if ((P0+P2)==1) { # binary case only
        
        S <- t(apply(SpecSens, 1, function(x) biomarkerGroups(x, binary=TRUE))) # each row is a set of Sens, Spec, etc. parameters
      
      } else { # trichotomous
        
        S <- t(apply(SpecSens, 1, function(x) biomarkerGroups(x, binary=FALSE))) # each row is a set of Sens, Spec, etc. parameters
      
      }
      
      # Now keep the S's in nPhase2 of the cases (deleting the rest) and in controlCaseRatio*nPhase2 controls
      casesinds <- c(1:N)[Y==1]
      keepcasesinds <- sample(casesinds,nPhase2,replace=FALSE)
      controlinds <- c(1:N)[Y==0]
      keepcontrolinds <- sample(controlinds,controlCaseRatio*nPhase2,replace=FALSE)
      keepinds <- sort(c(keepcasesinds,keepcontrolinds))
      
      # Those with data:
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
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
          pval <- round(min(2*(1-pnorm(abs(fit$coef[2]/sqrt(fit$covm[2,2])))),1.0),4)
          if (pval <= alpha & fit$coef[2] < 0) { powerstrinary[k,j] <- powerstrinary[k,j] + 1}
        }
      }
      
    }
    
    # Continuous biomarker
    
    if (!is.null(VElowest)){
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
        
        # Now keep the S's in nPhase2 of the cases (deleting the rest) and in controlCaseRatio*nPhase2 controls
        casesinds <- c(1:N)[Y==1]
        keepcasesinds <- sample(casesinds,nPhase2,replace=FALSE)
        controlinds <- c(1:N)[Y==0]
        keepcontrolinds <- sample(controlinds,controlCaseRatio*nPhase2,replace=FALSE)
        keepinds <- sort(c(keepcasesinds,keepcontrolinds))
        
        # Those with data:
        Ycc <- Y[keepinds]
        Scc <- t(apply(S,1, function(x) x[keepinds])) # nrow=length(rhos)
        
        for(k in 1:nrow(Scc)){
          fit <- tps(Ycc~Scc[k,],nn0=length(Y[Y==0]),nn1=length(Y[Y==1]),group=rep(1,length(Ycc)))
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
# FP1, FP2, FN1, FN2 (defined above)
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
# are relevant (FP1, FP2, FN1, FN2 are not used in the calculations)
##################################################

# Original function that requires Plat2 to be a scalar
computeSensSpecFPFN <- function(sigma2obs,rhos,Plat0,Plat2,P0,P2) {
  # sigma2tr = Var(X) = Var(Str) = rho*sigma2obs
  # If Plat0 + Plat2 = 1 then the method collapses to a binary biomarker,
  # and FP1, FP2, FN1, FN2 are irrelevant; this function simply returns 0's in that scenario
  
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
    FP1 <- rep(0,m)
    FP2 <- rep(0,m)
    FN1 <- rep(0,m)
    FN2 <- rep(0,m)
    tauhisolution <- rep(0,m)
    taulosolution <- rep(0,m)
    if (rhos[i] < 1) {  #*# if rho=1, then Sens=1, Spec=1, FP1=0, FP2=0, FN1=0, FN2=0
      # Stochastic integration
      X <- rnorm(20000,0,sqrt(sigma2tr[i]))
      S <- X + rnorm(20000,0,sqrt(sigma2e[i]))
      
      Phi <- sum(X>thetahiVE[i])/length(X)
      Plo <- sum(X<=thetaloVE[i])/length(X)
      Pmed <- 1 - Phi - Plo
      
      # Compute a kernel over a grid of tau values:
      n <- 10000
      tauhi <- seq(-2.5,2.5,len=n)
      taulo <- seq(-2.5,2.5,len=n)
      Sensvec <- rep(1,n)
      Specvec <- rep(1,n)
      FP1vec <- rep(0,n)
      FP2vec <- rep(0,n)
      FN1vec <- rep(0,n)
      FN2vec <- rep(0,n)
      
      for (j in 1:length(tauhi)) {
        Sensvec[j] <- (sum(S>tauhi[j] & X > thetahiVE[i])/length(S))/Phi
        Specvec[j] <- (sum(S<=taulo[j] & X <= thetaloVE[i])/length(S))/Plo
        if (Pmed==0) { #*# if binary biomarker, 0's for FP2, FP1, FN2, FN1
          FP2vec[j] <- 0
          FP1vec[j] <- 0
          FN1vec[j] <- 0
          FN2vec[j] <- 0 }
        if (Pmed > 0) {
          FP2vec[j] <- (sum(S>tauhi[j] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
          FP1vec[j] <- (sum(S>tauhi[j] & X <= thetaloVE[i])/length(S))/Plo
          FN1vec[j] <- (sum(S<=taulo[j] & X > thetahiVE[i])/length(S))/Phi
          FN2vec[j] <- (sum(S<=taulo[j] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        }
      }
      
      for (l in 1:m) {
        kernelhi <- Sensvec*Plat2 + FP2vec*Plat1 + FP1vec*Plat0 - P2[l]
        kernelhi <- kernelhi^2
        
        kernello <- Specvec*Plat0 + FN2vec*Plat1 + FN1vec*Plat2 - P0[l]
        kernello <- kernello^2
        
        indhi <- c(1:length(kernelhi))[kernelhi==min(kernelhi)] #*# equiv to which(kernelhi==min(kernelhi))
        tauhisol <- tauhi[indhi]
        indlo <- c(1:length(kernello))[kernello==min(kernello)]
        taulosol <- taulo[indlo]
        
        if (length(tauhisol)>1) {tauhisol <- tauhisol[1]} #*# if more than one tau given, choose first one
        if (length(taulosol)>1) {taulosol <- taulosol[1]}
        
        tauhisolution[l] <- tauhisol
        taulosolution[l] <- taulosol
        
        Sens[l] <- sum(S>tauhisolution[l] & X > thetahiVE[i])/sum(X>thetahiVE[i])
        Spec[l] <- sum(S<=taulosolution[l] & X <= thetaloVE[i])/sum(X<=thetaloVE[i])
        if (Pmed==0) {  #*# if binary biomarker, 0's for FP2, FP1, FN2, FN1
          FP2[l] <- 0
          FP1[l] <- 0
          FN2[l] <- 0
          FN1[l] <- 0 }
        if (Pmed > 0) {
          FP2[l] <- (sum(S>tauhisolution[l] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
          FP1[l] <- (sum(S>tauhisolution[l] & X <= thetaloVE[i])/length(S))/Plo
          FN1[l] <- (sum(S<=taulosolution[l] & X > thetahiVE[i])/length(S))/Phi
          FN2[l] <- (sum(S<=taulosolution[l] & X > thetaloVE[i] & X <= thetahiVE[i])/length(S))/Pmed
        }
      }
    }
    ans[[i]] <- cbind(rep(thetaloVE[i],m),rep(thetahiVE[i],m),rep(Plat0,m),rep(Plat1,m),rep(Plat2,m),P0,P2,
                      taulosolution,tauhisolution,Sens,Spec,FP1,FP2,FN1,FN2)
  }  
  return(ans)
}
